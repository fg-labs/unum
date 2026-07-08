//! Binary-level end-to-end tests for `unum extract --bam-mode`: proves that
//! `-i x.bam --bam-mode alignment` and `-b x.bam --bam-mode alignment` reach
//! the SAME coordinate-`alignment` extraction path, that `--bam-mode
//! no-alignment` routes by sort order (coordinate -> 2-pass, name-sorted ->
//! grouped one-pass) instead of erroring, and that every remaining
//! reserved/guard combination (missing `--bam-mode`, `alignment` on a
//! non-coordinate-sorted BAM) errors with a helpful, on-topic message.
//!
//! The coordinate `alignment` extraction itself is already byte-golden-gated
//! at the library level by `unum-core/tests/golden_bam_extract.rs`, so this
//! file does not re-freeze a new golden: its job is asserting the `-i`/`-b`
//! CLI invocations agree with EACH OTHER and are non-empty, which proves both
//! routes reach the real dispatcher rather than, say, two trivially-matching
//! empty files.
//!
//! BAMs are built programmatically via `rust_htslib::bam::Writer` in a
//! `tempfile::tempdir()` and never committed.
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{self, Header, Writer};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use unum_core::bam_extract::{self, CoordRecord};

/// A path under the workspace-level `fixtures/` directory (a sibling of
/// `crates/`), mirroring `extract_e2e.rs`'s `example`/`golden` helpers and
/// `unum-core/tests/golden_bam_extract.rs`'s `fixture` helper.
fn fixture(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures").join(rel)
}

/// The coord FASTA `-f` must point at for BAM mode: the SAME fixture
/// `unum-core/tests/golden_bam_extract.rs` builds its BAM header against
/// (`fixture("refbuild/golden/kir_rna_coord.fa")`), so `build_genes`'s
/// chrom-name resolution succeeds against the BAMs this file builds below
/// (both declare a single `chr19` `@SQ`, matching the coord FASTA's only
/// referenced chrom).
fn coord_fasta_path() -> PathBuf {
    fixture("refbuild/golden/kir_rna_coord.fa")
}

/// `chr19`'s KIR2DL1 gene-interval start used by
/// `unum-core/tests/golden_bam_extract.rs`'s BAM builder; reused here so the
/// on-target pair this file builds lands inside the real gene interval
/// `build_genes` resolves from the coord FASTA.
const KIR2DL1_CHROM: &str = "chr19";
const KIR2DL1_START: i64 = 54_769_793;

/// The `KIR2DL1*0010101` sequence bytes from the coord FASTA, mirroring
/// `unum-core/tests/golden_bam_extract.rs::kir2dl1_sequence` -- gives the
/// on-target pair this file's BAM builder writes real, filter-clearing KIR
/// bases instead of arbitrary noise.
fn kir2dl1_sequence() -> String {
    bam_extract::parse_coord_fa(&coord_fasta_path())
        .expect("parse coord fa")
        .into_iter()
        .find(|r: &CoordRecord| r.name == "KIR2DL1*0010101")
        .expect("KIR2DL1*0010101 present in coord FASTA")
        .seq
}

/// Writes a coordinate-sorted (`@HD SO:coordinate`) or name-sorted (`@HD
/// SO:queryname`) BAM containing a single on-target read pair, whose
/// sequence is a real substring of `KIR2DL1*0010101` long enough (100bp) to
/// clear the default paired `hitLenRequired` filter and produce non-empty
/// candidate output. Mirrors the builder pattern in
/// `unum-core/tests/golden_bam_extract.rs`'s `build_test_bam` (hand-set
/// tid/pos/mtid/mpos/flags), trimmed to the one on-target pair this file's
/// tests need -- the sort-order guard fires on the `@HD` tag before any
/// record is read, so a single pair (order-independent) is sufficient for
/// every test in this file, including the name-sorted guard test.
fn build_bam(dir: &Path, filename: &str, sort_order_tag: &str) -> PathBuf {
    let path = dir.join(filename);
    let kir_seq = kir2dl1_sequence();
    let kir_bytes = kir_seq.as_bytes();

    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", sort_order_tag);
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", KIR2DL1_CHROM);
    sq.push_tag(b"LN", 58_617_616);
    header.push_record(&sq);

    let mut writer = Writer::from_path(&path, &header, bam::Format::Bam).expect("bam writer");

    let mut r1 = bam::Record::new();
    r1.set(
        b"on_target_pair",
        Some(&CigarString(vec![Cigar::Match(100)])),
        &kir_bytes[0..100],
        &[30u8; 100],
    );
    r1.set_tid(0);
    r1.set_pos(KIR2DL1_START + 10);
    r1.set_mtid(0);
    r1.set_mpos(KIR2DL1_START + 210);
    r1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&r1).expect("write mate 1");

    let mut r2 = bam::Record::new();
    r2.set(
        b"on_target_pair",
        Some(&CigarString(vec![Cigar::Match(100)])),
        &kir_bytes[200..300],
        &[30u8; 100],
    );
    r2.set_tid(0);
    r2.set_pos(KIR2DL1_START + 210);
    r2.set_mtid(0);
    r2.set_mpos(KIR2DL1_START + 10);
    r2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&r2).expect("write mate 2");

    drop(writer);
    path
}

/// A coordinate-sorted (`@HD SO:coordinate`) BAM -- the input the Stage 2a
/// `alignment` path accepts.
fn build_coordinate_bam(dir: &Path) -> PathBuf {
    build_bam(dir, "coordinate.bam", "coordinate")
}

/// The SAME on-target pair as [`build_coordinate_bam`], but under `@HD
/// SO:queryname` -- exercises the sort-order guard in
/// `crate::stages::extract::run_bam_mode`, which rejects `alignment` on
/// anything but `SO:coordinate`.
fn build_name_sorted_bam(dir: &Path) -> PathBuf {
    build_bam(dir, "name_sorted.bam", "queryname")
}

/// The SAME on-target pair as [`build_coordinate_bam`], but under `@HD
/// SO:unsorted` -- exercises the distinct `SortOrder::Unsorted` branch of the
/// sort-order guard in `crate::stages::extract::run_bam_mode`, whose error
/// message differs from the name-sorted (`SO:queryname`) one. `SO:unsorted`
/// (like a missing `@HD SO`) maps to [`SortOrder::Unsorted`].
fn build_unsorted_bam(dir: &Path) -> PathBuf {
    build_bam(dir, "unsorted.bam", "unsorted")
}

/// `Path` -> `&str`, mirroring `extract_e2e.rs`'s `path_str` (named `s` here
/// to match the task brief's helper name).
fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

/// Runs `unum extract` with `args`, asserting a clean exit. Mirrors
/// `extract_e2e.rs`'s `run_extract`; renamed `run_extract_ok` here to pair
/// with [`run_extract_raw`], which the guard tests below need in order to
/// inspect a FAILING run's `stderr`.
fn run_extract_ok(args: &[&str]) {
    let status =
        Command::new(env!("CARGO_BIN_EXE_unum")).arg("extract").args(args).status().unwrap();
    assert!(status.success(), "`unum extract` exited non-zero for args {args:?}");
}

/// Runs `unum extract` with `args`, returning the full captured `Output`
/// (exit status + stdout/stderr) without asserting success -- for the
/// guard-path tests below, which assert on a non-zero exit plus a specific
/// `stderr` message.
fn run_extract_raw(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_unum")).arg("extract").args(args).output().unwrap()
}

/// The `{prefix}_1.fq` output of a paired `extract` run, mirroring
/// `extract_e2e.rs`'s `paired_out_1`. The BAMs this file builds are paired
/// (mate-paired flags on both records), so `alignment` extraction emits
/// `{prefix}_1.fq`/`{prefix}_2.fq`, not the single-end `{prefix}.fq`.
fn paired_out_1(prefix: &Path) -> PathBuf {
    PathBuf::from(format!("{}_1.fq", prefix.to_str().unwrap()))
}

/// The `{prefix}_2.fq` (mate-2) output of a paired `extract` run, mirroring
/// [`paired_out_1`]. Both mates are compared so a divergence on mate-2 output
/// between the `-i` and `-b` routes cannot slip through.
fn paired_out_2(prefix: &Path) -> PathBuf {
    PathBuf::from(format!("{}_2.fq", prefix.to_str().unwrap()))
}

#[test]
fn i_and_b_flag_coordinate_alignment_match_and_are_nonempty() {
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let bam = build_coordinate_bam(tmp.path());

    let out_i = tmp.path().join("via_i");
    run_extract_ok(&[
        "-c",
        s(&coord_fa),
        "-i",
        s(&bam),
        "--bam-mode",
        "alignment",
        "-o",
        s(&out_i),
    ]);

    let out_b = tmp.path().join("via_b");
    run_extract_ok(&[
        "-c",
        s(&coord_fa),
        "-b",
        s(&bam),
        "--bam-mode",
        "alignment",
        "-o",
        s(&out_b),
    ]);

    let a1 = std::fs::read(paired_out_1(&out_i)).unwrap();
    let b1 = std::fs::read(paired_out_1(&out_b)).unwrap();
    let a2 = std::fs::read(paired_out_2(&out_i)).unwrap();
    let b2 = std::fs::read(paired_out_2(&out_b)).unwrap();
    assert_eq!(a1, b1, "-i and -b coordinate alignment must be byte-identical (mate 1)");
    assert_eq!(a2, b2, "-i and -b coordinate alignment must be byte-identical (mate 2)");
    assert!(!a1.is_empty(), "coordinate alignment should emit mate-1 candidates");
    assert!(!a2.is_empty(), "coordinate alignment should emit mate-2 candidates");
}

#[test]
fn b_flag_without_bam_mode_errors() {
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let bam = build_coordinate_bam(tmp.path());
    let out = run_extract_raw(&["-f", s(&coord_fa), "-b", s(&bam), "-o", s(&tmp.path().join("o"))]);
    assert!(!out.status.success(), "-b without --bam-mode must fail");
    assert!(String::from_utf8_lossy(&out.stderr).contains("--bam-mode"));
}

#[test]
fn no_alignment_mode_on_coordinate_bam_runs_the_two_pass() {
    // `--bam-mode no-alignment` is no longer reserved (Stage 2b): on a
    // coordinate-sorted BAM it routes to the seekable 2-pass name-map
    // (`bam_extract::extract_from_bam_no_alignment`), same input this file's
    // `alignment` tests use, reusing the coord FASTA directly as `-f` (
    // no-alignment builds its `RefKmerFilter` straight from `-f`, same as the
    // FASTQ path -- no coord-FASTA/gene-interval parsing).
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let bam = build_coordinate_bam(tmp.path());
    let out_prefix = tmp.path().join("o");
    run_extract_ok(&[
        "-f",
        s(&coord_fa),
        "-i",
        s(&bam),
        "--bam-mode",
        "no-alignment",
        "-o",
        s(&out_prefix),
    ]);
    let out1 = std::fs::read(paired_out_1(&out_prefix)).unwrap();
    assert!(!out1.is_empty(), "no-alignment 2-pass should emit the on-target candidate pair");
}

#[test]
fn no_alignment_mode_on_name_sorted_bam_runs_the_grouped_one_pass() {
    // The SAME on-target pair, but `@HD SO:queryname` -- no-alignment routes
    // this to the stdin-capable grouped one-pass
    // (`bam_extract::extract_from_bam_no_alignment_grouped`) instead of the
    // 2-pass, unlike `alignment` (which still rejects this input as reserved
    // for Stage 2c -- see `alignment_on_name_sorted_bam_errors_with_sort_hint`).
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let bam = build_name_sorted_bam(tmp.path());
    let out_prefix = tmp.path().join("o");
    run_extract_ok(&[
        "-f",
        s(&coord_fa),
        "-i",
        s(&bam),
        "--bam-mode",
        "no-alignment",
        "-o",
        s(&out_prefix),
    ]);
    let out1 = std::fs::read(paired_out_1(&out_prefix)).unwrap();
    assert!(
        !out1.is_empty(),
        "no-alignment grouped one-pass should emit the on-target candidate pair"
    );
}

#[test]
fn alignment_on_name_sorted_bam_errors_with_sort_hint() {
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let bam = build_name_sorted_bam(tmp.path());
    let out = run_extract_raw(&[
        "-c",
        s(&coord_fa),
        "-i",
        s(&bam),
        "--bam-mode",
        "alignment",
        "-o",
        s(&tmp.path().join("o")),
    ]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("coordinate") || err.contains("samtools sort"));
}

#[test]
fn alignment_on_unsorted_bam_errors_with_sort_hint() {
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let bam = build_unsorted_bam(tmp.path());
    let out = run_extract_raw(&[
        "-f",
        s(&coord_fa),
        "-i",
        s(&bam),
        "--bam-mode",
        "alignment",
        "-o",
        s(&tmp.path().join("o")),
    ]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    // The `SortOrder::Unsorted` branch carries its own message distinct from
    // the name-sorted one: it names the unsorted/unstated state and hints at
    // `samtools sort`.
    assert!(err.contains("unsorted"), "expected unsorted-state message, got: {err}");
    assert!(
        err.contains("samtools sort") || err.contains("coordinate"),
        "expected a coordinate-sort hint, got: {err}"
    );
}

// -- `no-alignment ≡ FASTQ` CLI e2e (Stage 2b Task 7) -----------------------
//
// Proves `unum extract -i x.bam --bam-mode no-alignment` produces the SAME
// candidates as `unum extract -1/-2` (FASTQ mode) on the identical
// underlying reads, at the real CLI/process level -- the library-level
// equivalence (against `bam_extract::extract_from_bam_no_alignment(_grouped)`
// directly) is proven in `unum-core/tests/golden_bam_no_alignment.rs`; this
// file's job is asserting the CLI wiring (arg parsing, dispatcher routing,
// `FastqFileSink` naming) doesn't disturb that equivalence. Reuses this
// file's own `kir2dl1_sequence`/`coord_fasta_path` fixtures. Fixtures are
// CLEAN and homogeneous (all-paired, uniform read length, forward-strand
// only -- see `golden_bam_no_alignment.rs`'s module docs for why: mixed/
// orphan input is brief-sanctioned to diverge between BAM modes, and a
// `0x10`-flagged record's SEQ would be reverse-complemented by
// `Alignments::read_seq`, diverging from a raw-byte FASTQ fixture unless
// compensated for).

/// Off-reference noise, distinct from KIR2DL1's composition (mirrors
/// `golden_bam_no_alignment.rs`'s own `noise` helper).
fn cli_noise(len: usize) -> Vec<u8> {
    let pattern =
        b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
    (0..len).map(|i| pattern[i % pattern.len()]).collect()
}

const CLI_PAIR_LEN: usize = 100;

/// One paired-end template for the CLI equivalence tests: UNIFORM
/// `CLI_PAIR_LEN`-byte mate sequences, all-paired (no orphans).
struct CliPairTemplate {
    id: &'static str,
    seq1: Vec<u8>,
    seq2: Vec<u8>,
}

/// Three templates: OR-rescue (mate1 hits, mate2 is noise), both-hit, and
/// neither-hit (must be dropped) -- proves the CLI e2e isn't a vacuous
/// "everything passes"/"everything empty" comparison.
fn cli_pair_templates() -> Vec<CliPairTemplate> {
    let kir_seq = kir2dl1_sequence();
    let kir_bytes = kir_seq.as_bytes();
    vec![
        CliPairTemplate {
            id: "hit_rescue_pair",
            seq1: kir_bytes[0..CLI_PAIR_LEN].to_vec(),
            seq2: cli_noise(CLI_PAIR_LEN),
        },
        CliPairTemplate {
            id: "hit_both_pair",
            seq1: kir_bytes[200..200 + CLI_PAIR_LEN].to_vec(),
            seq2: kir_bytes[400..400 + CLI_PAIR_LEN].to_vec(),
        },
        CliPairTemplate {
            id: "fail_pair",
            seq1: cli_noise(CLI_PAIR_LEN),
            seq2: cli_noise(CLI_PAIR_LEN),
        },
    ]
}

/// Writes one FASTQ record with an all-`I` quality string.
fn write_cli_fastq_record(f: &mut std::fs::File, id: &str, seq: &[u8]) {
    let qual = "I".repeat(seq.len());
    writeln!(f, "@{id}\n{}\n+\n{qual}", std::str::from_utf8(seq).unwrap()).unwrap();
}

/// Writes `templates` as paired FASTQ files, mate1/mate2 in template order.
fn write_cli_paired_fastq(dir: &Path, templates: &[CliPairTemplate]) -> (PathBuf, PathBuf) {
    let r1 = dir.join("cli_r1.fq");
    let r2 = dir.join("cli_r2.fq");
    let mut f1 = std::fs::File::create(&r1).unwrap();
    let mut f2 = std::fs::File::create(&r2).unwrap();
    for t in templates {
        write_cli_fastq_record(&mut f1, t.id, &t.seq1);
        write_cli_fastq_record(&mut f2, t.id, &t.seq2);
    }
    (r1, r2)
}

/// Writes `templates` as a paired BAM (`@HD SO:<so_val>`), mate1 immediately
/// followed by mate2 per template, templates in the SAME order given -- the
/// order [`write_cli_paired_fastq`] itself iterates in, so a
/// `SO:queryname`-tagged write is a valid grouped fixture that matches the
/// FASTQ mate files' record order exactly. Every record is forward-strand
/// only (no `0x10`/`0x20`) and carries raw QUAL `40` (ASCII `'I'` after
/// `Alignments::qual`'s `+33` conversion), so the BAM no-alignment path's
/// emitted bytes can be byte-identical to this file's `'I'`-quality FASTQ
/// fixture (see the section doc comment above).
fn write_cli_paired_bam(path: &Path, templates: &[CliPairTemplate], so_val: &str) {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", so_val);
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", KIR2DL1_CHROM);
    sq.push_tag(b"LN", 58_617_616);
    header.push_record(&sq);

    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("bam writer");

    for (i, t) in templates.iter().enumerate() {
        let pos1 = KIR2DL1_START + 10 + i64::try_from(i).unwrap() * 1_000;
        let pos2 = pos1 + 300;
        let qual1 = vec![40u8; t.seq1.len()];
        let qual2 = vec![40u8; t.seq2.len()];

        let mut r1 = bam::Record::new();
        r1.set(
            t.id.as_bytes(),
            Some(&CigarString(vec![Cigar::Match(u32::try_from(t.seq1.len()).unwrap())])),
            &t.seq1,
            &qual1,
        );
        r1.set_tid(0);
        r1.set_pos(pos1);
        r1.set_mtid(0);
        r1.set_mpos(pos2);
        r1.set_flags(0x1 | 0x40);
        writer.write(&r1).unwrap();

        let mut r2 = bam::Record::new();
        r2.set(
            t.id.as_bytes(),
            Some(&CigarString(vec![Cigar::Match(u32::try_from(t.seq2.len()).unwrap())])),
            &t.seq2,
            &qual2,
        );
        r2.set_tid(0);
        r2.set_pos(pos2);
        r2.set_mtid(0);
        r2.set_mpos(pos1);
        r2.set_flags(0x1 | 0x80);
        writer.write(&r2).unwrap();
    }

    drop(writer);
}

/// Parses a candidate FASTQ's text into SORTED `(id, "seq\n+\nqual")` records,
/// for comparing the coordinate 2-pass's output to FASTQ mode's as a MULTISET
/// (pass 2 emits in pair-completion order, which need not match FASTQ's input
/// order). Sorting keeps the comparison order-independent while preserving
/// duplicate-record multiplicity -- a plain `id -> record` map would silently
/// collapse repeated ids, hiding a dropped or duplicated candidate.
fn fastq_records_sorted_by_id(text: &str) -> Vec<(String, String)> {
    let lines: Vec<&str> = text.lines().collect();
    let mut records: Vec<(String, String)> = lines
        .chunks(4)
        .map(|chunk| {
            let [id_line, seq, plus, qual] = chunk else {
                panic!("malformed FASTQ chunk: {chunk:?}")
            };
            let id = id_line.strip_prefix('@').expect("FASTQ id line must start with '@'");
            (id.to_string(), format!("{seq}\n{plus}\n{qual}"))
        })
        .collect();
    records.sort();
    records
}

#[test]
fn no_alignment_bam_equals_fastq_cli_coordinate() {
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let templates = cli_pair_templates();
    let (r1, r2) = write_cli_paired_fastq(tmp.path(), &templates);
    let bam = tmp.path().join("cli_coord.bam");
    write_cli_paired_bam(&bam, &templates, "coordinate");

    let fastq_prefix = tmp.path().join("fastq_out");
    run_extract_ok(&["-f", s(&coord_fa), "-1", s(&r1), "-2", s(&r2), "-o", s(&fastq_prefix)]);

    let bam_prefix = tmp.path().join("bam_out");
    run_extract_ok(&[
        "-f",
        s(&coord_fa),
        "-i",
        s(&bam),
        "--bam-mode",
        "no-alignment",
        "-o",
        s(&bam_prefix),
    ]);

    let fastq_out1 = std::fs::read_to_string(paired_out_1(&fastq_prefix)).unwrap();
    let fastq_out2 = std::fs::read_to_string(paired_out_2(&fastq_prefix)).unwrap();
    let bam_out1 = std::fs::read_to_string(paired_out_1(&bam_prefix)).unwrap();
    let bam_out2 = std::fs::read_to_string(paired_out_2(&bam_prefix)).unwrap();

    // Sanity/teeth: the expected pass/fail mix actually happened.
    assert!(fastq_out1.contains("hit_rescue_pair"));
    assert!(fastq_out1.contains("hit_both_pair"));
    assert!(!fastq_out1.contains("fail_pair"));
    assert!(!bam_out2.is_empty(), "no-alignment 2-pass should emit mate-2 candidates");

    // Coordinate 2-pass emits in pair-COMPLETION order -- compare EACH mate's
    // candidate FASTQ as a MULTISET keyed by id (sorted records), not
    // byte-for-byte: reordering is tolerated, but a dropped or duplicated
    // record is still caught (which an `id -> record` map would hide).
    assert_eq!(
        fastq_records_sorted_by_id(&fastq_out1),
        fastq_records_sorted_by_id(&bam_out1),
        "no-alignment coordinate BAM mate-1 candidates diverged from the FASTQ-mode CLI run"
    );
    assert_eq!(
        fastq_records_sorted_by_id(&fastq_out2),
        fastq_records_sorted_by_id(&bam_out2),
        "no-alignment coordinate BAM mate-2 candidates diverged from the FASTQ-mode CLI run"
    );
}

#[test]
fn no_alignment_bam_equals_fastq_cli_grouped() {
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let templates = cli_pair_templates();
    let (r1, r2) = write_cli_paired_fastq(tmp.path(), &templates);
    let bam = tmp.path().join("cli_grouped.bam");
    write_cli_paired_bam(&bam, &templates, "queryname");

    let fastq_prefix = tmp.path().join("fastq_out");
    run_extract_ok(&["-f", s(&coord_fa), "-1", s(&r1), "-2", s(&r2), "-o", s(&fastq_prefix)]);

    let bam_prefix = tmp.path().join("bam_out");
    run_extract_ok(&[
        "-f",
        s(&coord_fa),
        "-i",
        s(&bam),
        "--bam-mode",
        "no-alignment",
        "-o",
        s(&bam_prefix),
    ]);

    let fastq_out1 = std::fs::read(paired_out_1(&fastq_prefix)).unwrap();
    let fastq_out2 = std::fs::read(paired_out_2(&fastq_prefix)).unwrap();
    let bam_out1 = std::fs::read(paired_out_1(&bam_prefix)).unwrap();
    let bam_out2 = std::fs::read(paired_out_2(&bam_prefix)).unwrap();

    assert!(!fastq_out1.is_empty(), "should emit at least one candidate");

    // The name-sorted BAM preserves the FASTQ mate files' record order
    // exactly (mate1/mate2 adjacent per template, templates in input order),
    // so the grouped one-pass's output is byte-identical to FASTQ mode's.
    assert_eq!(
        fastq_out1, bam_out1,
        "no-alignment grouped BAM _1.fq diverged from the FASTQ-mode CLI run"
    );
    assert_eq!(
        fastq_out2, bam_out2,
        "no-alignment grouped BAM _2.fq diverged from the FASTQ-mode CLI run"
    );
}

/// The SAME grouped/name-sorted no-alignment equivalence as
/// [`no_alignment_bam_equals_fastq_cli_grouped`], but fed via STDIN (`-i -`)
/// instead of a file path. This exercises the separate stdin CLI branch
/// (`BamInputSpec::Stdin`) and the non-rewindable `Alignments::from_stdin`
/// one-pass, which the file-backed test does not reach. The grouped one-pass is
/// order-preserving, so the stdin output stays byte-identical to FASTQ mode's,
/// same as the file case.
#[test]
fn no_alignment_bam_equals_fastq_cli_grouped_stdin() {
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let templates = cli_pair_templates();
    let (r1, r2) = write_cli_paired_fastq(tmp.path(), &templates);
    let bam = tmp.path().join("cli_grouped_stdin.bam");
    write_cli_paired_bam(&bam, &templates, "queryname");

    let fastq_prefix = tmp.path().join("fastq_out");
    run_extract_ok(&["-f", s(&coord_fa), "-1", s(&r1), "-2", s(&r2), "-o", s(&fastq_prefix)]);

    // Pipe the BAM to `extract -i - --bam-mode no-alignment` on stdin.
    let bam_prefix = tmp.path().join("bam_out");
    let bam_file = std::fs::File::open(&bam).unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_unum"))
        .arg("extract")
        .args(["-f", s(&coord_fa), "-i", "-", "--bam-mode", "no-alignment", "-o", s(&bam_prefix)])
        .stdin(std::process::Stdio::from(bam_file))
        .status()
        .unwrap();
    assert!(
        status.success(),
        "`unum extract -i - --bam-mode no-alignment` (stdin) exited non-zero"
    );

    let fastq_out1 = std::fs::read(paired_out_1(&fastq_prefix)).unwrap();
    let fastq_out2 = std::fs::read(paired_out_2(&fastq_prefix)).unwrap();
    let bam_out1 = std::fs::read(paired_out_1(&bam_prefix)).unwrap();
    let bam_out2 = std::fs::read(paired_out_2(&bam_prefix)).unwrap();

    assert!(!bam_out1.is_empty(), "stdin grouped no-alignment should emit candidates");
    assert_eq!(
        fastq_out1, bam_out1,
        "no-alignment grouped BAM (stdin) _1.fq diverged from the FASTQ-mode CLI run"
    );
    assert_eq!(
        fastq_out2, bam_out2,
        "no-alignment grouped BAM (stdin) _2.fq diverged from the FASTQ-mode CLI run"
    );
}
