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
use std::process::{Command, Output, Stdio};
use std::time::{Duration, Instant};
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

/// A coordinate-sorted BAM with `copies` DUPLICATE on-target pairs (same
/// positions/sequences as [`build_coordinate_bam`]'s single pair, distinct
/// `QNAME`s), for the `run -b` fused-genotyping test below: a single read
/// pair's EM-estimated abundance falls under the genotyper's default
/// `filter_cov` (1.0) threshold (`Genotyper::select_alleles_for_genes_quality_scores`
/// zeroes `genotype_quality` whenever a rank's summed abundance is below
/// `filter_cov`), which in turn makes `output_representative_alleles` skip
/// writing ANY representative allele (`genotype_quality < 1` gate) --
/// leaving `_allele.tsv` empty and `analyze` erroring with "reference FASTA
/// contains no sequences matching the selected-allele list". A few
/// duplicate copies push the summed abundance comfortably past `filter_cov`
/// so a real representative allele gets selected and `analyze` has
/// something to work with, without changing [`build_bam`] (shared by other
/// tests in this file that rely on its single-pair shape).
fn build_coordinate_bam_replicated(dir: &Path, copies: u32) -> PathBuf {
    let path = dir.join("coordinate_replicated.bam");
    let kir_seq = kir2dl1_sequence();
    let kir_bytes = kir_seq.as_bytes();

    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", KIR2DL1_CHROM);
    sq.push_tag(b"LN", 58_617_616);
    header.push_record(&sq);

    let mut writer = Writer::from_path(&path, &header, bam::Format::Bam).expect("bam writer");
    for copy in 0..copies {
        let qname = format!("on_target_pair_{copy}");

        let mut r1 = bam::Record::new();
        r1.set(
            qname.as_bytes(),
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
            qname.as_bytes(),
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
    }
    drop(writer);
    path
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

// -- CRAM input (Stage 2c Task 4) -------------------------------------------
//
// Proves `unum extract` accepts CRAM the same way it accepts BAM -- routed
// through the SAME `--bam-mode` paths (it's just a codec) -- given an
// explicit `-r` reference, and rejects CRAM outright when `-r` is missing.
// Reuses the CRAM-writing + `.fai` fixture recipe from
// `unum-core::alignments::tests` (`write_fasta_with_fai`/`write_cram`,
// verified there against `samtools faidx`), duplicated here since it's a
// crate-private test helper there.

/// Writes `ref.fa` + `ref.fa.fai` (contig `chr1`, `len` bytes of a cycled
/// `ACGT` sequence) and returns the path plus the raw sequence bytes. Used as
/// BOTH `-r` (the CRAM writer/reader's decode reference) and `-f` (the
/// `RefKmerFilter` reference-sequence FASTA): the reads
/// [`write_cram_bam_reads`] builds are an exact substring of this same
/// sequence, so they trivially clear the filter's default similarity
/// threshold without needing a second, unrelated reference fixture.
///
/// The `.fai` is written by hand as `chr1\t<len>\t<offset>\t<len>\t<len+1>`
/// (offset 6 = `>chr1\n`; the whole sequence on one line) -- the exact
/// recipe `unum-core::alignments::tests::write_fasta_with_fai` uses,
/// independently verified there against `samtools faidx`.
fn write_ref_and_fai(dir: &Path, len: usize) -> (PathBuf, Vec<u8>) {
    let ref_fa = dir.join("ref.fa");
    let seq: Vec<u8> = b"ACGT".iter().copied().cycle().take(len).collect();

    let mut fasta = Vec::new();
    fasta.push(b'>');
    fasta.extend_from_slice(b"chr1");
    fasta.push(b'\n');
    fasta.extend_from_slice(&seq);
    fasta.push(b'\n');
    std::fs::write(&ref_fa, &fasta).unwrap();

    let offset = "chr1".len() + 2; // ">" + "chr1" + "\n"
    let fai_path = PathBuf::from(format!("{}.fai", ref_fa.display()));
    let fai_line = format!("chr1\t{}\t{offset}\t{}\t{}\n", seq.len(), seq.len(), seq.len() + 1);
    std::fs::write(&fai_path, fai_line).unwrap();

    (ref_fa, seq)
}

/// A single-contig `chr1` header (`LN = ref_len`), shared by the BAM and
/// CRAM fixtures below so their headers -- and therefore the records read
/// back from each -- agree exactly.
fn cram_bam_header(ref_len: usize) -> Header {
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", i64::try_from(ref_len).unwrap());
    header.push_record(&sq);
    header
}

/// RFC 1321 per-round left-rotation amounts, used by [`md5_process_block`].
const MD5_SHIFTS: [u32; 64] = [
    7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9,
    14, 20, 5, 9, 14, 20, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 6, 10, 15,
    21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21,
];

/// RFC 1321 per-round additive constants (`floor(abs(sin(round+1)) * 2^32)`),
/// used by [`md5_process_block`].
const MD5_CONSTANTS: [u32; 64] = [
    0xd76a_a478,
    0xe8c7_b756,
    0x2420_70db,
    0xc1bd_ceee,
    0xf57c_0faf,
    0x4787_c62a,
    0xa830_4613,
    0xfd46_9501,
    0x6980_98d8,
    0x8b44_f7af,
    0xffff_5bb1,
    0x895c_d7be,
    0x6b90_1122,
    0xfd98_7193,
    0xa679_438e,
    0x49b4_0821,
    0xf61e_2562,
    0xc040_b340,
    0x265e_5a51,
    0xe9b6_c7aa,
    0xd62f_105d,
    0x0244_1453,
    0xd8a1_e681,
    0xe7d3_fbc8,
    0x21e1_cde6,
    0xc337_07d6,
    0xf4d5_0d87,
    0x455a_14ed,
    0xa9e3_e905,
    0xfcef_a3f8,
    0x676f_02d9,
    0x8d2a_4c8a,
    0xfffa_3942,
    0x8771_f681,
    0x6d9d_6122,
    0xfde5_380c,
    0xa4be_ea44,
    0x4bde_cfa9,
    0xf6bb_4b60,
    0xbebf_bc70,
    0x289b_7ec6,
    0xeaa1_27fa,
    0xd4ef_3085,
    0x0488_1d05,
    0xd9d4_d039,
    0xe6db_99e5,
    0x1fa2_7cf8,
    0xc4ac_5665,
    0xf429_2244,
    0x432a_ff97,
    0xab94_23a7,
    0xfc93_a039,
    0x655b_59c3,
    0x8f0c_cc92,
    0xffef_f47d,
    0x8584_5dd1,
    0x6fa8_7e4f,
    0xfe2c_e6e0,
    0xa301_4314,
    0x4e08_11a1,
    0xf753_7e82,
    0xbd3a_f235,
    0x2ad7_d2bb,
    0xeb86_d391,
];

/// Runs the RFC 1321 MD5 compression function over one 64-byte `block`,
/// updating the running hash `state` (`[A, B, C, D]`) in place. Split out of
/// [`md5_hex`] purely to stay under this workspace's clippy
/// `too_many_lines` pedantic threshold.
fn md5_process_block(state: &mut [u32; 4], block: &[u8]) {
    let mut words = [0u32; 16];
    for (word, bytes) in words.iter_mut().zip(block.chunks_exact(4)) {
        *word = u32::from_le_bytes(bytes.try_into().unwrap());
    }

    let (mut aa, mut bb, mut cc, mut dd) = (state[0], state[1], state[2], state[3]);
    for round in 0..64 {
        let (mix, word_idx) = if round < 16 {
            ((bb & cc) | (!bb & dd), round)
        } else if round < 32 {
            ((dd & bb) | (!dd & cc), (5 * round + 1) % 16)
        } else if round < 48 {
            (bb ^ cc ^ dd, (3 * round + 5) % 16)
        } else {
            (cc ^ (bb | !dd), (7 * round) % 16)
        };
        let mix =
            mix.wrapping_add(aa).wrapping_add(MD5_CONSTANTS[round]).wrapping_add(words[word_idx]);
        aa = dd;
        dd = cc;
        cc = bb;
        bb = bb.wrapping_add(mix.rotate_left(MD5_SHIFTS[round]));
    }

    state[0] = state[0].wrapping_add(aa);
    state[1] = state[1].wrapping_add(bb);
    state[2] = state[2].wrapping_add(cc);
    state[3] = state[3].wrapping_add(dd);
}

/// Pads `data` per RFC 1321 (a `0x80` byte, zero-fill, then the original
/// bit-length as a little-endian `u64`) so the result's length is a
/// multiple of 64 bytes, ready for [`md5_process_block`].
fn md5_pad(data: &[u8]) -> Vec<u8> {
    let mut padded = data.to_vec();
    let bit_len = (data.len() as u64).wrapping_mul(8);
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_len.to_le_bytes());
    padded
}

/// A minimal, from-scratch MD5 (RFC 1321) implementation, needed ONLY to
/// compute a SAM-spec `M5` tag (`md5(uppercase(reference bases))`, no
/// whitespace) for the stdin-CRAM regression test below -- no MD5/digest
/// crate is otherwise a dependency of this workspace, and pulling one in
/// for a single test-fixture checksum isn't worth it. Self-checked against
/// the RFC 1321 test vectors in [`md5_hex_matches_rfc1321_vectors`].
fn md5_hex(data: &[u8]) -> String {
    use std::fmt::Write as _;

    let mut state = [0x6745_2301u32, 0xefcd_ab89u32, 0x98ba_dcfeu32, 0x1032_5476u32];
    for block in md5_pad(data).chunks_exact(64) {
        md5_process_block(&mut state, block);
    }

    state.iter().flat_map(|word| word.to_le_bytes()).fold(
        String::with_capacity(32),
        |mut hex, byte| {
            let _ = write!(hex, "{byte:02x}");
            hex
        },
    )
}

#[test]
fn md5_hex_matches_rfc1321_vectors() {
    // RFC 1321 sec. A.5 test suite (subset).
    assert_eq!(md5_hex(b""), "d41d8cd98f00b204e9800998ecf8427e");
    assert_eq!(md5_hex(b"abc"), "900150983cd24fb0d6963f7d28e17f72");
    assert_eq!(md5_hex(b"message digest"), "f96b697d7cb7938d525a2f31aaf161d0");
}

/// Same as [`cram_bam_header`], but with `@HD SO:queryname` and an explicit
/// `M5` tag on the `@SQ` record (SAM-spec `md5(uppercase(ref_seq))`). Both
/// are needed for the stdin-CRAM regression test below, for two independent
/// reasons:
///
/// - `SO:queryname`: `no-alignment`'s stdin arm rejects a coordinate/unsorted
///   (no-`@HD`) input outright, before ever reading a record (a pipe cannot
///   seek for the 2-pass name-map) -- so reaching actual CRAM record decode
///   requires an explicitly queryname-sorted header.
/// - `M5`: without an `M5` tag already present on `@SQ`, htslib's CRAM
///   *writer* cannot compute one at header-flush time either (rust-htslib's
///   `Writer::from_path` flushes the header at construction, before this
///   test's code gets a chance to call `set_reference`), so it silently
///   falls back to `embed_ref=2` (embedding the full reference in the CRAM
///   itself) -- which trivially decodes with NO external reference at all,
///   defeating the whole point of this regression test. A pre-supplied `M5`
///   skips that fallback, producing a genuinely externally-referenced CRAM.
///   `M5` is ALSO the field htslib's decode-side `cram_populate_ref` needs
///   in the header to even consider the `REF_PATH`/EBI network fallback in
///   the first place (no `M5` or `UR` tag -> immediate local error, never
///   reaching `REF_PATH` at all) -- so `M5` is what makes this fixture an
///   honest stand-in for a real-world externally-referenced production CRAM.
fn cram_bam_header_queryname_with_m5(ref_len: usize, ref_seq: &[u8]) -> Header {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "queryname");
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", i64::try_from(ref_len).unwrap());
    sq.push_tag(b"M5", md5_hex(ref_seq));
    header.push_record(&sq);
    header
}

/// Writes two mapped, single-end (unpaired) 100bp reads at reference
/// positions 0 and 100 to `path` in `format`, whose bases equal the
/// reference slice they align to (so CRAM's reference-vs-read diff encoding
/// is trivial and round-trips exactly). When `format` is `Format::Cram`,
/// `ref_fa` must be `Some` (the writer's own decode reference, distinct from
/// `-r`, which the READER uses); ignored for `Format::Bam`.
fn write_cram_bam_reads(
    path: &Path,
    header: &Header,
    format: bam::Format,
    ref_seq: &[u8],
    ref_fa: Option<&Path>,
) {
    let mut writer = Writer::from_path(path, header, format).expect("bam/cram writer");
    if let Some(ref_fa) = ref_fa {
        writer.set_reference(ref_fa).expect("set CRAM writer reference");
    }
    let qual = vec![30u8; 100];
    for (i, pos) in [0i64, 100].into_iter().enumerate() {
        let start = usize::try_from(pos).unwrap();
        let seq = &ref_seq[start..start + 100];
        let mut r = bam::Record::new();
        r.set(
            format!("read{i}").as_bytes(),
            Some(&CigarString(vec![Cigar::Match(100)])),
            seq,
            &qual,
        );
        r.set_tid(0);
        r.set_pos(pos);
        r.set_flags(0); // mapped, unpaired -- general_info's frag_stdev == 0 -> single-end output.
        writer.write(&r).unwrap();
    }
    drop(writer);
}

/// The single-end candidate output path (`{prefix}.fq`), matching
/// `FastqFileSink::create`'s non-paired naming -- our fixture reads are
/// unpaired (see [`write_cram_bam_reads`]), unlike this file's other CLI
/// e2e tests, which build paired fixtures.
fn single_end_out(prefix: &Path) -> PathBuf {
    PathBuf::from(format!("{}.fq", prefix.to_str().unwrap()))
}

#[test]
fn cram_no_alignment_matches_bam_byte_for_byte() {
    let tmp = tempfile::tempdir().unwrap();
    let (ref_fa, ref_seq) = write_ref_and_fai(tmp.path(), 300);
    let header = cram_bam_header(300);

    let bam_path = tmp.path().join("reads.bam");
    write_cram_bam_reads(&bam_path, &header, bam::Format::Bam, &ref_seq, None);

    let cram_path = tmp.path().join("reads.cram");
    write_cram_bam_reads(&cram_path, &header, bam::Format::Cram, &ref_seq, Some(&ref_fa));

    let bam_prefix = tmp.path().join("bam_out");
    run_extract_ok(&[
        "-f",
        s(&ref_fa),
        "-i",
        s(&bam_path),
        "--bam-mode",
        "no-alignment",
        "-o",
        s(&bam_prefix),
    ]);

    let cram_prefix = tmp.path().join("cram_out");
    run_extract_ok(&[
        "-f",
        s(&ref_fa),
        "-r",
        s(&ref_fa),
        "-i",
        s(&cram_path),
        "--bam-mode",
        "no-alignment",
        "-o",
        s(&cram_prefix),
    ]);

    let bam_out = std::fs::read(single_end_out(&bam_prefix)).unwrap();
    let cram_out = std::fs::read(single_end_out(&cram_prefix)).unwrap();
    assert!(
        !bam_out.is_empty(),
        "the fixture reads must clear the k-mer filter and emit candidates"
    );
    assert_eq!(bam_out, cram_out, "CRAM no-alignment output must be byte-identical to BAM's");
}

#[test]
fn cram_without_reference_errors_naming_cram_and_r() {
    let tmp = tempfile::tempdir().unwrap();
    let (ref_fa, ref_seq) = write_ref_and_fai(tmp.path(), 300);
    let header = cram_bam_header(300);
    let cram_path = tmp.path().join("reads.cram");
    write_cram_bam_reads(&cram_path, &header, bam::Format::Cram, &ref_seq, Some(&ref_fa));

    let out = run_extract_raw(&[
        "-b",
        s(&cram_path),
        "--bam-mode",
        "no-alignment",
        "-f",
        s(&ref_fa),
        "-o",
        s(&tmp.path().join("o")),
    ]);
    assert!(!out.status.success(), "CRAM without -r must fail");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("CRAM") && err.contains("-r"), "error must name CRAM and -r: {err}");
}

/// Regression test for the CRITICAL no-network finding: `unum extract -i -`
/// (stdin) cannot content-sniff BAM vs CRAM (a pipe can't be re-read), so a
/// CRAM arriving on stdin without `-r` used to fall through to htslib's
/// default CRAM reference chain, which defaults to a LIVE outbound fetch
/// against the EBI `REF_PATH` endpoint (`https://www.ebi.ac.uk/ena/cram/md5/...`)
/// -- a hard violation of this project's no-outbound-network rule. `main`'s
/// `neutralize_cram_ref_path_network_fallback` now points `REF_PATH` at an
/// inert local sentinel before `Cli::parse()` ever runs, so the SAME
/// `-r`-less stdin CRAM must instead fail LOCALLY (no reference resolves
/// against the sentinel), promptly, with no successful output.
///
/// The fixture header is queryname-sorted with an explicit `M5` tag (see
/// [`cram_bam_header_queryname_with_m5`]): `no-alignment`'s stdin arm rejects
/// a coordinate/unsorted stdin input before reading any record, which would
/// trivially "pass" this test for the wrong reason, and an absent `M5` would
/// (a) make the fixture writer silently embed the reference (defeating the
/// test) and (b) short-circuit htslib's decode-side reference resolution
/// before `REF_PATH` is even consulted (also defeating the test). With both
/// present, this fixture is a genuinely externally-referenced CRAM whose
/// decode reaches the exact htslib code path the CRITICAL finding flagged.
#[test]
fn stdin_cram_without_reference_fails_locally_not_over_network() {
    let tmp = tempfile::tempdir().unwrap();
    let (ref_fa, ref_seq) = write_ref_and_fai(tmp.path(), 300);
    let header = cram_bam_header_queryname_with_m5(300, &ref_seq);
    let cram_path = tmp.path().join("reads.cram");
    write_cram_bam_reads(&cram_path, &header, bam::Format::Cram, &ref_seq, Some(&ref_fa));

    let out_prefix = tmp.path().join("o");
    let cram_file = std::fs::File::open(&cram_path).expect("open fixture CRAM for stdin redirect");

    let start = Instant::now();
    let output = Command::new(env!("CARGO_BIN_EXE_unum"))
        .arg("extract")
        .args(["-i", "-", "--bam-mode", "no-alignment", "-f", s(&ref_fa), "-o", s(&out_prefix)])
        .stdin(Stdio::from(cram_file))
        .output()
        .expect("run unum extract with stdin redirected from the fixture CRAM");
    let elapsed = start.elapsed();

    // (c) completes promptly -- does not hang on a network round-trip
    // (DNS/connect/TLS to a real host would take far longer than a local
    // reference-resolution failure).
    assert!(
        elapsed < Duration::from_secs(20),
        "must fail locally/promptly rather than hang on an outbound network call: {elapsed:?}"
    );
    // (a) exits non-zero.
    assert!(
        !output.status.success(),
        "a -r-less stdin CRAM must fail, not silently succeed via a network reference fetch"
    );
    // (b) does not emit a successful extraction.
    assert!(
        !single_end_out(&out_prefix).exists(),
        "a failed run must not emit a successful extraction output file"
    );
}

// -- `run -b` fused BAM/CRAM extraction (Stage 2c Task 5) -------------------
//
// Proves `unum run -b <bam> --bam-mode no-alignment ...` -- the fused
// extract -> genotype -> analyze native path, with candidates handed to the
// genotyper entirely in memory (no `{prefix}_candidate_*.fq` on disk) --
// produces the SAME `_genotype.tsv` as running the two-step CLI path
// (`unum extract -b <bam> --bam-mode no-alignment` followed by a separate
// `unum genotype` on the emitted candidate FASTQs). `--bam-mode alignment`'s
// fused coordinate-alignment route is exercised indirectly by this same
// `build_coordinate_bam` fixture in the `extract`-only tests above; this
// file's job is proving the `run` in-memory hand-off itself, not
// re-proving extraction correctness (already golden-gated elsewhere).

/// Runs `unum run` with `args`, asserting a clean exit. Mirrors
/// [`run_extract_ok`].
fn run_run_ok(args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_unum")).arg("run").args(args).status().unwrap();
    assert!(status.success(), "`unum run` exited non-zero for args {args:?}");
}

/// Runs `unum genotype` with `args`, asserting a clean exit. Mirrors
/// [`run_extract_ok`].
fn run_genotype_ok(args: &[&str]) {
    let status =
        Command::new(env!("CARGO_BIN_EXE_unum")).arg("genotype").args(args).status().unwrap();
    assert!(status.success(), "`unum genotype` exited non-zero for args {args:?}");
}

/// Writes a single-sequence reference-sequence FASTA (`>KIR2DL1*0010101\n{seq}\n`, no header
/// comment) containing exactly [`kir2dl1_sequence`] -- the SAME allele [`build_coordinate_bam`]'s
/// on-target pair is a substring of. Used as `-f` for BOTH `--bam-mode no-alignment` extraction
/// (which only needs id+sequence, same as [`coord_fasta_path`]'s multi-gene coord FASTA already
/// proven to work there) AND `genotype`'s own reference load, which does NOT tolerate the coord
/// FASTA's `chrom start end strand` header comment: `genotype::load_reference` feeds that comment
/// to `AlleleRef::new`/`parse_exon_comment` as an exon-interval list, and `chrom start end`'s
/// digits parse as a nonsensical exon range (`start`..`end`, both far outside the sequence's own
/// length) that zeroes out every position's exon membership and panics downstream
/// (`get_seq_missing_base_coverage: allele has zero exon positions`). Omitting the comment
/// entirely sidesteps this: `AlleleRef::new(seq, None)` falls back to a single whole-sequence
/// exon, matching every OTHER `genotype`-mode reference fixture in this workspace (`fixtures/
/// example/ref/kir_rna_seq.fa` uses a real base-offset exon comment instead, which this minimal
/// fixture does not need since its one on-target pair is an exact, diff-free substring).
fn write_run_no_alignment_ref_fasta(dir: &Path) -> PathBuf {
    let path = dir.join("run_ref.fa");
    std::fs::write(&path, format!(">KIR2DL1*0010101\n{}\n", kir2dl1_sequence())).unwrap();
    path
}

#[test]
fn run_bam_no_alignment_matches_extract_then_genotype() {
    let tmp_dir = tempfile::tempdir().unwrap();
    let tmp = tmp_dir.path();
    let ref_fa = write_run_no_alignment_ref_fasta(tmp);
    // 5 duplicate on-target pairs -- see `build_coordinate_bam_replicated`'s doc comment for
    // why a single pair's abundance falls under the genotyper's default `filter_cov`.
    let bam = build_coordinate_bam_replicated(tmp, 5);

    // Fused: `unum run -b <bam> --bam-mode no-alignment -f <ref_fa>` -- extract, genotype,
    // and analyze all in one process.
    let run_prefix = tmp.join("via_run");
    run_run_ok(&[
        "-b",
        s(&bam),
        "--bam-mode",
        "no-alignment",
        "-f",
        s(&ref_fa),
        "-o",
        s(&run_prefix),
    ]);

    // Separately: `unum extract -b <bam> --bam-mode no-alignment` writes candidate FASTQs,
    // then `unum genotype` reads them back from disk -- the pre-Task-5 two-process path.
    let extract_prefix = tmp.join("via_extract");
    run_extract_ok(&[
        "-b",
        s(&bam),
        "--bam-mode",
        "no-alignment",
        "-f",
        s(&ref_fa),
        "-o",
        s(&extract_prefix),
    ]);
    let genotype_prefix = tmp.join("via_genotype");
    run_genotype_ok(&[
        "-f",
        s(&ref_fa),
        "-1",
        s(&paired_out_1(&extract_prefix)),
        "-2",
        s(&paired_out_2(&extract_prefix)),
        "-o",
        s(&genotype_prefix),
    ]);

    let run_genotype_tsv =
        std::fs::read_to_string(format!("{}_genotype.tsv", run_prefix.display())).unwrap();
    let separate_genotype_tsv =
        std::fs::read_to_string(format!("{}_genotype.tsv", genotype_prefix.display())).unwrap();
    assert_eq!(
        run_genotype_tsv, separate_genotype_tsv,
        "`run -b --bam-mode no-alignment`'s fused _genotype.tsv must match the separate \
         extract-then-genotype pipeline byte-for-byte"
    );

    // Teeth: a real genotype call must actually have happened against this fixture, not an
    // empty/no-call row (`Genotyper::get_allele_description`'s no-call sentinel is a literal
    // `.`, not an empty string -- see `fixtures/example/oracle_genotype.golden.tsv`'s
    // zero-`calledCnt` rows). `allele1` reports the "major allele" name (allele-nomenclature
    // group resolution, e.g. `KIR2DL1*001` for `KIR2DL1*0010101` -- `unum`/T1K report at this
    // resolution whenever finer-digit alleles aren't distinguished), so check for that prefix
    // rather than the full allele name the fixture's on-target pair was built from.
    let kir2dl1_line = run_genotype_tsv
        .lines()
        .find(|line| line.starts_with("KIR2DL1\t"))
        .expect("run_genotype_tsv must contain a KIR2DL1 row");
    assert!(
        kir2dl1_line.contains("KIR2DL1*001"),
        "expected a called KIR2DL1*001 allele, got: {kir2dl1_line}"
    );
    assert!(
        !kir2dl1_line.starts_with("KIR2DL1\t0\t"),
        "expected calledAlleleCnt > 0 (a real call), got: {kir2dl1_line}"
    );
}

/// Same parity check as [`run_bam_no_alignment_matches_extract_then_genotype`], but for the
/// OTHER Task 5 in-scope route: `run -b --bam-mode alignment` on a coordinate-sorted BAM. Unlike
/// `no-alignment`, `alignment` reads its gene intervals + k-mer seed reference from `-c` (the
/// coord FASTA) rather than `-f`; `-f` is still required by `RunArgs` (unlike `ExtractArgs`,
/// where it is optional and rejected in `alignment` mode) because [`genotype_args_for`] always
/// builds the genotyper's own reference from `args.ref_seq_fasta` regardless of `--bam-mode`, so
/// this test supplies both.
#[test]
fn run_bam_alignment_matches_extract_then_genotype() {
    let tmp_dir = tempfile::tempdir().unwrap();
    let tmp = tmp_dir.path();
    let coord_fa = coord_fasta_path();
    let ref_fa = write_run_no_alignment_ref_fasta(tmp);
    let bam = build_coordinate_bam_replicated(tmp, 5);

    // Fused: `unum run -b <bam> --bam-mode alignment -c <coord_fa> -f <ref_fa>`.
    let run_prefix = tmp.join("via_run");
    run_run_ok(&[
        "-b",
        s(&bam),
        "--bam-mode",
        "alignment",
        "-c",
        s(&coord_fa),
        "-f",
        s(&ref_fa),
        "-o",
        s(&run_prefix),
    ]);

    // Separately: `unum extract -b <bam> --bam-mode alignment -c <coord_fa>` then `unum
    // genotype -f <ref_fa>` on the emitted candidate FASTQs.
    let extract_prefix = tmp.join("via_extract");
    run_extract_ok(&[
        "-b",
        s(&bam),
        "--bam-mode",
        "alignment",
        "-c",
        s(&coord_fa),
        "-o",
        s(&extract_prefix),
    ]);
    let genotype_prefix = tmp.join("via_genotype");
    run_genotype_ok(&[
        "-f",
        s(&ref_fa),
        "-1",
        s(&paired_out_1(&extract_prefix)),
        "-2",
        s(&paired_out_2(&extract_prefix)),
        "-o",
        s(&genotype_prefix),
    ]);

    let run_genotype_tsv =
        std::fs::read_to_string(format!("{}_genotype.tsv", run_prefix.display())).unwrap();
    let separate_genotype_tsv =
        std::fs::read_to_string(format!("{}_genotype.tsv", genotype_prefix.display())).unwrap();
    assert_eq!(
        run_genotype_tsv, separate_genotype_tsv,
        "`run -b --bam-mode alignment`'s fused _genotype.tsv must match the separate \
         extract-then-genotype pipeline byte-for-byte"
    );

    let kir2dl1_line = run_genotype_tsv
        .lines()
        .find(|line| line.starts_with("KIR2DL1\t"))
        .expect("run_genotype_tsv must contain a KIR2DL1 row");
    assert!(
        kir2dl1_line.contains("KIR2DL1*001"),
        "expected a called KIR2DL1*001 allele, got: {kir2dl1_line}"
    );
    assert!(
        !kir2dl1_line.starts_with("KIR2DL1\t0\t"),
        "expected calledAlleleCnt > 0 (a real call), got: {kir2dl1_line}"
    );
}

/// `run -b --bam-mode alignment` on a grouped/name-sorted BAM is explicitly OUT of Task 5's
/// scope (the grouped-alignment one-pass extractor lands in a later task) -- proves it fails
/// with a clear deferral message rather than silently misrouting or panicking, mirroring
/// `alignment_on_name_sorted_bam_errors_with_sort_hint`'s coverage of the same guard on
/// `extract`.
#[test]
fn run_bam_alignment_on_grouped_bam_lands_in_later_task() {
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let ref_fa = write_run_no_alignment_ref_fasta(tmp.path());
    let bam = build_name_sorted_bam(tmp.path());

    let out = Command::new(env!("CARGO_BIN_EXE_unum"))
        .arg("run")
        .args([
            "-b",
            s(&bam),
            "--bam-mode",
            "alignment",
            "-c",
            s(&coord_fa),
            "-f",
            s(&ref_fa),
            "-o",
            s(&tmp.path().join("o")),
        ])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "run -b --bam-mode alignment on a grouped/name-sorted BAM must fail (deferred to a \
         later task)"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("later task") || stderr.contains("grouped-alignment"),
        "expected a later-task deferral message, got: {stderr}"
    );
}

/// `run -b` combined with `-1`/`-2`/`-u` must be rejected outright, mirroring
/// `extract::resolve_extract_input`'s own `-b` vs. `-1`/`-2`/`-u` mutual-exclusion check --
/// without this guard, `run`'s early `-b` branch would silently win and discard the FASTQ flags.
#[test]
fn run_bam_flag_with_fastq_flags_is_rejected() {
    let tmp = tempfile::tempdir().unwrap();
    let ref_fa = write_run_no_alignment_ref_fasta(tmp.path());
    let bam = build_coordinate_bam(tmp.path());

    let out = Command::new(env!("CARGO_BIN_EXE_unum"))
        .arg("run")
        .args([
            "-b",
            s(&bam),
            "--bam-mode",
            "no-alignment",
            "-f",
            s(&ref_fa),
            "-u",
            "reads.fq",
            "-o",
            s(&tmp.path().join("o")),
        ])
        .output()
        .unwrap();
    assert!(!out.status.success(), "-b combined with -u must fail, not silently win");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("mutually exclusive"),
        "expected a mutual-exclusion error, got: {stderr}"
    );
}

/// `-c` (coord FASTA) applies only to `--bam-mode alignment`; `run -b --bam-mode no-alignment`
/// must reject it, mirroring `extract::require_no_coord_fasta`'s guard on `extract`.
#[test]
fn run_bam_no_alignment_with_coord_fasta_is_rejected() {
    let tmp = tempfile::tempdir().unwrap();
    let ref_fa = write_run_no_alignment_ref_fasta(tmp.path());
    let coord_fa = coord_fasta_path();
    let bam = build_coordinate_bam(tmp.path());

    let out = Command::new(env!("CARGO_BIN_EXE_unum"))
        .arg("run")
        .args([
            "-b",
            s(&bam),
            "--bam-mode",
            "no-alignment",
            "-f",
            s(&ref_fa),
            "-c",
            s(&coord_fa),
            "-o",
            s(&tmp.path().join("o")),
        ])
        .output()
        .unwrap();
    assert!(!out.status.success(), "-c under --bam-mode no-alignment must fail, not be ignored");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("-c"), "expected a -c-specific error, got: {stderr}");
}

// -- `run -b` `RunArgs`-specific guard coverage (Stage 2c Task 5 review fix) -
//
// The four tests below cover `run_bam_fused`'s own error-path guards
// (`require_reference_for_cram_run`/`require_coord_fasta_run`/the
// `--bam-mode` presence check), which had no e2e coverage: Task 5's own
// `run_bam_*` tests above only exercise the success paths plus the
// grouped-alignment deferral and the `-b`/`-c` mutual-exclusion guards.
// Mirrors the sibling `extract`-side guard tests
// (`b_flag_without_bam_mode_errors`, `cram_without_reference_errors_naming_cram_and_r`)
// and reuses this file's BAM/CRAM fixture helpers.

/// `run -b` without `--bam-mode` must fail, naming `--bam-mode`, mirroring `extract`'s own
/// `b_flag_without_bam_mode_errors` guard test.
#[test]
fn run_bam_missing_bam_mode_errors() {
    let tmp = tempfile::tempdir().unwrap();
    let ref_fa = write_run_no_alignment_ref_fasta(tmp.path());
    let bam = build_coordinate_bam(tmp.path());

    let out = Command::new(env!("CARGO_BIN_EXE_unum"))
        .arg("run")
        .args(["-b", s(&bam), "-f", s(&ref_fa), "-o", s(&tmp.path().join("o"))])
        .output()
        .unwrap();
    assert!(!out.status.success(), "run -b without --bam-mode must fail");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("--bam-mode"), "expected a --bam-mode-specific error, got: {stderr}");
}

/// `run -b --bam-mode alignment` without `-c` must fail, naming `-c` (the gene coordinate
/// FASTA), mirroring `require_coord_fasta_run`'s doc-commented error and
/// `run_bam_no_alignment_with_coord_fasta_is_rejected`'s coverage of the sibling `-c` guard.
#[test]
fn run_bam_alignment_missing_coord_fasta_is_rejected() {
    let tmp = tempfile::tempdir().unwrap();
    let ref_fa = write_run_no_alignment_ref_fasta(tmp.path());
    let bam = build_coordinate_bam(tmp.path());

    let out = Command::new(env!("CARGO_BIN_EXE_unum"))
        .arg("run")
        .args([
            "-b",
            s(&bam),
            "--bam-mode",
            "alignment",
            "-f",
            s(&ref_fa),
            "-o",
            s(&tmp.path().join("o")),
        ])
        .output()
        .unwrap();
    assert!(!out.status.success(), "run -b --bam-mode alignment without -c must fail");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("-c"), "expected a -c-specific error, got: {stderr}");
}

/// `run -b <cram> --bam-mode no-alignment` without `-r` must fail, naming CRAM and `-r`,
/// mirroring `extract`'s own `cram_without_reference_errors_naming_cram_and_r` guard test.
#[test]
fn run_bam_cram_no_alignment_without_reference_is_rejected() {
    let tmp = tempfile::tempdir().unwrap();
    let (ref_fa, ref_seq) = write_ref_and_fai(tmp.path(), 300);
    let header = cram_bam_header(300);
    let cram_path = tmp.path().join("reads.cram");
    write_cram_bam_reads(&cram_path, &header, bam::Format::Cram, &ref_seq, Some(&ref_fa));

    let out = Command::new(env!("CARGO_BIN_EXE_unum"))
        .arg("run")
        .args([
            "-b",
            s(&cram_path),
            "--bam-mode",
            "no-alignment",
            "-f",
            s(&ref_fa),
            "-o",
            s(&tmp.path().join("o")),
        ])
        .output()
        .unwrap();
    assert!(!out.status.success(), "run -b <cram> without -r must fail");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("CRAM") && stderr.contains("-r"),
        "error must name CRAM and -r: {stderr}"
    );
}

/// `-r` with a BAM (not CRAM) is accepted-and-ignored:
/// `require_reference_for_cram_run` "returns `None` unconditionally" whenever `is_cram` is
/// false, even if `-r` was harmlessly also passed -- so `run -b <bam> --bam-mode no-alignment
/// -r <ref>` must NOT error on `-r` and must proceed through the full fused
/// extract -> genotype -> analyze pipeline, exactly as
/// `run_bam_no_alignment_matches_extract_then_genotype` does without `-r`.
#[test]
fn run_bam_no_alignment_reference_with_bam_is_accepted_and_ignored() {
    let tmp_dir = tempfile::tempdir().unwrap();
    let tmp = tmp_dir.path();
    let ref_fa = write_run_no_alignment_ref_fasta(tmp);
    let bam = build_coordinate_bam_replicated(tmp, 5);
    // A plausible-looking CRAM-decode-reference `-r` value (real file + `.fai` sibling) that
    // must be ignored entirely for BAM input -- irrelevant to the BAM's own contents.
    let (unused_reference, _ref_seq) = write_ref_and_fai(tmp, 300);

    let out_prefix = tmp.join("o");
    run_run_ok(&[
        "-b",
        s(&bam),
        "--bam-mode",
        "no-alignment",
        "-f",
        s(&ref_fa),
        "-r",
        s(&unused_reference),
        "-o",
        s(&out_prefix),
    ]);

    let genotype_tsv =
        std::fs::read_to_string(format!("{}_genotype.tsv", out_prefix.display())).unwrap();
    let kir2dl1_line = genotype_tsv
        .lines()
        .find(|line| line.starts_with("KIR2DL1\t"))
        .expect("genotype_tsv must contain a KIR2DL1 row");
    assert!(
        kir2dl1_line.contains("KIR2DL1*001"),
        "expected a called KIR2DL1*001 allele (proving the run proceeded past extraction \
         despite -r), got: {kir2dl1_line}"
    );
    assert!(
        !kir2dl1_line.starts_with("KIR2DL1\t0\t"),
        "expected calledAlleleCnt > 0 (a real call), got: {kir2dl1_line}"
    );
}
