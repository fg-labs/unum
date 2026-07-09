//! Binary-level end-to-end tests for `unum extract --bam-mode`: proves that
//! `-i x.bam --bam-mode alignment` and `-b x.bam --bam-mode alignment` reach
//! the SAME coordinate-`alignment` extraction path, and that every
//! reserved/guard combination the Stage 2a dispatcher is supposed to reject
//! (missing `--bam-mode`, `no-alignment`, a non-coordinate-sorted BAM) errors
//! with a helpful, on-topic message.
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
        "-f",
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
        "-f",
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
fn no_alignment_mode_errors_as_reserved() {
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let bam = build_coordinate_bam(tmp.path());
    let out = run_extract_raw(&[
        "-f",
        s(&coord_fa),
        "-i",
        s(&bam),
        "--bam-mode",
        "no-alignment",
        "-o",
        s(&tmp.path().join("o")),
    ]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("later release"));
}

#[test]
fn alignment_on_name_sorted_bam_errors_with_sort_hint() {
    let tmp = tempfile::tempdir().unwrap();
    let coord_fa = coord_fasta_path();
    let bam = build_name_sorted_bam(tmp.path());
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
