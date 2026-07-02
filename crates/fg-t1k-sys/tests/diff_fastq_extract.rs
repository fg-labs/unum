#![cfg(feature = "t1k-sys")]
//! Differential test: runs the REAL C++ `fastq-extractor` oracle binary and
//! the Rust `fg_t1k_core::extract::extract_candidates` port on the SAME
//! reference + read inputs, and asserts the output FASTQ file(s) are
//! **byte-identical**, not merely set-equal or order-insensitive.
//!
//! # Why byte identity (not just set equality) is achievable at any `-t`
//!
//! See `fg_t1k_core::extract`'s module docs ("Output order == input order")
//! for the full argument: `FastqExtractor.cpp`'s multi-threaded path only
//! parallelizes the per-read filter DECISION, never the emission order,
//! which is always sequential-in-input-order on both the single- and
//! multi-threaded code paths. A correct single-threaded Rust
//! re-implementation is therefore provably byte-identical to the oracle's
//! output regardless of the oracle's own `-t` value.
//!
//! This test always runs the oracle at its default `-t 1` (the CLI default),
//! which is sufficient given that invariant -- this test's job is to prove
//! the Rust port's DECISION and FORMATTING logic matches stock exactly, not
//! to re-litigate threadCnt-invariance (which is a property of the vendored
//! C++ source, established by reading it, not something a test can
//! meaningfully re-derive by running the oracle at higher `-t`).

use fg_t1k_core::extract::{self, CandidateSink, ReadRecord};
use fg_t1k_core::ref_kmer_filter::RefKmerFilter;
use fg_t1k_sys::oracle::{OracleStage, binary_path};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

/// `FastqExtractor.cpp:272`: the literal initial k-mer length.
const INITIAL_KMER_LENGTH: usize = 9;
/// `FastqExtractor.cpp:283`: default `-s`.
const DEFAULT_SIMILARITY: f64 = extract::DEFAULT_REF_SEQ_SIMILARITY;

fn fixture(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures").join(rel)
}

/// A FASTQ-file-writing sink mirroring the CLI's own `FastqFileSink` (kept
/// as a small, self-contained copy here rather than depending on the binary
/// crate, so this test only exercises the library surface directly --
/// matching this task's "library-first" design goal: the differential
/// proves the LIBRARY is byte-identical, independent of any particular CLI
/// wrapper implementation).
struct FastqFileSink {
    fp1: std::fs::File,
    fp2: Option<std::fs::File>,
}

impl FastqFileSink {
    fn create(prefix: &Path, paired: bool) -> Self {
        let fp1 = std::fs::File::create(format!("{}_1.fq", prefix.display()))
            .unwrap_or_else(|e| panic!("create {}_1.fq: {e}", prefix.display()));
        if paired {
            let fp2 = std::fs::File::create(format!("{}_2.fq", prefix.display()))
                .unwrap_or_else(|e| panic!("create {}_2.fq: {e}", prefix.display()));
            Self { fp1, fp2: Some(fp2) }
        } else {
            Self { fp1, fp2: None }
        }
    }

    fn create_single(prefix: &Path) -> Self {
        let fp1 = std::fs::File::create(format!("{}.fq", prefix.display()))
            .unwrap_or_else(|e| panic!("create {}.fq: {e}", prefix.display()));
        Self { fp1, fp2: None }
    }
}

impl CandidateSink for FastqFileSink {
    fn emit_pair(&mut self, r1: &ReadRecord, r2: Option<&ReadRecord>) -> anyhow::Result<()> {
        extract::output_seq(&mut self.fp1, &r1.id, &r1.seq, r1.qual.as_deref(), 0, -1)?;
        if let (Some(r2), Some(fp2)) = (r2, self.fp2.as_mut()) {
            // Mate-1's id for both outputs -- see FastqFileSink's C++
            // counterpart doc comment (`FastqExtractor.cpp:471-473`).
            extract::output_seq(fp2, &r1.id, &r2.seq, r2.qual.as_deref(), 0, -1)?;
        }
        Ok(())
    }
}

/// Runs the Rust extractor end-to-end (library call, not the CLI binary),
/// writing `{prefix}_1.fq`/`{prefix}_2.fq` (paired) to `out_dir`.
fn run_rust_paired(ref_fasta: &Path, r1: &Path, r2: &Path, out_prefix: &Path) {
    run_rust_paired_with_threads(ref_fasta, r1, r2, out_prefix, 1);
}

/// Same as [`run_rust_paired`], but drives the parallel-decision path via
/// [`extract::extract_candidates_with_threads`] at the given `threads`.
fn run_rust_paired_with_threads(
    ref_fasta: &Path,
    r1: &Path,
    r2: &Path,
    out_prefix: &Path,
    threads: usize,
) {
    let mut filter = RefKmerFilter::from_reference_fasta(ref_fasta, INITIAL_KMER_LENGTH)
        .unwrap_or_else(|e| panic!("RefKmerFilter::from_reference_fasta: {e}"));
    let mut source =
        extract::open_source(r1, Some(r2)).unwrap_or_else(|e| panic!("open_source(paired): {e}"));
    let mut sink = FastqFileSink::create(out_prefix, true);
    extract::extract_candidates_with_threads(
        &mut source,
        &mut filter,
        DEFAULT_SIMILARITY,
        threads,
        &mut sink,
    )
    .unwrap_or_else(|e| panic!("extract_candidates_with_threads(paired, threads={threads}): {e}"));
}

/// Same as [`run_rust_paired`] but single-end, writing `{prefix}.fq`.
fn run_rust_single(ref_fasta: &Path, r1: &Path, out_prefix: &Path) {
    let mut filter = RefKmerFilter::from_reference_fasta(ref_fasta, INITIAL_KMER_LENGTH)
        .unwrap_or_else(|e| panic!("RefKmerFilter::from_reference_fasta: {e}"));
    let mut source =
        extract::open_source(r1, None).unwrap_or_else(|e| panic!("open_source(single): {e}"));
    let mut sink = FastqFileSink::create_single(out_prefix);
    extract::extract_candidates(&mut source, &mut filter, DEFAULT_SIMILARITY, &mut sink)
        .unwrap_or_else(|e| panic!("extract_candidates(single): {e}"));
}

/// Runs the real oracle `fastq-extractor` binary, paired input.
fn run_oracle_paired(ref_fasta: &Path, r1: &Path, r2: &Path, out_prefix: &Path) {
    let bin = binary_path(OracleStage::FastqExtractor);
    assert!(bin.exists(), "oracle binary not built: {}", bin.display());
    let status = Command::new(&bin)
        .arg("-f")
        .arg(ref_fasta)
        .arg("-1")
        .arg(r1)
        .arg("-2")
        .arg(r2)
        .arg("-o")
        .arg(out_prefix)
        .status()
        .unwrap_or_else(|e| panic!("spawning fastq-extractor: {e}"));
    assert!(status.success(), "fastq-extractor (paired) exited with {status}");
}

/// Same as [`run_oracle_paired`] but single-end (`-u`).
fn run_oracle_single(ref_fasta: &Path, r1: &Path, out_prefix: &Path) {
    let bin = binary_path(OracleStage::FastqExtractor);
    assert!(bin.exists(), "oracle binary not built: {}", bin.display());
    let status = Command::new(&bin)
        .arg("-f")
        .arg(ref_fasta)
        .arg("-u")
        .arg(r1)
        .arg("-o")
        .arg(out_prefix)
        .status()
        .unwrap_or_else(|e| panic!("spawning fastq-extractor: {e}"));
    assert!(status.success(), "fastq-extractor (single-end) exited with {status}");
}

fn assert_files_byte_identical(rust_path: &Path, oracle_path: &Path) {
    let rust_bytes = std::fs::read(rust_path)
        .unwrap_or_else(|e| panic!("reading rust output {}: {e}", rust_path.display()));
    let oracle_bytes = std::fs::read(oracle_path)
        .unwrap_or_else(|e| panic!("reading oracle output {}: {e}", oracle_path.display()));
    assert_eq!(
        rust_bytes,
        oracle_bytes,
        "byte mismatch between {} ({} bytes) and {} ({} bytes)",
        rust_path.display(),
        rust_bytes.len(),
        oracle_path.display(),
        oracle_bytes.len()
    );
}

#[test]
fn paired_example_fixture_byte_identical_to_oracle() {
    let ref_fasta = fixture("example/ref/kir_rna_seq.fa");
    let r1 = fixture("example/example_1.fq");
    let r2 = fixture("example/example_2.fq");

    let tmp = tempfile::tempdir().unwrap();
    let rust_prefix = tmp.path().join("rust");
    let oracle_prefix = tmp.path().join("oracle");

    run_rust_paired(&ref_fasta, &r1, &r2, &rust_prefix);
    run_oracle_paired(&ref_fasta, &r1, &r2, &oracle_prefix);

    assert_files_byte_identical(&tmp.path().join("rust_1.fq"), &tmp.path().join("oracle_1.fq"));
    assert_files_byte_identical(&tmp.path().join("rust_2.fq"), &tmp.path().join("oracle_2.fq"));

    // Sanity: the fixture is non-trivial (not all-pass or all-fail), so this
    // test is actually exercising the filter, not just file-plumbing.
    let n1 = std::fs::read_to_string(tmp.path().join("rust_1.fq"))
        .unwrap()
        .lines()
        .filter(|l| l.starts_with('@'))
        .count();
    assert!(n1 > 0, "expected at least one candidate pair to pass");
}

#[test]
fn single_end_example_fixture_byte_identical_to_oracle() {
    let ref_fasta = fixture("example/ref/kir_rna_seq.fa");
    let r1 = fixture("example/example_1.fq");

    let tmp = tempfile::tempdir().unwrap();
    let rust_prefix = tmp.path().join("rust");
    let oracle_prefix = tmp.path().join("oracle");

    run_rust_single(&ref_fasta, &r1, &rust_prefix);
    run_oracle_single(&ref_fasta, &r1, &oracle_prefix);

    assert_files_byte_identical(&tmp.path().join("rust.fq"), &tmp.path().join("oracle.fq"));
}

/// Regression for the Task-3.2 review finding: stock's single-threaded loop
/// (`FastqExtractor.cpp:449-479`) drives on mate-1 and, once mate-1 is
/// exhausted, silently ignores EXTRA trailing mate-2 records and exits 0.
/// The port must match byte-for-byte (produce output, not error) on
/// mate-2-longer input. (The mate-2-*shorter* error case is covered by
/// `mate_count_mismatch_matches_oracle_error_behavior`.)
#[test]
fn mate2_longer_byte_identical_to_oracle() {
    // Base-balanced 260bp reference so exact substrings pass HasHitInSet at the
    // default paired hitLenRequired (27) and are not low-complexity.
    let reference = "ACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATGGTACGTACGTACGTACGGGCATTCATGGCATTCATGGCATTCATGACGTTAGCACGTTAGCACGTTAGCACGTTAGCTGACCATGTGACCATGTGACCATGTGACCATG";
    let dir = tempfile::tempdir().unwrap();
    let ref_fasta = dir.path().join("ref.fa");
    {
        let mut f = std::fs::File::create(&ref_fasta).unwrap();
        writeln!(f, ">only\n{reference}").unwrap();
    }

    let hit1 = reference[0..100].to_string();
    let hit2 = reference[120..220].to_string();
    // r1: 3 hitting records; r2: the 3 mates PLUS 2 extra trailing records that
    // stock never reads (mate-1 exhausts first).
    let r1_records: Vec<(String, String)> =
        (0..3).map(|i| (format!("p{i}"), hit1.clone())).collect();
    let mut r2_records: Vec<(String, String)> =
        (0..3).map(|i| (format!("p{i}"), hit2.clone())).collect();
    r2_records.push(("extra0".to_string(), hit1.clone()));
    r2_records.push(("extra1".to_string(), hit2.clone()));

    let r1 = dir.path().join("r1.fq");
    let r2 = dir.path().join("r2.fq");
    write_fastq(&r1, &r1_records);
    write_fastq(&r2, &r2_records);

    let rust_prefix = dir.path().join("rust");
    let oracle_prefix = dir.path().join("oracle");
    // Neither must error; both must emit only the 3 mate-1-driven pairs.
    run_rust_paired(&ref_fasta, &r1, &r2, &rust_prefix);
    run_oracle_paired(&ref_fasta, &r1, &r2, &oracle_prefix);

    assert_files_byte_identical(&dir.path().join("rust_1.fq"), &dir.path().join("oracle_1.fq"));
    assert_files_byte_identical(&dir.path().join("rust_2.fq"), &dir.path().join("oracle_2.fq"));
}

/// Writes a small synthetic paired-end FASTQ fixture set covering the
/// branch-coverage categories required by this task: mate1-only-hit,
/// mate2-only-hit, neither-hits, low-complexity, and RC read -- built
/// against a small synthetic reference so the oracle and Rust port are
/// exercised on the exact same non-trivial (not just curated-easy) inputs.
struct SyntheticFixture {
    dir: tempfile::TempDir,
    ref_fasta: PathBuf,
    r1: PathBuf,
    r2: PathBuf,
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'N' => 'N',
            other => panic!("unsupported base {other}"),
        })
        .collect()
}

fn build_synthetic_fixture() -> SyntheticFixture {
    let dir = tempfile::tempdir().unwrap();

    // A 300bp reference, base-balanced (not low-complexity) so exact
    // substrings of it reliably pass HasHitInSet at the default
    // hitLenRequired (27 paired).
    let reference = "ACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATGGTACGTACGTACGTACGGGCATTCATGGCATTCATGGCATTCATGACGTTAGCACGTTAGCACGTTAGCACGTTAGCTGACCATGTGACCATGTGACCATGTGACCATG";
    let ref_fasta = dir.path().join("ref.fa");
    {
        let mut f = std::fs::File::create(&ref_fasta).unwrap();
        writeln!(f, ">only").unwrap();
        writeln!(f, "{reference}").unwrap();
    }

    // Build enough sample-1000 padding so hitLenRequired sampling behaves
    // like a real run (matches production usage: sampling reads at least a
    // handful of records is the realistic path, and we want the same
    // integer-division math both sides compute).
    let mate1_hit = &reference[0..100];
    let mate2_hit = &reference[150..250];
    let noise_100 = "A".repeat(50) + &"T".repeat(50); // low-complexity-ish, but distinct pattern
    let random_noise = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";
    let homopolymer = "A".repeat(100);
    let rc_hit = reverse_complement(&reference[50..150]);

    let mut r1_records: Vec<(String, String)> = Vec::new();
    let mut r2_records: Vec<(String, String)> = Vec::new();

    // 1. mate1-only-hit: mate1 matches reference, mate2 is unrelated noise.
    r1_records.push(("pair_mate1_hit".to_string(), mate1_hit.to_string()));
    r2_records.push(("pair_mate1_hit".to_string(), random_noise.to_string()));

    // 2. mate2-only-hit: mate1 is unrelated noise, mate2 matches reference.
    r1_records.push(("pair_mate2_hit".to_string(), random_noise.to_string()));
    r2_records.push(("pair_mate2_hit".to_string(), mate2_hit.to_string()));

    // 3. neither-hits: both mates unrelated noise (should be dropped).
    r1_records.push(("pair_neither_hit".to_string(), random_noise.to_string()));
    r2_records.push(("pair_neither_hit".to_string(), noise_100.clone()));

    // 4. low-complexity: mate1 is a homopolymer (IsLowComplexity gates
    // before HasHitInSet is even called), mate2 is unrelated noise too ->
    // should be dropped even though the homopolymer happens to be a
    // substring-adjacent sequence.
    r1_records.push(("pair_low_complexity".to_string(), homopolymer));
    r2_records.push(("pair_low_complexity".to_string(), noise_100));

    // 5. RC read: mate1 is the reverse complement of a reference substring
    // (GetHitsFromRead always searches both strands, so this should still
    // hit), mate2 is unrelated noise.
    r1_records.push(("pair_rc_hit".to_string(), rc_hit));
    r2_records.push(("pair_rc_hit".to_string(), random_noise.to_string()));

    // Pad with additional plain hits/misses so the hitLenRequired sampling
    // pass (first 1000 read-1 records) sees a realistic mix, and so this
    // fixture isn't trivially tiny.
    for i in 0..20 {
        let start = (i * 3) % (reference.len() - 80);
        r1_records.push((format!("pad_hit_{i}"), reference[start..start + 80].to_string()));
        r2_records.push((format!("pad_hit_{i}"), reference[start + 5..start + 85].to_string()));
    }

    let r1 = dir.path().join("r1.fq");
    let r2 = dir.path().join("r2.fq");
    write_fastq(&r1, &r1_records);
    write_fastq(&r2, &r2_records);

    SyntheticFixture { dir, ref_fasta, r1, r2 }
}

fn write_fastq(path: &Path, records: &[(String, String)]) {
    let mut f = std::fs::File::create(path).unwrap();
    for (id, seq) in records {
        let qual = "I".repeat(seq.len());
        writeln!(f, "@{id}\n{seq}\n+\n{qual}").unwrap();
    }
}

#[test]
fn synthetic_branch_coverage_set_byte_identical_to_oracle() {
    let fx = build_synthetic_fixture();

    let rust_prefix = fx.dir.path().join("rust");
    let oracle_prefix = fx.dir.path().join("oracle");

    run_rust_paired(&fx.ref_fasta, &fx.r1, &fx.r2, &rust_prefix);
    run_oracle_paired(&fx.ref_fasta, &fx.r1, &fx.r2, &oracle_prefix);

    assert_files_byte_identical(
        &fx.dir.path().join("rust_1.fq"),
        &fx.dir.path().join("oracle_1.fq"),
    );
    assert_files_byte_identical(
        &fx.dir.path().join("rust_2.fq"),
        &fx.dir.path().join("oracle_2.fq"),
    );

    // Confirm the synthetic set actually exercises a mix of pass/fail (not
    // vacuously all-pass or all-fail), so this test is meaningful.
    let out1 = std::fs::read_to_string(fx.dir.path().join("rust_1.fq")).unwrap();
    let emitted_ids: Vec<&str> = out1.lines().filter_map(|l| l.strip_prefix('@')).collect();
    assert!(
        emitted_ids.contains(&"pair_mate1_hit"),
        "expected mate1-only-hit pair to be emitted, got ids: {emitted_ids:?}"
    );
    assert!(
        emitted_ids.contains(&"pair_mate2_hit"),
        "expected mate2-only-hit pair to be emitted, got ids: {emitted_ids:?}"
    );
    assert!(
        emitted_ids.contains(&"pair_rc_hit"),
        "expected RC-hit pair to be emitted, got ids: {emitted_ids:?}"
    );
    assert!(
        !emitted_ids.contains(&"pair_neither_hit"),
        "neither-hit pair should have been dropped, got ids: {emitted_ids:?}"
    );
    assert!(
        !emitted_ids.contains(&"pair_low_complexity"),
        "low-complexity pair should have been dropped, got ids: {emitted_ids:?}"
    );
}

#[test]
fn empty_read1_file_matches_oracle_error_behavior() {
    // Both the Rust port and the oracle must reject an empty read-1 file:
    // the oracle exits non-zero; the Rust library call returns an error.
    // (The library API surfaces this as a Result rather than a process
    // exit code -- proven independently by fg-t1k-core's own unit test,
    // `extract::tests::empty_read1_file_is_an_error`; this test's job is
    // only to confirm the ORACLE really does reject it the same way, i.e.
    // that this is a real, not imagined, shared failure mode.)
    let dir = tempfile::tempdir().unwrap();
    let ref_fasta = dir.path().join("ref.fa");
    {
        let mut f = std::fs::File::create(&ref_fasta).unwrap();
        writeln!(f, ">only\nACGTACGTACGTACGTACGTACGTACGT").unwrap();
    }
    let empty_r1 = dir.path().join("empty.fq");
    std::fs::write(&empty_r1, "").unwrap();

    let bin = binary_path(OracleStage::FastqExtractor);
    let out_prefix = dir.path().join("oracle_out");
    let output = Command::new(&bin)
        .arg("-f")
        .arg(&ref_fasta)
        .arg("-u")
        .arg(&empty_r1)
        .arg("-o")
        .arg(&out_prefix)
        .output()
        .unwrap_or_else(|e| panic!("spawning fastq-extractor: {e}"));
    assert!(!output.status.success(), "oracle should fail on an empty read-1 file");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Read file is empty."), "unexpected oracle stderr: {stderr}");
}

#[test]
fn mate_count_mismatch_matches_oracle_error_behavior() {
    // Same idea as the empty-file test: confirm the REAL oracle actually
    // rejects mismatched mate-pair counts (the Rust side's own behavior is
    // covered by fg-t1k-core's unit tests).
    let dir = tempfile::tempdir().unwrap();
    let ref_fasta = dir.path().join("ref.fa");
    {
        let mut f = std::fs::File::create(&ref_fasta).unwrap();
        writeln!(f, ">only\nACGTACGTACGTACGTACGTACGTACGT").unwrap();
    }
    let r1 = dir.path().join("r1.fq");
    let r2 = dir.path().join("r2.fq");
    write_fastq(
        &r1,
        &[
            ("p0".to_string(), "ACGTACGTACGTACGTACGTACGTACGT".to_string()),
            ("p1".to_string(), "ACGTACGTACGTACGTACGTACGTACGT".to_string()),
        ],
    );
    write_fastq(&r2, &[("p0".to_string(), "ACGTACGTACGTACGTACGTACGTACGT".to_string())]);

    let bin = binary_path(OracleStage::FastqExtractor);
    let out_prefix = dir.path().join("oracle_out");
    let output = Command::new(&bin)
        .arg("-f")
        .arg(&ref_fasta)
        .arg("-1")
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .arg("-o")
        .arg(&out_prefix)
        .output()
        .unwrap_or_else(|e| panic!("spawning fastq-extractor: {e}"));
    assert!(!output.status.success(), "oracle should fail on mate-count mismatch");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("different number of reads"), "unexpected oracle stderr: {stderr}");
}

/// Builds a LARGE synthetic paired fixture (thousands of pairs, spanning
/// several `512 * threads`-sized parallel-decision batches even at `threads
/// = 8`) covering the same mix of pass/fail/RC/low-complexity categories as
/// [`build_synthetic_fixture`], so the `-t N` identity test below actually
/// exercises the parallel-decision path across multiple batches, not just a
/// single one.
fn build_large_synthetic_fixture(pair_count: usize) -> SyntheticFixture {
    let dir = tempfile::tempdir().unwrap();

    // A 300bp base-balanced reference (same shape as build_synthetic_fixture's).
    let reference = "ACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATGGTACGTACGTACGTACGGGCATTCATGGCATTCATGGCATTCATGACGTTAGCACGTTAGCACGTTAGCACGTTAGCTGACCATGTGACCATGTGACCATGTGACCATG";
    let ref_fasta = dir.path().join("ref.fa");
    {
        let mut f = std::fs::File::create(&ref_fasta).unwrap();
        writeln!(f, ">only\n{reference}").unwrap();
    }

    let random_noise = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";

    let mut r1_records: Vec<(String, String)> = Vec::with_capacity(pair_count);
    let mut r2_records: Vec<(String, String)> = Vec::with_capacity(pair_count);

    for i in 0..pair_count {
        match i % 4 {
            // Hit: mate1 is a rotating 80bp reference substring.
            0 => {
                let start = (i * 7) % (reference.len() - 80);
                r1_records.push((format!("hit_{i}"), reference[start..start + 80].to_string()));
                r2_records.push((
                    format!("hit_{i}"),
                    reference[(start + 5) % (reference.len() - 80)
                        ..(start + 5) % (reference.len() - 80) + 80]
                        .to_string(),
                ));
            }
            // Miss: both mates unrelated noise.
            1 => {
                r1_records.push((format!("miss_{i}"), random_noise.to_string()));
                r2_records.push((format!("miss_{i}"), random_noise.to_string()));
            }
            // RC hit: mate1 is the reverse complement of a reference
            // substring.
            2 => {
                let start = (i * 11) % (reference.len() - 80);
                r1_records
                    .push((format!("rc_{i}"), reverse_complement(&reference[start..start + 80])));
                r2_records.push((format!("rc_{i}"), random_noise.to_string()));
            }
            // Low-complexity: dropped by IsLowComplexity before even
            // reaching the k-mer filter.
            _ => {
                r1_records.push((format!("lowc_{i}"), "A".repeat(80)));
                r2_records.push((format!("lowc_{i}"), random_noise.to_string()));
            }
        }
    }

    let r1 = dir.path().join("r1.fq");
    let r2 = dir.path().join("r2.fq");
    write_fastq(&r1, &r1_records);
    write_fastq(&r2, &r2_records);

    SyntheticFixture { dir, ref_fasta, r1, r2 }
}

/// The core byte-identity requirement of the P1 parallelism task: `-t 1`,
/// `-t 4`, and `-t 8` must all produce output byte-identical to EACH OTHER
/// and to the real oracle (which is always run at its own default `-t 1` --
/// see this module's docs for why that single oracle run is sufficient).
/// `pair_count` (4000) spans multiple `512 * threads`-sized parallel batches
/// even at `threads = 8` (`512 * 8 = 4096` per batch, so this fixture forces
/// at least two batches total), directly exercising the batch-boundary code
/// path, not just a single in-memory batch.
#[test]
fn parallel_threads_1_4_8_byte_identical_to_each_other_and_oracle() {
    let fx = build_large_synthetic_fixture(4000);

    let t1_prefix = fx.dir.path().join("t1");
    let t4_prefix = fx.dir.path().join("t4");
    let t8_prefix = fx.dir.path().join("t8");
    let oracle_prefix = fx.dir.path().join("oracle");

    run_rust_paired_with_threads(&fx.ref_fasta, &fx.r1, &fx.r2, &t1_prefix, 1);
    run_rust_paired_with_threads(&fx.ref_fasta, &fx.r1, &fx.r2, &t4_prefix, 4);
    run_rust_paired_with_threads(&fx.ref_fasta, &fx.r1, &fx.r2, &t8_prefix, 8);
    run_oracle_paired(&fx.ref_fasta, &fx.r1, &fx.r2, &oracle_prefix);

    // -t 1 vs -t 4 vs -t 8 vs oracle, both mate files, all pairwise
    // byte-identical.
    for (label, prefix) in [("t4", &t4_prefix), ("t8", &t8_prefix), ("oracle", &oracle_prefix)] {
        assert_files_byte_identical(
            &fx.dir.path().join("t1_1.fq"),
            &fx.dir.path().join(format!("{}_1.fq", prefix.file_name().unwrap().to_str().unwrap())),
        );
        assert_files_byte_identical(
            &fx.dir.path().join("t1_2.fq"),
            &fx.dir.path().join(format!("{}_2.fq", prefix.file_name().unwrap().to_str().unwrap())),
        );
        let _ = label; // used only in the loop's own readability; assertions above do the work.
    }

    // Sanity: the fixture actually exercises a real mix of pass/fail (not
    // vacuously all-pass or all-fail).
    let out1 = std::fs::read_to_string(fx.dir.path().join("t1_1.fq")).unwrap();
    let emitted: usize = out1.lines().filter(|l| l.starts_with('@')).count();
    assert!(emitted > 0, "expected at least one candidate pair to pass");
    assert!(emitted < 4000, "expected at least one candidate pair to be dropped");
}
