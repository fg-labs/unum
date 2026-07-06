//! Byte-golden test for `unum_core::extract::extract_candidates`, converted
//! from the retired T1K-oracle FFI differential (`diff_fastq_extract.rs`) FFI differential (see
//! `tests/common/mod.rs`). Runs the Rust extractor over the same inputs the
//! differential ran the C++ `fastq-extractor` oracle over, and asserts the
//! emitted FASTQ bytes match the committed byte-goldens (which are the Rust
//! port's output, byte-identical to the oracle's).
//!
//! # The golden locks "output order == input order"
//!
//! `FastqExtractor.cpp`'s multi-threaded path parallelizes only the per-read
//! filter DECISION, never the emission order, which is always
//! sequential-in-input-order. So a correct single-threaded Rust
//! re-implementation is byte-identical to the oracle at any `-t`; these
//! byte-goldens pin exactly that emitted stream.
//!
//! Dropped from the original: `empty_read1_file_matches_oracle_error_behavior`
//! and `mate_count_mismatch_matches_oracle_error_behavior` -- both only
//! asserted the C++ oracle's own error exit; the Rust side's Result-based
//! error behavior is covered by `unum_core::extract`'s own unit tests.

mod common;

use common::assert_byte_golden;
use unum_core::extract::{self, CandidateSink, ReadRecord};
use unum_core::ref_kmer_filter::RefKmerFilter;
use std::io::Write as _;
use std::path::{Path, PathBuf};

const INITIAL_KMER_LENGTH: usize = 9;
const DEFAULT_SIMILARITY: f64 = extract::DEFAULT_REF_SEQ_SIMILARITY;

fn fixture(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures").join(rel)
}

/// A FASTQ-writing sink mirroring the CLI's own `FastqFileSink`, so the test
/// exercises only the library surface.
struct FastqFileSink {
    fp1: std::fs::File,
    fp2: Option<std::fs::File>,
}

impl FastqFileSink {
    fn create(prefix: &Path, paired: bool) -> Self {
        let fp1 = std::fs::File::create(format!("{}_1.fq", prefix.display())).unwrap();
        let fp2 = paired.then(|| std::fs::File::create(format!("{}_2.fq", prefix.display())).unwrap());
        Self { fp1, fp2 }
    }
    fn create_single(prefix: &Path) -> Self {
        Self { fp1: std::fs::File::create(format!("{}.fq", prefix.display())).unwrap(), fp2: None }
    }
}

impl CandidateSink for FastqFileSink {
    fn emit_pair(&mut self, r1: &ReadRecord, r2: Option<&ReadRecord>) -> anyhow::Result<()> {
        extract::output_seq(&mut self.fp1, &r1.id, &r1.seq, r1.qual.as_deref(), 0, -1)?;
        if let (Some(r2), Some(fp2)) = (r2, self.fp2.as_mut()) {
            extract::output_seq(fp2, &r1.id, &r2.seq, r2.qual.as_deref(), 0, -1)?;
        }
        Ok(())
    }
}

fn run_rust_paired(ref_fasta: &Path, r1: &Path, r2: &Path, out_prefix: &Path) {
    let mut filter = RefKmerFilter::from_reference_fasta(ref_fasta, INITIAL_KMER_LENGTH).unwrap();
    let mut source = extract::open_source(r1, Some(r2)).unwrap();
    let mut sink = FastqFileSink::create(out_prefix, true);
    extract::extract_candidates(&mut source, &mut filter, DEFAULT_SIMILARITY, &mut sink).unwrap();
}

fn run_rust_single(ref_fasta: &Path, r1: &Path, out_prefix: &Path) {
    let mut filter = RefKmerFilter::from_reference_fasta(ref_fasta, INITIAL_KMER_LENGTH).unwrap();
    let mut source = extract::open_source(r1, None).unwrap();
    let mut sink = FastqFileSink::create_single(out_prefix);
    extract::extract_candidates(&mut source, &mut filter, DEFAULT_SIMILARITY, &mut sink).unwrap();
}

fn write_fastq(path: &Path, records: &[(String, String)]) {
    let mut f = std::fs::File::create(path).unwrap();
    for (id, seq) in records {
        let qual = "I".repeat(seq.len());
        writeln!(f, "@{id}\n{seq}\n+\n{qual}").unwrap();
    }
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

#[test]
fn paired_example_fixture_byte_golden() {
    let ref_fasta = fixture("example/ref/kir_rna_seq.fa");
    let r1 = fixture("example/example_1.fq");
    let r2 = fixture("example/example_2.fq");
    let tmp = tempfile::tempdir().unwrap();
    let prefix = tmp.path().join("rust");
    run_rust_paired(&ref_fasta, &r1, &r2, &prefix);

    let out1 = std::fs::read(tmp.path().join("rust_1.fq")).unwrap();
    let out2 = std::fs::read(tmp.path().join("rust_2.fq")).unwrap();
    // Sanity: non-trivial (some pair passed).
    assert!(String::from_utf8_lossy(&out1).lines().any(|l| l.starts_with('@')));
    assert_byte_golden("fastq_extract/example_paired_1.fq", &out1);
    assert_byte_golden("fastq_extract/example_paired_2.fq", &out2);
}

#[test]
fn single_end_example_fixture_byte_golden() {
    let ref_fasta = fixture("example/ref/kir_rna_seq.fa");
    let r1 = fixture("example/example_1.fq");
    let tmp = tempfile::tempdir().unwrap();
    let prefix = tmp.path().join("rust");
    run_rust_single(&ref_fasta, &r1, &prefix);
    let out = std::fs::read(tmp.path().join("rust.fq")).unwrap();
    assert_byte_golden("fastq_extract/example_single.fq", &out);
}

#[test]
fn mate2_longer_byte_golden() {
    let reference = "ACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATGGTACGTACGTACGTACGGGCATTCATGGCATTCATGGCATTCATGACGTTAGCACGTTAGCACGTTAGCACGTTAGCTGACCATGTGACCATGTGACCATGTGACCATG";
    let dir = tempfile::tempdir().unwrap();
    let ref_fasta = dir.path().join("ref.fa");
    {
        let mut f = std::fs::File::create(&ref_fasta).unwrap();
        writeln!(f, ">only\n{reference}").unwrap();
    }
    let hit1 = reference[0..100].to_string();
    let hit2 = reference[120..220].to_string();
    let r1_records: Vec<(String, String)> = (0..3).map(|i| (format!("p{i}"), hit1.clone())).collect();
    let mut r2_records: Vec<(String, String)> = (0..3).map(|i| (format!("p{i}"), hit2.clone())).collect();
    r2_records.push(("extra0".to_string(), hit1.clone()));
    r2_records.push(("extra1".to_string(), hit2.clone()));
    let r1 = dir.path().join("r1.fq");
    let r2 = dir.path().join("r2.fq");
    write_fastq(&r1, &r1_records);
    write_fastq(&r2, &r2_records);
    let prefix = dir.path().join("rust");
    run_rust_paired(&ref_fasta, &r1, &r2, &prefix);
    assert_byte_golden("fastq_extract/mate2_longer_1.fq", &std::fs::read(dir.path().join("rust_1.fq")).unwrap());
    assert_byte_golden("fastq_extract/mate2_longer_2.fq", &std::fs::read(dir.path().join("rust_2.fq")).unwrap());
}

fn build_synthetic_fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let reference = "ACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATGGTACGTACGTACGTACGGGCATTCATGGCATTCATGGCATTCATGACGTTAGCACGTTAGCACGTTAGCACGTTAGCTGACCATGTGACCATGTGACCATGTGACCATG";
    let ref_fasta = dir.join("ref.fa");
    {
        let mut f = std::fs::File::create(&ref_fasta).unwrap();
        writeln!(f, ">only\n{reference}").unwrap();
    }
    let mate1_hit = &reference[0..100];
    let mate2_hit = &reference[150..250];
    let noise_100 = "A".repeat(50) + &"T".repeat(50);
    let random_noise = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";
    let homopolymer = "A".repeat(100);
    let rc_hit = reverse_complement(&reference[50..150]);

    let mut r1_records: Vec<(String, String)> = Vec::new();
    let mut r2_records: Vec<(String, String)> = Vec::new();
    r1_records.push(("pair_mate1_hit".to_string(), mate1_hit.to_string()));
    r2_records.push(("pair_mate1_hit".to_string(), random_noise.to_string()));
    r1_records.push(("pair_mate2_hit".to_string(), random_noise.to_string()));
    r2_records.push(("pair_mate2_hit".to_string(), mate2_hit.to_string()));
    r1_records.push(("pair_neither_hit".to_string(), random_noise.to_string()));
    r2_records.push(("pair_neither_hit".to_string(), noise_100.clone()));
    r1_records.push(("pair_low_complexity".to_string(), homopolymer));
    r2_records.push(("pair_low_complexity".to_string(), noise_100));
    r1_records.push(("pair_rc_hit".to_string(), rc_hit));
    r2_records.push(("pair_rc_hit".to_string(), random_noise.to_string()));
    for i in 0..20 {
        let start = (i * 3) % (reference.len() - 80);
        r1_records.push((format!("pad_hit_{i}"), reference[start..start + 80].to_string()));
        r2_records.push((format!("pad_hit_{i}"), reference[start + 5..start + 85].to_string()));
    }
    let r1 = dir.join("r1.fq");
    let r2 = dir.join("r2.fq");
    write_fastq(&r1, &r1_records);
    write_fastq(&r2, &r2_records);
    (ref_fasta, r1, r2)
}

#[test]
fn synthetic_branch_coverage_byte_golden() {
    let dir = tempfile::tempdir().unwrap();
    let (ref_fasta, r1, r2) = build_synthetic_fixture(dir.path());
    let prefix = dir.path().join("rust");
    run_rust_paired(&ref_fasta, &r1, &r2, &prefix);

    let out1 = std::fs::read(dir.path().join("rust_1.fq")).unwrap();
    let out2 = std::fs::read(dir.path().join("rust_2.fq")).unwrap();

    // Sanity: the expected pass/fail mix (proves the test exercises the filter).
    let emitted: Vec<&str> =
        std::str::from_utf8(&out1).unwrap().lines().filter_map(|l| l.strip_prefix('@')).collect();
    assert!(emitted.contains(&"pair_mate1_hit"));
    assert!(emitted.contains(&"pair_mate2_hit"));
    assert!(emitted.contains(&"pair_rc_hit"));
    assert!(!emitted.contains(&"pair_neither_hit"));
    assert!(!emitted.contains(&"pair_low_complexity"));

    assert_byte_golden("fastq_extract/synthetic_1.fq", &out1);
    assert_byte_golden("fastq_extract/synthetic_2.fq", &out2);
}
