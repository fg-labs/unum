#![cfg(feature = "t1k-sys")]
//! Differential test: runs the REAL C++ `bam-extractor` oracle binary and
//! the Rust `fg_t1k_core::bam_extract::extract_from_bam` port on the SAME
//! `_coord.fa` + BAM inputs, and asserts the output FASTQ file(s) are
//! **byte-identical**, not merely set-equal or order-insensitive.
//!
//! # Subprocess, not in-process FFI -- symbol collision
//!
//! `vendor/t1k/samtools-0.1.19` (which `bam-extractor` links against
//! directly) and `rust-htslib`'s vendored modern htslib (which
//! `fg_t1k_core::alignments::Alignments` -- and therefore
//! `fg_t1k_core::bam_extract` -- links against) export a large overlapping
//! set of IDENTICALLY NAMED C symbols with INCOMPATIBLE ABIs (`bam_read1`,
//! `sam_read1`, `bam_write1`, every `bgzf_*`/`fai_*` function, etc.). Linking
//! both into one process is a genuine, unavoidable ABI collision (confirmed
//! in `crates/fg-t1k-sys/tests/diff_alignments.rs`'s module docs, which hit
//! the exact same issue for the Task-3.3a `Alignments` differential). This
//! test therefore runs the REAL `bam-extractor` oracle binary as a spawned
//! SUBPROCESS (`fg_t1k_sys::oracle::binary_path(OracleStage::BamExtractor)`
//! plus `std::process::Command`), never via in-process FFI -- this test
//! binary links only rust-htslib (via `fg_t1k_core`), and the spawned oracle
//! process links only samtools-0.1.19, so neither process has both.
//!
//! # `-t 1` on both sides makes output directly comparable
//!
//! `fg_t1k_core::bam_extract`'s module docs explain why this port only
//! reproduces `BamExtractor.cpp`'s `threadCnt == 1` code path: unlike
//! `FastqExtractor.cpp`, `BamExtractor.cpp`'s multi-threaded path changes
//! BATCHING/FLUSH TIMING for the unaligned-template-pair and single-end-
//! unmapped candidate paths, so output order/batching is NOT provably
//! thread-count-invariant there. This test always runs the oracle with
//! `-t 1` (matching the Rust port's only-ever-single-threaded semantics), so
//! both sides are directly, deterministically comparable byte-for-byte.
//!
//! # Fixture coverage
//!
//! [`build_test_bam`] builds a coordinate-sorted BAM covering every branch
//! documented in `fg_t1k_core::bam_extract`'s module docs:
//! - **On-target aligned pairs**: KIR2DL1-gene-coordinate-mapped read pairs
//!   whose SEQ is an exact substring of the KIR2DL1 reference sequence (from
//!   `fixtures/refbuild/golden/kir_rna_coord.fa`), landing inside that
//!   gene's `[start, end]` interval -- exercises the `tag`-advance/overlap
//!   gene-interval logic and paired-candidates-then-pass-2-completion path.
//! - **Unaligned-template pair**: two consecutive unmapped records (KIR-
//!   substring SEQ so they pass the candidate filter) -- exercises the
//!   direct pass-1 emission path for unaligned templates.
//! - **Alt-chrom aligned read pair**: mapped to a contig whose name contains
//!   `_` (`chr19_KI270938v1_alt`), KIR-substring SEQ -- exercises
//!   `ValidAlternativeChrom` + the paired-candidates path (same code path as
//!   an off-target aligned read, just reached via the chrom-name check
//!   instead of the gene-interval check).
//! - **Off-target (non-gene-overlapping) aligned pair on the primary
//!   chrom**: mapped to `chr19` but far outside every KIR gene interval --
//!   must be silently dropped (falls through the `tag`-advance `continue`
//!   before ever reaching the low-complexity/candidate check).
//! - **Low-complexity on-target pair**: mapped inside a KIR gene interval
//!   but a homopolymer SEQ -- must be dropped by `IsLowComplexity` despite
//!   overlapping a gene.
//!
//! Every category is paired (this port's fixture uses a paired BAM;
//! [`build_single_end_test_bam`] provides a separate single-end-only BAM for
//! the `frag_stdev == 0` code path, checked in a separate test).

use fg_t1k_core::alignments::Alignments;
use fg_t1k_core::bam_extract::{self, CoordRecord};
use fg_t1k_core::extract::{CandidateSink, ReadRecord, output_seq};
use fg_t1k_core::ref_kmer_filter::RefKmerFilter;
use fg_t1k_sys::oracle::{OracleStage, binary_path};
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{self, Header, Writer};
use std::path::{Path, PathBuf};
use std::process::Command;

/// `BamExtractor.cpp:480`: the literal initial k-mer length.
const INITIAL_KMER_LENGTH: usize = 9;

fn fixture(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures").join(rel)
}

fn coord_fa() -> PathBuf {
    fixture("refbuild/golden/kir_rna_coord.fa")
}

/// KIR2DL1's genomic interval per `kir_rna_coord.fa`'s header
/// (`>KIR2DL1*0010101 chr19 54769793 54784332 +`).
const KIR2DL1_CHROM: &str = "chr19";
const KIR2DL1_START: i64 = 54_769_793;

/// Loads `KIR2DL1*0010101`'s sequence out of `kir_rna_coord.fa` directly
/// (rather than hardcoding a copy here), so this fixture always tracks
/// whatever the golden coord FASTA actually contains.
fn kir2dl1_sequence() -> String {
    let records =
        bam_extract::parse_coord_fa(&coord_fa()).expect("parse kir_rna_coord.fa for fixture");
    records
        .into_iter()
        .find(|r: &CoordRecord| r.name == "KIR2DL1*0010101")
        .expect("KIR2DL1*0010101 record present in kir_rna_coord.fa")
        .seq
}

/// A base-balanced, non-low-complexity 100bp filler sequence used for
/// off-target/noise reads that must NOT be candidates.
fn noise_seq(len: usize) -> Vec<u8> {
    let pattern =
        b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
    (0..len).map(|i| pattern[i % pattern.len()]).collect()
}

/// Builds the paired-input coordinate-sorted BAM covering every branch
/// documented in this module's docs. Two contigs: `chr19` (length spanning
/// every KIR gene coordinate) and `chr19_KI270938v1_alt` (an alt contig,
/// name deliberately containing `_` to exercise `ValidAlternativeChrom`).
// Long by construction: each record is a distinct, deliberately documented
// coverage category (see module docs) -- splitting this into several tiny
// functions would scatter closely related per-record flag/CIGAR bookkeeping
// across the file for no readability benefit (same rationale as
// diff_alignments.rs's build_test_bam).
#[allow(clippy::too_many_lines)]
fn build_test_bam(path: &Path) {
    let kir_seq = kir2dl1_sequence();
    let kir_bytes = kir_seq.as_bytes();

    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);

    let mut sq_chr19 = HeaderRecord::new(b"SQ");
    sq_chr19.push_tag(b"SN", KIR2DL1_CHROM);
    sq_chr19.push_tag(b"LN", 58_617_616);
    header.push_record(&sq_chr19);

    let mut sq_alt = HeaderRecord::new(b"SQ");
    sq_alt.push_tag(b"SN", "chr19_KI270938v1_alt");
    sq_alt.push_tag(b"LN", 1_000_000);
    header.push_record(&sq_alt);

    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

    // --- Category 1: on-target aligned pair, mapped inside KIR2DL1's gene
    // interval on chr19, SEQ = exact 100bp substrings of KIR2DL1's own
    // reference sequence (guarantees HasHitInSet passes).
    let on_target_seq1 = &kir_bytes[0..100];
    let on_target_seq2 = &kir_bytes[200..300];
    let pos1 = KIR2DL1_START + 10; // inside [start, end]
    let pos2 = KIR2DL1_START + 210;

    let mut ot1 = bam::Record::new();
    ot1.set(
        b"on_target_pair",
        Some(&CigarString(vec![Cigar::Match(100)])),
        on_target_seq1,
        &[30u8; 100],
    );
    ot1.set_tid(0);
    ot1.set_pos(pos1);
    ot1.set_mtid(0);
    ot1.set_mpos(pos2);
    ot1.set_flags(0x1 | 0x2 | 0x20 | 0x40); // paired, proper, mate_reverse, first
    writer.write(&ot1).expect("write on-target mate 1");

    let mut ot2 = bam::Record::new();
    ot2.set(
        b"on_target_pair",
        Some(&CigarString(vec![Cigar::Match(100)])),
        on_target_seq2,
        &[30u8; 100],
    );
    ot2.set_tid(0);
    ot2.set_pos(pos2);
    ot2.set_mtid(0);
    ot2.set_mpos(pos1);
    ot2.set_flags(0x1 | 0x2 | 0x10 | 0x80); // paired, proper, reverse, last
    writer.write(&ot2).expect("write on-target mate 2");

    // --- Category 2: off-target aligned pair on chr19, FAR from every KIR
    // gene interval (position 1000, well before KIR2DP1's earliest start of
    // 54,755,023) -- must be silently dropped by the tag-advance/overlap
    // check before any candidate/low-complexity check runs.
    let off_target_seq = noise_seq(80);
    let mut off1 = bam::Record::new();
    off1.set(
        b"off_target_pair",
        Some(&CigarString(vec![Cigar::Match(80)])),
        &off_target_seq,
        &[30u8; 80],
    );
    off1.set_tid(0);
    off1.set_pos(1000);
    off1.set_mtid(0);
    off1.set_mpos(1200);
    off1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&off1).expect("write off-target mate 1");

    let mut off2 = bam::Record::new();
    off2.set(
        b"off_target_pair",
        Some(&CigarString(vec![Cigar::Match(80)])),
        &off_target_seq,
        &[30u8; 80],
    );
    off2.set_tid(0);
    off2.set_pos(1200);
    off2.set_mtid(0);
    off2.set_mpos(1000);
    off2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&off2).expect("write off-target mate 2");

    // --- Category 3: low-complexity pair, mapped INSIDE KIR2DL1's gene
    // interval (so it clears the gene-overlap check) but a homopolymer SEQ
    // that IsLowComplexity must reject.
    let homopolymer = vec![b'A'; 90];
    let lc_pos1 = KIR2DL1_START + 1000;
    let lc_pos2 = KIR2DL1_START + 1200;
    let mut lc1 = bam::Record::new();
    lc1.set(
        b"low_complexity_pair",
        Some(&CigarString(vec![Cigar::Match(90)])),
        &homopolymer,
        &[30u8; 90],
    );
    lc1.set_tid(0);
    lc1.set_pos(lc_pos1);
    lc1.set_mtid(0);
    lc1.set_mpos(lc_pos2);
    lc1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&lc1).expect("write low-complexity mate 1");

    let mut lc2 = bam::Record::new();
    lc2.set(
        b"low_complexity_pair",
        Some(&CigarString(vec![Cigar::Match(90)])),
        &homopolymer,
        &[30u8; 90],
    );
    lc2.set_tid(0);
    lc2.set_pos(lc_pos2);
    lc2.set_mtid(0);
    lc2.set_mpos(lc_pos1);
    lc2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&lc2).expect("write low-complexity mate 2");

    // --- Category 4: alt-chrom aligned pair, mapped to
    // "chr19_KI270938v1_alt" (name contains '_') with KIR-substring SEQ --
    // exercises ValidAlternativeChrom's paired-candidates path.
    let alt_seq1 = &kir_bytes[400..480];
    let alt_seq2 = &kir_bytes[600..680];
    let mut alt1 = bam::Record::new();
    alt1.set(b"alt_chrom_pair", Some(&CigarString(vec![Cigar::Match(80)])), alt_seq1, &[30u8; 80]);
    alt1.set_tid(1);
    alt1.set_pos(500);
    alt1.set_mtid(1);
    alt1.set_mpos(700);
    alt1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&alt1).expect("write alt-chrom mate 1");

    let mut alt2 = bam::Record::new();
    alt2.set(b"alt_chrom_pair", Some(&CigarString(vec![Cigar::Match(80)])), alt_seq2, &[30u8; 80]);
    alt2.set_tid(1);
    alt2.set_pos(700);
    alt2.set_mtid(1);
    alt2.set_mpos(500);
    alt2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&alt2).expect("write alt-chrom mate 2");

    // --- Category 5: unaligned-template pair -- two CONSECUTIVE unmapped
    // records, KIR-substring SEQ so they pass the candidate filter.
    // Coordinate-sorted BAMs conventionally place unmapped-template records
    // at the very end (tid=-1); this port's driver only requires the two
    // mates be ADJACENT, not any particular position relative to mapped
    // records, so placing them last here (matching real-world coordinate-
    // sort behavior) is both realistic and sufficient.
    let unmapped_seq1 = &kir_bytes[800..880];
    let unmapped_seq2 = &kir_bytes[900..980];
    let mut um1 = bam::Record::new();
    um1.set(b"unaligned_template_pair", None, unmapped_seq1, &[25u8; 80]);
    um1.set_tid(-1);
    um1.set_pos(-1);
    um1.set_mtid(-1);
    um1.set_mpos(-1);
    um1.set_flags(0x1 | 0x4 | 0x8 | 0x40); // paired, unmapped, mate_unmapped, first
    writer.write(&um1).expect("write unaligned-template mate 1");

    let mut um2 = bam::Record::new();
    um2.set(b"unaligned_template_pair", None, unmapped_seq2, &[25u8; 80]);
    um2.set_tid(-1);
    um2.set_pos(-1);
    um2.set_mtid(-1);
    um2.set_mpos(-1);
    um2.set_flags(0x1 | 0x4 | 0x8 | 0x80); // paired, unmapped, mate_unmapped, last
    writer.write(&um2).expect("write unaligned-template mate 2");

    // --- Category 6: an unaligned-template pair that does NOT pass the
    // candidate filter (both mates noise) -- must be silently dropped.
    let noise1 = noise_seq(70);
    let noise2 = noise_seq(70);
    let mut noise_um1 = bam::Record::new();
    noise_um1.set(b"unaligned_noise_pair", None, &noise1, &[25u8; 70]);
    noise_um1.set_tid(-1);
    noise_um1.set_pos(-1);
    noise_um1.set_mtid(-1);
    noise_um1.set_mpos(-1);
    noise_um1.set_flags(0x1 | 0x4 | 0x8 | 0x40);
    writer.write(&noise_um1).expect("write noise unaligned-template mate 1");

    let mut noise_um2 = bam::Record::new();
    noise_um2.set(b"unaligned_noise_pair", None, &noise2, &[25u8; 70]);
    noise_um2.set_tid(-1);
    noise_um2.set_pos(-1);
    noise_um2.set_mtid(-1);
    noise_um2.set_mpos(-1);
    noise_um2.set_flags(0x1 | 0x4 | 0x8 | 0x80);
    writer.write(&noise_um2).expect("write noise unaligned-template mate 2");

    drop(writer);
}

/// Builds a single-end-only BAM (no paired flag on any record) covering the
/// `frag_stdev == 0` code path: on-target aligned reads (KIR-substring SEQ,
/// mapped inside KIR2DL1's gene interval, with a DUPLICATE read id to
/// exercise `usedName` dedup), an off-target aligned read (dropped by the
/// gene-interval check), and a low-complexity on-target read (dropped by
/// `IsLowComplexity`).
fn build_single_end_test_bam(path: &Path) {
    let kir_seq = kir2dl1_sequence();
    let kir_bytes = kir_seq.as_bytes();

    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);
    let mut sq_chr19 = HeaderRecord::new(b"SQ");
    sq_chr19.push_tag(b"SN", KIR2DL1_CHROM);
    sq_chr19.push_tag(b"LN", 58_617_616);
    header.push_record(&sq_chr19);

    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

    let on_target_seq = &kir_bytes[0..100];
    let mut r1 = bam::Record::new();
    r1.set(
        b"single_on_target",
        Some(&CigarString(vec![Cigar::Match(100)])),
        on_target_seq,
        &[30u8; 100],
    );
    r1.set_tid(0);
    r1.set_pos(KIR2DL1_START + 10);
    r1.set_mtid(-1);
    r1.set_mpos(-1);
    r1.set_flags(0); // unpaired, forward
    writer.write(&r1).expect("write single on-target read");

    // Same read id, second alignment (e.g. multi-mapped to an alt-like
    // region) -- exercises usedName dedup on the single-end path.
    let mut r1_dup = bam::Record::new();
    r1_dup.set(
        b"single_on_target",
        Some(&CigarString(vec![Cigar::Match(100)])),
        on_target_seq,
        &[30u8; 100],
    );
    r1_dup.set_tid(0);
    r1_dup.set_pos(KIR2DL1_START + 20);
    r1_dup.set_mtid(-1);
    r1_dup.set_mpos(-1);
    r1_dup.set_flags(0x100); // secondary
    writer.write(&r1_dup).expect("write duplicate-name single read");

    let off_target_seq = noise_seq(80);
    let mut r2 = bam::Record::new();
    r2.set(
        b"single_off_target",
        Some(&CigarString(vec![Cigar::Match(80)])),
        &off_target_seq,
        &[30u8; 80],
    );
    r2.set_tid(0);
    r2.set_pos(2000);
    r2.set_mtid(-1);
    r2.set_mpos(-1);
    r2.set_flags(0);
    writer.write(&r2).expect("write single off-target read");

    let homopolymer = vec![b'T'; 90];
    let mut r3 = bam::Record::new();
    r3.set(
        b"single_low_complexity",
        Some(&CigarString(vec![Cigar::Match(90)])),
        &homopolymer,
        &[30u8; 90],
    );
    r3.set_tid(0);
    r3.set_pos(KIR2DL1_START + 500);
    r3.set_mtid(-1);
    r3.set_mpos(-1);
    r3.set_flags(0);
    writer.write(&r3).expect("write single low-complexity read");

    drop(writer);
}

/// A minimal [`CandidateSink`] writing to `{prefix}_1.fq`/`{prefix}_2.fq`
/// (paired) or `{prefix}.fq` (single-end), mirroring the CLI's own
/// `FastqFileSink` -- kept as a small self-contained copy so this
/// differential test exercises the LIBRARY surface directly (matching the
/// FASTQ-extractor differential's own pattern).
struct FastqFileSink {
    fp1: std::fs::File,
    fp2: Option<std::fs::File>,
}

impl FastqFileSink {
    fn create_paired(prefix: &Path) -> Self {
        let fp1 = std::fs::File::create(format!("{}_1.fq", prefix.display()))
            .unwrap_or_else(|e| panic!("create {prefix:?}_1.fq: {e}"));
        let fp2 = std::fs::File::create(format!("{}_2.fq", prefix.display()))
            .unwrap_or_else(|e| panic!("create {prefix:?}_2.fq: {e}"));
        Self { fp1, fp2: Some(fp2) }
    }

    fn create_single(prefix: &Path) -> Self {
        let fp1 = std::fs::File::create(format!("{}.fq", prefix.display()))
            .unwrap_or_else(|e| panic!("create {prefix:?}.fq: {e}"));
        Self { fp1, fp2: None }
    }
}

impl CandidateSink for FastqFileSink {
    fn emit_pair(&mut self, r1: &ReadRecord, r2: Option<&ReadRecord>) -> anyhow::Result<()> {
        output_seq(&mut self.fp1, &r1.id, &r1.seq, r1.qual.as_deref(), 0, -1)?;
        if let (Some(r2), Some(fp2)) = (r2, self.fp2.as_mut()) {
            output_seq(fp2, &r1.id, &r2.seq, r2.qual.as_deref(), 0, -1)?;
        }
        Ok(())
    }
}

/// Runs the Rust `bam_extract` driver end-to-end (library call, not the
/// CLI), writing `{prefix}_1.fq`/`{prefix}_2.fq` (paired) or `{prefix}.fq`
/// (single-end, detected from the BAM's own sampled `frag_stdev`) to
/// `out_prefix`'s directory.
fn run_rust(
    bam_path: &Path,
    coord_fasta: &Path,
    out_prefix: &Path,
) -> bam_extract::BamExtractMetrics {
    run_rust_with_threads(bam_path, coord_fasta, out_prefix, 1)
}

/// Same as [`run_rust`], but drives the parallel pass-1 decision path via
/// [`bam_extract::extract_from_bam_with_threads`] at the given `threads`.
fn run_rust_with_threads(
    bam_path: &Path,
    coord_fasta: &Path,
    out_prefix: &Path,
    threads: usize,
) -> bam_extract::BamExtractMetrics {
    let coord_records =
        bam_extract::parse_coord_fa(coord_fasta).unwrap_or_else(|e| panic!("parse_coord_fa: {e}"));
    let mut filter = RefKmerFilter::from_reference_fasta(coord_fasta, INITIAL_KMER_LENGTH)
        .unwrap_or_else(|e| panic!("RefKmerFilter::from_reference_fasta: {e}"));
    let mut alignments =
        Alignments::open(bam_path).unwrap_or_else(|e| panic!("Alignments::open: {e}"));
    let genes = bam_extract::build_genes(&alignments, &coord_records)
        .unwrap_or_else(|e| panic!("build_genes: {e}"));

    let single_end =
        alignments.general_info(true).unwrap_or_else(|e| panic!("general_info: {e}")).frag_stdev
            == 0;
    alignments.rewind().unwrap_or_else(|e| panic!("rewind: {e}"));

    let mut sink = if single_end {
        FastqFileSink::create_single(out_prefix)
    } else {
        FastqFileSink::create_paired(out_prefix)
    };

    bam_extract::extract_from_bam_with_threads(
        &mut alignments,
        &mut filter,
        &genes,
        false,
        -1,
        threads,
        &mut sink,
    )
    .unwrap_or_else(|e| panic!("extract_from_bam_with_threads(threads={threads}): {e}"))
}

/// Runs the real oracle `bam-extractor` binary, always at `-t 1` (see module
/// docs for why this makes output directly comparable to the Rust port).
fn run_oracle(bam_path: &Path, coord_fasta: &Path, out_prefix: &Path) {
    let bin = binary_path(OracleStage::BamExtractor);
    assert!(bin.exists(), "oracle binary not built: {bin:?}");
    let output = Command::new(&bin)
        .arg("-b")
        .arg(bam_path)
        .arg("-f")
        .arg(coord_fasta)
        .arg("-o")
        .arg(out_prefix)
        .arg("-t")
        .arg("1")
        .output()
        .unwrap_or_else(|e| panic!("spawning bam-extractor: {e}"));
    assert!(
        output.status.success(),
        "bam-extractor exited with {:?}; stderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
}

fn assert_files_byte_identical(rust_path: &Path, oracle_path: &Path) {
    let rust_bytes = std::fs::read(rust_path)
        .unwrap_or_else(|e| panic!("reading rust output {rust_path:?}: {e}"));
    let oracle_bytes = std::fs::read(oracle_path)
        .unwrap_or_else(|e| panic!("reading oracle output {oracle_path:?}: {e}"));
    assert_eq!(
        rust_bytes,
        oracle_bytes,
        "byte mismatch between {rust_path:?} ({} bytes) and {oracle_path:?} ({} bytes)",
        rust_bytes.len(),
        oracle_bytes.len()
    );
}

#[test]
fn paired_synthetic_fixture_byte_identical_to_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = dir.path().join("test.bam");
    build_test_bam(&bam_path);

    let coord_fasta = coord_fa();
    let rust_prefix = dir.path().join("rust");
    let oracle_prefix = dir.path().join("oracle");

    let metrics = run_rust(&bam_path, &coord_fasta, &rust_prefix);
    assert!(!metrics.single_end, "fixture BAM must be treated as paired");
    run_oracle(&bam_path, &coord_fasta, &oracle_prefix);

    assert_files_byte_identical(&dir.path().join("rust_1.fq"), &dir.path().join("oracle_1.fq"));
    assert_files_byte_identical(&dir.path().join("rust_2.fq"), &dir.path().join("oracle_2.fq"));

    // Sanity: confirm the expected categories actually made it through (not
    // vacuously all-pass or all-fail), proving this test exercises the
    // branches it claims to.
    let out1 = std::fs::read_to_string(dir.path().join("rust_1.fq")).unwrap();
    let emitted_ids: Vec<&str> = out1.lines().filter_map(|l| l.strip_prefix('@')).collect();
    assert!(
        emitted_ids.contains(&"on_target_pair"),
        "expected on-target aligned pair to be emitted, got: {emitted_ids:?}"
    );
    assert!(
        emitted_ids.contains(&"alt_chrom_pair"),
        "expected alt-chrom pair to be emitted, got: {emitted_ids:?}"
    );
    assert!(
        emitted_ids.contains(&"unaligned_template_pair"),
        "expected unaligned-template pair to be emitted, got: {emitted_ids:?}"
    );
    assert!(
        !emitted_ids.contains(&"off_target_pair"),
        "off-target pair (outside every gene interval) must be dropped, got: {emitted_ids:?}"
    );
    assert!(
        !emitted_ids.contains(&"low_complexity_pair"),
        "low-complexity on-target pair must be dropped, got: {emitted_ids:?}"
    );
    assert!(
        !emitted_ids.contains(&"unaligned_noise_pair"),
        "unaligned-template noise pair must be dropped, got: {emitted_ids:?}"
    );
}

#[test]
fn single_end_synthetic_fixture_byte_identical_to_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = dir.path().join("single.bam");
    build_single_end_test_bam(&bam_path);

    let coord_fasta = coord_fa();
    let rust_prefix = dir.path().join("rust");
    let oracle_prefix = dir.path().join("oracle");

    let metrics = run_rust(&bam_path, &coord_fasta, &rust_prefix);
    assert!(metrics.single_end, "fixture BAM must be treated as single-end");
    run_oracle(&bam_path, &coord_fasta, &oracle_prefix);

    assert_files_byte_identical(&dir.path().join("rust.fq"), &dir.path().join("oracle.fq"));

    let out = std::fs::read_to_string(dir.path().join("rust.fq")).unwrap();
    let emitted_ids: Vec<&str> = out.lines().filter_map(|l| l.strip_prefix('@')).collect();
    // usedName dedup: the on-target read appears twice in the BAM (primary +
    // secondary alignment, same QNAME) but must be emitted at most once.
    assert_eq!(
        emitted_ids.iter().filter(|&&id| id == "single_on_target").count(),
        1,
        "usedName dedup must emit the duplicate-QNAME aligned read only once, got: {emitted_ids:?}"
    );
    assert!(
        !emitted_ids.contains(&"single_off_target"),
        "off-target single-end read must be dropped, got: {emitted_ids:?}"
    );
    assert!(
        !emitted_ids.contains(&"single_low_complexity"),
        "low-complexity single-end read must be dropped, got: {emitted_ids:?}"
    );
}

/// Regression for the `-u`/`abnormal_unaligned_flag` path: WITHOUT `-u`, an
/// unmapped-template pair whose two mates are NOT adjacent records is a
/// hard error on both sides (`BamExtractor.cpp:657-672`). This test confirms
/// the REAL oracle rejects it (the Rust side's identical behavior is proven
/// by construction -- `handle_unaligned_template_pair` bails the same way --
/// but this confirms the oracle's error text/exit-code assumption this port
/// relies on is real, not imagined).
#[test]
fn missing_unaligned_mate_is_an_error_on_the_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = dir.path().join("missing_mate.bam");

    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", KIR2DL1_CHROM);
    sq.push_tag(b"LN", 58_617_616);
    header.push_record(&sq);

    {
        let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
        let seq = kir2dl1_sequence();
        let mut r = bam::Record::new();
        r.set(b"lonely_unmapped", None, &seq.as_bytes()[0..60], &[25u8; 60]);
        r.set_tid(-1);
        r.set_pos(-1);
        r.set_mtid(-1);
        r.set_mpos(-1);
        r.set_flags(0x1 | 0x4 | 0x8 | 0x40); // paired, unmapped, mate_unmapped, first -- NO mate follows
        writer.write(&r).unwrap();
    }

    let bin = binary_path(OracleStage::BamExtractor);
    let out_prefix = dir.path().join("oracle_out");
    let output = Command::new(&bin)
        .arg("-b")
        .arg(&bam_path)
        .arg("-f")
        .arg(coord_fa())
        .arg("-o")
        .arg(&out_prefix)
        .arg("-t")
        .arg("1")
        .output()
        .unwrap_or_else(|e| panic!("spawning bam-extractor: {e}"));
    assert!(!output.status.success(), "oracle should fail when the unaligned mate is missing");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Two reads from the unaligned fragment are not showing up together"),
        "unexpected oracle stderr: {stderr}"
    );

    // The Rust side must reject the exact same input the exact same way.
    let coord_records = bam_extract::parse_coord_fa(&coord_fa()).unwrap();
    let mut filter = RefKmerFilter::from_reference_fasta(&coord_fa(), INITIAL_KMER_LENGTH).unwrap();
    let mut alignments = Alignments::open(&bam_path).unwrap();
    let genes = bam_extract::build_genes(&alignments, &coord_records).unwrap();
    let mut sink = FastqFileSink::create_paired(&dir.path().join("rust_out"));
    let result =
        bam_extract::extract_from_bam(&mut alignments, &mut filter, &genes, false, -1, &mut sink);
    assert!(result.is_err(), "rust side must also reject a missing unaligned mate");
    assert!(
        result.unwrap_err().to_string().contains("not showing up together"),
        "rust error message should match the oracle's semantics"
    );
}

/// Builds a LARGE coordinate-sorted paired BAM (thousands of on-target
/// pairs, spanning several `512 * threads`-sized parallel pass-1-decision
/// batches even at `threads = 8`) plus one of every other pass-1 category
/// ([`build_test_bam`]'s categories), so the `-t N` identity test below
/// actually exercises the parallel-decision batch-boundary path, not just a
/// single in-memory batch.
#[allow(clippy::too_many_lines)]
fn build_large_test_bam(path: &Path, pair_count: usize) {
    let kir_seq = kir2dl1_sequence();
    let kir_bytes = kir_seq.as_bytes();

    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);
    let mut sq_chr19 = HeaderRecord::new(b"SQ");
    sq_chr19.push_tag(b"SN", KIR2DL1_CHROM);
    sq_chr19.push_tag(b"LN", 58_617_616);
    header.push_record(&sq_chr19);
    let mut sq_alt = HeaderRecord::new(b"SQ");
    sq_alt.push_tag(b"SN", "chr19_KI270938v1_alt");
    sq_alt.push_tag(b"LN", 1_000_000);
    header.push_record(&sq_alt);

    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

    // Many on-target pairs, at monotonically increasing coordinates within
    // KIR2DL1's interval, rotating through the reference sequence so they're
    // non-identical (still all exact substrings, so every one is a
    // candidate).
    let usable_len = kir_bytes.len() - 100;
    for i in 0..pair_count {
        let off1 = (i * 13) % usable_len;
        let off2 = (i * 13 + 40) % usable_len;
        let seq1 = &kir_bytes[off1..off1 + 60];
        let seq2 = &kir_bytes[off2..off2 + 60];
        let name = format!("on_target_{i}");
        let pos1 = KIR2DL1_START + 10 + i64::try_from(i).unwrap() * 4;
        let pos2 = pos1 + 300;

        let mut r1 = bam::Record::new();
        r1.set(name.as_bytes(), Some(&CigarString(vec![Cigar::Match(60)])), seq1, &[30u8; 60]);
        r1.set_tid(0);
        r1.set_pos(pos1);
        r1.set_mtid(0);
        r1.set_mpos(pos2);
        r1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&r1).unwrap();

        let mut r2 = bam::Record::new();
        r2.set(name.as_bytes(), Some(&CigarString(vec![Cigar::Match(60)])), seq2, &[30u8; 60]);
        r2.set_tid(0);
        r2.set_pos(pos2);
        r2.set_mtid(0);
        r2.set_mpos(pos1);
        r2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&r2).unwrap();
    }

    // One off-target pair (dropped before any candidate check).
    let off_target_seq = noise_seq(80);
    let mut off1 = bam::Record::new();
    off1.set(
        b"off_target_pair",
        Some(&CigarString(vec![Cigar::Match(80)])),
        &off_target_seq,
        &[30u8; 80],
    );
    off1.set_tid(0);
    off1.set_pos(1000);
    off1.set_mtid(0);
    off1.set_mpos(1200);
    off1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&off1).unwrap();
    let mut off2 = bam::Record::new();
    off2.set(
        b"off_target_pair",
        Some(&CigarString(vec![Cigar::Match(80)])),
        &off_target_seq,
        &[30u8; 80],
    );
    off2.set_tid(0);
    off2.set_pos(1200);
    off2.set_mtid(0);
    off2.set_mpos(1000);
    off2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&off2).unwrap();

    // One alt-chrom pair.
    let alt_seq1 = &kir_bytes[400..480];
    let alt_seq2 = &kir_bytes[600..680];
    let mut alt1 = bam::Record::new();
    alt1.set(b"alt_chrom_pair", Some(&CigarString(vec![Cigar::Match(80)])), alt_seq1, &[30u8; 80]);
    alt1.set_tid(1);
    alt1.set_pos(500);
    alt1.set_mtid(1);
    alt1.set_mpos(700);
    alt1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&alt1).unwrap();
    let mut alt2 = bam::Record::new();
    alt2.set(b"alt_chrom_pair", Some(&CigarString(vec![Cigar::Match(80)])), alt_seq2, &[30u8; 80]);
    alt2.set_tid(1);
    alt2.set_pos(700);
    alt2.set_mtid(1);
    alt2.set_mpos(500);
    alt2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&alt2).unwrap();

    // Several unaligned-template pairs (direct-emit path, consumes 2
    // records via a second alignments.next() call each).
    for i in 0..10 {
        let off1 = (i * 17) % usable_len;
        let off2 = (i * 17 + 20) % usable_len;
        let seq1 = &kir_bytes[off1..off1 + 70];
        let seq2 = &kir_bytes[off2..off2 + 70];
        let name = format!("unaligned_{i}");
        let mut um1 = bam::Record::new();
        um1.set(name.as_bytes(), None, seq1, &[25u8; 70]);
        um1.set_tid(-1);
        um1.set_pos(-1);
        um1.set_mtid(-1);
        um1.set_mpos(-1);
        um1.set_flags(0x1 | 0x4 | 0x8 | 0x40);
        writer.write(&um1).unwrap();
        let mut um2 = bam::Record::new();
        um2.set(name.as_bytes(), None, seq2, &[25u8; 70]);
        um2.set_tid(-1);
        um2.set_pos(-1);
        um2.set_mtid(-1);
        um2.set_mpos(-1);
        um2.set_flags(0x1 | 0x4 | 0x8 | 0x80);
        writer.write(&um2).unwrap();
    }

    drop(writer);
}

/// The core byte-identity requirement of the P1 parallelism task for BAM
/// input: `-t 1`, `-t 4`, and `-t 8` must all produce output byte-identical
/// to EACH OTHER and to the real oracle (always run at its own `-t 1` -- see
/// this module's docs for why that single oracle run is the correct
/// comparator). `pair_count` (3000) spans multiple `512 * threads`-sized
/// parallel pass-1 batches even at `threads = 8`.
#[test]
fn parallel_threads_1_4_8_byte_identical_to_each_other_and_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = dir.path().join("large.bam");
    build_large_test_bam(&bam_path, 3000);

    let coord_fasta = coord_fa();
    let t1_prefix = dir.path().join("t1");
    let t4_prefix = dir.path().join("t4");
    let t8_prefix = dir.path().join("t8");
    let oracle_prefix = dir.path().join("oracle");

    let m1 = run_rust_with_threads(&bam_path, &coord_fasta, &t1_prefix, 1);
    let m4 = run_rust_with_threads(&bam_path, &coord_fasta, &t4_prefix, 4);
    let m8 = run_rust_with_threads(&bam_path, &coord_fasta, &t8_prefix, 8);
    assert!(!m1.single_end && !m4.single_end && !m8.single_end, "fixture BAM must be paired");
    // Total emitted (pass1 + pass2) must be independent of threads too.
    assert_eq!(m1.pass1_emitted + m1.pass2_emitted, m4.pass1_emitted + m4.pass2_emitted);
    assert_eq!(m1.pass1_emitted + m1.pass2_emitted, m8.pass1_emitted + m8.pass2_emitted);

    run_oracle(&bam_path, &coord_fasta, &oracle_prefix);

    for suffix in ["_1.fq", "_2.fq"] {
        let t1_file = dir.path().join(format!("t1{suffix}"));
        assert_files_byte_identical(&t1_file, &dir.path().join(format!("t4{suffix}")));
        assert_files_byte_identical(&t1_file, &dir.path().join(format!("t8{suffix}")));
        assert_files_byte_identical(&t1_file, &dir.path().join(format!("oracle{suffix}")));
    }

    // Sanity: the fixture actually produced a non-trivial number of
    // candidates (not vacuously empty).
    let out1 = std::fs::read_to_string(dir.path().join("t1_1.fq")).unwrap();
    let emitted: usize = out1.lines().filter(|l| l.starts_with('@')).count();
    assert!(emitted > 1000, "expected a large number of candidate pairs, got {emitted}");
}
