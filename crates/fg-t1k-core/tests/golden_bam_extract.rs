//! Byte-golden test for `fg_t1k_core::bam_extract::extract_from_bam`,
//! converted from the retired `fg-t1k-sys` `diff_bam_extract.rs` FFI/subprocess
//! differential (see `tests/common/mod.rs`). Builds the SAME coordinate-sorted
//! BAMs programmatically (via `rust_htslib`), runs the Rust `bam_extract`
//! driver over them, and asserts the emitted FASTQ bytes match the committed
//! byte-goldens (the Rust port's output, byte-identical to the C++
//! `bam-extractor` oracle at `-t 1`). The Rust-only
//! `missing_unaligned_mate_is_rejected_by_the_port` test is kept verbatim.

mod common;

use common::assert_byte_golden;
use fg_t1k_core::alignments::Alignments;
use fg_t1k_core::bam_extract::{self, CoordRecord};
use fg_t1k_core::extract::{CandidateSink, ReadRecord, output_seq};
use fg_t1k_core::ref_kmer_filter::RefKmerFilter;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{self, Header, Writer};
use std::path::{Path, PathBuf};

const INITIAL_KMER_LENGTH: usize = 9;

fn fixture(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures").join(rel)
}

fn coord_fa() -> PathBuf {
    fixture("refbuild/golden/kir_rna_coord.fa")
}

const KIR2DL1_CHROM: &str = "chr19";
const KIR2DL1_START: i64 = 54_769_793;

fn kir2dl1_sequence() -> String {
    bam_extract::parse_coord_fa(&coord_fa())
        .expect("parse coord fa")
        .into_iter()
        .find(|r: &CoordRecord| r.name == "KIR2DL1*0010101")
        .expect("KIR2DL1*0010101 present")
        .seq
}

fn noise_seq(len: usize) -> Vec<u8> {
    let pattern = b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
    (0..len).map(|i| pattern[i % pattern.len()]).collect()
}

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

    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("writer");

    let mut ot1 = bam::Record::new();
    ot1.set(b"on_target_pair", Some(&CigarString(vec![Cigar::Match(100)])), &kir_bytes[0..100], &[30u8; 100]);
    ot1.set_tid(0);
    ot1.set_pos(KIR2DL1_START + 10);
    ot1.set_mtid(0);
    ot1.set_mpos(KIR2DL1_START + 210);
    ot1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&ot1).unwrap();
    let mut ot2 = bam::Record::new();
    ot2.set(b"on_target_pair", Some(&CigarString(vec![Cigar::Match(100)])), &kir_bytes[200..300], &[30u8; 100]);
    ot2.set_tid(0);
    ot2.set_pos(KIR2DL1_START + 210);
    ot2.set_mtid(0);
    ot2.set_mpos(KIR2DL1_START + 10);
    ot2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&ot2).unwrap();

    let off = noise_seq(80);
    let mut off1 = bam::Record::new();
    off1.set(b"off_target_pair", Some(&CigarString(vec![Cigar::Match(80)])), &off, &[30u8; 80]);
    off1.set_tid(0);
    off1.set_pos(1000);
    off1.set_mtid(0);
    off1.set_mpos(1200);
    off1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&off1).unwrap();
    let mut off2 = bam::Record::new();
    off2.set(b"off_target_pair", Some(&CigarString(vec![Cigar::Match(80)])), &off, &[30u8; 80]);
    off2.set_tid(0);
    off2.set_pos(1200);
    off2.set_mtid(0);
    off2.set_mpos(1000);
    off2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&off2).unwrap();

    let homopolymer = vec![b'A'; 90];
    let mut lc1 = bam::Record::new();
    lc1.set(b"low_complexity_pair", Some(&CigarString(vec![Cigar::Match(90)])), &homopolymer, &[30u8; 90]);
    lc1.set_tid(0);
    lc1.set_pos(KIR2DL1_START + 1000);
    lc1.set_mtid(0);
    lc1.set_mpos(KIR2DL1_START + 1200);
    lc1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&lc1).unwrap();
    let mut lc2 = bam::Record::new();
    lc2.set(b"low_complexity_pair", Some(&CigarString(vec![Cigar::Match(90)])), &homopolymer, &[30u8; 90]);
    lc2.set_tid(0);
    lc2.set_pos(KIR2DL1_START + 1200);
    lc2.set_mtid(0);
    lc2.set_mpos(KIR2DL1_START + 1000);
    lc2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&lc2).unwrap();

    let mut alt1 = bam::Record::new();
    alt1.set(b"alt_chrom_pair", Some(&CigarString(vec![Cigar::Match(80)])), &kir_bytes[400..480], &[30u8; 80]);
    alt1.set_tid(1);
    alt1.set_pos(500);
    alt1.set_mtid(1);
    alt1.set_mpos(700);
    alt1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&alt1).unwrap();
    let mut alt2 = bam::Record::new();
    alt2.set(b"alt_chrom_pair", Some(&CigarString(vec![Cigar::Match(80)])), &kir_bytes[600..680], &[30u8; 80]);
    alt2.set_tid(1);
    alt2.set_pos(700);
    alt2.set_mtid(1);
    alt2.set_mpos(500);
    alt2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&alt2).unwrap();

    let mut um1 = bam::Record::new();
    um1.set(b"unaligned_template_pair", None, &kir_bytes[800..880], &[25u8; 80]);
    um1.set_tid(-1);
    um1.set_pos(-1);
    um1.set_mtid(-1);
    um1.set_mpos(-1);
    um1.set_flags(0x1 | 0x4 | 0x8 | 0x40);
    writer.write(&um1).unwrap();
    let mut um2 = bam::Record::new();
    um2.set(b"unaligned_template_pair", None, &kir_bytes[900..980], &[25u8; 80]);
    um2.set_tid(-1);
    um2.set_pos(-1);
    um2.set_mtid(-1);
    um2.set_mpos(-1);
    um2.set_flags(0x1 | 0x4 | 0x8 | 0x80);
    writer.write(&um2).unwrap();

    let noise1 = noise_seq(70);
    let noise2 = noise_seq(70);
    let mut nu1 = bam::Record::new();
    nu1.set(b"unaligned_noise_pair", None, &noise1, &[25u8; 70]);
    nu1.set_tid(-1);
    nu1.set_pos(-1);
    nu1.set_mtid(-1);
    nu1.set_mpos(-1);
    nu1.set_flags(0x1 | 0x4 | 0x8 | 0x40);
    writer.write(&nu1).unwrap();
    let mut nu2 = bam::Record::new();
    nu2.set(b"unaligned_noise_pair", None, &noise2, &[25u8; 70]);
    nu2.set_tid(-1);
    nu2.set_pos(-1);
    nu2.set_mtid(-1);
    nu2.set_mpos(-1);
    nu2.set_flags(0x1 | 0x4 | 0x8 | 0x80);
    writer.write(&nu2).unwrap();

    drop(writer);
}

fn build_single_end_test_bam(path: &Path) {
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
    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("writer");

    let mut r1 = bam::Record::new();
    r1.set(b"single_on_target", Some(&CigarString(vec![Cigar::Match(100)])), &kir_bytes[0..100], &[30u8; 100]);
    r1.set_tid(0);
    r1.set_pos(KIR2DL1_START + 10);
    r1.set_mtid(-1);
    r1.set_mpos(-1);
    r1.set_flags(0);
    writer.write(&r1).unwrap();
    let mut r1_dup = bam::Record::new();
    r1_dup.set(b"single_on_target", Some(&CigarString(vec![Cigar::Match(100)])), &kir_bytes[0..100], &[30u8; 100]);
    r1_dup.set_tid(0);
    r1_dup.set_pos(KIR2DL1_START + 20);
    r1_dup.set_mtid(-1);
    r1_dup.set_mpos(-1);
    r1_dup.set_flags(0x100);
    writer.write(&r1_dup).unwrap();

    let off = noise_seq(80);
    let mut r2 = bam::Record::new();
    r2.set(b"single_off_target", Some(&CigarString(vec![Cigar::Match(80)])), &off, &[30u8; 80]);
    r2.set_tid(0);
    r2.set_pos(2000);
    r2.set_mtid(-1);
    r2.set_mpos(-1);
    r2.set_flags(0);
    writer.write(&r2).unwrap();

    let homopolymer = vec![b'T'; 90];
    let mut r3 = bam::Record::new();
    r3.set(b"single_low_complexity", Some(&CigarString(vec![Cigar::Match(90)])), &homopolymer, &[30u8; 90]);
    r3.set_tid(0);
    r3.set_pos(KIR2DL1_START + 500);
    r3.set_mtid(-1);
    r3.set_mpos(-1);
    r3.set_flags(0);
    writer.write(&r3).unwrap();

    drop(writer);
}

struct FastqFileSink {
    fp1: std::fs::File,
    fp2: Option<std::fs::File>,
}

impl FastqFileSink {
    fn create_paired(prefix: &Path) -> Self {
        Self {
            fp1: std::fs::File::create(format!("{}_1.fq", prefix.display())).unwrap(),
            fp2: Some(std::fs::File::create(format!("{}_2.fq", prefix.display())).unwrap()),
        }
    }
    fn create_single(prefix: &Path) -> Self {
        Self { fp1: std::fs::File::create(format!("{}.fq", prefix.display())).unwrap(), fp2: None }
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

fn run_rust(bam_path: &Path, coord_fasta: &Path, out_prefix: &Path) -> bam_extract::BamExtractMetrics {
    let coord_records = bam_extract::parse_coord_fa(coord_fasta).unwrap();
    let mut filter = RefKmerFilter::from_reference_fasta(coord_fasta, INITIAL_KMER_LENGTH).unwrap();
    let mut alignments = Alignments::open(bam_path).unwrap();
    let genes = bam_extract::build_genes(&alignments, &coord_records).unwrap();
    let single_end = alignments.general_info(true).unwrap().frag_stdev == 0;
    alignments.rewind().unwrap();
    let mut sink = if single_end {
        FastqFileSink::create_single(out_prefix)
    } else {
        FastqFileSink::create_paired(out_prefix)
    };
    bam_extract::extract_from_bam(&mut alignments, &mut filter, &genes, false, -1, &mut sink).unwrap()
}

#[test]
fn paired_synthetic_fixture_byte_golden() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = dir.path().join("test.bam");
    build_test_bam(&bam_path);
    let prefix = dir.path().join("rust");
    let metrics = run_rust(&bam_path, &coord_fa(), &prefix);
    assert!(!metrics.single_end);

    let out1 = std::fs::read(dir.path().join("rust_1.fq")).unwrap();
    let out2 = std::fs::read(dir.path().join("rust_2.fq")).unwrap();
    let emitted: Vec<&str> =
        std::str::from_utf8(&out1).unwrap().lines().filter_map(|l| l.strip_prefix('@')).collect();
    assert!(emitted.contains(&"on_target_pair"));
    assert!(emitted.contains(&"alt_chrom_pair"));
    assert!(emitted.contains(&"unaligned_template_pair"));
    assert!(!emitted.contains(&"off_target_pair"));
    assert!(!emitted.contains(&"low_complexity_pair"));
    assert!(!emitted.contains(&"unaligned_noise_pair"));

    assert_byte_golden("bam_extract/paired_1.fq", &out1);
    assert_byte_golden("bam_extract/paired_2.fq", &out2);
}

#[test]
fn single_end_synthetic_fixture_byte_golden() {
    let dir = tempfile::tempdir().unwrap();
    let bam_path = dir.path().join("single.bam");
    build_single_end_test_bam(&bam_path);
    let prefix = dir.path().join("rust");
    let metrics = run_rust(&bam_path, &coord_fa(), &prefix);
    assert!(metrics.single_end);

    let out = std::fs::read(dir.path().join("rust.fq")).unwrap();
    let emitted: Vec<&str> =
        std::str::from_utf8(&out).unwrap().lines().filter_map(|l| l.strip_prefix('@')).collect();
    assert_eq!(emitted.iter().filter(|&&id| id == "single_on_target").count(), 1, "usedName dedup");
    assert!(!emitted.contains(&"single_off_target"));
    assert!(!emitted.contains(&"single_low_complexity"));

    assert_byte_golden("bam_extract/single.fq", &out);
}

/// Rust-only regression (kept verbatim): a paired unmapped read whose mate
/// record is entirely absent must be rejected deterministically.
#[test]
fn missing_unaligned_mate_is_rejected_by_the_port() {
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
        r.set_flags(0x1 | 0x4 | 0x8 | 0x40);
        writer.write(&r).unwrap();
    }

    let coord_records = bam_extract::parse_coord_fa(&coord_fa()).unwrap();
    let mut filter = RefKmerFilter::from_reference_fasta(&coord_fa(), INITIAL_KMER_LENGTH).unwrap();
    let mut alignments = Alignments::open(&bam_path).unwrap();
    let genes = bam_extract::build_genes(&alignments, &coord_records).unwrap();
    let mut sink = FastqFileSink::create_paired(&dir.path().join("rust_out"));
    let result = bam_extract::extract_from_bam(&mut alignments, &mut filter, &genes, false, -1, &mut sink);
    assert!(result.is_err(), "must reject a missing unaligned mate");
    assert!(result.unwrap_err().to_string().contains("not showing up together"));
}
