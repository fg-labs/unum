//! Golden-file test for `fg_t1k_core::alignments::Alignments`, converted from
//! the retired `fg-t1k-sys` `diff_alignments.rs` FFI/subprocess differential
//! (see `tests/common/mod.rs`). Builds the SAME small coordinate-sorted BAMs
//! programmatically (via `rust_htslib::bam::Writer`), reads them with the Rust
//! `Alignments` port, and freezes each record's accessor fields + the global
//! `GeneralInfo` into the `alignments.txt` golden (the values were
//! byte-identical to the vendored C++ `Alignments` class when the oracle
//! existed). The reverse-complement record's hand-derived self-consistency
//! check is kept (its former oracle cross-check is dropped).

mod common;

use common::Golden;
use fg_t1k_core::alignments::Alignments;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{self, Header, Writer};
use std::path::Path;

#[allow(clippy::too_many_lines)]
fn build_test_bam(path: &Path) {
    let mut header = Header::new();

    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);

    let mut sq_chr19 = HeaderRecord::new(b"SQ");
    sq_chr19.push_tag(b"SN", "chr19");
    sq_chr19.push_tag(b"LN", 58_617_616);
    header.push_record(&sq_chr19);

    let mut sq_chr20 = HeaderRecord::new(b"SQ");
    sq_chr20.push_tag(b"SN", "chr20");
    sq_chr20.push_tag(b"LN", 64_444_167);
    header.push_record(&sq_chr20);

    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

    let mut r1 = bam::Record::new();
    let r1_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let r1_len = u32::try_from(r1_seq.len()).unwrap();
    r1.set(b"pair_fwd", Some(&CigarString(vec![Cigar::Match(r1_len)])), r1_seq, &vec![30u8; r1_seq.len()]);
    r1.set_tid(0);
    r1.set_pos(100);
    r1.set_mtid(0);
    r1.set_mpos(300);
    r1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&r1).expect("write r1");

    let mut r2 = bam::Record::new();
    let asym_seq = b"AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACC";
    let asym_len = u32::try_from(asym_seq.len()).unwrap();
    r2.set(
        b"pair_fwd",
        Some(&CigarString(vec![Cigar::Match(asym_len)])),
        asym_seq,
        &(0..asym_seq.len()).map(|i| 2 + u8::try_from(i % 40).unwrap()).collect::<Vec<u8>>(),
    );
    r2.set_tid(0);
    r2.set_pos(300);
    r2.set_mtid(0);
    r2.set_mpos(100);
    r2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&r2).expect("write r2");

    let mut r3 = bam::Record::new();
    let spliced_seq = b"TTGGCCAATTGGCCAATTGGCCAA";
    r3.set(
        b"spliced_read",
        Some(&CigarString(vec![Cigar::Match(12), Cigar::RefSkip(500), Cigar::Match(12)])),
        spliced_seq,
        &[25u8; 24],
    );
    r3.set_tid(0);
    r3.set_pos(1000);
    r3.set_mtid(-1);
    r3.set_mpos(-1);
    r3.set_flags(0x10);
    writer.write(&r3).expect("write r3");

    let mut r4 = bam::Record::new();
    let r4_seq = b"GATTACAGATTACAGATTACAGATTACAG";
    let r4_len = u32::try_from(r4_seq.len()).unwrap();
    r4.set(b"one_unmapped", Some(&CigarString(vec![Cigar::Match(r4_len)])), r4_seq, &vec![20u8; r4_seq.len()]);
    r4.set_tid(0);
    r4.set_pos(2000);
    r4.set_mtid(0);
    r4.set_mpos(2000);
    r4.set_flags(0x1 | 0x8 | 0x40);
    writer.write(&r4).expect("write r4");

    let mut r4_mate = bam::Record::new();
    let r4_mate_seq = b"CTGTAATCTGTAATCTGTAATCTGTAAT";
    r4_mate.set(b"one_unmapped", None, r4_mate_seq, &vec![20u8; r4_mate_seq.len()]);
    r4_mate.set_tid(0);
    r4_mate.set_pos(2000);
    r4_mate.set_mtid(0);
    r4_mate.set_mpos(2000);
    r4_mate.set_flags(0x1 | 0x4 | 0x80);
    writer.write(&r4_mate).expect("write r4_mate");

    let mut r5 = bam::Record::new();
    r5.set(b"both_unmapped", None, b"NNNNNNNNNNNNNNNNNNNN", &[2u8; 20]);
    r5.set_tid(-1);
    r5.set_pos(-1);
    r5.set_mtid(-1);
    r5.set_mpos(-1);
    r5.set_flags(0x1 | 0x4 | 0x8 | 0x40);
    writer.write(&r5).expect("write r5");

    let mut r6 = bam::Record::new();
    r6.set(b"both_unmapped", None, b"AAAAAAAAAAAAAAAAAAAA", &[3u8; 20]);
    r6.set_tid(-1);
    r6.set_pos(-1);
    r6.set_mtid(-1);
    r6.set_mpos(-1);
    r6.set_flags(0x1 | 0x4 | 0x8 | 0x80);
    writer.write(&r6).expect("write r6");

    let mut r7 = bam::Record::new();
    let three_exon_seq: Vec<u8> = (0..45).map(|i| b"ACGT"[i % 4]).collect();
    r7.set(
        b"three_exon",
        Some(&CigarString(vec![
            Cigar::Match(15),
            Cigar::RefSkip(200),
            Cigar::Match(15),
            Cigar::RefSkip(300),
            Cigar::Match(15),
        ])),
        &three_exon_seq,
        &vec![35u8; three_exon_seq.len()],
    );
    r7.set_tid(1);
    r7.set_pos(50_000);
    r7.set_mtid(-1);
    r7.set_mpos(-1);
    r7.set_flags(0);
    writer.write(&r7).expect("write r7");

    drop(writer);
}

/// Serializes the current record's every accessor field into one golden line.
fn serialize_record(a: &Alignments) -> String {
    let segs = a.segments().iter().map(|s| format!("{},{}", s.a, s.b)).collect::<Vec<_>>().join(";");
    format!(
        "seq={}|qual={}|id={}|first={}|rev={}|materev={}|aligned={}|tmplaligned={}|primary={}|chrom={}|segs={}",
        String::from_utf8_lossy(&a.read_seq()),
        a.qual().iter().map(u8::to_string).collect::<Vec<_>>().join(","),
        a.read_id(),
        u8::from(a.is_first_mate()),
        u8::from(a.is_reverse()),
        u8::from(a.is_mate_reverse()),
        u8::from(a.is_aligned()),
        u8::from(a.is_template_aligned()),
        u8::from(a.is_primary()),
        a.chrom_id(),
        segs,
    )
}

fn build_general_info_bam(path: &Path, n_pairs: i64, read_len: usize) {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr19");
    sq.push_tag(b"LN", 58_617_616);
    header.push_record(&sq);

    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("writer");
    let seq: Vec<u8> = (0..read_len).map(|i| b"ACGT"[i % 4]).collect();
    let cigar_len = u32::try_from(read_len).unwrap();
    let qual = vec![30u8; read_len];
    for i in 0..n_pairs {
        let pos1 = 1000 + i * 1000;
        let insert = 300 + i;
        let pos2 = pos1 + insert - i64::try_from(read_len).unwrap();

        let mut m1 = bam::Record::new();
        m1.set(format!("pair{i}").as_bytes(), Some(&CigarString(vec![Cigar::Match(cigar_len)])), &seq, &qual);
        m1.set_tid(0);
        m1.set_pos(pos1);
        m1.set_mtid(0);
        m1.set_mpos(pos2);
        m1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&m1).expect("write m1");

        let mut m2 = bam::Record::new();
        m2.set(format!("pair{i}").as_bytes(), Some(&CigarString(vec![Cigar::Match(cigar_len)])), &seq, &qual);
        m2.set_tid(0);
        m2.set_pos(pos2);
        m2.set_mtid(0);
        m2.set_mpos(pos1);
        m2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&m2).expect("write m2");
    }
    drop(writer);
}

fn build_single_end_bam(path: &Path) {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr19");
    sq.push_tag(b"LN", 58_617_616);
    header.push_record(&sq);

    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("writer");
    let seq: Vec<u8> = (0..75).map(|i| b"ACGT"[i % 4]).collect();
    let qual = vec![30u8; 75];
    for i in 0..20i64 {
        let mut r = bam::Record::new();
        r.set(format!("single{i}").as_bytes(), Some(&CigarString(vec![Cigar::Match(75)])), &seq, &qual);
        r.set_tid(0);
        r.set_pos(1000 + i * 200);
        r.set_mtid(-1);
        r.set_mpos(-1);
        r.set_flags(0);
        writer.write(&r).expect("write single-end record");
    }
    drop(writer);
}

#[test]
fn alignments_matches_golden() {
    let mut golden = Golden::open("alignments.txt");

    // Per-record fields across all eight coverage categories.
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("multi.bam");
    build_test_bam(&bam);
    let labels = [
        "r1_forward_first_mate",
        "r2_reverse_second_mate",
        "r3_spliced_reverse_unpaired",
        "r4_mapped_with_unmapped_mate",
        "r4_mate_unmapped_with_mapped_mate",
        "r5_both_unmapped_first_mate",
        "r6_both_unmapped_second_mate",
        "r7_three_exon_forward",
    ];
    let mut a = Alignments::open(&bam).expect("open");
    for label in labels {
        assert!(a.next().expect("next"), "expected a record for {label}");
        golden.record(format!("record/{label}"), serialize_record(&a));
    }
    assert!(!a.next().expect("next at end"), "unexpected extra records");

    // rewind() must return to the first record: record it again.
    a.rewind().expect("rewind");
    assert!(a.next().expect("next after rewind"));
    golden.record("rewind/first_record", serialize_record(&a));

    // General info: paired (frag_stdev > 0) and single-end (frag_stdev == 0).
    let gi_bam = dir.path().join("general_info.bam");
    build_general_info_bam(&gi_bam, 29, 50);
    let mut gi = Alignments::open(&gi_bam).expect("open general_info");
    let info = gi.general_info(false).expect("general_info");
    assert!(info.frag_stdev > 0, "paired dataset must exercise frag_stdev>0");
    golden.record("general_info/paired", format!("read_len={},frag_stdev={}", info.read_len, info.frag_stdev));

    let se_bam = dir.path().join("single_end.bam");
    build_single_end_bam(&se_bam);
    let mut se = Alignments::open(&se_bam).expect("open single_end");
    let se_info = se.general_info(false).expect("general_info single");
    assert_eq!(se_info.frag_stdev, 0, "single-end must yield frag_stdev==0");
    golden.record("general_info/single_end", format!("read_len={},frag_stdev={}", se_info.read_len, se_info.frag_stdev));

    golden.finish();
}

/// Directly proves the reverse-strand record (r2) genuinely exercises the
/// reverse-complement path (self-consistent; its former oracle cross-check is
/// dropped -- the RC output is now frozen in the golden above).
#[test]
fn reverse_strand_record_is_genuinely_reverse_complemented() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("multi.bam");
    build_test_bam(&bam);
    let mut a = Alignments::open(&bam).expect("open");

    assert!(a.next().unwrap()); // r1
    assert!(a.next().unwrap()); // r2
    assert!(a.is_reverse(), "expected r2 to be reverse-strand");

    let decoded = a.read_seq();
    let stored_forward = b"AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACC";
    assert_ne!(decoded.as_slice(), stored_forward.as_slice());
    assert_ne!(decoded, stored_forward.iter().rev().copied().collect::<Vec<u8>>());

    let mut expected = stored_forward.to_vec();
    expected.reverse();
    for b in &mut expected {
        *b = match *b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        };
    }
    assert_eq!(decoded, expected, "read_seq must be the true reverse complement");
}
