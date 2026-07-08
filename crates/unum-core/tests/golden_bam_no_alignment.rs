//! `no-alignment ≡ FASTQ` equivalence (Stage 2b Task 7): the headline
//! correctness proof for `--bam-mode no-alignment`.
//!
//! Builds the SAME reads two ways -- (a) FASTQ (paired or single-end files)
//! and (b) a BAM (coordinate-sorted OR grouped, `@HD GO:query`) -- and asserts
//! that running the FASTQ path
//! ([`unum_core::extract::extract_candidates`]) and the no-alignment BAM
//! paths ([`bam_extract::extract_from_bam_no_alignment`] /
//! [`bam_extract::extract_from_bam_no_alignment_grouped`]) over them into
//! `VecSink`s produces the SAME candidates. This includes a deliberately
//! UNALIGNED (`0x4`, `tid=-1`) but k-mer-matching read, proving
//! [`bam_extract::extract_from_bam_no_alignment`]'s pass-2 alignment-gate
//! bypass (Task 4) does not change the result FASTQ would have produced (the
//! FASTQ path never sees alignment information at all, so it is the natural
//! "no gate" reference).
//!
//! # Clean, homogeneous fixtures only
//!
//! `no-alignment ≡ FASTQ` provably holds only for CLEAN data: every template
//! either fully paired (both mates present, no orphans) or genuinely
//! single-end, and UNIFORM read length within a fixture --
//! `compute_hit_len_required_no_alignment`'s `sampled_len/(count*5)` bump is
//! sampling-order-independent only under uniform length, and the grouped
//! one-pass legitimately diverges from the coordinate 2-pass on mixed/orphan
//! input (brief-sanctioned per-mode behavior, not a bug -- see
//! `extract_from_bam_no_alignment_grouped`'s "Pairing rules" doc section).
//! Every fixture below is all-paired XOR all-single-end, never a mix, and
//! every read within a fixture shares one fixed length.
//!
//! # Comparison strategy (stated per test below too)
//!
//! - Coordinate 2-pass vs FASTQ: compared as a SET keyed by id (a `HashMap`
//!   id -> (seq1, seq2)) via [`as_set`], because pass 2 emits in
//!   PAIR-COMPLETION order (it walks the BAM a second time and emits whenever
//!   a candidate's second mate is found), which need not match FASTQ's input
//!   order.
//! - Grouped one-pass vs FASTQ: the BAM is built `@HD GO:query` with records
//!   in the SAME order as the FASTQ mate files (mate1, mate2 adjacent per
//!   template, templates in input order), so the grouped one-pass's
//!   emission order matches FASTQ's input order exactly -- compared as an
//!   ORDERED `Vec` via [`as_seq`].
//! - Single-end (both coordinate and grouped): single-end input never
//!   reaches pass 2 at all (a direct one-pass streaming emit in
//!   BAM-encounter order), so it too is compared as an ORDERED `Vec`.

use std::collections::HashMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{self, Header, Writer};
use unum_core::alignments::Alignments;
use unum_core::bam_extract;
use unum_core::extract::{self, CandidateSink, ReadRecord};
use unum_core::ref_kmer_filter::RefKmerFilter;

const INITIAL_KMER_LENGTH: usize = 9;
const DEFAULT_SIMILARITY: f64 = extract::DEFAULT_REF_SEQ_SIMILARITY;

/// The same 400bp base-balanced synthetic reference used by
/// `bam_extract.rs`'s internal `PARALLEL_TEST_REF` unit tests and
/// `golden_fastq_extract.rs`'s synthetic fixture -- long enough to carve
/// several distinct on-reference substrings that all clear the default
/// `hitLenRequired`/`-s` filter.
const REF_SEQ: &str = "ACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATGGTACGTACGTACGTACGGGCATTCATGGCATTCATGGCATTCATGACGTTAGCACGTTAGCACGTTAGCACGTTAGCTGACCATGTGACCATGTGACCATGTGACCATGGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATG";

/// Writes a single-contig reference FASTA from [`REF_SEQ`] for `-f`/
/// [`RefKmerFilter::from_reference_fasta`].
fn ref_fasta(dir: &Path) -> PathBuf {
    let path = dir.join("ref.fa");
    std::fs::write(&path, format!(">only\n{REF_SEQ}\n")).unwrap();
    path
}

/// Off-reference noise (mirrors the crate's other synthetic-BAM tests' own
/// `noise` helpers), distinct enough from [`REF_SEQ`]'s composition that it
/// never accidentally k-mer-matches.
fn noise(len: usize) -> Vec<u8> {
    let pattern =
        b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
    (0..len).map(|i| pattern[i % pattern.len()]).collect()
}

/// A [`CandidateSink`] that accumulates emitted pairs in memory, for
/// comparing FASTQ-path and BAM-no-alignment-path output directly (mirrors
/// `bam_extract.rs`'s own internal `VecSink` test helper).
#[derive(Default)]
struct VecSink {
    pairs: Vec<(ReadRecord, Option<ReadRecord>)>,
}

impl CandidateSink for VecSink {
    fn emit_pair(&mut self, r1: &ReadRecord, r2: Option<&ReadRecord>) -> anyhow::Result<()> {
        self.pairs.push((r1.clone(), r2.cloned()));
        Ok(())
    }
}

/// Reduces emitted pairs to a SET keyed by id -- order-independent, for
/// comparing the coordinate 2-pass (completion order) against FASTQ (input
/// order).
fn as_set(
    pairs: &[(ReadRecord, Option<ReadRecord>)],
) -> HashMap<String, (Vec<u8>, Option<Vec<u8>>)> {
    pairs
        .iter()
        .map(|(r1, r2)| (r1.id.clone(), (r1.seq.clone(), r2.as_ref().map(|r| r.seq.clone()))))
        .collect()
}

/// Reduces emitted pairs to an ORDERED sequence -- for comparing the grouped
/// one-pass (and single-end, either mode) against FASTQ, both of which
/// preserve input/encounter order exactly.
fn as_seq(pairs: &[(ReadRecord, Option<ReadRecord>)]) -> Vec<(String, Vec<u8>, Option<Vec<u8>>)> {
    pairs
        .iter()
        .map(|(r1, r2)| (r1.id.clone(), r1.seq.clone(), r2.as_ref().map(|r| r.seq.clone())))
        .collect()
}

/// Runs the FASTQ path ([`extract::extract_candidates`], `threads == 1`)
/// over `r1`/`r2` with a freshly-loaded filter, returning the emitted pairs.
fn run_fastq(fasta: &Path, r1: &Path, r2: Option<&Path>) -> Vec<(ReadRecord, Option<ReadRecord>)> {
    let mut filter = RefKmerFilter::from_reference_fasta(fasta, INITIAL_KMER_LENGTH).unwrap();
    let mut source = extract::open_source(r1, r2).unwrap();
    let mut sink = VecSink::default();
    extract::extract_candidates(&mut source, &mut filter, DEFAULT_SIMILARITY, &mut sink).unwrap();
    sink.pairs
}

/// Runs the coordinate/unsorted no-alignment 2-pass
/// ([`bam_extract::extract_from_bam_no_alignment`], `threads == 1`, default
/// `mate_id_len == -1`) over `bam_path` with a freshly-loaded filter,
/// returning the emitted pairs.
fn run_bam_coordinate(bam_path: &Path, fasta: &Path) -> Vec<(ReadRecord, Option<ReadRecord>)> {
    let mut filter = RefKmerFilter::from_reference_fasta(fasta, INITIAL_KMER_LENGTH).unwrap();
    let mut alignments = Alignments::open(bam_path).unwrap();
    let mut sink = VecSink::default();
    bam_extract::extract_from_bam_no_alignment(
        &mut alignments,
        &mut filter,
        DEFAULT_SIMILARITY,
        -1,
        1,
        &mut sink,
    )
    .unwrap();
    sink.pairs
}

/// Runs the grouped/name-sorted no-alignment one-pass
/// ([`bam_extract::extract_from_bam_no_alignment_grouped`], `threads == 1`,
/// default `mate_id_len == -1`) over `bam_path` with a freshly-loaded
/// filter, returning the emitted pairs.
fn run_bam_grouped(bam_path: &Path, fasta: &Path) -> Vec<(ReadRecord, Option<ReadRecord>)> {
    let mut filter = RefKmerFilter::from_reference_fasta(fasta, INITIAL_KMER_LENGTH).unwrap();
    let mut alignments = Alignments::open(bam_path).unwrap();
    let (_metrics, _single_end, sink) = bam_extract::extract_from_bam_no_alignment_grouped(
        &mut alignments,
        &mut filter,
        DEFAULT_SIMILARITY,
        -1,
        1,
        |_single_end| Ok::<VecSink, anyhow::Error>(VecSink::default()),
    )
    .unwrap();
    sink.pairs
}

/// Writes one FASTQ record with an all-`I` quality string (matching the
/// `[30u8; len]`/`[25u8; len]`-style fixed-quality convention the crate's
/// other synthetic-BAM tests use).
fn write_fastq_record(f: &mut std::fs::File, id: &str, seq: &[u8]) {
    let qual = "I".repeat(seq.len());
    writeln!(f, "@{id}\n{}\n+\n{qual}", std::str::from_utf8(seq).unwrap()).unwrap();
}

// -- Paired fixtures ---------------------------------------------------

const PAIR_LEN: usize = 90;

/// One paired-end template: UNIFORM `PAIR_LEN`-byte mate sequences,
/// all-paired (no orphans) -- see module docs on why the equivalence
/// requires exactly this shape.
struct PairTemplate {
    id: &'static str,
    seq1: Vec<u8>,
    seq2: Vec<u8>,
    /// Whether the BAM builders mark mate1 UNALIGNED (`0x4`, `tid=-1`;
    /// mate2 stays aligned, off-target). The FASTQ side never encodes
    /// alignment at all -- this field is the whole point of the
    /// gate-bypass proof: the BAM no-alignment path must include this
    /// template exactly as if this flag had never been set.
    mate1_unaligned: bool,
}

/// Four templates exercising every rule the paired no-alignment path must
/// agree with FASTQ on: OR-rescue (one mate hits, the other is noise),
/// both-hit, neither-hit (dropped), and the UNALIGNED-but-k-mer-matching
/// gate-bypass case.
fn paired_templates() -> Vec<PairTemplate> {
    let r = REF_SEQ.as_bytes();
    vec![
        PairTemplate {
            id: "hit_rescue_pair",
            seq1: r[0..PAIR_LEN].to_vec(),
            seq2: noise(PAIR_LEN),
            mate1_unaligned: false,
        },
        PairTemplate {
            id: "hit_both_pair",
            seq1: r[100..100 + PAIR_LEN].to_vec(),
            seq2: r[200..200 + PAIR_LEN].to_vec(),
            mate1_unaligned: false,
        },
        PairTemplate {
            id: "fail_pair",
            seq1: noise(PAIR_LEN),
            seq2: noise(PAIR_LEN),
            mate1_unaligned: false,
        },
        PairTemplate {
            id: "unaligned_gate_pair",
            seq1: r[300..300 + PAIR_LEN].to_vec(),
            seq2: noise(PAIR_LEN),
            mate1_unaligned: true,
        },
    ]
}

/// Writes `templates` as paired FASTQ files, mate1/mate2 in template order.
fn write_paired_fastq(dir: &Path, templates: &[PairTemplate]) -> (PathBuf, PathBuf) {
    let r1 = dir.join("r1.fq");
    let r2 = dir.join("r2.fq");
    let mut f1 = std::fs::File::create(&r1).unwrap();
    let mut f2 = std::fs::File::create(&r2).unwrap();
    for t in templates {
        write_fastq_record(&mut f1, t.id, &t.seq1);
        write_fastq_record(&mut f2, t.id, &t.seq2);
    }
    (r1, r2)
}

/// Writes `templates` as a paired BAM, mate1 immediately followed by mate2
/// per template (so a `GO:query`-tagged write is a valid grouped fixture),
/// templates in the SAME order given -- the order [`write_paired_fastq`]
/// itself iterates in, so a `GO:query` fixture built from the same
/// `templates` slice matches the FASTQ mate files' record order exactly.
/// `so_key`/`so_val` set the `@HD` sort-order-family tag (`(b"SO",
/// "coordinate")` or `(b"GO", "query")`).
fn write_paired_bam(path: &Path, templates: &[PairTemplate], so_key: &[u8], so_val: &str) {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(so_key, so_val);
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 1_000_000);
    header.push_record(&sq);

    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("bam writer");

    for (i, t) in templates.iter().enumerate() {
        let pos1 = 1_000 + i64::try_from(i).unwrap() * 1_000;
        let pos2 = pos1 + 200;
        let qual1 = vec![30u8; t.seq1.len()];
        let qual2 = vec![30u8; t.seq2.len()];

        let mut r1 = bam::Record::new();
        if t.mate1_unaligned {
            r1.set(t.id.as_bytes(), None, &t.seq1, &qual1);
            r1.set_tid(-1);
            r1.set_pos(-1);
            r1.set_mtid(0);
            r1.set_mpos(pos2);
            r1.set_flags(0x1 | 0x4 | 0x40);
        } else {
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
            // No `0x10`/`0x20` (reverse-strand) bits: `Alignments::read_seq`
            // reverse-complements a `0x10`-flagged record's SEQ (mirroring
            // real aligned-orientation BAMs), which would silently diverge
            // this raw-byte fixture from the FASTQ side unless the SEQ were
            // ALSO pre-reverse-complemented to compensate -- simpler to just
            // keep every record forward-strand, since no-alignment selection
            // never consults strand/position anyway.
            r1.set_flags(0x1 | 0x40);
        }
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
        if t.mate1_unaligned {
            r2.set_mtid(-1);
            r2.set_mpos(-1);
            r2.set_flags(0x1 | 0x8 | 0x80);
        } else {
            r2.set_mtid(0);
            r2.set_mpos(pos1);
            r2.set_flags(0x1 | 0x80); // forward-strand -- see r1's comment above.
        }
        writer.write(&r2).unwrap();
    }

    drop(writer);
}

#[test]
fn no_alignment_bam_equals_fastq_candidates_coordinate() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = ref_fasta(dir.path());
    let templates = paired_templates();
    let (r1, r2) = write_paired_fastq(dir.path(), &templates);
    let bam_path = dir.path().join("coord.bam");
    write_paired_bam(&bam_path, &templates, b"SO", "coordinate");

    let fastq_pairs = run_fastq(&fasta, &r1, Some(&r2));
    let bam_pairs = run_bam_coordinate(&bam_path, &fasta);

    // Sanity/teeth: the expected pass/fail mix actually happened on the
    // FASTQ side (proves the fixture and filter exercise real logic, not
    // just trivially agreeing on emptiness).
    let fastq_ids: Vec<&str> = fastq_pairs.iter().map(|(r1, _)| r1.id.as_str()).collect();
    assert!(fastq_ids.contains(&"hit_rescue_pair"));
    assert!(fastq_ids.contains(&"hit_both_pair"));
    assert!(fastq_ids.contains(&"unaligned_gate_pair"));
    assert!(!fastq_ids.contains(&"fail_pair"));

    // Coordinate 2-pass emits in pair-COMPLETION order -- compare as a SET
    // keyed by id (see module docs).
    assert_eq!(
        as_set(&fastq_pairs),
        as_set(&bam_pairs),
        "no-alignment coordinate 2-pass candidates diverged from the FASTQ path"
    );
}

#[test]
fn no_alignment_bam_equals_fastq_candidates_grouped() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = ref_fasta(dir.path());
    let templates = paired_templates();
    let (r1, r2) = write_paired_fastq(dir.path(), &templates);
    let bam_path = dir.path().join("grouped.bam");
    write_paired_bam(&bam_path, &templates, b"GO", "query");

    let fastq_pairs = run_fastq(&fasta, &r1, Some(&r2));
    let bam_pairs = run_bam_grouped(&bam_path, &fasta);

    // The grouped BAM is `@HD GO:query` with mate1/mate2 adjacent per
    // template, templates in the SAME order as the FASTQ mate files -- so
    // the grouped one-pass's emission order matches FASTQ's input order
    // exactly; compare as an ORDERED Vec (see module docs).
    assert_eq!(
        as_seq(&fastq_pairs),
        as_seq(&bam_pairs),
        "no-alignment grouped one-pass candidates diverged from the FASTQ path (order or content)"
    );
}

/// Dedicated, minimal proof of the pass-2 alignment-gate bypass: a
/// deliberately UNALIGNED (`0x4`, `tid=-1`) but k-mer-matching read must be
/// present in the no-alignment output, matching FASTQ (which never sees
/// alignment at all, so it is the natural "no gate" reference). A
/// non-matching pair is included alongside it to prove this isn't a
/// vacuous "everything passes" fixture.
#[test]
fn no_alignment_includes_kmer_passing_unaligned_reads() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = ref_fasta(dir.path());
    let templates = vec![
        PairTemplate {
            id: "unaligned_gate_pair",
            seq1: REF_SEQ.as_bytes()[300..300 + PAIR_LEN].to_vec(),
            seq2: noise(PAIR_LEN),
            mate1_unaligned: true,
        },
        PairTemplate {
            id: "fail_pair",
            seq1: noise(PAIR_LEN),
            seq2: noise(PAIR_LEN),
            mate1_unaligned: false,
        },
    ];
    let (r1, r2) = write_paired_fastq(dir.path(), &templates);

    let coord_bam_path = dir.path().join("gate_coord.bam");
    write_paired_bam(&coord_bam_path, &templates, b"SO", "coordinate");
    let grouped_bam_path = dir.path().join("gate_grouped.bam");
    write_paired_bam(&grouped_bam_path, &templates, b"GO", "query");

    let fastq_pairs = run_fastq(&fasta, &r1, Some(&r2));
    let coord_pairs = run_bam_coordinate(&coord_bam_path, &fasta);
    let grouped_pairs = run_bam_grouped(&grouped_bam_path, &fasta);

    let fastq_ids: Vec<&str> = fastq_pairs.iter().map(|(r1, _)| r1.id.as_str()).collect();
    let coord_ids: Vec<&str> = coord_pairs.iter().map(|(r1, _)| r1.id.as_str()).collect();
    let grouped_ids: Vec<&str> = grouped_pairs.iter().map(|(r1, _)| r1.id.as_str()).collect();

    assert!(
        fastq_ids.contains(&"unaligned_gate_pair"),
        "FASTQ never sees alignment -- must include it"
    );
    assert!(
        coord_ids.contains(&"unaligned_gate_pair"),
        "no-alignment coordinate 2-pass must include the k-mer-passing UNALIGNED read -- proves \
         run_pass2's !is_template_aligned() gate is bypassed under Selection::NoAlignment"
    );
    assert!(
        grouped_ids.contains(&"unaligned_gate_pair"),
        "no-alignment grouped one-pass must include the k-mer-passing UNALIGNED read too (it never \
         inspects alignment status at all)"
    );
    assert!(!fastq_ids.contains(&"fail_pair"));
    assert!(!coord_ids.contains(&"fail_pair"));
    assert!(!grouped_ids.contains(&"fail_pair"));

    // Coordinate 2-pass: SET compare (completion order). Grouped: ORDERED
    // Vec compare (matches FASTQ input order -- see module docs).
    assert_eq!(as_set(&fastq_pairs), as_set(&coord_pairs));
    assert_eq!(as_seq(&fastq_pairs), as_seq(&grouped_pairs));
}

// -- Single-end fixtures -------------------------------------------------

const SINGLE_LEN: usize = 80;

/// One single-end template: a UNIFORM `SINGLE_LEN`-byte read, no mate at
/// all (genuinely single-end, not an orphan).
struct SingleTemplate {
    id: &'static str,
    seq: Vec<u8>,
}

/// Three single-end templates: two on-reference (must be emitted) and one
/// noise (must be dropped).
fn single_templates() -> Vec<SingleTemplate> {
    let r = REF_SEQ.as_bytes();
    vec![
        SingleTemplate { id: "single_hit_a", seq: r[0..SINGLE_LEN].to_vec() },
        SingleTemplate { id: "single_fail", seq: noise(SINGLE_LEN) },
        SingleTemplate { id: "single_hit_b", seq: r[150..150 + SINGLE_LEN].to_vec() },
    ]
}

/// Writes `templates` as a single-end FASTQ file, in template order.
fn write_single_fastq(dir: &Path, templates: &[SingleTemplate]) -> PathBuf {
    let path = dir.join("single.fq");
    let mut f = std::fs::File::create(&path).unwrap();
    for t in templates {
        write_fastq_record(&mut f, t.id, &t.seq);
    }
    path
}

/// Writes `templates` as a single-end BAM (`0x1` UNSET on every record, no
/// mate fields), in template order. `so_key`/`so_val` set the `@HD`
/// sort-order-family tag, as in [`write_paired_bam`].
fn write_single_bam(path: &Path, templates: &[SingleTemplate], so_key: &[u8], so_val: &str) {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(so_key, so_val);
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", 1_000_000);
    header.push_record(&sq);

    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("bam writer");
    for (i, t) in templates.iter().enumerate() {
        let mut r = bam::Record::new();
        r.set(
            t.id.as_bytes(),
            Some(&CigarString(vec![Cigar::Match(u32::try_from(t.seq.len()).unwrap())])),
            &t.seq,
            &vec![30u8; t.seq.len()],
        );
        r.set_tid(0);
        r.set_pos(1_000 + i64::try_from(i).unwrap() * 500);
        r.set_mtid(-1);
        r.set_mpos(-1);
        r.set_flags(0);
        writer.write(&r).unwrap();
    }
    drop(writer);
}

#[test]
fn no_alignment_bam_equals_fastq_candidates_single_end_coordinate() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = ref_fasta(dir.path());
    let templates = single_templates();
    let r1 = write_single_fastq(dir.path(), &templates);
    let bam_path = dir.path().join("single_coord.bam");
    write_single_bam(&bam_path, &templates, b"SO", "coordinate");

    let fastq_pairs = run_fastq(&fasta, &r1, None);
    let bam_pairs = run_bam_coordinate(&bam_path, &fasta);

    assert!(fastq_pairs.iter().any(|(r, _)| r.id == "single_hit_a"));
    assert!(fastq_pairs.iter().any(|(r, _)| r.id == "single_hit_b"));
    assert!(!fastq_pairs.iter().any(|(r, _)| r.id == "single_fail"));

    // Single-end input skips pass 2 entirely on the coordinate path
    // (`extract_from_bam_no_alignment`'s doc comment) -- a direct pass-1
    // emit in BAM-encounter order, matching FASTQ's input order here since
    // both were built from the same `single_templates()` order. Compare as
    // an ORDERED Vec.
    assert_eq!(
        as_seq(&fastq_pairs),
        as_seq(&bam_pairs),
        "no-alignment coordinate single-end candidates diverged from the FASTQ path"
    );
}

#[test]
fn no_alignment_bam_equals_fastq_candidates_single_end_grouped() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = ref_fasta(dir.path());
    let templates = single_templates();
    let r1 = write_single_fastq(dir.path(), &templates);
    let bam_path = dir.path().join("single_grouped.bam");
    write_single_bam(&bam_path, &templates, b"GO", "query");

    let fastq_pairs = run_fastq(&fasta, &r1, None);
    let bam_pairs = run_bam_grouped(&bam_path, &fasta);

    // The grouped one-pass's genuine single-end (`0x1` UNSET) lone-emit
    // rule streams and flushes group-by-group in encounter order, matching
    // FASTQ's input order. Compare as an ORDERED Vec.
    assert_eq!(
        as_seq(&fastq_pairs),
        as_seq(&bam_pairs),
        "no-alignment grouped single-end candidates diverged from the FASTQ path"
    );
}
