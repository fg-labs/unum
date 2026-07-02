#![cfg(feature = "t1k-sys")]
//! Differential test: builds a Rust `RefKmerFilter` and a real C++ `SeqSet`
//! (via the opaque-handle FFI shim) from the SAME T1K-style reference FASTA,
//! then asserts `RefKmerFilter::get_overlaps_from_read` (Task 4b's port of
//! `SeqSet::GetOverlapsFromRead`, `SeqSet.hpp:1594-1912`) produces IDENTICAL
//! overlap lists to the real C++ oracle -- same count, and for each overlap
//! (in order): `seqIdx`/`strand`/`readStart`/`readEnd`/`seqStart`/`seqEnd`/
//! `matchCnt` exactly, and `similarity` compared with EXACT `f64` equality
//! (NOT an epsilon -- `similarity` is a plain, deterministic
//! `matchCnt / (seqSpan + readSpan)` ratio; see `ref_kmer_filter.rs`'s
//! `get_overlaps_from_read` doc comment for the FLOATS rule this enforces).
//!
//! # Adversarial coverage (learning from Task 4a's `diff_align_algo.rs`)
//!
//! Exact substrings alone would never exercise `AlignAlgo`'s gap-alignment
//! path inside `GetOverlapsFromRead` (SeqSet.hpp:1716-1733,1802-1818) at
//! all -- every consecutive hit pair would be perfectly colinear with no
//! gap, so the direct `matchCnt += 2*(a-aPrev)` branch would be the ONLY
//! branch ever taken. This test battery deliberately also covers: reads
//! with 1..K point mismatches at varied densities/positions (fragments the
//! hit chain, forcing gap alignment AND exercising the mismatch-counting
//! arithmetic); reads with small INSERTIONS and DELETIONS relative to the
//! reference (exercises the indel-counting arithmetic AND `AlignAlgo`'s
//! actual gap-open/gap-extend traceback, not just mismatch scoring -- this
//! is where Task 4a's own bugs hid); reads spanning widely-separated hits
//! (multi-gap, exercising the non-colinear hitCoords branches at
//! SeqSet.hpp:1770-1831); reads matching multiple alleles (both fixtures
//! have several highly similar KIR/HLA paralogs); reverse-complement reads;
//! near-`hitLenRequired`-threshold reads; reads expected to end up entirely
//! filtered out (similarity below threshold, or no overlap surviving at
//! all); and a large batch of seeded pseudo-random mutations sweeping
//! mismatch/indel counts.

use fg_t1k_core::ref_kmer_filter::{RefKmerFilter, Scratch};
use fg_t1k_sys::CppSeqSet;
use std::path::{Path, PathBuf};

/// `SeqSet::InferKmerLength()` evaluates to `8` for both `kir_rna_seq.fa`
/// (total ref length 8781) and `hla_rna_seq.fa` (12260), neither exceeding
/// `FastqExtractor.cpp`'s literal initial default of `9` -- so stock never
/// rebuilds via `UpdateKmerLength` for either fixture and ends up using `9`
/// unmodified. Matches `diff_refkmerfilter.rs`'s identical reasoning/value.
const KMER_LENGTH: usize = 9;

fn fixture_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/refbuild/golden").join(rel)
}

/// A loaded (Rust, C++) pair built from the same reference FASTA at
/// [`KMER_LENGTH`] with ctor-default `hitLenRequired=31`/
/// `refSeqSimilarity=0.8` on both sides (`RefKmerFilter::from_reference_fasta`
/// sets these defaults itself; `SeqSet(kmerLength)`'s constructor defaults
/// match -- see `SeqSet.hpp:763-768`), plus the raw per-record sequences
/// (for pulling exact substrings / building mutated reads).
struct FilterPair {
    rust: RefKmerFilter,
    cpp: CppSeqSet,
    records: Vec<(String, Vec<u8>)>,
}

fn load_pair(fasta_rel: &str) -> FilterPair {
    let path = fixture_path(fasta_rel);
    let rust = RefKmerFilter::from_reference_fasta(&path, KMER_LENGTH)
        .unwrap_or_else(|e| panic!("RefKmerFilter::from_reference_fasta({fasta_rel}): {e}"));

    let mut cpp = CppSeqSet::new(i32::try_from(KMER_LENGTH).unwrap());
    cpp.load_ref(&path);

    let records = parse_fasta_for_test(&path);
    assert_eq!(
        records.len(),
        rust.seq_count(),
        "test's own FASTA parse and RefKmerFilter disagree on record count for {fasta_rel}"
    );

    FilterPair { rust, cpp, records }
}

/// Trivial, independent FASTA parser -- deliberately NOT shared code with
/// `fg_t1k_core::ref_kmer_filter`'s own (private) FASTA parser, matching
/// `diff_refkmerfilter.rs`'s identical helper.
fn parse_fasta_for_test(path: &Path) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                records.push((id, std::mem::take(&mut current_seq)));
            }
            current_id = Some(rest.split_whitespace().next().unwrap_or("").to_string());
        } else {
            current_seq.extend_from_slice(line.trim_end().as_bytes());
        }
    }
    if let Some(id) = current_id {
        records.push((id, current_seq));
    }
    records
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            other => panic!("reverse_complement: unsupported base {other}"),
        })
        .collect()
}

/// Deterministic (seeded, reproducible) pseudo-random generator, matching
/// `diff_refkmerfilter.rs`'s identical xorshift64 helper (kept as a separate
/// copy rather than a shared dev-dependency, per this test crate's existing
/// convention of small, self-contained test files).
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        Self { state: seed | 1 }
    }

    fn next_u64(&mut self) -> u64 {
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }

    fn next_base(&mut self) -> u8 {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        BASES[(self.next_u64() % 4) as usize]
    }

    fn next_base_excluding(&mut self, exclude: u8) -> u8 {
        loop {
            let b = self.next_base();
            if b != exclude {
                return b;
            }
        }
    }

    fn next_range(&mut self, bound: usize) -> usize {
        usize::try_from(self.next_u64() % bound as u64).unwrap_or(usize::MAX)
    }
}

fn pseudo_random_acgt(seed: u64, len: usize) -> Vec<u8> {
    let mut rng = Xorshift64::new(seed);
    (0..len).map(|_| rng.next_base()).collect()
}

/// Introduces exactly `k` point substitutions at deterministic (seeded),
/// distinct positions, each guaranteed to actually change the base.
fn mutate_k_positions(seq: &[u8], k: usize, seed: u64) -> Vec<u8> {
    let mut rng = Xorshift64::new(seed);
    let mut out = seq.to_vec();
    let mut positions: Vec<usize> = (0..seq.len()).collect();
    for i in 0..k.min(positions.len()) {
        let j = i + rng.next_range(positions.len() - i);
        positions.swap(i, j);
    }
    for &pos in positions.iter().take(k.min(positions.len())) {
        out[pos] = rng.next_base_excluding(out[pos]);
    }
    out
}

/// Deletes `del_len` bases starting at `pos` (a small deletion relative to
/// the reference -- the resulting read is SHORTER than the reference
/// substring it was derived from).
fn delete_span(seq: &[u8], pos: usize, del_len: usize) -> Vec<u8> {
    let mut out = seq.to_vec();
    let end = (pos + del_len).min(out.len());
    out.drain(pos..end);
    out
}

/// Inserts `ins_len` pseudo-random bases at `pos` (a small insertion
/// relative to the reference -- the resulting read is LONGER than the
/// reference substring it was derived from).
fn insert_span(seq: &[u8], pos: usize, ins_len: usize, seed: u64) -> Vec<u8> {
    let mut out = seq[..pos].to_vec();
    out.extend(pseudo_random_acgt(seed, ins_len));
    out.extend_from_slice(&seq[pos..]);
    out
}

/// Asserts the Rust and C++ `GetOverlapsFromRead` results are IDENTICAL:
/// same overlap count, and every field of every overlap (in order) matches
/// exactly -- including bit-exact `f64` equality for `similarity`, per this
/// module's FLOATS rule.
fn assert_overlaps_agree(label: &str, pair: &FilterPair, read: &[u8]) {
    let mut scratch = Scratch::default();
    let rust_overlaps = pair.rust.get_overlaps_from_read(read, &mut scratch);
    let cpp_overlaps = pair.cpp.get_overlaps_from_read(read);

    assert_eq!(
        rust_overlaps.is_some(),
        cpp_overlaps.is_some(),
        "{label}: Some/None mismatch (rust_is_some={}, cpp_is_some={}, read={read:?})",
        rust_overlaps.is_some(),
        cpp_overlaps.is_some(),
    );

    match (&rust_overlaps, &cpp_overlaps) {
        (None, None) => (),
        (None, Some(_)) | (Some(_), None) => unreachable!("just asserted Some/None agree above"),
        (Some(rust), Some(cpp)) => {
            assert_eq!(
                rust.len(),
                cpp.len(),
                "{label}: overlap count mismatch (rust={}, cpp={}, read={read:?})\nrust={rust:?}\ncpp={cpp:?}",
                rust.len(),
                cpp.len(),
            );
            for (i, (r, c)) in rust.iter().zip(cpp.iter()).enumerate() {
                let c_seq_idx = u32::try_from(c.seq_idx).unwrap_or_else(|_| {
                    panic!("{label}: overlap[{i}] cpp seq_idx {} is negative", c.seq_idx)
                });
                assert_eq!(r.seq_idx, c_seq_idx, "{label}: overlap[{i}] seq_idx mismatch");
                assert_eq!(i32::from(r.strand), c.strand, "{label}: overlap[{i}] strand mismatch");
                assert_eq!(r.read_start, c.read_start, "{label}: overlap[{i}] read_start mismatch");
                assert_eq!(r.read_end, c.read_end, "{label}: overlap[{i}] read_end mismatch");
                assert_eq!(r.seq_start, c.seq_start, "{label}: overlap[{i}] seq_start mismatch");
                assert_eq!(r.seq_end, c.seq_end, "{label}: overlap[{i}] seq_end mismatch");
                assert_eq!(r.match_cnt, c.match_cnt, "{label}: overlap[{i}] match_cnt mismatch");
                // Exact f64 equality by design -- see module docs.
                #[allow(clippy::float_cmp)]
                let similarity_matches = r.similarity == c.similarity;
                assert!(
                    similarity_matches,
                    "{label}: overlap[{i}] similarity mismatch (rust={}, cpp={}, bits rust={:#x} cpp={:#x})",
                    r.similarity,
                    c.similarity,
                    r.similarity.to_bits(),
                    c.similarity.to_bits(),
                );
            }
        }
    }
}

// ---- exact substrings --------------------------------------------------

#[test]
fn kir_exact_substrings_agree() {
    let pair = load_pair("kir_rna_seq.fa");
    let cases: &[(usize, usize, usize)] = &[
        (0, 0, 60),
        (0, 200, 100),
        (1, 500, 150),
        (2, 50, 80),
        (3, 900, 120),
        (5, 10, 200),
        (6, 300, 45),
    ];
    for &(rec_idx, start, len) in cases {
        let (id, seq) = &pair.records[rec_idx];
        assert!(start + len <= seq.len(), "case out of bounds for {id}");
        let read = &seq[start..start + len];
        assert_overlaps_agree(&format!("kir substr {id}@{start}+{len}"), &pair, read);

        let rc = reverse_complement(read);
        assert_overlaps_agree(&format!("kir substr RC {id}@{start}+{len}"), &pair, &rc);
    }
}

#[test]
fn hla_exact_substrings_agree() {
    let pair = load_pair("hla_rna_seq.fa");
    let cases: &[(usize, usize, usize)] =
        &[(0, 0, 70), (1, 300, 90), (4, 700, 110), (8, 150, 60), (11, 700, 130)];
    for &(rec_idx, start, len) in cases {
        let (id, seq) = &pair.records[rec_idx];
        assert!(start + len <= seq.len(), "case out of bounds for {id}");
        let read = &seq[start..start + len];
        assert_overlaps_agree(&format!("hla substr {id}@{start}+{len}"), &pair, read);
    }
}

// ---- point mismatches: 1..K sweep ---------------------------------------

#[test]
fn kir_mismatch_sweep_agrees() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[2];
    let start = 100;
    let len = 150;
    let base_read = &seq[start..start + len];

    for k in 0..=40usize {
        let read = mutate_k_positions(base_read, k, 1000 + k as u64);
        assert_overlaps_agree(&format!("{id}@{start}+{len} mutated k={k}"), &pair, &read);
    }
}

#[test]
fn hla_mismatch_sweep_agrees() {
    let pair = load_pair("hla_rna_seq.fa");
    let (id, seq) = &pair.records[4];
    let start = 200;
    let len = 100;
    let base_read = &seq[start..start + len];

    for k in 0..=30usize {
        let read = mutate_k_positions(base_read, k, 2000 + k as u64);
        assert_overlaps_agree(&format!("{id}@{start}+{len} mutated k={k}"), &pair, &read);
    }
}

// ---- small insertions and deletions (AlignAlgo gap path) ----------------

#[test]
fn kir_small_deletions_agree() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[3];
    let start = 200;
    let len = 200;
    let base_read = &seq[start..start + len];

    for (del_pos, del_len) in [(50, 1), (50, 2), (50, 3), (100, 5), (150, 8), (30, 1), (170, 4)] {
        let read = delete_span(base_read, del_pos, del_len);
        assert_overlaps_agree(
            &format!("{id}@{start}+{len} deletion pos={del_pos} len={del_len}"),
            &pair,
            &read,
        );
    }
}

#[test]
fn kir_small_insertions_agree() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[3];
    let start = 200;
    let len = 200;
    let base_read = &seq[start..start + len];

    for (ins_pos, ins_len, seed) in
        [(50, 1, 11), (50, 2, 12), (50, 3, 13), (100, 5, 14), (150, 8, 15), (30, 1, 16)]
    {
        let read = insert_span(base_read, ins_pos, ins_len, seed);
        assert_overlaps_agree(
            &format!("{id}@{start}+{len} insertion pos={ins_pos} len={ins_len}"),
            &pair,
            &read,
        );
    }
}

#[test]
fn hla_small_indels_agree() {
    let pair = load_pair("hla_rna_seq.fa");
    let (id, seq) = &pair.records[6];
    let start = 150;
    let len = 180;
    let base_read = &seq[start..start + len];

    for (op, pos, l, seed) in [
        ("del", 40, 2usize, 0u64),
        ("del", 90, 4, 0),
        ("ins", 40, 2, 21),
        ("ins", 90, 4, 22),
        ("del", 20, 1, 0),
        ("ins", 160, 3, 23),
    ] {
        let read = if op == "del" {
            delete_span(base_read, pos, l)
        } else {
            insert_span(base_read, pos, l, seed)
        };
        assert_overlaps_agree(&format!("{id}@{start}+{len} {op} pos={pos} len={l}"), &pair, &read);
    }
}

#[test]
fn combined_mismatches_and_indels_agree() {
    // Exercises both the mismatch AND indel arithmetic in the same read:
    // a deletion followed by point mutations in the remaining sequence.
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[0];
    let start = 300;
    let len = 220;
    let base_read = &seq[start..start + len];

    let deleted = delete_span(base_read, 60, 3);
    let mutated = mutate_k_positions(&deleted, 5, 555);
    assert_overlaps_agree(&format!("{id}@{start}+{len} del+mutate"), &pair, &mutated);

    let inserted = insert_span(base_read, 120, 4, 777);
    let mutated2 = mutate_k_positions(&inserted, 4, 888);
    assert_overlaps_agree(&format!("{id}@{start}+{len} ins+mutate"), &pair, &mutated2);
}

// ---- multi-gap / widely-separated hits -----------------------------------

#[test]
fn widely_separated_hits_multi_gap_agrees() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id0, seq0) = &pair.records[0];
    let (id1, seq1) = &pair.records[1];

    // Two 75bp fragments, same sequence, far apart -- noncolinear on the
    // reference, exercising GetOverlapsFromRead's non-colinear hitCoords
    // branches (and possibly splitting into multiple overlaps/LIS chains).
    let frag_a = &seq0[50..125];
    let frag_b = &seq0[900..975];
    let mut same_seq_concat = frag_a.to_vec();
    same_seq_concat.extend_from_slice(frag_b);
    assert_overlaps_agree(
        &format!("{id0} noncolinear same-seq concat (50..125 + 900..975)"),
        &pair,
        &same_seq_concat,
    );

    let frag_c = &seq1[400..475];
    let mut cross_seq_concat = frag_a.to_vec();
    cross_seq_concat.extend_from_slice(frag_c);
    assert_overlaps_agree(
        &format!("{id0}/{id1} noncolinear cross-seq concat"),
        &pair,
        &cross_seq_concat,
    );

    // Small gap (simulated small deletion at the fragment boundary): should
    // still chain via the `adjustRadius` handling in most cases.
    let mut near_adjacent = seq0[200..280].to_vec();
    near_adjacent.extend_from_slice(&seq0[285..365]);
    assert_overlaps_agree(&format!("{id0} near-adjacent (5bp gap)"), &pair, &near_adjacent);

    // Three fragments concatenated -- multiple internal gaps in one read.
    let mut three_frags = seq0[10..90].to_vec();
    three_frags.extend_from_slice(&seq0[500..580]);
    three_frags.extend_from_slice(&seq0[950..1030.min(seq0.len())]);
    assert_overlaps_agree(&format!("{id0} three-fragment concat"), &pair, &three_frags);
}

// ---- multiple alleles / reads that may match >1 reference sequence ------

#[test]
fn reads_from_multiple_records_each_agree() {
    // A representative substring from EACH loaded record, in one batch --
    // if any allele-specific edge case in the port diverges (e.g. an
    // off-by-one only visible for a specific consensus length), this
    // sweeps across every record in the fixture rather than just a few.
    let pair = load_pair("kir_rna_seq.fa");
    for (idx, (id, seq)) in pair.records.iter().enumerate() {
        if seq.len() < 100 {
            continue;
        }
        let start = seq.len() / 4;
        let len = 70.min(seq.len() - start);
        let read = &seq[start..start + len];
        assert_overlaps_agree(&format!("record[{idx}]={id}@{start}+{len}"), &pair, read);
    }
}

// ---- reverse-complement reads ---------------------------------------------

#[test]
fn rc_mutated_and_indel_reads_agree() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[4];

    let base = &seq[80..200];
    let mutated = mutate_k_positions(base, 8, 555);
    let rc_mutated = reverse_complement(&mutated);
    assert_overlaps_agree(&format!("{id} RC mutated k=8"), &pair, &rc_mutated);

    let deleted = delete_span(base, 40, 3);
    let rc_deleted = reverse_complement(&deleted);
    assert_overlaps_agree(&format!("{id} RC deletion"), &pair, &rc_deleted);

    let inserted = insert_span(base, 40, 3, 999);
    let rc_inserted = reverse_complement(&inserted);
    assert_overlaps_agree(&format!("{id} RC insertion"), &pair, &rc_inserted);

    let frag_a = seq[10..90].to_vec();
    let frag_b = seq[600..680].to_vec();
    let mut concat = frag_a;
    concat.extend_from_slice(&frag_b);
    let rc_concat = reverse_complement(&concat);
    assert_overlaps_agree(&format!("{id} RC noncolinear concat"), &pair, &rc_concat);
}

// ---- near-hitLenRequired-threshold reads ---------------------------------

#[test]
fn boundary_length_reads_agree() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[3];

    for len in [35usize, 40, 44, 45, 46, 50, 60] {
        let start = 300;
        let read = &seq[start..start + len];
        assert_overlaps_agree(&format!("{id}@{start}+{len} boundary length"), &pair, read);
    }

    // Exactly kmer_length, and kmer_length+1: minimal possible reads.
    let exactly_k = &seq[50..50 + KMER_LENGTH];
    assert_overlaps_agree(&format!("{id} len=kmer_length exactly"), &pair, exactly_k);
    let just_over = &seq[50..=50 + KMER_LENGTH];
    assert_overlaps_agree(&format!("{id} len=kmer_length+1"), &pair, just_over);

    // Shorter than kmer_length: both sides must return None/-1.
    let too_short = &seq[0..5];
    assert_overlaps_agree(&format!("{id} len=5 (< kmer_length)"), &pair, too_short);
}

// ---- reads expected to be filtered out entirely --------------------------

#[test]
fn heavily_mutated_reads_may_be_filtered_and_still_agree() {
    // At a high mutation density, similarity should drop below
    // refSeqSimilarity=0.8, causing the final filter to drop the overlap
    // entirely (empty Vec on both sides) -- or, depending on how the hit
    // chain fragments, `GetOverlapsFromHits` may find no chain long enough
    // to produce ANY overlap. Either way, both sides must agree.
    let pair = load_pair("hla_rna_seq.fa");
    let (id, seq) = &pair.records[2];
    let start = 300;
    let len = 120;
    let base_read = &seq[start..start + len];

    let mut saw_empty = false;
    for k in [20usize, 30, 40, 50, 60] {
        let read = mutate_k_positions(base_read, k.min(len), 3000 + k as u64);
        assert_overlaps_agree(&format!("{id}@{start}+{len} heavy mutation k={k}"), &pair, &read);
        let mut scratch = Scratch::default();
        if let Some(overlaps) = pair.rust.get_overlaps_from_read(&read, &mut scratch) {
            if overlaps.is_empty() {
                saw_empty = true;
            }
        }
    }
    assert!(
        saw_empty,
        "expected at least one heavily-mutated read to produce zero surviving overlaps \
         (i.e. actually exercise the final similarity filter/no-chain-found path)"
    );
}

#[test]
fn pure_noise_reads_produce_no_overlaps_and_agree() {
    let pair = load_pair("kir_rna_seq.fa");
    for (seed, len) in [(1u64, 80), (2, 120), (3, 60), (4, 150), (5, 40)] {
        let read = pseudo_random_acgt(seed, len);
        assert_overlaps_agree(&format!("random seed={seed} len={len}"), &pair, &read);
    }
}

// ---- large seeded pseudo-random batch -------------------------------------

#[test]
fn large_random_read_batch_agrees() {
    let pair = load_pair("hla_rna_seq.fa");
    let seq_count = pair.records.len();

    let mut checked = 0usize;
    for seed in 0u64..150 {
        let mut rng = Xorshift64::new(seed ^ 0x5EED_C0DE);
        let len = 40 + rng.next_range(161); // 40..=200

        let category = rng.next_range(5);
        let read: Vec<u8> = match category {
            0 => pseudo_random_acgt(seed.wrapping_mul(31) + 7, len),
            1 => {
                let rec_idx = rng.next_range(seq_count);
                let (_, seq) = &pair.records[rec_idx];
                let len = len.min(seq.len());
                let max_start = seq.len() - len;
                let start = if max_start == 0 { 0 } else { rng.next_range(max_start + 1) };
                seq[start..start + len].to_vec()
            }
            2 => {
                let rec_idx = rng.next_range(seq_count);
                let (_, seq) = &pair.records[rec_idx];
                let len = len.min(seq.len());
                let max_start = seq.len() - len;
                let start = if max_start == 0 { 0 } else { rng.next_range(max_start + 1) };
                let base = &seq[start..start + len];
                let k = rng.next_range(len / 3 + 1);
                mutate_k_positions(base, k, seed.wrapping_mul(97) + 13)
            }
            3 => {
                let rec_idx = rng.next_range(seq_count);
                let (_, seq) = &pair.records[rec_idx];
                let len = len.min(seq.len()).max(20);
                let max_start = seq.len() - len;
                let start = if max_start == 0 { 0 } else { rng.next_range(max_start + 1) };
                let base = &seq[start..start + len];
                let del_pos = rng.next_range(base.len().saturating_sub(10).max(1));
                let del_len = 1 + rng.next_range(5);
                delete_span(base, del_pos, del_len)
            }
            _ => {
                let rec_idx = rng.next_range(seq_count);
                let (_, seq) = &pair.records[rec_idx];
                let len = len.min(seq.len()).max(20);
                let max_start = seq.len() - len;
                let start = if max_start == 0 { 0 } else { rng.next_range(max_start + 1) };
                let base = &seq[start..start + len];
                let ins_pos = rng.next_range(base.len().saturating_sub(10).max(1));
                let ins_len = 1 + rng.next_range(5);
                insert_span(base, ins_pos, ins_len, seed.wrapping_mul(211) + 3)
            }
        };

        assert_overlaps_agree(
            &format!("random batch seed={seed} category={category} len={len}"),
            &pair,
            &read,
        );
        checked += 1;
    }
    assert!(checked >= 150, "expected >= 150 random reads to be checked, got {checked}");
}
