#![cfg(feature = "t1k-sys")]
//! Differential test: builds a Rust `RefKmerFilter` and a real C++ `SeqSet`
//! (via the opaque-handle FFI shim) from the SAME T1K-style reference FASTA,
//! then asserts `is_low_complexity`/`has_hit_in_set`/`is_good_candidate`
//! agree bool-for-bool for a large, deliberately ADVERSARIAL set of reads --
//! not just curated exact-substring reads.
//!
//! # `has_hit_in_set` is now byte/value-identical to stock, not a superset
//!
//! `fg_t1k_core::ref_kmer_filter::RefKmerFilter::has_hit_in_set` now
//! implements BOTH of stock `SeqSet::HasHitInSet`'s gates: the bucket-count
//! gate (`SeqSet.hpp:1915-1964`, ported in Task 3.1) AND the
//! `GetOverlapsFromHits`-based (LIS hit-chaining) mismatch-threshold
//! confirmation (`SeqSet.hpp:1966-1990`, ported in this task -- see
//! `fg_t1k_core::overlap`'s module docs for why no `AlignAlgo`/Smith-Waterman
//! is needed on this path). This test's job is to prove that completion:
//! every read category below is chosen to actually exercise gate 2 (not just
//! gate 1), including reads that gate 1 alone would accept but gate 2
//! (correctly) rejects -- see `gate_2_rejects_reads_gate_1_alone_would_accept`
//! for the test that counts and asserts this gap is non-empty.

use fg_t1k_core::ref_kmer_filter::{RefKmerFilter, Scratch, is_low_complexity};
use fg_t1k_sys::CppSeqSet;
use std::path::{Path, PathBuf};

/// The k-mer length stock T1K would actually select for both fixture
/// references used here. `SeqSet::InferKmerLength()` (`SeqSet.hpp:2830-2845`)
/// evaluates to `8` for `kir_rna_seq.fa` (total reference length 8781) and
/// also `8` for `hla_rna_seq.fa` (total reference length 12260) -- neither
/// exceeds `FastqExtractor.cpp`'s initial literal default of `9`
/// (`FastqExtractor.cpp:272`), so stock never rebuilds via
/// `UpdateKmerLength` for either fixture, and ends up using `9` unmodified.
/// See `ref_kmer_filter`'s module doc comment for the full explanation of
/// why `kmer_length` is an explicit parameter here rather than replicated
/// two-stage inference logic.
const KMER_LENGTH: usize = 9;

fn fixture_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/refbuild/golden").join(rel)
}

/// A loaded (Rust, C++) pair built from the same reference FASTA at
/// [`KMER_LENGTH`], plus the raw per-record sequences (for pulling exact
/// substrings as "hit" test reads).
struct FilterPair {
    rust: RefKmerFilter,
    cpp: CppSeqSet,
    /// Concatenated (id, sequence) pairs in file order, parsed independently
    /// of both `RefKmerFilter`/`CppSeqSet` (a third, trivial parse) so that
    /// test-read construction does not depend on either port's internal
    /// FASTA handling.
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

/// Trivial, independent FASTA parser used only by this test to pull exact
/// reference substrings for "hit" reads -- deliberately NOT shared code with
/// `fg_t1k_core::ref_kmer_filter`'s own (private) FASTA parser, so a bug in
/// that parser cannot silently make this test agree with itself.
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

/// Deterministic (seeded, reproducible) pseudo-random ACGT sequence
/// generator -- NOT `rand`-crate-backed (no such dev-dependency here), just
/// a small xorshift64 so "random-looking" reads/mutations are reproducible
/// across runs without needing a new dependency.
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        Self { state: seed | 1 } // xorshift64 requires a nonzero seed
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

    /// A base guaranteed to differ from `exclude`.
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

/// Introduces exactly `k` point substitutions into `seq` at deterministic
/// (seeded), distinct positions, each guaranteed to actually change the
/// base (never a same-base "mutation").
fn mutate_k_positions(seq: &[u8], k: usize, seed: u64) -> Vec<u8> {
    let mut rng = Xorshift64::new(seed);
    let mut out = seq.to_vec();
    let mut positions: Vec<usize> = (0..seq.len()).collect();
    // Fisher-Yates partial shuffle to pick k distinct positions
    // deterministically from the seeded rng.
    for i in 0..k.min(positions.len()) {
        let j = i + rng.next_range(positions.len() - i);
        positions.swap(i, j);
    }
    for &pos in positions.iter().take(k.min(positions.len())) {
        out[pos] = rng.next_base_excluding(out[pos]);
    }
    out
}

/// The outcome of a single `assert_predicates_agree` check, used by the
/// gate-2-exercising census in `gate_2_rejects_reads_gate_1_alone_would_accept`.
struct AgreementResult {
    /// Whether Rust's gate-1-only decision (`passes_bucket_count_gate_only`)
    /// was `true` for this read.
    gate1_only_rust: bool,
    /// Whether the REAL C++ `HasHitInSet` (both gates) was `true`.
    cpp_hit: bool,
}

/// Asserts `is_low_complexity`, `has_hit_in_set` (via `is_good_candidate`
/// under `!is_low_complexity`), and `is_good_candidate` all agree
/// bool-for-bool between `pair.rust` and `pair.cpp` for `read`, tagging any
/// failure with `label` for a readable assertion message. Returns the
/// gate-1-only-vs-full-C++ comparison for the caller's own bookkeeping (used
/// by the gate-2-exercising census test).
fn assert_predicates_agree(label: &str, pair: &FilterPair, read: &[u8]) -> AgreementResult {
    let rust_low_complexity = is_low_complexity(read);
    let cpp_low_complexity = pair.cpp.is_low_complexity(read);
    assert_eq!(
        rust_low_complexity, cpp_low_complexity,
        "{label}: is_low_complexity mismatch (rust={rust_low_complexity}, cpp={cpp_low_complexity}, read={read:?})"
    );

    let rust_good = pair.rust.is_good_candidate(read);
    let cpp_good = pair.cpp.is_good_candidate(read);
    assert_eq!(
        rust_good, cpp_good,
        "{label}: is_good_candidate mismatch (rust={rust_good}, cpp={cpp_good}, read={read:?})"
    );

    let mut cpp_hit = cpp_good;
    if !rust_low_complexity {
        // is_good_candidate == has_hit_in_set here (is_low_complexity is
        // false on both sides, since we just asserted they agree), so this
        // is genuine has_hit_in_set-equivalent coverage.
        cpp_hit = pair.cpp.has_hit_in_set(read);
        assert_eq!(
            rust_good, cpp_hit,
            "{label}: has_hit_in_set mismatch (rust(via is_good_candidate)={rust_good}, cpp={cpp_hit}, read={read:?})"
        );
    }

    let mut scratch = Scratch::default();
    let gate1_only_rust = pair.rust.passes_bucket_count_gate_only(read, &mut scratch);

    AgreementResult { gate1_only_rust, cpp_hit }
}

#[test]
fn kir_reference_exact_substrings_are_good_candidates() {
    let pair = load_pair("kir_rna_seq.fa");
    assert!(pair.rust.seq_count() >= 7);

    // Exact contiguous substrings of loaded reference sequences, at varied
    // lengths/offsets/records -- trivially colinear and 100%-similar, so
    // stock's alignment-confirmation step cannot reject them.
    let cases: &[(usize, usize, usize)] = &[
        // (record index, start offset, length)
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
        assert_predicates_agree(&format!("kir substr {id}@{start}+{len}"), &pair, read);

        // Reverse-complement of the same substring: GetHitsFromRead always
        // searches both strands, so this must also agree.
        let rc = reverse_complement(read);
        assert_predicates_agree(&format!("kir substr RC {id}@{start}+{len}"), &pair, &rc);
    }
}

#[test]
fn hla_reference_exact_substrings_are_good_candidates() {
    let pair = load_pair("hla_rna_seq.fa");
    // `cases` below indexes records up to `pair.records[11]`, so require at
    // least 12 records; otherwise a smaller fixture would panic out-of-bounds
    // at the indexing step before the per-case bounds check below.
    assert!(pair.rust.seq_count() >= 12);

    let cases: &[(usize, usize, usize)] =
        &[(0, 0, 70), (1, 300, 90), (4, 700, 110), (8, 150, 60), (11, 700, 130)];
    for &(rec_idx, start, len) in cases {
        let (id, seq) = &pair.records[rec_idx];
        assert!(start + len <= seq.len(), "case out of bounds for {id}");
        let read = &seq[start..start + len];
        assert_predicates_agree(&format!("hla substr {id}@{start}+{len}"), &pair, read);
    }
}

/// Largest read the shim's fixed `char rcBuf[100001]` reverse-complement buffer
/// can hold: 100_000 read bytes plus the trailing NUL. Mirrors the private
/// `CppSeqSet::SHIM_RC_BUF_LEN` (`= 100_001`); a read this long is the boundary
/// that must be accepted, one byte more must be rejected.
const MAX_SHIM_READ_LEN: usize = 100_000;

// The shim reverse-complements each read into a fixed 100_001-byte stack buffer,
// so both FFI predicates must reject an oversized read (via the safe-wrapper
// length guard) rather than overflow the buffer inside C++. These pin that guard.

#[test]
#[should_panic(expected = "does not fit CppSeqSet::has_hit_in_set")]
fn has_hit_in_set_rejects_oversized_read() {
    let pair = load_pair("kir_rna_seq.fa");
    let too_long = vec![b'A'; MAX_SHIM_READ_LEN + 1];
    let _ = pair.cpp.has_hit_in_set(&too_long);
}

#[test]
#[should_panic(expected = "does not fit CppSeqSet::is_good_candidate")]
fn is_good_candidate_rejects_oversized_read() {
    let pair = load_pair("kir_rna_seq.fa");
    let too_long = vec![b'A'; MAX_SHIM_READ_LEN + 1];
    let _ = pair.cpp.is_good_candidate(&too_long);
}

/// A read of exactly the maximum length must pass the guard and run through the
/// shim's reverse-complement without panicking or overflowing -- proving the
/// boundary is `>= SHIM_RC_BUF_LEN`, not an over-strict off-by-one.
#[test]
fn maximum_length_read_is_accepted_without_overflow() {
    let pair = load_pair("kir_rna_seq.fa");
    let max_read = vec![b'A'; MAX_SHIM_READ_LEN];
    // `has_hit_in_set` always reverse-complements into `rcBuf`, exercising the
    // buffer at its exact capacity; the return value is irrelevant here.
    let _ = pair.cpp.has_hit_in_set(&max_read);
    let _ = pair.cpp.is_good_candidate(&max_read);
}

#[test]
fn substring_with_a_single_n_still_agrees() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[0];
    let start = 100;
    let len = 100;
    let mut read = seq[start..start + len].to_vec();
    // Corrupt a single base in the middle to 'N' -- reduces valid k-mer
    // windows locally but leaves plenty of exact matching windows on both
    // sides, so this should still comfortably pass both gates.
    let mid = len / 2;
    read[mid] = b'N';
    assert_predicates_agree(&format!("kir substr-with-N {id}@{start}+{len}"), &pair, &read);
}

#[test]
fn random_non_matching_reads_are_not_candidates() {
    let pair = load_pair("kir_rna_seq.fa");
    // Several different seeds/lengths of deterministic pseudo-random ACGT
    // sequence, verified (by this test passing) to not accidentally collect
    // enough reference hits to pass the bucket-count gate on either side.
    for (seed, len) in [(1u64, 80), (2, 120), (3, 60), (4, 150), (5, 40)] {
        let read = pseudo_random_acgt(seed, len);
        assert_predicates_agree(&format!("random seed={seed} len={len}"), &pair, &read);
    }
}

#[test]
fn low_complexity_reads_agree() {
    let pair = load_pair("kir_rna_seq.fa");

    let homopolymer_a = vec![b'A'; 100];
    assert_predicates_agree("homopolymer A", &pair, &homopolymer_a);

    let homopolymer_t = vec![b'T'; 60];
    assert_predicates_agree("homopolymer T", &pair, &homopolymer_t);

    let dinucleotide: Vec<u8> = (0..100).map(|i| if i % 2 == 0 { b'A' } else { b'T' }).collect();
    assert_predicates_agree("dinucleotide AT repeat", &pair, &dinucleotide);

    let dinucleotide_gc: Vec<u8> = (0..80).map(|i| if i % 2 == 0 { b'G' } else { b'C' }).collect();
    assert_predicates_agree("dinucleotide GC repeat", &pair, &dinucleotide_gc);

    // Mostly-N read: >= len/10 Ns triggers IsLowComplexity's N-fraction
    // check regardless of the rest of the composition.
    let mut mostly_n = pseudo_random_acgt(42, 100);
    for b in mostly_n.iter_mut().take(20) {
        *b = b'N';
    }
    assert_predicates_agree("mostly-N", &pair, &mostly_n);
}

#[test]
fn short_reads_agree() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[0];

    // Shorter than kmer_length (9): both sides gate on `len < kmerLength`
    // before ever computing a k-mer window.
    let too_short = &seq[0..5];
    assert_predicates_agree(&format!("{id} len=5 (< kmer_length)"), &pair, too_short);

    // Exactly kmer_length: exactly one possible k-mer window per strand.
    // A single window's hit count is very unlikely to reach the
    // hitLenRequired=31 gate (needs >= ceil(31/9)=4 hits in one bucket) for
    // an arbitrary reference position, so this is expected to NOT be a
    // candidate on either side -- itself a useful boundary case.
    let exactly_k = &seq[50..50 + KMER_LENGTH];
    assert_predicates_agree(&format!("{id} len=kmer_length exactly"), &pair, exactly_k);

    // One base longer than kmer_length: two overlapping windows.
    let just_over = &seq[50..=50 + KMER_LENGTH];
    assert_predicates_agree(&format!("{id} len=kmer_length+1"), &pair, just_over);
}

#[test]
fn high_frequency_kmer_skip_agrees_with_stock() {
    // Regression test for the SeqSet.hpp:1109/1176 `size >= 100`
    // high-frequency-kmer skip branch inside `GetHitsFromRead`, exercised
    // against the REAL C++ oracle (not just fg-t1k-core's own pure-Rust
    // unit tests for this exact branch -- see ref_kmer_filter.rs's
    // `high_frequency_kmer_skipped_when_not_at_read_boundary` and
    // `high_frequency_kmer_not_skipped_at_first_window`, which independently
    // hand-verify raw hit counts that are not observable through this
    // crate's is_good_candidate-only boundary). Neither KIR/HLA fixture is
    // large/repetitive enough for any k-mer to reach 100 index entries, so
    // this uses a synthetic single-sequence reference: `"ACGT"` repeated 150
    // times (600bp). At kmer_length=4, each of the 4 rotations of "ACGT"
    // ("ACGT"/"CGTA"/"GTAC"/"TACG") independently reaches ~150 index
    // entries -- all within this ONE sequence, so they concentrate into a
    // single (tag, seqIdx) bucket rather than spreading across many
    // sequences (spreading across many sequences, e.g. many distinct
    // single-copy records, would never pass the bucket-count gate at all,
    // since each bucket would only ever hold 1 hit -- a design pitfall this
    // comment flags explicitly since it is easy to get wrong). Critically,
    // this reference's periodic-but-4-base-balanced composition keeps
    // `IsLowComplexity` FALSE for both the reference and any exact
    // substring read of it (unlike a homopolymer, which would trip the
    // low-complexity gate before `HasHitInSet` is ever reached).
    use std::io::Write as _;
    let mut tmp = tempfile::NamedTempFile::new().unwrap();
    writeln!(tmp, ">only").unwrap();
    writeln!(tmp, "{}", "ACGT".repeat(150)).unwrap();
    tmp.flush().unwrap();

    let kmer_length = 4usize;
    let rust = RefKmerFilter::from_reference_fasta(tmp.path(), kmer_length).unwrap();
    let mut cpp = CppSeqSet::new(i32::try_from(kmer_length).unwrap());
    cpp.load_ref(tmp.path());

    // A 70bp exact substring of the repetitive reference: every window is
    // one of the four high-frequency rotations, so most (though not all --
    // the skip-limit logic still fires on non-boundary windows) get
    // skipped, yet the sheer repeat count still leaves thousands of
    // unskipped hits concentrated in the reference's single bucket --
    // comfortably over `hitLenRequired`.
    let reference = "ACGT".repeat(150);
    let hit_read = &reference.as_bytes()[10..80];
    assert!(!is_low_complexity(hit_read), "test read must not be low-complexity");
    let rust_hit = rust.is_good_candidate(hit_read);
    let cpp_hit = cpp.is_good_candidate(hit_read);
    assert_eq!(rust_hit, cpp_hit, "high-frequency-but-genuine-match mismatch");
    assert!(rust_hit, "a genuine repeat-region match must still be a candidate despite the skip");

    // An unrelated (non-matching, non-low-complexity) read: "TGCA" repeated
    // is a different permutation from "ACGT" (shares no rotation with it,
    // so no shared k-mers at all) but is equally base-balanced, so it stays
    // non-low-complexity too.
    let unrelated = "TGCA".repeat(20).into_bytes();
    assert!(!is_low_complexity(&unrelated));
    let rust_unrelated = rust.is_good_candidate(&unrelated);
    let cpp_unrelated = cpp.is_good_candidate(&unrelated);
    assert_eq!(rust_unrelated, cpp_unrelated, "unrelated-read mismatch");
    assert!(!rust_unrelated);
}

/// Reference substrings with a swept number `k` of point substitutions
/// (0..K), at varied positions -- directly exercises gate 2's mismatch
/// threshold (`SeqSet.hpp:1973`: `mismatchThreshold = int(len *
/// (1-refSeqSimilarity)) * kmerLength`). At `refSeqSimilarity=0.8`,
/// `kmerLength=9`, a 150bp read has `mismatchThreshold = int(150*0.2)*9 =
/// 30*9 = 270` -- but that threshold is compared against `len -
/// matchCnt/2`, and `matchCnt` comes from the LIS-chained hit length (a much
/// smaller quantity than `len` once mismatches fragment the k-mer hit
/// chain), so even a handful of mismatches can flip the gate. This test
/// sweeps k from 0 up to a value chosen per-length to comfortably cross that
/// flip point on both sides, asserting bool-for-bool agreement at every
/// step (not just at the extremes).
#[test]
fn mismatch_sweep_exercises_gate_2_threshold() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[2];
    let start = 100;
    let len = 150;
    let base_read = &seq[start..start + len];

    let mut saw_true = false;
    let mut saw_false = false;
    for k in 0..=40usize {
        let read = mutate_k_positions(base_read, k, 1000 + k as u64);
        let result =
            assert_predicates_agree(&format!("{id}@{start}+{len} mutated k={k}"), &pair, &read);
        if result.cpp_hit {
            saw_true = true;
        } else {
            saw_false = true;
        }
    }
    assert!(saw_true, "mismatch sweep never observed a C++ HasHitInSet==true case");
    assert!(
        saw_false,
        "mismatch sweep never observed a C++ HasHitInSet==false case (sweep didn't cross the threshold -- gate 2 not exercised)"
    );
}

/// A second mismatch sweep on a different record/fixture, at a different
/// base length, so the threshold-crossing behavior above is not an artifact
/// of one specific record.
#[test]
fn mismatch_sweep_hla_exercises_gate_2_threshold() {
    let pair = load_pair("hla_rna_seq.fa");
    let (id, seq) = &pair.records[4];
    let start = 200;
    let len = 100;
    let base_read = &seq[start..start + len];

    let mut saw_true = false;
    let mut saw_false = false;
    for k in 0..=30usize {
        let read = mutate_k_positions(base_read, k, 2000 + k as u64);
        let result =
            assert_predicates_agree(&format!("{id}@{start}+{len} mutated k={k}"), &pair, &read);
        if result.cpp_hit {
            saw_true = true;
        } else {
            saw_false = true;
        }
    }
    assert!(saw_true);
    assert!(saw_false, "hla mismatch sweep never crossed the gate-2 threshold");
}

/// Two disjoint reference fragments (from different, distant parts of the
/// SAME reference sequence, or different sequences entirely) concatenated
/// into one read: noncolinear on the reference, so the LIS-based chaining
/// in `GetOverlapsFromHits` can only pick ONE of the two fragments as its
/// chain (LIS is a chain in read AND seq coordinates simultaneously), never
/// both -- meaning gate 2's mismatch threshold sees roughly half the read
/// as "unmatched", which can flip the read from a bucket-gate pass to an
/// overall `HasHitInSet` failure. This directly targets the LIS
/// noncolinear-chain-selection logic (SeqSet.hpp:1360-1471), not just the
/// mismatch-count arithmetic the point-mutation sweep above targets.
#[test]
fn noncolinear_concatenated_fragments_exercise_lis_chaining() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id0, seq0) = &pair.records[0];
    let (id1, seq1) = &pair.records[1];

    // Two 75bp fragments from unrelated positions (and, for the "different
    // sequences" case, unrelated sequences): concatenating them creates a
    // read whose first half maps to one (seq, coordinate) region and second
    // half to a totally different one -- noncolinear on any single
    // reference sequence.
    let frag_a = &seq0[50..125];
    let frag_b = &seq0[900..975]; // same sequence, far-apart offset
    let mut same_seq_concat = frag_a.to_vec();
    same_seq_concat.extend_from_slice(frag_b);
    assert_predicates_agree(
        &format!("{id0} noncolinear same-seq concat (50..125 + 900..975)"),
        &pair,
        &same_seq_concat,
    );

    let frag_c = &seq1[400..475];
    let mut cross_seq_concat = frag_a.to_vec();
    cross_seq_concat.extend_from_slice(frag_c);
    assert_predicates_agree(
        &format!("{id0}/{id1} noncolinear cross-seq concat"),
        &pair,
        &cross_seq_concat,
    );

    // Near-adjacent fragments (small gap, still on the reference but not
    // perfectly contiguous -- e.g. simulating a small deletion): should
    // still chain successfully via LIS/adjustRadius handling in most cases,
    // exercising the "should still pass" side of the same machinery.
    let mut near_adjacent = seq0[200..280].to_vec();
    near_adjacent.extend_from_slice(&seq0[285..365]); // 5bp gap
    assert_predicates_agree(&format!("{id0} near-adjacent (5bp gap)"), &pair, &near_adjacent);
}

/// Reads at/near the `len`/`hitLenRequired` boundary: since gate 2's
/// mismatch threshold and gate 1's bucket-count gate both scale with `len`
/// and `kmerLength`, reads right at the edge of passing/failing are the
/// most likely to reveal an off-by-one in either port.
#[test]
fn boundary_length_reads_agree() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[3];

    // hitLenRequired=31 at kmer_length=9 needs >= ceil(31/9)=4 hits in the
    // winning bucket; a read of exactly kmer_length*4 - 1 = 35 has just
    // barely enough length for 4 overlapping windows (offsets 0..=26 give
    // 27 possible starting positions, way more than 4, so this is really
    // testing the gate-2 mismatch-threshold boundary, not gate 1's).
    for len in [35usize, 40, 44, 45, 46, 50, 60] {
        let start = 300;
        let read = &seq[start..start + len];
        assert_predicates_agree(&format!("{id}@{start}+{len} boundary length"), &pair, read);
    }
}

/// At least 200 seeded pseudo-random reads of varied lengths (40-200),
/// covering pure noise (unrelated to any reference), reference substrings,
/// and reference substrings with a random number of point mutations -- a
/// broad, non-curated adversarial sweep as required by this task.
#[test]
fn large_random_read_batch_agrees() {
    let pair = load_pair("hla_rna_seq.fa");
    let seq_count = pair.records.len();

    let mut agreements = 0usize;
    for seed in 0u64..220 {
        let mut rng = Xorshift64::new(seed ^ 0xABCD_1234);
        let len = 40 + rng.next_range(161); // 40..=200

        let category = rng.next_range(3);
        let read: Vec<u8> = match category {
            0 => {
                // Pure noise, unrelated to the reference.
                pseudo_random_acgt(seed.wrapping_mul(31) + 7, len)
            }
            1 => {
                // Exact reference substring.
                let rec_idx = rng.next_range(seq_count);
                let (_, seq) = &pair.records[rec_idx];
                let len = len.min(seq.len());
                let max_start = seq.len() - len;
                let start = if max_start == 0 { 0 } else { rng.next_range(max_start + 1) };
                seq[start..start + len].to_vec()
            }
            _ => {
                // Reference substring with a random number of mutations.
                let rec_idx = rng.next_range(seq_count);
                let (_, seq) = &pair.records[rec_idx];
                let len = len.min(seq.len());
                let max_start = seq.len() - len;
                let start = if max_start == 0 { 0 } else { rng.next_range(max_start + 1) };
                let base = &seq[start..start + len];
                let k = rng.next_range(len / 3 + 1);
                mutate_k_positions(base, k, seed.wrapping_mul(97) + 13)
            }
        };

        assert_predicates_agree(
            &format!("random batch seed={seed} category={category} len={len}"),
            &pair,
            &read,
        );
        agreements += 1;
    }
    assert!(agreements >= 200, "expected >= 200 random reads to be checked, got {agreements}");
}

/// Reverse-complement variants of a representative slice of the mutated /
/// noncolinear / random-batch categories above, plus reads containing `N`
/// and very short (but >= kmer_length) reads -- rounding out strand and
/// edge-case coverage specifically requested by the task brief.
#[test]
fn rc_variants_n_containing_and_short_reads_agree() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[4];

    // RC of a mutated substring.
    let base = &seq[80..200];
    let mutated = mutate_k_positions(base, 8, 555);
    let rc_mutated = reverse_complement(&mutated);
    assert_predicates_agree(&format!("{id} RC mutated k=8"), &pair, &rc_mutated);

    // RC of a noncolinear concatenation.
    let frag_a = seq[10..90].to_vec();
    let frag_b = seq[600..680].to_vec();
    let mut concat = frag_a;
    concat.extend_from_slice(&frag_b);
    let rc_concat = reverse_complement(&concat);
    assert_predicates_agree(&format!("{id} RC noncolinear concat"), &pair, &rc_concat);

    // Reads with N at varied positions/counts, still short of the
    // low-complexity N-fraction trigger (< len/10).
    let mut with_one_n = seq[300..400].to_vec();
    with_one_n[50] = b'N';
    assert_predicates_agree(&format!("{id} single N"), &pair, &with_one_n);

    let mut with_few_n = seq[300..400].to_vec();
    for &pos in &[10, 40, 70] {
        with_few_n[pos] = b'N';
    }
    assert_predicates_agree(&format!("{id} three Ns"), &pair, &with_few_n);

    // Short reads (>= kmer_length, but small): exact substring and mutated.
    let short_exact = &seq[500..500 + KMER_LENGTH + 5];
    assert_predicates_agree(&format!("{id} short exact"), &pair, short_exact);

    let short_mutated = mutate_k_positions(short_exact, 2, 777);
    assert_predicates_agree(&format!("{id} short mutated"), &pair, &short_mutated);
}

/// Reads engineered to create `std::sort` comparator ties in
/// `CompSortHitCoordDiff`/`CompSortPairBInc` -- specifically, reads with
/// REPEATED identical `(readOffset - seqOffset)` coordinate diffs, which
/// happens naturally whenever a read contains a tandem repeat / duplicated
/// motif that also occurs in the reference (each repeat unit produces hits
/// with the same coordinate diff as its neighbors). See `overlap.rs`'s
/// `comp_sort_hit_coord_diff` doc comment for why ties in that specific
/// comparator are provably unobservable (bit-identical tied entries), but
/// this test exercises the REAL C++ oracle on such inputs anyway, as the
/// task brief requires, rather than relying solely on that argument.
#[test]
fn tie_inducing_tandem_repeat_reads_agree() {
    let pair = load_pair("kir_rna_seq.fa");
    let (id, seq) = &pair.records[5];

    // A short reference motif, repeated several times back-to-back: each
    // repeat's k-mer hits (relative to ANY single reference occurrence of
    // the motif) share the same read-vs-seq coordinate diff pattern as its
    // neighbors, which is exactly the kind of repetitive structure that
    // produces duplicate/tied hitCoordDiff entries.
    let motif = &seq[40..70]; // 30bp motif, well over kmer_length=9
    let mut repeated = Vec::new();
    for _ in 0..4 {
        repeated.extend_from_slice(motif);
    }
    assert_predicates_agree(&format!("{id} tandem repeat x4 (30bp motif)"), &pair, &repeated);

    // Same idea, but the repeat unit itself is an exact reference substring
    // embedded amid random padding on both sides -- less likely to be
    // low-complexity (varied composition) while still causing internal
    // coordinate-diff duplication from the repeated core.
    let mut padded_repeat = pseudo_random_acgt(9001, 20);
    for _ in 0..3 {
        padded_repeat.extend_from_slice(motif);
    }
    padded_repeat.extend_from_slice(&pseudo_random_acgt(9002, 20));
    assert_predicates_agree(&format!("{id} padded tandem repeat x3"), &pair, &padded_repeat);

    // A longer motif repeated twice with a small number of mutations in the
    // second copy -- still produces (near-)duplicate coordinate diffs for
    // the unmutated windows, while also exercising the mismatch-threshold
    // machinery simultaneously.
    let motif2 = &seq[500..580]; // 80bp
    let mutated_copy = mutate_k_positions(motif2, 3, 424_242);
    let mut two_copies = motif2.to_vec();
    two_copies.extend_from_slice(&mutated_copy);
    assert_predicates_agree(&format!("{id} two motif copies, second mutated"), &pair, &two_copies);
}

/// **The core proof this task requires**: counts how many reads across a
/// dedicated adversarial batch have C++ `HasHitInSet == false` while the
/// Task-3.1 gate-1-only decision (`passes_bucket_count_gate_only`) would
/// have been `true` -- i.e. reads that would have INCORRECTLY passed under
/// the old (gate-1-only) port, and are now correctly rejected by gate 2.
/// This count MUST be `> 0`; if it were `0`, this whole differential
/// wouldn't actually be exercising gate 2's logic at all.
///
/// Deliberately reuses/re-derives reads from the categories above known to
/// approach or cross the gate-2 threshold (mismatch sweeps and noncolinear
/// concatenations), since those are exactly the categories expected to
/// produce a gate-1-true/gate-2-false split.
#[test]
fn gate_2_rejects_reads_gate_1_alone_would_accept() {
    let pair = load_pair("kir_rna_seq.fa");
    let hla_pair = load_pair("hla_rna_seq.fa");

    let mut gate1_true_cpp_false = 0usize;
    let mut total = 0usize;

    // Mismatch sweeps across several records/lengths.
    let mismatch_cases: &[(&FilterPair, usize, usize, usize, u64)] = &[
        (&pair, 2, 100, 150, 5000),
        (&pair, 4, 50, 120, 6000),
        (&pair, 6, 10, 90, 7000),
        (&hla_pair, 4, 200, 100, 8000),
        (&hla_pair, 8, 100, 130, 9000),
    ];
    for &(p, rec_idx, start, len, seed_base) in mismatch_cases {
        let (id, seq) = &p.records[rec_idx];
        let base_read = &seq[start..start + len];
        for k in 0..=(len / 3) {
            let read = mutate_k_positions(base_read, k, seed_base + k as u64);
            let result = assert_predicates_agree(
                &format!("{id}@{start}+{len} k={k} (gate2 census)"),
                p,
                &read,
            );
            total += 1;
            if result.gate1_only_rust && !result.cpp_hit {
                gate1_true_cpp_false += 1;
            }
        }
    }

    // Noncolinear concatenations across several fragment pairs/records.
    let (id0, seq0) = &pair.records[0];
    let (id1, seq1) = &pair.records[1];
    let (id6, seq6) = &pair.records[6];
    let noncolinear_cases: &[(&str, Vec<u8>)] = &[
        ("kir0 50..125+900..975", {
            let mut v = seq0[50..125].to_vec();
            v.extend_from_slice(&seq0[900..975]);
            v
        }),
        ("kir0/kir1 cross-seq", {
            let mut v = seq0[100..200].to_vec();
            v.extend_from_slice(&seq1[300..400]);
            v
        }),
        ("kir6 far-apart halves", {
            let mut v = seq6[0..60].to_vec();
            v.extend_from_slice(&seq6[300..360]);
            v
        }),
        ("kir0/kir6 cross-seq large", {
            let mut v = seq0[0..150].to_vec();
            v.extend_from_slice(&seq6[100..250]);
            v
        }),
    ];
    for (label, read) in noncolinear_cases {
        let result = assert_predicates_agree(&format!("{id0}/{id1}/{id6} {label}"), &pair, read);
        total += 1;
        if result.gate1_only_rust && !result.cpp_hit {
            gate1_true_cpp_false += 1;
        }
    }

    eprintln!(
        "gate_2_rejects_reads_gate_1_alone_would_accept: {gate1_true_cpp_false} / {total} reads \
         exercised gate 2 (gate-1-only would have been true, real C++ HasHitInSet is false)"
    );
    assert!(
        gate1_true_cpp_false > 0,
        "expected at least one read where gate 1 alone would accept but C++ HasHitInSet \
         rejects (i.e. gate 2 must be exercised) -- got 0 out of {total} reads; the adversarial \
         batch failed to reach the gate-2 mismatch/noncolinearity threshold"
    );
}
