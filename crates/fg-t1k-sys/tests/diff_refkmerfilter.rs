#![cfg(feature = "t1k-sys")]
//! Differential test: builds a Rust `RefKmerFilter` and a real C++ `SeqSet`
//! (via the opaque-handle FFI shim) from the SAME T1K-style reference FASTA,
//! then asserts `is_low_complexity`/`has_hit_in_set`/`is_good_candidate`
//! agree bool-for-bool for a varied set of reads.
//!
//! # Why `has_hit_in_set` can be compared against the REAL (full) C++ at all
//!
//! `fg_t1k_core::ref_kmer_filter::RefKmerFilter::has_hit_in_set` only
//! implements stock `SeqSet::HasHitInSet`'s bucket-count gate
//! (`SeqSet.hpp:1915-1964`), deliberately NOT the subsequent
//! `GetOverlapsFromHits`/`AlignAlgo`-based alignment confirmation
//! (`SeqSet.hpp:1966-1990`) -- see that module's doc comment for the full
//! rationale. This makes the Rust port a strict superset of stock's `true`
//! results: wherever stock returns `true`, so does this port, but this port
//! could in principle return `true` in rare cases stock's alignment step
//! would additionally reject (e.g. hits that are numerous but not
//! colinear/alignable).
//!
//! Every read below is deliberately chosen to avoid that gap: "hit" reads
//! are exact contiguous substrings of a loaded reference sequence (trivially
//! colinear and 100% similar, so stock's alignment step cannot reject them),
//! "non-hit" reads are either clearly unrelated to the reference or so short
//! they cannot possibly satisfy even the bucket-count gate. This is stated
//! up front so a future reader knows this test is not accidentally
//! vacuous -- it exercises the real C++ oracle's FULL `HasHitInSet`, not a
//! stub.

use fg_t1k_core::ref_kmer_filter::{RefKmerFilter, is_low_complexity};
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
/// a small xorshift64 so "random-looking" non-matching reads are
/// reproducible across runs without needing a new dependency.
fn pseudo_random_acgt(seed: u64, len: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut state = seed | 1; // xorshift64 requires a nonzero seed
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        out.push(bases[(state % 4) as usize]);
    }
    out
}

/// Asserts `is_low_complexity`, `has_hit_in_set`, and `is_good_candidate`
/// all agree bool-for-bool between `pair.rust` and `pair.cpp` for `read`,
/// tagging any failure with `label` for a readable assertion message.
///
/// `RefKmerFilter::has_hit_in_set` is `pub(crate)` inside `fg-t1k-core` (not
/// reachable from this crate), so it is not called directly here. Instead,
/// since `is_good_candidate == !is_low_complexity && has_hit_in_set` on BOTH
/// sides (`IsGoodCandidate`, `FastqExtractor.cpp:113-118`, and
/// `RefKmerFilter::is_good_candidate_with_scratch`'s identical formula),
/// checking `is_good_candidate` for a read where `is_low_complexity` is
/// `false` gives exactly the same coverage as checking `has_hit_in_set`
/// directly would.
fn assert_predicates_agree(label: &str, pair: &FilterPair, read: &[u8]) {
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

    if !rust_low_complexity {
        // is_good_candidate == has_hit_in_set here (is_low_complexity is
        // false on both sides, since we just asserted they agree), so this
        // is genuine has_hit_in_set-equivalent coverage.
        let cpp_hit = pair.cpp.has_hit_in_set(read);
        assert_eq!(
            rust_good, cpp_hit,
            "{label}: has_hit_in_set mismatch (rust(via is_good_candidate)={rust_good}, cpp={cpp_hit}, read={read:?})"
        );
    }
}

#[test]
fn kir_reference_exact_substrings_are_good_candidates() {
    let pair = load_pair("kir_rna_seq.fa");
    assert!(pair.rust.seq_count() >= 7);

    // Exact contiguous substrings of loaded reference sequences, at varied
    // lengths/offsets/records -- trivially colinear and 100%-similar, so
    // stock's excluded alignment-confirmation step cannot reject them (see
    // this file's module doc comment).
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
    // sides, so this should still comfortably pass the bucket-count gate
    // AND the (excluded, but real-on-the-C++-side) alignment confirmation.
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
