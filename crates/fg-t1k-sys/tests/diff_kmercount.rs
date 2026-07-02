#![cfg(feature = "t1k-sys")]
//! Differential test: drives a Rust `KmerCount` and a real C++ `KmerCount`
//! (via the opaque-handle FFI shim) with the same sequence of `AddCount`
//! calls, asserting the return values and every subsequent `GetCount` query
//! agree, plus a separate `GetCountSimilarityJaccard` comparison.

use fg_t1k_core::kmer_count::KmerCount;
use fg_t1k_sys::CppKmerCount;

/// A read set covering the scenarios called out in the task brief: reads
/// shorter than k, reads containing `N`, exact-duplicate reads (to exercise
/// count accumulation), and reverse-complement duplicates (to exercise
/// canonicalization).
fn read_set_a(k: usize) -> Vec<&'static [u8]> {
    let _ = k;
    vec![
        b"ACGTACGTACGTACGTACGTACGTACGTACGT", // long, all-ACGT
        b"TTTTAAAACCCCGGGGTTTTAAAACCCCGGGG",
        b"GATTACAGATTACAGATTACAGATTACAGATT",
        b"ACGTACGTACGTACGTACGTACGTACGTACGT", // exact duplicate of the first read
        b"ACGTNACGTACGTNACGTACGTACGTACGTAC", // contains N
        b"AC",                               // shorter than any tested k (2 bases)
        b"CCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA", // revcomp of the second read above
    ]
}

fn read_set_b(k: usize) -> Vec<&'static [u8]> {
    let _ = k;
    vec![
        b"CCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTT",
        b"AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
        b"ACGTACGTACGTACGTACGTACGTACGTACGT", // overlaps with read_set_a's first read
        b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        b"TG",
    ]
}

/// All length-`k` windows (as owned byte vectors) that occur in `reads`,
/// deduplicated, used to build a large set of "present" query k-mers.
fn all_windows(reads: &[&[u8]], k: usize) -> Vec<Vec<u8>> {
    let mut windows = Vec::new();
    for read in reads {
        if read.len() < k {
            continue;
        }
        for start in 0..=(read.len() - k) {
            windows.push(read[start..start + k].to_vec());
        }
    }
    windows.sort();
    windows.dedup();
    windows
}

/// A handful of k-mers guaranteed absent from either read set (all-A/all-T
/// homopolymers of the wrong k, plus alternating patterns unlikely to appear
/// in the constructed reads), plus explicit `N`-containing queries.
fn absent_and_n_queries(k: usize) -> Vec<Vec<u8>> {
    let mut queries = Vec::new();
    queries.push(vec![b'A'; k]); // homopolymer, unlikely to appear in test reads
    queries.push(vec![b'T'; k]);
    let mut alt = Vec::with_capacity(k);
    for i in 0..k {
        alt.push(if i % 2 == 0 { b'G' } else { b'C' });
    }
    queries.push(alt);

    // N-containing queries at a few different positions within the k-mer.
    let mut n_at_start = vec![b'A'; k];
    n_at_start[0] = b'N';
    queries.push(n_at_start);
    let mut n_at_end = vec![b'A'; k];
    n_at_end[k - 1] = b'N';
    queries.push(n_at_end);
    if k > 2 {
        let mut n_in_middle = vec![b'A'; k];
        n_in_middle[k / 2] = b'N';
        queries.push(n_in_middle);
    }

    queries
}

/// Feeds `reads` into both a Rust and a C++ `KmerCount` via `add_count`,
/// asserting the return values match at every step, then asserts `get_count`
/// agrees for every window present in `reads` plus a set of absent/`N`
/// queries. Returns the two populated counters so callers can additionally
/// exercise Jaccard similarity across two independently-built sets.
fn add_reads_and_check(k: usize, reads: &[&[u8]], label: &str) -> (KmerCount, CppKmerCount) {
    let mut rust = KmerCount::new(k);
    let mut cpp = CppKmerCount::new(i32::try_from(k).unwrap());

    for (i, read) in reads.iter().enumerate() {
        let rust_ret = rust.add_count(read);
        let cpp_ret = cpp.add_count(read);
        assert_eq!(
            rust_ret,
            cpp_ret,
            "{label} k={k}: add_count return mismatch on read #{i} ({:?})",
            String::from_utf8_lossy(read)
        );
    }

    for window in all_windows(reads, k) {
        let rust_count = rust.get_count(&window);
        let cpp_count = cpp.get_count(&window);
        assert_eq!(
            rust_count,
            cpp_count,
            "{label} k={k}: get_count mismatch for present window {:?}",
            String::from_utf8_lossy(&window)
        );
    }

    for query in absent_and_n_queries(k) {
        let rust_count = rust.get_count(&query);
        let cpp_count = cpp.get_count(&query);
        assert_eq!(
            rust_count,
            cpp_count,
            "{label} k={k}: get_count mismatch for absent/N query {:?}",
            String::from_utf8_lossy(&query)
        );
    }

    (rust, cpp)
}

fn run_add_count_and_get_count_diff(k: usize) {
    add_reads_and_check(k, &read_set_a(k), "read_set_a");
}

fn run_jaccard_diff(k: usize) {
    let (rust_a, cpp_a) = add_reads_and_check(k, &read_set_a(k), "jaccard/read_set_a");
    let (rust_b, cpp_b) = add_reads_and_check(k, &read_set_b(k), "jaccard/read_set_b");

    let rust_jaccard = rust_a.jaccard_similarity(&rust_b);
    let cpp_jaccard = cpp_a.jaccard(&cpp_b);
    assert_eq!(
        rust_jaccard.to_bits(),
        cpp_jaccard.to_bits(),
        "k={k}: jaccard similarity mismatch: rust={rust_jaccard} cpp={cpp_jaccard}"
    );
}

#[test]
fn add_count_and_get_count_match_t1k_k15() {
    run_add_count_and_get_count_diff(15);
}

#[test]
fn add_count_and_get_count_match_t1k_k31() {
    run_add_count_and_get_count_diff(31);
}

#[test]
fn jaccard_similarity_matches_t1k_k15() {
    run_jaccard_diff(15);
}

#[test]
fn jaccard_similarity_matches_t1k_k31() {
    run_jaccard_diff(31);
}

#[test]
fn kmercount_new_produces_a_usable_handle() {
    // Basic smoke test that the opaque handle itself is non-null-usable
    // (a NULL handle would panic inside CppKmerCount::new per its contract).
    let mut cpp = CppKmerCount::new(9);
    assert_eq!(cpp.add_count(b"ACGTACGTA"), 1);
    assert_eq!(cpp.get_count(b"ACGTACGTA"), 1);
}

#[test]
#[should_panic(expected = "get_count query must be at least k=9 bytes")]
fn get_count_rejects_query_shorter_than_k() {
    // The C++ GetCount reads exactly k bytes without calling strlen, so a
    // query slice shorter than k would read past the CString allocation.
    // The wrapper must reject it at the FFI boundary.
    let cpp = CppKmerCount::new(9);
    let _ = cpp.get_count(b"ACGT");
}
