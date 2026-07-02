#![cfg(feature = "t1k-sys")]
//! Differential test: drives a Rust `KmerIndex` and a real C++ `KmerIndex`
//! (via the opaque-handle FFI shim) with the same sequence of
//! `Insert`/`Remove`/`Search`/`BuildIndexFromRead` calls, asserting every
//! `Search` result matches element-for-element AND in order.
//!
//! Faithfulness traps this test specifically targets (see the module docs on
//! `fg_t1k_core::kmer_index` for the underlying C++ semantics):
//! - `Insert`/`Search` key on the *forward* `KmerCode::GetCode()`, not the
//!   canonical (min of forward/revcomp) code -- a reverse-complement pair of
//!   k-mers MUST land in different buckets.
//! - `strand` is accepted by `Insert`/`Remove` but never stored/consulted.
//! - `IsValid()` gates both `Insert` (drops silently) and `Search` (returns
//!   empty).
//! - `Search` results preserve insertion order (`SimpleVector::PushBack`
//!   append-only; `Remove` shifts left, doesn't swap-remove).
//! - `BuildIndexFromRead`'s consecutive-duplicate-k-mer dedup (`i == kl ||
//!   !kmerCode.IsEqual(prevKmerCode)`), including the edge case where the
//!   very first full k-mer's code coincidentally equals the freshly
//!   constructed `prevKmerCode`'s all-zero initial code (an all-`A` k-mer),
//!   which causes that first k-mer to be dropped.

use fg_t1k_core::kmer::KmerCode;
use fg_t1k_core::kmer_index::{IndexInfo, KmerIndex};
use fg_t1k_sys::{CppIndexInfo, CppKmerCode, CppKmerIndex};

/// Builds a matching pair of (Rust, C++) `KmerCode`s for k-mer length `k`,
/// fed the same `bytes` via `Append` (mirrors how `KmerCount::get_count` and
/// `KmerIndex::BuildIndexFromRead` build k-mer codes -- append, not prepend).
fn encode_pair(k: usize, bytes: &[u8]) -> (KmerCode, CppKmerCode) {
    let mut rust = KmerCode::new(k);
    let mut cpp = CppKmerCode::new(i32::try_from(k).unwrap());
    for &c in bytes {
        rust.append(c);
        cpp.append(c);
    }
    (rust, cpp)
}

/// Asserts a Rust `search` result (`&[IndexInfo]`) and a C++ `search` result
/// (`Vec<CppIndexInfo>`) contain the same `(idx, offset)` pairs IN ORDER
/// (order-sensitive: `SimpleVector`/`Vec` are append-only lists, not sets).
fn assert_search_matches(label: &str, rust_list: &[IndexInfo], cpp_list: &[CppIndexInfo]) {
    assert_eq!(
        rust_list.len(),
        cpp_list.len(),
        "{label}: search size mismatch (rust={rust_list:?}, cpp={cpp_list:?})"
    );
    for (i, (r, c)) in rust_list.iter().zip(cpp_list.iter()).enumerate() {
        assert_eq!(r.idx, c.idx, "{label}: entry #{i} idx mismatch (rust={r:?}, cpp={c:?})");
        assert_eq!(
            r.offset, c.offset,
            "{label}: entry #{i} offset mismatch (rust={r:?}, cpp={c:?})"
        );
    }
}

#[test]
fn insert_and_search_agree_on_basic_entries() {
    let k = 15;
    let mut rust = KmerIndex::new();
    let mut cpp = CppKmerIndex::new();

    let (kc_rust, kc_cpp) = encode_pair(k, b"ACGTACGTACGTACG");
    rust.insert(&kc_rust, 1, 0, 1);
    cpp.insert(&kc_cpp, 1, 0, 1);
    rust.insert(&kc_rust, 2, 5, 1);
    cpp.insert(&kc_cpp, 2, 5, 1);

    assert_search_matches("basic", rust.search(&kc_rust), &cpp.search(&kc_cpp));
}

#[test]
fn forward_keying_not_canonical_reverse_complement_pair_diverges() {
    // "AAAAAAAAAAAAAAC" (k=15) and its reverse complement
    // "GTTTTTTTTTTTTTT" are NOT palindromic, so their forward codes differ,
    // but their *canonical* codes (min of fwd/revcomp) would be identical.
    // KmerIndex keys on the forward code (GetCode()), so inserting distinct
    // entries under each must NOT collide -- proving forward, not canonical,
    // keying. If Insert/Search were mistakenly canonicalized, this test
    // would observe both entries merged into one bucket's search result.
    let k = 15;
    let fwd_bytes = b"AAAAAAAAAAAAAAC";
    let rc_bytes = b"GTTTTTTTTTTTTTT";

    // Sanity: confirm the two k-mers really are a non-identical
    // reverse-complement pair (same canonical code, different forward code).
    let (fwd_rust, fwd_cpp) = encode_pair(k, fwd_bytes);
    let (rc_rust, rc_cpp) = encode_pair(k, rc_bytes);
    assert_eq!(fwd_rust.get_canonical_kmer_code(), rc_rust.get_canonical_kmer_code());
    assert_ne!(fwd_rust.get_code(), rc_rust.get_code());
    assert_eq!(fwd_cpp.canonical(), rc_cpp.canonical());
    assert_ne!(fwd_cpp.get_code(), rc_cpp.get_code());

    let mut rust = KmerIndex::new();
    let mut cpp = CppKmerIndex::new();
    rust.insert(&fwd_rust, 100, 0, 1);
    cpp.insert(&fwd_cpp, 100, 0, 1);
    rust.insert(&rc_rust, 200, 0, 1);
    cpp.insert(&rc_cpp, 200, 0, 1);

    let rust_fwd_hits = rust.search(&fwd_rust).to_vec();
    let cpp_fwd_hits = cpp.search(&fwd_cpp);
    let rust_rc_hits = rust.search(&rc_rust).to_vec();
    let cpp_rc_hits = cpp.search(&rc_cpp);

    assert_search_matches("fwd bucket", &rust_fwd_hits, &cpp_fwd_hits);
    assert_search_matches("rc bucket", &rust_rc_hits, &cpp_rc_hits);

    // The two buckets must be DIFFERENT (each holds only its own entry),
    // proving Insert/Search do not collapse reverse-complement pairs.
    assert_eq!(rust_fwd_hits, vec![IndexInfo { idx: 100, offset: 0 }]);
    assert_eq!(rust_rc_hits, vec![IndexInfo { idx: 200, offset: 0 }]);
}

#[test]
fn invalid_kmer_is_dropped_by_insert_and_search() {
    let k = 15;
    // Contains 'N' in the middle -> IsValid() == false for the whole window.
    let (kc_rust, kc_cpp) = encode_pair(k, b"ACGTACGTNACGTAC");
    assert!(!kc_rust.is_valid());
    assert!(!kc_cpp.is_valid());

    let mut rust = KmerIndex::new();
    let mut cpp = CppKmerIndex::new();

    // Insert on an invalid k-mer must be a silent no-op on both sides.
    rust.insert(&kc_rust, 42, 7, 1);
    cpp.insert(&kc_cpp, 42, 7, 1);

    assert_search_matches("invalid-drop", rust.search(&kc_rust), &cpp.search(&kc_cpp));
    assert!(rust.search(&kc_rust).is_empty());

    // Remove on an invalid k-mer must also be a silent no-op (not a panic).
    rust.remove(&kc_rust, 42, 7, 1);
    cpp.remove(&kc_cpp, 42, 7, 1);
    assert_search_matches("invalid-remove-noop", rust.search(&kc_rust), &cpp.search(&kc_cpp));
}

#[test]
fn repeated_kmer_preserves_insertion_order() {
    let k = 15;
    let (kc_rust, kc_cpp) = encode_pair(k, b"GATTACAGATTACAG");

    let mut rust = KmerIndex::new();
    let mut cpp = CppKmerIndex::new();

    // Insert the SAME k-mer several times with different (idx, offset)
    // pairs, and once more with a differing `strand` value (which must have
    // zero effect on the stored entry -- the C++ `_indexInfo::strand` field
    // is commented out).
    let entries: [(u32, u32, i32); 5] = [(1, 0, 1), (2, 10, 1), (1, 0, -1), (3, 99, 0), (2, 10, 5)];
    for &(idx, offset, strand) in &entries {
        rust.insert(&kc_rust, idx, offset, strand);
        cpp.insert(&kc_cpp, idx, offset, strand);
    }

    assert_search_matches("insertion-order", rust.search(&kc_rust), &cpp.search(&kc_cpp));
    // Explicitly assert the order (not just size/set membership): five
    // entries in exactly the order inserted, `strand` values collapsed away.
    let expected: Vec<IndexInfo> =
        entries.iter().map(|&(idx, offset, _)| IndexInfo { idx, offset }).collect();
    assert_eq!(rust.search(&kc_rust), expected.as_slice());
}

#[test]
fn remove_matches_t1k_semantics() {
    let k = 15;
    let (kc_rust, kc_cpp) = encode_pair(k, b"TGCATGCATGCATGC");

    let mut rust = KmerIndex::new();
    let mut cpp = CppKmerIndex::new();

    for &(idx, offset) in &[(1u32, 0u32), (2, 0), (1, 0), (3, 7)] {
        rust.insert(&kc_rust, idx, offset, 1);
        cpp.insert(&kc_cpp, idx, offset, 1);
    }
    assert_search_matches("before-remove", rust.search(&kc_rust), &cpp.search(&kc_cpp));

    // Remove (1, 0): must remove only the FIRST matching entry, preserving
    // relative order of the rest (SimpleVector::Remove shifts left).
    rust.remove(&kc_rust, 1, 0, 1);
    cpp.remove(&kc_cpp, 1, 0, 1);
    assert_search_matches("after-first-remove", rust.search(&kc_rust), &cpp.search(&kc_cpp));

    // Remove a (idx, offset) pair that no longer/never matched -> no-op.
    rust.remove(&kc_rust, 99, 99, 1);
    cpp.remove(&kc_cpp, 99, 99, 1);
    assert_search_matches("noop-remove", rust.search(&kc_rust), &cpp.search(&kc_cpp));

    // Remove the remaining (1, 0) entry, then confirm it's actually gone.
    rust.remove(&kc_rust, 1, 0, 1);
    cpp.remove(&kc_cpp, 1, 0, 1);
    assert_search_matches("after-second-remove", rust.search(&kc_rust), &cpp.search(&kc_cpp));
}

/// A read long enough to cover k=15 and k=31, containing: a run of plain
/// non-repeating sequence, a long homopolymer run (to exercise
/// `BuildIndexFromRead`'s consecutive-duplicate dedup), and an `N` (to
/// exercise the `IsValid` gate mid-read).
fn build_index_read() -> &'static [u8] {
    b"GATTACAGATTACAGGGCTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTNNNNNACGTACGTACGTACGTGATTACAGATTACC"
}

/// Runs `BuildIndexFromRead`/`build_index_from_read` on both a Rust and a
/// C++ `KmerIndex` for k-mer length `k`, with the given `id`/`shift`, then
/// queries every length-`k` window of `read` (which, thanks to the read's
/// homopolymer run, includes windows that must have been *dropped* by the
/// dedup logic) and asserts the two indexes agree, in order, everywhere.
fn run_build_index_from_read_diff(k: usize, id: i32, shift: i32) {
    let read = build_index_read();
    assert!(read.len() >= k, "test read must be at least k={k} bases long");

    let mut rust_index = KmerIndex::new();
    let mut cpp_index = CppKmerIndex::new();
    let mut rust_kc = KmerCode::new(k);
    let mut cpp_kc = CppKmerCode::new(i32::try_from(k).unwrap());

    rust_index.build_index_from_read(&mut rust_kc, read, id, shift);
    cpp_index.build_index_from_read(&mut cpp_kc, read, id, shift);

    for start in 0..=(read.len() - k) {
        let window = &read[start..start + k];
        let (q_rust, q_cpp) = encode_pair(k, window);
        assert_search_matches(
            &format!(
                "build_index_from_read k={k} id={id} shift={shift} window@{start} ({window:?})"
            ),
            rust_index.search(&q_rust),
            &cpp_index.search(&q_cpp),
        );
    }
}

#[test]
fn build_index_from_read_matches_t1k_k15_positive_shift() {
    run_build_index_from_read_diff(15, 7, 3);
}

#[test]
fn build_index_from_read_matches_t1k_k31_positive_shift() {
    run_build_index_from_read_diff(31, 12, 100);
}

#[test]
fn build_index_from_read_matches_t1k_k15_zero_shift_and_id() {
    run_build_index_from_read_diff(15, 0, 0);
}

#[test]
fn build_index_from_read_matches_t1k_k15_negative_shift_and_id() {
    // Negative `shift`/`id`: BuildIndexFromRead computes `i - kl + 1 + shift`
    // as a plain (signed) `int` and implicitly converts it (and `id`) to
    // `index_t` (`uint32_t`) when calling `Insert` -- i.e. a two's-complement
    // wraparound, not a panic/clamp. This must match bit-for-bit between the
    // Rust port (`as u32`) and the real C++ implicit conversion.
    run_build_index_from_read_diff(15, -5, -20);
}

#[test]
fn build_index_from_read_matches_t1k_k31_negative_shift_and_id() {
    run_build_index_from_read_diff(31, -1, -1);
}

#[test]
fn kmerindex_new_produces_a_usable_handle() {
    // Basic smoke test that the opaque handle itself is non-null-usable (a
    // NULL handle would panic inside CppKmerIndex::new per its contract).
    let cpp = CppKmerIndex::new();
    let kc = CppKmerCode::new(9);
    assert!(cpp.search(&kc).is_empty());
}
