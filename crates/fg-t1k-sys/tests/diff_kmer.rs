#![cfg(feature = "t1k-sys")]
use fg_t1k_core::kmer::canonical_kmer;

/// Runs `canonical_kmer` (Rust port) and `fg_t1k_canonical_kmer` (C++ FFI
/// shim) on the same `(seq, k)` input and asserts the two implementations
/// agree bit-for-bit.
fn assert_kmer_matches(seq: &str, k: usize) {
    let bytes = seq.as_bytes();
    let rust = canonical_kmer(bytes, k);
    let cpp = unsafe {
        fg_t1k_sys::fg_t1k_canonical_kmer(
            bytes.as_ptr().cast::<std::os::raw::c_char>(),
            i32::try_from(bytes.len()).unwrap(),
            i32::try_from(k).unwrap(),
        )
    };
    assert_eq!(rust, cpp, "seq={seq} k={k}");
}

#[test]
fn canonical_kmer_matches_t1k() {
    let k = 9usize;
    // Original uppercase-ACGT cases, len == k.
    for seq in ["ACGTACGTA", "TTTTAAAAC", "GATTACAGA", "CCCGGGAAA"] {
        assert_kmer_matches(seq, k);
    }

    // `N`-containing sequence: exercises the `nucToNum[N] & 3 == 3` quirk
    // (an invalid table entry is treated the same as `T`).
    assert_kmer_matches("ACGTNACGT", k);

    // Rolling window: len > k, so only the last k bases contribute.
    assert_kmer_matches("ACGTACGTACGT", 9);

    // Short input: len < k.
    assert_kmer_matches("ACG", 9);

    // k boundaries: k = 1 (minimum useful window) and k = 31 (largest k that
    // still fits comfortably below the k = 32 packing limit).
    assert_kmer_matches("ACGTACGTA", 1);
    assert_kmer_matches("ACGTACGTACGTACGTACGTACGTACGTACGT", 31);
}
