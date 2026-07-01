#![cfg(feature = "t1k-sys")]
use fg_t1k_core::kmer::canonical_kmer;

#[test]
fn canonical_kmer_matches_t1k() {
    let k = 9usize;
    for seq in ["ACGTACGTA", "TTTTAAAAC", "GATTACAGA", "CCCGGGAAA"] {
        let bytes = seq.as_bytes();
        let rust = canonical_kmer(bytes, k);
        let cpp = unsafe {
            fg_t1k_sys::fg_t1k_canonical_kmer(
                bytes.as_ptr() as *const _,
                bytes.len() as i32,
                k as i32,
            )
        };
        assert_eq!(rust, cpp, "seq={seq}");
    }
}
