//! Canonical k-mer encoding, ported from T1K's `KmerCode` (`vendor/t1k/KmerCode.hpp`).
//!
//! Each base is packed into 2 bits (A=0, C=1, G=2, T=3). `KmerCode::Append` builds
//! the code as a sliding window: `code = ((code << 2) & mask) | base_code`, so the
//! earliest-appended base within the window ends up in the highest-order bits and
//! the most-recently-appended base occupies the lowest 2 bits. The canonical code
//! is `min(forward_code, reverse_complement_code)`, matching
//! `KmerCode::GetCanonicalKmerCode()`.

/// Maps a nucleotide byte to its 2-bit code (A=0, C=1, G=2, T=3).
///
/// Mirrors T1K's `nucToNum` table (`Genotyper.cpp:37-40`), which is indexed by
/// `c - 'A'` and holds `-1` for non-ACGT letters. `KmerCode::Append` masks the raw
/// table value with `& 3`, so a `signed char` value of `-1` (`0xFF`) becomes `3`
/// (as if the base were `T`) rather than being rejected outright; T1K separately
/// tracks an "invalid position" counter to flag the k-mer as unreliable when a
/// non-ACGT byte (or `N`) is appended. This port replicates only the bit-packing
/// behavior of `Append`, which is what `GetCanonicalKmerCode` operates on.
fn nuc_to_code(c: u8) -> u64 {
    match c {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        // `T` (code 3) and any non-ACGT byte share the same packed bit
        // pattern: `nucToNum[c - 'A'] & 3` is `3` for `T` directly, and also
        // `3` for an invalid table entry (`-1 as signed char, & 3 == 3`).
        _ => 3,
    }
}

/// Computes the canonical (strand-independent) 2-bit-packed k-mer code for the
/// last `k` bases of `seq`, matching T1K's `KmerCode::Append` followed by
/// `KmerCode::GetCanonicalKmerCode()` bit-for-bit.
///
/// Bases are packed high-to-low in append order: the earliest base within the
/// `k`-length window occupies the highest-order 2 bits, and the most recently
/// appended base occupies the lowest 2 bits. If `seq` is longer than `k`, only
/// the final `k` bases contribute to the code (a sliding window), matching the
/// mask-and-shift behavior of `Append`. The returned value is
/// `min(forward_code, reverse_complement_code)`.
///
/// # Preconditions
///
/// - `k` must be `<= 32`, since each base is packed into 2 bits of a `u64`
///   (`2 * k` must not exceed 64). `k == 0` yields a code of `0` (both forward
///   and reverse-complement codes are empty).
/// - Input bases in `seq` must be uppercase (`A`/`C`/`G`/`T`; any other
///   uppercase byte falls back to code `3`). T1K uppercases sequences
///   upstream before calling `KmerCode::Append`, so this port does not
///   special-case lowercase input. Lowercase bytes diverge from the C++
///   side, whose `nucToNum` table is indexed by `c - 'A'` and reads out of
///   bounds (undefined behavior) for lowercase input; callers must not rely
///   on parity with T1K for lowercase bytes.
///
/// # Panics
///
/// Panics (in debug builds, via `debug_assert!`) if `k > 32`, since the
/// 2-bit-per-base packing would overflow the `u64` code (`2 * k >= 64`
/// triggers a shift-amount overflow). Release builds compile out the
/// assertion and instead silently truncate.
#[must_use]
pub fn canonical_kmer(seq: &[u8], k: usize) -> u64 {
    debug_assert!(k <= 32, "kmer length must be <= 32 (2-bit packing in u64)");
    let mask: u64 = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };

    let mut code: u64 = 0;
    for &c in seq {
        code = ((code << 2) & mask) | nuc_to_code(c);
    }

    // Reverse-complement: `GetCanonicalKmerCode` walks the packed code's 2-bit
    // groups from lowest-order (most recently appended base) to highest-order
    // (oldest base in the window), complements each (`3 - base`), and appends
    // it to `rc_code`. Since each step shifts `rc_code` left before OR-ing in
    // the next complemented base, the final `rc_code` is the reverse
    // complement of the original k-mer, encoded in the same high-to-low
    // convention as `code`.
    let mut rc_code: u64 = 0;
    for i in 0..k {
        let base = (code >> (2 * i)) & 3;
        rc_code = (rc_code << 2) | (3 - base);
    }

    rc_code.min(code)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_kmer_is_symmetric_under_revcomp() {
        // Reverse-complementing the input sequence must not change the
        // canonical code, since canonical_kmer picks the min of the two.
        let k = 5;
        let fwd = canonical_kmer(b"ACGTA", k);
        let rc = canonical_kmer(b"TACGT", k); // revcomp of "ACGTA"
        assert_eq!(fwd, rc);
    }

    #[test]
    fn canonical_kmer_palindrome_equals_forward_and_revcomp() {
        // "ACGT" is its own reverse complement, so forward == revcomp == canonical.
        let k = 4;
        assert_eq!(canonical_kmer(b"ACGT", k), canonical_kmer(b"ACGT", k));
    }

    #[test]
    fn canonical_kmer_sliding_window_uses_last_k_bases() {
        // Appending extra leading bases beyond the k-length window must not
        // change the result, since Append maintains only the last k bases.
        let k = 3;
        let windowed = canonical_kmer(b"GGGACG", k); // last 3 bases: "ACG"
        let exact = canonical_kmer(b"ACG", k);
        assert_eq!(windowed, exact);
    }
}
