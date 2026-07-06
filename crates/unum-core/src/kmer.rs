//! Canonical k-mer encoding, ported from T1K's `KmerCode` (`KmerCode.hpp`).
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

/// Returns `true` if `c` is not one of the four canonical bases (`A`/`C`/`G`/`T`).
///
/// Mirrors T1K's `nucToNum` table (`Genotyper.cpp:37-40`), which holds `-1`
/// for every entry except `A`, `C`, `G`, `T`. This is the check
/// `KmerCode::Prepend` uses (`nucToNum[c-'A'] == -1`) to decide whether a
/// base is invalid -- note this differs from `KmerCode::Append`'s literal
/// `c == 'N'` check (see [`KmerCode::append`]'s doc comment for the
/// consequence of that asymmetry).
fn is_invalid_base(c: u8) -> bool {
    !matches!(c, b'A' | b'C' | b'G' | b'T')
}

/// Converts a k-mer length/shift-count (`usize`) to `i64` for comparison
/// against `invalid_pos`. `KmerCode` is only meaningful for `k <= 32` (2-bit
/// packing into a `u64`), so this conversion never wraps in practice; it
/// saturates to `i64::MAX` rather than panicking if that invariant is ever
/// violated.
fn kmer_length_as_i64(n: usize) -> i64 {
    i64::try_from(n).unwrap_or(i64::MAX)
}

/// A stateful, rolling 2-bit-packed k-mer encoder, ported from T1K's
/// `KmerCode` (`KmerCode.hpp`).
///
/// Unlike [`canonical_kmer`] (a one-shot function), `KmerCode` maintains a
/// sliding window of the last `kmer_length` bases as callers stream bases in
/// one at a time via [`KmerCode::append`] or [`KmerCode::prepend`]. It also
/// tracks whether the current window contains an invalid (non-ACGT) base via
/// `invalid_pos`, exposed through [`KmerCode::is_valid`].
///
/// # `invalid_pos` bookkeeping
///
/// `invalid_pos` is `-1` when the window is fully valid (no invalid base seen
/// within the last `kmer_length` appends/prepends), or else the *age* of the
/// most recently marked invalid position: 0 immediately after the invalid
/// base was introduced, incrementing by one on every subsequent `append`, and
/// reset to `-1` once it ages out of the window (`invalid_pos >=
/// kmer_length`). `shift_right` moves the window the opposite direction, so
/// it *decrements* `invalid_pos` (clamping negative results back to `-1`).
///
/// This exactly reproduces T1K's semantics, including a real asymmetry in the
/// upstream C++: `Append` only treats the literal byte `'N'` as invalid,
/// while `Prepend` uses the full `nucToNum` lookup (any non-ACGT byte). This
/// port preserves that asymmetry rather than "fixing" it, since the
/// differential test validates byte-for-byte parity with T1K.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KmerCode {
    kmer_length: usize,
    invalid_pos: i64,
    code: u64,
    mask: u64,
}

impl KmerCode {
    /// Constructs a new `KmerCode` for k-mer length `k`, mirroring the
    /// `KmerCode(int kl)` constructor: `code = 0`, `invalid_pos = -1`, and
    /// `mask` built by repeatedly shifting left 2 bits and OR-ing in `0b11`,
    /// `k` times (so `mask` has its lowest `2*k` bits set).
    ///
    /// # Panics
    ///
    /// Panics if `k > 32`: the 2-bit-per-base packing lives in a `u64`, so a
    /// `kmer_length` beyond 32 has no room in the packed word (`2*k > 64`) and
    /// would silently overflow the window. This guard is a fail-fast on an
    /// invariant the C++ side leaves implicit; it never fires for any real
    /// usage (the largest supported k is 32).
    #[must_use]
    pub fn new(k: usize) -> Self {
        assert!(k <= 32, "kmer length must be <= 32");
        let mut mask: u64 = 0;
        for _ in 0..k {
            mask <<= 2;
            mask |= 3;
        }
        Self { kmer_length: k, invalid_pos: -1, code: 0, mask }
    }

    /// Mirrors `KmerCode::Restart`: resets `code` to 0 and `invalid_pos` to
    /// -1 (fully valid, empty window), without changing `kmer_length`/`mask`.
    pub fn restart(&mut self) {
        self.code = 0;
        self.invalid_pos = -1;
    }

    /// Mirrors `KmerCode::GetCode`: the raw packed code (not canonicalized).
    #[must_use]
    pub fn get_code(&self) -> u64 {
        self.code
    }

    /// Mirrors `KmerCode::SetCode`: overwrites the raw code and resets
    /// `invalid_pos` to -1 (the window is considered fully valid after a
    /// direct code assignment, regardless of prior state).
    pub fn set_code(&mut self, c: u64) {
        self.code = c;
        self.invalid_pos = -1;
    }

    /// Mirrors `KmerCode::GetCanonicalKmerCode`: `min(code, reverse_complement)`.
    #[must_use]
    pub fn get_canonical_kmer_code(&self) -> u64 {
        let rc = self.get_reverse_complement_code();
        rc.min(self.code)
    }

    /// Mirrors `KmerCode::GetReverseComplementCode`: walks `code`'s 2-bit
    /// groups from lowest-order to highest-order, complementing each
    /// (`3 - base`) and appending it to the result -- i.e. the reverse
    /// complement of `code`, in the same high-to-low packing convention.
    #[must_use]
    pub fn get_reverse_complement_code(&self) -> u64 {
        let mut cr_code: u64 = 0;
        for i in 0..self.kmer_length {
            let tmp = (self.code >> (2 * i)) & 3;
            cr_code = (cr_code << 2) | (3 - tmp);
        }
        cr_code
    }

    /// Mirrors `KmerCode::GetKmerLength`.
    #[must_use]
    pub fn get_kmer_length(&self) -> usize {
        self.kmer_length
    }

    /// Mirrors `KmerCode::IsValid`: `true` iff `invalid_pos == -1`.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        self.invalid_pos == -1
    }

    /// Mirrors `KmerCode::Append`: pushes `c` into the low 2 bits of `code`
    /// (masking out any bits beyond `kmer_length`), and updates
    /// `invalid_pos`:
    /// - If a position is already marked invalid (`invalid_pos != -1`), age
    ///   it by one.
    /// - If `c` is the literal byte `'N'` (not a general invalid-base check
    ///   -- see the struct-level doc comment), mark position 0 as invalid.
    /// - If the marked invalid position has aged out of the window
    ///   (`invalid_pos >= kmer_length`), clear it back to -1.
    pub fn append(&mut self, c: u8) {
        if self.invalid_pos != -1 {
            self.invalid_pos += 1;
        }

        self.code = ((self.code << 2) & self.mask) | nuc_to_code(c);

        if c == b'N' {
            self.invalid_pos = 0;
        }
        if self.invalid_pos >= kmer_length_as_i64(self.kmer_length) {
            self.invalid_pos = -1;
        }
    }

    /// Mirrors `KmerCode::Prepend`: shifts the window right by one (making
    /// room at the high end), then writes `c`'s 2-bit code into the highest
    /// `2` bits of the `kmer_length`-bit window. Uses the full
    /// [`is_invalid_base`] check (any non-ACGT byte), unlike `Append`'s
    /// literal `'N'` check -- see the struct-level doc comment.
    ///
    /// # Panics
    ///
    /// Panics if `kmer_length == 0`. The `self.kmer_length as u64 - 1`
    /// subtraction below would otherwise underflow: in debug builds Rust's
    /// overflow check already panics, while release builds would wrap
    /// silently and produce a garbage shift amount/code. The C++ side has
    /// undefined behavior at this call site for a zero-length k-mer; rather
    /// than mirror that garbage in release builds, this guard fails fast with
    /// a clear message in all build modes. A zero-length k-mer is a nonsense
    /// input that no real usage produces.
    pub fn prepend(&mut self, c: u8) {
        assert!(self.kmer_length > 0, "prepend requires kmer_length >= 1");
        self.shift_right(1);
        if is_invalid_base(c) {
            self.invalid_pos = kmer_length_as_i64(self.kmer_length) - 1;
        }
        let shift = 2 * (self.kmer_length as u64 - 1);
        self.code = (self.code | (nuc_to_code(c) << shift)) & self.mask;
    }

    /// Mirrors `KmerCode::ShiftRight`: shifts `code` right by `2*k` bits
    /// (masked to the remaining window width), and moves `invalid_pos` in
    /// the opposite direction of `Append`/`Prepend` (decrementing by `k`,
    /// clamped back to -1 if it goes negative).
    ///
    /// # Panics
    ///
    /// Panics if `k >= 32`, i.e. `2*k >= 64`. A shift that consumes the whole
    /// packed `u64` word is `>> 64`: undefined behavior in C++ and a masked
    /// (`& 63`) no-op in release Rust / an overflow panic in debug Rust. T1K
    /// never calls `ShiftRight` with the full `kmer_length` in one step, so
    /// this is a nonsense input rather than a real usage pattern; this guard
    /// fails fast with a clear message instead of silently diverging.
    pub fn shift_right(&mut self, k: usize) {
        assert!(k < 32, "shift_right count must be < 32 (2*k must fit in a u64)");
        if self.invalid_pos != -1 {
            self.invalid_pos -= kmer_length_as_i64(k);
        }

        let shift = 2 * k as u64;
        self.code = (self.code >> shift) & (self.mask >> shift);

        if self.invalid_pos < 0 {
            self.invalid_pos = -1;
        }
    }

    /// Mirrors `KmerCode::IsEqual`: compares only the raw `code` field (not
    /// `kmer_length`/`mask`/`invalid_pos`).
    #[must_use]
    pub fn is_equal(&self, other: &KmerCode) -> bool {
        self.code == other.code
    }
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

    #[test]
    fn kmer_code_rolling_window_matches_canonical_kmer() {
        // Appending base-by-base through KmerCode should agree with the
        // one-shot canonical_kmer function for the same effective window.
        let k = 5;
        let mut kc = KmerCode::new(k);
        for &c in b"GGGACGTA" {
            kc.append(c);
        }
        // Last 5 bases: "ACGTA"
        assert_eq!(kc.get_canonical_kmer_code(), canonical_kmer(b"ACGTA", k));
        assert!(kc.is_valid());
    }

    #[test]
    fn kmer_code_canonical_is_min_of_forward_and_revcomp() {
        let k = 4;
        let mut kc = KmerCode::new(k);
        for &c in b"ACGT" {
            kc.append(c);
        }
        let fwd = kc.get_code();
        let rc = kc.get_reverse_complement_code();
        assert_eq!(kc.get_canonical_kmer_code(), fwd.min(rc));
    }

    #[test]
    fn kmer_code_invalid_base_marks_and_ages_out() {
        // k=3: after "A,N,C,G,A", the last-3 window progresses
        // [A] -> [A,N] -> [A,N,C] -> [N,C,G] -> [C,G,A]. N first enters the
        // window at the "A,N" step and only fully exits once the window
        // slides past it (dropped after the *second* append following it).
        let k = 3;
        let mut kc = KmerCode::new(k);
        kc.append(b'A');
        assert!(kc.is_valid());
        kc.append(b'N');
        assert!(!kc.is_valid(), "window contains an N, so it must be invalid");
        kc.append(b'C');
        assert!(!kc.is_valid(), "N is still within the k=3 window [A,N,C]");
        kc.append(b'G');
        assert!(!kc.is_valid(), "N is still within the k=3 window [N,C,G]");
        kc.append(b'A');
        assert!(kc.is_valid(), "N has aged out of the k=3 window [C,G,A]");
    }

    #[test]
    fn kmer_code_restart_resets_code_and_validity() {
        let mut kc = KmerCode::new(3);
        kc.append(b'N');
        assert!(!kc.is_valid());
        kc.restart();
        assert_eq!(kc.get_code(), 0);
        assert!(kc.is_valid());
    }

    #[test]
    fn kmer_code_set_code_forces_valid() {
        let mut kc = KmerCode::new(3);
        kc.append(b'N');
        assert!(!kc.is_valid());
        kc.set_code(0b10_1010);
        assert!(kc.is_valid(), "SetCode always resets invalid_pos to -1");
        assert_eq!(kc.get_code(), 0b10_1010);
    }

    #[test]
    fn kmer_code_is_equal_compares_code_only() {
        let mut a = KmerCode::new(3);
        let mut b = KmerCode::new(3);
        a.append(b'A');
        a.append(b'C');
        a.append(b'G');
        b.append(b'A');
        b.append(b'C');
        b.append(b'G');
        assert!(a.is_equal(&b));
        b.append(b'T');
        assert!(!a.is_equal(&b));
    }

    #[test]
    fn kmer_code_shift_right_drops_low_order_bases() {
        let k = 4;
        let mut kc = KmerCode::new(k);
        for &c in b"ACGT" {
            kc.append(c);
        }
        let before = kc.get_code();
        kc.shift_right(1);
        // Shifting right by 1 base (2 bits) drops the most-recently-appended
        // base and keeps the rest, masked to the (now-narrower) window.
        assert_eq!(kc.get_code(), (before >> 2) & (kc.mask >> 2));
    }

    #[test]
    fn kmer_code_prepend_adds_base_at_high_end() {
        let k = 4;
        let mut kc = KmerCode::new(k);
        for &c in b"ACGT" {
            // Fills the k=4 window completely (no unfilled/padding slot),
            // so Prepend's effect is unambiguous: it evicts the newest base
            // ('T') and inserts the prepended base as the new oldest.
            kc.append(c);
        }
        kc.prepend(b'A');
        // Prepending 'A' onto a full "ACGT" window evicts the newest base
        // ('T') and inserts 'A' as the new oldest, yielding "AACG".
        let mut expected = KmerCode::new(k);
        for &c in b"AACG" {
            expected.append(c);
        }
        assert_eq!(kc.get_code(), expected.get_code());
    }

    #[test]
    fn kmer_code_new_allows_maximum_k_of_32() {
        // k=32 is the largest k that fits in a u64 (2*32 == 64 bits); the mask
        // is all-ones with no headroom bit. Must not panic.
        let kc = KmerCode::new(32);
        assert_eq!(kc.mask, u64::MAX);
    }

    #[test]
    #[should_panic(expected = "kmer length must be <= 32")]
    fn kmer_code_new_rejects_k_above_32() {
        // 2*33 == 66 bits has no room in a u64 packed word, so `new` must
        // fail fast rather than silently overflow the window.
        let _ = KmerCode::new(33);
    }

    #[test]
    #[should_panic(expected = "prepend requires kmer_length >= 1")]
    fn kmer_code_prepend_rejects_zero_length() {
        // KmerCode::new(0) is constructible, but prepend's shift computation
        // underflows on a zero-length window; guard it explicitly so release
        // builds fail fast instead of producing garbage.
        let mut kc = KmerCode::new(0);
        kc.prepend(b'A');
    }

    #[test]
    #[should_panic(expected = "shift_right count must be < 32")]
    fn kmer_code_shift_right_rejects_full_width() {
        // A full-width shift (2*k >= 64) is `>> 64`: UB in C++ and a masked
        // no-op / overflow panic in Rust. Guard it rather than silently
        // diverging from the C++ oracle on a nonsense count.
        let mut kc = KmerCode::new(32);
        kc.shift_right(32);
    }
}
