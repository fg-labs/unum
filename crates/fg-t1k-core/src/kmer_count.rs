//! Canonical k-mer counting, ported from T1K's `KmerCount`
//! (`vendor/t1k/KmerCount.hpp`).
//!
//! The C++ class stores counts in `std::map<uint64_t, int> count[103]`,
//! bucketed by `GetHash(k) = k % 103` (`KmerCount.hpp:24-27`). That bucketing
//! is purely an implementation detail for `std::map` lookup performance under
//! the hood -- the count for a given canonical k-mer code `k` is always
//! `count[k % 103][k]`, so a single `HashMap<u64, i32>` keyed directly by the
//! canonical code yields identical per-k-mer counts without replicating the
//! 103-way bucket array. `Output`/`AddCountFromFile`
//! (`KmerCount.hpp:83-126`), whose serialized byte order depends on the
//! bucket iteration order, are intentionally NOT ported here -- byte-identical
//! `Output` is deferred until a consumer actually needs it.

use crate::kmer::KmerCode;
use std::collections::HashMap;

/// Counts canonical k-mers across reads, ported from T1K's `KmerCount`
/// (`vendor/t1k/KmerCount.hpp`).
///
/// Unlike the vendored C++ (which stores an owned `KmerCode` member reused
/// across calls via `Restart()`), this port constructs a fresh, local
/// [`KmerCode`] inside [`KmerCount::add_count`] and [`KmerCount::get_count`].
/// `KmerCode::new` and `KmerCode::Restart` produce byte-identical state
/// (`code = 0`, `invalid_pos = -1`), so this is behaviorally identical while
/// letting [`KmerCount::get_count`] take `&self` rather than `&mut self`.
#[derive(Debug, Clone)]
pub struct KmerCount {
    /// Per-canonical-k-mer counts. Equivalent to the union of all 103 buckets
    /// in the C++ `count[103]` array (see module docs).
    counts: HashMap<u64, i32>,
    /// The k-mer length; mirrors the C++ `kmerLength` field.
    kmer_length: usize,
}

impl KmerCount {
    /// Constructs an empty `KmerCount` for k-mer length `k`. Mirrors the
    /// `KmerCount(int k)` constructor (`KmerCount.hpp:29-35`), minus the
    /// `maxReadLen`/`c` scratch-buffer fields, which are internal
    /// bookkeeping not observable through `AddCount`/`GetCount`/Jaccard.
    #[must_use]
    pub fn new(k: usize) -> Self {
        Self { counts: HashMap::new(), kmer_length: k }
    }

    /// Ported from `KmerCount::AddCount` (`KmerCount.hpp:53-81`).
    ///
    /// Rolls a [`KmerCode`] across `read`: if `read` is shorter than the
    /// k-mer length, returns `0` without recording anything. Otherwise the
    /// first `kmer_length - 1` bases only fill the rolling window (no
    /// counting), and starting at index `kmer_length - 1` every full window
    /// is counted -- but only if [`KmerCode::is_valid`] (i.e. no `N` within
    /// the current window) -- by incrementing the entry for its canonical
    /// k-mer code. Returns `1` on success, matching the C++ return value
    /// semantics (`0` = read too short to yield any k-mer, `1` = processed).
    ///
    /// Uses `read.len()` as the read length, whereas T1K's C++ `AddCount`
    /// computes it via `strlen(read)`. The two diverge only if `read`
    /// contains an interior NUL byte (`strlen` would stop early; this port
    /// would keep going) -- which does not occur in real FASTQ/sequence
    /// data, and which the FFI wrapper (`CppKmerCount::add_count` in
    /// `fg-t1k-sys`) rejects up front via `CString::new`, since a NUL byte
    /// cannot be encoded in a NUL-terminated C string.
    pub fn add_count(&mut self, read: &[u8]) -> i32 {
        let len = read.len();
        if len < self.kmer_length {
            return 0;
        }

        let mut kmer_code = KmerCode::new(self.kmer_length);
        // Fill the first `kmer_length - 1` bases: these only build up the
        // rolling window and are never counted on their own (mirrors the
        // C++ `for (i = 0; i < kmerLength - 1; ++i) kmerCode.Append(...)`
        // loop, which runs zero times when `kmer_length == 0`).
        let fill = self.kmer_length.saturating_sub(1);
        let mut i = 0usize;
        while i < fill {
            kmer_code.append(read[i]);
            i += 1;
        }

        // From here on, every append completes a full k-length window;
        // count it if (and only if) the window is valid (no `N` inside it).
        while i < len {
            kmer_code.append(read[i]);
            if kmer_code.is_valid() {
                let kcode = kmer_code.get_canonical_kmer_code();
                *self.counts.entry(kcode).or_insert(0) += 1;
            }
            i += 1;
        }

        1
    }

    /// Ported from `KmerCount::GetCount` (`KmerCount.hpp:135-152`).
    ///
    /// Encodes exactly the first `kmer_length` bytes of `kmer` into a fresh
    /// [`KmerCode`] via [`KmerCode::append`] (NOT [`KmerCode::prepend`] --
    /// the C++ builds the query k-mer the same way `AddCount` builds each
    /// window, appending left-to-right), then looks up its canonical code.
    /// If the query k-mer is invalid (contains an `N`), or was never
    /// recorded via [`KmerCount::add_count`], returns `0`.
    ///
    /// # Panics
    ///
    /// Panics if `kmer.len() < kmer_length` (slice indexing out of bounds).
    /// The C++ has no such guard: `GetCount` unconditionally reads
    /// `kmer[0..kmerLength]` regardless of the buffer's actual length, which
    /// is undefined behavior for a too-short input. This port turns that UB
    /// into a deterministic panic instead of silently reading out of bounds.
    #[must_use]
    pub fn get_count(&self, kmer: &[u8]) -> i32 {
        let mut kmer_code = KmerCode::new(self.kmer_length);
        for &c in &kmer[..self.kmer_length] {
            kmer_code.append(c);
        }
        if kmer_code.is_valid() {
            let kcode = kmer_code.get_canonical_kmer_code();
            *self.counts.get(&kcode).unwrap_or(&0)
        } else {
            0
        }
    }

    /// Ported from `KmerCount::GetCountSimilarityJaccard`
    /// (`KmerCount.hpp:166-192`).
    ///
    /// `countA`/`countB` are the total k-mer occurrence counts in `self` and
    /// `other` respectively (summed over every distinct k-mer, not the
    /// number of distinct k-mers); `sharedCount` sums, for every k-mer
    /// present in both, `min(count_in_self, count_in_other)`. The result is
    /// `sharedCount / (countA + countB - sharedCount)`.
    ///
    /// Uses `i32` accumulators (matching the C++ `int countA, countB,
    /// sharedCount`) rather than a wider integer type, so the arithmetic --
    /// and hence the final `f64` division -- is bit-for-bit identical to the
    /// vendored C++, including in the (untested, and not expected in
    /// practice) case of `i32` overflow.
    #[must_use]
    pub fn jaccard_similarity(&self, other: &KmerCount) -> f64 {
        let mut count_a: i32 = 0;
        let mut count_b: i32 = 0;
        let mut shared_count: i32 = 0;

        for (kcode, &c) in &self.counts {
            count_a += c;
            if let Some(&other_c) = other.counts.get(kcode) {
                shared_count += c.min(other_c);
            }
        }
        for &c in other.counts.values() {
            count_b += c;
        }

        f64::from(shared_count) / f64::from(count_a + count_b - shared_count)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_count_returns_zero_for_short_read() {
        let mut kc = KmerCount::new(5);
        assert_eq!(kc.add_count(b"ACG"), 0);
        // No k-mers should have been recorded.
        assert_eq!(kc.get_count(b"ACGAC"), 0);
    }

    #[test]
    fn add_count_returns_one_and_counts_one_kmer_for_exact_length_read() {
        let mut kc = KmerCount::new(4);
        assert_eq!(kc.add_count(b"ACGT"), 1);
        // "ACGT" is its own reverse complement, so canonical == "ACGT".
        assert_eq!(kc.get_count(b"ACGT"), 1);
    }

    #[test]
    fn add_count_counts_every_rolling_window() {
        // len=6, k=4 => 3 windows: ACGT, CGTA, GTAC.
        let mut kc = KmerCount::new(4);
        assert_eq!(kc.add_count(b"ACGTAC"), 1);
        assert_eq!(kc.get_count(b"ACGT"), 1);
        assert_eq!(kc.get_count(b"CGTA"), 1);
        assert_eq!(kc.get_count(b"GTAC"), 1);
        // A window never seen must be absent (count 0).
        assert_eq!(kc.get_count(b"TTTT"), 0);
    }

    #[test]
    fn add_count_skips_windows_containing_n() {
        // k=3 over "ANCGA": windows are [A,N,C] (invalid), [N,C,G] (invalid),
        // [C,G,A] (valid, once N has aged out of the 3-length window).
        let mut kc = KmerCount::new(3);
        assert_eq!(kc.add_count(b"ANCGA"), 1);
        // Only the one valid window ("CGA") was recorded.
        assert_eq!(kc.get_count(b"CGA"), 1);
        // A different, never-seen window must be absent.
        assert_eq!(kc.get_count(b"TTT"), 0);
    }

    #[test]
    fn add_count_accumulates_across_repeated_reads() {
        let mut kc = KmerCount::new(4);
        kc.add_count(b"ACGT");
        kc.add_count(b"ACGT");
        kc.add_count(b"ACGT");
        assert_eq!(kc.get_count(b"ACGT"), 3);
    }

    #[test]
    fn add_count_canonicalizes_reverse_complement_duplicates() {
        // "AAAC" and its reverse complement "GTTT" must count as the same
        // canonical k-mer.
        let mut kc = KmerCount::new(4);
        kc.add_count(b"AAAC");
        kc.add_count(b"GTTT");
        assert_eq!(kc.get_count(b"AAAC"), 2);
        assert_eq!(kc.get_count(b"GTTT"), 2);
    }

    #[test]
    fn get_count_returns_zero_for_query_containing_n() {
        let mut kc = KmerCount::new(4);
        kc.add_count(b"ACGTACGT");
        assert_eq!(kc.get_count(b"ACNT"), 0);
    }

    #[test]
    // These hand-computed Jaccard values are exact integer ratios with
    // exactly-representable binary fractions (1.0, 0.25, 0.0), so strict
    // float equality is intentional here, not a rounding-tolerance bug.
    #[allow(clippy::float_cmp)]
    fn jaccard_similarity_of_identical_sets_is_one() {
        let mut a = KmerCount::new(4);
        a.add_count(b"ACGTAC");
        let mut b = KmerCount::new(4);
        b.add_count(b"ACGTAC");
        assert_eq!(a.jaccard_similarity(&b), 1.0);
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn jaccard_similarity_hand_computed() {
        // a: ACGT x2, CGTA x1 (counts: {ACGT:2, CGTA:1}), countA = 3.
        let mut a = KmerCount::new(4);
        a.add_count(b"ACGT");
        a.add_count(b"ACGT");
        a.add_count(b"CGTA");
        // b: ACGT x1, GGGG x1 (counts: {ACGT:1, GGGG:1}), countB = 2.
        let mut b = KmerCount::new(4);
        b.add_count(b"ACGT");
        b.add_count(b"GGGG");
        // shared = min(2,1) for ACGT = 1; CGTA/GGGG not shared.
        // jaccard = 1 / (3 + 2 - 1) = 1/4 = 0.25.
        assert_eq!(a.jaccard_similarity(&b), 0.25);
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn jaccard_similarity_disjoint_sets_is_zero() {
        let mut a = KmerCount::new(4);
        a.add_count(b"AAAA");
        let mut b = KmerCount::new(4);
        b.add_count(b"CCCC");
        assert_eq!(a.jaccard_similarity(&b), 0.0);
    }
}
