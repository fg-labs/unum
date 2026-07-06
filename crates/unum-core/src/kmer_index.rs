//! K-mer position index (k-mer code -> `(idx, offset)` list), ported from
//! T1K's `KmerIndex` (`vendor/t1k/KmerIndex.hpp`).
//!
//! # Forward vs. canonical keying
//!
//! Unlike [`crate::kmer_count::KmerCount`] (which keys on
//! `KmerCode::GetCanonicalKmerCode()`, the strand-independent `min(forward,
//! revcomp)` code), `KmerIndex::Insert`/`Search` key on the raw *forward*
//! `KmerCode::GetCode()` (`KmerIndex.hpp:66,97`). A k-mer and its reverse
//! complement therefore land in different buckets unless the k-mer happens
//! to be a palindrome. This is intentional in the upstream C++ (the caller,
//! `BuildIndexFromRead`, is responsible for indexing both strands if it
//! wants both -- and in fact it does not: see below).
//!
//! # `_indexInfo` has no `strand` field
//!
//! The C++ `struct _indexInfo` (`KmerIndex.hpp:12-17`) declares `idx` and
//! `offset` (both `index_t`); a `strand` field is present in a comment
//! (`//int strand ;`) but not compiled in. `Insert`/`Remove` both accept a
//! `strand: int` parameter, but neither reads nor stores it -- it is dead on
//! arrival. This port mirrors that exactly: [`IndexInfo`] has no `strand`
//! field, and `insert`/`remove`'s `strand` parameter is accepted (to keep
//! the call signature honest about what the C++ API looks like) and
//! ignored.
//!
//! # Bucketing is an implementation detail, not ported
//!
//! The C++ backing store is `std::map<uint64_t, SimpleVector<_indexInfo>>
//! index[KINDEX_HASH_MAX]` (`KINDEX_HASH_MAX == 1000003`,
//! `KmerIndex.hpp:20,25`), bucketed purely for `std::map` lookup
//! performance via `GetHash(k) = k % 1000003` (`KmerIndex.hpp:28-31`).
//! `Search`/`Insert`/`Remove` on a given code only ever touch that code's
//! single `std::map` key within its one bucket, so a single
//! `HashMap<u64, Vec<IndexInfo>>` keyed directly by the forward code is
//! behaviorally identical without replicating the 1000003-way array.
//!
//! # Not ported (out of scope for this port)
//!
//! `Clear`, `UpdateIndexFromRead`, and `RemoveIndexFromRead`
//! (`KmerIndex.hpp:50-56,133-191`) are not ported here; only
//! `Insert`/`Search`/`Remove`/`BuildIndexFromRead` (the primary
//! build/query/remove surface) are needed by callers so far.

use crate::kmer::KmerCode;
use std::collections::HashMap;

/// Matches T1K's `typedef uint32_t index_t` (`defs.h:15`), the declared type
/// of `_indexInfo::idx`/`_indexInfo::offset` and of the `idx`/`offset`
/// parameters to `Insert`/`Remove`/`BuildIndexFromRead`.
pub type IndexT = u32;

/// A single `(idx, offset)` entry, ported from T1K's `_indexInfo`
/// (`KmerIndex.hpp:12-17`).
///
/// Deliberately has no `strand` field: the C++ struct's `strand` field is
/// commented out (see module docs), so it is never stored by the vendored
/// implementation either.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexInfo {
    pub idx: IndexT,
    pub offset: IndexT,
}

/// K-mer position index, ported from T1K's `KmerIndex` (`KmerIndex.hpp:22-192`).
///
/// See the module-level docs for the forward-vs-canonical keying and
/// bucketing notes.
#[derive(Debug, Clone, Default)]
pub struct KmerIndex {
    /// Equivalent to the union of all `KINDEX_HASH_MAX` (1000003) buckets in
    /// the C++ `index[KINDEX_HASH_MAX]` array (see module docs): keyed
    /// directly by the forward `KmerCode::GetCode()`.
    index: HashMap<u64, Vec<IndexInfo>>,
}

impl KmerIndex {
    /// Constructs an empty `KmerIndex`. Mirrors the `KmerIndex()` constructor
    /// (`KmerIndex.hpp:33-36`), minus the `KINDEX_HASH_MAX`-bucket array
    /// allocation, which is an implementation detail not ported (see module
    /// docs).
    #[must_use]
    pub fn new() -> Self {
        Self { index: HashMap::new() }
    }

    /// Ported from `KmerIndex::Insert` (`KmerIndex.hpp:58-71`).
    ///
    /// If `!kmer_code.is_valid()`, returns without inserting anything
    /// (`KmerIndex.hpp:60-61`). Otherwise appends a new [`IndexInfo { idx,
    /// offset }`](IndexInfo) to the list keyed by `kmer_code.get_code()`
    /// (the *forward* code, not canonical -- see module docs), creating the
    /// list if this is the first entry for that code. Matches
    /// `SimpleVector::PushBack`'s append-only, order-preserving behavior.
    ///
    /// `strand` is accepted (matching the C++ signature) but ignored: see
    /// module docs.
    pub fn insert(&mut self, kmer_code: &KmerCode, idx: IndexT, offset: IndexT, _strand: i32) {
        if !kmer_code.is_valid() {
            return;
        }
        self.index.entry(kmer_code.get_code()).or_default().push(IndexInfo { idx, offset });
    }

    /// Ported from `KmerIndex::Search` (`KmerIndex.hpp:93-105`).
    ///
    /// Returns an empty slice if `!kmer_code.is_valid()` or the code has no
    /// entries (both cases return the C++ `nullHit` sentinel), otherwise the
    /// full entry list for `kmer_code.get_code()` (the *forward* code), in
    /// insertion order.
    #[must_use]
    pub fn search(&self, kmer_code: &KmerCode) -> &[IndexInfo] {
        if !kmer_code.is_valid() {
            return &[];
        }
        self.index.get(&kmer_code.get_code()).map_or(&[], Vec::as_slice)
    }

    /// Ported from `KmerIndex::Remove` (`KmerIndex.hpp:73-91`).
    ///
    /// If `!kmer_code.is_valid()`, returns without doing anything
    /// (`KmerIndex.hpp:75-76`, matching `Search`'s own gate rather than
    /// calling through to it). Otherwise scans the entry list for
    /// `kmer_code.get_code()` for the FIRST entry with matching `idx` AND
    /// `offset`, and removes only that one entry via an order-preserving
    /// shift-left (`SimpleVector::Remove`, `SimpleVector.hpp:224-239`, shifts
    /// every later element down by one rather than swap-removing with the
    /// last element). If no entry matches (or the code has no entries at
    /// all), this is a no-op.
    ///
    /// `strand` is accepted (matching the C++ signature) but ignored: see
    /// module docs.
    pub fn remove(&mut self, kmer_code: &KmerCode, idx: IndexT, offset: IndexT, _strand: i32) {
        if !kmer_code.is_valid() {
            return;
        }
        if let Some(list) = self.index.get_mut(&kmer_code.get_code()) {
            if let Some(pos) = list.iter().position(|e| e.idx == idx && e.offset == offset) {
                list.remove(pos);
            }
        }
    }

    /// Ported from `KmerIndex::BuildIndexFromRead` (`KmerIndex.hpp:107-130`).
    ///
    /// Rolls `kmer_code` (which callers must construct via
    /// [`KmerCode::new`] with the desired k-mer length beforehand) across
    /// `s`, mirroring the C++ loop structure exactly:
    ///
    /// - If `s.len() < kmer_code.get_kmer_length()`, returns immediately
    ///   without touching `kmer_code` (`KmerIndex.hpp:111-112`).
    /// - Otherwise `kmer_code.restart()`s (`KmerIndex.hpp:113`) and fills the
    ///   first `kmer_length - 1` bases with no insertion (`KmerIndex.hpp:
    ///   116-117`) -- the window isn't full yet.
    /// - From the first full window onward (index `kmer_length - 1` through
    ///   `s.len() - 1`), each `append` is followed by a **consecutive-
    ///   duplicate dedup check**: insert only if the window is valid AND
    ///   either this is the *second* full window (`i == kmer_length`,
    ///   unconditional -- `KmerIndex.hpp:121`) or the current code differs
    ///   from the immediately preceding window's code
    ///   (`!kmer_code.is_equal(&prev_kmer_code)`). This drops repeated
    ///   insertions for runs of identical consecutive k-mer codes (e.g. long
    ///   homopolymer stretches, where every shifted window is bit-for-bit
    ///   the same code).
    ///   - Edge case: `prev_kmer_code` starts as a **freshly constructed**
    ///     `KmerCode::new(kmer_length)` (`KmerIndex.hpp:115`), which has
    ///     `code == 0` (same as an all-`A` k-mer). At the very first full
    ///     window (`i == kmer_length - 1`, so the `i == kmer_length`
    ///     unconditional branch does NOT apply), if that first window's code
    ///     is *also* 0 (i.e. it is all-`A`), `is_equal` returns `true`
    ///     against this initial `prev_kmer_code`, so the first window is
    ///     silently dropped -- this port reproduces that exactly by
    ///     constructing `prev_kmer_code` the same way.
    /// - The offset passed to [`KmerIndex::insert`] is `i - kmer_length + 1 +
    ///   shift`, computed as a signed `i64` (to avoid the intermediate
    ///   arithmetic overflowing/panicking in Rust the way C++'s 32-bit `int`
    ///   silently wraps) and then truncated to `i32` before being cast to
    ///   `IndexT` (`u32`) -- mirroring C++'s implicit `int -> uint32_t`
    ///   conversion at the `Insert` call site (`KmerIndex.hpp:123`), which is
    ///   a well-defined two's-complement reinterpretation, not a
    ///   clamp/panic. `id` is converted to `IndexT` the same way. The strand
    ///   argument to `Insert` is always `1` (hardcoded, `KmerIndex.hpp:123`)
    ///   -- and, per [`KmerIndex::insert`], ignored regardless. The
    ///   commented-out reverse-complement insertion (`KmerIndex.hpp:125-126`)
    ///   is dead code in the vendored C++ and is not ported.
    pub fn build_index_from_read(
        &mut self,
        kmer_code: &mut KmerCode,
        s: &[u8],
        id: i32,
        shift: i32,
    ) {
        let len = s.len();
        let kl = kmer_code.get_kmer_length();
        if len < kl {
            return;
        }
        kmer_code.restart();
        let mut prev_kmer_code = KmerCode::new(kl);

        let mut i = 0usize;
        let fill = kl.saturating_sub(1);
        while i < fill {
            kmer_code.append(s[i]);
            i += 1;
        }

        while i < len {
            kmer_code.append(s[i]);
            if kmer_code.is_valid() && (i == kl || !kmer_code.is_equal(&prev_kmer_code)) {
                // `i`/`kl` are window positions bounded by `len` (a real
                // read length), so this never approaches i64::MAX; the
                // `as i64` casts just widen an already-small usize for
                // signed arithmetic (mirroring C++'s signed `int i`/`kl`).
                #[allow(clippy::cast_possible_wrap)]
                let offset_i64 = (i as i64) - (kl as i64) + 1 + i64::from(shift);
                #[allow(clippy::cast_possible_truncation)]
                let offset = offset_i64 as i32;
                // Mirrors C++'s implicit `int -> uint32_t` conversion when
                // passing `id`/`offset` (both signed `int` expressions) to
                // `Insert`'s `index_t` (`uint32_t`) parameters: a
                // well-defined two's-complement reinterpretation, not a
                // clamp -- so intentionally allowing sign loss here.
                #[allow(clippy::cast_sign_loss)]
                let idx = id as IndexT;
                #[allow(clippy::cast_sign_loss)]
                let offset = offset as IndexT;
                self.insert(kmer_code, idx, offset, 1);
            }
            prev_kmer_code = kmer_code.clone();
            i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn code_for(k: usize, bytes: &[u8]) -> KmerCode {
        let mut kc = KmerCode::new(k);
        for &c in bytes {
            kc.append(c);
        }
        kc
    }

    #[test]
    fn insert_and_search_preserve_order() {
        let mut idx = KmerIndex::new();
        let kc = code_for(4, b"ACGT");
        idx.insert(&kc, 1, 0, 1);
        idx.insert(&kc, 2, 10, 1);
        idx.insert(&kc, 3, 20, 1);
        assert_eq!(
            idx.search(&kc),
            &[
                IndexInfo { idx: 1, offset: 0 },
                IndexInfo { idx: 2, offset: 10 },
                IndexInfo { idx: 3, offset: 20 },
            ]
        );
    }

    #[test]
    fn insert_ignores_strand_parameter() {
        let mut idx = KmerIndex::new();
        let kc = code_for(4, b"ACGT");
        idx.insert(&kc, 5, 7, 1);
        idx.insert(&kc, 5, 7, -1);
        idx.insert(&kc, 5, 7, 999);
        // Three inserts of the identical (idx, offset), regardless of
        // strand, produce three identical entries (no dedup, no strand
        // stored/observed).
        assert_eq!(idx.search(&kc), &[IndexInfo { idx: 5, offset: 7 }; 3]);
    }

    #[test]
    fn forward_code_not_canonical_keys_distinctly() {
        // "AAAC" (k=4) and its reverse complement "GTTT" have the same
        // canonical code but different forward codes.
        let mut idx = KmerIndex::new();
        let fwd = code_for(4, b"AAAC");
        let rc = code_for(4, b"GTTT");
        assert_eq!(fwd.get_canonical_kmer_code(), rc.get_canonical_kmer_code());
        assert_ne!(fwd.get_code(), rc.get_code());

        idx.insert(&fwd, 1, 0, 1);
        idx.insert(&rc, 2, 0, 1);

        assert_eq!(idx.search(&fwd), &[IndexInfo { idx: 1, offset: 0 }]);
        assert_eq!(idx.search(&rc), &[IndexInfo { idx: 2, offset: 0 }]);
    }

    #[test]
    fn invalid_kmer_is_dropped_by_insert_and_search() {
        let mut idx = KmerIndex::new();
        let invalid = code_for(4, b"ACNT");
        assert!(!invalid.is_valid());
        idx.insert(&invalid, 1, 0, 1);
        assert!(idx.search(&invalid).is_empty());
    }

    #[test]
    fn remove_first_match_preserves_order_of_rest() {
        let mut idx = KmerIndex::new();
        let kc = code_for(4, b"ACGT");
        idx.insert(&kc, 1, 0, 1);
        idx.insert(&kc, 2, 0, 1); // duplicate idx/offset-distinct entry
        idx.insert(&kc, 1, 0, 1); // duplicate of the first entry
        idx.insert(&kc, 3, 9, 1);

        idx.remove(&kc, 1, 0, 1);
        // Only the FIRST (1, 0) is removed; the second (1, 0) later in the
        // list remains, and relative order of the rest is preserved.
        assert_eq!(
            idx.search(&kc),
            &[
                IndexInfo { idx: 2, offset: 0 },
                IndexInfo { idx: 1, offset: 0 },
                IndexInfo { idx: 3, offset: 9 },
            ]
        );
    }

    #[test]
    fn remove_nonexistent_entry_is_noop() {
        let mut idx = KmerIndex::new();
        let kc = code_for(4, b"ACGT");
        idx.insert(&kc, 1, 0, 1);
        idx.remove(&kc, 99, 99, 1);
        assert_eq!(idx.search(&kc), &[IndexInfo { idx: 1, offset: 0 }]);
    }

    #[test]
    fn remove_on_invalid_kmer_is_noop_not_panic() {
        let mut idx = KmerIndex::new();
        let invalid = code_for(4, b"ACNT");
        idx.remove(&invalid, 1, 0, 1); // must not panic
        assert!(idx.search(&invalid).is_empty());
    }

    #[test]
    fn build_index_from_read_hand_computed_offsets() {
        // k=3, read="ACGTAC" (len=6), shift=0, id=42.
        // Windows (0-based start): 0:"ACG" 1:"CGT" 2:"GTA" 3:"TAC".
        // First full window completes at i=2 (0-based) == kl-1=2. At this
        // point prev_kmer_code is the initial KmerCode::new(3) (code=0,
        // valid). "ACG" != all-A, so its code != 0, so !is_equal(prev) is
        // true -> inserted at offset = i - kl + 1 = 2-3+1 = 0.
        // i=3 ("CGT"): i == kl (3==3) -> unconditionally inserted at
        // offset = 3-3+1 = 1, regardless of equality to prev.
        // i=4 ("GTA"): i != kl, and code("GTA") != code("CGT") -> inserted
        // at offset = 4-3+1 = 2.
        // i=5 ("TAC"): i != kl, code("TAC") != code("GTA") -> inserted at
        // offset = 5-3+1 = 3.
        // So all four windows are inserted here (no duplicate consecutive
        // codes in this read), each exactly once, at offsets 0..=3.
        let mut idx = KmerIndex::new();
        let mut kc = KmerCode::new(3);
        idx.build_index_from_read(&mut kc, b"ACGTAC", 42, 0);

        for (window, expected_offset) in
            [(&b"ACG"[..], 0u32), (&b"CGT"[..], 1), (&b"GTA"[..], 2), (&b"TAC"[..], 3)]
        {
            let q = code_for(3, window);
            assert_eq!(
                idx.search(&q),
                &[IndexInfo { idx: 42, offset: expected_offset }],
                "window {window:?} (as utf8: {:?})",
                String::from_utf8_lossy(window)
            );
        }
    }

    #[test]
    fn build_index_from_read_dedups_homopolymer_run_but_drops_first_window() {
        // k=3, read="AAAAAA" (len=6, all-A homopolymer), shift=0, id=1.
        // Every full window ("AAA") has forward code 0 (A=0 in 2-bit
        // packing), so every window is bit-for-bit identical.
        //
        // i=2 (kl-1=2, first full window "AAA" at offset 0): prev_kmer_code
        // is the initial KmerCode::new(3), which ALSO has code=0. So
        // is_equal(prev) is true -> the `i == kl` branch doesn't apply
        // (2 != 3) -> NOT inserted. This is the documented edge case.
        // i=3 (i == kl == 3, second full window "AAA" at offset 1):
        // unconditionally inserted (the `i == kl` branch short-circuits the
        // equality check) -> ONE entry recorded, at offset 1.
        // i=4, i=5: code equals prev_kmer_code (still all-A) -> deduped,
        // not inserted.
        //
        // Net result: exactly ONE entry for the "AAA" code, at offset 1 (not
        // offset 0), despite there being 4 full windows in the read.
        let mut idx = KmerIndex::new();
        let mut kc = KmerCode::new(3);
        idx.build_index_from_read(&mut kc, b"AAAAAA", 1, 0);

        let q = code_for(3, b"AAA");
        assert_eq!(idx.search(&q), &[IndexInfo { idx: 1, offset: 1 }]);
    }

    #[test]
    fn build_index_from_read_returns_early_for_short_read() {
        let mut idx = KmerIndex::new();
        let mut kc = KmerCode::new(5);
        idx.build_index_from_read(&mut kc, b"ACGT", 1, 0); // len=4 < k=5
        let q = code_for(5, b"ACGTA");
        assert!(idx.search(&q).is_empty());
    }

    #[test]
    fn build_index_from_read_skips_invalid_windows() {
        // k=3, read="ACNGAC": windows are "ACN"(invalid), "CNG"(invalid),
        // "NGA"(invalid), "GAC"(valid, N aged out). Only "GAC" is inserted.
        let mut idx = KmerIndex::new();
        let mut kc = KmerCode::new(3);
        idx.build_index_from_read(&mut kc, b"ACNGAC", 7, 0);

        let gac = code_for(3, b"GAC");
        assert_eq!(idx.search(&gac), &[IndexInfo { idx: 7, offset: 3 }]);
        // No entries at all for any of the invalid windows' would-be codes
        // beyond what's asserted above; spot check one directly-computed
        // invalid query returns empty too (mirrors Search's IsValid gate).
        let acn = code_for(3, b"ACN");
        assert!(idx.search(&acn).is_empty());
    }

    #[test]
    fn build_index_from_read_applies_shift_to_offset() {
        let mut idx = KmerIndex::new();
        let mut kc = KmerCode::new(3);
        idx.build_index_from_read(&mut kc, b"ACGTAC", 1, 100);
        // Same windows as the hand-computed test above, but every offset is
        // shifted by +100.
        let q = code_for(3, b"ACG");
        assert_eq!(idx.search(&q), &[IndexInfo { idx: 1, offset: 100 }]);
    }

    #[test]
    fn build_index_from_read_negative_shift_wraps_like_cpp_implicit_conversion() {
        // shift=-10 on the first window (offset expression 0 + (-10) = -10)
        // must reinterpret as a u32 via two's-complement wraparound, exactly
        // like C++'s implicit `int -> uint32_t` conversion when passing a
        // negative int to an `index_t` parameter.
        let mut idx = KmerIndex::new();
        let mut kc = KmerCode::new(3);
        idx.build_index_from_read(&mut kc, b"ACGTAC", 1, -10);
        let q = code_for(3, b"ACG");
        // Two's-complement reinterpretation, matching the production code's
        // own `as IndexT` cast -- not a numeric truncation bug.
        #[allow(clippy::cast_sign_loss)]
        let expected_offset = (-10i32) as u32;
        assert_eq!(idx.search(&q), &[IndexInfo { idx: 1, offset: expected_offset }]);
    }
}
