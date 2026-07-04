//! Reference-k-mer read-candidate filter, ported from a slice of T1K's
//! `SeqSet` class (`SeqSet.hpp`) plus the free functions
//! `IsLowComplexity`/`IsGoodCandidate` (`FastqExtractor.cpp:89-118`,
//! byte-identical to `BamExtractor.cpp:144-166,214`).
//!
//! This is the *read-candidate-filtering* slice only -- the part of `SeqSet`
//! that `FastqExtractor`/`BamExtractor` use to decide whether a read is
//! worth keeping for downstream genotyping: reference load + index build
//! (`SeqSet::InputRefFa`), k-mer hit collection (`SeqSet::GetHitsFromRead`),
//! and BOTH gates of `SeqSet::HasHitInSet` -- the bucket-count gate AND the
//! `GetOverlapsFromHits`-based (LIS hit-chaining) mismatch-threshold
//! confirmation (ported in [`crate::overlap`]). It does **not** port
//! `AlignAlgo` (banded Smith-Waterman alignment) -- `GetOverlapsFromHits`,
//! as called from `HasHitInSet`, never invokes it (see [`crate::overlap`]'s
//! module docs); real alignment-based scoring is a separate, later-phase
//! concern unrelated to this gate.
//!
//! # `kmerLength` is a caller-supplied parameter here, not a fixed constant
//!
//! Stock T1K does not use a single hardcoded default k-mer length. In
//! `FastqExtractor.cpp:main` (`FastqExtractor.cpp:272-418`):
//! 1. `SeqSet refSet(9)` is constructed with a literal initial `kmerLength = 9`.
//! 2. The reference is loaded at that k-mer length (`refSet.InputRefFa(...)`,
//!    `FastqExtractor.cpp:302`).
//! 3. `SeqSet::InferKmerLength()` (`SeqSet.hpp:2830-2845`) computes a
//!    *data-dependent* k-mer length from the total loaded reference length
//!    (`totalLength`): `kmerLength = floor(log4(totalLength)) + 2` (the loop
//!    counts base-4 "digits" of `totalLength`, then adds one more).
//! 4. If that inferred value is *strictly greater* than the initial `9`, the
//!    `SeqSet` is rebuilt at the new k-mer length via `UpdateKmerLength`
//!    (`FastqExtractor.cpp:411-418`, `SeqSet.hpp:2847-2857`); otherwise it
//!    stays at `9`.
//!
//! [`RefKmerFilter::from_reference_fasta`] itself does not replicate this
//! two-stage default-then-infer-then-maybe-rebuild dance -- it takes
//! `kmer_length` as an explicit parameter, matching a single plain
//! `SeqSet(kmerLength)` construction followed by one `InputRefFa` call
//! (`SeqSet.hpp:760-772,872-904`, no `InferKmerLength`/`UpdateKmerLength`
//! call during construction itself). [`RefKmerFilter::infer_kmer_length`]
//! and [`RefKmerFilter::update_kmer_length`] (added for Task 3.2) now DO
//! port `SeqSet::InferKmerLength`/`SeqSet::UpdateKmerLength` faithfully; a
//! caller that wants stock's exact default-selection dance (e.g.
//! [`crate::extract`]) calls `from_reference_fasta` at an initial k, then
//! `infer_kmer_length`, then conditionally `update_kmer_length` -- see
//! [`crate::extract`]'s module docs for that full sequence. For
//! `fixtures/refbuild/golden/kir_rna_seq.fa` specifically (the fixture used
//! by this module's own differential tests, distinct from
//! `fixtures/example/ref/kir_rna_seq.fa` used by [`crate::extract`]'s
//! tests), `InferKmerLength()` evaluates to `8` for that fixture's total
//! reference length (8781 bases) -- NOT greater than the initial `9` -- so
//! stock T1K would actually use `kmerLength = 9` unmodified for this exact
//! fixture; the differential tests in THIS module use `kmer_length = 9`
//! directly (not via the infer/update dance) to match that real-world value
//! precisely.
//!
//! # `HasHitInSet`'s two gates, both now ported
//!
//! Stock `SeqSet::HasHitInSet` (`SeqSet.hpp:1915-1990`) does two things in
//! sequence:
//! 1. Collect hits from both strands (`GetHitsFromRead`), bucket them by
//!    `(strand-tag, seqIdx)`, and find the single largest bucket. If
//!    `kmerLength * <largest bucket size> < hitLenRequired`, return `false`
//!    immediately (`SeqSet.hpp:1929-1964`).
//! 2. Otherwise, run `GetOverlapsFromHits` (LIS-based colinear hit chaining;
//!    NOT `AlignAlgo`-based -- see [`crate::overlap`]'s module docs) on the
//!    winning bucket, and only return `true` if at least one resulting
//!    overlap's implied mismatch count is within `int(len *
//!    (1 - refSeqSimilarity)) * kmerLength` (`SeqSet.hpp:1966-1983`).
//!
//! [`RefKmerFilter::has_hit_in_set`] implements BOTH gates: step 1 (the
//! bucket-count gate) followed by step 2 (calling
//! [`crate::overlap::get_overlaps_from_hits`] on the winning bucket's hits,
//! reconstructed in original insertion order -- see that function's doc
//! comment for why a stable filter of `scratch.hits` is exactly equivalent
//! to stock's per-hit `buckets[tag][idx].PushBack(...)` construction). This
//! makes `has_hit_in_set` byte/value-identical to stock's `HasHitInSet` (not
//! merely a superset), confirmed by an adversarial differential
//! (`crates/unum-core/tests/golden_refkmerfilter.rs`) against the real C++
//! oracle across curated AND arbitrary/random/noncolinear/tie-inducing reads.
//!
//! # Reused Phase-2 primitives
//!
//! Reuses [`crate::kmer::KmerCode`] (rolling k-mer encoder) and
//! [`crate::kmer_index::KmerIndex`] (forward-code-keyed position index)
//! unmodified. `KmerIndex::search`'s forward-not-canonical keying and
//! `KmerIndex::build_index_from_read`'s consecutive-duplicate dedup are
//! exactly what `SeqSet::InputRefFa`'s `seqIndex.BuildIndexFromRead(...)`
//! call needs (`SeqSet.hpp:902`); see `kmer_index`'s module docs for that
//! dedup's exact semantics (already covered by Phase-2's differential
//! tests, not re-derived here). Also reuses [`crate::overlap`]'s
//! `get_overlaps_from_hits`/LIS port for `has_hit_in_set`'s gate 2.

use crate::kmer::KmerCode;
use crate::kmer_index::{IndexInfo, IndexT, for_each_kmer};
use crate::overlap::{self, OverlapHit};
use anyhow::{Context, Result};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

/// `SeqSet`'s default `hitLenRequired` (`SeqSet.hpp:764`), used by
/// [`RefKmerFilter::from_reference_fasta`] since it mirrors a plain
/// `SeqSet(kmerLength)` construction with no `SetHitLenRequired` call (see
/// module docs: `_load_ref`/`from_reference_fasta` do not replicate
/// `FastqExtractor.cpp`'s dynamic `SetHitLenRequired` call, which requires
/// sampling actual read lengths from a FASTQ file upstream of reference
/// loading).
const DEFAULT_HIT_LEN_REQUIRED: i32 = 31;

/// `SeqSet`'s default `radius` (`SeqSet.hpp:763`), used by
/// [`RefKmerFilter::has_hit_in_set`]'s gate 2 (`overlap::get_overlaps_from_hits`'s
/// `radius` parameter). Fixed at the `SeqSet(kmerLength)` constructor
/// default for the same reason as [`DEFAULT_HIT_LEN_REQUIRED`]: no
/// `SetRadius` call exists on the plain-construction path this port models.
const DEFAULT_RADIUS: i32 = 10;

/// `SeqSet`'s default `refSeqSimilarity` (`SeqSet.hpp:768`), used by
/// [`RefKmerFilter::has_hit_in_set`]'s gate 2 mismatch-threshold computation
/// (`SeqSet.hpp:1973`). Fixed at the `SeqSet(kmerLength)` constructor
/// default; see [`DEFAULT_HIT_LEN_REQUIRED`]'s doc comment for why this port
/// does not replicate any dynamic `Set*` call.
const DEFAULT_REF_SEQ_SIMILARITY: f64 = 0.8;

/// `SeqSet::isLongSeqSet`'s default (`SeqSet.hpp:765`). Always `false` on the
/// plain-construction path this port models (no public setter reachable from
/// `SeqSet(kmerLength)`; see `get_hits_from_read`'s doc comment for the same
/// point made about `downSample`). Passed straight through to
/// `overlap::get_overlaps_from_hits`'s `is_long_seq_set` parameter.
const IS_LONG_SEQ_SET: bool = false;

/// `SeqSet`'s default `novelSeqSimilarity` (`SeqSet.hpp:767`), used by
/// [`RefKmerFilter::get_overlaps_from_read`]'s final similarity filter for
/// non-`isRef` overlaps (`SeqSet.hpp:1898`). Fixed at the `SeqSet(kmerLength)`
/// constructor default; see [`DEFAULT_HIT_LEN_REQUIRED`]'s doc comment for
/// why this port does not replicate any dynamic `Set*` call. Every sequence
/// loaded via [`RefKmerFilter::from_reference_fasta`] is `isRef == true` (see
/// module docs), so this constant is currently unexercised by any reachable
/// overlap in this codebase -- ported anyway, per this port's "ported
/// generally, not overfit" discipline (see [`crate::overlap`]'s module docs
/// for the same principle applied to `GetOverlapsFromHits`'s `filter == 1`
/// branches).
const DEFAULT_NOVEL_SEQ_SIMILARITY: f64 = 0.9;

/// A single k-mer hit collected by [`RefKmerFilter::get_hits_from_read`],
/// mirroring the FFI-relevant fields of T1K's `struct _hit`
/// (`SeqSet.hpp:66-87`), including `repeats` (`_hit::repeats`,
/// `SeqSet.hpp:72`), populated exactly as stock does (`SeqSet.hpp:
/// 1122,1143,1189,1210`: `repeats = size` of the `seqIndex.Search(...)`
/// result, since `has_hit_in_set`'s `GetHitsFromRead` call always has
/// `barcode == -1` and `puse == NULL`, `SeqSet.hpp:1923`) so
/// [`overlap::get_overlaps_from_hits`]'s `filter == 1` branches (not
/// reachable from `has_hit_in_set`, which always passes `filter == 0`, but
/// ported generally -- see that module's docs) have a correct `repeats` to
/// read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Hit {
    /// The reference sequence index this hit maps to (`_indexInfo::idx`).
    pub idx: u32,
    /// The offset within that reference sequence (`_indexInfo::offset`).
    pub offset: u32,
    /// The offset within the (possibly reverse-complemented) read where this
    /// k-mer window starts (`_hit::readOffset`).
    pub read_offset: i32,
    /// `1` for a forward-strand hit, `-1` for a reverse-complement-strand
    /// hit (`_hit::strand`).
    pub strand: i8,
    /// How many times this hit's k-mer occurs across the whole reference
    /// index (`_hit::repeats`) -- see this struct's doc comment.
    pub repeats: i32,
}

impl Hit {
    /// Converts to [`overlap::OverlapHit`], the standalone hit type
    /// [`overlap::get_overlaps_from_hits`] operates on (see that module's
    /// docs for why it doesn't reuse this crate's `Hit` directly).
    fn to_overlap_hit(self) -> OverlapHit {
        OverlapHit {
            idx: self.idx,
            offset: self.offset,
            read_offset: self.read_offset,
            strand: self.strand,
            repeats: self.repeats,
        }
    }
}

/// Ported from `_hit::operator<` (`SeqSet.hpp:74-86`): ascending by
/// `strand`, then `indexHit.idx`, then `readOffset`, then `indexHit.offset`.
/// This is the comparator `SortHits`'s `std::sort` fallback branch used
/// (`SeqSet.hpp:1589`) before it was replaced by an equivalent radix sort on
/// [`pack_hit_key`]; see [`sort_hits`]'s doc comment for the full
/// bucket-sort-vs-radix dispatch this feeds into. It is retained as the
/// canonical strict-total-order specification that [`pack_hit_key`]'s
/// monotonicity is verified against (`pack_hit_key_matches_hit_less_than_order`);
/// production sorting no longer calls it, hence `#[allow(dead_code)]`.
#[allow(dead_code)]
#[must_use]
fn hit_less_than(a: &Hit, b: &Hit) -> bool {
    if a.strand != b.strand {
        return a.strand < b.strand;
    }
    if a.idx != b.idx {
        return a.idx < b.idx;
    }
    if a.read_offset != b.read_offset {
        return a.read_offset < b.read_offset;
    }
    a.offset < b.offset
}

/// Packs a [`Hit`]'s four ordered fields (`strand`, `idx`, `read_offset`,
/// `offset` -- exactly the fields [`hit_less_than`] compares, in that
/// priority) into a single `u128` whose unsigned ordering is monotonic in
/// [`hit_less_than`]: for any two hits `a`, `b`, `pack(a) < pack(b)` iff
/// `hit_less_than(a, b)`, and `pack(a) == pack(b)` iff `a` and `b` tie in all
/// four fields.
///
/// The packing is most-significant-field-first (`strand`, then `idx`, then
/// `read_offset`, then `offset`), each field occupying a full-width lane so
/// no reference size, read length, or offset value can overflow into a
/// neighbouring field:
///
/// * `strand` (`i8`, only ever `-1` or `+1`): sign-flipped via `^ 0x80` so
///   the two's-complement bit pattern orders `-1 < +1` as an unsigned byte
///   (`0x7f < 0x81`), matching `hit_less_than`'s `a.strand < b.strand`.
/// * `idx` (`u32`): used directly -- unsigned order already matches.
/// * `read_offset` (`i32`): sign-flipped via `^ 0x8000_0000` so negative
///   read offsets sort below non-negative ones as unsigned `u32`s.
/// * `offset` (`u32`): used directly.
///
/// Because [`hit_less_than`] is a strict total order (two hits tie only when
/// bit-identical in these four fields -- and each hit is a unique k-mer
/// match, so real hit lists never actually tie), a stable radix sort keyed on
/// this `u128` produces the SAME unique ordering as `sort_unstable_by(
/// hit_less_than)` -- byte-identical to the C++ oracle. Verified by the
/// `pack_hit_key_matches_hit_less_than_order` unit test.
#[must_use]
#[inline]
fn pack_hit_key(h: &Hit) -> u128 {
    #[allow(clippy::cast_sign_loss)]
    let strand_flip = (h.strand as u8) ^ 0x80;
    #[allow(clippy::cast_sign_loss)]
    let read_offset_flip = (h.read_offset as u32) ^ 0x8000_0000;
    (u128::from(strand_flip) << 96)
        | (u128::from(h.idx) << 64)
        | (u128::from(read_offset_flip) << 32)
        | u128::from(h.offset)
}

/// Ported from `SeqSet::SortHits` (`SeqSet.hpp:1558-1590`): reorders `hits`
/// in place, either via a `(strand-tag, seqIdx)` bucket sort (STABLE within
/// each bucket, since it's built by a single forward `PushBack` pass per
/// hit) or a plain `std::sort` using [`hit_less_than`] as the comparator.
///
/// The bucket-sort path is taken only when `hits.len() > 2 * seq_count &&
/// already_read_order` (mirrors `SeqSet.hpp:1561`); `GetOverlapsFromRead`
/// always passes `already_read_order = true` (`SeqSet.hpp:1605`), so this
/// port's only caller always has that half of the condition satisfied --
/// the `hits.len() > 2 * seq_count` half is genuinely data-dependent and
/// reproduced as-is.
///
/// # Bucket sort is a STABLE reordering by `(tag, idx)`, not a full sort
///
/// Unlike [`hit_less_than`] (which also orders by `readOffset`/`offset`
/// within a tied `(strand, idx)` group), the bucket-sort path does NOT sort
/// within each `(tag, idx)` bucket at all -- it simply partitions `hits`
/// into `2 * seqCnt` buckets (by ascending `tag` then ascending `idx`) and
/// concatenates them back together in bucket order, preserving each hit's
/// ORIGINAL relative order within its bucket (`PushBack` in encounter
/// order, `SeqSet.hpp:1570-1574`). This port reproduces that exact
/// semantics via a stable sort keyed only on `(tag, idx)` (Rust's
/// `sort_by_key`/`sort_by` are documented stable), NOT [`hit_less_than`]
/// (which would additionally reorder within a bucket by `readOffset`/
/// `offset` -- a different, and for this branch WRONG, result).
pub(crate) fn sort_hits(hits: &mut [Hit], already_read_order: bool, seq_count: usize) {
    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
    let two_seq_count = 2 * seq_count as i64;
    if i64::try_from(hits.len()).unwrap_or(i64::MAX) > two_seq_count && already_read_order {
        // Bucket sort: stable partition by (tag, idx), tag = (strand==1) as
        // 0/1, buckets visited tag=0 (idx ascending) THEN tag=1 (idx
        // ascending) -- SeqSet.hpp:1577-1583's `for k in 0..=1 { for i in
        // 0..seqCnt { ... } }` nesting.
        hits.sort_by_key(|h| (i8::from(h.strand == 1), h.idx));
    } else {
        // `hit_less_than` is a strict total order over every field of `Hit`
        // that participates in equality (strand, idx, read_offset, offset --
        // `repeats` is not compared), so two hits only tie here if they are
        // bit-identical in those four fields; the sorted order is therefore
        // UNIQUELY determined and independent of the sort's stability. We
        // exploit that to replace the O(n log n) comparison sort with an
        // O(n) LSD radix sort keyed on `pack_hit_key`, whose unsigned u128
        // order is monotonic in `hit_less_than` (see `pack_hit_key`'s doc and
        // the `pack_hit_key_matches_hit_less_than_order` test) -- this yields
        // the identical unique ordering, byte-identical to the C++ oracle,
        // while being the dominant cost of `GetOverlapsFromRead`. `radsort`'s
        // internal `unsafe` is the crate's own; `#![forbid(unsafe_code)]` on
        // this crate is unaffected.
        radsort::sort_by_key(hits, pack_hit_key);
    }
}

/// Reusable scratch buffers for [`RefKmerFilter`]'s per-read query path
/// (`has_hit_in_set`/`is_good_candidate`), so a caller processing many reads
/// in a loop (a future `extract`/`genotype` command) does not re-allocate a
/// reverse-complement buffer, hit list, and bucket map on every call.
#[derive(Debug, Default, Clone)]
pub struct Scratch {
    /// Reused buffer for the read's reverse complement (mirrors T1K's
    /// caller-owned `rcRead`/`buffer` scratch parameter to `HasHitInSet`/
    /// `GetHitsFromRead`, which is purely an output parameter in the C++ --
    /// see `get_hits_from_read`'s doc comment).
    rc_buf: Vec<u8>,
    /// Reused hit list, cleared and refilled by every `get_hits_from_read` call.
    hits: Vec<Hit>,
    /// Reused `(tag, seqIdx) -> count` bucket map for `has_hit_in_set`'s
    /// add02ca touched-buckets bucket sort (see that function's doc comment).
    buckets: HashMap<(i8, u32), u32>,
    /// Per-thread memoization cache for the DP alignment path
    /// (`global_alignment`). `assign_reads_parallel`'s `map_init` gives one
    /// `Scratch` -- hence one `DpCache` -- per rayon worker, which is the
    /// lock-free thread-local equivalent of the fork's `thread_local`
    /// `unordered_map` (folds fork commit `a35ed72`). See [`DpCache`].
    ///
    /// [`DpCache`]: crate::align_algo::DpCache
    pub dp_cache: crate::align_algo::DpCache,
}

/// Flat, search-only forward-code-keyed k-mer position index.
///
/// Behaviorally identical to [`crate::kmer_index::KmerIndex`] for the
/// `search` surface (the only surface `RefKmerFilter` needs) -- it returns a
/// byte-identical `&[IndexInfo]` slice for any given forward code -- but
/// stores all entries in ONE contiguous `positions` blob (each code mapping
/// to a `(start, len)` window into it) rather than a heap `Vec` per code.
/// That cuts allocator overhead and index RSS and keeps lookups
/// cache-friendly. Keying uses [`FxHashMap`] (rustc-hash) instead of std's
/// SipHash for faster build + lookup; rustc-hash is pure-safe, so
/// `#![forbid(unsafe_code)]` still holds.
///
/// # Byte-identity with `KmerIndex`
///
/// `KmerIndex::build_index_from_read` pushes each code's entries in (`id`
/// ascending, then `offset` ascending within an `id`) order, because
/// `from_reference_fasta`/`update_kmer_length` process sequence `id` 0 fully
/// before `id` 1 and roll offsets left->right within a sequence. This flat
/// build collects `(code, idx, offset)` tuples from every sequence via the
/// SAME shared [`for_each_kmer`] emitter, then `par_sort_unstable`s them
/// lexicographically. Since a given forward code occurs at most once per
/// `(seq, offset)`, the `(code, idx, offset)` tuples are all-distinct within
/// a code, so the sort order is unambiguous and reproduces the sequential
/// push order EXACTLY. Grouping consecutive-equal-`code` runs then yields,
/// for each code, a `positions` slice byte-identical to the old per-code
/// `Vec`.
#[derive(Debug, Default)]
struct FlatKmerIndex {
    /// Forward `KmerCode::get_code()` -> `(start, len)` window into
    /// [`positions`](Self::positions).
    map: FxHashMap<u64, (u32, u32)>,
    /// All `IndexInfo` entries across every code, contiguous and grouped by
    /// code, each group in (idx, offset)-ascending order (see type docs).
    positions: Vec<IndexInfo>,
}

impl FlatKmerIndex {
    /// Builds the flat index over `seqs` at k-mer length `kmer_length`, in
    /// parallel across the current rayon thread pool. Mirrors the result of
    /// running [`KmerIndex::build_index_from_read`] over `seqs` in 0-based
    /// load order (`id` = index), byte-for-byte (see type docs for the
    /// ordering proof).
    fn build(seqs: &[Vec<u8>], kmer_length: usize) -> Self {
        // 1. Emit every (code, idx, offset) tuple, in parallel per sequence,
        //    via the SAME shared emitter the sequential build uses.
        let mut tuples: Vec<(u64, IndexT, IndexT)> = seqs
            .par_iter()
            .enumerate()
            .flat_map_iter(|(seq_idx, seq)| {
                // `seqs.len()` is bounded by an `i32` at load time
                // (`from_reference_fasta`'s `i32::try_from`), so this cannot
                // truncate; matches the `id: i32` the sequential build uses.
                #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                let id = seq_idx as i32;
                let mut local: Vec<(u64, IndexT, IndexT)> = Vec::new();
                for_each_kmer(seq, id, kmer_length, 0, &mut |code, idx, offset| {
                    local.push((code, idx, offset));
                });
                local.into_iter()
            })
            .collect();

        // 2. Sort lexicographically by (code, idx, offset). This is exactly
        //    the sequential push order per code (see type docs). Unstable is
        //    safe because the tuples are all-distinct within a code.
        tuples.par_sort_unstable();

        // 3. Group consecutive-equal-`code` runs into `(start, len)` windows
        //    and flatten into the contiguous `positions` blob.
        let mut map: FxHashMap<u64, (u32, u32)> = FxHashMap::default();
        let mut positions: Vec<IndexInfo> = Vec::with_capacity(tuples.len());
        let mut run_start = 0usize;
        while run_start < tuples.len() {
            let code = tuples[run_start].0;
            let mut run_end = run_start + 1;
            while run_end < tuples.len() && tuples[run_end].0 == code {
                run_end += 1;
            }
            // Positions indices fit in u32: total entries are bounded by the
            // reference size (well under u32::MAX for any real T1K reference).
            #[allow(clippy::cast_possible_truncation)]
            let start = positions.len() as u32;
            #[allow(clippy::cast_possible_truncation)]
            let len = (run_end - run_start) as u32;
            map.insert(code, (start, len));
            for &(_, idx, offset) in &tuples[run_start..run_end] {
                positions.push(IndexInfo { idx, offset });
            }
            run_start = run_end;
        }

        Self { map, positions }
    }

    /// Ported search surface, byte-identical to
    /// [`KmerIndex::search`](crate::kmer_index::KmerIndex::search): returns an
    /// empty slice if `!kmer_code.is_valid()` or the forward code has no
    /// entries, otherwise the code's full entry slice in insertion order.
    #[must_use]
    fn search(&self, kmer_code: &KmerCode) -> &[IndexInfo] {
        if !kmer_code.is_valid() {
            return &[];
        }
        match self.map.get(&kmer_code.get_code()) {
            Some(&(start, len)) => {
                let start = start as usize;
                let end = start + len as usize;
                &self.positions[start..end]
            }
            None => &[],
        }
    }
}

/// Reference-k-mer read-candidate filter: the pure-Rust port of the
/// `SeqSet` slice `FastqExtractor`/`BamExtractor` use to decide whether a
/// read is a genotyping candidate. See the module docs for exact scope and
/// the deliberate `GetOverlapsFromHits`/`AlignAlgo` scope cut.
pub struct RefKmerFilter {
    /// Forward-code-keyed k-mer position index over every loaded reference
    /// sequence, built (in parallel) via the shared [`for_each_kmer`] emitter
    /// -- byte-identical to the sequential `KmerIndex::build_index_from_read`
    /// path (one emission per sequence, `id` = 0-based load order) -- mirrors
    /// `SeqSet::seqIndex` after `InputRefFa` (`SeqSet.hpp:220,872-904`). Stored
    /// flat (contiguous positions blob) rather than a heap `Vec` per code; see
    /// [`FlatKmerIndex`].
    seq_index: FlatKmerIndex,
    /// The k-mer length every sequence was indexed at (`SeqSet::kmerLength`,
    /// `SeqSet.hpp:221`). See the module docs for why this is an explicit
    /// caller-supplied parameter rather than a fixed default.
    kmer_length: usize,
    /// Number of reference sequences loaded (`SeqSet::seqs.size()` after
    /// `InputRefFa`; mirrors `SeqSet::Size()`, `SeqSet.hpp:795-798`).
    seq_count: usize,
    /// 0-based-load-order sequence names, in load order (index == the `idx`
    /// used to key `seq_index`). Not required by the ported slice itself,
    /// but cheap to keep and useful for a future caller that wants to report
    /// which reference sequence a candidate read matched.
    seq_names: Vec<String>,
    /// 0-based-load-order raw reference sequences (`SeqSet::seqs[i].consensus`
    /// equivalent), kept alongside `seq_names` so
    /// [`RefKmerFilter::update_kmer_length`] can rebuild `seq_index` at a new
    /// k-mer length without re-reading the source FASTA (mirrors
    /// `SeqSet::UpdateKmerLength` re-iterating `seqs[i].consensus`,
    /// `SeqSet.hpp:2847-2857`).
    seqs: Vec<Vec<u8>>,
    /// `SeqSet::hitLenRequired` (`SeqSet.hpp:223`), used by
    /// [`RefKmerFilter::has_hit_in_set`]'s bucket-count gate. Fixed at the
    /// `SeqSet(kmerLength)` constructor default (`SeqSet.hpp:764`); see
    /// [`DEFAULT_HIT_LEN_REQUIRED`]'s doc comment for why this is not
    /// dynamically computed the way `FastqExtractor.cpp`'s `main` does.
    hit_len_required: i32,
    /// `SeqSet::refSeqSimilarity` (`SeqSet.hpp:231`), used by
    /// [`RefKmerFilter::has_hit_in_set`]'s gate 2 mismatch-threshold
    /// computation. Fixed at the `SeqSet(kmerLength)` constructor default;
    /// see [`DEFAULT_REF_SEQ_SIMILARITY`]'s doc comment.
    ref_seq_similarity: f64,
}

impl RefKmerFilter {
    /// Loads a T1K-style reference FASTA (e.g. a `*_seq.fa` produced by
    /// `prepare-ref`/`t1k-build.pl`) and builds a k-mer index over every
    /// sequence, mirroring `SeqSet(kmerLength)` (`SeqSet.hpp:760-772`)
    /// followed by `InputRefFa(filename)` (`SeqSet.hpp:872-904`): each FASTA
    /// record is indexed with `id` = its 0-based position in the file
    /// (matching stock's `id = seqs.size()` at insertion time).
    ///
    /// Sequences must already be uppercase `A`/`C`/`G`/`T`/`N` (T1K assumes
    /// this; see [`crate::kmer::canonical_kmer`]'s doc comment for the same
    /// caveat on the reused `KmerCode`).
    ///
    /// # Errors
    ///
    /// Returns an error if `path` cannot be read, or if the file contains
    /// more sequences than fit in an `i32` (mirrors `KmerIndex`'s `id: i32`
    /// parameter, itself matching T1K's `int id`).
    ///
    /// # Panics
    ///
    /// Panics if `kmer_length == 0` (a zero-length k-mer is meaningless and
    /// would underflow arithmetic throughout this module, matching the
    /// undefined behavior a `kmerLength <= 0` would trigger on the C++
    /// side).
    pub fn from_reference_fasta(path: &Path, kmer_length: usize) -> Result<Self> {
        assert!(kmer_length >= 1, "kmer_length must be >= 1, got {kmer_length}");

        let text = fs::read_to_string(path)
            .with_context(|| format!("reading reference FASTA {}", path.display()))?;

        let mut seq_names = Vec::new();
        let mut seqs: Vec<Vec<u8>> = Vec::new();

        for record in parse_fasta(&text) {
            let seq_idx = seq_names.len();
            // Preserve the historical bound check (id must fit in i32) even
            // though the flat build derives ids from the `enumerate` index.
            i32::try_from(seq_idx).with_context(|| {
                format!("reference FASTA {} has more than i32::MAX sequences", path.display())
            })?;
            seq_names.push(record.id);
            seqs.push(record.seq.into_bytes());
        }

        // Build the flat k-mer index in parallel over all sequences (byte-
        // identical to the sequential per-sequence `build_index_from_read`
        // path; see `FlatKmerIndex` docs for the ordering proof).
        let seq_index = FlatKmerIndex::build(&seqs, kmer_length);

        let seq_count = seq_names.len();
        Ok(Self {
            seq_index,
            kmer_length,
            seq_count,
            seq_names,
            seqs,
            hit_len_required: DEFAULT_HIT_LEN_REQUIRED,
            ref_seq_similarity: DEFAULT_REF_SEQ_SIMILARITY,
        })
    }

    /// The k-mer length every reference sequence was indexed at.
    #[must_use]
    pub fn kmer_length(&self) -> usize {
        self.kmer_length
    }

    /// The number of reference sequences loaded (`SeqSet::Size()`).
    #[must_use]
    pub fn seq_count(&self) -> usize {
        self.seq_count
    }

    /// The name of the reference sequence at 0-based load-order index `idx`
    /// (`SeqSet::GetSeqName`, `SeqSet.hpp:810-813`).
    #[must_use]
    pub fn seq_name(&self, idx: usize) -> &str {
        &self.seq_names[idx]
    }

    /// `SeqSet::hitLenRequired` (`SeqSet.hpp:223`), used by
    /// [`RefKmerFilter::has_hit_in_set`]'s bucket-count gate.
    #[must_use]
    pub fn hit_len_required(&self) -> i32 {
        self.hit_len_required
    }

    /// Ported from `SeqSet::SetHitLenRequired` (`SeqSet.hpp:805-808`): sets
    /// `hitLenRequired`, the minimum `kmerLength * <hit count>` a candidate
    /// read's winning bucket must reach to pass gate 1 of `HasHitInSet`.
    /// `FastqExtractor.cpp:407` calls this with a data-dependent value
    /// computed from sampled read lengths (see [`crate::extract`]'s module
    /// docs for that computation) rather than relying on the
    /// [`DEFAULT_HIT_LEN_REQUIRED`] constructor default.
    pub fn set_hit_len_required(&mut self, hit_len_required: i32) {
        self.hit_len_required = hit_len_required;
    }

    /// Ported from `SeqSet::SetRefSeqSimilarity` (`SeqSet.hpp:835-838`): sets
    /// `refSeqSimilarity`, the minimum fraction of a candidate read's length
    /// that gate 2's LIS-chained mismatch-threshold check must confirm.
    /// `FastqExtractor.cpp:408` calls this with `filterAlignmentSimilarity`
    /// (the `-s` CLI flag, default 0.8) rather than relying on the
    /// [`DEFAULT_REF_SEQ_SIMILARITY`] constructor default.
    pub fn set_ref_seq_similarity(&mut self, ref_seq_similarity: f64) {
        self.ref_seq_similarity = ref_seq_similarity;
    }

    /// Ported from `SeqSet::GetRefSeqSimilarity` (`SeqSet.hpp:840-843`).
    /// Used by [`crate::genotyper::Genotyper::read_assignment_weight`]'s
    /// `segment = (1 - refSet.GetRefSeqSimilarity()) / 4.0` computation
    /// (`Genotyper.hpp:212`).
    #[must_use]
    pub fn ref_seq_similarity(&self) -> f64 {
        self.ref_seq_similarity
    }

    /// The raw reference sequence bytes at 0-based load-order index `idx`
    /// (`SeqSet::GetSeqConsensus`, `SeqSet.hpp:815-818`). Used by
    /// [`crate::genotyper::Genotyper::init_allele_info`], which needs each
    /// allele's own sequence to compute gene-level k-mer similarity and
    /// effective length -- mirroring `Genotyper::InitAlleleInfo`
    /// (`Genotyper.hpp:559-682`), which reads `refSet.GetSeqConsensus(...)`
    /// for exactly the same purpose.
    #[must_use]
    pub fn seq_consensus(&self, idx: usize) -> &[u8] {
        &self.seqs[idx]
    }

    /// Ported EXACTLY from `SeqSet::InferKmerLength` (`SeqSet.hpp:2830-2845`):
    /// sums every loaded reference sequence's length into `totalLength`, then
    /// repeatedly divides `totalLength` by 4 (integer division) counting each
    /// iteration until it reaches 0, and finally adds one more to that count.
    /// This is `floor(log4(totalLength)) + 2` for `totalLength > 0` (the loop
    /// counts base-4 "digits", i.e. `floor(log4(totalLength)) + 1`
    /// iterations, then one more `+= 1` after the loop) and exactly `1` for
    /// `totalLength == 0` (loop body never runs; the trailing `+= 1` still
    /// fires).
    ///
    /// `FastqExtractor.cpp:411` calls this AFTER `SetHitLenRequired`/
    /// `SetRefSeqSimilarity` but the total-length sum this reads
    /// (`seqs[i].consensusLen`, i.e. each loaded sequence's own length) is
    /// unaffected by either of those calls, so callers may invoke this at any
    /// point after [`RefKmerFilter::from_reference_fasta`] returns.
    #[must_use]
    pub fn infer_kmer_length(&self) -> usize {
        // Real reference FASTAs never approach i64::MAX total length, so
        // this sum cannot wrap in practice; matches the C++ `int
        // totalLength` (itself far smaller-range, `i32`) accumulation.
        #[allow(clippy::cast_possible_wrap)]
        let mut total_length: i64 = self.seqs.iter().map(|s| s.len() as i64).sum();
        let mut ret: i64 = 0;
        while total_length != 0 {
            ret += 1;
            total_length /= 4;
        }
        ret += 1;
        // `ret` is bounded by ~32 even for an astronomically large total
        // reference length (log4 of i64::MAX is well under 32), so this
        // never truncates in practice; matches the C++ `int` return type.
        usize::try_from(ret).unwrap_or(usize::MAX)
    }

    /// Ported from `SeqSet::UpdateKmerLength` (`SeqSet.hpp:2847-2857`):
    /// rebuilds `seq_index` from scratch at a new k-mer length `kl`, by
    /// clearing it and re-running [`KmerIndex::build_index_from_read`] over
    /// every stored reference sequence (in the same 0-based load order used
    /// by [`RefKmerFilter::from_reference_fasta`], so `id`s match exactly).
    ///
    /// `FastqExtractor.cpp:411-418` only calls this when `InferKmerLength()`
    /// returns a value strictly greater than the current k-mer length; this
    /// method itself does not gate on that condition (mirroring
    /// `UpdateKmerLength`'s own unconditional body) -- the caller decides
    /// whether to call it, exactly as `FastqExtractor.cpp:main` does.
    ///
    /// # Panics
    ///
    /// Panics if `kl == 0` (see [`RefKmerFilter::from_reference_fasta`]'s
    /// identical panic doc).
    pub fn update_kmer_length(&mut self, kl: usize) {
        assert!(kl >= 1, "kmer_length must be >= 1, got {kl}");
        self.kmer_length = kl;
        // Rebuild the flat index from scratch at the new k over the stored
        // sequences (same 0-based load order, so `id`s match exactly) --
        // byte-identical to the old sequential rebuild.
        self.seq_index = FlatKmerIndex::build(&self.seqs, kl);
    }

    /// Ergonomic entry point: `!is_low_complexity(read) &&
    /// has_hit_in_set(read)`, mirroring the free function `IsGoodCandidate`
    /// (`FastqExtractor.cpp:113-118`) exactly (including short-circuit
    /// order: `is_low_complexity` is always checked first, so `has_hit_in_set`
    /// is never called for a low-complexity read).
    ///
    /// Allocates a fresh [`Scratch`] internally on every call; for
    /// processing many reads in a loop, prefer
    /// [`RefKmerFilter::is_good_candidate_with_scratch`] with a
    /// caller-owned, reused `Scratch`.
    #[must_use]
    pub fn is_good_candidate(&self, read: &[u8]) -> bool {
        let mut scratch = Scratch::default();
        self.is_good_candidate_with_scratch(read, &mut scratch)
    }

    /// Same as [`RefKmerFilter::is_good_candidate`], but reuses a
    /// caller-owned [`Scratch`] instead of allocating one per call -- the
    /// intended hot-loop entry point for a future `extract`/`genotype`
    /// command processing many reads.
    #[must_use]
    pub fn is_good_candidate_with_scratch(&self, read: &[u8], scratch: &mut Scratch) -> bool {
        !is_low_complexity(read) && self.has_hit_in_set(read, scratch)
    }

    /// Diagnostic-only: evaluates ONLY `HasHitInSet`'s bucket-count gate
    /// (`SeqSet.hpp:1929-1964`, i.e. exactly what Task 3.1's
    /// `has_hit_in_set` computed before this task added gate 2) and returns
    /// whether it alone would accept `read`, WITHOUT running gate 2's
    /// `GetOverlapsFromHits`-based confirmation.
    ///
    /// This exists purely so a differential test can quantify how many reads
    /// gate 2 actually rejects that gate 1 alone would have accepted (i.e.
    /// how many reads exercise the fix this task makes) -- it is not used by
    /// [`RefKmerFilter::has_hit_in_set`]/[`RefKmerFilter::is_good_candidate`]
    /// itself (both call [`RefKmerFilter::bucket_count_gate`] directly, not
    /// this wrapper). See `crates/unum-core/tests/golden_refkmerfilter.rs`
    /// for the consumer.
    #[must_use]
    pub fn passes_bucket_count_gate_only(&self, read: &[u8], scratch: &mut Scratch) -> bool {
        self.bucket_count_gate(read, scratch).is_some()
    }

    /// Ported from the bucket-count gate at the top of `SeqSet::HasHitInSet`
    /// (`SeqSet.hpp:1929-1964`): collects hits from both strands, buckets
    /// them by `(strand-tag, seqIdx)`, and finds the single largest bucket.
    /// Returns `None` if that gate rejects the read (mirrors stock's early
    /// `return false`), otherwise `Some((winning_bucket, bucket_size))`.
    fn bucket_count_gate(&self, read: &[u8], scratch: &mut Scratch) -> Option<((i8, u32), u32)> {
        let len = read.len();
        if len < self.kmer_length {
            return None;
        }

        self.get_hits_from_read(read, scratch);
        if scratch.hits.is_empty() {
            return None;
        }

        // Bucket sort: (tag, seqIdx) -> hit count, mirroring stock's
        // `buckets[tag][idx].PushBack(...)` (SeqSet.hpp:1935-1939), except
        // we only need each bucket's SIZE to pick the winner (gate 2 below
        // reconstructs the winning bucket's actual MEMBER hits separately,
        // by filtering `scratch.hits` -- see the comment at that call site),
        // so a count map suffices here.
        scratch.buckets.clear();
        for hit in &scratch.hits {
            let tag = i8::from(hit.strand == 1);
            *scratch.buckets.entry((tag, hit.idx)).or_insert(0) += 1;
        }

        // add02ca touched-buckets selection -- see `select_best_bucket`'s
        // doc comment for why this is byte-identical in outcome to stock's
        // O(seqCnt) full-array scan (SeqSet.hpp:1941-1957).
        let (best_bucket, max) = select_best_bucket(&scratch.buckets)?;

        // SeqSet.hpp:1959: `if (kmerLength * max < hitLenRequired) return false;`
        // Widened to i64 to avoid any theoretical `i32` overflow from the
        // multiplication (stock uses plain `int` arithmetic here; `max` and
        // `kmerLength` are realistically small enough that this never
        // matters in practice, but i64 costs nothing and removes the
        // question).
        let kmer_length_i64 = i64::try_from(self.kmer_length).unwrap_or(i64::MAX);
        if kmer_length_i64 * i64::from(max) < i64::from(self.hit_len_required) {
            return None;
        }

        Some((best_bucket, max))
    }

    /// Ported from `SeqSet::HasHitInSet` in full (`SeqSet.hpp:1915-1990`):
    /// the bucket-count gate ([`RefKmerFilter::bucket_count_gate`]) followed
    /// by the `GetOverlapsFromHits`-based mismatch-threshold confirmation.
    /// See the module docs for the overall two-gate structure and
    /// [`crate::overlap`]'s module docs for why gate 2 needs no
    /// `AlignAlgo`/Smith-Waterman.
    ///
    /// Exposed (rather than fully private) so the golden test
    /// (`crates/unum-core/tests/golden_refkmerfilter.rs`) can drive it
    /// directly, matching the established pattern of exposing internals to
    /// the same-workspace golden tests without making them a documented part
    /// of the public API surface.
    #[must_use]
    pub(crate) fn has_hit_in_set(&self, read: &[u8], scratch: &mut Scratch) -> bool {
        let len = read.len();
        let Some((best_bucket, _max)) = self.bucket_count_gate(read, scratch) else {
            return false;
        };

        // Gate 2 (SeqSet.hpp:1966-1983): GetOverlapsFromHits on the winning
        // bucket, then accept if any resulting overlap's implied mismatch
        // count is within the similarity-derived threshold.
        //
        // Reconstructing the winning bucket's MEMBER hits: stock builds each
        // `buckets[tag][idx]` by a single forward pass over the full hit
        // list, `PushBack`-ing each hit into its `(tag, idx)` bucket in
        // encounter order (SeqSet.hpp:1935-1939) -- i.e. a stable filter of
        // the full hit list by `(tag, idx)`. `scratch.hits` (filled by
        // `get_hits_from_read` just above) is exactly that full hit list, in
        // the same forward-strand-then-reverse-strand encounter order stock
        // produces it in, so filtering it by `(tag, idx) == best_bucket`
        // here reproduces stock's bucket contents (and their order)
        // exactly, without needing to build all `2 * seqCount` buckets
        // up front.
        let (best_tag, best_idx) = best_bucket;
        let winning_bucket_hits: Vec<OverlapHit> = scratch
            .hits
            .iter()
            .filter(|hit| i8::from(hit.strand == 1) == best_tag && hit.idx == best_idx)
            .map(|hit| hit.to_overlap_hit())
            .collect();

        let overlaps = overlap::get_overlaps_from_hits(
            &winning_bucket_hits,
            self.hit_len_required,
            0,           // filter: HasHitInSet always passes filter=0 (SeqSet.hpp:1968).
            |_idx| true, // is_ref: every sequence loaded via from_reference_fasta is a reference (see module docs).
            DEFAULT_RADIUS,
            IS_LONG_SEQ_SET,
            #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
            {
                self.kmer_length as i32
            },
        );

        // SeqSet.hpp:1973: `int mismatchThreshold = int(len * (1 -
        // refSeqSimilarity)) * kmerLength;` -- truncates `len * (1 -
        // refSeqSimilarity)` to `int` FIRST, THEN multiplies by
        // `kmerLength` (NOT `int(len * (1 - refSeqSimilarity) *
        // kmerLength)`); reproduced in that exact operand order below.
        #[allow(clippy::cast_precision_loss)]
        let len_f64 = len as f64;
        #[allow(clippy::cast_possible_truncation)]
        let truncated = (len_f64 * (1.0 - self.ref_seq_similarity)) as i32;
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let kmer_length_i32 = self.kmer_length as i32;
        let mismatch_threshold = truncated * kmer_length_i32;

        // SeqSet.hpp:1974-1981: `valid = true` if ANY overlap satisfies
        // `len - overlaps[i].matchCnt / 2 <= mismatchThreshold` (integer
        // `matchCnt / 2`, i.e. `matchCnt` counted TWICE per hitLen and
        // halved back here -- see `overlap::Overlap::match_cnt`'s doc
        // comment).
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let read_len_i32 = len as i32;
        overlaps.iter().any(|o| read_len_i32 - o.match_cnt / 2 <= mismatch_threshold)
    }

    /// Ported from `SeqSet::GetOverlapsFromRead` (`SeqSet.hpp:1594-1912`):
    /// the read-to-allele alignment/scoring core the Genotyper (a later
    /// phase) depends on. Turns a single read into a list of scored,
    /// `AlignAlgo`-confirmed [`overlap::Overlap`]s -- unlike
    /// [`RefKmerFilter::has_hit_in_set`]'s gate 2 (which only ever reads
    /// `matchCnt` from the LIS-chained hit length, `overlap::Overlap::match_cnt`'s
    /// doc comment), this method fully reproduces stock's per-overlap
    /// `AlignAlgo`-based `matchCnt`/`similarity` refinement.
    ///
    /// Returns `None` if `read.len() < kmer_length` (mirrors stock's `return
    /// -1`, `SeqSet.hpp:1598-1599` -- surfaced here as `Option` rather than a
    /// sentinel `-1`, since this port's caller always has a real `Vec` to
    /// inspect either way).
    ///
    /// # `isRef` is always `true` for every seq loaded via `from_reference_fasta`
    ///
    /// Every reference sequence loaded by [`RefKmerFilter::from_reference_fasta`]
    /// is `isRef == true` (see module docs), and a ref-only `SeqSet` never
    /// populates `posWeight` -- so the `!isRef` / `GlobalAlignment_PosWeight`
    /// branches below are ported faithfully (not overfit away) but are
    /// PROVABLY UNREACHABLE from any overlap this method can actually
    /// produce in this codebase today: [`RefKmerFilter::is_ref`] (this
    /// type's per-seq `isRef` query, mirroring `_seqWrapper::isRef`) always
    /// returns `true`. A future caller that builds novel (non-reference)
    /// sequences -- e.g. a Genotyper port -- would need to extend
    /// [`RefKmerFilter`] with a real per-seq `isRef`/`posWeight` model
    /// (replacing [`RefKmerFilter::is_ref`]'s hardcoded `true` and the
    /// placeholder all-zero-count `PosWeight` this method substitutes below)
    /// before those branches become live.
    ///
    /// # FLOATS: `similarity` is a plain, deterministic `f64` ratio
    ///
    /// `similarity = (double)matchCnt / (seqSpan + readSpan)` (`SeqSet.hpp:
    /// 1838-1840`) is reproduced as `f64::from(match_cnt) /
    /// f64::from(seq_span + read_span)` in that exact operand order --
    /// bit-identical to the C++, not merely close. See [`overlap::overlap_less_than`]'s
    /// doc comment for the same exactness requirement applied to
    /// `_overlap::operator<`'s `similarity` comparison.
    #[must_use]
    // `int_plus_one`: this function's `... + kmerLength - 1 >= ...` guards are
    // deliberately kept in the C++'s own `x - 1 >= y` form (rather than
    // clippy's suggested `x > y`) throughout, to stay line-comparable to
    // `SeqSet.hpp:1704,1770-1786` for future maintainers cross-checking a fix
    // against the vendored source -- see `align_algo.rs`'s module docs for
    // the same "line-comparable over clippy's preferences" convention
    // applied elsewhere in this port.
    #[allow(clippy::too_many_lines, clippy::int_plus_one)]
    pub fn get_overlaps_from_read(
        &self,
        read: &[u8],
        scratch: &mut Scratch,
    ) -> Option<Vec<overlap::Overlap>> {
        let len = read.len();
        if len < self.kmer_length {
            return None;
        }

        self.get_hits_from_read(read, scratch);
        sort_hits(&mut scratch.hits, true, self.seq_count);

        let winning_hits: Vec<OverlapHit> =
            scratch.hits.iter().map(|hit| hit.to_overlap_hit()).collect();

        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let kmer_length_i32 = self.kmer_length as i32;
        let overlaps_with_coords = overlap::get_overlaps_from_hits_with_coords(
            &winning_hits,
            self.hit_len_required,
            0, // filter: GetOverlapsFromRead always passes filter=0 (SeqSet.hpp:1614).
            |idx| self.is_ref(idx),
            DEFAULT_RADIUS,
            IS_LONG_SEQ_SET,
            kmer_length_i32,
        );

        let mut overlaps: Vec<overlap::Overlap> =
            overlaps_with_coords.iter().map(|(o, _)| *o).collect();
        let mut overlaps_hit_coords: Vec<Vec<overlap::Pair>> =
            overlaps_with_coords.into_iter().map(|(_, coords)| coords).collect();

        // Strand filter (SeqSet.hpp:1619-1648): find bestOverlapIdx via
        // overlap::overlap_less_than (the FIRST index achieving the minimum
        // under that ordering -- `<` strictly, so later ties do not
        // displace an earlier winner, matching C++'s `if (overlaps[i] <
        // overlaps[bestOverlapIdx]) bestOverlapIdx = i` loop exactly), then
        // compact `overlaps`/`overlaps_hit_coords` in place to keep only
        // entries whose strand matches the winner's.
        if !overlaps.is_empty() {
            let mut best_overlap_idx = 0usize;
            for i in 1..overlaps.len() {
                if overlap::overlap_less_than(&overlaps[i], &overlaps[best_overlap_idx]) {
                    best_overlap_idx = i;
                }
            }
            let best_strand = overlaps[best_overlap_idx].strand;

            let mut k = 0usize;
            for i in 0..overlaps.len() {
                if overlaps[i].strand != best_strand {
                    continue;
                }
                if i != k {
                    overlaps[k] = overlaps[i];
                    overlaps_hit_coords[k] = std::mem::take(&mut overlaps_hit_coords[i]);
                }
                k += 1;
            }
            overlaps.truncate(k);
            overlaps_hit_coords.truncate(k);
        }

        reverse_complement_into(read, &mut scratch.rc_buf);
        let rc_read = scratch.rc_buf.clone();

        // Per-overlap AlignAlgo-based matchCnt/similarity refinement
        // (SeqSet.hpp:1660-1882). Stock also maintains `firstRef`/
        // `bestNovelOverlap`/`readOverlapRepresentatives` here, but -- as in
        // stock -- they never feed back into which overlaps are kept (only the
        // final refSeqSimilarity/novelSeqSimilarity filter below does that),
        // so their bookkeeping (a `Vec` allocation plus an O(overlaps^2)
        // containment scan) is omitted: it produced values this port's `let _
        // = ...` sinks proved were never read. See this crate's history for
        // the faithful-but-dead port that this removal replaces.
        for i in 0..overlaps.len() {
            let r: &[u8] = if overlaps[i].strand == 1 { read } else { &rc_read };
            let hit_coords = &overlaps_hit_coords[i];
            let hit_cnt = hit_coords.len();

            let mut match_cnt: i32 = 0;
            let mut mismatch_cnt: i32 = 0;
            let mut indel_cnt: i32 = 0;
            let mut similarity: f64 = 1.0;

            let is_ref = self.is_ref(overlaps[i].seq_idx);

            match_cnt += 2 * kmer_length_i32;

            for j in 1..hit_cnt {
                let prev = hit_coords[j - 1];
                let cur = hit_coords[j];

                if prev.b - prev.a == cur.b - cur.a {
                    // Same coordinate diff: colinear on both read and seq.
                    if prev.a + kmer_length_i32 - 1 >= cur.a {
                        match_cnt += 2 * (cur.a - prev.a);
                    } else {
                        match_cnt += 2 * kmer_length_i32;

                        let align = if is_ref {
                            let seq_slice = seq_gap_slice(
                                &self.seqs[overlaps[i].seq_idx as usize],
                                prev.b + kmer_length_i32,
                                cur.b - (prev.b + kmer_length_i32),
                            );
                            let read_slice = read_gap_slice(
                                r,
                                prev.a + kmer_length_i32,
                                cur.a - (prev.a + kmer_length_i32),
                            );
                            crate::align_algo::global_alignment_cached(
                                seq_slice,
                                read_slice,
                                crate::align_algo::DEFAULT_BAND,
                                &mut scratch.dp_cache,
                            )
                        } else {
                            // Novel-seq path: ported but unreachable from
                            // this codebase's overlaps (see this method's
                            // doc comment). `posWeight` is not modeled by
                            // `RefKmerFilter`, so this substitutes an
                            // all-zero-count (i.e. "no support", always
                            // `is_base_equal == true`) posWeight slice of
                            // the matching length -- the correct behavior
                            // once a real posWeight source exists is to
                            // replace this placeholder, not to change the
                            // control flow here.
                            let gap_len = cur.b - (prev.b + kmer_length_i32);
                            #[allow(clippy::cast_sign_loss)]
                            let weights = vec![
                                crate::align_algo::PosWeight::default();
                                gap_len.max(0) as usize
                            ];
                            let read_slice = read_gap_slice(
                                r,
                                prev.a + kmer_length_i32,
                                cur.a - (prev.a + kmer_length_i32),
                            );
                            crate::align_algo::global_alignment_pos_weight(&weights, read_slice)
                        };

                        let (mut m, mut mm, mut ind) = (0, 0, 0);
                        crate::align_algo::get_align_stats(
                            &align.align,
                            false,
                            &mut m,
                            &mut mm,
                            &mut ind,
                        );
                        match_cnt += 2 * m;
                        mismatch_cnt += mm;
                        indel_cnt += ind;

                        if (DEFAULT_RADIUS == 0 || !is_ref) && indel_cnt > 0 {
                            similarity = 0.0;
                            break;
                        }
                    }
                } else {
                    // Different coordinate diff: non-colinear hit pair.
                    if DEFAULT_RADIUS == 0 || !is_ref {
                        similarity = 0.0;
                        break;
                    }

                    if prev.a + kmer_length_i32 - 1 >= cur.a && prev.b + kmer_length_i32 - 1 < cur.b
                    {
                        match_cnt += 2 * (cur.a - prev.a);
                        indel_cnt += (cur.b - (prev.b + kmer_length_i32))
                            + (cur.a + kmer_length_i32 - prev.a);
                    } else if prev.a + kmer_length_i32 - 1 < cur.a
                        && prev.b + kmer_length_i32 - 1 >= cur.b
                    {
                        match_cnt += 2 * (cur.b - prev.b);
                        indel_cnt += (cur.a - (prev.a + kmer_length_i32))
                            + (cur.b + kmer_length_i32 - prev.b);
                    } else if prev.a + kmer_length_i32 - 1 >= cur.a
                        && prev.b + kmer_length_i32 - 1 >= cur.b
                    {
                        match_cnt += 2 * (cur.a - prev.a).min(cur.b - prev.b);
                        indel_cnt += ((cur.a - cur.b) - (prev.a - prev.b)).abs();
                    } else {
                        match_cnt += 2 * kmer_length_i32;

                        let align = if is_ref {
                            let seq_slice = seq_gap_slice(
                                &self.seqs[overlaps[i].seq_idx as usize],
                                prev.b + kmer_length_i32,
                                cur.b - (prev.b + kmer_length_i32),
                            );
                            let read_slice = read_gap_slice(
                                r,
                                prev.a + kmer_length_i32,
                                cur.a - (prev.a + kmer_length_i32),
                            );
                            crate::align_algo::global_alignment_cached(
                                seq_slice,
                                read_slice,
                                crate::align_algo::DEFAULT_BAND,
                                &mut scratch.dp_cache,
                            )
                        } else {
                            let gap_len = cur.b - (prev.b + kmer_length_i32);
                            #[allow(clippy::cast_sign_loss)]
                            let weights = vec![
                                crate::align_algo::PosWeight::default();
                                gap_len.max(0) as usize
                            ];
                            let read_slice = read_gap_slice(
                                r,
                                prev.a + kmer_length_i32,
                                cur.a - (prev.a + kmer_length_i32),
                            );
                            crate::align_algo::global_alignment_pos_weight(&weights, read_slice)
                        };

                        let (mut m, mut mm, mut ind) = (0, 0, 0);
                        crate::align_algo::get_align_stats(
                            &align.align,
                            false,
                            &mut m,
                            &mut mm,
                            &mut ind,
                        );
                        match_cnt += 2 * m;
                        mismatch_cnt += mm;
                        indel_cnt += ind;

                        if !is_ref && indel_cnt > 0 {
                            similarity = 0.0;
                            break;
                        }
                    }
                }
            }
            let _ = mismatch_cnt; // matches stock: computed, never read after the loop.

            overlaps[i].match_cnt = match_cnt;
            if (similarity - 1.0).abs() < f64::EPSILON {
                let seq_span = overlaps[i].seq_end - overlaps[i].seq_start + 1;
                let read_span = overlaps[i].read_end - overlaps[i].read_start + 1;
                overlaps[i].similarity = f64::from(match_cnt) / f64::from(seq_span + read_span);
            } else {
                overlaps[i].similarity = 0.0;
            }

            if is_overlap_low_complex(r, &overlaps[i]) {
                overlaps[i].similarity = 0.0;
            }
            overlaps[i].match_cnt = match_cnt;
        }

        // Final filter (SeqSet.hpp:1893-1908): keep only overlaps meeting
        // the isRef-vs-novel similarity threshold.
        overlaps.retain(|o| {
            let is_ref = self.is_ref(o.seq_idx);
            if is_ref {
                o.similarity >= self.ref_seq_similarity
            } else {
                o.similarity >= DEFAULT_NOVEL_SEQ_SIMILARITY
            }
        });

        Some(overlaps)
    }

    /// Every sequence loaded via [`RefKmerFilter::from_reference_fasta`] is a
    /// reference sequence (`_seqWrapper::isRef == true`, `SeqSet.hpp:883,911`)
    /// -- see module docs. `_idx` is accepted (rather than an always-`true`
    /// zero-arg closure) purely to match the `impl Fn(u32) -> bool` shape
    /// [`overlap::get_overlaps_from_hits`]/[`overlap::get_overlaps_from_hits_with_coords`]
    /// expect, and to document the extension point for a future novel-seq
    /// model (see [`RefKmerFilter::get_overlaps_from_read`]'s doc comment).
    #[allow(clippy::unused_self)]
    fn is_ref(&self, _idx: u32) -> bool {
        true
    }

    /// Ported from `SeqSet::GetHitsFromRead` (`SeqSet.hpp:1071-1229`),
    /// specialized to exactly the parameters `HasHitInSet` always passes
    /// (`SeqSet.hpp:1923`: `strand=0`, `barcode=-1`, `allowTotalSkip=false`,
    /// `puse=NULL`) -- both strands are searched, no barcode/`puse`
    /// filtering applies, and `allowTotalSkip`'s branch (dead under
    /// `allowTotalSkip=false`) is omitted. Fills `scratch.hits` (cleared
    /// first) and `scratch.rc_buf` (the read's reverse complement, computed
    /// as an output -- mirroring T1K's `rcRead` parameter to
    /// `GetHitsFromRead`/`HasHitInSet`, which callers pass as a pre-allocated
    /// scratch buffer purely to receive the computed reverse complement, not
    /// as an input).
    ///
    /// # `downSample` is always 1 in this port (dead branch omitted)
    ///
    /// Stock's `downSample` (`SeqSet.hpp:1084-1091`) is only ever `> 1` when
    /// `len > 200 && isLongSeqSet`. `isLongSeqSet` defaults to `false`
    /// (`SeqSet.hpp:765`) and has no public setter reachable from the plain
    /// `SeqSet(kmerLength)` construction this port's `from_reference_fasta`
    /// performs -- so `downSample` is unconditionally `1` for every
    /// `RefKmerFilter`, making the downsample-skip branch (`SeqSet.hpp:
    /// 1101-1102,1168-1169`: `if (downSample > 1 && ...) continue;`)
    /// unreachable dead code for this port's scope. It is omitted here
    /// (rather than implemented-but-inert) to keep the loop readable; if a
    /// future caller needs `isLongSeqSet` support, this is the first place
    /// to revisit.
    ///
    /// # `prevKmerCode` is NOT reset between the forward and reverse loops
    ///
    /// Stock only calls `kmerCode.Restart()` before the reverse-strand loop
    /// (`SeqSet.hpp:1160`) -- `prevKmerCode` (the consecutive-duplicate-
    /// window dedup state) carries over from wherever the forward loop left
    /// it. Since the reverse loop's very first search is always forced
    /// anyway (`i == kmerLength - 1` short-circuits the `IsEqual` check),
    /// this stale carry-over is only observable if that very first reverse
    /// window happens to be skipped by the high-frequency-kmer skip-limit
    /// branch below -- in which case `prevKmerCode` stays whatever the
    /// forward loop last left it, exactly matching stock. This port
    /// reproduces that by simply not resetting `prev_kmer_code` between the
    /// two loops (only `kmer_code` gets a fresh `KmerCode::new` +
    /// re-filled window, mirroring `Restart()`).
    fn get_hits_from_read(&self, read: &[u8], scratch: &mut Scratch) {
        scratch.hits.clear();
        let len = read.len();
        let kl = self.kmer_length;

        if len < kl {
            // Defensive guard not present in the C++ (which has no such
            // check inside `GetHitsFromRead` itself -- only `HasHitInSet`'s
            // caller-side `len < kmerLength` check protects it in practice,
            // SeqSet.hpp:1919-1920). Calling this directly with a short read
            // would index `read[i]` for `i >= len` in the fill loop below,
            // which Rust would panic on (a safe failure) where the C++
            // reads past the string's NUL terminator (undefined behavior).
            // `has_hit_in_set` (the only caller in this port) already
            // guards this exact condition before calling, so this is purely
            // a safety net for direct callers of this `pub(crate)` helper
            // (e.g. tests), not an observable behavior change on the tested
            // path.
            return;
        }

        let mut kmer_code = KmerCode::new(kl);
        let mut prev_kmer_code = KmerCode::new(kl);
        let skip_limit = kl / 2;
        let mut skip_cnt = 0usize;

        // Forward strand (SeqSet.hpp:1093-1154).
        let mut i = 0usize;
        while i < kl - 1 {
            kmer_code.append(read[i]);
            i += 1;
        }
        while i < len {
            kmer_code.append(read[i]);
            if i == kl - 1 || !prev_kmer_code.is_equal(&kmer_code) {
                let found = self.seq_index.search(&kmer_code);
                let size = found.len();
                if size >= 100 && i != kl - 1 && i != len - 1 && skip_cnt < skip_limit {
                    skip_cnt += 1;
                    i += 1;
                    continue; // matches C++: also skips the prev_kmer_code update below
                }
                skip_cnt = 0;
                // `i`/`kl` are window positions bounded by `len` (a real
                // read length), so this never approaches i32::MAX; matches
                // the equivalent cast in `kmer_index::build_index_from_read`.
                #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                let read_offset = (i - (kl - 1)) as i32;
                // `repeats = size` (SeqSet.hpp:1122): `HasHitInSet`'s
                // `GetHitsFromRead` call always has `puse == NULL` and
                // `barcode == -1` (SeqSet.hpp:1923), so neither the
                // `puse`-filtered recount (SeqSet.hpp:1123-1132) nor the
                // `barcode != -1` override (SeqSet.hpp:1134-1135) ever
                // applies here -- `repeats` is simply this k-mer's total
                // index hit count.
                #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                let repeats = size as i32;
                for entry in found {
                    scratch.hits.push(Hit {
                        idx: entry.idx,
                        offset: entry.offset,
                        read_offset,
                        strand: 1,
                        repeats,
                    });
                }
            }
            prev_kmer_code = kmer_code.clone();
            i += 1;
        }

        // Reverse complement (SeqSet.hpp:1156, always computed regardless
        // of which strand(s) are searched -- dead work in stock whenever
        // `strand == 1`, but this port only ever needs `strand == 0`, so
        // it's always used here).
        reverse_complement_into(read, &mut scratch.rc_buf);

        // Reverse strand (SeqSet.hpp:1158-1227). `kmer_code.Restart()`
        // equivalent: fresh code/invalid_pos at the same kmer_length.
        // `prev_kmer_code` is deliberately NOT reset -- see this function's
        // doc comment.
        kmer_code = KmerCode::new(kl);
        skip_cnt = 0; // stock explicitly re-zeroes skipCnt here (SeqSet.hpp:1164)
        let mut i = 0usize;
        while i < kl - 1 {
            kmer_code.append(scratch.rc_buf[i]);
            i += 1;
        }
        while i < len {
            kmer_code.append(scratch.rc_buf[i]);
            if i == kl - 1 || !prev_kmer_code.is_equal(&kmer_code) {
                let found = self.seq_index.search(&kmer_code);
                let size = found.len();
                if size >= 100 && i != kl - 1 && i != len - 1 && skip_cnt < skip_limit {
                    skip_cnt += 1;
                    i += 1;
                    continue;
                }
                skip_cnt = 0;
                // `i`/`kl` are window positions bounded by `len` (a real
                // read length), so this never approaches i32::MAX; matches
                // the equivalent cast in `kmer_index::build_index_from_read`.
                #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                let read_offset = (i - (kl - 1)) as i32;
                // `repeats = size` -- see the forward-strand loop's identical
                // comment above (SeqSet.hpp:1189, same reasoning).
                #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                let repeats = size as i32;
                for entry in found {
                    scratch.hits.push(Hit {
                        idx: entry.idx,
                        offset: entry.offset,
                        read_offset,
                        strand: -1,
                        repeats,
                    });
                }
            }
            prev_kmer_code = kmer_code.clone();
            i += 1;
        }
    }
}

/// Selects the single best `(tag, seqIdx)` bucket -- the one with the most
/// hits -- from a `(tag, seqIdx) -> count` map, replicating stock T1K's
/// `HasHitInSet` bucket-selection loop (`SeqSet.hpp:1941-1957`) bit-for-bit,
/// but WITHOUT stock's `O(seqCnt)` full-array scan. This is the "add02ca"
/// performance fold-in this port requires (named for the maintainer commit
/// that introduced the equivalent optimization): stock allocates
/// `SimpleVector<_hit> buckets[2][seqCnt]` (`SeqSet.hpp:1931-1933`) and scans
/// every one of those `2 * seqCnt` slots regardless of how many of them are
/// actually non-empty; a read only ever touches `hitCnt` distinct buckets
/// (`hitCnt` = the number of k-mer hits collected, independent of `seqCnt`),
/// so tracking only the touched buckets in a `HashMap` and scanning those is
/// `O(hitCnt)` instead of `O(seqCnt)`.
///
/// # Ordering and tie-breaking are what make this byte-identical to stock
///
/// Stock's selection loop is:
/// ```text
/// for k in 0..=1 {           // tag
///     for i in 0..seqCnt {   // seqIdx
///         if size(buckets[k][i]) > 0 && size(buckets[k][i]) > max {
///             maxTag = k; maxSeqIdx = i; max = size(buckets[k][i]);
///         }
///     }
/// }
/// ```
/// i.e. an ascending `(tag, seqIdx)` scan with a **strict `>`** update
/// condition, so the FIRST bucket (in `(tag, seqIdx)` order) to reach the
/// running maximum wins -- later buckets with an EQUAL count never displace
/// it. Every untouched `(tag, seqIdx)` pair has `size == 0` and can
/// therefore never satisfy `size > 0`, so it can never win regardless of how
/// `max` is initialized; restricting the scan to only the touched keys does
/// not change which bucket is selected, AS LONG AS the touched keys are
/// still visited in the SAME ascending `(tag, seqIdx)` order with the SAME
/// strict `>` tie-break. This function does exactly that: it collects the
/// touched keys, sorts them ascending by `(tag, seqIdx)` (Rust's default
/// tuple `Ord` is exactly the lexicographic `(tag, seqIdx)` order stock's
/// nested loop visits), and scans with `count > best_count` (strict).
/// The result -- which specific `(tag, seqIdx)` bucket is "the" winner among
/// ties -- is therefore identical to stock's, not just the winning count.
///
/// (This port's [`RefKmerFilter::has_hit_in_set`] only actually consumes the
/// winning COUNT, not the winning bucket's identity, since the identity is
/// only needed by the excluded `GetOverlapsFromHits` step -- but the
/// identity is computed and returned here anyway, both because it costs
/// nothing extra and because a future `GetOverlapsFromHits` port will need
/// it.)
fn select_best_bucket(buckets: &HashMap<(i8, u32), u32>) -> Option<((i8, u32), u32)> {
    let mut keys: Vec<(i8, u32)> = buckets.keys().copied().collect();
    keys.sort_unstable();

    let mut best: Option<((i8, u32), u32)> = None;
    for key in keys {
        let count = buckets[&key];
        let is_new_max = match best {
            None => true,
            Some((_, best_count)) => count > best_count,
        };
        if is_new_max {
            best = Some((key, count));
        }
    }
    best
}

/// Writes the reverse complement of `seq` into `out` (resizing as needed),
/// mirroring `SeqSet::ReverseComplement` (`SeqSet.hpp:2103-2114`) exactly:
/// each base is complemented (`A<->T`, `C<->G`) except `N`, which maps to
/// `N`. Only `A`/`C`/`G`/`T`/`N` are supported (matching this module's
/// general uppercase-ACGTN-only assumption, documented at
/// [`RefKmerFilter::from_reference_fasta`]); any other byte is out of scope
/// (the C++ original reads `nucToNum[c - 'A']` unconditionally for
/// non-`N` bytes, which is undefined behavior for bytes outside `'A'..='Z'`
/// and produces a garbage-but-defined base for other in-range letters --
/// neither of which this port replicates).
fn reverse_complement_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.resize(seq.len(), 0);
    for (out_i, &c) in seq.iter().rev().enumerate() {
        out[out_i] = complement_base(c);
    }
}

/// Complements a single base, matching `numToNuc[3 - nucToNum[c - 'A']]`
/// (with `N` bypassing the table via `SeqSet::ReverseComplement`'s explicit
/// `if (seq[...] != 'N')` check, `SeqSet.hpp:2108-2111`).
fn complement_base(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        other => panic!(
            "reverse_complement_into: unsupported base {other:?} (only A/C/G/T/N are supported; \
             see this module's doc comment)"
        ),
    }
}

/// Slices `seq_len` bytes out of `seq` starting at signed offset `start`,
/// mirroring the C++ pointer arithmetic `seqs[...].consensus + start` (a
/// `char*`) paired with a separately-computed length argument
/// (`GetOverlapsFromRead`'s `AlignAlgo::GlobalAlignment(consensus + ..., len,
/// ...)` call sites, `SeqSet.hpp:1716-1721` etc.). Both `start` and `seq_len`
/// are always non-negative for every reachable call site in
/// [`RefKmerFilter::get_overlaps_from_read`] (see that method's inline
/// comments deriving the non-negativity from each branch's guard condition);
/// this helper still uses checked `usize` conversions (rather than
/// `as usize`) so a violated invariant panics loudly instead of silently
/// wrapping.
///
/// # Panics
///
/// Panics if `start`/`seq_len` are negative, or if `start + seq_len` exceeds
/// `seq.len()`.
fn seq_gap_slice(seq: &[u8], start: i32, seq_len: i32) -> &[u8] {
    let start = usize::try_from(start).expect("seq_gap_slice: start must be non-negative");
    let seq_len = usize::try_from(seq_len).expect("seq_gap_slice: seq_len must be non-negative");
    &seq[start..start + seq_len]
}

/// Identical to [`seq_gap_slice`], but named separately for readability at
/// call sites that slice the read (`r + ...`) rather than the reference
/// sequence's consensus (`seqs[...].consensus + ...`) -- both C++ pointer
/// arithmetic idioms are structurally the same, but naming them apart keeps
/// [`RefKmerFilter::get_overlaps_from_read`]'s call sites line-comparable to
/// which C++ pointer (`consensus` vs. `r`) each one ports.
///
/// # Panics
///
/// Panics if `start`/`read_len` are negative, or if `start + read_len`
/// exceeds `r.len()`.
fn read_gap_slice(r: &[u8], start: i32, read_len: i32) -> &[u8] {
    let start = usize::try_from(start).expect("read_gap_slice: start must be non-negative");
    let read_len =
        usize::try_from(read_len).expect("read_gap_slice: read_len must be non-negative");
    &r[start..start + read_len]
}

/// Ported from `SeqSet::IsOverlapLowComplex` (`SeqSet.hpp:458-485`): flags an
/// overlap's read span (`r[o.readStart..=o.readEnd]`) as low-complexity by a
/// DIFFERENT (looser) rule than the whole-read [`is_low_complexity`]: counts
/// how many of the 4 bases appear `<= 2` times within just that span, and
/// how many total bases those low-count bases contribute. If that low-count
/// total is at least `1/7` of the span length, the span is considered
/// "not low-complex enough to reject" and this returns `false`
/// unconditionally (even if `lowCnt >= 2`); otherwise, `lowCnt >= 2` (at
/// least two bases each appearing `<= 2` times, contributing under `1/7` of
/// the span) flags it as low-complex.
///
/// Only `A`/`C`/`G`/`T`/`N` are supported (`N` positions are skipped
/// entirely, matching `SeqSet.hpp:464-465`'s `if (r[i]=='N') continue;`; any
/// other byte hits [`nuc_index`]'s same unsupported-base panic as
/// [`is_low_complexity`]).
///
/// # Panics
///
/// Panics if `o.read_start`/`o.read_end` are out of bounds for `r`, or if
/// `r` contains a byte other than `A`/`C`/`G`/`T`/`N` in that range.
fn is_overlap_low_complex(r: &[u8], o: &overlap::Overlap) -> bool {
    let mut cnt: [i32; 4] = [0; 4];
    #[allow(clippy::cast_sign_loss)]
    let (start, end) = (o.read_start as usize, o.read_end as usize);
    for &c in &r[start..=end] {
        if c == b'N' {
            continue;
        }
        cnt[nuc_index(c)] += 1;
    }
    let len = o.read_end - o.read_start + 1;
    let mut low_cnt = 0i32;
    let mut low_total_cnt = 0i32;
    for &c in &cnt {
        if c <= 2 {
            low_cnt += 1;
            low_total_cnt += c;
        }
    }
    if low_total_cnt * 7 >= len {
        return false;
    }
    low_cnt >= 2
}

/// Ported verbatim from the free function `IsLowComplexity`
/// (`FastqExtractor.cpp:89-111`, byte-identical to the copy at
/// `BamExtractor.cpp:144-166`): flags a read as low-complexity if any single
/// base makes up at least half the read, if `N`s make up at least a tenth of
/// the read, OR if at least two of the four bases each appear two times or
/// fewer in the whole read.
///
/// Only `A`/`C`/`G`/`T`/`N` are supported (see [`complement_base`]'s doc
/// comment for the same caveat, which applies identically here since both
/// functions mirror the same `nucToNum[c - 'A']` indexing pattern).
#[must_use]
pub fn is_low_complexity(seq: &[u8]) -> bool {
    let mut cnt: [i32; 5] = [0; 5];
    for &c in seq {
        if c == b'N' {
            cnt[4] += 1;
        } else {
            cnt[nuc_index(c)] += 1;
        }
    }

    // Read lengths never approach i32::MAX in practice; matches T1K's own
    // `int i` loop counter in the ported C++ (`FastqExtractor.cpp:92`).
    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
    let len = seq.len() as i32;
    if cnt[0] >= len / 2
        || cnt[1] >= len / 2
        || cnt[2] >= len / 2
        || cnt[3] >= len / 2
        || cnt[4] >= len / 10
    {
        return true;
    }

    let low_cnt = cnt[..4].iter().filter(|&&c| c <= 2).count();
    low_cnt >= 2
}

/// Maps `A`/`C`/`G`/`T` to their `nucToNum` index (0-3). Panics on any other
/// byte -- see [`is_low_complexity`]'s doc comment.
fn nuc_index(c: u8) -> usize {
    match c {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        other => panic!(
            "is_low_complexity: unsupported base {other:?} (only A/C/G/T/N are supported; \
             see this module's doc comment)"
        ),
    }
}

/// A single parsed FASTA record: `id` is the header token up to the first
/// whitespace (matching kseq's `name` field, which stock T1K reads via
/// `ReadFiles::Next` -- `ReadFiles.hpp:183`), `seq` is every subsequent line
/// up to (not including) the next `>` header or end of file, concatenated
/// with no separators (matching multi-line FASTA wrapping).
struct FastaRecord {
    id: String,
    seq: String,
}

/// Minimal FASTA parser sufficient for T1K-style reference files (plain
/// text, uppercase `A`/`C`/`G`/`T`/`N`, optionally multi-line-wrapped
/// sequences). Not a general-purpose FASTA/FASTQ reader -- T1K's own
/// `ReadFiles`/`kseq.h` (gzip-transparent, FASTA-and-FASTQ) is considerably
/// more general, but reference FASTAs in this pipeline are always plain-text
/// FASTA, so this suffices without adding a dependency.
fn parse_fasta(text: &str) -> Vec<FastaRecord> {
    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                records.push(FastaRecord { id, seq: std::mem::take(&mut current_seq) });
            }
            current_id = Some(rest.split_whitespace().next().unwrap_or("").to_string());
        } else {
            current_seq.push_str(line.trim_end());
        }
    }
    if let Some(id) = current_id {
        records.push(FastaRecord { id, seq: current_seq });
    }

    records
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_fasta(records: &[(&str, &str)]) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().expect("tempfile");
        for (id, seq) in records {
            writeln!(f, ">{id}").unwrap();
            writeln!(f, "{seq}").unwrap();
        }
        f.flush().unwrap();
        f
    }

    #[test]
    fn parse_fasta_splits_id_at_first_whitespace_and_concatenates_multiline_seq() {
        let text = ">seq1 some comment\nACGT\nACGT\n>seq2\nTTTT\n";
        let records = parse_fasta(text);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].seq, "ACGTACGT");
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].seq, "TTTT");
    }

    #[test]
    fn is_low_complexity_flags_homopolymer() {
        assert!(is_low_complexity(b"AAAAAAAAAAAAAAAAAAAA"));
    }

    #[test]
    fn is_low_complexity_flags_dinucleotide_repeat() {
        // "ATATAT..." has A and T each at ~50% (>= len/2), and C/G at 0
        // (<=2), so two independent triggers both fire on this input.
        assert!(is_low_complexity(b"ATATATATATATATATATAT"));
    }

    #[test]
    fn is_low_complexity_flags_mostly_n() {
        // 20 bases, 10 Ns: cnt[4]=10 >= len/10=2.
        assert!(is_low_complexity(b"NNNNNNNNNNACGTACGTAC"));
    }

    #[test]
    fn is_low_complexity_passes_balanced_sequence() {
        // Roughly balanced composition, no long homopolymer run collapsing
        // into a single dominant base, no more than one base at <=2 count.
        assert!(!is_low_complexity(b"ACGTACGTACGTACGTACGTACGTACGTACGT"));
    }

    #[test]
    fn from_reference_fasta_indexes_every_sequence() {
        let f =
            write_fasta(&[("a", "ACGTACGTACGTACGTACGTACGT"), ("b", "TTTTGGGGCCCCAAAATTTTGGGG")]);
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();
        assert_eq!(filter.seq_count(), 2);
        assert_eq!(filter.kmer_length(), 9);
        assert_eq!(filter.seq_name(0), "a");
        assert_eq!(filter.seq_name(1), "b");
    }

    #[test]
    fn is_good_candidate_true_for_exact_reference_substring() {
        let reference = "ACGTACGTACGTACGTACGTACGTGGATTACAGATTACAGATTACAGATTACAG\
                          CCCTGACGTGTGACGTGTGACGTGTGACGTGTGACGTGT";
        let f = write_fasta(&[("only", reference)]);
        // hitLenRequired defaults to 31; at kmer_length=9, this needs
        // max >= ceil(31/9) = 4 hits in a single bucket, easily satisfied
        // by a long exact substring.
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();

        // A 40bp exact substring of the loaded (and therefore indexed)
        // reference sequence.
        let read = &reference.as_bytes()[10..50];
        assert!(filter.is_good_candidate(read), "exact reference substring should be a candidate");
    }

    #[test]
    fn is_good_candidate_false_for_unrelated_sequence() {
        let reference = "ACGTACGTACGTACGTACGTACGTGGATTACAGATTACAGATTACAGATTACAG\
                          CCCTGACGTGTGACGTGTGACGTGTGACGTGTGACGTGT";
        let f = write_fasta(&[("only", reference)]);
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();

        // 40 bases chosen to avoid containing any 9-mer substring of
        // `reference` (heavy G/C alternation absent from the reference).
        let read = b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        assert!(!filter.is_good_candidate(read));
    }

    #[test]
    fn is_good_candidate_false_for_short_read() {
        let reference = "ACGTACGTACGTACGTACGTACGTGGATTACAGATTACA";
        let f = write_fasta(&[("only", reference)]);
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();
        assert!(!filter.is_good_candidate(b"ACGTACGT")); // len 8 < kmer_length 9
    }

    /// Builds a synthetic reference of 150 distinct 5bp `"AAAAA"` sequences.
    /// Building the k=4 index from a homopolymer sequence exactly `k` bases
    /// long would insert NOTHING (`KmerIndex::build_index_from_read`'s
    /// dedup drops the only window, since it coincidentally matches the
    /// all-zero initial `prevKmerCode` sentinel -- see `kmer_index`'s module
    /// docs); one extra base (5bp, not 4bp) gives a second window that IS
    /// unconditionally inserted (the `i == kl` branch), so each of the 150
    /// sequences contributes exactly one `"AAAA"` index entry, for exactly
    /// 150 entries total -- comfortably over the `size >= 100`
    /// high-frequency threshold this test exercises.
    fn high_frequency_reference() -> tempfile::NamedTempFile {
        let records: Vec<(String, &'static str)> =
            (0..150).map(|i| (format!("seq{i}"), "AAAAA")).collect();
        let record_refs: Vec<(&str, &str)> =
            records.iter().map(|(id, seq)| (id.as_str(), *seq)).collect();
        write_fasta(&record_refs)
    }

    #[test]
    fn high_frequency_kmer_skipped_when_not_at_read_boundary() {
        // Regression test for the SeqSet.hpp:1109/1176 `size >= 100`
        // high-frequency-kmer skip branch (GetHitsFromRead), which is
        // otherwise not naturally exercised by the small differential-test
        // fixture (kir_rna_seq.fa has too few sequences for any k-mer to
        // reach 100 index entries).
        let f = high_frequency_reference();
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 4).unwrap();

        // Sanity: the "AAAA" k-mer really does have >= 100 index entries.
        let mut kc = KmerCode::new(4);
        for &c in b"AAAA" {
            kc.append(c);
        }
        assert!(filter.seq_index.search(&kc).len() >= 100);

        // "AAAA" appears exactly once, at read positions 5-8 (0-based),
        // which is neither the read's first window (i == kmer_length - 1 ==
        // 3) nor its last window (i == len - 1 == 13) -- so it is eligible
        // for the high-frequency skip, and (with skip_cnt starting at 0 <
        // skip_limit = kmer_length/2 = 2) gets skipped. Every other 4-mer in
        // this read is unique (absent from a reference containing only
        // "AAAAA" sequences), so the ONLY possible source of hits is this
        // one "AAAA" window -- if it is correctly skipped, total hits must
        // be exactly 0.
        let read = b"GATTCAAAAGATTC";
        assert_eq!(read.len(), 14);
        let mut scratch = Scratch::default();
        filter.get_hits_from_read(read, &mut scratch);
        assert_eq!(
            scratch.hits.len(),
            0,
            "the sole high-frequency window is mid-read and must be skipped entirely"
        );
    }

    #[test]
    fn high_frequency_kmer_not_skipped_at_first_window() {
        // Same high-frequency reference as above, but this time "AAAA" is
        // the read's very FIRST window (i == kmer_length - 1), which is
        // exempt from the high-frequency skip regardless of `skip_cnt`
        // (SeqSet.hpp:1109: `i != kmerLength - 1 && ...`). Every hit
        // recorded must therefore come from that first window alone (all
        // other windows are unique, absent from the reference).
        let f = high_frequency_reference();
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 4).unwrap();

        let read = b"AAAAGATTC";
        let mut scratch = Scratch::default();
        filter.get_hits_from_read(read, &mut scratch);

        let forward_hits = scratch.hits.iter().filter(|h| h.strand == 1).count();
        assert_eq!(
            forward_hits, 150,
            "the first-window exemption must still record all 150 reference hits"
        );
    }

    #[test]
    fn set_hit_len_required_and_ref_seq_similarity_round_trip() {
        let f = write_fasta(&[("a", "ACGTACGTACGTACGTACGTACGT")]);
        let mut filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();
        assert_eq!(filter.hit_len_required(), DEFAULT_HIT_LEN_REQUIRED);

        filter.set_hit_len_required(42);
        assert_eq!(filter.hit_len_required(), 42);

        // ref_seq_similarity has no public getter; exercised indirectly via
        // is_good_candidate in the differential/extract-module tests. Here
        // we only confirm the setter does not panic and the filter remains
        // otherwise usable.
        filter.set_ref_seq_similarity(0.5);
        let read = b"ACGTACGTACGTACGTACGTACGT";
        let _ = filter.is_good_candidate(read);
    }

    #[test]
    fn infer_kmer_length_matches_seqset_formula_for_known_totals() {
        // InferKmerLength: ret = 0; while (total) { ret+=1; total/=4; } ret+=1;
        // total=0   -> loop never runs -> ret = 0 + 1 = 1
        // total=1   -> 1 iter (1/4=0)   -> ret = 1 + 1 = 2
        // total=8781 (kir_rna_seq.fa's real total, per ref_kmer_filter's own
        //   differential fixture) -> known real-world value: 8.
        let empty = tempfile::NamedTempFile::new().unwrap();
        let filter_empty = RefKmerFilter::from_reference_fasta(empty.path(), 9).unwrap();
        assert_eq!(filter_empty.seq_count(), 0);
        assert_eq!(filter_empty.infer_kmer_length(), 1);

        let f_one = write_fasta(&[("a", "A")]);
        let filter_one = RefKmerFilter::from_reference_fasta(f_one.path(), 9).unwrap();
        assert_eq!(filter_one.infer_kmer_length(), 2);

        // total=8781: 8781/4=2195, /4=548, /4=137, /4=34, /4=8, /4=2, /4=0 ->
        // 7 iterations -> ret = 7 + 1 = 8.
        let seq_8781 = "A".repeat(8781);
        let f_8781 = write_fasta(&[("a", &seq_8781)]);
        let filter_8781 = RefKmerFilter::from_reference_fasta(f_8781.path(), 9).unwrap();
        assert_eq!(filter_8781.infer_kmer_length(), 8);
    }

    #[test]
    fn update_kmer_length_rebuilds_index_at_new_k() {
        let reference = "ACGTACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTG";
        let f = write_fasta(&[("only", reference)]);
        let mut filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();
        assert_eq!(filter.kmer_length(), 9);

        // Before update: a 40bp exact substring is a candidate at k=9
        // (matches from_reference_fasta's own is_good_candidate test).
        let read40 = &reference.as_bytes()[0..40];
        assert!(filter.is_good_candidate(read40));

        filter.update_kmer_length(15);
        assert_eq!(filter.kmer_length(), 15);
        // Index rebuilt at k=15: an exact-substring read that is still
        // comfortably long relative to the new k should still be a
        // candidate (hitLenRequired defaults to 31, so needs >= ceil(31/15)
        // = 3 hits in one bucket; a 40bp exact substring at k=15 has 26
        // overlapping windows, easily satisfying this).
        let read = reference.as_bytes();
        assert!(filter.is_good_candidate(read));

        // A read shorter than the new k=15 can no longer produce any k-mer
        // window at all (previously fine at k=9).
        let too_short_for_new_k = &reference.as_bytes()[0..10];
        assert!(!filter.is_good_candidate(too_short_for_new_k));
    }

    // ---- hit_less_than (_hit::operator<) -----------------------------------

    fn h(idx: u32, offset: u32, read_offset: i32, strand: i8) -> Hit {
        Hit { idx, offset, read_offset, strand, repeats: 1 }
    }

    #[test]
    fn hit_less_than_orders_by_strand_first() {
        let a = h(5, 5, 5, -1);
        let b = h(0, 0, 0, 1);
        assert!(hit_less_than(&a, &b), "strand -1 < strand 1 regardless of other fields");
        assert!(!hit_less_than(&b, &a));
    }

    #[test]
    fn hit_less_than_ties_on_strand_fall_to_idx() {
        let a = h(1, 0, 0, 1);
        let b = h(2, 0, 0, 1);
        assert!(hit_less_than(&a, &b));
        assert!(!hit_less_than(&b, &a));
    }

    #[test]
    fn hit_less_than_ties_on_strand_and_idx_fall_to_read_offset() {
        let a = h(0, 0, 3, 1);
        let b = h(0, 0, 7, 1);
        assert!(hit_less_than(&a, &b));
        assert!(!hit_less_than(&b, &a));
    }

    #[test]
    fn hit_less_than_final_tiebreak_is_offset() {
        let a = h(0, 3, 0, 1);
        let b = h(0, 7, 0, 1);
        assert!(hit_less_than(&a, &b));
        assert!(!hit_less_than(&b, &a));
    }

    #[test]
    fn hit_less_than_bit_identical_hits_are_neither_less() {
        let a = h(1, 2, 3, 1);
        let b = a;
        assert!(!hit_less_than(&a, &b));
        assert!(!hit_less_than(&b, &a));
    }

    // ---- pack_hit_key (radix key monotonicity) -----------------------------

    #[test]
    fn pack_hit_key_matches_hit_less_than_order() {
        // For a diverse set of hits spanning both strands, small/large idx,
        // negative/zero/positive read_offset, and varying offset (including
        // extreme u32/i32 boundary values), the unsigned u128 ordering of
        // `pack_hit_key` must be monotonic in `hit_less_than`:
        //   pack(a) < pack(b)  <=>  hit_less_than(a, b)
        //   pack(a) == pack(b) <=>  neither less (bit-identical fields).
        // This is the invariant that makes a radix sort on `pack_hit_key`
        // byte-identical to `sort_unstable_by(hit_less_than)`.
        let hits = vec![
            h(0, 0, 0, -1),
            h(0, 0, 0, 1),
            h(0, 0, -5, 1),
            h(0, 0, i32::MIN, 1),
            h(0, 0, i32::MAX, 1),
            h(0, u32::MAX, 0, 1),
            h(u32::MAX, 0, 0, 1),
            h(3, 7, 3, -1),
            h(3, 7, 3, 1),
            h(3, 8, 3, 1),
            h(3, 7, 4, 1),
            h(3, 7, 3, 1), // duplicate of an earlier hit -> must tie
            h(1, u32::MAX, i32::MIN, 0),
        ];
        for a in &hits {
            for b in &hits {
                let ka = pack_hit_key(a);
                let kb = pack_hit_key(b);
                assert_eq!(
                    ka < kb,
                    hit_less_than(a, b),
                    "pack ordering must match hit_less_than for {a:?} vs {b:?}",
                );
                assert_eq!(
                    ka == kb,
                    !hit_less_than(a, b) && !hit_less_than(b, a),
                    "pack equality must match field-identity for {a:?} vs {b:?}",
                );
            }
        }

        // A full sort by pack_hit_key must equal a sort by hit_less_than.
        let mut by_pack = hits.clone();
        radsort::sort_by_key(&mut by_pack, pack_hit_key);
        let mut by_cmp = hits;
        by_cmp.sort_by(|a, b| {
            if hit_less_than(a, b) {
                std::cmp::Ordering::Less
            } else if hit_less_than(b, a) {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        });
        assert_eq!(by_pack, by_cmp);
    }

    // ---- sort_hits (SeqSet::SortHits) --------------------------------------

    #[test]
    fn sort_hits_uses_std_sort_path_when_below_bucket_threshold() {
        // hits.len() = 3, seq_count = 10 -> 3 > 2*10 is false, so this takes
        // the std::sort (hit_less_than) branch regardless of
        // already_read_order.
        let mut hits = vec![h(0, 0, 5, 1), h(0, 0, 1, 1), h(1, 0, 0, -1)];
        sort_hits(&mut hits, true, 10);
        // Expected order: strand -1 first, then strand 1 ascending by idx
        // then read_offset.
        assert_eq!(hits, vec![h(1, 0, 0, -1), h(0, 0, 1, 1), h(0, 0, 5, 1)]);
    }

    #[test]
    fn sort_hits_bucket_path_preserves_within_bucket_order() {
        // hits.len() = 7, seq_count = 2 -> 7 > 2*2 = 4, and
        // already_read_order = true -> bucket-sort path. Two sequences
        // (idx 0, 1), mixed strands; within each (tag, idx) bucket the
        // ORIGINAL encounter order must be preserved (NOT re-sorted by
        // read_offset/offset -- unlike the std::sort branch above).
        let mut hits = vec![
            h(1, 0, 30, 1),  // tag=1, idx=1
            h(0, 0, 20, -1), // tag=0, idx=0
            h(1, 0, 10, 1),  // tag=1, idx=1 (comes AFTER the read_offset=30 one above)
            h(0, 0, 5, -1),  // tag=0, idx=0
            h(0, 0, 99, 1),  // tag=1, idx=0
            h(1, 0, 1, -1),  // tag=0, idx=1
            h(0, 0, 2, 1),   // tag=1, idx=0
        ];
        let original = hits.clone();
        sort_hits(&mut hits, true, 2);

        // Bucket order: tag=0,idx=0 then tag=0,idx=1 then tag=1,idx=0 then
        // tag=1,idx=1. Within each bucket, original relative order.
        let expected = vec![
            original[1], // tag=0,idx=0 (read_offset=20)
            original[3], // tag=0,idx=0 (read_offset=5)  <- NOT resorted ahead of 20
            original[5], // tag=0,idx=1 (read_offset=1)
            original[4], // tag=1,idx=0 (read_offset=99)
            original[6], // tag=1,idx=0 (read_offset=2)
            original[0], // tag=1,idx=1 (read_offset=30)
            original[2], // tag=1,idx=1 (read_offset=10) <- stays AFTER the 30 one
        ];
        assert_eq!(hits, expected);
    }

    #[test]
    fn sort_hits_falls_back_to_std_sort_when_not_already_read_order() {
        // Even with hits.len() > 2*seq_count, already_read_order=false must
        // take the std::sort (hit_less_than) branch, not bucket sort.
        let mut hits = vec![h(0, 0, 30, 1), h(0, 0, 10, 1), h(0, 0, 20, 1)];
        sort_hits(&mut hits, false, 1); // 3 > 2*1 = 2, but already_read_order=false
        assert_eq!(hits, vec![h(0, 0, 10, 1), h(0, 0, 20, 1), h(0, 0, 30, 1)]);
    }

    // ---- is_overlap_low_complex (SeqSet::IsOverlapLowComplex) -------------

    fn overlap_span(read_start: i32, read_end: i32) -> overlap::Overlap {
        overlap::Overlap {
            seq_idx: 0,
            read_start,
            read_end,
            seq_start: 0,
            seq_end: read_end - read_start,
            strand: 1,
            match_cnt: 0,
            similarity: 0.0,
        }
    }

    #[test]
    fn is_overlap_low_complex_flags_homopolymer_span() {
        // All-A over a 30bp span: cnt=[30,0,0,0], lowCnt=3 (C/G/T each 0),
        // lowTotalCnt=0, 0*7 >= 30 is false, lowCnt>=2 -> true.
        let r = vec![b'A'; 30];
        let o = overlap_span(0, 29);
        assert!(is_overlap_low_complex(&r, &o));
    }

    #[test]
    fn is_overlap_low_complex_passes_balanced_span() {
        // Balanced ACGT repeat: each base ~equally represented, well above
        // the <=2 threshold for a 32bp span (8 of each).
        let r = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let o = overlap_span(0, 31);
        assert!(!is_overlap_low_complex(&r, &o));
    }

    #[test]
    fn is_overlap_low_complex_low_total_below_one_seventh_still_passes() {
        // 21bp span, mostly A with 2 C, 0 G, 0 T: cnt=[19,2,0,0].
        // lowCnt=3 (C=2,G=0,T=0 all <=2), lowTotalCnt=2+0+0=2.
        // 2*7=14 >= 21? No -> falls through to lowCnt>=2 -> true.
        // (This case is here mainly to hand-verify the arithmetic; see the
        // next test for a case where the 1/7 threshold actually flips it.)
        let mut r = vec![b'A'; 21];
        r[0] = b'C';
        r[1] = b'C';
        let o = overlap_span(0, 20);
        assert!(is_overlap_low_complex(&r, &o));
    }

    #[test]
    fn is_overlap_low_complex_ignores_n_positions() {
        // N positions are skipped entirely (not counted into any base's
        // tally): an all-N span has cnt=[0,0,0,0], lowCnt=4, lowTotalCnt=0,
        // 0*7 >= len is false (len>0) -> lowCnt>=2 -> true.
        let r = vec![b'N'; 10];
        let o = overlap_span(0, 9);
        assert!(is_overlap_low_complex(&r, &o));
    }

    // ---- get_overlaps_from_read: hand-computed exact-substring overlap ----

    #[test]
    fn get_overlaps_from_read_hand_computed_exact_substring() {
        // A single reference sequence; the read is an exact 60bp substring
        // of it. Exact substring -> a k-mer hit at EVERY read_offset in
        // 0..=(len-kmerLength), all sharing the same coord diff -> single
        // colinear hit chain -> no AlignAlgo gap branch is ever entered
        // (every consecutive hit pair has hitCoords[j-1].a + kmerLength - 1
        // >= hitCoords[j].a, i.e. the direct "matchCnt += 2*(a-aPrev)"
        // branch). Starting from `matchCnt = 2*kmerLength`
        // (SeqSet.hpp:1697) and accumulating `2*(a - aPrev)` for every
        // consecutive pair telescopes to `matchCnt = 2*kmerLength +
        // 2*((len-kmerLength) - 0) = 2*len` for a dense, single-step-offset
        // exact-match chain covering the whole read -- hand-verifiable
        // without any AlignAlgo call. `seqSpan = seqEnd-seqStart+1 = len`
        // and `readSpan = readEnd-readStart+1 = len` (both `+1`, NOT
        // `len-1` -- an inclusive span over `len` positions has length
        // `len`), so `similarity = matchCnt/(seqSpan+readSpan) = 2*len /
        // (2*len) = 1.0` exactly for this dense/exact-match case (a genuine
        // coincidence of this specific scenario -- the `similarity == 1`
        // check in the C++, `SeqSet.hpp:1677,1838`, is itself just a
        // SENTINEL flag meaning "no break-triggering condition was hit",
        // not a claim that the final ratio is always 1.0 in general).
        let reference = "ACGTACGTACGTACGTACGTACGTGGATTACAGATTACAGATTACAGATTACAG\
                          CCCTGACGTGTGACGTGTGACGTGTGACGTGTGACGTGT";
        let f = write_fasta(&[("only", reference)]);
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();

        let start = 10;
        let len = 60;
        let read = &reference.as_bytes()[start..start + len];
        let mut scratch = Scratch::default();
        let overlaps = filter
            .get_overlaps_from_read(read, &mut scratch)
            .expect("read is long enough for a k-mer window");

        assert_eq!(overlaps.len(), 1, "expected exactly one surviving overlap: {overlaps:?}");
        let o = &overlaps[0];
        assert_eq!(o.seq_idx, 0);
        assert_eq!(o.strand, 1);
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let len_i32 = len as i32;
        assert_eq!(o.read_start, 0);
        assert_eq!(o.read_end, len_i32 - 1);
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let start_i32 = start as i32;
        assert_eq!(o.seq_start, start_i32);
        assert_eq!(o.seq_end, start_i32 + len_i32 - 1);

        // matchCnt = 2*len (derivation above); seqSpan == readSpan == len.
        let expected_match_cnt = 2 * len_i32;
        assert_eq!(o.match_cnt, expected_match_cnt);
        let seq_span = o.seq_end - o.seq_start + 1;
        let read_span = o.read_end - o.read_start + 1;
        assert_eq!(seq_span, len_i32);
        assert_eq!(read_span, len_i32);
        let expected_similarity = f64::from(expected_match_cnt) / f64::from(seq_span + read_span);
        assert!(
            (expected_similarity - 1.0).abs() < f64::EPSILON,
            "sanity: hand-derivation should reduce to exactly 1.0"
        );
        assert!(
            (o.similarity - expected_similarity).abs() < f64::EPSILON,
            "expected similarity == {expected_similarity} exactly, got {}",
            o.similarity
        );
    }

    #[test]
    fn get_overlaps_from_read_returns_none_for_short_read() {
        let reference = "ACGTACGTACGTACGTACGTACGTGGATTACAGATTACA";
        let f = write_fasta(&[("only", reference)]);
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();
        let mut scratch = Scratch::default();
        assert!(filter.get_overlaps_from_read(b"ACGTACGT", &mut scratch).is_none());
    }
}
