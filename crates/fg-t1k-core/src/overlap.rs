//! LIS-based overlap/hit-chaining, ported from T1K's `SeqSet::GetOverlapsFromHits`
//! (`vendor/t1k/SeqSet.hpp:1232-1556`) and its helpers
//! `LongestIncreasingSubsequence`/`BinarySearch_LIS` (`SeqSet.hpp:352-436,327-348`),
//! `GetTotalHitLengthOnRead`/`GetTotalHitLengthOnSeq` (`SeqSet.hpp:1032-1069`), and
//! comparators `CompSortPairBInc`/`CompSortHitCoordDiff` (`SeqSet.hpp:233-239,266-274`).
//!
//! This is the "second gate" of `SeqSet::HasHitInSet` (`SeqSet.hpp:1966-1983`):
//! after the bucket-count gate ([`crate::ref_kmer_filter::RefKmerFilter::has_hit_in_set`]'s
//! Task-3.1 slice) picks a winning `(strand, seqIdx)` bucket, `GetOverlapsFromHits`
//! chains that bucket's k-mer hits into colinear runs via patience-sorting LIS,
//! and `HasHitInSet` accepts the read if any resulting overlap's implied
//! mismatch count is within a similarity-derived threshold.
//!
//! # No `AlignAlgo`/Smith-Waterman in this path
//!
//! `GetOverlapsFromHits` (as called from `HasHitInSet`, i.e. this module)
//! never invokes `AlignAlgo`/banded alignment -- `_overlap::matchCnt` is set
//! purely from the LIS-chained hit length (`matchCnt = 2 * hitLen`,
//! `SeqSet.hpp:1531`) and `similarity` is left at `0` (`SeqSet.hpp:1532`).
//! Real alignment-based `matchCnt`/`similarity` computation happens elsewhere
//! in stock T1K (`ExtendOverlap` and friends) and is out of scope here (a
//! later phase, unrelated to this gate).
//!
//! # Ported generally, not overfit to `filter == 0`
//!
//! `GetOverlapsFromHits` takes a `filter` parameter (`0` or `1`); only
//! `filter == 0` is reachable from `HasHitInSet`, but this port reproduces
//! BOTH branches faithfully (including the `filter == 1`
//! repeatability-based `novelMinHitRequired` recomputation and
//! `removeOnlyRepeats` filtering, `SeqSet.hpp:1256-1299,1320-1336,1408-1424`)
//! so a future caller (e.g. a `GetOverlapsFromRead`-equivalent port) can reuse
//! this function unmodified.
//!
//! # Every reference sequence has `isRef == true`
//!
//! [`crate::ref_kmer_filter::RefKmerFilter`] only ever loads sequences via
//! `InputRefFa`-equivalent logic, which stock always marks `isRef = true`
//! (`SeqSet.hpp:883,911`). This port therefore takes `is_ref: bool` as a
//! per-call parameter (mirroring what a caller would look up from
//! `seqs[idx].isRef`) rather than modeling a `_seqWrapper`-equivalent struct;
//! every caller in this codebase currently passes `true`.

/// Ported from `_pair` (`defs.h:18-21`): a generic `(a, b)` integer pair.
/// Used here for `(readOffset, seqOffset)` coordinate pairs (`concordantHitCoord`,
/// `hitCoordLIS`, `_overlap::hitCoords`-equivalent lists).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Pair {
    pub a: i32,
    pub b: i32,
}

/// Ported from `_triple` (`defs.h:40-43`): a generic `(a, b, c)` integer
/// triple. Used here for `hitCoordDiff` entries: `a` = `readOffset`, `b` =
/// `seqOffset`, `c` = `readOffset - seqOffset` (the coordinate diff used to
/// find concordant, colinear hit runs).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Triple {
    pub a: i32,
    pub b: i32,
    pub c: i32,
}

/// A single k-mer hit, ported from `_hit` (`SeqSet.hpp:66-87`). Field-identical
/// to [`crate::ref_kmer_filter::Hit`] (both carry `repeats`, `_hit::repeats`,
/// `SeqSet.hpp:72`); here `repeats` is consumed only by
/// [`get_overlaps_from_hits`]'s `filter == 1` branches. This is a standalone type (rather than
/// reusing `ref_kmer_filter::Hit` directly) so this module has no dependency
/// on `ref_kmer_filter`'s internals; callers convert at the boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct OverlapHit {
    /// The reference sequence index this hit maps to (`_indexInfo::idx`).
    pub idx: u32,
    /// The offset within that reference sequence (`_indexInfo::offset`).
    pub offset: u32,
    /// The offset within the (possibly reverse-complemented) read where this
    /// k-mer window starts (`_hit::readOffset`).
    pub read_offset: i32,
    /// `1` for a forward-strand hit, `-1` for a reverse-complement-strand hit
    /// (`_hit::strand`).
    pub strand: i8,
    /// How many times this hit's k-mer occurs across the whole index
    /// (`_hit::repeats`; only consumed by the `filter == 1` branches).
    pub repeats: i32,
}

/// A chained/colinear overlap between a read and a reference sequence,
/// ported from `_overlap` (`SeqSet.hpp:89-115`). Only the fields this port's
/// `AlignAlgo`-free `GetOverlapsFromHits` path actually populates are
/// included; `leftClip`/`rightClip`/`relaxedMatchCnt`/`align` are omitted
/// (stock's `GetOverlapsFromHits` never sets them -- they default to
/// whatever `_overlap`'s caller happens to leave, and no caller in the
/// `HasHitInSet` path reads them).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Overlap {
    /// `_overlap::seqIdx`.
    pub seq_idx: u32,
    /// `_overlap::readStart`. When `strand == -1`, this is the offset within
    /// the reverse-complemented read (matching `_hit::readOffset`'s own
    /// strand-dependent meaning).
    pub read_start: i32,
    /// `_overlap::readEnd`.
    pub read_end: i32,
    /// `_overlap::seqStart`.
    pub seq_start: i32,
    /// `_overlap::seqEnd`.
    pub seq_end: i32,
    /// `_overlap::strand`.
    pub strand: i8,
    /// `_overlap::matchCnt` -- the number of matched bases, counted TWICE
    /// (`= 2 * hitLen`, `SeqSet.hpp:1531`). NOT actually derived from
    /// alignment in this path; see module docs.
    pub match_cnt: i32,
    /// `_overlap::similarity`. Always `0.0` on this `AlignAlgo`-free path
    /// (`SeqSet.hpp:1532`).
    pub similarity: f64,
}

/// Ported from `SeqSet::CompSortPairBInc` (`SeqSet.hpp:233-239`): ascending
/// by `b`, ties broken ascending by `a`. This is a strict total order with no
/// ties possible for distinct `(a, b)` pairs, and for genuinely equal pairs
/// any total order agrees -- so this comparator poses no `std::sort`
/// stability hazard (see module docs on sort-tie handling below for the
/// comparator that DOES need care).
// By-reference parameters (rather than by-value, despite `Pair` being
// `Copy`/8 bytes) are required by `Vec::sort_unstable_by`'s `FnMut(&T, &T)`
// signature at this function's only call site below.
#[allow(clippy::trivially_copy_pass_by_ref)]
fn comp_sort_pair_b_inc(p1: &Pair, p2: &Pair) -> std::cmp::Ordering {
    p1.b.cmp(&p2.b).then_with(|| p1.a.cmp(&p2.a))
}

/// Ported from `SeqSet::CompSortHitCoordDiff` (`SeqSet.hpp:266-274`):
/// ascending by `c` (the coordinate diff), ties broken ascending by `b`, then
/// ascending by `a`. Like [`comp_sort_pair_b_inc`], this compares every field
/// of the key down to `a`, so two entries only compare equal if they are
/// bit-for-bit identical `(a, b, c)` triples -- a strict total order with no
/// observable ties. `hitCoordDiff` entries are built one per input hit
/// (`SeqSet.hpp:1339-1346`) and `a`/`b` come directly from that hit's
/// `(readOffset, seqOffset)`, so two entries CAN be genuinely identical
/// (e.g. two hits at the same read/seq coordinate from duplicate index
/// entries) -- in that case `std::sort`'s tie-order is unobservable anyway
/// (swapping two bit-identical elements changes nothing downstream), so
/// `sort_unstable_by` is safe here regardless of libstdc++'s introsort tie
/// behavior.
fn comp_sort_hit_coord_diff(a: &Triple, b: &Triple) -> std::cmp::Ordering {
    a.c.cmp(&b.c).then_with(|| a.b.cmp(&b.b)).then_with(|| a.a.cmp(&b.a))
}

/// Ported from `SeqSet::BinarySearch_LIS` (`SeqSet.hpp:327-348`): returns the
/// index (into `top`) of the last entry whose `hits[top[m]].a <= val_a`, or
/// `-1` if no such entry exists (all entries have `.a > val_a`). Reproduces
/// the exact `(l + r) / 2` midpoint and early-return-on-exact-match
/// semantics of the C++.
fn binary_search_lis(top: &[usize], size: usize, val_a: i32, hits: &[Pair]) -> i32 {
    let mut l: i32 = 0;
    let mut r: i32 = i32::try_from(size).unwrap_or(i32::MAX) - 1;
    while l <= r {
        let m = (l + r) / 2;
        #[allow(clippy::cast_sign_loss)]
        let hit_a = hits[top[m as usize]].a;
        match val_a.cmp(&hit_a) {
            std::cmp::Ordering::Equal => return m,
            std::cmp::Ordering::Less => r = m - 1,
            std::cmp::Ordering::Greater => l = m + 1,
        }
    }
    l - 1
}

/// Ported from `SeqSet::LongestIncreasingSubsequence` (`SeqSet.hpp:352-436`):
/// the O(n log n) patience-sorting LIS over `hits` (ascending on `.a`,
/// "biased towards left" per the C++'s own comment), followed by a
/// backtrace, an unconditional `Reverse()` of the reconstructed sequence, and
/// then a de-duplication pass that keeps only the FIRST element of any run of
/// consecutive elements sharing the same `.b` (`SeqSet.hpp:417-429`).
/// Returns the LIS (appended to `lis`, which is cleared first) -- the return
/// value is `lis.len()` after the dedup pass, matching the C++'s returned
/// `ret`.
///
/// `hits` must be non-empty (mirrors the C++, which is only ever called on a
/// non-empty `concordantHitCoord`, `SeqSet.hpp:1463`, and would read
/// out-of-bounds at `record[0] = 0`/`top[0] = 0` on an empty input).
///
/// # Panics
///
/// Panics if `hits` is empty.
fn longest_increasing_subsequence(hits: &[Pair], lis: &mut Vec<Pair>) -> usize {
    lis.clear();
    let size = hits.len();
    assert!(size > 0, "longest_increasing_subsequence: hits must be non-empty");

    // `record`: SeqSet.hpp:361-373. The C++ builds `record[0] = 0` then
    // `record[i] = i` for `i` in `1..size` (the commented-out
    // `if (hits[i].b == hits[i-1].b) continue;` dedup is dead code -- never
    // taken since it's commented out -- so `record` is simply `0..size`
    // unconditionally). `rcnt` therefore always equals `size`.
    let rcnt = size;

    // `top[k]`: index (into `hits`) of the hit with the smallest `.a` among
    // all increasing subsequences of length `k + 1` found so far.
    // `link[i]`: predecessor index (into `hits`) for backtrace, or `-1`.
    let mut top: Vec<usize> = vec![0; size];
    let mut link: Vec<i32> = vec![-1; size];

    top[0] = 0;
    link[0] = -1;
    let mut ret: usize = 1;

    for i in 1..rcnt {
        let record_i = i; // record[i] == i, see comment above.
        let tag: i32 = if hits[top[ret - 1]].a <= hits[record_i].a {
            i32::try_from(ret - 1).unwrap_or(i32::MAX)
        } else {
            binary_search_lis(&top, ret, hits[record_i].a, hits)
        };

        if tag == -1 {
            top[0] = record_i;
            link[record_i] = -1;
        } else {
            #[allow(clippy::cast_sign_loss)]
            let tag_usize = tag as usize;
            if hits[record_i].a > hits[top[tag_usize]].a {
                if tag_usize == ret - 1 {
                    top[ret] = record_i;
                    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                    let link_val = top[tag_usize] as i32;
                    link[record_i] = link_val;
                    ret += 1;
                } else if hits[record_i].a < hits[top[tag_usize + 1]].a {
                    top[tag_usize + 1] = record_i;
                    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                    let link_val = top[tag_usize] as i32;
                    link[record_i] = link_val;
                }
            }
        }
    }

    // Backtrace (SeqSet.hpp:407-413): walk `link` from `top[ret - 1]`,
    // pushing each visited hit, then reverse.
    let mut k: i32 = {
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let v = top[ret - 1] as i32;
        v
    };
    for _ in 0..ret {
        #[allow(clippy::cast_sign_loss)]
        let k_usize = k as usize;
        lis.push(hits[k_usize]);
        k = link[k_usize];
    }
    lis.reverse();

    // Remove elements with the same `.b` (SeqSet.hpp:417-429): keep index 0
    // unconditionally, then for each subsequent element, keep it only if its
    // `.b` differs from the LAST KEPT element's `.b`.
    if ret > 0 {
        let mut k = 1usize;
        for i in 1..ret {
            if lis[i].b == lis[k - 1].b {
                continue;
            }
            lis[k] = lis[i];
            k += 1;
        }
        lis.truncate(k);
        return k;
    }

    lis.len()
}

/// Ported from `SeqSet::GetTotalHitLengthOnRead` (`SeqSet.hpp:1032-1052`):
/// sums, over maximal runs of hits whose `read_offset`s are within
/// `kmer_length - 1` of the PREVIOUS hit in the run (i.e. their k-mer windows
/// overlap or abut on the read), `last.read_offset - first.read_offset +
/// kmer_length` per run. Requires `hits` to already be in ascending
/// `read_offset` order (as `finalHits`/`LIS` output always is here); no
/// sorting is performed internally (mirrors the C++, which also does not
/// sort here).
pub(crate) fn get_total_hit_length_on_read(hits: &[OverlapHit], kmer_length: i32) -> i32 {
    total_hit_length(hits, kmer_length, |h| h.read_offset)
}

/// Ported from `SeqSet::GetTotalHitLengthOnSeq` (`SeqSet.hpp:1054-1069`):
/// identical to [`get_total_hit_length_on_read`] but keyed on
/// `indexHit.offset` (the reference-sequence coordinate) instead of
/// `readOffset`.
pub(crate) fn get_total_hit_length_on_seq(hits: &[OverlapHit], kmer_length: i32) -> i32 {
    #[allow(clippy::cast_possible_wrap)]
    total_hit_length(hits, kmer_length, |h| h.offset as i32)
}

/// Shared implementation for [`get_total_hit_length_on_read`]/
/// [`get_total_hit_length_on_seq`]: both C++ functions have an identical
/// loop structure, differing only in which coordinate they key on.
fn total_hit_length(
    hits: &[OverlapHit],
    kmer_length: i32,
    coord: impl Fn(&OverlapHit) -> i32,
) -> i32 {
    let hit_size = hits.len();
    let mut ret: i32 = 0;
    let mut i = 0usize;
    while i < hit_size {
        let mut j = i + 1;
        while j < hit_size {
            if coord(&hits[j]) > coord(&hits[j - 1]) + kmer_length - 1 {
                break;
            }
            j += 1;
        }
        ret += coord(&hits[j - 1]) - coord(&hits[i]) + kmer_length;
        i = j;
    }
    ret
}

/// Ported from `SeqSet::GetOverlapsFromHits` (`SeqSet.hpp:1232-1556`) --
/// see module docs for scope notes (no `AlignAlgo`, both `filter` branches
/// ported, `is_ref` supplied per-call).
///
/// `hits` need not be pre-grouped by `(strand, idx)` -- like the C++, this
/// scans for maximal runs of consecutive same-`(strand, idx)` hits itself
/// (`SeqSet.hpp:1303-1307`). In practice (as called from `has_hit_in_set`'s
/// gate 2), `hits` already contains only a single `(strand, idx)` group (the
/// winning bucket), so this degenerates to one pass over the whole slice.
///
/// `radius`/`is_long_seq_set`/`kmer_length` mirror `SeqSet::radius`/
/// `SeqSet::isLongSeqSet`/`SeqSet::kmerLength` (constructor defaults `10`/
/// `false`, plus the caller-supplied k-mer length).
///
/// Returns the list of chained overlaps (mirrors the C++'s `overlaps` output
/// parameter, but returned by value here since this port doesn't need the
/// parallel `overlapsHitCoords` list). This port intentionally does NOT
/// return per-overlap hit-coordinate lists at all, since `HasHitInSet`'s
/// gate 2 (the only consumer in this codebase so far) never reads them -- it
/// only reads `matchCnt`. A future caller that needs them (e.g. a
/// `GetOverlapsFromRead` port) can add that alongside this function without
/// changing this function's core logic.
// Single-character loop/index variable names (`i`, `j`, `k`, `s`, `e`) are
// used deliberately throughout this function's body: they mirror
// `GetOverlapsFromHits`'s own C++ variable names one-for-one, which is
// intentional here to keep this a faithful, line-comparable port (per this
// module's doc comment) rather than a from-scratch reimplementation.
// Likewise, a couple of `for k in s..e { ... run[k] ... }`-style loops below
// are left in explicit index form (rather than converted to `run[s..e]`
// iterator chains) because they mix indexing into `run`/`hit_coord_lis`
// with indexing into OTHER same-length collections at the same position,
// and/or because the equivalent C++ is itself an explicit indexed loop this
// port is deliberately kept line-comparable to.
#[allow(clippy::too_many_lines, clippy::many_single_char_names, clippy::needless_range_loop)]
pub(crate) fn get_overlaps_from_hits(
    hits: &[OverlapHit],
    hit_len_required: i32,
    filter: i32,
    is_ref: impl Fn(u32) -> bool,
    radius: i32,
    is_long_seq_set: bool,
    kmer_length: i32,
) -> Vec<Overlap> {
    let hit_size = hits.len();
    let mut overlaps: Vec<Overlap> = Vec::new();
    if hit_size == 0 {
        return overlaps;
    }

    // `maxReadOffset`/`readOffsetUsed` (SeqSet.hpp:1242-1247): a scratch
    // array sized to the largest readOffset across ALL input hits (not just
    // the current [s,e) window), reused (and re-initialized per-entry, never
    // bulk-cleared) across the whole function -- mirrors the C++'s single
    // `new int[maxReadOffset + 1]` allocation shared by every iteration of
    // the outer/inner loops below.
    let max_read_offset = hits.iter().map(|h| h.read_offset).max().unwrap_or(-1);
    #[allow(clippy::cast_sign_loss)]
    let read_offset_used_len = (max_read_offset + 1).max(0) as usize;
    let mut read_offset_used: Vec<i32> = vec![0; read_offset_used_len];

    // `novelMinHitRequired`/`refMinHitRequired`/`removeOnlyRepeats`
    // (SeqSet.hpp:1252-1299): per-strand-tag thresholds/flags, only
    // recomputed away from their `{3, 3}`/`{false, false}` defaults when
    // `filter == 1`.
    let mut novel_min_hit_required: [i32; 2] = [3, 3];
    let ref_min_hit_required: [i32; 2] = [3, 3];
    let mut remove_only_repeats: [bool; 2] = [false, false];

    if filter == 1 {
        let mut longest_hits: [i32; 2] = [0, 0];
        let mut possible_overlap_cnt: [i32; 2] = [0, 0];

        let mut i = 0usize;
        while i < hit_size {
            let is_plus_strand = usize::from((1 + i32::from(hits[i].strand)) / 2 == 1);
            let mut j = i + 1;
            while j < hit_size {
                if hits[j].strand != hits[i].strand || hits[j].idx != hits[i].idx {
                    break;
                }
                j += 1;
            }

            if !is_ref(hits[i].idx) {
                #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                let run_len = (j - i) as i32;
                if run_len > novel_min_hit_required[is_plus_strand] {
                    possible_overlap_cnt[is_plus_strand] += 1;
                }
                if run_len > longest_hits[is_plus_strand] {
                    longest_hits[is_plus_strand] = run_len;
                }
            }

            if !remove_only_repeats[is_plus_strand] {
                let cnt = hits[i..j].iter().filter(|h| h.repeats <= 10000).count();
                if i32::try_from(cnt).unwrap_or(i32::MAX) >= novel_min_hit_required[is_plus_strand]
                {
                    remove_only_repeats[is_plus_strand] = true;
                }
            }

            i = j;
        }

        // filter based on the repeatability of overlaps (SeqSet.hpp:1287-1299).
        // Double->int truncation toward zero, computed in the exact operand
        // order the C++ uses (`longestHits[i] * 0.75` etc., THEN truncate).
        for i in 0..2 {
            if possible_overlap_cnt[i] > 100_000 {
                novel_min_hit_required[i] = trunc_i32(f64::from(longest_hits[i]) * 0.75);
            } else if possible_overlap_cnt[i] > 10_000 {
                novel_min_hit_required[i] = longest_hits[i] / 2;
            } else if possible_overlap_cnt[i] > 1_000 {
                novel_min_hit_required[i] = longest_hits[i] / 3;
            } else if possible_overlap_cnt[i] > 100 {
                novel_min_hit_required[i] = longest_hits[i] / 4;
            }
        }
    }

    let mut hit_coord_diff: Vec<Triple> = Vec::new();
    let mut concordant_hit_coord: Vec<Pair> = Vec::new();
    let mut hit_coord_lis: Vec<Pair> = Vec::new();
    let mut final_hits: Vec<OverlapHit> = Vec::new();

    let mut i = 0usize;
    while i < hit_size {
        let mut j = i + 1;
        while j < hit_size {
            if hits[j].strand != hits[i].strand || hits[j].idx != hits[i].idx {
                break;
            }
            j += 1;
        }

        let is_plus_strand = usize::from((1 + i32::from(hits[i].strand)) / 2 == 1);
        let min_hit_required = if is_ref(hits[i].idx) {
            ref_min_hit_required[is_plus_strand]
        } else {
            novel_min_hit_required[is_plus_strand]
        };

        // [i, j) holds the hits onto the same seq on the same strand.
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let run_len = (j - i) as i32;
        if run_len < min_hit_required {
            i = j;
            continue;
        }

        if remove_only_repeats[is_plus_strand] {
            let has_unique = hits[i..j].iter().any(|h| h.repeats <= 10000);
            if !has_unique {
                i = j;
                continue;
            }
        }

        hit_coord_diff.clear();
        for hit in &hits[i..j] {
            #[allow(clippy::cast_possible_wrap)]
            let offset_i32 = hit.offset as i32;
            hit_coord_diff.push(Triple {
                a: hit.read_offset,
                b: offset_i32,
                c: hit.read_offset - offset_i32,
            });
        }
        hit_coord_diff.sort_unstable_by(comp_sort_hit_coord_diff);

        // Pick the best concordant hits (SeqSet.hpp:1349-1551).
        let adjust_radius = if is_ref(hits[i].idx) { radius } else { 0 };

        let run = &hit_coord_diff[..];
        let run_len_usize = j - i;
        let mut s = 0usize;
        while s < run_len_usize {
            let mut current_coord_diff = run[s].c;
            let mut current_coord_diff_cnt: i32 = 1;
            let mut dominant_coord_diff = 0i32;
            let mut dominant_coord_diff_cnt: i32 = 0;
            #[allow(clippy::cast_sign_loss)]
            {
                read_offset_used[run[s].a as usize] = -1;
            }

            let mut e = s + 1;
            while e < run_len_usize {
                let mut diff = run[e].c - run[e - 1].c;
                if diff < 0 {
                    diff = -diff;
                }
                if diff > adjust_radius {
                    break;
                }

                if diff == 0 {
                    current_coord_diff_cnt += 1;
                } else {
                    if current_coord_diff_cnt > dominant_coord_diff_cnt {
                        dominant_coord_diff = current_coord_diff;
                        dominant_coord_diff_cnt = current_coord_diff_cnt;
                    }
                    current_coord_diff = run[e].c;
                    current_coord_diff_cnt = 1;
                }

                #[allow(clippy::cast_sign_loss)]
                {
                    read_offset_used[run[e].a as usize] = -1;
                }
                e += 1;
            }
            if current_coord_diff_cnt > dominant_coord_diff_cnt {
                dominant_coord_diff = current_coord_diff;
                // C++ SeqSet.hpp:1396 self-assigns `dominantCoordDiffCnt =
                // dominantCoordDiffCnt` here (a no-op typo in the vendored
                // source) -- deliberately not reproduced as a statement since
                // it has no observable effect (see this function's module
                // docs); `dominant_coord_diff_cnt` is not read again below.
            }

            if e - s < usize_from_i32_floor(min_hit_required)
                || i32::try_from(e - s).unwrap_or(i32::MAX) * kmer_length < hit_len_required
            {
                s = e;
                continue;
            }

            // SeqSet.hpp:1408-1424: NOTE this indexes the ORIGINAL, full
            // `hits` array (parameter-scoped indices `0..hitSize`) directly
            // with the post-hitCoordDiff-sort `s`/`e` bounds -- NOT
            // `hitCoordDiff` (whose entries were REORDERED by the
            // `CompSortHitCoordDiff` sort above) and NOT offset by the
            // current group's start `i`. This is almost certainly an
            // unintentional index-reuse bug in stock (comparing
            // sorted-hitCoordDiff-local positions against unrelated absolute
            // `hits` entries whenever `i > 0` or the sort actually
            // reordered anything), but it is EXACTLY what stock does, so it
            // is reproduced verbatim rather than "fixed". For this
            // codebase's only caller (`has_hit_in_set`'s gate 2, always
            // `filter == 0` on a single already-homogeneous `(strand, idx)`
            // bucket), `remove_only_repeats` is never set (`filter == 1`
            // only) and `i` is always `0`, so this quirk is provably inert
            // on the reachable path -- but it is still ported generally
            // (matching stock bit-for-bit) rather than special-cased away,
            // per this module's "ported generally" scope note.
            if remove_only_repeats[is_plus_strand] {
                let hi = s.min(hit_size);
                let he = e.min(hit_size);
                let has_unique = hits[hi..he].iter().any(|h| h.repeats <= 10000);
                if !has_unique {
                    s = e;
                    continue;
                }
            }

            // [s, e) holds the candidate in the array of hitCoordDiff.
            concordant_hit_coord.clear();
            for k in s..e {
                concordant_hit_coord.push(Pair { a: run[k].a, b: run[k].b });
            }

            if adjust_radius > 0 {
                for pair in &concordant_hit_coord {
                    #[allow(clippy::cast_sign_loss)]
                    let a_idx = pair.a as usize;
                    let dist = (pair.a - pair.b - dominant_coord_diff).abs();
                    if read_offset_used[a_idx] == -1 || read_offset_used[a_idx] > dist {
                        read_offset_used[a_idx] = dist;
                    }
                }
                let mut l = 0usize;
                for k in 0..concordant_hit_coord.len() {
                    let pair = concordant_hit_coord[k];
                    #[allow(clippy::cast_sign_loss)]
                    let a_idx = pair.a as usize;
                    let dist = (pair.a - pair.b - dominant_coord_diff).abs();
                    if dist == read_offset_used[a_idx] {
                        concordant_hit_coord[l] = pair;
                        l += 1;
                    }
                }
                concordant_hit_coord.truncate(l);
                concordant_hit_coord.sort_unstable_by(comp_sort_pair_b_inc);
            }

            // Compute the longest increasing subsequence.
            hit_coord_lis.clear();
            let lis_size = if concordant_hit_coord.is_empty() {
                0
            } else {
                longest_increasing_subsequence(&concordant_hit_coord, &mut hit_coord_lis)
            };
            #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
            let lis_size_i32 = lis_size as i32;
            if lis_size_i32 * kmer_length < hit_len_required {
                s = e;
                continue;
            }

            // Rebuild the hits.
            let mut lis_start = 0usize;
            let mut lis_end = lis_size - 1;
            // Ignore long insert gaps (SeqSet.hpp:1473-1498). Dead for this
            // codebase (`is_long_seq_set` is always `false` -- see module
            // docs / `RefKmerFilter`'s doc comment), but ported faithfully
            // and generally.
            if is_long_seq_set {
                let mut max_gap = 2 * hit_len_required + 3 * kmer_length;
                if filter == 0 {
                    max_gap *= 4;
                }
                if max_gap < 200 {
                    max_gap = 200;
                }
                let mut max_run: i32 = -1;
                let mut k = 0usize;
                while k < lis_size {
                    let mut l = k + 1;
                    while l < lis_size {
                        if hit_coord_lis[l].a - hit_coord_lis[l - 1].a > max_gap {
                            break;
                        }
                        l += 1;
                    }
                    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                    let run_len = (l - k) as i32;
                    if run_len > max_run {
                        max_run = run_len;
                        lis_start = k;
                        lis_end = l - 1;
                    }
                    k = l;
                }
            }

            final_hits.clear();
            for k in lis_start..=lis_end {
                let mut nh = hits[i];
                nh.read_offset = hit_coord_lis[k].a;
                #[allow(clippy::cast_sign_loss)]
                {
                    nh.offset = hit_coord_lis[k].b as u32;
                }
                final_hits.push(nh);
            }
            let lis_size = lis_end - lis_start + 1;

            let hit_len = get_total_hit_length_on_read(&final_hits, kmer_length);
            if hit_len < hit_len_required {
                s = e;
                continue;
            }
            if get_total_hit_length_on_seq(&final_hits, kmer_length) < hit_len_required {
                s = e;
                continue;
            }

            let seq_idx = hits[i].idx;
            let read_start = final_hits[0].read_offset;
            let read_end = final_hits[lis_size - 1].read_offset + kmer_length - 1;
            let strand = final_hits[0].strand;
            #[allow(clippy::cast_possible_wrap)]
            let seq_start = final_hits[0].offset as i32;
            #[allow(clippy::cast_possible_wrap)]
            let seq_end = final_hits[lis_size - 1].offset as i32 + kmer_length - 1;
            let match_cnt = 2 * hit_len;

            if !is_ref(seq_idx) && hit_len * 2 < seq_end - seq_start + 1 {
                s = e;
                continue;
            }

            overlaps.push(Overlap {
                seq_idx,
                read_start,
                read_end,
                seq_start,
                seq_end,
                strand,
                match_cnt,
                similarity: 0.0,
            });

            s = e;
        }

        i = j;
    }

    overlaps
}

/// Truncates a `f64` toward zero and casts to `i32`, matching C++'s
/// `int(double)` conversion (`SeqSet.hpp:1291`: `longestHits[i] * 0.75`).
fn trunc_i32(v: f64) -> i32 {
    #[allow(clippy::cast_possible_truncation)]
    let truncated = v.trunc() as i32;
    truncated
}

/// `min_hit_required` is always non-negative in practice (`{3, 3}` defaults,
/// or a `0.75`/`/2`/`/3`/`/4`-derived value of a non-negative `longestHits`
/// entry), but is typed `i32` throughout to mirror C++'s plain `int`; this
/// helper converts it to a `usize` for comparison against `e - s` (a `usize`
/// length) the same way C++'s signed/unsigned comparison
/// (`e - s < minHitRequired`) implicitly does for non-negative values.
fn usize_from_i32_floor(v: i32) -> usize {
    usize::try_from(v.max(0)).unwrap_or(usize::MAX)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(idx: u32, offset: u32, read_offset: i32, strand: i8) -> OverlapHit {
        OverlapHit { idx, offset, read_offset, strand, repeats: 1 }
    }

    #[test]
    fn lis_hand_computed_strictly_increasing() {
        // hits already strictly increasing in both a and b: LIS is the
        // whole sequence, unchanged.
        let hits = vec![Pair { a: 0, b: 100 }, Pair { a: 9, b: 109 }, Pair { a: 18, b: 118 }];
        let mut lis = Vec::new();
        let len = longest_increasing_subsequence(&hits, &mut lis);
        assert_eq!(len, 3);
        assert_eq!(lis, hits);
    }

    #[test]
    fn lis_hand_computed_with_out_of_order_element() {
        // a=[0, 30, 9, 18] with matching b offsets: the element at a=30
        // breaks the increasing run relative to what follows (18 < 30), so
        // the LIS should skip it and pick 0, 9, 18.
        let hits = vec![
            Pair { a: 0, b: 100 },
            Pair { a: 30, b: 130 },
            Pair { a: 9, b: 109 },
            Pair { a: 18, b: 118 },
        ];
        let mut lis = Vec::new();
        let len = longest_increasing_subsequence(&hits, &mut lis);
        assert_eq!(len, 3);
        assert_eq!(lis, vec![Pair { a: 0, b: 100 }, Pair { a: 9, b: 109 }, Pair { a: 18, b: 118 }]);
    }

    #[test]
    fn lis_dedups_consecutive_equal_b() {
        // Two entries with the SAME b (109) back to back in the LIS
        // reconstruction order: only the first-encountered (post-reverse)
        // survives the SeqSet.hpp:417-429 same-b dedup pass.
        let hits = vec![
            Pair { a: 0, b: 100 },
            Pair { a: 5, b: 109 },
            Pair { a: 9, b: 109 },
            Pair { a: 18, b: 118 },
        ];
        let mut lis = Vec::new();
        let len = longest_increasing_subsequence(&hits, &mut lis);
        // Full increasing run is length 4 (0,5,9,18 strictly increasing in
        // a); after same-b dedup, the second b=109 entry (a=9) is dropped
        // since it immediately follows another b=109 entry in LIS order.
        assert_eq!(len, 3);
        assert_eq!(lis[0], Pair { a: 0, b: 100 });
        assert_eq!(lis[1], Pair { a: 5, b: 109 });
        assert_eq!(lis[2], Pair { a: 18, b: 118 });
    }

    #[test]
    fn lis_single_element() {
        let hits = vec![Pair { a: 5, b: 50 }];
        let mut lis = Vec::new();
        let len = longest_increasing_subsequence(&hits, &mut lis);
        assert_eq!(len, 1);
        assert_eq!(lis, hits);
    }

    #[test]
    fn total_hit_length_on_read_merges_overlapping_windows() {
        // kmer_length = 9: hits at read_offset 0 and 5 have overlapping
        // 9-base windows ([0,9) and [5,14)), so they merge into a single run
        // contributing (5 - 0 + 9) = 14, not 9 + 9 = 18.
        let hits = vec![hit(0, 0, 0, 1), hit(0, 5, 5, 1)];
        assert_eq!(get_total_hit_length_on_read(&hits, 9), 14);
    }

    #[test]
    fn total_hit_length_on_read_separate_runs() {
        // kmer_length = 9: hits at read_offset 0 and 20 are far apart (20 >
        // 0 + 9 - 1 = 8), so they form two separate runs: 9 + 9 = 18.
        let hits = vec![hit(0, 0, 0, 1), hit(0, 20, 20, 1)];
        assert_eq!(get_total_hit_length_on_read(&hits, 9), 18);
    }

    #[test]
    fn total_hit_length_on_seq_keys_on_seq_offset() {
        let hits = vec![hit(0, 0, 0, 1), hit(0, 5, 100, 1)];
        // Same seq-offset adjacency logic as the read version, but keyed on
        // `offset` (5) instead of `read_offset` (100): 5-0+9=14.
        assert_eq!(get_total_hit_length_on_seq(&hits, 9), 14);
    }

    #[test]
    fn get_overlaps_hand_computed_single_colinear_chain() {
        // 5 perfectly colinear hits at kmer_length=9, read_offset stepping
        // by 1 (so seq_offset also steps by 1, coord diff constant at -100):
        // this should produce exactly one overlap spanning the whole chain.
        let hits: Vec<OverlapHit> =
            (0..5i32).map(|i| hit(0, u32::try_from(100 + i).unwrap(), i, 1)).collect();
        let overlaps =
            get_overlaps_from_hits(&hits, 9 /* hitLenRequired */, 0, |_| true, 10, false, 9);
        assert_eq!(overlaps.len(), 1);
        let o = overlaps[0];
        assert_eq!(o.seq_idx, 0);
        assert_eq!(o.strand, 1);
        assert_eq!(o.read_start, 0);
        assert_eq!(o.read_end, 4 + 9 - 1); // last readOffset + kmerLength - 1
        assert_eq!(o.seq_start, 100);
        assert_eq!(o.seq_end, 104 + 9 - 1);
        // hitLen on read: single run, (4 - 0) + 9 = 13. matchCnt = 2*hitLen.
        assert_eq!(o.match_cnt, 2 * 13);
        assert!((o.similarity - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn get_overlaps_rejects_run_below_min_hit_required() {
        // Only 2 hits on this (strand, idx): below the default
        // novelMinHitRequired/refMinHitRequired of 3, so the whole run is
        // skipped before even building hitCoordDiff.
        let hits = vec![hit(0, 100, 0, 1), hit(0, 101, 1, 1)];
        let overlaps = get_overlaps_from_hits(&hits, 9, 0, |_| true, 10, false, 9);
        assert!(overlaps.is_empty());
    }

    #[test]
    fn get_overlaps_two_disjoint_chains_picks_longer_lis_only() {
        // Two colinear chains at very different coord diffs (100 vs 500),
        // each individually long enough to pass thresholds. Since both are
        // "concordant" runs found via hitCoordDiff sorting by `c`
        // (readOffset - seqOffset), they should each become sorted into
        // adjacent hitCoordDiff regions with `diff` between them exceeding
        // `radius` (10), splitting them into TWO separate concordant runs
        // -> TWO overlaps, not one merged chain (this hand-verifies the
        // hitCoordDiff-based splitting, distinct from the LIS-noncolinear
        // case covered by the end-to-end differential test).
        let mut hits: Vec<OverlapHit> = Vec::new();
        for i in 0..5i32 {
            #[allow(clippy::cast_sign_loss)]
            hits.push(hit(0, (100 + i) as u32, i, 1));
        }
        for i in 0..5i32 {
            #[allow(clippy::cast_sign_loss)]
            hits.push(hit(0, (500 + i) as u32, 20 + i, 1));
        }
        let overlaps = get_overlaps_from_hits(&hits, 9, 0, |_| true, 10, false, 9);
        assert_eq!(overlaps.len(), 2);
    }
}
