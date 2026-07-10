//! Banded, integer-scored global alignment with traceback, ported from T1K's
//! `AlignAlgo` (`AlignAlgo.hpp`), plus the standalone
//! `SeqSet::GetAlignStats` helper (`SeqSet.hpp:438-453`) that consumes an
//! `align[]` op sequence.
//!
//! # Scope: only the genotyping-path slice is ported
//!
//! `AlignAlgo` (the C++ class) declares many more static methods than are
//! ported here. Cross-referencing every `AlignAlgo::` call site
//! (`grep -rn "AlignAlgo::" *.hpp *.cpp`) shows only
//! three are ever reached from the genotyping path (`SeqSet.hpp`,
//! `VariantCaller.hpp`):
//!
//! - [`global_alignment`] -- `AlignAlgo::GlobalAlignment` (banded affine-gap
//!   `char*` vs `char*`).
//! - [`global_alignment_pos_weight`] -- `AlignAlgo::GlobalAlignment_PosWeight`
//!   (banded, NON-affine `_posWeight[]` vs `char*`; only reachable when
//!   `lent == lenp` fails the free no-indel-alignment shortcut, see below).
//! - `AlignAlgo::VisualizeAlignment` -- debug-only `printf` visualization
//!   (`AlignAlgo.hpp:1187-1231`); its output is never captured or compared by
//!   any caller (both call sites in `SeqSet.hpp`/`VariantCaller.hpp` are
//!   inside `if (VERBOSE)`-style debug gates whose only observable effect is
//!   stdout text), so it is explicitly NOT ported here.
//!
//! `SemiGlobalAlignment`, `GlobalAlignment_PosWeight_Affine`,
//! `GlobalAlignment_OneEnd`, `GlobalAlignment_classic`, `LocalAlignment`,
//! `IsMateOverlap`, `LocatePartialSufPrefExactMatch`, and
//! `LocatePartialSufSufExactMatch` have zero call sites anywhere in the
//! vendored tree and are NOT ported (avoiding dead code per this port's
//! scope discipline).
//!
//! [`get_align_stats`] mirrors `SeqSet::GetAlignStats` (`SeqSet.hpp:438-453`)
//! -- technically a `SeqSet` *member* function (despite the brief's
//! expectation that it lived on `AlignAlgo`), but it touches no `SeqSet`
//! instance state and is the natural free-function counterpart that consumes
//! an `align[]` op sequence, so it lives in this module alongside its
//! producers.
//!
//! # All-integer scoring: score AND traceback are both byte-identical
//! # contracts
//!
//! Every DP cell and every accumulated score is a plain `i32` (mirroring
//! C++'s `int`); there is no floating point anywhere in the ported DP or
//! traceback. That means not just the final score, but the *exact sequence*
//! of `EDIT_*` ops emitted by the traceback, is a reproducible, checkable
//! byte-identity target -- this port's FFI differential test
//! (`crates/unum-core/tests/golden_align_algo.rs`) asserts both.
//!
//! # THE critical trap: traceback tie-breaking
//!
//! When multiple predecessor cells tie for the max score, the C++ resolves
//! the tie via a specific *last-write-wins* comparison chain, not a
//! `>=`-earliest-wins or symmetric priority order. Getting the *score*
//! right is comparatively easy (the recurrences are `MAX(...)` chains that
//! don't care about tie order); getting the *op sequence* right requires
//! reproducing this comparison order exactly:
//!
//! - [`global_alignment_pos_weight`]'s traceback (single `m` matrix, no
//!   affine gap): per step, `a` starts at a sentinel `0`==`EDIT_MATCH`, then
//!   is unconditionally overwritten (if the corresponding score-delta
//!   equality holds) in this order: DELETE-from-left, then
//!   INSERT-from-above, then diagonal (MATCH or MISMATCH) -- so the
//!   **diagonal move wins any tie**, INSERT beats DELETE, and DELETE is only
//!   the answer if neither INSERT nor the diagonal also matched.
//! - [`global_alignment`]'s traceback (three matrices `m`/`e`/`f`, banded
//!   affine gap) is a 3-state machine (`mat` = which matrix the backtrace is
//!   currently "in"): while `mat == 0` (in `m`), `a` starts as `EDIT_INSERT`,
//!   is overwritten to `EDIT_DELETE` if `f[i,j] >= e[i,j]` (**not** `>=
//!   m[i,j]` -- that variant only appears in the unported
//!   `GlobalAlignment_PosWeight_Affine`), then overwritten again to
//!   MATCH/MISMATCH if the diagonal-from-`m` equality holds -- so **the
//!   diagonal wins**, DELETE beats the INSERT default, and switching into
//!   `mat=1`/`mat=2` (rather than emitting an op immediately) defers the
//!   actual `EDIT_INSERT`/`EDIT_DELETE` emission to the *next* loop
//!   iteration. Within `mat==1`/`mat==2`, the "stay in the gap" vs.
//!   "gap-open back to `m`" choice is also a fixed-order `if`/`else`, not a
//!   `MAX` comparison: gap-open (`mat` returns to `0`) is checked FIRST and
//!   taken if it matches, extension is the `else` fallback.
//!
//! This port reproduces both comparison chains in the exact order written
//! above (see the `let mut a = ...; if ... { a = ...; }` chains in each
//! function body, which mirror the C++ `if` statements one-for-one, not
//! collapsed into `match`/`max_by_key`). The differential test's
//! `flip_tie_break_and_confirm_align_array_catches_it` case demonstrates
//! that comparing only the score would NOT catch a tie-break regression --
//! only comparing the full `align[]` array does.
//!
//! # Lint allowances: line-comparable naming over clippy's naming/length preferences
//!
//! `lent`/`lenp`, `tagi`/`tagj`, `t`/`p`, and the `m`/`e`/`f` matrix names
//! are the vendored C++'s OWN variable names (`AlignAlgo.hpp`'s `lent`
//! /`lenp` parameters, `tagi`/`tagj` traceback cursors, `t`/`p` sequence
//! pointers, `m`/`e`/`f` DP matrices) -- kept identical here deliberately so
//! this module stays line-comparable to the source it ports, which is more
//! valuable for future maintainers cross-checking a fix against the
//! original than clippy's `similar_names`/`many_single_char_names`
//! preferences. Likewise `global_alignment`/`global_alignment_pos_weight`
//! exceed `too_many_lines` because each is a single unbroken DP+traceback
//! block in the original (splitting it into helpers would only obscure the
//! line-for-line correspondence this port's byte-identity mandate depends
//! on), and the `lent`-vs-`lenp` band-widening `if`/`else if`/`else` chains
//! are kept as `if` chains (not collapsed into `match ... .cmp(...)`) to
//! mirror `AlignAlgo.hpp`'s own `if (lent > lenp) ... else if (lent < lenp)
//! ... ` structure verbatim.
#![allow(
    clippy::similar_names,
    clippy::many_single_char_names,
    clippy::too_many_lines,
    clippy::comparison_chain
)]

/// Mirrors `EDIT_MATCH` (`AlignAlgo.hpp:7`).
pub const EDIT_MATCH: i8 = 0;
/// Mirrors `EDIT_MISMATCH` (`AlignAlgo.hpp:8`).
pub const EDIT_MISMATCH: i8 = 1;
/// Mirrors `EDIT_INSERT` (`AlignAlgo.hpp:9`).
pub const EDIT_INSERT: i8 = 2;
/// Mirrors `EDIT_DELETE` (`AlignAlgo.hpp:10`).
pub const EDIT_DELETE: i8 = 3;

const SCORE_MATCH: i32 = 2;
const SCORE_MISMATCH: i32 = -2;
const SCORE_GAPOPEN: i32 = -4;
const SCORE_GAPEXTEND: i32 = -1;
const SCORE_INDEL: i32 = -4;

/// Default band half-width, mirroring `GlobalAlignment`'s `int band = 5`
/// default parameter (`AlignAlgo.hpp:215`) and
/// `GlobalAlignment_PosWeight`'s hardcoded `int leftBand = 5; int rightBand =
/// 5;` (`AlignAlgo.hpp:106-107`).
pub const DEFAULT_BAND: i32 = 5;

/// Maximum mismatch count for which [`global_alignment`]'s equal-length
/// diagonal fast path is provably the exact DP optimum (see that function's
/// doc comment for the score arithmetic: `SCORE_MATCH - SCORE_MISMATCH = 4`
/// gain per rescued mismatch vs. `2 * (SCORE_GAPOPEN + SCORE_GAPEXTEND) =
/// -10` for the cheapest gap pair, so `4 * 2 - 10 = -2 < 0` at 2 mismatches).
const DIAGONAL_FAST_PATH_MAX_MISMATCHES: u32 = 2;

/// Per-position base counts at one reference position, ported from
/// `struct _posWeight` (`AlignAlgo.hpp:21-44`).
///
/// Only [`PosWeight::sum`] (mirroring `_posWeight::Sum()`) is needed by
/// [`is_base_equal`]; `operator+=`/`Clear()` are not used by the
/// genotyping-path slice ported here and are omitted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PosWeight {
    /// Per-base observed counts, indexed by `nucToNum['A'..'T']` (A=0, C=1,
    /// G=2, T=3) -- matches `_posWeight::count[4]`.
    pub count: [i32; 4],
}

impl PosWeight {
    /// Mirrors `_posWeight::Sum()` (`AlignAlgo.hpp:35-38`).
    #[must_use]
    pub fn sum(&self) -> i32 {
        self.count[0] + self.count[1] + self.count[2] + self.count[3]
    }
}

/// Interior-mutable variant of [`PosWeight`] used for the SHARED per-allele
/// base-coverage counters in `AlleleRef::pos_weight`, so the base-coverage
/// marking step of `assign_read` can run across `rayon` workers in parallel
/// without a data race and WITHOUT per-thread copies of the (RSS-dominating)
/// coverage structures.
///
/// The only shared mutation `assign_read` performs is `count[code] += weight`
/// (an integer add). Integer addition is associative and commutative, so the
/// FINAL sum after all workers join is independent of the order the adds land
/// -- `fetch_add(_, Relaxed)` therefore yields byte-identical results to the
/// sequential `+=` loop while `Relaxed` avoids any unnecessary fences (we only
/// require the exact total once all threads have joined; no add ever depends
/// on the value another add produced, and no code reads `count` mid-pass, so
/// no `Acquire`/`Release` ordering is needed). Each `AtomicI32` is the same
/// size as an `i32`, so switching to this type does not grow peak RSS.
///
/// Reads (`sum`, `count`) are only ever performed AFTER the parallel
/// assignment pass has fully joined (`get_seq_missing_base_coverage`,
/// `exon_base_coverage`), where `Relaxed` loads observe every prior
/// `fetch_add`.
#[derive(Debug, Default)]
pub struct AtomicPosWeight {
    /// Per-base observed counts, indexed by `nucToNum['A'..'T']` (A=0, C=1,
    /// G=2, T=3) -- the atomic counterpart of [`PosWeight::count`].
    pub count: [std::sync::atomic::AtomicI32; 4],
}

impl AtomicPosWeight {
    /// Atomically adds `weight` to the `code`-th base count (mirrors
    /// `count[code] += weight`). `Relaxed` is sufficient -- see the type-level
    /// doc comment.
    pub fn add(&self, code: usize, weight: i32) {
        self.count[code].fetch_add(weight, std::sync::atomic::Ordering::Relaxed);
    }

    /// Loads the `code`-th base count. `Relaxed` is sufficient; only ever
    /// called after the parallel pass has joined.
    #[must_use]
    pub fn get(&self, code: usize) -> i32 {
        self.count[code].load(std::sync::atomic::Ordering::Relaxed)
    }
}

impl Clone for AtomicPosWeight {
    /// A fresh, independent copy of the current counts (loads each atomic with
    /// `Relaxed`). Only used off the hot path -- `AlleleRef::pos_weight` is
    /// built once and shared by reference during assignment, never cloned per
    /// read.
    fn clone(&self) -> Self {
        use std::sync::atomic::{AtomicI32, Ordering};
        Self {
            count: [
                AtomicI32::new(self.count[0].load(Ordering::Relaxed)),
                AtomicI32::new(self.count[1].load(Ordering::Relaxed)),
                AtomicI32::new(self.count[2].load(Ordering::Relaxed)),
                AtomicI32::new(self.count[3].load(Ordering::Relaxed)),
            ],
        }
    }
}

/// Maps `A`/`C`/`G`/`T` to their `nucToNum` index (0-3), matching the table
/// defined in `Genotyper.cpp:37-40` (the binary that actually links
/// `SeqSet`/`AlignAlgo`): every letter other than `A`/`C`/`G`/`T` maps to
/// `-1` in that table, INCLUDING `N` (unlike some other `nucToNum` copies
/// elsewhere in the vendored tree, e.g. `BamExtractor.cpp`, which map `N` to
/// `0`). [`is_base_equal`] never actually indexes with an `N` code (its `c ==
/// 'N'` check short-circuits first, exactly as the C++ does), so this
/// distinction is inert in practice.
///
/// Returns `None` for any byte the vendored table maps to `-1` (every
/// non-`A`/`C`/`G`/`T` letter, per `Genotyper.cpp:37-40` -- INCLUDING `N`,
/// though [`is_base_equal`] never actually reaches this function with `c ==
/// 'N'`) OR any byte outside `'A'..='Z'` (the C++ `nucToNum[c - 'A']`
/// indexing is undefined behavior there, which this port does not attempt
/// to replicate; only `'A'..='Z'` inputs are in scope, matching this port's
/// general uppercase-input assumption). Returns `Some(0..=3)` directly
/// (rather than the vendored table's raw `i8`, which is never negative for
/// any value this function actually returns) since the only consumer,
/// [`is_base_equal`], immediately uses the result as a `count` array index.
fn nuc_to_num(c: u8) -> Option<usize> {
    match c {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Mirrors the private `AlignAlgo::IsBaseEqual` (`AlignAlgo.hpp:49-55`):
/// decides whether a single reference position's base-count distribution
/// `w` is compatible with observed base `c`.
///
/// Three ways to be considered "equal": `w` has zero total support (`sum ==
/// 0`, an unobserved / all-gap column), `c == 'N'` (wildcard), or `c`'s own
/// count is less than a third of the total (i.e. `c` is not the dominant
/// call at this position, so treat it charitably as consistent with the
/// consensus rather than a hard mismatch). Note the `<` (not `<=`): a
/// position with `sum == 3 * count[c]` exactly (e.g. `count[c] == sum`, `c`
/// unanimous) does NOT count as a free match by this branch -- but a
/// unanimous position IS still a match by construction, since `sum <
/// 3*count[c]` reduces to `count[c] > sum/3`, which is satisfied whenever
/// `c` is not a strict minority (<=1/3) of the votes.
///
/// # Panics
///
/// Panics if `c` is not `'N'` and has no `nucToNum` mapping (see
/// [`nuc_to_num`]) -- mirrors the C++ side's undefined-behavior-on-garbage-
/// input, upgraded to an explicit panic since Rust has no silent UB
/// equivalent; only ever called with `'A'`/`'C'`/`'G'`/`'T'`/`'N'` in
/// practice.
#[must_use]
pub fn is_base_equal(w: &PosWeight, c: u8) -> bool {
    let sum = w.sum();
    if sum == 0 || c == b'N' {
        return true;
    }
    let idx = nuc_to_num(c).unwrap_or_else(|| panic!("is_base_equal: unmapped base {c:?}"));
    sum < 3 * w.count[idx]
}

/// A single reference base is compatible with an observed base under
/// [`global_alignment`]'s inline match test: exact equality, OR either side
/// is `'N'` (wildcard). Mirrors the repeated C++ expression `t[j-1] ==
/// p[i-1] || t[j-1] == 'N' || p[i-1] == 'N'` (`AlignAlgo.hpp:304-305,339,342`
/// etc.) and the `lent==1 && lenp==1` shortcut's `t[0] == p[0] || t[0] ==
/// 'N' || p[0] == 'N'` (`AlignAlgo.hpp:224`).
fn chars_match(t: u8, p: u8) -> bool {
    t == p || t == b'N' || p == b'N'
}

/// Result of a ported `AlignAlgo` alignment call: the integer score plus the
/// traceback op sequence (already reversed into left-to-right order, exactly
/// as the C++ `align[]` output array reads after its final in-place
/// reverse). No `-1` sentinel terminator is included -- unlike the C++
/// caller-allocated buffer, a `Vec<i8>`'s own length IS the terminator.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignResult {
    pub score: i32,
    pub align: Vec<i8>,
}

/// Per-thread memoization cache for the DP alignment path of
/// [`global_alignment`] (folds fork commit `a35ed72`).
///
/// # Why memoization is a large, byte-identical win here
///
/// During genotyping a read is aligned against ~212 highly similar candidate
/// HLA alleles, and there are tens of thousands of near-identical alleles in
/// the reference. The fork measured ~114x redundancy: 273.8M `GlobalAlignment`
/// calls resolved to only ~2.4M **distinct** `(t, p)` pairs. Because the DP is
/// a pure, deterministic function of `(t, p)` (same inputs => same score AND
/// same traceback `align[]`), returning a cached result for a repeated
/// `(t, p)` is *byte-identical* to recomputing it -- the exact property this
/// port's FFI differentials assert.
///
/// # Reused key buffer (the lesson from the matrix-reuse regression)
///
/// [`dp_key`] is a single `Vec<u8>` reused across every call
/// (`clear()` + rebuild), NOT allocated per call -- so a cache HIT costs one
/// hash of borrowed bytes and zero heap allocation. Only a cache MISS
/// allocates (the owned key to insert + the stored `align` copy), and misses
/// are ~1/114 of calls, so that allocation no longer dominates.
///
/// Construct one per worker thread (inside the per-rayon-worker [`Scratch`]),
/// which is the lock-free thread-local equivalent of the fork's
/// `thread_local unordered_map`.
///
/// [`Scratch`]: crate::ref_kmer_filter::Scratch
/// [`dp_key`]: DpCache::dp_key
#[derive(Debug, Default, Clone)]
pub struct DpCache {
    /// `(lent, lenp, t, p)`-keyed cache of `(score, align)` DP results. The key
    /// frames both lengths and both full sequences so two different `(t, p)`
    /// inputs can never alias to the same key (see [`DpCache::build_key`]).
    ///
    /// Keyed with [`rustc_hash::FxHashMap`] rather than std's SipHash: the key
    /// spans both full sequences (tens of bytes) and is hashed on *every* call
    /// (the ~99% cache-hit path), so per-byte hash speed dominates here. This is
    /// a pure memoization cache -- the hasher cannot affect which stored value a
    /// key maps to -- so the swap is byte-identical by construction.
    cache: rustc_hash::FxHashMap<Vec<u8>, (i32, Vec<i8>)>,
    /// Reused scratch buffer for the lookup key; cleared and rebuilt each call
    /// so a cache hit performs no heap allocation.
    dp_key: Vec<u8>,
    /// Reused DP score matrices (`m`/`e`/`f`) for the cache-miss path, so a miss
    /// re-zeroes and refills these buffers rather than allocating three fresh
    /// `(lenp+1)*(lent+1)` matrices per DP. Held here (per-worker, alongside the
    /// cache they back) since the miss path always has a `&mut DpCache`.
    dp_m: Matrix,
    dp_e: Matrix,
    dp_f: Matrix,
}

impl DpCache {
    /// Create an empty cache.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Rebuild [`Self::dp_key`] in place for the pair `(t, p)`.
    ///
    /// Layout: `lent` (8 bytes LE) `lenp` (8 bytes LE) `t[..]` `0xFF` `p[..]`.
    /// Including both explicit lengths AND a `0xFF` separator between the two
    /// sequences guarantees the byte string uniquely determines `(t, p)`:
    /// `t` and `p` are both plain uppercase nucleotide/`N` bytes (never `0xFF`),
    /// and the leading lengths pin the split even if `0xFF` ever appeared, so
    /// no two distinct `(t, p)` pairs can collide onto one key.
    fn build_key(&mut self, t: &[u8], p: &[u8]) {
        self.dp_key.clear();
        self.dp_key.extend_from_slice(&(t.len() as u64).to_le_bytes());
        self.dp_key.extend_from_slice(&(p.len() as u64).to_le_bytes());
        self.dp_key.extend_from_slice(t);
        self.dp_key.push(0xFF);
        self.dp_key.extend_from_slice(p);
    }
}

/// 2D DP matrix stored as a flat `Vec<i32>`, row-major with `bmax = lent + 1`
/// columns, mirroring the C++ side's manual `m[i * bmax + j]` flat indexing
/// exactly (kept flat, not `Vec<Vec<i32>>`, so index arithmetic is visibly
/// identical to the ported C++ at every call site).
#[derive(Debug, Default, Clone)]
struct Matrix {
    data: Vec<i32>,
    bmax: usize,
}

impl Matrix {
    fn new(lenp: usize, lent: usize) -> Self {
        Self { data: vec![0; (lenp + 1) * (lent + 1)], bmax: lent + 1 }
    }

    /// Reuse this matrix's backing allocation for a fresh `(lenp, lent)` DP,
    /// re-zeroing exactly `(lenp + 1) * (lent + 1)` cells. Byte-identical to a
    /// freshly [`Matrix::new`]-allocated matrix (same all-zero contents, same
    /// `bmax`), but reuses the heap buffer across calls instead of allocating
    /// and freeing a new one every DP -- the DP runs on the ~0.8% cache-miss
    /// path ~4x10^6 times per sample, each allocating three of these.
    fn reset(&mut self, lenp: usize, lent: usize) {
        self.data.clear();
        self.data.resize((lenp + 1) * (lent + 1), 0);
        self.bmax = lent + 1;
    }

    #[inline]
    fn get(&self, i: usize, j: usize) -> i32 {
        self.data[i * self.bmax + j]
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, v: i32) {
        self.data[i * self.bmax + j] = v;
    }

    /// Signed-cursor accessor used ONLY by
    /// [`global_alignment_pos_weight`]'s traceback, which (mirroring
    /// `AlignAlgo.hpp`'s `int tagi, tagj`) allows `i`/`j` to go transiently
    /// negative. C++ computes the flat offset as `tagi * bmax + tagj` using
    /// plain (possibly negative) `int` arithmetic and indexes the same
    /// backing array with it -- e.g. `tagj == -1` with `tagi >= 1` reads flat
    /// offset `tagi*bmax - 1`, which lands one element before row `tagi`,
    /// i.e. the LAST element of row `tagi - 1` (since each row is `bmax`
    /// elements wide). This is undefined behavior in C++ but is, in
    /// practice, an in-bounds read of that aliased cell for every reachable
    /// caller input (the loop terminates before `tagi`/`tagj` can run far
    /// enough negative to leave the backing allocation) -- reproduced here
    /// via the same flat-offset arithmetic so the traceback stays
    /// byte-identical instead of panicking on `usize` underflow.
    #[inline]
    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn get_signed(&self, i: isize, j: isize) -> i32 {
        let bmax = self.bmax as isize;
        let flat = i * bmax + j;
        self.data[flat as usize]
    }
}

/// Ported from `AlignAlgo::GlobalAlignment_PosWeight`
/// (`AlignAlgo.hpp:57-213`): banded, NON-affine (single linear-gap `m`
/// matrix, `SCORE_INDEL = -4` per gap base) global alignment of a
/// per-position base-count reference (`t_weights`) against a plain sequence
/// `p`.
///
/// # Early-exit shortcuts (ported verbatim, not simplifications)
///
/// - `lent == 0 || lenp == 0`: returns score `0` with an empty `align`.
/// - `lent == 1 && lenp == 1`: single-base compare via [`is_base_equal`],
///   bypassing the DP entirely.
/// - `lent == lenp`: first tries a **free no-indel alignment** (position-by-
///   position match/mismatch, no DP). If that no-indel score already beats
///   `lent * SCORE_MATCH + 2 * SCORE_INDEL` (i.e. is at least as good as the
///   best possible score achievable WITH at least one insertion and one
///   deletion), it is returned immediately without ever running the banded
///   DP -- the reasoning being that any real indel pair costs more than
///   `2 * SCORE_INDEL` worth of penalty, so it can never beat a
///   surprisingly-good no-indel alignment. This shortcut can change which
///   `align[]` is returned relative to running the DP unconditionally: the
///   DP might find an equally-scoring alignment with a different op
///   sequence, but the shortcut always wins ties in favor of the pure
///   diagonal walk.
///
/// # Panics
///
/// Panics if `t_weights.len() != lent` or `p.len() != lenp` (a caller
/// contract violation, not a data-dependent condition); mirrors the C++
/// side's undefined behavior on a length mismatch (both are separate
/// caller-supplied lengths in the original signature) upgraded to an
/// explicit panic.
#[must_use]
pub fn global_alignment_pos_weight(t_weights: &[PosWeight], p: &[u8]) -> AlignResult {
    let lent = t_weights.len();
    let lenp = p.len();

    if lent == 0 || lenp == 0 {
        return AlignResult { score: 0, align: Vec::new() };
    }
    if lent == 1 && lenp == 1 {
        return if is_base_equal(&t_weights[0], p[0]) {
            AlignResult { score: SCORE_MATCH, align: vec![EDIT_MATCH] }
        } else {
            AlignResult { score: SCORE_MISMATCH, align: vec![EDIT_MISMATCH] }
        };
    }

    if lent == lenp {
        let mut score = 0;
        let mut align = Vec::with_capacity(lent);
        for i in 0..lent {
            if is_base_equal(&t_weights[i], p[i]) {
                align.push(EDIT_MATCH);
                score += SCORE_MATCH;
            } else {
                align.push(EDIT_MISMATCH);
                score += SCORE_MISMATCH;
            }
        }
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let lent_i32 = lent as i32;
        if score >= lent_i32 * SCORE_MATCH + 2 * SCORE_INDEL {
            return AlignResult { score, align };
        }
    }

    let left_band: i32 = 5;
    let right_band: i32 = 5;
    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
    let (lent_i32, lenp_i32) = (lent as i32, lenp as i32);
    let (left_band, right_band) = if lent_i32 > lenp_i32 {
        (left_band, right_band + (lent_i32 - lenp_i32))
    } else if lent_i32 < lenp_i32 {
        (left_band + (lenp_i32 - lent_i32), right_band)
    } else {
        (left_band, right_band)
    };

    let neg_inf = (lent_i32 + 1) * (lenp_i32 + 1) * SCORE_INDEL;

    let mut m = Matrix::new(lenp, lent);
    m.set(0, 0, 0);
    for i in 1..=lenp {
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let i_i32 = i as i32;
        m.set(i, 0, SCORE_INDEL + i_i32 * SCORE_INDEL);
    }
    for j in 1..=lent {
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let j_i32 = j as i32;
        m.set(0, j, SCORE_INDEL + j_i32 * SCORE_INDEL);
    }

    for i in 1..=lenp {
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let i_i32 = i as i32;
        let start = if i_i32 - left_band < 1 { 1 } else { i_i32 - left_band };
        let end = if i_i32 + right_band > lent_i32 { lent_i32 } else { i_i32 + right_band };
        #[allow(clippy::cast_sign_loss)]
        let (start_u, end_u) = (start as usize, end as usize);

        if start > 1 {
            m.set(i, start_u - 1, neg_inf);
        }
        if end < lent_i32 {
            m.set(i, end_u + 1, neg_inf);
        }

        for j in start_u..=end_u {
            let diag = m.get(i - 1, j - 1)
                + if is_base_equal(&t_weights[j - 1], p[i - 1]) {
                    SCORE_MATCH
                } else {
                    SCORE_MISMATCH
                };
            let mut score = diag;
            score = score.max(m.get(i, j - 1) + SCORE_INDEL);
            score = score.max(m.get(i - 1, j) + SCORE_INDEL);
            m.set(i, j, score);
        }
    }

    let ret = m.get(lenp, lent);

    // Trace back. See module docs "THE critical trap" for the exact
    // last-write-wins tie-break order this reproduces: DELETE, then INSERT
    // (overwrites), then diagonal MATCH/MISMATCH (overwrites again and wins
    // any tie).
    //
    // Cursors are SIGNED (`isize`), mirroring C++'s `int tagi, tagj`
    // (`AlignAlgo.hpp:159`) exactly. This matters: when `tagj == 0 && tagi >
    // 0` and none of the DELETE/INSERT/diagonal branches fire (all three are
    // gated on `tagj > 0`/`tagi > 0`), the sentinel default `a == 0`
    // (`EDIT_MATCH`) survives and the final `else` arm decrements BOTH
    // cursors, sending `tagj` to `-1`. C++'s `int` tolerates this (the loop
    // still terminates once both cursors are `<= 0`); a `usize` cursor here
    // would panic with "attempt to subtract with overflow" on that same
    // input. Using `isize` + [`Matrix::get_signed`] reproduces the C++
    // pointer-arithmetic read (including its aliased-cell behavior) instead.
    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
    let mut tagi: isize = lenp as isize;
    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
    let mut tagj: isize = lent as isize;
    let mut align: Vec<i8> = Vec::new();

    while tagi > 0 || tagj > 0 {
        let max = m.get_signed(tagi, tagj);
        let mut a: i8 = 0; // sentinel default, mirrors C++ `int a = 0` (== EDIT_MATCH's value)

        if tagj > 0 && m.get_signed(tagi, tagj - 1) + SCORE_INDEL == max {
            a = EDIT_DELETE;
        }
        if tagi > 0 && m.get_signed(tagi - 1, tagj) + SCORE_INDEL == max {
            a = EDIT_INSERT;
        }
        if tagj > 0 && tagi > 0 {
            #[allow(clippy::cast_sign_loss)]
            let diff = if is_base_equal(&t_weights[(tagj - 1) as usize], p[(tagi - 1) as usize]) {
                SCORE_MATCH
            } else {
                SCORE_MISMATCH
            };
            if m.get_signed(tagi - 1, tagj - 1) + diff == max {
                a = if diff == SCORE_MATCH { EDIT_MATCH } else { EDIT_MISMATCH };
            }
        }

        align.push(a);
        if a == EDIT_DELETE {
            tagj -= 1;
        } else if a == EDIT_INSERT {
            tagi -= 1;
        } else {
            tagi -= 1;
            tagj -= 1;
        }
    }

    align.reverse();
    AlignResult { score: ret, align }
}

/// Whether [`global_alignment`]'s equal-length diagonal fast path applies to
/// `t`/`p`: both the same length AND differing by at most
/// [`DIAGONAL_FAST_PATH_MAX_MISMATCHES`] mismatches (the exact condition under
/// which the no-indel diagonal is provably the unique DP optimum -- see
/// [`global_alignment`]'s doc comment). Extracted so the `<= 2` cutoff is
/// directly unit-testable at its boundary: the fast path and the full DP are
/// byte-identical for every equal-length input at this mismatch count, so the
/// contract cannot otherwise be observed through `global_alignment`'s output
/// alone (a loosened cutoff would still produce identical results because the
/// bound is conservative, not tight).
fn diagonal_fast_path_applies(t: &[u8], p: &[u8]) -> bool {
    if t.len() != p.len() {
        return false;
    }
    let mut mismatches = 0u32;
    for k in 0..t.len() {
        if !chars_match(t[k], p[k]) {
            mismatches += 1;
            if mismatches > DIAGONAL_FAST_PATH_MAX_MISMATCHES {
                return false;
            }
        }
    }
    true
}

/// Ported from `AlignAlgo::GlobalAlignment` (`AlignAlgo.hpp:215-421`):
/// banded, affine-gap (`SCORE_GAPOPEN` + `SCORE_GAPEXTEND` per run) global
/// alignment of two plain sequences. `band` mirrors the C++ `int band = 5`
/// default parameter; pass [`DEFAULT_BAND`] to match every call site in the
/// vendored genotyping path (`SeqSet.hpp`/`VariantCaller.hpp` never pass a
/// non-default `band`).
///
/// Three DP matrices (`m` = match/mismatch "main" state, `e` = insertion
/// -to-the-text gap state, `f` = deletion-from-the-text gap state) are
/// tracked in lockstep, banded to `[i - left_band, i + right_band]` per row
/// (asymmetric when `lent != lenp`, widened on whichever side accommodates
/// the length difference -- see the C++ comments this mirrors verbatim).
///
/// # Diagonal fast path (perf, not a behavior change)
///
/// When `t.len() == p.len()` and the two sequences differ by at most 2
/// mismatches, the no-indel diagonal alignment is returned directly,
/// bypassing the DP below. This is provably byte-identical to running the
/// full DP -- not merely a good approximation -- because for equal-length
/// sequences any alignment using a gap needs at least one insertion AND one
/// deletion, so it is always strictly worse than the diagonal when there are
/// at most 2 mismatches to "fix" (see the inline comment at the shortcut for
/// the score arithmetic). Folds fork commit `b11027d`.
///
/// # Panics
///
/// Panics if `t.len()` or `p.len()` do not fit an `i32`, or if `band` is
/// negative (not a data-dependent condition; every real caller passes
/// [`DEFAULT_BAND`]).
#[must_use]
pub fn global_alignment(t: &[u8], p: &[u8], band: i32) -> AlignResult {
    global_alignment_impl(t, p, band, None)
}

/// [`global_alignment`] variant that memoizes the DP path through a per-thread
/// [`DpCache`]. The cheap early-exits (empty/single-base) and the diagonal
/// fast path run first (never cached -- they are already O(len) with no
/// allocation churn); only the expensive banded affine-gap DP is looked up in,
/// and stored to, `cache`. Byte-identical to [`global_alignment`] because the
/// DP is a pure function of `(t, p)` -- a cache hit returns the exact same
/// `(score, align)` a recompute would. Folds fork commit `a35ed72`.
#[must_use]
pub fn global_alignment_cached(t: &[u8], p: &[u8], band: i32, cache: &mut DpCache) -> AlignResult {
    global_alignment_impl(t, p, band, Some(cache))
}

/// Match/mismatch/indel tallies from an [`AlignResult::align`], as
/// [`get_align_stats`] with `update = false` would produce them.
pub type AlignStats = (i32, i32, i32);

/// Stats-only counterpart of [`global_alignment_cached`]: returns the
/// `(match, mismatch, indel)` tallies for the alignment of `(t, p)` WITHOUT
/// materializing an owned `align` `Vec` for the caller.
///
/// # Why this exists (the allocation it removes)
///
/// The read->allele overlap phase (~90% of genotype wall time) calls the
/// aligner ~5x10^8 times per sample, and every hot-path caller
/// ([`RefKmerFilter::get_overlaps_from_read`]) does exactly one thing with the
/// result: feed `align` to [`get_align_stats`] for the `(match, mismatch,
/// indel)` counts, then drop it. Going through [`global_alignment_cached`] there
/// allocates a fresh `Vec<i8>` on *every* call -- a diagonal-fast-path build, a
/// single/empty-case `vec!`, or (on the ~99% DP cache hit) a `clone()` of the
/// stored alignment -- purely to be counted and discarded. This entry point
/// computes the same tallies in place: the fast paths tally directly, and the
/// cache hit tallies the borrowed stored slice with no clone. Only a genuine
/// cache miss (~0.8% of calls) still runs the DP, and it moves (not clones) the
/// alignment into the cache.
///
/// # Byte-identical invariant
///
/// This MUST resolve `(t, p)` through the exact same branch ladder as
/// [`global_alignment_impl`] (empty/single -> diagonal fast path -> cache
/// hit/miss) so the tallies equal
/// `get_align_stats(&global_alignment_cached(t, p, band, cache).align, false, ..)`
/// for every input. The DP core, the cache, and the tally rule are all shared
/// (`global_alignment_dp`, `cache`, [`get_align_stats`]); only the branch
/// predicates (`diagonal_fast_path_applies`, the length checks, `chars_match`)
/// are mirrored here, and they mirror [`global_alignment_impl`]'s verbatim.
///
/// [`RefKmerFilter::get_overlaps_from_read`]: crate::ref_kmer_filter::RefKmerFilter::get_overlaps_from_read
#[must_use]
pub fn global_alignment_cached_stats(t: &[u8], p: &[u8], band: i32, cache: &mut DpCache) -> AlignStats {
    let lent = t.len();
    let lenp = p.len();

    // Empty: global_alignment_impl returns an empty `align` -> all-zero tallies.
    if lent == 0 || lenp == 0 {
        return (0, 0, 0);
    }
    // Single base: one EDIT_MATCH or one EDIT_MISMATCH.
    if lent == 1 && lenp == 1 {
        return if chars_match(t[0], p[0]) { (1, 0, 0) } else { (0, 1, 0) };
    }
    // Diagonal fast path: equal-length, <=2 mismatches, pure diagonal (no
    // indels). Mirrors the impl's `for k in 0..lent { push MATCH/MISMATCH }`.
    if diagonal_fast_path_applies(t, p) {
        let mut match_cnt = 0;
        let mut mismatch_cnt = 0;
        for k in 0..lent {
            if chars_match(t[k], p[k]) {
                match_cnt += 1;
            } else {
                mismatch_cnt += 1;
            }
        }
        return (match_cnt, mismatch_cnt, 0);
    }

    // Cache path: tally the stored alignment in place. A hit borrows it (no
    // clone); a miss runs the DP once, tallies the owned result, then moves it
    // into the cache under the (already-built) key.
    cache.build_key(t, p);
    if let Some((_, align)) = cache.cache.get(&cache.dp_key) {
        return stats_of(align);
    }
    let result = global_alignment_dp_with(t, p, band, &mut cache.dp_m, &mut cache.dp_e, &mut cache.dp_f);
    let stats = stats_of(&result.align);
    cache.cache.insert(cache.dp_key.clone(), (result.score, result.align));
    stats
}

/// `(match, mismatch, indel)` tally of an `align[]` op sequence -- the
/// allocation-free core shared by [`get_align_stats`] (`update = false`) and
/// [`global_alignment_cached_stats`].
#[inline]
fn stats_of(align: &[i8]) -> AlignStats {
    let mut match_cnt = 0;
    let mut mismatch_cnt = 0;
    let mut indel_cnt = 0;
    for &op in align {
        if op == EDIT_MATCH {
            match_cnt += 1;
        } else if op == EDIT_MISMATCH {
            mismatch_cnt += 1;
        } else {
            indel_cnt += 1;
        }
    }
    (match_cnt, mismatch_cnt, indel_cnt)
}

/// Shared implementation for [`global_alignment`] (no cache) and
/// [`global_alignment_cached`] (per-thread memoized DP path). See those
/// functions' docs.
#[must_use]
#[allow(clippy::too_many_lines)]
fn global_alignment_impl(
    t: &[u8],
    p: &[u8],
    band: i32,
    mut cache: Option<&mut DpCache>,
) -> AlignResult {
    let lent = t.len();
    let lenp = p.len();

    if lent == 0 || lenp == 0 {
        return AlignResult { score: 0, align: Vec::new() };
    }
    if lent == 1 && lenp == 1 {
        return if chars_match(t[0], p[0]) {
            AlignResult { score: SCORE_MATCH, align: vec![EDIT_MATCH] }
        } else {
            AlignResult { score: SCORE_MISMATCH, align: vec![EDIT_MISMATCH] }
        };
    }

    // Diagonal fast path (folds fork commit `b11027d`): for equal-length
    // sequences with at most 2 mismatches, the pure diagonal (no-indel)
    // alignment is PROVABLY the unique optimum, so it can be returned
    // directly without running the banded affine-gap DP below.
    //
    // Why it's byte-identical (not just "as good"): any gapped global
    // alignment of two equal-length sequences needs at least one insertion
    // AND at least one deletion (to stay equal-length end-to-end), costing
    // at least `2 * (SCORE_GAPOPEN + SCORE_GAPEXTEND)` = `2 * (-4 + -1)` =
    // `-10`. Each mismatch a gap pair could possibly rescue (turn into a
    // match) is worth at most `SCORE_MATCH - SCORE_MISMATCH` = `2 - (-2)` =
    // `+4`. With at most 2 mismatches, the best possible gain from adding a
    // gap pair is `4 * 2 - 10 = -2 < 0`, so the diagonal STRICTLY beats
    // every gapped alternative -- there is no tie for the DP to break
    // differently. The DP would therefore compute the exact same score and
    // the exact same all-`EDIT_MATCH`/`EDIT_MISMATCH` traceback, so skipping
    // it changes nothing observable.
    if diagonal_fast_path_applies(t, p) {
        let mut align = Vec::with_capacity(lent);
        let mut score: i32 = 0;
        for k in 0..lent {
            if chars_match(t[k], p[k]) {
                align.push(EDIT_MATCH);
                score += SCORE_MATCH;
            } else {
                align.push(EDIT_MISMATCH);
                score += SCORE_MISMATCH;
            }
        }
        return AlignResult { score, align };
    }

    // DP memoization (folds fork commit `a35ed72`): AFTER the diagonal fast
    // path, BEFORE allocating the DP matrices. A read is aligned against ~212
    // highly similar candidate alleles, so the same `(t, p)` recurs ~114x on
    // average; caching each distinct DP once is a pure, byte-identical win
    // (the DP is a deterministic function of `(t, p)`). The key buffer is
    // reused (no per-call allocation on a hit); only a miss allocates.
    let Some(cache) = cache.as_mut() else {
        return global_alignment_dp(t, p, band);
    };
    cache.build_key(t, p);
    if let Some((score, align)) = cache.cache.get(&cache.dp_key) {
        return AlignResult { score: *score, align: align.clone() };
    }

    // Cache miss: run the DP (reusing the per-worker score matrices), then
    // store it under the (already-built) key.
    let result = global_alignment_dp_with(t, p, band, &mut cache.dp_m, &mut cache.dp_e, &mut cache.dp_f);
    cache.cache.insert(cache.dp_key.clone(), (result.score, result.align.clone()));
    result
}

/// The banded affine-gap DP core of [`global_alignment`], factored out so both
/// the cached and uncached entry points share one identical implementation.
/// Callers must have already handled the empty/single-base/diagonal fast paths;
/// this is the expensive path the memoization cache wraps.
#[must_use]
fn global_alignment_dp(t: &[u8], p: &[u8], band: i32) -> AlignResult {
    // Uncached callers ([`global_alignment`], the no-cache branch): allocate
    // three throwaway matrices for this one DP.
    let (mut m, mut e, mut f) = (Matrix::default(), Matrix::default(), Matrix::default());
    global_alignment_dp_with(t, p, band, &mut m, &mut e, &mut f)
}

/// [`global_alignment_dp`] over caller-provided score matrices, so the
/// cache-miss hot path can reuse three per-worker [`Matrix`] buffers (held in
/// [`DpCache`]) instead of allocating and freeing three per call. `m`/`e`/`f`
/// are [`Matrix::reset`]-re-zeroed here before use, so their prior contents and
/// capacities are irrelevant on entry -- the result is byte-identical to
/// running on freshly allocated matrices.
#[must_use]
#[allow(clippy::too_many_lines)]
fn global_alignment_dp_with(
    t: &[u8],
    p: &[u8],
    band: i32,
    m: &mut Matrix,
    e: &mut Matrix,
    f: &mut Matrix,
) -> AlignResult {
    let lent = t.len();
    let lenp = p.len();

    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
    let (lent_i32, lenp_i32) = (lent as i32, lenp as i32);
    let (left_band, right_band) = if lent_i32 > lenp_i32 {
        (band, band + (lent_i32 - lenp_i32))
    } else if lent_i32 < lenp_i32 {
        (band + (lenp_i32 - lent_i32), band)
    } else {
        (band, band)
    };

    let neg_inf = (lent_i32 + 1) * (lenp_i32 + 1) * SCORE_GAPOPEN;

    m.reset(lenp, lent);
    e.reset(lenp, lent);
    f.reset(lenp, lent);

    m.set(0, 0, 0);
    e.set(0, 0, 0);
    f.set(0, 0, 0);

    // Boundary rows/columns mirror AlignAlgo.hpp:256-266. The col-0 (`i`) and
    // row-0 (`j`) gap-open/extend inits below are ordinary Gotoh boundaries,
    // with ONE subtlety in the row-0 `e` cell -- see `row0_e_sentinel` below.
    for i in 1..=lenp {
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let i_i32 = i as i32;
        e.set(i, 0, SCORE_GAPOPEN + i_i32 * SCORE_GAPEXTEND);
        f.set(i, 0, SCORE_GAPOPEN + i_i32 * SCORE_GAPOPEN);
        m.set(i, 0, SCORE_GAPOPEN + i_i32 * SCORE_GAPOPEN);
    }

    // Row-0 `e` is an UNREACHABLE gap-state boundary: with zero characters of
    // `p` consumed you cannot already be extending an `e`-direction gap, so this
    // cell must act as a -inf sentinel the traceback never enters from row 0.
    //
    // T1K reaches that sentinel BY ACCIDENT: `AlignAlgo.hpp:268` initializes
    // `e[0][j]` with `SCORE_GAPOPEN + i * SCORE_GAPOPEN`, reading the loop
    // variable `i` LEFT OVER from the preceding `for (i = 1; i <= lenp; ++i)`
    // loop (`i == lenp + 1`), NOT the column index `j`. The result is a single
    // large-negative constant for every column -- exactly the sentinel we want.
    // We keep that EXACT value: it is byte-identical to T1K (the goldens enforce
    // it) AND functionally correct.
    //
    // DO NOT "correct" this to a `j`-indexed value (`SCORE_GAPOPEN + j *
    // SCORE_GAPOPEN`, the "intended" formula): that makes the impossible state
    // cheaply reachable and REGRESSES real alignments. Verified 2026-07-06 --
    // the `j`-indexed form drove `golden_align_algo`'s `length_skew lent>>lenp
    // band=1` case from the correct score -22 to -63 (a lone match + all-deletes
    // instead of a grouped local alignment). This cell is live, not dead: the
    // `mat == 0` traceback reads it at `tagi == 0` (`AlignAlgo.hpp:332`).
    let row0_e_sentinel = SCORE_GAPOPEN + (lenp_i32 + 1) * SCORE_GAPOPEN;
    for j in 1..=lent {
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let j_i32 = j as i32;
        f.set(0, j, SCORE_GAPOPEN + j_i32 * SCORE_GAPEXTEND);
        e.set(0, j, row0_e_sentinel);
        m.set(0, j, SCORE_GAPOPEN + j_i32 * SCORE_GAPOPEN);
    }

    for i in 1..=lenp {
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let i_i32 = i as i32;
        let start = if i_i32 - left_band < 1 { 1 } else { i_i32 - left_band };
        let end = if i_i32 + right_band > lent_i32 { lent_i32 } else { i_i32 + right_band };
        #[allow(clippy::cast_sign_loss)]
        let (start_u, end_u) = (start as usize, end as usize);

        if start > 1 {
            let j = start_u - 1;
            e.set(i, j, neg_inf);
            f.set(i, j, neg_inf);
            m.set(i, j, neg_inf);
        }
        if end < lent_i32 {
            let j = end_u + 1;
            e.set(i, j, neg_inf);
            f.set(i, j, neg_inf);
            m.set(i, j, neg_inf);
        }

        // Slice the current (`i`) and previous (`i - 1`) rows of each matrix out
        // of the flat buffers ONCE per row, so the inner cell loop indexes into
        // fixed-window row slices -- the row base `i * bmax` is hoisted out of the
        // per-cell arithmetic instead of recomputing `i * bmax + j` and
        // bounds-checking against the entire matrix on every `Matrix::get`/`set`.
        // Each row slice is trimmed to `..=end_u`; since `j` ranges over
        // `start_u..=end_u` (and `j - 1 >= 0` because `start_u >= 1`), both `j`
        // and `j - 1` are provably in-bounds of a length-`(end_u + 1)` slice, so
        // LLVM can elide the slice bounds checks. This is byte-identical to the
        // `Matrix::get`/`set` form: identical reads, identical `max` chains,
        // identical writes -- only the access mechanics change.
        let bmax = m.bmax;
        let row = i * bmax;
        // `e` and `m` need both the previous and current row; `f` needs only the
        // current row (its recurrence reads `f[i][j-1]`, never row `i - 1`).
        let (e_head, e_tail) = e.data.split_at_mut(row);
        let e_prev = &e_head[row - bmax..][..=end_u];
        let e_cur = &mut e_tail[..=end_u];
        let (m_head, m_tail) = m.data.split_at_mut(row);
        let m_prev = &m_head[row - bmax..][..=end_u];
        let m_cur = &mut m_tail[..=end_u];
        let f_cur = &mut f.data[row..][..=end_u];

        for j in start_u..=end_u {
            // for e (insertion to the text): max(e[i-1][j]+GE, m[i-1][j]+GO+GE)
            let e_ij =
                (e_prev[j] + SCORE_GAPEXTEND).max(m_prev[j] + SCORE_GAPOPEN + SCORE_GAPEXTEND);
            e_cur[j] = e_ij;

            // for f (deletion to the text): max(f[i][j-1]+GE, m[i][j-1]+GO+GE)
            let f_ij = (f_cur[j - 1] + SCORE_GAPEXTEND)
                .max(m_cur[j - 1] + SCORE_GAPOPEN + SCORE_GAPEXTEND);
            f_cur[j] = f_ij;

            // for m: max(m[i-1][j-1]+match, e[i][j], f[i][j])
            let diag = m_prev[j - 1]
                + if chars_match(t[j - 1], p[i - 1]) { SCORE_MATCH } else { SCORE_MISMATCH };
            m_cur[j] = diag.max(e_ij).max(f_ij);
        }
    }

    let ret = m.get(lenp, lent);

    // Trace back: 3-state machine over which matrix ("m"=0, "e"=1, "f"=2)
    // the backtrace currently occupies. See module docs "THE critical trap"
    // for the exact tie-break order reproduced here.
    let mut tagi = lenp;
    let mut tagj = lent;
    let mut mat: u8 = 0;
    let mut align: Vec<i8> = Vec::new();

    while tagi > 0 || tagj > 0 {
        match mat {
            0 => {
                let max = e.get(tagi, tagj);
                let mut a = EDIT_INSERT;

                if f.get(tagi, tagj) >= max {
                    a = EDIT_DELETE;
                }
                if tagi > 0
                    && tagj > 0
                    && m.get(tagi - 1, tagj - 1)
                        + if chars_match(t[tagj - 1], p[tagi - 1]) {
                            SCORE_MATCH
                        } else {
                            SCORE_MISMATCH
                        }
                        == m.get(tagi, tagj)
                {
                    a = if chars_match(t[tagj - 1], p[tagi - 1]) {
                        EDIT_MATCH
                    } else {
                        EDIT_MISMATCH
                    };
                }

                if a == EDIT_MATCH || a == EDIT_MISMATCH {
                    align.push(a);
                    tagi -= 1;
                    tagj -= 1;
                } else if a == EDIT_INSERT {
                    mat = 1;
                } else if a == EDIT_DELETE {
                    mat = 2;
                }
            }
            1 => {
                // insertion to the text
                align.push(EDIT_INSERT);
                if tagi > 0 {
                    if m.get(tagi - 1, tagj) + SCORE_GAPOPEN + SCORE_GAPEXTEND == e.get(tagi, tagj)
                    {
                        tagi -= 1;
                        mat = 0;
                    } else {
                        tagi -= 1;
                        // mat stays 1
                    }
                } else {
                    mat = 2;
                }
            }
            _ => {
                // deletion to the text (mat == 2)
                align.push(EDIT_DELETE);
                if tagj > 0 {
                    if m.get(tagi, tagj - 1) + SCORE_GAPOPEN + SCORE_GAPEXTEND == f.get(tagi, tagj)
                    {
                        tagj -= 1;
                        mat = 0;
                    } else {
                        tagj -= 1;
                        // mat stays 2
                    }
                } else {
                    mat = 1;
                }
            }
        }
    }

    align.reverse();
    AlignResult { score: ret, align }
}

/// Ported from `SeqSet::GetAlignStats` (`SeqSet.hpp:438-453`): tallies
/// match/mismatch/indel counts from an `align[]` op sequence as produced by
/// [`global_alignment`]/[`global_alignment_pos_weight`].
///
/// If `update` is `false`, `match_cnt`/`mismatch_cnt`/`indel_cnt` are reset
/// to `0` before tallying (mirrors the C++ `if (!update) { matchCnt =
/// mismatchCnt = indelCnt = 0; }`); if `true`, the caller-supplied values are
/// accumulated onto (letting a caller sum stats across multiple `align[]`
/// sequences, matching `SeqSet.hpp`'s own `update=true` call sites at
/// `SeqSet.hpp:2057`, `VariantCaller.hpp`).
pub fn get_align_stats(
    align: &[i8],
    update: bool,
    match_cnt: &mut i32,
    mismatch_cnt: &mut i32,
    indel_cnt: &mut i32,
) {
    let (m, mm, ind) = stats_of(align);
    if update {
        *match_cnt += m;
        *mismatch_cnt += mm;
        *indel_cnt += ind;
    } else {
        *match_cnt = m;
        *mismatch_cnt = mm;
        *indel_cnt = ind;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pw(a: i32, c: i32, g: i32, t: i32) -> PosWeight {
        PosWeight { count: [a, c, g, t] }
    }

    // ---- is_base_equal --------------------------------------------------

    #[test]
    fn is_base_equal_empty_column_is_always_equal() {
        assert!(is_base_equal(&pw(0, 0, 0, 0), b'A'));
        assert!(is_base_equal(&pw(0, 0, 0, 0), b'T'));
    }

    #[test]
    fn is_base_equal_n_is_always_equal() {
        assert!(is_base_equal(&pw(10, 0, 0, 0), b'N'));
    }

    #[test]
    fn is_base_equal_dominant_base_is_equal() {
        // count[A] = 8, sum = 10 -> 10 < 3*8=24 -> true (dominant call).
        assert!(is_base_equal(&pw(8, 1, 1, 0), b'A'));
    }

    #[test]
    fn is_base_equal_minority_base_is_not_equal() {
        // count[T]=1, sum=10 -> 10 < 3*1=3? no -> false.
        assert!(!is_base_equal(&pw(8, 1, 1, 0), b'T'));
    }

    #[test]
    fn is_base_equal_unanimous_single_base_is_equal() {
        // count[A]=5, sum=5 -> 5 < 15 -> true.
        assert!(is_base_equal(&pw(5, 0, 0, 0), b'A'));
    }

    // ---- global_alignment: hand-computed cases ---------------------------

    #[test]
    fn global_alignment_empty_inputs() {
        let r = global_alignment(b"", b"ACGT", DEFAULT_BAND);
        assert_eq!(r.score, 0);
        assert!(r.align.is_empty());
        let r = global_alignment(b"ACGT", b"", DEFAULT_BAND);
        assert_eq!(r.score, 0);
        assert!(r.align.is_empty());
    }

    #[test]
    fn global_alignment_single_base_match() {
        let r = global_alignment(b"A", b"A", DEFAULT_BAND);
        assert_eq!(r.score, SCORE_MATCH);
        assert_eq!(r.align, vec![EDIT_MATCH]);
    }

    #[test]
    fn global_alignment_single_base_n_wildcard_match() {
        let r = global_alignment(b"N", b"A", DEFAULT_BAND);
        assert_eq!(r.score, SCORE_MATCH);
        assert_eq!(r.align, vec![EDIT_MATCH]);
        let r = global_alignment(b"A", b"N", DEFAULT_BAND);
        assert_eq!(r.score, SCORE_MATCH);
        assert_eq!(r.align, vec![EDIT_MATCH]);
    }

    #[test]
    fn global_alignment_single_base_mismatch() {
        let r = global_alignment(b"A", b"C", DEFAULT_BAND);
        assert_eq!(r.score, SCORE_MISMATCH);
        assert_eq!(r.align, vec![EDIT_MISMATCH]);
    }

    #[test]
    fn global_alignment_exact_match_multi_base() {
        let r = global_alignment(b"ACGTACGT", b"ACGTACGT", DEFAULT_BAND);
        assert_eq!(r.score, 8 * SCORE_MATCH);
        assert_eq!(r.align, vec![EDIT_MATCH; 8]);
    }

    #[test]
    fn global_alignment_single_mismatch_in_middle() {
        let r = global_alignment(b"ACGTACGT", b"ACGAACGT", DEFAULT_BAND);
        assert_eq!(r.score, 7 * SCORE_MATCH + SCORE_MISMATCH);
        let mut expected = vec![EDIT_MATCH; 8];
        expected[3] = EDIT_MISMATCH;
        assert_eq!(r.align, expected);
    }

    #[test]
    fn global_alignment_single_insertion() {
        // p has one extra base relative to t: t="ACGTACGT" (8), p="ACGTXACGT" (9).
        let r = global_alignment(b"ACGTACGT", b"ACGTAACGT", DEFAULT_BAND);
        // Best alignment: 8 matches + 1 inserted base (gap open+extend).
        assert_eq!(r.score, 8 * SCORE_MATCH + SCORE_GAPOPEN + SCORE_GAPEXTEND);
        let inserts = r.align.iter().filter(|&&a| a == EDIT_INSERT).count();
        let matches = r.align.iter().filter(|&&a| a == EDIT_MATCH).count();
        assert_eq!(inserts, 1);
        assert_eq!(matches, 8);
    }

    #[test]
    fn global_alignment_single_deletion() {
        // t has one extra base relative to p.
        let r = global_alignment(b"ACGTAACGT", b"ACGTACGT", DEFAULT_BAND);
        assert_eq!(r.score, 8 * SCORE_MATCH + SCORE_GAPOPEN + SCORE_GAPEXTEND);
        let deletes = r.align.iter().filter(|&&a| a == EDIT_DELETE).count();
        let matches = r.align.iter().filter(|&&a| a == EDIT_MATCH).count();
        assert_eq!(deletes, 1);
        assert_eq!(matches, 8);
    }

    #[test]
    fn global_alignment_align_ops_consume_expected_lengths() {
        // Sanity: number of non-insert ops == lent, number of non-delete ops == lenp.
        let t = b"ACGTACGTACGT";
        let p = b"ACGTAACGTACC";
        let r = global_alignment(t, p, DEFAULT_BAND);
        let t_consumed = r.align.iter().filter(|&&a| a != EDIT_INSERT).count();
        let p_consumed = r.align.iter().filter(|&&a| a != EDIT_DELETE).count();
        assert_eq!(t_consumed, t.len());
        assert_eq!(p_consumed, p.len());
    }

    // ---- global_alignment: diagonal fast path vs full DP (regression) -----
    //
    // These assert the fast-path (0/1/2 mismatches, equal length) returns
    // exactly the hand-computed diagonal score/align -- the same values the
    // full DP is required to produce per the byte-identity proof on
    // `global_alignment`'s doc comment. A 3-mismatch case is included to
    // prove the shortcut correctly declines (falls through to the DP) and
    // that the DP path still returns the correct result independently.

    #[test]
    fn global_alignment_fast_path_zero_mismatches() {
        let t = b"ACGTACGTACGT";
        let p = b"ACGTACGTACGT";
        let r = global_alignment(t, p, DEFAULT_BAND);
        assert_eq!(r.score, 12 * SCORE_MATCH);
        assert_eq!(r.align, vec![EDIT_MATCH; 12]);
    }

    #[test]
    fn global_alignment_fast_path_one_mismatch() {
        let t = b"ACGTACGTACGT";
        let p = b"ACCTACGTACGT"; // mismatch at position 2
        let r = global_alignment(t, p, DEFAULT_BAND);
        assert_eq!(r.score, 11 * SCORE_MATCH + SCORE_MISMATCH);
        let mut expected = vec![EDIT_MATCH; 12];
        expected[2] = EDIT_MISMATCH;
        assert_eq!(r.align, expected);
    }

    #[test]
    fn global_alignment_fast_path_two_mismatches() {
        let t = b"ACGTACGTACGT";
        let p = b"ACCTACCTACGT"; // mismatches at positions 2, 6
        let r = global_alignment(t, p, DEFAULT_BAND);
        assert_eq!(r.score, 10 * SCORE_MATCH + 2 * SCORE_MISMATCH);
        let mut expected = vec![EDIT_MATCH; 12];
        expected[2] = EDIT_MISMATCH;
        expected[6] = EDIT_MISMATCH;
        assert_eq!(r.align, expected);
    }

    #[test]
    fn global_alignment_fast_path_two_mismatches_with_n_wildcard() {
        // `N` on either side is a wildcard match (chars_match), so it must
        // NOT count toward the mismatch tally, and must score as a match.
        let t = b"ACGTACGTACGT";
        let p = b"ACCTACNTACGT"; // real mismatch at 2, N-wildcard "match" at 6
        let r = global_alignment(t, p, DEFAULT_BAND);
        assert_eq!(r.score, 11 * SCORE_MATCH + SCORE_MISMATCH);
        let mut expected = vec![EDIT_MATCH; 12];
        expected[2] = EDIT_MISMATCH;
        assert_eq!(r.align, expected);

        // Two `N`-wildcard positions (one in `t`, one in `p`, free matches)
        // PLUS two real mismatches: still only 2 scoring mismatches, so the
        // fast path must still trigger (4 differing positions total, but
        // only 2 count toward the `<= 2` mismatch tally).
        //   t2: A C N T A C G T A C C T
        //   p2: A C C T A C N T A A A T
        //       .  .  N  .  .  .  N  .  .  X  X  .   (N = wildcard match, X = real mismatch)
        let t2 = b"ACNTACGTACCT";
        let p2 = b"ACCTACNTAAAT";
        let r2 = global_alignment(t2, p2, DEFAULT_BAND);
        assert_eq!(r2.score, 10 * SCORE_MATCH + 2 * SCORE_MISMATCH);
        let mut expected2 = vec![EDIT_MATCH; 12];
        expected2[9] = EDIT_MISMATCH;
        expected2[10] = EDIT_MISMATCH;
        assert_eq!(r2.align, expected2);
    }

    #[test]
    fn global_alignment_three_mismatches_falls_through_to_dp() {
        // 3 mismatches (scattered, non-adjacent) is one past the fast
        // path's `<= 2` cutoff, so this must fall through to the full banded
        // DP below. Verified independently (see PR description / commit
        // message) that the DP still resolves this to the pure diagonal --
        // i.e. the fast path's `<= 2` threshold is conservative, not tight,
        // but this test locks in that the DP path alone (fast path
        // deliberately not taken) still returns the correct score/align.
        let t = b"ACGTACGTACGT";
        let p = b"ACCTACCTACCT"; // mismatches at positions 2, 6, 10
        let r = global_alignment(t, p, DEFAULT_BAND);
        assert_eq!(r.score, 9 * SCORE_MATCH + 3 * SCORE_MISMATCH);
        let mut expected = vec![EDIT_MATCH; 12];
        expected[2] = EDIT_MISMATCH;
        expected[6] = EDIT_MISMATCH;
        expected[10] = EDIT_MISMATCH;
        assert_eq!(r.align, expected);
    }

    #[test]
    fn diagonal_fast_path_cutoff_is_locked_at_two_mismatches() {
        // The `<= 2` cutoff cannot be observed through `global_alignment`'s
        // OUTPUT (the fast path and full DP are byte-identical for every
        // equal-length input at this many mismatches -- the bound is
        // conservative, so a loosened cutoff would still return the same
        // score/traceback). This test locks the literal contract by asserting
        // the eligibility predicate directly at its boundary.

        // Eligible: equal length, 0/1/2 mismatches.
        assert!(diagonal_fast_path_applies(b"ACGT", b"ACGT")); // 0 mismatches
        assert!(diagonal_fast_path_applies(b"ACGT", b"AGGT")); // 1 mismatch
        assert!(diagonal_fast_path_applies(b"ACGT", b"AGGA")); // 2 mismatches

        // NOT eligible: 3 mismatches is one past the cutoff -- must fall
        // through to the DP. (Same input as
        // `global_alignment_three_mismatches_falls_through_to_dp`.)
        assert!(!diagonal_fast_path_applies(b"ACGTACGTACGT", b"ACCTACCTACCT")); // 3 mismatches

        // NOT eligible: unequal length never takes the fast path, regardless
        // of mismatch count.
        assert!(!diagonal_fast_path_applies(b"ACGT", b"ACG"));
        assert!(!diagonal_fast_path_applies(b"ACG", b"ACGT"));
    }

    // ---- global_alignment_pos_weight: hand-computed cases -----------------

    #[test]
    fn pos_weight_empty_inputs() {
        let r = global_alignment_pos_weight(&[], b"ACGT");
        assert_eq!(r.score, 0);
        assert!(r.align.is_empty());
    }

    #[test]
    fn pos_weight_single_base_match_via_dominant_call() {
        let r = global_alignment_pos_weight(&[pw(5, 0, 0, 0)], b"A");
        assert_eq!(r.score, SCORE_MATCH);
        assert_eq!(r.align, vec![EDIT_MATCH]);
    }

    #[test]
    fn pos_weight_single_base_mismatch() {
        let r = global_alignment_pos_weight(&[pw(5, 0, 0, 0)], b"T");
        assert_eq!(r.score, SCORE_MISMATCH);
        assert_eq!(r.align, vec![EDIT_MISMATCH]);
    }

    #[test]
    fn pos_weight_exact_match_no_indel_shortcut() {
        let weights: Vec<PosWeight> = b"ACGTACGT"
            .iter()
            .map(|&c| match c {
                b'A' => pw(5, 0, 0, 0),
                b'C' => pw(0, 5, 0, 0),
                b'G' => pw(0, 0, 5, 0),
                b'T' => pw(0, 0, 0, 5),
                _ => unreachable!(),
            })
            .collect();
        let r = global_alignment_pos_weight(&weights, b"ACGTACGT");
        assert_eq!(r.score, 8 * SCORE_MATCH);
        assert_eq!(r.align, vec![EDIT_MATCH; 8]);
    }

    // ---- regression: Critical fix #1 (row-0 `e` init off-by-one) ----------

    /// `t` = 9x'G', `p` = "G", band = 5. Verified against the real C++
    /// oracle (`crates/unum-core/tests/golden_align_algo.rs`'s
    /// `regression_row0_e_init_off_by_one_lenp_plus_one`): score=-10,
    /// align=[MATCH, DELETE x8]. Before the fix, this Rust port's row-0 `e`
    /// boundary cell held the wrong value (`lenp * SCORE_GAPOPEN` instead of
    /// `(lenp + 1) * SCORE_GAPOPEN`), which the wide-band traceback reads at
    /// `tagi == 0` here, producing an incorrect score of -9.
    #[test]
    fn regression_row0_e_init_off_by_one_lenp_plus_one() {
        let t = b"GGGGGGGGG";
        let p = b"G";
        let r = global_alignment(t, p, DEFAULT_BAND);
        assert_eq!(r.score, -10, "score");
        let mut expected = vec![EDIT_DELETE; 9];
        expected[0] = EDIT_MATCH;
        assert_eq!(r.align, expected, "align[]");
    }

    // ---- regression: Critical fix #2 (traceback usize underflow panic) ----

    /// `lent=14, lenp=3, p="CGC"`, with the posweight matrix below. Frozen in
    /// the golden test (`crates/unum-core/tests/golden_align_algo.rs`'s
    /// `regression_traceback_usize_underflow_panic`):
    /// score=-42, align=[MATCH, DELETE x10, MATCH, MATCH, MATCH]. Before the
    /// fix, this Rust port's `usize` traceback cursors panicked with
    /// "attempt to subtract with overflow" on this input (the sentinel `a ==
    /// EDIT_MATCH` default survives at `tagj == 0, tagi > 0`, and the
    /// diagonal `else` arm decrements `tagj` below zero -- exactly what
    /// C++'s signed `int` cursors tolerate but a `usize` cannot).
    #[test]
    fn regression_traceback_usize_underflow_panic() {
        let counts: [[i32; 4]; 14] = [
            [0, 0, 0, 10],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 10],
            [0, 0, 10, 0],
            [5, 5, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 10],
            [0, 10, 0, 0],
            [0, 0, 5, 5],
            [5, 5, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 5, 5],
            [0, 10, 0, 0],
        ];
        let weights: Vec<PosWeight> =
            counts.iter().map(|&count| pw(count[0], count[1], count[2], count[3])).collect();
        let p = b"CGC";

        let r = global_alignment_pos_weight(&weights, p); // must not panic
        assert_eq!(r.score, -42, "score");
        let mut expected = vec![EDIT_DELETE; 14];
        expected[0] = EDIT_MATCH;
        expected[11] = EDIT_MATCH;
        expected[12] = EDIT_MATCH;
        expected[13] = EDIT_MATCH;
        assert_eq!(r.align, expected, "align[]");
    }

    // ---- get_align_stats --------------------------------------------------

    #[test]
    fn get_align_stats_counts_each_op_kind() {
        let align = [EDIT_MATCH, EDIT_MATCH, EDIT_MISMATCH, EDIT_INSERT, EDIT_DELETE, EDIT_MATCH];
        let (mut m, mut mm, mut i) = (0, 0, 0);
        get_align_stats(&align, false, &mut m, &mut mm, &mut i);
        assert_eq!((m, mm, i), (3, 1, 2));
    }

    #[test]
    fn get_align_stats_update_true_accumulates() {
        let align = [EDIT_MATCH, EDIT_MISMATCH];
        let (mut m, mut mm, mut i) = (10, 20, 30);
        get_align_stats(&align, true, &mut m, &mut mm, &mut i);
        assert_eq!((m, mm, i), (11, 21, 30));
    }

    #[test]
    fn get_align_stats_update_false_resets_before_counting() {
        let align = [EDIT_MATCH];
        let (mut m, mut mm, mut i) = (10, 20, 30);
        get_align_stats(&align, false, &mut m, &mut mm, &mut i);
        assert_eq!((m, mm, i), (1, 0, 0));
    }

    #[test]
    fn get_align_stats_empty_align_is_noop_when_updating() {
        let (mut m, mut mm, mut i) = (1, 2, 3);
        get_align_stats(&[], true, &mut m, &mut mm, &mut i);
        assert_eq!((m, mm, i), (1, 2, 3));
    }

    // ---- randomized differential fuzz of the banded-DP fill refactor ------
    //
    // The `global_alignment_dp_with` inner fill loop was refactored from
    // per-cell `Matrix::get`/`set` (flat `i * bmax + j` indexing) to per-row
    // slice windows. These tests fuzz that refactor against (a) an independent
    // naive reference DP written below with plain 2D `Vec` matrices (byte-for-
    // byte comparison of both `score` AND `align[]`), and (b) a self-consistency
    // check that the returned `align[]` op sequence, scored under the affine
    // scoring scheme, reproduces the returned `score`. All randomness is a
    // deterministic seeded LCG -- no system randomness or time.

    /// Deterministic 64-bit linear-congruential PRNG (SplitMix-style output
    /// mix) so the fuzz corpus is reproducible across runs and machines.
    struct Lcg {
        state: u64,
    }

    impl Lcg {
        fn new(seed: u64) -> Self {
            Self { state: seed }
        }

        fn next_u64(&mut self) -> u64 {
            // LCG step (Numerical Recipes constants), then a SplitMix64 output
            // finalizer for well-distributed low bits.
            self.state = self.state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            let mut z = self.state;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        }

        /// Uniform-ish integer in `[0, n)` for small `n` (fuzz corpus, not
        /// cryptographic uniformity).
        fn below(&mut self, n: usize) -> usize {
            (self.next_u64() % (n as u64)) as usize
        }
    }

    /// A random `ACGTN` byte string of length `len`.
    fn random_seq(rng: &mut Lcg, len: usize) -> Vec<u8> {
        const ALPHABET: &[u8; 5] = b"ACGTN";
        (0..len).map(|_| ALPHABET[rng.below(ALPHABET.len())]).collect()
    }

    /// Independent naive reference for [`global_alignment`]: the SAME banded
    /// affine-gap Gotoh DP, but transcribed with plain `Vec<Vec<i32>>` matrices
    /// and per-cell `[i][j]` indexing (no flat buffer, no row slicing) so it is
    /// a genuine cross-check of the sliced fill's matrix values. Reproduces the
    /// full public ladder (empty / single-base / diagonal fast path / banded
    /// DP) and the exact 3-state traceback, so it must be byte-identical to
    /// `global_alignment` for every input.
    fn naive_global_alignment(t: &[u8], p: &[u8], band: i32) -> AlignResult {
        let lent = t.len();
        let lenp = p.len();

        if lent == 0 || lenp == 0 {
            return AlignResult { score: 0, align: Vec::new() };
        }
        if lent == 1 && lenp == 1 {
            return if chars_match(t[0], p[0]) {
                AlignResult { score: SCORE_MATCH, align: vec![EDIT_MATCH] }
            } else {
                AlignResult { score: SCORE_MISMATCH, align: vec![EDIT_MISMATCH] }
            };
        }
        // Diagonal fast path (equal length, <= 2 mismatches).
        if lent == lenp {
            let mismatches = (0..lent).filter(|&k| !chars_match(t[k], p[k])).count();
            if mismatches <= DIAGONAL_FAST_PATH_MAX_MISMATCHES as usize {
                let mut align = Vec::with_capacity(lent);
                let mut score = 0;
                for k in 0..lent {
                    if chars_match(t[k], p[k]) {
                        align.push(EDIT_MATCH);
                        score += SCORE_MATCH;
                    } else {
                        align.push(EDIT_MISMATCH);
                        score += SCORE_MISMATCH;
                    }
                }
                return AlignResult { score, align };
            }
        }

        let (lent_i32, lenp_i32) = (lent as i32, lenp as i32);
        let (left_band, right_band) = if lent_i32 > lenp_i32 {
            (band, band + (lent_i32 - lenp_i32))
        } else if lent_i32 < lenp_i32 {
            (band + (lenp_i32 - lent_i32), band)
        } else {
            (band, band)
        };
        let neg_inf = (lent_i32 + 1) * (lenp_i32 + 1) * SCORE_GAPOPEN;

        let mut m = vec![vec![0i32; lent + 1]; lenp + 1];
        let mut e = vec![vec![0i32; lent + 1]; lenp + 1];
        let mut f = vec![vec![0i32; lent + 1]; lenp + 1];

        for i in 1..=lenp {
            let i_i32 = i as i32;
            e[i][0] = SCORE_GAPOPEN + i_i32 * SCORE_GAPEXTEND;
            f[i][0] = SCORE_GAPOPEN + i_i32 * SCORE_GAPOPEN;
            m[i][0] = SCORE_GAPOPEN + i_i32 * SCORE_GAPOPEN;
        }
        let row0_e_sentinel = SCORE_GAPOPEN + (lenp_i32 + 1) * SCORE_GAPOPEN;
        for j in 1..=lent {
            let j_i32 = j as i32;
            f[0][j] = SCORE_GAPOPEN + j_i32 * SCORE_GAPEXTEND;
            e[0][j] = row0_e_sentinel;
            m[0][j] = SCORE_GAPOPEN + j_i32 * SCORE_GAPOPEN;
        }

        for i in 1..=lenp {
            let i_i32 = i as i32;
            let start = if i_i32 - left_band < 1 { 1 } else { i_i32 - left_band };
            let end = if i_i32 + right_band > lent_i32 { lent_i32 } else { i_i32 + right_band };
            let (start_u, end_u) = (start as usize, end as usize);

            if start > 1 {
                let j = start_u - 1;
                e[i][j] = neg_inf;
                f[i][j] = neg_inf;
                m[i][j] = neg_inf;
            }
            if end < lent_i32 {
                let j = end_u + 1;
                e[i][j] = neg_inf;
                f[i][j] = neg_inf;
                m[i][j] = neg_inf;
            }

            for j in start_u..=end_u {
                let es = (e[i - 1][j] + SCORE_GAPEXTEND)
                    .max(m[i - 1][j] + SCORE_GAPOPEN + SCORE_GAPEXTEND);
                e[i][j] = es;
                let fs = (f[i][j - 1] + SCORE_GAPEXTEND)
                    .max(m[i][j - 1] + SCORE_GAPOPEN + SCORE_GAPEXTEND);
                f[i][j] = fs;
                let diag = m[i - 1][j - 1]
                    + if chars_match(t[j - 1], p[i - 1]) { SCORE_MATCH } else { SCORE_MISMATCH };
                m[i][j] = diag.max(es).max(fs);
            }
        }

        let ret = m[lenp][lent];

        let mut tagi = lenp;
        let mut tagj = lent;
        let mut mat: u8 = 0;
        let mut align: Vec<i8> = Vec::new();
        while tagi > 0 || tagj > 0 {
            match mat {
                0 => {
                    let max = e[tagi][tagj];
                    let mut a = EDIT_INSERT;
                    if f[tagi][tagj] >= max {
                        a = EDIT_DELETE;
                    }
                    if tagi > 0
                        && tagj > 0
                        && m[tagi - 1][tagj - 1]
                            + if chars_match(t[tagj - 1], p[tagi - 1]) {
                                SCORE_MATCH
                            } else {
                                SCORE_MISMATCH
                            }
                            == m[tagi][tagj]
                    {
                        a = if chars_match(t[tagj - 1], p[tagi - 1]) {
                            EDIT_MATCH
                        } else {
                            EDIT_MISMATCH
                        };
                    }
                    if a == EDIT_MATCH || a == EDIT_MISMATCH {
                        align.push(a);
                        tagi -= 1;
                        tagj -= 1;
                    } else if a == EDIT_INSERT {
                        mat = 1;
                    } else if a == EDIT_DELETE {
                        mat = 2;
                    }
                }
                1 => {
                    align.push(EDIT_INSERT);
                    if tagi > 0 {
                        if m[tagi - 1][tagj] + SCORE_GAPOPEN + SCORE_GAPEXTEND == e[tagi][tagj] {
                            tagi -= 1;
                            mat = 0;
                        } else {
                            tagi -= 1;
                        }
                    } else {
                        mat = 2;
                    }
                }
                _ => {
                    align.push(EDIT_DELETE);
                    if tagj > 0 {
                        if m[tagi][tagj - 1] + SCORE_GAPOPEN + SCORE_GAPEXTEND == f[tagi][tagj] {
                            tagj -= 1;
                            mat = 0;
                        } else {
                            tagj -= 1;
                        }
                    } else {
                        mat = 1;
                    }
                }
            }
        }
        align.reverse();
        AlignResult { score: ret, align }
    }

    #[test]
    fn fuzz_dp_matches_naive_reference() {
        let mut rng = Lcg::new(0x0DDC_0FFE_E123_4567);
        let mut dp_cases = 0u32; // count inputs that actually exercised the banded DP
        for _ in 0..20_000 {
            // Bias toward inputs that reach the banded DP: mostly unequal
            // lengths (or long enough to accumulate >2 mismatches), plus some
            // equal-length fast-path cases for full-ladder coverage.
            let lent = rng.below(81); // 0..=80
            let lenp = rng.below(81);
            let band = 1 + rng.below(8) as i32; // 1..=8
            let t = random_seq(&mut rng, lent);
            let p = random_seq(&mut rng, lenp);

            let got = global_alignment(&t, &p, band);
            let want = naive_global_alignment(&t, &p, band);
            assert_eq!(got.score, want.score, "score mismatch for t={t:?} p={p:?} band={band}");
            assert_eq!(got.align, want.align, "align mismatch for t={t:?} p={p:?} band={band}");

            // Cached and stats-only entry points share the refactored DP; they
            // must agree with the direct call.
            let mut cache = DpCache::new();
            let cached = global_alignment_cached(&t, &p, band, &mut cache);
            assert_eq!(cached, got, "cached != uncached for t={t:?} p={p:?} band={band}");
            let (m0, mm0, i0) = {
                let (mut a, mut b, mut c) = (0, 0, 0);
                get_align_stats(&got.align, false, &mut a, &mut b, &mut c);
                (a, b, c)
            };
            let mut cache2 = DpCache::new();
            let stats = global_alignment_cached_stats(&t, &p, band, &mut cache2);
            assert_eq!(stats, (m0, mm0, i0), "stats mismatch for t={t:?} p={p:?} band={band}");

            // Track how many cases fell through the fast paths into the DP so
            // the fuzz corpus is demonstrably exercising the refactored kernel.
            let is_trivial = lent == 0
                || lenp == 0
                || (lent == 1 && lenp == 1)
                || (lent == lenp
                    && (0..lent).filter(|&k| !chars_match(t[k], p[k])).count()
                        <= DIAGONAL_FAST_PATH_MAX_MISMATCHES as usize);
            if !is_trivial {
                dp_cases += 1;
            }
        }
        // Guard against the corpus silently degenerating into only fast-path
        // inputs (which would leave the refactored DP fill untested).
        assert!(dp_cases > 1000, "fuzz corpus exercised the banded DP only {dp_cases} times");
    }
}
