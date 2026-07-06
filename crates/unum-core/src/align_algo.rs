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

/// 2D DP matrix stored as a flat `Vec<i32>`, row-major with `bmax = lent + 1`
/// columns, mirroring the C++ side's manual `m[i * bmax + j]` flat indexing
/// exactly (kept flat, not `Vec<Vec<i32>>`, so index arithmetic is visibly
/// identical to the ported C++ at every call site).
struct Matrix {
    data: Vec<i32>,
    bmax: usize,
}

impl Matrix {
    fn new(lenp: usize, lent: usize) -> Self {
        Self { data: vec![0; (lenp + 1) * (lent + 1)], bmax: lent + 1 }
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
/// # Panics
///
/// Panics if `t.len()` or `p.len()` do not fit an `i32`, or if `band` is
/// negative (not a data-dependent condition; every real caller passes
/// [`DEFAULT_BAND`]).
#[must_use]
pub fn global_alignment(t: &[u8], p: &[u8], band: i32) -> AlignResult {
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

    let mut m = Matrix::new(lenp, lent);
    let mut e = Matrix::new(lenp, lent);
    let mut f = Matrix::new(lenp, lent);

    m.set(0, 0, 0);
    e.set(0, 0, 0);
    f.set(0, 0, 0);

    // Mirrors AlignAlgo.hpp:256-266 verbatim, including the row-0 `e`
    // initialization quirk: `e[0+j] = SCORE_GAPOPEN + i * SCORE_GAPOPEN` reads
    // the LOOP VARIABLE `i` LEFT OVER from the immediately preceding `for (i
    // = 1; i <= lenp; ++i)` loop, NOT the current `j`. In C++, a `for` loop's
    // control variable is NOT scoped to the loop body when declared outside
    // it (as it is here -- `i`/`j` are declared once, above both loops), so
    // after that loop exits its condition check fails at `i == lenp + 1`,
    // leaving `i == lenp + 1` for the second loop to read. So every `e[0][j]`
    // cell (all `j`) is set to the SAME value, `SCORE_GAPOPEN + (lenp + 1) *
    // SCORE_GAPOPEN`. This cell is NOT dead: the `mat == 0` traceback reads
    // `max = e[tagi * bmax + tagj]` unconditionally on entry to that state
    // (AlignAlgo.hpp:332), including when `tagi == 0` -- and widening the
    // band (e.g. `lent` much greater than `lenp`) can steer the traceback
    // into visiting `e` at `tagi == 0` before the loop terminates on `tagj ==
    // 0` too. It is reproduced exactly (leftover-loop-index value, not the
    // "intended" `j`-indexed value) per this port's byte-identity mandate.
    for i in 1..=lenp {
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let i_i32 = i as i32;
        e.set(i, 0, SCORE_GAPOPEN + i_i32 * SCORE_GAPEXTEND);
        f.set(i, 0, SCORE_GAPOPEN + i_i32 * SCORE_GAPOPEN);
        m.set(i, 0, SCORE_GAPOPEN + i_i32 * SCORE_GAPOPEN);
    }
    for j in 1..=lent {
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let j_i32 = j as i32;
        f.set(0, j, SCORE_GAPOPEN + j_i32 * SCORE_GAPEXTEND);
        e.set(0, j, SCORE_GAPOPEN + (lenp_i32 + 1) * SCORE_GAPOPEN);
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

        for j in start_u..=end_u {
            // for e (insertion to the text)
            let mut score = e.get(i - 1, j) + SCORE_GAPEXTEND;
            score = score.max(m.get(i - 1, j) + SCORE_GAPOPEN + SCORE_GAPEXTEND);
            e.set(i, j, score);

            // for f (deletion to the text)
            let mut score = f.get(i, j - 1) + SCORE_GAPEXTEND;
            score = score.max(m.get(i, j - 1) + SCORE_GAPOPEN + SCORE_GAPEXTEND);
            f.set(i, j, score);

            // for m
            let mut score = m.get(i - 1, j - 1)
                + if chars_match(t[j - 1], p[i - 1]) { SCORE_MATCH } else { SCORE_MISMATCH };
            score = score.max(e.get(i, j));
            score = score.max(f.get(i, j));
            m.set(i, j, score);
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
    if !update {
        *match_cnt = 0;
        *mismatch_cnt = 0;
        *indel_cnt = 0;
    }

    for &op in align {
        if op == EDIT_MATCH {
            *match_cnt += 1;
        } else if op == EDIT_MISMATCH {
            *mismatch_cnt += 1;
        } else {
            *indel_cnt += 1;
        }
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
}
