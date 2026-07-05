#![cfg(feature = "t1k-sys")]
//! Differential test: compares the Rust port of `AlignAlgo` in
//! `fg_t1k_core::align_algo` against the real, unmodified vendored C++
//! `AlignAlgo::GlobalAlignment`/`GlobalAlignment_PosWeight` (plus
//! `SeqSet::GetAlignStats`), via the free-function FFI shim
//! (`fg_t1k_sys::CppAlignAlgo`).
//!
//! # Why compare the FULL `align[]` array, not just the score
//!
//! Because scoring is all-integer, the DP score alone is comparatively easy
//! to get byte-identical -- the recurrences are `MAX(...)` chains, and tie
//! order never affects the numeric max. The traceback tie-break order,
//! however, decides WHICH of several equally-scoring alignments is emitted
//! as the op sequence. A scoring-only differential would happily pass a
//! Rust port whose traceback silently picks a different (but equally-scored)
//! alignment path than the C++ oracle -- e.g. reporting a run of 3
//! consecutive single-base mismatches where the real T1K would report a
//! 1-base insertion + 1-base deletion, changing downstream indel/mismatch
//! counts and every consumer that inspects individual ops. Every test below
//! therefore asserts `(score, align)` equality, not just `score` equality.
//!
//! `flip_tie_break_and_confirm_align_array_catches_it` (at the bottom) is
//! the falsifiability check the brief calls for: it patches a COPY of the
//! Rust traceback logic with one tie-break flipped, confirms score is
//! unaffected but `align[]` diverges from the real oracle on a case
//! specifically chosen to exercise a genuine multi-way score tie -- proving
//! the array-level assertion actually has teeth.

use fg_t1k_core::align_algo::{
    DEFAULT_BAND, EDIT_DELETE, EDIT_MATCH, PosWeight, get_align_stats, global_alignment,
    global_alignment_pos_weight,
};
use fg_t1k_sys::CppAlignAlgo;

/// Runs both the Rust port and the C++ oracle on `(t, p)` at the given
/// `band`, asserting `score` AND `align` (the full op sequence) are
/// identical, element-for-element.
fn assert_global_alignment_matches(label: &str, t: &[u8], p: &[u8], band: i32) {
    let rust = global_alignment(t, p, band);
    let cpp = CppAlignAlgo::global_alignment(t, p, band);

    assert_eq!(
        rust.score,
        cpp.score,
        "{label}: score mismatch (t={:?}, p={:?}, band={band}, rust_align={:?}, cpp_align={:?})",
        String::from_utf8_lossy(t),
        String::from_utf8_lossy(p),
        rust.align,
        cpp.align,
    );
    assert_eq!(
        rust.align,
        cpp.align,
        "{label}: align[] mismatch (t={:?}, p={:?}, band={band})",
        String::from_utf8_lossy(t),
        String::from_utf8_lossy(p),
    );
}

// ---- exact match, mismatches, indels, mixed -------------------------------

#[test]
fn exact_match() {
    assert_global_alignment_matches(
        "exact_match",
        b"ACGTACGTACGTACGT",
        b"ACGTACGTACGTACGT",
        DEFAULT_BAND,
    );
}

#[test]
fn single_mismatch() {
    assert_global_alignment_matches(
        "single_mismatch",
        b"ACGTACGTACGTACGT",
        b"ACGTAAGTACGTACGT",
        DEFAULT_BAND,
    );
}

#[test]
fn multiple_mismatches() {
    assert_global_alignment_matches(
        "multiple_mismatches",
        b"ACGTACGTACGTACGTACGT",
        b"ACCTACGAACGTACCTACGT",
        DEFAULT_BAND,
    );
}

#[test]
fn single_insertion() {
    // p has one extra base relative to t.
    assert_global_alignment_matches(
        "single_insertion",
        b"ACGTACGTACGTACGT",
        b"ACGTACCGTACGTACGTACGT"[..].as_ref(),
        DEFAULT_BAND,
    );
}

#[test]
fn single_deletion() {
    // t has one extra base relative to p.
    assert_global_alignment_matches(
        "single_deletion",
        b"ACGTACCGTACGTACGTACGT",
        b"ACGTACGTACGTACGT",
        DEFAULT_BAND,
    );
}

#[test]
fn mixed_indels_and_mismatches() {
    assert_global_alignment_matches(
        "mixed_indels_and_mismatches",
        b"ACGTACGTTTACGTACGTACGGGTACGT",
        b"ACGAACGTACGTACGCACGTACGT",
        DEFAULT_BAND,
    );
}

#[test]
fn n_wildcard_bases() {
    assert_global_alignment_matches(
        "n_wildcard_bases",
        b"ACGTNACGTACGT",
        b"ACGTAACGTNCGT",
        DEFAULT_BAND,
    );
}

// ---- band-edge cases --------------------------------------------------

#[test]
fn indel_at_band_edge_small_band() {
    // With band=2, an indel exactly at the edge of the allowed band forces
    // the DP to touch its negInf sentinel boundary cells.
    assert_global_alignment_matches(
        "indel_at_band_edge_small_band",
        b"ACGTACGTAACGTACGTACGTACGT",
        b"ACGTACGTACGTACGTACGTACGT",
        2,
    );
}

#[test]
fn indel_at_band_edge_band_one() {
    assert_global_alignment_matches(
        "indel_at_band_edge_band_one",
        b"AAAAAAAAAAAAAAAA",
        b"AAAAAAAAAAAAAAAAA",
        1,
    );
}

#[test]
fn differing_lengths_widens_band_left() {
    // lenp > lent: leftBand widens by (lenp - lent).
    assert_global_alignment_matches(
        "differing_lengths_widens_band_left",
        b"ACGTACGT",
        b"ACGTACGTACGTACGTACGT",
        DEFAULT_BAND,
    );
}

#[test]
fn differing_lengths_widens_band_right() {
    // lent > lenp: rightBand widens by (lent - lenp).
    assert_global_alignment_matches(
        "differing_lengths_widens_band_right",
        b"ACGTACGTACGTACGTACGT",
        b"ACGTACGT",
        DEFAULT_BAND,
    );
}

// ---- regression: row-0 `e` init off-by-one (Critical fix #1) --------------
//
// `AlignAlgo.hpp:266`'s `e[0 + j] = SCORE_GAPOPEN + i * SCORE_GAPOPEN` reads
// the LOOP VARIABLE `i` left over from the immediately preceding `for (i = 1;
// i <= lenp; ++i)` loop -- i.e. `lenp + 1`, not the current `j`. This cell is
// read by the `mat == 0` traceback state (`max = e[tagi * bmax + tagj]`) at
// `tagi == 0`, and a large `lent`-vs-`lenp` skew (band widened on the `lent`
// side) is exactly what steers the traceback into visiting it. `t` = 9×"G",
// `p` = "G", band = 5: verified against the real C++ oracle to score -10 with
// align = [MATCH, DELETE x8]; the pre-fix Rust port returned -9 (a different,
// wrong score) here because the dead-looking `e[0][j]` cell held the wrong
// value.
#[test]
fn regression_row0_e_init_off_by_one_lenp_plus_one() {
    assert_global_alignment_matches(
        "regression_row0_e_init_off_by_one_lenp_plus_one",
        b"GGGGGGGGG",
        b"G",
        5,
    );
}

// ---- short / empty ------------------------------------------------------

#[test]
fn empty_t() {
    assert_global_alignment_matches("empty_t", b"", b"ACGT", DEFAULT_BAND);
}

#[test]
fn empty_p() {
    assert_global_alignment_matches("empty_p", b"ACGT", b"", DEFAULT_BAND);
}

#[test]
fn both_empty() {
    assert_global_alignment_matches("both_empty", b"", b"", DEFAULT_BAND);
}

#[test]
fn single_base_match() {
    assert_global_alignment_matches("single_base_match", b"A", b"A", DEFAULT_BAND);
}

#[test]
fn single_base_mismatch() {
    assert_global_alignment_matches("single_base_mismatch", b"A", b"C", DEFAULT_BAND);
}

#[test]
fn two_base_sequences() {
    assert_global_alignment_matches("two_base_sequences", b"AC", b"AG", DEFAULT_BAND);
}

// ---- seeded pseudo-random pairs at several band values ---------------------

/// Minimal, deterministic xorshift PRNG (no external `rand` dependency
/// needed for this small test matrix) so the "seeded random" cases are
/// reproducible across runs and platforms.
struct XorShift64(u64);
impl XorShift64 {
    fn new(seed: u64) -> Self {
        Self(seed | 1)
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn next_base(&mut self) -> u8 {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        BASES[(self.next_u64() % 4) as usize]
    }
    #[allow(clippy::cast_possible_truncation)]
    fn next_range(&mut self, lo: usize, hi: usize) -> usize {
        lo + (self.next_u64() as usize) % (hi - lo + 1)
    }
}

fn random_seq(rng: &mut XorShift64, len: usize) -> Vec<u8> {
    (0..len).map(|_| rng.next_base()).collect()
}

/// Generates a `p` sequence derived from `t` by applying `n_mismatches`
/// substitutions and `n_indels` single-base insertions/deletions at random
/// positions, so the two sequences stay "close" (within the default band)
/// while still exercising a genuine mix of edit operations.
fn mutate(rng: &mut XorShift64, t: &[u8], n_mismatches: usize, n_indels: usize) -> Vec<u8> {
    let mut p: Vec<u8> = t.to_vec();
    for _ in 0..n_mismatches {
        if p.is_empty() {
            break;
        }
        let pos = rng.next_range(0, p.len() - 1);
        p[pos] = rng.next_base();
    }
    for _ in 0..n_indels {
        if rng.next_u64() % 2 == 0 {
            // insertion
            let pos = rng.next_range(0, p.len());
            p.insert(pos, rng.next_base());
        } else if !p.is_empty() {
            // deletion
            let pos = rng.next_range(0, p.len() - 1);
            p.remove(pos);
        }
    }
    p
}

#[test]
fn seeded_random_pairs_band_default() {
    let mut rng = XorShift64::new(0x00C0_FFEE);
    for trial in 0..30 {
        let t_len = rng.next_range(10, 60);
        let t = random_seq(&mut rng, t_len);
        let n_mismatches = rng.next_range(0, 5);
        let n_indels = rng.next_range(0, 3);
        let p = mutate(&mut rng, &t, n_mismatches, n_indels);
        assert_global_alignment_matches(&format!("seeded_random band=5 trial={trial}"), &t, &p, 5);
    }
}

#[test]
fn seeded_random_pairs_band_small() {
    let mut rng = XorShift64::new(0xDEAD_BEEF);
    for trial in 0..20 {
        let t_len = rng.next_range(10, 40);
        let t = random_seq(&mut rng, t_len);
        let n_mismatches = rng.next_range(0, 3);
        let n_indels = rng.next_range(0, 2);
        let p = mutate(&mut rng, &t, n_mismatches, n_indels);
        assert_global_alignment_matches(&format!("seeded_random band=2 trial={trial}"), &t, &p, 2);
    }
}

#[test]
fn seeded_random_pairs_band_large() {
    let mut rng = XorShift64::new(0x1234_5678_9ABC_DEF0);
    for trial in 0..20 {
        let t_len = rng.next_range(20, 80);
        let t = random_seq(&mut rng, t_len);
        let n_mismatches = rng.next_range(0, 6);
        let n_indels = rng.next_range(0, 5);
        let p = mutate(&mut rng, &t, n_mismatches, n_indels);
        assert_global_alignment_matches(
            &format!("seeded_random band=10 trial={trial}"),
            &t,
            &p,
            10,
        );
    }
}

#[test]
fn seeded_random_pairs_no_mutation() {
    // A batch of pure-random exact-match pairs (t == p), a degenerate but
    // real case that still exercises the "lent == lenp" fast paths.
    let mut rng = XorShift64::new(0x5EED_5EED);
    for trial in 0..10 {
        let t_len = rng.next_range(5, 30);
        let t = random_seq(&mut rng, t_len);
        assert_global_alignment_matches(
            &format!("seeded_random no-mutation trial={trial}"),
            &t,
            &t,
            5,
        );
    }
}

// ---- seeded pseudo-random pairs with LARGE length skew ---------------------
//
// Every seeded-random case above keeps `|lent - lenp|` small (mismatch/indel
// counts of at most a handful against similarly-sized sequences), which
// never exercises the wide-band, off-into-the-boundary-cells code paths that
// both Critical fixes live in. These generators instead pick `lent` and
// `lenp` independently (from very different ranges, including the `== 1`
// extremes) and fill each with unrelated random sequence -- no `mutate`
// relationship -- so the band is forced to widen asymmetrically by a large
// margin, exactly the condition needed to reach the row-0 `e` boundary cell
// (fix #1) and to stress the `tagj`/`tagi` traceback cursors far past their
// "normal" close-diagonal path (fix #2's `usize` underflow class of bug,
// even though the concrete panic reproducer lives in
// `global_alignment_pos_weight`, not here).
fn skewed_pair(
    rng: &mut XorShift64,
    short_len: usize,
    long_lo: usize,
    long_hi: usize,
) -> (Vec<u8>, Vec<u8>) {
    let long_len = rng.next_range(long_lo, long_hi);
    let short = random_seq(rng, short_len);
    let long = random_seq(rng, long_len);
    (long, short)
}

#[test]
fn seeded_random_pairs_length_skew_lent_much_greater_than_lenp() {
    let mut rng = XorShift64::new(0xA11C_E5CD);
    for &band in &[1, 2, 5, 10] {
        for trial in 0..8 {
            let lenp = rng.next_range(1, 3);
            let (t, p) = skewed_pair(&mut rng, lenp, 40, 90);
            assert_global_alignment_matches(
                &format!("length_skew lent>>lenp band={band} trial={trial}"),
                &t,
                &p,
                band,
            );
        }
    }
}

#[test]
fn seeded_random_pairs_length_skew_lenp_much_greater_than_lent() {
    let mut rng = XorShift64::new(0xB0BA_FE77);
    for &band in &[1, 2, 5, 10] {
        for trial in 0..8 {
            let lent = rng.next_range(1, 3);
            let (p, t) = skewed_pair(&mut rng, lent, 40, 90);
            assert_global_alignment_matches(
                &format!("length_skew lenp>>lent band={band} trial={trial}"),
                &t,
                &p,
                band,
            );
        }
    }
}

#[test]
fn tiny_p_against_long_t_all_bands() {
    // lenp == 1 against a long t, at every band width in scope.
    let mut rng = XorShift64::new(0x7EA5_1DE0);
    for &band in &[1, 2, 5, 10] {
        let t = random_seq(&mut rng, 70);
        let p = random_seq(&mut rng, 1);
        assert_global_alignment_matches(
            &format!("tiny_p_against_long_t band={band}"),
            &t,
            &p,
            band,
        );
    }
}

#[test]
fn tiny_t_against_long_p_all_bands() {
    // lent == 1 against a long p, at every band width in scope.
    let mut rng = XorShift64::new(0xC0DE_BABE);
    for &band in &[1, 2, 5, 10] {
        let p = random_seq(&mut rng, 70);
        let t = random_seq(&mut rng, 1);
        assert_global_alignment_matches(
            &format!("tiny_t_against_long_p band={band}"),
            &t,
            &p,
            band,
        );
    }
}

// ---- GlobalAlignment_PosWeight -----------------------------------------

/// Builds a `PosWeight` vector that unanimously "votes" for the given
/// reference sequence (each position has all its weight on the corresponding
/// base -- the common case for a consensus sequence with no observed
/// disagreement).
fn unanimous_pos_weights(seq: &[u8]) -> Vec<PosWeight> {
    seq.iter()
        .map(|&c| {
            let mut count = [0i32; 4];
            let idx = match c {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("unanimous_pos_weights: unsupported base {c:?}"),
            };
            count[idx] = 10;
            PosWeight { count }
        })
        .collect()
}

/// Flattens a `&[PosWeight]` into the `4 * lent`-element `i32` layout the
/// FFI shim expects (see `fg_t1k_alignalgo_global_alignment_pos_weight`'s
/// doc comment in `shim.h`).
fn flatten_weights(weights: &[PosWeight]) -> Vec<i32> {
    weights.iter().flat_map(|w| w.count).collect()
}

fn assert_pos_weight_matches(label: &str, weights: &[PosWeight], p: &[u8]) {
    let rust = global_alignment_pos_weight(weights, p);
    let flat = flatten_weights(weights);
    let cpp = CppAlignAlgo::global_alignment_pos_weight(&flat, p);

    assert_eq!(
        rust.score,
        cpp.score,
        "{label}: score mismatch (p={:?})",
        String::from_utf8_lossy(p)
    );
    assert_eq!(
        rust.align,
        cpp.align,
        "{label}: align[] mismatch (p={:?})",
        String::from_utf8_lossy(p)
    );
}

#[test]
fn pos_weight_exact_match() {
    let t = b"ACGTACGTACGTACGT";
    assert_pos_weight_matches("pos_weight_exact_match", &unanimous_pos_weights(t), t);
}

#[test]
fn pos_weight_with_mismatch() {
    let t = b"ACGTACGTACGTACGT";
    let p = b"ACGTAAGTACGTACGT";
    assert_pos_weight_matches("pos_weight_with_mismatch", &unanimous_pos_weights(t), p);
}

#[test]
fn pos_weight_with_indel() {
    let t = b"ACGTACGTACGTACGTACGT";
    let p = b"ACGTACGTAACGTACGTACGT";
    assert_pos_weight_matches("pos_weight_with_indel", &unanimous_pos_weights(t), p);
}

#[test]
fn pos_weight_ambiguous_column() {
    // A column with a genuine 50/50 split -- IsBaseEqual's minority-vote
    // leniency means BOTH bases at this position should count as "equal".
    let mut weights = unanimous_pos_weights(b"ACGTACGT");
    weights[2] = PosWeight { count: [0, 0, 5, 5] }; // position 2: G/T tie
    assert_pos_weight_matches("pos_weight_ambiguous_column vs G", &weights, b"ACGTACGT");
    assert_pos_weight_matches("pos_weight_ambiguous_column vs T", &weights, b"ACTTACGT");
}

#[test]
fn pos_weight_empty_column() {
    // sum == 0 at one position -> IsBaseEqual always true there regardless
    // of the observed base.
    let mut weights = unanimous_pos_weights(b"ACGTACGT");
    weights[4] = PosWeight { count: [0, 0, 0, 0] };
    assert_pos_weight_matches("pos_weight_empty_column vs A", &weights, b"ACGTACGT");
    assert_pos_weight_matches("pos_weight_empty_column vs T", &weights, b"ACGTTCGT");
}

// ---- regression: traceback `usize` underflow panic (Critical fix #2) ------
//
// C++'s `AlignAlgo::GlobalAlignment_PosWeight` traceback cursors (`int tagi,
// tagj`, `AlignAlgo.hpp:159`) are SIGNED: when `tagj == 0 && tagi > 0` and
// none of the DELETE/INSERT/diagonal branches fire (all three gated on
// `tagj > 0`/`tagi > 0`), the sentinel default `a == 0` (`EDIT_MATCH`)
// survives and the final `else` arm decrements BOTH cursors, sending `tagj`
// to `-1` -- harmless in C++, but a `usize` cursor here panicked with
// "attempt to subtract with overflow" on this exact input before the fix.
// `lent=14, lenp=3, p="CGC"` with the posweight matrix below (verified
// against the real C++ oracle): score=-42, align=[MATCH, DELETE x10, MATCH,
// MATCH, MATCH].
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
    let weights: Vec<PosWeight> = counts.iter().map(|&count| PosWeight { count }).collect();
    let p = b"CGC";

    let rust = global_alignment_pos_weight(&weights, p); // must not panic
    assert_eq!(rust.score, -42, "regression_traceback_usize_underflow_panic: score");
    assert_eq!(
        rust.align,
        vec![
            EDIT_MATCH,
            EDIT_DELETE,
            EDIT_DELETE,
            EDIT_DELETE,
            EDIT_DELETE,
            EDIT_DELETE,
            EDIT_DELETE,
            EDIT_DELETE,
            EDIT_DELETE,
            EDIT_DELETE,
            EDIT_DELETE,
            EDIT_MATCH,
            EDIT_MATCH,
            EDIT_MATCH,
        ],
        "regression_traceback_usize_underflow_panic: align[]"
    );

    assert_pos_weight_matches(
        "regression_traceback_usize_underflow_panic_vs_cpp_oracle",
        &weights,
        p,
    );
}

// ---- pos_weight: LARGE length skew (incl. empty/N/tie columns) ------------

#[test]
fn pos_weight_length_skew_lent_much_greater_than_lenp() {
    let mut rng = XorShift64::new(0x5CA1_AB1E);
    for trial in 0..8 {
        let t_len = rng.next_range(30, 70);
        let t = random_seq(&mut rng, t_len);
        let mut weights = unanimous_pos_weights(&t);
        // Sprinkle an empty (all-zero) column and a tied column into the
        // weights, so length-skew coverage overlaps with the "always equal"
        // / "minority leniency" special cases in `is_base_equal`.
        if weights.len() > 4 {
            weights[1] = PosWeight { count: [0, 0, 0, 0] }; // empty column
            weights[3] = PosWeight { count: [5, 5, 0, 0] }; // A/C tie column
        }
        let lenp = rng.next_range(1, 2);
        let p = random_seq(&mut rng, lenp);
        assert_pos_weight_matches(
            &format!("pos_weight_length_skew lent>>lenp trial={trial}"),
            &weights,
            &p,
        );
    }
}

#[test]
fn pos_weight_length_skew_lenp_much_greater_than_lent() {
    let mut rng = XorShift64::new(0xFACE_FEED);
    for trial in 0..8 {
        let lent = rng.next_range(1, 2);
        let t = random_seq(&mut rng, lent);
        let mut weights = unanimous_pos_weights(&t);
        if lent > 1 {
            weights[0] = PosWeight { count: [0, 0, 0, 0] }; // empty column
        }
        let p_len = rng.next_range(30, 70);
        let p = random_seq(&mut rng, p_len);
        assert_pos_weight_matches(
            &format!("pos_weight_length_skew lenp>>lent trial={trial}"),
            &weights,
            &p,
        );
    }
}

#[test]
fn pos_weight_n_wildcard_with_length_skew() {
    // 'N' in p is a wildcard at every position (IsBaseEqual short-circuits
    // on `c == 'N'`); combine with a big lenp>>lent skew.
    let t = b"ACG";
    let weights = unanimous_pos_weights(t);
    let p = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
    assert_pos_weight_matches("pos_weight_n_wildcard_with_length_skew", &weights, p);
}

#[test]
fn pos_weight_tiny_t_against_long_p() {
    let mut rng = XorShift64::new(0x0BAD_F00D);
    let t = random_seq(&mut rng, 1);
    let weights = unanimous_pos_weights(&t);
    let p = random_seq(&mut rng, 60);
    assert_pos_weight_matches("pos_weight_tiny_t_against_long_p", &weights, &p);
}

#[test]
fn pos_weight_seeded_random() {
    let mut rng = XorShift64::new(0xFEED_FACE);
    for trial in 0..15 {
        let t_len = rng.next_range(10, 40);
        let t = random_seq(&mut rng, t_len);
        let weights = unanimous_pos_weights(&t);
        let n_mismatches = rng.next_range(0, 4);
        let n_indels = rng.next_range(0, 2);
        let p = mutate(&mut rng, &t, n_mismatches, n_indels);
        assert_pos_weight_matches(&format!("pos_weight_seeded_random trial={trial}"), &weights, &p);
    }
}

// ---- GetAlignStats -------------------------------------------------------

fn assert_align_stats_matches(label: &str, align: &[i8], update: bool, initial: (i32, i32, i32)) {
    let (mut rm, mut rmm, mut ri) = initial;
    get_align_stats(align, update, &mut rm, &mut rmm, &mut ri);
    let cpp = CppAlignAlgo::get_align_stats(align, update, initial);
    assert_eq!(
        (rm, rmm, ri),
        cpp,
        "{label}: GetAlignStats mismatch (align={align:?}, update={update})"
    );
}

#[test]
fn align_stats_from_real_global_alignment_output() {
    let t = b"ACGTACGTTTACGTACGTACGGGTACGT";
    let p = b"ACGAACGTACGTACGCACGTACGT";
    let result = global_alignment(t, p, DEFAULT_BAND);
    assert_align_stats_matches("align_stats_from_real_output", &result.align, false, (0, 0, 0));
}

#[test]
fn align_stats_update_true_accumulates_on_top_of_initial() {
    let align = [
        fg_t1k_core::align_algo::EDIT_MATCH,
        fg_t1k_core::align_algo::EDIT_MISMATCH,
        fg_t1k_core::align_algo::EDIT_INSERT,
    ];
    assert_align_stats_matches("align_stats_update_true", &align, true, (5, 5, 5));
}

#[test]
fn align_stats_empty_align() {
    assert_align_stats_matches("align_stats_empty", &[], false, (0, 0, 0));
}

// ---- falsifiability check: confirm align[] comparison has teeth -----------

/// A (t, p) pair with a genuine `f == e` tie at some traceback step, found
/// by brute-force search over seeded-random pairs (see git history of this
/// file for the one-off search harness): `t` and `p` are close enough (one
/// substitution, one indel) that the DP's `e`/`f` gap matrices tie exactly
/// at the point the traceback decides DELETE-vs-INSERT-default. This
/// demonstrates that the FULL `align[]` comparison (not just `score`) is
/// necessary: a deliberately-mistracked traceback (see the local
/// `global_alignment_flipped_tie_break` copy below, which flips the
/// mat==0 DELETE-vs-INSERT-default check from `f >= e` to `f > e`) still
/// returns the IDENTICAL score, but a DIFFERENT `align[]`
/// (`[MATCH,INSERT,MISMATCH,MISMATCH,MATCH,MATCH,MATCH,DELETE,MATCH]` vs.
/// `[MATCH,DELETE,MISMATCH,MISMATCH,MATCH,MATCH,MATCH,INSERT,MATCH]`).
#[test]
fn flip_tie_break_and_confirm_align_array_catches_it() {
    let t = b"CACCACAG";
    let p = b"CGTACACG";

    let real = global_alignment(t, p, DEFAULT_BAND);
    let flipped = global_alignment_flipped_tie_break(t, p, DEFAULT_BAND);

    // The flipped-tie-break variant must still find the SAME optimal score
    // (tie-break order never affects the DP max, only which path is
    // reported) ...
    assert_eq!(real.score, flipped.score, "sanity: flipping a tie-break must not change the score");
    // ... but the two variants may emit different align[] sequences on
    // ties. On this trap case they do diverge, proving a score-only
    // differential would NOT have caught a traceback regression.
    assert_ne!(
        real.align, flipped.align,
        "expected the flipped tie-break to diverge on this deliberately tie-heavy case"
    );

    // And the real (unflipped) Rust port must still agree with the C++
    // oracle byte-for-byte on this same case.
    assert_global_alignment_matches("flip_tie_break_sanity_vs_cpp_oracle", t, p, DEFAULT_BAND);
}

/// A copy of [`fg_t1k_core::align_algo::global_alignment`]'s DP + traceback
/// with exactly ONE line changed from the real port: the `mat == 0`
/// DELETE-vs-INSERT-default tie-break condition `f.get(tagi, tagj) >= max`
/// is weakened to `f.get(tagi, tagj) > max` (strict). This is deliberately
/// NOT `pub` in the core crate -- it exists solely in this test file to
/// demonstrate that the differential's `align[]` comparison would catch such
/// a regression, without introducing a second maintained copy of the real
/// algorithm anywhere production code could accidentally call.
// Variable names (`t`/`p`/`lent`/`lenp`/`tagi`/`tagj`/`m`/`e`/`f`) and
// overall length/structure deliberately mirror `global_alignment` (and, in
// turn, AlignAlgo.hpp itself) line-for-line so the ONE flipped comparison is
// easy to spot by diffing against the real port; see
// fg_t1k_core::align_algo's module docs for the same naming rationale.
#[allow(
    clippy::similar_names,
    clippy::many_single_char_names,
    clippy::too_many_lines,
    clippy::comparison_chain,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss
)]
fn global_alignment_flipped_tie_break(
    t: &[u8],
    p: &[u8],
    band: i32,
) -> fg_t1k_core::align_algo::AlignResult {
    use fg_t1k_core::align_algo::{
        AlignResult, EDIT_DELETE, EDIT_INSERT, EDIT_MATCH, EDIT_MISMATCH,
    };

    const SCORE_MATCH: i32 = 2;
    const SCORE_MISMATCH: i32 = -2;
    const SCORE_GAPOPEN: i32 = -4;
    const SCORE_GAPEXTEND: i32 = -1;

    fn chars_match(t: u8, p: u8) -> bool {
        t == p || t == b'N' || p == b'N'
    }

    struct Matrix {
        data: Vec<i32>,
        bmax: usize,
    }
    impl Matrix {
        fn new(lenp: usize, lent: usize) -> Self {
            Self { data: vec![0; (lenp + 1) * (lent + 1)], bmax: lent + 1 }
        }
        fn get(&self, i: usize, j: usize) -> i32 {
            self.data[i * self.bmax + j]
        }
        fn set(&mut self, i: usize, j: usize, v: i32) {
            self.data[i * self.bmax + j] = v;
        }
    }

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
    for i in 1..=lenp {
        let i_i32 = i as i32;
        e.set(i, 0, SCORE_GAPOPEN + i_i32 * SCORE_GAPEXTEND);
        f.set(i, 0, SCORE_GAPOPEN + i_i32 * SCORE_GAPOPEN);
        m.set(i, 0, SCORE_GAPOPEN + i_i32 * SCORE_GAPOPEN);
    }
    for j in 1..=lent {
        let j_i32 = j as i32;
        f.set(0, j, SCORE_GAPOPEN + j_i32 * SCORE_GAPEXTEND);
        e.set(0, j, SCORE_GAPOPEN + lenp_i32 * SCORE_GAPOPEN);
        m.set(0, j, SCORE_GAPOPEN + j_i32 * SCORE_GAPOPEN);
    }
    for i in 1..=lenp {
        let i_i32 = i as i32;
        let start = if i_i32 - left_band < 1 { 1 } else { i_i32 - left_band };
        let end = if i_i32 + right_band > lent_i32 { lent_i32 } else { i_i32 + right_band };
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
            let mut score = e.get(i - 1, j) + SCORE_GAPEXTEND;
            score = score.max(m.get(i - 1, j) + SCORE_GAPOPEN + SCORE_GAPEXTEND);
            e.set(i, j, score);
            let mut score = f.get(i, j - 1) + SCORE_GAPEXTEND;
            score = score.max(m.get(i, j - 1) + SCORE_GAPOPEN + SCORE_GAPEXTEND);
            f.set(i, j, score);
            let mut score = m.get(i - 1, j - 1)
                + if chars_match(t[j - 1], p[i - 1]) { SCORE_MATCH } else { SCORE_MISMATCH };
            score = score.max(e.get(i, j));
            score = score.max(f.get(i, j));
            m.set(i, j, score);
        }
    }
    let ret = m.get(lenp, lent);

    let mut tagi = lenp;
    let mut tagj = lent;
    let mut mat: u8 = 0;
    let mut align: Vec<i8> = Vec::new();
    while tagi > 0 || tagj > 0 {
        match mat {
            0 => {
                let max = e.get(tagi, tagj);
                let mut a = EDIT_INSERT;
                // *** FLIPPED: `>` instead of `>=` (the real port's line). ***
                if f.get(tagi, tagj) > max {
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
                align.push(EDIT_INSERT);
                if tagi > 0 {
                    if m.get(tagi - 1, tagj) + SCORE_GAPOPEN + SCORE_GAPEXTEND == e.get(tagi, tagj)
                    {
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
                    if m.get(tagi, tagj - 1) + SCORE_GAPOPEN + SCORE_GAPEXTEND == f.get(tagi, tagj)
                    {
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
