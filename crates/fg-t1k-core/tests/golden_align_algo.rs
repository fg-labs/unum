//! Golden-file test for `fg_t1k_core::align_algo`, converted from the retired
//! `fg-t1k-sys` `diff_align_algo.rs` FFI differential (see
//! `tests/common/mod.rs` for the conversion rationale).
//!
//! # Why compare the FULL `align[]` array, not just the score
//!
//! Because scoring is all-integer, the DP score alone is comparatively easy
//! to get right -- the recurrences are `MAX(...)` chains, and tie order never
//! affects the numeric max. The traceback tie-break order, however, decides
//! WHICH of several equally-scoring alignments is emitted as the op sequence.
//! A score-only golden would happily pass a traceback that silently picks a
//! different (but equally-scored) alignment path -- e.g. reporting a run of 3
//! consecutive single-base mismatches where T1K would report a 1-base
//! insertion + 1-base deletion, changing downstream indel/mismatch counts.
//! Every case therefore freezes `(score, align)`, not just `score`.
//!
//! `flip_tie_break_and_confirm_align_array_catches_it` (at the bottom) is the
//! falsifiability check: it patches a COPY of the Rust traceback logic with
//! one tie-break flipped, confirms score is unaffected but `align[]` diverges
//! on a case chosen to exercise a genuine multi-way score tie -- proving the
//! array-level golden actually has teeth. (Its former C++-oracle arm is
//! dropped; the self-consistent flip check remains.)

mod common;

use common::Golden;
use fg_t1k_core::align_algo::{
    DEFAULT_BAND, EDIT_DELETE, EDIT_MATCH, PosWeight, get_align_stats, global_alignment,
    global_alignment_pos_weight,
};

/// Serializes an `(score, align)` result to the golden's per-case value:
/// `score;a0,a1,a2,...` (align ops as their raw `i8` codes).
fn serialize(score: i32, align: &[i8]) -> String {
    let ops = align.iter().map(|a| a.to_string()).collect::<Vec<_>>().join(",");
    format!("{score};{ops}")
}

/// Runs the Rust port on `(t, p)` at `band` and records `(score, align)`
/// under `label`.
fn record_global_alignment(golden: &mut Golden, label: &str, t: &[u8], p: &[u8], band: i32) {
    let r = global_alignment(t, p, band);
    golden.record(label, serialize(r.score, &r.align));
}

// ---- exact match, mismatches, indels, mixed / band edges / short-empty -----

fn record_curated_global_cases(golden: &mut Golden) {
    record_global_alignment(golden, "exact_match", b"ACGTACGTACGTACGT", b"ACGTACGTACGTACGT", DEFAULT_BAND);
    record_global_alignment(golden, "single_mismatch", b"ACGTACGTACGTACGT", b"ACGTAAGTACGTACGT", DEFAULT_BAND);
    record_global_alignment(
        golden,
        "multiple_mismatches",
        b"ACGTACGTACGTACGTACGT",
        b"ACCTACGAACGTACCTACGT",
        DEFAULT_BAND,
    );
    record_global_alignment(
        golden,
        "single_insertion",
        b"ACGTACGTACGTACGT",
        b"ACGTACCGTACGTACGTACGT",
        DEFAULT_BAND,
    );
    record_global_alignment(
        golden,
        "single_deletion",
        b"ACGTACCGTACGTACGTACGT",
        b"ACGTACGTACGTACGT",
        DEFAULT_BAND,
    );
    record_global_alignment(
        golden,
        "mixed_indels_and_mismatches",
        b"ACGTACGTTTACGTACGTACGGGTACGT",
        b"ACGAACGTACGTACGCACGTACGT",
        DEFAULT_BAND,
    );
    record_global_alignment(golden, "n_wildcard_bases", b"ACGTNACGTACGT", b"ACGTAACGTNCGT", DEFAULT_BAND);
    record_global_alignment(
        golden,
        "indel_at_band_edge_small_band",
        b"ACGTACGTAACGTACGTACGTACGT",
        b"ACGTACGTACGTACGTACGTACGT",
        2,
    );
    record_global_alignment(golden, "indel_at_band_edge_band_one", b"AAAAAAAAAAAAAAAA", b"AAAAAAAAAAAAAAAAA", 1);
    record_global_alignment(
        golden,
        "differing_lengths_widens_band_left",
        b"ACGTACGT",
        b"ACGTACGTACGTACGTACGT",
        DEFAULT_BAND,
    );
    record_global_alignment(
        golden,
        "differing_lengths_widens_band_right",
        b"ACGTACGTACGTACGTACGT",
        b"ACGTACGT",
        DEFAULT_BAND,
    );
    record_global_alignment(golden, "empty_t", b"", b"ACGT", DEFAULT_BAND);
    record_global_alignment(golden, "empty_p", b"ACGT", b"", DEFAULT_BAND);
    record_global_alignment(golden, "both_empty", b"", b"", DEFAULT_BAND);
    record_global_alignment(golden, "single_base_match", b"A", b"A", DEFAULT_BAND);
    record_global_alignment(golden, "single_base_mismatch", b"A", b"C", DEFAULT_BAND);
    record_global_alignment(golden, "two_base_sequences", b"AC", b"AG", DEFAULT_BAND);
}

// ---- regression: row-0 `e` init off-by-one (Critical fix #1) --------------
//
// `AlignAlgo.hpp:266`'s `e[0 + j] = SCORE_GAPOPEN + i * SCORE_GAPOPEN` reads
// the LOOP VARIABLE `i` left over from the preceding loop -- i.e. `lenp + 1`,
// not the current `j`. `t` = 9x"G", `p` = "G", band = 5: scores -10 with
// align = [MATCH, DELETE x8]; the pre-fix Rust port returned -9 here.
fn record_regression_row0(golden: &mut Golden) {
    let r = global_alignment(b"GGGGGGGGG", b"G", 5);
    assert_eq!(r.score, -10, "regression_row0: score");
    let mut expected = vec![EDIT_MATCH];
    expected.extend(std::iter::repeat_n(EDIT_DELETE, 8));
    assert_eq!(r.align, expected, "regression_row0: align");
    golden.record("regression_row0_e_init_off_by_one_lenp_plus_one", serialize(r.score, &r.align));
}

// ---- seeded pseudo-random pairs at several band values ---------------------

/// Minimal, deterministic xorshift PRNG (kept verbatim from the retired FFI
/// differential so the frozen goldens remain reproducible byte-for-byte).
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
            let pos = rng.next_range(0, p.len());
            p.insert(pos, rng.next_base());
        } else if !p.is_empty() {
            let pos = rng.next_range(0, p.len() - 1);
            p.remove(pos);
        }
    }
    p
}

fn record_seeded_random(golden: &mut Golden) {
    let mut rng = XorShift64::new(0x00C0_FFEE);
    for trial in 0..30 {
        let t_len = rng.next_range(10, 60);
        let t = random_seq(&mut rng, t_len);
        let n_mismatches = rng.next_range(0, 5);
        let n_indels = rng.next_range(0, 3);
        let p = mutate(&mut rng, &t, n_mismatches, n_indels);
        record_global_alignment(golden, &format!("seeded_random band=5 trial={trial}"), &t, &p, 5);
    }

    let mut rng = XorShift64::new(0xDEAD_BEEF);
    for trial in 0..20 {
        let t_len = rng.next_range(10, 40);
        let t = random_seq(&mut rng, t_len);
        let n_mismatches = rng.next_range(0, 3);
        let n_indels = rng.next_range(0, 2);
        let p = mutate(&mut rng, &t, n_mismatches, n_indels);
        record_global_alignment(golden, &format!("seeded_random band=2 trial={trial}"), &t, &p, 2);
    }

    let mut rng = XorShift64::new(0x1234_5678_9ABC_DEF0);
    for trial in 0..20 {
        let t_len = rng.next_range(20, 80);
        let t = random_seq(&mut rng, t_len);
        let n_mismatches = rng.next_range(0, 6);
        let n_indels = rng.next_range(0, 5);
        let p = mutate(&mut rng, &t, n_mismatches, n_indels);
        record_global_alignment(golden, &format!("seeded_random band=10 trial={trial}"), &t, &p, 10);
    }

    let mut rng = XorShift64::new(0x5EED_5EED);
    for trial in 0..10 {
        let t_len = rng.next_range(5, 30);
        let t = random_seq(&mut rng, t_len);
        record_global_alignment(golden, &format!("seeded_random no-mutation trial={trial}"), &t, &t, 5);
    }
}

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

fn record_length_skew(golden: &mut Golden) {
    let mut rng = XorShift64::new(0xA11C_E5CD);
    for &band in &[1, 2, 5, 10] {
        for trial in 0..8 {
            let lenp = rng.next_range(1, 3);
            let (t, p) = skewed_pair(&mut rng, lenp, 40, 90);
            record_global_alignment(
                golden,
                &format!("length_skew lent>>lenp band={band} trial={trial}"),
                &t,
                &p,
                band,
            );
        }
    }

    let mut rng = XorShift64::new(0xB0BA_FE77);
    for &band in &[1, 2, 5, 10] {
        for trial in 0..8 {
            let lent = rng.next_range(1, 3);
            let (p, t) = skewed_pair(&mut rng, lent, 40, 90);
            record_global_alignment(
                golden,
                &format!("length_skew lenp>>lent band={band} trial={trial}"),
                &t,
                &p,
                band,
            );
        }
    }

    let mut rng = XorShift64::new(0x7EA5_1DE0);
    for &band in &[1, 2, 5, 10] {
        let t = random_seq(&mut rng, 70);
        let p = random_seq(&mut rng, 1);
        record_global_alignment(golden, &format!("tiny_p_against_long_t band={band}"), &t, &p, band);
    }

    let mut rng = XorShift64::new(0xC0DE_BABE);
    for &band in &[1, 2, 5, 10] {
        let p = random_seq(&mut rng, 70);
        let t = random_seq(&mut rng, 1);
        record_global_alignment(golden, &format!("tiny_t_against_long_p band={band}"), &t, &p, band);
    }
}

// ---- GlobalAlignment_PosWeight -----------------------------------------

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

fn record_pos_weight(golden: &mut Golden, label: &str, weights: &[PosWeight], p: &[u8]) {
    let r = global_alignment_pos_weight(weights, p);
    golden.record(label, serialize(r.score, &r.align));
}

fn record_pos_weight_cases(golden: &mut Golden) {
    let t = b"ACGTACGTACGTACGT";
    record_pos_weight(golden, "pos_weight_exact_match", &unanimous_pos_weights(t), t);
    record_pos_weight(golden, "pos_weight_with_mismatch", &unanimous_pos_weights(t), b"ACGTAAGTACGTACGT");
    record_pos_weight(
        golden,
        "pos_weight_with_indel",
        &unanimous_pos_weights(b"ACGTACGTACGTACGTACGT"),
        b"ACGTACGTAACGTACGTACGT",
    );

    let mut weights = unanimous_pos_weights(b"ACGTACGT");
    weights[2] = PosWeight { count: [0, 0, 5, 5] };
    record_pos_weight(golden, "pos_weight_ambiguous_column vs G", &weights, b"ACGTACGT");
    record_pos_weight(golden, "pos_weight_ambiguous_column vs T", &weights, b"ACTTACGT");

    let mut weights = unanimous_pos_weights(b"ACGTACGT");
    weights[4] = PosWeight { count: [0, 0, 0, 0] };
    record_pos_weight(golden, "pos_weight_empty_column vs A", &weights, b"ACGTACGT");
    record_pos_weight(golden, "pos_weight_empty_column vs T", &weights, b"ACGTTCGT");

    // regression: traceback usize underflow panic (Critical fix #2). Must not
    // panic; score=-42, align=[MATCH, DELETE x10, MATCH x3].
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
    let r = global_alignment_pos_weight(&weights, b"CGC");
    assert_eq!(r.score, -42, "regression_traceback_usize_underflow: score");
    let mut expected = vec![EDIT_MATCH];
    expected.extend(std::iter::repeat_n(EDIT_DELETE, 10));
    expected.extend([EDIT_MATCH, EDIT_MATCH, EDIT_MATCH]);
    assert_eq!(r.align, expected, "regression_traceback_usize_underflow: align");
    golden.record("regression_traceback_usize_underflow_panic", serialize(r.score, &r.align));

    // length-skew posweight
    let mut rng = XorShift64::new(0x5CA1_AB1E);
    for trial in 0..8 {
        let t_len = rng.next_range(30, 70);
        let t = random_seq(&mut rng, t_len);
        let mut weights = unanimous_pos_weights(&t);
        if weights.len() > 4 {
            weights[1] = PosWeight { count: [0, 0, 0, 0] };
            weights[3] = PosWeight { count: [5, 5, 0, 0] };
        }
        let lenp = rng.next_range(1, 2);
        let p = random_seq(&mut rng, lenp);
        record_pos_weight(golden, &format!("pos_weight_length_skew lent>>lenp trial={trial}"), &weights, &p);
    }

    let mut rng = XorShift64::new(0xFACE_FEED);
    for trial in 0..8 {
        let lent = rng.next_range(1, 2);
        let t = random_seq(&mut rng, lent);
        let mut weights = unanimous_pos_weights(&t);
        if lent > 1 {
            weights[0] = PosWeight { count: [0, 0, 0, 0] };
        }
        let p_len = rng.next_range(30, 70);
        let p = random_seq(&mut rng, p_len);
        record_pos_weight(golden, &format!("pos_weight_length_skew lenp>>lent trial={trial}"), &weights, &p);
    }

    let t = b"ACG";
    let weights = unanimous_pos_weights(t);
    let p = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
    record_pos_weight(golden, "pos_weight_n_wildcard_with_length_skew", &weights, p);

    let mut rng = XorShift64::new(0x0BAD_F00D);
    let t = random_seq(&mut rng, 1);
    let weights = unanimous_pos_weights(&t);
    let p = random_seq(&mut rng, 60);
    record_pos_weight(golden, "pos_weight_tiny_t_against_long_p", &weights, &p);

    let mut rng = XorShift64::new(0xFEED_FACE);
    for trial in 0..15 {
        let t_len = rng.next_range(10, 40);
        let t = random_seq(&mut rng, t_len);
        let weights = unanimous_pos_weights(&t);
        let n_mismatches = rng.next_range(0, 4);
        let n_indels = rng.next_range(0, 2);
        let p = mutate(&mut rng, &t, n_mismatches, n_indels);
        record_pos_weight(golden, &format!("pos_weight_seeded_random trial={trial}"), &weights, &p);
    }
}

// ---- GetAlignStats -------------------------------------------------------

fn record_align_stats(
    golden: &mut Golden,
    label: &str,
    align: &[i8],
    update: bool,
    initial: (i32, i32, i32),
) {
    let (mut m, mut mm, mut i) = initial;
    get_align_stats(align, update, &mut m, &mut mm, &mut i);
    golden.record(label, format!("{m},{mm},{i}"));
}

fn record_align_stats_cases(golden: &mut Golden) {
    let t = b"ACGTACGTTTACGTACGTACGGGTACGT";
    let p = b"ACGAACGTACGTACGCACGTACGT";
    let result = global_alignment(t, p, DEFAULT_BAND);
    record_align_stats(golden, "align_stats_from_real_output", &result.align, false, (0, 0, 0));

    let align = [
        fg_t1k_core::align_algo::EDIT_MATCH,
        fg_t1k_core::align_algo::EDIT_MISMATCH,
        fg_t1k_core::align_algo::EDIT_INSERT,
    ];
    record_align_stats(golden, "align_stats_update_true", &align, true, (5, 5, 5));
    record_align_stats(golden, "align_stats_empty", &[], false, (0, 0, 0));
}

#[test]
fn align_algo_matches_golden() {
    let mut golden = Golden::open("align_algo.txt");
    record_curated_global_cases(&mut golden);
    record_regression_row0(&mut golden);
    record_seeded_random(&mut golden);
    record_length_skew(&mut golden);
    record_pos_weight_cases(&mut golden);
    record_align_stats_cases(&mut golden);
    golden.finish();
}

// ---- falsifiability check: confirm align[] comparison has teeth -----------

/// A `(t, p)` pair with a genuine `f == e` tie at some traceback step. The
/// flipped-tie-break variant returns the IDENTICAL score but a DIFFERENT
/// `align[]`, proving the full-array golden (not just score) is necessary.
/// (Its former C++-oracle arm is dropped; the self-consistent flip check
/// remains -- the byte-for-byte match against T1K is now frozen into the
/// `align_algo.txt` golden above.)
#[test]
fn flip_tie_break_and_confirm_align_array_catches_it() {
    let t = b"CACCACAG";
    let p = b"CGTACACG";

    let real = global_alignment(t, p, DEFAULT_BAND);
    let flipped = global_alignment_flipped_tie_break(t, p, DEFAULT_BAND);

    assert_eq!(real.score, flipped.score, "sanity: flipping a tie-break must not change the score");
    assert_ne!(
        real.align, flipped.align,
        "expected the flipped tie-break to diverge on this deliberately tie-heavy case"
    );
}

/// A copy of [`fg_t1k_core::align_algo::global_alignment`]'s DP + traceback
/// with exactly ONE line changed: the `mat == 0` DELETE-vs-INSERT-default
/// tie-break condition `f.get(tagi, tagj) >= max` weakened to `> max`. Not
/// `pub` in the core crate -- it exists solely here to demonstrate the
/// golden's `align[]` comparison would catch such a regression.
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
    use fg_t1k_core::align_algo::{AlignResult, EDIT_DELETE, EDIT_INSERT, EDIT_MATCH, EDIT_MISMATCH};

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
                    if m.get(tagi - 1, tagj) + SCORE_GAPOPEN + SCORE_GAPEXTEND == e.get(tagi, tagj) {
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
                    if m.get(tagi, tagj - 1) + SCORE_GAPOPEN + SCORE_GAPEXTEND == f.get(tagi, tagj) {
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
