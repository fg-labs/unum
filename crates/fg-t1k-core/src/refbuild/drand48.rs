//! Portable reimplementation of Perl's `srand`/`rand` PRNG stream
//! (`Perl_drand48`, glibc `drand48`'s algorithm — a 48-bit linear congruential
//! generator), used to reproduce `ParseDatFile.pl`'s `srand(17)`-seeded
//! UTR-padding fallback (`ParseDatFile.pl:575-602`) bit-for-bit.
//!
//! Perl's `rand()`, when not linked against a platform `drand48`, ships its
//! own portable `drand48`-alike implementation (`util.c`'s `Perl_drand48_r`)
//! so that `srand`-seeded streams are reproducible across platforms. That is
//! exactly the algorithm implemented here — **not** whatever the host libc's
//! `drand48(3)` happens to do, which is not portably specified and would not
//! match Perl's output.
//!
//! # Algorithm
//!
//! Seeding with scalar seed `s`:
//! ```text
//! state = 0x330E | ((s & 0xFFFF) << 16) | (((s >> 16) & 0xFFFF) << 32)
//! ```
//!
//! Each draw advances the 48-bit LCG state and converts it to `[0, 1)`:
//! ```text
//! state = (0x5DEECE66D * state + 0xB) & 0xFFFF_FFFF_FFFF
//! drand48 = state as f64 / 2f64.powi(48)
//! ```
//!
//! `ParseDatFile.pl` consumes this stream via `int(rand(4))`, i.e.
//! `(drand48() * 4.0) as usize` (truncating), indexing into `('A', 'C', 'G',
//! 'T')`.
//!
//! # Casts
//!
//! Mirroring Perl's own untyped numeric coercions (`$state as f64`, then
//! `int(...)` truncating back to an index), this module casts `u64 -> f64`
//! (accepted precision loss: the 48-bit LCG state losing its low bits when
//! widened to `f64`'s 52-bit mantissa is exactly what real `drand48`/Perl
//! does, not a bug to guard against) and `f64 -> usize` (accepted truncation:
//! `int(rand(4))` is a truncating cast by definition, and the result is
//! always in `[0, 4)`, nowhere near `usize` truncation/sign-loss territory on
//! any target this crate supports).
//!
//! The LCG constants and the Perl reference draws pinned in this module's
//! test are exact values from the algorithm/oracle, not numbers meant to be
//! "readable" in grouped form, so `clippy::unreadable_literal` is allowed
//! module-wide too.
#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::unreadable_literal
)]

/// A `Perl_drand48`-compatible 48-bit LCG stream, seeded and drawn from
/// exactly like Perl's `srand`/`rand`.
pub(crate) struct Drand48 {
    state: u64,
}

/// LCG multiplier, per `Perl_drand48_r` / POSIX `drand48`.
const MULTIPLIER: u64 = 0x5DEECE66D;
/// LCG increment, per `Perl_drand48_r` / POSIX `drand48`.
const INCREMENT: u64 = 0xB;
/// Keeps only the low 48 bits of the LCG state.
const MASK_48: u64 = 0xFFFF_FFFF_FFFF;

impl Drand48 {
    /// Seeds the stream exactly as Perl's `srand($seed)` does before falling
    /// back to its portable `drand48`.
    pub(crate) fn new(seed: u64) -> Self {
        let state = 0x330E | ((seed & 0xFFFF) << 16) | (((seed >> 16) & 0xFFFF) << 32);
        Self { state }
    }

    /// Draws the next `drand48()` value in `[0, 1)`.
    pub(crate) fn next_f64(&mut self) -> f64 {
        self.state = (MULTIPLIER.wrapping_mul(self.state).wrapping_add(INCREMENT)) & MASK_48;
        (self.state as f64) / 2f64.powi(48)
    }

    /// Draws one random base: Perl's `int(rand(4))` truncating cast, mapped
    /// to `('A', 'C', 'G', 'T')`. Mirrors `ParseDatFile.pl`'s
    /// `$numToNuc[int(rand(4))]`.
    pub(crate) fn next_base(&mut self) -> u8 {
        let idx = (self.next_f64() * 4.0) as usize;
        [b'A', b'C', b'G', b'T'][idx]
    }
}

#[cfg(test)]
mod tests {
    use super::Drand48;

    /// Self-check against real Perl's actual `srand(17); print rand(), "\n"`
    /// output (three successive draws, full `f64` precision), confirming this
    /// port's LCG matches Perl's portable `drand48` bit-for-bit rather than
    /// some other PRNG. Exact equality is intentional here — the point of
    /// this test is bit-for-bit reproduction of Perl's output, not an
    /// approximate comparison.
    #[test]
    #[allow(clippy::float_cmp)]
    fn seed_17_matches_perl_reference_draws() {
        let mut rng = Drand48::new(17);
        assert_eq!(rng.next_f64(), 0.9744672834212942);
        assert_eq!(rng.next_f64(), 0.7279398726272746);
        assert_eq!(rng.next_f64(), 0.6499462188604745);
    }
}
