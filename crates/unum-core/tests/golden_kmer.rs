//! Golden-file tests for `unum_core::kmer`, converted from the retired
//! the retired T1K-oracle FFI differential (`diff_kmer.rs`) and `diff_kmercode.rs` FFI differentials (see
//! `tests/common/mod.rs`). The former lockstep-drove the Rust `KmerCode` and a
//! real C++ `KmerCode` and asserted agreement after every operation; here the
//! same operation sequences are driven and each operation's full observable
//! state (`get_code`/`canonical`/`rc`/`is_valid`/`kmer_length`) is frozen into
//! the `kmer.txt` golden.
//!
//! Dropped from the original `diff_kmercode.rs`: the two `#[should_panic]`
//! FFI-guard tests (`cpp_shift_right_rejects_full_width`/`_negative_count`) --
//! they guarded the C++ shim's argument validation, which no longer exists.
//! `is_equal_matches_code_equality` (already Rust-only) is kept.

mod common;

use common::Golden;
use unum_core::kmer::{KmerCode, canonical_kmer};

// ---- canonical_kmer (hardcoded expected values) ---------------------------

/// `canonical_kmer` cases with their frozen expected canonical codes (the
/// Rust port's output, which was byte-identical to the T1K oracle
/// when the C++ shim still existed). Small enough to assert as inline
/// literals per the conversion brief.
#[test]
fn canonical_kmer_matches_frozen_values() {
    let cases: &[(&str, usize, u64)] = &[
        ("ACGTACGTA", 9, CANON_ACGTACGTA_9),
        ("TTTTAAAAC", 9, CANON_TTTTAAAAC_9),
        ("GATTACAGA", 9, CANON_GATTACAGA_9),
        ("CCCGGGAAA", 9, CANON_CCCGGGAAA_9),
        ("ACGTNACGT", 9, CANON_ACGTNACGT_9),
        ("ACGTACGTACGT", 9, CANON_ACGTACGTACGT_9),
        ("ACG", 9, CANON_ACG_9),
        ("ACGTACGTA", 1, CANON_ACGTACGTA_1),
        ("ACGTACGTACGTACGTACGTACGTACGTACGT", 31, CANON_LONG_31),
    ];
    for &(seq, k, expected) in cases {
        assert_eq!(canonical_kmer(seq.as_bytes(), k), expected, "canonical_kmer({seq:?}, {k})");
    }
}

// Frozen canonical k-mer codes (the Rust port's output, byte-identical to the
// retired T1K-oracle C++ shim). These constants are asserted, so a port change
// to `canonical_kmer` fails the test here rather than silently drifting.
const CANON_ACGTACGTA_9: u64 = 27756;
const CANON_TTTTAAAAC_9: u64 = 196_352;
const CANON_GATTACAGA_9: u64 = 146_504;
const CANON_CCCGGGAAA_9: u64 = 88704;
const CANON_ACGTNACGT_9: u64 = 27675;
const CANON_ACGTACGTACGT_9: u64 = 27756;
const CANON_LONG_31: u64 = 488_296_166_657_017_542;
const CANON_ACG_9: u64 = 6;
const CANON_ACGTACGTA_1: u64 = 0;

// ---- KmerCode operation sequences (frozen golden) -------------------------

/// Records the full observable state of `kc` under `label`:
/// `code,canonical,rc,is_valid,kmer_length`.
fn record_state(golden: &mut Golden, label: &str, kc: &KmerCode) {
    golden.record(
        label,
        format!(
            "{},{},{},{},{}",
            kc.get_code(),
            kc.get_canonical_kmer_code(),
            kc.get_reverse_complement_code(),
            u8::from(kc.is_valid()),
            kc.get_kmer_length()
        ),
    );
}

/// Drives one op sequence, recording state after `init` and after each op.
struct Driver<'g> {
    golden: &'g mut Golden,
    kc: KmerCode,
    prefix: String,
    step: usize,
}

impl<'g> Driver<'g> {
    fn new(golden: &'g mut Golden, k: usize, prefix: &str) -> Self {
        let kc = KmerCode::new(k);
        let mut d = Self { golden, kc, prefix: prefix.to_string(), step: 0 };
        d.snapshot("init");
        d
    }
    fn snapshot(&mut self, op: &str) {
        let label = format!("{}/{:03}/{op}", self.prefix, self.step);
        self.step += 1;
        record_state(self.golden, &label, &self.kc);
    }
    fn append(&mut self, c: u8) {
        self.kc.append(c);
        self.snapshot(&format!("append({})", c as char));
    }
    fn prepend(&mut self, c: u8) {
        self.kc.prepend(c);
        self.snapshot(&format!("prepend({})", c as char));
    }
    fn shift_right(&mut self, k: usize) {
        self.kc.shift_right(k);
        self.snapshot(&format!("shift_right({k})"));
    }
    fn set_code(&mut self, v: u64) {
        self.kc.set_code(v);
        self.snapshot(&format!("set_code({v})"));
    }
    fn restart(&mut self) {
        self.kc.restart();
        self.snapshot("restart");
    }
}

fn append_sequence(golden: &mut Golden, k: usize, seq: &[u8], label: &str) {
    let mut d = Driver::new(golden, k, label);
    for &c in seq {
        d.append(c);
    }
}

#[test]
fn kmer_code_op_sequences_match_golden() {
    let mut golden = Golden::open("kmer.txt");

    for seq in ["ACGTACGTA", "TTTTAAAAC", "GATTACAGA", "CCCGGGAAA"] {
        append_sequence(&mut golden, 9, seq.as_bytes(), seq);
    }
    append_sequence(&mut golden, 9, b"ACGTNACGTACGT", "ACGTNACGTACGT");
    append_sequence(&mut golden, 9, b"ACGTBACGTACGT", "ACGTBACGTACGT");
    append_sequence(&mut golden, 5, b"ACGTACGTNACGTACGT", "rolling-window-with-N");
    append_sequence(&mut golden, 5, b"AAAAAAAAAAAAAAAAAAAA", "rolling-window-homopolymer");

    {
        let mut d = Driver::new(&mut golden, 9, "shift_right");
        for &c in b"ACGTACGTA" {
            d.append(c);
        }
        d.shift_right(1);
        d.shift_right(3);
    }
    {
        let mut d = Driver::new(&mut golden, 9, "shift_right_invalid");
        for &c in b"ACGTNACGTA" {
            d.append(c);
        }
        d.shift_right(1);
        d.shift_right(1);
        d.shift_right(5);
    }
    {
        let mut d = Driver::new(&mut golden, 9, "set_code_restart");
        for &c in b"ACGTNACGTA" {
            d.append(c);
        }
        d.set_code(0x1234_5678);
        d.append(b'A');
        d.restart();
        d.append(b'C');
    }
    {
        let mut d = Driver::new(&mut golden, 9, "prepend");
        for &c in b"ACGTACGTA" {
            d.append(c);
        }
        d.prepend(b'T');
        d.prepend(b'G');
    }
    {
        let mut d = Driver::new(&mut golden, 9, "prepend_invalid");
        for &c in b"ACGTACGTA" {
            d.append(c);
        }
        d.prepend(b'B');
        d.prepend(b'A');
    }
    {
        let mut d = Driver::new(&mut golden, 9, "prepend_partial_window");
        for &c in b"ACG" {
            d.append(c);
        }
        d.prepend(b'T');
        d.prepend(b'A');
        d.prepend(b'N');
    }
    {
        let mut d = Driver::new(&mut golden, 1, "k=1");
        for &c in b"ACGTNACGT" {
            d.append(c);
        }
        d.shift_right(1);
        d.prepend(b'T');
    }
    {
        let mut d = Driver::new(&mut golden, 15, "k=15");
        for &c in b"ACGTACGTNACGTACGTACGTACGT" {
            d.append(c);
        }
        d.shift_right(2);
        d.prepend(b'G');
        d.set_code(0xAAAA);
    }
    {
        let mut d = Driver::new(&mut golden, 31, "k=31");
        for &c in b"ACGTACGTACGTACGTACGTACGTACGTACGTNACGTACGTACGTACGT" {
            d.append(c);
        }
        d.shift_right(4);
        d.prepend(b'A');
        d.restart();
        for &c in b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA" {
            d.append(c);
        }
    }
    {
        let mut d = Driver::new(&mut golden, 4, "k=4");
        for &c in b"ACGTNACGT" {
            d.append(c);
        }
        d.shift_right(1);
        d.prepend(b'T');
        d.set_code(0b1010_1010);
    }
    {
        let mut d = Driver::new(&mut golden, 8, "k=8");
        for &c in b"ACGTACGTNACGTACGT" {
            d.append(c);
        }
        d.shift_right(2);
        d.prepend(b'G');
        d.restart();
        for &c in b"TTTTGGGGCCCCAAAA" {
            d.append(c);
        }
    }
    {
        let mut d = Driver::new(&mut golden, 16, "k=16");
        for &c in b"ACGTACGTACGTACGTNACGTACGTACGTACGT" {
            d.append(c);
        }
        d.shift_right(3);
        d.prepend(b'A');
        d.set_code(0xDEAD_BEEF);
    }
    {
        let mut d = Driver::new(&mut golden, 30, "k=30");
        for &c in b"ACGTACGTACGTACGTACGTACGTACGTACNGTACGTACGTACGTACGTACGTACGTACGT" {
            d.append(c);
        }
        d.shift_right(4);
        d.prepend(b'A');
        d.restart();
        for &c in b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATT" {
            d.append(c);
        }
    }
    {
        let mut d = Driver::new(&mut golden, 32, "k=32");
        for &c in b"ACGTACGTACGTACGTACGTACGTACGTACGTNACGTACGTACGTACGTACGTACGTACGTACGT" {
            d.append(c);
        }
        d.shift_right(4);
        d.prepend(b'A');
        d.set_code(u64::MAX);
        d.restart();
        for &c in b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA" {
            d.append(c);
        }
    }

    golden.finish();
}

/// `is_equal` is Rust-only (no C++ counterpart) -- kept verbatim from the
/// original differential.
#[test]
fn is_equal_matches_code_equality() {
    let mut a = KmerCode::new(5);
    let mut b = KmerCode::new(5);
    for &c in b"ACGTA" {
        a.append(c);
    }
    for &c in b"ACGTA" {
        b.append(c);
    }
    assert!(a.is_equal(&b));
    b.append(b'C');
    assert!(!a.is_equal(&b));
}
