#![cfg(feature = "t1k-sys")]
//! Lockstep differential test: drives a Rust `KmerCode` and a real C++
//! `KmerCode` (via the opaque-handle FFI shim) through the same sequence of
//! operations, asserting `get_code`/`canonical`/`rc`/`is_valid`/`kmer_length`
//! agree after every single operation (not just at the end).

use fg_t1k_core::kmer::KmerCode;
use fg_t1k_sys::CppKmerCode;

/// Owns one Rust `KmerCode` and one C++ `CppKmerCode` side by side, and
/// asserts they agree after each driven operation.
struct Lockstep {
    rust: KmerCode,
    cpp: CppKmerCode,
    label: String,
}

impl Lockstep {
    fn new(k: usize, label: &str) -> Self {
        Self {
            rust: KmerCode::new(k),
            cpp: CppKmerCode::new(i32::try_from(k).unwrap()),
            label: label.to_string(),
        }
    }

    fn assert_equal(&self, step: &str) {
        assert_eq!(
            self.rust.get_code(),
            self.cpp.get_code(),
            "{}: get_code mismatch at {step}",
            self.label
        );
        assert_eq!(
            self.rust.get_canonical_kmer_code(),
            self.cpp.canonical(),
            "{}: canonical mismatch at {step}",
            self.label
        );
        assert_eq!(
            self.rust.get_reverse_complement_code(),
            self.cpp.rc(),
            "{}: rc mismatch at {step}",
            self.label
        );
        assert_eq!(
            self.rust.is_valid(),
            self.cpp.is_valid(),
            "{}: is_valid mismatch at {step}",
            self.label
        );
        assert_eq!(
            i32::try_from(self.rust.get_kmer_length()).unwrap(),
            self.cpp.kmer_length(),
            "{}: kmer_length mismatch at {step}",
            self.label
        );
    }

    fn append(&mut self, c: u8) {
        self.rust.append(c);
        self.cpp.append(c);
        self.assert_equal(&format!("append({})", c as char));
    }

    fn prepend(&mut self, c: u8) {
        self.rust.prepend(c);
        self.cpp.prepend(c);
        self.assert_equal(&format!("prepend({})", c as char));
    }

    fn shift_right(&mut self, k: usize) {
        self.rust.shift_right(k);
        self.cpp.shift_right(i32::try_from(k).unwrap());
        self.assert_equal(&format!("shift_right({k})"));
    }

    fn set_code(&mut self, v: u64) {
        self.rust.set_code(v);
        self.cpp.set_code(v);
        self.assert_equal(&format!("set_code({v})"));
    }

    fn restart(&mut self) {
        self.rust.restart();
        self.cpp.restart();
        self.assert_equal("restart");
    }
}

/// Appends `seq` base-by-base (including any invalid/`N` bytes), asserting
/// agreement after each individual base.
fn append_sequence(k: usize, seq: &[u8], label: &str) {
    let mut ls = Lockstep::new(k, label);
    ls.assert_equal("init");
    for &c in seq {
        ls.append(c);
    }
}

#[test]
fn append_matches_t1k_exact_length() {
    let k = 9;
    for seq in ["ACGTACGTA", "TTTTAAAAC", "GATTACAGA", "CCCGGGAAA"] {
        append_sequence(k, seq.as_bytes(), seq);
    }
}

#[test]
fn append_matches_t1k_with_invalid_base() {
    // 'N' is the canonical invalid base; exercises invalid_pos bookkeeping
    // through the full window (position enters, ages, and exits).
    append_sequence(9, b"ACGTNACGTACGT", "ACGTNACGTACGT");
}

#[test]
fn append_matches_t1k_with_non_n_invalid_base() {
    // Append's invalidity check is literally `c == 'N'`, not a general
    // nucToNum lookup -- confirm a different invalid byte ('B') does NOT
    // mark the position invalid (both sides must agree on this quirk).
    append_sequence(9, b"ACGTBACGTACGT", "ACGTBACGTACGT");
}

#[test]
fn append_rolling_window_longer_than_k() {
    // len(seq) > k: only the last k bases should contribute to the code,
    // but invalid_pos bookkeeping must also roll correctly across the window.
    append_sequence(5, b"ACGTACGTNACGTACGT", "rolling-window-with-N");
    append_sequence(5, b"AAAAAAAAAAAAAAAAAAAA", "rolling-window-homopolymer");
}

#[test]
fn shift_right_matches_t1k() {
    let mut ls = Lockstep::new(9, "shift_right");
    for &c in b"ACGTACGTA" {
        ls.append(c);
    }
    ls.shift_right(1);
    ls.shift_right(3);
}

#[test]
fn shift_right_with_invalid_position_matches_t1k() {
    let mut ls = Lockstep::new(9, "shift_right_invalid");
    for &c in b"ACGTNACGTA" {
        ls.append(c);
    }
    // invalid_pos is somewhere in-window; shifting must move it correctly
    // (including the `invalid_pos < 0 => -1` clamp).
    ls.shift_right(1);
    ls.shift_right(1);
    ls.shift_right(5);
}

#[test]
fn set_code_and_restart_match_t1k() {
    let mut ls = Lockstep::new(9, "set_code_restart");
    for &c in b"ACGTNACGTA" {
        ls.append(c);
    }
    // SetCode must reset invalid_pos to -1 (i.e. become valid) even though
    // an 'N' was previously appended.
    ls.set_code(0x1234_5678);
    ls.append(b'A');
    ls.restart();
    ls.append(b'C');
}

#[test]
fn prepend_matches_t1k() {
    let mut ls = Lockstep::new(9, "prepend");
    for &c in b"ACGTACGTA" {
        ls.append(c);
    }
    ls.prepend(b'T');
    ls.prepend(b'G');
}

#[test]
fn prepend_with_invalid_base_matches_t1k() {
    let mut ls = Lockstep::new(9, "prepend_invalid");
    for &c in b"ACGTACGTA" {
        ls.append(c);
    }
    // Prepend's invalidity check uses nucToNum (any non-ACGT byte), unlike
    // Append's literal `c == 'N'` check -- confirm both sides agree that a
    // non-'N' invalid byte still marks the position invalid via Prepend.
    ls.prepend(b'B');
    ls.prepend(b'A');
}

#[test]
fn boundary_k_equals_1() {
    let mut ls = Lockstep::new(1, "k=1");
    for &c in b"ACGTNACGT" {
        ls.append(c);
    }
    ls.shift_right(1);
    ls.prepend(b'T');
}

#[test]
fn boundary_k_equals_15() {
    let mut ls = Lockstep::new(15, "k=15");
    for &c in b"ACGTACGTNACGTACGTACGTACGT" {
        ls.append(c);
    }
    ls.shift_right(2);
    ls.prepend(b'G');
    ls.set_code(0xAAAA);
}

#[test]
fn boundary_k_equals_31() {
    let mut ls = Lockstep::new(31, "k=31");
    for &c in b"ACGTACGTACGTACGTACGTACGTACGTACGTNACGTACGTACGTACGT" {
        ls.append(c);
    }
    ls.shift_right(4);
    ls.prepend(b'A');
    ls.restart();
    for &c in b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA" {
        ls.append(c);
    }
}

#[test]
fn is_equal_matches_code_equality() {
    // is_equal is Rust-only (no direct C++ counterpart needed in the
    // lockstep loop above since it just compares `code` fields), but
    // validate it here against the Rust struct directly.
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
