//! Golden-file test for `unum_core::genotyper`'s deterministic scoring
//! slice (`parse_allele_name` / `read_assignment_weight`), converted from the
//! retired the retired T1K-oracle FFI differential (`diff_genotyper_model.rs`) FFI differential (see
//! `tests/common/mod.rs`).
//!
//! Both functions are fully deterministic with ENUMERATED inputs, so this is
//! the "hardcoded inputs -> unit test with captured values" conversion:
//! `parse_allele_name`'s `(gene, majorAllele)` strings and
//! `read_assignment_weight`'s exact `f64` bits (`to_bits()`) are frozen into
//! `genotyper_model.txt` (they were string-/bit-identical to the C++
//! `Genotyper` oracle). The long-allele buffer-overflow regression is kept as
//! a Rust-only self-consistency check.

mod common;

use common::Golden;
use unum_core::genotyper::{FragmentOverlap, Genotyper, parse_allele_name};

fn parse_allele_name_cases() -> Vec<(&'static str, i32, i32, u8)> {
    let mut cases = Vec::new();
    let kir_names = ["KIR2DL1*0010101", "KIR2DL1*00201", "KIR3DL1*09502", "KIR2DS4*0010102", "KIRUNKNOWN"];
    let hla_names = ["A*01:01:01:01", "A*01:01", "DRB1*03:01:01", "DQB1*05:01:01:02", "B*56:01:01", "C*07:01:01:14"];
    for &fields_type in &[0, 1] {
        for &allele in kir_names.iter().chain(hla_names.iter()) {
            cases.push((allele, fields_type, -1, b'\0'));
            cases.push((allele, fields_type, 2, b'\0'));
            cases.push((allele, fields_type, 4, b'\0'));
        }
    }
    for &allele in &["KIR2DL1*0010101", "A*01:01:01:01"] {
        cases.push((allele, 0, -1, b':'));
        cases.push((allele, 0, -1, b'*'));
    }
    cases
}

fn frag(similarity: f64, has_n: bool) -> FragmentOverlap {
    FragmentOverlap {
        seq_idx: 0,
        seq_start: 0,
        seq_end: 100,
        match_cnt: 100,
        relaxed_match_cnt: 100,
        similarity,
        has_mate_pair: true,
        o1_from_r2: false,
        qual: 1.0,
        has_n,
    }
}

#[test]
fn genotyper_model_matches_golden() {
    let mut golden = Golden::open("genotyper_model.txt");

    // ParseAlleleName
    let cases = parse_allele_name_cases();
    assert!(cases.len() > 50, "expected a substantial case battery, got {}", cases.len());
    for (allele, fields_type, digit_units, delimiter) in cases {
        let (gene, major) = parse_allele_name(allele, fields_type, digit_units, delimiter);
        let label = format!("parse/{allele}/ft{fields_type}/du{digit_units}/delim{delimiter}");
        golden.record(label, format!("{gene}|{major}"));
    }

    // ReadAssignmentWeight: sweep similarity around every segment boundary.
    let mut similarities: Vec<f64> = Vec::new();
    let mut s = 0.70_f64;
    while s <= 1.001 {
        similarities.push(s);
        s += 0.005;
    }
    for ref_seq_similarity in [0.8_f64, 0.99_f64] {
        let mut segment = (1.0 - ref_seq_similarity) / 4.0;
        if segment < 0.01 {
            segment = 0.01;
        }
        similarities.push(1.0 - 3.0 * segment);
        similarities.push(1.0 - 2.0 * segment);
        similarities.push(1.0 - segment);
    }

    for &ref_seq_similarity in &[0.8_f64, 0.99_f64, 0.5_f64] {
        for (si, &similarity) in similarities.iter().enumerate() {
            for &has_n in &[false, true] {
                let w = Genotyper::read_assignment_weight(&frag(similarity, has_n), ref_seq_similarity);
                let label = format!("weight/rss{ref_seq_similarity}/s{si:03}/n{}", u8::from(has_n));
                golden.record(label, w.to_bits().to_string());
            }
        }
    }

    golden.finish();
}

/// Rust-only regression (kept): a very long allele must not overflow any
/// fixed buffer and must parse consistently.
#[test]
fn parse_allele_name_long_allele_is_handled() {
    let allele = format!("KIR2DL1*{}", "0".repeat(5000));
    assert!(allele.len() > 4096);
    let (gene, major) = parse_allele_name(&allele, 0, -1, b'\0');
    assert_eq!(gene, "KIR2DL1");
    assert!(!major.is_empty());
}
