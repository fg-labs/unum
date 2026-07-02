#![cfg(feature = "t1k-sys")]
//! Differential test for Task 5a's deterministic `Genotyper` slice: asserts
//! [`fg_t1k_core::genotyper::parse_allele_name`] and
//! [`fg_t1k_core::genotyper::Genotyper::read_assignment_weight`] agree with
//! the real, unmodified vendored C++ `Genotyper::ParseAlleleName`/
//! `ReadAssignmentWeight` (via the opaque `fg_t1k_genotyper_*` FFI shim
//! entries) EXACTLY -- string-for-string for the allele-name parser, and
//! bit-for-bit (`f64` exact equality, not an epsilon -- both are pure
//! comparison/division ratio math with no `pow`/`log`) for the assignment
//! weight.
//!
//! # Coverage
//!
//! `ParseAlleleName`: a battery of KIR names (no delimiter,
//! `KIR2DL1*0010101`-style) and HLA names (colon-delimited,
//! `A*01:01:01:01`/`DRB1*03:01:01`-style) at every `fieldsType` (0 and 1),
//! crossed with every `(alleleDigitUnits, alleleDelimiter)` configuration
//! `Genotyper::SetAlleleNameStructure` can produce: the constructor default
//! (`-1`, `'\0'`), an explicit digit-units override (bypasses delimiter
//! auto-detection entirely -- see `genotyper.rs`'s own unit test for why
//! this is a genuinely surprising T1K behavior worth a dedicated
//! differential case), and an explicit delimiter override. Also covers
//! edge cases: a 2-field HLA name (field-truncation loop runs off the end
//! of the string), and a no-`'*'`-at-all name.
//!
//! `ReadAssignmentWeight`: sweeps `similarity` across and around every
//! segment-threshold boundary (both the default `refSeqSimilarity=0.8` and
//! a high value that exercises the `segment < 0.01` floor), crossed with
//! `hasN` true/false.

use fg_t1k_core::genotyper::{FragmentOverlap, Genotyper, parse_allele_name};
use fg_t1k_sys::CppGenotyper;

// --- ParseAlleleName ---

/// `(allele, fieldsType, alleleDigitUnits, alleleDelimiter)` cases run
/// against both the Rust port and the real C++ oracle.
fn parse_allele_name_cases() -> Vec<(&'static str, i32, i32, u8)> {
    let mut cases = Vec::new();

    let kir_names = [
        "KIR2DL1*0010101",
        "KIR2DL1*00201",
        "KIR3DL1*09502",
        "KIR2DS4*0010102",
        "KIRUNKNOWN", // no '*' at all
    ];
    let hla_names = [
        "A*01:01:01:01",
        "A*01:01",
        "DRB1*03:01:01",
        "DQB1*05:01:01:02",
        "B*56:01:01",
        "C*07:01:01:14",
    ];

    for &fields_type in &[0, 1] {
        for &allele in kir_names.iter().chain(hla_names.iter()) {
            // Constructor default: alleleDigitUnits=-1, alleleDelimiter='\0'.
            cases.push((allele, fields_type, -1, b'\0'));
            // Explicit alleleDigitUnits override (bypasses delimiter
            // auto-detection -- see module docs).
            cases.push((allele, fields_type, 2, b'\0'));
            cases.push((allele, fields_type, 4, b'\0'));
        }
    }

    // Explicit alleleDelimiter override, on both KIR- and HLA-shaped names.
    for &allele in &["KIR2DL1*0010101", "A*01:01:01:01"] {
        cases.push((allele, 0, -1, b':'));
        cases.push((allele, 0, -1, b'*'));
    }

    cases
}

#[test]
fn parse_allele_name_matches_cpp_oracle_exactly() {
    let cases = parse_allele_name_cases();
    assert!(cases.len() > 50, "expected a substantial case battery, got {}", cases.len());

    for (allele, fields_type, allele_digit_units, allele_delimiter) in cases {
        let (rust_gene, rust_major) =
            parse_allele_name(allele, fields_type, allele_digit_units, allele_delimiter);
        let (cpp_gene, cpp_major) = CppGenotyper::parse_allele_name(
            allele,
            fields_type,
            allele_digit_units,
            allele_delimiter,
        );

        assert_eq!(
            rust_gene, cpp_gene,
            "gene mismatch for allele={allele:?} fieldsType={fields_type} \
             alleleDigitUnits={allele_digit_units} alleleDelimiter={allele_delimiter:?}"
        );
        assert_eq!(
            rust_major, cpp_major,
            "majorAllele mismatch for allele={allele:?} fieldsType={fields_type} \
             alleleDigitUnits={allele_digit_units} alleleDelimiter={allele_delimiter:?}"
        );
    }
}

// --- ReadAssignmentWeight ---

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
fn read_assignment_weight_matches_cpp_oracle_exactly() {
    // Sweep similarity densely around and between every segment-threshold
    // boundary, for both a "normal" refSeqSimilarity (0.8, the SeqSet
    // constructor default) and a high one (0.99) that exercises the
    // `segment < 0.01` floor branch.
    let mut similarities: Vec<f64> = Vec::new();
    let mut s = 0.70_f64;
    while s <= 1.001 {
        similarities.push(s);
        s += 0.005;
    }
    // A few exact boundary values too (computed the same way the C++ does,
    // so any fp-rounding quirk at a boundary is exercised identically).
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
        for &similarity in &similarities {
            for &has_n in &[false, true] {
                let rust_w =
                    Genotyper::read_assignment_weight(&frag(similarity, has_n), ref_seq_similarity);
                let cpp_w =
                    CppGenotyper::read_assignment_weight(similarity, has_n, ref_seq_similarity);
                assert_eq!(
                    rust_w.to_bits(),
                    cpp_w.to_bits(),
                    "weight mismatch (exact f64 bits) for similarity={similarity} \
                     hasN={has_n} refSeqSimilarity={ref_seq_similarity}: rust={rust_w} cpp={cpp_w}"
                );
            }
        }
    }
}
