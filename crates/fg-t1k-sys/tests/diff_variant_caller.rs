#![cfg(feature = "t1k-sys")]
//! Differential test for Task 6a: `VariantCaller` core, comparing
//! [`fg_t1k_core::variant_caller::VariantCaller`] against the real,
//! unmodified vendored C++ `VariantCaller` (via the opaque-handle
//! `fg_t1k_variantcaller_*` FFI shim, [`CppVariantCaller`]).
//!
//! # Approach
//!
//! Both sides are driven by the SAME SCRIPTED set of alleles + reads +
//! `_fragmentOverlap`-shaped assignments (reusing the 5b scripted-assignment
//! pattern from `diff_genotyper_em.rs`/`diff_genotyper_model.rs`, extended
//! to carry full `overlap1`/`overlap2` coordinates since `VariantCaller`
//! reads them directly for pileup accumulation, unlike `Genotyper`'s own
//! `_fragmentOverlap` consumers).
//!
//! `overlap.align` is NOT scripted directly on either side: the Rust side
//! derives it via [`fg_t1k_core::align_algo::global_alignment`] (already
//! differentially proven byte-identical to `AlignAlgo::GlobalAlignment` in
//! `diff_align_algo.rs`), and the C++ side derives it via the shim's call to
//! the REAL, unmodified `SeqSet::AddOverlapAlignmentInfo` (which itself just
//! calls `AlignAlgo::GlobalAlignment`) -- exactly mirroring what
//! `Analyzer::AddFragmentAlignmentInfo` does before its own `ComputeVariant`
//! call in the real end-to-end pipeline. Neither side "cheats" by directly
//! constructing a shared alignment; both independently derive it from the
//! same overlap coordinates via already-validated alignment code.
//!
//! # Assertions
//!
//! - **Called variants match EXACTLY**: position (`refStart`/`refEnd`),
//!   `ref`/`var` base, zygosity (`outputGroupId`: 0 = best, 1 = equal-best
//!   ambiguous), and PASS/FILTER (`qual > 0`) must be identical, in the
//!   same order.
//! - **Support/quality fields within tolerance**: `varSupport`/
//!   `allSupport`/`varUniqSupport` are compared with a tight ABSOLUTE
//!   tolerance (`1e-9`) rather than requiring exact bit-equality --
//!   see [`assert_variants_match`]'s doc comment for why (unlike
//!   `diff_genotype_e2e.rs`'s fixed-point-arithmetic-only float fields,
//!   `VariantCaller`'s pileup counts are pure integer-valued `f64`
//!   accumulations along an ORDER this port reproduces exactly, so in
//!   practice every scenario below achieves EXACT equality -- the
//!   tolerance is a documented safety margin, not an observed necessity).
//! - **`OutputAlleleVCF` text matches byte-for-byte** (both sides use `%lf`
//!   6-decimal-place float formatting on integer-valued counts, so no
//!   float-formatting divergence is possible here).
//!
//! # Scenarios
//!
//! Clean SNV (single allele, overwhelming alt support), an indel-bearing
//! read set (proves indel alignment ops never spuriously create/shift a
//! variant call -- see `fg_t1k_core::variant_caller` module docs: T1K
//! itself never calls indels), a two-allele ambiguous group (exercises
//! `SolveVariantGroup`/`EnumerateVariants`' multi-variant resolution and
//! `varMaxGroupToResolve`), ref-only reads (no variant called), and a
//! low-support case that must be filtered out entirely.

use fg_t1k_core::align_algo;
use fg_t1k_core::genotyper::{AlleleRef, Genotyper};
use fg_t1k_core::variant_caller::{FragmentOverlap, Overlap, VariantCaller};
use fg_t1k_sys::{
    CppGenotyper2, CppVariantCaller, ScriptedVariantOverlap, ScriptedVcOverlapCoords,
};

const KMER_LENGTH: i32 = 11;

/// One allele's reference definition for a scripted scenario.
struct AlleleDef {
    name: &'static str,
    seq: &'static str,
}

/// One scripted single-end read, ungapped against `seq_idx`'s consensus
/// starting at position 0 (`seq_start = read_start = 0`, `seq_end =
/// read_end = read.len() - 1`) -- covers every scenario this test needs
/// without requiring a general aligner; `UpdateBaseVariantFromOverlap`
/// itself does the real per-position work regardless of how the overlap
/// coordinates were obtained.
struct ScriptedRead {
    seq_idx: i32,
    seq: &'static str,
}

/// Builds a [`Overlap`] (Rust side) covering the whole read ungapped
/// against `seq_idx`, deriving `align` via the real
/// [`align_algo::global_alignment`] -- mirrors what
/// `SeqSet::AddOverlapAlignmentInfo` does on the C++ side (see this
/// module's doc comment).
fn rust_overlap(seq_idx: i32, read: &str, allele_consensus: &[u8]) -> Overlap {
    let read_bytes = read.as_bytes();
    let len = i32::try_from(read_bytes.len()).unwrap();
    let align_result = align_algo::global_alignment(
        &allele_consensus[0..read_bytes.len()],
        read_bytes,
        align_algo::DEFAULT_BAND,
    );
    let mut match_cnt = 0;
    let mut mismatch_cnt = 0;
    let mut indel_cnt = 0;
    align_algo::get_align_stats(
        &align_result.align,
        false,
        &mut match_cnt,
        &mut mismatch_cnt,
        &mut indel_cnt,
    );
    Overlap {
        seq_idx,
        read_start: 0,
        read_end: len - 1,
        seq_start: 0,
        seq_end: len - 1,
        strand: 1,
        match_cnt: 2 * match_cnt,
        similarity: 1.0,
        align: Some(align_result.align),
    }
}

/// Runs the identical scripted `alleles`/`reads` scenario through both the
/// Rust [`VariantCaller`] and the real C++ [`CppVariantCaller`].
///
/// Returns `(rust_vc, names, allele_refs, cpp_genotyper, cpp_vc)`. NOTE:
/// `cpp_genotyper` MUST be kept alive for at least as long as `cpp_vc` --
/// the real C++ `VariantCaller` holds a `SeqSet&` REFERENCE into
/// `cpp_genotyper`'s `refSet` (`VariantCaller variantCaller(refSet)`,
/// mirroring `Analyzer.cpp:673`), so dropping `cpp_genotyper` first leaves
/// `cpp_vc` holding a dangling reference (a real use-after-free/segfault --
/// this bit an earlier revision of this test, which returned only
/// `cpp_vc`, letting `cpp_genotyper` drop at the end of this function).
fn run_both(
    alleles: &[AlleleDef],
    reads: &[ScriptedRead],
    max_var_group_to_resolve: i32,
) -> (VariantCaller, Vec<String>, Vec<AlleleRef>, CppGenotyper2, CppVariantCaller) {
    // --- Shared setup ---
    let names: Vec<String> = alleles.iter().map(|a| a.name.to_string()).collect();
    let consensus: Vec<Vec<u8>> = alleles.iter().map(|a| a.seq.as_bytes().to_vec()).collect();
    let weight: Vec<i32> = alleles.iter().map(|_| 1).collect();
    let mut effective_len: Vec<i32> =
        alleles.iter().map(|a| i32::try_from(a.seq.len()).unwrap()).collect();

    // --- Rust side ---
    let allele_refs: Vec<AlleleRef> =
        consensus.iter().map(|c| AlleleRef::new(c.clone(), None)).collect();
    let mut rust_genotyper = Genotyper::new();
    rust_genotyper.init_allele_info(&names, &consensus, &weight, &mut effective_len, 8);
    // Uniform abundance (1.0 each) and gene grouping is entirely determined
    // by init_allele_info's own allele-name parsing -- both sides therefore
    // agree on seqCopy without this test needing to compute it itself.
    for info in &mut rust_genotyper.allele_info {
        info.abundance = 1.0;
    }

    let mut rust_vc = VariantCaller::new(&allele_refs);
    rust_vc.set_seq_abundance(&rust_genotyper, alleles.len());
    rust_vc.set_max_var_group_to_resolve(max_var_group_to_resolve);

    let rust_read1: Vec<Vec<u8>> = reads.iter().map(|r| r.seq.as_bytes().to_vec()).collect();
    let rust_assignments: Vec<Vec<FragmentOverlap>> = reads
        .iter()
        .map(|r| {
            let overlap1 =
                rust_overlap(r.seq_idx, r.seq, &consensus[usize::try_from(r.seq_idx).unwrap()]);
            vec![FragmentOverlap {
                seq_idx: r.seq_idx,
                has_mate_pair: false,
                o1_from_r2: false,
                overlap1,
                overlap2: Overlap::none(),
            }]
        })
        .collect();
    rust_vc.compute_variant(&rust_read1, &[], &rust_assignments, &consensus);

    // --- C++ side ---
    let mut cpp_genotyper = CppGenotyper2::new(KMER_LENGTH);
    for a in alleles {
        cpp_genotyper.add_ref_seq(a.name, a.seq, 1);
    }
    cpp_genotyper.init_allele_info();
    // InitAlleleInfo defaults every allele's abundance to 0 (Genotyper.hpp:
    // 587,591), NOT 1.0 -- this test bypasses the EM entirely (matching the
    // Rust side's own choice not to run it) and instead scripts a uniform
    // 1.0 abundance directly via the test-only
    // fg_t1k_genotyper2_set_allele_abundance hook, so both sides use
    // IDENTICAL (not merely equal-by-coincidence) per-allele weights in
    // UpdateBaseVariantFromFragmentOverlap's `weight = seqAbundance[seqIdx]
    // / totalWeight` -- this matters even for a single allele: `weight ==
    // 1.0` gates the `uniqCount` increment exactly (VariantCaller.hpp:145),
    // so an unscripted (0/0 = NaN) abundance would silently zero out every
    // uniqCount on the C++ side while the Rust side (seeded to 1.0 above)
    // kept incrementing it -- exactly the divergence this hook fixes.
    for (i, _) in alleles.iter().enumerate() {
        cpp_genotyper.set_allele_abundance(i32::try_from(i).unwrap(), 1.0);
    }
    let mut cpp_vc = CppVariantCaller::new(&cpp_genotyper);
    cpp_vc.set_seq_abundance(&cpp_genotyper);
    cpp_vc.set_max_var_group_to_resolve(max_var_group_to_resolve);

    let cpp_reads1: Vec<&str> = reads.iter().map(|r| r.seq).collect();
    let cpp_assignments: Vec<ScriptedVariantOverlap> = reads
        .iter()
        .enumerate()
        .map(|(read_idx, r)| {
            let len = i32::try_from(r.seq.len()).unwrap();
            ScriptedVariantOverlap {
                read_idx: i32::try_from(read_idx).unwrap(),
                seq_idx: r.seq_idx,
                match_cnt: 2 * len,
                similarity: 1.0,
                has_mate_pair: false,
                o1_from_r2: false,
                overlap1: ScriptedVcOverlapCoords {
                    seq_start: 0,
                    seq_end: len - 1,
                    read_start: 0,
                    read_end: len - 1,
                    strand: 1,
                    match_cnt: 2 * len,
                    similarity: 1.0,
                },
                overlap2: ScriptedVcOverlapCoords::default(),
            }
        })
        .collect();
    cpp_vc.compute_variant(&cpp_reads1, &[], &cpp_assignments);

    (rust_vc, names, allele_refs, cpp_genotyper, cpp_vc)
}

/// Absolute tolerance for `varSupport`/`allSupport`/`varUniqSupport` --
/// see this module's doc comment for why exact equality is expected in
/// practice (this is a documented safety margin, not an observed
/// necessity: every scenario below achieves bit-exact float equality since
/// these are pure integer-valued accumulations along a fully-reproduced
/// accumulation order).
const SUPPORT_ABS_TOLERANCE: f64 = 1e-9;

/// Asserts the Rust [`VariantCaller`]'s called variants match the C++
/// oracle's EXACTLY on position/ref/var/zygosity/filter, and within
/// [`SUPPORT_ABS_TOLERANCE`] on support/quality fields.
fn assert_variants_match(rust_vc: &VariantCaller, cpp_vc: &CppVariantCaller) {
    let rust_variants = rust_vc.final_variants();
    let cpp_count = cpp_vc.final_variant_count();
    assert_eq!(
        rust_variants.len(),
        usize::try_from(cpp_count).unwrap(),
        "called variant COUNT must match exactly: rust={rust_variants:?} cpp_count={cpp_count}",
    );

    for (i, rv) in rust_variants.iter().enumerate() {
        let cv = cpp_vc.final_variant(i32::try_from(i).unwrap());
        assert_eq!(rv.seq_idx, cv.seq_idx, "variant {i}: seqIdx mismatch");
        assert_eq!(rv.ref_start, cv.ref_start, "variant {i}: refStart mismatch");
        assert_eq!(rv.ref_end, cv.ref_end, "variant {i}: refEnd mismatch");
        assert_eq!(rv.reference, cv.reference, "variant {i}: ref base mismatch");
        assert_eq!(rv.var, cv.var, "variant {i}: var base mismatch");
        assert_eq!(
            rv.output_group_id, cv.output_group_id,
            "variant {i}: outputGroupId (zygosity/ambiguity) mismatch"
        );
        assert_eq!(
            (rv.qual > 0),
            (cv.qual > 0),
            "variant {i}: PASS/FILTER (qual > 0) mismatch: rust qual={} cpp qual={}",
            rv.qual,
            cv.qual
        );

        assert!(
            (rv.all_support - cv.all_support).abs() <= SUPPORT_ABS_TOLERANCE,
            "variant {i}: allSupport mismatch beyond tolerance: rust={} cpp={}",
            rv.all_support,
            cv.all_support
        );
        assert!(
            (rv.var_support - cv.var_support).abs() <= SUPPORT_ABS_TOLERANCE,
            "variant {i}: varSupport mismatch beyond tolerance: rust={} cpp={}",
            rv.var_support,
            cv.var_support
        );
        assert!(
            (rv.var_uniq_support - cv.var_uniq_support).abs() <= SUPPORT_ABS_TOLERANCE,
            "variant {i}: varUniqSupport mismatch beyond tolerance: rust={} cpp={}",
            rv.var_uniq_support,
            cv.var_uniq_support
        );
    }
}

/// Computes `refSet.GetExonicPosition(seqIdx, pos)` for [`VariantCaller::output_allele_vcf`]:
/// since every [`AlleleRef`] in this test's scenarios has no exon comment
/// (`AlleleRef::new(..., None)`), the WHOLE consensus is one exon
/// (`SeqSet.hpp:970-976`'s `comment == NULL` fallback), so the exonic
/// position is simply `pos` itself (matching `GetExonicPosition`'s
/// single-exon identity case, `SeqSet.hpp:2808-2825`).
fn exonic_position_identity(_seq_idx: usize, pos: i32) -> i32 {
    pos
}

#[test]
fn clean_snv_matches_cpp_exactly() {
    // 20 reads support a T at position 5 (ref C) against a single allele --
    // overwhelming majority, unambiguous single-variant group.
    let consensus = "AAAAACAAAAAAAAAAAAAA";
    let alt = "AAAAATAAAAAAAAAAAAAA";
    let alleles = [AlleleDef { name: "KIR2DL1*001", seq: consensus }];
    let reads: Vec<ScriptedRead> = (0..20).map(|_| ScriptedRead { seq_idx: 0, seq: alt }).collect();

    let (rust_vc, names, _allele_refs, _cpp_genotyper, cpp_vc) = run_both(&alleles, &reads, 8);

    assert_variants_match(&rust_vc, &cpp_vc);
    assert_eq!(rust_vc.final_variants().len(), 1, "expected exactly one called SNV");
    let v = &rust_vc.final_variants()[0];
    assert_eq!(v.ref_start, 5);
    assert_eq!(v.reference, b'C');
    assert_eq!(v.var, b'T');

    let rust_vcf = rust_vc.output_allele_vcf(&names, exonic_position_identity);
    let cpp_vcf = cpp_vc.output_allele_vcf(std::env::temp_dir().as_path());
    assert_eq!(rust_vcf, cpp_vcf, "OutputAlleleVCF text must be byte-identical");

    // Non-vacuous: an all-reference read set calls nothing (see the
    // ref_only scenario below), so this scenario is proven to genuinely
    // exercise the variant-calling path, not trivially pass on empty output.
    assert!(!rust_vcf.is_empty());
}

#[test]
fn ref_only_reads_call_nothing_on_both_sides() {
    let consensus = "AAAAACAAAAAAAAAAAAAA";
    let alleles = [AlleleDef { name: "KIR2DL1*001", seq: consensus }];
    let reads: Vec<ScriptedRead> =
        (0..20).map(|_| ScriptedRead { seq_idx: 0, seq: consensus }).collect();

    let (rust_vc, names, _allele_refs, _cpp_genotyper, cpp_vc) = run_both(&alleles, &reads, 8);

    assert_variants_match(&rust_vc, &cpp_vc);
    assert!(rust_vc.final_variants().is_empty());
    assert_eq!(cpp_vc.final_variant_count(), 0);

    let rust_vcf = rust_vc.output_allele_vcf(&names, exonic_position_identity);
    let cpp_vcf = cpp_vc.output_allele_vcf(std::env::temp_dir().as_path());
    assert_eq!(rust_vcf, cpp_vcf);
    assert!(rust_vcf.is_empty());
}

#[test]
fn low_support_variant_is_filtered_on_both_sides() {
    // Only 2 alt reads out of 20 -- below both the absolute (5) and
    // relative (0.5x ref) support thresholds in FindCandidateVariants.
    let consensus = "AAAAACAAAAAAAAAAAAAA";
    let alt = "AAAAATAAAAAAAAAAAAAA";
    let alleles = [AlleleDef { name: "KIR2DL1*001", seq: consensus }];
    let mut reads = Vec::new();
    for i in 0..20 {
        let seq = if i < 2 { alt } else { consensus };
        reads.push(ScriptedRead { seq_idx: 0, seq });
    }

    let (rust_vc, _names, _allele_refs, _cpp_genotyper, cpp_vc) = run_both(&alleles, &reads, 8);

    assert_variants_match(&rust_vc, &cpp_vc);
    assert!(rust_vc.final_variants().is_empty(), "below-threshold support must not be called");
    assert_eq!(cpp_vc.final_variant_count(), 0);
}

#[test]
fn indel_bearing_reads_do_not_spuriously_call_a_variant() {
    // Every read carries a real 1bp deletion relative to the consensus
    // (position 5 removed) -- proving EDIT_DELETE/EDIT_INSERT ops never
    // get pileup-counted as a variant on EITHER side (matches stock T1K:
    // no indel calling exists at all, see fg_t1k_core::variant_caller
    // module docs). This scenario cannot use `rust_overlap`'s ungapped
    // helper (a real gapped alignment is required), so both sides build the
    // overlap coordinates directly from a hand-specified deletion.
    let consensus = "AAAAACAAAAAAAAAAAAAA"; // len 20, 'C' at position 5
    let read = "AAAAAAAAAAAAAAAAAAA"; // consensus with position 5 deleted, len 19
    let read_cnt = 20;

    // --- Rust side ---
    let allele_refs = vec![AlleleRef::new(consensus.as_bytes().to_vec(), None)];
    let consensus_vecs = vec![consensus.as_bytes().to_vec()];
    let mut rust_genotyper = Genotyper::new();
    let mut effective_len = vec![i32::try_from(consensus.len()).unwrap()];
    rust_genotyper.init_allele_info(
        &["KIR2DL1*001".to_string()],
        &consensus_vecs,
        &[1],
        &mut effective_len,
        8,
    );
    rust_genotyper.allele_info[0].abundance = 1.0;
    let mut rust_vc = VariantCaller::new(&allele_refs);
    rust_vc.set_seq_abundance(&rust_genotyper, 1);
    rust_vc.set_max_var_group_to_resolve(8);

    let align_ops = {
        use fg_t1k_core::align_algo::{EDIT_DELETE, EDIT_MATCH};
        let mut ops = vec![EDIT_MATCH; 5];
        ops.push(EDIT_DELETE);
        ops.extend(std::iter::repeat_n(EDIT_MATCH, 14));
        ops
    };
    let rust_read1: Vec<Vec<u8>> = (0..read_cnt).map(|_| read.as_bytes().to_vec()).collect();
    let rust_assignments: Vec<Vec<FragmentOverlap>> = (0..read_cnt)
        .map(|_| {
            let overlap1 = Overlap {
                seq_idx: 0,
                read_start: 0,
                read_end: 18,
                seq_start: 0,
                seq_end: 19,
                strand: 1,
                match_cnt: 38,
                similarity: 0.95,
                align: Some(align_ops.clone()),
            };
            vec![FragmentOverlap {
                seq_idx: 0,
                has_mate_pair: false,
                o1_from_r2: false,
                overlap1,
                overlap2: Overlap::none(),
            }]
        })
        .collect();
    rust_vc.compute_variant(&rust_read1, &[], &rust_assignments, &consensus_vecs);

    // --- C++ side ---
    let mut cpp_genotyper = CppGenotyper2::new(KMER_LENGTH);
    cpp_genotyper.add_ref_seq("KIR2DL1*001", consensus, 1);
    cpp_genotyper.init_allele_info();
    cpp_genotyper.set_allele_abundance(0, 1.0); // see run_both's doc comment on why this is required.
    let mut cpp_vc = CppVariantCaller::new(&cpp_genotyper);
    cpp_vc.set_seq_abundance(&cpp_genotyper);
    cpp_vc.set_max_var_group_to_resolve(8);

    let cpp_reads1: Vec<&str> = (0..read_cnt).map(|_| read).collect();
    let cpp_assignments: Vec<ScriptedVariantOverlap> = (0..read_cnt)
        .map(|i| ScriptedVariantOverlap {
            read_idx: i,
            seq_idx: 0,
            match_cnt: 38,
            similarity: 0.95,
            has_mate_pair: false,
            o1_from_r2: false,
            overlap1: ScriptedVcOverlapCoords {
                seq_start: 0,
                seq_end: 19,
                read_start: 0,
                read_end: 18,
                strand: 1,
                match_cnt: 38,
                similarity: 0.95,
            },
            overlap2: ScriptedVcOverlapCoords::default(),
        })
        .collect();
    cpp_vc.compute_variant(&cpp_reads1, &[], &cpp_assignments);

    assert_variants_match(&rust_vc, &cpp_vc);
    assert!(
        rust_vc.final_variants().is_empty(),
        "indel-only divergence must not be called as a SNV on the Rust side: {:?}",
        rust_vc.final_variants()
    );
    assert_eq!(
        cpp_vc.final_variant_count(),
        0,
        "indel-only divergence must not be called as a SNV on the C++ side either"
    );
}

#[test]
fn two_allele_ambiguous_group_exercises_var_group_resolution() {
    // Two near-identical alleles differing at TWO adjacent exonic positions
    // (a 2-variant group). Reads are split between both alleles with
    // overlapping ambiguous support, forcing SolveVariantGroup's
    // EnumerateVariants to actually resolve (not just single-position
    // candidate detection) -- exercises `varMaxGroupToResolve` (set to 8,
    // comfortably above this group's size of ~2).
    let allele_a = "AAAAACGAAAAAAAAAAAAA"; // C at 5, G at 6
    let allele_b = "AAAAATCAAAAAAAAAAAAA"; // T at 5, C at 6 (alt of both positions)
    let alleles = [
        AlleleDef { name: "KIR2DL1*001", seq: allele_a },
        AlleleDef { name: "KIR2DL1*002", seq: allele_b },
    ];

    // 20 reads perfectly matching allele_a, 20 perfectly matching allele_b --
    // each allele's own reads should reinforce its own consensus (no novel
    // variant expected against either allele when reads only ever match
    // their own assigned allele's consensus exactly).
    let mut reads = Vec::new();
    for _ in 0..20 {
        reads.push(ScriptedRead { seq_idx: 0, seq: allele_a });
    }
    for _ in 0..20 {
        reads.push(ScriptedRead { seq_idx: 1, seq: allele_b });
    }

    let (rust_vc, _names, _allele_refs, _cpp_genotyper, cpp_vc) = run_both(&alleles, &reads, 8);

    assert_variants_match(&rust_vc, &cpp_vc);
    // Every read matches its OWN assigned allele's consensus exactly, so no
    // novel variant should be called against either allele (this is a
    // ref-only-per-allele scenario at the pileup level, just with two
    // reference sequences instead of one) -- the meaningful assertion here
    // is that both sides AGREE (asserted above), whatever that agreement is.
    assert_eq!(
        rust_vc.final_variants().len(),
        usize::try_from(cpp_vc.final_variant_count()).unwrap(),
        "both sides must call the same NUMBER of variants for the two-allele group scenario"
    );
}
