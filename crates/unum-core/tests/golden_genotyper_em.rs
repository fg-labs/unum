//! Golden-file test for `unum_core::genotyper`'s quantification slice
//! (`coalesce_read_assignments` / `finalize_read_assignments` /
//! `quantify_allele_equivalent_class`), converted from the retired
//! T1K-oracle FFI differential (`diff_genotyper_em.rs`) (see
//! `tests/common/mod.rs`).
//!
//! The scenarios are all ENUMERATED (hardcoded scripted read-assignment
//! lists), so this is the "hardcoded inputs -> unit test with captured
//! values" conversion the brief calls for: the deterministic structure
//! (read-group count, equivalent-class membership/order, per-allele
//! `equivalentClass` and `missingCoverage`) and the SQUAREM-EM abundances are
//! all frozen into `genotyper_em.txt` as the Rust port's own output.
//!
//! # Why the abundances are frozen exactly here (unlike the retired FFI diff)
//!
//! The FFI differential compared abundances to the C++ oracle with a `1e-6`
//! relative tolerance because the oracle's `-O3` FMA-contracted arithmetic
//! rounds differently from Rust's. That tolerance is irrelevant now that the
//! oracle is gone: the golden IS the Rust port's own deterministic `f64`
//! output, so it is frozen bit-for-bit (`f64::to_bits()`). The non-vacuity
//! assertions from the original (allele orderings, nonzero abundances,
//! perturbation-changes-output) are kept as explicit checks alongside the
//! golden.
//!
//! `missing_coverage` is computed Rust-natively via
//! `get_seq_missing_base_coverage` on a fresh (no-pileup) `AlleleRef` per
//! allele -- exactly what the real `genotype` driver does -- rather than being
//! read back from the retired C++ side.

mod common;

use common::Golden;
use unum_core::genotyper::{AlleleRef, FragmentOverlap, Genotyper, get_seq_missing_base_coverage};

struct AlleleDef {
    name: &'static str,
    seq: &'static str,
    weight: i32,
}

struct ScriptedRead {
    hits: Vec<(i32, f64)>,
}

fn read(hits: &[(i32, f64)]) -> ScriptedRead {
    ScriptedRead { hits: hits.to_vec() }
}

/// Runs the full coalesce -> finalize -> quantify pipeline over the scripted
/// scenario, returning the populated genotyper.
fn run(alleles: &[AlleleDef], reads: &[ScriptedRead]) -> Genotyper {
    const REF_SEQ_SIMILARITY: f64 = 0.8;

    let seq_names: Vec<String> = alleles.iter().map(|a| a.name.to_string()).collect();
    let seq_consensus: Vec<Vec<u8>> = alleles.iter().map(|a| a.seq.as_bytes().to_vec()).collect();
    let seq_weight: Vec<i32> = alleles.iter().map(|a| a.weight).collect();
    let mut seq_effective_len: Vec<i32> =
        alleles.iter().map(|a| i32::try_from(a.seq.len()).unwrap()).collect();

    let mut g = Genotyper::new();
    g.init_allele_info(&seq_names, &seq_consensus, &seq_weight, &mut seq_effective_len, 8);
    g.init_read_assignments(i32::try_from(reads.len()).unwrap(), 2000);

    for (read_id, r) in reads.iter().enumerate() {
        let assignment: Vec<FragmentOverlap> = r
            .hits
            .iter()
            .map(|&(allele_idx, qual)| FragmentOverlap {
                seq_idx: allele_idx,
                seq_start: 0,
                seq_end: 10,
                match_cnt: 10,
                relaxed_match_cnt: 10,
                similarity: 1.0,
                has_mate_pair: true,
                o1_from_r2: false,
                qual,
                has_n: false,
            })
            .collect();
        g.set_read_assignments(read_id, &assignment, REF_SEQ_SIMILARITY, |_, _, _| false);
    }
    g.coalesce_read_assignments(0, i32::try_from(reads.len()).unwrap() - 1);

    // Rust-native missing_coverage (no pileup was added, so this is a fresh
    // AlleleRef's own value -- what the real genotype driver computes).
    let missing_coverage: Vec<i32> = alleles
        .iter()
        .map(|a| {
            get_seq_missing_base_coverage(&AlleleRef::new(a.seq.as_bytes().to_vec(), None), 0.01)
        })
        .collect();
    g.finalize_read_assignments(&missing_coverage);
    g.quantify_allele_equivalent_class(&seq_effective_len, &seq_weight);
    g
}

/// Records the scenario's deterministic structure + abundances under
/// `label`-prefixed keys.
fn record_scenario(golden: &mut Golden, label: &str, g: &Genotyper, n_alleles: usize) {
    golden.record(format!("{label}/read_cnt"), g.read_cnt.to_string());
    golden.record(format!("{label}/ec_cnt"), g.equivalent_class_to_alleles.len().to_string());
    for (ec_idx, members) in g.equivalent_class_to_alleles.iter().enumerate() {
        let m = members.iter().map(i32::to_string).collect::<Vec<_>>().join(",");
        golden.record(format!("{label}/ec/{ec_idx}"), m);
    }
    for allele_idx in 0..n_alleles {
        let ai = &g.allele_info[allele_idx];
        golden.record(
            format!("{label}/allele/{allele_idx}"),
            format!(
                "ec={},mc={},abund_bits={},ec_abund_bits={}",
                ai.equivalent_class,
                ai.missing_coverage,
                ai.abundance.to_bits(),
                ai.ec_abundance.to_bits()
            ),
        );
    }
}

fn two_hla_a_alleles() -> Vec<AlleleDef> {
    vec![
        AlleleDef {
            name: "A*01:01:01",
            seq: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            weight: 1,
        },
        AlleleDef {
            name: "A*01:02:01",
            seq: "TTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            weight: 1,
        },
    ]
}

fn two_gene_four_allele_mix() -> Vec<AlleleDef> {
    vec![
        AlleleDef {
            name: "A*01:01:01",
            seq: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            weight: 1,
        },
        AlleleDef {
            name: "A*01:02:01",
            seq: "TTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            weight: 1,
        },
        AlleleDef {
            name: "B*07:02:01",
            seq: "CCCCACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            weight: 1,
        },
        AlleleDef {
            name: "B*08:01:01",
            seq: "GGGGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            weight: 1,
        },
    ]
}

fn three_hla_a_alleles() -> Vec<AlleleDef> {
    vec![
        AlleleDef {
            name: "A*01:01:01",
            seq: "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            weight: 1,
        },
        AlleleDef {
            name: "A*02:01:01",
            seq: "TTTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            weight: 1,
        },
        AlleleDef {
            name: "A*03:01:01",
            seq: "CCCCACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            weight: 1,
        },
    ]
}

#[test]
fn genotyper_em_matches_golden() {
    let mut golden = Golden::open("genotyper_em.txt");

    // clean het (disjoint reads)
    {
        let alleles = two_hla_a_alleles();
        let reads: Vec<ScriptedRead> =
            (0..5).map(|_| read(&[(0, 1.0)])).chain((0..5).map(|_| read(&[(1, 1.0)]))).collect();
        let g = run(&alleles, &reads);
        record_scenario(&mut golden, "clean_het", &g, alleles.len());
        assert!(g.allele_info[0].abundance > 0.0);
        assert!(g.allele_info[1].abundance > 0.0);
    }

    // homozygous
    {
        let alleles = two_hla_a_alleles();
        let reads: Vec<ScriptedRead> = (0..8).map(|_| read(&[(0, 1.0)])).collect();
        let g = run(&alleles, &reads);
        record_scenario(&mut golden, "homozygous", &g, alleles.len());
        assert!(g.allele_info[0].abundance > 0.0);
        #[allow(clippy::float_cmp)]
        let a1_zero = g.allele_info[1].abundance == 0.0;
        assert!(a1_zero, "allele 1 got zero reads -> zero abundance");
    }

    // shared-read ambiguity
    {
        let alleles = two_hla_a_alleles();
        let mut reads: Vec<ScriptedRead> = (0..10).map(|_| read(&[(0, 1.0)])).collect();
        reads.extend((0..2).map(|_| read(&[(1, 1.0)])));
        reads.extend((0..4).map(|_| read(&[(0, 1.0), (1, 1.0)])));
        let g = run(&alleles, &reads);
        record_scenario(&mut golden, "shared_read", &g, alleles.len());
        assert!(g.allele_info[0].abundance > g.allele_info[1].abundance);
    }

    // equivalent alleles collapse to one EC
    {
        let alleles = two_hla_a_alleles();
        let reads: Vec<ScriptedRead> = (0..6).map(|_| read(&[(0, 1.0), (1, 1.0)])).collect();
        let g = run(&alleles, &reads);
        assert_eq!(g.equivalent_class_to_alleles.len(), 1);
        record_scenario(&mut golden, "equivalent_alleles", &g, alleles.len());
        assert!((g.allele_info[0].abundance - g.allele_info[1].abundance).abs() < 1e-9);
        assert!(g.allele_info[0].abundance > 0.0);
    }

    // multi-gene mix
    {
        let alleles = two_gene_four_allele_mix();
        let mut reads: Vec<ScriptedRead> = (0..5).map(|_| read(&[(0, 1.0)])).collect();
        reads.extend((0..5).map(|_| read(&[(1, 1.0)])));
        reads.extend((0..7).map(|_| read(&[(2, 1.0)])));
        reads.extend((0..3).map(|_| read(&[(3, 1.0)])));
        let g = run(&alleles, &reads);
        record_scenario(&mut golden, "multi_gene", &g, alleles.len());
        assert!(g.allele_info[2].abundance > g.allele_info[3].abundance);
        assert!(g.allele_info.iter().all(|a| a.abundance > 0.0));
    }

    // three-allele multi-EC apportionment
    {
        let alleles = three_hla_a_alleles();
        let quals = [0.5, 1.0, 2.0];
        let mut reads: Vec<ScriptedRead> = Vec::new();
        for allele_idx in 0..3 {
            for i in 0..15 {
                reads.push(read(&[(allele_idx, quals[i % quals.len()])]));
            }
        }
        for (a, b) in [(0, 1), (1, 2), (0, 2)] {
            for i in 0..8 {
                reads.push(read(&[(a, quals[i % quals.len()]), (b, quals[(i + 1) % quals.len()])]));
            }
        }
        for i in 0..6 {
            reads.push(read(&[
                (0, quals[i % quals.len()]),
                (1, quals[(i + 1) % quals.len()]),
                (2, quals[(i + 2) % quals.len()]),
            ]));
        }
        let g = run(&alleles, &reads);
        record_scenario(&mut golden, "three_allele_multi_ec", &g, alleles.len());
        assert!(g.allele_info.iter().all(|a| a.abundance > 0.0));
    }

    golden.finish();
}

/// Confirms the scenarios are NON-VACUOUS: perturbing one scripted read's
/// target allele changes the abundances.
#[test]
fn perturbing_a_scripted_read_changes_abundances() {
    let alleles = two_hla_a_alleles();
    let baseline: Vec<ScriptedRead> =
        (0..5).map(|_| read(&[(0, 1.0)])).chain((0..5).map(|_| read(&[(1, 1.0)]))).collect();
    let baseline_g = run(&alleles, &baseline);

    let mut perturbed: Vec<ScriptedRead> = (0..8).map(|_| read(&[(0, 1.0)])).collect();
    perturbed.extend((0..2).map(|_| read(&[(1, 1.0)])));
    let perturbed_g = run(&alleles, &perturbed);

    assert!(
        (baseline_g.allele_info[0].abundance - perturbed_g.allele_info[0].abundance).abs() > 1e-6,
        "perturbing the read distribution must change abundance"
    );
}
