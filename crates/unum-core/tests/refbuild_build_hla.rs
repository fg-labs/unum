//! Golden byte-identity test for the full reference-build orchestration
//! (`unum_core::refbuild::build_reference`) against the HLA subset fixture,
//! which (unlike the KIR fixture in `refbuild_build.rs`) exercises
//! `ParseDatFile.pl`'s `srand(17)`-seeded random-UTR-padding fallback for the
//! HLA-DRB2/HLA-DRB7 pseudogenes (RNA mode only) — see
//! `fixtures/refbuild/PINS.md`'s "Phase 1b, Task 1b.1" section and
//! `crates/unum-core/src/refbuild/dat.rs`'s `drand48` module docs.

use unum_core::refbuild;

fn fx(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/refbuild").join(rel)
}

#[test]
fn build_reference_hla_matches_golden() {
    let tmp = tempfile::tempdir().unwrap();
    refbuild::build_reference(&fx("hla_subset.dat"), &fx("hla_genes.gtf"), tmp.path(), "hla")
        .unwrap();
    for f in ["hla_dna_seq.fa", "hla_rna_seq.fa", "hla_dna_coord.fa", "hla_rna_coord.fa"] {
        let got = std::fs::read_to_string(tmp.path().join(f)).unwrap();
        let want = std::fs::read_to_string(fx(&format!("golden/{f}"))).unwrap();
        assert_eq!(got, want, "mismatch in {f}");
    }
}
