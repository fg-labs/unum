//! Golden byte-identity test for the full reference-build orchestration
//! (`unum_core::refbuild::build_reference`), ported from T1K's
//! `t1k-build.pl`.
//!
//! Runs the Rust `.dat` -> seq-FASTA -> coord-FASTA pipeline end-to-end on the
//! committed KIR subset fixture and asserts all four emitted files are
//! byte-for-byte identical to the golden files produced by the vendored Perl
//! oracle (see `fixtures/refbuild/PINS.md`).

use unum_core::refbuild;

fn fx(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/refbuild").join(rel)
}

#[test]
fn build_reference_matches_golden() {
    let tmp = tempfile::tempdir().unwrap();
    refbuild::build_reference(&fx("kir_subset.dat"), &fx("kir_genes.gtf"), tmp.path(), "kir")
        .unwrap();
    for f in ["kir_dna_seq.fa", "kir_rna_seq.fa", "kir_dna_coord.fa", "kir_rna_coord.fa"] {
        let got = std::fs::read_to_string(tmp.path().join(f)).unwrap();
        let want = std::fs::read_to_string(fx(&format!("golden/{f}"))).unwrap();
        assert_eq!(got, want, "mismatch in {f}");
    }
}
