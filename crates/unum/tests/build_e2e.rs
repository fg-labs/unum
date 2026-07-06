//! End-to-end test for the `unum build` subcommand: invokes the actual binary (not the
//! library function directly) on the pinned KIR fixtures and asserts all four emitted
//! reference files are byte-for-byte identical to the golden files produced by upstream T1K's
//! Perl `t1k-build.pl` (see `fixtures/refbuild/PINS.md`).
use std::process::Command;

/// Resolves a path under the workspace-level `fixtures/refbuild/` directory, relative to this
/// crate's manifest directory, so the test does not depend on the process's current working
/// directory.
fn fx(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/refbuild").join(rel)
}

#[test]
fn build_example_matches_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_unum"))
        .args([
            "build",
            "-d",
            fx("kir_subset.dat").to_str().unwrap(),
            "-g",
            fx("kir_genes.gtf").to_str().unwrap(),
            "--od",
            tmp.path().to_str().unwrap(),
            "--prefix",
            "kir",
        ])
        .status()
        .unwrap();
    assert!(status.success());

    for f in ["kir_dna_seq.fa", "kir_rna_seq.fa", "kir_dna_coord.fa", "kir_rna_coord.fa"] {
        let got = std::fs::read_to_string(tmp.path().join(f)).unwrap();
        let want = std::fs::read_to_string(fx(&format!("golden/{f}"))).unwrap();
        assert_eq!(got, want, "mismatch in {f}");
    }
}
