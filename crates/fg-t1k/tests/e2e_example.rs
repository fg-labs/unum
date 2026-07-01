#![cfg(feature = "t1k-sys")]
use std::process::Command;

/// Resolves a path under the workspace-level `fixtures/` directory, relative to this crate's
/// manifest directory, so the test does not depend on the process's current working directory.
fn fixture(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures").join(rel)
}

// Proves the `fg-t1k run` subprocess wrapper reproduces, byte-for-byte, the output of the
// vendored T1K oracle pipeline on the pinned KIR reference (wrapper == oracle).
#[test]
fn run_example_matches_oracle() {
    let tmp = tempfile::tempdir().unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_fg-t1k"))
        .args([
            "run",
            "-1",
            fixture("example/example_1.fq").to_str().unwrap(),
            "-2",
            fixture("example/example_2.fq").to_str().unwrap(),
            "-f",
            fixture("example/ref/kir_rna_seq.fa").to_str().unwrap(),
            "-o",
            "example",
            "--od",
            tmp.path().to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(tmp.path().join("example_genotype.tsv")).unwrap();
    let expected = std::fs::read_to_string(fixture("example/oracle_genotype.golden.tsv")).unwrap();
    assert_eq!(got, expected);
}
