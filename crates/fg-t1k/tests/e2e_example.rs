#![cfg(feature = "t1k-sys")]
use std::process::Command;

// Proves the `fg-t1k run` subprocess wrapper reproduces, byte-for-byte, the output of the
// vendored T1K oracle pipeline on the pinned KIR reference (wrapper == oracle).
#[test]
fn run_example_matches_oracle() {
    let tmp = tempfile::tempdir().unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_fg-t1k"))
        .args([
            "run",
            "-1",
            "fixtures/example/example_1.fq",
            "-2",
            "fixtures/example/example_2.fq",
            "-f",
            "fixtures/example/ref/kir_rna_seq.fa",
            "-o",
            "example",
            "--od",
            tmp.path().to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(tmp.path().join("example_genotype.tsv")).unwrap();
    let expected = std::fs::read_to_string("fixtures/example/oracle_genotype.golden.tsv").unwrap();
    assert_eq!(got, expected);
}
