//! Binary-level smoke test for the `unum extract` subcommand: invokes the
//! actual `unum` binary on the pinned KIR example and asserts the emitted
//! candidate FASTQs are byte-identical to the frozen extraction goldens shared
//! with the library-level `unum-core/tests/golden_fastq_extract.rs`. Guards the
//! real `stages::extract::run` CLI path (the CLI defaults `-s` /
//! initial-kmer-length match the parameters the golden was captured with).
use std::process::Command;

/// A path under the workspace-level `fixtures/example/` directory.
fn example(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/example").join(rel)
}

/// A path under `unum-core`'s frozen `fastq_extract` golden directory (the same
/// goldens the library-level extraction test asserts against).
fn golden(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../unum-core/tests/golden/fastq_extract")
        .join(rel)
}

#[test]
fn extract_example_matches_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let prefix = tmp.path().join("cand");
    let status = Command::new(env!("CARGO_BIN_EXE_unum"))
        .args([
            "extract",
            "-f",
            example("ref/kir_rna_seq.fa").to_str().unwrap(),
            "-1",
            example("example_1.fq").to_str().unwrap(),
            "-2",
            example("example_2.fq").to_str().unwrap(),
            "-o",
            prefix.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "`unum extract` exited non-zero");

    for (out, gold) in [("cand_1.fq", "example_paired_1.fq"), ("cand_2.fq", "example_paired_2.fq")]
    {
        let got = std::fs::read(tmp.path().join(out)).unwrap();
        let want = std::fs::read(golden(gold)).unwrap();
        assert_eq!(got, want, "`unum extract` {out} must be byte-identical to golden {gold}");
    }
}
