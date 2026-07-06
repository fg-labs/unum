//! Binary-level end-to-end smoke test for the `unum genotype` subcommand:
//! invokes the actual `unum` binary (not the library) on the pinned KIR
//! example candidate reads + reference, and asserts the emitted
//! `_genotype.tsv` is byte-identical to the pinned Phase-0 golden
//! (`fixtures/example/oracle_genotype.golden.tsv`).
//!
//! This guards the real `stages::genotype::run` CLI path end to end.
//! `unum-core/tests/golden_genotype_e2e.rs` freezes the same byte-identity at
//! the library level, but it must re-implement the driver (the `unum-core`
//! crate cannot depend on the `unum` CLI binary), so it cannot catch drift
//! between that re-implementation and the actual `stages::genotype::run`.
//! This test closes that gap the way `build_e2e.rs` does for `build`.
//!
//! # `#[ignore]` — run on demand
//!
//! Genotyping the 1000-pair KIR example through the **debug** (`opt-level = 0`)
//! `unum` binary takes ~90-120s (the alignment/EM inner loops are pathologically
//! slow unoptimized), so this test is `#[ignore]`d to keep default `cargo
//! test`/CI fast. Correctness/byte-identity of the genotyping is already
//! covered on every run by the fast library-level `golden_genotype_e2e`; this
//! test's job is to catch drift between that re-implementation and the real
//! `stages::genotype::run` CLI wiring, which is a rare event worth an on-demand
//! check: run it with `cargo test -p unum --test genotype_e2e -- --ignored`
//! (or against a release binary for speed).
use std::process::Command;

/// Resolves a path under the workspace-level `fixtures/example/` directory,
/// relative to this crate's manifest directory (CWD-independent).
fn fx(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/example").join(rel)
}

#[test]
#[ignore = "runs the genotyper through the debug binary (~90-120s at opt-level 0); \
            run with --ignored. Correctness is covered by the fast library-level \
            golden_genotype_e2e; this guards CLI-driver drift only."]
fn genotype_example_matches_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let prefix = tmp.path().join("out");
    let status = Command::new(env!("CARGO_BIN_EXE_unum"))
        .args([
            "genotype",
            "-f",
            fx("ref/kir_rna_seq.fa").to_str().unwrap(),
            "-1",
            fx("example_1.fq").to_str().unwrap(),
            "-2",
            fx("example_2.fq").to_str().unwrap(),
            "-o",
            prefix.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "`unum genotype` exited non-zero");

    let got = std::fs::read_to_string(format!("{}_genotype.tsv", prefix.display())).unwrap();
    let want = std::fs::read_to_string(fx("oracle_genotype.golden.tsv")).unwrap();
    assert_eq!(got, want, "`unum genotype` _genotype.tsv must match the pinned golden");
}
