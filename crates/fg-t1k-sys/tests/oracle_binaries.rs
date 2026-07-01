#![cfg(feature = "t1k-sys")]
use fg_t1k_sys::oracle::{OracleStage, binary_path};
use std::process::Command;

#[test]
fn genotyper_binary_runs() {
    let bin = binary_path(OracleStage::Genotyper);
    assert!(bin.exists(), "oracle binary not built: {bin:?}");
    let out = Command::new(&bin).output().expect("spawn genotyper");
    let text = String::from_utf8_lossy(&out.stderr);
    // Genotyper.cpp's `usage[]` literal starts with "./genotyper [OPTIONS]:" and never
    // contains the literal word "usage"/"Usage" (only the C variable is named `usage`),
    // so check for the actual printed marker instead of the word.
    assert!(
        text.to_lowercase().contains("[options]"),
        "no usage text: {text}"
    );
}
