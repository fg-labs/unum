//! Functional coverage for the opt-in `{prefix}_metrics.tsv` QC +
//! discriminative-quality panel (issues #19/#30). Runs the real `unum
//! genotype` binary on the pinned KIR example fixture (pure-Rust path, no
//! oracle needed) and asserts:
//!
//! - with `--emit-metrics`, a well-formed `{prefix}_metrics.tsv` is written
//!   (exact header, 15 columns per row, one row per called allele, values in
//!   their documented ranges);
//! - without the flag, NO `_metrics.tsv` is written;
//! - the byte-frozen `_genotype.tsv`/`_allele.tsv` are byte-identical with vs
//!   without the flag (the metrics panel is purely additive).

use std::process::Command;

const HEADER: &str = "gene\tallele_rank\tallele\tabundance\tbalance_ratio\tcov_min\tcov_p10\t\
cov_median\tfrac_bases_covered\tmissing_cov\tgt_quality\trunnerup_abundance\tq_gap\tq_min\t\
locus_gq_min";

/// Resolves a path under the workspace-level `fixtures/` directory relative to
/// this crate's manifest, independent of the process's working directory.
fn fixture(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures").join(rel)
}

/// Runs `unum genotype` on the KIR example fixture into `dir/out`, optionally
/// with `--emit-metrics`.
fn run_genotype(dir: &std::path::Path, emit_metrics: bool) {
    let prefix = dir.join("out");
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_unum"));
    cmd.args([
        "genotype",
        "-f",
        fixture("example/ref/kir_rna_seq.fa").to_str().unwrap(),
        "-1",
        fixture("example/example_1.fq").to_str().unwrap(),
        "-2",
        fixture("example/example_2.fq").to_str().unwrap(),
        "-o",
        prefix.to_str().unwrap(),
        "-t",
        "1",
    ]);
    if emit_metrics {
        cmd.arg("--emit-metrics");
    }
    let status = cmd.status().unwrap();
    assert!(status.success(), "unum genotype failed");
}

#[test]
fn emit_metrics_writes_well_formed_panel() {
    let tmp = tempfile::tempdir().unwrap();
    run_genotype(tmp.path(), true);

    let metrics_path = tmp.path().join("out_metrics.tsv");
    let content = std::fs::read_to_string(&metrics_path).expect("metrics file must exist");
    let mut lines = content.lines();

    assert_eq!(lines.next().unwrap(), HEADER, "unexpected metrics header");

    let mut row_cnt = 0usize;
    for line in lines {
        row_cnt += 1;
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 15, "each metrics row must have 15 columns: {line}");

        let allele_rank: i32 = cols[1].parse().unwrap();
        assert!(allele_rank == 0 || allele_rank == 1, "allele_rank must be 0 or 1");

        let balance_ratio: f64 = cols[4].parse().unwrap();
        assert!(balance_ratio > 0.0 && balance_ratio <= 1.0, "balance_ratio in (0,1]: {line}");

        let cov_min: i32 = cols[5].parse().unwrap();
        let cov_p10: i32 = cols[6].parse().unwrap();
        let cov_median: i32 = cols[7].parse().unwrap();
        assert!(cov_min >= 0, "cov_min >= 0");
        assert!(cov_p10 >= cov_min, "cov_p10 >= cov_min");
        assert!(cov_median >= cov_p10, "cov_median >= cov_p10");

        let frac: f64 = cols[8].parse().unwrap();
        assert!((0.0..=1.0).contains(&frac), "frac_bases_covered in [0,1]: {line}");

        let gt_quality: i32 = cols[10].parse().unwrap();
        let q_gap: i32 = cols[12].parse().unwrap();
        let q_min: i32 = cols[13].parse().unwrap();
        assert!((0..=60).contains(&gt_quality), "gt_quality in [0,60]: {line}");
        assert!((0..=60).contains(&q_gap), "q_gap in [0,60]: {line}");
        assert!((0..=60).contains(&q_min), "q_min in [0,60]: {line}");
        assert_eq!(q_min, gt_quality.min(q_gap), "q_min == min(gt_quality, q_gap): {line}");

        let locus_gq_min: i32 = cols[14].parse().unwrap();
        assert!(locus_gq_min <= gt_quality, "locus_gq_min <= this row's gt_quality: {line}");
    }
    assert!(row_cnt > 0, "expected at least one called-allele row");
}

#[test]
fn no_metrics_file_without_flag() {
    let tmp = tempfile::tempdir().unwrap();
    run_genotype(tmp.path(), false);
    assert!(
        !tmp.path().join("out_metrics.tsv").exists(),
        "no _metrics.tsv should be written without --emit-metrics"
    );
}

#[test]
fn frozen_outputs_byte_identical_with_and_without_flag() {
    let with_flag = tempfile::tempdir().unwrap();
    let without_flag = tempfile::tempdir().unwrap();
    run_genotype(with_flag.path(), true);
    run_genotype(without_flag.path(), false);

    for name in ["out_genotype.tsv", "out_allele.tsv"] {
        let a = std::fs::read(with_flag.path().join(name)).unwrap();
        let b = std::fs::read(without_flag.path().join(name)).unwrap();
        assert_eq!(a, b, "{name} must be byte-identical with vs without --emit-metrics");
    }
}
