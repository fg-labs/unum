//! `unum copy-number` -- port of T1K's `t1k-copynumber.py`.
//!
//! Reads the `-g` genotype TSV, infers per-allele copy number via
//! [`unum_core::copy_number`], and prints the result to stdout. Pure
//! post-processing over one genotype file -- no reads, no reference.

use crate::cli::CopyNumberArgs;
use anyhow::{Context, Result};
use std::fs;
use std::io::{self, Write};
use unum_core::copy_number::{self, CopyNumberConfig, GeneAlleles};

/// Run the `copy-number` stage: parse the genotype file, infer copy numbers,
/// and write the per-gene table to stdout.
pub fn run(args: &CopyNumberArgs) -> Result<()> {
    let content = fs::read_to_string(&args.genotype)
        .with_context(|| format!("reading genotype file {}", args.genotype))?;
    let genes: Vec<GeneAlleles> =
        content.lines().filter_map(copy_number::parse_genotype_line).collect();

    let nomissing_genes = if args.nomissing.is_empty() {
        Vec::new()
    } else {
        args.nomissing.split(',').map(str::to_owned).collect()
    };
    let cfg = CopyNumberConfig {
        nomissing_genes,
        upper_quantile: args.upper_quantile,
        lower_quantile: args.lower_quantile,
        adjust_var: args.adjust_var,
        min_quality: args.min_quality,
    };
    cfg.validate()?;

    let rows = copy_number::infer_copy_numbers(&genes, &cfg)?;
    io::stdout()
        .write_all(copy_number::format_rows(&rows).as_bytes())
        .context("writing copy-number table to stdout")
}
