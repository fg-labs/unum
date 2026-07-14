//! `unum combine` -- port of T1K's `t1k-merge.py`.
//!
//! Reads the `-l` list of `_genotype.tsv` paths, parses each into per-gene
//! calls, combines them into one allele-by-sample abundance matrix via
//! [`unum_core::combine`], and prints the matrix to stdout. Pure post-processing
//! over TSVs -- no reads, no reference.

use crate::cli::CombineArgs;
use anyhow::{Context, Result};
use std::fs;
use std::io::{self, Write};
use unum_core::combine::{self, CombineConfig, GeneCall};

/// Derive a sample name from a genotype-file path, matching `t1k-merge.py`:
/// take the basename, drop the final `.`-extension, then strip a trailing
/// `_genotype` (so `dir/HG00476_genotype.tsv` -> `HG00476`).
fn sample_name_from_path(path: &str) -> String {
    let base = path.rsplit('/').next().unwrap_or(path);
    let parts: Vec<&str> = base.split('.').collect();
    let stem = parts[..parts.len().saturating_sub(1)].join(".");
    stem.strip_suffix("_genotype").unwrap_or(&stem).to_owned()
}

/// Run the `combine` stage: parse every listed genotype TSV, combine, and write
/// the matrix to stdout.
pub fn run(args: &CombineArgs) -> Result<()> {
    let list = fs::read_to_string(&args.input_list)
        .with_context(|| format!("reading input list {}", args.input_list))?;

    let mut samples: Vec<(String, Vec<GeneCall>)> = Vec::new();
    for path in list.lines().map(str::trim).filter(|l| !l.is_empty()) {
        let content =
            fs::read_to_string(path).with_context(|| format!("reading genotype file {path}"))?;
        let calls: Vec<GeneCall> =
            content.lines().filter_map(combine::parse_genotype_line).collect();
        samples.push((sample_name_from_path(path), calls));
    }

    let cfg = CombineConfig {
        num_alleles_per_gene: args.num_alleles_per_gene,
        min_quality: args.min_quality,
        min_total_quality: args.min_total_quality,
    };
    let matrix = combine::combine(&samples, &cfg);
    io::stdout()
        .write_all(combine::format_matrix(&matrix).as_bytes())
        .context("writing combine matrix to stdout")
}

#[cfg(test)]
mod tests {
    use super::sample_name_from_path;

    #[test]
    fn sample_name_strips_dir_extension_and_genotype_suffix() {
        assert_eq!(sample_name_from_path("out/HG00476_genotype.tsv"), "HG00476");
        assert_eq!(sample_name_from_path("HG00476_genotype.tsv"), "HG00476");
        assert_eq!(sample_name_from_path("a/b/S1.tsv"), "S1");
        // Matches the Python: a basename with no extension yields an empty stem.
        assert_eq!(sample_name_from_path("noext"), "");
    }
}
