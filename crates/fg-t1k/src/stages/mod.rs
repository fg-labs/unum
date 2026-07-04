//! Pipeline stage orchestration: extract -> genotype -> analyze.
//!
//! Mirrors the vendored `run-t1k` perl wrapper's paired-end-FASTQ-plus-reference-sequence-FASTA
//! path (no `--preset`, no BAM input): each stage is dispatched, per the `--engine` overrides, to
//! either the vendored C++ oracle binary (subprocess, requires the `t1k-sys` feature) or the Rust
//! port (not yet implemented in Phase 0).
use crate::cli::RunArgs;
use crate::engine::{Engine, EngineOverrides};
use anyhow::Context;

pub mod build;
pub mod extract;

/// Runs the full extract -> genotype -> analyze pipeline for `args`.
///
/// Output file names are derived exactly as `run-t1k` derives them:
/// - prefix = `{output_dir}/{prefix}` when `--od` is given, else just `{prefix}`.
/// - extractor_prefix = `{prefix}_candidate`, producing `{extractor_prefix}_1.fq` / `_2.fq`.
/// - genotyper writes `{prefix}_allele.tsv`, `{prefix}_aligned_1.fa`, `{prefix}_aligned_2.fa`,
///   and `{prefix}_genotype.tsv`.
/// - analyzer reads those genotyper outputs and writes the final `{prefix}_genotype.tsv`
///   (post-analysis is not skipped in this path, matching `run-t1k`'s default).
pub fn run(args: &RunArgs) -> anyhow::Result<()> {
    let overrides = EngineOverrides::parse(&args.engine)?;
    let prefix = resolve_prefix(args)?;

    let extractor_prefix = format!("{prefix}_candidate");
    let candidate_1 = format!("{extractor_prefix}_1.fq");
    let candidate_2 = format!("{extractor_prefix}_2.fq");

    match overrides.engine_for("extract") {
        Engine::Rust => anyhow::bail!("rust engine for stage 'extract' is not yet implemented"),
        Engine::Cpp => cpp::run_extract(args, &extractor_prefix)?,
    }

    match overrides.engine_for("genotype") {
        Engine::Rust => anyhow::bail!("rust engine for stage 'genotype' is not yet implemented"),
        Engine::Cpp => cpp::run_genotype(args, &prefix, &candidate_1, &candidate_2)?,
    }

    let allele_tsv = format!("{prefix}_allele.tsv");
    let aligned_1 = format!("{prefix}_aligned_1.fa");
    let aligned_2 = format!("{prefix}_aligned_2.fa");

    match overrides.engine_for("analyze") {
        Engine::Rust => anyhow::bail!("rust engine for stage 'analyze' is not yet implemented"),
        Engine::Cpp => cpp::run_analyze(args, &prefix, &allele_tsv, &aligned_1, &aligned_2)?,
    }

    Ok(())
}

/// Combines `-o`/`--od` into the shared output-file prefix, creating `--od`'s directory if
/// needed — matching `run-t1k`'s `make_path($outputDirectory) ... $prefix = "$outputDirectory/$prefix"`.
fn resolve_prefix(args: &RunArgs) -> anyhow::Result<String> {
    match &args.output_dir {
        Some(dir) => {
            std::fs::create_dir_all(dir)
                .with_context(|| format!("creating output directory {dir}"))?;
            Ok(format!("{dir}/{}", args.prefix))
        }
        None => Ok(args.prefix.clone()),
    }
}

/// C++-oracle stage implementations. The `t1k-sys`-gated module shells out to the vendored
/// binaries; the fallback module (default build, no C++ compiled in) returns a clear error
/// instead of referencing any oracle-binary code.
mod cpp {
    #[cfg(feature = "t1k-sys")]
    mod imp {
        use crate::cli::RunArgs;
        use anyhow::{Context, Result, ensure};
        use fg_t1k_sys::oracle::{OracleStage, binary_path};
        use std::process::Command;

        pub fn run_extract(args: &RunArgs, extractor_prefix: &str) -> Result<()> {
            let status = Command::new(binary_path(OracleStage::FastqExtractor))
                .args(["-t", &args.threads.to_string()])
                .args(["-f", &args.ref_seq_fasta])
                .args(["-o", extractor_prefix])
                .args(["-1", &args.mate1])
                .args(["-2", &args.mate2])
                .status()
                .context("spawning fastq-extractor")?;
            ensure!(status.success(), "fastq-extractor failed: {status}");
            Ok(())
        }

        pub fn run_genotype(
            args: &RunArgs,
            prefix: &str,
            candidate_1: &str,
            candidate_2: &str,
        ) -> Result<()> {
            let status = Command::new(binary_path(OracleStage::Genotyper))
                .args(["-o", prefix])
                .args(["-t", &args.threads.to_string()])
                .args(["-f", &args.ref_seq_fasta])
                .args(["-1", candidate_1])
                .args(["-2", candidate_2])
                .status()
                .context("spawning genotyper")?;
            ensure!(status.success(), "genotyper failed: {status}");
            Ok(())
        }

        pub fn run_analyze(
            args: &RunArgs,
            prefix: &str,
            allele_tsv: &str,
            aligned_1: &str,
            aligned_2: &str,
        ) -> Result<()> {
            let status = Command::new(binary_path(OracleStage::Analyzer))
                .args(["-o", prefix])
                .args(["-t", &args.threads.to_string()])
                .args(["-f", &args.ref_seq_fasta])
                .args(["-a", allele_tsv])
                .args(["-1", aligned_1])
                .args(["-2", aligned_2])
                .status()
                .context("spawning analyzer")?;
            ensure!(status.success(), "analyzer failed: {status}");
            Ok(())
        }
    }

    #[cfg(not(feature = "t1k-sys"))]
    mod imp {
        use crate::cli::RunArgs;
        use anyhow::{Result, bail};

        const MISSING_FEATURE_MSG: &str =
            "the C++ oracle is not compiled in; rebuild with `--features t1k-sys`";

        pub fn run_extract(_args: &RunArgs, _extractor_prefix: &str) -> Result<()> {
            bail!(MISSING_FEATURE_MSG)
        }

        pub fn run_genotype(
            _args: &RunArgs,
            _prefix: &str,
            _candidate_1: &str,
            _candidate_2: &str,
        ) -> Result<()> {
            bail!(MISSING_FEATURE_MSG)
        }

        pub fn run_analyze(
            _args: &RunArgs,
            _prefix: &str,
            _allele_tsv: &str,
            _aligned_1: &str,
            _aligned_2: &str,
        ) -> Result<()> {
            bail!(MISSING_FEATURE_MSG)
        }
    }

    pub use imp::{run_analyze, run_extract, run_genotype};
}
