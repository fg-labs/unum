//! Reference-build orchestration: dispatches the `build` subcommand to
//! `fg_t1k_core::refbuild::build_reference`, T1K's Rust reference-build port
//! (`vendor/t1k/t1k-build.pl`).
//!
//! Unlike the `run` pipeline's per-stage `--engine cpp|rust` toggle, there is no C++ oracle
//! binary for reference-build to strangle: `t1k-build.pl` itself is a thin Perl orchestrator
//! around `ParseDatFile.pl`/`AddGeneCoord.pl`, already fully captured (byte-for-byte) by the
//! golden-file fixtures under `fixtures/refbuild/`. This stage is therefore always pure Rust.
use crate::cli::BuildArgs;

/// Builds the four T1K reference files for `args`, matching `t1k-build.pl`'s
/// `-d`/`-g`/`--prefix` invocation path (see [`crate::cli::BuildArgs`]'s doc comment for the one
/// deliberate flag-naming divergence, `--od` vs. `t1k-build.pl`'s own `-o`).
pub fn run(args: &BuildArgs) -> anyhow::Result<()> {
    fg_t1k_core::refbuild::build_reference(&args.dat, &args.gtf, &args.output_dir, &args.prefix)
}
