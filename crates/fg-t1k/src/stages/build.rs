//! Reference-build orchestration: dispatches the `build` subcommand to
//! `fg_t1k_core::refbuild::build_reference`, the Rust port of T1K's
//! `t1k-build.pl` reference-build tool.
//!
//! `t1k-build.pl` is a thin Perl orchestrator around `ParseDatFile.pl`/`AddGeneCoord.pl`, whose
//! output is captured (byte-for-byte) by the golden-file fixtures under `fixtures/refbuild/`
//! (see `crates/fg-t1k/tests/build_e2e.rs`). This stage is pure Rust.
use crate::cli::BuildArgs;

/// Builds the four T1K reference files for `args`, matching `t1k-build.pl`'s
/// `-d`/`-g`/`--prefix` invocation path (see [`crate::cli::BuildArgs`]'s doc comment for the one
/// deliberate flag-naming divergence, `--od` vs. `t1k-build.pl`'s own `-o`).
pub fn run(args: &BuildArgs) -> anyhow::Result<()> {
    fg_t1k_core::refbuild::build_reference(&args.dat, &args.gtf, &args.output_dir, &args.prefix)
}
