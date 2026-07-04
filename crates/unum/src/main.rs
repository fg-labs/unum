#![forbid(unsafe_code)]

mod cli;
mod stages;

use clap::Parser;
use cli::{Cli, Commands};

/// Use mimalloc as the global allocator. The genotyper's hot read-assignment
/// loop allocates many short-lived per-read/per-overlap buffers across rayon
/// workers; the system allocator (macOS libmalloc) serializes these on
/// per-size-class locks, capping multi-thread scaling. mimalloc's per-thread
/// heaps remove that contention. Allocator choice does not change any output
/// (byte-identical). The `#[global_allocator]` static is a safe declaration --
/// mimalloc's `GlobalAlloc` impl (with its own unsafe) lives inside the crate,
/// so `#![forbid(unsafe_code)]` on this crate still holds.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Build(args) => stages::build::run(&args),
        Commands::Extract(args) => stages::extract::run(&args),
        Commands::Genotype(args) => stages::genotype::run(&args),
        Commands::Analyze(args) => stages::analyze::run(&args),
    }
}
