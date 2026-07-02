#![forbid(unsafe_code)]

mod cli;
mod engine;
mod stages;

use clap::Parser;
use cli::{Cli, Commands};

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Run(args) => stages::run(&args),
        Commands::Build(args) => stages::build::run(&args),
        Commands::Extract(args) => stages::extract::run(&args),
        Commands::Genotype(args) => stages::genotype::run(&args),
    }
}
