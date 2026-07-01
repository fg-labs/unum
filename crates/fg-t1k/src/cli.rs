//! Command-line interface, mirroring the flag names of the vendored `run-t1k` perl wrapper.
use clap::{Args, Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(name = "fg-t1k", version, about = "Strangler-fig Rust port of T1K")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Run the T1K genotyping pipeline (extract -> genotype -> analyze).
    Run(RunArgs),
}

/// Arguments for the `run` subcommand, mirroring `run-t1k`'s flag names for the paired-end
/// FASTQ + reference-sequence-FASTA input path.
#[derive(Args, Debug)]
pub struct RunArgs {
    /// Path to the first-mate FASTQ file.
    #[arg(short = '1', value_name = "STRING")]
    pub mate1: String,

    /// Path to the second-mate FASTQ file.
    #[arg(short = '2', value_name = "STRING")]
    pub mate2: String,

    /// Path to the reference sequence FASTA file.
    #[arg(short = 'f', value_name = "STRING")]
    pub ref_seq_fasta: String,

    /// Path to the gene coordinate file (only required for BAM input; unused on this path).
    #[arg(short = 'c', value_name = "STRING")]
    pub ref_coord_fasta: Option<String>,

    /// Prefix of output files.
    #[arg(short = 'o', value_name = "STRING")]
    pub prefix: String,

    /// The directory for output files.
    #[arg(long = "od", value_name = "STRING")]
    pub output_dir: Option<String>,

    /// Number of threads.
    #[arg(short = 't', default_value_t = 1)]
    pub threads: u32,

    /// Per-stage engine override: `STAGE=cpp|rust` (repeatable).
    #[arg(long = "engine", value_name = "STAGE=cpp|rust")]
    pub engine: Vec<String>,
}
