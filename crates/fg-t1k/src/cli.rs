//! Command-line interface, mirroring the flag names of the vendored `run-t1k` perl wrapper.
use clap::{Args, Parser, Subcommand, ValueHint};
use std::path::PathBuf;

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

    /// Build T1K reference FASTAs from a `.dat` database and a genome GTF.
    Build(BuildArgs),

    /// Extract candidate reads from FASTQ input against a reference (the Rust port of
    /// `fastq-extractor`). NOT YET wired into the `run`/`--engine` strangler router -- see
    /// `crate::stages::extract`'s module docs for that follow-up.
    Extract(ExtractArgs),
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

/// Arguments for the `build` subcommand, mirroring `t1k-build.pl`'s `-d`/`-g`/`--prefix` flags
/// for the `.dat`-plus-GTF reference-build path (see [`crate::stages::build`]'s scope note for
/// which of `t1k-build.pl`'s flags are, and are not, implemented).
///
/// # `--od` vs. `t1k-build.pl`'s `-o`
///
/// `t1k-build.pl` itself uses `-o` for the output directory. We deliberately use `--od`
/// instead, for consistency with the `run` subcommand's own `--od` (output-directory)
/// convention established in Phase 0 — `-o` is reserved, project-wide, for an output
/// *file prefix* (as it is on `run`), not a directory.
#[derive(Args, Debug)]
pub struct BuildArgs {
    /// Path to the EMBL-style IPD-IMGT/HLA or IPD-KIR `.dat` database.
    #[arg(short = 'd', long = "dat", value_name = "FILE", value_hint = ValueHint::FilePath)]
    pub dat: PathBuf,

    /// Path to the genome GTF gene-annotation file.
    #[arg(short = 'g', long = "gtf", value_name = "FILE", value_hint = ValueHint::FilePath)]
    pub gtf: PathBuf,

    /// The directory for output files.
    #[arg(long = "od", value_name = "DIR", value_hint = ValueHint::DirPath)]
    pub output_dir: PathBuf,

    /// Prefix of output files.
    #[arg(long = "prefix", value_name = "STRING")]
    pub prefix: String,
}

/// Arguments for the `extract` subcommand, mirroring `fastq-extractor`'s flag names
/// (`FastqExtractor.cpp:12-33`'s `usage[]`) for the paired/single-end FASTQ-plus-reference-FASTA
/// candidate-extraction path. Barcode/single-cell (`--barcode*`) and interleaved (`-i`) input are
/// deliberately not exposed here -- see `crate::stages::extract`'s module docs.
#[derive(Args, Debug)]
pub struct ExtractArgs {
    /// Path to the reference sequence FASTA file.
    #[arg(short = 'f', value_name = "STRING")]
    pub ref_seq_fasta: String,

    /// Path to the first-mate FASTQ file (paired input; requires `-2`).
    #[arg(short = '1', value_name = "STRING")]
    pub mate1: Option<String>,

    /// Path to the second-mate FASTQ file (paired input; requires `-1`).
    #[arg(short = '2', value_name = "STRING")]
    pub mate2: Option<String>,

    /// Path to a single-end read file (mutually exclusive with `-1`/`-2`).
    #[arg(short = 'u', value_name = "STRING")]
    pub single: Option<String>,

    /// Prefix of the output file(s).
    #[arg(short = 'o', long = "prefix", default_value = "toassemble", value_name = "STRING")]
    pub prefix: String,

    /// Number of threads. Accepted for CLI compatibility; the extraction pass itself always
    /// runs single-threaded internally, since its output is provably threadCnt-invariant -- see
    /// `fg_t1k_core::extract`'s module docs ("Output order == input order").
    #[arg(short = 't', default_value_t = 1)]
    pub threads: u32,

    /// Filter alignments with alignment similarity less than the specified value.
    #[arg(short = 's', default_value_t = fg_t1k_core::extract::DEFAULT_REF_SEQ_SIMILARITY)]
    pub similarity: f64,
}
