//! Command-line interface, mirroring the flag names of the vendored `run-t1k` perl wrapper.
use clap::{Args, Parser, Subcommand, ValueHint};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "unum", version, about = "unum: a Rust port of the T1K HLA/KIR genotyper")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Build T1K reference FASTAs from a `.dat` database and a genome GTF.
    Build(BuildArgs),

    /// Extract candidate reads from FASTQ or BAM/CRAM input against a reference (the Rust port of
    /// `fastq-extractor` / `bam-extractor`).
    Extract(ExtractArgs),

    /// Call a genotype from candidate reads against a reference (the Rust port of `genotyper`).
    Genotype(GenotypeArgs),
}

/// Arguments for the `build` subcommand, mirroring `t1k-build.pl`'s `-d`/`-g`/`--prefix` flags
/// for the `.dat`-plus-GTF reference-build path (see [`crate::stages::build`]'s scope note for
/// which of `t1k-build.pl`'s flags are, and are not, implemented).
///
/// # `--od` vs. `t1k-build.pl`'s `-o`
///
/// `t1k-build.pl` itself uses `-o` for the output directory. We deliberately use `--od`
/// instead: `-o` is reserved, project-wide, for an output *file prefix* (as it is on
/// `extract`/`genotype`), not a directory.
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

/// Arguments for the `extract` subcommand. Mirrors `fastq-extractor`'s flag names
/// (`FastqExtractor.cpp:12-33`'s `usage[]`) for the paired/single-end FASTQ-plus-reference-FASTA
/// candidate-extraction path (the default, when `-b` is absent), OR `bam-extractor`'s flag names
/// (`BamExtractor.cpp:16-26`'s `usage[]`) for the BAM/CRAM-plus-coord-FASTA path (when `-b` is
/// given -- see `crate::stages::extract`'s module docs for how the two modes are dispatched).
/// Barcode/single-cell (`--barcode*`/`--UMI`) and interleaved (`-i`) input are deliberately not
/// exposed here in either mode -- see that module's docs.
#[derive(Args, Debug)]
pub struct ExtractArgs {
    /// Path to the reference sequence FASTA file (FASTQ mode) or the `_coord.fa` gene-coordinate
    /// reference (BAM mode, i.e. when `-b` is given).
    #[arg(short = 'f', value_name = "STRING")]
    pub ref_seq_fasta: String,

    /// Path to the first-mate FASTQ file (paired FASTQ-mode input; requires `-2`).
    #[arg(short = '1', value_name = "STRING")]
    pub mate1: Option<String>,

    /// Path to the second-mate FASTQ file (paired FASTQ-mode input; requires `-1`).
    #[arg(short = '2', value_name = "STRING")]
    pub mate2: Option<String>,

    /// Path to a single-end read file (FASTQ mode; mutually exclusive with `-1`/`-2` and `-b`).
    #[arg(short = 'u', value_name = "STRING")]
    pub single: Option<String>,

    /// Path to a BAM/CRAM file (switches to BAM mode; mutually exclusive with `-1`/`-2`/`-u`).
    #[arg(short = 'b', value_name = "STRING")]
    pub bam: Option<String>,

    /// BAM mode only: the flag or order of an unaligned read-pair is not ordinary (i.e. the two
    /// mates of an unaligned template are not guaranteed to be adjacent records). Mirrors
    /// `bam-extractor -u` / `BamExtractor.cpp`'s `abnormalUnalignedFlag`.
    #[arg(long = "abnormal-unmapped")]
    pub abnormal_unmapped: bool,

    /// BAM mode only: the suffix length in a read id to strip for mate matching (default: strip a
    /// trailing `/1` or `/2`). Mirrors `bam-extractor --mateIdSuffixLen`.
    #[arg(long = "mate-id-suffix-len", default_value_t = -1)]
    pub mate_id_suffix_len: i32,

    /// Prefix of the output file(s).
    #[arg(short = 'o', long = "prefix", default_value = "toassemble", value_name = "STRING")]
    pub prefix: String,

    /// Number of threads. Accepted for CLI compatibility; the extraction pass itself always
    /// runs single-threaded internally. For FASTQ mode, output is provably threadCnt-invariant --
    /// see `unum_core::extract`'s module docs ("Output order == input order"). For BAM mode,
    /// see `unum_core::bam_extract`'s module docs for why only `-t 1` oracle semantics are
    /// reproduced.
    #[arg(short = 't', default_value_t = 1)]
    pub threads: u32,

    /// FASTQ mode only: filter alignments with alignment similarity less than the specified
    /// value.
    #[arg(short = 's', default_value_t = unum_core::extract::DEFAULT_REF_SEQ_SIMILARITY)]
    pub similarity: f64,
}

/// Arguments for the `genotype` subcommand. Mirrors `genotyper`'s flag names (`Genotyper.cpp`'s
/// `usage[]`) for the paired/single-end candidate-FASTQ-plus-reference-sequence-FASTA path this
/// port targets (barcode/`-a`/`--alleleWhitelist`/`--outputReadAssignment` input are not exposed
/// here -- see `crate::stages::genotype`'s module docs).
#[derive(Args, Debug)]
pub struct GenotypeArgs {
    /// Path to the reference sequence FASTA file (e.g. `kir_rna_seq.fa`).
    #[arg(short = 'f', value_name = "STRING")]
    pub ref_seq_fasta: String,

    /// Path to the first-mate candidate-read FASTQ file (paired mode; requires `-2`).
    #[arg(short = '1', value_name = "STRING")]
    pub mate1: Option<String>,

    /// Path to the second-mate candidate-read FASTQ file (paired mode; requires `-1`).
    #[arg(short = '2', value_name = "STRING")]
    pub mate2: Option<String>,

    /// Path to a single-end candidate-read FASTQ file (mutually exclusive with `-1`/`-2`).
    #[arg(short = 'u', value_name = "STRING")]
    pub single: Option<String>,

    /// Prefix of the output files (`{prefix}_genotype.tsv`, `{prefix}_allele.tsv`).
    #[arg(short = 'o', long = "prefix", default_value = "t1k", value_name = "STRING")]
    pub prefix: String,

    /// Number of threads. Accepted for CLI compatibility -- this port always runs
    /// single-threaded internally (mirrors `Genotyper.cpp`'s `threadCnt <= 1` code path, the
    /// only one this port reproduces; both paths are deterministically byte-identical for a
    /// single genotyping run, so `-t 1` is the only value this CLI needs to accept for the
    /// end-to-end differential to compare cleanly).
    #[arg(short = 't', default_value_t = 1)]
    pub threads: u32,

    /// Maximal number of alleles per read.
    #[arg(short = 'n', default_value_t = 2000)]
    pub max_assign_cnt: i32,

    /// Filter alignments with alignment similarity less than the specified value.
    #[arg(short = 's', default_value_t = 0.8)]
    pub similarity: f64,

    /// Filter if abundance is less than the frac of the dominant allele.
    #[arg(long = "frac", default_value_t = 0.15)]
    pub filter_frac: f64,

    /// Filter genes with average coverage less than the specified value.
    #[arg(long = "cov", default_value_t = 1.0)]
    pub filter_cov: f64,

    /// The effect from other gene's expression.
    #[arg(long = "crossGeneRate", default_value_t = 0.04)]
    pub cross_gene_rate: f64,
}
