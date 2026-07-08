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

    /// Post-analyze a genotype call, restricting the reference to the CALLED alleles and
    /// emitting a novel-variant `_allele.vcf` (the Rust port of `analyzer`). Also runs as the
    /// final stage of the fused `run` pipeline.
    Analyze(AnalyzeArgs),

    /// Run the full extract -> genotype -> analyze pipeline in a single fused process (the Rust
    /// port of `run-t1k`), keeping candidate reads in memory between extract and genotype.
    Run(RunArgs),
}

/// Arguments for the `run` subcommand, mirroring `run-t1k`'s flag names for the paired-end
/// FASTQ + reference-sequence-FASTA input path.
#[derive(Args, Debug)]
pub struct RunArgs {
    /// Path to the first-mate FASTQ file (paired FASTQ input; requires `-2`).
    #[arg(short = '1', value_name = "STRING")]
    pub mate1: Option<String>,

    /// Path to the second-mate FASTQ file (paired FASTQ input; requires `-1`).
    #[arg(short = '2', value_name = "STRING")]
    pub mate2: Option<String>,

    /// Path to a single-end read file (mutually exclusive with `-1`/`-2` and `-b`).
    #[arg(short = 'u', value_name = "STRING")]
    pub single: Option<String>,

    /// Path to a BAM/CRAM file (mutually exclusive with `-1`/`-2`/`-u`). NOTE: BAM input is not
    /// yet wired into the native Rust fused path -- see `crate::stages::run`'s doc comment.
    #[arg(short = 'b', value_name = "STRING")]
    pub bam: Option<String>,

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

/// BAM/CRAM extraction-mode selector (`--bam-mode`). Required whenever the
/// input is a BAM/CRAM (there is no inferred default yet). Chooses the
/// candidate-*selection* criterion:
///
/// - `alignment` (Class B): gather reads by alignment position (on-target
///   ∪ alt/decoy ∪ unaligned), then k-mer-check — T1K's `bam-extractor`.
/// - `no-alignment` (Class A): pure k-mer selection on the read sequences,
///   identical to the FASTQ path (BAM as packaged reads).
///
/// (`alignment` is only implemented for a coordinate-sorted BAM as of this
/// stage; grouped/name-sorted or stdin `alignment` input is rejected with a
/// "later release" message. `no-alignment` is fully routed -- coordinate/
/// unsorted takes a 2-pass name-map, grouped/name-sorted (including stdin)
/// takes a one-pass -- see `crate::stages::extract::run_bam_no_alignment`.)
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum BamMode {
    /// Class B: coordinate/position-based selection (`bam-extractor` parity).
    Alignment,
    /// Class A: k-mer-only selection, identical to FASTQ (BAM as reads).
    NoAlignment,
}

/// Arguments for the `extract` subcommand. Mirrors `fastq-extractor`'s flag names
/// (`FastqExtractor.cpp:12-33`'s `usage[]`) for the paired/single-end FASTQ-plus-reference-FASTA
/// candidate-extraction path (the default, when `-b` is absent), OR `bam-extractor`'s flag names
/// (`BamExtractor.cpp:16-26`'s `usage[]`) for the BAM/CRAM-plus-coord-FASTA path (when `-b` is
/// given -- see `crate::stages::extract`'s module docs for how the two modes are dispatched).
/// Barcode/single-cell (`--barcode*`/`--UMI`) input is deliberately not exposed here in either
/// mode -- see that module's docs. `-i/--input` (unified single/interleaved/paired FASTQ input,
/// with content-based format detection and `-` for stdin) is mutually exclusive with
/// `-1`/`-2`/`-u`/`-b`.
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

    /// Unified input: 1 or 2 paths (`-` = stdin). Two paths = paired FASTQ;
    /// one path = single-end or interleaved FASTQ (auto-detected by content).
    /// Format is detected from content, not the file extension. Mutually
    /// exclusive with -1/-2/-u/-b.
    #[arg(short = 'i', long = "input", value_name = "PATH", num_args = 1..=2)]
    pub input: Vec<String>,

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

    /// BAM/CRAM extraction mode. REQUIRED whenever the input is a BAM/CRAM
    /// (via `-i` content-detection or `-b`); ignored for FASTQ input. See
    /// [`BamMode`].
    #[arg(long = "bam-mode", value_enum)]
    pub bam_mode: Option<BamMode>,

    /// Prefix of the output file(s).
    #[arg(short = 'o', long = "prefix", default_value = "toassemble", value_name = "STRING")]
    pub prefix: String,

    /// Number of threads used to parallelize the per-read candidate-filter decision (both FASTQ
    /// and BAM mode). Output is byte-identical at any `-t` -- see `unum_core::extract`'s and
    /// `unum_core::bam_extract`'s module docs ("Output order == input order") for why.
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
///
/// With `--emit-metrics`, an additional `{prefix}_metrics.tsv` per-call QC + discriminative-quality
/// panel is written (a unum extension, NOT part of T1K); it does not alter `{prefix}_genotype.tsv`
/// or `{prefix}_allele.tsv`.
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

    /// Number of threads used to parallelize the per-read `get_overlaps_from_read` computation
    /// in the read-assignment loop. Output is byte-identical at any `-t` -- see
    /// `unum_core::genotyper::assign_reads_parallel`'s doc comment for why (the shared-state
    /// `assign_read` mutation, and everything downstream, always runs sequentially in a fixed
    /// order regardless of `-t`).
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

    /// Also write a `{prefix}_metrics.tsv` per-call QC + discriminative-quality panel (a unum
    /// extension, NOT part of T1K): one row per called allele with coverage-distribution stats,
    /// allele balance, null-model and best-vs-second-best genotype qualities. Default off; does not
    /// alter `{prefix}_genotype.tsv` or `{prefix}_allele.tsv`.
    #[arg(long = "emit-metrics", default_value_t = false)]
    pub emit_metrics: bool,

    /// Path to an AFND-style allele-frequency TSV (`allele<TAB>frequency<TAB>count`, optional
    /// header) enabling the opt-in Hardy-Weinberg population-frequency prior on allele selection (a
    /// unum extension, NOT part of T1K). Default off: when this flag is absent -- or present but
    /// the prior is inactive for a gene (locus absent from the table, or all candidate alleles share
    /// the same effective frequency) -- selection is byte-identical to the T1K oracle. The prior is
    /// a bounded, coverage-vanishing tie-breaker that only tips near-ties; it can never override a
    /// clear coverage margin (see `--allele-freq-weight`). CAVEAT: the HWE model assumes per-locus
    /// allele independence, which is mis-specified under the strong linkage disequilibrium and
    /// admixture typical of HLA/KIR cohorts -- acceptable for a bounded opt-in tie-breaker, but not a
    /// population-conditioned frequency model.
    #[arg(long = "allele-freq", value_name = "STRING")]
    pub allele_freq: Option<String>,

    /// Weight `w` scaling the HWE log-prior into the same weighted-read units as the coverage
    /// objective (only meaningful with `--allele-freq`). The prior can flip a call only when the
    /// candidate coverage margin is within the prior span `w * |Δ ln P_HWE|`; a coverage margin
    /// above the span can never be overridden. At `w = 0` the span is 0 and the prior is inactive
    /// everywhere (byte-identical to off). Must be finite and non-negative. Default `2.0`.
    #[arg(
        long = "allele-freq-weight",
        default_value_t = 2.0,
        value_parser = parse_finite_non_negative_f64
    )]
    pub allele_freq_weight: f64,

    /// Fixed penalty (weighted-read units) biasing the caller away from asserting a
    /// null/low-expression allele (IMGT `N`/`L`/`S`/`C`/`A`/`Q` suffix) on marginal evidence. Unlike
    /// the HWE term this penalty is name-driven (it does not need `--allele-freq`) and does not vanish
    /// with coverage, so it is bounded: accepted range `[0, 16]` (half the ~32-weighted-read
    /// no-override span, since a homozygous-null call is penalized twice), so the null penalty ALONE
    /// can never flip a coverage margin beyond that span (the HWE term stacks on top of it). Default
    /// `0.0` (inert, byte-identical to off).
    #[arg(
        long = "allele-freq-null-penalty",
        default_value_t = 0.0,
        value_parser = parse_allele_freq_null_penalty
    )]
    pub allele_freq_null_penalty: f64,
}

/// Parses and range-validates `--allele-freq-null-penalty`: a finite value in
/// `[0, Genotyper::NULL_PENALTY_MAX]`. Rejecting out-of-range values up front
/// (rather than silently clamping) makes a mis-set penalty a hard, visible error.
fn parse_allele_freq_null_penalty(s: &str) -> Result<f64, String> {
    let p: f64 = s.parse().map_err(|_| format!("`{s}` is not a valid number"))?;
    let max = unum_core::genotyper::Genotyper::NULL_PENALTY_MAX;
    if p.is_finite() && (0.0..=max).contains(&p) {
        Ok(p)
    } else {
        Err(format!("must be a finite value in [0, {max}] (the prior no-override bound); got {p}"))
    }
}

/// Parses a `--allele-freq-weight` value, rejecting negatives and non-finite inputs.
///
/// A negative weight inverts the HWE prior, and NaN/±inf make the span comparisons
/// `w * |Δ ln P_HWE|` fail closed, so only finite, non-negative values are accepted.
fn parse_finite_non_negative_f64(s: &str) -> Result<f64, String> {
    let value: f64 = s.parse().map_err(|_| format!("`{s}` is not a valid number"))?;
    if !value.is_finite() {
        return Err(format!("must be a finite number, got `{s}`"));
    }
    if value < 0.0 {
        return Err(format!("must be non-negative, got `{s}`"));
    }
    Ok(value)
}

/// Arguments for the `analyze` subcommand. Mirrors `analyzer`'s flag names (`Analyzer.cpp`'s
/// `usage[]`) for the paired/single-end aligned-read-FASTA-plus-selected-allele-list path this
/// port targets (`--barcode`/`--relaxIntronAlign`/`--alleleDigitUnits`/`--alleleDelimiter` input
/// are not exposed here -- see `crate::stages::analyze`'s module docs).
#[derive(Args, Debug)]
pub struct AnalyzeArgs {
    /// Path to the reference sequence FASTA file (e.g. `kir_rna_seq.fa`).
    #[arg(short = 'f', value_name = "STRING")]
    pub ref_seq_fasta: String,

    /// Path to the selected-alleles list file (`{prefix}_allele.tsv`, `genotype`'s output).
    #[arg(short = 'a', value_name = "STRING")]
    pub allele_file: String,

    /// Path to the first-mate aligned-read FASTA file (paired mode; requires `-2`).
    #[arg(short = '1', value_name = "STRING")]
    pub mate1: Option<String>,

    /// Path to the second-mate aligned-read FASTA file (paired mode; requires `-1`).
    #[arg(short = '2', value_name = "STRING")]
    pub mate2: Option<String>,

    /// Path to a single-end aligned-read FASTA file (mutually exclusive with `-1`/`-2`).
    #[arg(short = 'u', value_name = "STRING")]
    pub single: Option<String>,

    /// Prefix of the output file (`{prefix}_allele.vcf`).
    #[arg(short = 'o', long = "prefix", default_value = "t1k", value_name = "STRING")]
    pub prefix: String,

    /// Number of threads used to parallelize the analyzer's three per-read loops: dedup
    /// `assign_read` (via `unum_core::genotyper::assign_reads_parallel`), slot-indexed fragment
    /// assembly (`compute_read_assignment` + `set_all_read_assignments`), and
    /// `AddFragmentAlignmentInfo`. Output is byte-identical at any `-t` -- each loop is either a
    /// pure per-slot computation or a per-read-independent mutation, with every order-dependent
    /// step (coalesce, quantification, `ComputeVariant`) still running sequentially in a fixed,
    /// thread-count-independent order (see `crate::stages::analyze`'s per-loop comments and
    /// `crate::stages::genotype`'s identical slot-index rationale).
    #[arg(short = 't', default_value_t = 1)]
    pub threads: u32,

    /// Maximal number of alleles per read.
    #[arg(short = 'n', default_value_t = 2000)]
    pub max_assign_cnt: i32,

    /// Filter alignments with alignment similarity less than the specified value.
    #[arg(short = 's', default_value_t = 0.8)]
    pub similarity: f64,

    /// The maximum variant group size to call a novel variant. `-1` for no limitation.
    #[arg(long = "varMaxGroup", default_value_t = 8)]
    pub var_max_group: i32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bam_mode_parses_alignment_and_no_alignment() {
        let cli = Cli::try_parse_from([
            "unum",
            "extract",
            "-f",
            "ref.fa",
            "-i",
            "x.bam",
            "--bam-mode",
            "alignment",
        ])
        .expect("alignment should parse");
        let Commands::Extract(args) = cli.command else { panic!("expected extract subcommand") };
        assert_eq!(args.bam_mode, Some(BamMode::Alignment));

        let cli = Cli::try_parse_from([
            "unum",
            "extract",
            "-f",
            "ref.fa",
            "-i",
            "x.bam",
            "--bam-mode",
            "no-alignment",
        ])
        .expect("no-alignment should parse");
        let Commands::Extract(args) = cli.command else { panic!("expected extract subcommand") };
        assert_eq!(args.bam_mode, Some(BamMode::NoAlignment));
    }

    #[test]
    fn bam_mode_absent_is_none_and_bad_value_errors() {
        let cli = Cli::try_parse_from(["unum", "extract", "-f", "ref.fa", "-u", "reads.fq"])
            .expect("no --bam-mode should parse");
        let Commands::Extract(args) = cli.command else { panic!("expected extract subcommand") };
        assert_eq!(args.bam_mode, None);

        assert!(
            Cli::try_parse_from([
                "unum",
                "extract",
                "-f",
                "ref.fa",
                "-i",
                "x.bam",
                "--bam-mode",
                "bogus",
            ])
            .is_err(),
            "an invalid --bam-mode value must be a clap error"
        );
    }

    #[test]
    fn parse_finite_non_negative_f64_accepts_valid_weights() {
        assert_eq!(parse_finite_non_negative_f64("0"), Ok(0.0));
        assert_eq!(parse_finite_non_negative_f64("2.0"), Ok(2.0));
        assert_eq!(parse_finite_non_negative_f64("0.5"), Ok(0.5));
    }

    #[test]
    fn parse_finite_non_negative_f64_rejects_negative() {
        assert!(parse_finite_non_negative_f64("-1.0").is_err());
        assert!(parse_finite_non_negative_f64("-0.0001").is_err());
    }

    #[test]
    fn parse_finite_non_negative_f64_rejects_non_finite() {
        assert!(parse_finite_non_negative_f64("nan").is_err());
        assert!(parse_finite_non_negative_f64("inf").is_err());
        assert!(parse_finite_non_negative_f64("-inf").is_err());
    }

    #[test]
    fn parse_finite_non_negative_f64_rejects_garbage() {
        assert!(parse_finite_non_negative_f64("abc").is_err());
    }

    #[test]
    fn genotype_cli_validates_allele_freq_weight() {
        // A complete, otherwise-valid genotype invocation parses when the weight is valid,
        // so the negative-weight failure below is attributable to the value_parser, not to
        // a missing required argument.
        let base = ["unum", "genotype", "-f", "ref.fa", "-u", "reads.fa"];
        assert!(
            Cli::try_parse_from(base.iter().chain(["--allele-freq-weight", "1.5"].iter())).is_ok()
        );
        assert!(
            Cli::try_parse_from(base.iter().chain(["--allele-freq-weight", "-1.0"].iter()))
                .is_err()
        );
    }
}
