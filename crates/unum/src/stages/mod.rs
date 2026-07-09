//! Pipeline stage implementations.
//!
//! Each stage is a self-contained, native-Rust port of the corresponding T1K
//! tool: [`build`] (reference build), [`extract`] (candidate-read extraction),
//! [`genotype`] (genotype calling), and [`analyze`] (variant calling / allele
//! VCF). [`run`] fuses extract -> genotype -> analyze into a single process.
//!
//! # Fused native Rust path (issue #28)
//!
//! [`run`] takes a FUSED path: extract's candidate reads are collected by an
//! in-memory [`unum_core::extract::InMemoryCandidateSink`] and handed straight
//! to the genotyper via [`genotype::run_with_candidate_reads`] -- **no
//! intermediate `{prefix}_candidate_*.fq` is ever written or re-read**. The
//! analyze stage then consumes the genotyper's on-disk `_aligned_*.fa`
//! (genotype -> analyze in-memory hand-off is a deliberate follow-up, out of
//! scope here).
//!
//! # Scope: FASTQ paired + single-end only; no native BAM (`-b`) yet
//!
//! `-1`/`-2` (paired) and `-u` (single-end) FASTQ input are supported. `-b`
//! (BAM/CRAM) input is NOT yet wired into the fused flow -- the extract
//! library's `bam_extract` path exists but plumbing it into `run` (coord-FASTA
//! resolution, gene intervals) is a follow-up; `-b` on `run` currently errors
//! with a clear message.
use crate::cli::{GenotypeArgs, RunArgs};
use anyhow::{Context, bail};

pub mod analyze;
pub mod build;
pub mod extract;
pub mod genotype;

/// The resolved FASTQ read source for the fused native path.
enum FastqSource<'a> {
    /// Paired: `(mate1, mate2)`.
    Paired(&'a str, &'a str),
    /// Single-end.
    Single(&'a str),
}

/// Resolves `-1`/`-2`/`-u`/`-b` into a [`FastqSource`], rejecting BAM input on the fused path and
/// the usual mutually-exclusive-flag misuse (mirroring [`extract`]'s own validation).
fn resolve_fastq_source(args: &RunArgs) -> anyhow::Result<FastqSource<'_>> {
    if args.bam.is_some() {
        bail!(
            "BAM/CRAM input (-b) is not yet supported by `unum run`; extract candidates with \
             `unum extract -b ...` and genotype separately"
        );
    }
    let paired = args.mate1.is_some() || args.mate2.is_some();
    if paired && args.single.is_some() {
        bail!("specify either -u (single-end) or -1/-2 (paired), not both");
    }
    if paired {
        let mate1 = args
            .mate1
            .as_deref()
            .context("paired input requires both -1 and -2 (got -2 without -1)")?;
        let mate2 = args
            .mate2
            .as_deref()
            .context("paired input requires both -1 and -2 (got -1 without -2)")?;
        Ok(FastqSource::Paired(mate1, mate2))
    } else if let Some(single) = args.single.as_deref() {
        Ok(FastqSource::Single(single))
    } else {
        bail!("must specify either -u (single-end) or -1/-2 (paired) read input")
    }
}

/// Runs the full extract -> genotype -> analyze pipeline for `args`, fused into a single process:
/// extract's candidate reads are kept in memory and handed straight to the genotyper (no
/// `{prefix}_candidate_*.fq` on disk), then the analyze stage runs on the genotyper's on-disk
/// `_aligned_*.fa`.
///
/// Output file names are derived exactly as `run-t1k` derives them:
/// - prefix = `{output_dir}/{prefix}` when `--od` is given, else just `{prefix}`.
/// - genotyper writes `{prefix}_allele.tsv`, `{prefix}_aligned_1.fa`, `{prefix}_aligned_2.fa`,
///   and `{prefix}_genotype.tsv`.
/// - analyzer reads those genotyper outputs and writes the final `{prefix}_allele.vcf`.
pub fn run(args: &RunArgs) -> anyhow::Result<()> {
    use std::path::Path;
    use unum_core::extract as core_extract;
    use unum_core::extract::InMemoryCandidateSink;
    use unum_core::ref_kmer_filter::RefKmerFilter;

    let prefix = resolve_prefix(args)?;
    let source_desc = resolve_fastq_source(args)?;

    // --- Extract candidates in memory (fastq-extractor port) ---
    let mut filter = RefKmerFilter::from_reference_fasta(
        Path::new(&args.ref_seq_fasta),
        self::extract::INITIAL_KMER_LENGTH,
    )
    .with_context(|| format!("loading reference FASTA {}", args.ref_seq_fasta))?;

    let (mate1_path, mate2_path): (&str, Option<&str>) = match &source_desc {
        FastqSource::Paired(m1, m2) => (m1, Some(*m2)),
        FastqSource::Single(u) => (u, None),
    };
    let mut source = core_extract::open_source(Path::new(mate1_path), mate2_path.map(Path::new))
        .context("opening read source")?;

    let mut sink = InMemoryCandidateSink::new();
    let threads = usize::try_from(args.threads).unwrap_or(usize::MAX).max(1);
    let metrics = core_extract::extract_candidates_with_threads(
        &mut source,
        &mut filter,
        core_extract::DEFAULT_REF_SEQ_SIMILARITY,
        threads,
        &mut sink,
    )
    .context("extracting candidate reads")?;
    let has_mate = mate2_path.is_some();
    eprintln!(
        "extracted {} / {} candidate {} (kmer_length={}, hit_len_required={})",
        metrics.candidates_emitted,
        metrics.total_reads,
        if has_mate { "pairs" } else { "reads" },
        metrics.kmer_length,
        metrics.hit_len_required,
    );

    let (reads1, reads2) = sink.into_reads();

    // Free the extract-phase k-mer index (~2 GB on the genomic HLA reference,
    // where `infer_kmer_length` picks k=15) before the genotyper builds its own
    // index -- otherwise two full reference indices are resident at once.
    // `filter` is not read after extraction, so dropping it is byte-neutral.
    drop(filter);

    // --- Genotype directly from the in-memory candidate reads ---
    let genotype_args = genotype_args_for(args, &prefix);
    genotype::run_with_candidate_reads(&genotype_args, reads1, reads2, has_mate)
        .context("genotyping candidate reads")?;

    // --- Analyze (native) on the genotyper's on-disk `_aligned_*.fa` ---
    let allele_tsv = format!("{prefix}_allele.tsv");
    let aligned_1 =
        if has_mate { format!("{prefix}_aligned_1.fa") } else { format!("{prefix}_aligned.fa") };
    let aligned_2 = format!("{prefix}_aligned_2.fa");
    run_analyze_native_source(args, &prefix, &allele_tsv, &aligned_1, &aligned_2, has_mate)?;

    Ok(())
}

/// Builds the [`GenotypeArgs`] the fused genotype stage runs with, mirroring `run-t1k`'s own
/// invocation (which passes only `-o`/`-t`/`-f`/`-1`/`-2` to the subprocess, so every other
/// genotyper parameter takes its default). Read-input fields are left `None` -- the fused path
/// supplies reads in memory via [`genotype::run_with_candidate_reads`], never from these paths.
fn genotype_args_for(args: &RunArgs, prefix: &str) -> GenotypeArgs {
    GenotypeArgs {
        ref_seq_fasta: args.ref_seq_fasta.clone(),
        mate1: None,
        mate2: None,
        single: None,
        prefix: prefix.to_string(),
        threads: args.threads,
        max_assign_cnt: 2000,
        similarity: 0.8,
        filter_frac: 0.15,
        filter_cov: 1.0,
        cross_gene_rate: 0.04,
        emit_metrics: false,
        allele_freq: None,
        allele_freq_weight: 2.0,
        allele_freq_null_penalty: 0.0,
    }
}

/// Runs the native (Rust) analyze stage, dispatching on `has_mate` to paired (`-1`/`-2`) or
/// single-end (`-u`) aligned-read input. Mirrors `run-t1k`'s analyzer invocation (only
/// `-o`/`-t`/`-f`/`-a` plus the aligned FASTA(s); every other analyzer parameter defaults).
fn run_analyze_native_source(
    args: &RunArgs,
    prefix: &str,
    allele_tsv: &str,
    aligned_1: &str,
    aligned_2: &str,
    has_mate: bool,
) -> anyhow::Result<()> {
    use crate::cli::AnalyzeArgs;
    let analyze_args = AnalyzeArgs {
        ref_seq_fasta: args.ref_seq_fasta.clone(),
        allele_file: allele_tsv.to_string(),
        mate1: if has_mate { Some(aligned_1.to_string()) } else { None },
        mate2: if has_mate { Some(aligned_2.to_string()) } else { None },
        single: if has_mate { None } else { Some(aligned_1.to_string()) },
        prefix: prefix.to_string(),
        threads: args.threads,
        max_assign_cnt: 2000,
        similarity: 0.8,
        var_max_group: 8,
    };
    analyze::run(&analyze_args).context("analyzing genotype call")
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
