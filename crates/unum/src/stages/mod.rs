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
//! # Scope: FASTQ paired/single-end, and BAM/CRAM (no-alignment + alignment)
//!
//! `-1`/`-2` (paired) and `-u` (single-end) FASTQ input are supported via [`run`]'s
//! original in-memory extract -> genotype fusion. `-b` (BAM/CRAM) input is ALSO fused
//! in-memory, via [`run_bam_fused`]: `--bam-mode no-alignment` (coordinate/unsorted 2-pass
//! or grouped/name-sorted one-pass, matching `extract::run_bam_no_alignment`'s routing)
//! and `--bam-mode alignment` (coordinate-sorted 2-pass or grouped/name-sorted one-pass,
//! matching `extract::run_bam_alignment`'s routing) are both wired through
//! [`unum_core::bam_extract`] into an [`unum_core::extract::InMemoryCandidateSink`], same
//! as the FASTQ path. `-r` is accepted (and required for CRAM); for BAM it is accepted and
//! ignored, matching `extract`'s own CRAM handling. Unlike the FASTQ path, `has_mate` for
//! BAM/CRAM input is derived from the extractor's own [`unum_core::bam_extract::BamExtractMetrics::single_end`]
//! (there is no `-2` flag for BAM). `-o` remains required for both input kinds -- no T1K
//! basename inference.
use crate::cli::{BamMode, GenotypeArgs, RunArgs};
use anyhow::{Context, bail, ensure};
use unum_core::bam_extract::BamExtractMetrics;
use unum_core::extract::InMemoryCandidateSink;

pub mod analyze;
pub mod build;
pub mod extract;
pub mod genotype;

/// `extract`'s `--abnormal-unmapped` flag isn't exposed by [`RunArgs`] (mirrors
/// `genotype_args_for`'s hardcoded defaults): [`run_bam_alignment_fused`] always passes the
/// same default `ExtractArgs` itself defaults to when the flag is omitted.
const RUN_BAM_ABNORMAL_UNALIGNED_FLAG: bool = false;

/// `extract`'s `--mate-id-suffix-len` flag isn't exposed by [`RunArgs`] (mirrors
/// [`RUN_BAM_ABNORMAL_UNALIGNED_FLAG`]'s doc comment): [`run_bam_no_alignment_fused`]/
/// [`run_bam_alignment_fused`] always pass the same default `ExtractArgs` itself defaults to
/// when the flag is omitted.
const RUN_BAM_MATE_ID_SUFFIX_LEN: i32 = -1;

/// The resolved FASTQ read source for the fused native path.
enum FastqSource<'a> {
    /// Paired: `(mate1, mate2)`.
    Paired(&'a str, &'a str),
    /// Single-end.
    Single(&'a str),
}

/// Resolves `-1`/`-2`/`-u` into a [`FastqSource`], enforcing the usual mutually-exclusive-flag
/// misuse (mirroring [`extract`]'s own validation). Never called when `-b` is given -- [`run`]
/// branches to [`run_bam_fused`] before reaching this function.
fn resolve_fastq_source(args: &RunArgs) -> anyhow::Result<FastqSource<'_>> {
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
/// - prefix = `{output_dir}/{prefix}` when `--output-dir` is given, else just `{prefix}`.
/// - genotyper writes `{prefix}_allele.tsv`, `{prefix}_aligned_1.fa`, `{prefix}_aligned_2.fa`,
///   and `{prefix}_genotype.tsv`.
/// - analyzer reads those genotyper outputs and writes the final `{prefix}_allele.vcf`.
pub fn run(args: &RunArgs) -> anyhow::Result<()> {
    use std::path::Path;
    use unum_core::extract as core_extract;
    use unum_core::ref_kmer_filter::RefKmerFilter;

    let prefix = resolve_prefix(args)?;

    if let Some(bam) = args.bam.as_deref() {
        // Mirrors `extract::resolve_extract_input`'s `-b` arm: `-b` is mutually exclusive
        // with the FASTQ flags, so a caller combining them gets a clear error instead of
        // `run_bam_fused` silently winning and the FASTQ flags being discarded.
        ensure!(
            args.mate1.is_none() && args.mate2.is_none() && args.single.is_none(),
            "-b (BAM/CRAM mode) is mutually exclusive with -1/-2/-u"
        );
        return run_bam_fused(args, &prefix, bam);
    }

    // `-c` (gene coordinate FASTA) is only meaningful for `-b --bam-mode alignment`; its doc
    // comment declares it rejected for FASTQ input. Without this guard the FASTQ path below would
    // silently ignore it, mirroring `run_bam_no_alignment_fused`'s reject rather than falling
    // through. (`-r`/`--bam-mode` are documented as *ignored* for FASTQ, so they are not rejected.)
    ensure!(
        args.ref_coord_fasta.is_none(),
        "-c (coord FASTA) applies only to run -b --bam-mode alignment, not FASTQ input"
    );

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

/// Runs the BAM/CRAM branch of the fused `run` pipeline: extracts candidates from `bam` via
/// [`unum_core::bam_extract`] straight into an in-memory
/// [`unum_core::extract::InMemoryCandidateSink`] (no `{prefix}_candidate_*.fq` on disk), then
/// hands them to the genotyper and analyzer exactly as [`run`]'s FASTQ branch does.
///
/// Routes on `--bam-mode` (`args.bam_mode`, required) and, for `no-alignment`, on the BAM's own
/// `@HD` sort order -- EXACTLY mirroring `extract::run_bam_no_alignment`/
/// `extract::run_bam_alignment`'s dispatch:
/// - `no-alignment` on a coordinate-sorted/unsorted file: the seekable 2-pass
///   ([`bam_extract::extract_from_bam_no_alignment`]).
/// - `no-alignment` on a grouped/name-sorted file: the one-pass
///   ([`bam_extract::extract_from_bam_no_alignment_grouped`]), fed a factory that always
///   returns a fresh [`InMemoryCandidateSink`] (there is no output-naming decision to make in
///   the fused path, unlike the CLI `extract` subcommand's `FastqFileSink`).
/// - `alignment` on a coordinate-sorted file: [`bam_extract::extract_from_bam_with_threads`],
///   fed gene intervals resolved from `-c` via [`bam_extract::parse_coord_fa`]/
///   [`bam_extract::build_genes`].
/// - `alignment` on a grouped/name-sorted file: [`bam_extract::extract_from_bam_alignment_grouped`],
///   fed the SAME gene intervals plus a factory that always returns a fresh
///   [`InMemoryCandidateSink`] (`run`'s `-b` is always a file/path -- there is no `-i -`
///   stdin route into `run`, so this is always called with `seekable = true`).
/// - `alignment` on an unsorted file: rejected (mirrors `extract::run_bam_alignment`).
///
/// Unlike the FASTQ branch, `has_mate` is derived from the extractor's own
/// `BamExtractMetrics::single_end` (there is no `-2` flag for BAM/CRAM input) -- see this
/// module's doc comment.
///
/// # Errors
///
/// Returns an error if: `--bam-mode` is not given; the input cannot be content-sniffed/opened
/// (including a missing `-r` for CRAM, or a mistaken FASTQ/FASTA `-b`); `--bam-mode alignment`
/// is given an unsorted BAM/CRAM; `--bam-mode alignment` is given without `-c`; `--bam-mode
/// no-alignment` is given `-c` (rejected -- see [`run_bam_no_alignment_fused`]); or the
/// underlying extraction/genotyping/analysis itself fails.
fn run_bam_fused(args: &RunArgs, prefix: &str, bam: &str) -> anyhow::Result<()> {
    use unum_core::read_input::{InputSpec, OpenedInput, open_input};

    let mode = args.bam_mode.context("run -b requires --bam-mode {alignment|no-alignment}")?;

    // Content-sniff to distinguish BAM from CRAM (and reject a mistaken FASTQ/FASTA `-b`),
    // mirroring `extract::resolve_extract_input`'s `-b` arm.
    let (opened, _fmt) = open_input(&InputSpec::Path(std::path::PathBuf::from(bam)))
        .with_context(|| format!("opening input {bam}"))?;
    let is_cram = match opened {
        // SAM is read by htslib exactly like BAM (no reference needed).
        OpenedInput::Sam | OpenedInput::Bam => false,
        OpenedInput::Cram => true,
        OpenedInput::Fastq(_) => {
            bail!("-b expects a SAM/BAM/CRAM file, but {bam} is FASTQ/FASTA input")
        }
    };
    let reference = require_reference_for_cram_run(args, is_cram)?;
    let threads = usize::try_from(args.threads).unwrap_or(usize::MAX).max(1);

    let (single_end, sink) = match mode {
        BamMode::NoAlignment => run_bam_no_alignment_fused(args, bam, reference, threads)?,
        BamMode::Alignment => run_bam_alignment_fused(args, bam, reference, threads)?,
    };

    let has_mate = !single_end;
    let (reads1, reads2) = sink.into_reads();

    // --- Genotype directly from the in-memory candidate reads ---
    let genotype_args = genotype_args_for(args, prefix);
    genotype::run_with_candidate_reads(&genotype_args, reads1, reads2, has_mate)
        .context("genotyping candidate reads")?;

    // --- Analyze (native) on the genotyper's on-disk `_aligned_*.fa` ---
    let allele_tsv = format!("{prefix}_allele.tsv");
    let aligned_1 =
        if has_mate { format!("{prefix}_aligned_1.fa") } else { format!("{prefix}_aligned.fa") };
    let aligned_2 = format!("{prefix}_aligned_2.fa");
    run_analyze_native_source(args, prefix, &allele_tsv, &aligned_1, &aligned_2, has_mate)?;

    Ok(())
}

/// `run_bam_fused`'s `--bam-mode no-alignment` route: builds the [`RefKmerFilter`] straight
/// from `-f` (same as the FASTQ path -- no coord-FASTA/gene-interval step) and routes by
/// `alignments`' own `@HD` sort order, exactly like `extract::run_bam_no_alignment`: a
/// coordinate/unsorted file takes the seekable 2-pass name-map
/// ([`bam_extract::extract_from_bam_no_alignment`]); a grouped/name-sorted file takes the
/// one-pass ([`bam_extract::extract_from_bam_no_alignment_grouped`]), fed a factory that always
/// returns a fresh [`InMemoryCandidateSink`] (there is no output-naming decision to make in the
/// fused path, unlike the CLI `extract` subcommand's `FastqFileSink`).
///
/// # Errors
///
/// Returns an error if `-c` was given (rejected: `no-alignment` has no aligned intervals to
/// build gene records from, mirroring `extract::require_no_coord_fasta`), the reference FASTA
/// cannot be opened/parsed, the BAM/CRAM cannot be opened, or the underlying extraction itself
/// fails.
fn run_bam_no_alignment_fused(
    args: &RunArgs,
    bam: &str,
    reference: Option<&str>,
    threads: usize,
) -> anyhow::Result<(bool, InMemoryCandidateSink)> {
    use std::path::Path;
    use unum_core::alignments::{Alignments, SortOrder};
    use unum_core::bam_extract;
    use unum_core::ref_kmer_filter::RefKmerFilter;

    ensure!(
        args.ref_coord_fasta.is_none(),
        "-c (coord FASTA) applies only to run -b --bam-mode alignment, not --bam-mode \
         no-alignment"
    );

    let mut filter = RefKmerFilter::from_reference_fasta(
        Path::new(&args.ref_seq_fasta),
        self::extract::INITIAL_KMER_LENGTH,
    )
    .with_context(|| format!("loading reference FASTA {}", args.ref_seq_fasta))?;
    let mut alignments = Alignments::open_with_reference(bam, reference)
        .with_context(|| format!("opening BAM/CRAM {bam}"))?;

    match alignments.sort_order() {
        SortOrder::QueryName | SortOrder::QueryGrouped => {
            let (metrics, single_end, sink) = bam_extract::extract_from_bam_no_alignment_grouped(
                &mut alignments,
                &mut filter,
                unum_core::extract::DEFAULT_REF_SEQ_SIMILARITY,
                RUN_BAM_MATE_ID_SUFFIX_LEN,
                threads,
                |_single_end| Ok(InMemoryCandidateSink::new()),
            )
            .context("extracting candidate reads from BAM (grouped no-alignment)")?;
            report_bam_extract_metrics(single_end, &metrics);
            Ok((single_end, sink))
        }
        SortOrder::Coordinate | SortOrder::Unsorted => {
            let mut sink = InMemoryCandidateSink::new();
            let metrics = bam_extract::extract_from_bam_no_alignment(
                &mut alignments,
                &mut filter,
                unum_core::extract::DEFAULT_REF_SEQ_SIMILARITY,
                RUN_BAM_MATE_ID_SUFFIX_LEN,
                threads,
                &mut sink,
            )
            .context("extracting candidate reads from BAM (no-alignment)")?;
            report_bam_extract_metrics(metrics.single_end, &metrics);
            Ok((metrics.single_end, sink))
        }
    }
}

/// `run_bam_fused`'s `--bam-mode alignment` route: parses `-c` as a coordinate FASTA (gene
/// intervals + k-mer seed reference, same as `extract::run_coordinate_alignment`), then routes
/// on `alignments`' `@HD` sort order, mirroring `extract::run_bam_alignment`'s dispatch --
/// coordinate-sorted takes [`bam_extract::extract_from_bam_with_threads`]; grouped/name-sorted
/// takes [`bam_extract::extract_from_bam_alignment_grouped`] with `seekable = true` (`run`'s
/// `-b` is always a file/path, never stdin -- there is no `-i -` route into `run`); unsorted is
/// rejected outright, mirroring `extract::run_bam_alignment`.
///
/// # Errors
///
/// Returns an error if `-c` is missing or cannot be opened/parsed, the BAM/CRAM cannot be
/// opened, `alignments` is unsorted, or the underlying extraction itself fails.
fn run_bam_alignment_fused(
    args: &RunArgs,
    bam: &str,
    reference: Option<&str>,
    threads: usize,
) -> anyhow::Result<(bool, InMemoryCandidateSink)> {
    use std::path::Path;
    use unum_core::alignments::{Alignments, SortOrder};
    use unum_core::bam_extract;
    use unum_core::ref_kmer_filter::RefKmerFilter;

    let coord = require_coord_fasta_run(args)?;
    let coord_records = bam_extract::parse_coord_fa(Path::new(coord))
        .with_context(|| format!("parsing coord FASTA {coord}"))?;
    let mut filter =
        RefKmerFilter::from_reference_fasta(Path::new(coord), self::extract::INITIAL_KMER_LENGTH)
            .with_context(|| format!("loading coord FASTA sequences {coord}"))?;
    let mut alignments = Alignments::open_with_reference(bam, reference)
        .with_context(|| format!("opening BAM/CRAM {bam}"))?;
    let genes = bam_extract::build_genes(&alignments, &coord_records)
        .context("resolving coord FASTA chroms to BAM header chrIds")?;

    match alignments.sort_order() {
        SortOrder::Coordinate => {
            let mut sink = InMemoryCandidateSink::new();
            let metrics = bam_extract::extract_from_bam_with_threads(
                &mut alignments,
                &mut filter,
                &genes,
                RUN_BAM_ABNORMAL_UNALIGNED_FLAG,
                RUN_BAM_MATE_ID_SUFFIX_LEN,
                threads,
                &mut sink,
            )
            .context("extracting candidate reads from BAM")?;
            report_bam_extract_metrics(metrics.single_end, &metrics);
            Ok((metrics.single_end, sink))
        }
        SortOrder::QueryName | SortOrder::QueryGrouped => {
            let (metrics, single_end, sink) = bam_extract::extract_from_bam_alignment_grouped(
                &mut alignments,
                &mut filter,
                &genes,
                RUN_BAM_ABNORMAL_UNALIGNED_FLAG,
                RUN_BAM_MATE_ID_SUFFIX_LEN,
                threads,
                true,
                |_single_end| Ok(InMemoryCandidateSink::new()),
            )
            .context("extracting candidate reads from BAM (grouped alignment)")?;
            report_bam_extract_metrics(single_end, &metrics);
            Ok((single_end, sink))
        }
        SortOrder::Unsorted => bail!(
            "run -b --bam-mode alignment requires a coordinate-sorted BAM (@HD SO:coordinate); \
             this input is unsorted/unstated -- run `samtools sort`"
        ),
    }
}

/// Shared BAM-extraction metrics `eprintln` for [`run_bam_fused`]'s two routes, mirroring
/// `extract::report_no_alignment_metrics`'s style/wording for the fused path.
fn report_bam_extract_metrics(single_end: bool, metrics: &BamExtractMetrics) {
    if single_end {
        eprintln!(
            "extracted {} candidate reads (single-end, kmer_length={}, hit_len_required={})",
            metrics.pass1_emitted, metrics.kmer_length, metrics.hit_len_required,
        );
    } else {
        eprintln!(
            "extracted {} + {} = {} candidate pairs (paired, kmer_length={}, \
             hit_len_required={}, candidates_recorded={})",
            metrics.pass1_emitted,
            metrics.pass2_emitted,
            metrics.pass1_emitted + metrics.pass2_emitted,
            metrics.kmer_length,
            metrics.hit_len_required,
            metrics.candidates_recorded,
        );
    }
}

/// [`RunArgs`] twin of `extract::require_reference_for_cram`: resolves the CRAM decoding
/// reference to pass into [`unum_core::alignments::Alignments::open_with_reference`] --
/// `Some(-r)` when `is_cram` is true (erroring if `-r` was not given), else `None`
/// unconditionally (a BAM never gets the reference threaded through, even if `-r` was
/// harmlessly also passed).
///
/// # Errors
///
/// Returns an error if `is_cram` is true and `-r` was not given, naming both CRAM and `-r`
/// plus the `.fai` sibling it requires.
fn require_reference_for_cram_run(args: &RunArgs, is_cram: bool) -> anyhow::Result<Option<&str>> {
    if is_cram {
        Ok(Some(args.reference.as_deref().context(
            "CRAM input requires -r <reference genome FASTA> (with a .fai sibling); the \
             reference is used exclusively -- no REF_PATH/REF_CACHE/network fallback",
        )?))
    } else {
        Ok(None)
    }
}

/// [`RunArgs`] twin of `extract::require_coord_fasta`: extracts the required `-c` gene-
/// coordinate FASTA path for `run -b --bam-mode alignment`, erroring if absent.
fn require_coord_fasta_run(args: &RunArgs) -> anyhow::Result<&str> {
    args.ref_coord_fasta
        .as_deref()
        .context("run -b --bam-mode alignment requires -c <gene coordinate FASTA (*_coord.fa)>")
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
        prefilter_frac: 0.0,
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

/// Combines `-o`/`--prefix` with `--output-dir` into the shared output-file prefix, creating
/// `--output-dir`'s directory if needed — matching `run-t1k`'s
/// `make_path($outputDirectory) ... $prefix = "$outputDirectory/$prefix"`.
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
