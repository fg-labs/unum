//! Thin CLI wrapper around `unum_core::extract` (the Rust port of
//! `fastq-extractor`, `FastqExtractor.cpp`) and
//! `unum_core::bam_extract` (the Rust port of `bam-extractor`,
//! `BamExtractor.cpp`). All extraction logic (data-dependent
//! setup, per-pair/per-pass filtering, `OutputSeq` formatting) lives in
//! `unum-core`; this module only:
//! 1. Routes the input via [`resolve_extract_input`]: `-b` or a
//!    content-sniffed BAM/CRAM under `-i` selects BAM mode (requires
//!    `--bam-mode`); `-1`/`-2`/`-u` (legacy) or FASTQ/FASTA under `-i`
//!    selects FASTQ mode (rejects `--bam-mode`). `-i` is mutually exclusive
//!    with `-1`/`-2`/`-u`/`-b`.
//! 2. FASTQ mode ([`run_fastq_with_source`]): constructs the initial `k=9`
//!    [`unum_core::ref_kmer_filter::RefKmerFilter`] from `-f`, and calls
//!    [`unum_core::extract::extract_candidates_with_threads`] with `-t`
//!    against the paired/single-end/interleaved [`ReadSource`] the router
//!    already resolved.
//! 3. BAM mode ([`run_bam_mode`] → [`run_coordinate_alignment`]): guards the
//!    mode/sort-order/stdin combination (Stage 2a implements only
//!    coordinate-sorted `alignment`; everything else errors as reserved for
//!    a later release), then parses `-f` as a `_coord.fa` (via
//!    [`unum_core::bam_extract::parse_coord_fa`]), builds the
//!    [`RefKmerFilter`] from its sequences and the sorted gene-interval list
//!    (via [`unum_core::bam_extract::build_genes`]), and calls
//!    [`unum_core::bam_extract::extract_from_bam_with_threads`] with `-t`.
//!
//! `-t` controls how many worker threads parallelize the per-read candidate
//! DECISION in both modes; output is byte-identical at any `-t` -- see
//! `unum_core::extract`'s and `unum_core::bam_extract`'s module docs for
//! why.
//!
//! Both modes share [`FastqFileSink`] (`{prefix}_1.fq`/`_2.fq` for paired,
//! `{prefix}.fq` for single-end -- `FastqExtractor.cpp:425-439` /
//! `BamExtractor.cpp:599-610`, identical naming convention).
use crate::cli::{BamMode, ExtractArgs};
use anyhow::{Context, Result, bail, ensure};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use unum_core::alignments::{Alignments, SortOrder};
use unum_core::bam_extract;
use unum_core::extract::{
    self, CandidateSink, ReadRecord, ReadSource, classify_single_input, open_source,
};
use unum_core::fastq::FastqReader;
use unum_core::read_input::{InputSpec, OpenedInput, open_fastq_reader, open_input};
use unum_core::ref_kmer_filter::RefKmerFilter;

/// `FastqExtractor.cpp:272` / `BamExtractor.cpp:480`: the literal initial
/// k-mer length the reference is first loaded at, before any data-dependent
/// `InferKmerLength`/`UpdateKmerLength` adjustment. Shared by both modes
/// (both vendored `main`s use the same literal `9`).
pub(crate) const INITIAL_KMER_LENGTH: usize = 9;

/// Resolves the CLI `-t` value into a worker-thread count for the core
/// extractors: widens `args.threads` into `usize` (saturating to `usize::MAX`
/// on the theoretical narrow-`usize` platform where `u32` doesn't fit), then
/// floors it at `1` (a `threads <= 1` value selects the sequential path in both
/// `unum_core::extract` and `unum_core::bam_extract`). Shared by
/// [`run_fastq_with_source`] and [`run_coordinate_alignment`] so the clamping
/// rule lives in one place.
fn resolve_threads(threads: u32) -> usize {
    usize::try_from(threads).unwrap_or(usize::MAX).max(1)
}

/// Runs the `extract` subcommand for `args`: routes the input via
/// [`resolve_extract_input`] (FASTQ vs. BAM/CRAM, enforcing the
/// `--bam-mode`-required-for-BAM / `--bam-mode`-forbidden-for-FASTQ rules),
/// then dispatches to [`run_fastq_with_source`] or [`run_bam_mode`].
///
/// # Errors
///
/// Returns an error if [`resolve_extract_input`] rejects the input/flag
/// combination, or the resolved runner ([`run_fastq_with_source`]/
/// [`run_bam_mode`]) itself fails.
pub fn run(args: &ExtractArgs) -> Result<()> {
    match resolve_extract_input(args)? {
        ResolvedExtractInput::Fastq(source) => run_fastq_with_source(args, source),
        ResolvedExtractInput::Bam { spec, mode } => run_bam_mode(args, &spec, mode),
    }
}

/// The reserved-for-a-later-release error for CRAM input. CRAM needs the
/// reference wiring (`-r`) that Stage 2c introduces; until then a detected
/// CRAM is rejected at routing time (in both the `-b` and `-i` arms) with
/// this message so a `BamInputSpec` only ever denotes a real BAM.
const CRAM_RESERVED_MSG: &str = "CRAM input requires reference wiring (-r) that is not yet available (later release); \
     convert to BAM (samtools view -b) or use a BAM input";

/// Where a BAM input comes from, once routing has decided it IS a BAM
/// (CRAM is rejected earlier, at routing time -- see [`CRAM_RESERVED_MSG`]).
enum BamInputSpec {
    Path(String),
    Stdin,
}

/// The resolved extraction input: FASTQ (a ready [`ReadSource`]) or BAM/CRAM
/// (a path/stdin spec plus the chosen [`BamMode`]).
enum ResolvedExtractInput {
    Fastq(ReadSource),
    Bam { spec: BamInputSpec, mode: BamMode },
}

/// Routes `-i`/`-b`/legacy flags to FASTQ or BAM, enforcing the
/// `--bam-mode`-required-for-BAM and `--bam-mode`-forbidden-for-FASTQ rules.
/// CRAM is content-sniffed and rejected as reserved for a later release (it
/// needs the `-r` reference wiring Stage 2c adds) in both the `-b` and `-i`
/// arms, so a resolved BAM is always a real BAM.
///
/// - `-b <path>` → BAM (requires `--bam-mode`); content-sniffed to reject
///   CRAM (reserved) and a mistaken FASTQ/FASTA `-b`.
/// - `-i` with 2 paths → paired FASTQ (rejects `--bam-mode`).
/// - `-i` with 1 path: `-` (stdin) routes by `--bam-mode` presence
///   (present → BAM stdin; absent → FASTQ stdin); a file is content-sniffed
///   via [`open_input`] (BAM → BAM, requires `--bam-mode`; CRAM → reserved
///   error; FASTQ/FASTA → FASTQ, rejects `--bam-mode`).
/// - legacy `-1`/`-2`/`-u` → FASTQ (rejects `--bam-mode`).
///
/// # Errors
///
/// Returns an error on a `--bam-mode`/input mismatch, a missing `--bam-mode`
/// for BAM, a detected CRAM (reserved), or any underlying open/detection
/// failure.
fn resolve_extract_input(args: &ExtractArgs) -> Result<ResolvedExtractInput> {
    // -b takes the BAM path directly.
    if let Some(bam_path) = args.bam.as_deref() {
        ensure!(
            args.input.is_empty()
                && args.mate1.is_none()
                && args.mate2.is_none()
                && args.single.is_none(),
            "-b (BAM mode) is mutually exclusive with -i and -1/-2/-u"
        );
        let mode = require_bam_mode(args)?;
        // Content-sniff (cheap magic-byte peek) to reject CRAM (reserved) and
        // a mistaken FASTQ/FASTA `-b` before `Alignments::open` is reached; a
        // `BamInputSpec::Path` must only ever denote a real BAM.
        let (opened, _fmt) = open_input(&InputSpec::Path(std::path::PathBuf::from(bam_path)))
            .with_context(|| format!("opening input {bam_path}"))?;
        match opened {
            OpenedInput::Bam => {}
            OpenedInput::Cram => bail!("{CRAM_RESERVED_MSG}"),
            OpenedInput::Fastq(_) => {
                bail!("-b expects a BAM/CRAM file, but {bam_path} is FASTQ/FASTA input")
            }
        }
        return Ok(ResolvedExtractInput::Bam {
            spec: BamInputSpec::Path(bam_path.to_string()),
            mode,
        });
    }

    // -i unified input.
    if !args.input.is_empty() {
        ensure!(
            args.mate1.is_none() && args.mate2.is_none() && args.single.is_none(),
            "-i/--input is mutually exclusive with -1/-2/-u/-b"
        );
        return resolve_i_input(args);
    }

    // Legacy -1/-2/-u FASTQ -- never BAM.
    ensure!(
        args.bam_mode.is_none(),
        "--bam-mode applies to BAM/CRAM input only, not -1/-2/-u FASTQ"
    );
    let source = resolve_legacy_fastq_source(args)?;
    Ok(ResolvedExtractInput::Fastq(source))
}

/// Resolves the 1-2 `-i` inputs.
fn resolve_i_input(args: &ExtractArgs) -> Result<ResolvedExtractInput> {
    match args.input.as_slice() {
        [a, b] => {
            // Two inputs = paired FASTQ; BAM is a single input.
            ensure!(
                args.bam_mode.is_none(),
                "--bam-mode applies to BAM/CRAM input, not paired FASTQ"
            );
            ensure!(
                !(a == "-" && b == "-"),
                "both paired inputs are stdin ('-'); use interleaved '-i -'"
            );
            let r1 = open_one_fastq(a)?;
            let r2 = open_one_fastq(b)?;
            Ok(ResolvedExtractInput::Fastq(ReadSource::Paired(r1, r2)))
        }
        [a] if a == "-" => {
            // stdin: route by --bam-mode presence (can't content-sniff a pipe for htslib).
            if let Some(mode) = args.bam_mode {
                Ok(ResolvedExtractInput::Bam { spec: BamInputSpec::Stdin, mode })
            } else {
                let reader = open_one_fastq("-")?;
                Ok(ResolvedExtractInput::Fastq(classify_single_input(reader)?))
            }
        }
        [a] => {
            // File: content-sniff.
            let (opened, _fmt) = open_input(&InputSpec::Path(std::path::PathBuf::from(a)))
                .with_context(|| format!("opening input {a}"))?;
            match opened {
                OpenedInput::Fastq(reader) => {
                    ensure!(
                        args.bam_mode.is_none(),
                        "--bam-mode given but {a} is FASTQ/FASTA input"
                    );
                    Ok(ResolvedExtractInput::Fastq(classify_single_input(reader)?))
                }
                OpenedInput::Bam => {
                    let mode = require_bam_mode(args)?;
                    Ok(ResolvedExtractInput::Bam { spec: BamInputSpec::Path(a.clone()), mode })
                }
                OpenedInput::Cram => bail!("{CRAM_RESERVED_MSG}"),
            }
        }
        _ => bail!("-i/--input takes 1 or 2 paths, got {}", args.input.len()),
    }
}

/// Extracts the required `--bam-mode`, erroring if absent.
fn require_bam_mode(args: &ExtractArgs) -> Result<BamMode> {
    args.bam_mode.context(
        "BAM/CRAM input requires --bam-mode {alignment|no-alignment}; \
         e.g. --bam-mode alignment (coordinate-sorted bam-extractor parity)",
    )
}

/// Opens one FASTQ/FASTA input (path or `-` for stdin) with content-based
/// format detection and niffler decompression, erroring on BAM/CRAM.
fn open_one_fastq(spec: &str) -> Result<FastqReader> {
    let input = if spec == "-" {
        InputSpec::Stdin
    } else {
        InputSpec::Path(std::path::PathBuf::from(spec))
    };
    let (reader, _fmt) =
        open_fastq_reader(&input).with_context(|| format!("opening input {spec}"))?;
    Ok(reader)
}

/// Resolves the legacy `-1`/`-2`/`-u` flags into a [`ReadSource`], mirroring
/// the pre-router `run_fastq`'s source-building exactly (same ensure-messages,
/// same [`open_source`] call) so the legacy FASTQ path stays byte-identical.
///
/// # Errors
///
/// Returns an error if: neither `-u` nor both `-1`/`-2` are given (or both
/// single-end and paired flags are given together); or the read file(s)
/// cannot be opened.
fn resolve_legacy_fastq_source(args: &ExtractArgs) -> Result<ReadSource> {
    let mate2 = args.mate2.as_deref();
    let paired = args.mate1.is_some() || mate2.is_some();
    let single = args.single.as_deref();

    ensure!(
        !(paired && single.is_some()),
        "specify either -u (single-end) or -1/-2 (paired), not both"
    );
    let (mate1_path, mate2_path): (&str, Option<&str>) = if paired {
        let mate1 = args
            .mate1
            .as_deref()
            .context("paired input requires both -1 and -2 (got -2 without -1)")?;
        let mate2 = mate2.context("paired input requires both -1 and -2 (got -1 without -2)")?;
        (mate1, Some(mate2))
    } else if let Some(single) = single {
        (single, None)
    } else {
        bail!("must specify either -u (single-end) or -1/-2 (paired) read input, or -b (BAM)");
    };

    open_source(Path::new(mate1_path), mate2_path.map(Path::new)).context("opening read source")
}

/// Runs FASTQ-mode extraction against an already-resolved [`ReadSource`] --
/// the shared executor for both `-i` FASTQ and legacy `-1`/`-2`/`-u` FASTQ
/// (the pre-existing `fastq-extractor` port). `paired` (for output-file
/// naming and the metrics log's "pairs"/"reads" label) is derived from
/// `source`'s variant, matching both predecessor paths' derivations exactly
/// ([`ReadSource::Paired`]/[`ReadSource::Interleaved`] ⇒ paired).
///
/// # Errors
///
/// Returns an error if the reference FASTA cannot be opened/parsed, or
/// [`unum_core::extract::extract_candidates`] itself fails (e.g. an empty
/// read-1 file or mismatched mate-pair counts -- see that function's doc
/// comment).
fn run_fastq_with_source(args: &ExtractArgs, mut source: ReadSource) -> Result<()> {
    let mut filter =
        RefKmerFilter::from_reference_fasta(Path::new(&args.ref_seq_fasta), INITIAL_KMER_LENGTH)
            .with_context(|| format!("loading reference FASTA {}", args.ref_seq_fasta))?;

    let paired = matches!(source, ReadSource::Paired(_, _) | ReadSource::Interleaved(_));

    let mut sink =
        FastqFileSink::create(&args.prefix, paired).context("creating output FASTQ file(s)")?;

    let threads = resolve_threads(args.threads);
    let metrics = extract::extract_candidates_with_threads(
        &mut source,
        &mut filter,
        args.similarity,
        threads,
        &mut sink,
    )
    .context("extracting candidate reads")?;

    eprintln!(
        "extracted {} / {} candidate {} (kmer_length={}, hit_len_required={})",
        metrics.candidates_emitted,
        metrics.total_reads,
        if paired { "pairs" } else { "reads" },
        metrics.kmer_length,
        metrics.hit_len_required,
    );

    sink.flush()?;
    Ok(())
}

/// Executes BAM extraction for the resolved `mode` (CRAM is already rejected
/// at routing time, so `spec` only ever denotes a BAM). Stage 2a implements
/// only coordinate-sorted `alignment`; everything else is a loud
/// reserved-for-a-later-release error.
///
/// # Errors
///
/// Returns an error for `no-alignment` (Stage 2b), grouped/name-sorted or
/// unsorted input, or stdin (Stage 2c), or any extraction failure.
fn run_bam_mode(args: &ExtractArgs, spec: &BamInputSpec, mode: BamMode) -> Result<()> {
    ensure!(
        mode == BamMode::Alignment,
        "--bam-mode no-alignment is not yet available (BAM-as-reads / Class-A \
         selection ships in a later release); use --bam-mode alignment"
    );

    // Stage 2a alignment == the existing coordinate 2-pass, which needs a
    // seekable file. stdin is reserved here (CRAM was already rejected at
    // routing time, so `spec` is always a BAM).
    let bam_path = match spec {
        BamInputSpec::Path(p) => p.as_str(),
        BamInputSpec::Stdin => bail!(
            "coordinate-sorted BAM extraction needs a seekable file; a BAM from stdin \
             requires the grouped/name-sorted one-pass (available in a later release) -- \
             redirect to a file, or pass a grouped BAM once supported"
        ),
    };

    let alignments =
        Alignments::open(bam_path).with_context(|| format!("opening BAM/CRAM {bam_path}"))?;

    // Sort-order guard: alignment (coordinate 2-pass) requires SO:coordinate.
    match alignments.sort_order() {
        SortOrder::Coordinate => {}
        SortOrder::QueryName | SortOrder::QueryGrouped => bail!(
            "--bam-mode alignment on a grouped/name-sorted BAM uses a one-pass interval \
             matcher that is not yet available (later release); coordinate-sort the BAM \
             (samtools sort) to use the current alignment path"
        ),
        SortOrder::Unsorted => bail!(
            "--bam-mode alignment requires a coordinate-sorted BAM (@HD SO:coordinate); \
             this input is unsorted/unstated -- run `samtools sort`"
        ),
    }

    run_coordinate_alignment(args, alignments)
}

/// The coordinate-sorted `alignment` extraction -- the former `run_bam` body,
/// now reached only after the mode/sort/stdin/CRAM guards pass. Byte-identical
/// to Stage 1's `-b` output: parses `-f` as a `_coord.fa`, builds the
/// [`RefKmerFilter`] and sorted gene-interval list, and calls
/// [`unum_core::bam_extract::extract_from_bam`]. Output file naming
/// (single-end vs. paired) is decided by the BAM's own sampled
/// `frag_stdev` (`BamExtractor.cpp:599-610`), not a CLI flag -- so this
/// function samples `general_info` itself up front (before opening
/// [`FastqFileSink`]) purely to pick the right filename(s); see this
/// function's body for why re-sampling here (rather than letting
/// `extract_from_bam`'s own internal sample decide) is harmless.
///
/// # Errors
///
/// Returns an error if the coord FASTA cannot be opened/parsed, or if
/// [`unum_core::bam_extract::extract_from_bam`] itself fails (e.g. an
/// unaligned-template mate-pairing error -- see that function's doc comment).
fn run_coordinate_alignment(args: &ExtractArgs, mut alignments: Alignments) -> Result<()> {
    let coord_records = bam_extract::parse_coord_fa(Path::new(&args.ref_seq_fasta))
        .with_context(|| format!("parsing coord FASTA {}", args.ref_seq_fasta))?;

    let mut filter =
        RefKmerFilter::from_reference_fasta(Path::new(&args.ref_seq_fasta), INITIAL_KMER_LENGTH)
            .with_context(|| format!("loading coord FASTA sequences {}", args.ref_seq_fasta))?;

    let genes = bam_extract::build_genes(&alignments, &coord_records)
        .context("resolving coord FASTA chroms to BAM header chrIds")?;

    // Output naming (`{prefix}.fq` vs. `{prefix}_1.fq`/`_2.fq`) depends on
    // frag_stdev (single-end vs. paired), matching BamExtractor.cpp:573-610:
    // `GetGeneralInfo` is called BEFORE the output files are opened there.
    // `extract_from_bam` computes this same statistic again internally (it
    // owns the full data-dependent setup sequence, including the
    // `hitLenRequired`/`InferKmerLength` steps that also depend on it) --
    // calling `general_info` a second time here is a cheap, harmless extra
    // BAM pass purely to decide a filename up front, not a semantic
    // divergence (both calls sample the same file from the same rewound
    // position and must agree).
    let single_end =
        alignments.general_info(true).context("sampling BAM for output naming")?.frag_stdev == 0;
    alignments.rewind().context("rewinding after output-naming sample")?;

    let mut sink = FastqFileSink::create(&args.prefix, !single_end)
        .with_context(|| format!("creating output FASTQ file(s) for prefix {}", args.prefix))?;

    let threads = resolve_threads(args.threads);
    let metrics = bam_extract::extract_from_bam_with_threads(
        &mut alignments,
        &mut filter,
        &genes,
        args.abnormal_unmapped,
        args.mate_id_suffix_len,
        threads,
        &mut sink,
    )
    .context("extracting candidate reads from BAM")?;
    sink.flush()?;

    if metrics.single_end {
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

    Ok(())
}

/// Writes candidate pairs/reads to FASTQ file(s), ported from
/// `FastqExtractor.cpp:425-439`'s output-file naming (`{prefix}_1.fq` /
/// `{prefix}_2.fq` for paired, `{prefix}.fq` for single-end) plus
/// `OutputSeq` (`FastqExtractor.cpp:120-153`, via
/// [`unum_core::extract::output_seq`]) for the actual record formatting.
///
/// # Both output mates use mate 1's id
///
/// `FastqExtractor.cpp:471-473` calls `OutputSeq(fp1, reads.id, reads.seq,
/// reads.qual, ...)` for mate 1 AND `OutputSeq(fp2, reads.id, mateReads.seq,
/// mateReads.qual, ...)` for mate 2 -- note `reads.id` (mate 1's id) is used
/// for BOTH calls; mate 2's own kseq-parsed id is discarded entirely on the
/// output path (even though it was read and is available). [`emit_pair`]
/// reproduces this exactly.
struct FastqFileSink {
    fp1: BufWriter<File>,
    fp2: Option<BufWriter<File>>,
    read1_start: i64,
    read1_end: i64,
    read2_start: i64,
    read2_end: i64,
}

impl FastqFileSink {
    /// Opens `{prefix}_1.fq` + `{prefix}_2.fq` (paired) or `{prefix}.fq`
    /// (single-end), matching `FastqExtractor.cpp:425-439` exactly. Trim
    /// parameters default to "no trim" (`start=0, end=-1`, matching
    /// `FastqExtractor.cpp`'s own `read1Start`/`read1End`/`read2Start`/
    /// `read2End` defaults, `FastqExtractor.cpp:286-289`) -- this CLI does
    /// not currently expose `--read1Start`/etc. (kept for a future
    /// generality pass; [`unum_core::extract::output_seq`] already
    /// supports them).
    fn create(prefix: &str, paired: bool) -> Result<Self> {
        if paired {
            let fp1 = File::create(format!("{prefix}_1.fq"))
                .with_context(|| format!("creating {prefix}_1.fq"))?;
            let fp2 = File::create(format!("{prefix}_2.fq"))
                .with_context(|| format!("creating {prefix}_2.fq"))?;
            Ok(Self {
                fp1: BufWriter::new(fp1),
                fp2: Some(BufWriter::new(fp2)),
                read1_start: 0,
                read1_end: -1,
                read2_start: 0,
                read2_end: -1,
            })
        } else {
            let fp1 = File::create(format!("{prefix}.fq"))
                .with_context(|| format!("creating {prefix}.fq"))?;
            Ok(Self {
                fp1: BufWriter::new(fp1),
                fp2: None,
                read1_start: 0,
                read1_end: -1,
                read2_start: 0,
                read2_end: -1,
            })
        }
    }

    fn flush(&mut self) -> Result<()> {
        self.fp1.flush().context("flushing mate-1 output file")?;
        if let Some(fp2) = &mut self.fp2 {
            fp2.flush().context("flushing mate-2 output file")?;
        }
        Ok(())
    }
}

impl CandidateSink for FastqFileSink {
    fn emit_pair(&mut self, r1: &ReadRecord, r2: Option<&ReadRecord>) -> Result<()> {
        extract::output_seq(
            &mut self.fp1,
            &r1.id,
            &r1.seq,
            r1.qual.as_deref(),
            self.read1_start,
            self.read1_end,
        )
        .context("writing mate-1 candidate record")?;

        if let (Some(r2), Some(fp2)) = (r2, self.fp2.as_mut()) {
            // r1.id (NOT r2.id) is used here -- see this struct's doc
            // comment.
            extract::output_seq(
                fp2,
                &r1.id,
                &r2.seq,
                r2.qual.as_deref(),
                self.read2_start,
                self.read2_end,
            )
            .context("writing mate-2 candidate record")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::BamMode;

    /// Builds an `ExtractArgs` with every field defaulted to a FASTQ-mode
    /// baseline (no `-b`, no `-i`, no `--bam-mode`); individual tests
    /// override just the fields they need.
    fn args() -> ExtractArgs {
        ExtractArgs {
            ref_seq_fasta: "ref.fa".into(),
            mate1: None,
            mate2: None,
            single: None,
            input: vec![],
            bam: None,
            abnormal_unmapped: false,
            mate_id_suffix_len: -1,
            prefix: "out".into(),
            threads: 1,
            similarity: unum_core::extract::DEFAULT_REF_SEQ_SIMILARITY,
            bam_mode: None,
        }
    }

    /// `resolve_extract_input`'s `Ok` variant contains a [`ReadSource`], which
    /// (via its `FastqReader`/`Box<dyn BufRead>`) does not implement `Debug`,
    /// so `Result::unwrap_err` can't be used directly. This extracts just
    /// the error message for assertions.
    fn resolve_extract_input_err_message(args: &ExtractArgs) -> String {
        match resolve_extract_input(args) {
            Ok(_) => panic!("expected resolve_extract_input to fail"),
            Err(e) => e.to_string(),
        }
    }

    #[test]
    fn bam_input_requires_bam_mode() {
        let mut a = args();
        a.bam = Some("x.bam".into()); // -b path, no --bam-mode
        let err = resolve_extract_input_err_message(&a);
        assert!(err.contains("--bam-mode"), "must demand --bam-mode: {err}");
    }

    #[test]
    fn b_flag_with_alignment_resolves_to_bam() {
        // `-b` now content-sniffs, so this needs a real BAM-magic file (the
        // path is opened for a magic-byte peek before it resolves).
        let (_dir, path) = write_temp_input(b"BAM\x01rest-of-a-fake-bam", "x.bam");
        let mut a = args();
        a.bam = Some(path);
        a.bam_mode = Some(BamMode::Alignment);
        let resolved = resolve_extract_input(&a).unwrap();
        assert!(matches!(resolved, ResolvedExtractInput::Bam { mode: BamMode::Alignment, .. }));
    }

    #[test]
    fn bam_mode_with_fastq_input_is_rejected() {
        let mut a = args();
        a.single = Some("reads.fq".into());
        a.bam_mode = Some(BamMode::Alignment); // --bam-mode on a FASTQ input
        let err = resolve_extract_input_err_message(&a);
        assert!(err.contains("--bam-mode"), "must reject --bam-mode on FASTQ: {err}");
    }

    #[test]
    fn no_alignment_mode_is_reserved() {
        let spec = BamInputSpec::Path("x.bam".into());
        let err = run_bam_mode(&args(), &spec, BamMode::NoAlignment).unwrap_err().to_string();
        assert!(
            err.contains("no-alignment") && err.contains("later release"),
            "no-alignment must be reserved: {err}"
        );
    }

    #[test]
    fn bam_from_stdin_is_reserved() {
        let err = run_bam_mode(&args(), &BamInputSpec::Stdin, BamMode::Alignment)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("stdin") || err.contains("seekable"),
            "BAM from stdin must be reserved with a file-redirect hint: {err}"
        );
    }

    /// Writes `bytes` to a fresh file named `name` under a tempdir, returning
    /// `(dir, path_string)`. The dir must outlive the path (it owns the file).
    fn write_temp_input(bytes: &[u8], name: &str) -> (tempfile::TempDir, String) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join(name);
        std::fs::write(&path, bytes).unwrap();
        let path_str = path.to_str().unwrap().to_string();
        (dir, path_str)
    }

    #[test]
    fn i_flag_cram_file_is_reserved() {
        // `CRAM\x03` + padding (>= 12 bytes so niffler doesn't FileTooShort).
        let (_dir, path) = write_temp_input(b"CRAM\x03padding-bytes", "x.cram");
        let mut a = args();
        a.input = vec![path];
        a.bam_mode = Some(BamMode::Alignment);
        let err = resolve_extract_input_err_message(&a);
        assert!(
            err.contains("CRAM") && err.contains("later release"),
            "CRAM via -i must be reserved: {err}"
        );
    }

    #[test]
    fn b_flag_cram_file_is_reserved() {
        let (_dir, path) = write_temp_input(b"CRAM\x03padding-bytes", "x.cram");
        let mut a = args();
        a.bam = Some(path);
        a.bam_mode = Some(BamMode::Alignment);
        let err = resolve_extract_input_err_message(&a);
        assert!(
            err.contains("CRAM") && err.contains("later release"),
            "CRAM via -b must be reserved: {err}"
        );
    }

    #[test]
    fn i_flag_bam_file_resolves_to_bam() {
        let (_dir, path) = write_temp_input(b"BAM\x01rest-of-a-fake-bam", "x.bam");
        let mut a = args();
        a.input = vec![path.clone()];
        a.bam_mode = Some(BamMode::Alignment);
        let resolved = resolve_extract_input(&a).unwrap();
        match resolved {
            ResolvedExtractInput::Bam { spec: BamInputSpec::Path(p), mode: BamMode::Alignment } => {
                assert_eq!(p, path);
            }
            _ => panic!("a BAM-magic -i file must resolve to a BAM path"),
        }
    }

    #[test]
    fn fastq_file_sink_uses_mate1_id_for_both_outputs() {
        let tmp = tempfile::tempdir().unwrap();
        let prefix = tmp.path().join("out").to_str().unwrap().to_string();
        let mut sink = FastqFileSink::create(&prefix, true).unwrap();

        let r1 = ReadRecord {
            id: "mate1_id".to_string(),
            seq: b"ACGT".to_vec(),
            qual: Some(b"IIII".to_vec()),
        };
        let r2 = ReadRecord {
            id: "totally_different_mate2_id".to_string(),
            seq: b"TTTT".to_vec(),
            qual: Some(b"JJJJ".to_vec()),
        };
        sink.emit_pair(&r1, Some(&r2)).unwrap();
        sink.flush().unwrap();

        let out1 = std::fs::read_to_string(format!("{prefix}_1.fq")).unwrap();
        let out2 = std::fs::read_to_string(format!("{prefix}_2.fq")).unwrap();
        assert_eq!(out1, "@mate1_id\nACGT\n+\nIIII\n");
        assert_eq!(out2, "@mate1_id\nTTTT\n+\nJJJJ\n", "mate-2 output must use mate-1's id");
    }

    #[test]
    fn fastq_file_sink_single_end_writes_one_file() {
        let tmp = tempfile::tempdir().unwrap();
        let prefix = tmp.path().join("out").to_str().unwrap().to_string();
        let mut sink = FastqFileSink::create(&prefix, false).unwrap();

        let r1 = ReadRecord {
            id: "s0".to_string(),
            seq: b"ACGT".to_vec(),
            qual: Some(b"IIII".to_vec()),
        };
        sink.emit_pair(&r1, None).unwrap();
        sink.flush().unwrap();

        assert!(!std::path::Path::new(&format!("{prefix}_1.fq")).exists());
        let out = std::fs::read_to_string(format!("{prefix}.fq")).unwrap();
        assert_eq!(out, "@s0\nACGT\n+\nIIII\n");
    }
}
