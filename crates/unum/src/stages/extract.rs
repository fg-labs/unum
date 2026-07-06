//! Thin CLI wrapper around `unum_core::extract` (the Rust port of
//! `fastq-extractor`, `vendor/t1k/FastqExtractor.cpp`) and
//! `unum_core::bam_extract` (the Rust port of `bam-extractor`,
//! `vendor/t1k/BamExtractor.cpp`). All extraction logic (data-dependent
//! setup, per-pair/per-pass filtering, `OutputSeq` formatting) lives in
//! `unum-core`; this module only:
//! 1. Dispatches on `-b`: BAM mode ([`run_bam`]) if given, else FASTQ mode
//!    ([`run_fastq`]).
//! 2. FASTQ mode: constructs the initial `k=9`
//!    [`unum_core::ref_kmer_filter::RefKmerFilter`] from `-f`, the
//!    paired/single-end read source from `-1`/`-2` or `-u`, and calls
//!    [`unum_core::extract::extract_candidates`].
//! 3. BAM mode: parses `-f` as a `_coord.fa` (via
//!    [`unum_core::bam_extract::parse_coord_fa`]), builds the
//!    [`RefKmerFilter`] from its sequences and the sorted gene-interval list
//!    (via [`unum_core::bam_extract::build_genes`]), opens `-b` as an
//!    [`unum_core::alignments::Alignments`], and calls
//!    [`unum_core::bam_extract::extract_from_bam`].
//!
//! Both modes share [`FastqFileSink`] (`{prefix}_1.fq`/`_2.fq` for paired,
//! `{prefix}.fq` for single-end -- `FastqExtractor.cpp:425-439` /
//! `BamExtractor.cpp:599-610`, identical naming convention).
use crate::cli::ExtractArgs;
use anyhow::{Context, Result, bail, ensure};
use unum_core::alignments::Alignments;
use unum_core::bam_extract;
use unum_core::extract::{self, CandidateSink, ReadRecord};
use unum_core::ref_kmer_filter::RefKmerFilter;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// `FastqExtractor.cpp:272` / `BamExtractor.cpp:480`: the literal initial
/// k-mer length the reference is first loaded at, before any data-dependent
/// `InferKmerLength`/`UpdateKmerLength` adjustment. Shared by both modes
/// (both vendored `main`s use the same literal `9`).
const INITIAL_KMER_LENGTH: usize = 9;

/// Runs the `extract` subcommand for `args`: dispatches to [`run_bam`] if
/// `-b` was given, else [`run_fastq`].
///
/// # Errors
///
/// See [`run_bam`]/[`run_fastq`].
pub fn run(args: &ExtractArgs) -> Result<()> {
    if let Some(bam_path) = args.bam.as_deref() {
        ensure!(
            args.mate1.is_none() && args.mate2.is_none() && args.single.is_none(),
            "-b (BAM mode) is mutually exclusive with -1/-2/-u (FASTQ mode)"
        );
        return run_bam(args, bam_path);
    }
    run_fastq(args)
}

/// Runs FASTQ-mode extraction (the pre-existing `fastq-extractor` port).
///
/// # Errors
///
/// Returns an error if: neither `-u` nor both `-1`/`-2` are given (or both
/// single-end and paired flags are given together); the reference or read
/// files cannot be opened/parsed; or [`unum_core::extract::extract_candidates`]
/// itself fails (e.g. an empty read-1 file or mismatched mate-pair counts --
/// see that function's doc comment).
fn run_fastq(args: &ExtractArgs) -> Result<()> {
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

    let mut filter =
        RefKmerFilter::from_reference_fasta(Path::new(&args.ref_seq_fasta), INITIAL_KMER_LENGTH)
            .with_context(|| format!("loading reference FASTA {}", args.ref_seq_fasta))?;

    let mut source = extract::open_source(Path::new(mate1_path), mate2_path.map(Path::new))
        .context("opening read source")?;

    let mut sink = FastqFileSink::create(&args.prefix, mate2_path.is_some())
        .context("creating output FASTQ file(s)")?;

    let metrics = extract::extract_candidates(&mut source, &mut filter, args.similarity, &mut sink)
        .context("extracting candidate reads")?;

    eprintln!(
        "extracted {} / {} candidate {} (kmer_length={}, hit_len_required={})",
        metrics.candidates_emitted,
        metrics.total_reads,
        if mate2_path.is_some() { "pairs" } else { "reads" },
        metrics.kmer_length,
        metrics.hit_len_required,
    );

    sink.flush()?;
    Ok(())
}

/// Runs BAM-mode extraction (the `bam-extractor` port): parses `-f` as a
/// `_coord.fa`, opens `bam_path`, builds the [`RefKmerFilter`] and sorted
/// gene-interval list, and calls
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
/// Returns an error if the coord FASTA or BAM/CRAM file cannot be
/// opened/parsed, or if [`unum_core::bam_extract::extract_from_bam`] itself
/// fails (e.g. an unaligned-template mate-pairing error -- see that
/// function's doc comment).
fn run_bam(args: &ExtractArgs, bam_path: &str) -> Result<()> {
    let coord_records = bam_extract::parse_coord_fa(Path::new(&args.ref_seq_fasta))
        .with_context(|| format!("parsing coord FASTA {}", args.ref_seq_fasta))?;

    let mut filter =
        RefKmerFilter::from_reference_fasta(Path::new(&args.ref_seq_fasta), INITIAL_KMER_LENGTH)
            .with_context(|| format!("loading coord FASTA sequences {}", args.ref_seq_fasta))?;

    let mut alignments =
        Alignments::open(bam_path).with_context(|| format!("opening BAM/CRAM {bam_path}"))?;

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

    let metrics = bam_extract::extract_from_bam(
        &mut alignments,
        &mut filter,
        &genes,
        args.abnormal_unmapped,
        args.mate_id_suffix_len,
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
