//! Thin CLI wrapper around `fg_t1k_core::extract`, the Rust port of
//! `fastq-extractor` (`vendor/t1k/FastqExtractor.cpp`). All extraction
//! logic (data-dependent setup, per-pair filtering, `OutputSeq` formatting)
//! lives in `fg-t1k-core`; this module only:
//! 1. Constructs the initial `k=9` [`fg_t1k_core::ref_kmer_filter::RefKmerFilter`]
//!    from `-f`.
//! 2. Constructs the paired/single-end read source from `-1`/`-2` or `-u`.
//! 3. Opens the output file(s) (`{prefix}_1.fq`/`_2.fq` for paired,
//!    `{prefix}.fq` for single-end -- `FastqExtractor.cpp:425-439`) and
//!    wraps them in a [`FastqFileSink`].
//! 4. Calls [`fg_t1k_core::extract::extract_candidates`].
//!
//! # Follow-up: wiring into the `--engine` strangler router
//!
//! `stages::run`'s `extract` stage currently only dispatches to the C++
//! oracle (`Engine::Cpp` in `stages::run`'s `match overrides.engine_for(...)`)
//! and returns an error for `Engine::Rust`. Wiring `extract_candidates` in as
//! that `Engine::Rust` implementation is a deliberate follow-up, not done in
//! this task -- this module currently only exposes a standalone `fg-t1k
//! extract` subcommand.
use crate::cli::ExtractArgs;
use anyhow::{Context, Result, bail, ensure};
use fg_t1k_core::extract::{self, CandidateSink, ReadRecord};
use fg_t1k_core::ref_kmer_filter::RefKmerFilter;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// `FastqExtractor.cpp:272`: the literal initial k-mer length the reference
/// is first loaded at, before any data-dependent `InferKmerLength`/
/// `UpdateKmerLength` adjustment.
const INITIAL_KMER_LENGTH: usize = 9;

/// Runs the `extract` subcommand for `args`.
///
/// # Errors
///
/// Returns an error if: neither `-u` nor both `-1`/`-2` are given (or both
/// single-end and paired flags are given together); the reference or read
/// files cannot be opened/parsed; or [`fg_t1k_core::extract::extract_candidates`]
/// itself fails (e.g. an empty read-1 file or mismatched mate-pair counts --
/// see that function's doc comment).
pub fn run(args: &ExtractArgs) -> Result<()> {
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
        bail!("must specify either -u (single-end) or -1/-2 (paired) read input");
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

/// Writes candidate pairs/reads to FASTQ file(s), ported from
/// `FastqExtractor.cpp:425-439`'s output-file naming (`{prefix}_1.fq` /
/// `{prefix}_2.fq` for paired, `{prefix}.fq` for single-end) plus
/// `OutputSeq` (`FastqExtractor.cpp:120-153`, via
/// [`fg_t1k_core::extract::output_seq`]) for the actual record formatting.
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
    /// generality pass; [`fg_t1k_core::extract::output_seq`] already
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
