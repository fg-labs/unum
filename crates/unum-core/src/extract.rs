//! Read-candidate extraction driver, ported from `FastqExtractor.cpp`'s
//! `main` (`FastqExtractor.cpp:259-628`): drives a configured
//! [`crate::ref_kmer_filter::RefKmerFilter`] over paired or single-end FASTQ
//! input and emits candidate reads (pairs where at least one mate is a good
//! candidate) to a [`CandidateSink`], in input order.
//!
//! # Library-first design
//!
//! This module is deliberately I/O-sink-agnostic: [`extract_candidates`]
//! takes a [`CandidateSink`] trait object rather than writing files itself,
//! so a future fused `genotype` command can reuse the exact same
//! candidate-selection logic while handing passing records directly to the
//! next pipeline stage in memory, without round-tripping through
//! intermediate FASTQ files. The CLI's `extract` subcommand (`unum`
//! binary crate) is a thin wrapper that constructs a FASTQ-file-writing
//! [`CandidateSink`] and calls this function -- see that crate's `extract`
//! stage module for the file-writing sink.
//!
//! # Data-dependent setup: `hitLenRequired` sampling, `InferKmerLength`,
//! `UpdateKmerLength`
//!
//! `FastqExtractor.cpp:main` does NOT use `RefKmerFilter`'s constructor
//! defaults; it performs a specific, data-dependent setup sequence BEFORE
//! any read is evaluated (`FastqExtractor.cpp:271-418`), reproduced exactly
//! by [`extract_candidates`]:
//!
//! 1. `kmerLength = 9` (a literal initial value, `FastqExtractor.cpp:272`);
//!    the reference is loaded at that k via
//!    [`crate::ref_kmer_filter::RefKmerFilter::from_reference_fasta`]
//!    (mirrors `SeqSet refSet(9); refSet.InputRefFa(...)`,
//!    `FastqExtractor.cpp:273,302`) -- this is the caller's responsibility
//!    (the filter is passed in already loaded), NOT this function's.
//! 2. `hitLenRequired = 27` if paired, else `23` if single-end
//!    (`FastqExtractor.cpp:390-392`).
//! 3. The FIRST 1000 records of the read-1 source are sampled (READ-1 ONLY,
//!    even when paired -- `FastqExtractor.cpp:394-399` calls `reads.Next()`,
//!    never touching `mateReads`), summing `len += strlen(seq)` and
//!    counting `i` = the number ACTUALLY read (may be `< 1000` if the file
//!    is shorter). If `i == 0` (the read-1 file is empty), this is a hard
//!    error, matching `FastqExtractor.cpp:400-404`'s
//!    `"Read file is empty."` + `EXIT_FAILURE`.
//! 4. `if (len / (i * 5) > hitLenRequired) hitLenRequired = len / (i * 5)`
//!    (`FastqExtractor.cpp:405-406`) -- BOTH the comparison and the
//!    assignment use INTEGER division (`len`/`i` are C `int`s); this is
//!    reproduced with Rust integer division on `i64` (to avoid the
//!    intermediate `i * 5` multiplication ever wrapping an `i32` on a
//!    pathological input), not floating-point division truncated
//!    afterward -- these can differ (e.g. `len=10, i=3`: `10/(3*5) =
//!    10/15 = 0` either way here since both operands are small, but the
//!    reproduction uses the same operation order as the C++ regardless).
//! 5. `SetHitLenRequired(hitLenRequired)` /
//!    `SetRefSeqSimilarity(filterAlignmentSimilarity)` are applied to the
//!    caller's filter (`FastqExtractor.cpp:407-408`), then the read-1 source
//!    is REWOUND (`FastqExtractor.cpp:409`) so the actual extraction pass
//!    below re-reads from the beginning (the 1000-record sample is
//!    otherwise consumed/discarded, not reused as real extraction input).
//! 6. `InferKmerLength()` is computed; if it is STRICTLY GREATER than the
//!    current `kmerLength` (9), the filter is rebuilt at that new k via
//!    `UpdateKmerLength`, and if the new k also exceeds `hitLenRequired`,
//!    `hitLenRequired` is bumped up to match it too
//!    (`FastqExtractor.cpp:411-418`). For every reference this port has
//!    been validated against (KIR/HLA RNA references), `InferKmerLength()`
//!    evaluates to `8`, which is NOT greater than the initial `9` -- so this
//!    branch is DEAD IN PRACTICE for those inputs, but it is ported fully
//!    and generally here (not hardcoded away), since a sufficiently large
//!    future reference set COULD push `InferKmerLength()` above 9.
//!
//! # Per-pair decision: exact short-circuit semantics
//!
//! The single-threaded reference path (`FastqExtractor.cpp:447-480`, the
//! ONLY semantics this module reproduces -- see "Output order == input
//! order" below for why the multi-threaded path's OUTPUT is provably
//! identical regardless) evaluates, for each pair/read: `IsGoodCandidate`
//! on read 1 first; read 2's `IsGoodCandidate` is evaluated ONLY if read 1
//! was NOT already a good candidate AND the input is paired
//! (`FastqExtractor.cpp:463-468`: `if (!goodCandidate && hasMate &&
//! IsGoodCandidate(mateReads.seq, ...))`) -- i.e. read 2 is never even
//! filter-evaluated when read 1 already passed. [`extract_candidates`]'s
//! inner loop reproduces this exact short-circuit via Rust's `||`
//! short-circuit evaluation.
//!
//! # Output order == input order, at ANY thread count -- why this port's
//! `-t N` is byte-identical to `-t 1` and to the oracle
//!
//! `FastqExtractor.cpp` has two code paths, gated on `threadCnt`:
//! - `threadCnt == 1` (`FastqExtractor.cpp:447-480`): a single `while
//!   (reads.Next())` loop that decides AND immediately writes each pair, in
//!   file order, by construction.
//! - `threadCnt > 1` (`FastqExtractor.cpp:481-567`): reads a batch of up to
//!   `512 * threadCnt` pairs in FILE ORDER via `GetBatch`
//!   (`FastqExtractor.cpp:521-524`, sequential, not parallelized), spawns
//!   `threadCnt` worker threads that each only handle indices `i` where `i %
//!   threadCnt == tid` (`FastqExtractor.cpp:205`, so the PARALLEL work is
//!   confined to the per-read `IsGoodCandidate` decision, writing the result
//!   into `readBatch[i].id[0] = '\0'` as a reject marker -- no output I/O
//!   happens inside the worker threads at all), joins all threads
//!   (`FastqExtractor.cpp:552-553`), and only THEN emits output with a
//!   single sequential `for (i = 0; i < batchSize; ++i)` loop over the
//!   now-fully-decided batch, in batch order (`FastqExtractor.cpp:555-566`)
//!   -- and batches themselves are read and processed strictly in file
//!   order, one after another. So multi-threading only parallelizes the
//!   per-read FILTER DECISION (a pure, side-effect-free boolean
//!   computation, order-independent by construction), never the emission
//!   order, which is ALWAYS sequential-in-input-order on both paths.
//!
//! [`extract_candidates_with_threads`] follows exactly the same pattern:
//! `threads <= 1` runs the plain sequential decide-and-emit loop
//! ([`run_sequential`], byte-identical to [`extract_candidates`]'s prior
//! single-threaded-only behavior); `threads > 1` reads a bounded batch (`512
//! * threads` pairs, mirroring stock's own batch size), evaluates every
//! pair's [`is_good_pair`] decision across a `rayon` thread pool (each worker
//! reusing its own [`Scratch`] via `map_init`, so no output-affecting state
//! is EVER shared across threads -- only the immutable [`RefKmerFilter`] is,
//! which is `Sync` because [`RefKmerFilter::is_good_candidate_with_scratch`]
//! takes `&self`), collects the decisions into an order-preserving `Vec`
//! (`rayon`'s `collect()` on an `IndexedParallelIterator` -- which
//! `Vec::par_iter()`'s `.map_init()` is -- preserves input order), then emits
//! passing pairs to `sink` SEQUENTIALLY in original batch order (see
//! [`run_batch_parallel`]). This makes `-t N` output byte-identical to `-t 1`
//! output (and therefore to the oracle) at ANY `N`, for exactly the same
//! reason stock's own `-t N` output is threadCnt-invariant.

use crate::fastq::{FastqReader, FastqRecord};
use crate::ref_kmer_filter::{RefKmerFilter, Scratch};
use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use std::path::Path;

/// Batch size used by the parallel (`threads > 1`) extraction path, mirroring
/// `FastqExtractor.cpp:521-524`'s `512 * threadCnt` batch-read size: enough
/// pairs per batch to keep every worker busy without unbounded memory growth
/// on very large inputs. Purely a throughput/memory knob -- it has NO effect
/// on output (see [`extract_candidates`]'s module docs: batches are still
/// processed strictly in file order, one after another, and decisions within
/// a batch are collected order-preserving before any emission happens).
const PARALLEL_BATCH_SIZE_PER_THREAD: usize = 512;

/// `FastqExtractor.cpp:390`: base `hitLenRequired` for paired input.
const HIT_LEN_REQUIRED_PAIRED: i32 = 27;
/// `FastqExtractor.cpp:392`: base `hitLenRequired` for single-end input.
const HIT_LEN_REQUIRED_SINGLE: i32 = 23;
/// `FastqExtractor.cpp:394`: the sampling pass reads at most this many
/// read-1 records to compute the data-dependent `hitLenRequired` bump.
const HIT_LEN_SAMPLE_SIZE: usize = 1000;
/// `FastqExtractor.cpp:283`: default `-s` (`filterAlignmentSimilarity`).
pub const DEFAULT_REF_SEQ_SIMILARITY: f64 = 0.8;

/// A single extracted read record: `id`, `seq` (raw bases, no line
/// wrapping), and `qual` (`None` for FASTA-sourced input). Mirrors
/// [`crate::fastq::FastqRecord`] but is re-exported under this module's own
/// name as the stable public type a future fused `genotype` command
/// consumes (decoupling downstream callers from `fastq`'s internal read
/// representation, even though the two happen to be structurally
/// equivalent today).
pub type ReadRecord = FastqRecord;

/// Where [`extract_candidates`] sends reads that pass the candidate filter.
/// Implementations decide what "emit" means: the CLI's `extract` subcommand
/// writes FASTQ files (`OutputSeq`-equivalent, see [`crate::extract`]'s
/// module docs and the `unum` binary crate's `extract` stage module); a
/// future fused `genotype` command could instead hand `r1`/`r2` directly to
/// the next in-memory pipeline stage.
///
/// `r2` is `None` for single-end input (`hasMate == false`); for paired
/// input, `r2` is always `Some` (a good candidate pair is only ever formed
/// from both mates being read together, matching
/// `FastqExtractor.cpp:472-473`'s `if (hasMate) OutputSeq(fp2, ...)`).
pub trait CandidateSink {
    /// Emits one passing pair (or single read, when `r2` is `None`).
    ///
    /// # Errors
    ///
    /// Implementations may fail (e.g. a file-writing sink on an I/O error).
    fn emit_pair(&mut self, r1: &ReadRecord, r2: Option<&ReadRecord>) -> Result<()>;
}

/// Whether the read source is paired or single-end -- determines the base
/// `hitLenRequired` (`FastqExtractor.cpp:390-392`) and the short-circuit
/// per-pair decision logic (`FastqExtractor.cpp:463-468`).
pub enum ReadSource {
    /// Single-end: one FASTQ file.
    Single(FastqReader),
    /// Paired: two FASTQ files (mate 1, mate 2).
    Paired(FastqReader, FastqReader),
}

/// Extraction run summary, returned by [`extract_candidates`] for caller
/// diagnostics/logging (not part of T1K's own output, which is silent on
/// this beyond the `PrintLog` timing lines this port does not replicate --
/// see the CLI's `extract` subcommand for user-facing logging built on top
/// of this).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExtractMetrics {
    /// Total pairs/reads read from the source.
    pub total_reads: u64,
    /// Number of pairs/reads that passed the candidate filter and were
    /// emitted.
    pub candidates_emitted: u64,
    /// The final `kmerLength` the filter ended up using (after the
    /// `InferKmerLength`/`UpdateKmerLength` step).
    pub kmer_length: usize,
    /// The final `hitLenRequired` the filter ended up using.
    pub hit_len_required: i32,
}

/// Runs the full `FastqExtractor.cpp:main` read-extraction driver: performs
/// the data-dependent setup (`hitLenRequired` sampling,
/// `InferKmerLength`/`UpdateKmerLength`, see module docs), then streams
/// every pair/read from `source`, applying the short-circuit
/// `IsGoodCandidate` decision, and emits passing pairs to `sink` in input
/// order.
///
/// `filter` must already be loaded (via
/// [`crate::ref_kmer_filter::RefKmerFilter::from_reference_fasta`]) at the
/// initial k-mer length (9, matching `FastqExtractor.cpp:272-273`'s literal
/// `SeqSet refSet(9)`) -- this function does not load the reference itself,
/// since [`RefKmerFilter`] construction requires a reference FASTA path the
/// caller already has (mirroring the CLI's `-f` flag).
///
/// `ref_seq_similarity` corresponds to the `-s` CLI flag (default
/// [`DEFAULT_REF_SEQ_SIMILARITY`]).
///
/// # Errors
///
/// Returns an error if the read-1 source is empty (`i == 0` after sampling,
/// matching `FastqExtractor.cpp:400-404`'s `"Read file is empty."`), if a
/// paired source's two mate files have different read counts (matching
/// `FastqExtractor.cpp:451-455`), or if the underlying FASTQ readers hit a
/// parse error.
pub fn extract_candidates(
    source: &mut ReadSource,
    filter: &mut RefKmerFilter,
    ref_seq_similarity: f64,
    sink: &mut impl CandidateSink,
) -> Result<ExtractMetrics> {
    extract_candidates_with_threads(source, filter, ref_seq_similarity, 1, sink)
}

/// Same as [`extract_candidates`], but runs the per-pair candidate DECISION
/// (`IsGoodCandidate` on mate 1, short-circuit-then-mate-2) across `threads`
/// worker threads via a scoped `rayon` thread pool, while keeping emission
/// strictly sequential in input order -- see the module docs ("Output order
/// == input order") for why this makes the output byte-identical to the
/// `threads == 1` path (and to the oracle) at ANY thread count.
///
/// `threads <= 1` takes the exact sequential fast path (no rayon pool is
/// built at all), matching [`extract_candidates`]'s prior behavior exactly.
/// `threads > 1` reads pairs into a bounded batch (up to
/// `threads * `[`PARALLEL_BATCH_SIZE_PER_THREAD`]` pairs per batch, mirroring
/// `FastqExtractor.cpp`'s own `512 * threadCnt` batching -- module docs),
/// evaluates every pair's candidate decision in parallel (each worker reusing
/// its own [`Scratch`] via `rayon`'s `map_init`, so `has_hit_in_set`'s hot
/// buffers are never reallocated per-read even under threading -- folding in
/// the per-thread scratch-buffer design from the single-threaded scratch-
/// reuse commits), collects the decisions in an order-preserving `Vec`
/// (`rayon`'s parallel iterators preserve input order on `collect()`), then
/// emits passing pairs to `sink` sequentially over the batch in original
/// order. Batches themselves are still read and processed strictly one after
/// another in file order, so no output can ever depend on `threads`.
///
/// # Errors
///
/// See [`extract_candidates`].
pub fn extract_candidates_with_threads(
    source: &mut ReadSource,
    filter: &mut RefKmerFilter,
    ref_seq_similarity: f64,
    threads: usize,
    sink: &mut impl CandidateSink,
) -> Result<ExtractMetrics> {
    let has_mate = matches!(source, ReadSource::Paired(_, _));

    // Step 2: base hitLenRequired (FastqExtractor.cpp:390-392).
    let mut hit_len_required =
        if has_mate { HIT_LEN_REQUIRED_PAIRED } else { HIT_LEN_REQUIRED_SINGLE };

    // Step 3: sample the first HIT_LEN_SAMPLE_SIZE read-1 records
    // (FastqExtractor.cpp:394-399) -- READ-1 ONLY, never read-2, even when
    // paired.
    let read1 = match source {
        ReadSource::Single(r) => r,
        ReadSource::Paired(r1, _) => r1,
    };
    let (sampled_count, sampled_len) = sample_hit_len_required_stats(read1)?;
    ensure!(sampled_count > 0, "Read file is empty.");

    // Step 4: integer-division bump (FastqExtractor.cpp:405-406).
    let candidate = sampled_len / (i64::try_from(sampled_count).unwrap_or(i64::MAX) * 5);
    if candidate > i64::from(hit_len_required) {
        // `candidate` is derived from real read-length sums over at most
        // 1000 reads, so it comfortably fits in an i32 in practice; clamp
        // defensively rather than panic on a pathological/adversarial input.
        hit_len_required = i32::try_from(candidate).unwrap_or(i32::MAX);
    }

    // Step 5: apply to the filter, then rewind read-1
    // (FastqExtractor.cpp:407-409).
    filter.set_hit_len_required(hit_len_required);
    filter.set_ref_seq_similarity(ref_seq_similarity);
    let read1 = match source {
        ReadSource::Single(r) => r,
        ReadSource::Paired(r1, _) => r1,
    };
    read1.rewind().context("rewinding read-1 source after hitLenRequired sampling")?;

    // Step 6: InferKmerLength / conditional UpdateKmerLength
    // (FastqExtractor.cpp:411-418). Dead branch for every reference this
    // port has been validated against (see module docs), but ported fully
    // and generally -- no hardcoding of "kmer length stays 9".
    let inferred = filter.infer_kmer_length();
    if inferred > filter.kmer_length() {
        filter.update_kmer_length(inferred);
        if inferred > usize::try_from(hit_len_required).unwrap_or(0) {
            hit_len_required = i32::try_from(inferred).unwrap_or(i32::MAX);
            filter.set_hit_len_required(hit_len_required);
        }
    }

    // Main extraction pass (FastqExtractor.cpp:447-480, single-threaded
    // reference semantics -- see module docs for why this is byte-identical
    // to the oracle's output at ANY -t). `threads <= 1` runs the plain
    // sequential decide-and-emit loop (unchanged from before this function
    // gained a `threads` parameter); `threads > 1` batches reads and
    // parallelizes only the decision (see `run_batch_parallel`'s doc
    // comment).
    let (total_reads, candidates_emitted) = if threads <= 1 {
        run_sequential(source, filter, sink)?
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .context("building rayon thread pool for parallel extraction")?;
        run_batch_parallel(source, filter, &pool, threads, sink)?
    };

    // DIVERGENCE from T1K (see docs/DIVERGENCES.md): stock loops on mate-1
    // (`FastqExtractor.cpp:449-479`) and errors only when mate-2 is *exhausted
    // early* (the in-loop check above). Once mate-1 is exhausted it stops,
    // silently ignoring any extra trailing mate-2 records and exiting 0 -- so a
    // truncated mate-1 file drops reads with no warning. `unum` adds the
    // symmetric check: after mate-1 is exhausted, mate-2 must be exhausted too,
    // otherwise the two files have unequal read counts and the input is
    // malformed -- the same hard error as the mate-2-exhausted-early case.
    if let ReadSource::Paired(_, r2) = source {
        if r2.next_record().context("reading mate-2 read")?.is_some() {
            bail!("The two mate-pair read files have different number of reads.");
        }
    }

    Ok(ExtractMetrics {
        total_reads,
        candidates_emitted,
        kmer_length: filter.kmer_length(),
        hit_len_required,
    })
}

/// Reads one pair/read from `source`. `Ok(None)` means mate-1 (or the
/// single-end source) is exhausted -- the caller should stop. An error is
/// returned if mate-2 is exhausted before mate-1 (matching
/// `FastqExtractor.cpp:451-455`'s "different number of reads" error).
fn read_next_pair(source: &mut ReadSource) -> Result<Option<(FastqRecord, Option<FastqRecord>)>> {
    match source {
        ReadSource::Single(r) => {
            let Some(rec1) = r.next_record().context("reading single-end read")? else {
                return Ok(None);
            };
            Ok(Some((rec1, None)))
        }
        ReadSource::Paired(r1, r2) => {
            let Some(rec1) = r1.next_record().context("reading mate-1 read")? else {
                return Ok(None);
            };
            let Some(rec2_val) = r2.next_record().context("reading mate-2 read")? else {
                bail!("The two mate-pair read files have different number of reads.");
            };
            Ok(Some((rec1, Some(rec2_val))))
        }
    }
}

/// Evaluates the per-pair candidate decision for a single pair/read,
/// reproducing `FastqExtractor.cpp:463-468`'s exact short-circuit: read 2 is
/// never filter-evaluated once read 1 already passed. Shared by both the
/// sequential and parallel paths so the decision logic itself has exactly
/// one implementation.
fn is_good_pair(
    filter: &RefKmerFilter,
    rec1: &FastqRecord,
    rec2: Option<&FastqRecord>,
    scratch: &mut Scratch,
) -> bool {
    filter.is_good_candidate_with_scratch(&rec1.seq, scratch)
        || rec2.is_some_and(|r2| filter.is_good_candidate_with_scratch(&r2.seq, scratch))
}

/// The plain sequential decide-and-emit loop: reads, decides, and (if
/// passing) emits each pair one at a time, in file order. Returns
/// `(total_reads, candidates_emitted)`.
fn run_sequential(
    source: &mut ReadSource,
    filter: &RefKmerFilter,
    sink: &mut impl CandidateSink,
) -> Result<(u64, u64)> {
    let mut scratch = Scratch::default();
    let mut total_reads: u64 = 0;
    let mut candidates_emitted: u64 = 0;

    while let Some((rec1, rec2)) = read_next_pair(source)? {
        total_reads += 1;
        if is_good_pair(filter, &rec1, rec2.as_ref(), &mut scratch) {
            sink.emit_pair(&rec1, rec2.as_ref())?;
            candidates_emitted += 1;
        }
    }

    Ok((total_reads, candidates_emitted))
}

/// The parallel-decision path: reads a bounded batch of pairs in file order,
/// evaluates every pair's [`is_good_pair`] decision across `pool`'s worker
/// threads (each worker reusing its own [`Scratch`] via `rayon::map_init`),
/// collects the decisions into an order-preserving `Vec` (`rayon`'s
/// `collect()` on an `IndexedParallelIterator` preserves input order), then
/// emits passing pairs to `sink` SEQUENTIALLY over the batch in original
/// order -- never in parallel. Batches are processed one after another,
/// strictly in file order, so this is byte-identical to [`run_sequential`]'s
/// output regardless of `threads`/batch size. Returns
/// `(total_reads, candidates_emitted)`.
fn run_batch_parallel(
    source: &mut ReadSource,
    filter: &RefKmerFilter,
    pool: &rayon::ThreadPool,
    threads: usize,
    sink: &mut impl CandidateSink,
) -> Result<(u64, u64)> {
    let batch_capacity = threads.saturating_mul(PARALLEL_BATCH_SIZE_PER_THREAD).max(1);
    let mut total_reads: u64 = 0;
    let mut candidates_emitted: u64 = 0;

    loop {
        // Read a batch in strict file order (sequential; only the DECISION
        // below is parallelized).
        let mut batch: Vec<(FastqRecord, Option<FastqRecord>)> = Vec::with_capacity(batch_capacity);
        while batch.len() < batch_capacity {
            match read_next_pair(source)? {
                Some(pair) => batch.push(pair),
                None => break,
            }
        }
        if batch.is_empty() {
            break;
        }
        total_reads += u64::try_from(batch.len()).unwrap_or(u64::MAX);

        // Parallel decision, order-preserving collect.
        let decisions: Vec<bool> = pool.install(|| {
            batch
                .par_iter()
                .map_init(Scratch::default, |scratch, (rec1, rec2)| {
                    is_good_pair(filter, rec1, rec2.as_ref(), scratch)
                })
                .collect()
        });

        // Sequential emit, in original batch (== file) order.
        for ((rec1, rec2), good) in batch.into_iter().zip(decisions) {
            if good {
                sink.emit_pair(&rec1, rec2.as_ref())?;
                candidates_emitted += 1;
            }
        }
    }

    Ok((total_reads, candidates_emitted))
}

/// Samples up to [`HIT_LEN_SAMPLE_SIZE`] records from `read1`, returning
/// `(count actually read, sum of sequence lengths)` -- ports
/// `FastqExtractor.cpp:394-399`'s sampling loop exactly (including stopping
/// early if the file has fewer than 1000 records, via `reads.Next()`
/// returning false).
fn sample_hit_len_required_stats(read1: &mut FastqReader) -> Result<(usize, i64)> {
    let mut count = 0usize;
    let mut total_len: i64 = 0;
    for _ in 0..HIT_LEN_SAMPLE_SIZE {
        let Some(rec) = read1.next_record().context("sampling read-1 for hitLenRequired")? else {
            break;
        };
        total_len += i64::try_from(rec.seq.len()).unwrap_or(i64::MAX);
        count += 1;
    }
    Ok((count, total_len))
}

/// Writes a candidate pair/read to a paired or single-end pair of FASTQ
/// (or, if `qual` is absent, FASTA) writers, ported from `OutputSeq`
/// (`FastqExtractor.cpp:120-153`).
///
/// # Default (no trim): `start == 0 && end == -1`
///
/// Writes the FULL sequence: `@{id}\n{seq}\n+\n{qual}\n` for FASTQ (a BARE
/// `+` line, no repeated id after it -- `FastqExtractor.cpp:125`), or
/// `>{id}\n{seq}\n` for FASTA (`qual == None`, `FastqExtractor.cpp:127`).
///
/// # Trim (`start`/`end` other than the default)
///
/// Writes the substring `seq[start..=end]` (INCLUSIVE end, matching C++'s
/// `for (i = s; i <= e; ++i)`, `FastqExtractor.cpp:135-136`), where `end ==
/// -1` means `seq.len() - 1` (i.e. "through the end", matching
/// `FastqExtractor.cpp:133`'s `end == -1 ? strlen(seq) - 1 : end`). Applies
/// the same substring range to `qual` if present
/// (`FastqExtractor.cpp:147-149`). This path is ported for generality (the
/// `--read1Start`/`--read1End`/`--read2Start`/`--read2End` CLI flags) but is
/// NOT exercised by the byte-identity differential, which always uses the
/// default (no trim).
///
/// # Errors
///
/// Returns an error if writing to `out` fails.
pub fn output_seq(
    out: &mut impl std::io::Write,
    id: &str,
    seq: &[u8],
    qual: Option<&[u8]>,
    start: i64,
    end: i64,
) -> Result<()> {
    if start == 0 && end == -1 {
        if let Some(qual) = qual {
            writeln!(out, "@{id}")?;
            out.write_all(seq)?;
            write!(out, "\n+\n")?;
            out.write_all(qual)?;
            writeln!(out)?;
        } else {
            writeln!(out, ">{id}")?;
            out.write_all(seq)?;
            writeln!(out)?;
        }
        return Ok(());
    }

    let s = usize::try_from(start).context("output_seq: negative trim start")?;
    let e = if end == -1 {
        seq.len().checked_sub(1).context("output_seq: empty sequence with end == -1")?
    } else {
        usize::try_from(end).context("output_seq: negative trim end")?
    };
    ensure!(
        s <= e && e < seq.len(),
        "output_seq: trim range [{s}, {e}] out of bounds for len {}",
        seq.len()
    );

    let seq_slice = &seq[s..=e];
    if let Some(qual) = qual {
        writeln!(out, "@{id}")?;
        out.write_all(seq_slice)?;
        write!(out, "\n+\n")?;
        out.write_all(&qual[s..=e])?;
        writeln!(out)?;
    } else {
        writeln!(out, ">{id}")?;
        out.write_all(seq_slice)?;
        writeln!(out)?;
    }
    Ok(())
}

/// Convenience: constructs a [`ReadSource`] for either a single-end file or
/// a paired pair of files.
///
/// # Errors
///
/// Returns an error if any file cannot be opened.
pub fn open_source(mate1: &Path, mate2: Option<&Path>) -> Result<ReadSource> {
    if let Some(mate2) = mate2 {
        let r1 = FastqReader::open(mate1)
            .with_context(|| format!("opening mate-1 file {}", mate1.display()))?;
        let r2 = FastqReader::open(mate2)
            .with_context(|| format!("opening mate-2 file {}", mate2.display()))?;
        Ok(ReadSource::Paired(r1, r2))
    } else {
        let r1 = FastqReader::open(mate1)
            .with_context(|| format!("opening single-end file {}", mate1.display()))?;
        Ok(ReadSource::Single(r1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn write_fastq(path: &Path, records: &[(&str, &str, &str)]) {
        let mut f = std::fs::File::create(path).unwrap();
        for (id, seq, qual) in records {
            writeln!(f, "@{id}\n{seq}\n+\n{qual}").unwrap();
        }
    }

    fn write_ref_fasta(path: &Path, records: &[(&str, &str)]) {
        let mut f = std::fs::File::create(path).unwrap();
        for (id, seq) in records {
            writeln!(f, ">{id}\n{seq}").unwrap();
        }
    }

    struct VecSink {
        pairs: Vec<(ReadRecord, Option<ReadRecord>)>,
    }

    impl CandidateSink for VecSink {
        fn emit_pair(&mut self, r1: &ReadRecord, r2: Option<&ReadRecord>) -> Result<()> {
            self.pairs.push((r1.clone(), r2.cloned()));
            Ok(())
        }
    }

    // A repetitive-but-balanced 200bp reference, long enough that exact
    // 100bp substrings comfortably clear both HasHitInSet gates at the
    // default hitLenRequired (27 for paired / 23 single-end) and default
    // ref_seq_similarity (0.8).
    const REF_SEQ: &str = "ACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATGGTACGTACGTACGTACG";

    fn ref_read1(len: usize, offset: usize) -> String {
        REF_SEQ[offset..offset + len].to_string()
    }

    #[test]
    fn hit_len_required_sampling_uses_integer_division_and_base_by_pairing() {
        // Directly exercises sample_hit_len_required_stats plus the
        // integer-division bump math, without going through a full
        // extract_candidates run (keeps this test fast and focused).
        let tmp = tempfile::tempdir().unwrap();
        let r1_path = tmp.path().join("r1.fq");
        // 10 reads of length 50 each: len=500, i=10 -> 500/(10*5)=500/50=10.
        // hitLenRequired base (paired)=27, candidate=10 -> stays 27.
        let records: Vec<(String, String, String)> =
            (0..10).map(|i| (format!("r{i}"), "A".repeat(50), "I".repeat(50))).collect();
        let record_refs: Vec<(&str, &str, &str)> =
            records.iter().map(|(a, b, c)| (a.as_str(), b.as_str(), c.as_str())).collect();
        write_fastq(&r1_path, &record_refs);

        let mut r1 = FastqReader::open(&r1_path).unwrap();
        let (count, len) = sample_hit_len_required_stats(&mut r1).unwrap();
        assert_eq!(count, 10);
        assert_eq!(len, 500);
        let candidate = len / (i64::try_from(count).unwrap() * 5);
        assert_eq!(candidate, 10);
        assert!(candidate <= i64::from(HIT_LEN_REQUIRED_PAIRED));
    }

    #[test]
    fn hit_len_required_bump_when_reads_are_long() {
        // Long reads push len/(i*5) above the paired base of 27: e.g. 10
        // reads of length 200 -> len=2000, i=10 -> 2000/50=40 > 27.
        let tmp = tempfile::tempdir().unwrap();
        let r1_path = tmp.path().join("r1.fq");
        let records: Vec<(String, String, String)> =
            (0..10).map(|i| (format!("r{i}"), "A".repeat(200), "I".repeat(200))).collect();
        let record_refs: Vec<(&str, &str, &str)> =
            records.iter().map(|(a, b, c)| (a.as_str(), b.as_str(), c.as_str())).collect();
        write_fastq(&r1_path, &record_refs);

        let mut r1 = FastqReader::open(&r1_path).unwrap();
        let (count, len) = sample_hit_len_required_stats(&mut r1).unwrap();
        let candidate = len / (i64::try_from(count).unwrap() * 5);
        assert_eq!(candidate, 40);
        assert!(candidate > i64::from(HIT_LEN_REQUIRED_PAIRED));
    }

    #[test]
    fn hit_len_required_sampling_stops_before_1000_on_short_file() {
        let tmp = tempfile::tempdir().unwrap();
        let r1_path = tmp.path().join("r1.fq");
        // Only 3 records, far fewer than the 1000-record sample cap.
        let records: [(&str, &str, &str); 3] =
            [("r0", "ACGT", "IIII"), ("r1", "ACGT", "IIII"), ("r2", "ACGT", "IIII")];
        write_fastq(&r1_path, &records);

        let mut r1 = FastqReader::open(&r1_path).unwrap();
        let (count, len) = sample_hit_len_required_stats(&mut r1).unwrap();
        assert_eq!(count, 3, "must stop at actual EOF, not require 1000 reads");
        assert_eq!(len, 12);
    }

    #[test]
    fn empty_read1_file_is_an_error() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        write_ref_fasta(&ref_path, &[("only", REF_SEQ)]);
        let r1_path = tmp.path().join("empty.fq");
        std::fs::write(&r1_path, "").unwrap();

        let mut filter = RefKmerFilter::from_reference_fasta(&ref_path, 9).unwrap();
        let mut source = open_source(&r1_path, None).unwrap();
        let mut sink = VecSink { pairs: Vec::new() };
        let result =
            extract_candidates(&mut source, &mut filter, DEFAULT_REF_SEQ_SIMILARITY, &mut sink);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Read file is empty."));
    }

    #[test]
    fn mate_count_mismatch_mate2_shorter_is_an_error() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        write_ref_fasta(&ref_path, &[("only", REF_SEQ)]);

        let r1_path = tmp.path().join("r1.fq");
        let r2_path = tmp.path().join("r2.fq");
        let r1_records: Vec<(String, String, String)> =
            (0..5).map(|i| (format!("p{i}"), ref_read1(60, 0), "I".repeat(60))).collect();
        let r1_refs: Vec<(&str, &str, &str)> =
            r1_records.iter().map(|(a, b, c)| (a.as_str(), b.as_str(), c.as_str())).collect();
        write_fastq(&r1_path, &r1_refs);

        let r2_records: Vec<(String, String, String)> =
            (0..3).map(|i| (format!("p{i}"), ref_read1(60, 20), "I".repeat(60))).collect();
        let r2_refs: Vec<(&str, &str, &str)> =
            r2_records.iter().map(|(a, b, c)| (a.as_str(), b.as_str(), c.as_str())).collect();
        write_fastq(&r2_path, &r2_refs);

        let mut filter = RefKmerFilter::from_reference_fasta(&ref_path, 9).unwrap();
        let mut source = open_source(&r1_path, Some(&r2_path)).unwrap();
        let mut sink = VecSink { pairs: Vec::new() };
        let result =
            extract_candidates(&mut source, &mut filter, DEFAULT_REF_SEQ_SIMILARITY, &mut sink);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("different number of reads"));
    }

    #[test]
    fn mate_count_mismatch_mate2_longer_is_an_error() {
        // DIVERGENCE from T1K (see docs/DIVERGENCES.md): stock loops on mate-1
        // (`FastqExtractor.cpp:449-479`) and silently ignores extra trailing
        // mate-2 records, exiting 0 -- so a truncated mate-1 file drops reads
        // with no warning. `unum` treats a mate-2-longer file as the same hard
        // error as mate-2-shorter: an unequal-read-count mate pair is malformed
        // input either way.
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        write_ref_fasta(&ref_path, &[("only", REF_SEQ)]);

        let r1_path = tmp.path().join("r1.fq");
        let r2_path = tmp.path().join("r2.fq");
        let r1_records: Vec<(String, String, String)> =
            (0..3).map(|i| (format!("p{i}"), ref_read1(60, 0), "I".repeat(60))).collect();
        let r1_refs: Vec<(&str, &str, &str)> =
            r1_records.iter().map(|(a, b, c)| (a.as_str(), b.as_str(), c.as_str())).collect();
        write_fastq(&r1_path, &r1_refs);

        let r2_records: Vec<(String, String, String)> =
            (0..5).map(|i| (format!("p{i}"), ref_read1(60, 20), "I".repeat(60))).collect();
        let r2_refs: Vec<(&str, &str, &str)> =
            r2_records.iter().map(|(a, b, c)| (a.as_str(), b.as_str(), c.as_str())).collect();
        write_fastq(&r2_path, &r2_refs);

        let mut filter = RefKmerFilter::from_reference_fasta(&ref_path, 9).unwrap();
        let mut source = open_source(&r1_path, Some(&r2_path)).unwrap();
        let mut sink = VecSink { pairs: Vec::new() };
        let result =
            extract_candidates(&mut source, &mut filter, DEFAULT_REF_SEQ_SIMILARITY, &mut sink);
        assert!(result.is_err(), "mate-2-longer must be a hard error, not a silent drop");
        assert!(result.unwrap_err().to_string().contains("different number of reads"));
    }

    #[test]
    fn short_circuit_skips_mate2_evaluation_when_mate1_already_good() {
        // If mate1 matches the reference but mate2 is total noise, the pair
        // must still be emitted (mate1 alone is sufficient) -- and per the
        // short-circuit semantics, mate2's IsGoodCandidate is never even
        // evaluated (though that's only externally observable via the
        // outcome here, not a side channel).
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        write_ref_fasta(&ref_path, &[("only", REF_SEQ)]);

        let r1_path = tmp.path().join("r1.fq");
        let r2_path = tmp.path().join("r2.fq");
        write_fastq(&r1_path, &[("p0", &ref_read1(100, 0), &"I".repeat(100))]);
        // mate2: low-complexity noise, would never be a candidate on its own.
        write_fastq(&r2_path, &[("p0", &"A".repeat(100), &"I".repeat(100))]);

        let mut filter = RefKmerFilter::from_reference_fasta(&ref_path, 9).unwrap();
        let mut source = open_source(&r1_path, Some(&r2_path)).unwrap();
        let mut sink = VecSink { pairs: Vec::new() };
        let metrics =
            extract_candidates(&mut source, &mut filter, DEFAULT_REF_SEQ_SIMILARITY, &mut sink)
                .unwrap();
        assert_eq!(metrics.total_reads, 1);
        assert_eq!(metrics.candidates_emitted, 1);
        assert_eq!(sink.pairs.len(), 1);
    }

    #[test]
    fn mate2_only_hit_still_emits_pair() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        write_ref_fasta(&ref_path, &[("only", REF_SEQ)]);

        let r1_path = tmp.path().join("r1.fq");
        let r2_path = tmp.path().join("r2.fq");
        write_fastq(&r1_path, &[("p0", &"A".repeat(100), &"I".repeat(100))]);
        write_fastq(&r2_path, &[("p0", &ref_read1(100, 40), &"I".repeat(100))]);

        let mut filter = RefKmerFilter::from_reference_fasta(&ref_path, 9).unwrap();
        let mut source = open_source(&r1_path, Some(&r2_path)).unwrap();
        let mut sink = VecSink { pairs: Vec::new() };
        let metrics =
            extract_candidates(&mut source, &mut filter, DEFAULT_REF_SEQ_SIMILARITY, &mut sink)
                .unwrap();
        assert_eq!(metrics.candidates_emitted, 1);
    }

    #[test]
    fn neither_mate_hits_pair_is_dropped() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        write_ref_fasta(&ref_path, &[("only", REF_SEQ)]);

        let r1_path = tmp.path().join("r1.fq");
        let r2_path = tmp.path().join("r2.fq");
        write_fastq(&r1_path, &[("p0", &"A".repeat(100), &"I".repeat(100))]);
        write_fastq(&r2_path, &[("p0", &"T".repeat(100), &"I".repeat(100))]);

        let mut filter = RefKmerFilter::from_reference_fasta(&ref_path, 9).unwrap();
        let mut source = open_source(&r1_path, Some(&r2_path)).unwrap();
        let mut sink = VecSink { pairs: Vec::new() };
        let metrics =
            extract_candidates(&mut source, &mut filter, DEFAULT_REF_SEQ_SIMILARITY, &mut sink)
                .unwrap();
        assert_eq!(metrics.total_reads, 1);
        assert_eq!(metrics.candidates_emitted, 0);
        assert!(sink.pairs.is_empty());
    }

    #[test]
    fn both_output_mates_use_read1_id() {
        // FastqExtractor.cpp:471-473: OutputSeq is called with `reads.id`
        // for BOTH fp1 AND fp2 -- mate2's own kseq-parsed id is never used
        // for output, even though it was read from mate2's own file (which
        // may have a different id string).
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        write_ref_fasta(&ref_path, &[("only", REF_SEQ)]);

        let r1_path = tmp.path().join("r1.fq");
        let r2_path = tmp.path().join("r2.fq");
        write_fastq(&r1_path, &[("mate1_id", &ref_read1(100, 0), &"I".repeat(100))]);
        write_fastq(
            &r2_path,
            &[("totally_different_mate2_id", &ref_read1(100, 20), &"I".repeat(100))],
        );

        let mut filter = RefKmerFilter::from_reference_fasta(&ref_path, 9).unwrap();
        let mut source = open_source(&r1_path, Some(&r2_path)).unwrap();
        let mut sink = VecSink { pairs: Vec::new() };
        extract_candidates(&mut source, &mut filter, DEFAULT_REF_SEQ_SIMILARITY, &mut sink)
            .unwrap();

        assert_eq!(sink.pairs.len(), 1);
        let (r1, r2) = &sink.pairs[0];
        assert_eq!(r1.id, "mate1_id");
        // The sink receives the raw records (mate2's own id intact); it is
        // the CALLER's (CLI's FASTQ-writing sink's) responsibility to use
        // r1.id for BOTH output files -- verified by the CLI-level
        // differential test and the CLI sink's own unit test, not here
        // (this test only proves the library hands both records through
        // unmodified so the sink can make that choice).
        assert_eq!(r2.as_ref().unwrap().id, "totally_different_mate2_id");
    }

    #[test]
    fn single_end_source_has_no_mate2() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        write_ref_fasta(&ref_path, &[("only", REF_SEQ)]);

        let r1_path = tmp.path().join("r1.fq");
        write_fastq(&r1_path, &[("s0", &ref_read1(100, 0), &"I".repeat(100))]);

        let mut filter = RefKmerFilter::from_reference_fasta(&ref_path, 9).unwrap();
        let mut source = open_source(&r1_path, None).unwrap();
        let mut sink = VecSink { pairs: Vec::new() };
        let metrics =
            extract_candidates(&mut source, &mut filter, DEFAULT_REF_SEQ_SIMILARITY, &mut sink)
                .unwrap();
        assert_eq!(metrics.candidates_emitted, 1);
        assert!(sink.pairs[0].1.is_none());
    }

    #[test]
    fn output_seq_default_writes_bare_plus_fastq() {
        let mut buf = Vec::new();
        output_seq(&mut buf, "r1", b"ACGT", Some(b"IIII"), 0, -1).unwrap();
        assert_eq!(std::str::from_utf8(&buf).unwrap(), "@r1\nACGT\n+\nIIII\n");
    }

    #[test]
    fn output_seq_default_fasta_when_no_qual() {
        let mut buf = Vec::new();
        output_seq(&mut buf, "r1", b"ACGT", None, 0, -1).unwrap();
        assert_eq!(std::str::from_utf8(&buf).unwrap(), ">r1\nACGT\n");
    }

    #[test]
    fn output_seq_trims_inclusive_range() {
        let mut buf = Vec::new();
        // seq "ACGTACGT" (len 8), trim [2, 5] inclusive -> "GTAC"
        output_seq(&mut buf, "r1", b"ACGTACGT", Some(b"IIIIJJJJ"), 2, 5).unwrap();
        assert_eq!(std::str::from_utf8(&buf).unwrap(), "@r1\nGTAC\n+\nIIJJ\n");
    }

    #[test]
    fn output_seq_trim_end_minus_one_means_through_end() {
        let mut buf = Vec::new();
        // start=2, end=-1 -> through the end of an 8bp seq -> [2..=7] = "GTACGT"
        output_seq(&mut buf, "r1", b"ACGTACGT", Some(b"IIIIJJJJ"), 2, -1).unwrap();
        assert_eq!(std::str::from_utf8(&buf).unwrap(), "@r1\nGTACGT\n+\nIIJJJJ\n");
    }

    #[test]
    fn output_seq_trim_fasta_no_qual() {
        let mut buf = Vec::new();
        output_seq(&mut buf, "r1", b"ACGTACGT", None, 1, 3).unwrap();
        assert_eq!(std::str::from_utf8(&buf).unwrap(), ">r1\nCGT\n");
    }
}
