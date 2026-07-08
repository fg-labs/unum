//! BAM/CRAM read-candidate extraction driver, ported from `BamExtractor.cpp`'s
//! `main` (`BamExtractor.cpp:464-949`): drives a configured
//! [`crate::ref_kmer_filter::RefKmerFilter`] plus a sorted list of gene
//! (`chrId`, `start`, `end`) intervals over a coordinate-sorted BAM/CRAM
//! ([`crate::alignments::Alignments`]) in a TWO-PASS scan, and emits
//! candidate reads (pairs where at least one mate is a good candidate, or an
//! unaligned template pair, or an on-target/alt-chrom aligned single-end
//! read) to a [`crate::extract::CandidateSink`], in BAM-encounter order.
//!
//! # Scope: barcode/UMI dropped, `-u`/`--mateIdSuffixLen` kept
//!
//! Stock `BamExtractor.cpp` also supports `--barcode`/`--UMI` aux-tag output
//! files (`fpBc`/`fpUMI`). Those paths are intentionally NOT ported here (no
//! `GetFieldZ` in [`crate::alignments::Alignments`] either, by the same
//! decision -- see that module's doc comment). `-u`
//! (`abnormalUnalignedFlag`) and `--mateIdSuffixLen` (`mateIdLen`) ARE
//! ported, since they affect which reads are extracted/how names are
//! matched, not just an optional side-output file.
//!
//! # Library-first design
//!
//! [`extract_from_bam`] is I/O-sink-agnostic (a [`crate::extract::CandidateSink`],
//! reused from the FASTQ-extractor port) and takes an already-open
//! [`crate::alignments::Alignments`] plus an already-loaded
//! [`crate::ref_kmer_filter::RefKmerFilter`] and `genes` interval list, so a
//! future fused genotype command can reuse the exact candidate-selection
//! logic. The CLI's `extract -b` mode is a thin wrapper that builds these
//! three inputs from `-f`/`-b` and a FASTQ-file-writing sink.
//!
//! # Threading: this port's own `-t N` is output-deterministic; the oracle's
//! `threadCnt == 1` semantics are what it reproduces
//!
//! `BamExtractor.cpp`'s OWN `threadCnt > 1` path changes BATCHING/OUTPUT-QUEUE
//! FLUSH TIMING for the unaligned-template-pair and single-end-unmapped
//! candidate paths (`ProcessUnmappedReads_Thread` +
//! `DistributeWork`/`AddWorkQueue`'s `workLoad = 2048` batching,
//! `BamExtractor.cpp:202-407`): candidates are queued up to 2048 at a time,
//! handed to a free worker thread, and flushed to `fp1`/`fp2` only once the
//! shared `outputQueue` exceeds `2 * candidates.size()`
//! (`BamExtractor.cpp:243`) -- i.e. the ORACLE's own output order/batching is
//! NOT provably thread-count-invariant the way `FastqExtractor.cpp`'s is.
//! This port's decision/emission SEMANTICS therefore only reproduce the
//! oracle's `threadCnt == 1` code path (`BamExtractor.cpp:675-696,754-778`:
//! the direct single-threaded `if`/`else` branches, never the
//! `AddWorkQueue`/`ProcessUnmappedReads_Thread` branches), and the golden
//! test always captures the `-t 1` output to match -- see
//! `crates/unum-core/tests/golden_bam_extract.rs`'s module docs.
//!
//! This port's OWN `-t N` (i.e. [`extract_from_bam_with_threads`]'s
//! `threads` parameter) is a SEPARATE, unrelated design: it parallelizes only
//! the per-read candidate-filter DECISION within pass 1 (the same
//! decision/emission-separation pattern [`crate::extract`] uses for
//! `FastqExtractor.cpp`), never pass 1's read-scan order, `tag`
//! gene-interval-pointer advance, `candidates`/`used_name` mutation, or
//! emission order -- ALL of which stay exactly as sequential (and therefore
//! exactly as BAM-encounter-order-deterministic) as the `threads == 1` path.
//! So this port's `-t N` output IS provably byte-identical to its own `-t 1`
//! output (and therefore to the oracle's `-t 1` output) at any `N`, even
//! though the ORACLE's own `-t N` is not thread-count-invariant against
//! itself. See [`run_pass1_with_threads`]'s doc comment for the
//! scan/evaluate/apply split that makes this possible.
//!
//! # Two-pass structure
//!
//! Both passes live in a single function, [`extract_from_bam`]. Pass 1
//! streams the BAM once, coordinate order assumed (`BamExtractor.cpp:631`'s
//! comment: "assuming the input is sorted by coordinate"):
//! - Unaligned-template pairs (paired input, `-u` NOT given): emitted
//!   directly, reading the mate as the very NEXT record
//!   (`BamExtractor.cpp:646-728`).
//! - Alt-chrom aligned reads and any other non-template-aligned reads
//!   (`ValidAlternativeChrom`): for paired input, their (trimmed) name is
//!   recorded in `candidates` (mate sequences filled in pass 2,
//!   `BamExtractor.cpp:732-748`); for single-end input, emitted directly
//!   with `usedName` dedup (`BamExtractor.cpp:749-778`).
//! - On-target aligned reads (`tag` monotonically advanced over sorted
//!   `genes`, `BamExtractor.cpp:805-822`): same candidates-vs-direct-emit
//!   split as above (`BamExtractor.cpp:824-850`).
//!
//! Pass 2 (paired input only) re-scans the BAM from the start,
//! filling in both mates of every `candidates` entry and emitting each pair
//! the moment BOTH mates have been seen (`BamExtractor.cpp:878-937`) --
//! output order is BAM re-encounter order of pair COMPLETION, not
//! `candidates` map-insertion order.

use crate::alignments::Alignments;
use crate::extract::{CandidateSink, HIT_LEN_REQUIRED_PAIRED, HIT_LEN_REQUIRED_SINGLE, ReadRecord};
use crate::ref_kmer_filter::{RefKmerFilter, Scratch, is_low_complexity};
use anyhow::{Context, Result, bail};
use rayon::prelude::*;
use std::collections::HashMap;

/// Batch size used by the parallel (`threads > 1`) pass-1 candidate-decision
/// path -- mirrors [`crate::extract::extract_candidates_with_threads`]'s own
/// batching knob and rationale: purely a throughput/memory bound, with no
/// effect on output (see [`evaluate_chunk`]'s doc comment for why the
/// decisions themselves are computed identically regardless of batch size).
const PARALLEL_BATCH_SIZE_PER_THREAD: usize = 512;

/// Max [`Pass1Site`]s scanned/evaluated/applied per chunk in
/// [`run_pass1_chunked`]'s scan-evaluate-apply loop, bounding pass-1 peak
/// memory to `O(chunk + candidates)` instead of `O(input)` (the whole point
/// of the chunk loop -- see that function's doc comment). Large enough to
/// amortize the parallel evaluate step (a multiple of
/// [`PARALLEL_BATCH_SIZE_PER_THREAD`]) while staying small relative to a WGS
/// BAM's read count.
const PASS1_CHUNK_SIZE: usize = PARALLEL_BATCH_SIZE_PER_THREAD * 64;

/// A single reference-coordinate gene interval, mirroring T1K's `_interval`
/// (`BamExtractor.cpp:49-63`): `chrId` plus an inclusive `[start, end]`
/// range. The `Ord`/`PartialOrd` derive matches `_interval::operator<`
/// EXACTLY: ascending `(chrId, start, end)` lexicographic order (the C++
/// operator compares `chrId` first, then `start`, then `end` -- precisely
/// Rust's derived tuple-field order for a 3-field struct).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct GeneInterval {
    pub chr_id: i32,
    pub start: i32,
    pub end: i32,
}

/// Parses a `_coord.fa`-style reference (e.g.
/// `fixtures/refbuild/golden/kir_rna_coord.fa`) into
/// `(name, chrom, start, end, strand)` tuples, one per FASTA record header,
/// mirroring the header-scan loop `BamExtractor.cpp:557-567`'s `fscanf( fpRef,
/// "%s %s %d %d %s", geneName, chrom, &start, &end, strand )` (five
/// whitespace-separated tokens on the header line, PAST the leading `>`) --
/// the record's SEQUENCE lines are read-but-discarded by this parse. Stock
/// reads the sequence with `fscanf( fpRef, "%s", buffer )` (`BamExtractor.cpp:566`)
/// -- exactly ONE whitespace-delimited token -- i.e. stock ASSUMES each
/// sequence is on a single unwrapped line and would desync (parse a wrapped
/// continuation line as the next header) on a multi-line sequence. This parse
/// is instead LINE-BASED: it accumulates every non-`>` line until the next
/// header, so it correctly handles multi-line-wrapped sequences that stock
/// cannot -- a deliberate DIVERGENCE in the safe direction. Every real
/// `_coord.fa` (machine-emitted by `unum build`) is single-line, so the two
/// agree byte-for-byte on all real inputs; the divergence only surfaces on
/// pathological wrapped input, where `unum` is the more robust of the two.
///
/// This function does NOT resolve `chrom` to a `chrId` -- that requires an
/// open [`Alignments`] (`GetChromIdFromName`), which the caller performs
/// separately (see [`build_genes`]), matching
/// `BamExtractor.cpp:560`'s `alignments.GetChromIdFromName( chrom )` call
/// inside the same loop (kept as a separate step here purely so this parse
/// step is independently unit-testable without an open BAM file).
///
/// # Errors
///
/// Returns an error if `path` cannot be read, or if any record's header line
/// does not have exactly five whitespace-separated tokens (name, chrom,
/// start, end, strand) or a non-integer start/end.
pub fn parse_coord_fa(path: &std::path::Path) -> Result<Vec<CoordRecord>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading coord FASTA {}", path.display()))?;

    let mut records = Vec::new();
    let mut current_header: Option<Vec<String>> = None;
    let mut current_seq = String::new();

    let flush =
        |header: &Option<Vec<String>>, seq: &str, records: &mut Vec<CoordRecord>| -> Result<()> {
            if let Some(tokens) = header {
                anyhow::ensure!(
                    tokens.len() == 5,
                    "coord FASTA header must have 5 whitespace-separated tokens (name chrom start \
                 end strand), got {}: {:?}",
                    tokens.len(),
                    tokens
                );
                let start: i32 = tokens[2].parse().with_context(|| {
                    format!("coord FASTA start not an integer: {:?}", tokens[2])
                })?;
                let end: i32 = tokens[3]
                    .parse()
                    .with_context(|| format!("coord FASTA end not an integer: {:?}", tokens[3]))?;
                records.push(CoordRecord {
                    name: tokens[0].clone(),
                    chrom: tokens[1].clone(),
                    start,
                    end,
                    strand: tokens[4].clone(),
                    seq: seq.to_string(),
                });
            }
            Ok(())
        };

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            flush(&current_header, &current_seq, &mut records)?;
            current_header = Some(rest.split_whitespace().map(str::to_string).collect());
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim_end());
        }
    }
    flush(&current_header, &current_seq, &mut records)?;

    Ok(records)
}

/// A single parsed `_coord.fa` record: the k-mer-reference sequence
/// (`name`/`seq`, fed to [`RefKmerFilter::from_reference_fasta`]-equivalent
/// loading) plus its genome interval (`chrom`/`start`/`end`/`strand`, fed to
/// [`build_genes`]). See [`parse_coord_fa`]'s doc comment for the exact
/// parse rules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoordRecord {
    pub name: String,
    pub chrom: String,
    pub start: i32,
    pub end: i32,
    pub strand: String,
    pub seq: String,
}

/// Resolves each [`CoordRecord`]'s `chrom` to a `chrId` via
/// `alignments.chrom_id_from_name` and returns the SORTED interval list,
/// mirroring `BamExtractor.cpp:559-570`: `ni.chrId =
/// alignments.GetChromIdFromName( chrom ); ... genes.push_back( ni ) ; ...
/// std::sort( genes.begin(), genes.end() ) ;` -- note `geneCnt = genes.size()`
/// is captured BEFORE the sort in stock, but since sorting never changes a
/// `Vec`'s length, capturing it after (via `.len()` on the returned `Vec`)
/// is equivalent; callers should use the returned `Vec`'s length as
/// `geneCnt`.
///
/// # Errors
///
/// Returns an error if any record's `chrom` cannot be resolved to a `chrId`
/// (mirrors `GetChromIdFromName`'s own `fprintf`+`exit(1)` on an unknown
/// name, surfaced here as a recoverable `Result`).
pub fn build_genes(alignments: &Alignments, records: &[CoordRecord]) -> Result<Vec<GeneInterval>> {
    let mut genes: Vec<GeneInterval> = records
        .iter()
        .map(|r| {
            let chr_id = alignments
                .chrom_id_from_name(&r.chrom)
                .with_context(|| format!("resolving coord FASTA chrom {:?}", r.chrom))?;
            Ok(GeneInterval { chr_id, start: r.start, end: r.end })
        })
        .collect::<Result<_>>()?;
    genes.sort_unstable();
    Ok(genes)
}

/// `BamExtractor.cpp:118-129`: `chrom` is a "valid alternative chromosome"
/// (i.e. treated the same as an off-target/alt-contig hit) if its name
/// contains `_`, `.`, OR `*` anywhere. EXACT reproduction, including the
/// dead/commented-out extra checks left in stock (not reproduced, since they
/// are unreachable in the C++ itself).
#[must_use]
pub fn valid_alternative_chrom(chrom: &str) -> bool {
    chrom.contains('_') || chrom.contains('.') || chrom.contains('*')
}

/// `BamExtractor.cpp:168-183`: trims a read name for mate-matching.
/// `trim_len == -1` (the default, matching `mateIdLen`'s stock default) means
/// "strip a trailing `/1` or `/2`" (last char `1` or `2` AND the char before
/// it is `/`); otherwise erases exactly the LAST `trim_len` characters,
/// unconditionally (mirrors `std::string::erase(len - trimLen, trimLen)`,
/// which -- like the C++ -- has no bounds check: a `trim_len` longer than
/// `name` panics here exactly where stock's `erase` would be undefined
/// behavior on an out-of-range position).
///
/// # Panics
///
/// Panics if `trim_len > 0` and `trim_len as usize > name.len()` (see above).
#[must_use]
pub fn trim_name(name: &str, trim_len: i32) -> String {
    if trim_len == -1 {
        let bytes = name.as_bytes();
        let len = bytes.len();
        if len >= 2 && (bytes[len - 1] == b'1' || bytes[len - 1] == b'2') && bytes[len - 2] == b'/'
        {
            return name[..len - 2].to_string();
        }
        return name.to_string();
    }

    let trim_len = usize::try_from(trim_len).unwrap_or_else(|_| {
        panic!("trim_name: negative trim_len {trim_len} (only -1 is a valid sentinel)")
    });
    let len = name.len();
    assert!(
        trim_len <= len,
        "trim_name: trim_len {trim_len} exceeds name length {len} for {name:?}"
    );
    name[..len - trim_len].to_string()
}

/// `BamExtractor.cpp:576-580`: computes `hitLenRequired` from the BAM's
/// sampled fragment/read statistics -- `21` for paired input, `17` for
/// single-end (`frag_stdev == 0`), then bumped up to `read_len / 5`
/// (INTEGER division) if that exceeds the base value.
#[must_use]
pub fn compute_hit_len_required(frag_stdev: i32, read_len: i32) -> i32 {
    let mut hit_len_required = if frag_stdev == 0 { 17 } else { 21 };
    if read_len / 5 > hit_len_required {
        hit_len_required = read_len / 5;
    }
    hit_len_required
}

/// Computes `--bam-mode no-alignment`'s `hitLenRequired` using the FASTQ
/// path's formula (`extract.rs`'s `extract_candidates_with_threads` setup,
/// ~lines 346-368: base [`HIT_LEN_REQUIRED_PAIRED`]/[`HIT_LEN_REQUIRED_SINGLE`]
/// (27/23), bumped to `sampled_len/(count*5)` when that exceeds the base) --
/// NOT [`compute_hit_len_required`] (the ALIGNMENT path's formula: 21/17,
/// `read_len/5`). Pinning to the FASTQ setup, byte-for-byte, is what makes
/// `no-alignment ≡ FASTQ` exact (an equivalence that holds only under uniform
/// read length, since `floor(sampled_len/(count*5))` is then
/// sampling-order-independent -- validated by adversarial review; see the
/// module docs and the `no_alignment ≡ FASTQ` equivalence test). Assumes
/// uniform read length: realistic sequencing runs are fixed-length, so this
/// assumption holds in practice even though it is not checked here.
///
/// `sampled_read1_len_sum` is the SUM (not mean) of sampled read-1 sequence
/// lengths -- the same quantity `extract.rs`'s `sample_head` returns and
/// [`crate::alignments::Alignments::sample_read1_len_sum`] computes over a
/// BAM. `sampled_count == 0` (an empty/no-read-1 sample) returns the base
/// value unchanged, guarding the division against a zero divisor.
///
/// Setup-only for now: no `no-alignment` driver wires this in yet (a later
/// task does), hence `#[allow(dead_code)]`.
#[allow(dead_code)]
fn compute_hit_len_required_no_alignment(
    sampled_read1_len_sum: i64,
    sampled_count: usize,
    has_mate: bool,
) -> i32 {
    let base = if has_mate { HIT_LEN_REQUIRED_PAIRED } else { HIT_LEN_REQUIRED_SINGLE };
    if sampled_count == 0 {
        return base;
    }
    let candidate = sampled_read1_len_sum / (i64::try_from(sampled_count).unwrap_or(i64::MAX) * 5);
    if candidate > i64::from(base) { i32::try_from(candidate).unwrap_or(i32::MAX) } else { base }
}

/// `BamExtractor.cpp:480`: the literal initial k-mer length the reference is
/// first loaded at, before `InferKmerLength`/`UpdateKmerLength`.
pub const INITIAL_KMER_LENGTH: usize = 9;

/// Extraction run summary, returned by [`extract_from_bam`] for caller
/// diagnostics/logging.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BamExtractMetrics {
    /// Whether the source was treated as single-end (`frag_stdev == 0`).
    pub single_end: bool,
    /// The final `hitLenRequired` the filter ended up using.
    pub hit_len_required: i32,
    /// The final `kmerLength` the filter ended up using.
    pub kmer_length: usize,
    /// Number of pass-1-direct-emitted pairs/reads (unaligned-template pairs
    /// for paired input; on-target/alt-chrom aligned OR unmapped reads for
    /// single-end input).
    pub pass1_emitted: u64,
    /// Number of paired candidates recorded in pass 1 (0 for single-end,
    /// which never populates `candidates`).
    pub candidates_recorded: u64,
    /// Number of paired candidates completed (both mates found) and emitted
    /// in pass 2 (0 for single-end).
    pub pass2_emitted: u64,
}

/// A single pass-1-recorded paired candidate awaiting its mate(s) in pass 2,
/// mirroring `_candidate` (`BamExtractor.cpp:65-72`): `mate1`/`mate2` start
/// `None` and are filled in independently as pass 2 encounters each mate.
#[derive(Debug, Default, Clone)]
struct PendingCandidate {
    mate1: Option<ReadRecord>,
    mate2: Option<ReadRecord>,
}

/// Runs the full two-pass `BamExtractor.cpp:main` BAM/CRAM read-extraction
/// driver (see module docs for the pass-1/pass-2 split), single-threaded
/// (`threadCnt == 1` semantics -- see module docs for why this is the ONLY
/// code path this port reproduces).
///
/// `filter` must already be loaded (via
/// [`RefKmerFilter::from_reference_fasta`]-equivalent, e.g. from the
/// `_coord.fa`'s sequences) at [`INITIAL_KMER_LENGTH`] (9). `genes` must
/// already be built (via [`build_genes`]) and sorted (which `build_genes`
/// does itself). `abnormal_unaligned_flag` corresponds to the `-u` CLI flag;
/// `mate_id_len` corresponds to `--mateIdSuffixLen` (default `-1`).
///
/// This function performs the data-dependent setup itself
/// (`compute_hit_len_required` + `infer_kmer_length`/`update_kmer_length`,
/// `BamExtractor.cpp:572-591`) from `alignments.general_info(true)`, so
/// `alignments` must be freshly opened (not yet advanced) when this is
/// called -- `general_info` consumes the reader's position, and this
/// function itself calls `alignments.rewind()` immediately after (mirroring
/// `BamExtractor.cpp:573-574`).
///
/// # Errors
///
/// Returns an error if: an unaligned-template pair's second record is
/// missing or has a mismatched trimmed name (mirrors
/// `BamExtractor.cpp:657-672`'s "Two reads from the unaligned fragment are
/// not showing up together" error); any underlying BAM read/rewind fails; or
/// [`build_genes`]'s chrom resolution fails (if not already done by the
/// caller).
pub fn extract_from_bam(
    alignments: &mut Alignments,
    filter: &mut RefKmerFilter,
    genes: &[GeneInterval],
    abnormal_unaligned_flag: bool,
    mate_id_len: i32,
    sink: &mut impl CandidateSink,
) -> Result<BamExtractMetrics> {
    extract_from_bam_with_threads(
        alignments,
        filter,
        genes,
        abnormal_unaligned_flag,
        mate_id_len,
        1,
        sink,
    )
}

/// Same as [`extract_from_bam`], but parallelizes pass 1's expensive per-read
/// candidate DECISION (`is_low_complexity` + `is_good_candidate_with_scratch`,
/// the only calls into the k-mer filter this whole two-pass driver makes --
/// pass 2 does none) across `threads` worker threads via a scoped `rayon`
/// thread pool, while keeping every OTHER aspect of pass 1 -- the `tag`
/// gene-interval-scan pointer advance, `used_name`/`candidates` map
/// mutation, and the sequential emit order -- exactly as sequential as
/// before. See [`run_pass1_with_threads`]'s doc comment for the
/// scan/evaluate/apply split this uses to make that possible.
///
/// `threads <= 1` takes the exact sequential fast path (no rayon pool is
/// built at all, and pass 1 runs as a single BAM-encounter-order loop,
/// identical to [`extract_from_bam`]'s prior behavior). `threads > 1`
/// produces byte-identical output to `threads <= 1` (module docs: this port,
/// unlike the oracle's own `threadCnt > 1` path, never changes batching/flush
/// timing for ANY candidate path -- pass 1's `apply` sub-phase always replays
/// every site's outcome in the exact BAM-encounter order the `scan`
/// sub-phase recorded it in, regardless of how many threads evaluated the
/// decisions).
///
/// # Errors
///
/// See [`extract_from_bam`].
pub fn extract_from_bam_with_threads(
    alignments: &mut Alignments,
    filter: &mut RefKmerFilter,
    genes: &[GeneInterval],
    abnormal_unaligned_flag: bool,
    mate_id_len: i32,
    threads: usize,
    sink: &mut impl CandidateSink,
) -> Result<BamExtractMetrics> {
    // Setup (BamExtractor.cpp:572-591).
    let general_info = alignments.general_info(true).context("computing general_info")?;
    alignments.rewind().context("rewinding after general_info")?;

    let mut hit_len_required =
        compute_hit_len_required(general_info.frag_stdev, general_info.read_len);
    filter.set_hit_len_required(hit_len_required);

    let inferred = filter.infer_kmer_length();
    if inferred > filter.kmer_length() {
        filter.update_kmer_length(inferred);
        if inferred > usize::try_from(hit_len_required).unwrap_or(0) {
            hit_len_required = i32::try_from(inferred).unwrap_or(i32::MAX);
            filter.set_hit_len_required(hit_len_required);
        }
    }

    let single_end = general_info.frag_stdev == 0;

    let (candidates, pass1_emitted) = run_pass1_with_threads(
        alignments,
        filter,
        genes,
        single_end,
        abnormal_unaligned_flag,
        mate_id_len,
        threads,
        sink,
        Selection::Alignment,
    )?;

    alignments.rewind().context("rewinding before pass 2 (or final close for single-end)")?;

    let candidates_recorded = u64::try_from(candidates.len()).unwrap_or(u64::MAX);
    if single_end {
        // BamExtractor.cpp:858-870: single-end terminates after pass 1.
        return Ok(BamExtractMetrics {
            single_end: true,
            hit_len_required,
            kmer_length: filter.kmer_length(),
            pass1_emitted,
            candidates_recorded: 0,
            pass2_emitted: 0,
        });
    }

    let pass2_emitted =
        run_pass2(alignments, candidates, abnormal_unaligned_flag, mate_id_len, sink)?;

    Ok(BamExtractMetrics {
        single_end: false,
        hit_len_required,
        kmer_length: filter.kmer_length(),
        pass1_emitted,
        candidates_recorded,
        pass2_emitted,
    })
}

/// Handles one unaligned-template pair within pass 1
/// (`BamExtractor.cpp:646-728`): the CURRENT record plus the NEXT record
/// (read here via a second `alignments.next()` call) are the two mates of an
/// unaligned template, which -- for a coordinate-sorted BAM without `-u`
/// (`abnormal_unaligned_flag`) -- must appear as two consecutive records.
/// Returns `true` if the pair was emitted (both mates pass the low-complexity
/// + candidate-filter gate), `false` if it was filtered out.
///
/// Only used by the `#[cfg(test)]`-only reference [`run_pass1`] now (see that
/// function's doc comment) -- [`run_pass1_with_threads`]'s production path
/// reproduces this exact logic inline, split across [`scan_pass1_chunk`] (record
/// reading) and [`evaluate_pass1_site`]'s [`Pass1Site::UnalignedPair`] arm
/// (the candidate decision).
///
/// # Errors
///
/// Returns an error if the second record is missing (EOF) or its trimmed
/// read id does not match the first record's (mirrors
/// `BamExtractor.cpp:657-672`'s "Two reads from the unaligned fragment are
/// not showing up together" error).
#[cfg(test)]
fn handle_unaligned_template_pair(
    alignments: &mut Alignments,
    filter: &mut RefKmerFilter,
    scratch: &mut Scratch,
    mate_id_len: i32,
    sink: &mut impl CandidateSink,
) -> Result<bool> {
    let seq1 = alignments.read_seq();
    let qual1 = alignments.qual();
    let name = alignments.read_id();

    if !alignments.next().context("pass 1: reading unaligned mate")? {
        bail!(
            "Two reads from the unaligned fragment are not showing up together. Please use \
             -u(--abnormalUnmapFlag from wrapper) option."
        );
    }
    let mate_name = alignments.read_id();
    let seq2 = alignments.read_seq();
    let qual2 = alignments.qual();
    // BamExtractor.cpp:681: `IsFirstMate()` is evaluated on the CURRENT
    // (second-read) record at this point in the C++, AFTER
    // `alignments.Next()` advanced to it -- not on the first record.
    let second_record_is_first_mate = alignments.is_first_mate();

    let trimmed_name = trim_name(&name, mate_id_len);
    let trimmed_mate_name = trim_name(&mate_name, mate_id_len);
    if trimmed_name != trimmed_mate_name {
        bail!(
            "Two reads from the unaligned fragment are not showing up together. Please use \
             -u(--abnormalUnmapFlag from wrapper) option."
        );
    }

    if (!is_low_complexity(&seq1) && !is_low_complexity(&seq2))
        && (filter.is_good_candidate_with_scratch(&seq1, scratch)
            || filter.is_good_candidate_with_scratch(&seq2, scratch))
    {
        let rec1 = ReadRecord { id: trimmed_name.clone(), seq: seq1, qual: Some(qual1) };
        let rec2 = ReadRecord { id: trimmed_name, seq: seq2, qual: Some(qual2) };
        // BamExtractor.cpp:681-690: order by the SECOND record's
        // IsFirstMate -- if the second record IS first-mate, it is mate1 and
        // the first record is mate2 (fp1=buffer, fp2=buffer2); else (second
        // record is mate2) the first record is mate1.
        if second_record_is_first_mate {
            sink.emit_pair(&rec2, Some(&rec1))?;
        } else {
            sink.emit_pair(&rec1, Some(&rec2))?;
        }
        return Ok(true);
    }
    Ok(false)
}

/// Which read-selection SEMANTICS pass 1 uses -- threaded through
/// [`scan_pass1_chunk`]/[`apply_pass1_chunk`]/[`run_pass1_chunked`]/
/// [`run_pass1_with_threads`] so both modes share the exact same bounded
/// scan/evaluate/apply chunk loop.
///
/// - `Alignment` is the EXISTING position/gene-interval-driven classification
///   (`BamExtractor.cpp:632-851`, unchanged by this enum's introduction --
///   [`extract_from_bam_with_threads`] passes this variant, so its output
///   stays byte-identical to before this enum existed).
/// - `NoAlignment` is the Class-A "BAM ~ FASTQ" selection: every PRIMARY read
///   is k-mer-tested on its own sequence, position/gene-interval/`tag`
///   entirely bypassed (see [`Pass1Site::KmerCandidate`]'s doc comment).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Selection {
    Alignment,
    NoAlignment,
}

/// One pass-1 site awaiting a candidate-filter DECISION, captured by
/// [`scan_pass1_chunk`] and resolved by [`evaluate_chunk`] -- see
/// [`run_pass1_with_threads`]'s doc comment for why splitting pass 1 into
/// scan/evaluate/apply sub-phases is what makes the expensive decision
/// (`is_low_complexity` + `is_good_candidate_with_scratch`) safely
/// parallelizable while every other aspect of pass 1 (the `tag`
/// gene-interval-scan pointer, `used_name`/`candidates` mutation, emit order)
/// stays exactly as sequential -- and therefore exactly as byte-identical to
/// the pre-existing single-threaded behavior -- as before.
///
/// Each variant carries only the raw sequence(s) [`evaluate_chunk`]
/// needs to test, plus whatever bookkeeping [`apply_pass1_chunk`] needs to
/// replay the ORIGINAL branch's side effect once the decision is known.
enum Pass1Site {
    /// `handle_unaligned_template_pair`'s candidate gate: both mates' raw
    /// sequences (evaluated as `!low_complexity(seq1) && !low_complexity(seq2)
    /// && (good(seq1) || good(seq2))`, matching that function's exact
    /// short-circuit-free `&&`/`||` combination -- see its doc comment),
    /// plus everything needed to reconstruct and emit the pair.
    UnalignedPair {
        seq1: Vec<u8>,
        qual1: Vec<u8>,
        seq2: Vec<u8>,
        qual2: Vec<u8>,
        trimmed_name: String,
        second_record_is_first_mate: bool,
    },
    /// The single-end not-template-aligned branch (`BamExtractor.cpp:749-778`,
    /// reached via `!template_aligned || chrom_alt`): evaluated as
    /// `!low_complexity(seq) && good(seq)`.
    SingleEndNotAligned { seq: Vec<u8>, qual: Vec<u8>, is_aligned: bool, name: String },
    /// The paired not-template-aligned/alt-chrom branch
    /// (`BamExtractor.cpp:732-748`): evaluated as `!low_complexity(seq) &&
    /// good(seq)`. Only the trimmed name is needed (candidates map key; no
    /// sequence/qual is stored here -- pass 2 re-reads it).
    PairedNotAlignedCandidate { seq: Vec<u8>, trimmed_name: String },
    /// The single-end on-target-aligned branch (`BamExtractor.cpp:824-850`
    /// single-end half): evaluated as `good(seq)`. `is_low_complexity` has
    /// ALREADY been checked by [`scan_pass1_chunk`] before this site is even
    /// created (`BamExtractor.cpp:812-814` runs that check unconditionally,
    /// before the single-end/paired split), so
    /// [`evaluate_pass1_site`] does not re-check it for this variant.
    SingleEndOnTarget { seq: Vec<u8>, qual: Vec<u8>, name: String },
    /// The paired on-target-aligned branch (`BamExtractor.cpp:824-850` paired
    /// half): resolves unconditionally to `true` (`is_low_complexity` already
    /// checked by [`scan_pass1_chunk`], same as [`Pass1Site::SingleEndOnTarget`]; no
    /// candidate-filter call on this path -- see [`evaluate_pass1_site`]'s
    /// doc comment on this variant). Only the trimmed name is needed
    /// (candidates map key; pass 2 re-reads the sequence).
    PairedOnTargetCandidate { trimmed_name: String },
    /// [`Selection::NoAlignment`]'s Class-A site: one PRIMARY read's raw
    /// sequence, tested REGARDLESS of alignment position/gene interval --
    /// evaluated as `!low_complexity(seq) && good(seq)`, the same gate as
    /// [`Pass1Site::PairedNotAlignedCandidate`]/[`Pass1Site::SingleEndNotAligned`].
    /// `qual` is only populated (and only needed) for single-end input, which
    /// emits a lone record directly on a pass; paired input only needs
    /// `trimmed_name` (the 2-pass name-map key -- pass 2 re-reads seq/qual),
    /// so it leaves `qual` empty. `is_first_mate` mirrors the record's own
    /// `0x40` flag for paired input (`true` for single-end, where mate order
    /// is moot) -- not consumed by [`apply_pass1_chunk`] today, carried for a
    /// future pass-2 consumer.
    KmerCandidate { seq: Vec<u8>, qual: Vec<u8>, trimmed_name: String, is_first_mate: bool },
}

/// Scans the BAM sequentially (BAM-encounter order, exactly as
/// [`run_pass1`]'s original single loop did), performing every CHEAP,
/// order-dependent decision itself (the `tag` gene-interval-scan pointer
/// advance, `template_aligned`/`chrom_alt`/`is_aligned` classification,
/// `is_low_complexity` pre-checks that gate whether a site is even created)
/// but DEFERRING every EXPENSIVE `is_good_candidate_with_scratch` call into a
/// [`Pass1Site`] pushed onto the returned `Vec`, in the exact order
/// encountered. This is the only sub-phase that touches `alignments`
/// (`Alignments::next()` is a stateful cursor -- see that method's doc
/// comment -- so record reading itself can never be parallelized; only the
/// DECISION on an already-read sequence can).
///
/// Stops once `limit` sites have been pushed (or at EOF), returning the
/// (possibly short/empty) chunk -- an empty return means EOF. `tag` (the
/// monotonic gene-interval-scan pointer) is threaded in/out via `&mut` so its
/// advance persists across chunks exactly as it would across one
/// unbounded scan: the caller owns `tag`'s storage and initializes it to `0`
/// once, before the first chunk. The `sites.len() < limit` check sits at the
/// LOOP TOP (not inside the unaligned-pair branch), so a chunk boundary can
/// never land between an unaligned pair's two records -- that branch reads
/// TWO records in one iteration but pushes exactly ONE site.
///
/// Errors (missing/mismatched unaligned mate) are raised here, at the exact
/// same point in the BAM scan stock would raise them at -- unaffected by
/// `threads` or `limit`, since no candidate-filter decision is needed to
/// detect them.
///
/// `selection` branches at the TOP of the loop: [`Selection::NoAlignment`]
/// skips non-primary reads (`!is_primary()`) and pushes every remaining
/// (primary) read as a [`Pass1Site::KmerCandidate`], entirely bypassing the
/// `tag`/gene-interval/`template_aligned` classification below --
/// [`Selection::Alignment`] takes the existing position-based branch,
/// unchanged.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn scan_pass1_chunk(
    alignments: &mut Alignments,
    genes: &[GeneInterval],
    single_end: bool,
    abnormal_unaligned_flag: bool,
    mate_id_len: i32,
    tag: &mut usize,
    limit: usize,
    selection: Selection,
) -> Result<Vec<Pass1Site>> {
    let gene_cnt = genes.len();
    let mut sites: Vec<Pass1Site> = Vec::with_capacity(limit.min(1024));

    while sites.len() < limit && alignments.next().context("pass 1: reading next BAM record")? {
        if selection == Selection::NoAlignment {
            // Class A: k-mer-test every PRIMARY read; ignore alignment
            // position entirely. Unlike the alignment path below (which does
            // NOT filter `0x900` in pass 1), non-primary (secondary/
            // supplementary) reads MUST be skipped here: a k-mer match on a
            // secondary alignment's (identical) sequence would double-count
            // the QNAME, and `no-alignment == FASTQ` requires each template's
            // primary read to be tested exactly once, matching the FASTQ
            // path's one-record-per-read semantics and pass 2's
            // primary-only fetch.
            if !alignments.is_primary() {
                continue;
            }
            let seq = alignments.read_seq();
            if single_end {
                let qual = alignments.qual();
                let name = alignments.read_id();
                sites.push(Pass1Site::KmerCandidate {
                    seq,
                    qual,
                    trimmed_name: name,
                    is_first_mate: true,
                });
            } else {
                let trimmed_name = trim_name(&alignments.read_id(), mate_id_len);
                sites.push(Pass1Site::KmerCandidate {
                    seq,
                    qual: Vec::new(),
                    trimmed_name,
                    is_first_mate: alignments.is_first_mate(),
                });
            }
            continue;
        }

        let template_aligned = alignments.is_template_aligned();
        let chrom_alt = alignments.is_aligned()
            && valid_alternative_chrom(&alignments.chrom_name(alignments.chrom_id()));

        if !template_aligned || chrom_alt {
            if !template_aligned && !single_end && !abnormal_unaligned_flag {
                let seq1 = alignments.read_seq();
                let qual1 = alignments.qual();
                let name = alignments.read_id();

                if !alignments.next().context("pass 1: reading unaligned mate")? {
                    bail!(
                        "Two reads from the unaligned fragment are not showing up together. \
                         Please use -u(--abnormalUnmapFlag from wrapper) option."
                    );
                }
                let mate_name = alignments.read_id();
                let seq2 = alignments.read_seq();
                let qual2 = alignments.qual();
                // BamExtractor.cpp:681: `IsFirstMate()` on the CURRENT
                // (second-read) record, after `Next()` advanced to it.
                let second_record_is_first_mate = alignments.is_first_mate();

                let trimmed_name = trim_name(&name, mate_id_len);
                let trimmed_mate_name = trim_name(&mate_name, mate_id_len);
                if trimmed_name != trimmed_mate_name {
                    bail!(
                        "Two reads from the unaligned fragment are not showing up together. \
                         Please use -u(--abnormalUnmapFlag from wrapper) option."
                    );
                }

                sites.push(Pass1Site::UnalignedPair {
                    seq1,
                    qual1,
                    seq2,
                    qual2,
                    trimmed_name,
                    second_record_is_first_mate,
                });
                continue;
            }

            if single_end {
                let seq = alignments.read_seq();
                let qual = alignments.qual();
                let is_aligned = alignments.is_aligned();
                let name = alignments.read_id();
                sites.push(Pass1Site::SingleEndNotAligned { seq, qual, is_aligned, name });
            } else {
                let seq = alignments.read_seq();
                let trimmed_name = trim_name(&alignments.read_id(), mate_id_len);
                sites.push(Pass1Site::PairedNotAlignedCandidate { seq, trimmed_name });
            }
            continue;
        }

        if !alignments.is_aligned() {
            // The unaligned mate of an aligned pair (BamExtractor.cpp:801-802).
            continue;
        }

        // Aligned reads reaching here (BamExtractor.cpp:804-850).
        let chr_id = alignments.chrom_id();
        let segments = alignments.segments();
        let Some(first_seg) = segments.first() else { continue };
        let start = first_seg.a;
        #[allow(clippy::cast_possible_truncation)]
        let start_i32 = start as i32;
        let end = segments[segments.len() - 1].b;
        #[allow(clippy::cast_possible_truncation)]
        let end_i32 = end as i32;

        while *tag < gene_cnt
            && (chr_id > genes[*tag].chr_id
                || (chr_id == genes[*tag].chr_id && start_i32 > genes[*tag].end))
        {
            *tag += 1;
        }

        if *tag >= gene_cnt {
            continue;
        }
        if chr_id < genes[*tag].chr_id
            || (chr_id == genes[*tag].chr_id && end_i32 <= genes[*tag].start)
        {
            continue;
        }

        let seq = alignments.read_seq();
        if is_low_complexity(&seq) {
            continue;
        }

        if single_end {
            let name = alignments.read_id();
            let qual = alignments.qual();
            sites.push(Pass1Site::SingleEndOnTarget { seq, qual, name });
        } else {
            let trimmed_name = trim_name(&alignments.read_id(), mate_id_len);
            sites.push(Pass1Site::PairedOnTargetCandidate { trimmed_name });
        }
    }

    Ok(sites)
}

/// Resolves every [`Pass1Site`]'s candidate-filter decision for one chunk,
/// either sequentially (`pool = None`, no rayon pool involved -- direct
/// index-order iteration) or across the given pre-built `rayon` pool's worker
/// threads (`Some(pool)`, each worker reusing its own [`Scratch`] via
/// `map_init`, batched [`PARALLEL_BATCH_SIZE_PER_THREAD`]`* threads` sites at
/// a time). The pool is built ONCE by the caller ([`run_pass1_chunked`]) and
/// reused across every chunk, not rebuilt per chunk. Returns a `Vec<bool>`
/// PARALLEL to `sites` (same length, same index order) -- `rayon` preserves
/// input order on `collect()` for an `IndexedParallelIterator`
/// (`Vec::par_iter()`'s `.map_init()` is one), and the sequential fallback is
/// index-order by construction, so both produce IDENTICAL `Vec<bool>`
/// contents regardless of `pool` -- the decision function itself
/// ([`Pass1Site`]'s per-variant boolean expression, documented on each
/// variant) has no cross-site dependency, so evaluating sites in any order
/// (or concurrently) cannot change any individual site's answer.
fn evaluate_chunk(
    pool: Option<&rayon::ThreadPool>,
    filter: &RefKmerFilter,
    sites: &[Pass1Site],
) -> Vec<bool> {
    match pool {
        None => {
            let mut scratch = Scratch::default();
            sites.iter().map(|site| evaluate_pass1_site(filter, site, &mut scratch)).collect()
        }
        Some(pool) => pool.install(|| {
            sites
                .par_iter()
                .with_min_len(PARALLEL_BATCH_SIZE_PER_THREAD.max(1))
                .map_init(Scratch::default, |scratch, site| {
                    evaluate_pass1_site(filter, site, scratch)
                })
                .collect()
        }),
    }
}

/// The per-[`Pass1Site`] boolean decision -- see each variant's doc comment
/// for the exact expression this reproduces from the original inline
/// `run_pass1` branches.
fn evaluate_pass1_site(filter: &RefKmerFilter, site: &Pass1Site, scratch: &mut Scratch) -> bool {
    match site {
        Pass1Site::UnalignedPair { seq1, seq2, .. } => {
            (!is_low_complexity(seq1) && !is_low_complexity(seq2))
                && (filter.is_good_candidate_with_scratch(seq1, scratch)
                    || filter.is_good_candidate_with_scratch(seq2, scratch))
        }
        Pass1Site::SingleEndNotAligned { seq, .. }
        | Pass1Site::PairedNotAlignedCandidate { seq, .. }
        | Pass1Site::KmerCandidate { seq, .. } => {
            !is_low_complexity(seq) && filter.is_good_candidate_with_scratch(seq, scratch)
        }
        // On-target sites resolve unconditionally to `true`: `BamExtractor.cpp:804-851`
        // emits/records on-target reads after ONLY the `IsLowComplexity` check (already
        // applied by `scan_pass1_chunk` before this site was created,
        // `BamExtractor.cpp:812-814`) -- it never calls `IsGoodCandidate` on the
        // on-target path. The untouched reference `run_pass1` (see its on-target arm)
        // reproduces this faithfully: no candidate-filter call. A prior version of this
        // function incorrectly gated on-target sites behind
        // `is_good_candidate_with_scratch`, which dropped on-target reads that overlap a
        // gene by coordinate but fail the kmer candidate filter (e.g. soft-clipped or
        // divergent-allele reads) -- a byte-identity regression vs. the oracle. Do not
        // reintroduce a candidate-filter call here; `scratch`/`filter` are accepted only
        // to keep this function's signature uniform across all `Pass1Site` variants.
        Pass1Site::SingleEndOnTarget { .. } | Pass1Site::PairedOnTargetCandidate { .. } => {
            let _ = (filter, scratch);
            true
        }
    }
}

/// Replays one chunk's [`Pass1Site`] pre-computed `decisions[i]` outcomes,
/// SEQUENTIALLY in the exact scan order [`scan_pass1_chunk`] recorded them in
/// -- reproducing `run_pass1`'s original mutation/emit side effects (`tag`
/// advance is NOT replayed here, since it was already fully resolved during
/// scanning and does not affect this sub-phase; `used_name`/`candidates`/
/// `pass1_emitted` ARE replayed here, in order, accumulating into the
/// caller-owned storage passed by `&mut` -- threaded across chunks by
/// [`run_pass1_chunked`] rather than recreated per chunk -- so their final
/// contents/emit order are identical to the original single-loop
/// implementation regardless of `threads` or chunk boundaries).
///
/// `selection` is accepted for signature uniformity across the scan/
/// evaluate/apply split (mirroring [`evaluate_pass1_site`]'s on-target arm,
/// which accepts `filter`/`scratch` it doesn't use for the same reason): the
/// [`Pass1Site::KmerCandidate`] arm below dispatches on `single_end` alone
/// (its site shape already encodes which [`Selection`] produced it), so this
/// function does not itself need to branch on `selection`.
#[allow(clippy::too_many_arguments)]
fn apply_pass1_chunk(
    sites: Vec<Pass1Site>,
    decisions: &[bool],
    single_end: bool,
    candidates: &mut HashMap<String, PendingCandidate>,
    used_name: &mut std::collections::HashSet<String>,
    pass1_emitted: &mut u64,
    sink: &mut impl CandidateSink,
    selection: Selection,
) -> Result<()> {
    let _ = selection;
    for (site, &good) in sites.into_iter().zip(decisions) {
        match site {
            Pass1Site::UnalignedPair {
                seq1,
                qual1,
                seq2,
                qual2,
                trimmed_name,
                second_record_is_first_mate,
            } => {
                if good {
                    let rec1 =
                        ReadRecord { id: trimmed_name.clone(), seq: seq1, qual: Some(qual1) };
                    let rec2 = ReadRecord { id: trimmed_name, seq: seq2, qual: Some(qual2) };
                    // BamExtractor.cpp:681-690: see Pass1Site::UnalignedPair's
                    // doc comment for the mate-order rationale.
                    if second_record_is_first_mate {
                        sink.emit_pair(&rec2, Some(&rec1))?;
                    } else {
                        sink.emit_pair(&rec1, Some(&rec2))?;
                    }
                    *pass1_emitted += 1;
                }
            }
            Pass1Site::SingleEndNotAligned { seq, qual, is_aligned, name } => {
                debug_assert!(single_end, "SingleEndNotAligned site produced for paired input");
                if is_aligned && used_name.contains(&name) {
                    continue;
                }
                if good {
                    if is_aligned {
                        used_name.insert(name.clone());
                    }
                    let rec = ReadRecord { id: name, seq, qual: Some(qual) };
                    sink.emit_pair(&rec, None)?;
                    *pass1_emitted += 1;
                }
            }
            Pass1Site::PairedNotAlignedCandidate { trimmed_name, .. } => {
                debug_assert!(
                    !single_end,
                    "PairedNotAlignedCandidate site produced for single-end input"
                );
                if good {
                    candidates.entry(trimmed_name).or_default();
                }
            }
            Pass1Site::SingleEndOnTarget { seq, qual, name } => {
                debug_assert!(single_end, "SingleEndOnTarget site produced for paired input");
                if used_name.contains(&name) {
                    continue;
                }
                if good {
                    used_name.insert(name.clone());
                    let rec = ReadRecord { id: name, seq, qual: Some(qual) };
                    sink.emit_pair(&rec, None)?;
                    *pass1_emitted += 1;
                }
            }
            Pass1Site::PairedOnTargetCandidate { trimmed_name, .. } => {
                debug_assert!(
                    !single_end,
                    "PairedOnTargetCandidate site produced for single-end input"
                );
                if good {
                    candidates.entry(trimmed_name).or_default();
                }
            }
            Pass1Site::KmerCandidate { seq, qual, trimmed_name, is_first_mate } => {
                if single_end {
                    if good {
                        let rec = ReadRecord { id: trimmed_name, seq, qual: Some(qual) };
                        sink.emit_pair(&rec, None)?;
                        *pass1_emitted += 1;
                    }
                } else {
                    // Paired input only needs the name for the 2-pass
                    // name-map -- pass 2 re-reads seq/qual; `is_first_mate` is
                    // carried on the site for a future pass-2 consumer, not
                    // needed here.
                    let _ = (qual, is_first_mate);
                    if good {
                        candidates.entry(trimmed_name).or_default();
                    }
                }
            }
        }
    }

    Ok(())
}

/// Pass 1, dispatching on `threads`, as a bounded scan -> evaluate -> apply
/// CHUNK LOOP: builds the `rayon` thread pool ONCE (`threads > 1`; `None`
/// otherwise, `evaluate_chunk`'s sequential fallback), then repeatedly (a)
/// scans up to `chunk_size` sites from wherever [`scan_pass1_chunk`] left off
/// (`Alignments` is a stateful cursor, and `tag`, the gene-interval-scan
/// pointer, is threaded across iterations via `&mut` so its monotonic advance
/// is identical to one unbounded scan), (b) resolves each chunk's candidate
/// decisions ([`evaluate_chunk`], parallel when `threads > 1`), and (c)
/// replays that chunk's outcomes in scan order ([`apply_pass1_chunk`]),
/// accumulating into `candidates`/`used_name`/`pass1_emitted` across chunks
/// rather than per-chunk-fresh maps -- until a chunk comes back empty (EOF).
/// This bounds pass-1 peak memory to `O(chunk + candidates)` instead of
/// `O(input)` (an unbounded `Vec<Pass1Site>` over the whole BAM), while
/// remaining byte-identical to a single unbounded scan/evaluate/apply pass at
/// any `threads`/`chunk_size` -- chunk boundaries only ever fall between
/// sites (see [`scan_pass1_chunk`]'s doc comment on why an unaligned pair can
/// never be split), so scan order, decisions, and apply-time side effects are
/// all unaffected by where those boundaries land. Returns the `candidates`
/// map recorded for pass 2 (empty for single-end input, which never
/// populates it) and the number of pairs/reads emitted directly during this
/// pass -- both IDENTICAL to [`run_pass1`]'s output for the same input, at
/// any `threads`/`chunk_size`.
#[allow(clippy::too_many_arguments)]
fn run_pass1_chunked(
    alignments: &mut Alignments,
    filter: &RefKmerFilter,
    genes: &[GeneInterval],
    single_end: bool,
    abnormal_unaligned_flag: bool,
    mate_id_len: i32,
    threads: usize,
    chunk_size: usize,
    sink: &mut impl CandidateSink,
    selection: Selection,
) -> Result<(HashMap<String, PendingCandidate>, u64)> {
    // Build the rayon pool ONCE (not per chunk).
    let pool = if threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .context("building rayon thread pool for parallel BAM pass-1 evaluation")?,
        )
    } else {
        None
    };

    let mut candidates: HashMap<String, PendingCandidate> = HashMap::new();
    let mut used_name: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut pass1_emitted: u64 = 0;
    let mut tag: usize = 0;

    loop {
        let sites = scan_pass1_chunk(
            alignments,
            genes,
            single_end,
            abnormal_unaligned_flag,
            mate_id_len,
            &mut tag,
            chunk_size,
            selection,
        )?;
        if sites.is_empty() {
            break; // EOF
        }
        let decisions = evaluate_chunk(pool.as_ref(), filter, &sites);
        apply_pass1_chunk(
            sites,
            &decisions,
            single_end,
            &mut candidates,
            &mut used_name,
            &mut pass1_emitted,
            sink,
            selection,
        )?;
    }

    Ok((candidates, pass1_emitted))
}

/// Pass 1, dispatching on `threads`: [`run_pass1_chunked`] with the
/// production [`PASS1_CHUNK_SIZE`]. See that function's doc comment for the
/// bounded scan/evaluate/apply chunk-loop structure and its byte-identity
/// argument. Tests call [`run_pass1_chunked`] directly with a small
/// `chunk_size` to exercise cross-chunk-boundary correctness cheaply, without
/// needing a `PASS1_CHUNK_SIZE`-sized fixture.
#[allow(clippy::too_many_arguments)]
fn run_pass1_with_threads(
    alignments: &mut Alignments,
    filter: &RefKmerFilter,
    genes: &[GeneInterval],
    single_end: bool,
    abnormal_unaligned_flag: bool,
    mate_id_len: i32,
    threads: usize,
    sink: &mut impl CandidateSink,
    selection: Selection,
) -> Result<(HashMap<String, PendingCandidate>, u64)> {
    run_pass1_chunked(
        alignments,
        filter,
        genes,
        single_end,
        abnormal_unaligned_flag,
        mate_id_len,
        threads,
        PASS1_CHUNK_SIZE,
        sink,
        selection,
    )
}

/// Pass 1 (`BamExtractor.cpp:632-851`): the ORIGINAL single-threaded,
/// single-loop implementation, kept unmodified (byte-for-byte, not merely
/// behaviorally) as [`run_pass1_with_threads`]'s `threads <= 1` reference
/// semantics and as a standalone regression fixture -- if
/// [`run_pass1_with_threads`]'s scan/evaluate/apply split ever drifts from
/// this, a differential between the two (same inputs, `threads = 1` on the
/// new path) would catch it immediately. See [`extract_from_bam`]'s and this
/// module's doc comments for the full branch structure. Returns the
/// `candidates` map recorded for pass 2 (empty for single-end input, which
/// never populates it) and the number of pairs/reads emitted directly during
/// this pass.
#[cfg(test)]
fn run_pass1(
    alignments: &mut Alignments,
    filter: &mut RefKmerFilter,
    genes: &[GeneInterval],
    single_end: bool,
    abnormal_unaligned_flag: bool,
    mate_id_len: i32,
    sink: &mut impl CandidateSink,
) -> Result<(HashMap<String, PendingCandidate>, u64)> {
    let gene_cnt = genes.len();
    let mut scratch = Scratch::default();
    let mut candidates: HashMap<String, PendingCandidate> = HashMap::new();
    let mut used_name: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut tag: usize = 0;
    let mut pass1_emitted: u64 = 0;

    while alignments.next().context("pass 1: reading next BAM record")? {
        let template_aligned = alignments.is_template_aligned();
        let chrom_alt = alignments.is_aligned()
            && valid_alternative_chrom(&alignments.chrom_name(alignments.chrom_id()));

        if !template_aligned || chrom_alt {
            if !template_aligned && !single_end && !abnormal_unaligned_flag {
                if handle_unaligned_template_pair(
                    alignments,
                    filter,
                    &mut scratch,
                    mate_id_len,
                    sink,
                )? {
                    pass1_emitted += 1;
                }
                continue;
            }

            if single_end {
                // Single-end (BamExtractor.cpp:749-778). Stock's
                // `threadCnt == 1 || alignments.IsAligned()` gate is always
                // true here since this port is threadCnt==1-only.
                let seq = alignments.read_seq();
                let qual = alignments.qual();
                let is_aligned = alignments.is_aligned();
                let name = alignments.read_id();
                if is_aligned && used_name.contains(&name) {
                    continue;
                }
                if !is_low_complexity(&seq)
                    && filter.is_good_candidate_with_scratch(&seq, &mut scratch)
                {
                    if is_aligned {
                        used_name.insert(name.clone());
                    }
                    let rec = ReadRecord { id: name, seq, qual: Some(qual) };
                    sink.emit_pair(&rec, None)?;
                    pass1_emitted += 1;
                }
            } else {
                // Paired, alt-chrom or otherwise-not-template-aligned-but-not-
                // the-unaligned-pair-case (BamExtractor.cpp:732-748).
                let seq = alignments.read_seq();
                if !is_low_complexity(&seq)
                    && filter.is_good_candidate_with_scratch(&seq, &mut scratch)
                {
                    let name = trim_name(&alignments.read_id(), mate_id_len);
                    candidates.entry(name).or_default();
                }
            }
            continue;
        }

        if !alignments.is_aligned() {
            // The unaligned mate of an aligned pair (BamExtractor.cpp:801-802).
            continue;
        }

        // Aligned reads reaching here (BamExtractor.cpp:804-850).
        let chr_id = alignments.chrom_id();
        let segments = alignments.segments();
        let Some(first_seg) = segments.first() else { continue };
        let start = first_seg.a;
        #[allow(clippy::cast_possible_truncation)]
        let start_i32 = start as i32;
        let end = segments[segments.len() - 1].b;
        #[allow(clippy::cast_possible_truncation)]
        let end_i32 = end as i32;

        while tag < gene_cnt
            && (chr_id > genes[tag].chr_id
                || (chr_id == genes[tag].chr_id && start_i32 > genes[tag].end))
        {
            tag += 1;
        }

        if tag >= gene_cnt {
            continue;
        }
        if chr_id < genes[tag].chr_id
            || (chr_id == genes[tag].chr_id && end_i32 <= genes[tag].start)
        {
            continue;
        }

        let seq = alignments.read_seq();
        if is_low_complexity(&seq) {
            continue;
        }

        if single_end {
            let name = alignments.read_id();
            if used_name.contains(&name) {
                continue;
            }
            used_name.insert(name.clone());
            let qual = alignments.qual();
            let rec = ReadRecord { id: name, seq, qual: Some(qual) };
            sink.emit_pair(&rec, None)?;
            pass1_emitted += 1;
        } else {
            let name = trim_name(&alignments.read_id(), mate_id_len);
            candidates.entry(name).or_default();
        }
    }

    Ok((candidates, pass1_emitted))
}

/// Pass 2 (`BamExtractor.cpp:878-937`, paired input only): re-scans the BAM
/// from the start (caller must have already `rewind()`-ed `alignments`),
/// filling in both mates of every `candidates` entry and emitting each pair
/// the moment BOTH mates have been seen. Returns the number of pairs
/// emitted.
fn run_pass2(
    alignments: &mut Alignments,
    mut candidates: HashMap<String, PendingCandidate>,
    abnormal_unaligned_flag: bool,
    mate_id_len: i32,
    sink: &mut impl CandidateSink,
) -> Result<u64> {
    let candidate_cnt = candidates.len();
    let mut output_cnt: u64 = 0;

    if candidate_cnt == 0 {
        return Ok(0);
    }

    while alignments.next().context("pass 2: reading next BAM record")? {
        if !alignments.is_primary() {
            continue;
        }
        if !alignments.is_template_aligned() && !abnormal_unaligned_flag {
            continue;
        }

        let name = trim_name(&alignments.read_id(), mate_id_len);
        let Some(entry) = candidates.get_mut(&name) else {
            continue;
        };

        let seq = alignments.read_seq();
        let qual = alignments.qual();
        if alignments.is_first_mate() {
            entry.mate1 = Some(ReadRecord { id: name.clone(), seq, qual: Some(qual) });
        } else {
            entry.mate2 = Some(ReadRecord { id: name.clone(), seq, qual: Some(qual) });
        }

        // BamExtractor.cpp:917: `it->second.mate1 != NULL &&
        // it->second.mate2 != NULL` -- only emit once BOTH mates have been
        // filled (across this call and any earlier one for the same
        // candidate); otherwise leave the entry as-is, still waiting on the
        // other mate.
        if let (Some(m1_ref), Some(m2_ref)) = (entry.mate1.as_ref(), entry.mate2.as_ref()) {
            sink.emit_pair(m1_ref, Some(m2_ref))?;
            entry.mate1 = None;
            entry.mate2 = None;
            output_cnt += 1;
            if output_cnt == u64::try_from(candidate_cnt).unwrap_or(u64::MAX) {
                break;
            }
        }
    }

    Ok(output_cnt)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// On-target padding reads added to the parallel-regression fixtures so
    /// the `threads > 1` path actually splits work. `evaluate_chunk`
    /// uses `with_min_len(PARALLEL_BATCH_SIZE_PER_THREAD)`, so `rayon` refuses
    /// to split a `par_iter` shorter than that -- the fixtures must carry
    /// MORE than `PARALLEL_BATCH_SIZE_PER_THREAD` sites for the parallel case
    /// to span more than one chunk. Kept comfortably under ~930 so every
    /// padded read's coordinates (`base_pos + i*3`, read length 80) stay
    /// inside the 1000..5000 gene interval and remain on-target.
    // `PARALLEL_BATCH_SIZE_PER_THREAD` (512) is a compile-time constant that
    // trivially fits in `u32`, so this cast can never truncate.
    #[allow(clippy::cast_possible_truncation)]
    const PARALLEL_PAD_COUNT: u32 = PARALLEL_BATCH_SIZE_PER_THREAD as u32 + 88; // 600

    #[test]
    fn valid_alternative_chrom_matches_underscore_dot_star() {
        assert!(valid_alternative_chrom("chr19_KI270938v1_alt"));
        assert!(valid_alternative_chrom("HLA-A*01:01"));
        assert!(valid_alternative_chrom("some.name"));
        assert!(!valid_alternative_chrom("chr19"));
        assert!(!valid_alternative_chrom("chrX"));
    }

    #[test]
    fn trim_name_default_strips_trailing_slash_1_or_2() {
        assert_eq!(trim_name("read001/1", -1), "read001");
        assert_eq!(trim_name("read001/2", -1), "read001");
        assert_eq!(trim_name("read001", -1), "read001");
        // Not a trailing /1 or /2: no-op.
        assert_eq!(trim_name("read001/3", -1), "read001/3");
        assert_eq!(trim_name("read0012", -1), "read0012");
    }

    #[test]
    fn trim_name_explicit_len_erases_last_n_chars_unconditionally() {
        assert_eq!(trim_name("read001.foo", 4), "read001");
        assert_eq!(trim_name("read001", 0), "read001");
    }

    #[test]
    #[should_panic(expected = "exceeds name length")]
    fn trim_name_explicit_len_longer_than_name_panics() {
        let _ = trim_name("ab", 5);
    }

    #[test]
    fn compute_hit_len_required_paired_base_21() {
        // frag_stdev != 0 (paired), short reads: base 21, read_len/5 doesn't
        // exceed it.
        assert_eq!(compute_hit_len_required(50, 100), 21);
    }

    #[test]
    fn compute_hit_len_required_single_end_base_17() {
        assert_eq!(compute_hit_len_required(0, 50), 17);
    }

    #[test]
    fn compute_hit_len_required_bumped_by_read_len_over_5() {
        // Paired base 21, but read_len=150 -> 150/5=30 > 21 -> bumped to 30.
        assert_eq!(compute_hit_len_required(50, 150), 30);
        // Single-end base 17, read_len=150 -> 150/5=30 > 17 -> bumped to 30.
        assert_eq!(compute_hit_len_required(0, 150), 30);
    }

    #[test]
    fn compute_hit_len_required_integer_division_does_not_round_up() {
        // read_len=104 -> 104/5=20 (integer division), which does NOT exceed
        // the paired base of 21, so it stays at 21.
        assert_eq!(compute_hit_len_required(50, 104), 21);
        // read_len=109 -> 109/5=21, exactly equal (not strictly greater), so
        // no bump.
        assert_eq!(compute_hit_len_required(50, 109), 21);
        // read_len=110 -> 110/5=22 > 21, bumped.
        assert_eq!(compute_hit_len_required(50, 110), 22);
    }

    #[test]
    fn no_alignment_hit_len_matches_fastq_formula() {
        // Uniform 150bp reads, paired: base 27, bumped to 150/5 = 30.
        assert_eq!(compute_hit_len_required_no_alignment(150 * 1000, 1000, true), 30);
        // Uniform 100bp reads, paired: 100/5 = 20 < 27 base -> stays 27.
        assert_eq!(compute_hit_len_required_no_alignment(100 * 1000, 1000, true), 27);
        // Single-end 90bp: base 23, 90/5 = 18 < 23 -> 23.
        assert_eq!(compute_hit_len_required_no_alignment(90 * 500, 500, false), 23);
    }

    #[test]
    fn gene_interval_ord_matches_chr_id_then_start_then_end() {
        let mut genes = vec![
            GeneInterval { chr_id: 1, start: 100, end: 200 },
            GeneInterval { chr_id: 0, start: 500, end: 600 },
            GeneInterval { chr_id: 1, start: 50, end: 80 },
            GeneInterval { chr_id: 1, start: 100, end: 150 },
        ];
        genes.sort_unstable();
        assert_eq!(
            genes,
            vec![
                GeneInterval { chr_id: 0, start: 500, end: 600 },
                GeneInterval { chr_id: 1, start: 50, end: 80 },
                GeneInterval { chr_id: 1, start: 100, end: 150 },
                GeneInterval { chr_id: 1, start: 100, end: 200 },
            ]
        );
    }

    #[test]
    fn parse_coord_fa_extracts_header_fields_and_ignores_sequence_wrapping() {
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("coord.fa");
        std::fs::write(
            &path,
            ">GENE1 chr19 100 200 +\nACGTACGT\n>GENE2 chr1 5000 6000 -\nTTTTGGGG\n",
        )
        .unwrap();

        let records = parse_coord_fa(&path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, "GENE1");
        assert_eq!(records[0].chrom, "chr19");
        assert_eq!(records[0].start, 100);
        assert_eq!(records[0].end, 200);
        assert_eq!(records[0].strand, "+");
        assert_eq!(records[0].seq, "ACGTACGT");
        assert_eq!(records[1].name, "GENE2");
        assert_eq!(records[1].chrom, "chr1");
        assert_eq!(records[1].strand, "-");
    }

    #[test]
    fn parse_coord_fa_rejects_malformed_header() {
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("coord.fa");
        std::fs::write(&path, ">GENE1 chr19 100\nACGT\n").unwrap();
        let result = parse_coord_fa(&path);
        assert!(result.is_err());
    }

    // -- Parallel pass-1 regression: threads in {1, 2, 4} must produce
    // byte-identical output to the original single-threaded reference
    // implementation (`run_pass1`), and to each other. This is the
    // unum-core-local half of the P1 parallelism task's byte-identity
    // requirement -- crates/unum-sys's diff_bam_extract.rs additionally
    // proves the `threads == 1` path matches the REAL oracle; this test
    // proves `run_pass1_with_threads` at threads>1 never diverges from that
    // already-oracle-matching `threads == 1`/original-`run_pass1` behavior.

    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::{self, Header, Writer};

    /// A 400bp base-balanced (non-low-complexity) synthetic reference,
    /// long enough to carve multiple distinct 80-100bp substrings that all
    /// pass `HasHitInSet` at the default `hitLenRequired`.
    const PARALLEL_TEST_REF: &str = "ACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATGGTACGTACGTACGTACGGGCATTCATGGCATTCATGGCATTCATGACGTTAGCACGTTAGCACGTTAGCACGTTAGCTGACCATGTGACCATGTGACCATGTGACCATGGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATG";

    fn parallel_test_ref_fasta(tmp: &std::path::Path) -> std::path::PathBuf {
        let path = tmp.join("ref.fa");
        std::fs::write(&path, format!(">only\n{PARALLEL_TEST_REF}\n")).unwrap();
        path
    }

    fn noise(len: usize) -> Vec<u8> {
        let pattern =
            b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
        (0..len).map(|i| pattern[i % pattern.len()]).collect()
    }

    /// Builds a coordinate-sorted paired BAM exercising: on-target aligned
    /// pair (candidates path), alt-chrom aligned pair (candidates path,
    /// `ValidAlternativeChrom`), off-target aligned pair (dropped, never
    /// reaches the candidate filter), low-complexity on-target pair
    /// (dropped by `IsLowComplexity`), and an unaligned-template pair
    /// (direct-emit path). Two contigs: `chr1` (gene interval
    /// `[1000, 5000]`) and `chr1_alt` (alt contig, name contains `_`).
    #[allow(clippy::too_many_lines)]
    fn build_paired_test_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq1 = HeaderRecord::new(b"SQ");
        sq1.push_tag(b"SN", "chr1");
        sq1.push_tag(b"LN", 1_000_000);
        header.push_record(&sq1);
        let mut sq2 = HeaderRecord::new(b"SQ");
        sq2.push_tag(b"SN", "chr1_alt");
        sq2.push_tag(b"LN", 1_000_000);
        header.push_record(&sq2);

        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        // On-target: inside [1000, 5000] on chr1, ref-substring SEQ.
        let ot_seq1 = &ref_bytes[0..90];
        let ot_seq2 = &ref_bytes[100..190];
        let mut ot1 = bam::Record::new();
        ot1.set(b"on_target", Some(&CigarString(vec![Cigar::Match(90)])), ot_seq1, &[30u8; 90]);
        ot1.set_tid(0);
        ot1.set_pos(1100);
        ot1.set_mtid(0);
        ot1.set_mpos(1300);
        ot1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&ot1).unwrap();
        let mut ot2 = bam::Record::new();
        ot2.set(b"on_target", Some(&CigarString(vec![Cigar::Match(90)])), ot_seq2, &[30u8; 90]);
        ot2.set_tid(0);
        ot2.set_pos(1300);
        ot2.set_mtid(0);
        ot2.set_mpos(1100);
        ot2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&ot2).unwrap();

        // Off-target: on chr1 but far outside [1000, 5000] -- dropped before
        // any candidate-filter call.
        let off_seq = noise(80);
        let mut off1 = bam::Record::new();
        off1.set(b"off_target", Some(&CigarString(vec![Cigar::Match(80)])), &off_seq, &[30u8; 80]);
        off1.set_tid(0);
        off1.set_pos(50_000);
        off1.set_mtid(0);
        off1.set_mpos(50_200);
        off1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&off1).unwrap();
        let mut off2 = bam::Record::new();
        off2.set(b"off_target", Some(&CigarString(vec![Cigar::Match(80)])), &off_seq, &[30u8; 80]);
        off2.set_tid(0);
        off2.set_pos(50_200);
        off2.set_mtid(0);
        off2.set_mpos(50_000);
        off2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&off2).unwrap();

        // Low-complexity on-target pair: inside the gene interval but a
        // homopolymer SEQ -- must be dropped by IsLowComplexity.
        let homopolymer = vec![b'A'; 90];
        let mut lc1 = bam::Record::new();
        lc1.set(
            b"low_complexity",
            Some(&CigarString(vec![Cigar::Match(90)])),
            &homopolymer,
            &[30u8; 90],
        );
        lc1.set_tid(0);
        lc1.set_pos(2000);
        lc1.set_mtid(0);
        lc1.set_mpos(2200);
        lc1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&lc1).unwrap();
        let mut lc2 = bam::Record::new();
        lc2.set(
            b"low_complexity",
            Some(&CigarString(vec![Cigar::Match(90)])),
            &homopolymer,
            &[30u8; 90],
        );
        lc2.set_tid(0);
        lc2.set_pos(2200);
        lc2.set_mtid(0);
        lc2.set_mpos(2000);
        lc2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&lc2).unwrap();

        // Alt-chrom pair: mapped to "chr1_alt" (name contains '_'), ref-
        // substring SEQ -- exercises ValidAlternativeChrom's candidates path.
        let alt_seq1 = &ref_bytes[200..280];
        let alt_seq2 = &ref_bytes[300..380];
        let mut alt1 = bam::Record::new();
        alt1.set(b"alt_chrom", Some(&CigarString(vec![Cigar::Match(80)])), alt_seq1, &[30u8; 80]);
        alt1.set_tid(1);
        alt1.set_pos(500);
        alt1.set_mtid(1);
        alt1.set_mpos(700);
        alt1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&alt1).unwrap();
        let mut alt2 = bam::Record::new();
        alt2.set(b"alt_chrom", Some(&CigarString(vec![Cigar::Match(80)])), alt_seq2, &[30u8; 80]);
        alt2.set_tid(1);
        alt2.set_pos(700);
        alt2.set_mtid(1);
        alt2.set_mpos(500);
        alt2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&alt2).unwrap();

        // Unaligned-template pair: two consecutive unmapped records with
        // ref-substring SEQ -- direct-emit path.
        let um_seq1 = &ref_bytes[0..70];
        let um_seq2 = &ref_bytes[100..170];
        let mut um1 = bam::Record::new();
        um1.set(b"unaligned", None, um_seq1, &[25u8; 70]);
        um1.set_tid(-1);
        um1.set_pos(-1);
        um1.set_mtid(-1);
        um1.set_mpos(-1);
        um1.set_flags(0x1 | 0x4 | 0x8 | 0x40);
        writer.write(&um1).unwrap();
        let mut um2 = bam::Record::new();
        um2.set(b"unaligned", None, um_seq2, &[25u8; 70]);
        um2.set_tid(-1);
        um2.set_pos(-1);
        um2.set_mtid(-1);
        um2.set_mpos(-1);
        um2.set_flags(0x1 | 0x4 | 0x8 | 0x80);
        writer.write(&um2).unwrap();

        // Padding: several more on-target pairs so the parallel batching
        // path (threads > 1) actually spans multiple sites, not just one.
        for i in 0..PARALLEL_PAD_COUNT {
            let off = (i as usize * 7) % (ref_bytes.len() - 80);
            let pad_seq1 = &ref_bytes[off..off + 80];
            let pad_seq2 = &ref_bytes
                [(off + 5) % (ref_bytes.len() - 80)..(off + 5) % (ref_bytes.len() - 80) + 80];
            let name = format!("pad_{i}");
            let p1 = 1500 + i64::from(i) * 3;
            let p2 = p1 + 200;
            let mut r1 = bam::Record::new();
            r1.set(
                name.as_bytes(),
                Some(&CigarString(vec![Cigar::Match(80)])),
                pad_seq1,
                &[30u8; 80],
            );
            r1.set_tid(0);
            r1.set_pos(p1);
            r1.set_mtid(0);
            r1.set_mpos(p2);
            r1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
            writer.write(&r1).unwrap();
            let mut r2 = bam::Record::new();
            r2.set(
                name.as_bytes(),
                Some(&CigarString(vec![Cigar::Match(80)])),
                pad_seq2,
                &[30u8; 80],
            );
            r2.set_tid(0);
            r2.set_pos(p2);
            r2.set_mtid(0);
            r2.set_mpos(p1);
            r2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
            writer.write(&r2).unwrap();
        }

        drop(writer);
    }

    #[derive(Default)]
    struct VecSink {
        pairs: Vec<(ReadRecord, Option<ReadRecord>)>,
    }

    impl CandidateSink for VecSink {
        fn emit_pair(&mut self, r1: &ReadRecord, r2: Option<&ReadRecord>) -> Result<()> {
            self.pairs.push((r1.clone(), r2.cloned()));
            Ok(())
        }
    }

    /// Runs the ORIGINAL single-threaded reference `run_pass1` end to end
    /// (setup + pass 1 + pass 2), returning the emitted pairs in order.
    fn run_full_reference(bam_path: &std::path::Path, ref_fasta: &std::path::Path) -> VecSink {
        let mut filter =
            RefKmerFilter::from_reference_fasta(ref_fasta, INITIAL_KMER_LENGTH).unwrap();
        let mut alignments = Alignments::open(bam_path).unwrap();
        let genes = build_genes(
            &alignments,
            &[CoordRecord {
                name: "only".to_string(),
                chrom: "chr1".to_string(),
                start: 1000,
                end: 5000,
                strand: "+".to_string(),
                seq: PARALLEL_TEST_REF.to_string(),
            }],
        )
        .unwrap();

        let general_info = alignments.general_info(true).unwrap();
        alignments.rewind().unwrap();
        let mut hit_len_required =
            compute_hit_len_required(general_info.frag_stdev, general_info.read_len);
        filter.set_hit_len_required(hit_len_required);
        let inferred = filter.infer_kmer_length();
        if inferred > filter.kmer_length() {
            filter.update_kmer_length(inferred);
            if inferred > usize::try_from(hit_len_required).unwrap_or(0) {
                hit_len_required = i32::try_from(inferred).unwrap_or(i32::MAX);
                filter.set_hit_len_required(hit_len_required);
            }
        }
        let single_end = general_info.frag_stdev == 0;

        let mut sink = VecSink::default();
        let (candidates, _pass1_emitted) =
            run_pass1(&mut alignments, &mut filter, &genes, single_end, false, -1, &mut sink)
                .unwrap();
        alignments.rewind().unwrap();
        if !single_end {
            run_pass2(&mut alignments, candidates, false, -1, &mut sink).unwrap();
        }
        sink
    }

    /// Runs [`extract_from_bam_with_threads`] end to end at the given
    /// `threads`, returning the emitted pairs in order.
    fn run_full_parallel(
        bam_path: &std::path::Path,
        ref_fasta: &std::path::Path,
        threads: usize,
    ) -> VecSink {
        let mut filter =
            RefKmerFilter::from_reference_fasta(ref_fasta, INITIAL_KMER_LENGTH).unwrap();
        let mut alignments = Alignments::open(bam_path).unwrap();
        let genes = build_genes(
            &alignments,
            &[CoordRecord {
                name: "only".to_string(),
                chrom: "chr1".to_string(),
                start: 1000,
                end: 5000,
                strand: "+".to_string(),
                seq: PARALLEL_TEST_REF.to_string(),
            }],
        )
        .unwrap();

        let mut sink = VecSink::default();
        extract_from_bam_with_threads(
            &mut alignments,
            &mut filter,
            &genes,
            false,
            -1,
            threads,
            &mut sink,
        )
        .unwrap();
        sink
    }

    /// Full per-record signature: every field the sink emits for each pair
    /// (`r1.id`/`seq`/`qual` and, when present, `r2.id`/`seq`/`qual`), so the
    /// parallel-vs-reference regression tests enforce byte-identity of the
    /// complete FASTQ output rather than a subset.
    #[allow(clippy::type_complexity)]
    fn sink_signature(
        sink: &VecSink,
    ) -> Vec<(String, Vec<u8>, Option<Vec<u8>>, Option<String>, Option<Vec<u8>>, Option<Vec<u8>>)>
    {
        sink.pairs
            .iter()
            .map(|(r1, r2)| {
                let r2 = r2.as_ref();
                (
                    r1.id.clone(),
                    r1.seq.clone(),
                    r1.qual.clone(),
                    r2.map(|r| r.id.clone()),
                    r2.map(|r| r.seq.clone()),
                    r2.and_then(|r| r.qual.clone()),
                )
            })
            .collect()
    }

    #[test]
    fn paired_parallel_pass1_matches_reference_at_threads_1_2_4() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("test.bam");
        build_paired_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let reference = run_full_reference(&bam_path, &ref_fasta);
        let reference_sig = sink_signature(&reference);
        assert!(!reference_sig.is_empty(), "fixture must produce at least one candidate");
        // Sanity: confirm expected categories are present/absent.
        let ids: Vec<&str> = reference.pairs.iter().map(|(r1, _)| r1.id.as_str()).collect();
        assert!(ids.contains(&"on_target"), "on-target pair should be emitted: {ids:?}");
        assert!(ids.contains(&"alt_chrom"), "alt-chrom pair should be emitted: {ids:?}");
        assert!(ids.contains(&"unaligned"), "unaligned-template pair should be emitted: {ids:?}");
        assert!(!ids.contains(&"off_target"), "off-target pair must be dropped: {ids:?}");
        assert!(!ids.contains(&"low_complexity"), "low-complexity pair must be dropped: {ids:?}");

        for threads in [1usize, 2, 4] {
            let parallel = run_full_parallel(&bam_path, &ref_fasta, threads);
            let parallel_sig = sink_signature(&parallel);
            assert_eq!(
                parallel_sig, reference_sig,
                "threads={threads} output diverged from the single-threaded reference"
            );
        }

        // Cross-check threads=2 vs threads=4 too, not just each-vs-reference.
        let p2 = sink_signature(&run_full_parallel(&bam_path, &ref_fasta, 2));
        let p4 = sink_signature(&run_full_parallel(&bam_path, &ref_fasta, 4));
        assert_eq!(p2, p4, "threads=2 and threads=4 outputs diverged from each other");
    }

    /// Builds a single-end BAM (frag_stdev == 0 path): on-target aligned
    /// reads (one with a duplicate QNAME, exercising `used_name` dedup) and
    /// an off-target aligned read (dropped).
    fn build_single_end_test_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);

        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        let seq = &ref_bytes[0..90];
        let mut r1 = bam::Record::new();
        r1.set(b"single_on_target", Some(&CigarString(vec![Cigar::Match(90)])), seq, &[30u8; 90]);
        r1.set_tid(0);
        r1.set_pos(1100);
        r1.set_mtid(-1);
        r1.set_mpos(-1);
        r1.set_flags(0);
        writer.write(&r1).unwrap();

        // Duplicate QNAME, second alignment -- usedName dedup.
        let mut r1_dup = bam::Record::new();
        r1_dup.set(
            b"single_on_target",
            Some(&CigarString(vec![Cigar::Match(90)])),
            seq,
            &[30u8; 90],
        );
        r1_dup.set_tid(0);
        r1_dup.set_pos(1200);
        r1_dup.set_mtid(-1);
        r1_dup.set_mpos(-1);
        r1_dup.set_flags(0x100);
        writer.write(&r1_dup).unwrap();

        let off_seq = noise(80);
        let mut r2 = bam::Record::new();
        r2.set(
            b"single_off_target",
            Some(&CigarString(vec![Cigar::Match(80)])),
            &off_seq,
            &[30u8; 80],
        );
        r2.set_tid(0);
        r2.set_pos(50_000);
        r2.set_mtid(-1);
        r2.set_mpos(-1);
        r2.set_flags(0);
        writer.write(&r2).unwrap();

        // A batch of additional on-target single-end reads so threads > 1
        // has multiple sites to parallelize over.
        for i in 0..PARALLEL_PAD_COUNT {
            let off = (i as usize * 7) % (ref_bytes.len() - 80);
            let seq = &ref_bytes[off..off + 80];
            let name = format!("single_pad_{i}");
            let mut r = bam::Record::new();
            r.set(name.as_bytes(), Some(&CigarString(vec![Cigar::Match(80)])), seq, &[30u8; 80]);
            r.set_tid(0);
            r.set_pos(1500 + i64::from(i) * 3);
            r.set_mtid(-1);
            r.set_mpos(-1);
            r.set_flags(0);
            writer.write(&r).unwrap();
        }

        drop(writer);
    }

    #[test]
    fn single_end_parallel_pass1_matches_reference_at_threads_1_2_4() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("single.bam");
        build_single_end_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let reference = run_full_reference(&bam_path, &ref_fasta);
        let reference_sig = sink_signature(&reference);
        assert!(!reference_sig.is_empty(), "fixture must produce at least one candidate");

        let ids: Vec<&str> = reference.pairs.iter().map(|(r1, _)| r1.id.as_str()).collect();
        assert_eq!(
            ids.iter().filter(|&&id| id == "single_on_target").count(),
            1,
            "usedName dedup must emit the duplicate-QNAME read only once: {ids:?}"
        );
        assert!(!ids.contains(&"single_off_target"), "off-target read must be dropped: {ids:?}");

        for threads in [1usize, 2, 4] {
            let parallel = run_full_parallel(&bam_path, &ref_fasta, threads);
            let parallel_sig = sink_signature(&parallel);
            assert_eq!(
                parallel_sig, reference_sig,
                "threads={threads} output diverged from the single-threaded reference"
            );
        }
    }

    /// Regression for the unaligned-template mate-mismatch error: both the
    /// original reference `run_pass1` AND `run_pass1_with_threads` must
    /// reject a lonely unmapped record (no adjacent mate) identically,
    /// regardless of `threads`.
    #[test]
    fn missing_unaligned_mate_errors_at_any_thread_count() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("missing_mate.bam");
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);
        {
            let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
            let seq = &PARALLEL_TEST_REF.as_bytes()[0..60];
            let mut r = bam::Record::new();
            r.set(b"lonely", None, seq, &[25u8; 60]);
            r.set_tid(-1);
            r.set_pos(-1);
            r.set_mtid(-1);
            r.set_mpos(-1);
            r.set_flags(0x1 | 0x4 | 0x8 | 0x40);
            writer.write(&r).unwrap();
        }

        for threads in [1usize, 2, 4] {
            let mut filter =
                RefKmerFilter::from_reference_fasta(&ref_fasta, INITIAL_KMER_LENGTH).unwrap();
            let mut alignments = Alignments::open(&bam_path).unwrap();
            let genes = build_genes(
                &alignments,
                &[CoordRecord {
                    name: "only".to_string(),
                    chrom: "chr1".to_string(),
                    start: 1000,
                    end: 5000,
                    strand: "+".to_string(),
                    seq: PARALLEL_TEST_REF.to_string(),
                }],
            )
            .unwrap();
            let mut sink = VecSink::default();
            let result = extract_from_bam_with_threads(
                &mut alignments,
                &mut filter,
                &genes,
                false,
                -1,
                threads,
                &mut sink,
            );
            assert!(result.is_err(), "threads={threads} must reject a missing unaligned mate");
            assert!(
                result.unwrap_err().to_string().contains("not showing up together"),
                "threads={threads}: unexpected error message"
            );
        }
    }

    // -- On-target/non-candidate regression (CRITICAL fix): `BamExtractor.cpp:804-851`
    // emits/records on-target reads UNCONDITIONALLY once they pass the `IsLowComplexity`
    // check -- there is no `IsGoodCandidate` gate on the on-target path. A prior version
    // of `evaluate_pass1_site`'s `SingleEndOnTarget | PairedOnTargetCandidate` arm
    // incorrectly called `is_good_candidate_with_scratch`, which dropped on-target reads
    // that overlap a gene by coordinate but are NOT kmer candidates (e.g. soft-clipped,
    // divergent-allele, or noisy reads) -- a byte-identity divergence from the oracle,
    // reproducible even at `threads == 1`. These tests build an on-target read whose SEQ
    // is `noise(..)` (not a substring of `PARALLEL_TEST_REF`, so it fails
    // `is_good_candidate_with_scratch`) and assert it is emitted anyway, matching both
    // the untouched reference `run_pass1` and `run_pass1_with_threads` at
    // threads in {1, 4, 8}.

    /// Paired variant of [`build_paired_test_bam`], adding one extra on-target pair
    /// (`ontarget_noise`) whose SEQ is low-complexity-clean but kmer-non-candidate noise.
    #[allow(clippy::too_many_lines)]
    fn build_paired_test_bam_with_ontarget_noise(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq1 = HeaderRecord::new(b"SQ");
        sq1.push_tag(b"SN", "chr1");
        sq1.push_tag(b"LN", 1_000_000);
        header.push_record(&sq1);

        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        // On-target, ref-substring SEQ (a good candidate) -- baseline sanity pair.
        let ot_seq1 = &ref_bytes[0..90];
        let ot_seq2 = &ref_bytes[100..190];
        let mut ot1 = bam::Record::new();
        ot1.set(b"on_target", Some(&CigarString(vec![Cigar::Match(90)])), ot_seq1, &[30u8; 90]);
        ot1.set_tid(0);
        ot1.set_pos(1100);
        ot1.set_mtid(0);
        ot1.set_mpos(1300);
        ot1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&ot1).unwrap();
        let mut ot2 = bam::Record::new();
        ot2.set(b"on_target", Some(&CigarString(vec![Cigar::Match(90)])), ot_seq2, &[30u8; 90]);
        ot2.set_tid(0);
        ot2.set_pos(1300);
        ot2.set_mtid(0);
        ot2.set_mpos(1100);
        ot2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&ot2).unwrap();

        // On-target, NOT a kmer candidate: `noise()` is neither low-complexity (it is a
        // 4-base repeating pattern, not a homopolymer) nor a substring/near-match of
        // `PARALLEL_TEST_REF`, so `is_good_candidate_with_scratch` returns false for it
        // (the same sequence generator is used for the existing `off_target` fixture,
        // which relies on this same non-candidacy). Positioned inside the gene interval
        // `[1000, 5000]` on `chr1` so it reaches the on-target branch.
        let noise_seq1 = noise(90);
        let noise_seq2 = noise(90);
        let mut n1 = bam::Record::new();
        n1.set(
            b"ontarget_noise",
            Some(&CigarString(vec![Cigar::Match(90)])),
            &noise_seq1,
            &[30u8; 90],
        );
        n1.set_tid(0);
        n1.set_pos(1600);
        n1.set_mtid(0);
        n1.set_mpos(1800);
        n1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&n1).unwrap();
        let mut n2 = bam::Record::new();
        n2.set(
            b"ontarget_noise",
            Some(&CigarString(vec![Cigar::Match(90)])),
            &noise_seq2,
            &[30u8; 90],
        );
        n2.set_tid(0);
        n2.set_pos(1800);
        n2.set_mtid(0);
        n2.set_mpos(1600);
        n2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&n2).unwrap();

        // Off-target, ref-substring SEQ (a good candidate but out of range) -- confirms
        // coordinate filtering still drops out-of-range reads regardless of candidacy.
        let off_seq1 = &ref_bytes[0..80];
        let off_seq2 = &ref_bytes[100..180];
        let mut off1 = bam::Record::new();
        off1.set(b"off_target", Some(&CigarString(vec![Cigar::Match(80)])), off_seq1, &[30u8; 80]);
        off1.set_tid(0);
        off1.set_pos(50_000);
        off1.set_mtid(0);
        off1.set_mpos(50_200);
        off1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&off1).unwrap();
        let mut off2 = bam::Record::new();
        off2.set(b"off_target", Some(&CigarString(vec![Cigar::Match(80)])), off_seq2, &[30u8; 80]);
        off2.set_tid(0);
        off2.set_pos(50_200);
        off2.set_mtid(0);
        off2.set_mpos(50_000);
        off2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&off2).unwrap();

        // Padding: several more on-target good-candidate pairs so the parallel batching
        // path (threads > 1) spans multiple sites, not just the two under test.
        for i in 0..PARALLEL_PAD_COUNT {
            let off = (i as usize * 7) % (ref_bytes.len() - 80);
            let pad_seq1 = &ref_bytes[off..off + 80];
            let alt_off = (off + 5) % (ref_bytes.len() - 80);
            let pad_seq2 = &ref_bytes[alt_off..alt_off + 80];
            let name = format!("pad_{i}");
            let p1 = 2000 + i64::from(i) * 3;
            let p2 = p1 + 200;
            let mut r1 = bam::Record::new();
            r1.set(
                name.as_bytes(),
                Some(&CigarString(vec![Cigar::Match(80)])),
                pad_seq1,
                &[30u8; 80],
            );
            r1.set_tid(0);
            r1.set_pos(p1);
            r1.set_mtid(0);
            r1.set_mpos(p2);
            r1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
            writer.write(&r1).unwrap();
            let mut r2 = bam::Record::new();
            r2.set(
                name.as_bytes(),
                Some(&CigarString(vec![Cigar::Match(80)])),
                pad_seq2,
                &[30u8; 80],
            );
            r2.set_tid(0);
            r2.set_pos(p2);
            r2.set_mtid(0);
            r2.set_mpos(p1);
            r2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
            writer.write(&r2).unwrap();
        }

        drop(writer);
    }

    #[test]
    fn paired_ontarget_non_candidate_read_is_emitted_unconditionally() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("ontarget_noise.bam");
        build_paired_test_bam_with_ontarget_noise(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        // The untouched reference `run_pass1` is the oracle-matching ground truth (see
        // its on-target arm: no `is_good_candidate_with_scratch` call at all) -- it MUST
        // emit `ontarget_noise` despite the sequence failing the kmer candidate filter.
        let reference = run_full_reference(&bam_path, &ref_fasta);
        let reference_ids: Vec<&str> =
            reference.pairs.iter().map(|(r1, _)| r1.id.as_str()).collect();
        assert!(
            reference_ids.contains(&"ontarget_noise"),
            "reference run_pass1 must emit an on-target non-candidate read unconditionally: {reference_ids:?}"
        );
        assert!(reference_ids.contains(&"on_target"));
        assert!(!reference_ids.contains(&"off_target"), "off-target must still be dropped");

        // `run_pass1_with_threads` (the production scan/evaluate/apply split) must match
        // the reference exactly at every thread count -- this is the regression this test
        // guards: before the fix, `evaluate_pass1_site`'s on-target arm gated on
        // `is_good_candidate_with_scratch`, dropping `ontarget_noise` at ALL thread
        // counts (the bug is not thread-count-dependent).
        let reference_sig = sink_signature(&reference);
        for threads in [1usize, 4, 8] {
            let parallel = run_full_parallel(&bam_path, &ref_fasta, threads);
            let parallel_ids: Vec<&str> =
                parallel.pairs.iter().map(|(r1, _)| r1.id.as_str()).collect();
            assert!(
                parallel_ids.contains(&"ontarget_noise"),
                "threads={threads}: on-target non-candidate read must be emitted, matching the oracle: {parallel_ids:?}"
            );
            let parallel_sig = sink_signature(&parallel);
            assert_eq!(
                parallel_sig, reference_sig,
                "threads={threads}: output diverged from the oracle-matching reference"
            );
        }
    }

    /// Single-end variant: `build_single_end_test_bam` plus one on-target read whose SEQ
    /// is kmer-non-candidate noise. Same regression as the paired test above, exercised
    /// through `Pass1Site::SingleEndOnTarget` instead of `PairedOnTargetCandidate`.
    fn build_single_end_test_bam_with_ontarget_noise(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);

        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        let seq = &ref_bytes[0..90];
        let mut r1 = bam::Record::new();
        r1.set(b"single_on_target", Some(&CigarString(vec![Cigar::Match(90)])), seq, &[30u8; 90]);
        r1.set_tid(0);
        r1.set_pos(1100);
        r1.set_mtid(-1);
        r1.set_mpos(-1);
        r1.set_flags(0);
        writer.write(&r1).unwrap();

        // On-target, NOT a kmer candidate.
        let noise_seq = noise(90);
        let mut n = bam::Record::new();
        n.set(
            b"single_ontarget_noise",
            Some(&CigarString(vec![Cigar::Match(90)])),
            &noise_seq,
            &[30u8; 90],
        );
        n.set_tid(0);
        n.set_pos(1600);
        n.set_mtid(-1);
        n.set_mpos(-1);
        n.set_flags(0);
        writer.write(&n).unwrap();

        let off_seq = &ref_bytes[0..80];
        let mut r2 = bam::Record::new();
        r2.set(
            b"single_off_target",
            Some(&CigarString(vec![Cigar::Match(80)])),
            off_seq,
            &[30u8; 80],
        );
        r2.set_tid(0);
        r2.set_pos(50_000);
        r2.set_mtid(-1);
        r2.set_mpos(-1);
        r2.set_flags(0);
        writer.write(&r2).unwrap();

        // Padding for the parallel batching path.
        for i in 0..PARALLEL_PAD_COUNT {
            let off = (i as usize * 7) % (ref_bytes.len() - 80);
            let seq = &ref_bytes[off..off + 80];
            let name = format!("single_pad_{i}");
            let mut r = bam::Record::new();
            r.set(name.as_bytes(), Some(&CigarString(vec![Cigar::Match(80)])), seq, &[30u8; 80]);
            r.set_tid(0);
            r.set_pos(2000 + i64::from(i) * 3);
            r.set_mtid(-1);
            r.set_mpos(-1);
            r.set_flags(0);
            writer.write(&r).unwrap();
        }

        drop(writer);
    }

    #[test]
    fn single_end_ontarget_non_candidate_read_is_emitted_unconditionally() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("single_ontarget_noise.bam");
        build_single_end_test_bam_with_ontarget_noise(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let reference = run_full_reference(&bam_path, &ref_fasta);
        let reference_ids: Vec<&str> =
            reference.pairs.iter().map(|(r1, _)| r1.id.as_str()).collect();
        assert!(
            reference_ids.contains(&"single_ontarget_noise"),
            "reference run_pass1 must emit an on-target non-candidate read unconditionally: {reference_ids:?}"
        );
        assert!(!reference_ids.contains(&"single_off_target"), "off-target must still be dropped");

        let reference_sig = sink_signature(&reference);
        for threads in [1usize, 4, 8] {
            let parallel = run_full_parallel(&bam_path, &ref_fasta, threads);
            let parallel_ids: Vec<&str> =
                parallel.pairs.iter().map(|(r1, _)| r1.id.as_str()).collect();
            assert!(
                parallel_ids.contains(&"single_ontarget_noise"),
                "threads={threads}: on-target non-candidate read must be emitted, matching the oracle: {parallel_ids:?}"
            );
            let parallel_sig = sink_signature(&parallel);
            assert_eq!(
                parallel_sig, reference_sig,
                "threads={threads}: output diverged from the oracle-matching reference"
            );
        }
    }

    // -- Bounded-chunk pass-1 regression (memory fix): `run_pass1_chunked` at a
    // SMALL `chunk_size` must match the untouched single-scan `run_pass1` oracle
    // exactly, proving the `tag`/`candidates`/`used_name`/`pass1_emitted`
    // cross-chunk threading in `run_pass1_chunked` (see its doc comment) is
    // correct -- not just at the production `PASS1_CHUNK_SIZE` (~32k, too slow
    // to exercise boundary-crossing in a unit test), but at a `chunk_size` that
    // forces several chunk boundaries over a small, cheap fixture.

    /// A small coordinate-sorted paired BAM (~20 records: 8 on-target pairs
    /// exercising the `candidates` accumulation path, one off-target pair
    /// dropped before any candidate-filter call, one unaligned-template pair
    /// exercising the direct-emit path) -- enough sites to cross several
    /// `chunk_size = 4` chunk boundaries without needing a
    /// [`PASS1_CHUNK_SIZE`]-scale fixture.
    fn build_small_paired_test_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);

        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        // 8 on-target pairs (16 records), each inside the gene interval
        // [1000, 5000] -- candidates path (`PairedOnTargetCandidate`).
        for i in 0..8u32 {
            let off1 = (i as usize * 11) % (ref_bytes.len() - 80);
            let off2 = (off1 + 13) % (ref_bytes.len() - 80);
            let seq1 = &ref_bytes[off1..off1 + 80];
            let seq2 = &ref_bytes[off2..off2 + 80];
            let name = format!("on_target_{i}");
            let p1 = 1200 + i64::from(i) * 100;
            let p2 = p1 + 50;
            let mut r1 = bam::Record::new();
            r1.set(name.as_bytes(), Some(&CigarString(vec![Cigar::Match(80)])), seq1, &[30u8; 80]);
            r1.set_tid(0);
            r1.set_pos(p1);
            r1.set_mtid(0);
            r1.set_mpos(p2);
            r1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
            writer.write(&r1).unwrap();
            let mut r2 = bam::Record::new();
            r2.set(name.as_bytes(), Some(&CigarString(vec![Cigar::Match(80)])), seq2, &[30u8; 80]);
            r2.set_tid(0);
            r2.set_pos(p2);
            r2.set_mtid(0);
            r2.set_mpos(p1);
            r2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
            writer.write(&r2).unwrap();
        }

        // Off-target pair (2 records): on chr1 but far outside [1000, 5000] --
        // dropped before any candidate-filter call.
        let off_seq = noise(80);
        let mut off1 = bam::Record::new();
        off1.set(b"off_target", Some(&CigarString(vec![Cigar::Match(80)])), &off_seq, &[30u8; 80]);
        off1.set_tid(0);
        off1.set_pos(50_000);
        off1.set_mtid(0);
        off1.set_mpos(50_200);
        off1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&off1).unwrap();
        let mut off2 = bam::Record::new();
        off2.set(b"off_target", Some(&CigarString(vec![Cigar::Match(80)])), &off_seq, &[30u8; 80]);
        off2.set_tid(0);
        off2.set_pos(50_200);
        off2.set_mtid(0);
        off2.set_mpos(50_000);
        off2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&off2).unwrap();

        // Unaligned-template pair (2 records): two consecutive unmapped
        // records with ref-substring SEQ -- direct-emit path
        // (`Pass1Site::UnalignedPair`), read as ONE site across TWO records.
        let um_seq1 = &ref_bytes[0..70];
        let um_seq2 = &ref_bytes[100..170];
        let mut um1 = bam::Record::new();
        um1.set(b"unaligned", None, um_seq1, &[25u8; 70]);
        um1.set_tid(-1);
        um1.set_pos(-1);
        um1.set_mtid(-1);
        um1.set_mpos(-1);
        um1.set_flags(0x1 | 0x4 | 0x8 | 0x40);
        writer.write(&um1).unwrap();
        let mut um2 = bam::Record::new();
        um2.set(b"unaligned", None, um_seq2, &[25u8; 70]);
        um2.set_tid(-1);
        um2.set_pos(-1);
        um2.set_mtid(-1);
        um2.set_mpos(-1);
        um2.set_flags(0x1 | 0x4 | 0x8 | 0x80);
        writer.write(&um2).unwrap();

        drop(writer);
    }

    /// Runs the ORIGINAL single-threaded, single-scan `run_pass1` oracle
    /// end-to-end (setup + pass 1 only, no pass 2), returning its
    /// `candidates` key set, `pass1_emitted` count, and the pairs it emitted
    /// directly during pass 1.
    fn run_pass1_oracle(
        bam_path: &std::path::Path,
        ref_fasta: &std::path::Path,
    ) -> (std::collections::HashSet<String>, u64, VecSink) {
        let mut filter =
            RefKmerFilter::from_reference_fasta(ref_fasta, INITIAL_KMER_LENGTH).unwrap();
        let mut alignments = Alignments::open(bam_path).unwrap();
        let genes = build_genes(
            &alignments,
            &[CoordRecord {
                name: "only".to_string(),
                chrom: "chr1".to_string(),
                start: 1000,
                end: 5000,
                strand: "+".to_string(),
                seq: PARALLEL_TEST_REF.to_string(),
            }],
        )
        .unwrap();

        let general_info = alignments.general_info(true).unwrap();
        alignments.rewind().unwrap();
        let mut hit_len_required =
            compute_hit_len_required(general_info.frag_stdev, general_info.read_len);
        filter.set_hit_len_required(hit_len_required);
        let inferred = filter.infer_kmer_length();
        if inferred > filter.kmer_length() {
            filter.update_kmer_length(inferred);
            if inferred > usize::try_from(hit_len_required).unwrap_or(0) {
                hit_len_required = i32::try_from(inferred).unwrap_or(i32::MAX);
                filter.set_hit_len_required(hit_len_required);
            }
        }
        let single_end = general_info.frag_stdev == 0;

        let mut sink = VecSink::default();
        let (candidates, pass1_emitted) =
            run_pass1(&mut alignments, &mut filter, &genes, single_end, false, -1, &mut sink)
                .unwrap();
        (candidates.into_keys().collect(), pass1_emitted, sink)
    }

    /// Runs [`run_pass1_chunked`] end-to-end (setup + pass 1 only) at the
    /// given `threads`/`chunk_size`, returning the same shape as
    /// [`run_pass1_oracle`] for direct comparison.
    fn run_pass1_chunked_test(
        bam_path: &std::path::Path,
        ref_fasta: &std::path::Path,
        threads: usize,
        chunk_size: usize,
        selection: Selection,
    ) -> (std::collections::HashSet<String>, u64, VecSink) {
        let mut filter =
            RefKmerFilter::from_reference_fasta(ref_fasta, INITIAL_KMER_LENGTH).unwrap();
        let mut alignments = Alignments::open(bam_path).unwrap();
        let genes = build_genes(
            &alignments,
            &[CoordRecord {
                name: "only".to_string(),
                chrom: "chr1".to_string(),
                start: 1000,
                end: 5000,
                strand: "+".to_string(),
                seq: PARALLEL_TEST_REF.to_string(),
            }],
        )
        .unwrap();

        let general_info = alignments.general_info(true).unwrap();
        alignments.rewind().unwrap();
        let mut hit_len_required =
            compute_hit_len_required(general_info.frag_stdev, general_info.read_len);
        filter.set_hit_len_required(hit_len_required);
        let inferred = filter.infer_kmer_length();
        if inferred > filter.kmer_length() {
            filter.update_kmer_length(inferred);
            if inferred > usize::try_from(hit_len_required).unwrap_or(0) {
                hit_len_required = i32::try_from(inferred).unwrap_or(i32::MAX);
                filter.set_hit_len_required(hit_len_required);
            }
        }
        let single_end = general_info.frag_stdev == 0;

        let mut sink = VecSink::default();
        let (candidates, pass1_emitted) = run_pass1_chunked(
            &mut alignments,
            &filter,
            &genes,
            single_end,
            false,
            -1,
            threads,
            chunk_size,
            &mut sink,
            selection,
        )
        .unwrap();
        (candidates.into_keys().collect(), pass1_emitted, sink)
    }

    #[test]
    fn bounded_chunk_pass1_matches_oracle_across_chunk_boundaries() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("chunked.bam");
        build_small_paired_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let (oracle_candidates, oracle_emitted, oracle_sink) =
            run_pass1_oracle(&bam_path, &ref_fasta);
        assert_eq!(oracle_candidates.len(), 8, "expected all 8 on-target pairs as candidates");
        assert_eq!(oracle_emitted, 1, "expected the single unaligned pair to be emitted directly");
        let oracle_sig = sink_signature(&oracle_sink);

        // `chunk_size = 4` over a 20-record fixture crosses several chunk
        // boundaries (5+ chunks), at both a sequential (`threads = 1`) and a
        // parallel (`threads = 4`) evaluate step.
        for threads in [1usize, 4] {
            let (candidates, pass1_emitted, sink) =
                run_pass1_chunked_test(&bam_path, &ref_fasta, threads, 4, Selection::Alignment);
            assert_eq!(
                candidates, oracle_candidates,
                "threads={threads}, chunk_size=4: candidates diverged from the oracle across chunk boundaries"
            );
            assert_eq!(
                pass1_emitted, oracle_emitted,
                "threads={threads}, chunk_size=4: pass1_emitted diverged from the oracle"
            );
            assert_eq!(
                sink_signature(&sink),
                oracle_sig,
                "threads={threads}, chunk_size=4: emitted pairs diverged from the oracle across chunk boundaries"
            );
        }
    }

    /// A small coordinate-sorted SINGLE-END BAM (~20 on-target aligned reads,
    /// single-end header/flags: `mtid = -1`, no paired/proper-pair flags) with
    /// a DUPLICATE QNAME (`shared_name`) appearing twice, at read indices 3 and
    /// 12 -- i.e. in different `chunk_size = 4` chunks (chunk 0 vs chunk 3), so
    /// the `used_name` dedup decision for the second occurrence depends on
    /// `used_name` state accumulated in an EARLIER chunk. This is the fixture
    /// that gives `used_name` cross-chunk threading teeth: a per-chunk
    /// `used_name` reset would let the second `shared_name` read slip through
    /// the dedup and be emitted a second time.
    fn build_small_single_end_test_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);

        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        // 20 on-target single-end reads inside the gene interval [1000, 5000].
        // Read indices 3 and 12 share the QNAME `shared_name` (the second one a
        // secondary alignment, flag 0x100) -- the `used_name` dedup path.
        for i in 0..20u32 {
            let off = (i as usize * 11) % (ref_bytes.len() - 80);
            let seq = &ref_bytes[off..off + 80];
            let name =
                if i == 3 || i == 12 { "shared_name".to_string() } else { format!("se_{i}") };
            // Second occurrence of `shared_name` marked secondary (0x100),
            // mirroring a realistic duplicate alignment; the dedup itself keys
            // on the QNAME, not the flag.
            let flags = if i == 12 { 0x100 } else { 0 };
            let mut r = bam::Record::new();
            r.set(name.as_bytes(), Some(&CigarString(vec![Cigar::Match(80)])), seq, &[30u8; 80]);
            r.set_tid(0);
            r.set_pos(1100 + i64::from(i) * 100);
            r.set_mtid(-1);
            r.set_mpos(-1);
            r.set_flags(flags);
            writer.write(&r).unwrap();
        }

        drop(writer);
    }

    #[test]
    fn bounded_chunk_pass1_single_end_matches_oracle_across_chunk_boundaries() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("chunked_single_end.bam");
        build_small_single_end_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let (oracle_candidates, oracle_emitted, oracle_sink) =
            run_pass1_oracle(&bam_path, &ref_fasta);
        // Single-end never populates `candidates`; all output is emitted
        // directly in pass 1. 20 reads, one QNAME duplicated -> 19 distinct
        // names emitted once each (the second `shared_name` is deduped by
        // `used_name`). If a per-chunk `used_name` reset regression existed,
        // the second `shared_name` (chunk 3) would not see the chunk-0 entry
        // and would be emitted a second time -> `oracle_emitted` would be 20
        // and this fixture-level sanity check pins the deduped count.
        assert!(oracle_candidates.is_empty(), "single-end input must not populate candidates");
        assert_eq!(
            oracle_emitted, 19,
            "expected 19 emitted reads (one QNAME deduped via used_name)"
        );
        let oracle_sig = sink_signature(&oracle_sink);

        // `chunk_size = 4` over 20 reads crosses several boundaries (5 chunks);
        // the duplicate `shared_name` at read indices 3 and 12 straddles ~3
        // chunk boundaries, so matching the oracle here proves `used_name`
        // accumulates across chunks. A per-chunk `used_name` reset would leave
        // the chunk-3 duplicate undeduped -> higher `pass1_emitted` and an
        // extra emitted read -> both asserts below would fail.
        for threads in [1usize, 4] {
            let (candidates, pass1_emitted, sink) =
                run_pass1_chunked_test(&bam_path, &ref_fasta, threads, 4, Selection::Alignment);
            assert_eq!(
                candidates, oracle_candidates,
                "threads={threads}, chunk_size=4: single-end candidates diverged from the oracle"
            );
            assert_eq!(
                pass1_emitted, oracle_emitted,
                "threads={threads}, chunk_size=4: pass1_emitted diverged from the oracle (used_name cross-chunk dedup)"
            );
            assert_eq!(
                sink_signature(&sink),
                oracle_sig,
                "threads={threads}, chunk_size=4: emitted reads diverged from the oracle across chunk boundaries"
            );
        }
    }

    // -- `Selection::NoAlignment` (Class-A k-mer selection): every PRIMARY
    // read is k-mer-tested on its own sequence, position/gene-interval/`tag`
    // entirely bypassed. `Selection::Alignment` stays byte-identical to its
    // pre-Task-3 behavior -- every test ABOVE this point exercises it
    // exclusively (via `run_pass1_chunked_test(..., Selection::Alignment)`),
    // so those tests double as the byte-identity regression guard.

    /// A 3-pair coordinate BAM for the no-alignment-vs-alignment selection
    /// differential: `on_target` (aligned inside the gene interval `[1000,
    /// 5000]`, ref-substring SEQ -- selected by BOTH modes), `off_target_kmer`
    /// (aligned to chr1 at position 50000, far outside the gene interval, but
    /// with a ref-substring SEQ that WOULD pass the k-mer candidate filter --
    /// [`Selection::Alignment`] drops this purely by position, via the
    /// `tag`/gene-interval-scan `continue` in `scan_pass1_chunk`, BEFORE any
    /// candidate-filter call; [`Selection::NoAlignment`] selects it via
    /// k-mer, position-independent), and `unaligned_kmer` (both mates
    /// unmapped, ref-substring SEQ -- `Selection::Alignment`'s direct-emit
    /// `UnalignedPair` path vs. `Selection::NoAlignment`'s per-read
    /// `KmerCandidate`/`candidates`-map path).
    fn build_kmer_selection_test_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);

        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        // on_target: inside [1000, 5000], ref-substring SEQ -- selected by
        // both modes.
        let ot_seq1 = &ref_bytes[0..90];
        let ot_seq2 = &ref_bytes[100..190];
        let mut ot1 = bam::Record::new();
        ot1.set(b"on_target", Some(&CigarString(vec![Cigar::Match(90)])), ot_seq1, &[30u8; 90]);
        ot1.set_tid(0);
        ot1.set_pos(1100);
        ot1.set_mtid(0);
        ot1.set_mpos(1300);
        ot1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&ot1).unwrap();
        let mut ot2 = bam::Record::new();
        ot2.set(b"on_target", Some(&CigarString(vec![Cigar::Match(90)])), ot_seq2, &[30u8; 90]);
        ot2.set_tid(0);
        ot2.set_pos(1300);
        ot2.set_mtid(0);
        ot2.set_mpos(1100);
        ot2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&ot2).unwrap();

        // off_target_kmer: aligned to chr1 far outside [1000, 5000] (written
        // AFTER on_target so the `tag` gene-interval pointer has already
        // correctly classified on_target before advancing past the only
        // gene), but a ref-substring SEQ that DOES pass the k-mer filter.
        let far_seq1 = &ref_bytes[200..280];
        let far_seq2 = &ref_bytes[300..380];
        let mut far1 = bam::Record::new();
        far1.set(
            b"off_target_kmer",
            Some(&CigarString(vec![Cigar::Match(80)])),
            far_seq1,
            &[30u8; 80],
        );
        far1.set_tid(0);
        far1.set_pos(50_000);
        far1.set_mtid(0);
        far1.set_mpos(50_200);
        far1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
        writer.write(&far1).unwrap();
        let mut far2 = bam::Record::new();
        far2.set(
            b"off_target_kmer",
            Some(&CigarString(vec![Cigar::Match(80)])),
            far_seq2,
            &[30u8; 80],
        );
        far2.set_tid(0);
        far2.set_pos(50_200);
        far2.set_mtid(0);
        far2.set_mpos(50_000);
        far2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
        writer.write(&far2).unwrap();

        // unaligned_kmer: two consecutive unmapped records, ref-substring SEQ.
        let um_seq1 = &ref_bytes[0..70];
        let um_seq2 = &ref_bytes[100..170];
        let mut um1 = bam::Record::new();
        um1.set(b"unaligned_kmer", None, um_seq1, &[25u8; 70]);
        um1.set_tid(-1);
        um1.set_pos(-1);
        um1.set_mtid(-1);
        um1.set_mpos(-1);
        um1.set_flags(0x1 | 0x4 | 0x8 | 0x40);
        writer.write(&um1).unwrap();
        let mut um2 = bam::Record::new();
        um2.set(b"unaligned_kmer", None, um_seq2, &[25u8; 70]);
        um2.set_tid(-1);
        um2.set_pos(-1);
        um2.set_mtid(-1);
        um2.set_mpos(-1);
        um2.set_flags(0x1 | 0x4 | 0x8 | 0x80);
        writer.write(&um2).unwrap();

        drop(writer);
    }

    #[test]
    fn no_alignment_pass1_selects_by_kmer_not_position() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("kmer_selection.bam");
        build_kmer_selection_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        // `Selection::NoAlignment`: every primary read is k-mer-tested on its
        // own sequence, position/gene-interval entirely bypassed -- all
        // three pairs' k-mer-matching sequences pass, so all three QNAMEs
        // land in `candidates` (paired input only records the name; nothing
        // is emitted directly in pass 1).
        let (no_alignment_candidates, no_alignment_emitted, no_alignment_sink) =
            run_pass1_chunked_test(&bam_path, &ref_fasta, 1, 1000, Selection::NoAlignment);
        assert_eq!(
            no_alignment_candidates,
            ["on_target", "off_target_kmer", "unaligned_kmer"]
                .into_iter()
                .map(str::to_string)
                .collect(),
            "no-alignment must record every k-mer-passing primary read's QNAME, regardless of \
             alignment position"
        );
        assert_eq!(
            no_alignment_emitted, 0,
            "no-alignment (paired input) records names via `candidates`, never emits directly"
        );
        assert!(no_alignment_sink.pairs.is_empty(), "no-alignment (paired) must not emit directly");

        // `Selection::Alignment`: unchanged, position-driven classification.
        // `off_target_kmer` is dropped purely by position (the `tag`/
        // gene-interval-scan `continue`, BEFORE any candidate-filter call) --
        // it never becomes a `Pass1Site` at all, k-mer match notwithstanding.
        let (alignment_candidates, alignment_emitted, alignment_sink) =
            run_pass1_chunked_test(&bam_path, &ref_fasta, 1, 1000, Selection::Alignment);
        assert_eq!(
            alignment_candidates,
            ["on_target"].into_iter().map(str::to_string).collect(),
            "alignment mode must record only the on-target pair -- the off-target-by-position \
             pair must NOT be recorded even though its sequence k-mer-matches"
        );
        assert_eq!(
            alignment_emitted, 1,
            "alignment mode emits the unaligned-template pair directly (UnalignedPair path)"
        );
        assert_eq!(
            alignment_sink.pairs.len(),
            1,
            "alignment mode must emit exactly the unaligned-template pair"
        );
        assert_eq!(alignment_sink.pairs[0].0.id, "unaligned_kmer");
    }

    /// A single-end BAM exercising `Selection::NoAlignment`'s PRIMARY-ONLY
    /// invariant: `good_read` is a real primary read whose SEQ k-mer-matches
    /// (must be selected); `primary_noise` is a PRIMARY read whose SEQ is
    /// non-matching noise (must NOT be selected on its own) immediately
    /// followed by a SECONDARY alignment (`0x100`) of the SAME QNAME whose
    /// SEQ DOES k-mer-match; `supp_noise` is likewise a non-matching PRIMARY
    /// read followed by a SUPPLEMENTARY alignment (`0x800`) of the same
    /// QNAME whose SEQ k-mer-matches. If `scan_pass1_chunk`'s
    /// `!is_primary()` skip were missing or wrong, the secondary/
    /// supplementary records' matching sequences would each independently
    /// pass the k-mer filter and get emitted under `primary_noise`/
    /// `supp_noise`'s QNAME -- this fixture gives that regression direct,
    /// unambiguous teeth.
    fn build_no_alignment_primary_only_test_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);

        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        // good_read: primary, k-mer-matching SEQ -- must be selected.
        let good_seq = &ref_bytes[0..90];
        let mut good = bam::Record::new();
        good.set(b"good_read", Some(&CigarString(vec![Cigar::Match(90)])), good_seq, &[30u8; 90]);
        good.set_tid(0);
        good.set_pos(1100);
        good.set_mtid(-1);
        good.set_mpos(-1);
        good.set_flags(0);
        writer.write(&good).unwrap();

        // primary_noise: PRIMARY, non-matching noise SEQ -- must NOT be
        // selected on its own.
        let noise_seq1 = noise(80);
        let mut noise1 = bam::Record::new();
        noise1.set(
            b"primary_noise",
            Some(&CigarString(vec![Cigar::Match(80)])),
            &noise_seq1,
            &[30u8; 80],
        );
        noise1.set_tid(0);
        noise1.set_pos(1200);
        noise1.set_mtid(-1);
        noise1.set_mpos(-1);
        noise1.set_flags(0);
        writer.write(&noise1).unwrap();

        // Same QNAME, SECONDARY (0x100), k-mer-matching SEQ -- must be
        // SKIPPED entirely (never tested) by the primary-only invariant.
        let sec_seq = &ref_bytes[100..190];
        let mut sec = bam::Record::new();
        sec.set(b"primary_noise", Some(&CigarString(vec![Cigar::Match(90)])), sec_seq, &[30u8; 90]);
        sec.set_tid(0);
        sec.set_pos(1210);
        sec.set_mtid(-1);
        sec.set_mpos(-1);
        sec.set_flags(0x100);
        writer.write(&sec).unwrap();

        // supp_noise: PRIMARY, non-matching noise SEQ -- must NOT be
        // selected on its own.
        let noise_seq2 = noise(80);
        let mut noise2 = bam::Record::new();
        noise2.set(
            b"supp_noise",
            Some(&CigarString(vec![Cigar::Match(80)])),
            &noise_seq2,
            &[30u8; 80],
        );
        noise2.set_tid(0);
        noise2.set_pos(1300);
        noise2.set_mtid(-1);
        noise2.set_mpos(-1);
        noise2.set_flags(0);
        writer.write(&noise2).unwrap();

        // Same QNAME, SUPPLEMENTARY (0x800), k-mer-matching SEQ -- must be
        // SKIPPED entirely by the primary-only invariant.
        let supp_seq = &ref_bytes[200..280];
        let mut supp = bam::Record::new();
        supp.set(b"supp_noise", Some(&CigarString(vec![Cigar::Match(80)])), supp_seq, &[30u8; 80]);
        supp.set_tid(0);
        supp.set_pos(1310);
        supp.set_mtid(-1);
        supp.set_mpos(-1);
        supp.set_flags(0x800);
        writer.write(&supp).unwrap();

        drop(writer);
    }

    #[test]
    fn no_alignment_pass1_skips_secondary_and_supplementary_reads() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("primary_only.bam");
        build_no_alignment_primary_only_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let (candidates, pass1_emitted, sink) =
            run_pass1_chunked_test(&bam_path, &ref_fasta, 1, 1000, Selection::NoAlignment);

        assert!(candidates.is_empty(), "single-end input must not populate candidates");
        assert_eq!(
            pass1_emitted, 1,
            "only `good_read` may be emitted -- the secondary/supplementary k-mer-matching \
             records must be skipped entirely, not tested independently of their (non-matching) \
             primary record"
        );
        let ids: Vec<&str> = sink.pairs.iter().map(|(r1, _)| r1.id.as_str()).collect();
        assert_eq!(ids, vec!["good_read"], "unexpected emitted read set: {ids:?}");
    }

    // -- No-alignment cross-chunk-boundary regression: Task 1's
    // cross-boundary tests (above) only exercise `Selection::Alignment`.
    // `run_pass1_chunked` at a SMALL `chunk_size` under `Selection::
    // NoAlignment` must match a single-chunk (large `chunk_size`)
    // `Selection::NoAlignment` run exactly, proving `KmerCandidate` sites and
    // `candidates`-map accumulation thread correctly across chunks under
    // this selection mode too.

    #[test]
    fn no_alignment_pass1_matches_across_chunk_boundaries() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("no_alignment_chunked.bam");
        build_small_paired_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        // Single-chunk (large `chunk_size`) `Selection::NoAlignment` run --
        // the reference this test's chunked runs must match.
        let (reference_candidates, reference_emitted, reference_sink) =
            run_pass1_chunked_test(&bam_path, &ref_fasta, 1, 1000, Selection::NoAlignment);

        // Concrete, non-tautological expectation (not just "chunked ==
        // itself"): every primary read is k-mer-tested regardless of
        // position, so all 8 on-target QNAMEs AND the unaligned-template
        // pair's QNAME (a real k-mer match) land in `candidates`. The
        // off-target pair's `noise(80)` SEQ fails the k-mer filter, so under
        // NoAlignment it IS tested (unlike under `Selection::Alignment`,
        // which drops it by position before ever reaching a candidate-filter
        // call) but not recorded, since it fails.
        let mut expected_candidates: std::collections::HashSet<String> =
            (0..8).map(|i| format!("on_target_{i}")).collect();
        expected_candidates.insert("unaligned".to_string());
        assert_eq!(reference_candidates, expected_candidates);
        assert_eq!(reference_emitted, 0, "paired no-alignment never emits directly");
        assert!(reference_sink.pairs.is_empty());

        // `chunk_size = 4` over this ~20-record fixture crosses several
        // chunk boundaries (5 chunks); both a sequential (`threads = 1`) and
        // a parallel (`threads = 4`) evaluate step must reproduce the
        // single-chunk result exactly.
        for threads in [1usize, 4] {
            let (candidates, pass1_emitted, sink) =
                run_pass1_chunked_test(&bam_path, &ref_fasta, threads, 4, Selection::NoAlignment);
            assert_eq!(
                candidates, reference_candidates,
                "threads={threads}, chunk_size=4: no-alignment candidates diverged from the \
                 single-chunk reference across chunk boundaries"
            );
            assert_eq!(
                pass1_emitted, reference_emitted,
                "threads={threads}, chunk_size=4: no-alignment pass1_emitted diverged from the \
                 single-chunk reference"
            );
            assert_eq!(
                sink_signature(&sink),
                sink_signature(&reference_sink),
                "threads={threads}, chunk_size=4: no-alignment emitted output diverged across \
                 chunk boundaries"
            );
        }
    }
}
