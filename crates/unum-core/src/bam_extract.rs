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
use crate::extract::{
    CandidateSink, HIT_LEN_REQUIRED_PAIRED, HIT_LEN_REQUIRED_SINGLE, HIT_LEN_SAMPLE_SIZE,
    ReadRecord,
};
use crate::ref_kmer_filter::{RefKmerFilter, Scratch, is_low_complexity};
use anyhow::{Context, Result, bail};
use rayon::prelude::*;
use std::collections::{HashMap, VecDeque};

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

/// Independent per-read gene-interval overlap test for non-coordinate-ordered
/// input (grouped/name-sorted alignment). Replaces the coordinate path's
/// monotonic `tag` cursor with a self-contained scan, preserving its EXACT
/// accept condition: `chr_id == g.chr_id && start <= g.end && end > g.start`.
/// A plain linear scan (NOT a binary search): `genes` is tiny (HLA/KIR), and
/// `end` is not monotonic within a chromosome under nested intervals, so a
/// binary search keyed on `end` would violate its precondition and drop valid
/// overlaps.
pub(crate) fn read_overlaps_any_gene(
    chr_id: i32,
    start: i32,
    end: i32,
    genes: &[GeneInterval],
) -> bool {
    genes.iter().any(|g| g.chr_id == chr_id && start <= g.end && end > g.start)
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
/// Wired into [`extract_from_bam_no_alignment`], the coordinate no-alignment
/// 2-pass entry point.
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

    let pass2_emitted = run_pass2(
        alignments,
        candidates,
        abnormal_unaligned_flag,
        mate_id_len,
        Selection::Alignment,
        sink,
    )?;

    Ok(BamExtractMetrics {
        single_end: false,
        hit_len_required,
        kmer_length: filter.kmer_length(),
        pass1_emitted,
        candidates_recorded,
        pass2_emitted,
    })
}

/// Runs the coordinate/unsorted `--bam-mode no-alignment` two-pass BAM/CRAM
/// read-extraction driver: the FASTQ-pinned `hitLenRequired` setup
/// ([`compute_hit_len_required_no_alignment`], matching
/// [`crate::extract::extract_candidates_with_threads`]'s own formula so
/// `no-alignment ≡ FASTQ`), a [`Selection::NoAlignment`] pass 1 (every
/// PRIMARY read k-mer-tested on its own sequence, position/gene-interval/
/// `tag` entirely bypassed -- see [`Selection`]'s doc comment), and a
/// gate-bypassed pass 2 ([`run_pass2`] never applies its
/// `!is_template_aligned()` gate under [`Selection::NoAlignment`], so every
/// recorded candidate's mates are fetched regardless of alignment position).
///
/// Unlike [`extract_from_bam_with_threads`], no `genes`/coordinate FASTA is
/// needed here -- selection is purely sequence-driven. `ref_seq_similarity`
/// (the `-s` CLI flag) is REQUIRED so the filter's similarity matches the
/// FASTQ path's value: `is_good_candidate_with_scratch` gates on
/// `filter.ref_seq_similarity`, and the FASTQ path sets it from
/// `args.similarity` (`extract.rs:373`) -- omitting this call here would
/// silently break the `no-alignment ≡ FASTQ` equivalence for any non-default
/// `-s`.
///
/// Same seekable-file assumption as [`extract_from_bam_with_threads`]:
/// `alignments` must be freshly opened when this is called. This function
/// rewinds `alignments` itself: once after `general_info` (mirroring
/// `BamExtractor.cpp:573-574`'s `GetGeneralInfo(true); Rewind();` pattern)
/// and again after [`crate::alignments::Alignments::sample_read1_len_sum`]
/// (whose doc comment notes it does NOT rewind, and which itself advances
/// the shared `total_read_cnt` counter via `next()` -- rewinding here keeps
/// that counter clean for any future caller, rather than leaving it
/// contaminated by the sample).
///
/// Single-end input skips pass 2 entirely (mirroring
/// [`extract_from_bam_with_threads`]): pass 1's
/// [`Pass1Site::KmerCandidate`] emits single-end reads directly, so
/// `candidates` is never populated for single-end input.
///
/// # Errors
///
/// See [`extract_from_bam_with_threads`].
pub fn extract_from_bam_no_alignment(
    alignments: &mut Alignments,
    filter: &mut RefKmerFilter,
    ref_seq_similarity: f64,
    mate_id_len: i32,
    threads: usize,
    sink: &mut impl CandidateSink,
) -> Result<BamExtractMetrics> {
    // Setup: general_info -> rewind -> sample_read1_len_sum -> rewind, to
    // keep `total_read_cnt` clean before pass 1 (see doc comment above).
    let info = alignments.general_info(true).context("computing general_info (no-alignment)")?;
    alignments.rewind().context("rewinding after general_info (no-alignment)")?;

    let single_end = !info.mate_paired;
    let has_mate = info.mate_paired;

    let (len_sum, count) = alignments
        .sample_read1_len_sum(HIT_LEN_SAMPLE_SIZE)
        .context("sampling read-1 lengths (no-alignment)")?;
    alignments.rewind().context("rewinding after sampling read-1 lengths (no-alignment)")?;

    // FASTQ-pinned hitLenRequired + similarity (the equivalence
    // prerequisites), then InferKmerLength / conditional UpdateKmerLength --
    // mirrors the FASTQ setup (extract.rs ~370-386) exactly.
    let mut hit_len_required = compute_hit_len_required_no_alignment(len_sum, count, has_mate);
    filter.set_hit_len_required(hit_len_required);
    filter.set_ref_seq_similarity(ref_seq_similarity);

    let inferred = filter.infer_kmer_length();
    if inferred > filter.kmer_length() {
        filter.update_kmer_length(inferred);
        if inferred > usize::try_from(hit_len_required).unwrap_or(0) {
            hit_len_required = i32::try_from(inferred).unwrap_or(i32::MAX);
            filter.set_hit_len_required(hit_len_required);
        }
    }

    let (candidates, pass1_emitted) = run_pass1_with_threads(
        alignments,
        filter,
        &[],
        single_end,
        // `abnormal_unaligned_flag` is dead under `Selection::NoAlignment`
        // (the unaligned-pair path is alignment-only), hence `false`.
        false,
        mate_id_len,
        threads,
        sink,
        Selection::NoAlignment,
    )?;

    let candidates_recorded = u64::try_from(candidates.len()).unwrap_or(u64::MAX);
    if single_end {
        // Single-end emits directly in pass 1; skip pass 2 entirely, same as
        // `extract_from_bam_with_threads`.
        return Ok(BamExtractMetrics {
            single_end: true,
            hit_len_required,
            kmer_length: filter.kmer_length(),
            pass1_emitted,
            candidates_recorded: 0,
            pass2_emitted: 0,
        });
    }

    alignments.rewind().context("rewinding before pass 2 (no-alignment)")?;
    let pass2_emitted =
        run_pass2(alignments, candidates, false, mate_id_len, Selection::NoAlignment, sink)?;

    Ok(BamExtractMetrics {
        single_end: false,
        hit_len_required,
        kmer_length: filter.kmer_length(),
        pass1_emitted,
        candidates_recorded,
        pass2_emitted,
    })
}

/// One buffered/live PRIMARY record captured by
/// [`extract_from_bam_no_alignment_grouped_with_head_limit`]'s head-then-live
/// group loop: enough of a record's state to re-derive both the setup
/// statistics (read-1 length sum, `0x1`-set count) and the per-group k-mer/
/// pairing decision, without holding on to the live [`Alignments`] cursor
/// itself (which a `VecDeque` of these can outlive across the head/live
/// boundary -- see [`next_buffered_or_live`]).
struct GroupedRecord {
    /// The QNAME with `mate_id_len` trimming already applied ([`trim_name`]),
    /// applied uniformly to every record regardless of its own `0x1` bit --
    /// this is the key the group loop compares consecutive records by, so a
    /// mate pair whose raw QNAMEs differ only by a trailing `/1`/`/2` (or
    /// whatever `mate_id_len` strips) still groups together. Also the id the
    /// PAIRED-arm emits (both mates share a QNAME, so the trimmed form is the
    /// mate-1 id every downstream sink expects -- matching the coordinate
    /// path's paired 2-pass name-map, which keys on the trimmed name too).
    trimmed_name: String,
    /// The RAW, UNtrimmed QNAME ([`Alignments::read_id`]). Emitted verbatim
    /// as the `ReadRecord.id` of a genuine single-end (0x1-UNSET) lone read,
    /// to match the coordinate no-alignment path's single-end emits
    /// (`Pass1Site::KmerCandidate`/`SingleEndNotAligned`/`SingleEndOnTarget`,
    /// which all carry `read_id()` UNTRIMMED for single-end input): under a
    /// positive `mate_id_len` or a single-end name bearing a `/1`/`/2`
    /// suffix, trimming would strip characters the coordinate/FASTQ
    /// reference keeps, diverging from the very path Task 7 proves this one
    /// equivalent to. Only the genuine single-end lone emit uses this field;
    /// the paired arm uses [`Self::trimmed_name`], and an orphan (`0x1`-SET,
    /// mate absent) is DROPPED rather than emitted (see
    /// [`flush_grouped_candidate`]).
    raw_name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    /// The `0x40` FLAG bit ([`Alignments::is_first_mate`]) -- used to order a
    /// 2-member group's `(mate1, mate2)` emission.
    is_first_mate: bool,
    /// The `0x1` FLAG bit ([`Alignments::is_paired`]) -- used only for the
    /// setup-phase majority-vote (`single_end` derivation); a group's
    /// pairing-rule dispatch itself is driven entirely by how many members
    /// it accumulates (0/1/2), not by this bit.
    is_paired: bool,
}

/// Runs the grouped/name-sorted (`GO:query`/`SO:queryname`) `--bam-mode
/// no-alignment` ONE-PASS BAM/CRAM read-extraction driver: unlike
/// [`extract_from_bam_no_alignment`] (coordinate/unsorted input, two passes
/// plus a `rewind()` between them), a grouped/name-sorted BAM keeps a
/// template's mate(s) adjacent in file order, so a SINGLE streaming pass can
/// both k-mer-test each primary read and reunite its QNAME group --
/// `O(group + candidates)` memory (never the whole-file `candidates` name
/// map the coordinate two-pass needs), and critically **stdin-capable**:
/// this function never calls [`Alignments::rewind`] or
/// [`Alignments::general_info`], so `alignments` may be backed by
/// [`Alignments::from_stdin`].
///
/// # Setup without `general_info` (B2)
///
/// The coordinate/unsorted entry points derive `single_end` from
/// `alignments.general_info(true).mate_paired`, which requires a `rewind()`
/// afterward -- impossible on a non-seekable stdin source. Instead, this
/// function buffers a bounded HEAD of up to [`HIT_LEN_SAMPLE_SIZE`] raw
/// PRIMARY records (see [`GroupedRecord`]) and derives BOTH the FASTQ-pinned
/// `hitLenRequired` (read-1 length sum/count -- the same rule
/// [`crate::alignments::Alignments::sample_read1_len_sum`] uses: read-1 is
/// the first mate of a paired-flag template, or any record of an unpaired
/// one) AND `single_end` from that SAME head, via the identical
/// `has_mate_cnt >= total / 2` majority-vote rule
/// [`crate::alignments::Alignments::general_info`] uses internally. Returns
/// `single_end` (alongside the sink `make_sink` was used to create -- see the
/// "Sink FACTORY" doc section below) so a caller (the dispatcher) can pick
/// the output filename from it directly, without ever calling
/// `general_info`.
///
/// # Buffering raw records, not pairs (M2)
///
/// The head buffers RAW records, not pre-paired templates -- unlike
/// `extract::sample_head` (FASTQ positional pairing; private, and typed on a
/// different reader), a grouped BAM pairs by QNAME ADJACENCY, and a QNAME
/// group can straddle the boundary between the sampled head and the live
/// stream that follows it. This function's group-accumulation loop pulls
/// from the head deque first, then the live [`Alignments`] cursor
/// ([`next_buffered_or_live`]), as ONE continuous stream -- so a group split
/// across that boundary (its first member sampled into the head, its second
/// member read live) is reunited exactly as if no head/live split existed at
/// all (see `grouped_no_alignment_reunites_group_straddling_head_boundary`).
///
/// # Pairing rules
///
/// Primary reads only (`!is_primary()` skipped, both while sampling the head
/// and while reading live -- matching [`extract_from_bam_no_alignment`]'s own
/// primary-only invariant, so a secondary/supplementary alignment of an
/// already-grouped QNAME can never masquerade as a second mate or a
/// duplicate single). Records are grouped by `trim_name`-trimmed QNAME
/// adjacency; a group is flushed (k-mer-tested and possibly emitted) the
/// moment a differently-named record is encountered, and once more at EOF
/// for whatever group is still open ([`flush_grouped_candidate`]):
/// - A 2-member group emits a PAIR iff EITHER member passes the k-mer gate
///   (`filter.is_good_candidate_with_scratch`, OR-rescue -- matching
///   [`Pass1Site::KmerCandidate`]'s own per-read-independent, OR-across-the-
///   template semantics in the coordinate no-alignment path, NOT
///   [`Pass1Site::UnalignedPair`]'s stricter "neither read low-complexity"
///   AND-gate), ordered `(mate1, mate2)` by each member's `is_first_mate` bit.
/// - A 1-member group with the `0x1` bit UNSET (a genuine single-end record)
///   emits a LONE read iff it passes, carrying its RAW (untrimmed) id --
///   matching the coordinate no-alignment path's single-end emits, so a
///   positive `mate_id_len`/`/1`-suffixed single-end name does not diverge
///   from the `≡FASTQ`/`≡coordinate` reference (see
///   [`GroupedRecord::raw_name`]).
/// - A 1-member group with `0x1` SET (an ORPHAN -- this template's other mate
///   never arrived, e.g. filtered upstream or genuinely missing from the file)
///   is DROPPED, exactly as the coordinate path drops it: a paired candidate
///   whose mate is absent is recorded in pass 1 but never emitted by pass 2
///   (which fetches both primary mates). Emitting it lone would break
///   set-equality with the coordinate path AND hand the fused `run` path an
///   unequal read-1/read-2 count.
/// - A group that grows past 2 primary members aborts with an error hint (a
///   QNAME group larger than 2 means the input is not actually
///   grouped/name-sorted the way this one-pass entry point requires).
///
/// `threads` is accepted for signature parity with the sibling
/// `extract_from_bam_no_alignment*` entry points but unused: a group here has
/// at most 2 members, so there is no useful unit of parallel work (unlike
/// [`run_pass1_with_threads`]'s whole-chunk-wide candidate-decision
/// parallelism).
///
/// # Sink FACTORY, not a concrete sink (the B2 chicken-and-egg)
///
/// A concrete sink (e.g. the CLI's `FastqFileSink`) often needs `single_end`
/// to pick its own output naming -- but `single_end` here is derived from
/// this function's OWN buffered head (see "Setup without `general_info`"
/// above), not knowable before this function runs. So this takes `make_sink`,
/// a `FnOnce(single_end) -> Result<S>` factory, and calls it itself the
/// MOMENT `single_end` is known (right after the head is buffered and the
/// filter is configured, before the group loop emits a single candidate) --
/// never before, and never via a `general_info` pre-sample (impossible on
/// stdin). The created sink is returned to the caller (rather than dropped
/// internally) so it can be flushed and its I/O errors observed, and so
/// tests can inspect what it collected.
///
/// # Errors
///
/// Returns an error if: `make_sink` itself fails; a QNAME group exceeds 2
/// primary records; or the underlying BAM read fails (a genuine parse error,
/// distinct from a clean EOF).
pub fn extract_from_bam_no_alignment_grouped<S, F>(
    alignments: &mut Alignments,
    filter: &mut RefKmerFilter,
    ref_seq_similarity: f64,
    mate_id_len: i32,
    threads: usize,
    make_sink: F,
) -> Result<(BamExtractMetrics, bool, S)>
where
    S: CandidateSink,
    F: FnOnce(bool) -> Result<S>,
{
    let _ = threads;
    extract_from_bam_no_alignment_grouped_with_head_limit(
        alignments,
        filter,
        ref_seq_similarity,
        mate_id_len,
        HIT_LEN_SAMPLE_SIZE,
        make_sink,
    )
}

/// [`extract_from_bam_no_alignment_grouped`]'s implementation, with the
/// sampled-head bound taken as an explicit parameter rather than the fixed
/// [`HIT_LEN_SAMPLE_SIZE`] constant -- so `#[cfg(test)]` can force a tiny head
/// (e.g. `1`) to directly exercise a QNAME group straddling the head/live
/// boundary (see the public wrapper's "Buffering raw records" doc section),
/// which no fixture small enough for a unit test could otherwise reach at the
/// real 1000-record bound.
fn extract_from_bam_no_alignment_grouped_with_head_limit<S, F>(
    alignments: &mut Alignments,
    filter: &mut RefKmerFilter,
    ref_seq_similarity: f64,
    mate_id_len: i32,
    head_limit: usize,
    make_sink: F,
) -> Result<(BamExtractMetrics, bool, S)>
where
    S: CandidateSink,
    F: FnOnce(bool) -> Result<S>,
{
    // -- (1) Bounded head of RAW PRIMARY records (no rewind -> stdin-safe).
    let mut head: VecDeque<GroupedRecord> = VecDeque::new();
    let mut read1_len_sum: i64 = 0;
    let mut read1_count: usize = 0;
    let mut has_mate_cnt: u64 = 0;
    let mut head_primary_cnt: u64 = 0;

    while head.len() < head_limit
        && alignments.next().context("reading BAM head (grouped no-alignment)")?
    {
        if !alignments.is_primary() {
            continue;
        }
        let record = read_grouped_record(alignments, mate_id_len);
        head_primary_cnt += 1;
        if record.is_paired {
            has_mate_cnt += 1;
        }
        let is_read1 = !record.is_paired || record.is_first_mate;
        if is_read1 {
            read1_len_sum += i64::try_from(record.seq.len()).unwrap_or(i64::MAX);
            read1_count += 1;
        }
        head.push_back(record);
    }

    // -- (2) Setup from the head alone (B2: never general_info).
    // `!(has_mate_cnt >= total / 2)`, simplified per clippy::nonminimal_bool
    // (equivalent for these integer counts): `single_end` is the negation of
    // `Alignments::general_info`'s own `mate_paired` majority-vote rule.
    let single_end = has_mate_cnt < head_primary_cnt / 2;
    let mut hit_len_required =
        compute_hit_len_required_no_alignment(read1_len_sum, read1_count, !single_end);
    filter.set_hit_len_required(hit_len_required);
    filter.set_ref_seq_similarity(ref_seq_similarity); // M1: required for ≡FASTQ.

    let inferred = filter.infer_kmer_length();
    if inferred > filter.kmer_length() {
        filter.update_kmer_length(inferred);
        if inferred > usize::try_from(hit_len_required).unwrap_or(0) {
            hit_len_required = i32::try_from(inferred).unwrap_or(i32::MAX);
            filter.set_hit_len_required(hit_len_required);
        }
    }

    // `single_end` is now known -- create the sink NOW, via the factory, and
    // BEFORE the group loop below emits a single candidate (the B2 fix: a
    // concrete sink like `FastqFileSink` needs `single_end` to pick its own
    // output naming, and that could not have been known any earlier than
    // this point on a non-seekable stdin source).
    let mut sink = make_sink(single_end)?;

    // -- (3) ONE continuous group-accumulation loop over head ++ live.
    let mut scratch = Scratch::default();
    let mut current_group: Vec<GroupedRecord> = Vec::new();
    let mut emitted: u64 = 0;

    while let Some(record) = next_buffered_or_live(&mut head, alignments, mate_id_len)? {
        if let Some(first) = current_group.first() {
            if first.trimmed_name == record.trimmed_name {
                if current_group.len() >= 2 {
                    bail!(
                        "more than two primary records share QNAME {:?} in grouped/name-sorted \
                         no-alignment mode; hint: the input's @HD GO:query/SO:queryname claim may \
                         not hold (mates are not actually adjacent) -- unum's grouped one-pass \
                         requires a genuinely grouped/name-sorted BAM",
                        record.trimmed_name
                    );
                }
                current_group.push(record);
                continue;
            }
            flush_grouped_candidate(&current_group, filter, &mut scratch, &mut sink, &mut emitted)?;
            current_group.clear();
        }
        current_group.push(record);
    }
    flush_grouped_candidate(&current_group, filter, &mut scratch, &mut sink, &mut emitted)?;

    let metrics = BamExtractMetrics {
        single_end,
        hit_len_required,
        kmer_length: filter.kmer_length(),
        pass1_emitted: emitted,
        candidates_recorded: 0,
        pass2_emitted: 0,
    };
    Ok((metrics, single_end, sink))
}

/// Reads the CURRENT record's [`GroupedRecord`] snapshot. Caller must have
/// already called `alignments.next()` and confirmed `is_primary()`.
fn read_grouped_record(alignments: &Alignments, mate_id_len: i32) -> GroupedRecord {
    let raw_name = alignments.read_id();
    GroupedRecord {
        trimmed_name: trim_name(&raw_name, mate_id_len),
        raw_name,
        seq: alignments.read_seq(),
        qual: alignments.qual(),
        is_first_mate: alignments.is_first_mate(),
        is_paired: alignments.is_paired(),
    }
}

/// Pulls the next PRIMARY record for
/// [`extract_from_bam_no_alignment_grouped_with_head_limit`]'s group loop:
/// drains `head` first (in original encounter order), then reads live from
/// `alignments` -- skipping non-primary records exactly as the head-buffering
/// loop does -- once `head` is empty. Returns `None` at EOF. This is what
/// makes a QNAME group split across the head/live boundary transparent to the
/// caller: from the group loop's point of view, the head is just a pushback
/// buffer in front of the live cursor, not a separately-paired structure (see
/// [`extract_from_bam_no_alignment_grouped_with_head_limit`]'s "Buffering raw
/// records" doc section).
///
/// # Errors
///
/// Returns an error if the underlying BAM read fails (a genuine parse error,
/// distinct from a clean EOF).
fn next_buffered_or_live(
    head: &mut VecDeque<GroupedRecord>,
    alignments: &mut Alignments,
    mate_id_len: i32,
) -> Result<Option<GroupedRecord>> {
    if let Some(record) = head.pop_front() {
        return Ok(Some(record));
    }
    while alignments.next().context("reading next BAM record (grouped no-alignment)")? {
        if !alignments.is_primary() {
            continue;
        }
        return Ok(Some(read_grouped_record(alignments, mate_id_len)));
    }
    Ok(None)
}

/// Resolves and emits (if it passes) one flushed QNAME group -- see
/// [`extract_from_bam_no_alignment_grouped`]'s "Pairing rules" doc section.
/// `group` is expected to have length 0 (a no-op, e.g. the very first call
/// before any record has been accumulated), 1, or 2; the caller's own
/// larger-than-2-member abort check in the group loop guarantees a 3rd
/// member is never pushed, so the `_` arm below is unreachable in practice
/// and exists only as a defensive fallback.
///
/// # Errors
///
/// Returns an error if `sink.emit_pair` fails, or (defensively) if `group`
/// somehow has more than 2 members.
fn flush_grouped_candidate(
    group: &[GroupedRecord],
    filter: &RefKmerFilter,
    scratch: &mut Scratch,
    sink: &mut impl CandidateSink,
    emitted: &mut u64,
) -> Result<()> {
    match group {
        [] => Ok(()),
        [lone] => {
            // Only a GENUINE single-end record (0x1 UNSET) emits lone -- with
            // its RAW, untrimmed id, matching the coordinate no-alignment path's
            // single-end emits (which carry `read_id()` untrimmed -- see
            // `GroupedRecord::raw_name`'s doc comment); this is the
            // `≡FASTQ`/`≡coordinate` (Task 7) requirement.
            //
            // An ORPHAN (0x1 SET, this template's mate never arrived in the
            // file) is DROPPED, matching the coordinate path exactly: a paired
            // candidate whose mate is absent is recorded in pass 1 but never
            // emitted by pass 2 (which fetches BOTH primary mates). Emitting it
            // lone would (a) break set-equality with the coordinate path on
            // orphan-containing input and (b) hand the fused `run` path an
            // unequal read-1/read-2 count, which `run_with_candidate_reads`
            // rejects. So the drop is required for both correctness and the
            // fused genotyper's mate-count invariant.
            if !lone.is_paired && filter.is_good_candidate_with_scratch(&lone.seq, scratch) {
                let read = ReadRecord {
                    id: lone.raw_name.clone(),
                    seq: lone.seq.clone(),
                    qual: Some(lone.qual.clone()),
                };
                sink.emit_pair(&read, None)?;
                *emitted += 1;
            }
            Ok(())
        }
        [a, b] => {
            // Order by `is_first_mate`; if the group is malformed (neither
            // or both members carry `0x40`), fall back to encounter order --
            // an arbitrary but deterministic choice for input that does not
            // actually distinguish its two mates.
            let (mate1, mate2) = if a.is_first_mate && !b.is_first_mate {
                (a, b)
            } else if b.is_first_mate && !a.is_first_mate {
                (b, a)
            } else {
                (a, b)
            };
            // OR-rescue: either mate passing the k-mer gate independently is
            // enough (matches `Pass1Site::KmerCandidate`'s per-read decision,
            // OR-combined via `candidates`-map insertion in the coordinate
            // no-alignment path -- see this function's caller's doc comment).
            let passes = filter.is_good_candidate_with_scratch(&mate1.seq, scratch)
                || filter.is_good_candidate_with_scratch(&mate2.seq, scratch);
            if passes {
                // Both output records carry mate1's id, matching every other
                // pair-emission site in this module (e.g.
                // `Pass1Site::UnalignedPair`'s apply arm) and
                // `InMemoryCandidateSink`'s documented byte-identity
                // requirement with the file-writing sink.
                let r1 = ReadRecord {
                    id: mate1.trimmed_name.clone(),
                    seq: mate1.seq.clone(),
                    qual: Some(mate1.qual.clone()),
                };
                let r2 = ReadRecord {
                    id: mate1.trimmed_name.clone(),
                    seq: mate2.seq.clone(),
                    qual: Some(mate2.qual.clone()),
                };
                sink.emit_pair(&r1, Some(&r2))?;
                *emitted += 1;
            }
            Ok(())
        }
        _ => bail!(
            "internal error: grouped no-alignment flush called with a {}-member group (expected \
             0, 1, or 2 -- the caller's >2-member abort check should have already fired)",
            group.len()
        ),
    }
}

/// One buffered/live record captured by [`extract_from_bam_alignment_grouped`]'s
/// group loop -- the ALIGNMENT analogue of [`GroupedRecord`]. Unlike the
/// no-alignment one, this captures EVERY record of a QNAME group (primary +
/// secondary + supplementary, no `is_primary()` skip), plus the per-record
/// alignment STATE ([`scan_pass1_chunk`]'s `Selection::Alignment` branch reads
/// off each record) needed to reproduce that path's status-conditioned
/// candidacy without re-reading the BAM: a secondary alignment overlapping a
/// gene makes the whole template a candidate, so its state must be carried,
/// even though only the PRIMARY sequences are ever emitted.
///
/// The multiple `bool` fields are independent per-record BAM FLAG/alignment
/// facts (a captured record snapshot), not a hidden state machine -- so
/// `clippy::struct_excessive_bools` does not apply; folding them into enums
/// would only obscure their direct correspondence to `scan_pass1_chunk`'s
/// per-record reads.
#[allow(clippy::struct_excessive_bools)]
struct AlignmentGroupedRecord {
    /// QNAME with `mate_id_len` trimming applied ([`trim_name`]) -- the group
    /// key, and the id a 2-primary pair emits (see [`GroupedRecord::trimmed_name`]).
    trimmed_name: String,
    /// Raw, untrimmed QNAME ([`Alignments::read_id`]) -- the id a genuine
    /// single-end (`0x1`-unset) lone read emits (see [`GroupedRecord::raw_name`]).
    raw_name: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    /// `0x40` (READ1) -- orders a 2-primary group's `(mate1, mate2)` emission.
    is_first_mate: bool,
    /// `0x1` (paired) -- a lone record's `0x1`-unset (genuine single-end, emit
    /// untrimmed id) vs `0x1`-set (orphan, DROPPED to match the coordinate
    /// path) dispatch, and the non-seekable head's `single_end` majority vote.
    is_paired: bool,
    /// `!(0x900)` -- only primary records are emitted; secondaries/
    /// supplementaries contribute to candidacy but never to output.
    is_primary: bool,
    /// [`Alignments::is_template_aligned`] -- the first branch key
    /// (`!template_aligned || chrom_alt`) of `scan_pass1_chunk`.
    is_template_aligned: bool,
    /// [`Alignments::is_aligned`] -- distinguishes an unaligned mate of an
    /// aligned pair (`template_aligned && !is_aligned`, not a candidate
    /// trigger) from a main-chrom on-target read.
    is_aligned: bool,
    /// `is_aligned && valid_alternative_chrom(chrom_name(chrom_id))` --
    /// captured here (guarded by `is_aligned` short-circuit so `chrom_name` is
    /// never called on an unmapped `tid == -1`), matching `scan_pass1_chunk`'s
    /// own `chrom_alt` computation.
    chrom_alt: bool,
    /// Reference id ([`Alignments::chrom_id`]) -- the overlap test's `chr_id`.
    chr_id: i32,
    /// `segments.first().a` / `segments.last().b`, meaningful only when
    /// [`Self::has_segments`]; the interval fed to [`read_overlaps_any_gene`].
    start: i32,
    end: i32,
    /// Whether the record has any CIGAR-derived segment. `scan_pass1_chunk`
    /// (`bam_extract.rs`'s `let Some(first_seg) = segments.first() else {
    /// continue }`) drops an aligned read with no segments before the overlap
    /// test; the on-target candidacy branch reproduces that by requiring this.
    has_segments: bool,
}

/// Reads the CURRENT record's [`AlignmentGroupedRecord`] snapshot. Caller must
/// have already called `alignments.next()`. Captures records regardless of
/// `is_primary()` (unlike the no-alignment [`read_grouped_record`]): the
/// alignment path classifies every member of a QNAME group.
fn read_alignment_grouped_record(
    alignments: &Alignments,
    mate_id_len: i32,
) -> AlignmentGroupedRecord {
    let raw_name = alignments.read_id();
    let is_aligned = alignments.is_aligned();
    // `is_aligned &&` short-circuits so `chrom_name` is never called on an
    // unmapped `tid == -1` (which would panic) -- exactly `scan_pass1_chunk`.
    let chrom_alt =
        is_aligned && valid_alternative_chrom(&alignments.chrom_name(alignments.chrom_id()));
    let segments = alignments.segments();
    let has_segments = !segments.is_empty();
    // `a`/`b` are `i64` reference coordinates; the same truncating cast
    // `scan_pass1_chunk` applies. Only consulted when `has_segments` (aligned).
    #[allow(clippy::cast_possible_truncation)]
    let (start, end) = match (segments.first(), segments.last()) {
        (Some(first), Some(last)) => (first.a as i32, last.b as i32),
        _ => (0, 0),
    };
    AlignmentGroupedRecord {
        trimmed_name: trim_name(&raw_name, mate_id_len),
        raw_name,
        seq: alignments.read_seq(),
        qual: alignments.qual(),
        is_first_mate: alignments.is_first_mate(),
        is_paired: alignments.is_paired(),
        is_primary: alignments.is_primary(),
        is_template_aligned: alignments.is_template_aligned(),
        is_aligned,
        chrom_alt,
        chr_id: alignments.chrom_id(),
        start,
        end,
        has_segments,
    }
}

/// Pulls the next record (primary OR non-primary) for
/// [`extract_from_bam_alignment_grouped`]'s group loop: drains `head` first (in
/// encounter order), then reads live from `alignments`. Unlike the
/// no-alignment [`next_buffered_or_live`], it does NOT skip `!is_primary()`
/// records -- the alignment path needs the whole QNAME group. `None` at EOF.
///
/// # Errors
///
/// Returns an error if the underlying BAM read fails (a genuine parse error,
/// distinct from a clean EOF).
fn next_alignment_buffered_or_live(
    head: &mut VecDeque<AlignmentGroupedRecord>,
    alignments: &mut Alignments,
    mate_id_len: i32,
) -> Result<Option<AlignmentGroupedRecord>> {
    if let Some(record) = head.pop_front() {
        return Ok(Some(record));
    }
    if alignments.next().context("reading next BAM record (grouped alignment)")? {
        return Ok(Some(read_alignment_grouped_record(alignments, mate_id_len)));
    }
    Ok(None)
}

/// The per-record candidacy decision for the NON-joint case, reproducing
/// [`scan_pass1_chunk`]'s `Selection::Alignment` status-conditioned branches
/// EXACTLY (`bam_extract.rs` lines ~1305-1405), for ONE record of a QNAME
/// group. A group (outside the genuine-unaligned-pair joint case) is a
/// candidate iff ANY of its records returns `true` here -- the set-union that
/// mirrors the coordinate path's `candidates.entry(trimmed_name).or_default()`
/// idempotent insertion across every record of a template.
///
/// The three branches, in `scan_pass1_chunk` order:
/// 1. `!is_template_aligned || chrom_alt` (ALT-contig, or `-u`/abnormal /
///    single-end unaligned reads that did NOT take the joint path): the
///    per-read predicate `!is_low_complexity(seq) && filter.good(seq)` --
///    mirrors `Pass1Site::PairedNotAlignedCandidate` / `SingleEndNotAligned`.
/// 2. `is_template_aligned && !is_aligned` (the unaligned mate of an aligned
///    pair): `false` -- `scan_pass1_chunk` `continue`s it (never a candidate
///    trigger), though its primary seq is still EMITTED if the template is a
///    candidate via another record.
/// 3. `is_template_aligned && is_aligned && !chrom_alt` (main-chrom on-target):
///    `has_segments && !is_low_complexity(seq) && read_overlaps_any_gene(...)`
///    -- position/gene-overlap only, with NO seed/k-mer test (matches
///    `Pass1Site::PairedOnTargetCandidate`/`SingleEndOnTarget`, which resolve
///    to `true` after only the low-complexity pre-check). The `has_segments`
///    guard reproduces `scan_pass1_chunk`'s `segments.first()` early-`continue`.
fn per_record_alignment_candidate(
    record: &AlignmentGroupedRecord,
    filter: &RefKmerFilter,
    genes: &[GeneInterval],
    scratch: &mut Scratch,
) -> bool {
    if !record.is_template_aligned || record.chrom_alt {
        return !is_low_complexity(&record.seq)
            && filter.is_good_candidate_with_scratch(&record.seq, scratch);
    }
    if !record.is_aligned {
        return false;
    }
    record.has_segments
        && !is_low_complexity(&record.seq)
        && read_overlaps_any_gene(record.chr_id, record.start, record.end, genes)
}

/// The COORDINATE sort key used to pick which of a single-end name's EMITTABLE
/// records the coordinate path would reach first (and therefore emit) --
/// ascending `(chrom_id, start)`, with UNALIGNED records (`tid == -1`) sorting
/// AFTER every aligned one (a coordinate-sorted BAM places unmapped records
/// last). The leading rank byte (`0` aligned / `1` unaligned) enforces that
/// unaligned-last ordering independently of the (negative) unaligned
/// `chrom_id`, so the raw `chr_id`/`start` never have to encode +infinity.
fn single_end_coordinate_key(record: &AlignmentGroupedRecord) -> (u8, i32, i32) {
    let aligned_rank = u8::from(!record.is_aligned);
    (aligned_rank, record.chr_id, record.start)
}

/// Resolves and emits (if it passes) one flushed QNAME group for
/// [`extract_from_bam_alignment_grouped`]. Candidacy is status-conditioned:
///
/// - **Genuine unaligned pair (joint):** exactly 2 primaries, both
///   `!is_template_aligned`, with `!single_end && !abnormal_unaligned_flag` --
///   the JOINT predicate `!lc(s1) && !lc(s2) && (good(s1) || good(s2))`,
///   mirroring [`Pass1Site::UnalignedPair`]. (Unmapped reads carry no
///   secondary alignments, so "2 primaries both unmapped" already implies a
///   2-member group; the explicit primary count keeps the rule robust.)
/// - **Otherwise:** the set-union of [`per_record_alignment_candidate`] over
///   ALL records (primary + secondary + supplementary). With `-u`, an
///   unaligned pair falls here (per-read), NOT the joint branch.
///
/// Emission is status-conditioned on `single_end`:
///
/// - **Single-end:** emit exactly ONE record's seq/qual -- the FIRST in
///   COORDINATE order (via [`single_end_coordinate_key`]) among the group's
///   EMITTABLE records ([`per_record_alignment_candidate`]) -- under its
///   UNTRIMMED id. This may be a NON-PRIMARY (e.g. an on-target `0x800`
///   supplementary), matching the coordinate path's `used_name`-deduped
///   "first emittable record in the scan wins" selection; emitting the
///   primary unconditionally would break set-equality when the primary is
///   off-target but a supplementary is on-target.
/// - **Paired:** PRIMARIES ONLY, with the same id/pairing rules as
///   [`flush_grouped_candidate`]: a 2-primary group emits `(mate1, mate2)`
///   ordered by `is_first_mate` on the trimmed name; a 1-primary `0x1`-unset
///   group emits a lone read with the UNTRIMMED id; a 1-primary `0x1`-set
///   orphan is DROPPED (matching the coordinate path, which emits only
///   complete pairs, and keeping the fused `run` path's mate counts equal).
///
/// # Errors
///
/// Returns an error if a QNAME group has more than 2 PRIMARY records (a
/// malformed grouped/name-sorted claim), or if `sink.emit_pair` fails.
#[allow(clippy::too_many_arguments)]
fn flush_alignment_grouped_candidate(
    group: &[AlignmentGroupedRecord],
    filter: &RefKmerFilter,
    genes: &[GeneInterval],
    single_end: bool,
    abnormal_unaligned_flag: bool,
    scratch: &mut Scratch,
    sink: &mut impl CandidateSink,
    emitted: &mut u64,
) -> Result<()> {
    if group.is_empty() {
        return Ok(());
    }

    let primaries: Vec<&AlignmentGroupedRecord> =
        group.iter().filter(|record| record.is_primary).collect();
    if primaries.len() > 2 {
        bail!(
            "more than two primary records share QNAME {:?} in grouped/name-sorted alignment \
             mode; hint: the input's @HD GO:query/SO:queryname claim may not hold (mates are not \
             actually adjacent) -- unum's grouped one-pass requires a genuinely grouped/\
             name-sorted BAM",
            group[0].trimmed_name
        );
    }

    let is_candidate = if !single_end
        && !abnormal_unaligned_flag
        && primaries.len() == 2
        && !primaries[0].is_template_aligned
        && !primaries[1].is_template_aligned
    {
        // JOINT unaligned-pair predicate (mirrors Pass1Site::UnalignedPair):
        // both mates non-low-complexity AND at least one a good candidate.
        let seq1 = &primaries[0].seq;
        let seq2 = &primaries[1].seq;
        (!is_low_complexity(seq1) && !is_low_complexity(seq2))
            && (filter.is_good_candidate_with_scratch(seq1, scratch)
                || filter.is_good_candidate_with_scratch(seq2, scratch))
    } else {
        // OR-union of the per-record decision over EVERY record of the group.
        group.iter().any(|record| per_record_alignment_candidate(record, filter, genes, scratch))
    };

    if !is_candidate {
        return Ok(());
    }

    // -- SINGLE-END emission (the faithful selection). The coordinate path
    // emits, per single-end name, the seq/qual of the record it reaches FIRST
    // in COORDINATE order among that name's EMITTABLE records -- which may be a
    // NON-PRIMARY: an on-target `0x800` supplementary retains real SEQ, so a
    // group whose primary is off-target but whose supplementary is on-target
    // emits the SUPPLEMENTARY's seq, not the primary's. (`scan_pass1_chunk` does
    // NOT filter `!is_primary()` on the alignment path, and its shared
    // `used_name` dedup lets the first emittable aligned record in coordinate
    // order win, blocking the rest.) Emitting the primary here (the previous
    // behavior) broke set-equality with the coordinate path on exactly that
    // shape. Reproduce the coordinate selection: among all records whose
    // per-record predicate holds -- identical to `per_record_alignment_candidate`,
    // the coordinate single-end EMITTABLE condition (on-target = `!lc && overlap`;
    // not-template-aligned/alt = `!lc && good`) -- pick the minimum by
    // `single_end_coordinate_key` (coordinate order, unaligned last). `min_by_key`
    // returns the FIRST of equal-key records, i.e. group/file order; two records
    // of one single-end read sharing an exact `(tid, pos)` does not arise in
    // practice (a single-end read's records are either all aligned at distinct
    // positions or a lone unmapped record), so this tie-break is deterministic and
    // never actually consulted. Emit that record's seq/qual under its UNTRIMMED id
    // (matching the coordinate single-end emit's `GetReadId()` / `raw_name`).
    if single_end {
        let selected = group
            .iter()
            .filter(|record| per_record_alignment_candidate(record, filter, genes, scratch))
            .min_by_key(|record| single_end_coordinate_key(record));
        if let Some(record) = selected {
            let read = ReadRecord {
                id: record.raw_name.clone(),
                seq: record.seq.clone(),
                qual: Some(record.qual.clone()),
            };
            sink.emit_pair(&read, None)?;
            *emitted += 1;
        }
        return Ok(());
    }

    // Emit PRIMARIES ONLY (same id/pairing rules as flush_grouped_candidate).
    match primaries.as_slice() {
        [] => Ok(()), // candidacy via a non-primary only, with no primary to emit
        [lone] => {
            // Reached only in paired mode (single_end == true returned above).
            // A lone paired primary is an ORPHAN (0x1 SET, mate absent from the
            // file): DROP it, exactly as the coordinate path does (pass 2 emits
            // only complete pairs). Emitting it lone would break set-equality
            // AND the fused `run` path's equal-mate-count invariant. A stray
            // 0x1-UNSET single-end record in an otherwise-paired file still
            // emits lone (its untrimmed id), matching the no-alignment flush.
            if !lone.is_paired {
                let read = ReadRecord {
                    id: lone.raw_name.clone(),
                    seq: lone.seq.clone(),
                    qual: Some(lone.qual.clone()),
                };
                sink.emit_pair(&read, None)?;
                *emitted += 1;
            }
            Ok(())
        }
        [a, b] => {
            let (mate1, mate2) = if a.is_first_mate && !b.is_first_mate {
                (a, b)
            } else if b.is_first_mate && !a.is_first_mate {
                (b, a)
            } else {
                (a, b)
            };
            let r1 = ReadRecord {
                id: mate1.trimmed_name.clone(),
                seq: mate1.seq.clone(),
                qual: Some(mate1.qual.clone()),
            };
            let r2 = ReadRecord {
                id: mate1.trimmed_name.clone(),
                seq: mate2.seq.clone(),
                qual: Some(mate2.qual.clone()),
            };
            sink.emit_pair(&r1, Some(&r2))?;
            *emitted += 1;
            Ok(())
        }
        // Guarded by the >2-primary bail above; defensive only.
        _ => bail!(
            "internal error: grouped alignment flush reached a {}-primary group after the \
             >2-primary abort check",
            primaries.len()
        ),
    }
}

/// Runs the grouped/name-sorted (`GO:query`/`SO:queryname`) `--bam-mode
/// alignment` ONE-PASS BAM/CRAM read-extraction driver: the ALIGNMENT analogue
/// of [`extract_from_bam_no_alignment_grouped`]. A grouped/name-sorted BAM
/// keeps a template's records (both mates, plus any secondary/supplementary
/// alignments) adjacent in file order, so a SINGLE streaming pass can both
/// classify every record and reunite its QNAME group -- reproducing the
/// coordinate two-pass path's ([`extract_from_bam_with_threads`]) candidate SET
/// without a `rewind()` (so `alignments` may be [`Alignments::from_stdin`] when
/// `seekable == false`).
///
/// # Set-equality with the coordinate path (the correctness contract)
///
/// This is a unum EXTENSION beyond T1K (T1K has no grouped-alignment path), so
/// its correctness is defined by SET-EQUALITY of the emitted candidates with
/// the coordinate path -- proven by Task 9's differential test. Every
/// per-record decision here reproduces [`scan_pass1_chunk`]'s
/// `Selection::Alignment` status-conditioned branch (see
/// [`per_record_alignment_candidate`] and [`flush_alignment_grouped_candidate`]
/// for the branch-by-branch correspondence), and a template is a candidate iff
/// ANY of its records triggers candidacy -- the union that mirrors the
/// coordinate path's idempotent `candidates` insertion per template. Only
/// PRIMARY sequences are emitted (the coordinate path's pass 2 fetches primary
/// mates by name), so a candidacy triggered by a secondary alignment still
/// emits the primary pair.
///
/// # Setup: `seekable` (file) vs `!seekable` (stdin)
///
/// - `seekable`: `general_info(true)` + `rewind()` give the bitwise-identical
///   `single_end`/`read_len` the coordinate path uses, and
///   [`compute_hit_len_required`] (21/17, `read_len/5`) the identical
///   `hitLenRequired`.
/// - `!seekable`: buffers a bounded HEAD of raw records (up to
///   [`HIT_LEN_SAMPLE_SIZE`], all records -- the group loop consumes them),
///   deriving `single_end` from the primary `0x1`-FLAG majority
///   (`has_mate_cnt < head_primary_cnt / 2`, the negation of
///   [`Alignments::general_info`]'s `mate_paired` rule) and `read_len` as the
///   max primary `l_qseq`, then `compute_hit_len_required(if single_end {0}
///   else {1}, read_len)`.
///
/// Both then apply `set_hit_len_required` + the `infer_kmer_length` /
/// conditional `update_kmer_length` block copied from
/// [`extract_from_bam_with_threads`]. Unlike the no-alignment grouped path,
/// this does NOT call `set_ref_seq_similarity` (the alignment mode never pins
/// similarity -- matching the coordinate alignment path).
///
/// # Sink FACTORY
///
/// Like [`extract_from_bam_no_alignment_grouped`], takes `make_sink`, a
/// `FnOnce(single_end) -> Result<S>` factory, and calls it the moment
/// `single_end` is known (after setup, before the group loop emits). Returns
/// `(metrics, single_end, sink)`.
///
/// `threads` is accepted for signature parity with the sibling entry points
/// but unused: a QNAME group is tiny, so there is no useful unit of parallel
/// work.
///
/// # Errors
///
/// Returns an error if: `make_sink` fails; a QNAME group has more than 2
/// PRIMARY records; the setup `general_info`/`rewind` fails (seekable); or an
/// underlying BAM read fails (a genuine parse error, distinct from clean EOF).
#[allow(clippy::too_many_arguments)]
pub fn extract_from_bam_alignment_grouped<S, F>(
    alignments: &mut Alignments,
    filter: &mut RefKmerFilter,
    genes: &[GeneInterval],
    abnormal_unaligned_flag: bool,
    mate_id_len: i32,
    threads: usize,
    seekable: bool,
    make_sink: F,
) -> Result<(BamExtractMetrics, bool, S)>
where
    S: CandidateSink,
    F: FnOnce(bool) -> Result<S>,
{
    let _ = threads;
    extract_from_bam_alignment_grouped_with_head_limit(
        alignments,
        filter,
        genes,
        abnormal_unaligned_flag,
        mate_id_len,
        seekable,
        HIT_LEN_SAMPLE_SIZE,
        make_sink,
    )
}

/// [`extract_from_bam_alignment_grouped`]'s implementation, with the
/// sampled-head bound taken as an explicit parameter (only consulted on the
/// `!seekable` path) so `#[cfg(test)]` could force a tiny head to exercise a
/// QNAME group straddling the head/live boundary, mirroring
/// [`extract_from_bam_no_alignment_grouped_with_head_limit`].
#[allow(clippy::too_many_arguments)]
fn extract_from_bam_alignment_grouped_with_head_limit<S, F>(
    alignments: &mut Alignments,
    filter: &mut RefKmerFilter,
    genes: &[GeneInterval],
    abnormal_unaligned_flag: bool,
    mate_id_len: i32,
    seekable: bool,
    head_limit: usize,
    make_sink: F,
) -> Result<(BamExtractMetrics, bool, S)>
where
    S: CandidateSink,
    F: FnOnce(bool) -> Result<S>,
{
    // -- (1) Setup: derive single_end + hitLenRequired. `head` is empty on the
    // seekable path (the group loop reads purely live from the rewound cursor)
    // and pre-filled on the non-seekable path (bounded head, no rewind).
    let mut head: VecDeque<AlignmentGroupedRecord> = VecDeque::new();
    let (single_end, mut hit_len_required) = if seekable {
        let general_info =
            alignments.general_info(true).context("computing general_info (grouped alignment)")?;
        alignments.rewind().context("rewinding after general_info (grouped alignment)")?;
        (
            general_info.frag_stdev == 0,
            compute_hit_len_required(general_info.frag_stdev, general_info.read_len),
        )
    } else {
        let mut head_primary_cnt: u64 = 0;
        let mut has_mate_cnt: u64 = 0;
        let mut max_read_len: i32 = 0;
        while head.len() < head_limit
            && alignments.next().context("reading BAM head (grouped alignment)")?
        {
            let record = read_alignment_grouped_record(alignments, mate_id_len);
            if record.is_primary {
                head_primary_cnt += 1;
                if record.is_paired {
                    has_mate_cnt += 1;
                }
                let len = i32::try_from(record.seq.len()).unwrap_or(i32::MAX);
                if len > max_read_len {
                    max_read_len = len;
                }
            }
            head.push_back(record);
        }
        // Negation of general_info's `has_mate_cnt >= total / 2` majority vote.
        let single_end = has_mate_cnt < head_primary_cnt / 2;
        // `compute_hit_len_required` keys only on `frag_stdev == 0`, so a
        // sentinel 0 (single-end) / 1 (paired) selects the 17/21 base exactly.
        let hlr = compute_hit_len_required(i32::from(!single_end), max_read_len);
        (single_end, hlr)
    };

    filter.set_hit_len_required(hit_len_required);
    // infer/update kmer length (copied from extract_from_bam_with_threads);
    // NO set_ref_seq_similarity (alignment mode does not pin similarity).
    let inferred = filter.infer_kmer_length();
    if inferred > filter.kmer_length() {
        filter.update_kmer_length(inferred);
        if inferred > usize::try_from(hit_len_required).unwrap_or(0) {
            hit_len_required = i32::try_from(inferred).unwrap_or(i32::MAX);
            filter.set_hit_len_required(hit_len_required);
        }
    }

    // `single_end` is now known -- create the sink via the factory before the
    // group loop emits (the same B2 chicken-and-egg fix the no-alignment
    // grouped path documents).
    let mut sink = make_sink(single_end)?;

    // -- (2) ONE continuous group-accumulation loop over head ++ live. Records
    // are grouped by trimmed-QNAME adjacency; a group is flushed on the first
    // differently-named record and once more at EOF. Unlike the no-alignment
    // loop, there is NO >2-member abort here (a group legitimately holds a
    // primary pair PLUS secondaries) -- the >2-PRIMARY guard lives in the flush.
    let mut scratch = Scratch::default();
    let mut current_group: Vec<AlignmentGroupedRecord> = Vec::new();
    let mut emitted: u64 = 0;

    while let Some(record) = next_alignment_buffered_or_live(&mut head, alignments, mate_id_len)? {
        if let Some(first) = current_group.first() {
            if first.trimmed_name != record.trimmed_name {
                flush_alignment_grouped_candidate(
                    &current_group,
                    filter,
                    genes,
                    single_end,
                    abnormal_unaligned_flag,
                    &mut scratch,
                    &mut sink,
                    &mut emitted,
                )?;
                current_group.clear();
            }
        }
        current_group.push(record);
    }
    flush_alignment_grouped_candidate(
        &current_group,
        filter,
        genes,
        single_end,
        abnormal_unaligned_flag,
        &mut scratch,
        &mut sink,
        &mut emitted,
    )?;

    let metrics = BamExtractMetrics {
        single_end,
        hit_len_required,
        kmer_length: filter.kmer_length(),
        pass1_emitted: emitted,
        candidates_recorded: 0,
        pass2_emitted: 0,
    };
    Ok((metrics, single_end, sink))
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
///
/// `selection` gates the alignment-position check below:
/// [`Selection::Alignment`] applies it unchanged (a template that is not
/// template-aligned, and `-u` was not given, is skipped -- it can never
/// complete a `candidates` entry recorded by the alignment-driven pass 1);
/// [`Selection::NoAlignment`] bypasses it entirely, since every recorded
/// candidate was selected purely by k-mer match on its own sequence
/// (pass 1's [`Pass1Site::KmerCandidate`]), independent of alignment
/// position -- an UNALIGNED mate of such a candidate must still be fetched
/// here.
fn run_pass2(
    alignments: &mut Alignments,
    mut candidates: HashMap<String, PendingCandidate>,
    abnormal_unaligned_flag: bool,
    mate_id_len: i32,
    selection: Selection,
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
        if selection == Selection::Alignment
            && !alignments.is_template_aligned()
            && !abnormal_unaligned_flag
        {
            continue; // alignment gate -- bypassed under Selection::NoAlignment
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
    fn read_overlaps_any_gene_matches_cursor_including_nested() {
        let genes = vec![
            GeneInterval { chr_id: 0, start: 100, end: 500 }, // contains the next
            GeneInterval { chr_id: 0, start: 200, end: 250 },
            GeneInterval { chr_id: 1, start: 50, end: 150 },
        ];
        // Nested case that a binary-search-on-end would wrongly drop:
        assert!(read_overlaps_any_gene(0, 260, 270, &genes)); // overlaps [100,500]
        assert!(read_overlaps_any_gene(0, 120, 180, &genes));
        assert!(read_overlaps_any_gene(0, 500, 550, &genes)); // start==end boundary: 500<=500 && 550>100
        assert!(!read_overlaps_any_gene(0, 40, 100, &genes)); // end==start boundary rejected: 100 > 100 false
        assert!(!read_overlaps_any_gene(0, 600, 700, &genes));
        assert!(!read_overlaps_any_gene(2, 120, 180, &genes)); // wrong chr
        assert!(read_overlaps_any_gene(1, 100, 120, &genes));
        assert!(!read_overlaps_any_gene(0, 120, 180, &[]));
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

    #[derive(Debug, Default)]
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
            run_pass2(&mut alignments, candidates, false, -1, Selection::Alignment, &mut sink)
                .unwrap();
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

    // -- Pass-2 alignment-gate bypass + `extract_from_bam_no_alignment`
    // entry point: proves both that `run_pass2`'s gate is only active under
    // `Selection::Alignment`, and that the coordinate no-alignment 2-pass
    // wires FASTQ-pinned setup + k-mer pass 1 + gate-bypassed pass 2
    // end-to-end.

    /// A 1-pair coordinate BAM for the pass-2 alignment-gate-bypass
    /// regression: mate1 is UNALIGNED (`0x4`, `tid == -1`) but its sequence
    /// k-mer-matches (the same `ref_bytes[0..70]` substring
    /// `build_kmer_selection_test_bam`'s `unaligned_kmer` uses, already
    /// proven to pass the default-similarity k-mer filter); mate2 is
    /// ALIGNED, off-target (position is irrelevant under no-alignment, which
    /// never consults gene intervals). Per `Alignments::is_template_aligned`'s
    /// doc comment, `tid < 0` alone fails template-alignment for mate1, so
    /// pass 2's `!is_template_aligned()` gate would drop it under
    /// `Selection::Alignment` -- only `Selection::NoAlignment`'s bypass lets
    /// both mates be fetched and the pair complete.
    fn build_no_alignment_gate_bypass_test_bam(path: &std::path::Path) {
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

        // mate1: UNALIGNED (tid = -1), k-mer-matching SEQ.
        let seq1 = &ref_bytes[0..70];
        let mut r1 = bam::Record::new();
        r1.set(b"gate_pair", None, seq1, &[25u8; 70]);
        r1.set_tid(-1);
        r1.set_pos(-1);
        r1.set_mtid(-1);
        r1.set_mpos(-1);
        r1.set_flags(0x1 | 0x4 | 0x40);
        writer.write(&r1).unwrap();

        // mate2: ALIGNED, off-target, k-mer-matching SEQ.
        let seq2 = &ref_bytes[100..170];
        let mut r2 = bam::Record::new();
        r2.set(b"gate_pair", Some(&CigarString(vec![Cigar::Match(70)])), seq2, &[25u8; 70]);
        r2.set_tid(0);
        r2.set_pos(50_000);
        r2.set_mtid(-1);
        r2.set_mpos(-1);
        r2.set_flags(0x1 | 0x8 | 0x80);
        writer.write(&r2).unwrap();

        drop(writer);
    }

    #[test]
    fn no_alignment_pass2_fetches_kmer_passing_unaligned_reads() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("no_alignment_gate_bypass.bam");
        build_no_alignment_gate_bypass_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        // Full no-alignment pipeline: FASTQ-pinned setup + k-mer pass 1 +
        // gate-bypassed pass 2.
        let mut filter =
            RefKmerFilter::from_reference_fasta(&ref_fasta, INITIAL_KMER_LENGTH).unwrap();
        let mut alignments = Alignments::open(&bam_path).unwrap();
        let mut sink = VecSink::default();
        let metrics =
            extract_from_bam_no_alignment(&mut alignments, &mut filter, 0.8, -1, 1, &mut sink)
                .unwrap();

        assert!(!metrics.single_end, "fixture is paired -- must take the pass-2 path");
        assert_eq!(
            metrics.pass1_emitted, 0,
            "no-alignment (paired input) never emits directly in pass 1"
        );
        assert_eq!(metrics.candidates_recorded, 1, "the one QNAME must be recorded as a candidate");
        assert_eq!(
            metrics.pass2_emitted, 1,
            "the pair must be completed in pass 2 despite mate1 being UNALIGNED -- proves the \
             alignment gate was bypassed"
        );

        assert_eq!(sink.pairs.len(), 1);
        let (r1, r2) = &sink.pairs[0];
        assert_eq!(r1.id, "gate_pair");
        let r2 = r2.as_ref().expect("mate2 (the ALIGNED, off-target mate) must be present");
        assert_eq!(r2.id, "gate_pair");

        // Differential: replaying the SAME candidate through `run_pass2`
        // under `Selection::Alignment` must NOT complete the pair -- mate1's
        // `!is_template_aligned()` (tid < 0) trips the (still-active)
        // alignment gate, so `entry.mate1` never fills and the candidate is
        // left incomplete. This is the concrete "the alignment gate would
        // drop it" check: without the `selection`-guarded bypass, this exact
        // fixture would silently lose mate1 forever.
        let mut alignment_candidates = HashMap::new();
        alignment_candidates.insert("gate_pair".to_string(), PendingCandidate::default());
        let mut reopened = Alignments::open(&bam_path).unwrap();
        let mut alignment_sink = VecSink::default();
        let alignment_emitted = run_pass2(
            &mut reopened,
            alignment_candidates,
            false,
            -1,
            Selection::Alignment,
            &mut alignment_sink,
        )
        .unwrap();
        assert_eq!(
            alignment_emitted, 0,
            "under Selection::Alignment the gate must drop the UNALIGNED mate1, leaving the \
             candidate incomplete"
        );
        assert!(alignment_sink.pairs.is_empty());
    }

    // -- Grouped/name-sorted no-alignment one-pass
    // (`extract_from_bam_no_alignment_grouped`): stdin-capable single-pass
    // extraction that k-mer-tests each primary read and reunites its QNAME
    // group by adjacency instead of a coordinate two-pass name map.

    /// Builds a name-sorted (`SO:queryname`) BAM exercising every grouped
    /// no-alignment pairing rule:
    /// - `pair_ok`: 2 adjacent primary records (mate1's SEQ k-mer-matches;
    ///   mate2's is `noise` and fails on its own) -- must be emitted as a
    ///   PAIR via OR-rescue (mate1 alone passing is enough).
    /// - `single_ok`: 1 primary record, `0x1` UNSET, k-mer-matching SEQ --
    ///   must be emitted LONE.
    /// - `orphan_ok`: 1 primary record, `0x1` SET (this template's other
    ///   mate never shows up in the file), k-mer-matching SEQ -- must be
    ///   DROPPED, matching the coordinate path (which emits only complete
    ///   pairs) and keeping the fused `run` path's mate counts equal.
    /// - `pair_fail`: 2 adjacent primary records, BOTH `noise` (fail
    ///   individually, so OR-rescue has nothing to rescue with) -- must NOT
    ///   be emitted at all.
    fn build_grouped_no_alignment_test_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "queryname");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);

        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        // pair_ok: mate1 passes on its own; mate2 is noise and fails on its
        // own -- OR-rescue must still emit the pair.
        let mut m1 = bam::Record::new();
        m1.set(b"pair_ok", None, &ref_bytes[0..90], &[30u8; 90]);
        m1.set_tid(-1);
        m1.set_pos(-1);
        m1.set_mtid(-1);
        m1.set_mpos(-1);
        m1.set_flags(0x1 | 0x40);
        writer.write(&m1).unwrap();

        let pair_ok_mate2_seq = noise(90);
        let mut m2 = bam::Record::new();
        m2.set(b"pair_ok", None, &pair_ok_mate2_seq, &[30u8; 90]);
        m2.set_tid(-1);
        m2.set_pos(-1);
        m2.set_mtid(-1);
        m2.set_mpos(-1);
        m2.set_flags(0x1 | 0x80);
        writer.write(&m2).unwrap();

        // single_ok: 0x1 UNSET, k-mer-matching SEQ -- lone emit.
        let mut single = bam::Record::new();
        single.set(b"single_ok", None, &ref_bytes[100..190], &[30u8; 90]);
        single.set_tid(-1);
        single.set_pos(-1);
        single.set_mtid(-1);
        single.set_mpos(-1);
        single.set_flags(0);
        writer.write(&single).unwrap();

        // orphan_ok: 0x1 SET, mate never shows up in the file,
        // k-mer-matching SEQ -- lone emit (orphan OR-rescue).
        let mut orphan = bam::Record::new();
        orphan.set(b"orphan_ok", None, &ref_bytes[200..290], &[30u8; 90]);
        orphan.set_tid(-1);
        orphan.set_pos(-1);
        orphan.set_mtid(-1);
        orphan.set_mpos(-1);
        orphan.set_flags(0x1 | 0x40);
        writer.write(&orphan).unwrap();

        // pair_fail: both mates are noise -- neither passes, so OR-rescue
        // has nothing to rescue with; must not be emitted.
        let pair_fail_mate1_seq = noise(90);
        let mut f1 = bam::Record::new();
        f1.set(b"pair_fail", None, &pair_fail_mate1_seq, &[30u8; 90]);
        f1.set_tid(-1);
        f1.set_pos(-1);
        f1.set_mtid(-1);
        f1.set_mpos(-1);
        f1.set_flags(0x1 | 0x40);
        writer.write(&f1).unwrap();

        let pair_fail_mate2_seq = noise(90);
        let mut f2 = bam::Record::new();
        f2.set(b"pair_fail", None, &pair_fail_mate2_seq, &[30u8; 90]);
        f2.set_tid(-1);
        f2.set_pos(-1);
        f2.set_mtid(-1);
        f2.set_mpos(-1);
        f2.set_flags(0x1 | 0x80);
        writer.write(&f2).unwrap();

        drop(writer);
    }

    #[test]
    fn grouped_no_alignment_pairs_singles_and_drops_orphans() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("grouped_no_alignment.bam");
        build_grouped_no_alignment_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut filter =
            RefKmerFilter::from_reference_fasta(&ref_fasta, INITIAL_KMER_LENGTH).unwrap();
        let mut alignments = Alignments::open(&bam_path).unwrap();

        let (metrics, single_end, sink) = extract_from_bam_no_alignment_grouped(
            &mut alignments,
            &mut filter,
            0.8,
            -1,
            1,
            |_single_end| Ok(VecSink::default()),
        )
        .unwrap();

        // 6 primary records total; 5 of them (pair_ok x2, orphan_ok,
        // pair_fail x2) carry the 0x1 (paired) FLAG bit -- the majority-vote
        // rule must classify this fixture as paired, not single-end.
        assert!(!single_end, "majority-paired fixture must not be classified single-end");
        assert_eq!(metrics.single_end, single_end, "returned tuple and metrics field must agree");
        assert_eq!(
            metrics.pass1_emitted, 2,
            "pair_ok + single_ok must be emitted; orphan_ok (0x1 SET, mate absent) must be DROPPED \
             to match the coordinate path (which emits only complete pairs); pair_fail must not"
        );
        assert_eq!(
            metrics.candidates_recorded, 0,
            "the one-pass path never populates a candidates map"
        );
        assert_eq!(metrics.pass2_emitted, 0, "the one-pass path has no pass 2");

        let mut by_id: HashMap<String, (ReadRecord, Option<ReadRecord>)> = HashMap::new();
        for (r1, r2) in &sink.pairs {
            by_id.insert(r1.id.clone(), (r1.clone(), r2.clone()));
        }
        assert_eq!(
            by_id.len(),
            2,
            "exactly 2 QNAMEs must be emitted, got {:?}",
            by_id.keys().collect::<Vec<_>>()
        );

        let (pair_r1, pair_r2) = by_id.get("pair_ok").expect("pair_ok must be emitted as a pair");
        let pair_r2 =
            pair_r2.as_ref().expect("pair_ok must be emitted WITH its mate2 (a pair, not lone)");
        assert_eq!(pair_r1.seq, ref_bytes[0..90], "mate1 (is_first_mate) must come first");
        assert_eq!(pair_r2.seq, noise(90), "mate2 (noise, the failing member) must come second");
        assert_eq!(pair_r2.id, "pair_ok", "mate2's emitted id must be mate1's id, not its own");

        let (single_r1, single_r2) = by_id.get("single_ok").expect("single_ok must be emitted");
        assert!(single_r2.is_none(), "single_ok (0x1 unset) must be emitted LONE, not as a pair");
        assert_eq!(single_r1.seq, ref_bytes[100..190]);

        assert!(
            !by_id.contains_key("orphan_ok"),
            "an orphan (0x1 SET, mate absent) must be DROPPED to match the coordinate path (a paired \
             candidate whose mate never arrives is not emitted); emitting it lone breaks set-equality \
             AND the fused `run` path's equal-mate-count invariant"
        );
        assert!(
            !by_id.contains_key("pair_fail"),
            "pair_fail (both mates noise) must not be emitted"
        );
    }

    /// Builds a name-sorted BAM with 3 primary records sharing one QNAME --
    /// a malformed "grouped/name-sorted" claim the one-pass entry point must
    /// reject rather than silently mis-grouping reads.
    fn build_grouped_no_alignment_three_primaries_test_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "queryname");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);

        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");
        for i in 0..3usize {
            let seq = &ref_bytes[i * 90..i * 90 + 90];
            let mut r = bam::Record::new();
            r.set(b"triple", None, seq, &[30u8; 90]);
            r.set_tid(-1);
            r.set_pos(-1);
            r.set_mtid(-1);
            r.set_mpos(-1);
            r.set_flags(0x1);
            writer.write(&r).unwrap();
        }
        drop(writer);
    }

    #[test]
    fn grouped_no_alignment_more_than_two_primaries_aborts() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("triple_qname.bam");
        build_grouped_no_alignment_three_primaries_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let mut filter =
            RefKmerFilter::from_reference_fasta(&ref_fasta, INITIAL_KMER_LENGTH).unwrap();
        let mut alignments = Alignments::open(&bam_path).unwrap();

        let result = extract_from_bam_no_alignment_grouped(
            &mut alignments,
            &mut filter,
            0.8,
            -1,
            1,
            |_single_end| Ok(VecSink::default()),
        );
        let err = result.expect_err("a >2-member QNAME group must abort with an error");
        let message = format!("{err:#}");
        assert!(
            message.contains("triple") && message.contains("two"),
            "error should hint at the offending QNAME and the >2-member rule, got: {message}"
        );
    }

    #[test]
    fn grouped_no_alignment_reunites_group_straddling_head_boundary() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("boundary_straddle.bam");
        // Reuses the same fixture as `grouped_no_alignment_pairs_orphans_and_singles`:
        // `pair_ok`'s two records are the FIRST two primary records in the
        // file.
        build_grouped_no_alignment_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let mut filter =
            RefKmerFilter::from_reference_fasta(&ref_fasta, INITIAL_KMER_LENGTH).unwrap();
        let mut alignments = Alignments::open(&bam_path).unwrap();

        // `head_limit = 1`: only `pair_ok`'s FIRST record (mate1) is
        // buffered into the head before the setup phase stops sampling;
        // mate2 is therefore read LIVE via `next_buffered_or_live`, in a
        // DIFFERENT sub-phase from mate1 despite sharing `pair_ok`'s QNAME.
        // If the group loop treated the head as a separately-paired buffer
        // (instead of draining head ++ live as one continuous stream), this
        // group would be split into two 1-member flushes instead of being
        // reunited into one 2-member group.
        let (metrics, _single_end, sink) = extract_from_bam_no_alignment_grouped_with_head_limit(
            &mut alignments,
            &mut filter,
            0.8,
            -1,
            1,
            |_single_end| Ok(VecSink::default()),
        )
        .unwrap();

        assert_eq!(
            metrics.pass1_emitted, 2,
            "boundary-straddled pair_ok must still merge into one 2-member group and emit \
             correctly, exactly like the unbounded-head run \
             (grouped_no_alignment_pairs_singles_and_drops_orphans: pair_ok + single_ok emitted, \
             orphan_ok dropped)"
        );
        let pair_entry = sink
            .pairs
            .iter()
            .find(|(r1, _)| r1.id == "pair_ok")
            .expect("pair_ok must be emitted despite straddling the head/live boundary");
        assert!(
            pair_entry.1.is_some(),
            "pair_ok must be emitted as a PAIR (both members reunited), not split into two \
             independent lone emissions"
        );
    }

    /// Builds a single-end (all `0x1` UNSET) name-sorted BAM. The read under
    /// test, `se_read/1`, carries a trailing `/1` suffix and a k-mer-matching
    /// SEQ; under the default `mate_id_len = -1`, `trim_name("se_read/1", -1)`
    /// strips the `/1` to `"se_read"` -- so an emit that (wrongly) used the
    /// TRIMMED name would drop the suffix, whereas the coordinate
    /// no-alignment path's single-end emit keeps the raw `read_id()` verbatim.
    /// Two additional plain single-end reads (`se_pad_a`/`se_pad_b`, distinct
    /// QNAMEs) are included so the `0x1`-FLAG majority vote clearly classifies
    /// the file single-end (`has_mate_cnt = 0 >= total/2` fails for
    /// `total = 3` -- the identical majority rule `general_info` uses; a
    /// single lone read would land in the `0 >= 0` degenerate tie that
    /// `general_info` itself treats as paired).
    fn build_grouped_no_alignment_single_end_suffix_test_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "queryname");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);

        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        // Name-sorted order: "se_pad_a", "se_pad_b", "se_read/1".
        let mut pad_a = bam::Record::new();
        pad_a.set(b"se_pad_a", None, &ref_bytes[100..190], &[30u8; 90]);
        pad_a.set_tid(-1);
        pad_a.set_pos(-1);
        pad_a.set_mtid(-1);
        pad_a.set_mpos(-1);
        pad_a.set_flags(0);
        writer.write(&pad_a).unwrap();

        let mut pad_b = bam::Record::new();
        pad_b.set(b"se_pad_b", None, &ref_bytes[200..290], &[30u8; 90]);
        pad_b.set_tid(-1);
        pad_b.set_pos(-1);
        pad_b.set_mtid(-1);
        pad_b.set_mpos(-1);
        pad_b.set_flags(0);
        writer.write(&pad_b).unwrap();

        let mut r = bam::Record::new();
        r.set(b"se_read/1", None, &ref_bytes[0..90], &[30u8; 90]);
        r.set_tid(-1);
        r.set_pos(-1);
        r.set_mtid(-1);
        r.set_mpos(-1);
        r.set_flags(0); // 0x1 UNSET -> genuine single-end.
        writer.write(&r).unwrap();
        drop(writer);
    }

    #[test]
    fn grouped_no_alignment_single_end_lone_emit_keeps_untrimmed_id() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("single_end_suffix.bam");
        build_grouped_no_alignment_single_end_suffix_test_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let mut filter =
            RefKmerFilter::from_reference_fasta(&ref_fasta, INITIAL_KMER_LENGTH).unwrap();
        let mut alignments = Alignments::open(&bam_path).unwrap();

        // Default `mate_id_len = -1` (strips a trailing `/1`/`/2` for the
        // grouping compare, but the single-end LONE emit must NOT trim its id).
        let (metrics, single_end, sink) = extract_from_bam_no_alignment_grouped(
            &mut alignments,
            &mut filter,
            0.8,
            -1,
            1,
            |_single_end| Ok(VecSink::default()),
        )
        .unwrap();

        assert!(single_end, "an all-0x1-unset fixture must be classified single-end");
        assert_eq!(metrics.pass1_emitted, 3, "all three single-end reads must be emitted lone");
        // Every single-end read is emitted lone (r2 == None).
        assert!(
            sink.pairs.iter().all(|(_, r2)| r2.is_none()),
            "single-end reads must be emitted LONE"
        );
        let suffix_read = sink
            .pairs
            .iter()
            .find(|(r1, _)| r1.seq == ref_bytes_first_90())
            .expect("the se_read/1 read (seq == ref[0..90]) must be emitted")
            .0
            .clone();
        // The load-bearing assertion: the emitted id is the RAW, untrimmed
        // name. Under the old (trimmed-name) behavior this would be
        // "se_read" and the assertion would FAIL -- matching the coordinate
        // no-alignment path, whose single-end emit keeps `read_id()` verbatim.
        assert_eq!(
            suffix_read.id, "se_read/1",
            "single-end lone emit must carry the UNTRIMMED read id (the coordinate no-alignment \
             path keeps read_id() verbatim for single-end); got the trimmed form instead"
        );
    }

    /// The first 90 bytes of [`PARALLEL_TEST_REF`] -- the `se_read/1`
    /// fixture read's SEQ, used to locate it among the emitted single-end
    /// reads without relying on emit order.
    fn ref_bytes_first_90() -> Vec<u8> {
        PARALLEL_TEST_REF.as_bytes()[0..90].to_vec()
    }

    // ---------------------------------------------------------------------
    // Task 7: one-pass grouped/name-sorted ALIGNMENT extractor
    // (`extract_from_bam_alignment_grouped`). Each test below exercises one
    // status-conditioned classification branch, asserting the emitted
    // candidate SET (via a `HashSet` derived from the collecting `VecSink`)
    // matches what the coordinate path's `scan_pass1_chunk` would select for
    // the same records -- the set-equality invariant Task 9's differential
    // test proves end-to-end. Fixtures are `@HD SO:queryname` (grouped/
    // name-sorted), two contigs: `chr1` (gene interval `[1000, 5000]`) and
    // `chr1_alt` (alt contig, name contains `_`).

    /// Builds the shared `@HD SO:queryname` header with the `chr1`/`chr1_alt`
    /// contigs every grouped-alignment fixture uses.
    fn grouped_alignment_header() -> Header {
        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "queryname");
        header.push_record(&hd);
        let mut sq1 = HeaderRecord::new(b"SQ");
        sq1.push_tag(b"SN", "chr1");
        sq1.push_tag(b"LN", 1_000_000);
        header.push_record(&sq1);
        let mut sq2 = HeaderRecord::new(b"SQ");
        sq2.push_tag(b"SN", "chr1_alt");
        sq2.push_tag(b"LN", 1_000_000);
        header.push_record(&sq2);
        header
    }

    /// Writes one record to a grouped-alignment fixture. `cigar = None` (with
    /// `tid = pos = -1`) writes an unmapped record; otherwise an aligned one.
    #[allow(clippy::too_many_arguments)]
    fn add_grouped_alignment_record(
        writer: &mut Writer,
        name: &[u8],
        seq: &[u8],
        flags: u16,
        tid: i32,
        pos: i64,
        mtid: i32,
        mpos: i64,
        cigar: Option<&CigarString>,
    ) {
        let qual = vec![30u8; seq.len()];
        let mut r = bam::Record::new();
        r.set(name, cigar, seq, &qual);
        r.set_tid(tid);
        r.set_pos(pos);
        r.set_mtid(mtid);
        r.set_mpos(mpos);
        r.set_flags(flags);
        writer.write(&r).unwrap();
    }

    /// The `chr1 [1000, 5000]` gene interval every grouped-alignment fixture
    /// selects on-target reads against.
    fn grouped_alignment_genes(alignments: &Alignments) -> Vec<GeneInterval> {
        build_genes(
            alignments,
            &[CoordRecord {
                name: "only".to_string(),
                chrom: "chr1".to_string(),
                start: 1000,
                end: 5000,
                strand: "+".to_string(),
                seq: PARALLEL_TEST_REF.to_string(),
            }],
        )
        .unwrap()
    }

    /// Runs [`extract_from_bam_alignment_grouped`] end to end over `bam_path`,
    /// returning the metrics, derived `single_end`, and the collecting sink.
    fn run_grouped_alignment(
        bam_path: &std::path::Path,
        ref_fasta: &std::path::Path,
        abnormal_unaligned_flag: bool,
        seekable: bool,
    ) -> (BamExtractMetrics, bool, VecSink) {
        let mut filter =
            RefKmerFilter::from_reference_fasta(ref_fasta, INITIAL_KMER_LENGTH).unwrap();
        let mut alignments = Alignments::open(bam_path).unwrap();
        let genes = grouped_alignment_genes(&alignments);
        extract_from_bam_alignment_grouped(
            &mut alignments,
            &mut filter,
            &genes,
            abnormal_unaligned_flag,
            -1,
            1,
            seekable,
            |_single_end| Ok(VecSink::default()),
        )
        .unwrap()
    }

    /// The SET of read-1 QNAMEs emitted into a [`VecSink`].
    fn emitted_ids(sink: &VecSink) -> std::collections::HashSet<String> {
        sink.pairs.iter().map(|(r1, _)| r1.id.clone()).collect()
    }

    /// On-target main-chrom aligned pair inside `[1000, 5000]` with
    /// ref-substring SEQs: `is_template_aligned && is_aligned && !chrom_alt`,
    /// overlaps a gene, not low-complexity -> CANDIDATE (no seed test needed).
    /// A second, off-target pair with the SAME ref-substring SEQs (so it WOULD
    /// seed-hit) confirms only the on-target one is selected.
    fn build_grouped_on_target_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();
        let header = grouped_alignment_header();
        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");
        let cig90 = CigarString(vec![Cigar::Match(90)]);

        // Name-sorted: "off_seed" < "on_target".
        // off_seed: main chrom, OUTSIDE the gene, ref-substring SEQs (seed-hit)
        // -- must be EXCLUDED (the main-chrom aligned branch has no seed test).
        // No reverse (0x10) flags: read_seq() then stores/emits SEQ forward, so
        // the emitted-SEQ assertions compare directly against ref substrings.
        add_grouped_alignment_record(
            &mut writer,
            b"off_seed",
            &ref_bytes[0..90],
            0x1 | 0x2 | 0x40,
            0,
            50_000,
            0,
            50_200,
            Some(&cig90),
        );
        add_grouped_alignment_record(
            &mut writer,
            b"off_seed",
            &ref_bytes[100..190],
            0x1 | 0x2 | 0x80,
            0,
            50_200,
            0,
            50_000,
            Some(&cig90),
        );
        // on_target: main chrom, INSIDE the gene -> CANDIDATE.
        add_grouped_alignment_record(
            &mut writer,
            b"on_target",
            &ref_bytes[0..90],
            0x1 | 0x2 | 0x40,
            0,
            1100,
            0,
            1300,
            Some(&cig90),
        );
        add_grouped_alignment_record(
            &mut writer,
            b"on_target",
            &ref_bytes[100..190],
            0x1 | 0x2 | 0x80,
            0,
            1300,
            0,
            1100,
            Some(&cig90),
        );
        drop(writer);
    }

    #[test]
    fn grouped_alignment_on_target_main_chrom_pair_is_candidate() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("on_target.bam");
        build_grouped_on_target_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let (metrics, single_end, sink) = run_grouped_alignment(&bam_path, &ref_fasta, false, true);

        assert!(!single_end, "paired fixture must not be single-end");
        assert_eq!(
            emitted_ids(&sink),
            ["on_target".to_string()].into_iter().collect(),
            "only the on-target pair is a candidate; the off-target pair (same seed-hitting \
             SEQs) must be excluded"
        );
        assert_eq!(metrics.pass1_emitted, 1);
        // The on-target pair is emitted as a PAIR of PRIMARY records.
        let (r1, r2) = sink.pairs.iter().find(|(r, _)| r.id == "on_target").unwrap();
        let r2 = r2.as_ref().expect("on-target must emit both mates");
        assert_eq!(r1.seq, ref_bytes_first_90(), "mate1 (is_first_mate) first");
        assert_eq!(r2.seq, PARALLEL_TEST_REF.as_bytes()[100..190]);
    }

    #[test]
    fn grouped_alignment_off_target_main_chrom_seed_hit_is_not_candidate() {
        // C-1: the main-chrom aligned branch selects purely by gene overlap;
        // there is NO seed/k-mer test on it. `build_grouped_on_target_bam`
        // includes an off-target pair whose SEQs DO seed-hit -- if the branch
        // (wrongly) applied a seed test, that pair would leak in.
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("off_seed.bam");
        build_grouped_on_target_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let (_metrics, _single_end, sink) =
            run_grouped_alignment(&bam_path, &ref_fasta, false, true);

        assert!(
            !emitted_ids(&sink).contains("off_seed"),
            "off-target main-chrom pair must NOT be a candidate even though its SEQs seed-hit"
        );
    }

    /// Genuine unaligned pair (both mates unmapped, `0x1` set): the JOINT
    /// predicate `!lc(s1) && !lc(s2) && (good(s1) || good(s2))`.
    /// - `unaligned_good`: both mates good, both non-lc -> JOINT passes.
    /// - `unaligned_onelc`: mate1 homopolymer (low-complexity), mate2 good --
    ///   JOINT's `!lc(s1)` fails, so the template is EXCLUDED. A per-read OR
    ///   would (wrongly) rescue it via mate2; asserting exclusion proves the
    ///   joint AND-gate.
    fn build_grouped_unaligned_pair_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();
        let header = grouped_alignment_header();
        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        // unaligned_good: both good.
        add_grouped_alignment_record(
            &mut writer,
            b"unaligned_good",
            &ref_bytes[0..90],
            0x1 | 0x4 | 0x8 | 0x40,
            -1,
            -1,
            -1,
            -1,
            None,
        );
        add_grouped_alignment_record(
            &mut writer,
            b"unaligned_good",
            &ref_bytes[100..190],
            0x1 | 0x4 | 0x8 | 0x80,
            -1,
            -1,
            -1,
            -1,
            None,
        );
        // unaligned_onelc: mate1 low-complexity, mate2 good.
        add_grouped_alignment_record(
            &mut writer,
            b"unaligned_onelc",
            &[b'A'; 90],
            0x1 | 0x4 | 0x8 | 0x40,
            -1,
            -1,
            -1,
            -1,
            None,
        );
        add_grouped_alignment_record(
            &mut writer,
            b"unaligned_onelc",
            &ref_bytes[100..190],
            0x1 | 0x4 | 0x8 | 0x80,
            -1,
            -1,
            -1,
            -1,
            None,
        );
        drop(writer);
    }

    #[test]
    fn grouped_alignment_genuine_unaligned_pair_uses_joint_predicate() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("unaligned_pair.bam");
        build_grouped_unaligned_pair_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let (_metrics, single_end, sink) =
            run_grouped_alignment(&bam_path, &ref_fasta, false, true);

        assert!(!single_end, "paired (0x1-set) unaligned fixture must classify paired");
        assert_eq!(
            emitted_ids(&sink),
            ["unaligned_good".to_string()].into_iter().collect(),
            "the joint predicate requires BOTH mates non-low-complexity; unaligned_onelc \
             (mate1 homopolymer) must be excluded, not rescued per-read"
        );
        // unaligned_good emits BOTH primaries as a pair.
        let (_r1, r2) = sink.pairs.iter().find(|(r, _)| r.id == "unaligned_good").unwrap();
        assert!(r2.is_some(), "an unaligned pair emits both mates");
    }

    /// ALT-contig pair (mapped to `chr1_alt`, `chrom_alt` true): the PER-READ
    /// predicate `!lc && good`, NOT the joint one. `alt_mixed`'s mate1 is a
    /// homopolymer (fails `!lc`), mate2 is good -- per-read OR rescues via
    /// mate2, so the template IS a candidate. A joint AND-gate would exclude
    /// it (mate1 low-complexity), so emission proves the per-read branch.
    fn build_grouped_alt_contig_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();
        let header = grouped_alignment_header();
        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");
        let cig90 = CigarString(vec![Cigar::Match(90)]);

        // alt_mixed on chr1_alt (tid 1): mate1 low-complexity, mate2 good.
        add_grouped_alignment_record(
            &mut writer,
            b"alt_mixed",
            &[b'A'; 90],
            0x1 | 0x2 | 0x20 | 0x40,
            1,
            500,
            1,
            700,
            Some(&cig90),
        );
        add_grouped_alignment_record(
            &mut writer,
            b"alt_mixed",
            &ref_bytes[200..290],
            0x1 | 0x2 | 0x10 | 0x80,
            1,
            700,
            1,
            500,
            Some(&cig90),
        );
        drop(writer);
    }

    #[test]
    fn grouped_alignment_alt_contig_pair_uses_per_read_predicate() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("alt_contig.bam");
        build_grouped_alt_contig_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let (_metrics, single_end, sink) =
            run_grouped_alignment(&bam_path, &ref_fasta, false, true);

        assert!(!single_end);
        assert_eq!(
            emitted_ids(&sink),
            ["alt_mixed".to_string()].into_iter().collect(),
            "ALT-contig pairs use the per-read `!lc && good` predicate; mate2 (good) rescues \
             even though mate1 is low-complexity -- a joint AND-gate would wrongly exclude it"
        );
        let (_r1, r2) = sink.pairs.iter().find(|(r, _)| r.id == "alt_mixed").unwrap();
        assert!(r2.is_some(), "an ALT-contig pair emits both primary mates");
    }

    /// Abnormal-unmapped (`-u`, `abnormal_unaligned_flag = true`): an unaligned
    /// pair falls to the PER-READ predicate instead of the joint one. Same
    /// mate1-lc/mate2-good shape as the ALT case: per-read rescues via mate2.
    fn build_grouped_abnormal_unmapped_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();
        let header = grouped_alignment_header();
        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

        // abn_mixed: unmapped pair, mate1 low-complexity, mate2 good.
        add_grouped_alignment_record(
            &mut writer,
            b"abn_mixed",
            &[b'A'; 90],
            0x1 | 0x4 | 0x8 | 0x40,
            -1,
            -1,
            -1,
            -1,
            None,
        );
        add_grouped_alignment_record(
            &mut writer,
            b"abn_mixed",
            &ref_bytes[100..190],
            0x1 | 0x4 | 0x8 | 0x80,
            -1,
            -1,
            -1,
            -1,
            None,
        );
        drop(writer);
    }

    #[test]
    fn grouped_alignment_abnormal_unmapped_uses_per_read_predicate() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("abnormal.bam");
        build_grouped_abnormal_unmapped_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        // With -u the unaligned pair takes the per-read branch.
        let (_metrics, single_end, sink) = run_grouped_alignment(&bam_path, &ref_fasta, true, true);

        assert!(!single_end);
        assert_eq!(
            emitted_ids(&sink),
            ["abn_mixed".to_string()].into_iter().collect(),
            "under -u an unaligned pair uses the per-read `!lc && good` predicate (mate2 \
             rescues); the joint AND-gate must NOT apply"
        );
    }

    /// Secondary alignment overlapping a gene makes the template a candidate,
    /// but ONLY the PRIMARY records are emitted. `sec_gene`'s two primaries are
    /// OFF-target (main chrom, outside the gene) with `noise` SEQs (candidates
    /// on neither their position nor a seed test), while a SECONDARY of mate1
    /// lands INSIDE the gene with a ref-substring SEQ -- so the template is a
    /// candidate via the secondary, and the emitted pair carries the PRIMARY
    /// (off-target `noise`) SEQs, never the secondary's.
    fn build_grouped_secondary_gene_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();
        let header = grouped_alignment_header();
        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");
        let cig90 = CigarString(vec![Cigar::Match(90)]);

        let prim1 = noise(90);
        let prim2 = noise(80);
        // Primary mate1: off-target, noise (no reverse flag -> forward SEQ).
        add_grouped_alignment_record(
            &mut writer,
            b"sec_gene",
            &prim1,
            0x1 | 0x2 | 0x40,
            0,
            50_000,
            0,
            50_200,
            Some(&cig90),
        );
        // Primary mate2: off-target, noise.
        add_grouped_alignment_record(
            &mut writer,
            b"sec_gene",
            &prim2,
            0x1 | 0x2 | 0x80,
            0,
            50_200,
            0,
            50_000,
            Some(&CigarString(vec![Cigar::Match(80)])),
        );
        // Secondary of mate1 (0x100): INSIDE the gene, ref-substring SEQ.
        add_grouped_alignment_record(
            &mut writer,
            b"sec_gene",
            &ref_bytes[0..90],
            0x1 | 0x2 | 0x100 | 0x40,
            0,
            2000,
            0,
            50_200,
            Some(&cig90),
        );
        drop(writer);
    }

    #[test]
    fn grouped_alignment_secondary_overlapping_gene_makes_template_candidate_emitting_primary_only()
    {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("secondary_gene.bam");
        build_grouped_secondary_gene_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let (_metrics, single_end, sink) =
            run_grouped_alignment(&bam_path, &ref_fasta, false, true);

        assert!(!single_end);
        assert_eq!(
            emitted_ids(&sink),
            ["sec_gene".to_string()].into_iter().collect(),
            "a secondary alignment overlapping a gene makes the template a candidate"
        );
        let (r1, r2) = sink.pairs.iter().find(|(r, _)| r.id == "sec_gene").unwrap();
        let r2 = r2.as_ref().expect("must emit both primary mates");
        // The load-bearing assertion: the emitted SEQs are the PRIMARY
        // (off-target noise) sequences, NOT the on-target secondary's.
        assert_eq!(r1.seq, noise(90), "mate1 must be the PRIMARY (noise) seq, not the secondary");
        assert_eq!(r2.seq, noise(80), "mate2 must be the PRIMARY (noise) seq");
        assert_ne!(
            r1.seq,
            ref_bytes_first_90(),
            "the on-target SECONDARY sequence must never be emitted"
        );
    }

    /// Single-end (`0x1` UNSET) on-target read: emitted LONE, carrying the
    /// UNTRIMMED read id. `se_read/1`'s trailing `/1` is stripped for the
    /// grouping compare (default `mate_id_len = -1`) but must be KEPT in the
    /// emitted id -- matching the coordinate path's single-end emit. Two
    /// off-target single-end pads (`se_pad_a`/`se_pad_b`) make the `0x1`-unset
    /// majority classify the file single-end.
    fn build_grouped_single_end_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();
        let header = grouped_alignment_header();
        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");
        let cig90 = CigarString(vec![Cigar::Match(90)]);

        // Name-sorted: "se_pad_a" < "se_pad_b" < "se_read/1".
        // Pads: off-target -> not candidates (keep the emitted set clean).
        add_grouped_alignment_record(
            &mut writer,
            b"se_pad_a",
            &noise(90),
            0,
            0,
            60_000,
            -1,
            -1,
            Some(&cig90),
        );
        add_grouped_alignment_record(
            &mut writer,
            b"se_pad_b",
            &noise(90),
            0,
            0,
            61_000,
            -1,
            -1,
            Some(&cig90),
        );
        // se_read/1: on-target single-end -> candidate, emitted lone.
        add_grouped_alignment_record(
            &mut writer,
            b"se_read/1",
            &ref_bytes[0..90],
            0,
            0,
            1100,
            -1,
            -1,
            Some(&cig90),
        );
        drop(writer);
    }

    #[test]
    fn grouped_alignment_single_end_on_target_emits_untrimmed_id() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("single_end.bam");
        build_grouped_single_end_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let (metrics, single_end, sink) = run_grouped_alignment(&bam_path, &ref_fasta, false, true);

        assert!(single_end, "an all-0x1-unset fixture must classify single-end");
        assert_eq!(
            emitted_ids(&sink),
            ["se_read/1".to_string()].into_iter().collect(),
            "only the on-target single-end read is a candidate, and its id keeps the UNTRIMMED \
             /1 suffix"
        );
        assert_eq!(metrics.pass1_emitted, 1);
        let (r1, r2) = sink.pairs.iter().find(|(r, _)| r.id == "se_read/1").unwrap();
        assert!(r2.is_none(), "a single-end read is emitted LONE");
        assert_eq!(r1.seq, ref_bytes_first_90());
    }

    #[test]
    fn grouped_alignment_non_seekable_head_derivation_matches_seekable() {
        // The non-seekable (stdin) path never calls general_info/rewind; it
        // derives single_end + hitLenRequired from a bounded head instead.
        // Over a file-backed reader the two setups must select the same set.
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("on_target_nonseekable.bam");
        build_grouped_on_target_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let (_m_seek, se_seek, sink_seek) =
            run_grouped_alignment(&bam_path, &ref_fasta, false, true);
        let (_m_head, se_head, sink_head) =
            run_grouped_alignment(&bam_path, &ref_fasta, false, false);

        assert_eq!(se_seek, se_head, "single_end derivation must agree across setup paths");
        assert_eq!(
            emitted_ids(&sink_seek),
            emitted_ids(&sink_head),
            "the non-seekable head-derived setup must select the same candidate set"
        );
    }

    // -- Task 7 single-end EMISSION-selection faithfulness (Defect A) + the
    // on-target-ignores-good faithfulness (the reviewer's "Defect B", verified
    // here to be a NON-defect: the coordinate path -- `evaluate_pass1_site`'s
    // `SingleEndOnTarget` arm and T1K `BamExtractor.cpp:836-850` -- emits
    // on-target single-end reads after ONLY the low-complexity check, with NO
    // `good()`/`HasHitInSet` gate; see the green
    // `single_end_ontarget_non_candidate_read_is_emitted_unconditionally`
    // regression). These prove the grouped single-end emitter selects the SAME
    // record the coordinate path would (the first EMITTABLE in coordinate order,
    // which may be a NON-PRIMARY), keyed on (id, seq, qual) -- Task 9's gate,
    // previewed here for single-end.

    /// Writes single-end records `(name, seq, flags, tid, pos)` to `path` in the
    /// EXACT order given (the caller pre-sorts: coordinate order for the
    /// coordinate extractor, name/QNAME-grouped order for the grouped one),
    /// under the shared `chr1`/`chr1_alt` header. `cigar` is `Match(seq.len())`
    /// for every aligned record and `None` (unmapped) when `tid < 0`.
    fn write_single_end_records_ordered(
        path: &std::path::Path,
        records: &[(&str, Vec<u8>, u16, i32, i64)],
    ) {
        let header = grouped_alignment_header();
        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");
        for (name, seq, flags, tid, pos) in records {
            #[allow(clippy::cast_possible_truncation)]
            let cig = CigarString(vec![Cigar::Match(seq.len() as u32)]);
            let cigar = if *tid < 0 { None } else { Some(&cig) };
            add_grouped_alignment_record(
                &mut writer,
                name.as_bytes(),
                seq,
                *flags,
                *tid,
                *pos,
                -1,
                -1,
                cigar,
            );
        }
        drop(writer);
    }

    /// The SET of `(id, seq, qual)` triples a [`VecSink`] emitted as read-1 (the
    /// downstream-genotyping key the grouped/coordinate set-equality is defined
    /// on). Single-end emits are all lone (read-2 is always `None`).
    #[allow(clippy::type_complexity)]
    fn single_end_emit_set(
        sink: &VecSink,
    ) -> std::collections::HashSet<(String, Vec<u8>, Option<Vec<u8>>)> {
        sink.pairs
            .iter()
            .map(|(r1, r2)| {
                assert!(r2.is_none(), "single-end emits must be lone");
                (r1.id.clone(), r1.seq.clone(), r1.qual.clone())
            })
            .collect()
    }

    /// Test 1 (the reviewer's `sup_x` counterexample): a single-end read whose
    /// PRIMARY is off-target (SEQ=A) but whose `0x800` SUPPLEMENTARY is on-target
    /// (SEQ=B). The coordinate path scans the on-target supplementary and emits
    /// `(sup_x, B)`; the grouped path MUST emit the supplementary's B, NOT the
    /// primary's A.
    fn build_grouped_single_end_supplementary_bam(path: &std::path::Path) {
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();
        // Name-sorted: "se_pad_a" < "se_pad_b" < "sup_x". Off-target single-end
        // pads (not candidates) keep the emitted set clean and make the
        // 0x1-unset single-end majority unambiguous.
        let records: Vec<(&str, Vec<u8>, u16, i32, i64)> = vec![
            ("se_pad_a", noise(90), 0, 0, 60_000),
            ("se_pad_b", noise(90), 0, 0, 61_000),
            // sup_x primary: OFF-target, SEQ=A=ref[0..90].
            ("sup_x", ref_bytes[0..90].to_vec(), 0, 0, 50_000),
            // sup_x supplementary (0x800): ON-target, SEQ=B=ref[100..190].
            ("sup_x", ref_bytes[100..190].to_vec(), 0x800, 0, 2000),
        ];
        write_single_end_records_ordered(path, &records);
    }

    #[test]
    fn grouped_alignment_single_end_emits_on_target_supplementary_not_off_target_primary() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("sup_x.bam");
        build_grouped_single_end_supplementary_bam(&bam_path);
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let (metrics, single_end, sink) = run_grouped_alignment(&bam_path, &ref_fasta, false, true);

        assert!(single_end, "an all-0x1-unset fixture must classify single-end");
        assert_eq!(metrics.pass1_emitted, 1, "only sup_x is a candidate");
        let (r1, r2) = sink.pairs.iter().find(|(r, _)| r.id == "sup_x").unwrap();
        assert!(r2.is_none(), "single-end read emits LONE");
        assert_eq!(
            r1.seq,
            PARALLEL_TEST_REF.as_bytes()[100..190].to_vec(),
            "must emit the ON-TARGET supplementary's SEQ (B), matching the coordinate path"
        );
        assert_ne!(r1.seq, ref_bytes_first_90(), "must NOT emit the off-target primary's SEQ (A)");
    }

    #[test]
    fn grouped_alignment_single_end_on_target_read_failing_good_is_still_emitted() {
        // "Defect B" verification: on-target single-end reads are emitted after
        // ONLY the low-complexity check -- NO good()/HasHitInSet gate (T1K
        // BamExtractor.cpp:836-850). A `noise(90)` read is non-low-complexity
        // but NOT a kmer candidate, yet -- placed on-target -- the coordinate
        // path emits it, so the grouped path must too. Gating on-target behind
        // good() here would break set-equality and the coordinate regression.
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("se_ontarget_noise.bam");
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        // Document that the on-target read's SEQ genuinely fails good().
        let filter = RefKmerFilter::from_reference_fasta(&ref_fasta, INITIAL_KMER_LENGTH).unwrap();
        assert!(
            !filter.is_good_candidate(&noise(90)),
            "the noise read must NOT be a good candidate (so this proves good() is not gated)"
        );

        let records: Vec<(&str, Vec<u8>, u16, i32, i64)> = vec![
            ("se_pad_a", noise(90), 0, 0, 60_000), // off-target: not a candidate
            ("noise_on", noise(90), 0, 0, 2000),   // on-target, fails good() -> STILL emitted
        ];
        write_single_end_records_ordered(&bam_path, &records);

        let (metrics, single_end, sink) = run_grouped_alignment(&bam_path, &ref_fasta, false, true);
        assert!(single_end);
        assert_eq!(metrics.pass1_emitted, 1);
        let (r1, _) = sink.pairs.iter().find(|(r, _)| r.id == "noise_on").unwrap();
        assert_eq!(r1.seq, noise(90), "the on-target non-candidate read must be emitted verbatim");
    }

    #[test]
    fn grouped_alignment_single_end_set_equals_coordinate_path() {
        // Test 3 (STRONGEST): build ONE single-end record set, write it BOTH
        // coordinate-sorted and QNAME-grouped, run the coordinate extractor on
        // the former and the grouped extractor on the latter, and assert the
        // emitted (id, seq, qual) SETS are EQUAL. Covers: a plain on-target read;
        // the off-target-primary / on-target-supplementary shape (`sup_x`); an
        // on-target read failing good() (`noise_on`); and `multi_on`, whose two
        // on-target records have DIFFERENT coordinate order vs. file order (so it
        // pins the "first EMITTABLE in COORDINATE order" selection, not file
        // order). This previews Task 9's set-equality gate for single-end.
        let tmp = tempfile::tempdir().unwrap();
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        // The shared record set (name, seq, flags, tid, pos), single-end (0x1
        // unset), all on `chr1` (tid 0). Insertion order below is the FILE order
        // WITHIN each QNAME group for the grouped BAM.
        let seq_a = ref_bytes[0..90].to_vec(); // sup_x primary (off-target)
        let seq_b = ref_bytes[100..190].to_vec(); // sup_x supplementary (on-target)
        let seq_c = ref_bytes[200..290].to_vec(); // multi_on primary (on-target, pos 3000)
        let seq_d = ref_bytes[0..90].to_vec(); // multi_on supplementary (on-target, pos 1100)
        let records: Vec<(&str, Vec<u8>, u16, i32, i64)> = vec![
            ("plain_on", ref_bytes[0..90].to_vec(), 0, 0, 1100),
            ("noise_on", noise(90), 0, 0, 2000),
            ("sup_x", seq_a.clone(), 0, 0, 50_000),
            ("sup_x", seq_b.clone(), 0x800, 0, 2500),
            ("multi_on", seq_c.clone(), 0, 0, 3000),
            ("multi_on", seq_d.clone(), 0x800, 0, 1100),
        ];

        // Coordinate BAM: records sorted by (tid, pos) ascending (stable).
        let mut coord_records = records.clone();
        coord_records.sort_by_key(|(_, _, _, tid, pos)| (*tid, *pos));
        let coord_bam = tmp.path().join("coord.bam");
        write_single_end_records_ordered(&coord_bam, &coord_records);

        // Grouped BAM: records sorted by name (stable -> preserves file order
        // within a QNAME group).
        let mut grouped_records = records.clone();
        grouped_records.sort_by_key(|(name, _, _, _, _)| *name);
        let grouped_bam = tmp.path().join("grouped.bam");
        write_single_end_records_ordered(&grouped_bam, &grouped_records);

        let coord_sink = run_full_parallel(&coord_bam, &ref_fasta, 1);
        let (_m, grouped_single_end, grouped_sink) =
            run_grouped_alignment(&grouped_bam, &ref_fasta, false, true);
        assert!(grouped_single_end, "the grouped fixture must classify single-end");

        let coord_set = single_end_emit_set(&coord_sink);
        let grouped_set = single_end_emit_set(&grouped_sink);
        assert_eq!(
            grouped_set, coord_set,
            "grouped single-end emitted (id, seq, qual) set must equal the coordinate path's"
        );

        // Spell out the expected (id, seq) pairs so a regression on EITHER path
        // is caught, not just a coincidental agreement. (Qual is uniform here and
        // already covered by the set-equality above; comparing (id, seq) keeps
        // this assertion free of the sink's phred-encoding detail.)
        let coord_id_seq: std::collections::HashSet<(String, Vec<u8>)> =
            coord_set.iter().map(|(id, seq, _)| (id.clone(), seq.clone())).collect();
        let expected: std::collections::HashSet<(String, Vec<u8>)> = [
            ("plain_on".to_string(), ref_bytes[0..90].to_vec()),
            ("noise_on".to_string(), noise(90)),
            ("sup_x".to_string(), seq_b), // on-target supplementary B, not off-target primary A
            ("multi_on".to_string(), seq_d), // earliest-coordinate on-target D, not file-first C
        ]
        .into_iter()
        .collect();
        assert_eq!(
            coord_id_seq, expected,
            "the coordinate path's emitted (id, seq) set must match expectations"
        );
    }

    // -- Task 9: the differential set-equality GATE for PAIRED input. The
    // single-end differential above previews the same invariant; this pins the
    // dominant PAIRED case across the branch-diverse classification ladder.

    /// Writes paired records `(name, seq, flags, tid, pos, mtid, mpos)` to `path`
    /// in the EXACT order given (the caller pre-sorts: coordinate order for the
    /// coordinate extractor, QNAME-grouped order for the grouped one), under the
    /// shared `chr1`/`chr1_alt` header. CIGAR is `Match(seq.len())` for an
    /// aligned record (`tid >= 0`) and `None` (unmapped) when `tid < 0`.
    #[allow(clippy::type_complexity)]
    fn write_paired_records_ordered(
        path: &std::path::Path,
        records: &[(&str, Vec<u8>, u16, i32, i64, i32, i64)],
    ) {
        let header = grouped_alignment_header();
        let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");
        for (name, seq, flags, tid, pos, mtid, mpos) in records {
            #[allow(clippy::cast_possible_truncation)]
            let cig = CigarString(vec![Cigar::Match(seq.len() as u32)]);
            let cigar = if *tid < 0 { None } else { Some(&cig) };
            add_grouped_alignment_record(
                &mut writer,
                name.as_bytes(),
                seq,
                *flags,
                *tid,
                *pos,
                *mtid,
                *mpos,
                cigar,
            );
        }
        drop(writer);
    }

    /// The SET of `(read-1 id, mate-1 seq, mate-2 seq)` triples a [`VecSink`]
    /// emitted -- the downstream-genotyping key the grouped/coordinate
    /// set-equality is defined on for paired input.
    #[allow(clippy::type_complexity)]
    fn paired_emit_set(
        sink: &VecSink,
    ) -> std::collections::HashSet<(String, Vec<u8>, Option<Vec<u8>>)> {
        sink.pairs
            .iter()
            .map(|(r1, r2)| (r1.id.clone(), r1.seq.clone(), r2.as_ref().map(|r| r.seq.clone())))
            .collect()
    }

    #[test]
    fn grouped_alignment_paired_set_equals_coordinate_path() {
        // Build ONE paired record set, write it BOTH coordinate-sorted and
        // QNAME-grouped, run the coordinate 2-pass on the former and the grouped
        // one-pass on the latter, and assert the emitted (id, seq1, seq2) SETS are
        // EQUAL. Branch coverage: `on_target` (main-chrom aligned, overlaps the
        // gene -> candidate by overlap, NO seed test); `off_seed` (main-chrom
        // aligned OUTSIDE the gene with ref-substring SEQs that WOULD seed-hit ->
        // MUST be excluded, proving no-seed-on-main-chrom); `alt_pair` (aligned to
        // the ALT contig `chr1_alt` -> per-read `!lc && good` seed branch);
        // `unaligned_joint` (genuine template-unaligned pair, non-abnormal -> the
        // JOINT `!lc(s1)&&!lc(s2)&&(good(s1)||good(s2))` predicate).
        let tmp = tempfile::tempdir().unwrap();
        let ref_bytes = PARALLEL_TEST_REF.as_bytes();
        let ref_fasta = parallel_test_ref_fasta(tmp.path());

        let m1: u16 = 0x1 | 0x2 | 0x40; // paired, proper, first mate
        let m2: u16 = 0x1 | 0x2 | 0x80; // paired, proper, second mate
        let u1: u16 = 0x1 | 0x40 | 0x4 | 0x8; // paired, first, unmapped, mate-unmapped
        let u2: u16 = 0x1 | 0x80 | 0x4 | 0x8; // paired, second, unmapped, mate-unmapped

        // Insertion order = FILE order within each QNAME group for the grouped BAM.
        #[allow(clippy::type_complexity)]
        let records: Vec<(&str, Vec<u8>, u16, i32, i64, i32, i64)> = vec![
            ("on_target", ref_bytes[0..90].to_vec(), m1, 0, 1100, 0, 1300),
            ("on_target", ref_bytes[100..190].to_vec(), m2, 0, 1300, 0, 1100),
            ("off_seed", ref_bytes[0..90].to_vec(), m1, 0, 50_000, 0, 50_200),
            ("off_seed", ref_bytes[100..190].to_vec(), m2, 0, 50_200, 0, 50_000),
            ("alt_pair", ref_bytes[0..90].to_vec(), m1, 1, 100, 1, 300),
            ("alt_pair", ref_bytes[100..190].to_vec(), m2, 1, 300, 1, 100),
            ("unaligned_joint", ref_bytes[0..90].to_vec(), u1, -1, -1, -1, -1),
            ("unaligned_joint", ref_bytes[100..190].to_vec(), u2, -1, -1, -1, -1),
            // orphan: 0x1 SET, on-target, but its mate NEVER appears in the file
            // (a single row). The coordinate path records it in pass 1 and drops
            // it in pass 2 (no mate to complete the pair); the grouped one-pass
            // must DROP it too (a lone paired primary), or set-equality breaks --
            // this row is the regression guard for the orphan-handling fix.
            ("orphan", ref_bytes[0..90].to_vec(), m1, 0, 1150, 0, 1350),
        ];

        // Coordinate BAM: (tid, pos) ascending with UNMAPPED (tid < 0) sorting
        // LAST (real coordinate-sort order), stable so a pair's members and the
        // unaligned pair's two mates stay adjacent.
        let mut coord_records = records.clone();
        coord_records.sort_by_key(|(_, _, _, tid, pos, _, _)| {
            (if *tid < 0 { i32::MAX } else { *tid }, *pos)
        });
        let coord_bam = tmp.path().join("coord.bam");
        write_paired_records_ordered(&coord_bam, &coord_records);

        // Grouped BAM: sorted by QNAME (mates adjacent), stable within a group.
        let mut grouped_records = records.clone();
        grouped_records.sort_by_key(|(name, _, _, _, _, _, _)| *name);
        let grouped_bam = tmp.path().join("grouped.bam");
        write_paired_records_ordered(&grouped_bam, &grouped_records);

        let coord_sink = run_full_parallel(&coord_bam, &ref_fasta, 1);
        let (_m, grouped_single_end, grouped_sink) =
            run_grouped_alignment(&grouped_bam, &ref_fasta, false, true);
        assert!(!grouped_single_end, "the paired fixture must classify paired");

        assert_eq!(
            paired_emit_set(&grouped_sink),
            paired_emit_set(&coord_sink),
            "grouped paired (id, seq1, seq2) set must equal the coordinate path's"
        );

        // Spell out the expected candidate ids so a regression on EITHER path is
        // caught, not just a coincidental agreement: `off_seed` excluded (no seed
        // on the main-chrom aligned branch), the other three present.
        let ids: std::collections::HashSet<String> =
            paired_emit_set(&coord_sink).into_iter().map(|(id, _, _)| id).collect();
        let expected: std::collections::HashSet<String> =
            ["on_target", "alt_pair", "unaligned_joint"].into_iter().map(String::from).collect();
        assert_eq!(ids, expected, "off_seed (main-chrom, off-gene, seed-hitting) must be excluded");
    }
}
