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
//! # Threading: this port is single-threaded-only, matching the oracle's `-t 1`
//!
//! Unlike `FastqExtractor.cpp` (where multi-threading only parallelizes the
//! per-read filter DECISION, not the emission order -- see
//! [`crate::extract`]'s module docs for why that makes single-threaded Rust
//! provably output-identical to the oracle at ANY `-t`),
//! `BamExtractor.cpp`'s `threadCnt > 1` path changes BATCHING/OUTPUT-QUEUE
//! FLUSH TIMING for the unaligned-template-pair and single-end-unmapped
//! candidate paths (`ProcessUnmappedReads_Thread` +
//! `DistributeWork`/`AddWorkQueue`'s `workLoad = 2048` batching,
//! `BamExtractor.cpp:202-407`): candidates are queued up to 2048 at a time,
//! handed to a free worker thread, and flushed to `fp1`/`fp2` only once the
//! shared `outputQueue` exceeds `2 * candidates.size()`
//! (`BamExtractor.cpp:243`) -- i.e. output order and batching are NOT
//! provably thread-count-invariant the way `FastqExtractor.cpp`'s is. This
//! port therefore only reproduces the `threadCnt == 1` code path
//! (`BamExtractor.cpp:675-696,754-778`: the direct single-threaded
//! `if`/`else` branches, never the `AddWorkQueue`/`ProcessUnmappedReads_Thread`
//! branches), and the differential test always runs the oracle at `-t 1` to
//! match -- see `crates/unum-core/tests/golden_bam_extract.rs`'s module docs.
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
use crate::extract::{CandidateSink, ReadRecord};
use crate::ref_kmer_filter::{RefKmerFilter, Scratch, is_low_complexity};
use anyhow::{Context, Result, bail};
use std::collections::HashMap;

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

    let (candidates, pass1_emitted) = run_pass1(
        alignments,
        filter,
        genes,
        single_end,
        abnormal_unaligned_flag,
        mate_id_len,
        sink,
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
/// # Errors
///
/// Returns an error if the second record is missing (EOF) or its trimmed
/// read id does not match the first record's (mirrors
/// `BamExtractor.cpp:657-672`'s "Two reads from the unaligned fragment are
/// not showing up together" error).
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

/// Pass 1 (`BamExtractor.cpp:632-851`): see [`extract_from_bam`]'s and this
/// module's doc comments for the full branch structure. Returns the
/// `candidates` map recorded for pass 2 (empty for single-end input, which
/// never populates it) and the number of pairs/reads emitted directly during
/// this pass.
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
}
