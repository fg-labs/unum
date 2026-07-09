//! Thin CLI wrapper around `unum_core::genotyper` (the Rust port of `genotyper`,
//! `Genotyper.cpp`/`Genotyper.hpp`/(the `ReadAssignmentToFragmentAssignment` slice
//! of) `SeqSet.hpp`). All genotyping logic (reference dedup/loading, the read-processing loop,
//! `SelectAllelesForGenes`, output formatting) is either in `unum_core::genotyper` or -- for
//! the reference-loading/read-processing DRIVER itself, which `Genotyper.cpp:main` owns directly
//! in the C++ rather than delegating to a `Genotyper` method -- right here, matching that same
//! split of responsibility.
//!
//! # Scope: paired/single-end FASTQ, no barcode/`-a`/whitelist
//!
//! This port reproduces `Genotyper.cpp:main`'s per-read logic exactly (`Genotyper.cpp:463-480,
//! 531-574`); it does not replicate stock's own `threadCnt > 1` batching mechanics, but achieves
//! the same output-determinism property via a different means. `-t`/`--threads` parallelizes two
//! read-independent phases: (1) the read-only `get_overlaps_from_read` phase of the read-
//! assignment loop (see `unum_core::genotyper::assign_reads_parallel`'s doc comment), and
//! (2) the fragment-assembly loop, which is SLOT-INDEXED -- each read `i` writes only its own
//! `all_read_assignments[i]` slot, computed by the pure
//! `read_assignment_to_fragment_assignment` + `Genotyper::compute_read_assignment` from read
//! `i`'s inputs plus immutable shared state (`max_assign_cnt`, per-allele `whitelist`, and the
//! static `read_assignment_weight`), with NO cross-read order dependence -- so the parallel map
//! produces byte-identical slots to the serial `-t 1` loop, which are then installed in one shot
//! via `Genotyper::set_all_read_assignments`. Every genuinely order-dependent shared-state
//! mutation still runs strictly sequentially in a fixed, thread-count-independent order:
//! `assign_read`'s `allele_refs` coverage marking (interleaved within phase 1's parallel pass but
//! order-invariant), and everything downstream -- `CoalesceReadAssignments`, the
//! EM/quantification, allele selection. So `-t N` output is byte-identical to `-t 1` (and to the
//! oracle, itself always run at `-t 1` for the end-to-end differential) at any `N`.
//! `--barcode`/`-a` (a precomputed abundance file, skipping
//! `QuantifyAlleleEquivalentClass` entirely)/`--alleleWhitelist`/`--outputReadAssignment` are not
//! exposed by [`GenotypeArgs`] -- see that struct's doc comment.
//!
//! # Reference loading mirrors `Genotyper::InitRefSet` (`Genotyper.hpp:706-727`)
//!
//! [`load_reference`] deduplicates identical reference sequences (an allele appearing more than
//! once in the FASTA collapses to ONE loaded allele with `weight` = duplicate count, matching
//! `InitRefSet`'s `usedSeq`/`UpdateSeqWeight` -- this is NOT the same as `RefKmerFilter::
//! from_reference_fasta`'s own no-dedup one-record-per-FASTA-line indexing, which this driver
//! does not use directly on the raw file for that reason), parses each record's exon-interval
//! header comment via [`unum_core::genotyper::AlleleRef::new`], and computes each allele's
//! `effectiveLen` via `ComputeEffectiveLen` (`SeqSet.hpp:747-758`: sequence length with
//! consecutive `N`s collapsed to 1) -- all BEFORE building the k-mer index, so the index is over
//! exactly the deduplicated allele set `Genotyper::init_allele_info` also sees.
//!
//! [`RefKmerFilter::from_reference_fasta`] only accepts a filesystem path (not an in-memory
//! sequence list), so [`load_reference`] writes the deduplicated sequences to a securely-created
//! scratch file under [`std::env::temp_dir`] (via `tempfile`, unlinked on drop) rather than
//! widening that Phase-3/3b/4 API -- this CLI crate is exactly the layer allowed to do ordinary
//! filesystem I/O as an implementation detail.
use crate::cli::GenotypeArgs;
use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write as _;
use std::path::Path;
use unum_core::allele_freq::AlleleFreqTable;
use unum_core::fastq::{FastqReader, FastqRecord};
use unum_core::genotyper::{self, AlleleRef, ExtendedOverlap, Genotyper, ReadAssignment};
use unum_core::ref_kmer_filter::RefKmerFilter;

/// The k-mer length the genotyper's reference index/alignment is built at
/// (`Genotyper genotyper(11)`, `Genotyper.cpp:207`) -- distinct from the
/// extractor's `k=9` (`crate::stages::extract::INITIAL_KMER_LENGTH`).
const GENOTYPER_KMER_LENGTH: usize = 11;

/// `KmerCount`'s own default-constructor k-mer length
/// (`unum_core::genotyper::Genotyper::init_allele_info`'s
/// `kmer_profile_k` parameter doc comment): the gene-similarity profile
/// k-mer length, independent of [`GENOTYPER_KMER_LENGTH`].
const GENE_SIMILARITY_KMER_LENGTH: usize = 31;

/// A loaded, deduplicated reference allele: everything
/// [`unum_core::genotyper::Genotyper::init_allele_info`] and this driver's
/// read-processing loop need per allele.
struct LoadedRef {
    names: Vec<String>,
    consensus: Vec<Vec<u8>>,
    weight: Vec<i32>,
    effective_len: Vec<i32>,
    allele_refs: Vec<AlleleRef>,
}

/// Ported from `Genotyper::InitRefSet(char *filename)` (`Genotyper.hpp:706-
/// 727`): reads `path` as a FASTA (id + sequence + header comment), collapses
/// exact-duplicate sequences into one allele with an incremented `weight`
/// (mirrors `usedSeq`/`refSet.UpdateSeqWeight`), and computes each surviving
/// allele's `effectiveLen` via `ComputeEffectiveLen` (`SeqSet.hpp:747-758`).
///
/// This port writes its own minimal FASTA-with-comment reader (rather than
/// reusing [`FastqReader`], which discards the header comment -- see this
/// module's doc comment) since the exon-interval comment is load-bearing for
/// [`AlleleRef::new`].
///
/// # Errors
///
/// Returns an error if `path` cannot be read or is not valid UTF-8/FASTA.
fn load_reference(path: &Path) -> Result<LoadedRef> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading reference FASTA {}", path.display()))?;

    let mut seq_to_idx: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut names = Vec::new();
    let mut consensus: Vec<Vec<u8>> = Vec::new();
    let mut weight: Vec<i32> = Vec::new();
    let mut comments: Vec<Option<String>> = Vec::new();

    let mut cur_id: Option<String> = None;
    let mut cur_comment: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();

    let flush = |id: Option<String>,
                 comment: Option<String>,
                 seq: Vec<u8>,
                 seq_to_idx: &mut HashMap<Vec<u8>, usize>,
                 names: &mut Vec<String>,
                 consensus: &mut Vec<Vec<u8>>,
                 weight: &mut Vec<i32>,
                 comments: &mut Vec<Option<String>>| {
        let Some(id) = id else { return };
        if seq.is_empty() {
            return;
        }
        if let Some(&idx) = seq_to_idx.get(&seq) {
            weight[idx] += 1;
        } else {
            let idx = names.len();
            seq_to_idx.insert(seq.clone(), idx);
            names.push(id);
            consensus.push(seq);
            weight.push(1);
            comments.push(comment);
        }
    };

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            flush(
                cur_id.take(),
                cur_comment.take(),
                std::mem::take(&mut cur_seq),
                &mut seq_to_idx,
                &mut names,
                &mut consensus,
                &mut weight,
                &mut comments,
            );
            let mut parts = rest.splitn(2, char::is_whitespace);
            cur_id = Some(parts.next().unwrap_or("").to_string());
            cur_comment =
                parts.next().map(str::trim_start).filter(|s| !s.is_empty()).map(String::from);
        } else {
            cur_seq.extend_from_slice(line.trim_end().as_bytes());
        }
    }
    flush(
        cur_id,
        cur_comment,
        cur_seq,
        &mut seq_to_idx,
        &mut names,
        &mut consensus,
        &mut weight,
        &mut comments,
    );

    ensure!(!names.is_empty(), "reference FASTA {} contains no sequences", path.display());

    let effective_len: Vec<i32> = consensus.iter().map(|s| compute_effective_len(s)).collect();
    let allele_refs: Vec<AlleleRef> = consensus
        .iter()
        .zip(&comments)
        .map(|(seq, comment)| AlleleRef::new(seq.clone(), comment.as_deref()))
        .collect();

    Ok(LoadedRef { names, consensus, weight, effective_len, allele_refs })
}

/// Ported from `SeqSet::ComputeEffectiveLen` (`SeqSet.hpp:747-758`): sequence
/// length with consecutive `N` runs collapsed to `1`.
fn compute_effective_len(seq: &[u8]) -> i32 {
    let mut ret = 0i32;
    for (i, &b) in seq.iter().enumerate() {
        if b != b'N' || (i > 0 && seq[i - 1] != b'N') {
            ret += 1;
        }
    }
    ret
}

/// One read (mate) loaded from a candidate FASTQ, mirroring the fields of
/// `_genotypeRead` (`Genotyper.cpp:59-81`) this driver actually needs. `id`
/// is written (matching `_genotypeRead::id`) to `{prefix}_aligned{,_1,_2}.fa`
/// for each assigned fragment (`Genotyper.cpp:688-703`); `seq` is used both
/// for read-end alignment and for that same aligned-FASTA output (emitted in
/// its original input orientation, never reverse-complemented).
struct GenotypeRead {
    id: String,
    seq: Vec<u8>,
    has_n: bool,
}

impl GenotypeRead {
    /// Builds a [`GenotypeRead`] from a parsed [`FastqRecord`], recomputing
    /// `has_n = seq.contains(&b'N')` (`Genotyper.cpp:59-81`). Used for both the
    /// file-based path (records read from `_candidate_*.fq`) and the fused
    /// in-memory path (records handed over by
    /// [`unum_core::extract::InMemoryCandidateSink`]) -- identical either way.
    fn from_record(record: FastqRecord) -> Self {
        let has_n = record.seq.contains(&b'N');
        Self { id: record.id, seq: record.seq, has_n }
    }
}

/// Runs the `genotype` subcommand for `args`.
///
/// # Errors
///
/// Returns an error if the reference/candidate-read files cannot be
/// opened/parsed, if neither `-u` nor both `-1`/`-2` are given, or if output
/// files cannot be created.
pub fn run(args: &GenotypeArgs) -> Result<()> {
    let mate2_path = args.mate2.as_deref();
    let paired = args.mate1.is_some() || mate2_path.is_some();
    let single_path = args.single.as_deref();
    ensure!(
        !(paired && single_path.is_some()),
        "specify either -u (single-end) or -1/-2 (paired), not both"
    );
    let (mate1_path, mate2_path): (&str, Option<&str>) = if paired {
        let mate1 = args
            .mate1
            .as_deref()
            .context("paired input requires both -1 and -2 (got -2 without -1)")?;
        let mate2 =
            mate2_path.context("paired input requires both -1 and -2 (got -1 without -2)")?;
        (mate1, Some(mate2))
    } else if let Some(single) = single_path {
        (single, None)
    } else {
        bail!("must specify either -u (single-end) or -1/-2 (paired) candidate-read input");
    };
    let has_mate = mate2_path.is_some();

    // --- Candidate FASTQ reading (Genotyper.cpp:362-443) ---
    // Each mate is read in its own file order into its own `Vec<FastqRecord>`,
    // then handed to [`run_with_candidate_reads`] (the shared driver from
    // reference-load onward), which is exactly what the fused in-memory `run`
    // path (issue #28) also feeds -- so the file-based and in-memory paths run
    // byte-identical downstream logic on byte-identical read vectors.
    let reads1 = read_records(Path::new(mate1_path))
        .with_context(|| format!("reading candidate read file {mate1_path}"))?;
    let reads2 = if let Some(mate2_path) = mate2_path {
        read_records(Path::new(mate2_path))
            .with_context(|| format!("reading candidate read file {mate2_path}"))?
    } else {
        Vec::new()
    };

    run_with_candidate_reads(args, reads1, reads2, has_mate)
}

/// Reads every record from the FASTA/FASTQ at `path` into a `Vec<FastqRecord>`
/// (id whitespace-truncated + `/1`/`/2`-stripped, seq/qual verbatim -- exactly
/// what the fused in-memory candidate sink also produces after extract has
/// already normalized ids on the way in; see [`run_with_candidate_reads`]).
fn read_records(path: &Path) -> Result<Vec<FastqRecord>> {
    let mut reader =
        FastqReader::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut out = Vec::new();
    while let Some(record) = reader.next_record()? {
        out.push(record);
    }
    Ok(out)
}

/// Runs the genotyper driver from reference-load onward on candidate reads
/// already present in memory as `Vec<FastqRecord>` -- the seam the fused
/// single-process `run` pipeline (issue #28) feeds directly from
/// [`unum_core::extract::InMemoryCandidateSink`], with NO intermediate
/// `{prefix}_candidate_*.fq` written and re-read.
///
/// `reads1`/`reads2` are the mate-1/mate-2 candidate records in strict input
/// order (index-aligned: `reads1[i]`/`reads2[i]` are the two ends of fragment
/// `i`). For paired input (`has_mate == true`) both vectors must have equal
/// length; for single-end input `reads2` is empty. Every mate-2 record must
/// already carry its pair's MATE-1 id (both the file path -- via extract's
/// `_candidate_2.fq` writing mate-1's id and genotype re-reading it -- and the
/// in-memory `InMemoryCandidateSink` guarantee this), so `_aligned_2.fa`'s
/// `>{id}` lines are byte-identical across both paths.
///
/// This is byte-identical to the previous inlined driver: it performs the same
/// reference load, builds the same `GenotypeRead`s (recomputing `has_n` as
/// `seq.contains(&b'N')`), and runs the identical alignment/fragment-assembly/
/// quantify/select/output sequence.
///
/// # Errors
///
/// Returns an error if the `has_mate`/`reads2` invariant is violated (paired
/// mate counts differ, or a single-end call carries mate-2 records), if the
/// reference cannot be loaded/parsed, or if output files cannot be created.
// Long by construction: mirrors Genotyper.cpp:main's single linear driver
// (reference load -> read load -> alignment -> fragment assembly ->
// quantify/select -> output) -- splitting this into several tiny functions
// would scatter closely-related sequential setup across the file for no
// readability benefit (same rationale as the vendored source's own
// single-function `main`).
#[allow(clippy::too_many_lines)]
pub fn run_with_candidate_reads(
    args: &GenotypeArgs,
    reads1: Vec<FastqRecord>,
    reads2: Vec<FastqRecord>,
    has_mate: bool,
) -> Result<()> {
    // Enforce the `has_mate`/`reads2` invariant at the seam, before any
    // reference load or shared-sequence processing: `has_mate == false` with a
    // non-empty `reads2` is an inconsistent single-end call whose mate-2
    // sequences would still feed shared read processing yet be dropped during
    // fragment assembly/output. In-tree callers always uphold this, but this is
    // a `pub` entrypoint, so guard it here.
    if has_mate {
        ensure!(
            reads1.len() == reads2.len(),
            "mate-1 ({}) and mate-2 ({}) candidate read counts differ",
            reads1.len(),
            reads2.len()
        );
    } else {
        ensure!(reads2.is_empty(), "single-end candidate reads must not include mate-2 records");
    }

    // --- Reference loading (Genotyper.hpp:706-727) ---
    let loaded = load_reference(Path::new(&args.ref_seq_fasta))?;
    let mut filter = build_ref_kmer_filter(&loaded.names, &loaded.consensus)?;
    filter.set_ref_seq_similarity(args.similarity);

    let mut genotyper = Genotyper::new();
    genotyper.set_filter_frac(args.filter_frac);
    genotyper.set_filter_cov(args.filter_cov);
    genotyper.set_cross_gene_rate(args.cross_gene_rate);
    let mut effective_len = loaded.effective_len.clone();
    genotyper.init_allele_info(
        &loaded.names,
        &loaded.consensus,
        &loaded.weight,
        &mut effective_len,
        GENE_SIMILARITY_KMER_LENGTH,
    );

    // Build the per-read `GenotypeRead`s (`has_n = seq.contains(&b'N')`,
    // `Genotyper.cpp:59-81`) from the provided records, mate 1 then mate 2.
    let reads1: Vec<GenotypeRead> = reads1.into_iter().map(GenotypeRead::from_record).collect();
    let reads2: Vec<GenotypeRead> = reads2.into_iter().map(GenotypeRead::from_record).collect();

    let allele_cnt = loaded.names.len();
    let read_cnt = reads1.len();
    let max_read_length =
        reads1.iter().chain(reads2.iter()).map(|r| r.seq.len()).max().unwrap_or(0);
    genotyper.set_read_length(i32::try_from(max_read_length).unwrap_or(0));

    // --- Read-end alignment, reusing identical sequences (Genotyper.cpp:445-480) ---
    // `overlaps_by_seq` mirrors `AssignReads_Thread`/the single-threaded loop
    // at `Genotyper.cpp:463-480`: sort ALL read ends (mate 1 + mate 2) by
    // sequence, and call AssignRead once per distinct sequence, reusing the
    // result for every read end sharing that exact sequence. `-t`/`--threads`
    // parallelizes the dominant `get_overlaps_from_read` cost across
    // `sorted_seqs`; `assign_read`'s `allele_refs` mutation always runs
    // sequentially in `sorted_seqs` order, so output is byte-identical at any
    // thread count -- see `genotyper::assign_reads_parallel`'s doc comment.
    // Not `mut`: `assign_reads_parallel` mutates each allele's `pos_weight`
    // base-coverage counters through interior mutability (`AtomicPosWeight`),
    // so it takes `&[AlleleRef]` -- the shared borrow is what lets the fused
    // pass mark coverage across `rayon` workers.
    let allele_refs = loaded.allele_refs;
    let mut all_seqs: Vec<&[u8]> = reads1.iter().map(|r| r.seq.as_slice()).collect();
    all_seqs.extend(reads2.iter().map(|r| r.seq.as_slice()));
    let mut sorted_seqs = all_seqs.clone();
    sorted_seqs.sort_unstable();
    sorted_seqs.dedup();

    let mut counted: HashMap<&[u8], i32> = HashMap::new();
    for &seq in &all_seqs {
        *counted.entry(seq).or_insert(0) += 1;
    }
    let threads = usize::try_from(args.threads).unwrap_or(usize::MAX).max(1);
    let extended_by_seq = genotyper::assign_reads_parallel(
        &filter,
        &sorted_seqs,
        &allele_refs,
        args.similarity,
        |seq| counted[seq],
        threads,
    );
    let mut overlaps_by_seq: HashMap<&[u8], Option<Vec<ExtendedOverlap>>> = HashMap::new();
    for (&seq, extended) in sorted_seqs.iter().zip(extended_by_seq) {
        overlaps_by_seq.insert(seq, extended);
    }

    // --- Fragment assembly + SetReadAssignments (Genotyper.cpp:531-574) ---
    genotyper.init_read_assignments(i32::try_from(read_cnt).unwrap_or(0), args.max_assign_cnt);
    let consensus_len_of = |idx: u32| {
        i32::try_from(allele_refs[usize::try_from(idx).unwrap_or(0)].consensus.len()).unwrap_or(0)
    };
    let hit_len_required = filter.hit_len_required();
    let ref_seq_similarity = filter.ref_seq_similarity();
    // `SeqSet::IsSeparatorInRange(s, e, seqIdx)` (`SeqSet.hpp:487-498`), fed
    // by each allele's `AlleleRef::separator` (built at load time by
    // `AlleleRef::new`, mirroring `_seqWrapper::separator`) -- see
    // `unum_core::genotyper::is_separator_in_range`'s doc comment (fixes #39:
    // this was previously stubbed as `|_, _, _| false`, over-assigning reads
    // whose alignment span crosses an interior reference `N`).
    let separator_in_range = |s: i32, e: i32, seq_idx: i32| {
        let idx = usize::try_from(seq_idx).unwrap_or(0);
        genotyper::is_separator_in_range(&allele_refs[idx].separator, s, e)
    };
    // `Genotyper::set_read_assignments`'s own `separator_lookup` closure is
    // called as `(seq_idx, s, e)` (see that function's doc comment) --
    // reordered here rather than at the call site.
    let separator_lookup_for_set_read_assignments = |seq_idx: i32, s: i32, e: i32| {
        let idx = usize::try_from(seq_idx).unwrap_or(0);
        genotyper::is_separator_in_range(&allele_refs[idx].separator, s, e)
    };

    // Fragment assembly is parallelized across reads and byte-identical to
    // `-t 1` because it is SLOT-INDEXED: each read `i`'s slot is computed by
    // the pure `read_assignment_to_fragment_assignment` +
    // `Genotyper::compute_read_assignment` from read `i`'s inputs plus
    // immutable shared state (`&genotyper`, `overlaps_by_seq`, the read
    // vectors, and the `Fn + Sync` lookup closures) -- there is NO cross-read
    // order dependence, so `into_par_iter().collect()` produces the same
    // `Vec<Vec<ReadAssignment>>` in the same slot order regardless of thread
    // count. The slots are then moved into `all_read_assignments` in one shot
    // (`set_all_read_assignments`). We run inside a scoped pool sized to `-t`
    // so `-t 1` is strictly single-threaded (and still byte-identical).
    let fragment_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("building rayon thread pool for parallel fragment assembly");
    // Alongside each read's fragment-assignment slot, capture whether the
    // fragment was assigned to >= 1 allele (`fragmentAssignment.size() > 0`,
    // `Genotyper.cpp:564-565`) for the `_aligned*.fa` output below. Like the
    // slots, `fragment_assigned[i]` is a pure function of read `i`'s inputs
    // (`!assignment.is_empty()`), so it is collected in the same deterministic
    // slot order regardless of thread count -- byte-identical to `-t 1`.
    let (fragment_assigned, slots): (Vec<bool>, Vec<Vec<ReadAssignment>>) =
        fragment_pool.install(|| {
            (0..read_cnt)
                .into_par_iter()
                .map(|i| {
                    let overlaps1 =
                        overlaps_by_seq.get(reads1[i].seq.as_slice()).and_then(Option::as_ref);
                    let overlaps2 = if has_mate {
                        overlaps_by_seq.get(reads2[i].seq.as_slice()).and_then(Option::as_ref)
                    } else {
                        None
                    };
                    let has_n = reads1[i].has_n || (has_mate && reads2[i].has_n);

                    let empty: Vec<ExtendedOverlap> = Vec::new();
                    let assignment = genotyper::read_assignment_to_fragment_assignment(
                        overlaps1.map_or(empty.as_slice(), Vec::as_slice),
                        if has_mate {
                            Some(overlaps2.map_or(empty.as_slice(), Vec::as_slice))
                        } else {
                            None
                        },
                        has_n,
                        hit_len_required,
                        consensus_len_of,
                        separator_in_range,
                    );
                    let assigned = !assignment.is_empty();
                    let slot = genotyper.compute_read_assignment(
                        &assignment,
                        ref_seq_similarity,
                        separator_lookup_for_set_read_assignments,
                    );
                    (assigned, slot)
                })
                .unzip()
        });
    genotyper.set_all_read_assignments(slots);

    genotyper.coalesce_read_assignments(0, i32::try_from(read_cnt.saturating_sub(1)).unwrap_or(0));

    // --- FinalizeReadAssignments + quantification + selection (Genotyper.cpp:634-650) ---
    // Per-allele independent (each reads only its own `pos_weight`, which is
    // fully accumulated once `assign_reads_parallel` above has joined), so this
    // is embarrassingly parallel and byte-identical to the serial map. Run it
    // inside `fragment_pool` (not Rayon's default global pool) so `-t` bounds
    // this stage's worker count exactly as it does the fragment pass above --
    // otherwise `-t 1` would silently parallelize this loop.
    let missing_coverage: Vec<i32> = fragment_pool.install(|| {
        (0..allele_cnt)
            .into_par_iter()
            .map(|i| genotyper::get_seq_missing_base_coverage(&allele_refs[i], 0.01))
            .collect()
    });
    genotyper.finalize_read_assignments(&missing_coverage);

    genotyper.quantify_allele_equivalent_class(&effective_len, &loaded.weight);
    genotyper.remove_low_likelihood_allele_in_equivalent_class(|idx| effective_len[idx]);

    // --- #29 opt-in Hardy-Weinberg population-frequency prior ---
    // Configured immediately before selection: when `--allele-freq` is absent
    // the table stays `None` and selection is byte-identical to the oracle. The
    // weight and null-penalty are always propagated (they are inert at their
    // `0.0`/default values; the null penalty is name-driven, so it can act even
    // without a table, but only once set > 0).
    if let Some(freq_path) = args.allele_freq.as_deref() {
        let table = AlleleFreqTable::from_tsv(Path::new(freq_path)).with_context(|| {
            format!("loading allele-frequency table for the HWE prior from {freq_path}")
        })?;
        genotyper.set_allele_freq(table);
    }
    genotyper.set_allele_freq_weight(args.allele_freq_weight);
    genotyper.set_allele_freq_null_penalty(args.allele_freq_null_penalty);

    genotyper.select_alleles_for_genes(|idx| loaded.weight[idx]);

    // --- Output (Genotyper.cpp:652-707) ---
    write_genotype_tsv(&genotyper, &args.prefix)?;
    genotyper
        .output_representative_alleles(Path::new(&format!("{}_allele.tsv", args.prefix)), |idx| {
            loaded.names[idx].clone()
        });
    write_aligned_fastas(&args.prefix, has_mate, &fragment_assigned, &reads1, &reads2)?;

    // Purely additive opt-in QC + discriminative-quality panel (issues #19/#30).
    // Reads only already-computed genotyper state + `allele_refs` coverage;
    // never touches the byte-frozen `_genotype.tsv`/`_allele.tsv` above.
    if args.emit_metrics {
        write_metrics_tsv(&genotyper, &allele_refs, &loaded.names, &args.prefix)?;
    }

    eprintln!(
        "genotyped {read_cnt} read fragments across {allele_cnt} alleles / {} genes (average {:.2} \
         alleles/read)",
        genotyper.gene_cnt,
        genotyper.get_average_read_assignment_cnt(),
    );

    Ok(())
}

/// Builds a [`RefKmerFilter`] over `names`/`consensus` (the deduplicated
/// allele set from [`load_reference`]) at [`GENOTYPER_KMER_LENGTH`].
/// [`RefKmerFilter::from_reference_fasta`] only accepts a path (see this
/// module's doc comment for why), so this writes a scratch FASTA to
/// [`std::env::temp_dir`] and removes it before returning.
fn build_ref_kmer_filter(names: &[String], consensus: &[Vec<u8>]) -> Result<RefKmerFilter> {
    // Use a securely-created temp file (random name + O_EXCL, via `tempfile`)
    // rather than a predictable `temp_dir()/pid-counter` path: a predictable
    // name in a shared temp dir can be pre-created or symlink-swapped by another
    // local user. `NamedTempFile` is unlinked on drop, so it cleans up on every
    // return path (including the error path below).
    let mut scratch = tempfile::Builder::new()
        .prefix("unum-genotype-ref-")
        .suffix(".fa")
        .tempfile()
        .context("creating scratch reference FASTA")?;
    for (name, seq) in names.iter().zip(consensus) {
        writeln!(scratch, ">{name}")?;
        scratch.write_all(seq)?;
        writeln!(scratch)?;
    }
    scratch.flush().context("flushing scratch reference FASTA")?;

    RefKmerFilter::from_reference_fasta(scratch.path(), GENOTYPER_KMER_LENGTH)
        .with_context(|| format!("building k-mer index over {} deduplicated alleles", names.len()))
}

/// Ported from `Genotyper.cpp:658-670`: writes `{prefix}_genotype.tsv`, one
/// row per gene, formatted as `{geneName}\t{calledAlleleCnt}\t{allele1}\t
/// {allele2}\t{secondaryAlleles}\n` (each of the three allele fields itself
/// already contains embedded `\t`/`;`-separated abundance/quality, from
/// [`Genotyper::get_allele_description`]).
fn write_genotype_tsv(genotyper: &Genotyper, prefix: &str) -> Result<()> {
    let path = format!("{prefix}_genotype.tsv");
    let mut out = std::fs::File::create(&path).with_context(|| format!("creating {path}"))?;
    let gene_cnt = usize::try_from(genotyper.gene_cnt).unwrap_or(0);
    for i in 0..gene_cnt {
        let (allele1, allele2, secondary, called_cnt) = genotyper.get_allele_description(i);
        writeln!(
            out,
            "{}\t{called_cnt}\t{allele1}\t{allele2}\t{secondary}",
            genotyper.gene_idx_to_name[i]
        )
        .with_context(|| format!("writing {path}"))?;
    }
    Ok(())
}

/// Writes the opt-in `{prefix}_metrics.tsv` per-call QC + discriminative-quality
/// panel (issues #19/#30 -- a unum extension, NOT part of T1K). One header
/// line plus one row per CALLED allele: the representative allele of rank 0 and
/// (for a heterozygous call) rank 1 of each called gene, matching the
/// representative-per-rank selection [`Genotyper::output_representative_alleles`]
/// uses for `_allele.tsv` (highest `ec_abundance` member of the rank, ties
/// broken by the LOWER allele index). A homozygous, rank-0-only call emits a
/// single row. Columns (tab-separated, in order):
///
/// `gene`, `allele_rank`, `allele`, `abundance`, `balance_ratio`, `cov_min`,
/// `cov_p10`, `cov_median`, `frac_bases_covered`, `missing_cov`, `gt_quality`,
/// `runnerup_abundance`, `q_gap`, `q_min`, `locus_gq_min`, `ec_set_size`,
/// `identifiability`, `ec_ambiguity_entropy`, `series_set`.
///
/// The final four columns (issue #33, report-only) describe the called
/// allele's post-prune equivalence class: `ec_set_size` is its member count
/// (1 = read-distinguishable), `identifiability` is `1 / ec_set_size`
/// (∈ (0,1]), `ec_ambiguity_entropy` is the Shannon entropy over the EC
/// members' `ec_abundance`, and `series_set` is the semicolon-joined
/// major-allele series of the EC members. Like the rest of this panel they are
/// read from already-retained genotyper state and never touch the byte-frozen
/// `_genotype.tsv`/`_allele.tsv`.
///
/// Every value is read from ALREADY-COMPUTED genotyper state or derived from
/// `allele_refs`' retained per-base coverage ([`genotyper::allele_coverage_stats`]);
/// nothing here is recomputed and nothing touches the byte-frozen
/// `_genotype.tsv`/`_allele.tsv`. `balance_ratio`, `runnerup_abundance`, and
/// `locus_gq_min` are per-GENE and repeated on each of the gene's rows.
/// `balance_ratio`/`runnerup_abundance` use the per-rank SUMMED abundances (the
/// same rank grouping as `select_alleles_for_genes_quality_scores`); the
/// `abundance` column is the representative allele's own abundance, formatted
/// with the same `%.6` precision as [`Genotyper::get_allele_description`].
fn write_metrics_tsv(
    genotyper: &Genotyper,
    allele_refs: &[AlleleRef],
    names: &[String],
    prefix: &str,
) -> Result<()> {
    let path = format!("{prefix}_metrics.tsv");
    let mut out = std::fs::File::create(&path).with_context(|| format!("creating {path}"))?;
    writeln!(
        out,
        "gene\tallele_rank\tallele\tabundance\tbalance_ratio\tcov_min\tcov_p10\tcov_median\t\
         frac_bases_covered\tmissing_cov\tgt_quality\trunnerup_abundance\tq_gap\tq_min\t\
         locus_gq_min\tec_set_size\tidentifiability\tec_ambiguity_entropy\tseries_set"
    )
    .with_context(|| format!("writing {path}"))?;

    let gene_cnt = usize::try_from(genotyper.gene_cnt).unwrap_or(0);
    for gene_idx in 0..gene_cnt {
        let selected = &genotyper.selected_alleles[gene_idx];

        // Per-gene aggregates over the selected alleles (same rank grouping as
        // `select_alleles_for_genes_quality_scores`'s `allele_rank_abund`):
        // rank-0/1 summed abundances (allele balance), rank-2 summed abundance
        // (the `a_second` runner-up used by `q_gap`). Also pick each rank's
        // REPRESENTATIVE allele (highest `ec_abundance`, ties -> lower allele
        // index), mirroring `output_representative_alleles`.
        let mut abund_r0 = 0.0f64;
        let mut abund_r1 = 0.0f64;
        let mut abund_r2 = 0.0f64;
        let mut has_r1 = false;
        let mut reps: [Option<i32>; 2] = [None, None];
        for &(allele_idx, rank) in selected {
            let info = &genotyper.allele_info[usize::try_from(allele_idx).unwrap()];
            match rank {
                0 => abund_r0 += info.abundance,
                1 => {
                    abund_r1 += info.abundance;
                    has_r1 = true;
                }
                2 => abund_r2 += info.abundance,
                _ => {}
            }
            if rank == 0 || rank == 1 {
                let slot = usize::try_from(rank).unwrap();
                // Exact `==` on `ec_abundance` mirrors `output_representative_alleles`'s
                // own exact tie-break (same representative-per-rank selection), not a
                // fuzzy float comparison.
                #[allow(clippy::float_cmp)]
                let better = match reps[slot] {
                    None => true,
                    Some(cur) => {
                        let cur_ec =
                            genotyper.allele_info[usize::try_from(cur).unwrap()].ec_abundance;
                        info.ec_abundance > cur_ec
                            || (info.ec_abundance == cur_ec && allele_idx < cur)
                    }
                };
                if better {
                    reps[slot] = Some(allele_idx);
                }
            }
        }

        // Het balance = min/max of the two rank abundances, in (0,1];
        // homozygous / single-rank call = 1.0.
        let balance_ratio = if has_r1 {
            let (lo, hi) =
                if abund_r0 <= abund_r1 { (abund_r0, abund_r1) } else { (abund_r1, abund_r0) };
            if hi > 0.0 { lo / hi } else { 1.0 }
        } else {
            1.0
        };
        let runnerup_abundance = abund_r2;

        // Per-locus minimum null-model quality across the present rank reps.
        let locus_gq = [reps[0], reps[1]]
            .into_iter()
            .flatten()
            .map(|idx| genotyper.allele_info[usize::try_from(idx).unwrap()].genotype_quality)
            .min()
            .unwrap_or(-1);

        for (rank, rep) in reps.iter().enumerate() {
            let Some(allele_idx) = *rep else { continue };
            let idx = usize::try_from(allele_idx).unwrap();
            let info = &genotyper.allele_info[idx];
            let stats = genotyper::allele_coverage_stats(&allele_refs[idx]);
            let q_gap = info.discriminative_quality;
            let q_min = info.genotype_quality.min(q_gap);
            // Issue #33 report-only identifiability columns, derived from the
            // called allele's already-retained post-prune equivalence class.
            let ec_set_size = genotyper.ec_set_size(idx);
            #[allow(clippy::cast_precision_loss)]
            let identifiability = 1.0 / ec_set_size as f64;
            let ec_ambiguity_entropy = genotyper.ec_ambiguity_entropy(idx);
            let series_set = genotyper.ec_series_set(idx);
            writeln!(
                out,
                "{gene}\t{rank}\t{allele}\t{abundance:.6}\t{balance_ratio:.6}\t{cov_min}\t\
                 {cov_p10}\t{cov_median}\t{frac:.6}\t{missing}\t{gt}\t{runnerup:.6}\t{q_gap}\t\
                 {q_min}\t{locus_gq}\t{ec_set_size}\t{identifiability:.6}\t\
                 {ec_ambiguity_entropy:.6}\t{series_set}",
                gene = genotyper.gene_idx_to_name[gene_idx],
                allele = names[idx],
                abundance = info.abundance,
                cov_min = stats.cov_min,
                cov_p10 = stats.cov_p10,
                cov_median = stats.cov_median,
                frac = stats.frac_covered,
                missing = info.missing_coverage,
                gt = info.genotype_quality,
                runnerup = runnerup_abundance,
            )
            .with_context(|| format!("writing {path}"))?;
        }
    }
    Ok(())
}

/// Ported from `Genotyper.cpp:680-707`: writes each assigned fragment's read
/// sequence(s) as FASTA so `unum analyze` (and T1K's `analyzer`) can consume
/// them. Paired input produces `{prefix}_aligned_1.fa` + `{prefix}_aligned_2.fa`;
/// single-end produces `{prefix}_aligned.fa`.
///
/// A fragment is emitted iff it was assigned to at least one allele
/// (`fragment_assigned[i]`, i.e. its `read_assignment_to_fragment_assignment`
/// was non-empty -- `Genotyper.cpp:564-565`). The SAME fragment-level flag
/// gates BOTH mate files (`reads1[i].fragmentAssigned` for mate 1 AND mate 2,
/// `Genotyper.cpp:688,701`), so mate-1/mate-2 records stay positionally paired.
fn write_aligned_fastas(
    prefix: &str,
    has_mate: bool,
    fragment_assigned: &[bool],
    reads1: &[GenotypeRead],
    reads2: &[GenotypeRead],
) -> Result<()> {
    let path1 =
        if has_mate { format!("{prefix}_aligned_1.fa") } else { format!("{prefix}_aligned.fa") };
    write_aligned_fasta(&path1, fragment_assigned, reads1)?;
    if has_mate {
        write_aligned_fasta(&format!("{prefix}_aligned_2.fa"), fragment_assigned, reads2)?;
    }
    Ok(())
}

/// Writes one aligned-read FASTA: for each `fragment_assigned[i]`, a
/// `>{id}\n{seq}\n` record (matching `fprintf(fpOutput, ">%s\n%s\n", ...)`,
/// `Genotyper.cpp:690`), in original input orientation.
fn write_aligned_fasta(
    path: &str,
    fragment_assigned: &[bool],
    reads: &[GenotypeRead],
) -> Result<()> {
    let mut f = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
    for (i, &assigned) in fragment_assigned.iter().enumerate() {
        if assigned {
            writeln!(f, ">{}", reads[i].id).with_context(|| format!("writing {path}"))?;
            f.write_all(&reads[i].seq).with_context(|| format!("writing {path}"))?;
            writeln!(f).with_context(|| format!("writing {path}"))?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compute_effective_len_collapses_consecutive_n_runs() {
        assert_eq!(compute_effective_len(b"ACGT"), 4);
        assert_eq!(compute_effective_len(b"ACNNNNGT"), 5); // NNNN collapses to 1
        assert_eq!(compute_effective_len(b"NNAC"), 2); // leading NN collapses to 1, but that 1 is never counted (i==0 case is always excluded when seq[0]=='N')
        assert_eq!(compute_effective_len(b"ACNNGTNN"), 6);
    }

    /// Builds a [`GenotypeArgs`] whose `ref_seq_fasta` points at a path that
    /// does not exist, so a test that reaches `load_reference` fails with a
    /// reference-load error rather than the seam guard we are exercising --
    /// letting the two be told apart by error text.
    fn args_with_missing_reference() -> GenotypeArgs {
        GenotypeArgs {
            ref_seq_fasta: "/nonexistent/reference/for/seam/guard/test.fa".to_string(),
            mate1: None,
            mate2: None,
            single: None,
            prefix: "t1k".to_string(),
            threads: 1,
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

    fn record(id: &str, seq: &[u8]) -> FastqRecord {
        FastqRecord { id: id.to_string(), seq: seq.to_vec(), qual: None }
    }

    #[test]
    fn run_with_candidate_reads_rejects_single_end_carrying_mate2() {
        // `has_mate == false` but a non-empty `reads2` is an inconsistent
        // single-end call: the mate-2 sequences would leak into shared read
        // processing yet be dropped from fragment assembly. The seam must
        // reject it up front, before any reference load is attempted.
        let args = args_with_missing_reference();
        let err = run_with_candidate_reads(
            &args,
            vec![record("r1", b"ACGT")],
            vec![record("r1", b"ACGT")],
            false,
        )
        .expect_err("single-end call with non-empty reads2 must error");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("single-end") && msg.contains("mate-2"),
            "expected the seam guard's single-end message, got: {msg}"
        );
    }

    #[test]
    fn run_with_candidate_reads_rejects_paired_count_mismatch() {
        // Mismatched paired mate counts must also be rejected at the seam,
        // before the reference is loaded.
        let args = args_with_missing_reference();
        let err = run_with_candidate_reads(
            &args,
            vec![record("r1", b"ACGT"), record("r2", b"ACGT")],
            vec![record("r1", b"ACGT")],
            true,
        )
        .expect_err("paired call with differing mate counts must error");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("candidate read counts differ"),
            "expected the seam guard's paired-mismatch message, got: {msg}"
        );
    }

    #[test]
    fn load_reference_dedups_identical_sequences_and_sums_weight() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            ">A*01:01 1 3\n\
             ACGT\n\
             >A*01:02 1 3\n\
             ACGT\n\
             >B*07:02 1 3\n\
             TTTT\n",
        )
        .unwrap();

        let loaded = load_reference(tmp.path()).unwrap();
        assert_eq!(loaded.names.len(), 2, "identical ACGT sequences must collapse to one allele");
        let a_idx = loaded.names.iter().position(|n| n == "A*01:01").unwrap();
        assert_eq!(loaded.weight[a_idx], 2, "duplicate sequence increments weight");
        let b_idx = loaded.names.iter().position(|n| n == "B*07:02").unwrap();
        assert_eq!(loaded.weight[b_idx], 1);
    }
}
