//! Thin CLI wrapper around `unum_core::genotyper`/`unum_core::variant_caller` (the Rust port
//! of `analyzer`, `vendor/t1k/Analyzer.cpp`). All post-analysis logic (restricted reference
//! loading, read reassignment, pileup/`ComputeVariant`) is either in those two modules or -- for
//! the read-I/O/reassignment DRIVER itself, which `Analyzer.cpp:main` owns directly rather than
//! delegating to a `Genotyper`/`VariantCaller` method -- right here, matching
//! `crate::stages::genotype`'s same split of responsibility.
//!
//! # Scope: `-t N`-parallel, paired/single-end aligned-read FASTA, no barcode
//!
//! This port reproduces `Analyzer.cpp:main`'s per-read logic exactly
//! (`Analyzer.cpp:467-484,528-565,623-634`); like [`crate::stages::genotype`], it does not
//! replicate stock's own `threadCnt > 1` batching mechanics but achieves the same output-
//! determinism property by parallelizing three read-independent per-read loops in a way that is
//! byte-identical to `-t 1` (and to the oracle, itself always run at `-t 1` for the end-to-end
//! differential) at any thread count: (0) the dedup `assign_read` loop (via
//! [`unum_core::genotyper::assign_reads_parallel`] with `weight = 0`, so its coverage marking is
//! a no-op -- there is no sequential coverage-marking hazard on the analyzer path); (A) the
//! slot-indexed fragment-assembly loop (each read `i` computes its own
//! `(fragment_overlaps_i, slot_i)` via the pure `&self` `compute_read_assignment`, installed in
//! one shot via `set_all_read_assignments`); and (B) `AddFragmentAlignmentInfo`, whose only
//! mutation is each read's own `fragment_assignments[i]`. Every genuinely order-dependent step
//! (`CoalesceReadAssignments`, EM/quantification, `ComputeVariant`) still runs strictly
//! sequentially in a fixed, thread-count-independent order.
//! `--barcode`/`--relaxIntronAlign`/`--alleleDigitUnits`/`--alleleDelimiter` are not exposed by
//! [`AnalyzeArgs`] -- see that struct's doc comment.
//!
//! # Reference loading mirrors `Genotyper::InitRefSet(char*, selectedAlleles)` (`Genotyper.hpp:732-757`)
//!
//! [`load_selected_reference`] is [`crate::stages::genotype::load_reference`]'s exact dedup/
//! effective-length logic, but restricted to allele names present in `selectedAlleles`
//! (`Analyzer.cpp:343-352`'s `-a` allele-list file, parsed by [`parse_selected_alleles`] --
//! ONE whitespace-delimited token per line, matching C++'s `sscanf(buffer, "%s", alleleName)`,
//! which silently ignores the trailing `genotypeQuality` column `_allele.tsv` also carries).
//!
//! # Abundances ARE re-quantified, restricted to the selected alleles
//!
//! `Analyzer.cpp:main` calls `genotyper.QuantifyAlleleEquivalentClass()` (`Analyzer.cpp:619`,
//! right after `FinalizeReadAssignments`/before `AddFragmentAlignmentInfo`) -- the SAME EM
//! abundance estimator `genotyper`'s own default (no `-a <abundance file>`) path runs, just over
//! the allele-RESTRICTED reference this driver's `InitRefSet(refFile, selectedAlleles)`-equivalent
//! loading built. `Genotyper::InitAlleleAbundance` (the `genotyper -a <file>` precomputed-
//! abundance mode) has exactly one caller (`Genotyper.cpp:641`) and is never used by `analyzer`
//! -- confirmed by `grep -n "InitAlleleAbundance" vendor/t1k/*.cpp`. Unlike `genotyper`'s own
//! driver, this one does NOT call `RemoveLowLikelihoodAlleleInEquivalentClass` or
//! `SelectAllelesForGenes` afterward (those are genotype-SELECTION steps; by the time `analyzer`
//! runs, the selection already happened and its result is exactly what `-a` fed in as
//! `selectedAlleles`) -- `VariantCaller::SetSeqAbundance` only reads `AlleleInfo::abundance`/
//! `geneIdx`, both of which `QuantifyAlleleEquivalentClass`/`init_allele_info` alone already
//! populate.
use crate::cli::AnalyzeArgs;
use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Arc;
use unum_core::fastq::FastqReader;
use unum_core::genotyper::{self, AlleleRef, ExtendedOverlap, Genotyper, ReadAssignment};
use unum_core::ref_kmer_filter::RefKmerFilter;
use unum_core::variant_caller::{self, FragmentOverlap, Overlap, VariantCaller};

/// The k-mer length the genotyper's reference index/alignment is built at
/// (`Genotyper genotyper(11)`, `Analyzer.cpp:231`) -- same value
/// `crate::stages::genotype::GENOTYPER_KMER_LENGTH` uses (both `Genotyper.cpp`
/// and `Analyzer.cpp` construct `Genotyper` at k=11).
const GENOTYPER_KMER_LENGTH: usize = 11;

/// `KmerCount`'s own default-constructor k-mer length -- see
/// `crate::stages::genotype::GENE_SIMILARITY_KMER_LENGTH`'s identical doc
/// comment.
const GENE_SIMILARITY_KMER_LENGTH: usize = 31;

/// A loaded, deduplicated, ALLELE-RESTRICTED reference: everything
/// [`unum_core::genotyper::Genotyper::init_allele_info`] and this driver's
/// read-reassignment loop need per (selected) allele.
struct LoadedRef {
    names: Vec<String>,
    consensus: Vec<Arc<[u8]>>,
    weight: Vec<i32>,
    effective_len: Vec<i32>,
    allele_refs: Vec<AlleleRef>,
}

/// Ported from `Analyzer.cpp:343-351`: reads the `-a` selected-allele-list
/// file (`{prefix}_allele.tsv`, `Genotyper::OutputRepresentativeAlleles`'s
/// output) into a set of allele NAMES. Only the FIRST whitespace-delimited
/// token of each line is kept (`sscanf(buffer, "%s", alleleName)`) -- the
/// trailing `genotypeQuality` column `_allele.tsv` also carries is parsed
/// into `buffer` by C++'s `fgets` but never extracted past that first
/// `sscanf` token, so it is correctly ignored here too.
///
/// # Errors
///
/// Returns an error if `path` cannot be read.
fn parse_selected_alleles(path: &Path) -> Result<HashSet<String>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading selected-allele file {}", path.display()))?;
    Ok(text.lines().filter_map(|line| line.split_whitespace().next()).map(String::from).collect())
}

/// Ported from `Genotyper::InitRefSet(char *filename, const
/// std::map<std::string, int> &selectedAlleles)` (`Genotyper.hpp:732-757`):
/// identical to [`crate::stages::genotype::load_reference`]'s dedup/
/// `effectiveLen` logic, except records whose id is not in `selected` are
/// skipped BEFORE the dedup step (`Genotyper.hpp:742-743`'s `if
/// (selectedAlleles.find(fa.id) == selectedAlleles.end()) continue ;`).
///
/// This port writes its own minimal FASTA-with-comment reader (matching
/// [`crate::stages::genotype::load_reference`]'s own rationale: exon-interval
/// comments are load-bearing for [`AlleleRef::new`], and [`FastqReader`]
/// discards header comments).
///
/// # Errors
///
/// Returns an error if `path` cannot be read, is not valid UTF-8/FASTA, or no
/// allele in `path` matches `selected`.
fn load_selected_reference(path: &Path, selected: &HashSet<String>) -> Result<LoadedRef> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading reference FASTA {}", path.display()))?;

    let mut seq_to_idx: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut names = Vec::new();
    let mut consensus: Vec<Arc<[u8]>> = Vec::new();
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
                 consensus: &mut Vec<Arc<[u8]>>,
                 weight: &mut Vec<i32>,
                 comments: &mut Vec<Option<String>>| {
        let Some(id) = id else { return };
        if seq.is_empty() || !selected.contains(&id) {
            return;
        }
        if let Some(&idx) = seq_to_idx.get(&seq) {
            weight[idx] += 1;
        } else {
            let idx = names.len();
            seq_to_idx.insert(seq.clone(), idx);
            names.push(id);
            consensus.push(Arc::from(seq));
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

    ensure!(
        !names.is_empty(),
        "reference FASTA {} contains no sequences matching the selected-allele list",
        path.display()
    );

    let effective_len: Vec<i32> = consensus.iter().map(|s| compute_effective_len(s)).collect();
    let allele_refs: Vec<AlleleRef> = consensus
        .iter()
        .zip(&comments)
        .map(|(seq, comment)| AlleleRef::new(seq.clone(), comment.as_deref()))
        .collect();

    Ok(LoadedRef { names, consensus, weight, effective_len, allele_refs })
}

/// Ported from `SeqSet::ComputeEffectiveLen` (`SeqSet.hpp:747-758`) -- see
/// `crate::stages::genotype::compute_effective_len`'s identical doc comment
/// (duplicated here rather than shared, matching this module's "additive,
/// not shared" convention with `genotype.rs`, since that function is
/// private).
fn compute_effective_len(seq: &[u8]) -> i32 {
    let mut ret = 0i32;
    for (i, &b) in seq.iter().enumerate() {
        if b != b'N' || (i > 0 && seq[i - 1] != b'N') {
            ret += 1;
        }
    }
    ret
}

/// Builds a [`RefKmerFilter`] over `names`/`consensus` (the deduplicated,
/// allele-restricted reference from [`load_selected_reference`]) at
/// [`GENOTYPER_KMER_LENGTH`]. See
/// `crate::stages::genotype::build_ref_kmer_filter` (a near-identical helper)
/// for how the `Arc<[u8]>` consensus buffers are shared with the k-mer index.
fn build_ref_kmer_filter(names: &[String], consensus: &[Arc<[u8]>]) -> RefKmerFilter {
    // Share the caller's `Arc<[u8]>` consensus buffers (see
    // `crate::stages::genotype::build_ref_kmer_filter`): no scratch FASTA, no
    // duplicate reference copy, byte-identical k-mer index.
    RefKmerFilter::from_consensus(names.to_vec(), consensus.to_vec(), GENOTYPER_KMER_LENGTH)
}

/// One aligned read (mate) loaded from an `_aligned_*.fa` FASTA, mirroring
/// the fields of `_genotypeRead` (`Analyzer.cpp:51-73`) this driver actually
/// needs (barcode/UMI/qual are dropped -- see this module's doc comment).
struct AnalyzeRead {
    seq: Vec<u8>,
    has_n: bool,
}

fn read_all(reader: &mut FastqReader) -> Result<Vec<AnalyzeRead>> {
    let mut reads = Vec::new();
    while let Some(record) = reader.next_record()? {
        let has_n = record.seq.contains(&b'N');
        reads.push(AnalyzeRead { seq: record.seq, has_n });
    }
    Ok(reads)
}

/// Converts one [`ExtendedOverlap`] into a [`variant_caller::Overlap`] with
/// `align: None` (populated later by [`variant_caller::add_fragment_alignment_info`]),
/// mirroring the `struct _overlap` fields `_fragmentOverlap::overlap1`/
/// `overlap2` carry forward from `ReadAssignmentToFragmentAssignment` before
/// `AddOverlapAlignmentInfo` fills in `align`.
fn to_overlap(o: &ExtendedOverlap) -> Overlap {
    Overlap {
        seq_idx: i32::try_from(o.seq_idx).unwrap(),
        read_start: o.read_start,
        read_end: o.read_end,
        seq_start: o.seq_start,
        seq_end: o.seq_end,
        strand: o.strand,
        match_cnt: o.match_cnt,
        similarity: o.similarity,
        align: None,
    }
}

/// Runs the `analyze` subcommand for `args`.
///
/// # Errors
///
/// Returns an error if the reference/selected-allele/aligned-read files
/// cannot be opened/parsed, if neither `-u` nor both `-1`/`-2` are given, or
/// if the output file cannot be created.
// Long by construction: mirrors Analyzer.cpp:main's single linear driver
// (selected-allele parse -> restricted reference load -> read load ->
// alignment -> fragment assembly -> AddFragmentAlignmentInfo -> ComputeVariant
// -> OutputAlleleVCF) -- same rationale as crate::stages::genotype::run's own
// too_many_lines allow.
#[allow(clippy::too_many_lines)]
pub fn run(args: &AnalyzeArgs) -> Result<()> {
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
        bail!("must specify either -u (single-end) or -1/-2 (paired) aligned-read input");
    };
    let has_mate = mate2_path.is_some();

    // --- Selected-allele parse + restricted reference loading (Analyzer.cpp:343-353, Genotyper.hpp:732-757) ---
    let selected = parse_selected_alleles(Path::new(&args.allele_file))?;
    let loaded = load_selected_reference(Path::new(&args.ref_seq_fasta), &selected)?;
    let allele_cnt = loaded.names.len();

    let mut filter = build_ref_kmer_filter(&loaded.names, &loaded.consensus);
    filter.set_ref_seq_similarity(args.similarity);

    let mut genotyper = Genotyper::new();
    let mut effective_len = loaded.effective_len.clone();
    genotyper.init_allele_info(
        &loaded.names,
        &loaded.consensus,
        &loaded.weight,
        &mut effective_len,
        GENE_SIMILARITY_KMER_LENGTH,
    );

    // --- Read aligned FASTA(s) (Analyzer.cpp:368-448) ---
    let mut reader1 = FastqReader::open(Path::new(mate1_path))
        .with_context(|| format!("opening aligned read file {mate1_path}"))?;
    let reads1 = read_all(&mut reader1)?;
    let reads2 = if let Some(mate2_path) = mate2_path {
        let mut reader2 = FastqReader::open(Path::new(mate2_path))
            .with_context(|| format!("opening aligned read file {mate2_path}"))?;
        read_all(&mut reader2)?
    } else {
        Vec::new()
    };
    if has_mate {
        ensure!(
            reads1.len() == reads2.len(),
            "mate-1 ({}) and mate-2 ({}) aligned read counts differ",
            reads1.len(),
            reads2.len()
        );
    }
    let read_cnt = reads1.len();
    let max_read_length =
        reads1.iter().chain(reads2.iter()).map(|r| r.seq.len()).max().unwrap_or(0);
    genotyper.set_read_length(i32::try_from(max_read_length).unwrap_or(0));

    // `-t`/`--threads`: honored byte-identically across the three per-read
    // loops below (Loop 0's `assign_reads_parallel`, Loop A's slot-indexed
    // fragment assembly, and Loop B's `AddFragmentAlignmentInfo`) -- see each
    // loop's comment and this module's doc comment for the per-loop byte-
    // identity argument. Mirror `crate::stages::genotype::run`'s clamp.
    let threads = usize::try_from(args.threads).unwrap_or(usize::MAX).max(1);

    // --- Read-end alignment, reusing identical sequences (Analyzer.cpp:455-511) ---
    // Not `mut`: `assign_read` takes `&[AlleleRef]` (coverage marking, when it
    // runs, goes through interior-mutable `AtomicPosWeight`); the analyzer path
    // passes `weight=0` so no coverage is marked at all here.
    let allele_refs = loaded.allele_refs;
    let mut all_seqs: Vec<&[u8]> = reads1.iter().map(|r| r.seq.as_slice()).collect();
    all_seqs.extend(reads2.iter().map(|r| r.seq.as_slice()));
    let mut sorted_seqs = all_seqs.clone();
    sorted_seqs.sort_unstable();
    sorted_seqs.dedup();

    // Loop 0 -- dedup `assign_read` (Analyzer.cpp:476). Delegates to the SAME
    // `assign_reads_parallel` helper `stages::genotype` uses, but with the
    // weight closure returning 0 (matching the analyzer's `assign_read(..., 0)`,
    // NOT an occurrence count): unlike the genotyper, the analyzer does NOT
    // accumulate per-base `pos_weight` coverage here (its `missing_coverage` is
    // never read on the analyzer path -- `em_update`/`set_allele_abundance` use
    // only ec length, and the missing-coverage consumers are the genotype-
    // selection steps analyzer never runs). With `weight == 0` the interior-
    // mutable `AtomicPosWeight` coverage marking is a no-op, so `allele_refs`'
    // end state -- and thus every downstream result -- is byte-identical at any
    // thread count and to the serial `-t 1` loop (see `assign_reads_parallel`'s
    // doc comment for the order-invariance argument).
    let extended_by_seq = genotyper::assign_reads_parallel(
        &filter,
        &sorted_seqs,
        &allele_refs,
        args.similarity,
        |_seq| 0,
        threads,
    );
    let mut overlaps_by_seq: HashMap<&[u8], Option<Vec<ExtendedOverlap>>> = HashMap::new();
    for (&seq, extended) in sorted_seqs.iter().zip(extended_by_seq) {
        overlaps_by_seq.insert(seq, extended);
    }

    // --- Fragment assembly + SetReadAssignments (Analyzer.cpp:452,528-565) ---
    // Fail loudly on an oversized input rather than masking the conversion to 0
    // (which would silently treat billions of reads as an empty run).
    let read_cnt_i32 = i32::try_from(read_cnt).with_context(|| {
        format!("read count {read_cnt} exceeds the supported maximum (i32::MAX)")
    })?;
    genotyper.init_read_assignments(read_cnt_i32, args.max_assign_cnt);
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

    // Loop A -- fragment assembly + read-assignment slot computation
    // (Analyzer.cpp:528-565). Parallelized across reads and byte-identical to
    // `-t 1` because it is SLOT-INDEXED (the exact pattern
    // `stages::genotype::run` uses): each read `i` produces a pure
    // `(fragment_overlaps_i, slot_i)` tuple from read `i`'s inputs plus
    // immutable shared state -- `overlaps_by_seq`, the read vectors, `&genotyper`
    // (read only via the `&self` `compute_read_assignment`), and the `Fn + Sync`
    // lookup closures (`consensus_len_of`/`separator_in_range`/
    // `separator_lookup_for_set_read_assignments`, each borrowing `allele_refs`
    // immutably) -- with NO cross-read order dependence. `.unzip()` on the
    // order-preserving `into_par_iter()` yields `fragment_assignments` and
    // `slots` in slot order regardless of thread count; the slots are then
    // installed in one shot via `set_all_read_assignments`.
    //
    // `fragment_assignments[i]` mirrors `Analyzer.cpp`'s own
    // `fragmentAssignments[i]` -- kept (not discarded like `genotype.rs`
    // does) since `AddFragmentAlignmentInfo`/`ComputeVariant` need it after
    // the read loop (`Analyzer.cpp:464,556,684`). Runs inside a pool sized to
    // `-t` so `-t 1` is strictly single-threaded (and still byte-identical).
    let fragment_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("building rayon thread pool for parallel fragment assembly");
    let (mut fragment_assignments, slots): (Vec<Vec<FragmentOverlap>>, Vec<Vec<ReadAssignment>>) =
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
                    let assembled = genotyper::read_assignment_to_fragment_assignment_with_overlaps(
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

                    let fragment_overlaps: Vec<FragmentOverlap> = assembled
                        .iter()
                        .map(|(fo, o1, o2)| FragmentOverlap {
                            seq_idx: fo.seq_idx,
                            has_mate_pair: fo.has_mate_pair,
                            o1_from_r2: fo.o1_from_r2,
                            overlap1: to_overlap(o1),
                            overlap2: o2.map_or_else(Overlap::none, |o2| to_overlap(&o2)),
                        })
                        .collect();

                    let plain_assignment: Vec<genotyper::FragmentOverlap> =
                        assembled.iter().map(|(fo, _, _)| *fo).collect();
                    let slot = genotyper.compute_read_assignment(
                        &plain_assignment,
                        ref_seq_similarity,
                        separator_lookup_for_set_read_assignments,
                    );

                    (fragment_overlaps, slot)
                })
                .unzip()
        });
    genotyper.set_all_read_assignments(slots);
    // Skip coalescing entirely when there are no reads: `coalesce_read_assignments`
    // would otherwise be handed `end = 0` and evaluate `total_read_cnt - 1` on an
    // empty assignment set.
    if read_cnt > 0 {
        genotyper.coalesce_read_assignments(0, read_cnt_i32 - 1);
    }

    // --- FinalizeReadAssignments + quantification (Analyzer.cpp:613,619) ---
    // No `RemoveLowLikelihoodAlleleInEquivalentClass`/`SelectAllelesForGenes`
    // call here -- those are genotype-SELECTION steps `analyzer` does not run
    // (see this module's doc comment).
    // Per-allele independent -> embarrassingly parallel, byte-identical to the
    // serial map (each allele reads only its own `pos_weight`).
    let missing_coverage: Vec<i32> = (0..allele_cnt)
        .into_par_iter()
        .map(|i| genotyper::get_seq_missing_base_coverage(&allele_refs[i], 0.01))
        .collect();
    genotyper.finalize_read_assignments(&missing_coverage);
    genotyper.quantify_allele_equivalent_class(&effective_len, &loaded.weight);

    // --- Base-level variant identification (Analyzer.cpp:671-686) ---
    let mut variant_caller = VariantCaller::new(&allele_refs);
    variant_caller.set_seq_abundance(&genotyper, allele_cnt);
    variant_caller.set_max_var_group_to_resolve(args.var_max_group);

    let read1_seqs: Vec<Vec<u8>> = reads1.iter().map(|r| r.seq.clone()).collect();
    let read2_seqs: Vec<Vec<u8>> =
        if has_mate { reads2.iter().map(|r| r.seq.clone()).collect() } else { Vec::new() };
    let allele_consensus: Vec<Vec<u8>> = allele_refs.iter().map(|a| a.consensus.to_vec()).collect();

    // Loop B -- AddFragmentAlignmentInfo (Analyzer.cpp:623-634): populate
    // `align` on every fragment overlap BEFORE ComputeVariant runs -- required,
    // not optional (see `variant_caller::add_fragment_alignment_info`'s doc
    // comment). Parallelized across reads: each call mutates ONLY its own
    // `fragment_assignments[i]` (zero shared mutable state --
    // `variant_caller.rs:469`), reading only the shared immutable
    // `read1_seqs`/`read2_seqs`/`allele_consensus` borrows, so
    // `par_iter_mut().enumerate()` is byte-identical to the serial `-t 1` loop
    // at any thread count. This is the only per-read loop with real DP cost
    // (`global_alignment` per overlap), so it is where any speedup comes from.
    // Runs inside a `-t`-sized pool so `-t 1` stays strictly single-threaded.
    let align_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("building rayon thread pool for parallel fragment alignment");
    align_pool.install(|| {
        fragment_assignments.par_iter_mut().enumerate().for_each(|(i, frag)| {
            let r2 = if has_mate { Some(read2_seqs[i].as_slice()) } else { None };
            variant_caller::add_fragment_alignment_info(
                &read1_seqs[i],
                r2,
                frag,
                &allele_consensus,
            );
        });
    });

    variant_caller.compute_variant(
        &read1_seqs,
        &read2_seqs,
        &fragment_assignments,
        &allele_consensus,
    );

    // --- Output (Analyzer.cpp:672,686) ---
    let vcf_path = format!("{}_allele.vcf", args.prefix);
    let vcf_text = variant_caller.output_allele_vcf(&loaded.names, |seq_idx, pos| {
        variant_caller::exonic_position(&allele_refs[seq_idx].exons, pos)
    });
    std::fs::write(&vcf_path, vcf_text).with_context(|| format!("writing {vcf_path}"))?;

    eprintln!(
        "post-analyzed {read_cnt} read fragments across {allele_cnt} selected alleles; wrote \
         {vcf_path}",
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_selected_alleles_keeps_only_the_first_token_per_line() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            tmp.path(),
            "KIR2DL1*0010101 60\n\
             KIR2DL3*0010101 60\n",
        )
        .unwrap();
        let selected = parse_selected_alleles(tmp.path()).unwrap();
        assert_eq!(selected.len(), 2);
        assert!(selected.contains("KIR2DL1*0010101"));
        assert!(selected.contains("KIR2DL3*0010101"));
    }

    #[test]
    fn parse_selected_alleles_ignores_blank_lines() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), "A*01:01 60\n\n\nB*07:02 60\n").unwrap();
        let selected = parse_selected_alleles(tmp.path()).unwrap();
        assert_eq!(selected.len(), 2);
    }

    #[test]
    fn load_selected_reference_restricts_to_selected_names() {
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

        let mut selected = HashSet::new();
        selected.insert("A*01:01".to_string());

        let loaded = load_selected_reference(tmp.path(), &selected).unwrap();
        assert_eq!(loaded.names, vec!["A*01:01".to_string()]);
        assert_eq!(loaded.weight, vec![1]);
        assert_eq!(loaded.consensus.len(), 1);
        assert_eq!(loaded.consensus[0].as_ref(), b"ACGT");
    }

    #[test]
    fn load_selected_reference_dedups_identical_selected_sequences() {
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

        let mut selected = HashSet::new();
        selected.insert("A*01:01".to_string());
        selected.insert("A*01:02".to_string());
        selected.insert("B*07:02".to_string());

        let loaded = load_selected_reference(tmp.path(), &selected).unwrap();
        assert_eq!(loaded.names.len(), 2, "identical ACGT sequences must collapse to one allele");
        let a_idx = loaded.names.iter().position(|n| n == "A*01:01").unwrap();
        assert_eq!(loaded.weight[a_idx], 2, "duplicate sequence increments weight");
    }

    #[test]
    fn load_selected_reference_errors_when_nothing_matches() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), ">A*01:01 1 3\nACGT\n").unwrap();
        let selected = HashSet::new();
        assert!(load_selected_reference(tmp.path(), &selected).is_err());
    }
}
