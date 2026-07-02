//! Thin CLI wrapper around `fg_t1k_core::genotyper`/`fg_t1k_core::variant_caller` (the Rust port
//! of `analyzer`, `vendor/t1k/Analyzer.cpp`). All post-analysis logic (restricted reference
//! loading, read reassignment, pileup/`ComputeVariant`) is either in those two modules or -- for
//! the read-I/O/reassignment DRIVER itself, which `Analyzer.cpp:main` owns directly rather than
//! delegating to a `Genotyper`/`VariantCaller` method -- right here, matching
//! `crate::stages::genotype`'s same split of responsibility.
//!
//! # Scope: single-threaded, paired/single-end aligned-read FASTA, no barcode
//!
//! This port only reproduces `Analyzer.cpp:main`'s `threadCnt <= 1` code path
//! (`Analyzer.cpp:467-484,528-565,623-634`) -- the multi-threaded path is a batching/parallelism
//! detail over the exact same per-read logic, not a different algorithm, so single-threaded
//! output is what an end-to-end differential against the real oracle (also run at `-t 1`) must
//! match. `--barcode`/`--relaxIntronAlign`/`--alleleDigitUnits`/`--alleleDelimiter` are not
//! exposed by [`AnalyzeArgs`] -- see that struct's doc comment.
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
use fg_t1k_core::fastq::FastqReader;
use fg_t1k_core::genotyper::{self, AlleleRef, ExtendedOverlap, Genotyper};
use fg_t1k_core::ref_kmer_filter::RefKmerFilter;
use fg_t1k_core::variant_caller::{self, FragmentOverlap, Overlap, VariantCaller};
use std::collections::{HashMap, HashSet};
use std::io::Write as _;
use std::path::{Path, PathBuf};

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
/// [`fg_t1k_core::genotyper::Genotyper::init_allele_info`] and this driver's
/// read-reassignment loop need per (selected) allele.
struct LoadedRef {
    names: Vec<String>,
    consensus: Vec<Vec<u8>>,
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
        if seq.is_empty() || !selected.contains(&id) {
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
/// `crate::stages::genotype::build_ref_kmer_filter`'s identical doc comment
/// for why a scratch file is used.
fn build_ref_kmer_filter(names: &[String], consensus: &[Vec<u8>]) -> Result<RefKmerFilter> {
    let scratch_path: PathBuf = std::env::temp_dir().join(format!(
        "fg-t1k-analyze-ref-{}-{}.fa",
        std::process::id(),
        SCRATCH_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    ));

    {
        let mut f = std::fs::File::create(&scratch_path).with_context(|| {
            format!("creating scratch reference FASTA {}", scratch_path.display())
        })?;
        for (name, seq) in names.iter().zip(consensus) {
            writeln!(f, ">{name}")?;
            f.write_all(seq)?;
            writeln!(f)?;
        }
    }

    let result = RefKmerFilter::from_reference_fasta(&scratch_path, GENOTYPER_KMER_LENGTH)
        .with_context(|| format!("building k-mer index over {} selected alleles", names.len()));
    let _ = std::fs::remove_file(&scratch_path);
    result
}

static SCRATCH_COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

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

    let mut filter = build_ref_kmer_filter(&loaded.names, &loaded.consensus)?;
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

    // --- Read-end alignment, reusing identical sequences (Analyzer.cpp:455-511) ---
    let mut allele_refs = loaded.allele_refs;
    let mut all_seqs: Vec<&[u8]> = reads1.iter().map(|r| r.seq.as_slice()).collect();
    all_seqs.extend(reads2.iter().map(|r| r.seq.as_slice()));
    let mut sorted_seqs = all_seqs.clone();
    sorted_seqs.sort_unstable();
    sorted_seqs.dedup();

    let mut overlaps_by_seq: HashMap<&[u8], Option<Vec<ExtendedOverlap>>> = HashMap::new();
    for &seq in &sorted_seqs {
        let raw_overlaps = filter
            .get_overlaps_from_read(seq, &mut fg_t1k_core::ref_kmer_filter::Scratch::default())
            .unwrap_or_default();
        // Analyzer.cpp:476 calls AssignRead(..., weight=0): unlike the genotyper,
        // the analyzer does NOT accumulate per-base pos_weight coverage here (its
        // `missing_coverage` is never read on the analyzer path -- em_update /
        // set_allele_abundance use only ec length, and the missing-coverage
        // consumers are the genotype-selection steps analyzer never runs). Pass 0
        // to match the oracle exactly rather than the genotyper's occurrence count.
        let extended =
            genotyper::assign_read(seq, &raw_overlaps, &mut allele_refs, args.similarity, 0);
        overlaps_by_seq.insert(seq, extended);
    }

    // --- Fragment assembly + SetReadAssignments (Analyzer.cpp:452,528-565) ---
    genotyper.init_read_assignments(i32::try_from(read_cnt).unwrap_or(0), args.max_assign_cnt);
    let consensus_len_of = |idx: u32| {
        i32::try_from(allele_refs[usize::try_from(idx).unwrap_or(0)].consensus.len()).unwrap_or(0)
    };
    let hit_len_required = filter.hit_len_required();
    let ref_seq_similarity = filter.ref_seq_similarity();

    // `fragment_assignments[i]` mirrors `Analyzer.cpp`'s own
    // `fragmentAssignments[i]` -- kept (not discarded like `genotype.rs`
    // does) since `AddFragmentAlignmentInfo`/`ComputeVariant` need it after
    // the read loop (`Analyzer.cpp:464,556,684`).
    let mut fragment_assignments: Vec<Vec<FragmentOverlap>> = Vec::with_capacity(read_cnt);
    for i in 0..read_cnt {
        let overlaps1 = overlaps_by_seq.get(reads1[i].seq.as_slice()).and_then(Option::as_ref);
        let overlaps2 = if has_mate {
            overlaps_by_seq.get(reads2[i].seq.as_slice()).and_then(Option::as_ref)
        } else {
            None
        };
        let has_n = reads1[i].has_n || (has_mate && reads2[i].has_n);

        let empty: Vec<ExtendedOverlap> = Vec::new();
        let assembled = genotyper::read_assignment_to_fragment_assignment_with_overlaps(
            overlaps1.map_or(empty.as_slice(), Vec::as_slice),
            if has_mate { Some(overlaps2.map_or(empty.as_slice(), Vec::as_slice)) } else { None },
            has_n,
            hit_len_required,
            consensus_len_of,
            |_, _, _| false, // no interior-N separators in this port's references (see genotyper.rs doc).
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
        genotyper.set_read_assignments(i, &plain_assignment, ref_seq_similarity, |_, _, _| false);

        fragment_assignments.push(fragment_overlaps);
    }
    genotyper.coalesce_read_assignments(0, i32::try_from(read_cnt.saturating_sub(1)).unwrap_or(0));

    // --- FinalizeReadAssignments + quantification (Analyzer.cpp:613,619) ---
    // No `RemoveLowLikelihoodAlleleInEquivalentClass`/`SelectAllelesForGenes`
    // call here -- those are genotype-SELECTION steps `analyzer` does not run
    // (see this module's doc comment).
    let missing_coverage: Vec<i32> = (0..allele_cnt)
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
    let allele_consensus: Vec<Vec<u8>> = allele_refs.iter().map(|a| a.consensus.clone()).collect();

    // AddFragmentAlignmentInfo (Analyzer.cpp:623-634): populate `align` on
    // every fragment overlap BEFORE ComputeVariant runs -- required, not
    // optional (see `variant_caller::add_fragment_alignment_info`'s doc
    // comment).
    for i in 0..read_cnt {
        let r2 = if has_mate { Some(read2_seqs[i].as_slice()) } else { None };
        variant_caller::add_fragment_alignment_info(
            &read1_seqs[i],
            r2,
            &mut fragment_assignments[i],
            &allele_consensus,
        );
    }

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
        assert_eq!(loaded.consensus, vec![b"ACGT".to_vec()]);
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
