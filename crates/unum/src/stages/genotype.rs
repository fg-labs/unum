//! Thin CLI wrapper around `unum_core::genotyper` (the Rust port of `genotyper`,
//! `Genotyper.cpp`/`Genotyper.hpp`/(the `ReadAssignmentToFragmentAssignment` slice
//! of) `SeqSet.hpp`). All genotyping logic (reference dedup/loading, the read-processing loop,
//! `SelectAllelesForGenes`, output formatting) is either in `unum_core::genotyper` or -- for
//! the reference-loading/read-processing DRIVER itself, which `Genotyper.cpp:main` owns directly
//! in the C++ rather than delegating to a `Genotyper` method -- right here, matching that same
//! split of responsibility.
//!
//! # Scope: single-threaded, paired/single-end FASTQ, no barcode/`-a`/whitelist
//!
//! This port only reproduces `Genotyper.cpp:main`'s `threadCnt <= 1` code path (`Genotyper.cpp:
//! 463-480,531-574`) -- the multi-threaded path is a batching/parallelism detail over the exact
//! same per-read logic, not a different algorithm, so single-threaded output is what an
//! end-to-end differential against the real oracle (also run at `-t 1`) must match.
//! `--barcode`/`-a` (a precomputed abundance file, skipping `QuantifyAlleleEquivalentClass`
//! entirely)/`--alleleWhitelist`/`--outputReadAssignment` are not exposed by [`GenotypeArgs`] --
//! see that struct's doc comment.
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
use std::collections::HashMap;
use std::io::Write as _;
use std::path::Path;
use unum_core::fastq::FastqReader;
use unum_core::genotyper::{self, AlleleRef, ExtendedOverlap, Genotyper};
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
/// is kept (matching `_genotypeRead::id`) even though this port's output
/// (`_genotype.tsv`/`_allele.tsv`) never reads it -- only `_aligned*.fa`
/// (`Genotyper.cpp:680-707`, which this port does not yet write) would need
/// it; kept for fidelity/future use rather than dropped.
#[allow(dead_code)]
struct GenotypeRead {
    id: String,
    seq: Vec<u8>,
    has_n: bool,
}

fn read_all(reader: &mut FastqReader) -> Result<Vec<GenotypeRead>> {
    let mut reads = Vec::new();
    while let Some(record) = reader.next_record()? {
        let has_n = record.seq.contains(&b'N');
        reads.push(GenotypeRead { id: record.id, seq: record.seq, has_n });
    }
    Ok(reads)
}

/// Runs the `genotype` subcommand for `args`.
///
/// # Errors
///
/// Returns an error if the reference/candidate-read files cannot be
/// opened/parsed, if neither `-u` nor both `-1`/`-2` are given, or if output
/// files cannot be created.
// Long by construction: mirrors Genotyper.cpp:main's single linear driver
// (reference load -> read load -> alignment -> fragment assembly ->
// quantify/select -> output) -- splitting this into several tiny functions
// would scatter closely-related sequential setup across the file for no
// readability benefit (same rationale as the vendored source's own
// single-function `main`).
#[allow(clippy::too_many_lines)]
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

    // --- Reference loading (Genotyper.hpp:706-727) ---
    let loaded = load_reference(Path::new(&args.ref_seq_fasta))?;
    let allele_cnt = loaded.names.len();

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

    // --- Read candidate FASTQ(s) (Genotyper.cpp:362-443) ---
    let mut reader1 = FastqReader::open(Path::new(mate1_path))
        .with_context(|| format!("opening candidate read file {mate1_path}"))?;
    let reads1 = read_all(&mut reader1)?;
    let reads2 = if let Some(mate2_path) = mate2_path {
        let mut reader2 = FastqReader::open(Path::new(mate2_path))
            .with_context(|| format!("opening candidate read file {mate2_path}"))?;
        read_all(&mut reader2)?
    } else {
        Vec::new()
    };
    if has_mate {
        ensure!(
            reads1.len() == reads2.len(),
            "mate-1 ({}) and mate-2 ({}) candidate read counts differ",
            reads1.len(),
            reads2.len()
        );
    }
    let read_cnt = reads1.len();
    let max_read_length =
        reads1.iter().chain(reads2.iter()).map(|r| r.seq.len()).max().unwrap_or(0);
    genotyper.set_read_length(i32::try_from(max_read_length).unwrap_or(0));

    // --- Read-end alignment, reusing identical sequences (Genotyper.cpp:445-480) ---
    // `overlaps_by_seq` mirrors `AssignReads_Thread`/the single-threaded loop
    // at `Genotyper.cpp:463-480`: sort ALL read ends (mate 1 + mate 2) by
    // sequence, and call AssignRead once per distinct sequence, reusing the
    // result for every read end sharing that exact sequence. `allele_refs`
    // must be `&mut` (AssignRead marks base coverage in place), so this is
    // a plain sequential loop -- see this module's "single-threaded" scope
    // note.
    let mut allele_refs = loaded.allele_refs;
    let mut all_seqs: Vec<&[u8]> = reads1.iter().map(|r| r.seq.as_slice()).collect();
    all_seqs.extend(reads2.iter().map(|r| r.seq.as_slice()));
    let mut sorted_seqs = all_seqs.clone();
    sorted_seqs.sort_unstable();
    sorted_seqs.dedup();

    let mut overlaps_by_seq: HashMap<&[u8], Option<Vec<ExtendedOverlap>>> = HashMap::new();
    let mut counted: HashMap<&[u8], i32> = HashMap::new();
    for &seq in &all_seqs {
        *counted.entry(seq).or_insert(0) += 1;
    }
    for &seq in &sorted_seqs {
        let weight = counted[seq];
        let raw_overlaps = filter
            .get_overlaps_from_read(seq, &mut unum_core::ref_kmer_filter::Scratch::default())
            .unwrap_or_default();
        let extended =
            genotyper::assign_read(seq, &raw_overlaps, &mut allele_refs, args.similarity, weight);
        overlaps_by_seq.insert(seq, extended);
    }

    // --- Fragment assembly + SetReadAssignments (Genotyper.cpp:531-574) ---
    genotyper.init_read_assignments(i32::try_from(read_cnt).unwrap_or(0), args.max_assign_cnt);
    let consensus_len_of = |idx: u32| {
        i32::try_from(allele_refs[usize::try_from(idx).unwrap_or(0)].consensus.len()).unwrap_or(0)
    };
    let hit_len_required = filter.hit_len_required();
    let ref_seq_similarity = filter.ref_seq_similarity();

    for i in 0..read_cnt {
        let overlaps1 = overlaps_by_seq.get(reads1[i].seq.as_slice()).and_then(Option::as_ref);
        let overlaps2 = if has_mate {
            overlaps_by_seq.get(reads2[i].seq.as_slice()).and_then(Option::as_ref)
        } else {
            None
        };
        let has_n = reads1[i].has_n || (has_mate && reads2[i].has_n);

        let empty: Vec<ExtendedOverlap> = Vec::new();
        let assignment = genotyper::read_assignment_to_fragment_assignment(
            overlaps1.map_or(empty.as_slice(), Vec::as_slice),
            if has_mate { Some(overlaps2.map_or(empty.as_slice(), Vec::as_slice)) } else { None },
            has_n,
            hit_len_required,
            consensus_len_of,
            |_, _, _| false, // no interior-N separators in this port's references (see genotyper.rs doc).
        );
        genotyper.set_read_assignments(i, &assignment, ref_seq_similarity, |_, _, _| false);
    }

    genotyper.coalesce_read_assignments(0, i32::try_from(read_cnt.saturating_sub(1)).unwrap_or(0));

    // --- FinalizeReadAssignments + quantification + selection (Genotyper.cpp:634-650) ---
    let missing_coverage: Vec<i32> = (0..allele_cnt)
        .map(|i| genotyper::get_seq_missing_base_coverage(&allele_refs[i], 0.01))
        .collect();
    genotyper.finalize_read_assignments(&missing_coverage);

    genotyper.quantify_allele_equivalent_class(&effective_len, &loaded.weight);
    genotyper.remove_low_likelihood_allele_in_equivalent_class(|idx| effective_len[idx]);
    genotyper.select_alleles_for_genes(|idx| loaded.weight[idx]);

    // --- Output (Genotyper.cpp:652-679) ---
    write_genotype_tsv(&genotyper, &args.prefix)?;
    genotyper
        .output_representative_alleles(Path::new(&format!("{}_allele.tsv", args.prefix)), |idx| {
            loaded.names[idx].clone()
        });

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
