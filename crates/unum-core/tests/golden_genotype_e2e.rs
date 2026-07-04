//! End-to-end golden-file test for the `genotype` capstone, converted from the retired T1K-oracle FFI/subprocess differential (`diff_genotype_e2e.rs`) (see
//! `tests/common/mod.rs`). Runs a from-scratch reimplementation of the
//! `genotype` driver (calling `unum_core::genotyper` directly -- the core
//! crate cannot depend on the `unum` CLI binary crate) over the
//! `fixtures/example` KIR candidate reads + reference, and asserts the
//! resulting `_genotype.tsv` / `_allele.tsv` are byte-identical to the
//! committed goldens (the Rust port's output, which was byte-identical to the
//! vendored `genotyper` oracle on all 17 gene rows).
//!
//! The `_genotype.tsv` golden is `fixtures/example/oracle_genotype.golden.tsv`
//! (the project's pinned Phase-0 golden, see `fixtures/example/PINS.md`); the
//! `_allele.tsv` golden is captured under `tests/golden/genotype_e2e/`. The
//! former oracle-subprocess test is dropped; its byte-identity to T1K is now
//! frozen into these goldens.

mod common;

use common::assert_byte_golden;
use std::collections::HashMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use unum_core::genotyper::{self, AlleleRef, ExtendedOverlap, Genotyper};
use unum_core::ref_kmer_filter::RefKmerFilter;

fn fixture(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures").join(rel)
}

/// A self-contained reimplementation of the `genotype` CLI driver (mirrors
/// `unum`'s `stages::genotype::run` and `Genotyper.cpp:main`), calling
/// `unum_core::genotyper` directly.
#[allow(clippy::too_many_lines, clippy::items_after_statements)]
fn run_rust(out_prefix: &Path) {
    const GENOTYPER_KMER_LENGTH: usize = 11;
    const GENE_SIMILARITY_KMER_LENGTH: usize = 31;
    const REF_SEQ_SIMILARITY: f64 = 0.8;

    let ref_path = fixture("example/ref/kir_rna_seq.fa");
    let text = std::fs::read_to_string(&ref_path).expect("read reference FASTA");

    let mut seq_to_idx: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut names: Vec<String> = Vec::new();
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

    fn compute_effective_len(seq: &[u8]) -> i32 {
        let mut ret = 0i32;
        for (i, &b) in seq.iter().enumerate() {
            if b != b'N' || (i > 0 && seq[i - 1] != b'N') {
                ret += 1;
            }
        }
        ret
    }
    let mut effective_len: Vec<i32> = consensus.iter().map(|s| compute_effective_len(s)).collect();
    let allele_refs: Vec<AlleleRef> = consensus
        .iter()
        .zip(&comments)
        .map(|(seq, comment)| AlleleRef::new(seq.clone(), comment.as_deref()))
        .collect();
    let allele_cnt = names.len();

    let scratch = tempfile::NamedTempFile::new().expect("scratch reference FASTA");
    {
        let mut f = scratch.reopen().expect("reopen scratch");
        for (name, seq) in names.iter().zip(&consensus) {
            writeln!(f, ">{name}").unwrap();
            f.write_all(seq).unwrap();
            writeln!(f).unwrap();
        }
    }
    let mut filter = RefKmerFilter::from_reference_fasta(scratch.path(), GENOTYPER_KMER_LENGTH)
        .expect("build RefKmerFilter");
    filter.set_ref_seq_similarity(REF_SEQ_SIMILARITY);

    let mut g = Genotyper::new();
    g.init_allele_info(
        &names,
        &consensus,
        &weight,
        &mut effective_len,
        GENE_SIMILARITY_KMER_LENGTH,
    );

    struct Read {
        seq: Vec<u8>,
        has_n: bool,
    }
    fn read_all(path: &Path) -> Vec<Read> {
        let mut reader = unum_core::fastq::FastqReader::open(path).expect("open FASTQ");
        let mut reads = Vec::new();
        while let Some(rec) = reader.next_record().expect("read record") {
            let has_n = rec.seq.contains(&b'N');
            reads.push(Read { seq: rec.seq, has_n });
        }
        reads
    }
    let reads1 = read_all(&fixture("example/example_1.fq"));
    let reads2 = read_all(&fixture("example/example_2.fq"));
    assert_eq!(reads1.len(), reads2.len(), "mate counts must match");
    let read_cnt = reads1.len();
    let max_read_length =
        reads1.iter().chain(reads2.iter()).map(|r| r.seq.len()).max().unwrap_or(0);
    g.set_read_length(i32::try_from(max_read_length).unwrap());

    let mut all_seqs: Vec<&[u8]> = reads1.iter().map(|r| r.seq.as_slice()).collect();
    all_seqs.extend(reads2.iter().map(|r| r.seq.as_slice()));
    let mut sorted_seqs = all_seqs.clone();
    sorted_seqs.sort_unstable();
    sorted_seqs.dedup();
    let mut counted: HashMap<&[u8], i32> = HashMap::new();
    for &seq in &all_seqs {
        *counted.entry(seq).or_insert(0) += 1;
    }
    let mut overlaps_by_seq: HashMap<&[u8], Option<Vec<ExtendedOverlap>>> = HashMap::new();
    let mut dp_cache = unum_core::align_algo::DpCache::new();
    for &seq in &sorted_seqs {
        let w = counted[seq];
        let raw = filter
            .get_overlaps_from_read(seq, &mut unum_core::ref_kmer_filter::Scratch::default())
            .unwrap_or_default();
        let extended =
            genotyper::assign_read(seq, &raw, &allele_refs, REF_SEQ_SIMILARITY, w, &mut dp_cache);
        overlaps_by_seq.insert(seq, extended);
    }

    g.init_read_assignments(i32::try_from(read_cnt).unwrap(), 2000);
    let consensus_len_of = |idx: u32| {
        i32::try_from(allele_refs[usize::try_from(idx).unwrap()].consensus.len()).unwrap()
    };
    let hit_len_required = filter.hit_len_required();
    let ref_seq_similarity = filter.ref_seq_similarity();
    let empty: Vec<ExtendedOverlap> = Vec::new();

    for i in 0..read_cnt {
        let o1 = overlaps_by_seq.get(reads1[i].seq.as_slice()).and_then(Option::as_ref);
        let o2 = overlaps_by_seq.get(reads2[i].seq.as_slice()).and_then(Option::as_ref);
        let has_n = reads1[i].has_n || reads2[i].has_n;
        let assignment = genotyper::read_assignment_to_fragment_assignment(
            o1.map_or(empty.as_slice(), Vec::as_slice),
            Some(o2.map_or(empty.as_slice(), Vec::as_slice)),
            has_n,
            hit_len_required,
            consensus_len_of,
            |_, _, _| false,
        );
        g.set_read_assignments(i, &assignment, ref_seq_similarity, |_, _, _| false);
    }
    g.coalesce_read_assignments(0, i32::try_from(read_cnt - 1).unwrap());

    let missing_coverage: Vec<i32> = (0..allele_cnt)
        .map(|i| genotyper::get_seq_missing_base_coverage(&allele_refs[i], 0.01))
        .collect();
    g.finalize_read_assignments(&missing_coverage);
    g.quantify_allele_equivalent_class(&effective_len, &weight);
    g.remove_low_likelihood_allele_in_equivalent_class(|idx| effective_len[idx]);
    g.select_alleles_for_genes(|idx| weight[idx]);

    let genotype_path = format!("{}_genotype.tsv", out_prefix.display());
    let mut out = std::fs::File::create(&genotype_path).expect("create genotype.tsv");
    let gene_cnt = usize::try_from(g.gene_cnt).unwrap();
    for i in 0..gene_cnt {
        let (allele1, allele2, secondary, called_cnt) = g.get_allele_description(i);
        writeln!(out, "{}\t{called_cnt}\t{allele1}\t{allele2}\t{secondary}", g.gene_idx_to_name[i])
            .unwrap();
    }

    let allele_path = format!("{}_allele.tsv", out_prefix.display());
    g.output_representative_alleles(Path::new(&allele_path), |idx| names[idx].clone());
}

#[test]
fn kir_example_genotype_and_allele_tsv_byte_golden() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("rust");
    run_rust(&prefix);

    // _genotype.tsv must be byte-identical to the pinned Phase-0 golden.
    let genotype = std::fs::read(format!("{}_genotype.tsv", prefix.display())).unwrap();
    let golden = std::fs::read(fixture("example/oracle_genotype.golden.tsv")).unwrap();
    assert_eq!(
        genotype, golden,
        "_genotype.tsv must be byte-identical to fixtures/example/oracle_genotype.golden.tsv"
    );

    // _allele.tsv frozen as a byte-golden under tests/golden/.
    let allele = std::fs::read(format!("{}_allele.tsv", prefix.display())).unwrap();
    assert_byte_golden("genotype_e2e/kir_example_allele.tsv", &allele);
}
