#![cfg(feature = "t1k-sys")]
//! End-to-end differential test for Task 6b (the `Analyzer` driver): runs
//! the REAL, unmodified vendored `analyzer` oracle binary AND a from-scratch
//! reimplementation of `fg-t1k analyze`'s driver (calling
//! `fg_t1k_core::genotyper`/`fg_t1k_core::variant_caller` directly -- this
//! crate cannot dev-depend on the `fg-t1k` CLI binary crate without a
//! dependency cycle, since `fg-t1k` itself depends on `fg-t1k-sys`; see
//! `diff_genotype_e2e.rs`'s [`run_rust`]-equivalent doc comment for the same
//! rationale) on the SAME analyzer inputs, both at `-t 1`, and compares the
//! resulting `_allele.vcf`.
//!
//! # How the analyzer inputs (`_aligned_*.fa` + `_allele.tsv`) are produced
//!
//! [`kir_example_analyzer_inputs`] runs the REAL, unmodified vendored
//! `genotyper` oracle binary (`OracleStage::Genotyper`) on
//! `fixtures/example`'s KIR candidate reads -- exactly `run-t1k`'s own first
//! pipeline stage (`vendor/t1k/run-t1k:430`) -- which writes
//! `{prefix}_aligned_1.fa`/`_aligned_2.fa`/`_allele.tsv` (`Genotyper.cpp:
//! 677-707`). Using the oracle (rather than the Rust `genotype` stage, which
//! does not yet write `_aligned_*.fa`) keeps this test's fixture generation
//! itself independently trustworthy, and isolates the differential to
//! EXACTLY the `analyzer` stage this task ports -- `diff_genotype_e2e.rs`
//! already separately proves the Rust `genotype` stage matches the oracle
//! byte-for-byte on this same fixture, so chaining Rust genotype output into
//! this test would not add coverage, only risk conflating two stages' bugs.
//!
//! # Two scenarios
//!
//! - [`kir_example_allele_vcf_matches_oracle`]: the REAL KIR fixture
//!   (`fixtures/example`, 1000 simulated read pairs across 20 called
//!   alleles). The simulated reads are error-free relative to their donor
//!   allele (`snps=0 indels=0` in every FASTQ header, confirmed by
//!   inspection), so BOTH the oracle and this port correctly produce an
//!   EMPTY `_allele.vcf` here -- an important wiring/parity check (restricted
//!   `InitRefSet`, re-quantified abundances, `AddFragmentAlignmentInfo`, all
//!   run end-to-end with zero crashes/mismatches on 1000 real fragments), but
//!   a weak content check on its own, hence the second scenario.
//! - [`synthetic_snv_allele_vcf_matches_oracle`]: a small SYNTHETIC scenario
//!   (one real KIR allele consensus reused verbatim as the single selected
//!   allele, no exon-interval header comment so the whole consensus counts
//!   as one exon -- see `fg_t1k_core::variant_caller`'s `exonic_position`
//!   doc comment) with a deliberately injected single-base substitution
//!   present in 8/10 single-end reads at one position -- well above
//!   `FindCandidateVariants`' both the absolute (`>=5`) and relative
//!   (`>= refCount * 0.5`) support thresholds, so a real novel-variant call
//!   is expected, exercising `_allele.vcf`'s actual variant-row formatting
//!   against the real oracle (not just an empty-file match).
use fg_t1k_core::fastq::FastqReader;
use fg_t1k_core::genotyper::{self, AlleleRef, ExtendedOverlap, Genotyper};
use fg_t1k_core::ref_kmer_filter::RefKmerFilter;
use fg_t1k_core::variant_caller::{self, FragmentOverlap, Overlap, VariantCaller};
use fg_t1k_sys::oracle::{OracleStage, binary_path};
use std::collections::{HashMap, HashSet};
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures").join(rel)
}

/// Runs the real oracle `genotyper` binary on the KIR example fixture to
/// produce `{prefix}_aligned_1.fa`/`_aligned_2.fa`/`_allele.tsv` -- the
/// `analyzer` stage's inputs. See this module's doc comment for why the
/// oracle (not the Rust `genotype` stage) is used to generate them.
fn kir_example_analyzer_inputs(prefix: &Path) {
    let bin = binary_path(OracleStage::Genotyper);
    assert!(bin.exists(), "oracle binary not built: {bin:?}");
    let output = Command::new(&bin)
        .arg("-f")
        .arg(fixture("example/ref/kir_rna_seq.fa"))
        .arg("-1")
        .arg(fixture("example/example_1.fq"))
        .arg("-2")
        .arg(fixture("example/example_2.fq"))
        .arg("-o")
        .arg(prefix)
        .arg("-t")
        .arg("1")
        .output()
        .unwrap_or_else(|e| panic!("spawning genotyper: {e}"));
    assert!(
        output.status.success(),
        "genotyper exited with {:?}; stderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Runs the real oracle `analyzer` binary on `{ref}`/`{allele_file}`/
/// `{aligned_1}`/`{aligned_2}`, writing `{out_prefix}_allele.vcf`.
fn run_oracle_analyzer(
    ref_fasta: &Path,
    allele_file: &Path,
    aligned_1: &Path,
    aligned_2: &Path,
    out_prefix: &Path,
) {
    let bin = binary_path(OracleStage::Analyzer);
    assert!(bin.exists(), "oracle binary not built: {bin:?}");
    let output = Command::new(&bin)
        .arg("-f")
        .arg(ref_fasta)
        .arg("-a")
        .arg(allele_file)
        .arg("-1")
        .arg(aligned_1)
        .arg("-2")
        .arg(aligned_2)
        .arg("-o")
        .arg(out_prefix)
        .arg("-t")
        .arg("1")
        .output()
        .unwrap_or_else(|e| panic!("spawning analyzer: {e}"));
    assert!(
        output.status.success(),
        "analyzer exited with {:?}; stderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
}

/// A self-contained reimplementation of `fg-t1k analyze`'s driver (see
/// `crates/fg-t1k/src/stages/analyze.rs`), calling
/// `fg_t1k_core::genotyper`/`fg_t1k_core::variant_caller` directly instead of
/// spawning the `fg-t1k` binary. See this module's doc comment ("How the
/// analyzer inputs...") for why this mirrors `diff_genotype_e2e.rs`'s
/// `run_rust` precedent rather than spawning a subprocess.
// Long by construction: mirrors Analyzer.cpp:main / stages::analyze::run's
// own too_many_lines allow for the same rationale. Test-local `struct
// Read`/helper fns are declared where first used (co-located with their one
// call site) rather than hoisted, matching diff_genotype_e2e.rs's own
// items_after_statements allow.
#[allow(clippy::too_many_lines, clippy::items_after_statements)]
fn run_rust_analyze(
    ref_fasta: &Path,
    allele_file: &Path,
    aligned_1: &Path,
    aligned_2: &Path,
    out_prefix: &Path,
) {
    const GENOTYPER_KMER_LENGTH: usize = 11;
    const GENE_SIMILARITY_KMER_LENGTH: usize = 31;
    const REF_SEQ_SIMILARITY: f64 = 0.8;
    const VAR_MAX_GROUP: i32 = 8;

    // --- Selected-allele parse (Analyzer.cpp:343-351) ---
    let selected: HashSet<String> = std::fs::read_to_string(allele_file)
        .expect("read allele file")
        .lines()
        .filter_map(|line| line.split_whitespace().next())
        .map(String::from)
        .collect();

    // --- Restricted reference loading (Genotyper::InitRefSet(file, selectedAlleles)) ---
    let text = std::fs::read_to_string(ref_fasta).expect("read reference FASTA");
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
    assert!(!names.is_empty(), "no selected alleles matched the reference FASTA");

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
    let mut allele_refs: Vec<AlleleRef> = consensus
        .iter()
        .zip(&comments)
        .map(|(seq, comment)| AlleleRef::new(seq.clone(), comment.as_deref()))
        .collect();
    let allele_cnt = names.len();

    // --- Build the RefKmerFilter over the selected, deduplicated sequences ---
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

    // --- Read aligned FASTA(s) ---
    struct Read {
        seq: Vec<u8>,
        has_n: bool,
    }
    fn read_all(path: &Path) -> Vec<Read> {
        let mut reader = FastqReader::open(path).expect("open FASTA");
        let mut reads = Vec::new();
        while let Some(rec) = reader.next_record().expect("read record") {
            let has_n = rec.seq.contains(&b'N');
            reads.push(Read { seq: rec.seq, has_n });
        }
        reads
    }
    let reads1 = read_all(aligned_1);
    let reads2 = read_all(aligned_2);
    assert_eq!(reads1.len(), reads2.len(), "mate counts must match");
    let read_cnt = reads1.len();
    let has_mate = true;
    let max_read_length =
        reads1.iter().chain(reads2.iter()).map(|r| r.seq.len()).max().unwrap_or(0);
    g.set_read_length(i32::try_from(max_read_length).unwrap());

    // --- Read-end alignment, deduplicated by exact sequence ---
    let mut all_seqs: Vec<&[u8]> = reads1.iter().map(|r| r.seq.as_slice()).collect();
    all_seqs.extend(reads2.iter().map(|r| r.seq.as_slice()));
    let mut sorted_seqs = all_seqs.clone();
    sorted_seqs.sort_unstable();
    sorted_seqs.dedup();
    let mut overlaps_by_seq: HashMap<&[u8], Option<Vec<ExtendedOverlap>>> = HashMap::new();
    for &seq in &sorted_seqs {
        let raw = filter
            .get_overlaps_from_read(seq, &mut fg_t1k_core::ref_kmer_filter::Scratch::default())
            .unwrap_or_default();
        // Analyzer.cpp:476 uses weight=0 (see stages/analyze.rs); mirror it.
        let extended = genotyper::assign_read(seq, &raw, &mut allele_refs, REF_SEQ_SIMILARITY, 0);
        overlaps_by_seq.insert(seq, extended);
    }

    // --- Fragment assembly + SetReadAssignments + AddFragmentAlignmentInfo inputs ---
    g.init_read_assignments(i32::try_from(read_cnt).unwrap(), 2000);
    let consensus_len_of = |idx: u32| {
        i32::try_from(allele_refs[usize::try_from(idx).unwrap()].consensus.len()).unwrap()
    };
    let hit_len_required = filter.hit_len_required();
    let ref_seq_similarity = filter.ref_seq_similarity();
    let empty: Vec<ExtendedOverlap> = Vec::new();

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

    let mut fragment_assignments: Vec<Vec<FragmentOverlap>> = Vec::with_capacity(read_cnt);
    for i in 0..read_cnt {
        let o1 = overlaps_by_seq.get(reads1[i].seq.as_slice()).and_then(Option::as_ref);
        let o2 = overlaps_by_seq.get(reads2[i].seq.as_slice()).and_then(Option::as_ref);
        let has_n = reads1[i].has_n || reads2[i].has_n;
        let assembled = genotyper::read_assignment_to_fragment_assignment_with_overlaps(
            o1.map_or(empty.as_slice(), Vec::as_slice),
            Some(o2.map_or(empty.as_slice(), Vec::as_slice)),
            has_n,
            hit_len_required,
            consensus_len_of,
            |_, _, _| false,
        );
        let fragment_overlaps: Vec<FragmentOverlap> = assembled
            .iter()
            .map(|(fo, ov1, ov2)| FragmentOverlap {
                seq_idx: fo.seq_idx,
                has_mate_pair: fo.has_mate_pair,
                o1_from_r2: fo.o1_from_r2,
                overlap1: to_overlap(ov1),
                overlap2: ov2.map_or_else(Overlap::none, |ov2| to_overlap(&ov2)),
            })
            .collect();
        let plain: Vec<genotyper::FragmentOverlap> =
            assembled.iter().map(|(fo, _, _)| *fo).collect();
        g.set_read_assignments(i, &plain, ref_seq_similarity, |_, _, _| false);
        fragment_assignments.push(fragment_overlaps);
    }
    g.coalesce_read_assignments(0, i32::try_from(read_cnt.saturating_sub(1)).unwrap());

    // --- FinalizeReadAssignments + quantification (Analyzer.cpp:613,619) ---
    let missing_coverage: Vec<i32> = (0..allele_cnt)
        .map(|i| genotyper::get_seq_missing_base_coverage(&allele_refs[i], 0.01))
        .collect();
    g.finalize_read_assignments(&missing_coverage);
    g.quantify_allele_equivalent_class(&effective_len, &weight);

    // --- VariantCaller: AddFragmentAlignmentInfo + ComputeVariant + OutputAlleleVCF ---
    let mut vc = VariantCaller::new(&allele_refs);
    vc.set_seq_abundance(&g, allele_cnt);
    vc.set_max_var_group_to_resolve(VAR_MAX_GROUP);

    let read1_seqs: Vec<Vec<u8>> = reads1.iter().map(|r| r.seq.clone()).collect();
    let read2_seqs: Vec<Vec<u8>> =
        if has_mate { reads2.iter().map(|r| r.seq.clone()).collect() } else { Vec::new() };
    let allele_consensus: Vec<Vec<u8>> = allele_refs.iter().map(|a| a.consensus.clone()).collect();

    for i in 0..read_cnt {
        let r2 = if has_mate { Some(read2_seqs[i].as_slice()) } else { None };
        variant_caller::add_fragment_alignment_info(
            &read1_seqs[i],
            r2,
            &mut fragment_assignments[i],
            &allele_consensus,
        );
    }
    vc.compute_variant(&read1_seqs, &read2_seqs, &fragment_assignments, &allele_consensus);

    let vcf_path = format!("{}_allele.vcf", out_prefix.display());
    let vcf_text = vc.output_allele_vcf(&names, |seq_idx, pos| {
        variant_caller::exonic_position(&allele_refs[seq_idx].exons, pos)
    });
    std::fs::write(&vcf_path, vcf_text).expect("write _allele.vcf");
}

#[test]
fn kir_example_allele_vcf_matches_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let geno_prefix = dir.path().join("geno");
    kir_example_analyzer_inputs(&geno_prefix);

    let ref_fasta = fixture("example/ref/kir_rna_seq.fa");
    let allele_file = PathBuf::from(format!("{}_allele.tsv", geno_prefix.display()));
    let aligned_1 = PathBuf::from(format!("{}_aligned_1.fa", geno_prefix.display()));
    let aligned_2 = PathBuf::from(format!("{}_aligned_2.fa", geno_prefix.display()));
    assert!(allele_file.exists() && aligned_1.exists() && aligned_2.exists());

    let oracle_prefix = dir.path().join("oracle");
    let rust_prefix = dir.path().join("rust");
    run_oracle_analyzer(&ref_fasta, &allele_file, &aligned_1, &aligned_2, &oracle_prefix);
    run_rust_analyze(&ref_fasta, &allele_file, &aligned_1, &aligned_2, &rust_prefix);

    let oracle_vcf =
        std::fs::read_to_string(format!("{}_allele.vcf", oracle_prefix.display())).unwrap();
    let rust_vcf =
        std::fs::read_to_string(format!("{}_allele.vcf", rust_prefix.display())).unwrap();

    // Confirmed by inspection (see this module's doc comment): the fixture's
    // simulated reads are error-free relative to their donor allele, so an
    // EMPTY `_allele.vcf` is the CORRECT output here, not a false negative --
    // both the oracle and this port must agree it is empty.
    assert_eq!(
        oracle_vcf, "",
        "sanity: oracle is expected to produce an empty VCF on this fixture"
    );
    assert_eq!(rust_vcf, oracle_vcf, "_allele.vcf must be byte-identical to the oracle");
}

/// Builds a synthetic single-allele, single-end analyzer input set (see this
/// module's doc comment for the injected-SNV design) directly on disk under
/// `dir`, returning `(ref_fasta, allele_file, aligned_1, aligned_2)` -- the
/// same 4-path shape [`run_oracle_analyzer`]/[`run_rust_analyze`] take.
/// `aligned_2` is a second copy of the same single-end reads (both sides are
/// invoked in single-end mode via `-u`/one-FASTA `read_all`, so this file is
/// unused by either driver -- it only exists so this helper can return the
/// same 4-tuple shape as the KIR scenario; see the call sites, which pass
/// `aligned_1` as `-u` and never touch `aligned_2`).
fn write_synthetic_snv_inputs(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    // Reuse a REAL KIR allele consensus (no header comment, so the whole
    // sequence is treated as one exon -- see `exonic_position`'s doc
    // comment) as the single selected allele.
    let ref_text = std::fs::read_to_string(fixture("example/ref/kir_rna_seq.fa")).unwrap();
    let first_seq = ref_text.lines().nth(1).expect("reference FASTA has a first sequence body");
    let consensus = first_seq.as_bytes();
    assert!(consensus.len() >= 400, "expected a realistically long KIR allele consensus");

    let name = "SYN1*01";
    let ref_fasta = dir.join("ref.fa");
    std::fs::write(&ref_fasta, format!(">{name}\n{first_seq}\n")).unwrap();

    let allele_file = dir.join("allele.tsv");
    std::fs::write(&allele_file, format!("{name} 60\n")).unwrap();

    // 10 single-end 100bp reads from a fixed consensus window; 8 of them
    // carry a single substitution at the window's midpoint (well clear of
    // both read ends, so it survives ExtendOverlap's overhang extension).
    let window_start = 200usize;
    let window_len = 100usize;
    let mut window = consensus[window_start..window_start + window_len].to_vec();
    let mut_pos = window_len / 2;
    let ref_base = window[mut_pos];
    let alt_base = if ref_base == b'A' { b'T' } else { b'A' };

    let mut records = String::new();
    for i in 0..10 {
        let mut read = window.clone();
        if i < 8 {
            read[mut_pos] = alt_base;
        }
        records.push_str(&format!(">read{i}\n{}\n", String::from_utf8(read).unwrap()));
    }
    let aligned_1 = dir.join("aligned_1.fa");
    std::fs::write(&aligned_1, &records).unwrap();
    // Not read by either driver in single-end (`-u`) mode; written only to
    // keep a stable on-disk artifact set for anyone inspecting the temp dir.
    let aligned_2 = dir.join("aligned_2.fa");
    std::fs::write(&aligned_2, &records).unwrap();

    // Suppress the unused-assignment warning for `window`'s post-mutation
    // reuse across iterations -- silence any accidental future refactor that
    // would otherwise mutate `window` itself in place.
    let _ = &mut window;

    (ref_fasta, allele_file, aligned_1)
}

#[test]
fn synthetic_snv_allele_vcf_matches_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let (ref_fasta, allele_file, aligned_1) = write_synthetic_snv_inputs(dir.path());

    let oracle_prefix = dir.path().join("oracle");
    let rust_prefix = dir.path().join("rust");

    // Single-end (`-u`): run the oracle directly (not via `run_oracle_analyzer`,
    // which is paired-only) and reimplement the single-end shape of
    // `run_rust_analyze` inline for the same reason `diff_genotype_e2e.rs`
    // reimplements the driver rather than spawning `fg-t1k`.
    let bin = binary_path(OracleStage::Analyzer);
    assert!(bin.exists(), "oracle binary not built: {bin:?}");
    let output = Command::new(&bin)
        .arg("-f")
        .arg(&ref_fasta)
        .arg("-a")
        .arg(&allele_file)
        .arg("-u")
        .arg(&aligned_1)
        .arg("-o")
        .arg(&oracle_prefix)
        .arg("-t")
        .arg("1")
        .output()
        .unwrap_or_else(|e| panic!("spawning analyzer: {e}"));
    assert!(
        output.status.success(),
        "analyzer exited with {:?}; stderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );

    run_rust_analyze_single_end(&ref_fasta, &allele_file, &aligned_1, &rust_prefix);

    let oracle_vcf =
        std::fs::read_to_string(format!("{}_allele.vcf", oracle_prefix.display())).unwrap();
    let rust_vcf =
        std::fs::read_to_string(format!("{}_allele.vcf", rust_prefix.display())).unwrap();

    assert!(!oracle_vcf.is_empty(), "sanity: the injected SNV must produce a non-empty oracle VCF");
    assert_eq!(rust_vcf, oracle_vcf, "_allele.vcf must be byte-identical to the oracle");
}

/// Single-end variant of [`run_rust_analyze`] (see that function's doc
/// comment for the overall shape) -- duplicated rather than parameterized
/// since threading an `Option`-everywhere single/paired split through
/// `run_rust_analyze` would obscure the paired-only KIR scenario's own
/// readability for the sake of one additional synthetic caller.
///
/// `needless_range_loop` is allowed on the `0..read_cnt` loop below: the
/// index `i` is used to look up THREE parallel collections
/// (`reads1`/`overlaps_by_seq`/`fragment_assignments`, the last built up
/// incrementally inside the loop itself), which an `.enumerate()` rewrite
/// cannot express as cleanly -- matches `run_rust_analyze`'s identical loop
/// shape (not flagged there only because it also indexes `read2_seqs`,
/// apparently enough to dodge this particular lint heuristic).
#[allow(clippy::too_many_lines, clippy::items_after_statements, clippy::needless_range_loop)]
fn run_rust_analyze_single_end(
    ref_fasta: &Path,
    allele_file: &Path,
    aligned_1: &Path,
    out_prefix: &Path,
) {
    const GENOTYPER_KMER_LENGTH: usize = 11;
    const GENE_SIMILARITY_KMER_LENGTH: usize = 31;
    const REF_SEQ_SIMILARITY: f64 = 0.8;
    const VAR_MAX_GROUP: i32 = 8;

    let selected: HashSet<String> = std::fs::read_to_string(allele_file)
        .expect("read allele file")
        .lines()
        .filter_map(|line| line.split_whitespace().next())
        .map(String::from)
        .collect();

    let text = std::fs::read_to_string(ref_fasta).expect("read reference FASTA");
    let mut names: Vec<String> = Vec::new();
    let mut consensus: Vec<Vec<u8>> = Vec::new();
    let mut weight: Vec<i32> = Vec::new();
    let mut comments: Vec<Option<String>> = Vec::new();
    let mut cur_id: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(id) = cur_id.take() {
                if selected.contains(&id) && !cur_seq.is_empty() {
                    names.push(id);
                    consensus.push(std::mem::take(&mut cur_seq));
                    weight.push(1);
                    comments.push(None);
                }
            }
            cur_seq.clear();
            cur_id = Some(rest.split_whitespace().next().unwrap_or("").to_string());
        } else {
            cur_seq.extend_from_slice(line.trim_end().as_bytes());
        }
    }
    if let Some(id) = cur_id {
        if selected.contains(&id) && !cur_seq.is_empty() {
            names.push(id);
            consensus.push(cur_seq);
            weight.push(1);
            comments.push(None);
        }
    }
    assert!(!names.is_empty(), "no selected alleles matched the reference FASTA");

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
    let mut allele_refs: Vec<AlleleRef> = consensus
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
        let mut reader = FastqReader::open(path).expect("open FASTA");
        let mut reads = Vec::new();
        while let Some(rec) = reader.next_record().expect("read record") {
            let has_n = rec.seq.contains(&b'N');
            reads.push(Read { seq: rec.seq, has_n });
        }
        reads
    }
    let reads1 = read_all(aligned_1);
    let read_cnt = reads1.len();
    g.set_read_length(
        i32::try_from(reads1.iter().map(|r| r.seq.len()).max().unwrap_or(0)).unwrap(),
    );

    let mut sorted_seqs: Vec<&[u8]> = reads1.iter().map(|r| r.seq.as_slice()).collect();
    sorted_seqs.sort_unstable();
    sorted_seqs.dedup();
    let mut overlaps_by_seq: HashMap<&[u8], Option<Vec<ExtendedOverlap>>> = HashMap::new();
    for &seq in &sorted_seqs {
        let raw = filter
            .get_overlaps_from_read(seq, &mut fg_t1k_core::ref_kmer_filter::Scratch::default())
            .unwrap_or_default();
        // Analyzer.cpp:476 uses weight=0 (see stages/analyze.rs); mirror it.
        let extended = genotyper::assign_read(seq, &raw, &mut allele_refs, REF_SEQ_SIMILARITY, 0);
        overlaps_by_seq.insert(seq, extended);
    }

    g.init_read_assignments(i32::try_from(read_cnt).unwrap(), 2000);
    let consensus_len_of = |idx: u32| {
        i32::try_from(allele_refs[usize::try_from(idx).unwrap()].consensus.len()).unwrap()
    };
    let hit_len_required = filter.hit_len_required();
    let ref_seq_similarity = filter.ref_seq_similarity();
    let empty: Vec<ExtendedOverlap> = Vec::new();

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

    let mut fragment_assignments: Vec<Vec<FragmentOverlap>> = Vec::with_capacity(read_cnt);
    for i in 0..read_cnt {
        let o1 = overlaps_by_seq.get(reads1[i].seq.as_slice()).and_then(Option::as_ref);
        let assembled = genotyper::read_assignment_to_fragment_assignment_with_overlaps(
            o1.map_or(empty.as_slice(), Vec::as_slice),
            None,
            reads1[i].has_n,
            hit_len_required,
            consensus_len_of,
            |_, _, _| false,
        );
        let fragment_overlaps: Vec<FragmentOverlap> = assembled
            .iter()
            .map(|(fo, ov1, ov2)| FragmentOverlap {
                seq_idx: fo.seq_idx,
                has_mate_pair: fo.has_mate_pair,
                o1_from_r2: fo.o1_from_r2,
                overlap1: to_overlap(ov1),
                overlap2: ov2.map_or_else(Overlap::none, |ov2| to_overlap(&ov2)),
            })
            .collect();
        let plain: Vec<genotyper::FragmentOverlap> =
            assembled.iter().map(|(fo, _, _)| *fo).collect();
        g.set_read_assignments(i, &plain, ref_seq_similarity, |_, _, _| false);
        fragment_assignments.push(fragment_overlaps);
    }
    g.coalesce_read_assignments(0, i32::try_from(read_cnt.saturating_sub(1)).unwrap());

    let missing_coverage: Vec<i32> = (0..allele_cnt)
        .map(|i| genotyper::get_seq_missing_base_coverage(&allele_refs[i], 0.01))
        .collect();
    g.finalize_read_assignments(&missing_coverage);
    g.quantify_allele_equivalent_class(&effective_len, &weight);

    let mut vc = VariantCaller::new(&allele_refs);
    vc.set_seq_abundance(&g, allele_cnt);
    vc.set_max_var_group_to_resolve(VAR_MAX_GROUP);

    let read1_seqs: Vec<Vec<u8>> = reads1.iter().map(|r| r.seq.clone()).collect();
    let read2_seqs: Vec<Vec<u8>> = Vec::new();
    let allele_consensus: Vec<Vec<u8>> = allele_refs.iter().map(|a| a.consensus.clone()).collect();

    for i in 0..read_cnt {
        variant_caller::add_fragment_alignment_info(
            &read1_seqs[i],
            None,
            &mut fragment_assignments[i],
            &allele_consensus,
        );
    }
    vc.compute_variant(&read1_seqs, &read2_seqs, &fragment_assignments, &allele_consensus);

    let vcf_path = format!("{}_allele.vcf", out_prefix.display());
    let vcf_text = vc.output_allele_vcf(&names, |seq_idx, pos| {
        variant_caller::exonic_position(&allele_refs[seq_idx].exons, pos)
    });
    std::fs::write(&vcf_path, vcf_text).expect("write _allele.vcf");
}
