#![cfg(feature = "t1k-sys")]
//! End-to-end differential test for Task 5c (the Genotyper capstone): runs
//! the REAL, unmodified vendored `genotyper` oracle binary AND a from-scratch
//! reimplementation of `fg-t1k genotype`'s driver (calling
//! `fg_t1k_core::genotyper` directly -- this crate cannot dev-depend on the
//! `fg-t1k` CLI binary crate without a dependency cycle, since `fg-t1k`
//! itself depends on `fg-t1k-sys`; see [`run_rust`]'s doc comment) on the
//! SAME `fixtures/example` KIR candidate reads + reference, both at `-t 1`,
//! and compares the resulting `_genotype.tsv`.
//!
//! # `_genotype.tsv` / `_allele.tsv` are byte-identical on this fixture
//!
//! Earlier revisions of this test tolerated a single-row (`KIR3DL2`)
//! abundance mismatch in the 6th decimal digit (`44.958453` vs
//! `44.958454`), attributed at the time to the same `-O3`/FMA-contraction
//! float drift documented in `diff_genotyper_em.rs`. That attribution was
//! WRONG: the actual cause was a logic bug in
//! [`fg_t1k_core::genotyper::read_assignment_to_fragment_assignment`] -- a
//! `relaxedTie` fragment-qual branch (`SeqSet.hpp:2514`) that C++ gates
//! behind `ignoreNonExonDiff` (default `false`, and never set by the
//! oracle invocation this test uses) was being applied
//! UNCONDITIONALLY in the port, admitting fragment assignments C++ would
//! have dropped. Once that branch was corrected to match C++'s
//! `ignoreNonExonDiff == false` behavior (i.e. never applied), the
//! KIR3DL2 abundance mismatch disappeared and `_genotype.tsv` became
//! byte-identical to the oracle on all 17 rows -- see that function's doc
//! comment for the full C++ citation.
//!
//! This does NOT contradict `diff_genotyper_em.rs`'s documented `~1e-8`
//! FMA-drift tolerance on randomized scripted EM inputs -- that drift is
//! real and unrelated to this bug; it just was not what was happening on
//! THIS fixture's one previously-mismatched row. This test therefore
//! asserts, per gene row of `_genotype.tsv`:
//! - **Called alleles (major-allele names at rank 0/1) must match EXACTLY.**
//!   This is the genotype-call tier: the actual clinical/scientific output.
//! - **`genotypeQuality` (the `int` field) must match EXACTLY.**
//! - **Abundance (the `%lf`-formatted float field) must match EXACTLY**
//!   (parsed as `f64` and compared for bit-for-bit equality) -- not a loose
//!   tolerance. See [`assert_genotype_tsv_matches`] for the exact per-field
//!   comparison, which additionally asserts every row's raw TSV line is
//!   byte-identical.
//!
//! # Also reproduces the Phase-0 golden's called alleles
//!
//! [`rust_reproduces_golden_called_alleles`] additionally confirms the Rust
//! `fg_t1k_core::genotyper` path reproduces
//! `fixtures/example/oracle_genotype.golden.tsv`'s called alleles (the
//! self-generated Phase-0 reference, `fixtures/example/PINS.md`) -- a second,
//! independent check that this port's calls are not merely self-consistent
//! with a freshly-run oracle binary but match the project's pinned golden
//! too.

use fg_t1k_core::genotyper::{self, AlleleRef, ExtendedOverlap, Genotyper};
use fg_t1k_core::ref_kmer_filter::RefKmerFilter;
use fg_t1k_sys::oracle::{OracleStage, binary_path};
use std::collections::HashMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures").join(rel)
}

/// Runs the real oracle `genotyper` binary on the KIR example fixture,
/// always at `-t 1` (matching this port's only-ever-single-threaded
/// semantics -- see `fg_t1k_core::genotyper`'s module docs).
fn run_oracle(out_prefix: &Path) {
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
        .arg(out_prefix)
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

/// A self-contained reimplementation of `fg-t1k genotype`'s driver (see
/// `crates/fg-t1k/src/stages/genotype.rs`), calling `fg_t1k_core::genotyper`
/// directly instead of spawning the `fg-t1k` binary.
///
/// # Why not spawn the `fg-t1k` binary as a subprocess
///
/// Every other `diff_*_e2e`-style test in this crate (e.g.
/// `diff_bam_extract.rs`) calls the Rust PORT'S LIBRARY FUNCTIONS directly
/// rather than spawning a CLI subprocess -- this crate (`fg-t1k-sys`) is a
/// dependency of the `fg-t1k` binary crate, so a dev-dependency in the other
/// direction would be a cycle. This function is therefore a deliberately
/// slim, test-local copy of `stages::genotype::run`'s driver logic (same
/// reference-loading/dedup, same read-processing loop, same output
/// assembly) -- if it ever drifts from the real CLI driver, that is a test
/// hygiene concern to fix, not a reason to doubt the underlying
/// `fg_t1k_core::genotyper` port itself (which is what both this test and
/// the CLI ultimately call).
// Long by construction: a linear driver mirroring Genotyper.cpp:main /
// stages::genotype::run (see that function's own too_many_lines allow for
// the same rationale) plus its own local helper fns/structs (which clippy's
// items_after_statements would otherwise also flag -- allowed below for the
// same "small test-local helpers co-located with their one call site"
// readability reason).
#[allow(clippy::too_many_lines, clippy::items_after_statements)]
fn run_rust(out_prefix: &Path) {
    const GENOTYPER_KMER_LENGTH: usize = 11;
    const GENE_SIMILARITY_KMER_LENGTH: usize = 31;
    const REF_SEQ_SIMILARITY: f64 = 0.8;

    // --- Reference loading: dedup identical sequences (InitRefSet) ---
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
    let mut allele_refs: Vec<AlleleRef> = consensus
        .iter()
        .zip(&comments)
        .map(|(seq, comment)| AlleleRef::new(seq.clone(), comment.as_deref()))
        .collect();
    let allele_cnt = names.len();

    // --- Build the RefKmerFilter over the deduplicated sequences ---
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

    // --- Read candidate FASTQs ---
    struct Read {
        seq: Vec<u8>,
        has_n: bool,
    }
    fn read_all(path: &Path) -> Vec<Read> {
        let mut reader = fg_t1k_core::fastq::FastqReader::open(path).expect("open FASTQ");
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

    // --- Read-end alignment, deduplicated by exact sequence ---
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
    for &seq in &sorted_seqs {
        let w = counted[seq];
        let raw = filter
            .get_overlaps_from_read(seq, &mut fg_t1k_core::ref_kmer_filter::Scratch::default())
            .unwrap_or_default();
        let extended = genotyper::assign_read(seq, &raw, &mut allele_refs, REF_SEQ_SIMILARITY, w);
        overlaps_by_seq.insert(seq, extended);
    }

    // --- Fragment assembly + SetReadAssignments ---
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

    // --- Finalize, quantify, select ---
    let missing_coverage: Vec<i32> = (0..allele_cnt)
        .map(|i| genotyper::get_seq_missing_base_coverage(&allele_refs[i], 0.01))
        .collect();
    g.finalize_read_assignments(&missing_coverage);
    g.quantify_allele_equivalent_class(&effective_len, &weight);
    g.remove_low_likelihood_allele_in_equivalent_class(|idx| effective_len[idx]);
    g.select_alleles_for_genes(|idx| weight[idx]);

    // --- Output ---
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

/// One parsed `_genotype.tsv` row: the raw line (for whole-line
/// byte-identity comparison), gene name, called-allele count, then the
/// three tab-split allele-description fields, each further split on
/// `\t`/`;` into `(majorAlleleNames, abundance, quality)` triples (rank 0,
/// rank 1, then zero or more secondary candidates).
#[derive(Debug)]
struct GenotypeRow {
    raw_line: String,
    gene: String,
    called_cnt: String,
    fields: Vec<String>, // the 3 raw allele-description fields (unparsed)
}

fn parse_genotype_tsv(path: &Path) -> Vec<GenotypeRow> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("reading {path:?}: {e}"));
    text.lines()
        .map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            // cols[0]=gene, cols[1]=calledCnt, cols[2..] = allele1/allele2/secondary
            // fields, but each of THOSE fields itself contains embedded tabs
            // (abundance/quality), so re-join everything after col 1 and
            // treat allele1/allele2/secondary as 3 comma-token-free ranges
            // is fragile -- instead, keep the raw remainder split by our
            // OWN knowledge of the format: each populated field is exactly
            // 3 tab-separated tokens (name, abundance, quality), and the
            // secondary field (if non-empty) is 1 token (semicolon-internal).
            let gene = cols[0].to_string();
            let called_cnt = cols[1].to_string();
            let rest = &cols[2..];
            // allele1 = rest[0..3], allele2 = rest[3..6], secondary = rest[6..] joined back with '\t'
            let allele1 = rest.get(0..3).map(|s| s.join("\t")).unwrap_or_default();
            let allele2 = rest.get(3..6).map(|s| s.join("\t")).unwrap_or_default();
            let secondary = rest.get(6..).map(|s| s.join("\t")).unwrap_or_default();
            GenotypeRow {
                raw_line: line.to_string(),
                gene,
                called_cnt,
                fields: vec![allele1, allele2, secondary],
            }
        })
        .collect()
}

/// Splits one allele-description field into `(names, abundance:
/// Option<f64>, quality: Option<i32>)`. The primary two fields
/// (allele1/allele2) use `\t` as the separator
/// (`"{names}\t{abundance}\t{quality}"`, or the literal `".\t0\t-1"`); the
/// secondary-candidates field uses `;` instead (`Genotyper::
/// GetAlleleDescription`'s `sep` switches to `';'` once `type > 1`,
/// `Genotyper.hpp:2130-2134`) -- this function tries `\t` first (matching
/// the fixture's own only-ever-single-secondary-group shape, where the `;`
/// field never itself contains an embedded tab) and falls back to `;` so
/// both field shapes parse correctly. Empty input (no secondary candidates)
/// returns all-`None`.
fn split_field(field: &str) -> (String, Option<f64>, Option<i32>) {
    if field.is_empty() {
        return (String::new(), None, None);
    }
    let mut parts: Vec<&str> = field.splitn(3, '\t').collect();
    if parts.len() != 3 {
        parts = field.splitn(3, ';').collect();
    }
    if parts.len() != 3 {
        return (field.to_string(), None, None);
    }
    let names = parts[0].to_string();
    let abund = parts[1].parse::<f64>().ok();
    let qual = parts[2].parse::<i32>().ok();
    (names, abund, qual)
}

/// Compares two parsed `_genotype.tsv` files: gene order, called-allele
/// COUNT, and each of the 3 allele-description fields' major-allele NAMES,
/// QUALITY, and ABUNDANCE, all EXACT (see this module's doc comment for why
/// this fixture is byte-identical against the oracle, not merely
/// within-tolerance). The per-field checks run first so a mismatch reports
/// precisely which field diverged; the whole-line `raw_line` byte-identity
/// assertion runs last as the umbrella check (it would also catch any
/// divergence the field-level parsing might not, e.g. whitespace).
fn assert_genotype_tsv_matches(rust_path: &Path, oracle_path: &Path) {
    let rust_rows = parse_genotype_tsv(rust_path);
    let oracle_rows = parse_genotype_tsv(oracle_path);
    assert_eq!(
        rust_rows.len(),
        oracle_rows.len(),
        "gene row count differs: rust={} oracle={}",
        rust_rows.len(),
        oracle_rows.len()
    );

    let mut any_row_checked = false;

    for (r, o) in rust_rows.iter().zip(&oracle_rows) {
        assert_eq!(r.gene, o.gene, "gene row order mismatch");
        assert_eq!(r.called_cnt, o.called_cnt, "gene {}: called-allele count mismatch", r.gene);

        for (label, rf, of) in [
            ("allele1", &r.fields[0], &o.fields[0]),
            ("allele2", &r.fields[1], &o.fields[1]),
            ("secondary", &r.fields[2], &o.fields[2]),
        ] {
            any_row_checked = true;
            let (r_names, r_abund, r_qual) = split_field(rf);
            let (o_names, o_abund, o_qual) = split_field(of);

            assert_eq!(
                r_names, o_names,
                "gene {} {label}: called major-allele name(s) mismatch",
                r.gene
            );
            assert_eq!(r_qual, o_qual, "gene {} {label}: genotypeQuality mismatch", r.gene);
            match (r_abund, o_abund) {
                (Some(ra), Some(oa)) => {
                    // Exact `f64` equality is intentional here (not a
                    // tolerance check) -- see this module's doc comment for
                    // why this fixture is byte-identical against the
                    // oracle, and `genotyper.rs`'s own `exact_tie` for the
                    // same `#[allow(clippy::float_cmp)]` pattern.
                    #[allow(clippy::float_cmp)]
                    let exact = ra == oa;
                    assert!(
                        exact,
                        "gene {} {label}: abundance mismatch (expected EXACT match; see this \
                         module's doc comment -- if this fires, it likely indicates a real \
                         regression, not float drift): rust={ra} oracle={oa}",
                        r.gene
                    );
                }
                (None, None) => {}
                _ => {
                    panic!("gene {} {label}: one side has an abundance, the other doesn't", r.gene)
                }
            }
        }

        assert_eq!(
            r.raw_line, o.raw_line,
            "gene {}: _genotype.tsv row is not byte-identical to the oracle",
            r.gene
        );
    }

    assert!(any_row_checked, "sanity: the fixture must produce at least one gene row");
}

#[test]
fn kir_example_called_alleles_and_quality_match_oracle() {
    let dir = tempfile::tempdir().unwrap();
    let rust_prefix = dir.path().join("rust");
    let oracle_prefix = dir.path().join("oracle");

    run_rust(&rust_prefix);
    run_oracle(&oracle_prefix);

    let rust_genotype = format!("{}_genotype.tsv", rust_prefix.display());
    let oracle_genotype = format!("{}_genotype.tsv", oracle_prefix.display());
    assert_genotype_tsv_matches(Path::new(&rust_genotype), Path::new(&oracle_genotype));

    // _allele.tsv (the representative-allele output) should also match on
    // called alleles + quality; this file's abundance is NOT printed (only
    // `{name} {quality}`), so a plain line-set comparison suffices and is a
    // stronger check (byte-identical) where it applies.
    let rust_allele =
        std::fs::read_to_string(format!("{}_allele.tsv", rust_prefix.display())).unwrap();
    let oracle_allele =
        std::fs::read_to_string(format!("{}_allele.tsv", oracle_prefix.display())).unwrap();
    assert_eq!(rust_allele, oracle_allele, "_allele.tsv must be byte-identical (no float fields)");
}

#[test]
fn rust_reproduces_golden_called_alleles() {
    let dir = tempfile::tempdir().unwrap();
    let rust_prefix = dir.path().join("rust");
    run_rust(&rust_prefix);

    let rust_genotype = format!("{}_genotype.tsv", rust_prefix.display());
    let golden_path = fixture("example/oracle_genotype.golden.tsv");

    let rust_rows = parse_genotype_tsv(Path::new(&rust_genotype));
    let golden_rows = parse_genotype_tsv(&golden_path);
    assert_eq!(rust_rows.len(), golden_rows.len());

    for (r, gd) in rust_rows.iter().zip(&golden_rows) {
        assert_eq!(r.gene, gd.gene);
        assert_eq!(
            r.called_cnt, gd.called_cnt,
            "gene {}: called-allele count mismatch vs golden",
            r.gene
        );
        for (label, rf, gf) in
            [("allele1", &r.fields[0], &gd.fields[0]), ("allele2", &r.fields[1], &gd.fields[1])]
        {
            let (r_names, _, r_qual) = split_field(rf);
            let (g_names, _, g_qual) = split_field(gf);
            assert_eq!(
                r_names, g_names,
                "gene {} {label}: called major-allele name(s) mismatch vs golden",
                r.gene
            );
            assert_eq!(
                r_qual, g_qual,
                "gene {} {label}: genotypeQuality mismatch vs golden",
                r.gene
            );
        }
    }
}
