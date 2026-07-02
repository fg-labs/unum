//! Golden byte-identity test for the seq-FASTA -> coord-FASTA port
//! (`fg_t1k_core::refbuild::gene_coord`), ported from T1K's
//! `AddGeneCoord.pl`.
//!
//! Runs the Rust GTF parser + header annotator on the committed golden seq
//! FASTAs (isolating this step from the `.dat` -> seq-FASTA build tested in
//! `refbuild_dat.rs`) and asserts the output is byte-for-byte identical to
//! the golden coord FASTAs produced by the vendored Perl oracle (see
//! `fixtures/refbuild/PINS.md`).

use fg_t1k_core::refbuild::gene_coord;

fn fixture(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/refbuild").join(rel)
}

#[test]
fn coord_fastas_match_golden() {
    let gtf = gene_coord::load_gtf(&fixture("kir_genes.gtf")).expect("parsing kir_genes.gtf");

    for (seq, coord) in
        [("kir_dna_seq.fa", "kir_dna_coord.fa"), ("kir_rna_seq.fa", "kir_rna_coord.fa")]
    {
        let out = gene_coord::annotate(&fixture(&format!("golden/{seq}")), &gtf)
            .unwrap_or_else(|e| panic!("annotating {seq}: {e:#}"));
        let golden = std::fs::read_to_string(fixture(&format!("golden/{coord}")))
            .unwrap_or_else(|e| panic!("reading golden {coord}: {e:#}"));
        assert_eq!(out, golden, "{coord} does not match golden byte-for-byte");
    }
}
