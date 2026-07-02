//! Golden byte-identity test for the `.dat` -> seq-FASTA port
//! (`fg_t1k_core::refbuild::dat`), ported from T1K's `ParseDatFile.pl`.
//!
//! Runs the Rust parser+emitter on the committed KIR subset fixture and asserts
//! the emitted `dna`/`rna` seq FASTAs are byte-for-byte identical to the golden
//! files produced by the vendored Perl oracle (see `fixtures/refbuild/PINS.md`).

use fg_t1k_core::refbuild::dat::{self, SeqKind};

fn fixture(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/refbuild").join(rel)
}

#[test]
// `golden_dna`/`golden_rna` (and `dna`/`rna`) are intentionally named to pair
// with `SeqKind::Dna`/`SeqKind::Rna` — the similarity is the point, not an
// accident clippy needs to flag.
#[allow(clippy::similar_names)]
fn dat_to_seq_fastas_match_golden() {
    let alleles = dat::parse_dat(&fixture("kir_subset.dat")).expect("parsing kir_subset.dat");

    let dna = dat::emit_seq_fasta(&alleles, SeqKind::Dna).expect("emitting dna seq fasta");
    let rna = dat::emit_seq_fasta(&alleles, SeqKind::Rna).expect("emitting rna seq fasta");

    let golden_dna = std::fs::read_to_string(fixture("golden/kir_dna_seq.fa"))
        .expect("reading golden dna fasta");
    let golden_rna = std::fs::read_to_string(fixture("golden/kir_rna_seq.fa"))
        .expect("reading golden rna fasta");

    assert_eq!(dna, golden_dna, "dna seq FASTA does not match golden byte-for-byte");
    assert_eq!(rna, golden_rna, "rna seq FASTA does not match golden byte-for-byte");
}
