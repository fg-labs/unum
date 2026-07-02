//! Port of T1K's `AddGeneCoord.pl`: rewrites each record's header in a seq
//! FASTA (`_dna_seq.fa` / `_rna_seq.fa`, produced by [`super::dat`]) to append
//! that allele's gene's genomic coordinates, as read from a GTF gene
//! annotation file. `t1k-build.pl` feeds this the output of [`super::dat`]
//! and a genome GTF to produce the final `_dna_coord.fa` / `_rna_coord.fa`
//! reference FASTAs.
//!
//! # Scope
//!
//! `AddGeneCoord.pl` supports a `--gtf-gene-name-mapping` CLI flag (a
//! comma-separated list of `GtfName:SeqName` pairs used to translate a GTF
//! `gene_name` into the gene prefix used by allele names in the seq FASTA,
//! e.g. because HLA-HFE's official gene symbol in most GTFs is `HFE`, not
//! `HLA-HFE`) and a `$hasChrPrefix` toggle for whether output chromosome
//! names carry a `chr` prefix. `t1k-build.pl` never passes
//! `--gtf-gene-name-mapping` (so the Perl's own default,
//! `"HFE:HLA-HFE"`, is always in effect) and `$hasChrPrefix` is a hardcoded
//! `1` with no CLI flag to change it at all. This port therefore hardcodes
//! both: [`GENE_NAME_MAPPING`] is always applied, and chromosome names always
//! gain a `chr` prefix if they lack one.
//!
//! # A note on the single-pass design
//!
//! `AddGeneCoord.pl` reads the seq FASTA *twice*: once to seed
//! `%geneCoord` with a `chr19 -1 -1 +` sentinel for every gene mentioned in
//! it, then again (after parsing the GTF) to rewrite headers. The seed pass
//! only matters so the GTF pass's `defined $geneCoord{$gname} && ... == -1`
//! guard can distinguish "first real match for a gene actually used by this
//! seq FASTA" from "a GTF gene irrelevant to this reference". Recording the
//! first `gene`-row coordinate for *every* GTF gene name (used or not, via
//! [`load_gtf`]) and looking each allele's gene up on demand in [`annotate`]
//! produces byte-identical output without the first pass: entries for genes
//! that never appear in the seq FASTA are simply never consulted, so keeping
//! them is harmless.
#![allow(clippy::doc_markdown)]

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};

/// Chromosome written for a gene with no matching `gene` row anywhere in the
/// GTF. Mirrors `AddGeneCoord.pl`'s `$defaultChr = "chr19"` (valid only
/// because every T1K KIR/HLA reference this script has ever been run against
/// lives on chr19/chr6; the Perl has no fallback for any other chromosome).
const DEFAULT_CHROM: &str = "chr19";

/// Start/end coordinate written for a gene with no matching `gene` row
/// anywhere in the GTF. Mirrors `AddGeneCoord.pl`'s unresolved-coordinate
/// sentinel (`-1`).
const UNRESOLVED_COORD: i64 = -1;

/// Strand written for a gene with no matching `gene` row anywhere in the
/// GTF. Mirrors `AddGeneCoord.pl`'s default strand (`+`).
const DEFAULT_STRAND: char = '+';

/// GTF `gene_name` -> seq-FASTA gene-prefix translations, applied while
/// parsing the GTF in [`load_gtf`]. Mirrors `AddGeneCoord.pl`'s
/// `--gtf-gene-name-mapping` default value (`"HFE:HLA-HFE"`), which
/// `t1k-build.pl` never overrides — see the module docs.
const GENE_NAME_MAPPING: &[(&str, &str)] = &[("HFE", "HLA-HFE")];

/// One gene's genomic coordinates, as found in a GTF `gene` row (or the
/// unresolved sentinel, when no such row exists — see [`GeneCoord::unresolved`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneCoord {
    /// Chromosome name, `chr`-prefixed (e.g. `chr19`).
    pub chrom: String,
    /// 1-based, inclusive start coordinate, taken verbatim from the GTF's
    /// `start` column.
    pub start: i64,
    /// 1-based, inclusive end coordinate, taken verbatim from the GTF's `end`
    /// column.
    pub end: i64,
    /// Strand (`+` or `-`), taken verbatim from the GTF's `strand` column.
    pub strand: char,
}

impl GeneCoord {
    /// The sentinel coordinate for a gene with no matching `gene` row
    /// anywhere in the GTF: mirrors `AddGeneCoord.pl`'s
    /// `$geneCoord{$gene} = "$defaultChr -1 -1 +"` seed, which is left
    /// untouched when the GTF pass never finds a match.
    #[must_use]
    pub fn unresolved() -> Self {
        Self {
            chrom: DEFAULT_CHROM.to_string(),
            start: UNRESOLVED_COORD,
            end: UNRESOLVED_COORD,
            strand: DEFAULT_STRAND,
        }
    }
}

impl std::fmt::Display for GeneCoord {
    /// Renders as `AddGeneCoord.pl`'s coordinate suffix does: space-joined
    /// `chrom start end strand`, e.g. `chr19 54769793 54784332 +`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {} {} {}", self.chrom, self.start, self.end, self.strand)
    }
}

/// Per-gene coordinate table parsed from a GTF, keyed by gene name (already
/// translated through [`GENE_NAME_MAPPING`]).
///
/// Built by [`load_gtf`]; see the module docs for why this holds an entry for
/// every GTF `gene` row rather than only the genes a particular seq FASTA
/// mentions.
#[derive(Debug, Default)]
pub struct GeneCoordTable {
    coords: HashMap<String, GeneCoord>,
}

impl GeneCoordTable {
    /// Looks up `gene`'s coordinates, if the GTF had a `gene` row for it.
    #[must_use]
    pub fn get(&self, gene: &str) -> Option<&GeneCoord> {
        self.coords.get(gene)
    }
}

/// Parses a GTF gene-annotation file into a [`GeneCoordTable`].
///
/// Only `gene`-feature rows (column 3, 0-based column 2) are consumed;
/// comment lines (starting with `#`) and every other feature type
/// (`transcript`, `exon`, `CDS`, ...) are skipped, mirroring
/// `AddGeneCoord.pl`'s `next if ($cols[2] ne "gene")` guard. For each `gene`
/// row: the gene name is extracted from the `gene_name "..."` attribute
/// (column 9), translated through [`GENE_NAME_MAPPING`], and used as the key
/// under which `(chrom, start, end, strand)` is recorded — but only the
/// *first* `gene` row seen for a given (post-mapping) gene name wins; later
/// rows for the same name are ignored, mirroring `AddGeneCoord.pl`'s
/// `... == -1` "not yet resolved" guard.
///
/// # Errors
///
/// Returns an error if `path` cannot be opened or read, if a `gene` row has
/// no `gene_name "..."` attribute (mirrors `AddGeneCoord.pl`'s
/// `die "No gene_name", $_, "\n"`), or if a `gene` row's start/end coordinate
/// column is not a valid integer (rather than silently recording the
/// `-1 -1` "unresolved" sentinel, which is reserved for genes with no `gene`
/// row at all).
pub fn load_gtf(path: &Path) -> Result<GeneCoordTable> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("opening GTF file {}", path.display()))?;
    load_gtf_reader(std::io::BufReader::new(file))
        .with_context(|| format!("parsing GTF file {}", path.display()))
}

/// State machine driving [`load_gtf`], factored out so it can be exercised on
/// any [`BufRead`] (a real file or an in-memory buffer in tests).
fn load_gtf_reader<R: BufRead>(reader: R) -> Result<GeneCoordTable> {
    let gene_name_mapping: HashMap<&str, &str> = GENE_NAME_MAPPING.iter().copied().collect();
    let mut coords: HashMap<String, GeneCoord> = HashMap::new();

    for line in reader.lines() {
        let line = line.context("reading a line")?;
        if line.starts_with('#') {
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 || cols[2] != "gene" {
            continue;
        }

        let gtf_gene_name = extract_gene_name(cols[8])
            .with_context(|| format!("no gene_name attribute in GTF row: {line}"))?;
        let gene_name = gene_name_mapping.get(gtf_gene_name).copied().unwrap_or(gtf_gene_name);

        // First `gene` row for this (post-mapping) name wins; later rows for
        // the same name are ignored. Parse the coordinates eagerly (only when
        // this is the first row for the gene) and surface a malformed
        // start/end as an error, rather than coercing it to the `-1`
        // sentinel that means "gene absent from GTF".
        if let std::collections::hash_map::Entry::Vacant(entry) =
            coords.entry(gene_name.to_string())
        {
            let start: i64 = cols[3]
                .parse()
                .with_context(|| format!("unparsable gene start coordinate in GTF row: {line}"))?;
            let end: i64 = cols[4]
                .parse()
                .with_context(|| format!("unparsable gene end coordinate in GTF row: {line}"))?;
            entry.insert(GeneCoord {
                chrom: add_chr_prefix(cols[0]),
                start,
                end,
                strand: cols[6].chars().next().unwrap_or(DEFAULT_STRAND),
            });
        }
    }

    Ok(GeneCoordTable { coords })
}

/// Extracts the value of a `gene_name "..."` attribute from a GTF row's
/// attribute column, e.g. `KIR2DL1` from
/// `gene_id "..."; gene_name "KIR2DL1"; ...`. Mirrors the Perl regex
/// `/gene_name \"(.*?)\"/` by taking everything up to the *first* closing
/// quote after `gene_name "`.
fn extract_gene_name(attributes: &str) -> Option<&str> {
    const KEY: &str = "gene_name \"";
    let after_key = &attributes[attributes.find(KEY)? + KEY.len()..];
    let end = after_key.find('"')?;
    Some(&after_key[..end])
}

/// Prefixes `chrom` with `chr` unless it already starts with `c`. Mirrors
/// `AddGeneCoord.pl`'s `if ($hasChrPrefix == 1 && !($cols[0] =~ /^c/))`
/// branch (the `$hasChrPrefix == 0` branch, which strips an existing `chr`
/// prefix, is hardcoded off — see the module docs).
fn add_chr_prefix(chrom: &str) -> String {
    if chrom.starts_with('c') { chrom.to_string() } else { format!("chr{chrom}") }
}

/// The gene/locus prefix of an allele ID (text before the first `*`), e.g.
/// `KIR2DL1` for `KIR2DL1*0010101`.
fn gene_of(allele_id: &str) -> &str {
    allele_id.split('*').next().unwrap_or(allele_id)
}

/// Rewrites every record header in the seq FASTA at `seq_fasta_path` to
/// append its gene's genomic coordinates from `gtf`, returning the resulting
/// FASTA text.
///
/// For each header line: the first whitespace-separated token (including its
/// leading `>`) is kept as-is and everything else on that header line (e.g.
/// the seq-FASTA's own exon-boundary integer list) is dropped; the allele's
/// gene (its ID's prefix before `*`) is looked up in `gtf`, falling back to
/// [`GeneCoord::unresolved`] if the gene has no GTF `gene` row; the looked-up
/// coordinate is appended as a space-separated suffix. Sequence lines are
/// passed through unchanged, with multi-line records concatenated onto a
/// single output line (mirroring `AddGeneCoord.pl`'s `$seq .= $_`
/// accumulation). Mirrors `AddGeneCoord.pl`'s final read-and-print pass.
///
/// # Errors
///
/// Returns an error if `seq_fasta_path` cannot be opened or read.
pub fn annotate(seq_fasta_path: &Path, gtf: &GeneCoordTable) -> Result<String> {
    let file = std::fs::File::open(seq_fasta_path)
        .with_context(|| format!("opening seq FASTA {}", seq_fasta_path.display()))?;
    annotate_reader(std::io::BufReader::new(file), gtf)
        .with_context(|| format!("annotating seq FASTA {}", seq_fasta_path.display()))
}

/// State machine driving [`annotate`], factored out so it can be exercised on
/// any [`BufRead`] (a real file or an in-memory buffer in tests).
fn annotate_reader<R: BufRead>(reader: R, gtf: &GeneCoordTable) -> Result<String> {
    let unresolved = GeneCoord::unresolved();
    let mut out = String::new();
    let mut seq = String::new();

    for line in reader.lines() {
        let line = line.context("reading a line")?;
        if line.starts_with('>') {
            // Flush the previous record's accumulated sequence line, if any
            // (mirrors `print $seq, "\n" if ($seq ne "")`).
            if !seq.is_empty() {
                out.push_str(&seq);
                out.push('\n');
                seq.clear();
            }

            let header_token = line.split_whitespace().next().unwrap_or(&line);
            let allele_id = header_token.strip_prefix('>').unwrap_or(header_token);
            let gene = gene_of(allele_id);
            let coord = gtf.get(gene).unwrap_or(&unresolved);
            out.push_str(header_token);
            out.push(' ');
            out.push_str(&coord.to_string());
            out.push('\n');
        } else {
            seq.push_str(&line);
        }
    }
    if !seq.is_empty() {
        out.push_str(&seq);
        out.push('\n');
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The `HFE:HLA-HFE` default mapping: a GTF `gene_name "HFE"` row must be
    /// looked up under `HLA-HFE`, not `HFE`, matching
    /// `AddGeneCoord.pl`'s default `--gtf-gene-name-mapping`.
    #[test]
    fn hfe_gene_name_is_mapped_to_hla_hfe() {
        let gtf = "\
chr6\tHAVANA\tgene\t26087281\t26098343\t.\t+\t.\tgene_id \"ENSG00000010704.18\"; gene_name \"HFE\";
";
        let table = load_gtf_reader(gtf.as_bytes()).expect("parsing synthetic GTF");
        assert!(table.get("HFE").is_none(), "raw GTF gene name must not be a lookup key");
        assert_eq!(
            table.get("HLA-HFE"),
            Some(&GeneCoord {
                chrom: "chr6".to_string(),
                start: 26_087_281,
                end: 26_098_343,
                strand: '+',
            })
        );
    }

    /// An allele whose gene has no matching GTF `gene` row falls back to the
    /// `chr19 -1 -1 +` sentinel, matching `AddGeneCoord.pl`'s untouched
    /// `%geneCoord` seed value for a gene the GTF pass never resolves.
    #[test]
    fn gene_not_found_in_gtf_falls_back_to_unresolved_sentinel() {
        let gtf = load_gtf_reader(std::io::empty()).expect("parsing empty GTF");
        let fasta = ">NOGENE*0010101 8 50 83\nACGTACGT\n";
        let out = annotate_reader(fasta.as_bytes(), &gtf).expect("annotating");
        assert_eq!(out, ">NOGENE*0010101 chr19 -1 -1 +\nACGTACGT\n");
    }

    /// A matching gene's coordinates are appended verbatim, and any trailing
    /// content on the original header line (e.g. the seq-FASTA's own
    /// exon-boundary list) is dropped.
    #[test]
    fn matching_gene_coordinate_replaces_trailing_header_content() {
        let gtf_text = "\
chr19\tHAVANA\tgene\t54769793\t54784332\t.\t+\t.\tgene_name \"KIR2DL1\";
";
        let gtf = load_gtf_reader(gtf_text.as_bytes()).expect("parsing synthetic GTF");
        let fasta = ">KIR2DL1*0010101 8 50 83 84 119\nACGT\n";
        let out = annotate_reader(fasta.as_bytes(), &gtf).expect("annotating");
        assert_eq!(out, ">KIR2DL1*0010101 chr19 54769793 54784332 +\nACGT\n");
    }

    /// Multi-line sequence records are concatenated onto a single output
    /// line, matching `AddGeneCoord.pl`'s `$seq .= $_` accumulation.
    #[test]
    fn multiline_sequence_is_concatenated_onto_one_output_line() {
        let gtf = load_gtf_reader(std::io::empty()).expect("parsing empty GTF");
        let fasta = ">GENE*0010101\nACGT\nTTTT\nGGGG\n";
        let out = annotate_reader(fasta.as_bytes(), &gtf).expect("annotating");
        assert_eq!(out, ">GENE*0010101 chr19 -1 -1 +\nACGTTTTTGGGG\n");
    }

    /// A GTF row with no `gene_name` attribute at all is a hard error,
    /// matching `AddGeneCoord.pl`'s `die "No gene_name", $_, "\n"`.
    #[test]
    fn gene_row_with_no_gene_name_attribute_errors() {
        let gtf = "chr19\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id \"X\";\n";
        assert!(load_gtf_reader(gtf.as_bytes()).is_err());
    }

    /// A `gene` row whose start/end coordinate column is not an integer is a
    /// hard error, not silently coerced to the `-1 -1` "unresolved" sentinel
    /// (which is reserved for genes with no `gene` row at all).
    #[test]
    fn gene_row_with_unparsable_coordinate_errors() {
        let gtf = "chr19\tHAVANA\tgene\tnope\t100\t.\t+\t.\tgene_name \"KIR2DL1\";\n";
        assert!(load_gtf_reader(gtf.as_bytes()).is_err());
    }
}
