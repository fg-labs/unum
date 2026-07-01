//! Rust port of T1K's reference-build pipeline (`vendor/t1k/t1k-build.pl` and the
//! scripts it drives), which turns an IPD-IMGT/HLA or IPD-KIR `.dat` database into
//! the FASTA references T1K's genotyper indexes.

pub mod dat;
pub mod gene_coord;

use std::path::Path;

use anyhow::{Context, Result};

use dat::SeqKind;

/// Builds the four T1K reference files T1K's genotyper indexes from an IPD-IMGT/HLA or
/// IPD-KIR `.dat` database and a genome GTF, reproducing `vendor/t1k/t1k-build.pl` byte-for-byte
/// (for the `-d`/`-g`/`--prefix` path it drives; see the scope note below).
///
/// Writes, into `out_dir`:
/// - `{prefix}_dna_seq.fa` / `{prefix}_rna_seq.fa` — via [`dat::parse_dat`] +
///   [`dat::emit_seq_fasta`], mirroring `t1k-build.pl`'s two `ParseDatFile.pl $dat --mode dna`
///   / `--mode rna` invocations. Both emissions share a single parse of `dat` (parsing is
///   mode-independent — `ParseDatFile.pl` re-parses the whole file per invocation, but the parsed
///   records it produces do not depend on `--mode`, so reusing one `parse_dat` call is
///   byte-identical to the Perl's two independent parses, just without the redundant I/O).
/// - `{prefix}_dna_coord.fa` / `{prefix}_rna_coord.fa` — via [`gene_coord::load_gtf`] +
///   [`gene_coord::annotate`] run against the two just-written seq FASTAs, mirroring
///   `t1k-build.pl`'s `AddGeneCoord.pl {rna,dna}_seq.fa $gtf` invocations (in that order — the
///   order does not affect either output file's contents).
///
/// `out_dir` is created (recursively) if it does not already exist, mirroring `t1k-build.pl`'s
/// `if (!-d $outputDirectory) { make_path $outputDirectory }` guard; unlike `make_path`, which
/// only calls `mkdir` when the check says the directory is missing, `create_dir_all` is already
/// a no-op (not an error) when the directory exists, so the pre-check is redundant in Rust and is
/// not reproduced.
///
/// # Scope
///
/// `t1k-build.pl` also supports `-f` (a pre-built plain-FASTA input, skipping `ParseDatFile.pl`
/// entirely), `--download`, `--target`, `--ignore-partial`, and `--partial-intron-noseq`. None of
/// those are implemented here: this function always takes the `-d`/`-g` `.dat`+GTF path with no
/// gene-target filter and `ParseDatFile.pl`'s hardcoded defaults (see [`dat`]'s module docs),
/// which is the only path the current pipeline (and its pinned KIR fixtures) exercises.
///
/// # Errors
///
/// Returns an error if `out_dir` cannot be created, `dat` or `gtf` cannot be read/parsed, any
/// seq/coord FASTA cannot be emitted (see [`dat::emit_seq_fasta`]'s synthetic-UTR-fallback caveat),
/// or an output file cannot be written.
pub fn build_reference(dat: &Path, gtf: &Path, out_dir: &Path, prefix: &str) -> Result<()> {
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("creating output directory {}", out_dir.display()))?;

    let alleles = dat::parse_dat(dat)?;
    let gtf_table = gene_coord::load_gtf(gtf)?;

    for kind in [SeqKind::Dna, SeqKind::Rna] {
        let tag = match kind {
            SeqKind::Dna => "dna",
            SeqKind::Rna => "rna",
        };

        let seq_fasta = dat::emit_seq_fasta(&alleles, kind)
            .with_context(|| format!("emitting {tag} seq FASTA"))?;
        let seq_path = out_dir.join(format!("{prefix}_{tag}_seq.fa"));
        std::fs::write(&seq_path, &seq_fasta)
            .with_context(|| format!("writing {}", seq_path.display()))?;

        let coord_fasta = gene_coord::annotate(&seq_path, &gtf_table)
            .with_context(|| format!("annotating {tag} seq FASTA with gene coordinates"))?;
        let coord_path = out_dir.join(format!("{prefix}_{tag}_coord.fa"));
        std::fs::write(&coord_path, &coord_fasta)
            .with_context(|| format!("writing {}", coord_path.display()))?;
    }

    Ok(())
}
