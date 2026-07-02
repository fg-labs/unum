//! Port of T1K's `ParseDatFile.pl`: parses EMBL-style IPD-IMGT/HLA or IPD-KIR
//! `.dat` records and emits the `_dna_seq.fa` (full padded genomic sequence) and
//! `_rna_seq.fa` (spliced, exon-only sequence) reference FASTAs that
//! `vendor/t1k/t1k-build.pl` feeds to the rest of the T1K reference-build
//! pipeline.
//!
//! # Scope
//!
//! `ParseDatFile.pl` supports a `--mode genome` in addition to `dna`/`rna`, plus a
//! handful of CLI flags (`-f`, `--gene`, `--ignorePartial`, `--partialInRnaMode`,
//! `--partialIntronHasNoSeq`, `--intronPadding`, `--dedup`). `t1k-build.pl` only
//! ever invokes it as `--mode dna` / `--mode rna` with none of those flags set, so
//! this port hardcodes the equivalent defaults (`$intronPaddingLength = 200`,
//! `$includePartialDiffLen = 0`, `$ignorePartial = 0`,
//! `$partialIntronHasNoSeq = 0`, `$dedup = 0`, no gene-prefix filter) and does not
//! implement `--mode genome`, since nothing in the current pipeline drives it.
//! [`SeqKind`] therefore only has `Dna`/`Rna` variants.
//!
//! One further behavior is intentionally *not* replicated bit-for-bit: when a
//! gene has **no** allele anywhere in the input with real 5'/3' flanking sequence
//! (every allele of that gene needs synthetic UTR padding), `ParseDatFile.pl`
//! falls back to `srand(17)`-seeded `rand()` calls to fabricate a random 50bp UTR.
//! Perl's `rand()`/`srand()` stream is not a portably-specified algorithm, so
//! reproducing it bit-for-bit in Rust cannot be verified without a fixture that
//! actually exercises the fallback (the pinned `kir_subset.dat` fixture does not:
//! every gene in it has at least one full-length record providing real UTR
//! sequence — see `fixtures/refbuild/PINS.md`). Rather than fabricate an
//! unverifiable RNG replica, [`emit_seq_fasta`] returns an error if this path
//! would ever be needed; revisit if/when a fixture requiring it is added.
//!
//! # A note on `i64`/`usize` casts
//!
//! This module intentionally represents every coordinate as `i64`, mirroring
//! Perl's untyped scalar arithmetic (`ParseDatFile.pl` routinely computes
//! transient negative values, e.g. a 5'-UTR start before it gets clamped to
//! `0`), while indexing into a `Vec`/`&str` naturally requires `usize`. Real
//! HLA/KIR `.dat` records are at most tens of kilobases, so every cast here is
//! nowhere near `i64`/`usize` truncation, sign-loss, or wraparound territory
//! on any target this crate supports; the alternative
//! (`i64::try_from(...).expect(...)` at every arithmetic site) would obscure
//! this port's fidelity to the original coordinate math without adding real
//! safety. `clippy::cast_possible_truncation` / `cast_sign_loss` /
//! `cast_possible_wrap` are allowed module-wide for this reason.
#![allow(clippy::cast_possible_truncation, clippy::cast_sign_loss, clippy::cast_possible_wrap)]

use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::path::Path;

use anyhow::{Context, Result};

/// Nominal length, in bases, of the synthetic 5'/3' UTR padding appended to every
/// allele's output sequence. Mirrors `ParseDatFile.pl`'s `$utrLength = 50`, which
/// is only ever overridden (to `0`) by `--mode genome` — not implemented here
/// (see module docs).
const UTR_LENGTH: i64 = 50;

/// Flanking padding (in bases) kept on each side of an intron when assembling
/// [`SeqKind::Dna`] output; introns short enough that two exons' padded windows
/// would overlap are kept whole instead (no `N` separator). Mirrors
/// `ParseDatFile.pl`'s `$intronPaddingLength`, whose default of `200` is never
/// overridden by `t1k-build.pl` (`--intronPadding` is not passed).
const INTRON_PADDING_LENGTH: i64 = 200;

/// How many bases short of a gene's modal effective length a partial allele may
/// be and still be rescued into the output. Mirrors `$includePartialDiffLen`;
/// `t1k-build.pl` never passes `--partialInRnaMode`, so this is always `0`.
const INCLUDE_PARTIAL_DIFF_LEN: i64 = 0;

/// Which spliced-vs-genomic view of an allele's sequence to build.
///
/// Mirrors `ParseDatFile.pl`'s `--mode` flag, restricted to the two modes
/// `t1k-build.pl` actually drives; see the module docs for why `--mode genome`
/// is out of scope.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqKind {
    /// Full genomic sequence: introns are kept, padded to
    /// [`INTRON_PADDING_LENGTH`] bases on each side of an exon, with adjacent
    /// exons merged (no `N` separator, full intron kept) whenever their padded
    /// windows would otherwise overlap.
    Dna,
    /// Spliced, exon-only sequence — introns are removed entirely.
    Rna,
}

/// One EMBL-style `.dat` record: a single allele's raw feature-table coordinates
/// and genomic sequence, parsed directly from the `ID`/`FT`/`SQ` lines of one
/// record. Independent of [`SeqKind`] — the same `RawAllele`s feed both the
/// `dna` and `rna` builds.
#[derive(Debug, Clone)]
pub struct RawAllele {
    /// The allele name from the `FT   /allele="..."` qualifier, e.g.
    /// `KIR2DL1*0010101`.
    pub name: String,
    /// 0-based, inclusive exon coordinate pairs in file (annotation) order.
    ///
    /// An exon immediately followed by a `/pseudo` qualifier is *not* included
    /// here: `ParseDatFile.pl`'s `elsif (/pseudo$/)` branch pops the
    /// just-pushed exon rather than emitting it, effectively folding a
    /// pseudo-exon back into its surrounding intron sequence.
    pub exons: Vec<(i64, i64)>,
    /// Whether any feature in this record carried a trailing `/partial`
    /// qualifier (`ParseDatFile.pl`'s `elsif (/partial$/)` branch). This is
    /// only one of several ways an allele can end up excluded/rescued — see
    /// [`emit_seq_fasta`]'s dna-mode adjacency check and its "exon runs past
    /// the end of the sequence" check, which are mode-specific and computed
    /// later, not stored on `RawAllele`.
    pub annotated_partial: bool,
    /// The concatenated, as-read genomic sequence from the `SQ` block (kept in
    /// its original case, typically lowercase; upper-cased on demand when
    /// building output, matching every `uc(substr(...))` call in the Perl).
    pub sequence: String,
}

impl RawAllele {
    /// The gene/locus prefix of [`Self::name`] — the text before the first
    /// `*`, e.g. `KIR2DL1` for `KIR2DL1*0010101`.
    #[must_use]
    pub fn gene(&self) -> &str {
        self.name.split('*').next().unwrap_or(&self.name)
    }
}

/// Parses every allele record out of an EMBL-style `.dat` file.
///
/// Records with no `allele="..."` qualifier or no `exon` features at all are
/// silently dropped, matching `ParseDatFile.pl`'s `last if ($allele eq "-1")`
/// / `last if (scalar(@exons) == 0)` guards (such a record contributes nothing
/// to any downstream hash in the Perl).
///
/// # Errors
///
/// Returns an error if `path` cannot be opened or read.
pub fn parse_dat(path: &Path) -> Result<Vec<RawAllele>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("opening .dat file {}", path.display()))?;
    parse_dat_reader(std::io::BufReader::new(file))
        .with_context(|| format!("parsing .dat file {}", path.display()))
}

/// State machine driving [`parse_dat`], factored out so it can be exercised on
/// any [`std::io::BufRead`] (a real file or an in-memory buffer in tests).
fn parse_dat_reader<R: std::io::BufRead>(reader: R) -> Result<Vec<RawAllele>> {
    /// Which section of the current `ID ... //` record we're inside.
    enum Section {
        /// Before the `SQ` line: `FT` feature-table lines are meaningful here.
        Header,
        /// After the `SQ` header line, before the closing `//`: every line is
        /// raw sequence data (the "skip the header" comment in the Perl).
        Sequence,
    }

    let mut out = Vec::new();
    let mut section = Section::Header;

    let mut allele: Option<String> = None;
    let mut exons: Vec<(i64, i64)> = Vec::new();
    let mut annotated_partial = false;
    let mut sequence = String::new();

    for line in reader.lines() {
        let line = line.context("reading a line")?;

        if line.starts_with("ID") {
            allele = None;
            exons.clear();
            annotated_partial = false;
            sequence.clear();
            section = Section::Header;
            continue;
        }

        match section {
            Section::Header => {
                if line.starts_with("FT") {
                    process_feature_table_line(
                        &line,
                        &mut allele,
                        &mut exons,
                        &mut annotated_partial,
                    );
                } else if line.starts_with("SQ") {
                    // The `SQ   Sequence NNN BP; ...` line itself carries no
                    // sequence data; only the lines that follow it do.
                    section = Section::Sequence;
                }
            }
            Section::Sequence => {
                if line.starts_with("//") {
                    if let (Some(name), false) = (allele.take(), exons.is_empty()) {
                        out.push(RawAllele {
                            name,
                            exons: std::mem::take(&mut exons),
                            annotated_partial,
                            sequence: std::mem::take(&mut sequence),
                        });
                    }
                    exons.clear();
                    annotated_partial = false;
                    sequence.clear();
                    section = Section::Header;
                } else {
                    append_sequence_line(&line, &mut sequence);
                }
            }
        }
    }

    Ok(out)
}

/// Updates per-record parse state from one `FT` line, mirroring the
/// `elsif (/^FT/) { if (/allele=.../) ... elsif (/\sexon\s/) ... }` chain in
/// `ParseDatFile.pl`. Branches are checked in the same order as the Perl and
/// are mutually exclusive (each real `.dat` line matches at most one).
///
/// The exon-coordinate adjustment `$start - 1 - $partialIntronLen` in the Perl
/// always reduces to `start - 1` here: `$partialIntronLen` only ever becomes
/// nonzero when `--partialIntronHasNoSeq` is passed, which `t1k-build.pl` never
/// does (see module docs).
fn process_feature_table_line(
    line: &str,
    allele: &mut Option<String>,
    exons: &mut Vec<(i64, i64)>,
    annotated_partial: &mut bool,
) {
    if let Some(name) = extract_allele_name(line) {
        *allele = Some(name);
    } else if line.contains(" exon ") {
        if let Some((start, end)) = parse_exon_coords(line) {
            exons.push((start - 1, end - 1));
        }
    } else if line.trim_end().ends_with("pseudo") {
        // A `/pseudo` qualifier immediately follows the exon it annotates;
        // fold it back out of the exon list entirely (it is not spliced/kept
        // as a distinct exon downstream).
        exons.pop();
    } else if line.contains(" intron ") {
        // Recognized (to keep branch order faithful to the Perl) but not
        // modeled further: `$hasIntron` only matters for `--mode genome`.
    } else if line.trim_end().ends_with("partial") {
        *annotated_partial = true;
    }
}

/// Extracts the value of an `allele="..."` qualifier from an `FT` line, e.g.
/// `KIR2DL1*0010101` from `FT                   /allele="KIR2DL1*0010101"`.
/// Mirrors the non-greedy Perl regex `/allele="(.*?)"/` by taking everything up
/// to the *first* closing quote after `allele="`.
fn extract_allele_name(line: &str) -> Option<String> {
    const KEY: &str = "allele=\"";
    let after_key = &line[line.find(KEY)? + KEY.len()..];
    let end = after_key.find('"')?;
    Some(after_key[..end].to_string())
}

/// Extracts the first two digit runs from the third whitespace-separated token
/// of an `exon` feature line (e.g. `269..302` from
/// `FT   exon            269..302`), mirroring
/// `($cols[2] =~ /(\d+)\.\.(\d+)/)` in the Perl. Scanning for digit runs
/// (rather than assuming an exact `NNN..MMM` shape) matches the Perl regex's
/// unanchored search within that token.
fn parse_exon_coords(line: &str) -> Option<(i64, i64)> {
    let token = line.split_whitespace().nth(2)?;
    let mut numbers = token.split(|c: char| !c.is_ascii_digit()).filter(|s| !s.is_empty());
    let start: i64 = numbers.next()?.parse().ok()?;
    let end: i64 = numbers.next()?.parse().ok()?;
    Some((start, end))
}

/// Appends one `SQ`-block data line's bases to `sequence`, dropping the
/// trailing running-length column (e.g. the `60` in
/// `     ggttcttctt gctgcagggg ...        60`). Mirrors
/// `foreach my $s (@cols[0..$#cols - 1]) { $seq .= $s }` after
/// `split /\s+/, $_`.
fn append_sequence_line(line: &str, sequence: &mut String) {
    let mut tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.is_empty() {
        return;
    }
    tokens.pop();
    for token in tokens {
        sequence.push_str(token);
    }
}

// ---------------------------------------------------------------------------
// FindMode: the tie-break rule shared by every "modal length" / "modal intron
// sequence" vote in `ParseDatFile.pl`.
// ---------------------------------------------------------------------------

/// Shared tie-break fold behind [`find_mode_str`] and [`find_mode_i64`].
///
/// Picks the entry with the highest count; ties are broken by comparing the
/// **string** form of the key (`$k ge $ret` in Perl — Perl hash keys are
/// always strings, so even a dist keyed by lengths breaks ties
/// lexicographically on the key's decimal text, not numerically). For example,
/// among two keys tied at the max count, `"9"` beats `"10"`, since `'9' >
/// '1'` as the first differing byte, even though `10 > 9` numerically. This
/// quirk is preserved faithfully because it can affect real output bytes
/// whenever a genuine multi-way tie occurs.
///
/// This fold is commutative/associative over its inputs (the running "best"
/// only ever moves to a strictly-higher count, or to a lexicographically
/// not-smaller key on an exact count tie), so it produces the same result
/// regardless of iteration order — unlike `ParseDatFile.pl`'s
/// partial-allele-rescue order (see [`emit_seq_fasta`]), this does not need a
/// separately-imposed deterministic order.
fn find_mode_by(entries: impl Iterator<Item = (String, i64)>) -> Option<String> {
    let mut best: Option<(String, i64)> = None;
    for (key, count) in entries {
        let replace = match &best {
            None => true,
            Some((best_key, best_count)) => {
                count > *best_count || (count == *best_count && key.as_str() >= best_key.as_str())
            }
        };
        if replace {
            best = Some((key, count));
        }
    }
    best.map(|(key, _)| key)
}

/// String-keyed mode, used for `%geneIntronSeq`'s "most common intron text"
/// vote, where hash keys are literal DNA sequence strings.
fn find_mode_str(counts: &HashMap<String, i64>) -> Option<String> {
    find_mode_by(counts.iter().map(|(key, &count)| (key.clone(), count)))
}

/// Integer-keyed mode, used for every other dist in `ParseDatFile.pl` (gene
/// effective length, exon count, per-exon length, per-intron length, per-allele
/// output length). See [`find_mode_by`] for why ties are broken on the keys'
/// *decimal string* form rather than their numeric value.
fn find_mode_i64(counts: &HashMap<i64, i64>) -> Option<i64> {
    find_mode_by(counts.iter().map(|(&key, &count)| (key.to_string(), count)))
        .map(|key| key.parse().expect("mode key was formatted from an i64"))
}

/// The gene/locus prefix of an allele name (text before the first `*`).
fn gene_of(allele: &str) -> String {
    allele.split('*').next().unwrap_or(allele).to_string()
}

/// `substr($seq, $offset, $length)` for the shapes this module needs
/// (`offset >= 0`, `length >= 0` by construction at every call site); clamps
/// to `seq`'s bounds defensively rather than panicking.
fn perl_substr(seq: &str, offset: i64, length: i64) -> &str {
    let len = seq.len() as i64;
    let start = offset.clamp(0, len) as usize;
    let end = (offset + length.max(0)).clamp(0, len) as usize;
    &seq[start..end.max(start)]
}

/// The last `n` bytes of `s`, mirroring the single call site that needs
/// two-argument `substr($s, -$n)` (taking a gene's 3'-UTR padding suffix).
fn perl_suffix(s: &str, n: i64) -> &str {
    let len = s.len() as i64;
    let n = n.clamp(0, len) as usize;
    &s[s.len() - n..]
}

/// The byte at `pos` in `s`, or `None` if `pos` is out of bounds (including
/// negative).
fn byte_at(s: &str, pos: i64) -> Option<u8> {
    usize::try_from(pos).ok().and_then(|pos| s.as_bytes().get(pos)).copied()
}

/// `end - start + 1` of the final exon pair in a flattened
/// `[start0, end0, start1, end1, ...]` coordinate list. Mirrors
/// `GetLastExonLength` in `ParseDatFile.pl`.
fn last_exon_length(exons: &[i64]) -> i64 {
    exons[exons.len() - 1] - exons[exons.len() - 2] + 1
}

// ---------------------------------------------------------------------------
// Per-allele build state, threaded through the whole pipeline in `.dat` file
// order (mirrors the package-level `%allele*`/`%gene*` hashes in the Perl).
// ---------------------------------------------------------------------------

/// Mutable state accumulated while building one [`SeqKind`]'s output,
/// mirroring `ParseDatFile.pl`'s package-level hashes (`%alleleSeq`,
/// `%alleleExonRegions`, `%allelePaddingLength`, `%gene5UTRPadding`, ...).
#[derive(Default)]
struct BuildState {
    /// Built output sequence per allele (`%alleleSeq`).
    seq: HashMap<String, String>,
    /// Flattened exon-boundary positions in the *output* sequence, per allele
    /// (`%alleleExonRegions`).
    exon_regions: HashMap<String, Vec<i64>>,
    /// Flattened *true* (pre-padding, pre-merge) exon coordinates, per allele;
    /// only populated for [`SeqKind::Dna`] (`%alleleTrueExonRegions`).
    true_exon_regions: HashMap<String, Vec<i64>>,
    /// Coding + UTR length used for gene-length-mode voting
    /// (`%alleleEffectiveLength`).
    effective_length: HashMap<String, i64>,
    /// `(five_prime_pad, three_prime_pad)` bases this allele still owes to its
    /// gene's shared UTR cache (`%allelePaddingLength`).
    padding: HashMap<String, (i64, i64)>,
    /// Alleles not flagged partial by the time their record finished
    /// building, in `.dat` file order (`@alleleOrder`, before rescue).
    order: Vec<String>,
    /// Alleles flagged partial for any reason (explicit `/partial`
    /// annotation, dna-mode touching-exon detection, or
    /// exon-runs-past-sequence-end), in first-flagged order.
    ///
    /// `ParseDatFile.pl`'s `%partialAlleles` is a plain Perl hash, so
    /// `foreach my $allele (keys %partialAlleles)` (used by the
    /// partial-allele-rescue step) iterates in Perl's randomized hash order.
    /// That only affects output *when two or more partial alleles of the same
    /// gene are simultaneously eligible for rescue* (see
    /// `fixtures/refbuild/PINS.md`). This port always imposes a fixed,
    /// deterministic order instead — first-flagged (i.e. `.dat` file) order,
    /// the same convention `@alleleOrder` already uses for non-partial
    /// alleles — rather than replaying any particular Perl hash-seed
    /// behavior. See the `rescue_order_is_deterministic_on_a_tied_input` unit
    /// test below.
    partial_order: Vec<String>,
    partial_set: HashSet<String>,
    /// First-seen 50bp of real sequence immediately before each gene's first
    /// exon (`%gene5UTRPadding`), and immediately after its last exon
    /// (`%gene3UTRPadding`).
    gene_5utr: HashMap<String, String>,
    gene_3utr: HashMap<String, String>,
    /// Distribution of true (genomic, pre-merge) last-exon lengths per gene,
    /// collected from *every* processed record regardless of partial status
    /// (`%geneLastExonLengthDist`).
    gene_last_exon_len_dist: HashMap<String, HashMap<i64, i64>>,
}

impl BuildState {
    /// Flags `allele` as partial, recording the first time this happens per
    /// allele (an ordered analogue of `$partialAlleles{$allele} = 1`).
    fn mark_partial(&mut self, allele: &str) {
        if self.partial_set.insert(allele.to_string()) {
            self.partial_order.push(allele.to_string());
        }
    }

    fn is_partial(&self, allele: &str) -> bool {
        self.partial_set.contains(allele)
    }
}

/// Builds one allele's output sequence and exon-region coordinates for `kind`,
/// updating `state` in place. Mirrors the body of `while (<FP>) { ... }`'s
/// `/^SQ/` branch in `ParseDatFile.pl` (lines building `$outputSeq` through the
/// final `push @alleleOrder, $allele` decision) for a single record.
fn build_one_allele(raw: &RawAllele, kind: SeqKind, state: &mut BuildState) {
    // `last if (scalar(@exons) == 0)`: `parse_dat` already drops such records.
    if raw.exons.is_empty() {
        return;
    }
    let gene = raw.gene().to_string();

    if raw.annotated_partial {
        state.mark_partial(&raw.name);
    }

    let exons: Vec<i64> = raw.exons.iter().flat_map(|&(s, e)| [s, e]).collect();
    let seq_len = raw.sequence.len() as i64;

    // --- 5' UTR before the first exon ---
    let mut padding = (0i64, 0i64);
    let mut prefix_start = exons[0] - UTR_LENGTH;
    let prefix_end = exons[0] - 1;
    if prefix_start < 0 {
        padding.0 = -prefix_start;
        prefix_start = 0;
    } else {
        state.gene_5utr.entry(gene.clone()).or_insert_with(|| {
            perl_substr(&raw.sequence, prefix_start, prefix_end - prefix_start + 1)
                .to_ascii_uppercase()
        });
    }
    let mut output = perl_substr(&raw.sequence, prefix_start, prefix_end - prefix_start + 1)
        .to_ascii_uppercase();

    // --- mode-specific exon/intron assembly ---
    let mut exon_offset = UTR_LENGTH;
    let mut exon_regions = Vec::new();
    match kind {
        SeqKind::Rna => {
            build_rna_body(&exons, &raw.sequence, &mut output, &mut exon_offset, &mut exon_regions);
        }
        SeqKind::Dna => {
            // Adjacent, touching exons (no real intron between them at all)
            // are always flagged partial in dna mode, independent of any
            // `/partial` annotation. Mirrors the Perl's own separate
            // pre-pass over `@exons` before it builds `$outputSeq`.
            let mut i = 2usize;
            while i < exons.len() {
                if exons[i] <= exons[i - 1] + 1 {
                    state.mark_partial(&raw.name);
                }
                i += 2;
            }
            build_dna_body(
                &exons,
                &raw.sequence,
                seq_len,
                &mut output,
                &mut exon_offset,
                &mut exon_regions,
            );
            state.true_exon_regions.insert(raw.name.clone(), exons.clone());
        }
    }

    *state
        .gene_last_exon_len_dist
        .entry(gene.clone())
        .or_default()
        .entry(last_exon_length(&exons))
        .or_insert(0) += 1;

    // --- 3' UTR after the last exon ---
    let after_start = exons[exons.len() - 1] + 1;
    if after_start > seq_len {
        // A partial allele whose annotated exon runs past the end of the
        // available sequence: nothing more is appended to `output`.
        state.mark_partial(&raw.name);
    } else {
        let mut after_end = after_start + UTR_LENGTH - 1;
        if after_end >= seq_len {
            padding.1 = after_end - seq_len + 1;
            after_end = seq_len - 1;
        } else {
            state.gene_3utr.entry(gene.clone()).or_insert_with(|| {
                perl_substr(&raw.sequence, after_start, after_end - after_start + 1)
                    .to_ascii_uppercase()
            });
        }
        output.push_str(
            &perl_substr(&raw.sequence, after_start, after_end - after_start + 1)
                .to_ascii_uppercase(),
        );
    }

    if !state.is_partial(&raw.name) {
        state.order.push(raw.name.clone());
    }

    let effective_length =
        2 * UTR_LENGTH + exons.chunks_exact(2).map(|pair| pair[1] - pair[0] + 1).sum::<i64>();

    state.padding.insert(raw.name.clone(), padding);
    state.seq.insert(raw.name.clone(), output);
    state.exon_regions.insert(raw.name.clone(), exon_regions);
    state.effective_length.insert(raw.name.clone(), effective_length);
}

/// Builds the `rna`-mode exon/intron assembly: append each exon's sequence
/// directly (introns dropped entirely). Mirrors the `if ($mode eq "rna")`
/// branch of `ParseDatFile.pl`'s per-record output builder.
fn build_rna_body(
    exons: &[i64],
    sequence: &str,
    output: &mut String,
    exon_offset: &mut i64,
    exon_regions: &mut Vec<i64>,
) {
    for pair in exons.chunks_exact(2) {
        let (start, end) = (pair[0], pair[1]);
        output.push_str(&perl_substr(sequence, start, end - start + 1).to_ascii_uppercase());
        exon_regions.push(*exon_offset);
        exon_regions.push(*exon_offset + end - start);
        *exon_offset += end - start + 1;
    }
}

/// Builds the `dna`-mode exon/intron assembly: introns are kept, padded to
/// [`INTRON_PADDING_LENGTH`] bases on each side, with runs of exons separated
/// only by a short intron (one whose two padded windows would overlap) merged
/// into a single continuous block (full intron kept, no `N` separator).
/// Mirrors the `elsif ($mode eq "dna")` branch of `ParseDatFile.pl`'s
/// per-record output builder (the intron-merge `for`/nested-`while` loop that
/// manually advances its outer index).
fn build_dna_body(
    exons: &[i64],
    sequence: &str,
    seq_len: i64,
    output: &mut String,
    exon_offset: &mut i64,
    exon_regions: &mut Vec<i64>,
) {
    let mut i = 0usize;
    while i < exons.len() {
        let block_start_idx = i;
        let mut block_start = exons[i];
        if i > 0 {
            block_start = (exons[i] - INTRON_PADDING_LENGTH).max(0);
            *exon_offset += 1 + INTRON_PADDING_LENGTH;
            output.push('N');
        }

        exon_regions.push(*exon_offset);
        exon_regions.push(*exon_offset + exons[i + 1] - exons[i]);

        let mut block_end = exons[i + 1];
        while i + 2 < exons.len() {
            block_end = (exons[i + 1] + INTRON_PADDING_LENGTH).min(seq_len - 1);
            if block_end >= exons[i + 2] - INTRON_PADDING_LENGTH {
                // Short intron: merge the next exon into this same block
                // instead of starting a new, separately-padded one.
                i += 2;
                block_end = exons[i + 1];
                exon_regions.push(*exon_offset + exons[i] - exons[block_start_idx]);
                exon_regions.push(*exon_offset + exons[i + 1] - exons[block_start_idx]);
            } else {
                break;
            }
        }

        output.push_str(
            &perl_substr(sequence, block_start, block_end - block_start + 1).to_ascii_uppercase(),
        );
        *exon_offset += exons[i + 1] - exons[block_start_idx] + 1;
        *exon_offset += INTRON_PADDING_LENGTH;
        i += 2;
    }
}

/// Top-level orchestrator: parses `alleles` (already produced by
/// [`parse_dat`]) into one [`SeqKind`]'s seq-FASTA text, matching
/// `ParseDatFile.pl`'s full pipeline: per-record build, gene shape statistics
/// (dna mode), partial-allele rescue, UTR-padding application, the dna-only
/// "exonization" length fix, the final `fixGeneLength` trim, and FASTA
/// rendering.
///
/// # Errors
///
/// Returns an error if some allele needs synthetic UTR padding but no allele
/// of its gene ever had real flanking sequence to source it from — see the
/// module docs on why this port does not implement `ParseDatFile.pl`'s
/// RNG-based fallback for that case.
pub fn emit_seq_fasta(alleles: &[RawAllele], kind: SeqKind) -> Result<String> {
    let mut state = BuildState::default();
    for raw in alleles {
        build_one_allele(raw, kind, &mut state);
    }

    let mut gene_exon_cnt_mode: HashMap<String, i64> = HashMap::new();
    let mut gene_exon_length_mode: HashMap<String, HashMap<i64, i64>> = HashMap::new();
    let mut gene_true_intron_length_mode: HashMap<String, HashMap<i64, i64>> = HashMap::new();
    let mut gene_length_mode: HashMap<String, i64> = HashMap::new();

    if kind == SeqKind::Dna {
        compute_dna_shape_stats(
            &state,
            &mut gene_exon_cnt_mode,
            &mut gene_exon_length_mode,
            &mut gene_true_intron_length_mode,
            &mut gene_length_mode,
        );
    }

    // `if (scalar(keys %geneLengthMode) == 0) { ... }`: always true for rna
    // mode (the stats block above only runs for dna mode).
    if gene_length_mode.is_empty() {
        gene_length_mode = effective_length_mode(&state);
    }

    rescue_partial_alleles(&mut state, kind, &gene_length_mode, &gene_exon_cnt_mode);

    apply_utr_padding(&mut state)?;

    if kind == SeqKind::Dna {
        fix_exonization(
            &mut state,
            &gene_exon_cnt_mode,
            &gene_exon_length_mode,
            &gene_true_intron_length_mode,
        );
    }

    trim_to_gene_length(&mut state);

    Ok(render_fasta(&state))
}

/// Computes the dna-mode-only gene shape statistics: modal exon count, modal
/// per-exon-position length, modal per-intron-position true length, and modal
/// effective (coding + UTR) length, each per gene. Mirrors
/// `ParseDatFile.pl`'s `if ($mode eq "dna") { ... }` statistics block.
fn compute_dna_shape_stats(
    state: &BuildState,
    gene_exon_cnt_mode: &mut HashMap<String, i64>,
    gene_exon_length_mode: &mut HashMap<String, HashMap<i64, i64>>,
    gene_true_intron_length_mode: &mut HashMap<String, HashMap<i64, i64>>,
    gene_length_mode: &mut HashMap<String, i64>,
) {
    let mut gene_eff_len_dist: HashMap<String, HashMap<i64, i64>> = HashMap::new();
    let mut gene_exon_cnt_dist: HashMap<String, HashMap<i64, i64>> = HashMap::new();
    for allele in &state.order {
        let gene = gene_of(allele);
        *gene_eff_len_dist
            .entry(gene.clone())
            .or_default()
            .entry(state.effective_length[allele])
            .or_insert(0) += 1;
        let exon_cnt = state.exon_regions[allele].len() as i64 / 2;
        *gene_exon_cnt_dist.entry(gene).or_default().entry(exon_cnt).or_insert(0) += 1;
    }
    for (gene, dist) in &gene_eff_len_dist {
        if let Some(mode) = find_mode_i64(dist) {
            gene_length_mode.insert(gene.clone(), mode);
        }
    }
    for (gene, dist) in &gene_exon_cnt_dist {
        if let Some(mode) = find_mode_i64(dist) {
            gene_exon_cnt_mode.insert(gene.clone(), mode);
        }
    }

    // Representative per-exon-position and per-intron-position lengths, using
    // only alleles whose exon count matches their gene's modal count.
    let mut gene_exon_len_dist: HashMap<String, HashMap<i64, HashMap<i64, i64>>> = HashMap::new();
    let mut gene_true_intron_dist: HashMap<String, HashMap<i64, HashMap<i64, i64>>> =
        HashMap::new();
    for allele in &state.order {
        let gene = gene_of(allele);
        let exon_regions = &state.exon_regions[allele];
        let true_exons = &state.true_exon_regions[allele];
        let exon_cnt = exon_regions.len() as i64 / 2;
        if Some(exon_cnt) != gene_exon_cnt_mode.get(&gene).copied() {
            continue;
        }
        for i in 0..exon_cnt {
            let idx = i as usize;
            let length = exon_regions[2 * idx + 1] - exon_regions[2 * idx] + 1;
            *gene_exon_len_dist
                .entry(gene.clone())
                .or_default()
                .entry(i)
                .or_default()
                .entry(length)
                .or_insert(0) += 1;
            if i < exon_cnt - 1 {
                let true_intron_len = true_exons[2 * idx + 2] - true_exons[2 * idx + 1] - 1;
                *gene_true_intron_dist
                    .entry(gene.clone())
                    .or_default()
                    .entry(i)
                    .or_default()
                    .entry(true_intron_len)
                    .or_insert(0) += 1;
            }
        }
    }
    for (gene, per_position) in &gene_exon_len_dist {
        let out = gene_exon_length_mode.entry(gene.clone()).or_default();
        for (position, dist) in per_position {
            if let Some(mode) = find_mode_i64(dist) {
                out.insert(*position, mode);
            }
        }
    }
    for (gene, per_position) in &gene_true_intron_dist {
        let out = gene_true_intron_length_mode.entry(gene.clone()).or_default();
        for (position, dist) in per_position {
            if let Some(mode) = find_mode_i64(dist) {
                out.insert(*position, mode);
            }
        }
    }
}

/// Computes each gene's modal effective length purely from
/// [`BuildState::effective_length`] (no exon-count/shape information needed).
/// Used to seed `gene_length_mode` for rna mode, where the dna-only shape
/// statistics never ran.
fn effective_length_mode(state: &BuildState) -> HashMap<String, i64> {
    let mut dist: HashMap<String, HashMap<i64, i64>> = HashMap::new();
    for allele in &state.order {
        let gene = gene_of(allele);
        *dist.entry(gene).or_default().entry(state.effective_length[allele]).or_insert(0) += 1;
    }
    dist.iter().filter_map(|(gene, d)| find_mode_i64(d).map(|mode| (gene.clone(), mode))).collect()
}

/// Appends alleles from [`BuildState::partial_order`] that meet the rescue
/// criteria back onto [`BuildState::order`] (mutating `state.seq` /
/// `state.exon_regions` in dna mode, to splice in synthetic intron
/// sequence). Mirrors the "Rescue partial alleles" block of
/// `ParseDatFile.pl`, which — given this port's hardcoded default options
/// (see module docs) — always runs (`$includePartialDiffLen >= 0 &&
/// $ignorePartial == 0` is always true).
///
/// Iterates `state.partial_order` (a fixed, deterministic order) rather than
/// Perl's `keys %partialAlleles` (randomized hash order); see
/// [`BuildState::partial_order`]'s docs.
fn rescue_partial_alleles(
    state: &mut BuildState,
    kind: SeqKind,
    gene_length_mode: &HashMap<String, i64>,
    gene_exon_cnt_mode: &HashMap<String, i64>,
) {
    let candidates = state.partial_order.clone();
    let mut rescued = Vec::new();

    match kind {
        SeqKind::Rna => {
            for allele in &candidates {
                let gene = gene_of(allele);
                let Some(&mode_len) = gene_length_mode.get(&gene) else { continue };
                if state.effective_length[allele] >= mode_len - INCLUDE_PARTIAL_DIFF_LEN {
                    rescued.push(allele.clone());
                }
            }
        }
        SeqKind::Dna => {
            let gene_intron_mode = build_intron_consensus(state, gene_exon_cnt_mode);
            for allele in &candidates {
                let gene = gene_of(allele);
                let Some(&mode_len) = gene_length_mode.get(&gene) else { continue };
                if state.effective_length[allele] < mode_len - INCLUDE_PARTIAL_DIFF_LEN {
                    continue;
                }
                let exon_cnt = state.exon_regions[allele].len() as i64 / 2;
                if Some(exon_cnt) != gene_exon_cnt_mode.get(&gene).copied() {
                    continue;
                }
                splice_synthetic_introns(state, allele, &gene, &gene_intron_mode);
                rescued.push(allele.clone());
            }
        }
    }

    state.order.extend(rescued);
}

/// Builds each gene's per-intron-position consensus sequence
/// (`%geneIntronSeqMode`), used to synthesize intron sequence for dna-mode
/// rescue. Only alleles already in `state.order` (i.e. non-partial,
/// pre-rescue) with an exon count matching their gene's mode contribute.
fn build_intron_consensus(
    state: &BuildState,
    gene_exon_cnt_mode: &HashMap<String, i64>,
) -> HashMap<String, HashMap<i64, String>> {
    let mut dist: HashMap<String, HashMap<i64, HashMap<String, i64>>> = HashMap::new();
    for allele in &state.order {
        let gene = gene_of(allele);
        let exon_regions = &state.exon_regions[allele];
        let exon_cnt = exon_regions.len() as i64 / 2;
        if Some(exon_cnt) != gene_exon_cnt_mode.get(&gene).copied() {
            continue;
        }
        let seq = &state.seq[allele];
        let mut i = 2usize;
        while (i as i64) < 2 * exon_cnt {
            let intron_seq = perl_substr(
                seq,
                exon_regions[i - 1] + 1,
                exon_regions[i] - exon_regions[i - 1] - 1,
            )
            .to_string();
            let position = i as i64 / 2 - 1;
            *dist
                .entry(gene.clone())
                .or_default()
                .entry(position)
                .or_default()
                .entry(intron_seq)
                .or_insert(0) += 1;
            i += 2;
        }
    }

    let mut mode: HashMap<String, HashMap<i64, String>> = HashMap::new();
    for (gene, per_position) in &dist {
        let out = mode.entry(gene.clone()).or_default();
        for (position, seq_dist) in per_position {
            if let Some(consensus) = find_mode_str(seq_dist) {
                out.insert(*position, consensus);
            }
        }
    }
    mode
}

/// Splices synthetic intron sequence (borrowed from `gene_intron_mode`) into
/// `allele`'s output at every exon boundary that is exactly touching (no gap
/// at all — the dna-mode signature of a record with no real intron
/// annotations), then rebases the exon-region coordinates for the newly
/// inserted bases. Mirrors the "Adding introns to the partial alleles" loop
/// of `ParseDatFile.pl`'s dna-mode rescue block.
fn splice_synthetic_introns(
    state: &mut BuildState,
    allele: &str,
    gene: &str,
    gene_intron_mode: &HashMap<String, HashMap<i64, String>>,
) {
    let mut exons = state.exon_regions[allele].clone();
    let exon_cnt = exons.len() as i64 / 2;
    let extra_5utr = state.padding[allele].0;

    // Exon-region coordinates assume the full 50bp 5' UTR is already present,
    // but UTR padding is applied later (see `apply_utr_padding`), so rebase
    // to the *current* (possibly still-short) output sequence first.
    for exon in &mut exons {
        *exon -= extra_5utr;
    }

    let mut output = state.seq[allele].clone();
    let mut exon_offset = 0i64;
    let mut i = 2usize;
    while (i as i64) < 2 * exon_cnt {
        if exons[i] + exon_offset == exons[i - 1] + 1 {
            let position = i as i64 / 2 - 1;
            if let Some(intron_seq) = gene_intron_mode.get(gene).and_then(|m| m.get(&position)) {
                let insert_at = (exons[i - 1] + 1) as usize;
                output.insert_str(insert_at, intron_seq);
                exon_offset += intron_seq.len() as i64;
            }
        }
        exons[i] += exon_offset;
        exons[i + 1] += exon_offset;
        i += 2;
    }

    for exon in &mut exons {
        *exon += extra_5utr;
    }

    state.exon_regions.insert(allele.to_string(), exons);
    state.seq.insert(allele.to_string(), output);
}

/// Prepends/appends each allele's owed 5'/3' UTR padding, borrowed from its
/// gene's cached flanking sequence. Mirrors `ParseDatFile.pl`'s final UTR-
/// application loop (the `srand(17)`-seeded random-fallback branches are not
/// implemented — see module docs — so this errors instead of silently
/// diverging from the Perl when a gene's cache is empty).
///
/// # Errors
///
/// See [`emit_seq_fasta`].
fn apply_utr_padding(state: &mut BuildState) -> Result<()> {
    for allele in state.order.clone() {
        let gene = gene_of(&allele);
        let (pad_5p, pad_3p) = state.padding[&allele];
        let mut seq = state.seq[&allele].clone();

        if pad_5p > 0 {
            let source = state
                .gene_5utr
                .get(&gene)
                .with_context(|| missing_utr_padding_message(&allele, &gene, "5'", pad_5p))?;
            seq = format!("{}{seq}", perl_substr(source, 0, pad_5p));
        }
        if pad_3p > 0 {
            let source = state
                .gene_3utr
                .get(&gene)
                .with_context(|| missing_utr_padding_message(&allele, &gene, "3'", pad_3p))?;
            seq.push_str(perl_suffix(source, pad_3p));
        }

        state.seq.insert(allele, seq);
    }
    Ok(())
}

fn missing_utr_padding_message(allele: &str, gene: &str, side: &str, bases: i64) -> String {
    format!(
        "allele {allele} (gene {gene}) needs {bases} bases of synthetic {side} UTR padding, but no \
         allele of this gene had real flanking sequence to borrow from; this port does not implement \
         ParseDatFile.pl's srand(17)-seeded random-UTR fallback for that case (see refbuild::dat module docs)"
    )
}

/// Trims mis-annotated "exonized" intron sequence from an exon whose *output*
/// length exceeds its gene's modal length for that exon position, when the
/// excess exactly matches the gap to the intron-length mode on one side and
/// that side is bounded by an `N` padding separator. Dna-mode only. Mirrors
/// `ParseDatFile.pl`'s "Fix the exonization" block.
///
/// Not exercised by the pinned `kir_subset.dat` golden fixture (every grouped
/// allele there already agrees exactly on exon-region shape), so this is unit
/// tested directly against synthetic `BuildState` data instead — see
/// `exonization_fix_trims_a_right_side_overhang` below.
fn fix_exonization(
    state: &mut BuildState,
    gene_exon_cnt_mode: &HashMap<String, i64>,
    gene_exon_length_mode: &HashMap<String, HashMap<i64, i64>>,
    gene_true_intron_length_mode: &HashMap<String, HashMap<i64, i64>>,
) {
    for allele in state.order.clone() {
        let gene = gene_of(&allele);
        let mut exons = state.exon_regions[&allele].clone();
        let exon_cnt = exons.len() as i64 / 2;
        if Some(exon_cnt) != gene_exon_cnt_mode.get(&gene).copied() {
            continue;
        }
        let Some(true_exons) = state.true_exon_regions.get(&allele).cloned() else { continue };
        let Some(exon_len_mode) = gene_exon_length_mode.get(&gene) else { continue };
        let Some(true_intron_mode) = gene_true_intron_length_mode.get(&gene) else { continue };

        let mut seq = state.seq[&allele].clone();
        let mut updated = false;

        for i in 0..(exon_cnt - 1) {
            let idx = i as usize;
            let exon_length = exons[2 * idx + 1] - exons[2 * idx] + 1;
            let Some(&mode_length) = exon_len_mode.get(&i) else { continue };
            if exon_length <= mode_length {
                continue;
            }
            let trim = exon_length - mode_length;

            let Some((trim_side, new_seq)) = exonization_trim_candidate(
                &seq,
                &exons,
                &true_exons,
                true_intron_mode,
                idx,
                i,
                trim,
            ) else {
                continue;
            };

            seq = new_seq;
            if trim > INTRON_PADDING_LENGTH {
                if trim_side == 1 {
                    exons[2 * idx + 1] -= trim - INTRON_PADDING_LENGTH;
                } else {
                    exons[2 * idx] += trim + INTRON_PADDING_LENGTH;
                }
            }
            if trim_side == -1 {
                exons[2 * idx] -= trim;
                exons[2 * idx + 1] -= trim;
            }
            for exon in exons.iter_mut().skip(2 * idx + 2) {
                *exon -= trim;
            }
            updated = true;
        }

        state.seq.insert(allele.clone(), seq);
        if updated {
            state.exon_regions.insert(allele, exons);
        }
    }
}

/// Decides whether exon `i` (flat index `idx = 2*i`) has an exonization
/// overhang trimmable on its right side (into the following intron) or left
/// side (into the preceding intron), and if so returns
/// `(trim_side, new_seq)` — `trim_side` is `1` for right, `-1` for left,
/// matching `ParseDatFile.pl`'s `$trimSide`. `ParseDatFile.pl`'s `$posN` is
/// only ever used to compute `$newSeq` (never read again afterwards), so it
/// is not part of this function's return value.
fn exonization_trim_candidate(
    seq: &str,
    exons: &[i64],
    true_exons: &[i64],
    true_intron_mode: &HashMap<i64, i64>,
    idx: usize,
    i: i64,
    trim: i64,
) -> Option<(i64, String)> {
    let right_true_intron = true_exons[2 * idx + 2] - true_exons[2 * idx + 1] - 1;
    let right_pos_n = exons[2 * idx + 1] + 1 + INTRON_PADDING_LENGTH;
    if true_intron_mode.get(&i).copied() == Some(right_true_intron + trim)
        && right_pos_n < seq.len() as i64
        && byte_at(seq, right_pos_n) == Some(b'N')
    {
        let new_seq =
            format!("{}{}", &seq[..(right_pos_n - trim) as usize], &seq[right_pos_n as usize..]);
        return Some((1, new_seq));
    }

    if i > 0 {
        let left_true_intron = true_exons[2 * idx] - true_exons[2 * idx - 1] - 1;
        let left_pos_n = exons[2 * idx] - 1 - INTRON_PADDING_LENGTH;
        if true_intron_mode.get(&(i - 1)).copied() == Some(left_true_intron + trim)
            && left_pos_n >= 0
            && byte_at(seq, exons[2 * idx - 1] - 1 - INTRON_PADDING_LENGTH) == Some(b'N')
        {
            let new_seq = format!(
                "{}{}",
                &seq[..(left_pos_n + 1) as usize],
                &seq[(left_pos_n + trim + 1) as usize..]
            );
            return Some((-1, new_seq));
        }
    }

    None
}

/// Trims each allele's final output sequence if it is longer than its gene's
/// modal output length *and* its last exon is longer than the gene's modal
/// last-exon length (using the *true*, pre-merge last-exon length collected
/// while building, from every processed record — not just alleles surviving
/// to `state.order`). Runs for both modes (`ParseDatFile.pl`'s
/// `$fixGeneLength` is `1` for both `dna` and `rna`; only `--mode genome`,
/// out of scope here, sets it to `0`).
fn trim_to_gene_length(state: &mut BuildState) {
    let mut gene_seq_len_dist: HashMap<String, HashMap<i64, i64>> = HashMap::new();
    for allele in &state.order {
        let gene = gene_of(allele);
        let length = state.seq[allele].len() as i64;
        *gene_seq_len_dist.entry(gene).or_default().entry(length).or_insert(0) += 1;
    }
    let gene_seq_length: HashMap<String, i64> = gene_seq_len_dist
        .iter()
        .filter_map(|(gene, dist)| find_mode_i64(dist).map(|mode| (gene.clone(), mode)))
        .collect();
    let gene_last_exon_length: HashMap<String, i64> = state
        .gene_last_exon_len_dist
        .iter()
        .filter_map(|(gene, dist)| find_mode_i64(dist).map(|mode| (gene.clone(), mode)))
        .collect();

    for allele in state.order.clone() {
        let gene = gene_of(&allele);
        let Some(&mode_seq_length) = gene_seq_length.get(&gene) else { continue };
        let Some(&mode_last_exon_length) = gene_last_exon_length.get(&gene) else { continue };
        let last_exon_len = last_exon_length(&state.exon_regions[&allele]);
        let trim = last_exon_len - mode_last_exon_length;

        let seq = state.seq.get_mut(&allele).expect("built earlier in the pipeline");
        if seq.len() as i64 > mode_seq_length && trim > 0 {
            let new_len = (seq.len() as i64 - trim).max(0) as usize;
            seq.truncate(new_len);
        }
    }
}

/// Renders `state.order`'s final sequences as seq-FASTA text: one
/// `>{allele} {exon_count} {region1} {region2} ...` header line followed by
/// one unwrapped sequence line, per allele. Mirrors `ParseDatFile.pl`'s final
/// `print` loop (with `--dedup` and `--gene` filtering omitted, since
/// `t1k-build.pl` never passes them — see module docs).
fn render_fasta(state: &BuildState) -> String {
    let mut out = String::new();
    for allele in &state.order {
        let seq = &state.seq[allele];
        if seq.is_empty() {
            continue;
        }
        let regions = &state.exon_regions[allele];
        let region_list = regions.iter().map(i64::to_string).collect::<Vec<_>>().join(" ");
        let _ = writeln!(out, ">{allele} {} {region_list}", regions.len() / 2);
        let _ = writeln!(out, "{seq}");
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------
    // find_mode: the string-based tie-break quirk
    // -----------------------------------------------------------------

    #[test]
    fn find_mode_i64_breaks_ties_lexicographically_not_numerically() {
        // Both keys tied at count 2: Perl's `$k ge $ret` compares the keys'
        // *string* form, so "9" (first byte '9') beats "10" (first byte
        // '1') even though 10 > 9 numerically.
        let mut counts = HashMap::new();
        counts.insert(9i64, 2i64);
        counts.insert(10i64, 2i64);
        assert_eq!(find_mode_i64(&counts), Some(9));
    }

    #[test]
    fn find_mode_i64_picks_the_strict_majority_regardless_of_key_text() {
        let mut counts = HashMap::new();
        counts.insert(10i64, 5i64);
        counts.insert(9i64, 1i64);
        assert_eq!(find_mode_i64(&counts), Some(10));
    }

    #[test]
    fn find_mode_str_breaks_ties_lexicographically() {
        let mut counts = HashMap::new();
        counts.insert("AAAA".to_string(), 3i64);
        counts.insert("TTTT".to_string(), 3i64);
        assert_eq!(find_mode_str(&counts), Some("TTTT".to_string()));
    }

    // -----------------------------------------------------------------
    // Deterministic partial-allele-rescue ordering (multi-way tie)
    // -----------------------------------------------------------------

    /// Builds a minimal in-memory `RawAllele` for a gene with one full-length
    /// "anchor" allele (8 exons, ample flanking sequence, providing the
    /// gene's UTR padding and effective-length mode) plus a set of partial
    /// alleles that are all simultaneously eligible for rna-mode rescue
    /// (same gene, same effective length, all `annotated_partial`). This is
    /// exactly the "multi-way tie" scenario `fixtures/refbuild/PINS.md`
    /// flagged as needing a designed (not Perl-hash-order-replayed)
    /// deterministic order.
    fn tied_rescue_fixture() -> Vec<RawAllele> {
        let exon_len = 30i64;
        let flank = "A".repeat(200);
        let make = |name: &str, partial: bool| {
            let mut seq = flank.clone();
            let exon_start = seq.len() as i64;
            seq.push_str(&"C".repeat(exon_len as usize));
            seq.push_str(&flank);
            RawAllele {
                name: name.to_string(),
                exons: vec![(exon_start, exon_start + exon_len - 1)],
                annotated_partial: partial,
                sequence: seq,
            }
        };
        vec![
            make("GENE*anchor", false),
            make("GENE*partialA", true),
            make("GENE*partialB", true),
            make("GENE*partialC", true),
        ]
    }

    #[test]
    fn rescue_order_is_deterministic_on_a_tied_input() {
        let alleles = tied_rescue_fixture();
        let first = emit_seq_fasta(&alleles, SeqKind::Rna).expect("first run");
        for _ in 0..20 {
            let again = emit_seq_fasta(&alleles, SeqKind::Rna).expect("repeat run");
            assert_eq!(first, again, "rescue order must be stable across repeated runs");
        }
    }

    #[test]
    fn rescue_order_matches_dat_file_encounter_order() {
        // All three partial alleles tie exactly on every rescue criterion
        // (same gene, same effective length); the deterministic tie-break
        // this port imposes is "first-encountered-in-the-.dat-file wins",
        // so they must appear after the anchor in that same order.
        let alleles = tied_rescue_fixture();
        let rna = emit_seq_fasta(&alleles, SeqKind::Rna).expect("emit rna");
        let headers: Vec<&str> = rna.lines().filter(|l| l.starts_with('>')).collect();
        let names: Vec<&str> =
            headers.iter().map(|h| h.trim_start_matches('>').split(' ').next().unwrap()).collect();
        assert_eq!(names, vec!["GENE*anchor", "GENE*partialA", "GENE*partialB", "GENE*partialC"]);
    }

    // -----------------------------------------------------------------
    // Parsing edge case: a `/pseudo`-qualified exon is folded back out
    // -----------------------------------------------------------------

    #[test]
    fn pseudo_qualified_exon_is_removed_from_the_exon_list() {
        let dat = "\
ID   TEST00001; SV 1; standard; DNA; HUM; 100 BP.
FT   source          1..100
FT   CDS             join(1..10,21..30,41..50)
FT                   /allele=\"TESTGENE*0010101\"
FT   exon            1..10
FT                   /number=\"1\"
FT   exon            21..30
FT                   /number=\"2\"
FT                   /pseudo
FT   exon            41..50
FT                   /number=\"3\"
SQ   Sequence 100 BP;
     aaaaaaaaaa cccccccccc gggggggggg aaaaaaaaaa cccccccccc gggggggggg
     aaaaaaaaaa cccccccccc gggggggggg aaaaaaaaaa                          100
//
";
        let alleles = parse_dat_reader(dat.as_bytes()).expect("parsing synthetic record");
        assert_eq!(alleles.len(), 1);
        let allele = &alleles[0];
        assert_eq!(allele.name, "TESTGENE*0010101");
        // The middle exon (21..30) was popped by the `/pseudo` qualifier, so
        // only the first and third annotated exons survive, 0-based.
        assert_eq!(allele.exons, vec![(0, 9), (40, 49)]);
        assert!(!allele.annotated_partial);
    }

    // -----------------------------------------------------------------
    // Exonization fix (dna mode only; not exercised by the golden fixture)
    // -----------------------------------------------------------------

    #[test]
    fn exonization_fix_trims_a_right_side_overhang() {
        // Two alleles of the same synthetic gene, both dna-mode, 2 exons
        // each. The first exon's *output* length disagrees between them by
        // 5 bases in a way that exactly matches "the boundary was drawn 5
        // bases into what should have been intron" (an over-called exon
        // eating into the adjacent N-padded intron on its right/3' side):
        // `allele_b`'s first exon is 5 bases longer, and its true intron is
        // exactly 5 bases shorter than `allele_a`'s, so trimming 5 bases off
        // the overhang and shifting the intron boundary reconciles them.
        let gene_exon_cnt_mode: HashMap<String, i64> = [("GENE".to_string(), 2i64)].into();

        let mut gene_exon_length_mode: HashMap<String, HashMap<i64, i64>> = HashMap::new();
        gene_exon_length_mode.entry("GENE".to_string()).or_default().insert(0, 20); // exon 0 mode length: 20

        let mut gene_true_intron_length_mode: HashMap<String, HashMap<i64, i64>> = HashMap::new();
        gene_true_intron_length_mode.entry("GENE".to_string()).or_default().insert(0, 300); // intron 0 mode length: 300

        let mut state = BuildState {
            order: vec!["GENE*a".to_string(), "GENE*b".to_string()],
            ..BuildState::default()
        };

        // allele_a: "clean" shape, exon 0 = 20 bases, true intron 0 = 300 bases.
        state.exon_regions.insert("GENE*a".to_string(), vec![50, 69, 570, 599]);
        state.true_exon_regions.insert("GENE*a".to_string(), vec![50, 69, 370, 399]);
        {
            let mut seq = "X".repeat(50);
            seq.push_str(&"E".repeat(20)); // exon 0 (positions 50..=69)
            seq.push_str(&"I".repeat(INTRON_PADDING_LENGTH as usize)); // padded intron
            seq.push('N'); // separator, at position 69 + 1 + INTRON_PADDING_LENGTH
            seq.push_str(&"E".repeat(30)); // exon 1 placeholder
            state.seq.insert("GENE*a".to_string(), seq);
        }

        // allele_b: exon 0 over-called by 5 bases (25 instead of 20); its
        // *true* intron 0 is 5 bases shorter (295) than allele_a's, exactly
        // accounting for the 5-base overhang. The byte right after
        // `exons[1] + 1 + INTRON_PADDING_LENGTH` must be 'N' for the
        // right-side trim branch to trigger.
        let exon0_len = 25i64;
        let exon0_end = 50 + exon0_len - 1; // 74
        let n_pos = exon0_end + 1 + INTRON_PADDING_LENGTH;
        state.exon_regions.insert("GENE*b".to_string(), vec![50, exon0_end, n_pos + 1, n_pos + 31]);
        state.true_exon_regions.insert(
            "GENE*b".to_string(),
            vec![50, exon0_end, 50 + exon0_len + 295, 50 + exon0_len + 295 + 30],
        );
        {
            let mut seq = "X".repeat(50);
            seq.push_str(&"E".repeat(exon0_len as usize)); // over-called exon 0
            seq.push_str(&"I".repeat(INTRON_PADDING_LENGTH as usize));
            seq.push('N'); // separator, at position `n_pos`
            seq.push_str(&"E".repeat(30));
            state.seq.insert("GENE*b".to_string(), seq);
        }

        let allele_a_seq_before = state.seq["GENE*a"].clone();
        let allele_b_len_before = state.seq["GENE*b"].len();
        fix_exonization(
            &mut state,
            &gene_exon_cnt_mode,
            &gene_exon_length_mode,
            &gene_true_intron_length_mode,
        );

        assert_eq!(
            state.seq["GENE*a"], allele_a_seq_before,
            "allele_a has no overhang and must be untouched"
        );
        assert_eq!(
            state.seq["GENE*b"].len(),
            allele_b_len_before - 5,
            "5-base overhang must be trimmed"
        );
        // allele_a's exon regions must be untouched (no overhang there).
        assert_eq!(state.exon_regions["GENE*a"], vec![50, 69, 570, 599]);
        // Since `trim` (5) does not exceed `INTRON_PADDING_LENGTH` (200), the
        // removed bytes come entirely from the intron padding *after* exon 0's
        // already-recorded end, not from the exon boundary itself — so
        // `ParseDatFile.pl` leaves exon 0's own region untouched and only
        // shifts every *later* exon left by `trim`.
        assert_eq!(state.exon_regions["GENE*b"][0], 50);
        assert_eq!(state.exon_regions["GENE*b"][1], exon0_end);
        assert_eq!(state.exon_regions["GENE*b"][2], n_pos + 1 - 5);
        assert_eq!(state.exon_regions["GENE*b"][3], n_pos + 31 - 5);
    }
}
