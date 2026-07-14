//! Combine genotype results from multiple samples into one allele-by-sample
//! abundance matrix -- the Rust port of T1K's `t1k-merge.py` (the `unum combine`
//! subcommand). Built for many-sample runs (e.g. smart-seq): given a list of
//! `{prefix}_genotype.tsv` files, it picks a set of representative alleles per
//! gene by quality-weighted voting, then reports each sample's abundance for
//! those alleles plus any calls that fell outside the representative set.
//!
//! # Algorithm (faithful to `t1k-merge.py`)
//!
//! 1. **Vote.** For every sample and every *active* allele slot (slot `k` is
//!    active iff `k < num_alleles`) whose quality is `> min_quality`, add that
//!    quality to a per-gene tally keyed on the slot's *first* allele (only the
//!    first of an equally-best set votes).
//! 2. **Select representatives.** Per gene, keep the top `num_alleles_per_gene`
//!    alleles by tallied quality whose tally is `>= min_total_quality`; ties are
//!    broken by first-voted order (the stable sort in `t1k-merge.py`). Their
//!    union (sorted) is the matrix's columns.
//! 3. **Fill.** For each sample and each active slot with quality `> min_quality`,
//!    if any allele in the slot's set is a representative, add the slot's
//!    abundance to that sample's count for it; otherwise record the slot as an
//!    inconsistency (`allele[_allele...]_abundance_quality`, raw strings).
//!
//! Cell formatting mirrors the Python: an untouched count prints as `0` (the
//! integer initializer), a touched one via Python `str(float)` semantics -- so
//! a whole-number abundance keeps its trailing `.0` (`50.0` -> `"50.0"`), unlike
//! Rust's bare `f64` `Display`. See [`format_count`].

use std::collections::BTreeMap;

/// One called allele slot from a genotype line: the (comma-joined) equally-best
/// allele set plus its abundance and genotype quality. The raw abundance/quality
/// strings are retained so inconsistency records reproduce the source verbatim.
#[derive(Debug, Clone, PartialEq)]
pub struct AlleleCall {
    /// The equally-best alleles for this slot (T1K's comma-separated set).
    pub alleles: Vec<String>,
    /// Parsed abundance (T1K column `i+1`).
    pub abundance: f64,
    /// Parsed genotype quality (T1K column `i+2`).
    pub quality: f64,
    /// Raw abundance string, verbatim from the file (for inconsistency output).
    pub abundance_raw: String,
    /// Raw quality string, verbatim from the file (for inconsistency output).
    pub quality_raw: String,
}

/// One gene's call, parsed from a `{prefix}_genotype.tsv` line: the gene name,
/// the number of alleles called (T1K column 2), and the *active* primary allele
/// slots (at most two, cols 3-5 and 6-8, gated by `num_alleles`).
#[derive(Debug, Clone, PartialEq)]
pub struct GeneCall {
    /// Gene name (column 1).
    pub gene: String,
    /// Number of alleles called (column 2); slots `k < num_alleles` are active.
    pub num_alleles: usize,
    /// Active primary slots, in file order.
    pub slots: Vec<AlleleCall>,
}

/// Configuration mirroring `t1k-merge.py`'s `-n`, `-q`, and `--tq` flags.
#[derive(Debug, Clone, Copy)]
pub struct CombineConfig {
    /// `-n`: representative alleles kept per gene (default 2).
    pub num_alleles_per_gene: usize,
    /// `-q`: a slot votes/contributes only when `quality > min_quality` (default 0).
    pub min_quality: f64,
    /// `--tq`: a candidate is representative only when its tallied quality
    /// `>= min_total_quality` (default 30).
    pub min_total_quality: f64,
}

impl Default for CombineConfig {
    fn default() -> Self {
        Self { num_alleles_per_gene: 2, min_quality: 0.0, min_total_quality: 30.0 }
    }
}

/// One output row: a sample, its abundance for each representative allele (in
/// column order), and the inconsistency records for calls outside the set.
#[derive(Debug, Clone, PartialEq)]
pub struct CombineRow {
    /// Sample name.
    pub sample: String,
    /// `Some(abundance)` if the allele was hit in this sample, else `None`
    /// (printed as the integer `0`, matching the Python initializer). Parallel
    /// to [`CombineMatrix::alleles`].
    pub counts: Vec<Option<f64>>,
    /// Inconsistency records (`allele[_allele...]_abundance_quality`).
    pub inconsistencies: Vec<String>,
}

/// The full allele-by-sample matrix: sorted representative alleles (columns) and
/// one [`CombineRow`] per input sample (in input order).
#[derive(Debug, Clone, PartialEq)]
pub struct CombineMatrix {
    /// Representative alleles, sorted -- the matrix columns.
    pub alleles: Vec<String>,
    /// Per-sample rows, in the order the samples were supplied.
    pub rows: Vec<CombineRow>,
}

/// Parse the active primary slots of one `{prefix}_genotype.tsv` line into a
/// [`GeneCall`], or `None` if the line is blank / lacks the gene+count columns.
///
/// Columns (tab-separated): `gene  num_alleles  a1 abund1 qual1  a2 abund2 qual2 ...`.
/// Slot `k` (at column `2 + 3*k`) is included only when `k < num_alleles` and the
/// three columns are present; the allele field is split on `,` into the set.
#[must_use]
pub fn parse_genotype_line(line: &str) -> Option<GeneCall> {
    let cols: Vec<&str> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
    if cols.len() < 2 {
        return None;
    }
    let gene = cols[0];
    if gene.is_empty() {
        return None;
    }
    let num_alleles: usize = cols[1].trim().parse().ok()?;
    let mut slots = Vec::new();
    for (k, i) in [2usize, 5].into_iter().enumerate() {
        if k >= num_alleles {
            break;
        }
        if cols.len() <= i + 2 {
            break;
        }
        let (abundance, quality) = (cols[i + 1], cols[i + 2]);
        slots.push(AlleleCall {
            alleles: cols[i].split(',').map(str::to_owned).collect(),
            abundance: abundance.parse().unwrap_or(0.0),
            quality: quality.parse().unwrap_or(0.0),
            abundance_raw: abundance.to_owned(),
            quality_raw: quality.to_owned(),
        });
    }
    Some(GeneCall { gene: gene.to_owned(), num_alleles, slots })
}

/// Combine per-sample [`GeneCall`] lists into a [`CombineMatrix`]. `samples` is
/// `(sample_name, gene_calls)` in the desired output row order. Pure: no I/O.
#[must_use]
pub fn combine(samples: &[(String, Vec<GeneCall>)], cfg: &CombineConfig) -> CombineMatrix {
    // Pass 1 -- quality-weighted voting on each slot's FIRST allele. Per gene we
    // keep both the tally (allele -> summed quality) and the order in which
    // alleles were first voted for; the latter is the tie-break key in Pass 2,
    // matching `t1k-merge.py` (see below).
    let mut gene_tally: BTreeMap<&str, BTreeMap<&str, f64>> = BTreeMap::new();
    let mut gene_vote_order: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for (_, calls) in samples {
        for gc in calls {
            for slot in &gc.slots {
                if slot.quality > cfg.min_quality {
                    if let Some(first) = slot.alleles.first() {
                        let tally = gene_tally.entry(&gc.gene).or_default();
                        if !tally.contains_key(first.as_str()) {
                            gene_vote_order.entry(&gc.gene).or_default().push(first.as_str());
                        }
                        *tally.entry(first.as_str()).or_insert(0.0) += slot.quality;
                    }
                }
            }
        }
    }

    // Pass 2 -- top-N per gene above the total-quality floor -> representatives.
    let mut representatives: BTreeMap<String, f64> = BTreeMap::new();
    for (gene, order) in &gene_vote_order {
        let tally = &gene_tally[gene];
        let mut ranked: Vec<(&str, f64)> = order.iter().map(|&a| (a, tally[a])).collect();
        // Highest tally first; ties keep first-voted order. `t1k-merge.py`'s
        // `sorted(..., reverse=True)` is a *stable* sort over the dict's
        // insertion order, so equal tallies retain the order they were first
        // voted for -- `slice::sort_by` is likewise stable and `ranked` is built
        // in first-voted order, so a tie-break-free comparator reproduces it.
        // `total_cmp` gives a panic-free total order (qualities are finite).
        ranked.sort_by(|a, b| b.1.total_cmp(&a.1));
        for (allele, total) in ranked.into_iter().take(cfg.num_alleles_per_gene) {
            if total >= cfg.min_total_quality {
                representatives.insert(allele.to_owned(), total);
            }
        }
    }
    let alleles: Vec<String> = representatives.keys().cloned().collect();
    let col_index: BTreeMap<&str, usize> =
        alleles.iter().enumerate().map(|(i, a)| (a.as_str(), i)).collect();

    // Pass 3 -- per-sample abundances + inconsistency records.
    let mut rows = Vec::with_capacity(samples.len());
    for (name, calls) in samples {
        let mut counts: Vec<Option<f64>> = vec![None; alleles.len()];
        let mut inconsistencies = Vec::new();
        for gc in calls {
            for slot in &gc.slots {
                if slot.quality <= cfg.min_quality {
                    continue;
                }
                if let Some(&idx) = slot.alleles.iter().find_map(|a| col_index.get(a.as_str())) {
                    counts[idx] = Some(counts[idx].unwrap_or(0.0) + slot.abundance);
                } else {
                    let mut parts = slot.alleles.clone();
                    parts.push(slot.abundance_raw.clone());
                    parts.push(slot.quality_raw.clone());
                    inconsistencies.push(parts.join("_"));
                }
            }
        }
        rows.push(CombineRow { sample: name.clone(), counts, inconsistencies });
    }

    CombineMatrix { alleles, rows }
}

/// Format a touched count exactly as Python's `str(float)` renders it, so
/// combine output stays byte-identical to `t1k-merge.py`. Rust's `f64` `Display`
/// drops the trailing `.0` on whole numbers (`50.0` -> `"50"`) whereas Python
/// keeps it (`str(50.0) == "50.0"`); reattach it for finite whole numbers.
/// Fractional values already agree (both use a shortest-round-trip algorithm),
/// as do non-finite values (`inf`/`NaN` never carry a `.0`), so they pass
/// through unchanged.
fn format_count(value: f64) -> String {
    let rendered = format!("{value}");
    if value.is_finite() && !rendered.contains(['.', 'e', 'E']) {
        format!("{rendered}.0")
    } else {
        rendered
    }
}

/// Render a [`CombineMatrix`] as the tab-separated text T1K prints to stdout:
/// header `sample <alleles...> inconsistency`, then one row per sample. An
/// untouched count prints as the integer `0`, a touched one via
/// [`format_count`] (Python `str(float)` semantics, e.g. `50.0` -> `"50.0"`).
#[must_use]
pub fn format_matrix(matrix: &CombineMatrix) -> String {
    let mut out = String::new();
    let mut header = vec!["sample".to_owned()];
    header.extend(matrix.alleles.iter().cloned());
    header.push("inconsistency".to_owned());
    out.push_str(&header.join("\t"));
    out.push('\n');
    for row in &matrix.rows {
        let mut fields = vec![row.sample.clone()];
        for count in &row.counts {
            fields.push(match count {
                Some(v) => format_count(*v),
                None => "0".to_owned(),
            });
        }
        fields.push(row.inconsistencies.join(","));
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(line: &str) -> GeneCall {
        parse_genotype_line(line).expect("parse")
    }

    #[test]
    fn parses_het_hom_and_nocall_lines() {
        let het = call("HLA-A\t2\tHLA-A*24:02\t61.5\t35\tHLA-A*11:01\t57.4\t32");
        assert_eq!(het.num_alleles, 2);
        assert_eq!(het.slots.len(), 2);
        assert_eq!(het.slots[0].alleles, vec!["HLA-A*24:02"]);
        assert!((het.slots[0].abundance - 61.5).abs() < 1e-9);
        assert!((het.slots[1].quality - 32.0).abs() < 1e-9);

        // Homozygous: only the first slot is active (k < num_alleles).
        let hom = call("HLA-B\t1\tHLA-B*40:01\t102.3\t60\t.\t0\t-1");
        assert_eq!(hom.slots.len(), 1);
        assert_eq!(hom.slots[0].alleles, vec!["HLA-B*40:01"]);

        // No-call: zero active slots.
        let nocall = call("HLA-R\t0\t.\t0\t-1\t.\t0\t-1");
        assert!(nocall.slots.is_empty());

        // Ambiguity set splits on comma.
        let amb = call("HLA-K\t1\tHLA-K*01:01,HLA-K*01:04\t29.5\t17\t.\t0\t-1");
        assert_eq!(amb.slots[0].alleles, vec!["HLA-K*01:01", "HLA-K*01:04"]);
    }

    #[test]
    fn only_first_allele_of_a_set_votes_and_topn_selects() {
        // Two samples both het A*01/A*02; a third, low-quality A*03 should lose.
        let s = |n: &str, l: &[&str]| (n.to_owned(), l.iter().map(|x| call(x)).collect::<Vec<_>>());
        let samples = vec![
            s("s1", &["HLA-A\t2\tHLA-A*01:01\t50\t40\tHLA-A*02:01\t48\t38"]),
            s("s2", &["HLA-A\t2\tHLA-A*01:01\t52\t42\tHLA-A*02:01\t45\t36"]),
            s("s3", &["HLA-A\t1\tHLA-A*03:01\t10\t5\t.\t0\t-1"]),
        ];
        let m = combine(&samples, &CombineConfig::default());
        // top-2 per gene => A*01:01 (82) and A*02:01 (74); A*03:01 (5) drops (< tq=30 and rank 3).
        assert_eq!(m.alleles, vec!["HLA-A*01:01", "HLA-A*02:01"]);
    }

    #[test]
    fn boundary_tie_keeps_first_voted_allele_matching_t1k_merge() {
        // When more than `n` alleles for a gene tie at the top-N boundary,
        // `t1k-merge.py` (a stable reverse sort over dict insertion order) keeps
        // the *first-voted* allele, not the alphabetically-first one. Here
        // A*02:01 (q50) wins outright; A*99:01 and A*01:01 both tally 40, and
        // A*99:01 is voted first (s1, before s2's A*01:01), so top-2 = {A*02:01,
        // A*99:01} and s2's A*01:01 falls out as an inconsistency. A name-based
        // tie-break would wrongly pick A*01:01 instead. Byte-for-byte identical
        // to `python3 t1k-merge.py -l list -n 2` on this input.
        let s = |n: &str, l: &str| (n.to_owned(), vec![call(l)]);
        let samples = vec![
            s("s1", "HLA-A\t2\tHLA-A*99:01\t10\t40\tHLA-A*02:01\t10\t50"),
            s("s2", "HLA-A\t1\tHLA-A*01:01\t10\t40\t.\t0\t-1"),
        ];
        let m = combine(&samples, &CombineConfig::default());
        // Columns are the sorted representative union: A*02:01 then A*99:01.
        assert_eq!(m.alleles, vec!["HLA-A*02:01", "HLA-A*99:01"]);
        // s2's A*01:01 lost the tie -> recorded as an inconsistency, not a column.
        assert_eq!(m.rows[1].inconsistencies, vec!["HLA-A*01:01_10_40"]);
    }

    #[test]
    fn total_quality_floor_excludes_low_vote_alleles() {
        let s = |n: &str, l: &str| (n.to_owned(), vec![call(l)]);
        // Single sample: A*02:01 quality 20 < tq(30) => excluded entirely.
        let samples = vec![s("s1", "HLA-A\t2\tHLA-A*01:01\t50\t40\tHLA-A*02:01\t48\t20")];
        let m = combine(&samples, &CombineConfig::default());
        assert_eq!(m.alleles, vec!["HLA-A*01:01"]);
        // s1's A*02:01 slot (qual 20 > min_quality 0) is not representative => inconsistency.
        assert_eq!(m.rows[0].inconsistencies, vec!["HLA-A*02:01_48_20"]);
        // A*01:01 count is its abundance; matrix cell is the float.
        assert_eq!(m.rows[0].counts, vec![Some(50.0)]);
    }

    #[test]
    fn min_quality_gates_voting_and_fill_entirely() {
        // A slot at quality <= min_quality neither votes nor appears as an
        // inconsistency (T1K's `quality > q` gate applies to both passes). This
        // is the path where t1k-merge.py crashes on a nonzero -q; unum treats it
        // as numeric.
        let s = |n: &str, l: &str| (n.to_owned(), vec![call(l)]);
        let cfg = CombineConfig { min_quality: 25.0, ..CombineConfig::default() };
        // A*01:01 qual 40 (passes), A*02:01 qual 20 (<= 25, gated out).
        let samples = vec![s("s1", "HLA-A\t2\tHLA-A*01:01\t50\t40\tHLA-A*02:01\t48\t20")];
        let m = combine(&samples, &cfg);
        assert_eq!(m.alleles, vec!["HLA-A*01:01"]);
        // The gated-out A*02:01 slot produces NO inconsistency record.
        assert!(m.rows[0].inconsistencies.is_empty());
        assert_eq!(m.rows[0].counts, vec![Some(50.0)]);
    }

    #[test]
    fn untouched_cell_is_integer_zero_touched_whole_number_keeps_trailing_zero() {
        // Faithful to `t1k-merge.py`, which initializes each cell to the integer
        // `0` and does `+= float(...)` on a hit: an untouched cell prints as `0`,
        // a touched one via Python `str(float)`. So a touched whole-number
        // abundance keeps its trailing `.0` (`str(50.0) == "50.0"`), while an
        // untouched cell stays `0`.
        let s = |n: &str, l: &str| (n.to_owned(), vec![call(l)]);
        let samples = vec![
            s("s1", "HLA-A\t2\tHLA-A*01:01\t50\t40\tHLA-A*02:01\t48\t38"),
            s("s2", "HLA-A\t1\tHLA-A*01:01\t99\t60\t.\t0\t-1"), // no A*02:01
        ];
        let m = combine(&samples, &CombineConfig::default());
        let text = format_matrix(&m);
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "sample\tHLA-A*01:01\tHLA-A*02:01\tinconsistency");
        // Both of s1's cells are touched whole numbers -> "50.0" and "48.0".
        assert_eq!(lines[1], "s1\t50.0\t48.0\t");
        // s2's A*01:01 is a touched whole number -> "99.0"; A*02:01 is untouched -> "0".
        assert_eq!(lines[2], "s2\t99.0\t0\t");
    }

    #[test]
    fn touched_fractional_abundance_matches_python_str_float() {
        // Non-whole abundances already agree between Rust `Display` and Python
        // `str(float)` (both use a shortest-round-trip representation). A summed
        // touched cell (61.5 + 57.4) must render identically to Python.
        let s = |n: &str, l: &str| (n.to_owned(), vec![call(l)]);
        let samples = vec![
            s("s1", "HLA-A\t2\tHLA-A*01:01\t61.5\t40\t.\t0\t-1"),
            s("s2", "HLA-A\t1\tHLA-A*01:01\t57.4\t60\t.\t0\t-1"),
        ];
        // Single-column matrix; both samples hit HLA-A*01:01.
        let m = combine(&samples, &CombineConfig::default());
        let text = format_matrix(&m);
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "sample\tHLA-A*01:01\tinconsistency");
        assert_eq!(lines[1], "s1\t61.5\t");
        assert_eq!(lines[2], "s2\t57.4\t");
    }
}
