//! Infer per-allele copy number from a genotype result -- the Rust port of
//! T1K's `t1k-copynumber.py` (the `unum copy-number` subcommand). Copy-number
//! inference is meaningful for the copy-number-variable KIR locus (genes present
//! in 0..N copies); for the fixed-two-copy HLA loci it is degenerate. The
//! algorithm itself is locus-agnostic -- it operates purely on the abundances in
//! a `{prefix}_genotype.tsv`.
//!
//! # Algorithm (faithful to `t1k-copynumber.py`)
//!
//! 1. Read the genotype file; per gene keep the alleles whose quality is
//!    `> min_quality`, recording each allele's abundance.
//! 2. Estimate the **one-copy** abundance distribution `N(mean, var)` from a
//!    `sqrt`-transformed sample: for `--nomissing` genes (assumed present on
//!    every chromosome) use every allele's `sqrt(abund)` (a homozygous such gene
//!    contributes `sqrt(abund)/2`, since its single allele spans two copies);
//!    then extend with the lower `[start:end)` quantile band of the sorted
//!    `sqrt(abund)` of all *heterozygous, non-nomissing* alleles. `var` is scaled
//!    by `--adjust-var`.
//! 3. For each allele, pick the copy number `c in 1..=8` maximizing the
//!    log-likelihood of `sqrt(abund)` under `N(mean*c, var*c)`; report `c` and a
//!    confidence `ratio` (log-likelihood gap to the runner-up copy number).
//!
//! The quantile band indices are computed from the *total* passing-allele count
//! (minus nomissing), exactly as the Python does, then applied to the
//! heterozygous list (clamped) -- preserved verbatim for byte-identity.

use anyhow::{Result, bail};
use std::collections::{HashMap, HashSet};

/// One allele's abundance and quality parsed from a genotype line.
#[derive(Debug, Clone, PartialEq)]
pub struct AlleleAbundance {
    /// Allele name, verbatim (the whole genotype-column field, comma-set included).
    pub allele: String,
    /// Abundance (genotype column `k+1`).
    pub abundance: f64,
    /// Genotype quality (column `k+2`); alleles with `quality <= min_quality` are dropped.
    pub quality: i64,
}

/// One gene's parsed alleles, in file order (before quality filtering).
#[derive(Debug, Clone, PartialEq)]
pub struct GeneAlleles {
    /// Gene name (column 1).
    pub gene: String,
    /// Alleles read from the line's active slots (gated by the genotype's own count).
    pub alleles: Vec<AlleleAbundance>,
}

/// Configuration mirroring `t1k-copynumber.py`'s flags.
#[derive(Debug, Clone)]
pub struct CopyNumberConfig {
    /// `--nomissing`: genes assumed present on every chromosome (used to anchor
    /// the one-copy parameters); order preserved.
    pub nomissing_genes: Vec<String>,
    /// `--upper-quantile` (default 0.3).
    pub upper_quantile: f64,
    /// `--lower-quantile` (default 0.0).
    pub lower_quantile: f64,
    /// `--adjust-var` (default 1.0).
    pub adjust_var: f64,
    /// `-q`: keep alleles with `quality > min_quality` (default 0).
    pub min_quality: i64,
}

impl Default for CopyNumberConfig {
    fn default() -> Self {
        Self {
            nomissing_genes: Vec::new(),
            upper_quantile: 0.3,
            lower_quantile: 0.0,
            adjust_var: 1.0,
            min_quality: 0,
        }
    }
}

impl CopyNumberConfig {
    /// Validate the numeric parameters at the CLI boundary, before any inference.
    ///
    /// Requires finite quantiles in `[0, 1]` with `lower_quantile <= upper_quantile`,
    /// and a finite `adjust_var > 0`. This is a deliberate divergence from
    /// `t1k-copynumber.py`, which does not validate these flags (see
    /// `docs/DIVERGENCES.md`): out-of-range values only ever select no anchors or
    /// drive `infer_copy_numbers` into its non-positive-variance error, so an
    /// up-front, actionable argument error is strictly friendlier.
    ///
    /// # Errors
    ///
    /// Returns an error naming the first out-of-range or non-finite parameter.
    pub fn validate(&self) -> Result<()> {
        if !(self.upper_quantile.is_finite() && (0.0..=1.0).contains(&self.upper_quantile)) {
            bail!(
                "copy-number: --upper-quantile must be a finite value in [0, 1] (got {})",
                self.upper_quantile
            );
        }
        if !(self.lower_quantile.is_finite() && (0.0..=1.0).contains(&self.lower_quantile)) {
            bail!(
                "copy-number: --lower-quantile must be a finite value in [0, 1] (got {})",
                self.lower_quantile
            );
        }
        if self.lower_quantile > self.upper_quantile {
            bail!(
                "copy-number: --lower-quantile ({}) must not exceed --upper-quantile ({})",
                self.lower_quantile,
                self.upper_quantile
            );
        }
        if !(self.adjust_var.is_finite() && self.adjust_var > 0.0) {
            bail!(
                "copy-number: --adjust-var must be a finite value greater than 0 (got {})",
                self.adjust_var
            );
        }
        Ok(())
    }
}

/// One output row: a gene, its passing-allele count, and up to two
/// `(allele, copy, ratio)` slots (missing slots render as `. -1 0`).
#[derive(Debug, Clone, PartialEq)]
pub struct CopyNumberRow {
    /// Gene name.
    pub gene: String,
    /// Number of alleles that passed the quality filter.
    pub allele_count: usize,
    /// Per passing allele (in file order): `(allele, copy_number, ratio)`.
    pub calls: Vec<(String, i64, f64)>,
}

fn abund_transform(x: f64) -> f64 {
    x.sqrt()
}

/// Log of the normal likelihood without the constant factor:
/// `-0.5 * ((x - mu)/sigma)^2 - ln(sigma)`, with `sigma = sqrt(var)`.
fn log_normal_likelihood(x: f64, mu: f64, var: f64) -> f64 {
    let sigma = var.sqrt();
    -0.5 * ((x - mu) / sigma).powi(2) - sigma.ln()
}

/// Parse one genotype line into a [`GeneAlleles`] (all active slots, unfiltered),
/// or `None` if the line lacks the gene + count columns. Fields are split on
/// whitespace (matching the Python `str.split()`); the allele field is kept
/// verbatim (comma-sets are not split).
#[must_use]
pub fn parse_genotype_line(line: &str) -> Option<GeneAlleles> {
    let cols: Vec<&str> = line.split_whitespace().collect();
    if cols.len() < 2 {
        return None;
    }
    let gene_copy: usize = cols[1].parse().ok()?;
    let mut alleles = Vec::new();
    for i in 0..gene_copy {
        let k = if i == 0 { 2 } else { 5 };
        if cols.len() <= k + 2 {
            break;
        }
        let Ok(abundance) = cols[k + 1].parse::<f64>() else { continue };
        let Ok(quality) = cols[k + 2].parse::<i64>() else { continue };
        alleles.push(AlleleAbundance { allele: cols[k].to_owned(), abundance, quality });
    }
    Some(GeneAlleles { gene: cols[0].to_owned(), alleles })
}

/// Infer copy numbers for every gene in `genes` (in the given order) under
/// `cfg`. Pure: no I/O. Returns one [`CopyNumberRow`] per gene.
///
/// # Errors
///
/// Returns an error when alleles need copy calls but the anchor sample cannot
/// yield a positive-variance one-copy distribution (empty, single, or all-equal
/// anchors) -- the point at which `t1k-copynumber.py` raises `ZeroDivisionError`
/// and produces no output; unum fails cleanly rather than emitting `NaN` copy
/// calls.
// Casts here are intentional and bounded: they mirror the Python's `int()`
// truncation of quantile indices and its float arithmetic. Allele counts are
// tiny and copy numbers are 1..=8, so no meaningful precision/sign loss occurs.
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub fn infer_copy_numbers(
    genes: &[GeneAlleles],
    cfg: &CopyNumberConfig,
) -> Result<Vec<CopyNumberRow>> {
    // Quality filter -> per-gene passing alleles (file order) + abundance lookup.
    let mut passing: Vec<(&str, Vec<&AlleleAbundance>)> = Vec::with_capacity(genes.len());
    let mut abund_of: HashMap<&str, f64> = HashMap::new();
    let mut total_passing = 0usize;
    for g in genes {
        let kept: Vec<&AlleleAbundance> =
            g.alleles.iter().filter(|a| a.quality > cfg.min_quality).collect();
        for a in &kept {
            abund_of.insert(a.allele.as_str(), a.abundance);
        }
        total_passing += kept.len();
        passing.push((g.gene.as_str(), kept));
    }
    let nomissing: HashSet<&str> = cfg.nomissing_genes.iter().map(String::as_str).collect();
    let gene_alleles: HashMap<&str, &Vec<&AlleleAbundance>> =
        passing.iter().map(|(gname, alls)| (*gname, alls)).collect();

    // --- estimate one-copy N(mean, var) ---
    let mut abundances: Vec<f64> = Vec::new();
    let mut used_alleles = 0usize;
    if !nomissing.is_empty() {
        for g in &cfg.nomissing_genes {
            let Some(alls) = gene_alleles.get(g.as_str()) else { continue };
            if alls.len() > 1 {
                abundances.extend(alls.iter().map(|a| abund_transform(a.abundance)));
            } else if alls.len() == 1 {
                abundances.push(abund_transform(alls[0].abundance) / 2.0);
            }
            used_alleles += alls.len();
        }
    }

    // Slice indices are computed from the TOTAL passing count (minus nomissing),
    // exactly as the Python -- then applied to the heterozygous list below.
    let span = total_passing.saturating_sub(used_alleles) as f64;
    let start = (span * cfg.lower_quantile).trunc() as usize;
    let end = (span * cfg.upper_quantile).trunc() as usize;

    let mut seen: HashSet<&str> = HashSet::new();
    let mut heter: Vec<f64> = Vec::new();
    for (gname, alls) in &passing {
        if nomissing.contains(gname) || alls.len() <= 1 {
            continue;
        }
        for a in alls {
            if seen.insert(a.allele.as_str()) {
                heter.push(abund_transform(a.abundance));
            }
        }
    }
    heter.sort_by(f64::total_cmp);
    let (s, e) = (start.min(heter.len()), end.min(heter.len()));
    if s < e {
        abundances.extend_from_slice(&heter[s..e]);
    }

    let n = abundances.len() as f64;
    let mean = abundances.iter().sum::<f64>() / n;
    let var = (abundances.iter().map(|a| a * a).sum::<f64>() / n - mean * mean) * cfg.adjust_var;

    // The one-copy `N(mean, var)` must be usable before any copy call: a
    // non-positive or non-finite `var` gives `sigma = sqrt(var) <= 0`, so
    // `log_normal_likelihood` divides by zero and returns `NaN`. This is exactly
    // where `t1k-copynumber.py` raises `ZeroDivisionError` and emits nothing --
    // an empty anchor sample (`n == 0` -> `var` is `NaN`), a single anchor, or
    // all-equal anchors (`var == 0`). Fail cleanly instead of emitting `NaN`
    // copy/ratio values. When no allele passed the quality filter
    // (`total_passing == 0`) there is nothing to call -- `copy_of` is never
    // invoked, so the params are unused and we fall through to sentinel rows.
    if total_passing > 0 && (var <= 0.0 || var.is_nan()) {
        bail!(
            "copy-number: cannot estimate a usable one-copy distribution from the anchor \
             alleles (variance is not positive); need at least two heterozygous \
             non-nomissing alleles with differing abundance, or present --nomissing gene(s)"
        );
    }

    // --- per-allele copy number (argmax over 1..=8) + confidence ratio ---
    let copy_of = |abund: f64| -> (i64, f64) {
        let x = abund_transform(abund);
        let mut lls: Vec<(i64, f64)> = (1..=8)
            .map(|c| (c, log_normal_likelihood(x, mean * c as f64, var * c as f64)))
            .collect();
        lls.sort_by(|a, b| b.1.total_cmp(&a.1));
        (lls[0].0, lls[0].1 - lls[1].1)
    };

    // --- assemble output rows (one per gene, in order) ---
    Ok(passing
        .iter()
        .map(|(gname, alls)| {
            let calls = alls
                .iter()
                .map(|a| {
                    let (copy, ratio) = copy_of(abund_of[a.allele.as_str()]);
                    (a.allele.clone(), copy, ratio)
                })
                .collect();
            CopyNumberRow { gene: (*gname).to_owned(), allele_count: alls.len(), calls }
        })
        .collect())
}

/// Render rows as the tab-separated text `t1k-copynumber.py` prints: per gene,
/// `gene  count  [allele copy ratio]x2`, missing slots as `. -1 0`, ratio to 2dp.
#[must_use]
pub fn format_rows(rows: &[CopyNumberRow]) -> String {
    let mut out = String::new();
    for row in rows {
        out.push_str(&format!("{}\t{}", row.gene, row.allele_count));
        for i in 0..2 {
            if let Some((allele, copy, ratio)) = row.calls.get(i) {
                out.push_str(&format!("\t{allele}\t{copy}\t{ratio:.2}"));
            } else {
                out.push_str("\t.\t-1\t0");
            }
        }
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genes(lines: &[&str]) -> Vec<GeneAlleles> {
        lines.iter().filter_map(|l| parse_genotype_line(l)).collect()
    }

    #[test]
    fn parses_active_slots_and_gates_are_slot_counted() {
        let g =
            parse_genotype_line("KIR2DL1\t2\tKIR2DL1*035\t42.9661\t15\tKIR2DL1*003\t42.0886\t14")
                .unwrap();
        assert_eq!(g.gene, "KIR2DL1");
        assert_eq!(g.alleles.len(), 2);
        assert_eq!(g.alleles[0].allele, "KIR2DL1*035");
        assert!((g.alleles[0].abundance - 42.9661).abs() < 1e-9);
        assert_eq!(g.alleles[1].quality, 14);

        // No-call: count 0 -> no alleles, but the gene is still present.
        let z = parse_genotype_line("KIR2DL4\t0\t.\t0\t-1\t.\t0\t-1").unwrap();
        assert_eq!(z.gene, "KIR2DL4");
        assert!(z.alleles.is_empty());
    }

    #[test]
    fn validate_accepts_defaults_and_valid_ranges() {
        CopyNumberConfig::default().validate().expect("defaults are valid");
        let cfg = CopyNumberConfig {
            lower_quantile: 0.1,
            upper_quantile: 0.1, // lower == upper is allowed
            adjust_var: 2.5,
            ..Default::default()
        };
        cfg.validate().expect("valid ranges accepted");
    }

    #[test]
    fn validate_rejects_out_of_range_and_non_finite_params() {
        // Each case perturbs exactly one field out of range; all must be rejected.
        let bad = [
            CopyNumberConfig { upper_quantile: 1.5, ..Default::default() },
            CopyNumberConfig { upper_quantile: -0.1, ..Default::default() },
            CopyNumberConfig { lower_quantile: 1.5, ..Default::default() },
            CopyNumberConfig { lower_quantile: -0.1, ..Default::default() },
            CopyNumberConfig { upper_quantile: f64::NAN, ..Default::default() },
            CopyNumberConfig { lower_quantile: f64::INFINITY, ..Default::default() },
            // lower > upper
            CopyNumberConfig { lower_quantile: 0.6, upper_quantile: 0.3, ..Default::default() },
            // adjust_var must be finite and strictly positive
            CopyNumberConfig { adjust_var: 0.0, ..Default::default() },
            CopyNumberConfig { adjust_var: -1.0, ..Default::default() },
            CopyNumberConfig { adjust_var: f64::NAN, ..Default::default() },
            CopyNumberConfig { adjust_var: f64::INFINITY, ..Default::default() },
        ];
        for cfg in bad {
            assert!(cfg.validate().is_err(), "expected rejection for {cfg:?}");
        }
    }

    #[test]
    fn homozygous_double_abundance_infers_two_copies() {
        // A het anchor near abundance ~42 sets one-copy mean; a lone allele at
        // ~2x that abundance should be called copy 2 (the KIR2DL3 pattern).
        let gs = genes(&[
            "GA\t2\tGA*01\t42.0\t20\tGA*02\t42.0\t20",
            "GB\t2\tGB*01\t40.0\t20\tGB*02\t44.0\t20",
            "GC\t1\tGC*01\t84.0\t38\t.\t0\t-1",
        ]);
        // Use a wide quantile band so the het alleles anchor one-copy params.
        let cfg =
            CopyNumberConfig { lower_quantile: 0.0, upper_quantile: 1.0, ..Default::default() };
        let rows = infer_copy_numbers(&gs, &cfg).expect("infer");
        let gc = rows.iter().find(|r| r.gene == "GC").unwrap();
        assert_eq!(gc.calls[0].1, 2, "lone allele at 2x one-copy abundance -> copy 2");
        let ga = rows.iter().find(|r| r.gene == "GA").unwrap();
        assert_eq!(ga.calls[0].1, 1, "het alleles at one-copy abundance -> copy 1");
    }

    #[test]
    fn quality_filter_drops_alleles_but_keeps_gene_row() {
        // G loses its low-quality second allele (count 2 -> 1 passing); a
        // separate het gene H anchors the one-copy params so G*01 still gets a
        // valid (finite) copy call and G's row is retained.
        let gs = genes(&[
            "G\t2\tG*01\t40.0\t20\tG*02\t5.0\t0", // G*02 qual 0 <= min(0) -> dropped
            "H\t2\tH*01\t40.0\t20\tH*02\t44.0\t20", // het anchor, both pass
        ]);
        let cfg =
            CopyNumberConfig { lower_quantile: 0.0, upper_quantile: 1.0, ..Default::default() };
        let rows = infer_copy_numbers(&gs, &cfg).expect("infer");
        let g = rows.iter().find(|r| r.gene == "G").expect("G row");
        assert_eq!(g.allele_count, 1); // only G*01 survives the filter
        assert_eq!(g.calls.len(), 1);
        assert!(g.calls[0].2.is_finite(), "copy-number ratio must be finite, not NaN");
    }

    #[test]
    fn nonzero_min_quality_filters_numerically() {
        // The path where t1k-copynumber.py crashes (`int <= str` on a passed -q);
        // unum treats -q as numeric. Both alleles (qual 15, 14) fall below 20.
        let gs = genes(&["KIR2DL1\t2\tKIR2DL1*035\t42.9\t15\tKIR2DL1*003\t42.0\t14"]);
        let cfg = CopyNumberConfig { min_quality: 20, ..Default::default() };
        let rows = infer_copy_numbers(&gs, &cfg).expect("infer");
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].allele_count, 0);
        assert!(rows[0].calls.is_empty());
    }

    #[test]
    fn no_anchor_alleles_errors_instead_of_nan() {
        // No --nomissing and every passing gene homozygous -> empty anchor
        // sample (var is NaN), but there are still alleles that need a copy call.
        // t1k-copynumber.py divides by zero here (ZeroDivisionError, no output);
        // unum errors cleanly rather than emitting NaN copy/ratio values.
        let gs = genes(&["GC\t1\tGC*01\t84.0\t30\t.\t0\t-1", "GD\t1\tGD*01\t42.0\t30\t.\t0\t-1"]);
        let err = infer_copy_numbers(&gs, &CopyNumberConfig::default()).unwrap_err();
        assert!(err.to_string().contains("one-copy distribution"), "unexpected error: {err}");
    }

    #[test]
    fn single_anchor_zero_variance_errors_instead_of_nan() {
        // A single anchor allele (or all-equal anchors) gives var == 0, so
        // sigma == 0 and the log-likelihood is NaN -- t1k-copynumber.py likewise
        // divides by zero (at the likelihood, not the mean). Here the default
        // upper-quantile band collapses to one anchor. unum must error, not emit
        // NaN, matching the empty-sample case (a sibling of the same defect).
        let gs = genes(&[
            "GA\t2\tGA*01\t40.0\t20\tGA*02\t44.0\t20", // het -> 2 candidate anchors
            "GB\t1\tGB*01\t84.0\t30\t.\t0\t-1",        // homozygous, still needs a call
        ]);
        // total_passing == 3, so the band [0 : trunc(3 * 0.5)] = [0:1] keeps
        // exactly one anchor abundance -> var == 0.
        let cfg =
            CopyNumberConfig { lower_quantile: 0.0, upper_quantile: 0.5, ..Default::default() };
        let err = infer_copy_numbers(&gs, &cfg).unwrap_err();
        assert!(err.to_string().contains("one-copy distribution"), "unexpected error: {err}");
    }

    #[test]
    fn all_alleles_filtered_emits_sentinel_rows_without_error() {
        // When the quality filter removes *every* allele, there is nothing to
        // call: copy_of is never invoked, so the missing anchors are harmless.
        // The gene row is still emitted with sentinel slots (no error, no NaN).
        let gs = genes(&["KIR2DL1\t2\tKIR2DL1*035\t42.9\t15\tKIR2DL1*003\t42.0\t14"]);
        let cfg = CopyNumberConfig { min_quality: 20, ..Default::default() };
        let rows = infer_copy_numbers(&gs, &cfg).expect("infer");
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].allele_count, 0);
        assert!(rows[0].calls.is_empty());
    }

    #[test]
    fn format_pads_missing_slots_and_rounds_ratio() {
        let rows = vec![
            CopyNumberRow {
                gene: "G".into(),
                allele_count: 1,
                calls: vec![("G*01".into(), 2, 32.4236)],
            },
            CopyNumberRow { gene: "H".into(), allele_count: 0, calls: vec![] },
        ];
        let text = format_rows(&rows);
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "G\t1\tG*01\t2\t32.42\t.\t-1\t0");
        assert_eq!(lines[1], "H\t0\t.\t-1\t0\t.\t-1\t0");
    }
}
