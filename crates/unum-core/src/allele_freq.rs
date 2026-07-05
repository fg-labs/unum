//! AFND-style population allele-frequency table for the opt-in Hardy-Weinberg
//! selection prior (#29).
//!
//! One responsibility: map an allele *series* name (e.g. `A*01:01`) to a
//! Jeffreys-smoothed, per-locus-normalized population frequency. Pure data; no
//! dependency on the genotyper.
//!
//! # Model
//!
//! The table is parsed from a tab-separated file with columns
//! `allele<TAB>frequency<TAB>count`, one row per allele at some field depth
//! (an optional one-line header is skipped). **`count` is the source of
//! truth** -- integer observation counts, robust to per-source normalization;
//! the `frequency` column is parsed-position-wise but its value is ignored.
//!
//! For allele `a` in locus `L` the frequency is Jeffreys-smoothed and
//! per-locus-normalized:
//!
//! ```text
//! f(a) = (count(a) + α) / (Σ_{b ∈ L} count(b) + α · |L|)
//! ```
//!
//! with α = [`AlleleFreqTable::ALPHA`] = 0.5 (Jeffreys) and `|L|` the number of
//! distinct observed alleles in locus `L`. An allele that is *absent* from the
//! DB but whose locus *is* present gets the strictly-positive floor
//! `f = α / (Σ + α · |L|)` (rare-allele safeguard, never 0). A locus entirely
//! absent from the DB yields `None` -- the caller treats that gene's prior as
//! inactive.
//!
//! # Resolution-aware matching
//!
//! A query series is matched to DB rows at field granularity: a row matches a
//! query when one is a field-prefix of the other (compared at the shorter
//! depth). So a coarse DB row `A*01:01` matches a finer query `A*01:01:01:01`
//! (pushdown), and a coarse query `A*01:01` sums all finer DB rows that share
//! its `A*01:01` prefix (rollup). All matching rows' counts are summed.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::Context;

/// A parsed AFND-style allele-frequency table, retaining per-allele observation
/// counts grouped by locus for Jeffreys-smoothed, per-locus normalization.
///
/// Construct with [`AlleleFreqTable::from_tsv`]; query with
/// [`AlleleFreqTable::frequency`].
#[derive(Debug, Clone, Default)]
pub struct AlleleFreqTable {
    /// Per-locus data: `locus -> LocusFreq`. The locus key is the text before
    /// `*` in each allele name.
    loci: HashMap<String, LocusFreq>,
}

/// Per-locus observation counts. Each entry is one DB row's parsed field list
/// plus its count; the locus total and distinct-allele count are precomputed
/// for the normalization denominator.
#[derive(Debug, Clone, Default)]
struct LocusFreq {
    /// One `(fields, count)` per DB row within this locus. `fields` is the
    /// allele's series split on `:` (the leading gene/first-field element
    /// keeps the `LOCUS*NN` form; subsequent elements are the numeric fields).
    rows: Vec<(Vec<String>, f64)>,
    /// Σ of all `count`s in this locus (the `Σ count` term).
    total_count: f64,
    /// `|L|` -- the number of distinct observed alleles (rows) in this locus.
    distinct_alleles: usize,
}

impl AlleleFreqTable {
    /// The Jeffreys pseudocount used for per-locus smoothing.
    pub const ALPHA: f64 = 0.5;

    /// Parses an AFND-style TSV: `allele<TAB>frequency<TAB>count`, one row per
    /// allele. An optional one-line header (a first row whose third column is
    /// not a parseable count) is skipped. `count` is authoritative; the
    /// `frequency` column is ignored. Blank lines are skipped.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened/read, or if a non-header
    /// data row is malformed (fewer than three tab-separated fields, or a
    /// third field that does not parse as a non-negative number).
    pub fn from_tsv<P: AsRef<Path>>(path: P) -> anyhow::Result<Self> {
        let path = path.as_ref();
        let file = std::fs::File::open(path)
            .with_context(|| format!("opening allele-frequency TSV {}", path.display()))?;
        let reader = std::io::BufReader::new(file);

        let mut loci: HashMap<String, LocusFreq> = HashMap::new();
        let mut seen_data_row = false;

        for (line_no, line) in reader.lines().enumerate() {
            let line =
                line.with_context(|| format!("reading {} line {}", path.display(), line_no + 1))?;
            let line = line.trim_end_matches(['\r', '\n']);
            if line.is_empty() {
                continue;
            }

            let cols: Vec<&str> = line.split('\t').collect();
            let allele = cols.first().copied().unwrap_or("").trim();
            let count_str = cols.get(2).copied().unwrap_or("").trim();

            // A first-row header is skipped: if we have not yet consumed a data
            // row and the third column does not parse as a count, treat the row
            // as a header. Any later unparseable row is a hard error.
            let parsed_count = count_str.parse::<f64>();
            if !seen_data_row && (cols.len() < 3 || parsed_count.is_err()) {
                continue;
            }

            anyhow::ensure!(
                cols.len() >= 3,
                "malformed allele-frequency row at {} line {}: expected at least 3 \
                 tab-separated columns (allele, frequency, count), got {}",
                path.display(),
                line_no + 1,
                cols.len()
            );
            let count = parsed_count.with_context(|| {
                format!(
                    "parsing count column at {} line {}: {count_str:?}",
                    path.display(),
                    line_no + 1
                )
            })?;
            anyhow::ensure!(
                count.is_finite() && count >= 0.0,
                "invalid count at {} line {}: {count} (must be finite and non-negative)",
                path.display(),
                line_no + 1
            );

            seen_data_row = true;

            let (locus, fields) = parse_series(allele);
            let entry = loci.entry(locus).or_default();
            entry.rows.push((fields, count));
            entry.total_count += count;
            entry.distinct_alleles += 1;
        }

        Ok(Self { loci })
    }

    /// Jeffreys-smoothed, per-locus-normalized frequency for `series` (e.g.
    /// `"A*01:01"`), matched at field granularity (rollup/pushdown/sum, see the
    /// module docs).
    ///
    /// Returns `None` iff the locus is absent from the table -- the caller
    /// treats that as an inactive prior for the gene. When the locus is present
    /// but the specific series is not observed, returns the strictly-positive
    /// Jeffreys floor.
    #[must_use]
    pub fn frequency(&self, series: &str) -> Option<f64> {
        let (locus, query_fields) = parse_series(series);
        let locus_freq = self.loci.get(&locus)?;

        let alpha = Self::ALPHA;
        #[allow(clippy::cast_precision_loss)]
        let denom = locus_freq.total_count + alpha * (locus_freq.distinct_alleles as f64);
        // `distinct_alleles >= 1` here (the locus exists only if a row created
        // it), so `denom >= alpha > 0`.

        let mut matched_count = 0.0;
        for (row_fields, count) in &locus_freq.rows {
            if fields_match(&query_fields, row_fields) {
                matched_count += count;
            }
        }

        Some((matched_count + alpha) / denom)
    }
}

/// Splits an allele series into `(locus, fields)`. The locus is the text before
/// the first `*` (or the whole string if there is no `*`). `fields` is the full
/// series split on `:` -- element 0 keeps the `LOCUS*NN` head, and each
/// subsequent element is one numeric field. Used for locus grouping and for
/// resolution-aware field matching.
fn parse_series(series: &str) -> (String, Vec<String>) {
    let trimmed = series.trim();
    let locus = match trimmed.find('*') {
        Some(idx) => trimmed[..idx].to_string(),
        None => trimmed.to_string(),
    };
    let fields = trimmed.split(':').map(str::to_string).collect();
    (locus, fields)
}

/// Returns `true` when two field lists match at field granularity: one is a
/// prefix of the other, compared at the shorter list's depth. This makes a
/// coarse DB row match a finer query (pushdown) and a coarse query sum all
/// finer DB rows that share its prefix (rollup).
fn fields_match(a: &[String], b: &[String]) -> bool {
    let depth = a.len().min(b.len());
    a[..depth] == b[..depth]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Writes `contents` to a uniquely-named temp file and returns its path.
    /// Uses the process id + a monotonic counter to avoid collisions without a
    /// temp-file dependency.
    fn write_temp_tsv(contents: &str) -> std::path::PathBuf {
        use std::sync::atomic::{AtomicU64, Ordering};
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, Ordering::Relaxed);
        let mut path = std::env::temp_dir();
        path.push(format!("unum_allele_freq_test_{}_{n}.tsv", std::process::id()));
        let mut file = std::fs::File::create(&path).expect("create temp tsv");
        file.write_all(contents.as_bytes()).expect("write temp tsv");
        path
    }

    /// (a) Per-locus normalization: the observed alleles' frequencies sum to
    /// ~1 when the Jeffreys floor mass on unobserved alleles is added back.
    /// Equivalently, Σ over observed alleles of `(count+α)/(Σ+α|L|)` equals
    /// `(Σ + α|L|)/(Σ + α|L|) = 1`.
    #[test]
    fn allele_freq_per_locus_normalization_sums_to_one() {
        let tsv = "allele\tfrequency\tcount\n\
                   A*01:01\t0.5\t30\n\
                   A*02:01\t0.3\t20\n\
                   A*03:01\t0.2\t10\n";
        let path = write_temp_tsv(tsv);
        let table = AlleleFreqTable::from_tsv(&path).expect("parse");

        let f1 = table.frequency("A*01:01").expect("locus present");
        let f2 = table.frequency("A*02:01").expect("locus present");
        let f3 = table.frequency("A*03:01").expect("locus present");
        // Σ = 60, |L| = 3, denom = 60 + 0.5*3 = 61.5.
        let sum = f1 + f2 + f3;
        // Observed-allele mass = (30+20+10 + 3*0.5)/61.5 = 61.5/61.5 = 1.0.
        assert!((sum - 1.0).abs() < 1e-12, "sum={sum}");

        std::fs::remove_file(&path).ok();
    }

    /// (b) Jeffreys: a single-observed-allele locus gives
    /// `f = (n+0.5)/(n+0.5·|L|)` with `|L| = 1`.
    #[test]
    fn allele_freq_jeffreys_single_observed_allele() {
        let tsv = "B*07:02\t0.9\t42\n";
        let path = write_temp_tsv(tsv);
        let table = AlleleFreqTable::from_tsv(&path).expect("parse");

        let f = table.frequency("B*07:02").expect("locus present");
        // n = 42, |L| = 1: f = (42+0.5)/(42+0.5) = 1.0.
        assert!((f - 1.0).abs() < 1e-12, "f={f}");

        std::fs::remove_file(&path).ok();
    }

    /// (c) An allele absent from a present locus gets the positive Jeffreys
    /// floor `0.5/(Σ + 0.5·|L|)`, never 0.
    #[test]
    fn allele_freq_absent_allele_gets_positive_floor() {
        let tsv = "A*01:01\t0.6\t30\n\
                   A*02:01\t0.4\t20\n";
        let path = write_temp_tsv(tsv);
        let table = AlleleFreqTable::from_tsv(&path).expect("parse");

        let f = table.frequency("A*99:99").expect("locus present");
        // Σ = 50, |L| = 2, denom = 50 + 1.0 = 51; floor = 0.5/51.
        let expected = 0.5 / 51.0;
        assert!(f > 0.0, "floor must be strictly positive, f={f}");
        assert!((f - expected).abs() < 1e-12, "f={f} expected={expected}");

        std::fs::remove_file(&path).ok();
    }

    /// (d) A locus entirely absent from the table yields `None`.
    #[test]
    fn allele_freq_absent_locus_is_none() {
        let tsv = "A*01:01\t0.6\t30\n";
        let path = write_temp_tsv(tsv);
        let table = AlleleFreqTable::from_tsv(&path).expect("parse");

        assert!(table.frequency("DRB1*03:01").is_none());

        std::fs::remove_file(&path).ok();
    }

    /// (e) Resolution-aware: a 4-field query matches a 2-field table row
    /// (pushdown), and a 2-field query sums finer table rows (rollup).
    #[test]
    fn allele_freq_resolution_aware_pushdown_and_rollup() {
        // Pushdown: the table has only a 2-field row; a 4-field query resolves
        // to it (its prefix equals the row).
        let tsv_pushdown = "A*01:01\t0.5\t40\n\
                            A*02:01\t0.5\t10\n";
        let path = write_temp_tsv(tsv_pushdown);
        let table = AlleleFreqTable::from_tsv(&path).expect("parse");
        let f_fine = table.frequency("A*01:01:01:01").expect("locus present");
        let f_coarse = table.frequency("A*01:01").expect("locus present");
        assert!(
            (f_fine - f_coarse).abs() < 1e-12,
            "4-field query should match the 2-field row: {f_fine} vs {f_coarse}"
        );
        std::fs::remove_file(&path).ok();

        // Rollup: the table has two finer rows under A*01:01; a 2-field query
        // sums both counts.
        let tsv_rollup = "A*01:01:01\t0.3\t25\n\
                          A*01:01:02\t0.2\t15\n\
                          A*02:01:01\t0.5\t10\n";
        let path = write_temp_tsv(tsv_rollup);
        let table = AlleleFreqTable::from_tsv(&path).expect("parse");
        let f_rollup = table.frequency("A*01:01").expect("locus present");
        // matched = 25 + 15 = 40; Σ = 50, |L| = 3, denom = 50 + 1.5 = 51.5.
        let expected = (40.0 + 0.5) / 51.5;
        assert!((f_rollup - expected).abs() < 1e-12, "rollup f={f_rollup} expected={expected}");
        std::fs::remove_file(&path).ok();
    }

    /// (f) The `frequency` column is ignored: two rows with identical counts
    /// but different freq columns produce identical frequencies.
    #[test]
    fn allele_freq_frequency_column_ignored() {
        let tsv_a = "A*01:01\t0.9\t30\n\
                     A*02:01\t0.1\t30\n";
        let tsv_b = "A*01:01\t0.1\t30\n\
                     A*02:01\t0.9\t30\n";
        let path_a = write_temp_tsv(tsv_a);
        let path_b = write_temp_tsv(tsv_b);
        let table_a = AlleleFreqTable::from_tsv(&path_a).expect("parse a");
        let table_b = AlleleFreqTable::from_tsv(&path_b).expect("parse b");

        let fa = table_a.frequency("A*01:01").expect("present");
        let fb = table_b.frequency("A*01:01").expect("present");
        assert!((fa - fb).abs() < 1e-12, "frequency column must be ignored: {fa} vs {fb}");

        std::fs::remove_file(&path_a).ok();
        std::fs::remove_file(&path_b).ok();
    }
}
