# Design: population-frequency prior (#29) and identifiability reporting (#33)

**Status:** proposed — implementable. Injection points pinned by function+line against
`crates/unum-core/src/genotyper.rs` @ this worktree (5975 lines).
**Issues:** #29 (population allele-frequency prior on selection), #33 (per-allele identifiability, report-only).
**Author:** Nils Homer

## Goal

Two opt-in, default-off, well-cited capabilities on `unum genotype`, byte-identical to the
vendored C++ T1K oracle when off:

1. **#29** — a Hardy–Weinberg population allele-frequency prior on allele selection, applied as a
   *bounded, coverage-vanishing* additive term to the selector's own objective that breaks
   near-ties and tips weak evidence but **never overrides strong read evidence**.
2. **#33** — a per-called-allele identifiability measure, **report-only**, written to the existing
   opt-in `{prefix}_metrics.tsv` (`--emit-metrics`), never changing the call.

Non-goals: changing default (oracle-parity) behavior; bundling frequency DBs; population-conditioned
frequencies; LD/haplotype priors (all future work).

## Hard byte-identity invariant

When `--allele-freq` is absent, `select_alleles_for_genes` and every frozen output
(`_genotype.tsv`, `_allele.tsv`, aligned FASTAs) are byte-identical to the current oracle — verified
by `cargo ci-test-sys` (`crates/unum-sys/tests/diff_genotype_e2e.rs`), which must stay green.
**Additionally byte-identical when `--allele-freq` is present but the effective prior is constant
across the compared candidates** (empty table, or a locus/loci absent from the table): the prior
term must be *skipped entirely*, not added-as-zero — see §29 "empty-table guarantee". #33 writes only
to `_metrics.tsv` and never touches selection, so it is byte-identical to frozen outputs even when on.

## Prior art / why

- **Polysolver** (Shukla 2015, *Nat Biotechnol*): the canonical explicit-prior HLA caller — a
  Bayesian posterior = read likelihood × HWE population-frequency prior, MAP over allele *pairs*.
  Motivates applying the prior at the genotype (allele-pair) level, where HWE lives.
- **arcasHLA** (Orenbuch 2020, *Bioinformatics*): EM-over-equivalence-classes (T1K's shape) that
  weights resolution by AFND two-field frequencies with Dirichlet smoothing, explicitly to *"break
  ties between alleles indistinguishable given the reads"* — never to override evidence.
- **HWE genotype prior** (Wigginton 2005): P(hom A)=f_A², P(het A,B)=2 f_A f_B.
- **fghla** (this author's caller): portable template — log-additive HWE prior, per-locus-normalized
  AFND *counts* (count is source of truth, not the freq column), Jeffreys 0.5 pseudocount (rare
  alleles floored not zeroed), prior magnitude bounded far below per-read evidence.
- **Identifiability**: kallisto/RSEM/Salmon (Bray 2016, Li 2011, Patro 2017) — references sharing
  *all* equivalence classes are non-identifiable; measure and report at the identifiable resolution.
  Paunić 2012 — Shannon entropy over the compatible allele set. HLA\*LA (Dilthey 2019) — G-group
  reporting. Frequency data: CIWD 3.0 (Hurley 2020), AFND (Gonzalez-Galarza).

---

## #29 — population-frequency prior

### Input

New `--allele-freq <TSV>` on `GenotypeArgs` (`crates/unum/src/cli.rs`, alongside `emit_metrics` at
:222). Default `None` → feature off. Format matches fghla's AFND export: tab-separated
`allele<TAB>frequency<TAB>count`, one row per allele at some field depth; optional one-line header
skipped. **`count` is the source of truth** (integer observation counts, robust to per-source
normalization); the `frequency` column is parsed but ignored. No data is bundled; docs show how to
build the TSV from CIWD 3.0 / AFND.

### Frequency model

- **Per-locus normalization + Jeffreys.** For allele `a` in locus `L`:
  `f(a) = (count(a) + α) / (Σ_{b∈L} count(b) + α·|L|)`, α=0.5. An allele absent from the DB but in
  locus `L` gets `f = α / (Σ + α·|L|)` — small but strictly positive (rare-allele safeguard).
- **Resolution-aware matching.** Match a called allele to DB rows at its own field depth via the
  existing `Genotyper::parse_allele_name` (genotyper.rs:2169; core parser :210). The called allele's
  **major-allele series** is the natural key: `major_allele_idx_to_name` (:1972) already holds the
  series string (e.g. `A*01:01`) that the caller resolves to. Sum DB counts of all rows whose
  parsed series (at the caller's `fields_type`) equals the called series. Coarser DB rows roll up;
  finer DB rows sum. Unit-tested.
- **Locus absent from DB** → prior is a no-op for that gene (constant f across its candidates → the
  skip rule in "empty-table guarantee" fires, byte-identical).
- **Null/expression variants** (`N`/`L`/`S`/`C`/`A`/`Q` suffix): a fixed log-penalty biasing the
  selector away from calling non-expressed alleles. **DEFERRED to a follow-up PR** — the original
  `--allele-freq-null-penalty` flag shipped here as plumbing only (never read in any objective), and
  the concrete wiring (detection granularity, injection points, span cap) was never specified. It is
  designed in `docs/2026-07-08-null-penalty-wiring-design.md` (workspace-level) and will land on top
  of the merged frequency prior. The flag/field/setter are removed here until then.

### Where the prior enters — BOTH selection paths (the injection-coverage call)

**Decision: inject into BOTH the >2-allele-type pair search AND the ≤2-type greedy zygosity
decision.** Scoping to only the pair search would be *indefensible*: the ≤2-type greedy path
(`select_alleles_for_genes`, :3561) is exactly where the majority of genes with a clean 1- or
2-allele picture are decided, including every hom-vs-het call for those genes (the pair search
`continue`s at `allele_type_cnt <= 2`, :3832). A frequency prior that fires only on the messy >2-type
genes would be the opposite of principled. So the prior is factored into ONE shared helper and called
from both paths at the two argmax/tie sites below.

**Unifying principle (governs BOTH paths, subsumes the activeness gate and the no-override bound).**
*The prior changes a call — including the candidate SET, the objective, and the zygosity — only when
the candidate coverage margin is within the prior span `w·|Δ ln P_HWE|`.* Equivalently: the prior is
**active for a gene** only when **(0) `w > 0`**, AND (i) the gene's locus is present in the table with
at least two *distinct* effective frequencies among the gene's candidates (so `Δ ln P_HWE ≠ 0`), AND
(ii) the top-two candidates' coverage margin is `≤ w·|Δ ln P_HWE|`. Condition (0) is enforced as an
explicit `w == 0` short-circuit in `prior_active_for_gene`, `covered_prior_bonus`, and
`reconcile_zygosity_with_prior`, so at `w=0` the span collapses to 0 *and* activeness is forced false
— an exact coverage tie (`margin = 0`) can never satisfy `0 ≤ 0` to switch the prior on. (The `≤` in
(ii) is deliberate: once `w > 0`, breaking an exact tie toward the more frequent allele is the prior's
intended job — the byte-identity guarantee rests on the `w > 0` gate, not on a strict `<`.) When the
prior is inactive for a gene, the selector runs *exactly* the default logic — same objective, same
candidate set — for that gene. At `w=0`, or empty/constant table, the span is 0 and the prior is
inactive everywhere → byte-identical.
This is what makes the prior a strict, bounded refinement of the *existing* selector on both the
>2-type and ≤2-type paths, not a different selector.

Shared helper (new, `--allele-freq`-guarded, returns `0.0` when the flag is off OR when the prior is
inactive for this gene per the unifying principle — see empty-table guarantee):

```rust
fn prior_active_for_gene(&self, gene_idx, candidate_series: &[SeriesIdx]) -> bool
    // false unless locus present AND ≥2 distinct effective f among candidates
fn hwe_log_prior(&self, gene_idx, rank_j_series_idx, rank_k_series_idx) -> f64
    // ln2 + ln f_j + ln f_k  (het, j≠k),  or  2·ln f_j  (hom, j==k)
fn covered_prior_bonus(&self, gene_idx, series_j, series_k) -> f64
    // w · hwe_log_prior(...), or literal 0.0 if disabled/inactive
```

**Path A — pair search (>2 allele types),** `select_alleles_for_genes_haplotype_search`, the
objective comparison at **:4004–4016**. `covered_read_cnt` (accumulated :3962–4001) is the data term.
Add the prior to the objective *before* the comparison:

```rust
// after covered_read_cnt is finalized (~:4001), before the tie test at :4004:
let objective = covered_read_cnt + self.covered_prior_bonus(i, series_of(j), series_of(k));
```
and compare `objective`/`max_cover` (rename `max_cover` to track the objective). The
abundance-product tiebreak `abundance_sum` (:4008/4011) is retained *unchanged, as the lower-priority
tiebreak* — the prior sits strictly between coverage and the abundance tiebreak. `series_of(rank)` is
the major-allele series of the representative allele at that rank in `selected` (the same rank→series
resolution the loop already does at :3921–3936).

**Path B — ≤2-type greedy path.** Here zygosity is decided by the `filter_frac`/`major_abund`
thresholds in the greedy `filter`/rescue logic (:3634–3739), NOT by the coverage functional — so a
naive "score both candidates on the coverage functional and pick the argmax" would *not* reproduce
the greedy call at `w→0` (different objective), breaking the strict-refinement property and letting a
tiny `--allele-freq-weight` flip a clean zygosity call. The principled, *minimal* injection that does
not perturb the frozen filter logic: after `select_alleles_for_genes_haplotype_search` returns
(:3741), run a **hom-vs-het reconciliation** over each gene with `allele_type_cnt ∈ {1,2}`, **gated on
`prior_active_for_gene`** (so it never runs on an inactive gene). It computes the hom-vs-het
*coverage* margin `Δcov` with the same coverage functional used in Path A, and applies the unifying
no-override rule: **the prior may overturn the greedy zygosity only when `Δcov ≤ w·|Δ ln P_HWE|`
(within the prior span); if `Δcov` exceeds the span, the greedy decision stands unchanged.** At
`w=0` the span is 0, so Path B never overrides and reproduces the greedy call exactly — the byte-
identity anchor. This makes Path B a strict, bounded, coverage-vanishing refinement identical in
spirit and in guarantee to Path A. When the flag is off, or the gene is inactive, reconciliation is
not called at all.

### The homozygous candidate (well-defined at both paths)

The classic bug is `ln2 + ln f_A + ln f_B > 2 ln f_A` whenever `f_B > f_A/2`, i.e. HWE can rightly
prefer a het of two common alleles over a common hom — this is *correct* HWE and must be bounded by
`w`, not suppressed. To make hom-vs-het a well-defined comparison, synthesize the hom candidate with
the **same coverage functional**:

- **Hom candidate `G=(A,A)`**: covered set = the reads covered by allele A's ranks alone (evaluate
  the Path-A `covered_from_a` accumulation, :3865–3888, with `k` set equal to `j` so the second
  allele contributes no *new* reads). Coverage functional value = `covered_read_cnt` of that single
  set (no `k`-side missing-coverage term; set `k_missing_coverage = j_missing_coverage`). HWE term =
  `2·ln f_A`.
- **Het candidate `G=(A,B)`**: the existing `k>j` evaluation. HWE term = `ln2 + ln f_A + ln f_B`.

At Path A the hom candidate is enumerated by extending the inner loop to include `k == j` (currently
`(j+1)..allele_type_cnt`, :3890) **only when the prior is active for this gene** (unifying principle
above) — NOT merely when the flag is present. This is a blocker if gated on the flag alone: with an
empty table (`w`-term = 0) a true-hom's `(A,A)` candidate TIES the het `(A,B)` on `covered_read_cnt`
(B, the second allele of a real homozygote, contributes no *new* reads) and would then WIN the
abundance-product tiebreak (:4008), since `abundance_A² ≥ abundance_A·abundance_B` for the
higher-abundance rank A — so the argmax could change even though `covered_prior_bonus` is literal
`0.0`. The literal-`0.0` guarantee covers the objective *term* but not the *candidate set*; gating
enumeration on `prior_active_for_gene` closes that hole (empty/constant/absent-locus table ⇒ inactive
⇒ no `k==j` candidate ⇒ candidate set identical to default). When active, the hom candidate is scored
and, per the unifying no-override rule, can only win if the hom-vs-het coverage margin is within the
prior span. At Path B the two candidates are (i) hom on the single surviving series, (ii) het on the
top-two surviving series; both scored identically, and Path B runs *only* when the prior is active for
that gene. The only hom-vs-het asymmetry is the HWE term, so the comparison is honest; the `w`-bound
(next) prevents a spurious het against clear hom coverage.

### w-scale and the no-override bound (in real units)

`covered_read_cnt` is a **weighted read count** — a sum of `adjust_weight` floats (:3964), each O(1)
per optimal read, so a gene's coverage margin between two candidate pairs is realistically **tens to
hundreds** (a few reads' difference) up to **thousands** (deep panels). The HWE term
`Δ ln P_HWE = (ln f_j + ln f_k) − (ln f_j' + ln f_k')` is O(1) nats: with Jeffreys smoothing and
realistic HLA frequencies (common ~0.1–0.3, rare ~1e-3–1e-4, DB total ~1e4), a 100× single-allele
frequency gap is `ln 100 ≈ 4.6` nats; the floor-to-common extreme is ~`2·ln(0.3/1e-4) ≈ 16` nats.

**Set `w = 2.0`.** Then the prior span `w·|Δ ln P_HWE|` is **≤ ~10 nats for a realistic single-allele
frequency gap, ≤ ~32 nats at the pathological floor-vs-common extreme**, expressed in the *same
weighted-read units* as `covered_read_cnt`. Stated guarantee: **the prior cannot flip a call whose
coverage margin exceeds the prior span**; concretely, a coverage margin of ≳ 32 weighted reads (≈ a
handful of reads at typical `adjust_weight`) can *never* be overridden by any frequency table. Below
that margin — exact ties, hom-vs-het boundaries, ~1-read differences — the prior *may* tip the call
(the intended behavior). `w` is a documented CLI knob (`--allele-freq-weight`, default 2.0) so the
span is auditable; the deferred null-penalty (above) is to be capped at the same span. As coverage → ∞ the
`covered_read_cnt` margin grows without bound while the HWE term stays O(1), so the call is
prior-independent in the limit (coverage-vanishing, tested).

### Empty-table / skip-when-constant guarantee (byte-identity with the flag present)

Byte-identity with the flag present requires TWO invariants, because the prior touches both the
objective *value* and the candidate *set*:

- **Objective invariant.** Adding a constant `w·lnP` to *both* sides of an exact-`f64` tie (:4005)
  can perturb the `==` via round-off. The guarantee: **`covered_prior_bonus` returns the exact literal
  `0.0` (not a computed near-zero) whenever the prior is inactive for the gene** — (a) `--allele-freq`
  absent, OR (b) the gene's locus absent from the table, OR (c) all compared series map to the
  identical effective `f` (e.g. all absent → same floor, so `Δ ln P_HWE = 0`). Then the objective is
  `covered_read_cnt + 0.0`, and since `x + 0.0 == x` for all finite non-NaN `x`, the exact `==` tie at
  :4005 and the whole comparison chain are bit-for-bit unchanged. Enforced by making the helper's
  first action a `prior_active_for_gene` check that early-returns the literal `0.0`; a unit test
  asserts the returned value's bits equal `0.0f64.to_bits()` on an empty table.
- **Candidate-set invariant.** The hom-candidate loop extension (`k==j`, Path A) and the entire
  Path-B reconciliation are gated on `prior_active_for_gene`, NOT on the flag. So when the prior is
  inactive the enumerated candidate set is *identical* to the default path — no phantom `(A,A)`
  candidate can win a tiebreak (the abundance-product hole in Refinement above). Term-zero alone is
  insufficient; this gate is what makes empty/constant/absent-locus tables provably byte-identical
  rather than merely "approximately."

Together these mean: flag present + empty table, and flag present + a table missing the called locus,
both reproduce the flag-off path exactly — objective, candidate set, and argmax.

### Default-off path

`--allele-freq` absent (or, per the unifying principle, the gene inactive) → `covered_prior_bonus`
returns literal `0.0`, the Path-A loop bound stays `(j+1)..allele_type_cnt` (no hom candidate),
Path-B reconciliation is not invoked → both selectors run verbatim on the default candidate set →
byte-identical. The abundance-product tiebreak ordering (:4008) is untouched.

### Modeling caveats (documented)

- HWE prior assumes per-locus allele independence; HLA loci are in strong LD and cohorts are
  admixed, so the prior is mis-specified for admixed samples — acceptable for a bounded opt-in
  tie-breaker, stated in user docs. Population-conditioning is future work.
- T1K's genotype/discriminative quality (Poisson/`alnorm`, `select_alleles_for_genes_quality_scores`)
  is a data-only score and is **not recalibrated** when the prior flips a call. Documented; the
  #33 metrics panel surfaces the prior's contribution so a prior-driven flip is visible.

---

## #33 — identifiability measure (report-only)

**Honest constraint discovered in the code.** The kallisto/RSEM identifiability notion is over the
data-driven EC partition (`equivalent_class_to_alleles` :2029; within one EC all members share the
same read set *by construction* — `is_assigned_read_the_same` :2553). Two consequences the prior
review's spec did not resolve, resolved here:

1. **`distinguishing_reads` as "reads separating an allele from its EC-mates" is degenerate** — EC
   mates share reads by definition, so intra-EC distinguishing reads are always 0. A *cross-EC*
   distinguishing count would need per-read allele-compatibility retained past selection, which it is
   not. **Dropped.**
2. **G-group is NOT retained** — the core model has no G-group annotation; only the major-allele
   *series* (`major_allele_idx`, `major_allele_idx_to_name` :1972). A "G-group compatible_set" would
   require threading new reference-annotation data through the whole pipeline. **Not worth it for a
   report-only column; dropped/deferred.** We report the *series* set instead, which IS retained and
   is an honest sample-independent resolution annotation.

`remove_low_mapq_allele_in_equivalent_class` / `remove_low_likelihood...` (:3460) *prune* EC
membership before selection, so all measures are computed on the **retained, post-prune EC state** —
this is the state the call was actually made from, so it is the correct basis for an identifiability
signal.

### What #33 ships (all computable from retained data, all report-only)

Computed post-selection in `write_metrics_tsv` (`crates/unum/src/stages/genotype.rs:572`), added as
new columns to the existing header (:582). Per representative called allele (the reps already picked
at :602–633):

- **`ec_set_size`** — number of member alleles in the called allele's post-prune equivalence class
  (`equivalent_class_to_alleles[allele_info[idx].equivalent_class].len()`). 1 ⇒ the allele is
  read-distinguishable in this sample; >1 ⇒ non-identifiable within the EC (the disciplined
  kallisto-style statement).
- **`ec_ambiguity_entropy`** — Shannon entropy `H = −Σ p_i ln p_i` over the EC members, with
  `p_i = ec_abundance_i / Σ ec_abundance` (a proper distribution over the sample-indistinguishable
  set; Paunić 2012 framing). Singleton → 0; uniform over k → `ln k`. Uses retained `ec_abundance`.
- **`identifiability`** — `1 / ec_set_size` ∈ (0,1]; a monotone, unit-free distinguishability
  fraction (1 = uniquely identified). Simple, honest, derived from the same retained set — replaces
  the degenerate distinguishing-reads ratio.
- **`series_set`** — the major-allele-*series* string(s) of the EC members (from
  `major_allele_idx_to_name`), semicolon-joined, plus their count. This is the "compatible set" at
  the resolution unum actually retains (series, not G-group), reported as an annotation distinct
  from the sample-data EC measures above.

**Deferred (not shipped):** cross-EC distinguishing reads, and G-group compatible sets — both need
extra retained data (per-read allele compatibility, and a G-group reference annotation respectively)
whose cost is not justified for report-only columns. Revisit if a G-group annotation is ever added to
the reference model for another reason.

### Output

Extends the `_metrics.tsv` header (:582) with `ec_set_size`, `identifiability`,
`ec_ambiguity_entropy`, `series_set`. Same `--emit-metrics` flag, same file, one row per rep as
today. No new flag; no change to any frozen file.

---

## Testing strategy

**Byte-identity gates (both features):**
- `cargo ci-test-sys` green with `--allele-freq` **absent** — proves the default path.
- `cargo ci-test-sys` green with `--allele-freq` **present + empty table** and **present + a table
  missing every called locus** — proves the skip-when-inactive guarantee on the flag-present path.
  Both must reproduce the flag-off path *exactly, candidate-set included* (not just equal outputs):
  assert no `k==j` hom candidate is enumerated and Path-B reconciliation never runs when inactive
  (e.g. a debug-instrumented candidate-set/enumeration count matches the flag-off run).
- `_metrics.tsv` columns present under `--emit-metrics` leave `_genotype.tsv`/`_allele.tsv`
  byte-identical (existing `crates/unum/tests/emit_metrics.rs` extended).

**#29 unit tests (in `genotyper.rs`, in-memory, no I/O):**
- freq loader: Jeffreys smoothing, per-locus normalization, resolution-aware series matching
  (coarse-DB↔fine-call rollup, fine-DB↔coarse-call sum), missing-allele floor, missing-locus no-op.
- `hwe_log_prior`: hom `2 ln f` vs het `ln2 + ln f + ln f` forms; the `f_B > f_A/2` crossover.
- `prior_active_for_gene`: false on empty table / absent locus / all-candidates-same-`f`; true only
  with ≥2 distinct effective frequencies among candidates. `covered_prior_bonus` returns bit-exact
  `0.0` whenever inactive (`.to_bits()` assertion).
- **Candidate-set invariant:** with an inactive gene, the Path-A `k==j` hom candidate is NOT
  enumerated and Path-B does NOT run — the enumerated candidate set equals the flag-off set
  (guards the abundance-product-tiebreak hole).
- **Path-B `w=0` reproduces greedy:** a ≤2-type gene whose greedy zygosity is decided by
  `filter_frac`/`major_abund` is left unchanged at `w=0`, and at `w>0` is overturned only when the
  hom-vs-het coverage margin is within `w·|Δ ln P_HWE|`.
- Five adversarial safeguard properties on synthetic `covered_read_cnt` + frequencies (uniform across
  both paths): (a) exact-tie coverage → higher-HWE candidate wins; (b) coverage margin > prior span
  (`w·Δlnf`) → prior cannot override; (c) coverage → ∞ (scale the read counts) → prior influence → 0;
  (d) single-observed-allele locus → hom forced, het branch never spuriously entered;
  (e) hom-vs-het adversarial (common-hom truth vs common/rare-het) at realistic coverage → prior does
  NOT manufacture the het.

**#29 functional:** a constructed near-tie fixture where a frequency table flips the call to the
common allele *only* with `--allele-freq`, and a clear-coverage-margin fixture that does not flip.

**#33 unit tests:** `ec_set_size`/`identifiability` on synthetic ECs, where `identifiability` is
defined as `1 / ec_set_size` (singleton → `ec_set_size` 1, `identifiability` 1/1.0 = 1.0; k-member →
`ec_set_size` k, `identifiability` 1/k); `ec_ambiguity_entropy` (singleton → 0; uniform-k → `ln k`;
skewed abundances < `ln k`);
`series_set` joins the member series independent of EC entropy.

## Implementation staging (two stacked PRs, on the current chain tip)

1. **#29** (`29/nh/frequency-prior`): freq loader + `--allele-freq`/`--allele-freq-weight`
   + `hwe_log_prior`/`covered_prior_bonus` + Path-A objective injection (:4004) + hom candidate
   (`k==j` guarded, collapses to 1-type on a hom win) + Path-B reconciliation. Larger; touches
   selection. (`--allele-freq-null-penalty` deferred — see the null/expression note above.)
2. **#33** (stacked on #29): four identifiability columns in `write_metrics_tsv`. Report-only; small.

## Citations

Polysolver (Shukla 2015); arcasHLA (Orenbuch 2020); T1K (Song 2023); HLA\*LA (Dilthey 2019);
Paunić 2012; CIWD 3.0 (Hurley 2020); AFND (Gonzalez-Galarza); Wigginton 2005 (HWE);
kallisto/RSEM/Salmon (Bray 2016, Li 2011, Patro 2017).
