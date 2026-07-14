# Divergences from T1K

`unum` began life as a byte-identical Rust port of [T1K](https://github.com/mourisl/T1K), validated by differential-testing against the original C++. Owning the port means we can fix T1K's latent bugs and improve its behavior where `unum` can do better. This document tracks every place `unum` **intentionally departs** from T1K, and — separately — the known T1K bugs/quirks `unum` still **reproduces on purpose** for parity (candidates for future divergence).

Each entry cites the T1K source it relates to (paths are relative to the [T1K repo](https://github.com/mourisl/T1K)) and the `unum` code that implements the divergence.

## Divergences (unum intentionally differs from T1K)

### Weighted pileup counts in the variant caller

- **unum:** `crates/unum-core/src/variant_caller.rs` — `VariantCaller::update_base_variant_from_overlap`
- **T1K:** `VariantCaller.hpp:147`
- **Affects real output:** yes.

T1K increments the per-base pileup `count[nucIdx] += 1` unconditionally, ignoring the `weight` argument. As a result `count` is always identical to the separate `unweightedCount` array (which also does `+= 1`), the align-info-only pass (`updateType == 1`, `weight = 0`) still contributes full support, and a fractional multi-allele assignment (`weight < 1`) counts as a whole read. The three parallel arrays (`count`, `uniqCount`, `unweightedCount`) and the field name `count` all indicate the author *intended* `count` to be the weighted tally — this is a latent T1K bug.

`unum` accumulates `bv.count[nuc_idx] += weight`, so `count` is the genuinely weighted support its name implies: the align-info pass adds nothing, and fractional assignments contribute their fraction. This changes which candidate variants clear the absolute `count >= 5` threshold in `find_candidate_variants` and the support values reported in the allele VCF.

### Non-ACGT reference bases skipped in candidate-variant discovery

- **unum:** `crates/unum-core/src/variant_caller.rs` — `VariantCaller::find_candidate_variants`
- **T1K:** `VariantCaller.hpp:323`
- **Affects real output:** only where the reference consensus contains non-ACGT bases (e.g. `N`), which real references do.

T1K's `nucToNum` table maps non-ACGT bases (such as `N`) to `-1` and then indexes `baseVariants[i][j].count[nucToNum[s[j] - 'A']]` — i.e. `count[-1]`, an out-of-bounds read (undefined behavior). Reference consensuses do contain `N`s, so a naive port would read undefined memory (or panic).

`unum` skips positions whose reference base is not A/C/G/T: no variant is ever called against a non-ACGT reference base. This is well-defined and avoids reading out of bounds.

### Trailing-insertion out-of-bounds read guarded in the fragment/variant graph

- **unum:** `crates/unum-core/src/variant_caller.rs` — `build_fragment_candidate_var_graph`
- **T1K:** `VariantCaller.hpp` (the `baseVariants[seqIdx][refPos]` walk over an overlap's alignment)
- **Affects real output:** no (behavior-equivalent on real data).

A trailing `EDIT_INSERT` in an overlap's alignment can advance `refPos` one position past the last consensus base. T1K then reads `baseVariants[seqIdx][refPos]` out of bounds (undefined behavior that, in practice, does not affect substitution counts because a trailing insertion contributes no base variant). This shape is confirmed reachable on the real `unum analyze` pipeline.

`unum` stops the alignment walk when `ref_pos` reaches the end of the consensus rather than reading out of bounds. Because a trailing insertion contributes no base variant, this is behavior-equivalent to T1K's (accidental) outcome while being memory-safe.

### T1K's null-pointer "shadowing" bug not reproduced

- **unum:** `crates/unum-core/src/variant_caller.rs` — `VariantCaller::update_base_variant_from_overlap` (the `align == None` recompute path)
- **T1K:** `VariantCaller.hpp` (the `if (align == NULL)` recompute branch)
- **Affects real output:** no (the branch is dead code on the real end-to-end path).

When an overlap arrives without a precomputed alignment, T1K declares a fresh `signed char *align = new signed char[...]` *inside* the `if (align == NULL)` block, shadowing the outer `align`. The recomputed alignment is never assigned back to the outer pointer, which stays `NULL`, so the subsequent loop dereferences a null pointer. This branch is dead code on the real pipeline (`align` is always populated there), so the bug never fires in practice.

`unum` uses the recomputed alignment directly and does not reproduce the shadowing bug. This only affects the test-only `align == None` path.

### Opt-in HWE frequency prior adds a homozygous pair-search candidate (#29)

- **unum:** `crates/unum-core/src/genotyper.rs` — `select_alleles_for_genes_haplotype_search` (the Path-A pair search) and `reconcile_zygosity_with_prior` (Path B)
- **T1K:** `Genotyper.hpp:1697-1996` (pair search enumerates only het pairs `k > j`)
- **Affects real output:** only when `--allele-freq` is supplied with a weight > 0 and the gene's locus is present with ≥2 distinct candidate frequencies; **off by default → byte-identical to T1K.**

With the opt-in population-frequency prior active for a gene, unum extends the Path-A pair search to also enumerate the **homozygous candidate `(A, A)`** (`k == j`), which T1K never considers (its inner loop is `k > j`, het pairs only). When that homozygous candidate wins the (coverage + HWE-prior) objective, unum **collapses the gene to a genuine 1-type homozygous call** — it drops the losing allele types and keeps only the winning type at rank 0 — so `get_gene_allele_types` returns 1 and the gene is treated identically to a naturally-homozygous gene. Without the collapse the winner remap would leave rank-1 in place, the gene would stay multi-typed, and the reconciliation loop would spin to `ITER_MAX`; the collapse mirrors Path B's `retain(rank != 1)` hom reduction. Covered by a dedicated regression test (`path_a_hom_winner_collapses_gene_to_one_type`).

### Opt-in null/expression-variant selection penalty (`--allele-freq-null-penalty`)

- **unum:** `crates/unum-core/src/genotyper.rs` — `is_null_expression_name`, `AlleleInfo::is_null_expression`, `Genotyper::null_penalty_term`, and the Path-A objective / Path-B margin
- **T1K:** no equivalent — T1K's selector has no notion of allele expression status
- **Affects real output:** only when `--allele-freq-null-penalty > 0`; **off by default (`0.0`) → byte-identical to T1K.**

A null / low-expression allele (IMGT `N`/`L`/`S`/`C`/`A`/`Q` suffix) typically differs from a common expressed allele by a single defining variant, so a caller can flip to the non-expressed allele on a handful of reads. `--allele-freq-null-penalty p` subtracts a fixed, coverage-independent penalty (`p` per null haplotype a candidate would call; `2p` for a homozygous-null call) from the selection objective in **both** Path A (the pair search) and Path B (hom-vs-het reconciliation), biasing the caller away from asserting a non-expressed allele on marginal evidence. It is **name-driven** — detection reads the full reference allele name (not the frequency table), so unlike the HWE term it acts even without `--allele-freq`. When the penalty is active for a gene (`p > 0` and the gene has a null/expression rank) it also **enables the `k == j` homozygous candidate (Path A) and runs Path B**, so the penalty can drive a null-het → expressed-hom call *without any frequency table* (the HWE term simply drops to 0). All of this activation is gated on `p > 0`, so at the default `p == 0` the candidate set and the whole selection are **byte-identical** to the oracle. The value is bounded to `[0, 16]` (`Genotyper::NULL_PENALTY_MAX`, half the ~32-weighted-read no-override span, since the max per-candidate penalty is `2p`) so the penalty alone can never flip a coverage margin beyond the span. Covered by detection, per-haplotype, Path-A demote/no-override, Path-B override, and no-table scope-full tests.

### `unum combine` treats `-q` as numeric (`t1k-merge.py` crashes on nonzero `-q`)

- **unum:** `crates/unum-core/src/combine.rs` / `crates/unum/src/cli.rs` — `CombineArgs::min_quality` is an `f64`.
- **T1K:** `t1k-merge.py:35`.
- **Affects real output:** only when `-q` is passed a nonzero value; the default (`-q 0`) is byte-identical.

`t1k-merge.py` compares `float(cols[i + 2]) > args.qual`, but `argparse` leaves a passed `-q` value as a **string** (only the unset default is the int `0`). So any `python3 t1k-merge.py -q 20 ...` raises `TypeError: '>' not supported between instances of 'float' and 'str'` and produces no output — the `-q` flag is effectively unusable upstream.

`unum combine` parses `-q` as `f64` and filters `quality > min_quality` as intended, so nonzero `-q` works. On every flag combination T1K *can* run (default, `-n`, `--total-quality`/`--tq`), `unum combine` is byte-identical to `t1k-merge.py` (validated on 15 1000G-DNA and 40 HPRC genotype outputs).

### `unum copy-number` treats `-q` as numeric (`t1k-copynumber.py` crashes on nonzero `-q`)

- **unum:** `crates/unum-core/src/copy_number.rs` / `crates/unum/src/cli.rs` — `CopyNumberArgs::min_quality` is an `i64`.
- **T1K:** `t1k-copynumber.py:60`.
- **Affects real output:** only when `-q` is passed a nonzero value; the default (`-q 0`) is byte-identical.

`t1k-copynumber.py` compares `quality <= args.qual` with `quality = int(...)`, but `argparse` leaves a passed `-q` value as a **string** (only the unset default is the int `0`). So `python3 t1k-copynumber.py -q 20 ...` raises `TypeError: '<=' not supported between instances of 'int' and 'str'` and produces no output — the same latent bug as `t1k-merge.py`.

`unum copy-number` parses `-q` as `i64` and filters `quality > min_quality` as intended. On every flag combination T1K *can* run (default, `--nomissing`, `--upper-quantile`/`--lower-quantile`, `--adjust-var`), `unum copy-number` is byte-identical to `t1k-copynumber.py` (validated on the T1K KIR example genotype).

### `unum copy-number` errors when the anchor sample can't yield a one-copy distribution (`t1k-copynumber.py` divides by zero / NaN)

- **unum:** `crates/unum-core/src/copy_number.rs` — `infer_copy_numbers` returns an error whenever at least one allele still needs a copy call but the anchor sample has non-positive variance: it is **empty** (no usable anchors), has a **single** value, or is **all-equal** (every anchor abundance identical). All three give `sigma = sqrt(var) <= 0`, so the one-copy `N(mean, var)` is unusable.
- **T1K:** `t1k-copynumber.py:95` (`mean = sum(abundances)/inspectAlleleCnt`) for the empty case; the per-allele likelihood's `sigma = sqrt(var)` division for the single / all-equal (`var == 0`) case.
- **Affects real output:** only on a degenerate anchor sample — no `--nomissing` gene present *and* every passing gene homozygous (empty), or the quantile band / `--nomissing` set collapsing to a single anchor or a run of identical abundances (`var == 0`); never on the KIR example, nor on any input with two or more heterozygous non-`--nomissing` alleles of differing abundance.

With an empty anchor sample, `t1k-copynumber.py` evaluates `sum(abundances)/len(abundances)` on an empty list and raises `ZeroDivisionError`; with a single or all-equal sample it computes `var == 0`, so `sigma = sqrt(var) == 0` and the per-allele likelihood divides by zero — both produce no usable output. `unum copy-number` cannot estimate the one-copy `N(mean, var)` in any of these cases either, but fails cleanly with a descriptive error rather than dividing to `NaN` and emitting garbage copy/ratio calls. When *no* allele passes the quality filter (nothing to call at all), unum instead emits each gene row with sentinel slots (`. -1 0`), since no copy call — and thus no `NaN` — is ever produced.

### `unum copy-number` validates its numeric flags (`t1k-copynumber.py` does not)

- **unum:** `crates/unum-core/src/copy_number.rs` — `CopyNumberConfig::validate` (called from the stage before inference) requires finite `--upper-quantile`/`--lower-quantile` in `[0, 1]` with `--lower-quantile <= --upper-quantile`, and a finite `--adjust-var > 0`.
- **T1K:** `t1k-copynumber.py` — `argparse` parses these as raw `float`s with no range check.
- **Affects real output:** only on invalid numeric flag values, including out-of-range/non-finite values and reversed quantile bounds. On every in-range combination (the documented usage) unum stays byte-identical to `t1k-copynumber.py`.

On an invalid flag, `t1k-copynumber.py` runs on regardless: an out-of-`[0, 1]` quantile or reversed bounds (`--lower-quantile > --upper-quantile`, even with both individually in `[0, 1]`) selects an empty, truncated, or wrapped anchor band from its unchecked Python slice, and a non-positive `--adjust-var` drives `var <= 0` so the one-copy `sigma` collapses to `0`/`NaN` — either way it emits a (possibly garbage) table instead of erroring. `unum copy-number` rejects these up front with an actionable argument error rather than reproducing that silent-garbage behavior. This is a deliberate divergence confined to invalid input; it never changes output for a run T1K would have produced meaningfully.

## Known T1K quirks/bugs still reproduced (candidates for future divergence)

These are places where `unum` currently matches T1K on purpose — a known bug, quirk, or UB in T1K that we reproduce for parity rather than fix. They are tracked here as candidates for future divergence: if any starts to matter for a real result, revisit it the way we did the weighted-count bug above. Each is marked in the source with a comment explaining it is a deliberate parity choice.

### Can affect real output

| `unum` source | T1K source | Quirk / bug reproduced |
| --- | --- | --- |
| `crates/unum-core/src/refbuild/dat.rs` (`find_mode`) | `ParseDatFile.pl` | Among keys tied for the max count, the winner is chosen by **string** comparison, so `"9"` beats `"10"` (`'9' > '1'`) even though `10 > 9` numerically. Can change chosen bytes on a genuine multi-way tie. |
| `crates/unum-core/src/refbuild/dat.rs` (5'UTR pick) | `ParseDatFile.pl` `substr($seq,0,$end)` | Perl `substr` off-by-one: the comparison uses `$end` bytes, one short of the full available range, when picking the best 5'UTR. Affects 5'UTR selection on ties/boundaries. |
| `crates/unum-core/src/align_algo.rs` (row-0 `e[]` init) | `AlignAlgo.hpp:256-266` | Row-0 gap-score init reads the **leftover loop variable `i`** from the preceding loop instead of the current `j` (`e[j] = SCORE_GAPOPEN + i*SCORE_GAPOPEN`). Affects alignment scores/traceback. Covered by a dedicated regression test. |
| `crates/unum-core/src/genotyper.rs` (read-assignment merge) | `Genotyper.hpp:893-894` | Merge compares `.end < .end` but then assigns `.end = ...start` (not `.end`) — a real stock quirk. Affects read-assignment `.end` values. |
| `crates/unum-core/src/genotyper.rs` (`GetGeneAlleleTypes`) | `Genotyper.hpp` | When a gene has more than one secondary allele-type group, only the group with the **highest** `type` value survives in the genotype string, silently dropping the others. |
| `crates/unum-core/src/genotyper.rs` (`passed` from similarity) | stock `SeqSet` | `passed` is computed from `extendedOverlap.similarity` **before** clip handling mutates that field, and never revisited. Affects which extended overlaps pass. |
| `crates/unum-core/src/kmer_index.rs` | `KmerIndex.hpp:115,123` | `prev_kmer_code` starts at `code == 0`, so a read whose **first full k-mer window is all-`A`** (`code == 0`) is silently dropped from the index. |
| `crates/unum-core/src/extract.rs` | `FastqExtractor.cpp:271-418` | Extra trailing mate-2 records (mate-2 file longer than mate-1) are silently ignored and the tool exits 0; the symmetric trailing check is deliberately omitted to match stock. |

### Latent, cosmetic, or pathological-input-only (reproduced but currently inert)

Reproduced verbatim for parity; each is either provably inert on the only reachable path or requires input that never occurs in the real pipeline. Listed briefly so they are not re-discovered as "new" bugs:

- `crates/unum-core/src/overlap.rs` — index-reuse across a re-sorted `hits` array (`SeqSet.hpp:1408-1424`); inert for its single caller (`filter == 0`, `i == 0`).
- `crates/unum-core/src/kmer.rs` — `Append` treats only literal `'N'` as invalid while `Prepend` uses the full `nucToNum` lookup; only matters for non-`N` non-ACGT bytes.
- `crates/unum-core/src/kmer_count.rs` — `jaccard_similarity` reproduces stock's `i32`-overflow arithmetic bit-for-bit (only on unexpected overflow).
- `crates/unum-core/src/alignments.rs` — CIGAR length switch adds reference length for any unlisted op, mirroring stock's `switch` (only exotic CIGAR ops).
- `crates/unum-core/src/bam_extract.rs` — coord-FASTA parse assumes each sequence is on a single unwrapped line (`fscanf %s`); a wrapped sequence would desync stock too.
- `crates/unum-core/src/genotyper.rs` (`parse_allele_name`) — with a non-default explicit `alleleDigitUnits`, name parsing walks raw bytes instead of `:`-delimited fields (not the pipeline default).
- `crates/unum/src/stages/copy_number.rs` — `--nomissing` entries are split on `,` **without** trimming whitespace, matching `t1k-copynumber.py:42` (`{g:1 for g in ...split(",")}`). A value with surrounding spaces (`--nomissing "A, B"`) never matches a whitespace-split gene name in either tool, so the whitespace-prefixed `B` anchor is silently dropped identically; unspaced lists (the documented usage) are unaffected.

### Structural numeric tolerance (not a discrete bug)

`unum`'s abundance/likelihood arithmetic targets the `-O3` C++ oracle within a tight **relative** tolerance (~`1e-6`), not exact `f64` bits, because the C++ oracle may FMA-contract multiply-add chains. This is an acknowledged, structural non-bit-identity across the genotyper rather than a single fixable site; see the module docs in `crates/unum-core/src/genotyper.rs`.

---

**Related work:** open PR #44 introduces a `docs/KNOWN-DIVERGENCES.md` documenting a floating-point/abundance-cutoff divergence ("Cluster A": T1K's absolute `ec_abundance <= 1e-6` cutoff in `select_alleles_for_genes`, `Genotyper.hpp:1519`, flips deep-resolution allele calls on 5 HPRC samples). That doc and this one should be consolidated into a single canonical divergence log once both land.
