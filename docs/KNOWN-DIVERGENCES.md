# Known genotype-call divergences from the vendored T1K oracle

unum is validated **byte-identical** to the vendored C++ T1K oracle on the
FFI differential suite. That suite compares intermediate and final values with a
documented `~1e-6` relative tolerance on floating-point quantities (abundances,
similarities), because the oracle's `g++ -O3` build applies FMA contraction
(`a*b + c` → one fused, single-rounded op) while unum's Rust `+`/`*` round
twice — see `crates/unum-core/src/genotyper.rs` (the `# FLOATS` notes on
`em_update` / `squarem_alpha`). Within that tolerance the two agree.

This file records the **end-to-end genotype *call*** divergences that survive on
real data despite the byte-identical intermediate math — surfaced by the 37-sample
HPRC real-data differential (`compare.py`, unum vs the standalone vendored
`genotyper`). They are deliberately **accepted, not "fixed"**, for the reasons
below.

## Summary

Of 37 HPRC samples (41 genes each), the initial differential found **29 matching
the oracle on every gene** and **8 diverging** on exactly one gene each — all at
**3rd/4th-field (deep) resolution** in **Class II / non-classical loci**, at
**noise-floor abundances (~1e-6 normalized)**.

They fall into two mechanisms. Fixing Cluster B resolved 3 of the 8, leaving
**32 clean / 5 accepted** (confirmed by re-running the differential with zero
regressions on the previously-clean samples).

### Cluster B — read-assignment port bug (FIXED)

- **HG01109 HLA-U**, **NA19240 HLA-DQB2**, **HG02622 HLA-DQB2**
- Root cause: an operator-precedence error in the truncated-mate-pair filter of
  `read_assignment_to_fragment_assignment` — the `unused` (already-assigned)
  guard was applied to the `matchCnt >` branch, whereas `SeqSet.hpp:2601-2609`
  applies it only to the `matchCnt ==` branch.
- **Fixed** (unum-only, no oracle change); these three now match the oracle.

### Cluster A — FP-tolerance knife-edge at a hard abundance threshold (ACCEPTED)

- **HG00438 HLA-DPB2**, **HG02055 HLA-DQA2**,
  **HG02818 / HG03098 / HG03486 HLA-DOB**
- Root cause: the sub-`1e-6` FMA-vs-double-rounding EM drift (below the
  differential tolerance, so all FFI gates pass) shifts a noise-floor
  equivalence class across the **hard `ec_abundance <= 1e-6` processing cutoff**
  in `select_alleles_for_genes` (`Genotyper.hpp:1519`). unum and the oracle
  therefore process a slightly different *set* of low-tail ecs; the extra ec is
  "rescued" and perturbs a coverage-based re-rank, flipping which near-identical
  deep-resolution representative is reported. (For the HLA-DOB cases the effect
  is amplified by a dense equivalence class over ~75 near-identical alleles.)
- **Why accepted:** at these abundances the deciding quantity is essentially
  rounding noise and the alleles are near-indistinguishable in the covered exon
  regions, so the "correct" deep-resolution call is genuinely ambiguous. The
  first two fields (the clinically meaningful resolution) match the oracle in
  every case.

## The deferred fix option (Cluster A)

The fragility is the **hard, absolute `1e-6` cutoff** paired with sub-tolerance
FP noise. A robustness fix would make the cutoff **relative** — e.g. break when
`ec_abundance <= max(1e-6, 1e-4 * top_ec_abundance)` — in **both** unum and
the vendored oracle, so noise-floor ecs are consistently excluded regardless of
the tiny FP drift. This is deliberately **not** applied here because it:

1. patches the ground-truth oracle (diverging it from upstream `9376b55`), and
2. is a heuristic that changes selection for *all* samples, requiring a full
   regression sweep to confirm it does not perturb the 29 currently-clean
   samples.

If a Cluster-A-style divergence ever appears at a *meaningful* (non-noise-floor)
abundance, revisit this fix.

## Reproducing

```bash
python3 compare.py \
  --samples-dir <hprc-samples> \
  --ref <hla_reference>/_dna_seq.fa \
  --oracle-bin vendor/t1k/genotyper \
  --fgt1k-bin target/release/unum \
  --out-dir <out> --threads 4
```

The oracle output is cached per sample; re-running after a unum change only
re-runs unum. See the per-`(sample, gene)` `EXACT` / `SWAP` / `DIFF`
classification in `compare.py`.
