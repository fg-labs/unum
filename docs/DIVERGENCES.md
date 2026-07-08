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

### Structural numeric tolerance (not a discrete bug)

`unum`'s abundance/likelihood arithmetic targets the `-O3` C++ oracle within a tight **relative** tolerance (~`1e-6`), not exact `f64` bits, because the C++ oracle may FMA-contract multiply-add chains. This is an acknowledged, structural non-bit-identity across the genotyper rather than a single fixable site; see the module docs in `crates/unum-core/src/genotyper.rs`.

---

**Related work:** open PR #44 introduces a `docs/KNOWN-DIVERGENCES.md` documenting a floating-point/abundance-cutoff divergence ("Cluster A": T1K's absolute `ec_abundance <= 1e-6` cutoff in `select_alleles_for_genes`, `Genotyper.hpp:1519`, flips deep-resolution allele calls on 5 HPRC samples). That doc and this one should be consolidated into a single canonical divergence log once both land.
