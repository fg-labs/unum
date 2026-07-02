# Reference-build golden fixtures: provenance and pins (Phase 1, Task 1.1)

This directory pins small, hand-picked KIR inputs and the byte-identical golden
output of the Perl oracle (`vendor/t1k/t1k-build.pl`) that a future Rust port of
`t1k-build` must reproduce exactly. No Rust code is involved in this task —
this is data pinning plus running the existing vendored Perl scripts.

## Source data

- **Repo**: `mourisl/T1K_manuscript_evaluation` (GitHub)
- **Commit that added the file**: `c95ddcd77be0c317763fedf09f9cd20ea3aa0a23`
  ("Add the reference sequences used in the manuscript", 2022-10-18)
- **Tarball**: `kiridx_2_10_0.tar.gz` (IPD-KIR v2.10.0), extracts to
  `kiridx/{KIR.dat, kir_dna_seq.fa, kir_rna_seq.fa, kir_dna_coord.fa, kir_rna_coord.fa}`
- **Local copy used**: `/tmp/t1k_eval/kiridx_2_10_0.tar.gz`, obtained from a
  prior Phase-0 clone of the repo. That local clone's git history is shallow/
  partially unreadable (`git log` fails to traverse parents), so provenance was
  re-verified independently: `git hash-object kiridx_2_10_0.tar.gz` locally
  gives `fdea8b794445db3f7fabdfc180322d7366fac8a3`, which matches the blob SHA
  returned by the GitHub Contents API for this exact path at commit
  `c95ddcd77be0c317763fedf09f9cd20ea3aa0a23` (size 2646869 bytes). The tarball
  content is therefore confirmed byte-identical to the one added in that
  commit, independent of the local clone's current HEAD (which has since moved
  to a later commit, `8a57d90a948aa9a4c80919153e3710cd37df334f`, 2024-03-13,
  that still carries the same file unmodified).
- `kiridx/KIR.dat` is 21,271,966 bytes, 1531 alleles (1531 `//` record
  terminators), spanning 17 KIR genes: KIR2DL1, KIR2DL2, KIR2DL3, KIR2DL4,
  KIR2DL5A, KIR2DL5B, KIR2DP1, KIR2DS1, KIR2DS2, KIR2DS3, KIR2DS4, KIR2DS5,
  KIR3DL1, KIR3DL2, KIR3DL3, KIR3DP1, KIR3DS1.

## Subset: `kir_subset.dat` (147,967 bytes)

Extracted verbatim (byte-exact record slices, each terminated by its original
`//\n`) from `kiridx/KIR.dat` using a script that splits on `//\n` boundaries
and looks up each wanted allele by its `allele="..."` FT qualifier. 8 alleles
across 4 genes:

| Allele | Gene | Exons | Introns | Partial? | Why chosen |
|---|---|---|---|---|---|
| `KIR2DL1*0010101` | KIR2DL1 | 9 | 8 | no | full-length genomic representative; also has a `/pseudo` exon-qualifier feature internally (exercises the pseudo-exon-trim code path in `ParseDatFile.pl`'s `FT ... pseudo$` handling even though the allele itself is not a pseudogene) |
| `KIR2DL1*0020101` | KIR2DL1 | 9 | 8 | no | second full-length allele for the same gene, to exercise mode/consensus logic (`FindMode`, exon-length mode, intron consensus) across >1 allele per gene |
| `KIR2DL4*00101` | KIR2DL4 | 8 | 0 | no | shorter/no-intron genomic record for KIR2DL4 (deliberately picked alongside a full intron-bearing record for the same gene to exercise `geneExonCntMode`/last-exon-length trimming across heterogeneous record shapes) |
| `KIR2DL4*0010201` | KIR2DL4 | 8 | 7 | no | full intron-bearing KIR2DL4 record |
| `KIR3DL1*0010101` | KIR3DL1 | 9 | 8 | no | full-length representative |
| `KIR3DL1*0010102` | KIR3DL1 | 9 | 8 | no | second full-length allele for the same gene |
| `KIR2DP1*0010201` | KIR2DP1 | 9 | 8 | no (has `/pseudo` exon feature) | **pseudogene** representative (KIR2DP1 is a KIR pseudogene locus) — full-length, becomes part of the emitted FASTA |
| `KIR2DP1*00101` | KIR2DP1 | 3 | 0 | **yes** (`FT ... partial` on the record) | genuinely partial pseudogene allele, included specifically to exercise `ParseDatFile.pl`'s partial-allele-rescue logic (`%partialAlleles`, `keys %partialAlleles` iteration — the exact hash-order-dependent code path flagged as a determinism risk). Kept in the subset `.dat` even though (see below) it does **not** end up in the emitted FASTA. |

Only real, verified-present records were used (looked up by exact allele name
in the source `.dat`, not invented).

## Genome annotation GTF: `kir_genes.gtf` (806 bytes, 4 lines)

- **Source**: GENCODE human, **release v50** (GRCh38, Ensembl 116),
  `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v50.annotation.gtf.gz`
  (annotation header confirms `version 50 (Ensembl 116)`, dated 2026-04-08 in
  the GENCODE file header)
- **Date fetched**: 2026-07-01
- Downloaded the full `gtf.gz` (~124.5 MB compressed) to scratch space,
  decompressed, filtered to `chr19` + feature type `gene` only (3543 rows),
  then grepped down to just the 4 rows matching `gene_name "KIR2DL1"`,
  `"KIR2DL4"`, `"KIR3DL1"`, `"KIR2DP1"` exactly. Nothing beyond this
  806-byte, 4-line file was written into the repo; the full GTF was deleted
  from scratch space after extraction.
- Read `vendor/t1k/AddGeneCoord.pl` before subsetting: it only consumes GTF
  rows where column 3 (`feature`) is exactly `gene`, reads gene name from the
  `gene_name "..."` attribute (col 9), and reads chromosome from column 1
  (prefixing `chr` only if missing — GENCODE already has the `chr` prefix, so
  no rewriting occurred). It applies a single hardcoded name alias
  (`HFE:HLA-HFE`, irrelevant to KIR) and otherwise expects an exact
  `gene_name` string match against the gene parsed from the FASTA header
  (text before `*`). All 4 KIR gene names in the subset matched directly with
  no aliasing needed.
- Verified: after running `t1k-build.pl -g`, all 7 emitted alleles (across all
  4 genes) got real chr19 coordinates in both `kir_dna_coord.fa` and
  `kir_rna_coord.fa` — none were left at the `-1 -1` "not found" default.

## Golden generation (`fixtures/refbuild/golden/`)

Command:
```
perl vendor/t1k/t1k-build.pl -d fixtures/refbuild/kir_subset.dat -g fixtures/refbuild/kir_genes.gtf \
  -o fixtures/refbuild/golden --prefix kir
```

All 4 expected outputs were generated (coord files were **not** deferred —
the GTF resolved cleanly for every gene in the subset):

| File | Bytes | Alleles emitted |
|---|---|---|
| `kir_dna_seq.fa` | 27,418 | 7 |
| `kir_rna_seq.fa` | 9,381 | 7 |
| `kir_dna_coord.fa` | 27,055 | 7 |
| `kir_rna_coord.fa` | 9,087 | 7 |

Only 7 of the 8 `.dat` records are emitted, in both RNA and DNA mode: the one
genuinely partial allele, `KIR2DP1*00101` (3 exons, no introns), is excluded
in both modes. Traced through `ParseDatFile.pl`: it enters
`%partialAlleles`, then the "rescue partial alleles" block runs (rescue is
active by default — `t1k-build.pl` never passes `--partialInRnaMode`, so
`$includePartialDiffLen` stays `0`, and `0 >= 0` is true) but its effective
length (3 exons + 2×50bp UTR) falls well short of `$geneLengthMode{KIR2DP1}`
(computed from the two 9-exon/8-intron KIR2DP1 records), so it fails the
`$len >= $geneLengthMode{$gene} - $includePartialDiffLen` rescue test in both
`rna` and `dna` mode and is silently dropped. This is expected, correct
behavior of the oracle (not a bug) and is exactly the kind of partial-allele
edge case Task 1.2's Rust port must reproduce faithfully.

## Determinism check (defines the Phase-1 invariant)

Ran `t1k-build.pl` with the same inputs into 5 independent output directories
(`/tmp/rb_a` .. `/tmp/rb_e`) plus the committed `golden/` directory (6 runs
total), then did pairwise `diff -r`:

```
diff -r /tmp/rb_a /tmp/rb_b   -> exit 0, 0 lines of diff
diff -rq /tmp/rb_a /tmp/rb_c  -> exit 0
diff -rq /tmp/rb_a /tmp/rb_d  -> exit 0
diff -rq /tmp/rb_a /tmp/rb_e  -> exit 0
diff -rq fixtures/refbuild/golden /tmp/rb_a -> exit 0
```

**Result: all 6 runs are byte-identical across all 4 output files.**

**Chosen invariant for Phase 1 (with this fixture): raw per-file
byte-identity.** A Rust port of `t1k-build`/`ParseDatFile`/`AddGeneCoord`
should reproduce `kir_{dna,rna}_{seq,coord}.fa` byte-for-byte given
`kir_subset.dat` + `kir_genes.gtf` as input.

### Important caveat for Task 1.2 (read before assuming this proves hash-order robustness)

`ParseDatFile.pl` contains a Perl-hash-order-dependent code path in its
partial-allele-rescue logic: `foreach my $allele (keys %partialAlleles)` (used
both in `rna` mode, line ~481, and `dna` mode, line ~524) iterates a hash
whose key order is randomized per-process by Perl >= 5.18 (hash-seed
randomization, unless `PERL_HASH_SEED` is pinned). The order in which rescued
alleles are appended to `@alleleOrder` — and therefore their order in the
output FASTA — depends on this iteration order **whenever two or more
partial alleles are simultaneously eligible for rescue**.

This fixture's subset has exactly **one** partial allele
(`KIR2DP1*00101`), and it never actually gets rescued (see above) — so
`%partialAlleles` never has more than one live entry undergoing the
order-sensitive iteration in a way that could produce a visible ordering
difference. **The 6-way byte-identical result above is real and correctly
observed, but it does not yet exercise or prove robustness against the
genuine multi-way tie case.** It only demonstrates that byte-identity holds
for a fixture with zero rescue-order ambiguity.

Separately, `ParseDatFile.pl`'s `FindMode()` function *also* iterates
`keys %dist` (mode/majority-vote length selection used throughout: gene
length mode, exon-count mode, per-exon length mode, intron-length mode), but
its tie-break (`$dist{$k} == $max && $k ge $ret`) is a proper commutative/
associative running-max fold over strings for any tied maximal keys — it
converges to the same lexicographically-largest tied key regardless of Perl's
hash iteration order. **`FindMode`'s tie-break is therefore already
order-independent and does not need a special deterministic-tiebreak port in
Task 1.2** — a straightforward "pick max by string comparison among ties"
reimplementation in Rust (iterating any collection in any order) will match it
exactly.

`AddGeneCoord.pl` has no hash-order-dependent output: gene coordinates are
looked up once per FASTA record and the output order is simply the input
FASTA's record order (first read of `$ARGV[0]` builds `%geneCoord`, second
read walks the GTF and fills it in with "first GTF row wins" semantics keyed
by exact gene-name match, third read re-emits records in original FASTA
order). It is fully deterministic as written.

**Recommendation for Task 1.2**: if/when a fixture with >=2 tied partial
alleles for the same gene (or partial alleles across >=2 genes, to test
whether rescue order affects file order) is needed to pin down the *actual*
tie-break the Rust port must replicate, it should NOT attempt to replay Perl's
specific hash-seed/iteration order (that's neither stable nor portable). It
should instead pick a designed deterministic order — the natural, already-used
candidate is: append rescued alleles from `%partialAlleles` in the same
allele-name order they were first encountered while parsing the `.dat` file
(i.e. preserve original file/allele order for rescued alleles too, the same
way `@alleleOrder` already does for non-partial alleles), rather than Perl's
hash order. This task defers pinning that case since our natural, real-allele
subset didn't happen to produce a multi-way tie; fabricating one artificially
was avoided per the "do not fabricate data" instruction for this task.

## Deviations / concerns

- None blocking. The only note-worthy caveat is the rescue-order-tie coverage
  gap documented above — flagged for a follow-up fixture in Task 1.2, not a
  blocker for this task's deliverable.

---

# HLA subset and golden (Phase 1b, Task 1b.1)

This section pins a small IPD-IMGT/HLA subset and its byte-identical Perl-oracle
golden, extending the Phase-1 KIR pinning to HLA. Its distinguishing purpose is
to **exercise the `srand(17)`/`rand()` UTR-padding path** in `ParseDatFile.pl`
(`:575-602`), which on HLA fires only for the pseudogenes HLA-DRB2 and HLA-DRB7
(RNA mode only). See the Phase-0 spike report (spike #8) for the full analysis.

## Source data

- **Database**: IPD-IMGT/HLA release **3.64.0** `hla.dat` (EMBL-flatfile format,
  a series of records separated by `//\n`).
- **Local source used**: `/Volumes/scratch-00001/t1k-hla-run/hlaidx/hla.dat`
  (338 MB, 46,201 allele records; release confirmed by the
  `IPD-IMGT/HLA Release Version 3.64.0` CC banner in each record). Not
  re-downloaded — reused verbatim.

## Subset: `hla_subset.dat` (168,729 bytes, 13 records)

Extracted verbatim (byte-exact whole records, each terminated by its original
`//\n`) by splitting the source `.dat` on the `//\n` record separator and
selecting records whose `/allele="..."` FT qualifier matched a wanted name.
Only real, verified-present records were used (looked up by exact allele name,
never invented).

| Allele | Gene | Mode role | Why chosen |
|---|---|---|---|
| `HLA-A*01:01:01:01`   | HLA-A    | typed | common class-I representative (2 alleles/gene exercises mode/consensus logic) |
| `HLA-A*02:01:01:01`   | HLA-A    | typed | second common HLA-A allele |
| `HLA-B*07:02:01:01`   | HLA-B    | typed | common class-I representative |
| `HLA-B*08:01:01:01`   | HLA-B    | typed | second common HLA-B allele |
| `HLA-C*07:01:01:01`   | HLA-C    | typed | common class-I representative |
| `HLA-C*07:02:01:01`   | HLA-C    | typed | second common HLA-C allele |
| `HLA-DRB1*01:01:01:01`| HLA-DRB1 | typed | common class-II representative |
| `HLA-DRB1*03:01:01:01`| HLA-DRB1 | typed | second common HLA-DRB1 allele |
| `HLA-DQB1*02:01:01:01`| HLA-DQB1 | typed | common class-II representative |
| `HLA-DQB1*06:02:01:01`| HLA-DQB1 | typed | second common HLA-DQB1 allele |
| `HLA-DRB2*01:01`      | HLA-DRB2 | **RNG pseudogene** | **Required.** Pseudogene whose only allele supplies **no** full 50 bp 3' flank, so the gene-level "no allele has a complete 50 bp 3'UTR" condition holds and the `srand(17)` 3'UTR random-padding path fires (RNA mode). Exercises Task 1b.2's drand48 port. |
| `HLA-DRB7*01:01:01`   | HLA-DRB7 | **RNG pseudogene** | **Required.** Same as DRB2 — its 3'UTR is filled with 50 random bytes from the seeded PRNG in RNA mode. |
| `HLA-DRB7*01:01:02`   | HLA-DRB7 | partial | second DRB7 allele, marked `partial` in the `.dat`; kept so DRB7 has >1 allele and to exercise `ParseDatFile.pl`'s `%partialAlleles` rescue path (analogous to KIR2DP1*00101 above). Like that KIR allele, it is **not** emitted to the FASTA (dropped by the length-mode rescue test). |

Note: the typed genes (A/B/C/DRB1/DQB1) never reach the RNG path — confirmed
below. Only the two pseudogenes do, and only in RNA mode.

## Genome annotation GTF: `hla_genes.gtf` (969 bytes, 5 lines)

- **Source**: GENCODE human, release **v50** (GRCh38), from
  `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v50.annotation.gtf.gz`
- **Date fetched**: 2026-07-01 (same GENCODE v50 release used for the KIR GTF above).
- Filtered to `chr6` + feature type `gene` + `gene_name` in
  {HLA-A, HLA-B, HLA-C, HLA-DRB1, HLA-DQB1, HLA-DRB2, HLA-DRB7}. The download
  was still in progress but had already fully covered chr6 and reached chr17,
  so the chr6 MHC region is complete; the large `gtf.gz` lives only in scratch
  and is not committed.
- **DRB2/DRB7 are absent from GENCODE v50** (they are pseudogenes GENCODE does
  not annotate), so only the 5 typed genes appear in the GTF. This is expected:
  per `AddGeneCoord.pl`, any FASTA gene not present in the GTF keeps the default
  not-found sentinel `chr19 -1 -1 +` (the script's `$defaultChr` is hardcoded to
  `chr19`, a KIR-era default it applies verbatim to HLA too). Reproducing this
  sentinel behavior for DRB2/DRB7 is part of the pinned oracle output.
- Read `vendor/t1k/AddGeneCoord.pl` before subsetting: it consumes only rows
  whose column 3 is exactly `gene`, reads the gene name from `gene_name "..."`
  in col 9 and the chromosome from col 1 (prefixing `chr` if absent; GENCODE
  already has it). It applies one hardcoded alias (`HFE:HLA-HFE`, irrelevant
  here) and otherwise requires an exact `gene_name` match against the gene
  parsed from the FASTA header (text before `*`). All 5 typed HLA gene names
  matched directly with no aliasing.

## Golden generation (`fixtures/refbuild/golden/`)

Command:
```
perl vendor/t1k/t1k-build.pl -d fixtures/refbuild/hla_subset.dat \
  -g fixtures/refbuild/hla_genes.gtf -o fixtures/refbuild/golden --prefix hla
```

| File | Bytes | Records | Notes |
|---|---|---|---|
| `hla_dna_seq.fa`   | 27,896 | 10 | DRB2/DRB7 pseudogenes not emitted in DNA mode |
| `hla_rna_seq.fa`   | 13,191 | 12 | includes DRB2*01:01 and DRB7*01:01:01 with RNG-padded 3'UTR |
| `hla_dna_coord.fa` | 27,502 | 10 | all 5 typed genes got real chr6 coords |
| `hla_rna_coord.fa` | 12,787 | 12 | 5 typed genes → chr6 coords; DRB2/DRB7 → `chr19 -1 -1 +` sentinel |

DRB7*01:01:02 (the `partial` allele) is dropped in both modes by the
length-mode rescue test — the same mechanism that drops KIR2DP1*00101 above.

## RNG-path confirmation (critical — the reason for this subset)

Confirmed empirically on **this subset** using a COPY of `ParseDatFile.pl` in
scratch (the vendored original was never edited; verified by `diff`):

1. **Instrumented copy** (`warn` added at both RNG blocks) on the subset:
   ```
   RNA mode warnings:  RNG5UTR gene=HLA-DRB2 bestlen=512
                       RNG3UTR gene=HLA-DRB2 bestlen=0   <-- 50 truly-random bytes
                       RNG5UTR gene=HLA-DRB7 bestlen=796
                       RNG3UTR gene=HLA-DRB7 bestlen=0   <-- 50 truly-random bytes
   DNA mode warnings:  (none)
   ```
   The instrumented output was byte-identical to the un-instrumented seed-17
   output (warn-only, no behavior change). `bestlen=0` on the two 3'UTR entries
   means no real flank exists, so all 50 padded bytes are pure PRNG output.

2. **Seed-perturbation test** — ran the subset through a `srand(999)` scratch
   copy and compared records order-independently against seed-17:
   ```
   RNA mode: exactly 2 records differ -> HLA-DRB2*01:01, HLA-DRB7*01:01:01
             (only their trailing 50 bp 3'UTR bytes)
   DNA mode: 0 records differ
   ```
   This proves the seeded RNG stream reaches output for exactly these two
   records, RNA mode only — the path Task 1b.2's drand48 port must reproduce.

3. **Seed-17 3'UTR random tails for THIS subset** (present verbatim in the
   committed `golden/hla_rna_seq.fa`; these ARE the drand48 output to match):
   ```
   HLA-DRB2*01:01     3'UTR tail: CTGTCGATGCTTCCACGGAAGATACGTGCCAGACAGTTCCGATAAATTTA
   HLA-DRB7*01:01:01  3'UTR tail: TCTGTAACAGACACGCTAGTCGGAAGCCGTGAACTCACTGTCCTGCGTAG
   ```
   These differ from the full-database spike-report tails **by design**: the
   `rand()` consumption order follows `@alleleOrder` (file-parse order of all
   parsed alleles), which differs between the 46k-record full `.dat` and this
   13-record subset. Each input has its own deterministic seed-17 stream; the
   golden pins the stream for THIS subset.

## Determinism check and chosen invariant

Ran the full golden generation (all 4 files, incl. `AddGeneCoord`) into two
independent temp dirs plus 5 additional seq-only runs. All comparisons:

```
hla_dna_seq  : RAW IDENTICAL   hla_rna_seq  : RAW IDENTICAL
hla_dna_coord: RAW IDENTICAL   hla_rna_coord: RAW IDENTICAL
(golden vs each temp run: RAW IDENTICAL for all 4 files)
```

**Chosen invariant for this HLA fixture: raw per-file byte-identity** (same as
Phase-1 KIR). A Rust port must reproduce `hla_{dna,rna}_{seq,coord}.fa`
byte-for-byte given `hla_subset.dat` + `hla_genes.gtf`.

**Caveat (from spike #8), important for Task 1b.2's oracle harness on the FULL
database:** `ParseDatFile.pl` output record ORDER is non-deterministic on large
inputs due to Perl per-process hash-key randomization in the partial-allele
rescue (`keys %partialAlleles`). On the full 46k-record `hla.dat`, raw file
bytes differ run-to-run while content sorted by header is identical. This small
13-record subset has only ONE partial allele (DRB7*01:01:02, never rescued), so
no rescue-order ambiguity arises and raw byte-identity holds here — but a
general HLA oracle should compare **records sorted by header**, not raw bytes.
The RNG-padded 3'UTR bytes themselves are deterministic in CONTENT (drand48),
independent of record order.

## Deviations / concerns

- None blocking. As with KIR, this subset's raw byte-identity does not exercise
  the multi-way partial-allele rescue-order tie (only one partial allele, never
  rescued) — the general full-database invariant is sort-by-header per spike #8.
- The GENCODE v50 `gtf.gz` download was still streaming when the chr6 rows were
  extracted; chr6 (through chr17) was fully present, so the 5 HLA gene rows are
  complete and correct. The partial archive was not committed (scratch only).
