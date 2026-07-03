//! Genotyper foundation, ported from T1K's `Genotyper`
//! (`Genotyper.hpp`).
//!
//! This is the Phase-5a slice ONLY: the deterministic pieces the EM
//! (`EMupdate`/`Quantify`, Phase 5b) and allele selection (`SelectAllele`/
//! output, Phase 5c) build on top of. Specifically:
//!
//! - [`parse_allele_name`] -- ported from `Genotyper::ParseAlleleName`
//!   (`Genotyper.hpp:63-131`): the KIR (no-delimiter)/HLA (delimiter)
//!   allele-name parser that yields `(gene, majorAllele)`.
//! - [`AlleleInfo`] + [`Genotyper::init_allele_info`] -- ported from the
//!   `_alleleInfo` struct (`Genotyper.hpp:16-31`) and
//!   `Genotyper::InitAlleleInfo` (`Genotyper.hpp:559-682`): the allele/gene
//!   data model and its initialization (gene/major-allele grouping, the
//!   gene-level k-mer-similarity matrix, and effective-length outlier
//!   correction).
//! - [`Genotyper::read_assignment_weight`] -- ported from
//!   `Genotyper::ReadAssignmentWeight` (`Genotyper.hpp:205-230`): the
//!   deterministic per-fragment-overlap weight (no `pow`/`log`, hence
//!   bit-identical to the C++ oracle).
//! - [`FragmentOverlap`]/[`ReadAssignment`]/[`ReadGroupInfo`] +
//!   [`Genotyper::init_read_assignments`]/[`Genotyper::set_read_assignments`]
//!   -- ported from the `_fragmentOverlap` (`SeqSet.hpp:146-173`) and
//!   `_readAssignment`/`_readGroupInfo` (`Genotyper.hpp:33-56`) structs and
//!   `Genotyper::InitReadAssignments`/`SetReadAssignments`
//!   (`Genotyper.hpp:759-832`): the per-read fragment-overlap storage and
//!   weight assignment that Phase 5b's EM consumes.
//!
//! This module also carries the Phase-5b slice: the quantification core built
//! on top of 5a's foundation.
//!
//! - [`Genotyper::coalesce_read_assignments`] -- ported from
//!   `Genotyper::CoalesceReadAssignments` (`Genotyper.hpp:841-908`):
//!   collapses per-read allele assignments with an identical (sorted)
//!   allele-index/qual fingerprint into shared `readAssignments` groups.
//!   Deterministic (a fingerprint hash + full tie-break comparison, no
//!   floating point beyond `==` on already-computed weights) -> byte-identical.
//! - [`Genotyper::finalize_read_assignments`] -- ported from
//!   `Genotyper::FinalizeReadAssignments` (`Genotyper.hpp:912-939`): builds
//!   `readsInAllele` from the coalesced `readAssignments`, then calls
//!   [`Genotyper::build_allele_equivalent_class`] and sets each allele's
//!   `missing_coverage` from a caller-supplied slice (mirrors
//!   `refSet.GetSeqMissingBaseCoverage`, out of this port's scope -- same
//!   "caller supplies the reference-derived input" pattern as
//!   [`Genotyper::init_allele_info`]'s `seq_effective_len`/`seq_weight`).
//! - [`Genotyper::build_allele_equivalent_class`] -- ported from
//!   `Genotyper::BuildAlleleEquivalentClass` (`Genotyper.hpp:1072-1139`):
//!   groups alleles that were hit by the identical (multi)set of coalesced
//!   reads into equivalence classes, via a fingerprint hash sorted
//!   descending (full tie-break, so `std::sort`'s non-stability never
//!   matters) then a linear adjacent-fingerprint scan with a real
//!   read-list-equality check. Deterministic -> byte-identical; this
//!   grouping (including EC index order) is load-bearing for the EM and for
//!   5c's allele selection.
//! - [`Genotyper::remove_low_mapq_allele_in_equivalent_class`] -- ported from
//!   `Genotyper::RemoveLowMAPQAlleleInEquivalentClass`
//!   (`Genotyper.hpp:1330-1368`): within each EC, keeps only the allele(s)
//!   with the maximum total per-read qual sum. Deterministic (`==` on
//!   already-computed sums) -> byte-identical.
//! - [`Genotyper::quantify_allele_equivalent_class`] -- ported from
//!   `Genotyper::QuantifyAlleleEquivalentClass` (`Genotyper.hpp:1142-1328`):
//!   builds `readGroupToAlleleEc`/[`EcInfo`], then runs the SQUAREM-
//!   accelerated EM (two [`Genotyper::em_update`] calls, one
//!   [`Genotyper::squarem_alpha`], the quadratic extrapolation, a
//!   stabilizing `EMupdate`, convergence check, periodic low-abundance
//!   masking via [`Genotyper::set_allele_abundance`]) to completion, then a
//!   final `SetAlleleAbundance`. See the FLOATS note on
//!   [`Genotyper::em_update`] for the accumulation-order discipline this
//!   preserves, and why that yields abundances matching the C++ oracle
//!   within a tight relative tolerance rather than exact `f64` bits.
//! - [`Genotyper::em_update`] / [`Genotyper::squarem_alpha`] /
//!   [`Genotyper::set_allele_abundance`] -- ported from `Genotyper::EMupdate`
//!   (`Genotyper.hpp:372-421`), `Genotyper::SQUAREMalpha`
//!   (`Genotyper.hpp:424-437`), and `Genotyper::SetAlleleAbundance`
//!   (`Genotyper.hpp:957-1014`, the `ecReadCount != NULL` branch only --
//!   5c's `InitAlleleAbundance` calls the `NULL` branch, out of this port's
//!   scope): the float kernels. `EMupdate`/`SQUAREMalpha` are pure
//!   `+`/`*`/`/`/`sqrt` (no `pow`/`log`/`exp`) over the same
//!   iteration/accumulation order as the C++ (no parallelism, no
//!   reordering) -- see their doc comments for why this port targets, and
//!   in practice achieves, per-EC/per-allele abundances that match the C++
//!   oracle within a tight relative tolerance (`1e-6`; observed divergence
//!   is far smaller, ~1e-8 on called alleles). This is NOT exact-`f64`-bits
//!   reproducibility: the vendored C++ oracle is compiled at `-O3`, which
//!   enables FMA contraction (fusing `a * b + c` into one, more-precisely-
//!   rounded operation), while Rust's `+`/`*` always round twice. That
//!   per-operation difference compounds over the SQUAREM loop's iterations,
//!   so exact bit-for-bit agreement is an inherently compiler/platform-
//!   dependent outcome, not a property of the algorithm -- it happens to
//!   hold in degenerate cases (e.g. a single dominant equivalence class)
//!   where FMA contraction has no opportunity to change the rounding
//!   outcome, but is not a general invariant. The deterministic structure
//!   this module also computes (read groups, EC membership/order,
//!   `equivalent_class`, `missing_coverage`) has no floating point in its
//!   computation and remains byte-identical to the C++ oracle. See
//!   `crates/unum-core/tests/golden_genotyper_em.rs` for the differential
//!   test (structure exact, abundances within tolerance) and its
//!   `assert_close_within_tolerance` doc comment for the full justification;
//!   the genotype CALL itself (Phase 5c, end-to-end) is the ultimate
//!   correctness gate, not bitwise abundance parity.
//!
//! This module also carries the Phase-5c slice: the capstone that turns
//! candidate reads into a called, quality-scored genotype and writes the
//! output TSVs. See `crates/unum/src/stages/genotype.rs` for the CLI/read-
//! processing driver (reference loading/dedup, [`extend_overlap`]/
//! [`assign_read`]/[`read_assignment_to_fragment_assignment`] wiring, the
//! main read loop) that ties this module's pieces together end-to-end.
//!
//! - [`alnorm`] -- ported from `Genotyper::alnorm` (`Genotyper.hpp:252-370`):
//!   the AS66 standard-normal-CDF polynomial approximation
//!   [`Genotyper::select_alleles_for_genes`]'s quality score is built from.
//! - [`extend_overlap`] -- ported from `SeqSet::ExtendOverlap`
//!   (`SeqSet.hpp:1994-2100`): extends a k-mer-chained [`overlap::Overlap`]
//!   to the read's/allele's full overhangs via `AlignAlgo::GlobalAlignment`.
//!   Only the `ignoreNonExonDiff == false` path is ported (see its doc
//!   comment) -- correct for every caller in this port's scope, since
//!   `Genotyper.cpp`'s CLI only sets `ignoreNonExonDiff` true via
//!   `--relaxIntronAlign`, which no invocation this port targets passes.
//! - [`assign_read`] -- ported from `SeqSet::AssignRead` (`SeqSet.hpp:2119-
//!   2303`), restricted to `weight == -1` (no base-coverage marking; see its
//!   doc comment for why `posWeight`/`GetSeqMissingBaseCoverage` are out of
//!   this port's scope) and `ignoreNonExonDiff == false` (same as
//!   [`extend_overlap`]): turns one (deduplicated) read sequence into its
//!   sorted, [`extend_overlap`]-refined list of [`overlap::Overlap`]s.
//! - [`read_assignment_to_fragment_assignment`] -- ported from
//!   `SeqSet::ReadAssignmentToFragmentAssignment` (`SeqSet.hpp:2310-2655`):
//!   combines a read pair's two mate-end overlap lists (as produced by
//!   [`assign_read`]) into per-allele [`FragmentOverlap`]s, the direct input
//!   to [`Genotyper::set_read_assignments`].
//! - [`assign_reads_parallel`] -- the P2 parallelism-campaign entry point:
//!   runs [`RefKmerFilter::get_overlaps_from_read`] then [`assign_read`] over
//!   every distinct read sequence, parallelizing only the read-only
//!   `get_overlaps_from_read` phase across a scoped `rayon` thread pool (each
//!   worker with its own [`Scratch`]) while keeping `assign_read`'s
//!   `allele_refs`-mutating phase strictly sequential in input order -- see
//!   its doc comment for why this makes `-t N` output byte-identical to
//!   `-t 1` at any `N`.
//! - [`Genotyper::remove_low_likelihood_allele_in_equivalent_class`] --
//!   ported from `Genotyper::RemoveLowLikelihoodAlleleInEquivalentClass`
//!   (`Genotyper.hpp:1371-1460`): a coverage-likelihood filter over each
//!   equivalence class's member alleles, run once (`Genotyper.cpp:647`)
//!   between quantification and allele selection.
//! - [`Genotyper::select_alleles_for_genes`] -- ported from
//!   `Genotyper::SelectAllelesForGenes` (`Genotyper.hpp:1462-2090`): the
//!   genotyping decision itself -- ranks equivalence classes by abundance,
//!   assigns alleles to genes/ranks (0/1/2+), runs the iterative
//!   best-haplotype-pair search, then scores each rank's `genotypeQuality`
//!   via [`alnorm`]. This is the capstone this whole module (5a/5b/5c)
//!   builds up to.
//! - [`Genotyper::get_gene_allele_types`]/
//!   [`Genotyper::is_reads_in_allele_idx_optimal`]/
//!   [`Genotyper::get_average_read_assignment_cnt`] -- small ported helpers
//!   [`Genotyper::select_alleles_for_genes`] and the output writers below
//!   depend on.
//! - [`Genotyper::get_allele_description`]/[`Genotyper::output_representative_alleles`]
//!   -- ported from `Genotyper::GetAlleleDescription` (`Genotyper.hpp:2103-
//!   2178`, main's `_genotype.tsv` row formatting) and
//!   `Genotyper::OutputRepresentativeAlleles` (`Genotyper.hpp:2180-2229`,
//!   `_allele.tsv`).

use std::collections::HashMap;

use crate::overlap;
use crate::ref_kmer_filter::{RefKmerFilter, Scratch};

use crate::kmer_count::KmerCount;
use rayon::prelude::*;

/// KIR gene type constant, ported from `GENETYPE_KIR` (`Genotyper.hpp:13`).
pub const GENE_TYPE_KIR: i32 = 0;
/// HLA gene type constant, ported from `GENETYPE_HLA` (`Genotyper.hpp:14`).
pub const GENE_TYPE_HLA: i32 = 1;

/// Ported from `Genotyper::ParseAlleleName` (`Genotyper.hpp:63-131`).
///
/// Splits an allele name (e.g. `"KIR2DL1*0010101"` or `"A*01:01:01:01"`)
/// into `(gene, majorAllele)`.
///
/// - `allele_digit_units` mirrors the `alleleDigitUnits` field
///   (`Genotyper.hpp:481`, set via `SetAlleleNameStructure`): `-1` (the
///   `Genotyper` constructor default, `Genotyper.hpp:506`) means "derive
///   `fieldsLength` from `fields_type` and whether `allele` contains a `:`
///   delimiter"; any other value is used as `fieldsLength` directly.
/// - `allele_delimiter` mirrors the `alleleDelimiter` field
///   (`Genotyper.hpp:482`, constructor default `'\0'`,
///   `Genotyper.hpp:507`): a non-NUL value forces `parseType = 2`
///   (delimiter-based/HLA-style parsing) using that exact byte as the
///   delimiter, overriding whatever delimiter auto-detection found.
/// - `fields_type` mirrors the `fieldsType` parameter (default `0` at every
///   call site in stock `Genotyper.hpp`/`.cpp` except `IsAlleleSameInExon`,
///   which passes `1`, `Genotyper.hpp:139-140`): `0` always yields
///   `fieldsLength = 3` (once `allele_digit_units == -1`); `>= 1` yields `5`
///   for no-delimiter (KIR) names and `3` for delimiter (HLA) names.
///
/// # Buffer semantics vs. the C++
///
/// The C++ writes into caller-owned `gene`/`majorAllele` `char*` buffers via
/// `strcpy` (so `gene`/`majorAllele` START as full copies of `allele`, then
/// get NUL-truncated in place). This port has no such aliasing: it directly
/// slices `allele` at the same byte offsets the C++ loops compute, which
/// yields byte-identical output without needing an initial full-string copy.
///
/// # Panics
///
/// Panics if `allele` is not valid UTF-8 sliceable at the computed byte
/// offsets -- not expected for real IMGT/HLA or KIR allele names, which are
/// pure ASCII.
#[must_use]
pub fn parse_allele_name(
    allele: &str,
    fields_type: i32,
    allele_digit_units: i32,
    allele_delimiter: u8,
) -> (String, String) {
    let bytes = allele.as_bytes();
    let n = bytes.len();

    let mut parse_type = 1; // 1 = no delimiter (KIR); 2 = delimiter (HLA).
    let mut fields_length = allele_digit_units;
    let mut delimiter: u8 = b'\0';

    if fields_length == -1 {
        fields_length = 3;
        for &b in bytes {
            if b == b':' {
                delimiter = b':';
                parse_type = 2;
            }
        }

        if fields_type == 0 {
            fields_length = 3;
        } else if fields_type >= 1 {
            if parse_type == 1 {
                fields_length = 5;
            } else if parse_type == 2 {
                fields_length = 3;
            }
        }
    }
    if allele_delimiter != b'\0' {
        delimiter = allele_delimiter;
        parse_type = 2;
    }

    if parse_type == 1 {
        // Ported from Genotyper.hpp:99-110.
        let mut i = 0usize;
        while i < n && bytes[i] != b'*' {
            i += 1;
        }
        let gene_end = i;

        // `for (j = 0 ; j <= fieldsLength && allele[i + j] ; ++j) ;` --
        // advances j while j <= fieldsLength AND allele[i+j] is a live
        // (non-NUL, i.e. in-bounds) byte. `j` only ever increments from 0,
        // so it is always non-negative -- the `as usize` cast never loses a
        // sign bit.
        let mut j: i32 = 0;
        #[allow(clippy::cast_sign_loss)]
        while j <= fields_length && i + (j as usize) < n {
            j += 1;
        }
        #[allow(clippy::cast_sign_loss)]
        let major_end = (i + (j as usize)).min(n);

        (
            String::from_utf8_lossy(&bytes[..gene_end]).into_owned(),
            String::from_utf8_lossy(&bytes[..major_end]).into_owned(),
        )
    } else {
        // parse_type == 2. Ported from Genotyper.hpp:111-130.
        let mut i = 0usize;
        while i < n && bytes[i] != b'*' {
            i += 1;
        }
        let gene_end = i;

        let mut k = 0;
        let mut j = i;
        while j < n {
            if bytes[j] == delimiter {
                k += 1;
                if k >= fields_length {
                    break;
                }
            }
            j += 1;
        }
        let major_end = j;

        (
            String::from_utf8_lossy(&bytes[..gene_end]).into_owned(),
            String::from_utf8_lossy(&bytes[..major_end]).into_owned(),
        )
    }
}

/// Ported from `Genotyper::alnorm` (`Genotyper.hpp:252-370`): Algorithm AS66
/// (Hill 1973), the polynomial approximation to the standard normal
/// cumulative distribution used by [`Genotyper::select_alleles_for_genes`]'s
/// `genotypeQuality` scoring (`Genotyper.hpp:2068`,
/// `-log(alnorm(...))/log(10)`).
///
/// `upper == true` integrates from `x` to `+Infinity`; `upper == false`
/// integrates from `-Infinity` to `x`. All coefficients/branches/constants
/// are copied verbatim from the vendored C++ (down to variable names in this
/// doc comment) -- this is a fixed numerical-methods polynomial, not
/// something to "clean up" or re-derive.
///
/// # FLOATS
///
/// Pure `+`/`*`/`/`/`exp`/`abs`/comparisons, no loops or data-dependent
/// accumulation order -- a straight-line polynomial evaluation. `exp` is
/// libm's `f64::exp`; like [`Genotyper::squarem_alpha`]'s `sqrt` note, this
/// single transcendental call is not itself a source of divergence, but the
/// C++ oracle's `-O3` build may still contract the surrounding
/// multiply-then-divide chains via FMA in ways Rust's non-fused `+`/`*` does
/// not -- so this targets (and the `diff_genotype_e2e` end-to-end test
/// verifies at the genotype-CALL tier) close agreement with the oracle, not
/// guaranteed bit-identical `f64` output.
#[must_use]
#[allow(clippy::many_single_char_names)] // matches the C++'s own single-letter variable names verbatim.
pub fn alnorm(x: f64, upper: bool) -> f64 {
    let a1 = 5.758_854_804_58;
    let a2 = 2.624_331_216_79;
    let a3 = 5.928_857_244_38;
    let b1 = -29.821_355_780_7;
    let b2 = 48.695_993_069_2;
    let c1 = -0.000_000_038_052;
    let c2 = 0.000_398_064_794;
    let c3 = -0.151_679_116_635;
    let c4 = 4.838_591_280_8;
    let c5 = 0.742_380_924_027;
    let c6 = 3.990_194_170_11;
    let con = 1.28;
    let d1 = 1.000_006_153_02;
    let d2 = 1.986_153_813_64;
    let d3 = 5.293_303_249_26;
    let d4 = -15.150_897_245_1;
    let d5 = 30.789_933_034;
    let ltone = 7.0;
    let p = 0.398_942_280_444;
    let q = 0.399_903_485_04;
    let r = 0.398_942_280_385;
    let utzero = 18.66;

    let mut up = upper;
    let mut z = x;

    if z < 0.0 {
        up = !up;
        z = -z;
    }

    if ltone < z && (!up || utzero < z) {
        return if up { 0.0 } else { 1.0 };
    }

    let y = 0.5 * z * z;

    let mut value = if z <= con {
        0.5 - z * (p - q * y / (y + a1 + b1 / (y + a2 + b2 / (y + a3))))
    } else {
        r * (-y).exp()
            / (z + c1
                + d1 / (z + c2 + d2 / (z + c3 + d3 / (z + c4 + d4 / (z + c5 + d5 / (z + c6))))))
    };

    if !up {
        value = 1.0 - value;
    }

    value
}

/// Ported from `_alleleInfo` (`Genotyper.hpp:16-31`): per-allele bookkeeping
/// populated by [`Genotyper::init_allele_info`] and consumed by later stages
/// (EM/selection, Phase 5b/5c).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlleleInfo {
    /// `_alleleInfo::majorAlleleIdx` -- index into
    /// [`Genotyper::major_allele_idx_to_name`].
    pub major_allele_idx: i32,
    /// `_alleleInfo::geneIdx` -- index into [`Genotyper::gene_idx_to_name`].
    pub gene_idx: i32,
    /// `_alleleInfo::alleleRank` -- `-1` not selected, `0` first allele, `1`
    /// second allele. Not touched by [`Genotyper::init_allele_info`] beyond
    /// its `-1` initialization (Phase 5c sets this).
    pub allele_rank: i32,
    /// `_alleleInfo::genotypeQuality` -- assignment quality. Not touched by
    /// [`Genotyper::init_allele_info`] beyond its `-1` initialization (Phase
    /// 5c sets this).
    pub genotype_quality: i32,
    /// `_alleleInfo::abundance`.  Not touched by [`Genotyper::init_allele_info`]
    /// beyond its `0` initialization (Phase 5b's EM sets this).
    pub abundance: f64,
    /// `_alleleInfo::equivalentClass` -- the class id for alleles sharing the
    /// same set of read alignments. Not populated by
    /// [`Genotyper::init_allele_info`] (Phase 5b's
    /// `BuildAlleleEquivalentClass` sets this); defaults to `0`, mirroring
    /// `_alleleInfo`'s C++ default-initialization (the field is never
    /// explicitly set in `InitAlleleInfo`).
    pub equivalent_class: i32,
    /// `_alleleInfo::ecAbundance`. Same "not populated here" note as
    /// `equivalent_class`; defaults to `0.0`.
    pub ec_abundance: f64,
    /// `_alleleInfo::missingCoverage`. Same "not populated here" note as
    /// `equivalent_class` (`FinalizeReadAssignments` sets this, Phase 5b);
    /// defaults to `0`.
    pub missing_coverage: i32,
    /// `_alleleInfo::whitelist` -- `true` unless
    /// `Genotyper::SetAlleleWhitelist` narrows it (not ported here; every
    /// allele stays whitelisted, matching `InitAlleleInfo`'s own
    /// `alleleInfo[i].whitelist = true` unconditional initialization,
    /// `Genotyper.hpp:593`).
    pub whitelist: bool,
}

impl Default for AlleleInfo {
    fn default() -> Self {
        Self {
            major_allele_idx: 0,
            gene_idx: 0,
            allele_rank: -1,
            genotype_quality: -1,
            abundance: 0.0,
            equivalent_class: 0,
            ec_abundance: 0.0,
            missing_coverage: 0,
            whitelist: true,
        }
    }
}

/// Ported from `_fragmentOverlap` (`SeqSet.hpp:146-173`): a fragment-level
/// (i.e. read-pair-level) overlap against one reference allele. Produced by
/// `SeqSet::ReadAssignmentToFragmentAssignment` (out of this port's 5a
/// scope, see module docs) and consumed by
/// [`Genotyper::read_assignment_weight`]/[`Genotyper::set_read_assignments`].
///
/// Fields are a subset of the full C++ struct: `overlap1`/`overlap2` (the
/// per-mate `_overlap`s) are omitted because neither
/// `ReadAssignmentWeight` nor `SetReadAssignments` reads them (the
/// commented-out `similarity` fallback at `Genotyper.hpp:210-211` that WOULD
/// read `o.overlap2.similarity` is dead code in stock T1K -- see
/// [`Genotyper::read_assignment_weight`]'s doc comment).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FragmentOverlap {
    /// `_fragmentOverlap::seqIdx` -- which allele this overlap is against.
    pub seq_idx: i32,
    /// `_fragmentOverlap::seqStart`.
    pub seq_start: i32,
    /// `_fragmentOverlap::seqEnd`.
    pub seq_end: i32,
    /// `_fragmentOverlap::matchCnt`.
    pub match_cnt: i32,
    /// `_fragmentOverlap::relaxedMatchCnt`.
    pub relaxed_match_cnt: i32,
    /// `_fragmentOverlap::similarity` -- consumed directly by
    /// [`Genotyper::read_assignment_weight`].
    pub similarity: f64,
    /// `_fragmentOverlap::hasMatePair`.
    pub has_mate_pair: bool,
    /// `_fragmentOverlap::o1FromR2`.
    pub o1_from_r2: bool,
    /// `_fragmentOverlap::qual` -- copied into
    /// [`ReadAssignment::qual`] by [`Genotyper::set_read_assignments`].
    pub qual: f64,
    /// `_fragmentOverlap::hasN` -- consumed by
    /// [`Genotyper::read_assignment_weight`]'s `ret /= 10.0` branch.
    pub has_n: bool,
}

/// Ported from `_readAssignment` (`Genotyper.hpp:44-56`): one allele
/// assignment for a single read, with its
/// [`Genotyper::read_assignment_weight`]-computed weight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReadAssignment {
    /// `_readAssignment::alleleIdx`.
    pub allele_idx: i32,
    /// `_readAssignment::start`.
    pub start: i32,
    /// `_readAssignment::end`.
    pub end: i32,
    /// `_readAssignment::weight` -- `f32` in the C++ struct (`Genotyper.hpp:48`,
    /// unlike the `f64` `ReadAssignmentWeight` return type it's assigned
    /// from at `Genotyper.hpp:827`); matched here for byte-identical
    /// storage (the narrowing f64->f32 truncation happens at the same
    /// assignment point as the C++).
    pub weight: f32,
    /// `_readAssignment::qual`.
    pub qual: f32,
    /// `_readAssignment::adjustWeight` -- the tie-break weight
    /// (`adjustFactor * weight`, see [`Genotyper::set_read_assignments`]).
    pub adjust_weight: f32,
}

/// Ported from `_readGroupInfo` (`Genotyper.hpp:33-36`): trivial per-read-
/// group bookkeeping (how many reads collapsed into this group). See the
/// module docs' "deliberately not ported" section -- this struct is defined
/// here for Phase 5b's `EMupdate` to reuse, but nothing in 5a populates it
/// (population happens inside `Quantify`, Phase 5b's territory).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReadGroupInfo {
    /// `_readGroupInfo::count`.
    pub count: f64,
}

/// Ported from `_ecInfo` (`Genotyper.hpp:38-42`): combined per-equivalent-
/// class bookkeeping the EM reads (never writes) during
/// [`Genotyper::quantify_allele_equivalent_class`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EcInfo {
    /// `_ecInfo::length` -- the MIN effective length across the EC's member
    /// alleles (`Genotyper.hpp:1203-1216`).
    pub length: i32,
    /// `_ecInfo::missingCoverage` -- the MIN missing-coverage across the
    /// EC's member alleles.
    pub missing_coverage: i32,
}

/// Ported from `_pairID` (`defs.h:28-32`): a generic `(int, double)` pair.
/// Used here as `readGroupToAlleleEc`'s element type -- `a` is an
/// equivalent-class index, `b` is that class's qual for the read group.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairId {
    /// `_pairID::a`.
    pub a: i32,
    /// `_pairID::b`.
    pub b: f64,
}

/// An [`overlap::Overlap`] refined by [`extend_overlap`]/[`assign_read`]:
/// carries the extra fields `_overlap` has that `overlap::Overlap` (Phase 4's
/// `GetOverlapsFromHits`-only port) does not need -- `relaxed_match_cnt`
/// (`_overlap::relaxedMatchCnt`) and `left_clip`/`right_clip`
/// (`_overlap::leftClip`/`rightClip`). Defined here rather than added to
/// [`overlap::Overlap`] itself, since that type is Phase 4's and this port's
/// scope is additive-only on top of it (see this module's doc comment).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExtendedOverlap {
    /// `_overlap::seqIdx`.
    pub seq_idx: u32,
    /// `_overlap::readStart`.
    pub read_start: i32,
    /// `_overlap::readEnd`.
    pub read_end: i32,
    /// `_overlap::seqStart`.
    pub seq_start: i32,
    /// `_overlap::seqEnd`.
    pub seq_end: i32,
    /// `_overlap::strand`.
    pub strand: i8,
    /// `_overlap::matchCnt` -- count TWICE, as with [`overlap::Overlap::match_cnt`].
    pub match_cnt: i32,
    /// `_overlap::similarity`.
    pub similarity: f64,
    /// `_overlap::relaxedMatchCnt`.
    pub relaxed_match_cnt: i32,
    /// `_overlap::leftClip`.
    pub left_clip: i32,
    /// `_overlap::rightClip`.
    pub right_clip: i32,
}

/// Ported from `_overlap::operator<` (`SeqSet.hpp:103-127`), specialized to
/// [`ExtendedOverlap`] -- identical field-by-field comparison to
/// [`overlap::overlap_less_than`] (both types share the same first eight
/// compared fields in the same order), duplicated here because
/// [`ExtendedOverlap`] is a distinct Rust type from [`overlap::Overlap`]. See
/// [`overlap::overlap_less_than`]'s doc comment for the full field-order
/// rationale (a strict weak ordering, not a total `Ord`).
#[must_use]
pub fn extended_overlap_less_than(a: &ExtendedOverlap, b: &ExtendedOverlap) -> bool {
    if a.match_cnt != b.match_cnt {
        return a.match_cnt > b.match_cnt;
    }
    #[allow(clippy::float_cmp)]
    let similarity_differs = a.similarity != b.similarity;
    if similarity_differs {
        return a.similarity > b.similarity;
    }
    let a_span = a.read_end - a.read_start;
    let b_span = b.read_end - b.read_start;
    if a_span != b_span {
        return a_span > b_span;
    }
    if a.seq_idx != b.seq_idx {
        return a.seq_idx < b.seq_idx;
    }
    if a.strand != b.strand {
        return a.strand < b.strand;
    }
    if a.read_start != b.read_start {
        return a.read_start < b.read_start;
    }
    if a.read_end != b.read_end {
        return a.read_end < b.read_end;
    }
    if a.seq_start != b.seq_start {
        return a.seq_start < b.seq_start;
    }
    a.seq_end < b.seq_end
}

/// Ported from `SeqSet::ExtendOverlap` (`SeqSet.hpp:1994-2100`): extends a
/// k-mer-chained `overlap` to the read's/allele's full available overhangs
/// on both ends via `AlignAlgo::GlobalAlignment`, stopping early at the
/// first `N` found in the reference overhang (mirrors the C++'s "locate the
/// boundary in the reference with 'N'" loops, `SeqSet.hpp:2008-2016,2043-
/// 2051`).
///
/// `read` is the (possibly already reverse-complemented, matching
/// `overlap.strand`) read sequence; `allele_consensus` is
/// `refSet.GetSeqConsensus(overlap.seq_idx)`.
///
/// Returns `Some(extended)` if the extended similarity meets
/// `ref_seq_similarity` (mirrors the C++ `ret == 1`, i.e.
/// `extendedOverlap.similarity >= refSeqSimilarity`, `SeqSet.hpp:2074-2075`),
/// `None` otherwise (`ret == 0`) -- callers still receive the extended
/// overlap either way via the
/// `Result`-independent output; unlike the C++ (which writes into
/// `extendedOverlap` unconditionally and returns a separate `int` status),
/// this port folds both into a single `Option` since every caller in this
/// port's scope (`assign_read`) only ever wants the value when it passed.
///
/// # `ignoreNonExonDiff` is not a parameter here
///
/// The C++ `ExtendOverlap` itself never reads `ignoreNonExonDiff` -- that
/// flag only affects the CALLER's (`AssignRead`'s) post-processing after
/// `ExtendOverlap` returns (see [`assign_read`]'s doc comment). This
/// function therefore has no `ignoreNonExonDiff`-shaped parameter to begin
/// with; it is included in this doc comment only to make that absence
/// explicit for a reader cross-checking against the module's other
/// `ignoreNonExonDiff == false`-only ports.
///
/// # Panics
///
/// Panics if `read.len()`/`allele_consensus.len()` do not fit an `i32`, or
/// if `overlap`'s coordinates are out of bounds for `read`/`allele_consensus`
/// (not expected: `overlap` is assumed to be a valid
/// [`crate::ref_kmer_filter::RefKmerFilter::get_overlaps_from_read`] result
/// against the same `read`/`allele_consensus`).
#[must_use]
#[allow(clippy::similar_names)]
pub fn extend_overlap(
    read: &[u8],
    allele_consensus: &[u8],
    overlap: &overlap::Overlap,
    ref_seq_similarity: f64,
    dp_cache: &mut crate::align_algo::DpCache,
) -> Option<ExtendedOverlap> {
    let len = i32::try_from(read.len()).expect("read length fits in i32");
    let consensus_len =
        i32::try_from(allele_consensus.len()).expect("consensus length fits in i32");

    // Extension to the 5' end (left).
    let mut left_overhang_size = overlap.read_start.min(overlap.seq_start);
    let mut left_clip = if overlap.read_start > overlap.seq_start {
        overlap.read_start - overlap.seq_start
    } else {
        0
    };
    for i in 0..left_overhang_size {
        let pos = usize::try_from(overlap.seq_start - i - 1).expect("in-bounds seq offset");
        if allele_consensus[pos] == b'N' {
            left_clip = left_overhang_size - i;
            left_overhang_size = i;
            break;
        }
    }

    let left_seq_start = usize::try_from(overlap.seq_start - left_overhang_size).unwrap();
    let left_seq_end = usize::try_from(overlap.seq_start).unwrap();
    let left_read_start = usize::try_from(overlap.read_start - left_overhang_size).unwrap();
    let left_read_end = usize::try_from(overlap.read_start).unwrap();
    let left_align = crate::align_algo::global_alignment_cached(
        &allele_consensus[left_seq_start..left_seq_end],
        &read[left_read_start..left_read_end],
        crate::align_algo::DEFAULT_BAND,
        dp_cache,
    );
    let (mut match_cnt, mut mismatch_cnt, mut indel_cnt) = (0, 0, 0);
    crate::align_algo::get_align_stats(
        &left_align.align,
        false,
        &mut match_cnt,
        &mut mismatch_cnt,
        &mut indel_cnt,
    );

    // Extension to the 3' end (right).
    let mut right_overhang_size =
        (len - 1 - overlap.read_end).min(consensus_len - 1 - overlap.seq_end);
    let mut right_clip = if len - 1 - overlap.read_end > consensus_len - 1 - overlap.seq_end {
        len - 1 - overlap.read_end - (consensus_len - 1 - overlap.seq_end)
    } else {
        0
    };
    for i in 0..right_overhang_size {
        let pos = usize::try_from(overlap.seq_end + 1 + i).expect("in-bounds seq offset");
        if allele_consensus[pos] == b'N' {
            right_clip = right_overhang_size - i;
            right_overhang_size = i;
            break;
        }
    }

    let right_seq_start = usize::try_from(overlap.seq_end + 1).unwrap();
    let right_seq_end = usize::try_from(overlap.seq_end + 1 + right_overhang_size).unwrap();
    let right_read_start = usize::try_from(overlap.read_end + 1).unwrap();
    let right_read_end = usize::try_from(overlap.read_end + 1 + right_overhang_size).unwrap();
    let right_align = crate::align_algo::global_alignment_cached(
        &allele_consensus[right_seq_start..right_seq_end],
        &read[right_read_start..right_read_end],
        crate::align_algo::DEFAULT_BAND,
        dp_cache,
    );
    // `update = true`: accumulate onto the left-extension's stats (SeqSet.hpp:2057).
    crate::align_algo::get_align_stats(
        &right_align.align,
        true,
        &mut match_cnt,
        &mut mismatch_cnt,
        &mut indel_cnt,
    );
    let _ = mismatch_cnt; // matches stock: computed, never read after this point.

    let mut extended = ExtendedOverlap {
        seq_idx: overlap.seq_idx,
        read_start: overlap.read_start - left_overhang_size,
        read_end: overlap.read_end + right_overhang_size,
        seq_start: overlap.seq_start - left_overhang_size,
        seq_end: overlap.seq_end + right_overhang_size,
        strand: overlap.strand,
        match_cnt: 2 * match_cnt + overlap.match_cnt,
        similarity: 0.0,
        relaxed_match_cnt: 0,
        left_clip,
        right_clip,
    };
    extended.similarity = f64::from(extended.match_cnt)
        / f64::from(
            extended.read_end - extended.read_start + 1 + extended.seq_end - extended.seq_start + 1,
        );
    extended.relaxed_match_cnt = extended.match_cnt;

    // `ret` (SeqSet.hpp:2074-2075): decided from the PRE-clip-adjustment
    // similarity and never revisited, even though `extendedOverlap.similarity`
    // itself IS mutated below when a clip is present -- ported verbatim
    // (not a bug in this port; matches stock exactly).
    let passed = extended.similarity >= ref_seq_similarity;

    if left_clip > 0 || right_clip > 0 {
        extended.match_cnt += 2 * left_clip + 2 * right_clip;
        extended.similarity = f64::from(extended.match_cnt)
            / f64::from(
                extended.read_end - extended.read_start + 1 + extended.seq_end - extended.seq_start
                    + 1
                    + 2 * left_clip
                    + 2 * right_clip,
            );
    }

    passed.then_some(extended)
}

/// Ported from the `nucToNum` table as linked into the `genotyper`/`analyzer`
/// binaries (`Genotyper.cpp:37-40`): maps `A`/`C`/`G`/`T` to a 0-3 index.
/// Duplicated from (rather than reusing) `align_algo`'s private, identically
/// defined `nuc_to_num` -- see this module's doc comment for why 5c stays
/// additive-only on top of Phase 3/3b/4 files rather than exposing their
/// internals.
fn nuc_to_num(c: u8) -> Option<usize> {
    match c {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Complements a single base: `A<->T`, `C<->G`, `N->N` -- mirrors
/// `numToNuc[3 - nucToNum[c - 'A']]` with `N` bypassing the table
/// (`SeqSet::ReverseComplement`, `SeqSet.hpp:2103-2114`). Duplicated from
/// (rather than reusing) `ref_kmer_filter`'s private, identically defined
/// `reverse_complement_into`/`complement_base` -- see this module's doc
/// comment for why 5c stays additive-only on top of Phase 3/3b/4 files.
fn complement_base(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// Returns the reverse complement of `seq` (see [`complement_base`]).
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| complement_base(c)).collect()
}

/// Per-allele reference bookkeeping this port's `AssignRead`/
/// `GetSeqMissingBaseCoverage` ports need beyond what
/// [`crate::ref_kmer_filter::RefKmerFilter`] tracks: exon intervals (for
/// `isValidDiff.exon`) and per-position base-coverage counts (`posWeight`).
/// Mirrors the fields of `_seqWrapper` (`SeqSet.hpp:25-42`) this port's 5c
/// slice actually reads/writes; a caller (the `unum genotype` CLI) builds
/// one of these per deduplicated allele sequence and keeps it alongside a
/// [`crate::ref_kmer_filter::RefKmerFilter`] built from the same
/// deduplicated sequences in the same order.
#[derive(Debug, Clone)]
pub struct AlleleRef {
    /// `_seqWrapper::consensus`.
    pub consensus: Vec<u8>,
    /// `_seqWrapper::exons` -- 0-based, inclusive `(start, end)` intervals,
    /// as parsed by [`parse_exon_comment`] from the reference FASTA header
    /// comment (or a single whole-sequence interval when no comment is
    /// present, matching `InputRefSeq`'s `comment == NULL` fallback,
    /// `SeqSet.hpp:970-976`).
    pub exons: Vec<(i32, i32)>,
    /// `_seqWrapper::posWeight` -- per-position base-observation counts,
    /// mutated in place by [`assign_read`]'s coverage-marking step (mirrors
    /// `SeqSet::AssignRead`'s `weight > 0` branch, `SeqSet.hpp:2253-2274`).
    /// Initialized to all-zero (`SetSeqExonInfo`'s `posWeight.SetZero`,
    /// `SeqSet.hpp:646`).
    ///
    /// Uses the interior-mutable [`crate::align_algo::AtomicPosWeight`] so the
    /// base-coverage marking (`count[code] += weight`, an order-independent
    /// integer add) can run across `rayon` workers via a SHARED `&[AlleleRef]`
    /// in [`assign_reads_parallel`]'s fused pass, byte-identically to the
    /// sequential loop (see `AtomicPosWeight`'s doc comment) and without
    /// per-thread copies of these RSS-dominating structures.
    pub pos_weight: Vec<crate::align_algo::AtomicPosWeight>,
    /// `_seqWrapper::separator` -- computed once at load time, identically in
    /// both `InputRefFa` (`SeqSet.hpp:~896`) and `InputRefSeq` (`SeqSet.hpp:
    /// ~913`): `[-1, <index of each 'N' byte>, seq_len]`, in that
    /// order. Consumed by [`is_separator_in_range`], which ports
    /// `SeqSet::IsSeparatorInRange` (`SeqSet.hpp:487-498`).
    pub separator: Vec<i32>,
}

impl AlleleRef {
    /// Builds a fresh [`AlleleRef`] for `consensus`, parsing `comment` (the
    /// FASTA header's post-id text, or `None`) into exon intervals via
    /// [`parse_exon_comment`] -- mirrors `InputRefSeq` + `SetSeqExonInfo`
    /// (`SeqSet.hpp:906-989,638-720`) for the `initExonInfo = true` call
    /// site `InitRefSet` always uses (`Genotyper.hpp:722`).
    #[must_use]
    pub fn new(consensus: Vec<u8>, comment: Option<&str>) -> Self {
        let len = consensus.len();
        let exons = match comment {
            Some(c) => parse_exon_comment(c, len),
            None => vec![(0, i32::try_from(len).unwrap_or(0) - 1)],
        };
        let separator = compute_separator(&consensus);
        Self {
            consensus,
            exons,
            pos_weight: (0..len).map(|_| crate::align_algo::AtomicPosWeight::default()).collect(),
            separator,
        }
    }

    /// `true` if position `pos` (0-based) falls within any parsed exon
    /// interval (`_validDiff::exon`, `SeqSet.hpp:657,670`).
    #[must_use]
    fn is_exon(&self, pos: i32) -> bool {
        self.exons.iter().any(|&(s, e)| pos >= s && pos <= e)
    }
}

/// Ported from the `separator`-building loop shared, IDENTICALLY, by
/// `SeqSet::InputRefFa` and `SeqSet::InputRefSeq` (`SeqSet.hpp:~896` and
/// `~913` respectively):
/// ```c
/// sw.separator.push_back(-1) ;
/// for (i = 0 ; sw.consensus[i] ; ++i)
///     if (sw.consensus[i] == 'N')
///         sw.separator.push_back(i) ;
/// sw.separator.push_back(i) ;   // i == seqLen here (loop exits at the NUL)
/// ```
/// i.e. `[-1, <index of each 'N' byte, in ascending order>, seq_len]`.
#[must_use]
pub fn compute_separator(consensus: &[u8]) -> Vec<i32> {
    // No pre-sized capacity: a `.filter().count()` bytecount pass just to
    // size this small, once-per-allele Vec is not worth the extra full scan
    // clippy's `naive_bytecount` lint is really guarding against (that lint
    // targets bytecount used as an END in itself, not an allocation hint).
    let mut separator = Vec::new();
    separator.push(-1);
    for (i, &b) in consensus.iter().enumerate() {
        if b == b'N' {
            separator.push(i32::try_from(i).unwrap_or(i32::MAX));
        }
    }
    separator.push(i32::try_from(consensus.len()).unwrap_or(i32::MAX));
    separator
}

/// Ported from `SeqSet::IsSeparatorInRange(int s, int e, int seqIdx)`
/// (`SeqSet.hpp:487-498`): `true` iff any entry of `separator` (the
/// `seqIdx`-th allele's [`AlleleRef::separator`], as built by
/// [`compute_separator`]) falls within the inclusive range `[s, e]`.
///
/// ```c
/// for each p in seqs[seqIdx].separator: if (p >= s && p <= e) return 1;
/// return 0;
/// ```
/// A plain linear scan matches the C++ exactly; `separator` is tiny (two
/// sentinels plus one entry per `N` byte -- including any at the sequence
/// endpoints) for every allele this port targets.
#[must_use]
pub fn is_separator_in_range(separator: &[i32], s: i32, e: i32) -> bool {
    separator.iter().any(|&p| p >= s && p <= e)
}

/// Ported from `InputRefSeq`'s comment-parsing loop (`SeqSet.hpp:936-961`):
/// splits `comment` into a sequence of non-negative integers (digit runs
/// separated by any non-digit byte), discards the FIRST number (an unused
/// leading marker in every reference this port targets), then pairs up the
/// rest into `(start, end)` inclusive exon intervals: `(nums[1], nums[2])`,
/// `(nums[3], nums[4])`, etc. -- i.e. `for i in (1..nums.len()).step_by(2)`.
///
/// Returns a single whole-sequence interval `[(0, seq_len - 1)]` if `comment`
/// yields fewer than 2 numbers (mirrors the C++'s `if (size > 0) {...} else
/// {whole-sequence interval}` -- `size` there is `nums.size()`, which is 0
/// only when `comment` has no digits at all; a `size` of exactly 1, i.e. only
/// the discarded leading number, still enters the `for (i = 1 ; i < size ...)`
/// loop as a zero-iteration range, yielding an EMPTY `exons` vec, not the
/// whole-sequence fallback -- ported exactly: `size <= 1` produces an empty
/// list here too, `size == 0` is the only case mapped to the fallback).
#[must_use]
pub fn parse_exon_comment(comment: &str, seq_len: usize) -> Vec<(i32, i32)> {
    let mut nums: Vec<i32> = Vec::new();
    let mut n: i32 = 0;
    let mut in_digits = false;
    for b in comment.bytes() {
        if b.is_ascii_digit() {
            n = n * 10 + i32::from(b - b'0');
            in_digits = true;
        } else if in_digits {
            nums.push(n);
            n = 0;
            in_digits = false;
        }
    }
    if in_digits {
        nums.push(n);
    }

    if nums.is_empty() {
        return vec![(0, i32::try_from(seq_len).unwrap_or(0) - 1)];
    }

    let mut exons = Vec::new();
    let mut i = 1;
    while i + 1 < nums.len() {
        exons.push((nums[i], nums[i + 1]));
        i += 2;
    }
    exons
}

/// Ported from `SeqSet::AssignRead` (`SeqSet.hpp:2119-2303`), restricted to
/// `ignoreNonExonDiff == false` (see this module's doc comment). Turns one
/// (deduplicated) read sequence into its sorted, extended overlap list,
/// mutating `allele_refs`' `pos_weight` in place when `weight > 0` (mirrors
/// the C++'s `weight > 0` base-coverage-marking gate,
/// `SeqSet.hpp:2253`) -- every caller in this port's scope
/// (`Genotyper.cpp:472`'s `weight = j - i`, the duplicate-sequence count)
/// passes `weight >= 1`, so that gate is always live here.
///
/// `overlaps` is `GetOverlapsFromRead`'s already-computed, already-similarity-
/// filtered result for `read` (i.e.
/// [`crate::ref_kmer_filter::RefKmerFilter::get_overlaps_from_read`]'s
/// output) -- this port takes it as a parameter rather than calling
/// `get_overlaps_from_read` itself, so it stays independent of exactly how a
/// caller obtains it (mirrors the "caller supplies reference-derived input"
/// pattern used throughout this module).
///
/// Returns `None` if `overlaps` is empty (mirrors the C++'s `overlapCnt <= 0`
/// early `return -1`, `SeqSet.hpp:2135-2138` -- `seqs.size() == 0` is
/// unreachable here since `allele_refs` is non-empty by construction in every
/// caller).
///
/// # `AssignRead`'s OWN two `IsSeparatorInRange` calls (fixes #40)
///
/// `AssignRead` (`SeqSet.hpp:2163-2169`) calls `IsSeparatorInRange` twice,
/// independently of [`Genotyper::set_read_assignments`]'s
/// `IsFragmentSpanSeparator` filter (see [`is_separator_in_range`]'s doc
/// comment, which #39 fixed): a per-overlap skip (`continue` when
/// `[seqStart, seqEnd]` spans a separator, `SeqSet.hpp:2163-2164`) and a
/// `needClip` computation (`SeqSet.hpp:2166-2169`) that feeds
/// `onlyConsiderClip`'s `similarity < 0.95` escape hatch
/// (`SeqSet.hpp:2170-2172`). This function ports BOTH, inline in the loop
/// below. For `kir_rna_seq.fa` (no interior `N` runs, so `separator` is
/// always exactly `[-1, len]`) both are a no-op, matching the oracle exactly.
/// For references with interior `N`s (e.g. the HLA DNA reference #39 was
/// diagnosed against) neither is a no-op: `similarity` for overlaps reaching
/// this loop is NOT always `0.0` (a real, `matchCnt`-derived ratio survives
/// `GetOverlapsFromRead`'s own `>= refSeqSimilarity` filter, `SeqSet.hpp:
/// 1836-1845,1896`, typically landing in `[0.8, 1.0]` for ref alleles since
/// `radius = 10` by default means per-overlap indels don't force it to `0`
/// either, `SeqSet.hpp:1754,1763`), so `needClip`'s `similarity < 0.95`
/// escape hatch is a genuinely reachable condition, not dead code -- both
/// calls are ported for fidelity.
///
/// # Panics
///
/// Panics if any `overlaps[i].seq_idx` is out of range for `allele_refs`
/// (not expected: `overlaps` is assumed to come from
/// [`crate::ref_kmer_filter::RefKmerFilter::get_overlaps_from_read`] against
/// the same reference set `allele_refs` was built from).
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn assign_read(
    read: &[u8],
    overlaps: &[overlap::Overlap],
    allele_refs: &[AlleleRef],
    ref_seq_similarity: f64,
    weight: i32,
    dp_cache: &mut crate::align_algo::DpCache,
) -> Option<Vec<ExtendedOverlap>> {
    if overlaps.is_empty() {
        return None;
    }

    let mut sorted_overlaps = overlaps.to_vec();
    sorted_overlaps.sort_by(|a, b| {
        if overlap::overlap_less_than(a, b) {
            std::cmp::Ordering::Less
        } else if overlap::overlap_less_than(b, a) {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        }
    });

    let rc: Vec<u8> = reverse_complement(read);
    let r: &[u8] = if sorted_overlaps[0].strand == -1 { &rc } else { read };

    let mut extended_overlaps: Vec<ExtendedOverlap> = Vec::new();
    let mut only_consider_clip = false;
    let mut good_match_cnt: i32 = -1;

    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
    let len = r.len() as i32;
    for o in &sorted_overlaps {
        let allele_idx = usize::try_from(o.seq_idx).expect("seq_idx is non-negative");
        // SeqSet.hpp:2163-2172 (fixes #40): see this function's doc comment
        // ("AssignRead's OWN two IsSeparatorInRange calls") for why both
        // calls are ported and why `needClip` is NOT dead code here.
        let separator = &allele_refs[allele_idx].separator;
        if is_separator_in_range(separator, o.seq_start, o.seq_end) {
            continue;
        }
        let need_clip = is_separator_in_range(
            separator,
            o.seq_start - o.read_start,
            o.seq_end + (len - o.read_end - 1),
        );
        if only_consider_clip && o.match_cnt < good_match_cnt && (!need_clip || o.similarity < 0.95)
        {
            continue;
        }
        let consensus = &allele_refs[allele_idx].consensus;
        if let Some(extended) = extend_overlap(r, consensus, o, ref_seq_similarity, dp_cache) {
            extended_overlaps.push(extended);
            if !only_consider_clip && (good_match_cnt == -1 || o.match_cnt > good_match_cnt) {
                good_match_cnt = o.match_cnt;
            }
        } else {
            only_consider_clip = true;
        }
    }

    // Adding the exonic and base support for good overlaps (SeqSet.hpp:2188-2286).
    if !extended_overlaps.is_empty() && weight >= 0 {
        let mut best = extended_overlaps[0];
        for &eo in &extended_overlaps {
            if extended_overlap_less_than(&eo, &best) {
                best = eo;
            }
        }

        for eo in &mut extended_overlaps {
            if eo.match_cnt < best.match_cnt - 10 {
                // "The assignment is very bad" (SeqSet.hpp:2261-2262): C++'s
                // `else` arm zeroes `relaxedMatchCnt` here (rather than
                // leaving it at `extend_overlap`'s `relaxedMatchCnt =
                // matchCnt` default) and skips the alignment/base-coverage
                // marking below entirely.
                eo.relaxed_match_cnt = 0;
                continue;
            }
            let allele_idx = usize::try_from(eo.seq_idx).expect("seq_idx is non-negative");
            let seq_start = usize::try_from(eo.seq_start).unwrap();
            let seq_end = usize::try_from(eo.seq_end).unwrap();
            let read_start = usize::try_from(eo.read_start).unwrap();
            let read_end = usize::try_from(eo.read_end).unwrap();
            let align = crate::align_algo::global_alignment_cached(
                &allele_refs[allele_idx].consensus[seq_start..=seq_end],
                &r[read_start..=read_end],
                crate::align_algo::DEFAULT_BAND,
                dp_cache,
            );

            // `ignoreNonExonDiff == false`: `relaxedMatchCnt = matchCnt`
            // (SeqSet.hpp:2247-2250), skipping the exon-walk that would
            // otherwise recompute it.
            eo.relaxed_match_cnt = eo.match_cnt;

            // Mark the base coverage (SeqSet.hpp:2252-2274). This is the ONLY
            // shared mutation `assign_read` performs and the ONLY reason
            // `assign_read` had to run sequentially. It is an order-independent
            // integer add (`count[code] += weight`), so `pos_weight` uses
            // `AtomicPosWeight` (`add` -> `fetch_add(_, Relaxed)`): the fused
            // parallel pass in `assign_reads_parallel` can mark coverage across
            // `rayon` workers via a shared `&[AlleleRef]` and still land the
            // exact same totals as the sequential loop (see `AtomicPosWeight`).
            // `assign_read` never READS `pos_weight` (of this or any allele)
            // during its own execution -- it only writes here -- so there is no
            // read-during-write race; all `pos_weight` reads happen after the
            // pass joins (`get_seq_missing_base_coverage`).
            if weight > 0 {
                let mut ref_pos = eo.seq_start;
                let mut read_pos = eo.read_start;
                for &op in &align.align {
                    if op == crate::align_algo::EDIT_MATCH {
                        let rb = r[usize::try_from(read_pos).unwrap()];
                        if rb != b'N' {
                            if let Some(code) = nuc_to_num(rb) {
                                let pos = usize::try_from(ref_pos).unwrap();
                                allele_refs[allele_idx].pos_weight[pos].add(code, weight);
                            }
                        }
                    }
                    if op != crate::align_algo::EDIT_INSERT {
                        ref_pos += 1;
                    }
                    if op != crate::align_algo::EDIT_DELETE {
                        read_pos += 1;
                    }
                }
            }
        }
    }

    // `extendedOverlaps.size() > 1000` truncation (SeqSet.hpp:2290-2298).
    if extended_overlaps.len() > 1000 {
        extended_overlaps.sort_by(|a, b| {
            if extended_overlap_less_than(a, b) {
                std::cmp::Ordering::Less
            } else if extended_overlap_less_than(b, a) {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        });
        let first_similarity = extended_overlaps[0].similarity;
        let mut j = 1usize;
        while j < extended_overlaps.len()
            && extended_overlaps[j].similarity >= first_similarity - 0.1
        {
            j += 1;
        }
        extended_overlaps.truncate(j);
    }

    Some(extended_overlaps)
}

/// Runs [`RefKmerFilter::get_overlaps_from_read`] followed by [`assign_read`]
/// for every distinct read sequence in `sorted_seqs` (mirrors
/// `Genotyper.cpp:463-480`/`AssignReads_Thread`'s per-distinct-sequence read
/// loop -- see `crates/unum/src/stages/genotype.rs`'s `run()` for the
/// caller that sorts/dedups `all_seqs` into `sorted_seqs` and computes
/// `weight_of`), returning one `Option<Vec<ExtendedOverlap>>` per
/// `sorted_seqs` entry, in the SAME order as `sorted_seqs` (index `i` of the
/// result corresponds to `sorted_seqs[i]`).
///
/// # Parallelism: one fused, order-preserving work-stealing pass
///
/// `threads <= 1` runs the plain sequential loop (no `rayon` pool is built at
/// all): for each `seq`, call `get_overlaps_from_read` then immediately
/// `assign_read`, exactly as the single-threaded reference driver always has.
///
/// `threads > 1` runs `get_overlaps_from_read` AND `assign_read` for each
/// `sorted_seqs[i]` inside ONE `par_iter().map_init(...)` pass, with no
/// barrier between them (an earlier design ran a parallel `get_overlaps`
/// phase, `collect()`ed at a barrier, then a fully-sequential `assign_read`
/// phase -- that barrier left ~24% of thread-time idle at `-t4`; fusing the
/// two halves into one work-stealing pass removes it). Details:
///
/// - **Per-worker state, never shared.** Each `rayon` worker threads its own
///   [`Scratch`] (for `get_overlaps_from_read`'s hit list / bucket map / RC
///   buffer) AND its own [`crate::align_algo::DpCache`] DP memo via
///   `map_init`. The DP cache stays a pure per-worker memo -- exactly as it
///   was when it lived in the old sequential phase -- so memoization changes
///   nothing observable => byte-identical.
/// - **`allele_refs` shared `&[AlleleRef]`.** The only shared mutation is
///   `assign_read`'s `pos_weight` base-coverage marking, done through
///   [`crate::align_algo::AtomicPosWeight`] (`add` ->
///   `fetch_add(_, Relaxed)`). That marking is an order-independent integer
///   add, so the FINAL accumulated coverage is identical no matter which
///   worker's add lands first (see `AtomicPosWeight`'s doc comment).
///   `assign_read` never READS `pos_weight` mid-pass (it only writes there),
///   and `get_overlaps_from_read` uses only local zero-count `PosWeight`
///   slices (never `allele_refs.pos_weight`), so there is no read-during-write
///   race anywhere in the fused body. Every real `pos_weight` read
///   (`get_seq_missing_base_coverage`) happens only after this pass joins.
/// - **Order-preserving.** `collect()` on this `IndexedParallelIterator`
///   (`(&[&[u8]]).par_iter()`) preserves input order, so the result `Vec`
///   lines up with `sorted_seqs` regardless of which worker computed which
///   entry or in what order they finished.
///
/// Together these make `-t N` output byte-identical to `-t 1` at any `N`
/// (proven by the `parallel_threads_1_4_8` differential): the coverage totals
/// are thread-count-invariant, and no other shared or global/thread-local
/// mutable state is touched (`get_overlaps_from_read`'s `AlignAlgo`-based
/// refinement allocates fresh local buffers per call -- there is no shared
/// alignment scratch object beyond the per-worker `Scratch`/`DpCache`).
///
/// # Panics
///
/// Panics if `threads > 1` and the OS refuses to spawn `rayon`'s worker
/// threads for the scoped pool (mirrors this module's existing convention of
/// `.expect()`-ing on invariants/environment preconditions rather than
/// threading `anyhow::Result` through a pure-computation module -- see e.g.
/// [`extend_overlap`]'s/[`assign_read`]'s own panicking-on-invariant-violation
/// doc comments).
#[must_use]
pub fn assign_reads_parallel(
    filter: &RefKmerFilter,
    sorted_seqs: &[&[u8]],
    allele_refs: &[AlleleRef],
    ref_seq_similarity: f64,
    weight_of: impl Fn(&[u8]) -> i32 + Sync,
    threads: usize,
) -> Vec<Option<Vec<ExtendedOverlap>>> {
    if threads <= 1 {
        let mut scratch = Scratch::default();
        let mut dp_cache = crate::align_algo::DpCache::new();
        return sorted_seqs
            .iter()
            .map(|&seq| {
                let raw = filter.get_overlaps_from_read(seq, &mut scratch).unwrap_or_default();
                assign_read(
                    seq,
                    &raw,
                    allele_refs,
                    ref_seq_similarity,
                    weight_of(seq),
                    &mut dp_cache,
                )
            })
            .collect();
    }

    // One fused, order-preserving parallel pass: each `rayon` worker computes
    // `get_overlaps_from_read` AND then `assign_read` for the sequences it is
    // handed, with NO barrier between the two (removes the ~24%-of-thread-time
    // rayon idle that the old two-phase `collect()`-then-sequential-loop
    // structure spent stalled at the barrier). Per-worker state (`Scratch` for
    // `get_overlaps_from_read`'s hit list/bucket map/RC buffer, plus the
    // `DpCache` DP memo) is threaded via `map_init`, so it is never shared
    // across workers -- exactly as the DP cache was per-worker before, keeping
    // it a pure per-worker memo => byte-identical results.
    //
    // `allele_refs` is now shared `&[AlleleRef]`: `assign_read`'s only shared
    // mutation is the `pos_weight` base-coverage marking, which is an
    // order-independent integer add done through `AtomicPosWeight` (see its doc
    // comment and `assign_read`'s marking loop). No other shared state is
    // touched by either half of the fused body, and no thread ever reads
    // `pos_weight` mid-pass, so the accumulated coverage is thread-count-
    // invariant -- `-t N` output is byte-identical to `-t 1` at any `N`.
    //
    // `collect()` on this `IndexedParallelIterator` (`(&[&[u8]]).par_iter()`)
    // preserves input order, so the result `Vec` lines up with `sorted_seqs`
    // regardless of which worker computed which entry or in what order they
    // finished (`extended_by_seq` stays in `sorted_seqs` order downstream).
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("building rayon thread pool for parallel read assignment");
    pool.install(|| {
        sorted_seqs
            .par_iter()
            .map_init(
                || (Scratch::default(), crate::align_algo::DpCache::new()),
                |(scratch, dp_cache), &seq| {
                    let raw = filter.get_overlaps_from_read(seq, scratch).unwrap_or_default();
                    assign_read(
                        seq,
                        &raw,
                        allele_refs,
                        ref_seq_similarity,
                        weight_of(seq),
                        dp_cache,
                    )
                },
            )
            .collect()
    })
}

/// Ported from `SeqSet::GetSeqMissingBaseCoverage` (`SeqSet.hpp:2717-2755`):
/// counts how many of `allele_ref`'s exon positions (sorted ascending by
/// observed base coverage) fall below `ratio` times the MEDIAN exon-position
/// coverage (floored at `1`) -- i.e. the number of "missing" (low-coverage)
/// exon bases, consumed as [`AlleleInfo::missing_coverage`] by
/// [`Genotyper::finalize_read_assignments`].
///
/// # Panics
///
/// Panics if `allele_ref` has zero exon positions (mirrors the C++'s
/// undefined behavior indexing `exonBaseCoverage[k / 2]` with `k == 0`) --
/// not expected for any real allele, which always has at least one exon
/// interval (see [`AlleleRef::new`]'s whole-sequence fallback).
#[must_use]
pub fn get_seq_missing_base_coverage(allele_ref: &AlleleRef, ratio: f64) -> i32 {
    let mut exon_base_coverage: Vec<i32> = Vec::new();
    for (i, &base) in allele_ref.consensus.iter().enumerate() {
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let i_i32 = i as i32;
        if allele_ref.is_exon(i_i32) {
            let code = nuc_to_num(base).expect("consensus base is A/C/G/T");
            exon_base_coverage.push(allele_ref.pos_weight[i].get(code));
        }
    }
    exon_base_coverage.sort_unstable();
    let k = exon_base_coverage.len();
    assert!(k > 0, "get_seq_missing_base_coverage: allele has zero exon positions");

    let mut cutoff = f64::from(exon_base_coverage[k / 2]) * ratio;
    if cutoff < 1.0 {
        cutoff = 1.0;
    }

    let mut i = 0i32;
    for &cov in &exon_base_coverage {
        if f64::from(cov) >= cutoff {
            break;
        }
        i += 1;
    }
    i
}

/// Ported from `SeqSet::IsOverlapIntersect` (`SeqSet.hpp:317-323`): `true` if
/// `a`/`b` are against the same allele and their `[seqStart, seqEnd]` ranges
/// overlap. Not currently reachable from
/// [`read_assignment_to_fragment_assignment`] (see that function's doc
/// comment for why -- the `ignoreNonExonDiff`-gated branch that would call
/// this is dead at this port's `ignoreNonExonDiff == false` default); ported
/// anyway for fidelity, matching this module's general "port the branch even
/// if currently dead at the default flags" convention (see
/// [`extend_overlap`]'s doc comment for the same rationale applied
/// elsewhere).
#[must_use]
pub fn is_overlap_intersect(a: &ExtendedOverlap, b: &ExtendedOverlap) -> bool {
    a.seq_idx == b.seq_idx
        && ((a.seq_start <= b.seq_start && a.seq_end >= b.seq_start)
            || (b.seq_start <= a.seq_start && b.seq_end >= a.seq_start))
}

/// Ported from `SeqSet::TruncatedMatePairOverlap` (`SeqSet.hpp:502-522`):
/// `o` is a mate-1-set overlap; tests whether `o`'s mate pair could be
/// missing because the reference allele sequence is too short to contain it
/// (i.e. the implied mate position falls off the end of `consensus_len`).
///
/// `separator_in_range` mirrors `IsSeparatorInRange` (see
/// [`is_separator_in_range`]) -- callers pass a closure over each allele's
/// [`AlleleRef::separator`] (built by [`compute_separator`]), matching
/// [`Genotyper::set_read_assignments`]'s own `separator_lookup` pattern. For
/// `kir_rna_seq.fa` (no interior `N`s) `separator` is always exactly
/// `[-1, len]`, so this is a no-op there; for references with interior `N`s
/// it is live.
///
/// # Panics
///
/// Panics if `o.seq_idx` does not fit an `i32` (not expected: allele indices
/// are always small and non-negative).
#[must_use]
pub fn truncated_mate_pair_overlap(
    o: &ExtendedOverlap,
    comp_mate1: &ExtendedOverlap,
    comp_mate2: &ExtendedOverlap,
    consensus_len: i32,
    mut separator_in_range: impl FnMut(i32, i32, i32) -> bool,
) -> bool {
    if o.strand == 1 {
        // Mate is on the right.
        consensus_len - 1 < o.seq_end + comp_mate2.seq_end - comp_mate1.seq_end
            || separator_in_range(
                o.seq_end,
                o.seq_end + comp_mate2.seq_end - comp_mate1.seq_end + 1,
                i32::try_from(o.seq_idx).expect("seq_idx fits in i32"),
            )
    } else {
        o.seq_start - (comp_mate1.seq_start - comp_mate2.seq_start) < 0
            || separator_in_range(
                o.seq_start - (comp_mate1.seq_start - comp_mate2.seq_start) - 1,
                o.seq_start,
                i32::try_from(o.seq_idx).expect("seq_idx fits in i32"),
            )
    }
}

/// Internal pairing of a [`FragmentOverlap`] with the raw `ExtendedOverlap`(s)
/// it was built from -- mirrors `_fragmentOverlap::overlap1`/`overlap2`
/// (`SeqSet.hpp:158-159`), which the C++ keeps attached to every
/// `_fragmentOverlap` all the way through
/// `ReadAssignmentToFragmentAssignment` (used by its
/// [`fragment_overlap_less_than`] tie-break and the truncated-mate-pair
/// filter) even though [`FragmentOverlap`] itself (5a/5b's struct, which
/// only stores what [`Genotyper::set_read_assignments`] needs downstream)
/// does not carry them. Kept as a transient, function-local pairing here
/// rather than widening [`FragmentOverlap`]'s public fields.
#[derive(Debug, Clone, Copy)]
struct Assembled {
    fragment: FragmentOverlap,
    overlap1: ExtendedOverlap,
    overlap2: Option<ExtendedOverlap>,
}

/// Ported from `SeqSet::ReadAssignmentToFragmentAssignment` (`SeqSet.hpp:
/// 2310-2655`): combines a read pair's two mate-end [`assign_read`] outputs
/// into per-allele [`FragmentOverlap`]s -- the direct input to
/// [`Genotyper::set_read_assignments`].
///
/// `overlaps1`/`overlaps2` mirror `pOverlaps1`/`pOverlaps2` (mate 1's / mate
/// 2's [`assign_read`] output; `overlaps2 = None` for single-end reads,
/// mirroring `pOverlaps2 == NULL`). `has_n` mirrors the C++'s
/// `reads1[i].hasN || reads2[i].hasN` (computed by the caller, since this
/// port has no `_genotypeRead`-equivalent read struct).
///
/// # `ignoreNonExonDiff == false`: the WHOLE relaxed-tie qual=1 branch is dead
///
/// C++'s tie-marking loop (`SeqSet.hpp:2492-2518`) has two `if`/`else if`
/// arms that can set `assign[k].qual = 1`: the `exactTie` arm (matchCnt AND
/// similarity both equal `bestAssign`'s, `SeqSet.hpp:2508`) and a second
/// `else if (ignoreNonExonDiff && assign[i].matchCnt >= bestAssign.matchCnt -
/// matchCntRelax && assign[i].relaxedMatchCnt == bestAssign.relaxedMatchCnt)`
/// arm (`SeqSet.hpp:2514`). That second arm's condition is gated on
/// `ignoreNonExonDiff` as a WHOLE -- not just the `matchCntRelax` widening
/// from `2` to `4` computed just above it (`SeqSet.hpp:2495-2505`, the
/// `IsOverlapIntersect` branch). At this port's `ignoreNonExonDiff == false`
/// default (see this module's doc comment; there is no `--relaxIntronAlign`
/// flag in this port's scope), the ENTIRE `else if` arm is unreachable in
/// C++ -- so only `exactTie` ever sets `qual = 1` here. This function
/// therefore does not compute a `relaxed_tie` term at all: doing so
/// unconditionally (i.e. not gated on a false `ignoreNonExonDiff`) would
/// admit fragment assignments C++ would drop, corrupting downstream
/// abundance/genotype calls. A future port of `--relaxIntronAlign` should
/// reintroduce this arm gated behind an explicit `ignore_non_exon_diff: bool`
/// parameter (mirroring `matchCntRelax`'s own dead `IsOverlapIntersect` arm,
/// which for the same reason is hardcoded to `2` below rather than computed).
///
/// `hit_len_required` mirrors `refSet.hitLenRequired` (the dangling-mate-pair
/// filter's span threshold, `SeqSet.hpp:2561`; see
/// [`crate::ref_kmer_filter::RefKmerFilter::hit_len_required`]).
/// `consensus_len_of(seq_idx)` mirrors `refSet.GetSeqConsensusLen`/
/// `seqs[seqIdx].consensusLen` (needed by [`truncated_mate_pair_overlap`]'s
/// filter, `SeqSet.hpp:2580-2653`); `separator_in_range` mirrors
/// `IsSeparatorInRange` (see [`is_separator_in_range`]) and is live for
/// references with interior `N`s (a no-op for `kir_rna_seq.fa`, which has
/// none).
///
/// # Panics
///
/// Panics if any `ExtendedOverlap::seq_idx` does not fit an `i32`/`u32` (not
/// expected: every allele index this port produces is small and
/// non-negative).
#[must_use]
pub fn read_assignment_to_fragment_assignment(
    overlaps1: &[ExtendedOverlap],
    overlaps2: Option<&[ExtendedOverlap]>,
    has_n: bool,
    hit_len_required: i32,
    consensus_len_of: impl Fn(u32) -> i32,
    separator_in_range: impl FnMut(i32, i32, i32) -> bool,
) -> Vec<FragmentOverlap> {
    read_assignment_to_fragment_assignment_impl(
        overlaps1,
        overlaps2,
        has_n,
        hit_len_required,
        consensus_len_of,
        separator_in_range,
    )
    .into_iter()
    .map(|a| a.fragment)
    .collect()
}

/// Same port as [`read_assignment_to_fragment_assignment`], but additionally
/// returns each kept fragment's underlying mate-end [`ExtendedOverlap`](s) it
/// was built from (`_fragmentOverlap::overlap1`/`overlap2`,
/// `SeqSet.hpp:158-159`) alongside the [`FragmentOverlap`] itself.
///
/// [`read_assignment_to_fragment_assignment`] (5a/5b's entry point) discards
/// this `overlap1`/`overlap2` state once fragment assembly finishes, since
/// [`Genotyper::set_read_assignments`] never reads it. Task 6b's `Analyzer`
/// driver DOES need it: `_fragmentOverlap::overlap1`/`overlap2` (the
/// per-mate-end coordinates/strand) are exactly what
/// `SeqSet::AddOverlapAlignmentInfo`/`AddFragmentAlignmentInfo`
/// (`SeqSet.hpp:2657-2680,2758-2778`) populate an alignment op-sequence onto
/// before `VariantCaller::ComputeVariant` runs -- see
/// [`crate::variant_caller::add_fragment_alignment_info`]. Rather than
/// widening [`FragmentOverlap`]'s public fields (which would ripple into
/// every 5a/5b caller that has no use for them), this is a second, additive
/// entry point over the exact same shared implementation.
#[must_use]
pub fn read_assignment_to_fragment_assignment_with_overlaps(
    overlaps1: &[ExtendedOverlap],
    overlaps2: Option<&[ExtendedOverlap]>,
    has_n: bool,
    hit_len_required: i32,
    consensus_len_of: impl Fn(u32) -> i32,
    separator_in_range: impl FnMut(i32, i32, i32) -> bool,
) -> Vec<(FragmentOverlap, ExtendedOverlap, Option<ExtendedOverlap>)> {
    read_assignment_to_fragment_assignment_impl(
        overlaps1,
        overlaps2,
        has_n,
        hit_len_required,
        consensus_len_of,
        separator_in_range,
    )
    .into_iter()
    .map(|a| (a.fragment, a.overlap1, a.overlap2))
    .collect()
}

/// Shared implementation behind [`read_assignment_to_fragment_assignment`]
/// and [`read_assignment_to_fragment_assignment_with_overlaps`] -- see the
/// former's doc comment for the full C++ citation; this function differs
/// from that doc comment only in its `Vec<Assembled>` return type (both
/// public wrappers project out exactly what they need from it).
#[allow(clippy::too_many_lines, clippy::similar_names)]
fn read_assignment_to_fragment_assignment_impl(
    overlaps1: &[ExtendedOverlap],
    overlaps2: Option<&[ExtendedOverlap]>,
    has_n: bool,
    hit_len_required: i32,
    consensus_len_of: impl Fn(u32) -> i32,
    mut separator_in_range: impl FnMut(i32, i32, i32) -> bool,
) -> Vec<Assembled> {
    // Build the (mate1Idx, mate2Idx) fragment index pairs (SeqSet.hpp:2314-2384).
    // `-1` (C++'s sentinel for "no pairing on this side") is modeled as `None`.
    let mut fragments: Vec<(Option<usize>, Option<usize>)> = Vec::new();
    match overlaps2 {
        None => {
            for i in 0..overlaps1.len() {
                fragments.push((Some(i), None));
            }
        }
        Some(overlaps2) if overlaps1.is_empty() || overlaps2.is_empty() => {
            for i in 0..overlaps1.len() {
                fragments.push((Some(i), None));
            }
            for i in 0..overlaps2.len() {
                fragments.push((None, Some(i)));
            }
        }
        Some(overlaps2) => {
            let mut seq_idx_to_overlap2: HashMap<u32, Vec<usize>> = HashMap::new();
            for (i, o2) in overlaps2.iter().enumerate() {
                seq_idx_to_overlap2.entry(o2.seq_idx).or_default().push(i);
            }
            for (i, o1) in overlaps1.iter().enumerate() {
                let Some(candidates) = seq_idx_to_overlap2.get(&o1.seq_idx) else { continue };
                for &j in candidates {
                    let o2 = &overlaps2[j];
                    // Compatible mate pairs: opposite strands, same allele,
                    // and o1 positioned "before" o2 on the strand.
                    if o1.strand == o2.strand || o1.seq_idx != o2.seq_idx {
                        continue;
                    }
                    if (o1.strand == 1 && o1.seq_start < o2.seq_start)
                        || (o1.strand == -1 && o1.seq_start > o2.seq_start)
                    {
                        fragments.push((Some(i), Some(j)));
                    }
                }
            }
        }
    }

    // For each seq idx, keep the best fragment (SeqSet.hpp:2386-2455).
    let mut assign: Vec<Assembled> = Vec::new();
    let mut seq_idx_to_assign_idx: HashMap<u32, usize> = HashMap::new();
    for &(a, b) in &fragments {
        let assembled = match (a, b) {
            (Some(a), b_opt) => {
                let o = overlaps1[a];
                let mut fo = FragmentOverlap {
                    seq_idx: i32::try_from(o.seq_idx).unwrap(),
                    seq_start: o.seq_start,
                    seq_end: o.seq_end,
                    match_cnt: o.match_cnt,
                    relaxed_match_cnt: o.relaxed_match_cnt,
                    similarity: o.similarity,
                    has_mate_pair: false,
                    o1_from_r2: false,
                    qual: 0.0,
                    has_n,
                };
                let mut overlap2 = None;
                if let Some(b) = b_opt {
                    let o2 = overlaps2.unwrap()[b];
                    fo.match_cnt += o2.match_cnt;
                    fo.relaxed_match_cnt += o2.relaxed_match_cnt;
                    if o.strand == 1 {
                        fo.seq_end = o2.seq_end;
                    } else {
                        fo.seq_start = o2.seq_start;
                    }
                    fo.similarity = f64::from(fo.match_cnt)
                        / f64::from(
                            (o.read_end - o.read_start + 1)
                                + (o2.read_end - o2.read_start + 1)
                                + (o.seq_end - o.seq_start + 1)
                                + (o2.seq_end - o2.seq_start + 1)
                                + 2 * o.left_clip
                                + 2 * o.right_clip
                                + 2 * o2.left_clip
                                + 2 * o2.right_clip,
                        );
                    fo.has_mate_pair = true;
                    overlap2 = Some(o2);
                }
                Assembled { fragment: fo, overlap1: o, overlap2 }
            }
            (None, Some(b)) => {
                // Dangling case: only mate 2 aligned; `overlap1` (the
                // "primary" overlap the C++ keeps for tie-breaking/filtering)
                // is set from mate 2's overlap here, matching
                // `fragmentOverlap.overlap1 = o` at `SeqSet.hpp:2439`.
                let o = overlaps2.unwrap()[b];
                let fo = FragmentOverlap {
                    seq_idx: i32::try_from(o.seq_idx).unwrap(),
                    seq_start: o.seq_start,
                    seq_end: o.seq_end,
                    match_cnt: o.match_cnt,
                    relaxed_match_cnt: o.relaxed_match_cnt,
                    similarity: o.similarity,
                    has_mate_pair: false,
                    o1_from_r2: true,
                    qual: 0.0,
                    has_n,
                };
                Assembled { fragment: fo, overlap1: o, overlap2: None }
            }
            (None, None) => continue,
        };

        let seq_idx_u32 = u32::try_from(assembled.fragment.seq_idx).unwrap();
        if let Some(&existing_idx) = seq_idx_to_assign_idx.get(&seq_idx_u32) {
            // Note `<` here is for ranking, so smaller has higher rank.
            if fragment_overlap_less_than(&assembled, &assign[existing_idx]) {
                assign[existing_idx] = assembled;
            }
        } else {
            assign.push(assembled);
            seq_idx_to_assign_idx.insert(seq_idx_u32, assign.len() - 1);
        }
    }

    // Pick the best assignment and mark qual=1 ties (SeqSet.hpp:2474-2545).
    let mut best_assign = FragmentOverlap {
        seq_idx: 0,
        seq_start: 0,
        seq_end: 0,
        match_cnt: -1,
        relaxed_match_cnt: 0,
        similarity: 0.0,
        has_mate_pair: false,
        o1_from_r2: false,
        qual: 0.0,
        has_n: false,
    };
    for a in &assign {
        #[allow(clippy::float_cmp)]
        let tie_similarity_better = a.fragment.match_cnt == best_assign.match_cnt
            && a.fragment.similarity > best_assign.similarity;
        if a.fragment.match_cnt > best_assign.match_cnt || tie_similarity_better {
            best_assign = a.fragment;
        }
    }

    // Only `exactTie` (SeqSet.hpp:2508) ever fires at `ignoreNonExonDiff ==
    // false` -- see this function's doc comment for why the `relaxedTie`
    // `else if` arm (SeqSet.hpp:2514) is not ported here at all.
    let mut kept: Vec<Assembled> = Vec::new();
    for a in &assign {
        #[allow(clippy::float_cmp)]
        let exact_tie = a.fragment.match_cnt == best_assign.match_cnt
            && a.fragment.similarity == best_assign.similarity;
        if exact_tie {
            let mut a = *a;
            a.fragment.qual = 1.0;
            kept.push(a);
        }
    }
    assign = kept;

    // Dangling-mate-pair filter (SeqSet.hpp:2553-2578): only applies when
    // paired-end AND the best assignment has NO mate pair.
    if !assign.is_empty() && overlaps2.is_some() && !assign[0].fragment.has_mate_pair {
        let mut i = 0usize;
        while i < assign.len() {
            let a = &assign[i];
            let span = a.fragment.seq_end - a.fragment.seq_start
                + 1
                + (a.overlap1.read_end - a.overlap1.read_start + 1);
            if a.fragment.similarity < 1.0
                || separator_in_range(a.fragment.seq_start, a.fragment.seq_end, a.fragment.seq_idx)
                || span < 3 * hit_len_required
            {
                break;
            }
            let span_range = 100;
            let consensus_len = consensus_len_of(u32::try_from(a.fragment.seq_idx).unwrap());
            if (a.overlap1.strand == 1 && a.fragment.seq_end + span_range < consensus_len)
                || (a.overlap1.strand == -1 && a.fragment.seq_start - span_range >= 0)
            {
                break;
            }
            i += 1;
        }
        if i < assign.len() {
            assign.clear();
        }
    }

    // Truncated-mate-pair filter (SeqSet.hpp:2580-2653): only applies when
    // paired-end AND the best assignment HAS a mate pair.
    if !assign.is_empty() && overlaps2.is_some() && assign[0].fragment.has_mate_pair {
        #[allow(clippy::float_cmp)] // exact C++ `== 1` on qual, not a fuzzy float comparison.
        let representative_idx = assign.iter().position(|a| a.fragment.qual == 1.0).unwrap_or(0);
        let representative = &assign[representative_idx];
        let rep_overlap1 = representative.overlap1;
        let rep_overlap2 = representative.overlap2.expect("has_mate_pair implies overlap2 is set");

        let mut filter = false;

        // "Better read 1": is there an unused mate-1 overlap that dominates
        // the representative's overlap1 AND looks truncated?
        for o in overlaps1 {
            if filter {
                break;
            }
            // SeqSet.hpp:2601-2609: the "unused" guard (`seqIdxToOverlapIdx.find
            // == end`) binds ONLY to the `matchCnt ==` branch, NOT the
            // `matchCnt >` branch -- an overlap that strictly dominates on
            // matchCnt triggers this filter even if it is already assigned.
            // (This crate previously gated the whole `dominates` on `unused`,
            // wrongly sparing already-assigned dominating overlaps and
            // over-keeping fragments -- e.g. HLA-U U*01:03 over-assignment.)
            let dominates = o.match_cnt > rep_overlap1.match_cnt
                || (o.match_cnt == rep_overlap1.match_cnt
                    && o.similarity > rep_overlap1.similarity
                    && !seq_idx_to_assign_idx.contains_key(&o.seq_idx));
            if dominates
                && (truncated_mate_pair_overlap(
                    o,
                    &rep_overlap1,
                    &rep_overlap2,
                    consensus_len_of(o.seq_idx),
                    &mut separator_in_range,
                ) || o.similarity > rep_overlap2.similarity + 0.1)
            {
                filter = true;
            }
        }

        // "Better read 2": same check against mate-2 overlaps.
        if !filter {
            if let Some(overlaps2) = overlaps2 {
                for o in overlaps2 {
                    if filter {
                        break;
                    }
                    // SeqSet.hpp:2627-2635: `unused` guard binds ONLY to the
                    // `matchCnt ==` branch, not the `matchCnt >` branch (see the
                    // "better read 1" block above for the full rationale).
                    let dominates = o.match_cnt > rep_overlap2.match_cnt
                        || (o.match_cnt == rep_overlap2.match_cnt
                            && o.similarity > rep_overlap2.similarity
                            && !seq_idx_to_assign_idx.contains_key(&o.seq_idx));
                    if dominates
                        && (truncated_mate_pair_overlap(
                            o,
                            &rep_overlap2,
                            &rep_overlap1,
                            consensus_len_of(o.seq_idx),
                            &mut separator_in_range,
                        ) || o.similarity > rep_overlap1.similarity + 0.1)
                    {
                        filter = true;
                    }
                }
            }
        }

        if filter {
            assign.clear();
        }
    }

    assign
}

/// Ported from `_fragmentOverlap::operator<` (`SeqSet.hpp:161-172`): the
/// "smaller is better" ranking [`read_assignment_to_fragment_assignment`]
/// uses to keep the best fragment overlap per allele. The final `overlap1 <
/// b.overlap1` tie-break (`SeqSet.hpp:169`) uses [`extended_overlap_less_than`]
/// on each [`Assembled`]'s `overlap1` (see that struct's doc comment for why
/// this state is tracked alongside [`FragmentOverlap`] rather than inside
/// it).
#[must_use]
fn fragment_overlap_less_than(a: &Assembled, b: &Assembled) -> bool {
    if a.fragment.match_cnt != b.fragment.match_cnt {
        return a.fragment.match_cnt > b.fragment.match_cnt;
    }
    #[allow(clippy::float_cmp)]
    let similarity_differs = a.fragment.similarity != b.fragment.similarity;
    if similarity_differs {
        return a.fragment.similarity > b.fragment.similarity;
    }
    extended_overlap_less_than(&a.overlap1, &b.overlap1)
}

/// The Phase-5a slice of T1K's `Genotyper` (`Genotyper.hpp`): the
/// allele/gene data model, its initialization, and read-assignment weight +
/// storage. See the module docs for exactly what is (and is not) ported.
#[derive(Debug, Clone)]
pub struct Genotyper {
    // --- allele/gene/majorAllele bookkeeping (Genotyper.hpp:449-461) ---
    /// `Genotyper::alleleInfo` (a `SimpleVector<_alleleInfo>` in C++; a plain
    /// `Vec` here since [`Genotyper::init_allele_info`] always builds the
    /// whole thing up front via `ExpandTo`, matching a `Vec::resize`).
    pub allele_info: Vec<AlleleInfo>,
    /// `Genotyper::majorAlleleNameToIdx`.
    pub major_allele_name_to_idx: HashMap<String, i32>,
    /// `Genotyper::geneNameToIdx`.
    pub gene_name_to_idx: HashMap<String, i32>,
    /// `Genotyper::majorAlleleSize` -- summed `refSet.GetSeqWeight(i)` per
    /// major allele.
    pub major_allele_size: Vec<i32>,
    /// `Genotyper::geneIdxToName`.
    pub gene_idx_to_name: Vec<String>,
    /// `Genotyper::majorAlleleIdxToName`.
    pub major_allele_idx_to_name: Vec<String>,
    /// `Genotyper::geneCnt`.
    pub gene_cnt: i32,
    /// `Genotyper::majorAlleleCnt`.
    pub major_allele_cnt: i32,
    /// `Genotyper::alleleCnt`.
    pub allele_cnt: i32,
    /// `Genotyper::geneSimilarity` -- a `geneCnt x geneCnt` matrix,
    /// `geneSimilarity[i][j] = kmerProfiles[i].GetCountSimilarity(kmerProfiles[j])`
    /// (asymmetric; `geneSimilarity[i][i] == 1.0` always).
    pub gene_similarity: Vec<Vec<f64>>,

    // --- allele-name parsing config (Genotyper.hpp:481-482, 506-507) ---
    /// `Genotyper::alleleDigitUnits`. Constructor default `-1`
    /// (`Genotyper.hpp:506`).
    pub allele_digit_units: i32,
    /// `Genotyper::alleleDelimiter`. Constructor default `'\0'`
    /// (`Genotyper.hpp:507`).
    pub allele_delimiter: u8,

    // --- read-assignment storage (Genotyper.hpp:439-447) ---
    /// `Genotyper::readCnt` -- the coalesced read-group count. NOT updated
    /// by [`Genotyper::set_read_assignments`] (only `CoalesceReadAssignments`,
    /// Phase 5b, advances it); [`Genotyper::init_read_assignments`] resets it
    /// to `0`, matching `InitReadAssignments` (`Genotyper.hpp:762`).
    pub read_cnt: i32,
    /// `Genotyper::totalReadCnt`.
    pub total_read_cnt: i32,
    /// `Genotyper::maxAssignCnt`. Constructor default `2000`
    /// (`Genotyper.hpp:493`).
    pub max_assign_cnt: i32,
    /// `Genotyper::allReadAssignments` -- per-raw-read (not yet coalesced)
    /// allele assignments, indexed by `readId`. Populated by
    /// [`Genotyper::set_read_assignments`].
    pub all_read_assignments: Vec<Vec<ReadAssignment>>,

    // --- Phase-5b: coalesced read assignments + equivalent classes
    // (Genotyper.hpp:442-446) ---
    /// `Genotyper::readsInAllele` -- per-allele list of `(readIdx,
    /// idxWithinReadAssignments[readIdx])` pairs (`_pair`'s `(a, b)`).
    /// Populated by [`Genotyper::finalize_read_assignments`].
    pub reads_in_allele: Vec<Vec<(i32, i32)>>,
    /// `Genotyper::readAssignments` -- the coalesced (deduplicated) per-
    /// read-group allele assignments, indexed `0..read_cnt`. Populated by
    /// [`Genotyper::coalesce_read_assignments`].
    pub read_assignments: Vec<Vec<ReadAssignment>>,
    /// `Genotyper::readAssignmentsFingerprintToIdx` -- maps a coalesced
    /// read-group's allele-index fingerprint to the `read_assignments`
    /// indices sharing that fingerprint (a hash-collision bucket;
    /// [`Genotyper::coalesce_read_assignments`] still does a full
    /// [`Genotyper::is_read_assignment_the_same`] check within the bucket).
    pub read_assignments_fingerprint_to_idx: HashMap<i32, Vec<i32>>,
    /// `Genotyper::equivalentClassToAlleles` -- per-equivalent-class list of
    /// member allele indices. Populated by
    /// [`Genotyper::build_allele_equivalent_class`]; EC index order is
    /// load-bearing (consumed by [`Genotyper::quantify_allele_equivalent_class`]
    /// and Phase 5c's allele selection).
    pub equivalent_class_to_alleles: Vec<Vec<i32>>,

    // --- Phase-5b: EM tuning knobs (Genotyper.hpp:472-477, 509) ---
    /// `Genotyper::filterFrac`. Constructor default `0.15`
    /// (`Genotyper.hpp:498`); read by
    /// [`Genotyper::quantify_allele_equivalent_class`]'s periodic
    /// low-abundance masking step.
    pub filter_frac: f64,
    /// `Genotyper::minSquaremAlpha`. Constructor default `0`
    /// (`Genotyper.hpp:509`); the C++ CLI (`Genotyper.cpp:main`) never calls
    /// `SetMinSquaremAlpha`, so this stays `0` in the stock pipeline -- the
    /// clamp that reads it (`Genotyper.hpp:1243-1244`,
    /// [`Genotyper::quantify_allele_equivalent_class`]) is dead code at that
    /// default (`minSquaremAlpha < 0` is false), ported here verbatim anyway
    /// for byte-identical behavior if a caller does set it non-default.
    pub min_squarem_alpha: f64,

    // --- Phase-5b: abundance outputs (Genotyper.hpp:464-466) ---
    /// `Genotyper::geneAbundance`. Populated by
    /// [`Genotyper::set_allele_abundance`].
    pub gene_abundance: Vec<f64>,
    /// `Genotyper::majorAlleleAbundance`. Populated by
    /// [`Genotyper::set_allele_abundance`].
    pub major_allele_abundance: Vec<f64>,
    /// `Genotyper::geneMaxMajorAlleleAbundance`. Populated by
    /// [`Genotyper::set_allele_abundance`].
    pub gene_max_major_allele_abundance: Vec<f64>,

    // --- Phase-5c: allele-selection tuning knobs (Genotyper.hpp:477-484,499-500,504) ---
    /// `Genotyper::filterCov`. Constructor default `1.0` (`Genotyper.hpp:499`);
    /// the minimum rank-level abundance [`Genotyper::select_alleles_for_genes`]
    /// requires before assigning a nonzero `genotypeQuality`.
    pub filter_cov: f64,
    /// `Genotyper::crossGeneRate`. Constructor default `0.04`
    /// (`Genotyper.hpp:500`); scales cross-gene `geneSimilarity`-weighted
    /// noise into each gene's null-hypothesis mean in
    /// [`Genotyper::select_alleles_for_genes`].
    pub cross_gene_rate: f64,
    /// `Genotyper::readLength`. Constructor default `0` (`Genotyper.hpp:504`,
    /// set via `SetReadLength`, `Genotyper.cpp:443`, from the max observed
    /// read length across the input FASTQ(s)); used by
    /// [`Genotyper::select_alleles_for_genes`]'s missing-coverage weight
    /// term.
    pub read_length: i32,

    // --- Phase-5c: allele-selection output (Genotyper.hpp:447) ---
    /// `Genotyper::selectedAlleles` -- per-gene list of `(alleleIdx,
    /// alleleRank)` pairs (`_pair`'s `(a, b)`). Populated by
    /// [`Genotyper::select_alleles_for_genes`]; consumed by
    /// [`Genotyper::get_gene_allele_types`]/[`Genotyper::get_allele_description`]/
    /// [`Genotyper::output_representative_alleles`].
    pub selected_alleles: Vec<Vec<(i32, i32)>>,
}

impl Default for Genotyper {
    fn default() -> Self {
        Self::new()
    }
}

impl Genotyper {
    /// Ported from the `Genotyper(int kmerLength)` constructor
    /// (`Genotyper.hpp:488-510`), restricted to the fields this Phase-5a
    /// slice actually uses. `kmerLength` itself is not stored here: it only
    /// feeds `refSet`'s construction in the C++ (`refSet(kmerLength)`,
    /// `Genotyper.hpp:488`), and `refSet`/`SeqSet` is out of this struct's
    /// scope (see [`Genotyper::init_allele_info`]'s doc comment for how
    /// callers supply the reference-derived inputs `InitAlleleInfo` would
    /// otherwise pull from `refSet` directly).
    #[must_use]
    pub fn new() -> Self {
        Self {
            allele_info: Vec::new(),
            major_allele_name_to_idx: HashMap::new(),
            gene_name_to_idx: HashMap::new(),
            major_allele_size: Vec::new(),
            gene_idx_to_name: Vec::new(),
            major_allele_idx_to_name: Vec::new(),
            gene_cnt: 0,
            major_allele_cnt: 0,
            allele_cnt: 0,
            gene_similarity: Vec::new(),
            allele_digit_units: -1,
            allele_delimiter: b'\0',
            read_cnt: 0,
            total_read_cnt: 0,
            max_assign_cnt: 2000,
            all_read_assignments: Vec::new(),
            reads_in_allele: Vec::new(),
            read_assignments: Vec::new(),
            read_assignments_fingerprint_to_idx: HashMap::new(),
            equivalent_class_to_alleles: Vec::new(),
            filter_frac: 0.15,
            min_squarem_alpha: 0.0,
            gene_abundance: Vec::new(),
            major_allele_abundance: Vec::new(),
            gene_max_major_allele_abundance: Vec::new(),
            filter_cov: 1.0,
            cross_gene_rate: 0.04,
            read_length: 0,
            selected_alleles: Vec::new(),
        }
    }

    /// Ported from `Genotyper::SetFilterFrac` (`Genotyper.hpp:528-531`).
    pub fn set_filter_frac(&mut self, f: f64) {
        self.filter_frac = f;
    }

    /// Ported from `Genotyper::SetFilterCov` (`Genotyper.hpp:533-536`).
    pub fn set_filter_cov(&mut self, c: f64) {
        self.filter_cov = c;
    }

    /// Ported from `Genotyper::SetCrossGeneRate` (`Genotyper.hpp:543-546`).
    pub fn set_cross_gene_rate(&mut self, r: f64) {
        self.cross_gene_rate = r;
    }

    /// Ported from `Genotyper::SetReadLength` (`Genotyper.hpp:538-541`).
    pub fn set_read_length(&mut self, rl: i32) {
        self.read_length = rl;
    }

    /// Ported from `Genotyper::SetMinSquaremAlpha` (`Genotyper.hpp:554-557`).
    pub fn set_min_squarem_alpha(&mut self, a: f64) {
        self.min_squarem_alpha = a;
    }

    /// Ported from `Genotyper::SetAlleleNameStructure` (`Genotyper.hpp:548-552`).
    pub fn set_allele_name_structure(&mut self, n: i32, d: u8) {
        self.allele_digit_units = n;
        self.allele_delimiter = d;
    }

    /// Calls [`parse_allele_name`] with this `Genotyper`'s configured
    /// `allele_digit_units`/`allele_delimiter` (mirrors the private member
    /// function `Genotyper::ParseAlleleName`, `Genotyper.hpp:63-131`, which
    /// reads those same two fields off `this`).
    #[must_use]
    pub fn parse_allele_name(&self, allele: &str, fields_type: i32) -> (String, String) {
        parse_allele_name(allele, fields_type, self.allele_digit_units, self.allele_delimiter)
    }

    /// Ported from `Genotyper::InitAlleleInfo` (`Genotyper.hpp:559-682`).
    ///
    /// Unlike the C++ (which reads `refSet.GetSeqName`/`GetSeqConsensus`/
    /// `GetSeqWeight`/`GetSeqEffectiveLen`/`SetSeqEffectiveLen` directly off
    /// its own `refSet: SeqSet` member), this port takes the same four
    /// per-sequence inputs as explicit parallel slices -- `seq_names`
    /// (`GetSeqName`), `seq_consensus` (`GetSeqConsensus`), `seq_weight`
    /// (`GetSeqWeight`), and `seq_effective_len` (`GetSeqEffectiveLen`/
    /// `SetSeqEffectiveLen`, mutated in place to match the C++ "adjust
    /// allele effective length" step, `Genotyper.hpp:641-681`) -- so this
    /// module has no dependency on a full `SeqSet`/`RefKmerFilter` port.
    /// `refSet.Size()` becomes `seq_names.len()`. A caller wiring this up
    /// against [`crate::ref_kmer_filter::RefKmerFilter`] supplies
    /// `seq_names[i] = ref_kmer_filter.seq_name(i)`,
    /// `seq_consensus[i] = ref_kmer_filter.seq_consensus(i)`, and an
    /// `effective_len`/`weight` computed the same way `SeqSet::InputRefSeq`
    /// does (`SeqSet.hpp:906-982`, itself out of this port's scope).
    ///
    /// The k-mer length used for the gene-similarity `KmerCount` profiles
    /// (`Genotyper.hpp:599`, `KmerCount kmerProfiles[geneCnt]`, using
    /// `KmerCount`'s own default constructor -- `KmerCount()`,
    /// `KmerCount.hpp:36-41` -- which fixes `kmerLength = 31`) is NOT the
    /// same as `refSet`'s alignment k-mer length; `kmer_profile_k` here
    /// mirrors that same `KmerCount()` default-constructor value (callers
    /// should pass `31`).
    ///
    /// # Panics
    ///
    /// Panics if the four input slices have different lengths, or if
    /// `seq_names[i]` is not valid UTF-8-sliceable (see
    /// [`parse_allele_name`]'s panic doc).
    pub fn init_allele_info(
        &mut self,
        seq_names: &[String],
        seq_consensus: &[Vec<u8>],
        seq_weight: &[i32],
        seq_effective_len: &mut [i32],
        kmer_profile_k: usize,
    ) {
        // Threshold for the "adjust allele effective length" step below
        // (Genotyper.hpp:672); hoisted to the top of the function (rather
        // than declared right before its one use site) only to satisfy
        // clippy::items_after_statements.
        const LARGE_DELETION: i32 = 500;

        assert!(
            seq_names.len() == seq_consensus.len()
                && seq_names.len() == seq_weight.len()
                && seq_names.len() == seq_effective_len.len(),
            "init_allele_info: seq_names/seq_consensus/seq_weight/seq_effective_len must all be \
             the same length"
        );

        self.allele_cnt = i32::try_from(seq_names.len()).expect("allele_cnt fits in i32");
        self.allele_info = vec![AlleleInfo::default(); seq_names.len()];

        for (i, name) in seq_names.iter().enumerate() {
            let (gene, major_allele) = self.parse_allele_name(name, 0);

            let gene_idx = *self.gene_name_to_idx.entry(gene.clone()).or_insert_with(|| {
                let idx = self.gene_cnt;
                self.gene_idx_to_name.push(gene.clone());
                self.gene_cnt += 1;
                idx
            });
            let major_allele_idx =
                *self.major_allele_name_to_idx.entry(major_allele.clone()).or_insert_with(|| {
                    let idx = self.major_allele_cnt;
                    self.major_allele_idx_to_name.push(major_allele.clone());
                    self.major_allele_size.push(0);
                    self.major_allele_cnt += 1;
                    idx
                });

            self.allele_info[i] =
                AlleleInfo { major_allele_idx, gene_idx, ..AlleleInfo::default() };
            let major_allele_idx_usize =
                usize::try_from(major_allele_idx).expect("major_allele_idx is non-negative");
            self.major_allele_size[major_allele_idx_usize] += seq_weight[i];
        }

        // Compute gene similarity (Genotyper.hpp:597-639): for each gene,
        // pick the lexicographically smallest consensus among its alleles
        // (as a UTF-8-agnostic BYTE comparison; identical to C++ `strcmp`
        // for the ASCII/DNA byte domain here -- consensus is uppercase ACGTN,
        // all < 128 with no interior NUL, so signed-char `strcmp` and Rust's
        // unsigned `&[u8]` ordering agree) to build that gene's k-mer profile, then
        // fill an asymmetric geneCnt x geneCnt similarity matrix.
        let gene_cnt = usize::try_from(self.gene_cnt).expect("gene_cnt is non-negative");
        let mut kmer_profiles: Vec<KmerCount> = Vec::with_capacity(gene_cnt);
        for gi in 0..gene_cnt {
            let gi32 = i32::try_from(gi).expect("gene index fits in i32");
            let mut min_tag: Option<usize> = None;
            for (j, info) in self.allele_info.iter().enumerate() {
                if info.gene_idx != gi32 {
                    continue;
                }
                min_tag = Some(match min_tag {
                    None => j,
                    Some(cur) => {
                        if seq_consensus[j] < seq_consensus[cur] {
                            j
                        } else {
                            cur
                        }
                    }
                });
            }
            let mut profile = KmerCount::new(kmer_profile_k);
            if let Some(tag) = min_tag {
                profile.add_count(&seq_consensus[tag]);
            }
            kmer_profiles.push(profile);
        }

        self.gene_similarity = vec![vec![0.0; gene_cnt]; gene_cnt];
        for i in 0..gene_cnt {
            for j in 0..gene_cnt {
                self.gene_similarity[i][j] =
                    if i == j { 1.0 } else { kmer_profiles[i].count_similarity(&kmer_profiles[j]) };
            }
        }

        // Adjust allele effective length (Genotyper.hpp:641-681): within
        // each gene, find the modal effective length among its alleles,
        // then bump any allele whose effective length is more than
        // `largeDeletion` (500) shorter than the mode up to the mode
        // itself (treats a large apparent deletion as a reference-building
        // artifact, not a real allele-defining feature).
        let mut gene_idx_to_allele_idx: Vec<Vec<usize>> = vec![Vec::new(); gene_cnt];
        for (i, info) in self.allele_info.iter().enumerate() {
            let gene_idx_usize = usize::try_from(info.gene_idx).expect("gene_idx is non-negative");
            gene_idx_to_allele_idx[gene_idx_usize].push(i);
        }
        for allele_ids in &gene_idx_to_allele_idx {
            let mut lens: Vec<i32> = allele_ids.iter().map(|&j| seq_effective_len[j]).collect();
            lens.sort_unstable();

            let size = lens.len();
            let mut len_mode = 0;
            let mut max = 0;
            let mut j = 0;
            while j < size {
                let mut k = j;
                while k < size && lens[k] == lens[j] {
                    k += 1;
                }
                if k - j > max {
                    max = k - j;
                    len_mode = lens[j];
                }
                j = k;
            }

            for &allele_idx in allele_ids {
                if seq_effective_len[allele_idx] < len_mode - LARGE_DELETION {
                    seq_effective_len[allele_idx] = len_mode;
                }
            }
        }
    }

    /// Ported from `Genotyper::ReadAssignmentWeight` (`Genotyper.hpp:205-230`).
    ///
    /// Pure deterministic ratio math (comparisons and one division), no
    /// `pow`/`log` -- bit-identical to the C++ oracle for the same inputs.
    /// `ref_seq_similarity` mirrors `refSet.GetRefSeqSimilarity()`
    /// (`Genotyper.hpp:212`); the commented-out `overlap2.similarity`
    /// fallback (`Genotyper.hpp:210-211`) and `pairedEndData`/`hasMatePair`
    /// adjustment (`Genotyper.hpp:223-224`) are dead code in stock T1K
    /// (`if`/`//` both commented out) and are correctly NOT ported here.
    #[must_use]
    pub fn read_assignment_weight(o: &FragmentOverlap, ref_seq_similarity: f64) -> f64 {
        let mut ret: f64 = 1.0;

        let similarity = o.similarity;
        let mut segment = (1.0 - ref_seq_similarity) / 4.0;
        if segment < 0.01 {
            segment = 0.01;
        }
        if similarity < 1.0 - 3.0 * segment {
            ret = 0.01;
        } else if similarity < 1.0 - 2.0 * segment {
            ret = 0.1;
        } else if similarity < 1.0 - segment {
            ret = 0.5;
        }

        if o.has_n {
            ret /= 10.0;
        }

        ret
    }

    /// Ported from `Genotyper::InitReadAssignments` (`Genotyper.hpp:759-777`).
    ///
    /// Also resets `readsInAllele`/`readAssignments`/
    /// `readAssignmentsFingerprintToIdx` (`Genotyper.hpp:764-767`),
    /// Phase-5b's `CoalesceReadAssignments`/`FinalizeReadAssignments`
    /// territory -- `readsInAllele` is sized to `allele_cnt` (matching
    /// `readsInAllele.resize(alleleCnt)`, `Genotyper.hpp:767`; requires
    /// [`Genotyper::init_allele_info`] to have already run, same
    /// precondition the C++ has).
    pub fn init_read_assignments(&mut self, total_read_cnt: i32, max_assign_cnt: i32) {
        self.max_assign_cnt = max_assign_cnt;
        self.read_cnt = 0;
        self.total_read_cnt = total_read_cnt;
        self.all_read_assignments = vec![Vec::new(); usize::try_from(total_read_cnt).unwrap_or(0)];
        self.read_assignments.clear();
        self.read_assignments_fingerprint_to_idx.clear();
        self.reads_in_allele = vec![Vec::new(); usize::try_from(self.allele_cnt).unwrap_or(0)];
    }

    /// Ported from `Genotyper::SetReadAssignments` (`Genotyper.hpp:778-832`).
    ///
    /// Stores `assignment`'s per-allele `FragmentOverlap`s (after
    /// whitelist-filtering and weight computation) into
    /// `self.all_read_assignments[read_id]`. Returns early (leaving that
    /// slot empty, matching the C++ bare `return;`) if:
    /// - `assignment.len() > max_assign_cnt` (when `max_assign_cnt > 0`,
    ///   `Genotyper.hpp:783-784`);
    /// - any entry's `(seq_idx, seq_start, seq_end)` span crosses a
    ///   reference-sequence separator (`separator_lookup`, `Genotyper.hpp:
    ///   798-802` via `refSet.IsFragmentSpanSeparator`).
    ///
    /// `separator_lookup(seq_idx, seq_start, seq_end)` mirrors
    /// `SeqSet::IsFragmentSpanSeparator`/`IsSeparatorInRange`
    /// (`SeqSet.hpp:2305-2308`, `487-498`): note the argument order here is
    /// `(seq_idx, s, e)` (matching `IsFragmentSpanSeparator`'s
    /// `IsSeparatorInRange(assign.seqStart, assign.seqEnd, assign.seqIdx)`
    /// call, just reordered to put `seq_idx` first) -- callers should build
    /// this closure from each allele's [`AlleleRef::separator`] (via
    /// [`is_separator_in_range`]), matching
    /// [`Genotyper::init_allele_info`]'s "caller supplies the
    /// reference-derived input" pattern.
    ///
    /// # Panics
    ///
    /// Panics if `read_id` is out of bounds for `self.all_read_assignments`
    /// (i.e. `>= total_read_cnt` as set by
    /// [`Genotyper::init_read_assignments`]) -- mirrors the C++'s undefined
    /// behavior (an out-of-bounds `std::vector` index) with a deterministic
    /// panic instead.
    pub fn set_read_assignments(
        &mut self,
        read_id: usize,
        assignment: &[FragmentOverlap],
        ref_seq_similarity: f64,
        mut separator_lookup: impl FnMut(i32, i32, i32) -> bool,
    ) {
        self.all_read_assignments[read_id].clear();

        let assignment_cnt = assignment.len();
        if self.max_assign_cnt > 0
            && i32::try_from(assignment_cnt).unwrap_or(i32::MAX) > self.max_assign_cnt
        {
            return;
        }

        for a in assignment {
            if separator_lookup(a.seq_idx, a.seq_start, a.seq_end) {
                return;
            }
        }

        let mut adjust_factor: f64 = 1.0;
        let mut max_similarity: f64 = 0.0;
        for a in assignment {
            if a.similarity > max_similarity {
                max_similarity = a.similarity;
            }
        }
        if max_similarity < 1.0 {
            adjust_factor = 0.25;
        }

        for a in assignment {
            let allele_idx = usize::try_from(a.seq_idx).expect("seq_idx must be non-negative");
            if !self.allele_info[allele_idx].whitelist {
                continue;
            }
            // `weight` is `f32` in the C++ struct (`_readAssignment::weight`),
            // so the `f64` `ReadAssignmentWeight` result is narrowed at
            // this exact assignment point, matching `Genotyper.hpp:827`.
            #[allow(clippy::cast_possible_truncation)]
            let weight = Self::read_assignment_weight(a, ref_seq_similarity) as f32;
            #[allow(clippy::cast_possible_truncation)]
            let na = ReadAssignment {
                allele_idx: a.seq_idx,
                start: a.seq_start,
                end: a.seq_end,
                weight,
                qual: a.qual as f32,
                adjust_weight: (adjust_factor as f32) * weight,
            };
            self.all_read_assignments[read_id].push(na);
        }
    }

    /// Ported from `Genotyper::IsReadAssignmentTheSame`
    /// (`Genotyper.hpp:181-196`): both vectors must have equal length and,
    /// element-by-element, equal `(alleleIdx, qual)` -- assumes both are
    /// already sorted the same way (by `alleleIdx`, the C++'s
    /// `_readAssignment::operator<`), matching the C++'s own "the read id in
    /// each vector should be sorted" precondition comment.
    #[must_use]
    #[allow(clippy::float_cmp)] // exact C++ `==` on qual, not a fuzzy float comparison.
    fn is_read_assignment_the_same(a1: &[ReadAssignment], a2: &[ReadAssignment]) -> bool {
        if a1.len() != a2.len() {
            return false;
        }
        a1.iter().zip(a2.iter()).all(|(x, y)| x.allele_idx == y.allele_idx && x.qual == y.qual)
    }

    /// Ported from `Genotyper::IsAssignedReadTheSame` (`Genotyper.hpp:164-179`):
    /// both `(readIdx, idxWithinReadAssignments)` lists must have equal
    /// length and, element-by-element, equal `readIdx` AND equal
    /// `read_assignments[readIdx][idxWithinReadAssignments].qual` --
    /// assumes both lists are already sorted by `readIdx` ("the read id in
    /// each vector should be sorted", matching the C++ comment).
    #[must_use]
    #[allow(clippy::float_cmp)] // exact C++ `==` on qual, not a fuzzy float comparison.
    fn is_assigned_read_the_same(&self, l1: &[(i32, i32)], l2: &[(i32, i32)]) -> bool {
        if l1.len() != l2.len() {
            return false;
        }
        l1.iter().zip(l2.iter()).all(|(&(a1, b1), &(a2, b2))| {
            if a1 != a2 {
                return false;
            }
            let qual = |read_idx: i32, within_idx: i32| {
                self.read_assignments[usize::try_from(read_idx).expect("readIdx is non-negative")]
                    [usize::try_from(within_idx).expect("within-idx is non-negative")]
                .qual
            };
            qual(a1, b1) == qual(a2, b2)
        })
    }

    /// Ported from `Genotyper::CoalesceReadAssignments(begin, end)`
    /// (`Genotyper.hpp:841-908`).
    ///
    /// For each read `i` in `[begin, end]` (also bounded by `total_read_cnt`,
    /// matching the C++'s own `&& i < totalReadCnt` loop guard) with a
    /// non-empty `all_read_assignments[i]`: sorts that read's assignments by
    /// `allele_idx` (mirrors `_readAssignment::operator<`, which orders
    /// solely on `allele_idx` -- see [`Genotyper::is_read_assignment_the_same`]'s
    /// doc comment for the sorted-input precondition this establishes for
    /// later `is_read_assignment_the_same`/`is_assigned_read_the_same`
    /// calls), computes a hash "fingerprint" of the sorted `allele_idx`
    /// sequence, then either merges into an existing read group with the
    /// identical (fingerprint-bucketed, then exactly verified) assignment
    /// set, or starts a new read group. Deterministic (integer fingerprint
    /// hash + `==` comparisons on already-computed `f32`/`f64` values, no
    /// new floating-point arithmetic) -> byte-identical to the C++ oracle.
    ///
    /// Frees `all_read_assignments[i]` for every `i` in range after
    /// processing (mirrors the C++'s `std::vector<...>().swap(...)` memory
    /// release, `Genotyper.hpp:903-906`) -- callers must not read
    /// `all_read_assignments` for coalesced reads afterward.
    ///
    /// Returns the count of reads in `[begin, end]` that had a non-empty
    /// assignment set (mirrors the C++ `ret`/`++ret`, `Genotyper.hpp:844,851`).
    ///
    /// # Panics
    ///
    /// Panics if any `allele_idx` is negative, or if `begin`/`end` cannot be
    /// converted to `usize` (both are expected to be non-negative read
    /// indices, matching the C++'s `int` loop bounds used as valid vector
    /// indices).
    pub fn coalesce_read_assignments(&mut self, begin: i32, end: i32) -> i32 {
        // `FINGERPRINT_MAX`, `Genotyper.hpp:847` (declared inside the C++
        // loop body, hoisted here since it is loop-invariant).
        const FINGERPRINT_MAX: i64 = 20_000_003;

        let mut ret = 0i32;
        let begin_usize = usize::try_from(begin).expect("begin is non-negative");
        // `end_inclusive` mirrors the C++ loop guard `i <= end && i <
        // totalReadCnt`; when the range is empty (`end_inclusive < begin`),
        // both the processing loop and the memory-release loop below
        // execute zero iterations, matching the C++ exactly.
        let end_inclusive = end.min(self.total_read_cnt - 1);
        let end_usize_exclusive = if end_inclusive < begin {
            begin_usize
        } else {
            usize::try_from(end_inclusive + 1).expect("end_inclusive + 1 is non-negative")
        };

        for i in begin_usize..end_usize_exclusive {
            let size = self.all_read_assignments[i].len();
            if size == 0 {
                continue;
            }
            ret += 1;

            self.all_read_assignments[i].sort_by_key(|a| a.allele_idx);

            let mut fingerprint: i64 = 0;
            for a in &self.all_read_assignments[i] {
                let k = i64::from(a.allele_idx);
                fingerprint = (fingerprint * i64::from(self.allele_cnt) + k) % FINGERPRINT_MAX;
            }
            // C++ stores `fingerprint` back into an `int` (implicit
            // truncating narrowing conversion at the `int fingerprint =`
            // declaration's subsequent assignments) before using it as the
            // `std::map<int, ...>` key -- replicate that exact truncation.
            #[allow(clippy::cast_possible_truncation)]
            let fingerprint = fingerprint as i32;

            let mut add_to: i32 = -1;
            if let Some(bucket) = self.read_assignments_fingerprint_to_idx.get(&fingerprint) {
                for &idx in bucket {
                    if Self::is_read_assignment_the_same(
                        &self.all_read_assignments[i],
                        &self.read_assignments[usize::try_from(idx).expect("idx is non-negative")],
                    ) {
                        add_to = idx;
                        break;
                    }
                }
            }

            if add_to == -1 {
                self.read_assignments.push(self.all_read_assignments[i].clone());
                self.read_assignments_fingerprint_to_idx
                    .entry(fingerprint)
                    .or_default()
                    .push(self.read_cnt);
                self.read_cnt += 1;
            } else {
                let add_to_usize = usize::try_from(add_to).expect("add_to is non-negative");
                for j in 0..size {
                    let src = self.all_read_assignments[i][j];
                    #[allow(clippy::float_cmp)]
                    let qual_is_one = src.qual == 1.0;
                    if qual_is_one {
                        if src.start < self.read_assignments[add_to_usize][j].start {
                            self.read_assignments[add_to_usize][j].start = src.start;
                        }
                        // Ported verbatim from `Genotyper.hpp:893-894`: the
                        // C++ compares `.end < .end` but then assigns
                        // `.end = allReadAssignments[i][j].start` (NOT
                        // `.end`) -- a real quirk in stock T1K, not a typo
                        // in this port. Preserved for byte-identical output.
                        if src.end < self.read_assignments[add_to_usize][j].end {
                            self.read_assignments[add_to_usize][j].end = src.start;
                        }
                    }
                    self.read_assignments[add_to_usize][j].weight += src.weight;
                    self.read_assignments[add_to_usize][j].adjust_weight += src.adjust_weight;
                }
            }
        }

        // Release the memory space for allReadAssignments (Genotyper.hpp:902-906).
        for i in begin_usize..end_usize_exclusive {
            self.all_read_assignments[i] = Vec::new();
        }

        ret
    }

    /// Ported from `Genotyper::BuildAlleleEquivalentClass`
    /// (`Genotyper.hpp:1072-1139`).
    ///
    /// For each allele, computes a fingerprint hash from the sorted-descending
    /// (by fingerprint, then ascending by allele index -- `CompSortPairByBDec`,
    /// `Genotyper.hpp:145-150`, a TOTAL order since it always falls back to
    /// `p1.a < p2.a`, so `std::sort`'s non-stability never matters here)
    /// `readsInAllele` read-index sequence, then scans backward from each
    /// position through same-fingerprint alleles doing a real
    /// [`Genotyper::is_assigned_read_the_same`] check to find (or start) its
    /// equivalence class. Deterministic (integer fingerprint hash + a full
    /// sort with total tie-break + `==` comparisons on already-computed
    /// values) -> byte-identical to the C++ oracle. This EC grouping AND
    /// its index order are load-bearing for
    /// [`Genotyper::quantify_allele_equivalent_class`] and Phase 5c's allele
    /// selection.
    ///
    /// Ends by calling [`Genotyper::remove_low_mapq_allele_in_equivalent_class`]
    /// (`Genotyper.hpp:1136`), matching the C++ unconditionally.
    ///
    /// Requires `self.reads_in_allele` to already be populated (by
    /// [`Genotyper::finalize_read_assignments`], which is this method's only
    /// caller in the stock pipeline).
    ///
    /// Returns the number of equivalence classes built (mirrors the C++
    /// `ecCnt`/return value, `Genotyper.hpp:1097,1138`) -- `0` if
    /// `allele_cnt == 0` or no allele has any assigned reads (mirrors the
    /// C++'s early `return 0` at `Genotyper.hpp:1100-1101`).
    ///
    /// # Panics
    ///
    /// Panics if `self.allele_cnt` is negative, or if `self.reads_in_allele`
    /// has fewer entries than `self.allele_cnt` (both indicate a caller
    /// skipped [`Genotyper::init_allele_info`]/[`Genotyper::init_read_assignments`]).
    pub fn build_allele_equivalent_class(&mut self) -> i32 {
        const FINGERPRINT_MAX: i64 = 1_000_003;

        let allele_cnt_usize =
            usize::try_from(self.allele_cnt).expect("allele_cnt is non-negative");
        let read_cnt_i64 = i64::from(self.read_cnt);

        // `alleleFingerprint`: (alleleIdx, fingerprint), fingerprint == -1
        // means "no assigned reads" (Genotyper.hpp:1079-1091).
        let mut allele_fingerprint: Vec<(i32, i64)> = Vec::with_capacity(allele_cnt_usize);
        for i in 0..allele_cnt_usize {
            self.allele_info[i].equivalent_class = -1;
            let reads = &self.reads_in_allele[i];
            let mut fp: i64 = if reads.is_empty() { -1 } else { 0 };
            if !reads.is_empty() {
                for &(read_idx, _) in reads {
                    // `(uint32_t)np.b * readCnt + readsInAllele[i][j].a`:
                    // the C++ casts the running fingerprint to `uint32_t`
                    // before the multiply (`np.b` is `int`/`_pair::b`), so
                    // replicate that exact 32-bit-unsigned-then-widen
                    // arithmetic rather than doing the multiply in i64
                    // throughout.
                    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
                    let fp_u32 = fp as u32;
                    fp = (i64::from(fp_u32) * read_cnt_i64 + i64::from(read_idx)) % FINGERPRINT_MAX;
                }
            }
            allele_fingerprint.push((i32::try_from(i).expect("allele idx fits in i32"), fp));
        }

        // CompSortPairByBDec: descending by fingerprint (`.b`), ascending by
        // allele idx (`.a`) on ties -- a total order, so `sort_by` (stable)
        // and `sort_unstable_by` agree bit-for-bit here.
        allele_fingerprint
            .sort_by(|p1, p2| if p1.1 == p2.1 { p1.0.cmp(&p2.0) } else { p2.1.cmp(&p1.1) });

        let mut ec_cnt = 0i32;
        self.equivalent_class_to_alleles.clear();

        if self.allele_cnt == 0 || allele_fingerprint[0].1 == -1 {
            return 0;
        }

        for i in 0..allele_cnt_usize {
            let (allele_idx, fp_i) = allele_fingerprint[i];
            if fp_i == -1 {
                break;
            }

            let mut new_ec = true;
            let mut matched_j: Option<usize> = None;
            let mut j = i;
            while j > 0 {
                j -= 1;
                let (allele_j, fp_j) = allele_fingerprint[j];
                if fp_i != fp_j {
                    break;
                }
                if self.is_assigned_read_the_same(
                    &self.reads_in_allele[usize::try_from(allele_idx).unwrap()],
                    &self.reads_in_allele[usize::try_from(allele_j).unwrap()],
                ) {
                    new_ec = false;
                    matched_j = Some(j);
                    break;
                }
            }

            if new_ec {
                self.equivalent_class_to_alleles.push(vec![allele_idx]);
                self.allele_info[usize::try_from(allele_idx).unwrap()].equivalent_class = ec_cnt;
                ec_cnt += 1;
            } else {
                let matched_allele_idx = allele_fingerprint[matched_j.unwrap()].0;
                let ec_idx =
                    self.allele_info[usize::try_from(matched_allele_idx).unwrap()].equivalent_class;
                self.equivalent_class_to_alleles[usize::try_from(ec_idx).unwrap()].push(allele_idx);
                self.allele_info[usize::try_from(allele_idx).unwrap()].equivalent_class = ec_idx;
            }
        }

        self.remove_low_mapq_allele_in_equivalent_class();

        ec_cnt
    }

    /// Ported from `Genotyper::RemoveLowMAPQAlleleInEquivalentClass`
    /// (`Genotyper.hpp:1330-1368`).
    ///
    /// For each allele, sums `qual` across every coalesced read assignment
    /// that hit it. Then, within each equivalence class, keeps only the
    /// member allele(s) whose qual sum equals the class's maximum (an exact
    /// `==` on already-computed `f64` sums -- deterministic, no new
    /// floating-point arithmetic beyond summation in a fixed iteration
    /// order matching the C++'s `readCnt`-then-per-read-assignment loop).
    fn remove_low_mapq_allele_in_equivalent_class(&mut self) {
        let allele_cnt_usize =
            usize::try_from(self.allele_cnt).expect("allele_cnt is non-negative");
        let mut allele_read_qual = vec![0.0f64; allele_cnt_usize];
        let read_cnt_usize = usize::try_from(self.read_cnt).expect("read_cnt is non-negative");
        for i in 0..read_cnt_usize {
            for a in &self.read_assignments[i] {
                allele_read_qual[usize::try_from(a.allele_idx).unwrap()] += f64::from(a.qual);
            }
        }

        for ec in &mut self.equivalent_class_to_alleles {
            let mut max_qual_sum = -1.0f64;
            for &allele_idx in ec.iter() {
                let q = allele_read_qual[usize::try_from(allele_idx).unwrap()];
                if q > max_qual_sum {
                    max_qual_sum = q;
                }
            }
            #[allow(clippy::float_cmp)]
            let kept: Vec<i32> = ec
                .iter()
                .copied()
                .filter(|&allele_idx| {
                    allele_read_qual[usize::try_from(allele_idx).unwrap()] == max_qual_sum
                })
                .collect();
            *ec = kept;
        }
    }

    /// Ported from `Genotyper::FinalizeReadAssignments` (`Genotyper.hpp:912-939`).
    ///
    /// Builds `reads_in_allele` from the coalesced `read_assignments` (one
    /// `(readIdx, idxWithinReadAssignments)` entry per per-read allele hit,
    /// in `read_assignments[i]`/inner-index order -- this determines the
    /// order [`Genotyper::build_allele_equivalent_class`]'s
    /// `is_assigned_read_the_same` compares against, so it must match the
    /// C++'s nested-loop order exactly), calls
    /// [`Genotyper::build_allele_equivalent_class`], then sets each allele's
    /// `missing_coverage` from `missing_coverage[alleleIdx]` (mirrors
    /// `refSet.GetSeqMissingBaseCoverage(i, 0.01)`, out of this port's scope
    /// -- same "caller supplies the reference-derived input" pattern as
    /// [`Genotyper::init_allele_info`]'s `seq_effective_len`/`seq_weight`).
    ///
    /// Returns the count of coalesced read groups with a non-empty
    /// assignment set (mirrors the C++ `ret`/`++ret`, `Genotyper.hpp:914-921`).
    ///
    /// # Panics
    ///
    /// Panics if `missing_coverage.len() != self.allele_info.len()`.
    #[allow(clippy::needless_range_loop)]
    pub fn finalize_read_assignments(&mut self, missing_coverage: &[i32]) -> i32 {
        assert_eq!(
            missing_coverage.len(),
            self.allele_info.len(),
            "finalize_read_assignments: missing_coverage must have one entry per allele"
        );

        let mut ret = 0i32;
        let read_cnt_usize = usize::try_from(self.read_cnt).expect("read_cnt is non-negative");
        for i in 0..read_cnt_usize {
            let assignment_cnt = self.read_assignments[i].len();
            if assignment_cnt > 0 {
                ret += 1;
            }
            for j in 0..assignment_cnt {
                let allele_idx = self.read_assignments[i][j].allele_idx;
                self.reads_in_allele[usize::try_from(allele_idx).unwrap()].push((
                    i32::try_from(i).expect("read idx fits in i32"),
                    i32::try_from(j).unwrap(),
                ));
            }
        }

        self.build_allele_equivalent_class();

        for (info, &mc) in self.allele_info.iter_mut().zip(missing_coverage) {
            info.missing_coverage = mc;
        }

        ret
    }

    /// Ported from `Genotyper::EMupdate` (`Genotyper.hpp:372-421`).
    ///
    /// One E-step/M-step iteration of the abundance EM. `ec_abundance0` is
    /// the previous iteration's per-EC abundance; `ec_abundance1` receives
    /// the updated abundance (write-only, matching the C++ out-parameter);
    /// `ec_read_count` is scratch space (write-only, zeroed at the start of
    /// this call, matching the C++ `memset`).
    ///
    /// # FLOATS
    ///
    /// Pure `+`/`*`/`/` (no `pow`/`log`/`exp`/`sqrt`) over a FIXED iteration
    /// order: read groups `0..rgCnt` (E-step), then equivalent classes
    /// `0..ecCnt` twice (the normalization-accumulation loop, then the
    /// per-EC-abundance loop, M-step) -- both loops match the C++'s
    /// left-to-right, no-reordering summation exactly (no `rayon`/parallel
    /// reduction). This gives close (within a `1e-6` relative tolerance) but
    /// NOT exact-`f64`-bits agreement with the C++ oracle in general: the
    /// vendored oracle is compiled at `-O3`, which may fuse this method's
    /// `a * b + c` accumulations (e.g. `psum += ec_abundance0[ec_idx] *
    /// adjust`, `normalization += ec_read_count[i] / ...`) via FMA
    /// contraction, a single more-precisely-rounded operation that Rust's
    /// `+`/`*` (always two roundings) does not perform implicitly. Verified
    /// by `crates/unum-core/tests/golden_genotyper_em.rs`, which asserts this
    /// tolerance (not bit equality) and documents the root cause in
    /// `assert_close_within_tolerance`'s doc comment.
    ///
    /// Note the two `let adjust = 1.0 / (...); adjust = 1;` C++ lines
    /// (`Genotyper.hpp:389-390,400-401`): the computed `missingCoverage`-
    /// based adjustment is immediately overwritten with the constant `1.0`
    /// (dead code -- the division is still PERFORMED and its NaN/Inf
    /// potential exists in the C++, but the RESULT is discarded before use).
    /// This port skips computing the discarded value entirely (using `1.0`
    /// directly) since it can never affect the final `f64` bit pattern of
    /// `adjust` -- the assignment fully overwrites it, so the discarded
    /// computation's result (even if it were NaN) is provably irrelevant to
    /// any downstream value.
    ///
    /// Returns `diffSum`, the M-step's sum of `|new - old|` abundance deltas
    /// across all ECs (mirrors the C++ return value, `Genotyper.hpp:407-420`).
    ///
    /// # Panics
    ///
    /// Panics if any `read_group_to_allele_ec` entry's `PairId::a` (EC
    /// index) is negative or out of range for `self.equivalent_class_to_alleles`/
    /// `ec_info`, or if the output slices (`ec_abundance1`/`ec_read_count`)
    /// are shorter than `self.equivalent_class_to_alleles.len()`.
    #[must_use]
    #[allow(clippy::needless_range_loop)]
    pub fn em_update(
        &self,
        ec_abundance0: &[f64],
        ec_abundance1: &mut [f64],
        ec_read_count: &mut [f64],
        read_group_to_allele_ec: &[Vec<PairId>],
        read_group_info: &[ReadGroupInfo],
        ec_info: &[EcInfo],
    ) -> f64 {
        let ec_cnt = self.equivalent_class_to_alleles.len();
        let rg_cnt = read_group_to_allele_ec.len();

        // E-step: find the expected number of reads (Genotyper.hpp:377-404).
        ec_read_count[..ec_cnt].fill(0.0);
        for i in 0..rg_cnt {
            let mut psum = 0.0f64;
            for pid in &read_group_to_allele_ec[i] {
                let ec_idx = usize::try_from(pid.a).expect("ec idx is non-negative");
                // `adjust` is unconditionally overwritten to `1.0` in the
                // C++ before use -- see this method's doc comment.
                let adjust = 1.0f64;
                psum += ec_abundance0[ec_idx] * adjust;
            }
            if psum == 0.0 {
                psum = 1.0;
            }
            for pid in &read_group_to_allele_ec[i] {
                let ec_idx = usize::try_from(pid.a).expect("ec idx is non-negative");
                let adjust = 1.0f64;
                ec_read_count[ec_idx] +=
                    read_group_info[i].count * (ec_abundance0[ec_idx] * adjust / psum);
            }
        }

        // M-step: recompute the abundance (Genotyper.hpp:406-420).
        let mut diff_sum = 0.0f64;
        let mut normalization = 0.0f64;
        for i in 0..ec_cnt {
            normalization += ec_read_count[i] / f64::from(ec_info[i].length);
        }

        for i in 0..ec_cnt {
            let tmp = ec_read_count[i] / f64::from(ec_info[i].length) / normalization;
            diff_sum += (tmp - ec_abundance0[i]).abs();
            ec_abundance1[i] = tmp;
        }

        diff_sum
    }

    /// Ported from `Genotyper::SQUAREMalpha` (`Genotyper.hpp:424-437`).
    ///
    /// Computes the SQUAREM acceleration coefficient from three consecutive
    /// EM iterates `t0`/`t1`/`t2` (each length `n`). Returns `-1.0` if the
    /// second-difference sum-of-squares (`sqrSumV`) is exactly `0.0`
    /// (mirrors the C++'s `if (sqrSumV == 0) return -1`).
    ///
    /// # FLOATS
    ///
    /// `+`/`*`/`sqrt`, over a fixed `0..n` accumulation order matching the
    /// C++ exactly. `sqrt` itself is libm's `f64::sqrt`, which is IEEE-754
    /// correctly-rounded (a single, exactly-specified operation with no
    /// platform-dependent multi-step approximation) on both Rust's std and
    /// C++'s libm on the platforms this port targets. However, the
    /// `sqr_sum_r`/`sqr_sum_v` accumulations feeding into it (`(t1[i] -
    /// t0[i]) * (t1[i] - t0[i])`, summed) are exactly the kind of `a * b + c`
    /// sequence the C++ oracle's `-O3` build may fuse via FMA contraction,
    /// which Rust's two-rounding `+`/`*` does not do implicitly -- so, like
    /// [`Genotyper::em_update`], this is close (within the differential
    /// test's `1e-6` relative tolerance) but not guaranteed exact-`f64`-bits
    /// reproducible across the full SQUAREM loop. See
    /// `crates/unum-core/tests/golden_genotyper_em.rs`'s
    /// `assert_close_within_tolerance` doc comment for the full
    /// justification.
    #[must_use]
    pub fn squarem_alpha(t0: &[f64], t1: &[f64], t2: &[f64], n: usize) -> f64 {
        let mut sqr_sum_r = 0.0f64;
        let mut sqr_sum_v = 0.0f64;
        for i in 0..n {
            sqr_sum_r += (t1[i] - t0[i]) * (t1[i] - t0[i]);
            let v = t2[i] - 2.0 * t1[i] + t0[i];
            sqr_sum_v += v * v;
        }
        if sqr_sum_v == 0.0 {
            return -1.0;
        }
        -sqr_sum_r.sqrt() / sqr_sum_v.sqrt()
    }

    /// Ported from `Genotyper::SetAlleleAbundance` (`Genotyper.hpp:957-1014`),
    /// restricted to the `ecReadCount != NULL` branch (the stock pipeline's
    /// only call pattern from
    /// [`Genotyper::quantify_allele_equivalent_class`]; the `NULL` branch
    /// belongs to `InitAlleleInfo`'s sibling `InitAlleleAbundance`, out of
    /// this port's 5b scope per the module docs).
    ///
    /// Sets each allele's `abundance`/`ec_abundance` from its EC's read
    /// count and length (FPK: fragments per kilobase), then recomputes
    /// `gene_abundance`/`major_allele_abundance`/`gene_max_major_allele_abundance`
    /// from scratch.
    ///
    /// # FLOATS
    ///
    /// Pure `+`/`*`/`/` over the C++'s exact iteration order (per-EC then
    /// per-allele-within-EC for the abundance loop; three separate
    /// `0..allele_cnt` passes for the gene/major-allele aggregation, matching
    /// `Genotyper.hpp:997-1013` exactly). This method's own arithmetic has
    /// little FMA-contraction opportunity in isolation, but its
    /// `ec_read_count` input is itself the (possibly FMA-diverged) output of
    /// [`Genotyper::em_update`]'s SQUAREM loop, so the abundances it
    /// produces inherit that same close-but-not-exact-`f64`-bits
    /// relationship to the C++ oracle -- see [`Genotyper::em_update`]'s
    /// FLOATS note.
    ///
    /// # Panics
    ///
    /// Panics if `ec_read_count`/`ec_info` are shorter than
    /// `self.equivalent_class_to_alleles.len()`, or if `self.gene_cnt`/
    /// `self.major_allele_cnt` are negative.
    #[allow(clippy::needless_range_loop)]
    pub fn set_allele_abundance(&mut self, ec_read_count: &[f64], ec_info: &[EcInfo]) {
        let ec_cnt = self.equivalent_class_to_alleles.len();
        let allele_cnt_usize =
            usize::try_from(self.allele_cnt).expect("allele_cnt is non-negative");

        for info in &mut self.allele_info {
            info.abundance = 0.0;
            info.ec_abundance = 0.0;
        }

        for i in 0..ec_cnt {
            let size = self.equivalent_class_to_alleles[i].len();
            let abund = ec_read_count[i];
            let abund = abund / f64::from(ec_info[i].length) * 1000.0; // FPK
            for j in 0..size {
                let k = usize::try_from(self.equivalent_class_to_alleles[i][j]).unwrap();
                #[allow(clippy::cast_precision_loss)]
                let size_f64 = size as f64;
                self.allele_info[k].abundance = abund / size_f64;
                self.allele_info[k].ec_abundance = abund;
            }
        }

        // Set major allele and gene abundances (Genotyper.hpp:989-1013).
        let gene_cnt_usize = usize::try_from(self.gene_cnt).expect("gene_cnt is non-negative");
        let major_allele_cnt_usize =
            usize::try_from(self.major_allele_cnt).expect("major_allele_cnt is non-negative");

        self.gene_abundance = vec![0.0; gene_cnt_usize];
        self.major_allele_abundance = vec![0.0; major_allele_cnt_usize];
        self.gene_max_major_allele_abundance = vec![0.0; gene_cnt_usize];

        for i in 0..allele_cnt_usize {
            let major_idx = usize::try_from(self.allele_info[i].major_allele_idx).unwrap();
            let gene_idx = usize::try_from(self.allele_info[i].gene_idx).unwrap();
            self.major_allele_abundance[major_idx] += self.allele_info[i].abundance;
            self.gene_abundance[gene_idx] += self.allele_info[i].abundance;
        }

        for i in 0..allele_cnt_usize {
            let major_idx = usize::try_from(self.allele_info[i].major_allele_idx).unwrap();
            let gene_idx = usize::try_from(self.allele_info[i].gene_idx).unwrap();
            let abund = self.major_allele_abundance[major_idx];
            if abund > self.gene_max_major_allele_abundance[gene_idx] {
                self.gene_max_major_allele_abundance[gene_idx] = abund;
            }
        }
    }

    /// Ported from `Genotyper::QuantifyAlleleEquivalentClass`
    /// (`Genotyper.hpp:1142-1328`): the SQUAREM-accelerated abundance EM
    /// driver.
    ///
    /// Unlike the C++ (which reads `refSet.GetSeqEffectiveLen`/
    /// `GetSeqWeight` directly), this port takes `seq_effective_len`/
    /// `seq_weight` as caller-supplied parallel slices (one entry per
    /// allele) -- same "caller supplies the reference-derived input"
    /// pattern as [`Genotyper::init_allele_info`].
    ///
    /// Builds `readGroupToAlleleEc` (per coalesced-read-group -> list of
    /// `(ecIdx, qual)`, deduplicated by EC within a group,
    /// `Genotyper.hpp:1165-1189`) and `readGroupInfo` (per-group weight,
    /// the MAX per-allele-assignment weight within the group,
    /// `Genotyper.hpp:1155-1164`), then `ecInfo` (per-EC MIN effective
    /// length / MIN missing coverage across member alleles,
    /// `Genotyper.hpp:1203-1217`), initializes `ecAbundance0` from summed
    /// `seq_weight` per EC (`Genotyper.hpp:1221-1232`), then runs up to
    /// `maxEMIterations` (1000) SQUAREM iterations:
    ///
    /// 1. `EMupdate` (ecAbundance0 -> ecAbundance1)
    /// 2. `EMupdate` (ecAbundance1 -> ecAbundance2)
    /// 3. `SQUAREMalpha(ecAbundance0, ecAbundance1, ecAbundance2)`, clamped
    ///    to `min_squarem_alpha` when `min_squarem_alpha < 0 && alpha <
    ///    min_squarem_alpha` (dead at the stock default `min_squarem_alpha
    ///    == 0`, see [`Genotyper::min_squarem_alpha`]'s doc comment)
    /// 4. Quadratic extrapolation into `ecAbundance3` (the exact formula
    ///    from `Genotyper.hpp:1250-1252`)
    /// 5. A stabilizing `EMupdate` (ecAbundance3 -> ecAbundance1)
    /// 6. `diffSum` convergence check: `< 1e-5` forces one more iteration
    ///    then stops (`Genotyper.hpp:1289-1290`)
    /// 7. Every `maskRound` (10) iterations (but not iteration 0): a
    ///    low-abundance masking pass via [`Genotyper::set_allele_abundance`]
    ///    and the `filter_frac` threshold, then `ecAbundance0` is reset from
    ///    the (possibly just-zeroed) `alleleInfo[...].ecAbundance`
    ///    (`Genotyper.hpp:1292-1313`)
    ///
    /// then a final [`Genotyper::set_allele_abundance`] call
    /// (`Genotyper.hpp:1316`).
    ///
    /// # FLOATS
    ///
    /// See [`Genotyper::em_update`]/[`Genotyper::squarem_alpha`]/
    /// [`Genotyper::set_allele_abundance`]'s own FLOATS notes: every kernel
    /// this driver calls is `+`/`*`/`/`/`sqrt` in a fixed accumulation order
    /// matching the C++ exactly (no parallelism, no reordering), but this
    /// does NOT make the final per-EC/per-allele abundances exact-`f64`-bits
    /// reproducible against the C++ oracle in general -- the vendored oracle
    /// is compiled at `-O3` with FMA contraction, which fuses `a * b + c`
    /// accumulations inside `em_update`/`squarem_alpha` into a single,
    /// differently-rounded operation that Rust does not perform implicitly,
    /// and this per-operation difference compounds over up to 1000 SQUAREM
    /// iterations. This port instead targets (and its differential test,
    /// `crates/unum-core/tests/golden_genotyper_em.rs`, verifies) agreement
    /// within a tight `1e-6` relative tolerance, while the deterministic
    /// structure this method depends on (equivalence classes, read groups)
    /// remains exact-`f64`-bits reproducible since it has no floating point
    /// in its own computation. The genotype CALL (Phase 5c, end-to-end) is
    /// the ultimate correctness gate, not bitwise abundance parity.
    ///
    /// Returns the number of EM iterations run (mirrors the C++ `ret`,
    /// incremented once per outer-loop pass regardless of early
    /// convergence, `Genotyper.hpp:1147,1219,1236`).
    ///
    /// # Panics
    ///
    /// Panics if `seq_effective_len.len() != self.allele_info.len()` or
    /// `seq_weight.len() != self.allele_info.len()`.
    #[allow(clippy::too_many_lines, clippy::needless_range_loop)]
    pub fn quantify_allele_equivalent_class(
        &mut self,
        seq_effective_len: &[i32],
        seq_weight: &[i32],
    ) -> i32 {
        // Max EM iterations (Genotyper.hpp:1195) and the low-abundance
        // masking cadence (Genotyper.hpp:1220) -- hoisted to the top of the
        // function (rather than declared right before their first use) only
        // to satisfy clippy::items_after_statements, matching
        // init_allele_info's own `LARGE_DELETION` hoist.
        const MAX_EM_ITERATIONS: i32 = 1000;
        const MASK_ROUND: i32 = 10;

        assert_eq!(
            seq_effective_len.len(),
            self.allele_info.len(),
            "quantify_allele_equivalent_class: seq_effective_len must have one entry per allele"
        );
        assert_eq!(
            seq_weight.len(),
            self.allele_info.len(),
            "quantify_allele_equivalent_class: seq_weight must have one entry per allele"
        );

        let ec_cnt = self.equivalent_class_to_alleles.len();
        let read_cnt_usize = usize::try_from(self.read_cnt).expect("read_cnt is non-negative");

        // Convert readassignment_to_allele to readassignment_to_alleleEquivalentClass
        // (Genotyper.hpp:1149-1189).
        let mut read_group_to_allele_ec: Vec<Vec<PairId>> = vec![Vec::new(); read_cnt_usize];
        let mut read_group_info: Vec<ReadGroupInfo> =
            vec![ReadGroupInfo { count: 0.0 }; read_cnt_usize];

        for i in 0..read_cnt_usize {
            let size = self.read_assignments[i].len();
            let mut count = f64::from(self.read_assignments[i][0].weight);
            for j in 1..size {
                let w = f64::from(self.read_assignments[i][j].weight);
                if w > count {
                    count = w;
                }
            }
            read_group_info[i].count = count;
        }
        for i in 0..read_cnt_usize {
            let size = self.read_assignments[i].len();
            let mut ec_used: HashMap<i32, usize> = HashMap::new();
            for j in 0..size {
                let allele_idx = self.read_assignments[i][j].allele_idx;
                let ec_idx =
                    self.allele_info[usize::try_from(allele_idx).unwrap()].equivalent_class;
                if let Some(&k) = ec_used.get(&ec_idx) {
                    // Should not happen though, the equivalent class makes
                    // sure the quality score is the same.
                    let qual = f64::from(self.read_assignments[i][j].qual);
                    if qual > read_group_to_allele_ec[i][k].b {
                        read_group_to_allele_ec[i][k].b = qual;
                    }
                } else {
                    let np = PairId { a: ec_idx, b: f64::from(self.read_assignments[i][j].qual) };
                    ec_used.insert(ec_idx, read_group_to_allele_ec[i].len());
                    read_group_to_allele_ec[i].push(np);
                }
            }
        }

        // Start the EM algorithm (Genotyper.hpp:1193-1232).
        let mut ec_abundance0 = vec![0.0f64; ec_cnt];
        let mut ec_abundance1 = vec![0.0f64; ec_cnt];
        let mut ec_abundance2 = vec![0.0f64; ec_cnt];
        let mut ec_abundance3 = vec![0.0f64; ec_cnt];
        let mut ec_read_count = vec![0.0f64; ec_cnt];
        let mut ec_info = vec![EcInfo { length: 0, missing_coverage: 0 }; ec_cnt];

        for i in 0..ec_cnt {
            let members = &self.equivalent_class_to_alleles[i];
            let size = members.len();
            let mut length = seq_effective_len[usize::try_from(members[0]).unwrap()];
            let mut missing_coverage =
                self.allele_info[usize::try_from(members[0]).unwrap()].missing_coverage;
            for j in 1..size {
                let len = seq_effective_len[usize::try_from(members[j]).unwrap()];
                if len < length {
                    length = len;
                }
                let mc = self.allele_info[usize::try_from(members[j]).unwrap()].missing_coverage;
                if mc < missing_coverage {
                    missing_coverage = mc;
                }
            }
            ec_info[i] = EcInfo { length, missing_coverage };
        }

        for i in 0..ec_cnt {
            let mut abund = 0.0f64;
            for &member in &self.equivalent_class_to_alleles[i] {
                abund += f64::from(seq_weight[usize::try_from(member).unwrap()]);
            }
            ec_abundance0[i] = abund;
        }

        // SQUAREM-accelerated EM loop (Genotyper.hpp:1234-1314).
        let mut ret = 0i32;
        let mut t = 0i32;
        while t < MAX_EM_ITERATIONS {
            ret += 1;

            // Both `diffSum` results are discarded here, matching the C++
            // (`Genotyper.hpp:1237-1240`: the first is unused entirely, the
            // second is immediately overwritten by `SQUAREMalpha`'s
            // unrelated `alpha` on the next line) -- `diff_sum` is
            // recomputed manually below (`Genotyper.hpp:1279-1287`) from the
            // POST-extrapolation stabilizing `EMupdate` call instead.
            let _ = self.em_update(
                &ec_abundance0,
                &mut ec_abundance1,
                &mut ec_read_count,
                &read_group_to_allele_ec,
                &read_group_info,
                &ec_info,
            );
            let _ = self.em_update(
                &ec_abundance1,
                &mut ec_abundance2,
                &mut ec_read_count,
                &read_group_to_allele_ec,
                &read_group_info,
                &ec_info,
            );

            let mut alpha =
                Self::squarem_alpha(&ec_abundance0, &ec_abundance1, &ec_abundance2, ec_cnt);
            if self.min_squarem_alpha < 0.0 && alpha < self.min_squarem_alpha {
                alpha = self.min_squarem_alpha;
            }

            for i in 0..ec_cnt {
                ec_abundance3[i] = ec_abundance0[i]
                    - 2.0 * alpha * (ec_abundance1[i] - ec_abundance0[i])
                    + alpha
                        * alpha
                        * (ec_abundance2[i] - 2.0 * ec_abundance1[i] + ec_abundance0[i]);
            }

            let _ = self.em_update(
                &ec_abundance3,
                &mut ec_abundance1,
                &mut ec_read_count,
                &read_group_to_allele_ec,
                &read_group_info,
                &ec_info,
            );

            let mut diff_sum = 0.0f64;
            for i in 0..ec_cnt {
                diff_sum += (ec_abundance1[i] - ec_abundance0[i]).abs();
                ec_abundance0[i] = ec_abundance1[i];
            }

            if diff_sum < 1e-5 && t < MAX_EM_ITERATIONS - 2 {
                t = MAX_EM_ITERATIONS - 2; // Force one more iteration
            }

            if t > 0 && t % MASK_ROUND == 0 {
                // Filter the low abundant ones (Genotyper.hpp:1294-1313).
                self.set_allele_abundance(&ec_read_count, &ec_info);
                for i in 0..usize::try_from(self.allele_cnt).unwrap() {
                    let major_idx = usize::try_from(self.allele_info[i].major_allele_idx).unwrap();
                    let gene_idx = usize::try_from(self.allele_info[i].gene_idx).unwrap();
                    if self.major_allele_abundance[major_idx]
                        < self.filter_frac * 0.5 * self.gene_max_major_allele_abundance[gene_idx]
                    {
                        self.allele_info[i].abundance = 0.0;
                        self.allele_info[i].ec_abundance = 0.0;
                    }
                }

                for i in 0..ec_cnt {
                    let k = usize::try_from(self.equivalent_class_to_alleles[i][0]).unwrap();
                    ec_abundance0[i] = self.allele_info[k].ec_abundance;
                }
            }

            t += 1;
        }

        self.set_allele_abundance(&ec_read_count, &ec_info);

        ret
    }

    /// Ported from `Genotyper::IsReadsInAlleleIdxOptimal` (`Genotyper.hpp:
    /// 198-203`): `true` if the `k`-th `(readIdx, idxWithinReadAssignments)`
    /// entry in `reads_in_allele` refers to a top-quality (`qual == 1.0`)
    /// read assignment.
    #[must_use]
    #[allow(clippy::float_cmp)] // exact C++ `== 1` on qual, not a fuzzy float comparison.
    fn is_reads_in_allele_idx_optimal(&self, reads_in_allele: &[(i32, i32)], k: usize) -> bool {
        let (read_idx, within_idx) = reads_in_allele[k];
        self.read_assignments[usize::try_from(read_idx).unwrap()]
            [usize::try_from(within_idx).unwrap()]
        .qual
            == 1.0
    }

    /// Ported from `Genotyper::GetAverageReadAssignmentCnt`
    /// (`Genotyper.hpp:941-955`): the mean number of allele assignments
    /// per (coalesced) read group, among groups with at least one
    /// assignment. Diagnostic-only (mirrors `Genotyper.cpp:635-636`'s log
    /// line) -- not consumed by [`Genotyper::select_alleles_for_genes`].
    #[must_use]
    pub fn get_average_read_assignment_cnt(&self) -> f64 {
        let mut sum = 0.0f64;
        let mut cnt = 0.0f64;
        for ra in &self.read_assignments {
            if !ra.is_empty() {
                #[allow(clippy::cast_precision_loss)]
                {
                    sum += ra.len() as f64;
                }
                cnt += 1.0;
            }
        }
        sum / cnt
    }

    /// Ported from `Genotyper::GetGeneAlleleTypes` (`Genotyper.hpp:1053-1069`):
    /// the number of distinct allele ranks (0, 1, 2, ...) selected for gene
    /// `gene_idx` so far -- `max(rank) + 1`, or `0` if no allele has been
    /// selected for this gene yet.
    #[must_use]
    pub fn get_gene_allele_types(&self, gene_idx: usize) -> i32 {
        let selected = &self.selected_alleles[gene_idx];
        if selected.is_empty() {
            return 0;
        }
        let mut ret = 0;
        for &(_, rank) in selected {
            if rank > ret {
                ret = rank;
            }
        }
        ret + 1
    }

    /// Ported from `Genotyper::RemoveLowLikelihoodAlleleInEquivalentClass`
    /// (`Genotyper.hpp:1371-1460`): within each equivalence class, keeps
    /// only the member allele(s) whose coverage-range-implied likelihood
    /// (`(effectiveLen / consensusLen) ^ ecAbundance`) is within a `0.05`
    /// relative cutoff of the class's maximum likelihood. `consensus_len_of`
    /// mirrors `refSet.GetSeqConsensusLen` (same "caller supplies
    /// reference-derived input" pattern as [`Genotyper::init_allele_info`]).
    /// Called once (`Genotyper.cpp:647`) between quantification and
    /// [`Genotyper::select_alleles_for_genes`].
    ///
    /// # FLOATS
    ///
    /// Uses `pow` (`f64::powf`) on already-computed `ecAbundance` -- a
    /// transcendental call, so (like [`alnorm`]) this targets close, not
    /// bit-identical, agreement with the C++ oracle; see [`alnorm`]'s FLOATS
    /// note for the same rationale.
    ///
    /// # Panics
    ///
    /// Panics if any equivalence class contains an out-of-range allele
    /// index, or if `self.reads_in_allele`/`self.read_assignments` have not
    /// been populated (i.e. [`Genotyper::finalize_read_assignments`] was not
    /// called first).
    pub fn remove_low_likelihood_allele_in_equivalent_class(
        &mut self,
        consensus_len_of: impl Fn(usize) -> i32,
    ) {
        // Relative-likelihood cutoff (Genotyper.hpp:1441); hoisted to the
        // top of the function only to satisfy clippy::items_after_statements
        // (see read_assignment_to_fragment_assignment's MATCH_CNT_RELAX for
        // the same convention).
        const CUTOFF: f64 = 0.05;

        for ec in &mut self.equivalent_class_to_alleles {
            let size = ec.len();
            let mut min_starts: Vec<i32> =
                ec.iter().map(|&idx| consensus_len_of(usize::try_from(idx).unwrap())).collect();
            let mut max_ends: Vec<i32> = vec![-1; size];
            let mut allele_idx_to_idx: HashMap<i32, usize> = HashMap::new();
            for (j, &allele_idx) in ec.iter().enumerate() {
                allele_idx_to_idx.insert(allele_idx, j);
            }

            // Setting up the range of coverage for each allele.
            let represent_allele_idx = usize::try_from(ec[0]).unwrap();
            let assigned_reads = self.reads_in_allele[represent_allele_idx].clone();
            for &(read_idx, _) in &assigned_reads {
                let read_idx_usize = usize::try_from(read_idx).unwrap();
                for assign in &self.read_assignments[read_idx_usize] {
                    let Some(&idx) = allele_idx_to_idx.get(&assign.allele_idx) else { continue };
                    if assign.start < min_starts[idx] {
                        min_starts[idx] = assign.start;
                    }
                    if assign.end > max_ends[idx] {
                        max_ends[idx] = assign.end;
                    }
                }
            }

            // Compute the likelihood.
            let mut max_likelihood = -1.0f64;
            let mut likelihoods = vec![0.0f64; size];
            for j in 0..size {
                let allele_idx = usize::try_from(ec[j]).unwrap();
                let len = consensus_len_of(allele_idx);
                let mut effective_len = max_ends[j] - min_starts[j] + 1;
                if effective_len > len {
                    effective_len = len;
                }
                let ll = (f64::from(effective_len) / f64::from(len))
                    .powf(self.allele_info[allele_idx].ec_abundance);
                if ll > max_likelihood {
                    max_likelihood = ll;
                }
                likelihoods[j] = ll;
            }

            // Store the kept alleles.
            let mut kept = Vec::new();
            for j in 0..size {
                #[allow(clippy::float_cmp)]
                let is_max = likelihoods[j] == max_likelihood;
                if likelihoods[j] / max_likelihood >= CUTOFF || is_max {
                    kept.push(ec[j]);
                }
            }
            *ec = kept;
        }
    }

    /// Ported from `Genotyper::SelectAllelesForGenes` (`Genotyper.hpp:1462-
    /// 2090`) -- the genotyping decision: ranks equivalence classes by
    /// abundance (descending), assigns each EC's member alleles to their
    /// gene's next available rank (0 = first allele, 1 = second, 2+ =
    /// "secondary" candidates) subject to an abundance-fraction filter (with
    /// a rescue pass for filtered alleles that share a selected allele's
    /// major-allele series), then iteratively searches for the best
    /// rank-0/rank-1 haplotype PAIR per gene (maximizing read coverage, up
    /// to `iterMax = 1000` passes), and finally scores each rank's
    /// `genotypeQuality` via [`alnorm`].
    ///
    /// `seq_weight_of` mirrors `refSet.GetSeqWeight` (`Genotyper.hpp:1926`,
    /// same "caller supplies reference-derived input" pattern used
    /// throughout this module).
    ///
    /// # FLOATS
    ///
    /// Uses [`alnorm`] (a polynomial CDF approximation, itself using `exp`)
    /// and `log`/`sqrt` on already-computed abundances -- transcendental
    /// calls, so (like [`alnorm`] and
    /// [`Genotyper::remove_low_likelihood_allele_in_equivalent_class`]'s
    /// `pow`) this targets close, not bit-identical, `genotypeQuality`
    /// agreement with the C++ oracle. The CALLED alleles themselves
    /// (`selected_alleles`' membership and rank assignment) are decided by
    /// integer/exact-`f64`-comparison logic with no transcendental
    /// functions, so they are expected to match the oracle exactly -- see
    /// `crates/unum-core/tests/golden_genotype_e2e.rs`.
    ///
    /// # Panics
    ///
    /// Panics if any `equivalent_class_to_alleles` entry is empty (would
    /// indicate [`Genotyper::build_allele_equivalent_class`] produced a
    /// malformed EC -- not expected).
    #[allow(clippy::too_many_lines, clippy::needless_range_loop)]
    pub fn select_alleles_for_genes(&mut self, seq_weight_of: impl Fn(usize) -> i32) {
        let gene_cnt_usize = usize::try_from(self.gene_cnt).expect("gene_cnt is non-negative");
        let read_cnt_usize = usize::try_from(self.read_cnt).expect("read_cnt is non-negative");

        let mut read_covered = vec![false; read_cnt_usize];
        self.selected_alleles = vec![Vec::new(); gene_cnt_usize];

        // Compute the abundance for equivalent classes, sorted descending
        // (CompSortPairIntDoubleBDec, Genotyper.hpp:152-156: descending `.b`
        // -- ecAbundance -- ties broken ascending `.a` -- ec index; a total
        // order, so `sort_by`'s stability never matters).
        let ec_cnt = self.equivalent_class_to_alleles.len();
        let mut ec_abundance_list: Vec<(usize, f64)> = (0..ec_cnt)
            .map(|i| {
                (
                    i,
                    self.allele_info
                        [usize::try_from(self.equivalent_class_to_alleles[i][0]).unwrap()]
                    .ec_abundance,
                )
            })
            .collect();
        ec_abundance_list.sort_by(|a, b| {
            // Descending by abundance (`.1`), ties broken ascending by EC
            // index (`.0`) -- CompSortPairIntDoubleBDec (Genotyper.hpp:152-
            // 156). Exact `!=` matches the C++'s own exact `!=` on already-
            // computed abundances, not a fuzzy comparison.
            #[allow(clippy::float_cmp)]
            let abundance_differs = a.1 != b.1;
            if abundance_differs { b.1.partial_cmp(&a.1).unwrap() } else { a.0.cmp(&b.0) }
        });

        let mut filtered_alleles: Vec<i32> = Vec::new();

        for &(ec, _) in &ec_abundance_list {
            let members = self.equivalent_class_to_alleles[ec].clone();
            let allele_idx0 = usize::try_from(members[0]).unwrap();

            if self.allele_info[allele_idx0].ec_abundance <= 1e-6 {
                break;
            }

            // Check whether there are uncovered reads.
            let read_list = self.reads_in_allele[allele_idx0].clone();
            let mut covered = 0.0f64;
            let mut total_assigned_weight = 0.0f64;
            for (j, &(read_idx, _)) in read_list.iter().enumerate() {
                if !self.is_reads_in_allele_idx_optimal(&read_list, j) {
                    continue;
                }
                let weight =
                    f64::from(self.read_assignments[usize::try_from(read_idx).unwrap()][0].weight);
                if read_covered[usize::try_from(read_idx).unwrap()] {
                    covered += weight;
                }
                total_assigned_weight += weight;
            }
            let _ = covered; // matches stock: computed, the "no uncovered reads" early-continue is commented out.

            // Add these alleles to the gene allele.
            let mut genes_to_add: Vec<i32> = Vec::new();
            let mut alleles_to_add: Vec<i32> = Vec::new();
            for &allele_idx in &members {
                let allele_idx_usize = usize::try_from(allele_idx).unwrap();
                let gene_idx =
                    usize::try_from(self.allele_info[allele_idx_usize].gene_idx).unwrap();
                let major_idx =
                    usize::try_from(self.allele_info[allele_idx_usize].major_allele_idx).unwrap();

                let ec_abund = self.allele_info[allele_idx_usize].ec_abundance;
                let gene_max = self.gene_max_major_allele_abundance[gene_idx];
                let major_abund = self.major_allele_abundance[major_idx];

                let mut filter = ec_abund < self.filter_frac * gene_max
                    && (ec_abund * 3.0 >= major_abund
                        || major_abund < 3.0 * self.filter_frac * gene_max);

                if !filter {
                    let selected = &self.selected_alleles[gene_idx];
                    #[allow(clippy::float_cmp)]
                    let covered_equals_total = covered == total_assigned_weight;
                    if covered_equals_total
                        && (ec_abund < 0.25 * gene_max
                            || selected.is_empty()
                            || ec_abund
                                < 0.5
                                    * self.allele_info
                                        [usize::try_from(selected.last().unwrap().0).unwrap()]
                                    .ec_abundance)
                    {
                        filter = true;
                    }
                }

                if filter {
                    filtered_alleles.push(allele_idx);
                    continue;
                }

                let gene_idx_i32 = i32::try_from(gene_idx).unwrap();
                if !genes_to_add.contains(&gene_idx_i32) {
                    genes_to_add.push(gene_idx_i32);
                }
                alleles_to_add.push(allele_idx);
            }

            let quality = if genes_to_add.len() > 1 { 0 } else { 60 };

            if !genes_to_add.is_empty() {
                for &(read_idx, within_idx) in &read_list {
                    let read_idx_usize = usize::try_from(read_idx).unwrap();
                    let within_idx_usize = usize::try_from(within_idx).unwrap();
                    #[allow(clippy::float_cmp)]
                    let qual_is_one =
                        self.read_assignments[read_idx_usize][within_idx_usize].qual == 1.0;
                    if qual_is_one {
                        read_covered[read_idx_usize] = true;
                    }
                }
            }

            let mut gene_allele_types: HashMap<usize, i32> = HashMap::new();
            for &allele_idx in &alleles_to_add {
                let allele_idx_usize = usize::try_from(allele_idx).unwrap();
                let gene_idx =
                    usize::try_from(self.allele_info[allele_idx_usize].gene_idx).unwrap();
                let major_idx = self.allele_info[allele_idx_usize].major_allele_idx;

                let mut allele_rank: i32 = -1;
                for &(a, b) in &self.selected_alleles[gene_idx] {
                    if self.allele_info[usize::try_from(a).unwrap()].major_allele_idx == major_idx {
                        allele_rank = b;
                        break;
                    }
                }
                if allele_rank == -1 {
                    allele_rank = *gene_allele_types
                        .entry(gene_idx)
                        .or_insert_with(|| self.get_gene_allele_types(gene_idx));
                }

                self.allele_info[allele_idx_usize].genotype_quality = quality;
                self.allele_info[allele_idx_usize].allele_rank = allele_rank;

                let ec_abund = self.allele_info[allele_idx_usize].ec_abundance;
                let gene_max = self.gene_max_major_allele_abundance[gene_idx];
                let major_idx_usize = usize::try_from(major_idx).unwrap();
                let major_abund = self.major_allele_abundance[major_idx_usize];
                if ec_abund < self.filter_frac * gene_max
                    && (ec_abund * 3.0 >= major_abund
                        || major_abund < 3.0 * self.filter_frac * gene_max)
                {
                    self.allele_info[allele_idx_usize].genotype_quality = 0;
                }

                self.selected_alleles[gene_idx].push((allele_idx, allele_rank));
            }
        }

        // Rescue some filtered alleles if they are in some valid major
        // allele series (Genotyper.hpp:1669-1695).
        for &allele_idx in &filtered_alleles {
            let allele_idx_usize = usize::try_from(allele_idx).unwrap();
            let gene_idx = usize::try_from(self.allele_info[allele_idx_usize].gene_idx).unwrap();
            if self.selected_alleles[gene_idx].is_empty() {
                continue;
            }
            let major_idx = self.allele_info[allele_idx_usize].major_allele_idx;
            let mut rank: i32 = -1;
            for &(a, b) in &self.selected_alleles[gene_idx] {
                if self.allele_info[usize::try_from(a).unwrap()].major_allele_idx == major_idx {
                    rank = b;
                    break;
                }
            }
            if rank != -1 {
                self.selected_alleles[gene_idx].push((allele_idx, rank));
            }
        }

        self.select_alleles_for_genes_haplotype_search(&mut read_covered, &seq_weight_of);
        self.select_alleles_for_genes_quality_scores();
    }

    /// The `iterMax = 1000` best-haplotype-pair search
    /// (`Genotyper.hpp:1697-1996`), factored out of
    /// [`Genotyper::select_alleles_for_genes`] purely for readability (this
    /// port has no single-huge-function constraint the C++ does). For each
    /// gene with more than 2 selected allele types, searches every rank
    /// pair `(j, k)` (`j <= 1`, `k > j`) for the pair that covers the most
    /// (adjust-weighted) previously-uncovered reads, rearranges the winning
    /// pair to ranks 0/1, and repeats until no gene changes (or `iterMax`
    /// iterations elapse).
    ///
    /// `coveredReads[key] |= 1`/`|= 2` bit tracking in the C++
    /// (`Genotyper.hpp:1830,1851`) is dead weight here: the ONLY use of
    /// those bits is the fully commented-out `it->second & ...` branches
    /// (`Genotyper.hpp:1889-1907`) -- the live code
    /// (`coveredReadCnt += readAssignments[it->first][0].adjustWeight`,
    /// `Genotyper.hpp:1893`) sums over every KEY in the map unconditionally,
    /// never reading `it->second`. This port therefore tracks the covered
    /// read-index SET directly (a plain integer set, via a small sorted
    /// `Vec` -- read counts per gene are small enough that this beats a
    /// `HashSet`'s overhead) rather than reproducing the unread bitmask.
    #[allow(clippy::needless_range_loop, clippy::too_many_lines)]
    fn select_alleles_for_genes_haplotype_search(
        &mut self,
        read_covered: &mut [bool],
        seq_weight_of: &impl Fn(usize) -> i32,
    ) {
        const ITER_MAX: i32 = 1000;
        let gene_cnt_usize = usize::try_from(self.gene_cnt).expect("gene_cnt is non-negative");

        let mut read_coverage = vec![0i32; read_covered.len()];
        let mut used_ec: std::collections::HashSet<i32> = std::collections::HashSet::new();
        for i in 0..gene_cnt_usize {
            used_ec.clear();
            for &(allele_idx, rank) in &self.selected_alleles[i].clone() {
                if rank > 1 {
                    continue;
                }
                let ec = self.allele_info[usize::try_from(allele_idx).unwrap()].equivalent_class;
                if used_ec.contains(&ec) {
                    continue;
                }
                used_ec.insert(ec);
                let reads = self.reads_in_allele[usize::try_from(allele_idx).unwrap()].clone();
                for r in 0..reads.len() {
                    if !self.is_reads_in_allele_idx_optimal(&reads, r) {
                        continue;
                    }
                    let read_idx = usize::try_from(reads[r].0).unwrap();
                    read_coverage[read_idx] += 1;
                }
            }
        }

        // The weight for each missing-coverage value, per gene
        // (Genotyper.hpp:1731-1770).
        let mut missing_coverage_allele_type_weight: Vec<HashMap<i32, f64>> =
            Vec::with_capacity(gene_cnt_usize);
        for i in 0..gene_cnt_usize {
            let allele_type_cnt = self.get_gene_allele_types(i);
            let allele_type_cnt_usize = usize::try_from(allele_type_cnt).unwrap_or(0);
            let mut allele_type_info: Vec<(i32, f64)> = vec![(-1, 0.0); allele_type_cnt_usize];
            for &(allele_idx, allele_type) in &self.selected_alleles[i] {
                let allele_idx_usize = usize::try_from(allele_idx).unwrap();
                let allele_type_usize = usize::try_from(allele_type).unwrap();
                allele_type_info[allele_type_usize].1 +=
                    self.allele_info[allele_idx_usize].abundance;
                let mc = self.allele_info[allele_idx_usize].missing_coverage;
                if allele_type_info[allele_type_usize].0 == -1
                    || mc < allele_type_info[allele_type_usize].0
                {
                    allele_type_info[allele_type_usize].0 = mc;
                }
            }
            let mut weight: HashMap<i32, f64> = HashMap::new();
            for &(mc, abund) in &allele_type_info {
                let entry = weight.entry(mc).or_insert(0.0);
                if abund > *entry {
                    *entry = abund;
                }
            }
            missing_coverage_allele_type_weight.push(weight);
        }

        for _iter in 0..ITER_MAX {
            let mut updated_gene_cnt = 0;
            for i in 0..gene_cnt_usize {
                let allele_type_cnt = self.get_gene_allele_types(i);
                if allele_type_cnt <= 2 {
                    continue;
                }
                let selected = self.selected_alleles[i].clone();
                let mut max_cover = 0.0f64;
                let mut max_cover_abundance = 0.0f64;

                // Remove the effects of the current gene.
                used_ec.clear();
                for &(allele_idx, rank) in &selected {
                    if rank > 1 {
                        continue;
                    }
                    let ec =
                        self.allele_info[usize::try_from(allele_idx).unwrap()].equivalent_class;
                    if used_ec.contains(&ec) {
                        continue;
                    }
                    used_ec.insert(ec);
                    let reads = self.reads_in_allele[usize::try_from(allele_idx).unwrap()].clone();
                    for r in 0..reads.len() {
                        if !self.is_reads_in_allele_idx_optimal(&reads, r) {
                            continue;
                        }
                        let read_idx = usize::try_from(reads[r].0).unwrap();
                        read_coverage[read_idx] -= 1;
                    }
                }

                let mut best_types: Vec<(i32, i32)> = Vec::new();
                let j_upper = (allele_type_cnt - 1).min(2); // `j < alleleTypeCnt - 1 && j <= 1`.
                for j in 0..j_upper {
                    used_ec.clear();
                    let mut covered_from_a: Vec<i32> = Vec::new();
                    let mut allele_j = 0usize;
                    for (l, &(allele_idx, rank)) in selected.iter().enumerate() {
                        if rank != j {
                            continue;
                        }
                        let allele_idx_usize = usize::try_from(allele_idx).unwrap();
                        let ec = self.allele_info[allele_idx_usize].equivalent_class;
                        if used_ec.contains(&ec) {
                            continue;
                        }
                        used_ec.insert(ec);
                        let reads = self.reads_in_allele[allele_idx_usize].clone();
                        for r in 0..reads.len() {
                            let read_idx = usize::try_from(reads[r].0).unwrap();
                            if read_coverage[read_idx] == 0
                                && self.is_reads_in_allele_idx_optimal(&reads, r)
                                && !covered_from_a.contains(&reads[r].0)
                            {
                                covered_from_a.push(reads[r].0);
                            }
                        }
                        allele_j = l;
                    }

                    for k in (j + 1)..allele_type_cnt {
                        let mut covered_reads = covered_from_a.clone();
                        used_ec.clear();
                        let mut allele_k = 0usize;
                        for (l, &(allele_idx, rank)) in selected.iter().enumerate() {
                            if rank != k {
                                continue;
                            }
                            let allele_idx_usize = usize::try_from(allele_idx).unwrap();
                            let ec = self.allele_info[allele_idx_usize].equivalent_class;
                            if used_ec.contains(&ec) {
                                continue;
                            }
                            used_ec.insert(ec);
                            let reads = self.reads_in_allele[allele_idx_usize].clone();
                            for r in 0..reads.len() {
                                let read_idx = usize::try_from(reads[r].0).unwrap();
                                if read_coverage[read_idx] == 0
                                    && self.is_reads_in_allele_idx_optimal(&reads, r)
                                    && !covered_reads.contains(&reads[r].0)
                                {
                                    covered_reads.push(reads[r].0);
                                }
                            }
                            allele_k = l;
                        }

                        let mut abundance_j = 0.0f64;
                        let mut abundance_k = 0.0f64;
                        let mut j_missing_coverage: i32 = -1;
                        let mut k_missing_coverage: i32 = -1;
                        for &(allele_idx, rank) in &selected {
                            let allele_idx_usize = usize::try_from(allele_idx).unwrap();
                            if rank == j {
                                abundance_j += self.allele_info[allele_idx_usize].abundance;
                                let mc = self.allele_info[allele_idx_usize].missing_coverage;
                                if j_missing_coverage == -1 || mc < j_missing_coverage {
                                    j_missing_coverage = mc;
                                }
                            } else if rank == k {
                                abundance_k += self.allele_info[allele_idx_usize].abundance;
                                let mc = self.allele_info[allele_idx_usize].missing_coverage;
                                if k_missing_coverage == -1 || mc < k_missing_coverage {
                                    k_missing_coverage = mc;
                                }
                            }
                        }
                        let abundance_sum = abundance_j * abundance_k;

                        // Summation-order note (Genotyper.hpp:1887-1893): C++
                        // accumulates `coveredReadCnt` by iterating a
                        // `std::map<int, int> coveredReads`, whose iteration
                        // order is always ascending by (integer) key --
                        // i.e. ascending read index -- regardless of insertion
                        // order. `covered_reads` here is a `Vec<i32>` built by
                        // concatenating each qualifying allele's
                        // `reads_in_allele` list (each individually ascending,
                        // since it is built by a single ascending pass over
                        // `read_assignments`) in `selected`-iteration order;
                        // when more than one allele contributes to the same
                        // rank (multiple equivalence classes at rank `j`/`k`),
                        // their read-index ranges are not guaranteed to be
                        // globally interleaved in ascending order, so this
                        // Vec's summation order can differ from C++'s map
                        // order. Because `f64` `+` is not associative, this
                        // CAN produce a sub-ULP difference in `covered_read_cnt`
                        // versus the oracle. This is not reordered to match
                        // C++ exactly (would require a `BTreeMap`/sort here,
                        // a larger change than this fidelity note warrants);
                        // the resulting drift is far below the coverage
                        // margins that drive the `> max_cover`/`== max_cover`
                        // choice below on every fixture observed so far.
                        let mut covered_read_cnt = 0.0f64;
                        for &read_idx in &covered_reads {
                            covered_read_cnt += f64::from(
                                self.read_assignments[usize::try_from(read_idx).unwrap()][0]
                                    .adjust_weight,
                            );
                        }

                        if allele_type_cnt > 3
                            || j_missing_coverage >= 10
                            || k_missing_coverage >= 10
                        {
                            let mut weight_j = missing_coverage_allele_type_weight[i]
                                .get(&j_missing_coverage)
                                .copied()
                                .unwrap_or(0.0);
                            let mut weight_k = missing_coverage_allele_type_weight[i]
                                .get(&k_missing_coverage)
                                .copied()
                                .unwrap_or(0.0);
                            if allele_type_cnt <= 3 {
                                if weight_j >= 1.0 {
                                    weight_j = weight_j.log10();
                                }
                                if weight_k >= 1.0 {
                                    weight_k = weight_k.log10();
                                }
                            }
                            let allele_j_idx = usize::try_from(selected[allele_j].0).unwrap();
                            covered_read_cnt = covered_read_cnt
                                - f64::from(j_missing_coverage)
                                    * weight_j
                                    * f64::from(self.read_length)
                                    / 150.0
                                - f64::from(k_missing_coverage)
                                    * weight_k
                                    * f64::from(self.read_length)
                                    / 150.0
                                + f64::from(seq_weight_of(allele_j_idx));
                        }
                        let _ = allele_k;

                        #[allow(clippy::float_cmp)]
                        let tie = covered_read_cnt == max_cover;
                        if best_types.is_empty()
                            || covered_read_cnt > max_cover
                            || (tie && abundance_sum > max_cover_abundance)
                        {
                            max_cover = covered_read_cnt;
                            max_cover_abundance = abundance_sum;
                            best_types.clear();
                            best_types.push((j, k));
                        } else if tie {
                            best_types.push((j, k));
                        }
                    }
                }

                let best_type = best_types[0];
                if best_type != (0, 1) {
                    updated_gene_cnt += 1;
                    for &(allele_idx, rank) in &selected {
                        let new_rank = if rank == best_type.0 {
                            0
                        } else if rank == best_type.1 {
                            1
                        } else if rank < best_type.0 {
                            rank + 2
                        } else if rank < best_type.1 {
                            rank + 1
                        } else {
                            continue;
                        };
                        let allele_idx_usize = usize::try_from(allele_idx).unwrap();
                        for entry in &mut self.selected_alleles[i] {
                            if entry.0 == allele_idx {
                                entry.1 = new_rank;
                            }
                        }
                        self.allele_info[allele_idx_usize].allele_rank = new_rank;
                    }
                }

                // Update read coverage.
                used_ec.clear();
                for &(allele_idx, rank) in &self.selected_alleles[i].clone() {
                    if rank > 1 {
                        continue;
                    }
                    let allele_idx_usize = usize::try_from(allele_idx).unwrap();
                    let ec = self.allele_info[allele_idx_usize].equivalent_class;
                    if used_ec.contains(&ec) {
                        continue;
                    }
                    used_ec.insert(ec);
                    let reads = self.reads_in_allele[allele_idx_usize].clone();
                    for r in 0..reads.len() {
                        if self.is_reads_in_allele_idx_optimal(&reads, r) {
                            let read_idx = usize::try_from(reads[r].0).unwrap();
                            read_coverage[read_idx] += 1;
                        }
                    }
                }
            }

            if updated_gene_cnt == 0 {
                break;
            }
        }
    }

    /// The `alnorm`-based `genotypeQuality` scoring pass
    /// (`Genotyper.hpp:2010-2088`), factored out of
    /// [`Genotyper::select_alleles_for_genes`] for readability (see
    /// [`Genotyper::select_alleles_for_genes_haplotype_search`]'s doc
    /// comment for the same rationale). For each gene, computes a
    /// null-hypothesis mean abundance per rank (self allele-balance noise
    /// plus `crossGeneRate`-weighted noise from every OTHER gene, scaled by
    /// `geneSimilarity`), then scores each rank via
    /// `-log10(alnorm(2*(sqrt(rankAbund) - sqrt(nullMean)), upper=true))`,
    /// clamped to `[0, 60]` and zeroed if the rank's abundance is below
    /// `filter_cov` -- overwriting each MEMBER allele's `genotype_quality`
    /// (but only when it was `> 0`, i.e. not already zeroed by the
    /// abundance-fraction filter in [`Genotyper::select_alleles_for_genes`]).
    #[allow(clippy::needless_range_loop)]
    fn select_alleles_for_genes_quality_scores(&mut self) {
        const CROSS_ALLELE_RATE: f64 = 0.01;
        let gene_cnt_usize = usize::try_from(self.gene_cnt).expect("gene_cnt is non-negative");

        let mut gene_abundances = vec![0.0f64; gene_cnt_usize];
        for i in 0..gene_cnt_usize {
            for &(allele_idx, _) in &self.selected_alleles[i] {
                gene_abundances[i] +=
                    self.allele_info[usize::try_from(allele_idx).unwrap()].abundance;
            }
        }

        for i in 0..gene_cnt_usize {
            let rank_cnt = usize::try_from(self.get_gene_allele_types(i)).unwrap_or(0);
            let mut allele_rank_abund = vec![0.0f64; rank_cnt];
            for &(allele_idx, rank) in &self.selected_alleles[i] {
                allele_rank_abund[usize::try_from(rank).unwrap()] +=
                    self.allele_info[usize::try_from(allele_idx).unwrap()].abundance;
            }

            let mut cross_gene_noise = 0.0f64;
            for j in 0..gene_cnt_usize {
                if i == j {
                    continue;
                }
                cross_gene_noise +=
                    self.cross_gene_rate * self.gene_similarity[j][i] * gene_abundances[j];
            }

            for j in 0..rank_cnt {
                let null_mean = (gene_abundances[i] - allele_rank_abund[j]) * CROSS_ALLELE_RATE
                    + cross_gene_noise;
                let mut score = 0.0f64;
                if allele_rank_abund[j] != 0.0 {
                    score = -(alnorm(2.0 * (allele_rank_abund[j].sqrt() - null_mean.sqrt()), true))
                        .ln()
                        / 10.0f64.ln();
                }
                // Not `.clamp(0.0, 60.0)`: matches the C++'s two independent
                // `if` checks exactly (`Genotyper.hpp:2071-2074`), which --
                // unlike `f64::clamp` -- do not panic/misbehave if `score`
                // is ever `NaN` (e.g. a degenerate `alnorm`/`ln` input).
                #[allow(clippy::manual_clamp)]
                {
                    if score > 60.0 {
                        score = 60.0;
                    }
                    if score < 0.0 {
                        score = 0.0;
                    }
                }
                if allele_rank_abund[j] < self.filter_cov {
                    score = 0.0;
                }

                let j_i32 = i32::try_from(j).unwrap();
                for &(allele_idx, rank) in &self.selected_alleles[i].clone() {
                    let allele_idx_usize = usize::try_from(allele_idx).unwrap();
                    if rank == j_i32 && self.allele_info[allele_idx_usize].genotype_quality > 0 {
                        #[allow(clippy::cast_possible_truncation)]
                        {
                            self.allele_info[allele_idx_usize].genotype_quality = score as i32;
                        }
                    }
                }
            }
        }
    }

    /// Ported from `Genotyper::IsAlleleSameInExon` (`Genotyper.hpp:133-142`):
    /// `true` if `name_a`/`name_b` share the same `fieldsType = 1`
    /// major-allele string (mirrors [`Genotyper::parse_allele_name`] with
    /// `fields_type = 1`; despite the name, this does NOT compare actual
    /// exon sequence -- it is purely a name-parsing comparison, matching the
    /// C++ exactly).
    #[must_use]
    fn is_allele_same_in_exon(&self, name_a: &str, name_b: &str) -> bool {
        self.parse_allele_name(name_a, 1).1 == self.parse_allele_name(name_b, 1).1
    }

    /// Ported from `Genotyper::GetAlleleDescription` (`Genotyper.hpp:2103-
    /// 2178`): formats gene `gene_idx`'s called alleles into the three
    /// `_genotype.tsv` fields (allele-1, allele-2, secondary-candidates),
    /// one call per gene from the CLI driver's output-writing loop
    /// (`Genotyper.cpp:660-670`). Returns `ret` (the C++'s `calledAlleleCnt`
    /// return value: `0`, `1`, or `2`, counting how many of rank 0/1 got a
    /// non-empty major-allele name).
    ///
    /// # The secondary-candidates field only ever shows the LAST `type >= 2`
    /// group
    ///
    /// The C++ resets `secondaryAlleles[0] = '\0'` unconditionally at the
    /// top of EVERY loop iteration whose `buffer` is `secondaryAlleles`
    /// (i.e. every `type >= 2` -- `Genotyper.hpp:2133`, `buffer[0] = '\0'`,
    /// reached regardless of `type`), so each `type >= 2` iteration
    /// OVERWRITES whatever the previous `type >= 2` iteration wrote. When
    /// `GetGeneAlleleTypes(gene_idx) > 3` (more than one secondary
    /// candidate group), only the group from the HIGHEST `type` value
    /// survives in the final string -- a genuine T1K quirk, not a typo in
    /// this port; preserved exactly.
    ///
    /// # Field formatting
    ///
    /// Each populated field is `{majorAlleleNames}<sep>{abundance:.6}<sep>{qual}`
    /// where `<sep>` is `\t` for the primary two fields (allele1/allele2)
    /// and `;` for the secondary field; multiple major alleles at the same
    /// rank are joined with `,` (same rank via a shared equivalence class)
    /// or `|` (first-vs-later DISTINCT major alleles at the secondary rank
    /// only -- see the C++'s `else // only happens for secondary alleles`
    /// comment, `Genotyper.hpp:2164`). An unpopulated rank-0/1 field
    /// (`localQual < 0`) is instead the literal `.\t0\t-1` (matching
    /// `sprintf(buffer + strlen(buffer), ".\t0\t-1")`, `Genotyper.hpp:2172`
    /// -- note this ignores `sep`, always emitting `\t` even though `sep`
    /// may already be `;` by that point).
    ///
    /// `%lf`'s default C `printf` precision is 6 fractional digits; `{:.6}`
    /// matches that exactly for the finite, non-huge magnitudes real
    /// abundances take (this port does not reproduce `%lf`'s
    /// arbitrary-magnitude/`inf`/`nan` formatting quirks, not reachable for
    /// a real abundance).
    ///
    /// # Panics
    ///
    /// Panics if `self.major_allele_cnt` is negative, or if `gene_idx` is
    /// out of range for `self.selected_alleles`.
    #[must_use]
    pub fn get_allele_description(&self, gene_idx: usize) -> (String, String, String, i32) {
        let major_allele_cnt_usize =
            usize::try_from(self.major_allele_cnt).expect("major_allele_cnt is non-negative");
        let mut used = vec![false; major_allele_cnt_usize];
        let mut ret = 0i32;

        let mut type_cnt = self.get_gene_allele_types(gene_idx);
        if type_cnt < 2 {
            type_cnt = 2;
        }

        let mut allele1 = String::new();
        let mut allele2 = String::new();
        let mut secondary_alleles = String::new();
        let mut qualities: [i32; 2] = [-1, -1];

        for type_ in 0..type_cnt {
            let mut abundance = 0.0f64;
            // `sep` for THIS iteration: '\t' for type <= 1, ';' for type > 1
            // (Genotyper.hpp:2119,2130-2134: `sep` is only ever reassigned
            // to ';' once `type > 1` and never reset back).
            let sep = if type_ > 1 { ';' } else { '\t' };
            let mut added = false;
            let mut buffer = String::new(); // Reset every iteration -- see doc comment for the type>=2 quirk.

            let mut local_qual: i32 = -1;
            if type_ == 1 && qualities[0] == 0 {
                used.fill(false);
            }
            for &(allele_idx, rank) in &self.selected_alleles[gene_idx] {
                if rank != type_ {
                    continue;
                }
                let allele_idx_usize = usize::try_from(allele_idx).unwrap();
                let major_allele_idx =
                    usize::try_from(self.allele_info[allele_idx_usize].major_allele_idx).unwrap();
                abundance += self.allele_info[allele_idx_usize].abundance;
                if !used[major_allele_idx] {
                    local_qual = self.allele_info[allele_idx_usize].genotype_quality;
                    if type_ <= 1 {
                        ret = type_ + 1;
                    }
                    if added {
                        buffer.push(',');
                        buffer.push_str(&self.major_allele_idx_to_name[major_allele_idx]);
                    } else {
                        if buffer.is_empty() {
                            buffer.push_str(&self.major_allele_idx_to_name[major_allele_idx]);
                        } else {
                            buffer.push('|');
                            buffer.push_str(&self.major_allele_idx_to_name[major_allele_idx]);
                        }
                        added = true;
                    }
                    used[major_allele_idx] = true;
                }
            }

            if local_qual >= 0 {
                buffer.push_str(&format!("{sep}{abundance:.6}{sep}{local_qual}"));
            } else if type_ <= 1 {
                buffer.push_str(".\t0\t-1");
            }
            if type_ <= 1 {
                qualities[usize::try_from(type_).unwrap()] = local_qual;
            }

            match type_ {
                0 => allele1 = buffer,
                1 => allele2 = buffer,
                _ => secondary_alleles = buffer,
            }
        }

        (allele1, allele2, secondary_alleles, ret)
    }

    /// Ported from `Genotyper::OutputRepresentativeAlleles` (`Genotyper.hpp:
    /// 2180-2229`): writes `_allele.tsv` -- one `{alleleName} {quality}` line
    /// per representative allele (up to 2 per gene: rank-0 and rank-1,
    /// picking the highest-`ecAbundance` member when a rank has several,
    /// ties broken by the LOWER allele index -- `Genotyper.hpp:2193-2199`),
    /// with a rescue for a homozygous-looking rank-0-only gene: if no rank-1
    /// allele was selected, but ANOTHER rank-0 allele exists that is neither
    /// in the SAME equivalence class as the chosen representative NOR
    /// [`Genotyper::is_allele_same_in_exon`] with it, that allele (the
    /// highest-`ecAbundance` such candidate, ties broken by LOWER allele
    /// index) becomes the rank-1 representative too (`Genotyper.hpp:2201-
    /// 2221`).
    ///
    /// # Panics
    ///
    /// Panics if `path` cannot be created/written.
    pub fn output_representative_alleles(
        &self,
        path: &std::path::Path,
        seq_name_of: impl Fn(usize) -> String,
    ) {
        use std::io::Write as _;

        let gene_cnt_usize = usize::try_from(self.gene_cnt).expect("gene_cnt is non-negative");
        let mut out = std::fs::File::create(path)
            .unwrap_or_else(|e| panic!("creating {}: {e}", path.display()));

        for i in 0..gene_cnt_usize {
            let mut representatives: [Option<i32>; 2] = [None, None];
            for &(allele_idx, tag) in &self.selected_alleles[i] {
                if tag > 1
                    || self.allele_info[usize::try_from(allele_idx).unwrap()].genotype_quality < 1
                {
                    continue;
                }
                let tag_usize = usize::try_from(tag).unwrap();
                let candidate_ec_abund =
                    self.allele_info[usize::try_from(allele_idx).unwrap()].ec_abundance;
                let better = match representatives[tag_usize] {
                    None => true,
                    Some(cur) => {
                        let cur_ec_abund =
                            self.allele_info[usize::try_from(cur).unwrap()].ec_abundance;
                        #[allow(clippy::float_cmp)]
                        let tie = cur_ec_abund == candidate_ec_abund;
                        cur_ec_abund < candidate_ec_abund || (tie && cur > allele_idx)
                    }
                };
                if better {
                    representatives[tag_usize] = Some(allele_idx);
                }
            }

            if representatives[1].is_none() {
                if let Some(rep0) = representatives[0] {
                    let rep0_usize = usize::try_from(rep0).unwrap();
                    let rep0_ec = self.allele_info[rep0_usize].equivalent_class;
                    let rep0_name = seq_name_of(rep0_usize).to_string();
                    let mut max = -1.0f64;
                    let mut max_allele_idx: Option<i32> = None;
                    for &(allele_idx, tag) in &self.selected_alleles[i] {
                        if tag != 0 {
                            continue;
                        }
                        let allele_idx_usize = usize::try_from(allele_idx).unwrap();
                        if self.allele_info[allele_idx_usize].equivalent_class == rep0_ec
                            || self
                                .is_allele_same_in_exon(&seq_name_of(allele_idx_usize), &rep0_name)
                        {
                            continue;
                        }
                        let ec_abund = self.allele_info[allele_idx_usize].ec_abundance;
                        #[allow(clippy::float_cmp)]
                        let tie = ec_abund == max;
                        let better = ec_abund > max
                            || (tie && max_allele_idx.is_some_and(|cur| allele_idx < cur));
                        if better {
                            max = ec_abund;
                            max_allele_idx = Some(allele_idx);
                        }
                    }
                    #[allow(clippy::float_cmp)]
                    // exact C++ `!= -1` sentinel check, not a fuzzy comparison.
                    let found_candidate = max != -1.0;
                    if found_candidate {
                        representatives[1] = max_allele_idx;
                    }
                }
            }

            for rep in representatives.into_iter().flatten() {
                let rep_usize = usize::try_from(rep).unwrap();
                writeln!(
                    out,
                    "{} {}",
                    seq_name_of(rep_usize),
                    self.allele_info[rep_usize].genotype_quality
                )
                .unwrap_or_else(|e| panic!("writing {}: {e}", path.display()));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- parse_allele_name ---

    #[test]
    fn parse_allele_name_kir_default_no_delimiter_config() {
        // alleleDigitUnits=-1, alleleDelimiter='\0', fieldsType=0: no ':' in
        // "KIR2DL1*0010101" -> parseType stays 1 (KIR/no-delimiter),
        // fieldsLength=3 (fieldsType==0 branch).
        let (gene, major) = parse_allele_name("KIR2DL1*0010101", 0, -1, b'\0');
        assert_eq!(gene, "KIR2DL1");
        assert_eq!(major, "KIR2DL1*001");
    }

    #[test]
    fn parse_allele_name_kir_fields_type_1_uses_5_digit_fields_length() {
        let (gene, major) = parse_allele_name("KIR2DL1*0010101", 1, -1, b'\0');
        assert_eq!(gene, "KIR2DL1");
        assert_eq!(major, "KIR2DL1*00101");
    }

    #[test]
    fn parse_allele_name_hla_default_delimiter_detected() {
        // "A*01:01:01:01" contains ':' -> parseType=2 (HLA), fieldsLength=3
        // (fieldsType==0), so majorAllele keeps 3 fields (3 colons -> stop
        // at k>=3).
        let (gene, major) = parse_allele_name("A*01:01:01:01", 0, -1, b'\0');
        assert_eq!(gene, "A");
        assert_eq!(major, "A*01:01:01");
    }

    #[test]
    fn parse_allele_name_hla_fields_type_1_still_3_fields_for_delimiter_names() {
        // parseType==2 (delimiter) at fieldsType>=1 still yields
        // fieldsLength=3 (Genotyper.hpp:89-90).
        let (gene, major) = parse_allele_name("DRB1*03:01:01", 1, -1, b'\0');
        assert_eq!(gene, "DRB1");
        assert_eq!(major, "DRB1*03:01:01");
    }

    #[test]
    fn parse_allele_name_hla_two_field_name_keeps_whole_string() {
        // Only 2 colons present; the field-truncation loop never reaches
        // k>=3, so it runs to the end of the string (majorAllele ==
        // allele).
        let (gene, major) = parse_allele_name("A*01:01", 0, -1, b'\0');
        assert_eq!(gene, "A");
        assert_eq!(major, "A*01:01");
    }

    #[test]
    fn parse_allele_name_explicit_allele_digit_units_overrides_auto_detect() {
        // alleleDigitUnits=2 (not -1) skips the ENTIRE `if (fieldsLength ==
        // -1)` block (Genotyper.hpp:71-92) -- including the ':' delimiter
        // auto-detection loop. With alleleDelimiter also '\0' (not
        // overriding), parseType stays at its 1 (KIR/no-delimiter) default
        // even though "A*01:01:01:01" contains colons: the name is parsed
        // as if it had no delimiter at all, walking `fieldsLength=2` RAW
        // BYTES (not fields) past the '*'. This is a genuinely surprising
        // T1K behavior, not a typo in this port.
        let (gene, major) = parse_allele_name("A*01:01:01:01", 0, 2, b'\0');
        assert_eq!(gene, "A");
        assert_eq!(major, "A*01");
    }

    #[test]
    fn parse_allele_name_explicit_delimiter_forces_parse_type_2() {
        // A KIR-style (no ':') name, but alleleDelimiter='*' forces
        // parseType=2 using '*' as the delimiter (fieldsLength=3, the
        // fieldsType==0 default). The scan for delimiters starts at `i`
        // (the FIRST '*', which is also the gene/allele separator), so
        // that first '*' itself counts as delimiter occurrence k=1; the
        // 3rd occurrence (k>=3) is the '*' before the final "01", where
        // the scan breaks -- one field short of consuming the whole name.
        let (gene, major) = parse_allele_name("KIR2DL1*001*01*01", 0, -1, b'*');
        assert_eq!(gene, "KIR2DL1");
        assert_eq!(major, "KIR2DL1*001*01");
    }

    #[test]
    fn parse_allele_name_no_star_treats_whole_string_as_gene() {
        // No '*' at all: `i` runs to the end of the string for both parse
        // types, so gene == the whole allele string.
        let (gene, major) = parse_allele_name("NOSTARNAME", 0, -1, b'\0');
        assert_eq!(gene, "NOSTARNAME");
        assert_eq!(major, "NOSTARNAME");
    }

    // --- Genotyper::init_allele_info ---

    fn effective_lens(n: usize, len: i32) -> Vec<i32> {
        vec![len; n]
    }

    #[test]
    fn init_allele_info_groups_kir_and_hla_alleles_by_gene_and_major_allele() {
        let names: Vec<String> = [
            "KIR2DL1*0010101",
            "KIR2DL1*0010102", // same major allele as above (KIR2DL1*001)
            "KIR2DL1*00201",   // different major allele, same gene
            "A*01:01:01:01",
            "A*01:01:01:02", // same major allele as above (A*01:01:01)
            "A*01:02:01:01", // different major allele, same gene
            "DRB1*03:01:01", // different gene
        ]
        .iter()
        .map(|s| (*s).to_string())
        .collect();
        let consensus: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGA".to_vec(),
            b"TTTTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"GGGGACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"GGGGACGTACGTACGTACGTACGTACGTACGTACGA".to_vec(),
            b"CCCCACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"AAAAACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        ];
        let weight = vec![1; names.len()];
        let mut eff_len = effective_lens(names.len(), 100);

        let mut g = Genotyper::new();
        g.init_allele_info(&names, &consensus, &weight, &mut eff_len, 8);

        assert_eq!(g.allele_cnt, 7);
        assert_eq!(g.gene_cnt, 3); // KIR2DL1, A, DRB1
        assert_eq!(g.gene_idx_to_name, vec!["KIR2DL1", "A", "DRB1"]);

        // 5 distinct major alleles: KIR2DL1*001, KIR2DL1*002, A*01:01:01,
        // A*01:02, DRB1*03:01:01.
        assert_eq!(g.major_allele_cnt, 5);
        assert_eq!(
            g.major_allele_idx_to_name,
            vec!["KIR2DL1*001", "KIR2DL1*002", "A*01:01:01", "A*01:02:01", "DRB1*03:01:01"]
        );

        // Alleles 0 and 1 share a gene and major allele.
        assert_eq!(g.allele_info[0].gene_idx, g.allele_info[1].gene_idx);
        assert_eq!(g.allele_info[0].major_allele_idx, g.allele_info[1].major_allele_idx);
        // Allele 2 shares the gene but not the major allele.
        assert_eq!(g.allele_info[2].gene_idx, g.allele_info[0].gene_idx);
        assert_ne!(g.allele_info[2].major_allele_idx, g.allele_info[0].major_allele_idx);
        // Alleles 3/4/5 (the A* alleles) share a gene, distinct from
        // KIR2DL1 and DRB1.
        assert_eq!(g.allele_info[3].gene_idx, g.allele_info[4].gene_idx);
        assert_eq!(g.allele_info[3].major_allele_idx, g.allele_info[4].major_allele_idx);
        assert_ne!(g.allele_info[5].major_allele_idx, g.allele_info[3].major_allele_idx);
        assert_eq!(g.allele_info[5].gene_idx, g.allele_info[3].gene_idx);
        // Allele 6 (DRB1) is its own gene.
        assert!(g.allele_info[..6].iter().all(|a| a.gene_idx != g.allele_info[6].gene_idx));

        // majorAlleleSize sums seq_weight per major allele: KIR2DL1*001 has
        // 2 alleles of weight 1 each.
        let kir001_idx = usize::try_from(g.major_allele_name_to_idx["KIR2DL1*001"])
            .expect("major allele idx is non-negative");
        assert_eq!(g.major_allele_size[kir001_idx], 2);

        // Gene similarity matrix is geneCnt x geneCnt with diagonal == 1.0
        // (an exact value, not a fuzzy float -- the C++ literally assigns
        // `geneSimilarity[i][i] = 1.0`, Genotyper.hpp:632).
        assert_eq!(g.gene_similarity.len(), 3);
        for row in &g.gene_similarity {
            assert_eq!(row.len(), 3);
        }
        for i in 0..3 {
            #[allow(clippy::float_cmp)]
            let diagonal_is_one = g.gene_similarity[i][i] == 1.0;
            assert!(diagonal_is_one);
        }

        // Every allele starts with the default bookkeeping fields.
        for info in &g.allele_info {
            assert_eq!(info.allele_rank, -1);
            assert_eq!(info.genotype_quality, -1);
            #[allow(clippy::float_cmp)]
            let abundance_is_zero = info.abundance == 0.0;
            assert!(abundance_is_zero);
            assert!(info.whitelist);
        }
    }

    #[test]
    fn init_allele_info_adjusts_large_deletion_outlier_effective_length() {
        // Three alleles of the same gene: two at effectiveLen=1000 (the
        // mode), one at 400 (a >500 deletion relative to the mode) -- the
        // outlier should be bumped up to 1000.
        let names: Vec<String> =
            ["A*01:01:01", "A*01:01:02", "A*01:02:01"].iter().map(|s| (*s).to_string()).collect();
        let consensus: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGA".to_vec(),
            b"TTTTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        ];
        let weight = vec![1; 3];
        let mut eff_len = vec![1000, 1000, 400];

        let mut g = Genotyper::new();
        g.init_allele_info(&names, &consensus, &weight, &mut eff_len, 8);

        assert_eq!(eff_len, vec![1000, 1000, 1000]);
    }

    #[test]
    fn init_allele_info_does_not_adjust_small_deletion() {
        // The same setup, but the outlier is only 300 shorter than the mode
        // (< largeDeletion=500) -- must NOT be adjusted.
        let names: Vec<String> =
            ["A*01:01:01", "A*01:01:02", "A*01:02:01"].iter().map(|s| (*s).to_string()).collect();
        let consensus: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGA".to_vec(),
            b"TTTTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        ];
        let weight = vec![1; 3];
        let mut eff_len = vec![1000, 1000, 700];

        let mut g = Genotyper::new();
        g.init_allele_info(&names, &consensus, &weight, &mut eff_len, 8);

        assert_eq!(eff_len, vec![1000, 1000, 700]);
    }

    // --- Genotyper::read_assignment_weight ---

    fn frag(similarity: f64, has_n: bool) -> FragmentOverlap {
        FragmentOverlap {
            seq_idx: 0,
            seq_start: 0,
            seq_end: 10,
            match_cnt: 10,
            relaxed_match_cnt: 10,
            similarity,
            has_mate_pair: true,
            o1_from_r2: false,
            qual: 1.0,
            has_n,
        }
    }

    // The `read_assignment_weight_*` tests below assert exact `f64` values
    // (0.01/0.1/0.5/1.0/0.05/etc): `ReadAssignmentWeight` is a pure
    // comparison-and-division ratio with no `pow`/`log` (see
    // `Genotyper::read_assignment_weight`'s doc comment), so its outputs
    // for these hand-picked inputs are exactly representable/reproducible
    // -- this is the same "exact deterministic value, not a fuzzy float"
    // rationale `kmer_count.rs`'s own `#[allow(clippy::float_cmp)]` tests
    // use.
    #[test]
    #[allow(clippy::float_cmp)]
    fn read_assignment_weight_perfect_similarity_is_one() {
        let w = Genotyper::read_assignment_weight(&frag(1.0, false), 0.8);
        assert_eq!(w, 1.0);
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn read_assignment_weight_default_ref_seq_similarity_segment_thresholds() {
        // ref_seq_similarity=0.8 -> segment = (1-0.8)/4 = 0.05 (up to fp
        // rounding). Thresholds (each an exclusive upper bound, `similarity
        // < threshold`): 1-3*segment≈0.85, 1-2*segment=0.90, 1-segment≈0.95.
        assert_eq!(Genotyper::read_assignment_weight(&frag(0.84, false), 0.8), 0.01);
        assert_eq!(Genotyper::read_assignment_weight(&frag(0.85, false), 0.8), 0.01);
        assert_eq!(Genotyper::read_assignment_weight(&frag(0.89, false), 0.8), 0.1);
        // similarity==0.90 is NOT < 0.90, so it falls through to the next
        // (< 0.95) branch, not the 0.1 branch.
        assert_eq!(Genotyper::read_assignment_weight(&frag(0.90, false), 0.8), 0.5);
        assert_eq!(Genotyper::read_assignment_weight(&frag(0.94, false), 0.8), 0.5);
        assert_eq!(Genotyper::read_assignment_weight(&frag(0.95, false), 0.8), 1.0);
        assert_eq!(Genotyper::read_assignment_weight(&frag(0.96, false), 0.8), 1.0);
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn read_assignment_weight_segment_floor_at_point_zero_one() {
        // ref_seq_similarity=0.99 -> (1-0.99)/4 = 0.0025, floored to 0.01.
        // Thresholds become the same as the 0.8 case: 0.97, 0.98, 0.99.
        assert_eq!(Genotyper::read_assignment_weight(&frag(0.96, false), 0.99), 0.01);
        assert_eq!(Genotyper::read_assignment_weight(&frag(0.975, false), 0.99), 0.1);
        assert_eq!(Genotyper::read_assignment_weight(&frag(0.985, false), 0.99), 0.5);
        assert_eq!(Genotyper::read_assignment_weight(&frag(1.0, false), 0.99), 1.0);
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn read_assignment_weight_has_n_divides_by_ten() {
        assert_eq!(Genotyper::read_assignment_weight(&frag(1.0, true), 0.8), 0.1);
        assert_eq!(Genotyper::read_assignment_weight(&frag(0.94, true), 0.8), 0.05);
    }

    // --- Genotyper::init_read_assignments / set_read_assignments ---

    fn small_genotyper() -> Genotyper {
        let names: Vec<String> =
            ["A*01:01:01", "A*01:02:01", "B*07:02:01"].iter().map(|s| (*s).to_string()).collect();
        let consensus: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"TTTTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"CCCCACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        ];
        let weight = vec![1; 3];
        let mut eff_len = vec![100, 100, 100];
        let mut g = Genotyper::new();
        g.init_allele_info(&names, &consensus, &weight, &mut eff_len, 8);
        g
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn set_read_assignments_stores_weighted_assignments() {
        let mut g = small_genotyper();
        g.init_read_assignments(2, 2000);

        let assignment = vec![
            FragmentOverlap {
                seq_idx: 0,
                seq_start: 0,
                seq_end: 10,
                match_cnt: 10,
                relaxed_match_cnt: 10,
                similarity: 1.0,
                has_mate_pair: true,
                o1_from_r2: false,
                qual: 1.0,
                has_n: false,
            },
            FragmentOverlap {
                seq_idx: 1,
                seq_start: 0,
                seq_end: 10,
                match_cnt: 9,
                relaxed_match_cnt: 9,
                similarity: 0.94,
                has_mate_pair: true,
                o1_from_r2: false,
                qual: 0.5,
                has_n: false,
            },
        ];
        g.set_read_assignments(0, &assignment, 0.8, |_, _, _| false);

        let stored = &g.all_read_assignments[0];
        assert_eq!(stored.len(), 2);
        assert_eq!(stored[0].allele_idx, 0);
        assert_eq!(stored[0].weight, 1.0);
        assert_eq!(stored[1].allele_idx, 1);
        assert_eq!(stored[1].weight, 0.5);
        // maxSimilarity == 1.0 (from seq_idx=0), so adjustFactor stays 1.0.
        assert_eq!(stored[0].adjust_weight, 1.0);
        assert_eq!(stored[1].adjust_weight, 0.5);

        // Unset read slot stays empty.
        assert!(g.all_read_assignments[1].is_empty());
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn set_read_assignments_applies_adjust_factor_when_max_similarity_below_one() {
        let mut g = small_genotyper();
        g.init_read_assignments(1, 2000);

        let assignment = vec![FragmentOverlap {
            seq_idx: 0,
            seq_start: 0,
            seq_end: 10,
            match_cnt: 9,
            relaxed_match_cnt: 9,
            similarity: 0.94,
            has_mate_pair: true,
            o1_from_r2: false,
            qual: 0.5,
            has_n: false,
        }];
        g.set_read_assignments(0, &assignment, 0.8, |_, _, _| false);

        let stored = &g.all_read_assignments[0];
        assert_eq!(stored.len(), 1);
        assert_eq!(stored[0].weight, 0.5);
        // maxSimilarity=0.94 < 1.0 -> adjustFactor=0.25 -> adjustWeight = 0.25*0.5.
        assert_eq!(stored[0].adjust_weight, 0.125);
    }

    #[test]
    fn set_read_assignments_skips_non_whitelisted_alleles() {
        let mut g = small_genotyper();
        g.allele_info[1].whitelist = false;
        g.init_read_assignments(1, 2000);

        let assignment = vec![
            FragmentOverlap {
                seq_idx: 0,
                seq_start: 0,
                seq_end: 10,
                match_cnt: 10,
                relaxed_match_cnt: 10,
                similarity: 1.0,
                has_mate_pair: true,
                o1_from_r2: false,
                qual: 1.0,
                has_n: false,
            },
            FragmentOverlap {
                seq_idx: 1,
                seq_start: 0,
                seq_end: 10,
                match_cnt: 10,
                relaxed_match_cnt: 10,
                similarity: 1.0,
                has_mate_pair: true,
                o1_from_r2: false,
                qual: 1.0,
                has_n: false,
            },
        ];
        g.set_read_assignments(0, &assignment, 0.8, |_, _, _| false);

        let stored = &g.all_read_assignments[0];
        assert_eq!(stored.len(), 1);
        assert_eq!(stored[0].allele_idx, 0);
    }

    #[test]
    fn set_read_assignments_returns_empty_when_exceeding_max_assign_cnt() {
        let mut g = small_genotyper();
        g.init_read_assignments(1, 1); // maxAssignCnt=1

        let assignment = vec![
            FragmentOverlap {
                seq_idx: 0,
                seq_start: 0,
                seq_end: 10,
                match_cnt: 10,
                relaxed_match_cnt: 10,
                similarity: 1.0,
                has_mate_pair: true,
                o1_from_r2: false,
                qual: 1.0,
                has_n: false,
            },
            FragmentOverlap {
                seq_idx: 1,
                seq_start: 0,
                seq_end: 10,
                match_cnt: 10,
                relaxed_match_cnt: 10,
                similarity: 1.0,
                has_mate_pair: true,
                o1_from_r2: false,
                qual: 1.0,
                has_n: false,
            },
        ];
        g.set_read_assignments(0, &assignment, 0.8, |_, _, _| false);

        assert!(g.all_read_assignments[0].is_empty());
    }

    #[test]
    fn set_read_assignments_returns_empty_when_spanning_separator() {
        let mut g = small_genotyper();
        g.init_read_assignments(1, 2000);

        let assignment = vec![FragmentOverlap {
            seq_idx: 0,
            seq_start: 5,
            seq_end: 15,
            match_cnt: 10,
            relaxed_match_cnt: 10,
            similarity: 1.0,
            has_mate_pair: true,
            o1_from_r2: false,
            qual: 1.0,
            has_n: false,
        }];
        // separator_lookup reports a separator inside [5, 15] for seq_idx 0.
        g.set_read_assignments(0, &assignment, 0.8, |seq_idx, s, e| {
            seq_idx == 0 && s <= 10 && 10 <= e
        });

        assert!(g.all_read_assignments[0].is_empty());
    }

    #[test]
    fn init_read_assignments_resets_state() {
        let mut g = small_genotyper();
        g.init_read_assignments(3, 500);
        assert_eq!(g.read_cnt, 0);
        assert_eq!(g.total_read_cnt, 3);
        assert_eq!(g.max_assign_cnt, 500);
        assert_eq!(g.all_read_assignments.len(), 3);
        assert!(g.all_read_assignments.iter().all(Vec::is_empty));
    }

    // --- Phase 5b: coalesce_read_assignments / build_allele_equivalent_class
    // / finalize_read_assignments / EM ---

    /// Sets `g.all_read_assignments[read_id]` directly (bypassing
    /// `set_read_assignments`'s weight computation/whitelist filtering) so
    /// tests can hand-construct exact `ReadAssignment` sets.
    fn set_raw_assignment(g: &mut Genotyper, read_id: usize, assignments: Vec<ReadAssignment>) {
        g.all_read_assignments[read_id] = assignments;
    }

    fn ra(allele_idx: i32, weight: f32, qual: f32) -> ReadAssignment {
        ReadAssignment { allele_idx, start: 0, end: 10, weight, qual, adjust_weight: weight }
    }

    #[test]
    fn coalesce_read_assignments_merges_identical_fingerprint_and_reassignment() {
        // Two reads that hit alleles {0, 1} with the same qual should
        // coalesce into ONE read group, with weight/adjustWeight summed.
        let mut g = small_genotyper();
        g.init_read_assignments(3, 2000);
        set_raw_assignment(&mut g, 0, vec![ra(0, 1.0, 1.0), ra(1, 0.5, 1.0)]);
        set_raw_assignment(&mut g, 1, vec![ra(0, 1.0, 1.0), ra(1, 0.5, 1.0)]);
        // A third read that hits only allele 2 -- a distinct fingerprint,
        // must become its own read group.
        set_raw_assignment(&mut g, 2, vec![ra(2, 1.0, 1.0)]);

        let ret = g.coalesce_read_assignments(0, 2);
        assert_eq!(ret, 3, "all three reads had a non-empty assignment set");
        assert_eq!(g.read_cnt, 2, "reads 0 and 1 coalesce; read 2 is separate");

        // The coalesced group's weight is the SUM of both reads' weights.
        let coalesced = g
            .read_assignments
            .iter()
            .find(|ra_vec| ra_vec.iter().any(|a| a.allele_idx == 0))
            .unwrap();
        assert_eq!(coalesced.len(), 2);
        #[allow(clippy::float_cmp)]
        let weight_summed = coalesced[0].weight == 2.0;
        assert!(weight_summed, "weight should be summed across the two coalesced reads");

        // allReadAssignments is freed after coalescing.
        assert!(g.all_read_assignments[0].is_empty());
        assert!(g.all_read_assignments[1].is_empty());
        assert!(g.all_read_assignments[2].is_empty());
    }

    #[test]
    fn coalesce_read_assignments_keeps_different_qual_as_separate_groups() {
        // Same allele set {0, 1} but DIFFERENT qual -> NOT the same read
        // assignment (IsReadAssignmentTheSame checks qual equality too).
        let mut g = small_genotyper();
        g.init_read_assignments(2, 2000);
        set_raw_assignment(&mut g, 0, vec![ra(0, 1.0, 1.0), ra(1, 0.5, 1.0)]);
        set_raw_assignment(&mut g, 1, vec![ra(0, 1.0, 0.5), ra(1, 0.5, 0.5)]);

        g.coalesce_read_assignments(0, 1);
        assert_eq!(g.read_cnt, 2, "differing qual must keep the two reads in separate groups");
    }

    #[test]
    fn coalesce_read_assignments_empty_assignment_is_skipped() {
        let mut g = small_genotyper();
        g.init_read_assignments(2, 2000);
        set_raw_assignment(&mut g, 0, vec![]);
        set_raw_assignment(&mut g, 1, vec![ra(0, 1.0, 1.0)]);

        let ret = g.coalesce_read_assignments(0, 1);
        assert_eq!(ret, 1, "only read 1 had a non-empty assignment set");
        assert_eq!(g.read_cnt, 1);
    }

    #[test]
    fn build_allele_equivalent_class_groups_identical_read_sets() {
        // Alleles 0 and 1 are hit by the exact same coalesced read group ->
        // same equivalence class. Allele 2 is hit by a different read group
        // -> its own class.
        let mut g = small_genotyper();
        g.init_read_assignments(2, 2000);
        set_raw_assignment(&mut g, 0, vec![ra(0, 1.0, 1.0), ra(1, 1.0, 1.0)]);
        set_raw_assignment(&mut g, 1, vec![ra(2, 1.0, 1.0)]);
        g.coalesce_read_assignments(0, 1);

        let missing_coverage = vec![0; g.allele_info.len()];
        g.finalize_read_assignments(&missing_coverage);

        assert_eq!(g.allele_info[0].equivalent_class, g.allele_info[1].equivalent_class);
        assert_ne!(g.allele_info[0].equivalent_class, g.allele_info[2].equivalent_class);
        assert_eq!(g.equivalent_class_to_alleles.len(), 2);
    }

    #[test]
    fn build_allele_equivalent_class_no_reads_returns_zero() {
        let mut g = small_genotyper();
        g.init_read_assignments(1, 2000);
        // No assignments at all set for the one read.
        g.coalesce_read_assignments(0, 0);
        let missing_coverage = vec![0; g.allele_info.len()];
        g.finalize_read_assignments(&missing_coverage);

        assert_eq!(g.equivalent_class_to_alleles.len(), 0);
        assert!(g.allele_info.iter().all(|a| a.equivalent_class == -1));
    }

    #[test]
    fn finalize_read_assignments_sets_missing_coverage_from_input() {
        let mut g = small_genotyper();
        g.init_read_assignments(1, 2000);
        set_raw_assignment(&mut g, 0, vec![ra(0, 1.0, 1.0)]);
        g.coalesce_read_assignments(0, 0);

        let missing_coverage = vec![5, 7, 9];
        g.finalize_read_assignments(&missing_coverage);

        assert_eq!(g.allele_info[0].missing_coverage, 5);
        assert_eq!(g.allele_info[1].missing_coverage, 7);
        assert_eq!(g.allele_info[2].missing_coverage, 9);
    }

    #[test]
    #[allow(clippy::float_cmp)]
    fn squarem_alpha_zero_second_difference_returns_negative_one() {
        // t2 - 2*t1 + t0 == 0 everywhere (a perfectly linear sequence) ->
        // sqrSumV == 0 -> the C++'s explicit -1 sentinel.
        let t0 = [1.0, 2.0];
        let t1 = [2.0, 3.0];
        let t2 = [3.0, 4.0];
        assert_eq!(Genotyper::squarem_alpha(&t0, &t1, &t2, 2), -1.0);
    }

    #[test]
    fn squarem_alpha_nonzero_case_matches_hand_computed_value() {
        let t0 = [0.0];
        let t1 = [1.0];
        let t2 = [1.5];
        // sqrSumR = (1-0)^2 = 1. sqrSumV = (1.5 - 2 + 0)^2 = 0.25.
        // alpha = -sqrt(1)/sqrt(0.25) = -1/0.5 = -2.0.
        let alpha = Genotyper::squarem_alpha(&t0, &t1, &t2, 1);
        assert!((alpha - (-2.0)).abs() < 1e-12, "alpha={alpha}");
    }

    #[test]
    fn em_update_two_ecs_shared_read_group_apportions_by_abundance() {
        // Two ECs (alleles 0 and 1 in separate ECs, both length 100), one
        // read group compatible with both at equal qual. E-step should
        // apportion the read group's count proportionally to
        // ecAbundance0, and the M-step normalizes.
        let g = small_genotyper();
        // Fake up equivalent_class_to_alleles as if BuildAlleleEquivalentClass
        // had produced two singleton classes (needed only for `ecCnt` via
        // equivalent_class_to_alleles.len() inside em_update).
        let mut g = g;
        g.equivalent_class_to_alleles = vec![vec![0], vec![1]];

        let read_group_to_allele_ec = vec![vec![PairId { a: 0, b: 1.0 }, PairId { a: 1, b: 1.0 }]];
        let read_group_info = vec![ReadGroupInfo { count: 10.0 }];
        let ec_info = vec![
            EcInfo { length: 100, missing_coverage: 0 },
            EcInfo { length: 100, missing_coverage: 0 },
        ];

        // Equal starting abundance -> the single read group's count should
        // split evenly (psum = 1+1 = 2, each EC gets count * (1/2) = 5).
        let ec_abundance0 = vec![1.0, 1.0];
        let mut ec_abundance1 = vec![0.0; 2];
        let mut ec_read_count = vec![0.0; 2];
        let _ = g.em_update(
            &ec_abundance0,
            &mut ec_abundance1,
            &mut ec_read_count,
            &read_group_to_allele_ec,
            &read_group_info,
            &ec_info,
        );
        assert!((ec_read_count[0] - 5.0).abs() < 1e-12);
        assert!((ec_read_count[1] - 5.0).abs() < 1e-12);
        // Equal length + equal read count -> equal normalized abundance.
        assert!((ec_abundance1[0] - ec_abundance1[1]).abs() < 1e-12);

        // Skewed starting abundance (EC0 favored 3:1) -> the read group
        // count should split 3:1 too (psum = 3+1 = 4).
        let ec_abundance0_skewed = vec![3.0, 1.0];
        let mut ec_abundance1_skewed = vec![0.0; 2];
        let mut ec_read_count_skewed = vec![0.0; 2];
        let _ = g.em_update(
            &ec_abundance0_skewed,
            &mut ec_abundance1_skewed,
            &mut ec_read_count_skewed,
            &read_group_to_allele_ec,
            &read_group_info,
            &ec_info,
        );
        assert!((ec_read_count_skewed[0] - 7.5).abs() < 1e-12);
        assert!((ec_read_count_skewed[1] - 2.5).abs() < 1e-12);
    }

    #[test]
    fn em_update_empty_read_group_gets_psum_floor_of_one() {
        // A read group with zero eligible ECs never occurs in practice
        // (readGroupToAlleleEc entries are only built from real assignments),
        // but an EC with abundance 0 across the board should hit the
        // `psum == 0 -> psum = 1` floor without panicking or dividing by
        // zero.
        let g = small_genotyper();
        let mut g = g;
        g.equivalent_class_to_alleles = vec![vec![0]];
        let read_group_to_allele_ec = vec![vec![PairId { a: 0, b: 1.0 }]];
        let read_group_info = vec![ReadGroupInfo { count: 4.0 }];
        let ec_info = vec![EcInfo { length: 100, missing_coverage: 0 }];

        let ec_abundance0 = vec![0.0];
        let mut ec_abundance1 = vec![0.0];
        let mut ec_read_count = vec![0.0];
        let _ = g.em_update(
            &ec_abundance0,
            &mut ec_abundance1,
            &mut ec_read_count,
            &read_group_to_allele_ec,
            &read_group_info,
            &ec_info,
        );
        // psum floors to 1, so ecReadCount[0] = count * (0 * 1 / 1) = 0.
        assert!((ec_read_count[0] - 0.0).abs() < 1e-12);
    }

    #[test]
    #[allow(clippy::similar_names)] // locus_a_gene_idx / locus_b_gene_idx are deliberately named.
    fn set_allele_abundance_computes_fpk_and_gene_aggregates() {
        let mut g = small_genotyper();
        g.equivalent_class_to_alleles = vec![vec![0], vec![1, 2]];
        let ec_read_count = vec![10.0, 20.0];
        let ec_info = vec![
            EcInfo { length: 100, missing_coverage: 0 },
            EcInfo { length: 200, missing_coverage: 0 },
        ];

        g.set_allele_abundance(&ec_read_count, &ec_info);

        // Allele 0: abund = 10/100*1000 = 100 (FPK), size=1 -> 100/1 = 100.
        assert!((g.allele_info[0].abundance - 100.0).abs() < 1e-9);
        assert!((g.allele_info[0].ec_abundance - 100.0).abs() < 1e-9);
        // Alleles 1,2 share EC1: abund = 20/200*1000 = 100, size=2 -> 50 each.
        assert!((g.allele_info[1].abundance - 50.0).abs() < 1e-9);
        assert!((g.allele_info[2].abundance - 50.0).abs() < 1e-9);
        assert!((g.allele_info[1].ec_abundance - 100.0).abs() < 1e-9);

        // Alleles 0 and 1 are A*01:01:01/A*01:02:01 (different major alleles,
        // same gene "A"); allele 2 is B*07:02:01 (different gene).
        assert_eq!(g.gene_abundance.len(), 2);
        let locus_a_gene_idx = usize::try_from(g.allele_info[0].gene_idx).unwrap();
        let locus_b_gene_idx = usize::try_from(g.allele_info[2].gene_idx).unwrap();
        assert!((g.gene_abundance[locus_a_gene_idx] - 150.0).abs() < 1e-9); // 100 + 50
        assert!((g.gene_abundance[locus_b_gene_idx] - 50.0).abs() < 1e-9);
    }

    /// Builds a `Genotyper` with 2 alleles of the SAME gene/major-allele
    /// group (so `filter_frac` masking, which compares
    /// `majorAlleleAbundance` against `geneMaxMajorAlleleAbundance`, is
    /// exercised meaningfully across the whole pipeline).
    fn two_allele_genotyper() -> Genotyper {
        let names: Vec<String> =
            ["A*01:01:01", "A*01:02:01"].iter().map(|s| (*s).to_string()).collect();
        let consensus: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"TTTTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        ];
        let weight = vec![1; 2];
        let mut eff_len = vec![100, 100];
        let mut g = Genotyper::new();
        g.init_allele_info(&names, &consensus, &weight, &mut eff_len, 8);
        g
    }

    #[test]
    fn quantify_allele_equivalent_class_clean_het_splits_evenly() {
        // A clean het: 5 reads uniquely hit allele 0, 5 uniquely hit allele
        // 1, both alleles the same effective length -> the EM should
        // converge to equal abundance for both.
        let mut g = two_allele_genotyper();
        g.init_read_assignments(10, 2000);
        for i in 0..5 {
            set_raw_assignment(&mut g, i, vec![ra(0, 1.0, 1.0)]);
        }
        for i in 5..10 {
            set_raw_assignment(&mut g, i, vec![ra(1, 1.0, 1.0)]);
        }
        g.coalesce_read_assignments(0, 9);
        let missing_coverage = vec![0; g.allele_info.len()];
        g.finalize_read_assignments(&missing_coverage);

        let seq_effective_len = vec![100, 100];
        let seq_weight = vec![1, 1];
        let iters = g.quantify_allele_equivalent_class(&seq_effective_len, &seq_weight);
        assert!(iters > 0);

        let a0 = g.allele_info[0].abundance;
        let a1 = g.allele_info[1].abundance;
        assert!(a0 > 0.0 && a1 > 0.0, "both alleles should retain nonzero abundance: {a0} {a1}");
        assert!((a0 - a1).abs() / a0.max(a1) < 1e-6, "expected near-equal abundance: {a0} vs {a1}");
    }

    #[test]
    fn quantify_allele_equivalent_class_homozygous_all_reads_one_allele() {
        let mut g = two_allele_genotyper();
        g.init_read_assignments(10, 2000);
        for i in 0..10 {
            set_raw_assignment(&mut g, i, vec![ra(0, 1.0, 1.0)]);
        }
        g.coalesce_read_assignments(0, 9);
        let missing_coverage = vec![0; g.allele_info.len()];
        g.finalize_read_assignments(&missing_coverage);

        let seq_effective_len = vec![100, 100];
        let seq_weight = vec![1, 1];
        g.quantify_allele_equivalent_class(&seq_effective_len, &seq_weight);

        assert!(g.allele_info[0].abundance > 0.0);
        #[allow(clippy::float_cmp)]
        let allele1_is_zero = g.allele_info[1].abundance == 0.0;
        assert!(allele1_is_zero, "allele 1 got no reads at all -> zero abundance");
    }

    #[test]
    fn quantify_allele_equivalent_class_shared_read_ambiguity_apportions() {
        // 10 reads unique to allele 0, 2 reads unique to allele 1, and 4
        // reads compatible with BOTH -- the EM should apportion the
        // ambiguous reads and still call allele 0 the dominant one.
        let mut g = two_allele_genotyper();
        g.init_read_assignments(16, 2000);
        for i in 0..10 {
            set_raw_assignment(&mut g, i, vec![ra(0, 1.0, 1.0)]);
        }
        for i in 10..12 {
            set_raw_assignment(&mut g, i, vec![ra(1, 1.0, 1.0)]);
        }
        for i in 12..16 {
            set_raw_assignment(&mut g, i, vec![ra(0, 1.0, 1.0), ra(1, 1.0, 1.0)]);
        }
        g.coalesce_read_assignments(0, 15);
        let missing_coverage = vec![0; g.allele_info.len()];
        g.finalize_read_assignments(&missing_coverage);

        let seq_effective_len = vec![100, 100];
        let seq_weight = vec![1, 1];
        g.quantify_allele_equivalent_class(&seq_effective_len, &seq_weight);

        assert!(
            g.allele_info[0].abundance > g.allele_info[1].abundance,
            "allele 0 (10 unique reads) should out-abund allele 1 (2 unique reads): {} vs {}",
            g.allele_info[0].abundance,
            g.allele_info[1].abundance
        );
    }

    #[test]
    fn quantify_allele_equivalent_class_equivalent_alleles_share_one_ec() {
        // Both alleles hit by the EXACT same coalesced read groups -> one
        // equivalence class, abundance split evenly between them by
        // set_allele_abundance's `abund / size`.
        let mut g = two_allele_genotyper();
        g.init_read_assignments(5, 2000);
        for i in 0..5 {
            set_raw_assignment(&mut g, i, vec![ra(0, 1.0, 1.0), ra(1, 1.0, 1.0)]);
        }
        g.coalesce_read_assignments(0, 4);
        let missing_coverage = vec![0; g.allele_info.len()];
        g.finalize_read_assignments(&missing_coverage);

        assert_eq!(
            g.equivalent_class_to_alleles.len(),
            1,
            "both alleles should collapse to one EC"
        );

        let seq_effective_len = vec![100, 100];
        let seq_weight = vec![1, 1];
        g.quantify_allele_equivalent_class(&seq_effective_len, &seq_weight);

        assert!((g.allele_info[0].abundance - g.allele_info[1].abundance).abs() < 1e-9);
        assert!(g.allele_info[0].abundance > 0.0);
    }

    #[test]
    fn quantify_allele_equivalent_class_multi_gene_mix_isolates_genes() {
        // Two genes (A and B), 2 alleles each. Reads never cross genes.
        // A's alleles should end up with abundance independent of B's.
        let names: Vec<String> = ["A*01:01:01", "A*01:02:01", "B*07:02:01", "B*08:01:01"]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
        let consensus: Vec<Vec<u8>> = vec![
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"TTTTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"CCCCACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
            b"GGGGACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(),
        ];
        let weight = vec![1; 4];
        let mut eff_len = vec![100; 4];
        let mut g = Genotyper::new();
        g.init_allele_info(&names, &consensus, &weight, &mut eff_len, 8);

        g.init_read_assignments(20, 2000);
        for i in 0..5 {
            set_raw_assignment(&mut g, i, vec![ra(0, 1.0, 1.0)]);
        }
        for i in 5..10 {
            set_raw_assignment(&mut g, i, vec![ra(1, 1.0, 1.0)]);
        }
        for i in 10..17 {
            set_raw_assignment(&mut g, i, vec![ra(2, 1.0, 1.0)]);
        }
        for i in 17..20 {
            set_raw_assignment(&mut g, i, vec![ra(3, 1.0, 1.0)]);
        }
        g.coalesce_read_assignments(0, 19);
        let missing_coverage = vec![0; g.allele_info.len()];
        g.finalize_read_assignments(&missing_coverage);

        let seq_effective_len = vec![100; 4];
        let seq_weight = vec![1; 4];
        g.quantify_allele_equivalent_class(&seq_effective_len, &seq_weight);

        // Gene A: alleles 0/1 near-equal (5 reads each).
        assert!(
            (g.allele_info[0].abundance - g.allele_info[1].abundance).abs()
                / g.allele_info[0].abundance.max(g.allele_info[1].abundance)
                < 1e-6
        );
        // Gene B: allele 2 (7 reads) should out-abund allele 3 (3 reads).
        assert!(g.allele_info[2].abundance > g.allele_info[3].abundance);
        // All four alleles retain nonzero abundance -- no gene contaminated
        // another gene's calls.
        assert!(g.allele_info.iter().all(|a| a.abundance > 0.0));
    }

    // --- alnorm ---

    #[test]
    fn alnorm_at_zero_is_one_half() {
        // Symmetric around 0: both upper and lower tails are exactly 0.5.
        assert!((alnorm(0.0, false) - 0.5).abs() < 1e-9);
        assert!((alnorm(0.0, true) - 0.5).abs() < 1e-9);
    }

    #[test]
    fn alnorm_matches_known_standard_normal_quantiles() {
        // Standard textbook values: Phi(1.96) ~= 0.975, Phi(-1.96) ~= 0.025.
        assert!((alnorm(1.96, false) - 0.975).abs() < 1e-3);
        assert!((alnorm(-1.96, false) - 0.025).abs() < 1e-3);
        // Upper tail is the complement of the lower tail.
        assert!((alnorm(1.96, true) - 0.025).abs() < 1e-3);
    }

    #[test]
    #[allow(clippy::float_cmp)] // exact literal-constant branch (`value = 0.0`/`1.0`), not a fuzzy comparison.
    fn alnorm_far_tail_saturates_to_bounds() {
        // z > ltone (7.0) and beyond utzero (18.66): exact 0.0/1.0 branch.
        assert_eq!(alnorm(20.0, true), 0.0);
        assert_eq!(alnorm(20.0, false), 1.0);
        assert_eq!(alnorm(-20.0, false), 0.0);
    }

    // --- parse_exon_comment ---

    #[test]
    fn parse_exon_comment_pairs_up_skipping_leading_number() {
        // Mirrors kir_rna_seq.fa's real header comment shape.
        let exons = parse_exon_comment("8 50 83 84 119 120 419 420", 500);
        assert_eq!(exons, vec![(50, 83), (84, 119), (120, 419)]);
    }

    #[test]
    fn parse_exon_comment_no_digits_falls_back_to_whole_sequence() {
        assert_eq!(parse_exon_comment("", 100), vec![(0, 99)]);
        assert_eq!(parse_exon_comment("no numbers here", 100), vec![(0, 99)]);
    }

    #[test]
    fn parse_exon_comment_single_number_yields_empty_exon_list() {
        // size == 1 (only the discarded leading number): the C++ for-loop
        // `for (i = 1 ; i < 1 ; i += 2)` runs zero iterations -- an EMPTY
        // list, not the whole-sequence fallback (see doc comment).
        assert_eq!(parse_exon_comment("42", 100), Vec::<(i32, i32)>::new());
    }

    // --- AlleleRef / get_seq_missing_base_coverage ---

    #[test]
    fn get_seq_missing_base_coverage_counts_low_coverage_exon_bases() {
        let allele_ref = AlleleRef::new(b"ACGTACGTAC".to_vec(), None); // whole-seq exon [0,9]
        // Give most positions high coverage, one position near-zero.
        for pw in &allele_ref.pos_weight {
            pw.add(0, 10); // every position "looks like" high-A coverage for this synthetic test
        }
        // one position with a T base but 0 T-count: leave pos_weight[3] at its
        // all-zero default (do NOT add the synthetic count[0]=10 there).
        {
            use std::sync::atomic::Ordering;
            for c in &allele_ref.pos_weight[3].count {
                c.store(0, Ordering::Relaxed);
            }
        }
        // Position 3 is 'T' (0-indexed "ACGTACGTAC"), whose own base-count
        // is 0 (not the synthetic count[0]=10 written above) -- so it reads
        // as a true low-coverage exon base.
        let missing = get_seq_missing_base_coverage(&allele_ref, 0.01);
        assert!(missing >= 1, "expected at least one low-coverage base, got {missing}");
    }

    #[test]
    fn allele_ref_is_exon_respects_parsed_intervals() {
        let allele_ref = AlleleRef::new(vec![b'A'; 20], Some("0 5 9 12 19"));
        // exons: (5,9), (12,19) (leading 0 discarded).
        assert!(!allele_ref.is_exon(0));
        assert!(allele_ref.is_exon(5));
        assert!(allele_ref.is_exon(9));
        assert!(!allele_ref.is_exon(10));
        assert!(allele_ref.is_exon(19));
    }

    // --- extend_overlap ---

    fn simple_overlap(
        seq_idx: u32,
        read_start: i32,
        read_end: i32,
        seq_start: i32,
        seq_end: i32,
        strand: i8,
    ) -> overlap::Overlap {
        overlap::Overlap {
            seq_idx,
            read_start,
            read_end,
            seq_start,
            seq_end,
            strand,
            match_cnt: 2 * (read_end - read_start + 1),
            similarity: 1.0,
        }
    }

    #[test]
    fn extend_overlap_perfect_match_extends_to_full_read() {
        // Allele consensus == read exactly; the k-mer-chained core overlap
        // already spans the read, so left/right overhang is 0 and the
        // extended overlap is unchanged, still perfect similarity.
        let consensus = b"ACGTACGTACGT";
        let read = b"ACGTACGTACGT";
        let core = simple_overlap(0, 0, 11, 0, 11, 1);
        let extended =
            extend_overlap(read, consensus, &core, 0.8, &mut crate::align_algo::DpCache::new())
                .expect("should pass ref_seq_similarity");
        assert_eq!(extended.read_start, 0);
        assert_eq!(extended.read_end, 11);
        #[allow(clippy::float_cmp)]
        let is_one = extended.similarity == 1.0;
        assert!(is_one);
    }

    #[test]
    fn extend_overlap_extends_partial_core_overlap_using_overhangs() {
        // Core overlap only covers the middle of both read and allele; the
        // full read/allele are identical, so extension should walk out to
        // the full span with perfect similarity.
        let consensus = b"AAAACGTACGTAAAA";
        let read = b"AAAACGTACGTAAAA";
        let core = simple_overlap(0, 4, 10, 4, 10, 1);
        let extended =
            extend_overlap(read, consensus, &core, 0.8, &mut crate::align_algo::DpCache::new())
                .unwrap();
        assert_eq!(extended.read_start, 0);
        assert_eq!(extended.read_end, 14);
        assert_eq!(extended.seq_start, 0);
        assert_eq!(extended.seq_end, 14);
    }

    #[test]
    fn extend_overlap_below_similarity_threshold_returns_none() {
        // Core overlap's own matchCnt reflects a fully-mismatched middle
        // (as a real LIS-chained core overlap would compute) -- extension
        // only walks the OUTER overhangs (a perfect match here), so the
        // combined similarity should still fall below a demanding
        // threshold thanks to the core's low matchCnt.
        let consensus = b"AAAATTTTTTTTAAAA";
        let read = b"AAAACCCCCCCCAAAA";
        let mut core = simple_overlap(0, 4, 11, 4, 11, 1); // the mismatching middle only
        core.match_cnt = 0; // no real matches in the core region
        let extended =
            extend_overlap(read, consensus, &core, 0.99, &mut crate::align_algo::DpCache::new());
        assert!(extended.is_none() || extended.unwrap().similarity < 0.99);
    }

    // --- assign_read ---

    /// Snapshots an [`AlleleRef`]'s `pos_weight` into plain `[i32; 4]` arrays
    /// so tests can compare coverage totals with `assert_eq!` (the atomics
    /// themselves are neither `Clone` into a comparable value nor `PartialEq`).
    fn pos_weight_counts(allele_ref: &AlleleRef) -> Vec<[i32; 4]> {
        allele_ref.pos_weight.iter().map(|w| [w.get(0), w.get(1), w.get(2), w.get(3)]).collect()
    }

    #[test]
    fn assign_read_returns_none_for_empty_overlaps() {
        let refs = vec![AlleleRef::new(b"ACGTACGT".to_vec(), None)];
        let result =
            assign_read(b"ACGTACGT", &[], &refs, 0.8, 1, &mut crate::align_algo::DpCache::new());
        assert!(result.is_none());
    }

    #[test]
    fn assign_read_marks_base_coverage_for_matched_positions() {
        let consensus = b"ACGTACGTACGT".to_vec();
        let refs = vec![AlleleRef::new(consensus.clone(), None)];
        let read = b"ACGTACGTACGT";
        let core = simple_overlap(0, 0, 11, 0, 11, 1);

        let result = assign_read(
            read,
            std::slice::from_ref(&core),
            &refs,
            0.8,
            1,
            &mut crate::align_algo::DpCache::new(),
        );
        assert!(result.is_some());
        // Every position should now have nonzero coverage for its own base.
        for (i, &base) in consensus.iter().enumerate() {
            let base_code = nuc_to_num(base).unwrap();
            assert!(
                refs[0].pos_weight[i].get(base_code) > 0,
                "position {i} (base {base}) should have coverage"
            );
        }
    }

    // --- assign_reads_parallel ---

    /// Writes a small reference FASTA to a tempfile and builds a
    /// [`RefKmerFilter`] over it at `kmer_length` -- [`RefKmerFilter::
    /// from_reference_fasta`] only accepts a path (see `ref_kmer_filter`'s
    /// module doc for why), so every `assign_reads_parallel` test needs this
    /// small setup helper.
    fn build_test_filter(records: &[(&str, &str)], kmer_length: usize) -> RefKmerFilter {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        {
            use std::io::Write as _;
            let mut f = tmp.reopen().unwrap();
            for (name, seq) in records {
                writeln!(f, ">{name}\n{seq}").unwrap();
            }
        }
        RefKmerFilter::from_reference_fasta(tmp.path(), kmer_length).unwrap()
    }

    #[test]
    fn assign_reads_parallel_threads_1_matches_manual_sequential_loop() {
        // `threads <= 1` must be EXACTLY the same as calling
        // get_overlaps_from_read + assign_read by hand in a loop (the prior
        // driver code, before this function existed) -- a regression here
        // would mean the "sequential fast path" silently changed behavior.
        let consensus = b"ACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTG".to_vec();
        let filter = build_test_filter(&[("only", std::str::from_utf8(&consensus).unwrap())], 11);

        let seq1 = b"ACGTACGTGGATTACAGATTACAGATTACAG".to_vec();
        let seq2 = b"ATTACAGCCCTGACGTGTGACGTGTGACGTG".to_vec();
        let mut sorted_seqs: Vec<&[u8]> = vec![&seq1, &seq2];
        sorted_seqs.sort_unstable();

        // Manual sequential loop (mirrors the pre-P2 driver code).
        let manual_refs = vec![AlleleRef::new(consensus.clone(), None)];
        let mut manual_dp_cache = crate::align_algo::DpCache::new();
        let manual: Vec<Option<Vec<ExtendedOverlap>>> = sorted_seqs
            .iter()
            .map(|&seq| {
                let raw = filter
                    .get_overlaps_from_read(seq, &mut crate::ref_kmer_filter::Scratch::default())
                    .unwrap_or_default();
                assign_read(seq, &raw, &manual_refs, 0.8, 1, &mut manual_dp_cache)
            })
            .collect();

        // assign_reads_parallel at threads=1.
        let parallel_refs = vec![AlleleRef::new(consensus, None)];
        let via_helper =
            assign_reads_parallel(&filter, &sorted_seqs, &parallel_refs, 0.8, |_| 1, 1);

        assert_eq!(manual, via_helper, "threads=1 must match the manual sequential loop exactly");
        assert_eq!(
            pos_weight_counts(&manual_refs[0]),
            pos_weight_counts(&parallel_refs[0]),
            "threads=1 must mutate pos_weight identically to the manual sequential loop"
        );
    }

    #[test]
    fn assign_reads_parallel_threads_4_matches_threads_1_byte_identical() {
        // The core P2 invariant: parallelizing get_overlaps_from_read across
        // several distinct sequences must not change either the returned
        // ExtendedOverlaps OR the sequentially-applied pos_weight mutation,
        // regardless of how many rayon workers computed which overlap.
        let consensus = b"ACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATGGTACGTACGTACGTACG".to_vec();
        let filter = build_test_filter(&[("only", std::str::from_utf8(&consensus).unwrap())], 11);

        // Several distinct 40bp windows of the reference, so sorted_seqs has
        // more entries than any thread count tested below.
        let windows: Vec<Vec<u8>> = (0..8)
            .map(|i| {
                let start = i * 15;
                consensus[start..start + 40].to_vec()
            })
            .collect();
        let mut sorted_seqs: Vec<&[u8]> = windows.iter().map(Vec::as_slice).collect();
        sorted_seqs.sort_unstable();
        sorted_seqs.dedup();

        let refs_t1 = vec![AlleleRef::new(consensus.clone(), None)];
        let result_t1 = assign_reads_parallel(&filter, &sorted_seqs, &refs_t1, 0.8, |_| 1, 1);

        let refs_t4 = vec![AlleleRef::new(consensus, None)];
        let result_t4 = assign_reads_parallel(&filter, &sorted_seqs, &refs_t4, 0.8, |_| 1, 4);

        assert_eq!(result_t1, result_t4, "threads=1 vs threads=4 ExtendedOverlaps must match");
        assert_eq!(
            pos_weight_counts(&refs_t1[0]),
            pos_weight_counts(&refs_t4[0]),
            "threads=1 vs threads=4 pos_weight mutations must match exactly"
        );
        assert!(
            result_t1.iter().any(Option::is_some),
            "sanity: at least one window should have produced overlaps"
        );
    }

    // --- read_assignment_to_fragment_assignment ---

    fn extended(
        seq_idx: u32,
        read_start: i32,
        read_end: i32,
        seq_start: i32,
        seq_end: i32,
        strand: i8,
    ) -> ExtendedOverlap {
        ExtendedOverlap {
            seq_idx,
            read_start,
            read_end,
            seq_start,
            seq_end,
            strand,
            match_cnt: 2 * (read_end - read_start + 1),
            similarity: 1.0,
            relaxed_match_cnt: 2 * (read_end - read_start + 1),
            left_clip: 0,
            right_clip: 0,
        }
    }

    #[test]
    fn read_assignment_to_fragment_assignment_single_end_passthrough() {
        let overlaps1 = vec![extended(0, 0, 9, 0, 9, 1)];
        let assign = read_assignment_to_fragment_assignment(
            &overlaps1,
            None,
            false,
            31,
            |_| 100,
            |_, _, _| false,
        );
        assert_eq!(assign.len(), 1);
        assert_eq!(assign[0].seq_idx, 0);
        assert!(!assign[0].has_mate_pair);
        #[allow(clippy::float_cmp)]
        let qual_is_one = assign[0].qual == 1.0;
        assert!(qual_is_one);
    }

    #[test]
    fn read_assignment_to_fragment_assignment_combines_compatible_mate_pair() {
        // Mate 1 forward strand at [0,9], mate 2 reverse strand starting
        // further along the allele -- a compatible, non-overlapping pair.
        let overlaps1 = vec![extended(0, 0, 9, 0, 9, 1)];
        let overlaps2 = vec![extended(0, 0, 9, 50, 59, -1)];
        let assign = read_assignment_to_fragment_assignment(
            &overlaps1,
            Some(&overlaps2),
            false,
            31,
            |_| 100,
            |_, _, _| false,
        );
        assert_eq!(assign.len(), 1);
        assert!(assign[0].has_mate_pair);
        assert_eq!(assign[0].seq_start, 0);
        assert_eq!(assign[0].seq_end, 59);
    }

    #[test]
    fn read_assignment_to_fragment_assignment_incompatible_strand_yields_no_fragment() {
        // Both mates on the SAME strand against the same allele: NOT a
        // compatible pair, and (unlike the "one side has zero overlaps"
        // branch) this "both sides non-empty" branch has no dangling
        // fallback -- zero fragments are produced for this allele, matching
        // `SeqSet::ReadAssignmentToFragmentAssignment`'s `else` branch
        // (`SeqSet.hpp:2348-2381`) exactly: it only ever pushes a fragment
        // when the strand/position compatibility check passes.
        let overlaps1 = vec![extended(0, 0, 9, 0, 9, 1)];
        let overlaps2 = vec![extended(0, 0, 9, 50, 59, 1)];
        let assign = read_assignment_to_fragment_assignment(
            &overlaps1,
            Some(&overlaps2),
            false,
            3,
            |_| 100,
            |_, _, _| false,
        );
        assert!(assign.is_empty());
    }

    #[test]
    fn read_assignment_to_fragment_assignment_one_sided_empty_falls_back_to_dangling() {
        // Mate 2 has NO overlaps at all (e.g. it failed to align anywhere):
        // this DOES take the dangling fallback branch (`SeqSet.hpp:2330-
        // 2347`), producing a single dangling (non-mate-paired,
        // `o1FromR2 == false`) fragment for mate 1's overlap.
        let overlaps1 = vec![extended(0, 0, 9, 0, 9, 1)];
        let overlaps2: Vec<ExtendedOverlap> = Vec::new();
        let assign = read_assignment_to_fragment_assignment(
            &overlaps1,
            Some(&overlaps2),
            false,
            3,
            |_| 100,
            |_, _, _| false,
        );
        assert_eq!(assign.len(), 1);
        assert!(!assign[0].has_mate_pair);
        assert!(!assign[0].o1_from_r2);
    }

    #[test]
    fn read_assignment_to_fragment_assignment_dangling_short_span_is_dropped() {
        // Same as above, but hit_len_required=31 (the real SeqSet default):
        // the dangling-mate-pair filter's `span >= 3 * hitLenRequired`
        // check must reject this short overlap entirely.
        let overlaps1 = vec![extended(0, 0, 9, 0, 9, 1)];
        let overlaps2 = vec![extended(0, 0, 9, 50, 59, 1)];
        let assign = read_assignment_to_fragment_assignment(
            &overlaps1,
            Some(&overlaps2),
            false,
            31,
            |_| 100,
            |_, _, _| false,
        );
        assert!(assign.is_empty());
    }

    #[test]
    fn read_assignment_to_fragment_assignment_dominating_matchcnt_filters_even_if_assigned() {
        // Regression (HLA-U U*01:03 over-assignment): SeqSet.hpp:2601-2609's
        // "better read 1" truncated-mate-pair filter applies the "unused"
        // (not-already-assigned) guard ONLY to the `matchCnt ==` branch, NOT
        // the `matchCnt >` branch -- an overlap that strictly dominates the
        // representative on matchCnt triggers the filter even when it is
        // already assigned. A prior port wrote `dominates && unused` (gating
        // the `>` branch on `unused` too), so an already-assigned dominating
        // overlap failed to fire the filter and the fragment was over-kept.
        //
        // allele 0: representative mate pair, fragment matchCnt = 42+42 = 84.
        // allele 1: its mate-1 overlap (matchCnt 62) STRICTLY DOMINATES allele
        //   0's mate-1 (42); allele 1 also forms its own (weaker, total 74)
        //   fragment so its seq_idx IS in the assign map ("used"); and its
        //   mate-1 projects beyond the small consensus_len (100) so it is
        //   "truncated". The correct port fires the filter and clears assign.
        let overlaps1 = vec![
            extended(0, 0, 20, 0, 20, 1), // allele 0 mate1, mc=42 (representative)
            extended(1, 0, 30, 0, 30, 1), // allele 1 mate1, mc=62 > 42, truncated
        ];
        let overlaps2 = vec![
            extended(0, 0, 20, 100, 120, -1), // allele 0 mate2, mc=42
            extended(1, 0, 5, 200, 205, -1),  // allele 1 mate2, mc=12 -> frag1 total 74 < 84
        ];
        let assign = read_assignment_to_fragment_assignment(
            &overlaps1,
            Some(&overlaps2),
            false,
            3,
            |_| 100, // small consensus_len => allele 1's mate1 is truncated
            |_, _, _| false,
        );
        assert!(
            assign.is_empty(),
            "an already-assigned overlap that strictly dominates the representative on \
             matchCnt must trigger the truncated-mate-pair filter (SeqSet.hpp:2601)"
        );
    }

    // --- select_alleles_for_genes / get_allele_description (small
    // end-to-end scenario) ---

    /// Builds a tiny two-allele-per-gene, single-gene Genotyper scenario
    /// with a clear heterozygous call, entirely through the public
    /// read-assignment API (no direct field poking) -- exercises
    /// [`Genotyper::select_alleles_for_genes`] and
    /// [`Genotyper::get_allele_description`] together.
    #[test]
    fn select_alleles_for_genes_calls_clear_heterozygous_pair() {
        let names: Vec<String> =
            ["A*01:01:01", "A*01:02:01"].iter().map(|s| (*s).to_string()).collect();
        let consensus: Vec<Vec<u8>> = vec![
            b"AAAACGTACGTACGTACGTACGTACGTACGTACGTAAAA".to_vec(),
            b"TTTTCGTACGTACGTACGTACGTACGTACGTACGTTTTT".to_vec(),
        ];
        let weight = vec![1; 2];
        let mut eff_len = vec![40, 40];

        let mut g = Genotyper::new();
        g.init_allele_info(&names, &consensus, &weight, &mut eff_len, 8);
        g.init_read_assignments(20, 2000);

        // 10 reads uniquely supporting allele 0, 10 uniquely supporting allele 1.
        for i in 0..10 {
            let frag = FragmentOverlap {
                seq_idx: 0,
                seq_start: 0,
                seq_end: 10,
                match_cnt: 22,
                relaxed_match_cnt: 22,
                similarity: 1.0,
                has_mate_pair: false,
                o1_from_r2: false,
                qual: 1.0,
                has_n: false,
            };
            g.set_read_assignments(i, &[frag], 0.8, |_, _, _| false);
        }
        for i in 10..20 {
            let frag = FragmentOverlap {
                seq_idx: 1,
                seq_start: 0,
                seq_end: 10,
                match_cnt: 22,
                relaxed_match_cnt: 22,
                similarity: 1.0,
                has_mate_pair: false,
                o1_from_r2: false,
                qual: 1.0,
                has_n: false,
            };
            g.set_read_assignments(i, &[frag], 0.8, |_, _, _| false);
        }

        g.coalesce_read_assignments(0, 19);
        g.finalize_read_assignments(&[0, 0]);
        g.quantify_allele_equivalent_class(&eff_len, &weight);
        g.remove_low_likelihood_allele_in_equivalent_class(|idx| eff_len[idx]);
        g.select_alleles_for_genes(|idx| weight[idx]);

        let (allele1, allele2, secondary, called_cnt) = g.get_allele_description(0);
        assert_eq!(called_cnt, 2, "both alleles should be called");
        assert!(allele1.starts_with("A*01:01") || allele1.starts_with("A*01:02"));
        assert!(allele2.starts_with("A*01:01") || allele2.starts_with("A*01:02"));
        assert_ne!(allele1, allele2, "the two called major alleles must differ");
        assert!(secondary.is_empty(), "no secondary candidates expected in this scenario");
    }

    #[test]
    fn get_gene_allele_types_reflects_max_selected_rank() {
        let mut g = Genotyper::new();
        g.gene_cnt = 1;
        g.selected_alleles = vec![vec![(0, 0), (1, 1), (2, 2)]];
        assert_eq!(g.get_gene_allele_types(0), 3);
    }

    #[test]
    fn get_gene_allele_types_zero_when_nothing_selected() {
        let mut g = Genotyper::new();
        g.gene_cnt = 1;
        g.selected_alleles = vec![Vec::new()];
        assert_eq!(g.get_gene_allele_types(0), 0);
    }
}
