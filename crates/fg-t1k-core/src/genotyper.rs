//! Genotyper foundation, ported from T1K's `Genotyper`
//! (`vendor/t1k/Genotyper.hpp`).
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
//!   `crates/fg-t1k-sys/tests/diff_genotyper_em.rs` for the differential
//!   test (structure exact, abundances within tolerance) and its
//!   `assert_close_within_tolerance` doc comment for the full justification;
//!   the genotype CALL itself (Phase 5c, end-to-end) is the ultimate
//!   correctness gate, not bitwise abundance parity.
//!
//! # Deliberately NOT ported here (deferred to 5c)
//!
//! - `SelectAllele`/genotype-quality scoring/output formatting -- Phase 5c.
//! - `ReadAssignmentToFragmentAssignment` (`SeqSet.hpp:2310-2556`, builds
//!   `_fragmentOverlap`s from per-mate `_overlap` lists) -- this belongs to
//!   the read-processing loop that calls `SetReadAssignments` in a loop
//!   (`Genotyper.cpp:160-192`); [`Genotyper::set_read_assignments`] here
//!   takes already-built [`FragmentOverlap`]s as input, matching the C++
//!   method signature exactly, so this port is agnostic to how its caller
//!   produces them.
//! - `RemoveLowLikelihoodAlleleInEquivalentClass`/`InitAlleleAbundance` --
//!   neither is called from `Genotyper.cpp:main`'s stock pipeline between
//!   `FinalizeReadAssignments` and `QuantifyAlleleEquivalentClass`
//!   (`Genotyper.cpp:634,644`), so they are out of the 5b slice.

use std::collections::HashMap;

use crate::kmer_count::KmerCount;

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

/// The Phase-5a slice of T1K's `Genotyper` (`vendor/t1k/Genotyper.hpp`): the
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
        }
    }

    /// Ported from `Genotyper::SetFilterFrac` (`Genotyper.hpp:528-531`).
    pub fn set_filter_frac(&mut self, f: f64) {
        self.filter_frac = f;
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
    /// (`SeqSet.hpp:2305-2308`, `487-498`): a caller-supplied predicate so
    /// this port stays independent of a full `SeqSet::_seqWrapper::separator`
    /// port; see [`Genotyper::init_allele_info`]'s doc comment for the same
    /// "caller supplies the reference-derived input" pattern.
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
    /// by `crates/fg-t1k-sys/tests/diff_genotyper_em.rs`, which asserts this
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
    /// `crates/fg-t1k-sys/tests/diff_genotyper_em.rs`'s
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
    /// `crates/fg-t1k-sys/tests/diff_genotyper_em.rs`, verifies) agreement
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
}
