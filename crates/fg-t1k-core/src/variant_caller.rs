//! Post-analysis novel-variant caller, ported from T1K's `VariantCaller`
//! (`vendor/t1k/VariantCaller.hpp`): builds a per-allele base pileup from
//! reads assigned to the CALLED alleles, detects SNV variants against each
//! allele's consensus, resolves ambiguous multi-variant groups via read-level
//! phasing, and formats the result as `_allele.vcf` text
//! ([`VariantCaller::output_allele_vcf`]).
//!
//! # Scope: the 6a slice
//!
//! This module ports the `VariantCaller` CORE: pileup accumulation
//! ([`VariantCaller::update_base_variant_from_overlap`], mirroring
//! `UpdateBaseVariantFromOverlap`, `VariantCaller.hpp:103-173`),
//! `SetSeqAbundance`/`SetMaxVarGroupToResolve`
//! ([`VariantCaller::set_seq_abundance`]/[`VariantCaller::set_max_var_group_to_resolve`],
//! `VariantCaller.hpp:249-270`), and the full `ComputeVariant` driver
//! ([`VariantCaller::compute_variant`], `VariantCaller.hpp:978-1145`):
//! candidate-variant discovery, iterative candidate-variant expansion via
//! fragment overlaps, variant-group construction, per-group phasing
//! (`SolveVariantGroup`/`EnumerateVariants`), and `OutputAlleleVCF`
//! (`VariantCaller.hpp:1202-1227`).
//!
//! `VariantCaller::compute_variant` still takes already-assigned
//! `read1`/`read2`/`fragmentAssignments` as plain parameters (mirroring
//! `ComputeVariant`'s own signature exactly), the same way
//! [`crate::genotyper`] takes already-computed overlaps rather than calling
//! `GetOverlapsFromRead` itself -- the DRIVER that produces those inputs
//! (reference loading restricted to the called alleles, read I/O, read
//! reassignment) lives in `crates/fg-t1k/src/stages/analyze.rs` (the CLI
//! layer), matching `stages/genotype.rs`'s own split of "core in
//! `fg-t1k-core`, driver in the CLI crate".
//!
//! # Scope: the 6b slice
//!
//! [`add_overlap_alignment_info`]/[`add_fragment_alignment_info`] -- ported
//! from `SeqSet::AddOverlapAlignmentInfo`/`AddFragmentAlignmentInfo`
//! (`SeqSet.hpp:2657-2680,2758-2778`): the REQUIRED pre-`ComputeVariant` step
//! that populates every fragment overlap's `align` (see
//! [`add_fragment_alignment_info`]'s doc comment for why this is not merely
//! an optimization). [`exonic_position`] -- ported from
//! `SeqSet::GetExonicPosition` (`SeqSet.hpp:2808-2828`): the genomic-to-exon-
//! relative coordinate remap [`VariantCaller::output_allele_vcf`]'s VCF `POS`
//! column needs.
//!
//! # No indel calling (matches stock T1K exactly)
//!
//! `UpdateBaseVariantFromOverlap`'s indel-handling branch is commented-out
//! dead code in the vendored source (`VariantCaller.hpp:161-163`, `//TODO:
//! handle indels`) -- `EDIT_INSERT`/`EDIT_DELETE` alignment ops only advance
//! `refPos`/`readPos` bookkeeping here, exactly as in C++; no pileup count is
//! ever incremented for them, and [`VariantCaller::find_candidate_variants`]
//! only ever detects single-base substitutions against the allele consensus.
//! This is a faithful port of stock behavior, not a gap in this port.
//!
//! # `_overlap::align` is always caller-provided in the real pipeline
//!
//! In stock T1K, `Analyzer::main` populates every fragment assignment's
//! `overlap1.align`/`overlap2.align` via
//! `SeqSet::AddFragmentAlignmentInfo`/`AddOverlapAlignmentInfo`
//! (`SeqSet.hpp:2657-2680,2758-2778`, `Analyzer.cpp:622-633`) BEFORE calling
//! `VariantCaller::ComputeVariant` -- so `UpdateBaseVariantFromOverlap`'s own
//! `align == NULL` fallback (`VariantCaller.hpp:116-123`, recomputing the
//! alignment via `AlignAlgo::GlobalAlignment` from `seqStart`/`seqEnd`/
//! `readStart`/`readEnd`) is dead code on that path (confirmed: every
//! `fragment_assignment` the real `fg-t1k analyze` driver builds has already
//! had [`add_fragment_alignment_info`] called on it -- see that function's
//! own doc comment for why this is a REQUIRED step, not just the lazy
//! optimization `update_base_variant_from_overlap`'s fallback name might
//! suggest). This port still implements both branches faithfully (see
//! [`VariantCaller::update_base_variant_from_overlap`]) since 6a's
//! [`Overlap::align`] is modeled as `Option<Vec<i8>>` for completeness and
//! testability without requiring the full `fg-t1k analyze` driver to
//! pre-populate it in unit tests.
//!
//! # `f64` counts, not `f32`
//!
//! `_baseVariant::count`/`uniqCount`/`unweightedCount` are `double` in C++
//! (`VariantCaller.hpp:24-26`) -- matched here as `f64`, unlike some other
//! ported structs in this crate that narrow to `f32` (e.g.
//! [`crate::genotyper::ReadAssignment::weight`]).

use crate::genotyper::{AlleleRef, Genotyper};

/// Ported from `_pairIntDouble` (`defs.h:34-38`): a generic `(int, double)`
/// pair. Used here for `_baseVariant::alignInfo` -- the best `(matchCnt,
/// similarity)` seen so far for each nucleotide at a position.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairIntDouble {
    pub a: i32,
    pub b: f64,
}

impl Default for PairIntDouble {
    fn default() -> Self {
        Self { a: 0, b: 0.0 }
    }
}

/// Ported from `_overlap` (`SeqSet.hpp:89-115`) as consumed by
/// `VariantCaller` -- i.e. [`crate::genotyper::ExtendedOverlap`]'s fields
/// plus the `align` op-sequence array C++'s `_overlap::align` carries (which
/// `ExtendedOverlap` omits; see that type's doc comment). `align = None`
/// mirrors `align == NULL`: [`VariantCaller::update_base_variant_from_overlap`]
/// recomputes it via [`crate::align_algo::global_alignment`], exactly as
/// `UpdateBaseVariantFromOverlap` does (`VariantCaller.hpp:116-123`).
#[derive(Debug, Clone, PartialEq)]
pub struct Overlap {
    /// `_overlap::seqIdx`. `-1` mirrors "no overlap"
    /// (`UpdateBaseVariantFromOverlap`'s `o.seqIdx == -1` early return).
    pub seq_idx: i32,
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
    /// `_overlap::matchCnt`.
    pub match_cnt: i32,
    /// `_overlap::similarity`.
    pub similarity: f64,
    /// `_overlap::align` -- the alignment op sequence (`EDIT_MATCH` /
    /// `EDIT_MISMATCH` / `EDIT_INSERT` / `EDIT_DELETE`), or `None` to mirror
    /// `align == NULL` (recomputed on demand).
    pub align: Option<Vec<i8>>,
}

impl Overlap {
    /// The sentinel "no overlap" value, mirroring a `_fragmentOverlap`'s
    /// unused `overlap2` when `hasMatePair == false` (never dereferenced by
    /// [`VariantCaller::update_base_variant_from_overlap`], which early-returns
    /// on `seqIdx == -1`).
    #[must_use]
    pub fn none() -> Self {
        Self {
            seq_idx: -1,
            read_start: 0,
            read_end: 0,
            seq_start: 0,
            seq_end: 0,
            strand: 1,
            match_cnt: 0,
            similarity: 0.0,
            align: None,
        }
    }
}

/// Ported from `_fragmentOverlap` (`SeqSet.hpp:146-173`) as consumed by
/// `VariantCaller` -- unlike [`crate::genotyper::FragmentOverlap`] (which
/// omits `overlap1`/`overlap2` since 5a's `ReadAssignmentWeight`/
/// `SetReadAssignments` never read them), this port's variant caller reads
/// `overlap1`/`overlap2` directly (they carry the per-mate alignment used
/// for pileup accumulation), so this is a distinct, additive type rather
/// than an extension of the 5a struct.
#[derive(Debug, Clone, PartialEq)]
pub struct FragmentOverlap {
    /// `_fragmentOverlap::seqIdx`.
    pub seq_idx: i32,
    /// `_fragmentOverlap::hasMatePair`.
    pub has_mate_pair: bool,
    /// `_fragmentOverlap::o1FromR2`.
    pub o1_from_r2: bool,
    /// `_fragmentOverlap::overlap1`.
    pub overlap1: Overlap,
    /// `_fragmentOverlap::overlap2`.
    pub overlap2: Overlap,
}

/// Ported from `_variant` (`VariantCaller.hpp:7-20`): one called novel
/// variant against an allele's consensus.
#[derive(Debug, Clone, PartialEq)]
pub struct Variant {
    /// `_variant::seqIdx` -- which allele this variant is against.
    pub seq_idx: i32,
    /// `_variant::refStart`.
    pub ref_start: i32,
    /// `_variant::refEnd`.
    pub ref_end: i32,
    /// `_variant::ref` -- always a single base in this port (see module
    /// docs: no indel calling).
    pub reference: u8,
    /// `_variant::var` -- always a single base.
    pub var: u8,
    /// `_variant::allSupport`.
    pub all_support: f64,
    /// `_variant::varSupport`.
    pub var_support: f64,
    /// `_variant::varUniqSupport`.
    pub var_uniq_support: f64,
    /// `_variant::varGroupId`.
    pub var_group_id: i32,
    /// `_variant::outputGroupId` -- `0` best variants, `1` equal-best
    /// variants.
    pub output_group_id: i32,
    /// `_variant::qual`.
    pub qual: i32,
}

/// Ported from `_baseVariant` (`VariantCaller.hpp:22-55`): per-(allele,
/// position) pileup accumulator.
#[derive(Debug, Clone)]
struct BaseVariant {
    count: [f64; 4],
    uniq_count: [f64; 4],
    unweighted_count: [f64; 4],
    align_info: [PairIntDouble; 4],
    exon: bool,
    /// `-1` not a variant candidate, otherwise the index into
    /// `candidate_variants`.
    candidate_id: i32,
    /// The id in the final variant table (`finalVariants`).
    final_variant_ids: Vec<usize>,
}

impl Default for BaseVariant {
    fn default() -> Self {
        Self {
            count: [0.0; 4],
            uniq_count: [0.0; 4],
            unweighted_count: [0.0; 4],
            align_info: [PairIntDouble::default(); 4],
            exon: false,
            candidate_id: -1,
            final_variant_ids: Vec::new(),
        }
    }
}

impl BaseVariant {
    /// Mirrors `_baseVariant::AllCountSum` (`VariantCaller.hpp:32-35`).
    fn all_count_sum(&self) -> f64 {
        self.count[0] + self.count[1] + self.count[2] + self.count[3]
    }

    /// Mirrors `_baseVariant::UnweightedCountSum` (`VariantCaller.hpp:42-45`).
    fn unweighted_count_sum(&self) -> f64 {
        self.unweighted_count[0]
            + self.unweighted_count[1]
            + self.unweighted_count[2]
            + self.unweighted_count[3]
    }

    /// Mirrors `_baseVariant::IsGoodAssignment` (`VariantCaller.hpp:47-54`).
    /// NOTE: the C++ return type is `double` but the body only ever returns
    /// `bool` literals (an implicit `bool`->`double` conversion at every
    /// `return` -- `false` becomes `0.0`, `true` becomes `1.0`); every call
    /// site uses the result in a boolean context (`if (...)` /
    /// `!IsGoodAssignment(...)`), so this port returns `bool` directly.
    fn is_good_assignment(&self, match_cnt: i32) -> bool {
        for align in &self.align_info {
            if match_cnt < align.a - 4 {
                return false;
            }
        }
        true
    }
}

/// Ported from `_adjBaseVariantToBaseVariant` (`VariantCaller.hpp:75-81`):
/// one node/edge in the candidate-variant adjacency list used to group
/// co-occurring candidate variants.
#[derive(Debug, Clone, Copy)]
struct AdjVarToVar {
    var_idx: i32,
    weight: f64,
    root_candidate: bool,
    next: i32,
}

/// Ported from `_adjFragmentToBaseVariant` (`VariantCaller.hpp:57-65`): one
/// edge from a fragment (read pair) to a candidate variant it supports, with
/// the nucleotide it observed there.
///
/// `nuc` is written (mirroring `strcpy(nFragToBaseVar.nuc, var)`,
/// `VariantCaller.hpp:664`) but never read anywhere in stock T1K either --
/// `grep -n "adjFrag\[.*\]\.nuc" vendor/t1k/VariantCaller.hpp` has zero
/// matches. Kept here (rather than dropped) for structural fidelity with the
/// C++ struct this is a 1:1 port of.
#[derive(Debug, Clone, Copy)]
struct AdjFragToVar {
    #[allow(dead_code)]
    nuc: u8,
    next: i32,
}

/// Ported from `_adjBaseVariantToFragment` (`VariantCaller.hpp:67-72`): the
/// reverse edge (candidate variant -> fragment).
#[derive(Debug, Clone, Copy)]
struct AdjVarToFrag {
    frag_idx: i32,
    nuc: u8,
    next: i32,
}

/// Ported from `_enumVarResult` (`VariantCaller.hpp:83-89`): the best (and,
/// if tied, second-best/"equal best") nucleotide assignment found by
/// [`VariantCaller::enumerate_variants`] for one variant group.
#[derive(Debug, Clone, Default)]
struct EnumVarResult {
    best_cover: f64,
    used_var_cnt: i32,
    best_enum_variants: Vec<u8>,
    equal_best_enum_variants: Vec<u8>,
}

/// Maps `A`/`C`/`G`/`T` to a 0-3 index, mirroring the `nucToNum` table
/// (`Genotyper.cpp:37-40` and friends) as used by `VariantCaller.hpp`.
/// Duplicated from (rather than reusing) `genotyper`'s private,
/// identically-defined `nuc_to_num` -- this module stays additive-only on
/// top of prior-phase files, matching this crate's established convention
/// (see e.g. `genotyper.rs`'s own doc comment on why it duplicates
/// `ref_kmer_filter`'s `nuc_to_num`).
fn nuc_to_num(c: u8) -> Option<usize> {
    match c {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// The inverse of [`nuc_to_num`], mirroring `numToNuc` (`Analyzer.cpp:39`
/// and friends): `{'A', 'C', 'G', 'T'}`.
const NUM_TO_NUC: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Complements a single base: `A<->T`, `C<->G`, `N->N` -- mirrors
/// `numToNuc[3 - nucToNum[c - 'A']]` with `N` bypassing the table
/// (`SeqSet::ReverseComplement`, `SeqSet.hpp:2103-2114`). Duplicated from
/// (rather than reusing) `genotyper`'s private, identically-defined
/// `complement_base`/`reverse_complement` -- see this module's doc comment
/// on why 6a stays additive-only on top of prior-phase files.
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

/// `true` if position `pos` (0-based) falls within any of `exons`' `(start,
/// end)` inclusive intervals -- mirrors `_validDiff::exon`
/// (`SeqSet.hpp:657,670`) via [`crate::genotyper::AlleleRef::exons`], which
/// is `pub` (unlike `AlleleRef::is_exon` itself, which this module does not
/// call -- see this module's doc comment on staying additive-only).
fn is_exon(exons: &[(i32, i32)], pos: i32) -> bool {
    exons.iter().any(|&(s, e)| pos >= s && pos <= e)
}

/// Ported from `SeqSet::GetExonicPosition` (`SeqSet.hpp:2808-2828`): maps a
/// 0-based genomic (full-consensus) position to its 0-based exon-relative
/// position -- the sum of every EARLIER exon's length plus the offset into
/// the exon `pos` itself falls in -- or `-1` if `pos` is not in any exon
/// (mirrors `!IsPosInExon(seqIdx, pos)`'s early return, via [`is_exon`]).
///
/// This is what [`VariantCaller::output_allele_vcf`]'s `exonic_position`
/// parameter is expected to compute (see that function's doc comment); it is
/// a free function here (not a [`VariantCaller`] method) since it depends
/// only on one allele's `exons`, matching [`is_exon`]'s own shape.
///
/// # `ConvertVariantsToExonCoord` is dead code -- this is the ONLY
/// `GetExonicPosition` call site on the `OutputAlleleVCF` path
///
/// `VariantCaller::ConvertVariantsToExonCoord` (`VariantCaller.hpp:1189-
/// 1200`) would call `GetExonicPosition` a SECOND time, mutating
/// `finalVariants[i].refStart`/`refEnd` in place before `OutputAlleleVCF`
/// runs -- but its body is `return ;` as its very first statement (dead
/// code), and `grep -n "ConvertVariantsToExonCoord" vendor/t1k/*.cpp
/// vendor/t1k/*.hpp` has exactly one match (the dead definition itself): it
/// is never called from `ComputeVariant` or any `main`. So `finalVariants[i]
/// .refStart`/`refEnd` are always genomic (never exonic) coordinates by the
/// time [`VariantCaller::output_allele_vcf`] runs, and this function's own
/// `GetExonicPosition` call (`VariantCaller.hpp:1215`) is the only exonic
/// remapping that ever actually happens -- exactly what this port's
/// `exonic_position` closure parameter is for.
#[must_use]
pub fn exonic_position(exons: &[(i32, i32)], pos: i32) -> i32 {
    if !is_exon(exons, pos) {
        return -1;
    }
    let mut psum = 0;
    for &(s, e) in exons {
        if pos >= s && pos <= e {
            return psum + pos - s;
        }
        psum += e - s + 1;
    }
    psum
}

/// Ported from `SeqSet::AddOverlapAlignmentInfo` (`SeqSet.hpp:2657-2680`):
/// populates `o.align` (unconditionally overwriting any previous value, NOT
/// gated on `align == NULL` -- unlike
/// [`VariantCaller::update_base_variant_from_overlap`]'s own `align == NULL`
/// fallback, which is a DIFFERENT code path with a different caller; see
/// that function's doc comment) via [`crate::align_algo::global_alignment`]
/// over exactly the overlap's aligned sub-ranges: the allele consensus slice
/// `[seq_start, seq_end]` against the (possibly reverse-complemented, when
/// `o.strand == -1`) read slice `[read_start, read_end]`.
///
/// A no-op (mirrors the C++'s `if (o.seqIdx == -1) return ;`) when `o` has
/// no overlap.
///
/// # Panics
///
/// Panics if `o`'s coordinates are out of bounds for `read`/`allele_consensus`
/// (not expected: `o` is assumed to be a valid overlap against the same
/// reference/read this is called with, e.g. one of
/// [`crate::genotyper::read_assignment_to_fragment_assignment_with_overlaps`]'s
/// own `ExtendedOverlap`s converted to an [`Overlap`]).
pub fn add_overlap_alignment_info(read: &[u8], o: &mut Overlap, allele_consensus: &[u8]) {
    if o.seq_idx == -1 {
        return;
    }
    let rc;
    let r: &[u8] = if o.strand == -1 {
        rc = reverse_complement(read);
        &rc
    } else {
        read
    };
    let seq_start = usize::try_from(o.seq_start).unwrap();
    let seq_end = usize::try_from(o.seq_end).unwrap();
    let read_start = usize::try_from(o.read_start).unwrap();
    let read_end = usize::try_from(o.read_end).unwrap();
    let result = crate::align_algo::global_alignment(
        &allele_consensus[seq_start..=seq_end],
        &r[read_start..=read_end],
        crate::align_algo::DEFAULT_BAND,
    );
    o.align = Some(result.align);
}

/// Ported from `SeqSet::AddFragmentAlignmentInfo` (`SeqSet.hpp:2758-2778`):
/// for each fragment overlap, populates `overlap1`'s (and, for mate-pair
/// fragments, `overlap2`'s) `align` via [`add_overlap_alignment_info`] --
/// this is the `Analyzer` driver step (`Analyzer.cpp:622-633`'s
/// `AddFragmentAlignmentInfo` call, gated on `reads1[i].fragmentAssigned`
/// at that call site, not here) that MUST run before
/// [`VariantCaller::compute_variant`]: unlike
/// [`VariantCaller::update_base_variant_from_overlap`]'s own lazy `align ==
/// NULL` recompute fallback, [`VariantCaller::expand_candidate_variants_from_fragment_overlap`]
/// and [`VariantCaller::build_fragment_candidate_var_graph`] both `continue`/
/// no-op on `align == None` rather than recomputing it -- so a
/// [`FragmentOverlap`] whose `align` was never populated silently
/// contributes nothing to variant-group construction on the real pipeline.
///
/// `read2` is `None` for single-end fragments (every `fragment_overlap` in
/// `fragment_assignment` is then expected to have `has_mate_pair == false`).
///
/// # Panics
///
/// Panics if any `fragment_assignment[i].seq_idx` is out of range for
/// `allele_consensus`, or if a mate-pair/`o1FromR2` fragment is passed with
/// `read2 == None` (not expected: every `fragment_assignment` here is
/// assumed to come from the same reference/read pair this is called with,
/// e.g. [`crate::genotyper::read_assignment_to_fragment_assignment_with_overlaps`]'s
/// own output).
pub fn add_fragment_alignment_info(
    read1: &[u8],
    read2: Option<&[u8]>,
    fragment_assignment: &mut [FragmentOverlap],
    allele_consensus: &[Vec<u8>],
) {
    for frag in fragment_assignment {
        let seq_idx = usize::try_from(frag.seq_idx).unwrap();
        let consensus = &allele_consensus[seq_idx];
        if frag.has_mate_pair {
            add_overlap_alignment_info(read1, &mut frag.overlap1, consensus);
            add_overlap_alignment_info(
                read2.expect("has_mate_pair implies read2 is present"),
                &mut frag.overlap2,
                consensus,
            );
        } else if !frag.o1_from_r2 {
            add_overlap_alignment_info(read1, &mut frag.overlap1, consensus);
        } else {
            add_overlap_alignment_info(
                read2.expect("o1_from_r2 implies read2 is present"),
                &mut frag.overlap1,
                consensus,
            );
        }
    }
}

/// Ported from `VariantCaller` (`VariantCaller.hpp:92-1312`). See module
/// docs for the exact 6a slice.
pub struct VariantCaller {
    /// `baseVariants[seqIdx][refPos]`.
    base_variants: Vec<Vec<BaseVariant>>,
    /// `seqAbundance[seqIdx]`.
    seq_abundance: Vec<f64>,
    /// `seqCopy[seqIdx]` -- 1-homozygous, 2-heterozygous (i.e. the number of
    /// alleles sharing that allele's gene).
    seq_copy: Vec<i32>,
    /// `candidateVariants[i] = (seqIdx, refPos)`.
    candidate_variants: Vec<(i32, i32)>,
    /// `candidateVariantGroupId[i]`.
    candidate_variant_group_id: Vec<i32>,
    /// `finalVariants`.
    final_variants: Vec<Variant>,
    /// `maxVarGroupToResolve`.
    max_var_group_to_resolve: i32,
}

impl VariantCaller {
    /// Ported from `VariantCaller::VariantCaller(SeqSet &inRefSeq)`
    /// (`VariantCaller.hpp:229-246`): allocates one zeroed [`BaseVariant`]
    /// per (allele, position), marking each position's `exon` flag from
    /// `allele_refs`. `names.len() == allele_refs.len()` is expected (both
    /// indexed by `seqIdx`), but `names` itself is not stored here -- see
    /// [`VariantCaller::output_allele_vcf`] for where it is threaded through
    /// (mirrors `refSet.GetSeqName` being a `SeqSet` accessor rather than
    /// `VariantCaller` state).
    #[must_use]
    pub fn new(allele_refs: &[AlleleRef]) -> Self {
        let base_variants = allele_refs
            .iter()
            .map(|a| {
                let len = a.consensus.len();
                (0..len)
                    .map(|pos| {
                        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                        let pos_i32 = pos as i32;
                        BaseVariant { exon: is_exon(&a.exons, pos_i32), ..BaseVariant::default() }
                    })
                    .collect()
            })
            .collect();
        Self {
            base_variants,
            seq_abundance: Vec::new(),
            seq_copy: Vec::new(),
            candidate_variants: Vec::new(),
            candidate_variant_group_id: Vec::new(),
            final_variants: Vec::new(),
            max_var_group_to_resolve: 8,
        }
    }

    /// Ported from `VariantCaller::SetSeqAbundance` (`VariantCaller.hpp:249-265`):
    /// copies each allele's abundance from `genotyper`, and derives
    /// `seqCopy[i]` -- the count of alleles sharing allele `i`'s gene (i.e.
    /// `geneAlleleCount[GetAlleleGeneIdx(i)]`, accumulated via a
    /// `std::map<int,int>` exactly as C++ does; iteration ORDER of that map
    /// never matters here since only final per-gene counts are read, so a
    /// `HashMap` is a safe, faithful substitute).
    pub fn set_seq_abundance(&mut self, genotyper: &Genotyper, allele_count: usize) {
        self.seq_abundance =
            (0..allele_count).map(|i| genotyper.allele_info[i].abundance).collect();

        let mut gene_allele_count: std::collections::HashMap<i32, i32> =
            std::collections::HashMap::new();
        for i in 0..allele_count {
            *gene_allele_count.entry(genotyper.allele_info[i].gene_idx).or_insert(0) += 1;
        }
        self.seq_copy = (0..allele_count)
            .map(|i| gene_allele_count[&genotyper.allele_info[i].gene_idx])
            .collect();
    }

    /// Ported from `VariantCaller::SetMaxVarGroupToResolve`
    /// (`VariantCaller.hpp:267-270`).
    pub fn set_max_var_group_to_resolve(&mut self, m: i32) {
        self.max_var_group_to_resolve = m;
    }

    /// Ported from `VariantCaller::UpdateBaseVariantFromOverlap`
    /// (`VariantCaller.hpp:103-173`): walks one read's alignment to one
    /// allele (from `o`), adding weighted pileup counts for each
    /// matched/mismatched position. `filter_low_qual` gates the
    /// `IsGoodAssignment` low-quality-position skip (C++'s `filterLowQual`
    /// parameter).
    ///
    /// `read` is the ORIGINAL (not reverse-complemented) read sequence --
    /// this function reverse-complements internally when `o.strand == -1`,
    /// exactly as C++'s `r = strdup(read); refSet.ReverseComplement(r, read,
    /// readLen);` does.
    ///
    /// # Panics
    ///
    /// Panics if `o.seq_idx` is out of range for `self.base_variants`, or if
    /// `o`'s coordinates are out of bounds for `read`/the allele consensus
    /// (not expected: `o` is assumed to be a valid overlap against the same
    /// reference this [`VariantCaller`] was constructed from).
    fn update_base_variant_from_overlap(
        &mut self,
        read: &[u8],
        weight: f64,
        filter_low_qual: bool,
        o: &Overlap,
        allele_consensus: &[u8],
    ) {
        if o.seq_idx == -1 {
            return;
        }
        let seq_idx = usize::try_from(o.seq_idx).expect("seq_idx is non-negative");

        let rc;
        let r: &[u8] = if o.strand == -1 {
            rc = reverse_complement(read);
            &rc
        } else {
            read
        };

        // `align == NULL` fallback (VariantCaller.hpp:116-123): recompute via
        // AlignAlgo::GlobalAlignment when the caller did not pre-populate it
        // (see module docs -- always populated in the real end-to-end
        // pipeline, but this port supports both for testability).
        //
        // Note: C++ actually has a latent bug here -- it declares `signed
        // char *align = new signed char[...]` INSIDE the `if (align ==
        // NULL)` block, shadowing the outer `align` in an inner scope. The
        // freshly computed alignment is never assigned back to the outer
        // `align`, which remains NULL, so C++'s subsequent loop dereferences
        // a null pointer. This branch is dead code on the real end-to-end
        // path (`align` is always populated there), so the bug is never hit
        // in practice. This port intentionally does NOT reproduce the
        // shadowing bug: it uses the recomputed `align` directly. This
        // divergence from C++ is benign/beneficial and only affects the
        // test-only `align == None` path.
        let computed;
        let align: &[i8] = if let Some(a) = &o.align {
            a
        } else {
            let seq_start = usize::try_from(o.seq_start).unwrap();
            let seq_end = usize::try_from(o.seq_end).unwrap();
            let read_start = usize::try_from(o.read_start).unwrap();
            let read_end = usize::try_from(o.read_end).unwrap();
            computed = crate::align_algo::global_alignment(
                &allele_consensus[seq_start..=seq_end],
                &r[read_start..=read_end],
                crate::align_algo::DEFAULT_BAND,
            );
            &computed.align
        };

        let mut ref_pos = o.seq_start;
        let mut read_pos = o.read_start;
        for &op in align {
            if op == crate::align_algo::EDIT_MATCH || op == crate::align_algo::EDIT_MISMATCH {
                let ref_pos_usize = usize::try_from(ref_pos).unwrap();
                let bv = &mut self.base_variants[seq_idx][ref_pos_usize];
                if filter_low_qual && !bv.is_good_assignment(o.match_cnt) {
                    // fallthrough to position advancement below (`continue`
                    // in C++ skips the count update but NOT the ref/read pos
                    // advancement, which happens after this if-block).
                } else {
                    let rb = r[usize::try_from(read_pos).unwrap()];
                    if rb != b'N' {
                        if let Some(nuc_idx) = nuc_to_num(rb) {
                            #[allow(clippy::float_cmp)] // exact C++ `weight == 1` check.
                            if weight == 1.0 {
                                bv.uniq_count[nuc_idx] += weight;
                            }
                            bv.count[nuc_idx] += 1.0;
                            bv.unweighted_count[nuc_idx] += 1.0;

                            if o.match_cnt > bv.align_info[nuc_idx].a {
                                bv.align_info[nuc_idx].a = o.match_cnt;
                                bv.align_info[nuc_idx].b = o.similarity;
                            } else if o.match_cnt == bv.align_info[nuc_idx].a
                                && o.similarity > bv.align_info[nuc_idx].b
                            {
                                bv.align_info[nuc_idx].b = o.similarity;
                            }
                        }
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

    /// Ported from `VariantCaller::UpdateBaseVariantFromFragmentOverlap`
    /// (`VariantCaller.hpp:272-305`): distributes one fragment's total
    /// weight (`seqAbundance[seqIdx] / totalWeight` across ALL alleles this
    /// fragment was assigned to) or `alignInfo`-only mode (`updateType ==
    /// 1`, `weight = 0`, `filterLowQual = false`) into the pileup.
    fn update_base_variant_from_fragment_overlap(
        &mut self,
        read1: &[u8],
        read2: Option<&[u8]>,
        update_type: i32,
        fragment_assignment: &[FragmentOverlap],
        allele_consensus: &[Vec<u8>],
    ) {
        let mut total_weight = 0.0;
        for frag in fragment_assignment {
            let idx = usize::try_from(frag.seq_idx).unwrap();
            total_weight += self.seq_abundance[idx];
        }

        for frag in fragment_assignment {
            let seq_idx = usize::try_from(frag.seq_idx).unwrap();
            let mut weight = self.seq_abundance[seq_idx] / total_weight;
            let mut filter_low_qual = true;
            if update_type == 1 {
                filter_low_qual = false;
                weight = 0.0;
            }
            let consensus = &allele_consensus[seq_idx];
            if frag.has_mate_pair {
                self.update_base_variant_from_overlap(
                    read1,
                    weight,
                    filter_low_qual,
                    &frag.overlap1,
                    consensus,
                );
                self.update_base_variant_from_overlap(
                    read2.expect("hasMatePair implies read2 is present"),
                    weight,
                    filter_low_qual,
                    &frag.overlap2,
                    consensus,
                );
            } else if !frag.o1_from_r2 {
                self.update_base_variant_from_overlap(
                    read1,
                    weight,
                    filter_low_qual,
                    &frag.overlap1,
                    consensus,
                );
            } else {
                self.update_base_variant_from_overlap(
                    read2.expect("o1FromR2 implies read2 is present"),
                    weight,
                    filter_low_qual,
                    &frag.overlap1,
                    consensus,
                );
            }
        }
    }

    /// Ported from `VariantCaller::FindCandidateVariants`
    /// (`VariantCaller.hpp:307-345`): scans every (allele, position) pileup
    /// for a non-reference nucleotide meeting BOTH the absolute
    /// (`countThreshold = 5`) and relative (`>= refCount * 0.5`) support
    /// thresholds, recording the FIRST such nucleotide found (in `A, C, G,
    /// T` order, `VariantCaller.hpp:324`'s `for k in 0..4`) as a candidate.
    #[allow(clippy::needless_range_loop)]
    fn find_candidate_variants(&mut self, consensus: &[Vec<u8>]) {
        // Hoisted above the `clear()` call to satisfy clippy::items_after_statements.
        const COUNT_THRESHOLD: f64 = 5.0;
        self.candidate_variants.clear();
        let seq_cnt = self.base_variants.len();
        for i in 0..seq_cnt {
            let len = self.base_variants[i].len();
            let s = &consensus[i];
            let factor = 0.5;
            for j in 0..len {
                let ref_nuc_idx = nuc_to_num(s[j]).expect("consensus base is A/C/G/T");
                let ref_count = self.base_variants[i][j].count[ref_nuc_idx];
                for k in 0..4 {
                    let count_k = self.base_variants[i][j].count[k];
                    if count_k >= COUNT_THRESHOLD
                        && count_k >= ref_count * factor
                        && k != ref_nuc_idx
                    {
                        let id = self.candidate_variants.len();
                        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                        self.candidate_variants.push((i as i32, j as i32));
                        self.base_variants[i][j].candidate_id = i32::try_from(id).unwrap();
                        self.candidate_variant_group_id.push(-1);
                        break;
                    }
                }
            }
        }
    }

    /// Ported from `VariantCaller::ComputeCandidateVarAccuCount`
    /// (`VariantCaller.hpp:214-227`): a running count, per allele, of how
    /// many candidate variants have been seen up to (and including) each
    /// position -- lets [`Self::contain_candidate_var`] answer "any
    /// candidate variant in `[start, end]`?" in O(1).
    fn compute_candidate_var_accu_count(&self, seq_idx: usize) -> Vec<i32> {
        let len = self.base_variants[seq_idx].len();
        let mut accu = vec![0i32; len + 1];
        for i in 0..len {
            accu[i + 1] = accu[i] + i32::from(self.base_variants[seq_idx][i].candidate_id != -1);
        }
        accu
    }

    /// Ported from `VariantCaller::containCandidateVar`
    /// (`VariantCaller.hpp:182-199`).
    fn contain_candidate_var(start: i32, end: i32, accu: &[i32]) -> bool {
        let s = usize::try_from(start).unwrap();
        let e = usize::try_from(end).unwrap();
        accu[s] != accu[e + 1]
    }

    /// Ported from `VariantCaller::SelectOverlapFromFragmentOverlap`
    /// (`VariantCaller.hpp:201-208`).
    fn select_overlap(k: i32, frag: &FragmentOverlap) -> &Overlap {
        if k == 0 { &frag.overlap1 } else { &frag.overlap2 }
    }

    /// Ported from `VariantCaller::ExpandCandidateVariantsFromFragmentOverlap`
    /// (`VariantCaller.hpp:347-571`): for each read end (`k = 0, 1`) of each
    /// fragment assignment, if ANY assigned allele's overlap region already
    /// contains a known candidate variant, walk every OTHER assigned
    /// allele's aligned position at that SAME read offset and register it as
    /// a NEW candidate variant too (so downstream group resolution can
    /// consider all alleles' evidence at once) -- and accumulate
    /// variant-to-variant co-occurrence weight (`adjVarToVar`) for
    /// [`Self::build_candidate_variant_group`].
    #[allow(clippy::too_many_lines, clippy::similar_names)]
    fn expand_candidate_variants_from_fragment_overlap(
        &mut self,
        read1: &[u8],
        read2: Option<&[u8]>,
        fragment_assignment: &[FragmentOverlap],
        adj_var_to_var: &mut Vec<AdjVarToVar>,
        seq_candidate_accu_count: &[Vec<i32>],
    ) {
        if fragment_assignment.is_empty() {
            return;
        }
        let assign_cnt = fragment_assignment.len();

        for k in 0..=1 {
            if k == 1 && !fragment_assignment[0].has_mate_pair {
                break;
            }

            // Check whether there is a candidate variant in ANY assigned
            // allele's overlap region for this read end.
            let mut has_candidate = false;
            for frag in fragment_assignment {
                let o = Self::select_overlap(k, frag);
                let seq_idx = usize::try_from(o.seq_idx).unwrap();
                if Self::contain_candidate_var(
                    o.seq_start,
                    o.seq_end,
                    &seq_candidate_accu_count[seq_idx],
                ) {
                    has_candidate = true;
                    break;
                }
            }
            if !has_candidate {
                continue;
            }

            let read: &[u8] = if k == 1 || (k == 0 && fragment_assignment[0].o1_from_r2) {
                read2.expect("k==1 or o1FromR2 implies read2 is present")
            } else {
                read1
            };
            let len = i32::try_from(read.len()).unwrap();

            let mut ref_pos: Vec<i32> = Vec::with_capacity(assign_cnt);
            let mut read_pos0: Vec<i32> = Vec::with_capacity(assign_cnt);
            for frag in fragment_assignment {
                let o = Self::select_overlap(k, frag);
                ref_pos.push(o.seq_start);
                read_pos0.push(o.read_start);
            }

            // They all should have the same start position in read position.
            if (1..assign_cnt).any(|i| read_pos0[i] != read_pos0[0]) {
                continue;
            }

            let mut read_pos = read_pos0.clone();
            let mut align_idx = vec![0usize; assign_cnt];
            let mut valid_assignment = vec![false; assign_cnt];

            for j in 0..len {
                // Expand the set of candidate variants.
                let mut first_candidate_id = -1i32;
                for i in 0..assign_cnt {
                    let o = Self::select_overlap(k, &fragment_assignment[i]);
                    let seq_idx = usize::try_from(o.seq_idx).unwrap();
                    valid_assignment[i] =
                        if ref_pos[i] < i32::try_from(self.base_variants[seq_idx].len()).unwrap() {
                            let p = usize::try_from(ref_pos[i]).unwrap();
                            self.base_variants[seq_idx][p].is_good_assignment(o.match_cnt)
                        } else {
                            false
                        };
                }
                for i in 0..assign_cnt {
                    if !valid_assignment[i] {
                        continue;
                    }
                    let o = Self::select_overlap(k, &fragment_assignment[i]);
                    let seq_idx = usize::try_from(o.seq_idx).unwrap();
                    if ref_pos[i] < i32::try_from(self.base_variants[seq_idx].len()).unwrap() {
                        let p = usize::try_from(ref_pos[i]).unwrap();
                        let cid = self.base_variants[seq_idx][p].candidate_id;
                        if cid != -1 {
                            first_candidate_id = cid;
                            break;
                        }
                    }
                }

                if first_candidate_id != -1 {
                    // Contains candidate variants: register every other
                    // valid, currently-uncandidate aligned position (that is
                    // a match/mismatch op) as a new candidate.
                    for i in 0..assign_cnt {
                        if !valid_assignment[i] {
                            continue;
                        }
                        let o = Self::select_overlap(k, &fragment_assignment[i]);
                        let seq_idx = usize::try_from(o.seq_idx).unwrap();
                        let p = usize::try_from(ref_pos[i]).unwrap();
                        let op = o.align.as_ref().and_then(|a| a.get(align_idx[i])).copied();
                        if self.base_variants[seq_idx][p].candidate_id == -1
                            && matches!(
                                op,
                                Some(
                                    crate::align_algo::EDIT_MATCH
                                        | crate::align_algo::EDIT_MISMATCH
                                )
                            )
                        {
                            let cid = self.candidate_variants.len();
                            self.candidate_variants.push((o.seq_idx, ref_pos[i]));
                            self.base_variants[seq_idx][p].candidate_id =
                                i32::try_from(cid).unwrap();
                            adj_var_to_var.push(AdjVarToVar {
                                var_idx: i32::try_from(cid).unwrap(),
                                weight: 0.0,
                                root_candidate: false,
                                next: -1,
                            });
                            self.candidate_variant_group_id.push(-1);
                        }
                        let cid = self.base_variants[seq_idx][p].candidate_id;
                        if cid != -1 {
                            self.candidate_variant_group_id[usize::try_from(cid).unwrap()] = -1;
                        }
                    }

                    // Update the var-to-var abundance.
                    for i in 0..assign_cnt {
                        if !valid_assignment[i] {
                            continue;
                        }
                        for l in 0..assign_cnt {
                            if i == l || !valid_assignment[l] {
                                continue;
                            }
                            let o_i = Self::select_overlap(k, &fragment_assignment[i]);
                            let o_l = Self::select_overlap(k, &fragment_assignment[l]);
                            let seq_i = usize::try_from(o_i.seq_idx).unwrap();
                            let seq_l = usize::try_from(o_l.seq_idx).unwrap();
                            let p_i = usize::try_from(ref_pos[i]).unwrap();
                            let p_l = usize::try_from(ref_pos[l]).unwrap();
                            let cid_i = self.base_variants[seq_i][p_i].candidate_id;
                            let cid_l = self.base_variants[seq_l][p_l].candidate_id;
                            if cid_i == -1 || cid_l == -1 {
                                continue;
                            }
                            let cid_i_usize = usize::try_from(cid_i).unwrap();
                            let mut p = adj_var_to_var[cid_i_usize].next;
                            let mut found = false;
                            while p != -1 {
                                let p_usize = usize::try_from(p).unwrap();
                                if adj_var_to_var[p_usize].var_idx == cid_l {
                                    adj_var_to_var[p_usize].weight += 1.0;
                                    found = true;
                                    break;
                                }
                                p = adj_var_to_var[p_usize].next;
                            }
                            if !found {
                                let na = AdjVarToVar {
                                    var_idx: cid_l,
                                    weight: 1.0,
                                    root_candidate: false,
                                    next: adj_var_to_var[cid_i_usize].next,
                                };
                                adj_var_to_var[cid_i_usize].next =
                                    i32::try_from(adj_var_to_var.len()).unwrap();
                                adj_var_to_var.push(na);
                            }
                        }
                    }
                }

                // Move to the next read position.
                for i in 0..assign_cnt {
                    let o = Self::select_overlap(k, &fragment_assignment[i]);
                    let Some(align) = o.align.as_ref() else { continue };
                    while align_idx[i] < align.len() && read_pos[i] <= j {
                        let aidx = align_idx[i];
                        let op = align[aidx];
                        if op != crate::align_algo::EDIT_INSERT {
                            ref_pos[i] += 1;
                        }
                        if op != crate::align_algo::EDIT_DELETE {
                            read_pos[i] += 1;
                        }
                        align_idx[i] += 1;
                    }
                }
            }
        }
    }

    /// Ported from `VariantCaller::BuildCandidateVariantGroup`
    /// (`VariantCaller.hpp:573-593`): DFS over `adj_var_to_var`, tagging
    /// every candidate variant reachable from `from` with `tag` -- an edge
    /// `from -> to` is only followed if its co-occurrence `weight` is at
    /// least 15% of EITHER endpoint's total unweighted pileup depth (a noise
    /// filter).
    fn build_candidate_variant_group(
        &mut self,
        from: i32,
        tag: i32,
        adj_var_to_var: &[AdjVarToVar],
    ) {
        let from_usize = usize::try_from(from).unwrap();
        if self.candidate_variant_group_id[from_usize] != -1 {
            return;
        }
        self.candidate_variant_group_id[from_usize] = tag;
        let mut p = adj_var_to_var[from_usize].next;
        while p != -1 {
            let p_usize = usize::try_from(p).unwrap();
            let to = adj_var_to_var[p_usize].var_idx;
            let (from_seq, from_pos) = self.candidate_variants[from_usize];
            let (to_seq, to_pos) = self.candidate_variants[usize::try_from(to).unwrap()];
            let from_sum = self.base_variants[usize::try_from(from_seq).unwrap()]
                [usize::try_from(from_pos).unwrap()]
            .unweighted_count_sum();
            let to_sum = self.base_variants[usize::try_from(to_seq).unwrap()]
                [usize::try_from(to_pos).unwrap()]
            .unweighted_count_sum();
            let weight = adj_var_to_var[p_usize].weight;
            if weight >= from_sum * 0.15 || weight >= to_sum * 0.15 {
                self.build_candidate_variant_group(to, tag, adj_var_to_var);
            }
            p = adj_var_to_var[p_usize].next;
        }
    }

    /// Ported from `VariantCaller::BuildFragmentCandidateVarGraph`
    /// (`VariantCaller.hpp:595-687`): for each read end of each fragment
    /// assignment, walks the alignment and records an edge `(candidate
    /// variant) <-> (fragment, observed nucleotide)` for every aligned
    /// position that IS a known candidate variant -- deduplicating repeated
    /// `(fragIdx, nuc)` edges to the same candidate.
    fn build_fragment_candidate_var_graph(
        &self,
        read1: &[u8],
        read2: Option<&[u8]>,
        frag_idx: i32,
        fragment_assignment: &[FragmentOverlap],
        adj_frag: &mut Vec<AdjFragToVar>,
        adj_var: &mut Vec<AdjVarToFrag>,
    ) {
        let assign_cnt = fragment_assignment.len();
        if assign_cnt == 0 {
            return;
        }
        let frag_idx_usize = usize::try_from(frag_idx).unwrap();

        for k in 0..=1 {
            if k == 1 && !fragment_assignment[0].has_mate_pair {
                break;
            }
            let read: &[u8] = if k == 1 || (k == 0 && fragment_assignment[0].o1_from_r2) {
                read2.expect("k==1 or o1FromR2 implies read2 is present")
            } else {
                read1
            };
            let rc = reverse_complement(read);

            for frag in fragment_assignment {
                let o = Self::select_overlap(k, frag);
                let r: &[u8] = if o.strand == 1 { read } else { &rc };
                let seq_idx = usize::try_from(o.seq_idx).unwrap();
                let mut ref_pos = o.seq_start;
                let mut read_pos = o.read_start;
                let Some(align) = o.align.as_ref() else { continue };
                for &op in align {
                    // A trailing EDIT_INSERT can advance `ref_pos` one past
                    // the last consensus position (C++ then reads
                    // `baseVariants[seqIdx][refPos]` as unchecked UB that
                    // happens to not affect substitution counts; here we
                    // stop the walk instead of reading out of bounds, which
                    // is behavior-equivalent since a trailing insertion
                    // contributes no base variant).
                    //
                    // CONFIRMED REACHABLE on the real `fg-t1k analyze`
                    // pipeline (Task 6b's IMP-1 follow-up), not just a
                    // hypothetical: `align` here always comes from
                    // [`add_overlap_alignment_info`]'s
                    // [`crate::align_algo::global_alignment`] call over `t =
                    // allele_consensus[seq_start..=seq_end]` (`lent`) vs `p =
                    // read[read_start..=read_end]` (`lenp`) -- and
                    // `SeqSet::ExtendOverlap` (`SeqSet.hpp:1994-2100`, this
                    // port's `extend_overlap`) does not guarantee `lent ==
                    // lenp`: the two overhang extensions each advance BOTH
                    // spans by the SAME clamped amount, but the underlying
                    // k-mer-CHAINED overlap they extend from can already have
                    // unequal `seqEnd-seqStart` vs `readEnd-readStart` spans
                    // whenever a real indel separates two k-mer hits
                    // (`GetOverlapsFromHits`, unrelated to this guard, is
                    // free to chain hits across such a gap). A scratch probe
                    // (`fg_t1k_core::align_algo::global_alignment(b"ACGTACGTAC",
                    // b"ACGTACGTACGG", DEFAULT_BAND)`, `lent=10 < lenp=12`)
                    // confirms `global_alignment` DOES place the resulting
                    // gap as a trailing `EDIT_INSERT` in that shape -- so this
                    // guard is load-bearing on real inputs, not merely a
                    // defensive no-op for a path C++'s own trace makes
                    // unreachable.
                    if usize::try_from(ref_pos).unwrap() >= self.base_variants[seq_idx].len() {
                        break;
                    }
                    let p = usize::try_from(ref_pos).unwrap();
                    let cid = self.base_variants[seq_idx][p].candidate_id;
                    if cid != -1 {
                        let nuc = r[usize::try_from(read_pos).unwrap()];
                        let cid_usize = usize::try_from(cid).unwrap();

                        // Check whether the edge has already been added.
                        let mut p_edge = adj_var[cid_usize].next;
                        let mut found = false;
                        while p_edge != -1 {
                            let pe = usize::try_from(p_edge).unwrap();
                            if adj_var[pe].frag_idx == frag_idx && adj_var[pe].nuc == nuc {
                                found = true;
                                break;
                            }
                            p_edge = adj_var[pe].next;
                        }

                        if !found {
                            let n_frag_to_var =
                                AdjFragToVar { nuc, next: adj_frag[frag_idx_usize].next };
                            adj_frag[frag_idx_usize].next = i32::try_from(adj_frag.len()).unwrap();
                            adj_frag.push(n_frag_to_var);

                            let n_var_to_frag =
                                AdjVarToFrag { frag_idx, nuc, next: adj_var[cid_usize].next };
                            adj_var[cid_usize].next = i32::try_from(adj_var.len()).unwrap();
                            adj_var.push(n_var_to_frag);
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

    /// Ported from `VariantCaller::EnumerateVariants`
    /// (`VariantCaller.hpp:689-820`): exhaustively enumerates every `4^n`
    /// nucleotide assignment (`n = vars.len()`) for a variant group, scoring
    /// each by how many fragments it "covers" (mirrors the C++'s exact
    /// scoring, including the `varCnt <= 1` single-variant noise-tolerance
    /// special case), and records the best-covering assignment
    /// (`result.best_enum_variants`) plus any assignment TIED with it
    /// (`result.equal_best_enum_variants`, OVERWRITTEN each time a tie is
    /// found -- matches C++'s `result.equalBestEnumVariants = choices`,
    /// which keeps only the LAST tie seen, not all of them).
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn enumerate_variants(
        &self,
        depth: usize,
        choices: &mut Vec<u8>,
        result: &mut EnumVarResult,
        frag_ids: &[i32],
        vars: &[i32],
        adj_var: &[AdjVarToFrag],
        consensus: &[Vec<u8>],
    ) {
        if depth == vars.len() {
            let var_cnt = vars.len();
            let mut used_var_cnt = 0i32;
            let max_frag_idx = frag_ids.iter().copied().max().unwrap_or(-1) + 1;
            let max_frag_idx_usize = usize::try_from(max_frag_idx.max(0)).unwrap();
            let mut frag_covered = vec![0u8; max_frag_idx_usize];

            for i in 0..var_cnt {
                let (seq_idx, ref_pos) = self.var_pos(vars[i]);
                if var_cnt <= 1
                    && self.seq_copy[usize::try_from(seq_idx).unwrap()] <= 1
                    && choices[i]
                        != consensus[usize::try_from(seq_idx).unwrap()]
                            [usize::try_from(ref_pos).unwrap()]
                {
                    continue;
                }
                let var_idx_usize = usize::try_from(vars[i]).unwrap();
                let mut p = adj_var[var_idx_usize].next;
                while p != -1 {
                    let pe = usize::try_from(p).unwrap();
                    let edge_frag_idx = adj_var[pe].frag_idx;
                    if edge_frag_idx < max_frag_idx && adj_var[pe].nuc == choices[i] {
                        frag_covered[usize::try_from(edge_frag_idx).unwrap()] = 1;
                    }
                    p = adj_var[pe].next;
                }
            }

            // Single-variant noise-tolerance pass (VariantCaller.hpp:732-776):
            // only runs when varCnt <= 1.
            for i in 0..var_cnt {
                if var_cnt > 1 {
                    break;
                }
                let (seq_idx, ref_pos) = self.var_pos(vars[i]);
                let seq_idx_usize = usize::try_from(seq_idx).unwrap();
                if self.seq_copy[seq_idx_usize] != 1 {
                    continue;
                }
                let ref_nuc = consensus[seq_idx_usize][usize::try_from(ref_pos).unwrap()];
                if choices[i] == ref_nuc {
                    continue;
                }

                let mut ref_contribution = 0i32;
                let mut alt_contribution = 0i32;
                let var_idx_usize = usize::try_from(vars[i]).unwrap();
                let mut p = adj_var[var_idx_usize].next;
                while p != -1 {
                    let pe = usize::try_from(p).unwrap();
                    if adj_var[pe].nuc == choices[i] {
                        alt_contribution += 1;
                    } else if ref_nuc == adj_var[pe].nuc {
                        ref_contribution += 1;
                    }
                    p = adj_var[pe].next;
                }

                let uniq_alt = self.base_variants[seq_idx_usize][usize::try_from(ref_pos).unwrap()]
                    .uniq_count[nuc_to_num(choices[i]).unwrap()];
                let include_alt = ((alt_contribution >= 2 && uniq_alt > 0.0)
                    || alt_contribution >= 10)
                    && f64::from(alt_contribution) > 0.15 * f64::from(ref_contribution);

                let mut p = adj_var[var_idx_usize].next;
                while p != -1 {
                    let pe = usize::try_from(p).unwrap();
                    if ref_nuc == adj_var[pe].nuc || (choices[i] == adj_var[pe].nuc && include_alt)
                    {
                        let edge_frag_idx = adj_var[pe].frag_idx;
                        let fi = usize::try_from(edge_frag_idx).unwrap();
                        if fi < frag_covered.len() && frag_covered[fi] == 0 {
                            frag_covered[fi] = 2;
                        }
                    }
                    p = adj_var[pe].next;
                }
            }

            let mut covered = 0.0;
            for &fi in frag_ids {
                let fi_usize = usize::try_from(fi).unwrap();
                if fi_usize < frag_covered.len() && frag_covered[fi_usize] != 0 {
                    covered += 1.0;
                }
            }

            for i in 0..var_cnt {
                let (seq_idx, ref_pos) = self.var_pos(vars[i]);
                let ref_nuc =
                    consensus[usize::try_from(seq_idx).unwrap()][usize::try_from(ref_pos).unwrap()];
                if ref_nuc != choices[i] {
                    used_var_cnt += 1;
                }
            }

            #[allow(clippy::float_cmp)]
            let covered_is_better = covered > result.best_cover;
            #[allow(clippy::float_cmp)]
            let covered_ties_and_fewer_vars =
                covered == result.best_cover && used_var_cnt < result.used_var_cnt;
            if covered_is_better || covered_ties_and_fewer_vars {
                result.best_cover = covered;
                result.used_var_cnt = used_var_cnt;
                result.best_enum_variants = choices.clone();
                result.equal_best_enum_variants.clear();
            } else {
                #[allow(clippy::float_cmp)]
                let exact_tie = covered == result.best_cover && used_var_cnt == result.used_var_cnt;
                if exact_tie {
                    result.equal_best_enum_variants = choices.clone();
                }
            }
            return;
        }

        for &nuc in &NUM_TO_NUC {
            choices[depth] = nuc;
            self.enumerate_variants(depth + 1, choices, result, frag_ids, vars, adj_var, consensus);
        }
    }

    /// Small helper: `(seqIdx, refPos)` for `candidate_variants[var_idx]`.
    fn var_pos(&self, var_idx: i32) -> (i32, i32) {
        self.candidate_variants[usize::try_from(var_idx).unwrap()]
    }

    /// Ported from `VariantCaller::SolveVariantGroup`
    /// (`VariantCaller.hpp:822-976`): resolves one variant group (a list of
    /// candidate-variant indices believed to co-occur), producing zero or
    /// more [`Variant`]s. Skipped entirely (matching C++'s early `return`)
    /// when: the group exceeds `max_var_group_to_resolve` (and that limit is
    /// non-negative), any allele appears more than once in the group, or NO
    /// variant in the group falls in an exon (`ComputeVariant`/
    /// `SolveVariantGroup` only ever calls variants inside exons).
    /// `_adj_frag` mirrors `SolveVariantGroup`'s own `adjFrag` parameter --
    /// accepted for call-site symmetry with the C++ signature (`adjFrag,
    /// adjVar` are gathered together at every call site) but, like the C++,
    /// never actually read inside this function's body.
    #[allow(clippy::too_many_lines)]
    fn solve_variant_group(
        &mut self,
        vars: &[i32],
        _adj_frag: &[AdjFragToVar],
        adj_var: &[AdjVarToFrag],
        consensus: &[Vec<u8>],
    ) {
        let var_cnt = vars.len();
        let var_cnt_i32 = i32::try_from(var_cnt).expect("variant group size fits in i32");
        if var_cnt_i32 > self.max_var_group_to_resolve && self.max_var_group_to_resolve >= 0 {
            return;
        }

        let mut in_exon = false;
        let mut skip = false;
        let mut seq_idx_used: std::collections::HashMap<i32, i32> =
            std::collections::HashMap::new();
        for &v in vars {
            let (seq_idx, ref_pos) = self.var_pos(v);
            if self.base_variants[usize::try_from(seq_idx).unwrap()]
                [usize::try_from(ref_pos).unwrap()]
            .exon
            {
                in_exon = true;
            }
            let e = seq_idx_used.entry(seq_idx).or_insert(0);
            *e += 1;
            if *e > 1 {
                skip = true;
                break;
            }
        }
        if skip || !in_exon {
            return;
        }

        // Obtain related fragments (VariantCaller.hpp:868-900). The
        // commented-out `varCnt > 1` group-boundary check in C++ is dead
        // code (`fragP` is always `-1`) -- not ported, matching stock
        // behavior exactly.
        let mut frag_used: std::collections::HashSet<i32> = std::collections::HashSet::new();
        let mut frag_ids: Vec<i32> = Vec::new();
        for &v in vars {
            let v_usize = usize::try_from(v).unwrap();
            let mut p = adj_var[v_usize].next;
            while p != -1 {
                let pe = usize::try_from(p).unwrap();
                let edge_frag_idx = adj_var[pe].frag_idx;
                if frag_used.insert(edge_frag_idx) {
                    frag_ids.push(edge_frag_idx);
                }
                p = adj_var[pe].next;
            }
        }

        let mut choices = vec![0u8; var_cnt];
        let mut result = EnumVarResult {
            best_cover: -1.0,
            used_var_cnt: i32::try_from(var_cnt).unwrap() + 1,
            ..EnumVarResult::default()
        };
        self.enumerate_variants(0, &mut choices, &mut result, &frag_ids, vars, adj_var, consensus);

        let uniq = result.equal_best_enum_variants.is_empty();

        for (i, &v) in vars.iter().enumerate() {
            let (seq_idx, ref_pos) = self.var_pos(v);
            let seq_idx_usize = usize::try_from(seq_idx).unwrap();
            let ref_pos_usize = usize::try_from(ref_pos).unwrap();
            if !self.base_variants[seq_idx_usize][ref_pos_usize].exon {
                continue;
            }
            let ref_nuc = consensus[seq_idx_usize][ref_pos_usize];
            let var_nuc = result.best_enum_variants[i];
            if ref_nuc == var_nuc {
                continue;
            }

            let nv = Variant {
                seq_idx,
                ref_start: ref_pos,
                ref_end: ref_pos,
                reference: ref_nuc,
                var: var_nuc,
                all_support: self.base_variants[seq_idx_usize][ref_pos_usize].all_count_sum(),
                var_support: self.base_variants[seq_idx_usize][ref_pos_usize].count
                    [nuc_to_num(var_nuc).unwrap()],
                var_uniq_support: self.base_variants[seq_idx_usize][ref_pos_usize].uniq_count
                    [nuc_to_num(var_nuc).unwrap()],
                var_group_id: self.candidate_variant_group_id[usize::try_from(v).unwrap()],
                output_group_id: 0,
                qual: if uniq { 60 } else { 0 },
            };
            self.final_variants.push(nv);
        }

        if !uniq {
            for (i, &v) in vars.iter().enumerate() {
                let (seq_idx, ref_pos) = self.var_pos(v);
                let seq_idx_usize = usize::try_from(seq_idx).unwrap();
                let ref_pos_usize = usize::try_from(ref_pos).unwrap();
                if !self.base_variants[seq_idx_usize][ref_pos_usize].exon {
                    continue;
                }
                let ref_nuc = consensus[seq_idx_usize][ref_pos_usize];
                let var_nuc = result.equal_best_enum_variants[i];
                if ref_nuc == var_nuc {
                    continue;
                }

                let nv = Variant {
                    seq_idx,
                    ref_start: ref_pos,
                    ref_end: ref_pos,
                    reference: ref_nuc,
                    var: var_nuc,
                    all_support: self.base_variants[seq_idx_usize][ref_pos_usize].all_count_sum(),
                    var_support: self.base_variants[seq_idx_usize][ref_pos_usize].count
                        [nuc_to_num(var_nuc).unwrap()],
                    var_uniq_support: self.base_variants[seq_idx_usize][ref_pos_usize].uniq_count
                        [nuc_to_num(var_nuc).unwrap()],
                    var_group_id: self.candidate_variant_group_id[usize::try_from(v).unwrap()],
                    output_group_id: 1,
                    qual: 0,
                };
                self.final_variants.push(nv);
            }
        }
    }

    /// Ported from `VariantCaller::ComputeVariant`
    /// (`VariantCaller.hpp:978-1145`): the full driver -- pileup
    /// accumulation (two passes: `alignInfo`-only, then weighted counts,
    /// matching C++'s `updateType 1` then `0` ORDER exactly, since
    /// `IsGoodAssignment` in the SECOND pass reads `alignInfo` populated by
    /// the FIRST), candidate-variant discovery, iterative candidate
    /// expansion to a fixed point, variant-group construction, and
    /// per-group resolution via [`Self::solve_variant_group`].
    ///
    /// `allele_consensus[seqIdx]` must be the same reference consensus this
    /// [`VariantCaller`] was constructed from (`VariantCaller::new`'s
    /// `allele_refs`), in the same order.
    ///
    /// No-ops (matching C++'s early `return`) when
    /// `max_var_group_to_resolve == 0`.
    ///
    /// # Panics
    ///
    /// Panics if `fragment_assignments`/`allele_consensus` are inconsistent
    /// with how this [`VariantCaller`] was constructed (e.g. a
    /// `FragmentOverlap::seq_idx` out of range for `allele_consensus`, or an
    /// index/count that does not fit its target integer type) -- not
    /// expected for any well-formed input built from the same reference this
    /// [`VariantCaller`] was constructed from (mirrors the C++'s own
    /// undefined behavior on such inputs; this port fails loudly instead).
    #[allow(clippy::too_many_lines, clippy::needless_range_loop)]
    pub fn compute_variant(
        &mut self,
        read1: &[Vec<u8>],
        read2: &[Vec<u8>],
        fragment_assignments: &[Vec<FragmentOverlap>],
        allele_consensus: &[Vec<u8>],
    ) {
        if self.max_var_group_to_resolve == 0 {
            return;
        }

        let frag_cnt = fragment_assignments.len();
        let seq_cnt = self.base_variants.len();
        let has_mate = !read2.is_empty();

        // Identify the preliminary set of candidate variants (alignInfo-only
        // pass, updateType = 1).
        for i in 0..frag_cnt {
            let r2 = if has_mate { Some(read2[i].as_slice()) } else { None };
            self.update_base_variant_from_fragment_overlap(
                &read1[i],
                r2,
                1,
                &fragment_assignments[i],
                allele_consensus,
            );
        }
        // Weighted-count pass (updateType = 0).
        for i in 0..frag_cnt {
            let r2 = if has_mate { Some(read2[i].as_slice()) } else { None };
            self.update_base_variant_from_fragment_overlap(
                &read1[i],
                r2,
                0,
                &fragment_assignments[i],
                allele_consensus,
            );
        }

        self.find_candidate_variants(allele_consensus);
        let mut candidate_var_cnt = self.candidate_variants.len();

        let mut total_seq_len = 0usize;
        for i in 0..seq_cnt {
            total_seq_len += self.base_variants[i].len();
        }
        let mut adj_var_to_var: Vec<AdjVarToVar> =
            Vec::with_capacity(total_seq_len + 2 * candidate_var_cnt);
        adj_var_to_var.resize(
            total_seq_len,
            AdjVarToVar { var_idx: -1, weight: 0.0, root_candidate: false, next: -1 },
        );
        for i in 0..candidate_var_cnt {
            adj_var_to_var[i] = AdjVarToVar {
                var_idx: i32::try_from(i).unwrap(),
                weight: 0.0,
                root_candidate: true,
                next: -1,
            };
        }

        loop {
            let prev_candidate_var_cnt = self.candidate_variants.len();
            for entry in adj_var_to_var.iter_mut().take(prev_candidate_var_cnt) {
                entry.next = -1;
            }
            adj_var_to_var.resize(
                total_seq_len,
                AdjVarToVar { var_idx: -1, weight: 0.0, root_candidate: false, next: -1 },
            );

            let seq_candidate_var_accu_count: Vec<Vec<i32>> =
                (0..seq_cnt).map(|i| self.compute_candidate_var_accu_count(i)).collect();

            for i in 0..frag_cnt {
                let r2 = if has_mate { Some(read2[i].as_slice()) } else { None };
                self.expand_candidate_variants_from_fragment_overlap(
                    &read1[i],
                    r2,
                    &fragment_assignments[i],
                    &mut adj_var_to_var,
                    &seq_candidate_var_accu_count,
                );
            }
            if prev_candidate_var_cnt == self.candidate_variants.len() {
                break;
            }
        }

        candidate_var_cnt = self.candidate_variants.len();
        let mut group_cnt = 0i32;
        for i in 0..candidate_var_cnt {
            if adj_var_to_var[i].root_candidate && self.candidate_variant_group_id[i] == -1 {
                self.build_candidate_variant_group(
                    i32::try_from(i).unwrap(),
                    group_cnt,
                    &adj_var_to_var,
                );
                group_cnt += 1;
            }
        }

        let mut adj_frag: Vec<AdjFragToVar> = vec![AdjFragToVar { nuc: 0, next: -1 }; frag_cnt];
        let mut adj_var: Vec<AdjVarToFrag> =
            vec![AdjVarToFrag { frag_idx: -1, nuc: 0, next: -1 }; candidate_var_cnt];

        for i in 0..frag_cnt {
            let r2 = if has_mate { Some(read2[i].as_slice()) } else { None };
            #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
            self.build_fragment_candidate_var_graph(
                &read1[i],
                r2,
                i as i32,
                &fragment_assignments[i],
                &mut adj_frag,
                &mut adj_var,
            );
        }

        let mut candidate_var_group: Vec<Vec<i32>> =
            vec![Vec::new(); usize::try_from(group_cnt).unwrap()];
        for i in 0..candidate_var_cnt {
            let gid = self.candidate_variant_group_id[i];
            if gid == -1 {
                continue;
            }
            candidate_var_group[usize::try_from(gid).unwrap()].push(i32::try_from(i).unwrap());
        }

        for group in &candidate_var_group.clone() {
            self.solve_variant_group(group, &adj_frag, &adj_var, allele_consensus);
        }

        let final_var_cnt = self.final_variants.len();
        for i in 0..final_var_cnt {
            let seq_idx = usize::try_from(self.final_variants[i].seq_idx).unwrap();
            let ref_start = usize::try_from(self.final_variants[i].ref_start).unwrap();
            self.base_variants[seq_idx][ref_start].final_variant_ids.push(i);
        }
    }

    /// Read-only access to the called variants (`finalVariants`).
    #[must_use]
    pub fn final_variants(&self) -> &[Variant] {
        &self.final_variants
    }

    /// Ported from `VariantCaller::OutputAlleleVCF` (`VariantCaller.hpp:1202-1227`):
    /// formats every called variant as one VCF-ish text line. `names[seqIdx]`
    /// is the allele name (`refSet.GetSeqName`); `exonic_position(seqIdx,
    /// pos)` mirrors `refSet.GetExonicPosition` (the 0-based exon-relative
    /// coordinate reported as the VCF `POS`, 1-based per the C++ comment
    /// `// the VCF file is 1-based`).
    ///
    /// Field order (space-separated, matching the C++ `fprintf` format
    /// string EXACTLY, including the literal `.` placeholder columns and the
    /// trailing `refStart`/`outputGroupId` extra columns C++ appends beyond
    /// standard VCF):
    /// `name pos . ref var . filter varSupport allSupport varUniqSupport refStart outputGroupId`
    ///
    /// # Panics
    ///
    /// Panics if any called variant's `seq_idx` is out of range for `names`
    /// (not expected: every variant's `seq_idx` comes from this
    /// [`VariantCaller`]'s own reference-derived state).
    #[must_use]
    pub fn output_allele_vcf(
        &self,
        names: &[String],
        exonic_position: impl Fn(usize, i32) -> i32,
    ) -> String {
        use std::fmt::Write as _;
        let mut out = String::new();
        for variant in &self.final_variants {
            let filter = if variant.qual > 0 { "PASS" } else { "FAIL" };
            let seq_idx = usize::try_from(variant.seq_idx).unwrap();
            let exon_ref_start = exonic_position(seq_idx, variant.ref_start);
            writeln!(
                out,
                "{} {} . {} {} . {} {:.6} {:.6} {:.6} {} {}",
                names[seq_idx],
                exon_ref_start + 1,
                variant.reference as char,
                variant.var as char,
                filter,
                variant.var_support,
                variant.all_support,
                variant.var_uniq_support,
                variant.ref_start,
                variant.output_group_id,
            )
            .unwrap();
        }
        out
    }
}

#[cfg(test)]
// Every float assertion below compares exact, deterministic integer-valued
// pileup counts (e.g. `count[idx] == 1.0` after exactly one read's worth of
// weight-1.0 accumulation) -- exact equality is the CORRECT check here, not
// a bug clippy::float_cmp should flag (matches this crate's established
// convention of scoping this allow to individual exact-by-construction float
// comparisons, e.g. genotyper.rs's `#[allow(clippy::float_cmp)]` call sites).
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use crate::genotyper::AlleleRef;

    fn allele_ref(consensus: &str) -> AlleleRef {
        AlleleRef::new(consensus.as_bytes().to_vec(), None)
    }

    fn overlap_no_align(
        seq_idx: i32,
        seq_start: i32,
        seq_end: i32,
        read_start: i32,
        read_end: i32,
        match_cnt: i32,
        similarity: f64,
    ) -> Overlap {
        Overlap {
            seq_idx,
            read_start,
            read_end,
            seq_start,
            seq_end,
            strand: 1,
            match_cnt,
            similarity,
            align: None,
        }
    }

    // ---- IsGoodAssignment -------------------------------------------------

    #[test]
    fn is_good_assignment_true_when_no_align_info_yet() {
        let bv = BaseVariant::default();
        assert!(bv.is_good_assignment(0));
        assert!(bv.is_good_assignment(100));
    }

    #[test]
    fn is_good_assignment_false_when_match_cnt_far_below_best_seen() {
        let mut bv = BaseVariant::default();
        bv.align_info[0] = PairIntDouble { a: 100, b: 0.9 };
        // matchCnt = 90 < 100 - 4 = 96 -> false.
        assert!(!bv.is_good_assignment(90));
    }

    #[test]
    fn is_good_assignment_true_within_tolerance_of_best_seen() {
        let mut bv = BaseVariant::default();
        bv.align_info[0] = PairIntDouble { a: 100, b: 0.9 };
        // matchCnt = 96 == 100 - 4 -> not < threshold -> true.
        assert!(bv.is_good_assignment(96));
    }

    #[test]
    fn is_good_assignment_checks_all_four_nucleotide_slots() {
        let mut bv = BaseVariant::default();
        bv.align_info[3] = PairIntDouble { a: 50, b: 0.5 };
        assert!(!bv.is_good_assignment(10));
    }

    // ---- UpdateBaseVariantFromOverlap --------------------------------------

    #[test]
    fn update_base_variant_perfect_match_increments_ref_nucleotide() {
        let refs = vec![allele_ref("ACGT")];
        let mut vc = VariantCaller::new(&refs);
        let read = b"ACGT";
        let o = overlap_no_align(0, 0, 3, 0, 3, 8, 1.0);
        let consensus = [b"ACGT".to_vec()];
        vc.update_base_variant_from_overlap(read, 1.0, true, &o, &consensus[0]);
        // Every position should have its consensus base incremented once.
        for (pos, &base) in b"ACGT".iter().enumerate() {
            let idx = nuc_to_num(base).unwrap();
            assert_eq!(vc.base_variants[0][pos].count[idx], 1.0);
            assert_eq!(vc.base_variants[0][pos].uniq_count[idx], 1.0);
            assert_eq!(vc.base_variants[0][pos].unweighted_count[idx], 1.0);
        }
    }

    #[test]
    fn update_base_variant_mismatch_increments_alt_nucleotide() {
        let refs = vec![allele_ref("ACGT")];
        let mut vc = VariantCaller::new(&refs);
        let read = b"ATGT"; // position 1: C -> T mismatch
        let o = overlap_no_align(0, 0, 3, 0, 3, 6, 0.9);
        let consensus = [b"ACGT".to_vec()];
        vc.update_base_variant_from_overlap(read, 1.0, true, &o, &consensus[0]);
        let t_idx = nuc_to_num(b'T').unwrap();
        assert_eq!(vc.base_variants[0][1].count[t_idx], 1.0);
        let c_idx = nuc_to_num(b'C').unwrap();
        assert_eq!(vc.base_variants[0][1].count[c_idx], 0.0);
    }

    #[test]
    fn update_base_variant_weight_below_one_does_not_bump_uniq_count() {
        let refs = vec![allele_ref("ACGT")];
        let mut vc = VariantCaller::new(&refs);
        let read = b"ACGT";
        let o = overlap_no_align(0, 0, 3, 0, 3, 8, 1.0);
        let consensus = [b"ACGT".to_vec()];
        vc.update_base_variant_from_overlap(read, 0.5, true, &o, &consensus[0]);
        let a_idx = nuc_to_num(b'A').unwrap();
        assert_eq!(vc.base_variants[0][0].count[a_idx], 1.0);
        assert_eq!(vc.base_variants[0][0].uniq_count[a_idx], 0.0);
    }

    #[test]
    fn update_base_variant_n_base_is_skipped() {
        let refs = vec![allele_ref("ACGT")];
        let mut vc = VariantCaller::new(&refs);
        let read = b"ANGT";
        let o = overlap_no_align(0, 0, 3, 0, 3, 6, 0.9);
        let consensus = [b"ACGT".to_vec()];
        vc.update_base_variant_from_overlap(read, 1.0, true, &o, &consensus[0]);
        assert_eq!(vc.base_variants[0][1].all_count_sum(), 0.0);
    }

    #[test]
    fn update_base_variant_no_seq_is_noop() {
        let refs = vec![allele_ref("ACGT")];
        let mut vc = VariantCaller::new(&refs);
        let read = b"ACGT";
        let mut o = overlap_no_align(0, 0, 3, 0, 3, 8, 1.0);
        o.seq_idx = -1;
        let consensus = [b"ACGT".to_vec()];
        vc.update_base_variant_from_overlap(read, 1.0, true, &o, &consensus[0]);
        assert_eq!(vc.base_variants[0][0].all_count_sum(), 0.0);
    }

    #[test]
    fn update_base_variant_filter_low_qual_skips_below_threshold() {
        let refs = vec![allele_ref("ACGT")];
        let mut vc = VariantCaller::new(&refs);
        // Seed alignInfo with a high matchCnt at position 1.
        vc.base_variants[0][1].align_info[nuc_to_num(b'C').unwrap()] =
            PairIntDouble { a: 100, b: 1.0 };
        let read = b"ATGT";
        let o = overlap_no_align(0, 0, 3, 0, 3, 10, 0.5); // matchCnt=10 << 100-4
        let consensus = [b"ACGT".to_vec()];
        vc.update_base_variant_from_overlap(read, 1.0, true, &o, &consensus[0]);
        // Position 1 should be entirely skipped (filtered as low-qual).
        assert_eq!(vc.base_variants[0][1].all_count_sum(), 0.0);
    }

    // ---- end-to-end ComputeVariant regression scenarios --------------------

    /// Builds a single-mate, ungapped, whole-read `FragmentOverlap` against
    /// `seq_idx`, with `align` PRE-POPULATED (all `EDIT_MATCH`, since these
    /// scenarios only ever construct reads that are ungapped vs. the
    /// consensus at every position) -- matching the real end-to-end
    /// pipeline, where `Analyzer::AddFragmentAlignmentInfo` always populates
    /// `overlap.align` before `ComputeVariant` runs (see this module's doc
    /// comment on `align == NULL` being dead code on that path). Both
    /// [`VariantCaller::expand_candidate_variants_from_fragment_overlap`] and
    /// [`VariantCaller::build_fragment_candidate_var_graph`] require a
    /// populated `align` to build the fragment/variant graph that
    /// [`VariantCaller::enumerate_variants`] needs to actually call (not
    /// just detect a candidate for) a variant.
    fn frag_overlap(seq_idx: i32, read: &[u8]) -> FragmentOverlap {
        let len = i32::try_from(read.len()).unwrap();
        let align = vec![crate::align_algo::EDIT_MATCH; read.len()];
        FragmentOverlap {
            seq_idx,
            has_mate_pair: false,
            o1_from_r2: false,
            overlap1: Overlap {
                seq_idx,
                read_start: 0,
                read_end: len - 1,
                seq_start: 0,
                seq_end: len - 1,
                strand: 1,
                match_cnt: 2 * len,
                similarity: 1.0,
                align: Some(align),
            },
            overlap2: Overlap::none(),
        }
    }

    #[test]
    fn compute_variant_clean_snv_is_called() {
        // 10 reads support a T at position 5 (ref C), overwhelming majority.
        let consensus = "AAAAACAAAAAAAAAAAAAA";
        let refs = vec![allele_ref(consensus)];
        let mut vc = VariantCaller::new(&refs);
        vc.seq_abundance = vec![1.0];
        vc.seq_copy = vec![1];

        let mut read1 = Vec::new();
        let mut assigns = Vec::new();
        let alt = "AAAAATAAAAAAAAAAAAAA";
        for _ in 0..10 {
            read1.push(alt.as_bytes().to_vec());
            assigns.push(vec![frag_overlap(0, alt.as_bytes())]);
        }
        let consensus_vecs = vec![consensus.as_bytes().to_vec()];
        vc.compute_variant(&read1, &[], &assigns, &consensus_vecs);

        let variants = vc.final_variants();
        assert_eq!(variants.len(), 1, "expected exactly one called variant: {variants:?}");
        let v = &variants[0];
        assert_eq!(v.ref_start, 5);
        assert_eq!(v.reference, b'C');
        assert_eq!(v.var, b'T');
        assert_eq!(v.qual, 60);
        assert_eq!(v.output_group_id, 0);
    }

    #[test]
    fn compute_variant_ref_only_reads_call_nothing() {
        let consensus = "AAAAACAAAAAAAAAAAAAA";
        let refs = vec![allele_ref(consensus)];
        let mut vc = VariantCaller::new(&refs);
        vc.seq_abundance = vec![1.0];
        vc.seq_copy = vec![1];

        let mut read1 = Vec::new();
        let mut assigns = Vec::new();
        for _ in 0..10 {
            read1.push(consensus.as_bytes().to_vec());
            assigns.push(vec![frag_overlap(0, consensus.as_bytes())]);
        }
        let consensus_vecs = vec![consensus.as_bytes().to_vec()];
        vc.compute_variant(&read1, &[], &assigns, &consensus_vecs);
        assert!(vc.final_variants().is_empty());
    }

    #[test]
    fn compute_variant_low_support_is_not_called() {
        // Only 2 alt reads out of 10 -- below both the absolute (5) and
        // relative (0.5x ref) thresholds.
        let consensus = "AAAAACAAAAAAAAAAAAAA";
        let refs = vec![allele_ref(consensus)];
        let mut vc = VariantCaller::new(&refs);
        vc.seq_abundance = vec![1.0];
        vc.seq_copy = vec![1];

        let mut read1 = Vec::new();
        let mut assigns = Vec::new();
        let alt = "AAAAATAAAAAAAAAAAAAA";
        for i in 0..10 {
            let seq = if i < 2 { alt } else { consensus };
            read1.push(seq.as_bytes().to_vec());
            assigns.push(vec![frag_overlap(0, seq.as_bytes())]);
        }
        let consensus_vecs = vec![consensus.as_bytes().to_vec()];
        vc.compute_variant(&read1, &[], &assigns, &consensus_vecs);
        assert!(vc.final_variants().is_empty());
    }

    #[test]
    fn compute_variant_indel_read_does_not_spuriously_call_a_variant() {
        // A read with a 1bp deletion relative to consensus should not create
        // a false SNV call at the deletion site or drift the pileup
        // position (matches stock: EDIT_INSERT/EDIT_DELETE are never
        // pileup-counted, only ref/read position advancement -- see module
        // docs "No indel calling").
        let consensus = "AAAAACAAAAAAAAAAAAAA"; // len 20
        let refs = vec![allele_ref(consensus)];
        let mut vc = VariantCaller::new(&refs);
        vc.seq_abundance = vec![1.0];
        vc.seq_copy = vec![1];

        // read = consensus with position 5 ('C') deleted -> len 19.
        let read = "AAAAAAAAAAAAAAAAAAA";
        let align_ops = {
            use crate::align_algo::{EDIT_DELETE, EDIT_MATCH};
            let mut ops = vec![EDIT_MATCH; 5];
            ops.push(EDIT_DELETE);
            ops.extend(std::iter::repeat_n(EDIT_MATCH, 14));
            ops
        };
        let mut read1 = Vec::new();
        let mut assigns = Vec::new();
        for _ in 0..10 {
            read1.push(read.as_bytes().to_vec());
            let o = Overlap {
                seq_idx: 0,
                read_start: 0,
                read_end: 18,
                seq_start: 0,
                seq_end: 19,
                strand: 1,
                match_cnt: 38,
                similarity: 0.95,
                align: Some(align_ops.clone()),
            };
            assigns.push(vec![FragmentOverlap {
                seq_idx: 0,
                has_mate_pair: false,
                o1_from_r2: false,
                overlap1: o,
                overlap2: Overlap::none(),
            }]);
        }
        let consensus_vecs = vec![consensus.as_bytes().to_vec()];
        vc.compute_variant(&read1, &[], &assigns, &consensus_vecs);
        assert!(
            vc.final_variants().is_empty(),
            "indel-only divergence must not be called as a SNV: {:?}",
            vc.final_variants()
        );
    }

    #[test]
    fn compute_variant_trailing_insert_does_not_panic() {
        // Regression for an out-of-bounds panic in
        // `build_fragment_candidate_var_graph` (VariantCaller.hpp:639's
        // `baseVariants[seqIdx][refPos]` read at the TOP of the align-op
        // loop, BEFORE `refPos` is advanced): when an overlap's alignment
        // consumes the full ref range (an internal 1bp deletion compensates
        // so total length still matches) and then ends in a TRAILING
        // EDIT_INSERT, `refPos` has already reached `consensus_len` by the
        // time that final op is processed, so the unguarded read indexes one
        // past the end of `base_variants[seqIdx]`.
        //
        // C++ performs the identical unchecked read as silent UB and
        // survives without corrupting the call (see the outcome assertion
        // below for why): the garbage `candidateId` it may read from beyond
        // the buffer does not happen to overturn the otherwise-unambiguous
        // SNV vote. This port must not panic on the same input, and must
        // reach the same called-variant outcome.
        let consensus = "AAAAACAAAAAAAAAAAAAA"; // len 20
        let refs = vec![allele_ref(consensus)];
        let mut vc = VariantCaller::new(&refs);
        vc.seq_abundance = vec![1.0];
        vc.seq_copy = vec![1];

        // A clean SNV at position 5 on most reads, so `find_candidate_variants`
        // registers a real candidate and the align-walk in
        // `build_fragment_candidate_var_graph` has a live `candidate_id` to
        // check against (matching the reviewer's repro shape, not just an
        // empty-candidate no-op walk).
        let alt = "AAAAATAAAAAAAAAAAAAA"; // len 20, position 5: C -> T
        let mut read1 = Vec::new();
        let mut assigns = Vec::new();
        for _ in 0..9 {
            read1.push(alt.as_bytes().to_vec());
            assigns.push(vec![frag_overlap(0, alt.as_bytes())]);
        }

        // The panic-triggering fragment: consensus with position 5 ('C')
        // deleted (19 ref-consuming ops covering all 20 ref positions via 19
        // matches + 1 internal delete), followed by a TRAILING EDIT_INSERT
        // that pushes `refPos` from 20 to 21 without a corresponding ref
        // base -- i.e. the alignment ends mid-insertion, exactly the shape
        // that overruns `base_variants[seqIdx]` (len 20) in the unguarded
        // C++-mirroring read.
        let read = "AAAAAAAAAAAAAAAAAAAG"; // len 20: 19 ref-matched bases + 1 trailing inserted base
        let align_ops = {
            use crate::align_algo::{EDIT_DELETE, EDIT_INSERT, EDIT_MATCH};
            let mut ops = vec![EDIT_MATCH; 5];
            ops.push(EDIT_DELETE); // consumes ref pos 5 ('C'), no read base
            ops.extend(std::iter::repeat_n(EDIT_MATCH, 14)); // ref 6..=19
            ops.push(EDIT_INSERT); // trailing insert: refPos already at 20
            ops
        };
        let o = Overlap {
            seq_idx: 0,
            read_start: 0,
            read_end: 19,
            seq_start: 0,
            seq_end: 19,
            strand: 1,
            match_cnt: 38,
            similarity: 0.95,
            align: Some(align_ops),
        };
        read1.push(read.as_bytes().to_vec());
        assigns.push(vec![FragmentOverlap {
            seq_idx: 0,
            has_mate_pair: false,
            o1_from_r2: false,
            overlap1: o,
            overlap2: Overlap::none(),
        }]);

        let consensus_vecs = vec![consensus.as_bytes().to_vec()];
        // Must not panic (the regression under test). The trailing-insert
        // fragment votes 'A' at position 5 (from its in-bounds EDIT_DELETE
        // step, processed BEFORE the out-of-bounds trailing insert is ever
        // reached) while the other 9 fragments vote 'T' -- 'T' still wins
        // uniquely, so the SNV call is unaffected by the guard, exactly
        // matching what C++ does too (its unguarded OOB read only happens
        // on the trailing insert, which contributes nothing either way).
        vc.compute_variant(&read1, &[], &assigns, &consensus_vecs);
        let variants = vc.final_variants();
        assert_eq!(variants.len(), 1, "expected exactly one called variant: {variants:?}");
        let v = &variants[0];
        assert_eq!(v.ref_start, 5);
        assert_eq!(v.reference, b'C');
        assert_eq!(v.var, b'T');
        assert_eq!(v.qual, 60);
    }

    // ---- exonic_position ---------------------------------------------------

    #[test]
    fn exonic_position_single_whole_sequence_exon_is_identity() {
        let exons = [(0, 19)];
        assert_eq!(exonic_position(&exons, 0), 0);
        assert_eq!(exonic_position(&exons, 19), 19);
        assert_eq!(exonic_position(&exons, 10), 10);
    }

    #[test]
    fn exonic_position_maps_across_multiple_exons() {
        // Exon 0: genomic [0,9] -> exonic [0,9]; exon 1: genomic [20,29] ->
        // exonic [10,19] (the intervening intron [10,19] is excluded).
        let exons = [(0, 9), (20, 29)];
        assert_eq!(exonic_position(&exons, 0), 0);
        assert_eq!(exonic_position(&exons, 9), 9);
        assert_eq!(exonic_position(&exons, 20), 10);
        assert_eq!(exonic_position(&exons, 25), 15);
        assert_eq!(exonic_position(&exons, 29), 19);
    }

    #[test]
    fn exonic_position_returns_negative_one_outside_any_exon() {
        let exons = [(0, 9), (20, 29)];
        assert_eq!(exonic_position(&exons, 10), -1, "intron position must return -1");
        assert_eq!(exonic_position(&exons, 15), -1);
    }

    // ---- add_overlap_alignment_info / add_fragment_alignment_info ---------

    #[test]
    fn add_overlap_alignment_info_populates_align_for_perfect_match() {
        let consensus = b"ACGTACGTAC";
        let read = b"ACGTACGTAC";
        let mut o = Overlap {
            seq_idx: 0,
            read_start: 0,
            read_end: 9,
            seq_start: 0,
            seq_end: 9,
            strand: 1,
            match_cnt: 20,
            similarity: 1.0,
            align: None,
        };
        add_overlap_alignment_info(read, &mut o, consensus);
        let align = o.align.expect("align must be populated");
        assert_eq!(align, vec![crate::align_algo::EDIT_MATCH; 10]);
    }

    #[test]
    fn add_overlap_alignment_info_reverse_strand_uses_reverse_complement() {
        // consensus and its reverse complement: a perfect match on strand -1
        // must reverse-complement the read before aligning.
        let consensus = b"AAAACCCC"; // len 8
        let read = b"GGGGTTTT"; // reverse complement is AAAACCCC
        let mut o = Overlap {
            seq_idx: 0,
            read_start: 0,
            read_end: 7,
            seq_start: 0,
            seq_end: 7,
            strand: -1,
            match_cnt: 16,
            similarity: 1.0,
            align: None,
        };
        add_overlap_alignment_info(read, &mut o, consensus);
        let align = o.align.expect("align must be populated");
        assert_eq!(
            align,
            vec![crate::align_algo::EDIT_MATCH; 8],
            "RC(read) must match consensus exactly"
        );
    }

    #[test]
    fn add_overlap_alignment_info_no_overlap_is_noop() {
        let consensus = b"ACGT";
        let read = b"ACGT";
        let mut o = Overlap::none();
        assert!(o.align.is_none());
        add_overlap_alignment_info(read, &mut o, consensus);
        assert!(o.align.is_none(), "seq_idx == -1 must leave align untouched (None)");
    }

    #[test]
    fn add_overlap_alignment_info_overwrites_existing_align_unconditionally() {
        // Unlike VariantCaller::update_base_variant_from_overlap's `align ==
        // NULL` fallback, AddOverlapAlignmentInfo always recomputes -- even
        // if `align` was already populated with something else.
        let consensus = b"ACGT";
        let read = b"ACGT";
        let mut o = Overlap {
            seq_idx: 0,
            read_start: 0,
            read_end: 3,
            seq_start: 0,
            seq_end: 3,
            strand: 1,
            match_cnt: 8,
            similarity: 1.0,
            align: Some(vec![crate::align_algo::EDIT_MISMATCH; 4]), // stale/wrong on purpose
        };
        add_overlap_alignment_info(read, &mut o, consensus);
        assert_eq!(
            o.align,
            Some(vec![crate::align_algo::EDIT_MATCH; 4]),
            "must overwrite, not keep stale align"
        );
    }

    #[test]
    fn add_fragment_alignment_info_populates_both_mates_for_mate_pair_fragment() {
        let consensus = vec![b"ACGTACGTAC".to_vec()];
        let read1 = b"ACGTA";
        let read2 = b"CGTAC";
        let mut assignment = vec![FragmentOverlap {
            seq_idx: 0,
            has_mate_pair: true,
            o1_from_r2: false,
            overlap1: Overlap {
                seq_idx: 0,
                read_start: 0,
                read_end: 4,
                seq_start: 0,
                seq_end: 4,
                strand: 1,
                match_cnt: 10,
                similarity: 1.0,
                align: None,
            },
            overlap2: Overlap {
                seq_idx: 0,
                read_start: 0,
                read_end: 4,
                seq_start: 5,
                seq_end: 9,
                strand: 1,
                match_cnt: 10,
                similarity: 1.0,
                align: None,
            },
        }];
        add_fragment_alignment_info(read1, Some(read2), &mut assignment, &consensus);
        assert!(assignment[0].overlap1.align.is_some());
        assert!(assignment[0].overlap2.align.is_some());
    }

    #[test]
    fn add_fragment_alignment_info_single_end_populates_overlap1_only() {
        let consensus = vec![b"ACGTACGTAC".to_vec()];
        let read1 = b"ACGTA";
        let mut assignment = vec![FragmentOverlap {
            seq_idx: 0,
            has_mate_pair: false,
            o1_from_r2: false,
            overlap1: Overlap {
                seq_idx: 0,
                read_start: 0,
                read_end: 4,
                seq_start: 0,
                seq_end: 4,
                strand: 1,
                match_cnt: 10,
                similarity: 1.0,
                align: None,
            },
            overlap2: Overlap::none(),
        }];
        add_fragment_alignment_info(read1, None, &mut assignment, &consensus);
        assert!(assignment[0].overlap1.align.is_some());
        assert!(assignment[0].overlap2.align.is_none(), "overlap2 is unused (no mate pair)");
    }

    #[test]
    fn add_fragment_alignment_info_o1_from_r2_aligns_against_read2() {
        let consensus = vec![b"ACGTACGTAC".to_vec()];
        let read1 = b"NNNNN"; // must not be touched
        let read2 = b"ACGTA";
        let mut assignment = vec![FragmentOverlap {
            seq_idx: 0,
            has_mate_pair: false,
            o1_from_r2: true,
            overlap1: Overlap {
                seq_idx: 0,
                read_start: 0,
                read_end: 4,
                seq_start: 0,
                seq_end: 4,
                strand: 1,
                match_cnt: 10,
                similarity: 1.0,
                align: None,
            },
            overlap2: Overlap::none(),
        }];
        add_fragment_alignment_info(read1, Some(read2), &mut assignment, &consensus);
        let align = assignment[0].overlap1.align.as_ref().expect("align must be populated");
        assert_eq!(
            align,
            &vec![crate::align_algo::EDIT_MATCH; 5],
            "must align against read2, not read1"
        );
    }
}
