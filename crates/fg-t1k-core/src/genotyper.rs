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
//!   weight assignment that Phase 5b's EM will consume.
//!
//! # Deliberately NOT ported here (deferred to 5b/5c)
//!
//! - `EMupdate`/`Quantify`/`SQUAREMalpha`/`CoalesceReadAssignments`/
//!   `FinalizeReadAssignments`/`BuildAlleleEquivalentClass` (the EM loop and
//!   its equivalent-class bookkeeping) -- Phase 5b.
//! - `SelectAllele`/genotype-quality scoring/output formatting -- Phase 5c.
//! - `ReadAssignmentToFragmentAssignment` (`SeqSet.hpp:2310-2556`, builds
//!   `_fragmentOverlap`s from per-mate `_overlap` lists) -- this belongs to
//!   the read-processing loop that calls `SetReadAssignments` in a loop
//!   (`Genotyper.cpp:160-192`); [`Genotyper::set_read_assignments`] here
//!   takes already-built [`FragmentOverlap`]s as input, matching the C++
//!   method signature exactly, so this port is agnostic to how its caller
//!   produces them.
//! - `_readGroupInfo`'s actual population (the `readCnt`-length vector built
//!   inside `Quantify`, `Genotyper.hpp:1151-1163`) -- `Quantify` is Phase
//!   5b's territory; [`ReadGroupInfo`] is defined here (it is trivial: a
//!   single `count: f64` field) purely so 5b's `EMupdate` port can reuse it
//!   without a follow-up struct addition.

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
        }
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
    /// `readsInAllele`/`readAssignments`/`readAssignmentsFingerprintToIdx`
    /// are Phase 5b's `CoalesceReadAssignments`/`FinalizeReadAssignments`
    /// territory (they are cleared here in the C++ but never touched by
    /// `SetReadAssignments`, only by those later methods) -- this port
    /// omits them from [`Genotyper`]'s fields entirely rather than storing
    /// always-empty placeholders; Phase 5b should add them alongside its
    /// own port of those methods.
    pub fn init_read_assignments(&mut self, total_read_cnt: i32, max_assign_cnt: i32) {
        self.max_assign_cnt = max_assign_cnt;
        self.read_cnt = 0;
        self.total_read_cnt = total_read_cnt;
        self.all_read_assignments = vec![Vec::new(); usize::try_from(total_read_cnt).unwrap_or(0)];
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
}
