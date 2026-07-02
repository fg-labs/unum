//! Reference-k-mer read-candidate filter, ported from a slice of T1K's
//! `SeqSet` class (`vendor/t1k/SeqSet.hpp`) plus the free functions
//! `IsLowComplexity`/`IsGoodCandidate` (`vendor/t1k/FastqExtractor.cpp:89-118`,
//! byte-identical to `vendor/t1k/BamExtractor.cpp:144-166,214`).
//!
//! This is the *read-candidate-filtering* slice only -- the part of `SeqSet`
//! that `FastqExtractor`/`BamExtractor` use to decide whether a read is
//! worth keeping for downstream genotyping: reference load + index build
//! (`SeqSet::InputRefFa`), k-mer hit collection (`SeqSet::GetHitsFromRead`),
//! and BOTH gates of `SeqSet::HasHitInSet` -- the bucket-count gate AND the
//! `GetOverlapsFromHits`-based (LIS hit-chaining) mismatch-threshold
//! confirmation (ported in [`crate::overlap`]). It does **not** port
//! `AlignAlgo` (banded Smith-Waterman alignment) -- `GetOverlapsFromHits`,
//! as called from `HasHitInSet`, never invokes it (see [`crate::overlap`]'s
//! module docs); real alignment-based scoring is a separate, later-phase
//! concern unrelated to this gate.
//!
//! # `kmerLength` is a caller-supplied parameter here, not a fixed constant
//!
//! Stock T1K does not use a single hardcoded default k-mer length. In
//! `FastqExtractor.cpp:main` (`FastqExtractor.cpp:272-418`):
//! 1. `SeqSet refSet(9)` is constructed with a literal initial `kmerLength = 9`.
//! 2. The reference is loaded at that k-mer length (`refSet.InputRefFa(...)`,
//!    `FastqExtractor.cpp:302`).
//! 3. `SeqSet::InferKmerLength()` (`SeqSet.hpp:2830-2845`) computes a
//!    *data-dependent* k-mer length from the total loaded reference length
//!    (`totalLength`): `kmerLength = floor(log4(totalLength)) + 2` (the loop
//!    counts base-4 "digits" of `totalLength`, then adds one more).
//! 4. If that inferred value is *strictly greater* than the initial `9`, the
//!    `SeqSet` is rebuilt at the new k-mer length via `UpdateKmerLength`
//!    (`FastqExtractor.cpp:411-418`, `SeqSet.hpp:2847-2857`); otherwise it
//!    stays at `9`.
//!
//! [`RefKmerFilter::from_reference_fasta`] intentionally does not replicate
//! this two-stage default-then-infer-then-maybe-rebuild dance -- it takes
//! `kmer_length` as an explicit parameter, matching a single plain
//! `SeqSet(kmerLength)` construction followed by one `InputRefFa` call
//! (`SeqSet.hpp:760-772,872-904`, no `InferKmerLength`/`UpdateKmerLength`
//! call). `InferKmerLength`/`UpdateKmerLength` are not ported here; a future
//! caller (a `fghla`-style CLI/pipeline layer) that wants stock's exact
//! default-selection behavior can compute it itself before calling
//! `from_reference_fasta`. For `fixtures/refbuild/golden/kir_rna_seq.fa`
//! specifically (the fixture used by this module's differential tests),
//! `InferKmerLength()` evaluates to `8` for that fixture's total reference
//! length (8781 bases) -- NOT greater than the initial `9` -- so stock T1K
//! would actually use `kmerLength = 9` unmodified for this exact fixture;
//! the differential tests use `kmer_length = 9` to match that real-world
//! value precisely.
//!
//! # `HasHitInSet`'s two gates, both now ported
//!
//! Stock `SeqSet::HasHitInSet` (`SeqSet.hpp:1915-1990`) does two things in
//! sequence:
//! 1. Collect hits from both strands (`GetHitsFromRead`), bucket them by
//!    `(strand-tag, seqIdx)`, and find the single largest bucket. If
//!    `kmerLength * <largest bucket size> < hitLenRequired`, return `false`
//!    immediately (`SeqSet.hpp:1929-1964`).
//! 2. Otherwise, run `GetOverlapsFromHits` (LIS-based colinear hit chaining;
//!    NOT `AlignAlgo`-based -- see [`crate::overlap`]'s module docs) on the
//!    winning bucket, and only return `true` if at least one resulting
//!    overlap's implied mismatch count is within `int(len *
//!    (1 - refSeqSimilarity)) * kmerLength` (`SeqSet.hpp:1966-1983`).
//!
//! [`RefKmerFilter::has_hit_in_set`] implements BOTH gates: step 1 (the
//! bucket-count gate) followed by step 2 (calling
//! [`crate::overlap::get_overlaps_from_hits`] on the winning bucket's hits,
//! reconstructed in original insertion order -- see that function's doc
//! comment for why a stable filter of `scratch.hits` is exactly equivalent
//! to stock's per-hit `buckets[tag][idx].PushBack(...)` construction). This
//! makes `has_hit_in_set` byte/value-identical to stock's `HasHitInSet` (not
//! merely a superset), confirmed by an adversarial differential
//! (`crates/fg-t1k-sys/tests/diff_refkmerfilter.rs`) against the real C++
//! oracle across curated AND arbitrary/random/noncolinear/tie-inducing reads.
//!
//! # Reused Phase-2 primitives
//!
//! Reuses [`crate::kmer::KmerCode`] (rolling k-mer encoder) and
//! [`crate::kmer_index::KmerIndex`] (forward-code-keyed position index)
//! unmodified. `KmerIndex::search`'s forward-not-canonical keying and
//! `KmerIndex::build_index_from_read`'s consecutive-duplicate dedup are
//! exactly what `SeqSet::InputRefFa`'s `seqIndex.BuildIndexFromRead(...)`
//! call needs (`SeqSet.hpp:902`); see `kmer_index`'s module docs for that
//! dedup's exact semantics (already covered by Phase-2's differential
//! tests, not re-derived here). Also reuses [`crate::overlap`]'s
//! `get_overlaps_from_hits`/LIS port for `has_hit_in_set`'s gate 2.

use crate::kmer::KmerCode;
use crate::kmer_index::KmerIndex;
use crate::overlap::{self, OverlapHit};
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs;
use std::path::Path;

/// `SeqSet`'s default `hitLenRequired` (`SeqSet.hpp:764`), used by
/// [`RefKmerFilter::from_reference_fasta`] since it mirrors a plain
/// `SeqSet(kmerLength)` construction with no `SetHitLenRequired` call (see
/// module docs: `_load_ref`/`from_reference_fasta` do not replicate
/// `FastqExtractor.cpp`'s dynamic `SetHitLenRequired` call, which requires
/// sampling actual read lengths from a FASTQ file upstream of reference
/// loading).
const DEFAULT_HIT_LEN_REQUIRED: i32 = 31;

/// `SeqSet`'s default `radius` (`SeqSet.hpp:763`), used by
/// [`RefKmerFilter::has_hit_in_set`]'s gate 2 (`overlap::get_overlaps_from_hits`'s
/// `radius` parameter). Fixed at the `SeqSet(kmerLength)` constructor
/// default for the same reason as [`DEFAULT_HIT_LEN_REQUIRED`]: no
/// `SetRadius` call exists on the plain-construction path this port models.
const DEFAULT_RADIUS: i32 = 10;

/// `SeqSet`'s default `refSeqSimilarity` (`SeqSet.hpp:768`), used by
/// [`RefKmerFilter::has_hit_in_set`]'s gate 2 mismatch-threshold computation
/// (`SeqSet.hpp:1973`). Fixed at the `SeqSet(kmerLength)` constructor
/// default; see [`DEFAULT_HIT_LEN_REQUIRED`]'s doc comment for why this port
/// does not replicate any dynamic `Set*` call.
const DEFAULT_REF_SEQ_SIMILARITY: f64 = 0.8;

/// `SeqSet::isLongSeqSet`'s default (`SeqSet.hpp:765`). Always `false` on the
/// plain-construction path this port models (no public setter reachable from
/// `SeqSet(kmerLength)`; see `get_hits_from_read`'s doc comment for the same
/// point made about `downSample`). Passed straight through to
/// `overlap::get_overlaps_from_hits`'s `is_long_seq_set` parameter.
const IS_LONG_SEQ_SET: bool = false;

/// A single k-mer hit collected by [`RefKmerFilter::get_hits_from_read`],
/// mirroring the FFI-relevant fields of T1K's `struct _hit`
/// (`SeqSet.hpp:66-87`), including `repeats` (`_hit::repeats`,
/// `SeqSet.hpp:72`), populated exactly as stock does (`SeqSet.hpp:
/// 1122,1143,1189,1210`: `repeats = size` of the `seqIndex.Search(...)`
/// result, since `has_hit_in_set`'s `GetHitsFromRead` call always has
/// `barcode == -1` and `puse == NULL`, `SeqSet.hpp:1923`) so
/// [`overlap::get_overlaps_from_hits`]'s `filter == 1` branches (not
/// reachable from `has_hit_in_set`, which always passes `filter == 0`, but
/// ported generally -- see that module's docs) have a correct `repeats` to
/// read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Hit {
    /// The reference sequence index this hit maps to (`_indexInfo::idx`).
    pub idx: u32,
    /// The offset within that reference sequence (`_indexInfo::offset`).
    pub offset: u32,
    /// The offset within the (possibly reverse-complemented) read where this
    /// k-mer window starts (`_hit::readOffset`).
    pub read_offset: i32,
    /// `1` for a forward-strand hit, `-1` for a reverse-complement-strand
    /// hit (`_hit::strand`).
    pub strand: i8,
    /// How many times this hit's k-mer occurs across the whole reference
    /// index (`_hit::repeats`) -- see this struct's doc comment.
    pub repeats: i32,
}

impl Hit {
    /// Converts to [`overlap::OverlapHit`], the standalone hit type
    /// [`overlap::get_overlaps_from_hits`] operates on (see that module's
    /// docs for why it doesn't reuse this crate's `Hit` directly).
    fn to_overlap_hit(self) -> OverlapHit {
        OverlapHit {
            idx: self.idx,
            offset: self.offset,
            read_offset: self.read_offset,
            strand: self.strand,
            repeats: self.repeats,
        }
    }
}

/// Reusable scratch buffers for [`RefKmerFilter`]'s per-read query path
/// (`has_hit_in_set`/`is_good_candidate`), so a caller processing many reads
/// in a loop (a future `extract`/`genotype` command) does not re-allocate a
/// reverse-complement buffer, hit list, and bucket map on every call.
#[derive(Debug, Default, Clone)]
pub struct Scratch {
    /// Reused buffer for the read's reverse complement (mirrors T1K's
    /// caller-owned `rcRead`/`buffer` scratch parameter to `HasHitInSet`/
    /// `GetHitsFromRead`, which is purely an output parameter in the C++ --
    /// see `get_hits_from_read`'s doc comment).
    rc_buf: Vec<u8>,
    /// Reused hit list, cleared and refilled by every `get_hits_from_read` call.
    hits: Vec<Hit>,
    /// Reused `(tag, seqIdx) -> count` bucket map for `has_hit_in_set`'s
    /// add02ca touched-buckets bucket sort (see that function's doc comment).
    buckets: HashMap<(i8, u32), u32>,
}

/// Reference-k-mer read-candidate filter: the pure-Rust port of the
/// `SeqSet` slice `FastqExtractor`/`BamExtractor` use to decide whether a
/// read is a genotyping candidate. See the module docs for exact scope and
/// the deliberate `GetOverlapsFromHits`/`AlignAlgo` scope cut.
pub struct RefKmerFilter {
    /// Forward-code-keyed k-mer position index over every loaded reference
    /// sequence, built via repeated `KmerIndex::build_index_from_read` calls
    /// (one per sequence, `id` = 0-based load order) -- mirrors
    /// `SeqSet::seqIndex` after `InputRefFa` (`SeqSet.hpp:220,872-904`).
    seq_index: KmerIndex,
    /// The k-mer length every sequence was indexed at (`SeqSet::kmerLength`,
    /// `SeqSet.hpp:221`). See the module docs for why this is an explicit
    /// caller-supplied parameter rather than a fixed default.
    kmer_length: usize,
    /// Number of reference sequences loaded (`SeqSet::seqs.size()` after
    /// `InputRefFa`; mirrors `SeqSet::Size()`, `SeqSet.hpp:795-798`).
    seq_count: usize,
    /// 0-based-load-order sequence names, in load order (index == the `idx`
    /// used to key `seq_index`). Not required by the ported slice itself,
    /// but cheap to keep and useful for a future caller that wants to report
    /// which reference sequence a candidate read matched.
    seq_names: Vec<String>,
    /// `SeqSet::hitLenRequired` (`SeqSet.hpp:223`), used by
    /// [`RefKmerFilter::has_hit_in_set`]'s bucket-count gate. Fixed at the
    /// `SeqSet(kmerLength)` constructor default (`SeqSet.hpp:764`); see
    /// [`DEFAULT_HIT_LEN_REQUIRED`]'s doc comment for why this is not
    /// dynamically computed the way `FastqExtractor.cpp`'s `main` does.
    hit_len_required: i32,
    /// `SeqSet::refSeqSimilarity` (`SeqSet.hpp:231`), used by
    /// [`RefKmerFilter::has_hit_in_set`]'s gate 2 mismatch-threshold
    /// computation. Fixed at the `SeqSet(kmerLength)` constructor default;
    /// see [`DEFAULT_REF_SEQ_SIMILARITY`]'s doc comment.
    ref_seq_similarity: f64,
}

impl RefKmerFilter {
    /// Loads a T1K-style reference FASTA (e.g. a `*_seq.fa` produced by
    /// `prepare-ref`/`t1k-build.pl`) and builds a k-mer index over every
    /// sequence, mirroring `SeqSet(kmerLength)` (`SeqSet.hpp:760-772`)
    /// followed by `InputRefFa(filename)` (`SeqSet.hpp:872-904`): each FASTA
    /// record is indexed with `id` = its 0-based position in the file
    /// (matching stock's `id = seqs.size()` at insertion time).
    ///
    /// Sequences must already be uppercase `A`/`C`/`G`/`T`/`N` (T1K assumes
    /// this; see [`crate::kmer::canonical_kmer`]'s doc comment for the same
    /// caveat on the reused `KmerCode`).
    ///
    /// # Errors
    ///
    /// Returns an error if `path` cannot be read, or if the file contains
    /// more sequences than fit in an `i32` (mirrors `KmerIndex`'s `id: i32`
    /// parameter, itself matching T1K's `int id`).
    ///
    /// # Panics
    ///
    /// Panics if `kmer_length == 0` (a zero-length k-mer is meaningless and
    /// would underflow arithmetic throughout this module, matching the
    /// undefined behavior a `kmerLength <= 0` would trigger on the C++
    /// side).
    pub fn from_reference_fasta(path: &Path, kmer_length: usize) -> Result<Self> {
        assert!(kmer_length >= 1, "kmer_length must be >= 1, got {kmer_length}");

        let text = fs::read_to_string(path)
            .with_context(|| format!("reading reference FASTA {}", path.display()))?;

        let mut seq_index = KmerIndex::new();
        let mut seq_names = Vec::new();
        let mut kmer_code = KmerCode::new(kmer_length);

        for record in parse_fasta(&text) {
            let seq_idx = seq_names.len();
            let id = i32::try_from(seq_idx).with_context(|| {
                format!("reference FASTA {} has more than i32::MAX sequences", path.display())
            })?;
            seq_index.build_index_from_read(&mut kmer_code, record.seq.as_bytes(), id, 0);
            seq_names.push(record.id);
        }

        let seq_count = seq_names.len();
        Ok(Self {
            seq_index,
            kmer_length,
            seq_count,
            seq_names,
            hit_len_required: DEFAULT_HIT_LEN_REQUIRED,
            ref_seq_similarity: DEFAULT_REF_SEQ_SIMILARITY,
        })
    }

    /// The k-mer length every reference sequence was indexed at.
    #[must_use]
    pub fn kmer_length(&self) -> usize {
        self.kmer_length
    }

    /// The number of reference sequences loaded (`SeqSet::Size()`).
    #[must_use]
    pub fn seq_count(&self) -> usize {
        self.seq_count
    }

    /// The name of the reference sequence at 0-based load-order index `idx`
    /// (`SeqSet::GetSeqName`, `SeqSet.hpp:810-813`).
    #[must_use]
    pub fn seq_name(&self, idx: usize) -> &str {
        &self.seq_names[idx]
    }

    /// Ergonomic entry point: `!is_low_complexity(read) &&
    /// has_hit_in_set(read)`, mirroring the free function `IsGoodCandidate`
    /// (`FastqExtractor.cpp:113-118`) exactly (including short-circuit
    /// order: `is_low_complexity` is always checked first, so `has_hit_in_set`
    /// is never called for a low-complexity read).
    ///
    /// Allocates a fresh [`Scratch`] internally on every call; for
    /// processing many reads in a loop, prefer
    /// [`RefKmerFilter::is_good_candidate_with_scratch`] with a
    /// caller-owned, reused `Scratch`.
    #[must_use]
    pub fn is_good_candidate(&self, read: &[u8]) -> bool {
        let mut scratch = Scratch::default();
        self.is_good_candidate_with_scratch(read, &mut scratch)
    }

    /// Same as [`RefKmerFilter::is_good_candidate`], but reuses a
    /// caller-owned [`Scratch`] instead of allocating one per call -- the
    /// intended hot-loop entry point for a future `extract`/`genotype`
    /// command processing many reads.
    #[must_use]
    pub fn is_good_candidate_with_scratch(&self, read: &[u8], scratch: &mut Scratch) -> bool {
        !is_low_complexity(read) && self.has_hit_in_set(read, scratch)
    }

    /// Diagnostic-only: evaluates ONLY `HasHitInSet`'s bucket-count gate
    /// (`SeqSet.hpp:1929-1964`, i.e. exactly what Task 3.1's
    /// `has_hit_in_set` computed before this task added gate 2) and returns
    /// whether it alone would accept `read`, WITHOUT running gate 2's
    /// `GetOverlapsFromHits`-based confirmation.
    ///
    /// This exists purely so a differential test can quantify how many reads
    /// gate 2 actually rejects that gate 1 alone would have accepted (i.e.
    /// how many reads exercise the fix this task makes) -- it is not used by
    /// [`RefKmerFilter::has_hit_in_set`]/[`RefKmerFilter::is_good_candidate`]
    /// itself (both call [`RefKmerFilter::bucket_count_gate`] directly, not
    /// this wrapper). See `crates/fg-t1k-sys/tests/diff_refkmerfilter.rs`
    /// for the consumer.
    #[must_use]
    pub fn passes_bucket_count_gate_only(&self, read: &[u8], scratch: &mut Scratch) -> bool {
        self.bucket_count_gate(read, scratch).is_some()
    }

    /// Ported from the bucket-count gate at the top of `SeqSet::HasHitInSet`
    /// (`SeqSet.hpp:1929-1964`): collects hits from both strands, buckets
    /// them by `(strand-tag, seqIdx)`, and finds the single largest bucket.
    /// Returns `None` if that gate rejects the read (mirrors stock's early
    /// `return false`), otherwise `Some((winning_bucket, bucket_size))`.
    fn bucket_count_gate(&self, read: &[u8], scratch: &mut Scratch) -> Option<((i8, u32), u32)> {
        let len = read.len();
        if len < self.kmer_length {
            return None;
        }

        self.get_hits_from_read(read, scratch);
        if scratch.hits.is_empty() {
            return None;
        }

        // Bucket sort: (tag, seqIdx) -> hit count, mirroring stock's
        // `buckets[tag][idx].PushBack(...)` (SeqSet.hpp:1935-1939), except
        // we only need each bucket's SIZE to pick the winner (gate 2 below
        // reconstructs the winning bucket's actual MEMBER hits separately,
        // by filtering `scratch.hits` -- see the comment at that call site),
        // so a count map suffices here.
        scratch.buckets.clear();
        for hit in &scratch.hits {
            let tag = i8::from(hit.strand == 1);
            *scratch.buckets.entry((tag, hit.idx)).or_insert(0) += 1;
        }

        // add02ca touched-buckets selection -- see `select_best_bucket`'s
        // doc comment for why this is byte-identical in outcome to stock's
        // O(seqCnt) full-array scan (SeqSet.hpp:1941-1957).
        let (best_bucket, max) = select_best_bucket(&scratch.buckets)?;

        // SeqSet.hpp:1959: `if (kmerLength * max < hitLenRequired) return false;`
        // Widened to i64 to avoid any theoretical `i32` overflow from the
        // multiplication (stock uses plain `int` arithmetic here; `max` and
        // `kmerLength` are realistically small enough that this never
        // matters in practice, but i64 costs nothing and removes the
        // question).
        let kmer_length_i64 = i64::try_from(self.kmer_length).unwrap_or(i64::MAX);
        if kmer_length_i64 * i64::from(max) < i64::from(self.hit_len_required) {
            return None;
        }

        Some((best_bucket, max))
    }

    /// Ported from `SeqSet::HasHitInSet` in full (`SeqSet.hpp:1915-1990`):
    /// the bucket-count gate ([`RefKmerFilter::bucket_count_gate`]) followed
    /// by the `GetOverlapsFromHits`-based mismatch-threshold confirmation.
    /// See the module docs for the overall two-gate structure and
    /// [`crate::overlap`]'s module docs for why gate 2 needs no
    /// `AlignAlgo`/Smith-Waterman.
    ///
    /// Exposed as `pub(crate)` (rather than fully private) so
    /// `crates/fg-t1k-sys`'s differential test can drive it directly if
    /// useful, matching the established pattern of exposing internals to
    /// same-workspace differential tests without making them part of the
    /// public API surface.
    #[must_use]
    pub(crate) fn has_hit_in_set(&self, read: &[u8], scratch: &mut Scratch) -> bool {
        let len = read.len();
        let Some((best_bucket, _max)) = self.bucket_count_gate(read, scratch) else {
            return false;
        };

        // Gate 2 (SeqSet.hpp:1966-1983): GetOverlapsFromHits on the winning
        // bucket, then accept if any resulting overlap's implied mismatch
        // count is within the similarity-derived threshold.
        //
        // Reconstructing the winning bucket's MEMBER hits: stock builds each
        // `buckets[tag][idx]` by a single forward pass over the full hit
        // list, `PushBack`-ing each hit into its `(tag, idx)` bucket in
        // encounter order (SeqSet.hpp:1935-1939) -- i.e. a stable filter of
        // the full hit list by `(tag, idx)`. `scratch.hits` (filled by
        // `get_hits_from_read` just above) is exactly that full hit list, in
        // the same forward-strand-then-reverse-strand encounter order stock
        // produces it in, so filtering it by `(tag, idx) == best_bucket`
        // here reproduces stock's bucket contents (and their order)
        // exactly, without needing to build all `2 * seqCount` buckets
        // up front.
        let (best_tag, best_idx) = best_bucket;
        let winning_bucket_hits: Vec<OverlapHit> = scratch
            .hits
            .iter()
            .filter(|hit| i8::from(hit.strand == 1) == best_tag && hit.idx == best_idx)
            .map(|hit| hit.to_overlap_hit())
            .collect();

        let overlaps = overlap::get_overlaps_from_hits(
            &winning_bucket_hits,
            self.hit_len_required,
            0,           // filter: HasHitInSet always passes filter=0 (SeqSet.hpp:1968).
            |_idx| true, // is_ref: every sequence loaded via from_reference_fasta is a reference (see module docs).
            DEFAULT_RADIUS,
            IS_LONG_SEQ_SET,
            #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
            {
                self.kmer_length as i32
            },
        );

        // SeqSet.hpp:1973: `int mismatchThreshold = int(len * (1 -
        // refSeqSimilarity)) * kmerLength;` -- truncates `len * (1 -
        // refSeqSimilarity)` to `int` FIRST, THEN multiplies by
        // `kmerLength` (NOT `int(len * (1 - refSeqSimilarity) *
        // kmerLength)`); reproduced in that exact operand order below.
        #[allow(clippy::cast_precision_loss)]
        let len_f64 = len as f64;
        #[allow(clippy::cast_possible_truncation)]
        let truncated = (len_f64 * (1.0 - self.ref_seq_similarity)) as i32;
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let kmer_length_i32 = self.kmer_length as i32;
        let mismatch_threshold = truncated * kmer_length_i32;

        // SeqSet.hpp:1974-1981: `valid = true` if ANY overlap satisfies
        // `len - overlaps[i].matchCnt / 2 <= mismatchThreshold` (integer
        // `matchCnt / 2`, i.e. `matchCnt` counted TWICE per hitLen and
        // halved back here -- see `overlap::Overlap::match_cnt`'s doc
        // comment).
        #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
        let read_len_i32 = len as i32;
        overlaps.iter().any(|o| read_len_i32 - o.match_cnt / 2 <= mismatch_threshold)
    }

    /// Ported from `SeqSet::GetHitsFromRead` (`SeqSet.hpp:1071-1229`),
    /// specialized to exactly the parameters `HasHitInSet` always passes
    /// (`SeqSet.hpp:1923`: `strand=0`, `barcode=-1`, `allowTotalSkip=false`,
    /// `puse=NULL`) -- both strands are searched, no barcode/`puse`
    /// filtering applies, and `allowTotalSkip`'s branch (dead under
    /// `allowTotalSkip=false`) is omitted. Fills `scratch.hits` (cleared
    /// first) and `scratch.rc_buf` (the read's reverse complement, computed
    /// as an output -- mirroring T1K's `rcRead` parameter to
    /// `GetHitsFromRead`/`HasHitInSet`, which callers pass as a pre-allocated
    /// scratch buffer purely to receive the computed reverse complement, not
    /// as an input).
    ///
    /// # `downSample` is always 1 in this port (dead branch omitted)
    ///
    /// Stock's `downSample` (`SeqSet.hpp:1084-1091`) is only ever `> 1` when
    /// `len > 200 && isLongSeqSet`. `isLongSeqSet` defaults to `false`
    /// (`SeqSet.hpp:765`) and has no public setter reachable from the plain
    /// `SeqSet(kmerLength)` construction this port's `from_reference_fasta`
    /// performs -- so `downSample` is unconditionally `1` for every
    /// `RefKmerFilter`, making the downsample-skip branch (`SeqSet.hpp:
    /// 1101-1102,1168-1169`: `if (downSample > 1 && ...) continue;`)
    /// unreachable dead code for this port's scope. It is omitted here
    /// (rather than implemented-but-inert) to keep the loop readable; if a
    /// future caller needs `isLongSeqSet` support, this is the first place
    /// to revisit.
    ///
    /// # `prevKmerCode` is NOT reset between the forward and reverse loops
    ///
    /// Stock only calls `kmerCode.Restart()` before the reverse-strand loop
    /// (`SeqSet.hpp:1160`) -- `prevKmerCode` (the consecutive-duplicate-
    /// window dedup state) carries over from wherever the forward loop left
    /// it. Since the reverse loop's very first search is always forced
    /// anyway (`i == kmerLength - 1` short-circuits the `IsEqual` check),
    /// this stale carry-over is only observable if that very first reverse
    /// window happens to be skipped by the high-frequency-kmer skip-limit
    /// branch below -- in which case `prevKmerCode` stays whatever the
    /// forward loop last left it, exactly matching stock. This port
    /// reproduces that by simply not resetting `prev_kmer_code` between the
    /// two loops (only `kmer_code` gets a fresh `KmerCode::new` +
    /// re-filled window, mirroring `Restart()`).
    fn get_hits_from_read(&self, read: &[u8], scratch: &mut Scratch) {
        scratch.hits.clear();
        let len = read.len();
        let kl = self.kmer_length;

        if len < kl {
            // Defensive guard not present in the C++ (which has no such
            // check inside `GetHitsFromRead` itself -- only `HasHitInSet`'s
            // caller-side `len < kmerLength` check protects it in practice,
            // SeqSet.hpp:1919-1920). Calling this directly with a short read
            // would index `read[i]` for `i >= len` in the fill loop below,
            // which Rust would panic on (a safe failure) where the C++
            // reads past the string's NUL terminator (undefined behavior).
            // `has_hit_in_set` (the only caller in this port) already
            // guards this exact condition before calling, so this is purely
            // a safety net for direct callers of this `pub(crate)` helper
            // (e.g. tests), not an observable behavior change on the tested
            // path.
            return;
        }

        let mut kmer_code = KmerCode::new(kl);
        let mut prev_kmer_code = KmerCode::new(kl);
        let skip_limit = kl / 2;
        let mut skip_cnt = 0usize;

        // Forward strand (SeqSet.hpp:1093-1154).
        let mut i = 0usize;
        while i < kl - 1 {
            kmer_code.append(read[i]);
            i += 1;
        }
        while i < len {
            kmer_code.append(read[i]);
            if i == kl - 1 || !prev_kmer_code.is_equal(&kmer_code) {
                let found = self.seq_index.search(&kmer_code);
                let size = found.len();
                if size >= 100 && i != kl - 1 && i != len - 1 && skip_cnt < skip_limit {
                    skip_cnt += 1;
                    i += 1;
                    continue; // matches C++: also skips the prev_kmer_code update below
                }
                skip_cnt = 0;
                // `i`/`kl` are window positions bounded by `len` (a real
                // read length), so this never approaches i32::MAX; matches
                // the equivalent cast in `kmer_index::build_index_from_read`.
                #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                let read_offset = (i - (kl - 1)) as i32;
                // `repeats = size` (SeqSet.hpp:1122): `HasHitInSet`'s
                // `GetHitsFromRead` call always has `puse == NULL` and
                // `barcode == -1` (SeqSet.hpp:1923), so neither the
                // `puse`-filtered recount (SeqSet.hpp:1123-1132) nor the
                // `barcode != -1` override (SeqSet.hpp:1134-1135) ever
                // applies here -- `repeats` is simply this k-mer's total
                // index hit count.
                #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                let repeats = size as i32;
                for entry in found {
                    scratch.hits.push(Hit {
                        idx: entry.idx,
                        offset: entry.offset,
                        read_offset,
                        strand: 1,
                        repeats,
                    });
                }
            }
            prev_kmer_code = kmer_code.clone();
            i += 1;
        }

        // Reverse complement (SeqSet.hpp:1156, always computed regardless
        // of which strand(s) are searched -- dead work in stock whenever
        // `strand == 1`, but this port only ever needs `strand == 0`, so
        // it's always used here).
        reverse_complement_into(read, &mut scratch.rc_buf);

        // Reverse strand (SeqSet.hpp:1158-1227). `kmer_code.Restart()`
        // equivalent: fresh code/invalid_pos at the same kmer_length.
        // `prev_kmer_code` is deliberately NOT reset -- see this function's
        // doc comment.
        kmer_code = KmerCode::new(kl);
        skip_cnt = 0; // stock explicitly re-zeroes skipCnt here (SeqSet.hpp:1164)
        let mut i = 0usize;
        while i < kl - 1 {
            kmer_code.append(scratch.rc_buf[i]);
            i += 1;
        }
        while i < len {
            kmer_code.append(scratch.rc_buf[i]);
            if i == kl - 1 || !prev_kmer_code.is_equal(&kmer_code) {
                let found = self.seq_index.search(&kmer_code);
                let size = found.len();
                if size >= 100 && i != kl - 1 && i != len - 1 && skip_cnt < skip_limit {
                    skip_cnt += 1;
                    i += 1;
                    continue;
                }
                skip_cnt = 0;
                // `i`/`kl` are window positions bounded by `len` (a real
                // read length), so this never approaches i32::MAX; matches
                // the equivalent cast in `kmer_index::build_index_from_read`.
                #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                let read_offset = (i - (kl - 1)) as i32;
                // `repeats = size` -- see the forward-strand loop's identical
                // comment above (SeqSet.hpp:1189, same reasoning).
                #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
                let repeats = size as i32;
                for entry in found {
                    scratch.hits.push(Hit {
                        idx: entry.idx,
                        offset: entry.offset,
                        read_offset,
                        strand: -1,
                        repeats,
                    });
                }
            }
            prev_kmer_code = kmer_code.clone();
            i += 1;
        }
    }
}

/// Selects the single best `(tag, seqIdx)` bucket -- the one with the most
/// hits -- from a `(tag, seqIdx) -> count` map, replicating stock T1K's
/// `HasHitInSet` bucket-selection loop (`SeqSet.hpp:1941-1957`) bit-for-bit,
/// but WITHOUT stock's `O(seqCnt)` full-array scan. This is the "add02ca"
/// performance fold-in this port requires (named for the maintainer commit
/// that introduced the equivalent optimization): stock allocates
/// `SimpleVector<_hit> buckets[2][seqCnt]` (`SeqSet.hpp:1931-1933`) and scans
/// every one of those `2 * seqCnt` slots regardless of how many of them are
/// actually non-empty; a read only ever touches `hitCnt` distinct buckets
/// (`hitCnt` = the number of k-mer hits collected, independent of `seqCnt`),
/// so tracking only the touched buckets in a `HashMap` and scanning those is
/// `O(hitCnt)` instead of `O(seqCnt)`.
///
/// # Ordering and tie-breaking are what make this byte-identical to stock
///
/// Stock's selection loop is:
/// ```text
/// for k in 0..=1 {           // tag
///     for i in 0..seqCnt {   // seqIdx
///         if size(buckets[k][i]) > 0 && size(buckets[k][i]) > max {
///             maxTag = k; maxSeqIdx = i; max = size(buckets[k][i]);
///         }
///     }
/// }
/// ```
/// i.e. an ascending `(tag, seqIdx)` scan with a **strict `>`** update
/// condition, so the FIRST bucket (in `(tag, seqIdx)` order) to reach the
/// running maximum wins -- later buckets with an EQUAL count never displace
/// it. Every untouched `(tag, seqIdx)` pair has `size == 0` and can
/// therefore never satisfy `size > 0`, so it can never win regardless of how
/// `max` is initialized; restricting the scan to only the touched keys does
/// not change which bucket is selected, AS LONG AS the touched keys are
/// still visited in the SAME ascending `(tag, seqIdx)` order with the SAME
/// strict `>` tie-break. This function does exactly that: it collects the
/// touched keys, sorts them ascending by `(tag, seqIdx)` (Rust's default
/// tuple `Ord` is exactly the lexicographic `(tag, seqIdx)` order stock's
/// nested loop visits), and scans with `count > best_count` (strict).
/// The result -- which specific `(tag, seqIdx)` bucket is "the" winner among
/// ties -- is therefore identical to stock's, not just the winning count.
///
/// (This port's [`RefKmerFilter::has_hit_in_set`] only actually consumes the
/// winning COUNT, not the winning bucket's identity, since the identity is
/// only needed by the excluded `GetOverlapsFromHits` step -- but the
/// identity is computed and returned here anyway, both because it costs
/// nothing extra and because a future `GetOverlapsFromHits` port will need
/// it.)
fn select_best_bucket(buckets: &HashMap<(i8, u32), u32>) -> Option<((i8, u32), u32)> {
    let mut keys: Vec<(i8, u32)> = buckets.keys().copied().collect();
    keys.sort_unstable();

    let mut best: Option<((i8, u32), u32)> = None;
    for key in keys {
        let count = buckets[&key];
        let is_new_max = match best {
            None => true,
            Some((_, best_count)) => count > best_count,
        };
        if is_new_max {
            best = Some((key, count));
        }
    }
    best
}

/// Writes the reverse complement of `seq` into `out` (resizing as needed),
/// mirroring `SeqSet::ReverseComplement` (`SeqSet.hpp:2103-2114`) exactly:
/// each base is complemented (`A<->T`, `C<->G`) except `N`, which maps to
/// `N`. Only `A`/`C`/`G`/`T`/`N` are supported (matching this module's
/// general uppercase-ACGTN-only assumption, documented at
/// [`RefKmerFilter::from_reference_fasta`]); any other byte is out of scope
/// (the C++ original reads `nucToNum[c - 'A']` unconditionally for
/// non-`N` bytes, which is undefined behavior for bytes outside `'A'..='Z'`
/// and produces a garbage-but-defined base for other in-range letters --
/// neither of which this port replicates).
fn reverse_complement_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.resize(seq.len(), 0);
    for (out_i, &c) in seq.iter().rev().enumerate() {
        out[out_i] = complement_base(c);
    }
}

/// Complements a single base, matching `numToNuc[3 - nucToNum[c - 'A']]`
/// (with `N` bypassing the table via `SeqSet::ReverseComplement`'s explicit
/// `if (seq[...] != 'N')` check, `SeqSet.hpp:2108-2111`).
fn complement_base(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        other => panic!(
            "reverse_complement_into: unsupported base {other:?} (only A/C/G/T/N are supported; \
             see this module's doc comment)"
        ),
    }
}

/// Ported verbatim from the free function `IsLowComplexity`
/// (`FastqExtractor.cpp:89-111`, byte-identical to the copy at
/// `BamExtractor.cpp:144-166`): flags a read as low-complexity if any single
/// base makes up at least half the read, if `N`s make up at least a tenth of
/// the read, OR if at least two of the four bases each appear two times or
/// fewer in the whole read.
///
/// Only `A`/`C`/`G`/`T`/`N` are supported (see [`complement_base`]'s doc
/// comment for the same caveat, which applies identically here since both
/// functions mirror the same `nucToNum[c - 'A']` indexing pattern).
#[must_use]
pub fn is_low_complexity(seq: &[u8]) -> bool {
    let mut cnt: [i32; 5] = [0; 5];
    for &c in seq {
        if c == b'N' {
            cnt[4] += 1;
        } else {
            cnt[nuc_index(c)] += 1;
        }
    }

    // Read lengths never approach i32::MAX in practice; matches T1K's own
    // `int i` loop counter in the ported C++ (`FastqExtractor.cpp:92`).
    #[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
    let len = seq.len() as i32;
    if cnt[0] >= len / 2
        || cnt[1] >= len / 2
        || cnt[2] >= len / 2
        || cnt[3] >= len / 2
        || cnt[4] >= len / 10
    {
        return true;
    }

    let low_cnt = cnt[..4].iter().filter(|&&c| c <= 2).count();
    low_cnt >= 2
}

/// Maps `A`/`C`/`G`/`T` to their `nucToNum` index (0-3). Panics on any other
/// byte -- see [`is_low_complexity`]'s doc comment.
fn nuc_index(c: u8) -> usize {
    match c {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        other => panic!(
            "is_low_complexity: unsupported base {other:?} (only A/C/G/T/N are supported; \
             see this module's doc comment)"
        ),
    }
}

/// A single parsed FASTA record: `id` is the header token up to the first
/// whitespace (matching kseq's `name` field, which stock T1K reads via
/// `ReadFiles::Next` -- `ReadFiles.hpp:183`), `seq` is every subsequent line
/// up to (not including) the next `>` header or end of file, concatenated
/// with no separators (matching multi-line FASTA wrapping).
struct FastaRecord {
    id: String,
    seq: String,
}

/// Minimal FASTA parser sufficient for T1K-style reference files (plain
/// text, uppercase `A`/`C`/`G`/`T`/`N`, optionally multi-line-wrapped
/// sequences). Not a general-purpose FASTA/FASTQ reader -- T1K's own
/// `ReadFiles`/`kseq.h` (gzip-transparent, FASTA-and-FASTQ) is considerably
/// more general, but reference FASTAs in this pipeline are always plain-text
/// FASTA, so this suffices without adding a dependency.
fn parse_fasta(text: &str) -> Vec<FastaRecord> {
    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                records.push(FastaRecord { id, seq: std::mem::take(&mut current_seq) });
            }
            current_id = Some(rest.split_whitespace().next().unwrap_or("").to_string());
        } else {
            current_seq.push_str(line.trim_end());
        }
    }
    if let Some(id) = current_id {
        records.push(FastaRecord { id, seq: current_seq });
    }

    records
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_fasta(records: &[(&str, &str)]) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().expect("tempfile");
        for (id, seq) in records {
            writeln!(f, ">{id}").unwrap();
            writeln!(f, "{seq}").unwrap();
        }
        f.flush().unwrap();
        f
    }

    #[test]
    fn parse_fasta_splits_id_at_first_whitespace_and_concatenates_multiline_seq() {
        let text = ">seq1 some comment\nACGT\nACGT\n>seq2\nTTTT\n";
        let records = parse_fasta(text);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].seq, "ACGTACGT");
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].seq, "TTTT");
    }

    #[test]
    fn is_low_complexity_flags_homopolymer() {
        assert!(is_low_complexity(b"AAAAAAAAAAAAAAAAAAAA"));
    }

    #[test]
    fn is_low_complexity_flags_dinucleotide_repeat() {
        // "ATATAT..." has A and T each at ~50% (>= len/2), and C/G at 0
        // (<=2), so two independent triggers both fire on this input.
        assert!(is_low_complexity(b"ATATATATATATATATATAT"));
    }

    #[test]
    fn is_low_complexity_flags_mostly_n() {
        // 20 bases, 10 Ns: cnt[4]=10 >= len/10=2.
        assert!(is_low_complexity(b"NNNNNNNNNNACGTACGTAC"));
    }

    #[test]
    fn is_low_complexity_passes_balanced_sequence() {
        // Roughly balanced composition, no long homopolymer run collapsing
        // into a single dominant base, no more than one base at <=2 count.
        assert!(!is_low_complexity(b"ACGTACGTACGTACGTACGTACGTACGTACGT"));
    }

    #[test]
    fn from_reference_fasta_indexes_every_sequence() {
        let f =
            write_fasta(&[("a", "ACGTACGTACGTACGTACGTACGT"), ("b", "TTTTGGGGCCCCAAAATTTTGGGG")]);
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();
        assert_eq!(filter.seq_count(), 2);
        assert_eq!(filter.kmer_length(), 9);
        assert_eq!(filter.seq_name(0), "a");
        assert_eq!(filter.seq_name(1), "b");
    }

    #[test]
    fn is_good_candidate_true_for_exact_reference_substring() {
        let reference = "ACGTACGTACGTACGTACGTACGTGGATTACAGATTACAGATTACAGATTACAG\
                          CCCTGACGTGTGACGTGTGACGTGTGACGTGTGACGTGT";
        let f = write_fasta(&[("only", reference)]);
        // hitLenRequired defaults to 31; at kmer_length=9, this needs
        // max >= ceil(31/9) = 4 hits in a single bucket, easily satisfied
        // by a long exact substring.
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();

        // A 40bp exact substring of the loaded (and therefore indexed)
        // reference sequence.
        let read = &reference.as_bytes()[10..50];
        assert!(filter.is_good_candidate(read), "exact reference substring should be a candidate");
    }

    #[test]
    fn is_good_candidate_false_for_unrelated_sequence() {
        let reference = "ACGTACGTACGTACGTACGTACGTGGATTACAGATTACAGATTACAGATTACAG\
                          CCCTGACGTGTGACGTGTGACGTGTGACGTGTGACGTGT";
        let f = write_fasta(&[("only", reference)]);
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();

        // 40 bases chosen to avoid containing any 9-mer substring of
        // `reference` (heavy G/C alternation absent from the reference).
        let read = b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";
        assert!(!filter.is_good_candidate(read));
    }

    #[test]
    fn is_good_candidate_false_for_short_read() {
        let reference = "ACGTACGTACGTACGTACGTACGTGGATTACAGATTACA";
        let f = write_fasta(&[("only", reference)]);
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 9).unwrap();
        assert!(!filter.is_good_candidate(b"ACGTACGT")); // len 8 < kmer_length 9
    }

    /// Builds a synthetic reference of 150 distinct 5bp `"AAAAA"` sequences.
    /// Building the k=4 index from a homopolymer sequence exactly `k` bases
    /// long would insert NOTHING (`KmerIndex::build_index_from_read`'s
    /// dedup drops the only window, since it coincidentally matches the
    /// all-zero initial `prevKmerCode` sentinel -- see `kmer_index`'s module
    /// docs); one extra base (5bp, not 4bp) gives a second window that IS
    /// unconditionally inserted (the `i == kl` branch), so each of the 150
    /// sequences contributes exactly one `"AAAA"` index entry, for exactly
    /// 150 entries total -- comfortably over the `size >= 100`
    /// high-frequency threshold this test exercises.
    fn high_frequency_reference() -> tempfile::NamedTempFile {
        let records: Vec<(String, &'static str)> =
            (0..150).map(|i| (format!("seq{i}"), "AAAAA")).collect();
        let record_refs: Vec<(&str, &str)> =
            records.iter().map(|(id, seq)| (id.as_str(), *seq)).collect();
        write_fasta(&record_refs)
    }

    #[test]
    fn high_frequency_kmer_skipped_when_not_at_read_boundary() {
        // Regression test for the SeqSet.hpp:1109/1176 `size >= 100`
        // high-frequency-kmer skip branch (GetHitsFromRead), which is
        // otherwise not naturally exercised by the small differential-test
        // fixture (kir_rna_seq.fa has too few sequences for any k-mer to
        // reach 100 index entries).
        let f = high_frequency_reference();
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 4).unwrap();

        // Sanity: the "AAAA" k-mer really does have >= 100 index entries.
        let mut kc = KmerCode::new(4);
        for &c in b"AAAA" {
            kc.append(c);
        }
        assert!(filter.seq_index.search(&kc).len() >= 100);

        // "AAAA" appears exactly once, at read positions 5-8 (0-based),
        // which is neither the read's first window (i == kmer_length - 1 ==
        // 3) nor its last window (i == len - 1 == 13) -- so it is eligible
        // for the high-frequency skip, and (with skip_cnt starting at 0 <
        // skip_limit = kmer_length/2 = 2) gets skipped. Every other 4-mer in
        // this read is unique (absent from a reference containing only
        // "AAAAA" sequences), so the ONLY possible source of hits is this
        // one "AAAA" window -- if it is correctly skipped, total hits must
        // be exactly 0.
        let read = b"GATTCAAAAGATTC";
        assert_eq!(read.len(), 14);
        let mut scratch = Scratch::default();
        filter.get_hits_from_read(read, &mut scratch);
        assert_eq!(
            scratch.hits.len(),
            0,
            "the sole high-frequency window is mid-read and must be skipped entirely"
        );
    }

    #[test]
    fn high_frequency_kmer_not_skipped_at_first_window() {
        // Same high-frequency reference as above, but this time "AAAA" is
        // the read's very FIRST window (i == kmer_length - 1), which is
        // exempt from the high-frequency skip regardless of `skip_cnt`
        // (SeqSet.hpp:1109: `i != kmerLength - 1 && ...`). Every hit
        // recorded must therefore come from that first window alone (all
        // other windows are unique, absent from the reference).
        let f = high_frequency_reference();
        let filter = RefKmerFilter::from_reference_fasta(f.path(), 4).unwrap();

        let read = b"AAAAGATTC";
        let mut scratch = Scratch::default();
        filter.get_hits_from_read(read, &mut scratch);

        let forward_hits = scratch.hits.iter().filter(|h| h.strand == 1).count();
        assert_eq!(
            forward_hits, 150,
            "the first-window exemption must still record all 150 reference hits"
        );
    }
}
