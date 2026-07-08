//! BAM/CRAM reader ported from the slice of T1K's `Alignments` class
//! (`alignments.hpp`) that `BamExtractor.cpp` uses, over
//! **rust-htslib** (`bam::Reader`) rather than the bundled samtools-0.1.19 /
//! htslib-1.15.1 C library T1K links against directly. Since both are
//! libhts-family readers operating on the same BAM/CRAM record layout, a
//! record parses identically either way; this module's job is to reproduce
//! `Alignments`' *interpretation* of that record (orientation handling,
//! CIGAR-to-segments splitting, the `fragStdev`/`readLen` pre-scan) bit-for-
//! bit, not the low-level parsing itself.
//!
//! # `GetReadSeq`/`GetQual`: reverse-strand records are returned in ORIGINAL
//! sequencing orientation, not reference-forward
//!
//! BAM always stores `SEQ`/`QUAL` in reference-forward orientation
//! (regardless of which strand the read originally sequenced from) --
//! standard SAM spec behavior. T1K's `Alignments::GetReadSeq`
//! (`alignments.hpp:527-563`) and `GetQual` (`alignments.hpp:565-580`) each
//! branch on `IsReverse()`: when the record is NOT reverse-strand, bases/
//! quals are copied out in stored (= original sequencing) order; when it IS
//! reverse-strand, `GetReadSeq` walks the stored SEQ back-to-front AND
//! substitutes each base for its complement (`A<->T`, `C<->G`) -- i.e. a full
//! reverse-complement, undoing BAM's reference-forward storage to recover the
//! original read -- while `GetQual` walks `QUAL` back-to-front with NO value
//! transformation (quality scores have no "complement", only a positional
//! reversal that must track the base reversal). [`Record::read_seq`] and
//! [`Record::qual`] reproduce this exactly. Getting either the RC condition
//! or the base-substitution table wrong silently corrupts ~50% of records
//! (every reverse-strand read) -- this is the single highest-risk trap in
//! this port.
//!
//! # Base decoding: nt16 codes other than A/C/G/T collapse to `'N'`
//!
//! `GetReadSeq`'s inner `switch` on `bam1_seqi(...)` only special-cases the 4
//! canonical nt16 codes (`1`=A, `2`=C, `4`=G, `8`=T per htslib's
//! `seq_nt16_str` table); every other code -- including real IUPAC ambiguity
//! codes (`M`,`R`,`S`,`V`,`W`,`Y`,`H`,`K`,`D`,`B`) that htslib's own decode
//! table (`rust_htslib::bam::record::Seq`) WOULD decode faithfully -- falls
//! to `default: 'N'`. This module therefore does NOT use
//! `rust_htslib::bam::record::Seq::as_bytes()` (which decodes the full IUPAC
//! table), and instead re-implements the 4-case-else-N decode directly from
//! `Seq::encoded_base`.
//!
//! # `GetGeneralInfo`: `fragStdev`/`readLen` pre-scan
//!
//! See [`Alignments::general_info`] for the ported sampling/statistics
//! formula (`alignments.hpp:597-690`).
//!
//! # `GetFieldZ` (`Z:` aux tags) is intentionally NOT ported
//!
//! `BamExtractor.cpp`'s barcode/UMI paths (which are the only callers of
//! `GetFieldZ`) are being dropped in the downstream BamExtractor port (Task
//! 3.3b), so this module does not implement it.

use anyhow::{Context, Result, bail, ensure};
use rust_htslib::bam::{self, Read as _, record::Cigar};
use std::path::{Path, PathBuf};

/// The `@HD` sort/group order of a BAM/CRAM, as it bears on extraction
/// strategy: `Coordinate` needs a 2-pass name-map (seekable file);
/// `QueryName`/`QueryGrouped` keep mates adjacent (1-pass, pipe-able);
/// `Unsorted` (or a missing `@HD SO`/`GO`) guarantees nothing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SortOrder {
    /// `@HD SO:coordinate`.
    Coordinate,
    /// `@HD SO:queryname` (fully name-sorted).
    QueryName,
    /// `@HD GO:query` (grouped by QNAME but not fully sorted).
    QueryGrouped,
    /// `SO:unsorted`/`SO:unknown`, or no `@HD`/`SO`/`GO` at all.
    Unsorted,
}

/// A single reference-coordinate span of an alignment block, mirroring T1K's
/// `_pair64` (`defs.h:23-26`, fields `a`/`b`) as used for `Alignments::segments`.
/// Both bounds are inclusive 0-based reference coordinates (`b == a + len -
/// 1` for a span of `len` reference-consuming bases), matching
/// `alignments.hpp:264-265`'s `segments[segCnt].a = start; segments[segCnt].b
/// = start + len - 1;`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Segment {
    pub a: i64,
    pub b: i64,
}

/// The current record's derived state, recomputed on every
/// [`Alignments::next`] call -- mirrors the public fields T1K's `Alignments`
/// exposes directly (`segments`/`segCnt`) plus the record itself.
struct CurrentRecord {
    record: bam::Record,
    segments: Vec<Segment>,
}

/// Global read/fragment statistics computed by [`Alignments::general_info`],
/// mirroring the public fields `Alignments::fragStdev`/`readLen` (plus
/// `fragLen`/`matePaired`, exposed here too since `BamExtractor.cpp` reads
/// `fragStdev`/`readLen` directly off the struct).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct GeneralInfo {
    /// Representative read length: the MAXIMUM sampled `l_qseq`
    /// (`alignments.hpp:658-659`: `qsort` ascending, then `lens[lensCnt-1]`),
    /// NOT a mean/median.
    pub read_len: i32,
    /// Mean fragment length (`(mateDiff + readLen)` averaged over the
    /// smallest 70% of sampled mate-distance values), or `read_len` itself
    /// when the file is treated as single-end (`alignments.hpp:679-680`).
    pub frag_len: i32,
    /// Fragment-length standard deviation. `0` is T1K's own single-end
    /// sentinel (`alignments.hpp:681`; also forced to `0` at construction,
    /// `alignments.hpp:107`); note the ADDITIONAL `alignments.hpp:685-686`
    /// rule that bumps a genuinely-computed `0` up to `1` for a paired file
    /// (so `frag_stdev == 0` unambiguously means "single-end" to callers,
    /// even in the edge case of a paired file whose sampled stdev rounds
    /// down to exactly zero).
    pub frag_stdev: i32,
    /// Whether the sampled majority of records carry the paired flag
    /// (`alignments.hpp:661`: `hasMateCnt >= totalReadCnt / 2`).
    pub mate_paired: bool,
}

/// A BAM/CRAM reader reproducing the slice of T1K's `Alignments` class that
/// `BamExtractor.cpp` uses (`alignments.hpp`), over
/// **rust-htslib** instead of T1K's bundled htslib/samtools. See the module
/// docs for the two highest-risk semantic traps this port must (and does)
/// reproduce exactly: `GetReadSeq`/`GetQual` orientation, and the
/// `GetGeneralInfo` sampling formula.
pub struct Alignments {
    /// `Some(path)` for a file-backed reader (rewindable via reopen);
    /// `None` for a stdin-backed reader ([`Alignments::from_stdin`]), which
    /// has no path to reopen -- [`Alignments::rewind`] checks this and
    /// returns a clear error rather than attempting anything stdin cannot
    /// support (stdin is a non-seekable, single-consumption stream).
    path: Option<PathBuf>,
    reader: bam::Reader,
    header: bam::HeaderView,
    /// `Some(path)` when an explicit CRAM reference was supplied
    /// ([`Alignments::open_with_reference`] /
    /// [`Alignments::from_stdin_with_reference`]); `None` otherwise. Held so
    /// [`Alignments::rewind`] can re-apply it after reopening the file (a
    /// reopen rebuilds the reader from scratch and would otherwise decode
    /// CRAM against htslib's default reference chain, including the forbidden
    /// EBI network fallback). See [`apply_reference`] for the no-network
    /// coverage preflight this reference goes through on every apply.
    reference: Option<PathBuf>,
    current: Option<CurrentRecord>,
    /// Mirrors `Alignments::totalReadCnt`: counts every PRIMARY record
    /// (`(flag & 0x900) == 0`) seen so far by [`Alignments::next`] since the
    /// last [`Alignments::rewind`]/construction (`alignments.hpp:184-207`).
    /// [`Alignments::general_info`] also increments this field via its own
    /// scan (`alignments.hpp:652`), matching stock's shared counter.
    total_read_cnt: u64,
}

impl Alignments {
    /// Opens `path` for reading, mirroring `Alignments::Open(char*)`
    /// (`alignments.hpp:124-128`) -- rust-htslib's `bam::Reader::from_path`
    /// auto-detects BAM/SAM/CRAM from the file's magic bytes, matching
    /// htslib's own `sam_open` dispatch that the C++ side relies on.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or its header cannot be
    /// parsed.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        Self::open_with_reference(path, None::<&Path>)
    }

    /// Opens `path` like [`Alignments::open`], but additionally binds an
    /// explicit reference FASTA for decoding CRAM (via
    /// `rust_htslib::bam::Reader::set_reference`). BAM/SAM inputs ignore the
    /// reference (they carry their own sequence), but supplying it is
    /// harmless.
    ///
    /// When `reference` is `Some`, decoding goes through a HARD no-network
    /// preflight ([`apply_reference`]): the reference's `.fai` sibling must
    /// exist and must cover every `@SQ` contig in the header, or opening
    /// fails outright -- htslib's default CRAM reference chain (which can fall
    /// back to an `M5`/`REF_PATH` lookup against the EBI network endpoint) is
    /// never consulted. The reference is remembered so [`Alignments::rewind`]
    /// can re-apply it after reopening the file.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened, its header cannot be
    /// parsed, or (when `reference` is `Some`) the no-network coverage
    /// preflight in [`apply_reference`] fails.
    pub fn open_with_reference<P: AsRef<Path>, Q: AsRef<Path>>(
        path: P,
        reference: Option<Q>,
    ) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let reference = reference.map(|r| r.as_ref().to_path_buf());
        let mut reader = bam::Reader::from_path(&path)
            .with_context(|| format!("opening alignment file {}", path.display()))?;
        let header = bam::HeaderView::from_header(&bam::Header::from_template(reader.header()));
        apply_reference(&mut reader, &header, reference.as_deref())?;
        Ok(Self { path: Some(path), reader, header, reference, current: None, total_read_cnt: 0 })
    }

    /// Opens standard input as a BAM/CRAM stream, for a pipe-fed source
    /// (`-` on the CLI). Uses `rust_htslib::bam::Reader::from_stdin()`
    /// directly -- **not** [`Alignments::open`]`("-")`, which would call
    /// rust-htslib's `bam::Reader::from_path` and therefore
    /// `path_as_bytes(path, must_exist=true)`; `Path::new("-").exists()` is
    /// `false` (there is no file literally named `-`), so that call fails
    /// with a `FileNotFound` error instead of reading the pipe.
    ///
    /// The returned reader is **not rewindable**: [`Alignments::rewind`]
    /// returns an error rather than silently reopening a real file, since
    /// stdin is a non-seekable, single-consumption stream -- there is
    /// nothing to reopen. This also means [`Alignments::general_info`] is
    /// effectively unusable here in practice: every caller of it in this
    /// crate immediately calls `rewind()` afterward (mirroring
    /// `BamExtractor.cpp:573-574`'s `GetGeneralInfo(true); Rewind();`
    /// pattern), and that follow-up `rewind()` call fails on a stdin-backed
    /// reader -- so a stdin source must instead go through a one-pass entry
    /// point that derives its setup statistics from a bounded, non-rewound
    /// head buffer (e.g.
    /// [`crate::bam_extract::extract_from_bam_no_alignment_grouped`]) rather
    /// than `general_info`.
    ///
    /// # Errors
    ///
    /// Returns an error if stdin cannot be opened as a BAM/CRAM stream or its
    /// header cannot be parsed.
    pub fn from_stdin() -> Result<Self> {
        Self::from_stdin_with_reference(None::<&Path>)
    }

    /// Opens standard input like [`Alignments::from_stdin`], but additionally
    /// binds an explicit reference FASTA for decoding a CRAM stream, going
    /// through the same HARD no-network preflight as
    /// [`Alignments::open_with_reference`] (see [`apply_reference`]). As with
    /// [`Alignments::from_stdin`], the returned reader is **not rewindable**:
    /// [`Alignments::rewind`] returns an error rather than reopening.
    ///
    /// # Errors
    ///
    /// Returns an error if stdin cannot be opened as a BAM/CRAM stream, its
    /// header cannot be parsed, or (when `reference` is `Some`) the no-network
    /// coverage preflight in [`apply_reference`] fails.
    pub fn from_stdin_with_reference<Q: AsRef<Path>>(reference: Option<Q>) -> Result<Self> {
        let reference = reference.map(|r| r.as_ref().to_path_buf());
        let mut reader =
            bam::Reader::from_stdin().context("opening alignment stream from stdin")?;
        let header = bam::HeaderView::from_header(&bam::Header::from_template(reader.header()));
        apply_reference(&mut reader, &header, reference.as_deref())?;
        Ok(Self { path: None, reader, header, reference, current: None, total_read_cnt: 0 })
    }

    /// Returns the BAM/CRAM's declared sort/group order from its `@HD` line.
    /// `SO:coordinate`/`SO:queryname` take precedence; else `GO:query` →
    /// [`SortOrder::QueryGrouped`]; else [`SortOrder::Unsorted`]. Never does
    /// I/O — reads the already-parsed header text.
    #[must_use]
    pub fn sort_order(&self) -> SortOrder {
        let text = self.header.as_bytes();
        for line in text.split(|&b| b == b'\n') {
            if !line.starts_with(b"@HD") {
                continue;
            }
            let mut so: Option<&[u8]> = None;
            let mut go: Option<&[u8]> = None;
            for field in line.split(|&b| b == b'\t') {
                if let Some(v) = field.strip_prefix(b"SO:") {
                    so = Some(v);
                } else if let Some(v) = field.strip_prefix(b"GO:") {
                    go = Some(v);
                }
            }
            return match so {
                Some(b"coordinate") => SortOrder::Coordinate,
                Some(b"queryname") => SortOrder::QueryName,
                _ => match go {
                    Some(b"query") => SortOrder::QueryGrouped,
                    _ => SortOrder::Unsorted,
                },
            };
        }
        SortOrder::Unsorted
    }

    /// Mirrors `Alignments::IsOpened`. Always `true` for a live [`Alignments`]
    /// (construction fails outright on open error, unlike the C++ side which
    /// separates `Open`'s `exit(1)`-on-failure from a boolean flag) -- kept
    /// for API symmetry with callers ported from `BamExtractor.cpp` that
    /// check it.
    #[must_use]
    pub fn is_opened(&self) -> bool {
        true
    }

    /// Mirrors `Alignments::Rewind` (`alignments.hpp:130-137`): closes and
    /// reopens the underlying file, resetting iteration to the beginning.
    /// rust-htslib readers are not seekable-to-start without an index, so
    /// reopen-from-path is the direct equivalent of T1K's own
    /// close-then-reopen rewind (T1K's `Rewind` is likewise a full
    /// close+reopen, not a `gzrewind`-style in-place seek, despite the name).
    /// Also resets `total_read_cnt` to `0`, matching `Next`'s own
    /// `atBegin`-gated reset (`alignments.hpp:184-185`) -- the first `next()`
    /// call after a fresh open/rewind always starts counting from zero.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be reopened, or if this reader
    /// was constructed via [`Alignments::from_stdin`] (stdin is not
    /// seekable/reopenable -- see that method's doc comment).
    pub fn rewind(&mut self) -> Result<()> {
        let Some(path) = self.path.clone() else {
            bail!(
                "cannot rewind: this Alignments reader was opened from stdin (Alignments::from_stdin), \
                 which is not seekable/reopenable; use a one-pass entry point that never rewinds for \
                 a stdin source"
            );
        };
        let mut reader = bam::Reader::from_path(&path)
            .with_context(|| format!("rewinding alignment file {}", path.display()))?;
        let header = bam::HeaderView::from_header(&bam::Header::from_template(reader.header()));
        // Reopening builds a fresh reader with no reference; re-apply the
        // explicit CRAM reference (if any) so pass 2 of a coordinate-CRAM
        // 2-pass extract decodes THROUGH the supplied reference, never
        // htslib's default (network-capable) reference chain.
        apply_reference(&mut reader, &header, self.reference.as_deref())?;
        self.header = header;
        self.reader = reader;
        self.current = None;
        // C++ `Rewind` does NOT reset `totalReadCnt`; instead the first `Next`
        // after rewind resets it via an `atBegin` gate. We zero it here, which
        // is behavior-equivalent for the ported consumer: `BamExtractor` never
        // reads the counter between `Rewind` and the first `Next`, and after
        // that first `Next` both sides have counted exactly one record.
        self.total_read_cnt = 0;
        Ok(())
    }

    /// Mirrors `Alignments::Next` (`alignments.hpp:178-314`): advances to the
    /// next record in file order (ALL records, primary and
    /// secondary/supplementary alike -- the `(b->core.flag & 0x900) == 0`
    /// check only gates the `totalReadCnt` counter, never filters which
    /// records `Next` yields; see the commented-out `if (b->core.flag &
    /// 0xC) continue;` immediately below it in the vendored source, dead
    /// code that confirms no read is ever skipped here), recomputing
    /// `segments`/`segCnt` from its CIGAR. Returns `false` at EOF (mirrors
    /// `Next` returning `0`).
    ///
    /// # Errors
    ///
    /// Returns an error if the underlying record read fails (a genuine parse
    /// error, distinct from a clean EOF).
    // Named `next` (not e.g. `advance`) deliberately, mirroring T1K's own
    // `Alignments::Next()` name exactly for port-fidelity/searchability; its
    // `Result<bool>` signature does not match `Iterator::next`'s
    // `Option<Self::Item>`, so there is no real risk of confusing the two.
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Result<bool> {
        let mut record = bam::Record::new();
        match self.reader.read(&mut record) {
            None => {
                self.current = None;
                return Ok(false);
            }
            Some(Err(e)) => return Err(e).context("reading BAM/CRAM record"),
            Some(Ok(())) => {}
        }

        // alignments.hpp:206-207: only PRIMARY records advance totalReadCnt.
        if (record.flags() & 0x900) == 0 {
            self.total_read_cnt += 1;
        }

        let segments = cigar_to_segments(record.pos(), &record.cigar());
        self.current = Some(CurrentRecord { record, segments });
        Ok(true)
    }

    /// Panics with a clear message if [`Alignments::next`] has not yet been
    /// called (or the reader is at EOF) -- every accessor below requires a
    /// current record, matching T1K's own precondition (the C++ side simply
    /// dereferences `b`, which would be a null-pointer read before the first
    /// `Next()`).
    fn current(&self) -> &CurrentRecord {
        self.current.as_ref().expect("Alignments accessor called with no current record (call next() first, and check its return value)")
    }

    /// Mirrors `Alignments::segCnt` / `Alignments::segments` -- the
    /// reference-coordinate spans of the current record's alignment blocks,
    /// split from its CIGAR (see [`cigar_to_segments`] for the exact
    /// per-op rules). Empty for an unmapped record with no CIGAR ops.
    #[must_use]
    pub fn segments(&self) -> &[Segment] {
        &self.current().segments
    }

    /// Mirrors `Alignments::GetChromId` (`alignments.hpp:317-320`): the
    /// current record's reference ID (`-1` if unmapped/no reference).
    #[must_use]
    pub fn chrom_id(&self) -> i32 {
        self.current().record.tid()
    }

    /// Mirrors `Alignments::GetChromName` (`alignments.hpp:322-325`): the
    /// target name for reference id `tid`.
    ///
    /// # Panics
    ///
    /// Panics if `tid` is out of range for the header's target list (matches
    /// the C++ side's unchecked `bHdr->target_name[tid]` array access, which
    /// is undefined behavior on an out-of-range `tid` -- this port fails loud
    /// instead of silently reading garbage).
    #[must_use]
    pub fn chrom_name(&self, tid: i32) -> String {
        let tid = u32::try_from(tid).unwrap_or_else(|_| panic!("chrom_name: negative tid {tid}"));
        let name = self.header.tid2name(tid);
        String::from_utf8_lossy(name).into_owned()
    }

    /// Mirrors `Alignments::GetChromIdFromName` (`alignments.hpp:327-348`):
    /// looks up a reference id by name, with T1K's specific 3-way fallback
    /// chain (exact match; then, if `s` is at least 4 characters, the
    /// substring starting at its 4th character, i.e. stripping a `"chr"`-like
    /// 3-character prefix REGARDLESS of what those 3 characters actually are;
    /// then `"chr" + s`).
    ///
    /// # Errors
    ///
    /// Returns an error if none of the three lookups match any header target
    /// (mirrors the C++ side's `fprintf`+`exit(1)` on an unknown name, but as
    /// a recoverable `Result` rather than a hard process exit).
    pub fn chrom_id_from_name(&self, s: &str) -> Result<i32> {
        if let Some(tid) = self.header.tid(s.as_bytes()) {
            return Ok(i32::try_from(tid).unwrap_or(i32::MAX));
        }
        if s.len() >= 4 {
            let stripped = &s[3..];
            if let Some(tid) = self.header.tid(stripped.as_bytes()) {
                return Ok(i32::try_from(tid).unwrap_or(i32::MAX));
            }
        }
        let prefixed = format!("chr{s}");
        if let Some(tid) = self.header.tid(prefixed.as_bytes()) {
            return Ok(i32::try_from(tid).unwrap_or(i32::MAX));
        }
        bail!("Unknown genome name: {s}");
    }

    /// Mirrors `Alignments::GetChromLength` (`alignments.hpp:350-353`).
    ///
    /// # Panics
    ///
    /// Panics if `tid` is out of range (mirrors the C++ side's unchecked
    /// array access, same rationale as [`Alignments::chrom_name`]).
    #[must_use]
    pub fn chrom_length(&self, tid: i32) -> u64 {
        let tid = u32::try_from(tid).unwrap_or_else(|_| panic!("chrom_length: negative tid {tid}"));
        self.header
            .target_len(tid)
            .unwrap_or_else(|| panic!("chrom_length: tid {tid} out of range"))
    }

    /// Mirrors `Alignments::GetChromCount` (`alignments.hpp:355-358`).
    #[must_use]
    pub fn chrom_count(&self) -> i32 {
        i32::try_from(self.header.target_count()).unwrap_or(i32::MAX)
    }

    /// Mirrors `Alignments::IsFirstMate` (`alignments.hpp:405-410`): the
    /// `0x40` (READ1) flag bit.
    #[must_use]
    pub fn is_first_mate(&self) -> bool {
        self.current().record.is_first_in_template()
    }

    /// Mirrors `Alignments::IsReverse` (`alignments.hpp:412-417`): the `0x10`
    /// flag bit.
    #[must_use]
    pub fn is_reverse(&self) -> bool {
        self.current().record.is_reverse()
    }

    /// Mirrors `Alignments::IsMateReverse` (`alignments.hpp:419-424`): the
    /// `0x20` flag bit.
    #[must_use]
    pub fn is_mate_reverse(&self) -> bool {
        self.current().record.is_mate_reverse()
    }

    /// Mirrors `Alignments::IsTemplateAligned` (`alignments.hpp:426-432`)
    /// EXACTLY: `false` if `(flag & 0xd) == 0xd` (paired + unmapped +
    /// mate-unmapped all set) OR `(flag & 0x5) == 0x4` (unpaired AND
    /// unmapped) OR `tid < 0`; `true` otherwise. Note this is subtly
    /// different from "both mates aligned" -- e.g. a PAIRED record that is
    /// itself unmapped but whose MATE is mapped (`flag & 0xd == 0x9`, not
    /// `0xd`) is still considered template-aligned by this predicate.
    #[must_use]
    pub fn is_template_aligned(&self) -> bool {
        let flag = self.current().record.flags();
        let tid = self.current().record.tid();
        if (flag & 0xd) == 0xd || (flag & 0x5) == 0x4 || tid < 0 {
            return false;
        }
        true
    }

    /// Mirrors `Alignments::IsAligned` (`alignments.hpp:434-439`): `false` if
    /// the unmapped flag (`0x4`) is set OR `tid < 0`; `true` otherwise.
    #[must_use]
    pub fn is_aligned(&self) -> bool {
        let flag = self.current().record.flags();
        let tid = self.current().record.tid();
        !(flag & 0x4 != 0 || tid < 0)
    }

    /// Mirrors `Alignments::IsPrimary` (`alignments.hpp:458-464`): `true`
    /// unless either the secondary (`0x100`) or supplementary (`0x800`) flag
    /// bit is set.
    #[must_use]
    pub fn is_primary(&self) -> bool {
        (self.current().record.flags() & 0x900) == 0
    }

    /// The `0x1` (paired) FLAG bit of the current record: whether it is part
    /// of a template with an expected mate. Not a direct port of a named
    /// T1K `Alignments` accessor (T1K reads this bit inline off the raw
    /// `bam1_t` where needed); exposed here for
    /// [`crate::bam_extract::extract_from_bam_no_alignment_grouped`], which
    /// -- unable to call [`Alignments::general_info`] on a non-rewindable
    /// (possibly stdin) reader -- derives `single_end` and its per-group
    /// pairing rule directly from this bit, sampled over a bounded head,
    /// using the same `has_mate_cnt >= total / 2` majority-vote rule
    /// [`Alignments::general_info`] uses internally.
    #[must_use]
    pub fn is_paired(&self) -> bool {
        (self.current().record.flags() & 0x1) != 0
    }

    /// Mirrors `Alignments::GetReadId` (`alignments.hpp:441-444`): the QNAME.
    #[must_use]
    pub fn read_id(&self) -> String {
        String::from_utf8_lossy(self.current().record.qname()).into_owned()
    }

    /// Mirrors `Alignments::GetReadSeq` (`alignments.hpp:527-563`) EXACTLY,
    /// including the reverse-strand reverse-complement (see module docs) and
    /// the "any non-ACGT nt16 code decodes to `'N'`" rule (module docs).
    #[must_use]
    pub fn read_seq(&self) -> Vec<u8> {
        let record = &self.current().record;
        let seq = record.seq();
        let len = seq.len();
        if record.is_reverse() {
            // alignments.hpp:548: `for (i=0, j=l_qseq-1; j>=0; ++i,--j)`
            (0..len).map(|i| decode_reverse_complement(seq.encoded_base(len - 1 - i))).collect()
        } else {
            (0..len).map(|i| decode_forward(seq.encoded_base(i))).collect()
        }
    }

    /// Mirrors `Alignments::GetQual` (`alignments.hpp:565-580`) EXACTLY:
    /// phred-encoded (raw, 0-based) quality bytes converted to ASCII
    /// (phred+33), reversed (but NOT complemented -- qualities have no base
    /// identity) on a reverse-strand record, matching `GetReadSeq`'s base
    /// reversal position-for-position.
    #[must_use]
    pub fn qual(&self) -> Vec<u8> {
        let record = &self.current().record;
        let raw = record.qual();
        let len = raw.len();
        if record.is_reverse() {
            (0..len).map(|i| raw[len - 1 - i].wrapping_add(33)).collect()
        } else {
            raw.iter().map(|q| q.wrapping_add(33)).collect()
        }
    }

    /// Mirrors `Alignments::GetGeneralInfo` (`alignments.hpp:597-690`)
    /// EXACTLY: samples records (via repeated internal `next()`-equivalent
    /// reads, consuming the reader's current position -- callers must
    /// `rewind()` afterward if they need to re-read from the start, same as
    /// `BamExtractor.cpp:573-574`'s `GetGeneralInfo(true); Rewind();`
    /// pattern) to compute `read_len`/`frag_len`/`frag_stdev`/`mate_paired`.
    ///
    /// # Sampling loop
    ///
    /// Reads records ONE AT A TIME (not via [`Alignments::next`] -- this
    /// method has its own internal read loop, mirroring the C++ side's
    /// separate `while` loop over `sam_read1`, `alignments.hpp:609-655`),
    /// skipping non-primary records entirely (they are read but neither
    /// counted nor sampled -- `alignments.hpp:629-630`: `if ((flag & 0x900)
    /// == 0) break;` inside the inner `while`, i.e. the outer loop body only
    /// runs for primary records). For each PRIMARY record:
    /// - Always samples `l_qseq` into the length list, up to `sampleMax =
    ///   100_000` entries (`alignments.hpp:635-639`).
    /// - Samples a mate-distance value ONLY when: same-chromosome pair
    ///   (`tid == mtid`), this record's `pos < mate's pos` (avoids
    ///   double-counting each pair), AND the two mates are on OPPOSITE
    ///   strands (`IsReverse() != IsMateReverse()`, `alignments.hpp:641-647`
    ///   -- an explicit anti-chimeric guard: same-strand "pairs" are
    ///   excluded from the distance sample even though they still count
    ///   toward `hasMateCnt` below).
    /// - Counts toward `hasMateCnt` whenever the paired flag (`0x1`) is set,
    ///   REGARDLESS of the same-strand/mate-distance-sampling guard above
    ///   (`alignments.hpp:649-650`) -- this is what later decides
    ///   single-vs-paired, independent of whether any distance sample was
    ///   actually collected.
    /// - Increments `total_read_cnt` every primary record
    ///   (`alignments.hpp:652`, the SAME counter [`Alignments::next`]
    ///   maintains).
    ///
    /// If `stop_early` is `true`, the loop additionally breaks once
    /// `total_read_cnt >= sampleMax` (`alignments.hpp:653-654`) -- but ONLY
    /// after that record has already been fully sampled/counted, so this
    /// never truncates mid-record; it just caps how many records are visited
    /// after the 100,000th sample slot fills.
    ///
    /// # Statistics
    ///
    /// - `read_len = max(sampled l_qseq values)` (`alignments.hpp:658-659`:
    ///   ascending `qsort` then take the last element -- the maximum, NOT a
    ///   median/mean despite superficially resembling a median computation).
    /// - `mate_paired = hasMateCnt >= total_read_cnt / 2` (INTEGER division,
    ///   `alignments.hpp:661`).
    /// - If `mate_paired`: sort the sampled mate-distance values ascending,
    ///   average `(distance + read_len)` over the smallest `70%` count
    ///   `k = ceil(mateDiffCnt * 0.7)`. C's `for (i=0; i < mateDiffCnt*0.7; ++i)`
    ///   is a loop-bound comparison of an `int` counter against a `double`
    ///   bound, so it counts every `i` strictly below the bound -- i.e. `k` is
    ///   the bound itself when it is an exact integer (multiples of 10), else
    ///   `floor + 1`. This is NOT a truncating cast (a `floor` would be off by
    ///   one for every non-multiple-of-10). Then `frag_len` = integer division
    ///   of the sum by `k`, and `frag_stdev = floor(sqrt(sumsq/k - fragLen^2))`
    ///   (`alignments.hpp:669-676`, all integer arithmetic until the final `sqrt`).
    /// - Else (not `mate_paired`): `frag_len = read_len`, `frag_stdev = 0`
    ///   (`alignments.hpp:679-682`).
    /// - Special case: if `mate_paired` AND the computed `frag_stdev` is
    ///   exactly `0`, it is bumped to `1` (`alignments.hpp:685-686`) -- so a
    ///   `frag_stdev == 0` result unambiguously means "single-end" to
    ///   callers (`BamExtractor.cpp` branches on exactly this: `if
    ///   (alignments.fragStdev == 0)`), never "paired but zero-variance".
    ///
    /// If the sampled-length list ends up empty (an empty file), this
    /// mirrors the C++ side's own undefined behavior at `lens[lensCnt - 1]`
    /// on `lensCnt == 0` by returning a zeroed [`GeneralInfo`] instead
    /// (a defined, if divergent-from-UB, choice for a degenerate input T1K
    /// itself does not handle safely).
    ///
    /// # Errors
    ///
    /// Returns an error if the underlying record read fails (a genuine parse
    /// error).
    pub fn general_info(&mut self, stop_early: bool) -> Result<GeneralInfo> {
        const SAMPLE_MAX: usize = 100_000;
        let mut lens: Vec<i32> = Vec::new();
        let mut mate_diff: Vec<i64> = Vec::new();
        let mut has_mate_cnt: u64 = 0;

        loop {
            let mut record = bam::Record::new();
            match self.reader.read(&mut record) {
                None => break,
                Some(Err(e)) => return Err(e).context("reading BAM/CRAM record in general_info"),
                Some(Ok(())) => {}
            }
            // alignments.hpp:629-630: only primary records reach the sampling
            // logic below; non-primary records are read but otherwise
            // ignored (not even counted).
            if (record.flags() & 0x900) != 0 {
                continue;
            }

            if lens.len() < SAMPLE_MAX {
                // alignments.hpp:637: `lens[lensCnt] = b->core.l_qseq;` --
                // the raw stored sequence length (`record.seq_len()` mirrors
                // `l_qseq` directly).
                lens.push(i32::try_from(record.seq_len()).unwrap_or(i32::MAX));
            }

            if mate_diff.len() < SAMPLE_MAX
                && record.tid() == record.mtid()
                && record.pos() < record.mpos()
                && record.is_reverse() != record.is_mate_reverse()
            {
                mate_diff.push(record.mpos() - record.pos());
            }

            if (record.flags() & 0x1) != 0 {
                has_mate_cnt += 1;
            }

            self.total_read_cnt += 1;
            if stop_early && self.total_read_cnt >= SAMPLE_MAX as u64 {
                break;
            }
        }

        if lens.is_empty() {
            return Ok(GeneralInfo::default());
        }

        lens.sort_unstable();
        let read_len = lens[lens.len() - 1];

        let mate_paired = has_mate_cnt >= self.total_read_cnt / 2;
        let (frag_len, mut frag_stdev) = if mate_paired {
            mate_diff.sort_unstable();
            // alignments.hpp:669: `for (i = 0; i < mateDiffCnt * 0.7; ++i) ; k = i;`
            // This is a loop-bound comparison, NOT a truncating cast: it counts
            // every integer `i` strictly less than the f64 bound
            // `mateDiffCnt * 0.7`, so the resulting `k` is `ceil(mateDiffCnt*0.7)`
            // -- equal to the bound when it lands exactly on an integer (e.g.
            // multiples of 10), else `floor + 1`. Using `as usize` (floor) here
            // would be off by one for every `mateDiffCnt` that is not a multiple
            // of 10, corrupting both `frag_len` and `frag_stdev`. Rust computes
            // the same `len as f64 * 0.7` product the C++ loop compares against,
            // so `ceil` reproduces the loop's iteration count bit-for-bit.
            #[allow(
                clippy::cast_precision_loss,
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss
            )]
            let k = (mate_diff.len() as f64 * 0.7).ceil() as usize;
            let k = k.min(mate_diff.len());
            if k == 0 {
                // C++ divides by `k` unconditionally here; k==0 would be a
                // division-by-zero/UB on the C++ side for a pathological
                // input with hasMateCnt>=total/2 but zero valid distance
                // samples. Match T1K's actual observable behavior for a
                // realistic mixed input by falling back to the single-end
                // formula rather than panicking/dividing by zero.
                (read_len, 0)
            } else {
                // Accumulate in i128, not i64: `v` is a mate distance plus the
                // read length, and `sumsq` sums `v * v` over up to `SAMPLE_MAX`
                // terms. On a whole-genome BAM with large mate distances that
                // sum can exceed `i64::MAX` -- a debug-build overflow panic (and
                // a silent wrap in release). i128 keeps the accumulation
                // overflow-free for any realistic input. (The vendored C++
                // squares `(mateDiff[i] + readLen)` in 32-bit `int` before
                // widening to `long long`, so it already diverges from this port
                // once a per-term value exceeds ~46340; the fixtures exercise
                // only short fragments, where every accumulator width agrees.)
                let mut sum: i128 = 0;
                let mut sumsq: i128 = 0;
                for &d in &mate_diff[..k] {
                    let v = i128::from(d) + i128::from(read_len);
                    sum += v;
                    sumsq += v * v;
                }
                let k_i128 = i128::try_from(k).unwrap_or(i128::MAX);
                let frag_len_i128 = sum / k_i128;
                #[allow(clippy::cast_precision_loss)]
                let variance = (sumsq / k_i128 - frag_len_i128 * frag_len_i128) as f64;
                #[allow(clippy::cast_possible_truncation)]
                let stdev = variance.sqrt() as i32;
                let frag_len = i32::try_from(frag_len_i128).unwrap_or(i32::MAX);
                (frag_len, stdev)
            }
        } else {
            (read_len, 0)
        };

        // alignments.hpp:685-686.
        if mate_paired && frag_stdev == 0 {
            frag_stdev = 1;
        }

        Ok(GeneralInfo { read_len, frag_len, frag_stdev, mate_paired })
    }

    /// Mirrors `Alignments::totalReadCnt` as observed by callers reading the
    /// field directly (`BamExtractor.cpp` does not do so today, but this is
    /// exposed for completeness/testability).
    #[must_use]
    pub fn total_read_cnt(&self) -> u64 {
        self.total_read_cnt
    }

    /// Samples read-1 sequence lengths from this BAM/CRAM's PRIMARY records,
    /// for [`crate::bam_extract::compute_hit_len_required_no_alignment`] to
    /// reproduce the FASTQ path's `hitLenRequired` formula
    /// (`crate::extract`'s `sample_head`) over BAM input instead of a FASTQ
    /// reader.
    ///
    /// Iterates records via [`Alignments::next`], skipping non-primary
    /// records (`is_primary()`) WITHOUT counting them against `limit` --
    /// mirroring [`Alignments::general_info`]'s own primary-only gate, where
    /// secondary/supplementary records are read but neither counted nor
    /// sampled. Stops once `limit` PRIMARY records have been seen (or at
    /// EOF); a limit of [`crate::extract::HIT_LEN_SAMPLE_SIZE`] mirrors the
    /// FASTQ path's own "first 1000" sample size.
    ///
    /// For each primary record, `read_seq().len()` is summed into the
    /// returned total ONLY for a "read-1" record: the first mate
    /// (`is_first_mate()`) of a paired-flag (`0x1`) template, or ANY record
    /// of an unpaired (paired-flag unset) template -- generalizing the FASTQ
    /// path's "read-1 lengths only" rule on a per-record basis, since a BAM
    /// (unlike a FASTQ pair of files) can in principle mix paired and
    /// unpaired records. Returns `(sum_of_read1_seq_lengths, count)`, the
    /// exact two values [`crate::bam_extract::compute_hit_len_required_no_alignment`]
    /// takes as `sampled_read1_len_sum`/`sampled_count`.
    ///
    /// Does NOT rewind: the caller must call [`Alignments::rewind`]
    /// afterward if it needs to re-read from the start, matching how
    /// [`Alignments::general_info`] and `sample_head` both leave rewinding to
    /// their callers.
    ///
    /// # Errors
    ///
    /// Returns an error if the underlying record read fails (a genuine parse
    /// error, distinct from a clean EOF).
    pub fn sample_read1_len_sum(&mut self, limit: usize) -> Result<(i64, usize)> {
        let mut sum: i64 = 0;
        let mut count: usize = 0;
        let mut primary_seen: usize = 0;

        while primary_seen < limit && self.next()? {
            if !self.is_primary() {
                continue;
            }
            primary_seen += 1;

            let is_read1 = !self.current().record.is_paired() || self.is_first_mate();
            if is_read1 {
                // `seq_len()` returns `l_qseq` directly (O(1), no allocation);
                // `read_seq()` would decode every base into a `Vec<u8>` just to
                // take its length. `general_info` samples lengths the same way.
                sum += i64::try_from(self.current().record.seq_len()).unwrap_or(i64::MAX);
                count += 1;
            }
        }

        Ok((sum, count))
    }
}

/// Parses the contig NAME column (first tab-separated field per line) from
/// the `.fai` index sibling of `ref_path` (i.e. `<ref_path>.fai`), returning
/// the set of contig names it declares.
///
/// # Errors
///
/// Returns an error if `<ref_path>.fai` does not exist (naming the missing
/// index and suggesting `samtools faidx`) or cannot be read. htslib's
/// `set_reference` loads the reference THROUGH this `.fai` and does NOT build
/// a missing one, so its absence is a hard error here rather than something to
/// paper over.
fn fai_contig_names(ref_path: &Path) -> Result<std::collections::HashSet<String>> {
    let fai = PathBuf::from(format!("{}.fai", ref_path.display()));
    ensure!(
        fai.exists(),
        "CRAM reference {} has no .fai index ({}); run `samtools faidx {}`",
        ref_path.display(),
        fai.display(),
        ref_path.display()
    );
    let text =
        std::fs::read_to_string(&fai).with_context(|| format!("reading {}", fai.display()))?;
    Ok(text
        .lines()
        .filter_map(|line| line.split('\t').next())
        .filter(|name| !name.is_empty())
        .map(str::to_owned)
        .collect())
}

/// Parses the `SN:` (sequence name) values of every `@SQ` line in `header`,
/// using the same raw byte-line scan [`Alignments::sort_order`] uses on the
/// already-parsed header text (no I/O).
fn header_sq_names(header: &bam::HeaderView) -> Vec<String> {
    let mut out = Vec::new();
    for line in header.as_bytes().split(|&b| b == b'\n') {
        if !line.starts_with(b"@SQ") {
            continue;
        }
        for field in line.split(|&b| b == b'\t') {
            if let Some(v) = field.strip_prefix(b"SN:") {
                out.push(String::from_utf8_lossy(v).into_owned());
            }
        }
    }
    out
}

/// Binds an explicit CRAM decoding reference onto `reader`, gated by a HARD
/// no-network preflight. When `reference` is `Some(path)`:
///
/// 1. `<path>.fai` must exist (see [`fai_contig_names`]).
/// 2. Every `@SQ SN:` contig in `header` must appear in the `.fai`'s contig
///    set; if ANY is absent, this errors instead of letting htslib fall back
///    to its default CRAM reference chain -- which can resolve an `M5`
///    checksum against the EBI `REF_PATH` network endpoint. THIS coverage
///    check is the no-network guarantee (the crate is `#![forbid(unsafe_code)]`
///    and edition 2024, so `std::env::set_var` to neutralize `REF_PATH` is not
///    available; the guarantee is enforced purely here).
/// 3. `reader.set_reference(path)` is called (which appends `.fai` and loads
///    the reference through that index).
///
/// A `None` reference is a no-op (BAM/SAM carry their own sequence, and a CRAM
/// opened without an explicit reference is the pre-existing behavior).
///
/// # Errors
///
/// Returns an error if the `.fai` is missing/unreadable, if any header `@SQ`
/// contig is not covered by the reference, or if `set_reference` itself fails.
fn apply_reference(
    reader: &mut bam::Reader,
    header: &bam::HeaderView,
    reference: Option<&Path>,
) -> Result<()> {
    if let Some(path) = reference {
        let names = fai_contig_names(path)?;
        for sn in header_sq_names(header) {
            ensure!(
                names.contains(&sn),
                "CRAM reference {} does not cover @SQ contig `{sn}` (no REF_PATH/network fallback \
                 is used); pass the matching reference genome via -r",
                path.display()
            );
        }
        reader
            .set_reference(path)
            .with_context(|| format!("setting CRAM reference {}", path.display()))?;
    }
    Ok(())
}

/// Decodes an htslib nt16 code the way `GetReadSeq`'s forward-strand branch
/// does (`alignments.hpp:532-544`): only `1`/`2`/`4`/`8` map to `A`/`C`/`G`/`T`;
/// every other code (including real IUPAC ambiguity codes) maps to `N`.
fn decode_forward(nt16: u8) -> u8 {
    match nt16 {
        1 => b'A',
        2 => b'C',
        4 => b'G',
        8 => b'T',
        _ => b'N',
    }
}

/// Decodes an htslib nt16 code AND complements it, the way `GetReadSeq`'s
/// reverse-strand branch does (`alignments.hpp:546-561`): `1`(A)->`T`,
/// `2`(C)->`G`, `4`(G)->`C`, `8`(T)->`A`; every other code maps to `N`
/// (complementing an already-ambiguous/unknown base is still just `N`, same
/// as the forward path).
fn decode_reverse_complement(nt16: u8) -> u8 {
    match nt16 {
        1 => b'T',
        2 => b'G',
        4 => b'C',
        8 => b'A',
        _ => b'N',
    }
}

/// Splits a CIGAR into reference-coordinate [`Segment`]s the way
/// `Alignments::Next`'s inline CIGAR walk does (`alignments.hpp:223-287`).
///
/// # Per-op reference-length contribution
///
/// - `Match`(M)/`Del`(D)/`Equal`(=)/`Diff`(X)/any-unrecognized-op: adds the
///   op's length to the pending segment's running `len` (`M`/`D` are
///   explicit `case`s; `=`/`X` fall through the switch's `default:` arm,
///   which ALSO adds -- see `alignments.hpp:271-272` -- since T1K's switch
///   only special-cases `M`/`D` explicitly, not the full "reference-
///   consuming" op set. This is not a bug in the port: it is a faithful
///   reproduction of what the vendored `switch` actually does for any op
///   code it does not explicitly list, including future/exotic ones).
/// - `Ins`(I)/`SoftClip`(S)/`HardClip`(H)/`Pad`(P): contributes `0` (S/H/P
///   fall through into I's `num = 0;` body via the switch's intentional lack
///   of `break` between the clip cases and `case Ins`,
///   `alignments.hpp:247-261`).
/// - `RefSkip`(N): flushes the current pending segment as `[start, start +
///   len - 1]`, advances `start` past both the flushed length AND the skip
///   itself (`start = start + len + num`), and resets `len` to `0` --
///   `alignments.hpp:262-270`.
///
/// After the CIGAR walk, if any pending (non-zero) `len` remains, it is
/// flushed as one final segment (`alignments.hpp:281-287`) -- this is what
/// makes a simple (no-`N`) CIGAR produce exactly ONE segment spanning the
/// whole alignment, and a spliced (`N`-containing) CIGAR produce one segment
/// per exon block.
fn cigar_to_segments(
    start_pos: i64,
    cigar: &rust_htslib::bam::record::CigarStringView,
) -> Vec<Segment> {
    let mut segments = Vec::new();
    let mut start = start_pos;
    let mut len: i64 = 0;

    for op in cigar {
        match *op {
            Cigar::Match(n) | Cigar::Del(n) => {
                len += i64::from(n);
            }
            Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) | Cigar::Ins(_) => {
                // alignments.hpp:247-261: clip ops fall through into Ins's
                // `num = 0;` -- all four contribute nothing to `len`.
            }
            Cigar::RefSkip(n) => {
                segments.push(Segment { a: start, b: start + len - 1 });
                start = start + len + i64::from(n);
                len = 0;
            }
            Cigar::Equal(n) | Cigar::Diff(n) => {
                // alignments.hpp:271-272: falls to the switch's `default:`
                // arm, same as Match/Del.
                len += i64::from(n);
            }
        }
    }

    if len > 0 {
        segments.push(Segment { a: start, b: start + len - 1 });
    }

    segments
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    fn cigar_view(ops: &[Cigar], pos: i64) -> rust_htslib::bam::record::CigarStringView {
        CigarString(ops.to_vec()).into_view(pos)
    }

    #[test]
    fn simple_match_produces_one_segment() {
        let cigar = cigar_view(&[Cigar::Match(100)], 1000);
        let segs = cigar_to_segments(1000, &cigar);
        assert_eq!(segs, vec![Segment { a: 1000, b: 1099 }]);
    }

    #[test]
    fn match_with_deletion_extends_single_segment() {
        // 50M10D50M: deletion is reference-consuming but does not split
        // into a new segment (only N does that).
        let cigar = cigar_view(&[Cigar::Match(50), Cigar::Del(10), Cigar::Match(50)], 0);
        let segs = cigar_to_segments(0, &cigar);
        assert_eq!(segs, vec![Segment { a: 0, b: 109 }]);
    }

    #[test]
    fn ref_skip_splits_into_two_segments() {
        // 50M1000N50M: spliced alignment, exactly the BamExtractor.cpp
        // scenario (RNA-seq reads spanning an intron).
        let cigar = cigar_view(&[Cigar::Match(50), Cigar::RefSkip(1000), Cigar::Match(50)], 100);
        let segs = cigar_to_segments(100, &cigar);
        assert_eq!(segs, vec![Segment { a: 100, b: 149 }, Segment { a: 1150, b: 1199 }]);
    }

    #[test]
    fn multiple_ref_skips_produce_multiple_segments() {
        // 30M100N30M100N30M: three exons.
        let cigar = cigar_view(
            &[
                Cigar::Match(30),
                Cigar::RefSkip(100),
                Cigar::Match(30),
                Cigar::RefSkip(100),
                Cigar::Match(30),
            ],
            0,
        );
        let segs = cigar_to_segments(0, &cigar);
        assert_eq!(
            segs,
            vec![Segment { a: 0, b: 29 }, Segment { a: 130, b: 159 }, Segment { a: 260, b: 289 },]
        );
    }

    #[test]
    fn soft_clip_at_head_and_tail_does_not_extend_segment() {
        // 10S80M10S: clips contribute 0 to the reference span.
        let cigar = cigar_view(&[Cigar::SoftClip(10), Cigar::Match(80), Cigar::SoftClip(10)], 500);
        let segs = cigar_to_segments(500, &cigar);
        assert_eq!(segs, vec![Segment { a: 500, b: 579 }]);
    }

    #[test]
    fn hard_clip_and_pad_do_not_extend_segment() {
        let cigar =
            cigar_view(&[Cigar::HardClip(5), Cigar::Match(40), Cigar::Pad(3), Cigar::Match(40)], 0);
        let segs = cigar_to_segments(0, &cigar);
        // Pad contributes 0 but does NOT split (only RefSkip splits), so
        // the two Match runs merge into one segment.
        assert_eq!(segs, vec![Segment { a: 0, b: 79 }]);
    }

    #[test]
    fn insertion_does_not_extend_reference_span() {
        // 40M5I40M: insertion is not reference-consuming.
        let cigar = cigar_view(&[Cigar::Match(40), Cigar::Ins(5), Cigar::Match(40)], 0);
        let segs = cigar_to_segments(0, &cigar);
        assert_eq!(segs, vec![Segment { a: 0, b: 79 }]);
    }

    #[test]
    fn equal_and_diff_ops_extend_segment_like_match() {
        // 30=5X30=: '=' and 'X' both fall to the switch's default arm,
        // which adds to len just like M/D.
        let cigar = cigar_view(&[Cigar::Equal(30), Cigar::Diff(5), Cigar::Equal(30)], 10);
        let segs = cigar_to_segments(10, &cigar);
        assert_eq!(segs, vec![Segment { a: 10, b: 74 }]);
    }

    #[test]
    fn empty_cigar_produces_no_segments() {
        let cigar = cigar_view(&[], 0);
        let segs = cigar_to_segments(0, &cigar);
        assert!(segs.is_empty());
    }

    // -- RC/qual-reverse hand-built-record tests -----------------------

    fn make_record(seq: &[u8], qual: &[u8], reverse: bool) -> bam::Record {
        let mut record = bam::Record::new();
        let cigar = CigarString(vec![Cigar::Match(u32::try_from(seq.len()).unwrap())]);
        record.set(b"read1", Some(&cigar), seq, qual);
        record.set_pos(0);
        record.set_tid(0);
        if reverse {
            record.set_reverse();
        } else {
            record.unset_reverse();
        }
        record
    }

    /// Directly exercises the forward-strand decode path (no RC) using the
    /// same decode function [`Alignments::read_seq`] uses, hand-verifying
    /// against a manually RC'd expectation to prove the orientation logic
    /// independent of any file I/O.
    #[test]
    fn forward_strand_seq_and_qual_are_unmodified() {
        let seq = b"ACGTACGT";
        let qual = [10u8, 20, 30, 40, 5, 15, 25, 35];
        let record = make_record(seq, &qual, false);

        let decoded: Vec<u8> =
            (0..record.seq().len()).map(|i| decode_forward(record.seq().encoded_base(i))).collect();
        assert_eq!(decoded, seq);

        let qual_out: Vec<u8> = record.qual().iter().map(|&q| q + 33).collect();
        let expected: Vec<u8> = qual.iter().map(|&q| q + 33).collect();
        assert_eq!(qual_out, expected);
    }

    /// Directly exercises the reverse-strand path: bases must be
    /// reverse-COMPLEMENTED (not just reversed), and quals must be
    /// REVERSED ONLY (no value transformation), tracking the same
    /// position-for-position reversal as the bases.
    #[test]
    fn reverse_strand_seq_is_reverse_complemented_and_qual_is_reversed_only() {
        let seq = b"ACGTACGT"; // as stored in BAM (reference-forward)
        let qual = [10u8, 20, 30, 40, 5, 15, 25, 35];
        let record = make_record(seq, &qual, true);

        let len = record.seq().len();
        let rc_decoded: Vec<u8> = (0..len)
            .map(|i| {
                let j = len - 1 - i;
                decode_reverse_complement(record.seq().encoded_base(j))
            })
            .collect();
        // "ACGTACGT" reverse-complemented is "ACGTACGT" is a palindrome
        // coincidentally at this length/composition; use an asymmetric
        // sequence to make the assertion meaningful.
        assert_eq!(rc_decoded, b"ACGTACGT"); // sanity: RC of this specific palindrome-ish seq

        let asym_seq = b"AACCGGTT";
        let asym_record = make_record(asym_seq, &qual, true);
        let asym_len = asym_record.seq().len();
        let asym_rc: Vec<u8> = (0..asym_len)
            .map(|i| {
                let j = asym_len - 1 - i;
                decode_reverse_complement(asym_record.seq().encoded_base(j))
            })
            .collect();
        // Reverse of "AACCGGTT" is "TTGGCCAA"; complement each base:
        // T->A, T->A, G->C, G->C, C->G, C->G, A->T, A->T = "AACCGGTT"...
        // compute directly instead of hand-deriving to avoid a transcription error:
        let mut expected = asym_seq.to_vec();
        expected.reverse();
        for b in &mut expected {
            *b = match *b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => other,
            };
        }
        assert_eq!(asym_rc, expected);

        // Qual: reversed ONLY, no value change.
        let qual_out: Vec<u8> = (0..len)
            .map(|i| {
                let j = len - 1 - i;
                asym_record.qual()[j] + 33
            })
            .collect();
        let mut expected_qual: Vec<u8> = qual.iter().map(|&q| q + 33).collect();
        expected_qual.reverse();
        assert_eq!(qual_out, expected_qual);
    }

    #[test]
    fn non_acgt_nt16_codes_decode_to_n_not_iupac_ambiguity() {
        // htslib's own Seq::as_bytes() would decode nt16 code 3 ('M', an
        // ambiguity code) faithfully; T1K's GetReadSeq collapses it (and
        // every other non-A/C/G/T code) to 'N'. Directly test the decode
        // functions against every non-canonical nt16 value (0-15 minus
        // 1/2/4/8).
        for code in 0u8..16 {
            if [1, 2, 4, 8].contains(&code) {
                continue;
            }
            assert_eq!(decode_forward(code), b'N', "nt16 code {code} should decode to N");
            assert_eq!(
                decode_reverse_complement(code),
                b'N',
                "nt16 code {code} should RC-decode to N"
            );
        }
    }

    #[test]
    fn sort_order_reads_hd_so_and_go() {
        use rust_htslib::bam::header::HeaderRecord;
        use rust_htslib::bam::{Format, Header, Writer};

        fn write_header_bam(dir: &std::path::Path, hd_tags: &[(&str, &str)]) -> std::path::PathBuf {
            let path = dir.join("hdr.bam");
            let mut header = Header::new();
            let mut hd = HeaderRecord::new(b"HD");
            hd.push_tag(b"VN", "1.6");
            for (k, v) in hd_tags {
                hd.push_tag(k.as_bytes(), v);
            }
            header.push_record(&hd);
            let mut sq = HeaderRecord::new(b"SQ");
            sq.push_tag(b"SN", "chr1");
            sq.push_tag(b"LN", 1000);
            header.push_record(&sq);
            // Writer with no records: header-only BAM is enough for sort_order.
            let _writer = Writer::from_path(&path, &header, Format::Bam).unwrap();
            path
        }

        let dir = tempfile::tempdir().unwrap();

        let coord = write_header_bam(dir.path(), &[("SO", "coordinate")]);
        assert_eq!(Alignments::open(&coord).unwrap().sort_order(), SortOrder::Coordinate);

        let qname = write_header_bam(dir.path(), &[("SO", "queryname")]);
        assert_eq!(Alignments::open(&qname).unwrap().sort_order(), SortOrder::QueryName);

        let grouped = write_header_bam(dir.path(), &[("GO", "query")]);
        assert_eq!(Alignments::open(&grouped).unwrap().sort_order(), SortOrder::QueryGrouped);

        let unsorted = write_header_bam(dir.path(), &[("SO", "unsorted")]);
        assert_eq!(Alignments::open(&unsorted).unwrap().sort_order(), SortOrder::Unsorted);

        // No SO/GO tag at all → Unsorted.
        let bare = write_header_bam(dir.path(), &[]);
        assert_eq!(Alignments::open(&bare).unwrap().sort_order(), SortOrder::Unsorted);
    }

    /// Builds a small coordinate-sorted paired BAM: 4 pairs, read-1/read-2
    /// lengths 100bp/80bp respectively (so a read-1-only sampler must ignore
    /// the 80bp mates), plus one secondary alignment (flag `0x100`) of an
    /// existing read-1 that must NOT be double-counted (mirrors
    /// `general_info`'s own primary-only gate).
    fn write_paired_test_bam(dir: &std::path::Path) -> std::path::PathBuf {
        use rust_htslib::bam::header::HeaderRecord;
        use rust_htslib::bam::record::{Cigar, CigarString};
        use rust_htslib::bam::{Format, Header, Writer};

        let path = dir.join("sample_read1.bam");
        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);

        let mut writer = Writer::from_path(&path, &header, Format::Bam).unwrap();
        let seq100 = vec![b'A'; 100];
        let seq80 = vec![b'C'; 80];
        let qual100 = vec![30u8; 100];
        let qual80 = vec![30u8; 80];

        for i in 0..4i64 {
            let mut r1 = bam::Record::new();
            r1.set(
                format!("pair{i}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(100)])),
                &seq100,
                &qual100,
            );
            r1.set_tid(0);
            r1.set_pos(1000 + i * 200);
            r1.set_mtid(0);
            r1.set_mpos(1100 + i * 200);
            r1.set_flags(0x1 | 0x2 | 0x40);
            writer.write(&r1).unwrap();

            let mut r2 = bam::Record::new();
            r2.set(
                format!("pair{i}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(80)])),
                &seq80,
                &qual80,
            );
            r2.set_tid(0);
            r2.set_pos(1100 + i * 200);
            r2.set_mtid(0);
            r2.set_mpos(1000 + i * 200);
            r2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
            writer.write(&r2).unwrap();
        }

        // A secondary alignment of pair0's read-1 -- must be skipped (not
        // primary), so it must not inflate the sampled sum/count.
        let mut secondary = bam::Record::new();
        secondary.set(b"pair0", Some(&CigarString(vec![Cigar::Match(100)])), &seq100, &qual100);
        secondary.set_tid(0);
        secondary.set_pos(5000);
        secondary.set_mtid(0);
        secondary.set_mpos(5100);
        secondary.set_flags(0x1 | 0x2 | 0x40 | 0x100);
        writer.write(&secondary).unwrap();

        drop(writer);
        path
    }

    #[test]
    fn sample_read1_len_sum_sums_only_primary_read1_records() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_paired_test_bam(dir.path());
        let mut alignments = Alignments::open(&path).unwrap();

        let (sum, count) = alignments.sample_read1_len_sum(1000).unwrap();

        // 4 read-1 records at 100bp each; the 4 read-2 (80bp) records and
        // the secondary read-1 alignment must be excluded.
        assert_eq!(count, 4);
        assert_eq!(sum, 400);
    }

    #[test]
    fn sample_read1_len_sum_respects_limit() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_paired_test_bam(dir.path());
        let mut alignments = Alignments::open(&path).unwrap();

        // limit=2 sites (any record, primary or not, counts against the
        // limit): the first 2 records scanned are pair0's read-1 (100bp)
        // and pair0's read-2 (80bp, not sampled since not first-mate).
        let (sum, count) = alignments.sample_read1_len_sum(2).unwrap();
        assert_eq!(count, 1);
        assert_eq!(sum, 100);
    }

    #[test]
    fn sample_read1_len_sum_on_single_end_samples_every_primary_record() {
        use rust_htslib::bam::header::HeaderRecord;
        use rust_htslib::bam::record::{Cigar, CigarString};
        use rust_htslib::bam::{Format, Header, Writer};

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("single_end.bam");

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1_000_000);
        header.push_record(&sq);

        {
            let mut writer = Writer::from_path(&path, &header, Format::Bam).unwrap();
            for i in 0..3i64 {
                let seq = vec![b'A'; 60];
                let qual = vec![30u8; 60];
                let mut r = bam::Record::new();
                r.set(
                    format!("se{i}").as_bytes(),
                    Some(&CigarString(vec![Cigar::Match(60)])),
                    &seq,
                    &qual,
                );
                r.set_tid(0);
                r.set_pos(1000 + i * 100);
                r.set_mtid(-1);
                r.set_mpos(-1);
                r.set_flags(0);
                writer.write(&r).unwrap();
            }
        }

        let mut alignments = Alignments::open(&path).unwrap();
        let (sum, count) = alignments.sample_read1_len_sum(1000).unwrap();
        assert_eq!(count, 3);
        assert_eq!(sum, 180);
    }

    /// Covers the loop's post-scan fallback (`sort_order`'s trailing
    /// `SortOrder::Unsorted` after the `for line in text.split(...)` loop),
    /// which is distinct from the in-loop `_ => Unsorted` arm that
    /// `sort_order_reads_hd_so_and_go`'s "bare" case already exercises for a
    /// header that HAS an `@HD` line but no `SO`/`GO` tag. Here the header has
    /// NO `@HD` record at all -- only `@SQ` -- so the loop's `if
    /// !line.starts_with(b"@HD")` guard skips every line and control falls
    /// through to the function's final `SortOrder::Unsorted` statement.
    #[test]
    fn sort_order_no_hd_line_is_unsorted() {
        use rust_htslib::bam::header::HeaderRecord;
        use rust_htslib::bam::{Format, Header, Writer};

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("no_hd.bam");

        // Deliberately do NOT push an `@HD` record -- only `@SQ`.
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 1000);
        header.push_record(&sq);

        // A header-only (zero-record) BAM without `@HD` is not reliably
        // readable on every htslib version, so write one minimal unmapped
        // record after the header to keep the file unambiguously valid.
        {
            let mut writer = Writer::from_path(&path, &header, Format::Bam).unwrap();
            let mut record = bam::Record::new();
            record.set(b"read1", None, b"ACGT", &[30u8, 30, 30, 30]);
            writer.write(&record).unwrap();
        }

        assert_eq!(Alignments::open(&path).unwrap().sort_order(), SortOrder::Unsorted);
    }

    // -- CRAM explicit-reference wiring --------------------------------

    /// Writes a single-contig FASTA (`>{contig}\n{seq}\n`) plus a hand-built
    /// `.fai` sibling. The `.fai` columns are NAME, LENGTH, OFFSET (byte
    /// offset of the first base = `>` + contig + `\n` = `contig.len() + 2`),
    /// LINEBASES (`seq.len()`, single-line sequence) and LINEWIDTH
    /// (`seq.len() + 1`, bases plus the trailing newline). Only the NAME
    /// column is consulted by the preflight, but the offset/linebase columns
    /// are written correctly so htslib can actually load the reference through
    /// the index (`set_reference` appends `.fai` and reads THROUGH it; it does
    /// not build a missing index).
    fn write_fasta_with_fai(ref_fa: &std::path::Path, contig: &str, seq: &[u8]) {
        let mut fasta = Vec::new();
        fasta.push(b'>');
        fasta.extend_from_slice(contig.as_bytes());
        fasta.push(b'\n');
        fasta.extend_from_slice(seq);
        fasta.push(b'\n');
        std::fs::write(ref_fa, &fasta).unwrap();

        let offset = contig.len() + 2; // ">" + contig + "\n"
        let linebases = seq.len();
        let linewidth = seq.len() + 1; // bases + "\n"
        let fai_path = std::path::PathBuf::from(format!("{}.fai", ref_fa.display()));
        let fai_line = format!("{contig}\t{}\t{offset}\t{linebases}\t{linewidth}\n", seq.len());
        std::fs::write(&fai_path, fai_line).unwrap();
    }

    /// Builds a CRAM with the given single-contig header and two aligned
    /// 100bp reads (at reference positions 0 and 100), reference-compressed
    /// against `ref_fa` (which must already have its `.fai` sibling written).
    fn write_cram(cram: &std::path::Path, ref_fa: &std::path::Path, contig: &str, ref_seq: &[u8]) {
        use rust_htslib::bam::header::HeaderRecord;
        use rust_htslib::bam::record::{Cigar, CigarString};
        use rust_htslib::bam::{Format, Header, Writer};

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        header.push_record(&hd);
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", contig);
        sq.push_tag(b"LN", i64::try_from(ref_seq.len()).unwrap());
        header.push_record(&sq);

        let mut writer = Writer::from_path(cram, &header, Format::Cram).unwrap();
        writer.set_reference(ref_fa).expect("set CRAM writer reference");

        let qual = vec![30u8; 100];
        for (i, pos) in [0i64, 100].into_iter().enumerate() {
            // Read bases equal to the reference slice they align to, so the
            // CRAM reference-vs-read features are trivial and round-trip.
            let start = usize::try_from(pos).unwrap();
            let seq = &ref_seq[start..start + 100];
            let mut r = bam::Record::new();
            r.set(
                format!("read{i}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(100)])),
                seq,
                &qual,
            );
            r.set_tid(0);
            r.set_pos(pos);
            r.set_flags(0); // mapped, unpaired
            writer.write(&r).unwrap();
        }
        drop(writer);
    }

    /// Writes `ref.fa` + `ref.fa.fai` (contig `chr1`, 200bp) and a CRAM of two
    /// aligned reads against it. The `.fai` is written by hand as
    /// `chr1\t200\t6\t200\t201` (offset 6 = `>chr1\n`; 200 bases on one line;
    /// linewidth 201).
    fn write_ref_fai_and_cram(ref_fa: &std::path::Path, cram: &std::path::Path) {
        let ref_seq: Vec<u8> = b"ACGT".iter().copied().cycle().take(200).collect();
        write_fasta_with_fai(ref_fa, "chr1", &ref_seq);
        write_cram(cram, ref_fa, "chr1", &ref_seq);
    }

    /// Builds a CRAM whose `@SQ` is `chr1` but returns a DIFFERENT reference
    /// (contig `chr2`) whose `.fai` does not cover `chr1`, for the no-network
    /// preflight test. The CRAM itself is built against a correct `chr1`
    /// reference (so writing succeeds); the returned reference is the
    /// mismatched one the preflight must reject.
    fn write_cram_and_mismatched_reference(
        dir: &std::path::Path,
    ) -> (std::path::PathBuf, std::path::PathBuf) {
        let build_ref = dir.join("build_ref.fa");
        let cram = dir.join("x.cram");
        write_ref_fai_and_cram(&build_ref, &cram);

        // A reference covering only `chr2` -- its `.fai` contig set does not
        // contain the CRAM's `@SQ SN:chr1`.
        let mismatch_ref = dir.join("mismatch_ref.fa");
        let mismatch_seq: Vec<u8> = b"TTTT".iter().copied().cycle().take(200).collect();
        write_fasta_with_fai(&mismatch_ref, "chr2", &mismatch_seq);

        (mismatch_ref, cram)
    }

    #[test]
    fn cram_with_reference_reads_and_survives_rewind() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_fa = tmp.path().join("ref.fa");
        let cram = tmp.path().join("x.cram");
        write_ref_fai_and_cram(&ref_fa, &cram);

        let mut a = Alignments::open_with_reference(&cram, Some(&ref_fa)).expect("open CRAM");
        let mut n = 0;
        while a.next().unwrap() {
            n += 1;
        }
        assert_eq!(n, 2);

        a.rewind().expect("rewind");
        let mut m = 0;
        while a.next().unwrap() {
            m += 1;
        }
        assert_eq!(m, 2);
    }

    #[test]
    fn cram_reference_missing_contig_errors_without_network() {
        // A ref.fa/.fai whose contig set does NOT cover the CRAM's @SQ must
        // error at open (preflight), never attempting htslib's M5/REF_PATH
        // network fallback.
        let tmp = tempfile::tempdir().unwrap();
        let (ref_fa, cram) = write_cram_and_mismatched_reference(tmp.path());
        // `.map(|_| ())` discards the `Ok(Alignments)` (which is not `Debug`)
        // so `unwrap_err` can format the error path on failure.
        let err = Alignments::open_with_reference(&cram, Some(&ref_fa))
            .map(|_| ())
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("reference") && (err.contains("contig") || err.contains("@SQ")),
            "unexpected error message: {err}"
        );
    }
}
