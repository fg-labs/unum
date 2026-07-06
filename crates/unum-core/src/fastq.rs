//! Minimal FASTA/FASTQ reader with kseq-parity semantics, ported to cover
//! exactly what [`crate::extract`] needs from T1K's `ReadFiles`/`kseq.h`
//! (`vendor/t1k/ReadFiles.hpp`, `vendor/t1k/kseq.h`).
//!
//! # `id` is the header token up to the first whitespace; the rest is a
//! dropped comment
//!
//! `kseq_read` (`kseq.h:187-188`) reads the header name via
//! `ks_getuntil(ks, 0, &seq->name, &c)` -- separator mode `0` is
//! `KS_SEP_SPACE` (`kseq.h:35`: "isspace(): \t, \n, \v, \f, \r" -- i.e. any
//! whitespace, matching `str::split_whitespace`'s definition, NOT just the
//! ASCII space character). Whatever follows on the same header line (after
//! that first whitespace run) is read separately into `seq->comment`
//! (`kseq.h:196`) and never appended back onto `name`/`id` on the read path
//! this module exercises (`ReadFiles::Next`, `ReadFiles.hpp:183`, uses
//! `strdup(inSeq[...]->name.s)` directly -- the commented-out block at
//! `ReadFiles.hpp:168-173` that WOULD reattach the comment is dead code,
//! never compiled in). So `id` here is exactly the first whitespace-delimited
//! token of the header line, and everything after it is parsed but
//! discarded.
//!
//! # Trailing `/1`/`/2` mate-suffix stripping
//!
//! `ReadFiles::Next` (`ReadFiles.hpp:184-189`) strips a trailing `/1` or `/2`
//! from `id` (checking the last two characters are `/` followed by `1` or
//! `2`) after reading it from kseq. [`FastqReader::next_record`] reproduces this
//! exactly.
//!
//! # Multi-line sequence/quality wrapping
//!
//! `kseq_read`'s body loop (`kseq.h:203-204`) concatenates every non-header,
//! non-`+`-line line into `seq`/`qual` with no separator -- i.e. arbitrary
//! line wrapping is supported and collapsed away. This reader does the same.
//!
//! # Format detection: FASTA (`>`) vs. FASTQ (`@`)
//!
//! kseq dispatches purely on which sentinel character (`>` or `@`) starts
//! the record (`kseq.h:190-191`), not on file extension. This reader mirrors
//! that: the first non-whitespace byte of the (decompressed) stream
//! determines the format for the whole file, matching `ReadFiles::AddReadFile`
//! setting a per-file `FILE_TYPE` once from the first `kseq_read` call
//! (`ReadFiles.hpp:94-107`).
//!
//! # Gzip transparency
//!
//! T1K opens every read file via `gzopen`/`gzread` (`ReadFiles.hpp:9,90-91`),
//! which transparently reads both plain and gzip-compressed input (zlib
//! auto-detects the gzip magic bytes). This reader mirrors that by sniffing
//! the gzip magic (`0x1f 0x8b`) at the start of the file and wrapping in a
//! [`flate2::read::MultiGzDecoder`] only when present, otherwise reading the
//! file directly -- functionally equivalent auto-detection, without needing
//! to inspect the file extension.

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

/// A single parsed FASTA/FASTQ record: `id` is the header token up to the
/// first whitespace (see module docs), `seq` is the concatenated sequence
/// (no line wrapping), and `qual` is `Some(quality)` for FASTQ records or
/// `None` for FASTA records (mirroring `ReadFiles::qual` being `NULL` when
/// `inSeq[...]->qual.l == 0`, `ReadFiles.hpp:192-195`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

/// A rewindable FASTA/FASTQ reader over a single file path, supporting
/// plain-text and gzip-compressed input (see module docs). Mirrors the
/// single-file slice of T1K's `ReadFiles` that [`crate::extract`] needs:
/// sequential [`FastqReader::next_record`] plus [`FastqReader::rewind`] (re-open
/// from the start, matching `gzrewind`/`kseq_rewind`,
/// `ReadFiles.hpp:132-140`).
///
/// Does not support multiple concatenated input files or interleaved input
/// (`ReadFiles::AddReadFile`'s `fileHasMate`/`interleavedId` generality) --
/// [`crate::extract`] only ever needs one file per mate, and `-i`
/// interleaved input is explicitly deferred (see [`crate::extract`]'s module
/// docs).
pub struct FastqReader {
    path: PathBuf,
    inner: Box<dyn BufRead>,
}

impl FastqReader {
    /// Opens `path` for reading, auto-detecting gzip compression (see module
    /// docs).
    ///
    /// # Errors
    ///
    /// Returns an error if `path` cannot be opened.
    pub fn open(path: &Path) -> Result<Self> {
        let inner = open_reader(path)?;
        Ok(Self { path: path.to_path_buf(), inner })
    }

    /// Re-opens the file from the beginning, matching `ReadFiles::Rewind`'s
    /// `gzrewind`/`kseq_rewind` (`ReadFiles.hpp:132-140`) -- functionally a
    /// fresh read position at byte 0, which for a gzip stream also means a
    /// fresh decompressor state (gzip files cannot be seeked-and-resumed
    /// mid-stream the way plain files can, so re-opening is the correct
    /// general-purpose equivalent, not merely a convenient one).
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be re-opened.
    pub fn rewind(&mut self) -> Result<()> {
        self.inner = open_reader(&self.path)?;
        Ok(())
    }

    /// Reads the next record, matching `ReadFiles::Next` (`ReadFiles.hpp:
    /// 155-204`): the header token becomes `id` (with a trailing `/1`/`/2`
    /// mate suffix stripped, see module docs), the concatenated
    /// (possibly-multi-line) sequence becomes `seq`, and (for FASTQ input)
    /// the concatenated quality string becomes `qual`. Returns `Ok(None)` at
    /// end of file (matching `Next`'s `return 0`).
    ///
    /// # Errors
    ///
    /// Returns an error on malformed input (e.g. a FASTQ record whose
    /// quality string length does not match its sequence length, mirroring
    /// `kseq_read`'s `-2` return, or a record with no header sentinel).
    ///
    /// # Panics
    ///
    /// Does not panic: the internal `expect` in the sequence-line loop is
    /// only ever reached immediately after [`peek_line`] has already
    /// confirmed a line is present, so the subsequent
    /// [`read_line_trimmed`] call in the same iteration cannot return
    /// `None`.
    pub fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        let Some(mut header) = read_line_trimmed(&mut self.inner)? else {
            return Ok(None);
        };
        loop {
            if header.is_empty() {
                let Some(next_header) = read_line_trimmed(&mut self.inner)? else {
                    return Ok(None);
                };
                header = next_header;
                continue;
            }
            break;
        }

        let is_fastq = match header.as_bytes().first() {
            Some(b'@') => true,
            Some(b'>') => false,
            _ => anyhow::bail!(
                "malformed record in {}: expected '@' or '>' header, got {header:?}",
                self.path.display()
            ),
        };
        let header_body = &header[1..];
        // kseq's `KS_SEP_SPACE` splits on ANY whitespace, not just ' ' (see
        // module docs) -- `split_whitespace` matches that exactly.
        let mut id = header_body.split_whitespace().next().unwrap_or("").to_string();
        strip_mate_suffix(&mut id);

        let mut seq: Vec<u8> = Vec::new();
        let mut qual: Option<Vec<u8>> = None;

        // Sequence lines, until a '+' (FASTQ separator) or the next header
        // sentinel ('@'/'>') or EOF.
        loop {
            let Some(peeked) = peek_line(&mut self.inner)? else { break };
            if peeked.is_empty() {
                // Blank line: consume and skip (kseq.h:200: "skip empty
                // lines").
                let _ = read_line_trimmed(&mut self.inner)?;
                continue;
            }
            let first = peeked.as_bytes()[0];
            if first == b'+' || (is_fastq && (first == b'@')) || (!is_fastq && first == b'>') {
                break;
            }
            let line = read_line_trimmed(&mut self.inner)?.expect("just peeked Some");
            seq.extend_from_slice(line.as_bytes());
        }

        if is_fastq {
            let Some(plus_line) = read_line_trimmed(&mut self.inner)? else {
                anyhow::bail!(
                    "truncated FASTQ record in {} (missing '+' line for {id})",
                    self.path.display()
                );
            };
            anyhow::ensure!(
                plus_line.starts_with('+'),
                "malformed FASTQ record in {}: expected '+' separator for {id}, got {plus_line:?}",
                self.path.display()
            );

            let mut q: Vec<u8> = Vec::new();
            while q.len() < seq.len() {
                let Some(line) = read_line_trimmed(&mut self.inner)? else {
                    anyhow::bail!(
                        "truncated quality string in {} for read {id}",
                        self.path.display()
                    );
                };
                q.extend_from_slice(line.as_bytes());
            }
            anyhow::ensure!(
                q.len() == seq.len(),
                "quality/sequence length mismatch in {} for read {id} (seq={}, qual={})",
                self.path.display(),
                seq.len(),
                q.len()
            );
            qual = Some(q);
        }

        Ok(Some(FastqRecord { id, seq, qual }))
    }
}

/// Strips a trailing `/1` or `/2` mate suffix from `id` in place, matching
/// `ReadFiles::Next`'s exact check (`ReadFiles.hpp:184-189`): the LAST
/// character must be `1` or `2` AND the SECOND-TO-LAST character must be
/// `/`.
fn strip_mate_suffix(id: &mut String) {
    let bytes = id.as_bytes();
    let len = bytes.len();
    if len >= 2 {
        let last = bytes[len - 1];
        let second_last = bytes[len - 2];
        if (last == b'1' || last == b'2') && second_last == b'/' {
            id.truncate(len - 2);
        }
    }
}

/// Opens `path`, auto-detecting gzip compression by sniffing the gzip magic
/// bytes (`0x1f 0x8b`) at the start of the file (see module docs).
fn open_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file =
        File::open(path).with_context(|| format!("opening read file {}", path.display()))?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic).with_context(|| format!("reading {}", path.display()))?;
    let is_gzip = n == 2 && magic == [0x1f, 0x8b];

    // Re-open from the start rather than trying to "un-read" the sniffed
    // bytes through the same handle -- simplest correct approach, and this
    // is only paid once per `open`/`rewind` call, not per record.
    let file = File::open(path).with_context(|| format!("opening read file {}", path.display()))?;
    if is_gzip {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Reads one line, stripping a trailing `\n` and (if present) `\r`, into
/// `String`. Returns `Ok(None)` at EOF with nothing read.
///
/// # Errors
///
/// Returns an error if the line is not valid UTF-8 (T1K's `char*`-based
/// kseq has no such restriction, but reference/read FASTA/FASTQ content in
/// this pipeline is always plain ASCII, so this is not a practical
/// limitation).
fn read_line_trimmed(r: &mut Box<dyn BufRead>) -> Result<Option<String>> {
    let mut buf = Vec::new();
    let n = r.read_until(b'\n', &mut buf).context("reading line")?;
    if n == 0 {
        return Ok(None);
    }
    if buf.last() == Some(&b'\n') {
        buf.pop();
        if buf.last() == Some(&b'\r') {
            buf.pop();
        }
    }
    let s = String::from_utf8(buf).context("read file line is not valid UTF-8")?;
    Ok(Some(s))
}

/// Peeks at the next line without consuming it (needed to decide, before
/// consuming, whether the next line starts a new record). Implemented via
/// `BufRead::fill_buf`, which fills its internal buffer with at least one
/// byte (or reports EOF) without consuming anything; since every line this
/// reader handles is short relative to the default `BufReader` capacity
/// (8 KiB), a single `fill_buf` call is sufficient to see the whole line in
/// practice for this pipeline's inputs (reference/read FASTA/FASTQ, not
/// arbitrary attacker-controlled data).
fn peek_line(r: &mut Box<dyn BufRead>) -> Result<Option<String>> {
    let buf = r.fill_buf().context("peeking line")?;
    if buf.is_empty() {
        return Ok(None);
    }
    let end = buf.iter().position(|&b| b == b'\n').unwrap_or(buf.len());
    let mut line = buf[..end].to_vec();
    if line.last() == Some(&b'\r') {
        line.pop();
    }
    let s = String::from_utf8(line).context("read file line is not valid UTF-8")?;
    Ok(Some(s))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp(contents: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn parses_id_up_to_first_whitespace_dropping_comment() {
        let f = write_temp("@read1 some comment here\nACGT\n+\nIIII\n");
        let mut r = FastqReader::open(f.path()).unwrap();
        let rec = r.next_record().unwrap().unwrap();
        assert_eq!(rec.id, "read1");
        assert_eq!(rec.seq, b"ACGT");
        assert_eq!(rec.qual, Some(b"IIII".to_vec()));
    }

    #[test]
    fn parses_id_split_on_tab_too() {
        let f = write_temp("@read1\tsome comment\nACGT\n+\nIIII\n");
        let mut r = FastqReader::open(f.path()).unwrap();
        let rec = r.next_record().unwrap().unwrap();
        assert_eq!(rec.id, "read1");
    }

    #[test]
    fn strips_trailing_slash_1_and_slash_2() {
        let f = write_temp("@read1/1\nACGT\n+\nIIII\n@read1/2\nTTTT\n+\nJJJJ\n");
        let mut r = FastqReader::open(f.path()).unwrap();
        let rec1 = r.next_record().unwrap().unwrap();
        assert_eq!(rec1.id, "read1");
        let rec2 = r.next_record().unwrap().unwrap();
        assert_eq!(rec2.id, "read1");
    }

    #[test]
    fn does_not_strip_non_mate_suffix() {
        let f = write_temp("@read_a1\nACGT\n+\nIIII\n");
        let mut r = FastqReader::open(f.path()).unwrap();
        let rec = r.next_record().unwrap().unwrap();
        assert_eq!(rec.id, "read_a1");
    }

    #[test]
    fn multiple_records_and_eof() {
        let f = write_temp("@r1\nACGT\n+\nIIII\n@r2\nGGCC\n+\nJJJJ\n");
        let mut r = FastqReader::open(f.path()).unwrap();
        assert_eq!(r.next_record().unwrap().unwrap().id, "r1");
        assert_eq!(r.next_record().unwrap().unwrap().id, "r2");
        assert!(r.next_record().unwrap().is_none());
    }

    #[test]
    fn fasta_input_has_no_qual() {
        let f = write_temp(">seq1 desc\nACGTACGT\n>seq2\nTTTT\n");
        let mut r = FastqReader::open(f.path()).unwrap();
        let rec1 = r.next_record().unwrap().unwrap();
        assert_eq!(rec1.id, "seq1");
        assert_eq!(rec1.seq, b"ACGTACGT");
        assert_eq!(rec1.qual, None);
        let rec2 = r.next_record().unwrap().unwrap();
        assert_eq!(rec2.id, "seq2");
    }

    #[test]
    fn multiline_sequence_and_quality_are_concatenated() {
        let f = write_temp("@r1\nACGT\nACGT\n+\nIIII\nJJJJ\n");
        let mut r = FastqReader::open(f.path()).unwrap();
        let rec = r.next_record().unwrap().unwrap();
        assert_eq!(rec.seq, b"ACGTACGT");
        assert_eq!(rec.qual, Some(b"IIIIJJJJ".to_vec()));
    }

    #[test]
    fn rewind_resets_to_first_record() {
        let f = write_temp("@r1\nACGT\n+\nIIII\n@r2\nGGCC\n+\nJJJJ\n");
        let mut r = FastqReader::open(f.path()).unwrap();
        assert_eq!(r.next_record().unwrap().unwrap().id, "r1");
        assert_eq!(r.next_record().unwrap().unwrap().id, "r2");
        r.rewind().unwrap();
        assert_eq!(r.next_record().unwrap().unwrap().id, "r1");
    }

    #[test]
    fn gzip_input_is_auto_detected() {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        {
            let mut enc = GzEncoder::new(&mut tmp, Compression::default());
            enc.write_all(b"@r1\nACGT\n+\nIIII\n").unwrap();
            enc.finish().unwrap();
        }
        let mut r = FastqReader::open(tmp.path()).unwrap();
        let rec = r.next_record().unwrap().unwrap();
        assert_eq!(rec.id, "r1");
        assert_eq!(rec.seq, b"ACGT");
    }

    #[test]
    fn quality_length_mismatch_is_an_error() {
        let f = write_temp("@r1\nACGT\n+\nII\n");
        let mut r = FastqReader::open(f.path()).unwrap();
        // Since qual is shorter than seq, the reader keeps reading lines
        // until it hits EOF looking for more quality bases -- this should
        // surface as an error (truncated quality), not silently succeed.
        let result = r.next_record();
        assert!(result.is_err(), "expected an error for truncated quality, got {result:?}");
    }
}
