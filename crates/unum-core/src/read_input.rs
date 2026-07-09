//! Content-based input format detection (magic bytes, never the file
//! extension) plus reader-opening helpers shared by `extract`/`run`.
//!
//! Compression is handled by niffler before detection (see
//! [`crate::fastq`]); this module distinguishes the decompressed container.

use anyhow::{Context, Result, bail};
use std::io::Read;
use std::path::PathBuf;

use crate::fastq::FastqReader;

/// The container format detected from an input's leading bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectedFormat {
    /// FASTQ (`@` sentinel).
    Fastq,
    /// FASTA (`>` sentinel).
    Fasta,
    /// SAM (a text `@HD`/`@SQ`/`@RG`/... header line: `@` + two uppercase
    /// letters + TAB). Read by htslib exactly like BAM/CRAM.
    Sam,
    /// BAM (`BAM\1` magic).
    Bam,
    /// CRAM (`CRAM` magic).
    Cram,
}

/// Where an input comes from.
#[derive(Debug, Clone)]
pub enum InputSpec {
    /// A filesystem path.
    Path(PathBuf),
    /// Standard input (`-`).
    Stdin,
}

/// The number of leading bytes needed to classify every container this
/// pipeline recognizes: the `BAM\x01` and `CRAM` magics are both 4 bytes long,
/// while FASTQ (`@`) and FASTA (`>`) are distinguished by their first byte.
/// [`read_magic_prefix`] buffers up to this many bytes before classifying.
const MAGIC_PREFIX_LEN: usize = 4;

/// The result of opening an input for extraction: either a ready FASTQ/FASTA
/// reader, or a marker that the input is a BAM/CRAM (which the caller routes
/// to htslib by its original path/stdin -- we never hand a peeked Rust stream
/// to htslib; see the stdin routing note in the design doc).
pub enum OpenedInput {
    /// FASTQ/FASTA: a reader positioned at the first record.
    Fastq(FastqReader),
    /// SAM (`@`-header text) -- route the original input to the BAM extractor
    /// (htslib reads SAM like BAM; no reference needed, same as BAM).
    Sam,
    /// BAM (`BAM\1` magic) -- route the original input to the BAM extractor.
    Bam,
    /// CRAM (`CRAM` magic) -- route the original input to the BAM extractor.
    Cram,
}

/// Classifies an input from its leading bytes (`prefix`, up to
/// [`MAGIC_PREFIX_LEN`] bytes as produced by [`read_magic_prefix`]).
///
/// # Errors
///
/// Returns an error if `prefix` is empty (empty stream) or starts with a byte
/// that is not one of the recognized signatures.
pub fn detect_format(prefix: &[u8]) -> Result<DetectedFormat> {
    if prefix.is_empty() {
        bail!("input is empty; cannot detect format");
    }
    if prefix.starts_with(b"CRAM") {
        return Ok(DetectedFormat::Cram);
    }
    if prefix.starts_with(b"BAM\x01") {
        return Ok(DetectedFormat::Bam);
    }
    // A SAM header line is `@` + a two-uppercase-letter record type (HD/SQ/RG/
    // PG/CO) + TAB. FASTQ also begins with `@`, but its id line never has this
    // exact shape (an id is not two uppercase letters followed by a tab), so
    // this reliably distinguishes a headered SAM from FASTQ. (Headerless SAM --
    // a bare record line -- is not detected here; `samtools view -h` always
    // emits a header, so this covers the practical case.)
    if prefix.len() >= 4
        && prefix[0] == b'@'
        && prefix[1].is_ascii_uppercase()
        && prefix[2].is_ascii_uppercase()
        && prefix[3] == b'\t'
    {
        return Ok(DetectedFormat::Sam);
    }
    match prefix[0] {
        b'@' => Ok(DetectedFormat::Fastq),
        b'>' => Ok(DetectedFormat::Fasta),
        other => bail!(
            "unrecognized input format: leading byte {other:#04x} is not one of \
             '@' (FASTQ), '>' (FASTA), 'BAM\\1', or 'CRAM'"
        ),
    }
}

/// Reads up to [`MAGIC_PREFIX_LEN`] leading bytes of `reader` into a `Vec`,
/// looping until the full prefix is buffered or EOF is reached. A single
/// `Read::read` may legally return fewer bytes than requested (a streamed
/// BAM/CRAM can deliver 1-3 bytes on the first read), so the magic must not be
/// matched until the full 4-byte prefix is in hand -- otherwise a short first
/// read misclassifies BAM/CRAM as unrecognized. These bytes are consumed from
/// `reader`; [`open_input`] chains them back in front of the remainder
/// so no input is lost.
///
/// # Errors
///
/// Returns an error on an I/O failure while reading.
fn read_magic_prefix(reader: &mut impl Read) -> Result<Vec<u8>> {
    let mut prefix = Vec::with_capacity(MAGIC_PREFIX_LEN);
    while prefix.len() < MAGIC_PREFIX_LEN {
        let mut chunk = [0u8; MAGIC_PREFIX_LEN];
        let want = MAGIC_PREFIX_LEN - prefix.len();
        let got = reader.read(&mut chunk[..want]).context("reading input for format detection")?;
        if got == 0 {
            break; // EOF: fewer than MAGIC_PREFIX_LEN bytes available.
        }
        prefix.extend_from_slice(&chunk[..got]);
    }
    Ok(prefix)
}

/// Opens `spec`, transparently decompresses (niffler), detects the format,
/// and returns a FASTQ reader or a BAM/CRAM marker. Never errors on BAM/CRAM
/// -- that routing is the caller's decision (`--bam-mode`).
///
/// Streams strictly shorter than 5 bytes (including empty input) can never
/// be gzip, so niffler cannot sniff them and returns
/// `niffler::Error::FileTooShort`; this is treated the same as
/// [`crate::fastq`]'s `open_reader` equivalent fallback -- the bytes are
/// read as plain (uncompressed) input rather than surfaced as a
/// compression-detection error, so short/empty input reaches
/// [`detect_format`] and gets the clean "input is empty" message instead of
/// an opaque niffler error.
///
/// `niffler::Error::FileTooShort` does not hand back the (fewer than 5)
/// bytes it already consumed while sniffing, so recovering them requires
/// re-reading from the start: for [`InputSpec::Path`] that means re-opening
/// the file (byte-exact, same as `fastq::open_reader`); for
/// [`InputSpec::Stdin`] a pipe's already-consumed bytes cannot be
/// re-obtained at all, so a short *non-empty* stdin input degrades to being
/// treated as empty (the same outcome as a genuinely empty stdin). This is
/// an inherent limitation of sniff-then-consume decoding against an
/// unbuffered, unseekable pipe, not something introduced here; closing it
/// would require pre-buffering all of stdin before any format detection,
/// which isn't worth it (real FASTQ/FASTA/BAM/CRAM input is never 1-4
/// bytes).
///
/// # Errors
///
/// Returns an error if the input cannot be opened, compression detection
/// fails, or the format is unrecognized.
pub fn open_input(spec: &InputSpec) -> Result<(OpenedInput, DetectedFormat)> {
    let (label, decoded): (String, Box<dyn std::io::Read>) = match spec {
        InputSpec::Path(p) => {
            let label = p.display().to_string();
            let file =
                std::fs::File::open(p).with_context(|| format!("opening {}", p.display()))?;
            match niffler::get_reader(Box::new(file)) {
                Ok((decoded, _compression)) => (label, decoded),
                Err(niffler::Error::FileTooShort) => {
                    let file = std::fs::File::open(p)
                        .with_context(|| format!("opening {}", p.display()))?;
                    (label, Box::new(file))
                }
                Err(e) => {
                    return Err(e).with_context(|| format!("detecting compression of {label}"));
                }
            }
        }
        InputSpec::Stdin => {
            let label = "<stdin>".to_string();
            match niffler::get_reader(Box::new(std::io::stdin())) {
                Ok((decoded, _compression)) => (label, decoded),
                Err(niffler::Error::FileTooShort) => (label, Box::new(std::io::empty())),
                Err(e) => {
                    return Err(e).with_context(|| format!("detecting compression of {label}"));
                }
            }
        }
    };
    // Peek the magic prefix by consuming up to 4 bytes (robust to a short
    // first read), classify, then chain those bytes back in front of the rest
    // so the FastqReader sees the byte-exact original stream.
    let mut decoded = decoded;
    let prefix = read_magic_prefix(&mut decoded)?;
    let fmt = detect_format(&prefix)?;
    let opened = match fmt {
        DetectedFormat::Fastq | DetectedFormat::Fasta => {
            let rejoined = std::io::BufReader::new(std::io::Cursor::new(prefix).chain(decoded));
            OpenedInput::Fastq(FastqReader::from_bufread(label, Box::new(rejoined)))
        }
        DetectedFormat::Sam => OpenedInput::Sam,
        DetectedFormat::Bam => OpenedInput::Bam,
        DetectedFormat::Cram => OpenedInput::Cram,
    };
    Ok((opened, fmt))
}

/// Opens `spec` as FASTQ/FASTA, erroring if it is BAM/CRAM (which needs
/// `--bam-mode`). Preserved for the pure-FASTQ callers/tests from Stage 1.
///
/// # Errors
///
/// As [`open_input`], plus an error if the detected format is BAM/CRAM.
pub fn open_fastq_reader(spec: &InputSpec) -> Result<(FastqReader, DetectedFormat)> {
    let label = match spec {
        InputSpec::Path(p) => p.display().to_string(),
        InputSpec::Stdin => "<stdin>".to_string(),
    };
    let (opened, fmt) = open_input(spec)?;
    match opened {
        OpenedInput::Fastq(reader) => Ok((reader, fmt)),
        OpenedInput::Sam | OpenedInput::Bam | OpenedInput::Cram => bail!(
            "{label} is a {fmt:?} file; SAM/BAM/CRAM input requires --bam-mode, not the FASTQ path"
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn detect(bytes: &[u8]) -> DetectedFormat {
        detect_format(bytes).unwrap()
    }

    #[test]
    fn detects_fastq_by_at_sign() {
        assert_eq!(detect(b"@r1\nACGT\n+\nIIII\n"), DetectedFormat::Fastq);
    }

    #[test]
    fn detects_sam_header_but_not_fastq_at_ids() {
        // SAM header lines: `@` + two uppercase letters + TAB.
        assert_eq!(detect(b"@HD\tVN:1.6\n"), DetectedFormat::Sam);
        assert_eq!(detect(b"@SQ\tSN:chr1\n"), DetectedFormat::Sam);
        assert_eq!(detect(b"@RG\tID:x\n"), DetectedFormat::Sam);
        // FASTQ ids that superficially start with `@` + letters are NOT SAM:
        assert_eq!(detect(b"@r1\n"), DetectedFormat::Fastq);
        assert_eq!(detect(b"@HD\n"), DetectedFormat::Fastq); // no tab after @HD
        assert_eq!(detect(b"@Hd\tx"), DetectedFormat::Fastq); // lowercase -> not a SAM type code
    }

    #[test]
    fn detects_fasta_by_gt_sign() {
        assert_eq!(detect(b">seq1\nACGT\n"), DetectedFormat::Fasta);
    }

    #[test]
    fn detects_bam_by_magic() {
        assert_eq!(detect(b"BAM\x01rest"), DetectedFormat::Bam);
    }

    #[test]
    fn detects_cram_by_magic() {
        assert_eq!(detect(b"CRAM\x03rest"), DetectedFormat::Cram);
    }

    #[test]
    fn read_magic_prefix_assembles_magic_split_across_one_byte_reads() {
        // A reader that hands back a single byte per `read` call must still
        // yield the full 4-byte magic: BAM/CRAM streamed one byte at a time
        // must not be misclassified (regression for the single-read bug where
        // a short first read fell through to the unrecognized-format path).
        struct OneByte(Cursor<Vec<u8>>);
        impl std::io::Read for OneByte {
            fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
                let mut one = [0u8; 1];
                let n = std::io::Read::read(&mut self.0, &mut one)?;
                if n == 1 {
                    buf[0] = one[0];
                }
                Ok(n)
            }
        }
        let mut r = OneByte(Cursor::new(b"BAM\x01rest".to_vec()));
        let prefix = read_magic_prefix(&mut r).unwrap();
        assert_eq!(prefix, b"BAM\x01", "the full 4-byte magic must be assembled");
        assert_eq!(detect_format(&prefix).unwrap(), DetectedFormat::Bam);
    }

    #[test]
    fn read_magic_prefix_stops_at_eof_below_prefix_len() {
        // A stream shorter than the magic prefix returns just its bytes (no
        // error, no over-read), so detection still sees the leading byte.
        let mut r = Cursor::new(b">s".to_vec());
        let prefix = read_magic_prefix(&mut r).unwrap();
        assert_eq!(prefix, b">s");
        assert_eq!(detect_format(&prefix).unwrap(), DetectedFormat::Fasta);
    }

    fn write_temp(contents: &[u8]) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        std::io::Write::write_all(&mut f, contents).unwrap();
        std::io::Write::flush(&mut f).unwrap();
        f
    }

    /// `open_fastq_reader`'s `Ok` variant contains `FastqReader`, which
    /// intentionally does not implement `Debug` (it wraps a `Box<dyn
    /// BufRead>`), so `Result::unwrap_err` can't be used directly. This
    /// extracts just the error message for assertions.
    fn open_fastq_reader_err_message(spec: &InputSpec) -> String {
        match open_fastq_reader(spec) {
            Ok(_) => panic!("expected open_fastq_reader to fail for {spec:?}"),
            Err(e) => e.to_string(),
        }
    }

    #[test]
    fn open_fastq_reader_opens_a_plain_fastq_path() {
        let f = write_temp(b"@r1\nACGT\n+\nIIII\n");
        let (mut reader, fmt) =
            open_fastq_reader(&InputSpec::Path(f.path().to_path_buf())).unwrap();
        assert_eq!(fmt, DetectedFormat::Fastq);
        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.id, "r1");
    }

    #[test]
    fn open_fastq_reader_opens_a_fasta_path() {
        let f = write_temp(b">seq1\nACGT\n");
        let (mut reader, fmt) =
            open_fastq_reader(&InputSpec::Path(f.path().to_path_buf())).unwrap();
        assert_eq!(fmt, DetectedFormat::Fasta);
        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.id, "seq1");
    }

    #[test]
    fn open_fastq_reader_errors_cleanly_on_an_empty_path() {
        // niffler::get_reader errors with FileTooShort on inputs under 5
        // bytes (including empty), which open_fastq_reader falls back on
        // (mirroring fastq::open_reader's identical fallback) so this
        // reaches detect_format's clean "input is empty" error rather than
        // an opaque niffler error.
        let f = write_temp(b"");
        let err = open_fastq_reader_err_message(&InputSpec::Path(f.path().to_path_buf()));
        assert!(err.contains("input is empty"), "expected an 'input is empty' error, got: {err}");
    }

    #[test]
    fn open_fastq_reader_detects_a_short_non_empty_fasta_path() {
        // Regression test for the FileTooShort fallback: ">s\nA" is 4 bytes,
        // strictly under niffler's 5-byte sniffing floor (a complete FASTQ
        // record cannot fit in <5 bytes, so a FASTA record is used here). The
        // fallback must recover the FULL original content (by re-opening the
        // file), not silently treat it as empty -- proving it re-reads from
        // the start rather than discarding the short input.
        let payload = b">s\nA";
        assert!(payload.len() < 5, "payload must be under niffler's 5-byte sniffing floor");
        let f = write_temp(payload);
        let (mut reader, fmt) =
            open_fastq_reader(&InputSpec::Path(f.path().to_path_buf())).unwrap();
        assert_eq!(fmt, DetectedFormat::Fasta);
        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.id, "s");
        assert_eq!(rec.seq, b"A");
    }

    #[test]
    fn open_fastq_reader_errors_on_bam_magic_pointing_at_bam_mode() {
        let f = write_temp(b"BAM\x01rest-of-a-fake-bam-payload");
        let err = open_fastq_reader_err_message(&InputSpec::Path(f.path().to_path_buf()));
        assert!(
            err.contains("--bam-mode"),
            "expected the error to point at --bam-mode, got: {err}"
        );
    }

    #[test]
    fn open_fastq_reader_errors_on_cram_magic_pointing_at_bam_mode() {
        let f = write_temp(b"CRAM\x03rest-of-a-fake-cram-payload");
        let err = open_fastq_reader_err_message(&InputSpec::Path(f.path().to_path_buf()));
        assert!(
            err.contains("--bam-mode"),
            "expected the error to point at --bam-mode, got: {err}"
        );
    }

    #[test]
    fn open_fastq_reader_errors_on_a_missing_path() {
        let missing = std::path::PathBuf::from("/nonexistent/path/does-not-exist.fastq");
        assert!(open_fastq_reader(&InputSpec::Path(missing)).is_err());
    }

    #[test]
    fn open_input_returns_fastq_reader_for_fastq() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("r.fq");
        std::fs::write(&p, b"@r1\nACGT\n+\nIIII\n").unwrap();
        let (opened, fmt) = open_input(&InputSpec::Path(p)).unwrap();
        assert_eq!(fmt, DetectedFormat::Fastq);
        let OpenedInput::Fastq(mut reader) = opened else { panic!("expected Fastq") };
        assert_eq!(reader.next_record().unwrap().unwrap().id, "r1");
    }

    #[test]
    fn open_input_marks_bam_and_cram_without_erroring() {
        let dir = tempfile::tempdir().unwrap();

        let bam = dir.path().join("x.bam");
        std::fs::write(&bam, b"BAM\x01rest-of-header-bytes").unwrap();
        let (opened, fmt) = open_input(&InputSpec::Path(bam)).unwrap();
        assert_eq!(fmt, DetectedFormat::Bam);
        assert!(matches!(opened, OpenedInput::Bam));

        let cram = dir.path().join("x.cram");
        std::fs::write(&cram, b"CRAM\x03rest").unwrap();
        let (opened, fmt) = open_input(&InputSpec::Path(cram)).unwrap();
        assert_eq!(fmt, DetectedFormat::Cram);
        assert!(matches!(opened, OpenedInput::Cram));
    }

    #[test]
    fn open_fastq_reader_still_errors_on_bam_pointing_at_bam_mode() {
        // The existing contract must be preserved (Stage 1 tests already cover
        // the message; this re-asserts it after the open_input refactor).
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("x.bam");
        std::fs::write(&bam, b"BAM\x01rest").unwrap();
        let err = open_fastq_reader_err_message(&InputSpec::Path(bam));
        assert!(err.contains("--bam-mode"), "message should still point at --bam-mode: {err}");
    }
}
