//! Standalone dump helper for `diff_alignments.rs`'s FFI differential test.
//!
//! # Why a separate PROCESS, not just an in-process FFI call
//!
//! Historically, T1K's vendored `Alignments` linked the bundled legacy
//! `vendor/t1k/samtools-0.1.19` while `rust-htslib`'s vendored modern htslib
//! (`hts-sys`, pulled in transitively by `fg-t1k-core`'s
//! `alignments::Alignments`) was a DIFFERENT, ABI-incompatible library
//! generation -- both exported identically named C symbols (`bam_read1`,
//! `sam_read1`, every `bgzf_*`/`fai_*` function, etc.) with incompatible
//! struct layouts (`bam1_t` was 56 bytes in samtools-0.1.19 vs. a larger,
//! different layout in modern htslib). Linking both into one process let the
//! platform linker silently resolve each duplicate symbol to exactly one
//! winner, corrupting the heap the moment the "wrong" `Alignments` bound to
//! the other generation's `bam_read1`.
//!
//! As of the htslib unification (see `crates/fg-t1k-sys/build.rs`), the
//! vendored `Alignments` class here is compiled `-DHTSLIB` against the SAME
//! htslib build (via hts-sys) that `rust-htslib` itself links, so that
//! original hard ABI blocker no longer exists -- an in-process
//! `CppAlignments` FFI call would be safe today. This binary is kept as a
//! separate process anyway purely to avoid touching a large, already-correct
//! test harness (`diff_alignments.rs`) as a side effect of the htslib
//! unification; it depends on `fg_t1k_sys` (the shim + `Alignments`,
//! compiled against hts-sys's htslib, never linking `rust-htslib`
//! directly), reads a BAM path from `argv[1]`, and prints every record's
//! `Alignments` fields plus the `GetGeneralInfo` summary to stdout in a
//! simple, line-oriented, tab-separated format `diff_alignments.rs` parses.
//!
//! # Output format
//!
//! ```text
//! RECORD<TAB>read_seq<TAB>qual<TAB>read_id<TAB>is_first_mate<TAB>is_reverse<TAB>is_mate_reverse<TAB>is_aligned<TAB>is_template_aligned<TAB>is_primary<TAB>chrom_id<TAB>seg_count<TAB>a0,b0;a1,b1;...
//! ...(one RECORD line per Alignments::Next() call, in file order)...
//! GENERAL_INFO<TAB>frag_stdev<TAB>read_len
//! ```
//!
//! `GENERAL_INFO` is computed via a FRESH `Alignments` instance (re-opened on
//! the same path), matching how `diff_alignments.rs`'s general_info tests
//! open a dedicated reader rather than reusing one already advanced by
//! `Next()` -- avoids coupling the RECORD dump's stream position to the
//! GetGeneralInfo scan's own internal consumption of the file.

use std::process::ExitCode;

#[cfg(not(feature = "t1k-sys"))]
fn main() -> ExitCode {
    eprintln!(
        "alignments_oracle_dump: built without the `t1k-sys` feature (no-op stub); \
         rebuild with `--features t1k-sys` to use this binary."
    );
    ExitCode::FAILURE
}

#[cfg(feature = "t1k-sys")]
fn main() -> ExitCode {
    run()
}

#[cfg(feature = "t1k-sys")]
fn run() -> ExitCode {
    use fg_t1k_sys::CppAlignments;
    use std::env;
    use std::path::PathBuf;

    let args: Vec<String> = env::args().collect();
    let Some(path) = args.get(1) else {
        eprintln!("usage: alignments_oracle_dump <path.bam>");
        return ExitCode::FAILURE;
    };
    let path = PathBuf::from(path);

    let mut alignments = CppAlignments::new();
    alignments.open(&path);

    while alignments.next() {
        let seq = String::from_utf8_lossy(&alignments.get_read_seq()).into_owned();
        let qual: String = alignments
            .get_qual()
            .iter()
            .map(|&q| {
                // Encode each qual byte as a 3-digit zero-padded decimal
                // joined by commas -- qual bytes are already phred+33 ASCII
                // per Alignments::GetQual, but some may not be printable/
                // TSV-safe (e.g. embedded whitespace at very low quality),
                // so a numeric encoding is unambiguous and trivially parsed.
                format!("{q}")
            })
            .collect::<Vec<_>>()
            .join(",");
        let read_id = alignments.get_read_id();
        let segments = alignments
            .segments()
            .iter()
            .map(|s| format!("{},{}", s.a, s.b))
            .collect::<Vec<_>>()
            .join(";");

        println!(
            "RECORD\t{seq}\t{qual}\t{read_id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{segments}",
            b2i(alignments.is_first_mate()),
            b2i(alignments.is_reverse()),
            b2i(alignments.is_mate_reverse()),
            b2i(alignments.is_aligned()),
            b2i(alignments.is_template_aligned()),
            b2i(alignments.is_primary()),
            alignments.get_chrom_id(),
            alignments.segments().len(),
        );
    }

    // Fresh instance for GetGeneralInfo, matching the doc comment above.
    let mut gi_alignments = CppAlignments::new();
    gi_alignments.open(&path);
    let info = gi_alignments.general_info(false);
    println!("GENERAL_INFO\t{}\t{}", info.frag_stdev, info.read_len);

    ExitCode::SUCCESS
}

#[cfg(feature = "t1k-sys")]
fn b2i(b: bool) -> i32 {
    i32::from(b)
}
