//! Standalone dump helper for `diff_alignments.rs`'s FFI differential test.
//!
//! # Why a separate PROCESS, not just an in-process FFI call
//!
//! `vendor/t1k/samtools-0.1.19` (the bundled, legacy htslib T1K's `Alignments`
//! links against here) and `rust-htslib`'s vendored modern htslib
//! (`hts-sys`, pulled in transitively by `fg-t1k-core`'s
//! `alignments::Alignments`) export a large overlapping set of IDENTICALLY
//! NAMED C symbols with INCOMPATIBLE ABIs -- `bam_read1`, `sam_read1`,
//! `bam_write1`, every `bgzf_*`/`fai_*` function, etc. (confirmed by
//! disassembling a test binary that linked both: `bam1_t`'s struct layout
//! differs between the two library generations, e.g. samtools-0.1.19's
//! `bam1_t` is 56 bytes vs. modern htslib's larger layout). When both
//! libraries are linked into ONE process, the platform linker resolves each
//! duplicate symbol to exactly ONE winner (macOS `ld`'s flat-namespace
//! archive-member resolution silently picks one, with no diagnostic) -- and
//! T1K's own `Alignments::Next`/`GetGeneralInfo` (compiled WITHOUT
//! `-DHTSLIB`, so it calls the LEGACY `samread`->`bam_read1` chain) ends up
//! silently bound to modern htslib's INCOMPATIBLE `bam_read1`, corrupting the
//! heap on the very first record read (verified via `lldb`: the crash is a
//! `bam1_t`-sized `calloc`/`free` mismatch inside `bam_destroy1`, caused by
//! `bam_read1` writing through the wrong struct layout).
//!
//! Since this ABI collision is a genuine, unavoidable consequence of linking
//! both library generations into one address space (not a bug in either
//! port), the fix is to run the C++ oracle in a SEPARATE PROCESS from the
//! Rust reader (which needs `rust-htslib`/modern htslib to read the same
//! fixture BAM). This binary is that separate process: it depends ONLY on
//! `fg_t1k_sys` (which pulls in samtools-0.1.19 + the shim, never
//! rust-htslib), reads a BAM path from `argv[1]`, and prints every record's
//! `Alignments` fields plus the `GetGeneralInfo` summary to stdout in a
//! simple, line-oriented, tab-separated format `diff_alignments.rs` parses.
//! `diff_alignments.rs` itself never links `fg_t1k_sys`'s FFI/shim code
//! directly (only spawns this binary via `std::process::Command`), so its
//! own process only ever links ONE htslib generation (rust-htslib's), same
//! as this binary only ever links the OTHER (samtools-0.1.19's) -- neither
//! process has the collision.
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
