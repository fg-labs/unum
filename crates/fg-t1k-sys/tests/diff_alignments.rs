#![cfg(feature = "t1k-sys")]
//! Differential test: builds a small, coordinate-sorted BAM programmatically
//! (via `rust_htslib::bam::Writer`), reads it with the Rust
//! `fg_t1k_core::alignments::Alignments` port IN-PROCESS, and separately runs
//! the REAL C++ `vendor/t1k/alignments.hpp` `Alignments` class over the SAME
//! BAM via a spawned subprocess (`alignments_oracle_dump`, see below),
//! parsing its output and asserting per-record + global agreement.
//!
//! # Why a subprocess, not an in-process FFI call via `CppAlignments`
//!
//! `vendor/t1k/samtools-0.1.19` (the bundled, legacy htslib the vendored
//! `Alignments` class links against) and `rust-htslib`'s vendored MODERN
//! htslib (`hts-sys`, pulled in transitively by `fg-t1k-core`'s
//! `alignments::Alignments`, which this test also needs to build/read the
//! fixture) export a large overlapping set of IDENTICALLY NAMED C symbols
//! with INCOMPATIBLE ABIs -- `bam_read1`, `sam_read1`, `bam_write1`, every
//! `bgzf_*`/`fai_*` function, etc. Linking BOTH into one process causes the
//! platform linker to silently resolve each duplicate symbol to exactly ONE
//! winner (confirmed via `lldb` + `otool -tv` disassembly: the final test
//! binary's `bam_read1`/`sam_read1` bind to MODERN htslib's implementation
//! everywhere, including inside the vendored `Alignments::Next`, which is
//! compiled WITHOUT `-DHTSLIB` and therefore expects the LEGACY
//! samtools-0.1.19 `bam1_t` layout -- e.g. samtools-0.1.19's `bam1_t` is 56
//! bytes vs. modern htslib's differently-laid-out, larger struct). The
//! result is a `calloc`/`free` size mismatch and an immediate SIGSEGV inside
//! `bam_destroy1` on the very first record read -- a genuine, unavoidable
//! ABI collision from linking two htslib generations into one address space,
//! not a bug in either port.
//!
//! The fix: run the C++ oracle in a SEPARATE PROCESS
//! (`alignments_oracle_dump`, `src/bin/alignments_oracle_dump.rs`) that
//! depends ONLY on `fg_t1k_sys` (samtools-0.1.19 + the shim, never
//! rust-htslib). This test file spawns that binary via
//! `std::process::Command` and parses its line-oriented stdout output --
//! it never calls `fg_t1k_sys::CppAlignments`/the FFI module directly, so
//! THIS process links only rust-htslib (via `fg_t1k_core`), and the spawned
//! process links only samtools-0.1.19 -- neither process has both, so
//! neither has the collision. This is still a byte-for-byte differential
//! against the REAL, unmodified vendored C++ `Alignments` class on the exact
//! same BAM bytes; only the process boundary differs from a more typical
//! in-process opaque-handle FFI comparison.
//!
//! # Coverage
//!
//! The test BAM (`build_test_bam`) covers every category the task brief
//! requires:
//! - A forward-strand, first-mate, primary, mapped record.
//! - A reverse-strand record (exercises `GetReadSeq`/`GetQual`'s
//!   reverse-complement/reversal path -- the single highest-risk trap this
//!   port must get right).
//! - A spliced (N-CIGAR) record, producing multiple `segments`.
//! - A first-mate/second-mate pair (same QNAME, `0x40`/`0x80` flags).
//! - An unmapped record whose mate IS mapped (paired, `0x4` set, `0x8`
//!   unset) -- exercises `IsAligned() == false` while `IsTemplateAligned()`
//!   can still be `true` (its mate is aligned).
//! - An unmapped-TEMPLATE pair: both mates unmapped (`0x4` AND `0x8` both
//!   set on each record) -- exercises `IsTemplateAligned() == false` via the
//!   `(flag & 0xd) == 0xd` branch.
//!
//! Every record is checked for: `get_read_seq`, `get_qual`, `get_read_id`,
//! `is_first_mate`, `is_reverse`, `is_mate_reverse`, `is_aligned`,
//! `is_template_aligned`, `is_primary`, `get_chrom_id`, segment count, and
//! each segment's `(a, b)`. The global `frag_stdev`/`read_len` from
//! `GetGeneralInfo` are checked separately (`general_info_agrees_with_cpp`).

use fg_t1k_core::alignments::Alignments;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{self, Header, Writer};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Builds a small coordinate-sorted BAM at `path` with two contigs
/// (`chr19`, `chr20`) and the record categories documented in the module
/// docs above.
// Long by construction: each of the 8 records is a distinct, deliberately
// documented coverage category (see module docs) -- splitting this into
// several tiny functions would scatter closely related per-record flag/
// CIGAR bookkeeping across the file for no readability benefit.
#[allow(clippy::too_many_lines)]
fn build_test_bam(path: &Path) {
    let mut header = Header::new();

    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);

    let mut sq_chr19 = HeaderRecord::new(b"SQ");
    sq_chr19.push_tag(b"SN", "chr19");
    sq_chr19.push_tag(b"LN", 58_617_616);
    header.push_record(&sq_chr19);

    let mut sq_chr20 = HeaderRecord::new(b"SQ");
    sq_chr20.push_tag(b"SN", "chr20");
    sq_chr20.push_tag(b"LN", 64_444_167);
    header.push_record(&sq_chr20);

    let mut writer = Writer::from_path(path, &header, bam::Format::Bam).expect("create writer");

    // --- Record 1: forward-strand, first-mate, primary, mapped, paired.
    let mut r1 = bam::Record::new();
    let r1_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let r1_len = u32::try_from(r1_seq.len()).unwrap();
    r1.set(
        b"pair_fwd",
        Some(&CigarString(vec![Cigar::Match(r1_len)])),
        r1_seq,
        &vec![30u8; r1_seq.len()],
    );
    r1.set_tid(0);
    r1.set_pos(100);
    r1.set_mtid(0);
    r1.set_mpos(300);
    // paired(0x1) + proper_pair(0x2) + mate_reverse(0x20) + first_in_template(0x40)
    r1.set_flags(0x1 | 0x2 | 0x20 | 0x40);
    writer.write(&r1).expect("write r1");

    // --- Record 2: second mate of pair_fwd, REVERSE-strand -- exercises
    // GetReadSeq/GetQual's RC path directly on a real paired record.
    let mut r2 = bam::Record::new();
    let asym_seq = b"AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACC"; // 37bp, deliberately asymmetric composition
    let asym_len = u32::try_from(asym_seq.len()).unwrap();
    r2.set(
        b"pair_fwd",
        Some(&CigarString(vec![Cigar::Match(asym_len)])),
        asym_seq,
        &(0..asym_seq.len()).map(|i| 2 + u8::try_from(i % 40).unwrap()).collect::<Vec<u8>>(),
    );
    r2.set_tid(0);
    r2.set_pos(300);
    r2.set_mtid(0);
    r2.set_mpos(100);
    // paired + proper_pair + reverse(0x10) + last_in_template(0x80)
    r2.set_flags(0x1 | 0x2 | 0x10 | 0x80);
    writer.write(&r2).expect("write r2");

    // --- Record 3: spliced (N-CIGAR) single-end-style record on chr19,
    // reverse strand too (covers RC + multi-segment simultaneously).
    let mut r3 = bam::Record::new();
    let spliced_seq = b"TTGGCCAATTGGCCAATTGGCCAA"; // 24bp
    r3.set(
        b"spliced_read",
        Some(&CigarString(vec![Cigar::Match(12), Cigar::RefSkip(500), Cigar::Match(12)])),
        spliced_seq,
        &[25u8; 24],
    );
    r3.set_tid(0);
    r3.set_pos(1000);
    r3.set_mtid(-1);
    r3.set_mpos(-1);
    r3.set_flags(0x10); // reverse only, unpaired
    writer.write(&r3).expect("write r3");

    // --- Record 4: mapped mate with an UNMAPPED partner -- exercises
    // IsAligned()==false while IsTemplateAligned() can still be true.
    let mut r4 = bam::Record::new();
    let r4_seq = b"GATTACAGATTACAGATTACAGATTACAG";
    let r4_len = u32::try_from(r4_seq.len()).unwrap();
    r4.set(
        b"one_unmapped",
        Some(&CigarString(vec![Cigar::Match(r4_len)])),
        r4_seq,
        &vec![20u8; r4_seq.len()],
    );
    r4.set_tid(0);
    r4.set_pos(2000);
    r4.set_mtid(0);
    r4.set_mpos(2000); // mate "coordinate" for sort purposes; mate itself is unmapped
    // paired + first_in_template (mate_unmapped=0x8 set on THIS record since ITS mate is unmapped)
    r4.set_flags(0x1 | 0x8 | 0x40);
    writer.write(&r4).expect("write r4");

    let mut r4_mate = bam::Record::new();
    let r4_mate_seq = b"CTGTAATCTGTAATCTGTAATCTGTAAT";
    r4_mate.set(b"one_unmapped", None, r4_mate_seq, &vec![20u8; r4_mate_seq.len()]);
    r4_mate.set_tid(0);
    r4_mate.set_pos(2000);
    r4_mate.set_mtid(0);
    r4_mate.set_mpos(2000);
    // paired + unmapped(0x4) + last_in_template(0x80); mate (r4) IS mapped so
    // mate_unmapped(0x8) is NOT set here.
    r4_mate.set_flags(0x1 | 0x4 | 0x80);
    writer.write(&r4_mate).expect("write r4_mate");

    // --- Record 5/6: an UNMAPPED-TEMPLATE pair -- both mates unmapped,
    // exercising IsTemplateAligned()'s `(flag & 0xd) == 0xd` branch.
    let mut r5 = bam::Record::new();
    r5.set(b"both_unmapped", None, b"NNNNNNNNNNNNNNNNNNNN", &[2u8; 20]);
    r5.set_tid(-1);
    r5.set_pos(-1);
    r5.set_mtid(-1);
    r5.set_mpos(-1);
    // paired + unmapped(0x4) + mate_unmapped(0x8) + first_in_template(0x40)
    r5.set_flags(0x1 | 0x4 | 0x8 | 0x40);
    writer.write(&r5).expect("write r5");

    let mut r6 = bam::Record::new();
    r6.set(b"both_unmapped", None, b"AAAAAAAAAAAAAAAAAAAA", &[3u8; 20]);
    r6.set_tid(-1);
    r6.set_pos(-1);
    r6.set_mtid(-1);
    r6.set_mpos(-1);
    r6.set_flags(0x1 | 0x4 | 0x8 | 0x80);
    writer.write(&r6).expect("write r6");

    // --- Record 7: a second spliced FORWARD-strand read on chr20, three
    // exons, for extra segment-count coverage on the second contig. Sequence
    // generated programmatically (3x15bp = 45bp) rather than hand-counted,
    // so the CIGAR-implied length always matches the literal's actual
    // length.
    let mut r7 = bam::Record::new();
    let three_exon_seq: Vec<u8> = (0..45).map(|i| b"ACGT"[i % 4]).collect();
    r7.set(
        b"three_exon",
        Some(&CigarString(vec![
            Cigar::Match(15),
            Cigar::RefSkip(200),
            Cigar::Match(15),
            Cigar::RefSkip(300),
            Cigar::Match(15),
        ])),
        &three_exon_seq,
        &vec![35u8; three_exon_seq.len()],
    );
    r7.set_tid(1);
    r7.set_pos(50_000);
    r7.set_mtid(-1);
    r7.set_mpos(-1);
    r7.set_flags(0); // forward, unpaired
    writer.write(&r7).expect("write r7");

    drop(writer);
}

fn fixture_bam() -> PathBuf {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("diff_alignments_test.bam");
    build_test_bam(&path);
    // Leak the tempdir so the path stays valid for the caller (test-only,
    // process-lifetime leak is fine here).
    std::mem::forget(dir);
    path
}

/// A single record's fields as parsed from `alignments_oracle_dump`'s
/// `RECORD` line -- see that binary's module docs for the exact format. The
/// six `bool` fields are six INDEPENDENT `Alignments` flag predicates (not
/// interacting state-machine flags), directly mirroring the six separate
/// C++ methods `assert_record_agrees` checks -- a set of enums would not
/// simplify this.
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Clone, PartialEq, Eq)]
struct OracleRecord {
    read_seq: Vec<u8>,
    qual: Vec<u8>,
    read_id: String,
    is_first_mate: bool,
    is_reverse: bool,
    is_mate_reverse: bool,
    is_aligned: bool,
    is_template_aligned: bool,
    is_primary: bool,
    chrom_id: i32,
    segments: Vec<(i64, i64)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct OracleGeneralInfo {
    frag_stdev: i32,
    read_len: i32,
}

/// Locates the `alignments_oracle_dump` binary built alongside this test
/// binary. Cargo places `[[bin]]`/`src/bin/*.rs` targets in the same
/// `target/<profile>/` directory as the test harness executable itself, one
/// level up from `target/<profile>/deps/` where the test binary lives.
fn oracle_dump_bin_path() -> PathBuf {
    let test_exe = std::env::current_exe().expect("current_exe");
    // test_exe is .../target/<profile>/deps/diff_alignments-<hash>
    let deps_dir = test_exe.parent().expect("deps dir");
    let profile_dir = deps_dir.parent().expect("profile dir");
    let candidate = profile_dir.join("alignments_oracle_dump");
    assert!(
        candidate.exists(),
        "alignments_oracle_dump binary not found at {candidate:?} -- expected it to be built \
         alongside the fg-t1k-sys test binaries (src/bin/alignments_oracle_dump.rs)"
    );
    candidate
}

/// Runs the `alignments_oracle_dump` subprocess against `bam_path` and
/// parses its stdout into `(per-record dump, general info)`.
fn run_oracle_dump(bam_path: &Path) -> (Vec<OracleRecord>, OracleGeneralInfo) {
    let bin = oracle_dump_bin_path();
    let output = Command::new(&bin)
        .arg(bam_path)
        .output()
        .unwrap_or_else(|e| panic!("spawning {bin:?}: {e}"));
    assert!(
        output.status.success(),
        "alignments_oracle_dump exited with {:?}; stderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("oracle dump stdout must be UTF-8");

    let mut records = Vec::new();
    let mut general_info = None;
    for line in stdout.lines() {
        let mut fields = line.split('\t');
        match fields.next() {
            Some("RECORD") => {
                let read_seq = fields.next().expect("read_seq field").as_bytes().to_vec();
                let qual: Vec<u8> = fields
                    .next()
                    .expect("qual field")
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .map(|s| s.parse::<u8>().expect("qual byte"))
                    .collect();
                let read_id = fields.next().expect("read_id field").to_string();
                let is_first_mate = fields.next().expect("is_first_mate") == "1";
                let is_reverse = fields.next().expect("is_reverse") == "1";
                let is_mate_reverse = fields.next().expect("is_mate_reverse") == "1";
                let is_aligned = fields.next().expect("is_aligned") == "1";
                let is_template_aligned = fields.next().expect("is_template_aligned") == "1";
                let is_primary = fields.next().expect("is_primary") == "1";
                let chrom_id: i32 = fields.next().expect("chrom_id").parse().expect("chrom_id int");
                let _seg_count: usize =
                    fields.next().expect("seg_count").parse().expect("seg_count int");
                let segments_field = fields.next().unwrap_or("");
                let segments: Vec<(i64, i64)> = if segments_field.is_empty() {
                    Vec::new()
                } else {
                    segments_field
                        .split(';')
                        .map(|pair| {
                            let mut it = pair.split(',');
                            let a: i64 = it.next().unwrap().parse().unwrap();
                            let b: i64 = it.next().unwrap().parse().unwrap();
                            (a, b)
                        })
                        .collect()
                };
                records.push(OracleRecord {
                    read_seq,
                    qual,
                    read_id,
                    is_first_mate,
                    is_reverse,
                    is_mate_reverse,
                    is_aligned,
                    is_template_aligned,
                    is_primary,
                    chrom_id,
                    segments,
                });
            }
            Some("GENERAL_INFO") => {
                let frag_stdev: i32 = fields.next().expect("frag_stdev").parse().unwrap();
                let read_len: i32 = fields.next().expect("read_len").parse().unwrap();
                general_info = Some(OracleGeneralInfo { frag_stdev, read_len });
            }
            Some(other) => panic!("unexpected oracle dump line prefix: {other:?}"),
            None => {}
        }
    }

    (records, general_info.expect("GENERAL_INFO line missing from oracle dump output"))
}

/// Asserts every per-record accessor [`Alignments`] exposes agrees with the
/// corresponding parsed [`OracleRecord`] for the CURRENT rust record (must
/// already have had `next()` called and returned success for this index).
fn assert_record_agrees(label: &str, rust: &Alignments, oracle: &OracleRecord) {
    assert_eq!(rust.read_seq(), oracle.read_seq, "{label}: read_seq mismatch");
    assert_eq!(rust.qual(), oracle.qual, "{label}: qual mismatch");
    assert_eq!(rust.read_id(), oracle.read_id, "{label}: read_id mismatch");
    assert_eq!(rust.is_first_mate(), oracle.is_first_mate, "{label}: is_first_mate mismatch");
    assert_eq!(rust.is_reverse(), oracle.is_reverse, "{label}: is_reverse mismatch");
    assert_eq!(rust.is_mate_reverse(), oracle.is_mate_reverse, "{label}: is_mate_reverse mismatch");
    assert_eq!(rust.is_aligned(), oracle.is_aligned, "{label}: is_aligned mismatch");
    assert_eq!(
        rust.is_template_aligned(),
        oracle.is_template_aligned,
        "{label}: is_template_aligned mismatch"
    );
    assert_eq!(rust.is_primary(), oracle.is_primary, "{label}: is_primary mismatch");
    assert_eq!(rust.chrom_id(), oracle.chrom_id, "{label}: chrom_id mismatch");

    let rust_segs: Vec<(i64, i64)> = rust.segments().iter().map(|s| (s.a, s.b)).collect();
    assert_eq!(rust_segs, oracle.segments, "{label}: segments mismatch");
}

#[test]
fn per_record_fields_agree_across_all_categories() {
    let path = fixture_bam();

    let mut rust = Alignments::open(&path).expect("open rust Alignments");
    let (oracle_records, _) = run_oracle_dump(&path);

    let labels = [
        "r1 forward first-mate",
        "r2 reverse second-mate",
        "r3 spliced reverse unpaired",
        "r4 mapped-with-unmapped-mate",
        "r4_mate unmapped-with-mapped-mate",
        "r5 both-unmapped first-mate",
        "r6 both-unmapped second-mate",
        "r7 three-exon forward",
    ];
    assert_eq!(
        oracle_records.len(),
        labels.len(),
        "oracle dump record count mismatch (fixture BAM changed?)"
    );

    let mut seen = 0usize;
    for (label, oracle_record) in labels.iter().zip(oracle_records.iter()) {
        let rust_has_next = rust.next().expect("rust next");
        assert!(rust_has_next, "{label}: expected a record but rust reader hit EOF");
        assert_record_agrees(label, &rust, oracle_record);
        seen += 1;
    }

    // Rust side must also reach EOF at exactly the same point (no extra
    // records the oracle dump didn't see).
    let rust_next = rust.next().expect("rust next at end");
    assert!(!rust_next, "rust reader has unexpected extra records");
    assert_eq!(seen, labels.len());
}

/// Directly proves the reverse-strand record (r2) genuinely exercises the
/// reverse-complement path: its Rust-decoded `read_seq` must differ from a
/// naive (non-RC'd) decode of the record's raw stored SEQ, and must equal
/// the true reverse complement of the stored (reference-forward) bases --
/// AND must match the real C++ oracle's output for that same record.
#[test]
fn reverse_strand_record_is_genuinely_reverse_complemented() {
    let path = fixture_bam();
    let mut rust = Alignments::open(&path).expect("open rust Alignments");

    // r1 (forward), then r2 (the reverse-strand record).
    assert!(rust.next().unwrap());
    assert!(rust.next().unwrap());
    assert!(rust.is_reverse(), "expected r2 to be the reverse-strand record");

    let decoded = rust.read_seq();
    let stored_forward = b"AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACC"; // what was written to SEQ
    assert_ne!(
        decoded.as_slice(),
        stored_forward.as_slice(),
        "reverse-strand read_seq must not equal the raw stored (reference-forward) bases"
    );
    assert_ne!(
        decoded,
        stored_forward.iter().rev().copied().collect::<Vec<u8>>(),
        "reverse-strand read_seq must not be a plain reversal without complementing"
    );

    let mut expected = stored_forward.to_vec();
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
    assert_eq!(decoded, expected, "read_seq must be the true reverse complement");

    // Also cross-check against the real C++ oracle directly (not just the
    // hand-derived expectation above), so this test's "genuinely exercises
    // RC" claim is itself validated against stock, not just self-consistent
    // Rust-side logic.
    let (oracle_records, _) = run_oracle_dump(&path);
    let oracle_r2 = &oracle_records[1];
    assert!(oracle_r2.is_reverse);
    assert_eq!(decoded, oracle_r2.read_seq, "rust RC output must match real C++ Alignments");
}

/// `frag_stdev`/`read_len` from `GetGeneralInfo` must agree between the two
/// sides. Uses a SEPARATE, larger synthetic BAM (paired, same-chrom,
/// opposite-strand mates at varied insert sizes) since the small
/// multi-category BAM above is not representative of a real fragment-size
/// distribution.
#[test]
fn general_info_agrees_with_cpp() {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("general_info_test.bam");

    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr19");
    sq.push_tag(b"LN", 58_617_616);
    header.push_record(&sq);

    {
        let mut writer = Writer::from_path(&path, &header, bam::Format::Bam).expect("writer");
        // 29 pairs, insert sizes varying, read length 50, mates on opposite
        // strands (proper orientation), same chromosome, first mate's pos <
        // mate's pos (matches GetGeneralInfo's sampling gate). NOTE: 29 is
        // deliberately NOT a multiple of 10 -- 29*0.7 = 20.3 (non-integer), so
        // the smallest-70% count is ceil(20.3)=21, exercising the loop-count
        // (ceil) semantics of alignments.hpp:669. A 30-pair fixture (30*0.7=21.0
        // exact) would mask a floor-vs-ceil off-by-one in that count.
        let seq: Vec<u8> = (0..50).map(|i| b"ACGT"[i % 4]).collect();
        let qual = vec![30u8; 50];
        for i in 0..29i64 {
            let pos1 = 1000 + i * 1000;
            let insert = 300 + i;
            let pos2 = pos1 + insert - 50;

            let mut m1 = bam::Record::new();
            m1.set(
                format!("pair{i}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(50)])),
                &seq,
                &qual,
            );
            m1.set_tid(0);
            m1.set_pos(pos1);
            m1.set_mtid(0);
            m1.set_mpos(pos2);
            m1.set_flags(0x1 | 0x2 | 0x20 | 0x40); // paired, proper, mate_reverse, first
            writer.write(&m1).expect("write m1");

            let mut m2 = bam::Record::new();
            m2.set(
                format!("pair{i}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(50)])),
                &seq,
                &qual,
            );
            m2.set_tid(0);
            m2.set_pos(pos2);
            m2.set_mtid(0);
            m2.set_mpos(pos1);
            m2.set_flags(0x1 | 0x2 | 0x10 | 0x80); // paired, proper, reverse, last
            writer.write(&m2).expect("write m2");
        }
    }

    let mut rust = Alignments::open(&path).expect("open rust");
    let rust_info = rust.general_info(false).expect("rust general_info");

    let (_, oracle_info) = run_oracle_dump(&path);

    assert_eq!(
        rust_info.read_len, oracle_info.read_len,
        "read_len mismatch: rust={rust_info:?} oracle={oracle_info:?}"
    );
    assert_eq!(
        rust_info.frag_stdev, oracle_info.frag_stdev,
        "frag_stdev mismatch: rust={rust_info:?} oracle={oracle_info:?}"
    );
    // This dataset is deliberately paired with real insert-size variance, so
    // frag_stdev must be > 0 on both sides -- proves the test actually
    // exercises the paired branch of GetGeneralInfo, not just the
    // single-end/frag_stdev==0 shortcut.
    assert!(oracle_info.frag_stdev > 0, "test data must exercise the paired frag_stdev>0 branch");
}

/// A single-end-only BAM (no paired flag on any record) must yield
/// `frag_stdev == 0` on both sides -- the single-end sentinel branch of
/// `GetGeneralInfo`.
#[test]
fn general_info_single_end_yields_zero_frag_stdev_on_both_sides() {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("single_end_test.bam");

    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "coordinate");
    header.push_record(&hd);
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr19");
    sq.push_tag(b"LN", 58_617_616);
    header.push_record(&sq);

    {
        let mut writer = Writer::from_path(&path, &header, bam::Format::Bam).expect("writer");
        let seq: Vec<u8> = (0..75).map(|i| b"ACGT"[i % 4]).collect();
        let qual = vec![30u8; 75];
        for i in 0..20i64 {
            let mut r = bam::Record::new();
            r.set(
                format!("single{i}").as_bytes(),
                Some(&CigarString(vec![Cigar::Match(75)])),
                &seq,
                &qual,
            );
            r.set_tid(0);
            r.set_pos(1000 + i * 200);
            r.set_mtid(-1);
            r.set_mpos(-1);
            r.set_flags(0); // unpaired, forward
            writer.write(&r).expect("write single-end record");
        }
    }

    let mut rust = Alignments::open(&path).expect("open rust");
    let rust_info = rust.general_info(false).expect("rust general_info");

    let (_, oracle_info) = run_oracle_dump(&path);

    assert_eq!(rust_info.read_len, oracle_info.read_len);
    assert_eq!(rust_info.frag_stdev, 0, "rust: single-end must yield frag_stdev==0");
    assert_eq!(oracle_info.frag_stdev, 0, "cpp: single-end must yield frag_stdev==0");
}

/// `rewind()` must reset the Rust reader to the beginning and reproduce the
/// FIRST record's fields again, matching `BamExtractor.cpp`'s own
/// `GetGeneralInfo(true); Rewind();` pattern -- checked against the oracle
/// dump's first record (a fresh oracle process is inherently "rewound").
#[test]
fn rewind_resets_to_first_record() {
    let path = fixture_bam();

    let mut rust = Alignments::open(&path).expect("open rust");
    let (oracle_records, _) = run_oracle_dump(&path);
    let first = &oracle_records[0];

    assert!(rust.next().unwrap());
    assert_eq!(rust.read_id(), first.read_id);

    // Advance a couple more records.
    assert!(rust.next().unwrap());
    assert!(rust.next().unwrap());

    rust.rewind().expect("rust rewind");

    assert!(rust.next().unwrap());
    assert_eq!(rust.read_id(), first.read_id, "rust rewind must return to the first record");
    assert_record_agrees("post-rewind first record", &rust, first);
}
