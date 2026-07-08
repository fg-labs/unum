//! Binary-level smoke test for the `unum extract` subcommand: invokes the
//! actual `unum` binary on the pinned KIR example and asserts the emitted
//! candidate FASTQs are byte-identical to the frozen extraction goldens shared
//! with the library-level `unum-core/tests/golden_fastq_extract.rs`. Guards the
//! real `stages::extract::run` CLI path (the CLI defaults `-s` /
//! initial-kmer-length match the parameters the golden was captured with).
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

/// A path under the workspace-level `fixtures/example/` directory.
fn example(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/example").join(rel)
}

/// A path under `unum-core`'s frozen `fastq_extract` golden directory (the same
/// goldens the library-level extraction test asserts against).
fn golden(rel: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../unum-core/tests/golden/fastq_extract")
        .join(rel)
}

#[test]
fn extract_example_matches_golden() {
    let tmp = tempfile::tempdir().unwrap();
    let prefix = tmp.path().join("cand");
    let status = Command::new(env!("CARGO_BIN_EXE_unum"))
        .args([
            "extract",
            "-f",
            example("ref/kir_rna_seq.fa").to_str().unwrap(),
            "-1",
            example("example_1.fq").to_str().unwrap(),
            "-2",
            example("example_2.fq").to_str().unwrap(),
            "-o",
            prefix.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success(), "`unum extract` exited non-zero");

    for (out, gold) in [("cand_1.fq", "example_paired_1.fq"), ("cand_2.fq", "example_paired_2.fq")]
    {
        let got = std::fs::read(tmp.path().join(out)).unwrap();
        let want = std::fs::read(golden(gold)).unwrap();
        assert_eq!(got, want, "`unum extract` {out} must be byte-identical to golden {gold}");
    }
}

// -- `-i/--input` interleaved-vs-paired-split equivalence -------------------
//
// Programmatic fixtures (no committed data files): a small repetitive-but-
// balanced reference plus exact-substring paired reads, long/unique enough
// to clear the default candidate filter (`hitLenRequired`=27 paired,
// `-s`=0.8) and actually produce non-empty candidate output, proving the
// `-i` interleaved path runs the real extraction pipeline (not just two
// empty files trivially matching).

/// A 165bp repetitive-but-balanced reference sequence, long enough that
/// 60bp exact substrings comfortably clear the default paired
/// `hitLenRequired` (27) and `-s` (0.8) filter thresholds.
const REF_SEQ: &str = "ACGTACGTGGATTACAGATTACAGATTACAGATTACAGCCCTGACGTGTGACGTGTGACGTGTGACGTGTGGATCAGATCAGATCAGATCAGGATCCATGGATCCATGGATCCATGACTGACTGACTGACTGCATGCATGCATGCATGGTACGTACGTACGTACG";

/// Writes a single-contig reference FASTA derived from [`REF_SEQ`].
fn write_reference(tmp: &tempfile::TempDir) -> PathBuf {
    let path = tmp.path().join("ref.fa");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, ">contig1\n{REF_SEQ}").unwrap();
    path
}

/// Writes a small paired FASTQ fixture: 10 pairs, each mate a 60bp exact
/// substring of [`REF_SEQ`] (so both mates stay on-reference and every pair
/// is emitted as a candidate). The same id (`r{i}`) is used for both mates,
/// matching `fixtures/example`'s convention and satisfying the mate-name
/// guard / interleaved auto-detection (which compares suffix-stripped ids).
///
/// Mate 1 and mate 2 use DISTINCT offsets (mate 2 mirrors mate 1's offset:
/// `off2 = max_off - off1`), so `_1.fq` and `_2.fq` carry different
/// sequences per pair. This is deliberate: it makes a within-pair mate
/// role-swap (emitting rec2 as mate 1) observable in the byte-comparison of
/// [`i_flag_interleaved_matches_paired_split`]. If both mates shared the same
/// sequence, a swap would be invisible and the test — the only oracle-free
/// proof of interleaved correctness — could pass on broken output.
fn write_paired_reads(tmp: &tempfile::TempDir) -> (PathBuf, PathBuf) {
    let r1_path = tmp.path().join("r1.fq");
    let r2_path = tmp.path().join("r2.fq");
    let mut f1 = std::fs::File::create(&r1_path).unwrap();
    let mut f2 = std::fs::File::create(&r2_path).unwrap();
    let read_len = 60;
    let max_off = REF_SEQ.len() - read_len;
    for i in 0..10 {
        let off1 = (i * 10) % max_off;
        // Mirror offset: distinct from `off1` for every i in 0..10, and still
        // within [0, max_off] so mate 2 stays a valid on-reference substring.
        let off2 = max_off - off1;
        let seq1 = &REF_SEQ[off1..off1 + read_len];
        let seq2 = &REF_SEQ[off2..off2 + read_len];
        let qual = "I".repeat(read_len);
        writeln!(f1, "@r{i}\n{seq1}\n+\n{qual}").unwrap();
        writeln!(f2, "@r{i}\n{seq2}\n+\n{qual}").unwrap();
    }
    (r1_path, r2_path)
}

/// Splits a FASTQ file's contents into its 4-line records (as joined
/// strings), preserving order.
fn fastq_records(path: &Path) -> Vec<String> {
    let contents = std::fs::read_to_string(path).unwrap();
    let lines: Vec<&str> = contents.lines().collect();
    lines.chunks(4).map(|rec| rec.join("\n")).collect()
}

/// NEW helper (not reused from elsewhere): interleaves two mate FASTQ files
/// record-by-record (`R1[0],R2[0],R1[1],R2[1],...`) into a single output
/// file, for exercising `-i`'s single-input interleaved auto-detection.
fn write_interleaved_from(tmp: &tempfile::TempDir, r1: &Path, r2: &Path) -> PathBuf {
    let recs1 = fastq_records(r1);
    let recs2 = fastq_records(r2);
    assert_eq!(recs1.len(), recs2.len(), "mate files must have the same record count");

    let out_path = tmp.path().join("interleaved.fq");
    let mut f = std::fs::File::create(&out_path).unwrap();
    for (rec1, rec2) in recs1.iter().zip(recs2.iter()) {
        writeln!(f, "{rec1}").unwrap();
        writeln!(f, "{rec2}").unwrap();
    }
    out_path
}

/// Runs `unum extract` with `args`, asserting a clean exit.
fn run_extract(args: &[&str]) {
    let status =
        Command::new(env!("CARGO_BIN_EXE_unum")).arg("extract").args(args).status().unwrap();
    assert!(status.success(), "`unum extract` exited non-zero for args {args:?}");
}

fn path_str(p: &Path) -> &str {
    p.to_str().unwrap()
}

/// The `{prefix}_1.fq` output of a paired/interleaved `extract` run.
fn paired_out_1(prefix: &Path) -> PathBuf {
    PathBuf::from(format!("{}_1.fq", prefix.to_str().unwrap()))
}

/// The `{prefix}_2.fq` output of a paired/interleaved `extract` run.
fn paired_out_2(prefix: &Path) -> PathBuf {
    PathBuf::from(format!("{}_2.fq", prefix.to_str().unwrap()))
}

#[test]
fn i_flag_interleaved_matches_paired_split() {
    let tmp = tempfile::tempdir().unwrap();
    let reference = write_reference(&tmp);
    let (r1, r2) = write_paired_reads(&tmp);
    let interleaved = write_interleaved_from(&tmp, &r1, &r2);

    // Paired split via legacy flags.
    let out_paired = tmp.path().join("paired");
    run_extract(&[
        "-f",
        path_str(&reference),
        "-1",
        path_str(&r1),
        "-2",
        path_str(&r2),
        "-o",
        path_str(&out_paired),
    ]);

    // Interleaved via -i (single input, auto-detected).
    let out_i = tmp.path().join("viai");
    run_extract(&[
        "-f",
        path_str(&reference),
        "-i",
        path_str(&interleaved),
        "-o",
        path_str(&out_i),
    ]);

    // The two candidate FASTQs must be byte-identical -- interleaved has no
    // separate oracle, so functional equivalence to the paired split is the
    // only correctness proof available.
    assert_eq!(
        std::fs::read(paired_out_1(&out_paired)).unwrap(),
        std::fs::read(paired_out_1(&out_i)).unwrap(),
        "-i interleaved mate-1 output must match the -1/-2 paired split"
    );
    assert_eq!(
        std::fs::read(paired_out_2(&out_paired)).unwrap(),
        std::fs::read(paired_out_2(&out_i)).unwrap(),
        "-i interleaved mate-2 output must match the -1/-2 paired split"
    );
}
