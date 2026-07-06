//! Golden-file test for `unum_core::ref_kmer_filter::get_overlaps_from_read`,
//! converted from the retired T1K-oracle FFI differential (`diff_get_overlaps_from_read.rs`) FFI
//! differential (see `tests/common/mod.rs`). Drives the same adversarial read
//! battery (exact substrings, mismatch sweeps, small indels, multi-gap /
//! noncolinear concatenations, RC reads, boundary lengths, heavy mutation,
//! pure noise, large seeded random batch) and freezes each overlap's every
//! field into the `get_overlaps.txt` golden.
//!
//! # Freeze sensitivity: `similarity` is a bit-exact `f64`
//!
//! `similarity` is a plain, deterministic `matchCnt / (seqSpan + readSpan)`
//! ratio -- the golden freezes it as `f64::to_bits()`, locking the port's
//! exact float output (it was byte-identical to the C++ `SeqSet` oracle). This
//! is the highest-freeze-sensitivity field in the suite; a change to any of
//! the overlap arithmetic will flip these bits and fail the golden.

mod common;

use common::Golden;
use unum_core::ref_kmer_filter::{RefKmerFilter, Scratch};
use std::path::{Path, PathBuf};

const KMER_LENGTH: usize = 9;

fn fixture_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/refbuild/golden").join(rel)
}

struct FilterPair {
    rust: RefKmerFilter,
    records: Vec<(String, Vec<u8>)>,
}

fn load_pair(fasta_rel: &str) -> FilterPair {
    let path = fixture_path(fasta_rel);
    let rust = RefKmerFilter::from_reference_fasta(&path, KMER_LENGTH)
        .unwrap_or_else(|e| panic!("from_reference_fasta({fasta_rel}): {e}"));
    let records = parse_fasta_for_test(&path);
    assert_eq!(records.len(), rust.seq_count(), "record-count disagreement for {fasta_rel}");
    FilterPair { rust, records }
}

fn parse_fasta_for_test(path: &Path) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                records.push((id, std::mem::take(&mut current_seq)));
            }
            current_id = Some(rest.split_whitespace().next().unwrap_or("").to_string());
        } else {
            current_seq.extend_from_slice(line.trim_end().as_bytes());
        }
    }
    if let Some(id) = current_id {
        records.push((id, current_seq));
    }
    records
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            other => panic!("reverse_complement: unsupported base {other}"),
        })
        .collect()
}

struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        Self { state: seed | 1 }
    }
    fn next_u64(&mut self) -> u64 {
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }
    fn next_base(&mut self) -> u8 {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        BASES[(self.next_u64() % 4) as usize]
    }
    fn next_base_excluding(&mut self, exclude: u8) -> u8 {
        loop {
            let b = self.next_base();
            if b != exclude {
                return b;
            }
        }
    }
    fn next_range(&mut self, bound: usize) -> usize {
        usize::try_from(self.next_u64() % bound as u64).unwrap_or(usize::MAX)
    }
}

fn pseudo_random_acgt(seed: u64, len: usize) -> Vec<u8> {
    let mut rng = Xorshift64::new(seed);
    (0..len).map(|_| rng.next_base()).collect()
}

fn mutate_k_positions(seq: &[u8], k: usize, seed: u64) -> Vec<u8> {
    let mut rng = Xorshift64::new(seed);
    let mut out = seq.to_vec();
    let mut positions: Vec<usize> = (0..seq.len()).collect();
    for i in 0..k.min(positions.len()) {
        let j = i + rng.next_range(positions.len() - i);
        positions.swap(i, j);
    }
    for &pos in positions.iter().take(k.min(positions.len())) {
        out[pos] = rng.next_base_excluding(out[pos]);
    }
    out
}

fn delete_span(seq: &[u8], pos: usize, del_len: usize) -> Vec<u8> {
    let mut out = seq.to_vec();
    let end = (pos + del_len).min(out.len());
    out.drain(pos..end);
    out
}

fn insert_span(seq: &[u8], pos: usize, ins_len: usize, seed: u64) -> Vec<u8> {
    let mut out = seq[..pos].to_vec();
    out.extend(pseudo_random_acgt(seed, ins_len));
    out.extend_from_slice(&seq[pos..]);
    out
}

/// Serializes the `get_overlaps_from_read` result: `None`, or
/// `count|o0|o1|...` where each `o` is
/// `seq_idx:strand:read_start:read_end:seq_start:seq_end:match_cnt:sim_bits`
/// (similarity as exact `f64::to_bits()`).
fn record_overlaps(golden: &mut Golden, label: &str, pair: &FilterPair, read: &[u8]) {
    let mut scratch = Scratch::default();
    let value = match pair.rust.get_overlaps_from_read(read, &mut scratch) {
        None => "None".to_string(),
        Some(overlaps) => {
            let parts = overlaps
                .iter()
                .map(|o| {
                    format!(
                        "{}:{}:{}:{}:{}:{}:{}:{}",
                        o.seq_idx,
                        o.strand,
                        o.read_start,
                        o.read_end,
                        o.seq_start,
                        o.seq_end,
                        o.match_cnt,
                        o.similarity.to_bits()
                    )
                })
                .collect::<Vec<_>>();
            format!("{}|{}", overlaps.len(), parts.join("|"))
        }
    };
    golden.record(label, value);
}

#[test]
fn get_overlaps_matches_golden() {
    let mut golden = Golden::open("get_overlaps.txt");
    let pair = load_pair("kir_rna_seq.fa");
    let hla_pair = load_pair("hla_rna_seq.fa");

    // kir exact substrings + RC.
    for &(rec_idx, start, len) in
        &[(0, 0, 60), (0, 200, 100), (1, 500, 150), (2, 50, 80), (3, 900, 120), (5, 10, 200), (6, 300, 45)]
    {
        let (id, seq) = &pair.records[rec_idx];
        let read = &seq[start..start + len];
        record_overlaps(&mut golden, &format!("kir_substr_{id}@{start}+{len}"), &pair, read);
        record_overlaps(&mut golden, &format!("kir_substr_RC_{id}@{start}+{len}"), &pair, &reverse_complement(read));
    }

    // hla exact substrings.
    for &(rec_idx, start, len) in &[(0, 0, 70), (1, 300, 90), (4, 700, 110), (8, 150, 60), (11, 700, 130)] {
        let (id, seq) = &hla_pair.records[rec_idx];
        let read = &seq[start..start + len];
        record_overlaps(&mut golden, &format!("hla_substr_{id}@{start}+{len}"), &hla_pair, read);
    }

    // mismatch sweeps.
    {
        let (id, seq) = &pair.records[2];
        let base = &seq[100..250];
        for k in 0..=40usize {
            record_overlaps(&mut golden, &format!("kir_mm_{id}_k{k}"), &pair, &mutate_k_positions(base, k, 1000 + k as u64));
        }
    }
    {
        let (id, seq) = &hla_pair.records[4];
        let base = &seq[200..300];
        for k in 0..=30usize {
            record_overlaps(&mut golden, &format!("hla_mm_{id}_k{k}"), &hla_pair, &mutate_k_positions(base, k, 2000 + k as u64));
        }
    }

    // small deletions / insertions.
    {
        let (id, seq) = &pair.records[3];
        let base = &seq[200..400];
        for (del_pos, del_len) in [(50, 1), (50, 2), (50, 3), (100, 5), (150, 8), (30, 1), (170, 4)] {
            record_overlaps(&mut golden, &format!("kir_del_{id}_p{del_pos}_l{del_len}"), &pair, &delete_span(base, del_pos, del_len));
        }
        for (ins_pos, ins_len, seed) in [(50, 1, 11), (50, 2, 12), (50, 3, 13), (100, 5, 14), (150, 8, 15), (30, 1, 16)] {
            record_overlaps(&mut golden, &format!("kir_ins_{id}_p{ins_pos}_l{ins_len}"), &pair, &insert_span(base, ins_pos, ins_len, seed));
        }
    }
    {
        let (id, seq) = &hla_pair.records[6];
        let base = &seq[150..330];
        for (op, pos, l, seed) in [("del", 40, 2usize, 0u64), ("del", 90, 4, 0), ("ins", 40, 2, 21), ("ins", 90, 4, 22), ("del", 20, 1, 0), ("ins", 160, 3, 23)] {
            let read = if op == "del" { delete_span(base, pos, l) } else { insert_span(base, pos, l, seed) };
            record_overlaps(&mut golden, &format!("hla_{op}_{id}_p{pos}_l{l}"), &hla_pair, &read);
        }
    }
    // combined mismatches + indels.
    {
        let (id, seq) = &pair.records[0];
        let base = &seq[300..520];
        record_overlaps(&mut golden, &format!("kir_del_mutate_{id}"), &pair, &mutate_k_positions(&delete_span(base, 60, 3), 5, 555));
        record_overlaps(&mut golden, &format!("kir_ins_mutate_{id}"), &pair, &mutate_k_positions(&insert_span(base, 120, 4, 777), 4, 888));
    }

    // widely-separated / multi-gap.
    {
        let (id0, seq0) = &pair.records[0];
        let (id1, seq1) = &pair.records[1];
        let frag_a = &seq0[50..125];
        let mut same = frag_a.to_vec();
        same.extend_from_slice(&seq0[900..975]);
        record_overlaps(&mut golden, &format!("{id0}_noncolinear_same"), &pair, &same);
        let mut cross = frag_a.to_vec();
        cross.extend_from_slice(&seq1[400..475]);
        record_overlaps(&mut golden, &format!("{id0}_{id1}_noncolinear_cross"), &pair, &cross);
        let mut near = seq0[200..280].to_vec();
        near.extend_from_slice(&seq0[285..365]);
        record_overlaps(&mut golden, &format!("{id0}_near_adjacent"), &pair, &near);
        let mut three = seq0[10..90].to_vec();
        three.extend_from_slice(&seq0[500..580]);
        three.extend_from_slice(&seq0[950..1030.min(seq0.len())]);
        record_overlaps(&mut golden, &format!("{id0}_three_frags"), &pair, &three);
    }

    // reads from each record.
    for (idx, (id, seq)) in pair.records.iter().enumerate() {
        if seq.len() < 100 {
            continue;
        }
        let start = seq.len() / 4;
        let len = 70.min(seq.len() - start);
        record_overlaps(&mut golden, &format!("record{idx}_{id}@{start}+{len}"), &pair, &seq[start..start + len]);
    }

    // RC mutated / indel.
    {
        let (id, seq) = &pair.records[4];
        let base = &seq[80..200];
        record_overlaps(&mut golden, &format!("{id}_rc_mutated"), &pair, &reverse_complement(&mutate_k_positions(base, 8, 555)));
        record_overlaps(&mut golden, &format!("{id}_rc_del"), &pair, &reverse_complement(&delete_span(base, 40, 3)));
        record_overlaps(&mut golden, &format!("{id}_rc_ins"), &pair, &reverse_complement(&insert_span(base, 40, 3, 999)));
        let mut concat = seq[10..90].to_vec();
        concat.extend_from_slice(&seq[600..680]);
        record_overlaps(&mut golden, &format!("{id}_rc_concat"), &pair, &reverse_complement(&concat));
    }

    // boundary lengths.
    {
        let (id, seq) = &pair.records[3];
        for len in [35usize, 40, 44, 45, 46, 50, 60] {
            record_overlaps(&mut golden, &format!("{id}_boundary_len{len}"), &pair, &seq[300..300 + len]);
        }
        record_overlaps(&mut golden, &format!("{id}_len_k"), &pair, &seq[50..50 + KMER_LENGTH]);
        record_overlaps(&mut golden, &format!("{id}_len_k+1"), &pair, &seq[50..=50 + KMER_LENGTH]);
        record_overlaps(&mut golden, &format!("{id}_len5"), &pair, &seq[0..5]);
    }

    // heavily mutated (may filter to empty).
    {
        let (id, seq) = &hla_pair.records[2];
        let base = &seq[300..420];
        for k in [20usize, 30, 40, 50, 60] {
            record_overlaps(&mut golden, &format!("{id}_heavy_k{k}"), &hla_pair, &mutate_k_positions(base, k.min(120), 3000 + k as u64));
        }
    }

    // pure noise.
    for (seed, len) in [(1u64, 80), (2, 120), (3, 60), (4, 150), (5, 40)] {
        record_overlaps(&mut golden, &format!("noise_seed{seed}_len{len}"), &pair, &pseudo_random_acgt(seed, len));
    }

    // large seeded random batch (hla).
    {
        let seq_count = hla_pair.records.len();
        for seed in 0u64..150 {
            let mut rng = Xorshift64::new(seed ^ 0x5EED_C0DE);
            let len = 40 + rng.next_range(161);
            let category = rng.next_range(5);
            let read: Vec<u8> = match category {
                0 => pseudo_random_acgt(seed.wrapping_mul(31) + 7, len),
                1 => {
                    let rec_idx = rng.next_range(seq_count);
                    let (_, seq) = &hla_pair.records[rec_idx];
                    let len = len.min(seq.len());
                    let max_start = seq.len() - len;
                    let start = if max_start == 0 { 0 } else { rng.next_range(max_start + 1) };
                    seq[start..start + len].to_vec()
                }
                2 => {
                    let rec_idx = rng.next_range(seq_count);
                    let (_, seq) = &hla_pair.records[rec_idx];
                    let len = len.min(seq.len());
                    let max_start = seq.len() - len;
                    let start = if max_start == 0 { 0 } else { rng.next_range(max_start + 1) };
                    let base = &seq[start..start + len];
                    let k = rng.next_range(len / 3 + 1);
                    mutate_k_positions(base, k, seed.wrapping_mul(97) + 13)
                }
                3 => {
                    let rec_idx = rng.next_range(seq_count);
                    let (_, seq) = &hla_pair.records[rec_idx];
                    let len = len.min(seq.len()).max(seq.len().min(20));
                    let max_start = seq.len() - len;
                    let start = if max_start == 0 { 0 } else { rng.next_range(max_start + 1) };
                    let base = &seq[start..start + len];
                    let del_pos = rng.next_range(base.len().saturating_sub(10).max(1));
                    let del_len = 1 + rng.next_range(5);
                    delete_span(base, del_pos, del_len)
                }
                _ => {
                    let rec_idx = rng.next_range(seq_count);
                    let (_, seq) = &hla_pair.records[rec_idx];
                    let len = len.min(seq.len()).max(seq.len().min(20));
                    let max_start = seq.len() - len;
                    let start = if max_start == 0 { 0 } else { rng.next_range(max_start + 1) };
                    let base = &seq[start..start + len];
                    let ins_pos = rng.next_range(base.len().saturating_sub(10).max(1));
                    let ins_len = 1 + rng.next_range(5);
                    insert_span(base, ins_pos, ins_len, seed.wrapping_mul(211) + 3)
                }
            };
            record_overlaps(&mut golden, &format!("random_batch_seed{seed}_cat{category}"), &hla_pair, &read);
        }
    }

    golden.finish();
}
