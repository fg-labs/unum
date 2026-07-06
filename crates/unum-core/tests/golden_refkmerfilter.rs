//! Golden-file test for `unum_core::ref_kmer_filter`'s
//! `is_low_complexity`/`is_good_candidate`/`passes_bucket_count_gate_only`,
//! converted from the retired T1K-oracle FFI differential (`diff_refkmerfilter.rs`) FFI
//! differential (see `tests/common/mod.rs`). Drives the same large
//! adversarial read battery (exact substrings, mismatch sweeps, indels,
//! noncolinear concatenations, RC reads, low-complexity, seeded random) and
//! freezes the three per-read booleans into the `refkmerfilter.txt` golden
//! (they were byte-identical to the C++ `SeqSet` oracle).
//!
//! Dropped from the original: the `#[should_panic]` oversized-read guards and
//! the max-length-read acceptance test (all guarded the C++ shim's fixed
//! `rcBuf[100001]` buffer, which no longer exists).

mod common;

use common::Golden;
use unum_core::ref_kmer_filter::{RefKmerFilter, Scratch, is_low_complexity};
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

/// Records the three per-read predicates under `label`:
/// `low_complexity,good_candidate,gate1_only`.
fn record_predicates(golden: &mut Golden, label: &str, pair: &FilterPair, read: &[u8]) {
    let low = is_low_complexity(read);
    let good = pair.rust.is_good_candidate(read);
    let mut scratch = Scratch::default();
    let gate1 = pair.rust.passes_bucket_count_gate_only(read, &mut scratch);
    golden.record(label, format!("{},{},{}", u8::from(low), u8::from(good), u8::from(gate1)));
}

#[test]
fn refkmerfilter_matches_golden() {
    let mut golden = Golden::open("refkmerfilter.txt");
    let pair = load_pair("kir_rna_seq.fa");
    let hla_pair = load_pair("hla_rna_seq.fa");

    // kir exact substrings + RC.
    assert!(pair.rust.seq_count() >= 7);
    for &(rec_idx, start, len) in
        &[(0, 0, 60), (0, 200, 100), (1, 500, 150), (2, 50, 80), (3, 900, 120), (5, 10, 200), (6, 300, 45)]
    {
        let (id, seq) = &pair.records[rec_idx];
        let read = &seq[start..start + len];
        record_predicates(&mut golden, &format!("kir_substr_{id}@{start}+{len}"), &pair, read);
        let rc = reverse_complement(read);
        record_predicates(&mut golden, &format!("kir_substr_RC_{id}@{start}+{len}"), &pair, &rc);
    }

    // hla exact substrings.
    assert!(hla_pair.rust.seq_count() >= 12);
    for &(rec_idx, start, len) in &[(0, 0, 70), (1, 300, 90), (4, 700, 110), (8, 150, 60), (11, 700, 130)] {
        let (id, seq) = &hla_pair.records[rec_idx];
        let read = &seq[start..start + len];
        record_predicates(&mut golden, &format!("hla_substr_{id}@{start}+{len}"), &hla_pair, read);
    }

    // substring with a single N.
    {
        let (id, seq) = &pair.records[0];
        let (start, len) = (100, 100);
        let mut read = seq[start..start + len].to_vec();
        read[len / 2] = b'N';
        record_predicates(&mut golden, &format!("kir_substr_N_{id}@{start}+{len}"), &pair, &read);
    }

    // random non-matching reads.
    for (seed, len) in [(1u64, 80), (2, 120), (3, 60), (4, 150), (5, 40)] {
        let read = pseudo_random_acgt(seed, len);
        record_predicates(&mut golden, &format!("random_seed{seed}_len{len}"), &pair, &read);
    }

    // low-complexity reads.
    record_predicates(&mut golden, "homopolymer_A", &pair, &vec![b'A'; 100]);
    record_predicates(&mut golden, "homopolymer_T", &pair, &vec![b'T'; 60]);
    record_predicates(
        &mut golden,
        "dinucleotide_AT",
        &pair,
        &(0..100).map(|i| if i % 2 == 0 { b'A' } else { b'T' }).collect::<Vec<_>>(),
    );
    record_predicates(
        &mut golden,
        "dinucleotide_GC",
        &pair,
        &(0..80).map(|i| if i % 2 == 0 { b'G' } else { b'C' }).collect::<Vec<_>>(),
    );
    {
        let mut mostly_n = pseudo_random_acgt(42, 100);
        for b in mostly_n.iter_mut().take(20) {
            *b = b'N';
        }
        record_predicates(&mut golden, "mostly_N", &pair, &mostly_n);
    }

    // short reads.
    {
        let (id, seq) = &pair.records[0];
        record_predicates(&mut golden, &format!("{id}_len5"), &pair, &seq[0..5]);
        record_predicates(&mut golden, &format!("{id}_len_k"), &pair, &seq[50..50 + KMER_LENGTH]);
        record_predicates(&mut golden, &format!("{id}_len_k+1"), &pair, &seq[50..=50 + KMER_LENGTH]);
    }

    // mismatch sweep kir.
    {
        let (id, seq) = &pair.records[2];
        let (start, len) = (100, 150);
        let base = &seq[start..start + len];
        for k in 0..=40usize {
            let read = mutate_k_positions(base, k, 1000 + k as u64);
            record_predicates(&mut golden, &format!("mm_sweep_kir_{id}_k{k}"), &pair, &read);
        }
    }
    // mismatch sweep hla.
    {
        let (id, seq) = &hla_pair.records[4];
        let (start, len) = (200, 100);
        let base = &seq[start..start + len];
        for k in 0..=30usize {
            let read = mutate_k_positions(base, k, 2000 + k as u64);
            record_predicates(&mut golden, &format!("mm_sweep_hla_{id}_k{k}"), &hla_pair, &read);
        }
    }

    // noncolinear concatenations.
    {
        let (id0, seq0) = &pair.records[0];
        let (id1, seq1) = &pair.records[1];
        let frag_a = &seq0[50..125];
        let mut same = frag_a.to_vec();
        same.extend_from_slice(&seq0[900..975]);
        record_predicates(&mut golden, &format!("{id0}_noncolinear_same"), &pair, &same);
        let mut cross = frag_a.to_vec();
        cross.extend_from_slice(&seq1[400..475]);
        record_predicates(&mut golden, &format!("{id0}_{id1}_noncolinear_cross"), &pair, &cross);
        let mut near = seq0[200..280].to_vec();
        near.extend_from_slice(&seq0[285..365]);
        record_predicates(&mut golden, &format!("{id0}_near_adjacent"), &pair, &near);
    }

    // high-frequency-kmer skip against a synthetic repetitive reference.
    {
        use std::io::Write as _;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, ">only").unwrap();
        writeln!(tmp, "{}", "ACGT".repeat(150)).unwrap();
        tmp.flush().unwrap();
        let kmer_length = 4usize;
        let rust = RefKmerFilter::from_reference_fasta(tmp.path(), kmer_length).unwrap();
        let reference = "ACGT".repeat(150);
        let hit_read = &reference.as_bytes()[10..80];
        assert!(!is_low_complexity(hit_read));
        golden.record("highfreq_hit", u8::from(rust.is_good_candidate(hit_read)).to_string());
        let unrelated = "TGCA".repeat(20).into_bytes();
        golden.record("highfreq_unrelated", u8::from(rust.is_good_candidate(&unrelated)).to_string());
    }

    // boundary length reads.
    {
        let (id, seq) = &pair.records[3];
        for len in [35usize, 40, 44, 45, 46, 50, 60] {
            record_predicates(&mut golden, &format!("{id}_boundary_len{len}"), &pair, &seq[300..300 + len]);
        }
    }

    // tandem repeat / tie-inducing reads.
    {
        let (id, seq) = &pair.records[5];
        let motif = &seq[40..70];
        let mut repeated = Vec::new();
        for _ in 0..4 {
            repeated.extend_from_slice(motif);
        }
        record_predicates(&mut golden, &format!("{id}_tandem_x4"), &pair, &repeated);
        let mut padded = pseudo_random_acgt(9001, 20);
        for _ in 0..3 {
            padded.extend_from_slice(motif);
        }
        padded.extend_from_slice(&pseudo_random_acgt(9002, 20));
        record_predicates(&mut golden, &format!("{id}_padded_tandem_x3"), &pair, &padded);
        let motif2 = &seq[500..580];
        let mut two_copies = motif2.to_vec();
        two_copies.extend_from_slice(&mutate_k_positions(motif2, 3, 424_242));
        record_predicates(&mut golden, &format!("{id}_two_motif_copies"), &pair, &two_copies);
    }

    // rc / N-containing / short reads.
    {
        let (id, seq) = &pair.records[4];
        let base = &seq[80..200];
        let rc_mutated = reverse_complement(&mutate_k_positions(base, 8, 555));
        record_predicates(&mut golden, &format!("{id}_rc_mutated"), &pair, &rc_mutated);
        let mut concat = seq[10..90].to_vec();
        concat.extend_from_slice(&seq[600..680]);
        record_predicates(&mut golden, &format!("{id}_rc_concat"), &pair, &reverse_complement(&concat));
        let mut with_one_n = seq[300..400].to_vec();
        with_one_n[50] = b'N';
        record_predicates(&mut golden, &format!("{id}_single_N"), &pair, &with_one_n);
        let mut with_few_n = seq[300..400].to_vec();
        for &pos in &[10, 40, 70] {
            with_few_n[pos] = b'N';
        }
        record_predicates(&mut golden, &format!("{id}_three_N"), &pair, &with_few_n);
        let short_exact = &seq[500..500 + KMER_LENGTH + 5];
        record_predicates(&mut golden, &format!("{id}_short_exact"), &pair, short_exact);
        record_predicates(&mut golden, &format!("{id}_short_mutated"), &pair, &mutate_k_positions(short_exact, 2, 777));
    }

    // large random batch (hla).
    {
        let seq_count = hla_pair.records.len();
        for seed in 0u64..220 {
            let mut rng = Xorshift64::new(seed ^ 0xABCD_1234);
            let len = 40 + rng.next_range(161);
            let category = rng.next_range(3);
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
                _ => {
                    let rec_idx = rng.next_range(seq_count);
                    let (_, seq) = &hla_pair.records[rec_idx];
                    let len = len.min(seq.len());
                    let max_start = seq.len() - len;
                    let start = if max_start == 0 { 0 } else { rng.next_range(max_start + 1) };
                    let base = &seq[start..start + len];
                    let k = rng.next_range(len / 3 + 1);
                    mutate_k_positions(base, k, seed.wrapping_mul(97) + 13)
                }
            };
            record_predicates(&mut golden, &format!("random_batch_seed{seed}_cat{category}"), &hla_pair, &read);
        }
    }

    golden.finish();
}
