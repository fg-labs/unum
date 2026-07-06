//! Golden-file test for `fg_t1k_core::kmer_count`, converted from the retired
//! `fg-t1k-sys` `diff_kmercount.rs` FFI differential (see `tests/common/mod.rs`).
//! Drives the same `add_count` sequences and `get_count`/`jaccard_similarity`
//! queries the differential drove, freezing the Rust output (which was
//! byte-identical to the C++ `KmerCount` oracle) into `kmer_count.txt`.
//!
//! Dropped from the original: `kmercount_new_produces_a_usable_handle` (an
//! FFI-handle smoke test) and the `#[should_panic]` `get_count_rejects_query_
//! shorter_than_k` FFI-boundary guard.

mod common;

use common::Golden;
use fg_t1k_core::kmer_count::KmerCount;

fn read_set_a() -> Vec<&'static [u8]> {
    vec![
        b"ACGTACGTACGTACGTACGTACGTACGTACGT",
        b"TTTTAAAACCCCGGGGTTTTAAAACCCCGGGG",
        b"GATTACAGATTACAGATTACAGATTACAGATT",
        b"ACGTACGTACGTACGTACGTACGTACGTACGT",
        b"ACGTNACGTACGTNACGTACGTACGTACGTAC",
        b"AC",
        b"CCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA",
    ]
}

fn read_set_b() -> Vec<&'static [u8]> {
    vec![
        b"CCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTT",
        b"AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
        b"ACGTACGTACGTACGTACGTACGTACGTACGT",
        b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        b"TG",
    ]
}

fn all_windows(reads: &[&[u8]], k: usize) -> Vec<Vec<u8>> {
    let mut windows = Vec::new();
    for read in reads {
        if read.len() < k {
            continue;
        }
        for start in 0..=(read.len() - k) {
            windows.push(read[start..start + k].to_vec());
        }
    }
    windows.sort();
    windows.dedup();
    windows
}

fn absent_and_n_queries(k: usize) -> Vec<Vec<u8>> {
    let mut queries = Vec::new();
    queries.push(vec![b'A'; k]);
    queries.push(vec![b'T'; k]);
    let mut alt = Vec::with_capacity(k);
    for i in 0..k {
        alt.push(if i % 2 == 0 { b'G' } else { b'C' });
    }
    queries.push(alt);
    let mut n_at_start = vec![b'A'; k];
    n_at_start[0] = b'N';
    queries.push(n_at_start);
    let mut n_at_end = vec![b'A'; k];
    n_at_end[k - 1] = b'N';
    queries.push(n_at_end);
    if k > 2 {
        let mut n_in_middle = vec![b'A'; k];
        n_in_middle[k / 2] = b'N';
        queries.push(n_in_middle);
    }
    queries
}

/// Builds a `KmerCount` for `reads`, recording every `add_count` return value
/// and every present-window / absent-query `get_count` result under
/// `label`-prefixed keys. Returns the populated counter for jaccard use.
fn add_reads_and_record(golden: &mut Golden, k: usize, reads: &[&[u8]], label: &str) -> KmerCount {
    let mut kc = KmerCount::new(k);
    for (i, read) in reads.iter().enumerate() {
        let ret = kc.add_count(read);
        golden.record(format!("{label}/k{k}/add/{i:03}"), ret.to_string());
    }
    for (i, window) in all_windows(reads, k).into_iter().enumerate() {
        let c = kc.get_count(&window);
        golden.record(format!("{label}/k{k}/present/{i:04}"), c.to_string());
    }
    for (i, query) in absent_and_n_queries(k).into_iter().enumerate() {
        let c = kc.get_count(&query);
        golden.record(format!("{label}/k{k}/absent/{i:02}"), c.to_string());
    }
    kc
}

#[test]
fn kmer_count_matches_golden() {
    let mut golden = Golden::open("kmer_count.txt");

    for k in [15usize, 31] {
        add_reads_and_record(&mut golden, k, &read_set_a(), "read_set_a");
    }

    for k in [15usize, 31] {
        let a = add_reads_and_record(&mut golden, k, &read_set_a(), "jaccard/read_set_a");
        let b = add_reads_and_record(&mut golden, k, &read_set_b(), "jaccard/read_set_b");
        let jaccard = a.jaccard_similarity(&b);
        golden.record(format!("jaccard/k{k}"), jaccard.to_bits().to_string());
    }

    golden.finish();
}
