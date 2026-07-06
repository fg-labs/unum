//! Tests for `unum_core::kmer_index`, converted from the retired
//! the retired T1K-oracle FFI differential (`diff_kmerindex.rs`) FFI differential (see `tests/common/mod.rs`).
//!
//! The `insert`/`search`/`remove` behavioral tests already asserted explicit
//! Rust-side expected values in the original differential (the C++ oracle was
//! only a cross-check on top of those explicit assertions), so those are kept
//! verbatim with their C++ arm dropped. The `build_index_from_read` cases,
//! which had NO standalone expected values (only Rust-vs-C++), have their Rust
//! output frozen into the `kmer_index.txt` golden.
//!
//! Dropped: `kmerindex_new_produces_a_usable_handle` (an FFI-handle smoke test).

mod common;

use common::Golden;
use unum_core::kmer::KmerCode;
use unum_core::kmer_index::{IndexInfo, KmerIndex};

fn encode(k: usize, bytes: &[u8]) -> KmerCode {
    let mut kc = KmerCode::new(k);
    for &c in bytes {
        kc.append(c);
    }
    kc
}

#[test]
fn insert_and_search_agree_on_basic_entries() {
    let k = 15;
    let mut index = KmerIndex::new();
    let kc = encode(k, b"ACGTACGTACGTACG");
    index.insert(&kc, 1, 0, 1);
    index.insert(&kc, 2, 5, 1);
    assert_eq!(
        index.search(&kc),
        &[IndexInfo { idx: 1, offset: 0 }, IndexInfo { idx: 2, offset: 5 }]
    );
}

#[test]
fn forward_keying_not_canonical_reverse_complement_pair_diverges() {
    // "AAAAAAAAAAAAAAC" and its RC "GTTTTTTTTTTTTTT" share a canonical code
    // but differ in forward code; KmerIndex keys on the forward code, so the
    // two must NOT collide.
    let k = 15;
    let fwd = encode(k, b"AAAAAAAAAAAAAAC");
    let rc = encode(k, b"GTTTTTTTTTTTTTT");
    assert_eq!(fwd.get_canonical_kmer_code(), rc.get_canonical_kmer_code());
    assert_ne!(fwd.get_code(), rc.get_code());

    let mut index = KmerIndex::new();
    index.insert(&fwd, 100, 0, 1);
    index.insert(&rc, 200, 0, 1);

    assert_eq!(index.search(&fwd), &[IndexInfo { idx: 100, offset: 0 }]);
    assert_eq!(index.search(&rc), &[IndexInfo { idx: 200, offset: 0 }]);
}

#[test]
fn invalid_kmer_is_dropped_by_insert_and_search() {
    let k = 15;
    let kc = encode(k, b"ACGTACGTNACGTAC");
    assert!(!kc.is_valid());

    let mut index = KmerIndex::new();
    index.insert(&kc, 42, 7, 1);
    assert!(index.search(&kc).is_empty());
    index.remove(&kc, 42, 7, 1); // silent no-op, must not panic
    assert!(index.search(&kc).is_empty());
}

#[test]
fn repeated_kmer_preserves_insertion_order() {
    let k = 15;
    let kc = encode(k, b"GATTACAGATTACAG");
    let mut index = KmerIndex::new();
    let entries: [(u32, u32, i32); 5] = [(1, 0, 1), (2, 10, 1), (1, 0, -1), (3, 99, 0), (2, 10, 5)];
    for &(idx, offset, strand) in &entries {
        index.insert(&kc, idx, offset, strand);
    }
    let expected: Vec<IndexInfo> =
        entries.iter().map(|&(idx, offset, _)| IndexInfo { idx, offset }).collect();
    assert_eq!(index.search(&kc), expected.as_slice());
}

#[test]
fn remove_matches_t1k_semantics() {
    let k = 15;
    let kc = encode(k, b"TGCATGCATGCATGC");
    let mut index = KmerIndex::new();
    for &(idx, offset) in &[(1u32, 0u32), (2, 0), (1, 0), (3, 7)] {
        index.insert(&kc, idx, offset, 1);
    }
    assert_eq!(
        index.search(&kc),
        &[
            IndexInfo { idx: 1, offset: 0 },
            IndexInfo { idx: 2, offset: 0 },
            IndexInfo { idx: 1, offset: 0 },
            IndexInfo { idx: 3, offset: 7 },
        ]
    );

    // Remove the FIRST matching (1, 0), preserving the rest's relative order.
    index.remove(&kc, 1, 0, 1);
    assert_eq!(
        index.search(&kc),
        &[
            IndexInfo { idx: 2, offset: 0 },
            IndexInfo { idx: 1, offset: 0 },
            IndexInfo { idx: 3, offset: 7 },
        ]
    );

    index.remove(&kc, 99, 99, 1); // no-op
    assert_eq!(
        index.search(&kc),
        &[
            IndexInfo { idx: 2, offset: 0 },
            IndexInfo { idx: 1, offset: 0 },
            IndexInfo { idx: 3, offset: 7 },
        ]
    );

    index.remove(&kc, 1, 0, 1);
    assert_eq!(
        index.search(&kc),
        &[IndexInfo { idx: 2, offset: 0 }, IndexInfo { idx: 3, offset: 7 }]
    );
}

// ---- build_index_from_read (frozen golden) --------------------------------

fn build_index_read() -> &'static [u8] {
    b"GATTACAGATTACAGGGCTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTNNNNNACGTACGTACGTACGTGATTACAGATTACC"
}

fn build_index_read_homopolymer_leading() -> &'static [u8] {
    b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTACGTACGTGATTACAGATTACCGGGCT"
}

fn serialize_hits(hits: &[IndexInfo]) -> String {
    if hits.is_empty() {
        return "-".to_string();
    }
    hits.iter().map(|h| format!("{}:{}", h.idx, h.offset)).collect::<Vec<_>>().join(",")
}

/// Runs `build_index_from_read` and records the search result of every
/// length-`k` window (in order) under a `label`-prefixed key.
fn record_build_index(golden: &mut Golden, read: &[u8], k: usize, id: i32, shift: i32, label: &str) {
    let mut index = KmerIndex::new();
    let mut kc = KmerCode::new(k);
    index.build_index_from_read(&mut kc, read, id, shift);
    for start in 0..=(read.len() - k) {
        let window = &read[start..start + k];
        let q = encode(k, window);
        golden.record(format!("{label}/w{start:03}"), serialize_hits(index.search(&q)));
    }
}

#[test]
fn build_index_from_read_matches_golden() {
    let mut golden = Golden::open("kmer_index.txt");
    record_build_index(&mut golden, build_index_read(), 15, 7, 3, "read/k15/id7/shift3");
    record_build_index(&mut golden, build_index_read(), 31, 12, 100, "read/k31/id12/shift100");
    record_build_index(&mut golden, build_index_read(), 15, 0, 0, "read/k15/id0/shift0");
    record_build_index(&mut golden, build_index_read(), 15, -5, -20, "read/k15/id-5/shift-20");
    record_build_index(&mut golden, build_index_read(), 31, -1, -1, "read/k31/id-1/shift-1");
    record_build_index(
        &mut golden,
        build_index_read_homopolymer_leading(),
        15,
        3,
        0,
        "homopoly/k15/id3/shift0",
    );
    record_build_index(
        &mut golden,
        build_index_read_homopolymer_leading(),
        31,
        -7,
        12,
        "homopoly/k31/id-7/shift12",
    );
    golden.finish();
}
