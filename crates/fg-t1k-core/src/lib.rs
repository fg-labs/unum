#![forbid(unsafe_code)]
//! fg-t1k core: the pure-Rust port of T1K's pipeline. No C++ dependency.

pub mod alignments;
pub mod extract;
pub mod fastq;
pub mod kmer;
pub mod kmer_count;
pub mod kmer_index;
pub(crate) mod overlap;
pub mod ref_kmer_filter;
pub mod refbuild;
