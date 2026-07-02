#![forbid(unsafe_code)]
//! fg-t1k core: the pure-Rust port of T1K's pipeline. No C++ dependency.

pub mod align_algo;
pub mod alignments;
pub mod bam_extract;
pub mod extract;
pub mod fastq;
pub mod genotyper;
pub mod kmer;
pub mod kmer_count;
pub mod kmer_index;
pub mod overlap;
pub mod ref_kmer_filter;
pub mod refbuild;
