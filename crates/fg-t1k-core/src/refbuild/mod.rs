//! Rust port of T1K's reference-build pipeline (`vendor/t1k/t1k-build.pl` and the
//! scripts it drives), which turns an IPD-IMGT/HLA or IPD-KIR `.dat` database into
//! the FASTA references T1K's genotyper indexes.

pub mod dat;
pub mod gene_coord;
