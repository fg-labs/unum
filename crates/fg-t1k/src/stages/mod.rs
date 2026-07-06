//! Pipeline stage implementations.
//!
//! Each stage is a self-contained, native-Rust port of the corresponding T1K
//! tool: [`build`] (reference build), [`extract`] (candidate-read extraction),
//! and [`genotype`] (genotype calling). There is no cross-stage orchestrator
//! here -- the vendored `run-t1k` wrapper's extract -> genotype -> analyze
//! pipeline was oracle-only (it shelled out to the vendored C++ binaries) and
//! was removed along with the oracle; each subcommand is invoked directly.

pub mod build;
pub mod extract;
pub mod genotype;
