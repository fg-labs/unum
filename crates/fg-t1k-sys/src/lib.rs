#![deny(warnings)]
#![allow(clippy::must_use_candidate)]
//! fg-t1k-sys: all C++ contact. FFI shims + vendored T1K build. Dev/test only.
//! `unsafe` is permitted here (FFI); nowhere else in the workspace.

// oracle module and FFI decls are feature-gated (they depend on build-script env vars).
#[cfg(feature = "t1k-sys")]
pub mod oracle;
