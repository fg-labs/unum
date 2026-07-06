// This module is a submodule of every golden test crate, but each individual
// test crate uses only the subset of helpers it needs, so the unused ones are
// legitimately dead in that crate's compilation.
#![allow(dead_code)]
//! Shared helpers for the self-contained golden-file tests that replaced the
//! retired T1K-oracle FFI differential harness.
//!
//! Each former differential drove the Rust port and a real C++ T1K oracle
//! over identical, deterministically-generated inputs and asserted they
//! agreed. When the oracle was retired, the Rust port's output (which was, at
//! that point, byte-identical to the oracle -- every differential passed) was
//! frozen into the committed golden files this module reads. The tests keep
//! their original deterministic input generators verbatim; only the "compare
//! against the C++ oracle" step is replaced by "compare against the frozen
//! golden".
//!
//! # Update mode
//!
//! Setting `UPDATE_GOLDENS=1` in the environment makes [`Golden::finish`]
//! (re)write the golden file from the current Rust output instead of
//! asserting against it. This is the ONLY supported way to regenerate a
//! golden, and must only be used deliberately (a golden change is a
//! port-behavior change and should be reviewed as such). Under normal `cargo
//! test`/`cargo nextest` runs the variable is unset and the goldens are
//! asserted read-only.

use std::collections::BTreeMap;
use std::path::PathBuf;

/// A keyed golden file: an ordered `label -> serialized-value` map persisted
/// as a simple `label<TAB>value` text file (one record per line, labels
/// unique and sorted). Tests record each case's Rust output via [`record`]
/// keyed by the SAME deterministic label the input generator produced, then
/// call [`finish`] to either assert every recorded value matches the frozen
/// golden (normal mode) or rewrite the golden (`UPDATE_GOLDENS=1`).
///
/// [`record`]: Golden::record
/// [`finish`]: Golden::finish
pub struct Golden {
    path: PathBuf,
    recorded: BTreeMap<String, String>,
}

impl Golden {
    /// Opens the golden named `name` (resolved under this crate's
    /// `tests/golden/` directory, so it is independent of the process's
    /// current working directory).
    pub fn open(name: &str) -> Self {
        let path =
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests").join("golden").join(name);
        Self { path, recorded: BTreeMap::new() }
    }

    /// Records `value` under `label`. Panics if the same label is recorded
    /// twice (a duplicate label means the input generator produced a
    /// non-unique case key, which would silently drop a golden entry).
    pub fn record(&mut self, label: impl Into<String>, value: impl Into<String>) {
        let label = label.into();
        assert!(
            !label.contains('\t') && !label.contains('\n'),
            "golden label must not contain a tab or newline: {label:?}"
        );
        let value = value.into();
        assert!(!value.contains('\n'), "golden value must not contain a newline: {value:?}");
        let prev = self.recorded.insert(label.clone(), value);
        assert!(prev.is_none(), "duplicate golden label recorded: {label:?}");
    }

    /// In `UPDATE_GOLDENS=1` mode, (re)writes the golden file from the
    /// recorded values. Otherwise, asserts the recorded values are IDENTICAL
    /// to the frozen golden (same labels, same serialized value for each) --
    /// panicking with a precise diff on any mismatch.
    pub fn finish(self) {
        let serialized =
            self.recorded.iter().map(|(k, v)| format!("{k}\t{v}")).collect::<Vec<_>>().join("\n");

        if matches!(std::env::var("UPDATE_GOLDENS").as_deref(), Ok("1")) {
            let serialized = format!("{serialized}\n");
            std::fs::write(&self.path, serialized)
                .unwrap_or_else(|e| panic!("writing golden {:?}: {e}", self.path));
            return;
        }

        let contents = std::fs::read_to_string(&self.path).unwrap_or_else(|e| {
            panic!("reading golden {:?}: {e}\n(run with UPDATE_GOLDENS=1 to create it)", self.path)
        });
        let mut golden: BTreeMap<&str, &str> = BTreeMap::new();
        for line in contents.lines() {
            if line.is_empty() {
                continue;
            }
            let (label, value) = line
                .split_once('\t')
                .unwrap_or_else(|| panic!("malformed golden line (no tab): {line:?}"));
            let prev = golden.insert(label, value);
            assert!(
                prev.is_none(),
                "duplicate label {label:?} in golden {:?} (corrupt fixture; run with UPDATE_GOLDENS=1 to regenerate)",
                self.path
            );
        }

        // Every recorded case must be present in the golden with the same value.
        for (label, value) in &self.recorded {
            match golden.get(label.as_str()) {
                Some(&golden_value) => assert_eq!(
                    value.as_str(),
                    golden_value,
                    "golden mismatch for label {label:?} in {:?}\n  recorded: {value}\n  golden:   {golden_value}",
                    self.path
                ),
                None => panic!(
                    "recorded label {label:?} is missing from golden {:?} \
                     (run with UPDATE_GOLDENS=1 to regenerate)",
                    self.path
                ),
            }
        }
        // And the golden must not contain stale labels the test no longer produces.
        for label in golden.keys() {
            assert!(
                self.recorded.contains_key(*label),
                "golden {:?} contains stale label {label:?} not produced by the test \
                 (run with UPDATE_GOLDENS=1 to regenerate)",
                self.path
            );
        }
    }
}

/// Reads a byte-golden file (used by the extract/bam-extract tests, whose
/// golden is the emitted FASTQ bytes) under `tests/golden/`. Returns the
/// resolved path plus, in normal mode, the frozen bytes to compare against.
pub fn byte_golden_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests").join("golden").join(name)
}

/// Compares `produced` bytes against the byte-golden `name`. In
/// `UPDATE_GOLDENS=1` mode, (re)writes the golden instead.
pub fn assert_byte_golden(name: &str, produced: &[u8]) {
    let path = byte_golden_path(name);
    if matches!(std::env::var("UPDATE_GOLDENS").as_deref(), Ok("1")) {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).unwrap_or_else(|e| panic!("mkdir {parent:?}: {e}"));
        }
        std::fs::write(&path, produced).unwrap_or_else(|e| panic!("writing golden {path:?}: {e}"));
        return;
    }
    let golden = std::fs::read(&path).unwrap_or_else(|e| {
        panic!("reading byte golden {path:?}: {e}\n(run with UPDATE_GOLDENS=1 to create it)")
    });
    assert_eq!(
        produced,
        golden.as_slice(),
        "byte-golden mismatch for {path:?} (produced {} bytes, golden {} bytes)",
        produced.len(),
        golden.len()
    );
}
