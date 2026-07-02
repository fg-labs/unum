#![cfg(feature = "t1k-sys")]
//! Guards against silent drift between the shim's hand-copied
//! `nucToNum`/`numToNuc` tables (`crates/fg-t1k-sys/shim/shim.cpp`) and the
//! vendored source of truth (`vendor/t1k/Genotyper.cpp:37-42`).
//!
//! The lockstep differential tests in `diff_kmer.rs`, `diff_kmercode.rs`,
//! `diff_kmercount.rs`, and `diff_kmerindex.rs` only validate the Rust port
//! against the *shim's copy* of these tables --
//! if the shim's copy ever drifted from the real T1K source (e.g. a future
//! T1K vendor update changes the table and nobody updates the shim), those
//! tests would keep passing while silently testing against a stale table.
//! This test closes that gap by parsing the vendored `.cpp` file directly at
//! test time and comparing every entry against what the shim's FFI
//! (`fg_t1k_nuc_to_num`/`fg_t1k_num_to_nuc`) actually returns.

use std::path::PathBuf;

/// Reads `vendor/t1k/Genotyper.cpp` and extracts the brace-delimited,
/// comma-separated initializer values for the C array declared as
/// `array_name[...] = { ... };`, tolerating arbitrary whitespace/newlines
/// inside the braces (the source wraps `nucToNum`'s initializer across
/// several lines).
fn parse_vendor_array_values(source: &str, array_name: &str) -> Vec<i64> {
    // Find the declaration `array_name[` followed eventually by `=` then `{`.
    let decl_pos = source
        .find(&format!("{array_name}["))
        .unwrap_or_else(|| panic!("could not find `{array_name}[` declaration in vendor source"));
    let after_decl = &source[decl_pos..];

    let eq_pos = after_decl
        .find('=')
        .unwrap_or_else(|| panic!("could not find `=` after `{array_name}[` declaration"));
    let after_eq = &after_decl[eq_pos + 1..];

    let open_brace = after_eq
        .find('{')
        .unwrap_or_else(|| panic!("could not find `{{` opening `{array_name}`'s initializer"));
    let close_brace = after_eq[open_brace..]
        .find('}')
        .unwrap_or_else(|| panic!("could not find `}}` closing `{array_name}`'s initializer"));
    let body = &after_eq[open_brace + 1..open_brace + close_brace];

    body.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.parse::<i64>()
                .unwrap_or_else(|e| panic!("could not parse `{array_name}` entry {s:?}: {e}"))
        })
        .collect()
}

/// Path to the vendored T1K source that defines `nucToNum`/`numToNuc`
/// (relative to this crate's manifest, resolved via `CARGO_MANIFEST_DIR` so
/// the test works regardless of the invocation's current directory).
fn vendor_genotyper_cpp_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../vendor/t1k/Genotyper.cpp")
}

#[test]
fn shim_nuc_to_num_matches_vendor_table() {
    let path = vendor_genotyper_cpp_path();
    let source = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("could not read {}: {e}", path.display()));
    let vendor_values = parse_vendor_array_values(&source, "nucToNum");
    assert_eq!(
        vendor_values.len(),
        26,
        "expected 26 entries in vendor `nucToNum[26]`, got {}",
        vendor_values.len()
    );

    for (i, &expected) in vendor_values.iter().enumerate() {
        let shim_value = unsafe { fg_t1k_sys::fg_t1k_nuc_to_num(i32::try_from(i).unwrap()) };
        assert_eq!(
            i64::from(shim_value),
            expected,
            "nucToNum[{i}] mismatch: shim={shim_value} vendor={expected} \
             (shim.cpp's hand-copied table has drifted from {})",
            path.display()
        );
    }
}

#[test]
// `numToNuc` only ever holds ASCII letters ('A'/'C'/'G'/'T'), all < 0x80, so
// reinterpreting the FFI `c_char` byte pattern as `u8` never loses
// information here -- it's a byte reinterpretation, not a numeric truncation
// (same rationale as the `cast_possible_wrap` allow in fg-t1k-sys/src/lib.rs).
#[allow(clippy::cast_sign_loss)]
fn shim_num_to_nuc_matches_vendor_table() {
    let path = vendor_genotyper_cpp_path();
    let source = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("could not read {}: {e}", path.display()));

    // numToNuc's initializer holds character literals (e.g. 'A'), not plain
    // integers, so it needs its own tiny parser rather than
    // parse_vendor_array_values's integer parsing.
    let decl_pos =
        source.find("numToNuc[").expect("could not find `numToNuc[` declaration in vendor source");
    let after_decl = &source[decl_pos..];
    let eq_pos = after_decl.find('=').expect("could not find `=` after `numToNuc[` declaration");
    let after_eq = &after_decl[eq_pos + 1..];
    let open_brace =
        after_eq.find('{').expect("could not find `{` opening `numToNuc`'s initializer");
    let close_brace = after_eq[open_brace..]
        .find('}')
        .expect("could not find `}` closing `numToNuc`'s initializer");
    let body = &after_eq[open_brace + 1..open_brace + close_brace];

    let vendor_chars: Vec<char> = body
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| {
            // Each entry looks like `'A'`: strip the surrounding single quotes.
            let stripped = s.strip_prefix('\'').and_then(|s| s.strip_suffix('\''));
            let stripped =
                stripped.unwrap_or_else(|| panic!("could not parse `numToNuc` entry {s:?}"));
            stripped
                .chars()
                .next()
                .unwrap_or_else(|| panic!("empty `numToNuc` char literal: {s:?}"))
        })
        .collect();

    assert_eq!(
        vendor_chars.len(),
        4,
        "expected 4 entries in vendor `numToNuc[4]`, got {}",
        vendor_chars.len()
    );

    for (i, &expected) in vendor_chars.iter().enumerate() {
        let shim_value = unsafe { fg_t1k_sys::fg_t1k_num_to_nuc(i32::try_from(i).unwrap()) };
        assert_eq!(
            shim_value as u8 as char,
            expected,
            "numToNuc[{i}] mismatch: shim={shim_value} vendor={expected} \
             (shim.cpp's hand-copied table has drifted from {})",
            path.display()
        );
    }
}
