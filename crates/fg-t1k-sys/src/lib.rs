#![deny(warnings)]
#![allow(clippy::must_use_candidate)]
//! fg-t1k-sys: all C++ contact. FFI shims + vendored T1K build. Dev/test only.
//! `unsafe` is permitted here (FFI); nowhere else in the workspace.

// oracle module and FFI decls are feature-gated (they depend on build-script env vars).
#[cfg(feature = "t1k-sys")]
pub mod oracle;

/// Raw FFI declarations for the header-only shim (`shim/shim.cpp`).
///
/// The shim includes T1K headers (e.g. `KmerCode.hpp`) but links zero T1K
/// `.cpp` files, proving that individual header-only components can be
/// exercised for differential testing without pulling in T1K's `main()`.
#[cfg(feature = "t1k-sys")]
mod ffi {
    use std::os::raw::{c_char, c_void};
    unsafe extern "C" {
        pub fn fg_t1k_canonical_kmer(seq: *const c_char, len: i32, k: i32) -> u64;

        pub fn fg_t1k_kmercode_new(k: i32) -> *mut c_void;
        pub fn fg_t1k_kmercode_free(p: *mut c_void);
        pub fn fg_t1k_kmercode_restart(p: *mut c_void);
        pub fn fg_t1k_kmercode_append(p: *mut c_void, c: c_char);
        pub fn fg_t1k_kmercode_prepend(p: *mut c_void, c: c_char);
        pub fn fg_t1k_kmercode_shift_right(p: *mut c_void, k: i32);
        pub fn fg_t1k_kmercode_set_code(p: *mut c_void, v: u64);
        pub fn fg_t1k_kmercode_get_code(p: *mut c_void) -> u64;
        pub fn fg_t1k_kmercode_canonical(p: *mut c_void) -> u64;
        pub fn fg_t1k_kmercode_rc(p: *mut c_void) -> u64;
        pub fn fg_t1k_kmercode_is_valid(p: *mut c_void) -> i32;
        pub fn fg_t1k_kmercode_kmer_length(p: *mut c_void) -> i32;

        /// Returns the shim's own `nucToNum[i]` entry (`i` in `0..26`). Bound
        /// as `i8` (not `c_char`) because the C signature is explicitly
        /// `signed char`, regardless of the platform's default char
        /// signedness.
        pub fn fg_t1k_nuc_to_num(i: i32) -> i8;
        /// Returns the shim's own `numToNuc[i]` entry (`i` in `0..4`).
        pub fn fg_t1k_num_to_nuc(i: i32) -> c_char;
    }
}
#[cfg(feature = "t1k-sys")]
pub use ffi::fg_t1k_canonical_kmer;
#[cfg(feature = "t1k-sys")]
pub use ffi::{fg_t1k_nuc_to_num, fg_t1k_num_to_nuc};

/// Safe Rust wrapper around the opaque C++ `KmerCode*` handle.
///
/// Owns the handle for its lifetime: [`CppKmerCode::new`] allocates the C++
/// object via `fg_t1k_kmercode_new`, and `Drop` calls `fg_t1k_kmercode_free`
/// exactly once. All methods forward to the corresponding FFI function; the
/// `unsafe` calls are sound because `self.handle` is a non-null pointer
/// created by `fg_t1k_kmercode_new` and never freed until `Drop`.
#[cfg(feature = "t1k-sys")]
pub struct CppKmerCode {
    handle: *mut std::os::raw::c_void,
}

#[cfg(feature = "t1k-sys")]
impl CppKmerCode {
    /// Constructs a new C++ `KmerCode` for k-mer length `k`.
    #[must_use]
    pub fn new(k: i32) -> Self {
        let handle = unsafe { ffi::fg_t1k_kmercode_new(k) };
        Self { handle }
    }

    /// Mirrors `KmerCode::Restart`.
    pub fn restart(&mut self) {
        unsafe { ffi::fg_t1k_kmercode_restart(self.handle) }
    }

    /// Mirrors `KmerCode::Append`.
    // `c as c_char` reinterprets the raw byte pattern to match C++'s `char`
    // parameter exactly (whether `c_char` is signed or unsigned is
    // platform-dependent, same as C++'s `char`); this is not a numeric
    // truncation, just crossing the FFI boundary.
    #[allow(clippy::cast_possible_wrap)]
    pub fn append(&mut self, c: u8) {
        unsafe { ffi::fg_t1k_kmercode_append(self.handle, c as std::os::raw::c_char) }
    }

    /// Mirrors `KmerCode::Prepend`.
    #[allow(clippy::cast_possible_wrap)]
    pub fn prepend(&mut self, c: u8) {
        unsafe { ffi::fg_t1k_kmercode_prepend(self.handle, c as std::os::raw::c_char) }
    }

    /// Mirrors `KmerCode::ShiftRight`.
    pub fn shift_right(&mut self, k: i32) {
        unsafe { ffi::fg_t1k_kmercode_shift_right(self.handle, k) }
    }

    /// Mirrors `KmerCode::SetCode`.
    pub fn set_code(&mut self, v: u64) {
        unsafe { ffi::fg_t1k_kmercode_set_code(self.handle, v) }
    }

    /// Mirrors `KmerCode::GetCode`.
    #[must_use]
    pub fn get_code(&self) -> u64 {
        unsafe { ffi::fg_t1k_kmercode_get_code(self.handle) }
    }

    /// Mirrors `KmerCode::GetCanonicalKmerCode`.
    #[must_use]
    pub fn canonical(&self) -> u64 {
        unsafe { ffi::fg_t1k_kmercode_canonical(self.handle) }
    }

    /// Mirrors `KmerCode::GetReverseComplementCode`.
    #[must_use]
    pub fn rc(&self) -> u64 {
        unsafe { ffi::fg_t1k_kmercode_rc(self.handle) }
    }

    /// Mirrors `KmerCode::IsValid`.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        unsafe { ffi::fg_t1k_kmercode_is_valid(self.handle) != 0 }
    }

    /// Mirrors `KmerCode::GetKmerLength`.
    #[must_use]
    pub fn kmer_length(&self) -> i32 {
        unsafe { ffi::fg_t1k_kmercode_kmer_length(self.handle) }
    }
}

#[cfg(feature = "t1k-sys")]
impl Drop for CppKmerCode {
    fn drop(&mut self) {
        unsafe { ffi::fg_t1k_kmercode_free(self.handle) }
    }
}
