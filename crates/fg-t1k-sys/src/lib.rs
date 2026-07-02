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

        /// Constructs a C++ `KmerCount(k)`. Returns NULL if construction
        /// throws (the shim catches every exception at this boundary); the
        /// caller MUST check for a NULL return before using the handle.
        pub fn fg_t1k_kmercount_new(k: i32) -> *mut c_void;
        pub fn fg_t1k_kmercount_free(p: *mut c_void);
        pub fn fg_t1k_kmercount_add_count(p: *mut c_void, read: *const c_char) -> i32;
        pub fn fg_t1k_kmercount_get_count(p: *mut c_void, kmer: *const c_char) -> i32;
        pub fn fg_t1k_kmercount_jaccard(a: *mut c_void, b: *mut c_void) -> f64;

        /// Constructs a C++ `KmerIndex()`. Returns NULL if construction
        /// throws (the shim catches every exception at this boundary); the
        /// caller MUST check for a NULL return before using the handle.
        pub fn fg_t1k_kmerindex_new() -> *mut c_void;
        pub fn fg_t1k_kmerindex_free(p: *mut c_void);
        pub fn fg_t1k_kmerindex_insert(
            idxp: *mut c_void,
            kcp: *mut c_void,
            idx: u32,
            offset: u32,
            strand: i32,
        );
        pub fn fg_t1k_kmerindex_remove(
            idxp: *mut c_void,
            kcp: *mut c_void,
            idx: u32,
            offset: u32,
            strand: i32,
        );
        pub fn fg_t1k_kmerindex_search_size(idxp: *mut c_void, kcp: *mut c_void) -> i32;
        pub fn fg_t1k_kmerindex_search_get(
            idxp: *mut c_void,
            kcp: *mut c_void,
            i: i32,
            out_idx: *mut u32,
            out_offset: *mut u32,
        );
        pub fn fg_t1k_kmerindex_build_index_from_read(
            idxp: *mut c_void,
            kcp: *mut c_void,
            s: *const c_char,
            len: i32,
            id: i32,
            shift: i32,
        );
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

/// Safe Rust wrapper around the opaque C++ `KmerCount*` handle.
///
/// Owns the handle for its lifetime: [`CppKmerCount::new`] allocates the C++
/// object via `fg_t1k_kmercount_new` and `Drop` calls `fg_t1k_kmercount_free`
/// exactly once. Construction checks the returned handle for NULL: the C++
/// constructor allocates heavy state (`new std::map<uint64_t,int>[103]`),
/// which can throw `std::bad_alloc`; the shim catches that at the `extern
/// "C"` boundary (an exception may not unwind across it) and returns NULL
/// instead, so this wrapper must -- and does -- refuse to proceed with a
/// NULL handle rather than silently dereferencing it later.
#[cfg(feature = "t1k-sys")]
pub struct CppKmerCount {
    handle: *mut std::os::raw::c_void,
}

#[cfg(feature = "t1k-sys")]
impl CppKmerCount {
    /// Constructs a new C++ `KmerCount` for k-mer length `k`.
    ///
    /// # Panics
    ///
    /// Panics if the underlying `fg_t1k_kmercount_new` call returns NULL
    /// (i.e. the C++ constructor threw an exception, most plausibly
    /// `std::bad_alloc`).
    #[must_use]
    pub fn new(k: i32) -> Self {
        let handle = unsafe { ffi::fg_t1k_kmercount_new(k) };
        assert!(
            !handle.is_null(),
            "fg_t1k_kmercount_new({k}) returned NULL: C++ KmerCount construction failed"
        );
        Self { handle }
    }

    /// Mirrors `KmerCount::AddCount`. `read` must not contain an interior NUL
    /// byte (mirrors the C++ side's `strlen`-based length calculation, which
    /// would silently truncate at the first NUL).
    ///
    /// # Panics
    ///
    /// Panics if `read` contains an interior NUL byte.
    pub fn add_count(&mut self, read: &[u8]) -> i32 {
        let c_read =
            std::ffi::CString::new(read).expect("read must not contain an interior NUL byte");
        unsafe { ffi::fg_t1k_kmercount_add_count(self.handle, c_read.as_ptr()) }
    }

    /// Mirrors `KmerCount::GetCount`. Only the first `k` bytes of `kmer` are
    /// read (matching the C++ side, which never calls `strlen` on the query
    /// k-mer); `kmer` must still not contain an interior NUL byte since it is
    /// passed across the FFI boundary as a NUL-terminated C string.
    ///
    /// # Panics
    ///
    /// Panics if `kmer` contains an interior NUL byte.
    pub fn get_count(&self, kmer: &[u8]) -> i32 {
        let c_kmer =
            std::ffi::CString::new(kmer).expect("kmer must not contain an interior NUL byte");
        unsafe { ffi::fg_t1k_kmercount_get_count(self.handle, c_kmer.as_ptr()) }
    }

    /// Mirrors `KmerCount::GetCountSimilarityJaccard`.
    #[must_use]
    pub fn jaccard(&self, other: &CppKmerCount) -> f64 {
        unsafe { ffi::fg_t1k_kmercount_jaccard(self.handle, other.handle) }
    }
}

#[cfg(feature = "t1k-sys")]
impl Drop for CppKmerCount {
    fn drop(&mut self) {
        unsafe { ffi::fg_t1k_kmercount_free(self.handle) }
    }
}

/// A single `(idx, offset)` entry, mirroring the FFI-exposed fields of the
/// C++ `_indexInfo` struct (`KmerIndex.hpp:12-17`). Used only by
/// [`CppKmerIndex::search`] to give differential tests an owned,
/// order-preserving snapshot of a C++ `Search` result.
#[cfg(feature = "t1k-sys")]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CppIndexInfo {
    pub idx: u32,
    pub offset: u32,
}

/// Safe Rust wrapper around the opaque C++ `KmerIndex*` handle.
///
/// Owns the handle for its lifetime: [`CppKmerIndex::new`] allocates the C++
/// object via `fg_t1k_kmerindex_new` and `Drop` calls `fg_t1k_kmerindex_free`
/// exactly once. Construction checks the returned handle for NULL: the C++
/// constructor allocates heavy state (`new std::map<uint64_t,
/// SimpleVector<_indexInfo>>[1000003]`), which can throw `std::bad_alloc`;
/// the shim catches that at the `extern "C"` boundary (an exception may not
/// unwind across it) and returns NULL instead, so this wrapper must -- and
/// does -- refuse to proceed with a NULL handle rather than silently
/// dereferencing it later.
///
/// Every method takes a [`CppKmerCode`] (the same opaque-handle wrapper used
/// for the KmerCode differential tests), so a test can drive one KmerCode
/// instance and feed it to both the Rust and C++ KmerIndex sides in
/// lockstep, guaranteeing identical `GetCode()`/`IsValid()` inputs.
#[cfg(feature = "t1k-sys")]
pub struct CppKmerIndex {
    handle: *mut std::os::raw::c_void,
}

#[cfg(feature = "t1k-sys")]
impl CppKmerIndex {
    /// Constructs a new C++ `KmerIndex`.
    ///
    /// # Panics
    ///
    /// Panics if the underlying `fg_t1k_kmerindex_new` call returns NULL
    /// (i.e. the C++ constructor threw an exception, most plausibly
    /// `std::bad_alloc`).
    #[must_use]
    pub fn new() -> Self {
        let handle = unsafe { ffi::fg_t1k_kmerindex_new() };
        assert!(
            !handle.is_null(),
            "fg_t1k_kmerindex_new() returned NULL: C++ KmerIndex construction failed"
        );
        Self { handle }
    }

    /// Mirrors `KmerIndex::Insert`. `strand` is accepted but ignored by the
    /// C++ side (the `strand` field of `_indexInfo` is commented out).
    pub fn insert(&mut self, kmer_code: &CppKmerCode, idx: u32, offset: u32, strand: i32) {
        unsafe { ffi::fg_t1k_kmerindex_insert(self.handle, kmer_code.handle, idx, offset, strand) }
    }

    /// Mirrors `KmerIndex::Remove`. `strand` is accepted but ignored by the
    /// C++ side, same as [`CppKmerIndex::insert`].
    pub fn remove(&mut self, kmer_code: &CppKmerCode, idx: u32, offset: u32, strand: i32) {
        unsafe { ffi::fg_t1k_kmerindex_remove(self.handle, kmer_code.handle, idx, offset, strand) }
    }

    /// Mirrors `KmerIndex::Search`, returning an owned, order-preserving
    /// snapshot of the matched `SimpleVector<_indexInfo>` (or an empty `Vec`
    /// if the k-mer is invalid or absent from the index -- matching the C++
    /// `nullHit` sentinel).
    ///
    /// # Panics
    ///
    /// Panics if `fg_t1k_kmerindex_search_size` returns a negative size
    /// (which would indicate a shim-level bug; `SimpleVector::Size()` cannot
    /// itself return a negative value).
    #[must_use]
    pub fn search(&self, kmer_code: &CppKmerCode) -> Vec<CppIndexInfo> {
        let size_i32 = unsafe { ffi::fg_t1k_kmerindex_search_size(self.handle, kmer_code.handle) };
        let size = usize::try_from(size_i32).unwrap_or_else(|_| {
            panic!("fg_t1k_kmerindex_search_size returned a negative size: {size_i32}")
        });
        let mut out = Vec::with_capacity(size);
        for i in 0..size_i32 {
            let mut idx: u32 = 0;
            let mut offset: u32 = 0;
            unsafe {
                ffi::fg_t1k_kmerindex_search_get(
                    self.handle,
                    kmer_code.handle,
                    i,
                    &mut idx,
                    &mut offset,
                );
            }
            out.push(CppIndexInfo { idx, offset });
        }
        out
    }

    /// Mirrors `KmerIndex::BuildIndexFromRead`. `s` need not be
    /// NUL-terminated; exactly `s.len()` bytes are read.
    ///
    /// # Panics
    ///
    /// Panics if `s.len()` does not fit in an `i32` (mirrors the C++ side's
    /// `int len` parameter).
    pub fn build_index_from_read(
        &mut self,
        kmer_code: &mut CppKmerCode,
        s: &[u8],
        id: i32,
        shift: i32,
    ) {
        let len = i32::try_from(s.len()).expect("read length must fit in an i32");
        unsafe {
            ffi::fg_t1k_kmerindex_build_index_from_read(
                self.handle,
                kmer_code.handle,
                s.as_ptr().cast::<std::os::raw::c_char>(),
                len,
                id,
                shift,
            );
        }
    }
}

#[cfg(feature = "t1k-sys")]
impl Default for CppKmerIndex {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(feature = "t1k-sys")]
impl Drop for CppKmerIndex {
    fn drop(&mut self) {
        unsafe { ffi::fg_t1k_kmerindex_free(self.handle) }
    }
}
