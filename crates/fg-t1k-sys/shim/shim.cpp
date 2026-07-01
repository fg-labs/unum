#include "KmerCode.hpp"   // header-only; declares extern nucToNum/numToNuc
#include "shim.h"

// nucToNum/numToNuc are extern in the headers, defined only in the T1K .cpp files
// (which we cannot link -- duplicate main/symbols). Define them here verbatim from
// Genotyper.cpp:37-42 so the shim links standalone.
signed char nucToNum[26] = { 0, -1, 1, -1, -1, -1, 2,
    -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, 3,
    -1, -1, -1, -1, -1, -1 };
char numToNuc[4] = {'A', 'C', 'G', 'T'};

uint64_t fg_t1k_canonical_kmer(const char* seq, int len, int k) {
    KmerCode kc(k);
    for (int i = 0; i < len; ++i) {
        kc.Append(seq[i]);
    }
    return kc.GetCanonicalKmerCode();
}

// Opaque-handle FFI for KmerCode. Each function casts the opaque `void*`
// back to `KmerCode*` and forwards to the corresponding method, so Rust can
// drive a real C++ KmerCode instance step-by-step (see shim.h).
void* fg_t1k_kmercode_new(int k) { return new KmerCode(k); }
void fg_t1k_kmercode_free(void* p) { delete static_cast<KmerCode*>(p); }
void fg_t1k_kmercode_restart(void* p) { static_cast<KmerCode*>(p)->Restart(); }
void fg_t1k_kmercode_append(void* p, char c) { static_cast<KmerCode*>(p)->Append(c); }
void fg_t1k_kmercode_prepend(void* p, char c) { static_cast<KmerCode*>(p)->Prepend(c); }
void fg_t1k_kmercode_shift_right(void* p, int k) { static_cast<KmerCode*>(p)->ShiftRight(k); }
void fg_t1k_kmercode_set_code(void* p, uint64_t v) { static_cast<KmerCode*>(p)->SetCode(v); }
uint64_t fg_t1k_kmercode_get_code(void* p) { return static_cast<KmerCode*>(p)->GetCode(); }
uint64_t fg_t1k_kmercode_canonical(void* p) {
    return static_cast<KmerCode*>(p)->GetCanonicalKmerCode();
}
uint64_t fg_t1k_kmercode_rc(void* p) {
    return static_cast<KmerCode*>(p)->GetReverseComplementCode();
}
int fg_t1k_kmercode_is_valid(void* p) { return static_cast<KmerCode*>(p)->IsValid() ? 1 : 0; }
int fg_t1k_kmercode_kmer_length(void* p) {
    return static_cast<KmerCode*>(p)->GetKmerLength();
}

// Exposes the shim's own nucToNum/numToNuc tables (defined above) so Rust
// tests can assert them against the vendored copy in Genotyper.cpp, guarding
// against undetected drift between the two.
signed char fg_t1k_nuc_to_num(int i) { return nucToNum[i]; }
char fg_t1k_num_to_nuc(int i) { return numToNuc[i]; }
