#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
uint64_t fg_t1k_canonical_kmer(const char* seq, int len, int k);

// Opaque-handle FFI for KmerCode: lets Rust drive a real C++ KmerCode
// instance step-by-step for lockstep differential testing. The handle is an
// opaque `void*` (really a `KmerCode*`); callers must free it exactly once
// via fg_t1k_kmercode_free.
void* fg_t1k_kmercode_new(int k);
void fg_t1k_kmercode_free(void* p);
void fg_t1k_kmercode_restart(void* p);
void fg_t1k_kmercode_append(void* p, char c);
void fg_t1k_kmercode_prepend(void* p, char c);
void fg_t1k_kmercode_shift_right(void* p, int k);
void fg_t1k_kmercode_set_code(void* p, uint64_t v);
uint64_t fg_t1k_kmercode_get_code(void* p);
uint64_t fg_t1k_kmercode_canonical(void* p);
uint64_t fg_t1k_kmercode_rc(void* p);
int fg_t1k_kmercode_is_valid(void* p);
int fg_t1k_kmercode_kmer_length(void* p);

// Exposes the shim's hand-copied nucToNum/numToNuc tables (see shim.cpp) so
// differential tests can assert they still match the vendored table in
// vendor/t1k/Genotyper.cpp, rather than only checking the shim's copy
// against itself. `i` must be in 0..26 for fg_t1k_nuc_to_num and 0..4 for
// fg_t1k_num_to_nuc.
signed char fg_t1k_nuc_to_num(int i);
char fg_t1k_num_to_nuc(int i);

// Opaque-handle FFI for KmerCount: lets Rust drive a real C++ KmerCount
// instance (AddCount/GetCount/Jaccard) for lockstep differential testing.
// The handle is an opaque `void*` (really a `KmerCount*`); callers must free
// it exactly once via fg_t1k_kmercount_free. `fg_t1k_kmercount_new` returns
// NULL if construction throws (e.g. std::bad_alloc from the heavy internal
// `std::map<uint64_t,int>[103]` allocation) -- callers MUST check for a NULL
// handle before using it.
void* fg_t1k_kmercount_new(int k);
void fg_t1k_kmercount_free(void* p);
// `read`/`kmer` must be NUL-terminated C strings (KmerCount::AddCount uses
// strlen(), and GetCount reads exactly `k` bytes starting at the pointer).
int fg_t1k_kmercount_add_count(void* p, const char* read);
int fg_t1k_kmercount_get_count(void* p, const char* kmer);
double fg_t1k_kmercount_jaccard(void* a, void* b);
#ifdef __cplusplus
}
#endif
