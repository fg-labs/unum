#include "KmerCode.hpp"   // header-only; declares extern nucToNum/numToNuc
#include "KmerCount.hpp"  // header-only; depends on KmerCode.hpp above
#include "KmerIndex.hpp"  // header-only; depends on KmerCode.hpp + SimpleVector.hpp
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

// Opaque-handle FFI for KmerCount. KmerCount's constructor allocates heavy
// state (`new std::map<uint64_t,int>[103]`), which can throw std::bad_alloc;
// letting an exception unwind through this extern "C" boundary is undefined
// behavior, so construction is wrapped in try/catch and returns NULL on
// failure. Every other KmerCount operation used here (AddCount/GetCount/
// GetCountSimilarityJaccard) only touches std::map insert/lookup, which can
// also throw (bad_alloc, or an out_of_range from map::at -- not used here),
// but per the same exception-safety rule we do not let those propagate
// across the boundary either.
void* fg_t1k_kmercount_new(int k) {
    try {
        return new KmerCount(k);
    } catch (...) {
        return nullptr;
    }
}

void fg_t1k_kmercount_free(void* p) { delete static_cast<KmerCount*>(p); }

int fg_t1k_kmercount_add_count(void* p, const char* read) {
    try {
        // KmerCount::AddCount takes `char*` (not `const char*`) but never
        // mutates the buffer; const_cast is safe here.
        return static_cast<KmerCount*>(p)->AddCount(const_cast<char*>(read));
    } catch (...) {
        return -1;
    }
}

int fg_t1k_kmercount_get_count(void* p, const char* kmer) {
    try {
        return static_cast<KmerCount*>(p)->GetCount(const_cast<char*>(kmer));
    } catch (...) {
        return -1;
    }
}

double fg_t1k_kmercount_jaccard(void* a, void* b) {
    try {
        return static_cast<KmerCount*>(a)->GetCountSimilarityJaccard(
            *static_cast<const KmerCount*>(b));
    } catch (...) {
        return -1.0;
    }
}

// Opaque-handle FFI for KmerIndex. KmerIndex's constructor allocates heavy
// state (`new std::map<uint64_t, SimpleVector<_indexInfo>>[1000003]`), which
// can throw std::bad_alloc; letting an exception unwind through this extern
// "C" boundary is undefined behavior, so construction is wrapped in
// try/catch and returns NULL on failure (mirrors fg_t1k_kmercount_new
// above). Insert/Remove/Search only touch std::map/SimpleVector operations
// that can also throw (bad_alloc), so those are guarded too, even though a
// throw there is not expected in practice.
void* fg_t1k_kmerindex_new(void) {
    try {
        return new KmerIndex();
    } catch (...) {
        return nullptr;
    }
}

void fg_t1k_kmerindex_free(void* p) { delete static_cast<KmerIndex*>(p); }

void fg_t1k_kmerindex_insert(void* idxp, void* kcp, uint32_t idx, uint32_t offset, int strand) {
    try {
        static_cast<KmerIndex*>(idxp)->Insert(*static_cast<KmerCode*>(kcp), idx, offset, strand);
    } catch (...) {
        // No return value to signal failure through; Insert has no
        // observable failure mode in the differential test's usage.
    }
}

void fg_t1k_kmerindex_remove(void* idxp, void* kcp, uint32_t idx, uint32_t offset, int strand) {
    try {
        static_cast<KmerIndex*>(idxp)->Remove(*static_cast<KmerCode*>(kcp), idx, offset, strand);
    } catch (...) {
    }
}

int fg_t1k_kmerindex_search_size(void* idxp, void* kcp) {
    try {
        SimpleVector<struct _indexInfo>* list =
            static_cast<KmerIndex*>(idxp)->Search(*static_cast<KmerCode*>(kcp));
        return list->Size();
    } catch (...) {
        return -1;
    }
}

void fg_t1k_kmerindex_search_get(void* idxp, void* kcp, int i, uint32_t* out_idx,
                                  uint32_t* out_offset) {
    SimpleVector<struct _indexInfo>* list =
        static_cast<KmerIndex*>(idxp)->Search(*static_cast<KmerCode*>(kcp));
    struct _indexInfo entry = (*list)[i];
    *out_idx = entry.idx;
    *out_offset = entry.offset;
}

void fg_t1k_kmerindex_build_index_from_read(void* idxp, void* kcp, const char* s, int len, int id,
                                             int shift) {
    try {
        static_cast<KmerIndex*>(idxp)->BuildIndexFromRead(
            *static_cast<KmerCode*>(kcp), const_cast<char*>(s), len, id, shift);
    } catch (...) {
    }
}
