// Standard C headers the vendored T1K headers rely on but do not include
// themselves (they only compile in T1K because a .cpp includes these first):
// KmerCount.hpp uses strlen/atoi/fopen/fscanf/printf. Include them BEFORE the
// T1K headers so the shim compiles standalone on gcc/libstdc++ (Linux), not just
// clang/libc++ (macOS), which pulls them in transitively.
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "KmerCode.hpp"   // header-only; declares extern nucToNum/numToNuc
#include "KmerCount.hpp"  // header-only; depends on KmerCode.hpp above
#include "KmerIndex.hpp"  // header-only; depends on KmerCode.hpp + SimpleVector.hpp
#include "SeqSet.hpp"     // header-only; depends on KmerIndex.hpp + ReadFiles.hpp (zlib) + AlignAlgo.hpp
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
void* fg_t1k_kmercode_new(int k) {
    try {
        return new KmerCode(k);
    } catch (...) {
        return nullptr;
    }
}
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
    try {
        SimpleVector<struct _indexInfo>* list =
            static_cast<KmerIndex*>(idxp)->Search(*static_cast<KmerCode*>(kcp));
        struct _indexInfo entry = (*list)[i];
        *out_idx = entry.idx;
        *out_offset = entry.offset;
    } catch (...) {
        // No return value to signal failure through; on exception, write a
        // clear sentinel rather than leaving *out_idx/*out_offset
        // uninitialized. Defense-in-depth only: the Rust caller only ever
        // calls this for `i in 0..search_size` on an unmutated index (see
        // CppKmerIndex::search), so `Search`/operator[] are not expected to
        // throw in that usage.
        *out_idx = 0;
        *out_offset = 0;
    }
}

void fg_t1k_kmerindex_build_index_from_read(void* idxp, void* kcp, const char* s, int len, int id,
                                             int shift) {
    try {
        static_cast<KmerIndex*>(idxp)->BuildIndexFromRead(
            *static_cast<KmerCode*>(kcp), const_cast<char*>(s), len, id, shift);
    } catch (...) {
    }
}

// IsLowComplexity is a free function defined identically in
// FastqExtractor.cpp:89-111 and BamExtractor.cpp:144-166 (verbatim
// byte-for-byte copies of each other) -- it is NOT a SeqSet method, and
// neither of those two .cpp files is header-only (each has its own main()),
// so we cannot #include either one here (same reasoning as the
// nucToNum/numToNuc verbatim copies above: duplicate-symbol/duplicate-main
// linkage). This is a verbatim copy of that function body, not a
// reimplementation.
static bool IsLowComplexity(char* seq) {
    int cnt[5] = { 0, 0, 0, 0, 0 };
    int i;
    for (i = 0; seq[i]; ++i) {
        if (seq[i] == 'N')
            ++cnt[4];
        else
            ++cnt[nucToNum[seq[i] - 'A']];
    }

    if (cnt[0] >= i / 2 || cnt[1] >= i / 2 || cnt[2] >= i / 2 || cnt[3] >= i / 2 || cnt[4] >= i / 10)
        return true;

    int lowCnt = 0;
    for (i = 0; i < 4; ++i)
        if (cnt[i] <= 2)
            ++lowCnt;
    if (lowCnt >= 2)
        return true;
    return false;
}

// Opaque-handle FFI for SeqSet, scoped to the FastqExtractor/BamExtractor
// read-candidate-filtering slice (reference load, GetHitsFromRead,
// HasHitInSet, IsLowComplexity/IsGoodCandidate) -- see shim.h for the
// documented scope. SeqSet's constructor is lightweight (no large fixed-size
// allocation like KmerIndex's), but InputRefFa can throw (e.g. bad_alloc, or
// a ReadFiles/zlib failure on a missing/corrupt FASTA), so both the
// constructor and the reference-loading call are wrapped in try/catch,
// matching the exception-safety rule used throughout this shim: no C++
// exception may unwind across an `extern "C"` boundary.
void* fg_t1k_seqset_new(int kmerLength) {
    try {
        return new SeqSet(kmerLength);
    } catch (...) {
        return nullptr;
    }
}

void fg_t1k_seqset_free(void* p) { delete static_cast<SeqSet*>(p); }

int fg_t1k_seqset_load_ref(void* p, const char* fastaPath) {
    try {
        // InputRefFa takes `char*` (not `const char*`) but never mutates the
        // buffer itself (it only opens `filename` for reading via
        // ReadFiles::AddReadFile); const_cast is safe here.
        static_cast<SeqSet*>(p)->InputRefFa(const_cast<char*>(fastaPath));
        return 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_seqset_is_low_complexity(void* /*p*/, const char* seq) {
    return IsLowComplexity(const_cast<char*>(seq)) ? 1 : 0;
}

int fg_t1k_seqset_has_hit_in_set(void* p, const char* read) {
    try {
        // HasHitInSet writes the read's reverse complement into `rcRead`
        // (an output scratch buffer, not an input) via GetHitsFromRead's
        // internal ReverseComplement call; T1K's own callers pass a
        // 100001-byte stack buffer (FastqExtractor.cpp/BamExtractor.cpp), so
        // a same-sized buffer here is more than sufficient for any realistic
        // test read. Deliberately NOT `static`: this shim may be called
        // concurrently from multiple test threads, and a `static` buffer
        // would be shared mutable state across those calls (a data race);
        // a plain local array gives each call its own stack-allocated copy.
        char rcBuf[100001];
        // Guard the fixed-size scratch buffer: HasHitInSet reverse-complements
        // `read` into `rcBuf`, writing strlen(read) bytes plus a NUL terminator
        // (SeqSet::ReverseComplement), so a read that leaves no room for that
        // terminator would overflow the stack. Reject it via the same -1
        // sentinel the catch(...) path uses.
        if (strlen(read) >= sizeof(rcBuf)) {
            return -1;
        }
        return static_cast<SeqSet*>(p)->HasHitInSet(const_cast<char*>(read), rcBuf) ? 1 : 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_seqset_is_good_candidate(void* p, const char* read) {
    try {
        // Not `static` -- see fg_t1k_seqset_has_hit_in_set above.
        char rcBuf[100001];
        // Guard the fixed-size scratch buffer -- see fg_t1k_seqset_has_hit_in_set above.
        if (strlen(read) >= sizeof(rcBuf)) {
            return -1;
        }
        SeqSet* seqSet = static_cast<SeqSet*>(p);
        if (!IsLowComplexity(const_cast<char*>(read)) &&
            seqSet->HasHitInSet(const_cast<char*>(read), rcBuf)) {
            return 1;
        }
        return 0;
    } catch (...) {
        return -1;
    }
}
