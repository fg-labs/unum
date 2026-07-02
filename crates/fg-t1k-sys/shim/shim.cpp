// Standard C headers the vendored T1K headers rely on but do not include
// themselves (they only compile in T1K because a .cpp includes these first):
// KmerCount.hpp uses strlen/atoi/fopen/fscanf/printf. Include them BEFORE the
// T1K headers so the shim compiles standalone on gcc/libstdc++ (Linux), not just
// clang/libc++ (macOS), which pulls them in transitively.
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "KmerCode.hpp"   // header-only; declares extern nucToNum/numToNuc
#include "KmerCount.hpp"  // header-only; depends on KmerCode.hpp above
#include "KmerIndex.hpp"  // header-only; depends on KmerCode.hpp + SimpleVector.hpp

// SeqSet::GetAlignStats (SeqSet.hpp:438-453) is declared `private`, but this
// shim needs to call it unmodified (not reimplemented) to serve as a real
// differential oracle for fg_t1k_core::align_algo::get_align_stats.
// Genotyper::ParseAlleleName/ReadAssignmentWeight (Genotyper.hpp:63,205) are
// likewise `private`, needed for the Task-5a differential (diff_genotyper_
// model.rs). `#define private public` around ONLY these includes is the
// standard, minimally-invasive technique for a test/differential harness to
// reach a private member without editing the vendored header at all -- it
// only affects how THIS translation unit's compiler front-end parses
// SeqSet's/Genotyper's access specifiers; it does not change
// vendor/t1k/SeqSet.hpp or vendor/t1k/Genotyper.hpp on disk, and every other
// translation unit (the real T1K binaries, any other shim TU) still sees
// `private` as private.
#define private public
#include "SeqSet.hpp"     // header-only; depends on KmerIndex.hpp + ReadFiles.hpp (zlib) + AlignAlgo.hpp
                           // (AlignAlgo.hpp is pulled in transitively by SeqSet.hpp; both
                           // AlignAlgo::* and struct _posWeight are usable directly)
#include "Genotyper.hpp"  // header-only; depends on SeqSet.hpp + KmerCount.hpp + SimpleVector.hpp above
#undef private

#include "alignments.hpp" // header-only; depends on samtools-0.1.19/sam.h (NOT built with -DHTSLIB,
                           // matching how build.rs compiles the bam-extractor oracle binary against
                           // the bundled samtools-0.1.19 -lbam, not htslib-1.15.1)
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
        return static_cast<SeqSet*>(p)->HasHitInSet(const_cast<char*>(read), rcBuf) ? 1 : 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_seqset_is_good_candidate(void* p, const char* read) {
    try {
        // Not `static` -- see fg_t1k_seqset_has_hit_in_set above.
        char rcBuf[100001];
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

int fg_t1k_seqset_get_overlaps_from_read(void* p, const char* read, int out_capacity,
                                          int* out_seq_idx, int* out_strand,
                                          int* out_read_start, int* out_read_end,
                                          int* out_seq_start, int* out_seq_end,
                                          int* out_match_cnt, double* out_similarity,
                                          int* out_count) {
    try {
        SeqSet* seqSet = static_cast<SeqSet*>(p);
        std::vector<struct _overlap> overlaps;
        // GetOverlapsFromRead takes non-const char* but never mutates
        // `read` itself (it only reads it and writes into its own
        // internally-allocated rcRead buffer); const_cast is safe here.
        // strand=0 (search both strands), barcode=-1 (no barcode
        // filtering) -- matching RefKmerFilter::get_overlaps_from_read's
        // own fixed parameters (see shim.h).
        int overlapCnt =
            seqSet->GetOverlapsFromRead(const_cast<char*>(read), 0, -1, overlaps);
        *out_count = overlapCnt;
        if (overlapCnt <= 0) {
            return 0;
        }
        int n = overlapCnt < out_capacity ? overlapCnt : out_capacity;
        for (int i = 0; i < n; ++i) {
            out_seq_idx[i] = overlaps[i].seqIdx;
            out_strand[i] = overlaps[i].strand;
            out_read_start[i] = overlaps[i].readStart;
            out_read_end[i] = overlaps[i].readEnd;
            out_seq_start[i] = overlaps[i].seqStart;
            out_seq_end[i] = overlaps[i].seqEnd;
            out_match_cnt[i] = overlaps[i].matchCnt;
            out_similarity[i] = overlaps[i].similarity;
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

// Opaque-handle FFI for Alignments (vendor/t1k/alignments.hpp), scoped to
// the BamExtractor.cpp slice -- see shim.h for the documented scope. Every
// entry point is wrapped in try/catch per this shim's exception-safety rule
// (no C++ exception may unwind across an `extern "C"` boundary), even though
// Alignments' constructor itself is lightweight (no heap allocation that can
// throw, unlike KmerCount/KmerIndex above).
void* fg_t1k_alignments_new(void) {
    try {
        return new Alignments();
    } catch (...) {
        return nullptr;
    }
}

void fg_t1k_alignments_free(void* p) { delete static_cast<Alignments*>(p); }

int fg_t1k_alignments_open(void* p, const char* path) {
    try {
        // Alignments::Open(char*) takes a non-const char* but only ever
        // strcpy()s it into its own internal fileName buffer; const_cast is
        // safe here. NOTE: the vendored Open() calls exit(1) (not a C++
        // throw) on a hard open failure, which this try/catch CANNOT
        // intercept -- callers must only pass paths known to exist and
        // parse as BAM/SAM (see shim.h).
        static_cast<Alignments*>(p)->Open(const_cast<char*>(path));
        return 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_rewind(void* p) {
    try {
        static_cast<Alignments*>(p)->Rewind();
        return 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_next(void* p) {
    try {
        return static_cast<Alignments*>(p)->Next();
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_get_read_seq(void* p, char* buffer) {
    try {
        static_cast<Alignments*>(p)->GetReadSeq(buffer);
        return 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_get_qual(void* p, char* buffer) {
    try {
        static_cast<Alignments*>(p)->GetQual(buffer);
        return 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_get_read_id(void* p, char* buffer) {
    try {
        strcpy(buffer, static_cast<Alignments*>(p)->GetReadId());
        return 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_is_first_mate(void* p) {
    try {
        return static_cast<Alignments*>(p)->IsFirstMate() ? 1 : 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_is_reverse(void* p) {
    try {
        return static_cast<Alignments*>(p)->IsReverse() ? 1 : 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_is_mate_reverse(void* p) {
    try {
        return static_cast<Alignments*>(p)->IsMateReverse() ? 1 : 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_is_aligned(void* p) {
    try {
        return static_cast<Alignments*>(p)->IsAligned() ? 1 : 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_is_template_aligned(void* p) {
    try {
        return static_cast<Alignments*>(p)->IsTemplateAligned() ? 1 : 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_is_primary(void* p) {
    try {
        return static_cast<Alignments*>(p)->IsPrimary() ? 1 : 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_get_chrom_id(void* p) {
    try {
        return static_cast<Alignments*>(p)->GetChromId();
    } catch (...) {
        return INT32_MIN;
    }
}

int fg_t1k_alignments_seg_count(void* p) {
    try {
        return static_cast<Alignments*>(p)->segCnt;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_seg(void* p, int i, int64_t* out_a, int64_t* out_b) {
    try {
        Alignments* a = static_cast<Alignments*>(p);
        *out_a = a->segments[i].a;
        *out_b = a->segments[i].b;
        return 0;
    } catch (...) {
        return -1;
    }
}

int fg_t1k_alignments_general_info(void* p, int stop_early, int* out_frag_stdev,
                                    int* out_read_len) {
    try {
        Alignments* a = static_cast<Alignments*>(p);
        a->GetGeneralInfo(stop_early != 0);
        *out_frag_stdev = a->fragStdev;
        *out_read_len = a->readLen;
        return 0;
    } catch (...) {
        return -1;
    }
}

// Free-function FFI for AlignAlgo (vendor/t1k/AlignAlgo.hpp) -- see shim.h
// for the documented scope (GlobalAlignment, GlobalAlignment_PosWeight, plus
// SeqSet::GetAlignStats). AlignAlgo::GlobalAlignment*  are `static`, so they
// are called directly with no instance; GetAlignStats is a non-static SeqSet
// member that touches no instance state, so a throwaway SeqSet is
// constructed to call it (SeqSet's constructor is lightweight, matching
// fg_t1k_seqset_new above).
int fg_t1k_alignalgo_global_alignment(const char* t, int lent, const char* p, int lenp, int band,
                                       signed char* out_align, int* out_len) {
    try {
        // GlobalAlignment takes non-const char* but never mutates t/p;
        // const_cast is safe here.
        int score = AlignAlgo::GlobalAlignment(const_cast<char*>(t), lent, const_cast<char*>(p),
                                                lenp, out_align, band);
        int len = 0;
        while (out_align[len] != -1) {
            ++len;
        }
        *out_len = len;
        return score;
    } catch (...) {
        *out_len = 0;
        return INT32_MIN;
    }
}

int fg_t1k_alignalgo_global_alignment_pos_weight(const int* t_weights, int lent, const char* p,
                                                   int lenp, signed char* out_align,
                                                   int* out_len) {
    try {
        std::vector<struct _posWeight> weights(lent);
        for (int i = 0; i < lent; ++i) {
            weights[i].count[0] = t_weights[4 * i + 0];
            weights[i].count[1] = t_weights[4 * i + 1];
            weights[i].count[2] = t_weights[4 * i + 2];
            weights[i].count[3] = t_weights[4 * i + 3];
        }
        // GlobalAlignment_PosWeight returns double but only ever produces
        // integer values (see align_algo.rs module docs); the shim surfaces
        // it as an int score to match GlobalAlignment's return type and the
        // Rust side's AlignResult::score: i32.
        double scoreD = AlignAlgo::GlobalAlignment_PosWeight(
            weights.data(), lent, const_cast<char*>(p), lenp, out_align);
        int len = 0;
        while (out_align[len] != -1) {
            ++len;
        }
        *out_len = len;
        return static_cast<int>(scoreD);
    } catch (...) {
        *out_len = 0;
        return INT32_MIN;
    }
}

int fg_t1k_seqset_get_align_stats(const signed char* align, int align_len, int update,
                                   int* out_match_cnt, int* out_mismatch_cnt,
                                   int* out_indel_cnt) {
    try {
        // GetAlignStats scans until it hits a -1 sentinel; build a
        // NUL(-1)-terminated copy since the Rust side passes an
        // un-terminated slice + explicit length.
        std::vector<signed char> buf(align, align + align_len);
        buf.push_back(-1);

        SeqSet seqSet(0);  // kmerLength is irrelevant; GetAlignStats touches no instance state.
        int matchCnt = *out_match_cnt;
        int mismatchCnt = *out_mismatch_cnt;
        int indelCnt = *out_indel_cnt;
        seqSet.GetAlignStats(buf.data(), update != 0, matchCnt, mismatchCnt, indelCnt);
        *out_match_cnt = matchCnt;
        *out_mismatch_cnt = mismatchCnt;
        *out_indel_cnt = indelCnt;
        return 0;
    } catch (...) {
        return -1;
    }
}

// Free-function FFI for Genotyper's Task-5a deterministic slice (see shim.h
// for the documented scope): ParseAlleleName and ReadAssignmentWeight. Both
// are private members reached via the `#define private public` trick above
// (see that comment for why this is safe/standard). Genotyper's constructor
// (Genotyper(int kmerLength)) is lightweight enough that a fresh, throwaway
// instance is constructed per call rather than exposing an opaque handle --
// neither method touches any Genotyper state beyond
// alleleDigitUnits/alleleDelimiter (ParseAlleleName) or refSet's
// refSeqSimilarity (ReadAssignmentWeight), both of which this shim sets
// explicitly from the caller's arguments before calling.
int fg_t1k_genotyper_parse_allele_name(const char* allele, int fieldsType, int alleleDigitUnits,
                                        char alleleDelimiter, char* out_gene, char* out_major) {
    try {
        Genotyper genotyper(11);
        genotyper.SetAlleleNameStructure(alleleDigitUnits, alleleDelimiter);
        // ParseAlleleName takes non-const char* but never mutates `allele`
        // itself (SetAlleleNameStructure it is not); const_cast is safe
        // here, matching the pattern used throughout this shim.
        genotyper.ParseAlleleName(const_cast<char*>(allele), out_gene, out_major, fieldsType);
        return 0;
    } catch (...) {
        return -1;
    }
}

double fg_t1k_genotyper_read_assignment_weight(double similarity, int hasN,
                                                double refSeqSimilarity) {
    try {
        Genotyper genotyper(11);
        genotyper.refSet.SetRefSeqSimilarity(refSeqSimilarity);

        struct _fragmentOverlap o;
        memset(&o, 0, sizeof(o));
        o.similarity = similarity;
        o.hasN = hasN != 0;

        return genotyper.ReadAssignmentWeight(o);
    } catch (...) {
        return -1.0;  // Not a valid ReadAssignmentWeight result (always in [0, 1]);
                       // an unambiguous "the call threw" sentinel for this f64-returning entry.
    }
}
