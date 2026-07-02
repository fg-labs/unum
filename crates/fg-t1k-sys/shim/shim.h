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

// Opaque-handle FFI for KmerIndex: lets Rust drive a real C++ KmerIndex
// instance (Insert/Search/Remove/BuildIndexFromRead) for lockstep
// differential testing. The handle is an opaque `void*` (really a
// `KmerIndex*`); callers must free it exactly once via
// fg_t1k_kmerindex_free. `fg_t1k_kmerindex_new` returns NULL if construction
// throws (KmerIndex's constructor does `new std::map<...>[1000003]`, a large
// allocation that can throw std::bad_alloc) -- callers MUST check for a NULL
// handle before using it.
//
// `idx`/`offset` are `uint32_t` to match T1K's `index_t` (`defs.h:15`:
// `typedef uint32_t index_t`), the declared parameter/field type of
// `_indexInfo` and `Insert`/`Remove`/`BuildIndexFromRead`'s `idx` and
// `offset` parameters. `strand` stays `int`, matching `Insert`/`Remove`'s
// `int strand` parameter -- it is accepted but never stored (the `strand`
// field of `_indexInfo` is commented out in the vendored header).
//
// The KmerCode arguments are the same opaque handles produced by
// fg_t1k_kmercode_new (see above), so a differential test can drive one
// KmerCode instance and feed it to both the Rust and C++ KmerIndex sides in
// lockstep, guaranteeing identical GetCode()/IsValid() inputs.
void* fg_t1k_kmerindex_new(void);
void fg_t1k_kmerindex_free(void* p);
void fg_t1k_kmerindex_insert(void* idxp, void* kcp, uint32_t idx, uint32_t offset, int strand);
void fg_t1k_kmerindex_remove(void* idxp, void* kcp, uint32_t idx, uint32_t offset, int strand);
int fg_t1k_kmerindex_search_size(void* idxp, void* kcp);
// Fills `*out_idx`/`*out_offset` with the `i`-th entry (0-based) of the most
// recent Search result. `i` must be in `0..fg_t1k_kmerindex_search_size(...)`.
void fg_t1k_kmerindex_search_get(void* idxp, void* kcp, int i, uint32_t* out_idx,
                                  uint32_t* out_offset);
// `s` need not be NUL-terminated; exactly `len` bytes are read (matches
// `BuildIndexFromRead(KmerCode&, char* s, int len, int id, int shift = 0)`,
// which never calls strlen on `s`).
void fg_t1k_kmerindex_build_index_from_read(void* idxp, void* kcp, const char* s, int len,
                                              int id, int shift);

// Opaque-handle FFI for SeqSet (the FastqExtractor/BamExtractor read-
// candidate-filtering slice: reference load + GetHitsFromRead + HasHitInSet +
// IsLowComplexity/IsGoodCandidate). These entries call the real, unmodified
// stock C++ (e.g. HasHitInSet runs IN FULL, including its
// GetOverlapsFromHits/AlignAlgo mismatch-threshold confirmation) so they serve
// as the differential oracle for fg-t1k-core's Rust port. The handle is an
// opaque `void*` (really a `SeqSet*`);
// callers must free it exactly once via fg_t1k_seqset_free.
// `fg_t1k_seqset_new` returns NULL if construction throws -- callers MUST
// check for a NULL handle before using it.
void* fg_t1k_seqset_new(int kmerLength);
void fg_t1k_seqset_free(void* p);
// Mirrors `SeqSet::InputRefFa(char* filename)` (SeqSet.hpp:872-904): loads
// every FASTA record in `fastaPath`, indexing each sequence's k-mers with
// `id` = its 0-based load order (matching stock's `id = seqs.size()` at
// insertion time). Returns 0 on success, -1 if the underlying call threw
// (e.g. the file could not be opened).
int fg_t1k_seqset_load_ref(void* p, const char* fastaPath);
// Mirrors the free function `IsLowComplexity` (FastqExtractor.cpp:89-111,
// byte-identical to BamExtractor.cpp:144-166) -- NOT a SeqSet method; `p` is
// unused (kept for call-site symmetry with the other fg_t1k_seqset_* FFI
// entries). `seq` must be a NUL-terminated C string. Returns 1/0.
int fg_t1k_seqset_is_low_complexity(void* p, const char* seq);
// Mirrors `SeqSet::HasHitInSet(char* read, char* rcRead)` (SeqSet.hpp:
// 1915-1990) IN FULL, including the GetOverlapsFromHits/AlignAlgo-based
// mismatch-threshold confirmation this port's Rust side deliberately does
// NOT implement (see ref_kmer_filter module docs) -- this is the real,
// unmodified stock C++ oracle. `read` must be a NUL-terminated C string;
// `rcRead`'s scratch buffer is allocated internally. Returns 1/0.
int fg_t1k_seqset_has_hit_in_set(void* p, const char* read);
// Mirrors the free function `IsGoodCandidate` (FastqExtractor.cpp:113-118):
// `!IsLowComplexity(read) && refSet->HasHitInSet(read, buffer)`. Returns
// 1/0.
int fg_t1k_seqset_is_good_candidate(void* p, const char* read);

// Opaque-handle FFI for Alignments (vendor/t1k/alignments.hpp), scoped to
// exactly the slice BamExtractor.cpp uses. The handle is an opaque `void*`
// (really an `Alignments*`); callers must free it exactly once via
// fg_t1k_alignments_free. Alignments' constructor is lightweight (no heap
// allocation that can throw), but Open/Next/GetGeneralInfo can still throw
// via a corrupt/unreadable file (or, on the real vendored source, `exit(1)`
// on a hard open failure -- that path cannot be caught by try/catch at all,
// so callers should only ever point fg_t1k_alignments_open at a file they
// already know exists and parses); every entry point here is still wrapped
// in try/catch per this shim's exception-safety rule (no C++ exception may
// unwind across an `extern "C"` boundary) for defense-in-depth.
void* fg_t1k_alignments_new(void);
void fg_t1k_alignments_free(void* p);
// Returns 0 on success, -1 if the underlying call threw. `path` must be a
// NUL-terminated C string.
int fg_t1k_alignments_open(void* p, const char* path);
// Returns 0 on success, -1 if the underlying call threw.
int fg_t1k_alignments_rewind(void* p);
// Mirrors `Alignments::Next`. Returns 1 if a record was read, 0 at EOF, -1
// if the underlying call threw.
int fg_t1k_alignments_next(void* p);
// Mirrors `Alignments::GetReadSeq`/`GetQual`: writes a NUL-terminated string
// into `buffer` (caller-owned, must be at least `record.l_qseq + 1` bytes --
// the shim itself does not know l_qseq, so callers should size generously,
// e.g. the same 100001-byte convention used elsewhere in this shim). Returns
// 0 on success, -1 if the underlying call threw.
int fg_t1k_alignments_get_read_seq(void* p, char* buffer);
int fg_t1k_alignments_get_qual(void* p, char* buffer);
// Mirrors `Alignments::GetReadId`. Writes a NUL-terminated string into
// `buffer` (caller-owned, must be large enough for any QNAME the test BAM
// contains). Returns 0 on success, -1 if the underlying call threw.
int fg_t1k_alignments_get_read_id(void* p, char* buffer);
// Flag predicates. Return 1/0, or -1 if the underlying call threw.
int fg_t1k_alignments_is_first_mate(void* p);
int fg_t1k_alignments_is_reverse(void* p);
int fg_t1k_alignments_is_mate_reverse(void* p);
int fg_t1k_alignments_is_aligned(void* p);
int fg_t1k_alignments_is_template_aligned(void* p);
int fg_t1k_alignments_is_primary(void* p);
// Mirrors `Alignments::GetChromId`. Returns the tid, or INT32_MIN if the
// underlying call threw (a tid is never that value in practice, so this
// remains an unambiguous sentinel without needing a separate out-param).
int fg_t1k_alignments_get_chrom_id(void* p);
// Mirrors `Alignments::segCnt`. Returns the count, or -1 if the underlying
// call threw.
int fg_t1k_alignments_seg_count(void* p);
// Mirrors `Alignments::segments[i]`. Writes `.a`/`.b` into the out-params.
// `i` must be in `0..fg_t1k_alignments_seg_count(p)`. Returns 0 on success,
// -1 if the underlying call threw.
int fg_t1k_alignments_seg(void* p, int i, int64_t* out_a, int64_t* out_b);
// Mirrors `Alignments::GetGeneralInfo(stopEarly)` then reads the resulting
// `fragStdev`/`readLen` public fields into the out-params. Returns 0 on
// success, -1 if the underlying call threw.
int fg_t1k_alignments_general_info(void* p, int stop_early, int* out_frag_stdev, int* out_read_len);
#ifdef __cplusplus
}
#endif
