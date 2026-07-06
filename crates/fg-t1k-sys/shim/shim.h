#include <stddef.h>
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
// Mirrors `SeqSet::GetOverlapsFromRead(char* read, int strand, int barcode,
// std::vector<_overlap>&)` (SeqSet.hpp:1594-1912) IN FULL -- the real,
// unmodified stock C++ read-to-allele alignment/scoring core (Task 4b's
// differential oracle). Always called with `strand=0, barcode=-1` (matching
// fg_t1k_core::ref_kmer_filter::RefKmerFilter::get_overlaps_from_read's own
// fixed parameters -- see that method's doc comment). Writes up to
// `out_capacity` overlaps into the caller-allocated `out_*` arrays (count-
// then-get is unnecessary here since GetOverlapsFromRead's result size is
// already known after one call); `*out_count` receives the actual number of
// overlaps returned by the real C++ call (which may exceed `out_capacity`,
// in which case only the first `out_capacity` are written -- callers should
// size `out_capacity` generously, e.g. 10000, and treat a `*out_count >
// out_capacity` return as a test-fixture sizing bug, not a real result).
// `read` must be a NUL-terminated C string. Returns 0 on success (including
// the legitimate `read.len() < kmerLength` case, which mirrors stock's `-1`
// return by setting `*out_count = -1`), -1 if the underlying call threw.
int fg_t1k_seqset_get_overlaps_from_read(void* p, const char* read, int out_capacity,
                                          int* out_seq_idx, int* out_strand,
                                          int* out_read_start, int* out_read_end,
                                          int* out_seq_start, int* out_seq_end,
                                          int* out_match_cnt, double* out_similarity,
                                          int* out_count);

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
// into `buffer` (caller-owned, `buffer_size` bytes). The record needs
// `l_qseq + 1` bytes; if it would not fit, the shim returns -1 instead of
// overflowing. Returns 0 on success, -1 on overflow or if the underlying call
// threw.
int fg_t1k_alignments_get_read_seq(void* p, char* buffer, size_t buffer_size);
int fg_t1k_alignments_get_qual(void* p, char* buffer, size_t buffer_size);
// Mirrors `Alignments::GetReadId`. Writes the NUL-terminated QNAME into
// `buffer` (caller-owned, `buffer_size` bytes); if the id would not fit, the
// shim returns -1 instead of overflowing. Returns 0 on success, -1 on overflow
// or if the underlying call threw.
int fg_t1k_alignments_get_read_id(void* p, char* buffer, size_t buffer_size);
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

// Free-function FFI for AlignAlgo (vendor/t1k/AlignAlgo.hpp), scoped to
// exactly the slice the genotyping path calls (see
// crates/fg-t1k-core/src/align_algo.rs module docs for the full
// used-vs-unused method inventory): GlobalAlignment,
// GlobalAlignment_PosWeight, plus SeqSet::GetAlignStats (a SeqSet member
// function that touches no instance state, so it is exposed here as a free
// function too). AlignAlgo is header-only and pure C++ (no htslib/samtools),
// so these entries call the real static AlignAlgo:: methods directly -- no
// opaque handle needed, unlike SeqSet/KmerIndex above.
//
// `out_align` is a caller-allocated buffer; both AlignAlgo methods can emit
// at most `lent + lenp` ops (every DP step consumes at least one base from
// t or p), so callers must size `out_align` to at least `lent + lenp` bytes
// (matching the 100001-byte convention used elsewhere in this shim is more
// than sufficient for any realistic test size). `*out_len` is set to the
// number of ops written (NOT including the C++ side's internal `-1`
// sentinel, which this shim strips before returning -- Rust callers get the
// op count directly, no sentinel scanning required). Returns the alignment
// score, or INT32_MIN if the underlying call threw.
int fg_t1k_alignalgo_global_alignment(const char* t, int lent, const char* p, int lenp, int band,
                                       signed char* out_align, int* out_len);
// `t_weights` is a flat array of `4 * lent` ints (four counts per reference
// position, matching `struct _posWeight { int count[4]; }` laid out
// contiguously); same `out_align`/`out_len` sizing/semantics as
// fg_t1k_alignalgo_global_alignment above.
int fg_t1k_alignalgo_global_alignment_pos_weight(const int* t_weights, int lent, const char* p,
                                                  int lenp, signed char* out_align, int* out_len);
// Mirrors SeqSet::GetAlignStats: `align` must be a NUL-terminated-by-`-1`
// signed-char op sequence (i.e. exactly what
// fg_t1k_alignalgo_global_alignment*'s C++ side produces internally before
// this shim strips the sentinel -- this entry re-appends a `-1` sentinel
// itself from `align`/`align_len` so it can call the real, unmodified
// GetAlignStats unchanged). `update` mirrors the C++ `bool update` parameter
// (0/1). Writes match/mismatch/indel counts into the out-params. Returns 0
// on success, -1 if the underlying call threw.
int fg_t1k_seqset_get_align_stats(const signed char* align, int align_len, int update,
                                   int* out_match_cnt, int* out_mismatch_cnt, int* out_indel_cnt);

// Free-function FFI for Genotyper's Task-5a deterministic slice (vendor/t1k/
// Genotyper.hpp), scoped to exactly ParseAlleleName and ReadAssignmentWeight
// -- see crates/fg-t1k-core/src/genotyper.rs module docs for the full port
// scope. Both are private C++ methods reached via a `#define private
// public` trick in shim.cpp (see that file for the full rationale); each
// call constructs a fresh, throwaway `Genotyper(11)` (matching the
// `Genotyper genotyper(11)` construction in Genotyper.cpp:207's real
// `main`), so these are free functions rather than opaque-handle FFI like
// fg_t1k_seqset_* above -- no persistent state needs to survive between
// calls.
//
// Mirrors `Genotyper::ParseAlleleName(char* allele, char* gene, char*
// majorAllele, int fieldsType)` (Genotyper.hpp:63-131) after calling
// `SetAlleleNameStructure(alleleDigitUnits, alleleDelimiter)`
// (Genotyper.hpp:548-552) to configure the same two fields the real
// `Genotyper.cpp:336` sets from CLI flags before ever calling
// ParseAlleleName. `allele` must be a NUL-terminated C string; `out_gene`/
// `out_major` are caller-allocated buffers, each must be at least
// `strlen(allele) + 1` bytes (ParseAlleleName never writes a byte beyond
// `strlen(allele)` into either -- it only NUL-truncates a copy of `allele`
// at an earlier offset). `alleleDelimiter` is `char` (NOT `unsigned char`)
// to exactly match the C++ `char alleleDelimiter` field's signedness on the
// build platform; pass `'\0'` for "not set" (`SetAlleleNameStructure`'s own
// "no override" sentinel). Returns 0 on success, -1 if the underlying call
// threw.
int fg_t1k_genotyper_parse_allele_name(const char* allele, int fieldsType, int alleleDigitUnits,
                                        char alleleDelimiter, char* out_gene, char* out_major);
// Mirrors `Genotyper::ReadAssignmentWeight(const _fragmentOverlap& o)`
// (Genotyper.hpp:205-230) after calling
// `genotyper.refSet.SetRefSeqSimilarity(refSeqSimilarity)`
// (SeqSet.hpp:835-838) so `o.similarity`'s comparison thresholds (which
// derive from `refSet.GetRefSeqSimilarity()`, Genotyper.hpp:212) use the
// caller's chosen value rather than `SeqSet`'s constructor default (`0.8`,
// SeqSet.hpp:768). Every other `_fragmentOverlap` field is zero-initialized
// (`ReadAssignmentWeight` reads only `.similarity`/`.hasN`, per this port's
// scope -- see genotyper.rs's `FragmentOverlap` doc comment for why the
// other fields are irrelevant here). Returns the real C++ weight (always in
// `[0, 1]`), or `-1.0` (an otherwise-impossible return value for this
// function) if the underlying call threw.
double fg_t1k_genotyper_read_assignment_weight(double similarity, int hasN,
                                                double refSeqSimilarity);

// Opaque-handle FFI for Genotyper's Task-5b quantification slice
// (CoalesceReadAssignments, BuildAlleleEquivalentClass via
// FinalizeReadAssignments, QuantifyAlleleEquivalentClass) -- see shim.cpp for
// the full rationale (a persistent handle, unlike the Task-5a free functions
// above, since a scripted differential test builds up state across many
// calls). The handle is an opaque `void*` (really a `Genotyper*`); callers
// must free it exactly once via fg_t1k_genotyper2_free.
// `fg_t1k_genotyper2_new` returns NULL if construction throws -- callers
// MUST check for a NULL handle before using it.
void* fg_t1k_genotyper2_new(int kmerLength);
void fg_t1k_genotyper2_free(void* p);
// Mirrors `genotyper.refSet.InputRefSeq(id, seq, weight,
// /*initExonInfo=*/true)`. `id`/`seq` must be NUL-terminated C strings.
// Returns the 0-based seqIdx (== alleleIdx) on success, -1 if the underlying
// call threw.
int fg_t1k_genotyper2_add_ref_seq(void* p, const char* id, const char* seq, int weight);
// Mirrors `Genotyper::InitAlleleInfo()`. Call after all ref seqs are loaded.
// Returns 0 on success, -1 if the underlying call threw.
int fg_t1k_genotyper2_init_allele_info(void* p);
// Mirrors `Genotyper::SetAlleleNameStructure(alleleDigitUnits,
// alleleDelimiter)`. Returns 0 on success, -1 if the underlying call threw.
int fg_t1k_genotyper2_set_allele_name_structure(void* p, int alleleDigitUnits,
                                                 char alleleDelimiter);
// Mirrors `Genotyper::InitReadAssignments(totalReadCnt, maxAssignCnt)`.
// Returns 0 on success, -1 if the underlying call threw.
int fg_t1k_genotyper2_init_read_assignments(void* p, int totalReadCnt, int maxAssignCnt);
// Mirrors `Genotyper::SetReadAssignments(readId, assignment)`, building the
// `std::vector<_fragmentOverlap>` from `n`-length flat parallel arrays (one
// entry per `_fragmentOverlap`); `refSeqSimilarity` is applied via
// `SetRefSeqSimilarity` before the call. Returns 0 on success, -1 if the
// underlying call threw.
int fg_t1k_genotyper2_set_read_assignments(void* p, int readId, const int* seqIdx,
                                            const int* seqStart, const int* seqEnd,
                                            const int* matchCnt, const int* relaxedMatchCnt,
                                            const double* similarity, const int* hasMatePair,
                                            const int* o1FromR2, const double* qual,
                                            const int* hasN, int n, double refSeqSimilarity);
// Mirrors `Genotyper::CoalesceReadAssignments(begin, end)`. Returns the C++
// method's own return value, or INT32_MIN if the underlying call threw.
int fg_t1k_genotyper2_coalesce_read_assignments(void* p, int begin, int end);
// Mirrors `Genotyper::FinalizeReadAssignments()` (also runs the real
// `GetSeqMissingBaseCoverage`/`BuildAlleleEquivalentClass`). Returns the C++
// method's own return value, or INT32_MIN if the underlying call threw.
int fg_t1k_genotyper2_finalize_read_assignments(void* p);
// Mirrors `Genotyper::QuantifyAlleleEquivalentClass()` -- the SQUAREM
// abundance EM, run to completion. Returns the number of EM iterations run,
// or INT32_MIN if the underlying call threw.
int fg_t1k_genotyper2_quantify(void* p);

// Read-back accessors (all real Genotyper state, no shim-side recomputation).
// The index-taking accessors bounds-check their index and return the sentinel
// noted below (rather than reading out of range, which would be UB): the
// `int`-returning ones return -1, the `double`-returning ones return -1.0.
// `ecIdx` is bounds `0..fg_t1k_genotyper2_ec_count(p)`, `memberIdx` is bounds
// `0..fg_t1k_genotyper2_ec_member_count(p, ecIdx)`, and `alleleIdx` is bounds
// the Genotyper's allele count (`refSet.Size()` for the `seq_*` entries).
int fg_t1k_genotyper2_read_cnt(void* p);
int fg_t1k_genotyper2_ec_count(void* p);
int fg_t1k_genotyper2_ec_member_count(void* p, int ecIdx);
int fg_t1k_genotyper2_ec_member(void* p, int ecIdx, int memberIdx);
int fg_t1k_genotyper2_allele_equivalent_class(void* p, int alleleIdx);
double fg_t1k_genotyper2_allele_abundance(void* p, int alleleIdx);
double fg_t1k_genotyper2_allele_ec_abundance(void* p, int alleleIdx);
int fg_t1k_genotyper2_allele_missing_coverage(void* p, int alleleIdx);
int fg_t1k_genotyper2_seq_effective_len(void* p, int alleleIdx);
int fg_t1k_genotyper2_seq_weight(void* p, int alleleIdx);
#ifdef __cplusplus
}
#endif
