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
