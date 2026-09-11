#pragma once
// ============================================================================
// LexiPair-BioRange | include/lexipair/lexipair_index.h
// ----------------------------------------------------------------------------
// Static cache-conscious index (Slides 12-13, 16 Algorithm BUILDINDEX/RANGEQUERY):
//
//   BUILDINDEX(H_1..H_N):
//     E <- [] ; for each hap i, for each 16-mer window: encode 2-bit code,
//          record (code,start,end,i); sort E by (code,start) with 64-byte
//          aligned alloc; build run table R by scanning E (run start,count).
//     return (E,R).  Complexity O(N log N), N = total occurrences (240,901
//     for primary Arabidopsis dataset, U=51,138 distinct).
//
//   RANGEQUERY((E,R), query_code, L, R):
//     1. binary search run table -> (s,c)                        O(log U)
//     2. lower_bound start>=L inside E[s..s+c-1]                  O(log m)
//     3. linear scan while start<=R, collect hits                 O(k)
//     Total O(log U + log m + k) vs naive O(m) dictionary scan.
//
// Validated against Python-dict baseline (Slide 11 Step 7).
// ============================================================================
#include <cstdint>
#include <string>
#include <vector>
#include "lexipair/kmer_entry.h"

namespace lexipair {

class LexiPairIndex {
public:
    std::vector<KmerEntry> entries;  // E, sorted by (code,start)
    std::vector<RunEntry>  runs;     // R, sorted by code

    LexiPairIndex() = default;

    // Build from haplotype strings. K=16, step=1. Deterministic.
    void build(const std::vector<std::string>& haplotypes, int K = 16);

    // Number of distinct k-mers U, total records N.
    size_t num_records() const { return entries.size(); }
    size_t num_distinct() const { return runs.size(); }

    // Static range query. Returns matching entries (copies).
    std::vector<KmerEntry> range_query(uint32_t query_code, int L, int R) const;

    // Helpers (public for testing / teaching):
    // binary search run table; -1 if absent.
    int find_run(uint32_t code) const;
    // lower bound of start>=L inside run [s, s+c).
    int lower_bound_in_run(int s, int c, int L) const;

    // Memory estimate: N*16 + U*12 bytes (structural).
    size_t structural_bytes() const {
        return entries.size()*sizeof(KmerEntry) + runs.size()*sizeof(RunEntry);
    }
    // Cache lines touched for k hits (worst): ceil(k/4)+1.
    static int cache_lines_for_hits(int k) {
        if (k <= 0) return 1;
        return (k + (int)ENTRIES_PER_CACHE_LINE - 1) / (int)ENTRIES_PER_CACHE_LINE + 1;
    }
};

}  // namespace lexipair
