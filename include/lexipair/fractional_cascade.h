#pragma once
// ============================================================================
// LexiPair-BioRange | include/lexipair/fractional_cascade.h
// ----------------------------------------------------------------------------
// Fractional-cascading multi-haplotype search (Slides 15-16):
//   Naive: H independent binary searches -> O(H log m).
//   Ours:  one binary search on top level + O(1) bridge hops per level
//          + tiny local correction -> O(log m + H).
//
// Teaching model (matches slide diagram):
//   - Level 0 (top): sampled catalog.
//   - Levels 1..H: per-haplotype sorted start arrays for the queried k-mer's
//     run (in full system: per-haplotype sorted positions).
//   - Bridges: every other element of level i+1 promoted into level i with a
//     pointer; build bottom-up in O(sum |level_i|).
//   - Query: binary search top, then follow bridge pointers down, correcting
//     by <=2 steps locally.
//
// The class below is a clean, testable implementation of that idea operating
// on integer position arrays. LexiPairIndex::range_query handles the single-
// k-mer run case; this class accelerates the H-level fan-out.
// ============================================================================
#include <vector>

namespace lexipair {

class FractionalCascade {
public:
    // Build from H sorted levels (each sorted ascending). Copies data and
    // builds bridge catalogs bottom-up.
    void build(const std::vector<std::vector<int>>& levels);
    // Query value L (lower bound) on every level. Returns per-level lower-
    // bound indices. Cost O(log m + H). Also returns total steps for teaching.
    std::vector<int> query_all(int L, long long* steps_out = nullptr) const;
    // Naive baseline for comparison: H binary searches.
    std::vector<int> query_naive(int L, long long* steps_out = nullptr) const;

    size_t num_levels() const { return levels_.size(); }

private:
    std::vector<std::vector<int>> levels_;   // original sorted levels
    // bridges_[i] = sampled catalog for level i (level i + every-other of i+1...)
    // For simplicity + fidelity we store full merged catalog with pointers:
    std::vector<std::vector<int>> catalogs_;
    std::vector<std::vector<int>> ptrs_;     // ptrs_[i][j] -> index in catalogs_[i+1]
};

}  // namespace lexipair
