#pragma once
// ============================================================================
// LexiPair-BioRange | include/lexipair/lsm_index.h
// ----------------------------------------------------------------------------
// LSM-style dynamic update protocol (Slides 14, 16 Algorithm DYNAMICUPDATE):
//   - Index = {Base (large, sorted)} + {S1..Sk active segments (small, sorted)}
//   - Insert new haplotype H_delta: encode delta_N records, sort ONLY the new
//     batch O(deltaN log deltaN) independent of existing size, append as new
//     segment S_{k+1}.
//   - Range query merges across segments (cost per segment O(logU+logm+k)).
//   - Compaction: when |S|>THRESHOLD, background merge into single sorted
//     array O(Nmerged log Nmerged) off query-critical path.
//
// This beats full rebuild O(N log N) / O(N) on dense panels (E.coli 95.3%
// faster at H=64, H.sapiens Chr21 up to 30.5x), with honest overhead on
// low-depth huge genomes (O.sativa -52.7%, D.mela -10.0%) — Slide 19/22.
// ============================================================================
#include <cstdint>
#include <string>
#include <vector>
#include "lexipair/lexipair_index.h"

namespace lexipair {

class LsmIndex {
public:
    explicit LsmIndex(int K = 16, int compaction_threshold = 8)
        : K_(K), threshold_(compaction_threshold) {}

    // Initial build becomes Base segment.
    void build_base(const std::vector<std::string>& haplotypes, int first_hap_id = 0);
    // Insert one new haplotype sequence with given hap id -> new segment.
    // Returns ms spent (measured with high_resolution_clock by caller if needed).
    void insert_haplotype(const std::string& seq, int hap_id);
    // Query across all segments, merged.
    std::vector<KmerEntry> range_query(uint32_t code, int L, int R) const;
    // Force compaction now (normally background when size>threshold).
    void compact();

    size_t num_segments() const { return segments_.size(); }
    size_t total_records() const;
    bool needs_compaction() const { return (int)segments_.size() > threshold_; }

private:
    int K_;
    int threshold_;
    std::vector<LexiPairIndex> segments_;  // segments_[0] = base
};

}  // namespace lexipair
