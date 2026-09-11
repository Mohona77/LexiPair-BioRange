// LexiPair-BioRange | src/lsm_index.cpp
#include "lexipair/lsm_index.h"
#include "lexipair/kmer_codec.h"
#include <algorithm>

namespace lexipair {

void LsmIndex::build_base(const std::vector<std::string>& haps, int first_hap_id) {
    segments_.clear();
    LexiPairIndex base;
    // Temporarily build with local hap ids then shift to global ids.
    base.build(haps, K_);
    if (first_hap_id != 0)
        for (auto &e : base.entries) e.hap += first_hap_id;
    segments_.push_back(std::move(base));
}

void LsmIndex::insert_haplotype(const std::string& seq, int hap_id) {
    // Encode only delta_N, sort only new batch.
    LexiPairIndex seg;
    seg.build(std::vector<std::string>{seq}, K_);
    for (auto &e : seg.entries) e.hap = hap_id;  // build() sets hap=0; fix
    segments_.push_back(std::move(seg));
    // NOTE: compaction is explicit via compact() / needs_compaction() so
    // query-critical path stays fast (background merge in production).
}

std::vector<KmerEntry> LsmIndex::range_query(uint32_t code, int L, int R) const {
    std::vector<KmerEntry> out;
    for (auto &seg : segments_) {
        auto part = seg.range_query(code, L, R);
        out.insert(out.end(), part.begin(), part.end());
    }
    std::sort(out.begin(), out.end(), [](const KmerEntry& a, const KmerEntry& b){
        if (a.start != b.start) return a.start < b.start;
        return a.hap < b.hap;
    });
    return out;
}

void LsmIndex::compact() {
    if (segments_.size() <= 1) return;
    std::vector<KmerEntry> all;
    for (auto &seg : segments_)
        all.insert(all.end(), seg.entries.begin(), seg.entries.end());
    std::sort(all.begin(), all.end(), [](const KmerEntry& a, const KmerEntry& b){
        if (a.code != b.code) return a.code < b.code;
        if (a.start != b.start) return a.start < b.start;
        return a.hap < b.hap;
    });
    LexiPairIndex merged;
    merged.entries = std::move(all);
    merged.runs.clear();
    for (size_t i = 0; i < merged.entries.size();) {
        size_t j = i + 1;
        while (j < merged.entries.size() && merged.entries[j].code == merged.entries[i].code) ++j;
        merged.runs.push_back({merged.entries[i].code, (int32_t)i, (int32_t)(j - i)});
        i = j;
    }
    segments_.clear();
    segments_.push_back(std::move(merged));
}

size_t LsmIndex::total_records() const {
    size_t n = 0;
    for (auto &s : segments_) n += s.num_records();
    return n;
}

}  // namespace lexipair
