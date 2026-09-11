// LexiPair-BioRange | src/lexipair_index.cpp
// Implements BUILDINDEX + RANGEQUERY (Slide 16, Algorithm 1 & 2).
#include "lexipair/lexipair_index.h"
#include "lexipair/kmer_codec.h"
#include <algorithm>

namespace lexipair {

void LexiPairIndex::build(const std::vector<std::string>& haplotypes, int K) {
    entries.clear(); runs.clear();
    for (size_t h = 0; h < haplotypes.size(); ++h) {
        const auto& seq = haplotypes[h];
        auto coded = encode_sliding(seq, K);
        for (auto &ck : coded)
            entries.push_back({ck.code, ck.start, ck.end, (int32_t)h});
    }
    std::sort(entries.begin(), entries.end(), [](const KmerEntry& a, const KmerEntry& b){
        if (a.code != b.code) return a.code < b.code;
        if (a.start != b.start) return a.start < b.start;
        return a.hap < b.hap;
    });
    // Run table by scanning sorted array.
    for (size_t i = 0; i < entries.size();) {
        size_t j = i + 1;
        while (j < entries.size() && entries[j].code == entries[i].code) ++j;
        runs.push_back({entries[i].code, (int32_t)i, (int32_t)(j - i)});
        i = j;
    }
}

int LexiPairIndex::find_run(uint32_t code) const {
    int lo = 0, hi = (int)runs.size() - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        if (runs[mid].code == code) return mid;
        else if (runs[mid].code < code) lo = mid + 1;
        else hi = mid - 1;
    }
    return -1;
}

int LexiPairIndex::lower_bound_in_run(int s, int c, int L) const {
    int lo = s, hi = s + c;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (entries[mid].start < L) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

std::vector<KmerEntry> LexiPairIndex::range_query(uint32_t q, int L, int R) const {
    std::vector<KmerEntry> out;
    int pos = find_run(q);
    if (pos < 0) return out;
    int s = runs[pos].run_start, c = runs[pos].run_count;
    int i = lower_bound_in_run(s, c, L);
    while (i < s + c && entries[i].start <= R) { out.push_back(entries[i]); ++i; }
    return out;
}

}  // namespace lexipair
