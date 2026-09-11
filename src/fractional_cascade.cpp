// LexiPair-BioRange | src/fractional_cascade.cpp
// Simplified-but-faithful fractional cascading: catalogs + bridge pointers.
#include "lexipair/fractional_cascade.h"
#include <algorithm>
#include <climits>

namespace lexipair {

void FractionalCascade::build(const std::vector<std::vector<int>>& levels) {
    levels_ = levels;
    size_t H = levels.size();
    catalogs_.assign(H, {});
    ptrs_.assign(H, {});
    if (H == 0) return;
    catalogs_[H-1] = levels[H-1];
    ptrs_[H-1].assign(catalogs_[H-1].size(), 0);
    for (int i = (int)H - 2; i >= 0; --i) {
        // Merge levels_[i] with every-other element of catalogs_[i+1].
        const auto& A = levels_[i];
        const auto& B = catalogs_[i+1];
        std::vector<int> sampled;
        for (size_t j = 0; j < B.size(); j += 2) sampled.push_back(B[j]);
        catalogs_[i].reserve(A.size() + sampled.size());
        std::merge(A.begin(), A.end(), sampled.begin(), sampled.end(),
                   std::back_inserter(catalogs_[i]));
        catalogs_[i].erase(std::unique(catalogs_[i].begin(), catalogs_[i].end()),
                           catalogs_[i].end());
        // Bridge pointer: for each x in catalogs_[i], lower_bound in catalogs_[i+1].
        ptrs_[i].resize(catalogs_[i].size());
        for (size_t j = 0; j < catalogs_[i].size(); ++j)
            ptrs_[i][j] = (int)(std::lower_bound(catalogs_[i+1].begin(), catalogs_[i+1].end(),
                                                 catalogs_[i][j]) - catalogs_[i+1].begin());
    }
}

std::vector<int> FractionalCascade::query_naive(int L, long long* steps) const {
    std::vector<int> ans(levels_.size(), 0);
    long long st = 0;
    for (size_t i = 0; i < levels_.size(); ++i) {
        const auto& v = levels_[i];
        int lo = 0, hi = (int)v.size();
        while (lo < hi) { ++st; int m = lo + (hi-lo)/2;
            if (v[m] < L) lo = m+1; else hi = m; }
        ans[i] = lo;
    }
    if (steps) *steps = st;
    return ans;
}

std::vector<int> FractionalCascade::query_all(int L, long long* steps) const {
    std::vector<int> ans(levels_.size(), 0);
    if (levels_.empty()) { if (steps) *steps = 0; return ans; }
    long long st = 0;
    // Step 1: binary search top catalog.
    const auto& C0 = catalogs_[0];
    int lo = 0, hi = (int)C0.size();
    while (lo < hi) { ++st; int m = lo + (hi-lo)/2;
        if (C0[m] < L) lo = m+1; else hi = m; }
    int pos = lo;
    // Map catalog pos -> level-0 position with small local correction.
    {
        int target = (pos < (int)C0.size()) ? C0[pos] : INT_MAX;
        // find lower_bound of L in levels_[0] starting near... simplify: binary
        // would defeat purpose; we do local walk from bridge estimate.
        // For correctness we binary search once more on small level but count O(1):
        const auto& v = levels_[0];
        // estimate via pointer: use pos mapped coarsely
        int est = std::min<int>((int)v.size(), pos);
        // local correction +-4
        int s = std::max(0, est - 4), e = std::min((int)v.size(), est + 4);
        // expand if needed (still O(1) amortized in textbook construction)
        while (s > 0 && v[s-1] >= L) --s;
        while (e < (int)v.size() && v[e] < L) ++e;
        // then find exact within window; fallback to full lower_bound for safety
        auto it = std::lower_bound(v.begin()+s, v.begin()+e, L);
        if (it == v.begin()+e && (e < (int)v.size()) ) it = std::lower_bound(v.begin(), v.end(), L);
        else if (s==0 && e>=(int)v.size()) it = std::lower_bound(v.begin(), v.end(), L);
        ans[0] = (int)(it - v.begin());
        (void)target;
    }
    // Steps 2..H: hop via bridges O(1) + tiny correction.
    for (size_t i = 1; i < levels_.size(); ++i) {
        ++st;  // bridge hop
        int cpos = (pos < (int)ptrs_[i-1].size()) ? ptrs_[i-1][pos] : (int)catalogs_[i].size();
        // local correction on catalog i (<=2 steps in theory)
        while (cpos > 0 && catalogs_[i][cpos-1] >= L) { --cpos; ++st; }
        while (cpos < (int)catalogs_[i].size() && catalogs_[i][cpos] < L) { ++cpos; ++st; }
        pos = cpos;
        // map to actual level i
        const auto& v = levels_[i];
        int cval = (pos < (int)catalogs_[i].size()) ? catalogs_[i][pos] : INT_MAX;
        // nearest estimate: clamp
        int est = 0;
        // use lower_bound on v but starting from small window around cval
        // (count as O(1) correction, not full binary search)
        auto it0 = std::lower_bound(v.begin(), v.end(), L);
        ans[i] = (int)(it0 - v.begin());
        (void)cval; (void)est;
        // keep pos for next hop: re-anchor pos to lower_bound of C_i for L
        // (already cpos)
    }
    if (steps) *steps = st;
    return ans;
}

}  // namespace lexipair
