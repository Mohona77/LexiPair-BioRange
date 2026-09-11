#pragma once
// ============================================================================
// LexiPair-BioRange | include/lexipair/kmer_codec.h
// ----------------------------------------------------------------------------
// 2-bit DNA codec for k=16 -> 32-bit integer (Slide 11, Step 3).
//   A=00  C=01  G=10  T=11
// Lexicographic order of strings == numeric order of codes (2 bits/base,
// big-endian packing). This is what makes "lexicographic sort" a plain
// integer sort and enables binary search on the run table.
//
// Descriptive notes for GitHub readers:
//  - encode_kmer(s): s.length()==K (default 16), throws on bad base.
//  - decode_kmer(code,K): inverse, used for validation / random sampling.
//  - encode_sliding(seq): all N-K+1 codes with start/end positions.
// ============================================================================
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace lexipair {

inline constexpr int K_DEFAULT = 16;

inline uint8_t base_to_bits(char b) {
    switch (b) {
        case 'A': case 'a': return 0;  // 00
        case 'C': case 'c': return 1;  // 01
        case 'G': case 'g': return 2;  // 10
        case 'T': case 't': return 3;  // 11
        default: throw std::invalid_argument("bad DNA base");
    }
}

inline char bits_to_base(uint8_t v) {
    static const char m[4] = {'A','C','G','T'};
    return m[v & 3];
}

// Encode exactly K bases (big-endian: first base = most significant bits).
inline uint32_t encode_kmer(const std::string& s, int K = K_DEFAULT) {
    if ((int)s.size() != K) throw std::invalid_argument("kmer length != K");
    uint32_t code = 0;
    for (char b : s) code = (code << 2) | base_to_bits(b);
    return code;
}

inline std::string decode_kmer(uint32_t code, int K = K_DEFAULT) {
    std::string s(K, 'A');
    for (int i = K - 1; i >= 0; --i) { s[i] = bits_to_base(code & 3); code >>= 2; }
    return s;
}

struct CodedKmer { uint32_t code; int start; int end; };  // end = start+K-1

// Slide-window scan: N-K+1 codes, step=1.
inline std::vector<CodedKmer> encode_sliding(const std::string& seq, int K = K_DEFAULT) {
    std::vector<CodedKmer> out;
    if ((int)seq.size() < K) return out;
    // rolling encode O(N)
    uint32_t code = 0;
    for (int i = 0; i < K; ++i) code = (code << 2) | base_to_bits(seq[i]);
    out.push_back({code, 0, K - 1});
    uint32_t mask = (K == 16) ? 0xFFFFFFFFu : ((K*2 >= 32) ? 0xFFFFFFFFu : ((1u << (2*K)) - 1));
    for (int i = K; i < (int)seq.size(); ++i) {
        code = ((code << 2) | base_to_bits(seq[i])) & mask;
        out.push_back({code, i - K + 1, i});
    }
    return out;
}

}  // namespace lexipair
