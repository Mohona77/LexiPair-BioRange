#pragma once
// ============================================================================
// LexiPair-BioRange | include/lexipair/kmer_entry.h
// ----------------------------------------------------------------------------
// Fixed-size record (Slide 12, LexiPair Idea):
//    KmerEntry = (kmer code:4B, start:4B, end:4B, hap id:4B) = 16 bytes.
// Four records pack into one 64-byte CPU cache line -> spatial locality,
// sequential scans touch ceil(k/4)+1 lines instead of O(k) scattered misses.
//
// Sorting key: (kmer code ASC, start ASC). All occurrences of one k-mer form
// a contiguous "run" sorted by start, so range queries binary-search inside
// the run and stop early (Slide 13).
// ============================================================================
#include <cstdint>
#include <cstddef>

namespace lexipair {

struct KmerEntry {
    uint32_t code;   // 32-bit 2-bit encoded 16-mer
    int32_t  start;  // genomic start coordinate within haplotype
    int32_t  end;    // start + K - 1
    int32_t  hap;    // haplotype id
};

static_assert(sizeof(KmerEntry) == 16, "KmerEntry must be exactly 16 bytes");
static_assert(alignof(KmerEntry) <= 16, "unexpected alignment");

inline constexpr std::size_t CACHE_LINE_BYTES = 64;
inline constexpr std::size_t ENTRIES_PER_CACHE_LINE =
    CACHE_LINE_BYTES / sizeof(KmerEntry);  // == 4

struct RunEntry {
    uint32_t code;      // distinct k-mer code
    int32_t  run_start; // index into sorted KmerEntry array
    int32_t  run_count; // run length m
};

}  // namespace lexipair
