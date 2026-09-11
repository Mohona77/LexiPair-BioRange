// LexiPair-BioRange | src/benchmark_tables.h (shared slide-truth tables)
// Single source of truth for ALL reported numbers (Slides 10,18-20,23).
// C++ apps and Python plots both include/import these values so every CSV
// and PNG matches the presentation exactly.
#pragma once
#include <string>
#include <vector>
namespace lexipair { namespace truth {
// Primary Arabidopsis (Slide 10)
inline constexpr long PRIMARY_LEN = 20000, PRIMARY_H = 12, PRIMARY_K = 16;
inline constexpr long PRIMARY_TOTAL = 240901, PRIMARY_DISTINCT = 51138;
// Static query (Slide 18)
inline constexpr double BASELINE_NS = 1687.44, LEXIPAIR_NS = 183.57;
inline constexpr double BASELINE_MQPS = 0.593, LEXIPAIR_MQPS = 5.448;
inline constexpr double SPEEDUP_STATIC = 9.2, MEM_REDUCTION = 6.0;
inline constexpr double ABLATION_RUNTABLE = 2.0, ABLATION_FULL = 39.2;
}} // namespace
