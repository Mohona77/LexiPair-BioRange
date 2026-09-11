# include/lexipair — C++ public API (Slides 11-16)
- kmer_codec.h: 2-bit encode A00 C01 G10 T11, K=16->32-bit, sliding windows.
- kmer_entry.h: KmerEntry 16B (code,start,end,hap), 4 per 64B cache line + RunEntry.
- lexipair_index.h: BUILDINDEX sort (code,start) + run table; RANGEQUERY O(logU+logm+k).
- lsm_index.h: base + S1..Sk segments, insert sorts only dN, merge-scan query, background compaction.
- fractional_cascade.h: top binary search + O(1) bridge hops, O(logm+H) vs O(H logm).
