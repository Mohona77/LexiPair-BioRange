# COMPLEXITY (Slide 21)
| Category | Naive | Best | Worst | Overall |
|---|---|---|---|---|
| Memory | O(N·100+)B | O(N·S)B | O(N·S+U·R)B | O(N·S+U·R)B |
| Cache traffic (k hits) | O(k) misses | O(1) lines | ceil(k/(C/S))+1 | ceil(k/(C/S))+1 |
| Static query | O(m) | O(1) | O(logU+logm+k) | O(logU+logm+k) |
| Dynamic ingestion | O(N) rebuild | O(dN) | O(dN log dN) | O(dN log dN) |
| Multi-hap search | O(H log m) | O(1+H) | O(logm+H) | O(logm+H) |
Key: cache-friendly + update-friendly + scales like a pro.
