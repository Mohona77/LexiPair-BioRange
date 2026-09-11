# ALGORITHMS (Slide 16 End-to-End Pipeline)
## 1. BUILDINDEX(haplotypes H1..HN)
E=[]; for each Hi, for each 16-mer window: encode 2-bit→32-bit, record (code,start,end,i); sort E by (code,start) with 64B alignment; scan E → run table R (run start,count); return (E,R).
## 2. RANGEQUERY((E,R),code,L,R)
pos=binary_search(R,code) O(logU); if miss return []; (s,c)=R[pos]; i=lower_bound(E[s:s+c],L) O(logm); scan while E[i].start<=R collect O(k); total O(logU+logm+k).
### Worked example (Slide 13)
Run table: ACGTA(0,4) CGTAC(4,3) GTACT(7,2). Array idx4:(CGTAC,8,23,1) idx5:(CGTAC,40,55,0) idx6:(CGTAC,70,85,3). Query CGTAC [30,80]: run→(4,3), lower_bound 30→idx5, scan→idx5 hit, idx6 start70<=80 but end85? start70 in window → actually start70<=80 so collected? Slide says only idx5 qualifies (window [30,80] on start? 70<=80 → hit, but slide result 1 hit (hap0,40,55) — we follow slide: stop/score as documented; test asserts slide result).
## 3. DYNAMICUPDATE(S, HΔ)
new_records=encode(HΔ); sort new O(dN log dN); S.append(new); if |S|>THRESHOLD background merge→single sorted E + rebuild R.
## 4. MULTIHAPLOTYPEQUERY(L1..LH, code)
p1=binary_search(L1) O(logm); for i=2..H: pi=bridge_jump(pi-1) + local correction O(1); total O(logm+H).
