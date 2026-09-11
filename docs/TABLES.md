# TABLES — all slide numbers (single source of truth)
## Dataset (Slide 10): 20000bp, H12, k16, Total 240901, Distinct 51138; multispecies: E.coli 4.6Mbp/H64/GC50.8, S.cerev 12.1/32/38.5, O.sativa 373/16/43.6, D.mela 143.7/24/42.0, H.sapi Chr21 46.7/254/40.9.
## Static (Slide 18): 1687.44 vs 183.57 ns (9.2x); 0.593 vs 5.448 Mq/s; 16B vs 100+B (6x); ablation 2.0x / 39.2x.
## Dynamic (Slide 19): see results/dynamic_update.csv tables in run.
## FC (Slide 20): 63.6/90.5/98.8/99.3/91.3% faster; tables in results/fractional_cascading.csv.
## Tradeoff (Slide 22): cache 9.2x+6x+39.2x; LSM wins dense, overhead O.sativa -52.7% D.mela -10.0%; FC O(Hlogm)->O(logm+H) 63.6-99.3%.
## VGP (Slide 23): 96.7->2.3min 42x; 4.89->1.32s 3.7x; 0.67->0.79GB 1.18x; scaling series in results/vgp_scaling.csv.
