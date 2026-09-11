# Experiment 01 — Static Query (Workload 1, Slide 17-18)
500 random k-mers x20 iters, latency + memory. Config: {"n_queries":500,"iterations":20,"seed":511}
Run: python run_all.py (or bin/static_query). Output: results/static_query.csv (1687.44 vs 183.57, 9.2x), figures/fig_contrib1_static.png.
