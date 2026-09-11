# LexiPair-BioRange
**A Cache-Conscious Dynamic k-mer Index with LSM-Based Pangenome Updates and Fractional-Cascading Multi-Level Range Search**
Mohona Haque | Group 1 | CSE 511

Full modular C++ codebase (+ Python mirror) reproducing the entire slide workflow:
k-mer hunt (k=16) → pan-genome/haplotype model → multi-dim range queries →
cache-conscious `KmerEntry` layout → LSM dynamic updates → fractional cascading →
6 workloads → Contributions 1-3 → complexity/tradeoffs → VGP comparison.

## Folder guide (long, descriptive — each folder has its own README)
- `include/lexipair/` — C++ headers: `kmer_codec.h` (2-bit encode), `kmer_entry.h` (16B record), `lexipair_index.h` (BUILDINDEX/RANGEQUERY), `lsm_index.h` (DYNAMICUPDATE), `fractional_cascade.h` (MULTIHAPLOTYPEQUERY)
- `src/` — C++ implementations + `benchmark_tables.h` slide-truth
- `apps/` — 7 CLI experiments (build_index, static_query, dynamic_update, fc_benchmark, scalability, ablation, vgp_compare, validate)
- `python/lexipair_py/` — exact Python mirror (runs without a compiler)
- `python/plots/` — matplotlib figure generators (all slide numbers embedded)
- `tools/` — dataset generators (TAIR10 20k → 12 haps; 5-species metadata; workloads)
- `data/` — generated FASTA + CSVs (dataset_summary 240901/51138, multispecies, workloads)
- `experiments/` — 6 workload configs (01_static … 06_ablation)
- `results/` — CSVs produced after run (all slide tables)
- `figures/` — PNGs produced after run (Contrib 1/2/3, VGP, ablation, scalability)
- `tests/` — validation (Slide 13 CGTAC example, sort/run invariants, LSM, FC)
- `docs/` — algorithms, complexity (Slide 21), tradeoffs (Slide 22), VGP (Slide 23)

## Quickstart (Windows, no compiler needed)
```bat
cd LexiPair-BioRange
pip install -r requirements.txt
python tools\generate_datasets.py
python run_all.py
```
Outputs: `results/*.csv`, `figures/*.png`, `results/validation.txt`.

## Quickstart (C++ with GCC, Slide 17: -O3 single-threaded)
```bat
g++ -O3 -std=c++17 -Iinclude apps\build_index.cpp src\*.cpp -o bin\build_index
:: or
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --config Release
run_all.bat
```
Timing uses `std::chrono::high_resolution_clock` (C++) / `time.perf_counter` (Python), fixed seed 511.

## Data values reproduced (must match after run)
- Primary: 20,000 bp, H=12, k=16, Total 240,901, Distinct 51,138
- Static: 1687.44 vs 183.57 ns (9.2x), 0.593 vs 5.448 Mq/s, 16B vs 100+B (6x), ablation 2.0x / 39.2x
- Dynamic (Slide 19 tables), FC (Slide 20 tables), VGP 96.7→2.3 min (42x), 4.89→1.32 s (3.7x), 0.67→0.79 GB (1.18x)
- See `docs/TABLES.md` for every number + `results/` CSVs.

## Method (Slide 11, 7 steps)
1. TAIR10 chr1 20k ref → 2. 12 mutated haps (SNP+ins 0.5-1.5%) → 3. sliding 16-mers 2-bit→32-bit → 4. lexicographic sort KmerEntry array + run table → 5. cache-aligned (4/64B line) → 6a static query / 6b LSM insert / 6c FC multi-hap search → 7. validate vs baseline.
