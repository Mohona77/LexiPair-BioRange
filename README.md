# 🧬 LexiPair-BioRange

> **Cache-Conscious Indexing for Dynamic Genomic k-mer Coordinate Range Queries in Pan-Genomic Variant Graphs**

<p align="center">
  <img src="https://img.shields.io/badge/Status-Research%20Proposal-blue?style=for-the-badge" />
  <img src="https://img.shields.io/badge/Language-C%2B%2B17-00599C?style=for-the-badge&logo=cplusplus" />
  <img src="https://img.shields.io/badge/Domain-Computational%20Genomics-06D6A0?style=for-the-badge" />
  <img src="https://img.shields.io/badge/Course-Advanced%20Algorithms-FFD166?style=for-the-badge" />
  <img src="https://img.shields.io/badge/License-MIT-green?style=for-the-badge" />
</p>

---

## 📖 Table of Contents

- [Overview](#overview)
- [Motivation](#motivation)
- [Problem Statement](#problem-statement)
- [Research Gap](#research-gap)
- [Proposed Method](#proposed-method)
- [Literature Review](#literature-review)
- [Algorithm Design](#algorithm-design)
- [Course Relevance — Advanced Algorithms](#course-relevance--advanced-algorithms)
- [Publishability](#publishability)
- [Project Timeline](#project-timeline)
- [Repository Structure](#repository-structure)
- [Getting Started](#getting-started)
- [Datasets](#datasets)
- [Baselines & Evaluation](#baselines--evaluation)
- [Expected Results](#expected-results)
- [Team](#team)
- [References](#references)
- [License](#license)

---

## Overview

**LexiPair-BioRange** is a research project proposing the first data structure that jointly addresses three open requirements in pan-genomic indexing:

1. **Exact coordinate range semantics** — given a k-mer, return every genomic interval `[chr:start–end]` across all haplotypes where it occurs.
2. **Cache-conscious memory layout** — organize index data in cache-line-aligned (64-byte) blocks to minimize CPU cache misses and memory stalls during range lookups.
3. **Dynamic updates** — incrementally insert or delete k-mers and coordinate entries as new haplotypes are added to a growing pangenome, without full index reconstruction.

No existing tool addresses all three simultaneously. This project bridges **string indexing theory** (BWT, FM-index, run-length compression) with **geometric range searching** (range trees, fractional cascading) and **cache-oblivious algorithm design** — producing a practically deployable C++ library evaluated on real human and plant pangenomes.

---

## Motivation

### Practical Drivers

Pangenomes now encode hundreds to thousands of haplotypes — the Human Pangenome Reference Consortium (HPRC) alone covers 94 assemblies, and plant pangenomes routinely exceed 200 genomes. The fundamental operation that downstream tools need — *"where, precisely, does this k-mer occur across all coordinate spaces?"* — is not answered by any existing index.

| Use Case | Why It Needs LexiPair-BioRange |
|---|---|
| **Clinical genomics** | Finding all genomic regions matching a disease-associated k-mer across 1,000 patient haplotypes in real time |
| **Population genetics** | Querying billions of k-mers per study; cache-inefficient indexes waste 40–70% of time on memory stalls |
| **Metagenomics** | Identifying strain-level k-mer coordinates across hundreds of reference genomes simultaneously |
| **Dynamic pangenomes** | Clinical pipelines cannot afford full index rebuilds when new assemblies are added weekly |

### Theoretical Drivers

- Formalizes a **new problem class**: lexicographic range queries over genomic coordinate spaces — bridging string indexing with geometric data structures.
- Cache-oblivious/conscious analysis of BWT-based pangenome indexes is **largely unexplored** — this project opens a new subfield of memory-hierarchy-aware genomic data structures.
- Provides **theoretical lower bounds** for dynamic k-mer coordinate range queries, contributing to fundamental compressed data structure theory.

---

## Problem Statement

**Formal Definition:**

Given a pan-genomic variant graph $G$ encoding $H$ haplotypes, and a k-mer query $q$ of length $k$, define a *coordinate range query* as:

$$\text{RangeQuery}(q, L, R) = \{ (h, \text{start}, \text{end}) \mid q \text{ occurs in haplotype } h \text{ at position } p \in [L, R] \}$$

where $h \in \{1 \ldots H\}$, and $[L, R]$ is a genomic window on a reference coordinate system.

**Requirements the solution must satisfy:**

- **Exactness** — no false positives or false negatives in the coordinate output.
- **Output-sensitivity** — query time proportional to the output size, not the full index.
- **Cache efficiency** — $O(|q|/B + \log_B n + \text{output}/B)$ block transfers in the I/O model, where $B$ is the block size.
- **Dynamicity** — support $O(\text{poly}\log\, n)$ amortized update time per inserted haplotype coordinate.

---

## Research Gap

Current tools each address a subset of the problem:

| Tool | Pattern Matching | Coordinate Output | Range Semantics | Dynamic | Cache-Optimized |
|---|:---:|:---:|:---:|:---:|:---:|
| GCSA2 / GBWT | ✅ | ⚠️ node-level | ❌ | ❌ | ❌ |
| MEMO (2025) | ✅ | ❌ membership only | ❌ | ❌ | ❌ |
| Tag Arrays (2025) | ✅ | ✅ one-to-all | ❌ no bounded range | ❌ | ❌ |
| KMC3 / KmerKeys | ✅ | ⚠️ point lookup | ❌ | ❌ | ⚠️ partial |
| BEDTools interval trees | ❌ no k-mer | ✅ | ✅ | ⚠️ | ❌ |
| **LexiPair-BioRange** | ✅ | ✅ exact chr:start–end | ✅ | ✅ | ✅ |

> **Key finding from the literature (2025–2026):** Five recent papers — MEMO, Tag Arrays, Varigraph, Roberts et al., and the 2026 pangenome graph survey — each explicitly identify exact coordinate range queries as an open problem. LexiPair-BioRange is designed to close this gap.

---

## Proposed Method

LexiPair-BioRange is built on four stacked components:

### 1. LexiPair Layer

Each k-mer $q$ is mapped to a **lexicographic pair** $(\ell, h_q)$, where:
- $\ell$ is the BWT-run interval produced by backward search on the pangenome FM-index.
- $h_q$ is a canonical hash of $q$ used for disambiguation in runs with multiple k-mers.

This enables $O(|q|)$ lookup into a sorted coordinate list, inheriting the compression of run-length encoded BWT while adding coordinate awareness.

### 2. Cache-Aligned Coordinate Arrays

Genomic coordinates associated with each k-mer are stored in **64-byte cache-line-aligned sorted arrays**. The layout is designed so that:
- Sequential range scans hit consecutive cache lines.
- Binary search for a boundary $[L, R]$ accesses at most $O(\log n / B)$ cache lines.
- Haplotype coordinate arrays are interleaved by run membership, not by haplotype ID, to maximize spatial locality.

### 3. Dynamic Insertion Protocol

When a new haplotype is added to the pangenome:
1. Identify the BWT runs affected by the new sequences using the existing r-index.
2. Compute only the new coordinate entries for affected k-mers.
3. Merge new entries into the coordinate arrays using a **log-structured merge strategy** (analogous to LSM-trees), deferring compaction to amortize the cost over many insertions.

This achieves $O(\log^2 n)$ amortized update time per new coordinate entry, with no full rebuild required.

### 4. Range Query Engine

Range queries are answered using **fractional cascading** across haplotype layers:
- Pre-process each haplotype's coordinate array with cross-pointers to adjacent arrays.
- For a query $(q, L, R)$: locate $q$ in $O(|q|)$ via the LexiPair layer, then execute a single cascaded binary search across all $H$ haplotype arrays in $O(\log n + H + \text{output})$ time.
- Total I/O complexity: $O(|q|/B + \log_B n + \text{output}/B)$ — optimal for this problem class.

---

## Literature Review

The following recent (2025–2026) papers define the landscape and confirm the research gap:

### Paper 1 — MEMO (2025)
**MEM-based Pangenome Indexing for k-mer Queries**
Hwang, Brown, Ahmed, Jenike, Kovaka, Schatz, Langmead — *Algorithms for Molecular Biology* 20:3, March 2025
🔗 https://doi.org/10.1186/s13015-025-00272-y

| | |
|---|---|
| **Dataset** | HPRC (89 human autosomal haplotypes), HLA locus, Arabidopsis (69 genomes) |
| **Achievement** | MEM-based index; 8.8× smaller than KMC3; 2.5× faster conservation queries; handles arbitrary-length k-mers |
| **Limitation** | Supports membership/conservation queries only — NOT coordinate range queries. Static index; no cache analysis |
| **Relevance** | Sets the compression and speed baseline that LexiPair-BioRange must surpass; MEM anchors could serve as LexiPair entry points |
| **Our gap solution** | Add a cache-optimized range index layer on top of MEM anchors; support exact `[chr:start–end]` output and dynamic updates |

---

### Paper 2 — Tag Arrays (2025/2026)
**Lossless Pangenome Indexing Using Tag Arrays**
Eskandar, Paten, Sirén — *WABI 2025* (LIPIcs Vol. 344) + *bioRxiv* 2026
🔗 https://doi.org/10.4230/LIPIcs.WABI.2025.8

| | |
|---|---|
| **Dataset** | HPRC pangenome graphs (human chromosomes) |
| **Achievement** | BWT tag-array index enables haplotype-aware query; lossless; one-to-all coordinate translation; run-length compressed |
| **Limitation** | No explicit range query API; coordinate translation is one-to-all (not range-bounded); cache efficiency unanalyzed; static |
| **Relevance** | Most directly adjacent to our work — tag arrays provide graph coordinates; LexiPair-BioRange adds range semantics |
| **Our gap solution** | Extend tag-array coordinate maps with a cache-conscious B⁺-tree range index; add a dynamic tag insertion protocol |

---

### Paper 3 — Varigraph (2025)
**Varigraph: An Accurate and Widely Applicable Pangenome Graph-Based Variant Genotyper**
Du, He, Xiao, Hu, Yang, Jiao — *Molecular Plant* 18(9):1587–1601, September 2025
🔗 https://doi.org/10.1016/j.molp.2025.08.001

| | |
|---|---|
| **Dataset** | Human genome, Maize, Wheat, Potato, Pummelo, Rapeseed, Rice (252 genomes) |
| **Achievement** | k-mer genotyping with Bloom filter + bitmap; handles polyploid and repetitive regions; GPU-accelerated; supports 200+ genome graphs |
| **Limitation** | Genotyping-only tool — no coordinate range query. k-mer length capped at 28. No cache analysis. Static graph index |
| **Relevance** | Demonstrates demand for k-mer indexes over very large graphs; bitmap data structures are relevant to our coordinate storage design |
| **Our gap solution** | Replace Bloom filter with cache-aligned interval structure; generalize k-mer length; add dynamic genome insertion support |

---

### Paper 4 — Roberts et al. (2025)
**k-mer-based Approaches to Bridging Pangenomics and Population Genetics**
Roberts, Davis, Josephs, Williamson — *Molecular Biology and Evolution* 42(3): msaf047, March 2025
🔗 https://doi.org/10.1093/molbev/msaf047

| | |
|---|---|
| **Dataset** | Comparative study across plant and human pangenomes |
| **Achievement** | Comprehensive review of k-mer methods for population genomics; distinguishes tools with known relative positions from position-agnostic ones |
| **Limitation** | Review paper — no new index proposed. Explicitly identifies coordinate precision as the key open problem in the field |
| **Relevance** | Directly confirms our research gap: most k-mer tools cannot give exact coordinate ranges across pangenomes |
| **Our gap solution** | LexiPair-BioRange is the direct algorithmic answer to the gap identified by this review |

---

### Paper 5 — Pangenome Graph Survey (2026)
**Pangenome Graphs: Concepts, Tools, and Emerging Trends in Genomic Analysis**
*Journal of Genome Biotechnology and Genetics* 1(1):5, 2026
🔗 https://www.mdpi.com/3042-8424/1/1/5

| | |
|---|---|
| **Dataset** | Review of HPRC, VGP, plant pangenomes; architecture comparison of dBG, VG, Cactus graphs |
| **Achievement** | Comprehensive taxonomy of pangenome graph tools and indexing architectures; identifies challenges for coordinate-level queries at scale |
| **Limitation** | Survey only; no algorithmic contribution. Notes cache-conscious coordinate range indexing is absent from all current architectures |
| **Relevance** | Provides institutional context that this gap exists across the entire field, not just in specific tools |
| **Our gap solution** | LexiPair-BioRange fills the architectural gap identified here with a concrete, deployable cache-conscious range index |

---

## Algorithm Design

### Complexity Summary

| Operation | Time Complexity | I/O Complexity | Space |
|---|---|---|---|
| Build index | $O(n \log n)$ | $O((n/B) \log_{M/B}(n/B))$ | $O(n)$ |
| k-mer lookup (LexiPair) | $O(\|q\|)$ | $O(\|q\|/B)$ | — |
| Coordinate range query | $O(\log n + H + \text{out})$ | $O(\log_B n + \text{out}/B)$ | — |
| Dynamic insert (amortized) | $O(\log^2 n)$ | $O(\log^2_B n)$ | $O(\log n)$ extra |

where $n$ = total coordinate entries, $H$ = number of haplotypes, $B$ = cache block size, $M$ = cache size.

### Data Structure Stack

```
┌─────────────────────────────────────────────────────┐
│              Range Query Engine                      │
│   Fractional cascading across haplotype layers      │
├─────────────────────────────────────────────────────┤
│           Cache-Aligned Coordinate Arrays            │
│   64B-aligned sorted arrays, one per BWT run        │
├─────────────────────────────────────────────────────┤
│               LexiPair Layer                         │
│   BWT-run interval × canonical hash → coord list    │
├─────────────────────────────────────────────────────┤
│         Pangenome FM-index / r-index (base)          │
│   Multi-string BWT over all haplotypes              │
└─────────────────────────────────────────────────────┘
```

### Dynamic Update Flow

```
New haplotype added
        │
        ▼
Identify affected BWT runs
        │
        ▼
Compute new coordinate entries
        │
        ▼
Append to per-run log buffer
        │
        ▼ (when buffer threshold reached)
Log-structured merge into main
coordinate arrays
        │
        ▼
Rebuild fractional cascade pointers
for affected runs only
```

---

## Course Relevance — Advanced Algorithms

This project is submitted as a research project for the **Advanced Algorithms** course at **North South University**. It is deeply grounded in the core theoretical topics of the course:

### 1. Cache-Oblivious & Cache-Conscious Algorithms
The I/O model, cache-line alignment, and the design of cache-oblivious data structures are central course topics. LexiPair-BioRange applies these directly — the coordinate array layout and the fractional cascading query engine are designed with explicit cache complexity analysis, providing a non-trivial real-world application of the theory.

### 2. Compressed & Succinct Data Structures (BWT / FM-index)
The Burrows-Wheeler Transform, FM-index, rank/select structures, and run-length encoded BWT (r-index) are advanced data structure topics covered in the course. LexiPair-BioRange extends these to a new application domain, showing how succinct string indexing interacts with geometric range data structures.

### 3. Range Trees & Fractional Cascading
The range query engine is a direct application of 2D range searching and fractional cascading — a canonical geometric data structure result taught in every advanced algorithms course. This project provides a biologically-motivated implementation of that theory.

### 4. Dynamic Data Structures & Amortized Analysis
The log-structured merge update protocol requires careful amortized analysis. This is a core advanced algorithms technique, and the project provides a realistic, non-trivial setting in which to apply and extend it beyond standard textbook examples.

### 5. Algorithm Engineering & Empirical Evaluation
Benchmarking cache miss rates, query throughput, index size, and dynamic update latency teaches the crucial skill of bridging complexity theory with empirical performance — an explicit learning objective of applied algorithm courses.

---

## Publishability

This study is targeted for submission to top-tier venues in algorithms and bioinformatics. The case for publication rests on five pillars:

**1. Novel Problem Formulation**
No prior work jointly addresses dynamic updates + cache-conscious layout + exact coordinate range semantics for genomic k-mers. The problem is formally new — a necessary condition for acceptance at any venue.

**2. Strong Venue Fit**
Primary targets:
- **WABI** (Workshop on Algorithms in Bioinformatics) — published MEMO (2025) and Tag Arrays (2025)
- **RECOMB** (Research in Computational Molecular Biology)
- **Bioinformatics** (Oxford University Press) — high-impact bioinformatics journal
- **Algorithms for Molecular Biology** (Springer) — published MEMO (2025)
- **Journal of Computational Biology**

**3. Three-Pillar Contribution Structure**
Reviewers at top venues expect: (a) a new data structure or algorithm, (b) theoretical analysis (complexity bounds), and (c) empirical evaluation on real data. LexiPair-BioRange delivers all three.

**4. Demonstrated & Citable Research Gap**
Five 2025–2026 papers explicitly leave the exact coordinate range query problem open. Reviewers can independently verify the gap is real, current, and unsolved — a strong signal for timeliness.

**5. Open Science & Reproducibility**
Full open-source C++ library + benchmark suite on public datasets (HPRC, VGP) meets the reproducibility standards now required by Nature Methods, Genome Biology, and all OUP journals.

---

## Project Timeline

| Phase | Weeks | Tasks |
|---|---|---|
| **Literature & Formalization** | 1–2 | Deep review of MEMO, Tag Arrays, GCSA2, GBWT; formal definition of dynamic coordinate range query; select evaluation datasets |
| **Algorithm Design — LexiPair Layer** | 3–4 | Design lexicographic pairing structure; theoretical analysis (time/space/cache); proof of correctness for range semantics |
| **Cache-Conscious Index Construction** | 5–6 | Implement cache-aligned coordinate arrays (64B-aligned); design BWT-run to coordinate memory layout; C++ prototype + unit tests |
| **Dynamic Update Protocol** | 7–8 | Implement log-structured merge; benchmark insert/delete latency vs. static rebuild baseline; correctness testing on simulated updates |
| **Range Query Engine & Experiments** | 9–10 | Implement fractional cascading; evaluate on HPRC (89 haplotypes) and plant pangenomes; compare vs. MEMO, Tag Arrays, KmerKeys |
| **Evaluation, Writing & Submission** | 11–12 | Full experimental evaluation + ablation studies; write paper; prepare code release, figures, supplementary material |

**Key Deliverables:**
- ✅ Open-source C++ library (this repository)
- ✅ Peer-reviewed paper draft
- ✅ Benchmark suite vs. state-of-the-art
- ✅ Full documentation and dataset scripts

---

## Repository Structure

```
LexiPair-BioRange/
│
├── README.md                        # This file
│
├── src/                             # Core C++ source
│   ├── lexipair/
│   │   ├── lexipair_layer.hpp       # BWT-run interval × canonical hash mapping
│   │   ├── lexipair_layer.cpp
│   │   └── canonical_hash.hpp       # k-mer canonical hash utilities
│   │
│   ├── index/
│   │   ├── coord_array.hpp          # Cache-aligned coordinate arrays
│   │   ├── coord_array.cpp
│   │   ├── bwt_index.hpp            # FM-index / r-index interface
│   │   └── pangenome_index.hpp      # Top-level pangenome index
│   │
│   ├── dynamic/
│   │   ├── lsm_merge.hpp            # Log-structured merge update protocol
│   │   └── lsm_merge.cpp
│   │
│   ├── query/
│   │   ├── range_engine.hpp         # Fractional cascading range query
│   │   ├── range_engine.cpp
│   │   └── query_result.hpp         # Result types: (haplotype, start, end)
│   │
│   └── main.cpp                     # CLI entry point
│
├── include/
│   └── lexipair_biorange.hpp        # Public API header
│
├── tests/
│   ├── unit/
│   │   ├── test_lexipair.cpp
│   │   ├── test_coord_array.cpp
│   │   ├── test_range_query.cpp
│   │   └── test_dynamic_update.cpp
│   └── integration/
│       ├── test_hprc.cpp            # End-to-end on HPRC graph
│       └── test_plant_pangenome.cpp
│
├── benchmarks/
│   ├── bench_query_throughput.cpp   # Queries per second vs. MEMO/Tag Arrays
│   ├── bench_cache_misses.cpp       # Cache miss rate profiling (perf/valgrind)
│   ├── bench_index_size.cpp         # Index size comparison
│   └── bench_dynamic_update.cpp     # Insert latency vs. static rebuild
│
├── scripts/
│   ├── download_hprc.sh             # Download HPRC pangenome graphs
│   ├── download_vgp.sh              # Download VGP data
│   ├── build_index.sh               # End-to-end index construction script
│   └── run_benchmarks.sh            # Full benchmark suite
│
├── data/
│   └── .gitkeep                     # Data downloaded via scripts (not tracked)
│
├── docs/
│   ├── algorithm.md                 # Detailed algorithm description
│   ├── api.md                       # Public API reference
│   ├── benchmarks.md                # Benchmark methodology and results
│   └── figures/                     # Diagrams and plots
│
├── presentation/
│   └── LexiPair-BioRange_Proposal.pptx   # Research proposal slides
│
├── CMakeLists.txt                   # Build system
├── .github/
│   └── workflows/
│       └── ci.yml                   # GitHub Actions CI
└── LICENSE                          # MIT License
```

---

## Getting Started

### Prerequisites

- C++17-compatible compiler (GCC ≥ 10, Clang ≥ 12)
- CMake ≥ 3.18
- [sdsl-lite](https://github.com/simongog/sdsl-lite) — succinct data structures library
- [libdivsufsort](https://github.com/y-256/libdivsufsort) — suffix array construction
- Python ≥ 3.9 (for benchmark scripts and dataset preparation)
- `zlib`, `libbz2` (compression support)

### Build

```bash
git clone https://github.com/[your-username]/LexiPair-BioRange.git
cd LexiPair-BioRange

# Install dependencies (Ubuntu/Debian)
sudo apt-get install cmake g++ zlib1g-dev libbz2-dev

# Build sdsl-lite
git clone https://github.com/simongog/sdsl-lite deps/sdsl-lite
cd deps/sdsl-lite && cmake -DCMAKE_BUILD_TYPE=Release . && make install && cd ../..

# Configure and build
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

### Quick Start

```bash
# 1. Download a small test pangenome (GRCh38 + CHM13)
bash scripts/download_hprc.sh --small

# 2. Build the LexiPair-BioRange index
./build/lexipair_biorange build \
    --graph data/hprc/hprc_v1.1_mc_chm13.gfa \
    --output data/index/hprc.lbr

# 3. Run a coordinate range query
./build/lexipair_biorange query \
    --index data/index/hprc.lbr \
    --kmer ACGTACGTACGTACGT \
    --range chr1:1000000-2000000

# Output format: haplotype_id  chr  start  end  strand
# HG00438.1    chr1  1045823  1045839  +
# HG00621.2    chr1  1387291  1387307  -
# ...
```

### API Usage (C++)

```cpp
#include "lexipair_biorange.hpp"

// Build index from pangenome graph
LexiPairIndex idx;
idx.build("path/to/pangenome.gfa");
idx.save("path/to/output.lbr");

// Load existing index
LexiPairIndex idx;
idx.load("path/to/output.lbr");

// Coordinate range query
std::string kmer = "ACGTACGTACGTACGT";
GenomicRange window{"chr1", 1000000, 2000000};

auto results = idx.range_query(kmer, window);
for (const auto& hit : results) {
    std::cout << hit.haplotype << "\t"
              << hit.chr << "\t"
              << hit.start << "\t"
              << hit.end << "\n";
}

// Dynamic update: add a new haplotype
idx.insert_haplotype("path/to/new_assembly.fa");
```

---

## Datasets

All datasets are publicly available and downloaded via the provided scripts.

| Dataset | Description | Size | Download |
|---|---|---|---|
| **HPRC v1.1** | Human Pangenome Reference Consortium — 89 human autosomal haplotypes (44 individuals + CHM13) | ~100 GB raw | `scripts/download_hprc.sh` |
| **HPRC HLA locus** | Subset of HPRC focused on the Human Leukocyte Antigen region (high variation, benchmark standard) | ~2 GB | Included in HPRC script |
| **VGP (Vertebrate Genomes Project)** | 16 vertebrate species pangenome | ~50 GB | `scripts/download_vgp.sh` |
| **Arabidopsis thaliana pangenome** | 69 genomes; used as MEMO benchmark | ~8 GB | `scripts/download_arabidopsis.sh` |
| **Rice pangenome (252 genomes)** | Used by Varigraph (2025) as stress test | ~120 GB | `scripts/download_rice.sh` |

> **Note:** Raw data is not tracked in this repository. Run the download scripts to populate `data/`.

---

## Baselines & Evaluation

LexiPair-BioRange is evaluated against the following baselines:

| Baseline | Version | What we measure |
|---|---|---|
| [MEMO](https://github.com/StephenHwang/MEMO) | v1.0 (2025) | Index size, k-mer query throughput |
| [Tag Arrays](https://github.com/ucsc-genomics/tag-arrays) | WABI 2025 release | Coordinate retrieval latency |
| [KMC3](https://github.com/refresh-bio/KMC) | v3.2.4 | Index construction time, index size |
| [PanKmer](https://github.com/aylward/PanKmer) | v1.0.3 | Index size, membership query speed |
| [BEDTools](https://github.com/arq5x/bedtools2) | v2.31 | Range query throughput (non-k-mer aware baseline) |

### Metrics

- **Query throughput** — k-mer range queries per second
- **Cache miss rate** — L1/L2/L3 cache misses per query (measured via `perf stat` and Valgrind/Cachegrind)
- **Index size** — bytes per k-mer coordinate entry
- **Construction time** — wall-clock time and peak memory to build the index
- **Update latency** — time to insert one new haplotype (LexiPair-BioRange) vs. full rebuild (static baselines)
- **Correctness** — recall and precision on synthetically generated ground-truth coordinate sets

---

## Expected Results

Based on the theoretical analysis and the limitations documented in baselines, we expect:

- **Query throughput:** 3–5× higher than MEMO on coordinate range queries (cache-aligned access vs. unstructured traversal).
- **Index size:** Comparable to MEMO (2–3 GB for HPRC 89-haplotype), larger than pure membership indexes but with full coordinate output.
- **Cache miss rate:** 40–60% reduction compared to hash-table-based approaches on large-range queries.
- **Update latency:** Log-structured merge achieves sub-second amortized insertion vs. 30–60 minute full rebuilds for static indexes.
- **Exactness:** 100% recall and precision on coordinate outputs by construction (no approximation).

---

## Team

| Role | Name | Institution |
|---|---|---|
| **Student / Researcher** | Mohona Haque | North South University |
| **Supervisor** | Dr. Ahsanur Rahman | North South University |

**Course:** Advanced Algorithms — North South University
**Submission Year:** 2026

---

## References

1. Hwang, S., Brown, N.K., Ahmed, O.Y., Jenike, K.M., Kovaka, S., Schatz, M.C., Langmead, B. (2025). MEM-based pangenome indexing for k-mer queries. *Algorithms for Molecular Biology* 20:3. https://doi.org/10.1186/s13015-025-00272-y

2. Eskandar, P., Paten, B., Sirén, J. (2025). Lossless Pangenome Indexing Using Tag Arrays. *25th International Conference on Algorithms in Bioinformatics (WABI 2025)*, LIPIcs Vol. 344. https://doi.org/10.4230/LIPIcs.WABI.2025.8

3. Du, Z.Z., He, J.B., Xiao, P.X., Hu, J., Yang, N., Jiao, W.B. (2025). Varigraph: An accurate and widely applicable pangenome graph-based variant genotyper for diploid and polyploid genomes. *Molecular Plant* 18(9):1587–1601. https://doi.org/10.1016/j.molp.2025.08.001

4. Roberts, M.D., Davis, O., Josephs, E.B., Williamson, R.J. (2025). k-mer-based approaches to bridging pangenomics and population genetics. *Molecular Biology and Evolution* 42(3): msaf047. https://doi.org/10.1093/molbev/msaf047

5. (2026). Pangenome Graphs: Concepts, Tools, and Emerging Trends in Genomic Analysis. *Journal of Genome Biotechnology and Genetics* 1(1):5. https://www.mdpi.com/3042-8424/1/1/5

6. Sirén, J., Garrison, E., Novak, A.M., Paten, B., Durbin, R. (2020). Haplotype-aware graph indexes. *Bioinformatics* 36(2):400–407.

7. Frigo, M., Leiserson, C.E., Prokop, H., Ramachandran, S. (1999). Cache-oblivious algorithms. *FOCS 1999*, pp. 285–298.

8. Chazelle, B. (1988). Functional approach to data structures and its use in multidimensional searching. *SIAM Journal on Computing* 17(3):427–462. (Fractional cascading)

9. Boucher, C., Gagie, T., Kuhnle, A., Langmead, B., Manzini, G., Mun, T. (2019). Prefix-free parsing for building big BWTs. *WABI 2019*.

10. Liao, W.W., et al. (2023). A draft human pangenome reference. *Nature* 617:312–324. (HPRC)

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2026 Mohona Haque, North South University

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
```

---

<p align="center">
  Made with ❤️ for computational genomics · North South University · 2026
</p>
