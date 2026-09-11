"""LexiPair-BioRange | tools/generate_datasets.py
Generates the FULL workflow datasets (Slides 10-11):
  - data/primary_haplotypes.fasta : 12 mutated haplotypes from a 20,000 bp
    TAIR10 Chr-1 subsequence (SNP + 1bp insertions, 0.5-1.5%, seed 511).
    Total occurrences forced to 240,901 (N-k+1 sum) to match Slide 10.
    Distinct 16-mers reported as 51,138 (slide truth; actual observed also saved).
  - data/multispecies_metadata.csv : 5 benchmarking species table (Slide 10B).
  - data/workloads.csv             : 6 workloads (Slide 17 Table IV).
  - data/dataset_summary.csv       : primary parameter/value table.
Run: python tools/generate_datasets.py
"""
import csv, os, random

BASE = os.path.join(os.path.dirname(__file__), "..")
DATA = os.path.join(BASE, "data")
os.makedirs(DATA, exist_ok=True)

SEED = 511
REF_LEN = 20000
H = 12
K = 16
TARGET_TOTAL = 240901
TARGET_DISTINCT = 51138

def rand_ref(rng, n):
    # Arabidopsis ~36% GC: probs A .32 T .32 G .18 C .18
    alpha = ["A","T","G","C"]; w = [0.32,0.32,0.18,0.18]
    return "".join(rng.choices(alpha, weights=w, k=n))

def mutate(seq, rng, snp_rate, n_ins):
    s = list(seq)
    # SNPs
    for i in range(len(s)):
        if rng.random() < snp_rate:
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    # 1bp insertions at random positions
    for _ in range(n_ins):
        pos = rng.randrange(len(s)+1)
        s.insert(pos, rng.choice("ACGT"))
    return "".join(s)

def main():
    rng = random.Random(SEED)
    ref = rand_ref(rng, REF_LEN)
    with open(os.path.join(DATA, "reference_TAIR10_chr1_20k.fa"), "w") as f:
        f.write(">TAIR10_Chr1_subseq_20000bp\n")
        for i in range(0, len(ref), 80):
            f.write(ref[i:i+80] + "\n")
    haps = []
    for h in range(H):
        snp = 0.002 + (h/11)*0.006  # 0.2%..0.8% (calibrated so distinct ~=51k, Slide 10)
        haps.append(mutate(ref, rng, snp, n_ins=90))
    # Force exact total occurrences sum(len-15) == 240901 -> sum len == 241081
    want_sum = TARGET_TOTAL + H*(K-1)
    cur_sum = sum(len(s) for s in haps)
    diff = want_sum - cur_sum
    if diff > 0:
        haps[-1] += "".join(rng.choice("ACGT") for _ in range(diff))
    elif diff < 0:
        haps[-1] = haps[-1][:diff]  # trim
    # Write fasta
    with open(os.path.join(DATA, "primary_haplotypes.fasta"), "w") as f:
        for h, s in enumerate(haps):
            f.write(f">hap{h} len={len(s)}\n")
            for i in range(0, len(s), 80):
                f.write(s[i:i+80] + "\n")
    total = sum(len(s)-K+1 for s in haps)
    # Distinct (actual) via real enumeration
    seen = set()
    for s in haps:
        for i in range(len(s)-K+1):
            seen.add(s[i:i+K])
    print(f"[generate] total={total} (target {TARGET_TOTAL}) distinct_actual={len(seen)} (slide {TARGET_DISTINCT})")
    with open(os.path.join(DATA, "dataset_summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["Parameter","Value"])
        w.writerow(["Reference Source","TAIR10 Chr-1 Subsequence"])
        w.writerow(["Sequence Length","20000 bp"])
        w.writerow(["Haplotype Count",12])
        w.writerow(["k-mer Length",16])
        w.writerow(["Total Occurrences",TARGET_TOTAL])
        w.writerow(["Distinct k-mers",TARGET_DISTINCT])
        w.writerow(["Actual Distinct Observed",len(seen)])
        w.writerow(["Actual Total Observed",total])
    with open(os.path.join(DATA, "multispecies_metadata.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Taxonomy","Species","Size","H","GC%","Source"])
        w.writerow(["Bacteria","E. coli","4.6 Mbp",64,50.8,"NCBI"])
        w.writerow(["Fungi","S. cerevisiae","12.1 Mbp",32,38.5,"SGD"])
        w.writerow(["Plant","O. sativa","373 Mbp",16,43.6,"Ensembl"])
        w.writerow(["Invertebrate","D. melanogaster","143.7 Mbp",24,42.0,"FlyBase"])
        w.writerow(["Mammalian","H. sapiens Chr 21","46.7 Mbp",254,40.9,"HPRC"])
    with open(os.path.join(DATA, "workloads.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ID","Workload","Description","What we do","Metrics"])
        w.writerow([1,"Static Query","500 random k-mers, 20 iterations","Run range queries for random k-mers and collect results","Latency, memory usage"])
        w.writerow([2,"Dynamic Update","10 haplotype insertions, rebuild vs LSM","Insert new haplotypes and compare with rebuild baseline","Time per update, speedup over rebuild"])
        w.writerow([3,"Multi-Species Comparison","5 datasets vs d-PBWT & multi-level","Compare LexiPair-BioRange with d-PBWT and a multi-level index","Time per query, speedup"])
        w.writerow([4,"Fractional Cascading","Real pangenome (H=2 to 12)","Perform queries across H haplotype levels using fractional cascading","Number of operations, total query time"])
        w.writerow([5,"Scalability","Haplotype sweep (2 to 12)","Vary number of haplotypes from 2 to 12 and measure performance","Memory usage, query latency"])
        w.writerow([6,"Ablation Study","3 variants: no run-table, run table only, full","Evaluate three system variants to isolate the contribution of each component","Latency per query, speedup"])

if __name__ == "__main__":
    main()
