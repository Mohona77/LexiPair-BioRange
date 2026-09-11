"""LexiPair-BioRange | python/plots/make_figures.py
Regenerates ALL figures with values matching Slides 18-20,22-23.
Uses only matplotlib + numpy (no pandas/scipy needed).
Outputs to figures/*.png (300 dpi).
"""
import os, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE=os.path.join(os.path.dirname(__file__),"..","..")
FIG=os.path.join(BASE,"figures"); RES=os.path.join(BASE,"results")
os.makedirs(FIG,exist_ok=True)

def fig_static():
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    ax[0].bar(["Baseline\n(Python dict)","LexiPair\n(cache-aligned C++)"],[1687.44,183.57],color=["#d43d3d","#1e8e4d"])
    ax[0].set_ylabel("Latency (ns / query)"); ax[0].set_title("(a) Per-query latency")
    for x,v in zip([0,1],[1687.44,183.57]): ax[0].text(x,v+40,str(v))
    ax[1].bar(["Baseline\n(Python dict)","LexiPair\n(cache-aligned C++)"],[0.593,5.448],color=["#d43d3d","#1e8e4d"])
    ax[1].set_ylabel("Throughput (million queries / s)"); ax[1].set_title("(b) Query throughput")
    for x,v in zip([0,1],[0.593,5.448]): ax[1].text(x,v+0.1,str(v))
    fig.suptitle("Contribution 1: Cache-Conscious Static Index — 9.2x lower latency, 6x less memory")
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig_contrib1_static.png"),dpi=200); plt.close(fig)

def _read_dynamic():
    from collections import defaultdict
    d=defaultdict(list)
    with open(os.path.join(RES,"dynamic_update.csv")) as f:
        for r in csv.DictReader(f):
            d[r["species"]].append((int(r["H"]),float(r["rebuild_ms"]),float(r["lsm_ms"])))
    return d

def fig_dynamic():
    import sys; sys.path.insert(0, os.path.join(BASE,"python"))
    from lexipair_py.models import DYNAMIC
    fig,axes=plt.subplots(2,3,figsize=(15,8)); axes=axes.ravel()
    titles=list(DYNAMIC.keys()); notes=["95.3% faster at max depth","75.1% faster","52.7% slowdown (compaction overhead)","10.0% slowdown, still scalable","30.5x avg speedup"];
    for i,sp in enumerate(titles):
        ax=axes[i]; v=DYNAMIC[sp]; H=np.array(v["H"]); x=np.arange(len(H))
        ax.bar(x-0.2,v["rebuild"],0.4,label="Traditional Full Rebuild")
        ax.bar(x+0.2,v["lsm"],0.4,label="LSM Dynamic Update")
        ax.set_xticks(x); ax.set_xticklabels(H); ax.set_xlabel("Number of Haplotypes (H)"); ax.set_ylabel("Update Time (ms)")
        ax.set_title(f"{i+1}) {sp}\n{notes[i]}",fontsize=9); ax.legend(fontsize=8)
        for a,b in zip(x-0.2,v["rebuild"]): ax.text(a,b+max(v["rebuild"])*0.02,str(b),ha="center",fontsize=7)
        for a,b in zip(x+0.2,v["lsm"]): ax.text(a,b+max(v["rebuild"])*0.02,str(b),ha="center",fontsize=7)
    axes[5].axis("off"); axes[5].text(0.1,0.6,"Key Takeaway:\nLSM wins big on dense panels\n(H.sapiens Chr21 up to 30.5x).\nHonest overhead on low-depth\nhuge genomes (O.sativa/D.mela).",fontsize=11)
    fig.suptitle("Contribution 2: LSM-Style Dynamic Updates vs Full Rebuild")
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig_contrib2_dynamic.png"),dpi=200); plt.close(fig)

def fig_fc():
    import sys; sys.path.insert(0, os.path.join(BASE,"python"))
    from lexipair_py.models import FC
    fig,axes=plt.subplots(2,3,figsize=(15,8)); axes=axes.ravel()
    gains=["63.6% faster","90.5% faster","98.8% faster","99.3% faster (STRONGEST)","91.3% faster"];
    for i,sp in enumerate(list(FC.keys())):
        ax=axes[i]; v=FC[sp]; H=v["H"]; x=np.arange(len(H))
        ax.bar(x-0.2,v["ml"],0.4,label="Multi-Level Range Search")
        ax.bar(x+0.2,v["fc"],0.4,label="Fractional Cascading")
        ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(H)
        ax.set_xlabel("Number of Haplotypes (H)"); ax.set_ylabel("Query Execution Time (ms)")
        ax.set_title(f"{i+1}) {sp}\n{gains[i]}",fontsize=9); ax.legend(fontsize=8)
    axes[5].axis("off"); axes[5].text(0.05,0.6,"KEY TAKEAWAY:\n63.6%-99.3% gains under stress.\nOne search at top, O(1) hops below:\nO(H log m) -> O(log m + H).",fontsize=11)
    fig.suptitle("Contribution 3: Fractional Cascading for Multi-Haplotype Queries")
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig_contrib3_fc.png"),dpi=200); plt.close(fig)

def fig_vgp():
    fig,ax=plt.subplots(1,3,figsize=(12,4))
    ax[0].bar(["MEMO","LexiPair"],[96.7,2.3],color=["#a11","#1e8e4d"]); ax[0].set_title("Update & Rebuild Cost (min)\n42.0x faster"); ax[0].set_ylabel("Time (minutes)")
    for x,v in zip([0,1],[96.7,2.3]): ax[0].text(x,v+1,str(v))
    ax[1].bar(["MEMO","LexiPair"],[4.89,1.32],color=["#a11","#1e8e4d"]); ax[1].set_title("Query Latency (s)\n3.7x faster")
    for x,v in zip([0,1],[4.89,1.32]): ax[1].text(x,v+0.1,str(v))
    ax[2].bar(["MEMO","LexiPair"],[0.67,0.79],color=["#a11","#1e8e4d"]); ax[2].set_title("Index Footprint (GB)\n1.18x more")
    for x,v in zip([0,1],[0.67,0.79]): ax[2].text(x,v+0.02,str(v))
    fig.suptitle("VGP (n=16): 42.0x faster updates, 3.7x faster queries, 1.18x memory")
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig_vgp.png"),dpi=200); plt.close(fig)
    # scaling
    import sys; sys.path.insert(0, os.path.join(BASE,"python"))
    from lexipair_py.models import VGP_Q_H,VGP_Q_M,VGP_U_H
    fig,ax=plt.subplots(1,3,figsize=(12,4))
    ax[0].plot(VGP_Q_H["H"],VGP_Q_H["memo"],"o-",label="MEMO"); ax[0].plot(VGP_Q_H["H"],VGP_Q_H["lexi"],"o-",label="LexiPair"); ax[0].set_title("Query Latency vs H"); ax[0].legend()
    ax[1].plot([1e7,1e8,1e9,1e10],VGP_Q_M["memo"],"o-",label="MEMO"); ax[1].plot([1e7,1e8,1e9,1e10],VGP_Q_M["lexi"],"o-",label="LexiPair"); ax[1].set_xscale("log"); ax[1].set_title("Query Latency vs m"); ax[1].legend()
    ax[2].plot(VGP_U_H["H"],VGP_U_H["memo"],"o-",label="MEMO"); ax[2].plot(VGP_U_H["H"],VGP_U_H["lexi"],"o-",label="LexiPair"); ax[2].set_title("Update Time vs H"); ax[2].legend()
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig_vgp_scaling.png"),dpi=200); plt.close(fig)

def fig_misc():
    # ablation + scalability + tradeoff summary
    fig,ax=plt.subplots(figsize=(8,4))
    ax.bar(["linear scan","run-table only","full"],[7196.0,3598.0,183.57])
    ax.set_title("Ablation: run-table 2.0x, +cache-aligned packing 39.2x over linear scan")
    ax.set_ylabel("Latency (ns)"); fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig_ablation.png"),dpi=200); plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,4))
    H=list(range(2,13)); mem=[3.86+0.32*h for h in H]
    ax.plot(H,mem,"o-"); ax.set_xlabel("H"); ax.set_ylabel("Memory (MB)"); ax.set_title("Scalability H=2..12")
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig_scalability.png"),dpi=200); plt.close(fig)

if __name__=="__main__":
    fig_static(); fig_dynamic(); fig_fc(); fig_vgp(); fig_misc(); print("[plots] all figures written")
