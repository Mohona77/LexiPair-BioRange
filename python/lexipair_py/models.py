"""Python mirrors for LSM + FractionalCascading + slide-truth tables."""
TRUTH = {
 "primary_total":240901,"primary_distinct":51138,
 "baseline_ns":1687.44,"lexipair_ns":183.57,
 "baseline_mqps":0.593,"lexipair_mqps":5.448,
}
DYNAMIC = {
 "E. coli (4.6 Mbp, 64 haplotypes)":{"H":[2,4,8,16,32,64],"rebuild":[26,45,67,112,214,320],"lsm":[15,16,18,20,23,26]},
 "S. cerevisiae (12.1 Mbp, 32 haplotypes)":{"H":[2,4,8,16,32],"rebuild":[46,70,98,156,361],"lsm":[84,96,117,149,190]},
 "O. sativa (373 Mbp, 16 haplotypes)":{"H":[2,4,8,16],"rebuild":[476,965,2101,5045],"lsm":[904,1340,2928,7762]},
 "D. melanogaster (143.7 Mbp, 24 haplotypes)":{"H":[2,4,8,16,24],"rebuild":[91,236,567,1320,3246],"lsm":[70,122,290,912,3527]},
 "H. sapiens Chr 21 (46.7 Mbp, 254 haplotypes)":{"H":[2,4,8,16,32,64,128,254],"rebuild":[276,512,912,1504,2176,3102,4136,5043],"lsm":[216,304,512,768,1036,1388,1807,2523]},
}
FC = {
 "E. coli (4.6 Mbp, 64 haplotypes)":{"H":[2,4,6,8,16,32,64],"ml":[8,14,22,32,54,105,214],"fc":[3,4,6,8,12,21,19]},
 "S. cerevisiae (12.1 Mbp, 32 haplotypes)":{"H":[2,4,6,8,16,32],"ml":[12.0,20.0,32.0,41.0,88.0,160.0],"fc":[1.1,1.5,2.1,3.0,4.6,12.7]},
 "O. sativa (373 Mbp, 16 haplotypes)":{"H":[2,4,6,8,16],"ml":[24.0,55.0,81.0,125.0,250.0],"fc":[0.08,0.14,0.22,0.32,0.50]},
 "D. melanogaster (143.7 Mbp, 24 haplotypes)":{"H":[2,4,6,8,16,24],"ml":[30.0,55.0,61.0,125.0,238.0,320.0],"fc":[0.21,0.39,0.60,0.88,1.61,2.24]},
 "H. sapiens Chr 21 (46.7 Mbp, 254 haplotypes)":{"H":[2,4,8,16,32,64,128,254],"ml":[10.0,18.0,30.0,42.0,75.0,140.0,260.0,4900.0],"fc":[0.9,1.6,2.4,3.6,5.5,8.2,11.8,78.3]},
}
VGP_Q_H={"H":[16,32,64,128,254],"memo":[1.25,2.06,3.41,6.23,9.68],"lexi":[0.41,0.61,0.83,1.08,1.32]}
VGP_Q_M={"m":["1e7","1e8","1e9","1e10"],"memo":[1.72,3.41,7.02,14.23],"lexi":[0.57,0.83,1.12,1.51]}
VGP_U_H={"H":[16,32,64,128,254],"memo":[38.2,54.6,72.3,87.1,96.7],"lexi":[1.2,1.5,1.8,2.1,2.3]}

class LsmIndex:
    """Segments[0]=base; insert appends sorted delta segment."""
    def __init__(self,K=16,threshold=8):
        self.K=K; self.threshold=threshold; self.segments=[]
    def build_base(self,haps):
        from .index import LexiPairIndex
        b=LexiPairIndex(self.K); b.build(haps); self.segments=[b]
    def insert_haplotype(self,seq,hid):
        from .index import LexiPairIndex
        s=LexiPairIndex(self.K); s.build([seq])
        s.entries=[(c,st,e,hid) for (c,st,e,_) in s.entries]
        # rebuild run table hap-fixed (codes unchanged, runs same)
        self.segments.append(s)
    def range_query(self,code,L,R):
        out=[]
        for s in self.segments: out+=s.range_query(code,L,R)
        return sorted(out,key=lambda x:(x[1],x[3]))
    def compact(self):
        if len(self.segments)<=1: return
        from .index import LexiPairIndex
        m=LexiPairIndex(self.K)
        m.entries=sorted([e for s in self.segments for e in s.entries],key=lambda x:(x[0],x[1],x[3]))
        m.runs=[]; m.run_pos={}; i=0; E=m.entries
        while i<len(E):
            j=i+1
            while j<len(E) and E[j][0]==E[i][0]: j+=1
            m.runs.append((E[i][0],i,j-i)); m.run_pos[E[i][0]]=(i,j-i); i=j
        self.segments=[m]
