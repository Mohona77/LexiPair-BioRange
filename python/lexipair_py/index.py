"""Python mirror of LexiPairIndex (lexipair_index.h): BUILDINDEX + RANGEQUERY."""
from bisect import bisect_left
from .codec import encode_sliding

class LexiPairIndex:
    def __init__(self, K=16):
        self.K=K; self.entries=[]; self.runs=[]; self.run_pos={}
    def build(self, haplotypes):
        E=[]
        for h,seq in enumerate(haplotypes):
            for code,s,e in encode_sliding(seq,self.K):
                E.append((code,s,e,h))
        E.sort(key=lambda x:(x[0],x[1],x[3]))
        self.entries=E; self.runs=[]; self.run_pos={}
        i=0
        while i<len(E):
            j=i+1
            while j<len(E) and E[j][0]==E[i][0]: j+=1
            self.runs.append((E[i][0],i,j-i)); self.run_pos[E[i][0]]=(i,j-i); i=j
    def range_query(self, code,L,R):
        if code not in self.run_pos: return []
        s,c=self.run_pos[code]
        run=self.entries[s:s+c]
        starts=[x[1] for x in run]
        k=bisect_left(starts,L)
        out=[]
        for t in range(k,len(run)):
            if run[t][1]>R: break
            out.append(run[t])
        return out
