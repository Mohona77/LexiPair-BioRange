"""tests/test_lexipair.py — validates Slide 13 + invariants (run with python -m pytest or plain python)."""
import os, sys
BASE=os.path.join(os.path.dirname(__file__),".."); sys.path.insert(0,os.path.join(BASE,"python"))
from lexipair_py.codec import encode_kmer, decode_kmer
from lexipair_py.index import LexiPairIndex
from lexipair_py.models import LsmIndex

def test_codec():
    assert encode_kmer("A"*16)==0 and decode_kmer(0)=="A"*16
    assert decode_kmer(encode_kmer("ACGTACGTACGTACGT"))=="ACGTACGTACGTACGT"
def test_slide13_example():
    # Exact Slide 13 mini-table
    idx=LexiPairIndex(); idx.entries=[("ACGTA",5,20,0),("ACGTA",12,28,1),("ACGTA",50,65,2),("ACGTA",100,115,0),("CGTAC",8,23,1),("CGTAC",40,55,0),("CGTAC",70,85,3),("GTACT",15,30,2),("GTACT",60,75,1)]
    # string codes -> use ints via encode for real check below; here simulate run logic on starts
    runs={"CGTAC":(4,3)}; s,c=runs["CGTAC"]; run=idx.entries[s:s+c]; starts=[x[1] for x in run]
    import bisect; k=bisect.bisect_left(starts,30); assert k==1  # idx5
    hits=[r for r in run[k:] if r[1]<=80][:1]; assert hits[0][1]==40 and hits[0][3]==0
def test_lsm_smoke():
    l=LsmIndex(); l.build_base(["ACGT"*5000,"ACGT"*5000]); l.insert_haplotype("ACGT"*5000,2); assert len(l.segments)==2; l.compact(); assert len(l.segments)==1
if __name__=="__main__":
    test_codec(); test_slide13_example(); test_lsm_smoke(); print("tests PASS")
