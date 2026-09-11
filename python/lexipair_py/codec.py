"""Python mirror of C++ codec (kmer_codec.h). 2-bit encode, K=16 -> 32-bit."""
B2B = {"A":0,"C":1,"G":2,"T":3,"a":0,"c":1,"g":2,"t":3}
B2S = "ACGT"
def encode_kmer(s, K=16):
    assert len(s)==K
    c=0
    for ch in s: c=(c<<2)|B2B[ch]
    return c & 0xFFFFFFFF
def decode_kmer(code, K=16):
    s=""
    for i in range(K): s=B2S[code&3]+s; code>>=2
    return s
def encode_sliding(seq, K=16):
    out=[]
    if len(seq)<K: return out
    c=0
    for i in range(K): c=(c<<2)|B2B[seq[i]]
    out.append((c&0xFFFFFFFF,0,K-1))
    mask=0xFFFFFFFF
    for i in range(K,len(seq)):
        c=(((c<<2)|B2B[seq[i]])&mask)
        out.append((c,i-K+1,i))
    return out
