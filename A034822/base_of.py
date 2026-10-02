#!/usr/bin/env python3
"""Map a root n (with n^2 a palindrome of even length L) to the GPU/CPU base index that
must find it: the canonical representative of the 8 square roots of n^2 mod 10^h
(unit 1..5, unit 5 => tens < 5, sigma2 digit < 5, top digit < 5)."""
import sys
def params(L):
    D=(L+1)//2
    h=17 if L>=68 else None
    return h
def ib_for(L,H):
    e=L-3*H
    for ib in range(2,H):
        lg=(H+e)/2+3*(1-ib)-1.2041199826559248
        if lg<-8 and H-ib<=7: return ib
    return H-1
def canonical(B,h):
    M=10**h; H5=M//2
    m2,m5=2**h,5**h
    e=1
    while e%m2!=m2-1: e+=m5
    E=e%M
    cands=set()
    for x in (B, B*E%M):
        for y in (x, (x+H5)%M):
            cands.add(y); cands.add((-y)%M)
    good=[]
    for x in sorted(cands):
        u=x%10
        if u<1 or u>5: continue
        d=[(x//10**i)%10 for i in range(h)]
        if u==5 and d[1]>4: continue
        if d[h-1]>4: continue
        ok=True
        if u!=5:
            for i in range(1,h-1):
                Bi=x%10**i
                if Bi and (Bi & -Bi).bit_length()-1==i-1 and d[i]>4: ok=False
        if ok: good.append(x)
    return good
if __name__=="__main__":
    L=int(sys.argv[1]); n=int(sys.argv[2]); h=int(sys.argv[3]) if len(sys.argv)>3 else L//4
    IB=ib_for(L,h)
    for x in canonical(n%10**h,h):
        u=x%10; rest=0
        for i in range(1,IB): rest+=((x//10**i)%10)*10**(i-1)
        print(L,h,IB,"canonical",x,"base",(u-1)*10**(IB-1)+rest)
