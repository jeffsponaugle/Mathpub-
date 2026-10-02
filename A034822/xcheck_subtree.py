#!/usr/bin/env python3
"""Independent check of the GPU filter on one subtree: for L=70, h=17, enumerate all
canonical leaves B whose low digits equal PREFIX (fixed), all 8 roots, exact windows,
t = class mod 11, and list (t, root) with S digits h..h+2 mirroring (3-digit filter truth)."""
import sys, math
sys.path.insert(0, '.')
from base_of import canonical
L, h = 70, 17
M = 10**h; H5 = M // 2; e = L - 3*h
m2, m5 = 2**h, 5**h
E_ = 1
while E_ % m2 != m2 - 1: E_ += m5
E_ %= M
def roots8(B):
    out = set()
    for x in (B, B*E_ % M):
        for y in (x, (x+H5) % M):
            out.add(y); out.add((-y) % M)
    return out
def ctz(x): return (x & -x).bit_length() - 1
prefix = int(sys.argv[1]); plen = int(sys.argv[2])   # low plen digits of B fixed
res = set(); leaves = 0
digs0 = [(prefix // 10**i) % 10 for i in range(plen)]
unit5 = digs0[0] == 5
def rec(i, B, s2):
    global leaves
    if i == h:
        leaves += 1
        P = int(str(B*B % M).zfill(h)[::-1])
        tlo = math.isqrt(P*10**e); thi = math.isqrt((P+1)*10**e)
        rs = roots8(B) if s2 else {B, (B+H5) % M, (-B) % M, (H5-B) % M}
        for r in rs:
            for t in range(tlo, thi+1):
                if (t*10**h + r) % 11: continue
                s = str((t*10**h + r)**2)
                if len(s) == L and s[:h] == s[::-1][:h] and s[h:h+3] == s[::-1][h:h+3]:
                    res.add((t, r))
        return
    dmax = 9
    if i == h-1: dmax = 4
    if i == 1 and unit5: dmax = 4
    r = (not unit5) and 1 <= i <= h-2 and B != 0 and ctz(B) == i-1
    if r: dmax = 4
    ds = [digs0[i]] if i < plen else range(dmax+1)
    for d in ds:
        if d > dmax: return
        rec(i+1, B + d*10**i, s2 or r)
rec(0, 0, False)
print("leaves", leaves, "hits", len(res))
with open(sys.argv[3], 'w') as f:
    for t, r in sorted(res): f.write(f"{t} {r}\n")
