#!/usr/bin/env python3
"""Compare per-edge planar-graph counts from the naive pipeline
(geng ... | planarg -q | countg --e) with geng_coplanar_lr -v output.

usage: compare_naive.py n NPIECES dir
  dir/naive_<r>.txt : countg --e output for piece r (ends with 'altogether')
  dir/ours_<r>.txt  : geng_coplanar_lr -v output for piece r (has a '>Z' line)
"""
import re, sys, os

n, npieces, d = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
E = n*(n-1)//2
naive, ours, ndone, odone = {}, {}, 0, 0
for r in range(npieces):
    f = os.path.join(d, f"naive_{r}.txt")
    if os.path.exists(f):
        t = open(f).read()
        if "altogether" in t:
            ndone += 1
            for c, k in re.findall(r"(\d+) graphs : e=(\d+)", t):
                naive[int(k)] = naive.get(int(k), 0) + int(c)
    f = os.path.join(d, f"ours_{r}.txt")
    if os.path.exists(f):
        t = open(f).read()
        if ">Z" in t:
            odone += 1
            for c, k in re.findall(r">C (\d+) graphs with (\d+) edges", t):
                ours[E - int(k)] = ours.get(E - int(k), 0) + int(c)
print(f"pieces complete: naive {ndone}/{npieces}, ours {odone}/{npieces}")
ok = ndone == npieces and odone == npieces
for k in sorted(set(naive) | set(ours)):
    a, b = naive.get(k, 0), ours.get(k, 0)
    flag = "OK" if a == b else "MISMATCH"
    if a != b: ok = False
    print(f"  k={k:2d}: naive {a:>12d}   ours {b:>12d}   {flag}")
print(f"totals: naive {sum(naive.values())}  ours {sum(ours.values())}")
print("RESULT:", "PASS" if ok else ("INCOMPLETE" if ndone < npieces or odone < npieces else "FAIL"))
