#!/usr/bin/env python3
"""Compare ./lego -SV symmetric counts with Deleuran's table.

usage: sym_vs_lasse.py SYM_FILE DELEURAN_SUM_FOR_SIZE_PY

Deleuran lists, for every bottleneck-free refinement, the number of classes
fixed by a 180-degree rotation; ours is (S180 + S90) / 2 per profile.
"""
import re, sys

sym = {}
for line in open(sys.argv[1]):
    m = re.match(r'\s+<(\d+)> S180=(\d+) S90=(\d+)', line)
    if m:
        sym[m.group(1)] = (int(m.group(2)) + int(m.group(3))) // 2

src = open(sys.argv[2]).read()
ok = bad = missing = 0
for k, tot, s in re.findall(r"'(\d+)',\s*(\d+),\s*\(?(\d+)\)?", src):
    if len(k) < 2:
        continue
    n = sum(int(c) for c in k)
    ours = sym.get(k, 0)
    if not any(p for p in sym if sum(int(c) for c in p) == n):
        missing += 1
        continue
    if ours == int(s):
        ok += 1
    else:
        bad += 1
        print('MISMATCH <%s> (n=%d): ours %d, Deleuran %s' % (k, n, ours, s))
print('symmetric counts: %d agree, %d differ, %d not compared (size not in SYM_FILE)' % (ok, bad, missing))
