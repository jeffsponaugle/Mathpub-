#!/usr/bin/env python3
"""Format results as OEIS b-file lines.

usage: bfile.py rows_with_nN.json N

Prints the A049334 b-file lines for row N (terms are numbered sequentially:
rows 1..13 occupy terms 1..212), the A003094 line for n=N, and the A005470
value for n=N via euler.py's transform.
"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from euler import euler_bivariate

rows = json.load(open(sys.argv[1])); N = int(sys.argv[2])
missing = [n for n in range(1, N + 1) if str(n) not in rows]
if missing:
    sys.exit(f"bfile.py: rows {missing} missing from {sys.argv[1]}; need every row 1..{N}")
start = 1 + sum(len(rows[str(n)]) for n in range(1, N))
print(f"# A049334 row {N}: terms {start}..{start + len(rows[str(N)]) - 1}")
for i, v in enumerate(rows[str(N)]):
    print(start + i, v)
print(f"# A003094: {N} {sum(rows[str(N)])}")
T = {int(n): {k: v for k, v in enumerate(r) if v} for n, r in rows.items() if int(n) <= N}
U = euler_bivariate(T, N)
print(f"# A005470: {N} {sum(U[N].values())}")
kmax = max(U[N])
print(f"# A039735 row {N} (k=0..{kmax}): {[U[N].get(k, 0) for k in range(kmax + 1)]}")
