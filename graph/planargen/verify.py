#!/usr/bin/env python3
"""Check geng_coplanar / geng_direct '-v' output against OEIS A049334.

usage: verify.py n {coplanar|direct} stderr_file [rows_json]

A049334 row n lists connected planar graphs with n nodes by edge count k.
For the coplanar variant the generated graphs are complements, so an output
line for k' edges corresponds to k = C(n,2) - k' edges of the planar graph.
Rows through n=13 are known (Georg Grasegger's b-file); for n=14 only
k <= 22 is known, and k = 3n-6 must equal the number of triangulations
A000109(n).
"""
import json, re, sys, os

A000109 = {4: 1, 5: 1, 6: 2, 7: 5, 8: 14, 9: 50, 10: 233, 11: 1249,
           12: 7595, 13: 49566, 14: 339722, 15: 2406841, 16: 17490241}
A003094 = [1, 1, 1, 2, 6, 20, 99, 646, 5974, 71885, 1052805, 17449299,
           313372298, 5942258308]
ROW14_PARTIAL = {13: 3159, 14: 39260, 15: 300748, 16: 1799700, 17: 9123403,
                 18: 40216577, 19: 153876251, 20: 505912342,
                 21: 1416544333, 22: 3360456524}

def main():
    n = int(sys.argv[1]); variant = sys.argv[2]; path = sys.argv[3]
    rows_json = sys.argv[4] if len(sys.argv) > 4 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "a049334_rows.json")
    rows = json.load(open(rows_json)) if os.path.exists(rows_json) else {}
    text = open(path).read()
    counts = {}
    for c, k in re.findall(r">C (\d+) graphs with (\d+) edges", text):
        counts[int(k)] = counts.get(int(k), 0) + int(c)
    total = None
    m = re.search(r">Z (\d+) graphs generated in ([\d.]+) sec", text)
    if m: total = int(m.group(1))
    if variant == "coplanar":
        e = n*(n-1)//2
        counts = {e - k: v for k, v in counts.items()}
    expected = {}
    if str(n) in rows:
        expected = {k: v for k, v in enumerate(rows[str(n)])}
    elif n == 14:
        expected = dict(ROW14_PARTIAL)
    if n in A000109:
        expected[3*n - 6] = A000109[n]
    ok = True
    for k in sorted(set(counts) | set(expected)):
        got, exp = counts.get(k, 0), expected.get(k)
        if exp is None:
            print(f"  k={k:2d}: {got:>14d}   (no reference)")
        else:
            flag = "OK" if got == exp else "MISMATCH"
            if got != exp: ok = False
            print(f"  k={k:2d}: {got:>14d}   expected {exp:>14d}  {flag}")
    s = sum(counts.values())
    ref = A003094[n] if n < len(A003094) else None
    print(f"n={n} {variant}: total={s} (geng >Z says {total})"
          + (f" expected {ref} {'OK' if s == ref else 'MISMATCH'}" if ref else ""))
    if ref is not None and s != ref: ok = False
    if total is not None and total != s: ok = False
    m = re.search(r"in ([\d.]+) sec", text)
    if m: print(f"  cpu {m.group(1)} sec")
    for line in text.splitlines():
        if line.startswith(">C preprune") or line.startswith(">C ") and "us cpu" in line:
            print("  " + line)
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
