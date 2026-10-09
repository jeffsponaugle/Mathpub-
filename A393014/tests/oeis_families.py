#!/usr/bin/env python3
"""Search several base pairs up to the OEIS b-file range and require identical term lists.

This is the independent check (the brute-force self-tests share code with the engines):
power-of-two bases (2,4) (2,8) (4,8) (4,10) (8,10) (10,16), coprime pairs (2,9) (3,10) (9,10)
(7,10) (2,3) (4,7) (4,9), and (2,10) up to 10^30.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import oeis  # noqa: E402
from campaign import roles  # noqa: E402

TESTS = [("A029804", 8, 10, None), ("A029961", 4, 10, None), ("A029731", 10, 16, None), ("A097856", 2, 4, None),
         ("A259380", 2, 8, None), ("A259382", 4, 8, None), ("A259387", 4, 9, None), ("A259378", 4, 7, None),
         ("A259385", 2, 9, None), ("A007633", 3, 10, 26), ("A029965", 9, 10, None), ("A029964", 7, 10, None),
         ("A060792", 2, 3, 71), ("A007632", 2, 10, 30)]


def ndig(n, b):
    d = 0
    while n:
        n //= b
        d += 1
    return d


def main():
    bad = 0
    for anum, b1, b2, Lcap in TESTS:
        known = sorted(set(oeis.all_terms(oeis.entry(anum))))
        P, Q = roles(b1, b2)
        L = ndig(max(known), P) - 1
        if Lcap:
            L = min(L, Lcap)
        lim = P ** L
        dense = len(known) > 1000
        exe = os.path.join(ROOT, "dbpal-cpu" if dense else "dbpal")
        cmd = [exe, "search", "-b", f"{b1},{b2}", "--Lmin", "1", "--Lmax", str(L)] + (["--engine", "cpu"] if dense else [])
        out = subprocess.run(cmd, capture_output=True, text=True)
        mine = sorted(int(l.split("\t")[0]) for l in out.stdout.splitlines() if l[:1].isdigit())
        kb = [n for n in known if n < lim]
        ok = out.returncode == 0 and mine == kb
        bad += not ok
        print(f"   {anum} bases {b1},{b2} below {P}^{L}: {len(mine)} found, OEIS {len(kb)}  {'OK' if ok else 'DIFFERENT'}")
        if not ok:
            print("      missing", sorted(set(kb) - set(mine))[:8], "extra", sorted(set(mine) - set(kb))[:8])
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
