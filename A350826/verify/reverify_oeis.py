#!/usr/bin/env python3
"""
reverify_oeis.py - re-verify OEIS data related to prime sextuplets with ../a350826.

Checks (results in reverify_oeis.txt):
  * A022008 b-file (n = 1..10000) and A271000 b-file (5940 terms) against `a350826 list`;
  * A200503 / A200504 / A233426 b-files (record gaps between sextuplets) for all records
    whose end lies below 1e15, from the sorted list of all 12,432,129 sextuplets < 1e15;
  * A343636(n) for 0 <= n <= 23 (smallest (n+1)-digit initial member);
  * no sextuplet straddles 10^n (3 <= n <= 24) or 2^k (60 <= k <= 79): no initial member
    in [x - 16, x).
Run from this directory after `make` in the parent; needs the b-files b*.txt here.
"""
import subprocess
import sys

TOOL = '../a350826'
OUT = open('reverify_oeis.txt', 'w')


def log(s):
    print(s, flush=True)
    OUT.write(s + '\n')
    OUT.flush()


def bfile(name):
    vals = []
    for line in open(name):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        n, v = line.split()[:2]
        vals.append((int(n), int(v)))
    return vals


def tool_list(lo, hi):
    r = subprocess.run(['nice', '-n', '10', TOOL, 'list', str(lo), str(hi)], capture_output=True, text=True, check=True)
    return [int(x) for x in r.stdout.split()]


def tool_count(lo, hi):
    r = subprocess.run([TOOL, 'count', str(lo), str(hi), '-q'], capture_output=True, text=True, check=True)
    for line in r.stdout.splitlines():
        if line.startswith('RESULT [') and ' count ' in line:
            return int(line.split(' count ')[1].split()[0])
    raise RuntimeError(r.stdout)


OFF = (0, 4, 6, 10, 12, 16)

# --- A022008 and A271000 ---------------------------------------------------------------
b = bfile('b022008.txt')
lst = tool_list(0, b[-1][1] + 1)
ok = [v for _, v in b] == lst[:len(b)] and len(lst) == len(b)
log(f"A022008 b-file n=1..{len(b)}: {'OK' if ok else 'MISMATCH'} ({len(lst)} sextuplets below {b[-1][1] + 1})")
b = bfile('b271000.txt')
members = [p + d for p in lst for d in OFF]
ok = [v for _, v in b] == members[:len(b)]
log(f"A271000 b-file n=1..{len(b)}: {'OK' if ok else 'MISMATCH'}")

# --- record gaps below 1e15 --------------------------------------------------------------
LIM = 10 ** 15
lst = tool_list(0, LIM)
log(f"sextuplets below 1e15: {len(lst)} (A063501(15) = 12432129)")
rec = []                     # (gap, start, end)
best = 0
for a, c in zip(lst, lst[1:]):
    if c - a > best:
        best = c - a
        rec.append((best, a, c))
for name, idx in (('A200503', 0), ('A200504', 1), ('A233426', 2)):
    b = bfile(f'b{name[1:]}.txt')
    k = min(len(rec), len(b))
    ok = all(rec[i][idx] == b[i][1] for i in range(k))
    log(f"{name} b-file n=1..{k}: {'OK' if ok else 'MISMATCH'} ({len(rec)} records with end < 1e15, b-file has {len(b)} terms)")
log(f"largest record gap below 1e15: {rec[-1][0]} from {rec[-1][1]} to {rec[-1][2]}")

# --- A343636 -------------------------------------------------------------------------------
b = dict(bfile('b343636.txt'))
bad = 0
for n in range(0, 24):
    lo, hi = 10 ** n, 10 ** (n + 1)
    first, w = None, 10 ** 9
    while first is None:
        top = min(hi, lo + w)
        got = tool_list(lo, top)
        if got:
            first = got[0]
        elif top == hi:
            break
        w *= 4
    a = first - lo if first is not None else 0
    good = a == b.get(n)
    bad += not good
    log(f"A343636({n}) = {a} {'OK' if good else 'MISMATCH (b-file ' + str(b.get(n)) + ')'}")
log(f"A343636 n=0..23: {'OK' if bad == 0 else str(bad) + ' MISMATCHES'}")

# --- straddling ------------------------------------------------------------------------------
bad = 0
for n in range(3, 25):
    x = 10 ** n
    c = tool_count(x - 16, x)
    bad += c != 0
for k in range(60, 80):
    x = 2 ** k
    c = tool_count(x - 16, x)
    bad += c != 0
log(f"no sextuplet with initial member in [x-16, x) for x = 10^3..10^24 and 2^60..2^79: {'OK' if bad == 0 else 'FAILED'}")
