#!/usr/bin/env python3
"""
Independent check of the a335406 tool.

Sieves the primes up to N with numpy (odd-only Eratosthenes), computes the prime
gaps, their runs and everything the tool reports -- pi(N), E (primes whose gap
equals the previous gap), the number of runs, the run-length histogram, the
(length, gap) table and the first occurrence of every run length -- with plain
Python/numpy logic that shares no code with the C program.  It checks the OEIS
terms (A335406, A089180, A006560) and, when the tool is present, compares its
own totals with the output of 'a335406 naive N -H' and 'a335406 scan N -H'.

usage: verify_a335406.py [N] [--tool PATH]        (default N = 2e8, tool ./a335406)
"""
import re
import subprocess
import sys

import numpy as np

# n: (a(n), gap index i, first prime, gap)  from A335406, A089180, A006560
KNOWN = {1: (1, 1, 2, 1), 2: (2, 2, 3, 2), 3: (49, 54, 251, 6),
         4: (633353, 654926, 9843019, 30), 5: (6706139, 6904737, 121174811, 30)}
LMAX, DMAX = 16, 4096
fails = 0


def check(ok, msg):
    global fails
    print(("ok   " if ok else "FAIL ") + msg)
    if not ok:
        fails += 1


def parse_num(s):
    s = s.replace(",", "")
    m = re.fullmatch(r"(\d+)e(\d+)", s)
    if m:
        return int(m.group(1)) * 10 ** int(m.group(2))
    m = re.fullmatch(r"(\d+)\^(\d+)", s)
    if m:
        return int(m.group(1)) ** int(m.group(2))
    return int(s)


def primes_upto(n):
    """all primes <= n, odd-only sieve of Eratosthenes"""
    sieve = np.ones(n // 2 + 1, dtype=np.bool_)          # index i <-> 2i+1
    sieve[0] = False                                      # 1 is not prime
    for i in range(1, int(n ** 0.5) // 2 + 1):
        if sieve[i]:
            p = 2 * i + 1
            sieve[p * p // 2::p] = False
    odd = 2 * np.flatnonzero(sieve) + 1
    return np.concatenate(([2], odd[odd <= n])).astype(np.int64)


def analyse(N):
    extra = 20000                                   # room to finish the last run past N
    pr = primes_upto(N + extra)
    npi = int(np.searchsorted(pr, N, side="right"))         # pi(N)
    gaps = np.diff(pr)                                       # gaps[k-1] = g_k = p_(k+1) - p_k
    same = gaps[1:] == gaps[:-1]                             # same[k-2]: g_(k-1) == g_k  (k >= 2)
    E = int(np.count_nonzero(same[: npi - 1]))               # primes p_k <= N (k >= 2) with g_(k-1) = g_k
    starts = np.concatenate(([0], np.flatnonzero(~same) + 1))  # 0-based gap indices where runs begin
    nruns = int(np.count_nonzero(starts < npi))              # runs beginning at a prime <= N
    if nruns >= len(starts):
        raise SystemExit("increase extra")
    lengths = (starts[1:] - starts[:-1])[:nruns]
    d = gaps[starts[:nruns]]
    res = {"pi": npi, "E": E, "runs": nruns, "maxlen": int(lengths.max()),
           "hist": {}, "ld": {}, "first": {}}
    hl = np.bincount(np.minimum(lengths, LMAX), minlength=LMAX + 1)
    for L in range(1, LMAX + 1):
        if hl[L]:
            res["hist"][L] = int(hl[L])
    key = np.minimum(lengths, LMAX) * DMAX + np.minimum(d, DMAX - 1)
    for k, c in zip(*np.unique(key, return_counts=True)):
        res["ld"][(int(k) // DMAX, int(k) % DMAX)] = int(c)
    for L in np.unique(lengths):
        r = int(np.flatnonzero(lengths == L)[0])            # 0-based run index
        res["first"][int(L)] = (r + 1, int(starts[r]) + 1, int(pr[starts[r]]), int(d[r]))
    return res


def run_tool(tool, args):
    out = subprocess.run([tool] + args, capture_output=True, text=True, check=True).stdout
    res = {"hist": {}, "ld": {}, "first": {}}
    for line in out.splitlines():
        m = re.match(r"primes <= [\d,]+ \(pi\)\s+([\d,]+)", line)
        if m:
            res["pi"] = parse_num(m.group(1))
        m = re.match(r"equal adjacent gaps \(E\)\s+([\d,]+)", line)
        if m:
            res["E"] = parse_num(m.group(1))
        m = re.match(r"runs \(= pi - E\)\s+([\d,]+)", line)
        if m:
            res["runs"] = parse_num(m.group(1))
        m = re.match(r"longest run\s+(\d+)", line)
        if m:
            res["maxlen"] = int(m.group(1))
        m = re.match(r"runs by length\s+(.*)", line)
        if m:
            for L, c in re.findall(r"(\d+): ([\d,]+)", m.group(1)):
                res["hist"][int(L)] = parse_num(c)
        m = re.match(r"a\((\d+)\) = ([\d,]+)\s+gap index i = ([\d,]+)\s+first prime ([\d,]+)\s+gap (\d+)", line)
        if m:
            res["first"][int(m.group(1))] = tuple(parse_num(m.group(k)) for k in range(2, 6))
        m = re.match(r"ld (\d+) (\d+) (\d+) (\d+)", line)
        if m:
            res["ld"][(int(m.group(2)), int(m.group(3)))] = int(m.group(4))
    return res


def main():
    N, tool = 200000000, "./a335406"
    args = sys.argv[1:]
    while args:
        a = args.pop(0)
        if a == "--tool":
            tool = args.pop(0)
        else:
            N = parse_num(a)

    mine = analyse(N)
    print(f"python: pi({N}) = {mine['pi']}, E = {mine['E']}, runs = {mine['runs']}, "
          f"longest run {mine['maxlen']}, runs by length {mine['hist']}")
    check(mine["runs"] == mine["pi"] - mine["E"], "runs = pi - E")
    for n, known in KNOWN.items():
        got = mine["first"].get(n)
        check(got == known, f"a({n}) = {got} (OEIS: a, gap index, first prime, gap = {known})")
    for L in range(1, mine["maxlen"] + 1):
        r, i, p, d = mine["first"][L]
        print(f"      first run of length {L}: run {r}, gap index {i}, primes {p} + {d}k")
    print("      runs of length >= 3 by gap:",
          {k: v for k, v in sorted(mine["ld"].items()) if k[0] >= 3})

    for cmd in (["naive", str(N), "-H"], ["scan", str(N), "-q", "-H"]):
        try:
            tr = run_tool(tool, cmd)
        except (OSError, subprocess.CalledProcessError) as e:
            check(False, f"could not run {tool} {' '.join(cmd)}: {e}")
            continue
        what = f"{tool} {' '.join(cmd)}"
        for k in ("pi", "E", "runs", "maxlen", "hist", "first", "ld"):
            check(tr.get(k) == mine[k], f"{what}: {k} agrees" if tr.get(k) == mine[k]
                  else f"{what}: {k} differs: tool {tr.get(k)} vs python {mine[k]}")

    print("verification " + ("FAILED" if fails else "passed"))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
