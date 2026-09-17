#!/usr/bin/env python3
"""
Near-tail model for A053686 from the histogram of an a053686 scan.

Reads the 'hist N GAP COUNT' lines that 'a053686 scan ... -H' prints (gaps of
size GAP after primes strictly inside record interval N, GAP within 198 of
the record g_N), and estimates how many repeats of g_N one would expect in
each interval from the gaps just below the record:

  * The count of gaps of size g in an interval is modelled as
        c(g) = c'(g) * A_N * exp(b * (g_N - g)),
    i.e. exponential in the deficit d = g_N - g (smaller gaps are more common;
    b ~ 1/20 means the counts double every 14 below the record), times the
    Hardy-Littlewood factor
        c'(g) = prod_{p | g, p > 2} (p-1)/(p-2),
    which makes gaps divisible by 3 about twice as common as their neighbours
    (and by 5, 7 ... a little more).  The slope b is fitted once on the
    normalised counts c(g)/c'(g) summed over all decided intervals for
    deficits 2..DMAX (default 100), each interval then gets its own amplitude
    A_N from its normalised count over the same deficits, and the expected
    number of repeats is
        E_N = c'(g_N) * A_N.
  * Because the record g_N is by construction the largest gap up to P_(N+1),
    the extrapolation to d = 0 is where the model is least certain; the sum of
    E_N over intervals with a known outcome, against the observed number of
    repeats, calibrates it.

usage: tail_a053686.py SCANFILE [--from N] [--to N] [--dmax D]
"""
import math
import re
import sys

REC = {}        # n -> record gap, from the hist lines' own record intervals (gap of the k = 0 entry)
KNOWN_TERMS = {2, 4, 6, 14, 34, 36, 52, 86, 132, 154, 250, 336}
RECORD_GAPS = [1, 2, 4, 6, 8, 14, 18, 20, 22, 34, 36, 44, 52, 72, 86, 96, 112, 114, 118, 132, 148, 154, 180, 210,
               220, 222, 234, 248, 250, 282, 288, 292, 320, 336, 354, 382, 384, 394, 456, 464, 468, 474, 486, 490,
               500, 514, 516, 532, 534, 540, 582, 588, 602, 652, 674, 716, 766, 778, 804, 806, 906, 916, 924, 1132,
               1184, 1198, 1220, 1224, 1248, 1272, 1328, 1356, 1370, 1442, 1476, 1488, 1510, 1526, 1530, 1550]


def hl(g):
    """prod_{p | g, p > 2} (p-1)/(p-2)"""
    if g < 2:
        return 1.0
    f, m, p = 1.0, g, 3
    while m % 2 == 0:
        m //= 2
    while p * p <= m:
        if m % p == 0:
            f *= (p - 1) / (p - 2)
            while m % p == 0:
                m //= p
        p += 2
    if m > 1:
        f *= (m - 1) / (m - 2)
    return f


def main():
    args = sys.argv[1:]
    lo, hi = 1, 10 ** 9
    if "--from" in args:
        k = args.index("--from"); lo = int(args[k + 1]); del args[k:k + 2]
    if "--to" in args:
        k = args.index("--to"); hi = int(args[k + 1]); del args[k:k + 2]
    dmax = 100
    if "--dmax" in args:
        k = args.index("--dmax"); dmax = int(args[k + 1]); del args[k:k + 2]
    if len(args) != 1:
        print(__doc__); sys.exit(2)
    hist = {}                                   # (n, gap) -> count
    decided = set()
    for line in open(args[0]):
        m = re.match(r"hist (\d+) (\d+) (\d+)$", line)
        if m:
            hist[(int(m.group(1)), int(m.group(2)))] = int(m.group(3))
        m = re.match(r"DECIDED\s+#(\d+)", line)
        if m:
            decided.add(int(m.group(1)))
        m = re.match(r"\s*(\d+)\s+\d+\s+\d+\s+\S+\s+\d+\s+.*\s(TERM|not a term),", line)
        if m:
            decided.add(int(m.group(1)))
    ns = sorted({n for n, g in hist if lo <= n <= hi and n in decided and n >= 10})
    if not ns:
        print("no decided record intervals with histogram data in the file"); sys.exit(1)

    # normalised counts by deficit, summed over intervals
    agg = {}
    for n in ns:
        G = RECORD_GAPS[n - 1]
        for d in range(2, min(dmax, G - 2) + 1, 2):
            c = hist.get((n, G - d), 0)
            agg[d] = agg.get(d, 0.0) + c / hl(G - d)
    # least squares on log counts, deficits with at least 5 gaps
    pts = [(d, math.log(v)) for d, v in agg.items() if v >= 5]
    if len(pts) < 5:
        print("too few gaps near the records to fit a slope"); sys.exit(1)
    mx = sum(d for d, _ in pts) / len(pts)
    my = sum(y for _, y in pts) / len(pts)
    b = sum((d - mx) * (y - my) for d, y in pts) / sum((d - mx) ** 2 for d, _ in pts)
    print(f"records {ns[0]}..{ns[-1]} ({len(ns)} decided intervals), deficits 2..{dmax}: fitted tail slope b = {b:.4f} "
          f"per unit of gap (counts double every {math.log(2) / b:.1f} below the record; 1/b = {1 / b:.1f})")
    # per-interval amplitude from its normalised count over deficits 2..dmax
    print(f"\n  n   gap   gaps within {dmax:3d}   HL(g)  expected repeats  observed")
    tot_e, tot_o = 0.0, 0
    for n in ns:
        G = RECORD_GAPS[n - 1]
        dm = min(dmax, G - 2)
        T = sum(hist.get((n, G - d), 0) / hl(G - d) for d in range(2, dm + 1, 2))
        within = sum(hist.get((n, G - d), 0) for d in range(2, dm + 1, 2))
        A = T / sum(math.exp(b * d) for d in range(2, dm + 1, 2))
        E = hl(G) * A
        obs = hist.get((n, G), 0)
        tot_e += E
        tot_o += obs
        tag = ""
        if n <= 34:
            tag = "  (OEIS term)" if G in KNOWN_TERMS else ""
        print(f" {n:3d}  {G:5d}  {within:14d}   {hl(G):5.2f}  {E:16.3f}  {obs:8d}{tag}")
    p_none = math.exp(-tot_e)
    print(f"\n  total expected repeats {tot_e:.2f}, observed {tot_o}; Poisson P(0 | {tot_e:.2f}) = {p_none:.3f}")
    if tot_o:
        print(f"  observed / expected = {tot_o / tot_e:.2f}")


if __name__ == "__main__":
    main()
