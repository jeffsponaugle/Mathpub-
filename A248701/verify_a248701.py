#!/usr/bin/env python3
"""Independent check of a248701.c (OEIS A248701-A248704).

Same definitions, different implementation: a segmented numpy sieve, run
lengths computed with vectorised prefix operations, and for single windows a
pure-Python deterministic Miller-Rabin test.  Nothing here is shared with the
C program apart from the definitions:

  the prime P sits between the gap into it and the gap out of it;
  peak depth          = min(non-decreasing run into P, non-increasing run out of P)
  valley depth        = min(non-increasing run into P, non-decreasing run out of P)
  strict peak depth   = min(strictly increasing run in,  strictly decreasing run out)
  strict valley depth = min(strictly decreasing run in,  strictly increasing run out)
  A248701(n) = smallest P with peak depth >= n            (offset 1)
  A248702(n) = smallest P with valley depth >= n          (offset 0, a(0) = 2)
  A248703(n) = smallest P with strict peak depth >= n+1   (offset 1)
  A248704(n) = smallest P with strict valley depth >= n   (offset 1)

usage: verify_a248701.py scan LIMIT [SEGMENT] [--procs=N]   exhaustive scan of the primes below LIMIT
       verify_a248701.py window P [N]              the N primes on each side of P, the gaps,
                                                   and P's depth in the four modes
Numbers may be written as 1e12 or 10**12.
"""
import sys
import time

import numpy as np

MODES = ("peak", "valley", "strict peak", "strict valley")
SEQ = ("A248701", "A248702", "A248703", "A248704")
SHIFT = (0, 0, 1, 0)
FIRSTN = (1, 0, 1, 1)
KNOWN = {
    "A248701": [3, 7, 359, 7853, 96401, 2812099, 294276293],
    "A248702": [2, 3, 19, 43, 2687, 179819, 1107791],
    "A248703": [23, 1439, 21433, 1130863, 19881311, 331542583],
    "A248704": [3, 19, 1429, 25243, 340577, 1107791],
}
DMAX = 16


def num(s):
    return int(eval(s.replace("^", "**"), {}, {})) if any(c in s for c in "e*^") else int(s)


# ---------------------------------------------------------------- primality (pure Python)

_SMALL = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n):
    """Deterministic for n < 3.3e24 (Miller-Rabin to the first 12 prime bases)."""
    if n < 2:
        return False
    for p in _SMALL:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _SMALL:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prev_prime(n):
    n -= 1
    while not is_prime(n):
        n -= 1
    return n


def next_prime(n):
    n += 1
    while not is_prime(n):
        n += 1
    return n


def runs_from_gaps(gaps_in, gaps_out):
    """gaps_in[0] is the gap into P (then going backwards), gaps_out[0] the gap out of P."""
    def run(seq, cont):
        r = 1 if seq else 0
        while r < len(seq) and cont(seq[r], seq[r - 1]):
            r += 1
        return r
    # going backwards from P, "earlier" is seq[r], "later" is seq[r-1]
    L = (run(gaps_in, lambda earlier, later: earlier <= later),   # non-decreasing into P
         run(gaps_in, lambda earlier, later: earlier >= later),   # non-increasing into P
         run(gaps_in, lambda earlier, later: earlier < later),    # strictly increasing into P
         run(gaps_in, lambda earlier, later: earlier > later))    # strictly decreasing into P
    M = (run(gaps_out, lambda nxt, cur: nxt <= cur),              # non-increasing out of P
         run(gaps_out, lambda nxt, cur: nxt >= cur),              # non-decreasing out of P
         run(gaps_out, lambda nxt, cur: nxt < cur),               # strictly decreasing out of P
         run(gaps_out, lambda nxt, cur: nxt > cur))               # strictly increasing out of P
    return [min(a, b) for a, b in zip(L, M)]


def cmd_window(P, n):
    if not is_prime(P):
        sys.exit(f"{P} is not prime")
    W = max(n, 64)
    before, after = [], []
    q = P
    while len(before) < W and q > 2:
        q = prev_prime(q)
        before.append(q)                         # before[0] is the prime just below P
    q = P
    while len(after) < W:
        q = next_prime(q)
        after.append(q)
    gaps_in = [P - before[0]] + [before[i] - before[i + 1] for i in range(len(before) - 1)] if before else []
    gaps_out = [after[0] - P] + [after[i + 1] - after[i] for i in range(len(after) - 1)]
    depth = runs_from_gaps(gaps_in, gaps_out)
    print(f"{P}: " + ", ".join(f"{MODES[m]} depth {depth[m]} ({SEQ[m]} index {depth[m] - SHIFT[m]})" for m in range(4)))
    nb = min(n, len(before))
    primes = list(reversed(before[:nb])) + [P] + after[:n]
    # every number between consecutive listed primes must be composite (the list is consecutive primes)
    for a, b in zip(primes, primes[1:]):
        assert all(not is_prime(x) for x in range(a + 1, b)), (a, b)
    print("  primes:", " ".join(str(p) for p in primes))
    print("  gaps:  ", " ".join(str(g) for g in reversed(gaps_in[:nb])), "|", " ".join(str(g) for g in gaps_out[:n]))
    print("  (checked with Miller-Rabin: all listed numbers are prime and the numbers between them are composite)")


# ---------------------------------------------------------------- segmented scan (numpy)

def simple_sieve(limit):
    s = np.ones(limit // 2 + 1, dtype=np.bool_)
    s[0] = False
    for i in range(1, int(limit ** 0.5) // 2 + 1):
        if s[i]:
            p = 2 * i + 1
            s[p * p // 2::p] = False
    return 2 * np.nonzero(s)[0] + 1            # odd primes <= limit


def segment_primes(lo, hi, base):
    """Primes in [lo, hi) using the odd base primes (all base primes with p*p < hi)."""
    lo_odd = lo | 1
    size = (hi - lo_odd + 1) // 2
    if size <= 0:
        return np.array([2], dtype=np.int64) if lo <= 2 < hi else np.zeros(0, dtype=np.int64)
    s = np.ones(size, dtype=np.bool_)
    for p in base[base.astype(np.float64) ** 2 < hi]:
        p = int(p)
        m = max(p * p, ((lo_odd + p - 1) // p) * p)
        if m % 2 == 0:
            m += p
        if m < hi:
            s[(m - lo_odd) // 2::p] = False
    if lo_odd == 1:
        s[0] = False
    out = lo_odd + 2 * np.nonzero(s)[0]
    if lo <= 2 < hi:
        out = np.concatenate(([2], out))
    return out.astype(np.int64)


def back_run(cond, j, m):
    reset = np.ones(m, dtype=np.bool_)
    reset[1:] = ~cond
    return j - np.maximum.accumulate(np.where(reset, j, 0)) + 1


def fwd_run(cond, j, m):
    reset = np.ones(m, dtype=np.bool_)
    reset[:-1] = ~cond
    return np.minimum.accumulate(np.where(reset, j, m - 1)[::-1])[::-1] - j + 1


_BASE = None


def _init_worker(limit):
    global _BASE
    _BASE = simple_sieve(int((limit + 40000) ** 0.5) + 1)


def _scan_segment(bounds):
    """Records and depth counts of the centres P in [lo, hi); independent of every other segment."""
    lo, hi = bounds
    margin = 20000
    primes = segment_primes(max(lo - margin, 0), hi + margin, _BASE)
    g = np.diff(primes)
    m = len(g)
    j = np.arange(m)
    ge, le, gt, lt = g[1:] >= g[:-1], g[1:] <= g[:-1], g[1:] > g[:-1], g[1:] < g[:-1]
    k = np.arange(1, m)
    P = primes[k]
    sel = (P >= lo) & (P < hi)
    depths = (np.minimum(back_run(ge, j, m)[k - 1], fwd_run(le, j, m)[k]),
              np.minimum(back_run(le, j, m)[k - 1], fwd_run(ge, j, m)[k]),
              np.minimum(back_run(gt, j, m)[k - 1], fwd_run(lt, j, m)[k]),
              np.minimum(back_run(lt, j, m)[k - 1], fwd_run(gt, j, m)[k]))
    first = [[0] * (DMAX + 1) for _ in range(4)]
    count = [[0] * (DMAX + 1) for _ in range(4)]
    for mode in range(4):
        d = depths[mode]
        for t in range(1, DMAX + 1):
            w = np.nonzero((d >= t) & sel)[0]
            if len(w) == 0:
                break
            count[mode][t] = len(w)
            first[mode][t] = int(P[w[0]])
    inseg = (primes >= lo) & (primes < hi)
    nprimes = int(np.count_nonzero(inseg))
    gi = np.nonzero(inseg[:-1])[0]
    maxgap = int(g[gi].max()) if len(gi) else 0
    return first, count, nprimes, maxgap


def cmd_scan(limit, seg, procs):
    t0 = time.time()
    first = [[0] * (DMAX + 1) for _ in range(4)]
    count = [[0] * (DMAX + 1) for _ in range(4)]
    nprimes = 0
    maxgap = 0
    segments = [(lo, min(lo + seg, limit)) for lo in range(0, limit, seg)]
    done = 0

    def merge(res):
        nonlocal nprimes, maxgap, done
        f, c, n, mg = res
        for mode in range(4):
            for t in range(1, DMAX + 1):
                count[mode][t] += c[mode][t]
                if f[mode][t] and (not first[mode][t] or f[mode][t] < first[mode][t]):
                    first[mode][t] = f[mode][t]
        nprimes += n
        maxgap = max(maxgap, mg)
        done += 1
        if sys.stderr.isatty() or done == len(segments) or done % 64 == 0:
            el = time.time() - t0
            print(f"\r{done}/{len(segments)} segments ({100 * done / len(segments):5.1f}%) {el:.0f} s",
                  end="" if sys.stderr.isatty() else "\n", file=sys.stderr, flush=True)

    if procs <= 1:
        _init_worker(limit)
        for b in segments:
            merge(_scan_segment(b))
    else:
        import multiprocessing as mp
        with mp.Pool(procs, initializer=_init_worker, initargs=(limit,)) as pool:
            for res in pool.imap_unordered(_scan_segment, segments):
                merge(res)
    print(file=sys.stderr)
    print(f"# scanned [0, {limit}) in {time.time() - t0:.1f} s ({procs} process{'es' if procs > 1 else ''}, segment {seg}): {nprimes} primes, largest gap {maxgap}")
    print("# centres of depth >= d:")
    for mode in range(4):
        print(f"#   {MODES[mode]:14s}", " ".join(f"{count[mode][t]:12d}" for t in range(1, DMAX + 1) if count[mode][t]))
    for mode in range(4):
        terms = []
        if mode == 1:
            terms.append("a(0) = 2")
        for t in range(1, DMAX + 1):
            n = t - SHIFT[mode]
            if n < FIRSTN[mode]:
                continue
            if first[mode][t]:
                known = KNOWN[SEQ[mode]]
                idx = n - FIRSTN[mode]
                tag = ""
                if idx < len(known):
                    tag = "  (OEIS)" if known[idx] == first[mode][t] else f"  ** OEIS has {known[idx]} **"
                terms.append(f"a({n}) = {first[mode][t]}{tag}")
            else:
                terms.append(f"a({n}) > {limit - 1}")
                break
        print(f"{SEQ[mode]} ({MODES[mode]}):")
        for t in terms:
            print("  " + t)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--procs")]
    procs = 1
    for a in sys.argv[1:]:
        if a.startswith("--procs="):
            procs = int(a.split("=", 1)[1])
    if len(args) >= 2 and args[0] == "scan":
        cmd_scan(num(args[1]), num(args[2]) if len(args) > 2 else 1 << 27, procs)
    elif len(args) >= 2 and args[0] == "window":
        cmd_window(num(args[1]), int(args[2]) if len(args) > 2 else 8)
    else:
        sys.exit(__doc__)
