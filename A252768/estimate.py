#!/usr/bin/env python3
"""
estimate.py [SCANNED]

Where should the next terms of A252768 be?  A heuristic model of the density
of primes of depth >= n, calibrated on the exhaustive scan below 10^12:

  * a prime p is the upper twin (S_1 = p - 2 prime) with probability
    2 C_2 / ln p  (C_2 = 0.66016, the twin prime constant);
  * given S_1..S_(k-1) prime, S_k is prime with probability B_k * 2 / ln S_k,
    where B_k = prod_{3 <= q <= k, q prime} q/(q-1) is the congruence boost
    (g^k mod q depends only on k mod q-1, so S_k is coprime to every prime
    q <= k once the earlier sums are prime);
  * S_(k+1)/S_k ~ 0.8 (k+1) ln p, fitted to the actual sums at a(9).

The model reproduces the measured per-k pass rates below 10^12 within 1-3 %.
SCANNED (default 1e12) is the bound already searched exhaustively; the table
gives the expected number of new primes of depth >= n in [SCANNED, X] and
P(a(n) < X) = 1 - exp(-E) for a Poisson count.
"""
import math
import sys

MEASURED_EXACT = {1: 1747205024, 2: 112645225, 3: 9907858, 4: 756501, 5: 65240,
                  6: 4934, 7: 401, 8: 32, 9: 5}          # exact depth counts below 1e12


def boost(k):
    b = 1.0
    for q in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        if q <= k:
            b *= q / (q - 1)
    return b


def ln_S(k, lnp):
    return lnp + (k - 1) * (math.log(0.8) + math.log(lnp)) + math.lgamma(k + 1)


def pass_rate(k, lnp):
    return min(1.0, boost(k) * 2.0 / ln_S(k, lnp))


def density(n, p):
    lnp = math.log(p)
    q = 1.3203 / lnp
    for k in range(2, n + 1):
        q *= pass_rate(k, lnp)
    return q / lnp


def expected(n, a, b, steps=4000):
    la, lb = math.log(max(a, 1e6)), math.log(b)
    h = (lb - la) / steps
    tot = 0.0
    for i in range(steps + 1):
        x = math.exp(la + i * h)
        tot += (0.5 if i in (0, steps) else 1.0) * density(n, x) * x * h
    return tot


def quantile(n, scanned, prob):
    lo, hi = math.log10(scanned), 22.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if 1 - math.exp(-expected(n, scanned, 10 ** mid)) < prob:
            lo = mid
        else:
            hi = mid
    return 10 ** ((lo + hi) / 2)


def main():
    scanned = float(sys.argv[1]) if len(sys.argv) > 1 else 1e12
    cum = {n: sum(v for m, v in MEASURED_EXACT.items() if m >= n) for n in range(1, 10)}
    print("model check: per-k pass rate below 1e12, measured vs model (evaluated at p = 3e11)")
    for k in range(2, 10):
        print(f"  k={k:2d}  measured {cum[k] / cum[k - 1]:.4f}   model {pass_rate(k, math.log(3e11)):.4f}")
    print("model check: primes of depth >= n below 1e12, model vs measured")
    for n in range(4, 10):
        print(f"  n={n:2d}  model {expected(n, 1e6, 1e12):9.1f}   measured {cum[n]}")
    print(f"\nexpected new primes of depth >= n in [{scanned:.0e}, X] and P(a(n) < X):")
    print("  X       " + "".join(f"        n={n:<10d}" for n in (10, 11, 12, 13)))
    for X in (1e13, 1e14, 1e15, 1e16, 1e17, 1e18):
        if X <= scanned:
            continue
        row = f"  {X:8.0e}"
        for n in (10, 11, 12, 13):
            E = expected(n, scanned, X)
            row += f"   E={E:7.2f} P={1 - math.exp(-E):6.1%}"
        print(row)
    print("\nquantiles of a(n) given a(n) > SCANNED:")
    for n in (10, 11, 12, 13):
        print(f"  a({n}): 50% by {quantile(n, scanned, 0.5):.1e}, 90% by {quantile(n, scanned, 0.9):.1e}")
    print("\nbits of S_k (128 is the tool's limit):")
    for k in (11, 12, 13):
        print(f"  k={k}: " + "  ".join(f"{p:.0e}: {ln_S(k, math.log(p)) / math.log(2):5.1f}" for p in (1e13, 1e14, 3e14, 1e15, 1e16)))


if __name__ == "__main__":
    main()
