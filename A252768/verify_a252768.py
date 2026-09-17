#!/usr/bin/env python3
"""
verify_a252768.py P [KMAX]

Independent check of a candidate term of OEIS A252768.  Shares no code with
a252768.c: the primes up to P come from a numpy segmented sieve (not
primesieve), the gap power sums S_k = sum over gaps g of count(g) * g^k are
formed exactly with Python integers, and each S_k is tested with Miller-Rabin
(the 13 prime bases 2..41, deterministic below 3.3e24, plus 40 random bases
for larger values).  Prints S_1..S_KMAX with verdicts and the depth of P:
the largest n such that S_1..S_n are all prime.
"""
import math
import random
import sys
import time

import numpy as np


def small_primes(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]


def gap_histogram(P, seg_odds=1 << 25):
    """Return (histogram dict gap -> count, number of primes <= P, last prime)."""
    base = [p for p in small_primes(int(math.isqrt(P)) + 1) if p > 2]
    hist = np.zeros(1 << 16, dtype=np.int64)
    prev = 2                      # the prime before the current segment
    nprimes = 1                   # the prime 2
    lo = 3                        # odd numbers lo, lo+2, ... per segment
    t0 = time.time()
    while lo <= P:
        n = min(seg_odds, (P - lo) // 2 + 1)
        hi = lo + 2 * n           # exclusive, odd numbers in [lo, hi)
        s = np.ones(n, dtype=bool)
        for p in base:
            if p * p >= hi:
                break
            start = max(p * p, ((lo + p - 1) // p) * p)
            if start % 2 == 0:
                start += p
            if start >= hi:
                continue
            s[(start - lo) // 2::p] = False
        primes = np.nonzero(s)[0] * 2 + lo
        if len(primes):
            gaps = np.diff(np.concatenate(([prev], primes)))
            hist += np.bincount(gaps, minlength=len(hist))[: len(hist)]
            prev = int(primes[-1])
            nprimes += len(primes)
        lo = hi
        if time.time() - t0 > 10:
            print(f"  ... sieved to {lo:,} ({100.0 * lo / P:.1f}%)", file=sys.stderr, flush=True)
            t0 = time.time()
    return {int(g): int(c) for g, c in enumerate(hist) if c}, nprimes, prev


def is_probable_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
    if n >= 3317044064679887385961981:
        bases += [random.randrange(2, n - 1) for _ in range(40)]
    for a in bases:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    P = int(float(sys.argv[1])) if "e" in sys.argv[1] else int(sys.argv[1])
    kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    t0 = time.time()
    hist, nprimes, last = gap_histogram(P)
    if last != P:
        print(f"{P} is not prime (largest prime <= P is {last})")
        sys.exit(1)
    print(f"p = {P}: {nprimes} primes <= p, {sum(hist.values())} gaps, largest gap {max(hist)} ({time.time() - t0:.1f} s)")
    assert sum(g * c for g, c in hist.items()) == P - 2
    depth = 0
    for k in range(1, kmax + 1):
        S = sum(c * g ** k for g, c in hist.items())
        prime = is_probable_prime(S)
        print(f"S_{k:<2d} = {S}  ({S.bit_length()} bits)  {'prime' if prime else 'composite'}")
        if prime and depth == k - 1:
            depth = k
    print(f"depth({P}) = {depth}: S_1..S_{depth} are prime")


if __name__ == "__main__":
    main()
