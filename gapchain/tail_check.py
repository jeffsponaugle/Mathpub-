#!/usr/bin/env python3
"""Close the gap between gapsieve's clamped --end and 2^64.

Lists every prime in [2^64 - W - 1000, 2^64 + 1000] with a Miller-Rabin test
that is deterministic below 3.3e24 (bases = the first 13 primes), then checks
every start p in [2^64 - W, 2^64) for an increasing chain of >= L primes
(gaps 2, 4, ..., 2(L-1)) and a decreasing L-prime pattern (gaps 2(L-1), ...,
4, 2).  Chains that start below 2^64 may end above it; that is the point.
usage: tail_check.py [L] [W]"""
import subprocess
import sys

L = int(sys.argv[1]) if len(sys.argv) > 1 else 17
W = int(sys.argv[2]) if len(sys.argv) > 2 else 10**6
T = 2**64
BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)


def is_prime(n):
    if n < 2:
        return False
    for b in BASES:
        if n % b == 0:
            return n == b
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in BASES:
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


lo, hi = T - W - 1000, T + 1000
# trial-sieve the window by small primes, then test the survivors
small = [q for q in range(3, 20000, 2) if all(q % r for r in range(3, int(q**0.5) + 1, 2))]
comp = bytearray(hi - lo + 1)
for q in small:
    for x in range((-lo) % q, hi - lo + 1, q):
        comp[x] = 1
primes = [lo + i for i in range(0, hi - lo + 1)
          if not comp[i] and (lo + i) % 2 and is_prime(lo + i)]
print(f"{len(primes)} primes in [2^64 - {T - lo}, 2^64 + {hi - T}]")
below = [p for p in primes if p < T]
print(f"largest below 2^64: 2^64 - {T - below[-1]}; "
      f"smallest above: 2^64 + {primes[len(below)] - T}")

# primesieve lists the part below 2^64 exactly: the two lists must agree
ps = subprocess.run(["primesieve", str(lo), str(T - 1), "-p"],
                    capture_output=True, text=True, check=True).stdout
assert [int(x) for x in ps.split()] == below, "Miller-Rabin disagrees with primesieve"
print(f"primesieve agrees on all {len(below)} primes below 2^64")

best_inc = best_dec = (0, None)
n = len(primes)
for i, p in enumerate(primes):
    if not T - W <= p < T:
        continue
    m = 1                       # increasing chain from p
    while i + m < n and primes[i + m] - primes[i + m - 1] == 2 * m:
        m += 1
    assert i + m < n
    best_inc = max(best_inc, (m, p))
    # longest decreasing run p, p+2(L-1), ... matching the L-prime pattern
    j = 0
    while j < L - 1 and primes[i + j + 1] - primes[i + j] == 2 * (L - 1 - j):
        j += 1
    best_dec = max(best_dec, (j + 1, p))
print(f"starts in [2^64 - {W}, 2^64): longest increasing chain {best_inc[0]} "
      f"primes (at 2^64 - {T - best_inc[1]}); longest prefix of the "
      f"{L}-prime decreasing pattern {best_dec[0]} (at 2^64 - {T - best_dec[1]})")
print("no chain of length", L, "of either shape" if best_inc[0] < L and best_dec[0] < L
      else "FOUND A CHAIN")
