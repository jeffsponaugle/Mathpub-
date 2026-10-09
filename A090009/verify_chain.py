#!/usr/bin/env python3
"""verify_chain.py -- independent check of a digit-sum prime chain.

Usage:  python3 verify_chain.py P [P ...]

For each starting number P it follows the chain P -> P + digitsum(P) -> ...
with real digit sums, and for every term it
  * runs Miller-Rabin with the 13 smallest prime bases, which is a proof of
    primality for every n below 3,317,044,064,679,887,385,961,981
    (Sorenson and Webster, 2015), and
  * builds a Pratt certificate from scratch: n-1 is factored completely
    (trial division, then Pollard-Brent) and a witness a is found with
    a^(n-1) = 1 and a^((n-1)/q) != 1 (mod n) for every prime q | n-1, each q
    proven prime the same way, down to trial division below 10^6.
It then checks that P is a true chain start: no prime r has r + digitsum(r) = P.

Pure Python, no dependencies, and no code shared with dschain or dspattern.
"""
import math
import random
import sys

BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
PSI13 = 3317044064679887385961981


def digitsum(n):
    return sum(map(int, str(n)))


def miller_rabin(n):
    """Deterministic below PSI13."""
    if n < 2:
        return False
    for p in BASES:
        if n % p == 0:
            return n == p
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


def pollard_brent(n):
    if n % 2 == 0:
        return 2
    while True:
        y, c, m = random.randrange(1, n), random.randrange(1, n), 128
        g = r = q = 1
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += m
            r *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
        if g != n:
            return g


def factor(n, out):
    """Append the prime factors of n (with multiplicity) to out."""
    if n == 1:
        return
    for p in range(2, 10000):
        if n % p == 0:
            out.append(p)
            factor(n // p, out)
            return
        if p * p > n:
            out.append(n)
            return
    if miller_rabin(n):
        out.append(n)
        return
    d = pollard_brent(n)
    factor(d, out)
    factor(n // d, out)


def pratt(n, lines=None):
    """True if a Pratt (Lucas n-1) certificate for n is found and checks out.
    lines, if given, collects the top-level certificate for printing."""
    if n < 10**6:
        return n > 1 and all(n % d for d in range(2, math.isqrt(n) + 1))
    fs = []
    factor(n - 1, fs)
    if math.prod(fs) != n - 1:
        return False
    qs = sorted(set(fs))
    if not all(pratt(q) for q in qs):
        return False
    for a in range(2, 10000):
        if pow(a, n - 1, n) != 1:
            return False                    # Fermat witness: composite
        if all(pow(a, (n - 1) // q, n) != 1 for q in qs):
            if lines is not None:
                lines.append("n-1 = " + " * ".join(map(str, fs)) + f", witness {a}")
            return True
    return False


def verify(p):
    print(f"chain from {p}:")
    if p >= PSI13:
        print("  (above 3.3e24: the Miller-Rabin result is not a proof; "
              "the Pratt certificate still is)")
    chain, n = [], p
    ok = True
    while True:
        cert = []
        mr, pr = miller_rabin(n), pratt(n, cert)
        if not (mr and pr):
            if mr != pr:
                print(f"  {n}: DISAGREEMENT, Miller-Rabin {mr}, Pratt {pr}")
                ok = False
            break
        chain.append(n)
        d = digitsum(n)
        print(f"  {len(chain):2d}. {n}  prime (Miller-Rabin and Pratt certificate: "
              f"{cert[0] if cert else 'trial division'})")
        print(f"      digit sum {d}" + ("  -- odd, so the next number is even" if d % 2 else ""))
        if d % 2:
            break
        n += d
    if chain and n != chain[-1]:
        print(f"      next {n} is composite")
    preds = [p - d for d in range(1, 9 * len(str(p)) + 1)
             if p - d > 1 and digitsum(p - d) == d and miller_rabin(p - d) and pratt(p - d)]
    print(f"  length {len(chain)}; "
          + ("a true chain start (no prime r has r + digitsum(r) = "
             f"{p})" if not preds else f"NOT a start: preceded by {preds}"))
    return ok


def main():
    if len(sys.argv) < 2:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        sys.exit(2)
    ok = True
    for a in sys.argv[1:]:
        ok &= verify(int(a.replace("_", "").replace(",", "")))
        print()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
