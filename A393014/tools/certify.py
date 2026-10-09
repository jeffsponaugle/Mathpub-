#!/usr/bin/env python3
"""Primality certificates (Pocklington / Lucas-Pratt) for the new prime terms.

  certify.py N [N ...]      prints a certificate for each N (or says why it could not make one)

Pocklington: if N - 1 = F * R with F fully factored, F > sqrt(N), and for every prime q | F there
is a base a with a^(N-1) == 1 (mod N) and gcd(a^((N-1)/q) - 1, N) = 1, then N is prime.  The prime
factors q of F are themselves certified recursively (small ones by trial division).
"""
import math
import random
import sys
import time

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from verify import is_prime  # BPSW, used only to decide which factors to certify recursively

SMALL = 10**12  # below this, primality is checked by trial division in the certificate


def brent(n, deadline):
    if n % 2 == 0:
        return 2
    while time.time() < deadline:
        y, c, m = random.randrange(1, n), random.randrange(1, n), 512
        g = r = q = 1
        x = ys = y
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
            if time.time() > deadline:
                return None
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
        if 1 < g < n:
            return g
    return None


def factor(n, deadline, found):
    """Split n into prime factors as far as possible before the deadline; returns the unfactored rest."""
    for p in range(2, 100000):
        while n % p == 0:
            found[p] = found.get(p, 0) + 1
            n //= p
        if p * p > n:
            break
    stack, rest = [n] if n > 1 else [], 1
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            found[m] = found.get(m, 0) + 1
            continue
        d = brent(m, deadline)
        if d is None:
            rest *= m
            continue
        stack += [d, m // d]
    return rest


def certify(N, depth=0, budget=120.0):
    """Return a nested certificate dict, or None."""
    if N < SMALL:
        return {"N": N, "method": "trial division"} if all(N % p for p in range(2, math.isqrt(N) + 1)) else None
    deadline = time.time() + budget
    found = {}
    rest = factor(N - 1, deadline, found)
    F = (N - 1) // rest
    if F * F <= N:
        return None
    sub = []
    for q in sorted(found):
        a = 2
        while a < 1000:
            if pow(a, N - 1, N) == 1 and math.gcd(pow(a, (N - 1) // q, N) - 1, N) == 1:
                break
            a += 1
        else:
            return None
        cq = certify(q, depth + 1, budget / 2)
        if cq is None:
            return None
        sub.append({"q": q, "e": found[q], "a": a, "cert": cq})
    return {"N": N, "method": "Pocklington", "F_digits": len(str(F)), "R": rest, "factors": sub}


def show(c, indent=0):
    pad = "  " * indent
    if c["method"] == "trial division":
        print(f"{pad}{c['N']}: prime (trial division)")
        return
    print(f"{pad}{c['N']}: prime by Pocklington; N-1 = F*R with F > sqrt(N) fully factored "
          f"({c['F_digits']} digits), R = {c['R']}")
    for f in c["factors"]:
        print(f"{pad}  q = {f['q']}^{f['e']}: witness a = {f['a']}")
        if f["cert"]["method"] != "trial division":
            show(f["cert"], indent + 2)


def main():
    for s in sys.argv[1:]:
        N = int(s)
        t0 = time.time()
        c = certify(N)
        if c is None:
            print(f"{N}: no certificate within the time budget (N-1 not factored far enough)")
        else:
            show(c)
            print(f"   ({time.time() - t0:.1f}s)")
        print()


if __name__ == "__main__":
    main()
