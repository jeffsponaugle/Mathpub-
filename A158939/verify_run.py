#!/usr/bin/env python3
"""Independent check of an A158939 term: verify_run.py P [N]

Prints the prime gaps following the prime P and its run length L(P), the
largest n with g_1 < g_2 < ... < g_n.  Uses its own deterministic Miller-Rabin
(valid below 2^64 with the bases 2, 325, 9375, 28178, 450775, 9780504,
1795265022) and a plain odd-number next-prime search, so it shares no code
with a158939.c or primesieve.  With N given, exits 1 unless L(P) == N.
"""
import sys

BASES = (2, 325, 9375, 28178, 450775, 9780504, 1795265022)
SMALL = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n):
    if n < 2:
        return False
    for p in SMALL:
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in BASES:
        a %= n
        if a == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def next_prime(n):
    m = n + 1
    if m <= 2:
        return 2
    if m % 2 == 0:
        m += 1
    while not is_prime(m):
        m += 2
    return m


def run(p):
    """Return (gaps, L): the L increasing gaps after p followed by the gap that ends the run."""
    gaps, prev, prevgap = [], p, 0
    while True:
        q = next_prime(prev)
        g = q - prev
        gaps.append(g)
        if g <= prevgap:
            return gaps, len(gaps) - 1
        prev, prevgap = q, g


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    p = int(eval(argv[1].replace('^', '**'), {"__builtins__": {}}, {}))
    if not is_prime(p):
        print(f"{p} is not prime")
        return 1
    gaps, L = run(p)
    print(f"{p}: gaps {' '.join(map(str, gaps[:-1]))} | {gaps[-1]}  ->  L(p) = {L}")
    if len(argv) > 2:
        n = int(argv[2])
        print("OK: run length is", n) if L == n else print(f"MISMATCH: expected {n}, got {L}")
        return 0 if L == n else 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
