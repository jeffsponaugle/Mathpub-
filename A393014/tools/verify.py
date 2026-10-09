#!/usr/bin/env python3
"""Independent verification of dbpal results (pure Python, no GMP).

  verify.py RESULTS.tsv            check every line: palindromic in both bases, primality flag
  verify.py --prime N...           BPSW probable-prime test of the given numbers
  verify.py --pal B1,B2 N...       show base representations / palindromicity

Primality: Baillie-PSW (strong base-2 Miller-Rabin + strong Lucas, Selfridge parameters),
plus Miller-Rabin to 20 fixed prime bases.  No BPSW pseudoprime is known.
"""
import sys
from math import isqrt


def digits(n, b):
    if n == 0:
        return [0]
    d = []
    while n:
        n, r = divmod(n, b)
        d.append(r)
    return d[::-1]


def is_pal(n, b):
    d = digits(n, b)
    return d == d[::-1]


def to_str(n, b):
    al = "0123456789abcdefghijklmnopqrstuvwxyz"
    return "".join(al[x] for x in digits(n, b))


def _mr(n, a):
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    x = pow(a, d, n)
    if x in (1, n - 1):
        return True
    for _ in range(s - 1):
        x = x * x % n
        if x == n - 1:
            return True
    return False


def _jacobi(a, n):
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def _strong_lucas(n):
    D = 5
    while True:
        j = _jacobi(D, n)
        if j == -1:
            break
        if j == 0 and abs(D) != n:
            return False
        D = -D - 2 if D > 0 else -D + 2
        if D == 17 and isqrt(n) ** 2 == n:
            return False
    P, Q = 1, (1 - D) // 4
    d, s = n + 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    # Lucas sequence U_d, V_d via binary method
    U, V, Qk = 1, P, Q % n
    inv2 = (n + 1) // 2
    for bit in bin(d)[3:]:
        U, V = U * V % n, (V * V - 2 * Qk) % n
        Qk = Qk * Qk % n
        if bit == "1":
            U, V = (P * U + V) * inv2 % n, (D * U + P * V) * inv2 % n
            Qk = Qk * Q % n
    if U == 0 or V == 0:
        return True
    for _ in range(s - 1):
        V = (V * V - 2 * Qk) % n
        Qk = Qk * Qk % n
        if V == 0:
            return True
    return False


SMALL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73]


def is_prime(n):
    if n < 2:
        return False
    for p in SMALL:
        if n % p == 0:
            return n == p
    if not _mr(n, 2) or not _strong_lucas(n):
        return False
    return all(_mr(n, a) for a in SMALL[1:])


def verify_file(path):
    bad = n_ok = 0
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or not p[0].isdigit():
                continue
            n, flag, b1, b2 = int(p[0]), p[1] == "1", int(p[2]), int(p[3])
            ok = is_pal(n, b1) and is_pal(n, b2) and is_prime(n) == flag
            if len(p) >= 6:
                ok = ok and p[4] == to_str(n, b1) and p[5] == to_str(n, b2)
            if ok:
                n_ok += 1
            else:
                bad += 1
                print("FAILED:", line.strip())
    print(f"{path}: {n_ok} lines verified, {bad} failures")
    return 1 if bad else 0


def main():
    a = sys.argv[1:]
    if not a:
        print(__doc__)
        return 1
    if a[0] == "--prime":
        for s in a[1:]:
            n = int(s)
            print(n, "prime (BPSW)" if is_prime(n) else "composite")
        return 0
    if a[0] == "--pal":
        b1, b2 = map(int, a[1].split(","))
        for s in a[2:]:
            n = int(s)
            print(n, f"base{b1}={to_str(n, b1)} {'pal' if is_pal(n, b1) else 'not pal'}",
                  f"base{b2}={to_str(n, b2)} {'pal' if is_pal(n, b2) else 'not pal'}")
        return 0
    rc = 0
    for path in a:
        rc |= verify_file(path)
    return rc


if __name__ == "__main__":
    sys.exit(main())
