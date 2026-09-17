#!/usr/bin/env python3
"""
verify_run.py L R

Independent check (pure Python, no third-party modules) that L and R are
consecutive members of A073491 (numbers with no prime gap in their
factorization), i.e. that L+1 .. R-1 is a maximal run of exactly R-L-1
consecutive numbers with at least one prime gap (A073492).  If it is the
first such run, L+1 = A137723(R-L-1).

The verdict never depends on a complete factorization.  A number m is split
into its prime factors below TRIAL (found by trial division) and a cofactor C
whose prime factors all exceed TRIAL.  Then:
  * if small primes were found, m is gap-free iff they are consecutive and the
    chain of following primes divides C down to 1 (a few divisions);
  * otherwise m is gap-free iff C is prime, or C is a product of consecutive
    primes p1 < ... < ps with t factors counted with multiplicity: then
    p1 < C^(1/t) < ps, so p1 is one of the t-1 largest primes below the t-th
    root of C, and each candidate is checked by chain division.
A full factorization (Pollard-Brent rho with an iteration budget) is only
attempted for display.

Example:  python3 verify_run.py 18014398509481951 18014398509481984
"""
import sys
from math import isqrt, gcd

TRIAL = 100000
SMALL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

_sieve = bytearray([1]) * (TRIAL + 1)
_sieve[0] = _sieve[1] = 0
for _i in range(2, isqrt(TRIAL) + 1):
    if _sieve[_i]:
        _sieve[_i * _i::_i] = bytearray(len(_sieve[_i * _i::_i]))
TRIAL_PRIMES = [i for i in range(2, TRIAL + 1) if _sieve[i]]


def is_prime(n):
    """Miller-Rabin: deterministic (first 12 prime bases) for n < 3.3e24,
    strong probable prime to 25 prime bases above that."""
    if n < 2:
        return False
    for p in SMALL:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in SMALL:
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


def next_prime(p):
    q = p + 1
    while not is_prime(q):
        q += 1
    return q


def prev_prime(p):
    q = p - 1
    while q >= 2 and not is_prime(q):
        q -= 1
    return q


def iroot(n, t):
    """floor(n ** (1/t))"""
    if t == 1:
        return n
    if t == 2:
        return isqrt(n)
    lo, hi = 1, 1 << (n.bit_length() // t + 1)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** t <= n:
            lo = mid
        else:
            hi = mid
    return lo


def chain(x, p):
    """Does x equal p^e * nextprime(p)^f * ... exactly?  Returns the factor list or None."""
    out = []
    while x > 1:
        if x % p:
            return None
        e = 0
        while x % p == 0:
            x //= p
            e += 1
        out.append((p, e))
        p = next_prime(p)
    return out


def split_small(m):
    """(small prime factors [(p, e)], cofactor C with all primes > TRIAL)"""
    f = []
    exhausted = True
    for p in TRIAL_PRIMES:
        if p * p > m:
            exhausted = False          # no prime factor <= sqrt(m): m is 1 or prime
            break
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            f.append((p, e))
    if m > 1 and not exhausted:
        f.append((m, 1))
        m = 1
    return f, m


def gap_free(m):
    """(verdict, [(p, e), ...] as far as known, unfactored cofactor or 1)"""
    if m == 1:
        return True, [], 1
    f, C = split_small(m)
    if C == 1:
        ps = [p for p, _ in f]
        ok = all(next_prime(a) == b for a, b in zip(ps, ps[1:]))
        return ok, f, 1
    if f:
        ps = [p for p, _ in f]
        if not all(next_prime(a) == b for a, b in zip(ps, ps[1:])):
            return False, f, C
        rest = chain(C, next_prime(ps[-1]))
        if rest is None:
            return False, f, C
        return True, f + rest, 1
    # no factor below TRIAL
    if is_prime(C):
        return True, [(C, 1)], 1
    t = 2
    while TRIAL ** t < C:
        c = iroot(C, t)
        exact = c ** t == C
        if exact and is_prime(c):
            return True, [(c, t)], 1
        q = prev_prime(c) if exact else (c if is_prime(c) else prev_prime(c))
        for _ in range(t - 1):
            if q <= TRIAL:
                break
            got = chain(C, q)
            if got is not None:
                return True, got, 1
            q = prev_prime(q)
        t += 1
    return False, [], C


def brent(n, budget):
    """Pollard-Brent rho with an iteration budget; a nontrivial factor or None."""
    if n % 2 == 0:
        return 2
    used = 0
    c = 1
    while used < budget:
        y, r, q, g = 2, 1, 1, 1
        while g == 1 and used < budget:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(128, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = gcd(q, n)
                k += 128
                used += 128
            r *= 2
        if g == 1:
            return None
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = gcd(abs(x - ys), n)
        if g != n:
            return g
        c += 1
    return None


def factor_display(m, budget=300000):
    """Best-effort full factorization for display: (sorted [(p, e)], unfactored cofactor or 1)."""
    f = {}
    small, C = split_small(m)
    for p, e in small:
        f[p] = f.get(p, 0) + e
    stack = [C] if C > 1 else []
    unf = 1
    while stack:
        x = stack.pop()
        if x == 1:
            continue
        if is_prime(x):
            f[x] = f.get(x, 0) + 1
            continue
        r = isqrt(x)
        if r * r == x:
            stack += [r, r]
            continue
        d = brent(x, budget)
        if d is None:
            unf *= x
        else:
            stack += [d, x // d]
    return sorted(f.items()), unf


def show(fac, unf=1):
    out, prev = [], None
    for p, e in fac:
        if prev is not None:
            out.append(" | " if next_prime(prev) != p else "*")
        out.append(f"{p}^{e}" if e > 1 else str(p))
        prev = p
    if unf > 1:
        out.append(("*" if out else "") + f"[C{unf}]")
    return "".join(out) if out else "1"


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    L, R = int(sys.argv[1]), int(sys.argv[2])
    ok = True
    for name, v in (("L", L), ("R", R)):
        gf, fac, C = gap_free(v)
        dfac, unf = factor_display(v)
        print(f"{name} = {v} = {show(dfac, unf)}  ->  {'gap-free' if gf else '** NOT gap-free **'}")
        ok &= gf
    n = R - L - 1
    print(f"checking the {n} numbers {L+1} .. {R-1}:")
    for m in range(L + 1, R):
        gf, fac, C = gap_free(m)
        dfac, unf = factor_display(m)
        note = "" if not gf else "   ** gap-free: run is NOT maximal **"
        if unf > 1:
            note += "   (cofactor not factored; verdict from the consecutive-prime test)"
        print(f"  {m} = {show(dfac, unf)}  ->  {'gap' if not gf else 'no gap'}{note}")
        ok &= not gf
    if ok:
        print(f"OK: {L} and {R} are consecutive gap-free numbers; {L+1} starts a maximal run of exactly {n}.")
    else:
        print("FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
