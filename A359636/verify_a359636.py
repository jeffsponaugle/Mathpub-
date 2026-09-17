#!/usr/bin/env python3
"""
Independent cross-check for the A359636 tool (a359636.c).

  verify_a359636.py verify N P            check that the prime gap after P qualifies at level N
                                          (own Miller-Rabin, own Pollard-rho factoring)
  verify_a359636.py brute N LIMIT         brute-force a(N) by sieving omega for every integer
                                          up to LIMIT (numpy) and walking every prime gap directly;
                                          no use of the m = 3 (mod 6) reduction of the C tool
  verify_a359636.py check-log FILE ...    re-verify every "SOLUTION n=.. p=.." line of a tool log
  verify_a359636.py selftest              brute-force a(1..5) and verify the OEIS terms a(1..8)

Requires numpy only for 'brute' and 'selftest'.
"""
import sys, math, random, re

KNOWN = {1: 7, 2: 19, 3: 643, 4: 51427, 5: 8083633, 6: 1077940147,
         7: 75582271489, 8: 34710483181813}

SMALL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]

def is_prime(n):
    if n < 2:
        return False
    for p in SMALL:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in SMALL:                      # deterministic below 3.3e24
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

def pollard_rho(n):
    if n % 2 == 0:
        return 2
    while True:
        c = random.randrange(1, n)
        x = y = random.randrange(0, n)
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = math.gcd(abs(x - y), n)
        if d != n:
            return d

def factor(n):
    """sorted list of (prime, exponent)"""
    fs = {}
    for p in SMALL + list(range(41, 1000, 2)):
        while n % p == 0:
            fs[p] = fs.get(p, 0) + 1
            n //= p
    stack = [n] if n > 1 else []
    while stack:
        x = stack.pop()
        if x == 1:
            continue
        if is_prime(x):
            fs[x] = fs.get(x, 0) + 1
            continue
        d = pollard_rho(x)
        stack += [d, x // d]
    return sorted(fs.items())

def fstr(fs):
    return "*".join(f"{p}^{e}" if e > 1 else str(p) for p, e in fs)

def next_prime(x):
    x += 1
    while not is_prime(x):
        x += 1
    return x

def verify(n, p, quiet=False):
    if not is_prime(p):
        print(f"{p} is not prime")
        return False
    q = next_prime(p)
    ok = q - p >= 4
    if not quiet:
        print(f"p = {p}, next prime {q}, gap {q - p}" + ("" if ok else "  (twin prime, excluded)"))
    for x in range(p + 1, q):
        fs = factor(x)
        w = len(fs)
        if not quiet:
            print(f"  {x} = {fstr(fs)}  omega={w}" + ("" if w >= n else "  <-- fails"))
        if w < n:
            ok = False
    if not quiet:
        print(f"gap after {p} {'QUALIFIES' if ok else 'does not qualify'} at level {n}")
    return ok

def brute(n, limit):
    """a(n) by direct definition for p <= limit, or None."""
    import numpy as np
    N = limit + 2000                              # room for the prime after limit
    isp = np.ones(N + 1, dtype=bool)
    isp[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if isp[i]:
            isp[i * i::i] = False
    omega = np.zeros(N + 1, dtype=np.uint8)
    primes = np.nonzero(isp)[0]
    cut = N // 1000
    for p in primes[primes <= cut]:               # small primes: strided slices
        omega[p::p] += 1
    big = primes[primes > cut]                    # large primes: at most 1000 multiples each
    k = 1
    while k * big[0] <= N:
        idx = big[k * big <= N] * k
        omega[idx] += 1
        k += 1
    bad = (~isp) & (omega < n)                    # composites that fail the level
    bad[:2] = False
    cb = np.concatenate(([0], np.cumsum(bad, dtype=np.int64)))   # cb[i] = #bad in [0, i)
    ps = primes[primes <= limit]
    qs = primes[1:len(ps) + 1]
    gaps = qs - ps
    good = (gaps >= 4) & (cb[qs] - cb[ps + 1] == 0)
    hits = ps[good]
    return int(hits[0]) if len(hits) else None

def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    cmd = argv[1]
    if cmd == "verify":
        n, p = int(argv[2]), int(argv[3])
        return 0 if verify(n, p) else 1
    if cmd == "brute":
        n, limit = int(argv[2]), int(float(argv[3]))
        r = brute(n, limit)
        print(f"brute force: a({n}) = {r}" if r else f"brute force: a({n}) > {limit}")
        return 0
    if cmd == "check-log":
        bad = 0
        for fn in argv[2:]:
            for line in open(fn):
                m = re.match(r"SOLUTION n=(\d+) p=(\d+) q=(\d+)", line)
                if not m:
                    continue
                n, p, q = map(int, m.groups())
                ok = verify(n, p, quiet=True) and next_prime(p) == q
                print(f"{fn}: n={n} p={p} q={q} -> {'ok' if ok else 'FAIL'}")
                bad += not ok
        return 1 if bad else 0
    if cmd == "selftest":
        fails = 0
        for n in range(1, 6):
            r = brute(n, KNOWN[n])
            ok = r == KNOWN[n]
            print(f"brute a({n}) = {r} expected {KNOWN[n]} -> {'ok' if ok else 'FAIL'}")
            fails += not ok
        for n in range(1, 9):
            ok = verify(n, KNOWN[n], quiet=True)
            print(f"verify gap after a({n}) = {KNOWN[n]} at level {n} -> {'ok' if ok else 'FAIL'}")
            fails += not ok
            if n >= 2:                             # the gap must not already qualify one level up
                ok2 = not verify(n + 1, KNOWN[n], quiet=True)
                print(f"   ... and not at level {n + 1} -> {'ok' if ok2 else 'FAIL'}")
                fails += not ok2
        print("selftest", "FAILED" if fails else "passed")
        return 1 if fails else 0
    print(__doc__)
    return 2

if __name__ == "__main__":
    sys.exit(main(sys.argv))
