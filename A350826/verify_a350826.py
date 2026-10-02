#!/usr/bin/env python3
"""
verify_a350826.py - independent check of a350826.c / cuda/a350826_cuda.cu.

Counts prime sextuplets (p, p+4, p+6, p+10, p+12, p+16) with LO <= p < HI and
the checksum used by the C and CUDA tools,

    cks = sum of mix64(lo64(p) ^ hi64(p) * 0x9E3779B97F4A7C15) mod 2^64,

with code that shares nothing with them: candidates p = 97 + 210 k are sieved
with numpy (one byte per k, every prime 11 <= q <= 2^16 striking the six
residues -d (mod q)), and each survivor is tested with Python's own big-integer
pow(): a base-2 strong test of each member, then Miller-Rabin with the first 13
prime bases, which is deterministic below 3.3e24 (Sorenson & Webster).  Initial
members below 2^16 (where a member could equal a sieving prime) are tested
directly, which covers p = 7.

Usage:
    verify_a350826.py count LO HI       count and checksum (also prints the p's if few)
    verify_a350826.py oeis N            A350826(1..N) by the same method
    verify_a350826.py check P           is P the start of a prime sextuplet?
Numbers may be written as 10^k, 2^k, 1e18, or sums/differences of those.
"""
import sys
import numpy as np

OFF = (0, 4, 6, 10, 12, 16)
MASK = (1 << 64) - 1
GOLD = 0x9E3779B97F4A7C15
BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
BOUND = 1 << 16
SEG = 1 << 24                      # candidates per numpy segment


def mix64(z):
    z = (z + GOLD) & MASK
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
    return z ^ (z >> 31)


def hashp(p):
    return mix64((p & MASK) ^ (((p >> 64) * GOLD) & MASK))


def strong(n, a):
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(s - 1):
        x = x * x % n
        if x == n - 1:
            return True
    return False


def is_prime(n):
    if n < 2:
        return False
    for b in BASES:
        if n == b:
            return True
        if n % b == 0:
            return False
    if n < 43 * 43:
        return True
    assert n < 3317044064679887385961981
    return all(strong(n, b) for b in BASES)


def is_sext(p):
    return all(is_prime(p + d) for d in OFF)


def small_primes(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return [int(x) for x in np.nonzero(s)[0]]


SP = [q for q in small_primes(BOUND) if q > 7]


def count(lo, hi, show=False):
    """count and checksum of sextuplets with lo <= p < hi"""
    n = cks = 0
    found = []
    # direct part: p < BOUND (also catches p = 7)
    for p in range(max(lo, 0), min(hi, BOUND)):
        if p % 30 == 7 and is_sext(p):
            n += 1
            cks = (cks + hashp(p)) & MASK
            found.append(p)
    a = max(lo, BOUND)
    if a < hi:
        k0 = (a - 97 + 209) // 210            # first k with 97 + 210 k >= a
        k1 = (hi - 97 + 209) // 210           # first k with 97 + 210 k >= hi
        inv = {q: pow(210, -1, q) for q in SP}
        for s0 in range(k0, k1, SEG):
            m = min(SEG, k1 - s0)
            bad = np.zeros(m, dtype=bool)
            base = 97 + 210 * s0
            for q in SP:
                r = base % q
                for d in OFF:
                    # 97 + 210 (s0 + j) + d = 0 (mod q)  <=>  j = -(r + d) / 210 (mod q)
                    j = (-(r + d) * inv[q]) % q
                    bad[j::q] = True
            for j in np.nonzero(~bad)[0]:
                p = base + 210 * int(j)
                if all(strong(p + d, 2) for d in OFF) and is_sext(p):
                    n += 1
                    cks = (cks + hashp(p)) & MASK
                    if show and len(found) < 100:
                        found.append(p)
    return n, cks, found


def parse(s):
    s = s.strip()
    for i in range(len(s) - 1, 0, -1):
        if s[i] in '+-' and s[i - 1] not in 'eE':
            a, b = parse(s[:i]), parse(s[i + 1:])
            return a + b if s[i] == '+' else a - b
    if '^' in s:
        b, e = s.split('^')
        return int(b) ** int(e)
    if 'e' in s or 'E' in s:
        m, e = s.lower().split('e')
        if '.' in m:
            ip, fp = m.split('.')
            return int(ip + fp) * 10 ** (int(e) - len(fp))
        return int(m) * 10 ** int(e)
    return int(s)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    cmd = sys.argv[1]
    if cmd == 'count' and len(sys.argv) == 4:
        lo, hi = parse(sys.argv[2]), parse(sys.argv[3])
        n, cks, found = count(lo, hi, show=True)
        if len(found) <= 20:
            for p in found:
                print(p)
        print(f"[{lo}, {hi}) count {n} cks {cks:016x}")
    elif cmd == 'oeis':
        N = int(sys.argv[2])
        for k in range(1, N + 1):
            n, cks, _ = count(10 ** (k - 1), 10 ** k)
            print(f"A350826({k}) = {n} cks {cks:016x}", flush=True)
    elif cmd == 'check':
        p = parse(sys.argv[2])
        ok = is_sext(p)
        print(f"{p} {'IS' if ok else 'is NOT'} the start of a prime sextuplet")
        sys.exit(0 if ok else 1)
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == '__main__':
    main()
