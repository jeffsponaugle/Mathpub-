#!/usr/bin/env python3
"""
verify_wieferich.py -- independent checks for wieferich.c, using Python integers only.

  verify_wieferich.py range LO HI [BASES]   prime count and per-base checksums of [LO, HI),
                                            printed like a line of the tool's chunk log (-L)
  verify_wieferich.py check P [BASES]        Fermat quotients q_P(b) and A = q or q - P
  verify_wieferich.py known LIMIT [BASES]    every solution p < LIMIT by brute force
  verify_wieferich.py cksheader              regenerate wieferich_cks.h (selftest checksums)

Primes come from a small-prime sieve of the window plus deterministic Miller-Rabin
(the first 12 prime bases, exact below 3.18e23), not from primesieve.  The checksum
is the tool's definition:

    cks_b = sum over primes p >= 3 with p not dividing b of
            mix64(q ^ (p * 0x9E3779B97F4A7C15 mod 2^64))  mod 2^64,
    q     = (b^(p-1) mod p^2 - 1) / p      (the Fermat quotient),

with mix64 the splitmix64 finalizer.  BASES: comma list / ranges, 'oeis' or 'known'.
"""
import sys

M64 = (1 << 64) - 1
GOLD = 0x9E3779B97F4A7C15
OEIS = [3, 5, 6, 7, 10, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 26, 30]
KNOWN = [2, 3, 5, 6, 7, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 26, 28, 29, 30]
MR_BASES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]


def mix64(z):
    z = (z + GOLD) & M64
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
    return z ^ (z >> 31)


def parse_num(s):
    s = s.strip()
    for i in range(len(s) - 1, 0, -1):
        if s[i] in '+-' and s[i - 1] not in 'eE':
            a, b = parse_num(s[:i]), parse_num(s[i + 1:])
            return a + b if s[i] == '+' else a - b
    if '^' in s:
        b, e = s.split('^')
        return int(b) ** int(e)
    if 'e' in s.lower():
        m, e = s.lower().split('e')
        if '.' in m:
            ip, fp = m.split('.')
            return int(ip + fp) * 10 ** (int(e) - len(fp))
        return int(m) * 10 ** int(e)
    return int(s)


def parse_bases(s):
    out = set()
    for tok in s.split(','):
        if tok == 'oeis':
            out.update(OEIS)
        elif tok == 'known':
            out.update(KNOWN)
        elif '-' in tok:
            lo, hi = tok.split('-')
            out.update(range(int(lo), int(hi) + 1))
        else:
            out.add(int(tok))
    return sorted(out)


def small_primes(n):
    """primes < n by a plain sieve"""
    s = bytearray([1]) * n
    s[0:2] = b'\0\0'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(range(i * i, n, i)))
    return [i for i in range(n) if s[i]]


def is_prime(n):
    if n < 2:
        return False
    for q in MR_BASES:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


SMALL = small_primes(1 << 16)


def primes_in(lo, hi):
    """primes in [lo, hi): sieve by primes < 2^16, then Miller-Rabin"""
    n = hi - lo
    alive = bytearray([1]) * n
    for q in SMALL:
        if q * q >= hi:
            break
        first = max(q * q, (lo + q - 1) // q * q)
        if first < hi:
            alive[first - lo::q] = bytearray(len(range(first - lo, n, q)))
    for i in range(n):
        if alive[i]:
            p = lo + i
            if p >= 2 and (p < (1 << 32) or is_prime(p)):
                yield p


def quotient(b, p):
    r = pow(b, p - 1, p * p)
    assert (r - 1) % p == 0
    return (r - 1) // p


def window(lo, hi, bases):
    cks = {b: 0 for b in bases}
    count = 0
    for p in primes_in(lo, hi):
        count += 1
        if p == 2:
            continue
        ph = (p * GOLD) & M64
        for b in bases:
            if b % p == 0:
                continue
            cks[b] = (cks[b] + mix64(quotient(b, p) ^ ph)) & M64
    return count, cks


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    cmd = argv[1]
    if cmd == 'range':
        lo, hi = parse_num(argv[2]), parse_num(argv[3])
        bases = parse_bases(argv[4] if len(argv) > 4 else 'oeis')
        count, cks = window(lo, hi, bases)
        print('C %d %d %d %s' % (lo, hi, count, ' '.join('%d:%016x' % (b, cks[b]) for b in bases)))
    elif cmd == 'check':
        p = parse_num(argv[2])
        bases = parse_bases(argv[3] if len(argv) > 3 else 'known')
        print('p = %d (%s)' % (p, 'prime' if is_prime(p) else 'NOT prime'))
        for b in bases:
            if b % p == 0:
                print('base %3d: p divides b' % b)
                continue
            q = quotient(b, p)
            A = q if q <= p // 2 else q - p
            print('base %3d: q = %d  A = %d%s' % (b, q, A, '  <-- WIEFERICH' if q == 0 else ''))
    elif cmd == 'known':
        limit = parse_num(argv[2])
        bases = parse_bases(argv[3] if len(argv) > 3 else 'known')
        ps = small_primes(limit)
        for b in bases:
            sol = [p for p in ps if b % p and pow(b, p - 1, p * p) == 1]
            print('base %d: %s' % (b, ' '.join(map(str, sol)) or 'none'))
    elif cmd == 'cksheader':
        wins = [(10 ** 15, 10 ** 15 + 4000000), (10 ** 18, 10 ** 18 + 2000000), (2 ** 63 - 2000000, 2 ** 63)]
        print('/* generated by verify_wieferich.py cksheader: {lo, hi, primes, cks for bases oeis} */')
        for lo, hi in wins:
            count, cks = window(lo, hi, OEIS)
            print('{%dULL, %dULL, %dULL, {%s}},' % (lo, hi, count, ', '.join('0x%016xULL' % cks[b] for b in OEIS)))
            sys.stdout.flush()
    else:
        print(__doc__)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
