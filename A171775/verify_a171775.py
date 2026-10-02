#!/usr/bin/env python3
"""
Independent checks for OEIS A171775 (pure Python, no shared code with a171775.c).

a(n) = smallest M such that for every k = 2..n there is a base b_k in which M
is a k-digit palindrome.

  verify_a171775.py check M [n]     digits of M in the smallest base of each length 2..n
  verify_a171775.py small           a(1..6) by testing M = 1, 2, 3, ... directly
  verify_a171775.py brute n HI      enumerate every n-digit palindrome <= HI in every base,
                                    test the (n-1)-digit condition by trying every base, and print
                                    the same statistics as `a171775 search n 1 HI`
                                    ((n-1)-palindrome pairs, canonical survivors, checksum,
                                    failures by length, solutions)
Numbers accept 2^72 style powers.
"""
import sys


def parse(s):
    s = s.strip()
    total = 0
    for term in s.replace('-', '+-').split('+'):
        if not term:
            continue
        sign = -1 if term.startswith('-') else 1
        term = term.lstrip('-')
        if '^' in term:
            a, e = term.split('^')
            v = int(a) ** int(e)
        elif 'e' in term:
            a, e = term.split('e')
            v = int(a) * 10 ** int(e)
        else:
            v = int(term)
        total += sign * v
    return total


def iroot(x, k):
    """floor(x^(1/k)) exactly"""
    if x < 2:
        return x
    r = int(round(x ** (1.0 / k)))
    while r ** k > x:
        r -= 1
    while (r + 1) ** k <= x:
        r += 1
    return r


def digits(M, b):
    d = []
    while M:
        M, r = divmod(M, b)
        d.append(r)
    return d


def pal_bases(M, L):
    """all bases b >= 2 in which M is an L-digit palindrome"""
    if L == 1:
        return []
    if L == 2:
        # M = d(b+1), 1 <= d < b
        out = []
        d = 1
        while d * (d + 2) <= M:          # need b = M/d - 1 > d
            if M % d == 0:
                b = M // d - 1
                if b > d:
                    out.append(b)
            d += 1
            if d > 10 ** 6:              # only existence matters for huge M: b = M-1 always works
                break
        if M >= 3 and (M - 1) not in out:
            out.append(M - 1)
        return sorted(out)
    lo = max(2, iroot(M, L) + 1)
    hi = iroot(M, L - 1)
    out = []
    for b in range(lo, hi + 1):
        r = M % b
        if r == 0 or M // b ** (L - 1) != r:
            continue
        d = digits(M, b)
        if len(d) == L and d == d[::-1]:
            out.append(b)
    return out


def first_pal_base(M, L):
    if L == 2:
        return M - 1 if M >= 3 else 0
    lo = max(2, iroot(M, L) + 1)
    hi = iroot(M, L - 1)
    for b in range(lo, hi + 1):
        r = M % b
        if r == 0 or M // b ** (L - 1) != r:
            continue
        d = digits(M, b)
        if len(d) == L and d == d[::-1]:
            return b
    return 0


def surv_hash(M):
    m64 = (1 << 64) - 1
    a = ((M & m64) * 0x9E3779B97F4A7C15) & m64
    b = (((M >> 64) & m64) * 0xC2B2AE3D27D4EB4F) & m64
    a ^= a >> 29
    return (a * 0xBF58476D1CE4E5B9 + b) & m64


def cmd_check(M, n):
    ok = True
    print("M =", M)
    for L in range(2, n + 1):
        b = first_pal_base(M, L)
        if not b:
            print("  L=%2d none" % L)
            ok = False
            continue
        d = digits(M, b)[::-1]
        print("  L=%2d base %d digits %s" % (L, b, " ".join(map(str, d)) if L > 2 else "1 1"))
    print("all lengths 2..%d: %s" % (n, "yes" if ok else "NO"))
    return ok


def cmd_small():
    known = [1, 3, 5, 52, 130, 1885]
    for n in range(1, 7):
        M = 1
        while not all(first_pal_base(M, L) for L in range(2, n + 1)):
            M += 1
        print("a(%d) = %d %s" % (n, M, "ok" if M == known[n - 1] else "WRONG"))


def palindromes(B, n, hi):
    """all n-digit palindromes in base B that are <= hi, ascending"""
    h = n // 2
    D = (n + 1) // 2
    w = [B ** (n - 1 - i) + B ** i for i in range(h)]
    if n % 2:
        w.append(B ** h)
    rest = [(B - 1) * sum(w[i + 1:]) for i in range(D)]   # max contribution of the digits after i

    def rec(i, acc):
        if i == D:
            yield acc
            return
        start = 1 if i == 0 else 0
        for d in range(start, B):
            v = acc + d * w[i]
            if v > hi:
                return
            if v + rest[i] < 1:
                continue
            yield from rec(i + 1, v)
    yield from rec(0, 0)


def cmd_brute(n, hi):
    st2 = 0
    canon = {}
    B = 2
    npal = 0
    while B ** (n - 1) <= hi:
        for M in palindromes(B, n, hi):
            npal += 1
            cs = pal_bases(M, n - 1)
            st2 += len(cs)
            if cs and M not in canon:
                canon[M] = (B, cs[0])
        B += 1
    fail = {}
    sols = []
    for M in sorted(canon):
        for L in range(n - 2, 2, -1):
            if not first_pal_base(M, L):
                fail[L] = fail.get(L, 0) + 1
                break
        else:
            sols.append(M)
    csum = 0
    for M in canon:
        csum = (csum + surv_hash(M)) & ((1 << 64) - 1)
    for M in sols:
        print("SOLUTION n=%d M=%d bases(L=3..%d)=%s" % (n, M, n, ",".join(str(first_pal_base(M, L)) for L in range(3, n + 1))))
    print("# brute n=%d [1, %d]: %d n-digit palindromes, (n-1)-palindromes %d, canonical %d, checksum %016x, solutions %d"
          % (n, hi, npal, st2, len(canon), csum, len(sols)))
    print("# failed at length:" + "".join(" L=%d:%d" % (L, fail.get(L, 0)) for L in range(n - 2, 2, -1)))


def main():
    a = sys.argv[1:]
    if not a:
        print(__doc__)
        sys.exit(1)
    if a[0] == 'check':
        sys.exit(0 if cmd_check(parse(a[1]), int(a[2]) if len(a) > 2 else 12) else 2)
    if a[0] == 'small':
        cmd_small()
        return
    if a[0] == 'brute':
        cmd_brute(int(a[1]), parse(a[2]))
        return
    print(__doc__)
    sys.exit(1)


if __name__ == '__main__':
    main()
