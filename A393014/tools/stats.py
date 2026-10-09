#!/usr/bin/env python3
"""Observed vs. heuristically expected counts of double-base palindromes, per length.

  stats.py results/b2_9.tsv [--json]

Heuristic: a number in an interval I of base-P length L and base-Q length LQ is a base-b
palindrome with probability ~ b^(ceil(len/2) - len) (independently), so the expected count is
    E = sum_LQ |I| * P^(ceil(L/2)-L) * Q^(ceil(LQ/2)-LQ),
with a parity correction when one base is 2 (binary palindromes are odd).
Expected primes: E * 2/ln(N) for the odd numbers (every double palindrome here is odd when one
base is even), E * 1/ln(N) otherwise.  Lengths ruled out by the parity rules contribute 0.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from campaign import roles  # noqa: E402
from report import bound_for  # noqa: E402


def ndig(n, b):
    d = 0
    while n:
        n //= b
        d += 1
    return d


def expected(P, Q, L):
    """(expected double palindromes, expected primes) with base-P length L."""
    if L % 2 == 0 and (P + 1) % Q == 0:
        return 0.0, 0.0
    lo, hi = P ** (L - 1), P ** L
    e = ep = 0.0
    for LQ in range(ndig(lo, Q), ndig(hi - 1, Q) + 1):
        if LQ % 2 == 0 and (Q + 1) % P == 0:
            continue
        a, b = max(lo, Q ** (LQ - 1)), min(hi, Q ** LQ)
        if b <= a:
            continue
        x = (b - a) * float(P) ** (math.ceil(L / 2) - L) * float(Q) ** (math.ceil(LQ / 2) - LQ)
        # base 2: binary palindromes are odd, so only the odd base-P palindromes can match, and among
        # odd numbers binary palindromes are twice as dense
        if Q == 2 or P == 2:
            B = P if Q == 2 else Q
            Lb = L if Q == 2 else LQ
            if B % 2 == 1:
                frac = ((B - 1) // 2) / B if Lb % 2 == 1 else 0.0  # odd centre digit
            else:
                frac = (B // 2) / (B - 1)  # odd leading (= last) digit
            x *= 2 * frac
        e += x
        # primes need odd lengths in both bases (else divisible by P+1 or Q+1)
        if L % 2 == 1 and LQ % 2 == 1:
            mid = math.log(a) / 2 + math.log(b) / 2
            odd = 2.0 if (P % 2 == 0 or Q % 2 == 0) else 1.0
            ep += x * odd / mid
    return e, ep


def main():
    path = sys.argv[1]
    b1, b2 = map(int, os.path.basename(path)[1:-4].split("_"))
    P, Q = roles(b1, b2)
    P_, Lmax, lim, _ = bound_for(b1, b2)
    obs, obsp = {}, {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0].isdigit():
            n = int(p[0])
            if 0 < n < lim:
                L = ndig(n, P)
                obs[L] = obs.get(L, 0) + 1
                if p[1] == "1":
                    obsp[L] = obsp.get(L, 0) + 1
    rows = []
    ce = cep = co = cop = 0.0
    for L in range(1, Lmax + 1):
        e, ep = expected(P, Q, L)
        ce += e
        cep += ep
        co += obs.get(L, 0)
        cop += obsp.get(L, 0)
        rows.append({"L": L, "observed": obs.get(L, 0), "expected": e, "primes": obsp.get(L, 0),
                     "expectedPrimes": ep, "cumObserved": co, "cumExpected": ce, "cumPrimes": cop,
                     "cumExpectedPrimes": cep})
    if "--json" in sys.argv:
        print(json.dumps({"bases": [b1, b2], "P": P, "Q": Q, "rows": rows}))
        return 0
    print(f"bases {b1},{b2} (P={P}, Q={Q}); searched below {P}^{Lmax}")
    print(f"{'L':>3} {'obs':>4} {'exp':>7} {'primes':>6} {'exp.pr':>7} | {'cum obs':>7} {'cum exp':>8} {'cum pr':>6} {'cum exp.pr':>10}")
    for r in rows:
        if r["expected"] == 0 and r["observed"] == 0:
            continue
        print(f"{r['L']:>3} {r['observed']:>4} {r['expected']:>7.2f} {r['primes']:>6} {r['expectedPrimes']:>7.3f} | "
              f"{r['cumObserved']:>7.0f} {r['cumExpected']:>8.1f} {r['cumPrimes']:>6.0f} {r['cumExpectedPrimes']:>10.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
