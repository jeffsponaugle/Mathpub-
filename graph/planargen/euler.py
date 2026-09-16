#!/usr/bin/env python3
"""Derived sequences from the connected-planar-by-edges table (A049334).

Given rows T(n,k) = number of connected planar graphs with n nodes and k
edges (n = 1..N), compute by the bivariate Euler transform (disjoint unions):
  * A039735 rows: all planar graphs with n nodes and k edges,
  * A005470: all planar graphs with n nodes (row sums),
and print A003094 (row sums of the input) for reference.

usage: euler.py rows.json [N]      (rows.json as written by verify/aggregate:
                                    {"1": [..], "2": [..], ...})
"""
import json, sys
from collections import defaultdict

def euler_bivariate(T, N):
    """T: dict n -> dict k -> count of connected structures.
    Returns U: dict n -> dict k -> count of all structures (multisets of
    connected ones), using the standard recurrence with weights on n only:
      n*U(n,k) = sum_{j=1..n} sum_{d|j} d * T(d, .) ... done as a formal power
    series in two variables via the multiset (Euler) product."""
    # Represent series as dict (n,k) -> coefficient; product over all
    # connected types c of 1/(1 - x^n y^k)^{T(n,k)}.  Use the logarithmic
    # derivative form: n*U(n,k) = sum over (n1,k1) with n1<=n of
    #   [sum_{d | (n1,k1)-multiples} ...] -- simpler: iterate over each
    #   connected type (n1,k1,count) and multiply by (1-x^n1 y^k1)^(-count)
    #   using the binomial series, truncated to n<=N.
    from math import comb
    U = defaultdict(int); U[(0, 0)] = 1
    for n1 in range(1, N + 1):
        for k1, cnt in T.get(n1, {}).items():
            if cnt == 0: continue
            # multiply U by sum_{m>=0} C(cnt+m-1, m) x^{m n1} y^{m k1}
            newU = defaultdict(int)
            for (n, k), v in U.items():
                m = 0
                while n + m * n1 <= N:
                    newU[(n + m * n1, k + m * k1)] += v * comb(cnt + m - 1, m)
                    m += 1
            U = newU
    out = defaultdict(dict)
    for (n, k), v in U.items():
        if n >= 1: out[n][k] = v
    return out

def main():
    rows = json.load(open(sys.argv[1]))
    N = int(sys.argv[2]) if len(sys.argv) > 2 else max(int(x) for x in rows)
    missing = [n for n in range(1, N + 1) if str(n) not in rows]
    if missing:
        sys.exit(f"euler.py: rows {missing} missing from {sys.argv[1]}; the Euler transform "
                 f"needs every row 1..{N} (add them to a049334_rows.json, e.g. from results_n14/)")
    T = {int(n): {k: v for k, v in enumerate(r) if v} for n, r in rows.items() if int(n) <= N}
    U = euler_bivariate(T, N)
    print("A003094 (connected planar):", [sum(T[n].values()) for n in range(1, N + 1)])
    print("A005470 (all planar)      :", [sum(U[n].values()) for n in range(1, N + 1)])
    for n in range(1, N + 1):
        kmax = max(U[n]) if U[n] else 0
        row = [U[n].get(k, 0) for k in range(0, kmax + 1)]
        print(f"A039735 row {n} (k=0..{kmax}): {row}")

if __name__ == "__main__":
    main()
