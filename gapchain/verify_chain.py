#!/usr/bin/env python3
"""Independent check of a chain of L primes with decreasing gaps (A263049,
the default) or increasing gaps (A016045, with "inc"): every number from p to
its last member is tested with 13-base Miller-Rabin (deterministic below
3.3e24) and with Baillie-PSW; the primes must be exactly the L members.
Usage: verify_chain.py p L [inc|dec]"""
import math, sys
def mr(n, bases=(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)):
    if n < 2: return False
    for q in bases:
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in bases:
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True
def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0
def lucas_strong(n):
    D = 5
    while True:
        j = jacobi(D, n)
        if j == -1: break
        if j == 0 and abs(D) != n: return False
        D = -D - 2 if D > 0 else -D + 2
    P, Q = 1, (1 - D) // 4
    d, s = n + 1, 0
    while d % 2 == 0: d //= 2; s += 1
    U, V, Qk = 0, 2, 1          # U_0, V_0, Q^0
    for bit in bin(d)[2:]:      # left to right
        U, V = U * V % n, (V * V - 2 * Qk) % n
        Qk = Qk * Qk % n
        if bit == "1":
            U, V = (P * U + V), (D * U + P * V)
            U = (U + n) * pow(2, -1, n) % n if U % 2 else U // 2 % n
            V = (V + n) * pow(2, -1, n) % n if V % 2 else V // 2 % n
            Qk = Qk * Q % n
    if U == 0 or V == 0: return True
    for _ in range(s - 1):
        V = (V * V - 2 * Qk) % n
        Qk = Qk * Qk % n
        if V == 0: return True
    return False
def bpsw(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    if not mr(n, (2,)): return False
    if math.isqrt(n) ** 2 == n: return False
    return lucas_strong(n)
p = int(sys.argv[1]); L = int(sys.argv[2])
inc = len(sys.argv) > 3 and sys.argv[3] == "inc"
span = (L - 1) * L                      # gaps 2, 4, ..., 2(L-1) in either order
primes = []
for n in range(p, p + span + 1):
    a, b = mr(n), bpsw(n)
    if a != b: print("TESTS DISAGREE at", n); sys.exit(1)
    if a: primes.append(n)
gaps = [b - a for a, b in zip(primes, primes[1:])]
want = list(range(2, 2 * L, 2)) if inc else list(range(2 * (L - 1), 0, -2))
print(f"{len(primes)} primes from {p} to {p + span}; gaps {gaps}")
print("VERIFIED" if primes[0] == p and gaps == want else "NOT A CHAIN")
# is the chain longer?  An increasing chain can only grow at its end (a gap
# of 2L after its last member), a decreasing one only at its start (2L before p)
if inc:
    q = p + span + 1
    while not mr(q): q += 1
    print(f"next prime {q} (gap {q - p - span}); a gap of {2 * L} there would make a longer chain: {'yes' if q - p - span == 2 * L else 'no'}")
else:
    q = p - 1
    while not mr(q): q -= 1
    print(f"previous prime {q} (gap {p - q}); a gap of {2 * L} there would make a longer chain: {'yes' if p - q == 2 * L else 'no'}")
