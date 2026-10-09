#!/usr/bin/env python3
"""Literal big-integer check of the A061955 definition for small n.

V(n) = binary string rev(n) rev(n-1) ... rev(1), where rev(k) is k's binary
expansion written least-significant digit first (trailing zeros of k kept as
leading zeros).  a(n) lists the n with n | V(n).
"""
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000

s = ""
hits = []
for n in range(1, N + 1):
    s = bin(n)[2:][::-1] + s          # prepend reversed binary of n
    if int(s, 2) % n == 0:
        hits.append(n)

assert int("11101110100111011", 2) == 122171
print(hits)
