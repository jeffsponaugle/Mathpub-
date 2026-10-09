#!/usr/bin/env python3
"""Independent exhaustive check of an nrepeat result (A331881/A331882).

usage: verify_exhaustive.py PIFILE N CLAIMED_END [CLAIMED_SUBSTRING]

Counts every N-digit window in the first CLAIMED_END fractional digits of the
file with numpy and verifies that exactly one substring reaches N occurrences
there (and, if given, that it equals CLAIMED_SUBSTRING with its N-th
occurrence ending exactly at CLAIMED_END).

Memory: ~16 bytes per digit of prefix (numpy int64 x2), so n=9's ~481M-digit
prefix needs ~8 GB; fine for n<=9, too big for n>=10 (use the tool's own
verification pass there).
"""
import sys
import numpy as np

path, N, END = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
claimed = sys.argv[4] if len(sys.argv) > 4 else None

with open(path) as f:
    s = f.read(END + 2)
if s.startswith("3."):
    s = s[2:2 + END]
else:
    s = s[:END]
assert s.startswith("14159265358979323846"), "not fractional pi digits"
assert len(s) == END, f"file has only {len(s)} digits, need {END}"

d = np.frombuffer(s.encode(), dtype=np.uint8).astype(np.int64) - ord("0")
assert d.min() >= 0 and d.max() <= 9, "non-digit bytes in prefix"
W = END - N + 1
vals = np.zeros(W, dtype=np.int64)
for k in range(N):
    vals += d[k:k + W] * 10 ** (N - 1 - k)
uniq, cnt = np.unique(vals, return_counts=True)
hot = uniq[cnt >= N]
print(f"windows checked: {W}")
print(f"max count in prefix: {cnt.max()}")
print(f"substrings with >= {N} occurrences:", [f"{v:0{N}d}" for v in hot])

ok = len(hot) == 1
if claimed:
    ok = ok and f"{hot[0]:0{N}d}" == claimed
    starts = np.nonzero(vals == hot[0])[0] + 1  # 1-based
    end = int(starts[N - 1]) + N - 1
    print(f"occurrence starts: {list(map(int, starts[:N]))}")
    print(f"{N}-th occurrence end: {end}")
    ok = ok and end == END
print("VERIFIED" if ok else "MISMATCH")
sys.exit(0 if ok else 1)
