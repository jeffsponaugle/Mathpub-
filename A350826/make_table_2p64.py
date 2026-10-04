#!/usr/bin/env python3
"""
make_table_2p64.py - table of pi_6(x) for x in [10^19, 2^64] from the two halves of
`a350826_cuda count 1e19 2^64 -w 37 -c 16384 -b $(cat runs/bounds_2p64.txt)` (runs/t64a.*, runs/t64b.*).

Adds the per-bin PART counts of both halves, checks the total against the a(20) run
(19907602564 sextuplets in [1e19, 2^64)) and, chunk by chunk, against the [1e19, 2^64) recount
(runs/v64a.chunks, runs/v64b.chunks), and writes runs/pi6_table_1e19_2p64.txt and .md.
"""
import re

PI6_1E19 = 28722086297
TOTAL = 19907602564                      # [1e19, 2^64) from the a(20) run and the recount
T = 2 ** 64


def part_lines(path):
    bins = {}
    for line in open(path):
        m = re.match(r'PART \[(\d+), (\d+)\) count (\d+) cks ([0-9a-f]+)$', line.strip())
        if m:
            a, b, c, h = int(m[1]), int(m[2]), int(m[3]), m[4]
            if (a, b) != (10 ** 19, T) or len(bins) == 0:
                bins[(a, b)] = (c, int(h, 16))
    return bins


def label(x):
    if x == T:
        return '2^64'
    if x % 10 ** 17 == 0:
        return f'{x // 10 ** 17 / 100:.2f}e19'
    if x % 2 ** 58 == 0:
        return f'{x // 2 ** 58}*2^58'
    d = T - x
    if d & (d - 1) == 0:
        return f'2^64-2^{d.bit_length() - 1}'
    return str(x)


def chunk_totals(paths, nfields=None):
    tot = {}
    for path in paths:
        for raw in open(path):
            f = raw.split()
            if len(f) >= 6 and f[0] == 'C':
                tot[int(f[1])] = sum(int(x.split(':')[0]) for x in f[5:])
    return tot


a = part_lines('runs/t64a.out')
b = part_lines('runs/t64b.out')
keys = sorted(k for k in a if k != (10 ** 19, T))
assert keys == sorted(k for k in b if k != (10 ** 19, T)), "the two halves have different bins"
rows, cum, total = [], PI6_1E19, 0
for k in keys:
    c = a[k][0] + b[k][0]
    total += c
    cum += c
    rows.append((k[0], k[1], c, cum))
print(f"{len(rows)} intervals, total {total} (a(20) run: {TOTAL}) {'OK' if total == TOTAL else 'MISMATCH'}")
print(f"pi_6(2^64) = {cum} (table at pzktupel.de: 48629687343, difference {cum - 48629687343})")

# chunk by chunk against the recount (sum over the bins of each chunk line)
new = chunk_totals(['runs/t64a.chunks', 'runs/t64b.chunks'])
old = chunk_totals(['runs/v64a.chunks', 'runs/v64b.chunks'])
common = set(new) & set(old)
bad = [k for k in common if new[k] != old[k]]
print(f"chunk totals vs the [1e19, 2^64) recount: {len(common)} chunks compared, {len(bad)} mismatches")

with open('runs/pi6_table_1e19_2p64.txt', 'w') as f:
    f.write("# pi_6(x): number of prime sextuplets (p, p+4, p+6, p+10, p+12, p+16) with p < x\n")
    f.write("# (no sextuplet straddles any x below), from pi_6(10^19) = 28722086297 (Desfontaines; confirmed)\n")
    f.write("# x  label  sextuplets_in_[previous_x,x)  pi_6(x)\n")
    f.write(f"{10 ** 19} 1.00e19 - {PI6_1E19}\n")
    for lo, hi, c, cum in rows:
        f.write(f"{hi} {label(hi)} {c} {cum}\n")
with open('runs/pi6_table_1e19_2p64.md', 'w') as f:
    f.write("| x | label | sextuplets in [prev x, x) | pi_6(x) |\n|---|---|---|---|\n")
    f.write(f"| {10 ** 19} | 1.00e19 | | {PI6_1E19} |\n")
    for lo, hi, c, cum in rows:
        f.write(f"| {hi} | {label(hi)} | {c} | {cum} |\n")
print("wrote runs/pi6_table_1e19_2p64.txt and .md")
