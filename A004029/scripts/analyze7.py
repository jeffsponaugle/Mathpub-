#!/usr/bin/env python3
"""Combine the 7-dimensional run results (q_*.tsv / z_*.tsv) into the lower
bounds, and list candidate duplicate Z-classes for a GL_7(Z)-conjugacy check.

usage: analyze7.py RESULTDIR [RESULTDIR...]
"""
import sys
from collections import defaultdict

import glob
qrows, zrows = [], []
seen_q, seen_z = set(), set()
for d in sys.argv[1:]:
    for f in sorted(glob.glob(f'{d}/q_*.tsv')):
        for l in open(f):
            r = l.rstrip('\n').split('\t')
            if l.strip() and r[0] not in seen_q:   # a Q-class may have been run twice
                seen_q.add(r[0])
                qrows.append(r)
    for f in sorted(glob.glob(f'{d}/z_*.tsv')):
        for l in open(f):
            r = l.rstrip('\n').split('\t')
            if l.strip() and r[0] not in seen_z:   # Z-class file names embed the Q-class id
                seen_z.add(r[0])
                zrows.append(r)

status = defaultdict(int)
for r in qrows:
    status[r[1]] += 1
qids = [r[0] for r in qrows]
dupq = len(qids) - len(set(qids))

good, bad = [], []
for r in zrows:
    (good if len(r) >= 8 and r[1] not in ('ERROR', 'TIMEOUT') else bad).append(r)
aff = sum(int(r[3]) for r in good)
prop = sum(int(r[4]) for r in good)

print(f'Q-classes processed: {len(qrows)} (duplicates in lists: {dupq}); status: {dict(status)}')
print(f'Z-classes counted:   {len(good)}   (errors/timeouts: {len(bad)})')
print(f'space groups (affine types):          {aff}')
print(f'space groups incl. enantiomorphs:     {prop}')
print(f'enantiomorphic pairs:                 {prop - aff}')
for r in bad[:10]:
    print('   bad:', r[:3])

# biggest contributors
good.sort(key=lambda r: -int(r[3]))
cum = 0
for i, r in enumerate(good):
    cum += int(r[3])
    if cum >= 0.9 * aff:
        print(f'top {i + 1} Z-classes carry 90% of the affine total')
        break
with open('results/z7_top.tsv', 'w') as f:
    for r in good[:20000]:
        f.write('\t'.join(r) + '\n')
for r in good[:8]:
    print('  ', r[0], r[1], r[3], r[7])
