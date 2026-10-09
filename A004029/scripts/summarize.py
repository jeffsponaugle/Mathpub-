#!/usr/bin/env python3
"""Summarize zverify results per dimension and compare with the OEIS terms.

Reads results/zverify_dim{1..6}.tsv (one line per Z-class:
file |G| H^1 affine proper pairs improper method) and prints

  n  Q-classes  Z-classes  chiral-Z  affine  proper  pairs

A Z-class is chiral (its GL_n(Z)-class splits into two SL_n(Z)-classes) iff
proper == 2*affine: the split extension is fixed by the whole normalizer, so it
splits into an enantiomorphic pair exactly when N <= SL_n(Z).
"""
import os
import sys

OEIS = {
    'A004027 (Z-classes)':                 [1, 2, 13, 73, 710, 6079, 85308],
    'A307288 (Z-classes incl. enant.)':    [1, 2, 13, 73, 780, 6079, None],
    'A004029 (space groups)':              [1, 2, 17, 219, 4783, 222018, 28927915],
    'A006227 (space groups incl. enant.)': [1, 2, 17, 230, 4894, 222097, None],
    'A395859 (enantiomorphic pairs)':      [0, 0, 0, 11, 111, 79, None],
}

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
rows = {0: dict(Q=1, Z=1, chiral=0, affine=1, proper=1, pairs=0, maxsg=1, maxH='1')}
for n in range(1, 7):
    path = os.path.join(root, 'results', f'zverify_dim{n}.tsv')
    if not os.path.exists(path):
        continue
    qs = set()
    r = dict(Q=0, Z=0, chiral=0, affine=0, proper=0, pairs=0, maxsg=0, maxH='')
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        qs.add(f[0].rsplit('__Z', 1)[0])
        a, p = int(f[3]), int(f[4])
        r['Z'] += 1
        r['affine'] += a
        r['proper'] += p
        r['chiral'] += (p == 2 * a)
        if a > r['maxsg']:
            r['maxsg'], r['maxH'] = a, f[2]
    r['Q'] = len(qs)
    r['pairs'] = r['proper'] - r['affine']
    rows[n] = r

print(f"{'n':>2} {'Q':>6} {'Z':>7} {'chiralZ':>8} {'Z+chiral':>9} {'affine':>11} {'proper':>11} {'pairs':>6}  largest single Z-class")
for n, r in sorted(rows.items()):
    h = r['maxH']
    hs = h if len(h) < 14 else f"2^{h.count('2')}" if set(h.split('x')) == {'2'} else h[:12] + '..'
    print(f"{n:>2} {r['Q']:>6} {r['Z']:>7} {r['chiral']:>8} {r['Z'] + r['chiral']:>9} {r['affine']:>11} "
          f"{r['proper']:>11} {r['pairs']:>6}  {r['maxsg']} groups (H^1 = {hs})")

print()
ok = True
for name, vals in OEIS.items():
    key = {'A004027': 'Z', 'A307288': None, 'A004029': 'affine', 'A006227': 'proper', 'A395859': 'pairs'}[name[:7]]
    got = []
    for n in range(7):
        if n not in rows:
            got.append(None)
            continue
        r = rows[n]
        got.append(r['Z'] + r['chiral'] if key is None else r[key])
    status = []
    for n, (want, have) in enumerate(zip(vals, got)):
        if have is None:
            status.append('?')
        elif want is None:
            status.append(f'NEW a({n})={have}')
        elif want == have:
            status.append('ok')
        else:
            status.append(f'MISMATCH a({n}): oeis {want} vs {have}')
            ok = False
    print(f"{name:38s} " + ', '.join(status))
sys.exit(0 if ok else 1)
