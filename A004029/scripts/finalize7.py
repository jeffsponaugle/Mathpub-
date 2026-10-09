#!/usr/bin/env python3
"""Final numbers of the 7-dimensional computation.

Reads, from each given results directory:
  q_*.tsv       per Q-class status (run7.py)
  z_*.tsv       per Z-class counts (zverify via run7.py)
  zdups_*.tsv   Z-classes to drop (zdedup7.py)
  renorm_*.tsv  normalizer cross-check (renorm7.py)
and prints the lower bounds for A004028(7), A004027(7), A004029(7),
A006227(7) and A395859(7).

usage: finalize7.py REPLIST RESULTDIR [RESULTDIR...]
"""
import glob
import sys

replist = [l.split('\t')[0].strip() for l in open(sys.argv[1]) if l.strip()]
dirs = sys.argv[2:]

q, z, drops, renorm = {}, {}, {}, {}
for d in dirs:
    for f in sorted(glob.glob(f'{d}/q_*.tsv')):
        for l in open(f):
            r = l.rstrip('\n').split('\t')
            q.setdefault(r[0], r)
    for f in sorted(glob.glob(f'{d}/z_*.tsv')):
        for l in open(f):
            r = l.rstrip('\n').split('\t')
            z.setdefault(r[0], r)
    for f in sorted(glob.glob(f'{d}/zdups_*.tsv')):
        for l in open(f):
            r = l.rstrip('\n').split('\t')
            if r[0]:
                drops[r[0]] = r[1] if len(r) > 1 else ''
    for f in sorted(glob.glob(f'{d}/renorm_*.tsv')):
        for l in open(f):
            r = l.rstrip('\n').split('\t')
            renorm[r[0].split('/')[-1]] = r

reps = {r.split('/')[-1] for r in replist}
ok = {k for k, r in q.items() if r[1] == 'ok'}
excluded = sorted(reps - ok)

bad = [k for k, r in z.items() if len(r) < 8 or r[1] in ('ERROR', 'TIMEOUT')]
counted = {k: r for k, r in z.items() if k not in drops and k not in bad}
# every counted Z-class must belong to a completed Q-class
orphans = [k for k in counted if k[3:].rsplit('__Z', 1)[0] not in ok]
for k in orphans:
    del counted[k]

aff = sum(int(r[3]) for r in counted.values())
prop = sum(int(r[4]) for r in counted.values())

checked = [k for k, r in renorm.items() if k in counted and r[5] == 'ok']
diffs = [k for k in checked if renorm[k][1:3] != renorm[k][3:5]]
aff_checked = sum(int(counted[k][3]) for k in checked)

print(f'7-dim Q-classes with a 1-dim rational constituent (distinct fingerprints): {len(reps)}')
print(f'  split into Z-classes:        {len(ok)}')
print(f'  excluded (CARAT failure or stopped): {len(excluded)}')
for e in excluded:
    print(f'      {e}  ({q[e][1] if e in q else "not run"})')
print(f'Z-classes counted:             {len(counted)}   '
      f'(dropped: {len(drops)} unclear Z_equiv, {len(bad)} errors, {len(orphans)} orphans)')
print(f'space-group types (affine):    {aff}')
print(f'incl. enantiomorphs (proper):  {prop}')
print(f'enantiomorphic pairs:          {prop - aff}')
print(f'normalizer cross-check: {len(checked)} Z-classes recomputed, {len(diffs)} differences, '
      f'covering {100.0 * aff_checked / aff:.2f}% of the affine total')
print()
print('Lower bounds:')
print(f'  A004028(7) >= {len(reps)}')
print(f'  A004027(7) >= {len(counted)}')
print(f'  A004029(7) >= {aff}')
print(f'  A006227(7) >= {prop}')
print(f'  A395859(7) >= {prop - aff}')
