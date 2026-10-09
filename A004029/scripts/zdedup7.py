#!/usr/bin/env python3
"""Certify that the Z-classes of each Q-class are pairwise non-conjugate.

Input: a Z-fingerprint table (extend1 -z output: path \\t |G| \\t zfp ...),
paths of the form OUTDIR/<qid>/<zfile>.  Z-classes of one Q-class with
different fingerprints are certainly distinct (the fingerprint is a
GL_n(Z)-conjugacy invariant).  Within each group of equal fingerprints every
member is tested with CARAT's Z_equiv against the representatives found so
far.  A member counts as new only if Z_equiv answers "not conjugated" for all
of them; otherwise (conjugate, or any unclear answer) it is listed as a
duplicate and must be dropped from the counts.

usage: zdedup7.py ZFPTABLE DUPFILE [-j JOBS] [--limit S]
output DUPFILE: zfile \\t reason   (one line per Z-class to drop)
"""
import argparse
import multiprocessing as mp
import os
import subprocess
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZEQ = next(p for p in (os.path.join(ROOT, 'Z_equiv'), os.path.join(ROOT, 'carat', 'bin', 'Z_equiv'))
           if os.path.exists(p))


def zequiv(a, b, limit):
    try:
        r = subprocess.run([ZEQ, a, b], capture_output=True, text=True, timeout=limit)
    except subprocess.TimeoutExpired:
        return 'timeout'
    out = r.stdout
    if 'not conjugated' in out:
        return 'distinct'
    if 'conjugates the group' in out:
        return 'conjugate'
    return 'unclear:' + (out + r.stderr).strip().replace('\n', ' ')[:120]


def work(args):
    members, limit = args
    reps, drops, calls = [members[0]], [], 0
    for m in members[1:]:
        verdict = 'distinct'
        for r in reps:
            calls += 1
            v = zequiv(r, m, limit)
            if v != 'distinct':
                verdict = f'{v} with {os.path.basename(r)}'
                break
        if verdict == 'distinct':
            reps.append(m)
        else:
            drops.append((m, verdict))
    return drops, calls


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('table')
    ap.add_argument('dupfile')
    ap.add_argument('-j', type=int, default=os.cpu_count())
    ap.add_argument('--limit', type=float, default=600)
    a = ap.parse_args()
    groups = defaultdict(list)
    for l in open(a.table):
        f = l.rstrip('\n').split('\t')
        if len(f) < 3:
            continue
        q = os.path.basename(os.path.dirname(f[0]))
        groups[(q, f[2])].append(f[0])
    todo = [(sorted(v), a.limit) for v in groups.values() if len(v) > 1]
    print(f'{sum(len(v) for v in groups.values())} Z-classes, {len(groups)} fingerprint groups, '
          f'{len(todo)} groups to check', flush=True)
    ndrop = ncalls = 0
    with mp.Pool(a.j) as pool, open(a.dupfile, 'w') as out:
        for drops, calls in pool.imap_unordered(work, todo, chunksize=16):
            ncalls += calls
            for m, why in drops:
                out.write(f'{os.path.basename(m)}\t{why}\n')
                ndrop += 1
    print(f'Z_equiv calls: {ncalls}; Z-classes to drop: {ndrop}', flush=True)


if __name__ == '__main__':
    main()
