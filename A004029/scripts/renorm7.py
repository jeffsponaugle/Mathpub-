#!/usr/bin/env python3
"""Normalizer cross-check for 7-dimensional Z-classes.

For each listed Z-class file (path relative to OUTDIR, i.e. <qid>/<zfile>),
recompute the normalizer from scratch with CARAT's standalone Normalizer,
then count again with zverify using the union of both generating sets
("A+B").  Both sets consist of verified normalizer elements, so the union
can only be closer to the full normalizer: a smaller count would reveal an
incomplete normalizer in the original run.

usage: renorm7.py LIST OUTDIR RESULTFILE [-j JOBS] [--limit S]
output: zfile \t original-affine \t original-proper \t union-affine \t union-proper \t status
"""
import argparse
import multiprocessing as mp
import os
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def find(name):
    for c in (os.path.join(ROOT, name), os.path.join(ROOT, 'carat', 'bin', name)):
        if os.path.exists(c):
            return c
    raise SystemExit(f'{name} not found')


NORMALIZER = find('Normalizer')
ZVERIFY = find('zverify_v2') if os.path.exists(os.path.join(ROOT, 'zverify_v2')) else find('zverify')


def work(args):
    rel, outdir, limit = args
    d, z = os.path.split(os.path.join(outdir, rel))
    try:
        r = subprocess.run([NORMALIZER, z], cwd=d, capture_output=True, text=True, timeout=limit)
        if r.returncode != 0 or not r.stdout.strip():
            return rel, None, None, 'normalizer-failed'
        with open(os.path.join(d, z + '.N2'), 'w') as f:
            f.write(r.stdout)
        r = subprocess.run([ZVERIFY, z, z + '+' + z + '.N2'], cwd=d, capture_output=True, text=True, timeout=limit)
        lines = [l.split('\t') for l in r.stdout.splitlines()]
        if len(lines) != 2 or 'ERROR' in (lines[0][1], lines[1][1]):
            return rel, None, None, 'zverify-failed: ' + r.stdout.strip()[:200]
        return rel, lines[0][3:5], lines[1][3:5], 'ok'
    except subprocess.TimeoutExpired:
        return rel, None, None, 'timeout'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('list')
    ap.add_argument('outdir')
    ap.add_argument('result')
    ap.add_argument('-j', type=int, default=os.cpu_count())
    ap.add_argument('--limit', type=float, default=1800)
    a = ap.parse_args()
    items = [l.split('\t')[0].strip() for l in open(a.list) if l.strip()]
    with mp.Pool(a.j) as pool, open(a.result, 'w') as out:
        for rel, o, u, st in pool.imap_unordered(work, [(i, a.outdir, a.limit) for i in items]):
            o = o or ['', '']
            u = u or ['', '']
            out.write(f'{rel}\t{o[0]}\t{o[1]}\t{u[0]}\t{u[1]}\t{st}\n')
            out.flush()


if __name__ == '__main__':
    main()
