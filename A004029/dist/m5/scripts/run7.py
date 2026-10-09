#!/usr/bin/env python3
"""Split 7-dimensional Q-classes into Z-classes and count space groups.

For every Q-class representative listed in REPLIST:
  1. `enantio -k` (CARAT q2z) splits it into Z-classes with normalizers and
     dumps them to OUT/<qid>/ ;
  2. `zverify` counts space-group types (affine and orientation-preserving)
     for each Z-class.
Results are appended to RESULTS/q_<tag>.tsv (one line per Q-class) and
RESULTS/z_<tag>.tsv (one line per Z-class, zverify format).  Q-classes that
already appear in q_<tag>.tsv are skipped, so the run can be restarted.

usage: run7.py REPLIST OUTDIR RESULTSDIR TAG [-j JOBS] [--split-limit S] [--count-limit S]
"""
import argparse
import multiprocessing as mp
import os
import subprocess
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def work(args):
    rep, out, split_limit, count_limit = args
    qid = os.path.basename(rep)
    d = os.path.join(out, qid)
    os.makedirs(d, exist_ok=True)
    env = dict(os.environ, ENANTIO_DUMP_DIR=d, CARAT_DIR=os.path.join(ROOT, 'carat'))
    t0 = time.time()
    try:
        r = subprocess.run([os.path.join(ROOT, 'enantio'), '-k', rep], env=env,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=split_limit)
        status = 'ok' if r.returncode == 0 else f'rc{r.returncode}'
    except subprocess.TimeoutExpired:
        status = 'split-timeout'
    t1 = time.time()
    zfiles = sorted(f for f in os.listdir(d) if '__Z' in f)
    zlines = []
    if zfiles and status == 'ok':
        try:
            r = subprocess.run([os.path.join(ROOT, 'zverify')] + zfiles, cwd=d,
                               capture_output=True, text=True, timeout=count_limit)
            zlines = [l for l in r.stdout.splitlines() if l.strip()]
        except subprocess.TimeoutExpired:
            # salvage file by file
            for z in zfiles:
                try:
                    r = subprocess.run([os.path.join(ROOT, 'zverify'), z], cwd=d,
                                       capture_output=True, text=True, timeout=count_limit)
                    zlines += [l for l in r.stdout.splitlines() if l.strip()]
                except subprocess.TimeoutExpired:
                    zlines.append(f'{z}\tTIMEOUT')
    t2 = time.time()
    return qid, status, len(zfiles), t1 - t0, t2 - t1, zlines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('replist')
    ap.add_argument('outdir')
    ap.add_argument('resultsdir')
    ap.add_argument('tag')
    ap.add_argument('-j', type=int, default=os.cpu_count())
    ap.add_argument('--split-limit', type=float, default=4 * 3600)
    ap.add_argument('--count-limit', type=float, default=4 * 3600)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(a.resultsdir, exist_ok=True)
    qfile = os.path.join(a.resultsdir, f'q_{a.tag}.tsv')
    zfile = os.path.join(a.resultsdir, f'z_{a.tag}.tsv')
    done = set()
    if os.path.exists(qfile):
        done = {l.split('\t')[0] for l in open(qfile)}
    reps = [l.split('\t')[0].strip() for l in open(a.replist) if l.strip()]
    todo = [r for r in reps if os.path.basename(r) not in done]
    print(f'{len(reps)} Q-classes listed, {len(done)} done, {len(todo)} to do, {a.j} jobs', flush=True)
    t_start = time.time()
    with mp.Pool(a.j) as pool, open(qfile, 'a') as qf, open(zfile, 'a') as zf:
        n = 0
        for qid, status, nz, ts, tc, zlines in pool.imap_unordered(
                work, [(r, a.outdir, a.split_limit, a.count_limit) for r in todo]):
            aff = prop = bad = 0
            for l in zlines:
                f = l.split('\t')
                if len(f) >= 8 and f[1] != 'ERROR':
                    aff += int(f[3])
                    prop += int(f[4])
                else:
                    bad += 1
                zf.write(l + '\n')
            bad += nz - len(zlines)
            qf.write(f'{qid}\t{status}\t{nz}\t{aff}\t{prop}\t{bad}\t{ts:.2f}\t{tc:.2f}\n')
            qf.flush()
            zf.flush()
            n += 1
            if n % 500 == 0:
                print(f'{n}/{len(todo)} done, {time.time() - t_start:.0f} s', flush=True)
    print('finished', flush=True)


if __name__ == '__main__':
    main()
