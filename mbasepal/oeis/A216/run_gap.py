#!/usr/bin/env python3
"""Run gappal over a range of base gaps d for one palindrome length L.

Each (L, d) is an independent exact search ("smallest L-digit palindrome in
bases b and b+d"). Jobs run in parallel as separate processes; finished
terms are appended to results/L<L>.tsv as "d value b1 b2 seconds", so an
interrupted run resumes where it left off. Every result is compared with
known/L<L>.txt (OEIS data or b-file) when that term is known.

    ./run_gap.py -L 8 --dmax 200 -j 14
"""
import argparse, os, re, subprocess, sys, threading, time
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
RES = re.compile(r'RESULT L=(\d+) d=(\d+) v=(\d+) b1=(\d+) b2=(\d+)')


def load_tsv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            p = line.split()
            if len(p) >= 4 and not line.startswith('#'):
                out[int(p[0])] = (int(p[1]), int(p[2]), int(p[3]))
    return out


def load_known(L):
    out = {}
    path = os.path.join(HERE, 'known', f'L{L}.txt')
    if os.path.exists(path):
        for line in open(path):
            p = line.split()
            if len(p) >= 2 and not line.startswith('#'):
                out[int(p[0])] = int(p[1])
    return out


def run_one(L, d, threads, extra):
    cmd = [os.path.join(HERE, 'gappal'), '-n', str(L), '-D', str(d),
           '-t', str(threads), '-q'] + (['-M'] if L >= 8 else []) + extra
    t0 = time.time()
    with open(os.path.join(HERE, 'logs', f'L{L}_d{d}.log'), 'w') as lf:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=lf, text=True)
    secs = time.time() - t0
    m = RES.search(p.stdout)
    return d, (int(m.group(3)), int(m.group(4)), int(m.group(5))) if m else None, secs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-L', type=int, required=True)
    ap.add_argument('--dmin', type=int, default=0, help='default: offset (1 odd L, 2 even L)')
    ap.add_argument('--dmax', type=int, required=True)
    ap.add_argument('-j', '--jobs', type=int, default=os.cpu_count())
    ap.add_argument('-t', '--threads', type=int, default=1, help='threads per gappal process')
    ap.add_argument('--extra', default='', help='extra gappal flags, e.g. "-T 700000000"')
    a = ap.parse_args()

    L = a.L
    dmin = a.dmin or (1 if L % 2 else 2)
    os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    tsv = os.path.join(HERE, 'results', f'L{L}.tsv')
    done = load_tsv(tsv)
    known = load_known(L)
    todo = [d for d in range(dmin, a.dmax + 1) if d not in done]
    print(f'L={L}: d={dmin}..{a.dmax}, {len(done)} already done, {len(todo)} to run, '
          f'{a.jobs} jobs x {a.threads} threads', flush=True)

    lock = threading.Lock()
    bad = []
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run_one, L, d, a.threads, a.extra.split()) for d in todo]
        for n, f in enumerate(as_completed(futs), 1):
            d, r, secs = f.result()
            if r is None:
                print(f'  d={d}: NO RESULT (see logs/L{L}_d{d}.log)', flush=True)
                bad.append(d)
                continue
            v, b1, b2 = r
            with lock:
                with open(tsv, 'a') as out:
                    out.write(f'{d} {v} {b1} {b2} {secs:.2f}\n')
            flag = ''
            if d in known:
                flag = '  [matches OEIS]' if known[d] == v else f'  [MISMATCH: OEIS has {known[d]}]'
                if known[d] != v:
                    bad.append(d)
            print(f'  [{n}/{len(todo)} {time.time()-t0:7.1f}s] d={d}: {v} (bases {b1},{b2}) '
                  f'{secs:.1f}s{flag}', flush=True)

    done = load_tsv(tsv)
    have = [d for d in range(dmin, a.dmax + 1) if d in done]
    checked = [d for d in have if d in known]
    mism = [d for d in checked if known[d] != done[d][0]]
    print(f'L={L} summary: {len(have)}/{a.dmax - dmin + 1} terms for d={dmin}..{a.dmax}; '
          f'{len(checked)} checked against OEIS, {len(mism)} mismatches {mism}; '
          f'failures {bad}; wall {time.time()-t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
