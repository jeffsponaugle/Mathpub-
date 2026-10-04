#!/usr/bin/env python3
"""Independent rows of A112390 (buildings by height) from ./lego.

usage: rows_a112390.py SYM_FILE NMAX [DELEURAN_SUM_FOR_SIZE_PY]

SYM_FILE is the output of `./lego -SV N` (S180/S90 for every profile).  Every
bottleneck-free profile of size <= NMAX is counted with ./lego (layer
aggregation when it enumerates <= 5 bricks, plain enumeration otherwise);
profiles with an interior layer of size 1 are composed:
    F(<R1 1 R2>) = F(<R1 1>) F(<1 R2>) / 2,  S180 likewise, S90 = 0.
Prints each row T(n,2..n) and, if given, compares profile values with
Deleuran's table.
"""
import itertools, re, subprocess, sys

symfile, nmax = sys.argv[1], int(sys.argv[2])
lasse = {}
if len(sys.argv) > 3:
    src = open(sys.argv[3]).read()
    for k, tot, sym in re.findall(r"'(\d+)',\s*(\d+),\s*\(?(\d+)\)?", src):
        lasse[k] = (int(tot), int(sym))

S180, S90 = {}, {}
for line in open(symfile):
    m = re.match(r'\s+<(\d+)> S180=(\d+) S90=(\d+)', line)
    if m:
        S180[m.group(1)] = int(m.group(2))
        S90[m.group(1)] = int(m.group(3))


def comps(n):
    for cuts in itertools.product([0, 1], repeat=n - 1):
        prof, run = [], 1
        for c in cuts:
            if c:
                prof.append(run); run = 1
            else:
                run += 1
        prof.append(run)
        yield ''.join(map(str, prof))


def plan(p):
    z = [int(c) for c in p]
    h = len(z)
    best = None
    for m in range(1, (1 << h) - 1):
        if m & (m >> 1) or any(z[l] > 4 for l in range(h) if m >> l & 1):
            continue
        r = sum(z[l] for l in range(h) if not m >> l & 1)
        k = max(z[l] for l in range(h) if m >> l & 1)
        if best is None or (r, k) < best[0]:
            best = ((r, k), m)
    return best


F = {}


def fat(p):
    return all(c != '1' for c in p[1:-1])


def run_lego(p):
    z = [int(c) for c in p]
    if len(z) == 2:
        args = ['./lego', '-b', p]
    else:
        pl = plan(p)
        if pl and pl[0][0] <= 5:
            agg = ''.join(str(l) for l in range(len(z)) if pl[1] >> l & 1)
            args = ['./lego', '-r', p, '-a', agg]
        else:
            args = ['./lego', '-b', p]
    out = subprocess.run(args, capture_output=True, text=True).stdout
    m = re.search(r'F=(\d+)', out)
    if not m:
        raise SystemExit('failed: %s\n%s' % (' '.join(args), out))
    return int(m.group(1)), ' '.join(args[1:])


for n in range(1, nmax + 1):
    row = {}
    for p in comps(n):
        h = len(p)
        if h == 1:
            F[p] = 2 if n == 1 else 0
            how = 'single layer'
        elif fat(p):
            F[p], how = run_lego(p)
        else:
            b = next(i for i in range(1, h - 1) if p[i] == '1')
            F[p] = F[p[:b + 1]] * F[p[b:]] // 2
            how = 'composed'
        s180, s90 = S180.get(p, 0), S90.get(p, 0)
        num = F[p] + s180 + 2 * s90
        assert num % 4 == 0 and (s180 + s90) % 2 == 0, p
        cnt, sym = num // 4, (s180 + s90) // 2
        row[h] = row.get(h, 0) + cnt
        if fat(p) and h > 1:
            tag = ''
            if p in lasse or p[::-1] in lasse:
                tag = 'OK' if lasse.get(p, lasse.get(p[::-1])) == (cnt, sym) else 'MISMATCH %s' % (lasse.get(p, lasse.get(p[::-1])),)
            print('  n=%d <%s> %d (%d)  [%s] %s' % (n, p, cnt, sym, how, tag), flush=True)
    if n >= 2:
        print('ROW %d: %s   total %d' % (n, ', '.join(str(row.get(k, 0)) for k in range(2, n + 1)),
                                       sum(row.values())), flush=True)
