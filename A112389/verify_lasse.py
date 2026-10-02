#!/usr/bin/env python3
"""Recompute refinements from Deleuran's sum-for-size.py with ./lego -r and compare.

usage: verify_lasse.py SUM_FOR_SIZE_PY SIZE [max_seconds_per_run] [profiles...]

Runs each listed profile (default: all bottleneck-free profiles of SIZE whose
aggregated plan enumerates at most 4 bricks) and prints OK / MISMATCH.
"""
import re, subprocess, sys

src = open(sys.argv[1]).read()
size = int(sys.argv[2])
tmax = int(sys.argv[3]) if len(sys.argv) > 3 else 1800
want = sys.argv[4:]

i = src.index('# Size %d:' % size)
j = src.index('# Size %d:' % (size + 1)) if ('# Size %d:' % (size + 1)) in src else src.index('# Impossible')
known = {}
for k, tot, sym in re.findall(r"'(\d+)',\s*(\d+),\s*\(?(\d+)\)?", src[i:j]):
    known[k] = (int(tot), int(sym))


def plan(prof):
    h = len(prof)
    best = None
    for m in range(1, 1 << h):
        if m & (m >> 1) or m == (1 << h) - 1:
            continue
        if any(prof[l] > 4 for l in range(h) if m >> l & 1):
            continue
        rb = sum(prof[l] for l in range(h) if not m >> l & 1)
        key = (rb, max(prof[l] for l in range(h) if m >> l & 1))
        if best is None or key < best[0]:
            best = (key, m)
    return best


todo = want or [k for k in known if plan([int(c) for c in k]) and plan([int(c) for c in k])[0][0] <= 4]
for k in todo:
    prof = [int(c) for c in k]
    (rb, kk), m = plan(prof)
    agg = ''.join(str(l) for l in range(len(prof)) if m >> l & 1)
    try:
        out = subprocess.run(['perl', '-e', 'alarm shift; exec @ARGV', str(tmax), './lego', '-r', k, '-a', agg],
                             capture_output=True, text=True).stdout.strip()
    except Exception as e:  # noqa
        out = ''
    mm = re.match(r'<(\d+)> (\d+) \((\d+)\)', out)
    if not mm:
        print('%-8s agg=%-4s TIMEOUT/FAIL' % (k, agg), flush=True)
        continue
    got = (int(mm.group(2)), int(mm.group(3)))
    exp = known.get(k)
    tag = 'OK' if exp == got else 'MISMATCH (Deleuran %s)' % (exp,)
    t = re.search(r'\[F ([\d.]+)s', out)
    print('%-8s agg=%-4s %s (%s)  %s  [%ss]' % (k, agg, got[0], got[1], tag, t.group(1) if t else '?'), flush=True)
