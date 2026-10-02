#!/usr/bin/env python3
"""
compare_bin0.py - compare the [1e19, 2^64) bin of the a(20) chunk logs with recount logs.

The a(20) runs (1e19..1e20, wheel 37, chunk 16384, bins starting with [1e19, 2^64)) wrote one
line per chunk; field 6 is "count:cks" of the [1e19, 2^64) bin.  A recount of `count 1e19 2^64
-w 37 -c 16384` over the same classes writes lines whose field 6 must be identical.

Usage: compare_bin0.py ORIG_SEGMENTS... -- RECOUNT_LOGS...
  ORIG_SEGMENTS as for combine_a20.py (LOG:LO:HI); RECOUNT_LOGS are plain chunk logs.
"""
import sys


def lines(path):
    d = {}
    for raw in open(path):
        f = raw.split()
        if len(f) >= 6 and f[0] == 'C':
            d[int(f[1])] = f[5]
    return d


def main():
    if '--' not in sys.argv:
        print(__doc__)
        sys.exit(2)
    i = sys.argv.index('--')
    orig = {}
    for seg in sys.argv[1:i]:
        path, lo, hi = seg.rsplit(':', 2)
        for k, v in lines(path).items():
            if int(lo) <= k < int(hi):
                orig[k] = v
    rec = {}
    for path in sys.argv[i + 1:]:
        rec.update(lines(path))
    common = sorted(set(orig) & set(rec))
    bad = [k for k in common if orig[k] != rec[k]]
    print(f"original chunks {len(orig)}, recounted chunks {len(rec)}, compared {len(common)}, mismatches {len(bad)}")
    for k in bad[:20]:
        print(f"  chunk {k}: original {orig[k]}  recount {rec[k]}")
    if common:
        so = sum(int(orig[k].split(':')[0]) for k in common)
        sr = sum(int(rec[k].split(':')[0]) for k in common)
        print(f"sum over compared chunks: original {so}, recount {sr}, difference {sr - so}")


if __name__ == '__main__':
    main()
