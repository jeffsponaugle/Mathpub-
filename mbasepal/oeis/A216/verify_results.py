#!/usr/bin/env python3
"""Independent check of every gappal result: v has exactly L digits and is a
palindrome in both bases, the bases differ by d, and known OEIS terms match."""
import glob, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))

def digits(v, b):
    d = []
    while v:
        v, r = divmod(v, b)
        d.append(r)
    return d

def known(L):
    p = os.path.join(HERE, 'known', f'L{L}.txt')
    out = {}
    if os.path.exists(p):
        for line in open(p):
            s = line.split()
            if len(s) >= 2 and not line.startswith('#'):
                out[int(s[0])] = int(s[1])
    return out

dirs = sys.argv[1:] or ['results']
for dname in dirs:
    for path in sorted(glob.glob(os.path.join(HERE, dname, 'L*.tsv')), key=lambda p: int(os.path.basename(p)[1:-4])):
        L = int(os.path.basename(path)[1:-4])
        kn = known(L)
        rows = {}
        bad = []
        for line in open(path):
            s = line.split()
            d, v, b1, b2 = map(int, s[:4])
            rows[d] = v
            ok = b2 - b1 == d
            for b in (b1, b2):
                ds = digits(v, b)
                ok = ok and len(ds) == L and ds == ds[::-1]
            if d in kn and kn[d] != v:
                ok = False
            if not ok:
                bad.append(d)
        off = 1 if L % 2 else 2
        ds = sorted(rows)
        contiguous = ds == list(range(off, off + len(ds)))
        print(f'{dname}/L{L}: {len(rows)} terms d={ds[0]}..{ds[-1]}, contiguous={contiguous}, '
              f'{sum(1 for d in rows if d in kn)} vs OEIS, failures={bad}')
