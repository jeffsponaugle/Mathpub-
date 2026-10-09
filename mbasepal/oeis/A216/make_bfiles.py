#!/usr/bin/env python3
"""Write OEIS b-files from gappal results (results/ and results_studio/).

Only the contiguous prefix of terms starting at the sequence offset is
written, so a b-file never has gaps. Every term is re-checked (exactly L
digits and a palindrome in both bases, bases differ by d) before writing.
"""
import glob, os

HERE = os.path.dirname(os.path.abspath(__file__))
ANUM = {2: 216841, 3: 216840, 4: 216843, 5: 216899, 6: 216900, 7: 216901, 8: 216902,
        9: 216903, 10: 216904, 11: 216905, 12: 216906, 13: 216907, 14: 216908,
        15: 216909, 17: 216910}


def digits(v, b):
    d = []
    while v:
        v, r = divmod(v, b)
        d.append(r)
    return d


def main():
    rows = {}
    for dname in ('results', 'results_studio'):
        for path in glob.glob(os.path.join(HERE, dname, 'L*.tsv')):
            L = int(os.path.basename(path)[1:-4])
            for line in open(path):
                d, v, b1, b2 = map(int, line.split()[:4])
                prev = rows.setdefault(L, {}).get(d)
                if prev and prev[0] != v:
                    raise SystemExit(f'conflict L={L} d={d}: {prev[0]} vs {v}')
                rows[L][d] = (v, b1, b2)
    os.makedirs(os.path.join(HERE, 'bfiles'), exist_ok=True)
    for L in sorted(rows):
        off = 1 if L % 2 else 2
        n = off
        while n in rows[L]:
            n += 1
        last = n - 1
        if last < off:
            continue
        for d in range(off, last + 1):
            v, b1, b2 = rows[L][d]
            for b in (b1, b2):
                ds = digits(v, b)
                assert len(ds) == L and ds == ds[::-1] and b2 - b1 == d, (L, d)
        a = ANUM.get(L)
        name = f'b{a}.txt' if a else f'b_new_L{L}.txt'
        title = (f'A{a}: Smallest palindromic number of length {L} in two bases differing by n.'
                 if a else f'(new) Smallest palindromic number of length {L} in two bases differing by n.')
        with open(os.path.join(HERE, 'bfiles', name), 'w') as f:
            f.write(f'# {title}\n')
            f.write(f'# n = {off}..{last}: exhaustive search over base pairs (b, b+n) in increasing\n')
            f.write('# larger base, meet-in-the-middle with exact 192-bit arithmetic (gappal.c).\n')
            f.write('# Computed 2026 by Jeff Sponaugle.\n')
            for d in range(off, last + 1):
                f.write(f'{d} {rows[L][d][0]}\n')
        print(f'L={L:2d}: {name}  n={off}..{last} ({last - off + 1} terms)')


if __name__ == '__main__':
    main()
