#!/usr/bin/env python3
"""Independent checks of the proposed b-files (run from this directory)."""
def rb(f):
    return [tuple(map(int, l.split())) for l in open(f) if l.strip() and not l.startswith('#')]
pal = lambda x: str(x) == str(x)[::-1]
ok = True
def check(cond, msg):
    global ok
    print(("ok   " if cond else "FAIL ") + msg); ok &= cond
r = rb('b002778.txt'); s = rb('b002779.txt'); e = rb('b016113.txt'); q = rb('b027829.txt'); c = rb('b263618.txt')
check([n for n, _ in r] == list(range(1, len(r) + 1)), "A002778 indices 1..%d" % len(r))
check(all(pal(v * v) for _, v in r), "A002778 every root has a palindromic square")
check(all(r[i][1] < r[i + 1][1] for i in range(len(r) - 1)), "A002778 strictly increasing")
check([v for _, v in s] == [v * v for _, v in r], "A002779 = A002778^2 termwise")
check([v for _, v in e] == [v for _, v in r if v > 0 and len(str(v * v)) % 2 == 0], "A016113 = even-length subset of A002778")
check([v for _, v in q] == [v * v for _, v in e], "A027829 = A016113^2 termwise")
cnt = {}
for _, v in r: cnt[len(str(v * v))] = cnt.get(len(str(v * v)), 0) + 1
check(all(v == cnt.get(n, 0) for n, v in c), "A263618 = length histogram of A002779 (n = 1..%d)" % len(c))
check(c[-2:] == [(68, 1), (69, 1443)], "A263618 new terms a(68)=1, a(69)=1443")
check(e[-1] == (23, 7256171055736382499839982033184475), "A016113 a(23) = 7256171055736382499839982033184475")
for name, f in [('A002778', 'orig_b002778.txt'), ('A002779', 'orig_b002779.txt'), ('A016113', 'orig_b016113.txt'),
                ('A027829', 'orig_b027829.txt'), ('A263618', 'orig_b263618.txt')]:
    o = rb(f); n = {'A002778': r, 'A002779': s, 'A016113': e, 'A027829': q, 'A263618': c}[name]
    check(n[:len(o)] == o, f"{name} extends the current OEIS b-file ({len(o)} -> {len(n)} terms)")
zeros = [n for n, v in c if v == 0 and n % 2 == 0]
check(zeros == [2, 4, 8, 10, 14, 18, 20, 24, 30, 38, 40, 46, 54, 56, 62, 64],
      "zeros of A263618 (n <= 69) are exactly A034822: " + ",".join(map(str, zeros)))
print("ALL OK" if ok else "SOME CHECKS FAILED")
