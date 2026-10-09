#!/usr/bin/env python3
"""Independent exact check of one Z-class: H^1(G, Q^n/Z^n), the action of the
normalizer on it, and the number of space-group types (affine and
orientation-preserving).  Shares no code with CARAT; only reads a CARAT
bravais_TYP file (generators + normalizer generators).

  t: cocycle values on the generators, t in Q^(n*k).
  t(g) for every g in G is T_g t, built along a BFS tree of the Cayley graph;
  every non-tree edge gives the constraint (T_g + g E_i - T_{g g_i}) t in Z^n.
  H^1 = {t : C t integral} / (ker C + Z^(nk)) = (+) Z/S_jj  (SNF of C).
  n in N acts by (n.t)(g_i) = n t(n^-1 g_i n).

usage: verify_zclass.py file [file...]
"""
import sys
from fractions import Fraction


def parse_bravais(path):
    toks = open(path).read().split('\n')
    lines = [l.split('%')[0].rstrip() for l in toks]
    pos = 0

    def nextline():
        nonlocal pos
        while pos < len(lines) and lines[pos].strip() == '':
            pos += 1
        l = lines[pos]
        pos += 1
        return l

    counts = dict(g=0, f=0, z=0, n=0, c=0)
    first = nextline().strip()
    if first.startswith('#'):
        import re
        for key, val in re.findall(r'([gfznc])(\d+)', first):
            counts[key] = int(val)
    else:
        counts['g'] = 1
        pos -= 1
    # remaining numbers stream, read matrices one by one
    def read_matrix():
        head = nextline().strip()
        kgv = 1
        if '/' in head:
            head, k = head.split('/')
            kgv = int(k)
        rows = int(''.join(ch for ch in head.split('x')[0].split('d')[0] if ch.isdigit()))
        sym = diag = scalar = False
        cols = rows
        if 'x' in head:
            c = int(head.split('x')[1])
            if c == 0:
                sym = True
            else:
                cols = c
        elif 'd' in head:
            c = int(head.split('d')[1])
            diag = True
            scalar = (c == 0)
        M = [[Fraction(0)] * cols for _ in range(rows)]
        def nums(cnt):
            out = []
            nonlocal pos
            while len(out) < cnt:
                l = nextline()
                out += l.split()
            if len(out) != cnt:
                raise ValueError('matrix row length mismatch in %s' % path)
            return [Fraction(x) for x in out]
        if scalar:
            v = nums(1)[0]
            for i in range(rows):
                M[i][i] = v
        elif diag:
            v = nums(rows)
            for i in range(rows):
                M[i][i] = v[i]
        else:
            for i in range(rows):
                v = nums(i + 1 if sym else cols)
                for j, x in enumerate(v):
                    M[i][j] = x
                    if sym:
                        M[j][i] = x
        return tuple(tuple(x / kgv for x in r) for r in M)

    mats = {}
    for key in 'gfznc':
        mats[key] = [read_matrix() for _ in range(counts[key])]
    return mats


def mul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)) for i in range(n))


def inv_unimodular(A):
    n = len(A)
    M = [list(map(Fraction, r)) + [Fraction(int(i == j)) for j in range(n)] for i, r in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(tuple(r[n:]) for r in M)


def det(A):
    n = len(A)
    M = [list(map(Fraction, r)) for r in A]
    d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return 0
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return d


def row_lattice_basis(rows, ncols):
    """Integer row reduction (Hermite) -> basis of the row lattice."""
    rows = [list(r) for r in rows if any(r)]
    basis = []
    col = 0
    while rows and col < ncols:
        nz = [r for r in rows if r[col] != 0]
        zr = [r for r in rows if r[col] == 0]
        if not nz:
            col += 1
            continue
        # gcd-combine all nz rows into one pivot row
        piv = nz[0]
        rest = []
        for r in nz[1:]:
            a, b = piv[col], r[col]
            while b != 0:
                q = a // b
                piv, r = r, [x - q * y for x, y in zip(piv, r)]
                a, b = b, a - q * b
            rest.append(r)  # r[col] == 0 now
        basis.append(piv)
        rows = zr + [r for r in rest if any(r)]
        col += 1
    return basis


def snf_with_V(W, ncols):
    """Return (diag, V) with U W V = diag (V unimodular, ncols x ncols)."""
    A = [list(r) for r in W]
    R = len(A)
    V = [[int(i == j) for j in range(ncols)] for i in range(ncols)]

    def colop_swap(i, j):
        for r in A:
            r[i], r[j] = r[j], r[i]
        for r in V:
            r[i], r[j] = r[j], r[i]

    def colop_add(dst, src, f):  # col dst += f * col src
        for r in A:
            r[dst] += f * r[src]
        for r in V:
            r[dst] += f * r[src]

    diag = []
    t = 0
    while t < min(R, ncols):
        # find nonzero entry with smallest abs value in submatrix
        best = None
        for i in range(t, R):
            for j in range(t, ncols):
                if A[i][j] != 0 and (best is None or abs(A[i][j]) < abs(A[best[0]][best[1]])):
                    best = (i, j)
        if best is None:
            break
        i, j = best
        A[t], A[i] = A[i], A[t]
        colop_swap(t, j)
        while True:
            done = True
            p = A[t][t]
            for j in range(t + 1, ncols):
                if A[t][j] != 0:
                    q = A[t][j] // p
                    colop_add(j, t, -q)
                    if A[t][j] != 0:
                        done = False
            for i in range(t + 1, R):
                if A[i][t] != 0:
                    q = A[i][t] // p
                    A[i] = [x - q * y for x, y in zip(A[i], A[t])]
                    if A[i][t] != 0:
                        done = False
            if done:
                # divisibility condition
                bad = None
                for i in range(t + 1, R):
                    for j in range(t + 1, ncols):
                        if A[i][j] % p != 0:
                            bad = i
                            break
                    if bad is not None:
                        break
                if bad is None:
                    break
                A[t] = [x + y for x, y in zip(A[t], A[bad])]
                continue
            # move the smallest entry of row t / col t to the pivot
            best = None
            for j in range(t, ncols):
                if A[t][j] != 0 and (best is None or abs(A[t][j]) < abs(A[t][best[1]] if best[0] == 'c' else A[best[1]][t])):
                    best = ('c', j)
            for i in range(t, R):
                if A[i][t] != 0:
                    cur = abs(A[t][best[1]]) if best[0] == 'c' else abs(A[best[1]][t])
                    if abs(A[i][t]) < cur:
                        best = ('r', i)
            if best[0] == 'c':
                colop_swap(t, best[1])
            else:
                A[t], A[best[1]] = A[best[1]], A[t]
        if A[t][t] < 0:
            A[t] = [-x for x in A[t]]
        diag.append(A[t][t])
        t += 1
    return diag, V


def analyze(path, verbose=True):
    M = parse_bravais(path)
    gens = [tuple(tuple(int(x) for x in r) for r in g) for g in M['g']]
    norms = [tuple(tuple(int(x) for x in r) for r in g) for g in M['n'] + M['c']]
    n = len(gens[0])
    k = len(gens)
    I = tuple(tuple(int(i == j) for j in range(n)) for i in range(n))
    # enumerate G with BFS tree; T[g] is n x nk integer matrix
    nk = n * k
    zero = tuple(tuple(0 for _ in range(nk)) for _ in range(n))
    T = {I: zero}
    order = [I]
    constraints = []
    qi = 0
    while qi < len(order):
        g = order[qi]
        qi += 1
        for i, gi in enumerate(gens):
            h = mul(g, gi)
            # value: T_g + g E_i
            val = [list(r) for r in T[g]]
            for a in range(n):
                for b in range(n):
                    val[a][i * n + b] += g[a][b]
            if h not in T:
                T[h] = tuple(tuple(r) for r in val)
                order.append(h)
            else:
                for a in range(n):
                    constraints.append([val[a][c] - T[h][a][c] for c in range(nk)])
    Gset = set(order)
    W = row_lattice_basis(constraints, nk)
    diag, V = snf_with_V(W, nk)
    r = len(diag)
    H = [(j, d) for j, d in enumerate(diag) if d > 1]
    mods = [d for _, d in H]
    Vinv = inv_unimodular(tuple(tuple(x for x in row) for row in V))
    # check normalizer elements
    for nm in norms:
        ninv = inv_unimodular(nm)
        for g in gens:
            h = mul(mul(ninv, g), nm)
            if tuple(tuple(int(x) for x in row) for row in h) not in Gset:
                raise RuntimeError('normalizer element does not normalize G')

    def action(nm):
        ninv = inv_unimodular(nm)
        # A_n: nk x nk, block row i = n * T_{n^-1 g_i n}
        A = []
        for gi in gens:
            h = tuple(tuple(int(x) for x in row) for row in mul(mul(ninv, gi), nm))
            blk = mul(nm, T[h])
            A += [list(row) for row in blk]
        # on H: columns
        res = []
        for (j, d) in H:
            # s = e_j / d ; t = V s ; t' = A t ; s' = Vinv t'
            t = [Fraction(V[a][j], d) for a in range(nk)]
            tp = [sum(A[a][b] * t[b] for b in range(nk)) for a in range(nk)]
            sp = [sum(Vinv[a][b] * tp[b] for b in range(nk)) for a in range(nk)]
            col = []
            for (jj, dd) in H:
                v = sp[jj] * dd
                if v.denominator != 1:
                    raise RuntimeError('non-integral coordinate')
                col.append(int(v) % dd)
            # coordinates jj < r with diag 1 must be integral
            for jj in range(r):
                if diag[jj] == 1 and sp[jj].denominator != 1:
                    raise RuntimeError('trivial coordinate not integral')
            res.append(col)
        # res[c][rr] = image of basis c in coordinate rr -> matrix rows rr
        m = len(H)
        return [[res[c][rr] for c in range(m)] for rr in range(m)]

    acts = [(action(nm), 1 if det(nm) > 0 else -1) for nm in norms]
    improper = any(det(g) < 0 for g in gens)
    if improper:
        m = len(H)
        acts.append(([[int(a == b) for b in range(m)] for a in range(m)], -1))
    # explicit orbit count on H x {+-1}
    m = len(H)
    size = 1
    for d in mods:
        size *= d

    def decode(v):
        x = []
        for d in mods:
            x.append(v % d)
            v //= d
        return x

    def encode(x):
        v, rad = 0, 1
        for xi, d in zip(x, mods):
            v += (xi % d) * rad
            rad *= d
        return v

    seen = set()
    orbH = orbHS = 0
    for v0 in range(size):
        for sgn in (1, -1):
            if (v0, sgn) in seen:
                continue
            if sgn == 1:
                orbH += 1
            orbHS += 1
            stack = [(v0, sgn)]
            seen.add((v0, sgn))
            while stack:
                v, s = stack.pop()
                x = decode(v)
                for A, ds in acts:
                    y = [sum(A[a][b] * x[b] for b in range(m)) for a in range(m)]
                    st = (encode(y), s * ds)
                    if st not in seen:
                        seen.add(st)
                        stack.append(st)
    if verbose:
        print(f'{path}: |G|={len(order)} H^1={"x".join("Z/%d" % d for d in mods) or "0"} '
              f'normalizer gens={len(norms)} affine={orbH} proper={orbHS} pairs={orbHS - orbH}')
        for idx, (A, ds) in enumerate(acts):
            print(f'   gen {idx} det {ds:+d}: {A}')
    return orbH, orbHS


if __name__ == '__main__':
    for p in sys.argv[1:]:
        analyze(p)
