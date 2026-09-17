#!/usr/bin/env python3
"""
Heuristic model for A335406(6): where does the first run of six equal prime
gaps (the smallest CPAP-7) begin, how large would a(6) be, and what would it
cost to compute?

The expected number of maximal runs of exactly L equal gaps d that begin at a
prime <= x is modelled as

    N_{L,d}(x) = integral_1000^x rho_{L,d}(t) dt,

with rho from a Cramer/Hardy-Littlewood type random model: integers coprime to
W = 2*3*5*7*11*13 are "prime" independently with probability (W/phi(W))/ln t.
For every residue class of the start t (mod W) the model gives exactly
  (i)   the probability that the L+1 members t, t+d, ..., t+Ld are all prime,
  (ii)  the probability that the L*(d-1) integers between them are composite,
  (iii) the probability that the gaps just before and after the progression
        are not d, so that the run is maximal.
Primes q > 13 enter through the Hardy-Littlewood singular-series factor for the
L+1 members and through a per-position correction of the interior primality
probabilities (an interior number congruent to a member mod q cannot be
divisible by q, the others are slightly more likely to be).

The model is checked against the exact counts written by
'a335406 scan 1e13 -m pow10 -H' for many (L, d) pairs -- including L = 1, 2
with d = 210, which have the same 209-number interiors as a CPAP-7.  It is
exact to about 1% for small gaps but over-predicts runs of large gaps, more so
the larger u = d/ln x is.  The deficit is regular: for single gaps
ln(observed/model) = -F(u)/ln x, and every further link of a chain costs
another factor exp(-G(u)/ln x).  F and G are fitted to the data (power laws in
u) and applied to L = 6, d = 210: six interiors and five extra links.  Both the
raw and the calibrated predictions are printed; the calibration changes the
density by a factor of order 2-3, which moves the predicted location of the
first CPAP-7 by much less than a decade because the expected count grows about
thirtyfold per decade there.

usage: model_a335406.py SCANFILE [RATE]
    SCANFILE  output of 'a335406 scan ... -m pow10 -H'
    RATE      sieve throughput in numbers/s (default 1e10), for the time estimates
"""
import math
import re
import sys

import numpy as np

SMALL = [2, 3, 5, 7, 11, 13]
W = 30030
WPHI = W / math.prod(q - 1 for q in SMALL)          # W/phi(W)
QMAX = 200000


def primes_upto(n):
    s = np.ones(n + 1, dtype=np.bool_)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)


BIGQ = [int(q) for q in primes_upto(QMAX) if q > 13]


def li(x):
    """logarithmic integral, li(x) = gamma + ln ln x + sum (ln x)^n / (n n!)"""
    u = math.log(x)
    s, term = 0.0, 1.0
    for n in range(1, 400):
        term *= u / n
        s += term / n
        if term / n < 1e-17 * s and n > u:
            break
    return 0.5772156649015329 + math.log(u) + s


class Model:
    def __init__(self, L, d):
        self.L, self.d = L, d
        members = [j * d for j in range(L + 1)]
        interior = np.array([j * d + r for j in range(L) for r in range(1, d)], dtype=np.int64)
        T = 1.0                                          # singular series, primes q > 13
        logc = np.zeros(len(interior))                   # interior corrections, primes q > 13
        for q in BIGQ:
            res = sorted(set(m % q for m in members))
            nu = len(res)
            if nu >= q:
                T = 0.0
                break
            T *= (1 - nu / q) / (1 - 1 / q) ** (L + 1)
            protected = np.isin(interior % q, res)
            logc += np.where(protected, -math.log1p(-1 / q),
                             math.log1p(-1 / (q - nu)) - math.log1p(-1 / q))
        self.T = T
        self.c = np.exp(logc)
        a = np.arange(W)
        ok = np.ones(W, dtype=np.bool_)
        for m in members:
            ok &= np.gcd(a + m, W) == 1
        cl = a[ok]
        self.nclasses = len(cl)
        self.cop = (np.gcd(cl[:, None] + interior[None, :], W) == 1).astype(np.float64)
        self.left_ok = (np.gcd(cl - d, W) == 1).astype(np.float64)
        self.left_n = (np.gcd(cl[:, None] + np.arange(-d + 1, 0)[None, :], W) == 1).sum(1)
        self.right_ok = (np.gcd(cl + (L + 1) * d, W) == 1).astype(np.float64)
        self.right_n = (np.gcd(cl[:, None] + np.arange(L * d + 1, (L + 1) * d)[None, :], W) == 1).sum(1)

    def rho(self, t):
        """probability that the integer t begins a maximal run of L gaps d"""
        if self.T == 0.0 or self.nclasses == 0:
            return 0.0
        pe = WPHI / math.log(t)                          # prime probability of a number coprime to W
        interior = np.exp(self.cop @ np.log1p(-pe * self.c))
        left = self.left_ok * pe * np.exp(self.left_n * math.log1p(-pe))
        right = self.right_ok * pe * np.exp(self.right_n * math.log1p(-pe))
        per_class = pe ** (self.L + 1) * self.T * interior * (1 - left) * (1 - right)
        return float(per_class.sum()) / W

    def N(self, xs, npts=1500):
        """expected number of runs beginning at primes <= x for each x in xs"""
        xs = np.asarray(xs, dtype=np.float64)
        u0, u1 = math.log(1000.0), math.log(xs.max())
        us = np.linspace(u0, u1, npts)
        f = np.array([self.rho(math.exp(u)) * math.exp(u) for u in us])
        cum = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(us))))
        return np.interp(np.log(xs), us, cum)


def parse_scan(path):
    """cumulative totals per mark: pos -> dict(pi, E, runs, ld{(L,d): count})"""
    marks = {}
    cur = None
    for line in open(path):
        m = re.match(r"--- (?:cumulative )?totals for the primes (?:<= |in \[[\d,]+, )([\d,]+)", line)
        if m:
            cur = int(m.group(1).replace(",", "").rstrip("]"))
            marks[cur] = {"ld": {}}
            continue
        if cur is None:
            continue
        m = re.match(r"primes <= [\d,]+ \(pi\)\s+([\d,]+)", line)
        if m:
            marks[cur]["pi"] = int(m.group(1).replace(",", ""))
        m = re.match(r"equal adjacent gaps \(E\)\s+([\d,]+)", line)
        if m:
            marks[cur]["E"] = int(m.group(1).replace(",", ""))
        m = re.match(r"ld (\d+) (\d+) (\d+) (\d+)", line)
        if m:
            marks[cur]["ld"][(int(m.group(2)), int(m.group(3)))] = int(m.group(4))
    return dict(sorted(marks.items()))


def fmt(x):
    return f"{x:.3g}" if x < 1e5 else f"{x:.3e}"


def fit_power(us, fs):
    """least squares fit f = c u^g in log-log; returns (c, g)"""
    us, fs = np.asarray(us), np.asarray(fs)
    A = np.vstack([np.ones(len(us)), np.log(us)]).T
    (lc, g), *_ = np.linalg.lstsq(A, np.log(fs), rcond=None)
    return math.exp(lc), g


def calibrate(marks, xs):
    """
    Fit the deficit of the raw model.  Single gaps (L = 1):
    F(u) = -ln(obs/model) * ln x, u = d/ln x.  Chains (L >= 2): each link beyond
    the first costs exp(-G(u)/ln x) more, G(u) = -ln((obs/model) / r1^L) * ln x / (L-1).
    Returns the fitted (F, G) as functions of u, plus the raw points.
    """
    Fu, Ff, Gu, Gf = [], [], [], []
    r1 = {}
    ds = list(range(6, 61, 6)) + list(range(66, 301, 6))
    for d in ds:
        pred = Model(1, d).N(xs, npts=600)
        for x, p in zip(xs, pred):
            obs = marks[int(x)]["ld"].get((1, d), 0)
            if obs < 1000 or p <= 0:
                continue
            u = d / math.log(x)
            r1[(d, x)] = obs / p
            if u >= 1.0:
                Fu.append(u)
                Ff.append(-math.log(obs / p) * math.log(x))
    for L in (2, 3, 4, 5):
        for d in (12, 18, 24, 30, 36, 42, 48, 60, 90, 120, 150, 180, 210):
            pred = Model(L, d).N(xs, npts=600)
            for x, p in zip(xs, pred):
                obs = marks[int(x)]["ld"].get((L, d), 0)
                if obs < 100 or p <= 0 or (d, x) not in r1:
                    continue
                e = (obs / p) / r1[(d, x)] ** L
                if e <= 0 or e >= 1.5:
                    continue
                Gu.append(d / math.log(x))
                Gf.append(-math.log(e) * math.log(x) / (L - 1))
    cF, gF = fit_power(Fu, Ff)
    keep = [i for i, g in enumerate(Gf) if g > 0]
    cG, gG = fit_power([Gu[i] for i in keep], [Gf[i] for i in keep])
    return (lambda u: cF * u ** gF), (lambda u: cG * u ** gG), (cF, gF, len(Fu)), (cG, gG, len(keep))


class Calibrated:
    """raw model density times exp(-(L F(u) + (L-1) G(u)) / ln t), u = d / ln t"""

    def __init__(self, model, F, G):
        self.m, self.F, self.G = model, F, G

    def N(self, xs, npts=1500):
        xs = np.asarray(xs, dtype=np.float64)
        u0, u1 = math.log(1000.0), math.log(xs.max())
        us = np.linspace(u0, u1, npts)
        f = np.empty(len(us))
        for k, uu in enumerate(us):
            t = math.exp(uu)
            u = self.m.d / uu
            corr = math.exp(-(self.m.L * self.F(u) + (self.m.L - 1) * self.G(u)) / uu)
            f[k] = self.m.rho(t) * t * corr
        cum = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(us))))
        return np.interp(np.log(xs), us, cum)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    marks = parse_scan(sys.argv[1])
    rate = float(sys.argv[2]) if len(sys.argv) > 2 else 1e10
    xs = np.array(sorted(marks), dtype=np.float64)
    if len(xs) == 0:
        sys.exit("no totals found in " + sys.argv[1])

    print("Model check: observed runs of exactly L gaps d beginning at primes <= x, vs the model")
    print("(L = 1, 2 with d = 210 have the same 209-number interiors as a CPAP-7)\n")
    pairs = [(1, 6), (1, 30), (1, 210), (2, 6), (2, 30), (2, 210), (3, 6), (3, 12), (3, 30),
             (4, 30), (4, 60), (5, 30), (5, 60), (2, 420), (1, 420)]
    print(f"{'L':>2} {'d':>4} " + " ".join(f"{'x=%.0e' % x:>16}" for x in xs))
    for L, d in pairs:
        mod = Model(L, d)
        pred = mod.N(xs)
        cells = []
        for x, p in zip(xs, pred):
            obs = marks[int(x)]["ld"].get((L, d), 0)
            cells.append(f"{obs:>8d}/{p:>7.3g}" if obs < 1e6 else f"{obs:>8.3g}/{p:>7.3g}")
        print(f"{L:>2} {d:>4} " + " ".join(f"{c:>16}" for c in cells))
        ratios = [marks[int(x)]["ld"].get((L, d), 0) / p if p > 0 else float('nan') for x, p in zip(xs, pred)]
        print(f"{'':>7} " + " ".join(f"{'obs/pred %.3f' % r:>16}" for r in ratios))
    print("\n(each cell is observed/predicted)\n")

    # calibration of the deficit
    xcal = xs[xs >= 1e9]
    F, G, (cF, gF, nF), (cG, gG, nG) = calibrate(marks, xcal)
    print(f"Calibration: single gaps  ln(obs/model) = -F(u)/ln x,  F(u) = {cF:.3f} u^{gF:.3f}   ({nF} points)")
    print(f"             chain links  extra factor exp(-G(u)/ln x) per link,  G(u) = {cG:.3f} u^{gG:.3f}   ({nG} points)")
    print("Calibrated model vs observed for the longer chains (each cell observed/calibrated, ratio):")
    for L, d in [(3, 30), (4, 30), (5, 30), (3, 60), (4, 60), (5, 60), (2, 120), (3, 120), (2, 150), (2, 180), (1, 210), (2, 210), (3, 210)]:
        cal = Calibrated(Model(L, d), F, G).N(xcal, npts=600)
        cells = []
        for x, p in zip(xcal, cal):
            obs = marks[int(x)]["ld"].get((L, d), 0)
            cells.append(f"{obs:>7d}/{p:>7.3g} {obs / p if p > 0 else float('nan'):5.2f}")
        print(f"{L:>2} {d:>4} " + " ".join(f"{c:>22}" for c in cells))
    print()

    # fraction of primes whose gap equals the previous gap:  E/pi ~ alpha/ln x + beta/ln^2 x
    xf = xs[xs >= 1e8]                                   # fit only where the asymptotics have set in
    lx = np.log(xf)
    frac = np.array([marks[int(x)]["E"] / marks[int(x)]["pi"] for x in xf])
    A = np.vstack([1 / lx, 1 / lx ** 2]).T
    (alpha, beta), *_ = np.linalg.lstsq(A, frac, rcond=None)
    print("Fraction E/pi of primes whose gap equals the previous gap (fit uses x >= 1e8):")
    for x, f in zip(xf, frac):
        print(f"   x = {x:.0e}: {100 * f:.4f}%   fit {100 * (alpha / math.log(x) + beta / math.log(x) ** 2):.4f}%")
    print(f"   fit E/pi = {alpha:.4f}/ln x + {beta:.4f}/ln^2 x\n")

    # L = 6: the first CPAP-7
    grid = np.logspace(15, 26, 221)
    m210, m420 = Model(6, 210), Model(6, 420)
    n210, n420 = m210.N(grid), m420.N(grid)
    c210, c420 = Calibrated(m210, F, G).N(grid), Calibrated(m420, F, G).N(grid)
    total_raw, total = n210 + n420, c210 + c420
    print("Runs of six equal gaps (maximal CPAP-7) beginning at primes <= x, expected number:")
    print(f"{'x':>8} {'raw d=210':>10} {'raw d=420':>10} {'calibrated':>11} {'P(first CPAP-7 < x)':>22}")
    for e in range(16, 26):
        x = 10.0 ** e
        i = int(round((e - 15) * 20))
        print(f"{'1e%d' % e:>8} {n210[i]:>10.3g} {n420[i]:>10.3g} {total[i]:>11.3g} {1 - math.exp(-total[i]):>22.4f}")
    qs, qs_raw = {}, {}
    for name, level in (("10%", 0.10), ("50%", 0.50), ("90%", 0.90)):
        target = -math.log(1 - level)
        qs[name] = float(np.interp(target, total, grid))
        qs_raw[name] = float(np.interp(target, total_raw, grid))
    print("\nFirst CPAP-7 start P, quantiles of the raw model:        " +
          ", ".join(f"{k}: {v:.2e}" for k, v in qs_raw.items()))
    print("First CPAP-7 start P, quantiles of the calibrated model: " +
          ", ".join(f"{k}: {v:.2e}" for k, v in qs.items()))
    zim = 71137654873189893604531
    print(f"Smallest known CPAP-7 (P. Zimmermann): {zim} = {zim:.3e}; expected number of CPAP-7 "
          f"below it: {float(np.interp(math.log10(zim), np.log10(grid), total)):.3g} (calibrated), "
          f"{float(np.interp(math.log10(zim), np.log10(grid), total_raw)):.3g} (raw)")

    print("\nImplied a(6) = runs below P = pi(P) - E(P) ~ li(P) (1 - E/pi), and the cost of sieving to P:")
    print(f"{'P':>10} {'pi(P) ~ li(P)':>16} {'a(6) estimate':>16} {'sieve time at %.2g numbers/s' % rate:>32}")
    for name, P in list(qs.items()) + [("known", float(zim))]:
        L = li(P)
        f = alpha / math.log(P) + beta / math.log(P) ** 2
        secs = P / rate
        print(f"{name + ' ' + '%.2e' % P:>10} {L:>16.4e} {L * (1 - f):>16.4e} "
              f"{secs / 3.15576e7:>22.3g} years")
    print("\nFor comparison, with the same rate:")
    for e in (13, 14, 15, 16, 18):
        secs = 10.0 ** e / rate
        print(f"   sieving to 1e{e}: {secs / 3600:.3g} h = {secs / 86400:.3g} days = {secs / 3.15576e7:.3g} years")


if __name__ == "__main__":
    main()
