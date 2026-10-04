# A002966 — can a(9) be computed?

A002966(n) = number of solutions of 1 = 1/x_1 + ... + 1/x_n with x_1 <= ... <= x_n.
Known: 1, 1, 3, 14, 147, 3462, 294314, 159330691 (a(8) by Dethridge & Le Normand, Jan 2004).
Siblings: A006585 (distinct denominators), A002967 (ordered), A156869 / A280520 (triangles whose
next term is a(9)). Nobody has computed a(9); Crawford's 2019 VT master's thesis
(hdl.handle.net/10919/90573) estimates a(9) ~ 3.0e12 - 8.5e12 and calls it infeasible.

## Verdict (2026-10-02): not feasible with known methods on our hardware

* Every exact method reduces to a scan over all "heads" (x_1..x_7 prefixes) with a
  divisor-class count for the last two terms. The head count for n=9 is
  **713,922,875,728,745** (reproduced exactly by `a002966_heads`, matching the thesis).
* Measured cost of the best per-head routine I built (`a002966_mitm`, meet-in-the-middle
  over divisors of (qx)^2): 6-25 us/head unoptimized on the M1 Pro, for heads of every size
  (q is ~1e12-1e13 for essentially all n=9 heads, so tau(q^2) is ~1e3 everywhere).
  Montgomery arithmetic might reach 2-4 us/head. That is **50-70 core-years**; the 96-core
  box plus Mac Studio would need 4-7 months, and a CUDA port is unproven for this kernel
  (hash lookups, 128-bit intervals, batch modular inverses).
* Sub-linear ideas I tried all leave a region that still needs the full scan:
  - quotient parametrization k = (u + p w)/a: covers small k, but tiny-d / huge-z
    solutions need per-head candidate checks;
  - gcd structure: ab - q^2 divides (pq gcd(x,y))^2, so x/gcd(x,y) | p gcd(x,z); the
    solutions are spread evenly over all scales of x/gcd(x,y) (measured, `a002966_gstat`);
  - small-d candidates d = (-qx mod a) + j a: the candidate count D/a is spread over 2^0..2^35;
  - factoring p u w' + q^2 (~1e26) per (u, w') is too slow for w' beyond a handful.
* Solutions per head are spread evenly over all dyadic ranges of a = p x - q (measured on
  n=8 heads with `a002966_tail3`), ~90% in the "big-w" regime, so no single regime can be
  special-cased away.

## Tools (all single-file C, build with `make`)

| file | purpose |
|---|---|
| `a002966_heads.c` | enumerate (n-3)-prefixes, count heads by x_{n-2} range (`./a002966_heads 9`, 17 s) |
| `a002966_baseline.c` | classic method (Le Normand 2004 style), reproduces a(3..8); a(8) in 30 min single-thread |
| `a002966_mitm.c` | faster per-head count (reduced fraction, h = gcd(x,q) structure, MITM); reproduces a(3..8) (a(8) in 462 s, 3.0 us/head, M1 Pro single thread); `./a002966_mitm N`, or `./a002966_mitm head p q xlo xhi` to time one head; `./a002966_mitm N dbg` cross-checks every head against brute force |
| `a002966_tail3.c` | histogram of 3-term tail solutions by a = p x - q for one head |
| `a002966_gstat.c` | histogram of solutions by x/gcd(x,y) and by D/a for one head |
| `heads9_range_ge_1e8.txt` | the 424,004 n=9 prefixes whose x_7 range is >= 1e8 (p, q, s, range) |

Key correctness detail (cost me a bug): the tail fraction a/(qx) must be reduced to lowest
terms a'/M before counting divisors d of M^2 with d ≡ -M (mod a'), (a' x - M) <= d <= M;
counting in the unreduced form overcounts (a(5) came out 152/155 instead of 147).

## Measured timings (M1 Pro, one core, `a002966_mitm head`)

| head (p, q) | q/p | region | us/x | solutions / 1e5 x |
|---|---|---|---|---|
| 1, 10650056950806 (Sylvester) | 1.07e13 | x = q+1.. | 60 | 7.2e6 |
| same | | a ~ q/2 | 12 | 1087 |
| same | | a ~ q | 12 | 1071 |
| same | | a ~ 1.8q | 12 | 0 |
| 37979, 1898951211612 | 5.0e7 | dense / sparse | 10 / 15 | 1696 / 0 |
| 5360, 2680005470961 | 5.0e8 | dense / sparse | 6 / 9 | 900 / 0 |
| 2131, 10657008082266 | 5.0e9 | dense / sparse | 7 / 10 | 4129 / 0 |
| 71, 3550249600170 | 5.0e10 | dense / sparse | 20 / 16 | 757227 / 20 |
| 1, 507148677126 | 5.1e11 | dense / sparse | 554 / 23 | 1.3e8 / 13400 |

## If someone wants to try anyway

The only route is a CUDA port of `count_x` in `a002966_mitm.c` (per-thread: sieve-factored x,
residue chains mod a', batch inversion, hash match) on both Sparks, with the 2.46e8 prefixes
as work units ordered by range. Expect months, and plan a second independent run for
verification. A genuinely new idea would be a sub-linear count of 1/x+1/y+1/z = p/q for
q/p ~ 1e13; I found none in the literature (Elsholtz-Tao, Browning-Elsholtz, Levaillant 2025,
Ghirlanda 2025 all stop at structure/asymptotics).
