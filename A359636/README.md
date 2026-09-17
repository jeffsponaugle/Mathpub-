# A359636 — prime gaps whose composites all have ≥ n distinct prime factors

Tool for computing and extending [OEIS A359636](https://oeis.org/A359636):

> a(n) is the least odd prime not in A001359 such that all subsequent composites
> in the gap up to the next prime have at least n distinct prime factors.

Take consecutive primes p < q with q − p ≥ 4 (so p is not the lesser member of a
twin prime pair, A001359). The gap *qualifies at level n* when ω(x) ≥ n for every
x in p+1 … q−1, where ω = A001221 counts distinct prime factors. a(n) is the
smallest such p.

The OEIS entry (Sep 2026, keyword `hard,more`) lists a(1..8):

```
7, 19, 643, 51427, 8083633, 1077940147, 75582271489, 34710483181813
```

with the bounds a(9) ≤ 76340177205657727 and a(10) ≤ 225096507194749219819
(David A. Corneth). Every known term and both bounds are prime gaps of length 4,
for example a(8) = 34710483181813:

```
34710483181814 = 2*11*17*29*47*229*409*727
34710483181815 = 3^2*5*7*19*31*101*263*7043
34710483181816 = 2^3*13*41*43*61*83*139*269      34710483181817 prime
```

## Results

Computed with this tool on an Apple M1 Pro (10 cores, Sep 16 2026).

* `a359636 selftest` reproduces a(1..7) in about two seconds and a(8) in
  10 min 56 s (exhaustive scan of all m ≤ 3.47·10¹³ at 53 Gm/s on 10 threads: 142,705,091 survivors, 7 triples, 1 qualifying gap; `repro_a8.txt`).
* `verify_a359636.py` (independent Python: numpy ω-sieve over every integer,
  every prime gap walked directly, own Miller–Rabin and Pollard rho) recomputes
  a(1..6) from the definition and re-verifies each known gap.
* Level 9, **new upper bound: a(9) ≤ 33955545649252303**, found by the
  structured hunts (see below) and verified independently by
  `verify_a359636.py`. This is 45% of Corneth's bound 76340177205657727.
  Its gap:

  ```
  33955545649252303 prime
  33955545649252304 = 2^4*7^2*11*17*29*37*163*1021*1297
  33955545649252305 = 3*5*13*19*61*71*127*251*66383
  33955545649252306 = 2*23*41*47*53*113*181*313*1129
  33955545649252307 prime
  ```

  The first hunt (`-L 1000`, m with all prime factors ≤ 1000) found six
  qualifying gaps below Corneth's bound, the smallest at p = 36189117287569243
  (m = 3·5·37·83·107·163·197·373·613); the free-largest-prime hunt
  (`-L 300 -Q 2e5`) then found the smaller one above, whose m carries the prime
  66383, the same shape as Corneth's solution.

  Whether it *is* a(9) is being settled by the exhaustive scan
  (`scan_n9.txt`/`.log`), which stops by itself once its frontier passes the
  smallest solution it finds; N9_STATUS

For comparison, the numbers of *triples* m ≡ 3 (mod 6) with ω(m−1), ω(m), ω(m+1)
all ≥ n found below a(n) (each is a potential term; it becomes one only when the
surrounding gap is bounded by primes and every composite in it qualifies):

| n | a(n) | triples with m < a(n)+2 | of which qualifying gaps |
|---|------|------|------|
| 5 | 8083633 | 10 | 1 |
| 6 | 1077940147 | 16 | 1 |
| 7 | 75582271489 | 3 | 1 |
| 8 | 34710483181813 | 7 | 1 |

## Method

**Reduction.** Every qualifying gap contains an odd multiple of 3, call it m,
such that m−1, m, m+1 are all composites of the gap: in a gap of 4 the prime
p is 1 (mod 3), so m = p+2; in a gap of 6, m = p+2 or p+4; a gap of 8 or more
holds three consecutive odd numbers, one of them 0 (mod 3). So it is enough to
find every m ≡ 3 (mod 6) with

    ω(m−1) ≥ n,  ω(m) ≥ n,  ω(m+1) ≥ n                          (*)

and then look at the prime gap around each of them. a(n) is the least
p = prevprime(m−1) over the m satisfying (*) whose whole gap qualifies.

**Three-target sieve.** Write m = 6k+3. For every prime 5 ≤ p ≤ T and each of
the three targets N = 6k+2, 6k+3, 6k+4, the k with p | N form one residue class
mod p. One byte counter per target, packed three to a 32-bit word per k, is
incremented along these classes on L1-sized segments (16384 k). The primes
5·7·11·13, 17·19·23 and 29·31·37 are applied as precomputed periodic patterns
(periods 5005, 7429, 33263); the remaining primes are marked with a strided
loop that handles the three targets of a prime together. The forced divisors
2 | m±1 and 3 | m are not counted; the thresholds are lowered by one instead.

**Why a small T is enough.** If N ≤ Bmax has n distinct prime factors, at most
rmax of them exceed T, where rmax is the largest r with

    (product of the first n−r primes) · (product of the first r primes above T) ≤ Bmax.

So N has c = n − rmax prime factors ≤ T, and a k whose three counters are all
≥ c−1 is a *survivor*. The tool picks T automatically as the smallest bound
that gives c ≥ 5; for n = 9 up to Corneth's bound this is T = 1000 (c = 5, i.e.
five prime factors below 1000 in each of m−1, m, m+1). Survivors occur at a
rate of about 2·10⁻⁵ per k near 7·10¹⁶, so the sieve — about 1.8–2.1 counter
increments per k — is essentially the whole cost.

**Exact verification.** Each survivor's three numbers are tested with trial
division using the same size bound as an early exit (a cofactor that still
needs r distinct primes but is smaller than the r-th power of the next prime
fails at once; a cofactor needing one more prime succeeds iff it exceeds 1; one
needing two is settled by Miller–Rabin plus a perfect-power test). For a k
satisfying (*), p = prevprime(m−1) and q = nextprime(m+1) are found with
deterministic Miller–Rabin and every composite in (p, q) is checked and printed
with its factorization.

**Hunt (upper bounds only).** `a359636 hunt N END -L L` enumerates every
m = 3^e · p₁⋯p_(N−1) (e = 1, 2; distinct primes 5 ≤ pᵢ ≤ L) with m ≤ END+2 and
tests m−1, m+1 and the gap exactly as above. It only sees the m whose own prime
factors are all small, so it can miss terms (a(8)'s m = 3²·5·7·19·31·101·263·7043
needs L ≥ 7043; a(6)'s needs L ≥ 1021), but for n = 9, L = 1000 it visits about
10¹¹ candidates below Corneth's bound instead of 1.3·10¹⁶ sieve positions, and
any hit is an upper bound for a(9) and a point where the exhaustive scan can stop.
Work is split over threads by the first two primes chosen.

With `-Q Q` the largest of the chosen primes may instead be any prime ≤ Q (this
is the shape of Corneth's solution, whose m = 3·7·13·17·23·137·163·211·151783).
Hunts run so far at level 9, all below the first hunt's bound 36189117287569243 unless noted:

| hunt | candidates | time | triples | qualifying gaps |
|---|---|---|---|---|
| `-L 1000` (below Corneth's 7.63·10¹⁶) | 1.47·10¹¹ | 2 h 12 min, 6 threads | 148 | 6, smallest p = 36189117287569243 |
| `-L 2000` | 2.22·10¹¹ | 3 h 39 min, 5 threads | 50 | none |
| `-L 300 -Q 2e5` | 1.69·10¹¹ | 3 h 08 min, 4 threads | 34 | 1, p = 33955545649252303 (m = 3·5·13·19·61·71·127·251·66383) |

**Parallel scan.** Chunks of 2²⁶ k are handed to worker threads through an
atomic counter; the lowest unfinished chunk (the *frontier*) drives the
progress line, the checkpoint file (`-S`, written every `-i` seconds and on
Ctrl-C) and the early stop: once the frontier has passed the smallest solution
found, that solution is a(n) and the scan ends (unless `-a` asks for all
solutions in the range).

## Performance and outlook

| range of m | k values | 10 threads (idle M1 Pro) |
|---|---|---|
| level 8 to a(8) = 3.47·10¹³ | 5.8·10¹² | 10 min 56 s |
| level 9 to 10¹⁶ | 1.7·10¹⁵ | ≈ 2 days |
| level 9 to Corneth's bound 7.63·10¹⁶ | 1.27·10¹⁶ | ≈ 15 days |

The sieve runs at ≈ 2.2 cycles per k on one performance core (≈ 57 Gm/s, i.e.
9.5·10⁹ k/s, on all ten cores). The ratios a(n+1)/a(n) so far are 80, 157,
133, 70 and 459, so a(9) is plausibly in the 10¹⁵–10¹⁷ range; whether the scan
finds it or only pushes a lower bound is a gamble on where it lies. a(10) is
beyond 64-bit arithmetic (its bound is 2.25·10²⁰) and out of reach for this
method.

## Usage

```
make                       # needs primesieve (brew install primesieve)
./a359636 selftest         # a(1..7); add -8 to include a(8) (~10 min)
./a359636 scan 9 0 76340177205657727 -t 10 -S scan_n9.state -i 120 > scan_n9.txt 2> scan_n9.log &
./a359636 scan 9 1e16 2e16 -a        # all solutions with 1e16 <= m <= 2e16+2
./a359636 hunt 9 76340177205657727 -L 1000   # structured upper-bound search (~1e11 candidates)
./a359636 hunt 9 36189117287569242 -L 300 -Q 2e5   # 7 primes <= 300 plus one free prime <= 2e5
./a359636 verify 9 76340177205657727 # check one gap, with factorizations
./a359636 omega 76340177205657729    # factor numbers / print omega
python3 verify_a359636.py selftest   # independent brute force a(1..5) + gap checks
python3 verify_a359636.py brute 6 1077940147
python3 verify_a359636.py check-log scan_n9.txt
```

Numbers may be written as decimal, `2^k`, `10^k`, `1e12`, or `X+Y` / `X-Y` of
those. `scan` searches m ≡ 3 (mod 6) with START ≤ m ≤ END+2, so every gap with
p ≤ END is covered. Options: `-t` threads, `-c` chunk size (k per chunk),
`-T` sieve prime bound (default automatic), `-S` state file, `-i` checkpoint
interval, `-a` scan the whole range instead of stopping at the first certain
a(n), `-q` quiet. Rerunning a `scan` with the same arguments and `-S` file
resumes from the checkpoint; the state file also records solutions found so far.

## Files

* `a359636.c`, `Makefile` — the tool (single C file, pthreads, primesieve).
* `verify_a359636.py` — independent Python cross-check (numpy for the brute force).
* `repro_a8.txt`, `repro_a8.log` — exhaustive level-8 scan reproducing a(8).
* `hunt_n9_L1000.txt`, `hunt_n9_L1000.log` — level-9 hunt (L = 1000) below Corneth's bound.
* `hunt_n9_L2000.txt`, `hunt_n9_L2000.log` — level-9 hunt (L = 2000) below 36189117287569243.
* `hunt_n9_L300_Q2e5.txt`, `.log` — level-9 hunt with a free largest prime (L = 300, Q = 2·10⁵).
* `OEIS_draft.md`, `b359636.txt` — draft submission text and b-file (pending confirmation).
* `scan_n9.txt`, `scan_n9.log`, `scan_n9.state` — the level-9 scan (output, progress, checkpoint).
