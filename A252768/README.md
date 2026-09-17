# A252768 — primes whose prime-gap power sums are all prime

Tool for computing and extending [OEIS A252768](https://oeis.org/A252768):

> Primes p with property that the sum of the k-th powers of the successive gaps
> between primes <= p are prime numbers for k = 1 to n.

Take the primes 2 = q₁ < q₂ < … < qₘ = p and their gaps gᵢ = qᵢ₊₁ − qᵢ, and let

    S_k(p) = Σ gᵢ^k        (k = 1, 2, …)

The *depth* of p is the largest n such that S₁(p), …, Sₙ(p) are all prime, and
a(n) is the smallest prime of depth ≥ n. Example: the primes up to 13 have gaps
1, 2, 2, 4, 2, so S₁ = 11, S₂ = 29, S₃ = 89 (all prime) and S₄ = 305 = 5·61, so
13 has depth 3 and a(3) = 13.

Since S₁(p) = p − 2, every term is the larger member of a twin prime pair
(A006512), and a(n) is non-decreasing in n.

The OEIS entry (Sep 2026, keyword `hard,more`) lists a(1..7):
5, 5, 13, 14593, 372313, 2315773, 541613713.

## Results

Everything below was computed with this tool on an Apple M1 Pro (10 cores, Sep 16 2026).
Both new terms were verified independently: by `a252768 verify` (a separate
single-threaded code path, every verdict cross-checked with GMP) and by
`verify_a252768.py` (numpy sieve, exact Python integers, its own Miller–Rabin),
which reproduce every power sum digit for digit.

**Two new terms:**

| n | a(n) | S₁..Sₙ prime, S_(n+1) composite |
|---|------|------|
| 8 | **7952072743** | S₉ = 38009162738943211881748481 = composite |
| 9 | **21814967833** | S₁₀ = 30245374383129007149396784129 = composite |

The exhaustive scan of all primes below 10¹² (37,607,912,018 primes, 111 s)
found no prime of depth ≥ 10, so **a(10) > 10¹²** (and a(11), a(12) > 10¹²;
sums were tested up to k = 12). The OEIS terms a(1..7) are reproduced by the
same scan. The sequence is now

```
5, 5, 13, 14593, 372313, 2315773, 541613713, 7952072743, 21814967833
```

The power sums at the new terms (all prime unless marked):

```
p = 7952072743   S_1 = 7952072741            S_2 = 306387204649           S_3 = 17232652308257
                 S_4 = 1275921083018209      S_5 = 117138327618559361     S_6 = 12826308390187298689
                 S_7 = 1629701772261995408897                 S_8 = 235364125518401969499649
                 S_9 = 38009162738943211881748481 (composite)

p = 21814967833  S_1 = 21814967831           S_2 = 883642617829           S_3 = 52322085406601
                 S_4 = 4080850722829393      S_5 = 394810145761684001     S_6 = 45578875566881080129
                 S_7 = 6111338676278533840001                 S_8 = 932919209682371038708993
                 S_9 = 159652836113962880019857921            S_10 = 30245374383129007149396784129 (composite)
```

### Depth statistics below 10¹²

Number of primes p < 10¹² of exact depth n (S₁..Sₙ prime, S_(n+1) not):

| n | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| count | 1747205024 | 112645225 | 9907858 | 756501 | 65240 | 4934 | 401 | 32 | 5 |

Each extra k costs a factor of only 8–13 (the conditions are positively
correlated, see below). The five primes of depth 9 below 10¹² are 21814967833,
60308763733, 73668807349, 331274803549 and 844836146389; all 438 primes of
depth ≥ 7 are listed in `scan_1e12.txt`.

### Where a(10) should be

About one in ten primes of depth 9 reaches depth 10, and the density of
depth-9 primes (3 below 10¹¹, 2 more below 10¹²) falls by roughly half per
decade, so the expected number of depth-10 primes is about 1.4 in
[10¹², 10¹³] and about 7 in [10¹³, 10¹⁴]: a(10) is probably below 10¹³
(~75 %) and almost certainly below 10¹⁴. a(11) is expected another factor of
~10 further out, i.e. around 10¹⁴–10¹⁵.

### Timings (M1 Pro, 10 threads)

| scan | time | rate |
|------|------|------|
| 0 .. 6·10⁸ (selftest) | 0.1 s | |
| 0 .. 10¹⁰ | 1.2 s | 8.0·10⁹ numbers/s |
| 0 .. 10¹¹ | 10.5 s | 9.5·10⁹ numbers/s |
| 0 .. 10¹² | 111 s | 9.0·10⁹ numbers/s |
| 0 .. 10¹³ (estimate) | ~20 min | |
| 0 .. 10¹⁴ (estimate) | ~3 h | |
| 0 .. 10¹⁵ (estimate) | ~1.3 days | |

primesieve alone counts primes 3–8× faster than this (1.4 s per 10¹¹ near 0,
4 s per 10¹¹ near 10¹⁵), so the tool is bound by its per-prime work — mostly
the S₂ tests of the 1.87·10⁹ twin primes below 10¹² — which shrinks relative
to the range as the primes thin out; the rate should stay near 9·10⁹/s up to
10¹⁵. `verify` (single-threaded, from scratch) takes 5 s for 2.2·10¹⁰.

## Building

Needs [primesieve](https://github.com/kimwalisch/primesieve) (`brew install primesieve`).
GMP (`brew install gmp`) is optional; when present the `verify` command and
the selftest cross-check every primality verdict with `mpz_probab_prime_p`.

```
make          # builds ./a252768 (adds -DHAVE_GMP -lgmp when gmp.h is found)
make test     # runs the selftest (about 1 s)
```

or by hand:

```
cc -O2 -std=gnu11 -pthread -I/opt/homebrew/include a252768.c \
   -L/opt/homebrew/lib -lprimesieve -lm -DHAVE_GMP -lgmp -o a252768
```

## Usage

```
a252768 scan [START] END [-t T] [-c CHUNK] [-k KMAX] [-r DEPTH] [-S STATE] [-i SECS] [-q]
a252768 verify P [-k KMAX] [-t T]
a252768 selftest [-t T]
```

Numbers may be written as `123`, `2^40`, `10^12`, `1e12` or `2^40-1`.

**scan** searches every prime in [START, END] (START defaults to 0). It prints
each prime of depth ≥ DEPTH (default 7) as it is found, with all its power
sums and the first composite one, then a table of a(n) for n ≤ KMAX (default
12, at most 16). With START > 0 the power sums up to START are first computed
in a fast totals-only pass, so any sub-range can be scanned on its own.

`-S FILE` writes a checkpoint every SECS seconds (default 60) and on Ctrl-C;
running the *same* command again resumes where it left off, and END may be
raised on resume. `-t` sets the thread count (default: all cores), `-c` the
chunk size (default: range/(32·threads), clamped to 10⁶..10⁹), `-q` silences
the status line on stderr.

```
./a252768 scan 1e12 -S a252768.state          # 2 minutes on an M1 Pro
./a252768 scan 1e13 -S a252768.state          # continues from the checkpoint
./a252768 scan 5e11 6e11 -r 8                 # only this sub-range, print depth >= 8
```

**verify P** recomputes S₁..S_KMAX at the prime P from scratch with a plain
single-threaded accumulation (no chunks, no offsets) and tests each sum, with
GMP as an independent second opinion. `-t T` with T > 1 uses the parallel
totals pipeline instead. With GMP, KMAX may exceed 16 here.

```
./a252768 verify 21814967833 -k 12
```

**selftest** checks the Miller–Rabin code against known (pseudo)primes and
20000 + 3000 random numbers versus GMP, recomputes the OEIS example and a(7),
reproduces a(1..7) with a chunked scan, checks that different chunk sizes and a
mid-range start give identical hits and sums, and round-trips a checkpoint.

`verify_a252768.py P [KMAX]` is a fully independent Python check (numpy
segmented sieve, exact Python integers, its own Miller–Rabin); it shares no
code with the C tool and needs about 1 minute per 10¹⁰.

## How it works

### The sequential problem, made parallel

S_k(p) is a running sum over *all* gaps below p, so a naive scan is inherently
sequential. The tool cuts the range into chunks and uses the fact that the sum
over a chunk depends only on the chunk's own gaps:

1. A worker takes the next chunk index from an atomic counter, sieves the
   chunk with primesieve (starting a little below the chunk to find the prime
   before it, so every prime's incoming gap is attributed to the chunk that
   contains the prime), and stores the gaps — one byte per prime (g/2, with
   an escape for the gap 1 and gaps ≥ 510) — together with a gap histogram.
2. From the histogram it forms the chunk totals T_k = Σ count(g)·g^k for
   k ≤ KMAX in 128-bit arithmetic.
3. **Ordered hand-off**: the worker waits until the running sums O_k at the
   start of its chunk have been published (chunk i publishes O_k(i+1) =
   O_k(i) + T_k(i) as soon as its *sieve* is done), then proceeds.
4. It walks its stored gaps a second time with the correct starting sums and
   tests candidates. Meanwhile the other workers do the same with their
   chunks, so the checking runs fully in parallel, the range is sieved exactly
   once, no barrier is needed and memory stays bounded (one chunk of gaps per
   thread; 10⁹ numbers ≈ 40 MB).

The bookkeeping ring is 4096 chunks deep; a completion frontier tracks the
lowest unfinished chunk, which is where checkpoints are written and where
`a(n) > X` statements are certified.

### Checking a chunk

S₁(p) = p − 2 is prime exactly when the gap into p is 2, so only twin primes
(about 1/ln p of all primes) are examined at all. S₂, S₃, S₄ are maintained
incrementally and tested in turn; the ~1 in 10⁵ primes surviving S₂..S₄ get
S₅, S₆, … evaluated from the running gap histogram (Σ count(g)·g^k over the
few hundred distinct gap sizes), so the cost per prime is a handful of adds.
Every sum is 128-bit with overflow detection; if some S_k would exceed 2¹²⁸
(S₁₂ does so around 10¹⁵), depths ≥ k stop being tested and the tool says so.

### Primality

* n < 2⁶⁴: Miller–Rabin to the bases 2, 325, 9375, 28178, 450775, 9780504,
  1795265022, which is deterministic for all 64-bit n. Trial division by the
  primes up to 47 first (two 32-bit remainders) removes two thirds of the
  candidates before the expensive part. Everything runs in Montgomery
  arithmetic; the base-2 exponentiation multiplies by the base with a modular
  doubling.
* n ≥ 2⁶⁴: 128-bit Montgomery Miller–Rabin to the 24 prime bases 2..89. The
  first 13 bases are deterministic below 3.3·10²⁴ (Sorenson–Webster); above
  that the verdict is a strong probable-prime test, and `verify` adds GMP's
  independent BPSW + Miller–Rabin verdict. The sums that matter here have
  at most ~110 bits.

### Why the depths are correlated

By Fermat, g^k mod q depends only on k mod (q − 1), so S_k ≡ S_k′ (mod q)
whenever k ≡ k′ (mod q − 1). Once S₁ and S₂ are prime, no later S_k is
divisible by 3 (odd k ≡ S₁, even k ≡ S₂ mod 3), S₅ ≡ S₁ and S₆ ≡ S₂ (mod 5),
S₇ ≡ S₁ (mod 7), and so on. This is why a(n) grows more slowly than
independent conditions would suggest: empirically each extra k costs a factor
of only about 8–12 in the count of qualifying primes (see the depth table
above), and S_k(p) ~ k!·p·(ln p)^(k−1) makes the per-k prime probability
~2·c_k/ln S_k(p) shrink only slowly with p.

## Verification

* `a252768 selftest` reproduces a(1..7) and cross-checks 23000 random
  Miller–Rabin verdicts against GMP.
* `a252768 verify P` recomputes the sums single-threaded from scratch (a code
  path without chunks or offsets) and shows the GMP verdict next to the
  tool's own for every k.
* `verify_a252768.py P` repeats the whole computation without primesieve or
  any of the C code.
* Inside every chunk the tool asserts that the published running sum S₁ equals
  (previous prime − 2), which pins the chunk hand-off to the actual primes.

Primality of the sums above 3.3·10²⁴ (S₈ and beyond for the terms found so
far) is probable-prime status from 24 Miller–Rabin bases plus GMP's BPSW test;
for an OEIS submission a certificate (e.g. PARI/GP `isprime`, which is APR-CL)
would make it rigorous. All the sums involved have at most ~110 bits.

## Files

* `a252768.c` — the tool (single file, C11, pthreads, primesieve, optional GMP)
* `Makefile` — `make`, `make test`
* `verify_a252768.py` — independent Python/numpy verification of one prime
* `scan_1e12.txt`, `scan_1e12.log` — output of the scan to 10¹² (all primes of depth ≥ 7)
* `a252768.state` — checkpoint of the scan (resume with the same `scan … -S a252768.state`)
* `b252768.txt`, `DATA.txt` — b-file and OEIS DATA line with the new terms
