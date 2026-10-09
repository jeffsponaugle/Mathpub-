# OEIS submission candidates from the A007508 project

Cross-check against OEIS done 2026-10-03 with the JSON search API
(`https://oeis.org/search?fmt=json&q=id:A......`). "Last term" is the last
value in the DATA field at that date, with the extension line that added it.

Our tools: the CPU twin sieve (`a007508.c`, verified to a(15), works to 10^20)
and the GPU twin sieve (`gpu/a007508_cuda.cu`, GB10, verified to a(15)).
Measured GPU rates: 3.8e11 numbers/s at 10^15, 2.2e11 at 10^17, 1.8e11 at
10^18. The time estimates below are for one DGX Spark unless noted.

## 1. Submittable now, no computation (derived from existing terms)

### A167874 — primes < 10^n that belong to a twin prime pair
- Now: 18 terms, last a(18) = 1617351777154871 (a(16) corrected by Falk
  Hüffner, May 2026). Formula in the entry: a(n) = 2*A007508(n) - 1
  (5 is counted once). Holds for all 18 terms.
- Submit: **a(19) = 14475036187469089** = 2*7237518093734545 - 1, from
  A007508(19) (Benjamin Chaffin, corrected Jul 08 2026). State the source.

### A350826 — prime sextuplets with an n-digit initial term
- Now: 17 terms, last a(17) = 488096844 (Hugo Pfoertner, Jan 2022).
- A063501 (sextuplets with largest member <= 10^n) now reaches a(20)
  (a(18)-a(19) Karl Desfontaines / confirmed Jeff Sponaugle, a(20) Jeff
  Sponaugle, Oct 2026). A sextuplet can straddle 10^n only if
  10^n-3, 10^n+1, ... are all prime, and 10^n+1 is composite for n >= 3, so
  for n >= 4: a(n) = A063501(n) - A063501(n-1). This reproduces a(4)..a(17).
- Submit: **a(18) = 3439443854, a(19) = 24711327817, a(20) = 180637585474**,
  with the formula and a note on the straddle argument.

### A350825 — prime 5-tuples with initial member between 10^(n-1) and 10^n
- Now: 9 terms, last a(9) = 5838 (from b-files of A022006/A022007).
- A125517 (5-tuples with largest member <= 10^n, both patterns) has 16 terms
  (a(16) = 3473135825, Norman Luhn, Nov 2021). Same straddle argument (only a
  (p, p+4, p+6, p+10, p+12) tuple with p = 10^n - 3 could straddle, needing
  10^n+1 prime), so for n >= 4: a(n) = A125517(n) - A125517(n-1); this
  reproduces a(4)..a(9).
- Submit: **a(10) = 33193, a(11) = 204898, a(12) = 1306911, a(13) = 8665337,
  a(14) = 59193016, a(15) = 415670619, a(16) = 2988054630**. Depends on
  Luhn's A125517 values; we could verify a(13)..a(16) with a 5-tuple sieve
  (section 4) if an editor asks.

## 2. Our main sequence

### A007508 — twin prime pairs below 10^n
- Now: 19 terms, a(19) = 7237518093734545 (Chaffin, Jun 2026, corrected Jul
  2026). a(17)-a(18) Oliveira e Silva; a(16) Sebah (corrected Jacobsen 2014).
- We reproduced a(1)..a(15) with two independent implementations (CPU, GPU).
  No new term is within reach: a(20) is ~100 GB10-GPU-months.
- Optional submissions: a link to the programs once they are on GitHub
  (the entry links Nicely's, Oliveira e Silva's and Richstein's programs and
  tables), or a short comment that a(14)-a(15) were independently confirmed.
  Editors treat such comments as low priority; the link is more useful.

## 3. New terms with the existing twin sieve plus a small patch

The sieve counts pairs with k < K, p = 30k + r. Counting below an arbitrary
bound N (2^n, the 10^n-th prime) needs one addition: finish the last partial
block of 30 by testing its at most three candidates directly.

### A033843 — twin prime pairs with lesser member < 2^n
- Now: 41 terms, a(41) = 3879202049 (a(39)-a(41) Hugo Pfoertner, Feb 2024).
- Submit a(42)..a(57) with a b-file. GPU time, cumulative from 2^41:
  2^50 (1.1e15) ~1 h, 2^53 ~8 h, 2^55 ~1.5 days, 2^57 (1.4e17) ~7 days.
  Also recompute a(1)..a(41) as a check (2^41 = 2.2e12: seconds).

### A369109 — twin pairs with p <= 10^n and p ≡ 1 (mod 4)
- Now: 12 terms, a(12) = 935286453 (Amiram Eldar, Jun 2024).
- p mod 4 follows from the parity of k and the class r, so the count splits
  by masking alternate k in the popcount. Submit a(13)..a(17):
  a(15) ~45 min, a(16) ~7 h, a(17) ~5 days cumulative; a(18) ~2 months.

### A049035 / A093683 — twin pairs with smaller member <= the 10^n-th prime
- Now: 13 terms each, a(13) = 408231310520 (Eduard Roure Perdices, Jan 2024).
- Bounds are the 10^n-th primes (A006988), known far beyond this. Submit
  a(14) (bound 3.5e15, ~3 h), a(15) (3.7e16, ~1.5 days), a(16) (3.9e17,
  ~3 weeks). The two sequences differ only in the boundary convention.

### A146214 — the 10^n-th lesser twin prime; A147797 — A001097(10^n)
- Now: A146214 has 11 terms, a(10) = 6100479510551; A147797 has 14 terms,
  a(13) = 4646528385435887 (Donovan Johnson, 2009).
- Needs locating the exact m-th survivor: count chunks until the running
  total passes m, then enumerate the last segment (half a day of work on the
  CPU tool). Submit A146214 a(11)..a(14) (the 10^14-th pair is near 1.1e17,
  ~6 days) and A147797 a(14)..a(15) (near 5.4e16 and 6e17: 2 days, 4 weeks).

## 4. New terms needing a generalized constellation sieve

The twin sieve is a "pairs (p, p+d) in residue classes R mod 30" sieve with
R = {11,17,29}, d = 2. Making R and the offset list parameters gives cousin
pairs ({7,13,19}, +4), sexy pairs ({1,7,11,13,17,23}, +6), and constellations
(one or two classes, 3 to 6 progressions per prime). The GPU kernels only
need the wheel tables and presieve patterns regenerated; roughly a day or
two of work and validation against brute force.

| sequence | counts | now | candidate terms | GPU time |
|---|---|---|---|---|
| A396644 | cousin pairs < 2^n (upper member) | 39 terms, 2^39 | a(40)..a(57) | ~1 week to 2^57 |
| A370158 | {p,p+2,p+6,p+12}, p < 10^n | a(14) (Martin Ehrenstein, Feb 2024) | a(15)..a(18) | 1 h, 8 h, 4 days, 5 weeks |
| A125517 | 5-tuples, largest <= 10^n | a(16) (Norman Luhn, Nov 2021) | a(17), a(18) | ~4 days, ~5 weeks |
| A050258 | quadruplets {p,p+2,p+6,p+8}, largest < 10^n | a(17) (Jonathan Webster, Jun 2018) | a(18), maybe a(19) | ~4 weeks; a(19) ~1 year |
| A055737 | triples {p,p+2,p+6} or {p,p+4,p+6} < 10^n | a(14) (Charles R Greathouse IV, Feb 2022) | a(15)..a(18) | 1 h, 8 h, 4 days, 5 weeks |
| A118552 | sum of twin prime pairs < 10^n | a(15) (Chaffin, Jun 2026) | a(16), a(17) | 7 h, 5 days (needs a 128-bit sum accumulator) |
| A152127 | sum of cousin primes < 10^n | a(15) (Chaffin, Jun 2026) | a(16), a(17) | same, with the cousin classes |

Times assume twin-sieve cost per number. Sparser constellations are cheaper
per number because their candidate classes are rarer: a quadruplet sieve with
a mod-30030 wheel has one candidate per 155 integers, and whatever produced
the A063501 terms to 10^20 should handle A050258 to 10^19 in comparable time.

## 5. Already extended by others; nothing to add below 10^20

| sequence | counts | state |
|---|---|---|
| A080840, A152052 | cousin primes < 10^n (upper / lower member) | a(19), Chaffin Jun 2026 (corrected Jul 2026) |
| A080841 | pairs (p, p+6) with q < 10^n | a(19), Chaffin Jun 2026 |
| A063501 | sextuplets, largest <= 10^n | a(20), Jeff Sponaugle Oct 2026 |
| A007508 | twin pairs below 10^n | a(19), Chaffin 2026 |

## 6. Suggested order

1. Submit the derived terms in section 1 (A167874, A350826, A350825). Each
   needs only the formula, the straddle remark and the source of the inputs.
2. Patch the bound handling and run A033843 to 2^53 overnight, then to 2^57;
   run A369109 to a(16) alongside (both reuse the twin kernels unchanged).
3. Generalize the residue classes and offsets; do A396644 and A370158 first
   (cheap, recent sequences with active editors), then A050258(18).
4. Keep A007508 itself as a program link or comment only.

OEIS mechanics: terms beyond the ~260-character DATA line go in a b-file
(`b######.txt`, "n a(n)" per line); add an EXTENSIONS line "a(k)-a(m) from
Jeff Sponaugle, <date>"; for derived terms state the formula in the FORMULA
section and the dependency in a COMMENT; for computed terms editors like a
note on the method and an independent check (our CPU/GPU pair gives one).
