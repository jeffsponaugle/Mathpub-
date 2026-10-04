# OEIS submission package — A137723 (runs of numbers with a prime gap in their factorization)

Prepared 2026-10-03 by Jeff Sponaugle (with Claude). Blocks are in OEIS field style so they can be pasted into the
edit form; sign new lines "_Jeff Sponaugle_, Oct 03 2026". Supporting files: `b137723.txt`, `structural_terms.txt`,
`primes_to_certify.txt`, `verify_run.py`, `a137723.c`.

## 0. OEIS status check (2026-10-03, via oeis.org internal format, revision history and draft page)

**The main result is already in the OEIS.** The Sep 16 2026 submission (revisions #9–#11) was reviewed by Michel
Marcus and Hugo Pfoertner and approved by Andrew Howroyd on Sep 17 2026 (revision #17); the b-file with n = 1..193
was installed as revision #18. No pending draft exists ("No pending changes").

| field | current content (revision #18) | note |
|---|---|---|
| DATA | a(1)..a(40), ending 33554394, 20810, 74960 | correct (an HTML rendering showed a stray leading "3"; the %S line starts with 10) |
| b-file | n = 1..193, our file with 5 `#` header lines | identical to `b137723.txt` |
| %C 1 | A073490(a(n)+k)>0 for 0<=k<n and A073490(a(n)-1)=A073490(a(n)+n)=0. | Zumkeller's original, fine |
| %C 2 | "Continuation after the missing a(14): 1832, 1261, ... missing, 15684, ..." | **obsolete**: all of these are now in the data |
| %C 3 | "a(32) > 10^11. - Lucas A. Brown, Oct 07 2024" | **obsolete**: a(32) is in the data |
| %C 4 | "a(194) > 2^64-194. - Jeff Sponaugle, Sep 16 2026" | added by H. Pfoertner from our note; still the best bound |
| %H | b-file; Lucas A. Brown, Python program | no link to our program yet |
| %e | a(5) = 84 example | fine; nothing about the large terms |
| %Y | Cf. A073492 | thin |
| %K | nonn,hard | "more" was dropped by the editors; nothing to do |
| %E | ... a(32) onward from Jeff Sponaugle, Sep 16 2026 | fine |

Related entries A073490, A073491, A073492 were not re-fetched today (last read Sep 16 2026); nothing in this package
edits them.

---

## 1. Follow-up edit to A137723 — ready now, no computation needed

Purpose: remove the two obsolete comments, explain *why* the even-n terms are huge (this is what makes a 90-digit term
credible to a reader), add an example for a(32), link the program, widen the cross-references.

**Delete** these two comment lines (both superseded by the data):
```
%C Continuation after the missing a(14): 1832, 1261, 1130, 1332, 1638, missing, 1952,4298, 4524, missing, 5120, 16385, 2972, 4832, 5352, 10801, 5592, missing, 8468, missing, 9552, missing, 39462, missing, 20810, missing, 38502, missing, 15684, ...
%C a(32) > 10^11. - _Lucas A. Brown_, Oct 07 2024
```
(Keep Brown's program link. If an editor prefers to keep the bound as history, it can stay; it is simply no longer informative.)

**Add** comments:
```
%C a(n)-1 and a(n)+n are consecutive terms of A073491 (numbers without a prime gap in their factorization), so a(n)-1 is the first term of A073491 that is followed by a difference of exactly n+1. - _Jeff Sponaugle_, Oct 03 2026
%C For even n, n+1 is odd, so exactly one of a(n)-1 and a(n)+n is even, and an even term of A073491 is of the form 2^e_1*3^e_2*...*p_k^e_k with all e_i >= 1 (the first k primes). If also 3 | n+1, i.e. n == 2 (mod 6), then 3 divides the even bound together with n+1, so a prime odd bound is only possible when the even bound is a power of 2. This is why these terms are so large: for every n == 2 (mod 6) with n <= 188, a(n) is adjacent to a power of 2, a(n) = 2^k+1 or 2^k-n (a(32) = 2^54-32, a(38) = 2^25-38, a(44) = 2^38-44, a(50) = 2^81-50, a(56) = 2^28-56, a(62) = 2^164-62, a(68) = 2^49+1, ...). The terms with a(n) > 10^11 were obtained by enumerating all even terms of A073491 below a bound (2^64, 2^127 or 2^1000) and testing both neighbors at distance n+1, which is exhaustive. - _Jeff Sponaugle_, Oct 03 2026
%C For odd n both bounds are odd and in most cases prime, so a(n)-1 is usually the first prime p such that the gap to the next prime is n+1 and no term of A073491 lies in between; e.g. a(33) = 8468 rather than 1328, because the first prime gap of 34, 1327..1361, contains 1331 = 11^3 and 1350 = 2*3^3*5^2. - _Jeff Sponaugle_, Oct 03 2026
```
If the editors find three comments too many, the middle one carries the essential information.

**Add** an example:
```
%e a(32) = 18014398509481952 = 2^54 - 32: 2^54 - 33 is prime, 2^54 is a prime power, and each of the 32 numbers in between has a prime gap in its factorization, e.g. 2^54 - 32 = 2^5 * 127 * 4432676798593 (3..113 are skipped) and 2^54 - 1 = 3^4 * 7 * 19 * 73 * 87211 * 262657. - _Jeff Sponaugle_, Oct 03 2026
```

**Add** links (upload `a137723.c` as an a-file; the OEIS accepts program files alongside Brown's a137723.py):
```
%H Jeff Sponaugle, <a href="/A137723/a137723.c.txt">C program</a> (exhaustive scan with primesieve plus the structural search for even n)
```
(The OEIS usually stores uploaded programs with a .txt extension; use whatever name the upload form assigns. If the
repository goes on GitHub, link that instead and keep the file local.)

**Widen** cross-references:
```
%Y Cf. A073490, A073491, A073492, A000230 (first occurrence of prime gaps), A002386 (maximal prime gaps).
```

No change to %S/%T/%U, %K or %E in this edit.

---

## 2. New terms — needs runs (nothing has been started; run these when convenient)

The b-file stops at a(193) because a(194) is unknown. Everything below is exhaustive up to the stated limit, so each
result is final. Times are for the M1 Pro, 10 threads.

| step | command | time | what it settles |
|---|---|---|---|
| 1 | `./a137723 search 200 206 212 218 230 236 242 248 260 266 272 278 290 296 -l 2^1000` | ~10 s total | even n with 3 \| n+1 and block {3}: each is 2^a ± d plus powers of 3 |
| 2 | `./a137723 search 194 224 254 284 -l 2^1000` | a few min each | block {3, 5}; ~45 M candidates each |
| 3 | `./a137723 search 204 214 234 264 274 294 -l 2^1000` | several min each | block {5}; the k <= 2 partners need Miller-Rabin |
| 4 | `./a137723 search 244 -l 2^400` (then 2^500 if empty) | minutes | block {5, 7}: 4*10^9 candidates at 2^1000, 10^8 at 2^400 |
| 5 | `./a137723 next -N 10^13 -l 2^64 -o next_1e13.txt` | ~5 min | odd n: first unknown is a(387); the first prime gap of 388 lies beyond 10^11 |

Then, for each new term: `python3 verify_run.py L R` (L, R printed by `search`), append `n a(n) form limit` to
`structural_terms.txt`, run `python3 make_bfile.py next_1e13.txt` (or with `next_1e11.txt` if step 5 was skipped),
and submit:
```
%H Jeff Sponaugle, <a href="/A137723/b137723.txt">Table of n, a(n) for n = 1..K</a>      (K = new consecutive range)
%C a(K+1) > ...  (replace the a(194) comment by the new first unknown term and its bound; bounds print as "a(n) > ..." )
%E a(194)-a(K) from _Jeff Sponaugle_, Oct xx 2026
```
Expected outcome: all open even n <= 296 resolved (the structure predicts they exist below 2^1000 with high
probability; if one is missing, its line becomes a bound) and, if the 10^13 scan fills the odd n up to 463, a b-file
reaching about n = 300 or beyond. Terms that come from the scan are unconditional; terms from `search` rely on a
probable-prime test for one endpoint, see section 3.

---

## 3. Primality certification (recommended; closes the only non-rigorous step)

The even-n terms above 2^64 depend on one endpoint being prime, established by strong probable-prime tests (24 bases in
the 128-bit code, 30 rounds of `mpz_probab_prime_p` in GMP mode) and re-checked by `verify_run.py` (another 25 bases).
All 16 such endpoints are in `primes_to_certify.txt` (n and the prime, 21 to 112 digits). They are far below the sizes
where a proof costs more than seconds:
```bash
brew install pari
grep -v '^#' primes_to_certify.txt | awk '{print $2}' > /tmp/a137723_primes.txt
gp -q -f <<'GP'
v = readvec("/tmp/a137723_primes.txt");
for (i = 1, #v, print(#digits(v[i]), " digits: ", isprime(v[i])))
GP
```
(`isprime` is APR-CL, a proof, unlike `ispseudoprime`.) Nothing in the OEIS entry needs to change afterwards unless a
certification fails, which would be a genuine surprise; it would mean the "prime" endpoint is composite and the run
is longer than stated.

Terms below 2^64 (a(32), a(44), a(68), a(80), a(86), a(92), a(94), a(110), a(114), a(124), a(128), a(132), a(134),
a(142), a(144), a(158), a(160), a(172), a(174), a(184)) use the deterministic 7-base test and need nothing.

---

## 4. Optional: two new sequences derived from the data (searched 2026-10-03: no match in the OEIS)

Records of run length by starting point, the analog of A005250/A002386 (maximal prime gaps) for A073491. Both lists
are exhaustive below 10^11 (49 terms). Interest: from 360653 on every record coincides with a maximal prime gap
(A002386), but the early records differ because a composite term of A073491 can split a prime gap: the maximal gap of
86 after 155921 contains 155925 = 3^4*5^2*7*11, so the record run there is only 81 and the first run of 85 waits
until 338034.

**A. Record lengths of runs of consecutive numbers with a prime gap in their factorization**
```
%N Record values in A137723 ordered by the start of the run: lengths n of maximal runs of consecutive numbers having a prime gap in their factorization (A073492) that are longer than every run starting earlier.
%S 1,3,4,5,6,9,10,11,13,17,18,19,21,27,28,29,31,33,35,43,51,71,81,85,95,111,113,117,131,147,153,179,209,219,221,233,247,249,281,287,291,319,335,353,381,383,393,455,463
%C The run of length a(n) starts at A?????(n)+1 (sequence B below). For a(n) >= 95 all runs found so far start at a maximal prime gap: a(n) = A005250(k)-1 for the corresponding k. Terms verified exhaustively below 10^11. - _Jeff Sponaugle_
%Y Cf. A137723, A073491, A073492, A005250, A002386.
%K nonn,hard
%O 1,2
```
**B. Where the record runs start: the term of A073491 preceding each record run**
```
%N Terms of A073491 that are followed by a record gap to the next term of A073491: a(n) = A137723(A?????(n)) - 1 with A????? = sequence A above.
%S 9,19,54,83,113,181,199,467,773,1129,1331,1637,1951,2971,4831,5351,5591,8467,9551,15683,19609,31397,155925,338033,360653,370261,492113,1349533,1357201,2010733,4652353,17051707,20831323,47326693,122164747,189695659,191912783,387096133,436273009,1294268491,1453168141,2300942549,3842610773,4302407359,10726904659,20678048297,22367084959,25056082087,42652618343
%C From a(25) = 360653 on, every term so far is also a term of A002386 (prime starting a maximal prime gap); 155925 = 3^4*5^2*7*11 and 1331 = 11^3 are examples of composite terms. Verified exhaustively below 10^11. - _Jeff Sponaugle_
%Y Cf. A137723, A073491, A002386, A005250.
%K nonn,hard
%O 1,1
```
To extend both: `./a137723 scan 10^13 -o scan_1e13.txt` and recompute the running maximum of a(n) over the table
(the record run below 10^11 is 463 starting at 42652618344, the maximal prime gap of 464). Submit only if you think
the interplay with maximal prime gaps is worth an entry; the data is solid either way.

---

## 5. Checklist

- [ ] Section 1 edit entered (delete 2 comments, add 3 comments + example + program link + crossrefs)
- [ ] Section 3 certification run (`primes_to_certify.txt`, 16 primes) — optional but cheap
- [ ] Section 2 runs done; `verify_run.py` on each new term; `structural_terms.txt` + `make_bfile.py` updated
- [ ] New b-file uploaded, a(194) comment replaced, %E line added
- [ ] Decide on section 4
