# OEIS submission draft — prime sextuplet counts (A350826, A063501)

Draft of 2026-10-01. Everything below has been computed or checked with the tools in this
directory (`a350826.c`, `cuda/a350826_cuda.cu`, `verify_a350826.py`, `verify/reverify_oeis.py`);
see README.md for methods and run logs (`runs/`).

## Summary

| item | status |
|---|---|
| A350826(18) = 3439443854 | **new to OEIS**, ready (own GPU run; equals the value implied by Desfontaines' pi_6) |
| A350826(19) = 24711327817 | **new to OEIS**, ready (own GPU run; equals the value implied by Desfontaines' pi_6) |
| A063501(18) = 4010758480 | **new to OEIS**, ready (= pi_6(10^18), Desfontaines; confirmed) |
| A063501(19) = 28722086297 | **new to OEIS**, ready (= pi_6(10^19), Desfontaines; confirmed) |
| A350826(20) = 180637585474 | **new**, computed 2026-10-02 (two DGX Sparks, 34.8 h + 25.6 h GPU); not in Luhn's table either |
| A063501(20) = 209359671771 | **new** (= pi_6(10^20)) |
| pi_6(2^64) | our value **48629688861** differs from Luhn's table (Desfontaines) **48629687343** by 1518; ours reproduced twice (see section 3) |
| A350826 comment "as far as we know ... same length" | can be replaced by a verified statement for 3 <= n <= 24 |
| Dead link in A350826 and A063501 (`PI_06.html`) | should point to `PI_06.php` |
| A350826(1..17), A063501(1..17), A022008 b-file, A271000 b-file, A200503/A200504/A233426 (terms 1..56), A343636 (n = 0..23) | re-verified, all agree |

Recommendation: submit a(18)-a(20) for both sequences in one edit each. Contact Norman Luhn about the
pi_6(2^64) discrepancy (section 3) first or in parallel; it does not affect the OEIS terms.

---

## 1. A350826 — Number of prime sextuplets with n-digit initial term

Current entry: #24, Jan 21 2022, data a(1..17), keywords `nonn,base,more,hard`.

### DATA (append three terms)

```
1,1,0,0,3,0,13,64,235,1296,7013,41782,253420,1607418,10520883,70785653,488096844,3439443854,24711327817,180637585474
```

### COMMENTS — replace the second comment

Current:

> For n = 1 and n = 2 (see Example), the last member of the sextuplet has one digit more than
> the initial member (so the count would be 0 for these two, if all terms of the sextuplet had to
> have the same length). As far as we know, for all n > 2, all members of the sextuplets have the
> same length. A sufficient condition for this is that A033874(n) > 16.

Proposed:

```
For n = 1 and n = 2 (see Example), the last member of the sextuplet has one digit more than the initial member (so the count would be 0 for these two, if all terms of the sextuplet had to have the same length). For 3 <= n <= 24, all members of each sextuplet have the same length: no initial member p satisfies 10^n - 16 <= p < 10^n. A sufficient condition for this is A033874(n) > 16, which fails for n = 3, 5, 7, 8, 12, 15, 17, 18, 20; for those n it was checked directly. Hence a(n) = A063501(n) - A063501(n-1) for 4 <= n <= 24.
```

(Checked with `a350826 count 10^n-16 10^n` for n = 3..24, `verify/reverify_oeis.txt`.)

### Optional comment (rigour of the counting)

```
A base-2 strong probable-prime test of the six members is not sufficient to count prime sextuplets: there are 5 numbers p in (10^17, 10^18) and 3 in (10^18, 10^19) for which p, p+4, p+6, p+10, p+12, p+16 are all base-2 strong probable primes but one of them is composite (a base-2 strong pseudoprime), e.g., p = 156502671571374757 with p+4 = 212630861*736029901.
```

The eight cases (all re-checked with Python's own primality test and Pollard rho):

| p | composite member |
|---|---|
| 156502671571374757 | p+4 = 212630861 · 736029901 |
| 174519900169098217 | p = 79150597 · 2204909461 |
| 273314385115345027 | p+10 = 426860153 · 640290229 |
| 483757682409253657 | p+6 = 169691191 · 2850811993 |
| 532522504106799457 | p+16 = 161832109 · 3290586197 |
| 1157159728167006787 | p+6 = 563142103 · 2054827231 |
| 4211906013907595737 | p+4 = 1510446653 · 2788516897 |
| 6526567117896650257 | p+16 = 154618367 · 42210813919 |

### Optional FORMULA (asymptotic, conjectural)

```
Conjecture (Hardy-Littlewood): a(n) ~ C_6 * Integral_{10^(n-1)..10^n} dt/log(t)^6, C_6 = 17.29861... . For n = 17..20 the integral differs from a(n) by less than 5*10^-5 a(n) (2*10^-6 for n = 20); it gives a(21) ~ 1.34117*10^12.
```

(C_6 = prod over primes q of (1 - w(q)/q)/(1 - 1/q)^6 with w(q) the number of residues of
{0,4,6,10,12,16} mod q; computed over q < 3*10^6 plus the tail correction.)

### LINKS

Fix the dead link (`.html` now returns "300 Multiple Choices"):

```
Norman Luhn, <a href="https://www.pzktupel.de/counting/PI_06.php">PI_6(10^n)</a>
```

Optional, only if the programs are published somewhere public:

```
Jeff Sponaugle, <a href="...">C and CUDA programs for counting prime sextuplets</a>
```

### CROSSREFS (add)

```
Cf. A350825 (5-tuples), A350827 (septuplets), A350828 (octuplets).
```

### EXTENSIONS (add)

Suggested wording (credit to Desfontaines, who computed pi_6(10^18) and pi_6(10^19) first):

```
a(18)-a(19) from Karl Desfontaines's counts in Norman Luhn's table (May 2026), independently confirmed by _Jeff Sponaugle_, Oct 01 2026
a(20) from _Jeff Sponaugle_, Oct 02 2026
```

(If the editors prefer, simply `a(18)-a(19) from _Jeff Sponaugle_, Oct 01 2026` with the
Desfontaines credit in the notes. Desfontaines has an OEIS b-file in A200503, so `_Karl Desfontaines_`
may work as a user link.)

Keyword `more` stays.

---

## 2. A063501 — Number of prime-sextuplets up to 10^n

Current entry: #25, Feb 10 2022, data a(1..17), offset 1,3, keywords `nonn,more`.

### DATA (append three terms)

```
0,1,2,2,5,5,18,82,317,1613,8626,50408,303828,1911246,12432129,83217782,571314626,4010758480,28722086297,209359671771
```

### COMMENTS (add)

```
A sextuplet (p, p+4, p+6, p+10, p+12, p+16) is counted when its largest member p+16 is at most 10^n; thus a(1) = 0 although p = 7 < 10.
```

### FORMULA (add)

```
a(n) = Sum_{k=1..n} A350826(k) for 3 <= n <= 24 (no prime sextuplet has members on both sides of 10^n for these n).
```

### LINKS

Same fix as in A350826:

```
Norman Luhn, <a href="https://www.pzktupel.de/counting/PI_06.php">PI_6(10^n)</a>
```

### CROSSREFS (add)

```
Cf. A343636, A350826.
```

### KEYWORDS

Consider adding `hard` (as in A350826).

### EXTENSIONS (add)

```
a(18)-a(19) from Karl Desfontaines (May 2026, see Luhn's table), independently confirmed by _Jeff Sponaugle_, Oct 01 2026
a(20) from _Jeff Sponaugle_, Oct 02 2026
```

---

## 3. a(20) — result and verification

`a350826_cuda count 1e19 1e20 -b 2^64,2e19,3e19,2^65,4e19,5e19,6e19,7e19,2^66,8e19,9e19` (wheel 37,
B = 2^16, 92569 chunks of 16384 classes) on two DGX Sparks: atom1 chunks 0..53417 (2026-10-01 00:55 to
2026-10-02 11:12) and 91858..92568, atom2 chunks 53416..91861 (from 2026-10-01 10:07); merged with
`combine_a20.py` (`runs/a20_combined.txt`): every chunk exactly once, the 10 chunks computed by both GPUs
have identical lines, candidates − count = 15 pseudo-sextuplets.

```
RESULT [1e19, 1e20) count 180637585474 cks 6fde7e53b8dbdda6  A350826(20)
RESULT A063501(20) = pi_6(10^20) = 209359671771
```

pi_6 at the bin boundaries (cumulative from pi_6(10^19) = 28722086297):

| x | sextuplets in [previous, x) | pi_6(x) |
|---|---|---|
| 2^64 | 19907602564 | 48629688861 |
| 2·10^19 | 3506083177 | 52135772038 |
| 3·10^19 | 21809123277 | 73944895315 |
| 2^65 | 14448377987 | 88393273302 |
| 4·10^19 | 6389745264 | 94783018566 |
| 5·10^19 | 20148907847 | 114931926413 |
| 6·10^19 | 19619158182 | 134551084595 |
| 7·10^19 | 19190757809 | 153741842404 |
| 2^66 | 7170897125 | 160912739529 |
| 8·10^19 | 11662053851 | 172574793380 |
| 9·10^19 | 18526297674 | 191101091054 |
| 10^20 | 18258580717 | 209359671771 |

**Discrepancy at 2^64.** Luhn's table lists pi_6(2^64) = 48629687343 (Desfontaines, Jul 2026); we get
48629688861, i.e. 1518 more sextuplets in [10^19, 2^64). Checks on our side:

* a full recount of [10^19, 2^64) on both GPUs with the halves swapped between the machines reproduced
  every chunk exactly (92569 of 92569 chunk lines identical, total 19907602564);
* the CPU program (independent marking code and primality tests) recomputed 40 random chunks of
  [10^19, 2^64): no mismatch;
* every counted sextuplet passed BPSW in all six members (no BPSW pseudoprime exists below 2^64), so our
  count cannot contain composites; the sieve only removes candidates; the same programs reproduce
  pi_6(10^18) and pi_6(10^19) of Desfontaines exactly;
* the 1518th sextuplet below 2^64 lies 6.82e11 below it (no natural cut-off there: [2^64 − 2^39, 2^64)
  holds 1211 sextuplets, [2^64 − 2^40, 2^64) holds 2453), so the difference is not simply a missing top
  segment of power-of-two size.

Our value therefore looks right; Luhn / Desfontaines should be told (it does not affect any OEIS term).
Pseudo-sextuplets excluded in [10^19, 10^20): 15 (12 listed in the README run notes, all verified).

## 4. Re-verified data (for the notes to the editors; no edits needed)

| sequence | what was checked | result |
|---|---|---|
| A350826 | all 17 terms a(1..17) (GPU tool; the CPU tool and the Python program for a(1..14) and a(1..10)), identical checksums wherever two programs ran | agree |
| A063501 | a(1..17) via the partial sums of A350826 and the straddling check | agree |
| Luhn PI_06 table | pi_6(10^18) = 4010758480, pi_6(10^19) = 28722086297 (Desfontaines, May 2026) | confirmed |
| A022008 | b-file n = 1..10000 (D. A. Corneth / Z. Seidov) against `a350826 list` | agree |
| A271000 | b-file n = 1..5940 (A. Wesolowski) | agree |
| A200503 | b-file terms 1..56 (all record gaps ending below 10^15, from the sorted list of all 12432129 sextuplets below 10^15) | agree |
| A200504 | b-file n = 1..56 (A. Kourbatov), same computation | agree |
| A233426 | b-file n = 1..56 (A. Kourbatov), same computation | agree |
| A343636 | a(0..23) of the b-file (M. F. Hasler / N. Luhn) | agree |

Largest record gap below 10^15: 1457965740, from 437804272277497 to 437805730243237
(= A200503(56), A200504(56), A233426(56)). Details: `verify/reverify_oeis.txt`.

---

## 5. Notes to the editors (paste into the edit form)

```
a(18)-a(20) were computed with C and CUDA programs: candidates p are restricted to the 1516640125 residue classes modulo 2*3*5*...*37 that are admissible for the pattern (0,4,6,10,12,16); each class is sieved as an arithmetic progression with all primes up to 2^16; every member of every survivor gets a base-2 strong probable-prime test, and every member of every candidate that passes all six gets a strong Lucas test (BPSW, which has no counterexample below 2^64). The values agree with PI_6(10^18) = 4010758480 and PI_6(10^19) = 28722086297 of Karl Desfontaines (May 2026) in Norman Luhn's table. Runs: a(18) 30.5 min, a(19) 4.0 h on one NVIDIA GB10 (DGX Spark), a(20) 60 GPU-hours on two. Checks: the programs reproduce a(1)..a(17); an independent Python program (numpy sieve, Python's own Miller-Rabin with the first 13 prime bases) agrees in count and checksum on 2*10^12-wide windows near 10^15, 10^17, 10^18, 10^19, 2^64, 10^20 and 10^22; the CPU and GPU programs give identical per-chunk counts and checksums, and randomly chosen chunks of both GPU runs were recomputed on the CPU. The Lucas test matters: 5 (in a(18)), 3 (in a(19)) and 15 (in a(20)) candidates pass the base-2 test in all six members but contain a base-2 strong pseudoprime. The link to Luhn's table is updated (the .html address is dead). a(20) = 180637585474 (new) was computed the same way on two DGX Sparks (60 GPU-hours) and its [10^19, 2^64) part was recounted with the machines swapped (identical per-chunk counts and checksums).
```

---

## 6. Possible further submissions (not ready yet)

* **New sequence: number of prime sextuplets below 2^n** (no such entry exists; searched
  "number of prime sextuplets", "sextuplets 2^n"). Needs pi_6(2^n) for all n up to 66: n <= 60 takes
  minutes, 2^61..2^63 about 4 GPU-hours, 2^64..2^66 come from the a(20) run.
* **New sequence: "pseudo-sextuplets"**, numbers p such that p, p+4, p+6, p+10, p+12, p+16 are all
  base-2 strong probable primes but not all prime. Eight terms in (10^17, 10^19) are known (table in
  section 1); [10^15, 10^16) has none; [0, 10^15) and (10^16, 10^17) still have to be listed
  (minutes on the GPU) before the sequence could start at its first term.
* **A200503 / A200504 / A233426**: nothing new (Desfontaines already extended A200503 to 91 terms).
* **Luhn's table**: pi_6(10^20) = 209359671771 and pi_6(2^65) = 88393273302, pi_6(2^66) = 160912739529,
  pi_6(k·10^19) (section 3) are new; pi_6(2^64) = 48629688861 would correct the table.

---

## 7. Before submitting

* Re-check https://oeis.org/draft/A350826 and https://oeis.org/draft/A063501 for pending edits
  (none on 2026-10-01).
* Re-check Luhn's table (https://www.pzktupel.de/counting/PI_06.php) for new pi_6 values, especially
  pi_6(10^20).
* Decide on the credit wording in EXTENSIONS (section 1).
* Confirm your OEIS user name for the `_Jeff Sponaugle_` links.
