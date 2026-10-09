# OEIS submission candidates

Generated Oct 06 2026 by `tools/submissions.py` from the dbpal results (`results/`), the search log
(`results/search_log.jsonl`, which certifies each bound length by length) and the current OEIS entries.
Attribution below uses _Jeff Sponaugle_; change it if your OEIS user name differs. Lines are in OEIS internal
format (`%S` data, `%C` comment, `%H` link, `%E` extension) so they can be pasted into the edit form.

General notes for every submission:

* **Completeness.** Every bound below means *all* numbers under it were searched: each length
  (in the enumerated base) is logged as done only after the whole length finished. Lengths split
  across two machines count only when all parts used an identical, fingerprinted search plan.
* **Verification.** Every reported number was re-checked independently (GMP and a separate
  pure-Python checker: palindromic in both bases, Baillie-PSW for primality). For each sequence,
  every term already in the OEIS below the bound was re-found (none missing). The engines are
  cross-checked against brute force and against 14 full OEIS b-files (`make test`).
* **Method** (for a comment or a program link, if wanted): exhaustive meet-in-the-middle search.
  Fix the outer digits of a palindrome in one base; that fixes the leading digits in the other
  base and hence (palindrome) its trailing digits; matching middles are found by table lookup,
  ~N^(1/4) work for coprime bases. GPU (Apple Metal and NVIDIA CUDA) + CPU implementation, 'dbpal'.
* **b-files** are in `results/bfiles/` (format `n a(n)`, with `#` header lines that you may delete).

How to submit (one edit per sequence; editors review each):

1. Sign in at oeis.org, open the sequence, choose **Edit**.
2. For a b-file: upload it first (the b-file link on the edit page), then change the `%H` line.
3. Paste the `%S` / `%C` / `%E` lines below; put a short note in the edit summary, e.g.
   "exhaustive search to X with a GPU meet-in-the-middle program; new term(s) verified, primality
   proven (Pocklington)".
4. Start with section 1 (new terms), then the bound comments, then the b-files. Several
   sequences of one family can be submitted the same day, but keep each edit self-contained.
5. If you publish the program (e.g. on GitHub), a `%H` link to it can be added to A393014 and the
   entries whose terms it found, e.g. `%H A393014 Jeff Sponaugle, <a href="https://...">dbpal</a>, GPU
   search for numbers palindromic in two bases.`

Bounds are those of the search log when this file was generated. Re-run
`python3 tools/submissions.py --refresh` after more searching or after OEIS edits, so that this
file proposes only what is still missing from the OEIS.

## 0. Edits already made on the OEIS

From the OEIS revision histories (fetched Oct 06 2026 14:50).
What is published or pending review is not proposed again below.

| Sequence | Your edit | OEIS status | Still to do |
|---|---|---|---|
| A046478 | a(9); comment `a(10) > 10^49, if it exists.` (Oct 06 2026) | proposed, awaiting review | nothing |
| A046479 | a(10)-a(11); comment `a(12) > 11^45 (about 7.29*10^46), if it exists.` (Oct 03 2026) | approved, published Oct 03 2026 (#26) | new bound comment `a(12) > 11^47 (about 8.82*10^48), if it exists.` replacing yours (section 2) |
| A046480 | a(9); comment `a(10) > 12^49 (about 7.58*10^52), if it exists.` (Oct 03 2026) | approved, published Oct 03 2026 (#20) | nothing |

## 1. New terms

### A060792 — Numbers that are palindromic in bases 2 and 3.

Currently 17 terms (largest 24014998963383302600955162866787153652444049); keywords `nonn,base,hard,nice`.
Existing bound comment: "a(18) (if it exists) is greater than 3^93. - _Ilya Nikulshin_, Feb 22 2016"

Proposed edits:

```
   DATA is already 252 characters: the new term(s) go in the b-file only
%E A060792 a(18)-a(21) from _Jeff Sponaugle_, Oct 06 2026
%C A060792 a(22) > 3^117 (about 6.66*10^55), if it exists. - _Jeff Sponaugle_, Oct 06 2026
%H A060792 Jeff Sponaugle, <a href="/A060792/b060792.txt">Table of n, a(n) for n = 1..21</a> (terms 1..17 from Ilya Nikulshin; terms a(11)..a(15) from Alan Grimes and _Keith F. Lynch_; a(16) from _Japheth Lim_; all terms < 3^117)
   (upload results/bfiles/b060792.txt)
```

* 93814833782752683486286194707262031368542962943489
  * base 2: `10000000011000011010010101000111110110101100100111010010001101100101110010111001110001110011101001110100110110001001011100100110101101111100010101001011000011000000001`
  * base 3: `202020002112102111210120212201001012122211221122012212210221122112221210100102212021012111201211200020202`
* 42468857786616486687891165709191634653913393785544903
  * base 2: `1110001100000100101111100100110010101100010001011111110000011001111010111000101010101111111010101010001110101111001100000111111101000100011010100110010011111010010000011000111`
  * base 3: `110120002221012001102221011212122120011112200020100222212222001020002211110021221212110122201100210122200021011`
* 676710174766987245264188686637424310057228994099521607
  * base 2: `11100010000101011111110010000010110001101001110011011111101110010111100101011110001000011111000010001111010100111101001110111111011001110010110001101000001001111111010100001000111`
  * base 3: `21102010110002110011211122200100001110010000000002002211111220020000000001001110000100222111211001120001101020112`
* 30622032620764000016137033607275848398941862554899405817
  * base 2: `10011111110110101011111001111011010111111011001001011101111111100110100100110111011001001110101110010011011101100100101100111111110111010010011011111101011011110011111010101101111111001`
  * base 3: `110102102001011210022110212010122001120222001021200012110212011210002120100222021100221010212011220012110100201201011`

The existing bound comment(s) are superseded by the new term; editors may delete or keep them.

Notable: first new term of A060792 since 2014 (a(18) > 3^93 was the last published progress, 2016).

## 2. New search bounds for prime sequences (no new term)

Each line is a proposed `%C` comment; where an older bound comment exists it can be replaced.
"(yours)" marks your own earlier comment: edit it in place rather than adding a second one.

| Sequence | Bases | Terms | Existing bound | Proposed comment |
|---|---|---|---|---|
| A046472 | 2, 10 | 6 | — | `a(7) > 10^55, if it exists. - _Jeff Sponaugle_, Oct 06 2026` |
| A046473 | 3, 10 | 7 | a(8) > 10^32 if it exists. | `a(8) > 10^49, if it exists. - _Jeff Sponaugle_, Oct 06 2026` |
| A046474 | 4, 10 | 8 | — | `a(9) > 10^54, if it exists. - _Jeff Sponaugle_, Oct 06 2026` |
| A046475 | 6, 10 | 9 | — | `a(10) > 10^55, if it exists. - _Jeff Sponaugle_, Oct 06 2026` |
| A046476 | 7, 10 | 9 | — | `a(10) > 10^47, if it exists. - _Jeff Sponaugle_, Oct 06 2026` |
| A046477 | 8, 10 | 10 | a(11) > 10^34 if it exists. | `a(11) > 10^55, if it exists. - _Jeff Sponaugle_, Oct 06 2026` |
| A046479 | 10, 11 | 11 | a(10) > 7.86174 * 10^18 using b-file at A029966.; a(10) > 10^36 if it exists.; a(12) > 11^45 (about 7.29*10^46), if it exists. (yours) | `a(12) > 11^47 (about 8.82*10^48), if it exists. - _Jeff Sponaugle_, Oct 06 2026` |
| A046481 | 10, 13 | 10 | — | `a(11) > 13^41 (about 4.70*10^45), if it exists. - _Jeff Sponaugle_, Oct 06 2026` |
| A046482 | 10, 14 | 8 | — | `a(9) > 14^46 (about 5.27*10^52), if it exists. - _Jeff Sponaugle_, Oct 06 2026` |
| A046483 | 10, 15 | 7 | — | `a(8) > 15^49 (about 4.25*10^57), if it exists. - _Jeff Sponaugle_, Oct 06 2026` |
| A046484 | 10, 16 | 13 | — | `a(14) > 10^51, if it exists. - _Jeff Sponaugle_, Oct 06 2026` |
| A393014 | 2, 9 | 9 | — | `a(10) > 9^54 (about 3.38*10^51), if it exists. - _Jeff Sponaugle_, Oct 06 2026` |

Already in the OEIS (published or pending), nothing to add: A046478 (`a(10) > 10^49, if it exists.`), A046480 (`a(10) > 12^49 (about 7.58*10^52), if it exists.`).

Notes:

* A393014: the largest known term is a(9) = 2073722573406568398837217962245683 (about 2.07*10^33).
  The search covered all numbers below 9^54 (every base-9 length up to 53; even base-9 lengths cannot
  occur, since an even-length base-9 palindrome is a multiple of 10). Worth adding to A393014 also:
  `%C A393014 a(10) > 9^54 (about 3.38*10^51), if it exists. - _Jeff Sponaugle_, Oct 06 2026`.
* A046472: the bound above is from this search. In addition, none of the 183 known terms of
  A007632 (its b-file, complete to 53 decimal digits according to E. Schacham, plus eight further
  terms with 54-55 digits) is prime, which would give a(7) > 10^53 if that b-file is complete.

## 3. Extended b-files (numbers palindromic in two bases)

Each b-file is complete up to the stated bound and starts with the existing terms (checked equal).
Upload the file from `results/bfiles/` and replace the `%H` b-file line with the one given.

| Sequence | Bases | Existing b-file | New b-file | New terms |
|---|---|---|---|---|
| A007632 | 2, 10 | 1..175 | 1..183 (all terms < 10^55) | 8 |
| A007633 | 3, 10 | 1..99 | 1..150 (all terms < 10^49) | 51 |
| A029731 | 10, 16 | 1..81 | 1..202 (all terms < 10^51) | 121 |
| A029804 | 8, 10 | 1..75 | 1..192 (all terms < 10^55) | 117 |
| A029961 | 4, 10 | 1..56 | 1..160 (all terms < 10^54) | 104 |
| A029962 | 5, 10 | 1..79 | 1..200 (all terms < 10^54) | 121 |
| A029963 | 6, 10 | 1..102 | 1..237 (all terms < 10^55) | 135 |
| A029964 | 7, 10 | 1..65 | 1..155 (all terms < 10^47) | 90 |
| A029965 | 9, 10 | 1..66 | 1..163 (all terms < 10^49) | 97 |
| A029966 | 10, 11 | 1..95 | 1..177 (all terms < 11^47) | 82 |
| A029967 | 10, 12 | 1..57 | 1..190 (all terms < 12^49) | 133 |
| A029968 | 10, 13 | 1..76 | 1..157 (all terms < 13^41) | 81 |
| A029969 | 10, 14 | 1..70 | 1..186 (all terms < 14^46) | 116 |
| A029970 | 10, 15 | 1..72 | 1..197 (all terms < 15^49) | 125 |
| A099165 | 10, 32 | 1..115 | 1..218 (all terms < 10^40) | 103 |
| A182232 | 2, 5 | 1..25 | 1..62 (all terms < 5^58) | 37 |
| A182233 | 2, 6 | 1..70 | 1..174 (all terms < 6^65) | 104 |
| A182234 | 2, 7 | 1..31 | 1..72 (all terms < 7^48) | 41 |
| A249155 | 6, 15 | 1..71 | 1..162 (all terms < 15^43) | 91 |
| A249156 | 5, 7 | 1..72 | 1..152 (all terms < 7^48) | 80 |
| A249157 | 11, 13 | 1..69 | 1..154 (all terms < 13^36) | 85 |
| A249158 | 7, 29 | 1..80 | 1..166 (all terms < 29^28) | 86 |
| A250408 | 10, 20 | 1..85 | 1..222 (all terms < 20^31) | 137 |
| A250409 | 10, 24 | 1..99 | 1..173 (all terms < 24^29) | 74 |
| A250411 | 10, 27 | 1..84 | 1..167 (all terms < 27^28) | 83 |
| A250412 | 10, 36 | 1..93 | 1..179 (all terms < 36^26) | 86 |
| A259374 | 3, 5 | 1..39 | 1..85 (all terms < 5^58) | 46 |
| A259375 | 3, 6 | 1..64 | 1..154 (all terms < 6^70) | 90 |
| A259376 | 4, 6 | 1..64 | 1..136 (all terms < 6^65) | 72 |
| A259377 | 3, 7 | 1..72 | 1..149 (all terms < 7^48) | 77 |
| A259378 | 4, 7 | 1..48 | 1..87 (all terms < 7^48) | 39 |
| A259381 | 3, 8 | 1..43 | 1..69 (all terms < 3^84) | 26 |
| A259383 | 5, 8 | 1..67 | 1..120 (all terms < 5^58) | 53 |
| A259384 | 6, 8 | 1..74 | 1..163 (all terms < 6^65) | 89 |
| A259385 | 2, 9 | 1..69 | 1..100 (all terms < 9^54) | 31 |
| A259387 | 4, 9 | 1..57 | 1..121 (all terms < 9^42) | 64 |
| A259388 | 5, 9 | 1..40 | 1..76 (all terms < 9^42) | 36 |
| A259389 | 6, 9 | 1..57 | 1..173 (all terms < 9^57) | 116 |
| A259390 | 7, 9 | 1..86 | 1..149 (all terms < 9^42) | 63 |

Proposed `%H` lines:

```
%H A007632 Jeff Sponaugle, <a href="/A007632/b007632.txt">Table of n, a(n) for n = 1..183</a> (terms 1..175 from Eshed Shaham; Terms 1..147 variously by Robert G. Wilson v, Charlton Harrison, Ilya Nikulshin, Andrey Astrelin; all terms < 10^55)
%H A007633 Jeff Sponaugle, <a href="/A007633/b007633.txt">Table of n, a(n) for n = 1..150</a> (terms 1..99 from Richard Jones; terms 1..41 from Robert Israel, terms 42..63 from Robert G. Wilson v, terms 64..65 from Patrick De Geest; all terms < 10^49)
%H A029731 Jeff Sponaugle, <a href="/A029731/b029731.txt">Table of n, a(n) for n = 1..202</a> (terms 1..81 from Robert G. Wilson v; all terms < 10^51)
%H A029804 Jeff Sponaugle, <a href="/A029804/b029804.txt">Table of n, a(n) for n = 1..192</a> (terms 1..75 from Robert G. Wilson v; all terms < 10^55)
%H A029961 Jeff Sponaugle, <a href="/A029961/b029961.txt">Table of n, a(n) for n = 1..160</a> (terms 1..56 from Robert G. Wilson v; all terms < 10^54)
%H A029962 Jeff Sponaugle, <a href="/A029962/b029962.txt">Table of n, a(n) for n = 1..200</a> (terms 1..79 from Robert G. Wilson v; all terms < 10^54)
%H A029963 Jeff Sponaugle, <a href="/A029963/b029963.txt">Table of n, a(n) for n = 1..237</a> (terms 1..102 from Robert G. Wilson v; all terms < 10^55)
%H A029964 Jeff Sponaugle, <a href="/A029964/b029964.txt">Table of n, a(n) for n = 1..155</a> (terms 1..65 from Robert G. Wilson v; all terms < 10^47)
%H A029965 Jeff Sponaugle, <a href="/A029965/b029965.txt">Table of n, a(n) for n = 1..163</a> (terms 1..66 from Robert G. Wilson v and Ray Chandler; terms < 10^18, first 52 terms from Robert G. Wilson v; all terms < 10^49)
%H A029966 Jeff Sponaugle, <a href="/A029966/b029966.txt">Table of n, a(n) for n = 1..177</a> (terms 1..95 from Max Alekseyev; terms for n = 1..79 from Ray Chandler and Robert G. Wilson v; all terms < 11^47)
%H A029967 Jeff Sponaugle, <a href="/A029967/b029967.txt">Table of n, a(n) for n = 1..190</a> (terms 1..57 from Robert G. Wilson v; all terms < 12^49)
%H A029968 Jeff Sponaugle, <a href="/A029968/b029968.txt">Table of n, a(n) for n = 1..157</a> (terms 1..76 from Robert G. Wilson v; all terms < 13^41)
%H A029969 Jeff Sponaugle, <a href="/A029969/b029969.txt">Table of n, a(n) for n = 1..186</a> (terms 1..70 from Robert G. Wilson v; first 63 terms from Ray Chandler; all terms < 14^46)
%H A029970 Jeff Sponaugle, <a href="/A029970/b029970.txt">Table of n, a(n) for n = 1..197</a> (terms 1..72 from Robert G. Wilson v; first 69 terms from Ray Chandler; all terms < 15^49)
%H A099165 Jeff Sponaugle, <a href="/A099165/b099165.txt">Table of n, a(n) for n = 1..218</a> (terms 1..115 from Ray Chandler and Robert G. Wilson v; terms a(88)-a(111) from Ray Chandler; all terms < 10^40)
%H A182232 Jeff Sponaugle, <a href="/A182232/b182232.txt">Table of n, a(n) for n = 1..62</a> (terms 1..25 from Ray Chandler; all terms < 5^58)
%H A182233 Jeff Sponaugle, <a href="/A182233/b182233.txt">Table of n, a(n) for n = 1..174</a> (terms 1..70 (those < 10^18) from Ray Chandler; all terms < 6^65)
%H A182234 Jeff Sponaugle, <a href="/A182234/b182234.txt">Table of n, a(n) for n = 1..72</a> (terms 1..31 (those < 10^18) from Ray Chandler; all terms < 7^48)
%H A249155 Jeff Sponaugle, <a href="/A249155/b249155.txt">Table of n, a(n) for n = 1..162</a> (terms 1..71 from Ray Chandler and Chai Wah Wu; terms < 6^28). First 65 terms from Ray Chandler; all terms < 15^43)
%H A249156 Jeff Sponaugle, <a href="/A249156/b249156.txt">Table of n, a(n) for n = 1..152</a> (terms 1..72 from Giovanni Resta; first 60 terms from Ray Chandler; all terms < 7^48)
%H A249157 Jeff Sponaugle, <a href="/A249157/b249157.txt">Table of n, a(n) for n = 1..154</a> (terms 1..69 (those < 10^18) from Ray Chandler; all terms < 13^36)
%H A249158 Jeff Sponaugle, <a href="/A249158/b249158.txt">Table of n, a(n) for n = 1..166</a> (terms 1..80 (those < 10^18) from Ray Chandler; all terms < 29^28)
%H A250408 Jeff Sponaugle, <a href="/A250408/b250408.txt">Table of n, a(n) for n = 1..222</a> (terms 1..85 from Robert G. Wilson v; all terms < 20^31)
%H A250409 Jeff Sponaugle, <a href="/A250409/b250409.txt">Table of n, a(n) for n = 1..173</a> (terms 1..99 from Robert G. Wilson v; all terms < 24^29)
%H A250411 Jeff Sponaugle, <a href="/A250411/b250411.txt">Table of n, a(n) for n = 1..167</a> (terms 1..84 from Robert G. Wilson v; all terms < 27^28)
%H A250412 Jeff Sponaugle, <a href="/A250412/b250412.txt">Table of n, a(n) for n = 1..179</a> (terms 1..93 from Robert G. Wilson v; all terms < 36^26)
%H A259374 Jeff Sponaugle, <a href="/A259374/b259374.txt">Table of n, a(n) for n = 1..85</a> (terms 1..39 from Giovanni Resta; all terms < 5^58)
%H A259375 Jeff Sponaugle, <a href="/A259375/b259375.txt">Table of n, a(n) for n = 1..154</a> (terms 1..64 from Giovanni Resta; all terms < 6^70)
%H A259376 Jeff Sponaugle, <a href="/A259376/b259376.txt">Table of n, a(n) for n = 1..136</a> (terms 1..64 from Giovanni Resta; all terms < 6^65)
%H A259377 Jeff Sponaugle, <a href="/A259377/b259377.txt">Table of n, a(n) for n = 1..149</a> (terms 1..72 from Giovanni Resta; all terms < 7^48)
%H A259378 Jeff Sponaugle, <a href="/A259378/b259378.txt">Table of n, a(n) for n = 1..87</a> (terms 1..48 from Giovanni Resta; all terms < 7^48)
%H A259381 Jeff Sponaugle, <a href="/A259381/b259381.txt">Table of n, a(n) for n = 1..69</a> (terms 1..43 from Giovanni Resta; all terms < 3^84)
%H A259383 Jeff Sponaugle, <a href="/A259383/b259383.txt">Table of n, a(n) for n = 1..120</a> (terms 1..67 from Giovanni Resta; all terms < 5^58)
%H A259384 Jeff Sponaugle, <a href="/A259384/b259384.txt">Table of n, a(n) for n = 1..163</a> (terms 1..74 from Giovanni Resta; all terms < 6^65)
%H A259385 Jeff Sponaugle, <a href="/A259385/b259385.txt">Table of n, a(n) for n = 1..100</a> (terms 1..69 from Max Alekseyev; terms for n = 1..44 from Giovanni Resta; all terms < 9^54)
%H A259387 Jeff Sponaugle, <a href="/A259387/b259387.txt">Table of n, a(n) for n = 1..121</a> (terms 1..57 from Giovanni Resta; all terms < 9^42)
%H A259388 Jeff Sponaugle, <a href="/A259388/b259388.txt">Table of n, a(n) for n = 1..76</a> (terms 1..40 from Giovanni Resta; all terms < 9^42)
%H A259389 Jeff Sponaugle, <a href="/A259389/b259389.txt">Table of n, a(n) for n = 1..173</a> (terms 1..57 from Giovanni Resta; all terms < 9^57)
%H A259390 Jeff Sponaugle, <a href="/A259390/b259390.txt">Table of n, a(n) for n = 1..149</a> (terms 1..86 from Giovanni Resta; all terms < 9^42)
```

Where the data section is short of the limit, a few more terms could also be appended to `%S`;
usually the b-file is enough.

## 4. Possible new sequences (optional)

Primes that are palindromic in two bases where the OEIS has no such sequence. An OEIS search for
the first terms of each list returned no match (Oct 2026). The (10, b) ones would continue the
family A046472-A046484 (which stops at b = 16); the others are the primes in an existing
"palindromic in bases b1 and b2" sequence (cross-reference it). Editors may consider some of these
too sparse or too arbitrary; submit only the ones you find worthwhile.

* **Primes that are palindromic in bases 10 and 20** (primes in A250408); complete below 20^31 (about 2.15*10^40): 2, 3, 5, 7, 11, 1917191, 145813396492606294693318541, 12345009773704925152940737790054321
* **Primes that are palindromic in bases 10 and 24** (primes in A250409); complete below 24^29 (about 1.06*10^40): 2, 3, 5, 7, 11, 1851581, 3673763, 7865687, 726674454476627, 37872782328727873, 95528099096969099082559
* **Primes that are palindromic in bases 10 and 27** (primes in A250411); complete below 27^28 (about 1.20*10^40): 2, 3, 5, 7, 11, 757, 919, 10301, 1975783979793875791, 92723610384148301632729
* **Primes that are palindromic in bases 10 and 32** (primes in A099165); complete below 10^40: 2, 3, 5, 7, 11, 19891, 10894711311749801, 13150353235305131, 17624701910742671, 16883812041873367076337814021838861
* **Primes that are palindromic in bases 10 and 36** (primes in A250412); complete below 36^26 (about 2.91*10^40): 2, 3, 5, 7, 11, 30803, 3315133, 111471214393412174111, 33567980330440491719404403308976533
* **Primes that are palindromic in bases 2 and 6** (primes in A182233); complete below 6^65 (about 3.80*10^50): 3, 5, 7, 85093, 129956651976247, 9027418797989294431, 34579344526412674075223, 540770568795997877440747003, 493729853550986902109639562979, 133957498046208178608172058508570751379
* **Primes that are palindromic in bases 2 and 7** (primes in A182234); complete below 7^48 (about 3.67*10^40): 3, 5, 107, 257, 1566399158893, 2059577471782805986901831467
* **Primes that are palindromic in bases 3 and 5** (primes in A259374); complete below 5^58 (about 3.47*10^40): 2, 1667, 635730140688854381713947926325539
* **Primes that are palindromic in bases 3 and 6** (primes in A259375); complete below 6^70 (about 2.96*10^54): 2, 78526951, 189403700304176275112370190723561706238418229
* **Primes that are palindromic in bases 3 and 8** (primes in A259381); complete below 3^84 (about 1.20*10^40): 2, 1554173, 9154421879941826706077226079
* **Primes that are palindromic in bases 4 and 6** (primes in A259376); complete below 6^65 (about 3.80*10^50): 2, 3, 5, 7517, 69313, 322509617, 6372205983137069544853
* **Primes that are palindromic in bases 4 and 7** (primes in A259378); complete below 7^48 (about 3.67*10^40): 2, 3, 5, 257, 424623193
* **Primes that are palindromic in bases 4 and 9** (primes in A259387); complete below 9^42 (about 1.20*10^40): 2, 3, 5, 373, 16319, 1181729, 22346657910193902193
* **Primes that are palindromic in bases 5 and 7** (primes in A249156); complete below 7^48 (about 3.67*10^40): 2, 3, 34853363, 4359505255951
* **Primes that are palindromic in bases 5 and 9** (primes in A259388); complete below 9^42 (about 1.20*10^40): 2, 3, 109, 701, 47653
* **Primes that are palindromic in bases 6 and 15** (primes in A249155); complete below 15^43 (about 3.73*10^50): 2, 3, 5, 7, 1627, 1777, 52201, 89731777, 4597931032867, 6953069365619265023329, 11397145769735915135407, 340081611711521780932476199, 355971209534794822452639919
* **Primes that are palindromic in bases 6 and 8** (primes in A259384); complete below 6^65 (about 3.80*10^50): 2, 3, 5, 7, 116755331881, 15813860880803, 6561870189686046211609, 36808172245371104292023
* **Primes that are palindromic in bases 6 and 9** (primes in A259389); complete below 9^57 (about 2.47*10^54): 2, 3, 5, 7, 191, 2295481, 61281793, 25237844535689857919, 7437875238052201437091, 70185808353033086807480636273678342496252011303749
* **Primes that are palindromic in bases 7 and 29** (primes in A249158); complete below 29^28 (about 8.85*10^40): 2, 3, 5, 13907, 1657740313, 2181959871705099791966033339

## Appendix A. Certified bounds per base pair

| Bases | Searched below | GPU hours logged |
|---|---|---|
| 2, 3 | 3^117 (about 6.66*10^55) | 98.54 |
| 2, 5 | 5^58 (about 3.47*10^40) | 0.01 |
| 2, 6 | 6^65 (about 3.80*10^50) | 0.05 |
| 2, 7 | 7^48 (about 3.67*10^40) | 0.01 |
| 2, 9 | 9^54 (about 3.38*10^51) | 8.66 |
| 2, 10 | 10^55 (about 1.00*10^55) | 11.48 |
| 3, 5 | 5^58 (about 3.47*10^40) | 0.09 |
| 3, 6 | 6^70 (about 2.96*10^54) | 0.01 |
| 3, 7 | 7^48 (about 3.67*10^40) | 0.05 |
| 3, 8 | 3^84 (about 1.20*10^40) | 0.01 |
| 3, 10 | 10^49 (about 1.00*10^49) | 15.55 |
| 4, 6 | 6^65 (about 3.80*10^50) | 0.06 |
| 4, 7 | 7^48 (about 3.67*10^40) | 0.02 |
| 4, 9 | 9^42 (about 1.20*10^40) | 0.02 |
| 4, 10 | 10^54 (about 1.00*10^54) | 5.11 |
| 5, 7 | 7^48 (about 3.67*10^40) | 0.34 |
| 5, 8 | 5^58 (about 3.47*10^40) | 0.04 |
| 5, 9 | 9^42 (about 1.20*10^40) | 0.04 |
| 5, 10 | 10^54 (about 1.00*10^54) | 0.00 |
| 6, 8 | 6^65 (about 3.80*10^50) | 0.06 |
| 6, 9 | 9^57 (about 2.47*10^54) | 0.01 |
| 6, 10 | 10^55 (about 1.00*10^55) | 16.46 |
| 6, 15 | 15^43 (about 3.73*10^50) | 0.08 |
| 7, 9 | 9^42 (about 1.20*10^40) | 0.11 |
| 7, 10 | 10^47 (about 1.00*10^47) | 6.00 |
| 7, 29 | 29^28 (about 8.85*10^40) | 0.16 |
| 8, 10 | 10^55 (about 1.00*10^55) | 13.92 |
| 9, 10 | 10^49 (about 1.00*10^49) | 7.02 |
| 10, 11 | 11^47 (about 8.82*10^48) | 7.34 |
| 10, 12 | 12^49 (about 7.58*10^52) | 8.09 |
| 10, 13 | 13^41 (about 4.70*10^45) | 2.80 |
| 10, 14 | 14^46 (about 5.27*10^52) | 6.03 |
| 10, 15 | 15^49 (about 4.25*10^57) | 0.27 |
| 10, 16 | 10^51 (about 1.00*10^51) | 0.61 |
| 10, 20 | 20^31 (about 2.15*10^40) | 0.00 |
| 10, 24 | 24^29 (about 1.06*10^40) | 0.00 |
| 10, 27 | 27^28 (about 1.20*10^40) | 0.10 |
| 10, 32 | 10^40 (about 1.00*10^40) | 0.00 |
| 10, 36 | 36^26 (about 2.91*10^40) | 0.00 |
| 11, 13 | 13^36 (about 1.26*10^40) | 0.04 |

## Appendix B. Primality certificates of the new prime terms

Pocklington's theorem: if N-1 = F*R with F > sqrt(N), and for every prime q dividing F some base a
satisfies a^(N-1) = 1 (mod N) and gcd(a^((N-1)/q) - 1, N) = 1, then N is prime. Here N-1 is factored
completely; prime factors above 10^12 are certified the same way (recursively), smaller ones by
trial division. Regenerate with `python3 tools/certify.py N`.

### A046478: 314419338131817706717607718131833914413 (submitted, pending review)

```
N = 314419338131817706717607718131833914413: N-1 fully factored, witnesses:
  q = 2^2  a = 2
  q = 727  a = 2
  q = 2939  a = 2
  q = 26573  a = 2
  q = 34471  a = 2
  q = 40162504844095109164597  a = 2
    N = 40162504844095109164597: N-1 fully factored, witnesses:
      q = 2^2  a = 2
      q = 3  a = 2
      q = 19  a = 2
      q = 188603  a = 2
      q = 933979507407119  a = 2
        N = 933979507407119: N-1 fully factored, witnesses:
          q = 2  a = 17
          q = 19  a = 2
          q = 2213  a = 2
          q = 11106375097  a = 2
```

