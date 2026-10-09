# OEIS submission opportunities — Pi prime search family

Cross-checked against live OEIS entries on 2026-10-03.

## Status of our results vs. OEIS

### A076106 (latest-appearing n-digit prime in Pi) — UP TO DATE

OEIS data (revision 22, 2026-08-29):
`7, 73, 373, 9337, 35569, 805289, 9271903, 43927427, 342263843, 7530425897, 70676925127`

Matches our computed terms exactly through a(11). Already present:
- Extension credit: "a(8)-a(11) from _Jeff Sponaugle_, Aug 29 2026"
- Comment: "a(12) estimated to need ~26 trillion digits of Pi. - _Jeff
  Sponaugle_, Aug 29 2026"
- Link: "Jeff Sponaugle, Verification for a(8)-a(11)" (/A076106/a076106.txt)
- Keywords still `hard,more,nonn,base` (correct: a(12) remains open).

### A076130 (position where that prime appears) — UP TO DATE

OEIS data (revision 30, 2026-08-29):
`13, 299, 5229, 75961, 715492, 11137824, 135224164, 1541659153, 20252853413, 238602469714, 2192839682323`

Matches our positions exactly through a(11). Same extension credit,
a(12)-estimate comment, and verification a-file as A076106.

**Conclusion: the a(8)-a(11) term submissions are done and approved.
No b-files needed for either** — at 11 terms, the OEIS-synthesized b-file
from the data lines suffices (explicit b-files that merely duplicate the
data are discouraged).

## Remaining submission opportunities

### 1. The C program (A076106 + A076130) — easy, high value

Both entries still carry only the 2021 Python program, which does a
`str.find` per prime over an in-memory digit string — impractical beyond
a(6)/a(7). Our `a076106.c` (streaming scan, segmented multithreaded sieve,
checkpointing) is what actually produced a(8)-a(11).

Options, in OEIS-preferred order:
- Upload the program as an a-file on A076106 (e.g. `a076106_1.c.txt`) via
  the entry's edit form, with a LINKS line:
  `Jeff Sponaugle, <a href="/A076106/a076106_1.c.txt">C program</a>`
  and a one-line PROG note: `(C) See Links section.` Cross-reference the
  same file from A076130's LINKS.
- Or host it (e.g. GitHub) and add an external link — OEIS prefers
  attached a-files since external links rot.

### 2. b-file for A076129 (position of FIRST n-digit prime in Pi) — easy gap-fill

Asymmetry found while cross-checking the family:
- A076094 (first n-digit prime) has a b-file to n=250 (Sean A. Irvine).
- A076129 (its position) has NO b-file and only 64 data terms.

A `b076129.txt` for n = 1..250 — the positions corresponding exactly to
Irvine's primes — is a clean, uncontroversial contribution. Cheap to
compute: scan the first few thousand digits of Pi with a provable
primality test (Pari/GP `isprime`, not just a pseudoprime test, since
candidates reach 250 digits). Format, one pair per line:

```
1 4
2 2
3 7
...
```

While at it, both A076094 and A076129 could be extended past n=250
(say to 500) in matching b-files; primality proofs of ~500-digit numbers
are still routine (APR-CL/ECPP).

### 3. a(12) / A076130(12) — the flagship open item

- Our own estimate comment (~26T digits) is already in both entries.
- Requirements: ~280 GB RAM (fits the >1TB servers), ~26T digits expected
  (the 100T file covers it, ~40T to be safe at 99%).
- Blocker is scan time: ~25 days single-threaded at the observed
  ~12M digits/s. Build the parallel scan first (worker blocks + ordered
  commit thread, est. 8-15x) to bring it to ~2-3 days.
- On completion: submit a(12) to both entries, update/remove our estimate
  comments, keep `more` (a(13) would need ~3*10^14 digits).

### 4. Possible new sequences: analogs for other constants — speculative

No OEIS sequence exists for "latest-appearing n-digit prime in the digits
of e" (searched 2026-10-03; the A076106 family is Pi-specific, from
Rivera's Puzzle 40). Our tool works on any digit stream, so e / sqrt(2) /
golden-ratio analogs (a value+position pair each) are computable to n~10
from existing digit files. Caveat: OEIS editors can be lukewarm on
"same idea, different constant" submissions — worth proposing only with a
few terms more than trivially computable (say through n=9-10), which we
can produce in a day per constant on the server.

## Adjacent family — checked, nothing for us to add

- A036903 / A080597 (scan until ALL n-digit strings seen): already
  extended through n=13 (2.94*10^14 digits) by others; beyond our data.
- A047658 (initial k digits of Pi's fractional part form a prime):
  pi-prime hunting via large-number primality proving, different problem.

## Suggested order

1. C program a-file upload to A076106/A076130 (minutes of effort).
2. b076129.txt for n=1..250 (an evening; needs Pari/GP or similar).
3. Parallel scan implementation, then the a(12) run on the 100T file.
4. Optionally, the e/sqrt(2) analogs as new sequences.
