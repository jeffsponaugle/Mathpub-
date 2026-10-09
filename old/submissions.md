# OEIS submission opportunities — status as of 2026-10-03

Cross-check of every sequence touched by the tools in `~/src/math/`
(`A057330/cchain`, `A119291/pdigits`) against the live OEIS entries.

---

## 1. Ready to submit now — Cunningham chain family

All three stem from Nenad Micic's March 2026 results published on Dirk
Augustin's records page (http://primerecords.dk/Cunningham_Chain_records.htm),
which OEIS has not yet absorbed. Credit Micic as discoverer; our `cchain`
tool independently verified every chain (50-round BPSW per member).

### A057330 — smallest prime starting a CC of the 2nd kind, length >= n
- **OEIS now**: 16 terms, last modified Apr 2023. a(17) missing.
- **Submit**: a(17) = 8058978143544782179441
- Basis: found AND proven minimal by Nenad Micic, 2026-03-06 (Augustin's
  records page, Part 2 + NEWSFLASH). 22 digits.
- Our verification (`cd ~/src/math/A057330 && ./cchain --verify 8058978143544782179441`):
  members m_0..m_16 all prime, m_17 composite — chain length exactly 17.
- Suggested %E: `a(17) found and proved minimal by Nenad Micic, Mar 06 2026.`

### A005603 — smallest prime starting a COMPLETE CC of the 2nd kind, length n
- **OEIS now**: 16 terms, last modified Nov 2025. a(17) missing. (A 2018
  editor comment already suggests adding terms from Augustin's page.)
- **Submit**: a(17) = 8058978143544782179441 (same chain as A057330)
- Completeness verified by us: forward, m_17 = 2*m_16 - 1 is composite;
  backward, (p+1)/2 = 4029489071772391089721 is composite — so the chain
  is complete with length exactly 17. Minimality: any complete length-17
  chain is a length->=17 chain, and Micic proved this is the smallest of
  those, so it is also the smallest complete one.

### A057331 — smallest p whose first n iterates of x->2x+1 are all prime
  (i.e. CC of the 1st kind, length >= n+1)
- **OEIS now**: 16 terms (through a(15), chain length 16), last modified
  Oct 2025. Prior terms were explicitly "added from A005602".
- **Submit**: a(16) = 2759832934171386593519
  - = A005602(17) (Wroblewski's chain, proven minimal by Micic Mar 2026;
    already in OEIS since Aug 2026). Chain length exactly 17 (we verified:
    m_17 composite), so first 16 iterates prime.
- **Likely also submittable**: a(17) = 776208125345634679522109
  - = A005602(18), added to OEIS Aug 2026 by Jonathan Pappas. Same
    derivation (chain length >= 18 => first 17 iterates prime). We
    verified the chain 2026-10-03: length exactly 18 (m_18 composite).
    Before submitting, confirm A005602's a(18) provenance traces to a
    proven-minimal complete CC18 of the 1st kind (check the entry's
    links/comments or Augustin's page for the proof claim).

### Not submittable in this family
- **A005602** (complete, 1st kind): a(17)-a(18) already added by Jonathan
  Pappas, Aug 22 2026. Done by others.
- **A057330 a(18)**: would need the smallest CC2nd of length >= 18 proven
  minimal; only unproven candidates exist (~1.17e26). Out of reach.

---

## 2. Already submitted and published — digit-count family (ours)

### A119291..A119300 — digit d (0..9) occurrences in the first 10^n primes
- **OEIS now**: 15 terms each.
  - a(13)-a(14): "_Jeff Sponaugle_, Aug 20 2026"
  - a(15): "_Jeff Sponaugle_, Sep 20 2026"
- Nothing pending. Reports preserved: `A119291/report13.txt`,
  `report14.txt` (a(15) run lived on the other machine).
- Optional polish: upload explicit b-files (n=1..15) for each of the ten
  entries if any lack one.

### Future term: a(16) for all ten
- Requires the N=16 campaign: first 10^16 primes, sieve to ~3.95e17.
- Tooling is ready (`pdigits -n 16 -s 12`, multi-machine `--chunk-min/max`
  split, `--benchmark` for fleet planning). ~180 days solo on the M1 Max
  (0.656 Gprimes/s); split across laptop + other machine + DGX Sparks
  (CPU) plausibly ~6 weeks. Last N with the independent A119290 digit-sum
  check.

---

## 3. Verified, nothing to submit (reference/cross-check sequences)

| Sequence | State | Our relationship |
|----------|-------|------------------|
| A091634..A091643 (primes < 10^n lacking digit d) | 19 terms each | independently re-verified n=1..15 during our runs; extending needs enumeration to 10^20 — out of reach |
| A119290 (total digits, first 10^n primes) | 22 terms (analytic) | used as independent check of our new terms; known far beyond our scope |
| A006988 (10^n-th prime) | known to n=18 | boundary check; nothing to add |
| A006880 (pi(10^n)) | known far beyond 10^18 | sieve-integrity check; nothing to add |

---

## 4. Submission mechanics reminder

- Draft edits at https://oeis.org/ (login required), one sequence per edit.
- For the Cunningham terms: cite Augustin's records page as the reference,
  credit Micic in the %E line, and mention independent verification
  (BPSW, 50 rounds) if asked. For a published announcement-grade check,
  certify members with PARI/GP `isprime` (actual proof at these sizes —
  22-28 digit members).
- A057331 terms are derivations from A005602; OEIS precedent is to note
  exactly that ("added from A005602"), as prior %E lines there do.
