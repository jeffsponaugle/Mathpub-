# OEIS submission status and remaining items

Verified against oeis.org on 2026-10-03.

## Already in OEIS — nothing to do

**A045875** (revision #54, Sep 26 2026):
- DATA ends `..., 356677212, 1973572199, 3238684956` — **a(18) and a(19) are in**.
- Extensions present: `a(18) from Jeff Sponaugle, Sep 24 2026` and
  `a(19) from Jeff Sponaugle, Sep 26 2026`.

**A215732** ("first digit to appear n times in succession in a power of 2",
revision #46, Sep 26 2026):
- DATA ends `..., 4, 9` — the run digits for n=18 ('4') and n=19 ('9') are in,
  with matching extensions credited.

**A215733 / A215734 / A215735** are the analogous sequences for powers of
3, 5, and 6 — unaffected by these results. (The tools here adapt readily:
replace doubling with multiply-by-small-base; a possible future project.)

## Remaining submittable item: the a(20) lower bound

A045875 traditionally records bounds for the next unknown term as comments
(existing examples: `a(9) > 20000`, `a(11) > 250000`, `a(15) > 10297974`,
`a(17) > 107000000`). The matching new comment, in OEIS draft format:

```
a(20) > 4033865312. - _Jeff Sponaugle_, Oct 03 2026
```

**Before submitting, use the freshest bound.** The search is currently
paused; 4,033,865,312 is the proven contiguous frontier at the 2026-09-27
pause (see RESUME-A20.md). If the search has resumed, the live bound is
atom1's... more precisely the lowest machine frontier once the
[4,033,865,312 – 4,046,000,000] gap is closed; one command prints it:

```bash
ssh jbs@10.1.30.36 'tr "\r" "\n" < A045875/a20-*.log | grep steps/s | tail -1'
```

(Replace the number in the comment with that m only when coverage below it
is contiguous — after the gap leg completes; otherwise keep 4033865312.)

## Supporting details (if an OEIS editor asks)

- **a(18) = 1973572199**: 2^m has 594,104,431 digits and contains exactly
  18 consecutive '4's at digit offsets 394,539,910..394,539,927 from the
  least significant digit. Context: `...04265984511(4^18)35479973487...`
- **a(19) = 3238684956**: 2^m has 974,941,319 digits and contains exactly
  19 consecutive '9's at offsets 875,306,593..875,306,611 from the least
  significant digit. Context: `...73365556368(9^19)4732931778...`
- The exact run lengths (18, not 19; 19, not 20) are what justify the
  free bounds a(19) > a(18) and a(20) > a(19) at the time of each find.
- **Method**: incremental doubling of the full decimal expansion in base
  10^9 (carries are local, so the pass parallelizes), with run detection
  via the aligned-repdigit-limb filter (any run of >= 17 digits contains a
  limb equal to d*111111111), on NVIDIA GB10 GPUs batching 28 doublings
  per memory pass. Every exponent in [a(17), a(19)] (and onward for the
  a(20) bound) was searched in logged, checkpointed runs; each find was
  re-verified by an independent direct GMP computation of 2^m. Code and
  provenance: see README.md and RESULTS.md in this repository.

## Bookkeeping

- When the a(20) comment is accepted, consider whether the stale historical
  bound comments need touching — OEIS convention is to leave them (they
  document search history), so: no.
- When a(20) itself is found: submit the term + extension line for A045875,
  the run digit for A215732, and a fresh `a(21) > ...` comment, same
  pattern as above.
