# OEIS submissions — pi_missing results

Cross-check of the `pi_missing` scan results (first 1.03×10^13 decimal digits
of Pi) against the live OEIS entries, as of **2026-10-03**.

## Status summary

| Sequence | Definition (digits of Pi) | OEIS last term | Our result | Status |
|---|---|---|---|---|
| [A228988](https://oeis.org/A228988) | smallest missing number in first 10^n digits | a(13) = 100000070502 | matches | ✅ submitted & accepted |
| [A260627](https://oeis.org/A260627) | largest n-digit number missing in first 10^n digits | a(13) = 9999999999998 | a(11) **mismatch** | ⚠️ **a(11) typo in OEIS — edit needed** |
| [A153221](https://oeis.org/A153221) | k found at position k−3 | a(12) = 4819181870297 | matches | ✅ submitted & accepted |
| [A153223](https://oeis.org/A153223) | k found at position k−4 | a(14) = 9029718698652 | matches | ✅ submitted & accepted |
| [A153224](https://oeis.org/A153224) | k found at position k−5 | a(10) = 5667529025936 | matches | ✅ submitted & accepted |

All five entries now carry extension credits to _Jeff Sponaugle_ (Aug 29 /
Sep 19–20, 2026), including the a(8) = 99999991 correction to A260627.

## 1. ACTION: fix the A260627 a(11) typo

The live entry reads:

```
..., 999999998, 9999999999, 9999999998, 999999999998, 9999999999998
         a(9)        a(10)       a(11)         a(12)          a(13)
```

a(11) is printed as **9999999998 — a 10-digit number**. By definition a(n)
is an n-digit number, and our run (independently reproduced on two machines
and two digit files) gives:

```
A260627(11) = 99999999998        (= 10^11 − 2, eleven digits)
```

A `9` was dropped somewhere between the run output and the entry. Submit an
edit replacing `9999999998` with `99999999998` at position 11, with the
one-line justification: *"a(11) must be an 11-digit number; the previous
value had 10 digits. 99999999998 = 10^11 − 2 does not occur in the first
10^11 digits, while 99999999999 and each of 10^11−1's other neighbors above
it do occur."*

## 2. Easy wins: b-files

None of the five entries has a b-file (`%H` is empty except A228988's
digit-source link). A b-file is `n a(n)` per line, named `b<Anumber>.txt`,
uploaded with a short edit. Ready-to-paste contents:

**b228988.txt** (offset 1)
```
1 0
2 12
3 103
4 1001
5 10000
6 14523
7 106945
8 1001823
9 10007363
10 100023783
11 1000020346
12 10000031203
13 100000070502
```

**b260627.txt** (offset 1, with the a(11) fix)
```
1 8
2 96
3 997
4 9997
5 99997
6 999997
7 9999997
8 99999991
9 999999998
10 9999999999
11 99999999998
12 999999999998
13 9999999999998
```

**b153221.txt**: terms 1–12 = 51, 875, 62843, 242424, 4308765, 9710721,
24747689, 126987778, 44159961449, 73514522243, 249973441618, 4819181870297.
**b153223.txt**: terms 1–14 = 9, 233, 1614, 9218, 27755974, 81259258,
120526011, 238732548, 651265325, 783609697, 4927138295, 29478491757,
6885005566681, 9029718698652.
**b153224.txt**: terms 1–10 = 26, 32, 41, 86, 2799, 469113068, 2544423437,
7730967174, 870163762036, 5667529025936.

## 3. Worthwhile additions to existing entries

- **Digit-source link (`%H`)** on A260627 and the three A15322x entries,
  matching the style of A228988's Peter Trüb link — e.g. a link to the
  Google Cloud 100-trillion-digit Pi computation the 10T file came from.
- **Search-bound comments**: for the self-locating sequences, a comment of
  the form *"No further terms below 1.03×10^13. - Jeff Sponaugle"* tells
  future extenders where the frontier is. Same for A228988/A260627 (next
  term requires 10^14 digits).
- **A260627 correction note**: the entry silently carries the a(8)
  correction in its `%E` line; a `%C` comment with the evidence (99999992
  occurs at position 3389380; 99999991 first occurs at position 209136149)
  would preempt future "discrepancy" reports against old literature.

## 4. Extension opportunities (same tool, new runs)

### Self-locating family, offsets −1 / 0 / +1 — needs a 6-line code change

Three classic "k at position k" sequences use the same machinery as
A153221/223/224; extending `g_selfseq[]` in `pi_missing.c` with offsets
−1, 0, +1 (counter k = i, i+1, i+2 at 0-based position i) covers them all
in one `-S` rerun:

| Sequence | Convention | Our offset | OEIS last term | 10^13 outlook |
|---|---|---|---|---|
| [A064810](https://oeis.org/A064810) | post-decimal digits 0-based | −1 | 8773143366618 (8.8×10^12) | frontier nearly reached; maybe one term |
| [A057680](https://oeis.org/A057680) | "3" counted as position 0 | 0 | 656430109694 (6.6×10^11) | ~1.2 decades open, ≈3 expected terms |
| [A057679](https://oeis.org/A057679) | "3" counted as first digit | +1 | 649661007154 (6.5×10^11) | ~1.2 decades open, ≈3 expected terms |

(Offset +2, "k at position k−2", has no OEIS entry — a candidate **new**
sequence that would complete the family from k−5 through k.)

### Other constants — A228988's siblings

The tool is constant-agnostic: any digit file works.

- **[A228989](https://oeis.org/A228989)** — smallest missing number in
  digits of **e**: stops at a(9); a y-cruncher run of e to 10^12–10^13
  digits extends it by 3–4 terms.
- **[A228990](https://oeis.org/A228990)** — same for **phi**: stops at a(8).
- "Largest n-digit missing" analogs of A260627 for e/phi do not exist —
  candidate new sequences, generated by the same `-L` pass.

### Longer Pi runs

A228988(14)/A260627(14) need 10^14 digits (~100 TB of file, ~1.5×10^12-wide
bit array ≈ 190 GB with `-m 1.5e12` — feasible on the 512 GB server; wall
time is the constraint).

## 5. Evidence appendix (for referee questions)

- Tool: `pi_missing.c` (this repo) — single-pass bit-array scan; windows
  never cross the boundary they are counted for; refuses to answer rather
  than guess when a tracking window saturates.
- Verification: reproduces every previously known term of all five
  sequences; digit file cross-checked byte-for-byte against an independent
  Chudnovsky binary-splitting computation (first 10^7 digits); key
  claims re-verified with independent substring search (`grep`, `pifind`).
- A260627(8) correction evidence: `99999992` at position 3389380;
  `99999991` absent before position 209136149; `99999993`–`99999999` all
  present within the first 10^8 digits (first occurrences at 1722776,
  14593770, 36356643, 20202676, 78088568, 15256174, 36356642).
- Self-locating matches carry their positions in the run logs
  (`missing-10t.txt`), each satisfying position = k − offset exactly.
