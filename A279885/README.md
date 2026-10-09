# pipal — A279885 search tool

Finds the terms of [OEIS A279885](https://oeis.org/A279885) ("First n-digit
palindrome in the decimal expansion of Pi") by streaming a file of pi digits.

## Build

```
make
```

(or `cc -O2 -o pipal pipal.c`)

## Usage

```
pipal [-o outfile] [-n maxlen] pifile
```

- `pifile` — decimal digits of pi. Non-digit bytes (the `3.` decimal point,
  newlines, spaces) are ignored, so `3.14159...` text files, wrapped-line
  files, and raw y-cruncher output all work. Designed for very large files
  (tested logic-equivalent down to a 200-byte window; memory use is constant
  at ~12 MB regardless of file size).
- `-o outfile` — write results there instead of stdout (lines are flushed as
  each term is found, so partial results survive an interrupted run).
- `-n maxlen` — longest palindrome length to track (default 64, max 128).
  If every length ≤ maxlen is found, the scan stops early.

Output is TSV: `n <TAB> start_index <TAB> palindrome`. Indices are 0-based
with the leading `3` as digit 0 (so `a(1) = 3` at index 0). Terms are digit
strings and may start with `0` (e.g. `a(10) = 0136776310`), per the OEIS
comment. Progress is reported to stderr every 10^9 digits when stderr is a
terminal.

## Algorithm

Single-threaded, one sequential pass. For every palindrome center (odd and
even), it expands outward while digits match, capped at maxlen/2. Since
random digits mismatch 90% of the time, expected work is ~2.2 comparisons
per digit independent of file size. Centers are processed left to right, so
the first center whose maximal palindrome reaches a still-missing length n
gives the first n-digit palindrome. Measured ~300 M digits/s from page
cache; a 10-trillion-digit scan is roughly 9 CPU-hours and will typically be
I/O-bound.

Verified against OEIS data a(1)–a(13) (values and positions) on the first
10^6 digits of pi, cross-checked with an independent brute-force script.

## Sequence status (OEIS entry as of revision 65, Sep 2026)

Known terms:

| n  | a(n)                    | n  | a(n)                    |
|----|-------------------------|----|-------------------------|
| 1  | 3                       | 12 | 450197791054            |
| 2  | 33                      | 13 | 9475082805749           |
| 3  | 141                     | 14 | 95715977951759          |
| 4  | 3993                    | 15 | 697763767367796         |
| 5  | 46264                   | 16 | 4855796336975584        |
| 6  | 999999                  | 17 | 42611063036011624       |
| 7  | 1736371                 | 18 | 700015327723510007      |
| 8  | 23911932                | 19 | 5851405951595041585     |
| 9  | 398989893               | 20 | 74054579700797545047    |
| 10 | 0136776310              | 21 | 483737413000314737384   |
| 11 | 21348884312             |    |                         |

- a(6) = 999999 is the Feynman point. a(10) begins with a leading zero.
- a(16)–a(20) were found by Michael S. Branicky (Jan 2022) within the first
  10^9 digits; the entry notes no further terms lie wholly within 10^9
  digits.
- a(21) = 483737413000314737384 was added by Jeff Sponaugle, Sep 19 2026.
- Cross-references: A000796 (digits of Pi), A002113 (palindromes), A226486,
  A280631.

## What could be found next

a(22) is the frontier. Treating pi's digits as random, a window of length n
is palindromic with probability 10^-floor(n/2), so a scan of N digits sees
about N · 10^-floor(n/2) n-digit palindromes, and the first one is expected
near digit position 10^floor(n/2). For a 10-trillion-digit (10^13) scan:

| n      | expected count in 10^13 digits | chance ≥ 1 exists | expected first position |
|--------|--------------------------------|--------------------|-------------------------|
| 22, 23 | ~100 each                      | ~certain           | ~10^11                  |
| 24, 25 | ~10 each                       | ~99.995%           | ~10^12                  |
| 26, 27 | ~1 each                        | ~63% each          | ~10^13                  |
| 28, 29 | ~0.1 each                      | ~10% each          | ~10^14                  |

So a full 10T-digit run should deliver a(22)–a(25) with near certainty,
has roughly even odds on each of a(26) and a(27), and would be lucky to see
anything longer. Note the first-occurrence positions are only known once
every earlier window has been checked, so terms must come from a scan that
covers all digits from the start — which is exactly what this tool does in
one pass. Pushing past a(27) will likely take ~10^14–10^15 digits of pi,
beyond any currently published computation (~3 × 10^14 digits as of 2025).
