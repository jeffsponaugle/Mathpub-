# OEIS submissions

Status as of 2026-10-06, checked against the live OEIS entries and their
draft (pending-edit) pages. Files ready to upload are in `A171702/`,
`A171703/`, `A171740/` and `A216/bfiles/`. All submissions are made by Jeff
Sponaugle.

## At a glance

| Sequence | Item | Status |
|---|---|---|
| A171701 | a(41) corrected 30160 → 20160; b-file n = 1..1152 | **live** (rev 19, Sep 30) |
| A171741 | a(3) corrected 154 → 121; a(8) upper bound | **live** (rev 28, Oct 2; editors also prepended a(1) = 1) |
| A171705 | a(4) ≤ 40968^5 | **live** (rev 20, Oct 4) |
| A171706 | a(4) ≤ 327696^6 | **live** (rev 21, Oct 1) |
| A171703 | a(10)–a(15) and the "confirmed by exhaustive search" comment | **live** (rev 16, Oct 1) |
| A171703 | a(16)–a(22), cube comment, b-file n = 1..22 | **pending review** (proposed Oct 1) |
| A171704 | a(5), a(6), and the a(7) > 1.244×10^17 comment | **pending review** (proposed Sep 30) |
| A171740 | b-file "n = 1..37" | **live but needs a fix**: the file skips n = 36, and its header has internal notes (section 1) |
| A171702 | b-file extension to n = 223 | **ready, not submitted** |
| A171742 | a(6) ≤ 40968^5 comment | **ready, not submitted** |
| A216902–A216910 | b-files for 9 base-gap sequences | **ready, not submitted** |
| new | length-16 base-gap sequence | **ready, not submitted** |

## 1. Needs attention: the A171740 b-file

The b-file published in rev 58 (Sep 29) is labeled "n = 1..37" but has rows
1..35 and 37, with no row 36, because a(36) is still unknown. OEIS b-files
need consecutive indices. Its header also carries internal working notes
("see repo README", "Pending for next submission: a(32) (final certificate
completing) ...", "(in search — gates the consecutive run)").

Suggested fix:
- **Replace the b-file** with `A171740/b171740.txt`: rows n = 1..35 only,
  with a clean three-line header. Link text: "Table of n, a(n) for n = 1..35".
- **Move a(37) into a comment** until a(36) is known: "a(37) =
  2405399997091326396338729111112302450399341, a 37-digit palindrome in bases
  14 and 15. - Jeff Sponaugle"
- **Add the a(36) lower bound**: "a(36) > 28^35 = 4.47*10^50: exhaustive
  search of all base pairs with larger base <= 27 (stages 23-27). - Jeff
  Sponaugle". Stage 28 is part-searched with no match; the paused search's
  resume points are in the main README, section 6.

## 2. Ready to submit

### A171702 (3-digit palindromes in >= n bases): b-file extension
- **New b-file** `A171702/b171702.txt`, n = 1..223. The current b-file
  (Resta, 2016) has n = 1..100; its terms are kept unchanged, with their
  credits, in the header.
- a(101)-a(223) come from an exhaustive CPU sieve of every 3-digit palindrome
  in every base for 13967553601 <= v < 10^12 (4.39*10^12 palindromes). They
  are valid terms because a(n) >= a(100) = 13967553601.
- Validation: the same program reproduces all 100 existing terms, and
  independent base counts confirm spot checks (a(101) = 16908091201 has 101
  bases; a(223) = 963761198401 has 223).
- Suggested comment: "b-file extended to n = 223 by exhaustive search to
  10^12. - Jeff Sponaugle"

### A171742 (first 4-fold intrinsically n-palindromic number): upper bound
- Its a(6) is the same number as A171705(4), whose entry now carries the bound.
- Suggested comment: "a(6) = A171705(4) <= 40968^5 = 115404784440424844525568,
  a 6-digit palindrome in bases 10241, 13655, 20483, 40967. - Jeff Sponaugle"

### Base-gap family: b-files
| Sequence | Length | Terms in OEIS now | New b-file |
|---|---|---|---|
| A216902 | 8 | 14 | `A216/bfiles/b216902.txt`, n = 2..1001 |
| A216903 | 9 | 10 | `A216/bfiles/b216903.txt`, n = 1..1000 |
| A216904 | 10 | 19 | `A216/bfiles/b216904.txt`, n = 2..1000 |
| A216905 | 11 | 14 | `A216/bfiles/b216905.txt`, n = 1..1000 |
| A216906 | 12 | 7 | `A216/bfiles/b216906.txt`, n = 2..200 |
| A216907 | 13 | 11 | `A216/bfiles/b216907.txt`, n = 1..200 |
| A216908 | 14 | 4 | `A216/bfiles/b216908.txt`, n = 2..200 |
| A216909 | 15 | 5 | `A216/bfiles/b216909.txt`, n = 1..200 |
| A216910 | 17 | 3 | `A216/bfiles/b216910.txt`, n = 1..60 |

Every term is independently re-checked as an L-digit palindrome in both
bases, with the right gap, and all previously known terms match OEIS
(1,184 known terms across the family, including Chai Wah Wu's b-files for
lengths 2-7).

### New sequence: length 16
- **Name**: Smallest palindromic number of length 16 in two bases differing by n.
- **b-file**: `A216/bfiles/b_new_L16.txt`, n = 2..60.
- **Data**: 530386561769238496, 3949277000413367696400, 235846132230130855735800,
  773759843075568329413720, 2100774733130395290456684,
  297357478494888295690227324, 4325221042711488170326944,
  73420354736530342140173040, 24039858092401995627330800,
  1164464285050781977724480388, 1162846680184119625009479424,
  611393260440971283573173652 (n = 2..13)
- **Offset**: 2. For n = 1 no solution exists: an even-length palindrome in
  base b is divisible by b+1, so in base b+1 it would end in digit 0.
- **Comment**: Fills the length-16 gap in the collection A216840, A216841,
  A216843, A216899-A216910 (the comment in A216910 notes that only the
  n = 2 term was known).
- **Example**: a(2) = 530386561769238496, the 16-digit palindrome in bases 13
  and 15, which is also A171740(16).
- **Crossrefs**: A171740, A216840, A216841, A216843, A216899-A216910.
- **Keywords**: nonn, base.

## 3. Pending review (submitted; awaiting editors)

### A171703 (4-digit palindromes in >= n bases), proposed Oct 1
The draft adds a(16) onward to the data, the cube comment, and the b-file
n = 1..22 (uploaded as `b171703_1.txt`).

| n | a(n) | Bases |
|---|---|---|
| 16 | 35158608576000000 = 327600^3 | 16, all of the form 327600/d - 1 |
| 17 | 40254862491648000 = 342720^3 | 17: 16 of that form plus the sporadic base 43349 |
| 18 | 205045005215232000 = 589680^3 | 18 |
| 19, 20 | 374368864117248000 = 720720^3 | 20 |
| 21 | 1168773800173056000 = 1053360^3 | 21 |
| 22 | 2994950912937984000 = 1441440^3 | 22 |

If editors ask for proof details:
- Exhaustive search of every 4-digit palindrome in every base below each
  term: CPU sieve below 8.4*10^15, then a contiguous GPU sieve to 1441440^3.
- Independent CPU sieves reproduce a(16)-a(18) (1.12*10^13 palindromes) and
  a(19)-a(21) (3.34*10^13 palindromes) exactly.
- An independent CPU sieve of a(21)..a(22) (4.71*10^13 palindromes) also
  reproduces a(22) exactly, so all 22 terms are confirmed on both CPU and
  GPU sieves.
- Every term's base count is re-checked separately in Python.
- a(23) <= 1965600^3 = 7594259452416000000 by the construction; its search is
  paused.

### A171704 (5-digit palindromes in >= n bases), proposed Sep 30
The draft adds a(5) and a(6) and the a(7) lower bound.

- a(5) = 4922057407205376 = 8376^4, a 5-digit palindrome in exactly 5 bases:
  2093, 2791, 4187 (each 8376/d - 1 with digits d^4 x (1,4,6,4,1)), 8375
  (1 4 6 4 1) and the sporadic base 5906.
- a(6) = 34777153514704896 = 13656^4, in exactly 6 bases: 3218, 3413, 4551,
  4828, 6827, 13655.
- a(7) > 124389107494560000.

If editors ask for proof details:
- Every 5-digit palindrome in every base below 1.244*10^17 was searched (CPU
  below 10^12, GPU above).
- Independent CPU sieves reproduce a(5) (3.7*10^12 palindromes) and a(6)
  (1.43*10^13 palindromes) exactly.

## 4. Already live (for the record)
- **A171701** rev 19: a(41) = 20160; b-file n = 1..1152 (exhaustive to 10^10).
- **A171741** rev 28: a(3) = 121 (232 in base 7, 171 in base 8, 121 in base
  10); a(8) <= 229644^7. Editors prepended a(1) = 1, so the offset is now 1.
- **A171705** rev 20: a(4) <= 40968^5 = 115404784440424844525568.
- **A171706** rev 21: a(4) <= 327696^6 = 1238302761306332849871154940215296.
- **A171703** rev 16: a(10)-a(15), all exactly Yamanouchi's 2014 bounds.
- **A171740** rev 58: b-file (see section 1 for the needed fix).
