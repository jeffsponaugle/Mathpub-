# a171740 — searching OEIS A171740

[OEIS A171740](https://oeis.org/A171740): **a(n) = the first number that is
palindromic with exactly n digits in more than one base.**

Examples: a(2) = 8, which is `22` in base 3 and `11` in base 7. a(4) = 624,
which is `4444` in base 5 and `1551` in base 7.

The OEIS entry (as of this writing) lists terms through a(17). The author's
2012 comment says he was stuck "waiting on the doubly 18-palindrome" — a(18)
was unknown — while already knowing a 19-digit candidate that could not be
added until a(18) was found. This tool was written to compute the missing
terms — **and it has: the sequence now extends through a(37)** (full table
in §4; a(32) is sealed, and a(36), the one remaining gap, is in its final
search stages), including
a(18) (open since 2009) and terms past 10^46 — far beyond 128-bit
arithmetic. Along the way came a structural theorem (§1.3) explaining why
even-indexed terms are so much harder than odd ones, a meet-in-the-middle
algorithm (§2.4) that cut the search exponent from W^(1/2) to W^(1/4), and
a 192-bit value type once the terms outgrew u128. All computed terms,
verified, are in `b171740.txt` (OEIS b-file format).

```
make
./a171740 -n 18 -t 14          # reproduce the a(18) discovery (~90 s)
./a171740 -n 19 -L 6411682614162861147   # reproduce the a(19) minimality proof (~2 s)
./a171740 -n 20 -t 14          # reproduce the a(20) discovery (~1 h from scratch)
```

- `-n N` — digit count (required)
- `-t T` — worker threads (default: all cores)
- `-L X` — search only values `< X` (exclusive). Used to confirm a known
  candidate: pass candidate+1; if the program finds exactly that value,
  nothing smaller exists.
- `-B B` — resume at stage b2 = B, asserting (without checking) that all
  earlier stages were already searched and found empty.
- `-H D` — with `-B`, resume mid-stage: D is the done-counter (the first
  number in parentheses) pasted from the last status line of the
  interrupted run. Because that counter sums completed work blocks rather
  than tracking a contiguous frontier, up to threads×block candidates
  below it may have been in flight when the run died; the program
  subtracts that margin automatically, so paste the number verbatim.
  Like `-B`, it asserts the skipped range was genuinely scanned — resuming
  past a would-be match makes the program report the *next* double
  palindrome instead of a(n). Works in MITM mode too: the counter resets
  per partner phase there, so on a multi-partner stage add `-P b1` naming
  the partner phase the counter came from (the last `MITM b1=X:` banner
  before the interruption); partner phases before it are assumed complete.
- `-q` — suppress periodic status lines

Status, progress (%done, rate, stage and overall ETA), and match
announcements go to **stderr**; the final result goes to **stdout**.

---

## 1. The mathematical structure

### 1.1 What "n digits in base b" pins down

A number N has exactly n digits in base b iff

```
b^(n-1) <= N < b^n
```

So for a *fixed* N and digit count n, the usable bases lie in a narrow band:
roughly `b ≈ N^(1/n)`. This is the key observation that makes the problem
finite and small: we are not looking for palindromes in "any two bases", but
in two bases that both squeeze N into an n-digit window.

### 1.2 Which base pairs are even possible

Suppose N is an n-digit palindrome in bases b1 < b2. Both windows must
overlap:

```
b2^(n-1) <= N < b1^n      =>      b2^(n-1) < b1^n
```

Taking logs: `b2 < b1^(n/(n-1))`. For large n the exponent tends to 1, so
**b2 must be very close to b1**. Concretely for n = 18 the condition
`(b1+1)^17 < b1^18` first holds at b1 = 8, and pairs with gap 2
(`(b1+2)^17 < b1^18`) first appear at b1 = 15. So the candidate pairs for
n = 18 are (8,9), (9,10), (10,11), ..., then (13,15), (14,15), (14,16), etc.
For small n the pairs can be far apart — a(2) = 8 uses bases 3 and 7 — but
then the numbers involved are tiny.

Each pair (b1, b2) defines a **search interval** of values:

```
[ b2^(n-1),  b1^n )
```

and any doubly-palindromic value for that pair must lie inside it.

### 1.3 The even-n theorem: adjacent bases are impossible

An even-length palindrome in base b is always divisible by b+1: modulo
b+1 we have b ≡ −1, so v ≡ Σ dᵢ(−1)ᵢ, and mirrored digit pairs (positions
i and n−1−i) cancel when n is even.

Now suppose n is even and v is an n-digit palindrome in both b1 and
b1+1. Palindromicity in base b1 forces (b1+1) | v. But v's *last* digit in
base b1+1 equals its *first* digit (it is a palindrome there too), which is
nonzero — so v is **not** divisible by b1+1. Contradiction:

> **For even n, no number is an n-digit palindrome in two adjacent bases.**

This is visible in the data — every known even-n term uses bases with gap
≥ 2 — and it has real consequences for the search: for even n all
adjacent-only stages are skipped entirely (raising the proven lower bound
for a(n) before any scanning), and within multi-partner stages the
adjacent partner is dropped. For n = 20, the first legal pair is (14, 16),
so a(20) ≥ 16^19 right out of the gate.

The divisibility also powers the even-n scanner below: candidates
enumerated in base b2 satisfy (b2+1) | v automatically, and the *extra*
requirement (b1+1) | v can be solved in closed form instead of tested.

### 1.4 Enumerating palindromes by their half

An n-digit palindrome in base b is completely determined by its leading
`m = ceil(n/2)` digits (the "half" h, an m-digit base-b number, leading
digit nonzero). With `k = floor(n/2)`:

- n even (m = k): `v(h) = h * b^k + rev(h)` where `rev` reverses all k
  digits of h.
- n odd (m = k+1): `v(h) = h * b^k + rev(h div b)` — the middle digit is
  h's last digit and is not mirrored.

Two facts make this enumeration ideal:

1. **Count.** There are only `(b-1) * b^(m-1)` n-digit palindromes in base
   b — about the *square root* of the interval width. An interval of 1e18
   values contains only ~1e9 palindromes of a given base.
2. **Monotonicity.** Since `rev(...) < b^k`, we have
   `h * b^k <= v(h) < (h+1) * b^k`, so v(h) is *strictly increasing in h*.
   Scanning halves in order scans values in order, which is what lets us
   stop the moment the current value exceeds the best match found.

### 1.5 The search strategy

For each pair we enumerate palindromes of the **larger** base b2 (fewer of
them per unit of value — count scales as `range / b2^k`) and test each value
for n-digit palindromicity in b1.

To find the *first* match across all pairs, pairs are processed grouped into
**stages**: stage b2 handles all valid partners b1 < b2 at once. Since a
stage's interval starts at `b2^(n-1)`, which is strictly increasing in b2,
processing stages in order b2 = 3, 4, 5, ... visits interval lower bounds in
increasing order. The loop terminates when the next stage's lower bound is
at or above the best match found so far (or the `-L` limit): no smaller
match can exist beyond that point. Within and across stages the running
best acts as a moving ceiling — every scan breaks as soon as its (monotone)
value stream reaches it — so the total work is bounded by the number of
candidate palindromes *below the answer*, and the minimum found is exact.

---

## 2. Implementation

Single C file, `a171740.c`. 128-bit arithmetic (`unsigned __int128`) for
values, 64-bit for halves; pthreads for parallelism.

### 2.1 Stage setup (`stage_setup`)

For stage b2:

- Collect all b1 in [2, b2) with `b1^n > b2^(n-1)` (non-empty interval),
  recording `p_lo = b1^(n-1)` and `p_hi = b1^n` for each. Powers are
  computed with saturation at 2^128−1 so overflow degrades safely.
- The stage's value window is `[b2^(n-1), min(max_b1 b1^n, limit))`.
- The half range is `h ∈ [b2^(m-1), h_end)` with
  `h_end = min(b2^m, hi/b2^k + 1)`; the per-candidate `v >= stop` check
  handles the exact boundary.

If the stage has exactly one partner and it is `b2 - 1` (the only shape that
occurs for large n), a specialized fast scanner is used; otherwise a generic
scanner handles the multi-partner / small-n cases.

### 2.2 The fast scanner (`scan_fast`)

This is where nearly all time is spent, so the hot loop contains **no
divisions at all**. Three layers of incremental computation:

**(a) Palindrome construction by parts.** Split the half as
`h = t*b2 + d0` (t = all but the last digit, d0 = last digit). Then with
`r = rev(t's m−1 digits)` computed once per t:

```
v(h) = t*b2*mulk + r + d0*step
  where mulk = b2^k
        step = mulk            (n odd:  d0 is the unmirrored middle digit)
        step = mulk + b2^(k-1) (n even: d0 also appears mirrored)
```

So the inner loop over d0 advances v by a single 128-bit addition, and the
reverse computation (m−1 divisions) is amortized over b2 candidates.

**(b) Incremental first/last-digit filter in b1.** A palindrome must have
equal first and last digits, and a random value fails this with probability
(b1−1)/b1 — a ~92% rejection rate for b1 = 12 at the cost of a couple of
adds. Both digits are maintained incrementally as v steps by `step`:

- last digit `bm = v mod b1`: maintained as `bm += (step mod b1)`, with a
  conditional subtract of b1 — no division.
- first digit `lead = v div b1^(n-1)`: maintained via the residual
  `t2 = v − lead*p_lo`; after `t2 += step`, a `while (t2 >= p_lo)` loop
  fires 0 or 1 times because `step ≪ p_lo`. Again no division.

Both are seeded once per t (one short division and at most b1 subtractions).

**(c) Full check only on survivors.** When `bm == lead` (~1/b1 of
candidates), `fullpal()` reverses the low k digits of v and compares with
the high k digits — the standard half-reversal palindrome test, k short
divisions.

**Making divisions cheap where they remain.** All division in the scanner
is by the small constants b1 or b2. The scanner body `scan_fast(B1, B2)` is
`static inline` and instantiated through a macro
(`DEFFAST(3) ... DEFFAST(30)`) as thirty specialized entry points, one per
adjacent pair. With the divisors compile-time constants, the compiler turns
every `/` and `%` — including those inside `udm128` — into
multiply-by-reciprocal sequences. `udm128` itself divides a 128-bit value by
a small constant portably using 64→32-bit limb long division (three 64-bit
constant divisions, no `__udivti3`, no hardware 128-bit divide needed —
relevant on arm64).

### 2.3 The even-n scanner (`scan_even`)

For even n ≥ 6 every stage uses this path (partner bases are never
adjacent, see §1.3, so the fast adjacent scanner does not apply). It uses
the same row decomposition `h = t*b2 + d0`, `v = base(t) + d0*step`, but
adds **divisibility skipping**: a base-b1 palindrome must satisfy
`v ≡ 0 (mod q)`, q = b1+1. Within a row, `v mod q` is linear in d0, so
the valid d0 form an arithmetic progression with common difference
`q / gcd(step, q)`. Per row and partner the scanner computes `base mod q`
(one short division), solves for the AP start with a precomputed modular
inverse, and then visits only every stride-th candidate — a factor
~q/gcd(step,q) (typically 7–23×) reduction before any filtering. Each
visited candidate still goes through the first/last-digit filter and only
then the full check. Rows where `gcd ∤ (−base mod q)` are skipped for that
partner with no visits at all.

A pleasing detail: for the adjacent partner q = b2 divides `step` exactly,
and the row condition degenerates to "t's leading digit ≡ 0" — impossible.
The scanner would discover the §1.3 theorem numerically on its own;
excluding those partners in setup just avoids paying for it.

### 2.4 The meet-in-the-middle scanner (`-M`, even n)

The scanners above all pay ~W^(1/2): they enumerate every palindrome of
the larger base below the bound W. The MITM scanner (mode 3, `-M` flag)
pays ~W^(1/4) by never enumerating whole palindromes at all.

For even n with m = n/2, a base-b2 palindrome is v = y·b2^m + rev(y).
Split the half y = yh·b2^t + yl into its top s and bottom t digits:

```
v = yh·b2^(m+t) + yl·b2^m + rev_t(yl)·b2^s + rev_s(yh)
```

— exactly separable into yh-terms and yl-terms. Now impose
b1-palindromicity: v mod b1^j must equal rev_j(p), where p is v's top j
base-b1 digits. Two observations close the trap:

1. p is determined by yh alone (the yl-terms are smaller than one yh
   step), up to a carry range of a couple of values.
2. Modulo Q = b1^j the palindrome equation splits into
   `B(yl) ≡ rev_j(p) − A(yh) (mod Q)` with B depending only on yl and A
   only on yh.

So: bucket all b2^t values of yl by B once (~W^(1/4) table entries), then
sweep yh in increasing order — increasing v, preserving first-match
semantics and the early-cutoff logic — and for each (yh, carry-candidate
p) look up the matching bucket and full-check the handful of survivors.
Completeness: a doubly palindromic v has a unique (yh, yl) split and its
true prefix p lies in the scanned carry range, so it cannot be missed.

Parameters are chosen per stage: t maximal with b2^t under the table cap
(`-T`, default 3×10^8 entries), j maximal with b1^j ≤ b2^t (buckets
average ~1 entry and the carry range stays ≤ 2–3). The bucket table is
built once per (b1, b2) pair and shared read-only by all threads.

Two refinements matter at depth. **Two-level filtering**: each table
entry also stores the residue at an extended prefix depth j+e (chosen
adaptively), so in-bucket candidates are rejected with one integer
compare before any wide arithmetic. This is what tames base pairs with a
common factor (24|26, 26|28), whose bucket keys otherwise collapse onto a
sublattice — mostly-empty buckets plus a few enormous ones that dominate
wall time. The table build itself is division-free: both residue terms
are maintained incrementally (arithmetic progression + digit odometer).
**w192 values**: beyond ~3.4×10^38 (mid-n=28) values no longer fit u128;
the MITM path and all window bookkeeping run on a 3×64-bit type with
exactly the operations the search needs, while the legacy scanners remain
u128 and refuse stages that would overflow.

Measured effect: n = 22 drops from 15 minutes (fleet) to **0.3 s**; the
n = 24 floor that took ~2 machine-days re-verifies in **7.8 s**; and
a(24) itself — which the fleet would have reached only after weeks —
falls in **50 s** on one machine. The scanner was validated by
reproducing all previously known terms exactly (even and odd) and by
agreeing with the half-enumeration fleet's stage-by-stage "empty"
verdicts for n = 24.

### 2.5 The generic scanner (`scan_gen`)

Used for small n (< 6) and for odd-n stages with several partner bases.
It reconstructs each palindrome from its half directly and tests the value
against *every* partner b1 whose window contains it, using the same
first/last-digit prefilter, but with runtime divisors.

### 2.5 Threading and termination

- Workers pull **blocks** from a shared atomic cursor (`g_next`): 4096
  t-values per grab in the fast path (~4096·b2 candidates), 65536 halves in
  the generic path. Work stealing is automatic; no partitioning decisions.
- The best match lives behind a mutex (`report_match` keeps the minimum and
  announces improvements on stderr). Workers re-read it once per block —
  cheap — and compare against a local `stop = min(best, stage_hi)`.
- Because values are monotone in the enumeration order, any worker seeing
  `v >= stop` sets a shared cutoff flag (`g_cut`); all workers drain out and
  the stage ends. A match found mid-stage therefore truncates both the
  current stage and the main stage loop naturally.
- Correctness of the minimum does not depend on scheduling: a match only
  *shrinks* the ceiling, and every candidate below the final answer is
  provably visited before its scan is allowed to stop.

### 2.6 Progress reporting

A status thread wakes every 2 s and prints to stderr: stage %done
(candidates processed / stage total), smoothed throughput, stage ETA, and —
once an upper bound exists (a match found, or `-L` given) — overall %done
and ETA, computed by summing the clamped candidate counts of all stages
below the bound (`stage_total_bounded`). Until a bound exists the total
workload is open-ended, so only stage-level progress is meaningful.

### 2.7 Result verification

The winning value is re-verified from scratch: the program extracts its
digits in every base from 2 to 128, and prints every base in which it is an
exactly-n-digit palindrome (so triple-base coincidences would be visible).
This check shares no code path with the search's incremental filters.

---

## 3. Validation

The program reproduces all sixteen computable known terms, a(2)–a(17),
including the non-adjacent-pair cases:

```
a(2)=8 (3,7)      a(3)=26 (3,5)     a(4)=624 (5,7)    a(5)=2293 (5,6)
a(6)=207702 (8,10)                  a(7)=186621 (6,7)
a(8)=342324801 (12,16)              a(9)=27924649 (7,8)
a(10)=260311602096 (15,17)          a(11)=1556085529 (7,8)
a(12)=248876637484140 (17,19)       a(13)=318713056300 (8,9)
a(14)=2544221971606336 (13,15)      a(15)=4712469842177 (7,8)
a(16)=530386561769238496 (13,15)    a(17)=1939137135947326 (8,9)
```

## 4. New results

### Complete table of terms (as of 2026-09-28)

| n | a(n) | bases | status |
|---|---|---|---|
| 1 | 1 | all b ≥ 2 | prior OEIS |
| 2 | 8 | 3, 7 | prior OEIS |
| 3 | 26 | 3, 5 | prior OEIS |
| 4 | 624 | 5, 7 | prior OEIS |
| 5 | 2293 | 5, 6 | prior OEIS |
| 6 | 207702 | 8, 10 | prior OEIS |
| 7 | 186621 | 6, 7 | prior OEIS |
| 8 | 342324801 | 12, 16 | prior OEIS |
| 9 | 27924649 | 7, 8 | prior OEIS |
| 10 | 260311602096 | 15, 17 | prior OEIS |
| 11 | 1556085529 | 7, 8 | prior OEIS |
| 12 | 248876637484140 | 17, 19 | prior OEIS |
| 13 | 318713056300 | 8, 9 | prior OEIS |
| 14 | 2544221971606336 | 13, 15 | prior OEIS |
| 15 | 4712469842177 | 7, 8 | prior OEIS |
| 16 | 530386561769238496 | 13, 15 | prior OEIS |
| 17 | 1939137135947326 | 8, 9 | prior OEIS |
| 18 | 73012066376764835227425 | 20, 22 | **new, proven** |
| 19 | 6411682614162861146 | 10, 11 | 2012 candidate, **confirmed** |
| 20 | 243177605553755014828397184 | 21, 23 | **new, proven** |
| 21 | 104618510424015816401 | 9, 10 | **new, proven** |
| 22 | 719621573294884539037143600 | 17, 19 | **new, proven** |
| 23 | 89403957605050675930498 | 10, 11 | **new, proven** |
| 24 | 1205824352934861911999927836207125 | 24, 26 | **new, proven** |
| 25 | 9986831781362631871386899 | 10, 11 | **new, proven** |
| 26 | 373535253009400410279563231242041 | 18, 20 | **new, proven** |
| 27 | 11729080397241840985487045425 | 11, 12 | **new, proven** |
| 28 | 5368922667686573493734112244355570527560 | 27, 29 | **new, proven** |
| 29 | 15794975400931374793240751579457 | 12, 13 | **new, proven** |
| 30 | 3019184969309064736498706779244850699507840 | 27, 29 | **new, proven** |
| 31 | 24683909945699206608160858458686717 | 13, 14 | **new, proven** |
| 32 | 48849507047691135113700496464447321657547821120 | 29, 31 | **new, proven** (sealed 2026-09-28) |
| 33 | 43490338130756441081401879246988772601 | 14, 15 | **new, proven** |
| 34 | 3807947712122126639487215545482127957525027225 | 22, 24 | **new, proven** |
| 35 | 955983290720152593554829802352458678189 | 13, 14 | **new, proven** |
| 36 | *(in search)* | — | stages 23–28 certified empty, so a(36) > 29³⁵ ≈ 1.53×10⁵¹ (§6) |
| 37 | 2405399997091326396338729111112302450399341 | 14, 15 | **new, proven** (publication gated on a(36) — OEIS needs consecutive terms) |

Every "proven" term carries an exhaustive minimality certificate (logged
stage-by-stage) and an independent from-scratch base-conversion
verification. Additionally, every term from a(33) on has been
independently re-derived end-to-end on a second GPU architecture
(CUDA-found terms re-discovered via Metal, 2026-09-29) — two
implementations, two silicon vendors, identical results. Patterns to date: all eleven odd-indexed new terms use
adjacent bases found near their range floor; even-indexed terms use
gap-2 pairs (odd-odd with the gcd = 2 boost, or even-even), with (27,29)
and (29,31) each appearing twice. The sequence is not monotone: a(34) <
a(32), a(26) < a(24), a(9) < a(8).

### Detailed discovery notes (chronological)

### a(18) = 73012066376764835227425 — new term, open since 2009

```
base 20:  5b830ba7227ab038b5
base 22:  1258i9fj77jf9i8521
```

Found by exhaustively scanning stages b2 = 15…22 (197 billion candidates,
88 s on 10 threads of an Apple M-series laptop). Stages 15–21 are empty;
the match appears early in the (20, 22) stage. Independently re-verified by
a from-scratch base-conversion scan over all bases 2–200, which finds
exactly the two representations above.

Why it stayed open for 14+ years: a(18) ≈ 7.3×10^22 is a **76-bit** number.
By the even-n theorem (§1.3) the search cannot even begin below the first
gap-2 pair (13, 15) at ~10^20, and the answer sits seven stages further
out. A scan at this depth needs 128-bit arithmetic and hundreds of billions
of candidates — far outside 2012-era reach, and the value is ~40000× larger
than the a(19) term that had already been found. Amusingly, the winning pair
(20, 22) enjoys no divisibility bonus (gcd(21, 23) = 1).

### a(19) = 6411682614162861146 — 2012 candidate confirmed minimal

```
base 10:  6411682614162861146
base 11:  11759a6746476a95711
```

Merickel's candidate was correct. The exhaustive scan below it (the only
two possible stages, (9, 10) and (10, 11): 712 million candidates, 2 s)
finds nothing smaller and hits exactly this value. With a(18) now known,
both terms are publishable to OEIS.

### a(20) = 243177605553755014828397184 — new term

```
base 21:  i79ce0fje44ejf0ec97i
base 23:  35m14bl46bb64lb41m53
```

An **88-bit** number, ≈ 2.43×10^26 — a thousand times larger than a(18).
An overnight run (pre-theorem binary) cleared stages through b2 = 22; the
match then appeared 82% of the way through the (23 vs {20,21}) stage, and
minimality was sealed by scanning the (24 vs {21,22}) window up to the
match value (the scan self-truncates there — 14% of that stage), 5.26
trillion candidates in ~17 minutes on 14 threads for the resumed portion.
Independently re-verified by a from-scratch base-conversion scan over
bases 2–300.

The winning pair (21, 23) fits the even-n pattern: both bases odd, so
gcd(b1+1, b2+1) = 2 doubles the coincidence probability — the same shape
as the known terms a(10) = (15,17), a(12) = (17,19), a(14)/a(16) = (13,15).

### a(21) = 104618510424015816401 — new term

```
base 9:   854013340414043310458
base 10:  104618510424015816401
```

The easy one of the batch: a 21-digit base-10 palindrome (readable
directly), found 2 seconds into the very first legal stage (10 vs 9) —
classic odd-n behavior, adjacent bases near the bottom of the range, just
like a(13), a(15), a(17), a(19).

### a(22) = 719621573294884539037143600 — new term

```
base 17:  a710g64f9b99b9f46g017a
base 19:  102dia3755885573aid201
```

A 90-bit number, ≈ 7.20×10^26 — only ~3× a(20), far below the cost
model's extrapolation, because the match landed just 1% into the
(19 vs {17}) stage after stages 17 and 18 came up empty (1.46T candidates,
~15 min). Another both-odd, gcd = 2 pair. Independently re-verified over
bases 2–300. (OEIS requires consecutive terms, so a(22) needed a(21) above
— the same rule that held Merickel's a(19) hostage to a(18) for 14 years.)

### a(23) = 89403957605050675930498 — new term

```
base 10:  89403957605050675930498
base 11:  110990a66a626a66a099011
```

Another base-10 palindrome, found 32 s into the first legal stage
(11 vs {10}) — the same base pair as a(19), and the fourth consecutive
odd-n term on adjacent bases. Independently re-verified over bases 2–300.

### a(24) = 1205824352934861911999927836207125 — new term, via MITM

```
base 24:  lgi53dd79enccne97dd35igl
base 26:  3bd6fe8i3cbbbbc3i8ef6db3
```

The term that motivated the meet-in-the-middle scanner (§2.4) — and
vindicated it: a(24) ≈ 1.21×10^33 (110 bits) sits in stage **26**.
Stages 23, 24 and 25 are all empty; the multi-machine half-enumeration
campaign, then grinding stage 23, was weeks away from it. The MITM run
found and proved it in 50 seconds on one laptop, itself re-verifying
stages 18–22 in agreement with the fleet. First even term on an
even-even base pair (gcd(25, 27) = 1 — no divisibility bonus).

### a(25) = 9986831781362631871386899 — new term

Bases 10 and 11 — the *third* appearance of that pair (a(19), a(23),
a(25)), readable directly as a base-10 palindrome. First stage, 26 s.

### a(26) = 373535253009400410279563231242041 — new term, via MITM

```
base 18:  f922hdh5e6a3443a6e5hdh229f
base 20:  1255f9f857g0770g758f9f5521
```

Found in **4 seconds**, second stage. Notably a(26) < a(24) — the
sequence is not monotonic (cf. a(9) < a(8)).

### a(27) = 11729080397241840985487045425 — new term, via odd-n MITM

Bases 11 and 12, first legal stage, found in **2 seconds** once the MITM
construction was extended to odd n (§2.4: shift exponent k instead of m,
low-side reverse over t−1 digits of yl div b2 — the unmirrored middle
digit is the only change). Validated by reproducing all nine known odd
terms. The brute-force scanner had spent an hour on the same stage and
was one-third through.

### a(28) = 5368922667686573493734112244355570527560 — new term, beyond u128

```
base 27:  c2ofqc7b327noaaon723b7cqfo2c
base 29:  1lspr8gbfakrpiiprkafbg8rpsl1
```

≈ 5.37×10^39 — a **132-bit** value, the first term that cannot be
represented in 128-bit arithmetic; finding it required the w192 wide-value
port. Stages 26–28 are empty (stage 26 certified by two independent runs);
the match sits 26% into the two-partner stage 29, where partner 26's whole
window (topping out below the match) was exhausted first and partner 27's
monotone sweep produced it. Bases 27 and 29: back to the classic odd-odd,
gcd = 2 pattern after a(24)/a(26)'s even-even interlude. Total search ~46
minutes on one Mac Studio; the pre-MITM cost model would have priced this
at ~10^20 brute-force candidates — thousands of years of fleet time.

### a(29) = 15794975400931374793240751579457 — new term

Bases 12 and 13, first legal stage, 2 seconds. The odd-n pattern is now
eight for eight: adjacent bases, found near the bottom of the range.

### a(30) = 3019184969309064736498706779244850699507840 — new term

```
base 27:  994kek89p60nm0bb0mn06p98kek499
base 29:  152pfcfjnfsk6dppd6ksfnjfcfp251
```

≈ 3.02×10^42 (142 bits) — and **the same base pair (27, 29) as a(28)**,
one digit-length deeper. Stages through 28 empty; the match sits early in
partner 27's sweep of stage 29, after partner 26's entire window (topping
out below the match at 26^30) was exhausted first. ~2 h 52 m on one Mac
Studio, 1.96T MITM steps.

### a(31) = 24683909945699206608160858458686717 — new term

Bases 13 and 14, first legal stage, 6 seconds. Odd-n pattern: nine for
nine (adjacent bases, bottom of the range).

With a(18)–a(31) all resolved, the sequence extends from 17 to 31 terms.
All new terms carry exhaustive minimality proofs and independent
from-scratch verification.

### The even-n adjacent-base theorem

Found while diagnosing why the n = 20 search was slow (§1.3): **for even n,
no number is an n-digit palindrome in two adjacent bases** — even-length
palindromes in base b1 are divisible by b1+1, contradicting a nonzero
trailing digit in base b1+1. It retroactively explains the shape of the
known data (every even-n term uses bases with gap ≥ 2, every odd-n term so
far uses adjacent bases), prunes provably-empty stages before any scanning,
and combined with AP-skipping (§2.3) sped the even-n scanner up by roughly
two orders of magnitude.

## 5. The tool family: three implementations, one algorithm

All three run the same MITM algorithm with identical parameter selection,
stage semantics, banners and certificate lines, and were each validated by
reproducing every known term (a(8)–a(31)) end to end. Any match from any
of them is re-verified in exact host arithmetic before being believed.

### 5.1 `a171740` — CPU (portable C, macOS/Linux)

`make a171740` (needs only cc + pthreads). The reference implementation:
all scanners (§2), all flags. MITM mode is `-M`. Work-stealing blocks mean
its status counters are *not* contiguous, so `-H` resume subtracts an
in-flight safety margin automatically. Measured ~45–300 M/s depending on
machine and stage depth (cache-bound on multi-GB tables).

### 5.2 `a171740g` — Apple-GPU (Metal, macOS)

`make a171740g` (builds `a171740g_host.m` + `kernel_src.h`; the kernel is
compiled at runtime — no Xcode Metal toolchain required). Table built by
CPU threads, shared with the GPU via unified memory; each GPU thread
sweeps a sequential yh-chunk, amortizing one exact 192-bit division seed.
Porting essentials: Apple GPUs have no integer divide, so every hot
division is an exact multiply-by-reciprocal (q = mulhi(x, ⌊2^64/d⌋) + one
correction, provably exact for all 64-bit x); no doubles exist and none
are needed. Dispatch slabs complete in ascending yh order, so status
counters are contiguous frontiers and `-H` resume is *exact* (no margin).
Flags: `-n -B -E -L -H -P -T -S -q` (`-S` = yh per dispatch slab).
Measured: **~2.3 G/s on M4 Max (32 cores), ~1.9–2.1 G/s on M2 Ultra** on
deep n=32 stages — one laptop GPU outruns the entire CPU fleet.
(Development history: `mitm_gpu.m` is the CPU-vs-GPU benchmark prototype
that validated the kernel; a(24)/a(28) rediscovered with byte-identical
filter statistics, 0.4×–17× vs 10 CPU cores depending on stage shape.)

### 5.3 `a171740c` — NVIDIA GPU (CUDA, Linux)

`nvcc -O3 -arch=native -o a171740c a171740c.cu -lpthread`. Direct port of
the Metal worker: same host logic, kernel translated to native CUDA
(`__umul64hi`, `atomicAdd`), tables in `cudaMallocManaged` (coherent
unified memory on Grace-Blackwell). Same flags and semantics as
`a171740g`, including exact `-H` resume. Compiled and validated first-try
on DGX Spark (GB10).

### 5.4 Measured performance by platform

Deep-stage MITM sweep rates (yh-steps/s on n=32–36 production stages;
same exact algorithm on every platform — rates vary with stage shape, so
ranges are shown; "dense band" = degenerate shared-factor bucket regions):

| machine | CPU | GPU | RAM | CPU rate | GPU rate | notes |
|---|---|---|---|---|---|---|
| MacBook M4 Max | 10P+4E | 32-core Apple | 36 GB | 33–117 M/s | **0.7–2.33 G/s** | GPU 4.7–17× its own CPU; 0.7–0.8 G/s in dense bands |
| Mac Studio M2 Ultra | 24-core | 60-core Apple | 64 GB | 155–297 M/s | **0.79–2.08 G/s** | fastest CPU; GPU rate dips in dense bands |
| DGX Spark ×2 (GB10) | 20-core Grace | Blackwell, 48 SMs | 121 GB | — (unused) | **1.2–2.8 G/s** | only boxes with RAM for t=7 tables (up to ~101 GB); 0.15–0.34 G/s in worst dense band |
| house1 (Apple, 10c) | 10-core | — | 64 GB | 45–169 M/s | — | rate strongly stage-dependent |
| mathf (Apple, 10c) | 10-core | — | ≥32 GB | 80–127 M/s | — | retired mid-campaign |
| 10.1.10.20 (Apple, 10c) | 10-core | — | 32 GB | 109–124 M/s | — | retired mid-campaign |
| 10.1.1.74 (Apple, 10c) | 10-core | — | 64 GB | 45–49 M/s | — | slowest fleet CPU |
| Jetson Orin Nano | 6× A78AE | 1024-core Ampere | 7.3 GB | 18.8 M/s | **~360 M/s** | RAM caps tables at t=5 (≈28× span inflation), so effective deep-stage throughput ÷28; great for odd terms / small slices |
| dual Xeon Gold 5315Y | 2×8-core Ice Lake | — | **1.5 TB** | 48–94 M/s | — | first x86 port (exact results); holds 133 GB+ t=7 tables no GPU box can; designated t=8 experiment platform |

Same-task wall-clock anchor — a(24) full search (`-M`, exhaustive through
its match): M4 Max CPU 57 s · Orin Nano GPU 1 m 57 s · Orin Nano CPU
4 m 46 s. Deep anchor — a(34) full campaign: Spark t=6 53 min, Spark
t=7 9 m 43 s, vs 2 h 52 m on the M2 Ultra-era CPU fleet.

Two structural lessons the table encodes: GPU latency-hiding beats CPU
cache hierarchies decisively once tables exceed L2 (every GPU outruns
every CPU, even a Jetson against a 24-core Mac Studio), and **table
memory is the real currency** — the Sparks' 121 GB (t=7 depth) is worth
more than raw FLOPS, while the Nano's 7 GB demotes an otherwise-fast GPU
to light duty. (Historical rates from the pre-MITM half-enumeration
scanners — e.g. "8.8 G candidates/s" on M4 Max CPU — measure a different
unit and are not comparable to these figures.)

### Fleet conventions (all versions)

Stage windows are disjoint, so machines combine by simple minimum over
matches once every lower stage/slice is certified. `-B`/`-E` fence a run
to owned stages; `-L` caps by value (slice boundaries); `-H`/`-P` resume
or start mid-stage. Keep every run's log: the stage banners and final
`done:`/`no ... found below` lines are the exhaustiveness certificates.
Campaign state lives in a per-target ledger (see `a32-ledger.md`).

## 6. Current search status: a(36) endgame

As of 2026-10-07 (a(32) was found, sealed and recorded — see §4; this
section tracks the one remaining gap):

- **a(36) is the only unresolved term below a(38)**: stages 23–28 of its
  campaign are **all certified empty**, putting **a(36) > 29³⁵ ≈
  1.53×10⁵¹** (~170 bits), the most expensive single term in the project's
  history (~650T MITM steps certified so far).
- **Stage 27 ((25,27) pair, 142T steps): certified empty** 2026-09-30
  05:28. It ran as a two-Spark pincer at t=7 with capped-offset tables
  (`-J 300000000`), ~650 M/s per Spark:
  - **Lower half** (atom1, `a36-s27lo.log`): no double palindrome below
    168497487175726784632929073454438737973189159843808.
  - **Upper half** (atom2, `a36-s27hi.log`): resumed exactly after atom2
    was lent out 20:35–22:35. At 03:46 its remaining 7.54×10¹² yh were split
    evenly:
    - atom2 swept to the fence
      209464258471269737150356012074732567339671682600348, exactly
      3,771,335,849,077 yh;
    - atom1 (`a36-s27top.log`) swept from counter 138471346183515 to the
      end of the stage (bound 28³⁵), exactly 3,771,335,849,076 yh.
  - Every segment boundary is an exact CUDA resume counter; no match
    anywhere.
- **Stage 28 ((26,28) pair, 279T steps): certified empty** 2026-10-07. It
  ran as a two-Spark pincer at t=7. The table was 28⁷ = 1.35×10¹⁰ entries
  (~103 GiB of a Spark's 121), with the bucket cap raised to
  `-J 310000000` so that Q = 26⁶ = 308,915,776. Speed varied by region from
  ~0.7 to 2.6 G/s per Spark, about 2–4× stage 27.
  - **Lower half** (atom1, `a36-s28lo.log`): no double palindrome below
    658135847746348224906843967211475118986813180477440. The final segment
    swept exactly 124,196,078,721,966 yh.
  - **Upper half** (atom2, `a36-s28hi.log`): no double palindrome below
    29³⁵. The final segment swept exactly 121,257,247,349,677 yh.
  - Both halves were paused 2026-09-30 08:05 and resumed 2026-10-06 17:15
    at their exact counters; no match anywhere.
- **Stage 29 ((27,29) pair, 511T steps per pass) in flight since
  2026-10-07 ~16:30**, using CUDA worker v6 (`a171740c_v6.cu`). Its t=7
  table, 29⁷ = 1.72×10¹⁰ entries (138 GB), does not fit a Spark, so v6
  adds two changes:
  - **Table split** (`-Y k/K`): build only part k of K of the yl table. K
    runs over the same yh range cover the stage. Here atom1 takes part 0/2
    (`a36-s29p0.log`) and atom2 part 1/2 (`a36-s29p1.log`), each ~72 GB and
    each sweeping all 510,712,773,402,664 yh.
  - **Variable residue width**: entries pack yl above a residue of
    bh = 64 − bits(b2^t − 1) bits (at most 30). v5's fixed 30 bits would
    have silently overflowed here, because 29⁷ > 2³⁴.

  Validation:
  - a(24) and a(26) are re-found exactly through the unchanged path.
  - n = 28, stage 29 at t=7 with the table split 0/2 + 1/2 (the 29-bit
    residue path): part 0 re-found a(28) =
    5368922667686573493734112244355570527560, and part 1 swept the whole
    stage (1,953,650,931 yh, both partners) with no match.

  A match found by one part is minimal once the other part has swept past
  it. Each part's certificate carries a "TABLE PART k/K ONLY" note; both
  are needed to certify the stage.
- **Outlook** (model of expected solutions per stage, calibrated against the
  known even terms n = 18–34, which it over-predicts by ~1.4×), given stages
  ≤ 28 are empty:
  - 19.5% chance a(36) is in stage 29;
  - 29% cumulative through stage 30, 46% through stage 31, 54% through
    stage 32;
  - on two Sparks, stages 29–31 take roughly 2 weeks to 2 months in total.
    Each needs multi-pass tables, because its t=7 table outgrows a Spark:
    stage 29's 138 GB table exceeds a Spark's 121 GB, so it needs splitting
    across two passes or denser entry packing;
  - past stage 33 (two partner bases per stage) the current method needs
    years.

  The Jetson Orin Nano holds
  a(39) stage 16 (odd, banked-term duty) but runs it at only ~14 M/s
  at t=7, so the stage goes back to a faster GPU when one frees up.
- **The stage-26 fight produced two permanent tool upgrades**: the
  *mega-bucket fix* (degenerate shared-factor pairs create buckets of
  millions of entries; sorted buckets + binary search turned observed
  1000× slowdowns into full-speed sweeps) and the *bucket-count cap*
  (`-J`), which shrinks offset arrays so the 121 GB Sparks can hold
  t=7 tables for big bases — turning stage 27 into a ~29-hour job on
  two Sparks instead of a multi-week one at t=6. Both fixes are validated in the CUDA
  worker (v4/v5) and the portable CPU tool, on known terms, at both
  table depths.
- **Big-memory x86 validated**: the 1.5 TB dual-Xeon reproduced a(34)
  exactly via a full t=7 campaign (54 m 32 s) and held stage 27's
  133 GB tables; its 48 M/s sweep rate benches it for now, but it is
  the designated t=8 experiment platform.
- Stage 27 held no term, so stage 28 (~2× larger) is running under the
  same pincer pattern. A match there completes the consecutive run
  a(1)–a(37). If stage 28 is also empty, stage 29's t=7 table no longer fits
  a Spark (see above). After a(36): full stop and records consolidation, per
  plan.

## 7. Cost model / frontier

Three algorithm generations, three cost regimes:

1. **Half-enumeration** (`scan_fast`/`scan_even`): ~W^(1/2) candidates
   below a bound W, ~1/(b1+1) of them visited for even n via divisibility
   skipping. Anchors: a(18) 197G candidates / 88 s; a(20) ~10.7T /
   overnight+17 min across a three-machine fleet. Practical wall:
   mid-a(20)s.
2. **Meet-in-the-middle** (`-M`, §2.4): ~W^(1/4). Anchors: a(22) 0.3 s
   (fleet: 15 min); a(24) 50 s (fleet projection: weeks); a(27) 2 s
   (brute scanner: 22+ h). Wall: the u128 value ceiling at ~3.4×10^38,
   i.e. mid-n=28.
3. **MITM on w192** (3×64-bit values): same W^(1/4) with ~1.5–2× constant
   overhead. Anchor: a(28) ≈ 5.4×10^39 (132 bits) in ~46 min on one Mac
   Studio — ~10^20 brute-force-equivalent candidates. Value range tops
   out near 6.3×10^57, enough for n ≈ 38–40.
4. **MITM on GPUs** (§5.2/§5.3): same algorithm, same exactness, ~10–50×
   the per-machine throughput on deep stages (2+ G/s per device). Anchor:
   a(30)'s 2h52m CPU discovery re-run in 24m39s on one laptop GPU; n=32
   stages that were CPU fleet-days each now run in ~2–4 GPU-hours. The
   wall moves from arithmetic width (fixed by w192) and per-machine
   compute (fixed by GPUs) to table memory and stage growth (~2× per
   stage) somewhere in the mid-to-high n=30s.

Odd-n terms remain near-free at any depth (adjacent-base stages at the
bottom of the range: a(29) took 2 s). Even-n terms set the pace: each one
starts several orders of magnitude deeper (§1.3) and costs roughly
W^(1/4) ~ b^(n/4) at its stage bases, so minutes at n=28 become hours in
the low 30s and the practical single-machine horizon — absent another
algorithmic idea — is somewhere in the mid 30s, now set by compute rather
than arithmetic width.

## 8. Follow-on work: A171740's cross-referenced sequences

A171740's 23 cross-references (all by James G. Merickel, 2009–2012) are
being extended in [`oeis/`](oeis/README.md), with new tools built on this
project's engine: `gappal` (this search restricted to one base gap) and
`multipal`, a CPU and GPU sieve that counts, for every number in a range,
how many bases make it a k-digit palindrome. Full details, validation and
the run log are in [`oeis/README.md`](oeis/README.md); submission drafts are
in [`oeis/SUBMISSIONS.md`](oeis/SUBMISSIONS.md). Headline results
(2026-09-29):

**A171704: first number that is a 5-digit palindrome in at least n bases.**
Stuck at four terms since 2009 (only "a(5) > 10¹⁴" was known):

| n | a(n) | Bases | Status |
|---|---|---|---|
| 5 | 4922057407205376 = 8376⁴ | 2093, 2791, 4187, 5906, 8375 | new, proven |
| 6 | 34777153514704896 = 13656⁴ | 3218, 3413, 4551, 4828, 6827, 13655 | new, proven |

8376⁴ is 1 4 6 4 1 in base 8375, and d⁴ × (1 4 6 4 1) in base 8376/d − 1 for
d = 2, 3, 4; its fifth base, 5906, is a coincidence. The exhaustive search
(every 5-digit palindrome in every base below each value, about 5×10¹³ in
total) is what finds such coincidences: the pure construction would have
bounded a(5) only at 1.24×10¹⁷.

**Other results**
- **A171703** (4-digit palindromes in ≥ n bases): from 9 proven terms to 22.
  All six of Yamanouchi's 2014 upper bounds, a(10)–a(15), are exact (CPU and
  GPU sieves agree below 8.4×10¹⁵). Five new terms:
  - a(16) = 327600³
  - a(17) = 342720³: 16 construction bases plus one coincidental base,
    which puts it three times below the construction's own bound
  - a(18) = 589680³
  - a(19) = a(20) = 720720³ = 374368864117248000, with 20 bases
  - a(21) = 1053360³ = 1168773800173056000, with 21 bases
  - a(22) = 1441440³ = 2994950912937984000, with 22 bases

  The search continues toward a(23).
- **A171702** (3-digit palindromes in ≥ n bases): from 100 terms to 223,
  exhaustive to 10¹².
- **Two OEIS data errors**: A171701's a(41) is 20160 (not 30160), and
  A171741's a(3) is 121 (not 154).
- **Base-gap family (A216xxx)**: complete. There are b-files of 199–1,000
  terms for lengths 8–15 and 60 terms for length 17, where OEIS had 3–19
  terms. All previous values are confirmed. **The new length-16
  sequence, the missing member of the family, is complete for n = 2..60.**
- **New upper bounds** for A171705(4), A171706(4) and A171741(8), whose
  entries give only lower bounds.
