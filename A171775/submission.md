# OEIS submission draft — A171775 (smallest number that is a k-digit palindrome for every k = 2..n)

Draft of 2026-10-03. Everything below was computed or checked with the tools in this directory
(`a171775.c`, `a171775_2d.c`, `cuda/a171775_cuda.cu`, `cuda/a171775_2d_cuda.cu`, `metal/`,
`verify_a171775.py`); see README.md for the methods and the run logs.

## Status on the OEIS (checked 2026-10-03 18:00 PDT)

* **A171775 — to do.** Unchanged since #64: a(9) by Max Alekseyev (#62–#63, Jun 05 2026),
  approved by Sean A. Irvine on Jun 06 2026. Data a(1..9), keywords `nonn,base,hard,more`, FORMULA
  "Conjecture: a(n) = 2^((n-1)*(n-2)) for n >= 7". No pending edit, no draft. Section 1 is the edit
  to submit.
* The cross-referenced entries (A007632, A029961–A029970, A171701–A171706, A171740–A171742,
  A253294, ...) are not affected by this work.
* The other `src/math` entries were checked the same day. All 16 earlier edits are approved and
  live; the appendix lists them and, per project, what is still open.

## Summary

| item | status |
|---|---|
| a(10) = 2^72 = 4722366482869645213696 | to submit (**new**; four complete searches, two methods) |
| a(11) = 2^90 = 1237940039285380274899124224 | to submit (**new**; three complete searches, plus a cross-method check) |
| comment with the survivor counts below 2^72 and 2^90 | to submit (recommended) |
| EXAMPLE lines for a(10), a(11) | optional |
| PROG (Python, reproduces a(1..6)) | optional |
| LINKS: program link to the public repository | recommended (the directory is already public) |
| Merickel's "equals it for n = 7 and 8 at least" | optional: change to "for 7 <= n <= 11" |
| a(12) | not proposed: only the trivial a(12) > 2^90 is known (see section 3) |

Next step: submit section 1 (DATA, COMMENTS, EXTENSIONS and the program link; the rest is optional).

---

## 1. A171775 — the edit to submit

Current entry: #64, Jun 06 2026, data a(1..9), keywords `nonn,base,hard,more`.

### DATA (append two terms)

```
1, 3, 5, 52, 130, 1885, 1073741824, 4398046511104, 72057594037927936, 4722366482869645213696, 1237940039285380274899124224
```

### COMMENTS (add)

```
There are 135382 numbers M <= 2^72 that are both a 10-digit and a 9-digit palindrome (in some bases); 2^72 is the only one of them that is also an 8-digit palindrome in some base. There are 524157 numbers M <= 2^90 that are both an 11-digit and a 10-digit palindrome; 2^90 is the only one of them that is also a 9-digit palindrome. Hence a(10) = 2^72 and a(11) = 2^90, as conjectured in the Formula section. - _Jeff Sponaugle_, Sep 30 2026
```

Optional, instead of a new line: Merickel's first comment says "(and equals it for n = 7 and 8 at
least)"; it could become "(and equals it for 7 <= n <= 11)". Or the conjecture in FORMULA could get
"(true for n <= 11)". The editors may prefer one of these; the comment above stands on its own.

### EXAMPLE (add, optional)

```
a(10) = 2^72 = 4722366482869645213696: it is 11 in base 2^72-1; in base 18669329 it is the 3-digit palindrome 13548844 15833512 13548844; in bases 2^19-1, 2^15-1, 2^13-1, 2^11-1, 2^10-1, 2^9-1, 2^8-1 its digits are 2^15*C(3,j), 2^12*C(4,j), 2^7*C(5,j), 2^6*C(6,j), 2^2*C(7,j), C(8,j), C(9,j), palindromes of 4 to 10 digits. These are the smallest bases for each length.
a(11) = 2^90 = 1237940039285380274899124224: it is 11 in base 2^90-1, and in bases 2^31-1, 2^23-1, 2^19-1, 2^16-1, 2^14-1, 2^12-1, 2^11-1, 2^10-1, 2^9-1 its digits are 2^28*C(2,j), 2^21*C(3,j), 2^14*C(4,j), 2^10*C(5,j), 2^6*C(6,j), 2^6*C(7,j), 2^2*C(8,j), C(9,j), C(10,j), palindromes of 3 to 11 digits. These are the smallest bases for each length.
```

(Checked digit by digit with `verify_a171775.py check 2^72 10` and `check 2^90 11`; output in README.)

### PROG (add, optional)

Plain brute force; it reproduces a(1..6) in a fraction of a second (tested with Python 3.14):

```
(Python)
def is_kpal(M, k): # is M a k-digit palindrome in some base?
    b = max(2, round(M**(1/k)) - 1)
    while b**(k-1) <= M:
        if M < b**k:
            d, x = [], M
            while x:
                d.append(x % b)
                x //= b
            if d == d[::-1]:
                return True
        b += 1
    return False
def a(n):
    M = 1
    while not all(is_kpal(M, k) for k in range(2, n+1)):
        M += 1
    return M
print([a(n) for n in range(1, 7)]) # _Jeff Sponaugle_, Oct 03 2026
```

### EXTENSIONS (add)

```
a(10)-a(11) from _Jeff Sponaugle_, Sep 30 2026
```

Keyword `more` stays: a(12) is open (conjecturally 2^110).

### LINKS (add, recommended)

```
%H A171775 Jeff Sponaugle, <a href="https://github.com/jeffsponaugle/Mathpub-/tree/main/A171775">C, CUDA and Metal programs, run logs and notes</a>
```

The directory is already public: commit e6f20f4 has the programs, `verify_a171775.py`, the b-file
and the a(10)/a(11) logs. Push the current README and this file before submitting, so the link
does not show the old `OEIS_notes.md`. Same link format as in `../A158939/submission.md`.

### Notes to the editors (paste into the discussion box)

Short version:

```
a(10) = 2^72 and a(11) = 2^90 are new and agree with the conjecture in the Formula section. Both come from exhaustive searches below these bounds: every number that is both an n-digit and an (n-1)-digit palindrome was found (135382 for n = 10, 524157 for n = 11) and then tested for the shorter lengths; only 2^72 and 2^90 pass. Each search was run several times, with different programs, methods and hardware, and gave identical counts and checksums.
```

Longer version, if the editors ask how the terms were computed:

```
Every candidate is an n-digit palindrome in some base B and an (n-1)-digit palindrome in some base c > B. The programs enumerate the palindromes of one of the two lengths and impose the other length by solving a linear congruence for the innermost free digits (the k lowest base-b digits of M must mirror its k highest), with all state updated by additions; matches are then checked digit by digit, and survivors are tested for every shorter length by scanning all admissible bases.
a(10): all 3.7*10^12 ten-digit palindromes <= 2^72 (bases 2..255) were enumerated, with the 9-digit condition solved for one digit: 5.1*10^12 steps; a C program on an Apple M1 Pro (29 min), a CUDA program on an NVIDIA GB10 (154 s) and the C program on a DGX Spark CPU gave identical statistics, and a second method that enumerates the 9-digit palindromes instead (and solves two digits with a lookup table) gave the same 135382 survivors and checksum (12 min).
a(11): the second method enumerated all 7.0*10^15 ten-digit palindromes <= 2^90 (bases 3..1024), with 5.4*10^12 table lookups; three complete runs (M1 Pro CPU 2.8 h, DGX Spark CPU 1.2 h, GB10 GPU with a separately written CUDA kernel 57 min) gave identical statistics (524157 survivors, 1.2*10^12 congruence matches, checksum 0a5d53d568af3f34). The first method gives the same survivors on the window [600^9, 600^9 + 3*600^8].
Checks: the programs reproduce a(1)..a(9); the congruence filters agree with naive enumeration on hundreds of random blocks for n = 7..12 and with a debug build that recomputes the state at every step; an independent Python brute force agrees on small ranges.
```

---

## 2. Verification chain (for reference)

| check | result |
|---|---|
| a(10), 1D method, M1 Pro CPU (`scan_n10.*`) | 1750 s, 5126601984459 steps, 135382 survivors, 1 solution (2^72) |
| a(10), 1D method, GB10 GPU (CUDA) | 154 s, identical, checksum e9eed6127c51a50f |
| a(10), 1D method, DGX Spark CPU | 653 s, identical, checksum e9eed6127c51a50f |
| a(10), 2D method, M1 Pro CPU | 715 s, 135382 survivors, checksum e9eed6127c51a50f |
| a(11), 2D method, M1 Pro CPU (`scan_n11.*`) | 10023 s, 5386659681911 lookups, 524157 survivors, checksum 0a5d53d568af3f34, 1 solution (2^90) |
| a(11), 2D method, DGX Spark CPU (gcc) | 4474 s, identical in every statistic |
| a(11), 2D method, GB10 GPU (CUDA port) | 3411 s, identical in every statistic |
| n = 11 window [600^9, 600^9 + 3·600^8], 1D method | 211 survivors, checksum 1dc0b33a97158c66 = 2D method |
| a(1..9) | reproduced by both methods (a(9): exhaustive in 32 s / 1.5 s) |
| 2^72, 2^90 | digits checked from the definition in Python for every length |

In every run, all survivors except the solution fail at the next length (8 for n = 10, 9 for n = 11).

## 3. Not proposed

* **a(12).** 2^90 is not a 12-digit palindrome in any base (`a171775 check 2^90 12`), so
  a(12) > 2^90, but that follows trivially from a(12) >= a(11). The conjectured value is 2^110. A search
  below 2^110 with the 2D method would need about 1.4·10¹⁶ lookups: roughly 100 GPU-days on one GB10
  at the rate of the a(11) run. It has not been attempted.
* **A b-file.** Not needed; the 11 terms fit in DATA.

---

## Appendix — OEIS status of the other `src/math` submissions (checked 2026-10-03 18:00 PDT)

Approved and live (credited to Jeff Sponaugle):

| entry | revision | content |
|---|---|---|
| A137723 | #17–#18, Sep 17 2026 | a(32) onward, b-file n = 1..193, comment a(194) > 2^64-194 |
| A252768 | #21, Sep 17 2026 | a(8)-a(9) and the depth-count comment |
| A053686 | #28, Sep 17 2026 | a(13) = 906 and the repeat comment |
| A133788 | #13, Sep 17 2026 | a(13) |
| A085237 | #44, Sep 17 2026 | a(80)-a(86), b-file n = 1..86 |
| A089180 | #28, Sep 17 2026 | example corrected ("corrected by Jeff Sponaugle, Sep 17 2026") |
| A248701 | #47, Sep 18 2026 | a(8)-a(10) |
| A248704 | #19, Sep 17 2026 | a(7)-a(9) |
| A248702 | #25, Sep 20 2026 | a(7)-a(11) |
| A248703 | #25, Sep 20 2026 | a(7)-a(9) |
| A359636 | #31, Sep 18 2026 | a(9) and a(10) > 31610555571634177 |
| A158939 | #38, Sep 22 2026 | a(16)-a(17) and a(18) > 5.30*10^16 |
| A229832 | #35, Sep 22 2026 | a(15)-a(16) |
| A133697 | #27, Sep 22 2026 | a(14)-a(15) |
| A306256 | #18, Oct 02 2026 | a(4) and the search-bound comment |
| A063501 | #29, Oct 03 2026 | a(18)-a(20) |

Still open. On 2026-10-03 each project got its own `submission.md`, and those files are
authoritative; the column below is only a pointer:

| project file | open items |
|---|---|
| `submission.md` (this file) | A171775 a(10)–a(11) |
| `../A350826/submission.md` | A350826 a(18)–a(20), a verified-range comment, the dead `PI_06.html` link |
| `../A248701/submission.md` | A248704(10) = 1958854030679863 (found Sep 19, never submitted); obsolete comment in A248702, `%E` indices in A248703, bounds |
| `../A001208/submission.md` | A053348(9) = A005344(8) = 6082, A001211(27) = 186942 plus a correction, literature terms (A001210, A053346), b-files (A001209, A084192, A084193, A196416), caveat comments |
| `../A034822/SUBMISSION.md` | b-file extensions of A263618, A002778, A002779, A016113, A027829 (De Geest's table); A034822 comment |
| `../A158939/submission.md` | program link, EXAMPLE for a(17); optional items; a(18) bound should wait |
| `../A359636/submission.md` | a(10) lower bound; optional EXAMPLE, comments, program link |
| `../A252768/submission.md` | a(10) lower bound, repair of the Sep 16 comment, program link; optional extras |
| `../A053686/submission.md` | comment listing the second occurrences, program link; optional extras |
| `../A137723/submission.md` | follow-up edit (no computation needed); new terms need runs |
| `../A001220/submission.md` | optional now: drop the old A306256 comment; bound updates when the search resumes (A039951 is blocked by R. Fischer's pending edit) |
| (no draft yet) | A055206(13) = 6470105925 (balanced primes below 10^13), from the A335406 scans (`../A335406/README.md`); A055206 still at #13 |

No submission: A335406 (a(6) is out of reach), A002966, A005234 and A007540 (feasibility
studies, verdict: not feasible), A112389 (in progress); A005115, A007508, A011541 and A101232
have no recorded results.
