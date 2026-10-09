# OEIS A007357 — proposed update: a(22) and a(23)

*Prepared 2026-10-01 against A007357 revision 64 (edited 2026-08-31).*

## Summary

Two new terms can be added. Both were already known infinitary perfect
numbers, found by Moews & Moews (1998); what is new is the proof that
nothing smaller was missed, so they are the 22nd and 23rd terms.

| term | value | factorization |
|---|---|---|
| a(22) | 6149822195556228360752824320000 | 2^18·3^5·5^4·7^2·11·41·79·83·157·313·331·65537 |
| a(23) | 9224733293334342541129236480000 | 2^17·3^6·5^4·7^2·11·41·79·83·157·313·331·65537 |

Proposed edits:

1. Add both terms (DATA and/or b-file — see the DATA length note).
2. Replace the "below 10^28" comment.
3. Add an EXTENSIONS line.
4. Replace the factorization file `a007357.txt` with one covering
   a(18)–a(23), and update its link title.
5. Optionally, upload a b-file and the search program.

## What A007357 says now (revision 64)

- DATA: 21 terms, ending at a(21) = 21623407345626345971712000.
- COMMENT: "No other infinitary perfect number below 10^28. - _Jeff
  Sponaugle_, Aug 30 2026"
- LINK: "Jeff Sponaugle, Factorization for a(18)-a(21)" → `/A007357/a007357.txt`
- EXTENSIONS: "a(18)-a(21) from _Jeff Sponaugle_, Aug 30 2026"
- b-file: none uploaded; OEIS shows one synthesized from DATA (21 terms).
- Keyword: `nonn`.

## Proposed field changes (ready to paste)

### DATA

Append:

```
,6149822195556228360752824320000,9224733293334342541129236480000
```

**Length note:** OEIS asks for DATA of about 260 characters at most. The
current DATA is already 282 characters; with these two terms it would be
346. Editors may prefer to leave DATA as it is and carry a(22)–a(23)
only in a b-file, so submit the b-file below either way and let the
editor decide.

### b-file (`b007357.txt`)

```
1 6
2 60
3 90
4 36720
5 12646368
6 22276800
7 126463680
8 4201148160
9 28770487200
10 287704872000
11 1446875426304
12 2548696550400
13 14468754263040
14 590325173932032
15 3291641594841600
16 8854877608980480
17 32916415948416000
18 1342989770695372800
19 20144846560430592000
20 1441560489708423064780800
21 21623407345626345971712000
22 6149822195556228360752824320000
23 9224733293334342541129236480000
```

### COMMENT

Replace the existing "No other infinitary perfect number below 10^28"
comment with:

```
a(18)-a(23) were known from earlier non-exhaustive searches (Cohen 1990; Moews & Moews 1998); an exhaustive search over Fermi-Dirac factorizations confirms there are no other infinitary perfect numbers up to a(23). The smallest known infinitary perfect number larger than a(23) is 92247332933343425411292364800000 = 10*a(23) (Moews & Moews 1998), but it is not known whether any lie in between. - _Jeff Sponaugle_, Oct 01 2026
```

If editors want a shorter version: "There are no other terms below
a(23). - _Jeff Sponaugle_, Oct 01 2026" is implied by DATA alone, so
the first sentence and the "next known" sentence carry the useful
information.

### EXTENSIONS

```
a(22)-a(23) from _Jeff Sponaugle_, Oct 01 2026
```

### LINK

Change the title to:

```
Jeff Sponaugle, <a href="/A007357/a007357.txt">Factorization for a(18)-a(23)</a>
```

## Replacement factorization file (`a007357.txt`)

The current file stops mid-sentence ("...(a(19), a(20))") and only
covers a(18)–a(21). Proposed replacement:

```
a(18) 1342989770695372800 2^12·3^5·5^2·7^2·11·13·17·41·43·257
a(19) 20144846560430592000 2^12·3^6·5^3·7^2·11·13·17·41·43·257
a(20) 1441560489708423064780800 2^16·3^5·5^2·7^3·11·13·41·83·331·65537
a(21) 21623407345626345971712000 2^16·3^6·5^3·7^3·11·13·41·83·331·65537
a(22) 6149822195556228360752824320000 2^18·3^5·5^4·7^2·11·41·79·83·157·313·331·65537
a(23) 9224733293334342541129236480000 2^17·3^6·5^4·7^2·11·41·79·83·157·313·331·65537

Sources: a(18), a(20), a(22), a(23) appear among the examples of Moews & Moews (1998); a(19) is from Cohen (1990). All are listed in David Moews' aliquot database (https://djm.cc/aliquot-database/aliquot-database-.1.txt).

These values were previously known as infinitary perfect numbers but not known to be consecutive terms. An exhaustive search over Fermi-Dirac factorizations (N = product of distinct p^(2^k), with prod(q+1) = 2*prod(q)) shows that no other infinitary perfect number exists up to a(23) = 9224733293334342541129236480000. In particular the intervals between a(17) and a(18), a(18) and a(19), ..., a(22) and a(23) are empty.

The smallest known infinitary perfect number above a(23) is 92247332933343425411292364800000 = 10*a(23) (Moews & Moews 1998); it is not known whether other terms lie between a(23) and that value.
```

## Optional PROG entry

A checker for listed terms (verifies σ∞(N) = 2N from the ordinary
factorization; trial division, so intended for checking terms, not for
searching). Tested: all 23 terms pass.

```
(Python)
def isok(n):  # sigma_infinity(n) == 2n
    s, m, p = 1, n, 2
    while p*p <= m:
        e = 0
        while m % p == 0: m //= p; e += 1
        k = 0
        while e:
            if e & 1: s *= p**(2**k) + 1
            e >>= 1; k += 1
        p += 1 if p == 2 else 2
    if m > 1: s *= m + 1
    return s == 2*n
```

Optionally also upload the search program (`ipn_search.c`, C with
primesieve and pthreads) as a linked file, e.g. `a007357.c`, so others
can reproduce the result.

## Evidence behind the claims

**Independent check of the values.** For each of the 23 terms,
σ∞(N) was computed from N's prime factorization — σ∞(p^e) is the
product of p^(2^k)+1 over the binary digits 2^k of e — and equals 2N.
This check does not use the search code.

**Exhaustiveness.** The search program enumerates every N up to a bound
X as a product of distinct Fermi-Dirac primes, with exact integer
arithmetic for the remaining ratio and pruning rules whose soundness is
documented at the top of `ipn_search.c` and in `README.md`. Any case
that would exceed the program's integer width is counted as a warning
instead of being skipped, so a run that reports `warnings 0` covered
its whole search space.

- **a(22):** exhaustive sweep of (10^28, 6149822195556228360752824320000],
  completed 2026-09-30 across two machines. All 308,792 subtree tasks
  are complete in both final checkpoints, with `warnings 0`; the
  solutions found are exactly a(1)–a(22).
- **a(23):** exhaustive sweep up to 9224733293334342541129236480000,
  completed 2026-10-01. The 311,775 subtree tasks were split across
  several runs on two machines; combining all checkpoints leaves 0 tasks
  unfinished, every run reported `warnings 0` throughout, and the
  solutions found are exactly a(1)–a(23). About 4.6×10^11 search nodes
  in total.

**Validation of the program.** Before any new claim, the program
reproduces the 17 terms published before 2026 exactly, and every
revision used in these runs reproduced the a(20) search with an
identical node count (1,340,384,551 nodes).

## Caveats to be upfront about

- **Single implementation.** The exhaustiveness proof rests on one
  program, written for this project. It has not yet been reproduced by
  an independent implementation. Offering the source (above) is the
  best answer if an editor asks.
- **Probabilistic primality for large candidates.** Candidate
  components above 2^64 are tested with 40-base Miller–Rabin. This
  cannot cause a missed term — Miller–Rabin never rejects a true prime —
  and any spurious "solution" it produced would fail the σ∞ check above.
  All 23 terms have components at most 65537.
- **Distributed endgame.** The final stretches of both searches were
  split across machines and merged by checkpoint. One deliberately
  stopped run counts only the tasks its checkpoint marks as finished. A
  development build with a checkpoint-accounting bug ran briefly during
  the a(22) search; its checkpoints were discarded and replaced with
  saved copies from before it ran.
- **Do not submit 92247332933343425411292364800000 as a(24).** It is
  the next *known* value, but the interval between a(23) and it has not
  been searched. Proving it is a(24) needs a sweep about 10 times
  larger, still within the program's checked arithmetic range (10^32).

## How to submit

1. Log in to oeis.org and open A007357 for editing.
2. Make the DATA, COMMENT, EXTENSIONS and LINK changes above.
3. Upload the new `a007357.txt` and the b-file (and optionally the
   program) through the entry's file upload.
4. In the notes to editors, mention that a(22)–a(23) were previously
   known (Moews & Moews 1998) and that the new contribution is the
   exhaustive confirmation; offer the source code.
