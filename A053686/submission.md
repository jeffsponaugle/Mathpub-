# A053686 — OEIS submission notes

Checked against the live entries on Oct 3 2026.

## Where the entries stand

The Sep 17 2026 submissions are all in and approved; nothing has changed since.

**A053686** revision #28 (Sep 17 2026):

```
%S 2,4,6,14,34,36,52,86,132,154,250,336,906
%C The record gap 906 (first after 218209405436543) occurs again after 543684371469023, whose
   next prime is 543684371469929, before the record 916 at 1189459969825483. The five other
   records (778, 804, 806, 916 and 924) do not repeat, so a(14) >= 1132. - Jeff Sponaugle, Sep 17 2026
%K more,nonn
%E a(13) from Jeff Sponaugle, Sep 17 2026
```

**A133788** revision #13 (Sep 17 2026): DATA ends `..., 3842610773, 218209405436543`,
`%E a(13) from Jeff Sponaugle, Sep 17 2026`.

**A085237** revision #44 (Sep 17 2026): b-file `b085237.txt` for n = 1..86 (terms 64..79
credited to Charles R Greathouse IV), ending 778, 804, 806, 906, 906, 916, 924, 1132;
`%E a(80)-a(86) from Jeff Sponaugle, Sep 17 2026`. The displayed DATA was trimmed by the
editors to the usual length (it now ends at 474); the b-file carries the rest. The Python
program's `while r < 778` is unchanged (harmless).

Not touched by this work and unchanged: A005250 (#290, Sep 30 2026) and A002386 (#185,
Jul 2026) — the 80-record table in the tool still matches them — and A053695.

An OEIS search for `5,13,31,293,8467,12853,25471` returns nothing, so the second
occurrences (item 1) are not in the database.

## What can still be submitted

Items 1 and 3 are cheap and worth doing now; 2 is a reasonable companion sequence; 4 and 5
are optional; 6 is the next real computation.

### 1. Comment in A053686 listing the second occurrences  (recommended)

The primes after which each term occurs for the second time are new data (from
`a053686 scan 0 1e10` for a(1..12), verified by the next-prime search in `selftest`, and
the run on the Mac Studio for a(13), verified four ways — see README):

```
%C A053686 The second occurrences of a(1)-a(13) are the gaps following the primes 5, 13, 31, 293, 8467, 12853, 25471, 338033, 1561919, 11113933, 428045491, 4275912661, 543684371469023. Only 4, 6, 14, 36 and 154 occur more than twice before the next record (4 also after 19; 6 also after 47, 53, 61, 73, 83; 14 also after 317; 36 also after 14107; 154 also after 15203977). - _Jeff Sponaugle_, Oct 03 2026
```

If the editors prefer a sequence over a list in a comment, see item 2a.

### 2. New companion sequences  (optional)

**2a. Second occurrences**, parallel to A133788 ("Primes where the record gaps in A053686
first appear"):

```
%N Primes where the record gaps in A053686 appear for the second time.
%S 5,13,31,293,8467,12853,25471,338033,1561919,11113933,428045491,4275912661,543684371469023
%O 1,1
%K nonn,hard,more
%C a(n) is the smallest prime p > A133788(n) such that the next prime is p + A053686(n); by definition of A053686 it lies before the start A002386(k+1) of the next record gap, where A005250(k) = A053686(n).
%e a(2) = 13: the record gap 4 first occurs between 7 and 11 (A133788(2) = 7) and again between 13 and 17, before the next record gap 6 between 23 and 29.
%Y Cf. A053686, A133788, A005250, A002386, A085237, A001223.
```

**2b. Multiplicity of each record gap** (how many times A005250(n) occurs as a gap between
consecutive primes before the next record; A053686 lists the A005250(n) with a(n) >= 2):

```
%N Number of times the n-th record prime gap A005250(n) occurs as a gap between consecutive primes below A002386(n+1), the start of the next record gap.
%S 1,2,3,7,1,3,1,1,1,2,3,1,2,1,2,1,1,1,1,2,1,3,1,1,1,1,1,1,2,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1
%O 1,2
%K nonn,hard,more
%C a(n) is the multiplicity of A005250(n) in A085237, so 1 + Sum_{k=1..n} a(k) is the index of A005250(n+1) in A085237 (e.g. A085237(86) = A005250(64) = 1132). A005250(n) is a term of A053686 iff a(n) >= 2. a(64) (for the record gap 1132) is not known; deciding it needs every prime gap up to A002386(65) = 43841547845541059.
%e a(4) = 7: the record gap 6 occurs between 23 and 29 and again after 31, 47, 53, 61, 73 and 83, seven times in all, before the record gap 8 between 89 and 97.
%Y Cf. A005250, A002386, A085237, A053686, A133788.
```

The 63 terms come from `scan_0_1e11.txt` (n <= 39), `scan_1e13.txt` (35..54) and the
Studio run (`scan_1.69e15_summary.txt`, 55..63); the "repeats" column of each table is
a(n) - 1.

### 3. Program link  (recommended)

The directory is public in the repository, so the LINKS entry (for A053686, and the same
line fits A133788 and A085237):

```
%H A053686 Jeff Sponaugle, <a href="https://github.com/jeffsponaugle/Mathpub-/tree/main/A053686">C program, Python verifier and notes</a>
```

### 4. EXAMPLE for A053686  (optional; the entry has none)

```
%e A053686 4 is a term: the record gap 4 between 7 and 11 occurs again between 13 and 17 before the next record gap 6 (23 to 29). 8 is not a term: the record gap 8 between 89 and 97 does not occur again before the record 14 between 113 and 127.
```

### 5. A085237 comment on the next term  (optional)

```
%C A085237 a(87) is 1132 or 1184, depending on whether the record gap 1132 (after A002386(64) = 1693182318746371) occurs again before the next record gap 1184 (after 43841547845541059); this needs every prime gap up to 4.4*10^16. - _Jeff Sponaugle_, Oct 03 2026
```

### 6. a(14): does 1132 repeat?  (the next computation)

Deciding record 64 means sieving [1693182318746371, 43841547845541059), 4.2*10^16
numbers: about 6 days on the Mac Studio (8*10^10/s, slowing beyond 10^16), or roughly
two days on the 96-core x86 box if it is available. The command continues from the
Studio's checkpoint (its `a053686.state` covers [4302407359, 1693182318746372)):

    caffeinate -i ./a053686 scan r35 43841547845541060 -S a053686.state -H -i 60 | tee -a a053686.log

Outcome either way is submittable: a `REPEAT` line gives a(14) = 1132 (and A133788(14) =
1693182318746371, A085237(87) = 1132); none gives "a(14) >= 1184" and A085237(87) = 1184.
The near-tail model (`tail_a053686.py`) rates 1132 (not divisible by 3) as a less likely
repeater than 906 was.

## Files backing the submission

* `scan_1e13.txt`, `scan_1e13.log` — M1 Pro leg [4302407359, 10^13): records 35..54, no repeats.
* `scan_1.69e15_summary.txt` — final table of the Mac Studio run to A002386(64): a(13) = 906.
* `scan_0_1e11.txt` — the known part [0, 10^11) with histogram lines (second occurrences, multiplicities).
* `a053686.c`, `verify_a053686.py`, `tail_a053686.py`, `README.md` — tool, independent verifier, model.
* The Studio's `a053686.log` / `a053686.state` hold the full run and the checkpoint for item 6.
