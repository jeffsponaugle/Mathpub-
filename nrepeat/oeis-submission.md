# OEIS edit text for A331881 / A331882 (extension a(9)-a(12))

Submit via the "edit" link on each entry (requires an OEIS account; changes
go through editorial review). Both entries carry keywords `base,nonn,hard,more`
— keep them ("more" is still true: a(13) is unknown).

Ratio sanity check of the new completion positions (A331882 grows ~9-12x per
term, matching the (n!)^(1/n)*10^(n-1) birthday estimate):
4519908297/481426671 = 9.39, 53506468331/4519908297 = 11.84,
579330690478/53506468331 = 10.83.

---

## A331881 — a(n) is the first n-digit substring to occur n times in the decimal expansion of the fractional part of Pi.

**DATA** (replace the term list):

```
1, 26, 446, 2796, 86538, 872117, 1591292, 66416662, 858275944, 5762246467, 47982672136, 100085029093
```

**EXAMPLE** (append, matching the entry's existing style — the n=9 and n=12
position lists below are from the verified run reports; fill in the n=10/n=11
lists from results.txt on house1 if desired, or include only a(9) as the
illustrative example to keep the section short):

```
a(9) = 858275944, which first appears 9 times within the first 481426671 digits (at positions 8364516, 108721072, 118923087, 209578974, 347121287, 365354466, 397701128, 426812266 and 481426663).
```

**EXTENSIONS**:

```
a(9)-a(12) from Jeff Sponaugle, Sep 01 2026
```

Suggested comment for the pink-box note to editors (not part of the entry):
computed by an exhaustive sliding-window scan over Pi digit files with exact
direct-indexed occurrence counters followed by an independent verification
pass recomputing all occurrence positions of each candidate; a(9) was
reproduced identically from two independently computed digit sources (MPFR
const_pi to 6e8 digits, and a 1e13-digit Pi file), and the n <= 8 terms were
reproduced from scratch as a self-test. Program available from the author.

---

## A331882 — a(n) is the number of digits in the decimal expansion of the fractional part of Pi needed to contain n occurrences of an n-digit substring.

**DATA** (replace the term list):

```
1, 22, 219, 1805, 25499, 168882, 3566679, 29325629, 481426671, 4519908297, 53506468331, 579330690478
```

**EXAMPLE** (append, matching the entry's existing style):

```
a(9) = 481426671, because we need 481426671 digits to have the first 9-digit substring ('858275944') appearing 9 times.
```

**EXTENSIONS**:

```
a(9)-a(12) from Jeff Sponaugle, Sep 01 2026
```

---

## Supporting data (occurrence start positions, 1-based, from the run reports)

n=9  (858275944):  8364516 108721072 118923087 209578974 347121287
                   365354466 397701128 426812266 481426663
n=10 (5762246467): 62083470 334368659 734076310 901367000 1611279195
                   2368347516 3645033154 4282821718 4472382985 4519908288
n=11 (47982672136): 10501601989 12403732092 13269952849 15713942497
                    16424571944 22902234759 25685796054 31751231100
                    42257680350 51109092955 53506468321
n=12 (100085029093): 51738833731 68671221941 75604597502 160932526320
                     231367395165 314009479133 324889612204 445320004727
                     447112780643 468314848852 526877014657 579330690467
```
