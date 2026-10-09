# OEIS submissions from the gapchain searches

Everything the gapchain searches can contribute to the OEIS, with the text for
each submission ready to paste. Checked against oeis.org on 2026-10-06.

The first six submissions are in the OEIS. The rest have not been submitted
yet. Every value below was computed and cross-checked as described in its
section, and [README.md](README.md) has the full search record behind each term.

## At a glance

| # | Entry | Submission | Status |
| --- | --- | --- | --- |
| 1 | A016045 | a(15) | Approved 2026-09-24 |
| 2 | A263049 | a(15) | Approved 2026-09-24 |
| 3 | A263049 | a(16), a comment, keyword `hard` | Approved 2026-09-28 |
| 4 | A016045 | a(16), a(17) | Approved 2026-09-27 |
| 5 | A349121 | New terms a(15)-a(17) | Approved 2026-10-06 |
| 6 | A090870 | New terms a(15)-a(17) | Approved 2026-10-05 |
| 7 | A349121 | Correct the PARI program, and the name's wording | Ready; run it in gp first |
| 8 | A016045 | Comments, cross-references, a Python program | Ready |
| 9 | A263049 | Comments, cross-references, a Python program | Ready |
| 10 | New | Initial primes of 10 consecutive primes with gaps 2, 4, ..., 18 | Ready |
| 11 | New | Initial primes of 10 consecutive primes with gaps 18, ..., 4, 2 | Ready |
| 12 | 14 related entries | Cross-references back to A016045 and A263049 | Optional |
| 13 | A016045, A263049 | The next terms, a(18) and a(17) | Needs new searches |

Of what's left, items 10 and 11 add new sequences and item 7 fixes an error,
so they matter most. Items 8 and 9 each fit in one edit of their entry.

## Done: in the OEIS

All six were submitted by Jeff Sponaugle. Revision numbers and times are
from each entry's history, in US Eastern time.

| Entry | Terms | Submitted | Approved |
| --- | --- | --- | --- |
| A016045 | a(15) = 1397398433200922807 | revision 32, Sep 24 00:26 | revision 39, Sep 24 17:23, by Michael De Vlieger |
| A263049 | a(15) = 253253149671986953 | revision 37, Sep 24 00:28 | revision 39, Sep 24 01:25, by Sean A. Irvine |
| A263049 | a(16) = 134491669675212201371, with a comment and the keyword `hard` | revision 40, Sep 27 12:59 | revision 50, Sep 28 00:16, by Sean A. Irvine |
| A016045 | a(16) = 254345176718302423991, a(17) = 307223174679659356577 | revision 40, Sep 27 20:20 | revision 43, Sep 27 23:30, by Michael De Vlieger |
| A349121 | a(15)-a(17), the same values as A016045's | revision 48, Oct 04 14:13 | revision 51, Oct 06 00:07, by Sean A. Irvine |
| A090870 | a(15) = 34287539062666038, a(16) = 5533761473972926789, a(17) = 6656852126755763340 | revision 16, Oct 04 16:11 | revision 19, Oct 05 08:32, by Michael De Vlieger |

The five chains, with the gaps on either side:

| Term | Primes | Gaps | Gap before | Gap after |
| --- | --- | --- | --- | --- |
| A016045 a(15) = 1397398433200922807 | 16 | 2, 4, ..., 30 | 10 | 60 |
| A016045 a(16) = 254345176718302423991 | 17 | 2, 4, ..., 32 | 114 | 36 |
| A016045 a(17) = 307223174679659356577 | 18 | 2, 4, ..., 34 | 24 | 30 |
| A263049 a(15) = 253253149671986953 | 16 | 30, 28, ..., 2 | 2 | 16 |
| A263049 a(16) = 134491669675212201371 | 17 | 32, 30, ..., 2 | 10 | 48 |

**How the terms were checked.** `gapsieve` proves each member prime with
Miller-Rabin on the 13 smallest prime bases, which is exact below 3.3×10²⁴, and
checks that every number between members is composite. `verify_chain.py`
retests every number from each chain's first prime to its last, separately with
the same Miller-Rabin and with Baillie-PSW. Both find exactly the members.

**Why each is the smallest.**

- **A016045 a(15):** the whole range from a(14) up was searched twice, on two
  different machines and builds. The sequence is non-decreasing, so nothing below
  a(14) can qualify.
- **A016045 a(16) and a(17):** every start below each chain was searched. The
  last range, 2×10²⁰ to 3×10²⁰, was covered in 27 pieces, checked end to end
  from their checkpoints.
- **A263049 a(15):** it is a(14) − 30. A 16-prime decreasing chain contains a
  15-prime one starting 30 above it, so no 16-prime chain can start below
  a(14) − 30, and one starts there.
- **A263049 a(16):** every start from a(15) − 32 up to the chain was searched.

**The A349121 and A090870 terms.** [A349121](https://oeis.org/A349121) is the
"exactly n gaps" version of A016045: the chain must stop after n gaps. A prime
followed by exactly n such gaps is followed by at least n, so A349121(n) >=
A016045(n), with equality when the chain for A016045(n) stops. The chains for
A016045's a(15)-a(17) are followed by gaps of 60, 36 and 30, not the 32, 34 and
36 that would continue them, so A349121's a(15)-a(17) are the same three
numbers. [A090870](https://oeis.org/A090870) is primePi(A016045(n)). Its new
terms were computed with primecount 8.7 using Gourdon's algorithm, then again
with its Deleglise-Rivat algorithm, and the two agree. The same program
reproduces the entry's a(11) and a(14).

During A349121's review, David A. Corneth noted that the name's "but not for
k>n" would be more precise as "but not for k=n+1". Jeff Sponaugle agreed in the
discussion. Item 7 can make that change.

## 7. A349121: correct the PARI program and the name

One edit to [A349121](https://oeis.org/A349121). As of revision 51, which
added the new terms, neither has changed.

The entry's PARI program computes A016045, not A349121. `isok()` checks that p
is followed by at least n gaps 2, 4, ..., 2n, but not that the chain stops
there. So it gives a(6) = 128981 instead of 665111.

Now:

```
(PARI) isok(p, n) = my(q=p); for (k=1, n, my(r = p+k+k^2); if (nextprime(q+1) != r, return (0)); q=r); return(1);
a(n) = my(p=2); while (!isok(p, n), p=nextprime(p+1)); p; \\ _Michel Marcus_, Nov 09 2021
```

Proposed: the last statement of `isok()` returns whether the next gap is
something other than 2n+2. After the loop, q is the chain's last prime,
p + n + n^2.

```
(PARI) isok(p, n) = my(q=p); for (k=1, n, my(r = p+k+k^2); if (nextprime(q+1) != r, return (0)); q=r); nextprime(q+1) != p+(n+1)+(n+1)^2;
a(n) = my(p=2); while (!isok(p, n), p=nextprime(p+1)); p; \\ _Michel Marcus_, Nov 09 2021
```

**Testing:** no PARI/GP was available, so both versions were tested as exact
Python transcriptions with sympy. The current program gives 3, 5, 17, 347,
13901, 128981, 128981 for n = 1..7. The corrected one gives 3, 5, 17, 347,
13901, 665111, 128981, matching the entry's data. **Run the corrected version in
gp before submitting.** Keep Michel Marcus's signature, and explain the change in
the discussion. The entry's Mathematica program already checks the stop and is
fine.

**NAME**, with David A. Corneth's wording from the review of revision 50:

```
a(n) is the smallest prime p, such that p + k + k^2 are consecutive primes for 0 <= k <= n, but not for k = n+1.
```

A note for the discussion:

> Corrected the PARI program: isok() did not check that the chain stops, so it
> returned A016045(n), e.g. 128981 for n = 6. Also changed "but not for k>n" to
> "but not for k = n+1" in the name, as David A. Corneth suggested in the review
> of revision 50.

## 8. A016045: comments, cross-references and a Python program

One edit to [A016045](https://oeis.org/A016045). The entry has no comments and
no programs yet. Its only cross-reference is A263049.

**COMMENTS:**

```
a(n) <= a(n+1), with equality exactly when the chain for a(n) continues with a gap of 2n+2, as for a(6) = a(7) = 128981. - _Jeff Sponaugle_, Oct 06 2026

The chains for a(15), a(16) and a(17) are followed by gaps of 60, 36 and 30, so each stops after exactly n gaps, and these are also terms of A349121. - _Jeff Sponaugle_, Oct 06 2026

a(n) is the first term of the sequence of initial primes of n+1 consecutive primes with gaps 2, 4, ..., 2n: A022004 (n = 2), A078847 (n = 3), A190814 (n = 4), A190817 (n = 5), A190819 (n = 6), A190838 (n = 7) and A281448 (n = 8). - _Jeff Sponaugle_, Oct 06 2026
```

If the new sequence in item 10 is accepted, add "and A?????? (n = 9)" to the
last comment.

**PROG:**

```
(Python)
from collections import deque
from sympy import nextprime
def A016045(n):
    d, p = deque([2], maxlen=n+1), 2
    while True:
        p = nextprime(p); d.append(p)
        if len(d) > n and all(d[i+1]-d[i] == 2*i+2 for i in range(n)): return d[0]
print([A016045(n) for n in range(1, 9)]) # _Jeff Sponaugle_, Oct 06 2026
```

Tested with sympy 1.12: it reproduces a(1)-a(8) exactly, in about 4 minutes.

**CROSSREFS**, replacing "Cf. A263049.":

```
Cf. A263049 (decreasing gaps).
Cf. A090870 (indices), A349121 (exactly n gaps).
Cf. A022004, A078847, A190814, A190817, A190819, A190838, A281448.
```

[A094749](https://oeis.org/A094749) (rows of n successive primes whose
differences rise by 2 from any first gap) and
[A084299](https://oeis.org/A084299) also cite A016045. They could be added too,
but they're further from it.

## 9. A263049: comments, cross-references and a Python program

One edit to [A263049](https://oeis.org/A263049). It already has the a(16)
comment from revision 40. Its only cross-reference is A016045.

**COMMENTS:**

```
a(n) >= a(n-1) - 2n, because the pattern for a(n) has the pattern for a(n-1) starting at its second prime. Equality holds exactly when the prime before a(n-1) is a(n-1) - 2n, as for a(10) = a(9) - 20 and a(15) = a(14) - 30. So the sequence is not monotonic. - _Jeff Sponaugle_, Oct 06 2026

a(n) is the first term of the sequence of initial primes of n+1 consecutive primes with gaps 2n, ..., 4, 2: A022005 (n = 2), A078855 (n = 3), A289907 (n = 4), A286891 (n = 5), A290161 (n = 6), A290162 (n = 7) and A290264 (n = 8). - _Jeff Sponaugle_, Oct 06 2026
```

The inequality could go in FORMULA instead, as "a(n) >= a(n-1) - 2n.", with the
rest kept as a comment.

**PROG:**

```
(Python)
from collections import deque
from sympy import nextprime
def A263049(n):
    d, p = deque([2], maxlen=n+1), 2
    while True:
        p = nextprime(p); d.append(p)
        if len(d) > n and all(d[i+1]-d[i] == 2*(n-i) for i in range(n)): return d[0]
print([A263049(n) for n in range(1, 9)]) # _Jeff Sponaugle_, Oct 06 2026
```

Tested with sympy 1.12: it reproduces a(1)-a(8) exactly, in about 1 minute.

**CROSSREFS**, replacing "Cf. A016045.":

```
Cf. A016045 (increasing gaps).
Cf. A022005, A078855, A289907, A286891, A290161, A290162, A290264.
```

## 10. New sequence: 10 consecutive primes with gaps 2, 4, ..., 18

The OEIS has this family for 3 to 9 primes (A022004, A078847, A190814,
A190817, A190819, A190838, A281448), but not for 10. A search for its first
terms, and for the gap list in the name, finds nothing.

**NAME:**

```
Initial primes of 10 consecutive primes with consecutive gaps 2, 4, 6, 8, 10, 12, 14, 16, 18.
```

**DATA** (21 terms, 1..21):

```
2426256797,6430890287,8518049207,55065405671,55373581421,60590486081,66945470477,76566117071,78067026071,95748657617,104217487301,111058349531,120546569177,135648105611,137168442221,137376420947,188343072341,209839688957,240441136151,292790049317,294026931551
```

**OFFSET:** 1,1

**COMMENTS:**

```
a(1) = A016045(9). Subsequence of A281448.
```

**EXAMPLE:**

```
a(1) = 2426256797: the ten consecutive primes 2426256797, 2426256799, 2426256803, 2426256809, 2426256817, 2426256827, 2426256839, 2426256853, 2426256869, 2426256887 have gaps 2, 4, 6, 8, 10, 12, 14, 16, 18.
```

**CROSSREFS:**

```
Cf. A016045, A078847, A190814, A190817, A190819, A190838, A281448.
```

**KEYWORD:** nonn. **AUTHOR:** _Jeff Sponaugle_, Oct 06 2026

**Data and checks:** `gapchain_ps --length 10 --end 1e13`, which enumerates the
actual primes with primesieve, finds 185 terms below 10¹³. Below 10¹²,
`gapsieve --length 10` finds the same 39 terms with an independent sieve.
The full list could be uploaded as a b-file ("Table of n, a(n) for n =
1..185"). The DATA line alone is enough for a first submission.

<details>
<summary>b-file: all 185 terms below 10¹³</summary>

```
1 2426256797
2 6430890287
3 8518049207
4 55065405671
5 55373581421
6 60590486081
7 66945470477
8 76566117071
9 78067026071
10 95748657617
11 104217487301
12 111058349531
13 120546569177
14 135648105611
15 137168442221
16 137376420947
17 188343072341
18 209839688957
19 240441136151
20 292790049317
21 294026931551
22 301595503481
23 354138215267
24 355076247257
25 361945174751
26 371248830761
27 378981593801
28 399020146697
29 473787509537
30 539313280547
31 540205177637
32 557149355507
33 562914638381
34 622255739321
35 632751717797
36 638368991081
37 710785340621
38 748172623397
39 799015982747
40 1055035747847
41 1063524061421
42 1091666916407
43 1092074843927
44 1117233465431
45 1131021834617
46 1293496877831
47 1427339937041
48 1431786828737
49 1431832975481
50 1439251952351
51 1491489049991
52 1544634887921
53 1555046019071
54 1673941031027
55 1715846996321
56 1720866650807
57 1920451125491
58 1984088238551
59 1992939202691
60 2003987481101
61 2131528031441
62 2170341748697
63 2172560720567
64 2172889244621
65 2202352794527
66 2229347130077
67 2236934638847
68 2262400428677
69 2269877834591
70 2313936247937
71 2364724681661
72 2535806437751
73 2547706912187
74 2616352005557
75 2635551337577
76 2641309357991
77 2656174523657
78 2739933923471
79 2760875495507
80 2814278347031
81 2889612975407
82 2914862511287
83 2923403124737
84 3007681931867
85 3015446431577
86 3032568069887
87 3090855886577
88 3132173481881
89 3176175919067
90 3253584760961
91 3481235927441
92 3513428103587
93 3548813907491
94 3552783611681
95 3630243449867
96 3695559614687
97 3717053283761
98 3721556957507
99 3736635776297
100 3759338377631
101 3788065254131
102 3820258647701
103 3844124810621
104 3967772025311
105 4059260464601
106 4115083623017
107 4235992318001
108 4331733220541
109 4394243968187
110 4404562432967
111 4446023043881
112 4549867661147
113 4588905521897
114 4623475760921
115 4650385577261
116 4656625081181
117 4793161871711
118 4840389162881
119 4867763134667
120 4876956073877
121 4880447140421
122 4942073842241
123 4947330421961
124 5055287057267
125 5267425948151
126 5291375776781
127 5364431251847
128 5527196538881
129 5536737361271
130 5605867192967
131 5626191372077
132 5694983051237
133 5742370877051
134 5758872759041
135 5851986024581
136 5867131427147
137 5893681737467
138 5962152574577
139 6008435924651
140 6275500600277
141 6342221970167
142 6378124212827
143 6461161331321
144 6474583441091
145 6508136293511
146 6559992810851
147 6645187171781
148 6774052899071
149 6917086110371
150 7007229783137
151 7321722214841
152 7358584496267
153 7527762257747
154 7535389177697
155 7578867540311
156 7613946616097
157 7661693319647
158 7686939420047
159 7735901468411
160 7816407366611
161 7826193182111
162 7931417733137
163 7940854727207
164 8249982984701
165 8257376897591
166 8286223666997
167 8374848737597
168 8418421533941
169 8521973368121
170 8630325519581
171 8764849284797
172 8919313672727
173 9022071533831
174 9215845939121
175 9238398616871
176 9274382590217
177 9669206777957
178 9714994219367
179 9871465665221
180 9882619372667
181 9923263355891
182 9930238716527
183 9949908106991
184 9954212767841
185 9969449552327
```

</details>

## 11. New sequence: 10 consecutive primes with gaps 18, 16, ..., 2

The mirror family exists for 3 to 9 primes (A022005, A078855, A289907,
A286891, A290161, A290162, A290264), but not for 10. A search for its first
terms finds nothing. As in A290264, a term is the start of the exact 10-prime
pattern, which can be the tail of a longer chain.

**NAME:**

```
Initial primes of 10 consecutive primes with 9 consecutive gaps 18, 16, 14, 12, 10, 8, 6, 4, 2.
```

**DATA** (21 terms, 1..21):

```
245333233,2636566333,24395129833,27797667559,43033031143,45000281449,50283448513,77128250509,89071927639,126454871293,133092893959,152703854539,189815716333,201216712699,227045652859,263021970313,273193446199,304022491873,354831179323,367968734233,372720867073
```

**OFFSET:** 1,1

**COMMENTS:**

```
a(1) = A263049(9). For every term p, p + 18 is in A290264.
```

**EXAMPLE:**

```
a(1) = 245333233: the ten consecutive primes 245333233, 245333251, 245333267, 245333281, 245333293, 245333303, 245333311, 245333317, 245333321, 245333323 have gaps 18, 16, 14, 12, 10, 8, 6, 4, 2.
```

**CROSSREFS:**

```
Cf. A263049, A078855, A289907, A286891, A290161, A290162, A290264.
```

**KEYWORD:** nonn. **AUTHOR:** _Jeff Sponaugle_, Oct 06 2026

**Data and checks:** `gapchain_ps --length 10 --end 1e13 --gaps dec` finds
208 terms below 10¹³. Below 10¹², `gapsieve --length 10 --gaps dec`
finds the same 43.

<details>
<summary>b-file: all 208 terms below 10¹³</summary>

```
1 245333233
2 2636566333
3 24395129833
4 27797667559
5 43033031143
6 45000281449
7 50283448513
8 77128250509
9 89071927639
10 126454871293
11 133092893959
12 152703854539
13 189815716333
14 201216712699
15 227045652859
16 263021970313
17 273193446199
18 304022491873
19 354831179323
20 367968734233
21 372720867073
22 385556696599
23 397788175423
24 399304194259
25 474789504499
26 480957565723
27 495996354943
28 532374878359
29 544895144533
30 570222789889
31 604487514643
32 632716074559
33 704010689989
34 714152709139
35 785260249303
36 788287290529
37 814937653549
38 850792939999
39 877503975493
40 914384386393
41 924363706333
42 961738788253
43 964171318849
44 1017034082089
45 1025477420419
46 1045589373799
47 1058148579133
48 1060024112053
49 1075248763339
50 1167119860723
51 1222319795023
52 1242300642949
53 1243749430753
54 1244815622503
55 1283197201519
56 1288444060993
57 1424294747059
58 1447695738673
59 1486453162723
60 1490620710793
61 1518332364859
62 1534462641493
63 1631933425939
64 1632512557519
65 1650088845703
66 1686301252759
67 1750915191583
68 1823110820143
69 1829860812799
70 1900498928143
71 1909866050419
72 1917636519103
73 1929883664629
74 2021832569209
75 2065916900539
76 2066921573953
77 2166109560253
78 2268973684429
79 2271540402583
80 2315814062083
81 2360213727583
82 2596207860229
83 2604881511853
84 2681911383253
85 2687551926589
86 2690006224603
87 2770543405723
88 2820474455353
89 2894510815033
90 2898286272439
91 3004715312233
92 3038890002703
93 3089516179489
94 3147003166969
95 3271429878469
96 3328006391833
97 3572326617403
98 3591317023033
99 3593970985903
100 3595831172803
101 3734837227819
102 3767312548309
103 3791131114783
104 4047890954209
105 4068549822229
106 4088833624873
107 4147454863333
108 4158447439633
109 4161883482703
110 4169683944253
111 4174955005699
112 4237330221763
113 4424362736533
114 4457560132699
115 4458650279053
116 4493603254753
117 4570069748773
118 4639969536799
119 4730803131199
120 4739138007979
121 4899586655329
122 4917064911133
123 5035651900753
124 5039142594043
125 5043245718523
126 5102550476839
127 5160165974989
128 5253943270963
129 5316194746249
130 5367814560859
131 5522420678989
132 5531013171499
133 5572847135329
134 5615907503659
135 5687711384899
136 5727543831193
137 5768709347503
138 5889062331583
139 5926543474303
140 5942986370983
141 5977232261023
142 5980877226349
143 6001680237019
144 6015336698473
145 6036192670813
146 6046864985089
147 6052536198013
148 6121051543213
149 6132874211293
150 6197085110059
151 6303405234733
152 6331811948929
153 6449657851249
154 6472966143019
155 6551625407563
156 6600437612689
157 6624672843709
158 6624992667529
159 6643215371959
160 6672890199649
161 6770953722349
162 6779361872203
163 6968763976963
164 6981303376213
165 7002869634373
166 7078700337343
167 7152502445503
168 7202952618523
169 7211280487303
170 7218525705883
171 7318925218663
172 7400521980343
173 7526579682283
174 7559724522109
175 7567321161523
176 7572641494819
177 7614495664603
178 7822566799693
179 7864819177693
180 7883315820223
181 7965442771129
182 8058260766253
183 8206538119873
184 8228698540969
185 8415797264803
186 8584306743763
187 8683844330653
188 8777641562449
189 8883475300009
190 8884868970343
191 8915476213669
192 8943717606199
193 9216572864683
194 9248183060209
195 9271188566623
196 9278729831263
197 9307539303529
198 9336610857313
199 9344885516089
200 9494105037013
201 9528923380219
202 9528967352833
203 9590802317563
204 9626375743723
205 9663030156253
206 9698185348849
207 9896129609389
208 9943877415859
```

</details>

## 12. Optional: cross-references back from the related entries

None of these cites A016045 or A263049. Each could get one cross-reference to
the sequence whose term is its first term. This is small and low priority, and
OEIS editors may prefer it done a few at a time.

| Entry | Gaps | a(1) | Add |
| --- | --- | --- | --- |
| A022004 | 2, 4 | 5 | Cf. A016045 |
| A078847 | 2, 4, 6 | 17 | Cf. A016045 |
| A190814 | 2, ..., 8 | 347 | Cf. A016045 |
| A190817 | 2, ..., 10 | 13901 | Cf. A016045 |
| A190819 | 2, ..., 12 | 128981 | Cf. A016045 |
| A190838 | 2, ..., 14 | 128981 | Cf. A016045 |
| A281448 | 2, ..., 16 | 113575727 | Cf. A016045 |
| A022005 | 4, 2 | 7 | Cf. A263049 |
| A078855 | 6, 4, 2 | 31 | Cf. A263049 |
| A289907 | 8, ..., 2 | 1979 | Cf. A263049 |
| A286891 | 10, ..., 2 | 41203 | Cf. A263049 |
| A290161 | 12, ..., 2 | 752251 | Cf. A263049 |
| A290162 | 14, ..., 2 | 5647457 | Cf. A263049 |
| A290264 | 16, ..., 2 | 32465047 | Cf. A263049 |

## 13. Later: the next terms

These need new searches. When they're found, they also extend A349121 and
A090870 (A016045's next term) and the comments above.

- **A016045 a(18)** needs a chain of 19 primes. It is larger than a(17) =
  307223174679659356577, and nothing above a(17) has been searched end to end.
  The fleet's 11 machines run at about 6,450 T/s together, so reaching 10²¹ would
  take about 30 hours, and each further 10²¹ about 43 hours.
- **A263049 a(17)** needs a decreasing chain of 18 primes, and it is larger
  than a(16) = 134491669675212201371. Such a chain has a 17-prime pattern
  starting 34 above its first prime, so it can't start below a(16) − 34.
  - It doesn't start at a(16) − 34, because the gap before a(16) is 10, not 34.
  - Starting between a(16) − 34 and a(16) would put that 17-prime pattern at
    a(16) + 32, the only prime in reach, where the next gap is 30 rather than 32.
  - It doesn't start at a(16) itself, whose first gap is 32, not 34.

  Some ranges above a(16) were searched partway in September, but not end to
  end from a(16).

## Submitting

- Edit an entry from its page while logged in. Put each change in its own field
  (DATA, COMMENTS, PROG, CROSSREFS, EXTENSIONS), and explain it in the discussion
  box. Then set the edit to "proposed".
- Use one edit per entry: item 7 is one A349121 edit, and items 8 and 9 are
  one edit each.
- New sequences are submitted from the OEIS "Contribute" page.
- From the A263049 discussion: an author b-file is only needed when the terms
  don't fit in DATA. Items 7 to 9 need none. For items 10 and 11, a b-file is
  optional.
- Signatures in comments and programs take the form `_Jeff Sponaugle_, Oct 06
  2026`. The drafts here use 2026-10-06; change it to the date of the edit.
- To answer an editor's comment on a draft under review, add a discussion
  comment rather than editing the draft. An edit sends a reviewed draft back to
  "editing", as the A349121 history shows for an earlier edit in 2021.
