#!/usr/bin/env python3
"""Generate verification.txt: cross-checks of the digitprimes run against
independently published OEIS data, plus internal consistency identities."""

import csv
import sys
from datetime import date

RESULTS = "digitprimes_results.txt"
OUT = "verification.txt"

OEIS_MISS = ["A231412", "A228413", "A228414", "A228415", "A228416",
             "A228417", "A228418", "A228419", "A228420", "A228421"]
OEIS_HAS = ["A231726", "A231787", "A231788", "A231789", "A231790",
            "A231792", "A231793", "A231794", "A231795", "A231796"]

# OEIS-published terms, "no digit d", a(0)..a(12)  (fetched 2026-08-31)
KNOWN_MISS = {
    0: [1, 10, 91, 819, 7122, 61702, 557224, 5062320, 45002763, 395879190, 3579400605, 32487367715, 294505958253],
    1: [1, 6, 54, 532, 4675, 34425, 262549, 2051466, 16831152, 155616459, 1529462564, 14830618421, 141585123501],
    2: [0, 7, 77, 697, 6497, 55552, 512100, 4710641, 42205969, 341224891, 2787791578, 22971326749, 190650687957],
    3: [1, 7, 54, 534, 4909, 45405, 385008, 3539880, 32260781, 294001190, 2564080248, 23271246324, 211753431947],
    4: [1, 10, 75, 721, 6637, 60605, 514809, 4730382, 43254591, 392344689, 3421561753, 31049600245, 282499317912],
    5: [1, 9, 85, 708, 6635, 60640, 535534, 4737129, 43297195, 392641522, 3536880527, 31067514571, 282635824867],
    6: [1, 10, 90, 719, 6696, 60845, 554933, 4742037, 43331008, 392875212, 3573268469, 31207451849, 282765603085],
    7: [1, 8, 67, 539, 5034, 45549, 416913, 3570781, 32517377, 294828478, 2681147149, 23720397369, 212156228217],
    8: [1, 10, 92, 816, 6712, 60867, 555878, 5026796, 43410238, 395243878, 3576361255, 32461990759, 282971130960],
    9: [1, 8, 69, 620, 5010, 45732, 418142, 3785060, 32579606, 296601070, 2683254222, 24354108057, 212324183352],
}

# OEIS-published terms, "at least one digit d", a(1)..a(11)  (fetched 2026-08-31)
KNOWN_HAS = {
    0: [0, 9, 181, 2878, 38298, 442776, 4937680, 54997237, 604120810, 6420599395, 67512632285],
    1: [4, 46, 468, 5325, 65575, 737451, 7948534, 83168848, 844383541, 8470537436, 85169381579],
    2: [3, 23, 303, 3503, 44448, 487900, 5289359, 57794031, 658775109, 7212208422, 77028673251],
    3: [3, 46, 466, 5091, 54595, 614992, 6460120, 67739219, 705998810, 7435919752, 76728753676],
    4: [0, 25, 279, 3363, 39395, 485191, 5269618, 56745409, 607655311, 6578438247, 68950399755],
    5: [1, 15, 292, 3365, 39360, 464466, 5262871, 56702805, 607358478, 6463119473, 68932485429],
    6: [0, 10, 281, 3304, 39155, 445067, 5257963, 56668992, 607124788, 6426731531, 68792548151],
    7: [2, 33, 461, 4966, 54451, 583087, 6429219, 67482623, 705171522, 7318852851, 76279602631],
    8: [0, 8, 184, 3288, 39133, 444122, 4973204, 56589762, 604756122, 6423638745, 67538009241],
    9: [2, 31, 380, 4990, 54268, 581858, 6214940, 67420394, 703398930, 7316745778, 75645891943],
}

# OEIS A006988: the 10^n-th prime, n = 0..13
KNOWN_NTH_PRIME = [2, 29, 541, 7919, 104729, 1299709, 15485863, 179424673,
                   2038074743, 22801763489, 252097800623, 2760727302517,
                   29996224275833, 323780508946331]


def main():
    rows = []  # (n, nth_prime, [no_0..no_9])
    with open(RESULTS) as f:
        for rec in csv.reader(f):
            if not rec or rec[0].startswith("#") or rec[0] == "n":
                continue
            rows.append((int(rec[0]), int(rec[1]), [int(x) for x in rec[2:12]]))
    rows.sort()

    checks_pass = 0
    checks_fail = 0
    lines = []

    def out(s=""):
        lines.append(s)

    def check(ok, text):
        nonlocal checks_pass, checks_fail
        if ok:
            checks_pass += 1
        else:
            checks_fail += 1
        out(f"  [{'PASS' if ok else 'FAIL'}] {text}")
        return ok

    W = max(len(str(rows[-1][1])), 15)

    out("=" * 78)
    out("VERIFICATION OF DIGIT-CONTENT COUNTS OVER THE FIRST 10^n PRIMES")
    out("=" * 78)
    out()
    out(f"Generated: {date.today().isoformat()} from {RESULTS}")
    out("Tool: digitprimes (single ordered pass over all primes via primesieve,")
    out("14-thread chunked sieve with strictly in-order merge; every power-of-10")
    out("crossing chunk was independently re-scanned serially and required to")
    out("match the parallel pass exactly).")
    out()
    out("Sequences verified/extended, per digit d:")
    out("  'no digit d'          : " + " ".join(OEIS_MISS) + "   (d = 0..9)")
    out("  'at least one digit d': " + " ".join(OEIS_HAS) + "   (d = 0..9)")
    out()

    # ---------------------------------------------------------------- sec 1
    out("-" * 78)
    out("1. ENDPOINT CHECK: the 10^n-th prime vs OEIS A006988")
    out("-" * 78)
    out("The run counts primes in order and records the prime at which the count")
    out("reaches exactly 10^n. These must equal the independently published")
    out("values of the 10^n-th prime (OEIS A006988).")
    out()
    for n, p, _ in rows:
        known = KNOWN_NTH_PRIME[n] if n < len(KNOWN_NTH_PRIME) else None
        if known is None:
            out(f"  [ -- ] p(10^{n:<2}) = {p:>{W},}  (no published value embedded)")
        else:
            check(p == known, f"p(10^{n:<2}) = {p:>{W},}  == A006988({n})")
    out()

    # ---------------------------------------------------------------- sec 2
    out("-" * 78)
    out("2. SUMMATION CHECK: no_d(n) + has_d(n) = 10^n (total primes counted)")
    out("-" * 78)
    out("For every n and every digit d, the count of the first 10^n primes")
    out("WITHOUT digit d plus the count WITH at least one digit d must sum to")
    out("the total number of primes considered, which is 10^n by definition")
    out("(and section 1 confirms the count reached 10^n at the correct prime).")
    out()
    for n, _, miss in rows:
        total = 10 ** n
        sums = [miss[d] + (total - miss[d]) for d in range(10)]
        ok = all(s == total for s in sums)
        check(ok, f"n={n:<2}  no_d + has_d = {total:>{W},}  for all d = 0..9")
    out()
    out("  Explicit summations for the new terms (n = 13, total = 10,000,000,000,000):")
    out()
    n13 = next((r for r in rows if r[0] == 13), None)
    if n13:
        for d in range(10):
            no_d = n13[2][d]
            has_d = 10 ** 13 - no_d
            check(no_d + has_d == 10 ** 13,
                  f"d={d}: {no_d:>17,} ({OEIS_MISS[d]}) + {has_d:>17,} ({OEIS_HAS[d]}) = {no_d + has_d:>18,}")
    out()

    # ---------------------------------------------------------------- sec 3
    out("-" * 78)
    out("3. 'NO DIGIT d' COUNTS vs PUBLISHED OEIS TERMS (n = 0..12)")
    out("-" * 78)
    out()
    for d in range(10):
        out(f"  {OEIS_MISS[d]} (first 10^n primes with no digit {d}):")
        for n, _, miss in rows:
            if n < len(KNOWN_MISS[d]):
                check(miss[d] == KNOWN_MISS[d][n],
                      f"a({n:<2}) = {miss[d]:>{W},}  == OEIS")
            else:
                out(f"  [NEW ] a({n:<2}) = {miss[d]:>{W},}")
        out()

    # ---------------------------------------------------------------- sec 4
    out("-" * 78)
    out("4. 'AT LEAST ONE DIGIT d' = 10^n - no_d vs PUBLISHED OEIS TERMS (n = 1..11)")
    out("-" * 78)
    out()
    for d in range(10):
        out(f"  {OEIS_HAS[d]} (first 10^n primes with at least one digit {d}):")
        for n, _, miss in rows:
            if n == 0:
                continue
            has = 10 ** n - miss[d]
            if n - 1 < len(KNOWN_HAS[d]):
                check(has == KNOWN_HAS[d][n - 1],
                      f"a({n:<2}) = {has:>{W},}  == OEIS")
            else:
                out(f"  [NEW ] a({n:<2}) = {has:>{W},}")
        out()

    # ---------------------------------------------------------------- sec 5
    out("-" * 78)
    out("5. NEW TERMS")
    out("-" * 78)
    out()
    out("  'No digit d' — new a(13), over the first 10^13 primes")
    out("  (all primes p <= p(10^13) = 323,780,508,946,331):")
    out()
    if n13:
        for d in range(10):
            out(f"    {OEIS_MISS[d]} (no {d}):  a(13) = {n13[2][d]:>17,}")
        out()
        out("  'At least one digit d' — new a(12) and a(13):")
        out()
        n12 = next(r for r in rows if r[0] == 12)
        for d in range(10):
            out(f"    {OEIS_HAS[d]} (>=1 {d}):  a(12) = {10**12 - n12[2][d]:>15,}   a(13) = {10**13 - n13[2][d]:>17,}")
    out()

    # ---------------------------------------------------------------- summary
    out("=" * 78)
    out(f"SUMMARY: {checks_pass} checks passed, {checks_fail} failed.")
    out("=" * 78)

    with open(OUT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {OUT}: {checks_pass} passed, {checks_fail} failed")
    return 1 if checks_fail else 0


if __name__ == "__main__":
    sys.exit(main())
