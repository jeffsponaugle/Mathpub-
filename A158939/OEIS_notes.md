# Proposed OEIS updates from the A158939 computation (not yet submitted)

## A158939  First primes followed by sequences of exactly n monotonic increasing prime gaps.

New term: **a(16) = 17293451238695141**.

DATA (offset 0): 7, 3, 2, 17, 347, 2903, 15373, 128981, 1319407, 17797517, 94097537,
6927837557, 48486712783, 968068681511, 1472840004017, 129001208165717, 17293451238695141

%e a(16) = 17293451238695141 is followed by the 16 increasing gaps
   2, 4, 6, 8, 22, 26, 30, 34, 36, 44, 46, 54, 56, 64, 86, 108 (the next gap is 46).
%C a(17) > 2*10^16. - Jeff Sponaugle, Sep 19 2026 (exhaustive scan to 2*10^16; prime count = pi(2*10^16) = 547863431950008)
%E a(16) from Jeff Sponaugle, Sep 19 2026 (GPU sieve; every prime below a(16) examined,
   prime count checked against pi(10^16)).

b-file: b158939.txt (n = 0..16).

## A229832  First term of smallest sequence of n consecutive weak primes.

New term: **a(15) = 17293451238695143** (the prime following A158939(16); the 15 consecutive
weak primes are 17293451238695143, ...147, ...153, ...161, ...183, ...209, ...239, ...273,
...309, ...353, ...399, ...453, ...509, ...573, ...659).

## A133697  Smallest k such that P(k)/P(k+1) > ... > P(k+n+1)/P(k+n+2).

New term: **a(14) = 475618519121221** = pi(17293451238695141), the index of A158939(16)
among the primes (independently computed with primecount 8.7; the GPU scan's running
prime count gives the same value).

## Verification chain

* Run length 16 at 17293451238695141 re-derived from the definition by three independent
  implementations: a158939 verify (deterministic Miller-Rabin next-prime search), GMP
  mpz_nextprime, verify_run.py (pure Python).
* Minimality: GPU scan of [0, 10^16) (prime count = pi(10^16) = 279238341033925 exactly),
  CPU scan of [10^16, 10^16 + 41615360) (no run longer than 9), GPU scan of
  [10^16 + 41615360, 17293451238695141].
* Both tools reproduce a(1)..a(15) and pi(a(n)) = A133697(n-2) for n = 2..15, and agree
  exactly (prime counts, first occurrences, run-length histograms) on [0, 10^13) and on a
  window around a(15).
