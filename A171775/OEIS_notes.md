# Proposed OEIS update from the A171775 computation (not yet submitted)

## A171775  a(n) = smallest number M such that there exist bases b_2, ..., b_n with M
##          written in base b_k a k-digit palindrome for all k = 2..n.

New terms: **a(10) = 4722366482869645213696 = 2^72** and **a(11) = 1237940039285380274899124224 = 2^90**.

DATA (offset 1):

    1, 3, 5, 52, 130, 1885, 1073741824, 4398046511104, 72057594037927936, 4722366482869645213696,
    1237940039285380274899124224

(a(11) has 28 digits; the DATA line stays within the usual limit.)

%e a(10) = 2^72 = 4722366482869645213696: in bases 255, 511, 1023, 2047, 8191, 32767,
   524287 it is the palindrome C(9,j), C(8,j), 4*C(7,j), 2^6*C(6,j), 2^7*C(5,j),
   2^12*C(4,j), 2^15*C(3,j) (lengths 10 down to 4), it is 13548844 15833512 13548844 in
   base 18669329 and 11 in base 2^72-1. These are the smallest bases of each length.

%e a(11) = 2^90: in bases 511, 1023, 2047, 4095, 16383, 65535, 524287, 8388607 it is the
   palindrome C(10,j), C(9,j), 4*C(8,j), 2^6*C(7,j), 2^6*C(6,j), 2^10*C(5,j), 2^14*C(4,j),
   2^21*C(3,j) (lengths 11 down to 4) and 2^28*(1,2,1) in base 2^31-1 (smallest bases).

%F (existing conjecture) a(n) = 2^((n-1)(n-2)) for n >= 7: now verified for n <= 11.

%E a(10) from Jeff Sponaugle, Sep 29 2026 (exhaustive search: every 10-digit palindrome
   <= 2^72 in every base, 3.7*10^12 numbers, was tested for the 9-digit condition by a
   modular method; 135382 numbers are both 10- and 9-digit palindromes and all except
   2^72 fail at length 8).
%E a(11) from Jeff Sponaugle, Sep 30 2026 (exhaustive search: every 10-digit palindrome
   <= 2^90 in every base, 7.0*10^15 numbers, was tested for the 11-digit condition by
   solving for its two innermost digits with a lookup table; 524157 numbers are both 11- and
   10-digit palindromes and all except 2^90 fail at length 9).

b-file: b171775.txt (n = 1..11).

## Verification chain

* Two implementations (C/pthreads and CUDA) of the search, run over [1, 2^72], agree on
  every statistic (5126601984459 inner steps, 31468592 three-digit matches, 135382
  (n-1)-palindromes, order-independent checksum of the survivors) and find exactly one
  solution, 2^72.
* Both reproduce a(7), a(8), a(9) = 2^30, 2^42, 2^56 by exhaustive searches from 1.
* The fast filter was compared with a naive scan on 400+ random blocks (n = 7..12) and
  checked state-by-state in a debug build; an independent pure-Python enumeration
  (verify_a171775.py) gives identical statistics on n = 7..10 test ranges.
* The digits of 2^72 and 2^90 in the listed bases were checked in Python from the definition.
* a(11): three complete runs over [1, 2^90] give identical statistics (5386659681911 lookups,
  1235755652260 matches, 524157 survivors, checksum 0a5d53d568af3f34, one solution 2^90):
  the swapped-2D tool on an Apple M1 Pro (clang) and on a DGX Spark CPU (gcc), and its
  CUDA port on the Spark's GB10 GPU (separate kernel code with 32-bit limb arithmetic).
* The swapped-2D method reproduces the 1D a(10) search exactly (same 135382 survivors and
  checksum) and a(8), a(9); it agrees with a naive scan on random blocks for n = 8..12 and
  with a state-by-state debug build; the 1D tool (the other enumeration) gives the same 211
  survivors and checksum on the n = 11 window [600^9, 600^9 + 3*600^8].
