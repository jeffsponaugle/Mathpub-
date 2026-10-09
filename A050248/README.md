# A050248 explorer

[OEIS A050248](https://oeis.org/A050248): integer averages of the first k primes —
the values S(k)/k whenever k divides S(k) = p₁ + p₂ + … + p_k.
The k values are [A045345](https://oeis.org/A045345), the sums S(k) are
[A050247](https://oeis.org/A050247).

## Build

Requires [primesieve](https://github.com/kimwalisch/primesieve)
(`brew install primesieve`; set `PREFIX` if it isn't in `/opt/homebrew`).

```
make
```

## Run

```
./a050248 [-q] [prime_limit]
```

- `prime_limit` — sieve all primes p ≤ limit; accepts `1e12` notation. Default 1e10.
- `-q` — suppress the progress lines printed to stderr every 250M primes.

Each term is printed as it is found, with n, k, p_k, S(k), and the average.

## Notes

- The running sum overflows 64 bits near p ≈ 6×10⁹, so the accumulator is
  `unsigned __int128` (safe past the sum of all primes below 10¹⁹).
- Throughput is roughly 175M primes/s single-threaded on an Apple Silicon Mac:
  p ≤ 1e11 (terms 1–10) takes ~25 s, p ≤ 1e12 on the order of 4–5 min.
- All 16 known terms lie below p ≈ 2.5×10¹⁷ and a(17) > 1.25×10¹⁷ (search
  frontier as of 2022), so finding a *new* term needs a distributed effort —
  this program is for exploration and verification of the known range.
