# a000230

Computes terms of [OEIS A000230](https://oeis.org/A000230): a(0) = 2; for
n ≥ 1, a(n) is the smallest prime p such that the gap between p and the
next prime is exactly 2n.

## Build

Requires [primesieve](https://github.com/kimwalisch/primesieve)
(`brew install primesieve`) and pkg-config.

```
make
```

## Usage

```
./a000230 -n <terms> [-t <threads>] [-v]
```

- `-n <terms>` — number of terms to produce, a(0) … a(terms−1)
- `-t <threads>` — worker threads (default: all CPUs)
- `-v` — progress on stderr

Output is `n a(n)` per line, the same format as the OEIS b-file, so it can
be diffed directly against <https://oeis.org/A000230/b000230.txt>.

## How it works

Primes are enumerated over expanding blocks (starting at 10⁷, doubling
each round) until every requested gap has been seen. Each block is split
into contiguous chunks, one per thread; a thread walks its chunk with a
`primesieve_iterator` and records the first prime in the chunk that
begins each gap of interest. A gap is attributed to the prime that starts
it, and every prime lies in exactly one chunk, so chunk and block
boundaries are handled exactly. Per-thread results are merged in
ascending chunk order, so the smallest prime for each gap wins.

On a 14-core machine: 150 terms in ~0.15 s, 200 terms in ~2.2 s
(sieving to 1.6×10¹¹). Runtime is dominated by the largest term, which
grows rapidly with n.
