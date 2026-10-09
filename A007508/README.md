# A007508 — Number of twin prime pairs below 10^n

    a(n) = #{ p : p and p+2 are both prime, p+2 < 10^n }

    2, 8, 35, 205, 1224, 8169, 58980, 440312, 3424506, 27412679,
    224376048, 1870585220, 15834664872, 135780321665, 1177209242304,
    10304195697298, 90948839353159, 808675888577436, 7237518093734545

Source: https://oeis.org/A007508. a(20) is not known.

## Build and run

    brew install primesieve
    make
    ./a007508 13                          # a(1)..a(13), all cores
    ./a007508 -f 15 -a 135780321665 15    # a(14) known: sieve only [10^14, 10^15)
    ./a007508 -l run.log 16               # checkpoint every chunk; rerun to resume

Options:

    -t threads   worker threads (default: all online CPUs)
    -k kernel    twin (default) | ps | iter          see "Kernels"
    -s chunk     chunk size in numbers (default 10^max(9, n-4), capped at 10^13)
    -S log2      twin-kernel segment size, log2 of k-units (default 22 = 1.5 MB)
    -f k -a A    start at decade k with a(k-1) = A
    -l file      append-only chunk log; a rerun with the same log resumes
    -r A:B       only chunks A <= idx < B of one decade (split across machines)

Output is one line per decade: `n a(n) [pairs in [10^(n-1),10^n)] [seconds]`.

## How it works

For p > 5 a twin pair (p, p+2) has p = 30k + r with r in {11, 17, 29}: those
are the only residues mod 30 where both p and p+2 are coprime to 2, 3 and 5.
The program indexes candidates by k = p/30 and a class j in {0,1,2}, and since
10^n = 30 K_n + 10, "p < 10^n" is exactly "k < K_n" with K_n = (10^n - 10)/30.
K_20 is about 3.3e18, so everything up to 10^20 fits in 64-bit integers even
though 10^20 itself does not. The pairs (3,5) and (5,7) are added by hand.

Each decade [K_(n-1), K_n) is cut into chunks. Worker threads take chunk
indices from an atomic counter (dynamic load balancing). A pair is attributed
to the chunk containing k, so nothing is counted twice at chunk boundaries.

### Kernels

* **twin** (default): a segmented sieve over the three twin classes only, so
  the bitmap holds 3 bits per 30 numbers, 2.7x denser than a normal mod-30
  sieve. A bit for candidate p survives only if neither p nor p+2 has a prime
  factor, so for each sieving prime q both progressions p ≡ 0 and p ≡ -2
  (mod q) are cleared: six arithmetic progressions per prime per period, versus
  eight for a full prime sieve. Sieving primes are split in tiers:
  - 7..97: precomputed periodic patterns ANDed into each segment (presieve).
  - 101..2048: unrolled 8-hit loops run over 32 KB sub-blocks so the
    read-modify-writes stay in L1.
  - up to segment/64: the same unrolled loops over the whole segment.
  - up to one segment: a wheel walk (one record per prime, six deltas of the
    form a·(q/30)+b), unrolled by one full period.
  - larger: a bucket sieve. Each prime is an 8-byte record (offset, class,
    wheel index, q/30) filed in the bucket of the segment it hits next.
  Works to 10^20; only this kernel can do n = 20.
* **ps**: `primesieve_count_twins` per chunk. 64-bit only, so n <= 19.
* **iter**: walk primes with `primesieve_iterator` and count gaps of 2.
  Slowest, kept as an independent check.

The sieving primes for a chunk are generated with primesieve up to
sqrt(30·K1 + 1). That setup is repeated per chunk, which is why the default
chunk grows with n (10^9 at n <= 13, 10^(n-4) up to 10^13): at 10^19 a chunk
narrower than 10^10 spends more time on setup than on sieving.

### Correctness checks

All three kernels agree with OEIS for a(1)..a(13), and the twin kernel
reproduced a(14) and a(15) in a 4.6 h run on 14 cores. The twin and ps kernels
give identical counts on identical chunks at 10^13, 10^15, 10^17, 10^18 and
1.83·10^19 (just under 2^64), across segment sizes 2^16..2^24 and chunk sizes
1000..10^11. Above 2^64 only the twin kernel runs; on the last 3·10^10
numbers below 10^20 it counts 18,683,392 pairs against a Hardy–Littlewood
prediction of 18.68 million, agreement to 2·10^-4.

### Performance (Apple M-series, 14 cores)

Single thread, per 10^10 numbers, chunk width 10^10:

| range | twin  | ps    |
|-------|-------|-------|
| 10^13 | 1.9 s | 1.7 s |
| 10^17 | 2.5 s | 3.7 s |
| 10^18 | 3.1 s | 4.2 s |

14 threads on 1.4·10^12 numbers at 10^14: twin 20.6 s, ps 28.7 s, i.e.
6.8·10^10 numbers per second. Per-decade wall time from that rate:

| n  | numbers | time (14 cores)      |
|----|---------|----------------------|
| 13 | 9e12    | ~2.5 min             |
| 14 | 9e13    | 20.4 min (measured)  |
| 15 | 9e14    | 4.2 h (measured)     |
| 16 | 9e15    | ~1.7 days            |
| 17 | 9e16    | ~3 weeks             |
| 18 | 9e17    | ~7 months            |
| 19 | 9e18    | ~6 years             |
| 20 | 9e19    | ~70 years            |

So a(20) is roughly 1,000 core-years with this code, or a few weeks on a
1,000-core cluster. Memory for the twin kernel is about 8 bytes per sieving
prime per thread: 0.1 GB at 10^17, 1.2 GB at 10^19, and 2.8 GB measured for
a chunk at 10^20. Per-thread memory times threads must fit in RAM, so a(20)
wants machines with a few GB per core.

## GPU version

`gpu/a007508_cuda.cu` is a CUDA port for the DGX Spark (GB10); see
[gpu/README.md](gpu/README.md). Same options and log format, so GPU and CPU
chunks merge. Measured on the GB10: 3.8e11 numbers/s at 10^15 and about
1.8e11 at 10^18, roughly 3 to 5 times this laptop's 14 cores.

## Resuming, checkpointing, splitting

Each decade only depends on a(n-1), and chunks are independent, so work is
never repeated:

    ./a007508 -f 14 -a 15834664872 14      # a(13) known
    ./a007508 -l run.log 15                # log every finished chunk and decade
    ./a007508 -l run.log 15                # after a crash: same command resumes

The log is append-only text. `k n idx chunk count` records one finished chunk
(chunk in k-units), `d n a(n)` a finished decade. On start the program reads
the log, skips finished decades, and inside the current decade skips finished
chunks. The chunk size must match; a mismatch is rejected. A partial last line
from an interrupted write is trimmed.

To split one decade across machines, give each a disjoint chunk range and its
own log, then concatenate the logs and run once more without `-r`:

    machine1$ ./a007508 -f 16 -a 1177209242304 -r 0:4500    -l m1.log 16
    machine2$ ./a007508 -f 16 -a 1177209242304 -r 4500:9000 -l m2.log 16
    $ cat m1.log m2.log > all.log && ./a007508 -f 16 -a 1177209242304 -l all.log 16

A run with `-r` prints a partial count and writes the `d` line only when every
chunk of the decade is in the log.

## Are there shortcuts that avoid enumerating every prime?

No exact one is known. For pi(x) there are combinatorial
(Meissel–Lehmer–Deléglise–Rivat–Gourdon) and analytic (Lagarias–Odlyzko,
Platt/Büthe/Franke) methods that run in about O(x^(2/3)) or O(x^(1/2+eps))
without touching every prime. Nothing analogous exists for twin primes: the
count depends on the joint primality of p and p+2, and no inclusion–exclusion
over "numbers with no small factors" collapses into a recursion the way the
Legendre/Meissel identity does for pi(x). Every record value of this sequence
(Brent, Nicely, Sebah, Oliveira e Silva, Chaffin's 2026 a(19)) came from
sieving all the way up. The Hardy–Littlewood estimate

    pi_2(x) ~ 2 C_2 * integral_2^x dt / (ln t)^2,   C_2 = 0.6601618...

is accurate to about seven digits at 10^18 but is a heuristic, not a proof,
and cannot give the exact integer.

What does exist are constant-factor shortcuts, all still sieves, and the twin
kernel uses them: sieve only the 3 of 8 residue classes that can host a twin,
clear "p and p+2" jointly, presieve the small primes, keep segments in cache,
bucket the large primes, and split the range across cores and machines.

Partial sieving followed by primality tests on the survivors does not help:
survivor density only falls like 1/(ln B)^2 in the sieve bound B, so at 10^20
testing survivors costs several times more than finishing the sieve.

## Verification advice for a(20)

a(19) was corrected once on OEIS, so redundancy matters. Cheap checks: log the
count per chunk and compare each with the Hardy–Littlewood prediction for that
interval (a corrupted chunk stands out by hundreds of standard deviations),
recompute a random sample of chunks with `-k ps` or on another machine, and
use ECC memory. The `-k ps` kernel cannot be used at 10^20 itself, so
cross-check there with a second independent twin-sieve implementation.
