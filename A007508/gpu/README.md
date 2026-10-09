# A007508 on the GPU (DGX Spark, GB10)

`a007508_cuda.cu` is the CUDA version of the twin prime counter in
`../cpu/a007508.c` (the CPU tool is `../a007508.c` in the repository root).
Same k-space layout, same decades, same options and the same checkpoint log
format, so chunks computed on the GPU and on CPUs can be merged with `cat`.

## Build and run (Spark: sm_121, CUDA 13, primesieve in ~/local)

    make                                   # or: make DEFS="-DNSEG=32 -DNSUB=32"
    ./a007508_cuda -v -f 16 -a 1177209242304 16     # a(16) from a(15), timing per chunk
    ./a007508_cuda -l gpu.log 17                    # checkpointed; rerun to resume
    ./a007508_cuda -f 18 -a 808675888577436 -r 0:450 -l m1.log 18   # half a decade

Options: `-s chunk` (numbers, default 10^max(10,n-4) capped at 10^14),
`-f k -a A`, `-l file`, `-r A:B`, `-v`. Decades 1..4 are done on the host by
trial division; the GPU needs the chunk start to be above sqrt(chunk end).

## Design

The bitmap is the packed twin layout (bit 3k+j for p = 30k + r_j, r in
{11,17,29}), with 1 = composite so that marking is an `atomicOr`. A chunk
[K0, K1) is cut into super-chunks of NSUB sub-chunks of NSEG segments of
SEG_K = 2^18 k-units (a segment is 96 KB of shared memory, 7.9M numbers;
a sub-chunk is a 6.3 MB bitmap that stays in the 24 MB L2).

Sieving primes are tiered by size:

| primes | where | how |
|---|---|---|
| 7..97 | kernel_segment | seven periodic presieve patterns ORed into the segment |
| 101..QSPLIT (2^21) | kernel_segment | shared-memory atomics; warp per prime below 4096 (30 lanes = 6 progressions x 5 strides), thread per prime above; first hit per segment by 32-bit Barrett reduction |
| QSPLIT..QL2 (2^26) | kernel_large, per sub-chunk | one thread per prime walks its wheel hits into the L2 bitmap; 8-byte state (next hit, wheel index) |
| above QL2 | kernel_generate, per super-chunk | one thread per prime appends 4-byte hit records to a list per sub-chunk through shared-memory bins flushed in bulk; kernel_apply later ORs each sub-chunk's records into the L2 bitmap |

Two streams overlap the marking kernels of sub-chunk t+1 with the segment
kernel of sub-chunk t (double-buffered bitmaps). The segment kernel loads
the L2 bitmap words into shared memory, ORs the patterns, sieves the small
primes, then counts the zero bits below the chunk end.

Why the tiers: shared-memory atomics are nearly free; L2 atomics run at
about 2.3e10/s on the GB10 (about 43 ps per hit); per-segment record lists
cost about 16 bytes of DRAM traffic per hit (about 100 ps at the bandwidth
we reach), so lists are used only to defer hits of the largest primes to
the right sub-chunk, after which they are applied with L2 atomics. The
generator visits each large prime once per super-chunk (2^28 k-units), so
the per-prime state traffic stays under 1 byte per number even at 10^20.

Setting NSEG (power of two), NSUB, QSPLIT, QL2, BLOCKDIM at compile time
(`-D...`) changes the geometry; `-DAPPLY_L2=0` switches the large primes to
per-segment lists (kernel_distribute, slower on the GB10).

## Correctness

Counts identical to the CPU tool on identical chunks at 10^13, 10^15, 10^17,
10^18, 1.01e19 and on the last 3e10 numbers below 10^20 (the CPU tool itself
matches OEIS a(1)..a(15) and primesieve up to 2^64).

## Throughput (GB10, 10^12-number chunks, numbers per second)

| range | GPU (this) | Spark CPU, 20 threads | M3 Max, 14 threads |
|---|---|---|---|
| 10^13 | 7.0e11 | | |
| 10^15 | 3.8e11 | 3.8e10 | 6.8e10 |
| 10^17 | 2.1e11 | | |
| 10^18 | 1.6-1.8e11 | | ~5e10 |

The GPU is bound by L2 atomic throughput plus the record traffic of the
largest primes, both of which scale with the hit density 6 sum 1/q over
the sieving primes, so throughput falls slowly with the height of the range.

Time for a full decade at these rates:

| n | numbers | GPU time |
|---|---|---|
| 16 | 9e15 | ~7 h |
| 17 | 9e16 | ~5 days |
| 18 | 9e17 | ~2 months |
| 19 | 9e18 | ~2 years |
| 20 | 9e19 | ~25 years |

So one Spark reaches a(17) comfortably and a(18) with patience; a(20) needs
on the order of 100 GB10-class GPUs for a few months, or fewer, faster
cards: the design is bandwidth-bound, and a desktop RTX 5090 (1.8 TB/s,
96 MB L2) should run it several times faster than the GB10.

## Splitting across machines

Same log format as the CPU tool: give each machine a `-r A:B` chunk range of
the decade and its own log, concatenate the logs, run once more without `-r`.
