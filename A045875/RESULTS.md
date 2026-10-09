# A045875 a(18) and a(19) — search results and provenance

## a(19) = 3238684956  (found 2026-09-26)

- 2^3238684956 has **974,941,319 decimal digits**
- Contains a run of **exactly 19 consecutive '9' digits**
- Run position: digits 875,306,593..875,306,611 from the least significant
  digit (99,634,707 digits in from the most significant digit)
- Digit context (MSD side): `...73365556368 9999999999999999999 4732931778...`
- Exact run length 19 establishes **a(20) > 3238684956**.

Minimality (contiguous, no gaps; below 1973572200 covered by the a(18) proof
since that run is exactly 18):

| range (m) | machine | log |
|---|---|---|
| 1,973,572,200 – 2,460,000,000 | atom1 | a19.log |
| 2,460,000,000 – 2,920,000,000 | atom2 | a19-r2460.log |
| 2,920,000,000 – 3,238,684,956 | atom1 (hit) | a19-r2920.log |

Independently verified by `verify_a18.c` (direct GMP recompute): digit count,
all 19 run digits, and exact run length confirmed. Bonus coverage for a
future a(20) hunt: atom2 also cleared [3,725,000,000 – ~3,897,000,000) of
19-runs (a19-r3725.log) — note a(20) coverage must still re-scan that range
for 20-runs only if desired; the a20 search simply starts at 3238684957.
Campaign wall time: ~2.2 days on two DGX Sparks; the hit landed at the ~63rd
percentile of the predicted distribution (ratio a(19)/a(18) ≈ 1.64).

Suggested OEIS lines:
- Term: `3238684956`
- Comment: `a(20) > 3238684956. - Jeff Sponaugle, Sep 26 2026`
- Extension: `a(19) from Jeff Sponaugle, Sep 26 2026`

# A045875 a(18) — search result and provenance

**a(18) = 1973572199** — the smallest m such that the decimal representation
of 2^m contains 18 consecutive identical digits.

Found 2026-09-24 by Jeff Sponaugle, using CPU/CUDA search tools in this
repository, on two NVIDIA DGX Spark (GB10) machines plus one x86 host.

## The hit

- 2^1973572199 has **594,104,431 decimal digits**
- Contains a run of **exactly 18 consecutive '4' digits**
- Run position: digits 394,539,910..394,539,927 from the least significant
  digit (199,564,503 digits in from the most significant digit)
- Digit context (MSD side): `...04265984511 444444444444444444 35479973487...`
- Since the run has exact length 18, this also establishes **a(19) > 1973572199**.

## Verification

1. Found by the CUDA tool (`a045875_gpu.cu`), which was validated bit-exactly
   against the independent CPU tool (`a045875.c`): md5-identical state
   checkpoints after identical doubling ranges, and exact reproduction of
   a(16) = 106892452 and a(17) = 356677212 (digit, run length, and position).
2. Independently confirmed by `verify_a18.c`: a direct GMP computation of
   2^1973572199 (binary shift + divide-and-conquer radix conversion — no code
   shared with the search), checking digit count, all 18 digits of the run,
   and that the run is exactly 18 (bounded by '1' and '3').
3. Detection pipeline sanity-checked at scale mid-search: an n=17 probe from
   m=1,760,864,744 found a 17-run within 4.2M steps, consistent with density
   statistics at 530M digits.

## Minimality: complete coverage of [a(17), a(18))

a(18) >= a(17) = 356677212 because an 18-run contains a 17-run.

| range (m) | machine / tool | log |
|---|---|---|
| 356,677,212 – 430,610,150 | house1, CPU tool | house1 session |
| 430,610,150 – 900,000,000 | atom1 GPU | a18.log |
| 900,000,000 – 1,161,000,000 | atom2 GPU | a18-r900.log |
| 1,161,000,000 – 1,500,000,000 | atom1 GPU | a18-r1161.log |
| 1,500,000,000 – 1,872,500,000 | atom2 GPU | a18-r1500.log |
| 1,872,500,000 – 1,965,000,000 | atom1 GPU | a18-gap.log |
| 1,965,000,000 – 1,973,572,199 | atom1 GPU (hit) | a18-r1965.log |

All ranges completed with "no run of 18 found" (or the final RESULT), with
checkpointed, resumable state; every start state self-checked against an
independently computed 2^m mod 10^9. Logs and checkpoints reside in
/home/jbs/A045875 on atom1 (10.1.30.36) and atom2 (10.1.30.37).

## Method (summary)

State: decimal expansion of 2^m as base-10^9 uint32 limbs, doubled in place.
Doubling in base 10^9 has purely local carries (carry into limb i is
`old a[i-1] >= 5e8`), and after j doublings limb i depends only on limbs
i-j..i. The GPU kernel exploits this: each thread takes a 64-limb span plus
a 28-limb halo and performs 28 doublings entirely in registers (carries as
bitmasks), checking every intermediate value; memory traffic per doubling
drops 28x. Any run of >= 17 identical digits must contain an aligned 9-digit
repdigit limb (value d*111111111), so detection is one multiply-and-compare
per limb per step, with rare candidates confirmed exactly on the host by
recomputing a small limb window at the exact intermediate step. Limbs are
stored tile-transposed in GPU memory for fully coalesced access.

Throughput: ~38,000 doublings/s at 107M digits, ~7,300/s at 590M digits per
GB10 (vs ~770/s for the 16-thread CPU tool on the original x86 host).
Total search wall time: about one day across the two Sparks.

## Suggested OEIS submission lines

- Term: `1973572199` (extends the DATA line after 356677212)
- Comment: `a(19) > 1973572199. - Jeff Sponaugle, Sep 24 2026`
- Extension: `a(18) from Jeff Sponaugle, Sep 24 2026`
