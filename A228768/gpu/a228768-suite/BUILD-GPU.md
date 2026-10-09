# Building the GPU sieve on DGX Spark

## Required library: the CUDA Toolkit, and nothing else

No Thrust, no CUB, no third-party dependency. Compaction is a warp-aggregated
atomic (a dozen lines), which avoids version-coupling the build to a template
library. DGX OS 7 (Ubuntu 24.04) ships CUDA; check with `nvcc --version`.

GB10 is compute capability **12.1** (`sm_121`), which requires **CUDA 12.9 or
newer**. Older toolkits fail with:

    nvcc fatal : Value 'sm_121' is not defined for option 'gpu-architecture'

Diagnose with `nvcc --version`, `ls -d /usr/local/cuda*`, `nvidia-smi`. Often
a newer toolkit is already installed and PATH points at an older one.

Everything is aarch64 — the CPU library uses the portable prefilter path
automatically (no x86 inline asm).

## Build

```sh
cc   -O3 -c -DA228768_LIB a228768x.c -o a228768lib.o
nvcc -O3 -arch=sm_121 a228768_gpu.cu a228768lib.o -o a228768gpu -lpthread -lm
```

### If the toolkit is older than 12.9

Emit PTX only and let the driver JIT it for the real device at load time. This
needs no toolkit knowledge of sm_121 at all, and costs a one-time JIT of a few
seconds on first run (cached afterwards in ~/.nv/ComputeCache):

```sh
nvcc -O3 -gencode arch=compute_120,code=compute_120 \
     a228768_gpu.cu a228768lib.o -o a228768gpu -lpthread -lm
```

A cubin built for `sm_120` will also load on a 12.1 device (cubins are
compatible with higher minor revisions of the same major architecture), so
`-arch=sm_120` usually works too. The PTX form above is safer because the
driver targets the actual hardware.

### Portable build line

```sh
ARCH=$(nvcc --list-gpu-code 2>/dev/null | grep -q sm_121 \
       && echo "-arch=sm_121" \
       || echo "-gencode arch=compute_120,code=compute_120")
nvcc -O3 $ARCH a228768_gpu.cu a228768lib.o -o a228768gpu -lpthread -lm
```

## Build with no GPU (algorithm check on any machine)

The kernels are written as `(block, thread)` functions so they also compile as
plain C++ and run serially. Results are bit-for-bit identical to the CUDA path.

```sh
cc  -O3 -c -DA228768_LIB a228768x.c -o a228768lib.o
c++ -O3 -DEMULATE -x c++ a228768_gpu.cu -x none a228768lib.o -o a228768emu -lpthread -lm
```

## Verify before trusting it

```sh
./a228768gpu --verify --blocks 1 --block-log 21 --sieve-limit 65536
```

Recomputes one block's survivors by honest trial division and compares as sets.
Must print `match`. Then the end-to-end check against a known term:

```sh
./a228768gpu --nbase 14 --start 19775643791000000 --limit 19775643792000000
```

Must print `HIT: 19775643791090909`.

## Run

```sh
./a228768gpu --nbase 15 --start 3e19 --limit 5*15^16 \
             --block-log 28 --gpu-bases 8 --threads 20 --resume w1.ckpt
```

`--resume` writes a watermark every 10 s and reads it at startup, so a reboot
costs at most a few blocks. The watermark only advances past a block whose CPU
cascade has finished, so nothing is skipped.

### Tuning

| flag | meaning | note |
|---|---|---|
| `--block-log K` | GPU block = 2^K numbers | 28 measured best on GB10 |
| `--gpu-bases N` | cascade stages on the GPU | 8 balances GPU and CPU |
| `--gpu-primes N` | small primes per stage | 10 default |
| `--gpu-threads T` | threads per CUDA block | 512 default |
| `--threads T` | CPU cascade pool | 20 on GB10 |
| `--sieve-only` | skip the cascade | raw sieve throughput |
| `--no-presieve` | disable the 3..19 pattern | A/B timing |

The end-of-run report states which side is the wall. If it says CPU-BOUND,
raise `--gpu-bases`; if GPU-BOUND with the cascade nearly free, lower it.

`--threads` sizes the CPU-side cascade pool (pthreads, not OpenMP —
nvcc's frontend strips `#pragma omp`, which would silently serialize it).
Defaults to the online CPU count; on GB10 that is 20.
