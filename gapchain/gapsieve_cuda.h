/*
 * gapsieve_cuda.h -- the interface between gapsieve_cuda_host.c (the host
 * side, which includes gapsieve.c) and gapsieve_cuda.cu (the CUDA kernels and
 * the device memory they use).  C linkage, so either side compiles alone.
 */
#ifndef GAPSIEVE_CUDA_H
#define GAPSIEVE_CUDA_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CU_MAXSLOT 8            /* units in flight, at most */
#define CU_MAXK 1024            /* member offsets per shape, at most */

/* one unit of one shape, as the kernels see it */
struct cu_params {
    uint64_t base_lo, base_hi;  /* the unit's start */
    uint64_t lim;               /* the largest candidate offset in range */
    uint64_t segsize;           /* a unit is nseg segments of this length */
    uint32_t W, bmodW, A, ngroups, nps, k;
    uint32_t cap;               /* capacity of each candidate list */
    uint32_t nseg;
    uint32_t dense;             /* pattern groups sieved before compacting */
};

/* select the device; its name goes in name.  0 on success */
int cu_init(char *name, size_t n);

/* the tables every unit shares: pattern groups (P, first, n, 64-bit word
   offset for each), the sieve primes with W^-1 mod q and CRT coefficients */
int cu_tables(uint32_t ngroups, const uint32_t *ginfo, uint32_t nps,
              const uint32_t *q, const uint32_t *qinv, const uint32_t *crt);

/* one shape's patterns (npat 32-bit words), admissible residues and member
   offsets */
int cu_shape(int d, const uint32_t *pat, size_t npat, const uint32_t *res,
             uint32_t A, const uint32_t *offs, uint32_t k);

/* buffers for one unit in flight: up to unit segments, lists of cap entries */
int cu_slot(int slot, uint32_t unit, uint32_t cap);

/* reallocate a slot's candidate lists with cap entries */
int cu_grow(int slot, uint32_t cap);

/* queue a unit: bq holds each segment's base mod each sieve prime; p[d] is
   NULL for a shape the unit doesn't need; split rounds of one member each,
   then one kernel for the rest.  Units run one after another, in order */
int cu_submit(int slot, const uint32_t *bq, uint32_t nseg, uint32_t nps,
              const struct cu_params *p[2], int split);

/* wait for a slot: 0 done, -1 failed.  stat[d] is the sieve's raw count
   (above cap: the lists overflowed), nfin[d] the survivors */
int cu_wait(int slot, uint32_t stat[2], uint32_t nfin[2]);

/* shape d's survivors after cu_wait: n of them */
const uint64_t *cu_list(int slot, int d, uint32_t n);

/* the base-2 strong test on count numbers (lo, hi pairs): for checking the
   kernels' arithmetic against the CPU's */
int cu_sprp_batch(const uint64_t *n, size_t count, uint8_t *ok);

#ifdef __cplusplus
}
#endif

#endif
