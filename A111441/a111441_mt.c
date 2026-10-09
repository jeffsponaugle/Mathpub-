/* A111441 parallel searcher: numbers k such that the sum of the squares of
 * the first k primes is divisible by k.
 *
 * Strategy: process prime-value windows [lo, hi).  Each window is split into
 * one subrange per thread.  Pass 1: every thread independently computes the
 * prime count and 192-bit sum of p^2 over its subrange.  A prefix sum over
 * subranges then yields the exact global prime index n and running sum S at
 * each subrange start.  Pass 2: every thread re-iterates its subrange doing
 * the divisibility checks with correct global state.
 *
 * Only n == 1 or 5 (mod 6) can be a term (no term is divisible by 2 or 3;
 * see OEIS comments), so ~1/3 of indices are tested.
 *
 * Checkpoints (value boundary, n, S) are written after every window so a
 * long search can be killed and resumed.
 *
 * Usage: ./a111441_mt <stop_prime_value> [ckpt_file] [threads] [window_width]
 *   resumes from ckpt_file if it exists, else starts from 0.
 *
 * Build: cc -O3 -o a111441_mt a111441_mt.c -lprimesieve -lpthread
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>
#include <pthread.h>
#include <primesieve.h>

typedef unsigned __int128 u128;
typedef struct { u128 lo; uint64_t hi; } u192;   /* S = hi*2^128 + lo */

static inline void u192_add_u128(u192 *s, u128 x)
{
    s->lo += x;
    if (s->lo < x) s->hi++;
}

static inline void u192_add(u192 *a, const u192 *b)
{
    a->lo += b->lo;
    a->hi += b->hi + (a->lo < b->lo);
}

static inline uint64_t u192_mod(const u192 *s, uint64_t n)
{
    uint64_t rlo = (uint64_t)(s->lo % n);
    if (s->hi == 0) return rlo;
    uint64_t t64  = (uint64_t)((((u128)(UINT64_MAX % n)) + 1) % n);
    uint64_t t128 = (uint64_t)(((u128)t64 * t64) % n);
    u128 r = ((u128)(s->hi % n) * t128) % n;
    return (uint64_t)((r + rlo) % n);
}

typedef struct {
    /* subrange [start, stop) of prime values */
    uint64_t start, stop;
    /* pass 1 outputs */
    uint64_t count;
    u192     sum;
    /* pass 2 inputs: global state just before this subrange */
    uint64_t n0;
    u192     s0;
} chunk_t;

static void *pass1(void *arg)
{
    chunk_t *c = (chunk_t *)arg;
    primesieve_iterator it;
    primesieve_init(&it);
    primesieve_jump_to(&it, c->start, c->stop);
    c->count = 0;
    c->sum.lo = 0; c->sum.hi = 0;
    for (;;) {
        uint64_t p = primesieve_next_prime(&it);
        if (p >= c->stop) break;
        c->count++;
        u192_add_u128(&c->sum, (u128)p * p);
    }
    primesieve_free_iterator(&it);
    return NULL;
}

static void *pass2(void *arg)
{
    chunk_t *c = (chunk_t *)arg;
    primesieve_iterator it;
    primesieve_init(&it);
    primesieve_jump_to(&it, c->start, c->stop);
    uint64_t n = c->n0;
    u192 s = c->s0;
    for (;;) {
        uint64_t p = primesieve_next_prime(&it);
        if (p >= c->stop) break;
        n++;
        u192_add_u128(&s, (u128)p * p);
        uint64_t m = n % 6;
        if ((m == 1 || m == 5) && u192_mod(&s, n) == 0) {
            printf("HIT  n = %" PRIu64 "  (p_n = %" PRIu64 ")\n", n, p);
            fflush(stdout);
        }
    }
    primesieve_free_iterator(&it);
    return NULL;
}

static void save_ckpt(const char *path, uint64_t next_lo, uint64_t n, const u192 *s)
{
    char tmp[1024];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { perror("ckpt"); return; }
    fprintf(f, "%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
            next_lo, n, s->hi,
            (uint64_t)(s->lo >> 64), (uint64_t)s->lo);
    fclose(f);
    rename(tmp, path);
}

static int load_ckpt(const char *path, uint64_t *next_lo, uint64_t *n, u192 *s)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    uint64_t hi, lohi, lolo;
    if (fscanf(f, "%" SCNu64 " %" SCNu64 " %" SCNu64 " %" SCNu64 " %" SCNu64,
               next_lo, n, &hi, &lohi, &lolo) != 5) {
        fclose(f);
        fprintf(stderr, "bad checkpoint file %s\n", path);
        exit(1);
    }
    fclose(f);
    s->hi = hi;
    s->lo = ((u128)lohi << 64) | lolo;
    return 1;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s <stop_prime_value> [ckpt_file] [threads] [window_width]\n", argv[0]);
        return 1;
    }
    uint64_t stop_value = strtoull(argv[1], NULL, 10);
    const char *ckpt = (argc > 2) ? argv[2] : "a111441.ckpt";
    int nthreads = (argc > 3) ? atoi(argv[3]) : 12;
    uint64_t window = (argc > 4) ? strtoull(argv[4], NULL, 10) : 120000000000ULL;

    uint64_t lo = 0, n = 0;
    u192 s = {0, 0};
    if (load_ckpt(ckpt, &lo, &n, &s))
        fprintf(stderr, "resumed: value=%" PRIu64 " n=%" PRIu64 "\n", lo, n);

    chunk_t *chunks = calloc((size_t)nthreads, sizeof *chunks);
    pthread_t *tids = calloc((size_t)nthreads, sizeof *tids);
    time_t t0 = time(NULL);
    uint64_t n_start = n;

    while (lo < stop_value) {
        uint64_t hi = lo + window;
        if (hi > stop_value) hi = stop_value;
        uint64_t width = (hi - lo) / (uint64_t)nthreads;

        for (int i = 0; i < nthreads; i++) {
            chunks[i].start = lo + width * (uint64_t)i;
            chunks[i].stop  = (i == nthreads - 1) ? hi : lo + width * (uint64_t)(i + 1);
        }
        for (int i = 0; i < nthreads; i++)
            pthread_create(&tids[i], NULL, pass1, &chunks[i]);
        for (int i = 0; i < nthreads; i++)
            pthread_join(tids[i], NULL);

        /* prefix sums -> global state at each subrange start */
        uint64_t np = n;
        u192 sp = s;
        for (int i = 0; i < nthreads; i++) {
            chunks[i].n0 = np;
            chunks[i].s0 = sp;
            np += chunks[i].count;
            u192_add(&sp, &chunks[i].sum);
        }

        for (int i = 0; i < nthreads; i++)
            pthread_create(&tids[i], NULL, pass2, &chunks[i]);
        for (int i = 0; i < nthreads; i++)
            pthread_join(tids[i], NULL);

        n = np;
        s = sp;
        lo = hi;
        save_ckpt(ckpt, lo, n, &s);

        long el = (long)(time(NULL) - t0);
        double rate = el > 0 ? (double)(n - n_start) / (double)el : 0;
        fprintf(stderr, "[%lds] value <= %" PRIu64 "  n = %" PRIu64 "  (%.0fM primes/s)\n",
                el, lo, n, rate / 1e6);
    }

    free(chunks);
    free(tids);
    return 0;
}
