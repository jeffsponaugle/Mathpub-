/*
 * bbp.c — multithreaded BBP hex digit extraction for pi.
 *
 *   pi = sum_{k>=0} 16^-k * ( 4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6) )
 *
 * Computes frac(16^n * pi) as an exact 128-bit fixed-point value: each
 * term is floor(2^128 * (16^(n-k) mod d) / d), accumulated mod 2^128
 * (natural unsigned wraparound). Accumulated floor error is < 8(n+32)
 * ulps of 2^-128, i.e. < 2^-83 even at n = 10^12, so the top digits are
 * reliable (output is capped at 16 hex digits = the top 64 bits).
 *
 * Modular exponentiation uses Montgomery multiplication. The even
 * moduli are reduced to odd ones first:
 *     16^e mod (8k+4) = 4 * (2^(4e-2) mod (2k+1))     (e >= 1)
 *     16^e mod (8k+6) = 2 * (2^(4e-1) mod (4k+3))     (e >= 1)
 *
 * Terms are independent, so threads pull chunks of the k-range from an
 * atomic counter — scaling is essentially linear in cores.
 *
 * This file is the engine only; bbp_main.c wraps it as a standalone CLI
 * and piverify.c uses it for verification.
 */

#include <inttypes.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "bbp.h"

/* inverse of odd u mod 2^64 (Newton doubling) */
static u64 inv64(u64 u) {
    u64 x = u;
    for (int i = 0; i < 6; i++)
        x *= 2 - u * x;
    return x;
}

static inline u64 mont_mul(u64 a, u64 b, u64 u, u64 nprime) {
    u128 t = (u128)a * b;
    u64 m = (u64)t * nprime;
    u128 s = t + (u128)m * u;
    u64 r = (u64)(s >> 64);
    return r >= u ? r - u : r;
}

/* 2^e mod u for odd u */
static u64 pow2_mod(u64 e, u64 u) {
    if (u == 1)
        return 0;
    u64 nprime = 0 - inv64(u);
    u64 r = (u64)((((u128)1) << 64) % u);      /* Montgomery form of 1 */
    u64 base = (u64)((((u128)2) << 64) % u);   /* Montgomery form of 2 */
    while (e) {
        if (e & 1)
            r = mont_mul(r, base, u, nprime);
        base = mont_mul(base, base, u, nprime);
        e >>= 1;
    }
    return mont_mul(r, 1, u, nprime);
}

/* floor(2^128 * r / d) for r < d */
static inline u128 fixdiv(u64 r, u64 d) {
    u128 a = (u128)r << 64;
    u64 qhi = (u64)(a / d);
    u128 b = (u128)(u64)(a % d) << 64;
    u64 qlo = (u64)(b / d);
    return ((u128)qhi << 64) | qlo;
}

/* contribution of series index k to frac(16^n * pi), mod 2^128 */
static u128 term_sum(u64 n, u64 k) {
    u64 e = n - k;
    u128 acc = 0;

    u64 d1 = 8 * k + 1;
    acc += 4 * fixdiv(pow2_mod(4 * e, d1), d1);

    u64 d4 = 8 * k + 4;
    u64 r4 = (e == 0) ? 1 : 4 * pow2_mod(4 * e - 2, 2 * k + 1);
    acc -= 2 * fixdiv(r4, d4);

    u64 d5 = 8 * k + 5;
    acc -= fixdiv(pow2_mod(4 * e, d5), d5);

    u64 d6 = 8 * k + 6;
    u64 r6 = (e == 0) ? 1 : 2 * pow2_mod(4 * e - 1, 4 * k + 3);
    acc -= fixdiv(r6, d6);

    return acc;
}

/* terms with k > n: 16^(n-k) is a plain right shift */
static u128 tail_sum(u64 n) {
    u128 acc = 0;
    for (u64 i = 1; 4 * i < 128; i++) {   /* sh=128 would be UB; term is 0 */
        u64 k = n + i;
        int sh = (int)(4 * i);
        static const int js[4] = {1, 4, 5, 6};
        static const int cf[4] = {4, -2, -1, -1};
        for (int t = 0; t < 4; t++) {
            u64 d = 8 * k + (u64)js[t];
            u128 x = (u128)d << sh;
            if ((u64)(x >> sh) != d)
                continue;                    /* d*16^i >= 2^128: term is 0 */
            u128 term = (~(u128)0) / x;      /* == floor(2^128 / x) here */
            acc += (cf[t] > 0) ? (u128)cf[t] * term
                               : (u128)0 - (u128)(-cf[t]) * term;
        }
    }
    return acc;
}

#define CHUNK 16384

typedef struct {
    u64 n;
    _Atomic(u64) *next;
    u128 *out;
    pthread_mutex_t *mu;
} ctx_t;

static void *worker(void *arg) {
    ctx_t *c = (ctx_t *)arg;
    u128 acc = 0;
    for (;;) {
        u64 k0 = atomic_fetch_add(c->next, (u64)CHUNK);
        if (k0 > c->n)
            break;
        u64 k1 = k0 + CHUNK;
        if (k1 > c->n + 1)
            k1 = c->n + 1;
        for (u64 k = k0; k < k1; k++)
            acc += term_sum(c->n, k);
    }
    pthread_mutex_lock(c->mu);
    *c->out += acc;
    pthread_mutex_unlock(c->mu);
    return NULL;
}

typedef struct {
    _Atomic(u64) *next;
    _Atomic(int) *done;
    u64 n;
} mon_t;

static void *monitor(void *arg) {
    mon_t *m = (mon_t *)arg;
    struct timespec ts = {2, 0};
    while (!atomic_load(m->done)) {
        u64 at = atomic_load(m->next);
        if (at > m->n)
            at = m->n;
        fprintf(stderr, "\r  bbp: %5.1f%% of %" PRIu64 " terms",
                m->n ? 100.0 * (double)at / (double)m->n : 100.0, m->n + 1);
        nanosleep(&ts, NULL);
    }
    fprintf(stderr, "\r%40s\r", "");
    return NULL;
}

u128 bbp_frac(u64 n, int threads, int verbose) {
    u128 total = 0;
    _Atomic(u64) next = 0;
    _Atomic(int) done = 0;
    pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
    ctx_t ctx = {n, &next, &total, &mu};

    if (threads < 1)
        threads = 1;
    pthread_t *tid = malloc((size_t)threads * sizeof(pthread_t));
    pthread_t mon;
    mon_t m = {&next, &done, n};
    if (verbose)
        pthread_create(&mon, NULL, monitor, &m);
    for (int i = 0; i < threads; i++)
        pthread_create(&tid[i], NULL, worker, &ctx);
    for (int i = 0; i < threads; i++)
        pthread_join(tid[i], NULL);
    atomic_store(&done, 1);
    if (verbose)
        pthread_join(mon, NULL);
    free(tid);

    return total + tail_sum(n);
}

void bbp_hex_digits(u64 pos0, int t, int threads, int verbose, char out[17]) {
    if (t < 1)
        t = 1;
    if (t > 16)
        t = 16;
    u128 f = bbp_frac(pos0, threads, verbose);
    snprintf(out, 17, "%016" PRIx64, (u64)(f >> 64));
    out[t] = '\0';
}
