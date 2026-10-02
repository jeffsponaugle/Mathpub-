/*
 * a005115.c
 *
 * Compute and extend OEIS A005115:
 *
 *   "Let i, i+d, i+2d, ..., i+(n-1)d be an n-term arithmetic progression of
 *    primes; choose the one which minimizes the last term; then a(n) = last
 *    term i+(n-1)d."
 *
 * Known terms (n = 1..26):
 *   2, 3, 7, 23, 29, 157, 907, 1669, 1879, 2089, 249037, 262897, 725663,
 *   36850999, 173471351, 198793279, 4827507229, 17010526363, 83547839407,
 *   572945039351, 6269243827111, 35742689530423, 449924511422857,
 *   1217585417914253, 5773236905395679, 12783396861134173
 * a(23)..a(26) were proved minimal by Joris Perrenet and Dmitry Petukhov
 * (2023-2026).  a(27) is open: the smallest known AP-27 (Rob Gahan /
 * PrimeGrid, 23 Sep 2019) is 224584605939537911 + 81292139*23#*j, j = 0..26,
 * so a(27) <= 696112717486210091.
 *
 * This file is the CPU reference implementation (any n = k >= 3).  It is
 * exact but far too slow for k = 27; the production search is the CUDA port
 * in cuda/, which runs the same algorithm and is validated against this file.
 *
 * Structure of an AP-k of primes
 * ------------------------------
 * Write the AP as a, a+D, ..., a+(k-1)D and L = a+(k-1)D for its last term.
 * If a > k then every prime p <= k divides D (otherwise the k terms run
 * through all residues mod p and one of them, >= a > p, is divisible by p).
 * So D = K*P with P = k# (product of the primes <= k).  The only other shape
 * is a = k with k prime and D = K'*(k-1)#, k not dividing K' (e.g. 7 + 150j,
 * a(7) = 907); a < k is impossible for k >= 3.
 *
 * For a prime q > k not dividing D the AP avoids q iff a mod q is not one
 * of the k residues -iD (i = 0..k-1), i.e. a/D mod q lies in [1, q-k].  For
 * q | D the condition is just a mod q != 0.  A term can equal a sieving
 * prime only if a itself is small, so the main sieve covers a > SMALL_A
 * (4096) and a separate direct pass covers a <= SMALL_A.
 *
 * The sieve (Jaroslaw Wroblewski's AP26 method, reorganised by last term)
 * ----------------------------------------------------------------------
 * Fix K.  Stage 1: pick a modulus MOD = product of a set Q1 of primes (2, 3,
 * 5 and the primes just above k with q not dividing K, chosen greedily by
 * ln(good/q)/ln q under MOD <= MODMAX).  The residues a mod MOD that are good
 * for every q in Q1 form a product set; with a CRT idempotent e_q per prime
 * they are R = sum_q (1+d_q) s_q (mod MOD), s_q = (D mod q) e_q, digits
 * 0 <= d_q < q-k (q-1 if q | D), enumerated by an odometer using only
 * additions.  For k = 27 and a normal K, Q1 = {2,3,5,29,...,59}, MOD =
 * 258559632607830 and there are 2,385,510,400 good residues.
 *
 * Stage 2 (Wroblewski's second trick): one residue R stands for the 64
 * numbers a = R + b*MOD, b = 64w..64w+63 (word w).  Whether a + iD avoids a
 * prime q depends only on R mod q, so a precomputed 64-bit mask per (q, R mod
 * q) tells which of the 64 are still alive.  AND-ing the masks for the next
 * N2 primes (61, 67, 71, ... for k = 27) kills a word after about 12 lookups
 * on average.  The residues R mod q for the first NTRACK primes are carried
 * along the odometer with additions; later primes use R % q when needed.
 *
 * Stage 3 checks the few surviving a against every other prime <= 4096
 * (including the primes <= k outside Q1, which only need a mod p != 0), and
 * stage 4 runs a base-2 strong probable prime test on all k terms followed
 * by a deterministic Miller-Rabin (bases 2..37) on each term.
 *
 * Order of the search
 * -------------------
 * The pairs (a, K) with L <= X form a triangle.  A work item is (K, w): it
 * covers a in [64w*MOD, 64(w+1)*MOD), and its smallest possible last term
 * is L_lo = (k-1)KP + 64w*MOD.  Items are processed in levels of L_lo
 * ([l*DELTA, (l+1)*DELTA)).  Once levels 0..l are complete, every AP-k with
 * L < (l+1)*DELTA has been seen, so the first AP found, once its level is
 * finished, is a(k) with a certificate.  Finding an AP lowers the cap X so
 * later items are skipped.
 *
 * Usage
 * -----
 *   a005115 selftest [KMAX]          reproduce a(3..KMAX) (default 20) and
 *                                    check them against the OEIS data and the
 *                                    progressions on Luhn's record page
 *   a005115 search k X [options]     find a(k) if a(k) <= X (certified)
 *   a005115 all k X [options]        list every AP-k with last term <= X
 *   a005115 range k X K0 K1          every AP-k with L <= X and K in [K0,K1]
 *                                    (for cross-checking the GPU tool)
 *   a005115 bench [k X K W SECS]     time the stage 1+2 kernel on items (K, W)
 *   a005115 verify a d k             check a, a+d, ..., a+(k-1)d
 * Options: -t T threads, -m MODMAX, -l DELTA level width, -n N2 stage-2
 * primes, -S FILE checkpoint (resume with the same command), -i SECS
 * checkpoint interval, -q quiet (no status line).
 */
#define _GNU_SOURCE
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <signal.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef unsigned __int128 u128;

/* ------------------------------------------------------------------ */
/* Parameters                                                          */
/* ------------------------------------------------------------------ */

#define SMALL_A     4096u  /* main sieve covers a > SMALL_A; all sieving primes <= SMALL_A */
#define TRIAL_MAX   4096u  /* stage 3 trial primes (must be <= SMALL_A) */
#define N1MAX       24
#define N2MAX       64
#define NTRACK      12     /* stage-2 primes whose residues ride along the odometer */
#define Q1_CAND_MAX 211u   /* candidate primes for Q1 */

static int      KLEN;        /* k, the AP length */
static uint64_t PRIM;        /* k# */
static uint64_t MODMAX;      /* bound for the stage-1 modulus */
static int      N2 = 40;     /* stage-2 primes */
static int      nthreads = 0;
static int      quiet = 0;

/* ------------------------------------------------------------------ */
/* Small primes                                                        */
/* ------------------------------------------------------------------ */

#define SP_LIMIT 65536u
static uint32_t smallp[6600];
static int nsmallp;

static void init_small_primes(void)
{
    static uint8_t comp[SP_LIMIT + 1];
    for (uint32_t i = 2; i <= SP_LIMIT; i++) {
        if (comp[i]) continue;
        smallp[nsmallp++] = i;
        for (uint64_t j = (uint64_t)i * i; j <= SP_LIMIT; j += i) comp[j] = 1;
    }
}

static uint64_t primorial(int k)
{
    uint64_t r = 1;
    for (int i = 0; smallp[i] <= (uint32_t)k; i++) r *= smallp[i];
    return r;
}

/* inverse of a mod m (gcd(a, m) = 1, m < 2^31) */
static uint32_t inv_mod32(uint32_t a, uint32_t m)
{
    int64_t t = 0, nt = 1, r = m, nr = a % m;
    while (nr) {
        int64_t q = r / nr, tmp;
        tmp = t - q * nt; t = nt; nt = tmp;
        tmp = r - q * nr; r = nr; nr = tmp;
    }
    if (t < 0) t += m;
    return (uint32_t)t;
}

/* ------------------------------------------------------------------ */
/* Primality: Montgomery strong probable prime test                    */
/* ------------------------------------------------------------------ */

static inline uint64_t mont_ninv(uint64_t n)
{
    uint64_t x = n;                       /* n*x = 1 mod 2^3 for odd n */
    for (int i = 0; i < 5; i++) x *= 2 - n * x;
    return x;
}

static inline uint64_t mont_mul(uint64_t a, uint64_t b, uint64_t n, uint64_t ninv)
{
    u128 t = (u128)a * b;
    uint64_t tl = (uint64_t)t, th = (uint64_t)(t >> 64);
    uint64_t m = tl * ninv;
    uint64_t mh = (uint64_t)(((u128)m * n) >> 64);
    uint64_t r = th - mh;
    return th < mh ? r + n : r;
}

/* strong probable prime to base 'base'; n odd, n > 2 */
static int sprp(uint64_t n, uint64_t base)
{
    uint64_t d = n - 1;
    int s = __builtin_ctzll(d);
    d >>= s;
    uint64_t ninv = mont_ninv(n);
    uint64_t one = (uint64_t)((((u128)1) << 64) % n), mone = n - one;
    uint64_t b = (uint64_t)((((u128)(base % n)) << 64) % n);
    if (b == 0) return 1;
    uint64_t x = one;
    for (int i = 63 - __builtin_clzll(d); i >= 0; i--) {
        x = mont_mul(x, x, n, ninv);
        if ((d >> i) & 1) x = mont_mul(x, b, n, ninv);
    }
    if (x == one || x == mone) return 1;
    for (int i = 1; i < s; i++) {
        x = mont_mul(x, x, n, ninv);
        if (x == mone) return 1;
        if (x == one) return 0;
    }
    return 0;
}

/* deterministic for n < 3.3e24 (bases 2..37) */
static int is_prime_u64(uint64_t n)
{
    static const uint32_t bases[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return 0;
    for (int i = 0; i < 12; i++) {
        if (n == bases[i]) return 1;
        if (n % bases[i] == 0) return 0;
    }
    if (n < 41 * 41) return 1;
    for (int i = 0; i < 12; i++)
        if (!sprp(n, bases[i])) return 0;
    return 1;
}

static _Atomic uint64_t n_sprp_calls;

/* all of a, a+d, ..., a+(k-1)d prime? (odd a, even d; base-2 screen first) */
static int ap_is_prime(uint64_t a, uint64_t d, int k)
{
    for (int i = 0; i < k; i++) {
        uint64_t n = a + (uint64_t)i * d;
        atomic_fetch_add_explicit(&n_sprp_calls, 1, memory_order_relaxed);
        if (n < 3 || !(n & 1)) { if (n != 2) return 0; continue; }
        if (!sprp(n, 2)) return 0;
    }
    for (int i = 0; i < k; i++)
        if (!is_prime_u64(a + (uint64_t)i * d)) return 0;
    return 1;
}

/* ------------------------------------------------------------------ */
/* Results                                                             */
/* ------------------------------------------------------------------ */

typedef struct { uint64_t a, d, L; } AP;

static pthread_mutex_t res_lock = PTHREAD_MUTEX_INITIALIZER;
static AP *found;
static size_t nfound, capfound;
static _Atomic uint64_t Xcap;          /* current cap on L */
static int search_min_mode;            /* lower Xcap to each AP found */
static FILE *results_fp;               /* optional copy of AP lines */
static int print_aps = 1;              /* print each AP as it is found */

static void report_ap(uint64_t a, uint64_t d)
{
    uint64_t L = a + (uint64_t)(KLEN - 1) * d;
    pthread_mutex_lock(&res_lock);
    for (size_t i = 0; i < nfound; i++)
        if (found[i].a == a && found[i].d == d) { pthread_mutex_unlock(&res_lock); return; }
    if (nfound == capfound) {
        capfound = capfound ? 2 * capfound : 1024;
        found = realloc(found, capfound * sizeof(AP));
    }
    found[nfound++] = (AP){a, d, L};
    if (!print_aps) ;
    else if (d % PRIM == 0)
        printf("AP%d a=%" PRIu64 " d=%" PRIu64 "*%d# L=%" PRIu64 "\n", KLEN, a, d / PRIM, KLEN, L);
    else
        printf("AP%d a=%" PRIu64 " d=%" PRIu64 " L=%" PRIu64 "\n", KLEN, a, d, L);
    fflush(stdout);
    if (results_fp) {
        fprintf(results_fp, "%d %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", KLEN, a, d, L);
        fflush(results_fp);
    }
    if (search_min_mode) {
        uint64_t cur = atomic_load(&Xcap);
        while (L < cur && !atomic_compare_exchange_weak(&Xcap, &cur, L)) {}
    }
    pthread_mutex_unlock(&res_lock);
}

static int cmp_ap(const void *x, const void *y)
{
    const AP *p = x, *q = y;
    if (p->L != q->L) return p->L < q->L ? -1 : 1;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Per-K context                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    uint64_t K, D, MOD, n0, nR;
    int n1;
    uint32_t q1[N1MAX], g1[N1MAX];
    uint64_t s1[N1MAX];              /* digit steps mod MOD */
    uint64_t c1[N1MAX];              /* odometer carry increments mod MOD */
    int n2, ntr;
    uint32_t q2[N2MAX], modq2[N2MAX], n0q2[N2MAX];
    uint32_t inc[N1MAX][NTRACK], incw[N1MAX][NTRACK];
    uint32_t toff[N2MAX + 1];        /* table offsets */
    uint8_t *bad;                    /* bad residues, per stage-2 prime at toff */
    uint64_t *tab;                   /* masks for the current word */
    int n3;
    uint32_t *q3, *mdinv3;           /* mdinv3 = -D^-1 mod q, or 0 when q | D */
    uint8_t *qd3;
} Kctx;

typedef struct {
    Kctx c;
    uint64_t words, cand2, cand3, aps;
    void *eng;                       /* GPU engine state (metal/, cuda/) */
    int dbg;                         /* record stage-2 survivors instead of testing them */
    uint64_t *dbg_a;
    size_t dbg_n, dbg_cap;
} Worker;

static void dbg_push(Worker *W, uint64_t a)
{
    if (W->dbg_n == W->dbg_cap) {
        W->dbg_cap = W->dbg_cap ? 2 * W->dbg_cap : 4096;
        W->dbg_a = realloc(W->dbg_a, W->dbg_cap * sizeof(uint64_t));
    }
    W->dbg_a[W->dbg_n++] = a;
}

static void kctx_alloc(Kctx *c)
{
    c->bad = malloc(70000);
    c->tab = malloc(70000 * sizeof(uint64_t));
    c->q3 = malloc(sizeof(uint32_t) * 700);
    c->mdinv3 = malloc(sizeof(uint32_t) * 700);
    c->qd3 = malloc(700);
}

/* choose Q1 for this K (greedy by ln(g/q)/ln q), returns MOD */
static uint64_t choose_q1(uint64_t K, uint32_t *q1, uint32_t *g1, int *n1)
{
    uint32_t cq[64], cg[64];
    double ce[64];
    int nc = 0;
    for (int i = 0; smallp[i] <= Q1_CAND_MAX; i++) {
        uint32_t q = smallp[i];
        uint32_t g = (q <= (uint32_t)KLEN || K % q == 0) ? q - 1 : q - KLEN;
        cq[nc] = q; cg[nc] = g;
        ce[nc] = log((double)g / q) / log((double)q);
        nc++;
    }
    /* insertion sort by efficiency (most negative first), ties by q */
    for (int i = 1; i < nc; i++) {
        uint32_t tq = cq[i], tg = cg[i]; double te = ce[i];
        int j = i - 1;
        while (j >= 0 && (ce[j] > te || (ce[j] == te && cq[j] > tq))) {
            cq[j + 1] = cq[j]; cg[j + 1] = cg[j]; ce[j + 1] = ce[j]; j--;
        }
        cq[j + 1] = tq; cg[j + 1] = tg; ce[j + 1] = te;
    }
    uint64_t MOD = 1;
    int n = 0;
    for (int i = 0; i < nc && n < N1MAX; i++) {
        if ((u128)MOD * cq[i] > MODMAX) continue;
        MOD *= cq[i];
        q1[n] = cq[i]; g1[n] = cg[i]; n++;
    }
    /* innermost digit = largest radix: sort ascending by g */
    for (int i = 1; i < n; i++) {
        uint32_t tq = q1[i], tg = g1[i];
        int j = i - 1;
        while (j >= 0 && g1[j] > tg) { q1[j + 1] = q1[j]; g1[j + 1] = g1[j]; j--; }
        q1[j + 1] = tq; g1[j + 1] = tg;
    }
    *n1 = n;
    return MOD;
}

static uint64_t mod_for_K(uint64_t K)
{
    uint32_t q1[N1MAX], g1[N1MAX];
    int n1;
    return choose_q1(K, q1, g1, &n1);
}

static void kctx_build(Kctx *c, uint64_t K)
{
    c->K = K;
    c->D = K * PRIM;
    c->MOD = choose_q1(K, c->q1, c->g1, &c->n1);
    uint64_t MOD = c->MOD, D = c->D;
    uint64_t s[N1MAX];
    c->nR = 1;
    uint64_t n0 = 0;
    for (int j = 0; j < c->n1; j++) {
        uint32_t q = c->q1[j];
        uint64_t Mq = MOD / q;
        uint64_t e = (uint64_t)((u128)Mq * inv_mod32((uint32_t)(Mq % q), q) % MOD);
        uint32_t delta = (D % q == 0) ? 1 : (uint32_t)(D % q);
        s[j] = (uint64_t)((u128)delta * e % MOD);
        c->s1[j] = s[j];
        n0 = (n0 + s[j]) % MOD;
        c->nR *= c->g1[j];
    }
    c->n0 = n0;
    for (int j = 0; j < c->n1; j++) {
        u128 v = s[j];
        for (int i = j + 1; i < c->n1; i++)
            v += (u128)(MOD - s[i]) * (c->g1[i] - 1);
        c->c1[j] = (uint64_t)(v % MOD);
    }

    /* stage-2 primes: the first N2 primes > k outside Q1 */
    uint8_t inq1[Q1_CAND_MAX + 1];
    memset(inq1, 0, sizeof inq1);
    for (int j = 0; j < c->n1; j++) inq1[c->q1[j]] = 1;
    c->n2 = 0;
    uint32_t off = 0;
    for (int i = 0; i < nsmallp && c->n2 < N2; i++) {
        uint32_t q = smallp[i];
        if (q <= (uint32_t)KLEN) continue;
        if (q <= Q1_CAND_MAX && inq1[q]) continue;
        int t = c->n2++;
        c->q2[t] = q;
        c->toff[t] = off;
        uint8_t *bad = c->bad + off;
        memset(bad, 0, q);
        uint32_t dq = (uint32_t)(D % q);
        for (int ii = 0; ii < KLEN; ii++)
            bad[(q - (uint32_t)((uint64_t)ii * dq % q)) % q] = 1;
        c->modq2[t] = (uint32_t)(MOD % q);
        c->n0q2[t] = (uint32_t)(c->n0 % q);
        off += q;
    }
    c->toff[c->n2] = off;
    c->ntr = c->n2 < NTRACK ? c->n2 : NTRACK;
    for (int j = 0; j < c->n1; j++)
        for (int t = 0; t < c->ntr; t++) {
            uint32_t q = c->q2[t];
            uint32_t v = (uint32_t)(c->c1[j] % q);
            c->inc[j][t] = v;
            c->incw[j][t] = (v + q - c->modq2[t]) % q;
        }

    /* stage-3 primes: everything else up to TRIAL_MAX, strongest first */
    uint8_t used[TRIAL_MAX + 1];
    memset(used, 0, sizeof used);
    for (int j = 0; j < c->n1; j++) used[c->q1[j]] = 1;
    for (int t = 0; t < c->n2; t++) if (c->q2[t] <= TRIAL_MAX) used[c->q2[t]] = 1;
    int n3 = 0;
    for (int i = 0; smallp[i] <= TRIAL_MAX; i++) {
        uint32_t q = smallp[i];
        if (used[q]) continue;
        uint32_t dq = (uint32_t)(D % q);
        c->q3[n3] = q;
        if (dq == 0) { c->qd3[n3] = 1; c->mdinv3[n3] = 0; }
        else { c->qd3[n3] = 0; c->mdinv3[n3] = (q - inv_mod32(dq, q)) % q; }
        n3++;
    }
    /* order: q not dividing D ascending (rejects k/q), then q | D (rejects 1/q) */
    {
        uint32_t tq[700], tm[700]; uint8_t td[700];
        int m = 0;
        for (int pass = 0; pass < 2; pass++)
            for (int i = 0; i < n3; i++)
                if (c->qd3[i] == pass) { tq[m] = c->q3[i]; tm[m] = c->mdinv3[i]; td[m] = c->qd3[i]; m++; }
        memcpy(c->q3, tq, n3 * sizeof(uint32_t));
        memcpy(c->mdinv3, tm, n3 * sizeof(uint32_t));
        memcpy(c->qd3, td, n3);
    }
    c->n3 = n3;
}

/* masks for word w: bit j of tab[toff[t] + r] = a = R + (64w+j)*MOD good mod q, r = R mod q */
static void build_tables(Kctx *c, uint64_t w)
{
    for (int t = 0; t < c->n2; t++) {
        uint32_t q = c->q2[t], m = c->modq2[t];
        uint32_t o = (uint32_t)((uint64_t)((64 * w) % q) * m % q);
        const uint8_t *bad = c->bad + c->toff[t];
        uint64_t *T = c->tab + c->toff[t];
        for (uint32_t r = 0; r < q; r++) {
            uint32_t x = r + o; if (x >= q) x -= q;
            uint64_t bits = 0;
            for (int j = 0; j < 64; j++) {
                if (!bad[x]) bits |= 1ULL << j;
                x += m; if (x >= q) x -= q;
            }
            T[r] = bits;
        }
    }
}

static inline int stage3_ok(const Kctx *c, uint64_t a)
{
    for (int t = 0; t < c->n3; t++) {
        uint32_t q = c->q3[t];
        uint32_t r = (uint32_t)(a % q);
        if (c->qd3[t]) { if (r == 0) return 0; }
        else if ((uint64_t)r * c->mdinv3[t] % q < (uint64_t)KLEN) return 0;
    }
    return 1;
}

/* stages 3 and 4 for one candidate first term a (SMALL_A < a <= A_hi) */
static void test_candidate(Worker *W, uint64_t a)
{
    Kctx *c = &W->c;
    if (!stage3_ok(c, a)) return;
    W->cand3++;
    if (ap_is_prime(a, c->D, KLEN)) { W->aps++; report_ap(a, c->D); }
}

static void handle_mask(Worker *W, uint64_t R, uint64_t w, uint64_t mask, uint64_t A_hi)
{
    Kctx *c = &W->c;
    while (mask) {
        int j = __builtin_ctzll(mask);
        mask &= mask - 1;
        uint64_t a = R + (64 * w + (uint64_t)j) * c->MOD;
        if (W->dbg) { dbg_push(W, a); continue; }
        W->cand2++;
        if (a <= SMALL_A || a > A_hi) continue;
        test_candidate(W, a);
    }
}

/* one work item: all good R, word w; APs with a <= A_hi */
static void process_item(Worker *W, uint64_t w, uint64_t A_hi)
{
    Kctx *c = &W->c;
    build_tables(c, w);
    const int n1 = c->n1, n2 = c->n2, ntr = c->ntr;
    const uint64_t MOD = c->MOD;
    const uint64_t *tab = c->tab;
    uint32_t d[N1MAX] = {0};
    uint32_t rr[NTRACK];
    const uint64_t *T[N2MAX];
    for (int t = 0; t < n2; t++) T[t] = tab + c->toff[t];
    for (int t = 0; t < ntr; t++) rr[t] = c->n0q2[t];
    uint64_t R = c->n0;
    for (;;) {
        uint64_t m = ~0ULL;
        int t = 0;
        for (; t < ntr; t++) { m &= T[t][rr[t]]; if (!m) goto next; }
        for (; t < n2; t++) { m &= T[t][R % c->q2[t]]; if (!m) goto next; }
        handle_mask(W, R, w, m, A_hi);
    next:;
        int j = n1 - 1;
        while (j >= 0 && d[j] == c->g1[j] - 1) { d[j] = 0; j--; }
        if (j < 0) break;
        d[j]++;
        R += c->c1[j];
        int wrap = R >= MOD;
        if (wrap) R -= MOD;
        const uint32_t *inc = wrap ? c->incw[j] : c->inc[j];
        for (t = 0; t < ntr; t++) {
            uint32_t v = rr[t] + inc[t];
            rr[t] = v >= c->q2[t] ? v - c->q2[t] : v;
        }
    }
    W->words += c->nR;
}

/* item engine (CPU by default; the GPU ports install their own) */
static void (*item_engine)(Worker *W, uint64_t w, uint64_t A_hi) = process_item;
static void (*worker_init_hook)(Worker *W) = NULL;

/* ------------------------------------------------------------------ */
/* Small first terms: a <= SMALL_A, tested directly                     */
/* ------------------------------------------------------------------ */

typedef struct {
    uint64_t Xc, K0, K1;
    _Atomic int next;
    int na;
    uint32_t *as;     /* first terms; case B (a = k prime) flagged by as == KLEN */
} SmallJob;

static void small_one(uint32_t a, uint64_t Xc, uint64_t K0, uint64_t K1)
{
    int caseB = (a == (uint32_t)KLEN);
    uint64_t unit = caseB ? PRIM / KLEN : PRIM;
    if (Xc < a) return;
    uint64_t Kmax = (Xc - a) / ((uint64_t)(KLEN - 1) * unit);
    if (!caseB && Kmax > K1) Kmax = K1;
    uint64_t Kmin = caseB ? 1 : K0;
    if (Kmax < Kmin) return;
    /* bad K residues for the first primes r > k, r != a */
    enum { NR = 24 };
    uint32_t r[NR];
    uint8_t *badk[NR];
    int nr = 0;
    for (int i = 0; i < nsmallp && nr < NR; i++) {
        uint32_t q = smallp[i];
        if (q <= (uint32_t)KLEN || q == a) continue;
        r[nr] = q;
        badk[nr] = calloc(q, 1);
        uint32_t um = (uint32_t)(unit % q), am = a % q;
        for (int ii = 1; ii < KLEN; ii++) {
            /* a + ii*K*unit = 0 mod q  <=>  K = -a / (ii*unit) */
            uint32_t den = (uint32_t)((uint64_t)ii * um % q);
            uint32_t x = (uint32_t)((uint64_t)(q - am) % q * inv_mod32(den, q) % q);
            badk[nr][x] = 1;
        }
        nr++;
    }
    uint32_t rmax = r[nr - 1];
    uint32_t k0 = (uint32_t)(Kmin % r[0]);
    for (uint64_t K = Kmin; K <= Kmax; K++, k0 = (k0 + 1 == r[0]) ? 0 : k0 + 1) {
        if (caseB && K % KLEN == 0) continue;
        /* the sieve is only valid once every term after a exceeds the sieving primes */
        if ((uint64_t)a + K * unit > rmax) {
            if (badk[0][k0]) continue;
            int ok = 1;
            for (int t = 1; t < nr; t++)
                if (badk[t][K % r[t]]) { ok = 0; break; }
            if (!ok) continue;
        }
        uint64_t d = K * unit;
        if ((uint64_t)a + (uint64_t)(KLEN - 1) * d > atomic_load(&Xcap)) break;
        if (ap_is_prime(a, d, KLEN)) report_ap(a, d);
    }
    for (int t = 0; t < nr; t++) free(badk[t]);
}

static void *small_thread(void *arg)
{
    SmallJob *J = arg;
    for (;;) {
        int i = atomic_fetch_add(&J->next, 1);
        if (i >= J->na) break;
        small_one(J->as[i], J->Xc, J->K0, J->K1);
    }
    return NULL;
}

static void small_pass(uint64_t Xc, uint64_t K0, uint64_t K1)
{
    SmallJob J;
    J.Xc = Xc; J.K0 = K0; J.K1 = K1;
    atomic_init(&J.next, 0);
    J.as = malloc(sizeof(uint32_t) * 700);
    J.na = 0;
    for (int i = 0; smallp[i] <= SMALL_A; i++) {
        uint32_t p = smallp[i];
        if (p > (uint32_t)KLEN || (p == (uint32_t)KLEN && K0 == 1)) J.as[J.na++] = p;
    }
    pthread_t th[256];
    for (int i = 0; i < nthreads; i++) pthread_create(&th[i], NULL, small_thread, &J);
    for (int i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
    free(J.as);
}

/* ------------------------------------------------------------------ */
/* Level driver                                                        */
/* ------------------------------------------------------------------ */

static volatile sig_atomic_t stop_flag;
static void on_sigint(int sig) { (void)sig; stop_flag = 1; }

typedef struct {
    uint64_t lev_lo, lev_hi;       /* L_lo range of this level [lo, hi) */
    uint64_t Klo, Khi;             /* K range */
    uint64_t chunk;
    _Atomic uint64_t nextK;
    /* frontier: chunks complete */
    uint8_t *done;
    uint64_t nchunks, frontier;    /* frontier = first incomplete chunk */
    pthread_mutex_t lock;
    pthread_cond_t cond;
    int nexited;
    Worker *workers;
    int all_words;                 /* range mode: every word, not just this level */
} LevelJob;

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void do_K(Worker *W, LevelJob *J, uint64_t K)
{
    uint64_t xc = atomic_load(&Xcap);
    uint64_t base = (uint64_t)(KLEN - 1) * K * PRIM;
    if (base + SMALL_A + 1 > xc) return;
    kctx_build(&W->c, K);
    uint64_t Wd = 64 * W->c.MOD;
    uint64_t wlo, whi;          /* words with L_lo in [lev_lo, lev_hi) */
    if (J->all_words) { wlo = 0; whi = UINT64_MAX; }
    else {
        wlo = J->lev_lo > base ? (J->lev_lo - base + Wd - 1) / Wd : 0;
        if (J->lev_hi <= base) return;
        whi = (J->lev_hi - base + Wd - 1) / Wd;      /* exclusive */
    }
    for (uint64_t w = wlo; w < whi; w++) {
        xc = atomic_load(&Xcap);
        if (base + w * Wd > xc) break;
        uint64_t A_hi = xc - base;
        item_engine(W, w, A_hi);
        if (stop_flag && J->all_words) break;
    }
}

static void *level_thread(void *arg)
{
    void **pa = arg;
    LevelJob *J = pa[0];
    Worker *W = pa[1];
    for (;;) {
        if (stop_flag) break;
        uint64_t K0 = atomic_fetch_add(&J->nextK, J->chunk);
        if (K0 > J->Khi) break;
        uint64_t K1 = K0 + J->chunk - 1;
        if (K1 > J->Khi) K1 = J->Khi;
        for (uint64_t K = K0; K <= K1; K++) do_K(W, J, K);
        if (stop_flag) break;                 /* chunk may be incomplete */
        uint64_t ci = (K0 - J->Klo) / J->chunk;
        pthread_mutex_lock(&J->lock);
        J->done[ci] = 1;
        while (J->frontier < J->nchunks && J->done[J->frontier]) J->frontier++;
        pthread_mutex_unlock(&J->lock);
    }
    pthread_mutex_lock(&J->lock);
    J->nexited++;
    pthread_cond_signal(&J->cond);
    pthread_mutex_unlock(&J->lock);
    return NULL;
}

/* checkpoint: k X DELTA MODMAX N2 level frontierK */
static const char *state_file;
static int ckpt_secs = 60;

static void write_state(uint64_t X, uint64_t delta, uint64_t level, uint64_t frontierK)
{
    if (!state_file) return;
    char tmp[1024];
    snprintf(tmp, sizeof tmp, "%s.tmp", state_file);
    FILE *f = fopen(tmp, "w");
    if (!f) { perror(tmp); return; }
    fprintf(f, "%d %" PRIu64 " %" PRIu64 " %" PRIu64 " %d %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
            KLEN, X, delta, MODMAX, N2, level, frontierK, atomic_load(&Xcap));
    pthread_mutex_lock(&res_lock);
    for (size_t i = 0; i < nfound; i++)
        fprintf(f, "ap %" PRIu64 " %" PRIu64 "\n", found[i].a, found[i].d);
    pthread_mutex_unlock(&res_lock);
    fclose(f);
    rename(tmp, state_file);
}

static int read_state(uint64_t X, uint64_t delta, uint64_t *level, uint64_t *frontierK)
{
    if (!state_file) return 0;
    FILE *f = fopen(state_file, "r");
    if (!f) return 0;
    int k, n2;
    uint64_t x, dl, mm, lv, fk, xc;
    if (fscanf(f, "%d %" SCNu64 " %" SCNu64 " %" SCNu64 " %d %" SCNu64 " %" SCNu64 " %" SCNu64,
               &k, &x, &dl, &mm, &n2, &lv, &fk, &xc) != 8) { fclose(f); return 0; }
    if (k != KLEN || x != X || dl != delta || mm != MODMAX || n2 != N2) {
        fprintf(stderr, "state file %s does not match these parameters\n", state_file);
        exit(1);
    }
    uint64_t a, d;
    int saved = quiet; quiet = 1;
    while (fscanf(f, " ap %" SCNu64 " %" SCNu64, &a, &d) == 2) report_ap(a, d);
    quiet = saved;
    fclose(f);
    atomic_store(&Xcap, xc);
    *level = lv; *frontierK = fk;
    return 1;
}

static void run_level(LevelJob *J, Worker *workers, uint64_t X, uint64_t delta, uint64_t level)
{
    J->nchunks = (J->Khi - J->Klo) / J->chunk + 1;
    J->done = calloc(J->nchunks, 1);
    J->frontier = 0;
    atomic_init(&J->nextK, J->Klo);
    pthread_mutex_init(&J->lock, NULL);
    pthread_cond_init(&J->cond, NULL);
    J->nexited = 0;
    pthread_t th[256];
    void *args[256][2];
    uint64_t w0 = 0, c20 = 0;
    for (int i = 0; i < nthreads; i++) { w0 += workers[i].words; c20 += workers[i].cand2; }
    double t0 = now_sec(), tlast = t0, tck = t0;
    for (int i = 0; i < nthreads; i++) {
        args[i][0] = J; args[i][1] = &workers[i];
        pthread_create(&th[i], NULL, level_thread, args[i]);
    }
    /* monitor */
    for (;;) {
        pthread_mutex_lock(&J->lock);
        if (J->nexited < nthreads) {
            struct timespec ts;
            clock_gettime(CLOCK_REALTIME, &ts);
            ts.tv_nsec += 250000000;
            if (ts.tv_nsec >= 1000000000) { ts.tv_sec++; ts.tv_nsec -= 1000000000; }
            pthread_cond_timedwait(&J->cond, &J->lock, &ts);
        }
        int finished = J->nexited == nthreads;
        pthread_mutex_unlock(&J->lock);
        uint64_t nk = atomic_load(&J->nextK);
        double t = now_sec();
        if (!quiet && t - tlast >= 2.0) {
            uint64_t wsum = 0, c3 = 0;
            for (int i = 0; i < nthreads; i++) { wsum += workers[i].words; c3 += workers[i].cand3; }
            fprintf(stderr, "\rlevel %" PRIu64 " L<%.4g  K %" PRIu64 "/%" PRIu64 "  %.3g words/s  cand3 %" PRIu64 "  APs %zu  cap %.6g   ",
                    level, (double)J->lev_hi, nk > J->Khi ? J->Khi : nk, J->Khi,
                    (wsum - w0) / (t - t0), c3, nfound, (double)atomic_load(&Xcap));
            tlast = t;
        }
        if (state_file && t - tck >= ckpt_secs) {
            pthread_mutex_lock(&J->lock);
            uint64_t fk = J->Klo + J->frontier * J->chunk;
            pthread_mutex_unlock(&J->lock);
            write_state(X, delta, level, fk);
            tck = t;
        }
        if (finished || stop_flag) break;
    }
    for (int i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
    (void)c20;
    if (stop_flag) {
        uint64_t fk = J->Klo + J->frontier * J->chunk;
        if (fk > J->Khi + 1) fk = J->Khi + 1;
        write_state(X, delta, level, fk);
        if (!quiet) fprintf(stderr, "\ninterrupted: level %" PRIu64 " complete below K = %" PRIu64 "\n", level, fk);
    }
    free(J->done);
    pthread_mutex_destroy(&J->lock);
    pthread_cond_destroy(&J->cond);
}

static void setup_k(int k)
{
    KLEN = k;
    PRIM = primorial(k);
}

static uint64_t default_modmax(uint64_t X)
{
    uint64_t m = X / 1024;
    return m < 30 ? 30 : m;
}

/* mode 0: find a(k) <= X (min, certified); mode 1: all APs with L <= X */
static int run_search(int mode, uint64_t X, uint64_t delta, uint64_t *answer)
{
    search_min_mode = (mode == 0);
    atomic_store(&Xcap, X);
    nfound = 0;
    Worker *workers = calloc(nthreads, sizeof(Worker));
    for (int i = 0; i < nthreads; i++) { kctx_alloc(&workers[i].c); if (worker_init_hook) worker_init_hook(&workers[i]); }
    uint64_t Wnom = 64 * mod_for_K(1);
    if (!delta) delta = Wnom;
    uint64_t level = 0, frontierK = 1;
    int resumed = read_state(X, delta, &level, &frontierK);
    if (!quiet)
        fprintf(stderr, "k=%d P=%" PRIu64 " X=%" PRIu64 " MODMAX=%" PRIu64 " MOD(K=1)=%" PRIu64
                " DELTA=%" PRIu64 " N2=%d threads=%d%s\n",
                KLEN, PRIM, X, MODMAX, Wnom / 64, delta, N2, nthreads, resumed ? " (resumed)" : "");
    if (!resumed || (level == 0 && frontierK == 1)) {
        small_pass(X, 1, UINT64_MAX);
        if (!quiet) fprintf(stderr, "small first terms (a <= %u) done: %zu APs\n", SMALL_A, nfound);
    }
    uint64_t certified = 0;
    int done = 0;
    for (; !done && !stop_flag; level++) {
        uint64_t xc = atomic_load(&Xcap);
        LevelJob J;
        memset(&J, 0, sizeof J);
        J.lev_lo = level * delta;
        J.lev_hi = (level + 1) * delta;
        if (J.lev_lo > xc) { certified = xc + 1; break; }
        uint64_t top = J.lev_hi - 1 < xc ? J.lev_hi - 1 : xc;
        J.Klo = frontierK;
        J.Khi = top > SMALL_A ? (top - SMALL_A - 1) / ((uint64_t)(KLEN - 1) * PRIM) : 0;
        frontierK = 1;
        if (J.Khi >= J.Klo) {
            uint64_t nK = J.Khi - J.Klo + 1;
            uint64_t ch = nK / ((uint64_t)nthreads * 64);
            J.chunk = ch < 1 ? 1 : ch > 4096 ? 4096 : ch;
            run_level(&J, workers, X, delta, level);
        }
        if (stop_flag) break;
        certified = J.lev_hi;
        xc = atomic_load(&Xcap);
        if (mode == 0 && nfound && xc < certified) done = 1;
        if (certified > xc) done = 1;
        write_state(X, delta, level + 1, 1);
    }
    if (!quiet) fprintf(stderr, "\n");
    qsort(found, nfound, sizeof(AP), cmp_ap);
    uint64_t wsum = 0, c2 = 0, c3 = 0;
    for (int i = 0; i < nthreads; i++) { wsum += workers[i].words; c2 += workers[i].cand2; c3 += workers[i].cand3; }
    if (!quiet)
        fprintf(stderr, "words %" PRIu64 "  stage-2 survivors %" PRIu64 "  stage-3 survivors %" PRIu64 "  sprp %" PRIu64 "\n",
                wsum, c2, c3, (uint64_t)atomic_load(&n_sprp_calls));
    int ok = !stop_flag;
    if (ok && mode == 0 && print_aps) {
        if (nfound) {
            *answer = found[0].L;
            printf("a(%d) = %" PRIu64 "  (a = %" PRIu64 ", d = %" PRIu64 "%s)\n", KLEN, found[0].L, found[0].a,
                   found[0].d % PRIM == 0 ? found[0].d / PRIM : found[0].d, found[0].d % PRIM == 0 ? "*k#" : "");
        } else {
            *answer = 0;
            printf("no AP-%d with last term <= %" PRIu64 "\n", KLEN, X);
        }
    }
    if (ok && mode == 0 && nfound) *answer = found[0].L;
    if (ok && mode == 0 && !nfound) *answer = 0;
    if (ok && mode == 1 && print_aps) printf("%zu AP-%d with last term <= %" PRIu64 "\n", nfound, KLEN, X);
    for (int i = 0; i < nthreads; i++) {
        free(workers[i].c.bad);
        if (!workers[i].eng) free(workers[i].c.tab);
        free(workers[i].c.q3); free(workers[i].c.mdinv3); free(workers[i].c.qd3);
    }
    free(workers);
    return ok;
}

/* every AP-k with L <= X and K in [K0, K1] (range mode) */
static void run_range(uint64_t X, uint64_t K0, uint64_t K1)
{
    search_min_mode = 0;
    atomic_store(&Xcap, X);
    Worker *workers = calloc(nthreads, sizeof(Worker));
    for (int i = 0; i < nthreads; i++) { kctx_alloc(&workers[i].c); if (worker_init_hook) worker_init_hook(&workers[i]); }
    small_pass(X, K0, K1);
    LevelJob J;
    memset(&J, 0, sizeof J);
    J.all_words = 1;
    J.Klo = K0;
    uint64_t kmax = X > SMALL_A ? (X - SMALL_A - 1) / ((uint64_t)(KLEN - 1) * PRIM) : 0;
    J.Khi = K1 < kmax ? K1 : kmax;
    J.chunk = 1;
    if (J.Khi >= J.Klo) run_level(&J, workers, X, 0, 0);
    if (!quiet) fprintf(stderr, "\n");
    qsort(found, nfound, sizeof(AP), cmp_ap);
    uint64_t wsum = 0, c2 = 0, c3 = 0;
    for (int i = 0; i < nthreads; i++) { wsum += workers[i].words; c2 += workers[i].cand2; c3 += workers[i].cand3; }
    printf("range K=[%" PRIu64 ",%" PRIu64 "] X=%" PRIu64 ": %zu AP-%d; words %" PRIu64 " stage2 %" PRIu64 " stage3 %" PRIu64 "\n",
           K0, K1, X, nfound, KLEN, wsum, c2, c3);
}

/* ------------------------------------------------------------------ */
/* Known data                                                          */
/* ------------------------------------------------------------------ */

static const uint64_t A005115[27] = {0,
    2, 3, 7, 23, 29, 157, 907, 1669, 1879, 2089, 249037, 262897, 725663,
    36850999, 173471351, 198793279, 4827507229ULL, 17010526363ULL, 83547839407ULL,
    572945039351ULL, 6269243827111ULL, 35742689530423ULL, 449924511422857ULL,
    1217585417914253ULL, 5773236905395679ULL, 12783396861134173ULL};

/* first term and difference (A113827, A093364 / Luhn's record page) */
static const struct { uint64_t a, d; } KNOWN_AP[27] = {
    {0, 0}, {2, 0}, {2, 1}, {3, 2}, {5, 6}, {5, 6}, {7, 30}, {7, 150}, {199, 210}, {199, 210},
    {199, 210}, {110437, 13860}, {110437, 13860}, {4943, 60060},
    {1045ULL * 30030 + 4189, 14ULL * 30030}, {3844ULL * 30030 + 18071, 138ULL * 30030},
    {1774ULL * 30030 + 24709, 323ULL * 30030}, {6720ULL * 510510 + 124669, 171ULL * 510510},
    {9418ULL * 510510 + 333163, 1406ULL * 510510}, {855ULL * 9699690 + 4409437, 431ULL * 9699690},
    {22151ULL * 9699690 + 3750431, 1943ULL * 9699690}, {592714ULL * 9699690 + 4390651, 2681ULL * 9699690},
    {1985821ULL * 9699690 + 1159033, 80910ULL * 9699690}, {1807252ULL * 223092870 + 181107397, 9523ULL * 223092870},
    {2310638ULL * 223092870 + 83578883, 136831ULL * 223092870}, {2526681ULL * 223092870 + 128097689, 972979ULL * 223092870},
    {15626261ULL * 223092870 + 59138353, 1666981ULL * 223092870}};

static int selftest(int kmax)
{
    int fails = 0;
    /* primality sanity */
    static const uint64_t pr[] = {2, 3, 5, 1000000007ULL, 18446744073709551557ULL, 4611686018427387847ULL};
    static const uint64_t cp[] = {1, 4, 561, 3215031751ULL, 341550071728321ULL, 3825123056546413051ULL};
    for (size_t i = 0; i < sizeof pr / sizeof pr[0]; i++)
        if (!is_prime_u64(pr[i])) { printf("FAIL prime %" PRIu64 "\n", pr[i]); fails++; }
    for (size_t i = 0; i < sizeof cp / sizeof cp[0]; i++)
        if (is_prime_u64(cp[i])) { printf("FAIL composite %" PRIu64 "\n", cp[i]); fails++; }
    /* the known progressions are all prime and end at a(k) */
    for (int k = 2; k <= 26; k++) {
        uint64_t a = KNOWN_AP[k].a, d = KNOWN_AP[k].d;
        int ok = 1;
        for (int i = 0; i < k; i++) if (!is_prime_u64(a + (uint64_t)i * d)) ok = 0;
        if (!ok || a + (uint64_t)(k - 1) * d != A005115[k]) {
            printf("FAIL known AP-%d data\n", k); fails++;
        }
    }
    /* the Gahan AP-27 */
    {
        uint64_t a = 224584605939537911ULL, d = 81292139ULL * 223092870ULL;
        int ok = 1;
        for (int i = 0; i < 27; i++) if (!is_prime_u64(a + (uint64_t)i * d)) ok = 0;
        if (!ok || a + 26 * d != 696112717486210091ULL) { printf("FAIL Gahan AP-27\n"); fails++; }
    }
    int saved_quiet = quiet;
    print_aps = 0;
    for (int k = 3; k <= kmax; k++) {
        setup_k(k);
        uint64_t X = A005115[k];
        MODMAX = default_modmax(X);
        quiet = 1;
        double t0 = now_sec();
        uint64_t ans = 0;
        nfound = 0;
        n_sprp_calls = 0;
        run_search(0, X, 0, &ans);
        double t = now_sec() - t0;
        int ok = (ans == A005115[k]);
        int apok = ok && nfound >= 1 && found[0].a == KNOWN_AP[k].a && found[0].d == KNOWN_AP[k].d;
        /* a second run with a looser cap must stop at the same place */
        uint64_t ans2 = 0;
        if (ok && k <= 18) {
            nfound = 0;
            run_search(0, 3 * X, 0, &ans2);
            if (ans2 != ans) ok = 0;
        }
        quiet = saved_quiet;
        printf("k=%2d  a(k) = %" PRIu64 "  %s  AP %s  (%.2f s)\n", k, ans, ok ? "ok" : "FAIL", apok ? "ok" : "FAIL", t);
        fflush(stdout);
        if (!ok || !apok) fails++;
    }
    printf("%s (%d failures)\n", fails ? "SELFTEST FAILED" : "selftest passed", fails);
    return fails ? 1 : 0;
}

/* ------------------------------------------------------------------ */
/* Bench                                                               */
/* ------------------------------------------------------------------ */

static void bench(int k, uint64_t X, uint64_t K, uint64_t w, double secs)
{
    setup_k(k);
    if (!MODMAX) MODMAX = default_modmax(X);
    atomic_store(&Xcap, X);
    Worker W;
    memset(&W, 0, sizeof W);
    kctx_alloc(&W.c);
    kctx_build(&W.c, K);
    printf("k=%d K=%" PRIu64 " MOD=%" PRIu64 " n1=%d residues=%" PRIu64 " N2=%d (q2 %u..%u) n3=%d\n", k, K, W.c.MOD, W.c.n1,
           W.c.nR, W.c.n2, W.c.q2[0], W.c.q2[W.c.n2 - 1], W.c.n3);
    printf("Q1:");
    for (int j = 0; j < W.c.n1; j++) printf(" %u(%u)", W.c.q1[j], W.c.g1[j]);
    printf("\n");
    uint64_t base = (uint64_t)(k - 1) * K * PRIM;
    uint64_t A_hi = X > base ? X - base : 0;
    double t0 = now_sec();
    int items = 0;
    do { process_item(&W, w, A_hi); items++; } while (now_sec() - t0 < secs);
    double t = now_sec() - t0;
    printf("%d items, %" PRIu64 " words in %.2f s: %.3g words/s per thread; stage-2 survivors/word %.4g, stage-3 %" PRIu64 "\n",
           items, W.words, t, W.words / t, (double)W.cand2 / W.words, W.cand3);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

#ifndef A005115_NO_MAIN
static uint64_t parse_u64(const char *s)
{
    char *end;
    double dv = strtod(s, &end);
    if (strchr(s, 'e') || strchr(s, 'E') || strchr(s, '.')) return (uint64_t)llround(dv);
    return strtoull(s, NULL, 10);
}

static void usage(void)
{
    fprintf(stderr,
            "usage: a005115 selftest [KMAX]\n"
            "       a005115 search k X [-t T] [-m MODMAX] [-l DELTA] [-n N2] [-S FILE] [-i SECS] [-o FILE] [-q]\n"
            "       a005115 all k X [options]\n"
            "       a005115 range k X K0 K1 [options]\n"
            "       a005115 bench [k X K W SECS]\n"
            "       a005115 verify a d k\n");
    exit(2);
}

int main(int argc, char **argv)
{
    init_small_primes();
    if (argc < 2) usage();
    const char *cmd = argv[1];
    uint64_t pos[8];
    int npos = 0;
    uint64_t delta = 0;
    const char *outfile = NULL;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) nthreads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-m") && i + 1 < argc) MODMAX = parse_u64(argv[++i]);
        else if (!strcmp(argv[i], "-l") && i + 1 < argc) delta = parse_u64(argv[++i]);
        else if (!strcmp(argv[i], "-n") && i + 1 < argc) N2 = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-S") && i + 1 < argc) state_file = argv[++i];
        else if (!strcmp(argv[i], "-i") && i + 1 < argc) ckpt_secs = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i + 1 < argc) outfile = argv[++i];
        else if (!strcmp(argv[i], "-q")) quiet = 1;
        else if (npos < 8) pos[npos++] = parse_u64(argv[i]);
        else usage();
    }
    if (N2 < 1 || N2 > N2MAX) { fprintf(stderr, "N2 must be 1..%d\n", N2MAX); return 2; }
    if (nthreads <= 0) {
        long n = sysconf(_SC_NPROCESSORS_ONLN);
        nthreads = n > 0 ? (int)n : 4;
    }
    if (nthreads > 256) nthreads = 256;
    signal(SIGINT, on_sigint);
    signal(SIGTERM, on_sigint);
    if (outfile) results_fp = fopen(outfile, "a");

    if (!strcmp(cmd, "selftest")) {
        int kmax = npos >= 1 ? (int)pos[0] : 20;
        return selftest(kmax);
    }
    if (!strcmp(cmd, "verify")) {
        if (npos < 3) usage();
        setup_k((int)pos[2]);
        int k = (int)pos[2];
        int bad = -1;
        for (int i = 0; i < k; i++) if (!is_prime_u64(pos[0] + (uint64_t)i * pos[1])) { bad = i; break; }
        if (bad < 0) printf("AP-%d of primes: %" PRIu64 " + %" PRIu64 "*j, last term %" PRIu64 "\n", k, pos[0], pos[1], pos[0] + (uint64_t)(k - 1) * pos[1]);
        else printf("term %d = %" PRIu64 " is composite\n", bad, pos[0] + (uint64_t)bad * pos[1]);
        return bad < 0 ? 0 : 1;
    }
    if (!strcmp(cmd, "bench")) {
        int k = npos >= 1 ? (int)pos[0] : 27;
        uint64_t X = npos >= 2 ? pos[1] : 696112717486210091ULL;
        uint64_t K = npos >= 3 ? pos[2] : 1000003;
        uint64_t w = npos >= 4 ? pos[3] : 3;
        double secs = npos >= 5 ? (double)pos[4] : 10;
        bench(k, X, K, w, secs);
        return 0;
    }
    if (!strcmp(cmd, "search") || !strcmp(cmd, "all")) {
        if (npos < 2) usage();
        setup_k((int)pos[0]);
        if (KLEN < 3 || KLEN > 40) { fprintf(stderr, "k must be 3..40\n"); return 2; }
        uint64_t X = pos[1];
        if (!MODMAX) MODMAX = default_modmax(X);
        uint64_t ans;
        int ok = run_search(!strcmp(cmd, "all"), X, delta, &ans);
        return ok ? 0 : 3;
    }
    if (!strcmp(cmd, "range")) {
        if (npos < 4) usage();
        setup_k((int)pos[0]);
        uint64_t X = pos[1];
        if (!MODMAX) MODMAX = default_modmax(X);
        run_range(X, pos[2], pos[3]);
        return 0;
    }
    usage();
    return 2;
}
#endif /* A005115_NO_MAIN */
