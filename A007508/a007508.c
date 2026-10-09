/*
 * a007508.c - Count twin prime pairs below 10^n  (OEIS A007508)
 *
 *   a(n) = #{ p : p and p+2 both prime, p+2 < 10^n }
 *
 * Everything is done in "k-space": for p > 5 a twin pair (p, p+2) has
 * p = 30k + r with r in {11, 17, 29}, so we index candidates by k = p/30 and
 * a class index j in {0,1,2}.  10^n = 30*K_n + 10, so "p < 10^n" is exactly
 * "k < K_n" with K_n = (10^n - 10)/30, which fits in 64 bits for n = 20.
 * The pairs (3,5) and (5,7) are added by hand.
 *
 * Kernels (-k):
 *   twin  own segmented sieve over the 3 twin residue classes mod 30, with a
 *         presieve for 7..19, unrolled loops for small primes, a wheel walk
 *         for medium primes and a bucket sieve for large primes.  Works to
 *         10^20.  (default)
 *   ps    primesieve_count_twins per chunk (64-bit only, n <= 19).
 *   iter  primesieve_iterator, count consecutive primes 2 apart (n <= 19).
 *
 * Parallelism: pthreads pull chunks of a decade from an atomic counter.
 * Resume/checkpoint: -l log, -f/-a, -r A:B  (see usage()).
 */

#define _GNU_SOURCE
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include <primesieve.h>

typedef unsigned __int128 u128;

/* ======================================================================= */
/*  Wheel tables for the twin sieve                                        */
/* ======================================================================= */

static const int R3[3] = {11, 17, 29};                 /* twin classes mod 30 */
static const int W8[8] = {1, 7, 11, 13, 17, 19, 23, 29};
static int8_t WIDX[30];

/* For a sieving prime q = 30*qd + w and class r_j, the candidates n = 30k+r_j
 * with n ≡ t (mod q), t in {0,-2}, form the progression k ≡ k0 (mod q) with
 * k0 = WM*qd + WC (the hit whose cofactor (n-t)/q equals WM < 30).
 * Sorted by WM the six progressions give the wheel order u = 0..5.        */
static int32_t WM[8][6], WC[8][6];
static uint8_t WJ[8][6];
static int8_t  WT[8][6];
static int32_t DKM[8][6], DKC[8][6];   /* k-delta hit u -> u+1: DKM*qd + DKC  */
static int32_t DBM[8][6], DBC[8][6];   /* same delta in bit units (3k + j)    */
static int     MAXDM;

typedef struct { int m, c, j, t; } wentry_t;

static int modinv30(int w)
{
    for (int x = 1; x < 30; x++) if ((w * x) % 30 == 1) return x;
    return -1;
}

static void wheel_init(void)
{
    memset(WIDX, -1, sizeof WIDX);
    for (int i = 0; i < 8; i++) WIDX[W8[i]] = (int8_t)i;
    MAXDM = 0;
    for (int wi = 0; wi < 8; wi++) {
        int w = W8[wi], winv = modinv30(w);
        wentry_t e[6];
        int cnt = 0;
        for (int j = 0; j < 3; j++)
            for (int ti = 0; ti < 2; ti++) {
                int t = ti ? -2 : 0;
                int d = R3[j] - t;                 /* 11,13,17,19,29,31 */
                int m = (d * winv) % 30;
                int num = m * w - d;               /* divisible by 30 */
                if (num % 30 != 0) { fprintf(stderr, "wheel bug\n"); exit(1); }
                e[cnt].m = m; e[cnt].c = num / 30; e[cnt].j = j; e[cnt].t = t;
                cnt++;
            }
        for (int a = 1; a < 6; a++)
            for (int b = a; b > 0 && e[b - 1].m > e[b].m; b--) {
                wentry_t tmp = e[b]; e[b] = e[b - 1]; e[b - 1] = tmp;
            }
        for (int u = 0; u < 6; u++) {
            WM[wi][u] = e[u].m; WC[wi][u] = e[u].c;
            WJ[wi][u] = (uint8_t)e[u].j; WT[wi][u] = (int8_t)e[u].t;
        }
        for (int u = 0; u < 6; u++) {
            int v = (u + 1) % 6;
            int dm = e[v].m - e[u].m + (v == 0 ? 30 : 0);
            int dc = e[v].c - e[u].c + (v == 0 ? w : 0);
            DKM[wi][u] = dm; DKC[wi][u] = dc;
            DBM[wi][u] = 3 * dm; DBC[wi][u] = 3 * dc + (e[v].j - e[u].j);
            if (dm > MAXDM) MAXDM = dm;
        }
    }
}

/* ======================================================================= */
/*  Presieve pattern for 7, 11, 13, 17, 19                                 */
/* ======================================================================= */

#define PS_NP     5
static const int PSP[PS_NP] = {7, 11, 13, 17, 19};
#define PS_PERIOD 323323ULL                 /* 7*11*13*17*19, in k-units    */
#define PS_BYTES  (3 * PS_PERIOD)           /* byte period of the bit pattern */
static uint8_t *PAT;                        /* bit 3k+j for absolute k < 8P  */

static void presieve_init(void)
{
    PAT = malloc(PS_BYTES);
    if (!PAT) { perror("malloc"); exit(1); }
    memset(PAT, 0xFF, PS_BYTES);
    for (int pi = 0; pi < PS_NP; pi++) {
        int q = PSP[pi];
        for (int j = 0; j < 3; j++)
            for (int k0 = 0; k0 < q; k0++) {
                int n = (30 * k0 + R3[j]) % q;
                if (n != 0 && n != q - 2) continue;
                for (uint64_t k = (uint64_t)k0; k < 8 * PS_PERIOD; k += (uint64_t)q) {
                    uint64_t b = 3 * k + (uint64_t)j;
                    PAT[b >> 3] &= (uint8_t)~(1u << (b & 7));
                }
            }
    }
}

/* ======================================================================= */
/*  Twin sieve kernel                                                      */
/* ======================================================================= */

/* Stacked presieve patterns for the primes 23..97, on top of PAT (7..19). */
#define NPAT 6
static const int PATP[NPAT][3] = {{23,29,31},{37,41,43},{47,53,59},{61,67,71},{73,79,83},{89,97,0}};
static uint64_t PATPER[NPAT];          /* period in k-units                    */
static uint64_t PATLEN[NPAT];          /* byte period = 3*period               */
static uint8_t *PATB[NPAT];

static void presieve2_init(void)
{
    for (int g = 0; g < NPAT; g++) {
        uint64_t per = 1;
        for (int i = 0; i < 3; i++) if (PATP[g][i]) per *= (uint64_t)PATP[g][i];
        PATPER[g] = per; PATLEN[g] = 3 * per;
        uint8_t *b = malloc(PATLEN[g]);
        if (!b) { perror("malloc"); exit(1); }
        memset(b, 0xFF, PATLEN[g]);
        for (int i = 0; i < 3; i++) {
            int q = PATP[g][i];
            if (!q) continue;
            for (int j = 0; j < 3; j++)
                for (int k0 = 0; k0 < q; k0++) {
                    int n = (30 * k0 + R3[j]) % q;
                    if (n != 0 && n != q - 2) continue;
                    for (uint64_t k = (uint64_t)k0; k < 8 * per; k += (uint64_t)q) {
                        uint64_t bb = 3 * k + (uint64_t)j;
                        b[bb >> 3] &= (uint8_t)~(1u << (bb & 7));
                    }
                }
        }
        PATB[g] = b;
    }
}

/* generic re-alignment of a byte-periodic pattern to chunk start K0 */
static void pat_align(uint8_t *dst, const uint8_t *src, uint64_t len, uint64_t K0)
{
    uint64_t B = 3 * K0;
    unsigned ph = (unsigned)(B & 7);
    uint64_t off = (B >> 3) % len;
    for (uint64_t i = 0; i < len; i++) {
        uint64_t i0 = off + i; if (i0 >= len) i0 -= len;
        uint64_t i1 = i0 + 1;  if (i1 >= len) i1 -= len;
        dst[i] = ph ? (uint8_t)((src[i0] >> ph) | (src[i1] << (8 - ph))) : src[i0];
    }
}

/* dst[0..n) = (op==0 ? src : dst & src), src taken cyclically from pattern */
static void pat_apply(uint8_t *dst, uint64_t n, const uint8_t *pat, uint64_t len, uint64_t start, int op)
{
    uint64_t pi = start % len, done = 0;
    while (done < n) {
        uint64_t m = len - pi;
        if (m > n - done) m = n - done;
        if (op == 0) memcpy(dst + done, pat + pi, m);
        else { uint8_t *d = dst + done; const uint8_t *s = pat + pi; for (uint64_t i = 0; i < m; i++) d[i] &= s[i]; }
        done += m; pi += m; if (pi >= len) pi = 0;
    }
}

static int is_prime_small(uint64_t n)
{
    if (n < 2) return 0;
    for (uint64_t d = 2; d * d <= n; d++) if (n % d == 0) return 0;
    return 1;
}

typedef struct { uint64_t pos[6]; uint32_t q; uint32_t pad; } sprime_t;   /* tiny/small */
typedef struct { uint64_t pos; uint32_t qd; uint8_t wi, u; uint16_t pad; } mprime_t; /* medium */

#define BLK 1023
typedef struct block { struct block *next; uint64_t rec[BLK]; } block_t;   /* 8192 bytes */
typedef struct { uint64_t *cur, *end; block_t *head, *tail; } bucket_t;

#define TINY_MAX  2048
#define SUB_BITS  (32 * 1024 * 8)      /* L1-sized sub-block for tiny primes */

typedef struct {
    uint64_t K0, K1, len;
    uint64_t seg_k, seg_bits, seg_bytes;
    int      seg_log;
    uint8_t *seg, *cp0, *cpx[NPAT];
    sprime_t *tp; size_t ntp, ctp;     /* tiny  : 101 <= q <= TINY_MAX        */
    sprime_t *sp; size_t nsp, csp;     /* small : TINY_MAX < q <= seg_k/64    */
    mprime_t *mp; size_t nmp, cmp;     /* medium: seg_k/64 < q <= seg_k       */
    bucket_t *bk; uint64_t nb, nbm;    /* large : q > seg_k                   */
    block_t *freelist;
    void   **arenas; size_t narenas, carenas;
} tw_t;

static block_t *blk_get(tw_t *T)
{
    if (!T->freelist) {
        size_t n = 512;
        block_t *a = malloc(n * sizeof(block_t));
        if (!a) { fprintf(stderr, "out of memory (bucket arena)\n"); exit(1); }
        if (T->narenas == T->carenas) {
            T->carenas = T->carenas ? 2 * T->carenas : 64;
            T->arenas = realloc(T->arenas, T->carenas * sizeof(void *));
            if (!T->arenas) { perror("realloc"); exit(1); }
        }
        T->arenas[T->narenas++] = a;
        for (size_t i = 0; i < n; i++) { a[i].next = T->freelist; T->freelist = &a[i]; }
    }
    block_t *b = T->freelist;
    T->freelist = b->next;
    b->next = NULL;
    return b;
}

static void bk_grow(tw_t *T, bucket_t *B)
{
    block_t *nb = blk_get(T);
    if (B->tail) B->tail->next = nb; else B->head = nb;
    B->tail = nb; B->cur = nb->rec; B->end = nb->rec + BLK;
}

static inline void bk_push(tw_t *T, uint64_t idx, uint64_t rec)
{
    bucket_t *B = &T->bk[idx];
    if (B->cur == B->end) bk_grow(T, B);
    *B->cur++ = rec;
}

static inline uint64_t pack_rec(uint64_t koff, unsigned j, unsigned u, unsigned wi, uint64_t qd)
{
    return koff | ((uint64_t)j << 24) | ((uint64_t)u << 26) | ((uint64_t)wi << 29) | (qd << 32);
}

static void push_sprime(sprime_t **arr, size_t *n, size_t *cap, uint64_t q, const uint64_t *kk, uint64_t K0, int wi)
{
    if (*n == *cap) {
        *cap = *cap ? 2 * *cap : 1024;
        *arr = realloc(*arr, *cap * sizeof(sprime_t));
        if (!*arr) { perror("realloc"); exit(1); }
    }
    sprime_t *S = &(*arr)[(*n)++];
    S->q = (uint32_t)q;
    for (int u = 0; u < 6; u++) S->pos[u] = 3 * (kk[u] - K0) + WJ[wi][u];
}

/* Place sieving prime q (>= 101, coprime to 30) into its tier for this chunk. */
static void tw_add_prime(tw_t *T, uint64_t q)
{
    uint64_t qd = q / 30;
    int wi = WIDX[q % 30];
    uint64_t K0 = T->K0;
    uint64_t kk[6];

    if (K0 >= q) {
        /* fast path: the cofactor-1 hits (n = q, n = q-2) lie below the chunk */
        uint64_t base = K0 % q;
        for (int u = 0; u < 6; u++) {
            int64_t v = (int64_t)WM[wi][u] * (int64_t)qd + WC[wi][u];
            if (v >= (int64_t)q) v %= (int64_t)q;
            uint64_t rho = (uint64_t)v;
            uint64_t off = rho >= base ? rho - base : rho + q - base;
            kk[u] = K0 + off;
        }
    } else {
        /* slow path (chunk starts below q): require n - t >= 7q so the
         * cofactor-1 hit is excluded                                    */
        for (int u = 0; u < 6; u++) {
            int64_t nb = 7 * (int64_t)q + WT[wi][u] - R3[WJ[wi][u]];
            uint64_t kb = (uint64_t)((nb + 29) / 30);
            if (kb < K0) kb = K0;
            int64_t v = (int64_t)WM[wi][u] * (int64_t)qd + WC[wi][u];
            v %= (int64_t)q; if (v < 0) v += (int64_t)q;
            uint64_t rho = (uint64_t)v;
            uint64_t base = kb % q;
            uint64_t off = rho >= base ? rho - base : rho + q - base;
            kk[u] = kb + off;
        }
    }

    if (q <= TINY_MAX)       { push_sprime(&T->tp, &T->ntp, &T->ctp, q, kk, K0, wi); return; }
    if (q <= T->seg_k / 64)  { push_sprime(&T->sp, &T->nsp, &T->csp, q, kk, K0, wi); return; }

    int um = 0;
    for (int u = 1; u < 6; u++) if (kk[u] < kk[um]) um = u;
    uint64_t krel = kk[um] - K0;
    if (krel >= T->len) return;

    if (q <= T->seg_k) {
        if (T->nmp == T->cmp) {
            T->cmp = T->cmp ? 2 * T->cmp : 4096;
            T->mp = realloc(T->mp, T->cmp * sizeof(mprime_t));
            if (!T->mp) { perror("realloc"); exit(1); }
        }
        mprime_t *M = &T->mp[T->nmp++];
        M->pos = 3 * krel + WJ[wi][um];
        M->qd = (uint32_t)qd; M->wi = (uint8_t)wi; M->u = (uint8_t)um;
        return;
    }
    uint64_t sg = krel >> T->seg_log;
    bk_push(T, sg & T->nbm, pack_rec(krel & (T->seg_k - 1), WJ[wi][um], (unsigned)um, (unsigned)wi, qd));
}

/* Clear bits sp, sp+step, ... below bend (bit indices into seg); returns the
 * first position >= bend.                                                 */
static inline uint64_t clear_prog(uint8_t *seg, uint64_t sp, uint64_t step, uint64_t bend)
{
    if (sp + 7 * step < bend) {
        unsigned ph = (unsigned)(sp & 7);
        uint8_t *p = seg + (sp >> 3);
        uint32_t d[8]; uint8_t m[8];
        for (int t = 0; t < 8; t++) {
            uint64_t bt = ph + (uint64_t)t * step;
            d[t] = (uint32_t)(bt >> 3); m[t] = (uint8_t)~(1u << (bt & 7));
        }
        do {
            p[d[0]] &= m[0]; p[d[1]] &= m[1]; p[d[2]] &= m[2]; p[d[3]] &= m[3];
            p[d[4]] &= m[4]; p[d[5]] &= m[5]; p[d[6]] &= m[6]; p[d[7]] &= m[7];
            p += step; sp += 8 * step;
        } while (sp + 7 * step < bend);
    }
    for (; sp < bend; sp += step) seg[sp >> 3] &= (uint8_t)~(1u << (sp & 7));
    return sp;
}

#ifdef PROFILE
static double PT[8];
static double pnow(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
#define PMARK(i) do { double _t = pnow(); PT[i] += _t - _pt; _pt = _t; } while (0)
#else
#define PMARK(i) do { } while (0)
#endif

static uint64_t tw_run(tw_t *T)
{
    uint64_t total = 0;
    uint8_t *seg = T->seg;
    const uint64_t SB = T->seg_bits;
    const uint64_t nseg = (T->len + T->seg_k - 1) >> T->seg_log;
#ifdef PROFILE
    double _pt = pnow();
#endif

    for (uint64_t s = 0; s < nseg; s++) {
        const uint64_t sbase = s << T->seg_log;          /* chunk-relative k */
        const uint64_t segbit0 = 3 * sbase;
        const uint64_t segend = segbit0 + SB;
        const uint64_t sb0 = segbit0 >> 3;               /* segment byte offset */

        /* 1. presieve 7..97 */
        pat_apply(seg, T->seg_bytes, T->cp0, PS_BYTES, sb0, 0);
        for (int g = 0; g < NPAT; g++) pat_apply(seg, T->seg_bytes, T->cpx[g], PATLEN[g], sb0, 1);
        PMARK(0);

        /* 2. tiny primes, in L1-sized sub-blocks */
        for (uint64_t b0 = 0; b0 < SB; b0 += SUB_BITS) {
            uint64_t bend = b0 + SUB_BITS < SB ? b0 + SUB_BITS : SB;
            for (size_t i = 0; i < T->ntp; i++) {
                sprime_t *S = &T->tp[i];
                const uint64_t step = 3 * (uint64_t)S->q;
                for (int u = 0; u < 6; u++) {
                    uint64_t sp = S->pos[u] - segbit0;
                    if (sp >= bend) continue;
                    S->pos[u] = segbit0 + clear_prog(seg, sp, step, bend);
                }
            }
        }
        PMARK(1);

        /* 3. small primes over the whole segment */
        for (size_t i = 0; i < T->nsp; i++) {
            sprime_t *S = &T->sp[i];
            const uint64_t step = 3 * (uint64_t)S->q;
            for (int u = 0; u < 6; u++) {
                uint64_t pos = S->pos[u];
                if (pos >= segend) continue;
                S->pos[u] = segbit0 + clear_prog(seg, pos - segbit0, step, SB);
            }
        }
        PMARK(2);

        /* 4. medium primes: wheel walk, unrolled one full period (6 hits) */
        for (size_t i = 0; i < T->nmp; i++) {
            mprime_t *M = &T->mp[i];
            uint64_t pos = M->pos;
            if (pos >= segend) continue;
            unsigned u = M->u;
            const int64_t qd = M->qd;
            const int32_t *dbm = DBM[M->wi], *dbc = DBC[M->wi];
            uint64_t E[6], period = 0;
            for (int v = 0; v < 6; v++) {
                E[v] = (uint64_t)((int64_t)dbm[(u + v) % 6] * qd + dbc[(u + v) % 6]);
                period += E[v];
            }
            uint8_t *sg = seg - (segbit0 >> 3);      /* index directly by pos */
            while (pos + period < segend) {
                sg[pos >> 3] &= (uint8_t)~(1u << (pos & 7)); pos += E[0];
                sg[pos >> 3] &= (uint8_t)~(1u << (pos & 7)); pos += E[1];
                sg[pos >> 3] &= (uint8_t)~(1u << (pos & 7)); pos += E[2];
                sg[pos >> 3] &= (uint8_t)~(1u << (pos & 7)); pos += E[3];
                sg[pos >> 3] &= (uint8_t)~(1u << (pos & 7)); pos += E[4];
                sg[pos >> 3] &= (uint8_t)~(1u << (pos & 7)); pos += E[5];
            }
            for (int v = 0; pos < segend; v++) {
                sg[pos >> 3] &= (uint8_t)~(1u << (pos & 7));
                pos += E[v];
                u = (u == 5) ? 0 : u + 1;
            }
            M->pos = pos; M->u = (uint8_t)u;
        }
        PMARK(3);

        /* 5. large primes: bucket for this segment */
        {
            bucket_t *B = &T->bk[s & T->nbm];
            block_t *blk = B->head, *tail = B->tail;
            uint64_t *cur = B->cur;
            B->head = B->tail = NULL; B->cur = B->end = NULL;
            const uint64_t segend_k = sbase + T->seg_k;
            while (blk) {
                const uint64_t *r = blk->rec, *re = (blk == tail) ? cur : blk->rec + BLK;
                for (; r < re; r++) {
                    uint64_t rec = *r;
                    uint32_t a = (uint32_t)rec;
                    const int64_t qd = (int64_t)(rec >> 32);
                    uint64_t k = sbase + (a & 0xFFFFFFu);
                    unsigned j = (a >> 24) & 3, u = (a >> 26) & 7, wi = a >> 29;
                    const int32_t *dkm = DKM[wi], *dkc = DKC[wi];
                    const uint8_t *wj = WJ[wi];
                    do {
                        uint64_t b = 3 * (k - sbase) + j;
                        seg[b >> 3] &= (uint8_t)~(1u << (b & 7));
                        k = (uint64_t)((int64_t)k + (int64_t)dkm[u] * qd + dkc[u]);
                        u = (u == 5) ? 0 : u + 1;
                        j = wj[u];
                    } while (k < segend_k);
                    if (k < T->len)
                        bk_push(T, (k >> T->seg_log) & T->nbm,
                                pack_rec(k & (T->seg_k - 1), j, u, wi, (uint64_t)qd));
                }
                block_t *nx = blk->next;
                blk->next = T->freelist; T->freelist = blk;
                blk = nx;
            }
        }
        PMARK(4);

        /* 6. exact bits for n < 120 (presieve kills 11, 17, 29, ...) */
        if (T->K0 == 0 && s == 0)
            for (uint64_t k = 0; k < 4; k++)
                for (int j = 0; j < 3; j++) {
                    uint64_t n = 30 * k + (uint64_t)R3[j], b = 3 * k + (uint64_t)j;
                    if (is_prime_small(n) && is_prime_small(n + 2)) seg[b >> 3] |= (uint8_t)(1u << (b & 7));
                    else seg[b >> 3] &= (uint8_t)~(1u << (b & 7));
                }

        /* 7. count survivors */
        {
            uint64_t nbits = (s == nseg - 1) ? 3 * (T->len - sbase) : SB;
            const uint64_t *W = (const uint64_t *)seg;
            uint64_t nw = nbits / 64, rem = nbits % 64, c = 0;
            for (uint64_t i = 0; i < nw; i++) c += (uint64_t)__builtin_popcountll(W[i]);
            if (rem) c += (uint64_t)__builtin_popcountll(W[nw] & ((1ULL << rem) - 1));
            total += c;
        }
        PMARK(5);
    }
    return total;
}

static uint64_t isqrt128(u128 x)
{
    uint64_t r = (uint64_t)sqrtl((long double)x);
    while ((u128)r * r > x) r--;
    while ((u128)(r + 1) * (r + 1) <= x) r++;
    return r;
}

/* Count twin pairs (p, p+2) with p = 30k + r, r in {11,17,29}, K0 <= k < K1. */
static uint64_t count_twins_twin(uint64_t K0, uint64_t K1, uint64_t seg_k)
{
    if (K1 <= K0) return 0;
    tw_t T; memset(&T, 0, sizeof T);
    T.K0 = K0; T.K1 = K1; T.len = K1 - K0;
    T.seg_k = seg_k; T.seg_log = __builtin_ctzll(seg_k);
    T.seg_bits = 3 * seg_k; T.seg_bytes = T.seg_bits / 8;
    T.seg = aligned_alloc(64, T.seg_bytes);
    T.cp0 = malloc(PS_BYTES);
    if (!T.seg || !T.cp0) { perror("malloc"); exit(1); }
    pat_align(T.cp0, PAT, PS_BYTES, K0);
    for (int g = 0; g < NPAT; g++) {
        T.cpx[g] = malloc(PATLEN[g]);
        if (!T.cpx[g]) { perror("malloc"); exit(1); }
        pat_align(T.cpx[g], PATB[g], PATLEN[g], K0);
    }

    uint64_t qmax = isqrt128((u128)30 * K1 + 1);
    uint64_t nb = 1;
    while (nb < 2 * (qmax / seg_k) + 4) nb <<= 1;
    T.nb = nb; T.nbm = nb - 1;
    T.bk = calloc(nb, sizeof(bucket_t));
    if (!T.bk) { perror("calloc"); exit(1); }

#ifdef PROFILE
    double _pt = pnow();
#endif
    if (qmax >= 101) {
        primesieve_iterator it;
        primesieve_init(&it);
        primesieve_jump_to(&it, 101, qmax);
        uint64_t q;
        while ((q = primesieve_next_prime(&it)) <= qmax) tw_add_prime(&T, q);
        primesieve_free_iterator(&it);
    }
    PMARK(6);

    uint64_t total = tw_run(&T);
#ifdef PROFILE
    fprintf(stderr, "# profile: setup %.2f  presieve %.2f  tiny %.2f  small %.2f  medium %.2f  large %.2f  count %.2f  (ntp=%zu nsp=%zu nmp=%zu nb=%" PRIu64 ")\n",
            PT[6], PT[0], PT[1], PT[2], PT[3], PT[4], PT[5], T.ntp, T.nsp, T.nmp, T.nb);
    memset(PT, 0, sizeof PT);
#endif

    for (size_t i = 0; i < T.narenas; i++) free(T.arenas[i]);
    free(T.arenas); free(T.bk); free(T.tp); free(T.sp); free(T.mp); free(T.cp0);
    for (int g = 0; g < NPAT; g++) free(T.cpx[g]);
    free(T.seg);
    return total;
}

/* ======================================================================= */
/*  primesieve-based kernels (64-bit only), on the same k-ranges           */
/* ======================================================================= */

/* pairs with p in [lo, hi) via a prime iterator */
static uint64_t count_twins_iter_n(uint64_t lo, uint64_t hi)
{
    primesieve_iterator it;
    primesieve_init(&it);
    primesieve_jump_to(&it, lo, hi + 2);
    uint64_t count = 0;
    uint64_t prev = primesieve_next_prime(&it);
    while (prev < hi) {
        uint64_t cur = primesieve_next_prime(&it);
        if (cur - prev == 2) count++;
        prev = cur;
    }
    primesieve_free_iterator(&it);
    return count;
}

/* k-range -> number range.  k_lo = 0 starts at 7 so (3,5),(5,7) stay out. */
static uint64_t count_twins_ps_k(uint64_t k_lo, uint64_t k_hi, int use_iter)
{
    uint64_t lo = 30 * k_lo, hi = 30 * k_hi;
    if (lo < 7) lo = 7;
    if (hi <= lo) return 0;
    return use_iter ? count_twins_iter_n(lo, hi) : primesieve_count_twins(lo, hi + 1);
}

/* ======================================================================= */
/*  Job control, logging, main                                             */
/* ======================================================================= */

enum { K_TWIN = 0, K_PS = 1, K_ITER = 2 };
static const char *KNAME[3] = {"twin", "ps", "iter"};

typedef struct {
    int n;                          /* decade                                */
    uint64_t K0, K1;                /* absolute k range of the decade         */
    uint64_t chunk_k, nchunks;
    uint64_t first, last;           /* chunk index range [first, last)       */
    _Atomic uint64_t next;
    _Atomic uint64_t total;
    _Atomic uint64_t ndone;         /* chunks finished (resumed + computed)   */
    int kernel;
    uint64_t seg_k;
    uint8_t *done;
    FILE *log;
    pthread_mutex_t mu;
} job_t;

static inline int  bit_get(const uint8_t *b, uint64_t i) { return (b[i >> 3] >> (i & 7)) & 1; }
static inline void bit_set(uint8_t *b, uint64_t i)       { b[i >> 3] |= (uint8_t)(1u << (i & 7)); }

static void *worker(void *arg)
{
    job_t *job = arg;
    uint64_t local = 0, ndone = 0;
    for (;;) {
        uint64_t i = atomic_fetch_add(&job->next, 1);
        if (i >= job->last) break;
        if (job->done && bit_get(job->done, i)) continue;
        uint64_t k_lo = job->K0 + i * job->chunk_k;
        uint64_t k_hi = k_lo + job->chunk_k;
        if (k_hi > job->K1 || k_hi < k_lo) k_hi = job->K1;
        uint64_t c;
        if (job->kernel == K_TWIN) c = count_twins_twin(k_lo, k_hi, job->seg_k);
        else                       c = count_twins_ps_k(k_lo, k_hi, job->kernel == K_ITER);
        local += c; ndone++;
        if (job->log) {
            pthread_mutex_lock(&job->mu);
            fprintf(job->log, "k %d %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", job->n, i, job->chunk_k, c);
            fflush(job->log);
            pthread_mutex_unlock(&job->mu);
        }
    }
    atomic_fetch_add(&job->total, local);
    atomic_fetch_add(&job->ndone, ndone);
    return NULL;
}

static double now(void)
{
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static u128 pow10_128(int n) { u128 r = 1; while (n-- > 0) r *= 10; return r; }

/* K_n = (10^n - 10)/30: number of k-blocks below 10^n.  K_0 = K_1 = 0. */
static uint64_t Kbound(int n)
{
    if (n < 1) return 0;
    return (uint64_t)((pow10_128(n) - 10) / 30);
}

static uint64_t default_chunk_k(int n)
{
    int e = n - 4; if (e < 9) e = 9; if (e > 13) e = 13;
    return (uint64_t)(pow10_128(e) / 30);
}

/* If the log does not end in a newline, cut it back to the last full line. */
static void trim_partial_line(const char *path)
{
    int fd = open(path, O_RDWR);
    if (fd < 0) return;
    struct stat st;
    if (fstat(fd, &st) != 0 || st.st_size == 0) { close(fd); return; }
    char c;
    if (pread(fd, &c, 1, st.st_size - 1) == 1 && c == '\n') { close(fd); return; }
    off_t keep = 0, pos = st.st_size;
    char buf[4096];
    while (pos > 0 && keep == 0) {
        size_t len = (size_t)(pos < (off_t)sizeof buf ? pos : (off_t)sizeof buf);
        pos -= (off_t)len;
        if (pread(fd, buf, len, pos) != (ssize_t)len) break;
        for (size_t i = len; i-- > 0; )
            if (buf[i] == '\n') { keep = pos + (off_t)i + 1; break; }
    }
    fprintf(stderr, "# log %s: dropping %lld bytes of partial last line\n", path, (long long)(st.st_size - keep));
    if (ftruncate(fd, keep) != 0) perror("ftruncate");
    close(fd);
}

static uint64_t load_chunks(const char *path, job_t *job)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    uint64_t resumed = 0;
    char line[256];
    while (fgets(line, sizeof line, f)) {
        int n; uint64_t idx, ck, cnt;
        if (!strchr(line, '\n')) continue;
        if (sscanf(line, "k %d %" SCNu64 " %" SCNu64 " %" SCNu64, &n, &idx, &ck, &cnt) != 4 || n != job->n) continue;
        if (ck != job->chunk_k) {
            fprintf(stderr, "log %s: decade %d was run with chunk %" PRIu64 " k-units, now %" PRIu64 "\n", path, n, ck, job->chunk_k);
            exit(1);
        }
        if (idx >= job->nchunks || bit_get(job->done, idx)) continue;
        bit_set(job->done, idx);
        atomic_fetch_add(&job->total, cnt);
        resumed++;
    }
    fclose(f);
    return resumed;
}

static int load_decades(const char *path, uint64_t *a_out)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    int best = 0; uint64_t a = 0; char line[256];
    while (fgets(line, sizeof line, f)) {
        int n; uint64_t v;
        if (!strchr(line, '\n')) continue;
        if (sscanf(line, "d %d %" SCNu64, &n, &v) == 2 && n > best) { best = n; a = v; }
    }
    fclose(f);
    *a_out = a;
    return best;
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "usage: %s [options] N\n"
        "  N            compute a(1)..a(N), twin prime pairs below 10^n (N <= 20)\n"
        "  -t threads   worker threads (default: online CPUs)\n"
        "  -k kernel    twin (default, works to 10^20) | ps | iter (primesieve, <= 10^19)\n"
        "  -c           same as -k ps\n"
        "  -s chunk     chunk size in numbers (default 10^max(9,n-4), capped at 10^13)\n"
        "  -S log2      twin-kernel segment size, log2 of k-units (16..24, default 22)\n"
        "  -f k         start at decade k\n"
        "  -a A         known a(k-1), used with -f\n"
        "  -l file      append-only chunk log; rerun with the same log to resume\n"
        "  -r A:B       only chunks A <= idx < B of a single decade (split across machines)\n", prog);
    exit(2);
}

int main(int argc, char **argv)
{
    int nthreads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    int kernel = K_TWIN, seg_log = 22, from = 0, have_a = 0, have_r = 0;
    uint64_t chunk_numbers = 0, a_known = 0, r_first = 0, r_last = UINT64_MAX;
    const char *logpath = NULL;
    int opt;

    while ((opt = getopt(argc, argv, "t:k:cs:S:f:a:l:r:h")) != -1) {
        switch (opt) {
        case 't': nthreads = atoi(optarg); break;
        case 'k':
            if (!strcmp(optarg, "twin")) kernel = K_TWIN;
            else if (!strcmp(optarg, "ps")) kernel = K_PS;
            else if (!strcmp(optarg, "iter")) kernel = K_ITER;
            else usage(argv[0]);
            break;
        case 'c': kernel = K_PS; break;
        case 's': chunk_numbers = strtoull(optarg, NULL, 10); break;
        case 'S': seg_log = atoi(optarg); break;
        case 'f': from = atoi(optarg); break;
        case 'a': a_known = strtoull(optarg, NULL, 10); have_a = 1; break;
        case 'l': logpath = optarg; break;
        case 'r':
            if (sscanf(optarg, "%" SCNu64 ":%" SCNu64, &r_first, &r_last) != 2 || r_last <= r_first) usage(argv[0]);
            have_r = 1; break;
        default: usage(argv[0]);
        }
    }
    if (optind != argc - 1) usage(argv[0]);
    int N = atoi(argv[optind]);
    if (N < 1 || N > 20) { fprintf(stderr, "N must be 1..20\n"); return 2; }
    if (nthreads < 1) nthreads = 1;
    if (seg_log < 16 || seg_log > 24) { fprintf(stderr, "-S must be 16..24\n"); return 2; }
    uint64_t seg_k = 1ULL << seg_log;

    wheel_init();
    presieve_init();
    presieve2_init();
    primesieve_set_num_threads(1);

    FILE *logf = NULL;
    if (logpath) {
        trim_partial_line(logpath);
        logf = fopen(logpath, "a");
        if (!logf) { perror(logpath); return 1; }
    }

    uint64_t cumulative = 0;
    int start = 1;
    if (from > 0) {
        if (from > 1 && !have_a) {
            uint64_t v; int k = logpath ? load_decades(logpath, &v) : 0;
            if (k == from - 1) { a_known = v; have_a = 1; }
            else { fprintf(stderr, "-f %d needs -a a(%d)\n", from, from - 1); return 2; }
        }
        start = from;
        cumulative = (from > 1) ? a_known : 0;
    } else if (logpath) {
        uint64_t v; int k = load_decades(logpath, &v);
        if (k > 0) { start = k + 1; cumulative = v; }
    }
    if (start > N) {
        fprintf(stderr, "nothing to do: a(%d) = %" PRIu64 " already known\n", start - 1, cumulative);
        if (logf) fclose(logf);
        return 0;
    }
    if (have_r && start != N) {
        fprintf(stderr, "-r needs a single decade: use -f %d (with -a) and N = %d\n", N, N);
        return 2;
    }

    pthread_t *tid = malloc(sizeof(pthread_t) * (size_t)nthreads);
    if (!tid) { perror("malloc"); return 1; }

    printf("# A007508: twin prime pairs below 10^n  (%d threads, kernel=%s, seg=2^%d k)\n", nthreads, KNAME[kernel], seg_log);
    printf("# n  a(n)  [decade count]  [seconds]\n");
    if (start > 1) printf("# resuming at decade %d with a(%d) = %" PRIu64 "\n", start, start - 1, cumulative);
    fflush(stdout);
    double t_start = now();

    for (int n = start; n <= N; n++) {
        job_t job; memset(&job, 0, sizeof job);
        job.n = n; job.K0 = Kbound(n - 1); job.K1 = Kbound(n);
        job.chunk_k = chunk_numbers ? chunk_numbers / 30 : default_chunk_k(n);
        if (job.chunk_k < 64) job.chunk_k = 64;
        uint64_t span = job.K1 - job.K0;
        job.nchunks = (span + job.chunk_k - 1) / job.chunk_k;
        job.first = 0; job.last = job.nchunks;
        if (have_r) {
            job.first = r_first < job.nchunks ? r_first : job.nchunks;
            job.last  = r_last  < job.nchunks ? r_last  : job.nchunks;
        }
        if (kernel != K_TWIN) {
            uint64_t top = job.last ? job.K0 + job.last * job.chunk_k : job.K0;
            if (top > job.K1) top = job.K1;
            if (top > UINT64_MAX / 30) {
                fprintf(stderr, "chunk range reaches beyond 2^64: only -k twin can do it\n");
                return 2;
            }
        }
        atomic_init(&job.next, job.first);
        atomic_init(&job.total, 0);
        atomic_init(&job.ndone, 0);
        job.kernel = kernel; job.seg_k = seg_k; job.log = logf;
        pthread_mutex_init(&job.mu, NULL);

        uint64_t resumed = 0;
        if (logpath) {
            job.done = calloc((job.nchunks + 7) / 8 + 1, 1);
            if (!job.done) { perror("calloc"); return 1; }
            resumed = load_chunks(logpath, &job);
            if (resumed) printf("# decade %d: %" PRIu64 " of %" PRIu64 " chunks resumed from %s\n", n, resumed, job.nchunks, logpath);
        }

        double t0 = now();
        int nt = nthreads;
        uint64_t todo = job.last > job.first ? job.last - job.first : 0;
        if ((uint64_t)nt > todo) nt = (int)todo;
        for (int i = 0; i < nt; i++)
            if (pthread_create(&tid[i], NULL, worker, &job) != 0) { perror("pthread_create"); return 1; }
        for (int i = 0; i < nt; i++) pthread_join(tid[i], NULL);

        uint64_t dec = atomic_load(&job.total) + (n == 1 ? 2 : 0);   /* (3,5),(5,7) */
        int complete = (resumed + atomic_load(&job.ndone)) == job.nchunks;
        if (complete) {
            cumulative += dec;
            printf("%2d %" PRIu64 "  [%" PRIu64 "]  [%.2f]\n", n, cumulative, dec, now() - t0);
            if (logf) { fprintf(logf, "d %d %" PRIu64 "\n", n, cumulative); fflush(logf); }
        } else {
            printf("%2d partial  [%" PRIu64 " in chunks %" PRIu64 "..%" PRIu64 " of %" PRIu64 "; %" PRIu64 " chunks logged]  [%.2f]\n",
                   n, dec, job.first, job.last, job.nchunks, resumed + atomic_load(&job.ndone), now() - t0);
        }
        fflush(stdout);
        free(job.done);
        pthread_mutex_destroy(&job.mu);
    }
    fprintf(stderr, "# total time %.2f s\n", now() - t_start);
    if (logf) fclose(logf);
    free(tid);
    return 0;
}
