/*
 * gapsearch.c -- search for prime chains following an arbitrary gap pattern.
 *
 *   --gaps 2^n   ->  2,4,8,16,...      p, p+2, p+6, p+14, ...
 *   --gaps 2n    ->  2,4,6,8,...       p, p+2, p+6, p+12, ...
 *   --gaps 2p    ->  4,6,10,14,22,...  p, p+4, p+10, p+20, ...  (2 * prime)
 *   --gaps 2p@3  ->  6,10,14,22,...    (2 * prime, starting from 3)
 *   --gaps 6,10,4,2,...                explicit list
 *
 * Before searching, the program computes the largest length for which the
 * pattern is ADMISSIBLE -- i.e. for which no small prime q covers every
 * residue class. If the offsets hit all q residues mod q, then one chain
 * member is always divisible by q and no chains of that length exist at all
 * (beyond the trivial ones where a member IS q). That check is cheap and it
 * saves searching for something that provably is not there.
 *
 * Method: constellation (k-tuple) sieve in a compressed wheel space. Wheel
 * residues surviving the pattern are enumerated by CRT, each small prime then
 * strikes one arithmetic progression per offset per residue row, and only the
 * survivors reach Miller-Rabin.
 *
 * Build:  gcc   -O3 -march=native -o gapsearch gapsearch.c -lpthread -lm
 *         clang -O3 -mcpu=native  -o gapsearch gapsearch.c -lpthread -lm
 */

#if !defined(__APPLE__)          /* Apple hides sysconf extensions in strict POSIX mode */
#define _POSIX_C_SOURCE 200809L
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdatomic.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#if defined(__APPLE__) || defined(__FreeBSD__)
#include <sys/sysctl.h>
#endif

#define MAXLEN 32

/* ------------------------------------------------------------------ */
/* globals (read-only after setup)                                     */
/* ------------------------------------------------------------------ */

static uint64_t OFF[MAXLEN + 1];      /* cumulative offsets, OFF[0] = 0 */
static uint64_t GAPS[MAXLEN];
static int      NOFF;                 /* offsets available = max chain length */
static int      ADMISSIBLE;           /* largest provably-possible length */
static char     GAPSPEC[128] = "2^n";
static int      MINLEN     = 6;
static int      REQ_CONSEC = 0;
static int      BRUTE      = 0;
static int      QUIET      = 0;
static uint64_t RSTART, REND;
static uint32_t SIEVE_B;
static int      NTHREADS   = 1;
static uint32_t JLEN;                 /* bits per residue class per block */
static size_t   BM_TARGET  = 4u << 20;/* bitmap bytes per thread */

/* wheel */
static uint64_t  WMOD;
static uint32_t  NRES;
static uint32_t *RES;

/* sieving primes (all > largest wheel prime, <= SIEVE_B) */
static uint32_t  NSP;
static uint32_t *SP;        /* prime            */
static uint32_t *SPINV;     /* W^-1 mod p       */
static uint32_t *DK;        /* NSP * MAXLEN: distinct (-off*W^-1) mod p */
static uint8_t  *DKN;       /* how many are distinct                   */

static uint64_t SPAN;       /* integers covered by one block */
static uint64_t BLK_FIRST, BLK_LAST;
static uint64_t SIEVE_LO, SIEVE_HI;

static atomic_ullong next_block;
static atomic_ullong stat_survivors, stat_chains, stat_blocks;

/* ------------------------------------------------------------------ */
/* number parsing: 12, 1.5T, 250B, 3e15, 1_000_000                      */
/* ------------------------------------------------------------------ */

static int is_small_prime(uint64_t n)
{
    if (n < 2) return 0;
    for (uint64_t d = 2; d * d <= n; d++) if (n % d == 0) return 0;
    return 1;
}

/* fill GAPS/OFF from a pattern spec; returns number of offsets */
static void gen_offsets(const char *spec)
{
    int ng = 0;
    snprintf(GAPSPEC, sizeof GAPSPEC, "%s", spec);

    static const char *KT[19] = { 0,0,
        "2","2,4","2,4,2","2,4,2,4","4,2,4,2,4","2,4,2,4,6,2","2,4,2,4,6,2,6",
        "2,4,2,4,6,2,6,4","2,4,2,4,6,2,6,4,2","2,4,2,4,6,2,6,4,2,4",
        "2,4,2,4,6,2,6,4,2,4,6","2,4,2,4,6,2,6,4,2,4,6,6",
        "2,4,2,4,6,2,6,4,2,4,6,6,2","2,4,2,4,6,2,6,4,2,4,6,6,2,6",
        "2,4,6,2,6,6,4,2,4,6,2,6,4,2,4","2,4,6,2,6,4,2,4,6,6,2,6,4,2,6,4",
        "4,2,4,6,2,6,4,2,4,6,6,2,6,4,2,6,4" };

    if (!strncmp(spec, "kt", 2)) {                  /* densest k-tuple */
        int k = atoi(spec + 2);
        if (k < 2 || k > 18) { fprintf(stderr, "kt: k must be 2..18\n"); exit(1); }
        gen_offsets(KT[k]);
        snprintf(GAPSPEC, sizeof GAPSPEC, "kt%d", k);
        return;
    }
    if (strchr(spec, 'x') && isdigit((unsigned char)spec[0])) {  /* e.g. 30x6 */
        char *e;
        unsigned long long g = strtoull(spec, &e, 10);
        int rep = atoi(e + 1);
        if (!g || rep < 1 || rep > MAXLEN - 1) { fprintf(stderr, "bad NxM spec\n"); exit(1); }
        while (ng < rep) GAPS[ng++] = g;
    } else if (strlen(spec) > 2 && !strcmp(spec + strlen(spec) - 2, "^n")
               && isdigit((unsigned char)spec[0])) {      /* r^n -> 2*r^k */
        uint64_t r = strtoull(spec, NULL, 10);
        if (r < 2) { fprintf(stderr, "bad r^n spec\n"); exit(1); }
        uint64_t v = 2;
        while (ng < MAXLEN - 1 && v < (1ULL << 40)) { GAPS[ng++] = v; v *= r; }
    } else if (!strcmp(spec, "2n")) {
        while (ng < MAXLEN - 1) { GAPS[ng] = 2 * (uint64_t)(ng + 1); ng++; }
    } else if (spec[0] == 'p' || (spec[0] == '2' && spec[1] == 'p')) {
        uint64_t mul = spec[0] == '2' ? 2 : 1;          /* 2p / p, opt. @start */
        const char *at = strchr(spec, '@');
        uint64_t from = at ? strtoull(at + 1, NULL, 10) : 2;
        if (from < 2) from = 2;
        for (uint64_t q = from; ng < MAXLEN - 1; q++)
            if (is_small_prime(q)) GAPS[ng++] = mul * q;
    } else if (strspn(spec, "0123456789,") == strlen(spec) && isdigit((unsigned char)spec[0])) {
        const char *p = spec;                            /* explicit gap list */
        while (*p && ng < MAXLEN - 1) {
            char *e;
            unsigned long long v = strtoull(p, &e, 10);
            if (e == p || v == 0) break;
            GAPS[ng++] = v;
            if (*e != ',') { p = e; break; }
            p = e + 1;
        }
        if (*p || ng == 0) { fprintf(stderr, "bad --gaps list: %s\n", spec); exit(1); }
    } else {
        fprintf(stderr, "unknown --gaps spec: %s\n"
                        "  try: 2^n | 2n | 2p | 2p@3 | p | 4,6,10,14\n", spec);
        exit(1);
    }

    OFF[0] = 0;
    for (int k = 0; k < ng; k++) {
        if (GAPS[k] > (1ULL << 40) || OFF[k] > (1ULL << 40)) { ng = k; break; }
        OFF[k + 1] = OFF[k] + GAPS[k];
    }
    NOFF = ng + 1;
}

/* largest L for which no prime q <= L covers all residues mod q */
static int max_admissible(int *bad_q)
{
    for (int L = 2; L <= NOFF; L++)
        for (int q = 2; q <= L && q <= MAXLEN; q++) {
            if (!is_small_prime(q)) continue;
            int seen[MAXLEN + 1] = {0}, s = 0;
            for (int k = 0; k < L; k++) {
                int r = (int)(OFF[k] % q);
                if (!seen[r]) { seen[r] = 1; s++; }
            }
            if (s == q) { *bad_q = q; return L - 1; }
        }
    *bad_q = 0;
    return NOFF;
}

static int cpu_count(void)
{
#if defined(__APPLE__) || defined(__FreeBSD__)
    int n = 0; size_t sz = sizeof n;
    if (sysctlbyname("hw.logicalcpu", &n, &sz, NULL, 0) == 0 && n > 0) return n;
#endif
#ifdef _SC_NPROCESSORS_ONLN
    long m = sysconf(_SC_NPROCESSORS_ONLN);
    if (m > 0) return (int)m;
#endif
    return 1;
}

static int parse_qty(const char *s, uint64_t *out)
{
    char buf[160];
    size_t n = 0;
    for (const char *p = s; *p && n < sizeof buf - 1; ++p)
        if (*p != ',' && *p != '_' && !isspace((unsigned char)*p)) buf[n++] = *p;
    buf[n] = 0;
    if (!n) return -1;

    unsigned long long mult = 1;
    size_t len = strlen(buf);
    switch (toupper((unsigned char)buf[len - 1])) {
        case 'K': mult = 1000ULL;                   break;
        case 'M': mult = 1000000ULL;                break;
        case 'G': case 'B': mult = 1000000000ULL;   break;
        case 'T': mult = 1000000000000ULL;          break;
        case 'P': mult = 1000000000000000ULL;       break;
        case 'E': mult = 1000000000000000000ULL;    break;
        default:  mult = 0;
    }
    if (mult) buf[--len] = 0; else mult = 1;
    if (!len) return -1;

    char *endp;
    if (strpbrk(buf, ".eE")) {                 /* fractional / scientific */
        errno = 0;
        long double v = strtold(buf, &endp);
        if (*endp || errno || v < 0) return -1;
        long double r = v * (long double)mult;
        if (r >= 1.8446744073709551615e19L) return -1;
        *out = (uint64_t)r;
    } else {
        errno = 0;
        unsigned long long v = strtoull(buf, &endp, 10);
        if (*endp || errno) return -1;
        unsigned __int128 r = (unsigned __int128)v * mult;
        if (r > UINT64_MAX) return -1;
        *out = (uint64_t)r;
    }
    return 0;
}

static void fmt_u64(uint64_t v, char *out)      /* 1234567 -> 1,234,567 */
{
    char t[24]; int n = snprintf(t, sizeof t, "%llu", (unsigned long long)v);
    int o = 0;
    for (int i = 0; i < n; i++) {
        if (i && (n - i) % 3 == 0) out[o++] = ',';
        out[o++] = t[i];
    }
    out[o] = 0;
}

/* ------------------------------------------------------------------ */
/* primality: deterministic Miller-Rabin for the full 64-bit range      */
/* ------------------------------------------------------------------ */

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m)
{ return (uint64_t)((__uint128_t)a * b % m); }

static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m)
{
    uint64_t r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}

static const uint32_t TRIALP[] = {
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
    101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
    193,197,199,211,223,227,229,233,239,241,251
};
#define NTRIALP (sizeof TRIALP / sizeof TRIALP[0])

static int is_prime_u64(uint64_t n)
{
    if (n < 2) return 0;
    for (size_t i = 0; i < NTRIALP; i++) {
        uint32_t p = TRIALP[i];
        if ((uint64_t)p * p > n) return 1;
        if (n % p == 0) return n == p;
    }
    uint64_t d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    static const uint64_t base[] = {2,325,9375,28178,450775,9780504,1795265022};
    for (int i = 0; i < 7; i++) {
        uint64_t a = base[i] % n;
        if (!a) continue;
        uint64_t x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int r = 1; r < s; r++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { ok = 1; break; }
        }
        if (!ok) return 0;
    }
    return 1;
}

/* how many chain members starting at c are prime (stops at first failure) */
static int chain_len_at(uint64_t c)
{
    int n = 0;
    for (int k = 0; k < NOFF; k++) {
        if (OFF[k] > UINT64_MAX - c) break;
        if (!is_prime_u64(c + OFF[k])) break;
        n++;
    }
    return n;
}

/* longest prefix whose members are *consecutive* primes (nothing between) */
static int consec_len(uint64_t c, int n)
{
    for (int k = 0; k + 1 < n; k++) {
        uint64_t hi = c + OFF[k + 1];
        for (uint64_t x = c + OFF[k] + 2; x < hi; x += 2)
            if (is_prime_u64(x)) return k + 1;
    }
    return n;
}

/* ------------------------------------------------------------------ */
/* small prime table                                                    */
/* ------------------------------------------------------------------ */

static uint32_t *make_primes(uint32_t lim, uint32_t *cnt)
{
    uint8_t *c = calloc(lim + 1, 1);
    if (!c) { perror("calloc"); exit(1); }
    for (uint64_t i = 2; i * i <= lim; i++)
        if (!c[i]) for (uint64_t j = i * i; j <= lim; j += i) c[j] = 1;
    uint32_t n = 0;
    for (uint32_t i = 2; i <= lim; i++) if (!c[i]) n++;
    uint32_t *p = malloc((size_t)n * sizeof *p);
    n = 0;
    for (uint32_t i = 2; i <= lim; i++) if (!c[i]) p[n++] = i;
    free(c);
    *cnt = n;
    return p;
}

/* ------------------------------------------------------------------ */
/* wheel construction                                                   */
/* ------------------------------------------------------------------ */

static uint32_t max_wheel_prime;

static int cmp_u32(const void *a, const void *b)
{
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return (x > y) - (x < y);
}

static void build_wheel(const uint32_t *pr, uint32_t npr,
                        uint32_t maxres, uint64_t maxW)
{
    uint64_t  M   = 1;
    uint32_t  cnt = 1;
    uint32_t *cur = malloc(sizeof *cur);
    cur[0] = 0;
    max_wheel_prime = 1;

    for (uint32_t pi = 0; pi < npr; pi++) {
        uint32_t q = pr[pi];
        if (M > maxW / q) break;

        /* which residues mod q survive the pattern? */
        uint8_t *bad = calloc(q, 1);
        for (int k = 0; k < MINLEN; k++) {
            uint32_t r = (uint32_t)(OFF[k] % q);      /* kill r == -off */
            bad[(q - r) % q] = 1;
        }
        uint32_t good = 0;
        for (uint32_t r = 0; r < q; r++) if (!bad[r]) good++;
        if (good == 0) {
            fprintf(stderr, "pattern is inadmissible mod %u -- no chains exist\n", q);
            exit(1);
        }
        if ((uint64_t)cnt * good > maxres) { free(bad); break; }

        uint32_t *nw = malloc((size_t)cnt * good * sizeof *nw);
        uint32_t  nn = 0;
        for (uint32_t t = 0; t < q; t++)
            for (uint32_t i = 0; i < cnt; i++) {
                uint64_t v = cur[i] + M * t;
                if (!bad[v % q]) nw[nn++] = (uint32_t)v;
            }
        free(bad); free(cur);
        cur = nw; cnt = nn; M *= q;
        max_wheel_prime = q;
    }

    qsort(cur, cnt, sizeof *cur, cmp_u32);

    WMOD = M; NRES = cnt; RES = cur;
}

/* ------------------------------------------------------------------ */
/* per-prime setup: W^-1 mod p and the distinct offset shifts           */
/* ------------------------------------------------------------------ */

static void build_prime_tables(const uint32_t *pr, uint32_t npr)
{
    NSP = 0;
    for (uint32_t i = 0; i < npr; i++) if (pr[i] > max_wheel_prime) NSP++;

    SP    = malloc((size_t)NSP * sizeof *SP);
    SPINV = malloc((size_t)NSP * sizeof *SPINV);
    DK    = malloc((size_t)NSP * MINLEN * sizeof *DK);
    DKN   = malloc((size_t)NSP);

    uint32_t n = 0;
    for (uint32_t i = 0; i < npr; i++) {
        uint32_t q = pr[i];
        if (q <= max_wheel_prime) continue;
        uint32_t wm  = (uint32_t)(WMOD % q);
        uint32_t inv = (uint32_t)powmod(wm, q - 2, q);   /* q is prime */
        SP[n] = q; SPINV[n] = inv;

        uint32_t *d = DK + (size_t)n * MINLEN;
        uint8_t   c = 0;
        for (int k = 0; k < MINLEN; k++) {
            uint32_t r  = (uint32_t)(OFF[k] % q);
            uint32_t nr = r ? q - r : 0;                       /* -off mod q */
            uint32_t v  = (uint32_t)((uint64_t)nr * inv % q);  /* * W^-1     */
            int dup = 0;
            for (uint8_t j = 0; j < c; j++) if (d[j] == v) { dup = 1; break; }
            if (!dup) d[c++] = v;
        }
        DKN[n] = c;
        n++;
    }
}

/* ------------------------------------------------------------------ */
/* singular series + expected-count estimate (informational)            */
/* ------------------------------------------------------------------ */

static double singular_series(const uint32_t *pr, uint32_t npr, int L)
{
    double S = 1.0;
    for (uint32_t i = 0; i < npr; i++) {
        uint32_t q = pr[i];
        uint8_t *seen = calloc(q, 1);
        uint32_t s = 0;
        for (int k = 0; k < L; k++) {
            uint32_t r = (uint32_t)(OFF[k] % q);
            if (!seen[r]) { seen[r] = 1; s++; }
        }
        free(seen);
        S *= (1.0 - (double)s / q) / pow(1.0 - 1.0 / q, L);
    }
    return S;
}

static double expected_count(double a, double b, int L, double S)
{
    if (a < 1000) a = 1000;          /* the heuristic is meaningless below this */
    if (b <= a) return 0;
    /* integrate in log space: x = e^t, dx = e^t dt -- uniform steps in x
       badly overweight the small-x end where 1/(ln x)^L blows up */
    const int N = 4000;
    double la = log(a), lb = log(b), h = (lb - la) / N, sum = 0;
    for (int i = 0; i <= N; i++) {
        double t = la + h * i;
        double f = exp(t) / pow(t, L);
        sum += (i == 0 || i == N) ? f : (i & 1 ? 4 * f : 2 * f);
    }
    return S * sum * h / 3.0;
}

/* ------------------------------------------------------------------ */
/* ordered output                                                       */
/* ------------------------------------------------------------------ */

typedef struct Res { uint64_t blk; char *txt; struct Res *next; } Res;
static Res *pending;
static uint64_t next_emit;
static pthread_mutex_t outmx = PTHREAD_MUTEX_INITIALIZER;

static void submit(uint64_t blk, char *txt)
{
    pthread_mutex_lock(&outmx);
    Res *nd = malloc(sizeof *nd);
    nd->blk = blk; nd->txt = txt; nd->next = pending; pending = nd;
    for (;;) {
        Res **pp = &pending, *hit = NULL;
        while (*pp) {
            if ((*pp)->blk == next_emit) { hit = *pp; *pp = hit->next; break; }
            pp = &(*pp)->next;
        }
        if (!hit) break;
        if (hit->txt && *hit->txt) { fputs(hit->txt, stdout); fflush(stdout); }
        free(hit->txt); free(hit);
        next_emit++;
    }
    pthread_mutex_unlock(&outmx);
}

/* ------------------------------------------------------------------ */
/* reporting a hit                                                      */
/* ------------------------------------------------------------------ */

typedef struct { char *b; size_t n, cap; } sbuf;

static void sb_add(sbuf *s, const char *fmt, ...)
{
    va_list ap;
    for (;;) {
        size_t room = s->cap - s->n;
        va_start(ap, fmt);
        int need = vsnprintf(s->b + s->n, room, fmt, ap);
        va_end(ap);
        if (need >= 0 && (size_t)need < room) { s->n += need; return; }
        s->cap = s->cap ? s->cap * 2 : 4096;
        if (s->cap < s->n + (size_t)need + 1) s->cap = s->n + need + 1;
        s->b = realloc(s->b, s->cap);
    }
}

static void report(sbuf *s, uint64_t c, int n, int m)
{
    char num[32]; fmt_u64(c, num);
    if (REQ_CONSEC)
        sb_add(s, "len %2d (consec %2d)  p = %s\n", n, m, num);
    else
        sb_add(s, "len %2d  p = %s\n", n, num);
    sb_add(s, "        ");
    for (int k = 0; k < n; k++)
        sb_add(s, "%llu%s", (unsigned long long)(c + OFF[k]), k + 1 < n ? " " : "\n");
}

/* ------------------------------------------------------------------ */
/* worker                                                               */
/* ------------------------------------------------------------------ */

typedef struct {
    uint64_t *bm;
    size_t    words;      /* per residue class */
} worker_t;

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static void *worker(void *arg)
{
    worker_t *w = arg;
    const uint32_t stride = (uint32_t)w->words;
    const int      L      = MINLEN;

    uint64_t *hits = NULL; size_t nhits, caphits = 0;

    for (;;) {
        uint64_t blk = atomic_fetch_add(&next_block, 1);
        if (blk > BLK_LAST) break;
        uint64_t base = blk * SPAN;

        memset(w->bm, 0, (size_t)NRES * stride * sizeof *w->bm);

        /* ---- constellation sieve ------------------------------------ */
        for (uint32_t pi = 0; pi < NSP; pi++) {
            const uint32_t  q    = SP[pi];
            const uint64_t  inv  = SPINV[pi];
            const uint32_t *dk   = DK + (size_t)pi * L;
            const uint8_t   nd   = DKN[pi];
            const uint32_t  bmod = (uint32_t)(base % q);

            for (uint32_t i = 0; i < NRES; i++) {
                uint32_t rq = (uint32_t)(((uint64_t)bmod + RES[i]) % q);
                uint32_t nb = rq ? q - rq : 0;
                uint32_t b  = (uint32_t)((uint64_t)nb * inv % q);
                uint64_t *row = w->bm + (size_t)i * stride;
                for (uint8_t k = 0; k < nd; k++) {
                    uint32_t j = b + dk[k];
                    if (j >= q) j -= q;
                    for (; j < JLEN; j += q) row[j >> 6] |= 1ULL << (j & 63);
                }
            }
        }

        /* ---- collect survivors and confirm -------------------------- */
        nhits = 0;
        uint64_t surv = 0;
        for (uint32_t i = 0; i < NRES; i++) {
            const uint64_t *row = w->bm + (size_t)i * stride;
            uint64_t cbase = base + RES[i];
            for (uint32_t wi = 0; wi < stride; wi++) {
                uint64_t m = ~row[wi];
                while (m) {
                    uint32_t j = (uint32_t)(wi * 64 + __builtin_ctzll(m));
                    m &= m - 1;
                    uint64_t c = cbase + WMOD * (uint64_t)j;
                    if (c < SIEVE_LO || c > SIEVE_HI) continue;
                    surv++;
                    int n = chain_len_at(c);
                    if (n < L) continue;
                    int mm = n;
                    if (REQ_CONSEC) {
                        mm = consec_len(c, n);
                        if (mm < L) continue;
                    }
                    if (nhits == caphits) {
                        caphits = caphits ? caphits * 2 : 64;
                        hits = realloc(hits, caphits * 3 * sizeof *hits);
                    }
                    hits[nhits * 3 + 0] = c;
                    hits[nhits * 3 + 1] = n;
                    hits[nhits * 3 + 2] = mm;
                    nhits++;
                }
            }
        }

        /* sort hits by p (classes are scanned out of order) */
        if (nhits > 1) {
            for (size_t a = 1; a < nhits; a++) {          /* insertion sort */
                uint64_t t0 = hits[a*3], t1 = hits[a*3+1], t2 = hits[a*3+2];
                size_t b = a;
                while (b && hits[(b-1)*3] > t0) {
                    memcpy(hits + b*3, hits + (b-1)*3, 3 * sizeof *hits);
                    b--;
                }
                hits[b*3] = t0; hits[b*3+1] = t1; hits[b*3+2] = t2;
            }
        }

        sbuf s = {0};
        for (size_t a = 0; a < nhits; a++)
            report(&s, hits[a*3], (int)hits[a*3+1], (int)hits[a*3+2]);

        atomic_fetch_add(&stat_survivors, surv);
        atomic_fetch_add(&stat_chains, nhits);
        atomic_fetch_add(&stat_blocks, 1);
        submit(blk, s.b);
    }
    free(hits);
    (void)cmp_u64;
    return NULL;
}

/* ------------------------------------------------------------------ */

static void brute_range(uint64_t lo, uint64_t hi)
{
    if (lo > hi) return;
    if (lo < 2) lo = 2;
    sbuf s = {0};
    for (uint64_t c = lo; c <= hi; c++) {
        if (c > 2 && !(c & 1)) continue;
        int n = chain_len_at(c);
        if (n < MINLEN) continue;
        int m = n;
        if (REQ_CONSEC) { m = consec_len(c, n); if (m < MINLEN) continue; }
        report(&s, c, n, m);
        atomic_fetch_add(&stat_chains, 1);
        if (c == hi) break;
    }
    if (s.b) { fputs(s.b, stdout); fflush(stdout); free(s.b); }
}

static void usage(const char *p)
{
    fprintf(stderr,
"usage: %s [options] <start> <end>\n"
"\n"
"  Finds primes p such that p, p+g1, p+g1+g2, ... are all prime, for a\n"
"  chosen gap pattern.\n"
"\n"
"  start/end accept K M B/G T P E suffixes and decimals: 250B, 1.5T, 3e15\n"
"\n"
"  -g, --gaps SPEC    gap pattern (default 2^n):\n"
"                       2^n        2,4,8,16,32,...   (doubling)\n"
"                       2n         2,4,6,8,10,...\n"
"                       2p         4,6,10,14,22,...  (twice each prime)\n"
"                       2p@3       6,10,14,22,...    (twice each prime >= 3)\n"
"                       p          2,3,5,7,...\n"
"                       3^n        2,6,18,54,...     (2*r^k for any r)\n"
"                       30x6       six gaps of 30    (AP / CPAP)\n"
"                       kt12       densest admissible 12-tuple (kt2..kt18)\n"
"                       4,6,10,14  explicit list\n"
"  -l, --len N        minimum chain length to report (default 6)\n"
"  -t, --threads N    worker threads (default: online CPUs)\n"
"  -c, --consecutive  require the members to be CONSECUTIVE primes\n"
"                     (nothing prime in between -- see notes in README)\n"
"  -b, --bound N      small-prime sieve bound (default: auto)\n"
"  -m, --mem N        bitmap bytes per thread (default 4M)\n"
"      --brute        naive verification mode (no sieve) -- for testing\n"
"  -q, --quiet        suppress the header/progress lines\n"
"  -h, --help         this message\n", p);
}

int main(int argc, char **argv)
{
    uint64_t args[2]; int nargs = 0;
    NTHREADS = cpu_count();
    uint32_t user_bound = 0;
    const char *gapspec = "2^n";

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        #define NEXT() (i + 1 < argc ? argv[++i] : (usage(argv[0]), exit(1), (char*)0))
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) { usage(argv[0]); return 0; }
        else if (!strcmp(a, "-g") || !strcmp(a, "--gaps"))    gapspec = NEXT();
        else if (!strcmp(a, "-l") || !strcmp(a, "--len"))     MINLEN = atoi(NEXT());
        else if (!strcmp(a, "-t") || !strcmp(a, "--threads")) NTHREADS = atoi(NEXT());
        else if (!strcmp(a, "-c") || !strcmp(a, "--consecutive")) REQ_CONSEC = 1;
        else if (!strcmp(a, "-q") || !strcmp(a, "--quiet"))   QUIET = 1;
        else if (!strcmp(a, "--brute")) BRUTE = 1;
        else if (!strcmp(a, "-b") || !strcmp(a, "--bound")) {
            uint64_t v; if (parse_qty(NEXT(), &v) || v < 1000 || v > (1u<<30)) {
                fprintf(stderr, "bad --bound\n"); return 1; }
            user_bound = (uint32_t)v;
        }
        else if (!strcmp(a, "-m") || !strcmp(a, "--mem")) {
            uint64_t v; if (parse_qty(NEXT(), &v) || v < (1u<<16)) {
                fprintf(stderr, "bad --mem\n"); return 1; }
            BM_TARGET = (size_t)v;
        }
        else if (a[0] == '-' && a[1]) { usage(argv[0]); return 1; }
        else {
            if (nargs == 2) { usage(argv[0]); return 1; }
            if (parse_qty(a, &args[nargs])) {
                fprintf(stderr, "cannot parse number: %s\n", a); return 1;
            }
            nargs++;
        }
    }
    if (nargs != 2) { usage(argv[0]); return 1; }
    RSTART = args[0]; REND = args[1];
    if (MINLEN < 2 || MINLEN > MAXLEN) {
        fprintf(stderr, "--len must be 2..%d\n", MAXLEN); return 1;
    }
    if (NTHREADS < 1) NTHREADS = 1;
    if (RSTART > REND) { fprintf(stderr, "start > end\n"); return 1; }

    gen_offsets(gapspec);
    int bad_q = 0;
    ADMISSIBLE = max_admissible(&bad_q);

    if (!QUIET) {
        fprintf(stderr, "# gaps       %s ->", GAPSPEC);
        for (int k = 0; k + 1 < NOFF && k < 12; k++)
            fprintf(stderr, " %llu", (unsigned long long)GAPS[k]);
        fprintf(stderr, "%s\n", NOFF > 13 ? " ..." : "");
        if (bad_q)
            fprintf(stderr, "# ADMISSIBLE only to length %d: at length %d the offsets\n"
                            "#            cover every residue mod %d, so one member is\n"
                            "#            always divisible by %d.\n",
                    ADMISSIBLE, ADMISSIBLE + 1, bad_q, bad_q);
        else
            fprintf(stderr, "# admissible at every length up to %d\n", NOFF);
    }

    if (MINLEN > ADMISSIBLE) {
        fprintf(stderr,
"# no chain of length %d can exist for p > %d, since one member is always\n"
"# divisible by %d. The only possible chains are those in which a member\n"
"# equals %d itself; checking those directly:\n", MINLEN, bad_q, bad_q, bad_q);
        uint64_t hi = REND < 64 ? REND : 64;
        brute_range(RSTART, hi);
        if (!stat_chains) fprintf(stderr, "# none.\n");
        return 0;
    }

    if (MINLEN > NOFF) {
        fprintf(stderr, "--len exceeds the %d offsets this pattern provides\n", NOFF);
        return 1;
    }
    uint64_t maxoff = OFF[NOFF - 1];
    if (REND > UINT64_MAX - maxoff - 2) {
        fprintf(stderr, "end too close to 2^64\n"); return 1;
    }

    /* sieve bound: small L leaves many survivors, so sieve harder there */
    if (user_bound) SIEVE_B = user_bound;
    else SIEVE_B = MINLEN <= 4 ? (1u << 21)
                 : MINLEN <= 6 ? (1u << 20)
                 : MINLEN <= 8 ? (1u << 18)
                               : (1u << 17);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    if (BRUTE) {
        if (!QUIET) fprintf(stderr, "# brute-force mode\n");
        brute_range(RSTART, REND);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (!QUIET)
            fprintf(stderr, "# %llu chains, %.2fs\n",
                    (unsigned long long)stat_chains,
                    (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec));
        return 0;
    }

    uint32_t npr, *pr = make_primes(SIEVE_B, &npr);
    build_wheel(pr, npr, 8192, 2000000000ULL);
    build_prime_tables(pr, npr);

    /* block geometry */
    uint64_t lo = RSTART > (uint64_t)SIEVE_B ? RSTART : (uint64_t)SIEVE_B + 1;
    SIEVE_LO = lo; SIEVE_HI = REND;

    size_t bits_per_block = BM_TARGET * 8;
    uint64_t jl = bits_per_block / NRES;
    if (jl < 1024) jl = 1024;
    if (jl > (1u << 20)) jl = 1u << 20;
    /* enough blocks to keep every thread busy */
    if (SIEVE_HI >= SIEVE_LO) {
        uint64_t want = (SIEVE_HI - SIEVE_LO) / WMOD / ((uint64_t)NTHREADS * 4) + 1;
        if (jl > want) jl = want < 1024 ? 1024 : want;
    }
    jl = (jl + 63) & ~63ULL;
    JLEN = (uint32_t)jl;
    SPAN = WMOD * (uint64_t)JLEN;

    if (!QUIET) {
        char a[32], b[32], c[32];
        fmt_u64(RSTART, a); fmt_u64(REND, b); fmt_u64(WMOD, c);
        double S = singular_series(pr, npr, MINLEN);
        double E = expected_count((double)RSTART, (double)REND, MINLEN, S);
        fprintf(stderr, "# range      %s .. %s\n", a, b);
        fprintf(stderr, "# pattern    len %d, offsets 0", MINLEN);
        for (int k = 1; k < MINLEN; k++)
            fprintf(stderr, ",%llu", (unsigned long long)OFF[k]);
        fprintf(stderr, "%s\n", REQ_CONSEC ? "  (consecutive primes required)" : "");
        fprintf(stderr, "# wheel      W = %s (primes <= %u), %u residues, density 1/%.0f\n",
                c, max_wheel_prime, NRES, (double)WMOD / NRES);
        fprintf(stderr, "# sieve      primes <= %u (%u of them), %u bits/class, "
                        "%.1f MB/thread, span %.3g/block\n",
                SIEVE_B, NSP, JLEN,
                (double)NRES * JLEN / 8 / 1048576.0, (double)SPAN);
        fprintf(stderr, "# estimate   singular series %.4g, expected chains %s%.4g\n",
                S, REQ_CONSEC ? "(ignoring consecutivity) " : "", E);
        fprintf(stderr, "# threads    %d\n", NTHREADS);
    }

    /* small numbers: any chain member could equal a sieving prime, so brute */
    if (RSTART <= (uint64_t)SIEVE_B)
        brute_range(RSTART, REND < (uint64_t)SIEVE_B ? REND : (uint64_t)SIEVE_B);

    if (SIEVE_LO <= SIEVE_HI) {
        BLK_FIRST = SIEVE_LO / SPAN;
        BLK_LAST  = SIEVE_HI / SPAN;
        next_block = BLK_FIRST;
        next_emit  = BLK_FIRST;

        pthread_t *th = malloc((size_t)NTHREADS * sizeof *th);
        worker_t  *wk = calloc((size_t)NTHREADS, sizeof *wk);
        size_t words = JLEN / 64;
        for (int i = 0; i < NTHREADS; i++) {
            wk[i].words = words;
            wk[i].bm = malloc((size_t)NRES * words * sizeof(uint64_t));
            if (!wk[i].bm) { perror("malloc"); return 1; }
            pthread_create(&th[i], NULL, worker, &wk[i]);
        }

        uint64_t total = BLK_LAST - BLK_FIRST + 1;
        while (!QUIET) {
            struct timespec ts = {1, 0};
            nanosleep(&ts, NULL);
            uint64_t done = atomic_load(&stat_blocks);
            if (done >= total) break;
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double el = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
            fprintf(stderr, "\r# %llu/%llu blocks  %.1f%%  %.0fs elapsed, ~%.0fs left   ",
                    (unsigned long long)done, (unsigned long long)total,
                    100.0 * done / total, el,
                    done ? el * (total - done) / done : 0.0);
        }
        for (int i = 0; i < NTHREADS; i++) pthread_join(th[i], NULL);
        if (!QUIET) fprintf(stderr, "\r%60s\r", "");
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (!QUIET) {
        double el = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
        fprintf(stderr, "# %llu chains found, %llu candidates confirmed, %.2fs\n",
                (unsigned long long)stat_chains,
                (unsigned long long)stat_survivors, el);
    }
    return 0;
}
