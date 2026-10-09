/*
 * cpap.c -- search for CPAP-k: k consecutive primes in arithmetic progression.
 *
 *   p, p+d, p+2d, ..., p+(k-1)d   all prime, and NOTHING prime in between.
 *
 * Target: beat the known upper bound for A006560(7),
 *   a(7) <= 71137654873189893604531   (P. Zimmermann)
 * This finds *a* CPAP-7 below that bound. It does NOT prove minimality --
 * that needs exhaustive search and is far out of reach (see README).
 *
 * Two phases:
 *
 * 1. CLASS GENERATION (`cpap classes`)
 *    The hard part of a CPAP is not the k primes, it is the ~1250 interior
 *    integers that must ALL be composite. So we don't search p uniformly:
 *    we pick a residue class p = r (mod M), M = product of primes <= Q, that
 *    forces as many interior positions as possible to be divisible by a small
 *    prime. Beam search over the primes of M, maximising kills.
 *
 *    For k=7,d=210,Q=37: a random admissible class leaves ~200 interior
 *    positions "alive" (needing a primality test and needing to be composite);
 *    the best classes leave ~168. Each position removed multiplies the hit
 *    rate by 1/(1-u) ~ 1.15, so this is worth many orders of magnitude.
 *
 * 2. RANGE SEARCH (`cpap search`)
 *    Within a class, p = r + M*t. Segmented sieve in t-space: for each prime
 *    q and each member i, the t killed form one arithmetic progression,
 *    t = -(r + i*d) * M^-1 (mod q). Survivors get deterministic 128-bit
 *    Miller-Rabin on the members (cheapest rejection first), and only then
 *    the interior positions are checked.
 *
 * Build: gcc -O3 -march=native -o cpap cpap.c -lpthread -lm
 */

#define _GNU_SOURCE
#include "u128.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#if defined(__APPLE__) || defined(__FreeBSD__)
#include <sys/sysctl.h>
#endif

#define MAXK    16
#define MAXSPAN 4096
#define WORDS   (MAXSPAN / 64)

/* ---------------- configuration ---------------- */
static int      K = 7;
static uint64_t Dd = 210;
static uint32_t QMOD = 37;          /* primes <= QMOD form the class modulus */
static uint32_t SIEVE_B = 10000000;
static int      NTHREADS = 1;
static int      QUIET = 0;
static u128     RLO, RHI;
static uint32_t SPAN;               /* (K-1)*d */

static u128     MOD;                /* primorial(QMOD) */
static uint32_t MODP[32], NMODP;

/* ---------------- misc ---------------- */

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

static int parse_u128(const char *s, u128 *out)
{
    char buf[64]; size_t n = 0;
    for (const char *p = s; *p && n < sizeof buf - 1; p++)
        if (*p != ',' && *p != '_') buf[n++] = *p;
    buf[n] = 0;
    if (!n) return -1;
    u128 mult = 1; size_t len = strlen(buf);
    switch (toupper((unsigned char)buf[len-1])) {
        case 'K': mult = 1000ULL; break;
        case 'M': mult = 1000000ULL; break;
        case 'B': case 'G': mult = 1000000000ULL; break;
        case 'T': mult = 1000000000000ULL; break;
        case 'P': mult = 1000000000000000ULL; break;
        case 'E': mult = 1000000000000000000ULL; break;
        default: mult = 0;
    }
    if (mult) buf[--len] = 0; else mult = 1;
    if (!len) return -1;
    if (strchr(buf, '.') || strchr(buf, 'e') || strchr(buf, 'E')) {
        long double v = strtold(buf, NULL) * (long double)mult;
        if (v < 0 || v > 1.7e38L) return -1;
        /* build via two halves to keep precision sane */
        long double hi = floorl(v / 1e18L);
        *out = (u128)(unsigned long long)hi * 1000000000000000000ULL
             + (u128)(unsigned long long)(v - hi * 1e18L);
        return 0;
    }
    u128 v = 0;
    for (size_t i = 0; i < len; i++) {
        if (!isdigit((unsigned char)buf[i])) return -1;
        v = v * 10 + (buf[i] - '0');
    }
    *out = v * mult;
    return 0;
}

static uint32_t *primes_upto(uint32_t lim, uint32_t *cnt)
{
    uint8_t *c = calloc((size_t)lim + 1, 1);
    for (uint64_t i = 2; i * i <= lim; i++)
        if (!c[i]) for (uint64_t j = i * i; j <= lim; j += i) c[j] = 1;
    uint32_t n = 0;
    for (uint32_t i = 2; i <= lim; i++) if (!c[i]) n++;
    uint32_t *p = malloc((size_t)n * sizeof *p); n = 0;
    for (uint32_t i = 2; i <= lim; i++) if (!c[i]) p[n++] = i;
    free(c);
    *cnt = n;
    return p;
}

/* ---------------- bitset over interior positions ---------------- */

static inline int bs_get(const uint64_t *b, uint32_t i) { return (b[i>>6] >> (i&63)) & 1; }
static inline void bs_clr(uint64_t *b, uint32_t i)      { b[i>>6] &= ~(1ULL << (i&63)); }
static inline uint32_t bs_count(const uint64_t *b, uint32_t nw)
{
    uint32_t c = 0;
    for (uint32_t i = 0; i < nw; i++) c += __builtin_popcountll(b[i]);
    return c;
}

/* ---------------- class generation (beam search) ---------------- */

/* Bounded beam expansion.
 *
 * The naive version builds every child (ncur * allowed-residues) and then
 * sorts -- at beam 1e6 and 30 residues that is 3e7 states, several GB. Instead
 * we keep a max-heap of size `beam` keyed on surviving-interior count, so a
 * child is only stored if it beats the current worst. Memory is O(beam), and
 * states live in a flat buffer whose stride follows the actual span rather
 * than MAXSPAN.
 */
static uint32_t ST_NW, ST_STRIDE;

#define ST_ALIVE(buf,i) ((uint64_t *)((char *)(buf) + (size_t)(i) * ST_STRIDE))
#define ST_R(buf,i)     (*(u128 *)((char *)(buf) + (size_t)(i) * ST_STRIDE + ST_NW * 8))
#define ST_S(buf,i)     (*(uint32_t *)((char *)(buf) + (size_t)(i) * ST_STRIDE + ST_NW * 8 + 16))

static void heap_sift_down(void *buf, uint32_t n, uint32_t i, void *tmp)
{
    for (;;) {
        uint32_t l = 2*i+1, r = l+1, big = i;
        if (l < n && ST_S(buf,l) > ST_S(buf,big)) big = l;
        if (r < n && ST_S(buf,r) > ST_S(buf,big)) big = r;
        if (big == i) return;
        memcpy(tmp, (char*)buf + (size_t)i*ST_STRIDE, ST_STRIDE);
        memcpy((char*)buf + (size_t)i*ST_STRIDE, (char*)buf + (size_t)big*ST_STRIDE, ST_STRIDE);
        memcpy((char*)buf + (size_t)big*ST_STRIDE, tmp, ST_STRIDE);
        i = big;
    }
}

static int st_cmp(const void *a, const void *b)
{
    uint32_t x = *(const uint32_t *)((const char *)a + ST_NW*8 + 16);
    uint32_t y = *(const uint32_t *)((const char *)b + ST_NW*8 + 16);
    return (int)x - (int)y;
}

static void *gen_classes(uint32_t beam, uint32_t *out_n)
{
    ST_NW = (SPAN + 63) / 64;
    ST_NW = (ST_NW + 1) & ~1u;               /* even -> r stays 16-byte aligned */
    ST_STRIDE = (ST_NW * 8 + 16 + 8 + 15) & ~15u;   /* u128 needs alignment */

    void *cur = calloc(1, ST_STRIDE);
    uint32_t ncur = 1;
    for (uint32_t x = 1; x < SPAN; x++)
        if (x % Dd) ST_ALIVE(cur,0)[x >> 6] |= 1ULL << (x & 63);
    ST_S(cur,0) = bs_count(ST_ALIVE(cur,0), ST_NW);
    ST_R(cur,0) = 0;

    void *heap = malloc((size_t)beam * ST_STRIDE);
    void *tmp  = malloc(ST_STRIDE);
    void *kid  = malloc(ST_STRIDE);
    if (!heap || !tmp || !kid) { perror("malloc"); exit(1); }

    u128 mdone = 1;
    for (uint32_t pi = 0; pi < NMODP; pi++) {
        uint32_t q = MODP[pi];
        uint8_t bad[512] = {0};
        for (int i = 0; i < K; i++)
            bad[(uint32_t)((q - (uint64_t)Dd * i % q) % q)] = 1;

        uint32_t nh = 0;
        for (uint32_t si = 0; si < ncur; si++) {
            u128 r = ST_R(cur,si);
            for (uint32_t a = 0; a < q; a++) {
                if (bad[a]) continue;
                u128 rr = r;
                while (rr % q != a) rr += mdone;     /* CRT lift */
                memcpy(ST_ALIVE(kid,0), ST_ALIVE(cur,si), ST_NW * 8);
                for (uint32_t x = 1; x < SPAN; x++)
                    if (bs_get(ST_ALIVE(kid,0), x) && (uint32_t)((rr + x) % q) == 0)
                        bs_clr(ST_ALIVE(kid,0), x);
                uint32_t sc = bs_count(ST_ALIVE(kid,0), ST_NW);
                ST_R(kid,0) = rr; ST_S(kid,0) = sc;

                if (nh < beam) {                       /* push */
                    memcpy((char*)heap + (size_t)nh*ST_STRIDE, kid, ST_STRIDE);
                    uint32_t i = nh++;
                    while (i) {
                        uint32_t p = (i-1)/2;
                        if (ST_S(heap,p) >= ST_S(heap,i)) break;
                        memcpy(tmp,(char*)heap+(size_t)i*ST_STRIDE,ST_STRIDE);
                        memcpy((char*)heap+(size_t)i*ST_STRIDE,(char*)heap+(size_t)p*ST_STRIDE,ST_STRIDE);
                        memcpy((char*)heap+(size_t)p*ST_STRIDE,tmp,ST_STRIDE);
                        i = p;
                    }
                } else if (sc < ST_S(heap,0)) {         /* replace worst */
                    memcpy(heap, kid, ST_STRIDE);
                    heap_sift_down(heap, nh, 0, tmp);
                }
            }
        }
        mdone *= q;
        free(cur);
        cur = malloc((size_t)nh * ST_STRIDE);
        memcpy(cur, heap, (size_t)nh * ST_STRIDE);
        ncur = nh;
    }
    qsort(cur, ncur, ST_STRIDE, st_cmp);
    free(heap); free(tmp); free(kid);
    *out_n = ncur;
    return cur;
}

/* ---------------- search ---------------- */

static uint32_t *SP; static uint32_t NSP;
static uint32_t *MINV;              /* M^-1 mod q, precomputed once */

static atomic_ullong stat_cand, stat_deep, stat_hit;
static uint32_t SEGBITS = 1u << 20;

static void *CLASSES; static uint32_t NCLASS;
static _Atomic uint32_t class_next;

static uint32_t modinv32(uint64_t a, uint32_t q)
{
    /* q prime -> a^(q-2) */
    uint64_t r = 1, b = a % q, e = q - 2;
    while (e) { if (e & 1) r = r * b % q; b = b * b % q; e >>= 1; }
    return (uint32_t)r;
}

typedef struct { uint64_t *bm; uint32_t *off; } wctx;

static pthread_mutex_t hitmx = PTHREAD_MUTEX_INITIALIZER;

static void report_hit(u128 p, u128 clr)
{
    pthread_mutex_lock(&hitmx);
    char b[48];
    printf("\n*** CPAP-%d FOUND ***\n", K);
    u128_str(p, b);
    printf("p = %s   (d = %llu)\n", b, (unsigned long long)Dd);
    for (int i = 0; i < K; i++) {
        u128_str(p + (u128)Dd * i, b);
        printf("  %s\n", b);
    }
    u128_str(clr, b); printf("class r = %s\n", b);
    fflush(stdout);
    pthread_mutex_unlock(&hitmx);
}

static void *worker(void *arg)
{
    wctx *w = arg;
    uint32_t nw = (SPAN + 63) / 64;
    uint32_t words = SEGBITS / 64;

    for (;;) {
        uint32_t ci = atomic_fetch_add(&class_next, 1);
        if (ci >= NCLASS) break;
        uint64_t *cl_alive = ST_ALIVE(CLASSES, ci);
        u128 cl_r = ST_R(CLASSES, ci);

        /* interior positions still alive for this class */
        uint32_t alivep[MAXSPAN]; uint32_t nalive = 0;
        for (uint32_t x = 1; x < SPAN; x++)
            if (bs_get(cl_alive, x)) alivep[nalive++] = x;

        /* t range so that r + M*t lies in [RLO, RHI] */
        u128 lo = cl_r >= RLO ? 0 : (RLO - cl_r + MOD - 1) / MOD;
        if (cl_r > RHI) continue;
        u128 hi = (RHI - cl_r) / MOD;

        /* per (prime, member) start offsets in t-space */
        for (uint32_t pi = 0; pi < NSP; pi++) {
            uint32_t q = SP[pi];
            uint64_t rq = (uint64_t)(cl_r % q);
            for (int i = 0; i < K; i++) {
                uint64_t v = (rq + (uint64_t)Dd % q * i) % q;
                uint64_t nv = v ? q - v : 0;
                w->off[(size_t)pi * K + i] = (uint32_t)(nv * MINV[pi] % q);
            }
        }

        for (u128 base = lo; base <= hi; base += SEGBITS) {
            uint32_t len = SEGBITS;
            if (hi - base + 1 < SEGBITS) len = (uint32_t)(hi - base + 1);
            memset(w->bm, 0, (size_t)words * sizeof(uint64_t));

            for (uint32_t pi = 0; pi < NSP; pi++) {
                uint32_t q = SP[pi];
                uint32_t bmod = (uint32_t)(base % q);
                for (int i = 0; i < K; i++) {
                    uint32_t st = w->off[(size_t)pi * K + i];
                    uint32_t j = st >= bmod ? st - bmod : st + q - bmod;
                    for (; j < len; j += q) w->bm[j >> 6] |= 1ULL << (j & 63);
                }
            }

            uint64_t cand = 0, deep = 0;
            for (uint32_t wi = 0; wi * 64 < len; wi++) {
                uint64_t m = ~w->bm[wi];
                if (wi * 64 + 64 > len) m &= (len - wi * 64) >= 64 ? ~0ULL
                                          : ((1ULL << (len - wi * 64)) - 1);
                while (m) {
                    uint32_t j = wi * 64 + __builtin_ctzll(m);
                    m &= m - 1;
                    cand++;
                    u128 p = cl_r + MOD * (base + j);
                    /* members: reject as early as possible */
                    int ok = 1;
                    for (int i = 0; i < K; i++)
                        if (!is_prime_u128(p + (u128)Dd * i)) { ok = 0; break; }
                    if (!ok) continue;
                    deep++;
                    for (uint32_t a = 0; a < nalive; a++)
                        if (is_prime_u128(p + alivep[a])) { ok = 0; break; }
                    if (!ok) continue;
                    report_hit(p, cl_r);
                    atomic_fetch_add(&stat_hit, 1);
                }
            }
            atomic_fetch_add(&stat_cand, cand);
            atomic_fetch_add(&stat_deep, deep);
            (void)nw;
        }
    }
    return NULL;
}

/* ---------------- main ---------------- */

static void usage(const char *p)
{
    fprintf(stderr,
"usage: %s [options]\n"
"\n"
"  Searches for CPAP-k: k consecutive primes in arithmetic progression.\n"
"  Default targets A006560(7): k=7, d=210, below Zimmermann's bound.\n"
"\n"
"  -k N          number of primes in the progression (default 7)\n"
"  -d N          common difference (default 210; must be divisible by\n"
"                every prime <= k, else no CPAP-k exists)\n"
"  -Q N          class modulus uses primes <= N (default 37)\n"
"  -w N          beam width / number of classes to search (default 4096)\n"
"  -b N          sieve bound (default 10M)\n"
"  -s N          segment size in t, power of 2 (default 1M)\n"
"  -t N          threads (default: online CPUs)\n"
"      --from X  lower end of the p range (default 0)\n"
"      --to X    upper end (default 71137654873189893604531, the record)\n"
"      --classes just print the top classes and exit\n"
"      --class-start N / --class-count N   split the class list across\n"
"                machines: node j of J runs --class-start j*C --class-count C\n"
"  -q            quiet\n"
"\n"
"  X accepts K M B/G T P E suffixes: 1.5T, 250B, 7.1e22\n", p);
}

int main(int argc, char **argv)
{
    uint32_t beam = 4096;
    uint32_t CLS_SKIP = 0, CLS_TAKE = 0;
    int only_classes = 0;
    NTHREADS = cpu_count();
    RLO = 0;
    parse_u128("71137654873189893604531", &RHI);

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        #define NEXT() (i+1 < argc ? argv[++i] : (usage(argv[0]), exit(1), (char*)0))
        if (!strcmp(a,"-h")||!strcmp(a,"--help")) { usage(argv[0]); return 0; }
        else if (!strcmp(a,"-k")) K = atoi(NEXT());
        else if (!strcmp(a,"-d")) Dd = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(a,"-Q")) QMOD = (uint32_t)atoi(NEXT());
        else if (!strcmp(a,"-w")) beam = (uint32_t)strtoul(NEXT(), NULL, 10);
        else if (!strcmp(a,"-b")) { u128 v; parse_u128(NEXT(), &v); SIEVE_B = (uint32_t)v; }
        else if (!strcmp(a,"-s")) { u128 v; parse_u128(NEXT(), &v); SEGBITS = (uint32_t)v; }
        else if (!strcmp(a,"-t")) NTHREADS = atoi(NEXT());
        else if (!strcmp(a,"--from")) { if (parse_u128(NEXT(), &RLO)) return 1; }
        else if (!strcmp(a,"--to"))   { if (parse_u128(NEXT(), &RHI)) return 1; }
        else if (!strcmp(a,"--classes")) only_classes = 1;
        else if (!strcmp(a,"--class-start")) CLS_SKIP = (uint32_t)strtoul(NEXT(), NULL, 10);
        else if (!strcmp(a,"--class-count")) CLS_TAKE = (uint32_t)strtoul(NEXT(), NULL, 10);
        else if (!strcmp(a,"-q")) QUIET = 1;
        else { usage(argv[0]); return 1; }
    }
    if (K < 3 || K > MAXK) { fprintf(stderr, "-k out of range\n"); return 1; }
    SPAN = (uint32_t)(Dd * (K - 1));
    if (SPAN >= MAXSPAN) { fprintf(stderr, "span too large (raise MAXSPAN)\n"); return 1; }
    SEGBITS = (SEGBITS + 63) & ~63u;
    if (NTHREADS < 1) NTHREADS = 1;

    /* d must be divisible by every prime <= K, else CPAP-k is impossible */
    for (int q = 2; q <= K; q++) {
        int isp = 1;
        for (int j = 2; j * j <= q; j++) if (q % j == 0) { isp = 0; break; }
        if (isp && Dd % q) {
            fprintf(stderr,
              "d=%llu is not divisible by %d, so one of the %d terms is always\n"
              "divisible by %d -- no CPAP-%d with this difference exists.\n",
              (unsigned long long)Dd, q, K, q, K);
            return 1;
        }
    }

    /* class modulus */
    uint32_t nq; uint32_t *qs = primes_upto(QMOD, &nq);
    MOD = 1; NMODP = 0;
    for (uint32_t i = 0; i < nq; i++) { MODP[NMODP++] = qs[i]; MOD *= qs[i]; }
    free(qs);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    uint32_t ncl;
    CLASSES = gen_classes(beam, &ncl);
    NCLASS = ncl;

    if (!QUIET || only_classes) {
        char b[48]; u128_str(MOD, b);
        uint32_t interior = 0;
        for (uint32_t x = 1; x < SPAN; x++) if (x % Dd) interior++;
        fprintf(stderr, "# target     CPAP-%d, d=%llu, span %u, %u interior positions\n",
                K, (unsigned long long)Dd, SPAN, interior);
        fprintf(stderr, "# class mod  M = %s (primes <= %u)\n", b, QMOD);
        fprintf(stderr, "# classes    %u kept; best leaves %u alive, worst kept %u\n",
                NCLASS, ST_S(CLASSES,0), ST_S(CLASSES,NCLASS-1));
        u128 per = MOD ? (RHI - RLO) / MOD : 0;
        char b2[48]; u128_str(per, b2);
        fprintf(stderr, "# per class  ~%s candidates in range\n", b2);
    }
    if (CLS_SKIP) {
        if (CLS_SKIP >= NCLASS) { fprintf(stderr, "--class-start beyond class count\n"); return 1; }
        CLASSES = (char *)CLASSES + (size_t)CLS_SKIP * ST_STRIDE; NCLASS -= CLS_SKIP;
    }
    if (CLS_TAKE && CLS_TAKE < NCLASS) NCLASS = CLS_TAKE;

    if (only_classes) {
        for (uint32_t i = 0; i < NCLASS && i < 40; i++) {
            char b[48]; u128_str(ST_R(CLASSES,i), b);
            printf("alive=%3u  r=%s\n", ST_S(CLASSES,i), b);
        }
        return 0;
    }

    /* sieving primes: those not dividing M */
    uint32_t nall; uint32_t *all = primes_upto(SIEVE_B, &nall);
    NSP = 0;
    for (uint32_t i = 0; i < nall; i++) if (all[i] > QMOD) NSP++;
    SP = malloc((size_t)NSP * sizeof *SP);
    MINV = malloc((size_t)NSP * sizeof *MINV);
    uint32_t n = 0;
    for (uint32_t i = 0; i < nall; i++) {
        if (all[i] <= QMOD) continue;
        SP[n] = all[i];
        MINV[n] = modinv32((uint64_t)(MOD % all[i]), all[i]);
        n++;
    }
    free(all);

    if (!QUIET) {
        fprintf(stderr, "# sieve      primes <= %u (%u), segment %u, %d threads\n",
                SIEVE_B, NSP, SEGBITS, NTHREADS);
        fprintf(stderr, "# searching...\n");
    }

    pthread_t *th = malloc((size_t)NTHREADS * sizeof *th);
    wctx *wk = calloc((size_t)NTHREADS, sizeof *wk);
    for (int i = 0; i < NTHREADS; i++) {
        wk[i].bm  = malloc((size_t)(SEGBITS / 64) * sizeof(uint64_t));
        wk[i].off = malloc((size_t)NSP * K * sizeof(uint32_t));
        pthread_create(&th[i], NULL, worker, &wk[i]);
    }
    if (!QUIET) {
        for (;;) {
            struct timespec ts = {5, 0};
            nanosleep(&ts, NULL);
            uint32_t done = atomic_load(&class_next);
            if (done >= NCLASS) break;
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double el = (t1.tv_sec-t0.tv_sec)+1e-9*(t1.tv_nsec-t0.tv_nsec);
            fprintf(stderr, "\r# class %u/%u  %llu cand  %llu AP-%d  %.0fs, ~%.1fh left    ",
                    done, NCLASS, (unsigned long long)stat_cand,
                    (unsigned long long)stat_deep, K, el,
                    done ? el*(NCLASS-done)/done/3600.0 : 0.0);
        }
        fprintf(stderr, "\r%78s\r", "");
    }
    for (int i = 0; i < NTHREADS; i++) pthread_join(th[i], NULL);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double el = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
    if (!QUIET)
        fprintf(stderr,
            "# done: %llu candidates sieved, %llu full AP-%d, %llu CPAP found, %.1fs\n",
            (unsigned long long)stat_cand, (unsigned long long)stat_deep, K,
            (unsigned long long)stat_hit, el);
    return 0;
}
