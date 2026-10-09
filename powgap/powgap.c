/*
 * powgap.c -- search for chains of primes whose gaps double: 2,4,8,16,32,...
 *
 *   p, p+2, p+6, p+14, p+30, p+62, p+126, p+254, ...     offset(k) = 2^(k+1)-2
 *
 * Method: constellation (k-tuple) sieve in a compressed wheel space, exactly
 * like the increasing-gap chain search, but the admissible-residue structure
 * here is far tighter, so the wheel does much more of the work.
 *
 *   - mod 2  : p must be odd                      -> 1 of 2 residues
 *   - mod 3  : offsets are {0,2} mod 3            -> 1 of 3 residues
 *   - mod 5  : offsets cycle {0,2,1,4}            -> 1 of 5 residues (len>=4)
 *   - mod q  : offsets take ord_q(2) values       -> q-ord_q(2) residues
 *
 * Since ord_q(2) <= q-1 the pattern is admissible for every length, but the
 * survivor density drops fast, so we sieve only the wheel survivors with
 * small primes and BPSW/Miller-Rabin only what is left.
 *
 * With -c (consecutive primes required) this reproduces OEIS A090807 (and
 * A079014, the same sequence indexed differently): a(n) = first prime in
 * the earliest chain of n consecutive primes with gaps 2^1, 2^2, ...,
 * 2^(n-1).  All published terms a(2)..a(8) -- 3, 5, 1997, 2237, 6824897,
 * 1356705137, 3637803390827 -- are confirmed by exhaustive search with
 * this program:  powgap -l N -c 1 <bound>   (https://oeis.org/A090807)
 *
 * The next term, a(9) = 14014732040120297, was found with this program
 * (Aug 2026, now published in A090807/A079014): the sweep  powgap -l 9 -c -b 4096 3.6e12 1.5e16  shows it is
 * the only chain of 9 consecutive primes in [a(8), 1.40148e16], and no
 * 9-chain can start below a(8) since its 8-prefix would contradict a(8).
 *
 * Build:  gcc -O3 -march=native -o powgap powgap.c -lpthread -lm
 *         clang -O3 -mcpu=native -o powgap powgap.c -lpthread -lm   (Apple silicon)
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
#include <signal.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#if defined(__APPLE__) || defined(__FreeBSD__)
#include <sys/sysctl.h>
#endif

#define MAXLEN 24

/* ------------------------------------------------------------------ */
/* globals (read-only after setup)                                     */
/* ------------------------------------------------------------------ */

static uint64_t OFF[MAXLEN + 1];      /* OFF[k] = 2^(k+1) - 2 */
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
static uint64_t *RES;
static uint64_t  WHEEL_MAXRES = 8192;   /* --wheel: max residues */

/* sieving primes (all > largest wheel prime, <= SIEVE_B) */
static uint32_t  NSP;
static uint32_t *SP;        /* prime            */
static uint32_t *SPINV;     /* W^-1 mod p       */
static uint32_t *DK;        /* NSP * MAXLEN: distinct (-off*W^-1) mod p */
static uint8_t  *DKN;       /* how many are distinct                   */
static uint32_t *RI;        /* NRES * NSP: RES[i]*W^-1 mod p (read-only) */

static uint64_t SPAN;       /* integers covered by one block */
static uint64_t BLK_FIRST, BLK_LAST;
static uint64_t SIEVE_LO, SIEVE_HI;

static atomic_ullong next_block;
static atomic_ullong stat_survivors, stat_chains, stat_blocks;
static atomic_ullong HIST[MAXLEN + 1];   /* chains found, by full length  */
static atomic_ullong CHIST[MAXLEN + 1];  /* by consecutive-prefix length  */

static volatile sig_atomic_t stop_requested;
static const char *CKPT_PATH;            /* --checkpoint file, or NULL    */
static double   ELAPSED_PREV;            /* elapsed time from prior runs  */
static uint64_t RESUME_FROM;             /* first number not yet searched */
static uint64_t BRUTE_POS;               /* brute-force progress cursor   */
static struct timespec T0;
#define CKPT_EVERY 10.0                  /* checkpoint interval, seconds  */

/* ------------------------------------------------------------------ */
/* number parsing: 12, 1.5T, 250B, 3e15, 1_000_000                      */
/* ------------------------------------------------------------------ */

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

static void fmt_dur(double s, char *out)        /* 5025 -> "1h23m45s" */
{
    if (s < 0) s = 0;
    long v = (long)(s + 0.5);
    if (v >= 86400)     snprintf(out, 24, "%ldd%02ldh%02ldm", v / 86400, v % 86400 / 3600, v % 3600 / 60);
    else if (v >= 3600) snprintf(out, 24, "%ldh%02ldm%02lds", v / 3600, v % 3600 / 60, v % 60);
    else if (v >= 60)   snprintf(out, 24, "%ldm%02lds", v / 60, v % 60);
    else                snprintf(out, 24, "%lds", v);
}

static double elapsed_now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (t.tv_sec - T0.tv_sec) + 1e-9 * (t.tv_nsec - T0.tv_nsec);
}

static void on_signal(int sig)
{
    (void)sig;
    stop_requested = 1;
    static const char msg[] =
        "\n# interrupt -- finishing current blocks (ctrl-c again to abort)\n";
    ssize_t r = write(2, msg, sizeof msg - 1);
    (void)r;
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
    for (int k = 0; k < MAXLEN; k++) {
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

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static void build_wheel(const uint32_t *pr, uint32_t npr,
                        uint64_t maxres, uint64_t maxW)
{
    uint64_t  M   = 1;
    uint32_t  cnt = 1;
    uint64_t *cur = malloc(sizeof *cur);
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

        uint64_t *nw = malloc((size_t)cnt * good * sizeof *nw);
        uint32_t  nn = 0;
        for (uint32_t t = 0; t < q; t++)
            for (uint32_t i = 0; i < cnt; i++) {
                uint64_t v = cur[i] + M * t;
                if (!bad[v % q]) nw[nn++] = v;
            }
        free(bad); free(cur);
        cur = nw; cnt = nn; M *= q;
        max_wheel_prime = q;
    }

    qsort(cur, cnt, sizeof *cur, cmp_u64);

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

    /* RI[i*NSP+pi] = RES[i] * W^-1 mod p, so that the per-block sieve start
       for (residue i, prime p) is a single subtraction instead of two
       64-bit divisions: b = (-(base+RES[i]) * W^-1) = Bp - RI  (mod p) */
    RI = malloc((size_t)NRES * NSP * sizeof *RI);
    if (!RI) { perror("malloc"); exit(1); }
    for (uint32_t i = 0; i < NRES; i++)
        for (uint32_t p = 0; p < NSP; p++)
            RI[(size_t)i * NSP + p] = (uint32_t)(RES[i] % SP[p] * SPINV[p] % SP[p]);
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
    if (a < 3) a = 3;
    if (b <= a) return 0;
    const int N = 20000;
    double h = (b - a) / N, sum = 0;
    for (int i = 0; i <= N; i++) {
        double x = a + h * i;
        double f = 1.0 / pow(log(x), L);
        sum += (i == 0 || i == N) ? f : (i & 1 ? 4 * f : 2 * f);
    }
    return S * sum * h / 3.0;
}

/* ------------------------------------------------------------------ */
/* ordered output                                                       */
/* ------------------------------------------------------------------ */

typedef struct Res {
    uint64_t  blk;
    char     *txt;
    uint16_t *lens;      /* (len, consec-len) pairs, one per chain */
    uint32_t  nlens;
    uint64_t  surv;
    struct Res *next;
} Res;
static Res *pending;
static uint64_t next_emit;
static pthread_mutex_t outmx = PTHREAD_MUTEX_INITIALIZER;

static void submit(uint64_t blk, char *txt,
                   uint16_t *lens, uint32_t nlens, uint64_t surv)
{
    pthread_mutex_lock(&outmx);
    Res *nd = malloc(sizeof *nd);
    nd->blk = blk; nd->txt = txt; nd->lens = lens; nd->nlens = nlens;
    nd->surv = surv; nd->next = pending; pending = nd;
    for (;;) {
        Res **pp = &pending, *hit = NULL;
        while (*pp) {
            if ((*pp)->blk == next_emit) { hit = *pp; *pp = hit->next; break; }
            pp = &(*pp)->next;
        }
        if (!hit) break;
        if (hit->txt && *hit->txt) { fputs(hit->txt, stdout); fflush(stdout); }
        /* stats count only in-order-emitted blocks, so a checkpoint's
           counters always agree exactly with the output produced so far */
        for (uint32_t i = 0; i < hit->nlens; i++) {
            atomic_fetch_add(&HIST[hit->lens[2*i]], 1);
            atomic_fetch_add(&CHIST[hit->lens[2*i+1]], 1);
        }
        atomic_fetch_add(&stat_chains, hit->nlens);
        atomic_fetch_add(&stat_survivors, hit->surv);
        free(hit->txt); free(hit->lens); free(hit);
        next_emit++;
    }
    pthread_mutex_unlock(&outmx);
}

/* ------------------------------------------------------------------ */
/* checkpointing                                                        */
/* ------------------------------------------------------------------ */

/* first number of the lowest block not yet fully emitted */
static uint64_t sieve_resume(void)
{
    pthread_mutex_lock(&outmx);
    uint64_t ne = next_emit;
    pthread_mutex_unlock(&outmx);
    if (ne > REND / SPAN) return REND + 1;
    uint64_t r = ne * SPAN;
    return r < SIEVE_LO ? SIEVE_LO : r;
}

static void ckpt_write(uint64_t resume)
{
    if (!CKPT_PATH) return;
    char tmp[1024];
    snprintf(tmp, sizeof tmp, "%s.tmp", CKPT_PATH);
    FILE *f = fopen(tmp, "w");
    if (!f) {
        fprintf(stderr, "# cannot write checkpoint %s: %s\n", tmp, strerror(errno));
        return;
    }
    fprintf(f, "powgap 1\nstart %llu\nend %llu\nlen %d\nconsec %d\n"
               "resume %llu\nelapsed %.3f\nchains %llu\nsurvivors %llu\n",
            (unsigned long long)RSTART, (unsigned long long)REND,
            MINLEN, REQ_CONSEC,
            (unsigned long long)resume, ELAPSED_PREV + elapsed_now(),
            (unsigned long long)atomic_load(&stat_chains),
            (unsigned long long)atomic_load(&stat_survivors));
    fprintf(f, "hist");
    for (int k = 0; k <= MAXLEN; k++)
        fprintf(f, " %llu", (unsigned long long)atomic_load(&HIST[k]));
    fprintf(f, "\nchist");
    for (int k = 0; k <= MAXLEN; k++)
        fprintf(f, " %llu", (unsigned long long)atomic_load(&CHIST[k]));
    fprintf(f, "\n");
    if (fclose(f) == 0) rename(tmp, CKPT_PATH);
}

static void ckpt_load(void)
{
    FILE *f = fopen(CKPT_PATH, "r");
    if (!f) return;                                        /* fresh start */
    char tag[16]; int ver = 0, ln = 0, co = 0; double el = 0;
    unsigned long long st = 0, en = 0, res = 0, ch = 0, su = 0, v;
    int ok = fscanf(f, "%15s %d", tag, &ver) == 2
          && !strcmp(tag, "powgap") && ver == 1
          && fscanf(f, " start %llu end %llu len %d consec %d resume %llu "
                       "elapsed %lf chains %llu survivors %llu",
                    &st, &en, &ln, &co, &res, &el, &ch, &su) == 8
          && fscanf(f, "%15s", tag) == 1 && !strcmp(tag, "hist");
    for (int k = 0; ok && k <= MAXLEN; k++)
        if (fscanf(f, "%llu", &v) == 1) HIST[k] = v; else ok = 0;
    ok = ok && fscanf(f, "%15s", tag) == 1 && !strcmp(tag, "chist");
    for (int k = 0; ok && k <= MAXLEN; k++)
        if (fscanf(f, "%llu", &v) == 1) CHIST[k] = v; else ok = 0;
    fclose(f);
    if (!ok) {
        fprintf(stderr, "corrupt checkpoint %s -- delete it to start over\n",
                CKPT_PATH);
        exit(1);
    }
    if (st != RSTART || en != REND || ln != MINLEN || co != REQ_CONSEC) {
        fprintf(stderr, "checkpoint %s is for different parameters "
                        "(start %llu, end %llu, len %d%s) -- "
                        "delete it to start over\n",
                CKPT_PATH, st, en, ln, co ? ", consecutive" : "");
        exit(1);
    }
    RESUME_FROM    = res;
    ELAPSED_PREV   = el;
    stat_chains    = ch;
    stat_survivors = su;
}

/* ------------------------------------------------------------------ */
/* final statistics                                                     */
/* ------------------------------------------------------------------ */

static void print_stats(int completed)
{
    char cb[32], sb2[32], eb[24];
    fmt_u64(atomic_load(&stat_chains), cb);
    fmt_u64(atomic_load(&stat_survivors), sb2);
    fmt_dur(ELAPSED_PREV + elapsed_now(), eb);
    fprintf(stderr, "# %s -- %s chains found, %s candidates confirmed, %s elapsed\n",
            completed ? "finished" : "interrupted", cb, sb2, eb);
    uint64_t tot = 0;
    for (int k = 0; k <= MAXLEN; k++) tot += atomic_load(&HIST[k]);
    if (!tot) {
        fprintf(stderr, "# no chains of length >= %d found\n", MINLEN);
        return;
    }
    fprintf(stderr, "# chains by length:\n");
    for (int k = 0; k <= MAXLEN; k++) {
        uint64_t n = atomic_load(&HIST[k]);
        if (!n) continue;
        char b[32]; fmt_u64(n, b);
        fprintf(stderr, "#   len %2d : %s\n", k, b);
    }
    if (REQ_CONSEC) {
        fprintf(stderr, "# chains by consecutive-prime length:\n");
        for (int k = 0; k <= MAXLEN; k++) {
            uint64_t n = atomic_load(&CHIST[k]);
            if (!n) continue;
            char b[32]; fmt_u64(n, b);
            fprintf(stderr, "#   len %2d : %s\n", k, b);
        }
    }
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

static void *worker(void *arg)
{
    worker_t *w = arg;
    const uint32_t stride = (uint32_t)w->words;
    const int      L      = MINLEN;

    uint64_t *hits = NULL; size_t nhits, caphits = 0;
    uint32_t *Bp = malloc((size_t)NSP * sizeof *Bp);   /* (-base*W^-1) mod q */

    for (;;) {
        if (stop_requested) break;
        uint64_t blk = atomic_fetch_add(&next_block, 1);
        if (blk > BLK_LAST) break;
        uint64_t base = blk * SPAN;

        memset(w->bm, 0, (size_t)NRES * stride * sizeof *w->bm);

        /* ---- constellation sieve ------------------------------------ */
        /* residue-class outer loop: each bitmap row is finished while it
           is hot in cache, which matters once --wheel makes the bitmap
           much larger than L2 */
        for (uint32_t pi = 0; pi < NSP; pi++) {
            uint32_t q = SP[pi];
            uint32_t t = (uint32_t)(base % q * SPINV[pi] % q);
            Bp[pi] = t ? q - t : 0;
        }

        for (uint32_t i = 0; i < NRES; i++) {
            uint64_t *row = w->bm + (size_t)i * stride;
            const uint32_t *ri = RI + (size_t)i * NSP;
            for (uint32_t pi = 0; pi < NSP; pi++) {
                const uint32_t  q   = SP[pi];
                const uint32_t *dk  = DK + (size_t)pi * L;
                const uint8_t   nd  = DKN[pi];
                uint32_t b = Bp[pi] >= ri[pi] ? Bp[pi] - ri[pi] : Bp[pi] + q - ri[pi];
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

        uint16_t *lens = NULL;
        if (nhits) {
            lens = malloc(nhits * 2 * sizeof *lens);
            for (size_t a = 0; a < nhits; a++) {
                lens[2*a]   = (uint16_t)hits[a*3+1];
                lens[2*a+1] = (uint16_t)hits[a*3+2];
            }
        }
        atomic_fetch_add(&stat_blocks, 1);
        submit(blk, s.b, lens, (uint32_t)nhits, surv);
    }
    free(hits); free(Bp);
    return NULL;
}

/* ------------------------------------------------------------------ */

static void brute_range(uint64_t lo, uint64_t hi)
{
    if (lo > hi) { BRUTE_POS = hi + 1; return; }
    if (lo < 2) lo = 2;
    sbuf s = {0};
    double t_start = elapsed_now(), last_tick = t_start, last_ck = t_start;
    int interrupted_here = 0;
    uint64_t c;
    BRUTE_POS = lo;
    for (c = lo; c <= hi; c++) {
        if (!(c & 0x3fff)) {
            if (stop_requested) { interrupted_here = 1; break; }
            double now = elapsed_now();
            if (now - last_tick >= 1.0) {
                last_tick = now;
                if (s.b && s.n) {
                    fputs(s.b, stdout); fflush(stdout);
                    s.n = 0; s.b[0] = 0;
                }
                BRUTE_POS = c;
                if (BRUTE && !QUIET) {
                    double done = (double)(c - lo), tot = (double)(hi - lo) + 1;
                    double el = now - t_start;
                    char cb[32], eb[24], rb[24];
                    fmt_u64(atomic_load(&stat_chains), cb);
                    fmt_dur(ELAPSED_PREV + el, eb);
                    fmt_dur(done > 0 ? el * (tot - done) / done : 0, rb);
                    fprintf(stderr, "\r# %.2f%%  %s chains  %s elapsed, ~%s left   ",
                            100.0 * done / tot, cb, eb, rb);
                }
                if (now - last_ck >= CKPT_EVERY) { last_ck = now; ckpt_write(c); }
            }
        }
        if (c > 2 && !(c & 1)) { if (c == hi) break; else continue; }
        int n = chain_len_at(c);
        if (n >= MINLEN) {
            int m = n;
            if (!REQ_CONSEC || (m = consec_len(c, n)) >= MINLEN) {
                report(&s, c, n, m);
                atomic_fetch_add(&stat_chains, 1);
                atomic_fetch_add(&HIST[n], 1);
                atomic_fetch_add(&CHIST[m], 1);
            }
        }
        if (c == hi) break;
    }
    BRUTE_POS = interrupted_here ? c : hi + 1;
    if (s.b) { if (s.n) { fputs(s.b, stdout); fflush(stdout); } free(s.b); }
}

static void usage(const char *p)
{
    fprintf(stderr,
"usage: %s [options] <start> <end>\n"
"\n"
"  Finds primes p such that p, p+2, p+6, p+14, p+30, ... are all prime\n"
"  (gaps 2,4,8,16,32,... -- offset k is 2^(k+1)-2).\n"
"\n"
"  start/end accept K M B/G T P E suffixes and decimals: 250B, 1.5T, 3e15\n"
"\n"
"  -l, --len N        minimum chain length to report (default 6)\n"
"  -t, --threads N    worker threads (default: online CPUs)\n"
"  -c, --consecutive  require the members to be CONSECUTIVE primes\n"
"                     (nothing prime in between -- see notes in README)\n"
"  -b, --bound N      small-prime sieve bound (default: auto)\n"
"  -m, --mem N        bitmap bytes per thread (default 4M)\n"
"  -w, --wheel N      max wheel residues (default 8192); larger values fold\n"
"                     more primes into the wheel: faster, more memory/thread\n"
"  -C, --checkpoint F save progress to file F every %.0fs; if F exists and\n"
"                     matches the parameters, resume from where it left off\n"
"      --brute        naive verification mode (no sieve) -- for testing\n"
"  -q, --quiet        suppress the header/progress lines\n"
"  -h, --help         this message\n"
"\n"
"  ctrl-c stops cleanly: workers finish their current block, a checkpoint\n"
"  is saved (with -C), and chain-length statistics are printed.\n",
            p, CKPT_EVERY);
}

int main(int argc, char **argv)
{
    uint64_t args[2]; int nargs = 0;
    NTHREADS = cpu_count();
    uint32_t user_bound = 0;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        #define NEXT() (i + 1 < argc ? argv[++i] : (usage(argv[0]), exit(1), (char*)0))
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) { usage(argv[0]); return 0; }
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
        else if (!strcmp(a, "-C") || !strcmp(a, "--checkpoint")) CKPT_PATH = NEXT();
        else if (!strcmp(a, "-w") || !strcmp(a, "--wheel")) {
            uint64_t v; if (parse_qty(NEXT(), &v) || v < 1 || v > (1ULL << 28)) {
                fprintf(stderr, "bad --wheel\n"); return 1; }
            WHEEL_MAXRES = v;
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

    for (int k = 0; k <= MAXLEN; k++) OFF[k] = (1ULL << (k + 1)) - 2;
    uint64_t maxoff = OFF[MAXLEN - 1];
    if (REND > UINT64_MAX - maxoff - 2) {
        fprintf(stderr, "end too close to 2^64\n"); return 1;
    }

    /* sieve bound: small L leaves many survivors, so sieve harder there */
    if (user_bound) SIEVE_B = user_bound;
    else SIEVE_B = MINLEN <= 4 ? (1u << 21)
                 : MINLEN <= 6 ? (1u << 20)
                 : MINLEN <= 8 ? (1u << 18)
                               : (1u << 17);

    clock_gettime(CLOCK_MONOTONIC, &T0);

    if (CKPT_PATH) ckpt_load();
    if (RESUME_FROM > RSTART && !QUIET) {
        char rb[32]; fmt_u64(RESUME_FROM, rb);
        fprintf(stderr, "# resuming from checkpoint %s at %s\n", CKPT_PATH, rb);
    }
    if (RESUME_FROM > REND) {
        fprintf(stderr, "# checkpoint shows this range is already complete\n");
        print_stats(1);
        return 0;
    }

    struct sigaction sa = {0};
    sa.sa_handler = on_signal;
    sa.sa_flags   = SA_RESETHAND;      /* a second ctrl-c kills immediately */
    sigaction(SIGINT,  &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    uint64_t lo0 = RSTART > RESUME_FROM ? RSTART : RESUME_FROM;

    if (BRUTE) {
        if (!QUIET) fprintf(stderr, "# brute-force mode\n");
        brute_range(lo0, REND);
        if (!QUIET) fprintf(stderr, "\r%90s\r", "");
        int completed = BRUTE_POS > REND;
        if (CKPT_PATH) {
            if (completed) unlink(CKPT_PATH);
            else {
                ckpt_write(BRUTE_POS);
                fprintf(stderr, "# checkpoint saved to %s -- "
                                "rerun the same command to resume\n", CKPT_PATH);
            }
        }
        print_stats(completed);
        return completed ? 0 : 130;
    }

    uint32_t npr, *pr = make_primes(SIEVE_B, &npr);
    build_wheel(pr, npr, WHEEL_MAXRES, 1ULL << 40);
    build_prime_tables(pr, npr);

    /* block geometry */
    uint64_t lo = RSTART > (uint64_t)SIEVE_B ? RSTART : (uint64_t)SIEVE_B + 1;
    if (lo < RESUME_FROM) lo = RESUME_FROM;
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
                        "%.1f MB/thread + %.0f MB shared, span %.3g/block\n",
                SIEVE_B, NSP, JLEN,
                (double)NRES * JLEN / 8 / 1048576.0,
                (double)NRES * NSP * sizeof *RI / 1048576.0, (double)SPAN);
        fprintf(stderr, "# estimate   singular series %.4g, expected chains %s%.4g\n",
                S, REQ_CONSEC ? "(ignoring consecutivity) " : "", E);
        fprintf(stderr, "# threads    %d\n", NTHREADS);
        if (CKPT_PATH)
            fprintf(stderr, "# checkpoint %s (every %.0fs)\n", CKPT_PATH, CKPT_EVERY);
    }

    /* small numbers: any chain member could equal a sieving prime, so brute */
    if (lo0 <= (uint64_t)SIEVE_B) {
        brute_range(lo0, REND < (uint64_t)SIEVE_B ? REND : (uint64_t)SIEVE_B);
        if (stop_requested) {
            if (CKPT_PATH) {
                ckpt_write(BRUTE_POS);
                fprintf(stderr, "# checkpoint saved to %s -- "
                                "rerun the same command to resume\n", CKPT_PATH);
            }
            print_stats(0);
            return 130;
        }
    }

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
        double last_ck = elapsed_now();
        for (;;) {
            struct timespec ts = {1, 0};
            nanosleep(&ts, NULL);              /* wakes early on ctrl-c */
            uint64_t done = atomic_load(&stat_blocks);
            if (stop_requested || done >= total) break;
            double el = elapsed_now();
            if (!QUIET) {
                char cb[32], eb[24], rb[24];
                fmt_u64(atomic_load(&stat_chains), cb);
                fmt_dur(ELAPSED_PREV + el, eb);
                fmt_dur(done ? el * (double)(total - done) / done : 0.0, rb);
                fprintf(stderr, "\r# %llu/%llu blocks  %.1f%%  %s chains  "
                                "%s elapsed, ~%s left   ",
                        (unsigned long long)done, (unsigned long long)total,
                        100.0 * done / total, cb, eb, rb);
            }
            if (CKPT_PATH && el - last_ck >= CKPT_EVERY) {
                last_ck = el;
                ckpt_write(sieve_resume());
            }
        }
        for (int i = 0; i < NTHREADS; i++) pthread_join(th[i], NULL);
        if (!QUIET) fprintf(stderr, "\r%90s\r", "");
    }

    int completed = SIEVE_LO > SIEVE_HI || sieve_resume() > REND;
    if (CKPT_PATH) {
        if (completed) unlink(CKPT_PATH);
        else {
            ckpt_write(sieve_resume());
            fprintf(stderr, "# checkpoint saved to %s -- "
                            "rerun the same command to resume\n", CKPT_PATH);
        }
    }
    print_stats(completed);
    return completed ? 0 : 130;
}
