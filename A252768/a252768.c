/*
 * a252768.c
 *
 * Compute and extend OEIS A252768:
 *
 *   "Primes p with property that the sum of the k-th powers of the successive
 *    gaps between primes <= p are prime numbers for k = 1 to n."
 *
 * Known terms (n = 1..7):  5, 5, 13, 14593, 372313, 2315773, 541613713.
 *
 * Definitions
 * -----------
 * Let 2 = q_1 < q_2 < ... < q_m = p be the primes up to p and g_i = q_(i+1) - q_i
 * their gaps.  Write S_k(p) = sum_i g_i^k.  The "depth" of p is the largest n
 * such that S_1(p), ..., S_n(p) are all prime, and a(n) is the smallest prime of
 * depth >= n.  Since S_1(p) = p - 2, every term is the larger member of a twin
 * prime pair (A006512), and a(n) is non-decreasing in n.
 *
 * Example: p = 13, gaps 1,2,2,4,2: S_1 = 11, S_2 = 29, S_3 = 89 are prime and
 * S_4 = 1+16+16+256+16 = 305 = 5*61 is not, so 13 has depth 3 and a(3) = 13.
 *
 * Method
 * ------
 * The search range is cut into chunks.  Worker threads take chunks from an
 * atomic counter, sieve them with primesieve and keep the gaps (one byte per
 * prime).  The gap histogram of a chunk gives its power-sum totals T_k, and an
 * ordered hand-off between the threads (chunk i waits until chunk i-1 has
 * published its totals) turns them into the running sums O_k at each chunk
 * start.  After that hand-off the chunks are checked fully in parallel, so the
 * range is sieved exactly once, memory stays bounded and no barrier is needed.
 *
 * Checking a chunk: S_1(p) = p - 2 is prime exactly when the gap into p is 2,
 * so only twin primes are examined.  S_2, S_3, S_4 are maintained incrementally
 * and tested in turn; S_5, S_6, ... are evaluated from the running gap
 * histogram (sum over g of count(g)*g^k) only for the rare p where S_2..S_4 are
 * all prime.  All sums are 128-bit with overflow detection.
 *
 * Primality: below 2^64 deterministic Miller-Rabin (bases 2, 325, 9375, 28178,
 * 450775, 9780504, 1795265022) in Montgomery arithmetic; above 2^64 Miller-
 * Rabin to the 24 prime bases 2..89 in 128-bit Montgomery arithmetic
 * (deterministic below 3.3e24, a strong probable-prime test beyond that).
 * When built with GMP, 'verify' cross-checks each sum with mpz_probab_prime_p.
 *
 * Congruences: g^k mod q depends only on k mod (q-1), so S_k = S_k' (mod q)
 * whenever k = k' (mod q-1).  Once S_1 and S_2 are prime, no S_k is divisible
 * by 3, S_5 = S_1 and S_6 = S_2 (mod 5), S_7 = S_1 (mod 7), and so on.  The
 * conditions are therefore positively correlated; the density of primes of
 * depth >= n falls off roughly like prod_k 2 c_k / ln S_k(p) with
 * S_k(p) ~ k! p (ln p)^(k-1).  Calibrated on a(4..7) this predicts a(8) in
 * the 10^10..10^12 range and a(9) below about 10^13.
 *
 * Usage
 * -----
 *   a252768 scan [START] END [-t T] [-c CHUNK] [-k KMAX] [-r DEPTH] [-S STATE] [-i SECS] [-q]
 *       Exhaustive scan of the primes in [START, END] (START defaults to 0; with
 *       START > 0 the sums up to START are first computed in a totals-only pass).
 *       Prints every prime of depth >= DEPTH (default 7) as it is found, then a
 *       table of a(n).  -S FILE keeps a checkpoint (every SECS seconds, default
 *       60, and on Ctrl-C); rerunning the same command resumes from it.
 *   a252768 verify P [-k KMAX] [-t T]
 *       Recompute S_1..S_KMAX at the prime P from scratch and test each one
 *       (also with GMP when compiled in).  -t 1 (default) is a plain
 *       single-threaded accumulation; -t T uses the parallel totals pipeline.
 *   a252768 selftest [-t T]
 *
 * Numbers may be written as decimal, 2^k, 10^k, 1e12, or X+Y / X-Y of those.
 *
 * Build:  cc -O2 -std=gnu11 -pthread -I/opt/homebrew/include a252768.c \
 *            -L/opt/homebrew/lib -lprimesieve -lm [-DHAVE_GMP -lgmp] -o a252768
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <errno.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include <time.h>
#include <inttypes.h>
#include <primesieve.h>
#ifdef HAVE_GMP
#include <gmp.h>
#endif

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint8_t  u8;
typedef unsigned __int128 u128;

#define KLIM    16      /* largest k the power tables support */
#define MAXG    4096    /* every prime gap must be < MAXG (largest gap below 2^64 is 1550) */
#define MARGIN  4096    /* sieve this far below a chunk to find the preceding prime */
#define RING    4096    /* chunk bookkeeping ring; bounds how far threads may run ahead */
#define NKNOWN  7

static const u64 KNOWN[NKNOWN + 1] = { 0, 5, 5, 13, 14593, 372313, 2315773, 541613713 };

static int  g_kmax   = 12;
static int  g_report = 7;
static bool g_quiet  = false;
static bool stderr_tty;
static volatile sig_atomic_t g_stop = 0;
static pthread_mutex_t out_mu = PTHREAD_MUTEX_INITIALIZER;

/* ------------------------------------------------------------------ */
/* Generic helpers                                                     */
/* ------------------------------------------------------------------ */

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + ts.tv_nsec * 1e-9;
}

__attribute__((noreturn)) static void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fputs(stderr_tty ? "\r\033[Ka252768: " : "a252768: ", stderr);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

static void note(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    pthread_mutex_lock(&out_mu);
    if (stderr_tty) fputs("\r\033[K", stderr);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    pthread_mutex_unlock(&out_mu);
    va_end(ap);
}

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) die("out of memory");
    return p;
}

static void *xcalloc(size_t n, size_t sz)
{
    void *p = calloc(n ? n : 1, sz);
    if (!p) die("out of memory");
    return p;
}

static void *xrealloc(void *p, size_t n)
{
    p = realloc(p, n ? n : 1);
    if (!p) die("out of memory");
    return p;
}

/* decimal string of a 128-bit value; buf must hold >= 40 chars */
static char *u128_str(u128 v, char *buf)
{
    char tmp[48];
    int i = 0, j = 0;
    if (v == 0) tmp[i++] = '0';
    while (v) { tmp[i++] = (char)('0' + (int)(v % 10)); v /= 10; }
    while (i) buf[j++] = tmp[--i];
    buf[j] = 0;
    return buf;
}

static int bitlen128(u128 v)
{
    u64 hi = (u64)(v >> 64), lo = (u64)v;
    if (hi) return 128 - __builtin_clzll(hi);
    return lo ? 64 - __builtin_clzll(lo) : 0;
}

static bool parse_dec(const char *s, size_t len, u128 *out)
{
    u128 r = 0;
    if (len == 0) return false;
    for (size_t i = 0; i < len; i++) {
        if (s[i] < '0' || s[i] > '9') return false;
        unsigned dgt = (unsigned)(s[i] - '0');
        if (r > (((u128)0 - 1) - dgt) / 10) return false;
        r = r * 10 + dgt;
    }
    *out = r;
    return true;
}

/* accepts "123", "2^64", "10^12", "1e15", and X+Y / X-Y of those ("2^54-32") */
static bool parse_u128(const char *s, u128 *out)
{
    const char *pm = strpbrk(s, "+-");
    if (pm && pm != s) {
        char head[128];
        size_t hl = (size_t)(pm - s);
        u128 a, b;
        if (hl >= sizeof head) return false;
        memcpy(head, s, hl);
        head[hl] = 0;
        if (!parse_u128(head, &a) || !parse_u128(pm + 1, &b)) return false;
        if (*pm == '+') { if (a > ((u128)0 - 1) - b) return false; *out = a + b; }
        else            { if (b > a) return false; *out = a - b; }
        return true;
    }
    const char *caret = strchr(s, '^');
    const char *ee = strpbrk(s, "eE");
    if (caret) {
        u128 b, e, r = 1;
        if (!parse_dec(s, (size_t)(caret - s), &b)) return false;
        if (!parse_dec(caret + 1, strlen(caret + 1), &e) || e > 400) return false;
        for (unsigned i = 0; i < (unsigned)e; i++) {
            if (b != 0 && r > ((u128)0 - 1) / b) return false;
            r *= b;
        }
        *out = r;
        return true;
    }
    if (ee) {
        u128 m, e;
        if (!parse_dec(s, (size_t)(ee - s), &m)) return false;
        if (!parse_dec(ee + 1, strlen(ee + 1), &e) || e > 40) return false;
        for (unsigned i = 0; i < (unsigned)e; i++) {
            if (m > ((u128)0 - 1) / 10) return false;
            m *= 10;
        }
        *out = m;
        return true;
    }
    return parse_dec(s, strlen(s), out);
}

static u64 arg_u64(const char *s, const char *what)
{
    u128 v;
    if (!parse_u128(s, &v) || v > UINT64_MAX) die("bad %s: '%s'", what, s);
    return (u64)v;
}

static int arg_int(const char *s, const char *what, int lo, int hi)
{
    char *end;
    long v = strtol(s, &end, 10);
    if (*s == 0 || *end != 0 || v < lo || v > hi) die("bad %s: '%s' (allowed %d..%d)", what, s, lo, hi);
    return (int)v;
}

static const char *fmt_eng(double v, char *buf)
{
    if (v < 1e6) sprintf(buf, "%.0f", v);
    else sprintf(buf, "%.3g", v);
    return buf;
}

static const char *fmt_hms(double s, char *buf)
{
    long t = s > 0 ? (long)(s + 0.5) : 0;
    if (t >= 360000L) sprintf(buf, "%ldd%02ldh", t / 86400, (t / 3600) % 24);
    else sprintf(buf, "%02ld:%02ld:%02ld", t / 3600, (t / 60) % 60, t % 60);
    return buf;
}

/* ------------------------------------------------------------------ */
/* Power tables: PW[k][g] = g^k (PWOVF set where g^k >= 2^128)         */
/* ------------------------------------------------------------------ */

static u128 PW[KLIM + 1][MAXG];
static u8   PWOVF[KLIM + 1][MAXG];

static void init_pow_tables(void)
{
    for (unsigned g = 0; g < MAXG; g++) {
        PW[0][g] = 1;
        PWOVF[0][g] = 0;
        for (int k = 1; k <= KLIM; k++) {
            if (PWOVF[k - 1][g] || __builtin_mul_overflow(PW[k - 1][g], (u128)g, &PW[k][g])) {
                PWOVF[k][g] = 1;
                PW[k][g] = 0;
            } else {
                PWOVF[k][g] = 0;
            }
        }
    }
}

/* ------------------------------------------------------------------ */
/* Primality: Montgomery Miller-Rabin, 64-bit deterministic            */
/* ------------------------------------------------------------------ */

static const u32 SMALLP[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                              53, 59, 61, 67, 71, 73, 79, 83, 89, 97 };
#define NSMALLP ((int)(sizeof SMALLP / sizeof SMALLP[0]))

typedef struct { u64 n, ninv, one; } mont64;

static inline void m64_init(mont64 *m, u64 n)
{
    u64 x = n;                              /* Newton iteration: n^-1 mod 2^64, n odd */
    for (int i = 0; i < 6; i++) x *= 2 - n * x;
    m->n = n;
    m->ninv = (u64)0 - x;
    m->one = ((u64)0 - n) % n;              /* 2^64 mod n */
}

static inline u64 m64_redc(const mont64 *m, u128 T)   /* T < n * 2^64 */
{
    u64 mm = (u64)T * m->ninv;
    u128 mn = (u128)mm * m->n;
    u128 s = (T >> 64) + (mn >> 64) + ((u64)T != 0);
    if (s >= m->n) s -= m->n;
    return (u64)s;
}

static inline u64 m64_mul(const mont64 *m, u64 a, u64 b) { return m64_redc(m, (u128)a * b); }

static inline u64 m64_add(const mont64 *m, u64 a, u64 b)
{
    u64 s = a + b;
    if (s < a || s >= m->n) s -= m->n;
    return s;
}

static inline u64 m64_to(const mont64 *m, u64 x) { return (u64)(((u128)x << 64) % m->n); }

/* strong probable prime test to a base given in Montgomery form; n - 1 = d 2^s */
static bool m64_sprp(const mont64 *m, u64 am, u64 d, int s)
{
    u64 x = m->one, b = am, nm1 = m->n - m->one;
    for (u64 e = d; e; e >>= 1) {
        if (e & 1) x = m64_mul(m, x, b);
        b = m64_mul(m, b, b);
    }
    if (x == m->one || x == nm1) return true;
    for (int i = 1; i < s; i++) {
        x = m64_mul(m, x, x);
        if (x == nm1) return true;
        if (x == m->one) return false;
    }
    return false;
}

/* base 2, left to right: multiplying by the base is a modular doubling */
static bool m64_sprp2(const mont64 *m, u64 d, int s)
{
    u64 two = m64_add(m, m->one, m->one), nm1 = m->n - m->one;
    u64 x = two;
    for (int i = 62 - __builtin_clzll(d); i >= 0; i--) {
        x = m64_mul(m, x, x);
        if ((d >> i) & 1) x = m64_add(m, x, x);
    }
    if (x == m->one || x == nm1) return true;
    for (int i = 1; i < s; i++) {
        x = m64_mul(m, x, x);
        if (x == nm1) return true;
        if (x == m->one) return false;
    }
    return false;
}

static bool is_prime64(u64 n)
{
    if (n < 2) return false;
    if (n < 4) return true;
    if (!(n & 1)) return false;
    if (n < 1000) {
        for (int i = 1; i < NSMALLP; i++) {
            if (n == SMALLP[i]) return true;
            if (n % SMALLP[i] == 0) return false;
        }
        return true;
    }
    u32 r = (u32)(n % 3234846615u);         /* 3*5*7*11*13*17*19*23*29 */
    if (!(r % 3) || !(r % 5) || !(r % 7) || !(r % 11) || !(r % 13) ||
        !(r % 17) || !(r % 19) || !(r % 23) || !(r % 29)) return false;
    r = (u32)(n % 95041567u);               /* 31*37*41*43*47 */
    if (!(r % 31) || !(r % 37) || !(r % 41) || !(r % 43) || !(r % 47)) return false;

    mont64 m;
    m64_init(&m, n);
    u64 d = n - 1;
    int s = __builtin_ctzll(d);
    d >>= s;
    if (!m64_sprp2(&m, d, s)) return false;
    if (n < 2047) return true;
    /* with base 2 these bases are deterministic for all n < 2^64 */
    static const u64 B[6] = { 325, 9375, 28178, 450775, 9780504, 1795265022 };
    for (int i = 0; i < 6; i++) {
        u64 a = B[i] % n;
        if (a == 0) continue;
        if (!m64_sprp(&m, m64_to(&m, a), d, s)) return false;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* Primality: 128-bit Montgomery Miller-Rabin                          */
/* ------------------------------------------------------------------ */

typedef struct { u128 n, ninv, one, r2; } mont128;

static inline void mul128_full(u128 a, u128 b, u128 *hi, u128 *lo)
{
    u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64);
    u128 p00 = (u128)a0 * b0, p01 = (u128)a0 * b1, p10 = (u128)a1 * b0, p11 = (u128)a1 * b1;
    u128 mid = (p00 >> 64) + (u64)p01 + (u64)p10;       /* < 3 * 2^64 */
    *lo = (mid << 64) | (u64)p00;
    *hi = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
}

static inline u128 m128_add(const mont128 *m, u128 a, u128 b)
{
    u128 s = a + b;
    if (s < a || s >= m->n) s -= m->n;
    return s;
}

static inline u128 m128_redc(const mont128 *m, u128 hi, u128 lo)   /* hi < n */
{
    u128 mm = lo * m->ninv;
    u128 h2, l2, t;
    mul128_full(mm, m->n, &h2, &l2);
    (void)l2;                                           /* lo + l2 is 0 or 2^128 */
    bool c1 = __builtin_add_overflow(hi, h2, &t);
    bool c2 = __builtin_add_overflow(t, (u128)(lo != 0), &t);
    if (c1 || c2 || t >= m->n) t -= m->n;
    return t;
}

static inline u128 m128_mul(const mont128 *m, u128 a, u128 b)
{
    u128 hi, lo;
    mul128_full(a, b, &hi, &lo);
    return m128_redc(m, hi, lo);
}

static void m128_init(mont128 *m, u128 n)
{
    u128 x = n;
    for (int i = 0; i < 7; i++) x *= 2 - n * x;
    m->n = n;
    m->ninv = (u128)0 - x;
    m->one = ((u128)0 - n) % n;             /* 2^128 mod n */
    u128 r = m->one;
    for (int i = 0; i < 128; i++) r = m128_add(m, r, r);
    m->r2 = r;                              /* 2^256 mod n */
}

static inline u128 m128_to(const mont128 *m, u128 x) { return m128_mul(m, x, m->r2); }   /* x < n */

static inline int ctz128(u128 v)
{
    u64 lo = (u64)v;
    return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((u64)(v >> 64));
}

static bool m128_sprp(const mont128 *m, u128 am, u128 d, int s)
{
    u128 x = m->one, b = am, nm1 = m->n - m->one;
    for (u128 e = d; e; e >>= 1) {
        if (e & 1) x = m128_mul(m, x, b);
        b = m128_mul(m, b, b);
    }
    if (x == m->one || x == nm1) return true;
    for (int i = 1; i < s; i++) {
        x = m128_mul(m, x, x);
        if (x == nm1) return true;
        if (x == m->one) return false;
    }
    return false;
}

/* the first 13 primes are deterministic below 3.3e24 (Sorenson-Webster); 11 more beyond */
static const u32 BASES128[24] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
                                  43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89 };

static bool is_prime128(u128 n)
{
    if (n < ((u128)1 << 64)) return is_prime64((u64)n);
    if (!(n & 1)) return false;
    for (int i = 1; i < NSMALLP; i++)
        if (n % SMALLP[i] == 0) return false;
    mont128 m;
    m128_init(&m, n);
    u128 d = n - 1;
    int s = ctz128(d);
    d >>= s;
    for (int i = 0; i < 24; i++)
        if (!m128_sprp(&m, m128_to(&m, BASES128[i]), d, s)) return false;
    return true;
}

#ifdef HAVE_GMP
static void mpz_set_u128(mpz_t z, u128 v)
{
    mpz_set_ui(z, (unsigned long)(u64)(v >> 64));
    mpz_mul_2exp(z, z, 64);
    mpz_add_ui(z, z, (unsigned long)(u64)v);
}

static int gmp_isprime128(u128 n)   /* 2 prime, 1 probable prime, 0 composite */
{
    mpz_t z;
    mpz_init(z);
    mpz_set_u128(z, n);
    int r = mpz_probab_prime_p(z, 40);
    mpz_clear(z);
    return r;
}
#endif

/* ------------------------------------------------------------------ */
/* Chunks: sieve a range, keep its gaps and gap histogram              */
/* ------------------------------------------------------------------ */

typedef struct {
    u8  *gb;  size_t gb_cap, n;       /* gap bytes: g/2 for even g <= 508, 255 = escaped */
    u32 *esc; size_t esc_cap, esc_n;  /* escaped gaps (1, or >= 510) in order */
    u64  cnt[MAXG];                   /* gap histogram of the chunk */
    u64  lo, hi, p_prev, p_last, nprimes;
    unsigned maxg;
    bool has_prev;                    /* false only for the chunk that contains 2 */
} chunk_t;

static inline void chunk_push(chunk_t *c, u64 q, u64 g)
{
    if (g >= MAXG) die("prime gap %" PRIu64 " before %" PRIu64 " exceeds MAXG", g, q);
    if (c->n == c->gb_cap) {
        c->gb_cap = c->gb_cap * 2 + 65536;
        c->gb = xrealloc(c->gb, c->gb_cap);
    }
    if (!(g & 1) && g <= 508) {
        c->gb[c->n++] = (u8)(g >> 1);
    } else {
        if (c->esc_n == c->esc_cap) {
            c->esc_cap = c->esc_cap * 2 + 64;
            c->esc = xrealloc(c->esc, c->esc_cap * sizeof(u32));
        }
        c->esc[c->esc_n++] = (u32)g;
        c->gb[c->n++] = 255;
    }
    c->cnt[g]++;
    if (g > c->maxg) c->maxg = (unsigned)g;
}

/* primes q in [lo, hi) and the gap into each of them */
static void sieve_chunk(primesieve_iterator *it, chunk_t *c, u64 lo, u64 hi)
{
    c->lo = lo; c->hi = hi;
    c->n = 0; c->esc_n = 0; c->nprimes = 0; c->maxg = 0;
    memset(c->cnt, 0, sizeof c->cnt);

    double lnx = log((double)lo + (double)(hi - lo) * 0.5 + 3.0);
    size_t need = (size_t)((double)(hi - lo) / (lnx > 1.0 ? lnx : 1.0) * 1.25) + 65536;
    if (c->gb_cap < need) {
        c->gb_cap = need;
        c->gb = xrealloc(c->gb, need);
    }

    u64 margin = MARGIN;
    for (;;) {
        u64 start = lo > margin ? lo - margin : 0;
        primesieve_jump_to(it, start, hi);
        u64 pp = 0, q;
        bool have = false;
        while ((q = primesieve_next_prime(it)) < lo) { pp = q; have = true; }
        if (!have && lo > 2) {
            margin *= 4;
            if (margin > ((u64)1 << 24)) die("no prime within 2^24 below %" PRIu64 "?!", lo);
            continue;
        }
        c->has_prev = have;
        c->p_prev = pp;
        while (q < hi) {
            if (have) chunk_push(c, q, q - pp);
            else have = true;                   /* q = 2, the first prime, has no gap */
            pp = q;
            c->nprimes++;
            q = primesieve_next_prime(it);
        }
        if (it->is_error) die("primesieve error in [%" PRIu64 ", %" PRIu64 ")", lo, hi);
        c->p_last = pp;
        return;
    }
}

/* T_k = sum over the chunk's gaps of g^k, k = 1..g_kmax */
static void chunk_totals(const chunk_t *c, u128 T[], bool Tovf[])
{
    for (int k = 1; k <= g_kmax; k++) {
        u128 s = 0;
        bool ovf = false;
        for (unsigned g = 1; g <= c->maxg; g++) {
            u64 n = c->cnt[g];
            u128 prod;
            if (!n) continue;
            if (PWOVF[k][g] || __builtin_mul_overflow((u128)n, PW[k][g], &prod) ||
                __builtin_add_overflow(s, prod, &s)) { ovf = true; break; }
        }
        T[k] = ovf ? 0 : s;
        Tovf[k] = ovf;
    }
}

/* ------------------------------------------------------------------ */
/* Scan state shared by the worker threads                             */
/* ------------------------------------------------------------------ */

typedef struct { u64 p; int depth; } hit_t;

typedef struct {
    u64  idx;
    u64  primes;                      /* primes below this chunk */
    u128 O[KLIM + 1];
    bool ovf[KLIM + 1];
} rec_t;

typedef struct scan scan_t;

typedef struct {
    int id;
    pthread_t th;
    scan_t *sc;
    primesieve_iterator it;
    chunk_t c;
    u64 hist[MAXG];
    u64 best[KLIM + 1];               /* per chunk: smallest p of exact depth d */
    u64 depth_count[KLIM + 1];
    u64 twins;
    hit_t *hits; size_t nhits, hits_cap;
} worker_t;

struct scan {
    /* parameters */
    u64  origin, start, end, chunk, nchunks;
    bool check, print_hits, verbose;
    int  report;
    u128 O0[KLIM + 1]; bool O0ovf[KLIM + 1];
    u64  primes0;
    /* runtime */
    atomic_uint_fast64_t next;
    atomic_int finished;
    pthread_mutex_t mu;
    pthread_cond_t cv;
    u64  pub;                         /* chunk whose start offsets are in O */
    u128 O[KLIM + 1]; bool Oovf[KLIM + 1]; u64 ovf_at[KLIM + 1];
    u64  primes_pub;
    u64  cf;                          /* completion frontier */
    u8   done[RING];
    rec_t rec[RING];
    /* results */
    u64  best[KLIM + 1];
    u64  depth_count[KLIM + 1];
    u64  primes, twins, chunks_done;
    hit_t *hits; size_t nhits, hits_cap;
    double t0, t1;
};

static void report_hit(u64 p, int depth, const u128 S[], int known, bool capped)
{
    char buf[48];
    pthread_mutex_lock(&out_mu);
    if (stderr_tty) fputs("\r\033[K", stderr);
    printf("p = %" PRIu64 "  n = %d  S_1..S_%d =", p, depth, depth);
    for (int k = 1; k <= depth; k++) printf(" %s", u128_str(S[k], buf));
    if (known > depth) printf("  S_%d = %s (composite)", depth + 1, u128_str(S[depth + 1], buf));
    else if (capped) printf("  [S_%d exceeds 2^128, not tested]", depth + 1);
    else printf("  [k_max = %d reached]", depth);
    printf("\n");
    fflush(stdout);
    pthread_mutex_unlock(&out_mu);
}

static void check_chunk(worker_t *w, const u128 O[], const bool Oovf[])
{
    scan_t *sc = w->sc;
    chunk_t *c = &w->c;
    char b1[48];
    int keff = g_kmax;
    for (int k = 1; k <= g_kmax; k++) if (Oovf[k]) { keff = k - 1; break; }
    if (keff < 4) die("running sums overflowed 128 bits already for k <= 4 (unexpected)");
    if (c->has_prev && O[1] != (u128)(c->p_prev - 2))
        die("internal error: S_1 offset %s at chunk start %" PRIu64 " but previous prime is %" PRIu64,
            u128_str(O[1], b1), c->lo, c->p_prev);

    memset(w->hist, 0, sizeof w->hist);
    u128 S2 = O[2], S3 = O[3], S4 = O[4];
    u64 p = c->has_prev ? c->p_prev : 2;    /* the gaps start at the prime 2 */
    size_t ei = 0;
    unsigned maxg = 0;
    const u8 *gb = c->gb;
    const u32 *esc = c->esc;

    for (size_t i = 0, n = c->n; i < n; i++) {
        u64 g = gb[i] == 255 ? (u64)esc[ei++] : ((u64)gb[i] << 1);
        p += g;
        S2 += PW[2][g];
        S3 += PW[3][g];
        S4 += PW[4][g];
        w->hist[g]++;
        if (g > maxg) maxg = (unsigned)g;
        if (g != 2) continue;               /* S_1 = p - 2 is prime iff the gap into p is 2 */
        w->twins++;
        if (!is_prime128(S2)) {
            w->depth_count[1]++;
            if (!w->best[1]) w->best[1] = p;
            continue;
        }

        u128 S[KLIM + 2];
        int depth = 2, known = 3;
        bool capped = false;
        S[1] = p - 2; S[2] = S2; S[3] = S3;
        if (is_prime128(S3)) {
            depth = 3; S[4] = S4; known = 4;
            if (is_prime128(S4)) {
                depth = 4;
                for (int k = 5; k <= keff; k++) {
                    u128 s = O[k];
                    bool ovf = false;
                    for (unsigned gg = 1; gg <= maxg; gg++) {
                        u64 cnt = w->hist[gg];
                        u128 prod;
                        if (!cnt) continue;
                        if (PWOVF[k][gg] || __builtin_mul_overflow((u128)cnt, PW[k][gg], &prod) ||
                            __builtin_add_overflow(s, prod, &s)) { ovf = true; break; }
                    }
                    if (ovf) { capped = true; break; }
                    S[k] = s; known = k;
                    if (!is_prime128(s)) break;
                    depth = k;
                }
                if (depth == keff && keff < g_kmax) capped = true;
            }
        }
        w->depth_count[depth]++;
        if (!w->best[depth] || p < w->best[depth]) w->best[depth] = p;
        if (depth >= sc->report) {
            if (sc->print_hits) report_hit(p, depth, S, known, capped);
            if (w->nhits == w->hits_cap) {
                w->hits_cap = w->hits_cap * 2 + 64;
                w->hits = xrealloc(w->hits, w->hits_cap * sizeof(hit_t));
            }
            w->hits[w->nhits].p = p;
            w->hits[w->nhits].depth = depth;
            w->nhits++;
        }
    }
}

static void *worker_main(void *arg)
{
    worker_t *w = arg;
    scan_t *sc = w->sc;
    chunk_t *c = &w->c;

    for (;;) {
        if (g_stop) break;
        u64 i = atomic_fetch_add(&sc->next, 1);
        if (i >= sc->nchunks) break;

        pthread_mutex_lock(&sc->mu);
        while (i >= sc->cf + RING) pthread_cond_wait(&sc->cv, &sc->mu);
        pthread_mutex_unlock(&sc->mu);

        u64 lo = sc->start + i * sc->chunk;
        u64 hi = (i + 1 == sc->nchunks) ? sc->end + 1 : lo + sc->chunk;
        sieve_chunk(&w->it, c, lo, hi);
        u128 T[KLIM + 1];
        bool Tovf[KLIM + 1];
        chunk_totals(c, T, Tovf);

        /* ordered hand-off: wait for the offsets of this chunk, publish the next */
        u128 O[KLIM + 1];
        bool Oovf[KLIM + 1];
        pthread_mutex_lock(&sc->mu);
        while (sc->pub != i) pthread_cond_wait(&sc->cv, &sc->mu);
        memcpy(O, sc->O, sizeof O);
        memcpy(Oovf, sc->Oovf, sizeof Oovf);
        rec_t *r = &sc->rec[i % RING];
        r->idx = i;
        r->primes = sc->primes_pub;
        memcpy(r->O, O, sizeof O);
        memcpy(r->ovf, Oovf, sizeof Oovf);
        for (int k = 1; k <= g_kmax; k++) {
            if (sc->Oovf[k]) continue;
            if (Tovf[k] || __builtin_add_overflow(sc->O[k], T[k], &sc->O[k])) {
                sc->Oovf[k] = true;
                sc->ovf_at[k] = lo;
                if (sc->check)
                    note("note: S_%d exceeds 2^128 from %" PRIu64 " on; depths >= %d are no longer tested", k, lo, k);
            }
        }
        sc->primes_pub += c->nprimes;
        sc->pub = i + 1;
        pthread_cond_broadcast(&sc->cv);
        pthread_mutex_unlock(&sc->mu);

        w->nhits = 0;
        w->twins = 0;
        memset(w->best, 0, sizeof w->best);
        memset(w->depth_count, 0, sizeof w->depth_count);
        if (sc->check) check_chunk(w, O, Oovf);

        pthread_mutex_lock(&sc->mu);
        sc->done[i % RING] = 1;
        while (sc->cf < sc->nchunks && sc->done[sc->cf % RING]) {
            sc->done[sc->cf % RING] = 0;
            sc->cf++;
        }
        sc->primes += c->nprimes;
        sc->twins += w->twins;
        sc->chunks_done++;
        for (int d = 1; d <= KLIM; d++) {
            sc->depth_count[d] += w->depth_count[d];
            if (w->best[d] && (!sc->best[d] || w->best[d] < sc->best[d])) sc->best[d] = w->best[d];
        }
        if (w->nhits) {
            if (sc->nhits + w->nhits > sc->hits_cap) {
                sc->hits_cap = (sc->nhits + w->nhits) * 2 + 256;
                sc->hits = xrealloc(sc->hits, sc->hits_cap * sizeof(hit_t));
            }
            memcpy(sc->hits + sc->nhits, w->hits, w->nhits * sizeof(hit_t));
            sc->nhits += w->nhits;
        }
        pthread_cond_broadcast(&sc->cv);
        pthread_mutex_unlock(&sc->mu);
    }
    atomic_fetch_add(&sc->finished, 1);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Checkpoint files                                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    int  kmax;
    u64  origin, end, lo, primes;
    u128 O[KLIM + 1];
    bool ovf[KLIM + 1];
    u64  best[KLIM + 1];
} state_t;

static void write_state(const char *path, const state_t *st)
{
    char tmp[4096], b[48];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { note("warning: cannot write %s: %s", tmp, strerror(errno)); return; }
    fprintf(f, "A252768-state 1\nkmax %d\norigin %" PRIu64 "\nend %" PRIu64 "\nlo %" PRIu64 "\nprimes %" PRIu64 "\n",
            st->kmax, st->origin, st->end, st->lo, st->primes);
    for (int k = 1; k <= st->kmax; k++)
        fprintf(f, "O %d %s %d\n", k, u128_str(st->O[k], b), st->ovf[k] ? 1 : 0);
    for (int d = 1; d <= KLIM; d++)
        if (st->best[d]) fprintf(f, "best %d %" PRIu64 "\n", d, st->best[d]);
    if (fclose(f) != 0 || rename(tmp, path) != 0)
        note("warning: cannot finish writing %s: %s", path, strerror(errno));
}

static bool read_state(const char *path, state_t *st)
{
    FILE *f = fopen(path, "r");
    if (!f) return false;
    memset(st, 0, sizeof *st);
    char line[512];
    int ver = 0;
    if (!fgets(line, sizeof line, f) || sscanf(line, "A252768-state %d", &ver) != 1 || ver != 1)
        die("%s is not an a252768 state file", path);
    while (fgets(line, sizeof line, f)) {
        char key[32], v1[64], v2[64], v3[64];
        int n = sscanf(line, "%31s %63s %63s %63s", key, v1, v2, v3);
        u128 v;
        if (n < 2) continue;
        if (!strcmp(key, "kmax")) st->kmax = atoi(v1);
        else if (!strcmp(key, "origin")) { if (!parse_u128(v1, &v)) die("bad state"); st->origin = (u64)v; }
        else if (!strcmp(key, "end"))    { if (!parse_u128(v1, &v)) die("bad state"); st->end = (u64)v; }
        else if (!strcmp(key, "lo"))     { if (!parse_u128(v1, &v)) die("bad state"); st->lo = (u64)v; }
        else if (!strcmp(key, "primes")) { if (!parse_u128(v1, &v)) die("bad state"); st->primes = (u64)v; }
        else if (!strcmp(key, "O") && n >= 4) {
            int k = atoi(v1);
            if (k < 1 || k > KLIM || !parse_u128(v2, &st->O[k])) die("bad O line in %s", path);
            st->ovf[k] = atoi(v3) != 0;
        } else if (!strcmp(key, "best") && n >= 3) {
            int d = atoi(v1);
            if (d < 1 || d > KLIM || !parse_u128(v2, &v)) die("bad best line in %s", path);
            st->best[d] = (u64)v;
        }
    }
    fclose(f);
    if (st->kmax < 4 || st->kmax > KLIM) die("bad kmax in %s", path);
    return true;
}

static void checkpoint(scan_t *sc, const char *path)
{
    state_t st;
    memset(&st, 0, sizeof st);
    bool ok = true;
    pthread_mutex_lock(&sc->mu);
    u64 cf = sc->cf;
    if (cf >= sc->nchunks) {
        st.lo = sc->end + 1;
        memcpy(st.O, sc->O, sizeof st.O);
        memcpy(st.ovf, sc->Oovf, sizeof st.ovf);
        st.primes = sc->primes_pub;
    } else if (cf == 0) {
        st.lo = sc->start;
        memcpy(st.O, sc->O0, sizeof st.O);
        memcpy(st.ovf, sc->O0ovf, sizeof st.ovf);
        st.primes = sc->primes0;
    } else if (sc->rec[cf % RING].idx == cf) {
        st.lo = sc->start + cf * sc->chunk;
        memcpy(st.O, sc->rec[cf % RING].O, sizeof st.O);
        memcpy(st.ovf, sc->rec[cf % RING].ovf, sizeof st.ovf);
        st.primes = sc->rec[cf % RING].primes;
    } else if (sc->pub == cf) {               /* chunk cf never started: O holds its offsets */
        st.lo = sc->start + cf * sc->chunk;
        memcpy(st.O, sc->O, sizeof st.O);
        memcpy(st.ovf, sc->Oovf, sizeof st.ovf);
        st.primes = sc->primes_pub;
    } else {
        ok = false;
    }
    if (ok) {
        st.kmax = g_kmax;
        st.origin = sc->origin;
        st.end = sc->end;
        for (int d = 1; d <= KLIM; d++)
            if (sc->best[d] && sc->best[d] < st.lo) st.best[d] = sc->best[d];
    }
    pthread_mutex_unlock(&sc->mu);
    if (ok) write_state(path, &st);
}

/* ------------------------------------------------------------------ */
/* Running a scan                                                      */
/* ------------------------------------------------------------------ */

static void print_status(scan_t *sc, double t, bool final)
{
    pthread_mutex_lock(&sc->mu);
    u64 cf = sc->cf;
    u64 lo = cf >= sc->nchunks ? sc->end + 1 : sc->start + cf * sc->chunk;
    u64 primes = sc->primes, twins = sc->twins;
    size_t nh = sc->nhits;
    int bd = 0;
    u64 bp = 0;
    for (int d = KLIM; d >= 1; d--) if (sc->best[d]) { bd = d; bp = sc->best[d]; break; }
    pthread_mutex_unlock(&sc->mu);

    double done = (double)(lo - sc->start), total = (double)(sc->end - sc->start) + 1.0;
    double el = t - sc->t0, rate = el > 0 ? done / el : 0, eta = rate > 0 ? (total - done) / rate : 0;
    char b1[32], b2[32], b3[32], b4[32], b5[32], b6[32];
    pthread_mutex_lock(&out_mu);
    if (stderr_tty) fputs("\r\033[K", stderr);
    fprintf(stderr, "[%s] %5.1f%% at %s  %s/s  ETA %s", fmt_hms(el, b1), 100.0 * done / total,
            fmt_eng((double)lo, b2), fmt_eng(rate, b3), fmt_hms(eta, b4));
    if (sc->check)
        fprintf(stderr, "  primes %s  twins %s  hits %zu  best n=%d @ %" PRIu64,
                fmt_eng((double)primes, b5), fmt_eng((double)twins, b6), nh, bd, bp);
    else
        fprintf(stderr, "  (sums only)");
    if (!stderr_tty || final) fputc('\n', stderr);
    fflush(stderr);
    pthread_mutex_unlock(&out_mu);
}

static void scan_run(scan_t *sc, int nthreads, const char *statefile, double ckpt_secs)
{
    if (sc->end >= primesieve_get_max_stop()) die("END too large for primesieve");
    if (nthreads < 1) nthreads = 1;
    u64 range = sc->end >= sc->start ? sc->end - sc->start + 1 : 0;
    if (!sc->chunk) {
        u64 c = range / ((u64)nthreads * 32);
        if (c < 1000000) c = 1000000;
        if (c > 1000000000) c = 1000000000;
        sc->chunk = c;
    }
    sc->nchunks = range ? (range + sc->chunk - 1) / sc->chunk : 0;

    atomic_store(&sc->next, 0);
    atomic_store(&sc->finished, 0);
    pthread_mutex_init(&sc->mu, NULL);
    pthread_cond_init(&sc->cv, NULL);
    sc->pub = 0;
    sc->cf = 0;
    memcpy(sc->O, sc->O0, sizeof sc->O);
    memcpy(sc->Oovf, sc->O0ovf, sizeof sc->Oovf);
    memset(sc->ovf_at, 0, sizeof sc->ovf_at);
    sc->primes_pub = sc->primes0;
    memset(sc->done, 0, sizeof sc->done);
    for (int i = 0; i < RING; i++) sc->rec[i].idx = UINT64_MAX;
    sc->primes = sc->twins = sc->chunks_done = 0;
    sc->nhits = 0;
    sc->t0 = now();

    worker_t *ws = xcalloc((size_t)nthreads, sizeof(worker_t));
    for (int i = 0; i < nthreads; i++) {
        ws[i].id = i;
        ws[i].sc = sc;
        primesieve_init(&ws[i].it);
        if (pthread_create(&ws[i].th, NULL, worker_main, &ws[i]) != 0) die("pthread_create failed");
    }

    double last_status = sc->t0, last_ckpt = sc->t0;
    while (atomic_load(&sc->finished) < nthreads) {
        usleep(100000);
        double t = now();
        if (sc->verbose && ((stderr_tty && t - last_status >= 0.5) || (!stderr_tty && t - last_status >= 60))) {
            print_status(sc, t, false);
            last_status = t;
        }
        if (statefile && t - last_ckpt >= ckpt_secs) {
            checkpoint(sc, statefile);
            last_ckpt = t;
        }
    }
    for (int i = 0; i < nthreads; i++) {
        pthread_join(ws[i].th, NULL);
        primesieve_free_iterator(&ws[i].it);
        free(ws[i].c.gb);
        free(ws[i].c.esc);
        free(ws[i].hits);
    }
    free(ws);
    sc->t1 = now();
    if (sc->verbose) print_status(sc, sc->t1, true);
    if (statefile) checkpoint(sc, statefile);
}

static scan_t *scan_new(void)
{
    scan_t *sc = xcalloc(1, sizeof *sc);
    sc->report = g_report;
    return sc;
}

static void scan_free(scan_t *sc)
{
    free(sc->hits);
    free(sc);
}

static u64 a_from_best(const u64 best[], int n)
{
    u64 a = 0;
    for (int d = n; d <= KLIM; d++)
        if (best[d] && (!a || best[d] < a)) a = best[d];
    return a;
}

static int cmp_hit(const void *x, const void *y)
{
    const hit_t *a = x, *b = y;
    return a->p < b->p ? -1 : a->p > b->p;
}

/* ------------------------------------------------------------------ */
/* Direct (single-threaded) computation of the gap histogram up to P   */
/* ------------------------------------------------------------------ */

static u64 direct_hist(u64 P, u64 cnt[MAXG])
{
    primesieve_iterator it;
    primesieve_init(&it);
    primesieve_jump_to(&it, 0, P);
    memset(cnt, 0, sizeof(u64) * MAXG);
    u64 pp = 0, q, nprimes = 0;
    while ((q = primesieve_next_prime(&it)) <= P) {
        if (pp) {
            u64 g = q - pp;
            if (g >= MAXG) die("prime gap %" PRIu64 " exceeds MAXG", g);
            cnt[g]++;
        }
        pp = q;
        nprimes++;
    }
    if (it.is_error) die("primesieve error");
    primesieve_free_iterator(&it);
    if (pp != P) die("%" PRIu64 " is not prime (largest prime <= it is %" PRIu64 ")", P, pp);
    return nprimes;
}

static void sums_from_hist(const u64 cnt[MAXG], int kmax, u128 S[], bool ovf[])
{
    for (int k = 1; k <= kmax; k++) {
        u128 s = 0;
        bool o = false;
        for (unsigned g = 1; g < MAXG; g++) {
            u128 prod;
            if (!cnt[g]) continue;
            if (PWOVF[k][g] || __builtin_mul_overflow((u128)cnt[g], PW[k][g], &prod) ||
                __builtin_add_overflow(s, prod, &s)) { o = true; break; }
        }
        S[k] = o ? 0 : s;
        ovf[k] = o;
    }
}

/* ------------------------------------------------------------------ */
/* Commands                                                            */
/* ------------------------------------------------------------------ */

static void usage(void)
{
    fputs(
"usage: a252768 scan [START] END [-t T] [-c CHUNK] [-k KMAX] [-r DEPTH] [-S STATE] [-i SECS] [-q]\n"
"       a252768 verify P [-k KMAX] [-t T]\n"
"       a252768 selftest [-t T]\n"
"\n"
"  scan     exhaustive search of the primes in [START, END] (START defaults to 0);\n"
"           prints every prime p whose first DEPTH power sums S_1(p)..S_DEPTH(p) are\n"
"           prime (default DEPTH 7) and the table of a(n) at the end.  -k sets how\n"
"           many sums are tested (default 12, max 16).  -S FILE writes a checkpoint\n"
"           every SECS seconds (default 60) and on Ctrl-C; running the same command\n"
"           again resumes from it.\n"
"  verify   recompute S_1..S_KMAX at the prime P from scratch and test each one\n"
"           (independently with GMP when compiled in).  -t T > 1 uses the parallel\n"
"           totals pipeline instead of a plain single-threaded accumulation.\n"
"  numbers  may be written as 123, 2^40, 10^12, 1e12, 2^40-1\n", stderr);
    exit(2);
}

static void on_signal(int s)
{
    (void)s;
    g_stop = 1;
}

static void print_summary(scan_t *sc, bool exhaustive_from_zero)
{
    char b1[32];
    double secs = sc->t1 - sc->t0;
    u64 total_primes = sc->primes0 + sc->primes;
    printf("# scanned [%" PRIu64 ", %" PRIu64 "]%s in %.1f s (%s numbers/s): %" PRIu64 " primes in range"
           " (%" PRIu64 " below %" PRIu64 "), %" PRIu64 " twin primes, %zu hits with n >= %d\n",
           sc->start, g_stop ? sc->start + sc->cf * sc->chunk - 1 : sc->end, g_stop ? " (interrupted)" : "",
           secs, fmt_eng((double)(sc->end - sc->start + 1) / (secs > 0 ? secs : 1), b1),
           sc->primes, sc->primes0, sc->start, sc->twins, sc->nhits, sc->report);
    printf("# primes in range by depth n:");
    for (int d = 1; d <= KLIM; d++)
        if (sc->depth_count[d]) printf("  n=%d: %" PRIu64, d, sc->depth_count[d]);
    printf("\n");
    (void)total_primes;

    int keff = g_kmax;
    for (int k = 1; k <= g_kmax; k++) if (sc->Oovf[k]) { keff = k - 1; break; }
    u64 scanned_to = g_stop ? sc->start + sc->cf * sc->chunk - 1 : sc->end;
    for (int n = 1; n <= g_kmax; n++) {
        u64 a = a_from_best(sc->best, n);
        if (a) {
            const char *tag = "";
            if (n <= NKNOWN) tag = a == KNOWN[n] ? "   (matches OEIS)" : "   ** DIFFERS FROM OEIS **";
            else tag = "   NEW";
            printf("a(%d) = %" PRIu64 "%s\n", n, a, tag);
        } else if (n > keff) {
            printf("a(%d): not tested (S_%d exceeds 2^128 from %" PRIu64 " on)\n", n, keff + 1, sc->ovf_at[keff + 1]);
        } else if (exhaustive_from_zero) {
            printf("a(%d) > %" PRIu64 "\n", n, scanned_to);
        } else {
            printf("a(%d): no prime of depth >= %d in [%" PRIu64 ", %" PRIu64 "]\n", n, n, sc->start, scanned_to);
        }
    }
    if (sc->nhits && sc->nhits <= 500) {
        qsort(sc->hits, sc->nhits, sizeof(hit_t), cmp_hit);
        printf("# hits (p, n) in increasing order:");
        for (size_t i = 0; i < sc->nhits; i++) printf(" (%" PRIu64 ", %d)", sc->hits[i].p, sc->hits[i].depth);
        printf("\n");
    }
    fflush(stdout);
}

static int cmd_scan(int argc, char **argv)
{
    u64 pos[2] = { 0, 0 };
    int npos = 0, threads = 0;
    u64 chunk = 0;
    const char *statefile = NULL;
    double ckpt = 60;

    for (int i = 0; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-t") && i + 1 < argc) threads = arg_int(argv[++i], "thread count", 1, 1024);
        else if (!strcmp(a, "-c") && i + 1 < argc) chunk = arg_u64(argv[++i], "chunk size");
        else if (!strcmp(a, "-k") && i + 1 < argc) g_kmax = arg_int(argv[++i], "KMAX", 4, KLIM);
        else if (!strcmp(a, "-r") && i + 1 < argc) g_report = arg_int(argv[++i], "report depth", 1, KLIM);
        else if (!strcmp(a, "-S") && i + 1 < argc) statefile = argv[++i];
        else if (!strcmp(a, "-i") && i + 1 < argc) ckpt = atof(argv[++i]);
        else if (!strcmp(a, "-q")) g_quiet = true;
        else if (a[0] == '-' && !isdigit((unsigned char)a[1])) die("unknown option '%s'", a);
        else if (npos < 2) pos[npos++] = arg_u64(a, "bound");
        else usage();
    }
    if (npos == 0) usage();
    if (chunk && chunk < 1000) die("chunk size too small");
    if (ckpt < 1) ckpt = 1;
    if (threads <= 0) {
        long n = sysconf(_SC_NPROCESSORS_ONLN);
        threads = n > 0 ? (int)n : 4;
    }

    u64 start = npos == 2 ? pos[0] : 0, end = pos[npos - 1], origin = start;
    scan_t *sc = scan_new();
    sc->report = g_report;
    bool resumed = false;
    state_t st;

    if (statefile && read_state(statefile, &st)) {
        if (st.kmax != g_kmax) die("%s was written with -k %d; use the same value", statefile, st.kmax);
        if (npos == 2 && pos[0] != st.lo)
            note("note: START %" PRIu64 " ignored, resuming from checkpoint at %" PRIu64, pos[0], st.lo);
        if (end < st.lo) {
            if (end == st.end) { note("%s: scan of [%" PRIu64 ", %" PRIu64 "] already complete", statefile, st.origin, st.end); return 0; }
            die("END %" PRIu64 " is below the checkpoint position %" PRIu64, end, st.lo);
        }
        start = st.lo;
        origin = st.origin;
        memcpy(sc->O0, st.O, sizeof sc->O0);
        memcpy(sc->O0ovf, st.ovf, sizeof sc->O0ovf);
        memcpy(sc->best, st.best, sizeof sc->best);
        sc->primes0 = st.primes;
        resumed = true;
        note("resuming from %s at %" PRIu64 " (%" PRIu64 " primes done)", statefile, start, st.primes);
    } else if (start > 2) {
        /* totals-only pass to get the power sums at START */
        scan_t *pre = scan_new();
        pre->origin = 0; pre->start = 0; pre->end = start - 1;
        pre->check = false; pre->verbose = !g_quiet; pre->print_hits = false;
        pre->chunk = chunk;
        if (!g_quiet) fprintf(stderr, "computing the power sums below %" PRIu64 " ...\n", start);
        scan_run(pre, threads, NULL, 0);
        if (g_stop) die("interrupted");
        memcpy(sc->O0, pre->O, sizeof sc->O0);
        memcpy(sc->O0ovf, pre->Oovf, sizeof sc->O0ovf);
        sc->primes0 = pre->primes_pub;
        scan_free(pre);
    }

    if (!chunk) {
        u64 range = end >= start ? end - start + 1 : 0;
        chunk = range / ((u64)threads * 32);
        if (chunk < 1000000) chunk = 1000000;
        if (chunk > 1000000000) chunk = 1000000000;
    }
    sc->origin = origin;
    sc->start = start;
    sc->end = end;
    sc->chunk = chunk;
    sc->check = true;
    sc->verbose = !g_quiet;
    sc->print_hits = true;

    char b[48];
    printf("# A252768: a(n) = least prime p such that S_1(p)..S_n(p) are prime, S_k(p) = sum of k-th powers of the prime gaps up to p\n");
    printf("# scanning primes in [%" PRIu64 ", %" PRIu64 "] with %d threads, chunk %" PRIu64 ", testing k <= %d, reporting depth >= %d%s\n",
           start, end, threads, chunk, g_kmax, g_report, resumed ? " (resumed)" : "");
    if (start > 2) {
        printf("# power sums at start:");
        for (int k = 1; k <= g_kmax; k++) printf(" S_%d=%s", k, sc->O0ovf[k] ? "overflow" : u128_str(sc->O0[k], b));
        printf("\n");
    }
    fflush(stdout);

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_signal;
    sa.sa_flags = SA_RESETHAND;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    scan_run(sc, threads, statefile, ckpt);
    if (g_stop) note("interrupted; all chunks below %" PRIu64 " are complete%s", start + sc->cf * sc->chunk,
                     statefile ? ", checkpoint written" : "");
    print_summary(sc, origin == 0);
    scan_free(sc);
    return 0;
}

static int cmd_verify(int argc, char **argv)
{
    u64 P = 0;
    int kmax = g_kmax, threads = 1;
    bool have_p = false;
    for (int i = 0; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-k") && i + 1 < argc) kmax = arg_int(argv[++i], "KMAX", 1, 64);
        else if (!strcmp(a, "-t") && i + 1 < argc) threads = arg_int(argv[++i], "thread count", 1, 1024);
        else if (a[0] == '-' && !isdigit((unsigned char)a[1])) die("unknown option '%s'", a);
        else { P = arg_u64(a, "prime"); have_p = true; }
    }
    if (!have_p) usage();
    if (!is_prime64(P)) die("%" PRIu64 " is not prime", P);
#ifndef HAVE_GMP
    if (kmax > KLIM) die("without GMP at most k = %d is supported", KLIM);
#endif
    if (threads > 1 && kmax > KLIM) die("with -t > 1 at most k = %d is supported", KLIM);

    double t0 = now();
    u64 cnt[MAXG];
    u64 nprimes;
    u128 S[KLIM + 1];
    bool ovf[KLIM + 1];
    bool have_cnt = false;
    if (threads <= 1) {
        nprimes = direct_hist(P, cnt);
        have_cnt = true;
        sums_from_hist(cnt, kmax <= KLIM ? kmax : KLIM, S, ovf);
    } else {
        int save = g_kmax;
        g_kmax = kmax;
        scan_t *sc = scan_new();
        sc->start = 0; sc->end = P; sc->check = false; sc->verbose = !g_quiet;
        scan_run(sc, threads, NULL, 0);
        nprimes = sc->primes_pub;
        memcpy(S, sc->O, sizeof S);
        memcpy(ovf, sc->Oovf, sizeof ovf);
        scan_free(sc);
        g_kmax = save;
    }
    double secs = now() - t0;

    char b[48];
    printf("p = %" PRIu64 " is prime; %" PRIu64 " primes <= p, %" PRIu64 " gaps (%.1f s)\n", P, nprimes, nprimes - 1, secs);
    int depth = 0;
    bool stopped = false;
#ifdef HAVE_GMP
    mpz_t z, t;
    mpz_init(z); mpz_init(t);
#endif
    for (int k = 1; k <= kmax; k++) {
        bool own_known = k <= KLIM && !ovf[k];
        bool own_prime = own_known && is_prime128(S[k]);
        int gmp_res = -1;
        char *dec = NULL;
#ifdef HAVE_GMP
        if (have_cnt) {
            mpz_set_ui(z, 0);
            for (unsigned g = 1; g < MAXG; g++) {
                if (!cnt[g]) continue;
                mpz_ui_pow_ui(t, g, (unsigned long)k);
                mpz_mul_ui(t, t, (unsigned long)cnt[g]);
                mpz_add(z, z, t);
            }
        } else if (own_known) {
            mpz_set_u128(z, S[k]);
        } else {
            mpz_set_ui(z, 0);
        }
        if (have_cnt || own_known) {
            gmp_res = mpz_probab_prime_p(z, 40);
            dec = mpz_get_str(NULL, 10, z);
            if (have_cnt && own_known) {
                mpz_t w;
                mpz_init(w);
                mpz_set_u128(w, S[k]);
                if (mpz_cmp(w, z) != 0) die("internal error: 128-bit and GMP sums differ for k = %d", k);
                mpz_clear(w);
            }
        }
#endif
        if (!dec && own_known) { dec = xmalloc(48); u128_str(S[k], dec); }
        bool prime = own_known ? own_prime : gmp_res > 0;
        if (!dec) {
            printf("S_%-2d exceeds 2^128 (not computed)\n", k);
            stopped = true;
            break;
        }
        int bits = own_known ? bitlen128(S[k]) : 0;
#ifdef HAVE_GMP
        if (!own_known) bits = (int)mpz_sizeinbase(z, 2);
#endif
        printf("S_%-2d = %s  (%d bits)  %s", k, dec, bits, prime ? "prime" : "composite");
        if (own_known) printf("  [mr]");
        if (gmp_res >= 0) printf("  [gmp: %s]", gmp_res == 2 ? "prime" : gmp_res == 1 ? "probable prime" : "composite");
        if (own_known && gmp_res >= 0 && own_prime != (gmp_res > 0)) printf("  ** DISAGREEMENT **");
        printf("\n");
        free(dec);
        if (prime && depth == k - 1) depth = k;
    }
#ifdef HAVE_GMP
    mpz_clear(z); mpz_clear(t);
#endif
    printf("depth(%" PRIu64 ") = %d: S_1..S_%d are prime%s\n", P, depth, depth,
           depth == kmax ? (stopped ? "" : " (all tested k)") : "");
    (void)b;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Selftest                                                            */
/* ------------------------------------------------------------------ */

static int st_fail;

static void expect(bool ok, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    printf("  %s  ", ok ? "ok  " : "FAIL");
    vprintf(fmt, ap);
    printf("\n");
    va_end(ap);
    if (!ok) st_fail++;
}

static u64 rng_state = 0x9E3779B97F4A7C15ull;
static u64 rng(void)
{
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    return rng_state * 0x2545F4914F6CDD1Dull;
}

/* run a quiet scan and hand back the result structure (caller frees) */
static scan_t *quiet_scan(u64 start, u64 end, u64 chunk, int threads, int report, const u128 *O0, const bool *O0ovf, u64 primes0)
{
    scan_t *sc = scan_new();
    sc->origin = 0; sc->start = start; sc->end = end; sc->chunk = chunk;
    sc->check = true; sc->verbose = false; sc->print_hits = false; sc->report = report;
    if (O0) { memcpy(sc->O0, O0, sizeof sc->O0); memcpy(sc->O0ovf, O0ovf, sizeof sc->O0ovf); sc->primes0 = primes0; }
    scan_run(sc, threads, NULL, 0);
    return sc;
}

static scan_t *quiet_totals(u64 end, u64 chunk, int threads)
{
    scan_t *sc = scan_new();
    sc->origin = 0; sc->start = 0; sc->end = end; sc->chunk = chunk;
    sc->check = false; sc->verbose = false; sc->print_hits = false;
    scan_run(sc, threads, NULL, 0);
    return sc;
}

static bool same_hits(scan_t *a, scan_t *b)
{
    if (a->nhits != b->nhits) return false;
    qsort(a->hits, a->nhits, sizeof(hit_t), cmp_hit);
    qsort(b->hits, b->nhits, sizeof(hit_t), cmp_hit);
    for (size_t i = 0; i < a->nhits; i++)
        if (a->hits[i].p != b->hits[i].p || a->hits[i].depth != b->hits[i].depth) return false;
    return true;
}

static int cmd_selftest(int argc, char **argv)
{
    int threads = 0;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) threads = arg_int(argv[++i], "thread count", 1, 1024);
        else die("unknown option '%s'", argv[i]);
    }
    if (threads <= 0) {
        long n = sysconf(_SC_NPROCESSORS_ONLN);
        threads = n > 0 ? (int)n : 4;
    }
    g_kmax = 12;
    char b1[48], b2[48];
    double t0 = now();
    printf("a252768 selftest (%d threads)\n", threads);

    /* 1. power tables */
    expect(PW[3][2] == 8 && PW[12][2] == 4096 && PW[1][MAXG - 1] == MAXG - 1 &&
           !PWOVF[12][1550] && PWOVF[13][1550] && PWOVF[11][MAXG - 1], "power tables");

    /* 2. primality on known values */
    static const struct { const char *s; bool prime; } PT[] = {
        { "2", true }, { "3", true }, { "4", false }, { "5", true }, { "9", false }, { "561", false },
        { "997", true }, { "1009", true }, { "2047", false }, { "1000003", true }, { "4294967291", true },
        { "4294967297", false }, { "3215031751", false }, { "3825123056546413051", false },
        { "18446744073709551557", true }, { "18446744073709551615", false },
        { "18446744073709551629", true },                                    /* 2^64 + 13 */
        { "318665857834031151167461", false },                               /* spsp to bases 2..37 */
        { "3317044064679887385961981", false },                              /* spsp to bases 2..41 */
        { "618970019642690137449562111", true },                             /* 2^89 - 1 */
        { "1856910058928070412348686333", false },                           /* 3 (2^89 - 1) */
        { "162259276829213363391578010288127", true },                       /* 2^107 - 1 */
        { "5316911983139663487003542222693990401", false },                  /* (2^61 - 1)^2 */
        { "170141183460469231731687303715884105727", true },                 /* 2^127 - 1 */
    };
    bool okp = true;
    for (size_t i = 0; i < sizeof PT / sizeof PT[0]; i++) {
        u128 n;
        if (!parse_u128(PT[i].s, &n)) die("selftest: bad literal");
        if (is_prime128(n) != PT[i].prime) { okp = false; printf("      wrong verdict for %s\n", PT[i].s); }
    }
    expect(okp, "Miller-Rabin on %zu known primes/composites", sizeof PT / sizeof PT[0]);

#ifdef HAVE_GMP
    {
        int bad = 0, nprime = 0;
        for (int i = 0; i < 20000; i++) {
            u64 n = rng() | 1;
            n >>= (rng() % 40);
            bool a = is_prime64(n), g = gmp_isprime128(n) > 0;
            if (a != g) bad++;
            nprime += a;
        }
        expect(bad == 0, "64-bit Miller-Rabin agrees with GMP on 20000 random odd numbers (%d primes)", nprime);
        bad = 0; nprime = 0;
        for (int i = 0; i < 3000; i++) {
            int bits = 65 + (int)(rng() % 63);
            u128 n = ((u128)rng() << 64) | rng();
            n >>= (128 - bits);
            n |= ((u128)1 << (bits - 1)) | 1;
            bool a = is_prime128(n), g = gmp_isprime128(n) > 0;
            if (a != g) bad++;
            nprime += a;
        }
        expect(bad == 0, "128-bit Miller-Rabin agrees with GMP on 3000 random odd numbers of 65..127 bits (%d primes)", nprime);
        bad = 0;
        for (int i = 0; i < 300; i++) {
            u64 p = (rng() >> 20) | 1, q = (rng() >> 24) | 1;
            while (!is_prime64(p)) p += 2;
            while (!is_prime64(q)) q += 2;
            u128 n = (u128)p * q;
            if (is_prime128(n) || gmp_isprime128(n) > 0) bad++;
        }
        expect(bad == 0, "products of two random primes (40-44 bits each) are recognised as composite");
    }
#endif

    /* 3. the OEIS example and a(7) by direct accumulation */
    {
        u64 cnt[MAXG];
        u128 S[KLIM + 1];
        bool ovf[KLIM + 1];
        direct_hist(13, cnt);
        sums_from_hist(cnt, 4, S, ovf);
        expect(S[1] == 11 && S[2] == 29 && S[3] == 89 && S[4] == 305, "p = 13: S_1..S_4 = 11, 29, 89, 305");
        direct_hist(541613713, cnt);
        sums_from_hist(cnt, 12, S, ovf);
        int depth = 0;
        for (int k = 1; k <= 12; k++) { if (ovf[k] || !is_prime128(S[k])) break; depth = k; }
        expect(depth == 7, "p = 541613713 = a(7) has depth exactly 7 (S_8 = %s is composite)", u128_str(S[8], b1));
        printf("        S_1..S_8 =");
        for (int k = 1; k <= 8; k++) printf(" %s", u128_str(S[k], b1));
        printf("\n");

        /* the parallel totals pipeline must reproduce these sums */
        scan_t *tot = quiet_totals(541613713, 777777, threads);
        bool same = tot->primes_pub == direct_hist(541613713, cnt);
        for (int k = 1; k <= 12; k++) if (tot->Oovf[k] || tot->O[k] != S[k]) same = false;
        expect(same, "parallel totals pipeline (chunk 777777) reproduces S_1..S_12 at 541613713 and the prime count %" PRIu64, tot->primes_pub);
        scan_free(tot);
    }

    /* 4. the sequence itself */
    {
        scan_t *sc = quiet_scan(0, 600000000, 1000000, threads, 3, NULL, NULL, 0);
        bool ok = true;
        for (int n = 1; n <= NKNOWN; n++) {
            u64 a = a_from_best(sc->best, n);
            if (a != KNOWN[n]) { ok = false; printf("      a(%d) = %" PRIu64 ", OEIS has %" PRIu64 "\n", n, a, KNOWN[n]); }
        }
        expect(ok, "scan [0, 6e8] (chunk 1e6) reproduces a(1..7) = 5, 5, 13, 14593, 372313, 2315773, 541613713");
        printf("        %" PRIu64 " primes, %" PRIu64 " twins, %zu primes of depth >= 3, %.1f s\n", sc->primes, sc->twins, sc->nhits, sc->t1 - sc->t0);
        u64 a8 = a_from_best(sc->best, 8);
        printf("        a(8) %s 6e8\n", a8 ? "<=" : ">");
        if (a8) printf("        a(8) = %" PRIu64 "\n", a8);

        /* 5. chunking independence */
        scan_t *s2 = quiet_scan(0, 300000000, 300000000, 1, 3, NULL, NULL, 0);
        scan_t *s3 = quiet_scan(0, 300000000, 777777, threads, 3, NULL, NULL, 0);
        scan_t *s4 = quiet_scan(0, 300000000, 1000000, threads, 3, NULL, NULL, 0);
        bool same = same_hits(s2, s3) && same_hits(s2, s4) && s2->primes == s3->primes && s3->primes == s4->primes;
        for (int k = 1; k <= 12; k++) if (s2->O[k] != s3->O[k] || s3->O[k] != s4->O[k]) same = false;
        expect(same, "scan [0, 3e8] gives identical hits and sums with one chunk, chunk 777777 and chunk 1e6 (%zu hits of depth >= 3)", s2->nhits);

        /* 6. a scan that starts in the middle (offsets from a totals pass) */
        scan_t *tot = quiet_totals(299999999, 0, threads);
        scan_t *s5 = quiet_scan(300000000, 600000000, 1000000, threads, 3, tot->O, tot->Oovf, tot->primes_pub);
        size_t nabove = 0;
        bool sub = true;
        qsort(sc->hits, sc->nhits, sizeof(hit_t), cmp_hit);
        qsort(s5->hits, s5->nhits, sizeof(hit_t), cmp_hit);
        for (size_t i = 0; i < sc->nhits; i++) {
            if (sc->hits[i].p < 300000000) continue;
            if (nabove >= s5->nhits || s5->hits[nabove].p != sc->hits[i].p || s5->hits[nabove].depth != sc->hits[i].depth) sub = false;
            nabove++;
        }
        if (nabove != s5->nhits) sub = false;
        expect(sub && a_from_best(s5->best, 7) == KNOWN[7], "scan [3e8, 6e8] after a totals pass finds exactly the hits above 3e8 (%zu) including a(7)", s5->nhits);
        scan_free(sc); scan_free(s2); scan_free(s3); scan_free(s4); scan_free(tot); scan_free(s5);
    }

    /* 7. checkpoint round trip */
    {
        state_t st, st2;
        memset(&st, 0, sizeof st);
        st.kmax = 12; st.origin = 0; st.end = 12345678901ull; st.lo = 4000000000ull; st.primes = 189961812;
        for (int k = 1; k <= 12; k++) { st.O[k] = ((u128)rng() << 64) | rng(); st.ovf[k] = k == 12; }
        st.best[7] = 541613713; st.best[3] = 13;
        const char *path = "/tmp/a252768-selftest.state";
        write_state(path, &st);
        bool ok = read_state(path, &st2) && st2.kmax == 12 && st2.end == st.end && st2.lo == st.lo && st2.primes == st.primes &&
                  st2.best[7] == 541613713 && st2.best[3] == 13 && st2.best[5] == 0;
        for (int k = 1; k <= 12; k++) if (st2.O[k] != st.O[k] || st2.ovf[k] != st.ovf[k]) ok = false;
        unlink(path);
        expect(ok, "checkpoint file round trip");
    }

    (void)b2;
    printf("%s (%.1f s)\n", st_fail ? "SELFTEST FAILED" : "all tests passed", now() - t0);
    return st_fail ? 1 : 0;
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    setvbuf(stdout, NULL, _IOLBF, 0);
    stderr_tty = isatty(2);
    init_pow_tables();
    if (argc < 2) usage();
    const char *cmd = argv[1];
    if (!strcmp(cmd, "scan")) return cmd_scan(argc - 2, argv + 2);
    if (!strcmp(cmd, "verify") || !strcmp(cmd, "depth")) return cmd_verify(argc - 2, argv + 2);
    if (!strcmp(cmd, "selftest")) return cmd_selftest(argc - 2, argv + 2);
    usage();
    return 2;
}
