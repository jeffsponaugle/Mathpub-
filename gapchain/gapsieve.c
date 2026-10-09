/*
 * gapsieve.c -- constellation sieve for chains of consecutive primes with
 * gaps 2, 4, 6, ...  Same search and same CLI as gapchain_ps.c, but instead
 * of enumerating every prime it hunts the pattern directly:
 *
 *   A chain head p needs p + m(m+1) prime for m = 0..L-1.  Small primes
 *   forbid most residues of p outright (p = 2 mod 3 is forced, only 2 of 5
 *   residues survive mod 5, ...), so a wheel over the first several primes
 *   leaves only a fraction of integers as candidates.  Each wheel class is
 *   then sieved against every member for primes up to the pre-sieve limit,
 *   the rare survivors get batched base-2 strong tests, and a full hit is
 *   verified exactly: every member deterministically prime, no interloper
 *   prime between members, then extended to its true length.
 *
 * The sieve can only ever strike a candidate whose member p+o is a genuine
 * multiple of a sieve prime q < p+o, so no true chain is ever lost -- with
 * one exception: p+o equal to q itself.  That is only possible for p below
 * the sieve depth, and that strip is searched exactly by a direct scan over
 * a small local prime sieve (which also keeps tiny cases like 3,5 correct).
 *
 * The sieve runs as precomputed bit patterns, not bit-by-bit strikes: which
 * candidates a prime q kills is the same set in every class up to a
 * rotation, so each q (and, by CRT, each group of q's) is one pattern read
 * at a per-class offset and ANDed a whole word at a time.  That makes
 * sieving against all L members nearly free, so every small prime removes
 * up to L/q of candidates instead of the J/q a strike loop could afford.
 * Primes between --presieve and --depth, if any, are still struck bit by
 * bit against the first J = min(L,6) members; at the defaults there are
 * none, because once the patterns have run, striking costs more than
 * screening the few survivors it would remove.
 *
 * --gaps dec searches the mirror image, gaps 2(L-1), ..., 4, 2 (OEIS
 * A263049), and --gaps both does both.  A decreasing chain from p has members
 * p + m(2L-1-m), the increasing pattern reflected, so it rules out as many
 * residues at every wheel prime.  It shares the wheel, the segment grid and
 * the pattern groups, with its own residues, patterns and final check, which
 * extends backwards: the pattern ends in a gap of 2, so a longer decreasing
 * chain can only grow at its front.  Each shape costs the same, so both
 * together take as long as two runs.
 *
 * Checkpointing, --first ordering, progress and stats work as in
 * gapchain_ps.c: segments are handed out in ascending order, a checkpoint is
 * the count of leading finished segments, and --first lets in-flight lower
 * segments finish before declaring a winner.  One difference: gapsieve
 * finishes the segment holding a --first hit, so the checkpoint must also
 * record the hit, or a resume would skip past the answer.
 *
 * Positions are 128-bit, so the search runs past 2^64.  A segment whose
 * numbers all stay below 2^64 is screened with one-word Montgomery arithmetic;
 * any other uses two words.  The ceiling is PSI13, about 3.317e24, below
 * which Miller-Rabin with the 13 smallest prime bases is a proven primality
 * test; every reported chain is certified with it.
 *
 * Build (no dependencies):
 *   cc -O3 -pthread gapsieve.c -o gapsieve
 */

#define _GNU_SOURCE
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <pthread.h>
#include <signal.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <time.h>
#include <unistd.h>

#define KBLOCKS  (1u << 17)          /* candidates per class per segment */
#define NBUCKET  9                   /* stats buckets: len l .. l+8 and up */
#define A_CAP    24000               /* max admissible residues in wheel */
#define W_CAP    30000000ULL         /* max wheel modulus */
#define EXT_CAP  80                  /* extension steps past l before "+" */
#define CKPT_MAGIC "gapsieve-checkpoint 1"
#define NO_SEG UINT64_MAX

typedef unsigned __int128 u128;
#define U128_MAX (~(u128)0)

/*
 * PSI13 = 3317044064679887385961981: Miller-Rabin with the 13 smallest prime
 * bases is correct for every n below it (Sorenson and Webster, 2015), so no
 * number the search checks may reach it.
 */
#define PSI13 ((u128)179817 << 64 | 0x51adc5b22410a5fdULL)

static u128 g_lo, g_hi;              /* chain heads searched */
static uint64_t g_segsize, g_nseg, g_next_seg;
static u128 g_sieve_lo;              /* first value covered by the sieve */
static int g_k, g_J;                 /* chain length; sieved offsets */
static int g_progress, g_quiet, g_first, g_gpu;
static uint32_t g_depth = 199;

/* wheel */
static uint64_t g_W;                 /* wheel modulus */

/* sieve primes above the wheel, up to g_depth */
static uint32_t *g_q, *g_qinv;       /* q and W^-1 mod q */
static uint32_t g_nq;

/*
 * The two chain shapes.  INC has gaps 2, 4, ..., 2(L-1) (OEIS A016045) and
 * members p + m(m+1); DEC has gaps 2(L-1), ..., 4, 2 (OEIS A263049) and
 * members p + m(2L-1-m).  DEC is INC's mirror image, so it has exactly as
 * many wheel classes and sieve survivors and costs the same to search.
 * found/hist/best are guarded by g_lock; workers read best_seg lock-free.
 */
enum { INC, DEC };
struct shape {
    int on;                          /* searched in this run */
    const char *tag, *oeis;
    uint32_t offs[1000];             /* member offsets from the head */
    uint32_t *res, A;                /* admissible residues mod W */
    uint32_t *dtab;                  /* strike table, [i*J + oi] */
    uint64_t found, hist[NBUCKET];   /* counting: finished prefix only */
    uint64_t best_seg;               /* --first: best candidate so far */
    u128 best_p;
    int best_len, best_plus;
};
static struct shape g_sh[2] = {
    { .tag = "inc", .oeis = "A016045", .best_seg = NO_SEG },
    { .tag = "dec", .oeis = "A263049", .best_seg = NO_SEG },
};

/*
 * Pre-sieve: the first g_nps sieve primes are applied as precomputed bit
 * patterns against all L offsets, several primes per pattern (CRT), with
 * whole-word ANDs.  The remaining primes up to g_depth are struck bit by
 * bit against the first J offsets only.
 */
struct psgroup {
    uint64_t *pat[2];                /* per shape: 1 = allowed; period P,
                                        then +K bits */
    uint32_t P, first, n;            /* period; primes g_q[first..first+n) */
};
static struct psgroup *g_groups;
static uint32_t g_ngroups, g_nps;
static uint64_t *g_crt;              /* CRT coefficient of g_q[i] in its group */
static uint32_t g_presieve = 199;
static uint32_t g_pscap = 1u << 24;  /* max pattern period, bits */

/* frontier bookkeeping, identical to gapchain_ps */
static pthread_mutex_t g_lock = PTHREAD_MUTEX_INITIALIZER;
static uint64_t g_frontier;
struct pend { uint64_t seg, hist[2][NBUCKET]; };
static struct pend *g_pend;
static size_t g_npend, g_cpend;

static pthread_mutex_t g_out = PTHREAD_MUTEX_INITIALIZER;
static volatile sig_atomic_t g_stop;
static int g_active;

static const char *g_ckfile;
static char *g_cktmp;
static unsigned g_ckint = 60;
static double g_t0;
static uint64_t g_base;

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static FILE *g_log;                  /* --log: timestamped copy of messages */

/* one screen line; a longer status line would wrap and stop overwriting */
static int term_width(void)
{
    struct winsize ws;
    if (ioctl(STDERR_FILENO, TIOCGWINSZ, &ws) == 0 && ws.ws_col > 1)
        return ws.ws_col;
    return 80;
}

static void log_time(char *ts, size_t n)
{
    time_t now = time(NULL);
    struct tm tm;
    localtime_r(&now, &tm);
    strftime(ts, n, "%Y-%m-%d %H:%M:%S", &tm);
}

/*
 * Append text to the log, one timestamped line per non-empty line, so
 * messages that open with a blank line (the summaries) still log cleanly.
 * Caller holds g_out.
 */
static void log_locked(const char *text)
{
    if (!g_log)
        return;
    char ts[32];
    log_time(ts, sizeof ts);
    for (const char *s = text; *s; ) {
        size_t n = strcspn(s, "\n");
        if (n)
            fprintf(g_log, "%s  %.*s\n", ts, (int)n, s);
        s += n + (s[n] == '\n');
    }
    fflush(g_log);
}

/* a message for the terminal and the log, clearing any status line first */
static void note(const char *fmt, ...)
{
    char buf[4096];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    pthread_mutex_lock(&g_out);
    if (g_progress)
        fprintf(stderr, "\r%*s\r", term_width() - 1, "");
    fputs(buf, stderr);
    fflush(stderr);
    log_locked(buf);
    pthread_mutex_unlock(&g_out);
}

/* the log only */
static void log_note(const char *fmt, ...)
{
    char buf[4096];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    pthread_mutex_lock(&g_out);
    log_locked(buf);
    pthread_mutex_unlock(&g_out);
}

static void fmt_dur(char *buf, size_t n, double s)
{
    if (!(s >= 0.0) || s > 3.0e8) {
        snprintf(buf, n, "--:--:--");
        return;
    }
    uint64_t t = (uint64_t)s, d = t / 86400;
    t %= 86400;
    if (d)
        snprintf(buf, n, "%" PRIu64 "d%02" PRIu64 ":%02" PRIu64 ":%02" PRIu64,
                 d, t / 3600, (t / 60) % 60, t % 60);
    else
        snprintf(buf, n, "%02" PRIu64 ":%02" PRIu64 ":%02" PRIu64,
                 t / 3600, (t / 60) % 60, t % 60);
}

static void fmt_rate(char *buf, size_t n, double v)
{
    static const char *unit[] = { "", "K", "M", "G", "T", "P" };
    int i = 0;
    while (v >= 1000.0 && i < 5) {
        v /= 1000.0;
        i++;
    }
    snprintf(buf, n, "%.2f%s/s", v, unit[i]);
}

/*
 * Parse N exactly.  Plain digits, e-notation (2.5e15) and magnitude suffixes
 * (300T) all go through integer arithmetic, since floating point cannot hold
 * numbers this large exactly; a fractional value rounds half up.  Separators
 * , _ and space are ignored.  "max" is the largest value, which main lowers
 * to what the search can reach.
 */
static int parse_num(const char *arg, u128 *out)
{
    char buf[128];
    size_t n = 0;
    for (const char *s = arg; *s; s++) {
        if (*s == ',' || *s == '_' || *s == ' ')
            continue;
        if (n + 1 >= sizeof buf)
            return -1;
        buf[n++] = *s;
    }
    buf[n] = '\0';
    if (!n)
        return -1;
    if (!strcmp(buf, "max")) {
        *out = U128_MAX;
        return 0;
    }

    /* the mantissa's digits, with its decimal point moved into exp */
    u128 m = 0;
    int exp = 0, digits = 0, point = 0;
    const char *s = buf;
    for (; *s; s++) {
        if (*s == '.' && !point) {
            point = 1;
            continue;
        }
        if (!isdigit((unsigned char)*s))
            break;
        unsigned dg = (unsigned)(*s - '0');
        if (m > (U128_MAX - dg) / 10)
            return -1;
        m = m * 10 + dg;
        digits++;
        exp -= point;
    }
    if (!digits)
        return -1;
    if (*s == 'e' || *s == 'E') {
        int neg = 0, e = 0, ed = 0;
        s++;
        if (*s == '+' || *s == '-')
            neg = *s++ == '-';
        for (; isdigit((unsigned char)*s); s++, ed++)
            if ((e = e * 10 + (*s - '0')) > 1000)
                return -1;
        if (!ed)
            return -1;
        exp += neg ? -e : e;
    }
    if (*s) {
        switch (toupper((unsigned char)*s)) {
        case 'K': exp += 3;  break;
        case 'M': exp += 6;  break;
        case 'B':
        case 'G': exp += 9;  break;
        case 'T': exp += 12; break;
        case 'P': exp += 15; break;
        default:  return -1;
        }
        if (*++s)
            return -1;
    }
    for (; exp > 0; exp--) {
        if (m > U128_MAX / 10)
            return -1;
        m *= 10;
    }
    /* the last digit dropped is the most significant one: round on it */
    int up = 0;
    for (; exp < 0; exp++) {
        up = m % 10 >= 5;
        m /= 10;
    }
    *out = m + (u128)up;
    return 0;
}

static int parse_u64(const char *arg, uint64_t *out)
{
    u128 v;
    if (parse_num(arg, &v) || v > UINT64_MAX)
        return -1;
    *out = (uint64_t)v;
    return 0;
}

/* decimal text of v into buf, which needs 40 bytes */
static char *u128_str(char *buf, u128 v)
{
    const uint64_t e19 = 10000000000000000000ULL;
    if (v <= UINT64_MAX) {
        snprintf(buf, 40, "%" PRIu64, (uint64_t)v);
        return buf;
    }
    u128 hi = v / e19;                          /* below 3.5e19 */
    uint64_t lo = (uint64_t)(v % e19);
    if (hi <= UINT64_MAX) {
        snprintf(buf, 40, "%" PRIu64 "%019" PRIu64, (uint64_t)hi, lo);
    } else {                                    /* 39 digits, the first 1..3 */
        buf[0] = (char)('0' + (int)(hi / e19));
        snprintf(buf + 1, 39, "%019" PRIu64 "%019" PRIu64,
                 (uint64_t)(hi % e19), lo);
    }
    return buf;
}

/* u128_str into a temporary that lives until the end of the enclosing block */
#define U128S(v) u128_str((char[40]){ 0 }, (v))

/* ---------- primality: Montgomery arithmetic, deterministic MR ---------- */

struct mont { uint64_t n, ninv, r2, one; };

static uint64_t mont_mul(const struct mont *m, uint64_t a, uint64_t b)
{
    /* t + lo*n can need 129 bits when n > 2^63, so sum the halves: the low
     * 64 bits cancel by construction, leaving only their carry. */
    __uint128_t t = (__uint128_t)a * b;
    uint64_t lo = (uint64_t)t * m->ninv;
    __uint128_t mn = (__uint128_t)lo * m->n;
    __uint128_t r = (t >> 64) + (mn >> 64) + ((uint64_t)t != 0);
    if (r >= m->n)
        r -= m->n;
    return (uint64_t)r;
}

static void mont_init(struct mont *m, uint64_t n)
{
    m->n = n;
    uint64_t inv = n;                 /* Newton: inverse of n mod 2^64 */
    for (int i = 0; i < 5; i++)
        inv *= 2 - n * inv;
    m->ninv = (uint64_t)0 - inv;      /* -n^-1 mod 2^64 */
    m->one = (0 - n) % n;             /* unsigned: (2^64 - n) mod n = 2^64 mod n */
    __uint128_t r2 = ((__uint128_t)m->one << 64) % n;
    m->r2 = (uint64_t)r2;             /* 2^128 mod n */
}

static int sprp(const struct mont *m, uint64_t a)
{
    uint64_t n = m->n;
    if (a >= n) {
        a %= n;
        if (!a)
            return 1;
    }
    uint64_t d = n - 1;
    int s = __builtin_ctzll(d);
    d >>= s;

    uint64_t x = mont_mul(m, a, m->r2);          /* a in Montgomery form */
    uint64_t r = m->one, nm1 = mont_mul(m, n - 1, m->r2);
    for (uint64_t e = d; e; e >>= 1) {
        if (e & 1)
            r = mont_mul(m, r, x);
        x = mont_mul(m, x, x);
    }
    if (r == m->one || r == nm1)
        return 1;
    while (--s > 0) {
        r = mont_mul(m, r, r);
        if (r == nm1)
            return 1;
    }
    return 0;
}

static const uint32_t small_primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23,
                                         29, 31, 37 };

/* deterministic for all n < 2^64 (12-base bound is 3.18e23) */
static int is_prime64(uint64_t n)
{
    if (n < 2)
        return 0;
    for (int i = 0; i < 12; i++) {
        if (n % small_primes[i] == 0)
            return n == small_primes[i];
    }
    if (n < 41 * 41)
        return 1;
    struct mont m;
    mont_init(&m, n);
    for (int i = 0; i < 12; i++)
        if (!sprp(&m, small_primes[i]))
            return 0;
    return 1;
}

/*
 * Base-2 strong probable-prime test for odd n > 2, specialised: left-to-right
 * exponentiation where "multiply by the base" is a modular doubling, and no
 * R^2 (so no 128-bit division) because the base never enters Montgomery form
 * through a multiply.  A "composite" verdict is a proof.
 */
static inline uint64_t mmul(uint64_t a, uint64_t b, uint64_t n, uint64_t ninv)
{
    __uint128_t t = (__uint128_t)a * b;
    uint64_t lo = (uint64_t)t * ninv;
    __uint128_t mn = (__uint128_t)lo * n;
    __uint128_t r = (t >> 64) + (mn >> 64) + ((uint64_t)t != 0);
    if (r >= n)
        r -= n;
    return (uint64_t)r;
}

static int sprp2(uint64_t n)
{
    uint64_t inv = n;
    for (int i = 0; i < 5; i++)
        inv *= 2 - n * inv;
    uint64_t ninv = (uint64_t)0 - inv;
    uint64_t one = ((uint64_t)0 - n) % n;        /* 2^64 mod n */
    uint64_t nm1 = n - one;                      /* Montgomery form of -1 */
    uint64_t d = n - 1;
    int s = __builtin_ctzll(d);
    d >>= s;

    int b = 63 - __builtin_clzll(d);
    uint64_t r = one >= n - one ? one - (n - one) : one + one;   /* 2 */
    while (--b >= 0) {
        r = mmul(r, r, n, ninv);
        if ((d >> b) & 1)
            r = r >= n - r ? r - (n - r) : r + r;
    }
    if (r == one || r == nm1)
        return 1;
    while (--s > 0) {
        r = mmul(r, r, n, ninv);
        if (r == nm1)
            return 1;
    }
    return 0;
}

#ifndef LANES
#define LANES 4
#endif

/*
 * sprp2 on LANES independent odd n > 2 in lockstep.  One test is a chain of
 * ~60 dependent Montgomery squarings, so a single test leaves the multipliers
 * mostly idle; interleaved chains keep them busy (4 measured best).  Selects replace the
 * per-bit branch, which is unpredictable.  ok[l] = 1 iff n[l] is a base-2
 * strong probable prime -- exactly sprp2(n[l]).
 */
static void sprp2_xN(const uint64_t *nv, int *ok)
{
    uint64_t n[LANES], ninv[LANES], one[LANES], nm1[LANES], d[LANES], r[LANES];
    int s[LANES], top = 0;

    for (int l = 0; l < LANES; l++) {
        uint64_t inv = n[l] = nv[l];
        for (int i = 0; i < 5; i++)
            inv *= 2 - n[l] * inv;
        ninv[l] = (uint64_t)0 - inv;
        one[l] = ((uint64_t)0 - n[l]) % n[l];
        nm1[l] = n[l] - one[l];
        d[l] = n[l] - 1;
        s[l] = __builtin_ctzll(d[l]);
        d[l] >>= s[l];
        int b = 63 - __builtin_clzll(d[l]);
        top = b > top ? b : top;
        r[l] = one[l];                /* squaring 1 is harmless above top */
    }
    for (int b = top; b >= 0; b--) {
        for (int l = 0; l < LANES; l++)
            r[l] = mmul(r[l], r[l], n[l], ninv[l]);
        for (int l = 0; l < LANES; l++) {
            uint64_t y = n[l] - r[l];
            uint64_t dbl = r[l] >= y ? r[l] - y : r[l] + r[l];
            uint64_t m = (uint64_t)0 - ((d[l] >> b) & 1);
            r[l] = (dbl & m) | (r[l] & ~m);
        }
    }
    for (int l = 0; l < LANES; l++) {
        uint64_t x = r[l];
        int good = x == one[l] || x == nm1[l];
        for (int i = 1; i < s[l] && !good; i++) {
            x = mmul(x, x, n[l], ninv[l]);
            good = x == nm1[l];
        }
        ok[l] = good;
    }
}

/* ---------- two-word arithmetic, for numbers from 2^64 up to PSI13 ---------- */

static inline int ctz128(u128 x)                /* x != 0 */
{
    uint64_t lo = (uint64_t)x;
    return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((uint64_t)(x >> 64));
}

static inline int top128(u128 x)                /* top set bit; x != 0 */
{
    uint64_t hi = (uint64_t)(x >> 64);
    return hi ? 127 - __builtin_clzll(hi) : 63 - __builtin_clzll((uint64_t)x);
}

/* x mod q for q < 2^32, without a 128-bit division */
static inline uint32_t mod_small(u128 x, uint32_t q)
{
    uint64_t r64 = ((uint64_t)0 - q) % q;       /* 2^64 mod q */
    return (uint32_t)(((uint64_t)(x >> 64) % q * r64 + (uint64_t)x % q) % q);
}

/* -n^-1 mod 2^64 for odd n, from its low word */
static inline uint64_t neg_inv64(u128 n)
{
    uint64_t n0 = (uint64_t)n, inv = n0;         /* Newton, as in mont_init */
    for (int i = 0; i < 5; i++)
        inv *= 2 - n0 * inv;
    return (uint64_t)0 - inv;
}

/*
 * Montgomery product a*b/2^128 mod n, for odd n < 2^127 and a, b < n: CIOS
 * over two 64-bit words, ninv = -n^-1 mod 2^64.  Each pass adds one word of
 * a*b, then a multiple of n that clears the low word, and drops that word.
 * Every sum fits in 128 bits, and the result is below n.
 */
static inline u128 mmul2(u128 a, u128 b, u128 n, uint64_t ninv)
{
    uint64_t a0 = (uint64_t)a, a1 = (uint64_t)(a >> 64);
    uint64_t b0 = (uint64_t)b, b1 = (uint64_t)(b >> 64);
    uint64_t n0 = (uint64_t)n, n1 = (uint64_t)(n >> 64);
    uint64_t t0, t1, t2, m;
    u128 c;

    c = (u128)a0 * b0;                          /* t = a*b0 */
    t0 = (uint64_t)c;
    c = (u128)a1 * b0 + (c >> 64);
    t1 = (uint64_t)c;
    t2 = (uint64_t)(c >> 64);
    m = t0 * ninv;                              /* t = (t + m*n) / 2^64 */
    c = ((u128)m * n0 + t0) >> 64;
    c += (u128)m * n1 + t1;
    t0 = (uint64_t)c;
    c = (c >> 64) + t2;
    t1 = (uint64_t)c;
    t2 = (uint64_t)(c >> 64);

    c = (u128)a0 * b1 + t0;                     /* t += a*b1 */
    t0 = (uint64_t)c;
    c = (u128)a1 * b1 + t1 + (c >> 64);
    t1 = (uint64_t)c;
    c = (c >> 64) + t2;
    t2 = (uint64_t)c;
    m = t0 * ninv;                              /* t = (t + m*n) / 2^64 */
    c = ((u128)m * n0 + t0) >> 64;
    c += (u128)m * n1 + t1;
    t0 = (uint64_t)c;
    c = (c >> 64) + t2;
    t1 = (uint64_t)c;
    t2 = (uint64_t)(c >> 64);

    u128 t = (u128)t1 << 64 | t0;               /* t < 2n */
    return t2 || t >= n ? t - n : t;
}

/*
 * 2^128 mod n, the Montgomery form of 1, for odd n < 2^82.  A 128-bit
 * remainder is a slow library loop, so above 2^64 the quotient comes from
 * floating point: 2^128/n in double precision is within 2^13 of the truth,
 * the leftover multiple of n is small enough to estimate the same way, and
 * a final step either way makes the remainder exact.
 */
static inline u128 mont_one2(u128 n)
{
    uint64_t n1 = (uint64_t)(n >> 64);
    if (!n1)                                    /* only in the segment at 2^64 */
        return ((u128)0 - n) % n;
    double dn = (double)n1 * 0x1p64 + (double)(uint64_t)n;
    double qd = 0x1p128 / dn;
    uint64_t q = qd < 0x1p64 ? (uint64_t)qd : UINT64_MAX;
    /* 2^128 - q*n, exact as a signed number: |it| < 2^13 n < 2^95 */
    __int128 r = (__int128)((u128)0 - (u128)q * n), sn = (__int128)n;
    double rd = (double)(int64_t)(r >> 64) * 0x1p64 + (double)(uint64_t)r;
    r -= (__int128)(int64_t)(rd / dn) * sn;
    while (r < 0)
        r += sn;
    while (r >= sn)
        r -= sn;
    return (u128)r;
}

/* a + b mod n, for a, b < n < 2^127 */
static inline u128 addmod2(u128 a, u128 b, u128 n)
{
    u128 y = n - b;
    return a >= y ? a - y : a + b;
}

/* sprp2 for odd n > 2 below 2^127, in two-word Montgomery form */
static int sprp2w(u128 n)
{
    uint64_t ninv = neg_inv64(n);
    u128 one = mont_one2(n);                    /* 2^128 mod n */
    u128 nm1 = n - one;
    u128 d = n - 1;
    int s = ctz128(d);
    d >>= s;

    int b = top128(d);
    u128 r = addmod2(one, one, n);              /* 2 */
    while (--b >= 0) {
        r = mmul2(r, r, n, ninv);
        if ((d >> b) & 1)
            r = addmod2(r, r, n);
    }
    if (r == one || r == nm1)
        return 1;
    while (--s > 0) {
        r = mmul2(r, r, n, ninv);
        if (r == nm1)
            return 1;
    }
    return 0;
}

/* sprp2w on LANES numbers in lockstep, like sprp2_xN */
static void sprp2w_xN(const u128 *nv, int *ok)
{
    u128 n[LANES], one[LANES], nm1[LANES], d[LANES], r[LANES];
    uint64_t ninv[LANES];
    int s[LANES], top = 0;

    for (int l = 0; l < LANES; l++) {
        n[l] = nv[l];
        ninv[l] = neg_inv64(n[l]);
        one[l] = mont_one2(n[l]);
        nm1[l] = n[l] - one[l];
        d[l] = n[l] - 1;
        s[l] = ctz128(d[l]);
        d[l] >>= s[l];
        int b = top128(d[l]);
        top = b > top ? b : top;
        r[l] = one[l];                /* squaring 1 is harmless above top */
    }
    for (int b = top; b >= 0; b--) {
        for (int l = 0; l < LANES; l++)
            r[l] = mmul2(r[l], r[l], n[l], ninv[l]);
        for (int l = 0; l < LANES; l++) {
            u128 dbl = addmod2(r[l], r[l], n[l]);
            u128 m = (u128)0 - (u128)((d[l] >> b) & 1);
            r[l] = (dbl & m) | (r[l] & ~m);
        }
    }
    for (int l = 0; l < LANES; l++) {
        u128 x = r[l];
        int good = x == one[l] || x == nm1[l];
        for (int i = 1; i < s[l] && !good; i++) {
            x = mmul2(x, x, n[l], ninv[l]);
            good = x == nm1[l];
        }
        ok[l] = good;
    }
}

/* strong probable-prime test to base a < n, n odd, two-word Montgomery */
static int sprpw(u128 n, uint64_t ninv, u128 one, uint32_t a)
{
    u128 nm1 = n - one, d = n - 1;
    int s = ctz128(d);
    d >>= s;

    u128 x = 0;                                 /* a*2^128 mod n */
    for (int b = 31 - __builtin_clz(a); b >= 0; b--) {
        x = addmod2(x, x, n);
        if ((a >> b) & 1)
            x = addmod2(x, one, n);
    }
    u128 r = one;
    for (u128 e = d; e; e >>= 1) {
        if (e & 1)
            r = mmul2(r, x, n, ninv);
        x = mmul2(x, x, n, ninv);
    }
    if (r == one || r == nm1)
        return 1;
    while (--s > 0) {
        r = mmul2(r, r, n, ninv);
        if (r == nm1)
            return 1;
    }
    return 0;
}

/* deterministic for 41 < n < PSI13: the 13 smallest prime bases */
static int is_prime_wide(u128 n)
{
    static const uint32_t base[13] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                                       31, 37, 41 };
    for (int i = 0; i < 13; i++)
        if (mod_small(n, base[i]) == 0)
            return 0;
    uint64_t ninv = neg_inv64(n);
    u128 one = mont_one2(n);
    for (int i = 0; i < 13; i++)
        if (!sprpw(n, ninv, one, base[i]))
            return 0;
    return 1;
}

/* exact, below PSI13 */
static int is_prime(u128 n)
{
    return n <= UINT64_MAX ? is_prime64((uint64_t)n) : is_prime_wide(n);
}

/*
 * Keep the candidates c[0..cnt), offsets from base, whose member
 * base + c + off passes the base-2 strong test, in order; returns the new
 * count.  Pads the last batch by repeating an entry.  wide says some member
 * in the segment reaches 2^64, so the two-word test is needed.
 */
static uint32_t screen_member(uint64_t *c, uint32_t cnt, u128 base,
                              uint64_t off, int wide)
{
    uint32_t kept = 0;
    int ok[LANES];
    if (!wide) {
        uint64_t b = (uint64_t)base + off;
        for (uint32_t i = 0; i < cnt; i += LANES) {
            uint64_t nv[LANES];
            for (int l = 0; l < LANES; l++)
                nv[l] = b + c[i + l < cnt ? i + l : cnt - 1];
            sprp2_xN(nv, ok);
            for (int l = 0; l < LANES && i + l < cnt; l++)
                if (ok[l])
                    c[kept++] = c[i + l];
        }
        return kept;
    }
    u128 b = base + off;
    for (uint32_t i = 0; i < cnt; i += LANES) {
        u128 nv[LANES];
        for (int l = 0; l < LANES; l++)
            nv[l] = b + c[i + l < cnt ? i + l : cnt - 1];
        sprp2w_xN(nv, ok);
        for (int l = 0; l < LANES && i + l < cnt; l++)
            if (ok[l])
                c[kept++] = c[i + l];
    }
    return kept;
}

/* cheap screen for n > 37: tiny trial division, then one base-2 strong test */
static int probable_prime(u128 n)
{
    if (n <= UINT64_MAX) {
        uint64_t v = (uint64_t)n;
        for (int i = 1; i < 12; i++)
            if (v % small_primes[i] == 0)
                return 0;
        return sprp2(v);
    }
    for (int i = 1; i < 12; i++)
        if (mod_small(n, small_primes[i]) == 0)
            return 0;
    return sprp2w(n);
}

/* ---------- exact chain verification ---------- */

/*
 * INC: p has passed screening for the first g_k members.  Verify exactly and
 * return the true chain length (extending past g_k), or 0 if p is not a
 * chain head.  Sets *plus if the chain is still alive at the extension cap.
 *
 * The walk uses cheap tests whose composite verdicts are exact (trial
 * division, base-2 strong test); only a number that survives both gets the
 * full deterministic test.  Members are provisionally prime during the walk
 * and certified at the end, truncating at a base-2 pseudoprime if one ever
 * slips through -- so the reported chain is deterministically correct.
 */
static int chain_verify_inc(u128 p, int *plus)
{
    *plus = 0;
    int len = 1;
    for (int m = 1; m <= g_k + EXT_CAP; m++) {
        u128 prev = p + (uint64_t)(m - 1) * m;
        u128 next = p + (uint64_t)m * (m + 1);
        /* no interloper prime strictly between the members */
        int clean = 1;
        for (u128 x = prev + 2; x < next && clean; x += 2)
            if (probable_prime(x) && is_prime(x))
                clean = 0;
        if (!clean || !probable_prime(next))
            break;
        len++;
        if (m == g_k + EXT_CAP)
            *plus = 1;
    }
    int true_len = 0;
    while (true_len < len &&
           is_prime(p + (uint64_t)true_len * (true_len + 1)))
        true_len++;
    if (true_len < len)
        *plus = 0;
    return true_len >= g_k ? true_len : 0;
}

/*
 * DEC: p has passed screening as the head of the exact L-prime pattern
 * p, p+2(L-1), ..., p+(L-1)L.  Verify it exactly, then extend backwards: the
 * pattern ends in a gap of 2, so it can only be the tail of a longer
 * decreasing chain.  Returns the chain's true length, or 0 if p fails, and
 * the chain's first prime in *head.  Sets *plus if still alive at the cap.
 * Backward steps can reach small numbers, where only is_prime is exact.
 */
static int chain_verify_dec(u128 p, u128 *head, int *plus)
{
    const uint32_t *o = g_sh[DEC].offs;

    *plus = 0;
    *head = p;
    for (int m = 0; m + 1 < g_k; m++)
        for (u128 x = p + o[m] + 2; x < p + o[m + 1]; x += 2)
            if (probable_prime(x) && is_prime(x))
                return 0;
    for (int m = 0; m < g_k; m++)
        if (!is_prime(p + o[m]))
            return 0;

    u128 h = p;
    int len = g_k;
    while (len < g_k + EXT_CAP + 1) {
        uint64_t gap = 2 * (uint64_t)len;
        if (h < gap + 3)
            break;
        u128 prev = h - gap;
        int ok = is_prime(prev);
        for (u128 x = prev + 2; ok && x < h; x += 2)
            ok = !is_prime(x);
        if (!ok)
            break;
        h = prev;
        len++;
    }
    *plus = len == g_k + EXT_CAP + 1;
    *head = h;
    return len;
}

/* ---------- shared reporting ---------- */

/*
 * A hit prints as its whole chain, from its first prime p.  Its key is the
 * head of the L-prime pattern that was searched for: p itself for INC, but
 * a DEC pattern ends its chain, so a longer DEC chain starts before its key.
 * --first ranks by key, the head A016045 and A263049 record.
 */
struct hit { u128 p, key; int len, plus, dir; };

static int hit_cmp(const void *a, const void *b)
{
    const struct hit *x = a, *y = b;
    if (x->p != y->p)
        return x->p < y->p ? -1 : 1;
    return x->dir - y->dir;
}

static void add_hit(struct hit **hits, size_t *nh, size_t *ch, struct hit h)
{
    if (*nh == *ch) {
        *ch = *ch ? *ch * 2 : 64;
        *hits = realloc(*hits, *ch * sizeof **hits);
        if (!*hits)
            _exit(1);
    }
    (*hits)[(*nh)++] = h;
}

/* offset of member m in a chain of len primes of shape d */
static uint64_t member_off(int d, int len, int m)
{
    return d == INC ? (uint64_t)m * (m + 1)
                    : (uint64_t)m * (uint64_t)(2 * len - 1 - m);
}

/* lines carry a shape tag whenever DEC is searched; INC-only output is
   unchanged */
static void fprint_hit(FILE *f, const struct hit *h)
{
    char s[40];
    if (g_sh[DEC].on)
        fprintf(f, "%s ", g_sh[h->dir].tag);
    fprintf(f, "%2d: %s", h->len, u128_str(s, h->p));
    for (int m = 1; m < h->len; m++)
        fprintf(f, " %s", u128_str(s, h->p + member_off(h->dir, h->len, m)));
    if (h->plus)
        fputs(" +", f);
    fputc('\n', f);
}

static void print_hit(const struct hit *h)
{
    fprint_hit(stdout, h);
}

/* sort a segment's hits, then print and count them */
static void flush_hits(struct hit *hits, size_t n, uint64_t (*hist)[NBUCKET])
{
    if (!n)
        return;
    qsort(hits, n, sizeof *hits, hit_cmp);
    pthread_mutex_lock(&g_out);
    for (size_t i = 0; i < n; i++) {
        print_hit(&hits[i]);
        int b = hits[i].len - g_k;
        hist[hits[i].dir][b > NBUCKET - 1 ? NBUCKET - 1 : b]++;
    }
    pthread_mutex_unlock(&g_out);
}

static void seg_complete(uint64_t seg, uint64_t (*hist)[NBUCKET])
{
    pthread_mutex_lock(&g_lock);
    if (g_npend == g_cpend) {
        size_t c = g_cpend ? g_cpend * 2 : 32;
        struct pend *p = realloc(g_pend, c * sizeof *p);
        if (!p) {
            note("gapsieve: out of memory tracking segments\n");
            _exit(1);
        }
        g_pend = p;
        g_cpend = c;
    }
    g_pend[g_npend].seg = seg;
    memcpy(g_pend[g_npend].hist, hist, sizeof g_pend[g_npend].hist);
    g_npend++;

    for (int moved = 1; moved; ) {
        moved = 0;
        for (size_t i = 0; i < g_npend; i++) {
            if (g_pend[i].seg != g_frontier)
                continue;
            for (int d = 0; d < 2; d++)
                for (int b = 0; b < NBUCKET; b++) {
                    g_sh[d].hist[b] += g_pend[i].hist[d][b];
                    g_sh[d].found += g_pend[i].hist[d][b];
                }
            g_pend[i] = g_pend[--g_npend];
            g_frontier++;
            moved = 1;
            break;
        }
    }
    pthread_mutex_unlock(&g_lock);
}

/* key is the pattern head; segments partition keys, so lower seg wins */
static void record_candidate(int d, uint64_t seg, u128 key, int len,
                             int plus)
{
    struct shape *S = &g_sh[d];
    uint64_t frontier;
    int better;

    pthread_mutex_lock(&g_lock);
    better = seg < S->best_seg;
    frontier = g_frontier;
    if (better) {
        S->best_p = key;
        S->best_len = len;
        S->best_plus = plus;
        __atomic_store_n(&S->best_seg, seg, __ATOMIC_RELEASE);
    }
    pthread_mutex_unlock(&g_lock);

    if (better)
        note("%s%scandidate: %s (length %d, segment %" PRIu64 ")%s\n",
             g_sh[DEC].on ? S->tag : "", g_sh[DEC].on ? " " : "",
             U128S(key), len, seg,
             seg > frontier ? "; finishing earlier segments" : "");
}

/* ---------- direct scan for the strip below the sieve region ---------- */

/*
 * Search [lo, hi] for chain heads by brute-force sieve; used for the small
 * region where a chain member could coincide with a sieve prime.  hi stays
 * below g_sieve_lo, which is at most depth + one wheel turn, so this sieve
 * is tiny.  Consecutiveness is exact: we walk the real prime list.
 */
static void direct_scan(uint64_t lo, uint64_t hi, const int *need,
                        struct hit **hits, size_t *nh, size_t *ch)
{
    uint64_t slack = 2ULL * (uint64_t)(g_k - 1) * g_k +
                     (uint64_t)(g_k + EXT_CAP + 2) * (g_k + EXT_CAP + 2) + 4096;
    size_t end = (size_t)(hi + slack);
    unsigned char *comp = calloc(end + 1, 1);
    uint64_t *pr = NULL;
    size_t np = 0, cp = 0;
    if (!comp)
        _exit(1);
    for (size_t i = 2; i * i <= end; i++)
        if (!comp[i])
            for (size_t j = i * i; j <= end; j += i)
                comp[j] = 1;
    for (size_t i = 2; i <= end; i++) {
        if (comp[i])
            continue;
        if (np == cp) {
            cp = cp ? cp * 2 : 4096;
            pr = realloc(pr, cp * sizeof *pr);
            if (!pr)
                _exit(1);
        }
        pr[np++] = i;
    }
    for (size_t i = 0; i < np && pr[i] <= hi; i++) {
        uint64_t p = pr[i];
        if (p < lo)
            continue;
        if (need[INC]) {
            int len = 1, plus = 0;
            for (size_t j = i + 1; j < np; j++) {
                if (pr[j] - pr[j - 1] != (uint64_t)(2 * (int)(j - i)))
                    break;
                len++;
                if (len == g_k + EXT_CAP + 1) {
                    plus = 1;
                    break;
                }
            }
            if (len >= g_k)
                add_hit(hits, nh, ch, (struct hit){ p, p, len, plus, INC });
        }
        if (need[DEC]) {
            /* the exact pattern from p, then back while the gaps keep
               growing by two */
            int ok = i + (size_t)g_k - 1 < np;
            for (int j = 1; ok && j < g_k; j++)
                ok = pr[i + j] - pr[i + j - 1] == 2 * (uint64_t)(g_k - j);
            if (ok) {
                size_t h = i;
                int len = g_k;
                while (h > 0 && len < g_k + EXT_CAP + 1 &&
                       pr[h] - pr[h - 1] == 2 * (uint64_t)len) {
                    h--;
                    len++;
                }
                add_hit(hits, nh, ch, (struct hit){ pr[h], p, len,
                        len == g_k + EXT_CAP + 1, DEC });
            }
        }
    }
    free(pr);
    free(comp);
}

/*
 * Build the candidate bitmap of one class of shape d: bit k set iff
 * candidate p = base + s0 + k*W survives the sieve.  bq[i] = base mod g_q[i];
 * psrc/psh are per-thread scratch of g_ngroups entries.
 */
static void class_bitmap(int d, uint64_t *bits, const uint64_t *bq,
                         uint64_t s0, const uint64_t **psrc, unsigned *psh)
{
    const uint32_t K = KBLOCKS, nq = g_nq;
    const int J = g_J;

        /*
         * Pre-sieve.  Bit k of this class corresponds to pattern bit
         * k + R, where R = CRT of r_q = W^-1 (base+s0) mod q; see
         * presieve_init for the derivation.
         */
        for (uint32_t gi = 0; gi < g_ngroups; gi++) {
            const struct psgroup *G = &g_groups[gi];
            uint64_t R = 0;
            for (uint32_t i = G->first; i < G->first + G->n; i++) {
                uint32_t q = g_q[i];
                uint64_t t = bq[i] + s0 % q;
                if (t >= q)
                    t -= q;
                R += (uint64_t)g_qinv[i] * t % q * g_crt[i];
            }
            R %= G->P;
            psrc[gi] = G->pat[d] + (R >> 6);
            psh[gi] = (unsigned)(R & 63);
        }
        /*
         * AND every pattern into an 8-word accumulator that stays in
         * registers, so the bitmap is written once rather than
         * read-modified-written once per pattern.
         * (x << 1) << (63 - sh) is x << (64 - sh), and 0 at sh = 0.
         */
        for (uint32_t w = 0; w < K / 64; w += 8) {
            uint64_t acc[8];
            for (int j = 0; j < 8; j++)
                acc[j] = ~0ULL;
            for (uint32_t gi = 0; gi < g_ngroups; gi++) {
                const uint64_t *src = psrc[gi] + w;
                unsigned sh = psh[gi];
                for (int j = 0; j < 8; j++)
                    acc[j] &= (src[j] >> sh) |
                              ((src[j + 1] << 1) << (63 - sh));
            }
            for (int j = 0; j < 8; j++)
                bits[w + j] = acc[j];
        }

        for (uint32_t i = g_nps; i < nq; i++) {
            uint32_t q = g_q[i];
            uint64_t t = bq[i] + s0 % q;
            if (t >= q)
                t -= q;
            /* k0 for offset 0: qinv * (q - t) mod q */
            uint32_t c0 = (uint32_t)
                (((uint64_t)g_qinv[i] * (t ? q - t : 0)) % q);
            const uint32_t *dt = &g_sh[d].dtab[(size_t)i * J];
            for (int oi = 0; oi < J; oi++) {
                uint32_t k0 = c0 >= dt[oi] ? c0 - dt[oi]
                                           : c0 + q - dt[oi];
                for (uint32_t k = k0; k < K; k += q)
                    bits[k >> 6] &= ~(1ULL << (k & 63));
            }
        }
}

/* ---------- the sieve worker ---------- */

/*
 * Claim the next segment and decide which shapes it needs.  --first: a shape
 * only needs the segments below its best candidate.  Within one run, once no
 * shape needs this claim none needs any later one either, since segments go
 * out in ascending order; after a resume that leaves exactly the unfinished
 * segments below the recorded hits.  Returns 0 when there is no more work.
 * The GPU claims up to max consecutive segments at a time (*n of them), which
 * share the first one's need[].
 */
static int claim_segments(uint64_t *seg, uint64_t *n, uint64_t max,
                          int need[2])
{
    if (g_stop)
        return 0;
    uint64_t s = __atomic_fetch_add(&g_next_seg, max, __ATOMIC_RELAXED);
    if (s >= g_nseg)
        return 0;
    int any = 0;
    for (int d = 0; d < 2; d++) {
        need[d] = g_sh[d].on &&
                  (!g_first ||
                   s < __atomic_load_n(&g_sh[d].best_seg, __ATOMIC_ACQUIRE));
        any |= need[d];
    }
    *seg = s;
    *n = g_nseg - s < max ? g_nseg - s : max;
    return any;
}

static int claim_segment(uint64_t *seg, int need[2])
{
    uint64_t n;
    return claim_segments(seg, &n, 1, need);
}

/* per-thread scratch for searching segments on the CPU */
struct scratch {
    uint64_t *bits, *bq, *cand;
    const uint64_t **psrc;
    unsigned *psh;
};

static int scratch_init(struct scratch *x)
{
    x->bits = malloc(KBLOCKS / 8);
    x->bq = malloc((size_t)g_nq * sizeof *x->bq);
    x->psrc = malloc((g_ngroups + 1) * sizeof *x->psrc);
    x->psh = malloc((g_ngroups + 1) * sizeof *x->psh);
    x->cand = malloc((size_t)KBLOCKS * sizeof *x->cand);
    return x->bits && x->bq && x->psrc && x->psh && x->cand ? 0 : -1;
}

static void scratch_free(struct scratch *x)
{
    free(x->bits);
    free(x->bq);
    free(x->psrc);
    free(x->psh);
    free(x->cand);
}

/* segment 0 owns the direct strip below the sieve region, which is then
   tiny: g_sieve_lo is depth + 512 */
static void segment_strip(uint64_t s, const int *need, struct hit **hits,
                          size_t *nh, size_t *ch)
{
    if (s == 0 && g_lo < g_sieve_lo) {
        u128 dhi = g_hi < g_sieve_lo - 1 ? g_hi : g_sieve_lo - 1;
        direct_scan((uint64_t)g_lo, (uint64_t)dhi, need, hits, nh, ch);
    }
}

/*
 * Search segment s on the CPU for the shapes in need[], appending its hits.
 * Only a stop mid-segment leaves it unfinished, and then this returns 0.  A
 * shape dropped because a lower --first hit appeared is cleared from need[]:
 * it no longer needs this segment, so the segment is complete once the
 * shapes that do need it are.
 */
static int cpu_segment(struct scratch *x, uint64_t s, int need[2],
                       struct hit **hits, size_t *nh, size_t *ch)
{
    const uint32_t K = KBLOCKS, nq = g_nq;
    const uint64_t W = g_W;
    const int J = g_J;
    uint64_t *bits = x->bits, *bq = x->bq, *cand = x->cand;

    segment_strip(s, need, hits, nh, ch);

    /*
     * Candidates are offsets from base.  A segment whose members all
     * stay below 2^64 is screened in one word, as before; the rest in
     * two.
     */
    u128 base = g_sieve_lo + (u128)s * g_segsize;
    uint64_t bmodW = mod_small(base, (uint32_t)W);
    u128 room = base <= g_hi ? g_hi - base : 0;
    uint64_t lim = room > UINT64_MAX ? UINT64_MAX : (uint64_t)room;
    int wide = base + g_segsize + (uint64_t)(g_k - 1) * g_k > UINT64_MAX;
    if (base <= g_hi)
        for (uint32_t i = 0; i < nq; i++)
            bq[i] = mod_small(base, g_q[i]);

    for (int d = 0; d < 2 && base <= g_hi; d++) {
        const struct shape *S = &g_sh[d];
        if (!need[d])
            continue;
        for (uint32_t j = 0; j < S->A; j++) {
            if (g_stop)
                return 0;
            /* a lower hit for this shape turned up: this segment can't win */
            if (g_first &&
                __atomic_load_n(&S->best_seg, __ATOMIC_ACQUIRE) < s) {
                need[d] = 0;
                break;
            }
            /* candidates p = base + s0 + k*W, k = 0..K-1 */
            uint64_t s0 = S->res[j] >= bmodW
                        ? S->res[j] - bmodW
                        : S->res[j] + W - bmodW;
            class_bitmap(d, bits, bq, s0, x->psrc, x->psh);

            uint32_t nc = 0;
            for (uint32_t w = 0; w < K / 64; w++) {
                uint64_t word = bits[w];
                while (word) {
                    uint32_t k = w * 64 + (uint32_t)__builtin_ctzll(word);
                    word &= word - 1;
                    uint64_t off = s0 + (uint64_t)k * W;
                    if (off <= lim)
                        cand[nc++] = off;
                }
            }
            /*
             * Screen one member at a time across the whole class, four
             * tests interleaved.  Members are all rough past the
             * pre-sieve limit, so trial division would be wasted;
             * members not struck beyond it go first.
             */
            for (int m = J; m < g_k && nc; m++)
                nc = screen_member(cand, nc, base, S->offs[m], wide);
            for (int m = 0; m < J && nc; m++)
                nc = screen_member(cand, nc, base, S->offs[m], wide);
            for (uint32_t ci = 0; ci < nc; ci++) {
                u128 p = base + cand[ci], head = p;
                int plus, len = d == INC
                    ? chain_verify_inc(p, &plus)
                    : chain_verify_dec(p, &head, &plus);
                if (len)
                    add_hit(hits, nh, ch, (struct hit){ head, p, len, plus, d });
            }
        }
    }
    return 1;
}

/* a finished segment: report its hits, or under --first each shape's best
   candidate, and advance the frontier */
static void segment_done(uint64_t s, const int *need, struct hit *hits,
                         size_t nh)
{
    uint64_t hist[2][NBUCKET] = { { 0 } };
    if (g_first) {
        for (int d = 0; d < 2; d++) {
            const struct hit *best = NULL;
            for (size_t i = 0; need[d] && i < nh; i++)
                if (hits[i].dir == d && (!best || hits[i].key < best->key))
                    best = &hits[i];
            if (best)
                record_candidate(d, s, best->key, best->len, best->plus);
        }
    } else {
        flush_hits(hits, nh, hist);
    }
    seg_complete(s, hist);
}

static void *worker(void *arg)
{
    (void)arg;
    struct scratch x;
    struct hit *hits = NULL;
    size_t nh, ch = 0;
    uint64_t s;
    int need[2];

    if (!scratch_init(&x))
        while (claim_segment(&s, need)) {
            nh = 0;
            if (cpu_segment(&x, s, need, &hits, &nh, &ch))
                segment_done(s, need, hits, nh);
        }
    free(hits);
    scratch_free(&x);
    __atomic_sub_fetch(&g_active, 1, __ATOMIC_RELEASE);
    return NULL;
}

/* ---------- wheel and sieve-prime setup ---------- */

static void wheel_init(void)
{
    for (int m = 0; m < g_k; m++) {
        g_sh[INC].offs[m] = (uint32_t)m * (m + 1);
        g_sh[DEC].offs[m] = (uint32_t)m * (uint32_t)(2 * g_k - 1 - m);
    }
    g_J = g_k < 6 ? g_k : 6;

    /*
     * Sized from INC.  DEC is INC's mirror image, so it has the same class
     * count at every wheel prime and the same wheel suits it; the segment
     * grid is therefore identical in every mode.
     */
    const uint32_t *offs = g_sh[INC].offs;
    static const uint32_t wp[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23 };
    uint64_t W = 1;
    uint64_t A = 1;
    for (size_t i = 0; i < sizeof wp / sizeof *wp; i++) {
        uint32_t q = wp[i];
        char seen[23] = { 0 };
        uint32_t v = 0;
        for (int m = 0; m < g_k; m++) {
            uint32_t r = offs[m] % q;
            if (!seen[r]) {
                seen[r] = 1;
                v++;
            }
        }
        if (v >= q)                   /* pattern inadmissible mod q: cannot */
            v = q - 1;                /* happen for m(m+1), but stay safe   */
        if (W * q > W_CAP || A * (q - v) > A_CAP)
            break;
        W *= q;
        A *= q - v;
    }
    g_W = W;

    unsigned char *bad = malloc(W);
    if (!bad)
        _exit(1);
    for (int d = 0; d < 2; d++) {
        struct shape *S = &g_sh[d];
        if (!S->on)
            continue;
        memset(bad, 0, W);
        for (size_t i = 0; i < sizeof wp / sizeof *wp && wp[i] <= W; i++) {
            uint32_t q = wp[i];
            if (W % q)
                continue;
            for (int m = 0; m < g_k; m++) {
                uint32_t o = S->offs[m] % q;
                for (uint64_t r = (q - o) % q; r < W; r += q)
                    bad[r] = 1;
            }
        }
        uint32_t a = 0;
        for (uint64_t r = 0; r < W; r++)
            a += !bad[r];
        S->A = a;
        S->res = malloc((size_t)(a ? a : 1) * sizeof *S->res);
        if (!S->res)
            _exit(1);
        a = 0;
        for (uint64_t r = 0; r < W; r++)
            if (!bad[r])
                S->res[a++] = (uint32_t)r;
    }
    free(bad);
}

static uint32_t powmod32(uint64_t b, uint32_t e, uint32_t q)
{
    uint64_t r = 1;
    b %= q;
    while (e) {
        if (e & 1)
            r = r * b % q;
        b = b * b % q;
        e >>= 1;
    }
    return (uint32_t)r;
}

static void sieve_primes_init(void)
{
    uint32_t B = g_depth;
    unsigned char *comp = calloc(B + 1, 1);
    if (!comp)
        _exit(1);
    for (uint32_t i = 2; (uint64_t)i * i <= B; i++)
        if (!comp[i])
            for (uint32_t j = i * i; j <= B; j += i)
                comp[j] = 1;
    uint32_t n = 0;
    for (uint32_t i = 2; i <= B; i++)
        if (!comp[i] && g_W % i)
            n++;
    g_q = malloc((n ? n : 1) * sizeof *g_q);
    g_qinv = malloc((n ? n : 1) * sizeof *g_qinv);
    if (!g_q || !g_qinv)
        _exit(1);
    for (int d = 0; d < 2; d++)
        if (g_sh[d].on && !(g_sh[d].dtab =
                malloc((size_t)(n ? n : 1) * g_J * sizeof *g_sh[d].dtab)))
            _exit(1);
    n = 0;
    for (uint32_t i = 2; i <= B; i++) {
        if (comp[i] || g_W % i == 0)
            continue;
        g_q[n] = i;
        g_qinv[n] = powmod32(g_W % i, i - 2, i);
        for (int d = 0; d < 2; d++)
            for (int oi = 0; g_sh[d].on && oi < g_J; oi++)
                g_sh[d].dtab[(size_t)n * g_J + oi] = (uint32_t)
                    ((uint64_t)g_qinv[n] * (g_sh[d].offs[oi] % i) % i);
        n++;
    }
    g_nq = n;
    free(comp);
}

static uint32_t inv_mod(uint64_t a, uint32_t q)
{
    return powmod32(a % q, q - 2, q);
}

/*
 * Candidate p = base + s0 + k*W has member p+o divisible by q exactly when
 *   k = c0 - d_o (mod q),  c0 = W^-1 * -(base+s0),  d_o = W^-1 * o.
 * Write x = k + r with r = -c0 = W^-1 * (base+s0) mod q: the condition
 * becomes x = -d_o (mod q), which no longer depends on the class.  So one
 * pattern per prime -- bit x clear iff x = -d_o for some offset o -- serves
 * every class, read from position r.  Primes are grouped so one pattern of
 * period P = q1*q2*... covers several at once, read from R = CRT(r_q).
 */
static void presieve_init(void)
{
    uint32_t n = 0;
    while (n < g_nq && g_q[n] <= g_presieve)
        n++;
    g_nps = n;
    g_crt = calloc(n ? n : 1, sizeof *g_crt);
    g_groups = calloc(n ? n : 1, sizeof *g_groups);
    unsigned char *seen = calloc(n ? g_q[n - 1] + 1 : 1, 1);
    if (!g_crt || !g_groups || !seen)
        _exit(1);

    uint32_t i = 0;
    while (i < n) {
        struct psgroup *G = &g_groups[g_ngroups++];
        uint64_t P = 1;
        G->first = i;
        while (i < n && P * g_q[i] <= g_pscap) {
            P *= g_q[i];
            i++;
        }
        if (i == G->first)           /* a prime above the cap: alone */
            P = g_q[i++];
        G->P = (uint32_t)P;
        G->n = i - G->first;

        for (uint32_t j = G->first; j < i; j++) {
            uint32_t q = g_q[j];
            uint64_t Mq = P / q;
            g_crt[j] = Mq * inv_mod(Mq, q) % P;
        }

        uint64_t nbits = P + KBLOCKS + 128;
        size_t words = (size_t)(nbits + 63) / 64;
        for (int d = 0; d < 2; d++) {
            const struct shape *S = &g_sh[d];
            if (!S->on)
                continue;
            uint64_t *pat = malloc(words * sizeof *pat);
            if (!pat)
                _exit(1);
            memset(pat, 0xFF, words * sizeof *pat);
            for (uint32_t j = G->first; j < i; j++) {
                uint32_t q = g_q[j];
                /* long chains repeat residues mod q; clear each only once */
                for (int m = 0; m < g_k; m++) {
                    uint32_t r = S->offs[m] % q;
                    if (seen[r])
                        continue;
                    seen[r] = 1;
                    uint32_t dd = (uint32_t)((uint64_t)g_qinv[j] * r % q);
                    for (uint64_t x = (q - dd) % q; x < nbits; x += q)
                        pat[x >> 6] &= ~(1ULL << (x & 63));
                }
                for (int m = 0; m < g_k; m++)
                    seen[S->offs[m] % q] = 0;
            }
            G->pat[d] = pat;
        }
    }
    free(seen);
}

/* ---------- progress / checkpoint / stats (as in gapchain_ps) ---------- */

static u128 pos_now(uint64_t f)
{
    if (!f)
        return g_lo;
    /* segments below g_nseg start at or below g_hi */
    if (f >= g_nseg)
        return g_hi;
    return g_sieve_lo + (u128)f * g_segsize;
}

/*
 * The status text, most important field first.  With width > 0, fields that
 * would overflow are dropped whole, never cut, so a narrow terminal cannot
 * show a truncated prime.
 */
static void format_progress(char *buf, size_t n, int width)
{
    uint64_t f, found[2], bseg[2];
    u128 bp[2];
    pthread_mutex_lock(&g_lock);
    f = g_frontier;
    for (int d = 0; d < 2; d++) {
        found[d] = g_sh[d].found;
        bseg[d] = g_sh[d].best_seg;
        bp[d] = g_sh[d].best_p;
    }
    pthread_mutex_unlock(&g_lock);

    double el = now_s() - g_t0;
    double pct = g_nseg ? 100.0 * (double)f / (double)g_nseg : 100.0;
    double rate = el > 0.0 ? (double)(f - g_base) * (double)g_segsize / el
                           : 0.0;
    double eta = rate > 0.0 ? (double)(g_nseg - f) * (double)g_segsize / rate
                            : -1.0;

    char sr[32], se[32], sl[32];
    fmt_rate(sr, sizeof sr, rate);
    fmt_dur(se, sizeof se, eta);
    fmt_dur(sl, sizeof sl, el);

    char fld[7][64];
    int nf = 0, both = g_sh[INC].on && g_sh[DEC].on;
    snprintf(fld[nf++], sizeof fld[0], "%.2f%%", pct);
    snprintf(fld[nf++], sizeof fld[0], "eta %s", se);
    for (int d = 0; d < 2; d++) {
        const char *tag = g_sh[d].tag;
        if (!g_sh[d].on)
            continue;
        if (!g_first)
            snprintf(fld[nf++], sizeof fld[0], "%" PRIu64 "%s%s chains",
                     found[d], both ? " " : "", both ? tag : "");
        else if (both)
            snprintf(fld[nf++], sizeof fld[0], bseg[d] == NO_SEG
                     ? "%s: no hit yet" : "%s: %s", tag, U128S(bp[d]));
        else if (bseg[d] == NO_SEG)
            snprintf(fld[nf++], sizeof fld[0], "no hit yet");
        else
            snprintf(fld[nf++], sizeof fld[0], "candidate %s", U128S(bp[d]));
    }
    snprintf(fld[nf++], sizeof fld[0], "%s", sr);
    snprintf(fld[nf++], sizeof fld[0], "at %.6g", (double)pos_now(f));
    snprintf(fld[nf++], sizeof fld[0], "elapsed %s", sl);

    size_t len = 0;
    buf[0] = '\0';
    for (int i = 0; i < nf && len < n; i++) {
        size_t need = (len ? 2 : 0) + strlen(fld[i]);
        if (width > 0 && len + need > (size_t)width)
            continue;
        len += (size_t)snprintf(buf + len, n - len, "%s%s", len ? "  " : "",
                                fld[i]);
    }
}

static const char *mode_name(void)
{
    return g_sh[INC].on && g_sh[DEC].on ? "both" : g_sh[DEC].on ? "dec" : "inc";
}

static void write_checkpoint(void)
{
    uint64_t frontier;
    struct shape snap[2];

    /*
     * One snapshot for frontier and --first candidate: a segment completes
     * only after its candidate is recorded, so any frontier past a hit is
     * written together with that hit.
     */
    pthread_mutex_lock(&g_lock);
    frontier = g_frontier;
    memcpy(snap, g_sh, sizeof snap);
    pthread_mutex_unlock(&g_lock);

    pthread_mutex_lock(&g_out);
    fflush(stdout);
    pthread_mutex_unlock(&g_out);

    FILE *f = fopen(g_cktmp, "w");
    if (!f) {
        note("gapsieve: %s: %s\n", g_cktmp, strerror(errno));
        return;
    }
    fprintf(f, "%s\n", CKPT_MAGIC);
    fprintf(f, "k %d\n", g_k);
    fprintf(f, "lo %s\n", U128S(g_lo));
    fprintf(f, "hi %s\n", U128S(g_hi));
    fprintf(f, "segsize %" PRIu64 "\n", g_segsize);
    fprintf(f, "nseg %" PRIu64 "\n", g_nseg);
    fprintf(f, "depth %u\n", g_depth);
    fprintf(f, "sievelo %s\n", U128S(g_sieve_lo));
    fprintf(f, "first %d\n", g_first ? 1 : 0);
    fprintf(f, "gaps %s\n", mode_name());
    fprintf(f, "done %" PRIu64 "\n", frontier);
    /* INC keeps the original line names; DEC's carry a "dec" prefix */
    for (int d = 0; d < 2; d++) {
        const struct shape *S = &snap[d];
        const char *pre = d == DEC ? "dec" : "";
        if (!S->on)
            continue;
        fprintf(f, "%sfound %" PRIu64 "\n", pre, S->found);
        for (int b = 0; b < NBUCKET; b++)
            fprintf(f, "%slen %d %" PRIu64 "\n", pre, g_k + b, S->hist[b]);
        if (S->best_seg != NO_SEG)
            fprintf(f, "%sbest %" PRIu64 " %s %d %d\n", pre, S->best_seg,
                    U128S(S->best_p), S->best_len, S->best_plus);
    }

    if (fflush(f) || fsync(fileno(f)) || fclose(f)) {
        note("gapsieve: %s: %s\n", g_cktmp, strerror(errno));
        return;
    }
    if (rename(g_cktmp, g_ckfile))
        note("gapsieve: %s: %s\n", g_ckfile, strerror(errno));
}

static int read_checkpoint(void)
{
    FILE *f = fopen(g_ckfile, "r");
    if (!f) {
        if (errno == ENOENT)
            return -1;
        note("gapsieve: %s: %s\n", g_ckfile, strerror(errno));
        exit(1);
    }
    char line[256];
    uint64_t segsize = 0, nseg = 0, done = NO_SEG;
    u128 lo = 0, hi = 0, slo = 0;
    uint64_t found[2] = { 0 }, hist[2][NBUCKET] = { { 0 } };
    uint64_t bseg[2] = { NO_SEG, NO_SEG };
    u128 bp[2] = { 0 };
    unsigned depth = 0;
    int k = 0, magic = 0, first = 0, blen[2] = { 0 }, bplus[2] = { 0 };
    int bad = 0;
    char gaps[16] = "inc";               /* files from before "gaps" */

    while (fgets(line, sizeof line, f)) {
        uint64_t v;
        int len;
        char num[64];                    /* 128-bit fields, parsed below */
        if (!strncmp(line, CKPT_MAGIC, strlen(CKPT_MAGIC)))
            magic = 1;
        else if (sscanf(line, "k %d", &k) == 1)
            ;
        else if (sscanf(line, "lo %63s", num) == 1)
            bad |= parse_num(num, &lo);
        else if (sscanf(line, "hi %63s", num) == 1)
            bad |= parse_num(num, &hi);
        else if (sscanf(line, "segsize %" SCNu64, &segsize) == 1)
            ;
        else if (sscanf(line, "nseg %" SCNu64, &nseg) == 1)
            ;
        else if (sscanf(line, "depth %u", &depth) == 1)
            ;
        else if (sscanf(line, "sievelo %63s", num) == 1)
            bad |= parse_num(num, &slo);
        else if (sscanf(line, "first %d", &first) == 1)
            ;
        else if (sscanf(line, "done %" SCNu64, &done) == 1)
            ;
        else if (sscanf(line, "gaps %15s", gaps) == 1)
            ;
        else if (sscanf(line, "found %" SCNu64, &found[INC]) == 1)
            ;
        else if (sscanf(line, "decfound %" SCNu64, &found[DEC]) == 1)
            ;
        else if (sscanf(line, "best %" SCNu64 " %63s %d %d",
                        &bseg[INC], num, &blen[INC], &bplus[INC]) == 4)
            bad |= parse_num(num, &bp[INC]);
        else if (sscanf(line, "decbest %" SCNu64 " %63s %d %d",
                        &bseg[DEC], num, &blen[DEC], &bplus[DEC]) == 4)
            bad |= parse_num(num, &bp[DEC]);
        else if (sscanf(line, "len %d %" SCNu64, &len, &v) == 2) {
            if (len >= k && len - k < NBUCKET)
                hist[INC][len - k] = v;
        }
        else if (sscanf(line, "declen %d %" SCNu64, &len, &v) == 2) {
            if (len >= k && len - k < NBUCKET)
                hist[DEC][len - k] = v;
        }
    }
    fclose(f);

    if (!magic || done == NO_SEG || bad) {
        note("gapsieve: %s: not a gapsieve checkpoint\n", g_ckfile);
        exit(1);
    }
    if (k != g_k || lo != g_lo || hi != g_hi ||
        first != (g_first ? 1 : 0) || strcmp(gaps, mode_name())) {
        note("gapsieve: %s describes a different search\n"
                "  checkpoint: -l %d -s %s -e %s --gaps %s%s\n"
                "  requested:  -l %d -s %s -e %s --gaps %s%s\n",
                g_ckfile, k, U128S(lo), U128S(hi), gaps,
                first ? " --first" : "", g_k, U128S(g_lo), U128S(g_hi),
                mode_name(), g_first ? " --first" : "");
        exit(1);
    }
    /*
     * The segment grid starts at sieve_lo, which below the depth depends on
     * the depth the file was written with, so adopt the file's grid (files
     * from before the "sievelo" line derive it from their depth).  Only the
     * depth must stay below it: the sieve is exact only while every member
     * exceeds every sieve prime.
     */
    if (!slo)
        slo = (u128)depth + 512 > lo ? (u128)depth + 512 : lo;
    if (slo <= (u128)g_depth + 256) {
        note("gapsieve: %s was cut at %s; resume it with "
                "--depth %s or less\n", g_ckfile, U128S(slo),
                U128S(slo - 512));
        exit(1);
    }
    g_sieve_lo = slo;
    g_nseg = g_hi >= slo ? (uint64_t)((g_hi - slo) / g_segsize) + 1 : 1;
    if (segsize != g_segsize || nseg != g_nseg || done > nseg) {
        note("gapsieve: %s: segment grid mismatch (checkpoint "
                "segsize %" PRIu64 ", this build %" PRIu64 ") -- was it "
                "written by a different gapsieve version?\n",
                g_ckfile, segsize, g_segsize);
        exit(1);
    }
    (void)depth;
    /*
     * --first hits recorded before the stop.  Each hit's segment was fully
     * scanned for its shape, so it stands until one of the unfinished
     * segments below it beats it.
     */
    for (int d = 0; d < 2; d++) {
        struct shape *S = &g_sh[d];
        if (bseg[d] != NO_SEG) {
            if (bseg[d] >= nseg || bp[d] < g_lo || bp[d] > g_hi ||
                blen[d] < g_k) {
                note("gapsieve: %s: corrupt \"%sbest\" line\n", g_ckfile,
                     d == DEC ? "dec" : "");
                exit(1);
            }
            S->best_p = bp[d];
            S->best_len = blen[d];
            S->best_plus = bplus[d];
            __atomic_store_n(&S->best_seg, bseg[d], __ATOMIC_RELEASE);
        }
        S->found = found[d];
        memcpy(S->hist, hist[d], sizeof S->hist);
    }
    g_frontier = g_next_seg = done;
    return 0;
}

static void print_first_result(void)
{
    uint64_t f = g_frontier;
    double pct = g_nseg ? 100.0 * (double)f / (double)g_nseg : 100.0;
    char el[32];
    int open = 0;
    fmt_dur(el, sizeof el, now_s() - g_t0);

    for (int d = 0; d < 2; d++) {
        const struct shape *S = &g_sh[d];
        char what[64];
        if (!S->on)
            continue;
        if (g_sh[DEC].on)
            snprintf(what, sizeof what, "%s-gap chain (%s)",
                     d == INC ? "increasing" : "decreasing", S->oeis);
        else
            snprintf(what, sizeof what, "chain");
        /* INC counts longer chains; a DEC pattern has exactly L primes */
        const char *atleast = d == INC ? ">= " : "";

        if (S->best_seg == NO_SEG) {
            if (g_stop)
                note("\nstopped at %.2f%% of %s..%s with no %s of length "
                     "%s%d yet, %s elapsed\n", pct, U128S(g_lo),
                     U128S(g_hi), what, atleast, g_k, el);
            else
                note("\nno %s of length %s%d in %s..%s, %s elapsed\n", what,
                     atleast, g_k, U128S(g_lo), U128S(g_hi), el);
            open |= g_stop;
            continue;
        }

        /* DEC prints the L-prime pattern itself: its head is the term */
        struct hit h = { S->best_p, S->best_p,
                         d == INC ? S->best_len : g_k,
                         d == INC ? S->best_plus : 0, d };
        print_hit(&h);
        fflush(stdout);
        if (g_log) {
            char ts[32];
            log_time(ts, sizeof ts);
            pthread_mutex_lock(&g_out);
            fprintf(g_log, "%s  ", ts);
            fprint_hit(g_log, &h);
            fflush(g_log);
            pthread_mutex_unlock(&g_out);
        }
        if (f >= S->best_seg) {
            note("\nfirst %s of length %s%d is at %s, after searching %.2f%% "
                 "of the range in %s\n", what, atleast, g_k,
                 U128S(S->best_p), pct, el);
        } else {
            note("\ncandidate %s of length %s%d at %s is not final "
                 "yet: %" PRIu64 " segments below it are unfinished\n",
                 what, atleast, g_k, U128S(S->best_p), S->best_seg - f);
            open = 1;
        }
        if (d == DEC && S->best_len > g_k) {
            u128 head;
            int plus;
            chain_verify_dec(S->best_p, &head, &plus);
            note("(the pattern ends a longer decreasing chain of %d primes "
                 "from %s)\n", S->best_len, U128S(head));
        }
    }
    if (open && g_ckfile)
        note("resume with: -c %s -r\n", g_ckfile);
}

static void print_stats(void)
{
    uint64_t f;
    struct shape snap[2];

    if (g_first) {
        print_first_result();
        return;
    }
    pthread_mutex_lock(&g_lock);
    f = g_frontier;
    memcpy(snap, g_sh, sizeof snap);
    pthread_mutex_unlock(&g_lock);

    double pct = g_nseg ? 100.0 * (double)f / (double)g_nseg : 100.0;
    char el[32];
    fmt_dur(el, sizeof el, now_s() - g_t0);

    note("\n%s %.2f%% of %s..%s (%" PRIu64 "/%" PRIu64 " segments), "
            "%s elapsed\n", f >= g_nseg ? "completed" : "stopped at", pct,
            U128S(g_lo), U128S(g_hi), f, g_nseg, el);

    for (int d = 0; d < 2; d++) {
        const struct shape *S = &snap[d];
        if (!S->on)
            continue;
        if (g_sh[DEC].on)
            note("%s gaps (%s):\n", d == INC ? "increasing" : "decreasing",
                 S->oeis);
        if (!S->found) {
            note("no chains found\n");
            continue;
        }
        note("chains by length:\n");
        for (int b = 0; b < NBUCKET; b++) {
            if (!S->hist[b])
                continue;
            char tag[16];
            snprintf(tag, sizeof tag, "%d%s", g_k + b,
                     b == NBUCKET - 1 ? "+" : "");
            note("  %6s  %14" PRIu64 "  %6.2f%%\n", tag, S->hist[b],
                    100.0 * (double)S->hist[b] / (double)S->found);
        }
        note("  %6s  %14" PRIu64 "\n", "total", S->found);
    }
    if (f < g_nseg && g_ckfile)
        note("resume with: -c %s -r\n", g_ckfile);
    else if (f < g_nseg)
        note("(no checkpoint file; rerun with -c FILE to make "
                        "interrupted runs resumable)\n");
}

static void on_signal(int sig)
{
    (void)sig;
    if (g_stop)
        _exit(130);
    g_stop = 1;
}

static void monitor(void)
{
    double last_ck = now_s(), last_log = now_s();
    for (;;) {
        struct timespec ts = { 0, 250000000L };
        nanosleep(&ts, NULL);
        int live = __atomic_load_n(&g_active, __ATOMIC_ACQUIRE);
        double t = now_s();
        int due = t - last_log >= 60.0;
        if (g_progress || due) {
            char line[256], full[256];
            int w = term_width() - 1;
            if (g_progress)
                format_progress(line, sizeof line, w);
            if (due)
                format_progress(full, sizeof full, 0);
            pthread_mutex_lock(&g_out);
            if (g_progress)             /* one screen line, old text blanked */
                fprintf(stderr, "\r%-*.*s", w, w, line);
            else if (due && !g_quiet)
                fprintf(stderr, "%s\n", full);
            fflush(stderr);
            if (due)
                log_locked(full);
            pthread_mutex_unlock(&g_out);
        }
        if (due)
            last_log = t;
        if (g_ckfile && g_ckint && t - last_ck >= (double)g_ckint) {
            write_checkpoint();
            last_ck = t;
        }
        if (!live || g_stop)
            break;
    }
}

#ifdef GAPSIEVE_GPU
/* the Metal GPU search, in gapsieve_gpu.m, which includes this file */
static const char *gpu_setup(void);
static const char *gpu_name(void);
static void *gpu_worker(void *arg);
#endif

/* never returns, so option cases that call it need no break */
static __attribute__((noreturn)) void usage(FILE *f, const char *prog,
                                            int status)
{
    fprintf(f,
        "usage: %s -l LEN -e END [options]\n"
        "\n"
        "Constellation-sieve search for runs of consecutive primes whose\n"
        "gaps grow by two: 2, 4, 6, ... (OEIS A016045) or, mirrored,\n"
        "..., 6, 4, 2 (A263049).  Same search and options as gapchain_ps,\n"
        "but much faster at large heights.\n"
        "\n"
        "required:\n"
        "  -l, --length N      chain length: least number of primes in a run\n"
        "  -e, --end N         highest prime to consider as a chain head;\n"
        "                      max = as high as the proven primality test\n"
        "                      allows (about 3.3e24)\n"
        "\n"
        "options:\n"
        "  -s, --start N       lowest prime to consider (default 3)\n"
        "  -t, --threads N     worker threads (default: online CPUs)\n"
        "  -g, --gaps WHICH    inc (2, 4, 6, ...; default), dec (..., 6, 4, 2)\n"
        "                      or both; both costs the time of two runs\n"
        "  -f, --first         report only the smallest qualifying chain\n"
        "                      (of each shape), then stop\n"
        "  -c, --checkpoint F  checkpoint to F periodically and on SIGINT\n"
        "  -i, --interval S    checkpoint seconds, 0 = on exit (default 60)\n"
        "  -r, --resume        resume from F (required if F exists)\n"
        "  -d, --depth N       largest sieve prime (default %u)\n"
        "  -p, --presieve N    primes up to N sieved with all L offsets via\n"
        "                      patterns (default %u); -d and -p are tuning\n"
        "                      only -- results are identical at any values\n"
        "      --log FILE      also append every message, timestamped, to\n"
        "                      FILE, with a status line once a minute\n"
        "      --gpu           also search on the GPU (Mac and CUDA builds);\n"
        "                      -t then sets the CPU threads beside it\n"
        "                      (default all but two cores, 0 for the GPU\n"
        "                      alone)\n"
        "  -q, --quiet         suppress progress output\n"
        "  -h, --help          this message\n"
        "\n"
        "N accepts 1e14, 2.5e15, 300T, 1_000_000 or plain digits.\n",
        prog, g_depth, g_presieve);
    exit(status);
}

int main(int argc, char **argv)
{
    static const struct option longopts[] = {
        { "length",     required_argument, NULL, 'l' },
        { "start",      required_argument, NULL, 's' },
        { "end",        required_argument, NULL, 'e' },
        { "threads",    required_argument, NULL, 't' },
        { "first",      no_argument,       NULL, 'f' },
        { "gaps",       required_argument, NULL, 'g' },
        { "checkpoint", required_argument, NULL, 'c' },
        { "interval",   required_argument, NULL, 'i' },
        { "resume",     no_argument,       NULL, 'r' },
        { "depth",      required_argument, NULL, 'd' },
        { "presieve",   required_argument, NULL, 'p' },
        { "log",        required_argument, NULL, 'L' },   /* long form only */
        { "gpu",        no_argument,       NULL, 'G' },   /* long form only */
        { "quiet",      no_argument,       NULL, 'q' },
        { "help",       no_argument,       NULL, 'h' },
        { NULL,         0,                 NULL,  0  },
    };
    const char *prog = argv[0];
    const char *arg_end = NULL, *logfile = NULL;
    long k = 0, nthreads = -1;          /* -1: not given */
    int resume = 0, opt;

    char cmdline[1024];
    size_t cl = 0;
    for (int i = 0; i < argc && cl < sizeof cmdline; i++)
        cl += (size_t)snprintf(cmdline + cl, sizeof cmdline - cl, "%s%s",
                               i ? " " : "", argv[i]);

    g_lo = 3;
    g_sh[INC].on = 1;
    while ((opt = getopt_long(argc, argv, "l:s:e:t:fg:c:i:rd:p:qh",
                              longopts, NULL)) != -1) {
        switch (opt) {
        case 'l':
            k = strtol(optarg, NULL, 10);
            if (k < 2 || k > 1000) {
                fprintf(stderr, "gapsieve: --length must be 2..1000\n");
                return 2;
            }
            break;
        case 's':
            if (parse_num(optarg, &g_lo)) {
                fprintf(stderr, "gapsieve: bad --start value '%s'\n", optarg);
                return 2;
            }
            break;
        case 'e':
            arg_end = optarg;
            if (parse_num(optarg, &g_hi)) {
                fprintf(stderr, "gapsieve: bad --end value '%s'\n", optarg);
                return 2;
            }
            break;
        case 't':
            nthreads = strtol(optarg, NULL, 10);
            if (nthreads < 0) {
                fprintf(stderr, "gapsieve: --threads must be >= 0\n");
                return 2;
            }
            break;
        case 'G': g_gpu = 1; break;
        case 'f': g_first = 1; break;
        case 'g':
            g_sh[INC].on = !strcmp(optarg, "inc") || !strcmp(optarg, "both");
            g_sh[DEC].on = !strcmp(optarg, "dec") || !strcmp(optarg, "both");
            if (!g_sh[INC].on && !g_sh[DEC].on) {
                fprintf(stderr, "gapsieve: --gaps must be inc, dec or both\n");
                return 2;
            }
            break;
        case 'c': g_ckfile = optarg; break;
        case 'i': g_ckint = (unsigned)strtoul(optarg, NULL, 10); break;
        case 'r': resume = 1; break;
        case 'd': {
            uint64_t d;
            if (parse_u64(optarg, &d) || d < 32 || d > (1u << 27)) {
                fprintf(stderr, "gapsieve: --depth must be 32..%u\n",
                        1u << 27);
                return 2;
            }
            g_depth = (uint32_t)d;
            break;
        }
        case 'p': {
            uint64_t d;
            if (parse_u64(optarg, &d) || d > (1u << 16)) {
                fprintf(stderr, "gapsieve: --presieve must be 0..%u\n",
                        1u << 16);
                return 2;
            }
            g_presieve = (uint32_t)d;
            break;
        }
        case 'L': logfile = optarg; break;
        case 'q': g_quiet = 1; break;
        case 'h': usage(stdout, prog, 0);
        default:  usage(stderr, prog, 2);
        }
    }
    if (optind < argc) {
        fprintf(stderr, "gapsieve: unexpected argument '%s'\n", argv[optind]);
        return 2;
    }
    if (!k || !arg_end) {
        fprintf(stderr, "gapsieve: %s is required\n",
                !k ? "--length" : "--end");
        usage(stderr, prog, 2);
    }
    if (resume && !g_ckfile) {
        fprintf(stderr, "gapsieve: --resume needs --checkpoint FILE\n");
        return 2;
    }
    g_k = (int)k;
    if (g_lo < 3)
        g_lo = 3;
    /*
     * Everything checked for a head p lies below p + reach: members and
     * interlopers reach p + (L-1)L, and the increasing extension walk
     * p + (L+EXT_CAP)(L+EXT_CAP+1).  All of it must stay below PSI13, where
     * the primality test is proven.
     */
    uint64_t reach = g_sh[INC].on
        ? (uint64_t)(g_k + EXT_CAP) * (uint64_t)(g_k + EXT_CAP + 1)
        : (uint64_t)(g_k - 1) * (uint64_t)g_k;
    int lowered = g_hi > PSI13 - 1 - reach;
    if (lowered)
        g_hi = PSI13 - 1 - reach;
    if (g_hi < g_lo) {
        fprintf(stderr, "gapsieve: --end is below --start\n");
        return 2;
    }
    if (nthreads == 0 && !g_gpu) {
        fprintf(stderr, "gapsieve: --threads 0 needs --gpu\n");
        return 2;
    }
    if (nthreads < 0) {
        nthreads = sysconf(_SC_NPROCESSORS_ONLN);
        if (g_gpu)          /* leave room for the GPU's feeder thread */
            nthreads -= 2;
        if (nthreads < 1)
            nthreads = 1;
    }
#ifndef GAPSIEVE_GPU
    if (g_gpu) {
        fprintf(stderr, "gapsieve: this build has no GPU support; on a Mac, "
                        "build it with make, and with an NVIDIA GPU, with "
                        "make CUDA=1\n");
        return 2;
    }
#endif

    wheel_init();
    sieve_primes_init();
    presieve_init();

    /*
     * The sieve starts above the depth so a member can never coincide with a
     * sieving prime; the strip below is searched exactly by direct_scan.
     */
    u128 D = (u128)g_depth + 512;
    if (g_lo > D)
        D = g_lo;
    g_sieve_lo = D;
    g_segsize = (uint64_t)KBLOCKS * g_W;
    g_nseg = g_hi >= g_sieve_lo
           ? (uint64_t)((g_hi - g_sieve_lo) / g_segsize) + 1 : 1;
    g_progress = isatty(STDERR_FILENO) && !g_quiet;

    /* appended, so a resumed run continues the same log */
    if (logfile) {
        g_log = fopen(logfile, "a");
        if (!g_log) {
            fprintf(stderr, "gapsieve: %s: %s\n", logfile, strerror(errno));
            return 2;
        }
        log_note("started (pid %ld): %s\n", (long)getpid(), cmdline);
    }
    if (lowered)
        note("gapsieve: --end lowered to %s, the last start whose chain can "
             "be checked below 3.317e24, the limit of the proven primality "
             "test\n", U128S(g_hi));

    if (g_ckfile) {
        size_t n = strlen(g_ckfile) + 5;
        g_cktmp = malloc(n);
        if (!g_cktmp)
            return 1;
        snprintf(g_cktmp, n, "%s.tmp", g_ckfile);
        if (!resume && !access(g_ckfile, F_OK)) {
            note("gapsieve: %s exists; pass --resume to resume "
                 "from it or remove it first\n", g_ckfile);
            return 2;
        }
        if (resume && read_checkpoint() == 0) {
            /* a --first hit with nothing unfinished below it is final */
            int settled = g_first;
            for (int d = 0; d < 2; d++)
                if (g_sh[d].on && g_frontier < g_sh[d].best_seg)
                    settled = 0;
            if (g_frontier >= g_nseg || settled) {
                note("gapsieve: %s is already complete\n", g_ckfile);
                g_t0 = now_s();
                print_stats();
                log_note("exit 0\n");
                return 0;
            }
            char more[400] = "";
            size_t ml = 0;
            for (int d = 0; d < 2 && ml < sizeof more; d++) {
                const struct shape *S = &g_sh[d];
                const char *tag = g_sh[DEC].on ? S->tag : "";
                const char *sp = g_sh[DEC].on ? " " : "";
                if (!S->on)
                    continue;
                if (!g_first)
                    ml += (size_t)snprintf(more + ml, sizeof more - ml,
                            ", %" PRIu64 "%s%s chains so far", S->found,
                            sp, tag);
                else if (S->best_seg != NO_SEG && g_frontier < S->best_seg)
                    ml += (size_t)snprintf(more + ml, sizeof more - ml,
                            "; %s%scandidate %s recorded in segment "
                            "%" PRIu64 ", finishing the %" PRIu64 " segments "
                            "below it", tag, sp, U128S(S->best_p),
                            S->best_seg, S->best_seg - g_frontier);
                else if (S->best_seg != NO_SEG)
                    ml += (size_t)snprintf(more + ml, sizeof more - ml,
                            "; %s%sanswer %s already settled", tag,
                            sp, U128S(S->best_p));
            }
            note("resuming from %s at %" PRIu64 "/%" PRIu64 " segments "
                 "(%.2f%%)%s\n", g_ckfile, g_frontier, g_nseg,
                 100.0 * (double)g_frontier / (double)g_nseg, more);
        }
    }

    char who[200];
    snprintf(who, sizeof who, "%ld thread%s", nthreads, nthreads == 1 ? "" : "s");
#ifdef GAPSIEVE_GPU
    if (g_gpu) {
        const char *why = gpu_setup();
        if (why) {
            note("gapsieve: cannot use the GPU: %s\n", why);
            log_note("exit 1\n");
            return 1;
        }
        snprintf(who, sizeof who, "GPU (%s) plus %ld CPU thread%s", gpu_name(),
                 nthreads, nthreads == 1 ? "" : "s");
    }
#endif

    note(
        "consecutive primes%s, chain length >= %d, range %s..%s, %s%s\n"
        "(sieve: wheel %" PRIu64 " with %u classes%s; %u primes pre-sieved in "
        "%u patterns, %d offsets to depth %u; segment %.3g)\n",
        !g_sh[DEC].on ? "" : g_sh[INC].on
            ? ", increasing and decreasing gaps" : ", decreasing gaps",
        g_k, U128S(g_lo), U128S(g_hi), who,
        g_first ? ", first match only" : "",
        g_W, g_sh[g_sh[INC].on ? INC : DEC].A,
        g_sh[INC].on && g_sh[DEC].on ? " per shape" : "",
        g_nps, g_ngroups, g_J, g_depth, (double)g_segsize);

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_signal;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    sigset_t block, prev;
    sigemptyset(&block);
    sigaddset(&block, SIGINT);
    sigaddset(&block, SIGTERM);
    pthread_sigmask(SIG_BLOCK, &block, &prev);

    g_t0 = now_s();
    g_base = g_frontier;

    long nt = nthreads + (g_gpu ? 1 : 0);   /* the GPU is fed by one more */
    pthread_t *t = calloc((size_t)nt, sizeof *t);
    if (!t)
        return 1;
    long started = 0;
    __atomic_store_n(&g_active, (int)nt, __ATOMIC_RELEASE);
    for (long i = 0; i < nt; i++) {
        void *(*fn)(void *) = worker;
#ifdef GAPSIEVE_GPU
        if (g_gpu && i == 0)
            fn = gpu_worker;
#endif
        if (pthread_create(&t[i], NULL, fn, NULL)) {
            __atomic_sub_fetch(&g_active, (int)(nt - i), __ATOMIC_RELEASE);
            break;
        }
        started++;
    }
    if (!started) {
        note("gapsieve: could not start any worker threads\n");
        return 1;
    }
    pthread_sigmask(SIG_SETMASK, &prev, NULL);
    monitor();

    if (g_stop)
        note("interrupt: finishing segments in flight...\n");
    for (long i = 0; i < started; i++)
        pthread_join(t[i], NULL);
    free(t);

    if (g_ckfile)
        write_checkpoint();
    if (g_progress)
        fprintf(stderr, "\r%*s\r", term_width() - 1, "");
    fflush(stdout);
    print_stats();
    log_note("exit %d\n", g_stop ? 130 : 0);
    free(g_pend);
    free(g_cktmp);
    free(g_q);
    free(g_qinv);
    for (int d = 0; d < 2; d++) {
        free(g_sh[d].res);
        free(g_sh[d].dtab);
        for (uint32_t i = 0; i < g_ngroups; i++)
            free(g_groups[i].pat[d]);
    }
    free(g_groups);
    free(g_crt);
    return g_stop ? 130 : 0;
}
