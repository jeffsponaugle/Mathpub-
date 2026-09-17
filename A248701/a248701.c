/*
 * a248701.c
 *
 * Compute and extend OEIS A248701 and its three companions A248702-A248704:
 *
 *   A248701  "Smallest prime such that the preceding n prime gaps are increasing
 *             and the following n prime gaps are decreasing."
 *             Known terms (n = 1..7): 3, 7, 359, 7853, 96401, 2812099, 294276293.
 *   A248702  the same with decreasing gaps before and increasing gaps after (a valley);
 *             offset 0, a(0) = 2.  Known (n = 0..6): 2, 3, 19, 43, 2687, 179819, 1107791,
 *             "a(7) >= 8960453 if it exists".
 *   A248703  "Strict Peak Primes": strictly increasing gaps before, strictly decreasing
 *             after, and the index counts the strict steps, so a(n) needs n+1 gaps on
 *             each side.  Known (n = 1..6): 23, 1439, 21433, 1130863, 19881311, 331542583.
 *   A248704  strict valley, n gaps each side.  Known (n = 1..6): 3, 19, 1429, 25243,
 *             340577, 1107791.
 *
 * Definitions
 * -----------
 * Let p_1 = 2 < p_2 < ... be the primes and g_i = p_(i+1) - p_i.  The prime p_k
 * sits between the gap g_(k-1) into it and the gap g_k out of it.  Its peak depth
 * is min(L, M), where L is the length of the longest non-decreasing run of gaps
 * that ends with g_(k-1) and M the length of the longest non-increasing run that
 * starts with g_k.  "Increasing" in the OEIS names is non-strict (a(2) = 7 has
 * gaps 2, 2 | 4, 2) and nothing is required of g_(k-1) against g_k; this is what
 * the PARI and Maple programs in the entry test.  A248701(n) is the smallest
 * prime of peak depth >= n.  Example: 359 lies in 337, 347, 349, 353, 359, 367,
 * 373, 379, 383 with gaps 10, 2, 4, 6 | 8, 6, 6, 4: L = 3, M = 4, depth 3.
 * Valley depth swaps the two directions, the strict depths use strict
 * inequalities.  A248702(n) = smallest prime of valley depth >= n, A248703(n) =
 * smallest prime of strict peak depth >= n+1, A248704(n) = smallest prime of
 * strict valley depth >= n.  Each a(n) is non-decreasing in n.
 *
 * Method
 * ------
 * The range is cut into chunks that worker threads take from an atomic counter.
 * A thread sieves its chunk with primesieve, starting a few hundred primes below
 * it and running a few hundred primes past it so that the run lengths at the
 * chunk edges are exact, and feeds the gaps to an O(1)-per-gap state machine:
 * the lengths of the non-decreasing, non-increasing, strictly increasing and
 * strictly decreasing runs ending at the current gap, plus rings holding the
 * last 256 of those lengths and the last 256 primes.  When the run ending at gap
 * j has length r, every n <= r names the centre k = j - n + 1 whose n following
 * gaps are that run's tail, and the centre has depth >= n iff the opposite run
 * ending at gap k - 1 has length >= n.  Each (centre, n) pair is examined once,
 * so a chunk yields, for each mode and depth d, the number of centres of depth
 * >= d and the smallest one.  Only runs of length >= DMIN (default 4) are
 * examined above 2^30: shorter depths are all known and examining them would
 * cost as much as the sieve itself.  Chunk results are folded into the global
 * totals in increasing order (a completion frontier), so the first record folded
 * is the smallest, the totals always describe an initial segment of the range,
 * and a checkpoint is just the frontier plus the totals.
 *
 * Heuristics: for independent continuous gaps P(depth >= d) = 1/(d!)^2.  The
 * observed ratio count(>= d) / count(>= d+1) is about 0.8 d^2 (ties help the
 * non-strict runs).  Below 4e8 there are two centres of peak depth >= 7, so one
 * expects about one centre of depth >= 8 per 5e8 primes, >= 9 per 3e10, >= 10
 * per 2e12 and >= 11 per 2e14: a(8) ~ 1e10, a(9) ~ 1e12, a(10) ~ 1e14,
 * a(11) ~ 1e16, each uncertain by a factor of several either way.
 *
 * Usage
 * -----
 *   a248701 scan [START] END [-t T] [-c CHUNK] [-n N] [-d DMIN] [-r DEPTH] [-S STATE] [-i SECS] [-q]
 *       Exhaustive scan of the primes in [START, END] (START defaults to 0).
 *       Prints every new record of depth >= DEPTH (default 6) in any of the four
 *       modes as it is found, then the tables of a(n) and of the depth counts.
 *       -n N stops as soon as A248701(N) is known.  -S FILE keeps a checkpoint
 *       (every SECS seconds, default 60, and on Ctrl-C); rerunning the same
 *       command resumes from it.
 *   a248701 show P [-n N]
 *       The depths of the prime P in the four modes and the primes and gaps
 *       around it.
 *   a248701 selftest [-t T]
 *
 * Numbers may be written as decimal, 2^k, 10^k, 1e12, or X+Y / X-Y of those.
 *
 * Build:  cc -O2 -std=gnu11 -pthread -I/opt/homebrew/include a248701.c \
 *            -L/opt/homebrew/lib -lprimesieve -lm -o a248701
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

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint8_t  u8;
typedef unsigned __int128 u128;

#define NMAX    63              /* largest depth tracked (63 monotone gaps in a row never happen below 2^64) */
#define RUNCAP  250             /* run lengths saturate here; > NMAX and fits a byte */
#define RMASK   255             /* rings of the last 256 gaps */
#define WARM    256             /* primes sieved before and after a chunk so that its edge runs are exact */
#define SMALL   (1ull << 30)    /* chunks starting below this track every depth whatever -d says; holds all known terms */
#define RING    1024            /* chunk bookkeeping ring; bounds how far threads may run ahead of the frontier */
#define NMODE   4

enum { PEAK, VALLEY, SPEAK, SVALLEY };
static const char *SEQ[NMODE]    = { "A248701", "A248702", "A248703", "A248704" };
static const char *MODE[NMODE]   = { "peak", "valley", "strict peak", "strict valley" };
static const int   SHIFT[NMODE]  = { 0, 0, 1, 0 };     /* OEIS index n <-> depth n + SHIFT */
static const int   FIRSTN[NMODE] = { 1, 0, 1, 1 };     /* OEIS offsets */
#define KDEPTH 7
static const u64 KNOWN[NMODE][KDEPTH + 1] = {           /* known terms indexed by depth */
    { 0, 3, 7, 359, 7853, 96401, 2812099, 294276293 },
    { 2, 3, 19, 43, 2687, 179819, 1107791, 0 },
    { 0, 0, 23, 1439, 21433, 1130863, 19881311, 331542583 },
    { 0, 3, 19, 1429, 25243, 340577, 1107791, 0 },
};

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
    fputs(stderr_tty ? "\r\033[Ka248701: " : "a248701: ", stderr);
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

static void *xcalloc(size_t n, size_t sz)
{
    void *p = calloc(n ? n : 1, sz);
    if (!p) die("out of memory");
    return p;
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
/* Scanning one chunk                                                  */
/* ------------------------------------------------------------------ */

typedef struct {
    u64 first[NMODE][NMAX + 1];     /* smallest centre in the chunk of depth >= d, 0 if none */
    u64 count[NMODE][NMAX + 1];     /* centres in the chunk of depth >= d (d >= dmin) */
    u64 nprimes;                    /* primes in [lo, hi) */
    u64 maxgap;                     /* largest gap starting in [lo, hi) */
} result_t;

typedef struct {
    u8  nd[RMASK + 1], ni[RMASK + 1], si[RMASK + 1], sd[RMASK + 1];   /* run lengths ending at gap j */
    u64 p[RMASK + 1];                                                 /* prime at the start of gap j */
} rings_t;

/* The run of kind "fwd" ending at gap j has length run: for every n in [dmin, run]
 * the centre k = j - n + 1 is followed by n gaps of that run; it has depth >= n
 * in this mode iff the opposite run ending at gap k - 1 (ring bwd) is >= n. */
static inline void examine(result_t *r, const rings_t *R, int mode, const u8 *bwd,
                           unsigned run, u64 j, unsigned dmin, u64 lo, u64 hi)
{
    unsigned top = run < NMAX ? run : NMAX;
    for (unsigned n = dmin; n <= top; n++) {
        if (n > j) break;                       /* the centre needs a gap before it */
        u64 k = j - n + 1;
        if (bwd[(k - 1) & RMASK] < n) continue;
        u64 P = R->p[k & RMASK];
        if (P < lo || P >= hi) continue;
        r->count[mode][n]++;
        if (!r->first[mode][n]) r->first[mode][n] = P;
    }
}

/* primes in [lo, hi) as centres; sieving starts WARM primes below lo and ends WARM primes above hi */
static void scan_chunk(primesieve_iterator *it, rings_t *R, u64 lo, u64 hi, unsigned dmin, result_t *r)
{
    u64 margin = 1ull << 16;
    for (;;) {
        memset(r, 0, sizeof *r);
        u64 start = lo > margin ? lo - margin : 0;
        primesieve_jump_to(it, start, hi + margin);
        u64 prev = primesieve_next_prime(it);
        u64 nbefore = 0, nafter = 0, j = 0, pg = 0;
        unsigned nd = 1, ni = 1, si = 1, sd = 1;
        if (prev < lo) nbefore++;
        else if (prev < hi) r->nprimes++;
        for (;;) {
            u64 q = primesieve_next_prime(it);
            u64 g = q - prev;
            if (j) {
                nd = g >= pg ? nd + 1 : 1;
                ni = g <= pg ? ni + 1 : 1;
                si = g >  pg ? si + 1 : 1;
                sd = g <  pg ? sd + 1 : 1;
                if (nd > RUNCAP) nd = RUNCAP;
                if (ni > RUNCAP) ni = RUNCAP;
                if (si > RUNCAP) si = RUNCAP;
                if (sd > RUNCAP) sd = RUNCAP;
            }
            unsigned x = (unsigned)j & RMASK;
            R->nd[x] = (u8)nd; R->ni[x] = (u8)ni; R->si[x] = (u8)si; R->sd[x] = (u8)sd;
            R->p[x] = prev;
            if (ni >= dmin) examine(r, R, PEAK,    R->nd, ni, j, dmin, lo, hi);
            if (nd >= dmin) examine(r, R, VALLEY,  R->ni, nd, j, dmin, lo, hi);
            if (sd >= dmin) examine(r, R, SPEAK,   R->si, sd, j, dmin, lo, hi);
            if (si >= dmin) examine(r, R, SVALLEY, R->sd, si, j, dmin, lo, hi);
            if (g > r->maxgap && prev >= lo && prev < hi) r->maxgap = g;
            if (q < lo) nbefore++;
            else if (q < hi) r->nprimes++;
            else if (++nafter > WARM) break;
            pg = g; prev = q; j++;
        }
        if (it->is_error) die("primesieve error in [%" PRIu64 ", %" PRIu64 ")", lo, hi);
        if (start > 0 && nbefore < WARM) {
            margin *= 4;
            if (margin > (1ull << 32)) die("fewer than %d primes within 2^32 below %" PRIu64 "?!", WARM, lo);
            continue;
        }
        return;
    }
}

/* Straightforward reference implementation (arrays, one walk per centre); selftest only. */
static void reference_scan(u64 lo, u64 hi, result_t *r)
{
    size_t n;
    u64 margin = 1ull << 16;
    u64 *pr = primesieve_generate_primes(lo > margin ? lo - margin : 0, hi + margin, &n, UINT64_PRIMES);
    memset(r, 0, sizeof *r);
    if (n < 3) { primesieve_free(pr); return; }
    for (size_t i = 0; i < n; i++)
        if (pr[i] >= lo && pr[i] < hi) r->nprimes++;
    for (size_t i = 0; i + 1 < n; i++)
        if (pr[i] >= lo && pr[i] < hi && pr[i + 1] - pr[i] > r->maxgap) r->maxgap = pr[i + 1] - pr[i];
    for (size_t k = 1; k + 1 < n; k++) {
        u64 P = pr[k];
        if (P < lo || P >= hi) continue;
        /* runs ending at gap k-1 (gaps g_i = pr[i+1] - pr[i]): walk back while the run continues */
        int L[4] = { 1, 1, 1, 1 }, M[4] = { 1, 1, 1, 1 };
        for (size_t i = k - 1; i >= 1 && i + NMAX >= k; i--) {
            u64 a = pr[i] - pr[i - 1], b = pr[i + 1] - pr[i];   /* a = g_(i-1), b = g_i */
            bool cont[4] = { a <= b, a >= b, a < b, a > b };     /* nd, ni, si, sd */
            bool any = false;
            for (int m = 0; m < 4; m++) if (L[m] == (int)(k - i) && cont[m]) { L[m]++; any = true; }
            if (!any) break;
        }
        for (size_t i = k; i + 2 < n && i < k + NMAX; i++) {
            u64 a = pr[i + 1] - pr[i], b = pr[i + 2] - pr[i + 1]; /* a = g_i, b = g_(i+1) */
            bool cont[4] = { b <= a, b >= a, b < a, b > a };     /* ni, nd, sd, si */
            bool any = false;
            for (int m = 0; m < 4; m++) if (M[m] == (int)(i - k + 1) && cont[m]) { M[m]++; any = true; }
            if (!any) break;
        }
        int depth[NMODE] = { L[0] < M[0] ? L[0] : M[0], L[1] < M[1] ? L[1] : M[1],
                             L[2] < M[2] ? L[2] : M[2], L[3] < M[3] ? L[3] : M[3] };
        for (int m = 0; m < NMODE; m++)
            for (int d = 1; d <= depth[m] && d <= NMAX; d++) {
                r->count[m][d]++;
                if (!r->first[m][d]) r->first[m][d] = P;
            }
    }
    primesieve_free(pr);
}

/* ------------------------------------------------------------------ */
/* Windows of consecutive primes around P                              */
/* ------------------------------------------------------------------ */

/* w[0..nb'+na] = the nb primes before P (fewer if P is near 2), P, the na primes after;
 * returns the number nb' of primes found before P, or -1 if P is not prime */
static int prime_window(u64 P, int nb, int na, u64 *w)
{
    primesieve_iterator it;
    primesieve_init(&it);
    u64 *ring = xcalloc((size_t)nb + 1, sizeof(u64));
    u64 span = (u64)(nb + na + 2) * 2000;
    int got;
    for (;;) {
        u64 start = P > span ? P - span : 0;
        primesieve_jump_to(&it, start, P + span);
        u64 q;
        int cnt = 0;
        while ((q = primesieve_next_prime(&it)) < P) { ring[cnt % (nb + 1)] = q; cnt++; }
        if (it.is_error) die("primesieve error near %" PRIu64, P);
        if (q != P) { got = -1; break; }
        if (cnt < nb && start > 0) { span *= 4; continue; }
        int have = cnt < nb ? cnt : nb;
        for (int i = 0; i < have; i++) w[i] = ring[(cnt - have + i) % (nb + 1)];
        w[have] = P;
        for (int i = 0; i < na; i++) {
            q = primesieve_next_prime(&it);
            if (it.is_error) die("primesieve error near %" PRIu64, P);
            w[have + 1 + i] = q;
        }
        got = have;
        break;
    }
    free(ring);
    primesieve_free_iterator(&it);
    return got;
}

/* depths of the prime P in the four modes (0 for P = 2); returns false if P is not prime */
static bool depths_at(u64 P, int depth[NMODE])
{
    enum { W = RUNCAP + 2 };
    u64 *w = xcalloc(2 * W + 1, sizeof(u64));
    int nb = prime_window(P, W, W, w);
    if (nb < 0) { free(w); return false; }
    u64 bg[W], ag[W];                           /* bg[0] is the gap into P, ag[0] the gap out of it */
    for (int i = 0; i < nb; i++) bg[i] = w[nb - i] - w[nb - i - 1];
    for (int i = 0; i < W; i++) ag[i] = w[nb + 1 + i] - w[nb + i];
    int L[4] = { nb > 0, nb > 0, nb > 0, nb > 0 }, M[4] = { 1, 1, 1, 1 };
    for (int i = 0; i + 1 < nb; i++) {          /* going back: bg[i+1] is the earlier gap */
        if (L[0] == i + 1 && bg[i + 1] <= bg[i]) L[0]++;   /* non-decreasing into P */
        if (L[1] == i + 1 && bg[i + 1] >= bg[i]) L[1]++;   /* non-increasing into P */
        if (L[2] == i + 1 && bg[i + 1] <  bg[i]) L[2]++;   /* strictly increasing into P */
        if (L[3] == i + 1 && bg[i + 1] >  bg[i]) L[3]++;   /* strictly decreasing into P */
    }
    for (int i = 0; i + 1 < W; i++) {
        if (M[0] == i + 1 && ag[i + 1] <= ag[i]) M[0]++;   /* non-increasing out of P */
        if (M[1] == i + 1 && ag[i + 1] >= ag[i]) M[1]++;   /* non-decreasing out of P */
        if (M[2] == i + 1 && ag[i + 1] <  ag[i]) M[2]++;   /* strictly decreasing out of P */
        if (M[3] == i + 1 && ag[i + 1] >  ag[i]) M[3]++;   /* strictly increasing out of P */
    }
    for (int m = 0; m < NMODE; m++) depth[m] = L[m] < M[m] ? L[m] : M[m];
    free(w);
    return true;
}

/* "gaps [10] 2 4 6 | 8 6 6 [4]" with n gaps on each side and the next gap outside in brackets */
static void print_gap_window(u64 P, int n)
{
    if (n > NMAX) n = NMAX;
    u64 *w = xcalloc(2 * (size_t)n + 5, sizeof(u64));
    int nb = prime_window(P, n + 1, n + 1, w);
    if (nb < 0) { printf("  (%" PRIu64 " is not prime)", P); free(w); return; }
    printf("  primes %" PRIu64 " .. %" PRIu64 " .. %" PRIu64 "  gaps", w[0], P, w[nb + n + 1]);
    for (int i = 0; i < nb; i++) {
        u64 g = w[i + 1] - w[i];
        printf(i == 0 && nb == n + 1 ? " [%" PRIu64 "]" : " %" PRIu64, g);
    }
    printf(" |");
    for (int i = 0; i <= n; i++) {
        u64 g = w[nb + 1 + i] - w[nb + i];
        printf(i == n ? " [%" PRIu64 "]" : " %" PRIu64, g);
    }
    free(w);
}

/* ------------------------------------------------------------------ */
/* Scan state shared by the worker threads                             */
/* ------------------------------------------------------------------ */

typedef struct { int mode, d; u64 P; } rec_t;

typedef struct scan scan_t;

typedef struct {
    int id;
    pthread_t th;
    scan_t *sc;
    primesieve_iterator it;
    rings_t R;
    result_t r;
} worker_t;

struct scan {
    /* parameters */
    u64  origin, start, end, chunk, nchunks;
    unsigned dmin;
    int  report, target;
    bool verbose, print_records;
    /* runtime */
    atomic_uint_fast64_t next;
    atomic_int finished;
    pthread_mutex_t mu;
    pthread_cond_t cv;
    u64  cf;                        /* completion frontier: chunks [0, cf) are folded */
    u64  stop_chunk;                /* no chunk above this index is started (-n) */
    u8   done[RING];
    result_t res[RING];
    /* folded totals: everything in [origin, start + cf*chunk) */
    u64  first[NMODE][NMAX + 1];
    u64  count[NMODE][NMAX + 1];
    u64  primes, maxgap, chunks_done;
    double t0, t1, elapsed0;        /* elapsed0: time spent by earlier runs (checkpoint) */
};

static void print_record(int mode, int d, u64 P)
{
    pthread_mutex_lock(&out_mu);
    if (stderr_tty) fputs("\r\033[K", stderr);
    printf("%s(%d) = %" PRIu64 "  [%s depth >= %d]", SEQ[mode], d - SHIFT[mode], P, MODE[mode], d);
    print_gap_window(P, d);
    printf("\n");
    fflush(stdout);
    pthread_mutex_unlock(&out_mu);
}

/* under sc->mu, in increasing chunk order */
static void fold(scan_t *sc, const result_t *r, rec_t *nr, int *nn)
{
    for (int m = 0; m < NMODE; m++)
        for (int d = 1; d <= NMAX; d++) {
            sc->count[m][d] += r->count[m][d];
            if (r->first[m][d] && !sc->first[m][d]) {
                sc->first[m][d] = r->first[m][d];
                if (d >= sc->report) { nr[*nn].mode = m; nr[*nn].d = d; nr[*nn].P = r->first[m][d]; (*nn)++; }
            }
        }
    sc->primes += r->nprimes;
    if (r->maxgap > sc->maxgap) sc->maxgap = r->maxgap;
}

static void *worker_main(void *arg)
{
    worker_t *w = arg;
    scan_t *sc = w->sc;

    for (;;) {
        if (g_stop) break;
        u64 i = atomic_fetch_add(&sc->next, 1);
        if (i >= sc->nchunks) break;

        pthread_mutex_lock(&sc->mu);
        while (i >= sc->cf + RING && !g_stop) pthread_cond_wait(&sc->cv, &sc->mu);
        bool skip = g_stop || i > sc->stop_chunk;
        pthread_mutex_unlock(&sc->mu);
        if (skip) break;

        u64 lo = sc->start + i * sc->chunk;
        u64 hi = (i + 1 == sc->nchunks) ? sc->end + 1 : lo + sc->chunk;
        unsigned dmin = lo < SMALL ? 1 : sc->dmin;
        scan_chunk(&w->it, &w->R, lo, hi, dmin, &w->r);

        rec_t newrec[NMODE * (NMAX + 1)];
        int nnew = 0;
        pthread_mutex_lock(&sc->mu);
        if (sc->target && w->r.first[PEAK][sc->target] && i < sc->stop_chunk) sc->stop_chunk = i;
        sc->res[i % RING] = w->r;
        sc->done[i % RING] = 1;
        while (sc->cf < sc->nchunks && sc->done[sc->cf % RING]) {
            fold(sc, &sc->res[sc->cf % RING], newrec, &nnew);
            sc->done[sc->cf % RING] = 0;
            sc->cf++;
        }
        sc->chunks_done++;
        pthread_cond_broadcast(&sc->cv);
        pthread_mutex_unlock(&sc->mu);
        if (sc->print_records)
            for (int k = 0; k < nnew; k++) print_record(newrec[k].mode, newrec[k].d, newrec[k].P);
    }
    atomic_fetch_add(&sc->finished, 1);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Checkpoint files                                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    unsigned dmin;
    u64  origin, end, lo, primes, maxgap;
    double elapsed;
    u64  first[NMODE][NMAX + 1];
    u64  count[NMODE][NMAX + 1];
} state_t;

static void write_state(const char *path, const state_t *st)
{
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { note("warning: cannot write %s: %s", tmp, strerror(errno)); return; }
    fprintf(f, "A248701-state 1\ndmin %u\norigin %" PRIu64 "\nend %" PRIu64 "\nlo %" PRIu64 "\nprimes %" PRIu64
               "\nmaxgap %" PRIu64 "\nelapsed %.1f\n", st->dmin, st->origin, st->end, st->lo, st->primes, st->maxgap, st->elapsed);
    for (int m = 0; m < NMODE; m++)
        for (int d = 1; d <= NMAX; d++) {
            if (st->first[m][d]) fprintf(f, "first %d %d %" PRIu64 "\n", m, d, st->first[m][d]);
            if (st->count[m][d]) fprintf(f, "count %d %d %" PRIu64 "\n", m, d, st->count[m][d]);
        }
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
    if (!fgets(line, sizeof line, f) || sscanf(line, "A248701-state %d", &ver) != 1 || ver != 1)
        die("%s is not an a248701 state file", path);
    while (fgets(line, sizeof line, f)) {
        char key[32], v1[64], v2[64], v3[64];
        int n = sscanf(line, "%31s %63s %63s %63s", key, v1, v2, v3);
        u128 v;
        if (n < 2) continue;
        if (!strcmp(key, "dmin")) st->dmin = (unsigned)atoi(v1);
        else if (!strcmp(key, "origin")) { if (!parse_u128(v1, &v)) die("bad state"); st->origin = (u64)v; }
        else if (!strcmp(key, "end"))    { if (!parse_u128(v1, &v)) die("bad state"); st->end = (u64)v; }
        else if (!strcmp(key, "lo"))     { if (!parse_u128(v1, &v)) die("bad state"); st->lo = (u64)v; }
        else if (!strcmp(key, "primes")) { if (!parse_u128(v1, &v)) die("bad state"); st->primes = (u64)v; }
        else if (!strcmp(key, "maxgap")) { if (!parse_u128(v1, &v)) die("bad state"); st->maxgap = (u64)v; }
        else if (!strcmp(key, "elapsed")) st->elapsed = atof(v1);
        else if ((!strcmp(key, "first") || !strcmp(key, "count")) && n >= 4) {
            int m = atoi(v1), d = atoi(v2);
            if (m < 0 || m >= NMODE || d < 1 || d > NMAX || !parse_u128(v3, &v)) die("bad %s line in %s", key, path);
            if (key[0] == 'f') st->first[m][d] = (u64)v;
            else st->count[m][d] = (u64)v;
        }
    }
    fclose(f);
    if (st->dmin < 1 || st->dmin > 8) die("bad dmin in %s", path);
    return true;
}

static void checkpoint(scan_t *sc, const char *path)
{
    state_t st;
    memset(&st, 0, sizeof st);
    pthread_mutex_lock(&sc->mu);
    st.dmin = sc->dmin;
    st.origin = sc->origin;
    st.end = sc->end;
    st.lo = sc->cf >= sc->nchunks ? sc->end + 1 : sc->start + sc->cf * sc->chunk;
    st.primes = sc->primes;
    st.maxgap = sc->maxgap;
    st.elapsed = sc->elapsed0 + now() - sc->t0;
    memcpy(st.first, sc->first, sizeof st.first);
    memcpy(st.count, sc->count, sizeof st.count);
    pthread_mutex_unlock(&sc->mu);
    write_state(path, &st);
}

/* ------------------------------------------------------------------ */
/* Running a scan                                                      */
/* ------------------------------------------------------------------ */

static u64 scanned_to(const scan_t *sc)
{
    return sc->cf >= sc->nchunks ? sc->end : sc->start + sc->cf * sc->chunk - 1;
}

static void print_status(scan_t *sc, double t, bool final)
{
    pthread_mutex_lock(&sc->mu);
    u64 pos = scanned_to(sc) + 1;
    u64 primes = sc->primes, maxgap = sc->maxgap;
    int best[NMODE] = { 0, 0, 0, 0 };
    u64 bestp = 0;
    for (int m = 0; m < NMODE; m++)
        for (int d = NMAX; d >= 1; d--) if (sc->first[m][d]) { best[m] = d; if (m == PEAK) bestp = sc->first[m][d]; break; }
    pthread_mutex_unlock(&sc->mu);

    double done = (double)(pos - sc->start), total = (double)(sc->end - sc->start) + 1.0;
    double el = t - sc->t0, rate = el > 0 ? done / el : 0, eta = rate > 0 ? (total - done) / rate : 0;
    char b1[32], b2[32], b3[32], b4[32], b5[32];
    pthread_mutex_lock(&out_mu);
    if (stderr_tty) fputs("\r\033[K", stderr);
    fprintf(stderr, "[%s] %5.1f%% at %s  %s/s  ETA %s  primes %s  maxgap %" PRIu64 "  depth peak %d @ %" PRIu64 " valley %d strict %d/%d",
            fmt_hms(el, b1), 100.0 * done / total, fmt_eng((double)pos, b2), fmt_eng(rate, b3), fmt_hms(eta, b4),
            fmt_eng((double)primes, b5), maxgap, best[PEAK], bestp, best[VALLEY], best[SPEAK], best[SVALLEY]);
    if (!stderr_tty || final) fputc('\n', stderr);
    fflush(stderr);
    pthread_mutex_unlock(&out_mu);
}

static u64 default_chunk(u64 range, int threads)
{
    u64 c = range / ((u64)threads * 64);
    if (c < (1ull << 20)) c = 1ull << 20;
    if (c > (1ull << 34)) c = 1ull << 34;
    return c;
}

static void scan_run(scan_t *sc, int nthreads, const char *statefile, double ckpt_secs)
{
    if (sc->end > primesieve_get_max_stop() - (1ull << 33)) die("END too large for primesieve");
    if (nthreads < 1) nthreads = 1;
    u64 range = sc->end >= sc->start ? sc->end - sc->start + 1 : 0;
    if (!sc->chunk) sc->chunk = default_chunk(range, nthreads);
    sc->nchunks = range ? (range + sc->chunk - 1) / sc->chunk : 0;
    if (sc->dmin < 1) sc->dmin = 1;

    atomic_store(&sc->next, 0);
    atomic_store(&sc->finished, 0);
    pthread_mutex_init(&sc->mu, NULL);
    pthread_cond_init(&sc->cv, NULL);
    sc->cf = 0;
    sc->stop_chunk = UINT64_MAX;
    memset(sc->done, 0, sizeof sc->done);
    sc->chunks_done = 0;
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
        if (g_stop) { pthread_mutex_lock(&sc->mu); pthread_cond_broadcast(&sc->cv); pthread_mutex_unlock(&sc->mu); }
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
    }
    free(ws);
    sc->t1 = now();
    if (sc->verbose) print_status(sc, sc->t1, true);
    if (statefile) checkpoint(sc, statefile);
}

static scan_t *scan_new(void)
{
    scan_t *sc = xcalloc(1, sizeof *sc);
    sc->dmin = 4;
    sc->report = 6;
    return sc;
}

/* ------------------------------------------------------------------ */
/* Commands                                                            */
/* ------------------------------------------------------------------ */

static void usage(void)
{
    fputs(
"usage: a248701 scan [START] END [-t T] [-c CHUNK] [-n N] [-d DMIN] [-r DEPTH] [-S STATE] [-i SECS] [-q]\n"
"       a248701 show P [-n N]\n"
"       a248701 selftest [-t T]\n"
"\n"
"  scan     exhaustive search of the primes in [START, END] (START defaults to 0) for the\n"
"           smallest prime of every depth in the four modes: A248701 (peak), A248702\n"
"           (valley), A248703 (strict peak), A248704 (strict valley).  Prints each new\n"
"           record of depth >= DEPTH (default 6) as it is found, then the tables.\n"
"           -n N stops as soon as A248701(N) is known.  -d DMIN (default 4, max 8) is the\n"
"           smallest depth examined above 2^30; all depths are examined below 2^30, which\n"
"           holds every known term.  -S FILE writes a checkpoint every SECS seconds\n"
"           (default 60) and on Ctrl-C; running the same command again resumes from it.\n"
"  show     the depths of the prime P in the four modes and the gaps around it\n"
"           (-n N: show N gaps on each side).\n"
"  numbers  may be written as 123, 2^40, 10^12, 1e12, 2^40-1\n", stderr);
    exit(2);
}

static void on_signal(int s)
{
    (void)s;
    g_stop = 1;
}

static void print_summary(scan_t *sc)
{
    char b1[32];
    double secs = sc->t1 - sc->t0;
    u64 to = scanned_to(sc);
    bool complete = sc->cf >= sc->nchunks;
    printf("# scanned [%" PRIu64 ", %" PRIu64 "]%s in %.1f s (%s numbers/s%s): %" PRIu64 " primes and largest gap %" PRIu64
           " in [%" PRIu64 ", %" PRIu64 "]\n",
           sc->start, to, complete ? "" : g_stop ? " (interrupted)" : " (stopped early, -n)", secs,
           fmt_eng((double)(to + 1 - sc->start) / (secs > 0 ? secs : 1), b1),
           sc->elapsed0 > 0 ? ", this run" : "", sc->primes, sc->maxgap, sc->origin, to);
    if (sc->elapsed0 > 0) printf("# total time including earlier runs %.1f s\n", sc->elapsed0 + secs);

    int maxd = 0;
    for (int m = 0; m < NMODE; m++)
        for (int d = NMAX; d >= 1; d--) if (sc->count[m][d]) { if (d > maxd) maxd = d; break; }
    unsigned dlo = sc->dmin;
    printf("# centres of depth >= d in [%" PRIu64 ", %" PRIu64 "]%s:\n", sc->origin, to,
           dlo > 1 ? " (depths below -d are only examined in chunks starting below 2^30)" : "");
    printf("#   d            :");
    for (int d = (int)dlo; d <= maxd; d++) printf(" %12d", d);
    printf("\n");
    for (int m = 0; m < NMODE; m++) {
        printf("#   %-13s:", MODE[m]);
        for (int d = (int)dlo; d <= maxd; d++) printf(" %12" PRIu64, sc->count[m][d]);
        printf("\n");
    }

    for (int m = 0; m < NMODE; m++) {
        printf("%s (%s depth, a(n) = depth >= %s):", SEQ[m], MODE[m], SHIFT[m] ? "n+1" : "n");
        if (m == VALLEY) printf("\n  a(0) = 2   (by definition)");
        for (int d = 1; d <= NMAX; d++) {
            int n = d - SHIFT[m];
            if (n < FIRSTN[m]) continue;
            u64 a = sc->first[m][d];
            if (!a && (unsigned)d < sc->dmin && sc->origin >= SMALL) {
                printf("\n  a(%d): depth %d not examined (-d %u)", n, d, sc->dmin);
                continue;
            }
            if (a) {
                const char *tag = "NEW";
                if (d <= KDEPTH && KNOWN[m][d]) tag = a == KNOWN[m][d] ? "(matches OEIS)" : "** DIFFERS FROM OEIS **";
                printf("\n  a(%d) = %" PRIu64 "   %s", n, a, tag);
            } else {
                if (sc->origin == 0) printf("\n  a(%d) > %" PRIu64, n, to);
                else printf("\n  a(%d): no prime of depth >= %d in [%" PRIu64 ", %" PRIu64 "]", n, d, sc->origin, to);
                break;
            }
        }
        printf("\n");
    }
    fflush(stdout);
}

static int cmd_scan(int argc, char **argv)
{
    u64 pos[2] = { 0, 0 };
    int npos = 0, threads = 0, target = 0, report = 6, dmin = 4;
    u64 chunk = 0;
    const char *statefile = NULL;
    double ckpt = 60;

    for (int i = 0; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-t") && i + 1 < argc) threads = arg_int(argv[++i], "thread count", 1, 1024);
        else if (!strcmp(a, "-c") && i + 1 < argc) chunk = arg_u64(argv[++i], "chunk size");
        else if (!strcmp(a, "-n") && i + 1 < argc) target = arg_int(argv[++i], "target depth", 1, NMAX);
        else if (!strcmp(a, "-d") && i + 1 < argc) dmin = arg_int(argv[++i], "DMIN", 1, 8);
        else if (!strcmp(a, "-r") && i + 1 < argc) report = arg_int(argv[++i], "report depth", 1, NMAX);
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
    if (end < start) die("END is below START");
    scan_t *sc = scan_new();
    sc->dmin = (unsigned)dmin;
    sc->report = report;
    sc->target = target;
    bool resumed = false;
    state_t st;

    if (statefile && read_state(statefile, &st)) {
        if (st.dmin != (unsigned)dmin) die("%s was written with -d %u; use the same value", statefile, st.dmin);
        if (npos == 2 && pos[0] != st.lo)
            note("note: START %" PRIu64 " ignored, resuming from checkpoint at %" PRIu64, pos[0], st.lo);
        if (end < st.lo) {
            if (end == st.end) { note("%s: scan of [%" PRIu64 ", %" PRIu64 "] already complete", statefile, st.origin, st.end); return 0; }
            die("END %" PRIu64 " is below the checkpoint position %" PRIu64, end, st.lo);
        }
        start = st.lo;
        origin = st.origin;
        memcpy(sc->first, st.first, sizeof sc->first);
        memcpy(sc->count, st.count, sizeof sc->count);
        sc->primes = st.primes;
        sc->maxgap = st.maxgap;
        sc->elapsed0 = st.elapsed;
        resumed = true;
        note("resuming from %s at %" PRIu64 " (%" PRIu64 " primes done, %.0f s spent)", statefile, start, st.primes, st.elapsed);
    }

    if (!chunk) chunk = default_chunk(end - start + 1, threads);
    sc->origin = origin;
    sc->start = start;
    sc->end = end;
    sc->chunk = chunk;
    sc->verbose = !g_quiet;
    sc->print_records = true;

    printf("# A248701 / A248702 / A248703 / A248704: smallest prime with n monotone prime gaps on each side (peak, valley, strict peak, strict valley)\n");
    printf("# scanning primes in [%" PRIu64 ", %" PRIu64 "] with %d threads, chunk %" PRIu64 ", -d %d, reporting depth >= %d%s%s\n",
           start, end, threads, chunk, dmin, report, resumed ? " (resumed)" : "", target ? ", stopping at A248701 depth" : "");
    if (target) printf("# stop target: A248701(%d)\n", target);
    fflush(stdout);

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_signal;
    sa.sa_flags = SA_RESETHAND;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    scan_run(sc, threads, statefile, ckpt);
    if (g_stop) note("interrupted; all chunks below %" PRIu64 " are complete%s", scanned_to(sc) + 1,
                     statefile ? ", checkpoint written" : "");
    print_summary(sc);
    free(sc);
    return 0;
}

static int cmd_show(int argc, char **argv)
{
    u64 P = 0;
    int nshow = 0;
    bool have = false;
    for (int i = 0; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-n") && i + 1 < argc) nshow = arg_int(argv[++i], "N", 1, NMAX);
        else if (a[0] == '-' && !isdigit((unsigned char)a[1])) die("unknown option '%s'", a);
        else { P = arg_u64(a, "prime"); have = true; }
    }
    if (!have) usage();
    int depth[NMODE];
    if (!depths_at(P, depth)) die("%" PRIu64 " is not prime", P);
    printf("%" PRIu64 ":", P);
    int maxd = 0;
    for (int m = 0; m < NMODE; m++) {
        printf("%s %s depth %d (%s index %d)", m ? "," : "", MODE[m], depth[m], SEQ[m], depth[m] - SHIFT[m]);
        if (depth[m] > maxd) maxd = depth[m];
    }
    printf("\n");
    int n = nshow ? nshow : (maxd > 0 ? maxd : 1);
    print_gap_window(P, n);
    printf("\n");
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

static scan_t *quiet_scan(u64 start, u64 end, u64 chunk, int threads, unsigned dmin, int target)
{
    scan_t *sc = scan_new();
    sc->origin = start; sc->start = start; sc->end = end; sc->chunk = chunk;
    sc->dmin = dmin; sc->target = target;
    sc->verbose = false; sc->print_records = false; sc->report = NMAX + 1;
    scan_run(sc, threads, NULL, 0);
    return sc;
}

static bool same_results(const result_t *a, const result_t *b, unsigned dlo, bool *why_first)
{
    *why_first = false;
    if (a->nprimes != b->nprimes || a->maxgap != b->maxgap) return false;
    for (int m = 0; m < NMODE; m++)
        for (unsigned d = dlo; d <= NMAX; d++) {
            if (a->first[m][d] != b->first[m][d]) { *why_first = true; return false; }
            if (a->count[m][d] != b->count[m][d]) return false;
        }
    return true;
}

static void totals_of(const scan_t *sc, result_t *r)
{
    memset(r, 0, sizeof *r);
    memcpy(r->first, sc->first, sizeof r->first);
    memcpy(r->count, sc->count, sizeof r->count);
    r->nprimes = sc->primes;
    r->maxgap = sc->maxgap;
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
    double t0 = now();
    printf("a248701 selftest (%d threads)\n", threads);

    /* 1. the OEIS example and the small terms through depths_at */
    {
        int d[NMODE];
        bool ok = depths_at(359, d) && d[PEAK] == 3 && d[SPEAK] == 2 && d[VALLEY] == 1 && d[SVALLEY] == 1;
        expect(ok, "359: peak depth 3, strict peak depth 2, valley depths 1 (gaps 10, 2, 4, 6 | 8, 6, 6, 4)");
        ok = depths_at(7, d) && d[PEAK] == 2 && d[SPEAK] == 1 && d[VALLEY] == 1;
        expect(ok, "7: peak depth 2 but strict peak depth 1 (gaps 2, 2 | 4, 2)");
        ok = depths_at(43, d) && d[VALLEY] == 3 && d[SVALLEY] == 2;
        expect(ok, "43: valley depth 3, strict valley depth 2 (gaps 6, 4, 2 | 4, 6, 6)");
        ok = depths_at(1429, d) && d[SVALLEY] == 3 && depths_at(1439, d) && d[SPEAK] == 3;
        expect(ok, "1429: strict valley depth 3; 1439: strict peak depth 3");
        ok = depths_at(2, d) && d[PEAK] == 0 && d[VALLEY] == 0 && depths_at(3, d) && d[PEAK] == 1 && d[VALLEY] == 1;
        expect(ok, "2 has depth 0 (no gap before it), 3 has depth 1");
        ok = !depths_at(360, d) && !depths_at(1, d);
        expect(ok, "360 and 1 are rejected as non-primes");
    }

    /* 2. the chunk scanner against the reference implementation */
    {
        static const u64 LOS[] = { 0, 1000000000ull, 1000000000000ull, 123456789012345ull };
        primesieve_iterator it;
        primesieve_init(&it);
        rings_t *R = xcalloc(1, sizeof *R);
        result_t ra, rb;
        bool wf;
        for (size_t i = 0; i < sizeof LOS / sizeof LOS[0]; i++) {
            u64 lo = LOS[i], hi = lo + 3000000;
            reference_scan(lo, hi, &ra);
            scan_chunk(&it, R, lo, hi, 1, &rb);
            bool ok = same_results(&ra, &rb, 1, &wf);
            int maxd = 0;
            for (int d = NMAX; d >= 1; d--) if (ra.count[PEAK][d]) { maxd = d; break; }
            expect(ok, "chunk [%" PRIu64 ", %" PRIu64 ") agrees with the reference implementation (%" PRIu64 " primes, max peak depth %d)%s",
                   lo, hi, ra.nprimes, maxd, ok ? "" : wf ? " [records differ]" : " [counts differ]");
        }
        /* a chunk split in two must give the same totals as one piece */
        result_t r1, r2;
        scan_chunk(&it, R, 5000000000ull, 5003000000ull, 1, &ra);
        scan_chunk(&it, R, 5000000000ull, 5001500000ull, 1, &r1);
        scan_chunk(&it, R, 5001500000ull, 5003000000ull, 1, &r2);
        bool ok = ra.nprimes == r1.nprimes + r2.nprimes && ra.maxgap == (r1.maxgap > r2.maxgap ? r1.maxgap : r2.maxgap);
        for (int m = 0; m < NMODE; m++)
            for (int d = 1; d <= NMAX; d++) {
                u64 f = r1.first[m][d] ? r1.first[m][d] : r2.first[m][d];
                if (ra.count[m][d] != r1.count[m][d] + r2.count[m][d] || ra.first[m][d] != f) ok = false;
            }
        expect(ok, "chunk [5e9, 5e9+3e6) equals the sum of its two halves");
        primesieve_free_iterator(&it);
        free(R);
    }

    /* 3. the four sequences below 4e8 and the depth counts (checked against an independent numpy computation) */
    scan_t *sc = quiet_scan(0, 399999999ull, 1000000, threads, 1, 0);
    {
        bool ok = true;
        for (int m = 0; m < NMODE; m++)
            for (int d = 1; d <= KDEPTH; d++) {
                if (!KNOWN[m][d]) continue;
                if (sc->first[m][d] != KNOWN[m][d]) {
                    ok = false;
                    printf("      %s depth %d: %" PRIu64 ", OEIS has %" PRIu64 "\n", SEQ[m], d, sc->first[m][d], KNOWN[m][d]);
                }
            }
        expect(ok, "scan [0, 4e8) reproduces all 26 known terms of A248701-A248704 (%" PRIu64 " primes, %.1f s)",
               sc->primes, sc->t1 - sc->t0);
        static const u64 EXP_PEAK[8]  = { 0, 21336325, 5601745, 670557, 47070, 2154, 74, 2 };
        static const u64 EXP_VALL[8]  = { 0, 21336325, 5597514, 669692, 46218, 2078, 69, 0 };
        static const u64 EXP_SPK[8]   = { 0, 21336325, 5016394, 481047, 24392, 715, 21, 1 };
        static const u64 EXP_SVL[8]   = { 0, 21336325, 5020624, 483720, 24355, 755, 16, 0 };
        ok = true;
        for (int d = 1; d <= 7; d++)
            if (sc->count[PEAK][d] != EXP_PEAK[d] || sc->count[VALLEY][d] != EXP_VALL[d] ||
                sc->count[SPEAK][d] != EXP_SPK[d] || sc->count[SVALLEY][d] != EXP_SVL[d]) ok = false;
        expect(ok, "depth counts below 4e8 match numpy: peak %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64,
               sc->count[PEAK][1], sc->count[PEAK][2], sc->count[PEAK][3], sc->count[PEAK][4], sc->count[PEAK][5], sc->count[PEAK][6], sc->count[PEAK][7]);
        expect(sc->count[PEAK][8] == 0 && sc->count[VALLEY][8] == 0, "no centre of depth >= 8 below 4e8");
    }

    /* 4. chunking independence and a scan that starts in the middle */
    {
        result_t ta, tb, tc, td, te;
        bool wf;
        scan_t *s1 = quiet_scan(0, 299999999ull, 300000000ull, 1, 1, 0);
        scan_t *s2 = quiet_scan(0, 299999999ull, 777777, threads, 1, 0);
        scan_t *s3 = quiet_scan(0, 299999999ull, 1 << 20, threads, 1, 0);
        totals_of(s1, &ta); totals_of(s2, &tb); totals_of(s3, &tc);
        bool ok = same_results(&ta, &tb, 1, &wf) && same_results(&ta, &tc, 1, &wf);
        expect(ok, "scan [0, 3e8) gives identical records and counts as one chunk, chunk 777777 and chunk 2^20");

        scan_t *s4 = quiet_scan(300000000ull, 399999999ull, 1 << 20, threads, 1, 0);
        totals_of(s4, &td);
        totals_of(sc, &te);
        ok = te.nprimes == ta.nprimes + td.nprimes;
        for (int m = 0; m < NMODE && ok; m++)
            for (int d = 1; d <= NMAX; d++) {
                if (te.count[m][d] != ta.count[m][d] + td.count[m][d]) ok = false;
                u64 f = ta.first[m][d] ? ta.first[m][d] : td.first[m][d];
                if (te.first[m][d] != f) ok = false;
            }
        expect(ok && td.first[SPEAK][7] == 331542583ull, "scan [3e8, 4e8) adds up with [0, 3e8) and finds A248703(6) = 331542583 above 3e8");
        free(s1); free(s2); free(s3); free(s4);
    }

    /* 5. -d 4 above 2^30 only drops depths below 4 */
    {
        result_t ta, tb;
        bool wf;
        scan_t *s1 = quiet_scan(1ull << 30, (1ull << 30) + 199999999ull, 1 << 22, threads, 1, 0);
        scan_t *s2 = quiet_scan(1ull << 30, (1ull << 30) + 199999999ull, 1 << 22, threads, 4, 0);
        totals_of(s1, &ta); totals_of(s2, &tb);
        bool ok = same_results(&ta, &tb, 4, &wf) && tb.count[PEAK][3] == 0 && ta.count[PEAK][3] > 0;
        expect(ok, "scan [2^30, 2^30+2e8) with -d 1 and -d 4 agree on every depth >= 4 (%" PRIu64 " centres of depth >= 4, %" PRIu64 " of depth >= 5)",
               ta.count[PEAK][4], ta.count[PEAK][5]);
        free(s1); free(s2);
    }

    /* 6. early stop */
    {
        scan_t *s1 = quiet_scan(0, 999999999ull, 1000000, threads, 1, 7);
        bool ok = s1->first[PEAK][7] == 294276293ull && s1->cf < s1->nchunks && scanned_to(s1) >= 294276293ull &&
                  scanned_to(s1) < 294276293ull + (u64)(threads + RING) * 1000000ull;
        expect(ok, "-n 7 on [0, 1e9) stops soon after a(7) = 294276293 (scanned to %" PRIu64 ")", scanned_to(s1));
        free(s1);
    }

    /* 7. checkpoint round trip */
    {
        state_t st, st2;
        memset(&st, 0, sizeof st);
        st.dmin = 4; st.origin = 0; st.end = 12345678901ull; st.lo = 4000000000ull; st.primes = 189961812; st.maxgap = 354; st.elapsed = 123.5;
        st.first[PEAK][7] = 294276293; st.first[SPEAK][7] = 331542583; st.count[PEAK][4] = 47070; st.count[SVALLEY][7] = 1;
        const char *path = "/tmp/a248701-selftest.state";
        write_state(path, &st);
        bool ok = read_state(path, &st2) && st2.dmin == 4 && st2.end == st.end && st2.lo == st.lo && st2.primes == st.primes &&
                  st2.maxgap == 354 && st2.elapsed == 123.5 && st2.first[PEAK][7] == 294276293 && st2.first[SPEAK][7] == 331542583 &&
                  st2.first[PEAK][6] == 0 && st2.count[PEAK][4] == 47070 && st2.count[SVALLEY][7] == 1 && st2.count[PEAK][5] == 0;
        unlink(path);
        expect(ok, "checkpoint file round trip");
    }

    free(sc);
    printf("%s (%.1f s)\n", st_fail ? "SELFTEST FAILED" : "all tests passed", now() - t0);
    return st_fail ? 1 : 0;
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    setvbuf(stdout, NULL, _IOLBF, 0);
    stderr_tty = isatty(2);
    if (argc < 2) usage();
    const char *cmd = argv[1];
    if (!strcmp(cmd, "scan")) return cmd_scan(argc - 2, argv + 2);
    if (!strcmp(cmd, "show") || !strcmp(cmd, "depth")) return cmd_show(argc - 2, argv + 2);
    if (!strcmp(cmd, "selftest")) return cmd_selftest(argc - 2, argv + 2);
    usage();
    return 2;
}
