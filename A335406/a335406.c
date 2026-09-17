/*
 * a335406.c
 *
 * Compute and try to extend OEIS A335406:
 *
 *   "First position of n in the sequence of run-lengths of the sequence of
 *    prime gaps."
 *
 * Known terms (n = 1..5):  1, 2, 49, 633353, 6706139.
 *
 * Definitions
 * -----------
 * Let p_1 = 2 < p_2 = 3 < ... be the primes and g_i = p_(i+1) - p_i the prime
 * gaps (A001223).  Cut the gap sequence into maximal runs of equal values,
 * (1), (2,2), (4), (2), (4), (2), (4), (6), (2), (6), (4), (2), (4), (6,6), ...
 * The run lengths form A333254 = 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, ...
 * and a(n) is the index of the first run of length exactly n.
 *
 * A run of n equal gaps d beginning at gap g_i is the same thing as the n+1
 * consecutive primes p_i, p_i + d, ..., p_i + n*d in arithmetic progression
 * (a "CPAP-(n+1)").  So the first run of length n begins at the smallest
 * CPAP-(n+1), which is A006560(n+1); its gap index i is A089180(n-1); and
 * because a new run begins exactly at the gaps that differ from their
 * predecessor,
 *
 *     a(n) = i - E(i-1),     E(m) = #{ 1 <= j <= m : g_j = g_(j+1) }.
 *
 * Example: 251, 257, 263, 269 are the first four consecutive primes in
 * arithmetic progression.  251 = p_54, and among g_1..g_53 exactly five gaps
 * equal their successor (j = 2, 15, 36, 39, 46), so a(3) = 54 - 5 = 49.
 *
 * Hence a(6) is the run index of the smallest CPAP-7.  Its common difference
 * is a multiple of 210; the smallest known CPAP-7 starts at
 * 71137654873189893604531 ~ 7.1e22 (P. Zimmermann, 2022) and is not known to
 * be the smallest -- the calibrated model in model_a335406.py expects the
 * first one between 1e20 and 8e20 (median 3.5e20), i.e. a(6) ~ 2e18..2e19.
 * Evaluating a(6) needs E at that point, i.e. every prime gap up to ~1e20,
 * which is centuries of sieving on a workstation (see the README).  What this
 * tool does: reproduce a(1..5), measure the run statistics as far as one
 * cares to sieve, and give the rigorous lower bound a(n) > (number of runs
 * below N) for the lengths n not yet seen; a(6) > 339595430914 from N = 1e13.
 *
 * Method
 * ------
 * The range is cut into chunks [lo, hi).  Worker threads take chunks from an
 * atomic counter and walk the primes of a chunk with a primesieve iterator,
 * starting from the last prime below lo so that the gap into the first prime
 * is known.  A chunk owns the runs that begin at one of its primes.  It counts
 * its primes, the primes whose gap equals the previous gap (the pairs counted
 * by E), a histogram of run lengths (with -H also of (length, gap) pairs), and
 * remembers the first run of each length that begins in the chunk as a pair of
 * offsets (primes and equal pairs of the chunk before it).  A run that begins
 * in a chunk is followed past hi to its end; a run that begins before lo is
 * left to the previous chunk; so every run is counted exactly once.  Finished
 * chunks are folded into the global totals in order, which turns the offsets
 * into absolute indices and gives a(n) = i - E(i-1).  The fold also prints
 * the cumulative statistics at the requested marks (-m).
 *
 * Everything is exact integer counting; no primality testing is involved.
 *
 * Usage
 * -----
 *   a335406 scan [START] END [-t T] [-c CHUNK] [-m MARKS|pow10] [-H] [-q]
 *       Scan the primes in [START, END] (START defaults to 0).  Prints each
 *       a(n) as soon as its run has been folded in, then the totals: primes,
 *       equal adjacent gaps, runs, runs by length, and the lower bound
 *       a(n) > runs for the lengths not found.  With START > 2 only the range
 *       statistics are meaningful and no run indices are printed.
 *       -m M1,M2,...  print the cumulative totals for the primes <= M_k
 *                     (chunk boundaries are aligned to the marks);
 *       -m pow10      the same at every power of ten inside the range;
 *       -H            also print the (length, gap) table as lines
 *                     "ld POS LEN GAP COUNT".
 *   a335406 runs N
 *       Print the first N terms of A333254 (the run lengths), single pass.
 *   a335406 naive END [-H]
 *       Single-threaded single-pass reference implementation of 'scan'.
 *   a335406 selftest [-t T]
 *       Check a(1..5), A089180, A006560, the first 60 terms of A333254, and
 *       the agreement of the chunked scan with the reference implementation
 *       for several chunk sizes, thread counts, marks and split ranges.
 *
 * Numbers may be written as decimal, 2^k, 10^k, 1e12, or X+Y / X-Y of those.
 *
 * Build:  cc -O2 -std=gnu11 -pthread -I/opt/homebrew/include a335406.c \
 *            -L/opt/homebrew/lib -lprimesieve -lm -o a335406
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
#include <time.h>
#include <unistd.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdatomic.h>
#include <primesieve.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;

#define LMAX 16                       /* run lengths 1..LMAX are tracked separately, longer ones land in LMAX */
#define DMAX 4096                     /* gaps 1..DMAX-1 in the (length, gap) table, larger ones land in DMAX-1 */
#define RING 256                      /* finished chunks that may wait for the ordered fold */
#define NLD ((size_t)(LMAX + 1) * DMAX)

/* OEIS data used by the selftest */
static const u64 KNOWN_A[6] = {0, 1, 2, 49, 633353, 6706139};        /* A335406(n) */
static const u64 KNOWN_I[6] = {0, 1, 2, 54, 654926, 6904737};        /* gap index i = A089180(n-1), n >= 2 */
static const u64 KNOWN_P[6] = {0, 2, 3, 251, 9843019, 121174811};    /* first prime = A006560(n+1) */
static const u64 KNOWN_D[6] = {0, 1, 2, 6, 30, 30};                  /* common gap */
static const char *A333254_DATA =
    "1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, "
    "1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1";

/* ------------------------------------------------------------------ */
/* Utilities                                                           */
/* ------------------------------------------------------------------ */

static volatile sig_atomic_t g_stop = 0;
static void on_sigint(int s) { (void)s; g_stop = 1; }

static pthread_mutex_t out_mu = PTHREAD_MUTEX_INITIALIZER;
static bool stderr_tty = false;

static void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fputs("a335406: ", stderr);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
    exit(1);
}

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
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
    if (s < 0 || s > 1e11) { strcpy(buf, "--:--:--"); return buf; }
    long t = (long)(s + 0.5);
    sprintf(buf, "%02ld:%02ld:%02ld", t / 3600, (t / 60) % 60, t % 60);
    return buf;
}

/* 1234567 -> "1,234,567" */
static const char *fmt_commas(u64 v, char *buf)
{
    char tmp[32];
    int n = sprintf(tmp, "%" PRIu64, v), o = 0;
    for (int i = 0; i < n; i++) {
        if (i && (n - i) % 3 == 0) buf[o++] = ',';
        buf[o++] = tmp[i];
    }
    buf[o] = 0;
    return buf;
}

/* ------------------------------------------------------------------ */
/* One chunk of the range                                              */
/* ------------------------------------------------------------------ */

typedef struct {
    bool set;
    u64 off;        /* primes of the chunk before the start prime */
    u64 eqb;        /* equal pairs of the chunk before the start prime */
    u64 p;          /* start prime */
    u64 d;          /* common gap */
} first_t;

typedef struct {
    u64 lo, hi;
    u64 nprimes;                /* primes in [lo, hi) */
    u64 neq;                    /* primes p_k in [lo, hi) with g_(k-1) = g_k */
    u64 hist[LMAX + 1];         /* runs starting in [lo, hi), by length (>= LMAX in LMAX) */
    unsigned maxlen;            /* longest run starting in [lo, hi) */
    first_t first[LMAX + 1];    /* first run of each exact length starting in the chunk */
    u64 *ld;                    /* (length, gap) counts, NLD entries, or NULL */
} chunk_t;

static void finish_run(chunk_t *c, unsigned len, const first_t *st)
{
    unsigned L = len > LMAX ? LMAX : len;
    c->hist[L]++;
    if (len > c->maxlen) c->maxlen = len;
    if (len <= LMAX && !c->first[len].set) {
        c->first[len] = *st;
        c->first[len].set = true;
    }
    if (c->ld) c->ld[(size_t)L * DMAX + (st->d < DMAX ? st->d : DMAX - 1)]++;
}

/* largest prime < lo, or 0 if there is none */
static u64 prime_below(primesieve_iterator *it, u64 lo)
{
    if (lo <= 2) return 0;
    primesieve_jump_to(it, lo - 1, lo);
    u64 p = primesieve_prev_prime(it);
    if (it->is_error) die("primesieve error below %" PRIu64, lo);
    return p;
}

/*
 * Walk the primes p_k in [lo, hi).  For each one the gap g_k into the next
 * prime is compared with the gap g_(k-1) into p_k: equal means p_k continues
 * the current run (and is one of the pairs counted by E), different means a
 * new run begins at p_k.  Runs that begin in the chunk are the chunk's own and
 * are followed to their end, past hi if necessary.
 */
static void process_chunk(primesieve_iterator *it, chunk_t *c, u64 lo, u64 hi)
{
    c->lo = lo;
    c->hi = hi;
    c->nprimes = c->neq = 0;
    c->maxlen = 0;
    memset(c->hist, 0, sizeof c->hist);
    memset(c->first, 0, sizeof c->first);
    if (c->ld) memset(c->ld, 0, NLD * sizeof(u64));

    u64 pprev = prime_below(it, lo);
    primesieve_jump_to(it, lo, hi);
    u64 p = primesieve_next_prime(it);      /* first prime >= lo */
    u64 q = primesieve_next_prime(it);      /* the prime after it */
    u64 gprev = pprev ? p - pprev : 0;      /* gap into p; 0 only for p = 2 */
    bool ours = false;                      /* does the current run begin in this chunk? */
    unsigned len = 0;
    first_t st = {0};
    u64 idx = 0, eq = 0;

    while (p < hi) {
        u64 g = q - p;
        if (gprev && g == gprev) {
            eq++;
            len++;
        } else {
            if (ours) finish_run(c, len, &st);
            ours = true;
            len = 1;
            st.off = idx;
            st.eqb = eq;
            st.p = p;
            st.d = g;
        }
        idx++;
        gprev = g;
        p = q;
        q = primesieve_next_prime(it);
    }
    c->nprimes = idx;
    c->neq = eq;
    if (ours) {                             /* follow the last run to its end */
        while (q - p == gprev) {
            len++;
            p = q;
            q = primesieve_next_prime(it);
        }
        finish_run(c, len, &st);
    }
    if (it->is_error) die("primesieve error in [%" PRIu64 ", %" PRIu64 ")", lo, hi);
}

/* ------------------------------------------------------------------ */
/* Totals and their printing                                           */
/* ------------------------------------------------------------------ */

typedef struct { bool set; u64 a, i, p, d; } term_t;

typedef struct {
    u64 pos;                    /* primes <= pos are included */
    u64 primes, eq;
    u64 hist[LMAX + 1];
    unsigned maxlen;
    term_t A[LMAX + 1];
} totals_t;

static void print_term(FILE *f, unsigned n, const term_t *t)
{
    char b1[32], b2[32], b3[32];
    fprintf(f, "a(%u) = %-14s  gap index i = %-14s  first prime %-14s  gap %" PRIu64 "\n",
            n, fmt_commas(t->a, b1), fmt_commas(t->i, b2), fmt_commas(t->p, b3), t->d);
}

static void print_ld(FILE *f, u64 pos, const u64 *ld, unsigned maxlen)
{
    for (unsigned L = 1; L <= LMAX && L <= maxlen; L++)
        for (unsigned d = 1; d < DMAX; d++) {
            u64 n = ld[(size_t)L * DMAX + d];
            if (n) fprintf(f, "ld %" PRIu64 " %u %u %" PRIu64 "\n", pos, L, d, n);
        }
}

static void print_totals(FILE *f, const totals_t *t, bool full, const u64 *ld)
{
    char b1[32], b2[32];
    u64 runs = t->primes - t->eq;
    fprintf(f, "primes <= %s (pi)         %s\n", fmt_commas(t->pos, b1), fmt_commas(t->primes, b2));
    fprintf(f, "equal adjacent gaps (E)   %s  (%.4f%% of the primes)\n", fmt_commas(t->eq, b1),
            t->primes ? 100.0 * (double)t->eq / (double)t->primes : 0.0);
    fprintf(f, "runs (= pi - E)           %s\n", fmt_commas(runs, b1));
    fprintf(f, "longest run               %u\n", t->maxlen);
    fprintf(f, "runs by length           ");
    for (unsigned L = 1; L <= LMAX && L <= t->maxlen; L++) {
        fprintf(f, " %u: %s", L, fmt_commas(t->hist[L], b1));
        if (L == LMAX && t->maxlen > LMAX) fprintf(f, " (>= %u)", LMAX);
    }
    fputc('\n', f);
    if (full) {
        unsigned missing = 0;
        for (unsigned n = 1; n <= LMAX; n++) {
            if (t->A[n].set) print_term(f, n, &t->A[n]);
            else if (!missing) missing = n;
        }
        if (missing)
            fprintf(f, "a(%u) > %s  (no run of length %u begins at a prime <= %s)\n",
                    missing, fmt_commas(runs, b1), missing, fmt_commas(t->pos, b2));
    }
    if (ld) print_ld(f, t->pos, ld, t->maxlen);
    fflush(f);
}

/* ------------------------------------------------------------------ */
/* Parallel scan                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    /* configuration */
    u64 start, end;             /* primes in [start, end] */
    u64 *bounds;                /* chunk i covers [bounds[i], bounds[i+1]) */
    u64 nchunks;
    u64 *marks;                 /* chunk boundaries at which cumulative totals are reported */
    int nmarks;
    int *seg;                   /* seg[i] = segment 0..nmarks of chunk i (for the ld sums) */
    bool want_ld, verbose, full, print_marks;
    int nthreads;

    /* shared state */
    atomic_uint_fast64_t next;
    atomic_int finished;
    pthread_mutex_t mu;
    pthread_cond_t cv;
    u64 cf;                     /* chunks folded so far, in order */
    chunk_t ring[RING];
    bool done[RING];
    u64 chunks_done;

    /* folded totals (chunks 0 .. cf-1) */
    u64 primes, eq;
    u64 hist[LMAX + 1];
    unsigned maxlen;
    term_t A[LMAX + 1];
    u64 *ld;                    /* (nmarks + 1) segments x NLD */
    u64 *ldtmp;                 /* NLD, scratch for cumulative sums */
    int next_mark;
    totals_t *marksnap;         /* totals at each mark */
    double t0, t1;
} scan_t;

typedef struct {
    int id;
    pthread_t th;
    scan_t *sc;
    primesieve_iterator it;
    chunk_t c;
    u64 *ld;
} worker_t;

static void snapshot(const scan_t *sc, u64 pos, totals_t *t)
{
    t->pos = pos;
    t->primes = sc->primes;
    t->eq = sc->eq;
    memcpy(t->hist, sc->hist, sizeof t->hist);
    t->maxlen = sc->maxlen;
    memcpy(t->A, sc->A, sizeof t->A);
}

/* cumulative (length, gap) counts of the segments 0..s */
static const u64 *ld_cumulative(scan_t *sc, int s)
{
    memset(sc->ldtmp, 0, NLD * sizeof(u64));
    for (int k = 0; k <= s; k++)
        for (size_t j = 0; j < NLD; j++) sc->ldtmp[j] += sc->ld[(size_t)k * NLD + j];
    return sc->ldtmp;
}

/* called with sc->mu held; chunk c is the next one in order */
static void fold(scan_t *sc, const chunk_t *c)
{
    for (unsigned n = 1; n <= LMAX; n++) {
        if (!c->first[n].set || sc->A[n].set) continue;
        u64 i = sc->primes + c->first[n].off + 1;      /* absolute index of the start prime */
        u64 E = sc->eq + c->first[n].eqb;              /* E(i-1) */
        sc->A[n] = (term_t){true, i - E, i, c->first[n].p, c->first[n].d};
        if (sc->full && sc->verbose) {
            pthread_mutex_lock(&out_mu);
            if (stderr_tty) fputs("\r\033[K", stderr);
            print_term(stdout, n, &sc->A[n]);
            fflush(stdout);
            pthread_mutex_unlock(&out_mu);
        }
    }
    sc->primes += c->nprimes;
    sc->eq += c->neq;
    for (unsigned L = 1; L <= LMAX; L++) sc->hist[L] += c->hist[L];
    if (c->maxlen > sc->maxlen) sc->maxlen = c->maxlen;

    if (sc->next_mark < sc->nmarks && c->hi == sc->marks[sc->next_mark]) {
        int m = sc->next_mark++;
        snapshot(sc, c->hi - 1, &sc->marksnap[m]);
        if (sc->print_marks) {
            char b1[32];
            pthread_mutex_lock(&out_mu);
            if (stderr_tty) fputs("\r\033[K", stderr);
            printf("--- cumulative totals for the primes <= %s ---\n", fmt_commas(c->hi - 1, b1));
            print_totals(stdout, &sc->marksnap[m], sc->full, sc->want_ld ? ld_cumulative(sc, m) : NULL);
            pthread_mutex_unlock(&out_mu);
        }
    }
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
        pthread_mutex_unlock(&sc->mu);
        if (g_stop) break;

        w->c.ld = sc->want_ld ? w->ld : NULL;
        process_chunk(&w->it, &w->c, sc->bounds[i], sc->bounds[i + 1]);

        pthread_mutex_lock(&sc->mu);
        if (sc->want_ld) {
            u64 *dst = sc->ld + (size_t)sc->seg[i] * NLD;
            for (size_t j = 0; j < NLD; j++) dst[j] += w->ld[j];
        }
        sc->ring[i % RING] = w->c;
        sc->ring[i % RING].ld = NULL;
        sc->done[i % RING] = true;
        while (sc->cf < sc->nchunks && sc->done[sc->cf % RING]) {
            fold(sc, &sc->ring[sc->cf % RING]);
            sc->done[sc->cf % RING] = false;
            sc->cf++;
        }
        sc->chunks_done++;
        pthread_cond_broadcast(&sc->cv);
        pthread_mutex_unlock(&sc->mu);
    }
    atomic_fetch_add(&sc->finished, 1);
    return NULL;
}

static void print_status(scan_t *sc, double t, bool final)
{
    pthread_mutex_lock(&sc->mu);
    u64 cf = sc->cf;
    u64 pos = cf >= sc->nchunks ? sc->end + 1 : sc->bounds[cf];
    u64 primes = sc->primes, runs = sc->primes - sc->eq;
    unsigned maxlen = sc->maxlen, found = 0;
    while (found < LMAX && sc->A[found + 1].set) found++;
    pthread_mutex_unlock(&sc->mu);

    double done = (double)(pos - sc->start), total = (double)(sc->end - sc->start) + 1.0;
    double el = t - sc->t0, rate = el > 0 ? done / el : 0, eta = rate > 0 ? (total - done) / rate : 0;
    char b1[32], b2[32], b3[32], b4[32], b5[32], b6[32];
    pthread_mutex_lock(&out_mu);
    if (stderr_tty) fputs("\r\033[K", stderr);
    fprintf(stderr, "[%s] %5.1f%% at %s  %s/s  ETA %s  primes %s  runs %s  longest %u",
            fmt_hms(el, b1), 100.0 * done / total, fmt_eng((double)pos, b2), fmt_eng(rate, b3),
            fmt_hms(eta, b4), fmt_eng((double)primes, b5), fmt_eng((double)runs, b6), maxlen);
    if (sc->full) fprintf(stderr, "  a(1..%u) found", found);
    if (!stderr_tty || final) fputc('\n', stderr);
    fflush(stderr);
    pthread_mutex_unlock(&out_mu);
}

static int cmp_u64(const void *a, const void *b)
{
    u64 x = *(const u64 *)a, y = *(const u64 *)b;
    return x < y ? -1 : x > y;
}

/*
 * Chunk boundaries: start, start + chunk, ...; every mark (a position M+1 with
 * start < M+1 <= end) becomes a boundary as well, so that the totals at a mark
 * are exactly the fold of the chunks before it.
 */
static void make_bounds(scan_t *sc, u64 chunk, const u64 *marks_in, int nmarks_in)
{
    u64 *marks = xcalloc((size_t)nmarks_in + 1, sizeof(u64));
    int nm = 0;
    for (int k = 0; k < nmarks_in; k++)
        if (marks_in[k] > sc->start && marks_in[k] <= sc->end) marks[nm++] = marks_in[k];
    qsort(marks, (size_t)nm, sizeof(u64), cmp_u64);
    int u = 0;
    for (int k = 0; k < nm; k++) if (!u || marks[k] != marks[u - 1]) marks[u++] = marks[k];
    nm = u;

    u64 range = sc->end - sc->start + 1;
    u64 cap = range / chunk + (u64)nm + 2;
    sc->bounds = xcalloc((size_t)cap + 1, sizeof(u64));
    sc->seg = xcalloc((size_t)cap + 1, sizeof(int));
    u64 pos = sc->start, n = 0;
    int mk = 0;
    for (;;) {
        sc->bounds[n] = pos;
        if (pos > sc->end) break;
        sc->seg[n] = mk;
        u64 nxt = (sc->end - pos >= chunk) ? pos + chunk : sc->end + 1;
        if (mk < nm && marks[mk] < nxt) nxt = marks[mk];
        if (mk < nm && marks[mk] == nxt) mk++;
        pos = nxt;
        n++;
    }
    sc->nchunks = n;
    sc->marks = marks;
    sc->nmarks = nm;
}

static scan_t *scan_new(u64 start, u64 end, int nthreads, u64 chunk, const u64 *marks, int nmarks,
                        bool want_ld, bool verbose, bool print_marks)
{
    if (end < start) die("END must be >= START");
    if (end >= primesieve_get_max_stop() - 1000000) die("END too large for primesieve");
    if (nthreads < 1) nthreads = 1;
    scan_t *sc = xcalloc(1, sizeof *sc);
    sc->start = start;
    sc->end = end;
    sc->nthreads = nthreads;
    sc->full = start <= 2;
    sc->want_ld = want_ld;
    sc->verbose = verbose;
    sc->print_marks = print_marks;
    if (!chunk) {
        u64 c = (end - start + 1) / ((u64)nthreads * 32);
        if (c < 1000000) c = 1000000;
        if (c > ((u64)1 << 30)) c = (u64)1 << 30;
        chunk = c;
    }
    make_bounds(sc, chunk, marks, nmarks);
    sc->marksnap = xcalloc((size_t)sc->nmarks + 1, sizeof(totals_t));
    if (want_ld) {
        sc->ld = xcalloc(((size_t)sc->nmarks + 1) * NLD, sizeof(u64));
        sc->ldtmp = xcalloc(NLD, sizeof(u64));
    }
    pthread_mutex_init(&sc->mu, NULL);
    pthread_cond_init(&sc->cv, NULL);
    return sc;
}

static void scan_free(scan_t *sc)
{
    free(sc->bounds);
    free(sc->seg);
    free(sc->marks);
    free(sc->marksnap);
    free(sc->ld);
    free(sc->ldtmp);
    pthread_mutex_destroy(&sc->mu);
    pthread_cond_destroy(&sc->cv);
    free(sc);
}

static void scan_run(scan_t *sc)
{
    atomic_store(&sc->next, 0);
    atomic_store(&sc->finished, 0);
    sc->t0 = now();

    worker_t *ws = xcalloc((size_t)sc->nthreads, sizeof(worker_t));
    for (int i = 0; i < sc->nthreads; i++) {
        ws[i].id = i;
        ws[i].sc = sc;
        primesieve_init(&ws[i].it);
        if (sc->want_ld) ws[i].ld = xcalloc(NLD, sizeof(u64));
        if (pthread_create(&ws[i].th, NULL, worker_main, &ws[i]) != 0) die("pthread_create failed");
    }

    double last_status = sc->t0;
    while (atomic_load(&sc->finished) < sc->nthreads) {
        usleep(100000);
        double t = now();
        if (g_stop) {
            pthread_mutex_lock(&sc->mu);
            pthread_cond_broadcast(&sc->cv);
            pthread_mutex_unlock(&sc->mu);
        }
        if (sc->verbose && ((stderr_tty && t - last_status >= 0.5) || (!stderr_tty && t - last_status >= 60))) {
            print_status(sc, t, false);
            last_status = t;
        }
    }
    for (int i = 0; i < sc->nthreads; i++) {
        pthread_join(ws[i].th, NULL);
        primesieve_free_iterator(&ws[i].it);
        free(ws[i].ld);
    }
    free(ws);
    sc->t1 = now();
    if (sc->verbose) print_status(sc, sc->t1, true);
}

/* position up to which the folded totals are complete */
static u64 scan_folded_pos(const scan_t *sc)
{
    return sc->cf >= sc->nchunks ? sc->end : sc->bounds[sc->cf] - 1;
}

/* ------------------------------------------------------------------ */
/* Reference implementation: one pass, no chunks                       */
/* ------------------------------------------------------------------ */

typedef struct {
    totals_t t;
    u64 runs;                   /* runs counted directly */
    u64 *ld;
    unsigned *runlen;           /* the first nkeep run lengths (A333254) */
    u64 nkeep, nrl;
} naive_t;

static void naive_close(naive_t *nv, unsigned len, u64 runidx, u64 i, u64 p, u64 d)
{
    unsigned L = len > LMAX ? LMAX : len;
    nv->t.hist[L]++;
    if (len > nv->t.maxlen) nv->t.maxlen = len;
    if (len <= LMAX && !nv->t.A[len].set) nv->t.A[len] = (term_t){true, runidx, i, p, d};
    if (nv->ld) nv->ld[(size_t)L * DMAX + (d < DMAX ? d : DMAX - 1)]++;
    if (runidx <= nv->nkeep) nv->runlen[nv->nrl++] = len;
}

static void naive_scan(naive_t *nv, u64 end, bool want_ld, u64 nkeep)
{
    memset(nv, 0, sizeof *nv);
    nv->t.pos = end;
    nv->nkeep = nkeep;
    if (nkeep) nv->runlen = xcalloc((size_t)nkeep, sizeof(unsigned));
    if (want_ld) nv->ld = xcalloc(NLD, sizeof(u64));

    primesieve_iterator it;
    primesieve_init(&it);
    primesieve_jump_to(&it, 0, end);
    u64 p = primesieve_next_prime(&it), q = primesieve_next_prime(&it);
    u64 k = 1, gprev = 0, runidx = 0, ri = 0, rp = 0, rd = 0;
    unsigned len = 0;
    while (p <= end) {
        u64 g = q - p;
        if (gprev && g == gprev) {
            nv->t.eq++;
            len++;
        } else {
            if (runidx) naive_close(nv, len, runidx, ri, rp, rd);
            runidx++;
            len = 1;
            ri = k;
            rp = p;
            rd = g;
        }
        nv->t.primes++;
        k++;
        gprev = g;
        p = q;
        q = primesieve_next_prime(&it);
    }
    while (runidx && q - p == gprev) {          /* finish the last run */
        len++;
        p = q;
        q = primesieve_next_prime(&it);
    }
    if (runidx) naive_close(nv, len, runidx, ri, rp, rd);
    nv->runs = runidx;
    if (it.is_error) die("primesieve error");
    primesieve_free_iterator(&it);
}

/* ------------------------------------------------------------------ */
/* Commands                                                            */
/* ------------------------------------------------------------------ */

static void usage(void)
{
    fputs(
"usage: a335406 scan [START] END [-t T] [-c CHUNK] [-m MARKS|pow10] [-H] [-q]\n"
"       a335406 runs N\n"
"       a335406 naive END [-H]\n"
"       a335406 selftest [-t T]\n"
"  -t T       worker threads (default: all cores)\n"
"  -c CHUNK   chunk size (default: range / (32 T), clamped to [1e6, 2^30])\n"
"  -m MARKS   comma-separated positions at which to print the cumulative totals,\n"
"             or 'pow10' for every power of ten in the range\n"
"  -H         also print the (run length, gap) table as lines 'ld POS LEN GAP COUNT'\n"
"  -q         no status line\n"
"Numbers: decimal, 2^k, 10^k, 1e12, X+Y, X-Y.\n", stderr);
    exit(2);
}

static int parse_marks(const char *s, u64 start, u64 end, u64 **out)
{
    u64 *m = xcalloc(64, sizeof(u64));
    int n = 0, cap = 64;
    if (strcmp(s, "pow10") == 0) {
        for (u64 v = 10; v > start && v <= end; v *= 10) {
            m[n++] = v + 1;                     /* boundary after the primes <= 10^k */
            if (v > UINT64_MAX / 10) break;
        }
        *out = m;
        return n;
    }
    char *copy = strdup(s), *save = NULL;
    for (char *tok = strtok_r(copy, ",", &save); tok; tok = strtok_r(NULL, ",", &save)) {
        if (n == cap) { cap *= 2; m = realloc(m, (size_t)cap * sizeof(u64)); if (!m) die("out of memory"); }
        u64 v = arg_u64(tok, "mark");
        if (v == UINT64_MAX) die("mark too large");
        m[n++] = v + 1;
    }
    free(copy);
    *out = m;
    return n;
}

static int cmd_scan(int argc, char **argv)
{
    int threads = 0;
    u64 chunk = 0, pos[2];
    int npos = 0;
    const char *marks_s = NULL;
    bool want_ld = false, quiet = false;
    for (int k = 0; k < argc; k++) {
        const char *a = argv[k];
        if (!strcmp(a, "-t") && k + 1 < argc) threads = arg_int(argv[++k], "thread count", 1, 4096);
        else if (!strcmp(a, "-c") && k + 1 < argc) chunk = arg_u64(argv[++k], "chunk size");
        else if (!strcmp(a, "-m") && k + 1 < argc) marks_s = argv[++k];
        else if (!strcmp(a, "-H")) want_ld = true;
        else if (!strcmp(a, "-q")) quiet = true;
        else if (a[0] == '-' && !isdigit((unsigned char)a[1])) die("unknown option '%s'", a);
        else if (npos < 2) pos[npos++] = arg_u64(a, "bound");
        else usage();
    }
    if (npos == 0) usage();
    if (chunk && chunk < 1000) die("chunk size too small");
    if (threads <= 0) {
        long n = sysconf(_SC_NPROCESSORS_ONLN);
        threads = n > 0 ? (int)n : 4;
    }
    u64 start = npos == 2 ? pos[0] : 0, end = pos[npos - 1];
    u64 *marks = NULL;
    int nmarks = marks_s ? parse_marks(marks_s, start, end, &marks) : 0;

    scan_t *sc = scan_new(start, end, threads, chunk, marks, nmarks, want_ld, !quiet, true);
    free(marks);
    char b1[32], b2[32], b3[32];
    printf("A335406 scan of the primes in [%s, %s]: %d threads, %" PRIu64 " chunks of %s\n",
           fmt_commas(start, b1), fmt_commas(end, b2), threads, sc->nchunks,
           fmt_commas(sc->nchunks ? sc->bounds[1] - sc->bounds[0] : 0, b3));
    if (!sc->full) puts("(START > 2: range statistics only, run indices are not absolute)");
    fflush(stdout);

    signal(SIGINT, on_sigint);
    scan_run(sc);

    totals_t t;
    u64 fpos = scan_folded_pos(sc);
    snapshot(sc, fpos, &t);
    double el = sc->t1 - sc->t0;
    if (g_stop) printf("\n*** interrupted: totals below cover the primes <= %s only ***\n", fmt_commas(fpos, b1));
    printf("--- totals for the primes in [%s, %s]  (%.1f s, %s numbers/s) ---\n", fmt_commas(start, b1),
           fmt_commas(fpos, b2), el, fmt_eng(el > 0 ? (double)(fpos - start + 1) / el : 0, b3));
    print_totals(stdout, &t, sc->full, want_ld && !g_stop ? ld_cumulative(sc, sc->nmarks) : NULL);
    scan_free(sc);
    return g_stop ? 1 : 0;
}

static int cmd_naive(int argc, char **argv)
{
    bool want_ld = false;
    u64 end = 0;
    int npos = 0;
    for (int k = 0; k < argc; k++) {
        if (!strcmp(argv[k], "-H")) want_ld = true;
        else if (npos == 0) { end = arg_u64(argv[k], "bound"); npos++; }
        else usage();
    }
    if (!npos) usage();
    naive_t nv;
    double t0 = now();
    naive_scan(&nv, end, want_ld, 0);
    char b1[32], b2[32];
    printf("A335406 reference pass over the primes <= %s  (%.1f s)\n", fmt_commas(end, b1), now() - t0);
    print_totals(stdout, &nv.t, true, nv.ld);
    printf("runs counted directly       %s\n", fmt_commas(nv.runs, b2));
    free(nv.ld);
    return 0;
}

static int cmd_runs(int argc, char **argv)
{
    if (argc != 1) usage();
    u64 n = arg_u64(argv[0], "count");
    if (n > 100000000) die("N too large");
    /* the k-th run begins at the k-th prime at the latest, and prime(k) <= 2 k ln k for k >= 6 */
    double bound = (double)n * (log((double)n + 6.0) + 1.0) * 2.0 + 1000.0;
    naive_t nv;
    naive_scan(&nv, (u64)bound, false, n);
    if (nv.nrl < n) die("internal: only %" PRIu64 " runs below the bound", nv.nrl);
    for (u64 k = 0; k < n; k++) printf("%s%u", k ? ", " : "", nv.runlen[k]);
    putchar('\n');
    free(nv.runlen);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Selftest                                                            */
/* ------------------------------------------------------------------ */

static int g_fail = 0;

static void check(bool ok, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    printf("%s ", ok ? "ok  " : "FAIL");
    vprintf(fmt, ap);
    putchar('\n');
    va_end(ap);
    if (!ok) g_fail++;
}

static bool totals_equal(const totals_t *a, const totals_t *b, bool with_terms)
{
    if (a->primes != b->primes || a->eq != b->eq || a->maxlen != b->maxlen) return false;
    for (unsigned L = 1; L <= LMAX; L++) if (a->hist[L] != b->hist[L]) return false;
    if (with_terms)
        for (unsigned n = 1; n <= LMAX; n++) {
            if (a->A[n].set != b->A[n].set) return false;
            if (a->A[n].set && (a->A[n].a != b->A[n].a || a->A[n].i != b->A[n].i ||
                                a->A[n].p != b->A[n].p || a->A[n].d != b->A[n].d)) return false;
        }
    return true;
}

static int cmd_selftest(int argc, char **argv)
{
    int threads = 0;
    for (int k = 0; k < argc; k++) {
        if (!strcmp(argv[k], "-t") && k + 1 < argc) threads = arg_int(argv[++k], "thread count", 1, 4096);
        else usage();
    }
    if (threads <= 0) {
        long n = sysconf(_SC_NPROCESSORS_ONLN);
        threads = n > 0 ? (int)n : 4;
    }
    const u64 N = 200000000, N1 = 100000000, SPLIT = 123456789;
    double t0 = now();

    /* 1. reference pass: the OEIS terms */
    naive_t nv, nv1;
    naive_scan(&nv, N, true, 60);
    for (unsigned n = 1; n <= 5; n++) {
        const term_t *t = &nv.t.A[n];
        check(t->set && t->a == KNOWN_A[n] && t->i == KNOWN_I[n] && t->p == KNOWN_P[n] && t->d == KNOWN_D[n],
              "a(%u) = %" PRIu64 " at gap index %" PRIu64 ", prime %" PRIu64 ", gap %" PRIu64
              " (OEIS: %" PRIu64 ", %" PRIu64 ", %" PRIu64 ", %" PRIu64 ")",
              n, t->a, t->i, t->p, t->d, KNOWN_A[n], KNOWN_I[n], KNOWN_P[n], KNOWN_D[n]);
    }
    check(nv.t.maxlen == 5 && !nv.t.A[6].set, "no run of length 6 below %" PRIu64 " (longest %u)", N, nv.t.maxlen);
    check(nv.runs == nv.t.primes - nv.t.eq, "runs counted directly (%" PRIu64 ") = pi - E (%" PRIu64 ")",
          nv.runs, nv.t.primes - nv.t.eq);
    u64 hsum = 0;
    for (unsigned L = 1; L <= LMAX; L++) hsum += nv.t.hist[L];
    check(hsum == nv.runs, "run-length histogram sums to the run count");
    check(nv.t.primes == 11078937, "pi(2e8) = %" PRIu64 " (expected 11078937)", nv.t.primes);

    /* 2. the first 60 terms of A333254 */
    char buf[512], *o = buf;
    for (u64 k = 0; k < 60 && k < nv.nrl; k++) o += sprintf(o, "%s%u", k ? ", " : "", nv.runlen[k]);
    check(nv.nrl >= 60 && !strcmp(buf, A333254_DATA), "first 60 run lengths match A333254");

    /* 3. chunked scans against the reference */
    struct { int t; u64 chunk; const char *what; } cfg[] = {
        {threads, 0, "default chunking"},
        {3, 12345, "3 threads, chunk 12345"},
        {threads, 1000003, "chunk 1000003"},
        {1, 77777777, "1 thread, chunk 77777777"},
    };
    for (size_t c = 0; c < sizeof cfg / sizeof cfg[0]; c++) {
        scan_t *sc = scan_new(0, N, cfg[c].t, cfg[c].chunk, NULL, 0, true, false, false);
        scan_run(sc);
        totals_t t;
        snapshot(sc, N, &t);
        bool ok = totals_equal(&t, &nv.t, true);
        bool ldok = memcmp(ld_cumulative(sc, 0), nv.ld, NLD * sizeof(u64)) == 0;
        check(ok && ldok, "scan to %" PRIu64 " with %s (%" PRIu64 " chunks) matches the reference%s", N,
              cfg[c].what, sc->nchunks, ldok ? "" : " (ld table differs)");
        scan_free(sc);
    }

    /* 4. marks: the cumulative totals at 1e8 equal a reference pass to 1e8 */
    naive_scan(&nv1, N1, true, 0);
    {
        u64 marks[1] = {N1 + 1};
        scan_t *sc = scan_new(0, N, threads, 3000000, marks, 1, true, false, false);
        scan_run(sc);
        bool ok = sc->nmarks == 1 && totals_equal(&sc->marksnap[0], &nv1.t, true) &&
                  memcmp(ld_cumulative(sc, 0), nv1.ld, NLD * sizeof(u64)) == 0;
        totals_t t;
        snapshot(sc, N, &t);
        check(ok, "cumulative totals at the mark 1e8 match a reference pass to 1e8");
        check(totals_equal(&t, &nv.t, true) && !memcmp(ld_cumulative(sc, 1), nv.ld, NLD * sizeof(u64)),
              "totals with a mark still match the reference at 2e8");
        scan_free(sc);
    }

    /* 5. split range: [0, SPLIT] + [SPLIT+1, N] */
    {
        scan_t *s1 = scan_new(0, SPLIT, threads, 5000000, NULL, 0, true, false, false);
        scan_t *s2 = scan_new(SPLIT + 1, N, threads, 5000000, NULL, 0, true, false, false);
        scan_run(s1);
        scan_run(s2);
        bool ok = s1->primes + s2->primes == nv.t.primes && s1->eq + s2->eq == nv.t.eq &&
                  (s1->maxlen > s2->maxlen ? s1->maxlen : s2->maxlen) == nv.t.maxlen;
        for (unsigned L = 1; L <= LMAX; L++) ok = ok && s1->hist[L] + s2->hist[L] == nv.t.hist[L];
        const u64 *l1 = ld_cumulative(s1, 0);
        for (size_t j = 0; j < NLD && ok; j++) ok = l1[j] + s2->ld[j] == nv.ld[j];
        check(ok, "split scans [0, %" PRIu64 "] + [%" PRIu64 ", %" PRIu64 "] add up to the reference",
              SPLIT, SPLIT + 1, N);
        scan_free(s1);
        scan_free(s2);
    }

    free(nv.ld);
    free(nv.runlen);
    free(nv1.ld);
    printf("%s (%.1f s)\n", g_fail ? "SELFTEST FAILED" : "selftest passed", now() - t0);
    return g_fail ? 1 : 0;
}

int main(int argc, char **argv)
{
    stderr_tty = isatty(2);
    if (argc < 2) usage();
    const char *cmd = argv[1];
    if (!strcmp(cmd, "scan")) return cmd_scan(argc - 2, argv + 2);
    if (!strcmp(cmd, "naive")) return cmd_naive(argc - 2, argv + 2);
    if (!strcmp(cmd, "runs")) return cmd_runs(argc - 2, argv + 2);
    if (!strcmp(cmd, "selftest")) return cmd_selftest(argc - 2, argv + 2);
    usage();
    return 2;
}
