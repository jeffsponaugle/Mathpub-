/*
 * a158939.c
 *
 * Compute and extend OEIS A158939:
 *
 *   "First primes followed by sequences of exactly n monotonic increasing
 *    prime gaps."
 *
 * Known terms (offset 0, n = 0..15):
 *   7, 3, 2, 17, 347, 2903, 15373, 128981, 1319407, 17797517, 94097537,
 *   6927837557, 48486712783, 968068681511, 1472840004017, 129001208165717
 *
 * Definition
 * ----------
 * For a prime p let p = q_0 < q_1 < q_2 < ... be the consecutive primes from p
 * on and g_i = q_i - q_(i-1) the gaps that follow p.  The run length of p is
 *
 *   L(p) = largest n such that g_1 < g_2 < ... < g_n,
 *
 * i.e. g_(n+1) <= g_n is the gap that ends the run (this is the PARI program
 * by Charles Greathouse in the entry).  a(n) = smallest prime with L(p) = n.
 * Examples: L(2) = 2 (gaps 1, 2, 2), L(3) = 1 (2, 2), L(17) = 3 (2, 4, 6, 2),
 * so a(2) = 2, a(1) = 3, a(3) = 17.  a(14) = 1472840004017 has the run
 * 2 4 6 8 10 12 14 28 30 38 48 64 66 74 ended by 34.
 *
 * No prime has L(p) = 0 under this definition; the entry's a(0) = 7 (gaps 4,
 * 2: the first prime whose next gap shrinks) is a convention of the original
 * author and is not computed here.  Related: A229832(n) is the prime after
 * a(n+1) (runs of n weak primes) and A133697(n) is the index of a(n+2) among
 * the primes.  a(n) is not monotone in principle, so the scan keeps the first
 * prime of every run length instead of stopping at a target.
 *
 * Method
 * ------
 * "First prime with ..." means every prime up to the answer has to be looked
 * at, so the cost is exactly one pass of a prime sieve over the range.
 * primesieve does that at ~3e9 numbers/s per core near 10^15 on an M1 Pro,
 * and the sieve is the whole cost: the bookkeeping per prime is a compare, an
 * increment and a store.
 *
 * The range is cut into chunks handed out through an atomic counter.  A
 * worker sieves its chunk [lo, hi) with a primesieve iterator, reads the
 * primes straight out of the iterator's buffer and keeps
 *
 *   run    = number of gaps in the strictly increasing run in progress,
 *   ring[] = the last 256 primes.
 *
 * When a gap g_(k+1) <= g_k arrives, the run g_(j+1) < ... < g_k of length
 * L = k - j has ended at q_k, and the primes q_(k-1), q_(k-2), ..., q_j have
 * run lengths 1, 2, ..., L.  Within a chunk the run lengths already seen are
 * always {1, ..., M} (a run of length L settles every n <= L), so a run needs
 * attention only when L > M: one compare per run.  A worker keeps going past
 * hi until the run in progress has ended, so every prime in [lo, hi) has its
 * run length settled by the worker that owns it; a prime's run length depends
 * only on the gaps after it, so nothing from before lo is needed.  Per-chunk
 * minima are merged into the global table under a mutex.  a(n) = p found in
 * chunk c is confirmed once every chunk before c is complete (the frontier);
 * the frontier is also what the checkpoint file records.
 *
 * Expected size of a(16)
 * ----------------------
 * Treating gaps as i.i.d. exponential, P(L(p) = n) = n/(n+1)!, so the first
 * run of length n appears near pi(x) = (n+1)!/n.  The known terms sit a
 * factor ~3 above that (ties between even gaps make strict increase rarer).
 * For n = 16 this puts a(16) near 2e15 (roughly 10% chance below 2.5e14, 50%
 * below 1.6e15, 90% below 5.5e15, 99% below 1.1e16); a(17) is expected near
 * 4e16.  Resta's scan to 1.3e14 (comment in A158940) found no run of 16.
 *
 * Usage
 * -----
 *   a158939 scan [START] END [-t T] [-c CHUNK] [-s KIB] [-n N] [-r NMIN] [-S FILE] [-i SECS] [-q]
 *       Exhaustive scan of the primes in [START, END) (START defaults to 0;
 *       END is rounded up to a chunk boundary).  Prints each new smallest
 *       prime with run length >= NMIN (default 14) as it is found, a
 *       "confirmed" line once everything below it has been examined, and a
 *       table of a(n) at the end.  -n N stops as soon as a(N) is confirmed.
 *       -S FILE keeps a checkpoint (every SECS seconds, default 60, and on
 *       Ctrl-C); rerunning the same command resumes from it.  -t T threads
 *       (default: all cores), -c CHUNK chunk size (default 10^10), -s KIB
 *       primesieve sieve size (default: library choice).
 *   a158939 verify P
 *       Print the gaps after the prime P and its run length L(P), computed
 *       with a Miller-Rabin next-prime search (independent of primesieve;
 *       also cross-checked with GMP when compiled in).
 *   a158939 bench [N] [-d SPAN] [-t T] [-c CHUNK] [-s KIB]
 *       Time the scan of [N, N+SPAN) (defaults 10^15 and 10^11) and project
 *       the running time of full scans.
 *   a158939 selftest [LIMIT] [-t T]
 *       Scan [0, LIMIT) (default 10^12) and compare with the known terms.
 *
 * Numbers may be written as decimal, 2^k, 10^k, 1e12, 1.29e14, or X+Y / X-Y
 * of those.
 *
 * Build:  cc -O2 -std=gnu11 -pthread -I/opt/homebrew/include a158939.c \
 *            -L/opt/homebrew/lib -lprimesieve -lm [-DHAVE_GMP -lgmp] -o a158939
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
#include <limits.h>
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
typedef uint16_t u16;
typedef unsigned __int128 u128;

#define NMAX      64        /* run lengths 1..NMAX-1 are tracked */
#define RINGP     256       /* recent primes kept by a worker (power of two, > NMAX+1) */
#define PMASK     (RINGP - 1)
#define RINGC     4096      /* chunk completion ring: bounds how far workers run ahead */
#define LOOKAHEAD 1000000   /* sieve hint past a chunk end, for the run in progress */
#define NKNOWN    15

static const u64 KNOWN[NKNOWN + 1] = {
    0, 3, 2, 17, 347, 2903, 15373, 128981, 1319407, 17797517, 94097537,
    6927837557ULL, 48486712783ULL, 968068681511ULL, 1472840004017ULL,
    129001208165717ULL
};

static bool stderr_tty;
static volatile sig_atomic_t g_stop = 0;

/* ------------------------------------------------------------------ */
/* Generic helpers                                                     */
/* ------------------------------------------------------------------ */

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void die(const char *fmt, ...)
{
    va_list ap;
    if (stderr_tty) fputs("\r\033[K", stderr);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

static const char *fmt_dur(double s, char *buf, size_t n)
{
    if (s < 0 || s != s || s > 1e9) { snprintf(buf, n, "?"); return buf; }
    long t = (long)s;
    if (t >= 86400) snprintf(buf, n, "%ldd%02ldh%02ldm", t / 86400, (t / 3600) % 24, (t / 60) % 60);
    else if (t >= 3600) snprintf(buf, n, "%ldh%02ldm%02lds", t / 3600, (t / 60) % 60, t % 60);
    else if (t >= 60) snprintf(buf, n, "%ldm%02lds", t / 60, t % 60);
    else snprintf(buf, n, "%.1fs", s);
    return buf;
}

static bool mul_ok(u64 a, u64 b, u64 *r)
{
    u128 x = (u128)a * b;
    if (x >> 64) return false;
    *r = (u64)x;
    return true;
}

static bool all_digits(const char *s, size_t n)
{
    if (n == 0) return false;
    for (size_t i = 0; i < n; i++)
        if (!isdigit((unsigned char)s[i])) return false;
    return true;
}

/* decimal, B^E, MeE (M may have a decimal point) */
static bool parse_term(const char *s, size_t n, u64 *out)
{
    char buf[64], *c;
    if (n == 0 || n >= sizeof buf) return false;
    memcpy(buf, s, n);
    buf[n] = 0;
    if ((c = strchr(buf, '^'))) {
        *c = 0;
        if (!all_digits(buf, strlen(buf)) || !all_digits(c + 1, strlen(c + 1))) return false;
        u64 b = strtoull(buf, 0, 10), e = strtoull(c + 1, 0, 10), r = 1;
        while (e--) if (!mul_ok(r, b, &r)) return false;
        *out = r;
        return true;
    }
    if ((c = strpbrk(buf, "eE"))) {
        *c = 0;
        if (!all_digits(c + 1, strlen(c + 1))) return false;
        long e = strtol(c + 1, 0, 10);
        char mant[64];
        int md = 0, frac = -1;
        for (const char *q = buf; *q; q++) {
            if (*q == '.') { if (frac >= 0) return false; frac = 0; continue; }
            if (!isdigit((unsigned char)*q)) return false;
            mant[md++] = *q;
            if (frac >= 0) frac++;
        }
        if (md == 0 || md > 19) return false;
        mant[md] = 0;
        if (frac > 0) e -= frac;
        if (e < 0) return false;
        u64 m = strtoull(mant, 0, 10);
        while (e-- > 0) if (!mul_ok(m, 10, &m)) return false;
        *out = m;
        return true;
    }
    if (!all_digits(buf, n) || n > 20) return false;
    errno = 0;
    u64 v = strtoull(buf, 0, 10);
    if (errno) return false;
    *out = v;
    return true;
}

static u64 parse_num(const char *s)
{
    size_t n = strlen(s);
    for (size_t i = n; i-- > 1; ) {           /* split at the last top-level + or - */
        if (s[i] != '+' && s[i] != '-') continue;
        u64 a, b;
        if (!parse_term(s, i, &a) || !parse_term(s + i + 1, n - i - 1, &b)) die("bad number '%s'", s);
        if (s[i] == '+') {
            if (a + b < a) die("overflow in '%s'", s);
            return a + b;
        }
        if (b > a) die("negative result in '%s'", s);
        return a - b;
    }
    u64 v;
    if (!parse_term(s, n, &v)) die("bad number '%s'", s);
    return v;
}

static int parse_int(const char *s, int lo, int hi, const char *what)
{
    u64 v = parse_num(s);
    if (v < (u64)lo || v > (u64)hi) die("%s must be between %d and %d", what, lo, hi);
    return (int)v;
}

/* ------------------------------------------------------------------ */
/* Independent primality: deterministic Miller-Rabin below 2^64         */
/* ------------------------------------------------------------------ */

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)(((u128)a * b) % m); }

static u64 powmod(u64 a, u64 e, u64 m)
{
    u64 r = 1;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

static bool is_prime64(u64 n)
{
    static const u64 small[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
    static const u64 bases[] = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };
    if (n < 2) return false;
    for (size_t i = 0; i < sizeof small / sizeof *small; i++)
        if (n % small[i] == 0) return n == small[i];
    u64 d = n - 1;
    int r = 0;
    while (!(d & 1)) { d >>= 1; r++; }
    for (size_t i = 0; i < sizeof bases / sizeof *bases; i++) {
        u64 a = bases[i] % n;
        if (a == 0) continue;
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int j = 1; j < r && composite; j++) {
            x = mulmod(x, x, n);
            if (x == n - 1) composite = false;
        }
        if (composite) return false;
    }
    return true;
}

static u64 next_prime64(u64 n)
{
    if (n < 2) return 2;
    u64 m = n + 1;
    if (!(m & 1)) m++;
    while (!is_prime64(m)) m += 2;
    return m;
}

/* Run length of the prime p by next-prime search; fills gaps[0..L] (L+1 values,
 * the last one being the gap that ends the run).  Returns L, or -1 if the run
 * is longer than maxgaps - 1. */
static int run_length64(u64 p, u16 *gaps, int maxgaps)
{
    u64 prev = p, prevgap = 0;
    int L = 0;
    for (;;) {
        u64 q = next_prime64(prev), g = q - prev;
        if (L >= maxgaps) return -1;
        gaps[L] = (u16)g;
        if (g <= prevgap) return L;
        L++;
        prevgap = g;
        prev = q;
    }
}

/* ------------------------------------------------------------------ */
/* Scan engine                                                         */
/* ------------------------------------------------------------------ */

typedef struct {
    u64  p;                 /* smallest prime seen with this run length, 0 = none */
    u64  chunk;             /* chunk it was found in */
    u64  idx;               /* number of primes in [chunk start, p] */
    u64  pi;                /* number of primes in [START, p], known once the frontier passes chunk */
    u16  gaps[NMAX + 1];    /* the n gaps of the run, then the gap that ends it */
    bool announced;         /* "confirmed" line printed */
} cand_t;

/* A133697(n) = pi(A158939(n+2)) for n = 0..13 */
static const u64 A133697[14] = {
    1, 7, 69, 420, 1796, 12073, 101397, 1139211, 5440508, 320620306ULL,
    2058187481ULL, 36451609409ULL, 54594153615ULL, 4100904808215ULL
};

static struct {
    u64 start, end, chunk, nchunks;
    _Atomic u64 next;               /* next chunk index to hand out */
    _Atomic bool overflow;          /* a run longer than NMAX-1 gaps was seen */
    /* protected by g_mu: */
    u64 frontier;                   /* chunks [0, frontier) are complete */
    struct {                        /* completed chunks waiting for the frontier, by chunk % RINGC */
        bool done;
        u64 primes;
        u64 hist[NMAX];
    } slot[RINGC];
    u64 primes;                     /* primes in [start, frontier position) */
    u64 hist[NMAX];                 /* hist[n] = primes with run length n in [start, frontier position) */
    cand_t best[NMAX];
} S;
static pthread_mutex_t g_mu = PTHREAD_MUTEX_INITIALIZER;
static int g_report = 14;           /* print new minima for run lengths >= this */
static int g_sieve_kib = 0;         /* primesieve sieve size in KiB, 0 = library default */

typedef struct {
    int id;
    primesieve_iterator it;
} worker_t;

static void print_cand(int n, const cand_t *c, const char *tag)
{
    if (stderr_tty) fputs("\r\033[K", stderr);
    printf("%s a(%d) = %" PRIu64 "  gaps:", tag, n, c->p);
    for (int i = 0; i < n; i++) printf(" %u", c->gaps[i]);
    printf(" | %u  (chunk %" PRIu64, c->gaps[n], c->chunk);
    if (c->pi) printf(", pi(p) = %" PRIu64, c->pi);
    printf(")\n");
    fflush(stdout);
}

static void scan_chunk(worker_t *w, u64 idx)
{
    const u64 lo = S.start + idx * S.chunk;
    u64 hi = lo + S.chunk;
    if (hi > S.end || hi < lo) hi = S.end;

    u64 first[NMAX], first_idx[NMAX];
    u16 gaps[NMAX][NMAX + 1];
    u64 ring[RINGP];
    u64 runs[NMAX];                 /* runs[L] = run ends of length L with all L primes < hi */
    u64 extra[NMAX];                /* run lengths from the final run end, per prime < hi */
    unsigned ri = 0;
    int M = 0, run = 0;
    u64 prev, prevgap = 0, count = 0;

    memset(first, 0, sizeof first);
    memset(runs, 0, sizeof runs);
    memset(extra, 0, sizeof extra);
    primesieve_jump_to(&w->it, lo, hi + LOOKAHEAD);
    /* The primes are read straight out of the iterator's buffer, which
     * primesieve_generate_next_primes() fills with primes[0..size).  (Mixing in
     * primesieve_next_prime() for the first prime would discard the rest of
     * that first buffer.) */
    primesieve_generate_next_primes(&w->it);
    if (w->it.is_error || w->it.size == 0) die("primesieve error in chunk %" PRIu64, idx);
    const u64 *P = w->it.primes;
    size_t np = w->it.size, i = 1;
    prev = P[0];
    if (prev >= hi) goto finish;                    /* no prime in this chunk */
    ring[ri++ & PMASK] = prev;
    count = 1;

    for (;;) {
        for (; i < np; i++) {
            const u64 p = P[i];
            const u64 g = p - prev;
            if (g > prevgap) {
                run++;
            } else {
                /* The run of `run` gaps ended at prev = q_k (= ring[ri-1]);
                 * prime ring[ri-1-n] has run length n for n = 1..run. */
                if (run >= NMAX) { atomic_store(&S.overflow, true); run = NMAX - 1; }
                if (prev < hi) {
                    runs[run]++;
                } else {                            /* final run end: some primes belong to the next chunk */
                    for (int n = 1; n <= run; n++)
                        if (ring[(ri - 1 - n) & PMASK] < hi) extra[n]++;
                }
                if (run > M) {
                    for (int n = M + 1; n <= run; n++) {
                        const u64 q = ring[(ri - 1 - n) & PMASK];
                        if (q >= hi) continue;          /* belongs to the next chunk */
                        first[n] = q;
                        first_idx[n] = count - (u64)n;  /* count = primes in [lo, prev] */
                        for (int m = 0; m < n; m++)
                            gaps[n][m] = (u16)(ring[(ri - n + m) & PMASK] - ring[(ri - 1 - n + m) & PMASK]);
                        gaps[n][n] = (u16)g;
                    }
                    M = run;
                }
                run = 1;
                if (prev >= hi) goto finish;        /* every prime < hi is settled */
            }
            prevgap = g;
            prev = p;
            ring[ri++ & PMASK] = p;
            count++;
        }
        primesieve_generate_next_primes(&w->it);
        if (w->it.is_error || w->it.size == 0) die("primesieve error in chunk %" PRIu64, idx);
        P = w->it.primes;
        np = w->it.size;
        i = 0;
    }

finish:
    for (unsigned j = 1; j <= ri && j < RINGP; j++) {   /* don't count primes >= hi */
        if (ring[(ri - j) & PMASK] >= hi) count--;
        else break;
    }

    /* primes with run length exactly n: one per run end of length >= n */
    u64 hist[NMAX];
    {
        u64 acc = 0;
        for (int n = NMAX - 1; n >= 1; n--) {
            acc += runs[n];
            hist[n] = acc + extra[n];
        }
        hist[0] = 0;
    }

    pthread_mutex_lock(&g_mu);
    for (int n = 1; n <= M && n < NMAX; n++) {
        if (!first[n]) continue;
        cand_t *c = &S.best[n];
        if (c->p == 0 || first[n] < c->p) {
            c->p = first[n];
            c->chunk = idx;
            c->idx = first_idx[n];
            c->pi = 0;
            c->announced = false;
            memcpy(c->gaps, gaps[n], (size_t)(n + 1) * sizeof(u16));
            if (n >= g_report) print_cand(n, c, "found    ");
        }
    }
    {
        unsigned s = (unsigned)(idx % RINGC);
        S.slot[s].done = true;
        S.slot[s].primes = count;
        memcpy(S.slot[s].hist, hist, sizeof hist);
    }
    while (S.slot[S.frontier % RINGC].done) {       /* fold complete chunks in frontier order */
        unsigned s = (unsigned)(S.frontier % RINGC);
        S.slot[s].done = false;
        for (int n = 1; n < NMAX; n++)                  /* S.primes = primes in [start, chunk start) */
            if (S.best[n].p && S.best[n].chunk == S.frontier) S.best[n].pi = S.primes + S.best[n].idx;
        S.primes += S.slot[s].primes;
        for (int n = 1; n < NMAX; n++) S.hist[n] += S.slot[s].hist[n];
        S.frontier++;
    }
    pthread_mutex_unlock(&g_mu);
}

static void *worker(void *arg)
{
    worker_t *w = arg;
    if (g_sieve_kib) primesieve_set_sieve_size(g_sieve_kib);
    primesieve_init(&w->it);
    while (!g_stop) {
        u64 idx = atomic_fetch_add(&S.next, 1);
        if (idx >= S.nchunks) break;
        for (;;) {                                   /* stay within RINGC of the frontier */
            pthread_mutex_lock(&g_mu);
            u64 f = S.frontier;
            pthread_mutex_unlock(&g_mu);
            if (idx < f + RINGC - 1 || g_stop) break;
            usleep(2000);
        }
        if (g_stop) break;
        scan_chunk(w, idx);
    }
    primesieve_free_iterator(&w->it);
    return NULL;
}

/* ---- checkpoint file ---- */

static void save_state(const char *path)
{
    char tmp[PATH_MAX];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { fprintf(stderr, "\nwarning: cannot write %s: %s\n", tmp, strerror(errno)); return; }
    pthread_mutex_lock(&g_mu);
    fprintf(f, "A158939 checkpoint 2\nstart %" PRIu64 "\nchunk %" PRIu64 "\nend %" PRIu64
               "\nfrontier %" PRIu64 "\nprimes %" PRIu64 "\n",
            S.start, S.chunk, S.end, S.frontier, S.primes);
    for (int n = 1; n < NMAX; n++)
        if (S.hist[n]) fprintf(f, "hist %d %" PRIu64 "\n", n, S.hist[n]);
    for (int n = 1; n < NMAX; n++) {
        if (!S.best[n].p) continue;
        fprintf(f, "cand %d %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64, n, S.best[n].p, S.best[n].chunk,
                S.best[n].idx, S.best[n].pi);
        for (int i = 0; i <= n; i++) fprintf(f, " %u", S.best[n].gaps[i]);
        fputc('\n', f);
    }
    pthread_mutex_unlock(&g_mu);
    fputs("end\n", f);
    if (fclose(f) != 0 || rename(tmp, path) != 0)
        fprintf(stderr, "\nwarning: cannot update %s: %s\n", path, strerror(errno));
}

/* Returns false if there is no checkpoint; dies on a mismatch. */
static bool load_state(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        if (errno == ENOENT) return false;
        die("cannot read %s: %s", path, strerror(errno));
    }
    char line[2048];
    if (!fgets(line, sizeof line, f) || strncmp(line, "A158939 checkpoint 2", 20) != 0)
        die("%s is not an A158939 checkpoint", path);
    u64 start = 0, chunk = 0, end = 0, frontier = 0, primes = 0;
    bool complete = false;
    while (fgets(line, sizeof line, f)) {
        if (sscanf(line, "start %" SCNu64, &start) == 1) continue;
        if (sscanf(line, "chunk %" SCNu64, &chunk) == 1) continue;
        if (sscanf(line, "end %" SCNu64, &end) == 1) continue;
        if (sscanf(line, "frontier %" SCNu64, &frontier) == 1) continue;
        if (sscanf(line, "primes %" SCNu64, &primes) == 1) continue;
        if (strncmp(line, "hist ", 5) == 0) {
            int n = 0;
            u64 c = 0;
            if (sscanf(line, "hist %d %" SCNu64, &n, &c) != 2 || n < 1 || n >= NMAX) die("%s: bad hist line", path);
            S.hist[n] = c;
            continue;
        }
        if (strncmp(line, "cand ", 5) == 0) {
            char *s = line + 5, *e;
            long n = strtol(s, &e, 10);
            if (n < 1 || n >= NMAX) die("%s: bad candidate line", path);
            cand_t *c = &S.best[n];
            c->p = strtoull(e, &e, 10);
            c->chunk = strtoull(e, &e, 10);
            c->idx = strtoull(e, &e, 10);
            c->pi = strtoull(e, &e, 10);
            for (int i = 0; i <= n; i++) c->gaps[i] = (u16)strtoul(e, &e, 10);
            c->announced = false;
            continue;
        }
        if (strncmp(line, "end", 3) == 0) { complete = true; break; }
    }
    fclose(f);
    if (!complete) die("%s is truncated", path);
    if (start != S.start || chunk != S.chunk)
        die("%s was written for start=%" PRIu64 " chunk=%" PRIu64 ", not start=%" PRIu64 " chunk=%" PRIu64,
            path, start, chunk, S.start, S.chunk);
    if (frontier > S.nchunks) frontier = S.nchunks;
    S.frontier = frontier;
    S.primes = primes;
    return true;
}

static void on_sigint(int sig)
{
    (void)sig;
    g_stop = 1;
}

/* Announce confirmed terms; returns the largest confirmed n (0 if none). */
static int announce_locked(void)
{
    int maxc = 0;
    for (int n = 1; n < NMAX; n++) {
        cand_t *c = &S.best[n];
        if (!c->p || c->chunk >= S.frontier) continue;
        if (n > maxc) maxc = n;
        if (!c->announced) {
            c->announced = true;
            if (n >= g_report) print_cand(n, c, "CONFIRMED");
        }
    }
    return maxc;
}

typedef struct {
    double seconds;
    u64 primes;
    u64 numbers;
} scan_stats_t;

/* Scan [start, end) with the given number of threads.  stop_n > 0 ends the
 * scan as soon as a(stop_n) is confirmed. */
static scan_stats_t run_scan(u64 start, u64 end, u64 chunk, int threads, const char *state,
                             int interval, int stop_n, bool quiet)
{
    if (end <= start) die("END must be larger than START");
    if (chunk == 0) die("CHUNK must be positive");
    S.start = start;
    S.chunk = chunk;
    S.nchunks = (end - start + chunk - 1) / chunk;
    if (S.nchunks == 0 || (S.nchunks - 1) > (UINT64_MAX - start) / chunk) die("range too large");
    S.end = start + S.nchunks * chunk;
    if (S.end < start) die("range too large");
    if (S.end > primesieve_get_max_stop() - 2 * LOOKAHEAD) die("END beyond primesieve's limit");
    if (S.end != end)
        fprintf(stderr, "note: END rounded up to the chunk boundary %" PRIu64 "\n", S.end);
    S.frontier = 0;
    S.primes = 0;
    memset(S.slot, 0, sizeof S.slot);
    memset(S.hist, 0, sizeof S.hist);
    memset(S.best, 0, sizeof S.best);
    atomic_store(&S.overflow, false);

    if (state && load_state(state))
        fprintf(stderr, "resuming from %s: %" PRIu64 " of %" PRIu64 " chunks done (position %" PRIu64
                        ", %" PRIu64 " primes examined)\n",
                state, S.frontier, S.nchunks, S.start + S.frontier * S.chunk, S.primes);
    const u64 frontier0 = S.frontier;
    atomic_store(&S.next, S.frontier);

    if (!quiet)
        fprintf(stderr, "scanning [%" PRIu64 ", %" PRIu64 ") in %" PRIu64 " chunks of %" PRIu64
                        " with %d threads%s\n", S.start, S.end, S.nchunks, S.chunk, threads,
                stop_n ? " (stopping when a(N) is confirmed)" : "");

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_sigint;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    worker_t *ws = calloc((size_t)threads, sizeof *ws);
    pthread_t *th = calloc((size_t)threads, sizeof *th);
    if (!ws || !th) die("out of memory");
    for (int i = 0; i < threads; i++) {
        ws[i].id = i;
        if (pthread_create(&th[i], NULL, worker, &ws[i]) != 0) die("pthread_create failed");
    }

    /* rate window: one sample per second, ~2 minutes deep */
    enum { NS = 128 };
    double ts[NS];
    u64 ps[NS];
    int ns = 0;
    const double t0 = now();
    const u64 pos0 = S.start + frontier0 * S.chunk;
    double last_status = 0, last_ckpt = t0, last_sample = 0;
    ts[0] = t0; ps[0] = pos0; ns = 1;
    char b1[32], b2[32];

    for (;;) {
        usleep(200000);
        double t = now();
        pthread_mutex_lock(&g_mu);
        u64 f = S.frontier, primes = S.primes;
        int maxc = announce_locked();
        if (stop_n > 0 && S.best[stop_n].p && S.best[stop_n].chunk < f) g_stop = 1;
        pthread_mutex_unlock(&g_mu);
        u64 pos = S.start + f * S.chunk;
        if (f >= S.nchunks || g_stop) break;

        if (t - last_sample >= 1.0) {
            last_sample = t;
            if (ns < NS) { ts[ns] = t; ps[ns] = pos; ns++; }
            else { memmove(ts, ts + 1, (NS - 1) * sizeof *ts); memmove(ps, ps + 1, (NS - 1) * sizeof *ps); ts[NS - 1] = t; ps[NS - 1] = pos; }
        }
        double rate = (t - ts[0] > 0.5) ? (double)(pos - ps[0]) / (t - ts[0]) : 0;
        double eta = rate > 0 ? (double)(S.end - pos) / rate : -1;
        bool show = stderr_tty ? (t - last_status >= 1.0) : (t - last_status >= interval);
        if (show && !quiet) {
            last_status = t;
            fprintf(stderr, "%s%s  pos %.6g (%.2f%%)  %.3g/s  %.3g primes  ETA %s  confirmed to n=%d",
                    stderr_tty ? "\r\033[K" : "progress: ",
                    fmt_dur(t - t0, b1, sizeof b1), (double)pos,
                    100.0 * (double)(pos - S.start) / (double)(S.end - S.start),
                    rate, (double)primes, fmt_dur(eta, b2, sizeof b2), maxc);
            if (stop_n > 0) {
                if (S.best[stop_n].p) fprintf(stderr, "  a(%d) candidate %" PRIu64, stop_n, S.best[stop_n].p);
                else fprintf(stderr, "  a(%d) > %.4g", stop_n, (double)pos);
            }
            if (!stderr_tty) fputc('\n', stderr);
            fflush(stderr);
        }
        if (state && t - last_ckpt >= interval) {
            last_ckpt = t;
            save_state(state);
        }
    }
    for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
    double t1 = now();
    if (stderr_tty && !quiet) fputs("\r\033[K", stderr);

    pthread_mutex_lock(&g_mu);
    announce_locked();
    pthread_mutex_unlock(&g_mu);
    if (state) save_state(state);
    if (atomic_load(&S.overflow))
        fprintf(stderr, "warning: a run longer than %d gaps occurred; a(n) for n >= %d is not tracked\n", NMAX - 1, NMAX);
    if (g_stop && S.frontier < S.nchunks && !quiet)
        fprintf(stderr, "stopped at position %" PRIu64 " (%" PRIu64 " of %" PRIu64 " chunks complete)%s\n",
                S.start + S.frontier * S.chunk, S.frontier, S.nchunks,
                state ? ", checkpoint saved" : "");

    scan_stats_t st;
    st.seconds = t1 - t0;
    st.primes = S.primes;
    st.numbers = (S.frontier - frontier0) * S.chunk;
    free(ws);
    free(th);
    return st;
}

/* pi(10^k) for k = 6..16, for checking the prime count of a scan [0, 10^k) */
static const u64 PI10[17] = {
    0, 0, 0, 0, 0, 0, 78498ULL, 664579ULL, 5761455ULL, 50847534ULL, 455052511ULL,
    4118054813ULL, 37607912018ULL, 346065536839ULL, 3204941750802ULL,
    29844570422669ULL, 279238341033925ULL
};

static void print_hist(void)
{
    const u64 pos = S.start + S.frontier * S.chunk;
    printf("\n%" PRIu64 " primes in [%" PRIu64 ", %" PRIu64 ")", S.primes, S.start, pos);
    if (S.start == 0)
        for (int k = 6; k <= 16; k++)
            if (pos == (u64)pow(10, k)) printf(S.primes == PI10[k] ? "  = pi(10^%d)" : "  BUT pi(10^%d) = %" PRIu64, k, PI10[k]);
    printf("\nrun-length distribution (model: gaps i.i.d. exponential, P(L = n) = n/(n+1)!):\n");
    printf("  n  primes with L(p) = n     model        ratio\n");
    double fact = 1;                                 /* becomes (n+1)! */
    for (int n = 1; n < NMAX; n++) {
        fact *= (n + 1);
        double model = (double)S.primes * n / fact;
        if (S.hist[n] == 0 && model < 0.01) break;
        printf(" %2d  %-24" PRIu64 " %-12.4g %.3f\n", n, S.hist[n], model, model > 0 ? (double)S.hist[n] / model : 0);
    }
    fflush(stdout);
}

static void print_table(void)
{
    const u64 pos = S.start + S.frontier * S.chunk;
    int maxn = 0;
    for (int n = 1; n < NMAX; n++) if (S.best[n].p) maxn = n;
    printf("\n %2s  %-24s %-21s %s\n", "n", S.start ? "a(n) (>= START)" : "a(n)",
           S.start ? "pi(a(n)) (from START)" : "pi(a(n))", "gaps | gap that ends the run");
    for (int n = 1; n <= maxn; n++) {
        const cand_t *c = &S.best[n];
        if (!c->p) { printf(" %2d  (none found)\n", n); continue; }
        bool confirmed = c->chunk < S.frontier;
        printf(" %2d  %-24" PRIu64, n, c->p);
        if (c->pi) printf(" %-21" PRIu64, c->pi); else printf(" %-21s", "?");
        for (int i = 0; i < n; i++) printf(" %u", c->gaps[i]);
        printf(" | %u", c->gaps[n]);
        if (S.start == 0 && n <= NKNOWN)
            printf(c->p == KNOWN[n] ? "  = OEIS" : "  DIFFERS from OEIS %" PRIu64, KNOWN[n]);
        else if (S.start == 0)
            printf("  NEW");
        if (S.start == 0 && c->pi && n >= 2 && n - 2 < 14)
            printf(c->pi == A133697[n - 2] ? ", pi = A133697(%d)" : ", pi DIFFERS from A133697(%d) = %" PRIu64,
                   n - 2, A133697[n - 2]);
        else if (S.start == 0 && c->pi && n >= 2)
            printf(", pi = new A133697(%d)", n - 2);
        if (!confirmed) printf("  [unconfirmed: earlier chunks incomplete]");
        putchar('\n');
    }
    if (S.frontier > 0)
        printf(" %2d  > %" PRIu64 "  (no%s prime below this position has run length %d)\n",
               maxn + 1, pos, S.start ? " such" : "", maxn + 1);
    print_hist();
}

/* Brute-force reference: first prime of every run length below limit, from
 * the plain definition applied to an array of primes. */
static void reference_first(u64 limit, u64 *first, u64 *hist)
{
    size_t n;
    u64 *pr = primesieve_generate_primes(0, limit + 1000000, &n, UINT64_PRIMES);
    if (!pr) die("primesieve_generate_primes failed");
    memset(first, 0, NMAX * sizeof *first);
    memset(hist, 0, NMAX * sizeof *hist);
    for (size_t i = 0; i < n && pr[i] < limit; i++) {
        int L = 0;
        u64 prevgap = 0;
        size_t j;
        for (j = i + 1; j < n; j++) {
            u64 g = pr[j] - pr[j - 1];
            if (g <= prevgap) break;
            L++;
            prevgap = g;
        }
        if (j == n) die("reference_first: run at %" PRIu64 " not settled", pr[i]);
        if (L < NMAX && first[L] == 0) first[L] = pr[i];
        if (L < NMAX) hist[L]++;
    }
    primesieve_free(pr);
}

/* ------------------------------------------------------------------ */
/* Commands                                                            */
/* ------------------------------------------------------------------ */

static int default_threads(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n < 1 ? 1 : (n > 1024 ? 1024 : (int)n);
}

static void usage(void)
{
    fprintf(stderr,
        "usage: a158939 scan [START] END [-t T] [-c CHUNK] [-s KIB] [-n N] [-r NMIN] [-S FILE] [-i SECS] [-q]\n"
        "       a158939 verify P\n"
        "       a158939 bench [N] [-d SPAN] [-t T] [-c CHUNK] [-s KIB]\n"
        "       a158939 selftest [LIMIT] [-t T]\n");
    exit(2);
}

static int cmd_scan(int argc, char **argv)
{
    u64 pos[2];
    int npos = 0, threads = default_threads(), interval = 60, stop_n = 0;
    u64 chunk = 10000000000ULL;
    const char *state = NULL;
    bool quiet = false;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) threads = parse_int(argv[++i], 1, 1024, "threads");
        else if (!strcmp(argv[i], "-c") && i + 1 < argc) chunk = parse_num(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i + 1 < argc) g_sieve_kib = parse_int(argv[++i], 16, 65536, "sieve KiB");
        else if (!strcmp(argv[i], "-n") && i + 1 < argc) stop_n = parse_int(argv[++i], 1, NMAX - 1, "N");
        else if (!strcmp(argv[i], "-r") && i + 1 < argc) g_report = parse_int(argv[++i], 1, 1000, "NMIN");
        else if (!strcmp(argv[i], "-S") && i + 1 < argc) state = argv[++i];
        else if (!strcmp(argv[i], "-i") && i + 1 < argc) interval = parse_int(argv[++i], 1, 86400, "SECS");
        else if (!strcmp(argv[i], "-q")) quiet = true;
        else if (argv[i][0] == '-' && !isdigit((unsigned char)argv[i][1])) usage();
        else if (npos < 2) pos[npos++] = parse_num(argv[i]);
        else usage();
    }
    if (npos == 0) usage();
    u64 start = npos == 2 ? pos[0] : 0, end = pos[npos - 1];
    if (chunk < 1000000) die("CHUNK must be at least 10^6");
    scan_stats_t st = run_scan(start, end, chunk, threads, state, interval, stop_n, quiet);
    print_table();
    char b[32];
    fprintf(stderr, "%" PRIu64 " numbers, %" PRIu64 " primes in %s (%.3g numbers/s)\n",
            st.numbers, st.primes, fmt_dur(st.seconds, b, sizeof b),
            st.seconds > 0 ? (double)st.numbers / st.seconds : 0);
    return 0;
}

static int cmd_verify(int argc, char **argv)
{
    if (argc < 1) usage();
    u64 p = parse_num(argv[0]);
    if (!is_prime64(p)) die("%" PRIu64 " is not prime", p);
    u16 gaps[NMAX + 1];
    int L = run_length64(p, gaps, NMAX + 1);
    if (L < 0) die("run longer than %d gaps", NMAX);
    printf("%" PRIu64 ": gaps", p);
    for (int i = 0; i < L; i++) printf(" %u", gaps[i]);
    printf(" | %u  ->  L(p) = %d", gaps[L], L);
    for (int n = 1; n <= NKNOWN; n++)
        if (KNOWN[n] == p) printf(L == n ? "  (= OEIS a(%d))" : "  (OEIS has this as a(%d): MISMATCH)", n);
    putchar('\n');
    printf("primes: %" PRIu64, p);
    u64 q = p;
    for (int i = 0; i <= L; i++) { q += gaps[i]; printf(" %" PRIu64, q); }
    putchar('\n');
#ifdef HAVE_GMP
    {
        mpz_t z;
        mpz_init_set_ui(z, p);
        bool ok = true;
        u64 r = p;
        for (int i = 0; i <= L && ok; i++) {
            mpz_nextprime(z, z);
            r += gaps[i];
            if (mpz_cmp_ui(z, r) != 0) ok = false;
        }
        mpz_clear(z);
        printf("GMP mpz_nextprime %s\n", ok ? "agrees" : "DISAGREES");
        if (!ok) return 1;
    }
#endif
    return 0;
}

static int cmd_bench(int argc, char **argv)
{
    u64 N = 1000000000000000ULL, span = 100000000000ULL, chunk = 10000000000ULL;
    int threads = default_threads(), npos = 0;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) threads = parse_int(argv[++i], 1, 1024, "threads");
        else if (!strcmp(argv[i], "-d") && i + 1 < argc) span = parse_num(argv[++i]);
        else if (!strcmp(argv[i], "-c") && i + 1 < argc) chunk = parse_num(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i + 1 < argc) g_sieve_kib = parse_int(argv[++i], 16, 65536, "sieve KiB");
        else if (argv[i][0] == '-' && !isdigit((unsigned char)argv[i][1])) usage();
        else if (npos++ == 0) N = parse_num(argv[i]);
        else usage();
    }
    if (chunk < 1000000) die("CHUNK must be at least 10^6");
    scan_stats_t st = run_scan(N, N + span, chunk, threads, NULL, 60, 0, true);
    double rate = (double)st.numbers / st.seconds;
    printf("bench: [%" PRIu64 ", %" PRIu64 ") with %d threads, chunk %" PRIu64 ": %.2f s\n",
           N, N + span, threads, chunk, st.seconds);
    printf("       %.4g numbers/s, %.4g primes/s (%" PRIu64 " primes)\n",
           rate, (double)st.primes / st.seconds, st.primes);
    printf("       time for a scan [0, X) at this rate (the sieve is faster below %.3g and slower above):\n", (double)N);
    static const double X[] = { 1e15, 2e15, 5e15, 1e16 };
    char b[32];
    for (size_t i = 0; i < sizeof X / sizeof *X; i++)
        printf("         X = %.0e: %s\n", X[i], fmt_dur(X[i] / rate, b, sizeof b));
    print_hist();
    return 0;
}

static int cmd_selftest(int argc, char **argv)
{
    u64 limit = 1000000000000ULL;
    int threads = default_threads(), npos = 0, fails = 0;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) threads = parse_int(argv[++i], 1, 1024, "threads");
        else if (argv[i][0] == '-' && !isdigit((unsigned char)argv[i][1])) usage();
        else if (npos++ == 0) limit = parse_num(argv[i]);
        else usage();
    }
#define CHECK(cond, ...) do { if (cond) printf("ok    " __VA_ARGS__); else { printf("FAIL  " __VA_ARGS__); fails++; } putchar('\n'); } while (0)

    CHECK(parse_num("2^54-32") == 18014398509481952ULL && parse_num("1e12") == 1000000000000ULL &&
          parse_num("1.29e14") == 129000000000000ULL && parse_num("10^15+7") == 1000000000000007ULL,
          "number parsing");
    CHECK(is_prime64(2) && is_prime64(3) && !is_prime64(1) && !is_prime64(9) && is_prime64(1000000007ULL) &&
          !is_prime64(3215031751ULL) && is_prime64(18446744073709551557ULL) && !is_prime64(18446744073709551555ULL) &&
          is_prime64(129001208165717ULL) && !is_prime64(129001208165719ULL - 2 + 4),
          "Miller-Rabin spot checks");
    for (int n = 1; n <= NKNOWN; n++) {
        u16 gaps[NMAX + 1];
        int L = run_length64(KNOWN[n], gaps, NMAX + 1);
        CHECK(L == n, "run length of a(%d) = %" PRIu64 " by next-prime search is %d", n, KNOWN[n], L);
    }

    /* chunk boundaries: odd chunk sizes against a brute-force reference */
    int save_report = g_report;
    g_report = NMAX;
    {
        const u64 lim = 200000000ULL;
        u64 ref[NMAX], refhist[NMAX];
        static const u64 chunks[] = { 1234567ULL, 7000001ULL, 200000000ULL };
        for (size_t ci = 0; ci < sizeof chunks / sizeof *chunks; ci++) {
            run_scan(0, lim, chunks[ci], threads, NULL, 60, 0, true);
            reference_first(S.end, ref, refhist);       /* S.end = lim rounded up to a chunk boundary */
            bool same = true, samehist = true;
            int nref = 0;
            u64 tot = 0;
            for (int n = 1; n < NMAX; n++) {
                if (ref[n]) nref = n;
                if (S.best[n].p != ref[n]) same = false;
                if (S.hist[n] != refhist[n]) samehist = false;
                tot += S.hist[n];
            }
            CHECK(same, "chunked scan (chunk %" PRIu64 ", %d threads) matches brute force on [0, %" PRIu64 "), n = 1..%d",
                  chunks[ci], threads, S.end, nref);
            CHECK(samehist && tot == S.primes && (S.end != 200000000ULL || S.primes == 11078937ULL),
                  "run-length histogram matches brute force and sums to pi(%" PRIu64 ") = %" PRIu64, S.end, S.primes);
        }
    }

    printf("sieve scan of [0, %" PRIu64 ") with %d threads ...\n", limit, threads);
    fflush(stdout);
    u64 chunk = limit >= 100000000000ULL ? 10000000000ULL : (limit >= 10000000000ULL ? 1000000000ULL : 100000000ULL);
    scan_stats_t st = run_scan(0, limit, chunk, threads, NULL, 60, 0, true);
    g_report = save_report;
    char b[32];
    printf("      %" PRIu64 " primes in %s (%.3g numbers/s)\n", st.primes, fmt_dur(st.seconds, b, sizeof b),
           st.seconds > 0 ? (double)st.numbers / st.seconds : 0);
    for (int n = 1; n <= NKNOWN; n++) {
        if (KNOWN[n] < S.end) {
            CHECK(S.best[n].p == KNOWN[n], "a(%d) = %" PRIu64 " (sieve found %" PRIu64 ")", n, KNOWN[n], S.best[n].p);
            if (n >= 2)
                CHECK(S.best[n].pi == A133697[n - 2], "pi(a(%d)) = %" PRIu64 " = A133697(%d)", n, S.best[n].pi, n - 2);
        }
        else
            CHECK(S.best[n].p == 0, "a(%d) = %" PRIu64 " is above the limit and was not found below it", n, KNOWN[n]);
    }
    for (int n = NKNOWN + 1; n < NMAX; n++) {
        if (!S.best[n].p) continue;
        if (S.end <= 130000000000000ULL)
            CHECK(false, "unexpected run length %d at %" PRIu64 " (Resta found none below 1.3e14)", n, S.best[n].p);
        else
            printf("NEW   a(%d) = %" PRIu64 " (not in OEIS; verify with 'a158939 verify')\n", n, S.best[n].p);
    }
    /* cross-check every sieve result with the independent next-prime search */
    {
        bool all = true;
        for (int n = 1; n < NMAX; n++) {
            if (!S.best[n].p) continue;
            u16 gaps[NMAX + 1];
            int L = run_length64(S.best[n].p, gaps, NMAX + 1);
            if (L != n || memcmp(gaps, S.best[n].gaps, (size_t)(n + 1) * sizeof(u16)) != 0) {
                all = false;
                printf("      sieve and next-prime search disagree at %" PRIu64 "\n", S.best[n].p);
            }
        }
        CHECK(all, "sieve gaps agree with next-prime search for every a(n) found");
    }
#undef CHECK
    printf(fails ? "SELFTEST FAILED (%d)\n" : "selftest passed\n", fails);
    return fails ? 1 : 0;
}

int main(int argc, char **argv)
{
    stderr_tty = isatty(2);
    setvbuf(stdout, NULL, _IOLBF, 0);
    if (argc < 2) usage();
    if (!strcmp(argv[1], "scan")) return cmd_scan(argc - 2, argv + 2);
    if (!strcmp(argv[1], "verify")) return cmd_verify(argc - 2, argv + 2);
    if (!strcmp(argv[1], "bench")) return cmd_bench(argc - 2, argv + 2);
    if (!strcmp(argv[1], "selftest")) return cmd_selftest(argc - 2, argv + 2);
    usage();
    return 2;
}
