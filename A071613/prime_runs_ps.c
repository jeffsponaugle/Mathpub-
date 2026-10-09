/*
 * prime_runs_ps.c -- Runs of CONSECUTIVE primes sharing a digit sum,
 *                    powered by the primesieve library.
 *
 * The prime stream comes from primesieve's iterator (state-of-the-art,
 * memory-light, generates primes in order very fast).  On top of it sits the
 * same machinery as before:
 *   - mod-(base-1) lazy reject: digit sum computed only when a run could extend
 *   - table-based digit sum for base 10
 *   - pthreads with overlap-and-ownership so output is identical for any -t
 *
 * New in this revision:
 *   - -S/-E accept scientific notation (3e5, 1.5e12) as well as K/M/B/T
 *   - --checkpoint FILE: periodic crash-safe state; rerunning resumes
 *   - SIGINT/SIGTERM: graceful stop, summary printed, checkpoint flushed
 *   - run-length histogram: how many maximal runs of each shorter length
 *
 * Build:
 *   cc -O3 -march=native -pthread -o prime_runs_ps prime_runs_ps.c -lprimesieve -lm
 *
 * Usage:
 *   ./prime_runs_ps -n 8 -L 3e9 -t 8
 *   ./prime_runs_ps -n 3 -L 60000 --all
 *   ./prime_runs_ps -n 9 -E 1e13 -t 16 --checkpoint run9.ck
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <sys/types.h>
#include <primesieve.h>

typedef uint64_t u64;
typedef unsigned long long ull;

static inline double mono(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

/* Parse a count with optional suffix: T=1e12, B=1e9, M=1e6, K=1e3.
 * Scientific / decimal notation is accepted too: "3e5", "1.5e12", "2.5B".
 * Plain integers are read exactly; '_' and ',' are ignored as separators. */
static u64 parse_num(const char *s) {
    char buf[80];
    size_t n = 0;
    for (const char *p = s; *p && n + 1 < sizeof buf; p++)
        if (*p != '_' && *p != ',') buf[n++] = *p;
    buf[n] = '\0';
    if (!n) return 0;

    u64 mult = 1;
    switch (buf[n - 1]) {
        case 'T': case 't': mult = 1000000000000ULL; buf[--n] = '\0'; break;
        case 'B': case 'b': mult = 1000000000ULL;    buf[--n] = '\0'; break;
        case 'M': case 'm': mult = 1000000ULL;       buf[--n] = '\0'; break;
        case 'K': case 'k': mult = 1000ULL;          buf[--n] = '\0'; break;
    }
    if (!n) return 0;

    int sci = 0;
    for (size_t i = 0; i < n; i++)
        if (buf[i] == '.' || buf[i] == 'e' || buf[i] == 'E') { sci = 1; break; }

    if (sci) {                                   /* 3e5, 1.5e12, 2.5B, ... */
        double d = strtod(buf, NULL) * (double)mult;
        if (d <= 0) return 0;
        if (d >= 18446744073709549568.0) return UINT64_MAX;
        return (u64)(d + 0.5);
    }
    return strtoull(buf, NULL, 10) * mult;       /* exact integer path */
}

static int  g_L, g_base, g_find_all, g_quiet, g_emirps;
static int  g_stream_only, g_sieve_kb;
static u64  g_M, g_start, g_end;

#define OVERLAP (1u << 20)        /* boundary margin (integers) */
#define MAXH    65                /* run-length histogram buckets (last = ">=") */

/* ---- digit sum ---------------------------------------------------------
 * dsw[i] is the digit sum of i in base g_base, for i < g_WIN, where g_WIN is
 * the largest power of the base that fits in 1 MB.  Because the primes arrive
 * in increasing order we never compute a digit sum from scratch: the low
 * window slides by the prime gap, and the digit sum of everything above the
 * window only changes when the window carries -- about once in 25,000 primes
 * at base 10.  So the steady-state cost is one L2-resident byte load and an
 * add, and the mod-(base-1) pre-filter it replaces is gone entirely. */
static uint8_t *dsw;
static u64 g_WIN = 1;
static inline int dsum(u64 n) {
    int s = 0; while (n) { s += dsw[n % g_WIN]; n /= g_WIN; } return s;
}

/* ---- deterministic Miller-Rabin for all 64-bit n ---------------------- */
static inline u64 mulmod(u64 a, u64 b, u64 m) {
    return (u64)((__uint128_t)a * b % m);
}
static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
static int is_prime_u64(u64 n) {
    if (n < 2) return 0;
    static const int sm[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) {
        if (n == (u64)sm[i]) return 1;
        if (n % (u64)sm[i] == 0) return 0;
    }
    u64 d = n - 1; int r = 0;
    while ((d & 1) == 0) { d >>= 1; r++; }
    for (int i = 0; i < 12; i++) {                 /* these bases are proven  */
        u64 x = powmod((u64)sm[i], d, n);          /* sufficient for n < 2^64 */
        if (x == 1 || x == n - 1) continue;
        int composite = 1;
        for (int j = 1; j < r; j++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { composite = 0; break; }
        }
        if (composite) return 0;
    }
    return 1;
}

/* Reverse the digits of n in base g_base; returns 0 on overflow. */
static u64 revdigits(u64 n) {
    u64 b = (u64)g_base, r = 0;
    while (n) {
        u64 dg = n % b; n /= b;
        if (r > (UINT64_MAX - dg) / b) return 0;   /* would overflow */
        r = r * b + dg;
    }
    return r;
}

/* An emirp: prime whose digit-reversal is a DIFFERENT prime.
 * (Palindromic primes are excluded, per the usual convention.)
 * n is already known to be prime when this is called. */
static int is_emirp(u64 n) {
    u64 r = revdigits(n);
    if (r == 0 || r == n) return 0;                /* overflow or palindrome */
    return is_prime_u64(r);
}

/* ---- per-thread tracker + results ------------------------------------ */
typedef struct {
    u64 lo, hi, A, B, A0;   /* own [lo,hi); scan [A,B); A0 = scan start of a
                               fresh (non-resumed) run, used for progress */
    u64 stat_lo;            /* count maximal runs whose head is >= this */
    u64 emit0;              /* on resume: only emit windows starting > this */
    int tid;
    u64 *first; int *dsv; u64 *members; size_t count, cap;
    int len; int ds; u64 *ring; size_t rpos, rmask; u64 head; int done;
    u64 prev, low; int hids;      /* sliding digit-sum window */
    signed char *eflag;  /* per ring slot: -1 unknown, 0 not emirp, 1 emirp */
    u64 *win;            /* scratch: the window being emitted */
    /* thread-local statistics, merged into the globals at publish points */
    u64 lhist[MAXH], lfirst[MAXH];
    u64 lruns, lmax, lmax_at, lprimes;
} Targ;

/* ---- shared state (all touched under g_pmx) -------------------------- */
static pthread_mutex_t g_pmx = PTHREAD_MUTEX_INITIALIZER;
static u64 g_primes_done;                    /* this session, for the rate  */
static u64 g_primes_total;                   /* including resumed sessions  */
static double g_t_start, g_t_last, g_rate;   /* timing + windowed rate */
static u64 g_count_last;
static int g_nthreads;
static u64 *g_pos;                           /* per-thread value-progress   */
static u64 g_total_span;                     /* sum of (B - A0) over threads */

static u64 g_hist[MAXH], g_hfirst[MAXH];     /* run-length histogram        */
static u64 g_runs_total, g_maxrun, g_maxrun_at;
static u64 g_found;                          /* windows recorded, all runs  */

static const char *g_ckpath;                 /* checkpoint file, or NULL    */
static char  g_runpath[4096];                /* "<ckpath>.runs"             */
static FILE *g_runfp;                        /* append-only results stream  */
static double g_ck_secs = 30.0, g_ck_last;
static u64 *g_ck_resume, *g_ck_emit;
static volatile sig_atomic_t g_stop;         /* set by SIGINT / SIGTERM     */

static void on_signal(int sig) { (void)sig; g_stop = 1; }

/* ---- statistics bookkeeping ------------------------------------------ */
/* Called just before a maximal run is torn down.  Ownership is by the run's
 * HEAD, exactly like the window ownership rule, so overlapping thread ranges
 * (and re-scanned regions after a resume) never double-count. */
static inline void note_run_end(Targ *t) {
    if (t->len < 1) return;
    if (t->head < t->stat_lo || t->head >= t->hi) return;
    int b = t->len < MAXH ? t->len : MAXH - 1;
    t->lhist[b]++;
    if (!t->lfirst[b] || t->head < t->lfirst[b]) t->lfirst[b] = t->head;
    t->lruns++;
    /* Primes are tallied per completed run rather than per prime, so the count
     * is not inflated by overlap or by re-scanning after a resume, and the
     * histogram's sum of (length x count) reproduces it exactly. */
    t->lprimes += (u64)t->len;
    if ((u64)t->len > t->lmax) { t->lmax = (u64)t->len; t->lmax_at = t->head; }
}

static void publish_locked(Targ *t, u64 resume) {
    g_ck_resume[t->tid] = resume;
    for (int i = 0; i < MAXH; i++) {
        g_hist[i] += t->lhist[i]; t->lhist[i] = 0;
        if (t->lfirst[i] && (!g_hfirst[i] || t->lfirst[i] < g_hfirst[i]))
            g_hfirst[i] = t->lfirst[i];
        t->lfirst[i] = 0;
    }
    g_runs_total += t->lruns; t->lruns = 0;
    /* ties break on the smaller starting prime, so the record does not depend
     * on which thread happened to publish first */
    if (t->lmax && (t->lmax > g_maxrun ||
                    (t->lmax == g_maxrun &&
                     (!g_maxrun_at || t->lmax_at < g_maxrun_at)))) {
        g_maxrun = t->lmax; g_maxrun_at = t->lmax_at;
    }
    t->lmax = 0;
    g_primes_done += t->lprimes; g_primes_total += t->lprimes; t->lprimes = 0;
}

/* ---- checkpoint ------------------------------------------------------- */
/* Written atomically (temp file + rename).  The results stream is flushed
 * first and its byte offset recorded, so a resume truncates away anything
 * written after the snapshot -- the two files always agree. */
static void ckpt_write_locked(Targ *ta) {
    if (!g_ckpath) return;
    char tmp[4200];
    snprintf(tmp, sizeof tmp, "%s.tmp", g_ckpath);

    long long off = 0;
    if (g_runfp) { fflush(g_runfp); off = (long long)ftello(g_runfp); }

    FILE *f = fopen(tmp, "w");
    if (!f) { g_ck_last = mono(); return; }
    fprintf(f, "#prime_runs_ps checkpoint v1\n");
    fprintf(f, "length %d\nbase %d\nemirps %d\nfindall %d\n",
            g_L, g_base, g_emirps, g_find_all);
    fprintf(f, "start %llu\nend %llu\n", (ull)g_start, (ull)g_end);
    fprintf(f, "threads %d\nrunsoff %lld\nfound %llu\n",
            g_nthreads, off, (ull)g_found);
    fprintf(f, "primes %llu\nruns %llu\nmaxrun %llu %llu\n",
            (ull)g_primes_total, (ull)g_runs_total,
            (ull)g_maxrun, (ull)g_maxrun_at);
    for (int i = 0; i < g_nthreads; i++)
        fprintf(f, "T %d %llu %llu %llu %llu\n", i,
                (ull)ta[i].lo, (ull)ta[i].hi,
                (ull)g_ck_resume[i], (ull)g_ck_emit[i]);
    for (int i = 0; i < MAXH; i++)
        if (g_hist[i])
            fprintf(f, "H %d %llu %llu\n", i, (ull)g_hist[i], (ull)g_hfirst[i]);
    fflush(f);
    fsync(fileno(f));
    fclose(f);
    rename(tmp, g_ckpath);
    g_ck_last = mono();
}

/* ---- emitting a hit --------------------------------------------------- */
static void store_local(Targ *t) {
    if (t->count == t->cap) {
        t->cap = t->cap ? t->cap * 2 : 16;
        t->first   = realloc(t->first,   t->cap * sizeof(u64));
        t->dsv     = realloc(t->dsv,     t->cap * sizeof(int));
        t->members = realloc(t->members, t->cap * (size_t)g_L * sizeof(u64));
    }
    u64 *dst = &t->members[t->count * (size_t)g_L];
    memcpy(dst, t->win, (size_t)g_L * sizeof(u64));
    t->first[t->count] = dst[0];
    t->dsv[t->count]   = t->ds;
    t->count++;
}

static void record(Targ *t) {
    /* the current window is the last g_L primes pushed into the ring */
    size_t base = t->rpos - (size_t)g_L;
    for (int j = 0; j < g_L; j++)
        t->win[j] = t->ring[(base + (size_t)j) & t->rmask];

    if (g_runfp) {
        /* The write and the "emitted up to" publication happen in the same
         * critical section, so a checkpoint can never disagree with the file. */
        pthread_mutex_lock(&g_pmx);
        fprintf(g_runfp, "%d", t->ds);
        for (int j = 0; j < g_L; j++) fprintf(g_runfp, " %llu", (ull)t->win[j]);
        fputc('\n', g_runfp);
        g_ck_emit[t->tid] = t->win[0];
        g_found++;
        pthread_mutex_unlock(&g_pmx);
    } else {
        store_local(t);
        pthread_mutex_lock(&g_pmx);
        g_found++;
        pthread_mutex_unlock(&g_pmx);
    }
}

static inline void feed(Targ *t, u64 q) {
    if (t->done) return;
    /* Slide the low window forward by the gap.  t->prev starts at 0, so the
     * first prime of a thread naturally takes the carry path and initializes
     * both halves. */
    t->low += q - t->prev; t->prev = q;
    if (t->low >= g_WIN) { t->low %= g_WIN; t->hids = dsum(q / g_WIN); }
    int ds = t->hids + (int)dsw[t->low];

    int extended = 0;
    if (t->len >= 1 && ds == t->ds) {
        t->eflag[t->rpos & t->rmask] = -1;
        t->ring[t->rpos & t->rmask] = q; t->rpos++;
        t->len++; extended = 1;
    }
    if (!extended) {                                   /* start a new run */
        note_run_end(t);                               /* stats for the old one */
        t->len = 1; t->ds = ds; t->head = q;
        t->rpos = 0;
        t->ring[0] = q; t->eflag[0] = -1; t->rpos = 1;
    }

    /* Every window of g_L consecutive primes inside a maximal run counts, so
     * emit once for each prime that brings the run to length >= g_L.  A
     * maximal run of length M therefore yields M - g_L + 1 windows. */
    if (t->len >= g_L) {
        if (g_emirps) {
            /* Only now -- once L primes actually share a digit sum -- is the
             * expensive reversal-primality test worth doing.  Flags are cached
             * per ring slot, so overlapping windows never retest a prime, and
             * we scan newest-first since that slot is the only new one. */
            size_t wbase = t->rpos - (size_t)g_L;
            for (int j = g_L - 1; j >= 0; j--) {
                size_t sl = (wbase + (size_t)j) & t->rmask;
                if (t->eflag[sl] < 0)
                    t->eflag[sl] = is_emirp(t->ring[sl]) ? 1 : 0;
                if (!t->eflag[sl]) return;             /* not all emirps */
            }
        }
        u64 f = t->ring[(t->rpos - (size_t)g_L) & t->rmask];
        if (f >= t->lo && f < t->hi && f > t->emit0) {
            record(t);
            if (!g_find_all) t->done = 1;
        }
    }
}

/* ---- progress -------------------------------------------------------- */
static void fmt_hms(double secs, char *out, size_t n) {
    if (secs < 0) secs = 0;
    if (secs > 3.15e9) { snprintf(out, n, ">100y"); return; }
    long s = (long)(secs + 0.5), h = s/3600; s %= 3600; long m = s/60; s %= 60;
    if (h > 0) snprintf(out, n, "%ld:%02ld:%02ld", h, m, s);
    else       snprintf(out, n, "%ld:%02ld", m, s);
}

/* ---- worker: stream primes from primesieve over [A,B) ---------------- */
static Targ *g_ta;                            /* for ckpt_write_locked()    */

static void *worker(void *arg) {
    Targ *t = (Targ *)arg;
    t->len = 0; t->ds = -1; t->done = 0; t->rpos = 0; t->head = 0;
    t->prev = 0; t->low = 0; t->hids = 0;
    size_t rs = 1; while (rs < (size_t)g_L) rs <<= 1;   /* power-of-two ring */
    t->rmask = rs - 1;
    t->ring  = (u64 *)malloc(rs * sizeof(u64));
    t->win   = (u64 *)malloc((size_t)g_L * sizeof(u64));
    t->eflag = (signed char *)malloc(rs);
    memset(t->eflag, -1, rs);

    primesieve_iterator it;
    primesieve_init(&it);
    primesieve_jump_to(&it, t->A, t->B);     /* first prime >= A */

    u64 local = 0, p, lastp = t->A;
    int interrupted = 0;
    while ((p = primesieve_next_prime(&it)) < t->B) {
        lastp = p;
        if (!g_stream_only) feed(t, p);
        if (t->done) break;
        if ((++local & 0xFFFFu) == 0) {
            if (g_stop) { interrupted = 1; break; }
            if ((local & ((1u << 22) - 1)) == 0) {
                pthread_mutex_lock(&g_pmx);
                u64 mark = t->len ? t->head : p;
                publish_locked(t, mark);
                g_pos[t->tid] = p - t->A0;
                double now = mono(), since = now - g_t_last;
                if (since >= 0.2) {                   /* refresh windowed rate */
                    g_rate = (double)(g_primes_done - g_count_last) / since;
                    g_t_last = now; g_count_last = g_primes_done;
                }
                if (!g_quiet) {
                    double elapsed = now - g_t_start;
                    double avg = elapsed > 0 ? (double)g_primes_done / elapsed : 0;
                    double cur = g_rate > 0 ? g_rate : avg;
                    u64 covered = 0;
                    for (int i = 0; i < g_nthreads; i++) covered += g_pos[i];
                    double frac = g_total_span ? (double)covered / (double)g_total_span : 0;
                    if (frac > 1) frac = 1;
                    char el[32], eta[32];
                    fmt_hms(elapsed, el, sizeof el);
                    if (frac > 1e-6) fmt_hms(elapsed * (1.0 - frac) / frac, eta, sizeof eta);
                    else             snprintf(eta, sizeof eta, "--:--");
                    fprintf(stderr,
                        "\r%5.1f%% | %.1f M/s | elapsed %s | ETA %s | found %llu   ",
                        frac * 100.0, cur / 1e6, el, eta, (ull)g_found);
                    fflush(stderr);
                }
                if (g_ckpath && now - g_ck_last >= g_ck_secs)
                    ckpt_write_locked(g_ta);
                pthread_mutex_unlock(&g_pmx);
            }
        }
    }
    primesieve_free_iterator(&it);

    /* Final publish.  If we ran to the end of our range there is nothing left
     * to redo, so the resume mark is B; if we were interrupted mid-run the
     * mark is the head of that run, which makes the rescan reproduce both the
     * run state and its statistics contribution exactly. */
    if (!interrupted) note_run_end(t);
    pthread_mutex_lock(&g_pmx);
    publish_locked(t, interrupted ? (t->len ? t->head : lastp) : t->B);
    if (!interrupted) g_pos[t->tid] = t->B - t->A0;
    pthread_mutex_unlock(&g_pmx);
    return NULL;
}

static u64 align30(u64 x) { return x - (x % 30); }

typedef struct { u64 first; int ds; const u64 *mem; } Rec;
static int reccmp(const void *a, const void *b) {
    u64 x = ((const Rec *)a)->first, y = ((const Rec *)b)->first;
    return (x > y) - (x < y);
}

/* ---- checkpoint loading ---------------------------------------------- */
typedef struct { u64 lo, hi, resume, emit; } CkT;

static int ckpt_load(CkT **out, int *nt, long long *runsoff) {
    FILE *f = fopen(g_ckpath, "r");
    if (!f) return 0;
    char line[512];
    int L = -1, base = -1, em = -1, fa = -1, n = -1, seen = 0;
    u64 st = 0, en = 0;
    CkT *ts = NULL;
    *runsoff = 0;
    while (fgets(line, sizeof line, f)) {
        u64 a, b, c, d; int i;
        if      (sscanf(line, "length %d", &L) == 1) seen++;
        else if (sscanf(line, "base %d", &base) == 1) seen++;
        else if (sscanf(line, "emirps %d", &em) == 1) seen++;
        else if (sscanf(line, "findall %d", &fa) == 1) seen++;
        else if (sscanf(line, "start %llu", (ull *)&st) == 1) seen++;
        else if (sscanf(line, "end %llu", (ull *)&en) == 1) seen++;
        else if (sscanf(line, "threads %d", &n) == 1) {
            if (n < 1 || n > 100000) { fclose(f); free(ts); return -1; }
            ts = (CkT *)calloc((size_t)n, sizeof(CkT));
        }
        else if (sscanf(line, "runsoff %lld", runsoff) == 1) ;
        else if (sscanf(line, "found %llu", (ull *)&a) == 1) g_found = a;
        else if (sscanf(line, "primes %llu", (ull *)&a) == 1) g_primes_total = a;
        else if (sscanf(line, "runs %llu", (ull *)&a) == 1) g_runs_total = a;
        else if (sscanf(line, "maxrun %llu %llu", (ull *)&a, (ull *)&b) == 2) {
            g_maxrun = a; g_maxrun_at = b;
        }
        else if (sscanf(line, "T %d %llu %llu %llu %llu",
                        &i, (ull *)&a, (ull *)&b, (ull *)&c, (ull *)&d) == 5) {
            if (!ts || i < 0 || i >= n) { fclose(f); free(ts); return -1; }
            ts[i].lo = a; ts[i].hi = b; ts[i].resume = c; ts[i].emit = d;
        }
        else if (sscanf(line, "H %d %llu %llu",
                        &i, (ull *)&a, (ull *)&b) == 3) {
            if (i >= 0 && i < MAXH) { g_hist[i] = a; g_hfirst[i] = b; }
        }
    }
    fclose(f);
    if (seen < 6 || !ts) { free(ts); return -1; }
    if (L != g_L || base != g_base || em != g_emirps || fa != g_find_all ||
        st != g_start || en != g_end) {
        fprintf(stderr,
            "checkpoint '%s' was made with different parameters\n"
            "  file: -n %d -b %d -S %llu -E %llu%s%s\n"
            "  now : -n %d -b %d -S %llu -E %llu%s%s\n"
            "Use a different --checkpoint name, or --restart to overwrite it.\n",
            g_ckpath, L, base, (ull)st, (ull)en,
            em ? " --emirps" : "", fa ? " --all" : "",
            g_L, g_base, (ull)g_start, (ull)g_end,
            g_emirps ? " --emirps" : "", g_find_all ? " --all" : "");
        free(ts);
        return -1;
    }
    *out = ts; *nt = n;
    return 1;
}

/* ---- results file -> in-memory records ------------------------------- */
static Rec *load_runs_file(size_t *total_out, u64 **mem_out) {
    *total_out = 0; *mem_out = NULL;
    FILE *f = fopen(g_runpath, "r");
    if (!f) return NULL;
    size_t lines = 0;
    int c, prev = '\n';
    while ((c = fgetc(f)) != EOF) { if (c == '\n') lines++; prev = c; }
    if (prev != '\n' && lines) lines++;
    if (!lines) { fclose(f); return NULL; }
    if (lines > 20000000UL) {                 /* refuse to eat all of RAM */
        fclose(f);
        *total_out = lines;
        return NULL;
    }
    rewind(f);
    Rec *recs = (Rec *)malloc(lines * sizeof(Rec));
    u64 *mem  = (u64 *)malloc(lines * (size_t)g_L * sizeof(u64));
    if (!recs || !mem) { fclose(f); free(recs); free(mem); return NULL; }
    char *buf = NULL; size_t bcap = 0; ssize_t got; size_t k = 0;
    while (k < lines && (got = getline(&buf, &bcap, f)) > 0) {
        char *p = buf, *e;
        long ds = strtol(p, &e, 10);
        if (e == p) continue;
        p = e;
        u64 *dst = &mem[k * (size_t)g_L];
        int j = 0;
        for (; j < g_L; j++) {
            dst[j] = strtoull(p, &e, 10);
            if (e == p) break;
            p = e;
        }
        if (j < g_L) continue;
        recs[k].ds = (int)ds; recs[k].first = dst[0]; recs[k].mem = dst;
        k++;
    }
    free(buf);
    fclose(f);
    *total_out = k; *mem_out = mem;
    return recs;
}

/* ---- reporting -------------------------------------------------------- */
static void print_stats(void) {
    printf("\nMaximal runs of consecutive primes with equal digit sum (base %d)\n",
           g_base);
    printf("  primes counted : %llu\n", (ull)g_primes_total);
    printf("  maximal runs   : %llu\n", (ull)g_runs_total);
    if (g_maxrun)
        printf("  longest run    : %llu, starting at %llu\n",
               (ull)g_maxrun, (ull)g_maxrun_at);
    if (!g_find_all && g_found)
        printf("  note: without --all each thread stops at its first hit, so\n"
               "        these counts cover only the range scanned before that.\n");
    if (g_emirps)
        printf("  note: the histogram is over digit-sum runs; the emirp test is\n"
               "        applied only to windows that already reach length %d.\n", g_L);
    printf("\n  length        count      share    prev/this   smallest run starts at\n");
    u64 prev = 0;
    for (int i = 1; i < MAXH; i++) {
        if (!g_hist[i]) { if (i > (int)g_maxrun) break; prev = 0; continue; }
        double share = g_runs_total ? 100.0 * (double)g_hist[i] / (double)g_runs_total : 0;
        printf("  %s%-4d %12llu %8.4f%%", i == MAXH - 1 ? ">=" : "  ", i,
               (ull)g_hist[i], share);
        if (prev) printf("   %9.3f", (double)prev / (double)g_hist[i]);
        else      printf("   %9s", "--");
        printf("   %llu\n", (ull)g_hfirst[i]);
        prev = g_hist[i];
    }
}

static void usage(const char *prog) {
    fprintf(stderr,
"usage: %s -n LEN [options]\n"
"  -n, --length N     run length: N consecutive primes sharing a digit sum\n"
"  -S, --start V      range start (default 0)\n"
"  -E, --end V        range end, exclusive (default 1e9); -L is an alias\n"
"                     V accepts 1234567, 1e10, 1.5e12, 2.5B, 12K, 1_000_000\n"
"  -b, --base B       digit-sum base, 2..36 (default 10)\n"
"  -t, --threads T    worker threads (default 1)\n"
"      --all          report every run, not just the first\n"
"      --emirps       require every member of the run to be an emirp\n"
"  -c, --checkpoint F save state to F (and hits to F.runs); rerun to resume\n"
"      --ckpt-secs S  checkpoint interval in seconds (default 30)\n"
"      --restart      discard an existing checkpoint instead of resuming\n"
"      --no-stats     skip the run-length histogram\n"
"      --sieve-size K primesieve segment size in KB (default: auto)\n"
"      --stream-only  generate primes but skip the search (timing baseline)\n"
"      --quiet        no progress line\n"
"SIGINT/SIGTERM stop cleanly, print a summary, and flush the checkpoint.\n",
        prog);
}

int main(int argc, char **argv) {
    g_L = 0; g_base = 10; g_find_all = 0; g_quiet = 0; g_emirps = 0;
    g_start = 0; g_end = 1000000000ULL;
    int nthreads = 1, restart = 0, nostats = 0;
    for (int i = 1; i < argc; i++) {
        if      ((!strcmp(argv[i],"--length") ||!strcmp(argv[i],"-n")) && i+1<argc) g_L=atoi(argv[++i]);
        else if ((!strcmp(argv[i],"--start")  ||!strcmp(argv[i],"-S")) && i+1<argc) g_start=parse_num(argv[++i]);
        else if ((!strcmp(argv[i],"--end")    ||!strcmp(argv[i],"-E")) && i+1<argc) g_end=parse_num(argv[++i]);
        else if ((!strcmp(argv[i],"--limit")  ||!strcmp(argv[i],"-L")) && i+1<argc) g_end=parse_num(argv[++i]); /* alias for -E */
        else if ((!strcmp(argv[i],"--base")   ||!strcmp(argv[i],"-b")) && i+1<argc) g_base=atoi(argv[++i]);
        else if ((!strcmp(argv[i],"--threads")||!strcmp(argv[i],"-t")) && i+1<argc) nthreads=atoi(argv[++i]);
        else if ((!strcmp(argv[i],"--checkpoint")||!strcmp(argv[i],"-c")) && i+1<argc) g_ckpath=argv[++i];
        else if ( !strcmp(argv[i],"--ckpt-secs") && i+1<argc) g_ck_secs=atof(argv[++i]);
        else if ( !strcmp(argv[i],"--restart")) restart=1;
        else if ( !strcmp(argv[i],"--no-stats")) nostats=1;
        else if ( !strcmp(argv[i],"--all"))    g_find_all=1;
        else if ( !strcmp(argv[i],"--emirps")) g_emirps=1;
        else if ( !strcmp(argv[i],"--quiet")) g_quiet=1;
        else if ( !strcmp(argv[i],"--sieve-size") && i+1<argc) g_sieve_kb=atoi(argv[++i]);
        else if ( !strcmp(argv[i],"--stream-only")) g_stream_only=1;
        else if ( !strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")) { usage(argv[0]); return 0; }
        else { fprintf(stderr,"bad arg: %s\n", argv[i]); usage(argv[0]); return 2; }
    }
    if (g_L < 1) { fprintf(stderr,"need --length >= 1\n"); usage(argv[0]); return 2; }
    if (g_base < 2 || g_base > 36) { fprintf(stderr,"base 2..36\n"); return 2; }
    if (g_end <= g_start) { fprintf(stderr,"need end (-E) > start (-S)\n"); return 2; }
    if (nthreads < 1) nthreads = 1;
    if (g_ck_secs < 1) g_ck_secs = 1;
    g_M = (u64)g_base - 1; if (g_M == 0) g_M = 1;
    /* largest power of the base that fits in 1 MB; >= 35937 for every base,
     * far larger than any prime gap below 2^64, so one carry per slide */
    while (g_WIN * (u64)g_base <= (1u << 20)) g_WIN *= (u64)g_base;
    dsw = (uint8_t *)malloc((size_t)g_WIN);
    if (!dsw) { fprintf(stderr, "out of memory\n"); return 2; }
    dsw[0] = 0;
    for (u64 i = 1; i < g_WIN; i++)
        dsw[i] = (uint8_t)(dsw[i / (u64)g_base] + i % (u64)g_base);
    if (g_sieve_kb > 0) primesieve_set_sieve_size(g_sieve_kb);

    /* ---- partition: fresh, or restored from a checkpoint -------------- */
    CkT *ck = NULL; int cknt = 0; long long runsoff = 0; int resumed = 0;
    if (g_ckpath) {
        snprintf(g_runpath, sizeof g_runpath, "%s.runs", g_ckpath);
        if (restart) { remove(g_ckpath); remove(g_runpath); }
        else {
            int r = ckpt_load(&ck, &cknt, &runsoff);
            if (r < 0) return 3;
            resumed = (r > 0);
        }
    }

    Targ *ta;
    if (resumed) {
        nthreads = cknt;
        ta = (Targ *)calloc(nthreads, sizeof(Targ));
        for (int i = 0; i < nthreads; i++) {
            ta[i].lo = ck[i].lo; ta[i].hi = ck[i].hi; ta[i].tid = i;
            u64 rs = ck[i].resume;
            if (rs < ta[i].lo) rs = ta[i].lo;
            if (rs > ta[i].hi) rs = ta[i].hi;
            ta[i].stat_lo = rs;
            ta[i].emit0   = ck[i].emit;
            ta[i].A  = (rs > OVERLAP) ? align30(rs - OVERLAP) : 0;
            ta[i].A0 = (ta[i].lo > OVERLAP) ? align30(ta[i].lo - OVERLAP) : 0;
            ta[i].B  = (ta[i].hi + OVERLAP < g_end) ? ta[i].hi + OVERLAP : g_end;
        }
        free(ck);
    } else {
        u64 span = g_end - g_start;
        while (nthreads > 1 && span / (u64)nthreads < (u64)OVERLAP * 8) nthreads--;
        ta = (Targ *)calloc(nthreads, sizeof(Targ));
        u64 chunk = align30(span / (u64)nthreads); if (chunk < 30) chunk = 30;
        for (int i = 0; i < nthreads; i++) {
            u64 lo = g_start + (u64)i * chunk;
            u64 hi = (i == nthreads - 1) ? g_end : g_start + (u64)(i + 1) * chunk;
            ta[i].lo = lo; ta[i].hi = hi; ta[i].tid = i;
            ta[i].stat_lo = lo; ta[i].emit0 = 0;
            ta[i].A  = (lo > OVERLAP) ? align30(lo - OVERLAP) : 0;
            ta[i].A0 = ta[i].A;
            ta[i].B  = (hi + OVERLAP < g_end) ? hi + OVERLAP : g_end;
        }
    }
    g_nthreads = nthreads;
    g_ta = ta;
    g_total_span = 0;
    for (int i = 0; i < nthreads; i++) g_total_span += ta[i].B - ta[i].A0;
    g_pos      = (u64 *)calloc(nthreads, sizeof(u64));
    g_ck_resume = (u64 *)calloc(nthreads, sizeof(u64));
    g_ck_emit   = (u64 *)calloc(nthreads, sizeof(u64));
    for (int i = 0; i < nthreads; i++) {
        g_ck_resume[i] = ta[i].stat_lo;
        g_ck_emit[i]   = ta[i].emit0;
        g_pos[i]       = ta[i].A - ta[i].A0;      /* already-covered span */
    }

    /* ---- results stream ------------------------------------------------ */
    if (g_ckpath) {
        if (resumed) {
            if (truncate(g_runpath, (off_t)runsoff) != 0 && runsoff)
                fprintf(stderr, "warning: could not truncate %s\n", g_runpath);
            g_runfp = fopen(g_runpath, "a");
        } else {
            g_runfp = fopen(g_runpath, "w");
        }
        if (!g_runfp) {
            fprintf(stderr, "cannot open %s for results\n", g_runpath);
            return 4;
        }
        if (!g_quiet) {
            if (resumed)
                fprintf(stderr, "resuming from %s (%llu primes done, %llu found)\n",
                        g_ckpath, (ull)g_primes_total, (ull)g_found);
            else
                fprintf(stderr, "checkpointing to %s\n", g_ckpath);
        }
    }

    int run_it = !(resumed && !g_find_all && g_found);
    if (!run_it && !g_quiet)
        fprintf(stderr, "checkpoint already holds a hit; nothing left to do\n");

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_signal;
    sigaction(SIGINT,  &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    g_t_start = g_t_last = mono();
    g_count_last = 0; g_rate = 0; g_ck_last = g_t_start;

    pthread_t *th = (pthread_t *)malloc(nthreads * sizeof(pthread_t));
    if (run_it) {
        for (int i = 0; i < nthreads; i++) pthread_create(&th[i], NULL, worker, &ta[i]);
        for (int i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
    }

    pthread_mutex_lock(&g_pmx);
    if (g_ckpath) ckpt_write_locked(ta);
    pthread_mutex_unlock(&g_pmx);
    if (g_runfp) { fflush(g_runfp); fclose(g_runfp); g_runfp = NULL; }

    if (!g_quiet && run_it) {
        char el[32]; fmt_hms(mono() - g_t_start, el, sizeof el);
        fprintf(stderr, "\r%s in %s%40s\n",
                g_stop ? "interrupted" : "completed", el, "");
    }

    /* ---- collect and report -------------------------------------------- */
    size_t total = 0;
    Rec *recs = NULL; u64 *filemem = NULL;
    if (g_ckpath) {
        recs = load_runs_file(&total, &filemem);
        if (!recs && total) {
            printf("%llu runs recorded in %s (too many to sort in memory).\n",
                   (ull)total, g_runpath);
            total = 0;
        }
    } else {
        for (int i = 0; i < nthreads; i++) total += ta[i].count;
        recs = (Rec *)malloc((total ? total : 1) * sizeof(Rec));
        size_t k = 0;
        for (int i = 0; i < nthreads; i++)
            for (size_t j = 0; j < ta[i].count; j++)
                recs[k++] = (Rec){ ta[i].first[j], ta[i].dsv[j],
                                   &ta[i].members[j*(size_t)g_L] };
    }
    if (recs) qsort(recs, total, sizeof(Rec), reccmp);

    if (g_stop)
        printf("** Interrupted at the user's request -- results below are partial. **\n");

    if (total == 0) {
        printf("No run of %d consecutive %sprimes with equal digit sum "
               "(base %d) in [%llu, %llu).\n",
               g_L, g_emirps ? "emirp " : "", g_base, (ull)g_start, (ull)g_end);
    } else {
        size_t show = g_find_all ? total : 1;
        printf("%s of %d consecutive %sprimes sharing a digit sum (base %d), "
               "in [%llu, %llu):\n",
               g_find_all ? "All runs" : "First run", g_L,
               g_emirps ? "emirp " : "", g_base,
               (ull)g_start, (ull)g_end);
        int idxw = 1; for (size_t tt = total; tt >= 10; tt /= 10) idxw++;
        for (size_t r = 0; r < show; r++) {
            if (g_find_all) printf("  #%*zu  digit sum %d: [", idxw, r + 1, recs[r].ds);
            else            printf("  digit sum %d: [", recs[r].ds);
            for (int j = 0; j < g_L; j++)
                printf("%llu%s", (ull)recs[r].mem[j], j+1 < g_L ? ", " : "");
            printf("]\n");
        }
        if (g_find_all)
            printf("Total sequences discovered: %zu\n", total);
    }

    if (!nostats) print_stats();

    if (g_stop && g_ckpath)
        printf("\nState saved in %s -- rerun the same command line to continue.\n",
               g_ckpath);

    free(filemem);
    free(recs);
    return g_stop ? 130 : 0;
}
