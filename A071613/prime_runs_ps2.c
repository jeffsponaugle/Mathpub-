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
 * Build:
 *   cc -O3 -march=native -pthread -o prime_runs_ps prime_runs_ps.c -lprimesieve -lm
 *
 * Usage:
 *   ./prime_runs_ps -n 8 -L 3000000000 -t 8
 *   ./prime_runs_ps -n 3 -L 60000 --all
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <time.h>
#include <primesieve.h>

typedef uint64_t u64;
typedef unsigned long long ull;

static inline double mono(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

/* Parse a count with optional suffix: T=1e12, B=1e9, M=1e6, K=1e3.
 * With a suffix the numeric part may be fractional ("1.5B" -> 1500000000).
 * Without a suffix it is read as an exact integer. */
static u64 parse_num(const char *s) {
    size_t n = strlen(s);
    if (!n) return 0;
    double mult = 0;
    switch (s[n - 1]) {
        case 'T': case 't': mult = 1e12; break;
        case 'B': case 'b': mult = 1e9;  break;
        case 'M': case 'm': mult = 1e6;  break;
        case 'K': case 'k': mult = 1e3;  break;
    }
    if (mult > 0) return (u64)(strtod(s, NULL) * mult + 0.5);
    return strtoull(s, NULL, 10);
}

static int  g_L, g_base, g_find_all, g_quiet;
static u64  g_M, g_start, g_end;

#define OVERLAP (1u << 20)        /* boundary margin (integers) */

/* ---- digit sum: base-10 table, generic fallback ---------------------- */
static int ds5[100000];
static inline int dsum(u64 n) {
    if (g_base == 10) { int s = 0; while (n) { s += ds5[n % 100000]; n /= 100000; } return s; }
    int s = 0; while (n) { s += (int)(n % (u64)g_base); n /= (u64)g_base; } return s;
}

/* ---- per-thread tracker + results ------------------------------------ */
typedef struct {
    u64 lo, hi, A, B; int tid;
    u64 *first; int *dsv; u64 *members; size_t count, cap;
    int len; u64 res; int ds; u64 *run; int done;
} Targ;

static void record(Targ *t) {
    if (t->count == t->cap) {
        t->cap = t->cap ? t->cap * 2 : 16;
        t->first   = realloc(t->first,   t->cap * sizeof(u64));
        t->dsv     = realloc(t->dsv,     t->cap * sizeof(int));
        t->members = realloc(t->members, t->cap * (size_t)g_L * sizeof(u64));
    }
    t->first[t->count] = t->run[0];
    t->dsv[t->count]   = t->ds;
    memcpy(&t->members[t->count * (size_t)g_L], t->run, (size_t)g_L * sizeof(u64));
    t->count++;
}

static inline void feed(Targ *t, u64 q) {
    if (t->done) return;
    /* For base 10 the divisor is the literal 9, which the compiler turns into
     * a multiply-shift; the runtime g_M path (other bases) keeps a real mod. */
    u64 rq = (g_base == 10) ? (q % 9u) : (q % g_M);
    int extended = 0;
    if (t->len >= 1 && rq == t->res) {                 /* cheap reject first */
        if (t->ds < 0) t->ds = dsum(t->run[0]);
        if (dsum(q) == t->ds) {
            if (t->len < g_L) t->run[t->len] = q;
            t->len++; extended = 1;
        }
    }
    if (!extended) { t->len = 1; t->res = rq; t->ds = -1; t->run[0] = q; }

    if (t->len == g_L) {
        u64 f = t->run[0];
        if (f >= t->lo && f < t->hi) {
            if (t->ds < 0) t->ds = dsum(t->run[0]);
            record(t);
            if (!g_find_all) t->done = 1;
        }
    }
}

/* ---- progress -------------------------------------------------------- */
static pthread_mutex_t g_pmx = PTHREAD_MUTEX_INITIALIZER;
static u64 g_primes_done;
static double g_t_start, g_t_last, g_rate;   /* timing + windowed rate */
static u64 g_count_last;

/* ---- worker: stream primes from primesieve over [A,B) ---------------- */
static void *worker(void *arg) {
    Targ *t = (Targ *)arg;
    t->len = 0; t->res = 0; t->ds = -1; t->done = 0;
    t->run = (u64 *)malloc((size_t)g_L * sizeof(u64));

    primesieve_iterator it;
    primesieve_init(&it);
    primesieve_jump_to(&it, t->A, t->B);     /* first prime >= A */

    u64 local = 0, p;
    while ((p = primesieve_next_prime(&it)) < t->B) {
        feed(t, p);
        if (t->done) break;
        if (!g_quiet && ((++local & ((1u << 22) - 1)) == 0)) {
            pthread_mutex_lock(&g_pmx);
            g_primes_done += (1u << 22);
            double now = mono(), since = now - g_t_last;
            if (since >= 0.2) {                       /* refresh windowed rate */
                g_rate = (double)(g_primes_done - g_count_last) / since;
                g_t_last = now; g_count_last = g_primes_done;
            }
            double elapsed = now - g_t_start;
            double avg = elapsed > 0 ? (double)g_primes_done / elapsed : 0;
            double cur = g_rate > 0 ? g_rate : avg;
            fprintf(stderr,
                "\rprocessed ~%llu primes | %.1f M/s (avg %.1f M/s)        ",
                (ull)g_primes_done, cur / 1e6, avg / 1e6);
            fflush(stderr);
            pthread_mutex_unlock(&g_pmx);
        }
    }
    primesieve_free_iterator(&it);
    return NULL;
}

static u64 align30(u64 x) { return x - (x % 30); }

typedef struct { u64 first; int ds; const u64 *mem; } Rec;
static int reccmp(const void *a, const void *b) {
    u64 x = ((const Rec *)a)->first, y = ((const Rec *)b)->first;
    return (x > y) - (x < y);
}

int main(int argc, char **argv) {
    g_L = 0; g_base = 10; g_find_all = 0; g_quiet = 0;
    g_start = 0; g_end = 1000000000ULL;
    int nthreads = 1;
    for (int i = 1; i < argc; i++) {
        if      ((!strcmp(argv[i],"--length") ||!strcmp(argv[i],"-n")) && i+1<argc) g_L=atoi(argv[++i]);
        else if ((!strcmp(argv[i],"--start")  ||!strcmp(argv[i],"-S")) && i+1<argc) g_start=parse_num(argv[++i]);
        else if ((!strcmp(argv[i],"--end")    ||!strcmp(argv[i],"-E")) && i+1<argc) g_end=parse_num(argv[++i]);
        else if ((!strcmp(argv[i],"--limit")  ||!strcmp(argv[i],"-L")) && i+1<argc) g_end=parse_num(argv[++i]); /* alias for -E */
        else if ((!strcmp(argv[i],"--base")   ||!strcmp(argv[i],"-b")) && i+1<argc) g_base=atoi(argv[++i]);
        else if ((!strcmp(argv[i],"--threads")||!strcmp(argv[i],"-t")) && i+1<argc) nthreads=atoi(argv[++i]);
        else if ( !strcmp(argv[i],"--all"))   g_find_all=1;
        else if ( !strcmp(argv[i],"--quiet")) g_quiet=1;
        else { fprintf(stderr,"bad arg: %s\n", argv[i]); return 2; }
    }
    if (g_L < 1) { fprintf(stderr,"need --length >= 1\n"); return 2; }
    if (g_base < 2 || g_base > 36) { fprintf(stderr,"base 2..36\n"); return 2; }
    if (g_end <= g_start) { fprintf(stderr,"need end (-E) > start (-S)\n"); return 2; }
    if (nthreads < 1) nthreads = 1;
    g_M = (u64)g_base - 1; if (g_M == 0) g_M = 1;
    if (g_base == 10) for (int i = 1; i < 100000; i++) ds5[i] = ds5[i/10] + i%10;

    u64 span = g_end - g_start;
    while (nthreads > 1 && span / (u64)nthreads < (u64)OVERLAP * 8) nthreads--;

    Targ *ta = (Targ *)calloc(nthreads, sizeof(Targ));
    u64 chunk = align30(span / (u64)nthreads); if (chunk < 30) chunk = 30;
    for (int i = 0; i < nthreads; i++) {
        u64 lo = g_start + (u64)i * chunk;
        u64 hi = (i == nthreads - 1) ? g_end : g_start + (u64)(i + 1) * chunk;
        ta[i].lo = lo; ta[i].hi = hi; ta[i].tid = i;
        ta[i].A = (lo > OVERLAP) ? align30(lo - OVERLAP) : 0;
        ta[i].B = (hi + OVERLAP < g_end) ? hi + OVERLAP : g_end;
    }

    g_t_start = g_t_last = mono();
    g_count_last = 0; g_rate = 0;
    pthread_t *th = (pthread_t *)malloc(nthreads * sizeof(pthread_t));
    for (int i = 0; i < nthreads; i++) pthread_create(&th[i], NULL, worker, &ta[i]);
    for (int i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
    if (!g_quiet) fprintf(stderr, "\r%60s\r", "");

    size_t total = 0;
    for (int i = 0; i < nthreads; i++) total += ta[i].count;
    Rec *recs = (Rec *)malloc((total ? total : 1) * sizeof(Rec));
    size_t k = 0;
    for (int i = 0; i < nthreads; i++)
        for (size_t j = 0; j < ta[i].count; j++)
            recs[k++] = (Rec){ ta[i].first[j], ta[i].dsv[j], &ta[i].members[j*(size_t)g_L] };
    qsort(recs, total, sizeof(Rec), reccmp);

    if (total == 0) {
        printf("No run of %d consecutive primes with equal digit sum "
               "(base %d) in [%llu, %llu).\n",
               g_L, g_base, (ull)g_start, (ull)g_end);
    } else {
        size_t show = g_find_all ? total : 1;
        printf("%s of %d consecutive primes sharing a digit sum (base %d), "
               "in [%llu, %llu):\n",
               g_find_all ? "All runs" : "First run", g_L, g_base,
               (ull)g_start, (ull)g_end);
        for (size_t r = 0; r < show; r++) {
            printf("  digit sum %d: [", recs[r].ds);
            for (int j = 0; j < g_L; j++)
                printf("%llu%s", (ull)recs[r].mem[j], j+1 < g_L ? ", " : "");
            printf("]\n");
        }
    }
    return 0;
}
