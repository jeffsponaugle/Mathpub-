/*
 * cchain.c — exhaustive search for the smallest Cunningham chain of length k.
 *
 *   kind 1:  p -> 2p+1   (members m_i = 2^i*(p+1) - 1, i = 0..k-1)
 *   kind 2:  p -> 2p-1   (members m_i = 2^i*(p-1) + 1, i = 0..k-1)
 *
 * Target: OEIS A057330 a(17) (kind 2, "chain of length at least k" semantics).
 * Also handles kind 1 (cf. A005602 and derived "at least" values).
 *
 * Method:
 *   1. Wheel: p must avoid, for every small prime q, the residues that force
 *      some chain member to be divisible by q.  Forbidden residues mod q are
 *      r = eps*(2^-i - 1) mod q, i = 0..min(k, ord_2(q))-1  (eps=+1 kind 1,
 *      eps=-1 kind 2).  CRT over wheel primes {2,3,5,7,11,13,17,19,23} gives
 *      the admissible residue classes mod M = 223092870.
 *   2. Segmented sieve: for each class r, candidates are p = base + r + M*t.
 *      For each sieve prime q in (23, B], mark t where some member ≡ 0 mod q.
 *   3. Survivors get chain-prefix probable-prime tests (GMP).  Full-length
 *      hits are re-verified with many MR rounds and logged.
 *
 * Exhaustiveness: sieving is only valid when p > B (a chain member equal to
 * a sieve prime q would be prime, not composite; that requires p <= B).
 * The tool warns when start <= B.
 *
 * Build: make    (requires GMP; -lgmp -lpthread)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <pthread.h>
#include <stdatomic.h>
#include <signal.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>
#include <getopt.h>
#include <gmp.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef unsigned __int128 u128;

/* ------------------------------------------------------------------ */
/* small utilities                                                     */
/* ------------------------------------------------------------------ */

static char *u128_str(u128 v, char *buf /* >= 40 bytes */)
{
    char tmp[48];
    int i = 0;
    if (v == 0) tmp[i++] = '0';
    while (v > 0) { tmp[i++] = (char)('0' + (int)(v % 10)); v /= 10; }
    int j = 0;
    while (i > 0) buf[j++] = tmp[--i];
    buf[j] = 0;
    return buf;
}

static double u128_d(u128 v) { return (double)v; }

/* parse decimal or scientific ("3.2e18") into u128; returns 0 on error */
static int parse_u128(const char *s, u128 *out)
{
    u128 mant = 0;
    int frac = 0, seen = 0, exp10 = 0;
    const char *p = s;
    if (!*p) return 0;
    for (; *p && *p != 'e' && *p != 'E'; p++) {
        if (*p == '.') { if (frac) return 0; frac = 1; continue; }
        if (*p < '0' || *p > '9') return 0;
        if (mant > (((u128)-1) - 9) / 10) return 0;
        mant = mant * 10 + (u128)(*p - '0');
        if (frac) exp10--;
        seen = 1;
    }
    if (!seen) return 0;
    if (*p == 'e' || *p == 'E') {
        p++;
        if (*p == '+') p++;
        int e = 0;
        if (!*p) return 0;
        for (; *p; p++) {
            if (*p < '0' || *p > '9') return 0;
            e = e * 10 + (*p - '0');
            if (e > 60) return 0;
        }
        exp10 += e;
    }
    if (exp10 < 0) return 0;               /* must be an exact integer */
    while (exp10-- > 0) {
        if (mant > ((u128)-1) / 10) return 0;
        mant *= 10;
    }
    *out = mant;
    return 1;
}

static u64 modinv_u64(u64 a, u64 m)       /* gcd(a,m)=1 */
{
    int64_t t = 0, nt = 1;
    int64_t r = (int64_t)m, nr = (int64_t)(a % m);
    while (nr != 0) {
        int64_t qq = r / nr, tmp;
        tmp = t - qq * nt; t = nt; nt = tmp;
        tmp = r - qq * nr; r = nr; nr = tmp;
    }
    if (t < 0) t += (int64_t)m;
    return (u64)t;
}

static void mpz_set_u128(mpz_t z, u128 v)
{
    mpz_set_ui(z, (unsigned long)(u64)(v >> 64));
    mpz_mul_2exp(z, z, 64);
    mpz_add_ui(z, z, (unsigned long)(u64)v);
}

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ------------------------------------------------------------------ */
/* configuration                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    int   kind;            /* 1 or 2 */
    int   k;               /* chain length target */
    u128  start, end;
    int   threads;
    u32   B;               /* sieve limit */
    int   chunk_log2;      /* t-window bits per work unit */
    int   find_all;        /* keep going after first confirmed minimum */
    int   quiet;
    double status_iv, ckpt_iv;
    const char *ckpt_file;
    const char *results_file;
    int   resume;
} cfg_t;

/* forbidden residues mod odd prime q; returns count (<= min(k, ord_2(q))) */
static int forbidden_residues(u32 q, int k, int kind, u32 *out)
{
    u64 inv2 = ((u64)q + 1) / 2, x = 1;
    int n = 0;
    for (int i = 0; i < k; i++) {
        u64 f = (kind == 1) ? (x + q - 1) % q : (1 + (u64)q - x) % q;
        out[n++] = (u32)f;
        x = x * inv2 % q;
        if (x == 1) break;                 /* full 2-cycle covered */
    }
    return n;
}

/* ------------------------------------------------------------------ */
/* search state (one search at a time)                                 */
/* ------------------------------------------------------------------ */

#define MAXK 64
#define RING 4096                          /* superblock completion ring */
#define MAXFOUND 256

typedef struct { u32 q, minv, off, nres; } sprime_t;

static struct {
    cfg_t cfg;
    /* wheel */
    u64   M;
    u64  *cls;      u32 C;
    /* sieve primes */
    sprime_t *sp;   u32 nsp;
    u32  *fr;                              /* flattened forbidden residues */
    /* geometry */
    u128  base_p;                          /* p = base_p + cls[c] + M*t */
    u64   t_last;                          /* inclusive max t */
    u64   end_sb;                          /* inclusive max superblock */
    u64   tchunk;
    /* scheduling */
    _Atomic u64 next_unit;
    u64   total_units;
    _Atomic u32 done_ring[RING];
    pthread_mutex_t lock;                  /* frontier, found, best */
    u64   frontier_sb;                     /* all sb < frontier_sb complete */
    _Atomic u64 frontier_mirror;
    /* results */
    u128  found[MAXFOUND]; int nfound;
    u128  found_min; int have_found;
    int   best_len; u128 best_p;
    /* stats */
    _Atomic u64 units_done, cand, prp;
    double t_start, elapsed_prev;          /* elapsed_prev from checkpoint */
    _Atomic int stop;
    int   exhausted;
} S;

static volatile sig_atomic_t g_sigint = 0, g_siginfo = 0;
static void on_sigint(int sig) { (void)sig; g_sigint = 1; }
#ifdef SIGINFO
static void on_siginfo(int sig) { (void)sig; g_siginfo = 1; }
#endif

/* ------------------------------------------------------------------ */
/* table construction                                                  */
/* ------------------------------------------------------------------ */

static const u32 WHEEL[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23 };
#define NWHEEL 9

static int cmp_u64(const void *a, const void *b)
{
    u64 x = *(const u64 *)a, y = *(const u64 *)b;
    return x < y ? -1 : x > y;
}

static void build_wheel(void)
{
    int k = S.cfg.k, kind = S.cfg.kind;
    u64 M = 1;
    u64 *res = malloc(sizeof(u64));
    u64 nres = 1;
    res[0] = 0;
    u32 fbuf[MAXK], allow[64];

    for (int w = 0; w < NWHEEL; w++) {
        u32 q = WHEEL[w];
        int na = 0;
        if (q == 2) {
            allow[na++] = 1;               /* p must be odd */
        } else {
            int nf = forbidden_residues(q, k, kind, fbuf);
            for (u32 a = 0; a < q; a++) {
                int bad = 0;
                for (int i = 0; i < nf; i++) if (fbuf[i] == a) { bad = 1; break; }
                if (!bad) allow[na++] = a;
            }
        }
        if (na == 0) {
            fprintf(stderr, "error: no admissible residues mod %u for k=%d kind %d\n",
                    q, k, kind);
            exit(1);
        }
        u64 *nres_arr = malloc(sizeof(u64) * nres * na);
        u64 minv = modinv_u64(M % q, q);
        u64 cnt = 0;
        for (u64 i = 0; i < nres; i++) {
            for (int j = 0; j < na; j++) {
                u64 d = ((u64)allow[j] + q - res[i] % q) % q;
                nres_arr[cnt++] = res[i] + M * (d * minv % q);
            }
        }
        free(res);
        res = nres_arr;
        nres = cnt;
        M *= q;
    }
    qsort(res, nres, sizeof(u64), cmp_u64);
    S.M = M;
    S.cls = res;
    S.C = (u32)nres;
}

static void build_sieve_primes(void)
{
    u32 B = S.cfg.B;
    /* simple odd sieve of Eratosthenes up to B */
    u64 nbits = (u64)B / 2 + 1;
    unsigned char *comp = calloc((nbits + 7) / 8, 1);
    for (u64 i = 1; ; i++) {               /* odd number 2i+1 */
        u64 p = 2 * i + 1;
        if (p * p > B) break;
        if (comp[i >> 3] & (1 << (i & 7))) continue;
        for (u64 j = (p * p - 1) / 2; j <= (u64)B / 2; j += p)
            comp[j >> 3] |= (unsigned char)(1 << (j & 7));
    }
    /* count primes > 23 */
    u32 cap = 1024, n = 0;
    sprime_t *sp = malloc(cap * sizeof(sprime_t));
    u32 frcap = 4096, frn = 0;
    u32 *fr = malloc(frcap * sizeof(u32));
    u32 fbuf[MAXK];
    int k = S.cfg.k, kind = S.cfg.kind;
    for (u64 i = 1; 2 * i + 1 <= B; i++) {
        if (comp[i >> 3] & (1 << (i & 7))) continue;
        u32 q = (u32)(2 * i + 1);
        if (q <= 23) continue;
        int nf = forbidden_residues(q, k, kind, fbuf);
        if (n == cap) { cap *= 2; sp = realloc(sp, cap * sizeof(sprime_t)); }
        while (frn + (u32)nf > frcap) { frcap *= 2; fr = realloc(fr, frcap * sizeof(u32)); }
        sp[n].q = q;
        sp[n].minv = (u32)modinv_u64(S.M % q, q);
        sp[n].off = frn;
        sp[n].nres = (u32)nf;
        for (int j = 0; j < nf; j++) fr[frn++] = fbuf[j];
        n++;
    }
    free(comp);
    S.sp = sp; S.nsp = n; S.fr = fr;
}

/* ------------------------------------------------------------------ */
/* chain testing                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    mpz_t m;
    int id;
    unsigned char *bm;
} worker_t;

/* length of prime prefix of the chain from p (capped at cap) */
static int chain_prefix_len(mpz_t m, u128 p, int cap, int reps, int kind,
                            _Atomic u64 *prp_ctr)
{
    mpz_set_u128(m, p);
    int len = 0;
    for (int i = 0; i < cap; i++) {
        if (prp_ctr) atomic_fetch_add_explicit(prp_ctr, 1, memory_order_relaxed);
        if (!mpz_probab_prime_p(m, reps)) break;
        len++;
        mpz_mul_2exp(m, m, 1);
        if (kind == 1) mpz_add_ui(m, m, 1);
        else           mpz_sub_ui(m, m, 1);
    }
    return len;
}

static void print_chain(FILE *f, u128 p, int len, int kind)
{
    char buf[48];
    mpz_t m;
    mpz_init(m);
    mpz_set_u128(m, p);
    for (int i = 0; i < len; i++) {
        gmp_fprintf(f, "%s%Zd", i ? ", " : "", m);
        mpz_mul_2exp(m, m, 1);
        if (kind == 1) mpz_add_ui(m, m, 1);
        else           mpz_sub_ui(m, m, 1);
    }
    mpz_clear(m);
    (void)buf;
}

static void record_found(worker_t *w, u128 p)
{
    char buf[48];
    pthread_mutex_lock(&S.lock);
    int dup = 0;
    for (int i = 0; i < S.nfound; i++) if (S.found[i] == p) dup = 1;
    if (!dup && S.nfound < MAXFOUND) {
        S.found[S.nfound++] = p;
        if (!S.have_found || p < S.found_min) { S.found_min = p; S.have_found = 1; }
        FILE *f = fopen(S.cfg.results_file, "a");
        if (f) {
            time_t t = time(NULL);
            char ts[64];
            strftime(ts, sizeof ts, "%Y-%m-%d %H:%M:%S", localtime(&t));
            fprintf(f, "[%s] kind=%d k=%d p=%s chain: ", ts, S.cfg.kind, S.cfg.k,
                    u128_str(p, buf));
            print_chain(f, p, S.cfg.k, S.cfg.kind);
            fprintf(f, "\n");
            fclose(f);
        }
        printf("\n*** FOUND length-%d chain (kind %d): p = %s ***\n",
               S.cfg.k, S.cfg.kind, u128_str(p, buf));
        printf("    chain: ");
        print_chain(stdout, p, S.cfg.k, S.cfg.kind);
        printf("\n    (minimality confirmed once the frontier passes this p)\n");
        fflush(stdout);
    }
    pthread_mutex_unlock(&S.lock);
    (void)w;
}

/* ------------------------------------------------------------------ */
/* worker                                                              */
/* ------------------------------------------------------------------ */

static void *worker(void *arg)
{
    worker_t *w = arg;
    const u64 M = S.M;
    const u64 tchunk = S.tchunk;
    const u32 C = S.C;
    const int k = S.cfg.k, kind = S.cfg.kind;
    const u128 pmin = S.cfg.start, pmax = S.cfg.end;

    for (;;) {
        if (atomic_load_explicit(&S.stop, memory_order_relaxed)) break;
        u64 id = atomic_fetch_add(&S.next_unit, 1);
        if (id >= S.total_units) break;
        u64 sb = id / C;
        u32 ci = (u32)(id % C);
        u64 r = S.cls[ci];
        u64 t0 = sb * tchunk;
        u64 tlen = tchunk;
        if (t0 + tlen - 1 > S.t_last) tlen = S.t_last - t0 + 1;

        u64 words = (tlen + 63) / 64;
        memset(w->bm, 0, words * 8);
        unsigned char *bm = w->bm;

        /* sieve */
        const sprime_t *sp = S.sp;
        const u32 *fr = S.fr;
        for (u32 i = 0; i < S.nsp; i++) {
            u32 q = sp[i].q;
            u64 minv = sp[i].minv;
            u64 rmod = r % q;
            u64 t0mod = t0 % q;
            const u32 *f = fr + sp[i].off;
            for (u32 j = 0; j < sp[i].nres; j++) {
                u64 s = (((u64)f[j] + q - rmod) % q) * minv % q;
                u64 idx = (s + q - t0mod) % q;
                for (; idx < tlen; idx += q)
                    bm[idx >> 3] |= (unsigned char)(1u << (idx & 7));
            }
        }

        /* scan survivors */
        u64 *bw = (u64 *)bm;
        for (u64 wi = 0; wi < words; wi++) {
            u64 x = ~bw[wi];
            while (x) {
                u64 bit = (u64)__builtin_ctzll(x);
                x &= x - 1;
                u64 t = wi * 64 + bit;
                if (t >= tlen) break;
                u128 p = S.base_p + r + (u128)M * (t0 + t);
                if (p < pmin || p > pmax) continue;
                atomic_fetch_add_explicit(&S.cand, 1, memory_order_relaxed);
                int len = chain_prefix_len(w->m, p, k, 2, kind, &S.prp);
                if (len >= k) {
                    /* strong re-check before recording */
                    if (chain_prefix_len(w->m, p, k, 40, kind, NULL) >= k)
                        record_found(w, p);
                } else if (len > S.best_len) {
                    pthread_mutex_lock(&S.lock);
                    if (len > S.best_len) { S.best_len = len; S.best_p = p; }
                    pthread_mutex_unlock(&S.lock);
                }
            }
        }

        /* completion bookkeeping */
        atomic_fetch_add_explicit(&S.units_done, 1, memory_order_relaxed);
        u32 dc = atomic_fetch_add(&S.done_ring[sb % RING], 1) + 1;
        if (dc == C) {
            pthread_mutex_lock(&S.lock);
            while (S.frontier_sb <= S.end_sb &&
                   atomic_load(&S.done_ring[S.frontier_sb % RING]) == C) {
                atomic_store(&S.done_ring[S.frontier_sb % RING], 0);
                S.frontier_sb++;
            }
            atomic_store(&S.frontier_mirror, S.frontier_sb);
            pthread_mutex_unlock(&S.lock);
        }
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/* checkpointing                                                       */
/* ------------------------------------------------------------------ */

static u128 frontier_p(void)
{
    u64 fsb = atomic_load(&S.frontier_mirror);
    u64 t = fsb * S.tchunk;
    if (t > S.t_last + 1) t = S.t_last + 1;
    u128 p = S.base_p + (u128)S.M * t;
    if (p > S.cfg.end) p = S.cfg.end;
    return p;
}

static void write_checkpoint(void)
{
    if (!S.cfg.ckpt_file) return;
    char tmp[1024], buf[48];
    snprintf(tmp, sizeof tmp, "%s.tmp", S.cfg.ckpt_file);
    FILE *f = fopen(tmp, "w");
    if (!f) { fprintf(stderr, "warn: cannot write checkpoint %s\n", tmp); return; }
    pthread_mutex_lock(&S.lock);
    fprintf(f, "cchain-v1\n");
    fprintf(f, "kind=%d\nk=%d\n", S.cfg.kind, S.cfg.k);
    fprintf(f, "start=%s\n", u128_str(S.cfg.start, buf));
    fprintf(f, "end=%s\n", u128_str(S.cfg.end, buf));
    fprintf(f, "B=%u\nchunk_log2=%d\nM=%llu\nC=%u\n",
            S.cfg.B, S.cfg.chunk_log2, (unsigned long long)S.M, S.C);
    fprintf(f, "frontier_sb=%llu\n", (unsigned long long)S.frontier_sb);
    fprintf(f, "elapsed=%.1f\n", S.elapsed_prev + (now_s() - S.t_start));
    fprintf(f, "cand=%llu\nprp=%llu\n",
            (unsigned long long)atomic_load(&S.cand),
            (unsigned long long)atomic_load(&S.prp));
    fprintf(f, "best_len=%d\nbest_p=%s\n", S.best_len, u128_str(S.best_p, buf));
    fprintf(f, "nfound=%d\n", S.nfound);
    for (int i = 0; i < S.nfound; i++)
        fprintf(f, "found=%s\n", u128_str(S.found[i], buf));
    pthread_mutex_unlock(&S.lock);
    fclose(f);
    rename(tmp, S.cfg.ckpt_file);
}

static int read_checkpoint(void)
{
    FILE *f = fopen(S.cfg.ckpt_file, "r");
    if (!f) { fprintf(stderr, "error: cannot open checkpoint %s\n", S.cfg.ckpt_file); return 0; }
    char line[256], buf[48];
    int ok = 1, kind = 0, k = 0, chunk = 0, nf = 0;
    u32 B = 0, C = 0;
    u64 M = 0, fsb = 0;
    u128 start = 0, end = 0;
    double elapsed = 0;
    if (!fgets(line, sizeof line, f) || strncmp(line, "cchain-v1", 9) != 0) ok = 0;
    while (ok && fgets(line, sizeof line, f)) {
        char *nl = strchr(line, '\n'); if (nl) *nl = 0;
        if (!strncmp(line, "kind=", 5)) kind = atoi(line + 5);
        else if (!strncmp(line, "k=", 2)) k = atoi(line + 2);
        else if (!strncmp(line, "start=", 6)) ok = parse_u128(line + 6, &start);
        else if (!strncmp(line, "end=", 4)) ok = parse_u128(line + 4, &end);
        else if (!strncmp(line, "B=", 2)) B = (u32)strtoul(line + 2, NULL, 10);
        else if (!strncmp(line, "chunk_log2=", 11)) chunk = atoi(line + 11);
        else if (!strncmp(line, "M=", 2)) M = strtoull(line + 2, NULL, 10);
        else if (!strncmp(line, "C=", 2)) C = (u32)strtoul(line + 2, NULL, 10);
        else if (!strncmp(line, "frontier_sb=", 12)) fsb = strtoull(line + 12, NULL, 10);
        else if (!strncmp(line, "elapsed=", 8)) elapsed = atof(line + 8);
        else if (!strncmp(line, "cand=", 5)) atomic_store(&S.cand, strtoull(line + 5, NULL, 10));
        else if (!strncmp(line, "prp=", 4)) atomic_store(&S.prp, strtoull(line + 4, NULL, 10));
        else if (!strncmp(line, "best_len=", 9)) S.best_len = atoi(line + 9);
        else if (!strncmp(line, "best_p=", 7)) ok = parse_u128(line + 7, &S.best_p);
        else if (!strncmp(line, "nfound=", 7)) nf = atoi(line + 7);
        else if (!strncmp(line, "found=", 6)) {
            u128 p; if (parse_u128(line + 6, &p) && S.nfound < MAXFOUND) {
                S.found[S.nfound++] = p;
                if (!S.have_found || p < S.found_min) { S.found_min = p; S.have_found = 1; }
            }
        }
    }
    fclose(f);
    if (!ok || kind != S.cfg.kind || k != S.cfg.k || start != S.cfg.start ||
        end != S.cfg.end || B != S.cfg.B || chunk != S.cfg.chunk_log2 ||
        M != S.M || C != S.C || nf != S.nfound) {
        fprintf(stderr, "error: checkpoint mismatch with current configuration\n");
        fprintf(stderr, "  (kind/k/start/end/B/chunk must match; got kind=%d k=%d "
                        "start=%s B=%u chunk=%d)\n",
                kind, k, u128_str(start, buf), B, chunk);
        return 0;
    }
    S.frontier_sb = fsb;
    atomic_store(&S.frontier_mirror, fsb);
    S.elapsed_prev = elapsed;
    return 1;
}

/* ------------------------------------------------------------------ */
/* search driver                                                       */
/* ------------------------------------------------------------------ */

typedef enum { R_FOUND, R_EXHAUSTED, R_INTERRUPTED } result_t;

static result_t run_search(u128 *out_min)
{
    char b1[48], b2[48];
    cfg_t *cfg = &S.cfg;

    build_wheel();
    build_sieve_primes();

    S.base_p = (cfg->start / S.M) * S.M;
    S.t_last = (u64)((cfg->end - S.base_p) / S.M);
    S.tchunk = (u64)1 << cfg->chunk_log2;
    S.end_sb = S.t_last >> cfg->chunk_log2;
    S.total_units = (S.end_sb + 1) * S.C;
    pthread_mutex_init(&S.lock, NULL);
    for (int i = 0; i < RING; i++) atomic_store(&S.done_ring[i], 0);
    S.frontier_sb = 0;
    atomic_store(&S.frontier_mirror, 0);
    atomic_store(&S.stop, 0);
    S.exhausted = 0;

    if (cfg->resume) {
        if (!read_checkpoint()) exit(1);
    }
    atomic_store(&S.next_unit, S.frontier_sb * S.C);

    if (!cfg->quiet) {
        printf("cchain: kind %d (p -> 2p%+d), k=%d\n", cfg->kind,
               cfg->kind == 1 ? 1 : -1, cfg->k);
        printf("  range [%s, %s]\n", u128_str(cfg->start, b1), u128_str(cfg->end, b2));
        printf("  wheel M=%llu, %u admissible classes; sieve B=%u (%u primes)\n",
               (unsigned long long)S.M, S.C, cfg->B, S.nsp);
        printf("  chunk 2^%d, %llu superblocks x %u units, %d threads\n",
               cfg->chunk_log2, (unsigned long long)(S.end_sb + 1), S.C, cfg->threads);
        if (cfg->resume)
            printf("  resumed at superblock %llu (frontier p=%s)\n",
                   (unsigned long long)S.frontier_sb, u128_str(frontier_p(), b1));
        if (cfg->start <= cfg->B)
            printf("  WARNING: start <= B; exhaustiveness only guaranteed for p > %u\n",
                   cfg->B);
        fflush(stdout);
    }

    S.t_start = now_s();
    u64 sb0 = S.frontier_sb;

    pthread_t *th = malloc(sizeof(pthread_t) * (size_t)cfg->threads);
    worker_t *ws = calloc((size_t)cfg->threads, sizeof(worker_t));
    for (int i = 0; i < cfg->threads; i++) {
        ws[i].id = i;
        mpz_init(ws[i].m);
        ws[i].bm = malloc(((size_t)S.tchunk + 63) / 64 * 8);
        pthread_create(&th[i], NULL, worker, &ws[i]);
    }

    double last_status = now_s(), last_ckpt = now_s();
    result_t res = R_INTERRUPTED;
    for (;;) {
        usleep(200000);
        double t = now_s();
        if (g_sigint) { atomic_store(&S.stop, 1); res = R_INTERRUPTED; break; }

        pthread_mutex_lock(&S.lock);
        int have = S.have_found; u128 fmin = S.found_min;
        pthread_mutex_unlock(&S.lock);
        u128 fp = frontier_p();
        if (have && !cfg->find_all && fp > fmin) { res = R_FOUND; break; }
        if (atomic_load(&S.frontier_mirror) > S.end_sb) {
            res = have ? R_FOUND : R_EXHAUSTED;
            S.exhausted = 1;
            break;
        }

        if (!cfg->quiet && (t - last_status >= cfg->status_iv || g_siginfo)) {
            g_siginfo = 0;
            last_status = t;
            double el = S.elapsed_prev + (t - S.t_start);
            u64 fsb = atomic_load(&S.frontier_mirror);
            double span_done = (double)(fsb - sb0) * (double)S.tchunk * (double)S.M;
            double rate = span_done / (t - S.t_start + 1e-9);
            double togo = u128_d(cfg->end) - u128_d(fp);
            double eta = rate > 0 ? togo / rate : 0;
            pthread_mutex_lock(&S.lock);
            int bl = S.best_len; u128 bp = S.best_p;
            pthread_mutex_unlock(&S.lock);
            printf("[%8.0fs] p=%.4e  sb %llu/%llu  %.3e p/s  cand %llu  prp %llu  "
                   "best %d/%d @ %s  ETA %.1fd\n",
                   el, u128_d(fp), (unsigned long long)fsb,
                   (unsigned long long)(S.end_sb + 1), rate,
                   (unsigned long long)atomic_load(&S.cand),
                   (unsigned long long)atomic_load(&S.prp),
                   bl, cfg->k, u128_str(bp, b1), eta / 86400.0);
            fflush(stdout);
        }
        if (cfg->ckpt_file && t - last_ckpt >= cfg->ckpt_iv) {
            last_ckpt = t;
            write_checkpoint();
        }
    }
    atomic_store(&S.stop, 1);
    for (int i = 0; i < cfg->threads; i++) pthread_join(th[i], NULL);
    if (cfg->ckpt_file) write_checkpoint();

    for (int i = 0; i < cfg->threads; i++) { mpz_clear(ws[i].m); free(ws[i].bm); }
    free(ws); free(th);

    if (!cfg->quiet) {
        u128 fp = frontier_p();
        printf("\n");
        if (res == R_FOUND) {
            pthread_mutex_lock(&S.lock);
            u128 fmin = S.found_min;
            pthread_mutex_unlock(&S.lock);
            printf("RESULT: smallest length-%d chain (kind %d) in range starts at p = %s\n",
                   cfg->k, cfg->kind, u128_str(fmin, b1));
            printf("        exhaustively verified: no smaller start in [%s, %s]\n",
                   u128_str(cfg->start, b1), u128_str(fmin, b2));
            *out_min = fmin;
        } else if (res == R_EXHAUSTED) {
            printf("RESULT: no length-%d chain (kind %d) in [%s, %s]\n",
                   cfg->k, cfg->kind, u128_str(cfg->start, b1), u128_str(cfg->end, b2));
        } else {
            printf("INTERRUPTED: frontier at p = %s (checkpoint written)\n",
                   u128_str(fp, b1));
        }
        fflush(stdout);
    } else if (res == R_FOUND) {
        pthread_mutex_lock(&S.lock);
        *out_min = S.found_min;
        pthread_mutex_unlock(&S.lock);
    }

    free(S.cls); free(S.sp); free(S.fr);
    S.cls = NULL; S.sp = NULL; S.fr = NULL;
    return res;
}

/* ------------------------------------------------------------------ */
/* verify / selftest                                                   */
/* ------------------------------------------------------------------ */

static int verify_p(u128 p, int kind, int quiet)
{
    mpz_t m;
    mpz_init(m);
    mpz_set_u128(m, p);
    int len = 0;
    if (!quiet) printf("chain from p (kind %d, p -> 2p%+d):\n", kind, kind == 1 ? 1 : -1);
    for (int i = 0; i < MAXK; i++) {
        int pr = mpz_probab_prime_p(m, 50);
        if (!quiet) gmp_printf("  m_%-2d = %Zd  %s\n", i, m, pr ? "prime" : "COMPOSITE");
        if (!pr) break;
        len++;
        mpz_mul_2exp(m, m, 1);
        if (kind == 1) mpz_add_ui(m, m, 1);
        else           mpz_sub_ui(m, m, 1);
    }
    mpz_clear(m);
    if (!quiet) printf("chain length: %d\n", len);
    return len;
}

/* A057330 (kind 2, "length >= n"): */
static const char *KNOWN2[] = { "2","2","2","1531","1531","16651","16651",
    "15514861","857095381","205528443121","1389122693971","216857744866621",
    "758083947856951","69257563144280941","69257563144280941",
    "3203000719597029781" };
/* kind 1, "length >= n", derived from A005602 (complete chains): */
static const char *KNOWN1[] = { "2","2","2","2","2","89","1122659","19099919",
    "85864769","26089808579","554688278429","554688278429","4090932431513069",
    "95405042230542329","90616211958465842219","810433818265726529159" };

static void verify_known(void)
{
    char buf[48];
    for (int kind = 2; kind >= 1; kind--) {
        const char **tab = kind == 2 ? KNOWN2 : KNOWN1;
        printf("known smallest chain starts, kind %d (%s):\n", kind,
               kind == 2 ? "A057330" : "derived from A005602");
        for (int n = 1; n <= 16; n++) {
            u128 p;
            parse_u128(tab[n - 1], &p);
            int len = verify_p(p, kind, 1);
            printf("  a(%2d) = %-22s  chain length %2d  %s\n", n,
                   u128_str(p, buf), len, len >= n ? "OK" : "** FAIL **");
        }
    }
}

static int selftest(int threads)
{
    struct { int kind, k; const char *exp; } cases[] = {
        { 2,  8, "15514861"      },
        { 2, 10, "205528443121"  },
        { 2, 12, "216857744866621" },
        { 1,  8, "19099919"      },
        { 1, 10, "26089808579"   },
    };
    char b1[48], b2[48];
    int pass = 1;
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        u128 expct, endr, got = 0;
        parse_u128(cases[i].exp, &expct);
        endr = expct + 1000000;
        memset(&S, 0, sizeof S);
        S.cfg = (cfg_t){ .kind = cases[i].kind, .k = cases[i].k, .start = 3,
                         .end = endr, .threads = threads, .B = 1u << 16,
                         .chunk_log2 = 20, .find_all = 0, .quiet = 1,
                         .status_iv = 3600, .ckpt_iv = 1e9, .ckpt_file = NULL,
                         .results_file = "/dev/null", .resume = 0 };
        double t0 = now_s();
        result_t r = run_search(&got);
        double dt = now_s() - t0;
        int ok = (r == R_FOUND && got == expct);
        printf("selftest kind %d k=%-2d: expected %-18s got %-18s  %-4s (%.1fs)\n",
               cases[i].kind, cases[i].k, u128_str(expct, b1),
               r == R_FOUND ? u128_str(got, b2) : "(none)",
               ok ? "PASS" : "FAIL", dt);
        fflush(stdout);
        if (!ok) pass = 0;
    }
    printf(pass ? "selftest: ALL PASS\n" : "selftest: FAILURES\n");
    return pass;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

static void usage(void)
{
    printf(
"cchain — exhaustive search for the smallest Cunningham chain of length k\n"
"         (OEIS A057330: kind 2, p -> 2p-1; kind 1 is p -> 2p+1)\n\n"
"usage: cchain [options]\n"
"  -k, --length N        chain length target (default 17)\n"
"  -K, --kind 1|2        chain kind (default 2)\n"
"  -s, --start P         search start, inclusive (default 3; decimal or 3.2e18)\n"
"  -e, --end P           search end, inclusive (default 1e21)\n"
"  -t, --threads N       worker threads (default: online cores)\n"
"  -B, --sieve-limit Q   sieve primes up to Q (default 1024; larger is usually\n"
"                        slower — GMP PRP tests on survivors are cheap)\n"
"  -x, --chunk-bits N    log2 of t-window per work unit (default 24)\n"
"  -c, --checkpoint F    checkpoint file (default cchain.ckpt; 'none' disables)\n"
"  -i, --ckpt-interval S checkpoint every S seconds (default 60)\n"
"  -u, --status-interval S  status line every S seconds (default 10)\n"
"  -o, --results F       append found chains to F (default found.txt)\n"
"  -r, --resume          resume from checkpoint\n"
"  -A, --find-all        do not stop after first confirmed-minimal chain\n"
"  -q, --quiet           suppress status output\n"
"      --verify P        verify/print the chain starting at P, then exit\n"
"      --verify-known    check all known terms of A057330 (+ kind-1 table)\n"
"      --selftest        rediscover known small terms from scratch\n"
"  -h, --help            this help\n\n"
"Ctrl-C checkpoints and exits; Ctrl-T (SIGINFO) prints an immediate status.\n"
"Exhaustiveness requires start > sieve-limit (the tool warns otherwise).\n");
}

int main(int argc, char **argv)
{
    memset(&S, 0, sizeof S);
    cfg_t cfg = {
        .kind = 2, .k = 17, .start = 3, .end = 0,
        .threads = (int)sysconf(_SC_NPROCESSORS_ONLN),
        .B = 1024, .chunk_log2 = 24, .find_all = 0, .quiet = 0,
        .status_iv = 10, .ckpt_iv = 60,
        .ckpt_file = "cchain.ckpt", .results_file = "found.txt", .resume = 0,
    };
    parse_u128("1e21", &cfg.end);
    const char *verify_arg = NULL;
    int do_selftest = 0, do_known = 0;

    static struct option lo[] = {
        { "length", 1, 0, 'k' }, { "kind", 1, 0, 'K' }, { "start", 1, 0, 's' },
        { "end", 1, 0, 'e' }, { "threads", 1, 0, 't' }, { "sieve-limit", 1, 0, 'B' },
        { "chunk-bits", 1, 0, 'x' }, { "checkpoint", 1, 0, 'c' },
        { "ckpt-interval", 1, 0, 'i' }, { "status-interval", 1, 0, 'u' },
        { "results", 1, 0, 'o' }, { "resume", 0, 0, 'r' }, { "find-all", 0, 0, 'A' },
        { "quiet", 0, 0, 'q' }, { "verify", 1, 0, 1 }, { "verify-known", 0, 0, 2 },
        { "selftest", 0, 0, 3 }, { "help", 0, 0, 'h' }, { 0, 0, 0, 0 } };
    int c;
    while ((c = getopt_long(argc, argv, "k:K:s:e:t:B:x:c:i:u:o:rAqh", lo, NULL)) != -1) {
        switch (c) {
        case 'k': cfg.k = atoi(optarg); break;
        case 'K': cfg.kind = atoi(optarg); break;
        case 's': if (!parse_u128(optarg, &cfg.start)) goto badnum; break;
        case 'e': if (!parse_u128(optarg, &cfg.end)) goto badnum; break;
        case 't': cfg.threads = atoi(optarg); break;
        case 'B': cfg.B = (u32)strtoul(optarg, NULL, 10); break;
        case 'x': cfg.chunk_log2 = atoi(optarg); break;
        case 'c': cfg.ckpt_file = strcmp(optarg, "none") ? optarg : NULL; break;
        case 'i': cfg.ckpt_iv = atof(optarg); break;
        case 'u': cfg.status_iv = atof(optarg); break;
        case 'o': cfg.results_file = optarg; break;
        case 'r': cfg.resume = 1; break;
        case 'A': cfg.find_all = 1; break;
        case 'q': cfg.quiet = 1; break;
        case 1: verify_arg = optarg; break;
        case 2: do_known = 1; break;
        case 3: do_selftest = 1; break;
        case 'h': usage(); return 0;
        default: usage(); return 1;
        badnum: fprintf(stderr, "error: cannot parse number '%s'\n", optarg); return 1;
        }
    }
    if (cfg.kind != 1 && cfg.kind != 2) { fprintf(stderr, "error: kind must be 1 or 2\n"); return 1; }
    if (cfg.k < 4 || cfg.k > MAXK) { fprintf(stderr, "error: k must be in [4,%d]\n", MAXK); return 1; }
    if (cfg.threads < 1) cfg.threads = 1;
    if (cfg.chunk_log2 < 14 || cfg.chunk_log2 > 28) { fprintf(stderr, "error: chunk-bits in [14,28]\n"); return 1; }
    if (cfg.B < 100) { fprintf(stderr, "error: sieve-limit too small\n"); return 1; }
    if (cfg.end <= cfg.start) { fprintf(stderr, "error: end <= start\n"); return 1; }

    if (verify_arg) {
        u128 p;
        if (!parse_u128(verify_arg, &p)) { fprintf(stderr, "bad number\n"); return 1; }
        verify_p(p, cfg.kind, 0);
        return 0;
    }
    if (do_known) { verify_known(); return 0; }
    if (do_selftest) return selftest(cfg.threads) ? 0 : 1;

    signal(SIGINT, on_sigint);
    signal(SIGTERM, on_sigint);
#ifdef SIGINFO
    signal(SIGINFO, on_siginfo);
#endif

    S.cfg = cfg;
    u128 got = 0;
    result_t r = run_search(&got);
    return r == R_INTERRUPTED ? 130 : 0;
}
