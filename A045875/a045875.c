/*
 * a045875.c — search for terms of OEIS A045875:
 *   a(n) = smallest m such that the decimal representation of 2^m
 *          contains n consecutive identical digits.
 *
 * Strategy: keep the full decimal expansion of 2^m as an array of base-10^9
 * limbs (uint32), and double it in place once per step.  When doubling in
 * base 10^9 the carry into limb i is exactly (old a[i-1] >= 5*10^8),
 * independent of anything below, so the pass parallelizes perfectly.
 *
 * Run detection for n >= 17: any run of >= 17 identical digits must contain
 * at least one fully-aligned 9-digit limb, i.e. a limb whose value is
 * d * 111111111 (d = 0..9).  Testing "v % 111111111 == 0" is one multiply
 * and compare, fused into the doubling pass; the rare candidate gets an
 * exact run-length check into its neighbors.
 *
 * For n <= 16 (used to verify against the known terms) a generic digit
 * scan with block overlap is used instead.
 *
 * Init: GMP converts 2^start to decimal once (fast: bit shift + D&C radix
 * conversion), so a search can begin at any exponent.  Checkpoint/resume
 * supported for long runs.
 *
 * Build:  cc -O3 -o a045875 a045875.c -I/opt/homebrew/include \
 *             -L/opt/homebrew/lib -lgmp -lpthread
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdatomic.h>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <time.h>
#include <signal.h>
#include <errno.h>

#if defined(__has_include)
#  if __has_include(<gmp.h>)
#    include <gmp.h>
#    define HAVE_GMP_H 1
#  endif
#endif
#ifndef HAVE_GMP_H
/* Minimal declarations against libgmp's stable ABI (real __gmpz_* symbols;
 * struct layout unchanged for decades).  Lets the tool build on boxes with
 * only the runtime library (libgmp.so.10) and no dev package.  The init
 * self-check below verifies the conversion independently. */
typedef struct { int _mp_alloc; int _mp_size; void *_mp_d; } __mpz_struct;
typedef __mpz_struct mpz_t[1];
typedef unsigned long mp_bitcnt_t;
void __gmpz_init(mpz_t);
void __gmpz_clear(mpz_t);
void __gmpz_set_ui(mpz_t, unsigned long);
void __gmpz_mul_2exp(mpz_t, const mpz_t, mp_bitcnt_t);
char *__gmpz_get_str(char *, int, const mpz_t);
#define mpz_init      __gmpz_init
#define mpz_clear     __gmpz_clear
#define mpz_set_ui    __gmpz_set_ui
#define mpz_mul_2exp  __gmpz_mul_2exp
#define mpz_get_str   __gmpz_get_str
#endif

#define BASE 1000000000u
#define HALF  500000000u
#define REP   111111111u          /* nine ones */
#define REPART_STEPS 4096         /* rebalance partition this often */
#define CAND_CAP 1024
#define CHUNK 8192

/* ---------------- global state ---------------- */
static uint32_t *A;               /* limbs, A[0] least significant */
static size_t cap_limbs;
static _Atomic size_t L_shared;   /* current limb count */

static int T = 8;
static unsigned n_target = 18;
static uint64_t start_m = 0;      /* exponent of initial state */
static uint64_t end_m = 0;        /* 0 = unlimited */
static char *ckpt_path = NULL;
static double ckpt_interval = 600.0;

static _Atomic int stop_now;      /* set only in phase C1, read in C2 */
static _Atomic uint64_t m_last;   /* final exponent, for the exit summary */
static _Atomic int found_any;
static _Atomic int interrupted;
static pthread_mutex_t hit_mu = PTHREAD_MUTEX_INITIALIZER;
static struct {
    int valid;
    uint64_t m;
    unsigned digit;
    size_t len;
    size_t pos_lsd;               /* offset of run's lowest digit from LSD */
} best_hit;

static volatile sig_atomic_t g_sig;
static void on_sig(int s){ (void)s; g_sig = 1; }

static double now_sec(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ---------------- barrier (macOS has no pthread_barrier_t) ---------------- */
static struct { _Atomic int count; _Atomic int sense; char pad[56]; } BAR;
static void bar_wait(int *my_sense){
    int s = 1 - *my_sense; *my_sense = s;
    if (atomic_fetch_add(&BAR.count, 1) == T - 1){
        atomic_store(&BAR.count, 0);
        atomic_store(&BAR.sense, s);
    } else {
        int spins = 0;
        while (atomic_load(&BAR.sense) != s)
            if (++spins > 20000){ sched_yield(); spins = 0; }
    }
}

/* ---------------- helpers ---------------- */
static const uint32_t P10[9] = {1,10,100,1000,10000,100000,
                                1000000,10000000,100000000};

static inline unsigned ndig_u32(uint32_t v){
    unsigned n = 1; while (v >= 10){ v /= 10; n++; } return n;
}
static size_t total_digits(size_t L){
    return 9*(L-1) + ndig_u32(A[L-1]);
}

static void record_hit(uint64_t m, unsigned digit, size_t len, size_t pos_lsd){
    pthread_mutex_lock(&hit_mu);
    if (!best_hit.valid || pos_lsd < best_hit.pos_lsd){
        best_hit.valid = 1; best_hit.m = m; best_hit.digit = digit;
        best_hit.len = len; best_hit.pos_lsd = pos_lsd;
    }
    pthread_mutex_unlock(&hit_mu);
    atomic_store(&found_any, 1);
    atomic_store(&stop_now, 1);
}

/* Exact run length around a fully-repdigit limb i (A[i] == d*REP). */
static size_t run_around(size_t i, size_t L, unsigned d, size_t *pos_lo){
    size_t len = 9;
    size_t lo_digit = i * 9;
    /* extend toward less significant digits */
    size_t j = i;
    while (j > 0){
        uint32_t v = A[j-1];
        unsigned cnt = 0;
        for (int k = 8; k >= 0; k--){
            if ((v / P10[k]) % 10 == (uint32_t)d) cnt++;
            else break;
        }
        len += cnt; lo_digit -= cnt;
        if (cnt < 9) break;
        j--;
    }
    /* extend toward more significant digits */
    j = i + 1;
    while (j < L){
        uint32_t v = A[j];
        unsigned nd = (j == L-1) ? ndig_u32(v) : 9;
        unsigned cnt = 0;
        while (cnt < nd && v % 10 == d){ v /= 10; cnt++; }
        len += cnt;
        if (cnt < nd) break;
        j++;
    }
    *pos_lo = lo_digit;
    return len;
}

/* Generic digit-stream scan of limbs [s,e); records the first run >= n. */
static void scan_range(size_t s, size_t e, size_t L, unsigned n, uint64_t m){
    unsigned cur = 10; size_t len = 0, start_pos = 0;
    int over = 0;
    for (size_t i = s; i < e; i++){
        uint32_t v = A[i];
        unsigned nd = (i == L-1) ? ndig_u32(v) : 9;
        for (unsigned k = 0; k < nd; k++){
            unsigned dig = v % 10; v /= 10;
            if (dig == cur) len++;
            else {
                if (over){ record_hit(m, cur, len, start_pos); return; }
                cur = dig; len = 1; start_pos = i*9 + k;
            }
            if (len >= n) over = 1;
        }
    }
    if (over) record_hit(m, cur, len, start_pos);
}

/* ---------------- checkpoint ---------------- */
struct ckpt_hdr { char magic[8]; uint64_t version, m, L; };

static void write_ckpt(uint64_t m, size_t L){
    if (!ckpt_path) return;
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", ckpt_path);
    FILE *f = fopen(tmp, "wb");
    if (!f){ fprintf(stderr, "checkpoint: cannot open %s: %s\n", tmp, strerror(errno)); return; }
    struct ckpt_hdr h = {{'A','0','4','5','8','7','5','C'}, 1, m, L};
    if (fwrite(&h, sizeof h, 1, f) != 1 ||
        fwrite(A, sizeof(uint32_t), L, f) != L){
        fprintf(stderr, "checkpoint: write failed\n"); fclose(f); return;
    }
    fclose(f);
    if (rename(tmp, ckpt_path) != 0)
        fprintf(stderr, "checkpoint: rename failed: %s\n", strerror(errno));
}

static int load_ckpt(uint64_t *m, size_t *L){
    FILE *f = fopen(ckpt_path, "rb");
    if (!f) return 0;
    struct ckpt_hdr h;
    if (fread(&h, sizeof h, 1, f) != 1 ||
        memcmp(h.magic, "A045875C", 8) != 0 || h.version != 1){
        fprintf(stderr, "checkpoint: bad header in %s\n", ckpt_path);
        exit(1);
    }
    if (h.L > cap_limbs){ fprintf(stderr, "checkpoint larger than buffer\n"); exit(1); }
    if (fread(A, sizeof(uint32_t), h.L, f) != h.L){
        fprintf(stderr, "checkpoint: truncated\n"); exit(1);
    }
    fclose(f);
    *m = h.m; *L = h.L;
    return 1;
}

/* ---------------- worker threads ---------------- */
static void compute_part(int id, size_t part_L, size_t *lo, size_t *hi){
    size_t b = part_L / (size_t)T;
    *lo = b * (size_t)id;
    *hi = (id == T-1) ? part_L : b * (size_t)(id + 1);  /* top: caller overrides */
}

static void *worker(void *argp){
    int id = (int)(intptr_t)argp;
    int sense = 0;
    int is_top = (id == T-1);
    uint64_t m = start_m;
    size_t part_L = atomic_load(&L_shared);
    size_t top_L = part_L;                /* live limb count; top thread only */
    size_t lo, hi;
    compute_part(id, part_L, &lo, &hi);
    if (is_top) hi = top_L;
    uint32_t stash_cin = (lo > 0 && lo <= part_L) ? (A[lo-1] >= HALF) : 0;

    size_t cand[CAND_CAP]; int ncand = 0, cand_over = 0;
    uint64_t steps = 0;
    double t_ck = now_sec(), t_pr = now_sec();
    uint64_t m_pr = m;

    for (;;){
        bar_wait(&sense);                 /* B3: array holds 2^m, stashes valid */

        /* ---- phase U: double in place, collect repdigit-limb candidates ---- */
        size_t uhi = is_top ? top_L : hi;
        ncand = 0; cand_over = 0;
        uint32_t oldtop = is_top ? A[uhi-1] : 0;
        size_t i = uhi;
        while (i > lo + 1){
            size_t c0 = (i > lo + 1 + CHUNK) ? i - CHUNK : lo + 1;
            uint32_t any = 0;
            for (size_t j = i; j-- > c0; ){
                uint32_t nv = 2*A[j] + (A[j-1] >= HALF);
                if (nv >= BASE) nv -= BASE;
                A[j] = nv;
                any |= (nv % REP == 0);
            }
            if (any){
                for (size_t j = c0; j < i; j++)
                    if (A[j] % REP == 0){
                        if (ncand < CAND_CAP) cand[ncand++] = j;
                        else cand_over = 1;
                    }
            }
            i = c0;
        }
        if (uhi > lo){                    /* lowest limb of this block */
            uint32_t nv = 2*A[lo] + stash_cin;
            if (nv >= BASE) nv -= BASE;
            A[lo] = nv;
            if (nv % REP == 0){
                if (ncand < CAND_CAP) cand[ncand++] = lo;
                else cand_over = 1;
            }
        }
        if (is_top && oldtop >= HALF){    /* number grew by one limb */
            A[top_L] = 1;
            top_L++;
            atomic_store(&L_shared, top_L);
        }

        bar_wait(&sense);                 /* B1: array holds 2^(m+1) */

        /* ---- phase C1: examine candidates / scan; set flags ---- */
        m++;
        size_t Lc = atomic_load(&L_shared);
        size_t chi = is_top ? top_L : hi;
        if (n_target >= 17){
            if (cand_over){               /* buffer overflowed: rescan block */
                for (size_t j = lo; j < chi; j++)
                    if (A[j] % REP == 0){
                        size_t pos, len;
                        unsigned d = (unsigned)(A[j] / REP);
                        len = run_around(j, Lc, d, &pos);
                        if (len >= n_target) record_hit(m, d, len, pos);
                    }
            } else {
                for (int c = 0; c < ncand; c++){
                    size_t j = cand[c], pos, len;
                    unsigned d = (unsigned)(A[j] / REP);
                    len = run_around(j, Lc, d, &pos);
                    if (len >= n_target) record_hit(m, d, len, pos);
                }
            }
        } else {
            size_t ext = (size_t)n_target/9 + 1;
            size_t s = (lo > ext) ? lo - ext : 0;
            if (chi > lo) scan_range(s, chi, Lc, n_target, m);
        }
        if (id == 0){
            if (g_sig){ atomic_store(&interrupted, 1); atomic_store(&stop_now, 1); }
            if (end_m && m >= end_m) atomic_store(&stop_now, 1);
            if (Lc + 16 >= cap_limbs){
                fprintf(stderr, "limb buffer nearly full; stopping (restart to continue)\n");
                atomic_store(&interrupted, 1); atomic_store(&stop_now, 1);
            }
        }

        bar_wait(&sense);                 /* B2: flags stable */

        /* ---- phase C2: bookkeeping, decide exit, restash ---- */
        int stop = atomic_load(&stop_now);
        if (id == 0 && stop) atomic_store(&m_last, m);
        if (id == 0){
            double t = now_sec();
            if (t - t_pr >= 5.0){
                fprintf(stderr, "\rm=%llu digits=%zu  %.1f steps/s   ",
                        (unsigned long long)m, total_digits(Lc),
                        (double)(m - m_pr) / (t - t_pr));
                fflush(stderr);
                t_pr = t; m_pr = m;
            }
            if (ckpt_path && (stop || t - t_ck >= ckpt_interval)){
                write_ckpt(m, Lc);
                t_ck = now_sec();
            }
        }
        if (stop) break;                  /* uniform across threads */
        steps++;
        if (steps % REPART_STEPS == 0){
            part_L = Lc;
            compute_part(id, part_L, &lo, &hi);
        }
        stash_cin = (lo > 0) ? (A[lo-1] >= HALF) : 0;
    }
    return NULL;
}

/* ---------------- init ---------------- */
static uint64_t pow2_mod(uint64_t e, uint64_t mod){  /* 2^e mod `mod` */
    unsigned __int128 r = 1, b = 2;
    while (e){ if (e & 1) r = r * b % mod; b = b * b % mod; e >>= 1; }
    return (uint64_t)r;
}

static size_t limbs_for(uint64_t m){
    long double d = (long double)m * 0.30102999566398119521L;
    return (size_t)(d / 9.0L) + 4;
}

static size_t init_from_gmp(uint64_t s){
    fprintf(stderr, "initializing 2^%llu via GMP...\n", (unsigned long long)s);
    double t0 = now_sec();
    mpz_t z; mpz_init(z);
    mpz_set_ui(z, 1);
    mpz_mul_2exp(z, z, (mp_bitcnt_t)s);
    char *str = mpz_get_str(NULL, 10, z);
    mpz_clear(z);
    size_t len = strlen(str);
    size_t L = (len + 8) / 9;
    if (L > cap_limbs){ fprintf(stderr, "start exponent exceeds buffer\n"); exit(1); }
    size_t j = 0;
    for (size_t i = len; i > 0; ){
        size_t st = (i >= 9) ? i - 9 : 0;
        uint32_t v = 0;
        for (size_t k = st; k < i; k++) v = v*10 + (uint32_t)(str[k] - '0');
        A[j++] = v;
        i = st;
    }
    free(str);
    if (A[0] != (uint32_t)pow2_mod(s, BASE)){
        fprintf(stderr, "GMP init self-check FAILED (2^s mod 1e9 mismatch)\n");
        exit(1);
    }
    fprintf(stderr, "init done: %zu digits, %zu limbs (%.1fs), self-check ok\n",
            len, L, now_sec() - t0);
    return L;
}

static void usage(const char *p){
    fprintf(stderr,
"usage: %s [-n runlen] [-s start_exp] [-e end_exp] [-t threads]\n"
"          [-c checkpoint_file] [-i ckpt_interval_sec]\n"
"  Finds smallest m >= start with a run of >= runlen identical digits in 2^m.\n"
"  A045875(18): %s -n 18 -s 356677212 -t 8 -c a18.ckpt\n", p, p);
    exit(2);
}

int main(int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, "n:s:e:t:c:i:h")) != -1){
        switch (opt){
        case 'n': n_target = (unsigned)strtoul(optarg, NULL, 10); break;
        case 's': start_m = strtoull(optarg, NULL, 10); break;
        case 'e': end_m = strtoull(optarg, NULL, 10); break;
        case 't': T = atoi(optarg); break;
        case 'c': ckpt_path = optarg; break;
        case 'i': ckpt_interval = atof(optarg); break;
        default: usage(argv[0]);
        }
    }
    if (n_target < 1 || T < 1) usage(argv[0]);

    uint64_t cap_m = end_m ? end_m + 16 : (start_m > 5000000000ull ?
                                           start_m * 2 : 10000000000ull);
    cap_limbs = limbs_for(cap_m) + 16;
    A = malloc(cap_limbs * sizeof(uint32_t));
    if (!A){ fprintf(stderr, "cannot allocate %zu limbs\n", cap_limbs); return 1; }

    signal(SIGINT, on_sig);
    signal(SIGTERM, on_sig);

    size_t L;
    uint64_t m0 = start_m;
    if (ckpt_path && load_ckpt(&m0, &L)){
        fprintf(stderr, "resumed from %s: m=%llu, %zu digits\n",
                ckpt_path, (unsigned long long)m0, total_digits(L));
        start_m = m0;
    } else if (start_m == 0){
        A[0] = 1; L = 1;
    } else {
        L = init_from_gmp(start_m);
    }
    atomic_store(&L_shared, L);

    /* scan the initial value itself */
    scan_range(0, L, L, n_target, m0);
    if (atomic_load(&found_any)){
        printf("RESULT: smallest m >= %llu with %u consecutive identical digits "
               "in 2^m is m=%llu (digit %u, run length %zu, offset %zu from "
               "least significant digit; 2^m has %zu digits)\n",
               (unsigned long long)start_m, n_target,
               (unsigned long long)best_hit.m, best_hit.digit, best_hit.len,
               best_hit.pos_lsd, total_digits(L));
        return 0;
    }

    if (g_sig){ fprintf(stderr, "interrupted during init\n"); return 130; }
    atomic_store(&BAR.count, 0);
    atomic_store(&BAR.sense, 0);

    fprintf(stderr, "searching: n=%u start=%llu%s threads=%d digits=%zu\n",
            n_target, (unsigned long long)start_m,
            end_m ? " (bounded)" : "", T, total_digits(L));

    pthread_t th[256];
    if (T > 256) T = 256;
    double t_run0 = now_sec();
    for (int t = 1; t < T; t++)
        pthread_create(&th[t], NULL, worker, (void*)(intptr_t)t);
    worker((void*)(intptr_t)0);
    for (int t = 1; t < T; t++)
        pthread_join(th[t], NULL);
    fputc('\n', stderr);                  /* end the \r progress line */
    {
        double el = now_sec() - t_run0;
        uint64_t mf = atomic_load(&m_last);
        if (mf > start_m && el > 0)
            fprintf(stderr, "summary: %d threads, %llu doublings in %.1fs = "
                    "%.1f steps/s\n", T, (unsigned long long)(mf - start_m),
                    el, (double)(mf - start_m) / el);
    }

    if (atomic_load(&found_any)){
        size_t Lf = atomic_load(&L_shared);
        size_t nd = total_digits(Lf);
        double el = now_sec() - t_run0;
        uint64_t steps = best_hit.m - start_m;
        printf("RESULT: smallest m >= %llu with %u consecutive identical digits "
               "in 2^m is m=%llu\n",
               (unsigned long long)start_m, n_target,
               (unsigned long long)best_hit.m);
        printf("==================== FINAL REPORT ====================\n");
        printf("A045875 candidate:   a(%zu) = %llu\n",
               best_hit.len, (unsigned long long)best_hit.m);
        printf("run found:           %zu consecutive '%u' digits\n",
               best_hit.len, best_hit.digit);
        printf("number:              2^%llu  (%zu decimal digits)\n",
               (unsigned long long)best_hit.m, nd);
        printf("run position:        digits %zu..%zu from the least significant "
               "digit\n                     (%zu digits in from the most "
               "significant digit)\n",
               best_hit.pos_lsd, best_hit.pos_lsd + best_hit.len - 1,
               nd - best_hit.pos_lsd - best_hit.len);
        printf("search:              started at m=%llu, %llu doublings, "
               "%.1f hours, avg %.1f steps/s\n",
               (unsigned long long)start_m, (unsigned long long)steps,
               el / 3600.0, steps > 0 && el > 0 ? steps / el : 0.0);
        printf("======================================================\n");
        fflush(stdout);
        return 0;
    }
    if (atomic_load(&interrupted)){
        fprintf(stderr, "interrupted; %s\n",
                ckpt_path ? "state checkpointed" : "no checkpoint file given (-c)");
        return 130;
    }
    fprintf(stderr, "no run of %u found up to m=%llu\n",
            n_target, (unsigned long long)end_m);
    return 1;
}
