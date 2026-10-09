/*
 * a128344 - produce terms of OEIS A128344:
 *   Numbers k such that (7^k - 5^k)/2 is prime.
 *   3, 5, 7, 113, 397, 577, 7573, 14561, 58543, 100019, ...
 *
 * Every term is prime (if k = ab is composite, (7^a - 5^a)/2 > 1 divides
 * (7^k - 5^k)/2 properly), so only prime k are tested. Candidate primes
 * come from primesieve; bignum arithmetic and primality testing from GMP.
 *
 * Per-candidate pipeline:
 *   1. Sieved trial factoring. For prime k, any prime q dividing
 *      (7^k - 5^k)/2 satisfies q == 1 (mod 2k): the order of 7/5 mod q
 *      divides the prime k and cannot be 1 (that would force q | 2), so
 *      it is exactly k and k | q-1; q is odd, hence 2k | q-1. We scan
 *      q = 2km+1: a bitmap sieve over m first removes the ~90% of q with
 *      a prime factor < 2^16 (any relevant prime divisor is itself of the
 *      form 2km'+1 and gets tested at its own m'), then survivors get a
 *      64/128-bit check of 7^k == 5^k (mod q). The depth auto-scales with
 *      k: a bignum PRP test costs ~k^2, so deep trial factoring pays for
 *      itself many times over at large k.
 *   2. Progressive base-2 Fermat test (see fermat_base2) with live
 *      progress; the rare survivor is confirmed with GMP's BPSW test.
 *      Small terms are proven prime; large ones are strong probable
 *      primes (same status as the published OEIS terms).
 *
 * With -F, stage 2 is skipped and trial-factoring survivors are appended
 * to a PFGW-compatible ABC file instead, so the expensive PRP tests can
 * run in PFGW/gwnum on bigger hardware.
 *
 * Work is distributed dynamically: worker threads pull the next candidate
 * index from a shared queue, and the main thread prints confirmed terms
 * (and writes export lines / checkpoints) in sequence order as the
 * completed prefix advances. The checkpoint (-c) records the frontier k
 * below which every candidate is fully resolved, so a killed or crashed
 * run resumes there; work in flight above the frontier is redone.
 */

#include <inttypes.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <gmp.h>
#include <primesieve.h>

#define CAND_BLOCK 256          /* primes generated per queue refill */
#define TF_SIEVE_PMAX 65536     /* sieve q = 2km+1 by primes below this */
#define TF_BLOCK (1u << 21)     /* m values per trial-factor sieve block */
#define TF_MMAX_CAP (1ULL << 33)
#define FERMAT_WBITS 16
#define STATUS_INTERVAL_SEC 1       /* in-place redraw on a terminal   */
#define STATUS_INTERVAL_FILE_SEC 10 /* one line per report when logged */
#define CKPT_MIN_INTERVAL_SEC 5
#define LAST_KNOWN_TERM 636043      /* a(17); nothing else below 10^6 */

enum { ST_PENDING = 0, ST_COMPOSITE, ST_PRIME, ST_SURVIVOR };
enum { PH_TF = 0, PH_FERMAT, PH_CONFIRM };

/* ---- shared state, guarded by g_mu ---- */
static pthread_mutex_t g_mu = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  g_cv = PTHREAD_COND_INITIALIZER;
static uint64_t *g_cand;        /* candidate primes k, in increasing order */
static uint8_t  *g_status;      /* ST_* per candidate                      */
static uint32_t *g_digits;      /* decimal digits of (7^k-5^k)/2 if prime  */
static size_t    g_ncand, g_cap;
static size_t    g_next;        /* next candidate index to hand out        */
static uint64_t  g_last_prime;  /* largest candidate generated so far      */
static uint64_t  g_done_cnt;    /* candidates fully tested (this session)  */
static uint64_t  g_done_maxk;   /* largest candidate fully tested          */
static uint64_t  g_cnt_tf;      /* composites caught by trial factoring    */
static uint64_t  g_cnt_prp;     /* candidates that needed a PRP test       */
static uint64_t  g_surv_cnt;    /* survivors written to the export file    */
static uint64_t  g_frontier_next; /* everything below this k is resolved   */
static bool      g_exhausted;   /* candidate range (-e) fully generated    */
static atomic_bool g_stop;

/* ---- config / run info ---- */
static uint64_t opt_mmax;           /* 0 = auto-scale with k             */
static bool     opt_verbose;
static bool     opt_next;           /* -N: stop at first term found      */
static bool     opt_export;         /* -F: TF only, export survivors     */
static const char *opt_export_path;
static const char *opt_ckpt;        /* -c: checkpoint file               */
static uint64_t g_start = 2;        /* first candidate value (session)   */
static uint64_t g_end;              /* -e: last candidate value (0=none) */
static double   g_t0;               /* search start time                 */

/* cumulative totals carried in from a resumed checkpoint */
static uint64_t g_ck_orig_start = 2;
static uint64_t g_ck_found, g_ck_checked, g_ck_tf, g_ck_prp, g_ck_surv;

/* trial-factor sieve primes (3..TF_SIEVE_PMAX), shared read-only */
static uint32_t *g_sprimes;
static size_t    g_nsprimes;

/* export file state (written only by the printer thread, in order) */
static FILE    *g_export_fp;
static uint64_t g_export_max;   /* largest k already in the file */

/*
 * Per-thread progress of the check in flight, read lock-free by the
 * status thread. During trial factoring bits_done/bits_total is m/mmax;
 * during the Fermat test it is the real position inside the modular
 * exponentiation.
 */
typedef struct {
    _Atomic uint64_t k;          /* candidate being worked on; 0 = idle */
    _Atomic uint64_t bits_done;
    _Atomic uint64_t bits_total;
    _Atomic int      phase;      /* PH_* */
    _Atomic double   t_start;
} prp_prog_t;
static prp_prog_t *g_prog;
static long        g_nthreads;

/*
 * On a terminal the status line is redrawn in place (CR, no LF); anything
 * else printed must clear it first so output doesn't land mid-line.
 */
static atomic_bool g_status_active;

static void status_clear(void)
{
    if (atomic_exchange(&g_status_active, false))
        fputs("\r\x1b[K", stderr);
}

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* ---- 64-bit modular arithmetic (128-bit intermediates) ---- */

static inline uint64_t mulmod_u64(uint64_t a, uint64_t b, uint64_t m)
{
    return (uint64_t)((__uint128_t)a * b % m);
}

static uint64_t powmod_u64(uint64_t b, uint64_t e, uint64_t m)
{
    uint64_t r = 1;
    b %= m;
    while (e) {
        if (e & 1)
            r = mulmod_u64(r, b, m);
        b = mulmod_u64(b, b, m);
        e >>= 1;
    }
    return r;
}

/* inverse of a mod m (m odd prime here, gcd(a,m) = 1), extended Euclid */
static uint64_t invmod_u64(uint64_t a, uint64_t m)
{
    int64_t t = 0, nt = 1;
    uint64_t r = m, nr = a % m;
    while (nr) {
        uint64_t q = r / nr;
        int64_t tmp_t = t - (int64_t)q * nt;  t = nt;  nt = tmp_t;
        uint64_t tmp_r = r - q * nr;          r = nr;  nr = tmp_r;
    }
    return t < 0 ? (uint64_t)(t + (int64_t)m) : (uint64_t)t;
}

/*
 * Trial-factor depth. A PRP test costs ~k^2 bit operations while a
 * trial-division candidate costs O(log k), so the economical depth grows
 * ~quadratically in k. Anchored at m=65536 for small k (where PRP tests
 * are cheap anyway) and capped where marginal returns fade.
 */
static uint64_t tf_mmax(uint64_t k)
{
    if (opt_mmax)
        return opt_mmax;
    long double m = 4096.0L * (long double)k / 1000.0L
                            * (long double)k / 1000.0L;
    if (m < 65536)
        return 65536;
    if (m > (long double)TF_MMAX_CAP)
        return TF_MMAX_CAP;
    return (uint64_t)m;
}

/*
 * Trial-factor (7^k - 5^k)/2 by q = 2km+1, m = 1..mmax.
 * Returns a factor, or 0 if none found (or the search is stopping).
 *
 * m values where q has a prime factor p < TF_SIEVE_PMAX are removed with
 * a bitmap sieve first (q == 0 mod p  <=>  m == -1/(2k) mod p), which
 * eliminates ~90% of the range; only survivors get the two powmods.
 * Composite q are redundant anyway: any prime divisor of the target is
 * itself == 1 (mod 2k) and appears at a smaller m. (For very small k a
 * sieve prime can coincide with q itself and hide a factor; that only
 * costs a PRP test, never a wrong answer.)
 */
static uint64_t trial_factor(uint64_t k, uint64_t mmax, prp_prog_t *pr)
{
    if (k < 32)                 /* value fits easily; GMP handles it */
        return 0;
    uint64_t twok = 2 * k;

    static __thread uint64_t *bmp;    /* TF_BLOCK bits */
    static __thread uint64_t *nextm;  /* per sieve prime: next m to mark */
    if (!bmp) {
        bmp = malloc(TF_BLOCK / 8);
        nextm = malloc(g_nsprimes * sizeof *nextm);
        if (!bmp || !nextm) {
            fprintf(stderr, "a128344: out of memory\n");
            exit(1);
        }
    }

    for (size_t i = 0; i < g_nsprimes; i++) {
        uint32_t p = g_sprimes[i];
        uint64_t t = twok % p;
        /* smallest m >= 1 with 2km+1 == 0 (mod p); none if p | 2k */
        uint64_t m0 = t ? p - invmod_u64(t, p) : UINT64_MAX;
        if (m0 == 0)
            m0 = p;
        /*
         * If the first multiple is q = p itself, q is prime and may be a
         * genuine factor: leave it unsieved and start at the next one.
         */
        if (m0 != UINT64_MAX && twok * m0 + 1 == p)
            m0 += p;
        nextm[i] = m0;
    }

    for (uint64_t base = 1; base <= mmax; base += TF_BLOCK) {
        if (atomic_load(&g_stop))
            return 0;
        uint64_t blk = mmax - base + 1;
        if (blk > TF_BLOCK)
            blk = TF_BLOCK;
        memset(bmp, 0, TF_BLOCK / 8);
        for (size_t i = 0; i < g_nsprimes; i++) {
            uint64_t m = nextm[i];
            if (m == UINT64_MAX)
                continue;
            uint32_t p = g_sprimes[i];
            for (; m < base + blk; m += p)
                if (m >= base)
                    bmp[(m - base) >> 6] |= 1ULL << ((m - base) & 63);
            nextm[i] = m;
        }
        for (uint64_t j = 0; j < blk; j++) {
            if (bmp[j >> 6] & (1ULL << (j & 63)))
                continue;
            uint64_t m = base + j, q;
            if (__builtin_mul_overflow(twok, m, &q) ||
                __builtin_add_overflow(q, 1, &q))
                return 0;
            if (powmod_u64(7, k, q) == powmod_u64(5, k, q))
                return q;       /* q | 7^k - 5^k, and q is odd */
        }
        atomic_store(&pr->bits_done, base + blk - 1);
    }
    return 0;
}

/*
 * Progressive Fermat test: true iff 2^(v-1) == 1 (mod v).
 *
 * This replaces running many Miller-Rabin rounds up front: composites
 * essentially never pass a base-2 Fermat test, so one pass here settles
 * them, and the rare survivor is confirmed with a BPSW test afterwards.
 *
 * The exponent is processed left-to-right in windows, r <- r^(2^c) * 2^w
 * (mod v): the squaring run stays inside GMP's fast mpz_powm, the
 * multiply-by-2^w is a cheap shift because the base is 2, and bits_done
 * advances after every window so the status line can show real progress
 * through the test. 16-bit windows measured best: wide enough to
 * amortize the per-call mpz_powm setup, narrow enough that the shift
 * (up to 65536 bits) stays cheap to reduce.
 */
static bool fermat_base2(const mpz_t v, prp_prog_t *pr)
{
    mpz_t e, r, exp2;
    mpz_inits(e, r, exp2, NULL);
    mpz_sub_ui(e, v, 1);
    size_t nb = mpz_sizeinbase(e, 2);
    atomic_store(&pr->bits_total, nb);
    atomic_store(&pr->bits_done, 0);
    mpz_set_ui(r, 1);
    size_t i = nb;                  /* exponent bits still to process */
    while (i > 0) {
        size_t c = i < FERMAT_WBITS ? i : FERMAT_WBITS;
        if (mpz_cmp_ui(r, 1) != 0) {
            mpz_set_ui(exp2, 0);
            mpz_setbit(exp2, c);    /* exp2 = 2^c */
            mpz_powm(r, r, exp2, v);
        }
        unsigned long w = 0;
        for (size_t j = 0; j < c; j++)
            w = (w << 1) | mpz_tstbit(e, i - 1 - j);
        i -= c;
        if (w) {
            mpz_mul_2exp(r, r, w);
            mpz_mod(r, r, v);
        }
        atomic_store(&pr->bits_done, nb - i);
    }
    bool ok = mpz_cmp_ui(r, 1) == 0;
    mpz_clears(e, r, exp2, NULL);
    return ok;
}

/* ---- candidate queue ---- */

static void extend_candidates_locked(void)
{
    if (g_end && g_last_prime >= g_end) {
        g_exhausted = true;
        return;
    }
    uint64_t *block = primesieve_generate_n_primes(CAND_BLOCK,
                                                   g_last_prime + 1,
                                                   UINT64_PRIMES);
    if (!block) {
        fprintf(stderr, "a128344: primesieve failed\n");
        exit(1);
    }
    size_t take = CAND_BLOCK;
    if (g_end) {
        while (take > 0 && block[take - 1] > g_end)
            take--;
        if (take < CAND_BLOCK)
            g_exhausted = true;
    }
    size_t n = g_ncand + take;
    if (n > g_cap) {
        g_cap = g_cap ? g_cap * 2 : 1024;
        if (g_cap < n)
            g_cap = n;
        g_cand = realloc(g_cand, g_cap * sizeof *g_cand);
        g_status = realloc(g_status, g_cap * sizeof *g_status);
        g_digits = realloc(g_digits, g_cap * sizeof *g_digits);
        if (!g_cand || !g_status || !g_digits) {
            fprintf(stderr, "a128344: out of memory\n");
            exit(1);
        }
    }
    memcpy(g_cand + g_ncand, block, take * sizeof *block);
    memset(g_status + g_ncand, ST_PENDING, take);
    g_ncand = n;
    g_last_prime = block[CAND_BLOCK - 1];
    primesieve_free(block);
}

/* ---- checkpoint ---- */

static void ckpt_write_locked(uint64_t found_total)
{
    if (!opt_ckpt)
        return;
    char tmp[1024];
    snprintf(tmp, sizeof tmp, "%s.tmp", opt_ckpt);
    FILE *f = fopen(tmp, "w");
    if (!f) {
        fprintf(stderr, "a128344: cannot write checkpoint %s\n", tmp);
        return;
    }
    fprintf(f, "a128344-ckpt-v1\n"
            "orig_start=%" PRIu64 "\n"
            "next=%" PRIu64 "\n"
            "found=%" PRIu64 "\n"
            "checked=%" PRIu64 "\n"
            "tf=%" PRIu64 "\n"
            "prp=%" PRIu64 "\n"
            "surv=%" PRIu64 "\n",
            g_ck_orig_start, g_frontier_next, found_total,
            g_ck_checked + g_done_cnt, g_ck_tf + g_cnt_tf,
            g_ck_prp + g_cnt_prp, g_ck_surv + g_surv_cnt);
    fclose(f);
    rename(tmp, opt_ckpt);
}

static bool ckpt_read(const char *path, uint64_t *next)
{
    FILE *f = fopen(path, "r");
    if (!f)
        return false;
    char line[256];
    if (!fgets(line, sizeof line, f) ||
        strncmp(line, "a128344-ckpt-v1", 15) != 0) {
        fprintf(stderr, "a128344: %s is not an a128344 checkpoint\n", path);
        fclose(f);
        exit(1);
    }
    bool have_next = false;
    while (fgets(line, sizeof line, f)) {
        uint64_t val;
        if (sscanf(line, "orig_start=%" SCNu64, &val) == 1)
            g_ck_orig_start = val;
        else if (sscanf(line, "next=%" SCNu64, &val) == 1) {
            *next = val;
            have_next = true;
        } else if (sscanf(line, "found=%" SCNu64, &val) == 1)
            g_ck_found = val;
        else if (sscanf(line, "checked=%" SCNu64, &val) == 1)
            g_ck_checked = val;
        else if (sscanf(line, "tf=%" SCNu64, &val) == 1)
            g_ck_tf = val;
        else if (sscanf(line, "prp=%" SCNu64, &val) == 1)
            g_ck_prp = val;
        else if (sscanf(line, "surv=%" SCNu64, &val) == 1)
            g_ck_surv = val;
    }
    fclose(f);
    if (!have_next) {
        fprintf(stderr, "a128344: checkpoint %s has no next= field\n", path);
        exit(1);
    }
    return true;
}

/* ---- worker ---- */

static void *worker(void *arg)
{
    prp_prog_t *pr = &g_prog[(intptr_t)arg];
    mpz_t a, b, v;
    mpz_inits(a, b, v, NULL);

    for (;;) {
        pthread_mutex_lock(&g_mu);
        while (!atomic_load(&g_stop) && g_next >= g_ncand && !g_exhausted)
            extend_candidates_locked();
        if (atomic_load(&g_stop) || g_next >= g_ncand) {
            pthread_cond_broadcast(&g_cv);
            pthread_mutex_unlock(&g_mu);
            break;
        }
        size_t idx = g_next++;
        uint64_t k = g_cand[idx];
        pthread_mutex_unlock(&g_mu);

        double t0 = now_sec();
        uint8_t st;
        uint32_t digits = 0;
        bool prp_ran = false;

        uint64_t mm = tf_mmax(k);
        atomic_store(&pr->t_start, t0);
        atomic_store(&pr->bits_done, 0);
        atomic_store(&pr->bits_total, mm);
        atomic_store(&pr->phase, PH_TF);
        atomic_store(&pr->k, k);

        uint64_t q = trial_factor(k, mm, pr);
        if (atomic_load(&g_stop)) {
            atomic_store(&pr->k, 0);
            break;              /* TF was cut short; leave k pending */
        }
        if (q) {
            st = ST_COMPOSITE;
            if (opt_verbose) {
                status_clear();
                fprintf(stderr, "  k=%" PRIu64 ": factor %" PRIu64
                        " (%.2fs)\n", k, q, now_sec() - t0);
            }
        } else if (opt_export) {
            st = ST_SURVIVOR;
            if (opt_verbose) {
                status_clear();
                fprintf(stderr, "  k=%" PRIu64 ": TF survivor -> export "
                        "(%.2fs)\n", k, now_sec() - t0);
            }
        } else {
            prp_ran = true;
            mpz_ui_pow_ui(a, 7, k);
            mpz_ui_pow_ui(b, 5, k);
            mpz_sub(v, a, b);
            mpz_fdiv_q_2exp(v, v, 1);
            digits = (uint32_t)mpz_sizeinbase(v, 10);

            atomic_store(&pr->bits_done, 0);
            atomic_store(&pr->bits_total, 0);
            atomic_store(&pr->phase, PH_FERMAT);
            if (!fermat_base2(v, pr)) {
                st = ST_COMPOSITE;
            } else {
                /* rare: an actual prime or a base-2 pseudoprime */
                atomic_store(&pr->phase, PH_CONFIRM);
                st = mpz_probab_prime_p(v, 2) > 0 ? ST_PRIME
                                                  : ST_COMPOSITE;
            }
            if (opt_verbose) {
                status_clear();
                fprintf(stderr, "  k=%" PRIu64 ": %s after PRP test of "
                        "%u-digit number (%.2fs)\n", k,
                        st == ST_PRIME ? "PRIME" : "composite",
                        digits, now_sec() - t0);
            }
        }
        atomic_store(&pr->k, 0);

        pthread_mutex_lock(&g_mu);
        g_status[idx] = st;
        g_digits[idx] = digits;
        if (prp_ran)
            g_cnt_prp++;
        else if (st == ST_COMPOSITE)
            g_cnt_tf++;
        g_done_cnt++;
        if (k > g_done_maxk)
            g_done_maxk = k;
        pthread_cond_broadcast(&g_cv);
        pthread_mutex_unlock(&g_mu);
    }

    mpz_clears(a, b, v, NULL);
    return NULL;
}

/*
 * Periodic progress report on stderr: position in k, rates, what the
 * workers are doing, elapsed. On a terminal the line is redrawn in place
 * with CR (no LF); when stderr is redirected to a file each report is a
 * normal line (at a lower frequency).
 */
static void *status_thread(void *arg)
{
    (void)arg;
    bool tty = isatty(fileno(stderr));
    for (;;) {
        struct timespec ts = {
            tty ? STATUS_INTERVAL_SEC : STATUS_INTERVAL_FILE_SEC, 0
        };
        nanosleep(&ts, NULL);

        pthread_mutex_lock(&g_mu);
        if (atomic_load(&g_stop)) {
            pthread_mutex_unlock(&g_mu);
            break;
        }
        uint64_t done = g_done_cnt;
        uint64_t maxk = g_done_maxk;
        uint64_t disp = g_next ? g_cand[g_next - 1] : g_start;
        uint64_t surv = g_ck_surv + g_surv_cnt;
        pthread_mutex_unlock(&g_mu);

        double el = now_sec() - g_t0;
        double cand_rate = el > 0 ? (double)done / el : 0;
        double k_rate = (el > 0 && maxk > g_start)
                        ? (double)(maxk - g_start) / el : 0;

        /* what's in flight: prefer PRP slots for the lead, else TF */
        int n_tf = 0, n_prp = 0, leadphase = -1;
        double leadfrac = -1, leadstart = 0;
        uint64_t leadk = 0;
        for (long i = 0; i < g_nthreads; i++) {
            uint64_t kk = atomic_load(&g_prog[i].k);
            if (!kk)
                continue;
            int ph = atomic_load(&g_prog[i].phase);
            uint64_t bt = atomic_load(&g_prog[i].bits_total);
            double frac = ph == PH_CONFIRM ? 1.0 : (bt ?
                (double)atomic_load(&g_prog[i].bits_done) / (double)bt : 0);
            if (ph == PH_TF)
                n_tf++;
            else
                n_prp++;
            bool better = leadphase < 0 ||
                (ph != PH_TF && leadphase == PH_TF) ||
                ((ph != PH_TF) == (leadphase != PH_TF) && frac > leadfrac);
            if (better) {
                leadfrac = frac;
                leadk = kk;
                leadphase = ph;
                leadstart = atomic_load(&g_prog[i].t_start);
            }
        }
        char prp[128] = "";
        if (leadphase == PH_CONFIRM) {
            snprintf(prp, sizeof prp, ", TF: %d, PRP: %d, lead k=%" PRIu64
                     " confirming (BPSW)", n_tf, n_prp, leadk);
        } else if (leadphase >= 0) {
            double tel = now_sec() - leadstart;
            double eta = leadfrac > 0.02
                         ? tel * (1 - leadfrac) / leadfrac : -1;
            char eb[32];
            if (eta < 0)
                snprintf(eb, sizeof eb, "?");
            else if (eta < 120)
                snprintf(eb, sizeof eb, "%.0fs", eta);
            else if (eta < 7200)
                snprintf(eb, sizeof eb, "%.1fm", eta / 60);
            else
                snprintf(eb, sizeof eb, "%.1fh", eta / 3600);
            snprintf(prp, sizeof prp, ", TF: %d, PRP: %d, lead k=%" PRIu64
                     " %s%.0f%% (~%s left)", n_tf, n_prp, leadk,
                     leadphase == PH_TF ? "TF " : "", leadfrac * 100, eb);
        }
        char sv[48] = "";
        if (opt_export)
            snprintf(sv, sizeof sv, ", survivors: %" PRIu64, surv);

        char line[384];
        snprintf(line, sizeof line, "status: k=%" PRIu64 " (dispatched to %"
                 PRIu64 "), %" PRIu64 " checked, %.2f cand/s, %.1f k/s%s%s"
                 ", elapsed %.0fs",
                 maxk, disp, done, cand_rate, k_rate, prp, sv, el);
        if (tty) {
            fprintf(stderr, "\r%s\x1b[K", line);
            fflush(stderr);
            atomic_store(&g_status_active, true);
        } else {
            fprintf(stderr, "%s\n", line);
        }
    }
    return NULL;
}

/* ---- export file (-F): PFGW ABC format ---- */

static void export_open(void)
{
    g_export_fp = fopen(opt_export_path, "a+");
    if (!g_export_fp) {
        fprintf(stderr, "a128344: cannot open %s\n", opt_export_path);
        exit(1);
    }
    rewind(g_export_fp);
    char line[128];
    bool empty = true;
    while (fgets(line, sizeof line, g_export_fp)) {
        empty = false;
        if (line[0] >= '0' && line[0] <= '9') {
            uint64_t k = strtoull(line, NULL, 10);
            if (k > g_export_max)
                g_export_max = k;
        }
    }
    if (empty)
        fprintf(g_export_fp, "ABC (7^$a-5^$a)/2\n");
    fflush(g_export_fp);
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "usage: %s [options] nterms\n"
        "       %s -N [options]\n"
        "       %s -F file [options]\n"
        "  produce terms of OEIS A128344: k such that (7^k-5^k)/2 is prime\n"
        "  -N          search for the next term only: stop at the first\n"
        "              hit and print it as 'HIT: k' (default start is\n"
        "              %u, just past the last known term)\n"
        "  -F file     no PRP tests: append trial-factoring survivors to\n"
        "              file in PFGW ABC format (run: pfgw -f0 file)\n"
        "  -t threads  worker threads (default: online CPU count)\n"
        "  -s start    search candidates k >= start (default 2)\n"
        "  -e end      stop after all candidates k <= end are resolved\n"
        "  -c file     checkpoint file: progress is saved there and a\n"
        "              killed/crashed run resumes from it\n"
        "  -d mmax     trial-factor depth, q = 2km+1 for m <= mmax\n"
        "              (default: auto-scales with k, up to 2^33)\n"
        "  -v          log per-candidate progress to stderr\n",
        prog, prog, prog, LAST_KNOWN_TERM + 1);
    exit(2);
}

int main(int argc, char **argv)
{
    long nthreads = sysconf(_SC_NPROCESSORS_ONLN);
    uint64_t start = 2;
    bool start_given = false;
    int c;

    while ((c = getopt(argc, argv, "NF:t:d:s:e:c:v")) != -1) {
        switch (c) {
        case 'N': opt_next = true; break;
        case 'F': opt_export = true; opt_export_path = optarg; break;
        case 't': nthreads = strtol(optarg, NULL, 10); break;
        case 'd': opt_mmax = strtoull(optarg, NULL, 10); break;
        case 's': start = strtoull(optarg, NULL, 10);
                  start_given = true; break;
        case 'e': g_end = strtoull(optarg, NULL, 10); break;
        case 'c': opt_ckpt = optarg; break;
        case 'v': opt_verbose = true; break;
        default:  usage(argv[0]);
        }
    }

    uint64_t target;
    if (opt_export) {
        if (opt_next || optind != argc)
            usage(argv[0]);        /* -F takes no nterms and excludes -N */
        target = UINT64_MAX;
    } else if (opt_next) {
        if (optind != argc)
            usage(argv[0]);        /* nterms makes no sense with -N */
        target = 1;
        if (!start_given && !opt_ckpt)
            start = LAST_KNOWN_TERM + 1;
    } else {
        if (optind != argc - 1)
            usage(argv[0]);
        target = strtoull(argv[optind], NULL, 10);
        if (target == 0)
            return 0;
    }
    if (nthreads < 1)
        nthreads = 1;

    /* a checkpoint, if present, overrides -s (and -N's default start) */
    bool resumed = false;
    if (opt_ckpt) {
        uint64_t next;
        if (ckpt_read(opt_ckpt, &next)) {
            resumed = true;
            if (start_given && start != next)
                fprintf(stderr, "a128344: resuming from checkpoint "
                        "(next k=%" PRIu64 "), ignoring -s %" PRIu64 "\n",
                        next, start);
            start = next;
        }
    }
    if (opt_next && !start_given && !resumed)
        start = LAST_KNOWN_TERM + 1;
    if (start < 2)
        start = 2;
    if (!resumed)
        g_ck_orig_start = start;
    if (g_end && g_end < start) {
        fprintf(stderr, "a128344: range [%" PRIu64 ", %" PRIu64 "] is "
                "empty%s\n", start, g_end,
                resumed ? " (already completed per checkpoint)" : "");
        return 0;
    }
    if (!opt_export && g_ck_found >= target) {
        fprintf(stderr, "a128344: checkpoint says %" PRIu64 " term%s "
                "already found; nothing to do\n", g_ck_found,
                g_ck_found == 1 ? "" : "s");
        return 0;
    }

    if (resumed)
        fprintf(stderr, "a128344: resumed from %s: continuing at k=%"
                PRIu64 " (%" PRIu64 " checked, %" PRIu64 " found so far)\n",
                opt_ckpt, start, g_ck_checked, g_ck_found);
    if (opt_next)
        fprintf(stderr, "a128344: searching for the next term, starting "
                "at k=%" PRIu64 "%s\n", start,
                start == LAST_KNOWN_TERM + 1 ? " (past a(17) = 636043)"
                                             : "");
    else if (!opt_export && target > 17 && g_ck_orig_start <= 2)
        fprintf(stderr,
            "a128344: note: only 17 terms are known (a(17) = 636043; no "
            "others below 10^6);\nterms beyond that are unexplored "
            "territory and may take a very long time.\n");

    if (opt_export)
        export_open();

    /* sieve primes for trial factoring */
    size_t nsp;
    g_sprimes = primesieve_generate_primes(3, TF_SIEVE_PMAX, &nsp,
                                           UINT32_PRIMES);
    if (!g_sprimes) {
        fprintf(stderr, "a128344: primesieve failed\n");
        return 1;
    }
    g_nsprimes = nsp;

    g_start = start;
    g_last_prime = start - 1;
    g_frontier_next = start;
    g_t0 = now_sec();
    double t0 = g_t0;

    pthread_t *tid = malloc((size_t)nthreads * sizeof *tid);
    g_prog = calloc((size_t)nthreads, sizeof *g_prog);
    g_nthreads = nthreads;
    if (!tid || !g_prog) {
        fprintf(stderr, "a128344: out of memory\n");
        return 1;
    }
    for (long i = 0; i < nthreads; i++) {
        if (pthread_create(&tid[i], NULL, worker, (void *)(intptr_t)i)) {
            fprintf(stderr, "a128344: pthread_create failed\n");
            return 1;
        }
    }
    pthread_t status_tid;
    if (pthread_create(&status_tid, NULL, status_thread, NULL)) {
        fprintf(stderr, "a128344: pthread_create failed\n");
        return 1;
    }

    /*
     * Print confirmed terms (and write export lines / checkpoints) in
     * order as the completed prefix advances. The export file is written
     * strictly in k order, so on resume anything <= the largest k already
     * in the file is a duplicate and is skipped.
     */
    uint64_t printed = 0;
    size_t pos = 0;
    double last_ckpt = 0;
    pthread_mutex_lock(&g_mu);
    for (;;) {
        bool advanced = false;
        while (pos < g_ncand && g_status[pos] != ST_PENDING &&
               g_ck_found + printed < target) {
            uint64_t k = g_cand[pos];
            if (g_status[pos] == ST_PRIME) {
                printed++;
                status_clear();
                if (opt_next) {
                    double el = now_sec() - t0;
                    printf("HIT: %" PRIu64 "\n", k);
                    fflush(stdout);
                    fprintf(stderr,
                        "a128344: HIT at k=%" PRIu64 " after %.1fs: "
                        "(7^k-5^k)/2 is a %" PRIu32 "-digit "
                        "(probable) prime\n"
                        "  searched k=%" PRIu64 "..%" PRIu64 ": %" PRIu64
                        " candidates checked (%" PRIu64 " eliminated by "
                        "trial factoring, %" PRIu64 " PRP-tested)\n"
                        "  session rates: %.2f cand/s, %.1f k/s\n",
                        k, el, g_digits[pos],
                        g_ck_orig_start, k,
                        g_ck_checked + g_done_cnt, g_ck_tf + g_cnt_tf,
                        g_ck_prp + g_cnt_prp,
                        el > 0 ? (double)g_done_cnt / el : 0,
                        (el > 0 && k > g_start)
                            ? (double)(k - g_start) / el : 0);
                } else {
                    printf("%" PRIu64 "\n", k);
                }
                fflush(stdout);
            } else if (g_status[pos] == ST_SURVIVOR && k > g_export_max) {
                fprintf(g_export_fp, "%" PRIu64 "\n", k);
                fflush(g_export_fp);
                g_surv_cnt++;
            }
            pos++;
            g_frontier_next = k + 1;
            advanced = true;
        }
        if (advanced && now_sec() - last_ckpt > CKPT_MIN_INTERVAL_SEC) {
            ckpt_write_locked(g_ck_found + printed);
            last_ckpt = now_sec();
        }
        if (g_ck_found + printed >= target)
            break;
        if (g_exhausted && pos >= g_ncand) {
            g_frontier_next = g_end + 1;
            break;
        }
        pthread_cond_wait(&g_cv, &g_mu);
    }
    atomic_store(&g_stop, true);
    ckpt_write_locked(g_ck_found + printed);
    pthread_mutex_unlock(&g_mu);

    status_clear();
    double el = now_sec() - t0;
    if (opt_export)
        fprintf(stderr, "a128344: range done at k=%" PRIu64 ": %" PRIu64
                " survivor%s in %s (%" PRIu64 " total), %.2fs (%ld "
                "thread%s)\n", g_end, g_surv_cnt,
                g_surv_cnt == 1 ? "" : "s", opt_export_path,
                g_ck_surv + g_surv_cnt, el,
                nthreads, nthreads == 1 ? "" : "s");
    else if (g_ck_found + printed < target)
        fprintf(stderr, "a128344: range exhausted at k=%" PRIu64 ": "
                "%" PRIu64 " term%s this session (%" PRIu64 " total) in "
                "%.2fs (%ld thread%s)%s\n", g_end, printed,
                printed == 1 ? "" : "s", g_ck_found + printed, el,
                nthreads, nthreads == 1 ? "" : "s",
                opt_next && printed == 0 ? "; no HIT in range" : "");
    else
        fprintf(stderr, "a128344: %" PRIu64 " term%s in %.2fs (%ld "
                "thread%s)\n", printed, printed == 1 ? "" : "s", el,
                nthreads, nthreads == 1 ? "" : "s");

    /*
     * Workers may still be deep inside a long test on a candidate beyond
     * the frontier; those results can no longer change the output (and
     * the checkpoint frontier is already saved), so exit rather than
     * wait for them.
     */
    exit(0);
}
