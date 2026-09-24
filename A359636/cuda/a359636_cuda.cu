/*
 * a359636_cuda.cu  --  GPU (CUDA) version of the exhaustive scan of a359636.c
 *
 * Same problem, same reduction, same output format as a359636.c (read that
 * header first).  The sieve over m = 6k+3 (three targets m-1, m, m+1, one
 * byte counter each packed in a 32-bit word per k) runs on the GPU: one
 * thread block per segment of SEGK consecutive k, the counters in shared
 * memory, the primes 5..23 applied from two periodic pattern tables and the
 * primes 29..T by strided shared-memory atomicAdd (warp per progression for
 * p <= PSMALL, thread per progression above).  Positions whose three counters
 * all reach the threshold are appended to a survivor list, which CPU threads
 * verify exactly (omega(m-1), omega(m), omega(m+1) >= n by trial division with
 * the size-bound early exit, then the prime gap with Miller-Rabin) -- the same
 * host code as the CPU tool.  Because the GPU has cycles to spare, T defaults
 * to the smallest bound that forces six (not five) small prime factors into
 * each target, which makes survivors about ten times rarer.
 *
 * Usage
 *   a359636_cuda scan N [START] END [-t CPUTHREADS] [-T SIEVEMAX] [-b BLOCKS] [-S STATE] [-i SECS] [-a] [-q]
 *   a359636_cuda selftest [-8]
 *   a359636_cuda verify N P
 *
 * Build (DGX Spark / GB10, CUDA 13; no dependencies beyond CUDA):
 *   nvcc -O3 -std=c++17 -arch=sm_121 -Xcompiler -pthread a359636_cuda.cu -o a359636_cuda
 * Measured on the GB10: 810 Gm/s at level 9 (T = 2000), 904 Gm/s at level 10 (T = 1000),
 * versus 82 Gm/s for the CPU tool on the Spark's 20 Arm cores.
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cstdarg>
#include <cctype>
#include <cmath>
#include <csignal>
#include <cerrno>
#include <ctime>
#include <cinttypes>
#include <atomic>
#include <vector>
#include <deque>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <unistd.h>
#include <cuda_runtime.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint8_t  u8;
typedef unsigned __int128 u128;

#ifndef SEGK
#define SEGK     8192           /* k per segment (block): 32 KB of shared u32 (16384 is 2x slower: occupancy) */
#endif
#ifndef BLOCKDIM
#define BLOCKDIM 256
#endif
#define PSMALL   127            /* primes <= PSMALL: one warp per progression */
#define PAT_A    5005           /* 5*7*11*13 */
#define PAT_B    7429           /* 17*19*23  */
#define RING     (1 << 16)
#define MAXSOL   4096
#define NKNOWN   8

static const u64 KNOWN[NKNOWN + 1] = { 0, 7, 19, 643, 51427, 8083633, 1077940147ULL,
                                       75582271489ULL, 34710483181813ULL };

static bool g_quiet = false;
static bool stderr_tty;
static volatile sig_atomic_t g_stop = 0;
static std::mutex out_mu;

#define CUDA_CHECK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) die("CUDA error %s at %s:%d: %s", #x, __FILE__, __LINE__, cudaGetErrorString(e_)); } while (0)

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
    va_start(ap, fmt);
    fputs("a359636_cuda: ", stderr);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
    exit(1);
}

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) die("out of memory");
    return p;
}

static u64 parse_num(const char *s)
{
    const char *plus = strchr(s + 1, '+'), *minus = strchr(s + 1, '-');
    if (plus || minus) {
        const char *op = plus ? plus : minus;
        char *a = strndup(s, (size_t)(op - s));
        u64 x = parse_num(a), y = parse_num(op + 1);
        free(a);
        return plus ? x + y : x - y;
    }
    char *end;
    errno = 0;
    unsigned long long base = strtoull(s, &end, 10);
    if (end == s) die("bad number '%s'", s);
    if (*end == '\0') return base;
    if (*end == '^' || *end == 'e' || *end == 'E') {
        bool pw = (*end == '^');
        unsigned long long e = strtoull(end + 1, &end, 10);
        if (*end != '\0') die("bad number '%s'", s);
        u128 r = 1, b = pw ? base : 10;
        for (unsigned long long i = 0; i < e; i++) {
            r *= b;
            if (r >> 64) die("number '%s' does not fit in 64 bits", s);
        }
        if (!pw) r *= base;
        if (r >> 64) die("number '%s' does not fit in 64 bits", s);
        return (u64)r;
    }
    die("bad number '%s'", s);
    return 0;
}

static void fmt_commas(char *out, u64 v)
{
    char tmp[32];
    int n = snprintf(tmp, sizeof tmp, "%" PRIu64, v), o = 0;
    for (int i = 0; i < n; i++) {
        if (i && (n - i) % 3 == 0) out[o++] = ',';
        out[o++] = tmp[i];
    }
    out[o] = 0;
}

static void fmt_time(char *out, double s)
{
    if (s < 0 || s > 1e9) { strcpy(out, "?"); return; }
    int t = (int)s;
    if (t >= 86400) sprintf(out, "%dd%02dh%02dm", t / 86400, t % 86400 / 3600, t % 3600 / 60);
    else if (t >= 3600) sprintf(out, "%dh%02dm%02ds", t / 3600, t % 3600 / 60, t % 60);
    else if (t >= 60) sprintf(out, "%dm%02ds", t / 60, t % 60);
    else sprintf(out, "%ds", t);
}

static void on_sigint(int) { g_stop = 1; }

/* ------------------------------------------------------------------ */
/* Arithmetic (host)                                                   */
/* ------------------------------------------------------------------ */

static inline u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }

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

static bool is_prime(u64 n)
{
    static const u64 bases[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
    if (n < 2) return false;
    for (int i = 0; i < 12; i++) if (n % bases[i] == 0) return n == bases[i];
    if (n < 41 * 41) return true;
    u64 d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        u64 x = powmod(bases[i], d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 1; r < s && comp; r++) {
            x = mulmod(x, x, n);
            if (x == n - 1) comp = false;
        }
        if (comp) return false;
    }
    return true;
}

static u64 iroot(u64 n, int r)
{
    if (r == 1) return n;
    if (n < 2) return n;
    u64 x = (u64)pow((double)n, 1.0 / r);
    if (x == 0) x = 1;
    for (;;) {
        u128 p = 1;
        bool over = false;
        for (int i = 0; i < r; i++) { p *= x; if (p > n) { over = true; break; } }
        if (over) { x--; continue; }
        u128 q = 1;
        over = false;
        for (int i = 0; i < r; i++) { q *= (x + 1); if (q > n) { over = true; break; } }
        if (!over) { x++; continue; }
        return x;
    }
}

static bool is_perfect_power(u64 n)
{
    for (int j = 2; j <= 63; j++) {
        u64 x = iroot(n, j);
        if (x < 2) return false;
        u128 p = 1;
        for (int i = 0; i < j; i++) p *= x;
        if (p == n) return true;
    }
    return false;
}

static u32 *g_primes;
static size_t g_nprimes;
static u64 g_plimit;

static void build_primes(u64 limit)          /* plain Eratosthenes; limit is a few million at most */
{
    if (limit > (1ull << 31)) die("prime table limit too large");
    std::vector<u8> comp(limit + 1, 0);
    std::vector<u32> pr;
    for (u64 i = 2; i <= limit; i++) {
        if (comp[i]) continue;
        pr.push_back((u32)i);
        for (u64 j = i * i; j <= limit; j += i) comp[j] = 1;
    }
    free(g_primes);
    g_primes = (u32 *)xmalloc(pr.size() * sizeof(u32));
    memcpy(g_primes, pr.data(), pr.size() * sizeof(u32));
    g_nprimes = pr.size();
    g_plimit = limit;
}

static bool omega_at_least(u64 N, int n)
{
    if (n <= 0) return true;
    if (N < 2) return false;
    int cnt = 0;
    u64 lim = iroot(N, n);
    for (size_t i = 0; i < g_nprimes; i++) {
        u64 p = g_primes[i];
        int r = n - cnt;
        if (r == 1) return N > 1;
        if (p > lim) return false;
        if (r == 2 && p * p * p > N) {
            if (is_prime(N)) return false;
            return !is_perfect_power(N);
        }
        if (N % p == 0) {
            cnt++;
            N /= p;
            while (N % p == 0) N /= p;
            r = n - cnt;
            if (r <= 0) return true;
            if (r == 1) return N > 1;
            lim = iroot(N, r);
        }
    }
    if (N == 1) return cnt >= n;
    if (is_prime(N)) return cnt + 1 >= n;
    return cnt + (is_perfect_power(N) ? 1 : 2) >= n;
}

static u64 rho(u64 n)
{
    if (!(n & 1)) return 2;
    for (u64 c = 1;; c++) {
        u64 x = 2, y = 2, d = 1;
        while (d == 1) {
            x = (mulmod(x, x, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            u64 diff = x > y ? x - y : y - x;
            u64 a = diff, b = n;
            while (b) { u64 t = a % b; a = b; b = t; }
            d = a;
        }
        if (d != n) return d;
    }
}

static int factor_full(u64 n, u64 *pr, int *ex)
{
    int k = 0;
    for (size_t i = 0; i < g_nprimes && (u64)g_primes[i] * g_primes[i] <= n; i++) {
        u64 p = g_primes[i];
        if (n % p == 0) {
            pr[k] = p; ex[k] = 0;
            while (n % p == 0) { n /= p; ex[k]++; }
            k++;
        }
    }
    u64 stack[64];
    int sp = 0;
    if (n > 1) stack[sp++] = n;
    while (sp) {
        u64 x = stack[--sp];
        if (x == 1) continue;
        if (is_prime(x)) {
            int j;
            for (j = 0; j < k; j++) if (pr[j] == x) { ex[j]++; break; }
            if (j == k) { pr[k] = x; ex[k] = 1; k++; }
            continue;
        }
        u64 d = rho(x);
        stack[sp++] = d;
        stack[sp++] = x / d;
    }
    for (int i = 1; i < k; i++)
        for (int j = i; j > 0 && pr[j - 1] > pr[j]; j--) {
            u64 tp = pr[j]; pr[j] = pr[j - 1]; pr[j - 1] = tp;
            int te = ex[j]; ex[j] = ex[j - 1]; ex[j - 1] = te;
        }
    return k;
}

static void fmt_factorization(char *out, u64 n)
{
    u64 pr[64]; int ex[64];
    int k = factor_full(n, pr, ex);
    int o = 0;
    for (int i = 0; i < k; i++) {
        o += sprintf(out + o, "%s%" PRIu64, i ? "*" : "", pr[i]);
        if (ex[i] > 1) o += sprintf(out + o, "^%d", ex[i]);
    }
    if (k == 0) sprintf(out + o, "%" PRIu64, n);
}

static int omega_exact(u64 n) { u64 pr[64]; int ex[64]; return factor_full(n, pr, ex); }

static u64 prev_prime(u64 x) { x--; if (!(x & 1)) x--; while (!is_prime(x)) x -= 2; return x; }
static u64 next_prime(u64 x) { x++; if (!(x & 1)) x++; while (!is_prime(x)) x += 2; return x; }

static bool gap_qualifies(u64 p, int n, u64 *q_out)
{
    u64 q = next_prime(p);
    *q_out = q;
    if (q - p < 4) return false;
    for (u64 x = p + 1; x < q; x++) if (!omega_at_least(x, n)) return false;
    return true;
}

static u64 inv_mod(u64 a, u64 p) { return powmod(a % p, p - 2, p); }

static u32 target_residue(u64 p, int t)          /* k with p | 6k+3+t */
{
    u64 inv6 = inv_mod(6, p);
    u64 c = (u64)(3 + t) % p;
    return (u32)(((p - c) % p) * inv6 % p);
}

static int compute_c(int n, u32 T, u64 Bmax)
{
    int rmax = 0;
    for (int r = 0; r <= n; r++) {
        u128 prod = 1;
        bool over = false;
        for (int i = 0; i < n - r; i++) { prod *= g_primes[i]; if (prod > Bmax) { over = true; break; } }
        size_t j = 0;
        while (j < g_nprimes && g_primes[j] <= T) j++;
        for (int i = 0; i < r && !over; i++) {
            if (j + i >= g_nprimes) die("prime table too small in compute_c");
            prod *= g_primes[j + i];
            if (prod > Bmax) over = true;
        }
        if (over) break;
        rmax = r;
    }
    return n - rmax;
}

/* ------------------------------------------------------------------ */
/* GPU sieve                                                           */
/* ------------------------------------------------------------------ */

struct Prog { u32 p, r, inc; };            /* prime, residue of k, packed increment */

__global__ void __launch_bounds__(BLOCKDIM)
sieve_kernel(u64 kbase, u64 klimit, const Prog *__restrict__ prog, const u32 *__restrict__ kb,
             int nprog, int nsmall, const u32 *__restrict__ patA, const u32 *__restrict__ patB,
             u32 thradd, u32 thrmask, u64 *__restrict__ out, unsigned int *__restrict__ outcount, unsigned int outcap)
{
    extern __shared__ u32 W[];
    const int tid = threadIdx.x;
    const u64 k0 = kbase + (u64)blockIdx.x * SEGK;
    if (k0 >= klimit) return;

    /* pattern init for the primes 5..23: W[j] = patA[(k0+j) % 5005] + patB[(k0+j) % 7429]
       (coalesced: consecutive threads take consecutive positions) */
#ifdef NOPAT
    for (int j = tid; j < SEGK; j += BLOCKDIM) W[j] = 0;
#else
    {
        /* running table indices: thread tid takes positions tid, tid+BLOCKDIM, ... */
        u32 a = (u32)(k0 % PAT_A) + tid, b = (u32)(k0 % PAT_B) + tid;
        if (a >= PAT_A) a -= PAT_A;
        if (b >= PAT_B) b -= PAT_B;
        for (int j = tid; j < SEGK; j += BLOCKDIM) {
            W[j] = __ldg(patA + a) + __ldg(patB + b);
            a += BLOCKDIM; if (a >= PAT_A) a -= PAT_A;
            b += BLOCKDIM; if (b >= PAT_B) b -= PAT_B;
        }
    }
#endif
    __syncthreads();

#ifdef NOATOM
#define MARK(j, inc) (W[j] += (inc))
#else
#define MARK(j, inc) atomicAdd(&W[j], (inc))
#endif
#ifndef NOMARK
    /* phase A: small primes, one warp per progression (lanes hit distinct banks: p odd) */
    {
        const int warp = tid >> 5, lane = tid & 31, nwarps = BLOCKDIM >> 5;
        for (int q = warp; q < nsmall; q += nwarps) {
            const u32 p = prog[q].p, inc = prog[q].inc, r = prog[q].r;
            const u32 k0m = (kb[q] + (u32)blockIdx.x * (SEGK % p)) % p;      /* < 2^31 for blocks < 65536 */
            u32 j = r >= k0m ? r - k0m : r + p - k0m;
            for (j += lane * p; j < SEGK; j += 32 * p) MARK(j, inc);
        }
    }
    /* phase B: larger primes, one thread per progression, snake order for balance */
    {
        const int nlarge = nprog - nsmall;
        for (int s = 0; s * BLOCKDIM < nlarge; s++) {
            const int idx = s * BLOCKDIM + ((s & 1) ? (BLOCKDIM - 1 - tid) : tid);
            if (idx >= nlarge) continue;
            const int q = nsmall + idx;
            const u32 p = prog[q].p, inc = prog[q].inc, r = prog[q].r;
            const u32 k0m = (kb[q] + (u32)blockIdx.x * (SEGK % p)) % p;
            u32 j = r >= k0m ? r - k0m : r + p - k0m;
            for (; j < SEGK; j += p) MARK(j, inc);
        }
    }
#endif
    __syncthreads();

    /* scan: all three counters >= thr  <=>  (~(W + thradd)) & 0x808080 == 0 */
    for (int j = tid; j < SEGK; j += BLOCKDIM) {
        if ((((~(W[j] + thradd)) & thrmask)) == 0) {
            u64 k = k0 + j;
            if (k < klimit) {
                unsigned int idx = atomicAdd(outcount, 1u);
                if (idx < outcap) out[idx] = k;
            }
        }
    }
}

/* ------------------------------------------------------------------ */
/* Scan driver                                                         */
/* ------------------------------------------------------------------ */

struct Solution { u64 p, q, m; };

struct Batch { u64 idx; std::vector<u64> ks; };

struct Scan {
    int n;
    u64 k_lo, k_hi;
    u64 lsize;                  /* k per launch */
    u64 nlaunch;
    u32 T; int c, thr;
    std::vector<Prog> prog; int nsmall;
    /* device */
    Prog *d_prog; u32 *d_kb; u32 *d_patA, *d_patB; u64 *d_out; unsigned int *d_count;
    unsigned int outcap;
    /* stats */
    std::atomic<unsigned long long> survivors{0}, triples{0}, nsol{0};
    std::vector<Solution> sol; std::mutex sol_mu;
    u64 best_p = 0, best_m = 0;
    /* verification queue */
    std::deque<Batch> queue; std::mutex q_mu; std::condition_variable q_cv; bool q_done = false;
    size_t q_pending_k = 0;
    /* ordered completion of launches */
    u64 frontier = 0; std::vector<u8> done; std::mutex fr_mu;
    bool scan_all = false, stop_now = false;
    double t0 = 0;
    const char *state_file = nullptr;
    u64 resume_frontier = 0;
};

static void record_solution(Scan *S, u64 p, u64 q, u64 m)
{
    std::lock_guard<std::mutex> lk(S->sol_mu);
    for (auto &s : S->sol) if (s.p == p) return;
    if (S->sol.size() < MAXSOL) S->sol.push_back({p, q, m});
    S->nsol++;
    if (!S->best_p || p < S->best_p) { S->best_p = p; S->best_m = m; }
    std::lock_guard<std::mutex> lo(out_mu);
    if (stderr_tty && !g_quiet) fputs("\r\033[K", stderr);
    printf("SOLUTION n=%d p=%" PRIu64 " q=%" PRIu64 " gap=%" PRIu64 " (m=%" PRIu64 ")\n", S->n, p, q, q - p, m);
    char buf[256];
    for (u64 x = p + 1; x < q; x++) {
        fmt_factorization(buf, x);
        printf("    %" PRIu64 " = %s  omega=%d\n", x, buf, omega_exact(x));
    }
    fflush(stdout);
}

static void check_survivor(Scan *S, u64 k)
{
    u64 m = 6 * k + 3;
    int n = S->n;
    if (!omega_at_least(m - 1, n)) return;
    if (!omega_at_least(m, n)) return;
    if (!omega_at_least(m + 1, n)) return;
    S->triples++;
    u64 p = prev_prime(m - 1), q;
    {
        std::lock_guard<std::mutex> lo(out_mu);
        if (stderr_tty && !g_quiet) fputs("\r\033[K", stderr);
        printf("triple   n=%d m=%" PRIu64 "  omega(m-1,m,m+1) = %d %d %d  prevprime=%" PRIu64 " (%+" PRId64 ")\n",
               n, m, omega_exact(m - 1), omega_exact(m), omega_exact(m + 1), p, (int64_t)(p - m));
        fflush(stdout);
    }
    if (gap_qualifies(p, n, &q)) record_solution(S, p, q, m);
}

static void launch_done(Scan *S, u64 idx)
{
    std::lock_guard<std::mutex> lk(S->fr_mu);
    S->done[idx % RING] = 1;
    while (S->frontier < S->nlaunch && S->done[S->frontier % RING]) {
        S->done[S->frontier % RING] = 0;
        S->frontier++;
    }
    if (S->best_p && !S->scan_all) {
        u64 kbest = (S->best_m - 3) / 6;
        u64 lbest = (kbest - S->k_lo) / S->lsize;
        if (S->frontier > lbest) S->stop_now = true;
    }
}

static void verifier(Scan *S)
{
    for (;;) {
        Batch b;
        {
            std::unique_lock<std::mutex> lk(S->q_mu);
            S->q_cv.wait(lk, [&] { return !S->queue.empty() || S->q_done; });
            if (S->queue.empty()) return;
            b = std::move(S->queue.front());
            S->queue.pop_front();
            S->q_pending_k -= b.ks.size();
        }
        S->q_cv.notify_all();
        for (u64 k : b.ks) check_survivor(S, k);
        launch_done(S, b.idx);
    }
}

static void save_state(Scan *S)
{
    if (!S->state_file) return;
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", S->state_file);
    FILE *f = fopen(tmp, "w");
    if (!f) { fprintf(stderr, "a359636_cuda: cannot write %s\n", tmp); return; }
    u64 fr;
    { std::lock_guard<std::mutex> lk(S->fr_mu); fr = S->frontier; }
    fprintf(f, "A359636 cuda state v1\n");
    fprintf(f, "n %d\nk_lo %" PRIu64 "\nk_hi %" PRIu64 "\nlsize %" PRIu64 "\nfrontier %" PRIu64 "\n",
            S->n, S->k_lo, S->k_hi, S->lsize, fr);
    fprintf(f, "survivors %llu\ntriples %llu\n", (unsigned long long)S->survivors.load(), (unsigned long long)S->triples.load());
    {
        std::lock_guard<std::mutex> lk(S->sol_mu);
        for (auto &s : S->sol) fprintf(f, "solution %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", s.p, s.q, s.m);
    }
    fclose(f);
    rename(tmp, S->state_file);
}

static void load_state(Scan *S)
{
    FILE *f = fopen(S->state_file, "r");
    if (!f) return;
    char line[512];
    if (!fgets(line, sizeof line, f) || strncmp(line, "A359636 cuda state v1", 21)) die("bad state file");
    int n = -1; u64 klo = 0, khi = 0, lsize = 0, fr = 0;
    unsigned long long surv = 0, trip = 0;
    while (fgets(line, sizeof line, f)) {
        u64 p, q, m;
        if (sscanf(line, "n %d", &n) == 1) continue;
        if (sscanf(line, "k_lo %" SCNu64, &klo) == 1) continue;
        if (sscanf(line, "k_hi %" SCNu64, &khi) == 1) continue;
        if (sscanf(line, "lsize %" SCNu64, &lsize) == 1) continue;
        if (sscanf(line, "frontier %" SCNu64, &fr) == 1) continue;
        if (sscanf(line, "survivors %llu", &surv) == 1) continue;
        if (sscanf(line, "triples %llu", &trip) == 1) continue;
        if (sscanf(line, "solution %" SCNu64 " %" SCNu64 " %" SCNu64, &p, &q, &m) == 3) {
            S->sol.push_back({p, q, m});
            S->nsol++;
            if (!S->best_p || p < S->best_p) { S->best_p = p; S->best_m = m; }
        }
    }
    fclose(f);
    if (n != S->n || klo != S->k_lo)
        die("state file %s belongs to a different scan (n=%d, k_lo=%" PRIu64 ")", S->state_file, n, klo);
    if (khi != S->k_hi && !g_quiet)              /* END changed: fine, the frontier is relative to k_lo */
        fprintf(stderr, "state has k_hi %" PRIu64 ", now %" PRIu64 " (END changed on resume)\n", khi, S->k_hi);
    if (lsize != S->lsize) {                 /* different launch size: restart at or below the old k frontier */
        u64 kfr = fr * lsize;
        fr = kfr / S->lsize;
        if (!g_quiet) fprintf(stderr, "state has launch size %" PRIu64 ", now %" PRIu64 ": frontier converted\n", lsize, S->lsize);
    }
    S->frontier = fr;
    S->resume_frontier = fr;
    S->survivors = surv;
    S->triples = trip;
    if (!g_quiet) {
        char c1[32];
        fmt_commas(c1, 6 * (klo + fr * lsize) + 3);
        fprintf(stderr, "resuming from %s: m = %s (%.2f%%), %zu solution(s) on file\n",
                S->state_file, c1, 100.0 * (double)(fr * lsize) / (double)(khi - klo), S->sol.size());
    }
}

static void progress(Scan *S, bool final_)
{
    static double last = 0;
    double t = now();
    if (!final_ && t - last < (stderr_tty ? 1.0 : 60.0)) return;
    last = t;
    u64 fr;
    { std::lock_guard<std::mutex> lk(S->fr_mu); fr = S->frontier; }
    u64 kdone = fr * S->lsize;
    if (kdone > S->k_hi - S->k_lo) kdone = S->k_hi - S->k_lo;
    u64 mnow = 6 * (S->k_lo + kdone) + 3;
    double frac = (double)kdone / (double)(S->k_hi - S->k_lo);
    double el = t - S->t0;
    double rate = el > 0 ? (double)(kdone - S->resume_frontier * S->lsize) / el : 0;
    char c1[32], eta[32], elap[32], best[48] = "";
    fmt_commas(c1, mnow);
    fmt_time(eta, rate > 0 ? (double)((S->k_hi - S->k_lo) - kdone) / rate : -1);
    fmt_time(elap, el);
    if (S->best_p) snprintf(best, sizeof best, " best=%" PRIu64, S->best_p);
    fprintf(stderr, "%sm=%s (%.2f%%)  %.1f Gm/s  surv=%llu trip=%llu sol=%llu%s  %s elapsed, ETA %s%s",
            stderr_tty ? "\r\033[K" : "", c1, frac * 100, rate * 6 / 1e9,
            (unsigned long long)S->survivors.load(), (unsigned long long)S->triples.load(),
            (unsigned long long)S->nsol.load(), best, elap, eta, stderr_tty && !final_ ? "" : "\n");
    fflush(stderr);
}

static u64 run_scan(int n, u64 start, u64 end, int cputhreads, u32 T, int blocks, const char *state_file,
                    int interval, bool scan_all, bool *certain)
{
    Scan *S = new Scan();
    S->n = n;
    S->scan_all = scan_all;
    S->state_file = state_file;

    u64 m_lo = start < 9 ? 9 : start;
    while (m_lo % 6 != 3) m_lo++;
    u64 m_hi = end + 2;
    if (m_hi < m_lo) die("empty range");
    S->k_lo = (m_lo - 3) / 6;
    S->k_hi = (m_hi - 3) / 6 + 1;
    u64 Bmax = m_hi + 1;

    u64 plim = (u64)iroot(Bmax, 3) + 1000;
    if (plim < 1u << 20) plim = 1u << 20;
    if (plim < 60000) plim = 60000;
    if (!g_primes || g_plimit < plim) build_primes(plim);

    if (T == 0) {
        static const u32 cand[] = { 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000, 15000, 20000 };
        int want = n < 6 ? n : 6;
        T = 0;
        for (u32 t : cand) if (compute_c(n, t, Bmax) >= want) { T = t; break; }
        if (!T) { want = n < 5 ? n : 5; for (u32 t : cand) if (compute_c(n, t, Bmax) >= want) { T = t; break; } }
        if (!T) T = 1000;
    }
    S->T = T;
    S->c = compute_c(n, T, Bmax);
    S->thr = S->c - 1;
    u32 thradd = (S->thr > 0) ? 0x010101u * (u32)(128 - S->thr) : 0x808080u, thrmask = 0x808080u;

    /* progressions for primes 29..T, ascending */
    for (size_t i = 0; i < g_nprimes && g_primes[i] <= T; i++) {
        u32 p = g_primes[i];
        if (p < 29) continue;
        for (int t = -1; t <= 1; t++) S->prog.push_back({p, target_residue(p, t), 1u << (8 * (t + 1))});
    }
    S->nsmall = 0;
    for (auto &pr : S->prog) if (pr.p <= PSMALL) S->nsmall++;
    /* patterns 5..13 and 17..23 */
    std::vector<u32> patA(PAT_A, 0), patB(PAT_B, 0);
    static const u32 pa[4] = { 5, 7, 11, 13 }, pb[3] = { 17, 19, 23 };
    for (u32 p : pa) if (p <= T) for (int t = -1; t <= 1; t++) { u32 r = target_residue(p, t), inc = 1u << (8 * (t + 1)); for (u32 k = r; k < PAT_A; k += p) patA[k] += inc; }
    for (u32 p : pb) if (p <= T) for (int t = -1; t <= 1; t++) { u32 r = target_residue(p, t), inc = 1u << (8 * (t + 1)); for (u32 k = r; k < PAT_B; k += p) patB[k] += inc; }

    if (blocks <= 0) blocks = 65536;
    if ((u64)blocks * T >= (1ull << 32))          /* kernel computes k0 mod p in 32 bits: blockIdx * (SEGK % p) */
        die("blocks * T must stay below 2^32 (have %d * %u); lower -b or -T", blocks, T);
    S->lsize = (u64)blocks * SEGK;
    S->nlaunch = (S->k_hi - S->k_lo + S->lsize - 1) / S->lsize;
    S->done.assign(RING, 0);
    S->outcap = 1u << 24;

    if (state_file) load_state(S);

    /* device setup */
    CUDA_CHECK(cudaFuncSetAttribute(sieve_kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, SEGK * sizeof(u32)));
    CUDA_CHECK(cudaMalloc(&S->d_prog, S->prog.size() * sizeof(Prog)));
    CUDA_CHECK(cudaMalloc(&S->d_kb, S->prog.size() * sizeof(u32)));
    CUDA_CHECK(cudaMalloc(&S->d_patA, PAT_A * sizeof(u32)));
    CUDA_CHECK(cudaMalloc(&S->d_patB, PAT_B * sizeof(u32)));
    CUDA_CHECK(cudaMalloc(&S->d_out, (size_t)S->outcap * sizeof(u64)));
    CUDA_CHECK(cudaMalloc(&S->d_count, sizeof(unsigned int)));
    CUDA_CHECK(cudaMemcpy(S->d_prog, S->prog.data(), S->prog.size() * sizeof(Prog), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(S->d_patA, patA.data(), PAT_A * sizeof(u32), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(S->d_patB, patB.data(), PAT_B * sizeof(u32), cudaMemcpyHostToDevice));
    std::vector<u32> kb(S->prog.size());
    std::vector<u64> hostout(S->outcap);

    if (!g_quiet) {
        char c1[32], c2[32];
        fmt_commas(c1, m_lo); fmt_commas(c2, m_hi);
        double marks = 0;
        for (auto &pr : S->prog) marks += 1.0 / pr.p;
        int dev; cudaDeviceProp prop; CUDA_CHECK(cudaGetDevice(&dev)); CUDA_CHECK(cudaGetDeviceProperties(&prop, dev));
        fprintf(stderr, "A359636 level n=%d (CUDA on %s, %d SMs): m in [%s, %s] (m = 3 mod 6, %" PRIu64 " values), "
                "%d blocks x %d k per launch, %d CPU verifier threads\n", n, prop.name, prop.multiProcessorCount,
                c1, c2, S->k_hi - S->k_lo, blocks, SEGK, cputhreads);
        fprintf(stderr, "sieve primes 5..%u (%zu strided progressions, %d warp-level), each target needs >= %d prime "
                "factors <= %u (threshold %d after the forced 2/3), %.2f marks per k\n",
                S->T, S->prog.size(), S->nsmall, S->c, S->T, S->thr, marks);
    }

    S->t0 = now();
    std::vector<std::thread> vth;
    for (int i = 0; i < cputhreads; i++) vth.emplace_back(verifier, S);

    double last_save = now();
    u64 li = S->frontier;
    for (; li < S->nlaunch && !g_stop && !S->stop_now; li++) {
        /* bound how far the GPU may run ahead of verification */
        {
            std::unique_lock<std::mutex> lk(S->q_mu);
            S->q_cv.wait(lk, [&] { return S->queue.size() < 64 || g_stop; });
        }
        {
            std::lock_guard<std::mutex> lk(S->fr_mu);
            if (li >= S->frontier + RING - 2) { /* should not happen with queue bound */ }
        }
        u64 kbase = S->k_lo + li * S->lsize;
        u64 klimit = kbase + S->lsize;
        if (klimit > S->k_hi) klimit = S->k_hi;
        for (size_t i = 0; i < S->prog.size(); i++) kb[i] = (u32)(kbase % S->prog[i].p);
        CUDA_CHECK(cudaMemcpy(S->d_kb, kb.data(), kb.size() * sizeof(u32), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemset(S->d_count, 0, sizeof(unsigned int)));
        int nb = (int)((klimit - kbase + SEGK - 1) / SEGK);
        sieve_kernel<<<nb, BLOCKDIM, SEGK * sizeof(u32)>>>(kbase, klimit, S->d_prog, S->d_kb, (int)S->prog.size(),
                                                           S->nsmall, S->d_patA, S->d_patB, thradd, thrmask,
                                                           S->d_out, S->d_count, S->outcap);
        CUDA_CHECK(cudaGetLastError());
        unsigned int cnt;
        CUDA_CHECK(cudaMemcpy(&cnt, S->d_count, sizeof cnt, cudaMemcpyDeviceToHost));
        if (cnt > S->outcap) die("survivor buffer overflow (%u > %u): raise capacity or lower -b", cnt, S->outcap);
        Batch b; b.idx = li; b.ks.resize(cnt);
        if (cnt) CUDA_CHECK(cudaMemcpy(b.ks.data(), S->d_out, (size_t)cnt * sizeof(u64), cudaMemcpyDeviceToHost));
        S->survivors += cnt;
        {
            std::lock_guard<std::mutex> lk(S->q_mu);
            S->q_pending_k += cnt;
            S->queue.push_back(std::move(b));
        }
        S->q_cv.notify_all();
        if (!g_quiet) progress(S, false);
        if (state_file && now() - last_save >= interval) { save_state(S); last_save = now(); }
    }
    /* drain */
    {
        std::lock_guard<std::mutex> lk(S->q_mu);
        S->q_done = true;
    }
    S->q_cv.notify_all();
    for (auto &t : vth) t.join();
    if (state_file) save_state(S);
    if (!g_quiet) progress(S, true);

    u64 covered_k = S->frontier >= S->nlaunch ? S->k_hi : S->k_lo + S->frontier * S->lsize;
    u64 covered_m = covered_k ? 6 * covered_k - 3 : 0;
    u64 result = 0;
    *certain = false;
    if (S->best_p && S->best_m <= covered_m) { result = S->best_p; *certain = true; }
    if (!g_quiet) {
        char c1[32];
        fmt_commas(c1, covered_m);
        fprintf(stderr, "scanned m <= %s: %llu survivors, %llu triples, %llu solutions\n", c1,
                (unsigned long long)S->survivors.load(), (unsigned long long)S->triples.load(),
                (unsigned long long)S->nsol.load());
    }
    bool partial = m_lo > 9;                       /* START given: statements are about this range only */
    if (*certain && !partial) printf("a(%d) = %" PRIu64 "\n", n, result);
    else if (*certain) printf("smallest qualifying p with m in [%" PRIu64 ", %" PRIu64 "] is %" PRIu64 " (partial range: a(%d) <= this)\n", m_lo, covered_m, result, n);
    else if (S->best_p) printf("a(%d) <= %" PRIu64 " (scan incomplete below it)\n", n, S->best_p);
    else if (S->frontier >= S->nlaunch && !partial) printf("a(%d) > %" PRIu64 " (no qualifying gap with p <= %" PRIu64 ")\n", n, end, end);
    else if (S->frontier >= S->nlaunch) printf("no qualifying gap with m in [%" PRIu64 ", %" PRIu64 "] (partial range)\n", m_lo, m_hi);
    else if (!partial) printf("a(%d) > %" PRIu64 " (scan stopped early; no qualifying gap with p below that)\n", n, covered_m > 2 ? covered_m - 2 : 0);
    else printf("no qualifying gap with m in [%" PRIu64 ", %" PRIu64 "] (partial range, stopped early)\n", m_lo, covered_m);
    fflush(stdout);
    cudaFree(S->d_prog); cudaFree(S->d_kb); cudaFree(S->d_patA); cudaFree(S->d_patB); cudaFree(S->d_out); cudaFree(S->d_count);
    delete S;
    return result;
}

/* ------------------------------------------------------------------ */
/* Commands                                                            */
/* ------------------------------------------------------------------ */

static int cmd_selftest(int cputhreads, bool with8)
{
    int fails = 0, top = with8 ? 8 : 7;
    for (int n = 1; n <= top; n++) {
        bool certain;
        double t = now();
        u64 r = run_scan(n, 0, KNOWN[n], cputhreads, 0, n <= 6 ? 256 : 4096, nullptr, 60, false, &certain);
        bool ok = certain && r == KNOWN[n];
        printf("selftest n=%d: got %" PRIu64 "%s expected %" PRIu64 " -> %s  (%.1fs)\n", n, r,
               certain ? "" : " (uncertain)", KNOWN[n], ok ? "ok" : "FAIL", now() - t);
        fflush(stdout);
        if (!ok) fails++;
    }
    printf("selftest %s\n", fails ? "FAILED" : "passed");
    return fails ? 1 : 0;
}

static int cmd_verify(int n, u64 p)
{
    build_primes(1u << 21);
    char buf[256];
    if (!is_prime(p)) { printf("%" PRIu64 " is not prime\n", p); return 1; }
    u64 q = next_prime(p);
    printf("p = %" PRIu64 ", next prime q = %" PRIu64 ", gap %" PRIu64 "\n", p, q, q - p);
    bool ok = q - p >= 4;
    for (u64 x = p + 1; x < q; x++) {
        fmt_factorization(buf, x);
        int w = omega_exact(x);
        printf("  %" PRIu64 " = %s  omega=%d%s\n", x, buf, w, w >= n ? "" : "  <-- fails");
        if (w < n) ok = false;
    }
    printf("gap after %" PRIu64 " %s for level %d\n", p, ok ? "QUALIFIES" : "does not qualify", n);
    return ok ? 0 : 1;
}

static void usage(void)
{
    fputs("usage:\n"
          "  a359636_cuda scan N [START] END [-t CPUTHREADS] [-T SIEVEMAX] [-b BLOCKS] [-S STATE] [-i SECS] [-a] [-q]\n"
          "  a359636_cuda verify N P\n"
          "  a359636_cuda selftest [-t CPUTHREADS] [-8]\n", stderr);
    exit(2);
}

int main(int argc, char **argv)
{
    stderr_tty = isatty(2);
    if (argc < 2) usage();
    signal(SIGINT, on_sigint);
    signal(SIGTERM, on_sigint);
    const char *cmd = argv[1];
    if (!strcmp(cmd, "selftest")) {
        int th = 8; bool with8 = false;
        for (int i = 2; i < argc; i++) {
            if (!strcmp(argv[i], "-t") && i + 1 < argc) th = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-8")) with8 = true;
            else if (!strcmp(argv[i], "-q")) g_quiet = true;
            else usage();
        }
        return cmd_selftest(th, with8);
    }
    if (!strcmp(cmd, "verify")) {
        if (argc != 4) usage();
        return cmd_verify(atoi(argv[2]), parse_num(argv[3]));
    }
    if (!strcmp(cmd, "scan")) {
        if (argc < 4) usage();
        int n = atoi(argv[2]);
        if (n < 1 || n > 20) die("N must be 1..20");
        u64 pos[2]; int npos = 0;
        int th = 8, blocks = 0, interval = 60; u32 T = 0; const char *state = nullptr; bool all = false;
        for (int i = 3; i < argc; i++) {
            if (!strcmp(argv[i], "-t") && i + 1 < argc) th = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-T") && i + 1 < argc) T = (u32)parse_num(argv[++i]);
            else if (!strcmp(argv[i], "-b") && i + 1 < argc) blocks = (int)parse_num(argv[++i]);
            else if (!strcmp(argv[i], "-S") && i + 1 < argc) state = argv[++i];
            else if (!strcmp(argv[i], "-i") && i + 1 < argc) interval = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-a")) all = true;
            else if (!strcmp(argv[i], "-q")) g_quiet = true;
            else if (argv[i][0] == '-' && !isdigit((unsigned char)argv[i][1])) usage();
            else if (npos < 2) pos[npos++] = parse_num(argv[i]);
            else usage();
        }
        if (npos == 0) usage();
        u64 start = npos == 2 ? pos[0] : 0, end = pos[npos - 1];
        if (th < 1) th = 1;
        if (T && T < 30) T = 30;
        bool certain;
        run_scan(n, start, end, th, T, blocks, state, interval, all, &certain);
        return 0;
    }
    usage();
    return 2;
}
