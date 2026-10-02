/*
 * a350826_cuda.cu
 *
 * CUDA version of ../a350826.c: counts prime sextuplets (p, p+4, p+6, p+10, p+12,
 * p+16) by initial member p in [LO, HI), for OEIS A350826 / A063501.  It uses the
 * same wheel classes, chunks, per-bin counts and checksums, chunk-log lines and
 * checkpoint format as the CPU tool, so runs can be compared chunk by chunk (for the
 * same w, B and chunk) and a checkpoint written by one tool can be resumed by the
 * other.
 *
 * GPU layout
 * ----------
 * One thread block (NT = 512 threads) per residue class r mod M.  The class
 * progression p = p0 + M k is sieved in segments of SEGB bits held in shared memory:
 *
 *   - primes q < PATMAX (default 256): 64-bit word patterns as in the CPU tool
 *     (T_q[s] = bits i with (s + i) mod q in {0, k1..k5}); thread t owns the WPT words
 *     [t WPT, (t+1) WPT) and ORs the patterns into registers;
 *   - primes PATMAX <= q < QW (default 4096, "medium"): one warp per prime; lane l
 *     marks rounds l-1, l+31, ... (a round = one period of q = six bits) with
 *     shared-memory atomicOr;
 *   - primes QW <= q <= B ("large"): one thread per prime, phases in registers, the
 *     CPU tool's round loop with atomicOr.
 *
 * Survivors go to a shared list and are tested breadth first: all of them get the
 * base-2 strong test of p, the passing ones that of p+4, and so on, so warps stay
 * converged.  When all six members pass, each gets a strong Lucas test (BPSW) on
 * the GPU; failures are reported as SPSP lines and not counted.  Initial members
 * below 2^20 are tested on the host.  Arithmetic: Montgomery with R = 2^64 below
 * 2^64 and R = 2^128 (unsigned __int128) above.
 *
 * Usage
 * -----
 *   a350826_cuda count LO HI [-b B1,B2,..] [-w W] [-B B] [-c CH] [-C C0:C1 | -P I/N]
 *                            [-S FILE] [-i SECS] [-L FILE] [-q] [-d DEV]
 *   a350826_cuda bench LO HI [-n N] [-w W] [-B B] [-d DEV]
 *   a350826_cuda selftest [N] [-d DEV]
 * Same meaning as for a350826.c; -d selects the CUDA device.
 *
 * Build:  nvcc -O3 -std=c++17 -arch=sm_121 -o a350826_cuda a350826_cuda.cu
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
#include <climits>
#include <cinttypes>
#include <ctime>
#include <unistd.h>
#include <cuda_runtime.h>

typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t i64;
typedef unsigned __int128 u128;

#define GOLD    0x9E3779B97F4A7C15ULL
#define MAXBIN  32
#define MAXSPSP 1024
#define P0      ((u64)1 << 20)
#define BMAX    ((u32)1 << 20)
#define HIMAX   ((u128)1 << 80)

/* kernel geometry (override with -D at build time) */
#ifndef NT
#define NT      512                 /* threads per block */
#endif
#ifndef MINB
#define MINB    2                   /* blocks per SM (launch bounds) */
#endif
#ifndef SEGB
#define SEGB    (1u << 18)          /* bits per segment (dynamic shared memory) */
#endif
#define NWARP   (NT / 32)
#define SEGW    (SEGB / 64)         /* u64 words per segment */
#define WPT     (SEGW / NT)         /* words per thread in the pattern pass */
#if SEGW % NT
#error "SEGB/64 must be a multiple of NT"
#endif
#define MAXPAT  64                  /* pattern primes */
#define MAXMPL  4                   /* medium primes per lane (NWARP * 32 * MAXMPL total) */
#define MAXLPT  48                  /* large primes per thread (NT * MAXLPT total) */

static const u32 OFF[6] = {0, 4, 6, 10, 12, 16};

static bool stderr_tty;
static volatile sig_atomic_t g_stop = 0;

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) die("CUDA error %s at %s:%d", cudaGetErrorString(e_), __FILE__, __LINE__); } while (0)

/* ------------------------------------------------------------------ */
/* Host helpers (as in a350826.c)                                      */
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
    if (s < 0 || s != s || s > 1e12) { snprintf(buf, n, "?"); return buf; }
    long t = (long)s;
    if (t >= 86400) snprintf(buf, n, "%ldd%02ldh%02ldm", t / 86400, (t / 3600) % 24, (t / 60) % 60);
    else if (t >= 3600) snprintf(buf, n, "%ldh%02ldm%02lds", t / 3600, (t / 60) % 60, t % 60);
    else if (t >= 60) snprintf(buf, n, "%ldm%02lds", t / 60, t % 60);
    else snprintf(buf, n, "%.1fs", s);
    return buf;
}

static const char *u128s(u128 x)
{
    static char ring[16][44];
    static int k = 0;
    char *b = ring[k++ & 15], t[44];
    int n = 0;
    do { t[n++] = (char)('0' + (int)(x % 10)); x /= 10; } while (x);
    for (int i = 0; i < n; i++) b[i] = t[n - 1 - i];
    b[n] = 0;
    return b;
}

static bool mul_ok(u128 a, u128 b, u128 *r)
{
    if (a && b > ~(u128)0 / a) return false;
    *r = a * b;
    return true;
}

static bool all_digits(const char *s, size_t n)
{
    if (n == 0) return false;
    for (size_t i = 0; i < n; i++)
        if (!isdigit((unsigned char)s[i])) return false;
    return true;
}

static bool parse_term(const char *s, size_t n, u128 *out)
{
    char buf[80], *c;
    if (n == 0 || n >= sizeof buf) return false;
    memcpy(buf, s, n);
    buf[n] = 0;
    if ((c = strchr(buf, '^'))) {
        *c = 0;
        if (!all_digits(buf, strlen(buf)) || !all_digits(c + 1, strlen(c + 1))) return false;
        u128 b = strtoull(buf, 0, 10), r = 1;
        u64 e = strtoull(c + 1, 0, 10);
        while (e--) if (!mul_ok(r, b, &r)) return false;
        *out = r;
        return true;
    }
    if ((c = strpbrk(buf, "eE"))) {
        *c = 0;
        char *dot = strchr(buf, '.');
        size_t ip = dot ? (size_t)(dot - buf) : strlen(buf);
        const char *frac = dot ? dot + 1 : "";
        if ((ip && !all_digits(buf, ip)) || (*frac && !all_digits(frac, strlen(frac)))) return false;
        if (!all_digits(c + 1, strlen(c + 1))) return false;
        long e = strtol(c + 1, 0, 10) - (long)strlen(frac);
        u128 m = 0;
        for (size_t i = 0; i < ip; i++) { if (!mul_ok(m, 10, &m)) return false; m += (u128)(buf[i] - '0'); }
        for (const char *f = frac; *f; f++) { if (!mul_ok(m, 10, &m)) return false; m += (u128)(*f - '0'); }
        for (; e > 0; e--) if (!mul_ok(m, 10, &m)) return false;
        for (; e < 0; e++) { if (m % 10) return false; m /= 10; }
        *out = m;
        return true;
    }
    if (!all_digits(buf, n)) return false;
    u128 m = 0;
    for (size_t i = 0; i < n; i++) { if (!mul_ok(m, 10, &m)) return false; m += (u128)(buf[i] - '0'); }
    *out = m;
    return true;
}

static u128 parse_num(const char *s)
{
    size_t n = strlen(s);
    for (size_t i = n; i-- > 1; ) {
        if (s[i] != '+' && s[i] != '-') continue;
        if (s[i - 1] == 'e' || s[i - 1] == 'E') continue;
        u128 a, b;
        if (!parse_term(s, i, &a) || !parse_term(s + i + 1, n - i - 1, &b)) die("bad number '%s'", s);
        if (s[i] == '+') {
            if (a + b < a) die("overflow in '%s'", s);
            return a + b;
        }
        if (b > a) die("negative result in '%s'", s);
        return a - b;
    }
    u128 v;
    if (!parse_term(s, n, &v)) die("bad number '%s'", s);
    return v;
}

static u64 parse_u64(const char *s, u64 lo, u64 hi, const char *what)
{
    u128 v = parse_num(s);
    if (v < lo || v > hi) die("%s must be between %" PRIu64 " and %" PRIu64, what, lo, hi);
    return (u64)v;
}

static inline u64 mix64(u64 z)
{
    z += GOLD;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline u64 hashp(u128 p) { return mix64((u64)p ^ ((u64)(p >> 64) * GOLD)); }

/* reference primality (host) */
static u128 mulmod_ref(u128 a, u128 b, u128 n)
{
    if (!(n >> 64)) return ((a % n) * (b % n)) % n;
    u128 r = 0;
    a %= n;
    while (b) {
        if (b & 1) { r += a; if (r >= n) r -= n; }
        a += a; if (a >= n) a -= n;
        b >>= 1;
    }
    return r;
}

static u128 powmod_ref(u128 a, u128 e, u128 n)
{
    u128 r = 1 % n;
    a %= n;
    while (e) {
        if (e & 1) r = mulmod_ref(r, a, n);
        a = mulmod_ref(a, a, n);
        e >>= 1;
    }
    return r;
}

static bool mr_ref(u128 n, u64 a)
{
    u128 d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    u128 x = powmod_ref(a, d, n);
    if (x == 1 || x == n - 1) return true;
    for (int i = 1; i < s; i++) {
        x = mulmod_ref(x, x, n);
        if (x == n - 1) return true;
        if (x == 1) return false;
    }
    return false;
}

static bool is_prime_ref(u128 n)
{
    static const u32 sp[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41};
    if (n < 2) return false;
    for (int i = 0; i < 13; i++) {
        if (n == sp[i]) return true;
        if (n % sp[i] == 0) return false;
    }
    if (n < 43 * 43) return true;
    if (n >= HIMAX * 2) die("is_prime_ref: argument too large");
    for (int i = 0; i < 13; i++) if (!mr_ref(n, sp[i])) return false;
    return true;
}

static bool is_sext_ref(u128 p)
{
    for (int i = 0; i < 6; i++) if (!is_prime_ref(p + OFF[i])) return false;
    return true;
}

/* ------------------------------------------------------------------ */
/* Wheel, sieving primes, patterns (host; as in a350826.c)             */
/* ------------------------------------------------------------------ */

#define MAXWD 12
static const u32 SMALLP[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

static struct {
    u32 w;
    u64 M, C;
    int nd;
    u32 q[MAXWD], na[MAXWD];
    u8 a[MAXWD][48];
    u64 E[MAXWD];
    u64 rad[MAXWD];
    u64 r0;
} W;

static u64 invmod(u64 a, u64 m)
{
    i64 t = 0, nt = 1;
    i64 r = (i64)m, nr = (i64)(a % m);
    while (nr) {
        i64 qq = r / nr, x;
        x = t - qq * nt; t = nt; nt = x;
        x = r - qq * nr; r = nr; nr = x;
    }
    if (r != 1) die("invmod: not invertible");
    return (u64)(t < 0 ? t + (i64)m : t);
}

static bool wheel_ok(u32 w)
{
    if (w < 7 || w > 47) return false;
    for (int i = 0; i < 15; i++) if (SMALLP[i] == w) return true;
    return false;
}

static void wheel_setup(u32 w)
{
    if (!wheel_ok(w)) die("wheel prime must be one of 7, 11, 13, ..., 47");
    memset(&W, 0, sizeof W);
    W.w = w;
    W.M = 1;
    for (int i = 0; i < 15 && SMALLP[i] <= w; i++) W.M *= SMALLP[i];
    W.C = 1;
    for (int i = 4; i < 15 && SMALLP[i] <= w; i++) {
        const u32 q = SMALLP[i];
        const int j = W.nd++;
        W.q[j] = q;
        for (u32 a = 0; a < q; a++) {
            bool ok = true;
            for (int k = 0; k < 6; k++) if ((a + OFF[k]) % q == 0) ok = false;
            if (ok) W.a[j][W.na[j]++] = (u8)a;
        }
        W.C *= W.na[j];
        const u64 Mq = W.M / q;
        W.E[j] = (u64)((u128)Mq * invmod(Mq % q, q) % W.M);
    }
    const u64 M210 = W.M / 210;
    const u64 E210 = (u64)((u128)M210 * invmod(M210 % 210, 210) % W.M);
    W.r0 = (u64)((u128)97 * E210 % W.M);
    u64 rad = 1;
    for (int j = W.nd - 1; j >= 0; j--) { W.rad[j] = rad; rad *= W.na[j]; }
}

static inline u64 class_r(u64 c)
{
    u128 r = W.r0;
    for (int j = 0; j < W.nd; j++) r += (u128)W.a[j][(c / W.rad[j]) % W.na[j]] * W.E[j];
    return (u64)(r % W.M);
}

static u32 *small_primes(u32 n, int *count)
{
    u8 *s = (u8 *)calloc(n + 1, 1);
    if (!s) die("out of memory");
    int c = 0;
    for (u32 i = 2; i <= n; i++) {
        if (s[i]) continue;
        c++;
        for (u64 j = (u64)i * i; j <= n; j += i) s[j] = 1;
    }
    u32 *p = (u32 *)malloc(sizeof *p * (size_t)(c ? c : 1));
    if (!p) die("out of memory");
    c = 0;
    for (u32 i = 2; i <= n; i++) if (!s[i]) p[c++] = i;
    free(s);
    *count = c;
    return p;
}

static int cmp_u32(const void *a, const void *b)
{
    u32 x = *(const u32 *)a, y = *(const u32 *)b;
    return x < y ? -1 : x > y;
}

/* host copies of the device tables */
static struct {
    u32 B, patmax, qw;
    int np, npat, nmed, nlarge;
    u32 *q, *k, *minv, *bmod, *segmod, *pstep, *pstept, *poff;   /* k: 5 * np, SoA */
    u64 *tab;
    size_t ntab;
} SP;

static void tables_setup(u32 B, u32 patmax, u32 qw, u128 base)
{
    if (B <= W.w || B > BMAX) die("sieving bound must be in (w, 2^20]");
    int n;
    u32 *pr = small_primes(B, &n);
    free(SP.q); free(SP.k); free(SP.minv); free(SP.bmod); free(SP.segmod); free(SP.pstep); free(SP.pstept);
    free(SP.poff); free(SP.tab);
    SP.B = B; SP.patmax = patmax; SP.qw = qw;
    SP.q = (u32 *)malloc(sizeof(u32) * n);
    SP.k = (u32 *)malloc(sizeof(u32) * 5 * n);
    SP.minv = (u32 *)malloc(sizeof(u32) * n);
    SP.bmod = (u32 *)malloc(sizeof(u32) * n);
    SP.segmod = (u32 *)malloc(sizeof(u32) * n);
    SP.pstep = (u32 *)malloc(sizeof(u32) * MAXPAT);
    SP.pstept = (u32 *)malloc(sizeof(u32) * MAXPAT);
    SP.poff = (u32 *)malloc(sizeof(u32) * MAXPAT);
    if (!SP.q || !SP.k || !SP.minv || !SP.bmod || !SP.segmod || !SP.pstep || !SP.pstept || !SP.poff) die("out of memory");
    int np = 0;
    for (int i = 0; i < n; i++) if (pr[i] > W.w) np++;
    SP.np = np;
    u32 *kt = (u32 *)malloc(sizeof(u32) * 5 * np);
    if (!kt) die("out of memory");
    np = 0;
    for (int i = 0; i < n; i++) {
        const u32 q = pr[i];
        if (q <= W.w) continue;
        const u32 mi = (u32)invmod(W.M % q, q);
        u32 c[5];
        for (int k = 1; k < 6; k++) c[k - 1] = (u32)((q - (u64)OFF[k] * mi % q) % q);
        qsort(c, 5, sizeof *c, cmp_u32);
        SP.q[np] = q;
        for (int k = 0; k < 5; k++) kt[k * SP.np + np] = c[k];
        SP.minv[np] = mi;
        SP.bmod[np] = (u32)(base % q);
        SP.segmod[np] = SEGB % q;
        np++;
    }
    memcpy(SP.k, kt, sizeof(u32) * 5 * SP.np);
    free(kt);
    free(pr);
    /* classes of primes */
    SP.npat = 0;
    while (SP.npat < SP.np && SP.q[SP.npat] < patmax && SP.npat < MAXPAT) SP.npat++;
    int m = SP.npat;
    while (m < SP.np && SP.q[m] < qw) m++;
    SP.nmed = m - SP.npat;
    SP.nlarge = SP.np - m;
    if (SP.nmed > NWARP * 32 * MAXMPL) die("too many medium primes (raise MAXMPL or lower QW)");
    if (SP.nlarge > NT * MAXLPT) die("too many large primes for B = %u (max %d)", B, NT * MAXLPT);
    /* pattern tables */
    size_t words = 0;
    for (int i = 0; i < SP.npat; i++) words += SP.q[i];
    SP.tab = (u64 *)malloc(sizeof(u64) * (words ? words : 1));
    if (!SP.tab) die("out of memory");
    SP.ntab = words;
    size_t off = 0;
    for (int i = 0; i < SP.npat; i++) {
        const u32 q = SP.q[i];
        u8 bad[1024] = {0};
        bad[0] = 1;
        for (int j = 0; j < 5; j++) bad[SP.k[j * SP.np + i]] = 1;
        for (u32 s = 0; s < q; s++) {
            u64 w = 0;
            for (u32 b = 0; b < 64; b++) if (bad[(s + b) % q]) w |= 1ULL << b;
            SP.tab[off + s] = w;
        }
        SP.poff[i] = (u32)off;
        SP.pstep[i] = 64 % q;
        SP.pstept[i] = (u32)(((u64)64 * WPT) % q);
        off += q;
    }
}

/* ------------------------------------------------------------------ */
/* Device code                                                         */
/* ------------------------------------------------------------------ */

struct DevP {
    u64 M, bmM, Q0, rem0, base_lo, base_hi;
    int nbins, npat, nmed, nlarge, np;
    int skip;               /* bench ablation: 1 tests, 2 patterns, 4 medium, 8 large */
    u64 bnd_lo[MAXBIN + 1], bnd_hi[MAXBIN + 1];
};

struct Cnt {                    /* per sub-launch buffer counters */
    unsigned long long nb[7];   /* survivors in the buffer before round m (nb[0] from the sieve) */
};

struct Res {
    unsigned long long cnt[MAXBIN], cks[MAXBIN];
    unsigned long long surv, cand, spsp, bits;
    unsigned int ovf;           /* survivor buffer overflow: redo the launch in pieces */
    unsigned int nsp;
    unsigned long long sp_n[MAXSPSP][2], sp_p[MAXSPSP][2];
};

struct DevT {
    const u32 *q, *k, *minv, *bmod, *segmod, *pstep, *pstept, *poff;
    const u64 *tab;
};

__constant__ DevP P;
__constant__ u32 OFFD[6] = {0, 4, 6, 10, 12, 16};

__device__ __forceinline__ u64 d_inv64(u64 n)
{
    u64 x = n;
    #pragma unroll
    for (int i = 0; i < 5; i++) x *= 2 - n * x;
    return x;
}

__device__ __forceinline__ u64 d_mmul64(u64 a, u64 b, u64 n, u64 ni)
{
    const u64 tl = a * b, th = __umul64hi(a, b);
    const u64 m = tl * ni;
    const u64 mh = __umul64hi(m, n);
    const u64 r = th - mh;
    return th < mh ? r + n : r;
}

__device__ __forceinline__ u64 d_madd64(u64 a, u64 b, u64 n) { u64 s = a + b; return (s < a || s >= n) ? s - n : s; }
__device__ __forceinline__ u64 d_msub64(u64 a, u64 b, u64 n) { return a >= b ? a - b : a - b + n; }
__device__ __forceinline__ u64 d_mhalf64(u64 a, u64 n) { return (a & 1) ? (a >> 1) + (n >> 1) + 1 : a >> 1; }

/* 2^64 mod n for odd n > 1: with f = floor(2^64/n) < 256 estimated in float (error < 1),
   k = max(1, f_est - 1) <= f, r = 2^64 - k n, then subtract n until r < n */
__device__ __forceinline__ u64 d_r64(u64 n)
{
    if (n >> 56) {
        const u64 ke = (u64)(1.8446744e19f / (float)n);
        const u64 k = ke > 2 ? ke - 1 : 1;
        u64 r = 0 - k * n;
        while (r >= n) r -= n;
        return r;
    }
    return (u64)(0 - n) % n;
}

__device__ bool d_sprp2_64(u64 n)
{
    const u64 ni = d_inv64(n);
    const u64 one = d_r64(n), mone = n - one;
    u64 d = n - 1;
    const int s = __ffsll((long long)d) - 1;
    d >>= s;
    int b = 63 - __clzll((long long)d);
    u64 x = d_madd64(one, one, n);
    while (b-- > 0) {
        x = d_mmul64(x, x, n, ni);
        if ((d >> b) & 1) x = d_madd64(x, x, n);
    }
    if (x == one || x == mone) return true;
    for (int i = 1; i < s; i++) {
        x = d_mmul64(x, x, n, ni);
        if (x == mone) return true;
        if (x == one) return false;
    }
    return false;
}

__device__ int d_jacobi(i64 a, u128 n)
{
    int t = 1;
    u128 x;
    if (a < 0) {
        x = (u128)(-(a + 1)) + 1;
        if ((n & 3) == 3) t = -t;
    } else x = (u128)a;
    x %= n;
    while (x) {
        while (!(x & 1)) {
            x >>= 1;
            unsigned r = (unsigned)(n & 7);
            if (r == 3 || r == 5) t = -t;
        }
        u128 tmp = x; x = n; n = tmp;
        if ((x & 3) == 3 && (n & 3) == 3) t = -t;
        x %= n;
    }
    return n == 1 ? t : 0;
}

__device__ bool d_is_square(u128 n)
{
    u128 r = (u128)sqrt((double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r * r == n;
}

__device__ int d_selfridge(u128 n, i64 *Dout)
{
    if (d_is_square(n)) return 0;
    i64 D = 5;
    for (;;) {
        int j = d_jacobi(D, n);
        if (j == -1) break;
        if (j == 0 && (u128)(D < 0 ? -D : D) != n) return 0;
        D = D > 0 ? -(D + 2) : -(D - 2);
    }
    *Dout = D;
    return 1;
}

/* v * x mod n for small v >= 0 (double and add) */
__device__ u64 d_mulsmall64(u64 v, u64 x, u64 n)
{
    u64 r = 0;
    for (int b = 63 - __clzll((long long)(v | 1)); b >= 0; b--) {
        r = d_madd64(r, r, n);
        if ((v >> b) & 1) r = d_madd64(r, x, n);
    }
    return r;
}

__device__ bool d_slprp64(u64 n)
{
    i64 D;
    if (!d_selfridge(n, &D)) return false;
    const i64 Q = (1 - D) / 4;
    const u64 ni = d_inv64(n), one = d_r64(n);
    const u64 Qa = Q >= 0 ? d_mulsmall64((u64)Q, one, n) : n - d_mulsmall64((u64)(-Q), one, n);
    const u64 Qm = Qa == n ? 0 : Qa;
    const u64 Dm = D >= 0 ? d_mulsmall64((u64)D, one, n) : n - d_mulsmall64((u64)(-D), one, n);
    u64 d = n + 1;
    const int s = __ffsll((long long)d) - 1;
    d >>= s;
    u64 U = one, V = one, Qk = Qm;
    for (int b = 62 - __clzll((long long)d); b >= 0; b--) {
        U = d_mmul64(U, V, n, ni);
        V = d_msub64(d_mmul64(V, V, n, ni), d_madd64(Qk, Qk, n), n);
        Qk = d_mmul64(Qk, Qk, n, ni);
        if ((d >> b) & 1) {
            const u64 U2 = d_madd64(U, V, n);
            const u64 V2 = d_madd64(d_mmul64(Dm, U, n, ni), V, n);
            U = d_mhalf64(U2, n);
            V = d_mhalf64(V2, n);
            Qk = d_mmul64(Qk, Qm, n, ni);
        }
    }
    if (U == 0 || V == 0) return true;
    for (int r = 1; r < s; r++) {
        V = d_msub64(d_mmul64(V, V, n, ni), d_madd64(Qk, Qk, n), n);
        if (V == 0) return true;
        Qk = d_mmul64(Qk, Qk, n, ni);
    }
    return false;
}

/* 128-bit Montgomery (n < 2^127), as mmul128 in a350826.c */
__device__ __forceinline__ u128 d_mmul128(u128 a, u128 b, u128 n, u64 ni)
{
    const u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64);
    const u64 n0 = (u64)n, n1 = (u64)(n >> 64);
    const u128 p00 = (u128)a0 * b0, p01 = (u128)a0 * b1, p10 = (u128)a1 * b0, p11 = (u128)a1 * b1;
    const u64 t0 = (u64)p00;
    const u128 mid = (p00 >> 64) + (u64)p01 + (u64)p10;
    const u64 t1 = (u64)mid;
    const u128 hi = (mid >> 64) + (p01 >> 64) + (p10 >> 64) + p11;
    const u64 m = t0 * ni;
    const u128 q0 = (u128)m * n0, q1 = (u128)m * n1;
    const u128 s1 = (u128)t1 + (q0 >> 64) + (u64)q1 + (t0 != 0);
    const u64 u1 = (u64)s1;
    const u128 hi2 = hi + (s1 >> 64) + (q1 >> 64);
    const u64 m2 = u1 * ni;
    const u128 r0 = (u128)m2 * n0, r1 = (u128)m2 * n1;
    const u128 s2 = (u128)(u64)hi2 + (r0 >> 64) + (u64)r1 + (u1 != 0);
    const u128 top = (hi2 >> 64) + (s2 >> 64) + (r1 >> 64);
    const u128 res = (top << 64) | (u64)s2;
    return res >= n ? res - n : res;
}

__device__ __forceinline__ u128 d_madd128(u128 a, u128 b, u128 n) { u128 s = a + b; return s >= n ? s - n : s; }
__device__ __forceinline__ u128 d_msub128(u128 a, u128 b, u128 n) { return a >= b ? a - b : a - b + n; }
__device__ __forceinline__ u128 d_mhalf128(u128 a, u128 n) { return (a & 1) ? (a >> 1) + (n >> 1) + 1 : a >> 1; }

/* 2^128 mod n for 2^64 < n < 2^100.  The quotient Q = floor(2^128/n) < 2^64 is estimated in
   double (error < 2^13) and lowered by 2^21, so q <= Q and r = 2^128 - q n lies in
   [0, (2^21 + 2^13 + 1) n); one more double estimate brings r below 4n. */
__device__ u128 d_r128(u128 n)
{
    const double inv = 1.0 / (double)n;
    const double qd = 3.402823669209385e38 * inv;
    u64 q = qd >= 1.8446744073709552e19 ? ~0ULL : (u64)qd;
    q = q > (1ULL << 21) ? q - (1ULL << 21) : 1;
    u128 r = (u128)0 - (u128)q * n;
    u64 k = (u64)((double)r * inv);
    k = k > 2 ? k - 2 : 0;
    r -= (u128)k * n;
    while (r >= n) r -= n;
    return r;
}

__device__ bool d_sprp2_128(u128 n)
{
    const u64 ni = 0 - d_inv64((u64)n);
    const u128 one = d_r128(n), mone = n - one;
    u128 d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    int b = 127;
    while (!((d >> b) & 1)) b--;
    u128 x = d_madd128(one, one, n);
    while (b-- > 0) {
        x = d_mmul128(x, x, n, ni);
        if ((d >> b) & 1) x = d_madd128(x, x, n);
    }
    if (x == one || x == mone) return true;
    for (int i = 1; i < s; i++) {
        x = d_mmul128(x, x, n, ni);
        if (x == mone) return true;
        if (x == one) return false;
    }
    return false;
}

__device__ u128 d_mulsmall128(u64 v, u128 x, u128 n)
{
    u128 r = 0;
    for (int b = 63 - __clzll((long long)(v | 1)); b >= 0; b--) {
        r = d_madd128(r, r, n);
        if ((v >> b) & 1) r = d_madd128(r, x, n);
    }
    return r;
}

__device__ bool d_slprp128(u128 n)
{
    i64 D;
    if (!d_selfridge(n, &D)) return false;
    const i64 Q = (1 - D) / 4;
    const u64 ni = 0 - d_inv64((u64)n);
    const u128 one = d_r128(n);
    const u128 Qa = Q >= 0 ? d_mulsmall128((u64)Q, one, n) : n - d_mulsmall128((u64)(-Q), one, n);
    const u128 Qm = Qa == n ? 0 : Qa;
    const u128 Dm = D >= 0 ? d_mulsmall128((u64)D, one, n) : n - d_mulsmall128((u64)(-D), one, n);
    u128 d = n + 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    int top = 127;
    while (!((d >> top) & 1)) top--;
    u128 U = one, V = one, Qk = Qm;
    for (int b = top - 1; b >= 0; b--) {
        U = d_mmul128(U, V, n, ni);
        V = d_msub128(d_mmul128(V, V, n, ni), d_madd128(Qk, Qk, n), n);
        Qk = d_mmul128(Qk, Qk, n, ni);
        if ((d >> b) & 1) {
            const u128 U2 = d_madd128(U, V, n);
            const u128 V2 = d_madd128(d_mmul128(Dm, U, n, ni), V, n);
            U = d_mhalf128(U2, n);
            V = d_mhalf128(V2, n);
            Qk = d_mmul128(Qk, Qm, n, ni);
        }
    }
    if (U == 0 || V == 0) return true;
    for (int r = 1; r < s; r++) {
        V = d_msub128(d_mmul128(V, V, n, ni), d_madd128(Qk, Qk, n), n);
        if (V == 0) return true;
        Qk = d_mmul128(Qk, Qk, n, ni);
    }
    return false;
}

/* 96-bit Montgomery (three 32-bit limbs, R = 2^96, CIOS) for 2^64 <= n < 2^96; ni = -n^-1 mod 2^32 */
__device__ __forceinline__ u128 d_mmul96(u128 a, u128 b, u128 n, u32 ni)
{
    const u32 a0 = (u32)a, a1 = (u32)(a >> 32), a2 = (u32)(a >> 64);
    const u32 b0 = (u32)b, b1 = (u32)(b >> 32), b2 = (u32)(b >> 64);
    const u32 n0 = (u32)n, n1 = (u32)(n >> 32), n2 = (u32)(n >> 64);
    u32 t0 = 0, t1 = 0, t2 = 0, t3 = 0;
    u64 x;
    u32 c, m;
#define CIOS_STEP(ai)                                                                     \
    x = (u64)t0 + (u64)(ai) * b0;           t0 = (u32)x; c = (u32)(x >> 32);              \
    x = (u64)t1 + (u64)(ai) * b1 + c;       t1 = (u32)x; c = (u32)(x >> 32);              \
    x = (u64)t2 + (u64)(ai) * b2 + c;       t2 = (u32)x; c = (u32)(x >> 32);              \
    x = (u64)t3 + c;                        t3 = (u32)x; { const u32 t4 = (u32)(x >> 32); \
    m = t0 * ni;                                                                          \
    x = (u64)t0 + (u64)m * n0;              c = (u32)(x >> 32);                           \
    x = (u64)t1 + (u64)m * n1 + c;          t0 = (u32)x; c = (u32)(x >> 32);              \
    x = (u64)t2 + (u64)m * n2 + c;          t1 = (u32)x; c = (u32)(x >> 32);              \
    x = (u64)t3 + c;                        t2 = (u32)x; t3 = t4 + (u32)(x >> 32); }
    CIOS_STEP(a0)
    CIOS_STEP(a1)
    CIOS_STEP(a2)
#undef CIOS_STEP
    const u128 r = ((u128)t3 << 96) | ((u128)t2 << 64) | ((u128)t1 << 32) | t0;
    return r >= n ? r - n : r;
}

/* 2^96 mod n for 2^64 <= n < 2^96.  Q = floor(2^96/n) <= 2^32 is estimated in float (error
   < 2^9) and lowered by 2^10, so r = 2^96 - q n is in [0, 2^11 n); a second float estimate of
   r/n (error < 1) brings r below 3n. */
__device__ __forceinline__ u128 d_r96(u128 n)
{
    const float nf = (float)n;
    u64 q = (u64)__fdividef(7.9228163e28f, nf);
    q = q > 1024 ? q - 1024 : 0;
    u128 r = ((u128)1 << 96) - (u128)q * n;
    u64 k = (u64)__fdividef((float)r, nf);
    k = k > 1 ? k - 1 : 0;
    r -= (u128)k * n;
    while (r >= n) r -= n;
    return r;
}

__device__ bool d_sprp2_96(u128 n)
{
    const u32 ni = (u32)(0 - d_inv64((u64)n));
    const u128 one = d_r96(n), mone = n - one;
    u128 d = n - 1;
    int s;
    if ((u64)d) s = __ffsll((long long)(u64)d) - 1;
    else s = 63 + __ffsll((long long)(u64)(d >> 64));
    d >>= s;
    int b = (d >> 64) ? 127 - __clzll((long long)(u64)(d >> 64)) : 63 - __clzll((long long)(u64)d);
    u128 x = d_madd128(one, one, n);
    while (b-- > 0) {
        x = d_mmul96(x, x, n, ni);
        if ((d >> b) & 1) x = d_madd128(x, x, n);
    }
    if (x == one || x == mone) return true;
    for (int i = 1; i < s; i++) {
        x = d_mmul96(x, x, n, ni);
        if (x == mone) return true;
        if (x == one) return false;
    }
    return false;
}

__device__ __forceinline__ bool d_sprp2(u128 n)
{
    if (!(n >> 64)) return d_sprp2_64((u64)n);
    return (n >> 96) ? d_sprp2_128(n) : d_sprp2_96(n);
}
__device__ __forceinline__ bool d_slprp(u128 n) { return (n >> 64) ? d_slprp128(n) : d_slprp64((u64)n); }

__device__ __forceinline__ u64 d_mix64(u64 z)
{
    z += GOLD;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* all six members passed base 2: confirm with strong Lucas and record */
__device__ void d_candidate(u128 p, Res *res)
{
    atomicAdd(&res->cand, 1ULL);
    for (int m = 0; m < 6; m++) {
        const u128 n = p + OFFD[m];
        if (!d_slprp(n)) {
            atomicAdd(&res->spsp, 1ULL);
            const unsigned i = atomicAdd(&res->nsp, 1u);
            if (i < MAXSPSP) {
                res->sp_n[i][0] = (u64)n; res->sp_n[i][1] = (u64)(n >> 64);
                res->sp_p[i][0] = (u64)p; res->sp_p[i][1] = (u64)(p >> 64);
            }
            return;
        }
    }
    int b = 0;
    while (b + 1 < P.nbins && p >= (((u128)P.bnd_hi[b + 1] << 64) | P.bnd_lo[b + 1])) b++;
    atomicAdd(&res->cnt[b], 1ULL);
    atomicAdd(&res->cks[b], (unsigned long long)d_mix64((u64)p ^ ((u64)(p >> 64) * GOLD)));
}

/* phase of prime i at the first candidate p0 = base + delta: phi = -(p0 mod q) M^-1 mod q */
__device__ __forceinline__ u32 d_phase(const DevT &T, int i, u64 delta)
{
    const u32 q = T.q[i];
    u32 u = (u32)(delta % q) + T.bmod[i];
    if (u >= q) u -= q;
    return u ? (u32)((u64)(q - u) * T.minv[i] % q) : 0;
}

#ifdef NOATOMIC                             /* timing experiment only: racy, wrong results */
#define AOR(a, v) (*(a) |= (v))
#else
#define AOR(a, v) atomicOr((a), (v))
#endif
/* unchecked mark, and a branch-free checked one (out-of-range positions OR 0 into word 0) */
#define AMARK(x) AOR(bw + ((u32)(x) >> 5), 1u << ((u32)(x) & 31))
#define MARK(x) do { const u32 x_ = (u32)(x); const bool in_ = x_ < len; \
                     AOR(bw + (in_ ? x_ >> 5 : 0u), in_ ? 1u << (x_ & 31) : 0u); } while (0)

template <int LPT>
__global__ void __launch_bounds__(NT, MINB) sext_kernel(const u64 *__restrict__ cls_r, DevT T, Res *__restrict__ res,
                                                        Cnt *__restrict__ cnt, u64 *__restrict__ buf, u64 cap)
{
    extern __shared__ u64 bm[];             /* SEGW words */
    __shared__ u32 phs[MAXPAT];
    u32 *const bw = (u32 *)bm;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;

    const u64 M = P.M;
    const u64 r = cls_r[blockIdx.x];
    const u64 delta = r >= P.bmM ? r - P.bmM : r + (M - P.bmM);
    const u64 L = delta <= P.rem0 ? P.Q0 + 1 : P.Q0;
    if (L == 0) return;
    const u128 ps = (((u128)P.base_hi << 64) | P.base_lo) + delta;
    const bool fits64 = ps + (u128)M * (L - 1) + 16 < ((u128)1 << 64);

    /* initial phases */
    for (int i = tid; i < P.npat; i += NT) phs[i] = d_phase(T, i, delta);
    u32 mph[MAXMPL];
    #pragma unroll
    for (int i = 0; i < MAXMPL; i++) {
        const int m = warp + NWARP * (lane + 32 * i);
        mph[i] = m < P.nmed ? d_phase(T, P.npat + m, delta) : 0;
    }
    const int lb = P.npat + P.nmed;
    u32 lph[LPT];
    #pragma unroll
    for (int i = 0; i < LPT; i++) {
        const int g = tid + NT * i;
        lph[i] = g < P.nlarge ? d_phase(T, lb + g, delta) : 0;
    }
    unsigned long long nsurv = 0;
    __syncthreads();

    for (u64 k0 = 0; k0 < L; k0 += SEGB) {
        const u32 len = L - k0 < SEGB ? (u32)(L - k0) : SEGB;
        const bool last = k0 + SEGB >= L;

        /* 1. pattern primes, into registers */
        {
            u64 acc[WPT];
            #pragma unroll
            for (int w = 0; w < WPT; w++) acc[w] = 0;
            for (int i = 0; i < ((P.skip & 2) ? 0 : P.npat); i++) {
                const u32 q = T.q[i], st = T.pstep[i];
                const u64 *tab = T.tab + T.poff[i];
                const u32 f = phs[i];
                const u32 s0 = f ? q - f : 0;
                u32 s = (u32)((s0 + (u64)tid * T.pstept[i]) % q);
                #pragma unroll
                for (int w = 0; w < WPT; w++) {
                    acc[w] |= __ldg(tab + s);
                    s += st;
                    s = s >= q ? s - q : s;
                }
            }
            #pragma unroll
            for (int w = 0; w < WPT; w++) bm[tid * WPT + w] = acc[w];
        }
        __syncthreads();

        /* 2. medium primes: one warp per prime, lanes take rounds */
        #pragma unroll
        for (int i = 0; i < MAXMPL; i++) {
            for (int jl = 0; jl < 32; jl++) {
                const int m = warp + NWARP * (jl + 32 * i);
                if (m >= P.nmed || (P.skip & 4)) break;
                const int pi = P.npat + m;
                const u32 phi = __shfl_sync(0xffffffffu, mph[i], jl);
                const int q = (int)T.q[pi];
                const int k1 = (int)T.k[pi], k2 = (int)T.k[P.np + pi], k3 = (int)T.k[2 * P.np + pi];
                const int k4 = (int)T.k[3 * P.np + pi], k5 = (int)T.k[4 * P.np + pi];
                int b = (int)phi + (lane - 1) * q;
                if (b < 0) {                   /* lane 0, round -1 */
                    MARK(b + k1); MARK(b + k2); MARK(b + k3); MARK(b + k4); MARK(b + k5);
                    b += 32 * q;
                }
                for (const int lim = (int)len - k5; b < lim; b += 32 * q) {
                    AMARK(b); AMARK(b + k1); AMARK(b + k2); AMARK(b + k3); AMARK(b + k4); AMARK(b + k5);
                }
                if (b < (int)len) {
                    MARK(b); MARK(b + k1); MARK(b + k2); MARK(b + k3); MARK(b + k4);
                }
            }
            if (!last) {
                const int m = warp + NWARP * (lane + 32 * i);
                if (m < P.nmed) {
                    const int pi = P.npat + m;
                    const u32 q = T.q[pi], sm = T.segmod[pi];
                    mph[i] = mph[i] >= sm ? mph[i] - sm : mph[i] + q - sm;
                }
            }
        }

        /* 3. large primes: one thread per prime */
        #pragma unroll
        for (int i = 0; i < LPT; i++) {
            const int g = tid + NT * i;
            if (g < P.nlarge && !(P.skip & 8)) {
                const int pi = lb + g;
                const u32 q = T.q[pi];
                const u32 k1 = T.k[pi], k2 = T.k[P.np + pi], k3 = T.k[2 * P.np + pi];
                const u32 k4 = T.k[3 * P.np + pi], k5 = T.k[4 * P.np + pi];
                const u32 phi = lph[i];
                const u32 t = q - phi;
                if (k5 >= t) {
                    if (k1 >= t) MARK(k1 - t);
                    if (k2 >= t) MARK(k2 - t);
                    if (k3 >= t) MARK(k3 - t);
                    if (k4 >= t) MARK(k4 - t);
                    MARK(k5 - t);
                }
                u32 b = phi;
                if (len > k5) {
                    const u32 lim = len - k5;
                    for (; b < lim; b += q) {
                        AMARK(b); AMARK(b + k1); AMARK(b + k2); AMARK(b + k3); AMARK(b + k4); AMARK(b + k5);
                    }
                }
                if (b < len) {
                    MARK(b); MARK(b + k1); MARK(b + k2); MARK(b + k3); MARK(b + k4);
                }
                lph[i] = b >= len ? b - len : b + q - len;
            }
        }
        __syncthreads();

        /* 4. survivors of my words into the global buffer (one atomic per warp) */
        if (!(P.skip & 16)) {
            u64 x[WPT];
            u32 c = 0;
            #pragma unroll
            for (int w = 0; w < WPT; w++) {
                const u32 j = tid * WPT + w;
                u64 v = ~bm[j];
                if ((j + 1) * 64 > len) v = (j * 64 >= len) ? 0 : v & ((1ULL << (len - j * 64)) - 1);
                x[w] = v;
                c += __popcll(v);
            }
            u32 incl = c;
            #pragma unroll
            for (int o = 1; o < 32; o <<= 1) {
                const u32 y = __shfl_up_sync(0xffffffffu, incl, o);
                if (lane >= o) incl += y;
            }
            const u32 wtot = __shfl_sync(0xffffffffu, incl, 31);
            unsigned long long wbase = 0;
            if (lane == 31 && wtot) wbase = atomicAdd(&cnt->nb[0], (unsigned long long)wtot);
            wbase = __shfl_sync(0xffffffffu, wbase, 31);
            u64 pos = wbase + incl - c;
            nsurv += c;
            const u128 pk = ps + (u128)M * k0;
            #pragma unroll
            for (int w = 0; w < WPT; w++) {
                u64 v = x[w];
                const u32 j = tid * WPT + w;
                while (v) {
                    const u32 kk = j * 64 + (u32)(__ffsll((long long)v) - 1);
                    v &= v - 1;
                    if (pos < cap) {
                        const u128 p = pk + (u128)M * kk;
                        buf[2 * pos] = (u64)p;
                        buf[2 * pos + 1] = (u64)(p >> 64);
                    } else res->ovf = 1;
                    pos++;
                }
            }
        }

        /* pattern phases for the next segment */
        if (!last)
            for (int i = tid; i < P.npat; i += NT) {
                const u32 q = T.q[i], sm = T.segmod[i];
                phs[i] = phs[i] >= sm ? phs[i] - sm : phs[i] + q - sm;
            }
        __syncthreads();
    }
    /* block totals */
    for (int o = 16; o; o >>= 1) nsurv += __shfl_down_sync(0xffffffffu, nsurv, o);
    if (lane == 0) atomicAdd(&res->surv, nsurv);
    if (tid == 0) atomicAdd(&res->bits, (unsigned long long)L);
}

/* Round m of the breadth-first tests: keep the p in `in` (nb[m] of them) whose member
   p + OFF[m] is a base-2 strong probable prime. */
__global__ void __launch_bounds__(256, 4) round_kernel(int m, const u64 *__restrict__ in, u64 *__restrict__ out,
                                                       Cnt *__restrict__ cnt, u64 cap)
{
    const u64 n0 = cnt->nb[m];
    const u64 n = n0 < cap ? n0 : cap;
    const int lane = threadIdx.x & 31;
    const u64 stride = (u64)gridDim.x * blockDim.x;
    for (u64 base = (u64)blockIdx.x * blockDim.x; base < n; base += stride) {
        const u64 i = base + threadIdx.x;
        bool ok = false;
        u128 p = 0;
        if (i < n) {
            p = ((u128)in[2 * i + 1] << 64) | in[2 * i];
            const u128 v = p + OFFD[m];       /* v < 2^81 (HIMAX) */
            ok = (v >> 64) ? d_sprp2_96(v) : d_sprp2_64((u64)v);
        }
        const unsigned bal = __ballot_sync(0xffffffffu, ok);
        unsigned long long wb = 0;
        if (lane == 0 && bal) wb = atomicAdd(&cnt->nb[m + 1], (unsigned long long)__popc(bal));
        wb = __shfl_sync(0xffffffffu, wb, 0);
        if (ok) {
            const u64 j = wb + __popc(bal & ((1u << lane) - 1));
            out[2 * j] = (u64)p;
            out[2 * j + 1] = (u64)(p >> 64);
        }
    }
}

/* survivors of all six rounds: strong Lucas confirmation, binning */
__global__ void __launch_bounds__(256) cand_kernel(const u64 *__restrict__ in, const Cnt *__restrict__ cnt,
                                                   Res *__restrict__ res)
{
    const u64 n = cnt->nb[6];
    for (u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x; i < n; i += (u64)gridDim.x * blockDim.x)
        d_candidate(((u128)in[2 * i + 1] << 64) | in[2 * i], res);
}

/* self-check kernel: base-2 and Lucas tests on a list of numbers */
__global__ void test_kernel(const u64 *v, int n, u8 *out)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const u128 x = ((u128)v[2 * i + 1] << 64) | v[2 * i];
    out[i] = (u8)((d_sprp2(x) ? 1 : 0) | (d_slprp(x) ? 2 : 0));
}

/* ------------------------------------------------------------------ */
/* Host driver                                                         */
/* ------------------------------------------------------------------ */

typedef struct {
    u64 cnt[MAXBIN], cks[MAXBIN];
    u64 surv, cand, spsp, bits;
} acc_t;

static struct {
    u128 lo, hi, base;
    u64 bmM;
    int nbins;
    u128 bnd[MAXBIN + 1];
} R;

static struct {
    u64 chunk, nchunks, frontier, c0, c1;
    acc_t tot;
    u64 *sample;
    u64 nsample;
} D;

static FILE *g_log = NULL;
static bool g_print = true;
static int g_skip = 0;                       /* bench -X ablation mask */
static bool g_empty = false;                 /* nothing to sieve (HI <= 2^20) */
static u64 g_nspsp_printed = 0;

static void on_sigint(int sig) { (void)sig; g_stop = 1; }

static void acc_add(acc_t *t, const acc_t *a)
{
    for (int b = 0; b < MAXBIN; b++) { t->cnt[b] += a->cnt[b]; t->cks[b] += a->cks[b]; }
    t->surv += a->surv; t->cand += a->cand; t->spsp += a->spsp; t->bits += a->bits;
}

static void log_chunk(u64 idx, const acc_t *a)
{
    /* absolute chunk index when the class range is chunk-aligned (-P parts), else relative */
    const u64 aidx = D.c0 % D.chunk == 0 ? D.c0 / D.chunk + idx : idx;
    fprintf(g_log, "C %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64, aidx, a->surv, a->cand, a->spsp);
    for (int b = 0; b < R.nbins; b++) fprintf(g_log, " %" PRIu64 ":%016" PRIx64, a->cnt[b], a->cks[b]);
    fputc('\n', g_log);
    fflush(g_log);
}

static void small_region(acc_t *A)
{
    if (R.lo >= P0) return;
    const u64 a = (u64)R.lo, b = R.hi < P0 ? (u64)R.hi : P0;
    for (u64 p = a; p < b; p++) {
        if (p % 30 != 7) continue;
        if (is_sext_ref(p)) {
            int k = 0;
            while (k + 1 < R.nbins && p >= R.bnd[k + 1]) k++;
            A->cnt[k]++;
            A->cks[k] += hashp(p);
        }
    }
}

/* device state */
static struct {
    DevT T;
    u32 *dq, *dk, *dminv, *dbmod, *dsegmod, *dpstep, *dpstept, *dpoff;
    u64 *dtab;
    u64 *dcls[2], *hcls[2];
    Res *dres[2], *hres[2];
    u64 *dbuf[2][2];        /* survivor ping-pong buffers per stream (2 u64 per entry) */
    Cnt *dcnt[2];
    u64 bcap;               /* entries per survivor buffer */
    u64 sub;                /* classes per sub-launch (expected survivors fit the buffer) */
    cudaStream_t st[2];
    u64 cap;                /* classes per launch buffer */
} G;

static void dev_free(void)
{
    cudaFree(G.dq); cudaFree(G.dk); cudaFree(G.dminv); cudaFree(G.dbmod); cudaFree(G.dsegmod);
    cudaFree(G.dpstep); cudaFree(G.dpstept); cudaFree(G.dpoff); cudaFree(G.dtab);
    memset(&G.T, 0, sizeof G.T);
    G.dq = G.dk = G.dminv = G.dbmod = G.dsegmod = G.dpstep = G.dpstept = G.dpoff = NULL;
    G.dtab = NULL;
}

template <class T> static T *dev_copy(const T *h, size_t n)
{
    T *d;
    CK(cudaMalloc(&d, sizeof(T) * (n ? n : 1)));
    if (n) CK(cudaMemcpy(d, h, sizeof(T) * n, cudaMemcpyHostToDevice));
    return d;
}

static void dev_upload(void)
{
    dev_free();
    G.dq = dev_copy(SP.q, SP.np);
    G.dk = dev_copy(SP.k, 5 * (size_t)SP.np);
    G.dminv = dev_copy(SP.minv, SP.np);
    G.dbmod = dev_copy(SP.bmod, SP.np);
    G.dsegmod = dev_copy(SP.segmod, SP.np);
    G.dpstep = dev_copy(SP.pstep, MAXPAT);
    G.dpstept = dev_copy(SP.pstept, MAXPAT);
    G.dpoff = dev_copy(SP.poff, MAXPAT);
    G.dtab = dev_copy(SP.tab, SP.ntab);
    G.T.q = G.dq; G.T.k = G.dk; G.T.minv = G.dminv; G.T.bmod = G.dbmod; G.T.segmod = G.dsegmod;
    G.T.pstep = G.dpstep; G.T.pstept = G.dpstept; G.T.poff = G.dpoff; G.T.tab = G.dtab;
    DevP p;
    memset(&p, 0, sizeof p);
    p.M = W.M;
    p.bmM = R.bmM;
    const u128 span = R.hi > R.base ? R.hi - 1 - R.base : 0;
    p.Q0 = R.hi > R.base ? (u64)(span / W.M) : 0;
    p.rem0 = R.hi > R.base ? (u64)(span % W.M) : 0;
    if (R.hi <= R.base) { p.Q0 = 0; p.rem0 = 0; }
    p.base_lo = (u64)R.base; p.base_hi = (u64)(R.base >> 64);
    p.nbins = R.nbins;
    for (int b = 0; b <= R.nbins; b++) { p.bnd_lo[b] = (u64)R.bnd[b]; p.bnd_hi[b] = (u64)(R.bnd[b] >> 64); }
    p.npat = SP.npat; p.nmed = SP.nmed; p.nlarge = SP.nlarge; p.np = SP.np;
    p.skip = g_skip;
    CK(cudaMemcpyToSymbol(P, &p, sizeof p));
}

static void dev_init(int dev)
{
    CK(cudaSetDevice(dev));
    G.cap = 1 << 16;
    for (int i = 0; i < 2; i++) {
        CK(cudaStreamCreate(&G.st[i]));
        CK(cudaMalloc(&G.dcls[i], sizeof(u64) * G.cap));
        CK(cudaMallocHost(&G.hcls[i], sizeof(u64) * G.cap));
        CK(cudaMalloc(&G.dres[i], sizeof(Res)));
        CK(cudaMallocHost(&G.hres[i], sizeof(Res)));
    }
    G.bcap = 1ULL << 25;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) CK(cudaMalloc(&G.dbuf[i][j], sizeof(u64) * 2 * G.bcap));
        CK(cudaMalloc(&G.dcnt[i], sizeof(Cnt)));
    }
    CK(cudaFuncSetAttribute(sext_kernel<12>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)(sizeof(u64) * SEGW)));
    CK(cudaFuncSetAttribute(sext_kernel<24>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)(sizeof(u64) * SEGW)));
    CK(cudaFuncSetAttribute(sext_kernel<MAXLPT>, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)(sizeof(u64) * SEGW)));
}

/* ---- checkpoint (same format as a350826.c) ---- */

static void save_state(const char *path)
{
    char tmp[PATH_MAX];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { fprintf(stderr, "\nwarning: cannot write %s: %s\n", tmp, strerror(errno)); return; }
    fprintf(f, "a350826 checkpoint 1\nlo %s\nhi %s\nwheel %u\nbound %u\nchunk %" PRIu64 "\nclasses %" PRIu64
               " %" PRIu64 "\nbins %d\n", u128s(R.lo), u128s(R.hi), W.w, SP.B, D.chunk, D.c0, D.c1, R.nbins);
    for (int b = 0; b <= R.nbins; b++) fprintf(f, "bnd %s\n", u128s(R.bnd[b]));
    fprintf(f, "frontier %" PRIu64 "\nsurv %" PRIu64 "\ncand %" PRIu64 "\nspsp %" PRIu64 "\nbits %" PRIu64 "\n",
            D.frontier, D.tot.surv, D.tot.cand, D.tot.spsp, D.tot.bits);
    for (int b = 0; b < R.nbins; b++) fprintf(f, "bin %d %" PRIu64 " %016" PRIx64 "\n", b, D.tot.cnt[b], D.tot.cks[b]);
    fputs("end\n", f);
    if (fclose(f) != 0 || rename(tmp, path) != 0)
        fprintf(stderr, "\nwarning: cannot update %s: %s\n", path, strerror(errno));
}

static bool load_state(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        if (errno == ENOENT) return false;
        die("cannot read %s: %s", path, strerror(errno));
    }
    char line[512], slo[64] = "", shi[64] = "";
    if (!fgets(line, sizeof line, f) || strncmp(line, "a350826 checkpoint 1", 20) != 0)
        die("%s is not an a350826 checkpoint", path);
    u32 w = 0, B = 0;
    u64 chunk = 0, frontier = 0, c0 = 0, c1 = 0;
    int nbins = -1, nb = 0;
    u128 bnd[MAXBIN + 1];
    acc_t t;
    memset(&t, 0, sizeof t);
    bool complete = false;
    while (fgets(line, sizeof line, f)) {
        char sv[64];
        int b;
        unsigned long long c, k;
        if (sscanf(line, "lo %63s", slo) == 1) continue;
        if (sscanf(line, "hi %63s", shi) == 1) continue;
        if (sscanf(line, "wheel %u", &w) == 1) continue;
        if (sscanf(line, "bound %u", &B) == 1) continue;
        if (sscanf(line, "chunk %" SCNu64, &chunk) == 1) continue;
        if (sscanf(line, "classes %" SCNu64 " %" SCNu64, &c0, &c1) == 2) continue;
        if (sscanf(line, "bins %d", &nbins) == 1) continue;
        if (sscanf(line, "bnd %63s", sv) == 1) { if (nb <= MAXBIN) bnd[nb++] = parse_num(sv); continue; }
        if (sscanf(line, "frontier %" SCNu64, &frontier) == 1) continue;
        if (sscanf(line, "surv %" SCNu64, &t.surv) == 1) continue;
        if (sscanf(line, "cand %" SCNu64, &t.cand) == 1) continue;
        if (sscanf(line, "spsp %" SCNu64, &t.spsp) == 1) continue;
        if (sscanf(line, "bits %" SCNu64, &t.bits) == 1) continue;
        if (sscanf(line, "bin %d %llu %llx", &b, &c, &k) == 3) {
            if (b < 0 || b >= MAXBIN) die("%s: bad bin line", path);
            t.cnt[b] = c; t.cks[b] = k;
            continue;
        }
        if (!strncmp(line, "end", 3)) { complete = true; break; }
    }
    fclose(f);
    if (!complete) die("%s is truncated", path);
    if (parse_num(slo) != R.lo || parse_num(shi) != R.hi || w != W.w || B != SP.B || chunk != D.chunk ||
        c0 != D.c0 || c1 != D.c1)
        die("%s was written for lo=%s hi=%s w=%u B=%u chunk=%" PRIu64 " classes [%" PRIu64 ", %" PRIu64
            "); rerun with the same parameters", path, slo, shi, w, B, chunk, c0, c1);
    if (nbins != R.nbins || nb != nbins + 1) die("%s was written for different bins", path);
    for (int b = 0; b <= nbins; b++) if (bnd[b] != R.bnd[b]) die("%s was written for different bins", path);
    if (frontier > D.nchunks) die("%s: frontier beyond the last chunk", path);
    D.frontier = frontier;
    D.tot = t;
    return true;
}

/* classes of chunk idx into hcls; returns their number */
static u64 fill_chunk(u64 idx, u64 *h)
{
    const u64 n = D.sample ? D.nsample : D.c1 - D.c0;
    const u64 c0 = idx * D.chunk, c1 = c0 + D.chunk < n ? c0 + D.chunk : n;
    for (u64 c = c0; c < c1; c++) h[c - c0] = class_r(D.sample ? D.sample[c] : D.c0 + c);
    return c1 - c0;
}

static void launch_classes(int s, const u64 *dcls, u64 ncls)
{
    if (!ncls || g_empty) return;
    const size_t smem = sizeof(u64) * SEGW;
    for (u64 a = 0; a < ncls; a += G.sub) {
        const u64 n = a + G.sub < ncls ? G.sub : ncls - a;
        CK(cudaMemsetAsync(G.dcnt[s], 0, sizeof(Cnt), G.st[s]));
        if (SP.nlarge <= NT * 12)
            sext_kernel<12><<<(unsigned)n, NT, smem, G.st[s]>>>(dcls + a, G.T, G.dres[s], G.dcnt[s], G.dbuf[s][0], G.bcap);
        else if (SP.nlarge <= NT * 24)
            sext_kernel<24><<<(unsigned)n, NT, smem, G.st[s]>>>(dcls + a, G.T, G.dres[s], G.dcnt[s], G.dbuf[s][0], G.bcap);
        else
            sext_kernel<MAXLPT><<<(unsigned)n, NT, smem, G.st[s]>>>(dcls + a, G.T, G.dres[s], G.dcnt[s], G.dbuf[s][0],
                                                                  G.bcap);
        CK(cudaGetLastError());
        if (g_skip & 1) continue;
        for (int m = 0; m < 6; m++)
            round_kernel<<<48 * 8, 256, 0, G.st[s]>>>(m, G.dbuf[s][m & 1], G.dbuf[s][(m + 1) & 1], G.dcnt[s], G.bcap);
        cand_kernel<<<48 * 4, 256, 0, G.st[s]>>>(G.dbuf[s][0], G.dcnt[s], G.dres[s]);
        CK(cudaGetLastError());
    }
}

static void launch(int s, u64 ncls)
{
    CK(cudaMemcpyAsync(G.dcls[s], G.hcls[s], sizeof(u64) * ncls, cudaMemcpyHostToDevice, G.st[s]));
    CK(cudaMemsetAsync(G.dres[s], 0, sizeof(Res), G.st[s]));
    launch_classes(s, G.dcls[s], ncls);
    CK(cudaMemcpyAsync(G.hres[s], G.dres[s], sizeof(Res), cudaMemcpyDeviceToHost, G.st[s]));
}

/* a launch overflowed its survivor buffer: redo its classes in pieces, summing the results */
static void redo_in_pieces(int s, u64 ncls, Res *out)
{
    Res *acc = (Res *)calloc(1, sizeof(Res));
    if (!acc) die("out of memory");
    const u64 sub0 = G.sub;
    G.sub = G.sub / 8 ? G.sub / 8 : 1;
    const u64 piece = ncls / 8 ? ncls / 8 : 1;
    for (u64 a = 0; a < ncls; a += piece) {
        const u64 n = a + piece < ncls ? piece : ncls - a;
        CK(cudaMemsetAsync(G.dres[s], 0, sizeof(Res), G.st[s]));
        launch_classes(s, G.dcls[s] + a, n);
        CK(cudaMemcpyAsync(G.hres[s], G.dres[s], sizeof(Res), cudaMemcpyDeviceToHost, G.st[s]));
        CK(cudaStreamSynchronize(G.st[s]));
        const Res *r = G.hres[s];
        if (r->ovf) die("survivor buffer overflow even for %" PRIu64 " classes", n);
        for (int b = 0; b < MAXBIN; b++) { acc->cnt[b] += r->cnt[b]; acc->cks[b] += r->cks[b]; }
        acc->surv += r->surv; acc->cand += r->cand; acc->spsp += r->spsp; acc->bits += r->bits;
        for (unsigned i = 0; i < r->nsp && i < MAXSPSP && acc->nsp < MAXSPSP; i++, acc->nsp++) {
            memcpy(acc->sp_n[acc->nsp], r->sp_n[i], sizeof r->sp_n[i]);
            memcpy(acc->sp_p[acc->nsp], r->sp_p[i], sizeof r->sp_p[i]);
        }
    }
    memcpy(out, acc, sizeof(Res));
    free(acc);
    G.sub = sub0;
}

static void fold(u64 idx, const Res *r)
{
    acc_t a;
    memset(&a, 0, sizeof a);
    for (int b = 0; b < MAXBIN; b++) { a.cnt[b] = r->cnt[b]; a.cks[b] = r->cks[b]; }
    a.surv = r->surv; a.cand = r->cand; a.spsp = r->spsp; a.bits = r->bits;
    if (idx == 0 && !D.sample && D.c0 == 0) small_region(&a);
    const unsigned nsp = r->nsp < MAXSPSP ? r->nsp : MAXSPSP;
    for (unsigned i = 0; i < nsp && g_print; i++) {
        const u128 n = ((u128)r->sp_n[i][1] << 64) | r->sp_n[i][0], p = ((u128)r->sp_p[i][1] << 64) | r->sp_p[i][0];
        if (stderr_tty) fputs("\r\033[K", stderr);
        printf("SPSP %s (base-2 strong pseudoprime in candidate p = %s)\n", u128s(n), u128s(p));
        g_nspsp_printed++;
    }
    acc_add(&D.tot, &a);
    if (g_log) log_chunk(idx, &a);
}

typedef struct { double seconds; u64 bits; } run_stats_t;

static run_stats_t run_count(const char *state, int interval, bool quiet)
{
    D.nchunks = ((D.sample ? D.nsample : D.c1 - D.c0) + D.chunk - 1) / D.chunk;
    if (D.chunk > G.cap) die("chunk larger than the launch buffer (%" PRIu64 ")", G.cap);
    D.frontier = 0;
    memset(&D.tot, 0, sizeof D.tot);
    if (state && load_state(state) && !quiet)
        fprintf(stderr, "resuming from %s: %" PRIu64 " of %" PRIu64 " chunks done\n", state, D.frontier, D.nchunks);
    const u64 bits0 = D.tot.bits, f0 = D.frontier;
    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_sigint;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    const double t0 = now();
    double last_status = 0, last_ckpt = t0;
    char b1[32], b2[32];
    u64 next = D.frontier;
    int inflight = 0;
    u64 idxs[2] = {0, 0};
    int cur = 0;
    /* two launches in flight on alternating streams; fold in order */
    while ((next < D.nchunks && !g_stop) || inflight) {
        while (inflight < 2 && next < D.nchunks && !g_stop) {
            const int s = (cur + inflight) & 1;
            const u64 n = fill_chunk(next, G.hcls[s]);
            launch(s, n);
            idxs[s] = next++;
            inflight++;
        }
        CK(cudaStreamSynchronize(G.st[cur]));
        if (G.hres[cur]->ovf) {
            const u64 n = fill_chunk(idxs[cur], G.hcls[cur]);
            CK(cudaMemcpy(G.dcls[cur], G.hcls[cur], sizeof(u64) * n, cudaMemcpyHostToDevice));
            redo_in_pieces(cur, n, G.hres[cur]);
        }
        fold(idxs[cur], G.hres[cur]);
        D.frontier = idxs[cur] + 1;
        inflight--;
        cur ^= 1;
        const double t = now();
        if (!quiet && (stderr_tty ? t - last_status >= 1.0 : t - last_status >= interval)) {
            last_status = t;
            const double el = t - t0;
            const double rate = el > 0 ? (double)(D.tot.bits - bits0) / el : 0;
            const double crate = el > 0 ? (double)(D.frontier - f0) / el : 0;
            u64 found = 0;
            for (int b = 0; b < R.nbins; b++) found += D.tot.cnt[b];
            fprintf(stderr, "%s%s  chunks %" PRIu64 "/%" PRIu64 " (%.3f%%)  %.3g bits/s  found %" PRIu64 "%s  ETA %s",
                    stderr_tty ? "\r\033[K" : "progress: ", fmt_dur(el, b1, sizeof b1), D.frontier, D.nchunks,
                    100.0 * (double)D.frontier / (double)D.nchunks, rate, found, D.tot.spsp ? "  SPSP!" : "",
                    fmt_dur(crate > 0 ? (double)(D.nchunks - D.frontier) / crate : -1, b2, sizeof b2));
            if (!stderr_tty) fputc('\n', stderr);
            fflush(stderr);
        }
        if (state && t - last_ckpt >= interval) { last_ckpt = t; save_state(state); }
    }
    const double t1 = now();
    if (stderr_tty && !quiet) fputs("\r\033[K", stderr);
    if (state) save_state(state);
    if (g_stop && D.frontier < D.nchunks && !quiet)
        fprintf(stderr, "stopped after %" PRIu64 " of %" PRIu64 " chunks%s\n", D.frontier, D.nchunks,
                state ? ", checkpoint saved" : "");
    run_stats_t st = {t1 - t0, D.tot.bits - bits0};
    return st;
}

/* ------------------------------------------------------------------ */
/* Parameters and commands                                             */
/* ------------------------------------------------------------------ */

static double rho_wheel(u32 w)
{
    double r = 1.0 / 210;
    for (int i = 4; i < 15 && SMALLP[i] <= w; i++) r *= (double)(SMALLP[i] - 6) / SMALLP[i];
    return r;
}

/* GPU cost model, measured on GB10 (B = 2^16): ~1.1 us per class (phases of all sieving primes,
   block start) plus ~5.8 ps per candidate bit for B = 2^16 (scaled by the marks per bit). */
static u32 auto_wheel(u128 lo, u128 hi, u32 B)
{
    const u128 base = lo > P0 ? lo : P0;
    if (hi <= base) return 7;
    int np;
    u32 *pr = small_primes(B, &np);
    const double len = (double)(hi - base);
    u32 best = 7;
    double bestc = 1e300;
    for (int i = 3; i < 15; i++) {
        const u32 w = SMALLP[i];
        double M = 1, C = 1, marks = 0;
        for (int j = 0; j < 15 && SMALLP[j] <= w; j++) M *= SMALLP[j];
        for (int j = 4; j < 15 && SMALLP[j] <= w; j++) C *= SMALLP[j] - 6;
        for (int j = 0; j < np; j++) if (pr[j] > w) marks += 6.0 / pr[j];
        if (len / M < 1 && w > 7) break;
        const double cost = C * 1.1e-6 * (np / 6542.0) + len * rho_wheel(w) * 5.8e-12 * (marks / 6.45);
        if (cost < bestc) { bestc = cost; best = w; }
    }
    free(pr);
    return best;
}

static u64 auto_chunk(u64 C)
{
    u64 c = 1;
    while (c * 2 <= C / 4096 && c < 16384) c *= 2;
    return c;
}

typedef struct {
    u32 w, B, pmax, qw;
    u64 chunk, nsample, cr0, cr1, part, nparts;
    int interval, dev;
    const char *state, *log;
    bool quiet;
    int nb;
    u128 bnd[MAXBIN];
} opts_t;

static void usage(void)
{
    fputs("usage: a350826_cuda count LO HI [-b B1,B2,..] [-w W] [-B B] [-c CH] [-C C0:C1 | -P I/N] [-S FILE] [-i SECS]\n"
          "                          [-L FILE] [-q] [-d DEV]\n"
          "       a350826_cuda bench LO HI [-n N] [-w W] [-B B] [-d DEV]\n"
          "       a350826_cuda selftest [N] [-d DEV]\n", stderr);
    exit(2);
}

static void parse_opts(int argc, char **argv, int first, opts_t *o, int npos, u128 *pos)
{
    memset(o, 0, sizeof *o);
    o->B = 1u << 16;
    o->pmax = 256;
    o->qw = 4096;
    o->interval = 60;
    o->nsample = 20000;
    int np = 0;
    for (int i = first; i < argc; i++) {
        const char *a = argv[i];
        if (a[0] == '-' && a[1] && !isdigit((unsigned char)a[1]) && a[2] == 0) {
            const char c = a[1];
            if (c == 'q') { o->quiet = true; continue; }
            if (i + 1 >= argc) usage();
            const char *v = argv[++i];
            switch (c) {
            case 'b': {
                char buf[1024];
                snprintf(buf, sizeof buf, "%s", v);
                for (char *tok = strtok(buf, ","); tok; tok = strtok(NULL, ",")) {
                    if (o->nb >= MAXBIN - 1) die("too many bin boundaries");
                    o->bnd[o->nb++] = parse_num(tok);
                }
                break;
            }
            case 'w': o->w = (u32)parse_u64(v, 7, 47, "wheel prime"); break;
            case 'B': o->B = (u32)parse_u64(v, 11, BMAX, "sieving bound"); break;
            case 'p': o->pmax = (u32)parse_u64(v, 0, 1024, "pattern bound"); break;
            case 'Q': o->qw = (u32)parse_u64(v, 0, BMAX, "medium bound"); break;
            case 'c': o->chunk = parse_u64(v, 1, 1ULL << 16, "chunk"); break;
            case 'n': o->nsample = parse_u64(v, 1, 1ULL << 40, "sample"); break;
            case 'i': o->interval = (int)parse_u64(v, 1, 86400, "interval"); break;
            case 'd': o->dev = (int)parse_u64(v, 0, 64, "device"); break;
            case 'X': g_skip = (int)parse_u64(v, 0, 31, "skip mask"); break;
            case 'S': o->state = v; break;
            case 'L': o->log = v; break;
            case 'C': {
                const char *colon = strchr(v, ':');
                if (!colon) die("-C wants C0:C1");
                char tmp[64];
                snprintf(tmp, sizeof tmp, "%.*s", (int)(colon - v), v);
                o->cr0 = (u64)parse_num(tmp);
                o->cr1 = (u64)parse_num(colon + 1);
                if (o->cr1 <= o->cr0) die("-C wants C0 < C1");
                break;
            }
            case 'P': {
                const char *sl = strchr(v, '/');
                if (!sl) die("-P wants I/N");
                char tmp[64];
                snprintf(tmp, sizeof tmp, "%.*s", (int)(sl - v), v);
                o->part = (u64)parse_num(tmp);
                o->nparts = (u64)parse_num(sl + 1);
                if (o->nparts == 0 || o->part >= o->nparts) die("-P wants 0 <= I < N");
                break;
            }
            default: usage();
            }
            continue;
        }
        if (np >= npos) usage();
        pos[np++] = parse_num(a);
    }
    if (np != npos) usage();
}

static void configure(opts_t *o, u128 lo, u128 hi)
{
    if (hi <= lo) die("HI must be larger than LO");
    if (hi > HIMAX) die("HI must be <= 2^80");
    R.lo = lo; R.hi = hi;
    R.base = lo > P0 ? lo : P0;
    const u32 w = o->w ? o->w : auto_wheel(lo, hi, o->B);
    wheel_setup(w);
    if (o->B <= W.w) die("sieving bound must exceed the wheel prime");
    R.bmM = (u64)(R.base % W.M);
    R.nbins = 0;
    R.bnd[R.nbins++] = lo;
    for (int i = 0; i < o->nb; i++) {
        if (o->bnd[i] <= R.bnd[R.nbins - 1] || o->bnd[i] >= hi) die("bin boundaries must increase inside (LO, HI)");
        R.bnd[R.nbins++] = o->bnd[i];
    }
    R.bnd[R.nbins] = hi;
    tables_setup(o->B, o->pmax, o->qw, R.base);
    const u128 span = R.hi > R.base ? (R.hi - 1 - R.base) / W.M : 0;
    if (span >> 62) die("class too long; use a larger wheel (-w)");
    g_empty = R.hi <= R.base;
    {                                         /* expected survivors per class -> classes per sub-launch */
        double sr = 1;
        for (int i = 0; i < SP.np; i++) sr *= 1.0 - 6.0 / SP.q[i];
        const double per = ((double)span + 1) * sr * 1.3 + 64;
        const double sub = (double)G.bcap / per;
        G.sub = sub < 1 ? 1 : sub > 65536 ? 65536 : (u64)sub;
    }
    dev_upload();
    D.chunk = o->chunk ? o->chunk : auto_chunk(W.C);
    D.sample = NULL;
    D.c0 = 0;
    D.c1 = W.C;
    if (o->cr1) {
        if (o->cr1 > W.C) die("-C range beyond the %" PRIu64 " classes of wheel %u", W.C, W.w);
        D.c0 = o->cr0;
        D.c1 = o->cr1;
    } else if (o->nparts) {
        const u64 nch = (W.C + D.chunk - 1) / D.chunk;
        const u64 a = (u64)((u128)nch * o->part / o->nparts), b = (u64)((u128)nch * (o->part + 1) / o->nparts);
        D.c0 = a * D.chunk;
        D.c1 = b * D.chunk < W.C ? b * D.chunk : W.C;
        if (D.c1 <= D.c0) die("part %" PRIu64 "/%" PRIu64 " is empty", o->part, o->nparts);
    }
}

static const char *decade_label(u128 a, u128 b)
{
    static char buf[32];
    u128 x = 1;
    for (int n = 1; n <= 24; n++) {
        if (a == x && b == x * 10) { snprintf(buf, sizeof buf, "A350826(%d)", n); return buf; }
        x *= 10;
    }
    return NULL;
}

static void print_results(const char *tag)
{
    u64 total = 0, cks = 0;
    for (int b = 0; b < R.nbins; b++) { total += D.tot.cnt[b]; cks += D.tot.cks[b]; }
    for (int b = 0; b < R.nbins; b++) {
        const char *lab = decade_label(R.bnd[b], R.bnd[b + 1]);
        printf("%s [%s, %s) count %" PRIu64 " cks %016" PRIx64 "%s%s\n", tag, u128s(R.bnd[b]), u128s(R.bnd[b + 1]),
               D.tot.cnt[b], D.tot.cks[b], lab ? "  " : "", lab ? lab : "");
    }
    if (R.nbins > 1)
        printf("%s [%s, %s) count %" PRIu64 " cks %016" PRIx64 "\n", tag, u128s(R.lo), u128s(R.hi), total, cks);
    printf("%s survivors %" PRIu64 " candidates %" PRIu64 " spsp %" PRIu64 " bits %" PRIu64 "\n", tag,
           D.tot.surv, D.tot.cand, D.tot.spsp, D.tot.bits);
    fflush(stdout);
}

static void print_banner(const opts_t *o)
{
    cudaDeviceProp prop;
    CK(cudaGetDeviceProperties(&prop, o->dev));
    fprintf(stderr, "[%s, %s) on %s: wheel %u (%" PRIu64 " classes, ~%.3g candidates), B = %u (%d pattern, %d medium, "
                    "%d large primes), %" PRIu64 " chunks of %" PRIu64 "\n",
            u128s(R.lo), u128s(R.hi), prop.name, W.w, W.C, (double)(R.hi - R.base) * rho_wheel(W.w), SP.B, SP.npat,
            SP.nmed, SP.nlarge, (D.c1 - D.c0 + D.chunk - 1) / D.chunk, D.chunk);
    if (D.c0 != 0 || D.c1 != W.C)
        fprintf(stderr, "classes [%" PRIu64 ", %" PRIu64 ") of %" PRIu64 " only (partial count)\n", D.c0, D.c1, W.C);
}

static int cmd_count(int argc, char **argv)
{
    opts_t o;
    u128 pos[2];
    parse_opts(argc, argv, 2, &o, 2, pos);
    dev_init(o.dev);
    configure(&o, pos[0], pos[1]);
    if (o.log) {
        g_log = fopen(o.log, "a");
        if (!g_log) die("cannot open %s: %s", o.log, strerror(errno));
    }
    if (!o.quiet) print_banner(&o);
    run_stats_t st = run_count(o.state, o.interval, o.quiet);
    if (g_log) fclose(g_log);
    if (D.frontier < D.nchunks) { print_results("PARTIAL"); return 1; }
    print_results(D.c0 == 0 && D.c1 == W.C ? "RESULT" : "PART");
    char b1[32];
    if (!o.quiet)
        fprintf(stderr, "done in %s (%.3g bits/s)\n", fmt_dur(st.seconds, b1, sizeof b1),
                st.seconds > 0 ? (double)st.bits / st.seconds : 0);
    return 0;
}

static int cmd_bench(int argc, char **argv)
{
    opts_t o;
    u128 pos[2];
    parse_opts(argc, argv, 2, &o, 2, pos);
    dev_init(o.dev);
    configure(&o, pos[0], pos[1]);
    const u64 n = o.nsample < W.C ? o.nsample : W.C;
    u64 *smp = (u64 *)malloc(sizeof *smp * n);
    if (!smp) die("out of memory");
    for (u64 i = 0; i < n; i++) smp[i] = (u64)((u128)i * W.C / n);
    D.sample = smp;
    D.nsample = n;
    D.chunk = n < 4096 ? n : 4096;
    g_print = false;
    print_banner(&o);
    /* warm-up launch */
    {
        const u64 m = fill_chunk(0, G.hcls[0]);
        launch(0, m < 256 ? m : 256);
        CK(cudaStreamSynchronize(G.st[0]));
    }
    run_stats_t st = run_count(NULL, 60, true);
    const double scale = (double)W.C / (double)n;
    char b1[32];
    printf("BENCH w %u B %u: %.3f s for %" PRIu64 " classes, %.4g bits/s\n", W.w, SP.B, st.seconds, n,
           (double)st.bits / st.seconds);
    printf("BENCH per class: %.4g bits, %.4g survivors, %.4g candidates; survivor rate %.3e\n",
           (double)D.tot.bits / n, (double)D.tot.surv / n, (double)D.tot.cand / n,
           (double)D.tot.surv / (double)D.tot.bits);
    u64 found = 0;
    for (int b = 0; b < R.nbins; b++) found += D.tot.cnt[b];
    printf("BENCH projection for all %" PRIu64 " classes: %s (~%.4g sextuplets)\n", W.C,
           fmt_dur(st.seconds * scale, b1, sizeof b1), (double)found * scale);
    free(smp);
    D.sample = NULL;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Selftest                                                            */
/* ------------------------------------------------------------------ */

static u64 rng_state = 0x243F6A8885A308D3ULL;
static u64 rnd(void) { return mix64(rng_state += GOLD); }

static int test_arith_gpu(void)
{
    /* random odd numbers up to 2^80 (no factor <= 7), known pseudoprimes; compare with the host */
    enum { N = 1 << 16 };
    static u64 v[2 * N];
    static u8 out[N];
    static u128 xs[N];
    int n = 0;
    static const u64 fixed[] = {2047, 3277, 4033, 4681, 8321, 15841, 29341, 42799, 49141, 52633, 65281, 74665,
                                5459, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519, 75077, 97439,
                                18446744073709551557ULL, 18446744073709551533ULL, 18446744073709551521ULL};
    for (size_t i = 0; i < sizeof fixed / sizeof *fixed; i++) xs[n++] = fixed[i];
    while (n < N) {
        const int bits = 20 + (int)(rnd() % 61);
        u128 x = ((u128)rnd() << 64 | rnd()) >> (128 - bits);
        x |= ((u128)1 << (bits - 1)) | 1;
        if (x % 3 == 0 || x % 5 == 0 || x % 7 == 0 || x < 50) continue;
        xs[n++] = x;
    }
    if (getenv("A350826_NTEST")) {
        n = atoi(getenv("A350826_NTEST"));
        if (getenv("A350826_FROM")) {
            const int f = atoi(getenv("A350826_FROM"));
            for (int i = 0; i < n; i++) xs[i] = xs[f + i];
        }
    }
    for (int i = 0; i < n; i++) { v[2 * i] = (u64)xs[i]; v[2 * i + 1] = (u64)(xs[i] >> 64); }
    u64 *dv;
    u8 *dout;
    CK(cudaMalloc(&dv, sizeof v));
    CK(cudaMalloc(&dout, sizeof out));
    CK(cudaMemcpy(dv, v, sizeof v, cudaMemcpyHostToDevice));
    test_kernel<<<(n + 255) / 256, 256>>>(dv, n, dout);
    CK(cudaGetLastError());
    CK(cudaMemcpy(out, dout, sizeof out, cudaMemcpyDeviceToHost));
    cudaFree(dv);
    cudaFree(dout);
    int bad = 0, primes = 0;
    for (int i = 0; i < n; i++) {
        const u128 x = xs[i];
        const bool p = is_prime_ref(x), s2 = mr_ref(x, 2);
        primes += p;
        const bool g2 = out[i] & 1, gl = (out[i] >> 1) & 1;
        if (g2 != s2) { if (bad < 10) printf("FAIL gpu sprp2(%s) = %d\n", u128s(x), g2); bad++; }
        if (p && !gl) { if (bad < 10) printf("FAIL gpu slprp(%s) rejects a prime\n", u128s(x)); bad++; }
        if (!p && g2 && gl) { if (bad < 10) printf("FAIL gpu BPSW accepts composite %s\n", u128s(x)); bad++; }
    }
    /* known strong Lucas pseudoprimes must pass slprp */
    for (int i = 12; i < 24; i++) if (!((out[i] >> 1) & 1)) { printf("FAIL gpu slprp(%s) should pass\n", u128s(xs[i])); bad++; }
    printf("gpu arithmetic: %s (%d numbers, %d primes)\n", bad ? "FAILED" : "ok", n, primes);
    return bad;
}

static const u64 A350826[] = {0, 1, 1, 0, 0, 3, 0, 13, 64, 235, 1296, 7013, 41782, 253420, 1607418, 10520883,
                              70785653, 488096844};

/* windows computed by verify_a350826.py (files in verify/): count, cks */
typedef struct { const char *lo, *hi; u64 count, cks; } window_t;
static const window_t WINDOWS[] = {
    {"1e15", "1e15+2e12", 20102, 0x698702d7f0805d7dULL},
    {"1e17-1e12", "1e17+1e12", 9407, 0xb42353cfea760f5fULL},
    {"1e18-1e12", "1e18+1e12", 6983, 0x67b3c4c4b0bfba36ULL},
    {"1e19-1e12", "1e19+1e12", 4881, 0x7b9d093b40d59fb3ULL},
    {"2^64-1e12", "2^64+1e12", 4488, 0xfb6f2691fb75e673ULL},
    {"1e20-1e12", "1e20+1e12", 3644, 0x398796b7104319cbULL},
    {"1e22", "1e22+2e12", 2059, 0x004f163c8720cbb2ULL},
    {NULL, NULL, 0, 0}
};

static void run_quiet(u128 lo, u128 hi, u32 w, u32 B, u64 c0, u64 c1, u64 *count, u64 *cks)
{
    opts_t o;
    memset(&o, 0, sizeof o);
    o.w = w; o.B = B; o.pmax = 256; o.qw = 4096; o.interval = 60;
    o.cr0 = c0; o.cr1 = c1;
    configure(&o, lo, hi);
    run_count(NULL, 60, true);
    *count = 0; *cks = 0;
    for (int b = 0; b < R.nbins; b++) { *count += D.tot.cnt[b]; *cks += D.tot.cks[b]; }
}

static int cmd_selftest(int argc, char **argv)
{
    opts_t o;
    u128 pos[1] = {14};
    int npos = argc > 2 && isdigit((unsigned char)argv[2][0]) ? 1 : 0;
    parse_opts(argc, argv, 2, &o, npos, pos);
    const int N = (int)pos[0];
    if (N < 1 || N > 17) die("N must be 1..17");
    dev_init(o.dev);
    g_print = false;
    const double t0 = now();
    int bad = test_arith_gpu();
    u128 lo = 1;
    for (int n = 1; n <= N; n++) {
        u64 c, k;
        run_quiet(lo, lo * 10, 0, 1u << 16, 0, 0, &c, &k);
        const bool ok = c == A350826[n];
        printf("A350826(%d) = %" PRIu64 " %s (cks %016" PRIx64 ", w %u, %.1fs)\n", n, c, ok ? "ok" : "FAIL", k, W.w,
               now() - t0);
        fflush(stdout);
        bad += !ok;
        lo *= 10;
    }
    for (size_t i = 0; WINDOWS[i].lo; i++) {
        u64 c, k;
        run_quiet(parse_num(WINDOWS[i].lo), parse_num(WINDOWS[i].hi), 0, 1u << 16, 0, 0, &c, &k);
        const bool ok = c == WINDOWS[i].count && k == WINDOWS[i].cks;
        printf("window [%s, %s): %" PRIu64 " %016" PRIx64 " %s\n", WINDOWS[i].lo, WINDOWS[i].hi, c, k, ok ? "ok" : "FAIL");
        fflush(stdout);
        bad += !ok;
    }
    char b1[32];
    printf("selftest %s (%s)\n", bad ? "FAILED" : "passed", fmt_dur(now() - t0, b1, sizeof b1));
    return bad ? 1 : 0;
}

static int cmd_mrbench(int argc, char **argv)      /* round-kernel throughput on random odd numbers near X */
{
    if (argc < 3) usage();
    const u128 x0 = parse_num(argv[2]);
    const u64 n = argc > 3 ? (u64)parse_num(argv[3]) : (1ULL << 24);
    dev_init(0);
    u64 *h = (u64 *)malloc(sizeof(u64) * 2 * n);
    for (u64 i = 0; i < n; i++) {
        const u128 v = (x0 + ((u128)(rnd() >> 8) << 1)) | 1;
        h[2 * i] = (u64)v; h[2 * i + 1] = (u64)(v >> 64);
    }
    CK(cudaMemcpy(G.dbuf[0][0], h, sizeof(u64) * 2 * n, cudaMemcpyHostToDevice));
    Cnt c;
    memset(&c, 0, sizeof c);
    c.nb[0] = n;
    for (int grid = 48 * 4; grid <= 48 * 32; grid *= 2) {
        CK(cudaMemcpy(G.dcnt[0], &c, sizeof c, cudaMemcpyHostToDevice));
        CK(cudaDeviceSynchronize());
        const double t0 = now();
        round_kernel<<<grid, 256>>>(0, G.dbuf[0][0], G.dbuf[0][1], G.dcnt[0], G.bcap);
        CK(cudaDeviceSynchronize());
        const double t = now() - t0;
        Cnt r;
        CK(cudaMemcpy(&r, G.dcnt[0], sizeof r, cudaMemcpyDeviceToHost));
        printf("grid %d: %.3f s, %.3g tests/s, %" PRIu64 " passed\n", grid, t, (double)n / t, (u64)r.nb[1]);
    }
    free(h);
    return 0;
}

static int cmd_num(int argc, char **argv)          /* debug: GPU tests of single numbers */
{
    dev_init(0);
    for (int i = 2; i < argc; i++) {
        const u128 x = parse_num(argv[i]);
        u64 v[2] = {(u64)x, (u64)(x >> 64)}, *dv;
        u8 out = 0, *dout;
        CK(cudaMalloc(&dv, sizeof v));
        CK(cudaMalloc(&dout, 1));
        CK(cudaMemcpy(dv, v, sizeof v, cudaMemcpyHostToDevice));
        test_kernel<<<1, 1>>>(dv, 1, dout);
        CK(cudaGetLastError());
        CK(cudaMemcpy(&out, dout, 1, cudaMemcpyDeviceToHost));
        printf("%s: gpu sprp2 %d slprp %d | host mr2 %d prime %d\n", u128s(x), out & 1, (out >> 1) & 1,
               (int)mr_ref(x, 2), (int)is_prime_ref(x));
        fflush(stdout);
        cudaFree(dv); cudaFree(dout);
    }
    return 0;
}

int main(int argc, char **argv)
{
    stderr_tty = isatty(2);
    setvbuf(stdout, NULL, _IOLBF, 0);
    if (argc < 2) usage();
    if (!strcmp(argv[1], "count")) return cmd_count(argc, argv);
    if (!strcmp(argv[1], "bench")) return cmd_bench(argc, argv);
    if (!strcmp(argv[1], "selftest")) return cmd_selftest(argc, argv);
    if (!strcmp(argv[1], "num")) return cmd_num(argc, argv);
    if (!strcmp(argv[1], "mrbench")) return cmd_mrbench(argc, argv);
    usage();
    return 2;
}
