/*
 * a350826.c
 *
 * Counts prime sextuplets (p, p+4, p+6, p+10, p+12, p+16) by their initial
 * member p, to extend
 *
 *   A350826  number of prime sextuplets with n-digit initial term,
 *   A063501  number of prime sextuplets up to 10^n (A022008 = initial members).
 *
 * OEIS (Jan 2022) lists A350826 through a(17) = 488096844 (sextuplets with
 * 10^16 < p < 10^17).  N. Luhn's table (pzktupel.de/counting/PI_06.php) adds
 * pi_6(10^18) = 4010758480, pi_6(10^19) = 28722086297 and pi_6(2^64) =
 * 48629687343 (K. Desfontaines, May/Jul 2026), which imply a(18) = 3439443854
 * and a(19) = 24711327817; a(20) = pi_6(10^20) - pi_6(10^19) is unknown
 * (Hardy-Littlewood: ~1.806e11).
 *
 * Method
 * ------
 * Every sextuplet except (7, ..., 23) has p = 97 (mod 210), and for a prime
 * q >= 11 the residue p mod q must avoid the six values -d (mod q), d in
 * {0,4,6,10,12,16} (distinct for q >= 11).  A wheel M = 2*3*5*...*w (largest
 * wheel prime w, chosen from the range by a cost model) leaves
 *
 *     C(w) = prod_{11 <= q <= w} (q - 6)          residue classes r mod M
 *
 * (85085 for w = 23, 1.5e9 for w = 37), a fraction rho(w) = 2.04e-4 of all
 * integers for w = 37.  Each class is the progression p = p0 + M k; it is
 * sieved in equal segments of at most S bits (bit k = candidate p0 + M k) with
 * every prime w < q <= B (default B = 2^16).  For each sieving prime the six
 * bad k form six residues phi + K_j (mod q), where K_0 = 0 < K_1 < ... < K_5 < q
 * are the class-independent offsets -d M^-1 (mod q) and only the phase phi
 * depends on the class.  Primes q < 256 (-p) are applied as precomputed word
 * patterns (T_q[s] = the 64 bits starting at offset s of the q-periodic bad set,
 * four primes ORed per pass), which avoids the store-forwarding stalls of
 * marking a word several times; the larger primes mark one "round" (one period
 * of q, six bits) per loop iteration, and phi carries over between segments.
 *
 * Survivors (about 1.4e-3 of the candidates for w = 37, B = 2^16) are tested
 * member by member with a base-2 strong probable-prime test: all survivors of a
 * segment for p, those passing for p+4, and so on, four Montgomery
 * exponentiations in lock step (64-bit below 2^64, two-limb R = 2^128 above).
 * When all six pass, each member is confirmed with a strong Lucas test
 * (Selfridge parameters), i.e. BPSW, which has no counterexample below 2^64
 * (Feitsma-Galway list of base-2 pseudoprimes) and none known at all; a member
 * that passes base 2 but fails Lucas is reported as an SPSP line and not
 * counted.  Initial members below 2^20 are tested directly (this includes p = 7).
 *
 * The classes are cut into chunks handed out through an atomic counter;
 * chunk results are folded in chunk order, so the checkpoint is a frontier.
 * Each found p goes into the count and the order-independent checksum
 *
 *     cks = sum of mix64(lo64(p) ^ hi64(p) * 0x9E3779B97F4A7C15) mod 2^64
 *
 * (mix64 = splitmix64 finalizer) of the bin it falls in (bins are given with
 * -b), so runs with any wheel, bound, segment size or thread count - and the
 * CUDA version and verify_a350826.py - give identical totals.  With -L FILE
 * one line per chunk (survivors, candidates, counts and checksums per bin) is
 * appended in frontier order, for comparing runs with the same w, B, chunk.
 *
 * Usage
 * -----
 *   a350826 count LO HI [-b B1,B2,..] [-w W] [-B B] [-s S] [-c CH] [-t T]
 *                       [-S FILE] [-i SECS] [-L FILE] [-q]
 *       Count sextuplets with LO <= p < HI.  -b splits the range into bins at
 *       the given boundaries (e.g. -b 1e18,2^64); counts and checksums are
 *       reported per bin, and bins [10^(n-1), 10^n) are labelled A350826(n).
 *       -w largest wheel prime (7..47), -B largest sieving prime (<= 2^20),
 *       -s segment bits (default 2^19; 2^20 is best on Apple P-cores), -p
 *       pattern bound (256), -c classes per chunk, -C C0:C1 only classes
 *       [C0, C1), -P I/N only part I of N (chunk-aligned; for splitting a run
 *       over machines), -t threads (default: all cores), -S checkpoint (every
 *       SECS s, default 60, and on Ctrl-C; rerunning the same command resumes),
 *       -L chunk log, -q quiet.
 *   a350826 list LO HI           print the initial members in [LO, HI) (sorted)
 *   a350826 check P              test one p with the fast and reference code
 *   a350826 bench LO HI [-n N] [-w W] [-B B] [-s S] [-t T] [-x]
 *       Sieve N classes (default 2000) spread over all classes, over the full
 *       range, and project the time of the whole count; -x skips the
 *       primality tests (sieve only).
 *   a350826 recheck LO HI CHUNKLOG -w W -c CH [-B B] [-b B1,..] [-n K] [-t T]
 *       Recompute K random chunks (default: all) of a chunk log written by
 *       either tool for the same range, wheel, bound, chunk size and bins, and
 *       compare the lines (spot check of a long GPU run on the CPU).
 *   a350826 selftest [N] [-t T]
 *       Arithmetic against slow references and known pseudoprimes, small
 *       ranges against a direct search, A350826(1..N) (default N = 14)
 *       against the OEIS data, wheel/bound/segment invariance, and windows
 *       near 1e17, 1e18, 1e19, 2^64 and 1e20 against verify_a350826.py.
 *
 * Numbers may be written as decimal, 2^k, 10^k, 1e18, 1.5e19, or X+Y / X-Y
 * of those.
 *
 * Build:  cc -O3 -std=gnu11 -pthread a350826.c -lm -o a350826
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
#include <inttypes.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include <time.h>

typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t i64;
typedef unsigned __int128 u128;

#define GOLD    0x9E3779B97F4A7C15ULL
#define MAXBIN  32
#define RINGC   8192            /* chunk completion ring: bounds how far workers run ahead */
#define MAXSPSP 4096
#define P0      ((u64)1 << 20)  /* initial members below P0 are tested directly */
#define BMAX    ((u32)1 << 20)  /* largest sieving bound */
#define HIMAX   ((u128)1 << 80) /* < 3.3e24, where the reference test is deterministic */

static const u32 OFF[6] = {0, 4, 6, 10, 12, 16};

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
    if (s < 0 || s != s || s > 1e12) { snprintf(buf, n, "?"); return buf; }
    long t = (long)s;
    if (t >= 86400) snprintf(buf, n, "%ldd%02ldh%02ldm", t / 86400, (t / 3600) % 24, (t / 60) % 60);
    else if (t >= 3600) snprintf(buf, n, "%ldh%02ldm%02lds", t / 3600, (t / 60) % 60, t % 60);
    else if (t >= 60) snprintf(buf, n, "%ldm%02lds", t / 60, t % 60);
    else snprintf(buf, n, "%.1fs", s);
    return buf;
}

/* decimal string of a u128 (ring of buffers; main thread or under g_mu only) */
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

/* decimal, B^E, MeE (M may have a decimal point) */
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
        for (size_t i = 0; i < ip; i++) if (!mul_ok(m, 10, &m)) return false; else m += (u128)(buf[i] - '0');
        for (const char *f = frac; *f; f++) if (!mul_ok(m, 10, &m)) return false; else m += (u128)(*f - '0');
        for (; e > 0; e--) if (!mul_ok(m, 10, &m)) return false;
        for (; e < 0; e++) { if (m % 10) return false; m /= 10; }
        *out = m;
        return true;
    }
    if (!all_digits(buf, n)) return false;
    u128 m = 0;
    for (size_t i = 0; i < n; i++) if (!mul_ok(m, 10, &m)) return false; else m += (u128)(buf[i] - '0');
    *out = m;
    return true;
}

static u128 parse_num(const char *s)
{
    size_t n = strlen(s);
    for (size_t i = n; i-- > 1; ) {           /* split at the last top-level + or - */
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

static int default_threads(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? (int)n : 1;
}

static inline u64 mix64(u64 z)
{
    z += GOLD;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline u64 hashp(u128 p) { return mix64((u64)p ^ ((u64)(p >> 64) * GOLD)); }

static u64 rnd_r_state = 0;                  /* recheck sampling: seeded from the clock */
static u64 rnd_r(void)
{
    if (!rnd_r_state) rnd_r_state = (u64)time(NULL) * GOLD | 1;
    return mix64(rnd_r_state += GOLD);
}

/* ------------------------------------------------------------------ */
/* Reference arithmetic (slow and simple; checks and selftest only)    */
/* ------------------------------------------------------------------ */

/* a*b mod n for n < 2^127; the loop runs over the bits of b */
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

static bool mr_ref(u128 n, u64 a)            /* strong test to base a; n odd > a */
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

/* Deterministic for n < 3317044064679887385961981 (first 13 prime bases, Sorenson-Webster 2015). */
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
/* Fast primality, 64-bit: Montgomery with R = 2^64                    */
/* ------------------------------------------------------------------ */

static inline u64 inv64(u64 n)                /* n^-1 mod 2^64, n odd */
{
    u64 x = n;                                /* correct to 3 bits */
    for (int i = 0; i < 5; i++) x *= 2 - n * x;
    return x;
}

/* a*b/R mod n for a, b < n < 2^64, ni = n^-1 mod 2^64: (t - m n)/R lies in (-n, n) */
static inline u64 mmul64(u64 a, u64 b, u64 n, u64 ni)
{
    u128 t = (u128)a * b;
    u64 m = (u64)t * ni;
    u64 mh = (u64)(((u128)m * n) >> 64);
    u64 th = (u64)(t >> 64);
    u64 r = th - mh;
    return th < mh ? r + n : r;
}

static inline u64 madd64(u64 a, u64 b, u64 n) { u64 s = a + b; return (s < a || s >= n) ? s - n : s; }
static inline u64 msub64(u64 a, u64 b, u64 n) { return a >= b ? a - b : a - b + n; }
static inline u64 mhalf64(u64 a, u64 n) { return (a & 1) ? (a >> 1) + (n >> 1) + 1 : a >> 1; }

/* strong probable prime to base 2; n odd, n > 2 */
static inline bool sprp2_64(u64 n)
{
    const u64 ni = inv64(n);
    const u64 one = (u64)(-n) % n, mone = n - one;
    u64 d = n - 1;
    const int s = __builtin_ctzll(d);
    d >>= s;
    int b = 63 - __builtin_clzll(d);
    u64 x = madd64(one, one, n);              /* top bit of d */
    while (b-- > 0) {
        x = mmul64(x, x, n, ni);
        if ((d >> b) & 1) x = madd64(x, x, n);
    }
    if (x == one || x == mone) return true;
    for (int i = 1; i < s; i++) {
        x = mmul64(x, x, n, ni);
        if (x == mone) return true;
        if (x == one) return false;
    }
    return false;
}

/* Four base-2 strong tests in lock step (independent multiply chains overlap in the core).
   Lanes whose exponent d is shorter start with leading zero bits (x stays at R mod n). */
static void sprp2_64x4(const u64 *n, bool *res)
{
    u64 ni[4], one[4], mone[4], d[4], x[4];
    int s[4], bmax = 0;
    for (int l = 0; l < 4; l++) {
        ni[l] = inv64(n[l]);
        one[l] = (u64)(-n[l]) % n[l];
        mone[l] = n[l] - one[l];
        d[l] = n[l] - 1;
        s[l] = __builtin_ctzll(d[l]);
        d[l] >>= s[l];
        const int b = 63 - __builtin_clzll(d[l]);
        if (b > bmax) bmax = b;
        x[l] = one[l];
    }
    for (int b = bmax; b >= 0; b--) {
        for (int l = 0; l < 4; l++) {
            u64 y = mmul64(x[l], x[l], n[l], ni[l]);
            if ((d[l] >> b) & 1) y = madd64(y, y, n[l]);
            x[l] = y;
        }
    }
    for (int l = 0; l < 4; l++) {
        bool r = x[l] == one[l] || x[l] == mone[l];
        for (int i = 1; i < s[l] && !r; i++) {
            x[l] = mmul64(x[l], x[l], n[l], ni[l]);
            if (x[l] == mone[l]) r = true;
            else if (x[l] == one[l]) break;
        }
        res[l] = r;
    }
}

/* Jacobi symbol (a/n), n odd positive */
static int jacobi(i64 a, u128 n)
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

static bool is_square(u128 n)
{
    u128 r = (u128)sqrtl((long double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r * r == n;
}

/* Selfridge parameters: first D in 5, -7, 9, -11, ... with (D/n) = -1.  Returns 0 if n is
   composite for a trivial reason (square or a factor in common with D), else 1. */
static int selfridge(u128 n, i64 *Dout)
{
    if (is_square(n)) return 0;
    i64 D = 5;
    for (;;) {
        int j = jacobi(D, n);
        if (j == -1) break;
        if (j == 0 && (u128)(D < 0 ? -D : D) != n) return 0;
        D = D > 0 ? -(D + 2) : -(D - 2);
    }
    *Dout = D;
    return 1;
}

/* strong Lucas probable prime with P = 1, Q = (1 - D)/4; n odd, n > 13, not a square */
static bool slprp64(u64 n)
{
    i64 D;
    if (!selfridge(n, &D)) return false;
    const i64 Q = (1 - D) / 4;
    const u64 ni = inv64(n), one = (u64)(-n) % n;
#define TOM64(v) ((u64)(((u128)(v) * one) % n))
    const u64 Qm = TOM64(Q >= 0 ? (u64)Q % n : n - (u64)(-Q) % n);
    const u64 Dm = TOM64(D >= 0 ? (u64)D % n : n - (u64)(-D) % n);
#undef TOM64
    u64 d = n + 1;
    const int s = __builtin_ctzll(d);
    d >>= s;
    u64 U = one, V = one, Qk = Qm;            /* U_1 = 1, V_1 = P = 1, Q^1 */
    for (int b = 62 - __builtin_clzll(d); b >= 0; b--) {
        U = mmul64(U, V, n, ni);                                   /* U_2k = U_k V_k */
        V = msub64(mmul64(V, V, n, ni), madd64(Qk, Qk, n), n);     /* V_2k = V_k^2 - 2 Q^k */
        Qk = mmul64(Qk, Qk, n, ni);
        if ((d >> b) & 1) {
            const u64 U2 = madd64(U, V, n);                        /* (P U + V)/2 */
            const u64 V2 = madd64(mmul64(Dm, U, n, ni), V, n);     /* (D U + P V)/2 */
            U = mhalf64(U2, n);
            V = mhalf64(V2, n);
            Qk = mmul64(Qk, Qm, n, ni);
        }
    }
    if (U == 0 || V == 0) return true;
    for (int r = 1; r < s; r++) {
        V = msub64(mmul64(V, V, n, ni), madd64(Qk, Qk, n), n);
        if (V == 0) return true;
        Qk = mmul64(Qk, Qk, n, ni);
    }
    return false;
}

/* ------------------------------------------------------------------ */
/* Fast primality, 128-bit: Montgomery with R = 2^128, n < 2^127       */
/* ------------------------------------------------------------------ */

/* a*b/R mod n for a, b < n < 2^127; ni = -n^-1 mod 2^64 */
static inline u128 mmul128(u128 a, u128 b, u128 n, u64 ni)
{
    const u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64);
    const u64 n0 = (u64)n, n1 = (u64)(n >> 64);
    const u128 p00 = (u128)a0 * b0, p01 = (u128)a0 * b1, p10 = (u128)a1 * b0, p11 = (u128)a1 * b1;
    const u64 t0 = (u64)p00;
    const u128 mid = (p00 >> 64) + (u64)p01 + (u64)p10;
    const u64 t1 = (u64)mid;
    const u128 hi = (mid >> 64) + (p01 >> 64) + (p10 >> 64) + p11;        /* limbs 2,3 */
    const u64 m = t0 * ni;                                                /* kill limb 0 */
    const u128 q0 = (u128)m * n0, q1 = (u128)m * n1;
    const u128 s1 = (u128)t1 + (q0 >> 64) + (u64)q1 + (t0 != 0);
    const u64 u1 = (u64)s1;
    const u128 hi2 = hi + (s1 >> 64) + (q1 >> 64);                         /* value = hi2:u1 */
    const u64 m2 = u1 * ni;                                               /* kill limb u1 */
    const u128 r0 = (u128)m2 * n0, r1 = (u128)m2 * n1;
    const u128 s2 = (u128)(u64)hi2 + (r0 >> 64) + (u64)r1 + (u1 != 0);
    const u128 top = (hi2 >> 64) + (s2 >> 64) + (r1 >> 64);
    u128 res = (top << 64) | (u64)s2;
    return res >= n ? res - n : res;
}

static inline u128 madd128(u128 a, u128 b, u128 n) { u128 s = a + b; return s >= n ? s - n : s; }
static inline u128 msub128(u128 a, u128 b, u128 n) { return a >= b ? a - b : a - b + n; }
static inline u128 mhalf128(u128 a, u128 n) { return (a & 1) ? (a >> 1) + (n >> 1) + 1 : a >> 1; }

static inline u128 one128(u128 n) { return (~(u128)0 % n + 1) % n; }       /* 2^128 mod n */

static bool sprp2_128(u128 n)
{
    const u64 ni = -inv64((u64)n);
    const u128 one = one128(n), mone = n - one;
    u128 d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    int b = 127;
    while (!((d >> b) & 1)) b--;
    u128 x = madd128(one, one, n);
    while (b-- > 0) {
        x = mmul128(x, x, n, ni);
        if ((d >> b) & 1) x = madd128(x, x, n);
    }
    if (x == one || x == mone) return true;
    for (int i = 1; i < s; i++) {
        x = mmul128(x, x, n, ni);
        if (x == mone) return true;
        if (x == one) return false;
    }
    return false;
}

/* four 128-bit base-2 strong tests in lock step (as sprp2_64x4) */
static void sprp2_128x4(const u128 *n, bool *res)
{
    u64 ni[4];
    u128 one[4], mone[4], d[4], x[4];
    int s[4], bmax = 0;
    for (int l = 0; l < 4; l++) {
        ni[l] = -inv64((u64)n[l]);
        one[l] = one128(n[l]);
        mone[l] = n[l] - one[l];
        d[l] = n[l] - 1;
        s[l] = 0;
        while (!(d[l] & 1)) { d[l] >>= 1; s[l]++; }
        int b = 127;
        while (!((d[l] >> b) & 1)) b--;
        if (b > bmax) bmax = b;
        x[l] = one[l];
    }
    for (int b = bmax; b >= 0; b--) {
        for (int l = 0; l < 4; l++) {
            u128 y = mmul128(x[l], x[l], n[l], ni[l]);
            if ((d[l] >> b) & 1) y = madd128(y, y, n[l]);
            x[l] = y;
        }
    }
    for (int l = 0; l < 4; l++) {
        bool r = x[l] == one[l] || x[l] == mone[l];
        for (int i = 1; i < s[l] && !r; i++) {
            x[l] = mmul128(x[l], x[l], n[l], ni[l]);
            if (x[l] == mone[l]) r = true;
            else if (x[l] == one[l]) break;
        }
        res[l] = r;
    }
}

static bool slprp128(u128 n)
{
    i64 D;
    if (!selfridge(n, &D)) return false;
    const i64 Q = (1 - D) / 4;
    const u64 ni = -inv64((u64)n);
    const u128 one = one128(n);
    const u128 Qm = mulmod_ref(one, Q >= 0 ? (u128)Q : n - (u128)(-Q), n);
    const u128 Dm = mulmod_ref(one, D >= 0 ? (u128)D : n - (u128)(-D), n);
    u128 d = n + 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    int top = 127;
    while (!((d >> top) & 1)) top--;
    u128 U = one, V = one, Qk = Qm;
    for (int b = top - 1; b >= 0; b--) {
        U = mmul128(U, V, n, ni);
        V = msub128(mmul128(V, V, n, ni), madd128(Qk, Qk, n), n);
        Qk = mmul128(Qk, Qk, n, ni);
        if ((d >> b) & 1) {
            const u128 U2 = madd128(U, V, n);
            const u128 V2 = madd128(mmul128(Dm, U, n, ni), V, n);
            U = mhalf128(U2, n);
            V = mhalf128(V2, n);
            Qk = mmul128(Qk, Qm, n, ni);
        }
    }
    if (U == 0 || V == 0) return true;
    for (int r = 1; r < s; r++) {
        V = msub128(mmul128(V, V, n, ni), madd128(Qk, Qk, n), n);
        if (V == 0) return true;
        Qk = mmul128(Qk, Qk, n, ni);
    }
    return false;
}

static bool bpsw(u128 n)                     /* n odd, > 41, no factor <= 41 */
{
    if (!(n >> 64)) return sprp2_64((u64)n) && slprp64((u64)n);
    return sprp2_128(n) && slprp128(n);
}

/* ------------------------------------------------------------------ */
/* Wheel                                                               */
/* ------------------------------------------------------------------ */

#define MAXWD 12
static const u32 SMALLP[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

static struct {
    u32 w;                  /* largest wheel prime */
    u64 M, C;               /* modulus, number of classes */
    int nd;                 /* wheel primes > 7 */
    u32 q[MAXWD], na[MAXWD];
    u8 a[MAXWD][48];        /* allowed residues mod q */
    u64 E[MAXWD];           /* CRT idempotents: E = 1 (mod q), 0 (mod M/q) */
    u64 rad[MAXWD];         /* mixed-radix weights of the class index */
    u64 r0;                 /* 97 * E_210 mod M */
} W;

static u64 invmod(u64 a, u64 m)               /* a^-1 mod m, gcd(a, m) = 1 */
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

/* ------------------------------------------------------------------ */
/* Sieving primes                                                      */
/* ------------------------------------------------------------------ */

typedef struct { u32 q, k[5]; } sprime_t;   /* bad k: phi + {0, k1..k5} (mod q) */

static struct {
    u32 B;
    int np;
    sprime_t *p;
    u32 *minv;              /* M^-1 mod q */
    u32 *bmod;              /* base mod q (per run) */
} SP;

static u32 *small_primes(u32 n, int *count)  /* primes <= n */
{
    u8 *s = calloc(n + 1, 1);
    if (!s) die("out of memory");
    int c = 0;
    for (u32 i = 2; i <= n; i++) {
        if (s[i]) continue;
        c++;
        for (u64 j = (u64)i * i; j <= n; j += i) s[j] = 1;
    }
    u32 *p = malloc(sizeof *p * (size_t)(c ? c : 1));
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

static void sieve_setup(u32 B)
{
    if (B <= W.w || B > BMAX) die("sieving bound must be in (w, 2^20]");
    free(SP.p); free(SP.minv); free(SP.bmod);
    int n;
    u32 *pr = small_primes(B, &n);
    SP.B = B;
    SP.p = malloc(sizeof *SP.p * (size_t)n);
    SP.minv = malloc(sizeof *SP.minv * (size_t)n);
    SP.bmod = malloc(sizeof *SP.bmod * (size_t)n);
    if (!SP.p || !SP.minv || !SP.bmod) die("out of memory");
    SP.np = 0;
    for (int i = 0; i < n; i++) {
        const u32 q = pr[i];
        if (q <= W.w) continue;
        const u32 mi = (u32)invmod(W.M % q, q);
        u32 c[5];
        for (int k = 1; k < 6; k++) c[k - 1] = (u32)((q - (u64)OFF[k] * mi % q) % q);
        qsort(c, 5, sizeof *c, cmp_u32);
        sprime_t *s = &SP.p[SP.np];
        s->q = q;
        memcpy(s->k, c, sizeof c);
        SP.minv[SP.np] = mi;
        SP.np++;
    }
    free(pr);
}

/* Small sieving primes are applied as precomputed word patterns: T_q[s] has bit i set iff
   (s + i) mod q is one of the six bad offsets {0, k1..k5}, and word j of a segment whose
   phase is phi gets T_q[(64 j - phi) mod q].  Groups of four primes are ORed in one pass. */
typedef struct { const u64 *t; u32 q, step; } pat_t;

static struct {
    u32 qmax;               /* primes q < qmax use patterns */
    int n;                  /* they are SP.p[0 .. n) */
    pat_t *p;
    u64 *tab;
} PT;

static void pattern_setup(u32 qmax)
{
    free(PT.p); free(PT.tab);
    PT.qmax = qmax;
    PT.n = 0;
    size_t words = 0;
    while (PT.n < SP.np && SP.p[PT.n].q < qmax) words += SP.p[PT.n++].q;
    PT.p = malloc(sizeof *PT.p * (size_t)(PT.n ? PT.n : 1));
    PT.tab = malloc(sizeof *PT.tab * (words ? words : 1));
    if (!PT.p || !PT.tab) die("out of memory");
    u64 *t = PT.tab;
    for (int i = 0; i < PT.n; i++) {
        const sprime_t *sp = &SP.p[i];
        const u32 q = sp->q;
        u8 bad[1024] = {0};
        bad[0] = 1;
        for (int j = 0; j < 5; j++) bad[sp->k[j]] = 1;
        for (u32 s = 0; s < q; s++) {
            u64 w = 0;
            for (u32 b = 0; b < 64; b++) if (bad[(s + b) % q]) w |= 1ULL << b;
            t[s] = w;
        }
        PT.p[i].t = t;
        PT.p[i].q = q;
        PT.p[i].step = 64 % q;
        t += q;
    }
}

#define ADV(s, st, q) do { s += (st); s = s >= (q) ? s - (q) : s; } while (0)

/* Write the patterns of the small primes into bm[0 .. nw) (replaces clearing it). */
static void presieve(u64 *restrict bm, u32 nw, const u32 *phi)
{
    int i = 0;
    bool first = true;
    for (; i + 4 <= PT.n; i += 4) {
        const pat_t *a = &PT.p[i], *b = &PT.p[i + 1], *c = &PT.p[i + 2], *d = &PT.p[i + 3];
        const u64 *ta = a->t, *tb = b->t, *tc = c->t, *td = d->t;
        const u32 qa = a->q, qb = b->q, qc = c->q, qd = d->q;
        const u32 da = a->step, db = b->step, dc = c->step, dd = d->step;
        u32 sa = phi[i] ? qa - phi[i] : 0, sb = phi[i + 1] ? qb - phi[i + 1] : 0;
        u32 sc = phi[i + 2] ? qc - phi[i + 2] : 0, sd = phi[i + 3] ? qd - phi[i + 3] : 0;
        if (first) {
            for (u32 j = 0; j < nw; j++) {
                bm[j] = ta[sa] | tb[sb] | tc[sc] | td[sd];
                ADV(sa, da, qa); ADV(sb, db, qb); ADV(sc, dc, qc); ADV(sd, dd, qd);
            }
            first = false;
        } else {
            for (u32 j = 0; j < nw; j++) {
                bm[j] |= ta[sa] | tb[sb] | tc[sc] | td[sd];
                ADV(sa, da, qa); ADV(sb, db, qb); ADV(sc, dc, qc); ADV(sd, dd, qd);
            }
        }
    }
    for (; i < PT.n; i++) {
        const u64 *ta = PT.p[i].t;
        const u32 qa = PT.p[i].q, da = PT.p[i].step;
        u32 sa = phi[i] ? qa - phi[i] : 0;
        if (first) {
            for (u32 j = 0; j < nw; j++) { bm[j] = ta[sa]; ADV(sa, da, qa); }
            first = false;
        } else {
            for (u32 j = 0; j < nw; j++) { bm[j] |= ta[sa]; ADV(sa, da, qa); }
        }
    }
    if (first) memset(bm, 0, (size_t)nw * 8);
}

/* ------------------------------------------------------------------ */
/* Run configuration                                                   */
/* ------------------------------------------------------------------ */

static struct {
    u128 lo, hi;            /* initial members in [lo, hi) */
    u128 base;              /* max(lo, P0): start of the sieved part */
    u64 bmM;                /* base mod M */
    int nbins;
    u128 bnd[MAXBIN + 1];   /* bin b = [bnd[b], bnd[b+1]) */
    u32 S;                  /* segment bits */
    bool notest;            /* bench -x: sieve only */
} R;

static void run_setup(u128 lo, u128 hi)
{
    R.lo = lo;
    R.hi = hi;
    R.base = lo > P0 ? lo : P0;
    R.bmM = (u64)(R.base % W.M);
    for (int i = 0; i < SP.np; i++) SP.bmod[i] = (u32)(R.base % SP.p[i].q);
}

/* ------------------------------------------------------------------ */
/* Per-class sieve                                                     */
/* ------------------------------------------------------------------ */

typedef struct {
    u64 cnt[MAXBIN], cks[MAXBIN];
    u64 surv, cand, spsp, bits;
} acc_t;

typedef struct {
    u64 *bm;                /* segment bitmap, 1 = composite */
    u32 *phi;               /* phase per sieving prime */
    u64 *sv;                /* survivors of a segment (initial members) */
    u128 *sv2;              /* the same for segments reaching 2^64 */
    u32 svcap;
} tctx_t;

static pthread_mutex_t g_mu = PTHREAD_MUTEX_INITIALIZER;
static bool g_print = true;

static struct { u128 n, p; } g_spsp[MAXSPSP];
static int g_nspsp = 0;

static u128 *g_list = NULL;                   /* list mode */
static size_t g_nlist = 0, g_caplist = 0;

static void record_spsp(u128 n, u128 p)
{
    pthread_mutex_lock(&g_mu);
    if (g_nspsp < MAXSPSP) { g_spsp[g_nspsp].n = n; g_spsp[g_nspsp].p = p; g_nspsp++; }
    if (g_print) {
        if (stderr_tty) fputs("\r\033[K", stderr);
        printf("SPSP %s (base-2 strong pseudoprime in candidate p = %s)\n", u128s(n), u128s(p));
        fflush(stdout);
    }
    pthread_mutex_unlock(&g_mu);
}

static void found(u128 p, acc_t *A)
{
    int b = 0;
    while (b + 1 < R.nbins && p >= R.bnd[b + 1]) b++;
    A->cnt[b]++;
    A->cks[b] += hashp(p);
    if (g_caplist) {
        pthread_mutex_lock(&g_mu);
        if (g_nlist < g_caplist) g_list[g_nlist++] = p;
        pthread_mutex_unlock(&g_mu);
    }
}


/* Survivors p[0..n): keep those whose member p+OFF[m] is a base-2 strong probable prime,
   for m = 0..5 in turn (4 lanes at a time); confirm the rest with test64. */
static void test_batch64(tctx_t *T, u64 *p, u32 n, acc_t *A)
{
    (void)T;
    A->surv += n;
    for (int m = 0; m < 6 && n; m++) {
        u32 k = 0, i = 0;
        for (; i + 4 <= n; i += 4) {
            u64 v[4];
            bool r[4];
            for (int l = 0; l < 4; l++) v[l] = p[i + l] + OFF[m];
            sprp2_64x4(v, r);
            for (int l = 0; l < 4; l++) if (r[l]) p[k++] = p[i + l];
        }
        for (; i < n; i++) if (sprp2_64(p[i] + OFF[m])) p[k++] = p[i];
        n = k;
    }
    for (u32 i = 0; i < n; i++) {
        A->cand++;
        bool ok = true;
        for (int m = 0; m < 6 && ok; m++)
            if (!slprp64(p[i] + OFF[m])) { A->spsp++; record_spsp(p[i] + OFF[m], p[i]); ok = false; }
        if (ok) found(p[i], A);
    }
}

/* as test_batch64, for segments that reach 2^64 (two-limb arithmetic for members >= 2^64) */
static void test_batch128(u128 *p, u32 n, acc_t *A)
{
    A->surv += n;
    for (int m = 0; m < 6 && n; m++) {
        u32 k = 0, i = 0;
        for (; i + 4 <= n; i += 4) {
            u128 v[4];
            bool r[4];
            for (int l = 0; l < 4; l++) v[l] = p[i + l] + OFF[m];
            sprp2_128x4(v, r);
            for (int l = 0; l < 4; l++) if (r[l]) p[k++] = p[i + l];
        }
        for (; i < n; i++) {
            const u128 v = p[i] + OFF[m];
            if ((v >> 64) ? sprp2_128(v) : sprp2_64((u64)v)) p[k++] = p[i];
        }
        n = k;
    }
    for (u32 i = 0; i < n; i++) {
        A->cand++;
        bool ok = true;
        for (int m = 0; m < 6 && ok; m++) {
            const u128 v = p[i] + OFF[m];
            if (!((v >> 64) ? slprp128(v) : slprp64((u64)v))) { A->spsp++; record_spsp(v, p[i]); ok = false; }
        }
        if (ok) found(p[i], A);
    }
}

#define SETBIT(bm, x) ((bm)[(x) >> 6] |= 1ULL << ((x) & 63))

/* Mark the six progressions of one prime in bits [0, len) and advance its phase by len. */
static inline void mark_prime(u64 *restrict bm, u32 len, const sprime_t *sp, u32 *phi_io)
{
    const u32 q = sp->q, k1 = sp->k[0], k2 = sp->k[1], k3 = sp->k[2], k4 = sp->k[3], k5 = sp->k[4];
    const u32 phi = *phi_io;
    const u32 t = q - phi;                    /* round phi - q: positions k_j - t for k_j >= t */
    if (k5 >= t) {
        if (k1 >= t && k1 - t < len) SETBIT(bm, k1 - t);
        if (k2 >= t && k2 - t < len) SETBIT(bm, k2 - t);
        if (k3 >= t && k3 - t < len) SETBIT(bm, k3 - t);
        if (k4 >= t && k4 - t < len) SETBIT(bm, k4 - t);
        if (k5 - t < len) SETBIT(bm, k5 - t);
    }
    u32 b = phi;
    if (len > k5) {
        const u32 lim = len - k5;
        for (; b < lim; b += q) {
            SETBIT(bm, b);
            SETBIT(bm, b + k1);
            SETBIT(bm, b + k2);
            SETBIT(bm, b + k3);
            SETBIT(bm, b + k4);
            SETBIT(bm, b + k5);
        }
    }
    if (b < len) {                            /* partial round; offsets are sorted */
        SETBIT(bm, b);
        if (b + k1 < len) { SETBIT(bm, b + k1);
        if (b + k2 < len) { SETBIT(bm, b + k2);
        if (b + k3 < len) { SETBIT(bm, b + k3);
        if (b + k4 < len) { SETBIT(bm, b + k4); } } } }
    }
    *phi_io = b >= len ? b - len : b + q - len;
}

static void process_class(tctx_t *T, u64 c, acc_t *A)
{
    const u64 M = W.M;
    const u64 r = class_r(c);
    const u64 delta = r >= R.bmM ? r - R.bmM : r + (M - R.bmM);
    const u128 ps = R.base + delta;           /* first candidate >= base in this class */
    if (ps >= R.hi) return;
    const u128 span = (R.hi - 1 - ps) / M;
    if (span >> 62) die("class too long; use a larger wheel (-w)");
    const u64 L = (u64)span + 1;
    A->bits += L;
    for (int i = 0; i < SP.np; i++) {
        const u32 q = SP.p[i].q;
        const u32 u = (u32)((SP.bmod[i] + delta % q) % q);
        T->phi[i] = u ? (u32)((u64)(q - u) * SP.minv[i] % q) : 0;
    }
    const bool fits64 = ps + (u128)M * (L - 1) + 16 < ((u128)1 << 64);
    const int np = SP.np;
    const sprime_t *sp = SP.p;
    const u64 nseg = (L + R.S - 1) / R.S;
    const u32 seg = (u32)(((L + nseg - 1) / nseg + 63) & ~(u64)63);
    for (u64 k0 = 0; k0 < L; k0 += seg) {
        const u32 len = L - k0 < seg ? (u32)(L - k0) : seg;
        const u32 nw = (len + 63) >> 6;
        u64 *bm = T->bm;
        presieve(bm, nw, T->phi);
        for (int i = 0; i < PT.n; i++) {          /* advance the pattern phases by len */
            const u32 q = sp[i].q, m = len % q, f = T->phi[i];
            T->phi[i] = f >= m ? f - m : f + q - m;
        }
        for (int i = PT.n; i < np; i++) mark_prime(bm, len, &sp[i], &T->phi[i]);
        if (len & 63) bm[nw - 1] |= ~0ULL << (len & 63);
        if (R.notest) {
            for (u32 j = 0; j < nw; j++) A->surv += (u64)__builtin_popcountll(~bm[j]);
            continue;
        }
        const u128 pk = ps + (u128)M * k0;
        if (fits64 || pk + (u128)M * (len - 1) + 16 < ((u128)1 << 64)) {
            const u64 p0 = (u64)pk;
            u64 *sv = T->sv;
            u32 ns = 0;
            for (u32 j = 0; j < nw; j++) {
                u64 x = ~bm[j];
                while (x) {
                    const u32 bit = (u32)__builtin_ctzll(x);
                    x &= x - 1;
                    if (ns == T->svcap) { test_batch64(T, sv, ns, A); ns = 0; }
                    sv[ns++] = p0 + M * ((u64)j * 64 + bit);
                }
            }
            test_batch64(T, sv, ns, A);
        } else {
            u128 *sv = T->sv2;
            u32 ns = 0;
            for (u32 j = 0; j < nw; j++) {
                u64 x = ~bm[j];
                while (x) {
                    const u32 bit = (u32)__builtin_ctzll(x);
                    x &= x - 1;
                    if (ns == T->svcap) { test_batch128(sv, ns, A); ns = 0; }
                    sv[ns++] = pk + (u128)M * ((u64)j * 64 + bit);
                }
            }
            test_batch128(sv, ns, A);
        }
    }
}

/* initial members in [lo, min(hi, P0)), tested directly (includes p = 7) */
static void small_region(acc_t *A)
{
    if (R.lo >= P0) return;
    const u64 a = (u64)R.lo, b = R.hi < P0 ? (u64)R.hi : P0;
    for (u64 p = a; p < b; p++) {
        if (p % 30 != 7) continue;            /* 7 and 97 + 210k are 7 mod 30 */
        if (is_sext_ref(p)) found(p, A);
    }
}

/* ------------------------------------------------------------------ */
/* Chunked parallel driver                                             */
/* ------------------------------------------------------------------ */

typedef struct {
    bool done;
    acc_t a;
} slot_t;

static struct {
    u64 chunk, nchunks;     /* classes per chunk, number of chunks */
    atomic_ullong next;
    u64 frontier;           /* chunks [0, frontier) are complete */
    acc_t tot;              /* totals over complete chunks */
    slot_t slot[RINGC];
    atomic_ullong bits_live;
    u64 c0, c1;             /* class range [c0, c1) of this run */
    u64 *sample;            /* bench: explicit class list (NULL = the class range) */
    u64 nsample;
} D;

static FILE *g_log = NULL;

static void log_chunk(u64 idx, const acc_t *a)
{
    /* absolute chunk index when the class range is chunk-aligned (-P parts), else relative */
    const u64 aidx = D.c0 % D.chunk == 0 ? D.c0 / D.chunk + idx : idx;
    fprintf(g_log, "C %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64, aidx, a->surv, a->cand, a->spsp);
    for (int b = 0; b < R.nbins; b++) fprintf(g_log, " %" PRIu64 ":%016" PRIx64, a->cnt[b], a->cks[b]);
    fputc('\n', g_log);
    fflush(g_log);
}

static void acc_add(acc_t *t, const acc_t *a)
{
    for (int b = 0; b < MAXBIN; b++) { t->cnt[b] += a->cnt[b]; t->cks[b] += a->cks[b]; }
    t->surv += a->surv; t->cand += a->cand; t->spsp += a->spsp; t->bits += a->bits;
}

static void process_chunk(tctx_t *T, u64 idx)
{
    acc_t a;
    memset(&a, 0, sizeof a);
    if (idx == 0 && !D.sample && D.c0 == 0) small_region(&a);
    const u64 c0 = idx * D.chunk;
    const u64 n = D.sample ? D.nsample : D.c1 - D.c0;
    const u64 c1 = c0 + D.chunk < n ? c0 + D.chunk : n;
    for (u64 c = c0; c < c1; c++) {
        const u64 bits0 = a.bits;
        process_class(T, D.sample ? D.sample[c] : D.c0 + c, &a);
        atomic_fetch_add_explicit(&D.bits_live, a.bits - bits0, memory_order_relaxed);
        if (g_stop) return;                   /* chunk abandoned; redone on resume */
    }
    pthread_mutex_lock(&g_mu);
    slot_t *s = &D.slot[idx % RINGC];
    s->done = true;
    s->a = a;
    while (D.frontier < D.nchunks && D.slot[D.frontier % RINGC].done) {
        slot_t *f = &D.slot[D.frontier % RINGC];
        f->done = false;
        acc_add(&D.tot, &f->a);
        if (g_log) log_chunk(D.frontier, &f->a);
        D.frontier++;
    }
    pthread_mutex_unlock(&g_mu);
}

static void *worker(void *arg)
{
    (void)arg;
    tctx_t T;
    T.bm = aligned_alloc(64, (((size_t)R.S / 64 + 8) * 8 + 63) & ~(size_t)63);
    T.phi = malloc(sizeof *T.phi * (size_t)(SP.np ? SP.np : 1));
    T.svcap = 4096;
    T.sv = malloc(sizeof *T.sv * T.svcap);
    T.sv2 = malloc(sizeof *T.sv2 * T.svcap);
    if (!T.bm || !T.phi || !T.sv || !T.sv2) die("out of memory");
    while (!g_stop) {
        const u64 idx = atomic_fetch_add(&D.next, 1);
        if (idx >= D.nchunks) break;
        for (;;) {                            /* stay within RINGC of the frontier */
            pthread_mutex_lock(&g_mu);
            const u64 f = D.frontier;
            pthread_mutex_unlock(&g_mu);
            if (idx < f + RINGC - 1 || g_stop) break;
            usleep(2000);
        }
        if (g_stop) break;
        process_chunk(&T, idx);
    }
    free(T.bm);
    free(T.phi);
    free(T.sv);
    free(T.sv2);
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
    fprintf(f, "a350826 checkpoint 1\nlo %s\nhi %s\nwheel %u\nbound %u\nchunk %" PRIu64 "\nclasses %" PRIu64
               " %" PRIu64 "\nbins %d\n", u128s(R.lo), u128s(R.hi), W.w, SP.B, D.chunk, D.c0, D.c1, R.nbins);
    for (int b = 0; b <= R.nbins; b++) fprintf(f, "bnd %s\n", u128s(R.bnd[b]));
    fprintf(f, "frontier %" PRIu64 "\nsurv %" PRIu64 "\ncand %" PRIu64 "\nspsp %" PRIu64 "\nbits %" PRIu64 "\n",
            D.frontier, D.tot.surv, D.tot.cand, D.tot.spsp, D.tot.bits);
    for (int b = 0; b < R.nbins; b++) fprintf(f, "bin %d %" PRIu64 " %016" PRIx64 "\n", b, D.tot.cnt[b], D.tot.cks[b]);
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
        u64 c, k;
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
        if (sscanf(line, "bin %d %" SCNu64 " %" SCNx64, &b, &c, &k) == 3) {
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

static void on_sigint(int sig)
{
    (void)sig;
    g_stop = 1;
}

typedef struct { double seconds; u64 bits; } run_stats_t;

static run_stats_t run_count(int threads, const char *state, int interval, bool quiet)
{
    D.nchunks = ((D.sample ? D.nsample : D.c1 - D.c0) + D.chunk - 1) / D.chunk;
    D.frontier = 0;
    memset(&D.tot, 0, sizeof D.tot);
    memset(D.slot, 0, sizeof D.slot);
    atomic_store(&D.bits_live, 0);
    g_nspsp = 0;
    if (state && load_state(state) && !quiet)
        fprintf(stderr, "resuming from %s: %" PRIu64 " of %" PRIu64 " chunks done\n", state, D.frontier, D.nchunks);
    const u64 frontier0 = D.frontier, bits0 = D.tot.bits;
    atomic_store(&D.next, D.frontier);

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_sigint;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    pthread_t *th = calloc((size_t)threads, sizeof *th);
    if (!th) die("out of memory");
    for (int i = 0; i < threads; i++)
        if (pthread_create(&th[i], NULL, worker, NULL) != 0) die("pthread_create failed");

    enum { NS = 128 };      /* rate window: one sample per second, ~2 minutes deep */
    double ts[NS];
    u64 bs[NS], fs[NS];
    int ns = 1;
    const double t0 = now();
    double last_status = 0, last_ckpt = t0, last_sample = t0;
    ts[0] = t0; bs[0] = 0; fs[0] = frontier0;
    char b1[32], b2[32];
    for (;;) {
        usleep(100000);
        const double t = now();
        pthread_mutex_lock(&g_mu);
        const u64 f = D.frontier;
        u64 found = 0;
        for (int b = 0; b < R.nbins; b++) found += D.tot.cnt[b];
        const u64 spsp = D.tot.spsp;
        pthread_mutex_unlock(&g_mu);
        if (f >= D.nchunks || g_stop) break;
        const u64 bl = atomic_load(&D.bits_live);
        if (t - last_sample >= 1.0) {
            last_sample = t;
            if (ns == NS) {
                memmove(ts, ts + 1, (NS - 1) * sizeof *ts);
                memmove(bs, bs + 1, (NS - 1) * sizeof *bs);
                memmove(fs, fs + 1, (NS - 1) * sizeof *fs);
                ns--;
            }
            ts[ns] = t; bs[ns] = bl; fs[ns] = f; ns++;
        }
        const double dt = t - ts[0];
        const double brate = dt > 0.5 ? (double)(bl - bs[0]) / dt : 0;
        const double crate = dt > 0.5 ? (double)(f - fs[0]) / dt : 0;
        const double eta = crate > 0 ? (double)(D.nchunks - f) / crate : -1;
        const bool show = stderr_tty ? (t - last_status >= 1.0) : (t - last_status >= interval);
        if (show && !quiet) {
            last_status = t;
            fprintf(stderr, "%s%s  chunks %" PRIu64 "/%" PRIu64 " (%.3f%%)  %.3g bits/s  found %" PRIu64 "%s  ETA %s",
                    stderr_tty ? "\r\033[K" : "progress: ", fmt_dur(t - t0, b1, sizeof b1), f, D.nchunks,
                    100.0 * (double)f / (double)D.nchunks, brate, found, spsp ? "  SPSP!" : "",
                    fmt_dur(eta, b2, sizeof b2));
            if (!stderr_tty) fputc('\n', stderr);
            fflush(stderr);
        }
        if (state && t - last_ckpt >= interval) {
            last_ckpt = t;
            save_state(state);
        }
    }
    for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
    const double t1 = now();
    if (stderr_tty && !quiet) fputs("\r\033[K", stderr);
    if (state) save_state(state);
    if (g_stop && D.frontier < D.nchunks && !quiet)
        fprintf(stderr, "stopped after %" PRIu64 " of %" PRIu64 " chunks%s\n", D.frontier, D.nchunks,
                state ? ", checkpoint saved" : "");
    free(th);
    run_stats_t st = {t1 - t0, D.tot.bits - bits0};
    return st;
}

/* ------------------------------------------------------------------ */
/* Parameter choice                                                    */
/* ------------------------------------------------------------------ */

static double rho_wheel(u32 w)                /* fraction of integers left by the wheel */
{
    double r = 1.0 / 210;
    for (int i = 4; i < 15 && SMALLP[i] <= w; i++) r *= (double)(SMALLP[i] - 6) / SMALLP[i];
    return r;
}

/* Cost model (cycles), fitted to M1 Pro runs with B = 2^16: per class ~100 per sieving prime
   (phases, first segment), per further segment ~12 per prime, per candidate bit ~1.2 per mark
   plus ~1.3 (patterns, scan, tests).  Picks the cheapest wheel for the range (w = 31 for
   [1e17, 1e18), 37 for [1e18, 1e20)). */
static u32 auto_wheel(u128 lo, u128 hi, u32 B, u32 S)
{
    const u128 base = lo > P0 ? lo : P0;
    if (hi <= base) return 7;                 /* nothing to sieve */
    int np;
    u32 *pr = small_primes(B, &np);
    const double len = (double)(hi - base);
    u32 best = 7;
    double bestc = 1e300;
    for (int i = 3; i < 15; i++) {
        const u32 w = SMALLP[i];
        double M = 1, C = 1, marks = 0;
        int nsp = 0;
        for (int j = 0; j < 15 && SMALLP[j] <= w; j++) M *= SMALLP[j];
        for (int j = 4; j < 15 && SMALLP[j] <= w; j++) C *= SMALLP[j] - 6;
        for (int j = 0; j < np; j++) if (pr[j] > w) { marks += 6.0 / pr[j]; nsp++; }
        const double L = len / M;
        if (L < 1 && w > 7) break;
        const double segs = ceil(L / S > 0 ? L / S : 1);
        const double cost = C * (nsp * 100.0 + (segs - 1) * nsp * 12.0 + 2000) + len * rho_wheel(w) * (1.2 * marks + 1.3);
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

/* ------------------------------------------------------------------ */
/* Commands                                                            */
/* ------------------------------------------------------------------ */

static void usage(void)
{
    fputs("usage: a350826 count LO HI [-b B1,B2,..] [-w W] [-B B] [-s S] [-c CH] [-C C0:C1 | -P I/N] [-t T]\n"
          "                     [-S FILE] [-i SECS] [-L FILE] [-q]\n"
          "       a350826 list LO HI\n"
          "       a350826 check P\n"
          "       a350826 bench LO HI [-n N] [-w W] [-B B] [-s S] [-t T] [-x]\n"
          "       a350826 recheck LO HI CHUNKLOG -w W -c CH [-B B] [-b ..] [-n K] [-t T]\n"
          "       a350826 selftest [N] [-t T]\n", stderr);
    exit(2);
}

typedef struct {
    u128 lo, hi;
    u32 w, B, S, pmax;
    u64 chunk, nsample;
    u64 cr0, cr1;           /* -C class range (cr1 = 0: all) */
    u64 part, nparts;       /* -P part/nparts */
    int threads, interval;
    const char *state, *log;
    bool quiet, notest;
    int nb;
    u128 bnd[MAXBIN];
} opts_t;

static void parse_bins(const char *s, opts_t *o)
{
    char buf[1024];
    snprintf(buf, sizeof buf, "%s", s);
    for (char *tok = strtok(buf, ","); tok; tok = strtok(NULL, ",")) {
        if (o->nb >= MAXBIN - 1) die("too many bin boundaries");
        o->bnd[o->nb++] = parse_num(tok);
    }
}

static void parse_opts(int argc, char **argv, int first, opts_t *o, int npos, u128 *pos)
{
    memset(o, 0, sizeof *o);
    o->B = 1u << 16;
    o->S = 1u << 19;
    o->pmax = 256;
    o->threads = default_threads();
    o->interval = 60;
    o->nsample = 2000;
    int np = 0;
    for (int i = first; i < argc; i++) {
        const char *a = argv[i];
        if (a[0] == '-' && a[1] && !isdigit((unsigned char)a[1]) && a[2] == 0) {
            const char c = a[1];
            if (c == 'q') { o->quiet = true; continue; }
            if (c == 'x') { o->notest = true; continue; }
            if (i + 1 >= argc) usage();
            const char *v = argv[++i];
            switch (c) {
            case 'b': parse_bins(v, o); break;
            case 'w': o->w = (u32)parse_u64(v, 7, 47, "wheel prime"); break;
            case 'B': o->B = (u32)parse_u64(v, 11, BMAX, "sieving bound"); break;
            case 'p': o->pmax = (u32)parse_u64(v, 0, 1024, "pattern bound"); break;
            case 's': o->S = (u32)parse_u64(v, 64, 1u << 26, "segment bits"); break;
            case 'c': o->chunk = parse_u64(v, 1, 1ULL << 40, "chunk"); break;
            case 'n': o->nsample = parse_u64(v, 1, 1ULL << 40, "sample"); break;
            case 't': o->threads = (int)parse_u64(v, 1, 1024, "threads"); break;
            case 'i': o->interval = (int)parse_u64(v, 1, 86400, "interval"); break;
            case 'S': o->state = v; break;
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
            case 'L': o->log = v; break;
            default: usage();
            }
            continue;
        }
        if (np >= npos) usage();
        pos[np++] = parse_num(a);
    }
    if (np != npos) usage();
    o->S = (o->S + 63) & ~63u;
}

/* Configure wheel, primes, run and bins for [lo, hi). */
static void configure(opts_t *o, u128 lo, u128 hi)
{
    if (hi <= lo) die("HI must be larger than LO");
    if (hi > HIMAX) die("HI must be <= 2^80");
    o->lo = lo; o->hi = hi;
    const u32 w = o->w ? o->w : auto_wheel(lo, hi, o->B, o->S);
    wheel_setup(w);
    if (o->B <= W.w) die("sieving bound must exceed the wheel prime");
    sieve_setup(o->B);
    pattern_setup(o->pmax);
    R.S = o->S;
    R.notest = o->notest;
    run_setup(lo, hi);
    R.nbins = 0;
    R.bnd[R.nbins++] = lo;
    for (int i = 0; i < o->nb; i++) {
        if (o->bnd[i] <= R.bnd[R.nbins - 1] || o->bnd[i] >= hi) die("bin boundaries must increase inside (LO, HI)");
        R.bnd[R.nbins++] = o->bnd[i];
    }
    R.bnd[R.nbins] = hi;
    D.chunk = o->chunk ? o->chunk : auto_chunk(W.C);
    D.sample = NULL;
    D.c0 = 0;
    D.c1 = W.C;
    if (o->cr1) {
        if (o->cr1 > W.C) die("-C range beyond the %" PRIu64 " classes of wheel %u", W.C, W.w);
        D.c0 = o->cr0;
        D.c1 = o->cr1;
    } else if (o->nparts) {                   /* chunk-aligned parts */
        const u64 nch = (W.C + D.chunk - 1) / D.chunk;
        const u64 a = (u64)((u128)nch * o->part / o->nparts), b = (u64)((u128)nch * (o->part + 1) / o->nparts);
        D.c0 = a * D.chunk;
        D.c1 = b * D.chunk < W.C ? b * D.chunk : W.C;
        if (D.c1 <= D.c0) die("part %" PRIu64 "/%" PRIu64 " is empty", o->part, o->nparts);
    }
}

/* "A350826(n)" if [a, b) = [10^(n-1), 10^n) */
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

static int cmd_count(int argc, char **argv)
{
    opts_t o;
    u128 pos[2];
    parse_opts(argc, argv, 2, &o, 2, pos);
    configure(&o, pos[0], pos[1]);
    if (o.log) {
        g_log = fopen(o.log, "a");
        if (!g_log) die("cannot open %s: %s", o.log, strerror(errno));
    }
    if (!o.quiet)
        fprintf(stderr, "counting [%s, %s): wheel %u (M = %" PRIu64 ", %" PRIu64 " classes, ~%.3g candidates), "
                        "B = %u (%d primes), S = %u, %" PRIu64 " chunks of %" PRIu64 ", %d threads\n",
                u128s(o.lo), u128s(o.hi), W.w, W.M, W.C, (double)(o.hi - R.base) * rho_wheel(W.w), SP.B, SP.np,
                R.S, (D.c1 - D.c0 + D.chunk - 1) / D.chunk, D.chunk, o.threads);
    if (!o.quiet && (D.c0 != 0 || D.c1 != W.C))
        fprintf(stderr, "classes [%" PRIu64 ", %" PRIu64 ") of %" PRIu64 " only (partial count)\n", D.c0, D.c1, W.C);
    run_stats_t st = run_count(o.threads, o.state, o.interval, o.quiet);
    if (g_log) fclose(g_log);
    if (D.frontier < D.nchunks) { print_results("PARTIAL"); return 1; }
    print_results(D.c0 == 0 && D.c1 == W.C ? "RESULT" : "PART");
    char b1[32];
    if (!o.quiet)
        fprintf(stderr, "done in %s (%.3g bits/s)\n", fmt_dur(st.seconds, b1, sizeof b1),
                st.seconds > 0 ? (double)st.bits / st.seconds : 0);
    return 0;
}

static int cmp_u128(const void *a, const void *b)
{
    u128 x = *(const u128 *)a, y = *(const u128 *)b;
    return x < y ? -1 : x > y;
}

static int cmd_list(int argc, char **argv)
{
    opts_t o;
    u128 pos[2];
    parse_opts(argc, argv, 2, &o, 2, pos);
    configure(&o, pos[0], pos[1]);
    g_caplist = 1u << 24;
    g_list = malloc(sizeof *g_list * g_caplist);
    if (!g_list) die("out of memory");
    run_count(o.threads, NULL, o.interval, true);
    if (g_nlist == g_caplist) die("too many sextuplets for list mode");
    qsort(g_list, g_nlist, sizeof *g_list, cmp_u128);
    for (size_t i = 0; i < g_nlist; i++) printf("%s\n", u128s(g_list[i]));
    free(g_list);
    g_list = NULL; g_caplist = 0; g_nlist = 0;
    return 0;
}

static int cmd_check(int argc, char **argv)
{
    if (argc != 3) usage();
    const u128 p = parse_num(argv[2]);
    if (p >= HIMAX) die("P must be < 2^80");
    printf("p = %s\n", u128s(p));
    bool all = true;
    for (int i = 0; i < 6; i++) {
        const u128 n = p + OFF[i];
        const bool ref = is_prime_ref(n);
        bool fast = ref;
        if (n > 41 && (n & 1) && n % 3 && n % 5 && n % 7) fast = bpsw(n);
        printf("  p+%-2u = %s  %s%s\n", OFF[i], u128s(n), ref ? "prime" : "composite",
               fast == ref ? "" : "  (BPSW DISAGREES)");
        all &= ref;
    }
    printf("%s a prime sextuplet\n", all ? "IS" : "is NOT");
    return all ? 0 : 1;
}

static int cmd_bench(int argc, char **argv)
{
    opts_t o;
    u128 pos[2];
    parse_opts(argc, argv, 2, &o, 2, pos);
    configure(&o, pos[0], pos[1]);
    const u64 n = o.nsample < W.C ? o.nsample : W.C;
    u64 *smp = malloc(sizeof *smp * n);
    if (!smp) die("out of memory");
    for (u64 i = 0; i < n; i++) smp[i] = (u64)((u128)i * W.C / n);
    D.sample = smp;
    D.nsample = n;
    D.chunk = 1;
    g_print = false;
    fprintf(stderr, "bench [%s, %s): wheel %u (%" PRIu64 " classes), B = %u (%d primes), S = %u, %d threads, "
                    "%" PRIu64 " sample classes%s\n", u128s(o.lo), u128s(o.hi), W.w, W.C, SP.B, SP.np, R.S, o.threads, n,
            o.notest ? ", sieve only" : "");
    run_stats_t st = run_count(o.threads, NULL, o.interval, true);
    const double scale = (double)W.C / (double)n;
    const double total = st.seconds * scale;
    char b1[32];
    printf("BENCH w %u B %u S %u threads %d: %.3f s for %" PRIu64 " classes, %.4g bits/s, %.4g bits/s/thread\n",
           W.w, SP.B, R.S, o.threads, st.seconds, n, (double)st.bits / st.seconds,
           (double)st.bits / st.seconds / o.threads);
    printf("BENCH per class: %.4g bits, %.4g survivors, %.4g candidates; survivor rate %.3e\n",
           (double)D.tot.bits / n, (double)D.tot.surv / n, (double)D.tot.cand / n,
           (double)D.tot.surv / (double)D.tot.bits);
    u64 found = 0;
    for (int b = 0; b < R.nbins; b++) found += D.tot.cnt[b];
    printf("BENCH projection for all %" PRIu64 " classes: %s with %d threads (~%.4g sextuplets)\n", W.C,
           fmt_dur(total, b1, sizeof b1), o.threads, (double)found * scale);
    free(smp);
    D.sample = NULL;
    return 0;
}

/* Recompute chunks of a chunk log (from either tool) and compare the lines. */
static int cmd_recheck(int argc, char **argv)
{
    if (argc < 5) usage();
    const char *logpath = argv[4];
    opts_t o;
    u128 pos[2];
    char *args[256];
    int na = 0;
    for (int i = 0; i < argc && na < 255; i++) if (i != 4) args[na++] = argv[i];
    parse_opts(na, args, 2, &o, 2, pos);
    if (!o.w || !o.chunk) die("recheck needs the original -w and -c (and -B, -b if not default)");
    FILE *f = fopen(logpath, "r");
    if (!f) die("cannot read %s: %s", logpath, strerror(errno));
    size_t cap = 1024, n = 0;
    char **lines = malloc(sizeof *lines * cap);
    u64 *idx = malloc(sizeof *idx * cap);
    char buf[4096];
    while (fgets(buf, sizeof buf, f)) {
        unsigned long long k;
        if (sscanf(buf, "C %llu ", &k) != 1) continue;
        if (n == cap) { cap *= 2; lines = realloc(lines, sizeof *lines * cap); idx = realloc(idx, sizeof *idx * cap); }
        buf[strcspn(buf, "\n")] = 0;
        lines[n] = strdup(buf);
        idx[n] = k;
        n++;
    }
    fclose(f);
    if (!n) die("no chunk lines in %s", logpath);
    const u64 want = o.nsample && o.nsample < n ? o.nsample : n;
    /* pick `want` distinct lines at random (partial Fisher-Yates) */
    for (u64 i = 0; i < want; i++) {
        const u64 j = i + rnd_r() % (n - i);
        char *t = lines[i]; lines[i] = lines[j]; lines[j] = t;
        const u64 x = idx[i]; idx[i] = idx[j]; idx[j] = x;
    }
    g_print = false;
    int bad = 0;
    const double t0 = now();
    for (u64 i = 0; i < want; i++) {
        opts_t c = o;
        configure(&c, pos[0], pos[1]);
        const u64 a = idx[i] * D.chunk, b = a + D.chunk < W.C ? a + D.chunk : W.C;
        if (a >= W.C) die("chunk %" PRIu64 " is beyond the %" PRIu64 " classes", idx[i], W.C);
        D.c0 = a;
        D.c1 = b;
        run_count(o.threads, NULL, 60, true);
        char mine[4096];
        int len = snprintf(mine, sizeof mine, "C %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64, idx[i], D.tot.surv,
                           D.tot.cand, D.tot.spsp);
        for (int bb = 0; bb < R.nbins; bb++)
            len += snprintf(mine + len, sizeof mine - (size_t)len, " %" PRIu64 ":%016" PRIx64, D.tot.cnt[bb], D.tot.cks[bb]);
        const bool ok = !strcmp(mine, lines[i]);
        bad += !ok;
        printf("%s %s%s%s\n", ok ? "OK      " : "MISMATCH", lines[i], ok ? "" : "  recomputed: ", ok ? "" : mine);
        fflush(stdout);
    }
    char b1[32];
    printf("recheck: %" PRIu64 " of %zu chunks recomputed, %d mismatches (%s)\n", want, n, bad,
           fmt_dur(now() - t0, b1, sizeof b1));
    return bad ? 1 : 0;
}

/* ------------------------------------------------------------------ */
/* Selftest                                                            */
/* ------------------------------------------------------------------ */

static u64 rng_state = 0x243F6A8885A308D3ULL;
static u64 rnd(void) { return mix64(rng_state += GOLD); }

static int test_arith(void)
{
    int bad = 0;
    /* base-2 strong pseudoprimes (A001262) must pass sprp2 and fail BPSW */
    static const u64 spsp2[] = {2047, 3277, 4033, 4681, 8321, 15841, 29341, 42799, 49141, 52633, 65281, 74665,
                                80581, 85489, 88357, 90751, 104653, 130561, 196093, 220729, 233017, 252601,
                                253241, 256999, 271951, 280601, 314821, 357761, 390937, 458989, 476971, 486737};
    for (size_t i = 0; i < sizeof spsp2 / sizeof *spsp2; i++) {
        const u64 n = spsp2[i];
        if (!sprp2_64(n)) { printf("FAIL sprp2_64(%" PRIu64 ") should pass\n", n); bad++; }
        if (!sprp2_128(n)) { printf("FAIL sprp2_128(%" PRIu64 ") should pass\n", n); bad++; }
        if (slprp64(n)) { printf("FAIL slprp64(%" PRIu64 ") should fail\n", n); bad++; }
    }
    /* strong Lucas pseudoprimes (A217255) must pass slprp and fail sprp2 */
    static const u64 slpsp[] = {5459, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519, 75077, 97439,
                                100127, 113573, 115639, 130139, 155819, 158399, 161027, 162133, 176399, 176471,
                                189419, 192509, 197801, 224369, 230691, 231703, 243629, 253259, 268349, 288919};
    for (size_t i = 0; i < sizeof slpsp / sizeof *slpsp; i++) {
        const u64 n = slpsp[i];
        if (!slprp64(n)) { printf("FAIL slprp64(%" PRIu64 ") should pass\n", n); bad++; }
        if (!slprp128(n)) { printf("FAIL slprp128(%" PRIu64 ") should pass\n", n); bad++; }
        if (sprp2_64(n)) { printf("FAIL sprp2_64(%" PRIu64 ") should fail\n", n); bad++; }
    }
    /* every odd n in [3, 2e5): sprp2 against the reference strong test, BPSW against the truth */
    for (u64 n = 43; n < 200000; n += 2) {
        if (n % 3 == 0 || n % 5 == 0 || n % 7 == 0) continue;
        const bool p = is_prime_ref(n);
        if (sprp2_64(n) != mr_ref(n, 2)) { printf("FAIL sprp2_64 at %" PRIu64 "\n", n); bad++; }
        if (bpsw(n) != p) { printf("FAIL bpsw at %" PRIu64 "\n", n); bad++; }
        if (slprp128(n) != slprp64(n)) { printf("FAIL slprp128 != slprp64 at %" PRIu64 "\n", n); bad++; }
    }
    /* random odd numbers of all sizes up to 2^80, including above 2^63 and 2^64 */
    int primes = 0;
    for (int i = 0; i < 200000; i++) {
        const int bits = 20 + (int)(rnd() % 61);
        u128 n = ((u128)rnd() << 64 | rnd()) >> (128 - bits);
        n |= ((u128)1 << (bits - 1)) | 1;
        if (n % 3 == 0 || n % 5 == 0 || n % 7 == 0 || n < 50) continue;
        const bool p = is_prime_ref(n);
        primes += p;
        const bool s2 = (n >> 64) ? sprp2_128(n) : sprp2_64((u64)n);
        if (s2 != mr_ref(n, 2)) { printf("FAIL sprp2 at %s\n", u128s(n)); bad++; }
        if (bpsw(n) != p) { printf("FAIL bpsw at %s\n", u128s(n)); bad++; }
        if (!(n >> 64)) {
            if (sprp2_128(n) != sprp2_64((u64)n)) { printf("FAIL sprp2_128 != sprp2_64 at %s\n", u128s(n)); bad++; }
            if (slprp128(n) != slprp64((u64)n)) { printf("FAIL slprp128 != slprp64 at %s\n", u128s(n)); bad++; }
        }
    }
    /* products of two primes near 2^32 and 2^40 (hard composites), and primes near 2^64 */
    for (int i = 0; i < 2000; i++) {
        u64 a = (rnd() >> 32) | 1, b = (rnd() >> 24) | 1;
        while (!is_prime_ref(a)) a += 2;
        while (!is_prime_ref(b)) b += 2;
        const u128 n = (u128)a * b;
        if (n > 41 && bpsw(n)) { printf("FAIL bpsw accepts %s = %" PRIu64 " * %" PRIu64 "\n", u128s(n), a, b); bad++; }
    }
    u64 n = UINT64_MAX;
    for (int i = 0; i < 300; i++, n -= 2) {
        if (n % 3 == 0 || n % 5 == 0 || n % 7 == 0) continue;
        if (bpsw(n) != is_prime_ref(n)) { printf("FAIL bpsw near 2^64 at %" PRIu64 "\n", n); bad++; }
    }
    printf("arithmetic: %s (%d random primes seen)\n", bad ? "FAILED" : "ok", primes);
    return bad;
}

/* count by direct testing of every p = 7 (mod 30) */
static u64 brute_count(u64 lo, u64 hi, u64 *cks)
{
    u64 c = 0;
    *cks = 0;
    for (u64 p = lo + (u64)((37 - lo % 30) % 30); p < hi; p += 30)
        if (is_sext_ref(p)) { c++; *cks += hashp(p); }
    return c;
}

static bool run_quiet(u128 lo, u128 hi, u32 w, u32 B, u32 S, u64 c0, u64 c1, int threads, u64 *count, u64 *cks)
{
    opts_t o;
    memset(&o, 0, sizeof o);
    o.w = w; o.B = B; o.S = S ? S : 1u << 18; o.threads = threads; o.interval = 60;
    o.cr0 = c0; o.cr1 = c1;
    configure(&o, lo, hi);
    run_count(threads, NULL, 60, true);
    *count = 0; *cks = 0;
    for (int b = 0; b < R.nbins; b++) { *count += D.tot.cnt[b]; *cks += D.tot.cks[b]; }
    return D.tot.spsp == 0;
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

static int cmd_selftest(int argc, char **argv)
{
    opts_t o;
    u128 pos[1] = {14};
    int npos = argc > 2 && isdigit((unsigned char)argv[2][0]) ? 1 : 0;
    parse_opts(argc, argv, 2, &o, npos, pos);
    const int N = (int)pos[0];
    if (N < 1 || N > 17) die("N must be 1..17");
    g_print = false;
    int bad = test_arith();
    double t0 = now();

    /* small ranges against direct testing, several configurations */
    static const u64 rng[][2] = {{0, 2000000}, {0, 100}, {7, 8}, {8, 97}, {97, 98}, {1000000, 3000000},
                                 {1048576 - 30000, 1048576 + 300000}, {123456789, 223456789}};
    static const u32 cfg[][3] = {{7, 11, 64}, {7, 101, 512}, {13, 1000, 4096}, {17, 1u << 20, 1u << 16},
                                 {23, 1u << 12, 1u << 10}, {29, 1u << 16, 1u << 18}};
    for (size_t i = 0; i < sizeof rng / sizeof *rng; i++) {
        u64 bc, bk;
        const u64 want = brute_count(rng[i][0], rng[i][1], &bc);
        bk = bc;
        for (size_t j = 0; j < sizeof cfg / sizeof *cfg; j++) {
            u64 c, k;
            run_quiet(rng[i][0], rng[i][1], cfg[j][0], cfg[j][1], cfg[j][2], 0, 0, o.threads, &c, &k);
            if (c != want || k != bk) {
                printf("FAIL [%" PRIu64 ", %" PRIu64 ") w=%u B=%u S=%u: %" PRIu64 " %016" PRIx64 ", direct %" PRIu64
                       " %016" PRIx64 "\n", rng[i][0], rng[i][1], cfg[j][0], cfg[j][1], cfg[j][2], c, k, want, bk);
                bad++;
            }
        }
    }
    printf("small ranges vs direct search: %s\n", bad ? "FAILED" : "ok");

    /* A350826(1..N) */
    u128 lo = 1;
    for (int n = 1; n <= N; n++) {
        u64 c, k;
        run_quiet(lo, lo * 10, 0, 1u << 16, 0, 0, 0, o.threads, &c, &k);
        const bool ok = c == A350826[n];
        printf("A350826(%d) = %" PRIu64 " %s (cks %016" PRIx64 ", w %u, %.1fs)\n", n, c, ok ? "ok" : "FAIL", k, W.w,
               now() - t0);
        fflush(stdout);
        bad += !ok;
        lo *= 10;
    }

    /* One w = 23 class r against the classes of larger wheels that refine it: with the
       smallest wheel primes as the most significant digits of the class index, the classes
       r' = r (mod 223092870) of wheel w form a contiguous block of C(w)/C(23) indices. */
    {
        static const struct { u32 w, B, S; } inv[] = {{23, 1u << 16, 1u << 18}, {23, 1u << 12, 1u << 12},
                                                     {29, 1u << 14, 1u << 16}, {31, 1u << 17, 1u << 20},
                                                     {37, 1u << 16, 1u << 18}, {37, 1u << 15, 1u << 10},
                                                     {41, 1u << 12, 1u << 12}, {43, 1u << 11, 1u << 8}};
        static const char *ranges[][2] = {{"1e18", "1e18+4e15"}, {"2^64-2e15", "2^64+2e15"}};
        static const u64 cls23[] = {0, 4321, 85084};
        for (size_t ri = 0; ri < 2; ri++) {
            for (size_t ci = 0; ci < 3; ci++) {
                u64 c0 = 0, k0 = 0, sub = 1;
                for (size_t j = 0; j < sizeof inv / sizeof *inv; j++) {
                    wheel_setup(inv[j].w);
                    sub = W.C / 85085;
                    u64 c, k;
                    run_quiet(parse_num(ranges[ri][0]), parse_num(ranges[ri][1]), inv[j].w, inv[j].B, inv[j].S,
                              cls23[ci] * sub, (cls23[ci] + 1) * sub, o.threads, &c, &k);
                    if (j == 0) { c0 = c; k0 = k; }
                    else if (c != c0 || k != k0) {
                        printf("FAIL invariance [%s, %s) class %" PRIu64 " w=%u B=%u S=%u: %" PRIu64 " %016" PRIx64
                               " vs %" PRIu64 " %016" PRIx64 "\n", ranges[ri][0], ranges[ri][1], cls23[ci], inv[j].w,
                               inv[j].B, inv[j].S, c, k, c0, k0);
                        bad++;
                    }
                }
                printf("invariance [%s, %s) w23-class %" PRIu64 ": count %" PRIu64 " cks %016" PRIx64 "\n",
                       ranges[ri][0], ranges[ri][1], cls23[ci], c0, k0);
            }
        }
        fflush(stdout);
    }

    /* windows against verify_a350826.py */
    for (size_t i = 0; WINDOWS[i].lo; i++) {
        u64 c, k;
        run_quiet(parse_num(WINDOWS[i].lo), parse_num(WINDOWS[i].hi), 0, 1u << 16, 0, 0, 0, o.threads, &c, &k);
        const bool ok = c == WINDOWS[i].count && k == WINDOWS[i].cks;
        printf("window [%s, %s): %" PRIu64 " %016" PRIx64 " %s\n", WINDOWS[i].lo, WINDOWS[i].hi, c, k, ok ? "ok" : "FAIL");
        bad += !ok;
    }
    char b1[32];
    printf("selftest %s (%s)\n", bad ? "FAILED" : "passed", fmt_dur(now() - t0, b1, sizeof b1));
    return bad ? 1 : 0;
}

int main(int argc, char **argv)
{
    stderr_tty = isatty(2);
    setvbuf(stdout, NULL, _IOLBF, 0);
    if (argc < 2) usage();
    if (!strcmp(argv[1], "count")) return cmd_count(argc, argv);
    if (!strcmp(argv[1], "list")) return cmd_list(argc, argv);
    if (!strcmp(argv[1], "check")) return cmd_check(argc, argv);
    if (!strcmp(argv[1], "bench")) return cmd_bench(argc, argv);
    if (!strcmp(argv[1], "recheck")) return cmd_recheck(argc, argv);
    if (!strcmp(argv[1], "selftest")) return cmd_selftest(argc, argv);
    usage();
    return 2;
}
