/*
 * a045875_gpu.cu -- CUDA search for terms of OEIS A045875 (n >= 17 only):
 *   a(n) = smallest m such that the decimal representation of 2^m
 *          contains n consecutive identical digits.
 *
 * Same state and checkpoint format as the CPU tool (a045875.c): the decimal
 * expansion of 2^m as base-10^9 uint32 limbs, least significant first.
 * Checkpoints are interchangeable between the two programs.
 *
 * WHY A GPU HELPS -- and the batching trick that makes it big:
 * Doubling in base 10^9 has purely local carries (carry into limb i is
 * (old a[i-1] >= 5*10^8)), and after j doublings limb i depends only on
 * initial limbs i-j..i.  So each GPU thread takes a span of limbs plus a
 * K-limb "halo" below it, runs K=24 doublings entirely in registers
 * (carries as bitmasks), checks EVERY intermediate value for repdigit
 * limbs, and writes the K-step result back.  Memory traffic per doubling
 * drops by 24x versus the CPU tool, which is why this is worth a GPU.
 *
 * Detection: any run of >= 17 identical digits contains a fully-aligned
 * 9-digit repdigit limb (value d*111111111, d=0..9).  The kernel flags
 * those; the host confirms exact run lengths by recomputing a small window
 * of limbs at the exact intermediate step from the pre-batch state (the
 * same triangle dependency makes that a few hundred operations).
 * The aligned-limb filter is only complete for runs >= 17, hence n >= 17;
 * use the CPU tool for smaller n.
 *
 * LIBRARIES: CUDA Toolkit only.  GMP is used for -s initialization when
 * available; on boxes with only the runtime (libgmp.so.10, no gmp.h) the
 * needed five entry points are declared locally against the long-stable
 * ABI, and the conversion is verified against an independently computed
 * 2^s mod 10^9 before the search starts.
 *
 * Build on DGX Spark (GB10 = sm_121, aarch64):
 *     nvcc -O3 -arch=sm_121 a045875_gpu.cu -o a045875gpu -lgmp
 *
 * Usage mirrors the CPU tool:
 *     ./a045875gpu -n 18 -s 356677212 -c a18.ckpt
 *     (-e end_exp, -i ckpt_interval_sec; -t is accepted and ignored)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <signal.h>
#include <errno.h>
#include <cuda_runtime.h>

#if defined(__has_include)
#  if __has_include(<gmp.h>)
#    include <gmp.h>
#    define HAVE_GMP_H 1
#  endif
#endif
#ifndef HAVE_GMP_H
/* Minimal declarations against libgmp's stable ABI (symbols are the real
 * __gmpz_* names; the mpz struct layout has been unchanged for decades). */
extern "C" {
typedef struct { int _mp_alloc; int _mp_size; void *_mp_d; } __mpz_struct;
typedef __mpz_struct mpz_t[1];
void __gmpz_init(mpz_t);
void __gmpz_clear(mpz_t);
void __gmpz_set_ui(mpz_t, unsigned long);
void __gmpz_mul_2exp(mpz_t, const mpz_t, unsigned long);
char *__gmpz_get_str(char *, int, const mpz_t);
}
#define mpz_init      __gmpz_init
#define mpz_clear     __gmpz_clear
#define mpz_set_ui    __gmpz_set_ui
#define mpz_mul_2exp  __gmpz_mul_2exp
#define mpz_get_str   __gmpz_get_str
#endif

typedef uint32_t u32;
typedef uint64_t u64;

#define BASE 1000000000u
#define HALF  500000000u
#define REP   111111111u
#define KBATCH 28                 /* doublings per kernel pass (<= 29) */
#define SPAN   64                 /* limbs owned per GPU thread */
#define BT     128                /* threads per block */
#define TILE   ((u64)BT * SPAN)   /* limbs per block */

/* The limb array lives in global memory in a tile-transposed layout so that
 * warp accesses are fully coalesced: logical limb g = tile*TILE + t*SPAN + i
 * (thread t of the tile, i-th limb of its span) is stored at physical index
 * tile*TILE + i*BT + t.  Checkpoints stay in canonical (linear) order. */
static inline u64 phys(u64 g){
    u64 tb = g / TILE, w = g % TILE;
    return tb * TILE + (w % SPAN) * BT + (w / SPAN);
}
static inline u32 limb(const u32 *st, u64 g){ return st[phys(g)]; }
#define CAND_CAP 65536

#define CUCHK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    fprintf(stderr, "\nCUDA error %s at %s:%d\n", \
            cudaGetErrorString(e_), __FILE__, __LINE__); exit(1); } } while (0)

/* ---------------- globals ---------------- */
static u32 *bufA, *bufB;          /* managed, ping-pong */
static size_t cap_limbs;
struct Cand { u32 x, j; };
static Cand *cands;               /* managed */
static u32 *ncand;                /* managed */

static unsigned n_target = 18;
static u64 start_m = 0, end_m = 0;
static char *ckpt_path = NULL;
static double ckpt_interval = 600.0;

static volatile sig_atomic_t g_sig;
static void on_sig(int){ g_sig = 1; }

static double now_sec(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}
static unsigned ndig_u32(u32 v){ unsigned n = 1; while (v >= 10){ v /= 10; n++; } return n; }

/* ---------------- kernel ---------------- */
template<int KK>
__global__ void kdouble(const u32 *__restrict__ in, u32 *__restrict__ out,
                        u64 Lpad, u64 Ltop, Cand *cands, u32 *ncand)
{
    u64 tile = blockIdx.x;
    int  t   = threadIdx.x;
    u64 o = tile * TILE + (u64)t * SPAN;      /* first owned logical limb */
    if (o >= Lpad) return;
    int nown = (int)((Lpad - o < SPAN) ? (Lpad - o) : SPAN);

    /* strided pointers: element i of a span sits at p[i*BT] */
    const u32 *pin  = in  + tile * TILE + t;
    u32       *pout = out + tile * TILE + t;

    u32 carr = 0;                 /* bit j = carry into current limb at step j+1 */

    /* halo: the KK limbs below o = tail of the previous span (or zeros at g<0) */
    const u32 *ph = NULL;
    if (t > 0)            ph = in + tile * TILE + (t - 1);
    else if (tile > 0)    ph = in + (tile - 1) * TILE + (BT - 1);
    if (ph){
        for (int i = SPAN - KK; i < SPAN; i++){
            u32 v = ph[(u64)i * BT], newc = 0;
            #pragma unroll
            for (int j = 0; j < KK; j++){
                u32 myc = (v >= HALF);
                v = 2u*v + ((carr >> j) & 1u) - myc * BASE;
                newc |= myc << j;
            }
            carr = newc;
        }
    }
    for (int i = 0; i < nown; i++){
        u32 v = pin[(u64)i * BT], newc = 0;
        u64 g = o + i;
        #pragma unroll
        for (int j = 0; j < KK; j++){
            u32 myc = (v >= HALF);
            v = 2u*v + ((carr >> j) & 1u) - myc * BASE;
            newc |= myc << j;
            if (v % REP == 0u && !(v == 0u && g >= Ltop)){
                u32 slot = atomicAdd(ncand, 1u);
                if (slot < CAND_CAP){ cands[slot].x = (u32)g; cands[slot].j = j + 1; }
            }
        }
        carr = newc;
        pout[(u64)i * BT] = v;
    }
}

static void launch(int k, const u32 *in, u32 *out, u64 Lpad, u64 Ltop)
{
    u32 nblocks = (u32)((Lpad + TILE - 1) / TILE);
    if (k == KBATCH) kdouble<KBATCH><<<nblocks, BT>>>(in, out, Lpad, Ltop, cands, ncand);
    else if (k == 1) kdouble<1><<<nblocks, BT>>>(in, out, Lpad, Ltop, cands, ncand);
    else { fprintf(stderr, "bad k\n"); exit(1); }
    CUCHK(cudaGetLastError());
}

/* ---------------- host-side exact confirmation ----------------
 * Values of limbs [lo, hi] after j doublings of `st`, using the triangle
 * dependency: simulate [max(0,lo-j), hi]; errors from the unknown carry at
 * the window bottom stay below index lo.  Returns pointer to value of limb
 * `lo` inside buf, i.e. ret[i] = limb lo+i at step j.  hi-lo+j must be small.
 */
#define SIMMAX 128
static u32 simbuf[SIMMAX];
static const u32 *vals_at_step(const u32 *st, u64 Lpad, long lo, long hi, int j)
{
    long lo0 = lo - j; if (lo0 < 0) lo0 = 0;
    long len = hi - lo0 + 1;
    if (len > SIMMAX){ fprintf(stderr, "sim window too big\n"); exit(1); }
    for (long i = 0; i < len; i++){
        long g = lo0 + i;
        simbuf[i] = (g < (long)Lpad) ? limb(st, (u64)g) : 0;
    }
    for (int s = 1; s <= j; s++){
        for (long i = len - 1; i >= 1; i--){
            u32 v = 2u*simbuf[i] + (simbuf[i-1] >= HALF);
            if (v >= BASE) v -= BASE;
            simbuf[i] = v;
        }
        {   u32 v = 2u*simbuf[0];       /* carry-in only correct if lo0==0 */
            if (v >= BASE) v -= BASE;
            simbuf[0] = v;
        }
    }
    return simbuf + (lo - lo0);
}

/* digit count and top limb index of the number at step j (from state st) */
static void topinfo_at_step(const u32 *st, u64 Lpad, int j,
                            u64 *top_idx, u64 *digits)
{
    long hi = (long)Lpad - 1;
    long lo = hi - 2; if (lo < 0) lo = 0;
    const u32 *v = vals_at_step(st, Lpad, lo, hi, j);
    long t = hi;
    while (t > lo && v[t - lo] == 0) t--;
    *top_idx = (u64)t;
    *digits = 9ull * (u64)t + ndig_u32(v[t - lo]);
}

/* Exact run length around candidate limb x at step j; 0 if below n or bogus. */
struct Hit { int valid; u64 m; unsigned digit; u64 len, pos_lsd, digits; };

static void confirm(const u32 *st, u64 Lpad, u32 x, int j, u64 m_batch,
                    unsigned n, Hit *best)
{
    u64 top_idx, digits;
    topinfo_at_step(st, Lpad, j, &top_idx, &digits);
    if ((u64)x > top_idx) return;

    long lo = (long)x - 10; if (lo < 0) lo = 0;
    long hi = (long)x + 10; if (hi > (long)top_idx) hi = (long)top_idx;
    const u32 *v = vals_at_step(st, Lpad, lo, hi, j);

    u32 vx = v[x - lo];
    if (vx % REP != 0){
        fprintf(stderr, "\nwarning: GPU/host disagree at limb %u step %d\n", x, j);
        return;
    }
    unsigned d = vx / REP;
    if (vx == 0 && (u64)x == top_idx) return;      /* leading-zero padding */

    u64 len = 9, pos = (u64)x * 9;
    for (long y = (long)x - 1; y >= lo; y--){      /* extend down (toward LSD) */
        u32 w = v[y - lo]; unsigned cnt = 0;
        static const u32 P10[9] = {1,10,100,1000,10000,100000,1000000,10000000,100000000};
        for (int k = 8; k >= 0; k--){
            if ((w / P10[k]) % 10 == d) cnt++; else break;
        }
        len += cnt; pos -= cnt;
        if (cnt < 9) break;
    }
    for (long y = (long)x + 1; y <= hi; y++){      /* extend up (toward MSD) */
        u32 w = v[y - lo];
        unsigned nd = ((u64)y == top_idx) ? ndig_u32(w) : 9;
        unsigned cnt = 0;
        while (cnt < nd && w % 10 == d){ w /= 10; cnt++; }
        len += cnt;
        if (cnt < nd) break;
    }
    if (len < n) return;

    u64 m_hit = m_batch + j;
    if (!best->valid || m_hit < best->m || (m_hit == best->m && pos < best->pos_lsd)){
        best->valid = 1; best->m = m_hit; best->digit = d;
        best->len = len; best->pos_lsd = pos; best->digits = digits;
    }
}

/* ---------------- init / checkpoint (format-compatible with CPU tool) ---- */
struct ckpt_hdr { char magic[8]; u64 version, m, L; };

static void write_ckpt(const u32 *st, u64 m, u64 L)
{
    if (!ckpt_path) return;
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", ckpt_path);
    FILE *f = fopen(tmp, "wb");
    if (!f){ fprintf(stderr, "checkpoint: cannot open %s: %s\n", tmp, strerror(errno)); return; }
    ckpt_hdr h; memcpy(h.magic, "A045875C", 8); h.version = 1; h.m = m; h.L = L;
    int ok = fwrite(&h, sizeof h, 1, f) == 1;
    static u32 chunk[1 << 18];
    for (u64 g0 = 0; ok && g0 < L; g0 += 1 << 18){
        u64 nn = L - g0; if (nn > (1 << 18)) nn = 1 << 18;
        for (u64 i = 0; i < nn; i++) chunk[i] = limb(st, g0 + i);
        ok = fwrite(chunk, sizeof(u32), nn, f) == nn;
    }
    if (!ok){
        fprintf(stderr, "checkpoint: write failed\n"); fclose(f); return;
    }
    fclose(f);
    if (rename(tmp, ckpt_path) != 0)
        fprintf(stderr, "checkpoint: rename failed: %s\n", strerror(errno));
}

static int load_ckpt(u32 *st, u64 *m, u64 *L)
{
    FILE *f = fopen(ckpt_path, "rb");
    if (!f) return 0;
    ckpt_hdr h;
    if (fread(&h, sizeof h, 1, f) != 1 ||
        memcmp(h.magic, "A045875C", 8) != 0 || h.version != 1){
        fprintf(stderr, "checkpoint: bad header in %s\n", ckpt_path); exit(1);
    }
    if (h.L > cap_limbs){ fprintf(stderr, "checkpoint larger than buffer\n"); exit(1); }
    static u32 chunk[1 << 18];
    for (u64 g0 = 0; g0 < h.L; g0 += 1 << 18){
        u64 nn = h.L - g0; if (nn > (1 << 18)) nn = 1 << 18;
        if (fread(chunk, sizeof(u32), nn, f) != nn){
            fprintf(stderr, "checkpoint: truncated\n"); exit(1);
        }
        for (u64 i = 0; i < nn; i++) st[phys(g0 + i)] = chunk[i];
    }
    fclose(f);
    *m = h.m; *L = h.L;
    return 1;
}

static u64 pow2_mod(u64 e, u64 mod)   /* 2^e mod `mod`, for the init self-check */
{
    unsigned __int128 r = 1, b = 2;
    while (e){ if (e & 1) r = r * b % mod; b = b * b % mod; e >>= 1; }
    return (u64)r;
}

static u64 init_from_gmp(u32 *st, u64 s)
{
    fprintf(stderr, "initializing 2^%llu via GMP...\n", (unsigned long long)s);
    double t0 = now_sec();
    mpz_t z; mpz_init(z);
    mpz_set_ui(z, 1);
    mpz_mul_2exp(z, z, (unsigned long)s);
    char *str = mpz_get_str(NULL, 10, z);
    mpz_clear(z);
    size_t len = strlen(str);
    u64 L = (len + 8) / 9;
    if (L > cap_limbs){ fprintf(stderr, "start exponent exceeds buffer\n"); exit(1); }
    u64 j = 0;
    for (size_t i = len; i > 0; ){
        size_t st0 = (i >= 9) ? i - 9 : 0;
        u32 v = 0;
        for (size_t k = st0; k < i; k++) v = v*10 + (u32)(str[k] - '0');
        st[phys(j)] = v; j++;
        i = st0;
    }
    free(str);
    /* self-check: lowest limb must equal 2^s mod 10^9, digit count must match */
    if (limb(st, 0) != (u32)pow2_mod(s, BASE)){
        fprintf(stderr, "GMP init self-check FAILED (2^s mod 1e9 mismatch)\n"); exit(1);
    }
    fprintf(stderr, "init done: %zu digits, %llu limbs (%.1fs), self-check ok\n",
            len, (unsigned long long)L, now_sec() - t0);
    return L;
}

/* full generic scan of the initial state (any run length) */
static int scan_initial(const u32 *st, u64 L, unsigned n, u64 m, Hit *hit)
{
    unsigned cur = 10; u64 len = 0, start_pos = 0; int over = 0;
    u64 digits = 9*(L-1) + ndig_u32(limb(st, L-1));
    for (u64 i = 0; i < L; i++){
        u32 v = limb(st, i);
        unsigned nd = (i == L-1) ? ndig_u32(v) : 9;
        for (unsigned k = 0; k < nd; k++){
            unsigned dig = v % 10; v /= 10;
            if (dig == cur) len++;
            else {
                if (over){ hit->valid=1; hit->m=m; hit->digit=cur; hit->len=len;
                           hit->pos_lsd=start_pos; hit->digits=digits; return 1; }
                cur = dig; len = 1; start_pos = i*9 + k;
            }
            if (len >= n) over = 1;
        }
    }
    if (over){ hit->valid=1; hit->m=m; hit->digit=cur; hit->len=len;
               hit->pos_lsd=start_pos; hit->digits=digits; return 1; }
    (void)n; return 0;
}

static void report(const Hit *h, u64 from, double el, u64 steps)
{
    printf("RESULT: smallest m >= %llu with %u consecutive identical digits "
           "in 2^m is m=%llu\n",
           (unsigned long long)from, n_target, (unsigned long long)h->m);
    printf("==================== FINAL REPORT ====================\n");
    printf("A045875 candidate:   a(%llu) = %llu\n",
           (unsigned long long)h->len, (unsigned long long)h->m);
    printf("run found:           %llu consecutive '%u' digits\n",
           (unsigned long long)h->len, h->digit);
    printf("number:              2^%llu  (%llu decimal digits)\n",
           (unsigned long long)h->m, (unsigned long long)h->digits);
    printf("run position:        digits %llu..%llu from the least significant digit\n"
           "                     (%llu digits in from the most significant digit)\n",
           (unsigned long long)h->pos_lsd,
           (unsigned long long)(h->pos_lsd + h->len - 1),
           (unsigned long long)(h->digits - h->pos_lsd - h->len));
    printf("search:              started at m=%llu, %llu doublings, "
           "%.1f hours, avg %.1f steps/s (GPU)\n",
           (unsigned long long)from, (unsigned long long)steps,
           el / 3600.0, el > 0 ? steps / el : 0.0);
    printf("======================================================\n");
    fflush(stdout);
}

static void usage(const char *p)
{
    fprintf(stderr,
"usage: %s [-n runlen>=17] [-s start_exp] [-e end_exp]\n"
"          [-c checkpoint_file] [-i ckpt_interval_sec]\n"
"  A045875(18): %s -n 18 -s 356677212 -c a18.ckpt\n", p, p);
    exit(2);
}

int main(int argc, char **argv)
{
    int opt;
    while ((opt = getopt(argc, argv, "n:s:e:t:c:i:h")) != -1){
        switch (opt){
        case 'n': n_target = (unsigned)strtoul(optarg, NULL, 10); break;
        case 's': start_m = strtoull(optarg, NULL, 10); break;
        case 'e': end_m = strtoull(optarg, NULL, 10); break;
        case 't': fprintf(stderr, "(-t ignored: GPU version)\n"); break;
        case 'c': ckpt_path = optarg; break;
        case 'i': ckpt_interval = atof(optarg); break;
        default: usage(argv[0]);
        }
    }
    if (n_target < 17){
        fprintf(stderr, "n must be >= 17 (aligned-limb filter); use the CPU tool for n <= 16\n");
        return 2;
    }

    cudaDeviceProp prop; CUCHK(cudaGetDeviceProperties(&prop, 0));
    fprintf(stderr, "GPU: %s, %d SMs, cc %d.%d\n",
            prop.name, prop.multiProcessorCount, prop.major, prop.minor);

    u64 cap_m = end_m ? end_m + KBATCH + 16
                      : (start_m > 5000000000ull ? start_m * 2 : 10000000000ull);
    cap_limbs = (size_t)((long double)cap_m * 0.30102999566398119521L / 9.0L) + 64;
    cap_limbs = (cap_limbs + TILE - 1) / TILE * TILE;
    CUCHK(cudaMallocManaged(&bufA, cap_limbs * sizeof(u32)));
    CUCHK(cudaMallocManaged(&bufB, cap_limbs * sizeof(u32)));
    CUCHK(cudaMallocManaged(&cands, CAND_CAP * sizeof(Cand)));
    CUCHK(cudaMallocManaged(&ncand, sizeof(u32)));
    CUCHK(cudaMemset(bufA, 0, cap_limbs * sizeof(u32)));
    CUCHK(cudaMemset(bufB, 0, cap_limbs * sizeof(u32)));

    signal(SIGINT, on_sig);
    signal(SIGTERM, on_sig);

    u64 m, L;
    if (ckpt_path && load_ckpt(bufA, &m, &L)){
        start_m = m;
        fprintf(stderr, "resumed from %s: m=%llu, %llu digits\n", ckpt_path,
                (unsigned long long)m,
                (unsigned long long)(9*(L-1) + ndig_u32(limb(bufA, L-1))));
        if (limb(bufA, 0) != (u32)pow2_mod(m, BASE)){
            fprintf(stderr, "checkpoint self-check FAILED (2^m mod 1e9 mismatch)\n");
            return 1;
        }
    } else if (start_m == 0){
        bufA[phys(0)] = 1; L = 1; m = 0;
    } else {
        L = init_from_gmp(bufA, start_m); m = start_m;
    }

    Hit best; memset(&best, 0, sizeof best);
    if (scan_initial(bufA, L, n_target, m, &best) && best.len >= n_target){
        report(&best, start_m, 0, 0);
        return 0;
    }
    memset(&best, 0, sizeof best);

    fprintf(stderr, "searching: n=%u start=%llu%s digits=%llu K=%d\n",
            n_target, (unsigned long long)start_m, end_m ? " (bounded)" : "",
            (unsigned long long)(9*(L-1) + ndig_u32(limb(bufA, L-1))), KBATCH);

    u32 *in = bufA, *out = bufB;
    double t0 = now_sec(), t_pr = t0, t_ck = t0;
    u64 m_pr = m;
    int interrupted = 0;

    while (!best.valid){
        if (g_sig){ interrupted = 1; break; }
        if (end_m && m >= end_m) break;
        int k = KBATCH;
        if (end_m && m + k > end_m) k = 1;

        *ncand = 0;
        launch(k, in, out, L + 1, L);
        CUCHK(cudaDeviceSynchronize());

        u32 nc = *ncand;
        if (nc > CAND_CAP){ fprintf(stderr, "\ncandidate overflow (%u)\n", nc); return 1; }
        for (u32 c = 0; c < nc; c++)
            confirm(in, L + 1, cands[c].x, (int)cands[c].j, m, n_target, &best);

        if (limb(out, L) != 0) L++;
        m += k;
        u32 *tmp = in; in = out; out = tmp;

        double t = now_sec();
        if (t - t_pr >= 5.0){
            fprintf(stderr, "\rm=%llu digits=%llu  %.1f steps/s   ",
                    (unsigned long long)m,
                    (unsigned long long)(9*(L-1) + ndig_u32(limb(in, L-1))),
                    (double)(m - m_pr) / (t - t_pr));
            fflush(stderr);
            t_pr = t; m_pr = m;
        }
        if (ckpt_path && t - t_ck >= ckpt_interval){
            write_ckpt(in, m, L);
            t_ck = now_sec();
        }
    }
    fputc('\n', stderr);

    double el = now_sec() - t0;
    if (m > start_m && el > 0)
        fprintf(stderr, "summary: GPU, %llu doublings in %.1fs = %.1f steps/s\n",
                (unsigned long long)(m - start_m), el, (m - start_m) / el);

    if (ckpt_path) write_ckpt(in, m, L);

    if (best.valid){
        report(&best, start_m, el, best.m - start_m);
        fprintf(stderr, "note: checkpoint holds 2^%llu (first batch boundary at/after the hit)\n",
                (unsigned long long)m);
        return 0;
    }
    if (interrupted){
        fprintf(stderr, "interrupted; %s\n",
                ckpt_path ? "state checkpointed" : "no checkpoint file given (-c)");
        return 130;
    }
    fprintf(stderr, "no run of %u found up to m=%llu\n",
            n_target, (unsigned long long)end_m);
    return 1;
}
