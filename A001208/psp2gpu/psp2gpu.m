/*
 * psp2gpu.m -- Metal (Apple GPU) version of the exhaustive 2-basis search
 *              (postage stamp problem, h = 2, OEIS A001212).
 *
 * The CPU enumerates admissible prefixes a_1..a_D (same exact prunes as psp2.c);
 * each GPU thread then runs the depth-first search below one prefix with an
 * explicit stack: admissibility, element-wise upper bounds, Challis gaps test,
 * and the first-gap candidate rule for the last element.  Masks are 8 x 32-bit.
 *
 * Build: clang -O2 -fobjc-arc -framework Metal -framework Foundation -o psp2gpu psp2gpu.m
 * Usage: ./psp2gpu -k K -t T [-d prefixdepth] [-b batchsize] [-q]
 *        ./psp2gpu -selftest
 */
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define NW 8
#define MAXK 32

static const char *kernelSrc = R"MSL(
#include <metal_stdlib>
using namespace metal;
#define NW 8
#define MAXK 32
#define MAXDEPTH 8
struct Params { uint K, T, D, pad; uint UB[MAXK+2]; uint CAP[MAXK+2]; uint TMASK[NW]; };
struct Prefix { uint a[MAXK+2]; uint E[NW]; uint S[NW]; };
struct Found  { uint a[MAXK+2]; };

inline void shl_or(thread uint *dst, thread const uint *src, uint s) {
    if (s >= 32*NW) return;
    int ws = s >> 5; uint bs = s & 31;
    for (int i = NW-1; i >= ws; i--) {
        uint v = src[i-ws] << bs;
        if (bs != 0 && i-ws-1 >= 0) v |= src[i-ws-1] >> (32-bs);
        dst[i] |= v;
    }
}
inline uint popc_masked(thread const uint *s, constant uint *tm) {
    uint c = 0; for (int i = 0; i < NW; i++) c += popcount(s[i] & tm[i]); return c;
}
inline uint first_zero(thread const uint *s) {
    for (int i = 0; i < NW; i++) { uint v = ~s[i]; if (v) return i*32 + ctz(v); } return 32*NW;
}

kernel void search(device const Prefix *pre      [[buffer(0)]],
                   constant Params &P            [[buffer(1)]],
                   device atomic_uint *nfound    [[buffer(2)]],
                   device Found *found           [[buffer(3)]],
                   device ulong *nodecount       [[buffer(4)]],
                   constant uint &nprefix        [[buffer(5)]],
                   uint gid [[thread_position_in_grid]])
{
    if (gid >= nprefix) return;
    const uint K = P.K, T = P.T, D = P.D;
    uint a[MAXK+2];
    for (uint i = 0; i < MAXK+2; i++) a[i] = pre[gid].a[i];
    uint S[MAXDEPTH+1][NW], E[NW];                   /* S per level (index j-D); E shared, undone on pop */
    uint cur[MAXDEPTH+1], hi[MAXDEPTH+1];
    for (int i = 0; i < NW; i++) { S[0][i] = pre[gid].S[i]; E[i] = pre[gid].E[i]; }
    ulong nodes = 0;
    int j = D;                                     /* elements a_1..a_j placed */
    { uint fz = first_zero(S[0]); uint h = min(fz, P.UB[D+1]); hi[0] = h; cur[0] = a[D] + 1; }
    while (true) {
        int lv = j - D;
        if ((uint)j + 1 == K) {
            /* leaf: choose a_K covering all gaps */
            uint G[NW]; bool any = false;
            for (int i = 0; i < NW; i++) { G[i] = P.TMASK[i] & ~S[lv][i]; any = any || (G[i] != 0); }
            uint last = a[j];
            if (!any) {
                uint idx = atomic_fetch_add_explicit(nfound, 1u, memory_order_relaxed);
                if (idx < 1024) { for (uint i = 0; i < MAXK+2; i++) found[idx].a[i] = (i == K) ? last+1 : a[i]; }
            } else {
                uint g = 0; for (int i = 0; i < NW; i++) if (G[i]) { g = i*32 + ctz(G[i]); break; }
                uint h = min(g, P.UB[K]);
                for (uint p = 0; p <= (uint)j + 1; p++) {
                    uint z;
                    if (p <= (uint)j) z = g - a[p]; else { if (g & 1) break; z = g / 2; }
                    if (z <= last || z > h) continue;
                    nodes++;
                    uint M[NW]; for (int i = 0; i < NW; i++) M[i] = 0;
                    shl_or(M, E, z);
                    if (2*z <= T) M[(2*z) >> 5] |= 1u << ((2*z) & 31);
                    bool ok = true; for (int i = 0; i < NW; i++) if (G[i] & ~M[i]) { ok = false; break; }
                    if (ok) {
                        uint idx = atomic_fetch_add_explicit(nfound, 1u, memory_order_relaxed);
                        if (idx < 1024) { for (uint i = 0; i < MAXK+2; i++) found[idx].a[i] = (i == K) ? z : a[i]; }
                    }
                }
            }
            j--; if (j < (int)D) break; { uint ax = a[j+1]; E[ax >> 5] &= ~(1u << (ax & 31)); } continue;
        }
        uint x = cur[lv];
        if (x > hi[lv]) { j--; if (j < (int)D) break; { uint ax = a[j+1]; E[ax >> 5] &= ~(1u << (ax & 31)); } continue; }
        cur[lv] = x + 1;
        /* S[lv+1] = S[lv] | E[lv] << x | bit(2x) */
        for (int i = 0; i < NW; i++) S[lv+1][i] = S[lv][i];
        shl_or(S[lv+1], E, x);
        if (2*x <= T) S[lv+1][(2*x) >> 5] |= 1u << ((2*x) & 31);
        nodes++;
        uint gaps = (T + 1) - popc_masked(S[lv+1], P.TMASK);
        if (gaps > P.CAP[j+1]) continue;
        E[x >> 5] |= 1u << (x & 31);
        a[j+1] = x;
        j++;
        { uint fz = first_zero(S[lv+1]); uint h = min(fz, P.UB[j+1]); hi[lv+1] = h; cur[lv+1] = x + 1; }
    }
    nodecount[gid] = nodes;
}
)MSL";

typedef struct { uint32_t K, T, D, pad; uint32_t UB[MAXK+2]; uint32_t CAP[MAXK+2]; uint32_t TMASK[NW]; } Params;
typedef struct { uint32_t a[MAXK+2]; uint32_t E[NW]; uint32_t S[NW]; } Prefix;
typedef struct { uint32_t a[MAXK+2]; } Found;

static const int N2[] = {0, 2, 4, 8, 12, 16, 20, 26, 32, 40, 46, 54, 64, 72, 80,
                         92, 104, 116, 128, 140, 152, 164, 180, 196, 212};
#define N2_KNOWN 24

static int K, T, D, QUIET = 0;
static size_t BATCH = 1 << 20;
static Params P;

/* ---- CPU prefix enumeration (same prunes as psp2.c), streaming in batches ---- */
typedef struct { uint32_t E[NW], S[NW]; } lvl;
static inline void cset(uint32_t *m, int i) { m[i >> 5] |= 1u << (i & 31); }
static inline void cshl_or(uint32_t *dst, const uint32_t *src, int s) {
    int ws = s >> 5, bs = s & 31;
    for (int i = NW - 1; i >= ws; i--) { uint32_t v = src[i - ws] << bs; if (bs && i - ws - 1 >= 0) v |= src[i - ws - 1] >> (32 - bs); dst[i] |= v; }
}
static inline int cpopc(const uint32_t *s) { int c = 0; for (int i = 0; i < NW; i++) c += __builtin_popcount(s[i] & P.TMASK[i]); return c; }
static inline int cfz(const uint32_t *s) { for (int i = 0; i < NW; i++) { uint32_t v = ~s[i]; if (v) return i * 32 + __builtin_ctz(v); } return 32 * NW; }

static Prefix *batchbuf; static size_t nbatch;
static uint64_t total_nodes = 0, total_found = 0, total_prefixes = 0;
static double gpu_time = 0;
static id<MTLDevice> dev; static id<MTLComputePipelineState> pso; static id<MTLCommandQueue> queue;
static id<MTLBuffer> bufPre, bufParams, bufNfound, bufFound, bufNodes, bufN;
static Found *foundlist; static size_t nfoundlist, foundcap;

static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

/* run the GPU on batchbuf[off .. off+n); on a command-buffer error (macOS watchdog
   "Impacting Interactivity" when a dispatch runs too long) split the batch and retry */
static void run_batch(size_t off, size_t n, int depth) {
    if (n == 0) return;
    memcpy(bufPre.contents, batchbuf + off, n * sizeof(Prefix));
    *(uint32_t *)bufNfound.contents = 0;
    *(uint32_t *)bufN.contents = (uint32_t)n;
    id<MTLCommandBuffer> cb = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    [enc setComputePipelineState:pso];
    [enc setBuffer:bufPre offset:0 atIndex:0];
    [enc setBuffer:bufParams offset:0 atIndex:1];
    [enc setBuffer:bufNfound offset:0 atIndex:2];
    [enc setBuffer:bufFound offset:0 atIndex:3];
    [enc setBuffer:bufNodes offset:0 atIndex:4];
    [enc setBuffer:bufN offset:0 atIndex:5];
    NSUInteger tg = pso.maxTotalThreadsPerThreadgroup; if (tg > 256) tg = 256;
    [enc dispatchThreads:MTLSizeMake(n, 1, 1) threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
    [enc endEncoding];
    [cb commit];
    [cb waitUntilCompleted];
    if (cb.error) {
        fprintf(stderr, "GPU error on batch of %zu (depth %d): %s -- splitting and retrying\n", n, depth, cb.error.localizedDescription.UTF8String);
        if (n <= 256 || depth > 12) { fprintf(stderr, "GPU error persists at minimal batch size; giving up\n"); exit(1); }
        run_batch(off, n / 2, depth + 1);
        run_batch(off + n / 2, n - n / 2, depth + 1);
        return;
    }
    uint64_t *nc = (uint64_t *)bufNodes.contents;
    for (size_t i = 0; i < n; i++) total_nodes += nc[i];
    uint32_t nf = *(uint32_t *)bufNfound.contents;
    Found *f = (Found *)bufFound.contents;
    for (uint32_t i = 0; i < nf && i < 1024; i++) {
        if (nfoundlist == foundcap) { foundcap = foundcap ? 2 * foundcap : 64; foundlist = realloc(foundlist, foundcap * sizeof(Found)); }
        foundlist[nfoundlist++] = f[i];
        if (!QUIET) { printf("FOUND k=%d range>=%d:", K, T); for (int q = 1; q <= K; q++) printf(" %u", f[i].a[q]); printf("\n"); fflush(stdout); }
    }
    total_found += nf;
    total_prefixes += n;
}

static void flush_batch(void) {
    if (nbatch == 0) return;
    double t0 = now();
    run_batch(0, nbatch, 0);
    gpu_time += now() - t0;
    nbatch = 0;
    static int nb = 0; static double tstart = 0; if (!tstart) tstart = t0;
    if ((++nb % 16) == 0) { fprintf(stderr, "[progress] batches=%d prefixes=%llu nodes=%.4g found=%llu elapsed=%.0fs rate=%.0fM/s\n", nb, (unsigned long long)total_prefixes, (double)total_nodes, (unsigned long long)total_found, now() - tstart, total_nodes / (now() - tstart) / 1e6); fflush(stderr); }
}

static void gen(int j, uint32_t *a, lvl *L) {
    if (j == D) {
        Prefix *p = &batchbuf[nbatch++];
        memset(p, 0, sizeof *p);
        for (int i = 0; i <= D; i++) p->a[i] = a[i];
        memcpy(p->E, L[j].E, sizeof p->E); memcpy(p->S, L[j].S, sizeof p->S);
        if (nbatch == BATCH) flush_batch();
        return;
    }
    int last = a[j];
    int hi = cfz(L[j].S); if (hi > (int)P.UB[j + 1]) hi = P.UB[j + 1];
    for (int x = last + 1; x <= hi; x++) {
        memcpy(L[j + 1].S, L[j].S, sizeof L[j].S); cshl_or(L[j + 1].S, L[j].E, x); if (2 * x <= T) cset(L[j + 1].S, 2 * x);
        int gaps = (T + 1) - cpopc(L[j + 1].S);
        if (gaps > (int)P.CAP[j + 1]) continue;
        memcpy(L[j + 1].E, L[j].E, sizeof L[j].E); cset(L[j + 1].E, x);
        a[j + 1] = x;
        gen(j + 1, a, L);
    }
}

static void setup_metal(void) {
    dev = MTLCreateSystemDefaultDevice();
    if (!dev) { fprintf(stderr, "no Metal device\n"); exit(1); }
    NSError *err = nil;
    MTLCompileOptions *opts = [MTLCompileOptions new];
    id<MTLLibrary> lib = [dev newLibraryWithSource:[NSString stringWithUTF8String:kernelSrc] options:opts error:&err];
    if (!lib) { fprintf(stderr, "compile error: %s\n", err.localizedDescription.UTF8String); exit(1); }
    id<MTLFunction> fn = [lib newFunctionWithName:@"search"];
    pso = [dev newComputePipelineStateWithFunction:fn error:&err];
    if (!pso) { fprintf(stderr, "pso error: %s\n", err.localizedDescription.UTF8String); exit(1); }
    queue = [dev newCommandQueue];
    bufPre = [dev newBufferWithLength:BATCH * sizeof(Prefix) options:MTLResourceStorageModeShared];
    bufParams = [dev newBufferWithLength:sizeof(Params) options:MTLResourceStorageModeShared];
    bufNfound = [dev newBufferWithLength:4 options:MTLResourceStorageModeShared];
    bufFound = [dev newBufferWithLength:1024 * sizeof(Found) options:MTLResourceStorageModeShared];
    bufNodes = [dev newBufferWithLength:BATCH * sizeof(uint64_t) options:MTLResourceStorageModeShared];
    bufN = [dev newBufferWithLength:4 options:MTLResourceStorageModeShared];
    batchbuf = malloc(BATCH * sizeof(Prefix));
    if (!QUIET) printf("Metal device: %s, maxThreadsPerThreadgroup=%lu\n", dev.name.UTF8String, (unsigned long)pso.maxTotalThreadsPerThreadgroup);
}

static uint64_t run(int k, int t, int d, int verbose) {
    K = k; T = t; D = d;
    if (T > 32 * NW - 1) { fprintf(stderr, "T too large\n"); exit(1); }
    if (D >= K - 1) D = K - 2; if (D < 1) D = 1;
    if (K - D > 8) { fprintf(stderr, "K-D must be <= 8 (MAXDEPTH)\n"); exit(1); }
    memset(&P, 0, sizeof P); P.K = K; P.T = T; P.D = D;
    for (int i = 0; i <= T; i++) cset(P.TMASK, i);
    for (int j = 0; j <= K + 1; j++) {
        int ub = (j >= 1 && j - 1 <= N2_KNOWN) ? N2[j - 1] + 1 : T; if (ub > T) ub = T; P.UB[j] = ub;
        int c = 0; for (int i = j + 1; i <= K; i++) c += i + 1; P.CAP[j] = c;
    }
    memcpy(bufParams.contents, &P, sizeof P);
    total_nodes = total_found = total_prefixes = 0; gpu_time = 0; nfoundlist = 0; nbatch = 0;
    double t0 = now();
    uint32_t a[MAXK + 2] = {0}; lvl L[MAXK + 2]; memset(L, 0, sizeof L);
    a[1] = 1; cset(L[1].E, 0); cset(L[1].E, 1); cset(L[1].S, 0); cset(L[1].S, 1); cset(L[1].S, 2);
    gen(1, a, L);
    flush_batch();
    double t1 = now();
    if (verbose) {
        printf("RESULT k=%d T=%d D=%d prefixes=%llu found=%llu gpu_nodes=%llu time=%.2fs (gpu %.2fs) rate=%.2fM nodes/s\n",
               K, T, D, (unsigned long long)total_prefixes, (unsigned long long)total_found, (unsigned long long)total_nodes,
               t1 - t0, gpu_time, total_nodes / (t1 - t0) / 1e6);
        fflush(stdout);
    }
    return total_found;
}

static void selftest(void) {
    struct { int k, n, count; } cases[] = {
        {8, 32, 2}, {9, 40, 1}, {10, 46, 2}, {11, 54, 4}, {12, 64, 1}, {13, 72, 1}, {14, 80, 3}, {15, 92, 1}, {16, 104, 1}, {17, 116, 1}
    };
    int ok = 1, sq = QUIET; QUIET = 1;
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        int d = cases[i].k - 7; if (d < 2) d = 2;
        double t0 = now();
        uint64_t f = run(cases[i].k, cases[i].n, d, 0);
        uint64_t n1 = total_nodes;
        uint64_t g = run(cases[i].k, cases[i].n + 1, d, 0);
        int good = (f == (uint64_t)cases[i].count && g == 0);
        printf("selftest k=%2d n=%3d D=%d: found %llu (expect %d), at n+1: %llu (expect 0), nodes %llu/%llu, %.2fs %s\n", cases[i].k, cases[i].n, d,
               (unsigned long long)f, cases[i].count, (unsigned long long)g, (unsigned long long)n1, (unsigned long long)total_nodes, now() - t0, good ? "OK" : "FAIL");
        fflush(stdout);
        if (!good) ok = 0;
    }
    QUIET = sq;
    printf(ok ? "SELFTEST PASSED\n" : "SELFTEST FAILED\n");
}

int main(int argc, char **argv) {
    @autoreleasepool {
        int k = 0, t = 0, d = 0, st = 0;
        for (int i = 1; i < argc; i++) {
            if (!strcmp(argv[i], "-k")) k = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-t")) t = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-d")) d = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-b")) BATCH = (size_t)atol(argv[++i]);
            else if (!strcmp(argv[i], "-q")) QUIET = 1;
            else if (!strcmp(argv[i], "-selftest")) st = 1;
            else { fprintf(stderr, "bad arg %s\n", argv[i]); return 1; }
        }
        setup_metal();
        if (st) { selftest(); return 0; }
        if (k < 3 || t < 1) { fprintf(stderr, "usage: psp2gpu -k K -t T [-d prefixdepth] [-b batch] [-q] | -selftest\n"); return 1; }
        if (d == 0) d = k - 7; if (d < 2) d = 2;
        run(k, t, d, 1);
    }
    return 0;
}
