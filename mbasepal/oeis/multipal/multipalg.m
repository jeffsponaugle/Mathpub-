/*
 * multipalg.m — Metal (Apple GPU) version of multipal: exhaustive search
 * for numbers that are k-digit palindromes in many bases at once.
 *
 * Same method and output as multipal.c (sparse mode): value chunks
 * [S, S + 2^w) are processed one at a time. The host lists every k-digit
 * palindrome of every base in the chunk as arithmetic runs of the
 * innermost digit (pieces of <= 256 values); the GPU scatters each value
 * into one of 65536 buckets by its high bits (kernel gen), then each
 * bucket is bitonic-sorted in threadgroup memory and its runs of equal
 * values counted (kernel bsort). Values with >= 2 representations come
 * back to the host, which re-checks every one of them in exact arithmetic
 * before recording it.
 *
 *   ./multipalg -k 5 -l 1000000000000 -u 124389107494560001 [-r 5] [-N keys] [-q]
 */
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <pthread.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;

#define NB      (1u << 16)
#define BCAP    2048u
#define HMAX    64
#define PIECE   256u
#define PCCAP   (1u << 22)
#define HITCAP  (1u << 20)

static const char *KSRC =
"#include <metal_stdlib>\n"
"using namespace metal;\n"
"struct Piece { ulong v0; ulong step; uint len; uint pad; };\n"
"kernel void gen(device const Piece *pc [[buffer(0)]],\n"
"                device atomic_uint *bcnt [[buffer(1)]],\n"
"                device uint *bkeys [[buffer(2)]],\n"
"                device atomic_uint *flag [[buffer(3)]],\n"
"                constant uint *prm [[buffer(4)]],\n"
"                uint gid [[thread_position_in_grid]])\n"
"{\n"
"    uint np = prm[0], sh = prm[1], cap = prm[2];\n"
"    if (gid >= np) return;\n"
"    Piece p = pc[gid];\n"
"    ulong off = p.v0;\n"
"    ulong mask = (1ul << sh) - 1ul;\n"
"    for (uint j = 0; j < p.len; j++, off += p.step){\n"
"        uint b = (uint)(off >> sh);\n"
"        uint slot = atomic_fetch_add_explicit(&bcnt[b], 1u, memory_order_relaxed);\n"
"        if (slot < cap) bkeys[(ulong)b * cap + slot] = (uint)(off & mask);\n"
"        else atomic_store_explicit(flag, 1u, memory_order_relaxed);\n"
"    }\n"
"}\n"
"kernel void bsort(device const uint *bcnt [[buffer(0)]],\n"
"                  device const uint *bkeys [[buffer(1)]],\n"
"                  device atomic_uint *nhits [[buffer(2)]],\n"
"                  device uint4 *hits [[buffer(3)]],\n"
"                  device atomic_uint *hist [[buffer(4)]],\n"
"                  constant uint *prm [[buffer(5)]],\n"
"                  uint tg [[threadgroup_position_in_grid]],\n"
"                  uint lid [[thread_position_in_threadgroup]],\n"
"                  uint tsz [[threads_per_threadgroup]])\n"
"{\n"
"    threadgroup uint s[2048];\n"
"    threadgroup atomic_uint lh[64];\n"
"    uint cap = prm[0], hitcap = prm[1];\n"
"    uint n = min(bcnt[tg], cap);\n"
"    uint N = 32; while (N < n) N <<= 1;\n"
"    for (uint i = lid; i < 64; i += tsz) atomic_store_explicit(&lh[i], 0u, memory_order_relaxed);\n"
"    for (uint i = lid; i < N; i += tsz) s[i] = (i < n) ? bkeys[(ulong)tg * cap + i] : 0xFFFFFFFFu;\n"
"    threadgroup_barrier(mem_flags::mem_threadgroup);\n"
"    for (uint k = 2; k <= N; k <<= 1){\n"
"        for (uint j = k >> 1; j > 0; j >>= 1){\n"
"            for (uint i = lid; i < N; i += tsz){\n"
"                uint ixj = i ^ j;\n"
"                if (ixj > i){\n"
"                    uint a = s[i], b = s[ixj];\n"
"                    bool up = ((i & k) == 0);\n"
"                    if ((a > b) == up){ s[i] = b; s[ixj] = a; }\n"
"                }\n"
"            }\n"
"            threadgroup_barrier(mem_flags::mem_threadgroup);\n"
"        }\n"
"    }\n"
"    for (uint i = lid; i < n; i += tsz){\n"
"        if (i == 0 || s[i] != s[i - 1]){\n"
"            uint r = 1;\n"
"            while (i + r < n && s[i + r] == s[i]) r++;\n"
"            atomic_fetch_add_explicit(&lh[min(r, 63u)], 1u, memory_order_relaxed);\n"
"            if (r >= 2){\n"
"                uint h = atomic_fetch_add_explicit(nhits, 1u, memory_order_relaxed);\n"
"                if (h < hitcap) hits[h] = uint4(tg, s[i], r, 0u);\n"
"            }\n"
"        }\n"
"    }\n"
"    threadgroup_barrier(mem_flags::mem_threadgroup);\n"
"    for (uint i = lid; i < 64; i += tsz){\n"
"        uint c = atomic_load_explicit(&lh[i], memory_order_relaxed);\n"
"        if (c) atomic_fetch_add_explicit(&hist[i], c, memory_order_relaxed);\n"
"    }\n"
"}\n";

typedef struct { u64 v0, step; u32 len, pad; } Piece;

static int K, LV, REPORT = 0, QUIET = 0;
static u64 LO = 0, HI = 0;
static double TARGET = 4.5e7;          /* keys per chunk (~690 per bucket) */
static u64 g_min[HMAX + 1], g_hist[HMAX + 1], g_keys = 0;

static u128 powsat(u64 b, int e){
    const u128 CAP = (u128)1 << 100;
    u128 r = 1;
    while (e-- > 0){ if (r > CAP / b) return CAP; r *= b; }
    return r;
}
static u64 iroot(u64 x, int k){
    if (k <= 1) return x;
    u64 r = (u64)pow((double)x, 1.0 / k);
    while (r > 0 && powsat(r, k) > x) r--;
    while (powsat(r + 1, k) <= x) r++;
    return r;
}
static double density(u64 v){
    double x = (double)iroot(v, K) + 1.0, y = (double)iroot(v, K - 1);
    if (y < 2) y = 2;
    if (x < 2) x = 2;
    if (y < x) return 0.0;
    int m = K - LV;
    if (m <= 1) return log((y + 0.5) / (x - 0.5));
    return (pow(x - 0.5, 1 - m) - pow(y + 0.5, 1 - m)) / (m - 1);
}

/* exact count of bases in which v is a K-digit palindrome */
static int reps(u64 v){
    int c = 0;
    u64 blo = iroot(v, K) + 1, bhi = iroot(v, K - 1);
    if (blo < 2) blo = 2;
    for (u64 b = blo; b <= bhi; b++){
        u32 d[64]; int n = 0; u64 x = v;
        while (x){ d[n++] = (u32)(x % b); x /= b; }
        if (n != K) continue;
        int ok = 1;
        for (int i = 0; i < n / 2; i++) if (d[i] != d[n - 1 - i]){ ok = 0; break; }
        c += ok;
    }
    return c;
}

/* ---- host enumeration of a chunk's palindromes as arithmetic pieces ----
 * Bases are dealt round-robin to HT threads, each filling its own piece
 * list; the lists are then concatenated into the GPU buffer. */
static Piece *g_pc; static u64 g_npc; static int g_pcover;
static int HT = 4;

typedef struct {
    Piece *pc; u64 n, cap; int over;
    u64 P[8], R[8], b, S, E, blo, bhi;
    int tid; pthread_t th;
} HCtx;
static HCtx *g_hc;

static void add_ap(HCtx *C, u64 v0, u64 step, u64 m){
    while (m > 0){
        if (C->n >= C->cap){ C->over = 1; return; }
        u32 l = m > PIECE ? PIECE : (u32)m;
        Piece *p = &C->pc[C->n++];
        p->v0 = v0 - C->S; p->step = step; p->len = l; p->pad = 0;
        v0 += (u64)l * step; m -= l;
    }
}
static void hgen(HCtx *C, int i, u64 acc){
    if (C->over) return;
    const u64 P = C->P[i], R = C->R[i];
    u64 dmin = (i == 0) ? 1 : 0, dmax = C->b - 1;
    u128 reach = (u128)acc + R;
    if ((u128)C->S > reach){
        u64 t = (u64)(((u128)C->S - reach + P - 1) / P);
        if (t > dmin) dmin = t;
    }
    if (acc > C->E - 1) return;
    u64 t2 = (C->E - 1 - acc) / P;
    if (t2 < dmax) dmax = t2;
    if (dmin > dmax) return;
    if (i == LV - 1){ add_ap(C, acc + dmin * P, P, dmax - dmin + 1); return; }
    for (u64 d = dmin; d <= dmax && !C->over; d++) hgen(C, i + 1, acc + d * P);
}
static void *henum(void *arg){
    HCtx *C = arg;
    for (u64 b = C->blo + (u64)C->tid; b <= C->bhi && !C->over; b += (u64)HT){
        C->b = b;
        for (int i = 0; i < LV; i++){
            int j = K - 1 - i;
            u128 p = (i == j) ? powsat(b, i) : powsat(b, j) + powsat(b, i);
            C->P[i] = (u64)p;
        }
        u128 r = 0;
        for (int i = LV - 1; i >= 0; i--){ C->R[i] = (u64)r; r += (u128)(b - 1) * C->P[i]; }
        hgen(C, 0, 0);
    }
    return NULL;
}
static u64 build_pieces(u64 S, u64 E){
    u64 blo = iroot(S, K) + 1, bhi = iroot(E - 1, K - 1), keys = 0;
    if (blo < 2) blo = 2;
    g_npc = 0; g_pcover = 0;
    for (int t = 0; t < HT; t++){
        HCtx *C = &g_hc[t];
        C->n = 0; C->over = 0; C->S = S; C->E = E; C->blo = blo; C->bhi = bhi; C->tid = t;
        if (HT > 1) pthread_create(&C->th, NULL, henum, C); else henum(C);
    }
    for (int t = 0; t < HT; t++){
        HCtx *C = &g_hc[t];
        if (HT > 1) pthread_join(C->th, NULL);
        if (C->over || g_npc + C->n > PCCAP){ g_pcover = 1; continue; }
        memcpy(g_pc + g_npc, C->pc, C->n * sizeof(Piece));
        g_npc += C->n;
    }
    if (g_pcover) return 0;
    for (u64 i = 0; i < g_npc; i++) keys += g_pc[i].len;
    return keys;
}

static double now_s(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv){
    @autoreleasepool {
    int opt;
    while ((opt = getopt(argc, argv, "k:l:u:r:N:H:q")) != -1){
        switch (opt){
        case 'k': K = atoi(optarg); break;
        case 'l': LO = strtoull(optarg, NULL, 10); break;
        case 'u': HI = strtoull(optarg, NULL, 10); break;
        case 'r': REPORT = atoi(optarg); break;
        case 'N': TARGET = atof(optarg); break;
        case 'H': HT = atoi(optarg); break;
        case 'q': QUIET = 1; break;
        default:
            fprintf(stderr, "usage: %s -k digits -l lo -u hi [-r c] [-N keys/chunk] [-H host_threads] [-q]\n", argv[0]);
            return 2;
        }
    }
    if (K < 3 || K > 10 || LO < 1 || HI <= LO){
        fprintf(stderr, "need 3 <= k <= 10 and 1 <= lo < hi\n");
        return 2;
    }
    LV = (K + 1) / 2;
    for (int c = 0; c <= HMAX; c++){ g_min[c] = UINT64_MAX; g_hist[c] = 0; }

    id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
    if (!dev){ fprintf(stderr, "no Metal device\n"); return 1; }
    NSError *err = nil;
    id<MTLLibrary> lib = [dev newLibraryWithSource:[NSString stringWithUTF8String:KSRC] options:nil error:&err];
    if (!lib){ fprintf(stderr, "kernel compile failed: %s\n", err.localizedDescription.UTF8String); return 1; }
    id<MTLComputePipelineState> psoGen = [dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"gen"] error:&err];
    id<MTLComputePipelineState> psoSort = [dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"bsort"] error:&err];
    if (!psoGen || !psoSort){ fprintf(stderr, "pipeline failed\n"); return 1; }
    id<MTLCommandQueue> queue = [dev newCommandQueue];

    id<MTLBuffer> bPc   = [dev newBufferWithLength:(NSUInteger)PCCAP * sizeof(Piece) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bCnt  = [dev newBufferWithLength:NB * 4 options:MTLResourceStorageModeShared];
    id<MTLBuffer> bKeys = [dev newBufferWithLength:(NSUInteger)NB * BCAP * 4 options:MTLResourceStorageModePrivate];
    id<MTLBuffer> bFlag = [dev newBufferWithLength:4 options:MTLResourceStorageModeShared];
    id<MTLBuffer> bNh   = [dev newBufferWithLength:4 options:MTLResourceStorageModeShared];
    id<MTLBuffer> bHits = [dev newBufferWithLength:(NSUInteger)HITCAP * 16 options:MTLResourceStorageModeShared];
    id<MTLBuffer> bHist = [dev newBufferWithLength:HMAX * 4 options:MTLResourceStorageModeShared];
    id<MTLBuffer> bPg   = [dev newBufferWithLength:16 options:MTLResourceStorageModeShared];
    id<MTLBuffer> bPs   = [dev newBufferWithLength:16 options:MTLResourceStorageModeShared];
    g_pc = (Piece *)bPc.contents;
    if (HT < 1) HT = 1;
    if (HT > 64) HT = 64;
    g_hc = calloc(HT, sizeof *g_hc);
    for (int t = 0; t < HT; t++){
        g_hc[t].cap = PCCAP;
        g_hc[t].pc = malloc((size_t)PCCAP * sizeof(Piece));
        if (!g_hc[t].pc){ fprintf(stderr, "out of memory\n"); return 1; }
    }

    fprintf(stderr, "multipalg: k=%d, v in [%llu, %llu) on %s\n", K,
            (unsigned long long)LO, (unsigned long long)HI, dev.name.UTF8String);

    /* chunk stack: [S, S + 2^w); split on any overflow */
    typedef struct { u64 S; int w; } Ch;
    size_t scap = 1024, sn = 0;
    Ch *stk = malloc(scap * sizeof *stk);
    u64 S = LO, covered = 0, nchunks = 0, nsplit = 0, nhitsTot = 0;
    double t0 = now_s(), tlast = t0;
    for (;;){
        Ch ch;
        if (sn){ ch = stk[--sn]; }
        else {
            if (S >= HI) break;
            double d = density(S);
            int w = (int)floor(log2(TARGET / (d > 1e-30 ? d : 1e-30)));
            if (w > 48) w = 48;
            if (w < 17) w = 17;
            ch.S = S; ch.w = w;
            u64 W = (u64)1 << w;
            S = (HI - S > W) ? S + W : HI;
        }
        u64 cS = ch.S, cE = (HI - cS > ((u64)1 << ch.w)) ? cS + ((u64)1 << ch.w) : HI;
        u64 keys = build_pieces(cS, cE);
        u32 shift = (u32)(ch.w - 16);
        int overflow = g_pcover;
        u32 nh = 0;
        if (!overflow && g_npc){
            memset(bCnt.contents, 0, NB * 4);
            memset(bHist.contents, 0, HMAX * 4);
            *(u32 *)bFlag.contents = 0; *(u32 *)bNh.contents = 0;
            u32 *pg = (u32 *)bPg.contents; pg[0] = (u32)g_npc; pg[1] = shift; pg[2] = BCAP;
            u32 *ps = (u32 *)bPs.contents; ps[0] = BCAP; ps[1] = HITCAP;
            id<MTLCommandBuffer> cb = [queue commandBuffer];
            id<MTLComputeCommandEncoder> e1 = [cb computeCommandEncoder];
            [e1 setComputePipelineState:psoGen];
            [e1 setBuffer:bPc offset:0 atIndex:0]; [e1 setBuffer:bCnt offset:0 atIndex:1];
            [e1 setBuffer:bKeys offset:0 atIndex:2]; [e1 setBuffer:bFlag offset:0 atIndex:3];
            [e1 setBuffer:bPg offset:0 atIndex:4];
            NSUInteger tg1 = psoGen.maxTotalThreadsPerThreadgroup; if (tg1 > 256) tg1 = 256;
            [e1 dispatchThreads:MTLSizeMake((NSUInteger)g_npc, 1, 1) threadsPerThreadgroup:MTLSizeMake(tg1, 1, 1)];
            [e1 endEncoding];
            id<MTLComputeCommandEncoder> e2 = [cb computeCommandEncoder];
            [e2 setComputePipelineState:psoSort];
            [e2 setBuffer:bCnt offset:0 atIndex:0]; [e2 setBuffer:bKeys offset:0 atIndex:1];
            [e2 setBuffer:bNh offset:0 atIndex:2]; [e2 setBuffer:bHits offset:0 atIndex:3];
            [e2 setBuffer:bHist offset:0 atIndex:4]; [e2 setBuffer:bPs offset:0 atIndex:5];
            NSUInteger tg2 = psoSort.maxTotalThreadsPerThreadgroup; if (tg2 > 512) tg2 = 512;
            [e2 dispatchThreadgroups:MTLSizeMake(NB, 1, 1) threadsPerThreadgroup:MTLSizeMake(tg2, 1, 1)];
            [e2 endEncoding];
            [cb commit]; [cb waitUntilCompleted];
            if (cb.error){ fprintf(stderr, "GPU error: %s\n", cb.error.localizedDescription.UTF8String); return 1; }
            nh = *(u32 *)bNh.contents;
            if (*(u32 *)bFlag.contents || nh > HITCAP) overflow = 1;
        }
        if (overflow){
            if (ch.w <= 17){ fprintf(stderr, "cannot split chunk at %llu further\n", (unsigned long long)cS); return 1; }
            if (sn + 2 > scap){ scap *= 2; stk = realloc(stk, scap * sizeof *stk); }
            Ch hi = { cS + ((u64)1 << (ch.w - 1)), ch.w - 1 }, lo = { cS, ch.w - 1 };
            if (hi.S < cE) stk[sn++] = hi;
            stk[sn++] = lo;
            nsplit++;
            continue;
        }
        if (g_npc){
            u32 *hist = (u32 *)bHist.contents;
            for (int c = 1; c < HMAX; c++) g_hist[c] += hist[c];
            uint32_t (*hits)[4] = (uint32_t (*)[4])bHits.contents;
            for (u32 h = 0; h < nh; h++){
                u64 v = cS + ((u64)hits[h][0] << shift) + hits[h][1];
                int r = (int)hits[h][2];
                if (REPORT ? r >= REPORT : r >= 3){
                    int x = reps(v);
                    if (x != r){
                        fprintf(stderr, "VERIFY FAIL: v=%llu gpu=%d host=%d\n", (unsigned long long)v, r, x);
                        return 1;
                    }
                }
                int c = r < HMAX ? r : HMAX;
                if (v < g_min[c]) g_min[c] = v;
                if (REPORT && r >= REPORT){ printf("HIT %llu %d\n", (unsigned long long)v, r); fflush(stdout); }
            }
            nhitsTot += nh;
        }
        g_keys += keys; covered += cE - cS; nchunks++;
        double t = now_s();
        if (!QUIET && t - tlast >= 10.0){
            tlast = t;
            double frac = (double)covered / (double)(HI - LO), el = t - t0;
            int top = 0;
            for (int c = HMAX; c >= 2; c--) if (g_min[c] != UINT64_MAX){ top = c; break; }
            double rate = g_keys / el;
            double eta = frac > 0 ? el * (1 - frac) / frac : 0;
            fprintf(stderr, "[k=%d GPU] %.3f%% done below %llu | %.4g palindromes, %.0fM/s | max reps %d (least %llu) | %.0fs, ETA %.0fs\n",
                    K, 100 * frac, (unsigned long long)cE, (double)g_keys, rate / 1e6, top,
                    top ? (unsigned long long)g_min[top] : 0ULL, el, eta);
        }
    }
    double el = now_s() - t0;
    printf("# multipalg k=%d: all v in [%llu, %llu) searched, %llu palindromes, %llu chunks (%llu splits), %.1fs\n",
           K, (unsigned long long)LO, (unsigned long long)HI, (unsigned long long)g_keys,
           (unsigned long long)nchunks, (unsigned long long)nsplit, el);
    printf("# c  least_v_with_exactly_c  number_of_such_v\n");
    int top = 0;
    for (int c = HMAX - 1; c >= 1; c--) if (g_hist[c]){ top = c; break; }
    for (int c = 2; c <= top; c++)
        if (g_hist[c]) printf("EXACT %d %llu %llu\n", c, (unsigned long long)g_min[c], (unsigned long long)g_hist[c]);
    printf("# n  least_v_with_at_least_n in [lo, hi)\n");
    u64 run = UINT64_MAX;
    u64 atl[HMAX + 1];
    for (int c = top; c >= 2; c--){ if (g_min[c] < run) run = g_min[c]; atl[c] = run; }
    for (int c = 2; c <= top; c++) printf("ATLEAST %d %llu\n", c, (unsigned long long)atl[c]);
    }
    return 0;
}
