/*
 * dspattern_gpu.m -- dspattern with a Metal GPU search, for macOS.  This file
 * is the whole program: it includes dspattern.c, which declares the three
 * functions below when DSPATTERN_GPU is defined, and adds one worker thread
 * that takes work items from the same queue as the CPU threads, so
 * checkpoints, --first, the statistics and the status line work unchanged.
 *
 * The GPU takes a unit of items at a time and runs the enumeration and the
 * sieve over it (dspattern.metal, embedded at build time): list_kernel,
 * GPU_KLIST candidates to a SIMD group whose lanes share the high digits and
 * split the list of low ones, or with DSPATTERN_GPU_KERNEL=1 the simpler
 * enum_kernel, GPU_K candidates to a thread.  (list_kernel leaves items with
 * fewer than VS free digits, a few edge blocks, to the CPU.)  The rare starts that pass the sieve come back as
 * (item, rank) pairs; this thread rebuilds each one and screens and proves
 * it on the CPU, exactly as run_item does.  Two units are in flight, so the
 * GPU sieves one while the CPU finishes the other.  Every command buffer's
 * status is checked: a failed unit is redone, and after GPU_RETRIES
 * failures its items are searched on the CPU instead; a survivor list that
 * overflows is enlarged and the unit redone.
 *
 * Tuning: DSPATTERN_GPU_K (candidates per SIMD group, or per thread with
 * DSPATTERN_GPU_KERNEL=1), DSPATTERN_GPU_UNIT
 * (candidates per unit, in millions) and DSPATTERN_GPU_DEPTH (units in
 * flight).  Test hooks: DSPATTERN_GPU_FAIL=N treats every Nth command buffer
 * as failed, and DSPATTERN_GPU_CAP=N starts with N-entry survivor lists.
 *
 * Build (macOS): make.  The Makefile embeds dspattern.metal with xxd.
 */
#define DSPATTERN_GPU 1
#include "dspattern.c"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "dspattern_metal.h"

#define GPU_K        1024              /* enum_kernel: candidates per GPU thread */
#define GPU_KLIST    16384             /* list_kernel: candidates per SIMD group */
#define GPU_UNIT     ((uint64_t)1 << 31) /* candidates per unit, about */
#define GPU_ITEMS    65536             /* items per unit, at most */
#define GPU_CAP      (1u << 20)        /* survivors per unit, to start */
#define GPU_DEPTH    2                 /* units in flight */
#define GPU_MAXDEPTH 4
#define GPU_RETRIES  3
#define WS           (9 * NDIG + 1)
#define GPU_VS       7                 /* list_kernel: the low digits, V */
#define VLMAX        73                /* list directory entries (VS <= 8) */

/* the kernel's structures (dspattern.metal) */
struct gpat {
    uint64_t mask[NSP][MASKW];
    uint32_t wt[NMOD][NDIG + 1];
    uint32_t nine[NDIG + 1];
    uint32_t nineu[4][NDIG + 1];
    uint32_t pad;
};
#define PREFN ((NMOD + 1) & ~1)          /* even, so count needs no padding */
struct gitem {
    uint32_t pat, r, need, no9;
    uint32_t tv, pad;
    uint32_t pref[PREFN];
    uint64_t count, chunk0;
};
struct gvlist { uint32_t start[VLMAX], len[VLMAX]; };
struct gparams { uint32_t nitems, K, cap, pad; uint64_t nchunks; };
struct gconst { uint32_t sp[NSP], mods[NMOD], first[NMOD + 1]; };
struct gsurv { uint32_t item, pad; uint64_t rank; };

_Static_assert(sizeof(struct gpat) % 8 == 0, "gpat layout");
_Static_assert(offsetof(struct gitem, count) == 24 + 4 * PREFN, "gitem layout");
_Static_assert(sizeof(struct gparams) == 24, "gparams layout");
_Static_assert(sizeof(struct gsurv) == 16, "gsurv layout");

/* one unit in flight */
struct gpu_slot {
    id<MTLBuffer> params, items, out, nout;
    uint32_t cap;                       /* entries in out */
    id<MTLCommandBuffer> cb;
    uint64_t *idx;                      /* the items, as indices into g_items */
    struct iparams *ip;
    uint32_t n;                         /* items in the unit */
    uint64_t chunks, cand;              /* GPU threads, and candidates */
    int tries;
};

static id<MTLDevice> g_dev;
static id<MTLCommandQueue> g_queue;
static id<MTLComputePipelineState> g_pipe;
static id<MTLBuffer> g_bPats, g_bW0, g_bW1, g_bConst;
static id<MTLBuffer> g_bListV, g_bVList, g_bTvp;
static int g_kernel = 2;                /* 2: list_kernel, 1: enum_kernel */
static int VS = GPU_VS;                 /* list_kernel's V: its low digits */
static uint32_t VN;                     /* 10^VS */
static uint32_t *g_tvoff;               /* each pattern's V table, as an offset */
static struct gpu_slot g_slot[GPU_MAXDEPTH];
static char g_gpuname[128];
static uint32_t g_K = GPU_K;
static uint64_t g_unit = GPU_UNIT;
static int g_depth = GPU_DEPTH;
static uint64_t g_failevery, g_ncb;

static const char *gpu_name(void)
{
    return g_gpuname;
}

static id<MTLBuffer> gpu_buf(size_t len)
{
    return [g_dev newBufferWithLength:len ? len : 16
                              options:MTLResourceStorageModeShared];
}

/* a tuning override from the environment, within [lo, hi] */
static uint64_t gpu_env(const char *name, uint64_t def, uint64_t lo, uint64_t hi)
{
    const char *e = getenv(name);
    if (!e || !*e)
        return def;
    uint64_t v = strtoull(e, NULL, 10);
    return v < lo ? lo : v > hi ? hi : v;
}

static const char *gpu_setup(void)
{
    static char msg[512];
    g_dev = MTLCreateSystemDefaultDevice();
    if (!g_dev)
        return "no Metal device";
    snprintf(g_gpuname, sizeof g_gpuname, "%s", g_dev.name.UTF8String);

    NSError *err = nil;
    VS = (int)gpu_env("DSPATTERN_GPU_VS", GPU_VS, 4, 7);    /* 8 would shift words by 32 */
    VN = 1;
    for (int i = 0; i < VS; i++)
        VN *= 10;
    NSString *body = [[NSString alloc] initWithBytes:dspattern_metal
                                              length:dspattern_metal_len
                                            encoding:NSUTF8StringEncoding];
    NSString *src = [NSString stringWithFormat:@"#define VS %d\n#define VN %uu\n"
                                                "#define NSP %d\n#define NMOD %d\n"
                                                "#define MASKW %d\n%@",
                                               VS, VN, NSP, NMOD, MASKW, body];
    id<MTLLibrary> lib = [g_dev newLibraryWithSource:src options:nil error:&err];
    if (!lib) {
        snprintf(msg, sizeof msg, "kernel compile failed: %s",
                 err.localizedDescription.UTF8String);
        return msg;
    }
    g_kernel = (int)gpu_env("DSPATTERN_GPU_KERNEL", 2, 1, 2);
    const char *kname = g_kernel == 2 ? "list_kernel" : "enum_kernel";
    id<MTLFunction> f = [lib newFunctionWithName:@(kname)];
    g_pipe = f ? [g_dev newComputePipelineStateWithFunction:f error:&err] : nil;
    if (!g_pipe) {
        snprintf(msg, sizeof msg, "no pipeline for %s", kname);
        return msg;
    }
    g_queue = [g_dev newCommandQueue];

    /* the patterns, in the kernel's layout */
    g_bPats = gpu_buf((size_t)(g_npat ? g_npat : 1) * sizeof(struct gpat));
    struct gpat *gp = g_bPats.contents;
    for (int i = 0; i < g_npat; i++) {
        memcpy(gp[i].mask, g_pat[i].mask, sizeof gp[i].mask);
        memcpy(gp[i].wt, g_pat[i].wt, sizeof gp[i].wt);
        memcpy(gp[i].nine, g_pat[i].nine, sizeof gp[i].nine);
        for (int k = 0; k < 4; k++) {           /* 9s in U's low positions */
            uint64_t acc = 0;
            for (int a = 0; a <= NDIG; a++) {
                gp[i].nineu[k][a] = (uint32_t)acc;
                if (VS + a <= NDIG)
                    acc = (acc + 9ULL * g_pat[i].wt[k][VS + a]) % MODS[k];
            }
        }
    }

    /* list_kernel's tables.  The V (VS-digit strings) with each digit sum,
       in increasing order, all of them or those not ending in 9; and for
       each trailing-9 count t (which fixes the digit weights), each listed
       V's share of p mod 7, 11, 13, 17, 19 and 23, packed in 26 bits (3, 4,
       4, 5, 5 and 5), in list order */
    int smax = 9 * VS;
    g_bListV = gpu_buf(2 * (size_t)VN * sizeof(uint32_t));
    g_bVList = gpu_buf(2 * sizeof(struct gvlist));
    uint32_t *lv = g_bListV.contents;
    struct gvlist *vl = g_bVList.contents;
    uint8_t *dsv = malloc(VN);
    if (!g_bListV || !g_bVList || !dsv)
        return "out of memory for the GPU's tables";
    for (uint32_t v = 0; v < VN; v++)
        dsv[v] = (uint8_t)digitsum64(v);
    for (int no9 = 0; no9 < 2; no9++) {
        uint32_t *cnt = calloc((size_t)smax + 1, sizeof *cnt);
        for (uint32_t v = 0; v < VN; v++)
            if (!(no9 && v % 10 == 9))
                cnt[dsv[v]]++;
        for (int sg = 0, e = 0; sg <= smax; sg++) {
            vl[no9].start[sg] = (uint32_t)e;
            vl[no9].len[sg] = cnt[sg];
            e += (int)cnt[sg];
            cnt[sg] = vl[no9].start[sg];        /* now the next free slot */
        }
        for (uint32_t v = 0; v < VN; v++)       /* in increasing order */
            if (!(no9 && v % 10 == 9))
                lv[no9 * VN + cnt[dsv[v]]++] = v;
        free(cnt);
    }
    free(dsv);
    int tmaxp = -1;
    for (int i = 0; i < g_npat; i++)
        tmaxp = g_pat[i].t > tmaxp ? g_pat[i].t : tmaxp;
    int nclass = tmaxp + 2;                     /* t = -1 .. tmaxp */
    char *made = calloc((size_t)nclass, 1);
    g_tvoff = malloc((size_t)(g_npat ? g_npat : 1) * sizeof *g_tvoff);
    int nmade = 0;
    for (int i = 0; i < g_npat; i++)
        if (!made[g_pat[i].t + 1]) {
            made[g_pat[i].t + 1] = 1;
            nmade++;
        }
    memset(made, 0, (size_t)nclass);
    g_bTvp = gpu_buf((size_t)(nmade ? nmade : 1) * VN * sizeof(uint32_t));
    if (!g_bTvp || !g_tvoff)
        return "out of memory for the GPU's tables";
    uint32_t *tv = g_bTvp.contents, *slot = calloc((size_t)nclass, sizeof *slot);
    static const uint32_t q6[6] = { 7, 11, 13, 17, 19, 23 }, sh6[6] = { 0, 3, 7, 11, 16, 21 };
    for (int i = 0, next = 0; i < g_npat; i++) {
        int c = g_pat[i].t + 1, no9 = g_pat[i].t >= 0;
        if (!made[c]) {
            made[c] = 1;
            slot[c] = (uint32_t)next++;
            const uint32_t *w = g_pat[i].wt[0];
            uint32_t *out = tv + (size_t)slot[c] * VN;
            uint32_t total = vl[no9].start[smax] + vl[no9].len[smax];
            for (uint32_t e = 0; e < total; e++) {
                uint32_t v = lv[no9 * VN + e];
                uint64_t x = 0;
                for (int j = 0; j < VS; j++, v /= 10)
                    x += (uint64_t)(v % 10) * w[j];
                uint32_t pk = 0;
                for (int k = 0; k < 6; k++)
                    pk |= (uint32_t)(x % q6[k]) << sh6[k];
                out[e] = pk;
            }
        }
        g_tvoff[i] = slot[c] * VN;
    }
    free(slot);
    free(made);
    /* the counting tables for unranking, saturated at 2^62 (a thread's rank
       is far smaller, so a saturated count still compares correctly) */
    g_bW0 = gpu_buf((size_t)(NDIG + 1) * WS * 8);
    g_bW1 = gpu_buf((size_t)(NDIG + 1) * WS * 8);
    uint64_t *w0 = g_bW0.contents, *w1 = g_bW1.contents;
    const u128 sat = (u128)1 << 62;
    for (int j = 0; j <= NDIG; j++)
        for (int s = 0; s < WS; s++) {
            u128 a = ways[j][s], b = j ? strings(j, s, 1) : (u128)(s == 0);
            w0[j * WS + s] = (uint64_t)(a < sat ? a : sat);
            w1[j * WS + s] = (uint64_t)(b < sat ? b : sat);
        }
    g_bConst = gpu_buf(sizeof(struct gconst));
    struct gconst *gc = g_bConst.contents;
    for (int i = 0; i < NSP; i++)
        gc->sp[i] = SP[i];
    for (int k = 0; k < NMOD; k++)
        gc->mods[k] = MODS[k];
    for (int k = 0; k <= NMOD; k++)
        gc->first[k] = MOD_FIRST[k];

    g_K = (uint32_t)gpu_env("DSPATTERN_GPU_K", g_kernel == 2 ? GPU_KLIST : GPU_K, 1, 1u << 24);
    g_unit = gpu_env("DSPATTERN_GPU_UNIT", GPU_UNIT >> 20, 1, 1u << 20) << 20;
    g_depth = (int)gpu_env("DSPATTERN_GPU_DEPTH", GPU_DEPTH, 1, GPU_MAXDEPTH);
    g_failevery = gpu_env("DSPATTERN_GPU_FAIL", 0, 0, UINT64_MAX);
    uint32_t cap = (uint32_t)gpu_env("DSPATTERN_GPU_CAP", GPU_CAP, 1, 1u << 30);
    for (int k = 0; k < g_depth; k++) {
        struct gpu_slot *X = &g_slot[k];
        X->params = gpu_buf(sizeof(struct gparams));
        X->items = gpu_buf(GPU_ITEMS * sizeof(struct gitem));
        X->nout = gpu_buf(4);
        X->out = gpu_buf((size_t)cap * sizeof(struct gsurv));
        X->cap = cap;
        X->idx = malloc(GPU_ITEMS * sizeof *X->idx);
        X->ip = malloc(GPU_ITEMS * sizeof *X->ip);
        if (!X->params || !X->items || !X->nout || !X->out || !X->idx || !X->ip)
            return "out of memory for the GPU's buffers";
    }
    return NULL;
}

/* claim the next unit of items; 0 when there are none */
static int gpu_fill(struct gpu_slot *X)
{
    struct gitem *gi = X->items.contents;
    X->n = 0;
    X->chunks = 0;
    X->cand = 0;
    X->tries = 0;
    while (X->n < GPU_ITEMS && X->cand < g_unit) {
        uint64_t i;
        if (!claim_item(&i))
            break;
        struct iparams *ip = &X->ip[X->n];
        item_params(&g_items[i], ip);
        if (!ip->count) {
            item_finished(i, 0, 0, 0);
            continue;
        }
        if (g_kernel == 2 && ip->r < VS) {      /* a small block: the CPU's */
            run_item(i);
            continue;
        }
        struct gitem *g = &gi[X->n];
        g->pat = (uint32_t)g_items[i].pat;
        g->tv = g_tvoff[g_items[i].pat];
        g->r = (uint32_t)ip->r;
        g->need = (uint32_t)ip->need;
        g->no9 = (uint32_t)ip->no9;
        memcpy(g->pref, ip->pref, sizeof g->pref);
        g->count = ip->count;
        g->chunk0 = X->chunks;
        X->idx[X->n++] = i;
        X->chunks += (ip->count + g_K - 1) / g_K;
        X->cand += ip->count;
    }
    return X->n > 0;
}

static void gpu_encode(struct gpu_slot *X)
{
    *(struct gparams *)X->params.contents =
        (struct gparams){ X->n, g_K, X->cap, 0, X->chunks };
    *(uint32_t *)X->nout.contents = 0;
    id<MTLCommandBuffer> cb = [g_queue commandBuffer];
    id<MTLComputeCommandEncoder> e = [cb computeCommandEncoder];
    [e setComputePipelineState:g_pipe];
    [e setBuffer:X->params offset:0 atIndex:0];
    [e setBuffer:g_bPats offset:0 atIndex:1];
    [e setBuffer:X->items offset:0 atIndex:2];
    [e setBuffer:g_bW0 offset:0 atIndex:3];
    [e setBuffer:g_bW1 offset:0 atIndex:4];
    [e setBuffer:g_bConst offset:0 atIndex:5];
    [e setBuffer:X->out offset:0 atIndex:6];
    [e setBuffer:X->nout offset:0 atIndex:7];
    NSUInteger tg = g_pipe.maxTotalThreadsPerThreadgroup;
    if (tg > 256)
        tg = 256;
    if (g_kernel == 2) {                        /* a chunk per SIMD group */
        NSUInteger w = g_pipe.threadExecutionWidth, nsg = tg / w;
        [e setBuffer:g_bListV offset:0 atIndex:8];
        [e setBuffer:g_bVList offset:0 atIndex:9];
        [e setBuffer:g_bTvp offset:0 atIndex:10];
        [e dispatchThreadgroups:MTLSizeMake((X->chunks + nsg - 1) / nsg, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(nsg * w, 1, 1)];
    } else {
        [e dispatchThreads:MTLSizeMake(X->chunks, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
    }
    [e endEncoding];
    [cb commit];
    X->cb = cb;
}

/* wait for a slot: 0 done, -1 failed, 1 the survivor list overflowed */
static int gpu_wait(struct gpu_slot *X)
{
    [X->cb waitUntilCompleted];
    int failed = X->cb.status != MTLCommandBufferStatusCompleted;
    X->cb = nil;
    if (g_failevery && ++g_ncb % g_failevery == 0)
        failed = 1;
    if (failed)
        return -1;
    return *(uint32_t *)X->nout.contents > X->cap;
}

/* prove the unit's survivors on the CPU and finish its items */
static void gpu_results(struct gpu_slot *X)
{
    uint32_t ns = *(uint32_t *)X->nout.contents;
    const struct gsurv *sv = X->out.contents;
    uint64_t *sieved = calloc(X->n, sizeof *sieved), *tests = calloc(X->n, sizeof *tests);
    if (!sieved || !tests) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    for (uint32_t s = 0; s < ns; s++) {
        uint32_t j = sv[s].item;
        sieved[j]++;
        finish_survivor(&g_items[X->idx[j]], &X->ip[j], sv[s].rank, &tests[j]);
    }
    add64(&g_cand, X->cand);
    for (uint32_t j = 0; j < X->n; j++)
        item_finished(X->idx[j], X->ip[j].count, sieved[j], tests[j]);
    free(sieved);
    free(tests);
}

static void *gpu_worker(void *arg)
{
    (void)arg;
    int head = 0, nfly = 0, more = 1;
    for (;;) {
        @autoreleasepool {
            while (more && nfly < g_depth) {
                struct gpu_slot *X = &g_slot[(head + nfly) % g_depth];
                if (!gpu_fill(X)) {
                    more = 0;
                    break;
                }
                gpu_encode(X);
                nfly++;
            }
            if (!nfly)
                break;

            struct gpu_slot *X = &g_slot[head];
            int r = gpu_wait(X);
            if (r == 1) {                       /* enlarge the list and redo */
                uint64_t cap = (uint64_t)*(uint32_t *)X->nout.contents * 5 / 4 + 1024;
                id<MTLBuffer> b = cap < (1ULL << 31) ? gpu_buf(cap * sizeof(struct gsurv)) : nil;
                if (b) {
                    X->out = b;
                    X->cap = (uint32_t)cap;
                    gpu_encode(X);
                    continue;
                }
                r = -1;
                X->tries = GPU_RETRIES;
            }
            if (r == -1 && ++X->tries <= GPU_RETRIES) {
                gpu_encode(X);
                continue;
            }
            if (r == -1) {                      /* give the unit to the CPU */
                pthread_mutex_lock(&g_lock);
                fprintf(stderr, "%snote: the GPU failed a unit %d times; searching its "
                                "%u items on the CPU\n", CLR, X->tries, X->n);
                pthread_mutex_unlock(&g_lock);
                for (uint32_t j = 0; j < X->n; j++)
                    run_item(X->idx[j]);
            } else {
                gpu_results(X);
            }
            head = (head + 1) % g_depth;
            nfly--;
        }
    }
    __atomic_sub_fetch(&g_active, 1, __ATOMIC_RELEASE);
    return NULL;
}
