/*
 * gapsieve_gpu.m -- gapsieve with a Metal GPU search, for macOS.  This file
 * is the whole program: it includes gapsieve.c, which declares the three
 * functions below when GAPSIEVE_GPU is defined, and adds one worker thread
 * that feeds segments to the GPU from the same queue the CPU threads use,
 * so checkpoints, --first, statistics and the status line are unchanged.
 *
 * The GPU takes a unit of GPU_UNIT consecutive segments at a time: for each
 * shape it sieves every class of every segment and screens every candidate
 * with the base-2 strong test (gapsieve.metal, embedded at build time); the
 * CPU then proves the rare survivors with chain_verify, exactly as it does
 * its own, and completes the unit's segments one by one.  A segment alone
 * left the GPU waiting on dispatches and barriers for about half its time,
 * and on an M2 Ultra its speed swung with the load on the CPU; units give
 * each dispatch several times the work.  GPU_SPLIT screening rounds compact
 * the list member by member, then one kernel tests the members left.  Two
 * units stay in flight, so one's sieve overlaps the other's screening.
 * Every command buffer's status is checked.  A failed unit is redone, and
 * after GPU_RETRIES failures its segments are searched on the CPU instead;
 * a candidate list that overflows is enlarged and the unit redone.
 *
 * Tuning: GAPSIEVE_GPU_UNIT, GAPSIEVE_GPU_SPLIT, GAPSIEVE_GPU_DEPTH and
 * GAPSIEVE_GPU_DENSE override the four defaults below (GAPSIEVE_GPU_ORDER
 * too, see gpu_encode), and GAPSIEVE_GPU_PROF=1 reports where the GPU's time
 * went.  Test hooks: GAPSIEVE_GPU_FAIL=N treats every Nth
 * command buffer as failed, and GAPSIEVE_GPU_CAP=N starts with N-entry
 * candidate lists.
 *
 * Build (macOS): make.  The Makefile embeds gapsieve.metal with xxd.
 */
#define GAPSIEVE_GPU 1
#include "gapsieve.c"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "gapsieve_metal.h"

#define GPU_UNIT 16             /* segments per unit */
#define GPU_SPLIT 3             /* screening rounds before tail_kernel */
#define GPU_DEPTH 2             /* units in flight */
#define GPU_DENSE 3             /* pattern groups sieved before compacting */
#define GPU_MAXUNIT 64
#define GPU_LIST (1u << 23)     /* about the most entries a unit's lists need */
#define GPU_MAXDEPTH 8
#define GPU_RETRIES 3

struct gpu_params {             /* segparams in gapsieve.metal */
    uint64_t base_lo, base_hi, lim, segsize;
    uint32_t W, bmodW, A, ngroups, nps, k, cap, nseg;
    uint32_t dense, pad;
};

/* one unit in flight; each shape has its own lists and counters */
struct gpu_slot {
    id<MTLBuffer> bq;           /* each segment's base mod each sieve prime */
    id<MTLBuffer> own;          /* written in place of g_bOrder */
    id<MTLBuffer> sp[2], roff[2], l0[2], c0[2], l1[2], c1[2], args[2], stat[2];
    uint32_t cap;               /* entries in each list */
    id<MTLCommandBuffer> cb;    /* nil: no GPU work for this unit */
    id<MTLCommandBuffer> cbs;   /* GAPSIEVE_GPU_PROF: the sieves, apart */
    uint64_t s, ns;             /* segments s .. s+ns-1 */
    u128 base;                  /* segment s's base */
    int need[2], fin[2], tries; /* fin: the survivors are in list 1 */
};

static id<MTLDevice> g_dev;
static id<MTLCommandQueue> g_queue;
static id<MTLComputePipelineState> g_pRoffs, g_pSieve, g_pPrep, g_pPass, g_pTail;
static id<MTLBuffer> g_bGroups, g_bQ, g_bQinv, g_bCrt, g_bPat[2], g_bRes[2];
static id<MTLBuffer> g_bOrder;  /* written by every unit: see gpu_encode */
static int g_gorder = 1;
static struct gpu_slot g_slot[GPU_MAXDEPTH];
static int g_gunit = GPU_UNIT, g_gsplit = GPU_SPLIT, g_gdepth = GPU_DEPTH;
static uint32_t g_gdense = GPU_DENSE;
static char g_gpuname[128];
static uint64_t g_failevery, g_ncb;

/* GAPSIEVE_GPU_PROF: GPU seconds in the sieves and the screening, and the
   span from the first command buffer's start to the last one's end */
static int g_prof;
static double g_tsieve, g_tpass, g_tspan0 = -1, g_tspan1;
static uint64_t g_nprof, g_ncand;

static const char *gpu_name(void)
{
    return g_gpuname;
}

static id<MTLBuffer> gpu_buf(size_t len)
{
    return [g_dev newBufferWithLength:len ? len : 16
                              options:MTLResourceStorageModeShared];
}

static int gpu_lists(struct gpu_slot *X, uint32_t cap)
{
    for (int d = 0; d < 2; d++) {
        if (!g_sh[d].on)
            continue;
        X->l0[d] = gpu_buf((size_t)cap * 8);
        X->l1[d] = gpu_buf((size_t)cap * 8);
        if (!X->l0[d] || !X->l1[d])
            return -1;
    }
    X->cap = cap;
    return 0;
}

/* the sieve's expected candidates per segment for shape d: the classes,
   times the fraction of each pattern's period left standing */
static double gpu_expect(int d)
{
    double f = (double)g_sh[d].A * KBLOCKS;
    for (uint32_t g = 0; g < g_ngroups; g++) {
        const uint64_t *p = g_groups[g].pat[d];
        uint64_t P = g_groups[g].P, n = 0;
        for (uint64_t w = 0; w < P / 64; w++)
            n += (uint64_t)__builtin_popcountll(p[w]);
        if (P % 64)
            n += (uint64_t)__builtin_popcountll(p[P / 64] & ((1ULL << (P % 64)) - 1));
        f *= (double)n / (double)P;
    }
    return f;
}

/* a tuning override from the environment, within [lo, hi] */
static int gpu_env(const char *name, int def, int lo, int hi)
{
    const char *e = getenv(name);
    if (!e || !*e)
        return def;
    int v = atoi(e);
    return v < lo ? lo : v > hi ? hi : v;
}

static const char *gpu_setup(void)
{
    static char msg[512];
    if (g_nps != g_nq)
        return "the GPU sieve needs --depth equal to --presieve";
    g_dev = MTLCreateSystemDefaultDevice();
    if (!g_dev)
        return "no Metal device";
    snprintf(g_gpuname, sizeof g_gpuname, "%s", g_dev.name.UTF8String);

    NSError *err = nil;
    NSString *src = [[NSString alloc] initWithBytes:gapsieve_metal
                                             length:gapsieve_metal_len
                                           encoding:NSUTF8StringEncoding];
    id<MTLLibrary> lib = [g_dev newLibraryWithSource:src options:nil error:&err];
    if (!lib) {
        snprintf(msg, sizeof msg, "kernel compile failed: %s",
                 err.localizedDescription.UTF8String);
        return msg;
    }
    id<MTLComputePipelineState> __strong *pp[] = {
        &g_pRoffs, &g_pSieve, &g_pPrep, &g_pPass, &g_pTail };
    const char *names[] = {
        "roffs_kernel", "sieve_kernel", "prep_kernel", "pass_kernel", "tail_kernel" };
    for (int i = 0; i < 5; i++) {
        id<MTLFunction> f = [lib newFunctionWithName:@(names[i])];
        *pp[i] = f ? [g_dev newComputePipelineStateWithFunction:f error:&err] : nil;
        if (!*pp[i]) {
            snprintf(msg, sizeof msg, "no pipeline for %s", names[i]);
            return msg;
        }
    }
    g_queue = [g_dev newCommandQueue];

    /* the pattern groups; each shape's patterns back to back */
    uint32_t *gi = calloc((size_t)g_ngroups * 4, sizeof *gi);
    size_t words = 0;
    for (uint32_t g = 0; g < g_ngroups; g++) {
        gi[4 * g] = g_groups[g].P;
        gi[4 * g + 1] = g_groups[g].first;
        gi[4 * g + 2] = g_groups[g].n;
        gi[4 * g + 3] = (uint32_t)words;
        words += (size_t)(g_groups[g].P + KBLOCKS + 128 + 63) / 64;
    }
    g_bGroups = gpu_buf((size_t)g_ngroups * 16);
    memcpy(g_bGroups.contents, gi, (size_t)g_ngroups * 16);
    free(gi);
    g_bQ = gpu_buf(g_nps * 4);
    g_bQinv = gpu_buf(g_nps * 4);
    g_bCrt = gpu_buf(g_nps * 4);
    memcpy(g_bQ.contents, g_q, g_nps * 4);
    memcpy(g_bQinv.contents, g_qinv, g_nps * 4);
    for (uint32_t i = 0; i < g_nps; i++)
        ((uint32_t *)g_bCrt.contents)[i] = (uint32_t)g_crt[i];
    for (int d = 0; d < 2; d++) {
        if (!g_sh[d].on)
            continue;
        g_bPat[d] = gpu_buf(words * 8);
        for (uint32_t g = 0, w = 0; g < g_ngroups; g++) {
            uint32_t n = (uint32_t)((g_groups[g].P + KBLOCKS + 128 + 63) / 64);
            memcpy((uint64_t *)g_bPat[d].contents + w, g_groups[g].pat[d], (size_t)n * 8);
            w += n;
        }
        g_bRes[d] = gpu_buf(g_sh[d].A * 4);
        memcpy(g_bRes[d].contents, g_sh[d].res, g_sh[d].A * 4);
    }

    g_gunit = gpu_env("GAPSIEVE_GPU_UNIT", GPU_UNIT, 1, GPU_MAXUNIT);
    g_gsplit = gpu_env("GAPSIEVE_GPU_SPLIT", GPU_SPLIT, 0, 1000);
    g_gdepth = gpu_env("GAPSIEVE_GPU_DEPTH", GPU_DEPTH, 1, GPU_MAXDEPTH);
    g_gdense = (uint32_t)gpu_env("GAPSIEVE_GPU_DENSE", GPU_DENSE, 1, 1000);
    if (g_gdense > g_ngroups)
        g_gdense = g_ngroups;
    g_prof = gpu_env("GAPSIEVE_GPU_PROF", 0, 0, 1);
    g_gorder = gpu_env("GAPSIEVE_GPU_ORDER", 1, 0, 2);
    g_bOrder = gpu_buf(4);
    /* short chains leave tens of millions of candidates a segment, so a
       unit's lists are sized from the patterns, with 25% to spare, and the
       unit shrinks to keep them near GPU_LIST entries */
    double est = 0;
    for (int d = 0; d < 2; d++)
        if (g_sh[d].on && gpu_expect(d) > est)
            est = gpu_expect(d);
    double fit = GPU_LIST / (1.25 * est + 1);
    if (g_gunit > fit)
        g_gunit = fit < 1 ? 1 : (int)fit;
    const char *e = getenv("GAPSIEVE_GPU_CAP");
    uint32_t cap = e ? (uint32_t)strtoul(e, NULL, 10)
                     : (uint32_t)(1.25 * est * g_gunit) + 4096;
    if (g_prof)
        note("gapsieve: GPU expects %.0f candidates a segment; %d-segment units, "
             "lists of %u\n", est, g_gunit, cap);
    e = getenv("GAPSIEVE_GPU_FAIL");
    g_failevery = e ? strtoull(e, NULL, 10) : 0;
    for (int k = 0; k < g_gdepth; k++) {
        struct gpu_slot *X = &g_slot[k];
        X->bq = gpu_buf((size_t)g_gunit * g_nps * 4);
        X->own = gpu_buf(4);
        for (int d = 0; d < 2; d++) {
            if (!g_sh[d].on)
                continue;
            X->sp[d] = gpu_buf(sizeof(struct gpu_params));
            X->roff[d] = gpu_buf((size_t)g_gunit * g_sh[d].A * g_ngroups * 4);
            X->c0[d] = gpu_buf(4);
            X->c1[d] = gpu_buf(4);
            X->args[d] = gpu_buf(12);
            X->stat[d] = gpu_buf(4);
        }
        if (gpu_lists(X, cap ? cap : 1))
            return "out of GPU memory";
    }
    return NULL;
}

/*
 * Encode and commit one unit: for each shape it needs, the offsets, the
 * sieve and the screening, with no host round trip in between.  Units that
 * ran side by side slowed each other badly (an M2 Ultra fell from 735 to
 * 263 T/s), but a unit queued behind the running one keeps the GPU busy
 * while the host collects results.  So every unit writes g_bOrder, and
 * Metal's hazard tracking starts each one only once the one before it is
 * done.  (GAPSIEVE_GPU_ORDER=0 lets them overlap; 2 orders just the sieves.)
 */
static void gpu_encode(struct gpu_slot *X)
{
    u128 room = g_hi - X->base;                 /* base <= g_hi here */
    uint64_t span = X->ns * g_segsize - 1;
    uint64_t lim = room > span ? span : (uint64_t)room;
    for (uint64_t g = 0; g < X->ns; g++)
        for (uint32_t i = 0; i < g_nps; i++)
            ((uint32_t *)X->bq.contents)[g * g_nps + i] =
                mod_small(X->base + (u128)g * g_segsize, g_q[i]);

    id<MTLCommandBuffer> cb = [g_queue commandBuffer];
    id<MTLComputeCommandEncoder> e = [cb computeCommandEncoder];
    for (int d = 0; d < 2; d++) {
        if (!X->need[d])
            continue;
        const struct shape *S = &g_sh[d];
        *(struct gpu_params *)X->sp[d].contents = (struct gpu_params){
            (uint64_t)X->base, (uint64_t)(X->base >> 64), lim, g_segsize,
            (uint32_t)g_W, mod_small(X->base, (uint32_t)g_W), S->A, g_ngroups,
            g_nps, (uint32_t)g_k, X->cap, (uint32_t)X->ns, g_gdense, 0 };
        *(uint32_t *)X->c0[d].contents = 0;

        [e setComputePipelineState:g_pRoffs];
        [e setBuffer:X->sp[d] offset:0 atIndex:0];
        [e setBuffer:g_bRes[d] offset:0 atIndex:1];
        [e setBuffer:X->bq offset:0 atIndex:2];
        [e setBuffer:g_bQ offset:0 atIndex:3];
        [e setBuffer:g_bQinv offset:0 atIndex:4];
        [e setBuffer:g_bCrt offset:0 atIndex:5];
        [e setBuffer:g_bGroups offset:0 atIndex:6];
        [e setBuffer:X->roff[d] offset:0 atIndex:7];
        [e setBuffer:g_gorder ? g_bOrder : X->own offset:0 atIndex:8];
        [e dispatchThreads:MTLSizeMake(g_ngroups, S->A, X->ns)
            threadsPerThreadgroup:MTLSizeMake(g_ngroups, 256 / g_ngroups ? 256 / g_ngroups : 1, 1)];
        [e memoryBarrierWithScope:MTLBarrierScopeBuffers];

        [e setComputePipelineState:g_pSieve];
        [e setBuffer:X->sp[d] offset:0 atIndex:0];
        [e setBuffer:g_bRes[d] offset:0 atIndex:1];
        [e setBuffer:g_bPat[d] offset:0 atIndex:2];
        [e setBuffer:g_bGroups offset:0 atIndex:3];
        [e setBuffer:X->roff[d] offset:0 atIndex:4];
        [e setBuffer:X->l0[d] offset:0 atIndex:5];
        [e setBuffer:X->c0[d] offset:0 atIndex:6];
        [e dispatchThreads:MTLSizeMake(KBLOCKS / 32 / 4, S->A, X->ns)
            threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
        [e memoryBarrierWithScope:MTLBarrierScopeBuffers];
    }
    X->cbs = nil;
    if (g_prof || g_gorder == 2) {
        [e endEncoding];
        [cb commit];
        X->cbs = cb;
        cb = [g_queue commandBuffer];
        e = [cb computeCommandEncoder];
    }
    for (int d = 0; d < 2; d++) {
        if (!X->need[d])
            continue;
        const struct shape *S = &g_sh[d];
        id<MTLBuffer> li = X->l0[d], ci = X->c0[d], lo = X->l1[d], co = X->c1[d];
        uint32_t cap = X->cap;
        int rounds = 0;
        for (int m = 0; m < g_k; rounds++) {
            uint32_t first = m == 0;
            [e setComputePipelineState:g_pPrep];
            [e setBuffer:ci offset:0 atIndex:0];
            [e setBuffer:co offset:0 atIndex:1];
            [e setBuffer:X->args[d] offset:0 atIndex:2];
            [e setBytes:&cap length:4 atIndex:3];
            [e setBuffer:X->stat[d] offset:0 atIndex:4];
            [e setBytes:&first length:4 atIndex:5];
            [e setBuffer:g_gorder == 1 ? g_bOrder : X->own offset:0 atIndex:6];
            [e dispatchThreads:MTLSizeMake(1, 1, 1) threadsPerThreadgroup:MTLSizeMake(1, 1, 1)];
            [e memoryBarrierWithScope:MTLBarrierScopeBuffers];
            [e setBuffer:X->sp[d] offset:0 atIndex:0];
            [e setBuffer:li offset:0 atIndex:1];
            [e setBuffer:ci offset:0 atIndex:2];
            [e setBuffer:lo offset:0 atIndex:3];
            [e setBuffer:co offset:0 atIndex:4];
            if (m < g_gsplit) {                 /* one member */
                uint32_t moff = S->offs[m++];
                [e setComputePipelineState:g_pPass];
                [e setBytes:&moff length:4 atIndex:5];
            } else {                            /* all the rest */
                uint32_t nm = (uint32_t)(g_k - m);
                [e setComputePipelineState:g_pTail];
                [e setBytes:&S->offs[m] length:nm * 4 atIndex:5];
                [e setBytes:&nm length:4 atIndex:6];
                m = g_k;
            }
            [e dispatchThreadgroupsWithIndirectBuffer:X->args[d] indirectBufferOffset:0
                                threadsPerThreadgroup:MTLSizeMake(256, 1, 1)];
            [e memoryBarrierWithScope:MTLBarrierScopeBuffers];
            id<MTLBuffer> t = li; li = lo; lo = t;
            t = ci; ci = co; co = t;
        }
        X->fin[d] = rounds % 2;
    }
    [e endEncoding];
    [cb commit];
    X->cb = cb;
}

/* wait for a slot: 0 done, -1 failed, 1 a candidate list overflowed */
static int gpu_wait(struct gpu_slot *X)
{
    [X->cb waitUntilCompleted];
    int failed = X->cb.status != MTLCommandBufferStatusCompleted;
    if (X->cbs) {
        [X->cbs waitUntilCompleted];
        failed |= X->cbs.status != MTLCommandBufferStatusCompleted;
        if (!failed) {
            double s0 = X->cbs.GPUStartTime, s1 = X->cbs.GPUEndTime;
            double p0 = X->cb.GPUStartTime, p1 = X->cb.GPUEndTime;
            g_tsieve += s1 - s0;
            g_tpass += p1 - p0;
            if (g_tspan0 < 0 || s0 < g_tspan0)
                g_tspan0 = s0;
            if (p1 > g_tspan1)
                g_tspan1 = p1;
            g_nprof += X->ns;
            for (int d = 0; d < 2; d++)
                if (X->need[d])
                    g_ncand += *(uint32_t *)X->stat[d].contents;
        }
        X->cbs = nil;
    }
    if (g_failevery && ++g_ncb % g_failevery == 0)
        failed = 1;
    X->cb = nil;
    if (failed)
        return -1;
    for (int d = 0; d < 2; d++)
        if (X->need[d] && *(uint32_t *)X->stat[d].contents > X->cap)
            return 1;
    return 0;
}

/* prove the GPU's survivors in the unit's segment g on the CPU, as
   cpu_segment does its own */
static void gpu_results(struct gpu_slot *X, uint64_t g, struct hit **hits,
                        size_t *nh, size_t *ch)
{
    for (int d = 0; d < 2; d++) {
        if (!X->need[d])
            continue;
        id<MTLBuffer> l = X->fin[d] ? X->l1[d] : X->l0[d];
        id<MTLBuffer> c = X->fin[d] ? X->c1[d] : X->c0[d];
        uint32_t n = *(uint32_t *)c.contents;
        const uint64_t *off = l.contents;
        for (uint32_t i = 0; i < n; i++) {
            if (off[i] / g_segsize != g)
                continue;
            u128 p = X->base + off[i], head = p;
            int plus, len = d == INC ? chain_verify_inc(p, &plus)
                                     : chain_verify_dec(p, &head, &plus);
            if (len)
                add_hit(hits, nh, ch, (struct hit){ head, p, len, plus, d });
        }
    }
}

static void *gpu_worker(void *arg)
{
    (void)arg;
    struct scratch x;                   /* for segments the GPU gives up on */
    int have_x = !scratch_init(&x);
    struct hit *hits = NULL;
    size_t nh = 0, ch = 0;
    int head = 0, nfly = 0, more = 1;

    for (;;) {
        @autoreleasepool {
            while (more && nfly < g_gdepth) {
                struct gpu_slot *X = &g_slot[(head + nfly) % g_gdepth];
                if (!claim_segments(&X->s, &X->ns, (uint64_t)g_gunit, X->need)) {
                    more = 0;
                    break;
                }
                X->base = g_sieve_lo + (u128)X->s * g_segsize;
                X->tries = 0;
                X->cb = nil;
                if (X->base <= g_hi)
                    gpu_encode(X);
                nfly++;
            }
            if (!nfly)
                break;

            struct gpu_slot *X = &g_slot[head];
            int gpu = X->cb != nil;
            int r = gpu ? gpu_wait(X) : 0;
            if (r == 1) {
                uint32_t need = X->cap;
                for (int d = 0; d < 2; d++)
                    if (X->need[d] && *(uint32_t *)X->stat[d].contents > need)
                        need = *(uint32_t *)X->stat[d].contents;
                uint64_t cap = (uint64_t)need + need / 4 + 1024;
                if (cap > UINT32_MAX / 2 || gpu_lists(X, (uint32_t)cap)) {
                    note("gapsieve: GPU candidate lists cannot grow to %" PRIu64
                         "; searching segments %" PRIu64 "-%" PRIu64 " on the CPU\n",
                         cap, X->s, X->s + X->ns - 1);
                    r = -2;
                } else {
                    log_note("GPU candidate lists enlarged to %" PRIu64 "\n", cap);
                    gpu_encode(X);
                    continue;
                }
            }
            if (r == -1 && ++X->tries <= GPU_RETRIES) {
                log_note("GPU command buffer failed on segments %" PRIu64 "-%" PRIu64
                         "; redoing them (try %d)\n", X->s, X->s + X->ns - 1,
                         X->tries + 1);
                gpu_encode(X);
                continue;
            }
            if (r == -1)
                note("gapsieve: the GPU failed segments %" PRIu64 "-%" PRIu64
                     " %d times; searching them on the CPU\n",
                     X->s, X->s + X->ns - 1, X->tries);
            for (uint64_t g = 0; g < X->ns; g++) {
                int need[2] = { X->need[0], X->need[1] };
                int complete = 1;
                nh = 0;
                if (r < 0) {
                    complete = have_x && cpu_segment(&x, X->s + g, need, &hits, &nh, &ch);
                } else {
                    segment_strip(X->s + g, need, &hits, &nh, &ch);
                    if (gpu)
                        gpu_results(X, g, &hits, &nh, &ch);
                }
                if (complete)
                    segment_done(X->s + g, need, hits, nh);
            }
            head = (head + 1) % g_gdepth;
            nfly--;
        }
    }
    if (g_prof && g_nprof)
        note("gapsieve: GPU profile, %d-segment units, %d rounds then the tail, "
             "%d in flight: %" PRIu64 " segments; per segment %.3f ms sieving, "
             "%.3f ms screening, %.3f ms of span (GPU busy %.1f%%), %.0f candidates\n",
             g_gunit, g_gsplit, g_gdepth, g_nprof, 1e3 * g_tsieve / g_nprof,
             1e3 * g_tpass / g_nprof, 1e3 * (g_tspan1 - g_tspan0) / g_nprof,
             100 * (g_tsieve + g_tpass) / (g_tspan1 - g_tspan0),
             (double)g_ncand / g_nprof);
    free(hits);
    if (have_x)
        scratch_free(&x);
    __atomic_sub_fetch(&g_active, 1, __ATOMIC_RELEASE);
    return NULL;
}
