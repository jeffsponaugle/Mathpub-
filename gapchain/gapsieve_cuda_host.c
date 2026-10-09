/*
 * gapsieve_cuda_host.c -- gapsieve with a CUDA GPU search, for Linux machines
 * with an NVIDIA GPU (written for the DGX Spark).  Like gapsieve_gpu.m on a
 * Mac, this file is the whole program apart from the kernels: it includes
 * gapsieve.c, which declares the three functions below when GAPSIEVE_GPU is
 * defined, and adds one worker thread that feeds the GPU from the same queue
 * of segments as the CPU threads, so checkpoints, --first, statistics and
 * the status line are unchanged.  gapsieve_cuda.cu holds the kernels.
 *
 * The GPU takes a unit of GPU_UNIT consecutive segments at a time (fewer for
 * short chains, whose lists would be huge), sieves every class of every
 * segment for each shape, and screens every candidate with the base-2 strong
 * test; the CPU then proves the rare survivors with chain_verify and
 * completes the unit's segments one by one.  Units run one after another on
 * one stream, the next always queued.  A failed unit is redone, and after
 * GPU_RETRIES failures its segments are searched on the CPU; a candidate list
 * that overflows is enlarged and the unit redone.
 *
 * Tuning: GAPSIEVE_GPU_UNIT, GAPSIEVE_GPU_SPLIT, GAPSIEVE_GPU_DEPTH and
 * GAPSIEVE_GPU_DENSE override the defaults below.  Test hooks:
 * GAPSIEVE_GPU_FAIL=N treats every Nth unit as failed, and GAPSIEVE_GPU_CAP=N
 * starts with N-entry candidate lists.
 *
 * Build: make CUDA=1 (see the Makefile).
 */
#define GAPSIEVE_GPU 1
#include "gapsieve.c"
#include "gapsieve_cuda.h"

#define GPU_UNIT 16             /* segments per unit */
#define GPU_SPLIT 3             /* screening rounds before the tail kernel */
#define GPU_DEPTH 2             /* units in flight */
#define GPU_DENSE 3             /* pattern groups sieved before compacting */
#define GPU_MAXUNIT 64
#define GPU_LIST (1u << 23)     /* about the most entries a unit's lists need */
#define GPU_RETRIES 3

struct gpu_slot {
    uint32_t cap;
    uint64_t s, ns;             /* segments s .. s+ns-1 */
    u128 base;                  /* segment s's base */
    int need[2], tries, queued;
};

static struct gpu_slot g_slot[CU_MAXSLOT];
static int g_gunit = GPU_UNIT, g_gsplit = GPU_SPLIT, g_gdepth = GPU_DEPTH;
static uint32_t g_gdense = GPU_DENSE;
static char g_gpuname[128];
static uint64_t g_failevery, g_ncb;
static uint32_t *g_bq;          /* scratch: each segment's base mod each prime */

static const char *gpu_name(void)
{
    return g_gpuname;
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
    if (g_nps != g_nq)
        return "the GPU sieve needs --depth equal to --presieve";
    if (cu_init(g_gpuname, sizeof g_gpuname))
        return "no CUDA device";

    /* the pattern groups; each shape's patterns back to back */
    uint32_t *gi = calloc((size_t)g_ngroups * 4, sizeof *gi);
    uint32_t *crt = malloc(g_nps * sizeof *crt);
    if (!gi || !crt)
        return "out of memory";
    size_t words = 0;
    for (uint32_t g = 0; g < g_ngroups; g++) {
        gi[4 * g] = g_groups[g].P;
        gi[4 * g + 1] = g_groups[g].first;
        gi[4 * g + 2] = g_groups[g].n;
        gi[4 * g + 3] = (uint32_t)words;
        words += (size_t)(g_groups[g].P + KBLOCKS + 128 + 63) / 64;
    }
    for (uint32_t i = 0; i < g_nps; i++)
        crt[i] = (uint32_t)g_crt[i];
    int bad = cu_tables(g_ngroups, gi, g_nps, g_q, g_qinv, crt);
    free(gi);
    free(crt);
    if (bad)
        return "could not load the sieve tables onto the GPU";
    uint64_t *pat = malloc(words * 8);
    if (!pat)
        return "out of memory";
    for (int d = 0; d < 2; d++) {
        if (!g_sh[d].on)
            continue;
        for (uint32_t g = 0, w = 0; g < g_ngroups; g++) {
            uint32_t n = (uint32_t)((g_groups[g].P + KBLOCKS + 128 + 63) / 64);
            memcpy(pat + w, g_groups[g].pat[d], (size_t)n * 8);
            w += n;
        }
        if (cu_shape(d, (const uint32_t *)pat, words * 2, g_sh[d].res, g_sh[d].A,
                     g_sh[d].offs, (uint32_t)g_k)) {
            free(pat);
            return "could not load the patterns onto the GPU";
        }
    }
    free(pat);

    g_gunit = gpu_env("GAPSIEVE_GPU_UNIT", GPU_UNIT, 1, GPU_MAXUNIT);
    g_gsplit = gpu_env("GAPSIEVE_GPU_SPLIT", GPU_SPLIT, 0, 1000);
    g_gdepth = gpu_env("GAPSIEVE_GPU_DEPTH", GPU_DEPTH, 1, CU_MAXSLOT);
    g_gdense = (uint32_t)gpu_env("GAPSIEVE_GPU_DENSE", GPU_DENSE, 1, 1000);
    if (g_gdense > g_ngroups)
        g_gdense = g_ngroups;
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
    e = getenv("GAPSIEVE_GPU_FAIL");
    g_failevery = e ? strtoull(e, NULL, 10) : 0;
    g_bq = malloc((size_t)g_gunit * g_nps * sizeof *g_bq);
    if (!g_bq)
        return "out of memory";
    for (int k = 0; k < g_gdepth; k++) {
        if (cu_slot(k, (uint32_t)g_gunit, cap ? cap : 1))
            return "out of GPU memory";
        g_slot[k].cap = cap ? cap : 1;
    }
    return NULL;
}

/* queue one unit: for each shape it needs, the offsets, the sieve and the
   screening, with no round trip to the host in between */
static int gpu_encode(int k)
{
    struct gpu_slot *X = &g_slot[k];
    u128 room = g_hi - X->base;                 /* base <= g_hi here */
    uint64_t span = X->ns * g_segsize - 1;
    uint64_t lim = room > span ? span : (uint64_t)room;
    for (uint64_t g = 0; g < X->ns; g++)
        for (uint32_t i = 0; i < g_nps; i++)
            g_bq[g * g_nps + i] = mod_small(X->base + (u128)g * g_segsize, g_q[i]);
    struct cu_params p[2];
    const struct cu_params *pp[2] = { NULL, NULL };
    for (int d = 0; d < 2; d++) {
        if (!X->need[d])
            continue;
        p[d] = (struct cu_params){
            (uint64_t)X->base, (uint64_t)(X->base >> 64), lim, g_segsize,
            (uint32_t)g_W, mod_small(X->base, (uint32_t)g_W), g_sh[d].A, g_ngroups,
            g_nps, (uint32_t)g_k, X->cap, (uint32_t)X->ns, g_gdense };
        pp[d] = &p[d];
    }
    X->queued = 1;
    return cu_submit(k, g_bq, (uint32_t)X->ns, g_nps, pp, g_gsplit);
}

/* wait for a slot: 0 done, -1 failed, 1 a candidate list overflowed */
static int gpu_wait(int k, uint32_t stat[2], uint32_t nfin[2])
{
    struct gpu_slot *X = &g_slot[k];
    int failed = cu_wait(k, stat, nfin) != 0;
    X->queued = 0;
    if (g_failevery && ++g_ncb % g_failevery == 0)
        failed = 1;
    if (failed)
        return -1;
    for (int d = 0; d < 2; d++)
        if (X->need[d] && stat[d] > X->cap)
            return 1;
    return 0;
}

/* prove the GPU's survivors in the unit's segment g on the CPU, as
   cpu_segment does its own; lst[d] holds nfin[d] of them */
static void gpu_results(int k, uint64_t g, const uint32_t nfin[2],
                        const uint64_t *lst[2], struct hit **hits, size_t *nh,
                        size_t *ch)
{
    struct gpu_slot *X = &g_slot[k];
    for (int d = 0; d < 2; d++) {
        if (!X->need[d] || !nfin[d])
            continue;
        const uint64_t *off = lst[d];
        for (uint32_t i = 0; i < nfin[d]; i++) {
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
        while (more && nfly < g_gdepth) {
            int k = (head + nfly) % g_gdepth;
            struct gpu_slot *X = &g_slot[k];
            if (!claim_segments(&X->s, &X->ns, (uint64_t)g_gunit, X->need)) {
                more = 0;
                break;
            }
            X->base = g_sieve_lo + (u128)X->s * g_segsize;
            X->tries = 0;
            X->queued = 0;
            if (X->base <= g_hi && gpu_encode(k) != 0)
                X->queued = -1;         /* submit failed: counts as a failure */
            nfly++;
        }
        if (!nfly)
            break;

        int k = head;
        struct gpu_slot *X = &g_slot[k];
        uint32_t stat[2] = { 0, 0 }, nfin[2] = { 0, 0 };
        int gpu = X->queued != 0;
        int r = X->queued < 0 ? -1 : gpu ? gpu_wait(k, stat, nfin) : 0;
        if (r == 1) {
            uint32_t need = X->cap;
            for (int d = 0; d < 2; d++)
                if (X->need[d] && stat[d] > need)
                    need = stat[d];
            uint64_t cap = (uint64_t)need + need / 4 + 1024;
            if (cap > UINT32_MAX / 2 || cu_grow(k, (uint32_t)cap)) {
                note("gapsieve: GPU candidate lists cannot grow to %" PRIu64
                     "; searching segments %" PRIu64 "-%" PRIu64 " on the CPU\n",
                     cap, X->s, X->s + X->ns - 1);
                r = -2;
            } else {
                X->cap = (uint32_t)cap;
                log_note("GPU candidate lists enlarged to %" PRIu64 "\n", cap);
                if (gpu_encode(k) != 0)
                    X->queued = -1;
                continue;
            }
        }
        if (r == -1 && ++X->tries <= GPU_RETRIES) {
            log_note("GPU work failed on segments %" PRIu64 "-%" PRIu64
                     "; redoing them (try %d)\n", X->s, X->s + X->ns - 1, X->tries + 1);
            if (gpu_encode(k) != 0)
                X->queued = -1;
            continue;
        }
        if (r == -1)
            note("gapsieve: the GPU failed segments %" PRIu64 "-%" PRIu64
                 " %d times; searching them on the CPU\n", X->s, X->s + X->ns - 1,
                 X->tries);
        const uint64_t *lst[2] = { NULL, NULL };
        for (int d = 0; d < 2 && r >= 0 && gpu; d++)
            if (X->need[d] && nfin[d] && !(lst[d] = cu_list(k, d, nfin[d]))) {
                note("gapsieve: could not read the GPU's survivors\n");
                _exit(1);
            }
        for (uint64_t g = 0; g < X->ns; g++) {
            int need[2] = { X->need[0], X->need[1] };
            int complete = 1;
            nh = 0;
            if (r < 0) {
                complete = have_x && cpu_segment(&x, X->s + g, need, &hits, &nh, &ch);
            } else {
                segment_strip(X->s + g, need, &hits, &nh, &ch);
                if (gpu)
                    gpu_results(k, g, nfin, lst, &hits, &nh, &ch);
            }
            if (complete)
                segment_done(X->s + g, need, hits, nh);
        }
        head = (head + 1) % g_gdepth;
        nfly--;
    }
    free(hits);
    if (have_x)
        scratch_free(&x);
    __atomic_sub_fetch(&g_active, 1, __ATOMIC_RELEASE);
    return NULL;
}
