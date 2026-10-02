/*
 * gpu_leaf.h -- CUDA leaf for psph (postage stamp exhaustive search), C-linkage API.
 *
 * psph.c (built with -DGPU_LEAF) hands every surviving (k-1)-prefix to gpu_submit_prefix() instead of
 * running its own a_k candidate loop.  Prefixes are batched per CPU thread; a full batch runs
 *   kernel 1 (filter: Challis X test, deep-hole targets, reach target, top-block targets, TGT) and
 *   kernel 2 (exact streaming full check of the survivors)
 * on the thread's own CUDA stream, and every basis that passes the full check is handed back through the
 * report callback (psph.c re-verifies it with an independent exact DP in report_solution).
 *
 * Cell type must match psph.c: uint8_t (default) or uint16_t (-DWIDE) -- both sides see the same macro.
 */
#ifndef GPU_LEAF_H
#define GPU_LEAF_H
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIDE
typedef uint16_t gpu_cell_t;
#else
typedef uint8_t gpu_cell_t;
#endif

/* called (on the submitting CPU thread, from inside gpu_submit_prefix / gpu_flush) for every a_k that
   passed the exact full check; a[1..k-1] is the prefix, ak the new last element */
typedef void (*gpu_report_fn)(int thread_id, const int64_t *a, int64_t ak, void *user);

typedef struct {
    unsigned long long cand, xpass, rpass, s2pass, fullfail, direct, probes, sols;
    unsigned long long fhist[18];              /* full-check failure position, blocks below the top block */
    unsigned long long batches, prefixes, tab_bytes, fc_cells, filter_runs;
    double t_h2d, t_filter, t_full, t_host;    /* seconds: GPU events (h2d, kernel 1, kernel 2), host wall of the pipeline */
} gpu_stats_t;

/* one call per search (run_search); allocations persist across calls and are released by gpu_shutdown().
   ndeep_cap: maximum number of deep-hole targets accepted per prefix.  Returns 0 on success. */
int  gpu_init(int h, int k, int64_t tgt, int nthreads, int ndeep_cap, gpu_report_fn report, void *user);

/* append one (k-1)-prefix to thread tid's batch.  T must be valid for indices 0..E (E = min(TGT, h*a_{k-1})),
   candidates are a_k in [lo, hi] (lo <= hi), a[1..k-1] is the prefix, Z the prefix part of the Challis
   target, reach[0..h] the reach sequence r_m = n_m(A_{k-1}), deep_y[0..ndeep) the deep-hole remainders
   (positions y in [1,hi) with large T[y]; any set is exact, only y < a_k are used per candidate). */
void gpu_submit_prefix(int tid, const gpu_cell_t *T, int64_t E, int64_t lo, int64_t hi,
                       const int64_t *a, int64_t Z, const int64_t *reach,
                       const int64_t *deep_y, int ndeep);

/* run the pending batch of thread tid (if any) and report its solutions */
void gpu_flush(int tid);

/* statistics of thread tid accumulated since gpu_init (safe to call from the progress thread: counters are
   updated once per finished batch, so they may lag by a batch) */
void gpu_thread_stats(int tid, gpu_stats_t *out);

/* flush all threads, synchronize the device; then gpu_print_stats prints the GPU summary to f */
void gpu_finish(void);
void gpu_print_stats(FILE *f);

/* release all device/host memory (end of program) */
void gpu_shutdown(void);

#ifdef __cplusplus
}
#endif
#endif
