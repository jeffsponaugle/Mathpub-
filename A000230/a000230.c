/*
 * a000230 - compute terms of OEIS A000230
 *
 *   a(0) = 2; for n >= 1, a(n) = smallest prime p such that the gap
 *   between p and the next prime is exactly 2n.
 *
 * Usage: a000230 -n <terms> [-t <threads>] [-v]
 *
 * Strategy: sieve primes over expanding blocks [lo, hi) using
 * libprimesieve. Each block is split into <threads> contiguous chunks;
 * a worker walks the primes of its chunk with a primesieve_iterator and
 * records the first prime in the chunk that begins each gap of interest.
 * A gap "belongs" to the prime that starts it, and every prime lies in
 * exactly one chunk, so no gap is missed or double-counted at chunk or
 * block boundaries. Results are merged in ascending chunk order, so the
 * first value recorded for a gap is the smallest such prime.
 */

#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <primesieve.h>

typedef struct {
    uint64_t lo, hi;      /* chunk range [lo, hi) */
    uint64_t *local;      /* local[n] = first prime in chunk starting a gap of 2n, 0 if none */
    const uint64_t *found; /* global results (read-only): skip gaps already found */
    int nterms;
} job_t;

static void *worker(void *arg)
{
    job_t *job = (job_t *)arg;
    int wanted = job->nterms;

    primesieve_iterator it;
    primesieve_init(&it);
    /* next call to next_prime() returns the first prime >= lo */
    primesieve_jump_to(&it, job->lo, job->hi);

    uint64_t p = primesieve_next_prime(&it);
    while (p < job->hi) {
        uint64_t q = primesieve_next_prime(&it);
        uint64_t gap = q - p;
        if ((gap & 1) == 0) {
            uint64_t n = gap >> 1;
            if (n < (uint64_t)wanted && job->local[n] == 0 && job->found[n] == 0)
                job->local[n] = p;
        }
        p = q;
    }

    primesieve_free_iterator(&it);
    return NULL;
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s -n <terms> [-t <threads>] [-v]\n"
            "  -n <terms>    number of terms to produce (a(0) .. a(terms-1))\n"
            "  -t <threads>  worker threads (default: number of CPUs)\n"
            "  -v            print progress to stderr\n",
            prog);
}

int main(int argc, char **argv)
{
    long nterms = 0;
    long nthreads = 0;
    int verbose = 0;
    int opt;

    while ((opt = getopt(argc, argv, "n:t:vh")) != -1) {
        switch (opt) {
        case 'n': nterms = strtol(optarg, NULL, 10); break;
        case 't': nthreads = strtol(optarg, NULL, 10); break;
        case 'v': verbose = 1; break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
        }
    }
    if (nterms <= 0) {
        usage(argv[0]);
        return 1;
    }
    if (nthreads <= 0)
        nthreads = sysconf(_SC_NPROCESSORS_ONLN);
    if (nthreads < 1)
        nthreads = 1;
    if (nthreads > 1024)
        nthreads = 1024;

    /* found[n] = a(n) once discovered; 0 = still missing */
    uint64_t *found = calloc((size_t)nterms, sizeof(uint64_t));
    job_t *jobs = calloc((size_t)nthreads, sizeof(job_t));
    pthread_t *tids = calloc((size_t)nthreads, sizeof(pthread_t));
    if (!found || !jobs || !tids) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    for (long i = 0; i < nthreads; i++) {
        jobs[i].local = calloc((size_t)nterms, sizeof(uint64_t));
        if (!jobs[i].local) {
            fprintf(stderr, "out of memory\n");
            return 1;
        }
        jobs[i].found = found;
        jobs[i].nterms = (int)nterms;
    }

    long remaining = nterms;
    found[0] = 2; /* a(0) = 2 by definition (gap of 1 between 2 and 3) */
    remaining--;

    uint64_t lo = 2;
    uint64_t hi = 10 * 1000 * 1000; /* first block; doubles each round */

    while (remaining > 0) {
        uint64_t span = hi - lo;
        uint64_t chunk = span / (uint64_t)nthreads;
        if (chunk == 0)
            chunk = 1;

        long used = 0;
        uint64_t c = lo;
        for (long i = 0; i < nthreads && c < hi; i++) {
            jobs[i].lo = c;
            jobs[i].hi = (i == nthreads - 1) ? hi : c + chunk;
            if (jobs[i].hi > hi)
                jobs[i].hi = hi;
            memset(jobs[i].local, 0, (size_t)nterms * sizeof(uint64_t));
            if (pthread_create(&tids[i], NULL, worker, &jobs[i]) != 0) {
                fprintf(stderr, "pthread_create failed\n");
                return 1;
            }
            c = jobs[i].hi;
            used++;
        }
        for (long i = 0; i < used; i++)
            pthread_join(tids[i], NULL);

        /* merge in ascending chunk order: earliest chunk = smallest primes */
        for (long i = 0; i < used; i++) {
            for (long n = 1; n < nterms; n++) {
                if (found[n] == 0 && jobs[i].local[n] != 0) {
                    found[n] = jobs[i].local[n];
                    remaining--;
                }
            }
        }

        if (verbose)
            fprintf(stderr, "sieved up to %" PRIu64 ": %ld/%ld terms found\n",
                    hi, nterms - remaining, nterms);

        lo = hi;
        if (hi > UINT64_MAX / 2) {
            fprintf(stderr, "search limit of 2^64 reached with %ld terms missing\n",
                    remaining);
            return 1;
        }
        hi *= 2;
    }

    for (long n = 0; n < nterms; n++)
        printf("%ld %" PRIu64 "\n", n, found[n]);

    for (long i = 0; i < nthreads; i++)
        free(jobs[i].local);
    free(jobs);
    free(tids);
    free(found);
    return 0;
}
