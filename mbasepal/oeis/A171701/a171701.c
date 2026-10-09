/* A171701: least number that is a 2-digit palindrome in at least n bases.
 *
 * "aa" in base b is v = a*(b+1) with 1 <= a <= b-1.  Writing m = b+1,
 * v is "aa" in base m-1 iff a | v and v/a >= a+2, i.e. a*(a+2) <= v.
 * So count(v) = #{ divisors a of v : a*(a+2) <= v }, and a(n) is the
 * first v with count(v) >= n.
 *
 * Segmented, multithreaded divisor sieve over v < V; prints b-file lines
 * "n a(n)" on stdout as new records appear.
 *
 *   ./a171701 [-V limit] [-N maxterms] [-t threads]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <time.h>

typedef uint64_t u64;
typedef uint16_t u16;

typedef struct { u64 lo, hi; u16 *cnt; } Job;

static void *work(void *p){
    Job *J = p;
    u64 lo = J->lo, hi = J->hi;
    memset(J->cnt, 0, (size_t)(hi - lo) * sizeof(u16));
    for (u64 a = 1; a * (a + 2) < hi; a++){
        u64 v = a * (a + 2);
        if (v < lo) v = (lo + a - 1) / a * a;
        for (; v < hi; v += a) J->cnt[v - lo]++;
    }
    return NULL;
}

static double now_s(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv){
    u64 V = 10000000000ULL, N = 0, seg = 1ULL << 27;
    int T = (int)sysconf(_SC_NPROCESSORS_ONLN), opt;
    while ((opt = getopt(argc, argv, "V:N:t:")) != -1){
        switch (opt){
        case 'V': V = strtoull(optarg, NULL, 10); break;
        case 'N': N = strtoull(optarg, NULL, 10); break;
        case 't': T = atoi(optarg); break;
        default:
            fprintf(stderr, "usage: %s [-V limit] [-N maxterms] [-t threads]\n", argv[0]);
            return 2;
        }
    }
    if (T < 1) T = 1;
    if (T > 256) T = 256;
    u16 *cnt = malloc((size_t)seg * sizeof(u16));
    if (!cnt){ fprintf(stderr, "out of memory\n"); return 1; }
    Job jobs[256];
    pthread_t tid[256];
    u64 maxc = 0;
    double t0 = now_s();
    for (u64 lo = 1; lo < V; lo += seg){
        u64 hi = lo + seg < V ? lo + seg : V;
        u64 span = hi - lo, per = (span + T - 1) / T;
        int nj = 0;
        for (u64 s = lo; s < hi; s += per, nj++){
            jobs[nj].lo = s;
            jobs[nj].hi = s + per < hi ? s + per : hi;
            jobs[nj].cnt = cnt + (s - lo);
            pthread_create(&tid[nj], NULL, work, &jobs[nj]);
        }
        for (int j = 0; j < nj; j++) pthread_join(tid[j], NULL);
        for (u64 v = lo; v < hi; v++){
            u64 c = cnt[v - lo];
            if (c > maxc){
                for (u64 n = maxc + 1; n <= c; n++){
                    printf("%llu %llu\n", (unsigned long long)n, (unsigned long long)v);
                    if (N && n >= N){ fflush(stdout); goto out; }
                }
                maxc = c;
            }
        }
        fflush(stdout);
        fprintf(stderr, "\r  sieved to %.3e  (%llu terms, %.1fs)   ",
                (double)hi, (unsigned long long)maxc, now_s() - t0);
    }
out:
    fprintf(stderr, "\ndone: %llu terms, all v < %llu sieved, %.1fs\n",
            (unsigned long long)maxc, (unsigned long long)V, now_s() - t0);
    free(cnt);
    return 0;
}
