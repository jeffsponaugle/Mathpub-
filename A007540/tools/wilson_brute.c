/*
 * wilson_brute.c - brute-force (p-1)! mod p^2 for p < 2^39, multithreaded.
 *
 * Residues are kept as x = a*p + b with 0 <= a,b < p, so everything fits in
 * unsigned __int128 and the only division is floor(b*i/p), done with a
 * double-precision estimate corrected exactly. Each thread multiplies a
 * contiguous range of i; partial products are combined at the end.
 *
 * usage: wilson_brute [-t threads] p [p ...]
 * prints: p  (p-1)! mod p^2   k   where (p-1)! = -1 + k*p (mod p^2), |k| <= p/2
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
typedef unsigned long long u64; typedef unsigned __int128 u128; typedef __int128 i128;
typedef struct { u64 a, b; } res_t;            /* a*p + b */
static u64 P; static double Pd;
static inline u64 divmod_p(u128 x, u64 *rem) {  /* x < 2^79: q = floor(x/P), rem = x mod P */
    u64 q = (u64)((double)x / Pd);
    i128 r = (i128)x - (i128)((u128)q * P);
    while (r >= (i128)P) { r -= P; q++; }
    while (r < 0)        { r += P; q--; }
    *rem = (u64)r; return q;
}
static inline res_t mul_small(res_t x, u64 i) {  /* (a*p+b)*i mod p^2, i < p */
    u64 d, c = divmod_p((u128)x.b * i, &d);      /* b*i = c*p + d */
    u64 f; divmod_p((u128)x.a * i + c, &f);      /* (a*i + c) mod p */
    return (res_t){ f, d };
}
static res_t mul_res(res_t x, res_t y) {         /* general product mod p^2 */
    u64 d, c = divmod_p((u128)x.b * y.b, &d);
    u128 t = (u128)x.a * y.b + (u128)x.b * y.a + c;   /* < 2*2^78 + 2^39 */
    u64 f; divmod_p(t % P, &f);  /* t % P is exact 128-bit mod (rare, cheap here) */
    return (res_t){ f, d };
}
typedef struct { u64 lo, hi; res_t out; } job_t;
static void *worker(void *arg) { job_t *j = arg; res_t r = {0, 1}; for (u64 i = j->lo; i < j->hi; i++) r = mul_small(r, i); j->out = r; return NULL; }
int main(int argc, char **argv) {
    int nt = 8, a = 1;
    if (argc > 2 && !strcmp(argv[1], "-t")) { nt = atoi(argv[2]); a = 3; }
    for (; a < argc; a++) {
        P = strtoull(argv[a], 0, 10); Pd = (double)P;
        if (P >= (1ULL << 39)) { fprintf(stderr, "p too large (need p < 2^39)\n"); return 1; }
        pthread_t th[256]; job_t jb[256]; u64 n = P - 1;      /* i = 1..p-1 */
        for (int t = 0; t < nt; t++) { jb[t].lo = 1 + n * t / nt; jb[t].hi = 1 + n * (t + 1) / nt; pthread_create(&th[t], 0, worker, &jb[t]); }
        res_t r = {0, 1};
        for (int t = 0; t < nt; t++) { pthread_join(th[t], 0); r = mul_res(r, jb[t].out); }
        u128 R = (u128)r.a * P + r.b;                          /* (p-1)! mod p^2 */
        /* R = -1 + k p  =>  k = (R+1)/p mod p  */
        u128 R1 = (R + 1) % ((u128)P * P); long long k = (long long)(R1 / P); if (k > (long long)(P / 2)) k -= (long long)P;
        unsigned long long hi = (unsigned long long)(R / 1000000000000000000ULL), lo = (unsigned long long)(R % 1000000000000000000ULL);
        if (hi) printf("%llu %llu%018llu  k=%lld\n", P, hi, lo, k); else printf("%llu %llu  k=%lld\n", P, lo, k);
        fflush(stdout);
    }
    return 0;
}
