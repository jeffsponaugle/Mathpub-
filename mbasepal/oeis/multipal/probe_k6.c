/*
 * probe_k6.c — structured probe for A171705 (first number that is a
 * 6-digit palindrome in at least n bases).
 *
 * Scans v = c * L^5 (1 <= c <= CMAX, v <= VMAX) and counts every base b in
 * which v is a 6-digit palindrome. An even-length palindrome in base b is
 * divisible by b+1, so only divisors m = b+1 of v with b^5 <= v < b^6 are
 * tested; they are generated from the factorizations of c and L. Prints
 * every v with at least NMIN bases.
 *
 *   ./probe_k6 [-c CMAX] [-n NMIN] [-t threads]   (VMAX = 40968^5)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <math.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;

static u128 VMAX;
static int CMAX = 64, NMIN = 4, NT = 8;
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static u128 best = 0;

static void u128str(u128 v, char *buf){
    char t[64]; int n = 0;
    if (!v){ strcpy(buf, "0"); return; }
    while (v){ t[n++] = '0' + (int)(v % 10); v /= 10; }
    for (int i = 0; i < n; i++) buf[i] = t[n - 1 - i];
    buf[n] = 0;
}

static int factor(u64 x, u64 *p, int *e){
    int k = 0;
    for (u64 q = 2; q * q <= x; q++){
        if (x % q == 0){ p[k] = q; e[k] = 0; while (x % q == 0){ x /= q; e[k]++; } k++; }
    }
    if (x > 1){ p[k] = x; e[k] = 1; k++; }
    return k;
}

static int is_pal6(u128 v, u64 b){
    unsigned d[8]; int n = 0;
    while (v && n < 8){ d[n++] = (unsigned)(v % b); v /= b; }
    if (v || n != 6) return 0;
    return d[0] == d[5] && d[1] == d[4] && d[2] == d[3];
}

typedef struct { u64 p[32]; int e[32]; int k; u128 v; u64 mlo, mhi; int cnt; u64 bases[64]; } Ctx;

static void gen(Ctx *C, int i, u128 m){
    if (m > C->mhi) return;
    if (i == C->k){
        if (m >= C->mlo && is_pal6(C->v, (u64)(m - 1))){
            if (C->cnt < 64) C->bases[C->cnt] = (u64)(m - 1);
            C->cnt++;
        }
        return;
    }
    u128 pw = 1;
    for (int j = 0; j <= C->e[i]; j++){
        if (m * pw > C->mhi) break;
        gen(C, i + 1, m * pw);
        pw *= C->p[i];
    }
}

static int count_bases(u128 v, u64 c, u64 L, Ctx *C){
    u64 pc[16], pl[16]; int ec[16], el[16];
    int kc = factor(c, pc, ec), kl = factor(L, pl, el);
    C->k = 0;
    for (int i = 0; i < kl; i++){ C->p[C->k] = pl[i]; C->e[C->k] = 5 * el[i]; C->k++; }
    for (int i = 0; i < kc; i++){
        int f = -1;
        for (int j = 0; j < C->k; j++) if (C->p[j] == pc[i]) f = j;
        if (f >= 0) C->e[f] += ec[i];
        else { C->p[C->k] = pc[i]; C->e[C->k] = ec[i]; C->k++; }
    }
    /* b^5 <= v < b^6  <=>  b in (v^(1/6), v^(1/5)] ; m = b + 1 */
    long double vd = (long double)v;
    u64 blo = (u64)powl(vd, 1.0L / 6) - 2, bhi = (u64)powl(vd, 0.2L) + 2;
    if (blo < 2) blo = 2;
    C->v = v; C->mlo = blo + 1; C->mhi = bhi + 1; C->cnt = 0;
    gen(C, 0, 1);
    return C->cnt;
}

static void *work(void *arg){
    int tid = (int)(intptr_t)arg;
    Ctx C;
    for (int c = 1; c <= CMAX; c++){
        for (u64 L = 2 + (u64)tid; ; L += (u64)NT){
            u128 v = (u128)c * L * L * L * L * L;
            if (v > VMAX) break;
            int n = count_bases(v, (u64)c, L, &C);
            if (n >= NMIN){
                char buf[64]; u128str(v, buf);
                pthread_mutex_lock(&mu);
                printf("v=%s = %d*%llu^5  bases=%d:", buf, c, (unsigned long long)L, n);
                for (int i = 0; i < n && i < 64; i++) printf(" %llu", (unsigned long long)C.bases[i]);
                printf("\n");
                if (!best || v < best) best = v;
                pthread_mutex_unlock(&mu);
            }
        }
    }
    return NULL;
}

int main(int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, "c:n:t:")) != -1){
        if (opt == 'c') CMAX = atoi(optarg);
        else if (opt == 'n') NMIN = atoi(optarg);
        else if (opt == 't') NT = atoi(optarg);
    }
    VMAX = (u128)40968 * 40968 * 40968 * 40968 * 40968;
    pthread_t th[64];
    for (int i = 0; i < NT; i++) pthread_create(&th[i], NULL, work, (void *)(intptr_t)i);
    for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);
    char buf[64]; u128str(best, buf);
    printf("# smallest v with >= %d bases among c*L^5, c <= %d, v <= 40968^5: %s\n", NMIN, CMAX, best ? buf : "none");
    return 0;
}
