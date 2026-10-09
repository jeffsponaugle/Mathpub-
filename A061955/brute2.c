// brute2.c -- direct reference for the whole family of OEIS "concatenate 1..n" sequences.
// Computes V(n) mod n straight from the definition, word by word (O(n log n)).
//
//   brute2 FAM B LO HI             print "n residue" for every n in [LO, HI)
//   brute2 FAM B -v N [THREADS]    print "N residue" for one (large) N, using THREADS
//                                  threads (default: all CPUs); residue 0 means N is a term.
//                                  Splits 1..N into ranges whose digits are incremented in
//                                  place (no division per word); the LO HI mode uses the
//                                  plain per-word conversion, and the two are cross-checked.
//
// FAM: Ln/Rn = normal digits (A029471+ / A029447+),
//      Lk/Rk = reversed digits, all zeros kept   (A061955+ / A029495+),
//      Ld/Rd = reversed digits, least significant zeros of k dropped (A029519+ / A061931+).
// L = left concatenation (n ... 2 1), R = right concatenation (1 2 ... n).
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef unsigned long long ull;

static int g_b, g_left, g_rev, g_drop;

// Word for k: its value mod n and its length in digits.
static u64 word(u64 k, u64 n, int *len)
{
    int dig[70], d = 0;
    for (u64 t = k; t; t /= g_b) dig[d++] = (int)(t % g_b);   // dig[0] = least significant
    u64 val = 0;
    if (!g_rev) {
        val = k % n;
        *len = d;
    } else {                                // word = dig[0] dig[1] ... dig[d-1]
        int s = 0;
        if (g_drop) while (s < d - 1 && dig[s] == 0) s++;
        for (int i = s; i < d; i++) val = (u64)(((u128)val * g_b + dig[i]) % n);
        *len = d - s;
    }
    return val;
}

// Words lo..hi (inclusive) in concatenation order, as (V, X) = (value, b^length) mod n.
typedef struct { u64 n, lo, hi, V, X; const u64 *bpow; } part_t;

// Same as word(), but for consecutive k: the digits of k are kept in dig[0..*d-1] and
// incremented in place (with carry), so no division is needed per word.  n < 2^58,
// so every word value (< b^d <= b k) fits in 64 bits before the final reduction.
static u64 next_word(int *dig, int *d, u64 n, int *len)
{
    u64 val = 0;
    if (!g_rev) {
        for (int i = *d - 1; i >= 0; i--) val = val * g_b + dig[i];
        *len = *d;
    } else {
        int s = 0;
        if (g_drop) while (s < *d - 1 && dig[s] == 0) s++;
        for (int i = s; i < *d; i++) val = val * g_b + dig[i];
        *len = *d - s;
    }
    int i = 0;                                // k -> k+1
    while (i < *d && dig[i] == g_b - 1) dig[i++] = 0;
    if (i == *d) dig[(*d)++] = 1; else dig[i]++;
    return val % n;
}

static void *run_part(void *arg)
{
    part_t *p = arg;
    u64 n = p->n, V = 0, X = 1 % n;
    int dig[70], d = 0;
    for (u64 t = p->lo; t; t /= g_b) dig[d++] = (int)(t % g_b);
    for (u64 k = p->lo; k <= p->hi; k++) {
        int len;
        u64 val = next_word(dig, &d, n, &len), bl = p->bpow[len];
        if (g_left) V = (u64)(((u128)val * X + V) % n);   // w(k) goes to the left
        else V = (u64)(((u128)V * bl + val) % n);          // w(k) goes to the right
        X = (u64)((u128)X * bl % n);
    }
    p->V = V; p->X = X;
    return 0;
}

static u64 residue(u64 n, int nth)
{
    u64 bpow[70];
    bpow[0] = 1 % n;
    for (int i = 1; i < 70; i++) bpow[i] = (u64)((u128)bpow[i - 1] * g_b % n);
    if ((u64)nth > n) nth = (int)n;
    part_t *P = calloc(nth, sizeof *P);
    pthread_t *th = calloc(nth, sizeof *th);
    for (int i = 0; i < nth; i++) {
        P[i].n = n; P[i].bpow = bpow;
        P[i].lo = 1 + (u64)((u128)n * i / nth);
        P[i].hi = (u64)((u128)n * (i + 1) / nth);
        if (nth > 1) pthread_create(&th[i], 0, run_part, &P[i]);
        else run_part(&P[i]);
    }
    if (nth > 1) for (int i = 0; i < nth; i++) pthread_join(th[i], 0);
    u64 V = 0, X = 1 % n;                   // join the parts in string order
    for (int i = 0; i < nth; i++) {
        if (g_left) V = (u64)(((u128)P[i].V * X + V) % n);
        else V = (u64)(((u128)V * P[i].X + P[i].V) % n);
        X = (u64)((u128)X * P[i].X % n);
    }
    free(P); free(th);
    return V;
}

int main(int argc, char **argv)
{
    if (argc < 5 || strlen(argv[1]) != 2 || !strchr("LR", argv[1][0]) || !strchr("nkd", argv[1][1])) {
        fprintf(stderr, "usage: %s Ln|Rn|Lk|Rk|Ld|Rd BASE LO HI\n       %s FAM BASE -v N [THREADS]\n",
                argv[0], argv[0]);
        return 1;
    }
    g_left = argv[1][0] == 'L';
    g_rev = argv[1][1] != 'n';
    g_drop = argv[1][1] == 'd';
    g_b = atoi(argv[2]);
    if (!strcmp(argv[3], "-v")) {
        u64 n = strtoull(argv[4], 0, 0);
        int nth = argc > 5 ? atoi(argv[5]) : (int)sysconf(_SC_NPROCESSORS_ONLN);
        printf("%llu %llu\n", (ull)n, (ull)residue(n, nth < 1 ? 1 : nth));
        return 0;
    }
    u64 lo = strtoull(argv[3], 0, 0), hi = strtoull(argv[4], 0, 0);
    for (u64 n = lo ? lo : 1; n < hi; n++) {  // the plainest form: word() for each k
        u64 V = 0, X = 1 % n, bpow[70];
        bpow[0] = 1 % n;
        for (int i = 1; i < 70; i++) bpow[i] = (u64)((u128)bpow[i - 1] * g_b % n);
        for (u64 k = 1; k <= n; k++) {
            int len;
            u64 val = word(k, n, &len);
            if (g_left) V = (u64)(((u128)val * X + V) % n);
            else V = (u64)(((u128)V * bpow[len] + val) % n);
            X = (u64)((u128)X * bpow[len] % n);
        }
        printf("%llu %llu\n", (ull)n, (ull)V);
    }
    return 0;
}
