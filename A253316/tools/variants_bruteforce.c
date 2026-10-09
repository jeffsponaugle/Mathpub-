/*
 * variants_bruteforce.c -- brute-force counts of 2n x 2n binary grids with n
 * ones in every row and column and no three equal consecutive entries in any
 * row or column, split by whether rows / columns are distinct (n <= 4; n = 4
 * takes a few seconds).
 *
 *   M = all such grids (no distinctness rule)
 *   a = rows distinct and columns distinct   (A253316)
 *   X = columns distinct (= rows distinct, by transposition)
 *
 * usage: cc -O2 -o variants_bruteforce variants_bruteforce.c && ./variants_bruteforce 4
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static int N, n, nR;
static uint32_t R[300], rows[16];
static int cnt[16];
static uint64_t tot[2][2];          /* [rows distinct][columns distinct] */

static int valid(uint32_t x)
{
    if (__builtin_popcount(x) != n) return 0;
    for (int i = 0; i + 2 < N; i++) {
        int t = (x >> i) & 7;
        if (t == 0 || t == 7) return 0;
    }
    return 1;
}

static void dfs(int k)
{
    uint32_t mask = (1u << N) - 1;
    if (k == N) {
        int dr = 1, dc = 1;
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++) if (rows[i] == rows[j]) dr = 0;
        uint32_t col[16];
        for (int j = 0; j < N; j++) {
            col[j] = 0;
            for (int i = 0; i < N; i++) col[j] |= ((rows[i] >> j) & 1u) << i;
        }
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++) if (col[i] == col[j]) dc = 0;
        tot[dr][dc]++;
        return;
    }
    for (int r = 0; r < nR; r++) {
        uint32_t x = R[r];
        if (k >= 2 && (~(rows[k - 2] ^ rows[k - 1]) & ~(rows[k - 1] ^ x) & mask)) continue;
        int ok = 1;
        for (int j = 0; j < N; j++) {
            int c = cnt[j] + (int)((x >> j) & 1);
            if (c > n || c + (N - 1 - k) < n) { ok = 0; break; }
        }
        if (!ok) continue;
        for (int j = 0; j < N; j++) cnt[j] += (x >> j) & 1;
        rows[k] = x;
        dfs(k + 1);
        for (int j = 0; j < N; j++) cnt[j] -= (x >> j) & 1;
    }
}

int main(int argc, char **argv)
{
    n = argc > 1 ? atoi(argv[1]) : 4;
    if (n < 1 || n > 4) { fprintf(stderr, "n must be 1..4\n"); return 2; }
    N = 2 * n;
    for (uint32_t x = 0; x < (1u << N); x++) if (valid(x)) R[nR++] = x;
    dfs(0);
    unsigned long long M = tot[0][0] + tot[0][1] + tot[1][0] + tot[1][1];
    printf("n=%d  M=%llu  a=%llu  X=%llu  (rows distinct only: %llu, columns distinct only: %llu, neither: %llu)\n",
           n, M, (unsigned long long)tot[1][1], (unsigned long long)(tot[1][1] + tot[0][1]),
           (unsigned long long)tot[1][0], (unsigned long long)tot[0][1], (unsigned long long)tot[0][0]);
    return 0;
}
