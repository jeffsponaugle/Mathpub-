/*
 * verify.c - independent brute-force check for A399971 (small n).
 *
 * Shares no code with a399971.c and uses none of its shortcuts:
 *   - enumerates EVERY partition of the n X n grid into paths of >= 2 cells
 *     (cycles rejected with a relabelling component array);
 *   - decides compactness by trying every routing of the endpoint pairs
 *     (all simple paths, pairs in a fixed order), only pruning on plain
 *     reachability; a routing that leaves a cell empty makes it non-compact;
 *   - reduces by symmetry with Burnside's lemma (fixed-point counts) instead
 *     of canonical representatives, for solutions and for starting positions.
 *
 * Build: cc -O2 -o verify verify.c
 * Usage: verify n
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { UP = 1, RT = 2, DN = 4, LT = 8 };
static const int DR[4] = {-1, 0, 1, 0}, DC[4] = {0, 1, 0, -1}; /* UP RT DN LT */

static int n, cells;
static int dirs[64];     /* edge directions used at each cell */
static int comp[64];     /* component label of each cell */
static int partner[64];  /* other endpoint of the path, -1 for inner cells */

static long long all_covers, compact_lab, fixsol[8], fixpuz[8][16];

static void sym(int g, int r, int c, int *R, int *C)
{
    const int m = n - 1;
    const int t[8][2] = {{r, c}, {c, m - r}, {m - r, m - c}, {m - c, r},
                         {r, m - c}, {m - r, c}, {c, r}, {m - c, m - r}};
    *R = t[g][0];
    *C = t[g][1];
}

/* image under g of a direction bit, found by mapping a unit step */
static int symdir(int g, int d)
{
    int r0, c0, r1, c1;
    sym(g, 1, 1, &r0, &c0);
    sym(g, 1 + DR[d], 1 + DC[d], &r1, &c1);
    for (int e = 0; e < 4; e++)
        if (r1 - r0 == DR[e] && c1 - c0 == DC[e]) return e;
    abort();
}

/* ---- exhaustive routing ---- */
static int np, ps[32], pt[32];
static unsigned char used[64];
static long long routings;

static int reachable(int s, int t)
{
    int q[64], h = 0, e = 0;
    unsigned char seen[64] = {0};
    q[e++] = s;
    seen[s] = 1;
    while (h < e) {
        const int x = q[h++], r = x / n, c = x % n;
        for (int d = 0; d < 4; d++) {
            const int rr = r + DR[d], cc = c + DC[d];
            if (rr < 0 || rr >= n || cc < 0 || cc >= n) continue;
            const int y = rr * n + cc;
            if (y == t) return 1;
            if (!seen[y] && !used[y]) {
                seen[y] = 1;
                q[e++] = y;
            }
        }
    }
    return 0;
}

static int route_pair(int j);

/* extend pair j from x; returns 1 as soon as a routing with an empty cell is found */
static int extend(int j, int x)
{
    const int r = x / n, c = x % n;
    for (int d = 0; d < 4; d++) {
        const int rr = r + DR[d], cc = c + DC[d];
        if (rr < 0 || rr >= n || cc < 0 || cc >= n) continue;
        const int y = rr * n + cc;
        if (y == pt[j]) {
            if (route_pair(j + 1)) return 1;
        } else if (!used[y]) {
            used[y] = 1;
            const int f = extend(j, y);
            used[y] = 0;
            if (f) return 1;
        }
    }
    return 0;
}

static int route_pair(int j)
{
    if (j == np) {
        routings++;
        for (int i = 0; i < cells; i++)
            if (!used[i]) return 1;
        return 0;
    }
    for (int i = j; i < np; i++)
        if (!reachable(ps[i], pt[i])) return 0;
    return extend(j, ps[j]);
}

static void leaf(void)
{
    all_covers++;
    np = 0;
    memset(used, 0, sizeof used);
    for (int i = 0; i < cells; i++) partner[i] = -1;
    for (int i = 0; i < cells; i++) {
        if (__builtin_popcount(dirs[i]) != 1 || partner[i] >= 0) continue;
        int x = i, from = -1;
        for (;;) { /* walk to the other end */
            int nx = -1;
            for (int d = 0; d < 4; d++)
                if (dirs[x] & (1 << d)) {
                    const int y = (x / n + DR[d]) * n + x % n + DC[d];
                    if (y != from) nx = y;
                }
            if (nx < 0) break;
            from = x;
            x = nx;
        }
        partner[i] = x;
        partner[x] = i;
        ps[np] = i;
        pt[np] = x;
        np++;
        used[i] = used[x] = 1;
    }
    routings = 0;
    if (route_pair(0)) return;
    const long long m = routings; /* every routing is full: m solutions share this puzzle */
    compact_lab++;
    for (int g = 0; g < 8; g++) {
        int same_sol = 1, same_puz = 1;
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++) {
                int R, C, img = 0;
                sym(g, r, c, &R, &C);
                for (int d = 0; d < 4; d++)
                    if (dirs[r * n + c] & (1 << d)) img |= 1 << symdir(g, d);
                if (img != dirs[R * n + C]) same_sol = 0;
                const int p = partner[r * n + c];
                int pr = -1, pc = -1;
                if (p >= 0) sym(g, p / n, p % n, &pr, &pc);
                if ((p < 0 ? -1 : pr * n + pc) != partner[R * n + C]) same_puz = 0;
            }
        fixsol[g] += same_sol;
        if (same_puz) fixpuz[g][m < 15 ? m : 15]++;
    }
}

static void enumerate(int i)
{
    if (i == cells) {
        leaf();
        return;
    }
    const int r = i / n, c = i % n;
    for (int o = 0; o < 4; o++) {
        const int rt = o & 1, dn = o >> 1;
        if ((rt && c + 1 >= n) || (dn && r + 1 >= n)) continue;
        const int deg = __builtin_popcount(dirs[i]) + rt + dn;
        if (deg < 1 || deg > 2) continue;
        if (rt && comp[i] == comp[i + 1]) continue;
        if (rt && __builtin_popcount(dirs[i + 1]) >= 2) continue;
        int save[64];
        memcpy(save, comp, sizeof save);
        if (rt) {
            const int a = comp[i + 1], b = comp[i];
            for (int x = 0; x < cells; x++)
                if (comp[x] == a) comp[x] = b;
        }
        if (dn) {
            if (comp[i] == comp[i + n]) {
                memcpy(comp, save, sizeof save);
                continue;
            }
            const int a = comp[i + n], b = comp[i];
            for (int x = 0; x < cells; x++)
                if (comp[x] == a) comp[x] = b;
        }
        const int di = dirs[i], dr = rt ? dirs[i + 1] : 0, dd = dn ? dirs[i + n] : 0;
        if (rt) {
            dirs[i] |= RT;
            dirs[i + 1] |= LT;
        }
        if (dn) {
            dirs[i] |= DN;
            dirs[i + n] |= UP;
        }
        enumerate(i + 1);
        dirs[i] = di;
        if (rt) dirs[i + 1] = dr;
        if (dn) dirs[i + n] = dd;
        memcpy(comp, save, sizeof save);
    }
}

int main(int argc, char **argv)
{
    if (argc < 2 || (n = atoi(argv[1])) < 1 || n > 6) {
        fprintf(stderr, "usage: verify n (1..6)\n");
        return 1;
    }
    cells = n * n;
    for (int i = 0; i < cells; i++) comp[i] = i;
    enumerate(0);
    long long sol = 0, puz = 0;
    for (int g = 0; g < 8; g++) {
        sol += fixsol[g];
        for (int m = 1; m < 16; m++) {
            if (fixpuz[g][m] % m) fprintf(stderr, "warning: g=%d m=%d not divisible\n", g, m);
            puz += fixpuz[g][m] / m;
        }
    }
    printf("n=%d  all full solutions: %lld  compact (unreduced): %lld  "
           "compact (reduced): %lld  compact starting positions (reduced): %lld\n",
           n, all_covers, compact_lab, sol / 8, puz / 8);
    if (sol % 8 || puz % 8) fprintf(stderr, "warning: Burnside sums not divisible by 8\n");
    return 0;
}
