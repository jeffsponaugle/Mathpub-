/* biggraphs.c -- generate large test graphs with known planarity, in
 * sparse6 on stdout (loops and parallel edges allowed).
 *
 * usage: biggraphs FAMILY args...          (seed: optional last argument
 *                                            for the random families)
 *   planar:
 *     grid R C          R x C grid                     path N     path
 *     cylinder R C      path_R x cycle_C               cycle N    cycle
 *     ladder N          2 x N grid                     star N     K_{1,N}
 *     wheel N           hub + cycle_N                  k2n N      K_{2,N}
 *     tree N [seed]     random tree
 *     apollonian N [seed]        random Apollonian network (maximal planar, 3N-6 edges)
 *     apollodel N D [seed]       apollonian with D random edges deleted
 *     multi N [seed]             apollonian with parallel edges and loops added (planar)
 *     forest N [seed]            random forest with several components
 *   non-planar:
 *     torus R C         cycle_R x cycle_C (R,C >= 3)   mobius N   Moebius ladder, 2N vertices (N >= 3)
 *     k3n N             K_{3,N} (N >= 3)               k5sub N    K5 with every edge subdivided into paths (N vertices total)
 *     apollok5 N [seed]          apollonian + K5 on 5 random vertices
 *     apollok33 N [seed]         apollonian + K_{3,3} on 6 random vertices
 *     k5far N                    path on N vertices plus a separate K5 component
 *     multink33 N [seed]         k3n-like multigraph: K_{3,3} plus a long path with parallel edges and loops
 */

#include "gtools.h"
#include "nausparse.h"

static size_t cap, ne;
static int *eu, *ev;
static int nv;

static void
edge(int u, int v)
{
    if (ne == cap)
    {
        cap = cap ? 2*cap : 1024;
        eu = (int*)realloc(eu, cap*sizeof(int));
        ev = (int*)realloc(ev, cap*sizeof(int));
        if (!eu || !ev) { fprintf(stderr, "out of memory\n"); exit(1); }
    }
    eu[ne] = u; ev[ne] = v; ++ne;
}

static unsigned long long rng_state = 88172645463325252ULL;
static unsigned long long
rnd(void)   /* xorshift64 */
{
    rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17;
    return rng_state;
}
static int rndint(int n) { return (int)(rnd() % (unsigned long long)n); }

/* Build the sparsegraph from the edge list and write it as sparse6. */
static void
emit(void)
{
    SG_DECL(sg);
    size_t i, nde = 0;
    int *pos;

    for (i = 0; i < ne; ++i) nde += (eu[i] == ev[i]) ? 1 : 2;
    SG_ALLOC(sg, nv, nde, "biggraphs");
    sg.nv = nv; sg.nde = nde;
    for (i = 0; i < (size_t)nv; ++i) sg.d[i] = 0;
    for (i = 0; i < ne; ++i) { ++sg.d[eu[i]]; if (eu[i] != ev[i]) ++sg.d[ev[i]]; }
    sg.v[0] = 0;
    for (i = 1; i < (size_t)nv; ++i) sg.v[i] = sg.v[i-1] + sg.d[i-1];
    pos = (int*)calloc(nv, sizeof(int));
    for (i = 0; i < ne; ++i)
    {
        sg.e[sg.v[eu[i]] + pos[eu[i]]++] = ev[i];
        if (eu[i] != ev[i]) sg.e[sg.v[ev[i]] + pos[ev[i]]++] = eu[i];
    }
    writes6_sg(stdout, &sg);
    fprintf(stderr, "biggraphs: n=%d edges=%lu\n", nv, (unsigned long)ne);
}

static void
grid(int R, int C, int wrapR, int wrapC)
{
    int r, c;
    nv = R*C;
    for (r = 0; r < R; ++r)
        for (c = 0; c < C; ++c)
        {
            if (c + 1 < C) edge(r*C + c, r*C + c + 1);
            else if (wrapC && C > 2) edge(r*C + c, r*C);
            if (r + 1 < R) edge(r*C + c, (r+1)*C + c);
            else if (wrapR && R > 2) edge(r*C + c, c);
        }
}

/* Random Apollonian network: start with a triangle, repeatedly split a
 * random face by a new vertex.  Faces kept in an array of triples. */
static void
apollonian(int N)
{
    int *fa, *fb, *fc, nf, f, a, b, c, v;

    if (N < 3) N = 3;
    nv = N;
    fa = (int*)malloc((2*N+1)*sizeof(int)); fb = (int*)malloc((2*N+1)*sizeof(int));
    fc = (int*)malloc((2*N+1)*sizeof(int));
    edge(0,1); edge(1,2); edge(0,2);
    fa[0]=0; fb[0]=1; fc[0]=2;  fa[1]=0; fb[1]=1; fc[1]=2;   /* two sides of the triangle */
    nf = 2;
    for (v = 3; v < N; ++v)
    {
        f = rndint(nf);
        a = fa[f]; b = fb[f]; c = fc[f];
        edge(v,a); edge(v,b); edge(v,c);
        fa[f]=a; fb[f]=b; fc[f]=v;
        fa[nf]=a; fb[nf]=c; fc[nf]=v; ++nf;
        fa[nf]=b; fb[nf]=c; fc[nf]=v; ++nf;
    }
    free(fa); free(fb); free(fc);
}

static void
delete_random_edges(int D)
{
    int i, j;
    for (i = 0; i < D && ne > 0; ++i)
    {
        j = rndint((int)ne);
        eu[j] = eu[ne-1]; ev[j] = ev[ne-1]; --ne;
    }
}

static void
random_tree(int N, int first)   /* vertices first..first+N-1 */
{
    int v;
    for (v = 1; v < N; ++v) edge(first + v, first + rndint(v));
}

static int
pick_distinct(int *out, int k)  /* k distinct random vertices */
{
    int i, j, x, ok;
    for (i = 0; i < k; ++i)
    {
        do { x = rndint(nv); ok = 1; for (j = 0; j < i; ++j) if (out[j] == x) ok = 0; } while (!ok);
        out[i] = x;
    }
    return k;
}

int
main(int argc, char *argv[])
{
    const char *fam;
    int N, R, C, D, i, j, x[6];

    if (argc < 3) { fprintf(stderr, "usage: biggraphs FAMILY N [args] [seed]  (see source)\n"); return 2; }
    fam = argv[1];
    N = atoi(argv[2]);
    for (i = 1; i < argc; ++i) ;                     /* seed = last arg if there is one more than needed */

#define SEEDFROM(k) if (argc > (k)) rng_state ^= (unsigned long long)atoll(argv[k]) * 0x9E3779B97F4A7C15ULL

    if (strcmp(fam, "grid") == 0)          { R = N; C = atoi(argv[3]); grid(R, C, 0, 0); }
    else if (strcmp(fam, "cylinder") == 0) { R = N; C = atoi(argv[3]); grid(R, C, 0, 1); }
    else if (strcmp(fam, "torus") == 0)    { R = N; C = atoi(argv[3]); grid(R, C, 1, 1); }
    else if (strcmp(fam, "ladder") == 0)   { grid(2, N, 0, 0); }
    else if (strcmp(fam, "mobius") == 0)   /* 2N vertices: cycle 0..2N-1 plus rungs i -- i+N */
    {
        nv = 2*N;
        for (i = 0; i < 2*N; ++i) edge(i, (i+1) % (2*N));
        for (i = 0; i < N; ++i) edge(i, i + N);
    }
    else if (strcmp(fam, "path") == 0)     { nv = N; for (i = 0; i + 1 < N; ++i) edge(i, i+1); }
    else if (strcmp(fam, "cycle") == 0)    { nv = N; for (i = 0; i < N; ++i) edge(i, (i+1) % N); }
    else if (strcmp(fam, "star") == 0)     { nv = N + 1; for (i = 1; i <= N; ++i) edge(0, i); }
    else if (strcmp(fam, "wheel") == 0)    { nv = N + 1; for (i = 1; i <= N; ++i) { edge(0, i); edge(i, i % N + 1); } }
    else if (strcmp(fam, "k2n") == 0)      { nv = N + 2; for (i = 2; i < N + 2; ++i) { edge(0, i); edge(1, i); } }
    else if (strcmp(fam, "k3n") == 0)      { nv = N + 3; for (i = 3; i < N + 3; ++i) { edge(0, i); edge(1, i); edge(2, i); } }
    else if (strcmp(fam, "tree") == 0)     { SEEDFROM(3); nv = N; random_tree(N, 0); }
    else if (strcmp(fam, "forest") == 0)   /* ~N/1000 random trees */
    {
        int first = 0, size;
        SEEDFROM(3); nv = N;
        while (first < N)
        {
            size = 1 + rndint(2000); if (first + size > N) size = N - first;
            random_tree(size, first); first += size;
        }
    }
    else if (strcmp(fam, "apollonian") == 0) { SEEDFROM(3); apollonian(N); }
    else if (strcmp(fam, "apollodel") == 0)  { D = atoi(argv[3]); SEEDFROM(4); apollonian(N); delete_random_edges(D); }
    else if (strcmp(fam, "multi") == 0)      /* planar multigraph: duplicate ~N/10 edges, add ~N/10 loops */
    {
        SEEDFROM(3); apollonian(N);
        for (i = 0; i < N/10; ++i) { j = rndint((int)ne); edge(eu[j], ev[j]); }
        for (i = 0; i < N/10; ++i) { j = rndint(N); edge(j, j); }
    }
    else if (strcmp(fam, "apollok5") == 0)
    {
        SEEDFROM(3); apollonian(N); pick_distinct(x, 5);
        for (i = 0; i < 5; ++i) for (j = i+1; j < 5; ++j) edge(x[i], x[j]);
    }
    else if (strcmp(fam, "apollok33") == 0)
    {
        SEEDFROM(3); apollonian(N); pick_distinct(x, 6);
        for (i = 0; i < 3; ++i) for (j = 3; j < 6; ++j) edge(x[i], x[j]);
    }
    else if (strcmp(fam, "k5sub") == 0)    /* K5 on vertices 0..4, each of the 10 edges a path with (N-5)/10 inner vertices */
    {
        int inner = (N - 5) / 10, next = 5, prev, k;
        nv = 5 + 10*inner;
        for (i = 0; i < 5; ++i) for (j = i+1; j < 5; ++j)
        {
            prev = i;
            for (k = 0; k < inner; ++k) { edge(prev, next); prev = next++; }
            edge(prev, j);
        }
    }
    else if (strcmp(fam, "k5far") == 0)    /* disconnected: path on N vertices, K5 on 5 more */
    {
        nv = N + 5;
        for (i = 0; i + 1 < N; ++i) edge(i, i+1);
        for (i = N; i < N+5; ++i) for (j = i+1; j < N+5; ++j) edge(i, j);
    }
    else if (strcmp(fam, "multink33") == 0) /* non-planar multigraph: K3,3 on 0..5, path 5..N-1 with duplicates and loops */
    {
        SEEDFROM(3); nv = N;
        for (i = 0; i < 3; ++i) for (j = 3; j < 6; ++j) edge(i, j);
        for (i = 5; i + 1 < N; ++i) { edge(i, i+1); if (rndint(4) == 0) edge(i, i+1); if (rndint(8) == 0) edge(i, i); }
    }
    else { fprintf(stderr, "biggraphs: unknown family %s\n", fam); return 2; }

    emit();
    return 0;
}
