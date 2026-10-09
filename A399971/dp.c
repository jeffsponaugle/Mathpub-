/*
 * dp.c - transfer-matrix count of compact Numberlink solutions (OEIS A399971),
 *        not reduced for symmetry.
 *
 * The grid is swept cell by cell in row-major order.  A state describes the
 * processed part as seen from the frontier (the last n processed cells):
 *   P : the cover itself - degree and path label of each frontier cell, and for
 *       each open path whether one of its ends is already a finalized terminal;
 *   Q : the set of frontier configurations of every partial linkage of the same
 *       endpoint pairs (vertex-disjoint paths, cells may stay empty), each with
 *       a flag telling whether it already left a processed cell empty.
 * Which terminals form a pair is only known once P's path is complete, so a Q
 * configuration records, per component, the endpoint it is tied to: the
 * finalized terminal of an open P path (by P label), a token for a pair whose
 * two terminals are known (shared by the two components that must still meet),
 * or none; and per open P path a requirement "must merge with path b" when a
 * Q path already joins the terminals of two open P paths.
 *
 * Dominance: a flag-0 configuration with a flag-1 twin is dropped (any
 * completion of it completes the twin, which leaves a cell empty).  When the
 * cover's own configuration has a flag-1 twin, the cover is not compact and
 * the state is dropped.  After a virtual extra row that only lets the last
 * row leave, the surviving states are exactly the compact covers.
 *
 * Build: cc -O3 -march=native -o dp dp.c
 * Usage: dp [-P] n      count (-P: all path covers, ignoring compactness)
 *        dp -D n        check every cover's sweep against exhaustive routing
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef uint64_t u64;

static int n, ponly;
static int forceP;         /* debug: follow one given cover */
static int fH[64], fV[64]; /* its edges right of / below each cell */

enum { ABSENT = 3, TOK = 32, COMPLETE = 64, FRESH_P = 8, FRESH_U = 9, FRESH_L = 10 };

typedef struct {
    int deg[8], lab[8]; /* deg ABSENT: outside the grid */
    int term[16];       /* per label: 1 if one end is a finalized terminal */
} PS;

typedef struct {
    int deg[8], lab[8];
    int att[16]; /* per Q label: 0, P label, TOK+t, COMPLETE */
    int req[16]; /* per P label: the open P path it must merge with */
    int flag;
} QC;

typedef struct {
    u64 lo, hi;
} Q128;

static inline int qless(Q128 a, Q128 b) { return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
static inline int qeq(Q128 a, Q128 b) { return a.hi == b.hi && a.lo == b.lo; }

/* ---- packing (n <= 7) ---- */
static u64 packP(const PS *p)
{
    u64 x = 0;
    for (int j = 0; j < n; j++) x |= (u64)(p->deg[j] | p->lab[j] << 2) << (5 * j);
    for (int L = 1; L <= n; L++) x |= (u64)(p->term[L] & 1) << (35 + L);
    return x;
}

static void unpackP(u64 x, PS *p)
{
    memset(p, 0, sizeof *p);
    for (int j = 0; j < n; j++) {
        const int v = x >> (5 * j) & 31;
        p->deg[j] = v & 3;
        p->lab[j] = v >> 2;
    }
    for (int L = 1; L <= n; L++) p->term[L] = x >> (35 + L) & 1;
}

/* lo: flag, cells (5 bits each), att (4 bits per label: 1..7 P, 8..15 token)
   hi: req (3 bits per P label) */
static Q128 packQ(const QC *q)
{
    Q128 k = {q->flag & 1, 0};
    for (int j = 0; j < n; j++) k.lo |= (u64)(q->deg[j] | q->lab[j] << 2) << (1 + 5 * j);
    for (int L = 1; L <= n; L++) {
        const int a = q->att[L];
        k.lo |= (u64)((a >= TOK ? 8 + (a - TOK) : a) & 15) << (36 + 4 * (L - 1));
        k.hi |= (u64)(q->req[L] & 7) << (3 * (L - 1));
    }
    return k;
}

static void unpackQ(Q128 k, QC *q)
{
    memset(q, 0, sizeof *q);
    q->flag = k.lo & 1;
    for (int j = 0; j < n; j++) {
        const int v = k.lo >> (1 + 5 * j) & 31;
        q->deg[j] = v & 3;
        q->lab[j] = v >> 2;
    }
    for (int L = 1; L <= n; L++) {
        const int v = k.lo >> (36 + 4 * (L - 1)) & 15;
        q->att[L] = v >= 8 ? TOK + (v - 8) : v;
        q->req[L] = k.hi >> (3 * (L - 1)) & 7;
    }
}

/* relabel P by first appearance; map: old label -> new (0 = gone) */
static void canonP(PS *p, int *map)
{
    int term2[16] = {0}, next = 1;
    for (int i = 0; i < 16; i++) map[i] = 0;
    for (int j = 0; j < n; j++) {
        const int L = p->lab[j];
        if (L && !map[L]) {
            map[L] = next++;
            term2[map[L]] = p->term[L];
        }
        p->lab[j] = map[L];
    }
    memcpy(p->term, term2, sizeof term2);
}

static long long bad_stuck, bad_close, bad_own, bad_pref;

static Q128 canonQ(QC *q, const int *pmap)
{
    int qmap[16] = {0}, att2[16] = {0}, req2[16] = {0}, tok[64], next = 1, ntok = 0;
    for (int i = 0; i < 64; i++) tok[i] = -1;
    for (int j = 0; j < n; j++) {
        const int L = q->lab[j];
        if (!L) continue;
        if (!qmap[L]) {
            qmap[L] = next++;
            int a = q->att[L];
            if (a >= 1 && a < TOK) {
                a = pmap[a];
                if (!a) bad_pref++;
            } else if (a >= TOK) {
                const int t = a - TOK;
                if (tok[t] < 0) tok[t] = ntok++;
                a = TOK + tok[t];
            }
            att2[qmap[L]] = a;
        }
        q->lab[j] = qmap[L];
    }
    for (int L = 1; L < 16; L++)
        if (q->req[L]) {
            if (!pmap[L] || !pmap[q->req[L]]) bad_pref++;
            req2[pmap[L]] = pmap[q->req[L]];
        }
    memcpy(q->att, att2, sizeof att2);
    memcpy(q->req, req2, sizeof req2);
    return packQ(q);
}

/* ---- transitions ---- */
typedef struct {
    int c, up, left, xreal;
    int merged, A, B; /* P labels of the up / left cell; B merged into A */
    int uterm, ulab;  /* the up cell leaves as a terminal of P path ulab */
    int closeL;       /* P label whose two terminals are now both known */
    int pmap[16];
} Ev;

/* join the endpoints tied to two merging Q components; -1 if impossible */
static int comb(QC *t, int x, int y)
{
    if (!x) return y;
    if (!y) return x;
    if (x == COMPLETE || y == COMPLETE) return -1;
    if (x >= TOK || y >= TOK) return x == y ? COMPLETE : -1;
    if (x == y) return COMPLETE;
    /* terminals of two different open P paths: fine only if they merge later */
    if (t->req[x] || t->req[y]) return -1;
    t->req[x] = y;
    t->req[y] = x;
    return COMPLETE;
}

static int qstep(const QC *q0, int qu, int ql, const Ev *ev, Q128 *out)
{
    QC t = *q0;
    const int c = ev->c;
    if (ev->merged) { /* P path B was merged into A */
        const int A = ev->A, B = ev->B;
        if (t.req[A] == B) t.req[A] = t.req[B] = 0; /* the requirement is met */
        for (int L = 1; L < 16; L++) {
            if (t.att[L] == B) t.att[L] = A;
            if (t.req[L] == B) t.req[L] = A;
        }
        if (t.req[B]) {
            if (t.req[A]) return 0;
            t.req[A] = t.req[B];
            t.req[B] = 0;
        }
    }
    if (qu && t.deg[c] >= 2) return 0;
    if (ql && t.deg[c - 1] >= 2) return 0;
    int a = ev->up ? t.lab[c] : 0, b = ev->left ? t.lab[c - 1] : 0;
    if (qu && t.deg[c] == 0) {
        a = FRESH_U;
        t.lab[c] = a;
        t.att[a] = 0;
    }
    if (ql && t.deg[c - 1] == 0) {
        b = FRESH_L;
        t.lab[c - 1] = b;
        t.att[b] = 0;
    }
    if (qu && ql && a == b) return 0;
    int X = 0;
    if (qu && ql) {
        const int m = comb(&t, t.att[a], t.att[b]);
        if (m < 0) return 0;
        for (int j = 0; j < n; j++)
            if (t.lab[j] == b) t.lab[j] = a;
        t.att[a] = m;
        t.att[b] = 0;
        X = a;
    } else if (qu) X = a;
    else if (ql) X = b;
    if (ql) t.deg[c - 1]++;
    if (ev->up) { /* the up cell leaves the frontier */
        const int fq = t.deg[c] + qu;
        if (ev->uterm) {
            if (fq != 1) return 0;
            const int K = qu ? X : t.lab[c];
            const int m = comb(&t, t.att[K], ev->ulab);
            if (m < 0) return 0;
            t.att[K] = m;
        } else {
            if (fq == 1) return 0;
            if (fq == 0) t.flag = 1;
        }
    }
    t.deg[c] = ev->xreal ? qu + ql : ABSENT;
    t.lab[c] = ev->xreal && qu + ql ? X : 0;
    if (ev->closeL) { /* both terminals of this pair are now known */
        if (t.req[ev->closeL]) return 0;
        for (int L = 1; L < 16; L++)
            if (t.att[L] == ev->closeL) t.att[L] = TOK + 15;
    }
    for (int K = 1; K < 16; K++) { /* finished paths leave the frontier */
        if (t.att[K] != COMPLETE) continue;
        for (int j = 0; j < n; j++)
            if (t.lab[j] == K) {
                if (t.deg[j] != 2) return 0;
                t.lab[j] = 0;
            }
        t.att[K] = 0;
    }
    { /* a component tied to an endpoint must still reach the frontier */
        int seen[16] = {0};
        for (int j = 0; j < n; j++) seen[t.lab[j]] = 1;
        for (int K = 1; K < 16; K++)
            if (t.att[K] && !seen[K]) {
                bad_stuck++;
                return 0;
            }
    }
    *out = canonQ(&t, ev->pmap);
    return 1;
}

/* ---- hash map of states ---- */
typedef struct {
    u64 p, off, count;
    int nq;
} Ent;

typedef struct {
    Ent *e;
    size_t cap, used;
    Q128 *arena;
    size_t acap, aused;
} Map;

static u64 hashkey(u64 p, const Q128 *q, int nq)
{
    u64 h = p * 0x9E3779B97F4A7C15ULL;
    for (int i = 0; i < nq; i++) {
        h ^= q[i].lo + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
        h *= 0xBF58476D1CE4E5B9ULL;
        h ^= q[i].hi + (h >> 29);
    }
    return h ^ (h >> 31);
}

static void map_init(Map *m, size_t cap)
{
    m->cap = cap;
    m->used = 0;
    m->e = malloc(cap * sizeof(Ent));
    for (size_t i = 0; i < cap; i++) m->e[i].nq = -1;
    m->acap = 1 << 16;
    m->aused = 0;
    m->arena = malloc(m->acap * sizeof(Q128));
}

static void map_free(Map *m)
{
    free(m->e);
    free(m->arena);
}

static void map_grow(Map *m)
{
    const size_t cap = m->cap * 2;
    Ent *e = malloc(cap * sizeof(Ent));
    for (size_t i = 0; i < cap; i++) e[i].nq = -1;
    for (size_t i = 0; i < m->cap; i++) {
        if (m->e[i].nq < 0) continue;
        u64 h = hashkey(m->e[i].p, m->arena + m->e[i].off, m->e[i].nq) & (cap - 1);
        while (e[h].nq >= 0) h = (h + 1) & (cap - 1);
        e[h] = m->e[i];
    }
    free(m->e);
    m->e = e;
    m->cap = cap;
}

static void map_add(Map *m, u64 p, const Q128 *q, int nq, u64 cnt)
{
    if (2 * (m->used + 1) > m->cap) map_grow(m);
    u64 h = hashkey(p, q, nq) & (m->cap - 1);
    for (;;) {
        Ent *e = &m->e[h];
        if (e->nq < 0) {
            if (m->aused + nq > m->acap) {
                while (m->aused + nq > m->acap) m->acap *= 2;
                m->arena = realloc(m->arena, m->acap * sizeof(Q128));
                if (!m->arena) {
                    perror("realloc");
                    exit(1);
                }
            }
            memcpy(m->arena + m->aused, q, nq * sizeof(Q128));
            e->p = p;
            e->off = m->aused;
            e->nq = nq;
            e->count = cnt;
            m->aused += nq;
            m->used++;
            return;
        }
        if (e->p == p && e->nq == nq && !memcmp(m->arena + e->off, q, nq * sizeof(Q128))) {
            e->count += cnt;
            return;
        }
        h = (h + 1) & (m->cap - 1);
    }
}

static int cmpq(const void *a, const void *b)
{
    const Q128 x = *(const Q128 *)a, y = *(const Q128 *)b;
    return qless(x, y) ? -1 : qless(y, x);
}

static Q128 *buf;
static size_t bufcap;

/* Dominance.  Cells with label 0 are either empty (degree 0) or inside a
   finished Q path (degree 2, "done").  If two configurations differ only in
   such cells, and wherever they differ the first is empty and the second done,
   then every completion of the second also completes the first, which then
   leaves those cells empty.  So the second can be dropped; the same holds for
   equal configurations with flags 1 and 0.  (core: lo without the flag and
   with done cells written as empty.) */
static u64 CELLMASK[8]; /* degree bits of each frontier cell in lo */

static inline u64 qcore(Q128 q, u64 *emask)
{
    u64 lo = q.lo & ~1ULL, e = 0;
    for (int j = 0; j < n; j++) {
        const int v = lo >> (1 + 5 * j) & 31;
        if (v == 0) e |= 1u << j; /* empty */
        else if (v == 2) lo &= ~CELLMASK[j]; /* done: write as empty */
    }
    *emask = e;
    return lo;
}

typedef struct {
    u64 hi, core, e;
    Q128 q;
} DomKey;

static int cmpdom(const void *a, const void *b)
{
    const DomKey *x = a, *y = b;
    if (x->hi != y->hi) return x->hi < y->hi ? -1 : 1;
    if (x->core != y->core) return x->core < y->core ? -1 : 1;
    return 0;
}

static DomKey *dk;
static size_t dkcap;

/* reduce buf[0..nb) to its undominated configurations; returns the new size
   and sets *own_dom if the configuration own is dominated */
static size_t dominance(size_t nb, Q128 own, int *own_dom)
{
    if (dkcap < nb) {
        dkcap = 2 * nb + 16;
        dk = realloc(dk, dkcap * sizeof(DomKey));
    }
    for (size_t i = 0; i < nb; i++) {
        dk[i].hi = buf[i].hi;
        dk[i].core = qcore(buf[i], &dk[i].e);
        dk[i].q = buf[i];
    }
    qsort(dk, nb, sizeof(DomKey), cmpdom);
    u64 ownE;
    const u64 ownc = qcore(own, &ownE);
    *own_dom = 0;
    size_t m = 0;
    for (size_t g = 0; g < nb;) {
        size_t h = g;
        while (h < nb && dk[h].hi == dk[g].hi && dk[h].core == dk[g].core) h++;
        const int owngroup = dk[g].hi == own.hi && dk[g].core == ownc;
        for (size_t i = g; i < h; i++) {
            int dominated = 0;
            for (size_t j = g; j < h && !dominated; j++) {
                if (j == i) continue;
                const u64 ei = dk[i].e, ej = dk[j].e;
                if ((ei & ej) != ei) continue; /* j must be empty wherever i is */
                if (ej != ei) dominated = 1;
                else if ((dk[j].q.lo & 1) > (dk[i].q.lo & 1)) dominated = 1;
                else if (dk[j].q.lo == dk[i].q.lo && j < i) dominated = 1; /* duplicate */
            }
            if (!dominated) buf[m++] = dk[i].q;
        }
        if (owngroup)
            for (size_t j = g; j < h; j++) {
                const u64 ej = dk[j].e;
                if ((ownE & ej) == ownE && (ej != ownE || (dk[j].q.lo & 1))) *own_dom = 1;
            }
        g = h;
    }
    return m;
}
static QC *qcs;
static size_t qcscap;

static void printQ(Q128 k);

/* process one state at cell (r, c) */
static void step(Map *next, u64 pkey, const Q128 *qs, int nq, u64 count, int r, int c)
{
    PS p;
    unpackP(pkey, &p);
    const int xreal = r < n;
    const int up = p.deg[c] != ABSENT;
    const int left = xreal && c > 0;
    if (qcscap < (size_t)nq) {
        qcscap = 2 * nq;
        qcs = realloc(qcs, qcscap * sizeof(QC));
    }
    for (int i = 0; i < nq; i++) unpackQ(qs[i], &qcs[i]);

    for (int eu = 0; eu <= (up && xreal); eu++)
        for (int el = 0; el <= left; el++) {
            if (forceP && xreal) {
                if (up && eu != fV[(r - 1) * n + c]) continue;
                if (left && el != fH[r * n + c - 1]) continue;
            }
            PS t = p;
            if (eu && t.deg[c] >= 2) continue;
            if (el && t.deg[c - 1] >= 2) continue;
            const int A = up ? t.lab[c] : 0, B = left ? t.lab[c - 1] : 0;
            if (eu && el && A == B) continue;
            int uterm = 0;
            if (up) {
                const int fd = t.deg[c] + eu;
                if (fd == 0) continue;
                uterm = fd == 1;
            }
            Ev ev;
            memset(&ev, 0, sizeof ev);
            ev.c = c;
            ev.up = up;
            ev.left = left;
            ev.xreal = xreal;
            int X = 0;
            if (eu && el) {
                for (int j = 0; j < n; j++)
                    if (t.lab[j] == B) t.lab[j] = A;
                t.term[A] += t.term[B];
                t.term[B] = 0;
                ev.merged = 1;
                ev.A = A;
                ev.B = B;
                X = A;
            } else if (eu) X = A;
            else if (el) X = B;
            else if (xreal) {
                X = FRESH_P;
                t.term[X] = 0;
            }
            if (el) t.deg[c - 1]++;
            if (uterm) {
                t.term[A]++;
                ev.uterm = 1;
                ev.ulab = A;
            }
            t.deg[c] = xreal ? eu + el : ABSENT;
            t.lab[c] = xreal ? X : 0;
            for (int L = 1; L < 16; L++)
                if (t.term[L] >= 2) {
                    ev.closeL = L;
                    t.term[L] = 0;
                    for (int j = 0; j < n; j++)
                        if (t.lab[j] == L) {
                            if (t.deg[j] != 2) bad_close++;
                            t.lab[j] = 0;
                        }
                }
            canonP(&t, ev.pmap);
            const u64 pk = packP(&t);

            if (ponly) {
                Q128 z = {0, 0};
                map_add(next, pk, &z, 1, count);
                continue;
            }

            QC own; /* the cover's own configuration */
            memset(&own, 0, sizeof own);
            for (int j = 0; j < n; j++) {
                own.deg[j] = t.deg[j];
                own.lab[j] = t.deg[j] == 0 ? 0 : t.lab[j]; /* Q gives isolated cells no label */
            }
            for (int L = 1; L <= n; L++) own.att[L] = t.term[L] ? L : 0;
            int idmap[16];
            for (int i = 0; i < 16; i++) idmap[i] = i;
            const Q128 ownk = canonQ(&own, idmap);

            if (bufcap < (size_t)nq * 4 + 1) {
                bufcap = (size_t)nq * 8 + 16;
                buf = realloc(buf, bufcap * sizeof(Q128));
            }
            size_t nb = 0;
            for (int i = 0; i < nq; i++)
                for (int qu = 0; qu <= (up && xreal); qu++)
                    for (int ql = 0; ql <= left; ql++) {
                        Q128 o;
                        if (qstep(&qcs[i], qu, ql, &ev, &o)) buf[nb++] = o;
                    }
            int dominated = 0, present = 0;
            size_t m = dominance(nb, ownk, &dominated);
            if (dominated) continue;
            qsort(buf, m, sizeof(Q128), cmpq);
            for (size_t i = 0; i < m; i++)
                if (qeq(buf[i], ownk)) present = 1;
            if (!present) {
                if (!bad_own++ && forceP) {
                    printf("own config missing at cell (%d,%d); own =", r, c);
                    printQ(ownk);
                    printf("  set:\n");
                    for (size_t i = 0; i < m; i++) printQ(buf[i]);
                }
            }
            map_add(next, pk, buf, (int)m, count);
        }
}

static void init_map(Map *m)
{
    for (int j = 0; j < n; j++) CELLMASK[j] = 3ULL << (1 + 5 * j);
    map_init(m, 16);
    PS p0;
    memset(&p0, 0, sizeof p0);
    for (int j = 0; j < n; j++) p0.deg[j] = ABSENT;
    QC q0;
    memset(&q0, 0, sizeof q0);
    for (int j = 0; j < n; j++) q0.deg[j] = ABSENT;
    const Q128 qk = packQ(&q0);
    map_add(m, packP(&p0), &qk, 1, 1);
}

/* ---- debugging: single covers against exhaustive routing ---- */
static void printQ(Q128 k)
{
    QC q;
    unpackQ(k, &q);
    printf("    flag %d |", q.flag);
    for (int j = 0; j < n; j++) {
        if (q.deg[j] == ABSENT) printf(" --");
        else printf(" %d%c", q.deg[j], q.lab[j] ? 'a' + q.lab[j] - 1 : '.');
    }
    printf(" | att");
    for (int L = 1; L <= n; L++)
        if (q.att[L])
            printf(" %c:%s%d", 'a' + L - 1, q.att[L] >= TOK ? "T" : "P",
                   q.att[L] >= TOK ? q.att[L] - TOK : q.att[L]);
    printf(" | req");
    for (int L = 1; L <= n; L++)
        if (q.req[L]) printf(" P%d-P%d", L, q.req[L]);
    printf("\n");
}

static int sweep(int verbose)
{
    Map cur, nxt;
    init_map(&cur);
    for (int r = 0; r <= n; r++)
        for (int c = 0; c < n; c++) {
            map_init(&nxt, 16);
            for (size_t i = 0; i < cur.cap; i++) {
                const Ent *e = &cur.e[i];
                if (e->nq >= 0) step(&nxt, e->p, cur.arena + e->off, e->nq, e->count, r, c);
            }
            map_free(&cur);
            cur = nxt;
            if (!verbose) continue;
            printf("  after (%d,%d): %zu state(s)\n", r, c, cur.used);
            for (size_t i = 0; i < cur.cap; i++) {
                const Ent *e = &cur.e[i];
                if (e->nq < 0) continue;
                PS p;
                unpackP(e->p, &p);
                printf("   P:");
                for (int j = 0; j < n; j++) {
                    if (p.deg[j] == ABSENT) printf(" --");
                    else printf(" %d%c", p.deg[j], p.lab[j] ? 'a' + p.lab[j] - 1 : '.');
                }
                printf("  term");
                for (int L = 1; L <= n; L++)
                    if (p.term[L]) printf(" P%d", L);
                printf("\n");
                for (int k = 0; k < e->nq; k++) printQ(cur.arena[e->off + k]);
            }
        }
    const int alive = cur.used > 0;
    map_free(&cur);
    return alive;
}

static int nvp, nvs[32], nvt[32];
static unsigned char nvused[64];
static int nv_route(int j);

static int nv_extend(int j, int x)
{
    const int r = x / n, c = x % n;
    const int nb[4] = {r ? x - n : -1, c + 1 < n ? x + 1 : -1, r + 1 < n ? x + n : -1, c ? x - 1 : -1};
    for (int d = 0; d < 4; d++) {
        const int y = nb[d];
        if (y < 0) continue;
        if (y == nvt[j]) {
            if (nv_route(j + 1)) return 1;
        } else if (!nvused[y]) {
            nvused[y] = 1;
            const int f = nv_extend(j, y);
            nvused[y] = 0;
            if (f) return 1;
        }
    }
    return 0;
}

static int nv_route(int j)
{
    if (j == nvp) {
        for (int i = 0; i < n * n; i++)
            if (!nvused[i]) return 1;
        return 0;
    }
    return nv_extend(j, nvs[j]);
}

static void debug_all(void)
{
    const int NN = n * n;
    int ne = 0, eu[64], ev[64];
    for (int i = 0; i < NN; i++) {
        if (i % n + 1 < n) eu[ne] = i, ev[ne] = i + 1, ne++;
        if (i / n + 1 < n) eu[ne] = i, ev[ne] = i + n, ne++;
    }
    long long covers = 0, compact = 0, mism = 0;
    forceP = 1;
    for (long long mask = 0; mask < (1LL << ne); mask++) {
        int deg[64] = {0}, comp[64], ok = 1;
        for (int i = 0; i < NN; i++) comp[i] = i;
        memset(fH, 0, sizeof fH);
        memset(fV, 0, sizeof fV);
        for (int e = 0; e < ne && ok; e++)
            if (mask >> e & 1) {
                const int a = eu[e], b = ev[e];
                if (++deg[a] > 2 || ++deg[b] > 2) ok = 0;
                const int ca = comp[a], cb = comp[b];
                if (ca == cb) ok = 0;
                for (int i = 0; i < NN; i++)
                    if (comp[i] == cb) comp[i] = ca;
                if (b == a + 1) fH[a] = 1;
                else fV[a] = 1;
            }
        for (int i = 0; i < NN && ok; i++)
            if (!deg[i]) ok = 0;
        if (!ok) continue;
        covers++;
        nvp = 0;
        memset(nvused, 0, sizeof nvused);
        int seen[64] = {0};
        for (int i = 0; i < NN; i++) {
            if (deg[i] != 1 || seen[i]) continue;
            for (int j = 0; j < NN; j++)
                if (j != i && deg[j] == 1 && comp[j] == comp[i]) {
                    nvs[nvp] = i;
                    nvt[nvp] = j;
                    nvp++;
                    seen[i] = seen[j] = 1;
                    nvused[i] = nvused[j] = 1;
                }
        }
        const int truth = !nv_route(0), got = sweep(0);
        compact += truth;
        if (truth != got && ++mism == 1) {
            printf("MISMATCH (truth %s, dp %s):\n", truth ? "compact" : "not compact",
                   got ? "compact" : "not compact");
            for (int r = 0; r < n; r++) {
                for (int c = 0; c < n; c++) printf("o%s", c + 1 < n ? (fH[r * n + c] ? "-" : " ") : "");
                printf("\n");
                if (r + 1 < n) {
                    for (int c = 0; c < n; c++) printf("%s ", fV[r * n + c] ? "|" : " ");
                    printf("\n");
                }
            }
            sweep(1);
        }
    }
    printf("n=%d: %lld covers, %lld compact, %lld mismatches\n", n, covers, compact, mism);
    if (bad_stuck || bad_close || bad_own || bad_pref)
        printf("  invariant violations: stuck %lld close %lld own %lld pref %lld\n", bad_stuck,
               bad_close, bad_own, bad_pref);
}

int main(int argc, char **argv)
{
    int ai = 1, debug = 0;
    for (; ai < argc && argv[ai][0] == '-'; ai++) {
        if (!strcmp(argv[ai], "-P")) ponly = 1;
        else if (!strcmp(argv[ai], "-D")) debug = 1;
    }
    if (ai >= argc || (n = atoi(argv[ai])) < 1 || n > 7) {
        fprintf(stderr, "usage: dp [-P | -D] n   (1 <= n <= 7)\n");
        return 1;
    }
    if (debug) {
        debug_all();
        return 0;
    }
    Map cur, nxt;
    init_map(&cur);
    const clock_t t0 = clock();
    for (int r = 0; r <= n; r++)
        for (int c = 0; c < n; c++) {
            map_init(&nxt, 1 << 10);
            size_t maxq = 0, sumq = 0;
            for (size_t i = 0; i < cur.cap; i++) {
                const Ent *e = &cur.e[i];
                if (e->nq < 0) continue;
                if ((size_t)e->nq > maxq) maxq = e->nq;
                sumq += e->nq;
                step(&nxt, e->p, cur.arena + e->off, e->nq, e->count, r, c);
            }
            fprintf(stderr, "cell (%d,%d): %zu states, avg Q %.1f, max Q %zu -> %zu  [%.1fs]\n", r,
                    c, cur.used, cur.used ? (double)sumq / cur.used : 0, maxq, nxt.used,
                    (double)(clock() - t0) / CLOCKS_PER_SEC);
            map_free(&cur);
            cur = nxt;
        }
    u64 total = 0;
    for (size_t i = 0; i < cur.cap; i++)
        if (cur.e[i].nq >= 0) total += cur.e[i].count;
    printf("n=%d: %s %llu\n", n, ponly ? "path covers" : "compact solutions (not reduced)",
           (unsigned long long)total);
    if (bad_stuck || bad_close || bad_own || bad_pref)
        printf("  invariant violations: stuck %lld close %lld own %lld pref %lld\n", bad_stuck,
               bad_close, bad_own, bad_pref);
    return 0;
}
