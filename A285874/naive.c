/* naive.c - deliberately simple, independent perft for cross-checking perft.cpp.
 *
 * Shares no code or representation with perft.cpp: 0x88 mailbox board, pseudo-legal
 * move generation, legality tested by making each move and scanning for attacks on the
 * king. No bulk counting, no hashing. Slow but easy to audit.
 *
 * Usage: naive "<FEN>" depth [split]   (split: print per-move counts at ply 1)
 * Build: clang -O3 -o naive naive.c -lpthread
 */
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { EMPTY = 0, P = 1, N, B, R, Q, K };
#define WHITE 0x08
#define BLACK 0x10
#define COLOR(x) ((x) & 0x18)
#define TYPE(x) ((x) & 0x07)
#define ONBOARD(s) (!((s) & 0x88))

typedef struct {
  uint8_t sq[128];
  int stm;       /* WHITE or BLACK */
  int castle;    /* 1 WK, 2 WQ, 4 BK, 8 BQ */
  int ep;        /* 0x88 square or -1 */
  int king[2];   /* [0] white, [1] black */
} Board;

typedef struct { int from, to, promo, flags; } Mv; /* flags: 1 ep, 2 castle, 4 double */

static const int KN[8] = {33, 31, 18, 14, -33, -31, -18, -14};
static const int KG[8] = {1, -1, 16, -16, 17, 15, -17, -15};
static const int DG[4] = {17, 15, -17, -15};
static const int OR[4] = {1, -1, 16, -16};

static int attacked(const Board* b, int s, int by) {
  int i, t;
  int pd = by == WHITE ? -16 : 16; /* attacking pawn sits one rank "behind" s from its own view */
  for (i = -1; i <= 1; i += 2) {
    t = s + pd + i;
    if (ONBOARD(t) && b->sq[t] == (by | P)) return 1;
  }
  for (i = 0; i < 8; i++) {
    t = s + KN[i];
    if (ONBOARD(t) && b->sq[t] == (by | N)) return 1;
    t = s + KG[i];
    if (ONBOARD(t) && b->sq[t] == (by | K)) return 1;
  }
  for (i = 0; i < 4; i++) {
    for (t = s + DG[i]; ONBOARD(t); t += DG[i]) {
      if (b->sq[t]) {
        if (b->sq[t] == (by | B) || b->sq[t] == (by | Q)) return 1;
        break;
      }
    }
    for (t = s + OR[i]; ONBOARD(t); t += OR[i]) {
      if (b->sq[t]) {
        if (b->sq[t] == (by | R) || b->sq[t] == (by | Q)) return 1;
        break;
      }
    }
  }
  return 0;
}

static int addMove(Mv* list, int n, int from, int to, int promo, int flags) {
  list[n].from = from; list[n].to = to; list[n].promo = promo; list[n].flags = flags;
  return n + 1;
}

static int genPseudo(const Board* b, Mv* list) {
  int n = 0, s, i, t, us = b->stm, them = us ^ 0x18;
  int up = us == WHITE ? 16 : -16;
  int startRank = us == WHITE ? 1 : 6, lastRank = us == WHITE ? 7 : 0;
  for (s = 0; s < 128; s++) {
    if (!ONBOARD(s) || COLOR(b->sq[s]) != us) continue;
    int pc = TYPE(b->sq[s]);
    if (pc == P) {
      t = s + up;
      if (ONBOARD(t) && !b->sq[t]) {
        if ((t >> 4) == lastRank) {
          for (i = Q; i >= N; i--) n = addMove(list, n, s, t, i, 0);
        } else {
          n = addMove(list, n, s, t, 0, 0);
          if ((s >> 4) == startRank && !b->sq[t + up]) n = addMove(list, n, s, t + up, 0, 4);
        }
      }
      for (i = -1; i <= 1; i += 2) {
        t = s + up + i;
        if (!ONBOARD(t)) continue;
        if (COLOR(b->sq[t]) == them) {
          if ((t >> 4) == lastRank) {
            int pr;
            for (pr = Q; pr >= N; pr--) n = addMove(list, n, s, t, pr, 0);
          } else {
            n = addMove(list, n, s, t, 0, 0);
          }
        } else if (t == b->ep) {
          n = addMove(list, n, s, t, 0, 1);
        }
      }
    } else if (pc == N || pc == K) {
      const int* d = pc == N ? KN : KG;
      for (i = 0; i < 8; i++) {
        t = s + d[i];
        if (ONBOARD(t) && COLOR(b->sq[t]) != us) n = addMove(list, n, s, t, 0, 0);
      }
    } else {
      int k;
      if (pc == B || pc == Q)
        for (k = 0; k < 4; k++)
          for (t = s + DG[k]; ONBOARD(t); t += DG[k]) {
            if (COLOR(b->sq[t]) == us) break;
            n = addMove(list, n, s, t, 0, 0);
            if (b->sq[t]) break;
          }
      if (pc == R || pc == Q)
        for (k = 0; k < 4; k++)
          for (t = s + OR[k]; ONBOARD(t); t += OR[k]) {
            if (COLOR(b->sq[t]) == us) break;
            n = addMove(list, n, s, t, 0, 0);
            if (b->sq[t]) break;
          }
    }
  }
  /* castling: king and rook on home squares, path empty, king not in check and not crossing attack */
  {
    int home = us == WHITE ? 0x04 : 0x74, ks = us == WHITE ? 1 : 4, qs = us == WHITE ? 2 : 8;
    if (b->sq[home] == (us | K) && !attacked(b, home, them)) {
      if ((b->castle & ks) && b->sq[home + 3] == (us | R) && !b->sq[home + 1] && !b->sq[home + 2] &&
          !attacked(b, home + 1, them) && !attacked(b, home + 2, them))
        n = addMove(list, n, home, home + 2, 0, 2);
      if ((b->castle & qs) && b->sq[home - 4] == (us | R) && !b->sq[home - 1] && !b->sq[home - 2] &&
          !b->sq[home - 3] && !attacked(b, home - 1, them) && !attacked(b, home - 2, them))
        n = addMove(list, n, home, home - 2, 0, 2);
    }
  }
  return n;
}

static void clearRights(Board* b, int s) {
  if (s == 0x04) b->castle &= ~3;
  if (s == 0x74) b->castle &= ~12;
  if (s == 0x07) b->castle &= ~1;
  if (s == 0x00) b->castle &= ~2;
  if (s == 0x77) b->castle &= ~4;
  if (s == 0x70) b->castle &= ~8;
}

static void make(Board* b, const Mv* m) {
  int us = b->stm, pc = b->sq[m->from];
  b->sq[m->to] = m->promo ? (us | m->promo) : pc;
  b->sq[m->from] = EMPTY;
  if (m->flags & 1) b->sq[m->to + (us == WHITE ? -16 : 16)] = EMPTY;
  if (m->flags & 2) {
    if (m->to > m->from) { b->sq[m->to - 1] = b->sq[m->to + 1]; b->sq[m->to + 1] = EMPTY; }
    else { b->sq[m->to + 1] = b->sq[m->to - 2]; b->sq[m->to - 2] = EMPTY; }
  }
  if (TYPE(pc) == K) b->king[us == WHITE ? 0 : 1] = m->to;
  clearRights(b, m->from);
  clearRights(b, m->to);
  b->ep = (m->flags & 4) ? (m->from + m->to) / 2 : -1;
  b->stm = us ^ 0x18;
}

static uint64_t perft(const Board* b, int depth) {
  Mv list[256];
  int n = genPseudo(b, list), i;
  uint64_t total = 0;
  for (i = 0; i < n; i++) {
    Board c = *b;
    make(&c, &list[i]);
    if (attacked(&c, c.king[b->stm == WHITE ? 0 : 1], c.stm)) continue;
    total += depth <= 1 ? 1 : perft(&c, depth - 1);
  }
  return total;
}

static int parseFen(Board* b, const char* f) {
  int r = 7, file = 0;
  memset(b, 0, sizeof *b);
  b->ep = -1;
  for (; *f && *f != ' '; f++) {
    if (*f == '/') { r--; file = 0; }
    else if (*f >= '1' && *f <= '8') file += *f - '0';
    else {
      const char* pcs = " PNBRQK  pnbrqk";
      const char* q = strchr(pcs, *f);
      if (!q || *f == ' ') return 0;
      int idx = (int)(q - pcs), col = idx >= 8 ? BLACK : WHITE, t = idx & 7, s = r * 16 + file++;
      b->sq[s] = col | t;
      if (t == K) b->king[col == WHITE ? 0 : 1] = s;
    }
  }
  while (*f == ' ') f++;
  b->stm = *f == 'b' ? BLACK : WHITE;
  if (*f) f++;
  while (*f == ' ') f++;
  for (; *f && *f != ' '; f++) {
    if (*f == 'K') b->castle |= 1;
    if (*f == 'Q') b->castle |= 2;
    if (*f == 'k') b->castle |= 4;
    if (*f == 'q') b->castle |= 8;
  }
  while (*f == ' ') f++;
  if (*f && *f != '-') b->ep = (f[1] - '1') * 16 + (f[0] - 'a');
  return 1;
}

typedef struct { Board b; int depth; uint64_t result; Mv m; } Task;
static Task* tasks;
static int ntasks, nextTask;
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;

static void* worker(void* arg) {
  (void)arg;
  for (;;) {
    pthread_mutex_lock(&mu);
    int i = nextTask++;
    pthread_mutex_unlock(&mu);
    if (i >= ntasks) return NULL;
    tasks[i].result = tasks[i].depth == 0 ? 1 : perft(&tasks[i].b, tasks[i].depth);
  }
}

int main(int argc, char** argv) {
  Board root;
  if (argc < 3 || !parseFen(&root, argv[1])) {
    fprintf(stderr, "usage: naive FEN depth [threads]\n");
    return 1;
  }
  int depth = atoi(argv[2]), threads = argc > 3 ? atoi(argv[3]) : 8, i;
  Mv list[256];
  int n = genPseudo(&root, list);
  tasks = calloc(256, sizeof *tasks);
  for (i = 0; i < n; i++) {
    Board c = root;
    make(&c, &list[i]);
    if (attacked(&c, c.king[root.stm == WHITE ? 0 : 1], c.stm)) continue;
    tasks[ntasks].b = c;
    tasks[ntasks].depth = depth - 1;
    tasks[ntasks].m = list[i];
    ntasks++;
  }
  pthread_t th[64];
  for (i = 0; i < threads; i++) pthread_create(&th[i], NULL, worker, NULL);
  for (i = 0; i < threads; i++) pthread_join(th[i], NULL);
  uint64_t total = 0;
  for (i = 0; i < ntasks; i++) {
    Mv* m = &tasks[i].m;
    printf("%c%c%c%c%s: %llu\n", 'a' + (m->from & 7), '1' + (m->from >> 4), 'a' + (m->to & 7), '1' + (m->to >> 4),
           m->promo ? (m->promo == Q ? "q" : m->promo == R ? "r" : m->promo == B ? "b" : "n") : "",
           (unsigned long long)tasks[i].result);
    total += tasks[i].result;
  }
  printf("total: %llu\n", (unsigned long long)total);
  return 0;
}
