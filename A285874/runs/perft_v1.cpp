// perft.cpp - multithreaded perft with a lockless, 128-bit-keyed transposition table.
//
// Written to extend OEIS A285874 (chess games after n plies, both sides starting
// without rooks) and its siblings A285873, A285875-A285878, validated against A048987.
//
// Method: enumerate the unique positions at a split ply k (with path multiplicities),
// then compute perft(n-k) of each on worker threads sharing one hash table, and sum
// multiplicity * count in 128-bit arithmetic. Leaves are bulk-counted (depth-1 nodes
// count legal moves without making them). Table entries are located by one Zobrist
// key and verified by a second, independent 64-bit key plus depth, so a false hit
// needs a 64-bit collision within a bucket.
//
// Build: clang++ -O3 -std=c++20 -mcpu=native -o perft perft.cpp

#include <pthread.h>
#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using U64 = uint64_t;
using U128 = unsigned __int128;

static inline int lsb(U64 b) { return __builtin_ctzll(b); }
static inline int popcnt(U64 b) { return __builtin_popcountll(b); }
static inline U64 bit(int s) { return 1ULL << s; }
static inline int popLsb(U64& b) {
  int s = lsb(b);
  b &= b - 1;
  return s;
}

enum { WHITE, BLACK };
enum { PAWN, KNIGHT, BISHOP, ROOK, QUEEN, KING, NO_PIECE };
enum { A1 = 0, B1, C1, D1, E1, F1, G1, H1, A8 = 56, B8, C8, D8, E8, F8, G8, H8 };
constexpr int NO_SQ = 64;

constexpr U64 FILE_A = 0x0101010101010101ULL, FILE_H = FILE_A << 7;
constexpr U64 RANK_1 = 0xFFULL, RANK_3 = RANK_1 << 16, RANK_6 = RANK_1 << 40, RANK_8 = RANK_1 << 56;

// ---------------------------------------------------------------- attack tables

static U64 KnightAtt[64], KingAtt[64], PawnAtt[2][64], Between[64][64], Line[64][64];

struct Magic {
  U64 mask, magic;
  U64* att;
  unsigned shift;
  unsigned index(U64 occ) const { return unsigned(((occ & mask) * magic) >> shift); }
};
static Magic RookM[64], BishopM[64];
static U64 RookTable[0x19000], BishopTable[0x1480];

static inline U64 rookAtt(int s, U64 occ) { return RookM[s].att[RookM[s].index(occ)]; }
static inline U64 bishopAtt(int s, U64 occ) { return BishopM[s].att[BishopM[s].index(occ)]; }

static const int RookDir[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
static const int BishopDir[4][2] = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};

static U64 slowSlide(int s, U64 occ, const int dir[4][2]) {
  U64 a = 0;
  for (int d = 0; d < 4; d++) {
    int r = s / 8 + dir[d][0], f = s % 8 + dir[d][1];
    while (r >= 0 && r < 8 && f >= 0 && f < 8) {
      a |= bit(r * 8 + f);
      if (occ & bit(r * 8 + f)) break;
      r += dir[d][0];
      f += dir[d][1];
    }
  }
  return a;
}

struct Rng {  // xorshift64*
  U64 s;
  explicit Rng(U64 seed) : s(seed ? seed : 0x9E3779B97F4A7C15ULL) {}
  U64 next() {
    s ^= s >> 12;
    s ^= s << 25;
    s ^= s >> 27;
    return s * 0x2545F4914F6CDD1DULL;
  }
};

static void initMagics(Magic* M, U64* table, const int dir[4][2]) {
  std::vector<U64> occs(4096), refs(4096);
  std::vector<int> epoch(4096, 0);
  int cnt = 0;
  Rng rng(0x5DEECE66DULL);
  U64* t = table;
  for (int s = 0; s < 64; s++) {
    const U64 rankS = RANK_1 << (8 * (s / 8)), fileS = FILE_A << (s % 8);
    const U64 edges = ((RANK_1 | RANK_8) & ~rankS) | ((FILE_A | FILE_H) & ~fileS);
    Magic& m = M[s];
    m.mask = slowSlide(s, 0, dir) & ~edges;
    m.shift = 64 - popcnt(m.mask);
    m.att = t;
    int size = 0;
    U64 b = 0;
    do {
      occs[size] = b;
      refs[size] = slowSlide(s, b, dir);
      size++;
      b = (b - m.mask) & m.mask;
    } while (b);
    t += size;
    for (int i = 0; i < size;) {
      do m.magic = rng.next() & rng.next() & rng.next();
      while (popcnt((m.magic * m.mask) >> 56) < 6);
      ++cnt;
      for (i = 0; i < size; i++) {
        const unsigned idx = m.index(occs[i]);
        if (epoch[idx] < cnt) {
          epoch[idx] = cnt;
          m.att[idx] = refs[i];
        } else if (m.att[idx] != refs[i]) {
          break;
        }
      }
    }
  }
}

static uint8_t CastleMask[64];

static void initTables() {
  for (int s = 0; s < 64; s++) {
    const int r = s / 8, f = s % 8;
    auto add = [&](U64& bb, int dr, int df) {
      const int rr = r + dr, ff = f + df;
      if (rr >= 0 && rr < 8 && ff >= 0 && ff < 8) bb |= bit(rr * 8 + ff);
    };
    static const int kn[8][2] = {{1, 2}, {2, 1}, {-1, 2}, {-2, 1}, {1, -2}, {2, -1}, {-1, -2}, {-2, -1}};
    for (auto& d : kn) add(KnightAtt[s], d[0], d[1]);
    for (int dr = -1; dr <= 1; dr++)
      for (int df = -1; df <= 1; df++)
        if (dr || df) add(KingAtt[s], dr, df);
    add(PawnAtt[WHITE][s], 1, -1);
    add(PawnAtt[WHITE][s], 1, 1);
    add(PawnAtt[BLACK][s], -1, -1);
    add(PawnAtt[BLACK][s], -1, 1);
  }
  initMagics(RookM, RookTable, RookDir);
  initMagics(BishopM, BishopTable, BishopDir);
  for (int a = 0; a < 64; a++)
    for (int b = 0; b < 64; b++) {
      if (a == b) continue;
      for (const auto* dir : {RookDir, BishopDir}) {
        if (slowSlide(a, 0, dir) & bit(b)) {
          Line[a][b] = (slowSlide(a, 0, dir) & slowSlide(b, 0, dir)) | bit(a) | bit(b);
          Between[a][b] = slowSlide(a, bit(b), dir) & slowSlide(b, bit(a), dir);
        }
      }
    }
  for (int s = 0; s < 64; s++) CastleMask[s] = 15;
  CastleMask[E1] = 15 & ~3;
  CastleMask[H1] = 15 & ~1;
  CastleMask[A1] = 15 & ~2;
  CastleMask[E8] = 15 & ~12;
  CastleMask[H8] = 15 & ~4;
  CastleMask[A8] = 15 & ~8;
}

// ---------------------------------------------------------------- position

// Two independent Zobrist key sets: key 0 picks the table bucket, key 1 verifies.
static U64 Z[2][2][6][64], ZSide[2], ZCastle[2][16], ZEp[2][8];

static void initZobrist(U64 seed) {
  Rng r(seed * 0x9E3779B97F4A7C15ULL + 0x632BE59BD9B4E019ULL);
  for (int k = 0; k < 2; k++) {
    for (int c = 0; c < 2; c++)
      for (int p = 0; p < 6; p++)
        for (int s = 0; s < 64; s++) Z[k][c][p][s] = r.next();
    ZSide[k] = r.next();
    for (int i = 0; i < 16; i++) ZCastle[k][i] = r.next();
    for (int i = 0; i < 8; i++) ZEp[k][i] = r.next();
  }
}

struct Pos {
  U64 bb[6];   // by piece type, both colors
  U64 cbb[2];  // by color
  U64 k1, k2;  // Zobrist keys
  uint8_t stm, ep, castle;  // ep is set only when an en passant capture is pseudo-legal
};

static void computeKeys(Pos& p) {
  U64 k[2] = {0, 0};
  for (int i = 0; i < 2; i++) {
    for (int c = 0; c < 2; c++)
      for (int pc = 0; pc < 6; pc++)
        for (U64 b = p.bb[pc] & p.cbb[c]; b;) k[i] ^= Z[i][c][pc][popLsb(b)];
    if (p.stm == BLACK) k[i] ^= ZSide[i];
    k[i] ^= ZCastle[i][p.castle];
    if (p.ep != NO_SQ) k[i] ^= ZEp[i][p.ep & 7];
  }
  p.k1 = k[0];
  p.k2 = k[1];
}

static bool parseFen(const std::string& fen, Pos& p) {
  memset(&p, 0, sizeof p);
  const char* s = fen.c_str();
  int r = 7, f = 0;
  for (; *s && *s != ' '; s++) {
    if (*s == '/') {
      r--;
      f = 0;
    } else if (*s >= '1' && *s <= '8') {
      f += *s - '0';
    } else {
      static const char* pcs = "PNBRQKpnbrqk";
      const char* q = strchr(pcs, *s);
      if (!q || r < 0 || f > 7) return false;
      const int idx = int(q - pcs), sq = r * 8 + f++;
      p.bb[idx % 6] |= bit(sq);
      p.cbb[idx / 6] |= bit(sq);
    }
  }
  while (*s == ' ') s++;
  p.stm = (*s == 'b') ? BLACK : WHITE;
  if (*s) s++;
  while (*s == ' ') s++;
  for (; *s && *s != ' '; s++) {
    if (*s == 'K') p.castle |= 1;
    if (*s == 'Q') p.castle |= 2;
    if (*s == 'k') p.castle |= 4;
    if (*s == 'q') p.castle |= 8;
  }
  // Drop rights whose king or rook is not on its original square.
  const U64 wr = p.bb[ROOK] & p.cbb[WHITE], br = p.bb[ROOK] & p.cbb[BLACK];
  const U64 wk = p.bb[KING] & p.cbb[WHITE], bk = p.bb[KING] & p.cbb[BLACK];
  if (!(wk & bit(E1)) || !(wr & bit(H1))) p.castle &= ~1;
  if (!(wk & bit(E1)) || !(wr & bit(A1))) p.castle &= ~2;
  if (!(bk & bit(E8)) || !(br & bit(H8))) p.castle &= ~4;
  if (!(bk & bit(E8)) || !(br & bit(A8))) p.castle &= ~8;
  while (*s == ' ') s++;
  p.ep = NO_SQ;
  if (*s && *s != '-' && s[1]) {
    const int sq = (s[1] - '1') * 8 + (s[0] - 'a');
    if (sq >= 0 && sq < 64 && (PawnAtt[p.stm ^ 1][sq] & p.bb[PAWN] & p.cbb[p.stm])) p.ep = uint8_t(sq);
  }
  if (popcnt(wk) != 1 || popcnt(bk) != 1) return false;
  computeKeys(p);
  return true;
}

// ---------------------------------------------------------------- moves

using Move = uint32_t;
enum { F_NORMAL, F_DOUBLE, F_EP, F_CASTLE };

static inline Move mkMove(int from, int to, int pc, int promo = NO_PIECE, int flag = F_NORMAL) {
  return Move(from | (to << 6) | (pc << 12) | (promo << 15) | (flag << 18));
}

struct MoveList {
  int n = 0;
  Move m[256];
};

template <int Us>
static inline U64 up(U64 b) { return Us == WHITE ? b << 8 : b >> 8; }
template <int Us>
static inline U64 upWest(U64 b) { return Us == WHITE ? (b & ~FILE_A) << 7 : (b & ~FILE_A) >> 9; }
template <int Us>
static inline U64 upEast(U64 b) { return Us == WHITE ? (b & ~FILE_H) << 9 : (b & ~FILE_H) >> 7; }

// Legal move generation (Gen=true fills ml) or counting (Gen=false returns the count).
template <int Us, bool Gen>
static inline int legalMoves(const Pos& p, MoveList* ml) {
  constexpr int Them = Us ^ 1;
  constexpr int Up = Us == WHITE ? 8 : -8;
  constexpr int West = Us == WHITE ? 7 : -9;
  constexpr int East = Us == WHITE ? 9 : -7;
  constexpr U64 Promo = Us == WHITE ? RANK_8 : RANK_1;
  constexpr U64 Third = Us == WHITE ? RANK_3 : RANK_6;

  const U64 us = p.cbb[Us], them = p.cbb[Them], occ = us | them;
  const int ksq = lsb(p.bb[KING] & us);
  const U64 tP = p.bb[PAWN] & them, tN = p.bb[KNIGHT] & them, tK = p.bb[KING] & them;
  const U64 tBQ = (p.bb[BISHOP] | p.bb[QUEEN]) & them;
  const U64 tRQ = (p.bb[ROOK] | p.bb[QUEEN]) & them;

  int n = 0;
  auto push = [&](Move m) {
    if constexpr (Gen) ml->m[ml->n++] = m;
  };
  auto emit = [&](int from, U64 to, int pc) {
    if constexpr (Gen) {
      while (to) push(mkMove(from, popLsb(to), pc));
    } else {
      n += popcnt(to);
    }
  };
  auto attacked = [&](int s, U64 o) {
    return ((PawnAtt[Us][s] & tP) | (KnightAtt[s] & tN) | (KingAtt[s] & tK) | (bishopAtt(s, o) & tBQ) |
            (rookAtt(s, o) & tRQ)) != 0;
  };

  const U64 checkers =
      (PawnAtt[Us][ksq] & tP) | (KnightAtt[ksq] & tN) | (bishopAtt(ksq, occ) & tBQ) | (rookAtt(ksq, occ) & tRQ);

  {  // king moves: test each destination with the king lifted so sliders see through it
    const U64 o = occ ^ bit(ksq);
    U64 cand = KingAtt[ksq] & ~us, ok = 0;
    while (cand) {
      const int s = popLsb(cand);
      if (!attacked(s, o)) ok |= bit(s);
    }
    emit(ksq, ok, KING);
  }
  if (checkers & (checkers - 1)) return Gen ? ml->n : n;

  const U64 target = checkers ? (Between[ksq][lsb(checkers)] | checkers) : ~us;

  U64 pinned = 0;
  for (U64 sn = (rookAtt(ksq, them) & tRQ) | (bishopAtt(ksq, them) & tBQ); sn;) {
    const U64 b = Between[ksq][popLsb(sn)] & occ;
    if (b && !(b & (b - 1))) pinned |= b & us;
  }

  for (U64 b = p.bb[KNIGHT] & us & ~pinned; b;) {
    const int s = popLsb(b);
    emit(s, KnightAtt[s] & target, KNIGHT);
  }
  for (U64 b = (p.bb[BISHOP] | p.bb[QUEEN]) & us; b;) {
    const int s = popLsb(b);
    U64 a = bishopAtt(s, occ) & target;
    if (pinned & bit(s)) a &= Line[ksq][s];
    emit(s, a, (p.bb[QUEEN] & bit(s)) ? QUEEN : BISHOP);
  }
  for (U64 b = (p.bb[ROOK] | p.bb[QUEEN]) & us; b;) {
    const int s = popLsb(b);
    U64 a = rookAtt(s, occ) & target;
    if (pinned & bit(s)) a &= Line[ksq][s];
    emit(s, a, (p.bb[QUEEN] & bit(s)) ? QUEEN : ROOK);
  }

  const U64 pawns = p.bb[PAWN] & us, empty = ~occ;
  {  // unpinned pawns, set-wise
    const U64 fr = pawns & ~pinned;
    const U64 one = up<Us>(fr) & empty;
    const U64 two = up<Us>(one & Third) & empty & target;
    const U64 s1 = one & target;
    const U64 cw = upWest<Us>(fr) & them & target;
    const U64 ce = upEast<Us>(fr) & them & target;
    if constexpr (Gen) {
      auto gen = [&](U64 to, int delta, int flag) {
        while (to) {
          const int t = popLsb(to), f = t - delta;
          if (bit(t) & Promo) {
            for (int pr = QUEEN; pr >= KNIGHT; pr--) push(mkMove(f, t, PAWN, pr));
          } else {
            push(mkMove(f, t, PAWN, NO_PIECE, flag));
          }
        }
      };
      gen(s1, Up, F_NORMAL);
      gen(two, 2 * Up, F_DOUBLE);
      gen(cw, West, F_NORMAL);
      gen(ce, East, F_NORMAL);
    } else {
      n += popcnt(s1 & ~Promo) + 4 * popcnt(s1 & Promo) + popcnt(two) + popcnt(cw & ~Promo) +
           4 * popcnt(cw & Promo) + popcnt(ce & ~Promo) + 4 * popcnt(ce & Promo);
    }
  }
  if (!checkers) {  // pinned pawns move only along the pin line (and never when in check)
    for (U64 b = pawns & pinned; b;) {
      const int s = popLsb(b);
      const U64 one = up<Us>(bit(s)) & empty;
      const U64 two = up<Us>(one & Third) & empty;
      U64 to = (one | two | (PawnAtt[Us][s] & them)) & Line[ksq][s];
      if constexpr (Gen) {
        while (to) {
          const int t = popLsb(to);
          if (bit(t) & Promo) {
            for (int pr = QUEEN; pr >= KNIGHT; pr--) push(mkMove(s, t, PAWN, pr));
          } else {
            push(mkMove(s, t, PAWN, NO_PIECE, (t - s == 2 * Up) ? F_DOUBLE : F_NORMAL));
          }
        }
      } else {
        n += popcnt(to & ~Promo) + 4 * popcnt(to & Promo);
      }
    }
  }
  if (p.ep != NO_SQ) {  // en passant: verify by testing the king against the resulting occupancy
    const int ep = p.ep, cap = ep - Up;
    for (U64 b = PawnAtt[Them][ep] & pawns; b;) {
      const int s = popLsb(b);
      const U64 o = occ ^ bit(s) ^ bit(ep) ^ bit(cap);
      const bool exposed = ((PawnAtt[Us][ksq] & tP & ~bit(cap)) | (KnightAtt[ksq] & tN) |
                            (bishopAtt(ksq, o) & tBQ) | (rookAtt(ksq, o) & tRQ)) != 0;
      if (!exposed) {
        if constexpr (Gen) push(mkMove(s, ep, PAWN, NO_PIECE, F_EP));
        else n++;
      }
    }
  }
  if (!checkers && p.castle) {
    constexpr int KS = Us == WHITE ? 1 : 4, QS = Us == WHITE ? 2 : 8;
    constexpr int E = Us == WHITE ? E1 : E8;
    const U64 ourR = p.bb[ROOK] & us;
    if ((p.castle & KS) && ksq == E && (ourR & bit(E + 3)) && !(occ & (bit(E + 1) | bit(E + 2))) &&
        !attacked(E + 1, occ) && !attacked(E + 2, occ)) {
      if constexpr (Gen) push(mkMove(E, E + 2, KING, NO_PIECE, F_CASTLE));
      else n++;
    }
    if ((p.castle & QS) && ksq == E && (ourR & bit(E - 4)) && !(occ & (bit(E - 1) | bit(E - 2) | bit(E - 3))) &&
        !attacked(E - 1, occ) && !attacked(E - 2, occ)) {
      if constexpr (Gen) push(mkMove(E, E - 2, KING, NO_PIECE, F_CASTLE));
      else n++;
    }
  }
  return Gen ? ml->n : n;
}

template <int Us>
static inline void makeMove(const Pos& p, Pos& n, Move m) {
  constexpr int Them = Us ^ 1;
  constexpr int Up = Us == WHITE ? 8 : -8;
  const int from = m & 63, to = (m >> 6) & 63, pc = (m >> 12) & 7, promo = (m >> 15) & 7, flag = (m >> 18) & 3;
  const U64 fb = bit(from), tb = bit(to);
  n = p;
  U64 k1 = p.k1, k2 = p.k2;
  if (p.cbb[Them] & tb) {
    int cap = PAWN;
    while (!(p.bb[cap] & tb)) cap++;
    n.bb[cap] ^= tb;
    n.cbb[Them] ^= tb;
    k1 ^= Z[0][Them][cap][to];
    k2 ^= Z[1][Them][cap][to];
  }
  const int placed = promo != NO_PIECE ? promo : pc;
  n.bb[pc] ^= fb;
  n.bb[placed] ^= tb;
  n.cbb[Us] ^= fb | tb;
  k1 ^= Z[0][Us][pc][from] ^ Z[0][Us][placed][to];
  k2 ^= Z[1][Us][pc][from] ^ Z[1][Us][placed][to];
  if (flag == F_EP) {
    const int c = to - Up;
    n.bb[PAWN] ^= bit(c);
    n.cbb[Them] ^= bit(c);
    k1 ^= Z[0][Them][PAWN][c];
    k2 ^= Z[1][Them][PAWN][c];
  } else if (flag == F_CASTLE) {
    const int rf = to > from ? to + 1 : to - 2, rt = to > from ? to - 1 : to + 1;
    n.bb[ROOK] ^= bit(rf) | bit(rt);
    n.cbb[Us] ^= bit(rf) | bit(rt);
    k1 ^= Z[0][Us][ROOK][rf] ^ Z[0][Us][ROOK][rt];
    k2 ^= Z[1][Us][ROOK][rf] ^ Z[1][Us][ROOK][rt];
  }
  const int nc = p.castle & CastleMask[from] & CastleMask[to];
  if (nc != p.castle) {
    n.castle = uint8_t(nc);
    k1 ^= ZCastle[0][p.castle] ^ ZCastle[0][nc];
    k2 ^= ZCastle[1][p.castle] ^ ZCastle[1][nc];
  }
  if (p.ep != NO_SQ) {
    n.ep = NO_SQ;
    k1 ^= ZEp[0][p.ep & 7];
    k2 ^= ZEp[1][p.ep & 7];
  }
  if (flag == F_DOUBLE) {
    const int e = from + Up;
    if (PawnAtt[Us][e] & p.bb[PAWN] & p.cbb[Them]) {
      n.ep = uint8_t(e);
      k1 ^= ZEp[0][e & 7];
      k2 ^= ZEp[1][e & 7];
    }
  }
  n.stm = Them;
  n.k1 = k1 ^ ZSide[0];
  n.k2 = k2 ^ ZSide[1];
}

static std::string moveStr(Move m) {
  std::string s;
  const int from = m & 63, to = (m >> 6) & 63, promo = (m >> 15) & 7;
  s += char('a' + from % 8);
  s += char('1' + from / 8);
  s += char('a' + to % 8);
  s += char('1' + to / 8);
  if (promo != NO_PIECE) s += "pnbrqk"[promo];
  return s;
}

// ---------------------------------------------------------------- transposition table

struct TT {
  struct Entry {
    std::atomic<U64> lock, data;  // lock = k2 ^ data; data = count << 8 | depth
  };
  static constexpr int W = 8;  // entries per 128-byte bucket (Apple Silicon cache line)
  Entry* t = nullptr;
  U64 mask = 0;
  size_t bytes = 0;

  void init(size_t want) {
    U64 nb = 1;
    while (nb * 2 * W * sizeof(Entry) <= want) nb *= 2;
    bytes = nb * W * sizeof(Entry);
    void* mem = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
    if (mem == MAP_FAILED) {
      perror("mmap");
      exit(1);
    }
    t = static_cast<Entry*>(mem);
    mask = nb - 1;
  }
  ~TT() {
    if (t) munmap(t, bytes);
  }
  void prefetch(U64 k1) const { __builtin_prefetch(&t[(k1 & mask) * W]); }
  bool probe(U64 k1, U64 k2, int depth, U64& out) const {
    const Entry* b = &t[(k1 & mask) * W];
    for (int i = 0; i < W; i++) {
      const U64 d = b[i].data.load(std::memory_order_relaxed);
      const U64 l = b[i].lock.load(std::memory_order_relaxed);
      if ((l ^ d) == k2 && int(d & 0xFF) == depth) {
        out = d >> 8;
        return true;
      }
    }
    return false;
  }
  void store(U64 k1, U64 k2, int depth, U64 count) {
    if (count >> 56) return;
    const U64 nd = (count << 8) | U64(depth);
    Entry* b = &t[(k1 & mask) * W];
    int victim = 0;
    U64 best = ~0ULL;
    for (int i = 0; i < W; i++) {
      const U64 d = b[i].data.load(std::memory_order_relaxed);
      const U64 l = b[i].lock.load(std::memory_order_relaxed);
      if ((l ^ d) == k2 && int(d & 0xFF) == depth) return;
      const U64 score = ((d & 0xFF) << 56) | (d >> 8);  // replace shallowest, then smallest subtree
      if (score < best) {
        best = score;
        victim = i;
      }
    }
    b[victim].lock.store(k2 ^ nd, std::memory_order_relaxed);
    b[victim].data.store(nd, std::memory_order_relaxed);
  }
};

// ---------------------------------------------------------------- perft

template <int Us>
static U64 perft(const Pos& p, int depth, TT* tt, int ttMin) {
  constexpr int Them = Us ^ 1;
  if (depth == 1) return U64(legalMoves<Us, false>(p, nullptr));
  const bool useTT = tt && depth >= ttMin;
  U64 r;
  if (useTT && tt->probe(p.k1, p.k2, depth, r)) return r;
  MoveList ml;
  legalMoves<Us, true>(p, &ml);
  r = 0;
  if (depth == 2) {
    Pos c;
    for (int i = 0; i < ml.n; i++) {
      makeMove<Us>(p, c, ml.m[i]);
      r += U64(legalMoves<Them, false>(c, nullptr));
    }
  } else {
    Pos kids[256];
    const bool pf = tt && depth - 1 >= ttMin;
    for (int i = 0; i < ml.n; i++) {
      makeMove<Us>(p, kids[i], ml.m[i]);
      if (pf) tt->prefetch(kids[i].k1);
    }
    for (int i = 0; i < ml.n; i++) r += perft<Them>(kids[i], depth - 1, tt, ttMin);
  }
  if (useTT) tt->store(p.k1, p.k2, depth, r);
  return r;
}

static U64 perftAny(const Pos& p, int depth, TT* tt, int ttMin) {
  if (depth == 0) return 1;
  return p.stm == WHITE ? perft<WHITE>(p, depth, tt, ttMin) : perft<BLACK>(p, depth, tt, ttMin);
}

static bool kingAttacked(const Pos& p, int c) {
  const int k = lsb(p.bb[KING] & p.cbb[c]);
  const U64 occ = p.cbb[0] | p.cbb[1], them = p.cbb[c ^ 1];
  return (PawnAtt[c][k] & p.bb[PAWN] & them) || (KnightAtt[k] & p.bb[KNIGHT] & them) ||
         (KingAtt[k] & p.bb[KING] & them) || (bishopAtt(k, occ) & (p.bb[BISHOP] | p.bb[QUEEN]) & them) ||
         (rookAtt(k, occ) & (p.bb[ROOK] | p.bb[QUEEN]) & them);
}

// Reference perft: makes every move, no bulk counting, no table. Cross-checks the
// generator against the counter, incremental keys against recomputed ones, and that
// no generated move leaves the mover's king attacked.
template <int Us>
static U64 perftSlow(const Pos& p, int depth) {
  if (depth == 0) return 1;
  MoveList ml;
  legalMoves<Us, true>(p, &ml);
  if (legalMoves<Us, false>(p, nullptr) != ml.n) {
    fprintf(stderr, "FATAL: generator/counter mismatch\n");
    exit(2);
  }
  U64 r = 0;
  for (int i = 0; i < ml.n; i++) {
    Pos c;
    makeMove<Us>(p, c, ml.m[i]);
    Pos chk = c;
    computeKeys(chk);
    if (chk.k1 != c.k1 || chk.k2 != c.k2) {
      fprintf(stderr, "FATAL: incremental key mismatch after %s\n", moveStr(ml.m[i]).c_str());
      exit(2);
    }
    if (kingAttacked(c, Us)) {
      fprintf(stderr, "FATAL: illegal move generated: %s\n", moveStr(ml.m[i]).c_str());
      exit(2);
    }
    r += perftSlow<Us ^ 1>(c, depth - 1);
  }
  return r;
}

static U64 perftSlowAny(const Pos& p, int depth) {
  return p.stm == WHITE ? perftSlow<WHITE>(p, depth) : perftSlow<BLACK>(p, depth);
}

// ---------------------------------------------------------------- split-ply frontier

struct PosHash {
  size_t operator()(const Pos& p) const { return size_t(p.k1); }
};
struct PosEq {
  bool operator()(const Pos& a, const Pos& b) const {
    return !memcmp(a.bb, b.bb, sizeof a.bb) && a.cbb[0] == b.cbb[0] && a.cbb[1] == b.cbb[1] && a.stm == b.stm &&
           a.ep == b.ep && a.castle == b.castle;
  }
};

struct Frontier {
  std::vector<Pos> pos;
  std::vector<U64> mult;
  std::unordered_map<Pos, uint32_t, PosHash, PosEq> index;
};

template <int Us>
static void enumerate(const Pos& p, int ply, Frontier& F) {
  if (ply == 0) {
    auto [it, inserted] = F.index.try_emplace(p, uint32_t(F.pos.size()));
    if (inserted) {
      F.pos.push_back(p);
      F.mult.push_back(1);
    } else {
      F.mult[it->second]++;
    }
    return;
  }
  MoveList ml;
  legalMoves<Us, true>(p, &ml);
  for (int i = 0; i < ml.n; i++) {
    Pos c;
    makeMove<Us>(p, c, ml.m[i]);
    enumerate<Us ^ 1>(c, ply - 1, F);
  }
}

// ---------------------------------------------------------------- parallel driver

static std::string u128str(U128 v) {
  if (v == 0) return "0";
  std::string s;
  while (v) {
    s += char('0' + int(v % 10));
    v /= 10;
  }
  std::reverse(s.begin(), s.end());
  return s;
}

static double now() {
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}

struct Job {
  const std::vector<Pos>* pos = nullptr;
  int depth = 0, ttMin = 2;
  TT* tt = nullptr;
  std::vector<U64> result;
  std::unique_ptr<std::atomic<uint8_t>[]> done;
  std::atomic<size_t> next{0}, completed{0};
};

static void* workerMain(void* arg) {
  Job* j = static_cast<Job*>(arg);
  const size_t N = j->pos->size();
  for (;;) {
    const size_t i = j->next.fetch_add(1);
    if (i >= N) break;
    if (j->done[i].load(std::memory_order_acquire)) continue;
    j->result[i] = perftAny((*j->pos)[i], j->depth, j->tt, j->ttMin);
    j->done[i].store(1, std::memory_order_release);
    j->completed.fetch_add(1);
  }
  return nullptr;
}

struct RunCfg {
  int split = -1, threads = 0, ttMin = 2;
  double reportEvery = 10;
  bool quiet = false;
  std::string ckpt, fen;
  U64 seed = 1;
};

static U128 runPerft(const Pos& root, int depth, TT* tt, const RunCfg& cfg) {
  if (depth == 0) return 1;
  int split = cfg.split >= 0 ? cfg.split : std::clamp(depth - 6, 0, 6);
  split = std::min(split, depth - 1);
  const double t0 = now();

  Frontier F;
  if (root.stm == WHITE) enumerate<WHITE>(root, split, F);
  else enumerate<BLACK>(root, split, F);
  F.index.clear();
  const size_t N = F.pos.size();
  U128 paths = 0;
  for (U64 m : F.mult) paths += m;
  if (!cfg.quiet)
    printf("split ply %d: %zu unique positions (%s paths), %.1fs\n", split, N, u128str(paths).c_str(), now() - t0);

  Job job;
  job.pos = &F.pos;
  job.depth = depth - split;
  job.tt = tt;
  job.ttMin = cfg.ttMin;
  job.result.assign(N, 0);
  job.done.reset(new std::atomic<uint8_t>[N]);
  for (size_t i = 0; i < N; i++) job.done[i].store(0);

  // Resume from checkpoint: lines "index key1 count", key1 re-verified against the task.
  std::vector<uint8_t> written(N, 0);
  FILE* ck = nullptr;
  const std::string header = "# fen=" + cfg.fen + " depth=" + std::to_string(depth) + " split=" +
                             std::to_string(split) + " seed=" + std::to_string(cfg.seed);
  if (!cfg.ckpt.empty()) {
    size_t resumed = 0;
    if (FILE* in = fopen(cfg.ckpt.c_str(), "r")) {
      char line[1024];
      if (!fgets(line, sizeof line, in) || std::string(line) != header + "\n") {
        fprintf(stderr, "checkpoint %s does not match this run\n", cfg.ckpt.c_str());
        exit(1);
      }
      size_t i;
      U64 k, c;
      while (fscanf(in, "%zu %" SCNx64 " %" SCNu64, &i, &k, &c) == 3) {
        if (i >= N || F.pos[i].k1 != k) {
          fprintf(stderr, "checkpoint entry %zu does not match task list\n", i);
          exit(1);
        }
        job.result[i] = c;
        job.done[i].store(1);
        written[i] = 1;
        resumed++;
      }
      fclose(in);
      ck = fopen(cfg.ckpt.c_str(), "a");
    } else {
      ck = fopen(cfg.ckpt.c_str(), "w");
      fprintf(ck, "%s\n", header.c_str());
      fflush(ck);
    }
    job.completed = resumed;
    if (resumed && !cfg.quiet) printf("resumed %zu/%zu tasks from %s\n", resumed, N, cfg.ckpt.c_str());
  }

  const int T = cfg.threads > 0 ? cfg.threads : int(std::thread::hardware_concurrency());
  std::vector<pthread_t> th(T);
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setstacksize(&attr, 64 << 20);
  for (int i = 0; i < T; i++) pthread_create(&th[i], &attr, workerMain, &job);

  const double t1 = now();
  double lastReport = t1;
  auto flushCkpt = [&] {
    if (!ck) return;
    for (size_t i = 0; i < N; i++)
      if (!written[i] && job.done[i].load(std::memory_order_acquire)) {
        fprintf(ck, "%zu %" PRIx64 " %" PRIu64 "\n", i, F.pos[i].k1, job.result[i]);
        written[i] = 1;
      }
    fflush(ck);
    fsync(fileno(ck));
  };
  while (job.completed.load() < N) {
    usleep(200000);
    const double t = now();
    if (t - lastReport >= cfg.reportEvery) {
      lastReport = t;
      flushCkpt();
      if (!cfg.quiet) {
        const size_t c = job.completed.load();
        printf("  %zu/%zu tasks (%.2f%%)  %.0fs elapsed  ETA %.0fs\n", c, N, 100.0 * c / N, t - t0,
               c ? (t - t1) * (N - c) / c : 0.0);
        fflush(stdout);
      }
    }
  }
  for (int i = 0; i < T; i++) pthread_join(th[i], nullptr);
  flushCkpt();
  if (ck) fclose(ck);

  U128 total = 0;
  for (size_t i = 0; i < N; i++) total += U128(F.mult[i]) * job.result[i];
  return total;
}

// ---------------------------------------------------------------- self-test

struct Case {
  const char* name;
  const char* fen;
  std::vector<U64> expect;  // expect[d] = perft(d)
};

static const char* const FEN_START = "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1";

static const std::vector<Case> kVariants = {
    {"std (A048987)", FEN_START,
     {1, 20, 400, 8902, 197281, 4865609, 119060324, 3195901860ULL, 84998978956ULL, 2439530234167ULL,
      69352859712417ULL, 2097651003696806ULL, 62854969236701747ULL, 1981066775000396239ULL}},
    {"noqueens (A285873)", "rnb1kbnr/pppppppp/8/8/8/8/PPPPPPPP/RNB1KBNR w KQkq - 0 1",
     {1, 21, 441, 9872, 220447, 5247292, 124278971, 3113440755ULL, 77520962327ULL, 2024021927610ULL}},
    {"norooks (A285874)", "1nbqkbn1/pppppppp/8/8/8/8/PPPPPPPP/1NBQKBN1 w - - 0 1",
     {1, 20, 400, 8702, 188473, 4505624, 106770421, 2770746488ULL, 71151220765ULL, 1969755500063ULL}},
    {"noknights (A285875)", "r1bqkb1r/pppppppp/8/8/8/8/PPPPPPPP/R1BQKB1R w KQkq - 0 1",
     {1, 18, 324, 6572, 132640, 3030492, 68633066, 1733220521ULL, 43321058602ULL, 1182486223832ULL}},
    {"nobishops (A285876)", "rn1qk1nr/pppppppp/8/8/8/8/PPPPPPPP/RN1QK1NR w KQkq - 0 1",
     {1, 22, 484, 11248, 260904, 6434922, 158069690, 4126252938ULL, 107097735673ULL, 2940365284820ULL}},
    {"nopawns (A285877)", "rnbqkbnr/8/8/8/8/8/8/RNBQKBNR w KQkq - 0 1",
     {1, 50, 2125, 96062, 4200525, 191462298, 8509434855ULL, 390020597683ULL}},
    {"pawnsking (A285878)", "4k3/pppppppp/8/8/8/8/PPPPPPPP/4K3 w - - 0 1",
     {1, 18, 324, 5658, 98766, 1683597, 28677387, 479763588ULL, 8014917042ULL, 132060434889ULL,
      2170519816231ULL}},
};

static const std::vector<Case> kSuite = {
    {"kiwipete", "r3k2r/p1ppqpb1/bn2pnp1/3PN3/1p2P3/2N2Q1p/PPPBBPPP/R3K2R w KQkq - 0 1",
     {1, 48, 2039, 97862, 4085603, 193690690, 8031647685ULL}},
    {"cpw3", "8/2p5/3p4/KP5r/1R3p1k/8/4P1P1/8 w - - 0 1",
     {1, 14, 191, 2812, 43238, 674624, 11030083, 178633661, 3009794393ULL}},
    {"cpw4", "r3k2r/Pppp1ppp/1b3nbN/nP6/BBP1P3/q4N2/Pp1P2PP/R2Q1RK1 w kq - 0 1",
     {1, 6, 264, 9467, 422333, 15833292, 706045033}},
    {"cpw4m", "r2q1rk1/pP1p2pp/Q4n2/bbp1p3/Np6/1B3NBn/pPPP1PPP/R3K2R b KQ - 0 1",
     {1, 6, 264, 9467, 422333, 15833292, 706045033}},
    {"cpw5", "rnbq1k1r/pp1Pbppp/2p5/8/2B5/8/PPP1NnPP/RNBQK2R w KQ - 1 8", {1, 44, 1486, 62379, 2103487, 89941194}},
    {"cpw6", "r4rk1/1pp1qppp/p1np1n2/2b1p1B1/2B1P1b1/P1NP1N2/1PP1QPPP/R4RK1 w - - 0 10",
     {1, 46, 2079, 89890, 3894594, 164075551, 6923051137ULL, 287188994746ULL}},
};

static int selfTest(TT* tt, int threads, U64 maxNodes) {
  int fails = 0;
  RunCfg cfg;
  cfg.threads = threads;
  cfg.quiet = true;
  auto check = [&](const char* what, const char* name, int d, U128 got, U64 want) {
    const bool ok = got == U128(want);
    if (!ok) fails++;
    printf("  %-5s %-20s d=%-2d %s %22s%s\n", ok ? "ok" : "FAIL", name, d, what, u128str(got).c_str(),
           ok ? "" : ("  expected " + std::to_string(want)).c_str());
    fflush(stdout);
  };
  for (const auto* list : {&kSuite, &kVariants}) {
    for (const Case& c : *list) {
      Pos p;
      if (!parseFen(c.fen, p)) {
        printf("bad fen %s\n", c.fen);
        return 1;
      }
      for (int d = 1; d < int(c.expect.size()); d++) {
        const U64 want = c.expect[d];
        if (want > maxNodes) break;
        if (want <= 5000000) check("slow", c.name, d, perftSlowAny(p, d), want);
        if (want <= 2000000000ULL) {
          cfg.split = -1;
          check("noTT", c.name, d, runPerft(p, d, nullptr, cfg), want);
        }
        check("TT  ", c.name, d, runPerft(p, d, tt, cfg), want);
      }
    }
  }
  printf(fails ? "SELF-TEST FAILED: %d failures\n" : "SELF-TEST PASSED\n", fails);
  return fails ? 1 : 0;
}

// ---------------------------------------------------------------- main

static void usage() {
  printf(
      "usage: perft [--variant NAME | --fen FEN] --depth N [--to M] [options]\n"
      "       perft --test [--max NODES]\n"
      "variants: std noqueens norooks noknights nobishops nopawns pawnsking\n"
      "options: --threads T  --hash MB  --split K  --ttmin D  --seed S  --notable\n"
      "         --ckpt FILE  --divide  --report SECONDS\n");
}

int main(int argc, char** argv) {
  initTables();
  std::string fen = FEN_START;
  int depth = -1, depthTo = -1, threads = 0, split = -1, ttMin = 2;
  size_t hashMB = 1024;
  U64 seed = 1, maxNodes = 300000000000ULL;
  bool test = false, noTable = false, divide = false;
  double report = 10;
  std::string ckpt;
  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    auto val = [&]() -> std::string {
      if (i + 1 >= argc) {
        usage();
        exit(1);
      }
      return argv[++i];
    };
    if (a == "--fen") fen = val();
    else if (a == "--variant") {
      const std::string v = val();
      bool found = false;
      for (const Case& c : kVariants)
        if (std::string(c.name).rfind(v + " ", 0) == 0) {
          fen = c.fen;
          found = true;
        }
      if (!found) {
        fprintf(stderr, "unknown variant %s\n", v.c_str());
        return 1;
      }
    } else if (a == "--depth") depth = std::stoi(val());
    else if (a == "--to") depthTo = std::stoi(val());
    else if (a == "--threads") threads = std::stoi(val());
    else if (a == "--hash") hashMB = std::stoull(val());
    else if (a == "--split") split = std::stoi(val());
    else if (a == "--ttmin") ttMin = std::stoi(val());
    else if (a == "--seed") seed = std::stoull(val());
    else if (a == "--max") maxNodes = std::stoull(val());
    else if (a == "--ckpt") ckpt = val();
    else if (a == "--report") report = std::stod(val());
    else if (a == "--notable") noTable = true;
    else if (a == "--divide") divide = true;
    else if (a == "--test") test = true;
    else {
      usage();
      return 1;
    }
  }
  initZobrist(seed);

  TT table;
  TT* tt = nullptr;
  if (!noTable) {
    table.init(hashMB << 20);
    tt = &table;
  }
  if (test) return selfTest(tt, threads, maxNodes);
  if (depth < 0) {
    usage();
    return 1;
  }
  Pos root;
  if (!parseFen(fen, root)) {
    fprintf(stderr, "bad FEN\n");
    return 1;
  }
  printf("fen: %s\nthreads: %d  table: %s  seed: %" PRIu64 "  ttmin: %d\n", fen.c_str(),
         threads > 0 ? threads : int(std::thread::hardware_concurrency()),
         tt ? (std::to_string(table.bytes >> 20) + " MB").c_str() : "off", seed, ttMin);
  fflush(stdout);

  RunCfg cfg;
  cfg.split = split;
  cfg.threads = threads;
  cfg.ttMin = ttMin;
  cfg.reportEvery = report;
  cfg.fen = fen;
  cfg.seed = seed;
  if (divide) {
    MoveList ml;
    if (root.stm == WHITE) legalMoves<WHITE, true>(root, &ml);
    else legalMoves<BLACK, true>(root, &ml);
    cfg.quiet = true;
    U128 total = 0;
    for (int i = 0; i < ml.n; i++) {
      Pos c;
      if (root.stm == WHITE) makeMove<WHITE>(root, c, ml.m[i]);
      else makeMove<BLACK>(root, c, ml.m[i]);
      const U128 r = depth > 0 ? runPerft(c, depth - 1, tt, cfg) : 1;
      total += r;
      printf("%s: %s\n", moveStr(ml.m[i]).c_str(), u128str(r).c_str());
      fflush(stdout);
    }
    printf("total: %s\n", u128str(total).c_str());
    return 0;
  }
  for (int d = depth; d <= std::max(depth, depthTo); d++) {
    RunCfg c = cfg;
    if (!ckpt.empty()) c.ckpt = depthTo > depth ? ckpt + ".d" + std::to_string(d) : ckpt;
    const double t0 = now();
    const U128 r = runPerft(root, d, tt, c);
    const double dt = now() - t0;
    printf("perft(%d) = %s   (%.2fs, %.3g leaves/s)\n", d, u128str(r).c_str(), dt, double(r) / std::max(dt, 1e-9));
    fflush(stdout);
  }
  return 0;
}
