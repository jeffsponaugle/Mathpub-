// zverify.cpp -- independent check of space-group counts for single Z-classes.
//
// Shares no code with CARAT.  Reads a CARAT bravais_TYP file (generators of a
// finite G <= GL_n(Z) plus generators of its normalizer N in GL_n(Z)) and
// computes from scratch:
//
//   * G, by closure of the generators;
//   * H = H^1(G, Q^n/Z^n) from the full Cayley graph (no presentation needed):
//     a cocycle is t = (t_1..t_k) in Q^(nk) (values on the generators); along a
//     BFS tree t(g g_i) = t(g) + g t_i defines t(g) = T_g t, and every
//     non-tree edge g -> h = g g_i gives the condition (T_g + g E_i - T_h) t
//     in Z^n.  With C the matrix of all conditions and U C V = diag(S),
//         H = { t : C t integral } / (ker C + Z^(nk))  =  (+)_j  Z/S_jj ;
//   * the action of n in N on H:  (n.t)(g_i) = n t(n^-1 g_i n);
//   * the number of N-orbits on H (space-group types = affine classes) and of
//     N-orbits on H x {+1,-1} with n.(x,s) = (n.x, det(n) s), which equals the
//     number of orientation-preserving classes (A006227 counts);
//     pairs = proper - affine (A395859 counts).
//
// Checks performed along the way: every normalizer generator conjugates G
// into itself; every computed action is an automorphism of H; the point
// group itself acts trivially on H; Burnside sums are divisible.
//
// usage: zverify [-x=k] [-y=n] file...    (explicit orbits if |H| <= 2^k)
//        a file argument "A+B" takes the generators from A and the union of the
//        normalizer generators of A and B
// output: file  |G|  H  affine  proper  pairs  improper  method

#include <gmpxx.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using std::string;
using std::vector;
#include "carat_io.hpp"

static int EXPLICIT_LOG2 = 24;
static long MAX_IMAGE = 3000000;

// ---------------------------------------------------------------- lattice insertion (HNF)

struct RowLattice {
    int nc;
    vector<vector<i128>> piv;  // piv[c] = basis row with pivot at column c (empty if none)
    explicit RowLattice(int n) : nc(n), piv(n) {}
    static bool big(i128 x) { return x > ((i128)1 << 100) || x < -((i128)1 << 100); }
    void insert(vector<i128> v, const string &file)
    {
        for (int c = 0; c < nc; c++) {
            if (v[c] == 0) continue;
            if (piv[c].empty()) {
                if (v[c] < 0) for (auto &x : v) x = -x;
                piv[c] = v;
                reduce_row(c);
                return;
            }
            vector<i128> &b = piv[c];
            i128 p = b[c], q = v[c];
            if (q % p == 0) {
                i128 f = q / p;
                for (int j = c; j < nc; j++) { v[j] -= f * b[j]; if (big(v[j])) fail(file, "HNF overflow"); }
                continue;
            }
            // extended gcd: u p + w q = g
            i128 r0 = p, r1 = q, s0 = 1, s1 = 0, t0 = 0, t1 = 1;
            while (r1 != 0) {
                i128 qq = r0 / r1, tmp;
                tmp = r0 - qq * r1; r0 = r1; r1 = tmp;
                tmp = s0 - qq * s1; s0 = s1; s1 = tmp;
                tmp = t0 - qq * t1; t0 = t1; t1 = tmp;
            }
            i128 g = r0, u = s0, w = t0;
            if (g < 0) { g = -g; u = -u; w = -w; }
            i128 pg = p / g, qg = q / g;
            vector<i128> nb(nc), nv(nc);
            for (int j = c; j < nc; j++) {
                nb[j] = u * b[j] + w * v[j];
                nv[j] = pg * v[j] - qg * b[j];
                if (big(nb[j]) || big(nv[j])) fail(file, "HNF overflow");
            }
            b = nb;
            reduce_row(c);
            v = nv;
        }
    }
    // reduce entries of row c at later pivot columns modulo those pivots
    void reduce_row(int c)
    {
        vector<i128> &b = piv[c];
        for (int j = c + 1; j < nc; j++) {
            if (piv[j].empty() || b[j] == 0) continue;
            i128 p = piv[j][j];
            i128 f = b[j] / p;
            if (b[j] - f * p < 0) f -= 1;
            if (f == 0) continue;
            for (int t = j; t < nc; t++) b[t] -= f * piv[j][t];
        }
    }
};

// ---------------------------------------------------------------- SNF with V, V^-1 (GMP)

struct SNF {
    vector<mpz_class> diag;            // nonzero elementary divisors (rank r)
    vector<vector<mpz_class>> V, Vi;   // nk x nk, U W V = S, Vi = V^-1
};

static SNF snf(vector<vector<mpz_class>> A, int nc)
{
    int R = (int)A.size();
    SNF S;
    S.V.assign(nc, vector<mpz_class>(nc));
    S.Vi.assign(nc, vector<mpz_class>(nc));
    for (int i = 0; i < nc; i++) S.V[i][i] = S.Vi[i][i] = 1;
    // column op: col dst += f col src  => V same; Vi: row src -= f row dst
    auto col_add = [&](int dst, int src, const mpz_class &f) {
        for (int i = 0; i < R; i++) A[i][dst] += f * A[i][src];
        for (int i = 0; i < nc; i++) S.V[i][dst] += f * S.V[i][src];
        for (int j = 0; j < nc; j++) S.Vi[src][j] -= f * S.Vi[dst][j];
    };
    auto col_swap = [&](int a, int b) {
        if (a == b) return;
        for (int i = 0; i < R; i++) std::swap(A[i][a], A[i][b]);
        for (int i = 0; i < nc; i++) std::swap(S.V[i][a], S.V[i][b]);
        std::swap(S.Vi[a], S.Vi[b]);
    };
    int t = 0;
    while (t < R && t < nc) {
        // pivot: smallest nonzero |entry| in the remaining block
        int bi = -1, bj = -1;
        for (int i = t; i < R; i++)
            for (int j = t; j < nc; j++)
                if (A[i][j] != 0 && (bi < 0 || abs(A[i][j]) < abs(A[bi][bj]))) { bi = i; bj = j; }
        if (bi < 0) break;
        std::swap(A[t], A[bi]);
        col_swap(t, bj);
        for (;;) {
            bool clean = true;
            for (int j = t + 1; j < nc; j++) {
                if (A[t][j] == 0) continue;
                mpz_class q;
                mpz_fdiv_q(q.get_mpz_t(), A[t][j].get_mpz_t(), A[t][t].get_mpz_t());
                col_add(j, t, -q);
                if (A[t][j] != 0) clean = false;
            }
            for (int i = t + 1; i < R; i++) {
                if (A[i][t] == 0) continue;
                mpz_class q;
                mpz_fdiv_q(q.get_mpz_t(), A[i][t].get_mpz_t(), A[t][t].get_mpz_t());
                for (int j = t; j < nc; j++) A[i][j] -= q * A[t][j];
                if (A[i][t] != 0) clean = false;
            }
            if (clean) {
                int bad = -1;
                for (int i = t + 1; i < R && bad < 0; i++)
                    for (int j = t + 1; j < nc; j++)
                        if (A[i][j] % A[t][t] != 0) { bad = i; break; }
                if (bad < 0) break;
                for (int j = t; j < nc; j++) A[t][j] += A[bad][j];
                continue;
            }
            // bring the smallest nonzero entry of row t / column t to (t,t)
            int mi = t, mj = t;
            for (int j = t; j < nc; j++)
                if (A[t][j] != 0 && abs(A[t][j]) < abs(A[mi][mj])) { mi = t; mj = j; }
            for (int i = t; i < R; i++)
                if (A[i][t] != 0 && abs(A[i][t]) < abs(A[mi][mj])) { mi = i; mj = t; }
            if (mi != t) std::swap(A[t], A[mi]);
            if (mj != t) col_swap(t, mj);
        }
        if (A[t][t] < 0) for (int j = t; j < nc; j++) A[t][j] = -A[t][j];
        S.diag.push_back(A[t][t]);
        t++;
    }
    return S;
}

// ---------------------------------------------------------------- the finite abelian group H

struct HGrp {
    int m = 0;
    vector<int> d;
    double log2h = 0;
    int L = 1;
};

// |ker(A - 1)| = [Z^m : (A-1)Z^m + diag(d)Z^m]
static u128 fixcount(const HGrp &H, const vector<int> &A)
{
    int m = H.m;
    long long L = H.L;
    vector<vector<long long>> cols;
    for (int c = 0; c < m; c++) {
        vector<long long> v(m);
        for (int r = 0; r < m; r++) {
            long long x = A[r * m + c] - (r == c);
            x %= L; if (x < 0) x += L;
            v[r] = x;
        }
        cols.push_back(v);
    }
    for (int c = 0; c < m; c++) {
        vector<long long> v(m, 0);
        v[c] = H.d[c] % L;
        cols.push_back(v);
    }
    vector<char> active(cols.size(), 1);
    u128 index = 1;
    for (int r = 0; r < m; r++) {
        int piv = -1;
        for (size_t c = 0; c < cols.size(); c++) {
            if (!active[c] || cols[c][r] == 0) continue;
            if (piv < 0) { piv = (int)c; continue; }
            long long a = cols[piv][r], b = cols[c][r];
            long long r0 = a, r1 = b, s0 = 1, s1 = 0, t0 = 0, t1 = 1;
            while (r1) {
                long long q = r0 / r1, tmp;
                tmp = r0 - q * r1; r0 = r1; r1 = tmp;
                tmp = s0 - q * s1; s0 = s1; s1 = tmp;
                tmp = t0 - q * t1; t0 = t1; t1 = tmp;
            }
            long long g = r0, u = s0, v = t0, ag = a / g, bg = b / g;
            for (int k = r; k < m; k++) {
                i128 p = cols[piv][k], q = cols[c][k];
                i128 np = (u * p + v * q) % L, nq = (-bg * p + ag * q) % L;
                if (np < 0) np += L;
                if (nq < 0) nq += L;
                cols[piv][k] = (long long)np;
                cols[c][k] = (long long)nq;
            }
        }
        if (piv < 0) { index *= (u128)L; continue; }
        long long a = cols[piv][r];
        long long r0 = a, r1 = L;
        while (r1) { long long q = r0 / r1, tmp = r0 - q * r1; r0 = r1; r1 = tmp; }
        long long g = r0;
        index *= (u128)g;
        long long Lg = L / g;
        vector<long long> Q(m, 0);
        bool nz = false;
        for (int k = r + 1; k < m; k++) {
            i128 q = ((i128)(-Lg) * cols[piv][k]) % L;
            if (q < 0) q += L;
            Q[k] = (long long)q;
            if (q) nz = true;
        }
        active[piv] = 0;
        if (nz) { cols.push_back(Q); active.push_back(1); }
    }
    return index;
}

static bool is_aut(const HGrp &H, const vector<int> &A)
{
    int m = H.m;
    for (int r = 0; r < m; r++)
        for (int c = 0; c < m; c++)
            if (((long long)A[r * m + c] * H.d[c]) % H.d[r] != 0) return false;
    vector<int> B(A);
    for (int i = 0; i < m; i++) B[i * m + i] += 1;
    return fixcount(H, B) == 1;
}

struct Gen {
    vector<int> A;
    int s;
};

static bool burnside(const HGrp &H, const vector<Gen> &G, u128 &oH, u128 &oHS)
{
    int m = H.m, mm = m * m;
    struct Hs {
        size_t operator()(const string &s) const { return std::hash<string>()(s); }
    };
    std::unordered_map<string, int> idx;
    vector<string> elts;
    auto enc = [&](const vector<int> &A, int s) {
        string k(mm * 2 + 1, '\0');
        for (int i = 0; i < mm; i++) { k[2 * i] = (char)(A[i] & 0xff); k[2 * i + 1] = (char)(A[i] >> 8); }
        k[mm * 2] = (char)(s > 0);
        return k;
    };
    auto dec = [&](const string &k, vector<int> &A, int &s) {
        A.resize(mm);
        for (int i = 0; i < mm; i++) A[i] = (unsigned char)k[2 * i] | ((unsigned char)k[2 * i + 1] << 8);
        s = k[mm * 2] ? 1 : -1;
    };
    vector<int> I(mm, 0);
    for (int i = 0; i < m; i++) I[i * m + i] = 1 % H.d[i];
    elts.push_back(enc(I, 1));
    idx[elts[0]] = 0;
    vector<int> E, P(mm);
    int es;
    for (size_t k = 0; k < elts.size(); k++) {
        dec(elts[k], E, es);
        for (auto &g : G) {
            for (int r = 0; r < m; r++)
                for (int c = 0; c < m; c++) {
                    long long s = 0;
                    for (int t = 0; t < m; t++) s += (long long)g.A[r * m + t] * E[t * m + c];
                    P[r * m + c] = (int)(s % H.d[r]);
                }
            string key = enc(P, g.s * es);
            if (!idx.count(key)) {
                long memcap = (768L << 20) / (long)(2 * (2 * mm + 1) + 64);
                if ((long)elts.size() >= std::min(MAX_IMAGE, memcap)) return false;
                idx[key] = (int)elts.size();
                elts.push_back(key);
            }
        }
    }
    u128 sumAll = 0, sumPlus = 0;
    for (auto &k : elts) {
        dec(k, E, es);
        u128 f = fixcount(H, E);
        sumAll += f;
        if (es > 0) sumPlus += 2 * f;
    }
    u128 n = elts.size();
    if (sumAll % n || sumPlus % n) throw std::runtime_error("Burnside sums not divisible");
    oH = sumAll / n;
    oHS = sumPlus / n;
    return true;
}

static void explicit_orbits(const HGrp &H, const vector<Gen> &G, u128 &oH, u128 &oHS)
{
    int m = H.m;
    uint64_t h = 1;
    for (int d : H.d) h *= (uint64_t)d;
    vector<uint64_t> vis((2 * h + 63) / 64, 0);
    vector<uint64_t> radix(m);
    radix[0] = 1;
    for (int i = 1; i < m; i++) radix[i] = radix[i - 1] * H.d[i - 1];
    vector<uint64_t> Q;
    vector<int> x(m);
    oH = oHS = 0;
    for (uint64_t v0 = 0; v0 < h; v0++) {
        for (int pass = 0; pass < 2; pass++) {
            uint64_t st0 = 2 * v0 + pass;
            if ((vis[st0 >> 6] >> (st0 & 63)) & 1) continue;
            if (pass == 0) oH++;
            oHS++;
            Q.clear();
            vis[st0 >> 6] |= 1ULL << (st0 & 63);
            Q.push_back(st0);
            for (size_t qh = 0; qh < Q.size(); qh++) {
                uint64_t st = Q[qh], v = st >> 1;
                int neg = (int)(st & 1);
                for (int i = 0; i < m; i++) { x[i] = (int)(v % H.d[i]); v /= H.d[i]; }
                for (auto &g : G) {
                    uint64_t w = 0;
                    for (int r = 0; r < m; r++) {
                        long long s = 0;
                        for (int c = 0; c < m; c++) s += (long long)g.A[r * m + c] * x[c];
                        w += (uint64_t)(s % H.d[r]) * radix[r];
                    }
                    uint64_t nst = 2 * w + (uint64_t)(neg ^ (g.s < 0));
                    if (!((vis[nst >> 6] >> (nst & 63)) & 1)) {
                        vis[nst >> 6] |= 1ULL << (nst & 63);
                        Q.push_back(nst);
                    }
                }
            }
        }
    }
}

static string u128s(u128 x)
{
    if (!x) return "0";
    string s;
    while (x) { s += char('0' + (int)(x % 10)); x /= 10; }
    std::reverse(s.begin(), s.end());
    return s;
}

// ---------------------------------------------------------------- main analysis

static void analyze(const string &file)
{
    // "A+B": generators from A, normalizer generators from A and from B
    // (B must describe the same group, e.g. the output of CARAT's Normalizer)
    Bravais B;
    size_t plus = file.find('+');
    if (plus != string::npos) {
        B = parse_bravais(file.substr(0, plus));
        Bravais B2 = parse_bravais(file.substr(plus + 1));
        if (B2.gen.size() != B.gen.size()) fail(file, "the two files describe different generators");
        for (size_t i = 0; i < B.gen.size(); i++)
            if (B2.gen[i].a != B.gen[i].a) fail(file, "the two files describe different generators");
        for (auto &x : B2.normal) B.normal.push_back(x);
        for (auto &x : B2.cen) B.cen.push_back(x);
    } else {
        B = parse_bravais(file);
    }
    if (B.gen.empty()) fail(file, "no generators");
    const int n = B.gen[0].r;
    const int k = (int)B.gen.size();
    const int nk = n * k;
    IMat Id;
    Id.r = Id.c = n;
    Id.a.assign(n * n, 0);
    for (int i = 0; i < n; i++) Id.at(i, i) = 1;

    // enumerate G; T[g] is n x nk
    std::unordered_map<string, int> index;
    vector<IMat> elts;
    vector<vector<int>> T;
    elts.push_back(Id);
    index[key_of(Id)] = 0;
    T.push_back(vector<int>((size_t)n * nk, 0));
    RowLattice lat(nk);
    for (size_t q = 0; q < elts.size(); q++) {
        for (int i = 0; i < k; i++) {
            IMat h = mul(elts[q], B.gen[i]);
            string key = key_of(h);
            if (key.empty()) fail(file, "matrix entries too large");
            vector<int> val(T[q]);
            const IMat &g = elts[q];
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++) val[(size_t)a * nk + i * n + b] += (int)g.at(a, b);
            auto it = index.find(key);
            if (it == index.end()) {
                if (elts.size() > 2000000) fail(file, "group too large / infinite");
                index[key] = (int)elts.size();
                elts.push_back(h);
                T.push_back(val);
            } else {
                const vector<int> &Th = T[it->second];
                for (int a = 0; a < n; a++) {
                    vector<i128> row(nk);
                    bool nz = false;
                    for (int c = 0; c < nk; c++) {
                        row[c] = (i128)val[(size_t)a * nk + c] - Th[(size_t)a * nk + c];
                        if (row[c]) nz = true;
                    }
                    if (nz) lat.insert(row, file);
                }
            }
        }
    }
    const size_t order = elts.size();

    // SNF of the constraint lattice basis
    vector<vector<mpz_class>> W;
    for (int c = 0; c < nk; c++) {
        if (lat.piv[c].empty()) continue;
        vector<mpz_class> row(nk);
        for (int j = 0; j < nk; j++) {
            i128 v = lat.piv[c][j];
            bool neg = v < 0;
            u128 u = neg ? (u128)(-v) : (u128)v;
            mpz_class z = (unsigned long)(uint64_t)(u >> 64);
            z <<= 64;
            z += (unsigned long)(uint64_t)u;
            row[j] = neg ? mpz_class(-z) : z;
        }
        W.push_back(row);
    }
    SNF S = snf(W, nk);
    const int r = (int)S.diag.size();
    HGrp H;
    vector<int> hidx;  // which SNF index each H coordinate is
    for (int j = 0; j < r; j++) {
        if (S.diag[j] > 1) {
            if (S.diag[j] > 65535) fail(file, "elementary divisor too large");
            if (order % S.diag[j].get_ui() != 0) fail(file, "elementary divisor does not divide |G|");
            H.d.push_back((int)S.diag[j].get_si());
            hidx.push_back(j);
        }
    }
    H.m = (int)H.d.size();
    // order the cyclic factors so that d[i] | d[i+1] (SNF already does)
    for (int i = 0; i + 1 < H.m; i++)
        if (H.d[i + 1] % H.d[i]) fail(file, "SNF divisibility violated");
    for (int d : H.d) H.log2h += std::log2((double)d);
    H.L = H.m ? H.d[H.m - 1] : 1;

    // action of a matrix N normalizing G
    auto action = [&](const IMat &N, vector<int> &out) {
        IMat Ni = inv_unimodular(N, file);
        // A_N (nk x nk): block row i = N * T_{N^-1 g_i N}
        vector<vector<mpz_class>> A(nk, vector<mpz_class>(nk));
        for (int i = 0; i < k; i++) {
            IMat h = mul(mul(Ni, B.gen[i]), N);
            auto it = index.find(key_of(h));
            if (it == index.end()) fail(file, "normalizer element does not normalize G");
            const vector<int> &Th = T[it->second];
            for (int a = 0; a < n; a++)
                for (int c = 0; c < nk; c++) {
                    long long s = 0;
                    for (int b = 0; b < n; b++) s += N.at(a, b) * Th[(size_t)b * nk + c];
                    A[i * n + a][c] = (long)s;
                }
        }
        // M = Vi * A * V restricted to the needed columns
        int m = H.m;
        out.assign(m * m, 0);
        for (int cc = 0; cc < m; cc++) {
            int j = hidx[cc];
            vector<mpz_class> t(nk), tp(nk), sp(nk);
            for (int a = 0; a < nk; a++) t[a] = S.V[a][j];
            for (int a = 0; a < nk; a++) {
                mpz_class s = 0;
                for (int b = 0; b < nk; b++) if (t[b] != 0) s += A[a][b] * t[b];
                tp[a] = s;
            }
            for (int a = 0; a < r; a++) {
                mpz_class s = 0;
                for (int b = 0; b < nk; b++) if (tp[b] != 0) s += S.Vi[a][b] * tp[b];
                sp[a] = s;
            }
            mpz_class dj = S.diag[j];
            // coordinates with S_aa == 1 must be integral: sp[a] / dj in Z
            for (int a = 0; a < r; a++) {
                if (S.diag[a] == 1 && sp[a] % dj != 0) fail(file, "inconsistent action (trivial coordinate)");
            }
            for (int rr = 0; rr < m; rr++) {
                int a = hidx[rr];
                mpz_class num = sp[a] * S.diag[a];
                if (num % dj != 0) fail(file, "inconsistent action (non-integral)");
                mpz_class v = num / dj;
                mpz_class md = S.diag[a];
                v %= md;
                if (v < 0) v += md;
                out[rr * m + cc] = (int)v.get_si();
            }
        }
    };

    bool improper = false;
    for (auto &g : B.gen) if (det_ll(g) < 0) improper = true;

    vector<Gen> gens;
    if (H.m > 0) {
        // the point group must act trivially
        for (auto &g : B.gen) {
            vector<int> A;
            action(g, A);
            for (int a = 0; a < H.m; a++)
                for (int b = 0; b < H.m; b++)
                    if (A[a * H.m + b] != (a == b ? 1 % H.d[a] : 0)) fail(file, "point group acts nontrivially on H");
        }
        vector<const IMat *> NN;
        for (auto &x : B.cen) NN.push_back(&x);
        for (auto &x : B.normal) NN.push_back(&x);
        for (auto *N : NN) {
            long long dt = det_ll(*N);
            if (dt != 1 && dt != -1) fail(file, "normalizer generator not unimodular");
            vector<int> A;
            action(*N, A);
            if (!is_aut(H, A)) fail(file, "computed action is not an automorphism");
            bool isid = true;
            for (int a = 0; a < H.m; a++)
                for (int b = 0; b < H.m; b++)
                    if (A[a * H.m + b] != (a == b ? 1 % H.d[a] : 0)) isid = false;
            bool dup = isid && dt > 0;
            for (auto &g : gens) if (g.s == (int)dt && g.A == A) dup = true;
            if (!dup) gens.push_back({A, (int)dt});
        }
        if (improper) {
            vector<int> I(H.m * H.m, 0);
            for (int i = 0; i < H.m; i++) I[i * H.m + i] = 1 % H.d[i];
            gens.push_back({I, -1});
        }
    }

    u128 oH, oHS;
    const char *method;
    if (H.m == 0) {
        bool anyneg = improper;
        for (auto &x : B.cen) if (det_ll(x) < 0) anyneg = true;
        for (auto &x : B.normal) if (det_ll(x) < 0) anyneg = true;
        oH = 1;
        oHS = anyneg ? 1 : 2;
        method = "trivial";
    } else if (H.log2h <= 16.0) {
        explicit_orbits(H, gens, oH, oHS);
        method = "explicit";
    } else if (burnside(H, gens, oH, oHS)) {
        method = "burnside";
    } else if (H.log2h <= EXPLICIT_LOG2) {
        explicit_orbits(H, gens, oH, oHS);
        method = "explicit";
    } else {
        fail(file, "H too large and image too large");
    }
    string hs;
    for (int i = 0; i < H.m; i++) hs += (i ? "x" : "") + std::to_string(H.d[i]);
    if (hs.empty()) hs = "1";
    printf("%s\t%zu\t%s\t%s\t%s\t%s\t%d\t%s\n", file.c_str(), order, hs.c_str(), u128s(oH).c_str(),
           u128s(oHS).c_str(), u128s(oHS - oH).c_str(), improper ? 1 : 0, method);
    fflush(stdout);
}

int main(int argc, char **argv)
{
    setvbuf(stdout, nullptr, _IOLBF, 1 << 16);  // one write() per line
    vector<string> files;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a.rfind("-x=", 0) == 0) EXPLICIT_LOG2 = atoi(a.c_str() + 3);
        else if (a.rfind("-y=", 0) == 0) MAX_IMAGE = atol(a.c_str() + 3);
        else files.push_back(a);
    }
    if (files.empty()) {
        fprintf(stderr, "usage: zverify [-x=k] [-y=n] zclass_file...\n");
        return 1;
    }
    int bad = 0;
    for (auto &f : files) {
        try {
            analyze(f);
        } catch (const std::exception &e) {
            printf("%s\tERROR\t%s\n", f.c_str(), e.what());
            fflush(stdout);
            bad++;
        }
    }
    return bad ? 2 : 0;
}
