// extend1.cpp -- build (n+1)-dimensional finite unimodular groups that have a
// one-dimensional rational constituent from n-dimensional Q-class
// representatives, and fingerprint them.
//
// For a finite G <= GL_n(Z) with generators g_1..g_k it writes
//   <id>__T       generators g_i (+) 1                    (trivial on e_{n+1})
//   <id>__C       generators g_i (+) 1  and  I (+) -1     (G x C_2)
//   <id>__X<bits> generators g_i (+) chi(g_i)             (each nontrivial
//                                                          chi: G -> {+-1})
// Every finite subgroup of GL_{n+1}(Q) with a one-dimensional rational
// constituent is Q-conjugate to one of these, for a suitable Q-class
// representative G: the constituent is a character chi of order <= 2, and
// the group is either the graph of chi on its image in GL_n(Q), or that
// image times C_2.
//
// Fingerprint, an invariant of the Q-class (= isomorphism class of the group
// together with its rational character): |G|; the sorted multiset, over
// conjugacy classes, of (class size, element order, tr(g), ..., tr(g^dim)),
// each extended by the same data for the classes of g^2, g^3, g^5, g^7
// (power maps); and the sorted multiset, over the sign characters chi, of
// (multiplicity of chi in the representation, profile of ker chi).
// Q-equivalent groups always have equal fingerprints, so keeping one group
// per fingerprint never counts a Q-class twice (it may merge distinct ones,
// which only makes a lower bound weaker).
//
// usage: extend1 -o outdir file...    construct and fingerprint the extensions
//        extend1 -f file...           fingerprint the given groups only
//        extend1 -z file...           Z-class fingerprint (adds GL_n(Z)-conjugacy
//                                     invariants: elementary divisors of g-1, g+1,
//                                     the norm of <g> per class, coinvariants of
//                                     L and L^*); equal for Z-equivalent groups
// output: name \t |G| \t fingerprint-hash \t has1d \t #classes
//         (has1d = the group has a 1-dim rational constituent)

#include "carat_io.hpp"

#include <unordered_map>

using std::unordered_map;

static string factor_string(unsigned long long n)
{
    string s;
    unsigned long long m = n;
    for (unsigned long long p = 2; p * p <= m; p++) {
        int e = 0;
        while (m % p == 0) { m /= p; e++; }
        if (e) s += (s.empty() ? "" : " * ") + std::to_string(p) + "^" + std::to_string(e);
    }
    if (m > 1) s += (s.empty() ? "" : " * ") + std::to_string(m) + "^1";
    if (s.empty()) s = "1";
    return s + " = " + std::to_string(n);
}

struct Group {
    int n = 0;
    vector<IMat> gens;
    vector<IMat> elts;
    unordered_map<string, int> index;
    vector<int> parent, pgen;  // BFS tree (right multiplication by generators)
    vector<int> next;          // next[e*k + i] = index of elts[e] * gens[i]
};

static void enumerate(Group &G, const string &file)
{
    IMat I;
    I.r = I.c = G.n;
    I.a.assign((size_t)G.n * G.n, 0);
    for (int i = 0; i < G.n; i++) I.at(i, i) = 1;
    G.elts.assign(1, I);
    G.index.clear();
    G.index[key_of(I)] = 0;
    G.parent.assign(1, -1);
    G.pgen.assign(1, -1);
    G.next.clear();
    for (size_t q = 0; q < G.elts.size(); q++) {
        for (size_t i = 0; i < G.gens.size(); i++) {
            IMat h = mul(G.elts[q], G.gens[i]);
            string k = key_of(h);
            if (k.empty()) fail(file, "entries too large");
            auto it = G.index.find(k);
            if (it == G.index.end()) {
                if (G.elts.size() > 5000000) fail(file, "group too large");
                G.next.push_back((int)G.elts.size());
                G.index[k] = (int)G.elts.size();
                G.elts.push_back(h);
                G.parent.push_back((int)q);
                G.pgen.push_back((int)i);
            } else {
                G.next.push_back(it->second);
            }
        }
    }
}

// all homomorphisms G -> {+-1}, as sign vectors on the generators
static vector<vector<int>> sign_characters(const Group &G)
{
    vector<vector<int>> res;
    size_t k = G.gens.size();
    if (k > 20) return res;
    vector<int> chi(G.elts.size());
    for (unsigned long mask = 0; mask < (1UL << k); mask++) {
        vector<int> s(k);
        for (size_t i = 0; i < k; i++) s[i] = (mask >> i & 1) ? -1 : 1;
        chi[0] = 1;
        for (size_t e = 1; e < G.elts.size(); e++) chi[e] = chi[G.parent[e]] * s[G.pgen[e]];
        bool ok = true;
        for (size_t e = 0; e < G.elts.size() && ok; e++)
            for (size_t i = 0; i < k && ok; i++)
                if (chi[G.next[e * k + i]] != chi[e] * s[i]) ok = false;
        if (ok) res.push_back(s);
    }
    return res;
}

static long long trace(const IMat &A)
{
    long long t = 0;
    for (int i = 0; i < A.r; i++) t += A.at(i, i);
    return t;
}

struct FP {
    string hash;
    int nclasses;
    bool has1d;
};

static uint64_t fnv(const string &s, uint64_t h)
{
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
    return h;
}

// elementary divisors of a small integer matrix, as "d1.d2...dr/rank"
static string snf_string(vector<vector<long long>> A, const string &file)
{
    const long long LIM = 1LL << 61;
    int R = (int)A.size(), C = R ? (int)A[0].size() : 0;
    vector<long long> d;
    int t = 0;
    while (t < R && t < C) {
        int bi = -1, bj = -1;
        for (int i = t; i < R; i++)
            for (int j = t; j < C; j++)
                if (A[i][j] != 0 && (bi < 0 || std::llabs(A[i][j]) < std::llabs(A[bi][bj]))) { bi = i; bj = j; }
        if (bi < 0) break;
        std::swap(A[t], A[bi]);
        for (int i = 0; i < R; i++) std::swap(A[i][t], A[i][bj]);
        for (;;) {
            bool clean = true;
            long long p = A[t][t];
            for (int i = t + 1; i < R; i++) {
                if (!A[i][t]) continue;
                long long q = A[i][t] / p;
                for (int j = t; j < C; j++) {
                    __int128 v = (__int128)A[i][j] - (__int128)q * A[t][j];
                    if (v > LIM || v < -LIM) fail(file, "SNF overflow");
                    A[i][j] = (long long)v;
                }
                if (A[i][t]) clean = false;
            }
            for (int j = t + 1; j < C; j++) {
                if (!A[t][j]) continue;
                long long q = A[t][j] / p;
                for (int i = t; i < R; i++) {
                    __int128 v = (__int128)A[i][j] - (__int128)q * A[i][t];
                    if (v > LIM || v < -LIM) fail(file, "SNF overflow");
                    A[i][j] = (long long)v;
                }
                if (A[t][j]) clean = false;
            }
            if (!clean) {
                // move the smallest nonzero entry of row/column t to the pivot
                int mi = t, mj = t;
                for (int i = t; i < R; i++)
                    if (A[i][t] && std::llabs(A[i][t]) < std::llabs(A[mi][mj])) { mi = i; mj = t; }
                for (int j = t; j < C; j++)
                    if (A[t][j] && std::llabs(A[t][j]) < std::llabs(A[mi][mj])) { mi = t; mj = j; }
                std::swap(A[t], A[mi]);
                for (int i = 0; i < R; i++) std::swap(A[i][t], A[i][mj]);
                continue;
            }
            int bad = -1;
            for (int i = t + 1; i < R && bad < 0; i++)
                for (int j = t + 1; j < C; j++)
                    if (A[i][j] % p) { bad = i; break; }
            if (bad < 0) break;
            for (int j = t; j < C; j++) A[t][j] += A[bad][j];
        }
        d.push_back(std::llabs(A[t][t]));
        t++;
    }
    string s;
    for (long long x : d) s += std::to_string(x) + ".";
    return s + "/" + std::to_string(d.size());
}

static vector<vector<long long>> rows_of(const IMat &A, long long diag)
{
    vector<vector<long long>> M(A.r, vector<long long>(A.c));
    for (int i = 0; i < A.r; i++)
        for (int j = 0; j < A.c; j++) M[i][j] = A.at(i, j) + (i == j ? diag : 0);
    return M;
}

static FP fingerprint(Group &G, const string &file, bool zlevel = false)
{
    size_t N = G.elts.size();
    int n = G.n;
    vector<IMat> ginv;
    for (auto &g : G.gens) ginv.push_back(inv_unimodular(g, file));
    const string idkey = key_of(G.elts[0]);
    // conjugacy classes
    vector<int> cls(N, -1);
    vector<int> rep;
    vector<size_t> csize;
    for (size_t e = 0; e < N; e++) {
        if (cls[e] >= 0) continue;
        int c = (int)rep.size();
        vector<int> orbit(1, (int)e);
        cls[e] = c;
        for (size_t q = 0; q < orbit.size(); q++)
            for (size_t i = 0; i < G.gens.size(); i++) {
                IMat h = mul(mul(G.gens[i], G.elts[orbit[q]]), ginv[i]);
                int j = G.index.at(key_of(h));
                if (cls[j] < 0) { cls[j] = c; orbit.push_back(j); }
            }
        rep.push_back((int)e);
        csize.push_back(orbit.size());
    }
    const int nc = (int)rep.size();
    // level 1: class size, element order, traces of x^1..x^n
    vector<string> base(nc);
    vector<vector<int>> powcls(nc, vector<int>(4, 0));  // classes of x^2, x^3, x^5, x^7
    const int primes[4] = {2, 3, 5, 7};
    for (int c = 0; c < nc; c++) {
        const IMat &x = G.elts[rep[c]];
        string inv = std::to_string(csize[c]);
        vector<IMat> pw(1, x);  // pw[j] = x^(j+1)
        for (int j = 1; j < std::max(n, 7); j++) pw.push_back(mul(pw.back(), x));
        int ord = 1;
        IMat q = x;
        while (key_of(q) != idkey) { q = mul(q, x); ord++; }
        inv += "," + std::to_string(ord);
        for (int j = 0; j < n; j++) inv += "," + std::to_string(trace(pw[j]));
        base[c] = inv;
        for (int t = 0; t < 4; t++) powcls[c][t] = cls[G.index.at(key_of(pw[primes[t] - 1]))];
    }
    // level 2: add the classes of the prime powers (power maps)
    vector<string> pm(nc);
    for (int c = 0; c < nc; c++) {
        pm[c] = base[c];
        for (int t = 0; t < 4; t++) pm[c] += "|" + base[powcls[c][t]];
    }
    // Z-level (GL_n(Z)-conjugacy) invariants of each class: elementary
    // divisors of g - 1, g + 1 and of the norm 1 + g + ... + g^(ord-1)
    if (zlevel) {
        for (int c = 0; c < nc; c++) {
            const IMat &x = G.elts[rep[c]];
            IMat N = G.elts[0], p = G.elts[0];
            for (;;) {
                p = mul(p, x);
                if (key_of(p) == idkey) break;
                for (size_t i = 0; i < N.a.size(); i++) N.a[i] += p.a[i];
            }
            pm[c] += "|" + snf_string(rows_of(x, -1), file) + "|" + snf_string(rows_of(x, 1), file) + "|" +
                     snf_string(rows_of(N, 0), file);
        }
    }
    vector<string> parts(pm);
    std::sort(parts.begin(), parts.end());
    string canon = std::to_string(N) + ":";
    for (auto &x : parts) canon += x + ";";
    // sign characters: multiplicity in the representation and kernel profile
    FP f;
    f.has1d = false;
    vector<string> chis;
    vector<int> chi(N);
    for (auto &sv : sign_characters(G)) {
        chi[0] = 1;
        for (size_t e = 1; e < N; e++) chi[e] = chi[G.parent[e]] * sv[G.pgen[e]];
        long long sum = 0;
        for (size_t e = 0; e < N; e++) sum += chi[e] * trace(G.elts[e]);
        if (sum % (long long)N) fail(file, "non-integral character multiplicity");
        long long mult = sum / (long long)N;
        if (mult > 0) f.has1d = true;
        vector<string> ker;
        for (int c = 0; c < nc; c++) if (chi[rep[c]] == 1) ker.push_back(pm[c]);
        std::sort(ker.begin(), ker.end());
        string prof = std::to_string(mult) + "{";
        for (auto &x : ker) prof += x + ";";
        chis.push_back(prof + "}");
    }
    std::sort(chis.begin(), chis.end());
    canon += "#";
    for (auto &x : chis) canon += x + "/";
    if (zlevel) {
        // coinvariants of L and of its dual: Z^n / sum (g_i - 1) Z^n
        vector<vector<long long>> M(n), Mt(n);
        for (auto &g : G.gens)
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    M[i].push_back(g.at(i, j) - (i == j));
                    Mt[i].push_back(g.at(j, i) - (i == j));
                }
        canon += "#" + snf_string(M, file) + "#" + snf_string(Mt, file);
    }
    char buf[64];
    snprintf(buf, sizeof buf, "%016llx%016llx", (unsigned long long)fnv(canon, 1469598103934665603ULL),
             (unsigned long long)fnv(canon, 0x9e3779b97f4a7c15ULL));
    f.hash = buf;
    f.nclasses = nc;
    return f;
}

static void write_group(const string &path, const vector<IMat> &gens, size_t order)
{
    FILE *f = fopen(path.c_str(), "w");
    if (!f) fail(path, "cannot write");
    fprintf(f, "#g%zu\n", gens.size());
    for (auto &g : gens) {
        fprintf(f, "%d\n", g.r);
        for (int i = 0; i < g.r; i++) {
            for (int j = 0; j < g.c; j++) fprintf(f, "%lld ", g.at(i, j));
            fprintf(f, "\n");
        }
    }
    fprintf(f, "%s\n", factor_string(order).c_str());
    fclose(f);
}

static IMat extend(const IMat &g, int last)
{
    IMat h;
    h.r = h.c = g.r + 1;
    h.a.assign((size_t)h.r * h.c, 0);
    for (int i = 0; i < g.r; i++)
        for (int j = 0; j < g.c; j++) h.at(i, j) = g.at(i, j);
    h.at(g.r, g.r) = last;
    return h;
}

static string ident(const string &path)
{
    // .../dir.<family>/ordnung.<o>/<type>/<name>  ->  <family>__<name>
    string name = path.substr(path.find_last_of('/') + 1);
    string fam = "x";
    size_t p = path.find("/dir.");
    if (p != string::npos) {
        size_t e = path.find('/', p + 5);
        fam = path.substr(p + 5, e - p - 5);
        for (auto &ch : fam) ch = ch == ';' ? 's' : ch == ',' ? 'c' : ch == '\'' ? 'p' : ch;
    }
    return fam + "__" + name;
}

int main(int argc, char **argv)
{
    // one write() per line: several processes may share this stdout
    setvbuf(stdout, nullptr, _IOLBF, 1 << 16);
    string outdir;
    bool fponly = false, zfp = false;
    vector<string> files;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "-o" && i + 1 < argc) outdir = argv[++i];
        else if (a == "-f") fponly = true;
        else if (a == "-z") { fponly = true; zfp = true; }
        else files.push_back(a);
    }
    if (files.empty() || (!fponly && outdir.empty())) {
        fprintf(stderr, "usage: extend1 -o outdir file... | extend1 -f file...\n");
        return 1;
    }
    for (auto &file : files) try {
        Bravais B = parse_bravais(file);
        Group G;
        G.n = B.gen[0].r;
        G.gens = B.gen;
        enumerate(G, file);
        if (fponly) {
            FP f = fingerprint(G, file, zfp);
            printf("%s\t%zu\t%s\t%d\t%d\n", file.c_str(), G.elts.size(), f.hash.c_str(), f.has1d, f.nclasses);
            continue;
        }
        string id = ident(file);
        vector<std::pair<string, vector<IMat>>> cons;
        {
            vector<IMat> g;
            for (auto &x : G.gens) g.push_back(extend(x, 1));
            cons.push_back({"T", g});
            IMat e = extend(G.elts[0], -1);
            g.push_back(e);
            cons.push_back({"C", g});
        }
        for (auto &s : sign_characters(G)) {
            bool trivial = true;
            string bits;
            for (int v : s) { if (v < 0) trivial = false; bits += v < 0 ? '1' : '0'; }
            if (trivial) continue;
            vector<IMat> g;
            for (size_t i = 0; i < G.gens.size(); i++) g.push_back(extend(G.gens[i], s[i]));
            cons.push_back({"X" + bits, g});
        }
        for (auto &c : cons) {
            Group H;
            H.n = G.n + 1;
            H.gens = c.second;
            enumerate(H, file);
            FP f = fingerprint(H, file);
            string out = outdir + "/" + id + "__" + c.first;
            write_group(out, H.gens, H.elts.size());
            printf("%s\t%zu\t%s\t%d\t%d\n", out.c_str(), H.elts.size(), f.hash.c_str(), f.has1d, f.nclasses);
        }
        fflush(stdout);
    } catch (const std::exception &e) {
        fprintf(stderr, "extend1: %s\n", e.what());
    }
    return 0;
}
