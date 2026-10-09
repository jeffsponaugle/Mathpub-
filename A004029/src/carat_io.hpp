// carat_io.hpp -- shared helpers for zverify / extend1: reading CARAT
// bravais_TYP files and small exact integer matrix arithmetic.
#pragma once

#include <gmpxx.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using std::string;
using std::vector;

typedef __int128 i128;
typedef unsigned __int128 u128;


#include <stdexcept>

// errors in one input file are reported as exceptions so that a batch run
// can log them and continue with the next file
[[noreturn]] inline void fail(const string &file, const string &msg)
{
    throw std::runtime_error(file + ": " + msg);
}

// ---------------------------------------------------------------- parsing

struct IMat {
    int r = 0, c = 0;
    vector<long long> a;  // row-major
    long long &at(int i, int j) { return a[(size_t)i * c + j]; }
    long long at(int i, int j) const { return a[(size_t)i * c + j]; }
};

struct Bravais {
    vector<IMat> gen, form, zentr, normal, cen;
};

inline Bravais parse_bravais(const string &path)
{
    std::ifstream in(path);
    if (!in) fail(path, "cannot open");
    vector<string> lines;
    string l;
    while (std::getline(in, l)) {
        size_t p = l.find('%');
        if (p != string::npos) l = l.substr(0, p);
        lines.push_back(l);
    }
    size_t pos = 0;
    auto nextline = [&]() -> string {
        while (pos < lines.size() && lines[pos].find_first_not_of(" \t\r") == string::npos) pos++;
        if (pos >= lines.size()) fail(path, "unexpected end of file");
        return lines[pos++];
    };
    int cnt[5] = {0, 0, 0, 0, 0};  // g f z n c
    string first = nextline();
    size_t s0 = first.find_first_not_of(" \t");
    if (first[s0] == '#') {
        const char *keys = "gfznc";
        for (size_t i = s0 + 1; i < first.size(); i++) {
            const char *k = strchr(keys, first[i]);
            if (k && first[i]) cnt[k - keys] = atoi(first.c_str() + i + 1);
        }
    } else {
        cnt[0] = 1;
        pos--;
    }
    // token stream for matrix bodies
    auto read_matrix = [&]() -> IMat {
        string head = nextline();
        size_t a = head.find_first_not_of(" \t");
        head = head.substr(a);
        long long kgv = 1;
        size_t sl = head.find('/');
        if (sl != string::npos) {
            kgv = atoll(head.c_str() + sl + 1);
            head = head.substr(0, sl);
        }
        int rows = atoi(head.c_str());
        int cols = rows;
        bool sym = false, diag = false, scal = false;
        size_t px = head.find('x'), pd = head.find('d');
        if (px != string::npos) {
            int cc = atoi(head.c_str() + px + 1);
            if (cc == 0) sym = true; else cols = cc;
        } else if (pd != string::npos) {
            int cc = atoi(head.c_str() + pd + 1);
            diag = true;
            scal = (cc == 0);
        }
        if (kgv != 1) fail(path, "rational matrix not supported");
        vector<long long> toks;
        size_t need = scal ? 1 : diag ? rows : sym ? (size_t)rows * (rows + 1) / 2 : (size_t)rows * cols;
        while (toks.size() < need) {
            std::istringstream ss(nextline());
            string t;
            while (ss >> t) {
                if (t.find('/') != string::npos) fail(path, "rational entry not supported");
                toks.push_back(atoll(t.c_str()));
            }
        }
        if (toks.size() != need) fail(path, "matrix size mismatch");
        IMat M;
        M.r = rows;
        M.c = cols;
        M.a.assign((size_t)rows * cols, 0);
        size_t t = 0;
        if (scal) {
            for (int i = 0; i < rows; i++) M.at(i, i) = toks[0];
        } else if (diag) {
            for (int i = 0; i < rows; i++) M.at(i, i) = toks[i];
        } else if (sym) {
            for (int i = 0; i < rows; i++)
                for (int j = 0; j <= i; j++) { M.at(i, j) = toks[t]; M.at(j, i) = toks[t]; t++; }
        } else {
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++) M.at(i, j) = toks[t++];
        }
        return M;
    };
    Bravais B;
    vector<IMat> *dst[5] = {&B.gen, &B.form, &B.zentr, &B.normal, &B.cen};
    for (int k = 0; k < 5; k++)
        for (int i = 0; i < cnt[k]; i++) dst[k]->push_back(read_matrix());
    return B;
}

// ---------------------------------------------------------------- small matrix ops

inline IMat mul(const IMat &A, const IMat &B)
{
    IMat C;
    C.r = A.r;
    C.c = B.c;
    C.a.assign((size_t)C.r * C.c, 0);
    for (int i = 0; i < A.r; i++)
        for (int k = 0; k < A.c; k++) {
            long long x = A.at(i, k);
            if (!x) continue;
            for (int j = 0; j < B.c; j++) C.at(i, j) += x * B.at(k, j);
        }
    return C;
}

inline long long det_ll(const IMat &A)
{
    int n = A.r;
    vector<i128> M(A.a.begin(), A.a.end());
    i128 prev = 1;
    int sign = 1;
    for (int k = 0; k < n - 1; k++) {
        if (M[k * n + k] == 0) {
            int p = -1;
            for (int i = k + 1; i < n; i++) if (M[i * n + k] != 0) { p = i; break; }
            if (p < 0) return 0;
            for (int j = 0; j < n; j++) std::swap(M[k * n + j], M[p * n + j]);
            sign = -sign;
        }
        for (int i = k + 1; i < n; i++)
            for (int j = k + 1; j < n; j++)
                M[i * n + j] = (M[i * n + j] * M[k * n + k] - M[i * n + k] * M[k * n + j]) / prev;
        prev = M[k * n + k];
    }
    return sign * (long long)M[(n - 1) * n + n - 1];
}

// inverse of a unimodular matrix (Gauss-Jordan over Q with exact checks)
inline IMat inv_unimodular(const IMat &A, const string &file)
{
    int n = A.r;
    vector<mpq_class> M((size_t)n * 2 * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) M[i * 2 * n + j] = (long)A.at(i, j);
        M[i * 2 * n + n + i] = 1;
    }
    for (int c = 0; c < n; c++) {
        int p = -1;
        for (int r = c; r < n; r++) if (M[r * 2 * n + c] != 0) { p = r; break; }
        if (p < 0) fail(file, "singular matrix");
        for (int j = 0; j < 2 * n; j++) std::swap(M[c * 2 * n + j], M[p * 2 * n + j]);
        mpq_class pv = M[c * 2 * n + c];
        for (int j = 0; j < 2 * n; j++) M[c * 2 * n + j] /= pv;
        for (int r = 0; r < n; r++) {
            if (r == c || M[r * 2 * n + c] == 0) continue;
            mpq_class f = M[r * 2 * n + c];
            for (int j = 0; j < 2 * n; j++) M[r * 2 * n + j] -= f * M[c * 2 * n + j];
        }
    }
    IMat R;
    R.r = R.c = n;
    R.a.resize((size_t)n * n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            mpq_class v = M[i * 2 * n + n + j];
            if (v.get_den() != 1) fail(file, "matrix not unimodular");
            R.at(i, j) = v.get_num().get_si();
        }
    return R;
}

inline string key_of(const IMat &A)
{
    string s(A.a.size() * 2, '\0');
    for (size_t i = 0; i < A.a.size(); i++) {
        long long v = A.a[i];
        if (v < -32768 || v > 32767) return string();
        uint16_t u = (uint16_t)(int16_t)v;
        s[2 * i] = (char)(u & 0xff);
        s[2 * i + 1] = (char)(u >> 8);
    }
    return s;
}

