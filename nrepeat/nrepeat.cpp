// nrepeat.cpp — tool for extending OEIS A331881 / A331882
//
// A331881(n): the first n-digit substring to occur n times in the decimal
//             expansion of the fractional part of Pi.
// A331882(n): the number of digits of the fractional part of Pi needed to
//             contain those n occurrences (i.e. the end position, 1-based,
//             of the n-th occurrence of that substring).
//
// Method:
//   Phase 1 (detection): stream the digit file, slide an n-digit window over
//     every start position, and increment an exact counter tab[value] (direct
//     indexing over all 10^n values -- no hashing, no collisions).  The moment
//     any counter transitions to n, that substring is recorded as a candidate.
//     We finish the current block and stop.
//   Phase 2 (verification): rescan from the beginning, tracking only the
//     candidate substrings, collecting the n smallest occurrence end
//     positions of each.  The winner is the candidate whose n-th occurrence
//     ends first.  (Any substring whose n-th occurrence ends within the
//     scanned prefix necessarily has count >= n there, so it is a candidate;
//     hence the minimum over candidates is exact.)
//
// Counters are 8-bit by default, or 4-bit packed (-c 4) to halve memory
// (valid for n <= 15).  Table memory: 10^n bytes at 8-bit, 10^n/2 at 4-bit.
//
// Build: see Makefile.  Single translation unit, C++17, pthreads.

#include <atomic>
#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <csignal>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

using u8  = uint8_t;
using u64 = uint64_t;

static double wallNow() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

static std::string fmtHMS(double s) {
    if (s < 0 || s > 400 * 86400.0) return "--:--:--";
    long t = (long)s;
    char b[32];
    snprintf(b, sizeof b, "%ld:%02ld:%02ld", t / 3600, (t / 60) % 60, t % 60);
    return b;
}

static std::string fmtCount(double v) {
    char b[32];
    if (v >= 1e12) snprintf(b, sizeof b, "%.3fT", v / 1e12);
    else if (v >= 1e9) snprintf(b, sizeof b, "%.3fG", v / 1e9);
    else if (v >= 1e6) snprintf(b, sizeof b, "%.2fM", v / 1e6);
    else if (v >= 1e3) snprintf(b, sizeof b, "%.1fk", v / 1e3);
    else snprintf(b, sizeof b, "%.0f", v);
    return b;
}

static volatile sig_atomic_t gInterrupted = 0;
static void onSigint(int) { gInterrupted = 1; }

// ---------------------------------------------------------------------------
// Digit reader: streams a text file of pi digits, dropping every non-digit
// byte.  Auto-detects and skips the integer part ("3." or a leading "3"
// followed by "141...") so the stream is the fractional part: 1,4,1,5,9,...
// ---------------------------------------------------------------------------
struct DigitReader {
    int fd = -1;
    u64 fileSize = 0;
    u64 startOff = 0;          // raw byte offset where digit data begins
    std::string prefixNote;    // human description of what was skipped
    u64 maxDigits = UINT64_MAX;

    std::vector<u8> raw;
    size_t rawPos = 0, rawLen = 0;
    u64 rawConsumed = 0;       // raw bytes consumed since reset (for progress)
    u64 digitsOut = 0;
    bool eof = false;

    bool open(const char* path) {
        fd = ::open(path, O_RDONLY);
        if (fd < 0) return false;
        struct stat st{};
        if (fstat(fd, &st) != 0) return false;
        fileSize = (u64)st.st_size;
        raw.resize(1u << 22);  // 4 MB raw buffer

        // Sniff the head of the file to decide how much prefix to skip.
        char head[16] = {0};
        ssize_t got = pread(fd, head, sizeof head, 0);
        (void)got;
        if (head[0] == '3' && head[1] == '.') {
            startOff = 2;
            prefixNote = "leading \"3.\" skipped";
        } else if (memcmp(head, "3141", 4) == 0) {
            startOff = 1;
            prefixNote = "leading \"3\" (no decimal point) skipped";
        } else {
            startOff = 0;
            prefixNote = "no integer-part prefix detected";
        }
        reset();
        return true;
    }

    void reset() {
        lseek(fd, (off_t)startOff, SEEK_SET);
        rawPos = rawLen = 0;
        rawConsumed = startOff;
        digitsOut = 0;
        eof = false;
    }

    // Fill out[0..want) with digit values 0..9; returns number written.
    size_t fill(u8* out, size_t want) {
        if (digitsOut + want > maxDigits)
            want = (maxDigits > digitsOut) ? (size_t)(maxDigits - digitsOut) : 0;
        size_t k = 0;
        while (k < want) {
            if (rawPos == rawLen) {
                if (eof) break;
                ssize_t r = read(fd, raw.data(), raw.size());
                if (r <= 0) { eof = true; break; }
                rawLen = (size_t)r;
                rawPos = 0;
            }
            const u8* p = raw.data() + rawPos;
            size_t   m = rawLen - rawPos;
            size_t   i = 0;
            // filter digits
            for (; i < m && k < want; i++) {
                unsigned c = (unsigned)p[i] - '0';
                if (c <= 9) out[k++] = (u8)c;
            }
            rawPos += i;
            rawConsumed += i;
        }
        digitsOut += k;
        return k;
    }

    ~DigitReader() { if (fd >= 0) close(fd); }
};

// ---------------------------------------------------------------------------
// Options / globals
// ---------------------------------------------------------------------------
struct Options {
    std::string piFile;
    std::string outFile;
    int  n = 0;
    int  threads = 0;
    int  counterBits = 8;
    int  partitions = 1;
    u64  maxBlockDigits = 256ull << 20;   // 256 Mdigit blocks max
    u64  firstBlockDigits = 1ull << 20;   // adaptive: start at 1 Mdigit, double
    u64  maxDigits = UINT64_MAX;
    double progressInterval = 1.0;
    bool quiet = false;
    bool selftest = false;
};

static u64 POW10[20];
static void initPow10() {
    POW10[0] = 1;
    for (int i = 1; i < 20; i++) POW10[i] = POW10[i - 1] * 10;
}

// Shared progress state (written by workers, read by the reporter thread).
struct Progress {
    std::atomic<u64> digitsDone{0};   // window starts processed
    std::atomic<u64> rawBytes{0};     // raw file bytes consumed
    std::atomic<int> candCount{0};
    char phaseBuf[32] = "scan";
    const char* phase = phaseBuf;
    u64 fileSize = 0;
    u64 phaseGoalDigits = 0;          // 0 = unknown (use file size %)
    std::atomic<bool> stopFlag{false};
};

static void progressThread(Progress* pg, double interval) {
    double t0 = wallNow();
    u64 lastDigits = 0;
    double lastT = t0;
    double emaRate = 0;
    while (!pg->stopFlag.load(std::memory_order_relaxed)) {
        std::this_thread::sleep_for(std::chrono::milliseconds((int)(interval * 1000)));
        double t = wallNow();
        u64 d = pg->digitsDone.load(std::memory_order_relaxed);
        u64 rb = pg->rawBytes.load(std::memory_order_relaxed);
        if (d < lastDigits) { lastDigits = d; emaRate = 0; }  // phase change reset
        double inst = (d - lastDigits) / std::max(1e-3, t - lastT);
        emaRate = (emaRate == 0) ? inst : 0.7 * emaRate + 0.3 * inst;
        lastDigits = d; lastT = t;

        double pct, remainDigits;
        if (pg->phaseGoalDigits) {
            pct = 100.0 * d / (double)pg->phaseGoalDigits;
            remainDigits = (double)pg->phaseGoalDigits - (double)d;
        } else {
            pct = pg->fileSize ? 100.0 * rb / (double)pg->fileSize : 0.0;
            // estimate total digits from bytes->digits ratio so far
            double digitsPerByte = rb ? (double)d / (double)rb : 1.0;
            remainDigits = ((double)pg->fileSize - (double)rb) * digitsPerByte;
        }
        double eta = emaRate > 0 ? remainDigits / emaRate : -1;
        fprintf(stderr,
                "\r[%s] %6.2f%% | %sd done | %sd/s | elapsed %s | ETA%s %s | cands %d   ",
                pg->phase, pct,
                fmtCount((double)d).c_str(), fmtCount(emaRate).c_str(),
                fmtHMS(t - t0).c_str(),
                pg->phaseGoalDigits ? "" : "(full file)",
                fmtHMS(eta).c_str(),
                pg->candCount.load(std::memory_order_relaxed));
        fflush(stderr);
    }
    fprintf(stderr, "\n");
}

// ---------------------------------------------------------------------------
// Counter table with atomic saturating increments (8-bit or 4-bit packed).
// ---------------------------------------------------------------------------
struct CounterTable {
    u8* tab = nullptr;
    u64 bytes = 0;
    u64 base = 0;    // value of slot 0 (partition low bound)
    bool nibble = false;
    int target = 0;  // n: report transition to exactly `target`

    std::mutex candMx;
    std::vector<u64> cands;

    bool alloc(u64 slots, u64 lo, int bits) {
        base = lo;
        nibble = (bits == 4);
        bytes = nibble ? (slots + 1) / 2 : slots;
        void* p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                       MAP_PRIVATE | MAP_ANONYMOUS
#ifdef MAP_NORESERVE
                       | MAP_NORESERVE
#endif
                       , -1, 0);
        if (p == MAP_FAILED) return false;
#ifdef MADV_HUGEPAGE
        madvise(p, bytes, MADV_HUGEPAGE);
#endif
        tab = (u8*)p;
        return true;
    }
    void free() { if (tab) { munmap(tab, bytes); tab = nullptr; } }

    inline void recordCandidate(u64 v) {
        std::lock_guard<std::mutex> lk(candMx);
        cands.push_back(v);
    }

    // saturating atomic increment of slot `idx` (value = base + idx);
    // records candidate on transition to target
    inline void inc(u64 idx) {
        if (!nibble) {
            u8* p = tab + idx;
            u8 old = __atomic_load_n(p, __ATOMIC_RELAXED);
            while (old < 255) {
                if (__atomic_compare_exchange_n(p, &old, (u8)(old + 1), true,
                                                __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
                    if ((int)old == target - 1) recordCandidate(base + idx);
                    return;
                }
            }
        } else {
            u8* p = tab + (idx >> 1);
            int sh = (int)(idx & 1) * 4;
            u8 old = __atomic_load_n(p, __ATOMIC_RELAXED);
            for (;;) {
                int c = (old >> sh) & 0xF;
                if (c >= 15) return;
                u8 nw = (u8)(old + (1u << sh));
                if (__atomic_compare_exchange_n(p, &old, nw, true,
                                                __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
                    if (c == target - 1) recordCandidate(base + idx);
                    return;
                }
            }
        }
    }
};

// ---------------------------------------------------------------------------
// Block driver: streams the file in adaptive blocks (1 Md doubling up to the
// configured max), carries the (n-1)-digit tail between blocks, and hands
// each block's window-start range to `threads` workers.
//
//   runRange(dig, lo, hi, g0, tid): process window starts [lo,hi) of buffer
//     `dig`, whose index 0 has global (0-based) start position g0.
//   afterBlock(nextGlobalStart): called after each block with the exclusive
//     bound of processed global starts; return true to stop.
// ---------------------------------------------------------------------------
template <class RangeFn, class AfterFn>
static u64 driveBlocks(DigitReader& rd, int n, const Options& opt, Progress& pg,
                       RangeFn runRange, AfterFn afterBlock) {
    u64 blockSize = std::min(opt.firstBlockDigits, opt.maxBlockDigits);
    std::vector<u8> buf;
    buf.reserve((size_t)(opt.maxBlockDigits + n + 1));
    size_t carryLen = 0;
    u64 g0 = 0;  // global 0-based index of buf[0]
    int T = opt.threads;

    for (;;) {
        buf.resize(carryLen + (size_t)blockSize + 1);
        size_t got = rd.fill(buf.data() + carryLen, (size_t)blockSize);
        size_t L = carryLen + got;
        pg.rawBytes.store(rd.rawConsumed, std::memory_order_relaxed);
        if (L < (size_t)n) return g0;  // no (more) full windows
        size_t W = L - n + 1;
        buf[L] = 0;  // pad: rolling update reads one digit past the last window

        size_t chunk = (W + T - 1) / T;
        std::vector<std::thread> ths;
        for (int t = 0; t < T; t++) {
            size_t lo = (size_t)t * chunk;
            size_t hi = std::min(W, lo + chunk);
            if (lo >= hi) break;
            ths.emplace_back([&, lo, hi, t] { runRange(buf.data(), lo, hi, g0, t); });
        }
        for (auto& th : ths) th.join();

        u64 nextStart = g0 + W;
        bool ended = (got < blockSize);  // EOF or maxDigits cap
        if (afterBlock(nextStart) || ended || gInterrupted) return nextStart;

        size_t keep = (size_t)n - 1;
        memmove(buf.data(), buf.data() + L - keep, keep);
        carryLen = keep;
        g0 += L - keep;
        blockSize = std::min(blockSize * 2, opt.maxBlockDigits);
    }
}

// ---------------------------------------------------------------------------
// Result of one solve
// ---------------------------------------------------------------------------
struct Result {
    bool found = false;
    bool interrupted = false;
    bool partialBest = false;   // interrupted, but a best-so-far exists
    int n = 0;
    int partitions = 1;
    int winPartition = -1;
    u64 value = 0;              // A331881(n) as a number (print zero-padded)
    u64 endPos = 0;             // A331882(n)
    std::vector<u64> occStarts; // 1-based start positions of the n occurrences
    u64 scanDigits = 0;         // window starts processed in detection (all passes)
    u64 verifyDigits = 0;
    size_t candCount = 0;
    double scanSec = 0, verifySec = 0;
    std::string prefixNote;
    u64 fileSize = 0;
};

// merge-keepers for verification: per-candidate count + n smallest end positions
struct VStat {
    u64 count = 0;
    std::vector<u64> ends;  // sorted ascending, size <= n
    void add(u64 e, int n) {
        count++;
        if ((int)ends.size() < n) {
            ends.insert(std::lower_bound(ends.begin(), ends.end(), e), e);
        } else if (e < ends.back()) {
            ends.pop_back();
            ends.insert(std::lower_bound(ends.begin(), ends.end(), e), e);
        }
    }
    void merge(const VStat& o, int n) {
        count += o.count;
        for (u64 e : o.ends) {
            if ((int)ends.size() < n)
                ends.insert(std::lower_bound(ends.begin(), ends.end(), e), e);
            else if (e < ends.back()) {
                ends.pop_back();
                ends.insert(std::lower_bound(ends.begin(), ends.end(), e), e);
            }
        }
    }
};

static Result solve(const Options& opt, int n) {
    Result res;
    res.n = n;
    res.partitions = opt.partitions;

    // open once for metadata + head sanity check
    {
        DigitReader rd0;
        if (!rd0.open(opt.piFile.c_str())) {
            fprintf(stderr, "error: cannot open %s\n", opt.piFile.c_str());
            exit(1);
        }
        res.prefixNote = rd0.prefixNote;
        res.fileSize = rd0.fileSize;
        u8 first20[20];
        size_t g = rd0.fill(first20, 20);
        static const u8 want[20] = {1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4,6};
        if (g < 20 || memcmp(first20, want, 20) != 0)
            fprintf(stderr,
                    "warning: stream does not start with 14159265358979323846 -- "
                    "check the pi file format (%s)\n", rd0.prefixNote.c_str());
    }

    const int K = opt.partitions;
    const u64 slotsTotal = POW10[n];
    const u64 P = POW10[n - 1];

    Progress pg;
    pg.fileSize = res.fileSize;
    std::thread rep;
    if (!opt.quiet) rep = std::thread(progressThread, &pg, opt.progressInterval);

    bool haveBest = false;
    u64 bestVal = 0, bestEnd = 0;
    std::vector<u64> bestStarts;
    int bestPart = -1;

    for (int part = 0; part < K && !gInterrupted; part++) {
        const u64 pLo   = (u64)((__uint128_t)slotsTotal * (unsigned)part / (unsigned)K);
        const u64 pHi   = (u64)((__uint128_t)slotsTotal * (unsigned)(part + 1) / (unsigned)K);
        const u64 pSpan = pHi - pLo;
        if (pSpan == 0) continue;  // more partitions than values (tiny n)

        // A later partition only matters if it completes strictly before the
        // best end found so far, so cap its scan there.
        u64 cap = opt.maxDigits;
        if (haveBest && bestEnd - 1 < cap) cap = bestEnd - 1;
        if (cap < (u64)n) break;  // no window fits under the cap

        CounterTable ct;
        ct.target = n;
        if (!ct.alloc(pSpan, pLo, opt.counterBits)) {
            fprintf(stderr, "error: cannot allocate counter table (%sB). "
                            "Try -c 4, a larger -P, or a smaller n.\n",
                    fmtCount((double)(opt.counterBits == 4 ? (pSpan + 1) / 2 : pSpan)).c_str());
            exit(1);
        }
        if (!opt.quiet && part == 0)
            fprintf(stderr, "n=%d: counter table %sB (%d-bit)%s, threads=%d\n",
                    n, fmtCount((double)ct.bytes).c_str(), opt.counterBits,
                    K > 1 ? " per partition" : "", opt.threads);

        // ---- detection pass for this partition ----
        DigitReader rd;
        rd.open(opt.piFile.c_str());
        rd.maxDigits = cap;
        pg.digitsDone = 0;
        pg.candCount = 0;
        pg.phaseGoalDigits = (cap != UINT64_MAX) ? cap : 0;
        if (K > 1) snprintf(pg.phaseBuf, sizeof pg.phaseBuf, "scan %d/%d", part + 1, K);
        else       snprintf(pg.phaseBuf, sizeof pg.phaseBuf, "scan");

        double t0 = wallNow();
        u64 scanBound = driveBlocks(
            rd, n, opt, pg,
            [&](const u8* dig, size_t lo, size_t hi, u64 /*g0*/, int /*tid*/) {
                u64 v = 0;
                for (int k = 0; k < n; k++) v = v * 10 + dig[lo + k];
                constexpr size_t CH = 4096;
                u64 vals[CH];
                size_t i = lo;
                while (i < hi) {
                    size_t m = std::min(CH, hi - i);
                    size_t m2 = 0;
                    for (size_t j = 0; j < m; j++) {
                        u64 idx = v - pLo;  // unsigned wrap puts v<pLo out of range
                        v = (v - (u64)dig[i + j] * P) * 10 + dig[i + j + n];
                        if (idx < pSpan) {
                            vals[m2++] = idx;
                            __builtin_prefetch(ct.tab + (ct.nibble ? (idx >> 1) : idx), 1, 1);
                        }
                    }
                    for (size_t j = 0; j < m2; j++) ct.inc(vals[j]);
                    pg.digitsDone.fetch_add(m, std::memory_order_relaxed);
                    i += m;
                }
            },
            [&](u64 /*nextStart*/) {
                std::lock_guard<std::mutex> lk(ct.candMx);
                pg.candCount.store((int)ct.cands.size(), std::memory_order_relaxed);
                return !ct.cands.empty();
            });
        res.scanSec += wallNow() - t0;
        res.scanDigits += scanBound;

        std::vector<u64> cands;
        {
            std::lock_guard<std::mutex> lk(ct.candMx);
            cands = ct.cands;
        }
        std::sort(cands.begin(), cands.end());
        cands.erase(std::unique(cands.begin(), cands.end()), cands.end());
        res.candCount += cands.size();
        ct.free();  // release the big table before the verification pass

        if (gInterrupted || cands.empty()) continue;

        // ---- verification pass (exact positions for this partition) ----
        pg.digitsDone = 0;
        pg.phaseGoalDigits = scanBound;
        if (K > 1) snprintf(pg.phaseBuf, sizeof pg.phaseBuf, "verify %d/%d", part + 1, K);
        else       snprintf(pg.phaseBuf, sizeof pg.phaseBuf, "verify");
        t0 = wallNow();

        std::vector<std::vector<VStat>> tstats((size_t)opt.threads);
        for (auto& s : tstats) s.assign(cands.size(), VStat{});
        std::vector<VStat> global(cands.size());
        bool haveWinner = false;
        size_t winIdx = 0;
        u64 winEnd = 0;

        DigitReader rd2;
        rd2.open(opt.piFile.c_str());
        // Need every window with start < scanBound, i.e. digits up to scanBound+n-1.
        rd2.maxDigits = scanBound + (u64)n - 1;

        driveBlocks(
            rd2, n, opt, pg,
            [&](const u8* dig, size_t lo, size_t hi, u64 g0, int tid) {
                u64 v = 0;
                for (int k = 0; k < n; k++) v = v * 10 + dig[lo + k];
                auto& st = tstats[(size_t)tid];
                const u64* cb = cands.data();
                const u64* ce = cb + cands.size();
                for (size_t i = lo; i < hi; i++) {
                    const u64* it = std::lower_bound(cb, ce, v);
                    if (it != ce && *it == v)
                        st[(size_t)(it - cb)].add(g0 + i + (u64)n, n);
                    v = (v - (u64)dig[i] * P) * 10 + dig[i + n];
                }
                pg.digitsDone.fetch_add(hi - lo, std::memory_order_relaxed);
            },
            [&](u64 /*nextStart*/) {
                // fold per-thread stats into the global view, then test for a winner
                for (auto& st : tstats)
                    for (size_t c = 0; c < cands.size(); c++) {
                        if (st[c].count) global[c].merge(st[c], n);
                        st[c] = VStat{};
                    }
                haveWinner = false;
                for (size_t c = 0; c < cands.size(); c++) {
                    if (global[c].count >= (u64)n) {
                        u64 e = global[c].ends[(size_t)n - 1];
                        if (!haveWinner || e < winEnd) { haveWinner = true; winIdx = c; winEnd = e; }
                    }
                }
                return haveWinner;
            });
        res.verifySec += wallNow() - t0;
        res.verifyDigits += pg.digitsDone.load();

        if (haveWinner && (!haveBest || winEnd < bestEnd)) {
            haveBest = true;
            bestVal = cands[winIdx];
            bestEnd = winEnd;
            bestPart = part;
            bestStarts.clear();
            for (int k = 0; k < n; k++)
                bestStarts.push_back(global[winIdx].ends[(size_t)k] - (u64)n + 1);
        }
    }

    pg.stopFlag = true;
    if (rep.joinable()) rep.join();

    if (haveBest) {
        res.value = bestVal;
        res.endPos = bestEnd;
        res.occStarts = bestStarts;
        res.winPartition = bestPart;
    }
    if (gInterrupted) {
        res.interrupted = true;
        res.partialBest = haveBest;
        return res;
    }
    res.found = haveBest;
    return res;
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------
static void reportResult(const Options& opt, const Result& r, FILE* f) {
    fprintf(f, "==== nrepeat report ====\n");
    fprintf(f, "n                    : %d\n", r.n);
    if (r.found) {
        char l1[32], l2[32];
        snprintf(l1, sizeof l1, "A331881(%d)", r.n);
        snprintf(l2, sizeof l2, "A331882(%d)", r.n);
        fprintf(f, "%-21s: %0*" PRIu64 "   (first %d-digit substring occurring %d times)\n",
                l1, r.n, r.value, r.n, r.n);
        fprintf(f, "%-21s: %" PRIu64 "   (digits of frac(Pi) needed)\n", l2, r.endPos);
        fprintf(f, "occurrence starts    :");
        for (u64 s : r.occStarts) fprintf(f, " %" PRIu64, s);
        fprintf(f, "\n");
    } else if (r.interrupted) {
        fprintf(f, "result               : INTERRUPTED after %" PRIu64 " digits\n", r.scanDigits);
        if (r.partialBest) {
            fprintf(f, "best so far (NOT final; partitions incomplete):\n");
            fprintf(f, "  substring %0*" PRIu64 ", end position %" PRIu64 "\n",
                    r.n, r.value, r.endPos);
        }
    } else {
        fprintf(f, "result               : NOT FOUND within %" PRIu64 " digits scanned\n", r.scanDigits);
    }
    fprintf(f, "pi source            : %s (%sB, %s)\n", opt.piFile.c_str(),
            fmtCount((double)r.fileSize).c_str(), r.prefixNote.c_str());
    fprintf(f, "digits scanned       : %" PRIu64 " (detection), %" PRIu64 " (verification)\n",
            r.scanDigits, r.verifyDigits);
    fprintf(f, "candidates           : %zu\n", r.candCount);
    fprintf(f, "threads              : %d\n", opt.threads);
    u64 span = (POW10[r.n] + (u64)r.partitions - 1) / (u64)r.partitions;
    fprintf(f, "counter table        : %s slots, %d-bit (%sB)%s\n",
            fmtCount((double)span).c_str(), opt.counterBits,
            fmtCount((double)(opt.counterBits == 4 ? (span + 1) / 2 : span)).c_str(),
            r.partitions > 1 ? " per partition" : "");
    if (r.partitions > 1) {
        fprintf(f, "partitions           : %d", r.partitions);
        if (r.winPartition >= 0) fprintf(f, " (winner in partition %d)", r.winPartition + 1);
        fprintf(f, "\n");
    }
    double tot = r.scanSec + r.verifySec;
    fprintf(f, "time                 : scan %.1fs, verify %.1fs, total %.1fs (%sd/s avg)\n",
            r.scanSec, r.verifySec, tot,
            fmtCount(tot > 0 ? (r.scanDigits + r.verifyDigits) / tot : 0).c_str());
    fprintf(f, "\n");
}

// ---------------------------------------------------------------------------
// Self-test against the known terms
// ---------------------------------------------------------------------------
static int runSelftest(Options opt) {
    static const u64 kA331881[8] = {1, 26, 446, 2796, 86538, 872117, 1591292, 66416662};
    static const u64 kA331882[8] = {1, 22, 219, 1805, 25499, 168882, 3566679, 29325629};
    int fails = 0;
    for (int n = 1; n <= 8; n++) {
        Result r = solve(opt, n);
        bool ok = r.found && r.value == kA331881[n - 1] && r.endPos == kA331882[n - 1];
        printf("n=%d: %s  A331881=%" PRIu64 " (want %" PRIu64 ")  A331882=%" PRIu64 " (want %" PRIu64 ")\n",
               n, ok ? "PASS" : "FAIL",
               r.found ? r.value : 0, kA331881[n - 1],
               r.found ? r.endPos : 0, kA331882[n - 1]);
        if (!ok) fails++;
        if (gInterrupted) break;
    }
    printf(fails ? "SELFTEST FAILED (%d)\n" : "SELFTEST PASSED\n", fails);
    return fails ? 1 : 0;
}

// ---------------------------------------------------------------------------
static void usage() {
    fprintf(stderr,
        "nrepeat -- extend OEIS A331881 / A331882\n"
        "usage: nrepeat [options] <pi-digit-file>\n"
        "  -n N            compute term N (required unless --selftest)\n"
        "  -o FILE         append the summary report to FILE\n"
        "  -t T            worker threads (default: all hardware threads)\n"
        "  -c 4|8          counter bits (default 8; 4 halves table memory, needs n<=15)\n"
        "  -P K            split the value space into K partitions scanned in K\n"
        "                  passes; divides table memory by K at the cost of K\n"
        "                  reads of the digit prefix (default 1)\n"
        "  -b MDIGITS      max block size in Mdigits (default 256)\n"
        "  -m DIGITS       cap the number of fractional digits scanned\n"
        "  --selftest      recompute n=1..8 and compare with the known terms\n"
        "  --progress SEC  progress interval on stderr (default 1.0)\n"
        "  -q              quiet (no progress output)\n"
        "\nThe pi file is plain text; every non-digit byte is ignored and a\n"
        "leading \"3.\" (or bare \"3\") integer part is detected and skipped.\n"
        "Table memory: 10^n bytes at 8-bit counters (n=9: 1GB, n=10: 10GB,\n"
        "n=11: 100GB, n=12: 1TB), half that with -c 4.\n");
}

int main(int argc, char** argv) {
    initPow10();
    Options opt;

    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto val = [&]() -> const char* {
            if (i + 1 >= argc) { usage(); exit(1); }
            return argv[++i];
        };
        if (a == "-n") opt.n = atoi(val());
        else if (a == "-o") opt.outFile = val();
        else if (a == "-t") opt.threads = atoi(val());
        else if (a == "-c") opt.counterBits = atoi(val());
        else if (a == "-P") opt.partitions = atoi(val());
        else if (a == "-b") opt.maxBlockDigits = strtoull(val(), nullptr, 10) << 20;
        else if (a == "-m") opt.maxDigits = strtoull(val(), nullptr, 10);
        else if (a == "--selftest") opt.selftest = true;
        else if (a == "--progress") opt.progressInterval = atof(val());
        else if (a == "-q") opt.quiet = true;
        else if (a == "-h" || a == "--help") { usage(); return 0; }
        else if (a[0] == '-') { fprintf(stderr, "unknown option %s\n", a.c_str()); usage(); return 1; }
        else opt.piFile = a;
    }

    if (opt.piFile.empty()) { usage(); return 1; }
    if (!opt.selftest && (opt.n < 1 || opt.n > 19)) {
        fprintf(stderr, "error: -n must be 1..19\n");
        return 1;
    }
    if (opt.counterBits != 4 && opt.counterBits != 8) {
        fprintf(stderr, "error: -c must be 4 or 8\n");
        return 1;
    }
    if (opt.counterBits == 4 && opt.n > 15) {
        fprintf(stderr, "error: 4-bit counters support n<=15 only\n");
        return 1;
    }
    if (opt.partitions < 1 || opt.partitions > 100000) {
        fprintf(stderr, "error: -P must be 1..100000\n");
        return 1;
    }
    if (opt.threads <= 0)
        opt.threads = std::max(1u, std::thread::hardware_concurrency());

    signal(SIGINT, onSigint);

    if (opt.selftest) return runSelftest(opt);

    Result r = solve(opt, opt.n);
    reportResult(opt, r, stdout);
    if (!opt.outFile.empty()) {
        FILE* f = fopen(opt.outFile.c_str(), "a");
        if (!f) {
            fprintf(stderr, "error: cannot open output file %s\n", opt.outFile.c_str());
            return 1;
        }
        time_t t = time(nullptr);
        char ts[64];
        strftime(ts, sizeof ts, "%Y-%m-%d %H:%M:%S", localtime(&t));
        fprintf(f, "# run at %s\n", ts);
        reportResult(opt, r, f);
        fclose(f);
    }
    if (r.interrupted) return 130;
    return r.found ? 0 : 2;
}
