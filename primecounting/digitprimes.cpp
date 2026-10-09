// digitprimes: extend OEIS A231412, A228413-A228421 (first 10^n primes with no
// digit d) and A231726/A231787-A231790/A231792-A231796 (at least one digit d).
//
// Iterates all primes in order with primesieve, computes a 10-bit
// digit-presence mask per prime, and records the per-digit counts every time
// the running prime count crosses a power of 10. One pass produces all 20
// sequences at once. Multithreaded over fixed-size ranges, merged in order so
// crossings are exact; checkpoints so a multi-hour run can be resumed.
//
// Build:  make        (needs primesieve, e.g. `brew install primesieve`)
// Run:    ./digitprimes --target 13
//
// The next unknown term of the family is a(13): digit stats over the first
// 10^13 primes (all primes below ~3.24e14).

#include <primesieve.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using u64 = uint64_t;

// ---------------------------------------------------------------------------
// Digit-presence mask: bit d set iff decimal digit d appears in n.
// Table lookup 4 digits at a time; PAD includes leading zeros (for interior
// digit groups), TOP does not (for the most significant group).
// ---------------------------------------------------------------------------
static uint16_t PAD[10000], TOP[10000];

static void build_tables() {
    for (int v = 0; v < 10000; v++) {
        int a = v / 1000, b = v / 100 % 10, c = v / 10 % 10, d = v % 10;
        PAD[v] = (uint16_t)((1 << a) | (1 << b) | (1 << c) | (1 << d));
        uint16_t m = (uint16_t)(1 << d);
        if (v >= 10) m |= (uint16_t)(1 << c);
        if (v >= 100) m |= (uint16_t)(1 << b);
        if (v >= 1000) m |= (uint16_t)(1 << a);
        TOP[v] = m;
    }
}

static inline uint32_t digit_mask(u64 n) {
    uint32_t m = 0;
    while (n >= 10000) {
        m |= PAD[n % 10000];
        n /= 10000;
    }
    return m | TOP[n];
}

// ---------------------------------------------------------------------------
// Known OEIS values for verification: missing[d][n] = count of the first 10^n
// primes that do not contain digit d, for n = 0..12.
// ---------------------------------------------------------------------------
static const int KNOWN_N = 13;  // n = 0..12
static const u64 KNOWN[10][KNOWN_N] = {
    /*0 A231412*/ {1, 10, 91, 819, 7122, 61702, 557224, 5062320, 45002763, 395879190, 3579400605, 32487367715, 294505958253},
    /*1 A228413*/ {1, 6, 54, 532, 4675, 34425, 262549, 2051466, 16831152, 155616459, 1529462564, 14830618421, 141585123501},
    /*2 A228414*/ {0, 7, 77, 697, 6497, 55552, 512100, 4710641, 42205969, 341224891, 2787791578, 22971326749, 190650687957},
    /*3 A228415*/ {1, 7, 54, 534, 4909, 45405, 385008, 3539880, 32260781, 294001190, 2564080248, 23271246324, 211753431947},
    /*4 A228416*/ {1, 10, 75, 721, 6637, 60605, 514809, 4730382, 43254591, 392344689, 3421561753, 31049600245, 282499317912},
    /*5 A228417*/ {1, 9, 85, 708, 6635, 60640, 535534, 4737129, 43297195, 392641522, 3536880527, 31067514571, 282635824867},
    /*6 A228418*/ {1, 10, 90, 719, 6696, 60845, 554933, 4742037, 43331008, 392875212, 3573268469, 31207451849, 282765603085},
    /*7 A228419*/ {1, 8, 67, 539, 5034, 45549, 416913, 3570781, 32517377, 294828478, 2681147149, 23720397369, 212156228217},
    /*8 A228420*/ {1, 10, 92, 816, 6712, 60867, 555878, 5026796, 43410238, 395243878, 3576361255, 32461990759, 282971130960},
    /*9 A228421*/ {1, 8, 69, 620, 5010, 45732, 418142, 3785060, 32579606, 296601070, 2683254222, 24354108057, 212324183352},
};

static const char* OEIS_MISS[10] = {"A231412", "A228413", "A228414", "A228415", "A228416",
                                    "A228417", "A228418", "A228419", "A228420", "A228421"};
static const char* OEIS_HAS[10] = {"A231726", "A231787", "A231788", "A231789", "A231790",
                                   "A231792", "A231793", "A231794", "A231795", "A231796"};

// ---------------------------------------------------------------------------
// Shared state
// ---------------------------------------------------------------------------
struct ChunkResult {
    u64 nprimes = 0;
    u64 missing[10] = {0};
};

struct Row {           // counts at the moment the 10^n-th prime is reached
    int n;             // power of 10
    u64 prime;         // the 10^n-th prime itself
    u64 missing[10];   // cumulative "no digit d" counts
};

static u64 g_chunk = 1000000000ULL;  // numbers per chunk
static int g_target = 13;            // compute through the first 10^target primes
static std::string g_ckpt = "digitprimes.ckpt";
static std::string g_out = "digitprimes_results.txt";

static std::mutex g_mx;
static std::condition_variable g_cv_main, g_cv_space;
static std::map<u64, ChunkResult> g_results;
static std::atomic<u64> g_next_claim{0};
static std::atomic<bool> g_stop{false};
static u64 g_merge_next = 0;      // next chunk index to merge (protected by g_mx for claims)
static volatile sig_atomic_t g_sigint = 0;

static void on_sigint(int) { g_sigint = 1; }

static u64 pow10u(int n) {
    u64 r = 1;
    while (n--) r *= 10;
    return r;
}

// ---------------------------------------------------------------------------
// Worker: claim chunk indices, sieve (lo, hi], histogram digit masks.
// ---------------------------------------------------------------------------
static void worker(int max_ahead) {
    std::vector<u64> hist(1024);
    while (true) {
        u64 idx;
        {
            std::unique_lock<std::mutex> lk(g_mx);
            g_cv_space.wait(lk, [&] {
                return g_stop.load() || g_next_claim.load() < g_merge_next + (u64)max_ahead;
            });
            if (g_stop.load()) return;
            idx = g_next_claim.fetch_add(1);
        }
        u64 lo = idx * g_chunk;         // exclusive
        u64 hi = lo + g_chunk;          // inclusive
        std::fill(hist.begin(), hist.end(), 0);
        primesieve::iterator it(lo, hi);
        u64 p, np = 0;
        while ((p = it.next_prime()) <= hi) {
            hist[digit_mask(p)]++;
            np++;
        }
        ChunkResult r;
        r.nprimes = np;
        for (int m = 0; m < 1024; m++) {
            u64 c = hist[m];
            if (!c) continue;
            for (int d = 0; d < 10; d++)
                if (!(m >> d & 1)) r.missing[d] += c;
        }
        {
            std::lock_guard<std::mutex> lk(g_mx);
            g_results[idx] = r;
        }
        g_cv_main.notify_one();
    }
}

// ---------------------------------------------------------------------------
// Serial re-scan of one chunk to pin down exact counts at each 10^n crossing
// that falls inside it. Called rarely (once per crossing chunk).
// Returns the full-chunk totals so the caller can cross-check the parallel
// result. cum/cum_missing are the totals for everything before this chunk.
// ---------------------------------------------------------------------------
static ChunkResult rescan_chunk(u64 idx, u64 cum, const u64 cum_missing[10],
                                int& next_target, std::vector<Row>& rows) {
    u64 lo = idx * g_chunk, hi = lo + g_chunk;
    ChunkResult r;
    primesieve::iterator it(lo, hi);
    u64 p;
    while ((p = it.next_prime()) <= hi) {
        uint32_t m = digit_mask(p);
        r.nprimes++;
        for (int d = 0; d < 10; d++)
            if (!(m >> d & 1)) r.missing[d]++;
        if (next_target <= g_target && cum + r.nprimes == pow10u(next_target)) {
            Row row;
            row.n = next_target;
            row.prime = p;
            for (int d = 0; d < 10; d++) row.missing[d] = cum_missing[d] + r.missing[d];
            rows.push_back(row);
            printf("\n[crossing] pi = 10^%d  (p = %llu)\n", row.n, (unsigned long long)p);
            for (int d = 0; d < 10; d++) {
                const char* mark = "  (new)";
                if (row.n < KNOWN_N)
                    mark = (row.missing[d] == KNOWN[d][row.n]) ? "  OK" : "  *** MISMATCH ***";
                printf("  no %d (%s): %-15llu%s\n", d, OEIS_MISS[d],
                       (unsigned long long)row.missing[d], mark);
            }
            fflush(stdout);
            next_target++;
        }
    }
    return r;
}

// ---------------------------------------------------------------------------
// Checkpoint / results I/O
// ---------------------------------------------------------------------------
static void save_checkpoint(u64 merge_next, u64 cum, const u64 cum_missing[10],
                            int next_target, const std::vector<Row>& rows) {
    std::string tmp = g_ckpt + ".tmp";
    std::ofstream f(tmp);
    f << "digitprimes-ckpt-v1\n";
    f << "chunk " << g_chunk << "\n";
    f << "target " << g_target << "\n";
    f << "merge_next " << merge_next << "\n";
    f << "cum " << cum << "\n";
    f << "cum_missing";
    for (int d = 0; d < 10; d++) f << " " << cum_missing[d];
    f << "\nnext_target " << next_target << "\n";
    f << "rows " << rows.size() << "\n";
    for (auto& r : rows) {
        f << r.n << " " << r.prime;
        for (int d = 0; d < 10; d++) f << " " << r.missing[d];
        f << "\n";
    }
    f.close();
    std::rename(tmp.c_str(), g_ckpt.c_str());
}

static bool load_checkpoint(u64& merge_next, u64& cum, u64 cum_missing[10],
                            int& next_target, std::vector<Row>& rows) {
    std::ifstream f(g_ckpt);
    if (!f) return false;
    std::string tag;
    u64 chunk;
    int target;
    f >> tag;
    if (tag != "digitprimes-ckpt-v1") {
        fprintf(stderr, "checkpoint %s: unknown format, ignoring\n", g_ckpt.c_str());
        return false;
    }
    f >> tag >> chunk >> tag >> target >> tag >> merge_next >> tag >> cum >> tag;
    for (int d = 0; d < 10; d++) f >> cum_missing[d];
    f >> tag >> next_target >> tag;
    size_t nrows;
    f >> nrows;
    rows.resize(nrows);
    for (auto& r : rows) {
        f >> r.n >> r.prime;
        for (int d = 0; d < 10; d++) f >> r.missing[d];
    }
    if (!f) {
        fprintf(stderr, "checkpoint %s: parse error, ignoring\n", g_ckpt.c_str());
        rows.clear();
        return false;
    }
    if (chunk != g_chunk) {
        fprintf(stderr, "checkpoint chunk size %llu != --chunk %llu; use matching --chunk or --fresh\n",
                (unsigned long long)chunk, (unsigned long long)g_chunk);
        exit(1);
    }
    if (target > g_target) g_target = target;
    return true;
}

static void write_results(const std::vector<Row>& rows) {
    std::ofstream f(g_out);
    f << "# Count of the first 10^n primes by decimal digit content.\n";
    f << "# Generated by digitprimes (single ordered pass over all primes with primesieve).\n#\n";
    f << "# n, 10^n-th prime, then for d = 0..9 the count with no digit d.\n";
    f << "n,nth_prime";
    for (int d = 0; d < 10; d++) f << ",no_" << d;
    f << "\n";
    for (auto& r : rows) {
        f << r.n << "," << r.prime;
        for (int d = 0; d < 10; d++) f << "," << r.missing[d];
        f << "\n";
    }
    f << "\n# OEIS sequences, 'no digit d' (a(0)..a(" << (rows.empty() ? 0 : rows.back().n) << ")):\n";
    for (int d = 0; d < 10; d++) {
        f << "# " << OEIS_MISS[d] << " (no " << d << "):";
        for (auto& r : rows) f << " " << r.missing[d] << ",";
        f << "\n";
    }
    f << "\n# OEIS sequences, 'at least one digit d' = 10^n - no_d (a(1)..):\n";
    for (int d = 0; d < 10; d++) {
        f << "# " << OEIS_HAS[d] << " (>=1 " << d << "):";
        for (auto& r : rows)
            if (r.n >= 1) f << " " << (pow10u(r.n) - r.missing[d]) << ",";
        f << "\n";
    }
}

// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    int threads = (int)std::thread::hardware_concurrency();
    bool fresh = false;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() -> const char* {
            if (i + 1 >= argc) { fprintf(stderr, "%s needs a value\n", a.c_str()); exit(1); }
            return argv[++i];
        };
        if (a == "--target") g_target = atoi(next());
        else if (a == "--threads") threads = atoi(next());
        else if (a == "--chunk") g_chunk = strtoull(next(), nullptr, 10);
        else if (a == "--checkpoint") g_ckpt = next();
        else if (a == "--out") g_out = next();
        else if (a == "--fresh") fresh = true;
        else {
            fprintf(stderr,
                    "usage: %s [--target N] [--threads T] [--chunk C] "
                    "[--checkpoint FILE] [--out FILE] [--fresh]\n",
                    argv[0]);
            return a == "--help" || a == "-h" ? 0 : 1;
        }
    }
    if (g_target < 0 || g_target > 15) { fprintf(stderr, "--target must be 0..15\n"); return 1; }
    build_tables();

    u64 cum = 0, cum_missing[10] = {0};
    int next_target = 0;
    std::vector<Row> rows;
    if (!fresh && load_checkpoint(g_merge_next, cum, cum_missing, next_target, rows)) {
        printf("resumed from %s: %llu chunks merged, %llu primes, next crossing 10^%d\n",
               g_ckpt.c_str(), (unsigned long long)g_merge_next, (unsigned long long)cum,
               next_target);
    }
    g_next_claim = g_merge_next;

    u64 final_target = pow10u(g_target);
    printf("digitprimes: first 10^%d primes, %d threads, chunk %llu\n", g_target, threads,
           (unsigned long long)g_chunk);
    fflush(stdout);

    std::signal(SIGINT, on_sigint);
    std::signal(SIGTERM, on_sigint);

    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) pool.emplace_back(worker, threads * 4);

    auto t0 = std::chrono::steady_clock::now();
    u64 cum0 = cum;
    auto last_status = t0, last_ckpt = t0;
    bool done = next_target > g_target;

    while (!done && !g_sigint) {
        ChunkResult r;
        {
            std::unique_lock<std::mutex> lk(g_mx);
            if (!g_cv_main.wait_for(lk, std::chrono::seconds(1),
                                    [&] { return g_results.count(g_merge_next) > 0; }))
                continue;  // timeout: recheck g_sigint
            r = g_results[g_merge_next];
            g_results.erase(g_merge_next);
        }
        // Does a 10^n crossing fall inside this chunk?
        if (next_target <= g_target && cum + r.nprimes >= pow10u(next_target)) {
            ChunkResult chk = rescan_chunk(g_merge_next, cum, cum_missing, next_target, rows);
            if (chk.nprimes != r.nprimes || memcmp(chk.missing, r.missing, sizeof chk.missing)) {
                fprintf(stderr, "internal error: rescan of chunk %llu disagrees with parallel pass\n",
                        (unsigned long long)g_merge_next);
                return 1;
            }
            write_results(rows);
        }
        cum += r.nprimes;
        for (int d = 0; d < 10; d++) cum_missing[d] += r.missing[d];
        {
            std::lock_guard<std::mutex> lk(g_mx);
            g_merge_next++;
        }
        g_cv_space.notify_all();
        done = next_target > g_target;

        auto now = std::chrono::steady_clock::now();
        double since_status = std::chrono::duration<double>(now - last_status).count();
        if (since_status >= 15.0) {
            double el = std::chrono::duration<double>(now - t0).count();
            double rate = (double)(cum - cum0) / el;
            double eta = rate > 0 ? (double)(final_target - cum) / rate : 0;
            printf("progress: %.4f%%  %llu primes  up to %llu  %.1fM primes/s  eta %.1f h\n",
                   100.0 * (double)cum / (double)final_target, (unsigned long long)cum,
                   (unsigned long long)(g_merge_next * g_chunk), rate / 1e6, eta / 3600.0);
            fflush(stdout);
            last_status = now;
        }
        if (std::chrono::duration<double>(now - last_ckpt).count() >= 60.0 || done) {
            save_checkpoint(g_merge_next, cum, cum_missing, next_target, rows);
            last_ckpt = now;
        }
    }

    g_stop = true;
    g_cv_space.notify_all();
    for (auto& t : pool) t.join();
    save_checkpoint(g_merge_next, cum, cum_missing, next_target, rows);
    write_results(rows);

    if (g_sigint && !done) {
        printf("\ninterrupted; checkpoint saved to %s (rerun to resume)\n", g_ckpt.c_str());
        return 130;
    }

    printf("\ndone. results in %s\n\nno digit d ('missing') sequences:\n", g_out.c_str());
    for (int d = 0; d < 10; d++) {
        printf("%s (no %d):", OEIS_MISS[d], d);
        for (auto& r : rows) printf(" %llu,", (unsigned long long)r.missing[d]);
        printf("\n");
    }
    printf("\nat least one digit d sequences (10^n - no_d, from n=1):\n");
    for (int d = 0; d < 10; d++) {
        printf("%s (>=1 %d):", OEIS_HAS[d], d);
        for (auto& r : rows)
            if (r.n >= 1) printf(" %llu,", (unsigned long long)(pow10u(r.n) - r.missing[d]));
        printf("\n");
    }
    return 0;
}
