/*
 * a287994 - search the decimal expansion of pi for OEIS A287994 / A290977.
 *
 *   A287994(n) = position of the first time an n-digit number appears twice
 *                in a row after the decimal point of pi (1-based, counting
 *                from the first digit after the decimal point).
 *   A290977(n) = the number itself (first n-digit number to appear twice in
 *                a row).
 *
 * Example: pi = 3.14159265358979323846264(33)... -> "3" repeats at position
 * 24, so A287994(1) = 24 and A290977(1) = 3.
 *
 * A block with a leading zero ("012012") is arguably not an "n-digit number".
 * The programs that produced the published terms do not exclude leading
 * zeros (for every known term the two definitions coincide), so this tool
 * tracks BOTH: the first repeat of any n-digit block, and the first repeat
 * whose block does not start with '0'.  If they ever differ for some n it
 * reports both and says so.
 *
 * Search: for a candidate start p we need d[i] == d[i+n] for all i in
 * [p, p+n).  The window is checked back to front; a mismatch at i kills
 * every start in [i-n+1, i], so the scan jumps straight to i+1.  With
 * random digits that is ~1.24 byte-compares per n digits advanced, so a
 * whole range of lengths costs only a few compares per digit and the scan
 * runs at memory-bandwidth speed.
 *
 * Threads grab fixed-size chunks of start positions in increasing order.
 * A per-length minimum is confirmed (and reported) once every chunk before
 * it has completed, so results are exact regardless of thread scheduling.
 * Lengths drop out of the scan as soon as they are resolved.
 *
 * The pi file is mmap'd read-only and never read up front: only the first
 * 64MB are sniffed to classify it.  A pure-digit file (y-cruncher / MIT
 * format) is searched in place, each worker validating its chunk's bytes as
 * it scans, so the file streams from disk exactly once - multi-TB files
 * work in bounded RAM (pages behind the confirmed frontier are dropped).
 * A file with junk near the start (line wraps etc.) is compacted into RAM,
 * which requires it to fit.
 *
 * Build:  cc -O3 -march=native -pthread -o a287994 a287994.c -lm
 */
#include "pi_common.h"
#include <pthread.h>
#include <stdatomic.h>

#define MAX_N 64

static struct {
    const char *pifile, *outfile;
    int nlo, nhi;            /* length range to search                  */
    uint64_t limit;          /* max digits to scan (UINT64_MAX = all)   */
    int threads;
    double progress_secs;
    uint64_t chunk;          /* start positions per work unit           */
} opt = { NULL, NULL, 0, 0, UINT64_MAX, 0, 0.5, 8u << 20 };

/* ------------------------------------------------------------------ */
/* Digits                                                              */
/* ------------------------------------------------------------------ */
static const uint8_t *D;     /* ASCII digits, D[0] = first digit after "." */
static uint64_t ND;          /* number of digits available                 */
static void *g_map;          /* whole-file mapping (may be NULL)           */
static size_t g_map_len;
static uint8_t *g_clean;     /* compacted copy when the file has junk      */

/* ------------------------------------------------------------------ */
/* Results                                                             */
/* ------------------------------------------------------------------ */
#define NOPOS UINT64_MAX
typedef struct {
    uint64_t any;            /* best candidate start (0-based), NOPOS=none */
    uint64_t clean;          /* same, but block must not start with '0'    */
    int any_final, clean_final;
    int reported;
} nres;
static nres R[MAX_N + 1];

static pthread_mutex_t g_lock = PTHREAD_MUTEX_INITIALIZER;
static uint8_t *g_chunk_done;
static uint64_t g_nchunks;
static uint64_t g_frontier;          /* chunks 0..g_frontier-1 all done   */
static int g_unresolved;             /* lengths not yet fully reported    */
static _Atomic uint64_t g_next_chunk;
static _Atomic int g_all_done;
static _Atomic uint64_t g_searched;  /* confirmed frontier in digits      */

static FILE *g_out;                  /* -o file (NULL if none)            */
static double g_t0;
static int g_tty;
static int g_status_shown;

static void clear_status(void)
{
    if (g_tty && g_status_shown) { fprintf(stderr, "\r\033[K"); g_status_shown = 0; }
}

static void usage(const char *argv0)
{
    fprintf(stderr,
"Usage: %s -p <pifile> -l <len | lo-hi> [options]\n"
"\n"
"Searches pi for the first n-digit block that appears twice in a row:\n"
"OEIS A287994 (position) and A290977 (value).\n"
"\n"
"  -p FILE     pi source text file (3.14159... ; non-digits are ignored)\n"
"  -l SPEC     block length n, or an inclusive range: -l 9  or  -l 1-10\n"
"  -n COUNT    digits of pi to scan (500M, 10B, all)  [all]\n"
"  -t N        worker threads  [all cores]\n"
"  -o FILE     append confirmed results + summary to FILE\n"
"  -P SECS     progress interval  [0.5, min 30 if stderr is not a tty]\n"
"  -C COUNT    work chunk size in digits  [8M]\n"
"  -h          this help\n", argv0);
    exit(1);
}

/* ------------------------------------------------------------------ */
/* Loading                                                             */
/* ------------------------------------------------------------------ */
/*
 * No upfront pass over the file: only the first SNIFF_BYTES are checked to
 * classify it.  A clean file (pure digits after "3.", the usual y-cruncher /
 * MIT format) is searched straight off the read-only mapping, with each
 * worker validating its chunk's bytes as it goes - the file streams from
 * disk exactly once, during the timed search.  A file with junk near the
 * start (line wraps etc.) is compacted into RAM as before, which is only
 * possible when it fits.  Junk appearing later in a sniffed-clean file is a
 * hard error with the offending offset.
 */
#define SNIFF_BYTES (64u << 20)

static int g_lazy;            /* 1 = workers validate chunk bytes as they scan */
static uint64_t g_sniffed;    /* digits already validated up front             */
static size_t g_skip;         /* file offset of the first digit                */
static uint64_t g_dropped;    /* file offset below which pages were discarded  */

static void validate_digits(const uint8_t *d, uint64_t s, uint64_t e)
{
    for (uint64_t i = s; i < e; i++) {
        if ((unsigned)(d[i] - '0') >= 10u)
            die("non-digit byte 0x%02x at digit position %" PRIu64 " - "
                "pi source is not pure digits; clean it first (e.g. tr -cd 0-9)",
                d[i], i + 1);
    }
}

static int is_ws(unsigned char c)
{
    return c == '\n' || c == '\r' || c == ' ' || c == '\t';
}

static void load_digits(void)
{
    int fd = open(opt.pifile, O_RDONLY);
    if (fd < 0) die("cannot open pi source '%s': %s", opt.pifile, strerror(errno));
    struct stat st;
    if (fstat(fd, &st) != 0 || !S_ISREG(st.st_mode))
        die("'%s' is not a regular file", opt.pifile);
    if (st.st_size == 0) die("pi source '%s' is empty", opt.pifile);
    size_t size = (size_t)st.st_size;

    unsigned char head[4] = {0};
    ssize_t hn = pread(fd, head, sizeof head, 0);
    if (hn < 0) die("read error on '%s': %s", opt.pifile, strerror(errno));
    size_t skip = 0;
    if (hn >= 2 && head[0] == '3' && head[1] == '.') skip = 2;
    else if (hn >= 4 && head[0] == '3' && head[1] == '1' && head[2] == '4' && head[3] == '1') skip = 1;
    else if (hn >= 3 && head[0] == '1' && head[1] == '4' && head[2] == '1') skip = 0;
    else warn("pi source does not start with '3.1415' or '1415' - assuming it begins at the first digit after the decimal point");

    /* Trim trailing whitespace (e.g. a final newline) by reading backwards. */
    uint64_t end = size;
    while (end > skip) {
        unsigned char tb[4096];
        size_t want = end - skip < sizeof tb ? (size_t)(end - skip) : sizeof tb;
        ssize_t rn = pread(fd, tb, want, (off_t)(end - want));
        if (rn <= 0) die("read error near end of '%s': %s", opt.pifile, strerror(errno));
        size_t k = (size_t)rn;
        while (k > 0 && is_ws(tb[k - 1])) k--;
        end -= (uint64_t)((size_t)rn - k);
        if (k > 0) break;
    }
    if (end <= skip) die("no digits in '%s'", opt.pifile);

    g_map = mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (g_map == MAP_FAILED) die("mmap of '%s' failed: %s", opt.pifile, strerror(errno));
    g_map_len = size;
    g_skip = skip;
    close(fd);
#ifdef MADV_SEQUENTIAL
    madvise(g_map, size, MADV_SEQUENTIAL);
#endif

    D = (const uint8_t *)g_map + skip;
    ND = end - skip;

    /* Sniff the head; a wrapped/dirty file shows junk almost immediately. */
    uint64_t sniff = ND < SNIFF_BYTES ? ND : SNIFF_BYTES;
    uint64_t bad = NOPOS;
    for (uint64_t i = 0; i < sniff; i++)
        if ((unsigned)(D[i] - '0') >= 10u) { bad = i; break; }
    if (bad == NOPOS) {
        g_lazy = 1;
        g_sniffed = sniff;
        return;
    }

    /* Junk near the start (line-wrapped file etc.): compact into RAM. */
    char b[32];
    uint64_t cap = ND < opt.limit ? ND : opt.limit;   /* only what -n will search */
    fprintf(stderr, "pi source has non-digit bytes (first at digit offset %" PRIu64 "); compacting %s into RAM...\n",
            bad + 1, fmt_bytes((double)cap, b, sizeof b));
    check_fits_in_ram(cap, "compacted pi digits");
    g_clean = alloc_huge(cap, 0, "compacted pi digits");
    uint64_t out = 0, scanned = 0;
    for (uint64_t k = 0; k < ND && out < cap; k++, scanned++) {
        unsigned c = (unsigned)(D[k] - '0');
        if (c < 10u) g_clean[out++] = D[k];
    }
    if (out < scanned / 2)
        warn("more than half of the pi source was non-digit bytes - wrong file?");
    munmap(g_map, g_map_len);
    g_map = NULL;
    D = g_clean;
    ND = out;
}

/* ------------------------------------------------------------------ */
/* Reporting (g_lock held)                                             */
/* ------------------------------------------------------------------ */
static void print_result_line(FILE *f, int n)
{
    char pb[40], vb[MAX_N + 1], vb2[MAX_N + 1];
    if (R[n].any == NOPOS) {
        fprintf(f, "n=%-2d  no repeat found in first %s digits\n",
                n, fmt_u64(ND, pb, sizeof pb));
        return;
    }
    memcpy(vb, D + R[n].any, (size_t)n); vb[n] = 0;
    if (R[n].any == R[n].clean) {
        fprintf(f, "n=%-2d  A287994(%d) = %" PRIu64 "  A290977(%d) = %s   (\"%s%s\" at position %s)\n",
                n, n, R[n].any + 1, n, vb, vb, vb, fmt_u64(R[n].any + 1, pb, sizeof pb));
    } else {
        fprintf(f, "n=%-2d  first repeat \"%s%s\" at position %" PRIu64 " has a LEADING ZERO;",
                n, vb, vb, R[n].any + 1);
        if (R[n].clean == NOPOS) {
            fprintf(f, " no leading-zero-free repeat in first %s digits\n",
                    fmt_u64(ND, pb, sizeof pb));
        } else {
            memcpy(vb2, D + R[n].clean, (size_t)n); vb2[n] = 0;
            fprintf(f, " first n-digit-number repeat is %s at position %" PRIu64 "\n",
                    vb2, R[n].clean + 1);
        }
        fprintf(f, "      (published A287994/A290977 do not exclude leading zeros -> "
                   "a(%d) = %" PRIu64 " / %s by that convention)\n", n, R[n].any + 1, vb);
    }
}

static void report_n(int n)
{
    R[n].reported = 1;
    g_unresolved--;
    clear_status();
    print_result_line(stdout, n);
    fflush(stdout);
    if (g_out) { print_result_line(g_out, n); fflush(g_out); }
    if (g_unresolved == 0) atomic_store(&g_all_done, 1);
}

/* Confirm minima: a candidate is final once every chunk before it is done. */
static void resolve_check(void)
{
    uint64_t fpos = g_frontier >= g_nchunks ? NOPOS : g_frontier * opt.chunk;
    atomic_store(&g_searched, g_frontier >= g_nchunks ? ND : fpos);
    for (int n = opt.nlo; n <= opt.nhi; n++) {
        if (R[n].reported) continue;
        if (!R[n].any_final   && (fpos == NOPOS || R[n].any   < fpos)) R[n].any_final = 1;
        if (!R[n].clean_final && (fpos == NOPOS || R[n].clean < fpos)) R[n].clean_final = 1;
        if (R[n].any_final && R[n].clean_final) report_n(n);
    }
}

static void record(int n, uint64_t p, int is_clean)
{
    pthread_mutex_lock(&g_lock);
    if (p < R[n].any) R[n].any = p;
    if (is_clean && p < R[n].clean) R[n].clean = p;
    resolve_check();
    pthread_mutex_unlock(&g_lock);
}

/* Digits behind the confirmed frontier are never read again (workers only
 * read forward from their chunk start), so tell the kernel to drop those
 * pages - a multi-TB file then streams through RAM instead of filling it.
 * (The summary re-reads a few bytes at each result; they just fault back in.) */
static void drop_scanned_pages(void)     /* g_lock held */
{
    if (!g_map) return;
    uint64_t fpos = g_frontier * opt.chunk;
    if (fpos > ND) fpos = ND;
    uint64_t new_end = ((uint64_t)g_skip + fpos) & ~((1ull << 20) - 1);
    if (new_end > g_dropped + (256ull << 20)) {
#ifdef MADV_DONTNEED
        madvise((char *)g_map + g_dropped, (size_t)(new_end - g_dropped), MADV_DONTNEED);
#endif
        g_dropped = new_end;
    }
}

static void mark_done(uint64_t c)
{
    pthread_mutex_lock(&g_lock);
    g_chunk_done[c] = 1;
    while (g_frontier < g_nchunks && g_chunk_done[g_frontier]) g_frontier++;
    drop_scanned_pages();
    resolve_check();
    pthread_mutex_unlock(&g_lock);
}

/* ------------------------------------------------------------------ */
/* Search                                                              */
/* ------------------------------------------------------------------ */
static void scan_chunk_n(int n, uint64_t s, uint64_t e, int need_any, int need_clean)
{
    if (ND < (uint64_t)(2 * n)) return;
    uint64_t cap = ND - (uint64_t)(2 * n) + 1;   /* exclusive bound on starts */
    if (e < cap) cap = e;
    const uint8_t *d = D;
    uint64_t p = s;
    while (p < cap) {
        /* need d[i] == d[i+n] for all i in [p, p+n); check back to front */
        uint64_t j = p + (uint64_t)n;
        while (j > p && d[j - 1] == d[j - 1 + (uint64_t)n]) j--;
        if (j > p) { p = j; continue; }          /* mismatch at j-1 kills starts <= j-1 */
        if (need_any) { record(n, p, 0); need_any = 0; }
        if (d[p] != '0') {
            if (need_clean) record(n, p, 1);
            need_clean = 0;
        }
        if (!need_any && !need_clean) return;
        p++;
    }
}

static void *worker(void *arg)
{
    (void)arg;
    for (;;) {
        if (atomic_load(&g_all_done) || g_stop) break;
        uint64_t c = atomic_fetch_add(&g_next_chunk, 1);
        if (c >= g_nchunks) break;
        uint64_t s = c * opt.chunk;
        uint64_t e = s + opt.chunk;
        if (e > ND) e = ND;

        if (g_lazy) {
            /* Validate the bytes this chunk can read (windows reach 2n-1
             * digits past the chunk end) before trusting the comparisons. */
            uint64_t ve = e + (uint64_t)(2 * opt.nhi) - 1;
            if (ve > ND) ve = ND;
            uint64_t vs = s < g_sniffed ? g_sniffed : s;
            if (vs < ve) validate_digits(D, vs, ve);
        }

        uint64_t snap_any[MAX_N + 1], snap_clean[MAX_N + 1];
        pthread_mutex_lock(&g_lock);
        for (int n = opt.nlo; n <= opt.nhi; n++) {
            snap_any[n] = R[n].any;
            snap_clean[n] = R[n].clean;
        }
        pthread_mutex_unlock(&g_lock);

        for (int n = opt.nlo; n <= opt.nhi; n++) {
            int need_any = snap_any[n] > s;      /* could this chunk beat it? */
            int need_clean = snap_clean[n] > s;
            if (need_any || need_clean)
                scan_chunk_n(n, s, e, need_any, need_clean);
        }
        mark_done(c);
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Progress                                                            */
/* ------------------------------------------------------------------ */
static void print_status(void)
{
    uint64_t done = atomic_load(&g_searched);
    double t = now_sec() - g_t0;
    double rate = t > 0 ? (double)done / t : 0;
    double pct = ND ? 100.0 * (double)done / (double)ND : 100.0;
    char b1[32], b2[32], b3[32], d1[32];
    char left[128] = "";
    int ln = 0, first = 1;
    pthread_mutex_lock(&g_lock);
    int unresolved = g_unresolved;
    for (int n = opt.nlo; n <= opt.nhi && ln < (int)sizeof left - 8; n++) {
        if (R[n].reported) continue;
        ln += snprintf(left + ln, sizeof left - (size_t)ln, "%s%d", first ? "" : ",", n);
        first = 0;
    }
    pthread_mutex_unlock(&g_lock);
    double eta = rate > 0 ? (double)(ND - done) / rate : 0;
    if (g_tty) {
        fprintf(stderr, "\r\033[K%5.1f%%  %s / %s digits  %s/s  ETA %s  searching n=%s ",
                pct, fmt_u64(done, b1, sizeof b1), fmt_u64(ND, b2, sizeof b2),
                fmt_u64((uint64_t)rate, b3, sizeof b3),
                fmt_duration(eta, d1, sizeof d1), unresolved ? left : "-");
        g_status_shown = 1;
    } else {
        fprintf(stderr, "%5.1f%%  %s / %s digits  %s/s  ETA %s  searching n=%s\n",
                pct, fmt_u64(done, b1, sizeof b1), fmt_u64(ND, b2, sizeof b2),
                fmt_u64((uint64_t)rate, b3, sizeof b3),
                fmt_duration(eta, d1, sizeof d1), unresolved ? left : "-");
    }
    fflush(stderr);
}

/* ------------------------------------------------------------------ */
/* Main                                                                */
/* ------------------------------------------------------------------ */
static void parse_lenspec(const char *s)
{
    char *end;
    long lo = strtol(s, &end, 10);
    long hi = lo;
    if (*end == '-') hi = strtol(end + 1, &end, 10);
    if (*end || lo < 1 || hi < lo || hi > MAX_N)
        die("bad length spec '%s' (want N or LO-HI, 1..%d)", s, MAX_N);
    opt.nlo = (int)lo;
    opt.nhi = (int)hi;
}

int main(int argc, char **argv)
{
    int c;
    while ((c = getopt(argc, argv, "p:l:n:t:o:P:C:h")) != -1) {
        switch (c) {
        case 'p': opt.pifile = optarg; break;
        case 'l': parse_lenspec(optarg); break;
        case 'n': if (parse_count(optarg, &opt.limit)) die("bad -n '%s'", optarg); break;
        case 't': opt.threads = atoi(optarg); break;
        case 'o': opt.outfile = optarg; break;
        case 'P': opt.progress_secs = atof(optarg); break;
        case 'C': if (parse_count(optarg, &opt.chunk) || opt.chunk == 0) die("bad -C '%s'", optarg); break;
        default: usage(argv[0]);
        }
    }
    if (!opt.pifile || !opt.nlo) usage(argv[0]);
    if (opt.threads <= 0) {
        long nc = sysconf(_SC_NPROCESSORS_ONLN);
        opt.threads = nc > 0 ? (int)nc : 1;
    }
    g_tty = isatty(2);
    if (!g_tty && opt.progress_secs < 30) opt.progress_secs = 30;
    install_signals();

    load_digits();
    if (opt.limit < ND) ND = opt.limit;
    if (ND < (uint64_t)(2 * opt.nlo))
        die("only %" PRIu64 " digits available - need at least %d", ND, 2 * opt.nlo);

    if (opt.outfile) {
        g_out = fopen(opt.outfile, "a");
        if (!g_out) die("cannot open output file '%s': %s", opt.outfile, strerror(errno));
    }

    char b1[32], b2[32];
    fprintf(stderr, "a287994: scanning %s digits of pi from '%s' for n=%d..%d, %d thread%s, %s-digit chunks\n",
            fmt_u64(ND, b1, sizeof b1), opt.pifile, opt.nlo, opt.nhi,
            opt.threads, opt.threads == 1 ? "" : "s",
            fmt_u64(opt.chunk, b2, sizeof b2));
    if (g_out) {
        fprintf(g_out, "# a287994 -p %s -l %d-%d -n %" PRIu64 " -t %d\n",
                opt.pifile, opt.nlo, opt.nhi, ND, opt.threads);
        fflush(g_out);
    }

    for (int n = 0; n <= MAX_N; n++) { R[n].any = R[n].clean = NOPOS; }
    g_unresolved = opt.nhi - opt.nlo + 1;
    g_nchunks = (ND + opt.chunk - 1) / opt.chunk;
    if (g_nchunks == 0) g_nchunks = 1;
    g_chunk_done = calloc((size_t)g_nchunks, 1);
    if (!g_chunk_done) die("out of memory for %" PRIu64 " chunk flags", g_nchunks);

    g_t0 = now_sec();
    pthread_t *tid = calloc((size_t)opt.threads, sizeof *tid);
    if (!tid) die("out of memory");
    for (int t = 0; t < opt.threads; t++)
        if (pthread_create(&tid[t], NULL, worker, NULL))
            die("pthread_create failed");

    /* Tick fast so completion is noticed promptly; print at -P intervals. */
    struct timespec ts = { 0, 100 * 1000 * 1000 };
    double next_print = now_sec() + opt.progress_secs;
    while (!atomic_load(&g_all_done) && !g_stop &&
           atomic_load(&g_searched) < ND) {
        nanosleep(&ts, NULL);
        if (now_sec() >= next_print) {
            print_status();
            next_print = now_sec() + opt.progress_secs;
        }
    }
    for (int t = 0; t < opt.threads; t++) pthread_join(tid[t], NULL);
    clear_status();

    double dt = now_sec() - g_t0;
    uint64_t searched = atomic_load(&g_searched);
    if (g_stop) fprintf(stderr, "interrupted - partial results up to %s digits are exact\n",
                        fmt_u64(searched, b1, sizeof b1));

    /* Anything never reported (only possible on interrupt): show best guess. */
    pthread_mutex_lock(&g_lock);
    for (int n = opt.nlo; n <= opt.nhi; n++) {
        if (R[n].reported) continue;
        clear_status();
        if (R[n].any == NOPOS)
            printf("n=%-2d  UNCONFIRMED: no repeat found in first %s searched digits\n",
                   n, fmt_u64(searched, b1, sizeof b1));
        else
            printf("n=%-2d  UNCONFIRMED (interrupted): best candidate so far at position %" PRIu64 "\n",
                   n, R[n].any + 1);
    }
    pthread_mutex_unlock(&g_lock);

    printf("\nSummary (positions are 1-based after the decimal point):\n");
    printf("  n    A287994 (position)    A290977 (value)\n");
    for (int n = opt.nlo; n <= opt.nhi; n++) {
        char vb[MAX_N + 1];
        if (R[n].any == NOPOS) {
            printf(" %2d    %-20s  not found in %s digits\n", n, "-",
                   fmt_u64(searched, b1, sizeof b1));
            continue;
        }
        memcpy(vb, D + R[n].any, (size_t)n); vb[n] = 0;
        printf(" %2d    %-20" PRIu64 "  %s%s%s\n", n, R[n].any + 1, vb,
               R[n].reported ? "" : "  (unconfirmed)",
               vb[0] == '0' ? "  [leading zero!]" : "");
    }
    fprintf(stderr, "scanned %s digits in %.1fs (%s digits/s), peak RSS %.0f MB\n",
            fmt_u64(searched, b1, sizeof b1), dt,
            fmt_u64(dt > 0 ? (uint64_t)((double)searched / dt) : 0, b2, sizeof b2),
            peak_rss_mb());
    if (g_out) {
        fprintf(g_out, "# done: scanned %" PRIu64 " digits in %.1fs\n", searched, dt);
        fclose(g_out);
    }
    if (g_map) munmap(g_map, g_map_len);
    return g_stop ? 130 : 0;
}
