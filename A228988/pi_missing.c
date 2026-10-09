/*
 * pi_missing.c — scan the digits of Pi for terms of several OEIS sequences:
 *
 *   A228988 (always):  a(n) = smallest missing number in the first 10^n
 *                      digits after the decimal point of Pi.
 *   A260627 (-L):      a(n) = largest n-digit number missing in the first
 *                      10^n digits after the decimal point of Pi.
 *   A153221 (-S):      numbers k such that the string k is found at
 *                      (1-based) position k-3 after the decimal point.
 *   A153223/A153224 (-S): same with position k-4 / k-5.
 *
 * The input file is a plain text file of Pi's digits.  A leading "3." is
 * auto-detected and skipped; everything after that must be pure digit
 * characters (trailing whitespace at end of file is ignored).  If your file
 * has line breaks, strip them first:  tr -cd '0-9' < in.txt > out.txt
 *
 * Method: slide a window over the digit stream.  For every position that
 * does not start with '0' (leading zeros don't count), take the value of the
 * 1-, 2-, ... digit windows starting there and set the corresponding bit in
 * a shared bit array, stopping once the value reaches the tracking limit.
 * At each 10^k boundary, the smallest clear bit is A228988(k).  Marking is
 * fanned out across threads; bits are set with an atomic OR (guarded by a
 * plain read, since almost all bits are already set).
 *
 * For A260627 the same pass also maintains, for each window length l, a
 * small bit array covering only the top -W numbers of [10^(l-1), 10^l)
 * (empirically the answer is within 10 of 10^n - 1, so the default window
 * of 2^20 is enormous margin).  At each 10^k boundary the highest clear bit
 * of the length-k window is A260627(k); if the window has no clear bit the
 * program says to enlarge -W rather than guess.
 *
 * For A153221 each thread carries the decimal string of i+4 as a running
 * counter (position i is 0-based, so 1-based position i+1 = k-3 gives
 * k = i+4) and compares it against the stream at every position; matches
 * are printed as they are found.
 *
 * Once a boundary reports smallest-missing = s, every value below s is
 * provably present, so later segments skip the bit-array access entirely
 * for window values < s.  In the final segment of a large run this cuts
 * the random memory accesses per position by roughly 10x, which is where
 * nearly all the time goes.
 *
 * The A228988 tracking limit defaults to max(N/25, 2e6), ~4x the empirical
 * answer (a(n) is roughly 10^n / 100).  If every tracked number turns out
 * to be present the program refuses to answer and tells you to raise -m;
 * it never reports a wrong term.  Memory use is limit/8 bytes (-m 4e10 ->
 * 5 GB).
 *
 * Build:  cc -O3 -march=native -pthread -o pi_missing pi_missing.c
 * Usage:  ./pi_missing [-L] [-S] [-n digits] [-m limit] [-W window]
 *                      [-t threads] [-o skip] pifile
 */

#define _GNU_SOURCE
#include <ctype.h>
#include <fcntl.h>
#include <signal.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

static const unsigned char *g_dig;   /* digit stream, g_n chars            */
static uint64_t  g_n;                /* number of digits used              */
static uint64_t  g_maxv;             /* track numbers in [0, g_maxv)       */
static int       g_lmax;             /* longest window marked anywhere     */
static uint64_t *g_bits;             /* g_maxv bits (A228988)              */
static uint64_t  g_badpos;           /* first non-digit position, or ~0    */
static uint64_t  g_smin;             /* all values < g_smin known present  */
static pthread_mutex_t g_lock = PTHREAD_MUTEX_INITIALIZER;

/* A260627 state */
static int       g_do_top;
static int       g_kmax;             /* largest k with 10^k <= g_n         */
static uint64_t  g_topw;             /* window bits per length (mult of 64)*/
static uint64_t  g_toplo[24];        /* window start per length, ~0 if off */
static uint64_t *g_topbits;          /* g_kmax blocks of g_topw bits       */

/* Self-locating searches: k whose string occurs at 1-based position k-off */
static int       g_do_self;
#define NSELF 3
static const struct { int off; const char *name; } g_selfseq[NSELF] = {
    {3, "A153221"}, {4, "A153223"}, {5, "A153224"}
};

static uint64_t g_pow10[20];

/* Checkpointing */
static const char *g_ckpt;           /* checkpoint file, or NULL           */
static volatile sig_atomic_t g_stop_flag;

static void on_signal(int sig) { (void)sig; g_stop_flag = 1; }

struct ckpt {
    char     magic[8];
    uint64_t n, maxv, topw;
    int64_t  kmax, do_top, do_self, bi;
    uint64_t segstart, lo, done, cursor, smin;
    unsigned char digsig[32];        /* first digits, guards wrong file    */
    double   elapsed;
};
#define CKPT_MAGIC "PIMISS02"

static inline void mark(uint64_t v)
{
    uint64_t w = v >> 6, m = 1ULL << (v & 63);
    if (!(g_bits[w] & m))
        __atomic_fetch_or(&g_bits[w], m, __ATOMIC_RELAXED);
}

static inline void marktop(int len, uint64_t v)
{
    uint64_t b = (uint64_t)(len - 1) * g_topw + (v - g_toplo[len]);
    uint64_t w = b >> 6, m = 1ULL << (b & 63);
    if (!(g_topbits[w] & m))
        __atomic_fetch_or(&g_topbits[w], m, __ATOMIC_RELAXED);
}

/* Increment a decimal string in place. */
static inline void incbuf(char *buf, int *len)
{
    int p = *len - 1;
    while (p >= 0 && buf[p] == '9') buf[p--] = '0';
    if (p < 0) { memmove(buf + 1, buf, (size_t)*len); buf[0] = '1'; (*len)++; }
    else buf[p]++;
}

struct job { uint64_t lo, hi, wend, selfstart; };

static void *worker(void *arg)
{
    struct job *j = arg;
    const unsigned char *dig = g_dig;
    const uint64_t maxv = g_maxv, wend = j->wend, n = g_n;
    const uint64_t smin = g_smin;    /* values below this are known present */
    const int lmax = g_lmax;
    const int toplen = g_do_top ? g_kmax : 0;

    char self[NSELF][24];
    int selflen[NSELF] = {0};
    if (g_do_self) {
        /* counters start at the first position actually checked:
         * at 0-based i (1-based p = i+1), candidate k = p + off = i+1+off */
        uint64_t first = j->lo > j->selfstart ? j->lo : j->selfstart;
        for (int s = 0; s < NSELF; s++)
            selflen[s] = snprintf(self[s], sizeof self[s], "%" PRIu64,
                                  first + 1 + (uint64_t)g_selfseq[s].off);
    }

    for (uint64_t i = j->lo; i < j->hi; i++) {
        unsigned d = (unsigned)dig[i] - '0';
        if (d > 9) {
            pthread_mutex_lock(&g_lock);
            if (i < g_badpos) g_badpos = i;
            pthread_mutex_unlock(&g_lock);
            return NULL;
        }

        if (g_do_self && i >= j->selfstart) {
            for (int s = 0; s < NSELF; s++) {
                int len = selflen[s];
                if (i + (uint64_t)len <= n
                    && dig[i] == (unsigned char)self[s][0]) {
                    int m = 1;
                    while (m < len && dig[i + m] == (unsigned char)self[s][m])
                        m++;
                    if (m == len) {
                        pthread_mutex_lock(&g_lock);
                        printf("%s: %.*s   (at position %" PRIu64 ")\n",
                               g_selfseq[s].name, len, self[s], i + 1);
                        fflush(stdout);
                        pthread_mutex_unlock(&g_lock);
                    }
                }
                incbuf(self[s], &selflen[s]);
            }
        }

        if (d == 0) {                        /* leading zeros don't count */
            if (!smin) mark(0);
            if (toplen && g_toplo[1] == 0) marktop(1, 0);
            continue;
        }
        uint64_t v = 0;
        uint64_t lim = i + (uint64_t)lmax;
        if (lim > wend) lim = wend;
        int l = 0;
        for (uint64_t p = i; p < lim; p++) {
            v = v * 10 + ((unsigned)dig[p] - '0');
            l++;
            if (v < maxv && v >= smin) mark(v);
            if (v >= g_toplo[l]) marktop(l, v);       /* ~0 when inactive */
            if (v >= maxv && l >= toplen) break;      /* v only grows     */
        }
    }
    return NULL;
}

/* Mark all windows starting in [lo,hi), clamped to end before wend.
 * Positions below selfstart skip the A153221 check (already done). */
static void process_range(uint64_t lo, uint64_t hi, uint64_t wend,
                          uint64_t selfstart, int nthreads)
{
    if (hi <= lo) return;
    uint64_t span = hi - lo;
    int t = nthreads;
    if (span < (uint64_t)t * 4096) t = 1;

    pthread_t tid[t];
    struct job job[t];
    for (int k = 0; k < t; k++) {
        job[k].lo        = lo + span * k / t;
        job[k].hi        = lo + span * (k + 1) / t;
        job[k].wend      = wend;
        job[k].selfstart = selfstart;
        if (pthread_create(&tid[k], NULL, worker, &job[k])) {
            perror("pthread_create");
            exit(1);
        }
    }
    for (int k = 0; k < t; k++)
        pthread_join(tid[k], NULL);
}

/* Smallest clear bit >= from, or g_maxv if all of [from, g_maxv) is set. */
static uint64_t find_missing(uint64_t from)
{
    uint64_t nw = (g_maxv + 63) >> 6;
    uint64_t wi = from >> 6;
    if (wi >= nw) return g_maxv;
    uint64_t cur = g_bits[wi];
    if (from & 63) cur |= (1ULL << (from & 63)) - 1;   /* ignore bits < from */
    for (;;) {
        if (cur != ~0ULL) {
            uint64_t v = (wi << 6) + (uint64_t)__builtin_ctzll(~cur);
            return v < g_maxv ? v : g_maxv;
        }
        if (++wi >= nw) return g_maxv;
        cur = g_bits[wi];
    }
}

/* Highest clear bit in the length-k top window, or ~0 if all set. */
static uint64_t find_top_missing(int k)
{
    uint64_t lo = g_toplo[k];
    uint64_t nbits = g_pow10[k] - lo;
    uint64_t *bits = g_topbits + (uint64_t)(k - 1) * (g_topw >> 6);
    int64_t last = (int64_t)((nbits + 63) >> 6) - 1;
    for (int64_t wi = last; wi >= 0; wi--) {
        uint64_t w = bits[wi];
        if (wi == last && (nbits & 63))
            w |= ~((1ULL << (nbits & 63)) - 1);       /* pad past 10^k */
        if (w != ~0ULL)
            return lo + ((uint64_t)wi << 6)
                      + (uint64_t)(63 - __builtin_clzll(~w));
    }
    return UINT64_MAX;
}

static size_t ckpt_bits_words(void)    { return (g_maxv + 63) / 64; }
static size_t ckpt_top_words(void)     { return (size_t)g_kmax * (g_topw >> 6); }

/* Write header + bit arrays to <g_ckpt>.tmp, then rename atomically. */
static int ckpt_write(int bi, uint64_t segstart, uint64_t lo, uint64_t done,
                      uint64_t cursor, double elapsed)
{
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", g_ckpt);
    FILE *f = fopen(tmp, "wb");
    if (!f) { perror(tmp); return -1; }
    struct ckpt h;
    memset(&h, 0, sizeof h);
    memcpy(h.magic, CKPT_MAGIC, 8);
    h.n = g_n; h.maxv = g_maxv; h.topw = g_topw;
    h.kmax = g_kmax; h.do_top = g_do_top; h.do_self = g_do_self;
    h.bi = bi; h.segstart = segstart; h.lo = lo; h.done = done;
    h.cursor = cursor; h.smin = g_smin; h.elapsed = elapsed;
    memcpy(h.digsig, g_dig, g_n < 32 ? (size_t)g_n : 32);
    int ok = fwrite(&h, sizeof h, 1, f) == 1
          && fwrite(g_bits, 8, ckpt_bits_words(), f) == ckpt_bits_words();
    if (ok && g_do_top)
        ok = fwrite(g_topbits, 8, ckpt_top_words(), f) == ckpt_top_words();
    if (ok) { fflush(f); fsync(fileno(f)); }
    if (fclose(f) != 0) ok = 0;
    if (!ok || rename(tmp, g_ckpt) != 0) {
        fprintf(stderr, "\ncheckpoint write to %s failed\n", tmp);
        return -1;
    }
    return 0;
}

/* Load a checkpoint into the already-allocated arrays; dies on mismatch. */
static void ckpt_read(struct ckpt *h)
{
    FILE *f = fopen(g_ckpt, "rb");
    if (!f) { perror(g_ckpt); exit(1); }
    unsigned char sig[32] = {0};
    memcpy(sig, g_dig, g_n < 32 ? (size_t)g_n : 32);
    if (fread(h, sizeof *h, 1, f) != 1
        || memcmp(h->magic, CKPT_MAGIC, 8) != 0
        || h->n != g_n || h->maxv != g_maxv || h->topw != g_topw
        || h->kmax != g_kmax || h->do_top != g_do_top
        || h->do_self != g_do_self
        || memcmp(h->digsig, sig, 32) != 0) {
        fprintf(stderr, "checkpoint %s does not match this input file and "
                        "option set (-n/-m/-W/-L/-S must be identical)\n",
                g_ckpt);
        exit(1);
    }
    int ok = fread(g_bits, 8, ckpt_bits_words(), f) == ckpt_bits_words();
    if (ok && g_do_top)
        ok = fread(g_topbits, 8, ckpt_top_words(), f) == ckpt_top_words();
    fclose(f);
    if (!ok) { fprintf(stderr, "checkpoint %s truncated\n", g_ckpt); exit(1); }
}

static uint64_t parse_size(const char *s)
{
    char *end;
    if (strpbrk(s, "eE.")) {                 /* allow 1e12 style */
        double d = strtod(s, &end);
        if (*end || d < 0 || d > 9e18) goto bad;
        return (uint64_t)(d + 0.5);
    }
    uint64_t v = strtoull(s, &end, 10);
    switch (*end) {
    case 'k': case 'K': v *= 1000ULL; end++; break;
    case 'm': case 'M': v *= 1000000ULL; end++; break;
    case 'g': case 'G': v *= 1000000000ULL; end++; break;
    case 't': case 'T': v *= 1000000000000ULL; end++; break;
    }
    if (*end) goto bad;
    return v;
bad:
    fprintf(stderr, "bad size argument: %s\n", s);
    exit(1);
}

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void fmt_dur(char *buf, size_t sz, double secs)
{
    if (secs < 0) secs = 0;
    long s = (long)secs;
    if (s >= 86400)
        snprintf(buf, sz, "%ldd %02ld:%02ld", s / 86400,
                 s % 86400 / 3600, s % 3600 / 60);
    else
        snprintf(buf, sz, "%ld:%02ld:%02ld", s / 3600, s % 3600 / 60, s % 60);
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "usage: %s [-L] [-S] [-n digits] [-m limit] [-W window] [-t threads]\n"
        "          [-o skip] pifile\n"
        "  -L          also compute A260627 (largest n-digit number missing in\n"
        "              first 10^n digits) at each power-of-10 boundary\n"
        "  -S          also search for self-locating strings: A153221 (k at\n"
        "              position k-3), A153223 (k-4), A153224 (k-5)\n"
        "  -n digits   use only the first <digits> digits (default: whole file;\n"
        "              suffixes k/M/G/T and 1e12 notation accepted)\n"
        "  -m limit    A228988: track numbers 0..limit-1, memory = limit/8\n"
        "              bytes (default max(digits/25, 2e6))\n"
        "  -W window   A260627: per-length window below 10^n (default 2^20)\n"
        "  -t threads  worker threads (default: number of CPUs)\n"
        "  -o skip     skip <skip> leading bytes of the file (a leading \"3.\"\n"
        "              is skipped automatically)\n"
        "  -c file     write a checkpoint to <file> every -i seconds and on\n"
        "              SIGINT/SIGTERM (which then exit cleanly, status 3)\n"
        "  -i secs     checkpoint interval (default 1800)\n"
        "  -r          resume from the -c checkpoint file\n", prog);
    exit(1);
}

int main(int argc, char **argv)
{
    uint64_t want_n = 0, maxv = 0, skip = 0, topw = 1ULL << 20;
    int nthreads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (nthreads < 1) nthreads = 1;
    int do_top = 0, do_self = 0, resume = 0;
    uint64_t ck_interval = 1800;

    int opt;
    while ((opt = getopt(argc, argv, "LSrn:m:W:t:o:c:i:h")) != -1) {
        switch (opt) {
        case 'L': do_top = 1; break;
        case 'S': do_self = 1; break;
        case 'r': resume = 1; break;
        case 'n': want_n = parse_size(optarg); break;
        case 'm': maxv = parse_size(optarg); break;
        case 'W': topw = parse_size(optarg); break;
        case 't': nthreads = atoi(optarg); break;
        case 'o': skip = parse_size(optarg); break;
        case 'c': g_ckpt = optarg; break;
        case 'i': ck_interval = parse_size(optarg); break;
        default: usage(argv[0]);
        }
    }
    if (optind != argc - 1) usage(argv[0]);
    if (resume && !g_ckpt) {
        fprintf(stderr, "-r needs a checkpoint file (-c)\n");
        return 1;
    }
    if (nthreads < 1 || nthreads > 1024) {
        fprintf(stderr, "bad thread count\n");
        return 1;
    }

    g_pow10[0] = 1;
    for (int i = 1; i < 20; i++) g_pow10[i] = g_pow10[i - 1] * 10;

    /* Map the file. */
    int fd = open(argv[optind], O_RDONLY);
    if (fd < 0) { perror(argv[optind]); return 1; }
    struct stat st;
    if (fstat(fd, &st) < 0) { perror("fstat"); return 1; }
    if (st.st_size == 0) { fprintf(stderr, "empty file\n"); return 1; }
    unsigned char *map = mmap(NULL, (size_t)st.st_size, PROT_READ,
                              MAP_PRIVATE, fd, 0);
    if (map == MAP_FAILED) { perror("mmap"); return 1; }
    madvise(map, (size_t)st.st_size, MADV_SEQUENTIAL);

    /* Skip a leading "3." (digits *after* the decimal point are what count) */
    if (skip == 0 && st.st_size >= 2 && map[0] == '3' && map[1] == '.')
        skip = 2;
    if ((int64_t)skip >= st.st_size) {
        fprintf(stderr, "skip larger than file\n");
        return 1;
    }
    uint64_t avail = (uint64_t)st.st_size - skip;
    while (avail > 0 && isspace(map[skip + avail - 1]))
        avail--;                              /* trailing newline etc. */

    g_dig = map + skip;
    g_n = want_n ? want_n : avail;
    if (g_n > avail) {
        fprintf(stderr, "file has only %" PRIu64 " digits, %" PRIu64
                        " requested\n", avail, g_n);
        return 1;
    }

    g_maxv = maxv ? maxv : (g_n / 25 > 2000000 ? g_n / 25 : 2000000);
    if (g_maxv > 1000000000000000000ULL) {
        fprintf(stderr, "-m limit too large\n");
        return 1;
    }
    g_lmax = 1;
    for (uint64_t x = g_maxv - 1; x >= 10; x /= 10) g_lmax++;
    g_badpos = UINT64_MAX;
    g_do_self = do_self;

    /* A260627 setup: one top window per boundary length 1..kmax. */
    g_kmax = 0;
    while (g_kmax + 1 < 19 && g_pow10[g_kmax + 1] <= g_n) g_kmax++;
    g_do_top = do_top && g_kmax >= 1;
    for (int l = 0; l < 24; l++) g_toplo[l] = UINT64_MAX;
    if (g_do_top) {
        g_topw = (topw + 63) & ~63ULL;        /* round to whole words */
        if (g_topw < 64) g_topw = 64;
        for (int l = 1; l <= g_kmax; l++) {
            uint64_t lower = (l == 1) ? 0 : g_pow10[l - 1];
            uint64_t lo = g_pow10[l] > g_topw ? g_pow10[l] - g_topw : 0;
            g_toplo[l] = lo > lower ? lo : lower;
        }
        if (g_kmax > g_lmax) g_lmax = g_kmax; /* need k-digit windows too */
        g_topbits = calloc((size_t)g_kmax * (g_topw >> 6), 8);
        if (!g_topbits) { fprintf(stderr, "window allocation failed\n"); return 1; }
    }

    fprintf(stderr, "digits: %" PRIu64 "   tracking 0..%" PRIu64
                    " (%.2f GB)   threads: %d%s%s\n",
            g_n, g_maxv - 1, (double)g_maxv / 8 / 1e9, nthreads,
            g_do_top ? "   +A260627" : "", g_do_self ? "   +A153221" : "");

    g_bits = calloc((g_maxv + 63) / 64, 8);
    if (!g_bits) { fprintf(stderr, "bit array allocation failed\n"); return 1; }

    /* Boundaries: powers of 10 up to g_n, then g_n itself. */
    uint64_t bound[24];
    int nb = 0;
    for (uint64_t b = 10; b <= g_n && nb < 20; b *= 10) {
        bound[nb++] = b;
        if (b > g_n / 10) break;              /* avoid overflow */
    }
    if (nb == 0 || bound[nb - 1] != g_n) bound[nb++] = g_n;

    const uint64_t SLICE = 1ULL << 30;        /* progress/abort granularity */
    double t0 = now();
    uint64_t done = 0, cursor = 0, segstart = 0, resume_lo = 0;
    int bi0 = 0, resumed = 0;

    if (g_ckpt) {
        struct sigaction sa;
        memset(&sa, 0, sizeof sa);
        sa.sa_handler = on_signal;
        sa.sa_flags = SA_RESTART;
        sigaction(SIGINT, &sa, NULL);
        sigaction(SIGTERM, &sa, NULL);
    }
    if (resume) {
        struct ckpt h;
        ckpt_read(&h);
        if (h.bi < 0 || h.bi >= nb || h.lo > h.n) {
            fprintf(stderr, "checkpoint state out of range\n");
            return 1;
        }
        bi0 = (int)h.bi; segstart = h.segstart; resume_lo = h.lo;
        done = h.done; cursor = h.cursor; g_smin = h.smin;
        resumed = 1;
        t0 = now() - h.elapsed;
        /* Boundaries already reported: their top windows stay disabled. */
        if (g_do_top)
            for (int j = 0; j < bi0; j++) {
                uint64_t p = bound[j]; int k = 0;
                while (p % 10 == 0 && p > 1) { p /= 10; k++; }
                if (p == 1 && k <= g_kmax) g_toplo[k] = UINT64_MAX;
            }
        fprintf(stderr, "resumed from %s at position %" PRIu64
                        " (%.1f%%, %.0fs elapsed)\n",
                g_ckpt, resume_lo, 100.0 * done / (g_n + 0.0), h.elapsed);
    }
    double last_ck = now();

    /* Status-line state: EMA of recent throughput for rate/ETA display. */
    int    st_tty = isatty(fileno(stderr));
    double st_prev_t = now(), st_ema = 0, st_last_line = 0;
    uint64_t st_prev_done = done;

    for (int bi = bi0; bi < nb; bi++) {
        uint64_t B = bound[bi];
        /* Redo the last few positions of the previous segment: their
         * windows were clamped at the old boundary and may now extend. */
        uint64_t lo = segstart;
        if (lo > 0) lo -= (lo < (uint64_t)(g_lmax - 1)) ? lo : (uint64_t)(g_lmax - 1);
        if (resumed) { lo = resume_lo; resumed = 0; }

        while (lo < B) {
            uint64_t hi = lo + SLICE < B ? lo + SLICE : B;
            process_range(lo, hi, B, segstart, nthreads);
            if (g_badpos != UINT64_MAX) {
                fprintf(stderr, "\nnon-digit byte 0x%02x at file offset %"
                        PRIu64 " -- clean the file with: tr -cd '0-9'\n",
                        g_dig[g_badpos], skip + g_badpos);
                return 1;
            }
            done += hi - lo;
            lo = hi;
            if (g_n > SLICE) {
                double tn = now(), dt = tn - st_prev_t;
                if (dt > 0.01) {
                    double inst = (double)(done - st_prev_done) / dt;
                    st_ema = st_ema > 0 ? 0.7 * st_ema + 0.3 * inst : inst;
                    st_prev_t = tn;
                    st_prev_done = done;
                }
                /* tty: redraw in place; log file: one full line per minute */
                if (st_ema > 0 && (st_tty || tn - st_last_line >= 60)) {
                    char el[32], etab[32], etae[32], fin[32], blabel[24];
                    fmt_dur(el, sizeof el, tn - t0);
                    fmt_dur(etab, sizeof etab, (double)(B - lo) / st_ema);
                    fmt_dur(etae, sizeof etae, (double)(g_n - lo) / st_ema);
                    time_t ft = time(NULL)
                              + (time_t)((double)(g_n - lo) / st_ema);
                    strftime(fin, sizeof fin, "%b %d %H:%M", localtime(&ft));
                    uint64_t p = B; int k = 0;
                    while (p % 10 == 0 && p > 1) { p /= 10; k++; }
                    snprintf(blabel, sizeof blabel,
                             p == 1 ? "10^%d" : "end", k);
                    fprintf(stderr, "%s%5.1f%%  %" PRIu64 "/%" PRIu64
                            "G  %.1f Md/s  elapsed %s  %s in %s"
                            "  ETA %s (%s)%s",
                            st_tty ? "\r" : "",
                            100.0 * done / (g_n + 0.0),
                            lo / 1000000000, g_n / 1000000000,
                            st_ema / 1e6, el,
                            blabel, etab, etae, fin,
                            st_tty ? "      " : "\n");
                    st_last_line = tn;
                }
            }
            if (g_ckpt && (g_stop_flag || now() - last_ck >= ck_interval)) {
                if (ckpt_write(bi, segstart, lo, done, cursor,
                               now() - t0) == 0 && g_stop_flag) {
                    fprintf(stderr, "\ncheckpoint saved to %s at position %"
                            PRIu64 "; resume with -r\n", g_ckpt, lo);
                    return 3;
                }
                last_ck = now();
            }
        }
        if (g_n > SLICE && st_tty) fprintf(stderr, "\r%110s\r", "");

        cursor = find_missing(cursor);
        g_smin = cursor;   /* later segments can skip marking values < cursor */
        if (cursor >= g_maxv) {
            fprintf(stderr, "every number below %" PRIu64 " is present in the "
                    "first %" PRIu64 " digits;\nrerun with a larger -m limit\n",
                    g_maxv, B);
            return 2;
        }
        /* Label power-of-10 boundaries as a(k). */
        uint64_t p = B; int k = 0;
        while (p % 10 == 0 && p > 1) { p /= 10; k++; }
        if (p == 1) {
            printf("A228988(%d) = %" PRIu64 "   (first 10^%d digits, %.1fs)\n",
                   k, cursor, k, now() - t0);
            if (g_do_top && k <= g_kmax) {
                uint64_t top = find_top_missing(k);
                if (top == UINT64_MAX)
                    printf("A260627(%d): every number in [%" PRIu64 ", 10^%d) "
                           "is present -- rerun with a larger -W window\n",
                           k, g_toplo[k], k);
                else
                    printf("A260627(%d) = %" PRIu64 "\n", k, top);
                g_toplo[k] = UINT64_MAX;   /* boundary reported; stop marking */
            }
        } else {
            printf("A228988: smallest missing in first %" PRIu64 " digits: %"
                   PRIu64 "   (%.1fs)\n", B, cursor, now() - t0);
        }
        fflush(stdout);
        segstart = B;
    }
    return 0;
}
