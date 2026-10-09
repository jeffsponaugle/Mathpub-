/*
 * pi_common.h - shared helpers for pibloom / pisearch (A197123 tools)
 *
 * Conventions (match OEIS A197123 and the original tools):
 *   - Digit positions are 1-based, counting from the first digit AFTER the
 *     decimal point.  "1 appears at positions 1 and 3", "26 at 6 and 21".
 *   - The leading "3." in the source file is skipped and never part of a window.
 *   - A "window" of length L at position p is the L digits starting at p.
 */
#ifndef PI_COMMON_H
#define PI_COMMON_H

#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <signal.h>
#include <math.h>
#include <inttypes.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <stdarg.h>

#define MAX_SEQ_LEN 38          /* two uint64 decimal halves: 19 + 19 digits */
#define LO_DIGITS   19          /* digits held in the 'lo' word              */

/* ------------------------------------------------------------------ */
/* Error / logging                                                     */
/* ------------------------------------------------------------------ */
static void die(const char *fmt, ...) __attribute__((noreturn, format(printf, 1, 2)));
static void die(const char *fmt, ...)
{
    va_list ap;
    fflush(stdout);
    fprintf(stderr, "\nFATAL: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    exit(2);
}
static void warn(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
static void warn(const char *fmt, ...)
{
    va_list ap;
    fflush(stdout);
    fprintf(stderr, "WARNING: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
}

/* ------------------------------------------------------------------ */
/* Time                                                                */
/* ------------------------------------------------------------------ */
static inline double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static const char *fmt_duration(double s, char *buf, size_t n)
{
    if (s < 0) s = 0;
    long t = (long)s;
    long d = t / 86400, h = (t / 3600) % 24, m = (t / 60) % 60, sec = t % 60;
    if (d > 0)      snprintf(buf, n, "%ldd %02ld:%02ld:%02ld", d, h, m, sec);
    else if (h > 0) snprintf(buf, n, "%ld:%02ld:%02ld", h, m, sec);
    else            snprintf(buf, n, "%ld:%02ld", m, sec);
    return buf;
}

/* Thousands-separated integer: 1234567 -> "1,234,567" */
static const char *fmt_u64(uint64_t v, char *buf, size_t n)
{
    char tmp[32];
    int len = snprintf(tmp, sizeof tmp, "%" PRIu64, v);
    int commas = (len - 1) / 3;
    int out = len + commas;
    if ((size_t)out + 1 > n) { snprintf(buf, n, "%" PRIu64, v); return buf; }
    buf[out] = 0;
    int j = out - 1;
    for (int i = len - 1, c = 0; i >= 0; i--) {
        buf[j--] = tmp[i];
        if (++c == 3 && i > 0) { buf[j--] = ','; c = 0; }
    }
    return buf;
}

/* Human readable byte size */
static const char *fmt_bytes(double b, char *buf, size_t n)
{
    const char *u[] = {"B", "KiB", "MiB", "GiB", "TiB", "PiB"};
    int i = 0;
    while (b >= 1024.0 && i < 5) { b /= 1024.0; i++; }
    snprintf(buf, n, "%.2f %s", b, u[i]);
    return buf;
}

/* ------------------------------------------------------------------ */
/* Size parsing                                                        */
/* ------------------------------------------------------------------ */
/* Memory sizes: 1024-based.  "512M", "1G", "1.5T", "800000000" */
static int parse_mem_size(const char *s, uint64_t *out)
{
    char *end;
    errno = 0;
    double v = strtod(s, &end);
    if (end == s || errno != 0 || v < 0) return -1;
    double mult = 1;
    if (*end) {
        switch (*end) {
        case 'k': case 'K': mult = 1024.0; break;
        case 'm': case 'M': mult = 1024.0 * 1024; break;
        case 'g': case 'G': mult = 1024.0 * 1024 * 1024; break;
        case 't': case 'T': mult = 1024.0 * 1024 * 1024 * 1024; break;
        default: return -1;
        }
        end++;
        if (*end == 'B' || *end == 'b') end++;     /* allow "MB", "GiB" */
        if (*end == 'i') end++;
        if (*end == 'B' || *end == 'b') end++;
        if (*end) return -1;
    }
    double r = v * mult;
    if (r > 1.8e19) return -1;
    *out = (uint64_t)r;
    return 0;
}

/* Digit counts: 1000-based. "500M", "10B", "1T", "all" -> UINT64_MAX */
static int parse_count(const char *s, uint64_t *out)
{
    if (!strcmp(s, "all") || !strcmp(s, "ALL")) { *out = UINT64_MAX; return 0; }
    char *end;
    errno = 0;
    double v = strtod(s, &end);
    if (end == s || errno != 0 || v < 0) return -1;
    double mult = 1;
    if (*end) {
        switch (*end) {
        case 'k': case 'K': mult = 1e3; break;
        case 'm': case 'M': mult = 1e6; break;
        case 'b': case 'B': case 'g': case 'G': mult = 1e9; break;
        case 't': case 'T': mult = 1e12; break;
        default: return -1;
        }
        end++;
        if (*end) return -1;
    }
    double r = v * mult;
    if (r > 1.8e19) return -1;
    *out = (uint64_t)r;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Hashing                                                             */
/* ------------------------------------------------------------------ */
/* Key is the 128-bit pair (hi, lo): lo = last 19 digits as a decimal number,
 * hi = preceding (L-19) digits (0 when L <= 19).  This encoding is injective
 * for a fixed L, so a good 64-bit mixer gives near-ideal bloom behaviour for
 * every sequence length with no per-length special casing. */
static inline uint64_t mix64(uint64_t x)
{
    /* splitmix64 finalizer */
    x ^= x >> 30; x *= 0xBF58476D1CE4E5B9ULL;
    x ^= x >> 27; x *= 0x94D049BB133111EBULL;
    x ^= x >> 31;
    return x;
}

static inline void hash_window(uint64_t lo, uint64_t hi, uint64_t *h1, uint64_t *h2)
{
    uint64_t a = mix64(lo + 0x9E3779B97F4A7C15ULL);
    uint64_t b = mix64(hi ^ a ^ 0xD1B54A32D192ED03ULL);
    *h1 = a ^ (b << 1);
    *h2 = mix64(b + a) | 1;          /* odd: good stride for double hashing */
}

/* Map a uniform 64-bit value into [0, range) without a division (Lemire). */
static inline uint64_t fastrange64(uint64_t h, uint64_t range)
{
    return (uint64_t)(((__uint128_t)h * (__uint128_t)range) >> 64);
}

/* ------------------------------------------------------------------ */
/* Window (hi,lo) rolling state                                        */
/* ------------------------------------------------------------------ */
static const uint64_t POW10[20] = {
    1ULL, 10ULL, 100ULL, 1000ULL, 10000ULL, 100000ULL, 1000000ULL, 10000000ULL,
    100000000ULL, 1000000000ULL, 10000000000ULL, 100000000000ULL, 1000000000000ULL,
    10000000000000ULL, 100000000000000ULL, 1000000000000000ULL, 10000000000000000ULL,
    100000000000000000ULL, 1000000000000000000ULL, 10000000000000000000ULL };

typedef struct {
    int L;          /* sequence length                */
    int Lh;         /* digits in hi (0 if L <= 19)     */
    int Ll;         /* digits in lo (min(L,19))        */
    uint64_t hi_top;/* POW10[Lh-1] (0 if Lh == 0)      */
    uint64_t lo_top;/* POW10[Ll-1]                     */
} winspec;

static void winspec_init(winspec *w, int L)
{
    if (L < 1 || L > MAX_SEQ_LEN) die("sequence length must be 1..%d", MAX_SEQ_LEN);
    w->L = L;
    w->Ll = L < LO_DIGITS ? L : LO_DIGITS;
    w->Lh = L - w->Ll;
    w->lo_top = POW10[w->Ll - 1];
    w->hi_top = w->Lh ? POW10[w->Lh - 1] : 0;
}

/* Compute (hi,lo) from scratch for digits d[0..L-1] (values 0..9). */
static inline void win_compute(const winspec *w, const uint8_t *d, uint64_t *hi, uint64_t *lo)
{
    uint64_t h = 0, l = 0;
    int i = 0;
    for (; i < w->Lh; i++) h = h * 10 + d[i];
    for (; i < w->L; i++)  l = l * 10 + d[i];
    *hi = h; *lo = l;
}

/* Roll from window at d (i.e. d[0..L-1]) to window at d+1. */
static inline void win_roll(const winspec *w, const uint8_t *d, uint64_t *hi, uint64_t *lo)
{
    if (w->Lh) {
        *hi = (*hi - (uint64_t)d[0] * w->hi_top) * 10 + d[w->Lh];
    }
    *lo = (*lo - (uint64_t)d[w->Lh] * w->lo_top) * 10 + d[w->L];
}

/* Parse L ASCII digits into (hi,lo). Returns -1 on non-digit. */
static int win_parse_ascii(const winspec *w, const char *s, uint64_t *hi, uint64_t *lo)
{
    uint8_t d[MAX_SEQ_LEN];
    for (int i = 0; i < w->L; i++) {
        if (s[i] < '0' || s[i] > '9') return -1;
        d[i] = (uint8_t)(s[i] - '0');
    }
    win_compute(w, d, hi, lo);
    return 0;
}

/* Render (hi,lo) back to an L-digit ASCII string (NUL terminated). */
static void win_to_ascii(const winspec *w, uint64_t hi, uint64_t lo, char *out)
{
    for (int i = w->L - 1; i >= w->Lh; i--) { out[i] = (char)('0' + lo % 10); lo /= 10; }
    for (int i = w->Lh - 1; i >= 0; i--)   { out[i] = (char)('0' + hi % 10); hi /= 10; }
    out[w->L] = 0;
}

/* ------------------------------------------------------------------ */
/* Pi digit reader                                                     */
/* ------------------------------------------------------------------ */
/*
 * Reads the pi text file sequentially in large chunks, strips anything that
 * is not a decimal digit (newlines, spaces, the "3."), converts to values
 * 0..9, and hands back a buffer that always starts with the last (L-1)
 * digits of the previous chunk so that every window is contiguous.
 *
 * After pireader_fill():
 *    r->buf[0 .. r->ndigits-1]   cleaned digits (ndigits includes overlap)
 *    r->base                     1-based position of r->buf[0]
 *    windows available this chunk = r->ndigits - (L-1)  (if >= 1)
 * Windows at buffer index i correspond to absolute position r->base + i.
 */
typedef struct {
    int fd;
    int L;
    size_t chunk;           /* raw bytes read per fill            */
    uint8_t *raw;           /* raw read buffer                     */
    uint8_t *buf;           /* cleaned digits, with overlap prefix */
    size_t ndigits;         /* valid digits in buf                 */
    uint64_t base;          /* absolute 1-based position of buf[0] */
    uint64_t next_pos;      /* position of next digit to be read   */
    uint64_t total_raw;     /* raw bytes consumed                  */
    uint64_t total_skipped; /* non-digit bytes dropped             */
    uint64_t file_size;     /* bytes, or 0 if unknown              */
    int eof;
    int started;
    uint64_t limit;         /* max digits to deliver (UINT64_MAX=all) */
    uint64_t delivered;     /* digits delivered so far             */
} pireader;

static void pireader_open(pireader *r, const char *path, int L, size_t chunk, uint64_t limit_digits)
{
    memset(r, 0, sizeof *r);
    r->L = L;
    r->chunk = chunk;
    r->limit = limit_digits;
    r->fd = open(path, O_RDONLY);
    if (r->fd < 0) die("cannot open pi source '%s': %s", path, strerror(errno));
    struct stat st;
    if (fstat(r->fd, &st) == 0 && S_ISREG(st.st_mode)) r->file_size = (uint64_t)st.st_size;
#ifdef POSIX_FADV_SEQUENTIAL
    posix_fadvise(r->fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
#ifdef F_RDAHEAD
    fcntl(r->fd, F_RDAHEAD, 1);
#endif
    r->raw = malloc(chunk);
    r->buf = malloc(chunk + (size_t)L + 64);
    if (!r->raw || !r->buf) die("cannot allocate %zu byte read buffers", chunk);

    /* Peek the first bytes to skip a leading "3." (or "3") */
    unsigned char head[4];
    ssize_t n = read(r->fd, head, sizeof head);
    if (n < 0) die("read error on '%s': %s", path, strerror(errno));
    if (n == 0) die("pi source '%s' is empty", path);
    off_t skip = 0;
    if (n >= 2 && head[0] == '3' && head[1] == '.') skip = 2;
    else if (n >= 4 && head[0] == '3' && head[1] == '1' && head[2] == '4' && head[3] == '1') skip = 1;
    else if (n >= 3 && head[0] == '1' && head[1] == '4' && head[2] == '1') skip = 0;
    else warn("pi source does not start with '3.1415' or '1415' - assuming it begins at the first digit after the decimal point");
    if (lseek(r->fd, skip, SEEK_SET) < 0) die("lseek failed on '%s': %s", path, strerror(errno));
    r->total_raw = (uint64_t)skip;
    r->next_pos = 1;
    r->base = 1;
}

/* Returns number of new windows available (0 at EOF / limit). */
static size_t pireader_fill(pireader *r)
{
    const int L = r->L;
    size_t keep = 0;
    if (r->started && r->ndigits >= (size_t)(L - 1)) {
        keep = (size_t)(L - 1);
        memmove(r->buf, r->buf + r->ndigits - keep, keep);
        r->base = r->next_pos - keep;
    } else if (r->started) {
        keep = r->ndigits;
        r->base = r->next_pos - keep;
    }
    r->started = 1;
    r->ndigits = keep;

    while (r->ndigits - keep < r->chunk - 4096 && !r->eof) {
        if (r->delivered >= r->limit) { r->eof = 1; break; }
        size_t want = r->chunk - (r->ndigits - keep);
        if (want > r->chunk) want = r->chunk;
        ssize_t n = read(r->fd, r->raw, want);
        if (n < 0) {
            if (errno == EINTR) continue;
            die("read error on pi source: %s", strerror(errno));
        }
        if (n == 0) { r->eof = 1; break; }
        r->total_raw += (uint64_t)n;
        uint8_t *dst = r->buf + r->ndigits;
        const uint8_t *src = r->raw;
        size_t out = 0;
        uint64_t room = r->limit - r->delivered;
        for (ssize_t i = 0; i < n; i++) {
            unsigned c = (unsigned)src[i] - '0';
            if (c < 10u) {
                if ((uint64_t)out >= room) { r->eof = 1; break; }
                dst[out++] = (uint8_t)c;
            } else {
                r->total_skipped++;
            }
        }
        r->ndigits += out;
        r->delivered += out;
        r->next_pos += out;
        /* A buffer full of junk (e.g. wrong file) must not spin forever */
        if (out == 0 && r->total_skipped > (uint64_t)r->chunk * 4)
            die("pi source contains no digits after %" PRIu64 " bytes - wrong file?", r->total_raw);
    }
    if (r->ndigits < (size_t)L) return 0;
    return r->ndigits - (size_t)(L - 1);
}

static void pireader_close(pireader *r)
{
    if (r->fd >= 0) close(r->fd);
    free(r->raw); free(r->buf);
    r->raw = r->buf = NULL; r->fd = -1;
}

/* ------------------------------------------------------------------ */
/* Large memory allocation                                             */
/* ------------------------------------------------------------------ */
static void *alloc_huge(uint64_t bytes, int prefault, const char *what)
{
    if (bytes == 0) return NULL;
    void *p = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS
#ifdef MAP_NORESERVE
                   | (prefault ? 0 : MAP_NORESERVE)
#endif
                   , -1, 0);
    if (p == MAP_FAILED) {
        char b[32];
        die("cannot allocate %s for %s: %s", fmt_bytes((double)bytes, b, sizeof b), what, strerror(errno));
    }
#ifdef MADV_HUGEPAGE
    madvise(p, bytes, MADV_HUGEPAGE);
#endif
    if (prefault) {
        /* Touch every page now so page faults don't land in the hot loop and
         * so an over-committed allocation fails here, not 3 days in. */
        char b[32], b2[32], b3[32], d1[32];
        printf("Pre-faulting %s for %s ...", fmt_bytes((double)bytes, b, sizeof b), what);
        fflush(stdout);
        int tty = isatty(1);
        double t0 = now_sec(), tlast = t0;
        int started_progress = 0, last_decile = 0;
        const uint64_t step = 256ULL << 20;              /* 256 MiB per touch */
        for (uint64_t off = 0; off < bytes; ) {
            uint64_t n = bytes - off < step ? bytes - off : step;
            memset((char *)p + off, 0, n);
            off += n;
            double t = now_sec();
            /* only bother with live progress once this is taking a while */
            if (t - t0 < 1.0) continue;
            double pct = 100.0 * (double)off / (double)bytes;
            double rate = (double)off / (t - t0);
            double eta = rate > 0 ? (double)(bytes - off) / rate : 0;
            if (tty) {
                if (t - tlast < 0.25 && off < bytes) continue;
                if (!started_progress) { printf("\n"); started_progress = 1; }
                printf("\r\033[K  %5.1f%%  %s of %s   %s/s   ETA %s ", pct,
                       fmt_bytes((double)off, b, sizeof b), fmt_bytes((double)bytes, b2, sizeof b2),
                       fmt_bytes(rate, b3, sizeof b3), fmt_duration(eta, d1, sizeof d1));
                fflush(stdout);
                tlast = t;
            } else {
                int decile = (int)(pct / 10.0);
                if (decile > last_decile) {
                    if (!started_progress) { printf("\n"); started_progress = 1; }
                    printf("  %3.0f%%  %s   %s/s   ETA %s\n", 10.0 * decile,
                           fmt_bytes((double)off, b, sizeof b), fmt_bytes(rate, b3, sizeof b3),
                           fmt_duration(eta, d1, sizeof d1));
                    last_decile = decile;
                }
            }
        }
        double dt = now_sec() - t0;
        if (started_progress && tty) printf("\r\033[K");
        if (started_progress)
            printf("  pre-fault done: %s in %.1fs (%s/s)\n", fmt_bytes((double)bytes, b, sizeof b), dt,
                   fmt_bytes(dt > 0 ? (double)bytes / dt : 0, b2, sizeof b2));
        else
            printf(" done in %.1fs\n", dt);
    }
    return p;
}

static uint64_t phys_mem_bytes(void)
{
    long pages = sysconf(_SC_PHYS_PAGES), psz = sysconf(_SC_PAGESIZE);
    if (pages <= 0 || psz <= 0) return 0;
    return (uint64_t)pages * (uint64_t)psz;
}

/* Warn when an allocation is unlikely to stay resident. Random access into
 * swapped/compressed memory is 1000x slower than RAM - the run would never finish. */
static void check_fits_in_ram(uint64_t bytes, const char *what)
{
    uint64_t phys = phys_mem_bytes();
    if (!phys) return;
    char b1[32], b2[32];
    if (bytes > phys)
        die("%s (%s) is larger than physical RAM (%s)", what, fmt_bytes((double)bytes, b1, sizeof b1), fmt_bytes((double)phys, b2, sizeof b2));
    if (bytes > phys * 0.85)
        warn("%s (%s) is more than 85%% of physical RAM (%s); if anything else is using memory the run will swap and crawl",
             what, fmt_bytes((double)bytes, b1, sizeof b1), fmt_bytes((double)phys, b2, sizeof b2));
}

static double peak_rss_mb(void)
{
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    return (double)ru.ru_maxrss / (1024.0 * 1024.0);
#else
    return (double)ru.ru_maxrss / 1024.0;
#endif
}

/* ------------------------------------------------------------------ */
/* Signals                                                             */
/* ------------------------------------------------------------------ */
static volatile sig_atomic_t g_stop = 0;
static void on_signal(int sig)
{
    (void)sig;
    if (g_stop) _exit(130);   /* second Ctrl-C: hard exit */
    g_stop = 1;
}
static void install_signals(void)
{
    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_signal;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);
}

#endif /* PI_COMMON_H */
