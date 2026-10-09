/*
 * a076106 - compute OEIS A076106 and A076130 from a file of Pi digits.
 *
 * A076106(n): of all n-digit primes, the one whose first occurrence in the
 *             decimal digits of Pi (ignoring the initial 3) is latest.
 * A076130(n): the (1-indexed) position where that prime first appears.
 *
 * Known terms (a(8), a(9) computed by this tool, Aug 2026):
 *   A076106: 7, 73, 373, 9337, 35569, 805289, 9271903, 43927427, 342263843
 *   A076130: 13, 299, 5229, 75961, 715492, 11137824, 135224164,
 *            1541659153, 20252853413
 *
 * The digits file is streamed, never loaded into memory, so arbitrarily
 * large files (billions/trillions of digits) work. Non-digit characters
 * (whitespace, "3." prefix) are ignored; a leading 3 is dropped so
 * positions are counted from the first digit after the decimal point.
 *
 * All requested n are searched in one pass over the file, and reading
 * stops as soon as every n-digit prime has been seen for every n.
 *
 * Checkpointing (-c statefile): progress (per-n seen bitmaps + counters)
 * is saved periodically, at end of input, and on SIGINT/SIGTERM. A later
 * run with the same -c file resumes where the scan left off, fast-skipping
 * the digits file to the saved position -- so the scan can be continued
 * against a longer digits file without redoing the work. The digits file
 * used on resume must contain the same digit stream from the beginning
 * (formatting/whitespace may differ).
 *
 * Memory: two bitmaps of ~10^maxN bits (prime sieve + seen flags), about
 * 2.8 * 10^maxN / 8 bytes total: ~28 MB for maxN=8, ~2.8 GB for maxN=10,
 * ~28 GB for maxN=11, ~280 GB for maxN=12 (needs a large-memory server).
 *
 * Usage:
 *   a076106 [-n maxN] [-d maxdigits] [-t threads] [-c statefile]
 *           [-i ckptdigits] [-o outfile] pifile
 *
 * Build: cc -O3 -pthread -o a076106 a076106.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#include <time.h>

/* Results log (-o): discoveries as they happen, plus the final table. */
static FILE *outf;

static const char *timestamp(void)
{
    static char ts[32];
    time_t t = time(NULL);
    strftime(ts, sizeof ts, "%Y-%m-%d %H:%M:%S", localtime(&t));
    return ts;
}

#define MAXN 12

/* ---------------- Prime sieve (segmented, multithreaded bitmap) ---------- */

static uint8_t *sieve_bits;  /* bit i set => i is composite (or 0/1) */
static uint64_t sieve_limit; /* bits valid for 0 <= i < sieve_limit */
static uint32_t *base_primes;
static uint32_t nbase;

static inline int is_prime(uint64_t x)
{
    return !(sieve_bits[x >> 3] & (1u << (x & 7)));
}

/* Numbers per segment; a multiple of 8 so segments are byte-aligned. */
#define SEG_NUMS ((uint64_t)8 * (4u << 20)) /* 4 MB of bits per segment */

struct sieve_job {
    uint64_t lo, hi; /* [lo, hi) range of numbers, lo % SEG_NUMS == 0 */
};

static void *sieve_thread(void *arg)
{
    struct sieve_job *job = arg;
    for (uint64_t lo = job->lo; lo < job->hi; lo += SEG_NUMS) {
        uint64_t hi = lo + SEG_NUMS < job->hi ? lo + SEG_NUMS : job->hi;
        for (uint32_t b = 0; b < nbase; b++) {
            uint64_t p = base_primes[b];
            if (p * p >= hi)
                break;
            uint64_t j = p * p;
            if (j < lo)
                j = ((lo + p - 1) / p) * p;
            for (; j < hi; j += p)
                sieve_bits[j >> 3] |= (uint8_t)(1u << (j & 7));
        }
    }
    return NULL;
}

static void build_sieve(uint64_t limit, int nthreads)
{
    sieve_limit = limit;
    /* +16: count_primes reads the bitmap as 64-bit words (little-endian). */
    sieve_bits = calloc(limit / 8 + 16, 1);
    if (!sieve_bits) {
        fprintf(stderr, "error: out of memory for sieve to %" PRIu64 "\n", limit);
        exit(1);
    }

    /* Base primes up to sqrt(limit), by a small byte sieve. */
    uint64_t r = 2;
    while (r * r < limit)
        r++;
    uint8_t *small = calloc(r + 1, 1);
    base_primes = malloc((r / 2 + 1) * sizeof *base_primes);
    if (!small || !base_primes) {
        fprintf(stderr, "error: out of memory for base primes\n");
        exit(1);
    }
    for (uint64_t i = 2; i <= r; i++) {
        if (small[i])
            continue;
        base_primes[nbase++] = (uint32_t)i;
        for (uint64_t j = i * i; j <= r; j += i)
            small[j] = 1;
    }
    free(small);

    /* Split whole segments across threads. */
    uint64_t nsegs = (limit + SEG_NUMS - 1) / SEG_NUMS;
    uint64_t per = (nsegs + (uint64_t)nthreads - 1) / (uint64_t)nthreads;
    pthread_t tid[256];
    struct sieve_job jobs[256];
    int started = 0;
    for (int t = 0; t < nthreads; t++) {
        uint64_t s0 = (uint64_t)t * per, s1 = s0 + per;
        if (s0 >= nsegs)
            break;
        if (s1 > nsegs)
            s1 = nsegs;
        jobs[t].lo = s0 * SEG_NUMS;
        jobs[t].hi = s1 * SEG_NUMS < limit ? s1 * SEG_NUMS : limit;
        if (pthread_create(&tid[t], NULL, sieve_thread, &jobs[t])) {
            fprintf(stderr, "error: pthread_create failed\n");
            exit(1);
        }
        started++;
    }
    for (int t = 0; t < started; t++)
        pthread_join(tid[t], NULL);

    sieve_bits[0] |= 3; /* 0 and 1 are not prime */
}

/* Count primes in [lo, hi) by popcounting composite bits word-wise. */
static uint64_t count_primes(uint64_t lo, uint64_t hi)
{
    uint64_t composites = 0, w0 = lo / 64, w1 = hi / 64;
    const uint64_t *words = (const uint64_t *)sieve_bits;
    if (w0 == w1) {
        uint64_t m = ((hi - lo) == 64 ? ~0ULL : ((1ULL << (hi - lo)) - 1))
                     << (lo % 64);
        return (hi - lo) - (uint64_t)__builtin_popcountll(words[w0] & m);
    }
    if (lo % 64)
        composites += (uint64_t)__builtin_popcountll(words[w0++] & ~((1ULL << (lo % 64)) - 1));
    for (uint64_t w = w0; w < w1; w++)
        composites += (uint64_t)__builtin_popcountll(words[w]);
    if (hi % 64)
        composites += (uint64_t)__builtin_popcountll(words[w1] & ((1ULL << (hi % 64)) - 1));
    return (hi - lo) - composites;
}

static uint64_t popcount_bytes(const uint8_t *p, size_t nbytes)
{
    /* p is calloc'd with >= 8 bytes of zero padding past nbytes. */
    uint64_t total = 0;
    const uint64_t *w = (const uint64_t *)p;
    for (size_t k = 0; k < (nbytes + 7) / 8; k++)
        total += (uint64_t)__builtin_popcountll(w[k]);
    return total;
}

/* ---------------- Per-n search state ------------------------------------- */

struct search {
    int n;               /* digit length */
    uint64_t low;        /* smallest n-digit number (1 for n=1) */
    uint64_t p10n1;      /* 10^(n-1), weight of the outgoing digit */
    uint64_t v;          /* current window value */
    uint8_t *seen;       /* bitmap: prime already appeared */
    uint64_t seen_bytes; /* bitmap size (10^n/8 + 1) */
    uint64_t nprimes;    /* primes in [low, 10^n) */
    uint64_t remaining;  /* primes not yet seen */
    uint64_t last_prime; /* most recent first-occurrence */
    uint64_t last_pos;   /* its 1-indexed position */
};

static struct search S[MAXN];
static struct search *active[MAXN];
static int nactive;

/* Ring buffer of recent digits for the rolling windows (size > MAXN). */
#define RING 16
static uint8_t ring[RING];

/*
 * Candidate queue: windows that passed the cheap composite filters are
 * queued (with their sieve/seen bytes prefetched) and resolved in order
 * once the queue fills, hiding DRAM latency on the big bitmaps.
 */
#define QCAP 1024
static struct cand {
    struct search *s;
    uint64_t v, pos;
} q[QCAP];
static int qn;

static void flush_queue(void)
{
    for (int c = 0; c < qn; c++) {
        struct search *s = q[c].s;
        uint64_t v = q[c].v;
        if (!is_prime(v))
            continue;
        if (s->seen[v >> 3] & (1u << (v & 7)))
            continue;
        s->seen[v >> 3] |= (uint8_t)(1u << (v & 7));
        s->last_prime = v;
        s->last_pos = q[c].pos;
        if (--s->remaining == 0) {
            fprintf(stderr,
                    "n=%d done: a(%d)=%" PRIu64 " at position %" PRIu64
                    " (all %" PRIu64 " primes seen)\n",
                    s->n, s->n, s->last_prime, s->last_pos, s->nprimes);
            if (outf) {
                fprintf(outf,
                        "[%s] n=%d done: a(%d)=%" PRIu64 " at position %"
                        PRIu64 " (all %" PRIu64 " primes seen)\n",
                        timestamp(), s->n, s->n, s->last_prime, s->last_pos,
                        s->nprimes);
                fflush(outf);
            }
            for (int a = 0; a < nactive; a++)
                if (active[a] == s) {
                    memmove(&active[a], &active[a + 1],
                            (size_t)(nactive - a - 1) * sizeof active[0]);
                    nactive--;
                    break;
                }
        }
    }
    qn = 0;
}

/* ---------------- Checkpointing ------------------------------------------ */

#define CKPT_MAGIC "A076106S"
#define CKPT_VERSION 1

struct ckpt_header {
    char magic[8];
    uint32_t version;
    uint32_t maxN;
    uint32_t dropped3; /* did the original run drop a leading 3? */
    uint32_t pad;
    uint64_t digits;   /* digits scanned so far */
};

struct ckpt_search {
    uint64_t nprimes, remaining, last_prime, last_pos, seen_bytes;
};

static const char *ckpt_path;
static int g_dropped3;

static void save_checkpoint(int maxN, uint64_t digits)
{
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", ckpt_path);
    FILE *f = fopen(tmp, "wb");
    if (!f) {
        fprintf(stderr, "warning: cannot write checkpoint %s: %s\n",
                tmp, strerror(errno));
        return;
    }
    struct ckpt_header h;
    memset(&h, 0, sizeof h);
    memcpy(h.magic, CKPT_MAGIC, 8);
    h.version = CKPT_VERSION;
    h.maxN = (uint32_t)maxN;
    h.dropped3 = (uint32_t)g_dropped3;
    h.digits = digits;
    int ok = fwrite(&h, sizeof h, 1, f) == 1;
    for (int n = 1; ok && n <= maxN; n++) {
        struct search *s = &S[n - 1];
        struct ckpt_search cs = { s->nprimes, s->remaining, s->last_prime,
                                  s->last_pos, s->seen_bytes };
        ok = fwrite(&cs, sizeof cs, 1, f) == 1 &&
             fwrite(s->seen, 1, s->seen_bytes, f) == s->seen_bytes;
    }
    ok = !fclose(f) && ok;
    if (!ok || rename(tmp, ckpt_path)) {
        fprintf(stderr, "warning: checkpoint save to %s failed: %s\n",
                ckpt_path, strerror(errno));
        return;
    }
    fprintf(stderr, "checkpoint saved to %s at digit %" PRIu64 "\n",
            ckpt_path, digits);
}

/* Returns digits scanned; exits on a checkpoint inconsistent with args. */
static uint64_t load_checkpoint(int maxN)
{
    FILE *f = fopen(ckpt_path, "rb");
    if (!f)
        return 0; /* no checkpoint yet: fresh start */
    struct ckpt_header h;
    if (fread(&h, sizeof h, 1, f) != 1 || memcmp(h.magic, CKPT_MAGIC, 8) ||
        h.version != CKPT_VERSION) {
        fprintf(stderr, "error: %s is not a valid checkpoint file\n", ckpt_path);
        exit(1);
    }
    if (h.maxN != (uint32_t)maxN) {
        fprintf(stderr, "error: checkpoint %s was made with -n %u, not -n %d\n",
                ckpt_path, h.maxN, maxN);
        exit(1);
    }
    for (int n = 1; n <= maxN; n++) {
        struct search *s = &S[n - 1];
        struct ckpt_search cs;
        if (fread(&cs, sizeof cs, 1, f) != 1 ||
            cs.seen_bytes != s->seen_bytes ||
            fread(s->seen, 1, s->seen_bytes, f) != s->seen_bytes) {
            fprintf(stderr, "error: checkpoint %s is truncated or corrupt\n",
                    ckpt_path);
            exit(1);
        }
        if (cs.nprimes != s->nprimes) {
            fprintf(stderr, "error: checkpoint %s prime count mismatch for "
                    "n=%d (%" PRIu64 " vs %" PRIu64 ")\n",
                    ckpt_path, n, cs.nprimes, s->nprimes);
            exit(1);
        }
        if (popcount_bytes(s->seen, s->seen_bytes) !=
            s->nprimes - cs.remaining) {
            fprintf(stderr, "error: checkpoint %s seen-bitmap inconsistent "
                    "for n=%d\n", ckpt_path, n);
            exit(1);
        }
        s->remaining = cs.remaining;
        s->last_prime = cs.last_prime;
        s->last_pos = cs.last_pos;
    }
    fclose(f);
    g_dropped3 = (int)h.dropped3;
    fprintf(stderr, "resuming from %s: %" PRIu64 " digits already scanned\n",
            ckpt_path, h.digits);
    return h.digits;
}

/*
 * Fast-skip the digit stream to `target` digits (same filtering rules as
 * the scan), keeping the ring buffer filled with the trailing digits so
 * the rolling windows can be rebuilt. Returns 1 on success.
 */
static int skip_digits(FILE *f, uint64_t target)
{
    if (target == 0)
        return 1;
    static unsigned char buf[4 << 20];
    uint64_t i = 0, next_note = 100000000000ULL;
    int first = 1;
    size_t got;
    while ((got = fread(buf, 1, sizeof buf, f)) > 0) {
        for (size_t k = 0; k < got; k++) {
            unsigned d = (unsigned)buf[k] - '0';
            if (d > 9)
                continue;
            if (first) {
                first = 0;
                if (g_dropped3) {
                    if (d != 3)
                        fprintf(stderr, "warning: resume file does not start "
                                "with 3, but the checkpointed run's did\n");
                    else
                        continue;
                }
            }
            ring[i & (RING - 1)] = (uint8_t)d;
            i++;
            if (i == target) {
                /* Rewind the unconsumed tail of this buffer so the scan
                 * resumes at exactly the next byte. */
                if (fseeko(f, -(off_t)(got - k - 1), SEEK_CUR)) {
                    perror("fseeko");
                    return 0;
                }
                return 1;
            }
        }
        if (i >= next_note) {
            fprintf(stderr, "skipped %" PRIu64 " of %" PRIu64 " digits...\n",
                    i, target);
            next_note += 100000000000ULL;
        }
    }
    fprintf(stderr, "error: digits file ended at %" PRIu64 " digits, before "
            "the checkpoint position %" PRIu64 "\n", i, target);
    return 0;
}

/* ---------------- Signals ------------------------------------------------ */

static volatile sig_atomic_t g_stop;

static void on_signal(int sig)
{
    (void)sig;
    g_stop = 1;
}

/* ---------------- Main --------------------------------------------------- */

static void print_table(FILE *out, int maxN)
{
    fprintf(out, "%-3s %-16s %-18s %s\n", "n", "A076106(n)", "A076130(n)",
            "status");
    for (int n = 1; n <= maxN; n++) {
        struct search *s = &S[n - 1];
        if (s->remaining == 0)
            fprintf(out, "%-3d %-16" PRIu64 " %-18" PRIu64 " all %" PRIu64
                    " primes found\n",
                    n, s->last_prime, s->last_pos, s->nprimes);
        else
            fprintf(out, "%-3d %-16s %-18s %" PRIu64 " of %" PRIu64
                    " primes not yet seen -- need more digits\n",
                    n, "?", "?", s->remaining, s->nprimes);
    }
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "usage: %s [-n maxN] [-d maxdigits] [-t threads] [-c statefile]\n"
        "       %*s [-i ckptdigits] [-o outfile] pifile\n"
        "\n"
        "Computes OEIS A076106 and A076130: for each digit length n, the\n"
        "n-digit prime whose first appearance in the decimal digits of Pi\n"
        "(ignoring the initial 3) is latest, and the position where it appears.\n"
        "\n"
        "  -n maxN      compute a(1)..a(maxN), 1..%d (default 8)\n"
        "  -d maxdigits stop after scanning this many digits (default: whole file)\n"
        "  -t threads   threads for prime sieving (default: all CPUs)\n"
        "  -c statefile checkpoint file: resumed from if it exists; progress is\n"
        "               saved there periodically, at end of input, and on\n"
        "               SIGINT/SIGTERM. Resume must use the same -n, and a\n"
        "               digits file with the same digit stream from the start.\n"
        "  -i ckptdigits digits between periodic checkpoints (default 100e9,\n"
        "               0 = only at end/signal)\n"
        "  -o outfile   append discoveries (timestamped, as they happen) and\n"
        "               the final results table to this file\n"
        "  pifile       text file of Pi digits; streamed, so any size works.\n"
        "               A leading \"3.\" and all whitespace are ignored.\n"
        "\n"
        "Memory: ~2.8*10^maxN/8 bytes of bitmaps -- ~28 MB (maxN=8),\n"
        "~2.8 GB (maxN=10), ~28 GB (maxN=11), ~280 GB (maxN=12).\n"
        "\n"
        "example: %s -n 10 -c state10.ckpt pi-1t.txt\n",
        prog, (int)strlen(prog), "", MAXN, prog);
    exit(2);
}

int main(int argc, char **argv)
{
    int maxN = 8;
    int nthreads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    uint64_t maxdigits = UINT64_MAX;
    uint64_t ckpt_interval = 100000000000ULL;
    const char *pifile = NULL, *out_path = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-n") && i + 1 < argc)
            maxN = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-d") && i + 1 < argc)
            maxdigits = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-t") && i + 1 < argc)
            nthreads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-c") && i + 1 < argc)
            ckpt_path = argv[++i];
        else if (!strcmp(argv[i], "-i") && i + 1 < argc)
            ckpt_interval = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-o") && i + 1 < argc)
            out_path = argv[++i];
        else if (argv[i][0] != '-' && !pifile)
            pifile = argv[i];
        else
            usage(argv[0]);
    }
    if (!pifile)
        usage(argv[0]);
    if (maxN < 1 || maxN > MAXN) {
        fprintf(stderr, "error: maxN must be 1..%d\n", MAXN);
        return 2;
    }
    if (nthreads < 1)
        nthreads = 1;
    if (nthreads > 256)
        nthreads = 256;

    FILE *f = fopen(pifile, "rb");
    if (!f) {
        perror(pifile);
        return 1;
    }

    if (out_path) {
        outf = fopen(out_path, "a");
        if (!outf) {
            perror(out_path);
            return 1;
        }
        fprintf(outf, "[%s] run started: -n %d, digits file %s%s%s\n",
                timestamp(), maxN, pifile,
                ckpt_path ? ", checkpoint " : "", ckpt_path ? ckpt_path : "");
        fflush(outf);
    }

    uint64_t limit = 1;
    for (int i = 0; i < maxN; i++)
        limit *= 10;

    {
        uint64_t bytes = limit / 8 + 1, hi = 10;
        for (int n = 1; n <= maxN; n++, hi *= 10)
            bytes += hi / 8 + 1;
        fprintf(stderr, "estimated memory for bitmaps: %.1f GB\n",
                (double)bytes / 1e9);
    }

    fprintf(stderr, "sieving primes up to %" PRIu64 " (%d threads)...\n",
            limit, nthreads);
    build_sieve(limit + 1, nthreads);

    for (int n = 1; n <= maxN; n++) {
        struct search *s = &S[n - 1];
        memset(s, 0, sizeof *s);
        s->n = n;
        s->p10n1 = 1;
        for (int i = 1; i < n; i++)
            s->p10n1 *= 10;
        s->low = (n == 1) ? 1 : s->p10n1;
        uint64_t hi = s->p10n1 * 10;
        s->nprimes = count_primes(s->low, hi);
        s->remaining = s->nprimes;
        s->seen_bytes = hi / 8 + 1;
        /* +16: popcount_bytes reads the bitmap as 64-bit words. */
        s->seen = calloc(s->seen_bytes + 16, 1);
        if (!s->seen) {
            fprintf(stderr, "error: out of memory for n=%d\n", n);
            return 1;
        }
    }

    /* Resume from a checkpoint if one exists. */
    uint64_t i = 0; /* 0-based index into the fractional digit stream */
    int checked_prefix = 0;
    if (ckpt_path) {
        i = load_checkpoint(maxN);
        if (i > 0) {
            if (!skip_digits(f, i))
                return 1;
            checked_prefix = 1;
            /* Rebuild each rolling window from the trailing digits. */
            for (int n = 1; n <= maxN; n++) {
                struct search *s = &S[n - 1];
                uint64_t m = (i < (uint64_t)n) ? i : (uint64_t)n;
                s->v = 0;
                for (uint64_t k = i - m; k < i; k++)
                    s->v = s->v * 10 + ring[k & (RING - 1)];
            }
        }
    }

    nactive = 0;
    for (int n = 1; n <= maxN; n++) /* smallest n first */
        if (S[n - 1].remaining > 0)
            active[nactive++] = &S[n - 1];

    signal(SIGINT, on_signal);
    signal(SIGTERM, on_signal);

    static unsigned char buf[4 << 20];
    uint64_t next_report = (i / 1000000000ULL + 1) * 1000000000ULL;
    uint64_t next_ckpt = ckpt_interval ? i + ckpt_interval : UINT64_MAX;
    size_t got;

    while (nactive > 0 && i < maxdigits && !g_stop &&
           (got = fread(buf, 1, sizeof buf, f)) > 0) {
        for (size_t k = 0; k < got; k++) {
            unsigned d = (unsigned)buf[k] - '0';
            if (d > 9)
                continue; /* skip "3." dot, whitespace, etc. */

            /* Drop the integer-part 3 so position 1 is the first digit
             * after the decimal point (a raw fractional file starts 1415...,
             * so it has no leading 3 to drop). */
            if (!checked_prefix) {
                checked_prefix = 1;
                if (d == 3) {
                    g_dropped3 = 1;
                    continue;
                }
            }
            if (i < 8 && d != (unsigned)("14159265"[i] - '0'))
                fprintf(stderr, "warning: digit %" PRIu64 " is %u, but Pi "
                        "begins 3.14159265 -- is this really a Pi file?\n",
                        i + 1, d);

            for (int a = 0; a < nactive; a++) {
                struct search *s = active[a];
                if (i < (uint64_t)s->n) {
                    s->v = s->v * 10 + d;
                    if (i + 1 < (uint64_t)s->n)
                        continue;
                } else {
                    s->v = (s->v - ring[(i - s->n) & (RING - 1)] * s->p10n1)
                               * 10 + d;
                }
                uint64_t v = s->v;
                /*
                 * Cheap composite filters before touching the big bitmaps:
                 * an n>=2 prime must end in 1/3/7/9 and be coprime to 3 and
                 * 7. n=1 (primes 2,3,5,7) skips the filters.
                 */
                if (v < s->low)
                    continue;
                if (s->n > 1) {
                    if (!(d == 1 || d == 3 || d == 7 || d == 9))
                        continue;
                    if (v % 3 == 0 || v % 7 == 0)
                        continue;
                }
                __builtin_prefetch(&sieve_bits[v >> 3], 0, 0);
                __builtin_prefetch(&s->seen[v >> 3], 1, 0);
                q[qn].s = s;
                q[qn].v = v;
                q[qn].pos = i + 2 - (uint64_t)s->n;
                qn++;
            }
            ring[i & (RING - 1)] = (uint8_t)d;
            i++;

            if (qn > QCAP - MAXN) {
                flush_queue();
                if (nactive == 0)
                    break;
            }
            if (i >= next_report) {
                flush_queue();
                fprintf(stderr, "scanned %" PRIu64 " digits;", i);
                for (int a = 0; a < nactive; a++)
                    fprintf(stderr, " n=%d: %" PRIu64 " primes left",
                            active[a]->n, active[a]->remaining);
                fprintf(stderr, "\n");
                next_report += 1000000000ULL;
            }
            if (i >= next_ckpt && ckpt_path && nactive > 0) {
                flush_queue();
                save_checkpoint(maxN, i);
                next_ckpt = i + ckpt_interval;
            }
            if (nactive == 0 || i >= maxdigits)
                break;
        }
    }
    flush_queue();
    fclose(f);
    if (g_stop)
        fprintf(stderr, "interrupted at digit %" PRIu64 "\n", i);
    fprintf(stderr, "scanned %" PRIu64 " digits total\n", i);
    if (ckpt_path && nactive > 0)
        save_checkpoint(maxN, i);

    print_table(stdout, maxN);
    if (outf) {
        fprintf(outf, "[%s] run ended%s: scanned %" PRIu64 " digits total\n",
                timestamp(), g_stop ? " (interrupted)" : "", i);
        print_table(outf, maxN);
        fclose(outf);
    }

    for (int n = 1; n <= maxN; n++)
        free(S[n - 1].seen);
    free(sieve_bits);
    free(base_primes);
    return g_stop ? 130 : 0;
}
