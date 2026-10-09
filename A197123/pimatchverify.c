/*
 * pimatchverify - pass 3 of the minimizer/sort repeat finder (see pimatch.h).
 *
 * Reads the candidate groups from pimatchsort, collects the digits around
 * every candidate position (one sequential pass over pi, or direct seeks
 * with -S), and for every pair of positions sharing a k-mer extends the
 * match left and right to its full length.  Matches of at least the
 * requested length are reported with their maximal extent, sorted by second
 * occurrence, and the earliest second occurrence is called out for every
 * length the run covers - that is a(n) for each such n.
 *
 * Build: cc -O3 -march=native -pthread -o pimatchverify pimatchverify.c -lm
 */
#include "pimatch.h"            /* manifest, digit sources, formatting */
#include <pthread.h>            /* parallel seek workers */
#include <stdatomic.h>          /* batch counter shared by the seekers */

static struct {
    const char *src, *dir, *cand, *out;
    int Lmin, k, ctx, seek, seek_threads;
    double progress_secs;
    size_t chunk;
} opt = { NULL, NULL, NULL, NULL, 0, 0, 100, 0, 32, 30.0, 64u << 20 };

static void usage(const char *argv0)
{
    fprintf(stderr,
"pimatchverify v" PI_TOOLS_VERSION " - verify and extend candidate repeats (pass 3 of 3)\n"
"Usage: %s -d <scratchdir> [-p <source>] [options]\n"
"       %s -p <source> -c <candidates> -k <kmer> -L <minlen> [options]\n"
"\n"
"  -d DIR      scratch directory (supplies source, k, L and candidates.txt)\n"
"  -p SRC      pi digits (overrides the manifest's source)\n"
"  -c FILE     candidate file [DIR/candidates.txt]\n"
"  -k N        k-mer length of the candidates [from manifest]\n"
"  -L N        minimum match length to report [manifest L]; below the\n"
"              scan's guarantee the list is not exhaustive\n"
"  -E N        context digits kept either side of a k-mer [100]; matches\n"
"              longer than 2E+k are reported as truncated\n"
"  -S          seek to each position instead of streaming the whole source\n"
"              (ycd / pure-digit text; reads only the snippet bytes)\n"
"  -j N        parallel seek threads for -S [32]; on a disk array use enough\n"
"              to keep every spindle busy (32-128), on NFS 64+\n"
"  -o FILE     results file [DIR/matches.txt or matches.txt]\n"
"  -I SECS     progress interval [30]\n"
"  -V          print version and exit\n"
"  -h          this help\n", argv0, argv0);
    exit(1);
}

static void parse_args(int argc, char **argv)
{
    int c;
    while ((c = getopt(argc, argv, "d:p:c:k:L:E:Sj:o:I:Vh")) != -1) {
        switch (c) {
        case 'd': opt.dir = primary_dir(optarg); break;
        case 'p': opt.src = optarg; break;
        case 'c': opt.cand = optarg; break;
        case 'k': opt.k = atoi(optarg); break;
        case 'L': opt.Lmin = atoi(optarg); break;
        case 'E': opt.ctx = atoi(optarg); if (opt.ctx < 1 || opt.ctx > 100000) die("bad -E"); break;
        case 'S': opt.seek = 1; break;
        case 'j': opt.seek_threads = atoi(optarg); if (opt.seek_threads < 1 || opt.seek_threads > 512) die("-j must be 1..512"); break;
        case 'o': opt.out = optarg; break;
        case 'I': opt.progress_secs = atof(optarg); if (opt.progress_secs <= 0) die("bad -I"); break;
        case 'V': print_version("pimatchverify"); exit(0);
        default: usage(argv[0]);
        }
    }
    if (!opt.dir && (!opt.src || !opt.cand || !opt.k || !opt.Lmin)) { fprintf(stderr, "need -d, or all of -p -c -k -L\n"); usage(argv[0]); }
}

/* ---------------- candidates ---------------- */
typedef struct { uint64_t key; uint32_t first, n; } group;   /* positions gpos[first .. first+n) */
static group *groups; static size_t ngroups, gcap;
static uint64_t *gpos; static size_t ngpos, pcap;

static void load_candidates(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) die("cannot open candidate file '%s': %s", path, strerror(errno));
    char *line = NULL; size_t lcap = 0; ssize_t len;
    while ((len = getline(&line, &lcap, f)) > 0) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char *p = line, *e;
        uint64_t key = strtoull(p, &e, 10);
        if (e == p) die("bad candidate line: %s", line);
        if (ngroups == gcap) { gcap = gcap ? gcap * 2 : 4096; groups = realloc(groups, gcap * sizeof *groups); if (!groups) die("out of memory"); }
        group *g = &groups[ngroups];
        g->key = key; g->first = (uint32_t)ngpos; g->n = 0;
        p = e;
        for (;;) {
            uint64_t v = strtoull(p, &e, 10);
            if (e == p) break;
            if (ngpos == pcap) { pcap = pcap ? pcap * 2 : 8192; gpos = realloc(gpos, pcap * sizeof *gpos); if (!gpos) die("out of memory"); }
            gpos[ngpos++] = v; g->n++;
            p = e;
        }
        if (g->n >= 2) ngroups++; else ngpos = g->first;
    }
    free(line);
    fclose(f);
}

/* ---------------- snippets ---------------- */
static uint64_t *upos; static size_t nupos;      /* unique candidate positions, sorted */
static uint8_t *snip;                            /* nupos x snip_len digits, 0xFF = unavailable */
static int snip_len;                             /* 2E + k */

static int cmp_u64(const void *x, const void *y) { uint64_t a = *(const uint64_t *)x, b = *(const uint64_t *)y; return a < b ? -1 : a > b; }

static size_t upos_index(uint64_t p)
{
    size_t lo = 0, hi = nupos;
    while (lo < hi) { size_t mid = (lo + hi) / 2; if (upos[mid] < p) lo = mid + 1; else hi = mid; }
    return lo;
}

/* the snippet of position p covers digits [p - E, p - E + snip_len); digit d
 * of the snippet is at position p - E + d */
static inline uint64_t snip_start(uint64_t p) { return p > (uint64_t)opt.ctx ? p - (uint64_t)opt.ctx : 1; }

/*
 * collect_streaming - one sequential pass over the source; every window of
 * snip_len digits whose start is a wanted snippet start is copied out.
 */
static void collect_streaming(const char *src)
{
    pireader rd;
    pireader_open(&rd, src, snip_len, opt.chunk, UINT64_MAX);
    char desc[256]; pistream_describe(&rd.ps, desc, sizeof desc);
    printf("  source format  : %s\n", desc);
    size_t next = 0;                                  /* next snippet to fill */
    double t0 = now_sec(), tl = t0;
    int tty = isatty(1), shown = 0;
    size_t nwin;
    while (next < nupos && (nwin = pireader_fill(&rd)) > 0) {
        uint64_t base = rd.base;
        while (next < nupos) {
            uint64_t s = snip_start(upos[next]);
            if (s < base) { next++; continue; }        /* cannot happen: positions are sorted and we never go back */
            if (s >= base + nwin) break;
            memcpy(snip + next * (size_t)snip_len, rd.buf + (s - base), (size_t)snip_len);
            next++;
        }
        double t = now_sec();
        if (t - tl >= opt.progress_secs) {
            char c1[32], c2[32], d1[32];
            double pct = rd.digits_est ? 100.0 * (double)(rd.base + nwin - rd.first_pos) / (double)rd.digits_est : 0;
            printf("%s[%s] pos %s (%.1f%%)  snippets %s/%zu%s", tty ? "\r\033[K" : "", fmt_duration(t - t0, d1, sizeof d1),
                   fmt_u64(rd.base + nwin - 1, c1, sizeof c1), pct, fmt_u64(next, c2, sizeof c2), nupos, tty ? "" : "\n");
            fflush(stdout); shown = tty; tl = t;
        }
        if (g_stop) break;
    }
    /* tail: positions whose full snippet does not fit before EOF get what exists */
    if (next < nupos && rd.ndigits) {
        uint64_t base = rd.base;
        for (; next < nupos; next++) {
            uint64_t s = snip_start(upos[next]);
            if (s < base || s >= base + rd.ndigits) continue;
            size_t avail = (size_t)(base + rd.ndigits - s);
            if (avail > (size_t)snip_len) avail = (size_t)snip_len;
            memcpy(snip + next * (size_t)snip_len, rd.buf + (s - base), avail);
        }
    }
    if (shown) printf("\n");
    pireader_close(&rd);
}

/* ---------------- seek mode: direct parallel preads ---------------- */
/*
 * Seek mode reads only the bytes of each snippet with pread(): no file
 * reopen, no readahead hint, no sequential prefetcher (each of those cost
 * megabytes of disk traffic per 200-digit snippet in the pistream path and
 * capped a 21 TB text source at ~73 snippets/s).  -j threads keep many
 * reads in flight so a disk array delivers its aggregate IOPS; the sorted
 * position list is handed out in small batches, so the outstanding reads
 * are always close together on disk (elevator-friendly).
 */
#define SEEK_BATCH 64
static pistream *seek_ps;                        /* file table of the source */
static _Atomic size_t seek_next, seek_done;

typedef struct {
    int fd, cur;                                 /* open file, index in seek_ps->files */
    uint8_t *raw;                                /* pread buffer */
    uint8_t *dec;                                /* decoded ycd words */
    size_t rawcap;
} seeker;

/*
 * seek_find_file - index of the source file holding 0-based digit offset o
 * (the seeker's current file is tried first), or -1 past the end.
 */
static int seek_find_file(const seeker *sk, uint64_t o)
{
    const pistream *ps = seek_ps;
    if (sk->cur >= 0) {
        const pifile *f = &ps->files[sk->cur];
        if (o >= f->start && o < f->start + f->ndigits) return sk->cur;
    }
    int lo = 0, hi = ps->nfiles;
    while (lo < hi) { int mid = (lo + hi) / 2; if (ps->files[mid].start + ps->files[mid].ndigits <= o) lo = mid + 1; else hi = mid; }
    if (lo < ps->nfiles && o >= ps->files[lo].start) return lo;
    return -1;
}

/* seeker_use - make file idx the seeker's open file (random-access hint) */
static void seeker_use(seeker *sk, int idx)
{
    if (sk->fd >= 0 && sk->cur == idx) return;
    if (sk->fd >= 0) close(sk->fd);
    sk->fd = open(seek_ps->files[idx].path, O_RDONLY);
    if (sk->fd < 0) die("cannot open pi source '%s': %s", seek_ps->files[idx].path, strerror(errno));
#ifdef POSIX_FADV_RANDOM
    posix_fadvise(sk->fd, 0, 0, POSIX_FADV_RANDOM);
#endif
    sk->cur = idx;
}

/* pread_full - read up to n bytes at off; returns bytes read (short at EOF) */
static size_t pread_full(int fd, uint8_t *buf, size_t n, uint64_t off)
{
    size_t got = 0;
    while (got < n) {
        ssize_t r = pread(fd, buf + got, n - got, (off_t)(off + got));
        if (r < 0) { if (errno == EINTR) continue; die("pread failed: %s", strerror(errno)); }
        if (r == 0) break;
        got += (size_t)r;
    }
    return got;
}

/*
 * seek_read - fill dst[0..n) with the digit values starting at 0-based
 * stream offset o (text: byte offset after "3."; ycd: word offset space
 * including the lead digit when the block set carries it), spanning file
 * boundaries.  Digits past the end of the data are left as 0xFF.
 */
static void seek_read(seeker *sk, uint64_t o, uint8_t *dst, size_t n)
{
    while (n > 0) {
        int idx = seek_find_file(sk, o);
        if (idx < 0) return;
        seeker_use(sk, idx);
        const pifile *f = &seek_ps->files[idx];
        uint64_t local = o - f->start;
        size_t take = f->ndigits - local < n ? (size_t)(f->ndigits - local) : n;
        if (f->kind == PIF_TEXT) {
            size_t got = pread_full(sk->fd, sk->raw, take, f->data_off + local);
            for (size_t i = 0; i < got; i++) { unsigned c = (unsigned)sk->raw[i] - '0'; dst[i] = c < 10u ? (uint8_t)c : 0xFF; }
            if (got < take) return;
        } else {
            uint64_t w0 = local / YCD_DIGITS_PER_WORD;
            size_t drop = (size_t)(local % YCD_DIGITS_PER_WORD);
            size_t nw = (drop + take + YCD_DIGITS_PER_WORD - 1) / YCD_DIGITS_PER_WORD;
            size_t got = pread_full(sk->fd, sk->raw, nw * 8, f->data_off + w0 * 8) / 8;
            for (size_t w = 0; w < got; w++) ycd_decode_word(ycd_load_le64(sk->raw + w * 8), sk->dec + w * YCD_DIGITS_PER_WORD);
            size_t avail = got * YCD_DIGITS_PER_WORD > drop ? got * YCD_DIGITS_PER_WORD - drop : 0;
            if (avail > take) avail = take;
            memcpy(dst, sk->dec + drop, avail);
            if (avail < take) return;
        }
        dst += take; o += take; n -= take;
    }
}

/* seek_worker - thread body: claim batches of sorted positions, read snippets */
static void *seek_worker(void *arg)
{
    (void)arg;
    seeker sk = { -1, -1, NULL, NULL, 0 };
    sk.rawcap = (size_t)snip_len + 2 * YCD_DIGITS_PER_WORD + 64;
    sk.raw = malloc(sk.rawcap * 8 / YCD_DIGITS_PER_WORD + sk.rawcap);   /* enough for either format */
    sk.dec = malloc(sk.rawcap + 2 * YCD_DIGITS_PER_WORD);
    if (!sk.raw || !sk.dec) die("out of memory");
    uint64_t first = seek_ps->first_pos;
    for (;;) {
        size_t b = atomic_fetch_add(&seek_next, (size_t)SEEK_BATCH);
        if (b >= nupos || g_stop) break;
        size_t e = b + SEEK_BATCH < nupos ? b + SEEK_BATCH : nupos;
        for (size_t i = b; i < e; i++) {
            uint64_t s = snip_start(upos[i]);
            uint8_t *dst = snip + i * (size_t)snip_len;
            size_t n = (size_t)snip_len;
            if (s < first) {                              /* source starts after position 1: leave the head unavailable */
                size_t gap = (size_t)(first - s);
                if (gap >= n) continue;
                dst += gap; n -= gap; s = first;
            }
            uint64_t o = seek_ps->kind == PIF_YCD ? (s - 1) + (uint64_t)seek_ps->lead3 : s - 1;
            seek_read(&sk, o, dst, n);
        }
        atomic_fetch_add(&seek_done, e - b);
    }
    if (sk.fd >= 0) close(sk.fd);
    free(sk.raw); free(sk.dec);
    return NULL;
}

/*
 * collect_seek - read every snippet directly with -j parallel seekers
 * (ycd or pure-digit text sources only).
 */
static void collect_seek(const char *src)
{
    pistream ps;
    setenv("PI_NO_PREFETCH", "1", 1);               /* no sequential prefetcher for random access */
    pistream_open(&ps, src, 1u << 20);
    char desc[256]; pistream_describe(&ps, desc, sizeof desc);
    printf("  source format  : %s\n", desc);
    if (ps.kind == PIF_TEXT && !pistream_text_pure(&ps)) die("this source cannot be seeked (text with non-digit bytes) - drop -S");
    seek_ps = &ps;
    atomic_store(&seek_next, 0); atomic_store(&seek_done, 0);
    pthread_t *th = calloc((size_t)opt.seek_threads, sizeof *th);
    if (!th) die("out of memory");
    for (int t = 0; t < opt.seek_threads; t++)
        if (pthread_create(&th[t], NULL, seek_worker, NULL)) die("cannot create seek thread");
    double t0 = now_sec(), tl = t0;
    int tty = isatty(1), shown = 0;
    size_t last_done = 0;
    for (;;) {
        struct timespec ts = { 0, 100 * 1000 * 1000 };
        nanosleep(&ts, NULL);
        size_t done = atomic_load(&seek_done);
        double t = now_sec();
        if (done >= nupos || g_stop) break;
        if (t - tl >= opt.progress_secs) {
            char d1[32], d2[32], c1[32], c2[32];
            double rate = (double)(done - last_done) / (t - tl);
            double eta = rate > 0 ? (double)(nupos - done) / rate : 0;
            printf("%s[%s] snippets %s/%s (%.0f/s, %d threads)  ETA %s%s", tty ? "\r\033[K" : "", fmt_duration(t - t0, d1, sizeof d1),
                   fmt_u64(done, c1, sizeof c1), fmt_u64(nupos, c2, sizeof c2), rate, opt.seek_threads,
                   fmt_duration(eta, d2, sizeof d2), tty ? "" : "\n");
            fflush(stdout); shown = tty; tl = t; last_done = done;
        }
    }
    for (int t = 0; t < opt.seek_threads; t++) pthread_join(th[t], NULL);
    free(th);
    if (shown) printf("\n");
    pistream_close(&ps);
}

/* ---------------- matches ---------------- */
typedef struct { uint64_t p1, p2; uint32_t len; uint8_t truncated; } match;
static match *matches; static size_t nmatches, mcap;

static int cmp_match_second(const void *a, const void *b)
{
    const match *x = a, *y = b;
    if (x->p2 != y->p2) return x->p2 < y->p2 ? -1 : 1;
    if (x->p1 != y->p1) return x->p1 < y->p1 ? -1 : 1;
    return x->len < y->len ? -1 : x->len > y->len;
}
static int cmp_match_key(const void *a, const void *b)
{
    const match *x = a, *y = b;
    if (x->p1 != y->p1) return x->p1 < y->p1 ? -1 : 1;
    if (x->p2 != y->p2) return x->p2 < y->p2 ? -1 : 1;
    return x->len < y->len ? -1 : x->len > y->len;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);
    setvbuf(stdout, NULL, _IOLBF, 0);
    install_signals();

    manifest m; memset(&m, 0, sizeof m);
    char candname[2048], outname[2048];
    if (opt.dir) {
        manifest_read(opt.dir, &m);
        if (!opt.src) opt.src = m.source;
        if (!opt.k) opt.k = m.k;
        if (!opt.Lmin) opt.Lmin = m.L;
        if (!opt.cand) { snprintf(candname, sizeof candname, "%s/candidates.txt", opt.dir); opt.cand = candname; }
        if (!opt.out) { snprintf(outname, sizeof outname, "%s/matches.txt", opt.dir); opt.out = outname; }
    }
    if (!opt.out) opt.out = "matches.txt";
    if (opt.Lmin < opt.k) die("-L (%d) must be at least k (%d)", opt.Lmin, opt.k);
    snip_len = 2 * opt.ctx + opt.k;

    char b1[32], b2[32], b3[32];
    printf("pimatchverify v%s - verify and extend candidate repeats (pass 3 of 3)\n", PI_TOOLS_VERSION);
    printf("  source         : %s\n  candidates     : %s\n  k-mer / min len: %d / %d\n", opt.src, opt.cand, opt.k, opt.Lmin);
    if (opt.dir && m.L && opt.Lmin < m.L)
        warn("the scan guarantees repeats of >= %d digits; the list for lengths below that is not exhaustive", m.L);
    if (opt.dir && !m.complete) warn("the scan in '%s' did not complete; results cover only the scanned range", opt.dir);
    if (opt.dir && m.hpasses > 1)
        printf("  hash range     : pass %d of %d - these results cover one slice of the key space; the answer is the\n"
               "                   earliest second occurrence across all %d passes' matches.txt files\n", m.hpass, m.hpasses, m.hpasses);

    load_candidates(opt.cand);
    uint64_t npairs = 0; for (size_t g = 0; g < ngroups; g++) npairs += (uint64_t)groups[g].n * (groups[g].n - 1) / 2;
    printf("  candidate groups: %s, %s positions, %s pairs\n", fmt_u64(ngroups, b1, sizeof b1), fmt_u64(ngpos, b2, sizeof b2), fmt_u64(npairs, b3, sizeof b3));
    if (ngroups == 0) { printf("Nothing to verify.\n"); return 0; }

    /* unique sorted positions and snippet storage */
    upos = malloc(ngpos * sizeof *upos); if (!upos) die("out of memory");
    memcpy(upos, gpos, ngpos * sizeof *upos);
    qsort(upos, ngpos, sizeof *upos, cmp_u64);
    nupos = 0;
    for (size_t i = 0; i < ngpos; i++) if (nupos == 0 || upos[i] != upos[nupos - 1]) upos[nupos++] = upos[i];
    snip = malloc(nupos * (size_t)snip_len);
    if (!snip) die("cannot allocate %s of snippets", fmt_bytes((double)nupos * snip_len, b1, sizeof b1));
    memset(snip, 0xFF, nupos * (size_t)snip_len);
    printf("  snippets       : %s x %d digits = %s (%s mode)\n", fmt_u64(nupos, b1, sizeof b1), snip_len,
           fmt_bytes((double)nupos * snip_len, b2, sizeof b2), opt.seek ? "seek" : "streaming");

    double t0 = now_sec();
    if (opt.seek) collect_seek(opt.src); else collect_streaming(opt.src);
    double t_collect = now_sec() - t0;
    if (g_stop) { printf("*** interrupted ***\n"); return 130; }

    /* extend every pair */
    uint64_t checked = 0, missing = 0, kmismatch = 0, tooshort = 0, truncated = 0;
    for (size_t g = 0; g < ngroups; g++) {
        const group *G = &groups[g];
        for (uint32_t a = 0; a < G->n; a++) for (uint32_t b = a + 1; b < G->n; b++) {
            uint64_t pa = gpos[G->first + a], pb = gpos[G->first + b];
            if (pa == pb) continue;
            checked++;
            size_t ia = upos_index(pa), ib = upos_index(pb);
            const uint8_t *sa = snip + ia * (size_t)snip_len, *sb = snip + ib * (size_t)snip_len;
            int oa = (int)(pa - snip_start(pa)), ob = (int)(pb - snip_start(pb));   /* k-mer offset in snippet */
            /* k-mer must be present and equal */
            int ok = 1;
            for (int i = 0; i < opt.k; i++) if (sa[oa + i] > 9 || sb[ob + i] > 9 || sa[oa + i] != sb[ob + i]) { ok = 0; break; }
            if (!ok) { if (sa[oa] > 9 || sb[ob] > 9) missing++; else kmismatch++; continue; }
            int left = 0, right = 0, trunc = 0;
            while (oa - left - 1 >= 0 && ob - left - 1 >= 0 && sa[oa - left - 1] <= 9 && sa[oa - left - 1] == sb[ob - left - 1]) left++;
            if (oa - left == 0 || ob - left == 0) trunc = 1;
            while (oa + opt.k + right < snip_len && ob + opt.k + right < snip_len &&
                   sa[oa + opt.k + right] <= 9 && sa[oa + opt.k + right] == sb[ob + opt.k + right]) right++;
            if (oa + opt.k + right >= snip_len || ob + opt.k + right >= snip_len) trunc = 1;
            uint32_t len = (uint32_t)(left + opt.k + right);
            if (len < (uint32_t)opt.Lmin) { tooshort++; continue; }
            if (nmatches == mcap) { mcap = mcap ? mcap * 2 : 1024; matches = realloc(matches, mcap * sizeof *matches); if (!matches) die("out of memory"); }
            match *mt = &matches[nmatches++];
            mt->p1 = pa - (uint64_t)left; mt->p2 = pb - (uint64_t)left; mt->len = len; mt->truncated = (uint8_t)trunc;
            if (mt->p1 > mt->p2) { uint64_t t = mt->p1; mt->p1 = mt->p2; mt->p2 = t; }
        }
    }
    /* dedupe (the same repeat is found once per shared minimizer) */
    qsort(matches, nmatches, sizeof *matches, cmp_match_key);
    size_t u = 0;
    for (size_t i = 0; i < nmatches; i++)
        if (u == 0 || matches[i].p1 != matches[u - 1].p1 || matches[i].p2 != matches[u - 1].p2 || matches[i].len != matches[u - 1].len)
            matches[u++] = matches[i];
    nmatches = u;
    for (size_t i = 0; i < nmatches; i++) truncated += matches[i].truncated;
    qsort(matches, nmatches, sizeof *matches, cmp_match_second);

    /* results file */
    FILE *out = fopen(opt.out, "w");
    if (!out) die("cannot create '%s': %s", opt.out, strerror(errno));
    fprintf(out, "# pimatchverify v%s source=%s candidates=%s k=%d minlen=%d  format: len pos1 pos2 digits[T=truncated]\n",
            PI_TOOLS_VERSION, opt.src, opt.cand, opt.k, opt.Lmin);
    for (size_t i = 0; i < nmatches; i++) {
        const match *mt = &matches[i];
        size_t ia = upos_index(mt->p1 <= mt->p2 ? mt->p1 : mt->p2);
        /* the digits: the match starts at p1, which may precede the snippet's own position; find a snippet containing it */
        (void)ia;
        fprintf(out, "%u %" PRIu64 " %" PRIu64 " ", mt->len, mt->p1, mt->p2);
        /* locate a snippet whose window covers [p1, p1+len): search unique positions near p1 */
        size_t j = upos_index(mt->p1);
        int printed = 0;
        for (size_t q = (j > 4 ? j - 4 : 0); q < nupos && q < j + 8 && !printed; q++) {
            uint64_t s = snip_start(upos[q]);
            if (mt->p1 >= s && mt->p1 + mt->len <= s + (uint64_t)snip_len) {
                const uint8_t *sp = snip + q * (size_t)snip_len + (mt->p1 - s);
                for (uint32_t d = 0; d < mt->len; d++) fputc(sp[d] > 9 ? '?' : '0' + sp[d], out);
                printed = 1;
            }
        }
        if (!printed) fputs("(digits not in any snippet)", out);
        if (mt->truncated) fputs(" T", out);
        fputc('\n', out);
    }

    /* summary */
    char d1[32];
    printf("\nVerification complete\n");
    printf("  pairs checked   : %s  (%s no k-mer match, %s shorter than %d, %s unavailable)\n", fmt_u64(checked, b1, sizeof b1),
           fmt_u64(kmismatch, b2, sizeof b2), fmt_u64(tooshort, b3, sizeof b3), opt.Lmin, fmt_u64(missing, d1, sizeof d1));
    printf("  repeats >= %d   : %s distinct (%s truncated at the %d-digit context - raise -E to see their full length)\n",
           opt.Lmin, fmt_u64(nmatches, b1, sizeof b1), fmt_u64(truncated, b2, sizeof b2), snip_len);
    printf("  collect time    : %s (%s)\n", fmt_duration(t_collect, d1, sizeof d1), opt.seek ? "seek" : "streaming");
    if (nmatches) {
        uint32_t maxlen = 0;
        for (size_t i = 0; i < nmatches; i++) if (matches[i].len > maxlen) maxlen = matches[i].len;
        printf("\nEarliest second occurrence per length (the a(n) candidates this run establishes):\n");
        for (uint32_t n = (uint32_t)opt.Lmin; n <= maxlen && n <= (uint32_t)opt.Lmin + 20; n++) {
            const match *best = NULL;
            for (size_t i = 0; i < nmatches; i++) if (matches[i].len >= n) { best = &matches[i]; break; }   /* sorted by p2 */
            if (!best) break;
            /* digits of the n-digit prefix */
            size_t j = upos_index(best->p1); char digs[64] = "?";
            for (size_t q = (j > 4 ? j - 4 : 0); q < nupos && q < j + 8; q++) {
                uint64_t s = snip_start(upos[q]);
                if (best->p1 >= s && best->p1 + n <= s + (uint64_t)snip_len) {
                    const uint8_t *sp = snip + q * (size_t)snip_len + (best->p1 - s);
                    size_t dn = n < 60 ? n : 60;
                    for (size_t d = 0; d < dn; d++) digs[d] = sp[d] > 9 ? '?' : (char)('0' + sp[d]);
                    digs[dn] = 0; break;
                }
            }
            printf("  n=%-3u %s%s  at %s and %s%s\n", n, digs, n > 60 ? "..." : "", fmt_u64(best->p1, b1, sizeof b1), fmt_u64(best->p2, b2, sizeof b2),
                   (opt.dir && (int)n < m.L) ? "  (below the scan guarantee: not exhaustive)" : "");
        }
        fprintf(out, "# summary: %zu repeats >= %d digits; earliest second occurrence: len %u at %" PRIu64 " and %" PRIu64 "\n",
                nmatches, opt.Lmin, matches[0].len, matches[0].p1, matches[0].p2);
    } else printf("\nNo repeats of >= %d digits among the candidates.\n", opt.Lmin);
    if (fclose(out)) die("write error on '%s'", opt.out);
    printf("\nResults: %s\n", opt.out);
    return 0;
}
