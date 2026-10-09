/*
 * pimatchscan - pass 1 of the minimizer/sort repeat finder (see pimatch.h).
 *
 * Streams the digits of pi once.  Every window of w consecutive k-mers
 * (i.e. every L = k+w-1 digit window) selects its minimizer - the k-mer
 * with the smallest random rank - and a (k-mer value, position) record is
 * written for each distinct minimizer into the partition file chosen by a
 * hash of the value.  Identical L-digit strings anywhere in pi therefore
 * produce records with equal keys in the same partition.
 *
 * Pipeline: one reader thread (source -> chunk ring), N worker threads each
 * taking whole chunks (sliding-minimum over k-mer ranks, O(1) per digit),
 * per-worker per-partition buffers, and block-buffered partition writers
 * (large writes, no fsync except at checkpoints).  Checkpoints record the
 * resume position and every partition file's length; --resume truncates and
 * continues.  Duplicate records around chunk/resume boundaries are dropped
 * by pimatchsort.
 *
 * Build: cc -O3 -march=native -pthread -o pimatchscan pimatchscan.c -lm
 */
#include "pimatch.h"            /* records, manifest, checkpoint, digit sources */
#include <pthread.h>            /* reader, workers, ring synchronisation */
#include <stdatomic.h>          /* progress counters */

#define MAX_WORKERS 128
#define WBUF_RECS 2048          /* per-worker per-partition buffer (28 KB) */

static struct {
    const char *src, *dirspec;
    const char *dir;            /* primary scratch directory (manifest, checkpoint) */
    uint64_t ndigits;
    int L, k, w;
    int threads, nparts;
    int hpass, hpasses;         /* -H pass/npasses hash range [1/1] */
    uint64_t block_bytes;       /* partition write block */
    double ckpt_secs, progress_secs;
    int resume, dryrun, force;
    size_t chunk;
} opt = { NULL, NULL, NULL, UINT64_MAX, 0, KMER_MAX_EXACT, 0, 0, 64, 1, 1, 4u << 20, 600.0, 30.0, 0, 0, 0, 64u << 20 };
static scratch SC;              /* scratch directories and partition ranges */

static winspec W;               /* only for L bookkeeping */
static uint64_t POW10K;         /* 10^(k-1) */

/* ---------------- partition writers ---------------- */
typedef struct {
    pthread_mutex_t m;
    int fd;
    uint8_t *buf;
    size_t used, cap;
    uint64_t bytes;             /* bytes written to the file (incl. checkpointed) */
    uint64_t records;
} pwriter;
static pwriter *pw;

/*
 * pw_flush_locked - write the partition's block buffer out (caller holds m).
 */
static void pw_flush_locked(pwriter *p, int idx)
{
    size_t off = 0;
    while (off < p->used) {
        ssize_t n = write(p->fd, p->buf + off, p->used - off);
        if (n < 0) { if (errno == EINTR) continue; die("write error on partition %d: %s", idx, strerror(errno)); }
        off += (size_t)n;
    }
    p->bytes += p->used;
    p->used = 0;
}

/*
 * pw_append - append n records to partition idx (thread-safe, buffered).
 */
static void pw_append(int idx, const mrec *recs, size_t n)
{
    pwriter *p = &pw[idx];
    size_t bytes = n * MREC_BYTES;
    pthread_mutex_lock(&p->m);
    p->records += n;
    if (!opt.dryrun) {
        if (p->used + bytes > p->cap) pw_flush_locked(p, idx);
        if (bytes > p->cap) {                       /* larger than the block: write directly */
            size_t off = 0;
            while (off < bytes) {
                ssize_t r = write(p->fd, (const uint8_t *)recs + off, bytes - off);
                if (r < 0) { if (errno == EINTR) continue; die("write error on partition %d: %s", idx, strerror(errno)); }
                off += (size_t)r;
            }
            p->bytes += bytes;
        } else {
            memcpy(p->buf + p->used, recs, bytes);
            p->used += bytes;
        }
    }
    pthread_mutex_unlock(&p->m);
}

/* ---------------- chunk ring ---------------- */
typedef struct {
    uint8_t *buf;
    size_t nwin;
    uint64_t base, seq;
    int state;                  /* 0 free, 1 published, 2 taken */
} slot_t;
static slot_t *slots;
static int nslots;
static pthread_mutex_t ring_mtx = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t ring_can_fill = PTHREAD_COND_INITIALIZER, ring_can_take = PTHREAD_COND_INITIALIZER;
static uint64_t chunks_published, chunks_claimed;
static int producing_done;
static uint64_t g_windows, g_last_pos;
static uint64_t *chunk_end;     /* chunk_end[seq % ring] = base + nwin of that chunk */
static uint8_t  *chunk_done;    /* chunk_done[seq % ring] */
static uint64_t done_contig;    /* chunks 0..done_contig-1 fully processed */
static uint64_t resume_pos;     /* position where this run started */

/* ---------------- statistics ---------------- */
typedef struct { _Atomic uint64_t records, windows; char pad[48]; } wstat_t;
static wstat_t wstats[MAX_WORKERS];
static uint64_t base_records;   /* records from before a resume */

static void usage(const char *argv0)
{
    fprintf(stderr,
"pimatchscan v" PI_TOOLS_VERSION " - minimizer scan of pi into partition files (pass 1 of 3)\n"
"Usage: %s -p <source> -d <scratchdir> -L <len> [options]\n"
"\n"
"  -p SRC      pi digits: text or .ycd file(s)/directories, comma-separated\n"
"  -d DIRS     scratch directory (must exist), or a comma-separated list: the\n"
"              partitions are spread over the list in proportion to free space\n"
"              and the first directory holds the manifest and candidates\n"
"  -H K/R      hash-range pass K of R [1/1]: write only the keys whose hash\n"
"              falls in slice K, so a scan needing more scratch than the set\n"
"              holds can be run as R independent passes (sort and verify\n"
"              each; the answer is the earliest across the R results)\n"
"  -L N        guaranteed match length: every repeat of >= N digits is found\n"
"  -k N        k-mer length, 1..19 (exact 64-bit keys) [19]\n"
"  -w N        minimizer window (k-mers per window) [derived: L-k+1]\n"
"  -n COUNT    digits of pi to scan (500M, 10B, 1T, all) [all]\n"
"  -P N        number of partition files, 1..%d [64]; size them so one\n"
"              partition is <= 40%% of the sorting machine's RAM\n"
"  -t N        worker threads [min(8, CPUs)]\n"
"  -B SIZE     partition write block [4M]; writes are buffered, never synced\n"
"              except at checkpoints\n"
"  -C SECS     checkpoint interval [600]; 0 disables\n"
"  -R          resume from the checkpoint in DIR\n"
"  -f          allow a non-empty DIR without -R (files are overwritten)\n"
"  -T          dry run: scan and count records, write nothing\n"
"  -I SECS     progress interval [30]\n"
"  -c SIZE     read chunk size [64M]\n"
"  -V          print version and exit\n"
"  -h          this help\n", argv0, MAX_PARTS);
    exit(1);
}

static void parse_args(int argc, char **argv)
{
    int c;
    long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    opt.threads = ncpu >= 8 ? 8 : (int)(ncpu > 0 ? ncpu : 1);
    while ((c = getopt(argc, argv, "p:d:L:k:w:n:P:t:B:C:H:RfTI:c:Vh")) != -1) {
        switch (c) {
        case 'p': opt.src = optarg; break;
        case 'd': opt.dirspec = optarg; break;
        case 'H':
            if (sscanf(optarg, "%d/%d", &opt.hpass, &opt.hpasses) != 2 || opt.hpass < 1 || opt.hpasses < 1 || opt.hpass > opt.hpasses || opt.hpasses > 1024)
                die("-H must be K/R with 1 <= K <= R <= 1024");
            break;
        case 'L': opt.L = atoi(optarg); break;
        case 'k': opt.k = atoi(optarg); break;
        case 'w': opt.w = atoi(optarg); break;
        case 'n': if (parse_count(optarg, &opt.ndigits)) die("bad digit count '%s'", optarg); break;
        case 'P': opt.nparts = atoi(optarg); break;
        case 't': opt.threads = atoi(optarg); break;
        case 'B': if (parse_mem_size(optarg, &opt.block_bytes) || opt.block_bytes < MREC_BYTES * 64) die("bad -B"); break;
        case 'C': opt.ckpt_secs = atof(optarg); break;
        case 'R': opt.resume = 1; break;
        case 'f': opt.force = 1; break;
        case 'T': opt.dryrun = 1; break;
        case 'I': opt.progress_secs = atof(optarg); if (opt.progress_secs <= 0) die("bad -I"); break;
        case 'c': { uint64_t v; if (parse_mem_size(optarg, &v) || v < (1u << 20)) die("bad chunk size"); opt.chunk = (size_t)v; } break;
        case 'V': print_version("pimatchscan"); exit(0);
        default: usage(argv[0]);
        }
    }
    if (!opt.src || !opt.dirspec) { fprintf(stderr, "missing -p or -d\n"); usage(argv[0]); }
    if (opt.k < 1 || opt.k > KMER_MAX_EXACT) die("-k must be 1..%d", KMER_MAX_EXACT);
    if (opt.L == 0 && opt.w == 0) die("give -L (guaranteed match length) or -w");
    if (opt.w == 0) opt.w = opt.L - opt.k + 1;
    if (opt.L == 0) opt.L = opt.k + opt.w - 1;
    if (opt.w < 1) die("-L %d is shorter than -k %d: use a smaller k", opt.L, opt.k);
    if (opt.L != opt.k + opt.w - 1) die("-L, -k and -w disagree: L must equal k + w - 1 (%d + %d - 1 = %d)", opt.k, opt.w, opt.k + opt.w - 1);
    if (opt.L > MAX_SEQ_LEN) die("-L must be <= %d", MAX_SEQ_LEN);
    if (opt.nparts < 1 || opt.nparts > MAX_PARTS) die("-P must be 1..%d", MAX_PARTS);
    if (opt.threads < 1 || opt.threads > MAX_WORKERS) die("-t must be 1..%d", MAX_WORKERS);
    if (opt.ndigits == 0) die("-n must be > 0");
}

/* ---------------- reader thread ---------------- */
static pireader rd;

/*
 * reader_main - fill free ring slots with chunks from the digit source and
 * publish them; workers claim published chunks in order.
 */
static void *reader_main(void *arg)
{
    (void)arg;
    uint64_t seq = 0;
    size_t nwin;
    while (!g_stop && (nwin = pireader_fill(&rd)) > 0) {
        slot_t *s = &slots[seq % (uint64_t)nslots];
        pthread_mutex_lock(&ring_mtx);
        /* wait for a free slot AND keep every in-flight chunk within nslots of
         * done_contig, so chunk_done/chunk_end[seq % nslots] are unambiguous */
        while (s->state != 0 || seq >= done_contig + (uint64_t)nslots) pthread_cond_wait(&ring_can_fill, &ring_mtx);
        pthread_mutex_unlock(&ring_mtx);
        memcpy(s->buf, rd.buf, rd.ndigits);
        s->nwin = nwin; s->base = rd.base; s->seq = seq;
        pthread_mutex_lock(&ring_mtx);
        s->state = 1;
        chunk_end[seq % (uint64_t)nslots] = rd.base + nwin;
        chunk_done[seq % (uint64_t)nslots] = 0;
        chunks_published = ++seq;
        g_windows += nwin;
        g_last_pos = rd.base + nwin - 1;
        pthread_cond_broadcast(&ring_can_take);
        pthread_mutex_unlock(&ring_mtx);
    }
    pthread_mutex_lock(&ring_mtx);
    producing_done = 1;
    pthread_cond_broadcast(&ring_can_take);
    pthread_mutex_unlock(&ring_mtx);
    return NULL;
}

/* ---------------- workers ---------------- */
typedef struct { int id; } warg_t;

/*
 * worker_main - claim chunks in sequence order; for each chunk roll the
 * k-mer value, compute its rank, keep a monotonic deque of candidate
 * minimizers over the last w k-mers, and emit one record per distinct
 * minimizer position into per-partition buffers (flushed to the shared
 * partition writers when full and at the end of every chunk, so a
 * checkpoint taken after the chunk holds all its records).
 */
static void *worker_main(void *argp)
{
    const int id = ((warg_t *)argp)->id;
    wstat_t *st = &wstats[id];
    const int k = opt.k, w = opt.w, np = opt.nparts, hp = opt.hpass, hn = opt.hpasses;

    mrec *wb = malloc((size_t)np * WBUF_RECS * sizeof(mrec));
    uint32_t *wn = calloc((size_t)np, sizeof *wn);
    uint64_t *rank = NULL; uint64_t *val = NULL; size_t rank_cap = 0;
    int *dq = NULL;                     /* deque of k-mer indices */
    if (!wb || !wn) die("out of memory");

    uint64_t local_recs = 0, local_win = 0;
    for (;;) {
        pthread_mutex_lock(&ring_mtx);
        /* claim the next unclaimed chunk atomically AFTER waiting: reading the
         * counter before the wait let several woken workers claim one chunk */
        while (chunks_claimed >= chunks_published && !producing_done) pthread_cond_wait(&ring_can_take, &ring_mtx);
        if (chunks_claimed >= chunks_published && producing_done) { pthread_mutex_unlock(&ring_mtx); break; }
        uint64_t seq = chunks_claimed++;
        slot_t *s = &slots[seq % (uint64_t)nslots];
        s->state = 2;
        pthread_mutex_unlock(&ring_mtx);

        const uint8_t *d = s->buf;
        const size_t nwin = s->nwin;
        const uint64_t base = s->base;
        const size_t nk = nwin + (size_t)w - 1;          /* k-mers in this chunk */
        if (nk > rank_cap) {
            rank_cap = nk + 1024;
            rank = realloc(rank, rank_cap * sizeof *rank);
            val  = realloc(val,  rank_cap * sizeof *val);
            dq   = realloc(dq,   rank_cap * sizeof *dq);
            if (!rank || !val || !dq) die("out of memory");
        }
        /* k-mer values and ranks */
        uint64_t v = 0;
        for (int i = 0; i < k; i++) v = v * 10 + d[i];
        val[0] = v; rank[0] = kmer_rank(v);
        for (size_t i = 1; i < nk; i++) {
            v = (v - (uint64_t)d[i - 1] * POW10K) * 10 + d[i + (size_t)k - 1];
            val[i] = v; rank[i] = kmer_rank(v);
        }
        /* sliding minimum over windows of w k-mers; ties -> leftmost */
        int head = 0, tail = 0;             /* dq[head..tail) */
        long last_emitted = -1;
        for (size_t i = 0; i < nk; i++) {
            while (tail > head && rank[dq[tail - 1]] > rank[i]) tail--;
            dq[tail++] = (int)i;
            if (i + 1 >= (size_t)w) {                        /* window ending at i is complete */
                size_t j = i + 1 - (size_t)w;                /* window start */
                while ((size_t)dq[head] < j) head++;
                int m = dq[head];
                if (m != last_emitted) {
                    last_emitted = m;
                    uint64_t key = val[m];
                    int p = key_partition_pass(key, np, hp, hn);
                    if (p >= 0) {                                /* -1: another hash-range pass owns this key */
                        mrec *b = &wb[(size_t)p * WBUF_RECS + wn[p]];
                        mrec_pack(b, key, base + (uint64_t)m);
                        if (++wn[p] == WBUF_RECS) { pw_append(p, &wb[(size_t)p * WBUF_RECS], WBUF_RECS); wn[p] = 0; }
                        local_recs++;
                    }
                }
            }
        }
        /* end of chunk: push everything so a checkpoint after this chunk is complete */
        for (int p = 0; p < np; p++) if (wn[p]) { pw_append(p, &wb[(size_t)p * WBUF_RECS], wn[p]); wn[p] = 0; }
        local_win += nwin;
        atomic_store_explicit(&st->records, local_recs, memory_order_relaxed);
        atomic_store_explicit(&st->windows, local_win, memory_order_relaxed);

        pthread_mutex_lock(&ring_mtx);
        s->state = 0;
        chunk_done[seq % (uint64_t)nslots] = 1;
        /* chunks 0..done_contig-1 are all finished: advance over finished ones */
        while (done_contig < chunks_published && chunk_done[done_contig % (uint64_t)nslots]) done_contig++;
        pthread_cond_broadcast(&ring_can_fill);
        pthread_mutex_unlock(&ring_mtx);
    }
    free(wb); free(wn); free(rank); free(val); free(dq);
    return NULL;
}

/* ---------------- checkpoint ---------------- */
/*
 * do_checkpoint - flush and fsync every partition, then record the resume
 * position (start of the first chunk not yet fully processed) and the file
 * lengths.  Records from later chunks that were already flushed simply get
 * written again after a resume and are deduplicated by pimatchsort.
 */
static void do_checkpoint(int final)
{
    if (opt.dryrun) return;
    pthread_mutex_lock(&ring_mtx);
    uint64_t dc = done_contig;
    uint64_t next = dc == 0 ? resume_pos : chunk_end[(dc - 1) % (uint64_t)nslots];
    pthread_mutex_unlock(&ring_mtx);
    checkpoint c;
    c.nparts = opt.nparts;
    c.sizes = malloc((size_t)opt.nparts * sizeof *c.sizes);
    if (!c.sizes) die("out of memory");
    uint64_t recs = base_records;
    for (int i = 0; i < opt.nparts; i++) {
        pthread_mutex_lock(&pw[i].m);
        pw_flush_locked(&pw[i], i);
        if (fsync(pw[i].fd)) die("fsync failed on partition %d: %s", i, strerror(errno));
        c.sizes[i] = pw[i].bytes;
        recs += pw[i].records;
        pthread_mutex_unlock(&pw[i].m);
    }
    c.next_pos = next;
    c.records = recs;
    checkpoint_write(opt.dir, &c);
    free(c.sizes);
    char b[32];
    if (!final) printf("\n  checkpoint: resume position %s written\n", fmt_u64(next, b, sizeof b));
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);
    setvbuf(stdout, NULL, _IOLBF, 0);
    install_signals();
    winspec_init(&W, opt.L);
    POW10K = POW10[opt.k - 1];

    scratch_parse(opt.dirspec, &SC);
    opt.dir = SC.dirs[0];

    /* ---- source probe and summary ---- */
    printf("pimatchscan v%s - minimizer scan (pass 1 of 3)%s\n", PI_TOOLS_VERSION, opt.dryrun ? "  *** DRY RUN: nothing written ***" : "");
    pireader_open(&rd, opt.src, opt.L, opt.chunk, opt.ndigits);
    {
        char desc[256], b1[32], b2[32], b3[32];
        pistream_describe(&rd.ps, desc, sizeof desc);
        printf("  source         : %s\n  source format  : %s\n", opt.src, desc);
        uint64_t n = opt.ndigits == UINT64_MAX ? rd.digits_est : (opt.ndigits < rd.digits_est ? opt.ndigits : rd.digits_est);
        double frac = 2.0 / (opt.w + 1);
        double recs = (double)n * frac / opt.hpasses;
        double bytes = recs * MREC_BYTES;
        printf("  match length   : every repeat of >= %d digits is found (k=%d, w=%d); shorter ones only by luck\n", opt.L, opt.k, opt.w);
        printf("  digits to scan : %s\n", fmt_u64(n, b1, sizeof b1));
        printf("  sampling       : minimizer density ~2/(w+1) = %.1f%% of positions -> ~%.3g records%s\n", frac * 100, recs,
               opt.hpasses > 1 ? " in this pass" : "");
        if (opt.hpasses > 1)
            printf("  hash range     : pass %d of %d - keys in this slice only (1/%d of the records); run all %d passes\n",
                   opt.hpass, opt.hpasses, opt.hpasses, opt.hpasses);
        printf("  scratch        : ~%s in %d partitions of ~%s (%d-byte records)\n", fmt_bytes(bytes, b1, sizeof b1), opt.nparts,
               fmt_bytes(bytes / opt.nparts, b2, sizeof b2), MREC_BYTES);
        printf("  sort RAM       : pimatchsort needs ~3x the largest partition = %s\n", fmt_bytes(3.0 * bytes / opt.nparts, b1, sizeof b1));
        /* N^2/(2*10^k) for a uniform sample; minimizers prefer low-rank k-mers,
         * which collide with each other more often - measured ~3x (301 vs 99
         * on 200B digits at k=19), so scale the estimate accordingly */
        double false_pairs = 3.0 * recs * recs / (2.0 * pow(10.0, opt.k));
        printf("  false pairs    : ~%.3g pairs of equal %d-mers that are not %d-digit repeats (verified away in pass 3)\n",
               false_pairs, opt.k, opt.L);
        printf("  write rate     : %.1f MB of scratch per 100M digits/s of source\n", frac * MREC_BYTES * 100.0);
        printf("  threads        : %d workers + reader; write block %s per partition (%s buffered)\n", opt.threads,
               fmt_bytes((double)opt.block_bytes, b2, sizeof b2), fmt_bytes((double)opt.block_bytes * opt.nparts, b3, sizeof b3));
        /* scratch layout: on resume the manifest's split is authoritative */
        if (opt.resume) {
            manifest mm; manifest_read(opt.dir, &mm);
            if (strcmp(mm.sc.dirs[0], opt.dir)) die("resume: '%s' is not the primary directory of this scan ('%s')", opt.dir, mm.sc.dirs[0]);
            SC = mm.sc;
            manifest_check_dirs(&mm);
            printf("  scratch dirs   : %d (layout from the manifest)\n", SC.ndirs);
            for (int d = 0; d < SC.ndirs; d++)
                printf("                   %s: partitions %d..%d\n", SC.dirs[d], SC.first_part[d], SC.first_part[d + 1] - 1);
        } else {
            uint64_t freeb[MAX_DIRS]; double total_free = 0;
            for (int d = 0; d < SC.ndirs; d++) { freeb[d] = dir_free_bytes(SC.dirs[d]); total_free += (double)freeb[d]; }
            int fits = scratch_assign(&SC, opt.nparts, bytes / opt.nparts, freeb);
            printf("  scratch dirs   : %d, %s free in total, need ~%s (+5%% margin)%s\n", SC.ndirs, fmt_bytes(total_free, b1, sizeof b1),
                   fmt_bytes(bytes, b2, sizeof b2), fits ? "  *** DOES NOT FIT ***" : "");
            for (int d = 0; d < SC.ndirs; d++) {
                int np_d = SC.first_part[d + 1] - SC.first_part[d];
                printf("                   %s: %s free, %d partitions (%d..%d, ~%s)\n", SC.dirs[d], fmt_bytes((double)freeb[d], b1, sizeof b1),
                       np_d, SC.first_part[d], np_d ? SC.first_part[d + 1] - 1 : SC.first_part[d] - 1,
                       fmt_bytes(bytes / opt.nparts * np_d, b2, sizeof b2));
            }
            if (fits && !opt.dryrun) die("not enough free space: need ~%s across the scratch set, have %s",
                                         fmt_bytes(bytes * 1.05, b2, sizeof b2), fmt_bytes(total_free, b3, sizeof b3));
        }
        if (false_pairs > 1e9) warn("k=%d is small for this many records: pass 3 will have to check ~%.3g pairs", opt.k, false_pairs);
    }

    /* ---- manifest / resume ---- */
    manifest m;
    checkpoint ck;
    memset(&ck, 0, sizeof ck);
    int have_ck = 0;
    if (opt.resume) {
        manifest_read(opt.dir, &m);
        if (m.L != opt.L || m.k != opt.k || m.nparts != opt.nparts || strcmp(m.source, opt.src) || m.hpass != opt.hpass || m.hpasses != opt.hpasses)
            die("resume: options differ from the manifest in '%s' (L=%d k=%d P=%d H=%d/%d source=%s)", opt.dir, m.L, m.k, m.nparts, m.hpass, m.hpasses, m.source);
        have_ck = checkpoint_read(opt.dir, &ck, opt.nparts);
        if (!have_ck) die("resume requested but no checkpoint in '%s'", opt.dir);
        if (m.complete) die("the scan in '%s' already completed", opt.dir);
    } else if (!opt.dryrun) {
        char p[2048];
        manifest_path(opt.dir, p, sizeof p);
        if (file_size_of(p) && !opt.force) die("'%s' already holds a scan (manifest present) - use -R to resume or -f to overwrite", opt.dir);
        memset(&m, 0, sizeof m);
        snprintf(m.version, sizeof m.version, "%s", PI_TOOLS_VERSION);
        snprintf(m.source, sizeof m.source, "%s", opt.src);
        m.L = opt.L; m.k = opt.k; m.w = opt.w; m.ndigits = opt.ndigits; m.first_pos = rd.first_pos;
        m.nparts = opt.nparts; m.complete = 0;
        m.hpass = opt.hpass; m.hpasses = opt.hpasses; m.sc = SC;
        manifest_write(opt.dir, &m);
        checkpoint_path(opt.dir, p, sizeof p);
        unlink(p);                                   /* a fresh scan starts without a checkpoint */
    }

    /* ---- partition writers ---- */
    pw = calloc((size_t)opt.nparts, sizeof *pw);
    if (!pw) die("out of memory");
    for (int i = 0; i < opt.nparts; i++) {
        pthread_mutex_init(&pw[i].m, NULL);
        pw[i].fd = -1;
        if (!opt.dryrun) {
            char p[2048];
            part_path(&SC, i, p, sizeof p);
            pw[i].fd = open(p, O_WRONLY | O_CREAT | (opt.resume ? 0 : O_TRUNC), 0644);
            if (pw[i].fd < 0) die("cannot open partition file '%s': %s", p, strerror(errno));
            if (opt.resume) {
                uint64_t have = file_size_of(p);
                if (have < ck.sizes[i]) die("partition file '%s' is shorter (%" PRIu64 ") than the checkpoint recorded (%" PRIu64 ") - scratch is damaged, restart from scratch",
                                            p, have, ck.sizes[i]);
                if (ftruncate(pw[i].fd, (off_t)ck.sizes[i])) die("truncate failed on '%s': %s", p, strerror(errno));
                if (lseek(pw[i].fd, (off_t)ck.sizes[i], SEEK_SET) < 0) die("lseek failed on '%s'", p);
                pw[i].bytes = ck.sizes[i];
            }
            pw[i].cap = (size_t)opt.block_bytes;
            pw[i].buf = malloc(pw[i].cap);
            if (!pw[i].buf) die("cannot allocate partition write buffers (%d x %" PRIu64 " bytes)", opt.nparts, opt.block_bytes);
        }
    }

    /* ---- position the reader ---- */
    resume_pos = rd.first_pos;
    if (opt.resume) {
        resume_pos = ck.next_pos;
        base_records = ck.records;
        char b1[32], b2[32];
        printf("  resuming       : at position %s (%s records already written)\n", fmt_u64(resume_pos, b1, sizeof b1), fmt_u64(ck.records, b2, sizeof b2));
        if (pistream_seek(&rd.ps, resume_pos) != 0) die("this source cannot be seeked (text with non-digit bytes) - cannot resume");
        rd.base = rd.next_pos = resume_pos;
        if (rd.limit != UINT64_MAX) {
            uint64_t skipped = resume_pos - rd.first_pos;
            rd.limit = rd.limit > skipped ? rd.limit - skipped : 0;
        }
    }

    /* ---- ring ---- */
    nslots = 2 * opt.threads + 2;
    slots = calloc((size_t)nslots, sizeof *slots);
    chunk_end = calloc((size_t)nslots, sizeof *chunk_end);
    chunk_done = calloc((size_t)nslots, sizeof *chunk_done);
    if (!slots || !chunk_end || !chunk_done) die("out of memory");
    for (int i = 0; i < nslots; i++) {
        slots[i].buf = malloc(opt.chunk + (size_t)opt.L + 64);
        if (!slots[i].buf) die("cannot allocate chunk ring (%d x %zu bytes)", nslots, opt.chunk);
    }

    printf("Scanning ...\n");
    double t_start = now_sec(), t_last = t_start, t_ck = t_start;
    uint64_t last_wins = 0, last_raw = 0, last_recs = 0;
    const int tty = isatty(1);
    int status_shown = 0;

    pthread_t rth, wth[MAX_WORKERS];
    warg_t wargs[MAX_WORKERS];
    if (pthread_create(&rth, NULL, reader_main, NULL)) die("pthread_create(reader) failed");
    for (int i = 0; i < opt.threads; i++) {
        wargs[i].id = i;
        if (pthread_create(&wth[i], NULL, worker_main, &wargs[i])) die("pthread_create(worker) failed");
    }

    int done = 0;
    while (!done) {
        struct timespec ts = { 0, 200 * 1000 * 1000 };
        nanosleep(&ts, NULL);
        pthread_mutex_lock(&ring_mtx);
        done = producing_done && done_contig >= chunks_published;
        uint64_t wins = g_windows, pos = g_last_pos;
        pthread_mutex_unlock(&ring_mtx);
        double t = now_sec();
        if (!done && t - t_last >= opt.progress_secs) {
            uint64_t recs = base_records, raw = rd.ps.total_raw, done_wins = 0;
            for (int i = 0; i < opt.threads; i++) {
                recs += atomic_load_explicit(&wstats[i].records, memory_order_relaxed);
                done_wins += atomic_load_explicit(&wstats[i].windows, memory_order_relaxed);
            }
            char c1[32], c2[32], d1[32], d2[32];
            double dt = t - t_last, tot = t - t_start;
            double dig_rate = dt > 0 ? (double)(wins - last_wins) / dt : 0;
            double total = (double)(opt.ndigits == UINT64_MAX ? rd.digits_est : (opt.ndigits < rd.digits_est ? opt.ndigits : rd.digits_est));
            double scanned = (double)(pos + 1 - rd.first_pos);
            double pct = total > 0 ? 100.0 * scanned / total : 0;
            double eta = (dig_rate > 0 && total > scanned) ? (total - scanned) / dig_rate : -1;
            printf("%s[%s] pos %s (%.1f%% done)  records %s (%.1f%% of digits)  %.0f M digits/s (%.0f MB/s read, %.0f MB/s written)%s%s%s",
                   tty ? "\r\033[K" : "", fmt_duration(tot, d1, sizeof d1), fmt_u64(pos, c1, sizeof c1), pct,
                   fmt_u64(recs, c2, sizeof c2), done_wins ? 100.0 * (double)(recs - base_records) / (double)done_wins : 0,
                   dig_rate / 1e6, (double)(raw - last_raw) / dt / 1e6, (double)(recs - last_recs) * MREC_BYTES / dt / 1e6,
                   eta >= 0 ? "  ETA " : "", eta >= 0 ? fmt_duration(eta, d2, sizeof d2) : "", tty ? "" : "\n");
            fflush(stdout);
            status_shown = tty;
            t_last = t; last_wins = wins; last_raw = raw; last_recs = recs;
        }
        if (!done && opt.ckpt_secs > 0 && t - t_ck >= opt.ckpt_secs) {
            if (status_shown) { printf("\n"); status_shown = 0; }
            do_checkpoint(0);
            t_ck = now_sec();
        }
    }
    pthread_join(rth, NULL);
    for (int i = 0; i < opt.threads; i++) pthread_join(wth[i], NULL);
    int interrupted = g_stop != 0;
    if (status_shown) printf("\n");

    /* final flush / checkpoint / manifest */
    uint64_t recs = base_records, bytes = 0, minp = UINT64_MAX, maxp = 0;
    do_checkpoint(1);
    for (int i = 0; i < opt.nparts; i++) {
        recs += pw[i].records;
        uint64_t pb = opt.dryrun ? pw[i].records * MREC_BYTES : pw[i].bytes;
        bytes += pb;
        if (pb < minp) minp = pb;
        if (pb > maxp) maxp = pb;
        if (pw[i].fd >= 0) close(pw[i].fd);
    }
    if (!opt.dryrun) {
        pthread_mutex_lock(&ring_mtx);
        m.end_pos = done_contig == 0 ? resume_pos : chunk_end[(done_contig - 1) % (uint64_t)nslots];
        pthread_mutex_unlock(&ring_mtx);
        m.complete = !interrupted;
        manifest_write(opt.dir, &m);
    }

    double elapsed = now_sec() - t_start;
    char b1[32], b2[32], d1[32];
    printf("\n%s\n", interrupted ? "*** INTERRUPTED - checkpoint written; rerun with -R to resume ***" : "Scan complete");
    printf("  windows scanned : %s (positions %s..%s)\n", fmt_u64(g_windows, b1, sizeof b1), fmt_u64(resume_pos, d1, sizeof d1), fmt_u64(g_last_pos, b2, sizeof b2));
    printf("  records         : %s total (%.2f%% of positions; expected ~%.2f%%%s)\n", fmt_u64(recs, b1, sizeof b1),
           g_windows ? 100.0 * (double)(recs - base_records) / (double)g_windows : 0, 200.0 / (opt.w + 1) / opt.hpasses,
           opt.hpasses > 1 ? " for this hash pass" : "");
    printf("  scratch         : %s in %d partitions, %s .. %s each\n", fmt_bytes((double)bytes, b1, sizeof b1), opt.nparts,
           fmt_bytes((double)minp, b2, sizeof b2), fmt_bytes((double)maxp, d1, sizeof d1));
    printf("  elapsed         : %s  (%.0f M digits/s, %.0f MB/s of records)\n", fmt_duration(elapsed, d1, sizeof d1),
           elapsed > 0 ? g_windows / elapsed / 1e6 : 0, elapsed > 0 ? (double)(recs - base_records) * MREC_BYTES / elapsed / 1e6 : 0);
    printf("  peak RSS        : %.0f MiB\n", peak_rss_mb());
    if (!interrupted && !opt.dryrun) printf("\nNext: pimatchsort -d %s\n", opt.dir);
    pireader_close(&rd);
    return interrupted ? 130 : 0;
}
