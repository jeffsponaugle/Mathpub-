/*
 * pimatchsort - pass 2 of the minimizer/sort repeat finder (see pimatch.h).
 *
 * For each partition file written by pimatchscan: load it into RAM, radix
 * sort the 14-byte records by key, and write every key that occurs at more
 * than one distinct position as a candidate group "key pos1 pos2 ...".
 * Partitions larger than the RAM budget are processed in sub-passes that
 * each keep only the records whose hash falls in a slice (the file is read
 * once per slice).  Each partition's candidates go to cand_NNNN.txt with a
 * .done marker, so -R resumes after a crash; the final candidates.txt is the
 * concatenation.
 *
 * Loading overlaps sorting: three sort buffers rotate between a loader thread
 * (reads the next unit - partition or sub-pass - straight from disk into a
 * free buffer) and the main thread (radix sorts the loaded buffer into a
 * second one and emits the groups).  With the disk and the sort each taking
 * about half of the old per-partition time this roughly doubles pass 2.
 *
 * Build: cc -O3 -march=native -pthread -o pimatchsort pimatchsort.c -lm
 */
#include "pimatch.h"            /* records, manifest, digit helpers */
#include <pthread.h>            /* parallel radix sort */

#define MAX_THREADS 128

static struct {
    const char *dir, *out;
    uint64_t ram;               /* budget for the three sort buffers together */
    int threads, resume, delete_done;
    double progress_secs;
} opt = { NULL, NULL, 0, 0, 0, 0, 10.0 };

static void usage(const char *argv0)
{
    fprintf(stderr,
"pimatchsort v" PI_TOOLS_VERSION " - sort partitions, emit candidate groups (pass 2 of 3)\n"
"Usage: %s -d <scratchdir> [options]\n"
"\n"
"  -d DIR      scratch directory written by pimatchscan (the primary one; the\n"
"              manifest names the others of a multi-directory set)\n"
"  -D          delete each partition file once its candidates are written\n"
"              (frees the scratch for the next hash-range pass)\n"
"  -m SIZE     RAM budget for sorting [50%% of physical RAM]; a partition needs\n"
"              ~3x its size (two sort buffers + one being loaded for the next\n"
"              partition), larger ones are done in sub-passes (re-reading the file)\n"
"  -o FILE     candidate output file [DIR/candidates.txt]\n"
"  -t N        sort threads [min(8, CPUs)]\n"
"  -R          resume: skip partitions that already have a .done marker\n"
"  -I SECS     progress interval [10]\n"
"  -V          print version and exit\n"
"  -h          this help\n", argv0);
    exit(1);
}

static void parse_args(int argc, char **argv)
{
    int c;
    long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    opt.threads = ncpu >= 8 ? 8 : (int)(ncpu > 0 ? ncpu : 1);
    while ((c = getopt(argc, argv, "d:m:o:t:RDI:Vh")) != -1) {
        switch (c) {
        case 'd': opt.dir = primary_dir(optarg); break;
        case 'D': opt.delete_done = 1; break;
        case 'm': if (parse_mem_size(optarg, &opt.ram) || opt.ram < (64u << 20)) die("bad -m"); break;
        case 'o': opt.out = optarg; break;
        case 't': opt.threads = atoi(optarg); if (opt.threads < 1 || opt.threads > MAX_THREADS) die("-t must be 1..%d", MAX_THREADS); break;
        case 'R': opt.resume = 1; break;
        case 'I': opt.progress_secs = atof(optarg); if (opt.progress_secs <= 0) die("bad -I"); break;
        case 'V': print_version("pimatchsort"); exit(0);
        default: usage(argv[0]);
        }
    }
    if (!opt.dir) { fprintf(stderr, "missing -d\n"); usage(argv[0]); }
    if (!opt.ram) { uint64_t phys = phys_mem_bytes(); opt.ram = phys ? phys / 2 : (8ULL << 30); }
}

/* ---------------- parallel LSD radix sort on the 64-bit key ---------------- */
typedef struct {
    const mrec *src; mrec *dst;
    size_t from, to;
    int shift;
    uint64_t hist[65536];
    uint64_t start[65536];      /* this thread's write cursor per bucket */
} rjob;

static inline unsigned rdigit(const mrec *r, int shift)
{
    /* 16-bit digit of the little-endian key at byte offset shift/8 */
    return (unsigned)r->b[shift / 8] | ((unsigned)r->b[shift / 8 + 1] << 8);
}

static void *radix_hist(void *a)
{
    rjob *j = a;
    memset(j->hist, 0, sizeof j->hist);
    for (size_t i = j->from; i < j->to; i++) j->hist[rdigit(&j->src[i], j->shift)]++;
    return NULL;
}
static void *radix_scatter(void *a)
{
    rjob *j = a;
    for (size_t i = j->from; i < j->to; i++) {
        unsigned d = rdigit(&j->src[i], j->shift);
        j->dst[j->start[d]++] = j->src[i];
    }
    return NULL;
}

/*
 * radix_sort - sort n records by key using 4 passes of 16 bits, parallel
 * histogram + scatter across nt threads.  Passes whose digit is constant
 * (e.g. the top bits of a 19-digit key) are skipped.  Returns the buffer
 * holding the sorted data (a or b).
 */
static mrec *radix_sort(mrec *a, mrec *b, size_t n, int nt, rjob *jobs, pthread_t *th)
{
    mrec *src = a, *dst = b;
    for (int shift = 0; shift < 64; shift += 16) {
        for (int t = 0; t < nt; t++) {
            jobs[t].src = src; jobs[t].dst = dst; jobs[t].shift = shift;
            jobs[t].from = n * (size_t)t / (size_t)nt; jobs[t].to = n * (size_t)(t + 1) / (size_t)nt;
            pthread_create(&th[t], NULL, radix_hist, &jobs[t]);
        }
        for (int t = 0; t < nt; t++) pthread_join(th[t], NULL);
        /* skip a pass if every record has the same digit */
        int nonzero = 0;
        for (unsigned d = 0; d < 65536 && nonzero < 2; d++) {
            uint64_t tot = 0;
            for (int t = 0; t < nt; t++) tot += jobs[t].hist[d];
            if (tot) nonzero++;
        }
        if (nonzero < 2) continue;
        uint64_t pos = 0;
        for (unsigned d = 0; d < 65536; d++)
            for (int t = 0; t < nt; t++) { jobs[t].start[d] = pos; pos += jobs[t].hist[d]; }
        for (int t = 0; t < nt; t++) pthread_create(&th[t], NULL, radix_scatter, &jobs[t]);
        for (int t = 0; t < nt; t++) pthread_join(th[t], NULL);
        mrec *tmp = src; src = dst; dst = tmp;
    }
    return src;
}

/* ---------------- candidate emission ---------------- */
static uint64_t st_groups, st_pairs, st_dupes, st_maxgroup, st_records;

static int cmp_u64(const void *x, const void *y) { uint64_t a = *(const uint64_t *)x, b = *(const uint64_t *)y; return a < b ? -1 : a > b; }

/*
 * emit_groups - walk sorted records; every key with >= 2 distinct positions
 * becomes one line "key pos1 pos2 ..." (positions ascending, duplicates
 * removed).
 */
static void emit_groups(const mrec *s, size_t n, FILE *out)
{
    uint64_t *poss = NULL; size_t cap = 0;
    size_t i = 0;
    while (i < n) {
        uint64_t key = mrec_key(&s[i]);
        size_t j = i;
        while (j < n && mrec_key(&s[j]) == key) j++;
        size_t m = j - i;
        if (m >= 2) {
            if (m > cap) { cap = m + 16; poss = realloc(poss, cap * sizeof *poss); if (!poss) die("out of memory"); }
            for (size_t q = 0; q < m; q++) poss[q] = mrec_pos(&s[i + q]);
            qsort(poss, m, sizeof *poss, cmp_u64);
            size_t u = 0;
            for (size_t q = 0; q < m; q++) if (u == 0 || poss[q] != poss[u - 1]) poss[u++] = poss[q]; else st_dupes++;
            if (u >= 2) {
                fprintf(out, "%" PRIu64, key);
                for (size_t q = 0; q < u; q++) fprintf(out, " %" PRIu64, poss[q]);
                fputc('\n', out);
                st_groups++;
                st_pairs += u * (u - 1) / 2;
                if (u > st_maxgroup) st_maxgroup = u;
            }
        }
        i = j;
    }
    free(poss);
}


/* ---------------- work units and the loader thread ---------------- */
/*
 * A unit is one load+sort+emit step: a whole partition (subs == 0) or one
 * hash slice of a partition that does not fit in a buffer (pass 0..2^subs-1).
 * The loader thread fills the buffer the main thread hands it and publishes
 * the record count; the main thread meanwhile sorts and emits the previous
 * unit.  Exactly three buffers are needed: sort source, sort destination and
 * the one being loaded.
 */
typedef struct { int part, pass, subs; uint64_t nrec_file; } unit_t;

static struct {
    pthread_mutex_t mu;
    pthread_cond_t cv;
    unit_t *units; int nunits;
    const scratch *sc;          /* where the partition files live */
    size_t cap;                 /* records per buffer */
    mrec *target;               /* buffer to fill next (assigned by main; NULL = none) */
    int ready;                  /* a loaded unit is waiting to be sorted */
    mrec *ready_buf; size_t ready_n; double ready_secs;
    int finished;               /* loader exited (all units done, or stopped) */
    uint8_t *stage;             /* 64 MiB staging area for filtered (sub-pass) loads */
} ld;

#define STAGE_BYTES (64u << 20)

/*
 * load_unit - read one unit into buf.  Unfiltered units are read straight
 * from the file into the buffer (no staging copy); sub-pass units stream
 * through the staging area keeping only the records whose key hash falls in
 * this pass's slice.  Returns the record count, or SIZE_MAX if interrupted.
 */
static size_t load_unit(const unit_t *u, mrec *buf)
{
    char p[2048]; part_path(ld.sc, u->part, p, sizeof p);
    int fd = open(p, O_RDONLY);
    if (fd < 0) die("cannot open '%s': %s", p, strerror(errno));
#ifdef POSIX_FADV_SEQUENTIAL
    posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
#ifdef F_RDAHEAD
    fcntl(fd, F_RDAHEAD, 1);
#endif
    uint64_t left = u->nrec_file * MREC_BYTES;
    size_t n = 0;
    if (u->subs == 0) {
        if (u->nrec_file > ld.cap) die("internal: partition %d larger than the sort buffer", u->part);
        uint8_t *dst = (uint8_t *)buf;
        while (left > 0) {
            size_t want = left < STAGE_BYTES ? (size_t)left : STAGE_BYTES;
            ssize_t r = read(fd, dst, want);
            if (r < 0) { if (errno == EINTR) continue; die("read error on '%s': %s", p, strerror(errno)); }
            if (r == 0) die("'%s' ended early", p);
            dst += r; left -= (uint64_t)r;
            if (g_stop) { close(fd); return SIZE_MAX; }
        }
        n = u->nrec_file;
    } else {
        size_t carry = 0;
        unsigned mask = (1u << u->subs) - 1;
        while (left > 0) {
            size_t want = left < STAGE_BYTES - carry ? (size_t)left : STAGE_BYTES - carry;
            ssize_t r = read(fd, ld.stage + carry, want);
            if (r < 0) { if (errno == EINTR) continue; die("read error on '%s': %s", p, strerror(errno)); }
            if (r == 0) die("'%s' ended early", p);
            left -= (uint64_t)r;
            size_t have = carry + (size_t)r, nw = have / MREC_BYTES;
            const mrec *in = (const mrec *)ld.stage;
            for (size_t q = 0; q < nw; q++) {
                uint64_t h = mix64(mrec_key(&in[q]) + 0x1234567887654321ULL);
                if ((int)(h & mask) != u->pass) continue;
                if (n >= ld.cap) die("sub-pass %d of partition %d overflowed the buffer (uneven hash split) - rerun with more RAM or more partitions", u->pass, u->part);
                buf[n++] = in[q];
            }
            carry = have - nw * MREC_BYTES;
            if (carry) memmove(ld.stage, ld.stage + nw * MREC_BYTES, carry);
            if (g_stop) { close(fd); return SIZE_MAX; }
        }
    }
    close(fd);
    return n;
}

/*
 * loader_main - loader thread body: for every unit wait for main to assign a
 * free buffer, load into it, publish (ready_buf/ready_n/ready_secs).  Stops
 * early on g_stop; always sets finished on exit so main never waits forever.
 */
static void *loader_main(void *arg)
{
    (void)arg;
    for (int i = 0; i < ld.nunits; i++) {
        pthread_mutex_lock(&ld.mu);
        while (!ld.target && !g_stop) pthread_cond_wait(&ld.cv, &ld.mu);
        mrec *buf = ld.target; ld.target = NULL;
        pthread_mutex_unlock(&ld.mu);
        if (!buf) break;                                    /* g_stop */
        double t = now_sec();
        size_t n = load_unit(&ld.units[i], buf);
        if (n == SIZE_MAX) break;
        pthread_mutex_lock(&ld.mu);
        ld.ready = 1; ld.ready_buf = buf; ld.ready_n = n; ld.ready_secs = now_sec() - t;
        pthread_cond_broadcast(&ld.cv);
        pthread_mutex_unlock(&ld.mu);
    }
    pthread_mutex_lock(&ld.mu);
    ld.finished = 1;
    pthread_cond_broadcast(&ld.cv);
    pthread_mutex_unlock(&ld.mu);
    return NULL;
}

/* assign_target - hand the loader a free buffer for its next unit */
static void assign_target(mrec *buf)
{
    pthread_mutex_lock(&ld.mu);
    ld.target = buf;
    pthread_cond_broadcast(&ld.cv);
    pthread_mutex_unlock(&ld.mu);
}

/*
 * take_ready - wait for the loader's next unit; returns 0 and the buffer +
 * count, or -1 if the loader stopped without producing it (interrupted).
 */
static int take_ready(mrec **buf, size_t *n, double *secs)
{
    pthread_mutex_lock(&ld.mu);
    while (!ld.ready && !ld.finished) pthread_cond_wait(&ld.cv, &ld.mu);
    int ok = ld.ready;
    if (ok) { *buf = ld.ready_buf; *n = ld.ready_n; *secs = ld.ready_secs; ld.ready = 0; }
    pthread_mutex_unlock(&ld.mu);
    return ok ? 0 : -1;
}

int main(int argc, char **argv)
{
    parse_args(argc, argv);
    setvbuf(stdout, NULL, _IOLBF, 0);
    install_signals();

    manifest m;
    manifest_read(opt.dir, &m);
    manifest_check_dirs(&m);
    char outname[2048];
    if (!opt.out) { snprintf(outname, sizeof outname, "%s/candidates.txt", opt.dir); opt.out = outname; }

    char b1[32], b2[32], b3[32];
    printf("pimatchsort v%s - sort partitions, emit candidate groups (pass 2 of 3)\n", PI_TOOLS_VERSION);
    printf("  scratch dir    : %s (%d partitions, L=%d k=%d w=%d, scan %s)\n", opt.dir, m.nparts, m.L, m.k, m.w,
           m.complete ? "complete" : "INCOMPLETE - results will be partial");
    if (!m.complete) warn("the scan in '%s' did not complete; candidates cover only the scanned range", opt.dir);
    if (m.hpasses > 1) printf("  hash range     : pass %d of %d (candidates cover this slice of the key space only)\n", m.hpass, m.hpasses);
    if (m.sc.ndirs > 1) {
        printf("  scratch set    : %d directories\n", m.sc.ndirs);
        for (int d = 0; d < m.sc.ndirs; d++)
            printf("                   %s: partitions %d..%d\n", m.sc.dirs[d], m.sc.first_part[d], m.sc.first_part[d + 1] - 1);
    }

    /* sizes */
    uint64_t total = 0, maxsz = 0;
    for (int i = 0; i < m.nparts; i++) {
        char p[2048], dp[2048]; part_path(&m.sc, i, p, sizeof p); done_path(opt.dir, i, dp, sizeof dp);
        if (access(p, F_OK) != 0) {
            if (opt.resume && access(dp, F_OK) == 0) continue;      /* sorted earlier and deleted with -D */
            die("partition file '%s' is missing", p);
        }
        uint64_t sz = file_size_of(p);
        if (sz % MREC_BYTES) warn("partition %d has a partial record (%" PRIu64 " bytes); the tail is ignored", i, sz);
        total += sz; if (sz > maxsz) maxsz = sz;
    }
    printf("  scratch data   : %s total, largest partition %s, RAM budget %s (%d sort threads)\n",
           fmt_bytes((double)total, b1, sizeof b1), fmt_bytes((double)maxsz, b2, sizeof b2), fmt_bytes((double)opt.ram, b3, sizeof b3), opt.threads);
    uint64_t per_buf = opt.ram / 3;
    if (maxsz > per_buf)
        warn("largest partition (%s) exceeds a third of the RAM budget (%s): it will be sorted in %d sub-passes, re-reading the file each time",
             fmt_bytes((double)maxsz, b1, sizeof b1), fmt_bytes((double)per_buf, b2, sizeof b2),
             1 << (int)ceil(log2((double)maxsz / (double)per_buf)));
    uint64_t bufsz = maxsz < per_buf ? maxsz : per_buf;
    if (bufsz < (1u << 20)) bufsz = 1u << 20;
    bufsz = (bufsz / MREC_BYTES + 1) * MREC_BYTES;
    check_fits_in_ram(3 * bufsz, "sort buffers");
    mrec *bufs[3];
    bufs[0] = alloc_huge(bufsz, 0, "sort buffer A");
    bufs[1] = alloc_huge(bufsz, 0, "sort buffer B");
    bufs[2] = alloc_huge(bufsz, 0, "load buffer C");
    printf("  buffers        : 3 x %s (sort source, sort destination, next partition loading)\n", fmt_bytes((double)bufsz, b1, sizeof b1));
    size_t cap = bufsz / MREC_BYTES;
    rjob *jobs = calloc((size_t)opt.threads, sizeof *jobs);
    pthread_t *th = calloc((size_t)opt.threads, sizeof *th);
    if (!jobs || !th) die("out of memory");

    /* build the unit list (resume skips finished partitions) */
    int skipped = 0;
    uint64_t done_bytes = 0;
    unit_t *units = calloc((size_t)m.nparts * 64, sizeof *units);
    if (!units) die("out of memory");
    int nunits = 0;
    for (int i = 0; i < m.nparts; i++) {
        char p[2048], dp[2048];
        part_path(&m.sc, i, p, sizeof p); done_path(opt.dir, i, dp, sizeof dp);
        uint64_t sz = file_size_of(p);
        if (opt.resume && access(dp, F_OK) == 0) { skipped++; done_bytes += sz; continue; }
        uint64_t nrec_file = sz / MREC_BYTES;
        int subs = 0;
        while (((uint64_t)cap << subs) < nrec_file) subs++;          /* 2^subs passes */
        if (subs > 6) die("partition %d would need %d sub-passes; give more RAM (-m) or rescan with more partitions", i, 1 << subs);
        for (int pass = 0; pass < (1 << subs); pass++)
            units[nunits++] = (unit_t){ i, pass, subs, nrec_file };
    }

    /* start the loader on the first unit */
    pthread_mutex_init(&ld.mu, NULL); pthread_cond_init(&ld.cv, NULL);
    ld.units = units; ld.nunits = nunits; ld.cap = cap; ld.sc = &m.sc;
    ld.stage = malloc(STAGE_BYTES);
    if (!ld.stage) die("out of memory");
    mrec *freebuf[3] = { bufs[0], bufs[1], bufs[2] }; int nfree = 3;
    pthread_t loader;
    if (nunits > 0) assign_target(freebuf[--nfree]);
    if (pthread_create(&loader, NULL, loader_main, NULL)) die("cannot create loader thread");

    double t0 = now_sec();
    uint64_t run_bytes = 0;                                  /* bytes sorted by this run (rate/ETA) */
    int tty = isatty(1), interrupted = 0;
    FILE *out = NULL;
    uint64_t g0 = 0;
    double part_load = 0, part_sort = 0;
    for (int u = 0; u < nunits && !interrupted; u++) {
        const unit_t *un = &units[u];
        char cp[2048], dp[2048];
        cand_path(opt.dir, un->part, cp, sizeof cp); done_path(opt.dir, un->part, dp, sizeof dp);
        if (un->pass == 0) {
            out = fopen(cp, "w");
            if (!out) die("cannot create '%s': %s", cp, strerror(errno));
            setvbuf(out, NULL, _IOFBF, 1 << 20);
            g0 = st_groups; part_load = part_sort = 0;
        }
        mrec *src; size_t n; double lsecs;
        if (tty) { printf("\r\033[K  partition %d/%d%s: waiting for the loader ...", un->part + 1, m.nparts, un->subs ? " (sub-pass)" : ""); fflush(stdout); }
        if (take_ready(&src, &n, &lsecs) < 0) { interrupted = 1; break; }
        part_load += lsecs;
        mrec *dst = freebuf[--nfree];
        if (u + 1 < nunits) assign_target(freebuf[--nfree]);   /* loader fills the third buffer meanwhile */
        st_records += n;
        if (tty) { printf("\r\033[K  partition %d/%d%s: %s records loaded, sorting ...", un->part + 1, m.nparts,
                          un->subs ? " (sub-pass)" : "", fmt_u64(n, b1, sizeof b1)); fflush(stdout); }
        double ts = now_sec();
        mrec *sorted = radix_sort(src, dst, n, opt.threads, jobs, th);
        emit_groups(sorted, n, out);
        part_sort += now_sec() - ts;
        freebuf[nfree++] = src; freebuf[nfree++] = dst;
        if (g_stop) interrupted = 1;
        if (un->pass + 1 < (1 << un->subs)) continue;        /* more sub-passes of this partition */

        if (fflush(out) || fsync(fileno(out)) || fclose(out)) die("write error on '%s'", cp);
        out = NULL;
        if (interrupted) break;
        FILE *d = fopen(dp, "w"); if (!d) die("cannot create '%s'", dp); fclose(d);
        if (opt.delete_done) {
            char pp[2048]; part_path(&m.sc, un->part, pp, sizeof pp);
            if (unlink(pp)) warn("cannot delete '%s': %s", pp, strerror(errno));
        }
        uint64_t sz = un->nrec_file * MREC_BYTES;
        done_bytes += sz; run_bytes += sz;
        double t = now_sec(), el = t - t0;
        double rate = el > 0 ? (double)run_bytes / el : 0;
        double eta = rate > 0 ? (double)(total - done_bytes) / rate : 0;
        char d1[32], d2[32];
        printf("%s[%s] partition %d/%d: %s, %s groups (+%s)  load %.0fs sort %.0fs  %.0f MB/s  ETA %s%s",
               tty ? "\r\033[K" : "", fmt_duration(el, d1, sizeof d1), un->part + 1, m.nparts, fmt_bytes((double)sz, b1, sizeof b1),
               fmt_u64(st_groups, b2, sizeof b2), fmt_u64(st_groups - g0, b3, sizeof b3), part_load, part_sort,
               rate / 1e6, fmt_duration(eta, d2, sizeof d2), tty ? "" : "\n");
        fflush(stdout);
    }
    if (out) fclose(out);                                    /* interrupted mid-partition: no .done marker */
    if (interrupted) g_stop = 1;                             /* release a loader waiting for a buffer */
    pthread_mutex_lock(&ld.mu); pthread_cond_broadcast(&ld.cv); pthread_mutex_unlock(&ld.mu);
    pthread_join(loader, NULL);
    if (tty) printf("\n");
    if (interrupted) { printf("*** INTERRUPTED - rerun with -R to continue from the next partition ***\n"); return 130; }

    /* concatenate per-partition candidates */
    out = fopen(opt.out, "w");
    if (!out) die("cannot create '%s': %s", opt.out, strerror(errno));
    fprintf(out, "# pimatchsort v%s L=%d k=%d w=%d source=%s nparts=%d hash_pass=%d/%d  format: key pos1 pos2 ... (positions 1-based after the decimal point)\n",
            PI_TOOLS_VERSION, m.L, m.k, m.w, m.source, m.nparts, m.hpass, m.hpasses);
    uint64_t lines = 0;
    for (int i = 0; i < m.nparts; i++) {
        char cp[2048]; cand_path(opt.dir, i, cp, sizeof cp);
        FILE *in = fopen(cp, "r");
        if (!in) die("missing '%s'", cp);
        char line[4096];
        while (fgets(line, sizeof line, in)) { fputs(line, out); lines++; }
        fclose(in);
    }
    if (fclose(out)) die("write error on '%s'", opt.out);

    double el = now_sec() - t0;
    char d1[32];
    printf("\nSort complete\n");
    printf("  records sorted  : %s%s\n", fmt_u64(st_records, b1, sizeof b1), skipped ? " (this run; some partitions were resumed)" : "");
    printf("  candidate groups: %s (%s pairs to verify, largest group %s, %s duplicate records dropped)\n",
           fmt_u64(lines, b1, sizeof b1), fmt_u64(st_pairs, b2, sizeof b2), fmt_u64(st_maxgroup, b3, sizeof b3), fmt_u64(st_dupes, d1, sizeof d1));
    printf("  elapsed         : %s  (%.0f MB/s%s)\n", fmt_duration(el, d1, sizeof d1), el > 0 ? (double)run_bytes / el / 1e6 : 0,
           skipped ? ", this run only" : "");
    printf("  peak RSS        : %.0f MiB\n", peak_rss_mb());
    printf("\nCandidates: %s\nNext: pimatchverify -d %s\n", opt.out, opt.dir);
    return 0;
}
