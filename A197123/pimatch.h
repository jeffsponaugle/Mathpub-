/*
 * pimatch.h - shared definitions for the minimizer / partition / sort
 * repeat finder: pimatchscan, pimatchsort, pimatchverify.
 *
 * Idea: stream pi once; for every window of w consecutive k-mers keep the
 * one with the smallest (hashed) rank - the "minimizer" - and write a record
 * (k-mer value, position) into one of P partition files chosen by a hash of
 * the value.  Two identical strings of length >= L = k + w - 1 contain the
 * same window of k-mers and therefore the same minimizer, so their records
 * carry equal keys and land in the same partition.  Sorting each partition
 * (pimatchsort) makes equal keys adjacent; pimatchverify re-reads the digits
 * around each candidate pair and extends the match to its full length.
 *
 * Records are 14 bytes: an 8-byte key (the k-mer as a decimal integer - exact
 * for k <= 19, so there are no hash false positives) and a 6-byte position
 * (1-based after the decimal point, up to 2^48).
 */
#ifndef PIMATCH_H
#define PIMATCH_H

#include "pi_common.h"          /* digit sources, hashing, timing, formatting, versioning */
#include <sys/statvfs.h>        /* statvfs: free space check for the scratch directory */

#define MREC_BYTES 14
#define MREC_POS_BYTES 6
#define MAX_PARTS 4096
#define KMER_MAX_EXACT 19       /* a k-mer of up to 19 digits is an exact 64-bit key */

typedef struct { uint8_t b[MREC_BYTES]; } mrec;

/* mrec_pack / mrec_key / mrec_pos - little-endian packing of (key, position) */
static inline void mrec_pack(mrec *r, uint64_t key, uint64_t pos)
{
    for (int i = 0; i < 8; i++) r->b[i] = (uint8_t)(key >> (8 * i));
    for (int i = 0; i < MREC_POS_BYTES; i++) r->b[8 + i] = (uint8_t)(pos >> (8 * i));
}
static inline uint64_t mrec_key(const mrec *r)
{
    uint64_t k = 0;
    for (int i = 7; i >= 0; i--) k = (k << 8) | r->b[i];
    return k;
}
static inline uint64_t mrec_pos(const mrec *r)
{
    uint64_t p = 0;
    for (int i = MREC_POS_BYTES - 1; i >= 0; i--) p = (p << 8) | r->b[8 + i];
    return p;
}

/* Rank used for minimizer selection: a random-looking order of k-mers so
 * that minimizers are spread uniformly (density ~2/(w+1)). */
static inline uint64_t kmer_rank(uint64_t value) { return mix64(value ^ 0xA5A5F00DBEEF1234ULL); }

/* Partition of a key: top bits of an independent hash of the value. */
static inline uint32_t key_partition(uint64_t value, int nparts)
{
    return (uint32_t)fastrange64(mix64(value + 0x1234567887654321ULL), (uint64_t)nparts);
}

/*
 * key_partition_pass - partition of a key when the scan is split into
 * npasses hash-range passes (-H pass/npasses): the hash is spread over
 * nparts * npasses global partitions and pass p (1-based) owns global
 * partitions [(p-1)*nparts, p*nparts).  Returns the local partition, or -1
 * if the key belongs to another pass.  With npasses == 1 this is exactly
 * key_partition(), so single-pass scratch sets are unchanged, and pass p of
 * R at P partitions writes files identical to partitions (p-1)*P .. p*P-1 of
 * a single scan with R*P partitions.
 */
static inline int key_partition_pass(uint64_t value, int nparts, int pass, int npasses)
{
    uint64_t g = fastrange64(mix64(value + 0x1234567887654321ULL), (uint64_t)nparts * (uint64_t)npasses);
    uint64_t lo = (uint64_t)(pass - 1) * (uint64_t)nparts;
    if (g < lo || g >= lo + (uint64_t)nparts) return -1;
    return (int)(g - lo);
}

/* ------------------------------------------------------------------ */
/* Scratch set: one or more directories holding the partition files    */
/* ------------------------------------------------------------------ */
/*
 * The first directory is the primary: it holds manifest, checkpoint, the
 * per-partition candidate files and the results.  Partition i lives in the
 * directory d with first_part[d] <= i < first_part[d+1]; pimatchscan chooses
 * the split in proportion to free space and records it in the manifest.
 */
#define MAX_DIRS 32
typedef struct {
    int ndirs;
    char *dirs[MAX_DIRS];
    int first_part[MAX_DIRS + 1];   /* cumulative partition counts; [ndirs] = nparts */
} scratch;

/*
 * scratch_parse - split a comma-separated directory list; every entry must
 * be an existing directory.  Partition ranges are left unset (all zero).
 */
static void scratch_parse(const char *spec, scratch *sc)
{
    memset(sc, 0, sizeof *sc);
    char *copy = strdup(spec);
    if (!copy) die("out of memory");
    for (char *tok = strtok(copy, ","); tok; tok = strtok(NULL, ",")) {
        while (*tok == ' ') tok++;
        size_t n = strlen(tok);
        while (n > 1 && (tok[n - 1] == ' ' || tok[n - 1] == '/')) tok[--n] = 0;
        if (!n) continue;
        if (sc->ndirs >= MAX_DIRS) die("more than %d scratch directories", MAX_DIRS);
        struct stat st;
        if (stat(tok, &st) != 0 || !S_ISDIR(st.st_mode)) die("scratch directory '%s' does not exist", tok);
        for (int i = 0; i < sc->ndirs; i++) if (!strcmp(sc->dirs[i], tok)) die("scratch directory '%s' given twice", tok);
        sc->dirs[sc->ndirs++] = strdup(tok);
    }
    if (!sc->ndirs) die("no scratch directory given");
}

/* scratch_dir_of - directory index holding partition i */
static int scratch_dir_of(const scratch *sc, int i)
{
    for (int d = 0; d < sc->ndirs; d++) if (i < sc->first_part[d + 1]) return d;
    die("internal: partition %d has no directory", i);
    return 0;
}

/*
 * scratch_assign - split nparts partitions over the directories in
 * proportion to their free space (each partition is expected to hold
 * per_part bytes; a 5% margin is kept).  Returns -1 (after a warning) when
 * the set cannot hold the expected total.  Directories that cannot fit a
 * single partition get none.
 */
static int scratch_assign(scratch *sc, int nparts, double per_part, const uint64_t *freeb)
{
    double total_free = 0;
    int cap[MAX_DIRS], n[MAX_DIRS], capsum = 0;
    for (int d = 0; d < sc->ndirs; d++) {
        double usable = (double)freeb[d] / 1.05;
        cap[d] = per_part > 0 ? (int)(usable / per_part) : nparts;
        if (cap[d] > nparts) cap[d] = nparts;
        capsum += cap[d];
        total_free += (double)freeb[d];
        n[d] = 0;
    }
    if (capsum < nparts) {
        char b1[32], b2[32];
        warn("not enough free space across the scratch set: need ~%s (+5%%), the %d director%s can hold %d of %d partitions (%s free in total)",
             fmt_bytes(per_part * nparts, b1, sizeof b1), sc->ndirs, sc->ndirs == 1 ? "y" : "ies", capsum, nparts,
             fmt_bytes(total_free, b2, sizeof b2));
        sc->first_part[0] = 0;                        /* nominal layout so a dry run can still print */
        for (int d = 0; d < sc->ndirs; d++) sc->first_part[d + 1] = d == 0 ? nparts : nparts;
        return -1;
    }
    int assigned = 0;
    for (int d = 0; d < sc->ndirs; d++) {
        n[d] = total_free > 0 ? (int)((double)nparts * (double)freeb[d] / total_free) : nparts / sc->ndirs;
        if (n[d] > cap[d]) n[d] = cap[d];
        assigned += n[d];
    }
    while (assigned < nparts) {                       /* hand out the rounding remainder to the roomiest */
        int best = -1; double bestroom = -1;
        for (int d = 0; d < sc->ndirs; d++) {
            if (n[d] >= cap[d]) continue;
            double room = (double)freeb[d] - (double)n[d] * per_part;
            if (room > bestroom) { bestroom = room; best = d; }
        }
        if (best < 0) die("internal: cannot place all partitions");
        n[best]++; assigned++;
    }
    sc->first_part[0] = 0;
    for (int d = 0; d < sc->ndirs; d++) sc->first_part[d + 1] = sc->first_part[d] + n[d];
    return 0;
}

/* ------------------------------------------------------------------ */
/* Manifest: written by pimatchscan, read by the other tools            */
/* ------------------------------------------------------------------ */
typedef struct {
    char source[1024];
    int L, k, w;
    uint64_t ndigits;       /* digits requested (UINT64_MAX = all)   */
    uint64_t first_pos;     /* first position of the source          */
    uint64_t end_pos;       /* last position + 1 actually scanned    */
    int nparts;
    char version[32];
    int complete;           /* scan finished normally                */
    int hpass, hpasses;     /* hash-range pass of the scan (1/1 = whole key space) */
    scratch sc;             /* directories and partition ranges      */
} manifest;

/*
 * primary_dir - the first directory of a -d list (the one holding the
 * manifest); returns a malloc'd copy with any ",more,dirs" removed.
 */
static char *primary_dir(const char *spec)
{
    const char *c = strchr(spec, ',');
    size_t n = c ? (size_t)(c - spec) : strlen(spec);
    while (n > 1 && spec[n - 1] == '/') n--;
    char *p = malloc(n + 1);
    if (!p) die("out of memory");
    memcpy(p, spec, n); p[n] = 0;
    return p;
}

static void manifest_path(const char *dir, char *out, size_t n) { snprintf(out, n, "%s/manifest.txt", dir); }
static void checkpoint_path(const char *dir, char *out, size_t n) { snprintf(out, n, "%s/checkpoint.txt", dir); }
static void part_path(const scratch *sc, int i, char *out, size_t n) { snprintf(out, n, "%s/part_%04d.bin", sc->dirs[scratch_dir_of(sc, i)], i); }
static void cand_path(const char *dir, int i, char *out, size_t n) { snprintf(out, n, "%s/cand_%04d.txt", dir, i); }
static void done_path(const char *dir, int i, char *out, size_t n) { snprintf(out, n, "%s/cand_%04d.done", dir, i); }

/* manifest_write - atomically (temp file + rename) write the manifest. */
static void manifest_write(const char *dir, const manifest *m)
{
    char path[2048], tmp[2048 + 8];      /* room for the ".tmp" suffix */
    manifest_path(dir, path, sizeof path);
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) die("cannot write '%s': %s", tmp, strerror(errno));
    fprintf(f, "version=%s\nsource=%s\nL=%d\nk=%d\nw=%d\nndigits=%" PRIu64 "\nfirst_pos=%" PRIu64
               "\nend_pos=%" PRIu64 "\nnparts=%d\nrecord_bytes=%d\ncomplete=%d\n",
            m->version, m->source, m->L, m->k, m->w, m->ndigits, m->first_pos, m->end_pos, m->nparts,
            MREC_BYTES, m->complete);
    fprintf(f, "hash_pass=%d\nhash_passes=%d\nndirs=%d\n", m->hpass > 0 ? m->hpass : 1, m->hpasses > 0 ? m->hpasses : 1, m->sc.ndirs);
    for (int d = 0; d < m->sc.ndirs; d++)
        fprintf(f, "dir%d=%s\ndir%d_parts=%d\n", d, m->sc.dirs[d], d, m->sc.first_part[d + 1] - m->sc.first_part[d]);
    if (fclose(f)) die("cannot write '%s': %s", tmp, strerror(errno));
    if (rename(tmp, path)) die("cannot rename '%s': %s", tmp, strerror(errno));
}

/* manifest_read - parse the manifest; dies if missing or malformed. */
static void manifest_read(const char *dir, manifest *m)
{
    char path[2048], line[1200];
    manifest_path(dir, path, sizeof path);
    FILE *f = fopen(path, "r");
    if (!f) die("no manifest in '%s' (%s) - run pimatchscan first", dir, strerror(errno));
    memset(m, 0, sizeof *m);
    int rb = 0, ndirs_seen = 0, parts_seen[MAX_DIRS];
    memset(parts_seen, 0, sizeof parts_seen);
    while (fgets(line, sizeof line, f)) {
        char *nl = strchr(line, '\n'); if (nl) *nl = 0;
        char *eq = strchr(line, '='); if (!eq) continue;
        *eq = 0; const char *v = eq + 1;
        if (!strcmp(line, "version")) snprintf(m->version, sizeof m->version, "%s", v);
        else if (!strcmp(line, "source")) snprintf(m->source, sizeof m->source, "%s", v);
        else if (!strcmp(line, "L")) m->L = atoi(v);
        else if (!strcmp(line, "k")) m->k = atoi(v);
        else if (!strcmp(line, "w")) m->w = atoi(v);
        else if (!strcmp(line, "ndigits")) m->ndigits = strtoull(v, NULL, 10);
        else if (!strcmp(line, "first_pos")) m->first_pos = strtoull(v, NULL, 10);
        else if (!strcmp(line, "end_pos")) m->end_pos = strtoull(v, NULL, 10);
        else if (!strcmp(line, "nparts")) m->nparts = atoi(v);
        else if (!strcmp(line, "record_bytes")) rb = atoi(v);
        else if (!strcmp(line, "complete")) m->complete = atoi(v);
        else if (!strcmp(line, "hash_pass")) m->hpass = atoi(v);
        else if (!strcmp(line, "hash_passes")) m->hpasses = atoi(v);
        else if (!strcmp(line, "ndirs")) ndirs_seen = atoi(v);
        else if (!strncmp(line, "dir", 3)) {
            char *e; long d = strtol(line + 3, &e, 10);
            if (d < 0 || d >= MAX_DIRS) die("manifest in '%s': bad directory index in '%s'", dir, line);
            if (!*e) { m->sc.dirs[d] = strdup(v); if ((int)d + 1 > m->sc.ndirs) m->sc.ndirs = (int)d + 1; }
            else if (!strcmp(e, "_parts")) parts_seen[d] = atoi(v);
        }
    }
    fclose(f);
    if (m->nparts < 1 || m->nparts > MAX_PARTS || m->k < 1 || m->w < 1 || m->L != m->k + m->w - 1)
        die("manifest in '%s' is malformed", dir);
    if (rb != MREC_BYTES) die("manifest record size %d does not match this build (%d)", rb, MREC_BYTES);
    if (m->hpass < 1) m->hpass = 1;
    if (m->hpasses < 1) m->hpasses = 1;
    if (m->hpass > m->hpasses) die("manifest in '%s': hash pass %d of %d", dir, m->hpass, m->hpasses);
    if (ndirs_seen == 0) {                              /* pre-1.6 manifest: one directory, all partitions */
        m->sc.ndirs = 1; m->sc.dirs[0] = strdup(dir);
        m->sc.first_part[0] = 0; m->sc.first_part[1] = m->nparts;
    } else {
        if (ndirs_seen != m->sc.ndirs) die("manifest in '%s' lists %d directories but ndirs=%d", dir, m->sc.ndirs, ndirs_seen);
        m->sc.first_part[0] = 0;
        for (int d = 0; d < m->sc.ndirs; d++) {
            if (!m->sc.dirs[d]) die("manifest in '%s': dir%d missing", dir, d);
            m->sc.first_part[d + 1] = m->sc.first_part[d] + parts_seen[d];
        }
        if (m->sc.first_part[m->sc.ndirs] != m->nparts) die("manifest in '%s': partition counts per directory do not add up to %d", dir, m->nparts);
        if (strcmp(m->sc.dirs[0], dir)) warn("manifest names '%s' as its primary directory, reading it from '%s'", m->sc.dirs[0], dir);
    }
}

/*
 * manifest_check_dirs - verify every scratch directory of a manifest exists
 * (partition files may legitimately be gone after pimatchsort -D).
 */
static void manifest_check_dirs(const manifest *m)
{
    for (int d = 0; d < m->sc.ndirs; d++) {
        struct stat st;
        if (stat(m->sc.dirs[d], &st) != 0 || !S_ISDIR(st.st_mode))
            die("scratch directory '%s' (from the manifest) is not accessible", m->sc.dirs[d]);
    }
}

/* ------------------------------------------------------------------ */
/* Checkpoint (pimatchscan): resume position + partition file lengths  */
/* ------------------------------------------------------------------ */
typedef struct {
    uint64_t next_pos;      /* first position not yet fully processed */
    uint64_t records;       /* records written so far (for statistics) */
    int nparts;
    uint64_t *sizes;        /* byte length of each partition file      */
} checkpoint;

static void checkpoint_write(const char *dir, const checkpoint *c)
{
    char path[2048], tmp[2048 + 8];      /* room for the ".tmp" suffix */
    checkpoint_path(dir, path, sizeof path);
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) die("cannot write '%s': %s", tmp, strerror(errno));
    fprintf(f, "next_pos=%" PRIu64 "\nrecords=%" PRIu64 "\nnparts=%d\n", c->next_pos, c->records, c->nparts);
    for (int i = 0; i < c->nparts; i++) fprintf(f, "%" PRIu64 "\n", c->sizes[i]);
    if (fflush(f) || fsync(fileno(f)) || fclose(f)) die("cannot write '%s': %s", tmp, strerror(errno));
    if (rename(tmp, path)) die("cannot rename '%s': %s", tmp, strerror(errno));
}

/* returns 0 if no checkpoint exists */
static int checkpoint_read(const char *dir, checkpoint *c, int nparts)
{
    char path[2048], line[64];
    checkpoint_path(dir, path, sizeof path);
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    memset(c, 0, sizeof *c);
    if (!fgets(line, sizeof line, f) || sscanf(line, "next_pos=%" SCNu64, &c->next_pos) != 1) die("bad checkpoint");
    if (!fgets(line, sizeof line, f) || sscanf(line, "records=%" SCNu64, &c->records) != 1) die("bad checkpoint");
    if (!fgets(line, sizeof line, f) || sscanf(line, "nparts=%d", &c->nparts) != 1) die("bad checkpoint");
    if (c->nparts != nparts) die("checkpoint has %d partitions, manifest %d", c->nparts, nparts);
    c->sizes = calloc((size_t)nparts, sizeof *c->sizes);
    if (!c->sizes) die("out of memory");
    for (int i = 0; i < nparts; i++)
        if (!fgets(line, sizeof line, f) || sscanf(line, "%" SCNu64, &c->sizes[i]) != 1) die("bad checkpoint (partition %d)", i);
    fclose(f);
    return 1;
}

/* free bytes on the filesystem holding dir (0 if unknown) */
static uint64_t dir_free_bytes(const char *dir)
{
    struct statvfs sv;
    if (statvfs(dir, &sv) != 0) return 0;
    return (uint64_t)sv.f_bavail * (uint64_t)sv.f_frsize;
}

static uint64_t file_size_of(const char *path)
{
    struct stat st;
    if (stat(path, &st) != 0) return 0;
    return (uint64_t)st.st_size;
}

/* fmt helpers shared by the tools */
static const char *fmt_pct(double x, char *buf, size_t n) { snprintf(buf, n, "%.3g%%", x * 100.0); return buf; }

#endif /* PIMATCH_H */
