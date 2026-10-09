/*
 * ycd.c - read and write y-cruncher compressed digit files (.ycd)
 * See ycd.h for the API and ycd.md for the format.
 */
#define _POSIX_C_SOURCE 200809L
#define _DEFAULT_SOURCE             /* strcasecmp, strdup on glibc */
#define _DARWIN_C_SOURCE
#include "ycd.h"

#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

const uint8_t ycd_pairs[200] = {
    0,0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9, 1,0,1,1,1,2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,
    2,0,2,1,2,2,2,3,2,4,2,5,2,6,2,7,2,8,2,9, 3,0,3,1,3,2,3,3,3,4,3,5,3,6,3,7,3,8,3,9,
    4,0,4,1,4,2,4,3,4,4,4,5,4,6,4,7,4,8,4,9, 5,0,5,1,5,2,5,3,5,4,5,5,5,6,5,7,5,8,5,9,
    6,0,6,1,6,2,6,3,6,4,6,5,6,6,6,7,6,8,6,9, 7,0,7,1,7,2,7,3,7,4,7,5,7,6,7,7,7,8,7,9,
    8,0,8,1,8,2,8,3,8,4,8,5,8,6,8,7,8,8,8,9, 9,0,9,1,9,2,9,3,9,4,9,5,9,6,9,7,9,8,9,9 };

#define HEADER_MAX 16384            /* a header is ~200 bytes; this bounds the search for its NUL */

/* set_err - printf into an optional error buffer; always returns -1 */
static int set_err(char *err, size_t errlen, const char *fmt, ...)
{
    if (err && errlen) {
        va_list ap;
        va_start(ap, fmt);
        vsnprintf(err, errlen, fmt, ap);
        va_end(ap);
    }
    return -1;
}

/* ================================================================ */
/* Header                                                            */
/* ================================================================ */

int ycd_is_ycd(const char *path)
{
    int fd = open(path, O_RDONLY);
    if (fd < 0) return -1;
    char head[sizeof YCD_MAGIC];
    ssize_t n = read(fd, head, sizeof YCD_MAGIC - 1);
    close(fd);
    if (n < 0) return -1;
    return n == (ssize_t)(sizeof YCD_MAGIC - 1) && !memcmp(head, YCD_MAGIC, sizeof YCD_MAGIC - 1);
}

void ycd_file_free(ycd_file *f)
{
    if (f) { free(f->path); f->path = NULL; }
}

int ycd_read_header(const char *path, ycd_file *f, char *err, size_t errlen)
{
    memset(f, 0, sizeof *f);
    int fd = open(path, O_RDONLY);
    if (fd < 0) return set_err(err, errlen, "cannot open '%s': %s", path, strerror(errno));
    char *hdr = malloc(HEADER_MAX + 1);
    if (!hdr) { close(fd); return set_err(err, errlen, "out of memory"); }
    ssize_t n = 0;
    while (n < HEADER_MAX) {                        /* short reads are legal on pipes/NFS */
        ssize_t r = read(fd, hdr + n, (size_t)(HEADER_MAX - n));
        if (r < 0) { if (errno == EINTR) continue; break; }
        if (r == 0) break;
        n += r;
    }
    struct stat st;
    int st_ok = fstat(fd, &st) == 0;
    close(fd);
    if (n <= 0 || !st_ok) { free(hdr); return set_err(err, errlen, "cannot read '%s'", path); }
    hdr[n] = 0;
    if ((size_t)n < sizeof YCD_MAGIC - 1 || memcmp(hdr, YCD_MAGIC, sizeof YCD_MAGIC - 1)) {
        free(hdr);
        return set_err(err, errlen, "'%s' is not a .ycd file (no \"%s\" line)", path, YCD_MAGIC);
    }
    ssize_t nul = -1;
    for (ssize_t i = 0; i < n; i++) if (hdr[i] == 0) { nul = i; break; }
    if (nul < 0) { free(hdr); return set_err(err, errlen, "'%s': header not terminated by a NUL within %d bytes", path, HEADER_MAX); }

    int seen_end = 0;
    char *p = hdr;
    while (p && *p) {
        char *line = p;
        char *nl = strchr(p, '\n');
        if (nl) { *nl = 0; p = nl + 1; } else p = NULL;
        size_t len = strlen(line);
        while (len && (line[len - 1] == '\r' || line[len - 1] == ' ' || line[len - 1] == '\t')) line[--len] = 0;
        if (!len) continue;
        if (!strcmp(line, "EndHeader")) { seen_end = 1; break; }
        char *colon = strchr(line, ':');
        if (!colon) continue;                       /* the magic line, or anything unknown */
        *colon = 0;
        char *val = colon + 1;
        while (*val == '\t' || *val == ' ') val++;
        if (!strcasecmp(line, "FileVersion")) snprintf(f->version, sizeof f->version, "%s", val);
        else if (!strcasecmp(line, "Base")) f->base = atoi(val);
        else if (!strcasecmp(line, "FirstDigits")) snprintf(f->first_digits, sizeof f->first_digits, "%s", val);
        else if (!strcasecmp(line, "TotalDigits")) f->total_digits = strtoull(val, NULL, 10);
        else if (!strcasecmp(line, "Blocksize")) f->blocksize = strtoull(val, NULL, 10);
        else if (!strcasecmp(line, "BlockID")) f->blockid = strtoull(val, NULL, 10);
    }
    free(hdr);
    if (!seen_end) return set_err(err, errlen, "'%s': header has no EndHeader line", path);
    if (strcmp(f->version, "1.0.0") && strcmp(f->version, "1.1.0"))
        return set_err(err, errlen, "'%s': unsupported FileVersion '%s' (1.0.0 and 1.1.0 are known)", path, f->version);
    if (f->base == 10) f->digits_per_word = YCD_DEC_DIGITS_PER_WORD;
    else if (f->base == 16) f->digits_per_word = YCD_HEX_DIGITS_PER_WORD;
    else return set_err(err, errlen, "'%s': unsupported Base %d (10 and 16 are known)", path, f->base);
    if (f->blocksize == 0) return set_err(err, errlen, "'%s': missing or zero Blocksize", path);
    if (f->blockid > UINT64_MAX / f->blocksize - 1) return set_err(err, errlen, "'%s': BlockID x Blocksize overflows", path);

    f->data_offset = (uint64_t)nul + 1;
    f->file_size = (uint64_t)st.st_size;
    f->start = f->blockid * f->blocksize;
    uint64_t end = f->start + f->blocksize;
    if (f->total_digits) {
        if (f->total_digits <= f->start) return set_err(err, errlen, "'%s': BlockID %" PRIu64 " lies beyond TotalDigits %" PRIu64, path, f->blockid, f->total_digits);
        if (end > f->total_digits) end = f->total_digits;
    }
    f->ndigits = end - f->start;

    /* files are longer than their data, so only a shortfall is detectable */
    uint64_t have_words = f->file_size > f->data_offset ? (f->file_size - f->data_offset) / 8 : 0;
    uint64_t need_words = (f->ndigits + (uint64_t)f->digits_per_word - 1) / (uint64_t)f->digits_per_word;
    if (have_words < need_words) {
        if (f->total_digits)
            return set_err(err, errlen, "'%s' is truncated: %" PRIu64 " words present, %" PRIu64 " needed", path, have_words, need_words);
        /* TotalDigits unknown: accept a short file as the final block of a set */
        f->ndigits = have_words * (uint64_t)f->digits_per_word;
        f->length_uncertain = 1;
        if (f->ndigits == 0) return set_err(err, errlen, "'%s' holds no digit words", path);
    }
    f->path = strdup(path);
    if (!f->path) return set_err(err, errlen, "out of memory");
    return 0;
}

/* ================================================================ */
/* Set                                                               */
/* ================================================================ */

struct ycd_set {
    ycd_file *files;
    int nfiles;
    int base, dpw;
    int int_digits;             /* integer-part digits carried by the stream */
    int convention_checked;
    uint64_t first_pos, last_pos;
};

static int cmp_block(const void *a, const void *b)
{
    const ycd_file *x = a, *y = b;
    return x->blockid < y->blockid ? -1 : x->blockid > y->blockid;
}

/* add_path - append one path to a growing list */
static int add_path(char ***list, int *n, int *cap, const char *path)
{
    if (*n == *cap) {
        int nc = *cap ? *cap * 2 : 64;
        char **nl = realloc(*list, (size_t)nc * sizeof **list);
        if (!nl) return -1;
        *list = nl; *cap = nc;
    }
    (*list)[*n] = strdup(path);
    if (!(*list)[*n]) return -1;
    (*n)++;
    return 0;
}

/*
 * check_convention - compare block 0's first word with the FirstDigits
 * header to learn whether the stream starts after the radix point (the
 * y-cruncher layout) or carries the integer part as well.
 */
static void check_convention(ycd_set *s)
{
    s->int_digits = 0;
    s->convention_checked = 0;
    const ycd_file *f = &s->files[0];
    if (f->blockid != 0 || !f->first_digits[0] || f->ndigits < (uint64_t)s->dpw) return;
    const char *dot = strchr(f->first_digits, '.');
    if (!dot) return;
    size_t ilen = (size_t)(dot - f->first_digits);
    if (ilen == 0 || ilen >= (size_t)s->dpw || strlen(dot + 1) < (size_t)s->dpw) return;

    uint8_t raw[8];
    int fd = open(f->path, O_RDONLY);
    if (fd < 0) return;
    ssize_t n = pread(fd, raw, 8, (off_t)f->data_offset);
    close(fd);
    if (n != 8) return;
    uint8_t got[YCD_MAX_DIGITS_PER_WORD];
    if (s->base == 10) ycd_decode_dec(ycd_load_le64(raw), got); else ycd_decode_hex(ycd_load_le64(raw), got);

    /* expected digit values for the two layouts */
    char frac[YCD_MAX_DIGITS_PER_WORD + 1], with_int[YCD_MAX_DIGITS_PER_WORD + 1];
    memcpy(frac, dot + 1, (size_t)s->dpw); frac[s->dpw] = 0;
    memcpy(with_int, f->first_digits, ilen);
    memcpy(with_int + ilen, dot + 1, (size_t)s->dpw - ilen); with_int[s->dpw] = 0;
    int m_frac = 1, m_int = 1;
    for (int i = 0; i < s->dpw; i++) {
        char c = (char)(got[i] < 10 ? '0' + got[i] : 'a' + got[i] - 10);
        char a = frac[i], b = with_int[i];
        if (a >= 'A' && a <= 'F') a = (char)(a - 'A' + 'a');
        if (b >= 'A' && b <= 'F') b = (char)(b - 'A' + 'a');
        if (c != a) m_frac = 0;
        if (c != b) m_int = 0;
    }
    if (m_frac) { s->int_digits = 0; s->convention_checked = 1; }
    else if (m_int) { s->int_digits = (int)ilen; s->convention_checked = 1; }
}

ycd_set *ycd_open_paths(const char *const *paths, int npaths, char *err, size_t errlen)
{
    char **list = NULL; int n = 0, cap = 0;
    ycd_set *s = NULL;
    if (npaths <= 0) { set_err(err, errlen, "no .ycd files given"); return NULL; }

    for (int i = 0; i < npaths; i++) {
        struct stat st;
        if (stat(paths[i], &st) != 0) { set_err(err, errlen, "cannot stat '%s': %s", paths[i], strerror(errno)); goto fail; }
        if (S_ISDIR(st.st_mode)) {
            DIR *d = opendir(paths[i]);
            if (!d) { set_err(err, errlen, "cannot open directory '%s': %s", paths[i], strerror(errno)); goto fail; }
            struct dirent *de; int found = 0;
            while ((de = readdir(d)) != NULL) {
                size_t l = strlen(de->d_name);
                if (l < 5 || strcasecmp(de->d_name + l - 4, ".ycd")) continue;
                size_t pl = strlen(paths[i]) + 1 + l + 1;
                char *p = malloc(pl);
                if (!p) { closedir(d); set_err(err, errlen, "out of memory"); goto fail; }
                snprintf(p, pl, "%s/%s", paths[i], de->d_name);
                int rc = add_path(&list, &n, &cap, p);
                free(p);
                if (rc) { closedir(d); set_err(err, errlen, "out of memory"); goto fail; }
                found++;
            }
            closedir(d);
            if (!found) { set_err(err, errlen, "directory '%s' contains no .ycd files", paths[i]); goto fail; }
        } else if (add_path(&list, &n, &cap, paths[i])) { set_err(err, errlen, "out of memory"); goto fail; }
    }

    s = calloc(1, sizeof *s);
    if (s) s->files = calloc((size_t)n, sizeof *s->files);
    if (!s || !s->files) { set_err(err, errlen, "out of memory"); goto fail; }
    for (int i = 0; i < n; i++) {
        if (ycd_read_header(list[i], &s->files[s->nfiles], err, errlen)) goto fail;
        s->nfiles++;
    }
    qsort(s->files, (size_t)s->nfiles, sizeof *s->files, cmp_block);
    s->base = s->files[0].base; s->dpw = s->files[0].digits_per_word;
    for (int i = 1; i < s->nfiles; i++) {
        const ycd_file *a = &s->files[i - 1], *b = &s->files[i];
        if (b->base != a->base) { set_err(err, errlen, "'%s' and '%s' have different Base", a->path, b->path); goto fail; }
        if (b->blocksize != a->blocksize) { set_err(err, errlen, "'%s' and '%s' have different Blocksize", a->path, b->path); goto fail; }
        if (b->blockid == a->blockid) { set_err(err, errlen, "duplicate BlockID %" PRIu64 ": '%s' and '%s'", a->blockid, a->path, b->path); goto fail; }
        if (b->blockid != a->blockid + 1) {
            set_err(err, errlen, "blocks are not contiguous: BlockID %" PRIu64 " ('%s') is followed by %" PRIu64 " ('%s')",
                    a->blockid, a->path, b->blockid, b->path);
            goto fail;
        }
        if (a->ndigits != a->blocksize) { set_err(err, errlen, "'%s' is a short block but is not the last of the set", a->path); goto fail; }
    }
    check_convention(s);
    {
        const ycd_file *f0 = &s->files[0], *fl = &s->files[s->nfiles - 1];
        uint64_t idg = (uint64_t)s->int_digits;
        /* stream offset o holds position o + 1 - int_digits; offsets below int_digits are the integer part */
        s->first_pos = f0->start >= idg ? f0->start + 1 - idg : 1;
        uint64_t end_off = fl->start + fl->ndigits;                 /* one past the last digit */
        if (end_off <= idg) { set_err(err, errlen, "the set holds no fractional digits"); goto fail; }
        s->last_pos = end_off - idg;
    }
    for (int i = 0; i < n; i++) free(list[i]);
    free(list);
    return s;

fail:
    for (int i = 0; i < n; i++) free(list[i]);
    free(list);
    ycd_close(s);
    return NULL;
}

ycd_set *ycd_open(const char *spec, char *err, size_t errlen)
{
    char *copy = strdup(spec ? spec : "");
    if (!copy) { set_err(err, errlen, "out of memory"); return NULL; }
    const char **paths = NULL; int n = 0, cap = 0;
    for (char *tok = copy, *next; tok; tok = next) {
        char *comma = strchr(tok, ',');
        if (comma) { *comma = 0; next = comma + 1; } else next = NULL;
        if (!*tok) continue;
        if (n == cap) {
            cap = cap ? cap * 2 : 16;
            const char **np = realloc((void *)paths, (size_t)cap * sizeof *paths);
            if (!np) { free((void *)paths); free(copy); set_err(err, errlen, "out of memory"); return NULL; }
            paths = np;
        }
        paths[n++] = tok;
    }
    ycd_set *s = ycd_open_paths(paths, n, err, errlen);
    free((void *)paths);
    free(copy);
    return s;
}

void ycd_close(ycd_set *s)
{
    if (!s) return;
    for (int i = 0; i < s->nfiles; i++) ycd_file_free(&s->files[i]);
    free(s->files);
    free(s);
}

int ycd_nfiles(const ycd_set *s) { return s->nfiles; }
const ycd_file *ycd_get_file(const ycd_set *s, int i) { return i >= 0 && i < s->nfiles ? &s->files[i] : NULL; }
int ycd_base(const ycd_set *s) { return s->base; }
uint64_t ycd_first_pos(const ycd_set *s) { return s->first_pos; }
uint64_t ycd_last_pos(const ycd_set *s) { return s->last_pos; }
uint64_t ycd_ndigits(const ycd_set *s) { return s->last_pos - s->first_pos + 1; }
int ycd_integer_digits(const ycd_set *s) { return s->int_digits; }
int ycd_convention_checked(const ycd_set *s) { return s->convention_checked; }
int ycd_length_uncertain(const ycd_set *s) { return s->files[s->nfiles - 1].length_uncertain; }

/* ================================================================ */
/* Reader                                                            */
/* ================================================================ */

struct ycd_reader {
    const ycd_set *s;
    int flags;
    int fd, cur;                /* open file and its index (-1 = none) */
    uint64_t pos;               /* position of the next digit to deliver */
    uint64_t word;              /* index of the next word to decode in the current file */
    uint64_t file_left;         /* digits of the current file not yet decoded (from word's first digit) */
    int drop;                   /* digits of the next word to discard (after a seek) */
    uint8_t pend[YCD_MAX_DIGITS_PER_WORD];
    int pend_off, pend_len;     /* decoded digits not yet delivered */
    uint8_t *raw; size_t rawcap, rawlen, rawpos;   /* buffered words of the current file */
    int eof;
    char err[256];
};

ycd_reader *ycd_reader_new(const ycd_set *s, int flags)
{
    ycd_reader *r = calloc(1, sizeof *r);
    if (!r) return NULL;
    r->s = s; r->flags = flags; r->fd = -1; r->cur = -1;
    r->rawcap = 1u << 20;
    if (ycd_reader_seek(r, s->first_pos)) { free(r); return NULL; }
    return r;
}

void ycd_reader_free(ycd_reader *r)
{
    if (!r) return;
    if (r->fd >= 0) close(r->fd);
    free(r->raw);
    free(r);
}

int ycd_reader_set_buffer(ycd_reader *r, size_t bytes)
{
    if (bytes < 64) bytes = 64;
    bytes &= ~(size_t)7;
    if (r->rawpos < r->rawlen) {                    /* keep position: drop buffered words, re-read later */
        r->rawlen = r->rawpos = 0;
    }
    free(r->raw); r->raw = NULL;
    r->rawcap = bytes;
    return 0;
}

const char *ycd_reader_error(const ycd_reader *r) { return r->err; }
uint64_t ycd_reader_tell(const ycd_reader *r) { return r->pos; }

/* use_file - make file idx current, positioned at word index `word` */
static int use_file(ycd_reader *r, int idx, uint64_t word)
{
    const ycd_file *f = &r->s->files[idx];
    if (r->cur != idx) {
        if (r->fd >= 0) close(r->fd);
        r->fd = open(f->path, O_RDONLY);
        if (r->fd < 0) { r->cur = -1; return set_err(r->err, sizeof r->err, "cannot open '%s': %s", f->path, strerror(errno)); }
        r->cur = idx;
    }
    r->word = word;
    r->file_left = f->ndigits - word * (uint64_t)r->s->dpw;
    r->rawlen = r->rawpos = 0;
    return 0;
}

int ycd_reader_seek(ycd_reader *r, uint64_t pos)
{
    const ycd_set *s = r->s;
    if (pos < s->first_pos || pos > s->last_pos + 1)
        return set_err(r->err, sizeof r->err, "position %" PRIu64 " is outside the set (%" PRIu64 "..%" PRIu64 ")", pos, s->first_pos, s->last_pos);
    r->pos = pos; r->pend_len = r->pend_off = 0; r->drop = 0; r->eof = 0;
    if (pos == s->last_pos + 1) { r->eof = 1; return 0; }
    uint64_t o = pos - 1 + (uint64_t)s->int_digits;         /* 0-based stream offset */
    int lo = 0, hi = s->nfiles - 1;
    while (lo < hi) {                                        /* last file with start <= o */
        int mid = (lo + hi + 1) / 2;
        if (s->files[mid].start <= o) lo = mid; else hi = mid - 1;
    }
    uint64_t local = o - s->files[lo].start;
    if (use_file(r, lo, local / (uint64_t)s->dpw)) return -1;
    r->drop = (int)(local % (uint64_t)s->dpw);
    return 0;
}

/* fill_raw - buffer more words of the current file; -1 on error or short file */
static int fill_raw(ycd_reader *r)
{
    const ycd_file *f = &r->s->files[r->cur];
    if (!r->raw) {
        r->raw = malloc(r->rawcap);
        if (!r->raw) return set_err(r->err, sizeof r->err, "out of memory");
    }
    uint64_t words_left = (r->file_left + (uint64_t)r->s->dpw - 1) / (uint64_t)r->s->dpw;
    size_t want = r->rawcap;
    if (words_left * 8 < want) want = (size_t)(words_left * 8);
    size_t got = 0;
    uint64_t off = f->data_offset + r->word * 8;
    while (got < want) {
        ssize_t n = pread(r->fd, r->raw + got, want - got, (off_t)(off + got));
        if (n < 0) { if (errno == EINTR) continue; return set_err(r->err, sizeof r->err, "read error on '%s': %s", f->path, strerror(errno)); }
        if (n == 0) break;
        got += (size_t)n;
    }
    got &= ~(size_t)7;
    if (got == 0) return set_err(r->err, sizeof r->err, "'%s' ended before its block was complete", f->path);
    r->rawlen = got; r->rawpos = 0;
    return 0;
}

int64_t ycd_reader_read(ycd_reader *r, uint8_t *dst, size_t max)
{
    const ycd_set *s = r->s;
    const int dpw = s->dpw, dec = s->base == 10;
    size_t produced = 0;
    while (produced < max) {
        if (r->pend_len) {                                   /* leftovers of a split word */
            size_t take = (size_t)r->pend_len < max - produced ? (size_t)r->pend_len : max - produced;
            memcpy(dst + produced, r->pend + r->pend_off, take);
            produced += take; r->pend_off += (int)take; r->pend_len -= (int)take;
            continue;
        }
        if (r->eof) break;
        if (r->file_left == 0) {                             /* next block */
            if (r->cur + 1 >= s->nfiles) { r->eof = 1; break; }
            if (use_file(r, r->cur + 1, 0)) return -1;
        }
        if (r->rawpos >= r->rawlen && fill_raw(r)) return -1;

        /* bulk path: whole words straight into dst */
        if (r->drop == 0) {
            size_t nw = (r->rawlen - r->rawpos) / 8;
            uint64_t full = r->file_left / (uint64_t)dpw;
            if (nw > full) nw = (size_t)full;
            size_t room = (max - produced) / (size_t)dpw;
            if (nw > room) nw = room;
            const uint8_t *p = r->raw + r->rawpos;
            uint8_t *o = dst + produced;
            for (size_t i = 0; i < nw; i++, p += 8, o += dpw) {
                uint64_t w = ycd_load_le64(p);
                if (dec) {
                    if (w >= YCD_DEC_WORD_LIMIT)
                        return set_err(r->err, sizeof r->err, "'%s': word %" PRIu64 " is not a 19-digit value - corrupt file or wrong format",
                                       s->files[r->cur].path, r->word + i);
                    ycd_decode_dec(w, o);
                } else ycd_decode_hex(w, o);
            }
            if (nw) {
                r->rawpos += nw * 8; r->word += nw;
                r->file_left -= (uint64_t)nw * (uint64_t)dpw;
                produced += nw * (size_t)dpw;
                continue;
            }
        }
        /* one word through the pending buffer: after a seek, at a block's
         * partial last word, or when fewer than a word's digits are wanted */
        uint64_t w = ycd_load_le64(r->raw + r->rawpos);
        if (dec && w >= YCD_DEC_WORD_LIMIT)
            return set_err(r->err, sizeof r->err, "'%s': word %" PRIu64 " is not a 19-digit value - corrupt file or wrong format",
                           s->files[r->cur].path, r->word);
        if (dec) ycd_decode_dec(w, r->pend); else ycd_decode_hex(w, r->pend);
        int valid = r->file_left < (uint64_t)dpw ? (int)r->file_left : dpw;   /* partial last word: leading digits count */
        r->rawpos += 8; r->word++;
        r->file_left -= (uint64_t)valid;
        r->pend_off = r->drop < valid ? r->drop : valid;
        r->pend_len = valid - r->pend_off;
        r->drop = 0;
    }
    if (r->flags & YCD_ASCII)
        for (size_t i = 0; i < produced; i++) dst[i] = (uint8_t)(dst[i] < 10 ? '0' + dst[i] : 'a' + dst[i] - 10);
    r->pos += produced;
    return (int64_t)produced;
}

int64_t ycd_reader_read_at(ycd_reader *r, uint64_t pos, uint8_t *dst, size_t n)
{
    if (ycd_reader_seek(r, pos)) return -1;
    return ycd_reader_read(r, dst, n);
}

/* ================================================================ */
/* Writer                                                            */
/* ================================================================ */

struct ycd_writer {
    FILE *f;
    int base, dpw, flags;
    uint64_t limit, written;        /* digits allowed in this block / appended so far */
    int total_known;                /* TotalDigits was recorded */
    uint8_t word[YCD_MAX_DIGITS_PER_WORD];
    int nword;                      /* digits collected for the next word */
    uint64_t bytes;                 /* bytes written to the file */
    char err[256];
};

ycd_writer *ycd_writer_create(const char *path, const ycd_write_opts *o, int flags, char *err, size_t errlen)
{
    int base = o->base ? o->base : 10;
    if (base != 10 && base != 16) { set_err(err, errlen, "unsupported base %d", base); return NULL; }
    if (o->blocksize == 0) { set_err(err, errlen, "blocksize must be > 0"); return NULL; }
    uint64_t start = o->blockid * o->blocksize, limit = o->blocksize;
    if (o->total_digits) {
        if (o->total_digits <= start) { set_err(err, errlen, "block %" PRIu64 " lies beyond total_digits", o->blockid); return NULL; }
        if (start + limit > o->total_digits) limit = o->total_digits - start;
    }
    ycd_writer *w = calloc(1, sizeof *w);
    if (!w) { set_err(err, errlen, "out of memory"); return NULL; }
    w->f = fopen(path, "wb");
    if (!w->f) { set_err(err, errlen, "cannot create '%s': %s", path, strerror(errno)); free(w); return NULL; }
    setvbuf(w->f, NULL, _IOFBF, 1 << 20);
    w->base = base; w->dpw = base == 10 ? YCD_DEC_DIGITS_PER_WORD : YCD_HEX_DIGITS_PER_WORD;
    w->flags = flags; w->limit = limit; w->total_known = o->total_digits != 0;
    /* the exact layout y-cruncher writes: CRLF lines, a tab after each colon, NUL after the last CRLF */
    int n = fprintf(w->f, "%s\r\n\r\nFileVersion:\t%s\r\n\r\nBase:\t%d\r\n\r\n", YCD_MAGIC, o->version ? o->version : "1.1.0", base);
    if (o->first_digits && *o->first_digits) n += fprintf(w->f, "FirstDigits:\t%s\r\n\r\n", o->first_digits);
    n += fprintf(w->f, "TotalDigits:\t%" PRIu64 "\r\n\r\nBlocksize:\t%" PRIu64 "\r\nBlockID:\t%" PRIu64 "\r\n\r\nEndHeader\r\n\r\n",
                 o->total_digits, o->blocksize, o->blockid);
    fputc(0, w->f);
    w->bytes = (uint64_t)n + 1;
    if (ferror(w->f)) { set_err(err, errlen, "write error on '%s'", path); fclose(w->f); free(w); return NULL; }
    return w;
}

const char *ycd_writer_error(const ycd_writer *w) { return w->err; }

/* flush_word - emit the collected digits as one word (short words are left-aligned) */
static void flush_word(ycd_writer *w)
{
    for (int i = w->nword; i < w->dpw; i++) w->word[i] = 0;
    uint8_t out[8];
    ycd_store_le64(out, w->base == 10 ? ycd_encode_dec(w->word) : ycd_encode_hex(w->word));
    fwrite(out, 1, 8, w->f);
    w->bytes += 8;
    w->nword = 0;
}

int ycd_writer_append(ycd_writer *w, const uint8_t *digits, size_t n)
{
    if (w->written + n > w->limit)
        return set_err(w->err, sizeof w->err, "block overflow: %" PRIu64 " digits appended to a block of %" PRIu64, w->written + n, w->limit);
    for (size_t i = 0; i < n; i++) {
        unsigned v = digits[i];
        if (w->flags & YCD_ASCII) {
            if (v >= '0' && v <= '9') v -= '0';
            else if (v >= 'a' && v <= 'f') v = v - 'a' + 10;
            else if (v >= 'A' && v <= 'F') v = v - 'A' + 10;
            else return set_err(w->err, sizeof w->err, "byte 0x%02x is not a digit", v);
        }
        if (v >= (unsigned)w->base) return set_err(w->err, sizeof w->err, "digit value %u is out of range for base %d", v, w->base);
        w->word[w->nword++] = (uint8_t)v;
        if (w->nword == w->dpw) flush_word(w);
    }
    w->written += n;
    if (ferror(w->f)) return set_err(w->err, sizeof w->err, "write error: %s", strerror(errno));
    return 0;
}

int ycd_writer_close(ycd_writer *w)
{
    if (!w) return -1;
    if (w->nword) flush_word(w);
    while (w->bytes % YCD_FILE_ALIGN) { fputc(0, w->f); w->bytes++; }
    int bad = ferror(w->f);
    if (fclose(w->f)) bad = 1;
    /* a short block whose set does not record TotalDigits cannot be read back
     * correctly: no reader can tell its digits from the padding */
    if (w->written < w->limit && !w->total_known) bad = 1;
    free(w);
    return bad ? -1 : 0;
}
