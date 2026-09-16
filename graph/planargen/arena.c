/* arena.c -- trivial bump allocator used by planarity_arena.c.
 * The planarity tester allocates ~100 small blocks per call and frees them
 * all before returning, so a reset-per-call arena replaces malloc/free.
 * Every block is preceded by a 16-byte header holding its size (needed by
 * arena_realloc). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#define ARENA_SIZE ((size_t)4 << 20)
static unsigned char arena_buf[ARENA_SIZE] __attribute__((aligned(16)));
static size_t arena_pos;

void arena_reset(void) { arena_pos = 0; }

void *
arena_malloc(size_t sz)
{
    size_t need = (sz + 15) & ~(size_t)15;
    unsigned char *p;
    if (arena_pos + need + 16 > ARENA_SIZE)
    {
        fprintf(stderr, ">E arena overflow (%zu bytes requested)\n", sz);
        exit(2);
    }
    p = arena_buf + arena_pos;
    *(size_t *)p = need;
    arena_pos += 16 + need;
    return p + 16;
}

void *
arena_realloc(void *p, size_t sz)
{
    size_t old;
    void *q;
    if (p == NULL) return arena_malloc(sz);
    old = *(size_t *)((unsigned char *)p - 16);
    q = arena_malloc(sz);
    memcpy(q, p, old < sz ? old : sz);
    return q;
}

void arena_free(void *p) { (void)p; }
