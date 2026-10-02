# psph -- exhaustive search for extremal postage stamp bases
# macOS/Apple Silicon: make            (clang, NEON)
# Linux x86-64:        make linux      (gcc or clang, AVX2 via -march=native; add NOSIMD=1 for scalar)
# Linux aarch64:       make linux      (NEON)
CC      ?= clang
ARCHFLAGS ?= -mcpu=apple-m1
CFLAGS  = -O3 $(ARCHFLAGS) -std=c11 -pthread -Wall -Wextra -Wno-unused-parameter $(if $(NOSIMD),-DNO_SIMD,)
LDFLAGS = -lm

all: psph psph16

psph: psph.c
	$(CC) $(CFLAGS) -o psph psph.c $(LDFLAGS)

# 16-bit stamp counts: needed for h >= 255 (e.g. k=4, h=302/303); ~1.6x slower per candidate
psph16: psph.c
	$(CC) $(CFLAGS) -DWIDE -o psph16 psph.c $(LDFLAGS)

linux:
	$(MAKE) CC=$(or $(LCC),gcc) ARCHFLAGS=-march=native all

selftest: psph
	./psph -selftest quick

clean:
	rm -f psph psph16

.PHONY: all linux selftest clean
