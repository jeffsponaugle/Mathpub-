CC      ?= cc
CFLAGS  ?= -O3 -Wall -Wextra
# -mcpu=native on Apple clang / arm64, -march=native elsewhere
ARCH    := $(shell $(CC) -mcpu=native -x c -E /dev/null >/dev/null 2>&1 && echo -mcpu=native || echo -march=native)
LDLIBS  := -lpthread -lm

a253316: a253316.c
	$(CC) $(CFLAGS) $(ARCH) -o $@ $< $(LDLIBS)

test: a253316
	./a253316 --selftest
	./a253316 -n 6 --check-orbits 4

clean:
	rm -f a253316

.PHONY: test clean
