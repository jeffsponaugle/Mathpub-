# A171775 tool (single file, pthreads, no other dependencies).
CC      ?= cc
CFLAGS  ?= -O3 -std=gnu11 -Wall -Wextra -pthread
LDLIBS  ?= -lm

a171775: a171775.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

test: a171775
	./a171775 selftest

clean:
	rm -f a171775

.PHONY: test clean paranoid

# debug build: recomputes the incremental state at every step
a171775_paranoid: a171775.c
	$(CC) $(CFLAGS) -DPARANOID -o $@ $< $(LDLIBS)

paranoid: a171775_paranoid
	./a171775_paranoid paranoid

# swapped 2D method (enumerate (n-1)-digit palindromes, lookup two digits)
a171775_2d: a171775_2d.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

a171775_2d_paranoid: a171775_2d.c
	$(CC) $(CFLAGS) -DPARANOID -o $@ $< $(LDLIBS)
