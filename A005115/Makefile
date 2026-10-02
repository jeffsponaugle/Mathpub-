# A005115 tool (CPU reference). No library dependencies.
CC      ?= cc
CFLAGS  ?= -O3 -std=gnu11 -Wall -Wextra -pthread -march=native
LDLIBS  ?= -lm

a005115: a005115.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

test: a005115
	./a005115 selftest

bench: a005115
	./a005115 bench

clean:
	rm -f a005115

.PHONY: test bench clean
