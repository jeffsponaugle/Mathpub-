# A350826 tool (prime sextuplet counts). No dependencies.
CC      ?= cc
CFLAGS  ?= -O3 -std=gnu11 -Wall -Wextra -pthread
LDLIBS  ?= -lm

a350826: a350826.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

test: a350826
	./a350826 selftest

clean:
	rm -f a350826

.PHONY: test clean
