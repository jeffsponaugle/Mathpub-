# A335406 tool. Needs primesieve (brew install primesieve).
BREW_PREFIX ?= $(shell brew --prefix 2>/dev/null || echo /usr/local)
CC      ?= cc
CFLAGS  ?= -O2 -std=gnu11 -Wall -Wextra -pthread -I$(BREW_PREFIX)/include
LDFLAGS ?= -L$(BREW_PREFIX)/lib
LDLIBS  ?= -lprimesieve -lm

a335406: a335406.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

test: a335406
	./a335406 selftest

clean:
	rm -f a335406

.PHONY: test clean
