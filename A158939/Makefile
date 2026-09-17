# A158939 tool. Needs primesieve (brew install primesieve); GMP (brew install gmp) adds an independent check to 'verify'.
BREW_PREFIX ?= $(shell brew --prefix 2>/dev/null || echo /usr/local)
CC      ?= cc
CFLAGS  ?= -O2 -std=gnu11 -Wall -Wextra -pthread -I$(BREW_PREFIX)/include
LDFLAGS ?= -L$(BREW_PREFIX)/lib
LDLIBS  ?= -lprimesieve -lm

ifneq ($(wildcard $(BREW_PREFIX)/include/gmp.h),)
CFLAGS  += -DHAVE_GMP
LDLIBS  += -lgmp
endif

a158939: a158939.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

test: a158939
	./a158939 selftest

bench: a158939
	./a158939 bench

clean:
	rm -f a158939

.PHONY: test bench clean
