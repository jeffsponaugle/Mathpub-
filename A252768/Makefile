# A252768 tool. Needs primesieve (brew install primesieve); GMP (brew install gmp) adds independent primality checks.
BREW_PREFIX ?= $(shell brew --prefix 2>/dev/null || echo /usr/local)
CC      ?= cc
CFLAGS  ?= -O2 -std=gnu11 -Wall -Wextra -pthread -I$(BREW_PREFIX)/include
LDFLAGS ?= -L$(BREW_PREFIX)/lib
LDLIBS  ?= -lprimesieve -lm

ifneq ($(wildcard $(BREW_PREFIX)/include/gmp.h),)
CFLAGS  += -DHAVE_GMP
LDLIBS  += -lgmp
endif

a252768: a252768.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

test: a252768
	./a252768 selftest

clean:
	rm -f a252768

.PHONY: test clean
