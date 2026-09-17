# A137723 tool. Needs primesieve (brew install primesieve).
BREW_PREFIX ?= $(shell brew --prefix 2>/dev/null || echo /usr/local)
CC      ?= cc
CFLAGS  ?= -O2 -std=gnu11 -Wall -Wextra -pthread -I$(BREW_PREFIX)/include
LDFLAGS ?= -L$(BREW_PREFIX)/lib
LDLIBS  ?= -lprimesieve

# GMP (brew install gmp) enables 'search -l' above 2^127
ifneq ($(wildcard $(BREW_PREFIX)/include/gmp.h),)
CFLAGS  += -DHAVE_GMP
LDLIBS  += -lgmp
endif

a137723: a137723.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

test: a137723
	./a137723 selftest

clean:
	rm -f a137723

.PHONY: test clean
