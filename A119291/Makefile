BREW_PREFIX ?= $(shell brew --prefix 2>/dev/null || echo /usr/local)
CC      ?= cc
CFLAGS  += -O3 -march=native -Wall -Wextra -I$(BREW_PREFIX)/include
LDFLAGS += -L$(BREW_PREFIX)/lib -lprimesieve -lpthread

pdigits: pdigits.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

.PHONY: clean test
clean:
	rm -f pdigits

test: pdigits
	./pdigits --selftest
