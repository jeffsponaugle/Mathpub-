GMP_PREFIX ?= $(shell brew --prefix gmp 2>/dev/null || echo /usr/local)
CC      ?= cc
CFLAGS  += -O3 -march=native -Wall -Wextra -I$(GMP_PREFIX)/include
LDFLAGS += -L$(GMP_PREFIX)/lib -lgmp -lpthread

cchain: cchain.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

.PHONY: clean test
clean:
	rm -f cchain

test: cchain
	./cchain --verify-known
	./cchain --selftest
