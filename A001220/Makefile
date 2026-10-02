# wieferich: multi-base Wieferich prime search. Needs primesieve (brew install primesieve;
# on the Sparks: the system package or a build in ~/local).
BREW_PREFIX ?= $(shell brew --prefix 2>/dev/null || (test -f $(HOME)/local/include/primesieve.h && echo $(HOME)/local) || echo /usr/local)
CC      ?= cc
CFLAGS  ?= -O2 -std=gnu11 -Wall -Wextra -pthread -I$(BREW_PREFIX)/include
LDFLAGS ?= -L$(BREW_PREFIX)/lib
LDLIBS  ?= -lprimesieve -lm

# a static libprimesieve.a without a shared library next to it also needs the C++ runtime
ifeq ($(wildcard $(BREW_PREFIX)/lib/libprimesieve.so* $(BREW_PREFIX)/lib/libprimesieve*.dylib),)
ifneq ($(wildcard $(BREW_PREFIX)/lib/libprimesieve.a),)
LDLIBS  += -lstdc++
endif
endif

wieferich: wieferich.c wieferich_cks.h
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

wieferich_cks.h: verify_wieferich.py
	python3 verify_wieferich.py cksheader > $@

test: wieferich
	./wieferich selftest

bench: wieferich
	./wieferich bench

clean:
	rm -f wieferich

.PHONY: test bench clean
