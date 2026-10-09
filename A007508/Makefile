PRIMESIEVE ?= $(shell brew --prefix primesieve 2>/dev/null || echo /usr/local)
CC      ?= cc
CFLAGS  ?= -O3 -march=native -Wall -Wextra -std=c11
CPPFLAGS += -I$(PRIMESIEVE)/include
LDFLAGS  += -L$(PRIMESIEVE)/lib
LDLIBS   += -lprimesieve -lpthread -lm

a007508: a007508.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

clean:
	rm -f a007508

.PHONY: clean
