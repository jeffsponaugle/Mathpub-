PREFIX ?= /opt/homebrew
CC ?= cc
CFLAGS ?= -O3 -Wall -Wextra
CPPFLAGS += -I$(PREFIX)/include
LDFLAGS += -L$(PREFIX)/lib
LDLIBS += -lprimesieve

a050248: a050248.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

clean:
	rm -f a050248

.PHONY: clean
