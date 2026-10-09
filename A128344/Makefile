PREFIX ?= /opt/homebrew
CC ?= cc
CFLAGS ?= -O3 -Wall -Wextra
CPPFLAGS += -I$(PREFIX)/include
LDFLAGS += -L$(PREFIX)/lib
LDLIBS = -lprimesieve -lgmp -lpthread

a128344: a128344.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f a128344

.PHONY: clean
