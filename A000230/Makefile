CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra
PS_CFLAGS := $(shell pkg-config --cflags primesieve)
PS_LIBS   := $(shell pkg-config --libs primesieve)

a000230: a000230.c
	$(CC) $(CFLAGS) $(PS_CFLAGS) -o $@ $< $(PS_LIBS) -lpthread

clean:
	rm -f a000230

.PHONY: clean
