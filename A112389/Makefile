# A112389 tool (LEGO 2x4 buildings).  No external libraries.
CC      ?= cc
CFLAGS  ?= -O3 -march=native -std=gnu11 -Wall -Wextra -pthread
LDLIBS  ?= -lpthread

lego: lego.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

test: lego
	./lego -s

clean:
	rm -f lego

.PHONY: test clean
