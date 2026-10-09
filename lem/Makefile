CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -std=c11 -D_FILE_OFFSET_BITS=64
# Homebrew paths (harmless elsewhere)
CFLAGS  += -I/opt/homebrew/include -I/usr/local/include
LDFLAGS += -L/opt/homebrew/lib -L/usr/local/lib
LIBS     = -lmpfr -lgmp -lm

all: lemverify

lemverify: lemverify.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

clean:
	rm -f lemverify

.PHONY: all clean
