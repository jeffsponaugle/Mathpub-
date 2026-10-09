CC ?= cc
CFLAGS ?= -O3 -march=native -Wall -Wextra

pi_missing: pi_missing.c
	$(CC) $(CFLAGS) -pthread -o $@ $<

clean:
	rm -f pi_missing

.PHONY: clean
