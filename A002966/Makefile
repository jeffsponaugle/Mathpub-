CC=cc
CFLAGS=-O2
all: a002966_heads a002966_baseline a002966_mitm a002966_tail3 a002966_gstat
a002966_heads: a002966_heads.c
	$(CC) $(CFLAGS) -o $@ $<
a002966_baseline: a002966_baseline.c
	$(CC) $(CFLAGS) -o $@ $<
a002966_mitm: a002966_mitm.c
	$(CC) $(CFLAGS) -o $@ $<
a002966_tail3: a002966_tail3.c
	$(CC) $(CFLAGS) -o $@ $<
a002966_gstat: a002966_gstat.c
	$(CC) $(CFLAGS) -o $@ $<
selftest: a002966_mitm
	./a002966_mitm 5 && ./a002966_mitm 6 && ./a002966_mitm 7
clean:
	rm -f a002966_heads a002966_baseline a002966_mitm a002966_tail3 a002966_gstat
