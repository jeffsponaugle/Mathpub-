CC      ?= cc
CFLAGS  ?= -O3 -march=native -Wall -Wextra -Wno-unused-function
LDFLAGS ?= -lm

all: pibloom pibloomt pisearch pifind piself

pibloom: pibloom.c pi_common.h
	$(CC) $(CFLAGS) -o $@ pibloom.c $(LDFLAGS)

pisearch: pisearch.c pi_common.h
	$(CC) $(CFLAGS) -o $@ pisearch.c $(LDFLAGS)

pibloomt: pibloomt.c pi_common.h
	$(CC) $(CFLAGS) -pthread -o $@ pibloomt.c $(LDFLAGS)

pifind: pifind.c pi_common.h
	$(CC) $(CFLAGS) -o $@ pifind.c $(LDFLAGS)

piself: piself.c pi_common.h
	$(CC) $(CFLAGS) -o $@ piself.c $(LDFLAGS)

# quick correctness check against known OEIS terms (needs pi-1b.txt)
test: all
	./pibloom  -p pi-1b.txt -l 7 -n 10K -b 1M -o /tmp/_t7.txt > /dev/null
	./pisearch -p pi-1b.txt -l 7 -n 10K -c /tmp/_t7.txt -o /tmp/_r7.txt | grep -q '8530614  at positions 4,167 and 4,601' && echo "a(7) OK"
	./pibloomt -p pi-1b.txt -l 7 -n 10K -b 1M -o /tmp/_t7t.txt > /dev/null
	./pisearch -p pi-1b.txt -l 7 -n 10K -c /tmp/_t7t.txt -o /tmp/_r7t.txt | grep -q '8530614  at positions 4,167 and 4,601' && echo "a(7) OK (threaded)"
	./pifind -p pi-1b.txt -n 10K 8530614 | grep -q 'position 4,167' && echo "pifind OK"
	./piself pi-1b.txt -n 100K | grep -q '16470 begins at position 16,470' && echo "piself OK"

clean:
	rm -f pibloom pibloomt pisearch pifind piself
.PHONY: all test clean
