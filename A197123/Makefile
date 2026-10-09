CC      ?= cc
CFLAGS  ?= -O3 -march=native -Wall -Wextra -Wno-unused-function -pthread
LDFLAGS ?= -lm

all: pibloom pibloomt pisearch pifind piself pimatchscan pimatchsort pimatchverify ycdtool ycd_test

pibloom: pibloom.c pi_common.h
	$(CC) $(CFLAGS) -o $@ pibloom.c $(LDFLAGS)

pisearch: pisearch.c pi_common.h
	$(CC) $(CFLAGS) -o $@ pisearch.c $(LDFLAGS)

pibloomt: pibloomt.c pi_common.h
	$(CC) $(CFLAGS) -pthread -o $@ pibloomt.c $(LDFLAGS)

pifind: pifind.c pi_common.h
	$(CC) $(CFLAGS) -o $@ pifind.c $(LDFLAGS)

pimatchscan: pimatchscan.c pimatch.h pi_common.h
	$(CC) $(CFLAGS) -o $@ pimatchscan.c $(LDFLAGS)

pimatchsort: pimatchsort.c pimatch.h pi_common.h
	$(CC) $(CFLAGS) -o $@ pimatchsort.c $(LDFLAGS)

pimatchverify: pimatchverify.c pimatch.h pi_common.h
	$(CC) $(CFLAGS) -o $@ pimatchverify.c $(LDFLAGS)

# standalone .ycd library (ycd.h + ycd.c, no dependency on the pi tools) and its front end / self-test
ycdtool: ycdtool.c ycd.c ycd.h
	$(CC) -O2 -std=c99 -Wall -Wextra -o $@ ycdtool.c ycd.c

ycd_test: ycd_test.c ycd.c ycd.h
	$(CC) -O2 -std=c99 -Wall -Wextra -o $@ ycd_test.c ycd.c

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
	rm -rf /tmp/_pimatch && mkdir -p /tmp/_pimatch
	./pimatchscan -p pi-1b.txt -n 100M -d /tmp/_pimatch -L 10 -k 9 -P 4 -t 2 -C 0 > /dev/null
	./pimatchsort -d /tmp/_pimatch -t 2 > /dev/null
	./pimatchverify -d /tmp/_pimatch -E 40 | grep -q 'n=10  4392366484  at 182,105 and 240,479' && echo "pimatch OK"
	rm -rf /tmp/_pimatch /tmp/_pimatch2 && mkdir -p /tmp/_pimatch /tmp/_pimatch2
	./pimatchscan -p pi-1b.txt -n 100M -d /tmp/_pimatch,/tmp/_pimatch2 -L 10 -k 9 -P 4 -t 2 -C 0 -H 1/2 > /dev/null
	./pimatchsort -d /tmp/_pimatch -t 2 -D > /dev/null
	./pimatchverify -d /tmp/_pimatch -E 40 | grep -q 'n=10  4392366484  at 182,105 and 240,479' && echo "pimatch two-dir hash-pass OK"
	rm -rf /tmp/_pimatch /tmp/_pimatch2
	./ycd_test /tmp | grep -q 'all passed' && echo "ycd library OK"
	rm -rf /tmp/_ycd && mkdir -p /tmp/_ycd
	./ycdtool encode pi-1b.txt /tmp/_ycd -b 3000007 -n 10M -N "Pi - Dec - Test" > /dev/null
	./pifind -p /tmp/_ycd 4392366484 | grep -q 'position 182,105' && ./ycdtool dump /tmp/_ycd 182105 10 | grep -q '4392366484' && echo "ycdtool <-> pifind OK"
	rm -rf /tmp/_ycd

clean:
	rm -f pibloom pibloomt pisearch pifind piself pimatchscan pimatchsort pimatchverify ycdtool ycd_test
.PHONY: all test clean
