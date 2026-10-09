CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra

all: hexcount

hexcount: hexcount.c
	$(CC) $(CFLAGS) -o $@ $<

# Generate 1e5 digits of Pi and verify all 16 sequences against OEIS
test: hexcount
	python3 gen_pihex.py 100000 pihex_1e5.txt
	./hexcount --check-pi pihex_1e5.txt

clean:
	rm -f hexcount pihex_1e5.txt counts.csv

.PHONY: all test clean
