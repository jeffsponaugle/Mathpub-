# A004029 / A006227 / A395859: space-group counts via CARAT
CARAT   = carat
CC      = clang
CFLAGS  = -O2 -std=gnu11 -w -I$(CARAT)/include -I/opt/homebrew/include
LDLIBS  = -L/opt/homebrew/lib -lgmp -lm

all: enantio zverify extend1

carat:
	$(MAKE) -C $(CARAT) -f Makefile.local -j14 all

enantio: src/enantio.c $(CARAT)/libfunctions.a
	$(CC) $(CFLAGS) $< $(CARAT)/libfunctions.a $(LDLIBS) -o $@

clean:
	rm -f enantio

.PHONY: all carat clean

zverify: src/zverify.cpp src/carat_io.hpp
	clang++ -O2 -std=c++17 -Wall -I/opt/homebrew/include $< -L/opt/homebrew/lib -lgmpxx -lgmp -o $@

extend1: src/extend1.cpp src/carat_io.hpp
	clang++ -O2 -std=c++17 -Wall -I/opt/homebrew/include $< -L/opt/homebrew/lib -lgmpxx -lgmp -o $@
