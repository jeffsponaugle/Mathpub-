CC      ?= cc
CFLAGS  ?= -O3 -march=native -Wall -Wextra -Wno-unused-function
LDFLAGS ?= -lm

all: a287994

a287994: a287994.c pi_common.h
	$(CC) $(CFLAGS) -pthread -o $@ a287994.c $(LDFLAGS)

# quick correctness check against the known OEIS terms (needs ../A197123/pi-1b.txt)
PI1B ?= ../A197123/pi-1b.txt
test: a287994
	./a287994 -p $(PI1B) -l 1-6 -n 2M 2>/dev/null > /tmp/_a287994_test.txt
	grep -q 'A287994(1) = 24 '      /tmp/_a287994_test.txt && echo "a(1) position OK"
	grep -q 'A290977(1) = 3 '       /tmp/_a287994_test.txt && echo "a(1) value OK"
	grep -q 'A287994(2) = 413 '     /tmp/_a287994_test.txt && echo "a(2) OK"
	grep -q 'A290977(2) = 59 '      /tmp/_a287994_test.txt && echo "a(2) value OK"
	grep -q 'A287994(3) = 326 '     /tmp/_a287994_test.txt && echo "a(3) OK"
	grep -q 'A290977(3) = 209 '     /tmp/_a287994_test.txt && echo "a(3) value OK"
	grep -q 'A287994(4) = 8239 '    /tmp/_a287994_test.txt && echo "a(4) OK"
	grep -q 'A290977(4) = 9314 '    /tmp/_a287994_test.txt && echo "a(4) value OK"
	grep -q 'A287994(5) = 107472 '  /tmp/_a287994_test.txt && echo "a(5) OK"
	grep -q 'A290977(5) = 64015 '   /tmp/_a287994_test.txt && echo "a(5) value OK"
	grep -q 'A287994(6) = 1632152 ' /tmp/_a287994_test.txt && echo "a(6) OK"
	grep -q 'A290977(6) = 886287 '  /tmp/_a287994_test.txt && echo "a(6) value OK"

clean:
	rm -f a287994
.PHONY: all test clean
