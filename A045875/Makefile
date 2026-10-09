CC=cc
CFLAGS=-O3 -Wall -Wextra -I/opt/homebrew/include
LDFLAGS=-L/opt/homebrew/lib -lgmp -lpthread

a045875: a045875.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f a045875
