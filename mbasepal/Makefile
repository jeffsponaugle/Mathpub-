CC=cc
CFLAGS=-O3 -Wall -Wextra

a171740: a171740.c
	$(CC) $(CFLAGS) -o $@ $< -lpthread

clean:
	rm -f a171740

a171740g: a171740g_host.m kernel_src.h
	clang -O3 -fobjc-arc -framework Metal -framework Foundation -o a171740g a171740g_host.m
