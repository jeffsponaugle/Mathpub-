# Metal (Apple GPU) version of the swapped 2D A171775 search. The kernel source a171775_2d.metal is
# compiled at run time by the OS Metal compiler (no Xcode Metal toolchain needed); it must sit next
# to the executable or be given with -L. 'make metallib' precompiles it if the toolchain is installed.
CC      ?= clang
CFLAGS  ?= -O3 -Wall -Wextra -Wno-unused-function -fobjc-arc
FRAMEWORKS = -framework Metal -framework Foundation

a171775_metal: a171775_metal.m ../a171775_2d.c
	$(CC) $(CFLAGS) -o $@ a171775_metal.m $(FRAMEWORKS)

metallib: a171775_2d.metal
	xcrun -sdk macosx metal -O3 -c a171775_2d.metal -o a171775_2d.air
	xcrun -sdk macosx metallib a171775_2d.air -o a171775_2d.metallib

clean:
	rm -f a171775_metal a171775_2d.metallib a171775_2d.air
.PHONY: clean metallib
