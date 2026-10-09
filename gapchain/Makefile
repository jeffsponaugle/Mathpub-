# Build gapsieve and gapchain_ps -- see README.md.
#
#   make              build both programs
#   make gapsieve     build only the sieve (no dependencies; on a Mac it
#                     includes the Metal GPU search, --gpu)
#   make gapchain_ps  build only the enumeration tool (needs primesieve)
#   make clean        remove the binaries
#   make NATIVE=      build for a generic CPU instead of this one
#   make GPU=0        on a Mac, build gapsieve without the GPU search
#   make CUDA=1       on Linux with an NVIDIA GPU and the CUDA toolkit, build
#                     gapsieve with the CUDA GPU search (--gpu)
#
# Override on the command line, e.g.  make CC=gcc CFLAGS=-O2

CC         ?= cc
PKG_CONFIG ?= pkg-config

# Compile for the CPU doing the build.  On an x86 box this doubled the sieve's
# speed (55 -> 112 T/s), since compilers otherwise target a generic x86-64
# without AVX; on Apple Silicon it makes no difference.  The binary may not run
# on older CPUs, so build on the machine that runs the search.  Compilers name
# the flag differently, so use the first one this compiler accepts.
ifeq ($(origin NATIVE),undefined)
NATIVE := $(shell for f in -march=native -mcpu=native; do \
            $(CC) $$f -x c -c /dev/null -o /dev/null 2>/dev/null && \
            echo $$f && break; done)
endif

CFLAGS     ?= -O3 -Wall -Wextra $(NATIVE)

# primesieve, for gapchain_ps only: pkg-config if it knows the library,
# otherwise Homebrew's prefix, otherwise the compiler's default paths.
PS_CFLAGS := $(shell $(PKG_CONFIG) --cflags primesieve 2>/dev/null)
PS_LIBS   := $(shell $(PKG_CONFIG) --libs primesieve 2>/dev/null)
ifeq ($(strip $(PS_LIBS)),)
# (brew is often not on the PATH of a non-interactive ssh session)
BREW_PREFIX := $(shell brew --prefix 2>/dev/null || \
                 { test -d /opt/homebrew && echo /opt/homebrew; })
ifneq ($(BREW_PREFIX),)
PS_CFLAGS := -I$(BREW_PREFIX)/include
PS_LIBS   := -L$(BREW_PREFIX)/lib -lprimesieve
else
PS_LIBS   := -lprimesieve
endif
endif

all: gapsieve gapchain_ps

# On a Mac, gapsieve also gets the Metal GPU search: gapsieve_gpu.m includes
# gapsieve.c, and xxd embeds the kernels so the binary stands alone.
ifeq ($(shell uname -s),Darwin)
GPU ?= 1
endif

# With CUDA=1, gapsieve gets the CUDA GPU search instead: gapsieve_cuda_host.c
# includes gapsieve.c, and gapsieve_cuda.cu holds the kernels.  -arch=native
# compiles for the GPU in the machine doing the build (CUDA 11.5 and later),
# such as sm_121 on a DGX Spark or sm_87 on a Jetson Orin.  To build for a
# different GPU, set CUDA_ARCH, e.g. make CUDA=1 CUDA_ARCH=-arch=sm_87.
NVCC      ?= nvcc
CUDA_ARCH ?= -arch=native

ifeq ($(CUDA),1)
gapsieve: gapsieve.c gapsieve_cuda_host.c gapsieve_cuda.cu gapsieve_cuda.h
	$(CC) $(CFLAGS) -pthread -c -o gapsieve_cuda_host.o gapsieve_cuda_host.c
	$(NVCC) -O3 $(CUDA_ARCH) -o $@ gapsieve_cuda.cu gapsieve_cuda_host.o -lpthread
else ifeq ($(GPU),1)
gapsieve: gapsieve.c gapsieve_gpu.m gapsieve.metal
	xxd -i gapsieve.metal > gapsieve_metal.h
	$(CC) $(CFLAGS) -fobjc-arc -pthread -o $@ gapsieve_gpu.m \
	    -framework Metal -framework Foundation
else
gapsieve: gapsieve.c
	$(CC) $(CFLAGS) -pthread -o $@ gapsieve.c
endif

# -lm for sqrtl, which glibc keeps in libm
gapchain_ps: gapchain_ps.c
	$(CC) $(CFLAGS) -pthread $(PS_CFLAGS) -o $@ gapchain_ps.c $(PS_LIBS) -lm

clean:
	rm -f gapsieve gapchain_ps gapsieve_metal.h gapsieve_cuda_host.o

.PHONY: all clean
