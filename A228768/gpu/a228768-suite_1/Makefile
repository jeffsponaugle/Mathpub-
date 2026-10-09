# OEIS A228768 search -- CPU and CUDA builds
#
#   make            CPU search tool          -> a228768x
#   make gpu        CUDA build (DGX Spark)   -> a228768gpu
#   make emu        GPU code on the host     -> a228768emu   (no GPU needed)
#   make test       run every correctness check
#
# The GPU arch is auto-detected: sm_121 (GB10) when the toolkit knows it,
# otherwise PTX for compute_120 which the driver JITs at load.

CC     ?= cc
CXX    ?= c++
NVCC   ?= nvcc
CFLAGS ?= -O3 -pthread
ARCH   ?= $(shell $(NVCC) --list-gpu-code 2>/dev/null | grep -q sm_121 \
                  && echo -arch=sm_121 \
                  || echo -gencode arch=compute_120,code=compute_120)

all: a228768x
gpu: a228768gpu
emu: a228768emu

a228768x: a228768x.c
	$(CC) $(CFLAGS) -o $@ $< -lm

a228768lib.o: a228768x.c
	$(CC) -O3 -c -DA228768_LIB $< -o $@

a228768gpu: a228768_gpu.cu a228768lib.o
	$(NVCC) -O3 $(ARCH) a228768_gpu.cu a228768lib.o -o $@ -lpthread -lm

a228768emu: a228768_gpu.cu a228768lib.o
	$(CXX) -O3 -DEMULATE -x c++ a228768_gpu.cu -x none a228768lib.o -o $@ -lpthread -lm

test: a228768x a228768emu
	./a228768x --self-test
	./a228768emu --nbase 15 --start 3e19 --block-log 20 --sieve-limit 65536 --blocks 1 --verify
	./a228768emu --nbase 14 --start 19775643790000000 --limit 19775643800000000 --block-log 21

clean:
	rm -f a228768x a228768gpu a228768emu a228768lib.o

.PHONY: all gpu emu test clean
