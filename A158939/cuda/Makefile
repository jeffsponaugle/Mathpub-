# CUDA version of the A158939 scan (DGX Spark / GB10: sm_121, CUDA 13). No library dependencies.
NVCC    ?= /usr/local/cuda/bin/nvcc
ARCH    ?= sm_121
NVFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -Xcompiler -Wall

a158939_cuda: a158939_cuda.cu
	$(NVCC) $(NVFLAGS) -o $@ $<

test: a158939_cuda
	./a158939_cuda selftest

clean:
	rm -f a158939_cuda
.PHONY: test clean
