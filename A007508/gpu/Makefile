# GPU version for the DGX Spark (GB10, sm_121, CUDA 13). Needs primesieve in $(PS).
NVCC   ?= /usr/local/cuda/bin/nvcc
ARCH   ?= sm_121
PS     ?= $(HOME)/local
NVFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -Xcompiler -Wall
DEFS   ?=

a007508_cuda: a007508_cuda.cu
	$(NVCC) $(NVFLAGS) $(DEFS) -I$(PS)/include -o $@ $< -L$(PS)/lib -lprimesieve

clean:
	rm -f a007508_cuda
.PHONY: clean
