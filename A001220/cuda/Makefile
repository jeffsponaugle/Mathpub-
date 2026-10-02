# CUDA version of the multi-base Wieferich search (DGX Spark / GB10: sm_121, CUDA 13).
# Needs primesieve on the host side: system package or a build in ~/local.
NVCC      ?= /usr/local/cuda/bin/nvcc
ARCH      ?= sm_121
PS_PREFIX ?= $(shell test -f $(HOME)/local/include/primesieve.h && echo $(HOME)/local || echo /usr)
NVFLAGS   ?= -O3 -std=c++17 -arch=$(ARCH) -Xcompiler -Wall,-pthread -I$(PS_PREFIX)/include
LDLIBS    ?= -L$(PS_PREFIX)/lib -lprimesieve -lpthread

wieferich_cuda: wieferich_cuda.cu ../wieferich_cks.h
	$(NVCC) $(NVFLAGS) -o $@ $< $(LDLIBS)

regs: wieferich_cuda.cu
	$(NVCC) $(NVFLAGS) -Xptxas -v -c -o /dev/null $< 2>&1 | grep -A1 k_wieferich

test: wieferich_cuda
	./wieferich_cuda selftest

clean:
	rm -f wieferich_cuda
.PHONY: regs test clean
