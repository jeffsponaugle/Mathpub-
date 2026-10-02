# CUDA version of the A350826 tool (DGX Spark / GB10: sm_121). No dependencies beyond CUDA.
NVCC    ?= /usr/local/cuda/bin/nvcc
ARCH    ?= sm_121
NVFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -Xcompiler -Wall -Xptxas -v

a350826_cuda: a350826_cuda.cu
	$(NVCC) $(NVFLAGS) -o $@ $<

clean:
	rm -f a350826_cuda
.PHONY: clean
