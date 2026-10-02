# CUDA version of the A171775 search (DGX Spark / GB10: sm_121). No dependencies beyond CUDA.
NVCC    ?= /usr/local/cuda/bin/nvcc
ARCH    ?= sm_121
NVFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -Xcompiler -pthread -Xcompiler -Wall -Xptxas -v

a171775_cuda: a171775_cuda.cu
	$(NVCC) $(NVFLAGS) -o $@ $< 2> build.log || (cat build.log; false)
	@grep -E "registers|spill" build.log | sort | uniq -c | head -20

clean:
	rm -f a171775_cuda build.log
.PHONY: clean

# swapped 2D method (host logic from ../a171775_2d.c)
a171775_2d_cuda: a171775_2d_cuda.cu ../a171775_2d.c
	$(NVCC) -O3 -std=c++17 -arch=$(ARCH) -Xcompiler -pthread -Xcompiler -Wall -Xcompiler -Wno-unused-function -Xptxas -v -o $@ a171775_2d_cuda.cu 2> build2d.log || (cat build2d.log; false)
	@grep -E "registers|spill" build2d.log | sort | uniq -c | head
