# CUDA version of the A359636 scan (DGX Spark / GB10: sm_121). Needs libprimesieve-dev.
NVCC   ?= nvcc
ARCH   ?= sm_121
NVFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) -Xcompiler -pthread -Xcompiler -Wall

a359636_cuda: a359636_cuda.cu
	$(NVCC) $(NVFLAGS) -o $@ $< -lprimesieve

clean:
	rm -f a359636_cuda
.PHONY: clean
