# dbpal - double-base palindrome search (CPU + GPU: Apple Metal on macOS, NVIDIA CUDA on Linux)
CXX      ?= clang++
BREW     ?= /opt/homebrew
CUDA     ?= /usr/local/cuda
CXXFLAGS ?= -std=c++17 -O3 -mcpu=native -Wall -Wextra -Wno-unused-parameter
NLIMB    ?= 3
CPPFLAGS += -isystem $(BREW)/include -Isrc -DNLIMB=$(NLIMB)
LDLIBS   += -L$(BREW)/lib -lgmpxx -lgmp -lpthread

SRC  := src/main.cpp src/plan.cpp src/plan2.cpp src/table.cpp src/engine_cpu.cpp src/engine2_cpu.cpp
HDR  := $(wildcard src/*.h)

UNAME := $(shell uname -s)

ifeq ($(UNAME),Darwin)
  GPUSRC   := src/gpu_metal.mm
  GPUINC   := src/kernels_src.inc
  GPUFLAGS := -fobjc-arc
  LDLIBS   += -framework Metal -framework Foundation
else ifneq ($(wildcard $(CUDA)/include/nvrtc.h),)
  GPUSRC   := src/gpu_cuda.cpp
  GPUINC   := src/kernels_cuda_src.inc
  CPPFLAGS += -I$(CUDA)/include
  LDLIBS   += -L$(CUDA)/lib64 -Wl,-rpath,$(CUDA)/lib64 -lcuda -lnvrtc
else
  GPUSRC   := src/gpu_stub.cpp
endif
TARGETS := dbpal

all: $(TARGETS)

dbpal: $(SRC) $(GPUSRC) $(HDR) $(GPUINC)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(SRC) $(GPUSRC) -o $@.tmp $(LDLIBS) $(GPUFLAGS)
	mv -f $@.tmp $@   # atomic replace: running searches keep the old binary

# 4-limb build for numbers beyond the 3-limb range (~10^54); campaign.py switches to it by itself.
dbpal4: $(SRC) $(GPUSRC) $(HDR) $(GPUINC)
	$(CXX) $(CXXFLAGS) $(filter-out -DNLIMB=%,$(CPPFLAGS)) -DNLIMB=4 $(SRC) $(GPUSRC) -o $@.tmp $(LDLIBS) $(GPUFLAGS)
	mv -f $@.tmp $@

# The GPU kernels are compiled at run time (Metal: newLibraryWithSource, CUDA: NVRTC) with the
# constants of each search, so no offline GPU toolchain is needed: core.h + the kernel file are
# embedded as one raw string literal.
src/kernels_src.inc: src/core.h src/kernels.metal
	( echo 'R"METALSRC('; grep -v '#pragma once' src/core.h; grep -v '#include "core.h"' src/kernels.metal; echo ')METALSRC"' ) > $@

src/kernels_cuda_src.inc: src/core.h src/kernels.cu
	( echo 'R"CUDASRC('; grep -v '#pragma once' src/core.h; grep -v '#include "core.h"' src/kernels.cu; echo ')CUDASRC"' ) > $@

# offline-compiled library (optional; needs `xcodebuild -downloadComponent MetalToolchain`)
dbpal.metallib: src/kernels.metal src/core.h
	xcrun -sdk macosx metal -std=metal3.1 -O3 -Isrc -c src/kernels.metal -o kernels.air
	xcrun -sdk macosx metallib kernels.air -o $@
	rm -f kernels.air

cpu: dbpal-cpu
dbpal-cpu: $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(SRC) src/gpu_stub.cpp -o dbpal-cpu.tmp -L$(BREW)/lib -lgmpxx -lgmp -lpthread
	mv -f dbpal-cpu.tmp dbpal-cpu

test: all cpu dbpal4
	bash tests/run_all.sh

clean:
	rm -f dbpal dbpal4 dbpal-cpu dbpal.metallib kernels.air src/kernels_src.inc src/kernels_cuda_src.inc

.PHONY: all cpu test clean
