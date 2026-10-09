# Makefile for dschain and dspattern -- detects OpenMP and primesieve
# automatically (dschain only; dspattern needs nothing but pthreads).
#
#   make            build ./dschain and ./dspattern
#   make dspattern  build only the pattern search (no dependencies; on a
#                   Mac it includes the Metal GPU search, --gpu)
#   make GPU=0      on a Mac, build dspattern without the GPU search
#   make config     show what was detected, build nothing
#   make test       quick correctness check of both programs
#   make clean

CC      ?= cc
CFLAGS  ?= -O3 -march=native -Wall -Wextra
LDLIBS   = -lm
TARGET   = dschain

UNAME_S := $(shell uname -s)

# ---- primesieve ------------------------------------------------------
# Prefer pkg-config; fall back to probing for the header via brew or the
# usual system paths. If nothing turns up, use the built-in sieve.
HAVE_PS := $(shell pkg-config --exists primesieve 2>/dev/null && echo yes)

ifeq ($(HAVE_PS),yes)
  PS_CFLAGS := $(shell pkg-config --cflags primesieve)
  PS_LIBS   := $(shell pkg-config --libs primesieve)
else
  BREW_PS := $(shell brew --prefix primesieve 2>/dev/null)
  ifneq ($(wildcard $(BREW_PS)/include/primesieve.h),)
    HAVE_PS   := yes
    PS_CFLAGS := -I$(BREW_PS)/include
    PS_LIBS   := -L$(BREW_PS)/lib -lprimesieve
  else ifneq ($(wildcard /usr/include/primesieve.h)$(wildcard /usr/local/include/primesieve.h),)
    HAVE_PS   := yes
    PS_CFLAGS :=
    PS_LIBS   := -lprimesieve
  else
    HAVE_PS   := no
    PS_CFLAGS := -DNO_PRIMESIEVE
    PS_LIBS   :=
  endif
endif

# ---- OpenMP ----------------------------------------------------------
# Apple's stock clang rejects a bare -fopenmp and needs Homebrew's libomp.
ifeq ($(UNAME_S),Darwin)
  BREW_OMP := $(shell brew --prefix libomp 2>/dev/null)
  ifneq ($(wildcard $(BREW_OMP)/include/omp.h),)
    HAVE_OMP  := yes
    OMP_CFLAGS := -Xpreprocessor -fopenmp -I$(BREW_OMP)/include
    OMP_LIBS   := -L$(BREW_OMP)/lib -lomp
  else
    # maybe CC is a real gcc (brew install gcc) that supports -fopenmp directly
    HAVE_OMP := $(shell printf 'int main(void){return 0;}' > /tmp/_omp$$$$.c; \
                  $(CC) -fopenmp /tmp/_omp$$$$.c -o /tmp/_omp$$$$.out 2>/dev/null \
                  && echo yes || echo no; rm -f /tmp/_omp$$$$.c /tmp/_omp$$$$.out)
    ifeq ($(HAVE_OMP),yes)
      OMP_CFLAGS := -fopenmp
      OMP_LIBS   := -fopenmp
    else
      OMP_CFLAGS :=
      OMP_LIBS   :=
    endif
  endif
else
  HAVE_OMP := $(shell printf 'int main(void){return 0;}' > /tmp/_omp$$$$.c; \
                $(CC) -fopenmp /tmp/_omp$$$$.c -o /tmp/_omp$$$$.out 2>/dev/null \
                && echo yes || echo no; rm -f /tmp/_omp$$$$.c /tmp/_omp$$$$.out)
  ifeq ($(HAVE_OMP),yes)
    OMP_CFLAGS := -fopenmp
    OMP_LIBS   := -fopenmp
  else
    OMP_CFLAGS :=
    OMP_LIBS   :=
  endif
endif

ALL_CFLAGS := $(CFLAGS) $(OMP_CFLAGS) $(PS_CFLAGS)
ALL_LIBS   := $(OMP_LIBS) $(PS_LIBS) $(LDLIBS)

# ---- rules -----------------------------------------------------------

all: $(TARGET) dspattern

$(TARGET): dschain.c Makefile
	@echo "  OpenMP     : $(HAVE_OMP)"
	@echo "  primesieve : $(HAVE_PS)"
	$(CC) $(ALL_CFLAGS) -o $@ dschain.c $(ALL_LIBS)
ifneq ($(HAVE_OMP),yes)
	@echo ""
	@echo "  !! OpenMP not found -- --threads will have no effect."
ifeq ($(UNAME_S),Darwin)
	@echo "  !! Fix with:  brew install libomp   (then re-run make)"
else
	@echo "  !! Install a gcc/clang with OpenMP support, then re-run make."
endif
	@echo ""
endif
ifneq ($(HAVE_PS),yes)
	@echo "  (using the built-in segmented sieve; 'brew install primesieve'"
	@echo "   or 'apt install libprimesieve-dev' is 2-3x faster)"
endif

# dspattern is plain C with pthreads and __int128: any 64-bit gcc or clang.
# On a Mac it also gets the Metal GPU search: dspattern_gpu.m includes
# dspattern.c, and xxd embeds the kernel so the binary stands alone.
ifeq ($(UNAME_S),Darwin)
GPU ?= 1
endif

ifeq ($(GPU),1)
dspattern: dspattern.c dspattern_gpu.m dspattern.metal Makefile
	xxd -i dspattern.metal > dspattern_metal.h
	$(CC) $(CFLAGS) -fobjc-arc -pthread -o $@ dspattern_gpu.m \
	    -framework Metal -framework Foundation
else
dspattern: dspattern.c Makefile
	$(CC) $(CFLAGS) -pthread -o $@ dspattern.c
endif

config:
	@echo "uname      : $(UNAME_S)"
	@echo "CC         : $(CC)"
	@echo "OpenMP     : $(HAVE_OMP)   $(OMP_CFLAGS) $(OMP_LIBS)"
	@echo "primesieve : $(HAVE_PS)   $(PS_CFLAGS) $(PS_LIBS)"
	@echo "compile    : $(CC) $(ALL_CFLAGS) -o $(TARGET) dschain.c $(ALL_LIBS)"

test: $(TARGET) dspattern
	@echo "dschain: expecting 12 maximal chains of length >= 6 below 1e8"
	@./$(TARGET) -l 6 -s 2 -e 100M -q | grep -E "chains >="
	@echo "dspattern: expecting the same 12"
	@./dspattern -l 6 -s 2 -e 100M -q | grep -E "chains >="
	@echo "dspattern: expecting a(9) = 56676324799"
	@./dspattern -l 9 -s 2 -e 1e11 --first -q | grep "^earliest"
	@echo "dspattern: expecting a(10) = 63782998969989799429"
	@./dspattern -l 10 -s 56676324800 -e 1e20 --first -q | grep "^earliest"

clean:
	rm -f $(TARGET) dspattern dspattern_metal.h

.PHONY: all config test clean
