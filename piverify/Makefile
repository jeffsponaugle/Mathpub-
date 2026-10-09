CC ?= cc
CFLAGS ?= -O3 -Wall -Wextra
CFLAGS += -D_FILE_OFFSET_BITS=64

# Locate GMP: pkg-config first (Linux), then Homebrew (macOS), else default paths
GMP_LIBS := $(shell pkg-config --libs gmp 2>/dev/null)
ifneq ($(GMP_LIBS),)
CFLAGS += $(shell pkg-config --cflags gmp 2>/dev/null)
else
GMP_PREFIX := $(shell brew --prefix gmp 2>/dev/null)
ifneq ($(GMP_PREFIX),)
CFLAGS += -I$(GMP_PREFIX)/include
LDFLAGS += -L$(GMP_PREFIX)/lib
endif
GMP_LIBS := -lgmp
endif

all: piverify bbp

piverify: piverify.c bbp.c bbp.h
	@echo '#include <gmp.h>' | $(CC) $(CFLAGS) -E -xc - >/dev/null 2>&1 || { \
	  echo ""; \
	  echo "error: gmp.h not found. piverify needs the GMP development package:"; \
	  echo "  Debian/Ubuntu:  sudo apt install libgmp-dev"; \
	  echo "  Fedora/RHEL:    sudo dnf install gmp-devel"; \
	  echo "  Arch:           sudo pacman -S gmp"; \
	  echo "  macOS:          brew install gmp"; \
	  echo ""; \
	  echo "(\"make bbp\" builds the standalone digit extractor, which"; \
	  echo " needs no GMP)"; \
	  echo ""; \
	  exit 1; }
	$(CC) $(CFLAGS) -o piverify piverify.c bbp.c $(LDFLAGS) $(GMP_LIBS) -lpthread

bbp: bbp_main.c bbp.c bbp.h
	$(CC) $(CFLAGS) -o bbp bbp_main.c bbp.c -lpthread

clean:
	rm -f piverify bbp

.PHONY: all clean
