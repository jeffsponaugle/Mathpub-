PREFIX ?= /opt/homebrew
CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -march=native -Wall -Wextra
CPPFLAGS += -I$(PREFIX)/include
LDFLAGS += -L$(PREFIX)/lib
LDLIBS += -lprimesieve -lpthread

digitprimes: digitprimes.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f digitprimes

.PHONY: clean
