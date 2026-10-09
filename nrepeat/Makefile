CXX      ?= c++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra
LDFLAGS  ?= -pthread

nrepeat: nrepeat.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f nrepeat

.PHONY: clean
