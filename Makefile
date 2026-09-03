CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -pthread -Wall -Wextra

all: hamming verify

hamming: hamming.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

verify: verify.cpp hamming.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: all test run clean
test: verify
	./verify

run: hamming
	./hamming 1024 1M

clean:
	rm -f hamming verify
