CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra

hamming: hamming.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean run
run: hamming
	./hamming 1024 100000

clean:
	rm -f hamming
