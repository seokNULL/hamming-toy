CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -pthread -Wall -Wextra

hamming: hamming.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean run
run: hamming
	./hamming 1024 1M

clean:
	rm -f hamming
