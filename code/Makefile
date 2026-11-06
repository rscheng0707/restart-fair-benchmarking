# Makefile — time-fair benchmarking harness
#
# Toolchain is an experimental factor (paper §3.5). Build with each compiler
# under test and record the exact version in the output metadata:
#   make CXX=g++-13
#   make CXX=g++-14
#   make CXX=clang++
#
# -O2 vs -O3 is itself a factor worth reporting. Do NOT use -march=native for
# the main matrix unless it is recorded as a factor level: it changes generated
# code based on the build machine, which silently couples results to hardware
# and would confound the toolchain comparison.

CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -fopenmp -Wall -Wextra
TARGETS   = test_equiv pilot

all: $(TARGETS)

test_equiv: test_equivalence.cpp benchmarks.hpp pso.hpp
	$(CXX) $(CXXFLAGS) test_equivalence.cpp -o $@

pilot: pilot.cpp benchmarks.hpp pso.hpp harness.hpp
	$(CXX) $(CXXFLAGS) pilot.cpp -o $@

# The equivalence gate must pass before any pilot run.
check: test_equiv
	./test_equiv

version:
	@$(CXX) --version | head -1
	@echo "CXXFLAGS = $(CXXFLAGS)"

clean:
	rm -f $(TARGETS)

.PHONY: all check version clean
