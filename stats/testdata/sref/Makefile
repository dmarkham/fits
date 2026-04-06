# Siril/RawTherapee reference harness for validating the Go stats package
# against the real findMinMaxPercentile implementation.
#
#   make          -> build sref
#   make goldens  -> generate golden fixtures
#   make clean    -> remove binary + goldens

CXX     ?= c++
CXXFLAGS ?= -O2 -Wall -std=c++17
LDLIBS  := -lm

sref: sref.cc rt_algo_excerpt.cc
	$(CXX) $(CXXFLAGS) -o sref sref.cc $(LDLIBS)

.PHONY: goldens
goldens: sref
	mkdir -p golden
	./sref gen golden

.PHONY: clean
clean:
	rm -f sref
	rm -rf golden
