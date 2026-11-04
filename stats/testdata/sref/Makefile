# Siril reference harness for the Go stats package. Every expected value
# comes from Siril's own source, copied verbatim (extract_siril.py).
#
#   make                 -> build sref (gitignored); BIN=path builds elsewhere
#   make goldens         -> regenerate golden/siril_stats_golden.json
#   make check-excerpt   -> prove siril_excerpt.c / rt_algo_excerpt.cc are
#                           verbatim copies of the Siril source at SIRIL_SRC
#   make clean           -> remove the binary
#
# Flags follow Siril's meson release build (-O3 -fno-math-errno, no
# -march=native). -ffp-contract=off keeps float expressions unfused on
# every architecture, which is what an x86-64 baseline build does anyway
# (no FMA instructions available).

SIRIL_SRC ?= $(error set SIRIL_SRC to an unpacked Siril 1.4.0 source tree)
CC       ?= cc
CXX      ?= c++
FPFLAGS  := -O3 -fno-math-errno -ffp-contract=off
CFLAGS   := $(FPFLAGS) -std=gnu11 -Wall -Wno-maybe-uninitialized
CXXFLAGS := $(FPFLAGS) -std=c++17 -Wall
LDLIBS   := -lm

BIN ?= sref

# The Siril C code is compiled as C (not C++); the object is temporary.
$(BIN): sref.cc rt_algo_excerpt.cc siril_wrappers.c siril_excerpt.c siril_stubs.h sha256.h
	$(CC) $(CFLAGS) -c -o $(BIN).siril.o siril_wrappers.c
	$(CXX) $(CXXFLAGS) -o $(BIN) sref.cc $(BIN).siril.o $(LDLIBS)
	rm -f $(BIN).siril.o

.PHONY: goldens
goldens: $(BIN)
	mkdir -p golden
	./$(BIN) gen golden

.PHONY: check-excerpt
check-excerpt:
	python3 extract_siril.py --check $(SIRIL_SRC)

.PHONY: clean
clean:
	rm -f sref sref.siril.o
