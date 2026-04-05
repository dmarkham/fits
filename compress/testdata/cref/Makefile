# C reference harness for validating the Go port of fits_img_stats_float
# and fits_quantize_float against the real cfitsio library.
#
#   make            -> build cref
#   make goldens    -> run cref gen golden/, producing the fixtures the
#                      Go test suite diffs against.
#   make clean      -> remove binary and goldens
#
# Requires: libcfitsio headers + library (Arch: extra/cfitsio)

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -std=c99
LDLIBS  := -lcfitsio -lm

cref: cref.c
	$(CC) $(CFLAGS) -o cref cref.c $(LDLIBS)

.PHONY: goldens
goldens: cref
	mkdir -p golden
	./cref gen golden

.PHONY: clean
clean:
	rm -f cref
	rm -rf golden
