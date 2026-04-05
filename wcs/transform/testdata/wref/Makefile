# wcslib reference harness for validating the Go port of the projection
# math against the C implementation. Analogous to compress/testdata/cref.
#
#   make            -> build wref
#   make goldens    -> run wref gen golden/ to produce the canonical
#                      per-projection .json files the Go tests diff against
#   make clean      -> remove binary + goldens
#
# Requires: libwcs (Arch: extra/wcslib)

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -std=c99 -I/usr/include
LDLIBS  := -lwcs -lm

wref: wref.c
	$(CC) $(CFLAGS) -o wref wref.c $(LDLIBS)

.PHONY: goldens
goldens: wref
	mkdir -p golden
	./wref gen golden

.PHONY: clean
clean:
	rm -f wref
	rm -rf golden
