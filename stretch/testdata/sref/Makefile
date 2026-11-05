# Siril 1.4.0 stretch oracle (see PROVENANCE.md in ../).
#
#   make                 build sref, sref_libm and sref_clahe (gitignored binaries)
#   make goldens         rebuild + regenerate the Siril goldens in ../golden
#   make verify          regenerate into .verify-golden and diff with ../golden
#   make verify-excerpt  re-extract siril/ from the 1.4.0 source and diff
#   make libm-probe      glibc float libm vs correctly rounded (GHT tolerance)
#
# Flags: plain -O2 for x86-64 (no -march, no -ffast-math, no OpenMP), which is
# what Siril's release build uses for these files; -ffp-contract=off pins
# "no fused multiply-add" explicitly.
SHELL   := /bin/bash
CC      ?= cc
CXX     ?= c++
CFLAGS  := -O2 -std=gnu11 -ffp-contract=off -Ishim -I. -Isiril $(shell pkg-config --cflags glib-2.0)
CXXFLAGS:= -O2 -std=gnu++17 -ffp-contract=off -Ishim -I. -Isiril $(shell pkg-config --cflags glib-2.0 opencv5)
# Python with numpy: $(FITS_ASTROPY_PYTHON) if set (the variable the Go
# interop tests use), else the repository venv, else python3 on PATH.
PY      ?= $(firstword $(FITS_ASTROPY_PYTHON) $(wildcard $(abspath ../../../venv/bin/python)) python3)
SIRIL_SRC ?= $(error set SIRIL_SRC to the src/ directory of an unpacked Siril 1.4.0 source tree)

all: sref sref_libm sref_clahe

sref: sref.c glue.c siril/filters/mtf.c siril/filters/ght.c siril/filters/asinh_excerpt.c siril/stats_excerpt.c siril/rt_algo_excerpt.cc
	$(CC) $(CFLAGS) -c sref.c -o sref.o
	$(CC) $(CFLAGS) -c glue.c -o glue.o
	$(CC) $(CFLAGS) -c siril/filters/mtf.c -o mtf.o
	$(CC) $(CFLAGS) -c siril/filters/ght.c -o ght.o
	$(CC) $(CFLAGS) -c siril/filters/asinh_excerpt.c -o asinh.o
	$(CXX) $(CXXFLAGS) -c siril/rt_algo_excerpt.cc -o rt_algo.o
	$(CXX) -o $@ sref.o glue.o mtf.o ght.o asinh.o rt_algo.o $(shell pkg-config --libs glib-2.0) -lm
	rm -f *.o

# The same driver with the verbatim ght.c built against double-rounded
# logf/powf/expf (shim/libm_double.h). Only used to locate the pixels where
# glibc's float libm changes the GHT output; never a source of expected
# values.
sref_libm: sref.c glue.c siril/filters/mtf.c siril/filters/ght.c siril/filters/asinh_excerpt.c siril/stats_excerpt.c siril/rt_algo_excerpt.cc shim/libm_double.h
	$(CC) $(CFLAGS) -c sref.c -o sref_l.o
	$(CC) $(CFLAGS) -c glue.c -o glue_l.o
	$(CC) $(CFLAGS) -c siril/filters/mtf.c -o mtf_l.o
	$(CC) $(CFLAGS) -include shim/libm_double.h -c siril/filters/ght.c -o ght_l.o
	$(CC) $(CFLAGS) -c siril/filters/asinh_excerpt.c -o asinh_l.o
	$(CXX) $(CXXFLAGS) -c siril/rt_algo_excerpt.cc -o rt_algo_l.o
	$(CXX) -o $@ sref_l.o glue_l.o mtf_l.o ght_l.o asinh_l.o rt_algo_l.o $(shell pkg-config --libs glib-2.0) -lm
	rm -f *_l.o

sref_clahe: sref_clahe.cpp siril/opencv_excerpt.cpp
	$(CXX) $(CXXFLAGS) -o $@ sref_clahe.cpp $(shell pkg-config --libs glib-2.0) -lopencv_imgproc -lopencv_core

libm_probe: libm_probe.c
	$(CC) -O2 -o $@ libm_probe.c -lm

.PHONY: goldens verify verify-excerpt clean libm-probe
libm-probe: libm_probe
	./libm_probe

goldens: sref sref_libm sref_clahe
	@$(PY) -c 'import numpy' || { echo 'PY=$(PY) has no numpy; pass PY=<venv python>' >&2; exit 1; }
	cd .. && $(PY) gen_siril_goldens.py

# Byte-for-byte reproduce check of every file gen_siril_goldens.py writes
# (TestSirilGoldensReproduce does the same from go test); the gzipped
# sweep goldens are compared decompressed.
verify: sref sref_libm sref_clahe
	rm -rf .verify-golden
	cd .. && $(PY) gen_siril_goldens.py --out sref/.verify-golden >/dev/null
	cd .verify-golden && for f in $$(find . -type f | sort); do \
	  case $$f in *.gz) gzip -dc $$f | cmp -s - <(gzip -dc ../../golden/$$f) ;; *) cmp -s $$f ../../golden/$$f ;; esac || { echo "DIFFERS: $$f"; exit 1; }; \
	done
	cd ../golden && for f in sref/*.json sref/MANIFEST sweep/*.json.gz $$(ls stretch_*.json | grep -v clahe_cv2) $$(ls inputs/* | grep -v /cv2_); do test -f ../sref/.verify-golden/$$f || { echo "NOT REGENERATED: $$f"; exit 1; }; done
	rm -rf .verify-golden
	@echo "every Siril golden reproduces byte for byte"

verify-excerpt:
	rm -rf .verify && mkdir .verify && cp extract.sh .verify/ && cd .verify && SIRIL_SRC=$(SIRIL_SRC) ./extract.sh >/dev/null
	diff -r siril .verify/siril && echo "siril/ matches $(SIRIL_SRC)"
	rm -rf .verify

clean:
	rm -f sref sref_libm sref_clahe libm_probe *.o
