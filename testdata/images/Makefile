# Image / checksum truth samples and goldens for the Go fits test suite.
#
#   make            -> build ./imgref (gitignored) and ../cfref/cfref
#   make all        -> regenerate every sample and golden (gen.py)
#   make verify     -> regenerate into verify.out/ and byte-compare with
#                      the committed samples and goldens
#   make clean      -> remove the binary and verify.out/ ONLY
#
# Requires libcfitsio (4.7.0 for the committed files) and the repo venv
# (astropy 8.0.2, numpy 2.4.6). See PROVENANCE.md.

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -std=c99
LDLIBS  := -lcfitsio -lm
PYTHON  ?= ../../venv/bin/python
DIRS    := cfwrite astropy_written astropy_data framing refused checksum keywords duplicate golden

imgref: imgref.c
	$(CC) $(CFLAGS) -o imgref imgref.c $(LDLIBS)

.PHONY: cfref
cfref:
	$(MAKE) -C ../cfref cfref

.PHONY: all
all: imgref cfref
	$(PYTHON) gen.py

.PHONY: verify
verify: imgref cfref
	rm -rf verify.out
	mkdir -p verify.out
	$(PYTHON) gen.py $(CURDIR)/verify.out
	@set -e; for d in $(DIRS); do \
	  $(PYTHON) cmp_golden.py $$d verify.out/$$d; \
	done; echo "verify: every sample and golden reproduces byte-for-byte (gzipped goldens compared decompressed)"
	rm -rf verify.out

.PHONY: clean
clean:
	rm -f imgref
	rm -rf verify.out
