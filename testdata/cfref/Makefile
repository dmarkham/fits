# cfitsio truth harness for the Go fits test suite.
#
#   make            -> build ./cfref (gitignored build artifact)
#   make goldens    -> (re)write golden/<file>.json for every FITS file in
#                      ../ (testdata/), using the system libcfitsio, and
#                      golden/twins/<twin>.json for the decompressed twins of
#                      the compressed read corpus, which ../gen_twins.sh
#                      writes into twins.out (funpack; astropy at ZBITPIX 64:
#                      set FITS_ASTROPY_PYTHON); mid-size value arrays go to
#                      golden/pool once each (cfref.c jarr), and so do the
#                      header cards (card-<sha256>.json, cfref.c dump_cards);
#                      the pool files are committed gzipped (gzip -9n,
#                      <kind>-<sha256>.json.gz; verify compares the
#                      decompressed bytes; testutil.ReadGolden reads them);
#                      also ../comp_bad/golden/<file>.json for ../comp_bad
#   make verify     -> regenerate into a scratch dir and byte-compare with
#                      the committed goldens (reproducibility check)
#   make clean      -> remove the binary and scratch output ONLY; committed
#                      goldens are never deleted by any target
#
# Requires: libcfitsio headers + library (Arch: extra/cfitsio). Goldens were
# generated with cfitsio 4.7.0; see PROVENANCE.md.

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -std=c99
LDLIBS  := -lcfitsio -lm

# Every FITS sample committed in testdata/ (the directory above this one).
SAMPLES := $(sort $(wildcard ../*.fits ../*.fit ../*.std))
# The twins (not committed) are regenerated here by ../gen_twins.sh.
TWINS   := twins.out
# The refused layouts of the compressed read corpus; their goldens are
# ../comp_bad/golden/<file>.json (pooled like the others).
BAD     := $(sort $(wildcard ../comp_bad/*.fits))

cfref: cfref.c
	$(CC) $(CFLAGS) -o cfref cfref.c $(LDLIBS)

.PHONY: goldens
goldens: cfref
	rm -rf golden/pool golden/twins
	mkdir -p golden/pool golden/twins
	sh ../gen_twins.sh $(TWINS)
	@set -e; for f in $(SAMPLES); do \
	  b=$$(basename $$f); \
	  echo "cfref dump -pool golden/pool $$f > golden/$$b.json"; \
	  ./cfref dump -pool golden/pool $$f > golden/$$b.json.tmp; \
	  mv golden/$$b.json.tmp golden/$$b.json; \
	done; for f in $(TWINS)/*.fits; do \
	  b=$$(basename $$f); \
	  ./cfref dump -pool golden/pool $$f > golden/twins/$$b.json.tmp; \
	  mv golden/twins/$$b.json.tmp golden/twins/$$b.json; \
	done; for f in $(BAD); do \
	  b=$$(basename $$f); \
	  ./cfref dump -pool golden/pool $$f > ../comp_bad/golden/$$b.json.tmp; \
	  mv ../comp_bad/golden/$$b.json.tmp ../comp_bad/golden/$$b.json; \
	done
	gzip -9n golden/pool/*.json

.PHONY: verify
verify: cfref
	rm -rf verify.out
	mkdir -p verify.out/pool verify.out/twins verify.out/bad
	sh ../gen_twins.sh verify.out/twins.fits
	@set -e; for f in $(SAMPLES); do \
	  b=$$(basename $$f); \
	  ./cfref dump -pool verify.out/pool $$f > verify.out/$$b.json; \
	  cmp golden/$$b.json verify.out/$$b.json; \
	done; for f in verify.out/twins.fits/*.fits; do \
	  b=$$(basename $$f); \
	  ./cfref dump -pool verify.out/pool $$f > verify.out/twins/$$b.json; \
	done; for f in $(BAD); do \
	  b=$$(basename $$f); \
	  ./cfref dump -pool verify.out/pool $$f > verify.out/bad/$$b.json; \
	done; diff -r golden/twins verify.out/twins; \
	diff -r ../comp_bad/golden verify.out/bad; \
	cmp ../comp_fzun_same.txt verify.out/twins.fits/comp_fzun_same.txt; \
	n=0; for f in verify.out/pool/*.json; do \
	  gzip -dc golden/pool/$$(basename $$f).gz | cmp - $$f; n=$$((n+1)); \
	done; [ $$n -eq $$(ls golden/pool | wc -l) ] || { echo "pool: $$n regenerated, $$(ls golden/pool | wc -l) committed"; exit 1; }; \
	echo "verify: all $(words $(SAMPLES)) goldens, the twin goldens and the pool reproduce byte-for-byte"
	rm -rf verify.out

.PHONY: clean
clean:
	rm -f cfref
	rm -rf verify.out $(TWINS)
