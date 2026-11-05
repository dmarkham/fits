# Header-package samples and cfitsio goldens.
#
#   make samples   -> rewrite the hand-built samples (gen_samples.py: malformed,
#                     continue, adversarial) and the
#                     cfitsio-written sample cfwritten.fits (hdrref write)
#   make goldens   -> rewrite every golden in golden/ from the samples:
#                       golden/<s>.json          cfref dump <s>   (shared harness)
#                       golden/<s>.longstr.json  hdrref longstr <s>
#                       golden/testdata/<t>.longstr.json for every repo
#                                                testdata/ sample
#                       golden/testkey.json      hdrref testkey <KEYS>
#   make verify    -> regenerate samples and goldens into verify.out/ and
#                     byte-compare with the committed files
#   make clean     -> remove the hdrref binary and verify.out/ only
#
# Requires libcfitsio (4.7.0 for the committed goldens), a C compiler and
# python3 (stdlib only). The cfref harness is built from testdata/cfref.
# See PROVENANCE.md.

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -std=c99
LDLIBS  := -lcfitsio -lm
PYTHON  ?= python3
CFREF   := ../../testdata/cfref/cfref

SAMPLES  := malformed.fits continue.fits adversarial.fits cfwritten.fits
TDSAMPLES := $(sort $(wildcard ../../testdata/*.fits ../../testdata/*.fit ../../testdata/*.std))

# Keyword names whose legality cfitsio's fits_test_keyword decides
# (golden/testkey.json). Upper case: the Go header package upper-cases
# keyword names before validating them.
KEYS := SIMPLE NAXIS1 DATE-OBS A_B A-B 1234567 12345678 LOWER \
        'A B' 'A=B' 'A.B' 'A+B' 'A/B' "A'B" 'A&B' 'É' 'TAB	X' 'A*'

hdrref/hdrref: hdrref/hdrref.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

$(CFREF):
	$(MAKE) -C ../../testdata/cfref cfref

.PHONY: samples
samples: hdrref/hdrref
	$(PYTHON) gen_samples.py
	./hdrref/hdrref write cfwritten.fits

# out=<dir> selects where goldens go (default golden/).
out ?= golden
.PHONY: goldens
goldens: hdrref/hdrref $(CFREF)
	mkdir -p $(out)/testdata
	@set -e; for s in $(SAMPLES); do \
	  echo "golden $$s"; \
	  $(CFREF) dump $$s > $(out)/$$s.json.tmp; mv $(out)/$$s.json.tmp $(out)/$$s.json; \
	  ./hdrref/hdrref longstr $$s > $(out)/$$s.longstr.json.tmp; mv $(out)/$$s.longstr.json.tmp $(out)/$$s.longstr.json; \
	done
	@set -e; for f in $(TDSAMPLES); do \
	  b=$$(basename $$f); \
	  ./hdrref/hdrref longstr $$f > $(out)/testdata/$$b.longstr.json.tmp; \
	  mv $(out)/testdata/$$b.longstr.json.tmp $(out)/testdata/$$b.longstr.json; \
	done
	./hdrref/hdrref testkey $(KEYS) > $(out)/testkey.json.tmp
	mv $(out)/testkey.json.tmp $(out)/testkey.json

.PHONY: verify
verify: hdrref/hdrref $(CFREF)
	rm -rf verify.out && mkdir -p verify.out/samples
	cd verify.out/samples && cp ../../gen_samples.py . && $(PYTHON) gen_samples.py >/dev/null && rm gen_samples.py
	./hdrref/hdrref write verify.out/samples/cfwritten.fits
	@set -e; for s in $(SAMPLES); do test -f verify.out/samples/$$s || { echo "verify: $$s was not regenerated"; exit 1; }; cmp $$s verify.out/samples/$$s; done
	$(MAKE) --no-print-directory goldens out=verify.out/golden >/dev/null
	diff -r golden verify.out/golden
	rm -rf verify.out
	@echo "verify: samples and goldens reproduce byte-for-byte"

.PHONY: clean
clean:
	rm -f hdrref/hdrref
	rm -rf verify.out
