# Table samples and their cfitsio goldens.
#
#   make samples   -> (re)write the cf_*.fits samples with gen_tables.c
#                     (system libcfitsio), gaia_dr3_sample.fits with
#                     gen_gaia.py and vla_theap_gap.fits with gen_theap.py
#                     (astropy; PYTHON or FITS_ASTROPY_PYTHON picks the venv),
#                     strarray_astropy.fits with gen_strarray.py and
#                     astropy_only/astropy_special.fits with gen_special.py
#   make goldens   -> (re)write golden/<file>.json for every *.fits here
#                     with the shared cfitsio harness ../cfref/cfref, and
#                     the cfitsio + astropy goldens of refused/*.fits
#   make tfloat-goldens -> (re)write golden_tfloat/<file>.json, cfitsio's
#                     float32 (TFLOAT) read of every numeric column of every
#                     *.fits here (cf_tfloat.c)
#   make heap-goldens -> (re)write heap/fitsverify.txt (fitsverify -q on the
#                     heap/ boundary samples written by heap/gen_heap.py)
#   make astropy-goldens -> (re)write golden_astropy/<file>.json (astropy's
#                     read, astropy_dump.py --golden) for every table sample:
#                     *.fits here, astropy_only/*.fits and the cfitsio test
#                     files in ../ that the table tests read
#   make verify    -> regenerate samples and every golden set into
#                     verify.out/ and byte-compare with the committed files
#   make clean     -> remove build artifacts only (never samples/goldens)
#
# The generator binaries are named gen_tables.out and cf_tfloat.out so the
# repository's existing "*.out" ignore rule keeps them out of git.

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -std=c99 -D_DEFAULT_SOURCE
LDLIBS  := -lcfitsio -lm
PYTHON  ?= $(or $(FITS_ASTROPY_PYTHON),../../venv/bin/python)
CFREF   := ../cfref/cfref
# The cfitsio goldens are stored compact (no whitespace between tokens,
# every token byte for byte): compact_json.py is plain Python.
COMPACT := python3 compact_json.py

CF_SAMPLES := cf_bintable_fixed.fits cf_bintable_uint64.fits cf_vla_p.fits \
              cf_vla_q.fits cf_ascii.fits cf_ascii_implied.fits \
              cf_bintable_kzero.fits cf_bintable_substr.fits cf_ascii_scaled.fits \
              cf_bintable_ranges.fits cf_bintable_tzero_big.fits cf_ascii_exact.fits \
              cf_bintable_a3dtable.fits cf_bintable_3dtable.fits \
              cf_bintable_badlogical.fits cf_bintable_scaled_exact.fits \
              cf_ascii_scaled_exact.fits cf_bintable_strarray.fits \
              cf_ascii_scaled_edge.fits cf_bintable_scaled_special.fits \
              cf_ascii_tzero_round.fits cf_bintable_tzero_round.fits
PY_SAMPLES := gaia_dr3_sample.fits vla_theap_gap.fits strarray_astropy.fits
# Table samples outside this directory that the table tests also read.
TOP_SAMPLES := ../testprog.std ../testf77.std ../iter_a.fit ../iter_b.fit ../iter_c.fit
ASTROPY_SAMPLES = $(sort $(wildcard *.fits)) astropy_only/astropy_special.fits $(TOP_SAMPLES)

gen_tables.out: gen_tables.c
	$(CC) $(CFLAGS) -o $@ gen_tables.c $(LDLIBS)

cf_tfloat.out: cf_tfloat.c
	$(CC) $(CFLAGS) -o $@ cf_tfloat.c $(LDLIBS)

$(CFREF):
	$(MAKE) -C ../cfref

.PHONY: tfloat-goldens
tfloat-goldens: cf_tfloat.out
	mkdir -p golden_tfloat
	@set -e; for f in $(sort $(wildcard *.fits)); do \
	  ./cf_tfloat.out $$f > golden_tfloat/$$f.json.tmp; \
	  mv golden_tfloat/$$f.json.tmp golden_tfloat/$$f.json; \
	done

.PHONY: samples
samples: gen_tables.out
	mkdir -p refused
	./gen_tables.out .
	$(PYTHON) gen_gaia.py gaia_dr3_sample.fits
	$(PYTHON) gen_theap.py vla_theap_gap.fits
	$(PYTHON) gen_strarray.py strarray_astropy.fits
	$(PYTHON) astropy_only/gen_special.py astropy_only/astropy_special.fits
	$(PYTHON) heap/gen_heap.py heap

.PHONY: heap-goldens
heap-goldens:
	cd heap && for f in heap_*.fits; do fitsverify -q $$f; done > fitsverify.txt 2>&1

.PHONY: astropy-goldens
astropy-goldens:
	mkdir -p golden_astropy
	@set -e; for f in $(ASTROPY_SAMPLES); do \
	  b=$$(basename $$f); echo "astropy_dump.py --golden $$f > golden_astropy/$$b.json"; \
	  $(PYTHON) astropy_dump.py --golden $$f > golden_astropy/$$b.json.tmp; \
	  mv golden_astropy/$$b.json.tmp golden_astropy/$$b.json; \
	done

.PHONY: goldens
goldens: $(CFREF)
	$(CFREF) dump refused/cf_bintable_cplx_scaled.fits > refused/cplx.raw.tmp
	$(COMPACT) < refused/cplx.raw.tmp > refused/cf_bintable_cplx_scaled.fits.cfref.json
	rm -f refused/cplx.raw.tmp
	$(PYTHON) astropy_dump.py --golden refused/cf_bintable_cplx_scaled.fits > refused/cf_bintable_cplx_scaled.fits.astropy.json
	mkdir -p golden
	@set -e; for f in $(sort $(wildcard *.fits)); do \
	  echo "cfref dump $$f | compact > golden/$$f.json"; \
	  $(CFREF) dump $$f > golden/$$f.json.raw; \
	  $(COMPACT) < golden/$$f.json.raw > golden/$$f.json.tmp; rm -f golden/$$f.json.raw; \
	  mv golden/$$f.json.tmp golden/$$f.json; \
	done

.PHONY: verify
verify: gen_tables.out cf_tfloat.out $(CFREF)
	rm -rf verify.out
	mkdir -p verify.out/refused
	./gen_tables.out verify.out
	cmp refused/cf_bintable_cplx_scaled.fits verify.out/refused/cf_bintable_cplx_scaled.fits
	$(CFREF) dump refused/cf_bintable_cplx_scaled.fits > verify.out/refused/cplx.raw
	$(COMPACT) < verify.out/refused/cplx.raw > verify.out/refused/cplx.cfref.json
	cmp refused/cf_bintable_cplx_scaled.fits.cfref.json verify.out/refused/cplx.cfref.json
	$(PYTHON) astropy_dump.py --golden refused/cf_bintable_cplx_scaled.fits > verify.out/refused/cplx.astropy.json
	cmp refused/cf_bintable_cplx_scaled.fits.astropy.json verify.out/refused/cplx.astropy.json
	$(PYTHON) gen_gaia.py verify.out/gaia_dr3_sample.fits
	$(PYTHON) gen_theap.py verify.out/vla_theap_gap.fits
	$(PYTHON) gen_strarray.py verify.out/strarray_astropy.fits
	$(PYTHON) astropy_only/gen_special.py verify.out/astropy_special.fits
	mkdir -p verify.out/heap
	$(PYTHON) heap/gen_heap.py verify.out/heap
	@set -e; for f in heap/heap_*.fits; do cmp $$f verify.out/$$f; done
	cd verify.out/heap && for f in heap_*.fits; do fitsverify -q $$f; done > fitsverify.txt 2>&1
	cmp heap/fitsverify.txt verify.out/heap/fitsverify.txt
	@set -e; for f in $(CF_SAMPLES) $(PY_SAMPLES); do \
	  cmp $$f verify.out/$$f; \
	done; cmp astropy_only/astropy_special.fits verify.out/astropy_special.fits; \
	echo "verify: generated samples reproduce byte-for-byte"
	@set -e; for f in $(sort $(wildcard *.fits)); do \
	  $(CFREF) dump $$f > verify.out/$$f.raw; \
	  $(COMPACT) < verify.out/$$f.raw > verify.out/$$f.json; \
	  cmp golden/$$f.json verify.out/$$f.json; \
	done; echo "verify: goldens reproduce byte-for-byte"
	@set -e; for f in $(sort $(wildcard *.fits)); do \
	  ./cf_tfloat.out $$f > verify.out/$$f.tfloat.json; \
	  cmp golden_tfloat/$$f.json verify.out/$$f.tfloat.json; \
	done; test "$$(ls golden_tfloat | wc -l)" = "$$(ls *.fits | wc -l)"; \
	echo "verify: TFLOAT goldens reproduce byte-for-byte"
	@set -e; for f in $(ASTROPY_SAMPLES); do \
	  b=$$(basename $$f); \
	  $(PYTHON) astropy_dump.py --golden $$f > verify.out/$$b.astropy.json; \
	  cmp golden_astropy/$$b.json verify.out/$$b.astropy.json; \
	done; echo "verify: astropy goldens reproduce byte-for-byte"
	rm -rf verify.out

.PHONY: clean
clean:
	rm -f gen_tables.out cf_tfloat.out
	rm -rf verify.out
