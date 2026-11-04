# Siril FITS->float truth harness (see sirilref.c, PROVENANCE.md).
#   make          build ./sirilref (gitignored build artifact)
#   make goldens  write golden/<sample>.json for every ../*.fits, ../*.fz
#   make verify   regenerate into verify.out/ and cmp with the goldens
#   make clean    remove the binary and verify.out/ only
CC     ?= cc
CFLAGS ?= -O2 -Wall -std=c99
LDLIBS := -lcfitsio -lm
SAMPLES := $(sort $(wildcard ../*.fits ../*.fz))

sirilref: sirilref.c
	$(CC) $(CFLAGS) -o sirilref sirilref.c $(LDLIBS)

.PHONY: goldens verify clean
goldens: sirilref
	mkdir -p golden
	@set -e; for f in $(SAMPLES); do b=$$(basename $$f); ./sirilref $$f > golden/$$b.json; done
	@echo "goldens: $(words $(SAMPLES)) files"

verify: sirilref
	rm -rf verify.out; mkdir -p verify.out
	@set -e; for f in $(SAMPLES); do b=$$(basename $$f); ./sirilref $$f > verify.out/$$b.json; cmp golden/$$b.json verify.out/$$b.json; done
	@echo "verify: all $(words $(SAMPLES)) goldens reproduce byte-for-byte"
	rm -rf verify.out

clean:
	rm -f sirilref
	rm -rf verify.out

# Check that the VERBATIM functions in sirilref.c still equal Siril's source.
SIRIL_SRC ?= $(error set SIRIL_SRC to an unpacked Siril 1.4.0 source tree)
.PHONY: verbatim
verbatim:
	@set -e; for fn in 'int fit_stats(fitsfile *fptr, float *mini, float *maxi) {' \
	  'static void convert_floats(int bitpix, float *data, size_t nbdata, gboolean verbose) {' \
	  'void manage_bitpix(fitsfile *fptr, int *bitpix, int *orig_bitpix) {'; do \
	  awk -v f="$$fn" 'index($$0, f)==1 {p=1} p {print} p && /^}/ {exit}' $(SIRIL_SRC)/src/io/image_format_fits.c > verbatim.a; \
	  awk -v f="$$fn" 'index($$0, f)==1 {p=1} p {print} p && /^}/ {exit}' sirilref.c > verbatim.b; \
	  test -s verbatim.a; cmp verbatim.a verbatim.b; done; rm -f verbatim.a verbatim.b; \
	  echo "verbatim: fit_stats, convert_floats, manage_bitpix identical to $(SIRIL_SRC)"
