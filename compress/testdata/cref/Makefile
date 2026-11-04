# C reference harness: the Go compress package and the root package's
# compressed-image writer are diffed against the real libcfitsio.
#
#   make            -> build ./cref (gitignored build artifact)
#   make goldens    -> regenerate EVERY committed golden in place:
#                        ./cref gen    golden          quantizer/noise goldens
#                        ./cref encode golden          writer-parity goldens
#                                                     (rice/ hcomp/ plio/ img/
#                                                     imgw/ img64/ img2/ and
#                                                     golden/MANIFEST)
#                      then gzips every cfitsio-written FITS file
#                      (<case>.fits -> <case>.fits.gz, gzip -n -9: no name,
#                      no timestamp); the tests inflate them
#   make verify     -> regenerate into verify.out/ and byte-compare with the
#                      committed goldens, the .fits.gz files inflated first
#                      (so the check is on the bytes cfitsio wrote, whatever
#                      gzip produced the committed .gz)
#   make clean      -> remove the binary and verify.out/ ONLY. Goldens are
#                      committed test truth and are never deleted by any
#                      target; `git diff golden/` shows what a regeneration
#                      changed.
#
# Requires: libcfitsio 4.7.0 headers + library (Arch: extra/cfitsio), zlib.
# See PROVENANCE.md for the exact versions the committed goldens came from.

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -Wno-missing-field-initializers -std=c99
LDLIBS  := -lcfitsio -lm

ENCDIRS := rice hcomp plio img imgw img64 img2

cref: cref.c
	$(CC) $(CFLAGS) -o cref cref.c $(LDLIBS)

.PHONY: goldens
goldens: cref
	mkdir -p golden $(addprefix golden/,$(ENCDIRS))
	./cref gen golden
	./cref encode golden
	find golden -name '*.fits' -exec gzip -n -9 -f {} +

.PHONY: verify
verify: cref
	rm -rf verify.out verify.ref
	mkdir -p verify.out $(addprefix verify.out/,$(ENCDIRS))
	./cref gen verify.out > /dev/null
	./cref encode verify.out > /dev/null
	cp -r golden verify.ref
	find verify.ref -name '*.fits.gz' -exec gzip -d {} +
	diff -r verify.ref verify.out
	@echo "verify: every golden reproduces byte-for-byte"
	rm -rf verify.out verify.ref

.PHONY: clean
clean:
	rm -f cref
	rm -rf verify.out verify.ref
