# cfitsio truth for HDUByName (see cfmovnam.c): duplicate EXTNAME cards
# (dupext.fits) and the names where HDUByName deliberately differs from
# fits_movnam_hdu (names.fits).
#
#   make          -> build ./cfmovnam, ./mkdupext and ./mknames (gitignored)
#   make goldens  -> write dupext.fits, names.fits and golden/{dupext,names}.json
#   make verify   -> regenerate into verify.out/ and cmp with the committed files
#   make clean    -> remove the binaries and verify.out/ only
#
# NAMES / DIFFNAMES are the names looked up; the tests
# (TestHDUByNameDuplicateEXTNAMEMatchesCfitsio,
# TestHDUByNameDifferencesFromCfitsio in hdubyname_test.go) read every name
# from the golden.

CC     ?= cc
CFLAGS ?= -O2 -Wall -Wextra -std=c99
LDLIBS := -lcfitsio -lm

NAMES := SCI ERR DQ TAB LAST NONE
DIFFNAMES := SCI sci Sci ERR NUM 5 NONE

all: cfmovnam mkdupext mknames

cfmovnam: cfmovnam.c
	$(CC) $(CFLAGS) -o cfmovnam cfmovnam.c $(LDLIBS)

mkdupext: mkdupext.c
	$(CC) $(CFLAGS) -o mkdupext mkdupext.c $(LDLIBS)

mknames: mknames.c
	$(CC) $(CFLAGS) -o mknames mknames.c $(LDLIBS)

.PHONY: goldens
goldens: cfmovnam mkdupext mknames
	./mkdupext dupext.fits
	./mknames names.fits
	mkdir -p golden
	./cfmovnam dupext.fits $(NAMES) > golden/dupext.json
	./cfmovnam names.fits $(DIFFNAMES) > golden/names.json

.PHONY: verify
verify: cfmovnam mkdupext mknames
	rm -rf verify.out && mkdir -p verify.out
	./mkdupext verify.out/dupext.fits
	cmp dupext.fits verify.out/dupext.fits
	./mknames verify.out/names.fits
	cmp names.fits verify.out/names.fits
	./cfmovnam dupext.fits $(NAMES) > verify.out/dupext.json
	cmp golden/dupext.json verify.out/dupext.json
	./cfmovnam names.fits $(DIFFNAMES) > verify.out/names.json
	cmp golden/names.json verify.out/names.json
	@echo "verify: dupext.fits, names.fits and their goldens reproduce byte-for-byte"
	rm -rf verify.out

.PHONY: clean
clean:
	rm -f cfmovnam mkdupext mknames
	rm -rf verify.out
