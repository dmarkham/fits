# cfitsio truth for EditFile HDU copy / reorder / skip (see cfcopy.c).
#
#   make          -> build ./cfcopy and ./mkcomments (gitignored)
#   make goldens  -> write golden/<scenario>.fits with fits_copy_hdu
#   make verify   -> regenerate into verify.out/ and cmp with the goldens
#   make clean    -> remove the binary and verify.out/ only
#
# Each scenario is "name:sample:hdu,hdu,...". The HDU list is the order in
# which the 0-based HDUs of ../<sample> are copied; omitted HDUs are
# skipped. A sample under cfcopy/golden/ is the golden of an earlier
# scenario (cube3_1_0 and naxis0_1_0 start from an IMAGE extension that
# fits_copy_hdu made), so its scenario must come after that one. Keep this list in sync with editCopyScenarios in
# editfile_copy_test.go.
#
# reserve.fits (committed) is a scenario input made by cfitsio too: HDUs
# 0, 2 and 6 of testprog.std copied with morekeys = 40, so every header
# carries reserved blank cards before END.
#
# comments.fits (committed) is the input of comments_1_0, written by
# ./mkcomments (cfitsio): COMMENT cards that do and do not match the
# prefixes fits_copy_hdu uses to drop the standard FITS comments.
#
# cube3.fits, cube4.fits and cube6.fits (committed, 5760 bytes each) are the
# NAXIS 3/4/6 cube inputs: funpack -C (cfitsio 4.7.0) of the compressed cube
# samples ../comp_fz_cube16_rice_t3.fits, ../comp_fz_cube4d16_rice.fits and
# ../comp_fz_cube6d16_gzip2.fits. They are the bytes of the funpack twins
# that used to be committed as ../comp_fzun_cube*.fits.

CC     ?= cc
CFLAGS ?= -O2 -Wall -Wextra -std=c99
LDLIBS := -lcfitsio -lm

SCENARIOS := \
	testprog_2_0_1-7:testprog.std:2,0,1,3,4,5,6,7 \
	testprog_skip0:testprog.std:1,2,3,4,5,6,7 \
	testprog_6_2:testprog.std:6,2 \
	testprog_0_3_0:testprog.std:0,3,0 \
	iter_image_0_0:iter_image.fit:0,0 \
	vari_skip0:vari.fits:1 \
	reserve_1_0_2:cfcopy/reserve.fits:1,0,2 \
	comments_1_0:cfcopy/comments.fits:1,0 \
	comp_skip0:comp_fz_u16_rice.fits:1 \
	comp_1_0:comp_fz_u16_rice.fits:1,0 \
	cube3_0_0:cfcopy/cube3.fits:0,0 \
	cube4_0_0:cfcopy/cube4.fits:0,0 \
	cube6_0_0:cfcopy/cube6.fits:0,0 \
	cube3_1_0:cfcopy/golden/cube3_0_0.fits:1,0 \
	naxis0_0_0:vari.fits:0,0 \
	naxis0_1_0:cfcopy/golden/naxis0_0_0.fits:1,0

all: cfcopy mkcomments

cfcopy: cfcopy.c
	$(CC) $(CFLAGS) -o cfcopy cfcopy.c $(LDLIBS)

mkcomments: mkcomments.c
	$(CC) $(CFLAGS) -o mkcomments mkcomments.c $(LDLIBS)

RESERVE_SRC := ../testprog.std 0 2 6
CUBES := cube3:cube16_rice_t3 cube4:cube4d16_rice cube6:cube6d16_gzip2

define cubes
	@set -e; for c in $(CUBES); do \
	  rm -f $(1)/$${c%%:*}.fits; funpack -C -O $(1)/$${c%%:*}.fits ../comp_fz_$${c#*:}.fits; \
	done
endef

define run
	@set -e; for s in $(SCENARIOS); do \
	  name=$${s%%:*}; rest=$${s#*:}; sample=$${rest%%:*}; hdus=$${rest#*:}; \
	  echo "cfcopy $(1)/$$name.fits ../$$sample $$(echo $$hdus | tr , ' ')"; \
	  ./cfcopy $(1)/$$name.fits ../$$sample $$(echo $$hdus | tr , ' '); \
	done
endef

.PHONY: goldens
goldens: cfcopy mkcomments
	./cfcopy -r 40 reserve.fits $(RESERVE_SRC)
	./mkcomments comments.fits
	$(call cubes,.)
	mkdir -p golden
	$(call run,golden)

.PHONY: verify
verify: cfcopy mkcomments
	rm -rf verify.out && mkdir -p verify.out
	./cfcopy -r 40 verify.out/reserve.fits $(RESERVE_SRC)
	cmp reserve.fits verify.out/reserve.fits
	./mkcomments verify.out/comments.fits
	cmp comments.fits verify.out/comments.fits
	$(call cubes,verify.out)
	@set -e; for c in $(CUBES); do cmp $${c%%:*}.fits verify.out/$${c%%:*}.fits; done
	$(call run,verify.out)
	@set -e; for f in golden/*.fits; do cmp $$f verify.out/$$(basename $$f); done
	@echo "verify: every golden reproduces byte-for-byte"
	rm -rf verify.out

.PHONY: clean
clean:
	rm -f cfcopy mkcomments
	rm -rf verify.out
