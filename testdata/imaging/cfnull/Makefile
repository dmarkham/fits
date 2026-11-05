# cfitsio null/scaling truth for the imaging samples.
#   make          build ./cfnull (gitignored)
#   make goldens  write golden/<sample>.json for every ../comp_*.fz and
#                 cards/<sample>.json for every ../card_*.fits (the
#                 BSCALE / BZERO / BLANK card variants)
#   make verify   regenerate into verify.out/ and cmp with the goldens
CC     ?= cc
CFLAGS ?= -O2 -Wall -std=c99
LDLIBS := -lcfitsio -lm
SAMPLES := $(sort $(wildcard ../comp_*.fz))
CARDS   := $(sort $(wildcard ../card_*.fits))

cfnull: cfnull.c
	$(CC) $(CFLAGS) -o cfnull cfnull.c $(LDLIBS)

.PHONY: goldens verify clean
goldens: cfnull
	mkdir -p golden cards
	@set -e; for f in $(SAMPLES); do b=$$(basename $$f); ./cfnull $$f > golden/$$b.json; done
	@set -e; for f in $(CARDS); do b=$$(basename $$f); ./cfnull $$f > cards/$$b.json; done
	@echo "goldens: $(words $(SAMPLES)) + $(words $(CARDS)) files"

verify: cfnull
	rm -rf verify.out; mkdir -p verify.out/golden verify.out/cards
	@set -e; for f in $(SAMPLES); do b=$$(basename $$f); ./cfnull $$f > verify.out/golden/$$b.json; cmp golden/$$b.json verify.out/golden/$$b.json; done
	@set -e; for f in $(CARDS); do b=$$(basename $$f); ./cfnull $$f > verify.out/cards/$$b.json; cmp cards/$$b.json verify.out/cards/$$b.json; done
	@set -e; test "$$(ls cards | wc -l)" -eq $(words $(CARDS))
	@echo "verify: all $(words $(SAMPLES)) + $(words $(CARDS)) goldens reproduce byte-for-byte"
	rm -rf verify.out

clean:
	rm -f cfnull
	rm -rf verify.out
