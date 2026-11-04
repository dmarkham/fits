# wcslib reference harnesses for validating the Go port against the C
# implementation. Analogous to compress/testdata/cref.
#
#   make              -> build wref and hdr2sky
#   make goldens      -> regenerate ALL committed goldens (wref + hdr)
#   make wref-goldens -> ./wref gen golden   (bare-projection goldens,
#                        golden/<case>.json + golden/<case>.bin) and
#                        ./wref trig libm (libm/wcstrig_guards.json)
#   make verify       -> regenerate every wref and hdr golden into
#                        verify.out/ and require them byte for byte
#                        (the .bin compressed bytes need the zlib recorded
#                        in each .json; TestWrefGoldensReproduce compares
#                        the decoded payload when the zlib differs)
#   make hdr-goldens  -> ./hdr2sky hdrcases/<name>.case | gzip -9n >
#                        golden/hdr/<name>.json.gz for every case file (full
#                        header pipeline goldens; committed gzipped for the
#                        repository size budget, verify and
#                        TestHdrGoldensReproduce compare the decompressed
#                        bytes)
#   make clean        -> remove the BINARIES only. Goldens are committed
#                        test truth and are never deleted by clean; they are
#                        overwritten in place by the *-goldens targets, and
#                        `git diff golden/` shows what a regeneration changed.
#
# Requires: libwcs 8.10 (Arch: extra/wcslib) and zlib. See PROVENANCE.md.

CC      ?= cc
CFLAGS  ?= -O2 -Wall -Wextra -std=c99 -I/usr/include
LDLIBS  := -lwcs -lm
ZLIBS   := -lz

HDRCASES := $(sort $(wildcard hdrcases/*.case))

.PHONY: all
all: wref hdr2sky

wref: wref.c
	$(CC) $(CFLAGS) -o wref wref.c $(LDLIBS) $(ZLIBS)

hdr2sky: hdr2sky.c
	$(CC) $(CFLAGS) -o hdr2sky hdr2sky.c $(LDLIBS)

.PHONY: goldens wref-goldens hdr-goldens
goldens: wref-goldens hdr-goldens

wref-goldens: wref
	mkdir -p golden libm
	./wref gen golden
	./wref trig libm

# Each golden is written to a temp file and renamed only on success, so a
# failing case never leaves a truncated golden behind.
hdr-goldens: hdr2sky
	mkdir -p golden/hdr
	@set -e; for c in $(HDRCASES); do \
		n=$$(basename $$c .case); \
		echo "hdr2sky $$c -> golden/hdr/$$n.json.gz"; \
		./hdr2sky $$c > golden/hdr/$$n.json.tmp; \
		gzip -9n -c golden/hdr/$$n.json.tmp > golden/hdr/$$n.json.gz.tmp; \
		rm golden/hdr/$$n.json.tmp; \
		mv golden/hdr/$$n.json.gz.tmp golden/hdr/$$n.json.gz; \
	done

.PHONY: verify
verify: wref hdr2sky
	rm -rf verify.out && mkdir -p verify.out/hdr
	./wref gen verify.out >/dev/null
	./wref trig verify.out >/dev/null
	@set -e; for f in golden/*.json golden/*.bin libm/wcstrig_guards.json; do cmp $$f verify.out/$$(basename $$f); done
	@set -e; for c in $(HDRCASES); do \
		n=$$(basename $$c .case); \
		./hdr2sky $$c > verify.out/hdr/$$n.json; \
		gzip -dc golden/hdr/$$n.json.gz | cmp - verify.out/hdr/$$n.json; \
	done; \
	for g in golden/hdr/*.json.gz; do n=$$(basename $$g .json.gz); [ -f hdrcases/$$n.case ] || { echo "orphan $$g"; exit 1; }; done
	rm -rf verify.out
	@echo "verify: every wref golden and all $(words $(HDRCASES)) hdr goldens reproduced byte for byte"

.PHONY: clean
clean:
	rm -f wref hdr2sky
	rm -rf verify.out
