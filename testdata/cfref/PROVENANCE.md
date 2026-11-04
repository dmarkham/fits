# testdata/cfref: cfitsio truth harness

`cfref.c` is a small C program linked against the system libcfitsio. Its
`dump` mode writes a JSON description of a FITS file as cfitsio reads it.
The Go test suite treats that JSON as the independent truth: HDU byte
addresses, every header card, every image pixel and every table cell.
The field meanings are documented at the top of `cfref.c`, and the typed Go
decoder is `internal/testutil` (`LoadCfrefGolden`, `CfrefDump`).

## Golden format: explicit values, or an exact digest past 1024 values

An array of more than 1024 values (an image's `phys` / `i64` / `u64` /
`raw` pixels, or a fixed-width table column's values over all rows) is
not written element by element. It is stored as a digest object:

```json
{"n": 90000, "sha256": "<hex>", "nulls": [[start, count], ...],
 "sample": [[index, value], ...]}
```

- `sha256` is taken over a canonical big-endian encoding of every element
  in order (a flag byte for cfitsio's null flag, then the IEEE bits of a
  double with every NaN as `0x7ff8000000000000`, the 8-byte integer, or a
  length-prefixed string; the exact layout is in `cfref.c`, `jdigest`).
  A single differing element, null flag, sign of zero or NaN-vs-number
  changes it.
- `nulls` lists every null element as runs.
- `sample` keeps explicit values at a deterministic index set (the first
  and last 8, the image corners of the first and last plane, the first
  minimum and maximum, 64 strided indices), so a failure names the
  indices that differ.
- A large table column is written as `"cells": {"per_row": k, ...}`
  instead of `"rows"`, one digest per value kind over all rows in row
  order; a BITPIX -32 image with digested pixels also carries
  `phys_f32` (the same values rounded to float32 by a C cast).

Arrays of up to 1024 values stay explicit; that keeps the 960-pixel
synthetic compression corpus explicit, because its tests relate cfitsio's
decode to the funpack twin and the source element by element. Per-row
cells are always explicit.

The Go side is `internal/testutil` (`cfref_values.go`): `Values.Check` /
`Diff` compare a result exactly in either form (count, null runs, every
sample, sha256 over all elements), and `Values.Resolve` /
`CfrefColumn.Resolved` / `CfrefImage.PhysAsFloat32` hand a test the
golden's values one by one after proving the test's own reading equals
the digest exactly. `TestCfrefDigestCatchesOneChange` shows that one
changed pixel or table cell of a real decode fails;
`TestCfrefDigestsConsistent` checks every committed digest and that no
array over 1024 values is left explicit.

Size: the 267 goldens here went from 17.0 MB (every value explicit) to
5.1 MB. The same format applies to every other cfref dump in the
repository (`testdata/images/golden/*.fits.json`, `testdata/tables/golden`,
`testdata/comp_bad/golden`, `header/testdata/golden/*.fits.json`).

## Pool: mid-size arrays stored once (`make goldens` here only)

`make goldens` runs `cfref dump -pool golden/pool`: an array of 65..1024
values (the 40x24 compression corpus, small images) stays explicit but is
written once to `golden/pool/<kind>-<sha256>.json` (exactly the JSON
array the dump would print), and the golden holds
`{"n": N, "sha256": "<hex>", "pool": "<kind>"}`. The funpack twins and
sources hold the same pixels as the compressed files, and an integer
image's phys/i64/raw are often equal, so 162 files hold every such array
of the 265 goldens. `internal/testutil` loads the file and checks kind,
count and sha256 before a test sees the values (`TestCfrefPoolComplete`,
`TestCfrefPoolRefusesBadFiles`); `make verify` and
`TestCompressedGoldensReproduce` regenerate the pool byte for byte.
golden/ went from 4.92 MB to 3.51 MB. Other cfref dumps (images/,
tables/, comp_bad/, header/) are made without `-pool`.

## Tool versions used for the committed goldens

| tool       | version                                   |
|------------|-------------------------------------------|
| libcfitsio | 4.7.0 (Arch package `cfitsio 1:4.7.0-1`; `/usr/lib/libcfitsio.so.4.7.0`, sha256 `946cc8ae…bbab65d`) |
| compiler   | gcc 16.2.1 (Arch `gcc 16.2.1+r23+gd564253eb6c8-1`), `-O2 -Wall -Wextra -std=c99` |
| host       | Linux x86_64 (LP64)                       |

Every golden records the cfitsio version it was made with (`"cfitsio"`), and
`internal/testutil` checks that it is `4.7.0`.

## Commands

```sh
make -C testdata/cfref            # build ./cfref (gitignored)
make -C testdata/cfref goldens    # write golden/<sample>.json for every
                                  # testdata/*.fits, *.fit, *.std
make -C testdata/cfref verify     # regenerate into verify.out/ and cmp
                                  # against the committed goldens
```

Each golden is exactly `./cfref dump ../<sample> > golden/<sample>.json`.
`make clean` removes only the binary and `verify.out/`; no target deletes
committed goldens.

Reproducibility: `make verify` and the Go test
`internal/testutil.TestCfrefDumpReproducesGoldens` (tier 2, `go test -tags
interop`) both require the live harness output to match the committed
golden byte-for-byte.

Readings of files the tests write (not committed samples) are not goldens
here: they are recorded per test in `testdata/recorded` (see
`docs/VALIDATION.md`, "Test tiers"). Re-record a test's readings after a
deliberate change to the bytes it writes (or to `cfref.c`, whose sha256
keys every reading) with
`FITS_RECORD=1 FITS_ASTROPY_PYTHON=$PWD/venv/bin/python go test -tags interop -run '^TestName$' ./<pkg>`;
the readings are written only for a test that passes against the live
tools.

## Samples covered

All 185 FITS files at the top level of `testdata/` (provenance of the
samples themselves is in `testdata/README.md` and
`testdata/COMP_FZ_PROVENANCE.md`): `testprog.std`, `testf77.std`,
`iter_a.fit`, `iter_b.fit`, `iter_c.fit`, `iter_image.fit`, `vari.fits`, the
36 astropy-written `comp_<name>.fits` / `comp_raw_<name>.fits` files (the
dither fixtures are `comp_rice_f32_dither1/dither2`), and the fpack /
libcfitsio corpus `comp_fz_*`, `comp_fzun_*` (funpack twins) and
`comp_fzsrc_*` (sources). `make verify` reproduces all 185 goldens byte for
byte. `TestGoldensCoverEveryTestdataSample` fails if a sample has no golden
or a golden has no sample.

## What cfitsio itself reports (recorded in the goldens, not hidden)

- `comp_nocompress_i16.fits` (ZCMPTYPE = 'NOCOMPRESS', written by astropy):
  cfitsio 4.7.0 cannot read the pixels: status 414 "error uncompressing
  image" ("unknown compression algorithm"); `funpack` fails the same way.
  The golden has `"phys": null` and the error in `"errors"`.
- Tile-compressed integer images: cfitsio's decompressor rejects TLONGLONG
  (status 410), so `cfref` reads `"i64"` for compressed images as TLONG,
  which is 64-bit on LP64.
- `iter_image.fit` HDU 0 has a valid CHECKSUM/DATASUM (`"checksum"` 1/1); no
  other sample carries checksum keywords (0/0).
- fitsverify (not part of the goldens, for reference): `testprog.std` 1
  warning, `testf77.std` 3 warnings, `iter_c.fit` 1 error, others clean.
- libcfitsio 4.7.0 segfaults when an HCOMPRESS_1 image with BZERO is
  re-read through the same handle after `fits_set_bscale` (the tile cache
  is re-used). The harness reads the `raw` values of a compressed image
  through a second handle on the file, which reads the same stored values
  without the crash (`comp_fz_u16sat_hcomp_s2` crashed it before; every
  other golden is unchanged by this). The same TLONG
  `fits_read_imgnull` crash (inside `imcomp_decompress_tile`) hits
  cfitsio's own HCOMPRESS-with-BZERO goldens in
  `compress/testdata/cref/img` and `testdata/imaging/comp_i16_bscale_neg` /
  `comp_i16_bscale0`; those areas read TDOUBLE through their own harnesses
  (`cref decode`, `cfnull`) instead (see their PROVENANCE.md).
- `fits_read_img` returns GZIP_COMPRESSED_DATA / UNCOMPRESSED_DATA fallback
  tiles of a scaled or BLANK image unscaled and without null checks, while
  funpack and astropy apply the image scaling. The goldens record cfitsio's
  values; `TestCompressedCorpusRead` treats the funpack twin as the truth
  for those rows (rawRows) and fails if cfitsio stops differing.
- Uninitialized memory: libcfitsio 4.7.0's TDOUBLE `fits_read_img` of an
  integer image with 8-byte GZIP tiles (and of a quantized float image with
  such tiles) returns its tile buffer without writing it, status 0. What a
  caller gets is whatever the heap held: the goldens used to show zeros by
  chance, and a harness change that moved allocations made them heap
  addresses. `cfref` therefore re-executes itself once with
  `GLIBC_TUNABLES=glibc.malloc.tcache_count=0:glibc.malloc.perturb=90`, so
  glibc fills every malloc'd block with 0xA5 bytes: those pixels read as the
  double 0xA5A5A5A5A5A5A5A5 on every run, and the tests pin exactly that.
  Every other golden is byte-identical with and without the fill (checked
  over every sample), so nothing else reads uninitialized memory. The
  harness also zero-fills (calloc) its own read buffers.
- Card pool: with `dump -pool DIR` each header card whose cfitsio parse
  depends on its 80 bytes alone (all but a string card ending in `&`, whose
  `str` takes the CONTINUE cards after it) is stored once as
  `DIR/card-<sha256 of the card text>.json` (the parse: every field but `n`
  and `raw`), and the golden holds the card text alone; `n` is its
  position. A second, different parse for the same text stops the dump.
  internal/testutil reads it back; `TestCfrefDumpReproducesGoldens`
  requires the golden with its cards and arrays inlined to equal the plain
  dump byte for byte, and `TestCfrefPoolComplete` that every pool file is
  referenced and every reference has its file.
- The funpack/astropy twins of the compressed read corpus are not
  committed: `make goldens` and `make verify` regenerate them
  (`../gen_twins.sh` into `twins.out`; needs funpack and, for the ZBITPIX
  64 twins, `FITS_ASTROPY_PYTHON`), and their goldens are
  `golden/twins/comp_fzun_<case>.fits.json`. The `comp_bad` goldens are
  pooled the same way and made by the same targets.
