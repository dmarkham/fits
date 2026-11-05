# testdata/images: image, framing and checksum truth corpus

Samples and goldens for the image read/write, HDU layout, framing and
CHECKSUM tests of the root `fits` package (`image_*_truth_test.go`,
`image_write_readback_test.go`, `image_parsed_header_truth_test.go`,
`hdu_layout_truth_test.go`, `checksum_truth_test.go`) and of
`internal/checksum`. Every expected value comes from cfitsio or astropy;
none is produced by the Go library.

## Regenerating

```sh
make -C testdata/images                 # build ./imgref (gitignored) and ../cfref/cfref
make -C testdata/images all PYTHON=$PWD/venv/bin/python     # rewrite every sample + golden
make -C testdata/images verify PYTHON=$PWD/venv/bin/python  # regenerate into verify.out/, diff -r
```

`gen.py` is a pure function of its inputs: `make verify` regenerates
everything into a scratch directory and requires byte-identical output
(it did for the committed files).

## Tool versions

| tool       | version |
|------------|---------|
| libcfitsio | 4.7.0 (Arch `cfitsio 1:4.7.0-1`, `/usr/lib/libcfitsio.so.4.7.0` sha256 `946cc8ae…bbab65d`) |
| compiler   | gcc 16.2.1 (Arch `gcc 16.2.1+r23+gd564253eb6c8-1`), `-O2 -Wall -Wextra -std=c99` |
| astropy    | 8.0.2 (repo venv, pyenv Python 3.11.11) |
| numpy      | 2.4.6 |
| fitsverify | 4.22 (CFITSIO 4.7.0), used live by the tests, not by the generator |
| glibc      | 2.44 (`strtod` / `strtoll` behind cfitsio's `ffc2dd` / `ffc2jj`; `golden/strtod.txt`) |
| host       | Linux x86_64 (LP64) |

`golden/manifest.json` records the astropy/numpy versions, the sha256 of
every file copied from astropy's test data, cfitsio's write status for each
matrix case, and the list of samples the tests iterate.

## Harness

`imgref.c` (linked only against libcfitsio; usage at the top of the file):

* `write`: one image via `fits_create_img` (+ BSCALE/BZERO/BLANK cards and
  `fits_set_bscale`) and `fits_write_img`; reports cfitsio's status, so a
  NUM_OVERFLOW (412, cfitsio clamped a value) is recorded, not hidden.
* `chksum`: `fits_get_chksum` (data and HDU sums) and `fits_verify_chksum`
  per HDU (its failure is recorded as `"verify_error"`).
* `addchk`: copy a file and add DATASUM/CHECKSUM to every HDU with cfitsio's
  own sums (same steps as `fits_write_chksum`, fixed comments instead of a
  timestamp so the output is reproducible).
* `esum`: `fits_encode_chksum`.
* `pix`: `fits_read_imgnull` as TFLOAT for every image (each float32
  written as its exact double), as TULONGLONG for BITPIX 64 (cfitsio
  4.7.0's `fits_get_img_equivtype` reports LONGLONG_IMG for the uint64
  convention, so `cfref dump` only records an overflowing TLONGLONG read
  there), and, for integer images with null pixels, `"f64nn"`: TDOUBLE
  `fits_read_img` with `nulval` 0, i.e. null checking off, the scaled
  stored value of every pixel including the BLANK ones. An array of more
  than 1024 values is stored as an exact digest (sha256 over the canonical
  big-endian encoding, every null as a run, explicit samples), the same
  encoding and code as `../cfref/cfref.c` `jdigest` (copied verbatim);
  the tests resolve it with `testutil.Values`, so one differing element
  fails.
* `write`: when `fits_set_bscale` refuses the scaling it prints the status
  (`golden/write_refused.json`: ZERO_SCALE 322 for BSCALE 0 on BITPIX 8,
  16, 64 and -32). `golden/write_overflow.json` (`gen.py:WRITE_OVERFLOW`):
  each value of the overflow cases written ALONE with the case's
  BSCALE/BZERO, with cfitsio's status (412 NUM_OVERFLOW or 0) and the
  element it stored (hex); the truth for `ErrNumOverflow` in `WriteImage`
  and `OverwritePixels`.
* `strtod`: the C library's `strtod` and `strtoll` (base 10), which
  cfitsio's `ffc2dd` / `ffc2jj` call to convert BSCALE/BZERO/BLANK, on each
  line of stdin: `golden/strtod.txt` (3084 inputs: boundary cases plus
  seeded random tokens and numbers, `gen.py:gen_strtod`). glibc 2.44 on
  the generating host.
* `subset`: `fits_read_subset` for the cases in `gen.py:SUBSETS`.
* `overwrite`: copy a file and `fits_write_img` over its primary image, so
  cfitsio scales the values with the BSCALE/BZERO/BLANK it read from that
  header (`golden/keywords_overwrite.json`: every keywords/ sample, the
  written data unit in hex and cfitsio's status). The truth for
  `OverwritePixels`.
* `c2`: cfitsio's own `ffc2dd` / `ffc2jj` (fitsio2.h) on each line of stdin
  (`golden/c2.txt`: the strtod boundary inputs and random tokens plus
  wrapper cases: trailing blank, blank-separated junk, several `D`s, a `D`
  string over 72 characters).

Goldens are stored compact: `gen.py` removes the whitespace between JSON
tokens of the harness output (`compact`, which leaves every string and
number byte as the harness printed it and checks the parsed value is
unchanged) and writes its own JSON without indentation. A top-level
`testdata/` sample whose image HDUs are all tile-compressed (cfref golden
`"compressed"`) gets no `.pix.json`: no test reads those pixels here (the
compressed-read corpus proves them), and the corpus loader fails if such a
file reappears.

The shared `../cfref/cfref dump` supplies the per-file goldens
`golden/<dir>__<file>.json` (addresses, cards, TDOUBLE/TLONGLONG pixels,
null flags, checksum status). `golden/<dir>__<file>.chksum.json` and
`.pix.json` come from `imgref`; `golden/testdata__<file>.*.json` are the
same `imgref` reads for the top-level `testdata/` samples.

## Samples

| directory | what | source |
|-----------|------|--------|
| `cfwrite/` | the write-interop matrix: one file per case in `cases.json` (every BITPIX; int8/uint16/uint32/uint64 conventions; explicit BSCALE/BZERO with rounding up to the edges of the stored range incl. BITPIX 8 0..255 and lo-0.5+1 rounding to lo (`*_bzero_half`; the values beyond the edges are in `golden/write_overflow.json`); non-integer and negative BSCALE on uint16/uint32/uint64 (`u16_bscale_1p1`, `u16_bscale_neg2`, `u32_bscale_1p1`, `u32_bscale_neg3`, `u64_bscale_1p1`, `u64_bscale_neg1p5`) and caller scaling on int8 (`i8_bscale_1p1`, `i8_bscale_neg1`); BSCALE 1.1 values where dividing by BSCALE and multiplying by 1/BSCALE round differently (`*_bscale_1p1`); BITPIX 64 values where rounding to float32 directly and through float64 differ (`i64_f32_round`); BITPIX 64 with non-integer BSCALE incl. the int64 edges (`i64_bscale_half`, `i64_bscale_1p5`); negative BSCALE on BITPIX 16 and 32 (`i16_bscale_neg1`, `i16_bscale_neg2_bzero10`, `i32_bscale_neg3`); float32/float64 signaling and payload NaNs given by bit pattern (`"bits:0x..."` values, `f32_nanbits`, `f64_nanbits`; cfitsio writes them verbatim); a caller header with only BZERO or only BSCALE (`i16_bzero_only`, `i16_bscale_only`, `f64_bzero_only`; cfitsio writes both cards); BLANK; 1-D to 4-D) | written by cfitsio (`imgref write`) |
| `astropy_written/` | int8, uint8, uint16, uint32, uint64, int64, float32 NaN/Inf, a uint64 cube, an int16 BSCALE file from `ImageHDU.scale`, and `u16_stis_frame.fits`: the real 62x44 HST/STIS raw detector frame (uint16 counts) from `o4sp040b0_raw.fits[1]` re-written by astropy as a BZERO=32768 primary image | written by astropy 8.0.2 |
| `astropy_data/` | real files copied verbatim from astropy's `io/fits/tests/data` (BSD-3-Clause, (c) Astropy Developers): `arange.fits` (3-D int32), `blank.fits` (BITPIX 64 + BLANK), `checksum.fits` / `checksum_false.fits` (valid / wrong CHECKSUM, image + bintable), `memtest.fits` (CHECKSUM without DATASUM, bintable with wrong sums), `fixed-1890.fits` (uint16 100x100), `group.fits`, `random_groups.fits` (random groups), `o4sp040b0_raw.fits` (HST/STIS raw, real uint16 detector frames in 7 HDUs), `scale.fits` (non-integer BSCALE), `test0.fits` (multi-extension int16) | astropy 8.0.2 wheel, sha256 in `golden/manifest.json` |
| `framing/` | trailing zero block, trailing junk bytes, trailing junk block, last data unit without its final pad (alone and after a full HDU), a stray NAXIS3 card with NAXIS = 1 (once, and twice with different values), a malformed non-structural card (`DATE-OBS= 2020-01-01T00:00:00`, unquoted), a random-groups primary followed by an IMAGE extension, an `XTENSION = 'FOREIGN'` extension (NAXIS1 3000, PCOUNT 100) between the primary and an IMAGE (`xtension_foreign.fits`; cfitsio frames it and reads it as an IMAGE_HDU, the library refuses it with `ErrUnsupportedExtension`) | assembled byte by byte by `gen.py`; truth is cfitsio's reading (and astropy's, `golden/framing_astropy.json`, for the pixels cfitsio cannot read when the final pad is missing) |
| `refused/` | files the library refuses with a typed error: a data unit cut short (primary, and an extension after a full HDU), a primary without NAXIS2 / BITPIX / NAXIS, an extension without NAXIS1 / BITPIX, an extension whose XTENSION value is not a string | assembled byte by byte by `gen.py` (`gen_refused`); truth is cfitsio's status for each (`imgref open`, `golden/refused__<name>.open.json`): it refuses to open the primaries (222/223/224), will not move to the bad extensions (222/224), fails the read of the truncated data (107/108), and reads the XTENSION = 1 extension as an IMAGE (a documented difference) |
| `keywords/` | one int16 (or int32) image per BSCALE/BZERO/BLANK card form (`gen.py:KEYWORD_CASES`): D and d exponents, duplicates (also with a bad first or last card), quoted, logical, complex, Inf/NaN, over/underflowing, subnormal (exact and inexact) and hex values, `=` without the following blank, `BSCALE=2.0`, `BSCALE=1= 2` (value after the `= ` in columns 9-10), `0xDD` (only the first D replaced), no `=`, a `/ x=3` comment, an empty value, lowercase and HIERARCH names, BSCALE 0, values that round to or fall below the smallest normal, BLANK with `+`, a fraction, hex, overflow, outside int16, wrapping as int32, and cfitsio's NULL_UNDEFINED 1234554321. The stored values start 0, 1, so the BZERO/BSCALE cfitsio applies are its first two pixels | assembled byte by byte by `gen.py`; truth is cfitsio's reading (cfref dump, imgref pix); kept out of the corpus sample list because cfref's `fits_read_key` BSCALE differs from the scaling cfitsio applies on purpose |
| `duplicate/` | a second BITPIX, NAXIS, NAXIS1 (smaller and larger; also followed by a duplicated stray NAXIS3), SIMPLE in a primary; a second PCOUNT, GCOUNT, XTENSION in an IMAGE extension; a second TFIELDS in a BINTABLE. The library refuses them (`DuplicateKeywordError`); cfitsio opens them | assembled byte by byte; truth recorded: cfref dump |
| `checksum/` | DATASUM/CHECKSUM added by cfitsio to `cfwrite/u16_conv.fits`, `cfwrite/f64_plain.fits`, the STIS raw file, `testprog.std` (images + ASCII/binary tables) and a header-only HDU; `corrupt_data.fits` (one data byte flipped: both sums wrong) and `corrupt_header.fits` (one comment byte changed: only CHECKSUM wrong); `datasum_*.fits`: DATASUM values that cfitsio parses with `atof` and an `(unsigned long)` cast (leading blanks, a fraction, trailing junk, an exponent, hex, the sum + 2^32, unquoted; on a header-only HDU -0.5, -1 and -2^32, and a DATASUM / a CHECKSUM card with no value, where `fits_verify_chksum` fails with 204, recorded as `"verify_error"`; 2^32 on a zero sum); card order and name case from `f64_plain.fits` (`datasum_dup_*`: a second DATASUM after or before the CHECKSUM card, none after it, the one after it undefined, no CHECKSUM; `checksum_dup_empty_first`; lowercase `datasum` right and wrong: cfitsio's ffgcrd finds CHECKSUM from the top, DATASUM from the card after CHECKSUM, wrapping, names in upper case) | `imgref addchk`, then byte edits in `gen.py` (`gen_datasum`) |

## The scaling-card rule (`golden/keywords_rule.json`)

`gen.py gen_keyword_rule` records, for every `keywords/` sample, the
outcome of the library's one BSCALE/BZERO/BLANK rule (`KEYWORD_RULE`:
`accept`, or the reason for the error: `not_value_card`, `not_number`,
`cfitsio_range`, `hierarch`, `blank_range`, `blank_undefined`), astropy
8.0.2's scaled `getdata` of the sample (dtype and values as Python
`repr` strings, NaN as null; or the exception astropy raises) and whether
cfitsio's TDOUBLE reading (the `keywords__<name>.json` golden), cast to
astropy's dtype, equals it (`refs_agree`). The generator fails if an
accepted sample is read differently by the two references.
`TestScalingKeywordsRule` checks every accepted sample against both and
every other one for `ErrInvalidScaling` on every read path.

## What cfitsio does that the library deliberately does not

These are recorded so a reader of the goldens is not surprised; the tests
compare against the stored bytes where cfitsio substitutes values.

* Scaling cards cfitsio reads one way and astropy another (or refuses):
  cfitsio replaces an unconvertible BSCALE/BZERO/BLANK by the default,
  applies `HIERARCH BSCALE`, reads `BSCALE=2.0` and a strtod prefix,
  compares a BITPIX 32 BLANK as 32 bits and drops BLANK 1234554321; the
  library returns `ErrInvalidScaling` for all of them
  (`golden/keywords_rule.json`). The cfitsio goldens of those samples
  record cfitsio's reading.
* Null-checked float reads: cfitsio reports +-Inf as null and returns +0
  for -0 and subnormals (`fnan()` in `fitsio2.h`). The library returns the
  stored IEEE value and treats only NaN as undefined.
* TFLOAT reads of values beyond float32 range: cfitsio clamps to +-FLT_MAX
  with status 412 (`"f32_overflow": true` in the golden); Go's float32
  conversion gives +-Inf.
* Writing a value whose stored form leaves the BITPIX range: cfitsio
  stores the range limit and returns NUM_OVERFLOW (412); `WriteImage` and
  `OverwritePixels` return `ErrNumOverflow` and write nothing (until
  round 2 of the proof audit they clamped silently, the behaviour PR #2,
  feature/write-clamp, 0258acd, introduced). The refusal is exactly
  cfitsio's 412 per value (`golden/write_overflow.json`), on both edges,
  including a value of exactly lo-0.5 (`u8_bzero_half`, `i16_bzero_half`,
  `i32_bzero_half`), plus the values cfitsio stores wrongly without a
  status: see the next two items. The matrix (`cases.json`) keeps only
  storable values (cfitsio status 0 for every case), up to the edges.
* Float BITPIX with a scaling that sends a finite value to infinity
  (`f32_bscale_tiny`, `f64_bscale_tiny` in `write_overflow.json`): cfitsio
  stores the infinity with status 0; the library refuses it
  (`ErrNumOverflow`).
* Duplicate structural keywords (`duplicate/`): cfitsio reads BITPIX,
  NAXIS and NAXISn by card position and an image's PCOUNT/GCOUNT last;
  the library refuses the header (`DuplicateKeywordError`).
* Scaling cards cfitsio would read differently from the caller's card
  (duplicates, unconvertible values, BSCALE 0, BLANK 1234554321):
  `WriteImage` refuses them (`ErrInvalidScaling`); cfitsio's own writer
  refuses only BSCALE 0 (`golden/write_refused.json`).
* BITPIX 64, scaled value exactly 2^63 (`i64_bscale_half` in
  `write_overflow.json`: v = 2^62 and v = 2^62-1 with BSCALE 0.5, BZERO
  0.25). cfitsio's range test is `dvalue > DLONGLONG_MAX`, and
  DLONGLONG_MAX is 2^63 as a double, so 2^63 passes it and the
  out-of-range `(LONGLONG)` cast stores INT64_MIN (on x86-64) with no
  overflow status; likewise v = -2^63 with BZERO 1000 (`i64_bzero1000`),
  whose exact stored value -2^63-1000 cfitsio's double arithmetic rounds
  to -2^63. The library refuses these (`ErrNumOverflow`). The test asserts
  the golden really holds INT64_MIN there with status 0, and that it reads
  back further from v than one BSCALE step.
* `fits_get_img_equivtype` does not range-check BITPIX 64: it reports
  LONGLONG_IMG for any integral BZERO, including 2^63 (uint64) and 1000
  (whose physical range exceeds int64). The library's integer targets use
  the exact physical range instead.
* A trailing block that starts with `XTENSION=` but has no END card: cfitsio
  silently ends the HDU list there; the library returns
  ErrUnterminatedHeader (astropy also refuses such files).

## Storage of large goldens

`gen.py` (`put_golden`) writes a golden whose text is over 16384 bytes as
`<name>.gz` (Python `gzip.compress`, level 9, mtime 0); smaller ones stay
plain. The Go tests read both through `testutil.ReadGolden`, and
`make verify` compares the regenerated tree with the committed one through
`cmp_golden.py`, which decompresses `.gz` files first, so the harness
output is still compared byte for byte (the compressed bytes depend on the
zlib build and are not).
