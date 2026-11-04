# testdata

FITS samples and the truth (goldens) the root `fits` package is tested
against. Every claim the tests make about a sample is checked against an
independent tool (libcfitsio 4.7.0, fpack/funpack, fitsverify, astropy
8.0.2, wcslib); see `../docs/VALIDATION.md` for which test proves what.

Do NOT hand-edit samples or goldens. Each one is produced by a committed
generator, and regenerating must reproduce the committed bytes. Large
text goldens are committed gzipped as `<name>.gz` (level 9, no name or
time stamp); the tests read them through `testutil.ReadGolden` and the
reproduce checks compare the decompressed bytes. `TestRepositorySizeBudget`
keeps every tracked file under 1 MB and all testdata under 20 MB.

## Top-level samples

| files | origin | regenerate |
|-------|--------|------------|
| `testprog.std`, `testf77.std`, `iter_a.fit`, `iter_b.fit`, `iter_c.fit`, `iter_image.fit`, `vari.fits` | copied verbatim from the HEASARC cfitsio source tree (`../reference/cfitsio/`). Despite its name, `vari.fits` holds a plain binary table (7 scalar `1E` columns, 91 rows, PCOUNT 0): it has no variable-length arrays. VLA samples live in `tables/`. | not generated |
| `comp_<name>.fits`, `comp_raw_<name>.fits` (36) | astropy 8.0.2 / numpy 2.4.6: each compressed fixture and its uncompressed twin | `venv/bin/python testdata/generate_compressed.py` (reproducible: no gzip mtime) |
| `comp_fz_*.fits`, `comp_fzsrc_*.fits` (148 cases plus 12 ZBITPIX 64 cases), `comp_bad/*.fits` | fpack/funpack 1.7.0 (CFITSIO 4.7.0) and libcfitsio 4.7.0 (`gen_fpack.sh` + `gen_fpack.c`); `fzsrc` is the source image. The funpack twins `comp_fzun_*.fits` are not committed: `testdata/gen_twins.sh` regenerates them into the gitignored directory twins.out under testdata/cfref, `comp_fzun_same.txt` maps the cases whose twin is another case's, and their cfitsio goldens are `cfref/golden/twins/` | `FITS_ASTROPY_PYTHON=$PWD/venv/bin/python sh testdata/gen_fpack.sh`, then `FITS_ASTROPY_PYTHON=... make -C testdata/cfref goldens` (also writes `comp_bad/golden`) |
| `comp_fz_real_with_nan.fits`, `comp_fz_real_image.fits` | verbatim astropy 8.0.2 test data (BSD-3) | not generated |
| `wcs_crosscheck.json` | astropy 8.0.2 (wcslib 8.6), numpy 2.4.6; format `{generator, astropy, wcslib, numpy, fixtures:[{name, source, cards, tol_arcsec, tol_px, p2s, s2p}]}`. The fixtures are synthetic or transcribed headers; the real archive headers live in `wcs/transform/testdata/wref/hdrcases/real_*`. Checked by `TestCrossCheckFixtures`, reproduced by `TestCrossCheckFixturesReproduce` | `venv/bin/python testdata/generate_wcs_crosscheck.py > testdata/wcs_crosscheck.json` |

Provenance and corpus rules for the compressed samples: `COMP_FZ_PROVENANCE.md`.

## Subdirectories

| dir | what | provenance |
|-----|------|------------|
| `cfref/` | cfitsio truth harness (`cfref.c`) and `golden/<sample>.json` for every top-level sample: HDU addresses, every card, every pixel and table cell as libcfitsio reads it. Arrays over 1024 values are sha256 digests with samples; mid-size arrays (`<kind>-<sha256>.json`) and header cards (`card-<sha256>.json`) are stored once in `golden/pool`, committed gzipped (`.json.gz`) | `cfref/PROVENANCE.md` |
| `cfcopy/` | cfitsio `fits_copy_hdu` goldens for `EditFile` HDU conversion, plus the cfitsio-made inputs `reserve.fits` and `comments.fits` and the funpack-made cube inputs `cube3.fits`, `cube4.fits`, `cube6.fits` | `cfcopy/PROVENANCE.md` |
| `cfmovnam/` | cfitsio `fits_movnam_hdu` goldens for `HDUByName`: `dupext.fits`, `names.fits` (made by `mkdupext.c`, `mknames.c`), `golden/dupext.json`, `golden/names.json` (`cfmovnam.c`); `make -C testdata/cfmovnam verify` | `cfmovnam/PROVENANCE.md` |
| `images/` | image read/write truth corpus: cfitsio-written write matrix (`cfwrite/`), astropy-written and real samples, framing edge cases, `imgref` goldens (pixels, subsets, checksums) | `images/PROVENANCE.md` |
| `imaging/` | `ReadFloat32` / `ReadFrame` / `AsImage` sample matrix incl. BITPIX 32/64 rounding samples, a real ASI533MC Pro crop, compressed scaling/ZBLANK variants, the Siril harness (`sirilref/`) and the cfitsio null-flag harness (`cfnull/`) | `imaging/PROVENANCE.md` |
| `tables/` | binary and ASCII table samples (cfitsio `gen_tables.c`, astropy, real Gaia; incl. `cf_ascii_tzero_round.fits`, `cf_bintable_tzero_round.fits` and the heap layouts `heap/heap_neg_off*.fits`, `heap/heap_theap_float.fits`), cfitsio table goldens, and `golden_tfloat/` (cfitsio TFLOAT reads, `tables/cf_tfloat.c`) | `tables/PROVENANCE.md` |

Other packages keep their own samples next to their code: `header/testdata`,
`compress/testdata/cref`, `wcs/transform/testdata/wref`, `wcs/transform/testdata/frames`, `wcs/transform/testdata/tnx`,
`healpix/testdata`, `stats/testdata`, `stretch/testdata`.

## How the samples are used

- `TestCopyToByteIdentical` (`fits_fixtures_test.go`): `CopyTo` reproduces
  every top-level sample byte for byte. This checks framing only.
- `TestHDULayoutMatchesCfitsio`: HDU offsets, kinds and header card text
  equal cfitsio's for every sample.
- `TestGoldensCoverEveryTestdataSample`: every sample has a cfitsio golden
  and every golden a sample.
