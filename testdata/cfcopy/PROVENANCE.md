# testdata/cfcopy: cfitsio truth for EditFile HDU copy / reorder / skip

`cfcopy.c` copies a chosen list of HDUs from a committed sample into a new
file with cfitsio's `fits_copy_hdu`, which converts a primary array that
lands after position 0 into an IMAGE extension, an IMAGE extension that lands
at position 0 into a primary array, and writes an empty primary in front of
a table copied to position 0 (cfitsio `edithdu.c`, `ffcphd`).

`golden/<scenario>.fits` is exactly `./cfcopy golden/<scenario>.fits
../<sample> <hdu> ...` for each scenario in the `Makefile`. The Go tests in
`editfile_copy_test.go` require EditFile's output for the same HDU order to
match these files card for card and data byte for data byte, cfitsio
(`testdata/cfref`) to read both identically, and fitsverify to report the
same counts for both.

## reserve.fits

`reserve.fits` is a scenario input, also written by cfitsio: `./cfcopy -r 40
reserve.fits ../testprog.std 0 2 6` copies HDUs 0, 2 and 6 of testprog.std
with `fits_copy_hdu(..., morekeys = 40, ...)`, so each header ends with
reserved blank cards before END (30, 29 and 21 of them). The
`reserve_1_0_2` scenario moves its IMAGE extension to the primary position
and its primary to an extension, which shows that a converted header drops
the reserved cards (cfitsio's keyword count stops at the last non-blank
card) while an HDU that keeps its form keeps them.

## comments.fits

`comments.fits` is the input of the `comments_1_0` scenario, written by
cfitsio: `./mkcomments comments.fits` (source `mkcomments.c`). Its primary
array carries, after the standard cards `fits_create_img` writes, COMMENT
cards that share the first 58 (or 47) characters with the two standard
"FITS (Flexible Image Transport System) format ..." comments but differ
later, cards that differ inside those prefixes, and a keyword `EXTENDED`.
When `fits_copy_hdu` turns that primary into an IMAGE extension it drops
EXTEND and every card matching the prefixes (cfitsio `edithdu.c`, `ffcphd`,
`FSTRNCMP(card, ..., 58)` and `47`) and keeps the others; the golden
shows which. Its IMAGE extension carries the two standard comments, which
the reverse conversion keeps after the ones it inserts. fitsverify reports
no warnings or errors for it.

## Tool versions

| tool       | version |
|------------|---------|
| libcfitsio | 4.7.0 (Arch `cfitsio 1:4.7.0-1`, system `/usr/lib/libcfitsio.so.4.7.0`) |
| compiler   | gcc 16.2.1, `-O2 -Wall -Wextra -std=c99` |
| fitsverify | system `/usr/bin/fitsverify` (used at test time, not for the goldens) |

## Inputs

The scenarios read the committed samples `testdata/testprog.std`,
`testdata/iter_image.fit`, `testdata/vari.fits` (provenance in
`testdata/README.md`), `testdata/comp_fz_u16_rice.fits` (fpack-written
tile-compressed uint16, `testdata/COMP_FZ_PROVENANCE.md`; scenarios
`comp_skip0` and `comp_1_0`, whose copied pixels are also checked against
the source's cfref golden), `reserve.fits` and `comments.fits` (above),
and the funpack-written data cubes `cube3.fits` (NAXIS 3), `cube4.fits`
(NAXIS 4) and `cube6.fits` (NAXIS 6, see below) (scenarios `cube3_0_0`,
`cube4_0_0`, `cube6_0_0`: a cube primary to positions 0 and 1).

Two scenarios start from an earlier golden, so `make goldens` writes them
in Makefile order: `cube3_1_0` moves the NAXIS 3 IMAGE extension of
`golden/cube3_0_0.fits` back to the primary position, and `naxis0_1_0`
does the same for the NAXIS 0 IMAGE extension of `golden/naxis0_0_0.fits`
(vari.fits's NAXIS 0 primary copied to positions 0 and 1).

## Commands

```sh
make -C testdata/cfcopy           # build ./cfcopy and ./mkcomments (gitignored)
make -C testdata/cfcopy goldens   # rewrite reserve.fits, comments.fits and golden/*.fits
make -C testdata/cfcopy verify    # regenerate and cmp against the goldens
```

`TestCfcopyReproducesGoldens` runs the built harnesses (for `reserve.fits`,
`comments.fits` and every golden) and requires byte-identical output (tier 2,
`go test -tags interop`: a harness that is not built is a failure).

## cube3.fits, cube4.fits, cube6.fits

Inputs of the `cube3_0_0`, `cube4_0_0` and `cube6_0_0` scenarios (NAXIS 3, 4
and 6 int16 cubes). Each is `funpack -C` (cfitsio 4.7.0, `/usr/bin/funpack`)
of a committed compressed cube sample: `../comp_fz_cube16_rice_t3.fits`,
`../comp_fz_cube4d16_rice.fits`, `../comp_fz_cube6d16_gzip2.fits`. They are
byte-identical to the funpack twins that were committed as
`../comp_fzun_cube*.fits` before those twins became generated files
(`testdata/gen_twins.sh`). `make goldens` writes them and `make verify`
regenerates and compares them.
