# testdata/cfmovnam: cfitsio truth for HDUByName with duplicate EXTNAME cards

`dupext.fits` is written by cfitsio: `./mkdupext dupext.fits` (source
`mkdupext.c`). `fits_write_key` always appends a card, so HDUs 1, 3, 4 and 5
carry two EXTNAME cards each (HDU 1 `SCI` then `ERR`, HDU 2 `ERR` only, HDU 3
`DQ` then `SCI`, HDU 4 `TAB` twice, HDU 5 `LAST` then `DQ`; HDU 0 has none).

`golden/dupext.json` is exactly `./cfmovnam dupext.fits SCI ERR DQ TAB LAST
NONE` (source `cfmovnam.c`): for each name, the 0-based HDU that
`fits_movnam_hdu(ANY_HDU, name, 0)` moves to, or -1 for `BAD_HDU_NUM`.
cfitsio reads each HDU's EXTNAME with `fits_read_key_str`, which finds the
first EXTNAME card, so `SCI` is HDU 1, `ERR` HDU 2, `DQ` HDU 3 and `LAST` HDU 5.

`TestHDUByNameDuplicateEXTNAMEMatchesCfitsio` requires `HDUByName` to give
the same HDU (or an error for -1) for every name in the golden, for a
read-only file and in ModeEdit before and after `Header()` loads each header.
`TestCfmovnamReproducesGoldens` runs the built harnesses and requires
byte-identical output (tier 2, `go test -tags interop`: a harness that is
not built is a failure).

Those names avoid the cases where `HDUByName` deliberately differs from
`fits_movnam_hdu` (case-insensitive compare, the HDUNAME fallback, a
non-string EXTNAME); see the `HDUByName` doc. Those cases are in
`names.fits`, written by cfitsio with `./mknames names.fits` (source
`mknames.c`; every HDU is a NAXIS = 0 image, one block each): HDU 1 `EXTNAME = 'sci'`, HDU 2 `'SCI'`, HDU 3 only
`HDUNAME = 'ERR'`, HDU 4 `EXTNAME = 'ERR'`, HDU 5 `EXTNAME = 5` (an integer,
`fits_write_key` TINT) then `EXTNAME = 'NUM'`, HDU 6 `EXTNAME = 'NUM'`.
`golden/names.json` is exactly `./cfmovnam names.fits SCI sci Sci ERR NUM 5
NONE`. `TestHDUByNameDifferencesFromCfitsio` pins `HDUByName`'s answer for
each name together with whether cfitsio agrees, and requires the golden to
say the same: `SCI` (cfitsio 1, ours 2), `Sci` (cfitsio 1, ours none),
`ERR` (cfitsio 3 by HDUNAME, ours 4) and `5` (cfitsio 5, ours none)
differ; `sci`, `NUM` (HDU 5's first EXTNAME card is 5, so both give HDU 6)
and `NONE` agree.

## Tool versions

| tool       | version |
|------------|---------|
| libcfitsio | 4.7.0 (Arch `cfitsio 1:4.7.0-1`, system `/usr/lib/libcfitsio.so.4.7.0`) |
| compiler   | gcc 16.2.1, `-O2 -Wall -Wextra -std=c99` |

## Commands

```sh
make -C testdata/cfmovnam           # build ./cfmovnam, ./mkdupext and ./mknames (gitignored)
make -C testdata/cfmovnam goldens   # rewrite dupext.fits, names.fits and golden/*.json
make -C testdata/cfmovnam verify    # regenerate and cmp against the committed files
```
