# header/testdata: samples and cfitsio goldens for the header package

The header package's decoder and encoder are tested against cfitsio
reading the same bytes. Nothing here is produced by the Go code.

## Samples

| file             | made by                                  | what it holds |
|------------------|------------------------------------------|---------------|
| `malformed.fits` | `gen_samples.py` (plain Python bytes)    | value cards that are not valid FITS: unquoted string, `TRUE`, `Tx`, integers outside int64 (incl. 2^63 and 2^64-1), unterminated quote, float overflow / underflow / subnormal, `NaN`, `Inf`, hexadecimal, lower-case `d` exponent, two tokens, complex without `)`, `=` without the blank / in column 10, undefined values, comment spacing, string blanks and quotes, lower-case keyword, HIERARCH without `=` / with runs of blanks, `HIERARCH=` |
| `continue.fits`  | `gen_samples.py`                         | CONTINUE chains with comments, lone CONTINUE (no `&`, after an integer), dangling `&`, empty CONTINUE with and without comment, `CONTINUE=` with a value indicator, unquoted continuation, HIERARCH long string, full-width first card (comment joined without a blank), a quote split across cards |
| `adversarial.fits` | `gen_samples.py`                      | about 80 more edge cases, most first written by an adversarial review: keyword names cut by a blank or `=` inside columns 1-8 or longer than 8 characters (cfitsio `ffgknm`), `-0`, number forms outside FITS (`1E`, `++5`, `1_000`, `0x1p3`, `1.5E3D`, subnormal and boundary doubles, a 22-digit integer), complex forms (`( 1 , 2 )`, `(1.5D2, 3)`, `(1,2,3)`, `(1, )`, `(,1)`, `( )`, `(1 2)`), `HISTORY=` / `COMMENT =` with a value, a lone CONTINUE chain with an empty comment-only piece, and cards that reach cfitsio's buffer limits (value token over 70 characters, quoted string over 70, complex over 70, comment over 72) |
| `cfwritten.fits` | `hdrref write` (cfitsio's own writer)    | `fits_write_key_longstr` (200-char string with quotes and a 120-char comment, 68/69-char boundaries, leading blanks), HIERARCH long and short strings, `fits_write_history` / `fits_write_comment` of 200 / 150 characters, doubles (120.0, 1e-7, 1.5e300, 0.1), uint64 max, int64 min, logical, complex |

Each is one NAXIS = 0 primary HDU (2880-5760 bytes).

## Goldens (`golden/`)

| file                                  | command |
|---------------------------------------|---------|
| `<sample>.json`                       | `testdata/cfref/cfref dump <sample>` (the shared cfitsio harness: every card, cfitsio's value/comment split and typed reads) |
| `<sample>.longstr.json`               | `hdrref longstr <sample>`: `fits_read_key_longstr` value and comment, `TULONGLONG` and `TDBLCOMPLEX` reads |
| `testdata/<repo sample>.longstr.json` | `hdrref longstr` of every FITS file in the repository's `testdata/` |
| `testkey.json`                        | `hdrref testkey ...`: `fits_test_keyword` status for a list of names |

`hdrref readkey <file> <name>...` and `hdrref delkey <in> <out> <name>` are
used live by `TestLookupCaseAndFirstCardMatchCfitsio`,
`TestDeleteFirstCardMatchesCfitsio`,
`TestDuplicateKeywordFirstCardMatchesCfitsio` and
`TestDeleteCommentaryRemovesFirstMatchesCfitsio` (no committed golden):
cfitsio's `fits_read_keyword` and `fits_delete_key` from the top of a
header (one this package encoded, or one written byte for byte by the
test), with lower-case names, duplicated keywords and repeated HISTORY
cards. `hdrref writekeys <out>` is used live by
`TestFixedFormatMatchesCfitsioBytes`: cfitsio's `fits_write_key` writes
integers of every printed width (1 to 20 characters), uint64 max and both
logicals, and the test compares the card bytes with this package's; it
also passes the lower-case name `lowkey`, which cfitsio writes as
`LOWKEY` (`TestLowerCaseKeyWrittenUpperCasedLikeCfitsio`).

## Tool versions

| tool       | version |
|------------|---------|
| libcfitsio | 4.7.0 (Arch `cfitsio 1:4.7.0-1`); every longstr golden records it and the tests require `4.7.0` |
| compiler   | gcc 16.2.1, `-O2 -Wall -Wextra -std=c99` |
| python     | 3.11.11 (stdlib only for `gen_samples.py`) |
| fitsverify | 4.22 (CFITSIO V4.7.0), used at test time on encoder output |
| astropy    | 8.0.2 with numpy 2.4.6 (repo venv), used at test time as the second reader of encoder output |

## Commands

```sh
make -C header/testdata hdrref/hdrref   # build the harness (gitignored)
make -C header/testdata samples         # rewrite the four samples
make -C header/testdata goldens         # rewrite every golden
make -C header/testdata verify PYTHON=venv/bin/python
                                        # regenerate samples + goldens into
                                        # verify.out/ and compare byte for byte
```

`make clean` removes only the binary and `verify.out/`.

## cfitsio behaviour recorded here (and how the Go package treats it)

- cfitsio opens every sample; per-card problems surface only when that
  card is read. The Go decoder likewise never fails a header for a bad
  value: such a card is `TypeUnparsed` with cfitsio's value text.
- cfitsio reads `TRUE` and `Tx` as logical true, and a hexadecimal token
  (`0x10`) as the number 16 (C `strtod`). FITS has neither; the Go decoder
  keeps them as `TypeUnparsed`.
- cfitsio reports float underflow (`1.0E-400`) and subnormals (`4.9E-324`)
  as conversion overflow (status 412); so does the Go decoder.
- `fits_read_key_longstr` truncates a joined CONTINUE comment to 72
  characters; the Go parser keeps all of it (cfitsio's is its prefix).
- cfitsio does not count blank cards just before END as keywords; the Go
  parser keeps them as blank commentary cards (needed for byte-identical
  re-encoding).
- cfitsio names a non-HIERARCH card from column 1 up to the first blank
  or `=` (`ffgknm`, at most 74 characters): `B9  =    5` is `B9`,
  `KEYLONG10= 5` is `KEYLONG10`. The Go decoder does the same.
- For a HIERARCH card with runs of blanks in its path cfitsio returns the
  runs (`ESO   DBL   SPACE`); the Go Key has single blanks
  (`HIERARCH ESO DBL SPACE`). The test asserts exactly that form.
- cfitsio reads a complex value with an empty part (`(1, )`, `(,1)`) as
  having 0 there; FITS has no such value and the Go decoder keeps it as
  `TypeUnparsed`. Where cfitsio's complex read fails (`(1,2,3)`, `( )`,
  `(1 2)`) the Go card is `TypeUnparsed` too.
- cfitsio's `fits_read_key_dbl` reads the integer card `-0` as -0; so does
  `Header.Float`.
- The commentary text of a card without a value (columns 9-80) is at
  most 72 characters, so the 72-character comment limit modelled on
  cfitsio is reachable only through the value-card comment
  (`L5`/`L6` in `adversarial.fits`), which pins it.
- A string written by the Go encoder that needs CONTINUE cards loses
  its trailing blanks first. Otherwise a last piece of only blanks
  without a comment, `CONTINUE  '     '`, is read by cfitsio and by this
  package as an empty CONTINUE: the chain ends and the `&` stays in the
  value (astropy drops it, so the two readers disagree). Such a string
  ending in `&` is refused for the same reason: cfitsio keeps the `&`,
  astropy drops it.
- Joined CONTINUE comments. The Go decoder joins comment pieces as
  cfitsio's `fits_read_key_longstr` does: one blank between pieces,
  none after a card that is a full 80 columns wide. cfitsio then cuts
  the joined comment at 72 characters; Go keeps all of it. Past 72
  characters the truth is the writer's input: `cfwritten.fits` LONGSTR
  was written by `hdrref.c` with `c120 = rep("comment text ", 120)`, and
  `TestJoinedContinueCommentIsWholeComment` requires exactly that text.
  astropy 8.0.2 joins differently and is not used as the truth here: it
  reads that comment as `... com ment text ...` (a blank inserted where
  cfitsio split a word across cards) and the `FULL` comment of
  `continue.fits` as `c...c more` where cfitsio and Go read `c...cmore`.
- cfitsio 4.7.0 aborts (glibc fortify, buffer overflow in `ffc2ujj`) when
  `fits_read_key(TULONGLONG)` meets an integer outside uint64, e.g. the
  24-digit `BIGINT`; `hdrref` therefore reads `TULONGLONG` only for
  unsigned decimals in range and records `"skipped"` otherwise.
