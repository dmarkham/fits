# testdata

FITS fixtures used by the fits library test suite.

## Provenance

All `.fit`, `.fits`, and `.std` files in this directory are copied verbatim
from the HEASARC cfitsio reference source tree at
`../reference/cfitsio/`.

| file             | origin                | purpose                                                |
|------------------|-----------------------|--------------------------------------------------------|
| iter_image.fit   | cfitsio               | image HDU with multiple extensions                     |
| iter_a.fit       | cfitsio               | iterator example input (binary table)                  |
| iter_b.fit       | cfitsio               | iterator example input                                 |
| iter_c.fit       | cfitsio               | iterator example input                                 |
| vari.fits        | cfitsio               | variable-length-array binary table (critical for VLA)  |
| testprog.std     | cfitsio testprog.c    | end-to-end correctness oracle (~69 KB, many features)  |
| testf77.std      | cfitsio testf77.f     | Fortran-origin oracle (secondary)                      |

## Usage

Regression tests in `../fits_fixtures_test.go` open each fixture through
`fits.Open`, iterate every HDU, and verify that `(*File).CopyTo(w)` emits
byte-for-byte identical output (plan §"Done criteria" item 5).

Do NOT modify these files. If a test needs a variant (different BITPIX,
different keyword layout, etc.) generate a synthetic fixture under a new
name.
