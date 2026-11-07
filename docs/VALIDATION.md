# Validation ledger

This file lists every capability the library claims, the committed test that
proves it, the independent truth source (with its version) the test compares
against, the committed samples and goldens it uses, and the tolerance with
the reason for it. It ends with what is **not** supported or not yet
proven, and with how to regenerate every golden.

Rules every entry follows:

1. The expected values come from an external tool (or exact arithmetic),
   never from this library's own Go code.
2. Goldens and samples are committed next to the generator that made them,
   and re-running the generator reproduces them byte for byte. The
   `*Reproduce*` tests enforce it for: cfref (`TestCfrefDumpReproducesGoldens`),
   tables (samples and goldens: `TestTableSamplesReproduce`,
   `TestTableGoldensReproduce`), cfcopy, header (samples and goldens:
   `TestHeaderSamplesReproduce`), compress/cref
   (`TestCrefGoldensReproduce`), testdata/images (`TestImagesCorpusReproduce`),
   imaging samples / cfnull / sirilref / gen_expected
   (`TestImagingSamplesReproduce`, `TestCfnullGoldensReproduce`,
   `TestSirilrefGoldensReproduce`, `TestImagingExpectedReproduce`), healpix
   goldens and map samples (`TestHealpixGoldensReproduce`,
   `TestHealpixMapSamplesReproduce`), stretch sref and cv2
   (`TestSirilGoldensReproduce`, `TestCLAHECv2GoldensReproduce`), stats
   (`TestGoldensReproduce`) and every WCS set (wref, hdr, the generated
   case files `TestPrjCaseGeneratorReproduce` /
   `TestHighOrderCaseGeneratorReproduce`, frames, TNX goldens and synthetic
   headers `TestTNXSyntheticHeadersReproduce`, crosscheck), and the
   compressed read corpus `testdata/gen_fpack.sh` +
   `testdata/generate_compressed.py` (`TestCompressedCorpusReproduce`).
   Checked by hand only: the sources that need downloads (see the end of
   this file).
3. A missing committed golden or sample is a test failure, never a skip.
   Plain `go test ./...` (tier 1) runs no external tool at all; tests that
   run one live are tier 2 (`-tags interop`), where a missing tool is a
   failure that names it (see Test tiers).
4. Tolerances are as tight as the truth allows: exact for integers and
   codecs, a stated ULP / arcsec / pixel budget with its reason for floats.

`TestValidationLedger` (`validation_test.go`) parses this file and fails if
a test named here (any backticked `Test...` / `Fuzz...` name) has no
`func` in the repository, if a backticked repository path (or glob) under
`testdata/`, `header/`, `compress/`, `wcs/`, `healpix/`, `stats/`,
`stretch/`, `cmd/`, `internal/` or `docs/` matches no file, or if a ledger
row names no test.

## Test tiers

**Tier 1: `go test ./...`** needs nothing but Go and the files committed
in the repository: no C compiler, no make, no cfitsio / fpack / funpack /
fitsverify / wcslib / WCSTools program, no Python, no environment
variable, no network. Every check compares the library with committed
goldens, committed samples, or **recorded tool readings**. No Go file
that tier 1 compiles starts another program, except the test binary
itself (cross-process tests): `TestTier1StartsNoPrograms` scans every Go
file a plain build compiles on linux, darwin and windows (amd64, arm64)
and fails on `os/exec` (other than `exec.Command(os.Args[0], ...)`),
`os.StartProcess` or `syscall.Exec`; `TestTier1ScannerCatches` pins the
scanner. Its skips are platform ones only (editing on Windows, a file
system without POSIX ACLs, a user without a supplementary group).

**Tier 2: `go test -tags interop ./...`** compiles the
`*_interop_test.go` files and `internal/testutil/interop.go` and `internal/testutil/record_live.go`.
It runs the reference programs live: the cfitsio, wcslib and healpix
harnesses (built first with `make -C <dir>`), fitsverify, fpack, funpack,
`xy2sky`/`sky2xy`, `HPXcvt`, the Python in `FITS_ASTROPY_PYTHON` (astropy,
numpy, healpy, mpmath, cv2), and for the regeneration of the Siril and
HEALPix goldens the source trees in `FITS_SIRIL_SRC` and
`FITS_HEALPIX_CXX` and the executable `FITS_SIRIL_CLI`. Opting in means
every tool is required: a missing tool, harness or variable is a test
failure that names it, never a skip. The `*Reproduce*` tests (rule 2)
are tier 2; byte-exact regeneration is pinned to amd64
(`testutil.RequireGoldenPlatform`, a documented platform skip on arm64).

**The write-golden chain.** A test that judges a file the library wrote
(or a file a tool wrote) by an external reader calls `testutil.Run`
(`internal/testutil/record.go`; the helpers `CfrefDump`, `Fitsverify`,
`AstropyScript`, `AstropyCode` and the per-harness tools go through it).
A reading is identified by a key over the tool, the committed source of
the harness (e.g. `testdata/cfref/cfref.c`), the arguments and standard
input with file paths replaced by placeholders, and the base name and
sha256 of every input file. The readings are committed in
`testdata/recorded/<package>/<TestName>.jsonl.gz` (one gzipped bundle per
top-level test, every distinct output line stored once):

1. Tier 2 with `FITS_RECORD=1` runs the tools on the bytes the library
   writes, hands their live readings to the test's assertions, and records
   the readings of every test that **passes**: the tools verify the bytes
   first, the recording follows.
2. Tier 1 writes the same file, computes its sha256 and replays the
   reading recorded for exactly those bytes; the test's own assertions run
   on it, so every proof (cfitsio reads every pixel and card as intended,
   fitsverify finds no error, astropy agrees) still holds. A writer that
   produces other bytes finds no recording and fails, naming the file,
   its sha256 and the command that re-verifies and records it. So tier 1
   also asserts that the writer reproduces the verified bytes exactly.
3. Tier 2 without `FITS_RECORD` re-runs every tool live and, on amd64,
   requires each live reading to equal the committed one, so the
   recordings cannot drift from the tools.

A reading too large to commit value by value (the decoded pixels of the
compressed-writer matrix, `TestCompressedWriterInterop`) is recorded as an
exact sha256 digest plus a few sampled values (`testutil.FloatsDigest`),
and only when the live reading equals a candidate the test computes (this
package's own decode, or another reader's reading); tier 1 accepts the
candidate only if it matches the digest exactly. A committed sample that
has a cfitsio golden (`testdata/cfref/golden`) is read from its golden.
The output of a tool that is only handed to another tool (a funpacked
image read back by the cfitsio harness) is recorded by digest only.

Each ledger row's **Tier** column says where its proof runs: `1` (tier 1,
committed goldens or samples), `1 + 2 live` (tier 1 replays the recorded
readings of the exact written bytes; tier 2 re-runs the tools live), `2`
or `2: <tests>` (live tools or golden regeneration, tier 2 only).

Run everything:

```sh
go vet ./... && go vet -tags interop ./...
go test ./...                                   # tier 1: Go only
for d in $(git ls-files '*Makefile' | xargs -n1 dirname); do make -C $d; done   # harnesses (hpref is built by its test)
FITS_ASTROPY_PYTHON=$PWD/venv/bin/python \
FITS_SIRIL_SRC=<siril-1.4.0 source> FITS_HEALPIX_CXX=<healpix_cxx 8.3.0 source> \
FITS_SIRIL_CLI=<siril-cli 1.4.0> go test -tags interop ./...   # tier 2
FITS_RECORD=1 <the same> go test -tags interop -run '^TestName$' ./pkg   # re-record one test's readings
GOARCH=arm64 go test -exec qemu-aarch64-static ./...   # tier 1 on arm64 (fused multiply-adds)
```

The tier 2 tests `TestFMAGuards_GOAMD64v3` and `TestCrossArch_arm64`
re-run the `stretch` package's tier 1 tests with GOAMD64=v3 and built for
arm64 under qemu.

## Truth sources

| Tool | Version | Used for |
|------|---------|----------|
| libcfitsio | 4.7.0 (Arch `cfitsio 1:4.7.0-1`) | HDU framing, cards, pixels, tables, checksums, compression read/write, `fits_copy_hdu` |
| fpack / funpack / fitsverify | 1.7.0 / 1.7.0 (both linked to CFITSIO 4.7.0) / 4.22 | compressed corpus, decode twins, validity of every written file |
| astropy / numpy | 8.0.2 / 2.4.6 (Python 3.11.11 venv) | read-back of written files, WCS crosscheck, frames (ERFA), sigma clip, histograms, sample generation |
| wcslib | 8.10 (`hdr2sky.c`, `wref.c`) | every projection and the full header pipeline |
| WCSTools | 3.9.7 (`/usr/bin/xy2sky`, `sky2xy`) | IRAF TNX |
| healpy / healpix_cxx | 1.20.1 / 8.3.0 (`hpref`) | HEALPix |
| Siril | 1.4.0 release source, compiled verbatim (`stats/testdata/sref`, `stretch/testdata/sref`, `testdata/imaging/sirilref`) | stats, stretches, ReadFloat32 convention |
| OpenCV | 5.0.0 (opencv-python-headless 5.0.0.93, Arch opencv 5.0.0) | CLAHE and Lab colour conversion |
| mpmath | 1.4.1 | libm accuracy reference values, XPH and HPX paper equations, the correctly rounded cos/sin/atan2 fed to healpix_cxx |
| HPXcvt | wcslib 8.10 | real HPX/XPH image headers |
| gcc / glibc | 16.2.1 / 2.44, x86-64 | all C harnesses |

## Ledger

### File framing and HDUs (root package)

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| HDU header/data/end addresses, count and kind equal cfitsio's for every sample | `TestHDUAddressesMatchCfitsio`, `TestHDULayoutMatchesCfitsio`, `TestCfrefGoldenHDUAddresses` | cfitsio `fits_get_hduaddrll` / `fits_get_num_hdus` | `testdata/cfref/golden`, `testdata/images/golden`, `testdata/images/framing` | exact (integers) | 1 |
| An extension type the library does not decode (XTENSION = 'FOREIGN') is framed like cfitsio (the IMAGE after it is read) and refused with `*UnsupportedExtensionError` / `ErrUnsupportedExtension` naming the XTENSION value (cfitsio reads it as an IMAGE_HDU); fitsdump reports it in place | `TestUnsupportedExtensionFramedAndRefused`, `TestHDUAddressesMatchCfitsio`, `TestFramingToleranceMatchesCfitsio`, `TestDumpReportsUnsupportedExtension` | cfitsio `fits_get_hduaddrll` and pixels of the IMAGE after it | `testdata/images/framing/xtension_foreign.fits`, `testdata/images/golden/framing__xtension_foreign.fits.json` | exact | 1 |
| `HDUByName`: the first HDU in file order, the primary included, whose EXTNAME equals the name case-sensitively (cfitsio's `fits_movnam_hdu` also starts at the primary but compares case-insensitively); the EXTNAME compared is the one `Header()` holds now, so an unflushed rename counts (as for cfitsio, whose header holds the update) and a flushed one is followed | `TestHDUByNameFirstMatchCaseSensitive`, `TestHDUByNameSeesUnflushedRename`, `TestFlushedFileCopiesAndIndexesNewHeaders` | the documented contract; cfitsio ffmnhd loop order | — | exact | 1 + 2 live |
| `HDUByName` with two EXTNAME cards in one header uses the first, as `fits_movnam_hdu` (`fits_read_key_str`) does, read-only and in ModeEdit before and after `Header()` loads the header | `TestHDUByNameDuplicateEXTNAMEMatchesCfitsio`, `TestCfmovnamReproducesGoldens` | cfitsio 4.7.0 `fits_movnam_hdu` (cfmovnam) on a cfitsio-written file | `testdata/cfmovnam/dupext.fits`, `testdata/cfmovnam/golden/dupext.json` | exact (HDU index) | 1; 2: `TestCfmovnamReproducesGoldens` |
| `HDUByName` where it deliberately differs from `fits_movnam_hdu` (case-sensitive compare, no HDUNAME fallback, a first EXTNAME card with no string value is no name): each name's HDU is pinned with whether cfitsio agrees (`SCI` HDU 2 vs cfitsio 1, `Sci` none vs 1, `ERR` 4 vs 3 by HDUNAME, `5` none vs 5; `sci`, `NUM`, `NONE` agree), read-only and in ModeEdit | `TestHDUByNameDifferencesFromCfitsio`, `TestCfmovnamReproducesGoldens` | cfitsio 4.7.0 `fits_movnam_hdu` (cfmovnam) on a cfitsio-written file | `testdata/cfmovnam/names.fits`, `testdata/cfmovnam/golden/names.json` | exact (HDU index) | 1; 2: `TestCfmovnamReproducesGoldens` |
| A card named END (cfitsio's ffgknm: up to the first blank, `=` or NUL; for a HIERARCH card the text between `HIERARCH ` and `=`, blanks trimmed) that is not the END card (`END  = 1`, `END=1`, `END     junk`, NUL-padded, `HIERARCH END = 1`, `HIERARCH  END= 'x'`) is refused by `Open` with `ErrMalformedEND`; cfitsio also refuses it where a required keyword must be, and elsewhere ends the header there (with the END card in the next block it reads the data at 2880); `HIERARCH ENDX`, `HIERARCH ESO END`, `HIERARCH end` stay ordinary cards for both | `TestMalformedENDRefused`, `TestOpenedHeadersAlwaysParse` | cfitsio refusal and data offset (live cfref dump); FITS 4.0 §4.4.1.1 | — | exact | 1 + 2 live |
| Random groups sized per FITS eq. 4, then refused with `ErrRandomGroups` | `TestRandomGroupsAreRefusedNotMisread` | cfitsio addresses; astropy's real random-groups files | `testdata/images/astropy_data/random_groups.fits`, `testdata/images/astropy_data/group.fits`, `testdata/images/framing/random_groups_ext.fits` | exact | 1 |
| Trailing zero/junk blocks and bytes, missing final pad tolerated; `TrailingBytes` reports them; stray NAXISn (n > NAXIS) ignored, also when given twice (`stray_naxis3_twice.fits`, not a duplicate-keyword refusal) | `TestFramingToleranceMatchesCfitsio`, `TestZeroFillTo` | cfitsio; astropy pixels for the missing-pad files | `testdata/images/framing/*.fits` | exact | 1 |
| A truncated trailing XTENSION header is an error (deliberately stricter than cfitsio, like astropy) | `TestTruncatedExtensionHeaderIsAnError` | astropy ('Header missing END card') | `testdata/images/PROVENANCE.md` | exact | 1 |
| A malformed non-structural card does not make the file unopenable; `HeaderOf` returns the error | `TestMalformedCardStillOpens`, `TestOpenToleratesWhatCfitsioOpens` | cfitsio opens and reads the file | `testdata/images/framing/malformed_card.fits`, `header/testdata/malformed.fits` | exact | 1 |
| `CopyTo` reproduces every top-level sample byte for byte (framing only); addresses and CopyTo checked through every open path (`Open`, `OpenReader`, `OpenReadAll` on a non-seeking stream) | `TestCopyToByteIdentical`, `TestHDUAddressesMatchCfitsio` | the committed sample bytes; cfitsio addresses | `testdata/*.fits` | exact | 1 |
| `CopyTo` of a file with bytes after its last HDU drops them, and of a last data unit without its final padding writes the padding as zeros (the only two ways the output differs from the input) | `TestCopyToDropsTrailingBytesAndPads` | `testdata/testprog.std`; cfitsio's padded `f32_plain.fits` | `testdata/testprog.std`, `testdata/images/cfwrite/f32_plain.fits` | exact | 1 |
| Empty / non-FITS input: `ErrEmptyFile`, `ErrNotFITS` (checked before the END scan, also for a file starting with XTENSION) | `TestEmptyFile`, `TestNotFITS`, `TestReadNonexistent` | FITS 4.0 §3, §4.4.1.1 | — | exact | 1 |
| Truncated data unit: `*ErrTruncatedData` with HDU / declared / present bytes; missing BITPIX / NAXIS / NAXISn on a primary or extension, PCOUNT / GCOUNT / TFIELDS on a BINTABLE or TABLE extension, or a non-string XTENSION: `*ErrMissingRequiredKeyword` with HDU and keyword; an IMAGE extension may omit PCOUNT/GCOUNT (cfitsio accepts it) | `TestRefusedFramingMatchesCfitsio`, `TestFramingToleranceMatchesCfitsio` | cfitsio 4.7.0 status for each file (`imgref open`): open refused 222/223/224, move to the extension refused 222/224/228/229/230, read fails 107/108; IMAGE without PCOUNT/GCOUNT read by cfitsio | `testdata/images/refused/*.fits`, `testdata/images/golden/refused__*.open.json`, `testdata/images/framing/img_no_pcount.fits`, `testdata/images/framing/img_no_gcount.fits` | exact | 1 |
| A second SIMPLE, XTENSION, BITPIX, NAXIS, NAXISn (n <= NAXIS, also when a duplicated stray NAXISn with n > NAXIS follows it), PCOUNT, GCOUNT or TFIELDS card with a different value is refused at Open with `*DuplicateKeywordError` / `ErrDuplicateKeyword` naming HDU and keyword, through every open path (cfitsio opens these files: it reads BITPIX/NAXIS/NAXISn by position and an image's PCOUNT/GCOUNT last, so a stray second BITPIX would change the element size); the same keyword repeated with the same value opens and frames the HDUs as declared; a repeat is compared with its type (8 and 8.0 differ), XTENSION trimmed and upper-cased (astrometry.net 5200-series index files written by ESO-QFITS repeat BITPIX..TFIELDS in every kdtree header: the cut-down real index-5205-00, index-5206-38 and the control index-4107 frame byte-identically to cfitsio 4.7.0 fits_get_hduaddrll and read every cell as cfitsio and astropy do); a structural edit of such a header through Header.Set (first card only) is refused by Flush / CopyHDUWithHeader with the cause named and nothing written | `TestDuplicateStructuralKeywordsRefused`, `TestDuplicateStructuralKeywordsSameValue`, `TestDuplicateStructuralKeywordsDifferentValueInTable`, `TestDuplicateStructuralValueComparison`, `TestAstrometryIndexFilesMatchCfitsio`, `TestTablesMatchCfitsio`, `TestTablesMatchAstropy`, `TestEditOfRepeatedStructuralKeywordNamesTheCause` | cfitsio 4.7.0 cfref dump of each file (documented difference: cfitsio opens the different-value files); astropy 8.0.2 | `testdata/images/duplicate/*.fits`, `testdata/images/golden/duplicate__*.json`, `testdata/tables/anet_*.fits` (`testdata/tables/extract_anet.py`, sha256-pinned), `testdata/tables/golden/anet_*.json`, `testdata/tables/golden_astropy/anet_*.json` | exact | 1 |
| Legacy `XTENSION = 'A3DTABLE'` and `'3DTABLE'` read as a binary table, `'IUEIMAGE'` as an IMAGE extension (cfitsio's legacy names) | `TestA3DTABLEIsABinaryTable`, `TestTablesMatchCfitsio`, `TestTablesMatchAstropy`, `TestHDUAddressesMatchCfitsio`, `TestFramingToleranceMatchesCfitsio` | cfitsio 4.7.0 golden (BINARY_TBL / IMAGE_HDU); astropy 8.0.2 (A3DTABLE only; it does not know 3DTABLE) | `testdata/tables/cf_bintable_a3dtable.fits`, `testdata/tables/cf_bintable_3dtable.fits`, `testdata/images/framing/xtension_iueimage.fits` | exact | 1 |
| Out-of-range `HDU(i)`, missing `HDUByName`, every `ReadSubset` bound/stride/length check, `WriteImage`/`AppendImage` length (`ErrShapeMismatch`) and shape checks refuse without writing | `TestValidationBranchesRefuse` | the documented contract (a mutant reading block padding as data is caught) | — | exact | 1 |
| Context: `OpenContext`, `OpenForEditContext`, `EditFileContext` (before anything is touched and per HDU of the scan), `CreateContext`, `ReadPixelsContext` (checked once before the read) | `TestContextFunctions`, `TestReadPixelsCompressedContextCanceled`, `TestReadPixelsCompressedContextCanceledBetweenTiles` | `context.Canceled` | `testdata/testprog.std` | exact | 1 |

### Images (root package)

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| `ReadPixels` float32/float64 targets equal cfitsio TFLOAT/TDOUBLE for every image in the corpus; a BLANK pixel of an integer image reads as its scaled stored value (cfitsio's TDOUBLE read with null checking off) | `TestReadPixelsFloat64MatchesCfitsio`, `TestReadPixelsFloat32MatchesCfitsio` | cfitsio (cfref, imgref goldens incl. `"f64nn"`) | `testdata/images/golden`, `testdata/cfref/golden` | bit-exact; past float32 range ±Inf accepted where cfitsio gives ±FLT_MAX (documented difference) | 1 |
| One BSCALE/BZERO/BLANK rule (`scalingCards`) for every image path: `ReadPixels`, `ReadPixelsMasked`, `CanConvert`, `ImageHDU.BSCALE`/`BZERO`/`BLANK` (now returning an error), `ReadFloat32`, `ReadFrame`, `AsImage`/`image.Decode`/`DecodeConfig`, `OverwritePixels`, `WriteImage` caller cards and tile-compressed images; the same answer from the raw cards, after `HeaderOf` and after an in-memory edit (judged by the bytes Flush would write), and cfitsio reads the flushed file unchanged. The LAST of duplicate cards is used on a plain image, as cfitsio 4.7.0 and astropy 8.0.2 both use it (`dup_bscale`, `dup_bzero`, `dup_blank`, `dup_bscale_bad_first`); an accepted sample reads bit-exactly as cfitsio and as astropy read it. Every card the references read differently, refuse or drop is an error wrapping `ErrInvalidScaling` on every read and on `OverwritePixels` (file unchanged): a value that is not a FITS number (quoted, logical, complex, hex, Inf/NaN, 1E400, lowercase `d`, trailing junk, empty, BLANK `2.0`/`0x2`/overflow, a bad last card), a FITS real cfitsio rejects with ERANGE (subnormal, underflow, just below 2^-1022), no `= ` in columns 9-10 or a lowercase keyword (also wraps `header.ErrInvalidKeyword`), `HIERARCH BSCALE`/`BZERO`/`BLANK` (cfitsio applies, astropy ignores), a BLANK outside the BITPIX range (`blank_wrap16`, `blank_wrap32`: cfitsio compares 32 bits) and BLANK 1234554321 (cfitsio: no BLANK). The generator asserts that no accepted sample is read differently by cfitsio and astropy. cStrtod/ffc2dd/ffc2jj stay proven against glibc and cfitsio (used by the rule) | `TestScalingKeywordsRule`, `TestCStrtodMatchesGlibc`, `TestC2DDAndC2JJMatchCfitsio` | cfitsio 4.7.0 cfref dump + imgref pix of each sample, live cfref dump after Flush; astropy 8.0.2 scaled `getdata` (`gen.py gen_keyword_rule`, rule and both readings in `keywords_rule.json`); glibc 2.44 `strtod`/`strtoll` (`imgref strtod`, 3084 inputs); cfitsio's own `ffc2dd`/`ffc2jj` called directly (`imgref c2`, 1601 inputs) | `testdata/images/keywords/*.fits`, `testdata/images/golden/keywords__*.json`, `testdata/images/golden/keywords_rule.json`, `testdata/images/golden/strtod.txt`, `testdata/images/golden/c2.txt` | bit-exact; exact (error) | 1 + 2 live |
| Integer targets accepted only when they hold every physical value (else `ErrTypeMismatch`); int64/uint64 exact without float64 | `TestReadPixelsIntegerTargets`, `TestUint64ImagesReadWithoutHeaderParse` | cfitsio TLONGLONG/TULONGLONG; exact range in math/big | `testdata/images/cfwrite`, `testdata/images/astropy_written` | exact | 1 |
| BLANK/NaN mask equals cfitsio null flags | `TestReadPixelsMaskedMatchesCfitsio` | cfitsio `fits_read_imgnull` | `testdata/images/astropy_data/blank.fits`, `testdata/images/golden` | exact | 1 |
| Reads after an in-memory header edit use the edited BSCALE/BZERO/BLANK | `TestReadsAfterHeaderParseMatchCfitsio` | cfitsio goldens | `testdata/images/golden` | bit-exact | 1 |
| `ReadSubset` (N-D, strides, int64, BZERO unsigned, float32) | `TestReadSubsetMatchesCfitsio`, `TestReadSubset2D`, `TestReadSubsetStride` | cfitsio `fits_read_subset` | `testdata/images/golden/subsets.json` | bit-exact | 1 |
| `WriteImage` data unit byte-identical to cfitsio's for the 53-case matrix (every BITPIX, unsigned/int8 conventions, caller BSCALE/BZERO incl. a caller header with only BZERO or only BSCALE (the other defaults to 1 / 0; `i16_bzero_only`, `i16_bscale_only`, `f64_bzero_only`, required by name), negative and non-integer on every integer Go type incl. uint16/uint32/uint64 and int8, values up to both edges of the stored range incl. lo-0.5+1 rounding to lo on BITPIX 8/16/32 and the int64 edges under BSCALE 0.5, float32/float64 signaling and payload NaNs written verbatim; every case written by cfitsio with status 0) | `TestWriteImageMatchesCfitsioBytes`, `TestWriteImageMatrixInRange`, `TestWriteImage_IntegerBZERORoundTrip`, `TestWriteImage_FloatNoClamping` | cfitsio-written files `fits_create_img`+`fits_write_img` | `testdata/images/cases.json`, `testdata/images/cfwrite` | byte-exact; header bytes NOT compared | 1 |
| `WriteImage` and `OverwritePixels` never clamp: a value that cannot be stored returns `ErrNumOverflow` and nothing is written (file unchanged), exactly per value where cfitsio's `fits_write_img` reports NUM_OVERFLOW (412) on both edges of BITPIX 8/16/32/64 (incl. lo-0.5, unsigned/int8 conventions with negative and non-integer BSCALE), plus where cfitsio stores a wrong value with status 0 (BITPIX 64 scaled value exactly 2^63, stored INT64_MIN; -2^63 with BZERO 1000; float BITPIX scaled to ±Inf); every other value of the committed set stores cfitsio's bytes (not BITPIX 64 under an integral BZERO other than 0 and 2^63 where v - BZERO needs more than 53 bits: see Deliberate differences) | `TestWriteImageOverflowMatchesCfitsio`, `TestOverwritePixelsOverflowMatchesCfitsio`, `TestWriteImage_IntegerBZERORoundTrip` | each value written alone by cfitsio 4.7.0 (`imgref write`): status and stored element | `testdata/images/golden/write_overflow.json` | exact (status per value, stored bytes) | 1 |
| `WriteImage` with the header of each keywords/ sample the scaling-card rule accepts as the caller header (D exponents, `+2.5`, `.5`, `2.`, `1E1`, BLANK `+2`) stores exactly cfitsio's `fits_write_img` bytes under that header; a BSCALE under which the values overflow (2^-1022) is `ErrNumOverflow`; the headers the rule refuses are refused (`TestWriteImageRefusesInvalidScaling`) | `TestWriteImageKeywordsMatchCfitsio` | cfitsio `fits_write_img` over each keywords/ sample (`imgref overwrite`) | `testdata/images/golden/keywords_overwrite.json`, `testdata/images/keywords` | byte-exact | 1 + 2 live |
| Written images read back identically by cfitsio, fitsverify, astropy | `TestWriteImageCfitsioReadsBack`, `TestWriteImageFitsverify`, `TestWriteImageAstropyReadsBack`, `TestIntegrationWriteThenReadSemantic` | cfitsio 4.7.0, fitsverify, astropy 8.0.2 | `testdata/images/cases.json` | exact | 1 + 2 live |
| `OverwritePixels` applies the BSCALE/BZERO inverse like `WriteImage`, also after unflushed header edits; over every keywords/ header the scaling-card rule accepts it stores exactly what cfitsio's `fits_write_img` stores over the same HDU (last of duplicate cards), refusing the headers the rule refuses (`ErrInvalidScaling`), BSCALE 0 (`ErrInvalidScaling`) and values that overflow under the HDU's BSCALE (`ErrNumOverflow`, where cfitsio reports 412), file unchanged | `TestOverwritePixelsMatchesCfitsio`, `TestOverwritePixelsAfterInMemoryHeaderEdit`, `TestOverwritePixelsKeywordsMatchCfitsio`, `TestOverwritePixels`, `TestOverwritePixelsShapeMismatch` | cfitsio-written data units; cfitsio `fits_write_img` over each keywords/ sample (`imgref overwrite`); `keywords_rule.json` | `testdata/images/cfwrite`, `testdata/images/golden/keywords_overwrite.json`, `testdata/images/golden/keywords_rule.json` | byte-exact | 1 + 2 live |
| `WriteImage` refuses (`ErrInvalidScaling`, nothing written) BSCALE = 0 as a real or integer card (cfitsio `fits_set_bscale`: ZERO_SCALE 322), a second BSCALE/BZERO/BLANK card, every caller card the scaling-card rule refuses (keywords/ samples per `keywords_rule.json`) and a scaling card that cannot be encoded (e.g. a value edited to NaN in `Cards()`); every other keywords/ header writes and cfitsio's live reading of the file equals the library's; `OverwritePixels` refuses an HDU whose BSCALE is 0 (also from an in-memory edit) and an in-memory scaling card that cannot be encoded, file unchanged | `TestWriteImageRefusesInvalidScaling`, `TestOverwritePixelsRefusesZeroBSCALE`, `TestOverwritePixelsRefusesUnencodableScaling` | cfitsio 4.7.0 `fits_set_bscale` status; keywords/ goldens; live cfref dump | `testdata/images/golden/write_refused.json`, `testdata/images/keywords`, `testdata/images/golden/keywords_rule.json` | exact | 1 + 2 live |
| uint64 images (integer card BZERO = 9223372036854775808, parsed as TypeBigInt): header edit, `WriteChecksum` and `Flush`, after which cfitsio reads every value exactly, verifies the checksums and fitsverify is clean | `TestUint64HeaderEditFlushCfitsioReads`, `TestUint64ImagesReadWithoutHeaderParse`, `TestWriteChecksumCfitsioVerifies` | cfitsio TULONGLONG (imgref pix), `fits_verify_chksum`, fitsverify | `testdata/images/cfwrite/u64_conv.fits`, `testdata/images/astropy_written/u64.fits`, `testdata/images/astropy_written/u64_cube.fits` | exact | 1 + 2 live |
| `OverwritePixels` refuses every Go type other than the one `WriteImage` uses for the BITPIX (`*TypeMismatchError`) and a wrong length (`ErrShapeMismatch`), file unchanged | `TestOverwritePixelsTypeAndShapeMatrix` | the documented contract over every cfitsio-written BITPIX | `testdata/images/cfwrite` | exact | 1 |
| Read-only files (`Open`, `OpenReader`, `OpenReadAll`) refuse Flush after a header edit, OverwritePixels, WriteImage, AppendImage, AppendBinaryTable, AppendCompressedImage, WriteChecksum with `ErrReadOnly`, file unchanged; OverwritePixels also refuses ModeCreate | `TestReadOnlyRefusesEveryWrite` | the documented contract, every top-level sample | `testdata/*.fits` | exact | 1 |
| CHECKSUM/DATASUM: `ChecksumStatus`, `VerifyChecksum`, `WriteChecksum`; ModeCreate refused; DATASUM parsed like cfitsio's ffvcks (atof prefix, truncation toward zero: a fraction, trailing junk, hex, leading blanks, -0.5 on a zero sum match; sum + 2^32, -1, -2^32 do not); a DATASUM or CHECKSUM card with no value is `ErrUndefinedChecksum` (cfitsio: VALUE_UNDEFINED 204); the cards found as ffvcks/ffgcrd find them: the first CHECKSUM, then the first DATASUM after it, wrapping to the top (duplicate DATASUM before/after CHECKSUM, none after it, the read one undefined, a blank first CHECKSUM), names compared in upper case (lowercase `datasum` right and wrong), and 2^32 on a zero sum does not match | `TestChecksumMatchesCfitsio`, `TestWriteChecksumCfitsioVerifies`, `TestWriteChecksumRefusesCreateMode`, `TestComputeLargeBufferMatchesCfitsio`, `TestEncodeMatchesCfitsio` | cfitsio `fits_get_chksum` / `fits_verify_chksum` / `fits_encode_chksum`; fitsverify | `testdata/images/checksum`, `testdata/images/astropy_data/checksum.fits` | exact | 1 + 2 live |
| Real samples: real uint16 detector frame, astropy- and cfitsio-written int8/uint32/uint64 | `TestReadPixelsFloat64MatchesCfitsio`, `TestHDUAddressesMatchCfitsio` | cfitsio goldens | `testdata/images/astropy_data/o4sp040b0_raw.fits`, `testdata/images/astropy_data/fixed-1890.fits`, `testdata/images/astropy_written/u64.fits` | bit-exact | 1 |

### Headers (`fits/header`)

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| Every card of every sample decodes to cfitsio's key, value, type and comment (malformed, adversarial, HIERARCH, complex, -0, buffer limits) | `TestDecodeMatchesCfitsioHeaderSamples`, `TestDecodeMatchesCfitsioRepoSamples`, `TestLongstrGoldensCoverRepoSamples`, `TestHeaderGettersMatchCfitsio`, `TestLookupHistoryCommentsMatchCfitsio` | cfitsio `ffpsvc`, `fits_read_key*`, `fits_read_key_longstr` (hdrref, cfref) | `header/testdata/*.fits`, `header/testdata/golden`, `header/testdata/golden/testdata` | exact; documented divergences asserted exactly (see below) | 1 |
| CONTINUE joining as `fits_read_key_longstr`; joined comment kept whole past 72 characters | `TestJoinedContinueCommentIsWholeComment` | cfitsio; past 72 chars the writer input in `hdrref.c` | `header/testdata/cfwritten.fits`, `header/testdata/continue.fits` | exact | 1 |
| Unmodified headers re-encode byte for byte; `Flush` and `CopyHDUWithHeader` no-op are byte-identical | `TestReencodeUnmodifiedIsByteIdentical`, `TestFlushWithoutMutationIsByteIdentical`, `TestCopyHDUWithHeaderUnmodifiedIsByteIdentical`, `TestFlushWithoutChangesIsByteIdentical` | original file bytes; cfitsio HDU offsets | `testdata/*.std`, `header/testdata/adversarial.fits` | exact | 1 |
| In-memory cards read back identically by cfitsio and astropy (long strings, HIERARCH long strings, commentary splitting, floats, uint64 as `TypeBigInt`); a CONTINUE'd comment is split only at single blanks (one with only double blanks is refused) | `TestEncodeReadByCfitsio`, `TestEncodeReadByAstropy`, `TestEncodeStringsReadBackByGo`, `TestEncodeFloatAlwaysFloatSyntax`, `TestIntOverflowGetters`, `TestCommentWordTooLongToContinue`, `TestLongCommentWithDoubleBlanksRefused`, `FuzzEncodeString` | cfitsio 4.7.0 (hdrref), astropy 8.0.2, fitsverify | `header/testdata/golden` | exact (floats bit-exact via cfitsio / `float.hex()`) | 1 + 2 live |
| Write rules: keyword validation equals cfitsio `fits_test_keyword` (refusals are `ErrInvalidKeyword`); NaN/Inf, non-ASCII (`ErrInvalidText`), CONTINUE, bad HIERARCH spacing, commentary value+comment refused; RawValue accepted only when it reads back (an empty text and a `(` without `)` refused with their own reason); a CONTINUE card built in memory refused by AddCard, EncodeCard and Encode whatever its type and letter case; a `TypeBigInt` card with a nil value or one inside int64 refused (cfitsio cannot read the `<nil>` it would write; `Int` would overflow on an ordinary integer); a lower-case keyword written from fields is upper-cased by every path, byte-identical to cfitsio's `fits_write_key`; a card whose keyword upper-cases to END is written only as the bare END card (a lower-case `end` card, which cfitsio reads as an ordinary card, or END with a comment is refused by EncodeCard, EncodeCardWithContinuation, Encode and AddCard); Encode refuses joining cards and names the swallowed card; every encoder guard has a case (integer complex with fractional parts, HIERARCH key with no room for a string, END/commentary with a value, value/type mismatches, `EncodeCard` on a string needing CONTINUE, `Bool` on a non-logical card) | `TestKeywordValidationMatchesCfitsio`, `TestEncodeRejectsWhatFITSCannotHold`, `TestNonASCIIIsInvalidFITS`, `TestSetRefusesWhatWouldNotReadBack`, `TestRawValueRefusedUnlessItReadsBack`, `TestRawValueReadBackByGoAndCfitsio`, `TestEncodeRefusesJoiningCards`, `TestEncodeJoinErrorNamesTheSwallowedCard`, `TestInMemoryContinueRefusedOnEveryPath`, `TestBigIntCardValueRefusedUnlessOutsideInt64`, `TestLowerCaseKeyWrittenUpperCasedLikeCfitsio`, `TestRawValueRefusalNamesItsReason`, `TestEndKeywordNeverWrittenFromFields`, `TestAddCard`, `FuzzEncodeRawValue` | cfitsio `fits_test_keyword` (`header/testdata/golden/testkey.json`), cfitsio `fits_write_key` (`hdrref writekeys`), cfitsio read-back, fitsverify | `header/testdata/golden/testkey.json` | exact | 1 + 2 live |
| Writers keep cards decoded from a file as their source bytes and refuse invalid in-memory cards / non-ASCII column metadata | `TestWritersKeepCardsReadFromAFile`, `TestWritersRefuseInvalidUserCards` | the source card bytes (their meaning pinned to cfitsio above) | `header/testdata/malformed.fits` | exact | 1 |
| A copy of a decoded card whose Key, Type, Value or Comment was edited (Raw kept) is written from its fields, by Encode, AddCard and every writer (WriteImage scales with the BZERO it writes); its bytes are used only while they decode to the card, compared exactly: key letter case, float bits (sign of zero), both complex parts, a switch between integer and float complex (an integer complex with fractional parts is refused, an integer one switched to float is written as float and reads back so); a decoded `-0` integer edited to another value reads that value through `Float` (so `WriteImage` scales with it), one whose comment alone was edited stays -0 and is written `-0`; an edited commentary copy past 72 characters is split as cfitsio splits it; an edited copy whose Type no longer fits its Value is refused (one accepted difference: a float over an integer literal of exactly the same value and sign, the writers' BZERO = 2^63 / BSCALE = 1); an edited card Set could not build is an error; FromCards and AddCard store an edited copy under its upper-cased keyword, so a writer leaves out an edited lower-case structural copy (`naxis1`, `bitpix`, `naxis`, `naxis2`, `tfields`) as it leaves out an upper-case one | `TestEditedDecodedCardIsWrittenAsEdited`, `TestEditedComplexTypeSwitchReadsBackAsHeld`, `TestFromCardsUpperCasesEditedKeys`, `TestWritersLeaveOutEditedLowerCaseStructuralCopies`, `TestRawKeptOnlyWhenItReadsBackAsTheCard`, `TestEditedCardThatCannotBeWrittenIsAnError`, `TestReplacedRawIsNotWrittenWithOldContinueCards`, `TestWritersWriteEditedCardCopies`, `TestEditedNegZeroIntegerReadsBackAsHeld`, `TestWriteImageEditedNegZeroBZERO`, `TestEditedFloatZeroSignReadsBackAsHeld`, `TestEditedComplexPartReadsBackAsHeld`, `TestEditedKeyCaseIsWritten`, `TestEditedCommentaryCopyPastOneCard`, `TestEditedCopyWithMismatchedValueRefused` | cfitsio 4.7.0 reading the written bytes (cfref dump, `hdrref longstr`, TDOUBLE pixels), astropy 8.0.2 (image data), fitsverify | — | exact (floats bit-exact, integers exact) | 1 + 2 live |
| Unsupported, refused with `ErrInvalidKeyword`: copying a decoded card whose (non-HIERARCH) keyword is not upper case (`bzero`, `naxis1`, `object`, `end`) into another header as its bytes, by `AddCard` and so by WriteImage, AppendBinaryTable and AppendCompressedImage; readers disagree on such a card (cfitsio's pixel reader does not apply a `bzero` card, astropy does; fitsverify reports an error), so no copy is right. Upper-casing the copy's Key writes it from its fields; Encode of an unedited parsed header still writes it as its bytes; a lower-case HIERARCH path is copied (fitsverify accepts it) | `TestAddCardRefusesUneditedLowerCaseKeyword`, `TestWritersRefuseUneditedLowerCaseCards` | cfitsio 4.7.0 TDOUBLE pixels and astropy 8.0.2 data of the same source bytes, fitsverify | — | exact | 1 + 2 live |
| Integer and logical cards are byte-identical to cfitsio's `fits_write_key` (value right-justified to column 30, every printed width 1-20, with and without comment) | `TestFixedFormatMatchesCfitsioBytes` | cfitsio 4.7.0 (`hdrref writekeys`) | `header/testdata/hdrref/hdrref.c` | exact (card bytes) | 1 + 2 live |
| Header fuzzers compare decoded values, not counts; every seed sample is required by name | `FuzzParseHeader`, `FuzzCardDecode` | decode/encode agreement (meaning pinned by the goldens) | `testdata/testprog.std` | — | 1 |
| Case-insensitive lookup for Get/Has/String/Int/Float/Bool/Set/Delete; Set on a name in another case replaces that card; with a duplicated keyword every getter returns the FIRST card, for headers read from a file (FromCards), built with Add or copied with AddCard, and again after each Delete, which removes the first (also for HISTORY, including a HISTORY that is card 0); absent keys give `ErrKeyNotFound`, wrong types `ErrWrongType` | `TestLookupCaseAndFirstCardMatchCfitsio`, `TestDeleteFirstCardMatchesCfitsio`, `TestDuplicateKeywordFirstCardMatchesCfitsio`, `TestDeleteCommentaryRemovesFirstMatchesCfitsio`, `TestTypedGettersErrKeyNotFound` | cfitsio 4.7.0 `fits_read_keyword` / `fits_delete_key` from the top of the encoded header (`hdrref readkey` / `delkey`) | `header/testdata/hdrref/hdrref.c` | exact (card bytes) | 1 + 2 live |
| Clone independence, BlockCountFor | `TestCloneIsIndependent`, `TestBlockCountForMatchesEncode` | definitional / Encode | — | exact | 1 |

### Tables (root package, `internal/tform`)

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| Every fixed TFORM (L X B I J K A E D C M), scalar and vector, read correctly; TNULL masks; TSCAL/TZERO incl. unsigned conventions; a cell we call null is NaN / TNULL in astropy (±Inf stays a value) | `TestTablesMatchCfitsio`, `TestTablesMatchAstropy`, `TestTableReadExactnessErrors` | cfitsio 4.7.0 goldens; committed astropy 8.0.2 goldens (every column astropy can read; what it cannot is recorded with its error text and checked) | `testdata/tables/*.fits`, `testdata/tables/golden`, `testdata/tables/golden_astropy` | exact for integers; floats bit-exact | 1 |
| Astropy goldens reproduce from `astropy_dump.py --golden`; large columns as exact sha256 digests plus samples | `TestAstropyGoldensReproduce`, `TestTableSamplesReproduce` | astropy 8.0.2 | `testdata/tables/astropy_dump.py`, `testdata/tables/golden_astropy` | byte-exact | 2 |
| Every numeric column of every sample read as every Go target type (int8..uint64, float32, float64); integer targets never wrap on the exact paths; float32 targets get the float32 nearest the exact value on the exact paths (TSCAL 1, integer TZERO, unscaled: no double rounding through the float64) and the float32 of the float64 arithmetic otherwise | `TestTablesMatchCfitsio`, `TestTargetRangesSample`, `TestTablesFloat32MatchCfitsio` | cfitsio stored integers + TZERO literal, exact math/big; cfitsio 4.7.0 TFLOAT reads of every numeric column (`cf_tfloat.c`): bit-identical except cfitsio's flush/±Inf-null, its ±FLT_MAX overflow, and cells where its TFLOAT is the float32 of its own TDOUBLE (pinned per sample) | `testdata/tables/cf_bintable_ranges.fits`, `testdata/tables/cf_bintable_tzero_round.fits`, `testdata/tables/golden_tfloat` | exact; float32 bit-exact | 1 |
| Integer TZERO literals read exactly (incl. 2^63-1, 2^62+1, 2^63+1, 2^64, uint64 TZERO = 2^63 on K, 2^63 on J); TSCAL/TZERO literals with a D exponent (`4.611686018427387905D18`, `1.0D-1`, `-2.5D0`) read as exact decimals; a float literal near but not at 2^63 refused for integer targets; the writer writes TZERO 2^63 and every integral TScal/TZero as the exact integer literal and any other value as the shortest decimal that reads back as the float64 ("0.1", "0.3333333333333333"; read back as that decimal exactly), refuses an integral value whose literal does not fit a card (70 characters, about 1e70 and beyond) and refuses hdr TSCAL/TZERO literals it cannot write back | `TestTablesMatchCfitsio`, `TestAmbiguousTZERONear2p63`, `TestWriterScalingLiteralsExact`, `TestWriterNonIntegralScalingLiteral`, `TestWriterHugeScalingLiterals`, `TestWriteTableInteropCfitsio` | cfitsio raw values + literal text (math/big); cfitsio's value string of the written TSCAL/TZERO cards | `testdata/tables/cf_bintable_tzero_big.fits`, `testdata/tables/cf_bintable_kzero.fits`, `testdata/tables/cf_bintable_uint64.fits`, `testdata/tables/cf_bintable_scaled_exact.fits`, `testdata/tables/cf_bintable_scaled_special.fits` | exact | 1 + 2 live |
| Integer targets on the float64 scaling path (TSCAL != 1, fractional TZERO) get the exact value stored*TSCAL + TZERO (literals as exact decimals) or ErrTypeMismatch, never the rounded float64; float targets get the float64 arithmetic (cfitsio's TDOUBLE); binary and ASCII. E/D under TSCAL 1 and an integer TZERO no double holds (2^53+1, 2^64+2049) give stored + TZERO correctly rounded, not TZERO's double plus stored (cfitsio and astropy round twice; pinned). ±Inf stored in a scaled E/D column (scalar and vector) is ±Inf*TSCAL + TZERO for float targets (astropy) and ErrTypeMismatch for integer targets; NaN is null | `TestTablesMatchCfitsio`, `TestTablesMatchAstropy`, `TestTargetRangesSample` | cfitsio raw stored integers / E-D bytes / ASCII field text + card literals, math/big; cfitsio TDOUBLE for float targets; astropy 8.0.2 for the scaled ±Inf cells (cfitsio flags them null) | `testdata/tables/cf_bintable_scaled_exact.fits`, `testdata/tables/cf_ascii_scaled_exact.fits`, `testdata/tables/cf_bintable_ranges.fits`, `testdata/tables/cf_bintable_scaled_special.fits`, `testdata/tables/cf_bintable_tzero_round.fits` | exact | 1 |
| ASCII tables: exact integers at any size, implied decimals, blanks, TNULL per §7.2.2 (also a TNULL longer than the field), digitless fields, fractions under an integer TZERO (exact sum, correctly rounded; also TZERO beyond 2^64 at a rounding tie, fields of 802 digits, and exponents below -400 whose sign decides a TZERO tie: "1.5E-900" + (2^53+1) is 2^53+2), integers beyond int64 / uint64 under an int64 or bigger TZERO (exact sum, correctly rounded: 18446744073709553665 - 2 is 2^64, not 2^64+4096), negative zeros ("-0", "-0.0", "-", "-.") under integer, big, fractional TZERO and TSCAL 3 (exactly TZERO), fractions never read into integer targets at any length or exponent; TBCOL 0 and a TBCOL or field width near 2^63 refused (no int64 wrap past the NAXIS1 check) | `TestTablesMatchCfitsio`, `TestTablesMatchAstropy`, `TestASCIIRejectsBadFields` | FITS 4.0 §7.2 field text evaluated with big.Rat (cfitsio goes through a double and is cross-checked within 2 ulp of the result and operands) | `testdata/tables/cf_ascii_exact.fits`, `testdata/tables/cf_ascii.fits`, `testdata/tables/cf_ascii_implied.fits`, `testdata/tables/cf_ascii_scaled.fits`, `testdata/tables/cf_ascii_scaled_exact.fits`, `testdata/tables/cf_ascii_scaled_edge.fits`, `testdata/tables/cf_ascii_tzero_round.fits` | exact | 1 |
| `VarColumn.Rows()` and `Raw()` are the same zero-copy views as `At()`; `ColumnByName` (binary and ASCII) finds each cfitsio TTYPE exactly and is case-sensitive | `TestVarColumnViewsAgree`, `TestColumnByNameIsCaseSensitive` | `At()` pinned to cfitsio by `TestTablesMatchCfitsio`; cfitsio TTYPE cards | `testdata/tables/cf_vla_p.fits`, `testdata/tables/cf_vla_q.fits` | exact | 1 |
| P and Q VLAs of every element type; heap descriptor validation (negative count, negative offset, beyond heap, Q overflow) and the heap / THEAP limits at their boundaries (payload ending at the heap end vs one byte past, offset equal to the heap size, THEAP at / one past the data unit end, THEAP inside the rows, a THEAP gap, a non-integer THEAP); 0P/0Q columns | `TestTablesMatchCfitsio`, `TestVarColumnRejectsBadHeap`, `TestVLAHeapBoundaries`, `TestZeroRepeatVLA`, `TestPutDescriptorRange`, `TestVarColumnRoundTripFloat32`, `TestMixedFixedAndVarColumns` | cfitsio goldens; FITS 4.0 §7.3.5 arithmetic with fitsverify's committed verdicts (it agrees, negative offsets and THEAP = 12.0 included, except for the THEAP gap and a THEAP past the data unit, which it does not check; cfitsio reads past the heap silently); 0P/0Q: no reference tool reads them | `testdata/tables/cf_vla_p.fits`, `testdata/tables/cf_vla_q.fits`, `testdata/tables/vla_theap_gap.fits`, `testdata/tables/astropy_theap-gap.fits`, `testdata/tables/heap/*.fits`, `testdata/tables/heap/fitsverify.txt` | exact | 1 |
| E/D/C/M keep ±Inf, -0, subnormals (float32 and float64), the largest finite values; NaN is the E/D null; fixed and VLA; TSCAL/TZERO on complex is an error, because neither reference follows FITS 4.0 7.3.2 (real part TSCAL*re + TZERO, imaginary part TSCAL*im): cfitsio 4.7.0 adds TZERO to the imaginary part too and astropy 8.0.2 drops the imaginary part | `TestSpecialValuesSample`, `TestTablesMatchAstropy`, `TestComplexSpecialValuesAstropy`, `TestComplexScalingRejected`, `TestComplexScalingReferencesDisagree` | astropy 8.0.2 writer and committed astropy read (astropy-written sample); astropy read-back of Go-written cells; committed cfitsio and astropy reads of a scaled 1C/1M sample against the standard's formula on the stored bytes | `testdata/tables/astropy_only/astropy_special.fits`, `testdata/tables/golden_astropy/astropy_special.fits.json`, `testdata/tables/refused/cf_bintable_cplx_scaled.fits` | bit-exact | 1 + 2 live |
| Binary writer: unsigned TZERO convention, TSCAL/TZERO/TNULL/TUNIT/TDISP written (TNULL at the top of the B range and TNULL = 0 read back as null by cfitsio), conflicts and invalid input refused (TNULL one past each stored range; one reject case per guard, incl. no / 1000 columns, VLA repeat, Valid on a VLA column, TDIM on a VLA column, VarValid on a fixed column, complex nulls, unparsable hdr TDIM/TSCAL/TZERO, the 1P descriptor range) | `TestWriteTableInteropCfitsio`, `TestWriteTableInteropAstropy`, `TestWriteTableRejects`, `TestWriteTableRejectsTable`, `TestBinaryTableRoundTrip` | cfitsio read-back, fitsverify, astropy | — | exact | 1 + 2 live |
| ASCII writer: Fw.d, Ew.d, Dw.d (Fw.0 / Ew.0 with an explicit point, as printf `%#.0f` / `%#.0E`), -0, TNULL (one as wide as the field; one character wider refused), fields laid out with one blank between neighbours (cfitsio's TBCOLn / NAXIS1), TUNITn exactly when a Unit is given (also checked for the binary writer in `TestWriteTableInteropCfitsio`); error paths (one reject case per guard, incl. no / 1000 columns, zero or two Data fields, non-printable TNULL, a header-builder refusal, hdr TxxxN beyond the columns) | `TestWriteASCIITableInteropCfitsio`, `TestWriteASCIITableInteropAstropy`, `TestWriteASCIITableRejects` | fitsverify, cfitsio raw rows vs glibc printf, astropy | — | exact | 1 + 2 live |
| TFORM/TDIM grammar; a repeat count, TDIM product or row width whose byte arithmetic overflows int64 refused by `Columns` (not wrapped past the NAXIS1 check into a reader panic) | `TestParseBinaryGrammar`, `TestParseASCIIGrammar`, `TestParseDimGrammar`, `TestParseMatchesCfitsioColumns`, `FuzzTForm`, `TestColumnLayoutOverflow` | FITS 4.0 Tables 15/18; cfitsio; int64 arithmetic limits | `testdata/tables/golden` | exact | 1 |
| Invalid logical bytes (not `T`, `F` or 0) are an error from `ReadColumn[bool]`, `ReadColumnMasked[bool]` and `ReadLogicalColumn`, scalar and vector (deliberate difference: cfitsio reads 0x01 as true) | `TestInvalidLogicalByteIsAnError`, `TestTablesMatchCfitsio` | cfitsio 4.7.0 golden of the byte-patched cfitsio-written table | `testdata/tables/cf_bintable_badlogical.fits` | exact | 1 |
| String arrays: TDIM (w, d2, ...) on rA and rAw substring arrays read per element (`ReadStringArrayColumn`; NUL-padded astropy elements kept); a TDIM on an rAw column wins over w (8A2 with TDIM (8) is one string, with (4,2) two of 4; cfitsio and astropy agree); `ReadColumn[string]` refuses a multi-string cell; layouts the references disagree on (TDIM covering fewer characters, r not a multiple of w, TDIM string arrays in VLAs) and writing string arrays are errors; `RowBytes` returns the raw row | `TestTablesMatchCfitsio`, `TestTablesMatchAstropy`, `TestStringLayoutsRefused`, `TestWriteTableRejects` | cfitsio goldens (cfitsio-written and astropy-written string arrays, 20A5); committed astropy goldens | `testdata/tables/cf_bintable_strarray.fits`, `testdata/tables/strarray_astropy.fits`, `testdata/tables/cf_bintable_substr.fits` | exact | 1 |
| Real table samples, goldens complete and reproducible | `TestTableGoldensCoverSamples`, `TestTableGoldensReproduce` | cfitsio, astropy, Gaia DR3 | `testdata/tables/gaia_dr3_sample.fits`, `testdata/tables/astropy_chandra_time.fits` | exact | 1; 2: `TestTableGoldensReproduce` |

### Editing (root package, `cmd/`)

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| Header resize Flush produces a correct file, keeps mode and symlink; like an in-place Flush it keeps every byte after the edited header as stored: the bytes after the last HDU (`TrailingBytes`: zero block, junk block, junk bytes) and a last data unit without its final padding (before: the rewrite dropped the former and padded the latter) | `TestFlushHeaderResizeMatchesCfitsio`, `TestFlushResizeThenKeepEditing`, `TestFlushHeaderResizeKeepsModeAndSymlink`, `TestEditHeaderInPlace`, `TestFlushKeepsBytesAfterLastHDU` | cfitsio goldens + live cfref dump, fitsverify; byte equality with the original file after the header | `testdata/cfref/golden`, `testdata/images/framing/trail_*.fits`, `testdata/images/framing/missing_final_pad*.fits` | exact | 1 + 2 live |
| Process crash at every step of Flush/EditFile recovers to exactly original or result; `ErrPendingJournal` before recovery (`Open`, and `OpenReader` of an `*os.File`); legacy/torn/wrong-size journals | `TestFlushHeaderResizeCrashAtEveryStep`, `TestEditFileCrashAtEveryStep`, `TestOpenRefusesPendingJournal`, `TestOpenReaderRefusesPendingJournal`, `TestRecoverLegacyJournalRefuses`, `TestRecoverTornJournalRemoved`, `TestRecoverJournalSizeMismatchRefuses`, `TestRecoveryCrashIsRecoverable` | byte equality with the original and the fault-free result | `testdata/testprog.std` | exact | 1 |
| `OverwritePixels` (also on a data unit without its final padding, where rollback cuts the file back) and the edit-mode appends (image, binary table, ASCII table, tile-compressed integer and float32 image, also over trailing junk) run under the undo journal: a crash after every step and at every 80-byte point of the data write is refused by `Open` until recovery, and recovery gives exactly the original (or the result once committed); a real write/fsync/journal-removal error is rolled back before it is returned (also for the in-place header Flush, which can then be retried) | `TestInPlaceWriteCrashAtEveryStep`, `TestInPlaceWriteIOErrorRollsBack` | byte equality with the original and the fault-free result | `testdata/iter_image.fit`, `testdata/testprog.std`, `testdata/images/cfwrite/f32_plain.fits` (cut to 2924 bytes) | exact | 1 |
| A journal with flag bits this release does not know is refused by `OpenForEdit` and `EditFile`, file and journal untouched | `TestRecoverJournalUnknownFlagsRefuses` | the journal format; byte equality | `testdata/iter_a.fit` | exact | 1 |
| `OverwritePixels` writes the final zero padding when the file lacks it or holds other bytes there | `TestOverwritePixelsRestoresPadding` | cfitsio-written files (with the documented deviations of `expectedCaseData`) | `testdata/images/cfwrite` | byte-exact | 1 |
| The temp-file replacement (resizing Flush, EditFile) refuses a hard-linked file (`ErrHardLinked`), keeps group/owner, extended attributes and POSIX ACLs, does not add ones the temp file inherited (a directory default ACL), and refuses (`ErrCannotPreserve`) when one cannot be kept (Linux) | `TestReplaceRefusesHardLinks`, `TestReplaceKeepsGroupXattrsAndACLKernel`, `TestReplaceDropsInheritedACLKernel`, `TestReplaceKeepsGroupXattrsAndACL`, `TestReplaceDropsInheritedACL`, `TestReplaceRefusesWhenOwnerCannotBeKept` | the kernel's listxattr/getxattr (tier 1, attributes set with setxattr, ACLs in posix_acl_xattr form); getfattr / getfacl (attr, acl) in tier 2; stat; byte equality for refusals | `testdata/testprog.std` | exact | 1; 2: `TestReplaceKeepsGroupXattrsAndACL`, `TestReplaceDropsInheritedACL` |
| macOS: the same replacement refuses a hard-linked file, keeps the group and every listed extended attribute byte for byte (a user attribute, `com.apple.FinderInfo`, `com.apple.quarantine`, a 70000-byte `com.apple.ResourceFork`), removes attributes the source lacks, refuses an owner it cannot keep; an ACL is NOT carried over (pinned known gap) | `TestReplaceRefusesHardLinks`, `TestReplaceKeepsGroupAndXattrsDarwinKernel`, `TestCopyXattrsMakesDestinationEqualDarwinKernel`, `TestReplaceKeepsGroupAndXattrsDarwin`, `TestCopyXattrsMakesDestinationEqualDarwin`, `TestReplaceRefusesWhenOwnerCannotBeKept`, `TestReplaceDropsACLDarwin` | the kernel's listxattr/getxattr and stat (tier 1); /usr/bin/xattr -lx, stat -f %g, ls -le (tier 2); byte equality for refusals | `testdata/testprog.std` | exact; **pending Mac run** (the tier 1 twins are vetted for darwin, not yet run there) | 1; 2: `TestReplaceKeepsGroupAndXattrsDarwin`, `TestCopyXattrsMakesDestinationEqualDarwin`, `TestReplaceDropsACLDarwin` |
| One editor at a time: a second `OpenForEdit` / `EditFile` while one is open or running (also after a resizing Flush, and when it opened the file before the first editor's rename and locked it after) fails with `ErrLocked` before any recovery, so a running editor's journal and temp file are never touched; with no editor left it edits the file now at the path | `TestEditLockRefusesConcurrentEditors`, `TestEditLockFollowsReplacedFile` | the documented contract (flock(2) on Linux, an OFD record lock on macOS); byte equality; the bytes at the path | `testdata/testprog.std` | exact (Linux); macOS **pending Mac run** | 1 |
| An in-place write (Flush without resize, `OverwritePixels`, edit-mode appends) while the file is open through `Open` or `OpenReader(*os.File)` fails with `ErrOpenForReading` and writes nothing; after the reader closes it gives the fault-free result; a resizing Flush leaves the reader on the old, complete file | `TestInPlaceWriteRefusedWhileReaderOpen` | byte equality with the original and the fault-free result (OFD locks, Linux and macOS) | `testdata/testprog.std`, `testdata/iter_image.fit` | exact (Linux); macOS **pending Mac run** | 1 |
| The same locks across processes (the child is the test binary re-run as a helper): a second editor in the other process (`OpenForEdit`, `EditFile`, `Create`) gets `ErrLocked`, also after a resizing Flush, both ways round; a reader in one process makes the other's in-place Flush and `Create` fail with `ErrOpenForReading` without stopping its `OpenForEdit`, both ways round; an `Open` in the other process waits for a running in-place write and then reads the new header; closing another descriptor of the file (an `os.ReadFile`, the caller's file after `OpenReader`) keeps every lock this process holds through its other descriptors (classic POSIX record locks fail this: checked by switching the Linux commands to F_SETLK/F_SETLKW, which fails `reader-vs-writer` and `other-close-keeps-locks`) | `TestLocksAcrossProcesses` | the documented errors; byte equality after refusals; the header the waiting reader read | `testdata/testprog.std` | exact (Linux); macOS **pending Mac run** | 1 |
| On a file system without record locks (injected ENOLCK / EINVAL / ENOTSUP / EOPNOTSUPP on the read/write lock) `Open` still reads and every in-place write is refused, nothing written | `TestNoRecordLocksRefusesInPlaceWrites` | byte equality | `testdata/testprog.std` | exact (Linux); macOS **pending Mac run** | 1 |
| After a Flush (in place and resize) `CopyTo` and an EditFile copy out of the flushed `*File` emit the bytes on disk, which cfitsio reads with the new EXTNAMEs | `TestFlushedFileCopiesAndIndexesNewHeaders` | bytes on disk, live cfref dump | `testdata/testprog.std` | exact | 1 + 2 live |
| fsync order of the journal (Flush, OverwritePixels, appends) and rename protocols (on macOS each sync is F_FULLFSYNC, falling back to fsync where the file system does not support it: `syncFile`) | `TestFsyncOrder` | the protocol's ordering rules (no power cut simulated) | — | exact order (Linux); macOS **pending Mac run** | 1 |
| The journal of a file is found whatever path spelling reaches it: through a symlinked directory, the real directory, a symbolic link to the file in another directory, and the fully resolved spelling (macOS /var -> /private/var); on a case-insensitive file system through another case (skipped on a case-sensitive one) | `TestJournalFoundThroughSymlinkedDir`, `TestJournalFoundWithOtherCase` | `ErrPendingJournal` from every spelling, byte equality after recovery (naming the journal after the unresolved path fails `crash-via-file-link`) | `testdata/testprog.std` | exact (Linux; the case test skips there); macOS **pending Mac run** | 1 |
| Flush / CopyHDUWithHeader refuse structural changes and anything the encoder cannot represent exactly; each lossy-looking edit has a pinned outcome: CONTINUE long strings (kept and newly set) and integer-valued floats are written and read back by cfitsio, commentary text longer than one card is refused by the writer (card count), NaN and a keyword named END by `Set`, subnormal floats by the per-card read-back check (cfitsio cannot read them either: status 412) | `TestFlushRefusesStructuralChange`, `TestFlushRefusesLossyEncoding`, `TestEditFileCopyHDUWithHeaderRefusesLossyEncoding`, `TestFlushEditMatrixMatchesCfitsio` | cfitsio golden KEY_PKLS, live cfref dump (`fits_read_key_longstr`, value type F), fitsverify | `testdata/cfref/golden/testprog.std.json` | exact | 1 + 2 live |
| Flush (and `CopyHDUWithHeader`), and every writer (`WriteImage`, `AppendImage`, `AppendBinaryTable`, `AppendASCIITable`, `AppendCompressedImage`, `AppendCompressedFloat32Image`, `AppendCompressedFloat64Image`, on Create files and as ModeEdit appends: one rule, `refuseHierarchAlias`, wrapping `header.ErrInvalidKeyword`, no HDU added) refuse a caller HIERARCH card whose path is a standard keyword name (`HIERARCH BSCALE`, `BZERO`, `BLANK`, `END`, `NAXIS1`, `EXTNAME`, `OBJECT`, `TSCAL1`), file byte-identical: cfitsio names it by the short keyword (it scaled pixels 1..8 to 2..16 for `HIERARCH BSCALE = 2`, and read the data at 2880 instead of 5760 for `HIERARCH END`); HIERARCH cards with a real path are written and cfitsio reads the file as we do | `TestFlushRefusesHierarchAliasOfStandardKeyword`, `TestEditFileCopyHDUWithHeaderRefusesHierarchAlias`, `TestWritersRefuseHierarchAlias` | live cfref dump of the flushed file and of the would-be header (ffgknm) | — | exact (pixels, data offset) | 1 + 2 live |
| Cards `Set` refuses (a blank inside the name, non-ASCII, control and DEL bytes) in a header read from disk survive a Flush of another edit byte-identical, and cfitsio reads every earlier card exactly as before | `TestFlushKeepsCardsSetWouldRefuse` | live cfref dump before and after | — | exact | 1 + 2 live |
| Edits through Flush read back by cfitsio as set: HIERARCH string and float, `'abc&'` without CONTINUE, `/` in value and comment, quotes, 2^53+1, MinInt64, MaxInt64, a widened float32, extreme and -0.0 floats, an 80-character string (CONTINUE), `'=x'`, `' '`, `''`, T, F | `TestFlushEditMatrixMatchesCfitsio` | live cfref dump (typed value, comment), fitsverify | `testdata/testprog.std` | exact (float bits) | 1 + 2 live |
| Lock order and release: `Open` during an in-place write waits for it (its non-waiting request refused by the kernel; on Linux the blocked request is in /proc/locks), then reads the new bytes, or, if the writer died, fails with `ErrPendingJournal` (the journal is checked after the lock); `OpenReader`'s `Close`, and an `OpenReader` whose scan fails, release the lock with the caller's file open; `EditFile` keeps the path locked across its rename (`ErrLocked`); journal rollback waits for a reader's lock; `OpenReader` ignores the journal of a path that now names another file | `TestOpenWaitsForInPlaceWriteThenChecksJournal`, `TestOpenReaderReleasesLockWithCallerFileOpen`, `TestEditFileKeepsPathLockedAfterRename`, `TestRollbackWaitsForReaderLock`, `TestOpenReaderIgnoresJournalOfReplacedPath` | the documented contract (OFD locks, and flock(2) on Linux); the kernel's EAGAIN to the non-waiting request (`lockWaitHook`), and on Linux the blocked request in /proc/locks; byte equality | `testdata/testprog.std` | exact (Linux); macOS **pending Mac run** | 1 |
| `Create` does not rewrite a file another edit owns: a pending journal (`ErrPendingJournal`, so recovery cannot write stale bytes into the new file), an open editor (`ErrLocked`) or reader (`ErrOpenForReading`) leave it untouched; while Create's `*File` is open editors get `ErrLocked`; a path renamed over between open and lock is followed | `TestCreateRespectsJournalAndLocks` | byte equality; the pixels the next reader sees | — | exact | 1 |
| `Close` does not flush | `TestCloseDoesNotFlush` | cfref dump, byte equality | — | exact | 1 |
| `CloseAndFlush` whose Flush fails returns the error and leaves the file open: fixing the header and calling it again writes the file | `TestCloseAndFlushLeavesFileOpenOnFlushError` | cfref dump of the result; byte equality after the refused Flush | — | exact | 1 + 2 live |
| EditFile HDU conversion equals cfitsio `fits_copy_hdu` (EXTEND / COMMENT prefixes, reserved space, NAXISn order, NAXIS 0 and NAXIS 3/4/6 primaries to extensions and back); tile-compressed images copied to position 0 (empty primary in front) and before their primary, pixels still decoding to the source's cfitsio values | `TestEditFileCopyMatchesCfitsioCopyHDU`, `TestConvertForPositionRefusesBadAxisOrder`, `TestEditFileRefusesRandomGroupsAsExtension`, `TestCfcopyReproducesGoldens` | cfitsio 4.7.0 `fits_copy_hdu` (cfcopy); cfref golden of the source | `testdata/cfcopy/golden`, `testdata/cfcopy/reserve.fits`, `testdata/cfcopy/comments.fits`, `testdata/cfcopy/golden/comp_skip0.fits`, `testdata/cfcopy/golden/comp_1_0.fits`, `testdata/cfcopy/golden/cube*.fits`, `testdata/cfcopy/golden/naxis0_*.fits` | card for card, data byte for byte, pixels exact | 1 + 2 live; 2: `TestCfcopyReproducesGoldens` |
| New images through `WriteImage(out, ...)` at position 0 (primary) or after copies (IMAGE extension) | `TestEditFileWriteImageAfterCopy`, `TestEditFileWriteImageThenCopyPrimary` | live cfref dump, fitsverify | `testdata/cfref/golden/testprog.std.json` | exact | 1 + 2 live |
| EditFile safety: empty output, data size / kind changes refused, mode kept, symlinks followed, source intact on error/panic | `TestEditFileNoHDURefused`, `TestEditFileCopyHDUWithHeaderRejectsSizeChange`, `TestEditFileCopyHDUWithHeaderRejectsKindChange`, `TestEditFileCopyHDUWithHeader`, `TestEditFileKeepsPermissions`, `TestEditFileFollowsSymlink`, `TestEditFileFailureLeavesSourceIntact`, `TestEditFilePanicLeavesSourceIntact`, `TestEditFileDropHDU` | byte equality, file mode, Lstat, cfitsio | — | exact | 1 |
| Appends in edit mode, Create, CreateWriter (a non-empty stream is refused with `ErrStreamNotEmpty`, nothing written); an extension written first to a fresh Create (binary/ASCII table, compressed image) gets an empty primary in front, as cfitsio's ffcrtb does | `TestAppendInEdit`, `TestCreateTruncatesAndWritesCompleteHDUs`, `TestCreateWriter`, `TestCreateExtensionFirstGetsPrimary` | cfref dump, fitsverify, byte equality | — | exact | 1 + 2 live |
| `cmd/fitsdump` output equals cfitsio's view (also for a lone CONTINUE card on a primary, a BINTABLE and a tile-compressed HDU); undecodable HDU/columns reported in place, exit 1 | `TestDumpMatchesCfitsio`, `TestDumpExitStatus`, `TestDumpReportsColumnError`, `TestDumpLoneContinueMatchesCfitsio`, `TestDumpReportsHDUError` | cfref goldens / live dump | `testdata/cfref/golden` | exact | 1 + 2 live |
| Every header of a file `Open` accepts parses (`HeaderOf` never fails; `Header()`, `Columns()`, `CompressionType()` do not panic), and the parse ends at the END card the scan found (no earlier card named END survives `Open`) | `TestOpenedHeadersAlwaysParse`, `FuzzOpenedHeadersParse` | the property over crafted END/CONTINUE/NUL/non-ASCII cards and 3000 seeded mutations of real samples | `testdata/testprog.std`, `testdata/iter_a.fit`, `testdata/comp_rice_i16.fits`, `testdata/vari.fits` | exact | 1 |
| `cmd/fitscopy` byte copy, refuses copy onto itself, exit 3 on difference (also for trailing bytes and a missing final pad, which `CopyTo` does not reproduce) | `TestCopyEverySample`, `TestCopyOntoItselfRefused`, `TestCopyUsageAndErrors`, `TestCopyDifferenceReported`, `TestCopyReportsInputsCopyToChanges` | input bytes | `testdata/*.fits` | exact | 1 |

### Tile compression: read (root package, `fits/compress`)

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| RICE_1, GZIP_1, GZIP_2, HCOMPRESS_1 (SMOOTH, 64-bit lossy 32-bit), PLIO_1; ZBITPIX 8/16/32/-32/-64; BSCALE/BZERO incl. uint16/uint32/int8 conventions: every algorithm x ZBITPIX x convention cell cfitsio can write has a sample (148 files incl. 4-D and 6-D images, real M13 files; the cells cfitsio cannot write are listed with the reason); BSCALE/BZERO on quantized and fallback tiles and ZSCALE/ZZERO keywords on lossless tiles for both -32 and -64. The funpack twins are not committed (`testdata/gen_twins.sh` regenerates them from the samples); their cfitsio goldens are, once per distinct twin (`testdata/comp_fzun_same.txt`) | `TestCompressedCorpusRead`, `TestCompressedCorpusComplete`, `TestCompressedCorpusReproduce`, `TestCompressedGoldensReproduce` | libcfitsio 4.7.0 decode of each fpack file (`fits_read_subset` for more than 3 axes; funpack twin + astropy for 6-D, where that reader is wrong) + cfitsio's read of the funpack twin | `testdata/comp_fz_*.fits`, `testdata/cfref/golden`, `testdata/cfref/golden/twins`, `testdata/gen_fpack.sh`, `testdata/gen_twins.sh` | bit-exact | 1; 2: `TestCompressedCorpusReproduce`, `TestCompressedGoldensReproduce` |
| PLIO_1 uint16 (BZERO 32768: cfitsio stores value + 32768; BLANK compared with the PLIO value as `fits_read_img` does, also on the ReadFloat32 path); PLIO_1 at ZBITPIX 8 (incl. int8) and 32 | `TestReaderDecodesCfitsioWrittenGoldens`, `TestCompressedCorpusRead`, `TestReadFloat32NullsMatchCfitsioCorpus` | libcfitsio 4.7.0 (`fits_img_compress` + `fits_read_img`); fpack + funpack twin | `compress/testdata/cref/golden/img/plio_u16.fits.gz`, `testdata/comp_fz_u16_plio.fits`, `testdata/comp_fz_u8_plio.fits`, `testdata/comp_fz_s8_plio.fits`, `testdata/comp_fz_i32_plio.fits`, `testdata/comp_fz_u16blank_plio.fits` | exact | 1 |
| The reader decodes every image cfitsio's `fits_img_compress` wrote (91 files: every algorithm, tilings, ZDITHER0 1..10000, dither random-sequence wrap inside 10000-pixel -32 tiles, IEEE specials, NaN nulls, a caller BLANK on a float image, lossy HCOMPRESS integers, 1-D/3-D/4-D/6-D) exactly as cfitsio decodes them; ReadFloat32 on the 8/16-bit ones | `TestReaderDecodesCfitsioWrittenGoldens` | libcfitsio 4.7.0 `fits_read_img` TDOUBLE (the `dec` digest in `<case>.json`: sha256 over every value and null flag plus a sample) or the cfitsio-verified lossless input (regenerated by `internal/crefgen`, sha256-checked) | `compress/testdata/cref/golden/img/*.json`, `compress/testdata/cref/golden/img/*.fits.gz` | bit-exact (float32-rounded for ZBITPIX -32) | 1 |
| Integer targets of compressed reads follow the plain-image rule (ZBITPIX, BSCALE, BZERO); refused ones give `*TypeMismatchError` with its fields | `TestReadPixelsCompressedIntegerTargets` | cfitsio golden of every corpus file; exact range in math/big | `testdata/comp_fz_*.fits` | exact | 1 |
| Tile integers of another width than ZBITPIX, as cfitsio reads them: GZIP_1/GZIP_2 width from the inflated length (2/4/8 bytes on 8-bit, 1/4/8 on 16-bit, 1/2/8 on 32-bit, 1/4 on 64-bit), RICE_1 BYTEPIX 1, 2, else 4 (2/4 on 8-bit, 1/4 on 16-bit, 1/2 on 32-bit, 1/2/8 on quantized floats; ZVAL2 8 decodes as 4), quantized float GZIP tiles of 1 and 2 bytes, converted to the image type; the lowest and highest value of each ZBITPIX in wider tiles read, one past either refused (cfitsio status 412, funpack refuses); the null value cast to the tile type as cfitsio does (ZBLANK 406 is the unsigned byte 150 on 1-byte quantized tiles, BLANK 65541 is 5 on 1-byte tiles of a 32-bit image; a BLANK keyword outside the ZBITPIX range, e.g. 65536 on a 16-bit image, which cfitsio casts to 0 for 2-byte tiles and compares uncast with the 4-byte ints of HCOMPRESS_1/PLIO_1, is refused by the scaling-card rule, astropy reading it differently); real M13 (2006 SkyView) GZIP_1 4-byte tiles, RICE_1 without ZVAL2, HCOMPRESS_1 and PLIO_1; cfitsio's TDOUBLE read of 8-byte integer tiles returning its tile buffer unwritten (the harness's 0xA5 malloc fill), while its TLONG read is right, asserted | `TestCompressedCorpusRead`, `TestCompressedCorpusComplete`, `TestInflateTileWidths`, `TestCompressedInt64Read` | libcfitsio 4.7.0 decode (TLONG read where its TDOUBLE read of 8-byte tiles returns unwritten memory, and where funpack returns zeros: 8-byte tiles on 8-bit) + funpack twin; astropy 8.0.2 agrees | `testdata/comp_fz_real_m13_gzip.fits`, `testdata/comp_fz_real_m13_rice.fits`, `testdata/comp_fz_real_m13_hcomp.fits`, `testdata/comp_fz_real_m13_plio.fits`, `testdata/comp_fz_u8_gzip1_w4.fits`, `testdata/comp_fz_i16blank_gzip1_w8.fits`, `testdata/comp_fz_u8_rice_bp4.fits`, `testdata/comp_fz_u8_gzip2_w8.fits`, `testdata/comp_fz_i16small_rice_bp1.fits`, `testdata/comp_fz_i32small_gzip2_w2.fits`, `testdata/comp_fz_i16small_rice_blankcast.fits`, `testdata/comp_fz_f32_rice_bp1.fits`, `testdata/comp_fz_f64_rice_bp8.fits`, `testdata/comp_fz_f32_rice_bp1_zblank406.fits`, `testdata/comp_fz_i16small_hcomp_blank65536.fits`, `testdata/comp_fz_i16small_plio_blank65536.fits`, `testdata/comp_fz_f32_gzip1_q_w1.fits`, `testdata/comp_fz_f32_gzip2_q_w2.fits`, `testdata/comp_fz_f64_gzip1_q_w2.fits`, `testdata/comp_fz_u8_gzip2_w2_edge.fits`, `testdata/comp_fz_u16_gzip1_w4_edge.fits`, `testdata/comp_fz_i32_gzip1_w8_edge.fits` | bit-exact | 1 |
| ZBITPIX 64 (astropy writes it with GZIP_1/GZIP_2/NOCOMPRESS; fpack cannot): int64 over the full range and uint64 (BZERO 2^63), 8-byte and narrower GZIP tiles (1, 2, 4 bytes), read exactly as integers; BLANK in 8-, 2- and 1-byte tiles (masked and unmasked reads); RICE_1 (BYTEPIX 4), HCOMPRESS_1 and PLIO_1 at ZBITPIX 64, built from libcfitsio's own coders (`gen_fpack mk64`), which cfitsio decodes into C ints (both int32 extremes read) | `TestCompressedInt64Read`, `TestCompressedFixtures` | libcfitsio 4.7.0 TLONG `fits_read_img` (int64), its stored values + 2^63 (uint64: its TLONG read overflows, pinned), each equal to astropy 8.0.2's decode (NaN where cfitsio flags a BLANK); NOCOMPRESS: cfitsio's read of astropy's raw twin | `testdata/comp_fz_i64_rice.fits`, `testdata/comp_fz_i64_hcomp.fits`, `testdata/comp_fz_i64_plio.fits`, `testdata/comp_fz_i64small_gzip1_blank.fits`, `testdata/comp_fz_i64small_gzip1_w2_blank.fits`, `testdata/comp_fz_i64_gzip1.fits`, `testdata/comp_fz_u64_gzip2.fits`, `testdata/comp_fz_i64small_gzip1_w1.fits`, `testdata/comp_nocompress_i64.fits` | exact (integers) | 1 |
| Bad tiles and layouts refused, never decoded to wrong values: a GZIP tile of any other inflated length incl. leftover bytes (cfitsio status 414), tile values the ZBITPIX type cannot hold (cfitsio status 412; one past each edge at 8/16/32 bits), a null tile value the ZBITPIX type cannot hold in an unmasked read (a per-tile ZBLANK column, `i16_gzip1_w4_nullout2`: the masked reads and ReadFloat32 flag it, as cfitsio does; a BLANK keyword outside ZBITPIX, `i16_gzip1_w4_nullout`, is refused on every read: cfitsio flags it, astropy does not), BLANK / ZBLANK keyword or ZBLANK column value outside int32, with one past each end (2^31, -2^31 - 1) refused and both ends themselves read as null values (cfitsio ignores the keyword or fails 412, astropy refuses), a ZBLANK column winning over an out-of-range BLANK keyword, the first of two out-of-range null values named, quantized GZIP tiles of 8 bytes (cfitsio returns unwritten memory, funpack and astropy zeros), ZNAXIS outside 1..6 (cfitsio BAD_NAXIS; no panic in `Shape`, fitsdump reports it), a ZBITPIX that is no BITPIX and a negative ZNAXISn (both panicked; cfitsio reads the first and refuses the second with 414), a header that does not parse (`Header()` no longer panics), RICE_1 BYTEPIX 8 in the decoder (cfitsio has no 64-bit Rice coder) | `TestCompressedBadSamplesRefused`, `TestCompressedBadTilesFunpackRefuses`, `TestCompressedNullOutsideZBITPIX`, `TestCompressedBlankRangeAstropy`, `TestGZIPDecodersRefuseWrongLength`, `TestRiceDecoderRefusesBytepix8`, `TestCompressedReadRejects`, `TestDumpReportsCompressedMetadataError`, `TestCompressedHeaderParseErrorNoPanic`, `TestCompressedCorpusRead` | cfitsio goldens of the `comp_bad` samples; funpack 1.7.0 exit status; cfitsio imcompress.c | `testdata/comp_bad/u16_gzip1_extrabyte.fits`, `testdata/comp_bad/i16_gzip1_w4_overflow.fits`, `testdata/comp_bad/u16_rice_znaxis0.fits`, `testdata/comp_bad/u16_rice_znaxism1.fits`, `testdata/comp_bad/gzip1_7d.fits`, `testdata/comp_bad/u16_gzip1_w4_hi.fits`, `testdata/comp_bad/i16_gzip1_w4_nullout.fits`, `testdata/comp_bad/i32_rice_blank2p32.fits`, `testdata/comp_bad/i64_gzip1_blank2p40.fits`, `testdata/comp_bad/i16_rice_zblankcol1k.fits`, `testdata/comp_bad/f32_gzip1_q_w8.fits`, `testdata/comp_bad/i32_gzip1_blank2p31.fits`, `testdata/comp_bad/i32_gzip1_blankm2p31m1.fits`, `testdata/comp_bad/i32_gzip1_zblankcol2p31.fits`, `testdata/comp_bad/i32_gzip1_zblankcolm2p31m1.fits`, `testdata/comp_fz_i32_gzip1_blankmin.fits`, `testdata/comp_fz_i32_gzip1_blankmax.fits`, `testdata/comp_fz_i32_gzip1_zblankcolmin.fits`, `testdata/comp_fz_i16_zblankcol_blank2p32.fits`, `testdata/comp_bad/i16_gzip1_w4_nullout2.fits`, `testdata/comp_bad/u16_rice_zbitpix12.fits`, `testdata/comp_bad/u16_rice_zbitpix0.fits`, `testdata/comp_bad/u16_rice_znaxisneg.fits` | exact | 1 + 2 live |
| HCOMPRESS_1 and PLIO_1 integer tiles decode to C ints (null value compared uncast; BLANK 0 / 3 flag the zeros / 3s). Lossy HCOMPRESS_1 pixels reconstructed outside the ZBITPIX range: clipped for 16-bit images as funpack and cfitsio's TSHORT read do (cfitsio's TDOUBLE read returns them unclipped: pinned; its TUSHORT read of the uint16 sample fails with status 412, not pinned), refused for 8-bit (funpack 412), and refused at every ZBITPIX outside int32 (cfitsio's `(int)` cast wraps a saturated 32-bit pixel to about -2^31, also funpack; one past either int32 end at ZBITPIX 64). PLIO_1 values the ZBITPIX type cannot hold refused (cfitsio's read as the image type and funpack 412), also after the uint16 shift | `TestCompressedCorpusRead`, `TestCompressedBadSamplesRefused`, `TestCompressedBadTilesFunpackRefuses`, `TestReadFloat32NullsMatchCfitsioCorpus` | funpack twin (clip; cfitsio's TDOUBLE golden equal to it off the clipped pixels and outside the range on them); cfitsio goldens of the `comp_bad` samples; funpack exit status | `testdata/comp_fz_u16sat_hcomp_s2.fits`, `testdata/comp_fz_i16sat_hcomp_s2.fits`, `testdata/comp_fz_i16small_hcomp_blank0.fits`, `testdata/comp_fz_i16mask_plio_blank3.fits`, `testdata/comp_bad/u8sat_hcomp_s2.fits`, `testdata/comp_bad/i32sat_hcomp_s2.fits`, `testdata/comp_bad/i64_hcomp_2p31.fits`, `testdata/comp_bad/i64_hcomp_m2p31m1.fits`, `testdata/comp_bad/u16_plio_bzero0.fits`, `testdata/comp_bad/i32_plio_zbitpix8.fits`, `testdata/comp_bad/i32_plio_zbitpix16_bzero32768.fits` | exact | 1 + 2 live |
| ZCMPTYPE compared case-sensitively like cfitsio: `RICE_ONE` accepted, `rice_1` refused (`ErrUnknownCompression`) | `TestZCMPTYPEComparedLikeCfitsio`, `TestCompressedCorpusRead`, `TestParseAlgorithm` | cfitsio golden (refuses the HDU); funpack 1.7.0 (CFITSIO 4.7.0) status 414 | `testdata/comp_bad/u16_rice_lowercase.fits`, `testdata/comp_bad/golden/u16_rice_lowercase.fits.json`, `testdata/comp_fz_u16_rice_riceone.fits` | exact | 1 |
| Quantized floats: SUBTRACTIVE_DITHER_1/2 (exact zeros kept, -32 and -64), NO_DITHER (-32 and -64), absolute q, lossless ZQUANTIZ NONE, ZSCALE/ZZERO keywords, 3-D tiles; the dither random-sequence wrap inside one 10200-pixel -64 tile; ZQUANTIZ 'NONE' with ZSCALE/ZZERO columns at -32 and -64 (cfitsio NO_QUANTIZE: 4-byte tile integers taken as IEEE float32, widened for -64, and scaled by the tile's ZSCALE/ZZERO; astropy 8.0.2 dequantizes instead, the library follows cfitsio and funpack) | `TestCompressedCorpusRead`, `TestDequantizeSubtractive1`, `TestDequantizeSubtractive2PreservesZero` | libcfitsio golden + funpack twin | `testdata/comp_fz_f32_rice_sd2.fits`, `testdata/comp_fz_f64_rice_sd2.fits`, `testdata/comp_fz_f64_rice_nodith.fits`, `testdata/comp_fz_cubef32_rice_t3.fits`, `testdata/comp_fz_f64_rice_wrap.fits`, `testdata/comp_fz_f32_rice_zqnone.fits`, `testdata/comp_fz_f32_rice_zqnone_scaled.fits`, `testdata/comp_fz_f64_rice_zqnone.fits` | bit-exact | 1 |
| Masked read: ZBLANK keyword over BLANK, BLANK, per-tile ZBLANK column, NaN, including NaN inside GZIP_COMPRESSED_DATA and UNCOMPRESSED_DATA fallback tiles at -32 and -64; the unmasked reads return a null integer pixel as its stored value, scaled | `TestCompressedCorpusRead`, `TestCompressedUnmaskedNullsAreStoredValues` | cfitsio null flags (`fits_read_imgnull`); funpack twin where `fits_read_img` skips the null check (a NaN in a tile of an image without ZBLANK) | `testdata/comp_fz_i16_zblankcol.fits`, `testdata/comp_fz_i16blank_zblankkw.fits`, `testdata/comp_fz_f64flatnan_rice.fits`, `testdata/comp_fz_f32flatnan_uncomp.fits`, `testdata/comp_fz_f64flatnan_uncomp.fits` | exact | 1 |
| RICE BLOCKSIZE/BYTEPIX and HCOMPRESS SMOOTH taken from ZVAL1/ZVAL2 by position (NOISEBIT rule: ZNAME2 'NOISEBIT' with ZVAL2 = 1 on 4-byte Rice tiles reads as BYTEPIX 4, equal to the file it was made from; swap rule) | `TestCompressedCorpusRead`, `TestCompressedHDUHCompressMetadata` | libcfitsio golden + funpack twin | `testdata/comp_fz_u16_rice_znameswap.fits`, `testdata/comp_fz_f32_rice_noisebit.fits`, `testdata/comp_fz_u16_rice_zvalswap.fits` | bit-exact | 1 |
| Fallback tiles GZIP_COMPRESSED_DATA / UNCOMPRESSED_DATA in every column type; in scaled images we follow funpack/astropy | `TestCompressedCorpusRead` | funpack 1.7.0 (CFITSIO 4.7.0) twin (bit-exact) and astropy where `fits_read_img` differs (rawRows) | `testdata/comp_fz_f32flat_uncomp_scaled.fits`, `testdata/comp_fz_i16blank_uncomp.fits` | bit-exact | 1 |
| Invalid / unproven layouts refused loudly (unknown ZQUANTIZ on quantized; RICE_1 tiles of an unquantized float image other than 4-byte; HCOMPRESS_1/PLIO_1 tiles of an unquantized float image; one ZSCALE column, ZSCALE on integers; ZTILEn = 0, no COMPRESSED_DATA column, PLIO_1 with BZERO 32768 at ZBITPIX 8/32, a tile with no payload (cfitsio 415), an integer tile in GZIP_COMPRESSED_DATA (cfitsio 414), more table rows than tiles (cfitsio refuses the HDU), an UNCOMPRESSED_DATA tile of the wrong length (cfitsio reads it silently) on committed samples); tiled TABLE compression (ZTABLE = T) refused with `ErrTiledTableCompression` | `TestCompressedReadRejects`, `TestCompressedBadSamplesRefused`, `TestTiledTableCompressionRefused` | cfitsio behaviour (imcompress.c); cfitsio goldens of the `comp_bad` samples (cfitsio will not move to the ZTILE/COMPRESSED_DATA ones) | `testdata/comp_bad/*.fits`, `testdata/comp_bad/golden` | exact | 1 |
| Context cancellation before the read (before the metadata is parsed) and between tiles | `TestReadPixelsCompressedContextCanceled`, `TestReadPixelsCompressedContextCheckedBeforeRead`, `TestReadPixelsCompressedContextCanceledBetweenTiles` | `context.Canceled` | — | exact | 1 |
| astropy fixtures decode to the cfitsio-decoded truth; NOCOMPRESS at ZBITPIX 8/16/32/64 and the quantized -32/-64 astropy writes | `TestCompressedFixtures`, `TestCompressedFixtureSetComplete`, `TestReadFloat32_CompressedImage` | cfitsio golden of each astropy fixture; NOCOMPRESS: astropy 8.0.2 (input, or its own decode for floats; cfitsio cannot read NOCOMPRESS) | `testdata/comp_*.fits`, `testdata/comp_raw_*.fits` | bit-exact | 1 |
| Codec unit behaviour | `TestNOCOMPRESSIdentity`, `TestGZIP1Roundtrip`, `TestGZIP2Roundtrip`, `TestPLIO1Tile0`, `TestSelectHCompress`, `TestParseAlgorithm` | definitional / cfitsio streams | — | exact | 1 |

### Tile compression: write

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| Writer output is the same compressed image as `fits_img_compress` (96 cases: 87 identical incl. 4-D and 6-D, the quantizer range margins `rice_f32_margin_64x2` / `rice_f32_center_64` / `rice_f32_halfmargin_64`, a real -1e30 beside a NaN `rice_f32_nan_realnull_t40x5` / `rice_f32_nan_realnull_rows`, the dither-table re-seed in a tile with nulls `rice_f32_wrap_nan_z9000`; 1 identical except the documented FLOATNULLVALUE pixel; 5 cfitsio refusals refused incl. 7 axes; 3 deliberate refusals incl. a caller BLANK on float input); plus golden/img2: SUBTRACTIVE_DITHER_2 with NaNs and exact zeros in the same tiles (also BITPIX -64, `img64/rice_f64_sd2_nan_zeros_t16x8`) and a caller SIMPLE = F written as ZSIMPLE = F, as cfitsio writes it | `TestCompressedWriterMatchesCfitsio`, `TestCompressedWriterSimpleFMatchesCfitsio`, `TestReaderDecodesCfitsioImg2Goldens`, `TestHcompressDefaultTileMatchesCfitsio` | libcfitsio 4.7.0 `fits_img_compress` | `compress/testdata/cref/golden/img`, `compress/testdata/cref/golden/img2`, `compress/testdata/cref/golden/MANIFEST` | RICE/HCOMPRESS/PLIO byte-identical; GZIP inflated-identical (DEFLATE differs from zlib); ZSCALE/ZZERO bit-equal; keywords equal | 1 |
| Every written variant passes fitsverify and decodes identically in cfitsio, funpack, astropy (IEEE specials, BlockSize 16, oversized tiles, 4-px HCOMPRESS tiles, caller Z-keys, a real -1e30 in tiles with NaN). cfitsio's decode is `fits_read_img` with a NaN null value: its `fits_read_imgnull` flags real pixels null once several GZIP_COMPRESSED_DATA tiles hold NaNs, also on cfitsio's own file (327 stale flags in the golden; pinned) | `TestCompressedWriterInterop`, `TestFunpackMisdecodesNullTiledFloats`, `TestPLIOUint16AstropyMisreadsCfitsio`, `TestCfitsioImgnullStaleFallbackFlags` | fitsverify, cfitsio, funpack 1.7.0 (CFITSIO 4.7.0), astropy 8.0.2 | `compress/testdata/cref/golden/img/rice_f32_nan_realnull_rows.fits.gz` | exact; quantized floats within the half-step bound | 1 + 2 live |
| A real pixel equal to FLOATNULLVALUE (-9.11912e-36) or to the default NullValue (-1e30) in a tile with a NaN stays a value. In a GZIP_COMPRESSED_DATA (lossless) tile our reader, libcfitsio `fits_read_img` with a NaN null value (cref `decode`) and astropy return it exactly; funpack returns NaN for the FLOATNULLVALUE pixel (its output uses FLOATNULLVALUE as the null; pinned); `fits_img_compress` nulls that pixel (golden, the one byte range allowed to differ). In QUANTIZED tiles the tiles, ZSCALE and ZZERO equal cfitsio's when cfitsio writes the image with `fits_write_imgnull` and null value -1e30 (which quantizes the pixel as a value), and cfitsio's decode holds it as a value | `TestFloatNullValuePixelStaysValue`, `TestCompressedWriterMatchesCfitsio`, `TestCompressedWriterFloatNullValueQuantizedMatchesCfitsio` | libcfitsio 4.7.0 (`fits_img_compress`, `fits_write_imgnull`), astropy 8.0.2, funpack 1.7.0 | `compress/testdata/cref/golden/img/rice_f32_floatnullvalue_t40x5.fits.gz`, `compress/testdata/cref/golden/img/rice_f32_nan_realnull_t40x5.fits.gz`, `compress/testdata/cref/golden/imgw/rice_f32_fnv_quant_t40x5.fits.gz`, `compress/testdata/cref/golden/imgw/gzip2_f32_fnv_quant_nodither.fits.gz` | exact (lossless tile); quantized tiles byte-identical (RICE) / inflated-identical (GZIP), ZSCALE/ZZERO bit-equal | 1 + 2 live |
| A caller null value (Nullcheck with NullValue / NullValue64) is the null value of cfitsio's `fits_write_imgnull`: pixels equal to it, compared raw (a denormal NullValue nulls its pixels, another denormal stays a value), are nulls; the tiles, ZSCALE and ZZERO equal cfitsio's for -999 (float32, float64) and 1.04488e-40, and cfitsio's decode of its file has a null exactly at the NaNs and the null-value pixels. A null value of the other width (NullValue to the float64 writer, NullValue64 to the float32 writer) is refused instead of being ignored | `TestCompressedWriterNullValueMatchesCfitsio`, `TestCompressFloatNullValueOfOtherWidthRefused`, `TestCompressedWriterInterop` | libcfitsio 4.7.0 `fits_write_imgnull` into a compressing HDU | `compress/testdata/cref/golden/imgw/rice_f32_nullval_m999_t40x5.fits.gz`, `compress/testdata/cref/golden/imgw/rice_f32_nullval_denorm_t40x5.fits.gz`, `compress/testdata/cref/golden/imgw/rice_f64_nullval_m999_t40x5.fits.gz` | quantized tiles byte-identical, ZSCALE/ZZERO bit-equal, null positions exact | 1 + 2 live |
| Where cfitsio's quantizer overflows and writes wrong pixels with status 0 -- the int64 iqfactor for a tiny absolute QLevel (a tile of 100 at -1e-20 decodes as -0.0922, 1e10 at -1e-12 as -9223372.002), an int32 pixel when the centered zero point's float sum rounds low (14294966272 decodes as 9999998976), a centered minimum quantized onto NULL_VALUE (decodes as a null) -- the writer stores the tile losslessly and it reads back bit for bit in our reader, cfitsio, funpack and astropy; cfitsio's own wrong decode is pinned so the deliberate difference stays visible | `TestQuantizerOverflowStoredLosslessly`, `TestCompressedWriterInterop`, `TestReaderDecodesCfitsioImg2Goldens` | the input itself (lossless storage); libcfitsio 4.7.0 decode of its own file (the defect); fitsverify, funpack 1.7.0, astropy 8.0.2 | `compress/testdata/cref/golden/img2/rice_f32_zpt_overflow_abs.json`, `compress/testdata/cref/golden/img2/rice_f32_center_overflow.json`, `compress/testdata/cref/golden/img2/rice_f32_null_sentinel_hit.json`, `compress/testdata/cref/golden/img64/rice_f64_zpt_overflow_abs.json` | exact | 1 + 2 live |
| A caller BLANK on a float32 image is refused: cfitsio copies it into the ZBITPIX -32 header (no ZBLANK) and then cfitsio and astropy null the real pixels whose quantized integer equals it (64 of 1000 in the golden) | `TestWriteCompressedRefusals`, `TestFloatBLANKNullsRealPixels` | FITS 4.0 sect. 4.4.2.5; libcfitsio 4.7.0 decode of its own file; astropy 8.0.2 | `compress/testdata/cref/golden/img/rice_f32_blank_hdr.fits.gz`, `compress/testdata/cref/golden/img/rice_f32_blank_hdr.src.fits.gz` | exact (same null positions) | 1 + 2 live |
| RICE_1 BlockSize 16 is coded and declared: a one-tile image stores exactly `fits_rcomp*`'s blocksize-16 stream, ZVAL1 = 16, and reads back; a keyword is never written twice (caller EXTEND+ZEXTEND, CHECKSUM+ZHECKSUM, BLOCKED+ZBLOCKED, stale compression keys) | `TestCompressedWriterBlockSize16MatchesFitsRcomp`, `TestCompressedWriterNoDuplicateKeys` | `fits_rcomp`, `fits_rcomp_short`, `fits_rcomp_byte` (cref); fitsverify's duplicate-keyword rule (also run live in `TestCompressedWriterInterop`) | `compress/testdata/cref/golden/rice/rice_s_bs16.rice` | byte-exact | 1 + 2 live |
| RICE_1, HCOMPRESS_1, PLIO_1 encoders byte-identical to cfitsio (HCOMPRESS incl. quadrants with one side a power of two and the other larger: 8x10, 10x8, 16x18, 18x16, 8x9, 9x8); exact-buffer contract | `TestRiceEncodeMatchesCfitsio`, `TestHcompressEncodeMatchesCfitsio`, `TestPLIOEncodeMatchesCfitsio`, `TestPLIOEncodeRejectsOutOfRange`, `TestEncodersRefuseBytepix8`, `TestGzipEncodersBufferBoundary`, `TestEncodeGoldensComplete` | `fits_rcomp*`, `fits_hcompress[64]`, `pl_p2li` via cref | `compress/testdata/cref/golden/rice`, `compress/testdata/cref/golden/hcomp`, `compress/testdata/cref/golden/plio` | byte-exact | 1 |
| Float quantizer and noise estimator bit-exact (incl. cfitsio's float arithmetic for maxval - minval and minval + maxval, both N_RESERVED_VALUES range margins reached with a positive qlevel (the upper one also at 4294967286.75 steps, `margin_far_64x1`, whose centered extremes are clear of the sentinels), the centered zero point, qlevel 0, FnNoise1/FnNoise3 short rows; FnNoise5 rows with nvals 0, 1 with nvals2 0, 2 with nvals2 1; a FnNoise1 row with one difference and the 5-sigma clip on an exact tie; FnMeanSigma with one good pixel; an absolute qlevel with nulls; all-null tiles, the ngood == 0 dummy stdev, for qlevel 4/16/-0.5 and every dither). FnNoise1's `nx < 3` return and FnNoise3's second `nx < 5` test are equivalent mutants (a row of fewer than 2 differences is skipped and the median of no rows is 0; 4 flattened pixels give the same min/max/ngood through the row loop): the outputs for those shapes are pinned by `tiny_2x8` / `tiny_2x2_abs`; dither table | `TestQuantizeFloat32_CRef`, `TestImgStatsFloat32_CRef`, `TestCrefGoldensComplete`, `TestCrefRealNaNsInputHasNaNs`, `TestFitsRandomsMatchCfitsio` | `fits_quantize_float`, `fits_img_stats_float`, `fits_rand_value` | `compress/testdata/cref/golden/*.out.json` (quantized integers as an exact digest), `compress/testdata/cref/golden/fits_rand_value.json` | bit-exact | 1 |
| Refusals for their stated reason (int64/uint64, float32/float64 through the integer writer, float HCOMPRESS/PLIO, NOCOMPRESS, PLIO range/uint16, HCOMPRESS tile < 4, unknown dither, a caller BLANK on float input, a caller SIMPLE that is not logical, a QLevel of NaN or -Inf: no quantization step, every pixel would read back NaN; +Inf and tiny negative QLevels stay accepted: the tile's range or its zero point does not fit the integers and the tile is stored losslessly, `TestQuantizerOverflowStoredLosslessly`; a null value of the other width); shape refusals before anything is written (data length with `ErrShapeMismatch`, empty shape, non-positive axis, more than 6 axes, tile rank, non-positive tile axis) | `TestWriteCompressedRefusals`, `TestCompressedWriterShapeRefusals`, `TestCompressedFloatQLevelNaNInfRefused`, `TestCompressFloatNullValueOfOtherWidthRefused` | cfitsio limits cited per case; `fits_img_compress` status 212 for 7 axes (`img/rice_i16_7d`) | `compress/testdata/cref/golden/img/rice_i16_7d.json` | exact | 1 |
| float64 (ZBITPIX -64) writer `AppendCompressedFloat64Image` is the same compressed image as `fits_img_compress` (`fits_quantize_double`) on 19 BITPIX -64 images (the 20th, `rice_f64_zpt_overflow_abs`, is a cfitsio overflow stored losslessly instead; RICE_1/GZIP_1/GZIP_2, SUBTRACTIVE_DITHER_1/2, NO_DITHER, absolute qlevel, default qlevel, tiles 16x8 and 100x100 with ZDITHER0 9000, a cube, constant and huge rows (GZIP_COMPRESSED_DATA float64 fallback), IEEE specials incl. float64 denormals, NaN beside a real -1e30, a 1e8 background with 0.01 noise that a float32 quantizer cannot reproduce, every pixel above FLT_MAX, SUBTRACTIVE_DITHER_2 with NaNs and zeros; one documented difference, a real DOUBLENULLVALUE kept as a value); our float64 read of our file equals cfitsio's TDOUBLE decode of its own; fitsverify, cfitsio, funpack and astropy read 11 written variants (incl. Nullcheck with NullValue64, a caller header) within half a quantization step | `TestCompressedWriterFloat64MatchesCfitsio`, `TestCompressedWriterInterop` | libcfitsio 4.7.0 `fits_img_compress` / `fits_read_img` TDOUBLE, fitsverify, funpack 1.7.0, astropy 8.0.2 | `compress/testdata/cref/golden/img64` | RICE byte-identical, GZIP inflated-identical, ZSCALE/ZZERO bit-equal, decode bit-equal; interop: half a quantization step + 1 ulp | 1 + 2 live |
| Self-consistency only, not an independent proof: our writer then our reader, the half-step bound taken from the ZSCALE our writer wrote (ZSCALE itself is pinned bit-equal to cfitsio by `TestCompressedWriterMatchesCfitsio`) | `TestWriteCompressedFloat32_Roundtrip`, `TestWriteCompressedFloat32_GzipFallback`, `TestWriteCompressedRICE` | none (round trip) | — | half quantization step | 1 |

### WCS (`fits/wcs`, `fits/wcs/transform`)

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| All 28 projection codes (every one a degree-unit port of wcslib 8.10 prj.c; TSC and QSC were radian ports until round 2c), both directions, accept/reject boundaries pinned to 1e-7 deg, small-angle switches pinned: PCO's 1e-4 deg switches to 0.1% either side and at the threshold (`<` vs `<=`), and PCO's y = -0.0 at theta = 0 bit for bit; next to the native pole SIN forward (1e-5 rad, to 0.05%) and slant inverse (r^2 < 1e-10), SZP (r^2 < 1e-10), AIR (xi < 1e-4 rad, both directions) and XPH (sigma < 1e-4, both directions) to 0.1%, all at 4 ULP of wcslib plus the row's own libm sensitivity (every one of those switch mutants still fails), and the inverse SIN/SZP/AIR switches at an input exactly on the threshold; QSC's forward zeco < 1e-8 switch on all six faces to 0.05% (4 ULP plus libm sensitivity; `<=` for `<` survives, as before: no input ties); the final prjbchk of AZP (mu = 1, tilted: theta a few 1e-13 deg below -90 at its horizon) and XPH (phi past +-180 on its north-pole seam) pinned on both sides of the tolerance, and proven unreachable for TSC/CSC/QSC/AIR; QSC's inverse zeta range and CSC's float32-underflow guards proven equivalent (mutation sweep: every such threshold moved by 0.1%, 2%, 20%, 2x and 10x either way, every prjbchk tolerance to 1e-6, and `<=` for `<` wherever a tie is possible, is killed; `TestSmallAngleSwitchTies` and `TestQSCZecoSwitchTie` prove the other four untieable) | `TestWrefProjections`, `TestWrefGoldensReproduce`, `TestPCOSmallAngleBand`, `TestPoleSmallAngleBand`, `TestQuadcubeFaceCentreBand`, `TestSmallAngleSwitchTies`, `TestQSCZecoSwitchTie`, `TestQSCInverseZetaRange`, `TestPrjbchkRangeByConstruction`, `TestCSCUnderflowGuardPremises`, `TestQuadcubeFaceBound`, `TestSelectRejectsWhatPrjsetRejects`, `TestAllPaperIICodesImplemented`, `TestProjectionProofManifest`, `TestLibmPerturbHook`, `TestMolSolveMatchesBisection` | wcslib 8.10 `prjset`/`prjx2s`/`prjs2x` (`wref.c`) | `wcs/transform/testdata/wref/golden/*.json`, `wcs/transform/testdata/wref/golden/*.bin`, `wcs/transform/testdata/wref/golden/azp_mu1_gamma45.json` | 1e-11 deg within 1000 deg, relative 1e-14 beyond (only next to a forward divergence); boundaries ±1e-7 deg (above wcslib's 1e-12 prjbchk tolerance). The port calls Go's math functions where wcslib calls glibc's, so a row may exceed the budget by at most its own libm sensitivity (the change of our output when every sin/cos/atan2/log/... result moves by 4 ULP, `libm_sens_test.go`): needed by 333 of 171178 wref rows, all next to a pole, horizon or divergence (AIR, AZP, HPX, MOL, QSC, STG, XPH; STG at theta = -89.99 is 7e-9 relative, 1 ULP of sin(theta) magnified by 1/(1 + sin theta)), the measured difference never above 0.33 of the sensitivity, capped at 1% of the rows. Accept/reject decisions identical except at boundary-ambiguous points, where our own decision flips under a nudge of up to 4 ULP of either input or the 4-ULP libm perturbation: 2 rows (AZP mu = 1, gamma = 45 at (+-45, -165), its horizon), capped at 0.01% | 1; 2: `TestWrefGoldensReproduce` |
| Accuracy of the libm calls under the projections: Go's math.Sin/Cos/Sincos, math.Atan2/Atan, and the cancellation-free libAsin/libAcos (math.Asin/Acos are off by up to 973 / 2.5e7 ULP next to +-1 and are not used); the premise of the value budgets above | `TestLibmAccuracy`, `TestAtan2Accuracy`, `TestLibmGeneratorReproduces` | mpmath 1.4.1 at 200 bits | `wcs/transform/testdata/wref/libm/trig.json`, `wcs/transform/testdata/wref/libm/trig.bin`, `wcs/transform/testdata/wref/libm/atan2.json`, `wcs/transform/testdata/wref/libm/atan2.bin` (exact doubles, `wcs/transform/testdata/wref/f64pack.py`) | sin/cos 2 ULP (2^-100 absolute next to their zeros; measured 2), atan2/atan 1 ULP (measured 1), libAsin/libAcos 1 ULP (measured 1). Not correctly rounded: the proof audit's correctly rounded versions were removed (perf/trig) because they made the hot paths 3x to 17x slower | 1; 2: `TestLibmGeneratorReproduces` |
| asind/acosd out-of-range tolerance (WCSTRIG_TOL = 1e-10), every branch | `TestWcstrigGuards`, `TestWrefGoldensReproduce` | wcslib 8.10 `asind`/`acosd` (`wref trig`) | `wcs/transform/testdata/wref/libm/wcstrig_guards.json` | bit-exact at the guards; 2 ULP inside (-1, 1) (1 ULP of asin/acos between Go and glibc becomes up to 2 after * R2D: acosd(0.5) is 59.999999999999993 here, 60.000000000000007 in wcslib) | 1; 2: `TestWrefGoldensReproduce` |
| XPH against its defining paper, incl. the precision next to both poles | `TestXPHPaperEquations`, `TestXPHPolarPrecision`, `TestXPHPaperGeneratorReproduces` | Calabretta & Lowe 2013 eqs. (2)-(17) in mpmath | `wcs/transform/testdata/wref/paper/xph_paper.json`, `wcs/transform/testdata/wref/paper/xph_paper.bin`, `wcs/transform/testdata/wref/gen_xph_paper.py` | ~5e-14 deg, plus the stated sigma-rounding bound next to the poles; measured polar worst: forward north 2.6e-11 deg at theta = 89.99 (bound 3e-11; 2.4e-11 while the port's sin was correctly rounded: Go's sin(89.99 deg) is 1 ULP off, magnified by 1 - sin theta), inverse ~1e-14, forward SOUTH 4.08e-7 deg at theta = -89.9999997 (wcslib keeps eq. (9) there; the port follows wcslib) | 1; 2: `TestXPHPaperGeneratorReproduces` |
| HPX against its defining paper, incl. the precision next to the poles | `TestHPXPaperEquations`, `TestHPXPaperGeneratorReproduces` | Calabretta & Roukema 2007 sec. 6 in mpmath, (H, K) = (4, 3), (5, 5) | `wcs/transform/testdata/wref/paper/hpx_paper.json`, `wcs/transform/testdata/wref/paper/hpx_paper.bin`, `wcs/transform/testdata/wref/gen_hpx_paper.py` | 5e-14 deg (measured 1.24e-14; 1e-11 until round 2c let a 1 + 1e-13 scaled sigma pass), in the polar zones on top of the double-rounding bound of sigma, measured worst forward 2.98e-7 deg (theta = 90 - 3e-7), inverse 4.67e-7 deg (sigma ~ 1e-8), as in wcslib 8.10 | 1; 2: `TestHPXPaperGeneratorReproduces` |
| Full header pipeline (wcspih → wcsset → p2s/s2p incl. rejection status) for every projection and real headers | `TestHdrGoldens`, `TestHdrGoldensReproduce`, `TestHdrProjectionCoverage`, `TestPrjCaseGeneratorReproduce`, `TestSampleGeneratorsReproduce` | wcslib 8.10 (`hdr2sky.c`) | `wcs/transform/testdata/wref/hdrcases/prj_*.case`, `wcs/transform/testdata/wref/hdrcases/sample_1904-66_*.case`, `wcs/transform/testdata/wref/golden/hdr` | p2s 1e-8 arcsec (measured up to 3.2e-9, `prj_szp_cm40_mu1p5`); s2p 1e-9 px default, 1e-8 px for the iterative SIP/TPV inverses (real headers, 0.03-0.05 arcsec/px: wcslib's own round-off), 5e-8 px for COP next to its divergence (3 cases), 1e-7 px absolute for the `sample_1904-66_*` sky grids (28 cases, \|pixel\| up to ~3e7); each override is a `TOL_PX` line with its reason in the case file. An s2p point may exceed its budget by at most its own libm sensitivity (as for the projections): one point, `sample_1904-66_cop` at (90.315, 89.21), pixel ~3.3e7, 3.75e-4 px (1.2e-11 relative) | 1; 2: `TestHdrGoldensReproduce`, `TestPrjCaseGeneratorReproduce`, `TestSampleGeneratorsReproduce` |
| Real SIP / TPV headers: HST ACS/WFC3, JWST NIRCam (with AP/BP), Spitzer IRAC cryo/warm, astrometry.net, ZTF, PTF, ATLAS (TPV and its SIP copy, each alone) | `TestHdrGoldens`, `TestRealHdrCasesIntact` | wcslib 8.10 on verbatim archive cards (sha256-pinned) | `wcs/transform/testdata/wref/hdrcases/real_*.case`, `wcs/transform/testdata/wref/realhdr/make_cases.py` | as above | 1 |
| Real GLON survey headers with LONPOLE = 0 (ROSAT AIT, MSX CAR; CRVAL2 = 0, where 0 is the default); an explicit non-default LONPOLE = 0 honoured (TAN, CRVAL2 = 40) and refused where invalid (GLON-CAR/AIT at CRVAL2 < 0, as wcslib); HPXcvt HPX/XPH images | `TestHdrGoldens`, `TestHdrRejectCases`, `TestHPXcvtPixelCentres`, `TestHdrDerivedParameters` | wcslib 8.10; HPXcvt + healpy 1.20.1 | `wcs/transform/testdata/wref/hdrcases/real_rosat_glon_ait_lp0.case`, `wcs/transform/testdata/wref/hdrcases/real_msx_glon_car_lp0.case`, `wcs/transform/testdata/wref/hdrcases/lonpole_tan_dec40_lp0.case`, `wcs/transform/testdata/wref/hdrcases_reject/lonpole_glon_car_bm30_lp0.case`, `wcs/transform/testdata/wref/hdrcases_reject/lonpole_glon_ait_bm20_lp0.case`, `wcs/transform/testdata/wref/samples/hpxcvt_hpx_gal_n8.fits` | as above; HPXcvt vs healpy 8.5e-7 deg (CDELT written to 8 digits) | 1 + 2 live |
| LONPOLE/LATPOLE, PVi_0..4 (incl. PV1_0 = 1 with only one of PV1_1/PV1_2: no offset), CROTA2 (PC > CD > CROTA2), CUNIT, PC+CD, partial CD, swapped axes, cubes, WCSAXESa, alternate WCS, latitude-first headers | `TestHdrGoldens`, `TestHdrDerivedParameters`, `TestParseNAxisRules`, `TestParsePCAndCDBothAccepted`, `TestParseAltDescription` | wcslib 8.10 goldens | `wcs/transform/testdata/wref/hdrcases/lonpole_*.case`, `wcs/transform/testdata/wref/hdrcases/crota2_*.case`, `wcs/transform/testdata/wref/hdrcases/cunit_*.case`, `wcs/transform/testdata/wref/hdrcases/alt_a_*.case`, `wcs/transform/testdata/wref/hdrcases/prjswap_*.case`, `wcs/transform/testdata/wref/hdrcases/pvfid_*.case`, `wcs/transform/testdata/wref/hdrcases/pv1off_*.case` | as above | 1 |
| AIPS NCP (→ SIN, PV2_2 = cot CRVAL2, given PV2_1/2 overwritten; refused on the equator) and GLS (→ SFL with the fiducial offset) translated as wcslib `wcsset` does | `TestHdrGoldens`, `TestHdrRejectCases` | wcslib 8.10 (`hdr2sky`, wcsset alias translation) | `wcs/transform/testdata/wref/hdrcases/alias_*.case`, `wcs/transform/testdata/wref/hdrcases_reject/ncp_equator.case` | as above | 1 + 2 live |
| Defunct PC00i00j / CD00i00j / PROJPn keywords read as PCi_j / CDi_j / PVi_n of the latitude axis (wcslib `WCSHDR_PC00i00j`, `CD00i00j`, `PROJPn`; astropy relax=True), incl. a real Pan-STARRS GPC1 skycell header; a defunct card that disagrees with its modern twin is refused | `TestHdrGoldens`, `TestRealHdrCasesIntact`, `TestHdrStricterThanWcslib` | wcslib 8.10 `wcspih` (WCSHDR_all in `hdr2sky`) | `wcs/transform/testdata/wref/hdrcases/real_ps1_pc00i00j.case`, `wcs/transform/testdata/wref/hdrcases/legacy_*.case`, `wcs/transform/testdata/wref/hdrcases_stricter/legacy_*.case` | as above | 1 + 2 live |
| SIP keywords without the `-SIP` CTYPE suffix are applied as wcslib applies them (wcspih makes SIP the prior distortion from A_p_q/B_p_q alone) on TAN, SIN and every projection but TPV, on the primary only (an alternate description stays undistorted); `RA---TPV-SIP` (the suffix written out) applies SIP then TPV, as wcslib does (it accepts the CTYPE; proven on a synthetic header and on the real ATLAS TPV solution with a non-redundant SIP set, whose answer moves off TPV alone by exactly the SIP vector, 2.33-3.0 px, to 1.1e-10 px). A `RA---TPV` header without `-SIP` ignores its SIP keywords (next row, a deliberate difference). On a TNX header (with or without `-SIP`) SIP keywords are refused, by `wcs.Parse` and so by `New` and `NewTNX` (wcslib applies SIP with TAN and ignores the WAT correction when the SIP cards come last, and fails when they come first) | `TestHdrGoldens`, `TestParseSIPKeysWithoutSuffix`, `TestNewAppliesSIP`, `TestTPVSIPSuffixAppliesBoth`, `TestHdrStricterThanWcslib`, `TestNewTNXRefusesSIP`, `TestHdrRequiredCases` | wcslib 8.10 (`hdr2sky`) | `wcs/transform/testdata/wref/hdrcases/sip_nosuffix_tan.case`, `wcs/transform/testdata/wref/hdrcases/sip_nosuffix_sin.case`, `wcs/transform/testdata/wref/hdrcases/tpv_sip_suffix.case`, `wcs/transform/testdata/wref/hdrcases/tpv_sip_suffix_atlas.case`, `wcs/transform/testdata/wref/hdrcases/sip_nosuffix_alt_a.case`, `wcs/transform/testdata/wref/hdrcases/sip_nosuffix_primary_of_alt.case`, `wcs/transform/testdata/wref/hdrcases/sip_swapped_nosuffix.case`, `wcs/transform/testdata/wref/hdrcases/sip_swapped_suffix.case` (DEC on pixel axis 1: A stays on pixel axis 1), `wcs/transform/testdata/wref/hdrcases/sip_cube_radec_freq.case` (RA, DEC, FREQ cube), `wcs/transform/testdata/wref/hdrcases_stricter/sip_a_without_b_nosuffix.case`, `wcs/transform/testdata/wref/hdrcases_stricter/sip_nosuffix_cube_axes23.case`, `wcs/transform/testdata/wref/hdrcases_stricter/tnx_sip_nosuffix.case`, `wcs/transform/testdata/wref/hdrcases_stricter/tnx_sip_suffix.case` | as above (measured <= 3.4e-10 arcsec, 4.7e-10 px) | 1 + 2 live |
| `RA---TPV` / `DEC--TPV` without `-SIP` plus SIP keywords (the ATLAS form: a SIP copy of the TPV correction for SIP-only readers): the SIP keywords are ignored (`wcs.Header.SIPIgnored`; not read, so never an error), as astropy 8.0.2 ignores them, a deliberate difference from wcslib 8.10, which applies both. Evidence, the real ATLAS header of PR #4 (Daryl Manning; stamp 01a59742o0240o.diff, 100x100, 12 x 12 grid incl. the corners, and the centre): wcslib on its TPV part alone and on its SIP part alone (as TAN-SIP) agree to 2.31e-7 arcsec (1.2e-7 px; corner offsets (-12.28, 6.06) px at (1, 1), -12.03..-13.14 / 5.63..6.33 px over the four corners); wcslib on the header as written lands 24.45-26.84 arcsec, 13.28-14.59 px, from the TPV answer at every point (the correction counted twice); astropy on the header as written equals wcslib on TPV alone to 3.6e-10 arcsec; ours equals astropy to 3.7e-10 arcsec and 6.9e-11 px. Also the real PTF TPV + SIP pair (wcslib 72.7-133.4 arcsec off) and a synthetic non-redundant SIP set (0-8.3 arcsec) | `TestHdrDeliberateDifferences`, `TestATLASSIPCopiesTPV`, `TestATLASTPVPlusRedundantSIP`, `TestParseTPVIgnoresSIPWithoutSuffix`, `TestRealHdrCasesIntact`, `TestHdrRequiredCases` | astropy 8.0.2 `all_pix2world` / `all_world2pix` (`astropy_hdr.py`, recorded); wcslib 8.10 (`hdr2sky`: the TPV-only goldens, and its recorded reading of each header as written) | `wcs/transform/testdata/wref/hdrcases_differ/real_atlas_k2_tpv_sip.case`, `wcs/transform/testdata/wref/hdrcases/real_atlas_k2_tpv.case`, `wcs/transform/testdata/wref/hdrcases/real_atlas_k2_tansip.case`, `wcs/transform/testdata/wref/hdrcases_differ/real_ptf_tpv_plus_sip.case`, `wcs/transform/testdata/wref/hdrcases/real_ptf_tpv4.case`, `wcs/transform/testdata/wref/hdrcases_differ/tpv_sip_nosuffix.case`, `wcs/transform/testdata/wref/hdrcases/tpv_sip_nosuffix_tpvonly.case`, `wcs/transform/testdata/wref/astropy_hdr.py`, `testdata/recorded/wcs/transform/TestHdrDeliberateDifferences.jsonl.gz`, `testdata/recorded/wcs/transform/TestATLASTPVPlusRedundantSIP.jsonl.gz` | ours vs the TPV-only wcslib golden: that golden's budgets (1e-8 arcsec, 1e-9 px); ours and wcslib-TPV vs astropy 1e-8 arcsec, 1e-9 px (wcslib 8.6 in astropy vs 8.10, measured <= 3.7e-10 arcsec, 5.4e-10 px); SIP part vs TPV part 2.4e-7 arcsec, corner offsets 5e-6 px (the SIP copy is a fit, measured 2.31e-7 arcsec); wcslib's doubling pinned 24.4..26.9 arcsec, 13.2..14.6 px | 1 + 2 live |
| DSS plate solutions (CNPIX, PPO3/6, X/YPIXELSZ, PLTRA*/PLTDEC*, AMDX/AMDY 1..13) applied as wcslib 8.10 does (except that a PLATEID card after the plate keywords does not switch the solution off, see Deliberate differences): the plate keywords replace CRPIX/PC/CDELT/CRVAL/CTYPE and erase CD (wcspih_final), the plate polynomial is a sequent distortion in wcslib's TPD form (degree 3, or 5 with AMD13; dssset), sky->pixel by Newton iteration (Jacobian pinned against central differences of the forward polynomial; 6000 px outside the plate, where wcslib's disx2p fails, the returned pixel is checked against the source pixel of wcslib's pixel->sky); refused: partial plate keywords (wcslib reads uninitialised memory), AMDm with m > 13 non-zero (wcslib crashes), PLTDECSN other than +/- or an unpadded '-', zero-padded AMD indices, DSS with SIP, DSS on 3 axes, zero plate scale | `TestHdrGoldens`, `TestHdrRejectCases`, `TestHdrStricterThanWcslib`, `TestHdrLenientThanWcslib`, `TestDSSJacobian`, `TestParseDSSRefusals`, `TestParseDSSTranslation`, `TestRealHdrCasesIntact`, `TestHdrRequiredCases` | wcslib 8.10 (`hdr2sky`); astropy 8.0.2 agrees on the probe | `wcs/transform/testdata/wref/hdrcases/real_dss_plate_solution.case`, `wcs/transform/testdata/wref/hdrcases/dss_modern_cutout.case`, `wcs/transform/testdata/wref/hdrcases/dss_all_terms_deg5.case`, `wcs/transform/testdata/wref/hdrcases/dss_far_outside_plate.case`, `wcs/transform/testdata/wref/hdrcases_lenient/dss_plateid_last.case`, `wcs/transform/testdata/wref/hdrcases_reject/dss_with_sip.case`, `wcs/transform/testdata/wref/hdrcases_reject/dss_zero_scale.case`, `wcs/transform/testdata/wref/hdrcases_stricter/dss_naxis3.case`, `wcs/transform/testdata/wref/hdrcases_stricter/dss_pltdecsn_unpadded.case`, `wcs/transform/testdata/wref/hdrcases_stricter/dss_amd_leading_zero.case` | p2s 1e-8 arcsec (measured <= 2.2e-10); s2p 1e-9 px (measured <= 7.8e-10 px), 3e-9 px for `dss_far_outside_plate` (wcslib's own round trip is 1.3e-9 px off there) | 1 + 2 live |
| Keyword rules pinned one case each: per-axis CUNIT (mixed units), PV1_0 = 1 with only one of PV1_1/PV1_2 (no offset, the given phi_0/theta_0 still sets the pole), theta_0 just above 90 clamped within celset's 1e-10 (90 + 5e-11) and refused just beyond it (90 + 1.5e-10), every celset tolerance bracketed (z = cos latp cos lat0 at 1.6e-11 and 1.75e-10, a nearer native-pole latitude 1e-9 past +-90 passed over on either solution, both solutions 1e-9 out refused, slz = 1 + 5e-11 clamped and 1 + 2e-10 refused), sky points exactly at and 1e-7 deg from the native poles of CAR (wcslib `sphs2x` longitude-origin convention and its rearranged |x| < 1e-5 formula, ported in `sphS2X`), the celestial south pole at the native pole (celset z = 0, latp < 0) and the LATPOLE tie between two solutions, celset's z = 0 branch (theta_0 = 0, \|phi_p - phi_0\| = 90: LATPOLE 95 / -95 clamped to +-90, CRVAL2 != 0 refused), LONPOLE = 999 read as the default; refusals of PV2_40 (TPV latitude axis), PV2_5 on a longitude axis that is axis 2, two latitude axes; stricter refusals of TPV with a non-degree latitude CUNIT, matrix coupling into the latitude (and longitude) row, a numeric PSi_m | `TestHdrGoldens`, `TestHdrDerivedParameters`, `TestHdrRejectCases`, `TestHdrStricterThanWcslib`, `TestHdrRequiredCases` | wcslib 8.10 (`hdr2sky`; acceptance for the stricter cases) | `wcs/transform/testdata/wref/hdrcases/cunit_mixed_deg_arcsec.case`, `wcs/transform/testdata/wref/hdrcases/cunit_mixed_arcmin_deg.case`, `wcs/transform/testdata/wref/hdrcases/pvlon_car_offset_phi0_only.case`, `wcs/transform/testdata/wref/hdrcases/pvlon_car_offset_theta0_only.case`, `wcs/transform/testdata/wref/hdrcases/pvlon_car_theta0_tol.case`, `wcs/transform/testdata/wref/hdrcases/celpole_car_eq_south.case`, `wcs/transform/testdata/wref/hdrcases/celpole_car_eq_latpole_tie.case`, `wcs/transform/testdata/wref/hdrcases/celpole_z_above_tol.case`, `wcs/transform/testdata/wref/hdrcases/celpole_z_within_tol.case`, `wcs/transform/testdata/wref/hdrcases/celpole_latp1_tol_out.case`, `wcs/transform/testdata/wref/hdrcases/celpole_latp2_tol_out.case`, `wcs/transform/testdata/wref/hdrcases/celpole_slz_within_tol.case`, `wcs/transform/testdata/wref/hdrcases/celpole_car_native_pole_sky.case`, `wcs/transform/testdata/wref/hdrcases/celpole_z0_latpole95.case`, `wcs/transform/testdata/wref/hdrcases/celpole_z0_latpolem95.case`, `wcs/transform/testdata/wref/hdrcases/lonpole_999_default.case`, `wcs/transform/testdata/wref/hdrcases_reject/celpole_z0_lat30.case`, `wcs/transform/testdata/wref/hdrcases_reject/pvlon_car_theta0_tol_out.case`, `wcs/transform/testdata/wref/hdrcases_reject/celpole_slz_tol_out.case`, `wcs/transform/testdata/wref/hdrcases_reject/celpole_latp_tol_out.case`, `wcs/transform/testdata/wref/hdrcases_reject/tpv_pv2_40.case`, `wcs/transform/testdata/wref/hdrcases_reject/pvlon_axis2_pv5.case`, `wcs/transform/testdata/wref/hdrcases_reject/pvlon_car_theta0_90p5.case`, `wcs/transform/testdata/wref/hdrcases_reject/two_lat_axes.case`, `wcs/transform/testdata/wref/hdrcases_stricter/tpv_cunit2_arcsec.case`, `wcs/transform/testdata/wref/hdrcases_stricter/couple_lat_row.case`, `wcs/transform/testdata/wref/hdrcases_stricter/couple_lng_row.case`, `wcs/transform/testdata/wref/hdrcases_stricter/ps_numeric.case` | as above | 1 + 2 live |
| Every hand-written header case is required by name (deleting one, or its golden, fails) | `TestHdrRequiredCases` | — | `wcs/transform/hdr_required_test.go` | exact | 1 |
| SLON/SLAT and HLON/HLAT pairs; `-SIP` on a projection other than TAN | `TestHdrGoldens`, `TestParseCelestialPairRules` | wcslib 8.10 | `wcs/transform/testdata/wref/hdrcases/frame_slon_tan.case`, `wcs/transform/testdata/wref/hdrcases/frame_hlon_zea.case`, `wcs/transform/testdata/wref/hdrcases/sip_on_sin.case` | as above | 1 |
| Headers wcslib refuses are refused | `TestHdrRejectCases` | wcslib 8.10 refusals | `wcs/transform/testdata/wref/hdrcases_reject/*.case` | exact | 1 + 2 live |
| Headers we refuse on purpose although wcslib accepts them (incl. a singular CD: wcslib 8.10 sets it up and computes pixel→sky, only sky→pixel has no answer) | `TestHdrStricterThanWcslib`, `TestSingularLinearMatrixRejected`, `TestNewRefusesUnimplementedDistortions`, `TestNewRefusesCoupledMatrix`, `TestNewRefusesSIPOffAxes12`, `TestNewRefusesTPVWithNonDegreeUnits`, `TestDistortionCodeMustMatchOnBothAxes`, `TestParseUnsupportedDistortionListed` | wcslib 8.10 acceptance (first test); policy for the rest | `wcs/transform/testdata/wref/hdrcases_stricter/*.case` | exact | 1 + 2 live |
| Headers wcslib refuses that we accept on purpose, each reproducing a named wcslib golden: a SIP header without AP/BP or with one of them empty (wcslib needs an inverse coefficient on each axis although it only seeds the iteration), a DSS header with PLATEID after the plate keywords (wcslib's lexer resets the distortion type there) | `TestHdrLenientThanWcslib`, `TestHdrRequiredCases` | wcslib 8.10 refusal (as written) + the SAME-AS golden | `wcs/transform/testdata/wref/hdrcases_lenient/*.case` | the SAME-AS golden's budgets | 1 + 2 live |
| SIP/TPV inverse: non-convergence is an error (incl. 2000 targets just past each fold, where a slow creep must not be accepted); a round-off plateau is accepted only below `newtonStallTol`; AP/BP only a starting guess; relative convergence at large offsets | `TestSIPInverseNonConvergenceNotOK`, `TestTPVInverseNonConvergenceNotOK`, `TestNewtonStallAtRoundOffAccepted`, `TestSIPInverseConvergesAtLargeOffsets`, `TestSIPAPInverseIsStartingGuessOnly`, `TestTPVMissingPV1DefaultsToOne`, `TestNewtonRefusesNonFiniteJacobian`, `TestTPVRadialJacobian` | analytic folds; exact quantised-evaluation model; wcslib `wcss2p`; TPV registry; overflowing-determinant / overflowing-step inputs (exact); closed-form d(c r^p) | `wcs/transform/testdata/wref/hdrcases/sip4_*.case`, `wcs/transform/testdata/wref/hdrcases/tpv_*.case` | as above | 1 |
| SIP of every accepted order: orders 7 and 9 with every A/B term set (no AP/BP); TPV with every PVi_0..PVi_39 non-zero on both axes (CD and PC + CDELT forms); ZPN of degree 12 and with PV2_29 | `TestHdrGoldens`, `TestHighOrderCaseGeneratorReproduce`, `TestWrefProjections`, `TestProjectionProofManifest` | wcslib 8.10 (`hdr2sky`, `wref`) | `wcs/transform/testdata/wref/hdrcases/hiorder_*.case`, `wcs/transform/testdata/wref/gen_highorder_cases.py`, `wcs/transform/testdata/wref/hdrcases/prj_zpn_c0_deg12.case`, `wcs/transform/testdata/wref/golden/zpn_pv29.json` | p2s 1e-8 arcsec (measured 1.8e-10); s2p 1e-8 px (measured 4.6e-10) | 1; 2: `TestHighOrderCaseGeneratorReproduce` |
| `wcs.Parse` rules: SIP orders 0..9, SIP coefficient cards of degree > 9 refused even when zero (wcslib "Invalid powers", `wcs/transform/testdata/wref/hdrcases_reject/sip_zero_coeff_deg10.case`), partial SIP, non-numeric A/B/AP/BP coefficients, SIP coefficient keywords whose index is not a single digit (`A_01_1`, `B_0_02`, `A_1_10`: wcslib 8.10 discards the card silently, so they are refused, `wcs/transform/testdata/wref/hdrcases_stricter/sip_index_leading_zero.case`), wrong-typed keywords, celestial pairs, defunct PC00i00j/CD00i00j/PROJPn; every SIP refusal stricter than wcslib (which never reads A_ORDER/B_ORDER) has a wcslib-accepts case (`sip_beyond_order_nonzero`, `sip_order10`, `sip_neg_order`, `sip_order_float`, `sip_coef_string`, `sip_mixed_suffix`, `sip_alt_suffix`, `sip_unknown_dist_code`), as does wcslib's own `TPU` alias (`tpu_alias`), which we do not implement | `TestParseSIPOrderValidation`, `TestParseWrongTypeIsAnError`, `TestParseCelestialPairRules`, `TestParseSIPSuffix`, `TestParseDefunctKeywords`, `TestHdrStricterThanWcslib`, `TestHdrRejectCases` | SIP convention; wcslib TPD limit; wcslib 8.10 acceptance / refusal; policy | — | exact | 1 + 2 live |
| astropy crosscheck through `wcs.Parse` (cards, both directions) | `TestCrossCheckFixtures`, `TestCrossCheckFixturesReproduce` | astropy 8.0.2 (wcslib 8.6) `all_pix2world` / `all_world2pix` | `testdata/wcs_crosscheck.json`, `testdata/generate_wcs_crosscheck.py` | per-fixture `tol_arcsec` / `tol_px` (1e-9 px) | 1; 2: `TestCrossCheckFixturesReproduce` |
| End-to-end file path (WriteImage → Open → ParseAlt → New) | `TestEndToEndFileAgainstWcslib`, `TestEndToEndFITSWCS` | wcslib goldens; TAN closed form | `wcs/transform/testdata/wref/golden/hdr` | the header-pipeline budgets above (1e-8 arcsec, 1e-9 px default, per-case `TOL_PX`); TAN closed form 1e-10 deg / 1e-7 px | 1 |
| Sky frames: ICRS/FK5 ↔ galactic/supergalactic, ICRS ↔ IAU 2006 ecliptic of date, FK5 precession (IAU 1976), FK4 B1950 ↔ FK5 J2000 | `TestFramesAstropyGolden`, `TestFramesGoldenComplete`, `TestFramesGoldenReproduce`, `TestGalacticDefiningPoints`, `TestSupergalacticDefiningPoints`, `TestObliquityJ2000`, `TestPrecessMeeus21b`, `TestVec3ToLonLatWrap`, `TestRoundTrips` | astropy 8.0.2 / ERFA (`prec76` for precession); Meeus ex. 21.b | `wcs/transform/testdata/frames/frames_golden.json` | 1e-8 arcsec (144 points per case, 31 cases required) | 1; 2: `TestFramesGoldenReproduce` |
| IRAF TNX: WAT joining (gaps, spellings, duplicate keys, CONTINUE, 1999-character limit refused; non-numeric WATn_ keys ignored), surface types 1/2/3, cross terms, LONPOLE precedence (card <= 360, negative, WAT, 180), lngcor/latcor cross-axis fallback, RA/DEC, GLON/GLAT, ELON/ELAT, SLON/SLAT, HLON/HLAT axes, no-WAT TNX = TAN (without PVi_m), latitude-axis PVi_m ignored when a correction surface exists, PROJPn ignored, linear step (CD; PC + CDELT pair with PC1_1/PC001001; PC00i00j; 1 deg/pixel with no scaling), CUNIT deg, strict validation, Newton inverse; refused where WCSTools and wcslib disagree or WCSTools ignores a card (longitude-axis PVi_m, latitude-axis PVi_m with no surface, CD with PC, CD00i00j without CDi_j, one of CDELT1/CDELT2 alone or PC without CDELT1 or PC1_1 when there is no CD, SIP keywords with or without -SIP, lowercase CTYPE, LONPOLE > 360, CUNIT other than deg) | `TestTNXWCSToolsGolden`, `TestTNXRefusedHeaders`, `TestTNXGoldensReproduce`, `TestTNXSyntheticHeadersReproduce`, `TestTNXRealHeaderDetails`, `TestTNXRejects`, `TestConcatWATRefusesWhatWCSToolsDrops`, `TestConcatWATKeepsChunkBlanks`, `TestWATCardStringDoubledQuote`, `TestTNXBasisDerivatives`, `TestTNXSurfaceDerivativesRealHeaders`, `TestTNXInverseNonConvergence`, `TestNewRejectsTNX` | WCSTools 3.9.7 `xy2sky`/`sky2xy`; closed-form polynomials | `wcs/transform/testdata/tnx/headers/kpno_kp4m_mosaic1_ccd1.hdr`, `wcs/transform/testdata/tnx/headers/noao_ct4m_tnx1.hdr`, `wcs/transform/testdata/tnx/headers/registry_tnxdoc_tnx3.hdr`, `wcs/transform/testdata/tnx/golden` (32 cases required by name and point count), `wcs/transform/testdata/tnx/refused` (17 refused cases with their WCSTools answers) | p2s 1e-8 arcsec (12-decimal xy2sky output); s2p 1e-7 px vs source, 1e-4 px vs sky2xy (WCSTools stops its Newton iteration at 2.8e-8 deg; measured worst 1.4e-5 px) | 1; 2: `TestTNXGoldensReproduce`, `TestTNXSyntheticHeadersReproduce` |
| Bisection solver reports non-convergence | `TestSolveBisectionNonConvergenceIsAnError` | analytic | — | exact | 1 |

### HEALPix (`fits/healpix`)

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| Ang2Pix / Vec2Pix incl. the near-pole branch switch-over points, z exactly 2/3 (equatorial/polar split, RING and NESTED, Nside 2^0..2^29) and directions within an ULP of a boundary; the arithmetic and every branch choice after the libm calls (`ang2pixCS`, `vec2pixPhi`) fed the correctly rounded cos/sin/atan2 the references were made with | `TestAng2Pix_Healpy`, `TestAng2PixRandom_Healpy`, `TestAng2PixPolar_Healpy`, `TestBranchSwitch_Healpy`, `TestAng2PixTwoThird_Healpy`, `TestVec2Pix_Healpy`, `TestAng2PixPhiWrap`, `TestBoundaryAmbiguityDetects` | healpy 1.20.1; real healpix_cxx 8.3.0 (`hpref`) fed mpmath-rounded trig | `healpix/testdata/golden/ang2pix_*.json`, `healpix/testdata/golden/branch_switch.json` | exact (pixel indices) on the correctly rounded libm values and on the standard and random direction sets. With Go's math.Cos/Sin/Atan2 a different pixel is accepted only at a boundary-ambiguous direction (the pixel changes under a nudge of up to 4 ULP of the inputs or of cos/sin/atan2, and the reference is among the pixels the nudges give), counted and capped: 1536 of 56000 polar-boundary pixels (2.7%, cap 6%), 237 of 1946 switch-over inputs (12.2%, cap 20%), 914 of 38400 z = 2/3 inputs (2.4%, cap 5%), all constructed on boundaries | 1 |
| Vec2Pix for any non-zero finite vector length (every vector rescaled by a power of 2; phi from the unscaled components; a rescale threshold of 2^-482 or lower fails, 2^-481 and above are equivalent); healpy itself moves none of the near-pole vectors at 2^-480 / 2^-481, the first at 2^-482 (Nside 2^29 - 1 RING), tens at 2^-490; healpy's (1e-200, 0, 1e-200) is pixel 0 (RING) or another wrong pixel (NESTED) at every Nside > 1 | `TestVec2PixExtremeLength`, `TestVec2PixScaledPolar`, `TestVec2PixScaledPolar_HealpyMoves`, `TestVec2PixTinyVector_Healpy` | healpy at unit scale (boundary-ambiguous vectors as above: 481 of 28000) + exact power-of-2 invariance (bit for bit against the unscaled pixel); the polar-boundary vectors scaled by every 2^k for k = -512..-470 and 2^-1060..2^1020; healpy's scaled answers at every k = -490..-480, -495, -500 | `healpix/testdata/golden/vec2pix_extreme.json`, `healpix/testdata/golden/ang2pix_polar.json`, `healpix/testdata/golden/vec2pix_scaled_healpy.json` | exact | 1 |
| RING any nside ≤ 2^29, NESTED power of 2; Pix2Ang/Pix2Vec, xyf, ordering, neighbours | `TestPixFuncs_Healpy`, `TestPix2Ang_Healpy`, `TestPix2Vec_Healpy`, `TestXYF_Healpy`, `TestOrdering_Healpy`, `TestOrderingSampled_Healpy`, `TestNeighbors_Healpy`, `TestNeighboursOrdered_Healpy`, `TestAngleULPBudget_Measured`, `TestRoundTrip_Exhaustive`, `TestRingDecomposition_ConstantTime` | healpy 1.20.1 | `healpix/testdata/golden/pixfuncs_*.json` | indices exact; angles 1 ULP; vectors 2·2^-52 | 1 |
| MaxPixRad, PixelResol, PixelArea, Npix2Nside | `TestUtilities_Healpy`, `TestUtilitiesExt_Healpy` | healpy; mpmath exact definition | `healpix/testdata/golden/utilities_ext.json` | MaxPixRad 2 ULP of exact / 2^-52 abs vs healpy (healpy's own cancellation); others bit-exact | 1 |
| Invalid input panics with the package's own message | `TestInvalidInputPanics` | documented contract | — | exact | 1 |
| HEALPix maps in FITS binary tables | `TestRealMapHeader_WMAP`, `TestHealpyWrittenMap` | NASA LAMBDA header; healpy `write_map` | `healpix/testdata/maps/wmap_ilc_9yr_v5.hdr`, `healpix/testdata/maps/healpy_map_nside16_ring.fits` | exact | 1 |

### ReadFloat32, frames and `image.Image` (root package)

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| ReadFloat32 normalisation (float32 nearest the exact value; BITPIX 32/64 double-rounding cases; negative BSCALE; float BITPIX with BZERO and no BSCALE; BSCALE=0 and overflow errors; uint64 with BZERO 2^63 as a real-valued and as an integer-literal card) | `TestReadFloat32_ExpectedMatrix`, `TestReadFloat32_RangeInvariant`, `TestReadFloat32_UnsupportedHDU` | astropy raw values + `fractions.Fraction` | `testdata/imaging/readfloat32_expected.json`, `testdata/imaging/i32_round.fits`, `testdata/imaging/i64_round.fits`, `testdata/imaging/u64.fits`, `testdata/imaging/u64_intcard.fits`, `testdata/imaging/f32_bzero_2d.fits`, `testdata/imaging/f64_bzero_2d.fits` | bit-exact | 1 |
| A BSCALE / BZERO card that is present but not a finite FITS number (string, quoted number `'10.0'` / `'0.5'`, `2.0.0`, hex float `0x1p3`, `1E400`) is an error from ReadFloat32, ReadFrame, AsImage and DecodeConfig, plain and tile-compressed, on a fresh Open and after `Header()` alike (cfitsio uses the default, or strtod's value, silently); D-exponent values (`5.0D-1`, `1.0D1`) are read on both paths | `TestReadFloat32_InvalidScalingCards`, `TestReadFloat32_DExponentScaling`, `TestReadFloat32_CompressedInvalidScaling`, `TestReadFloat32_ExpectedMatrix`, `TestAsImage_Matrix` | FITS 4.0 §4.2.1/4.2.3/4.2.4 syntax; samples classified from their raw cards (`gen_expected.py` `scaling_cards`), D-exponent values from astropy 8.0.2 (cfitsio 4.7.0 agrees); fitsverify rejects the invalid ones; cfitsio 4.7.0 behaviour recorded by the Siril harness | `testdata/imaging/i16_bscale_str_2d.fits`, `testdata/imaging/i16_bscale_bad_2d.fits`, `testdata/imaging/f32_bzero_huge_2d.fits`, `testdata/imaging/f32_bzero_qnum_2d.fits`, `testdata/imaging/i16_bscale_qnum_2d.fits`, `testdata/imaging/f32_bzero_hex_2d.fits`, `testdata/imaging/i16_dexp_2d.fits`, `testdata/imaging/f32_bzero_dexp_2d.fits` | exact (error); D exponents bit-exact | 1 |
| The scaling-card rule (`scalingCards`, the one rule of every image path) on ReadFloat32, ReadFrame, AsImage and DecodeConfig, identical on a fresh Open and after `Header()` (an in-memory `Set` card is applied as it would be written): the LAST card of a keyword is used on a plain image (BZERO 10 then 20 reads 20; BSCALE 2 then -1; BLANK 1 then 100), the FIRST on a tile-compressed one (as cfitsio and astropy each read it there: `comp_card_i16_bscale_dup_2d`, both readings in `refs`); on a tile-compressed header `BSCALE = 'abc'` (cfitsio 1, astropy fails) and `HIERARCH BSCALE` (both apply it there, astropy ignores it on a plain image) are refused like on a plain image, by every compressed read and `BSCALE()`/`BZERO()`, with `Validate` nil and Shape/BITPIX still readable (`comp_card_i16_bscale_str`, `comp_card_i16_bscale_hierarch`), and so is a BLANK outside ZBITPIX (`i16small_rice_blankcast`, `i16small_hcomp_blank65536`, `i16small_plio_blank65536`, `i16_zblankcol_blank2p32`, `comp_bad/i16_gzip1_w4_nullout`: cfitsio casts it to the tile type, astropy refuses or compares the integer); that card must be a valid value card (upper-case keyword in columns 1-8, `= ` in columns 9-10, a FITS number cfitsio converts alike, lower-case `e` accepted, BLANK an integer in the BITPIX range), else an error: `BZERO   =10.0`, `BZERO=10.0` and `=` in column 10/11 (cfitsio 10, astropy 0), `bzero` and `1.0d1` (cfitsio 0, astropy 10), `1 0.0` / `10.0 x` (cfitsio's strtod prefix, astropy refuses the file), an empty value, BLANK `100.0` / `'abc'` (both references drop the null marker silently), `BLANK   =100` | `TestReadFloat32_ScalingCardVariants`, `TestCfnullGoldensCoverImagingSamples`, `TestCompressedCorpusRead`, `TestCompressedNullOutsideZBITPIX`, `TestReadFloat32_InMemoryScalingCard`, `TestReadFloat32_InMemoryBlankCard`, `TestReadFloat32_InMemoryBigIntScaling`, `TestReadFloat32_ExpectedMatrix`, `TestReadFloat32_SirilConvention`, `TestAsImage_Matrix` | libcfitsio 4.7.0 `fits_read_pixnull` TDOUBLE (`cfnull/cards`, checked by the test and by `gen_expected.py`) and astropy 8.0.2 scaled `getdata` (asserted equal to the documented values for every accepted sample in `gen_expected.py`; both readings recorded under `refs`); fitsverify 4.7.0 | `testdata/imaging/card_*.fits`, `testdata/imaging/comp_card_*.fits.fz`, `testdata/imaging/cfnull/cards` | bit-exact; exact (error) | 1 |
| BLANK → NaN; compressed null value: per-tile ZBLANK column, else ZBLANK over BLANK; compressed BSCALE/BZERO as cfitsio reads them | `TestReadFloat32_ExpectedMatrix`, `TestReadFloat32_CompressedFloat_Cfitsio` | cfitsio 4.7.0 null flags (`cfnull`) | `testdata/imaging/comp_*.fits.fz`, `testdata/imaging/cfnull/golden` | bit-exact | 1 |
| Siril convention: bit-identical for BITPIX 8/16 (BZERO 0/32768), every other difference asserted | `TestReadFloat32_SirilConvention` | Siril functions verbatim, linked against libcfitsio (`sirilref`) | `testdata/imaging/sirilref/golden` | bit-exact | 1 |
| ReadFrame / ReadFrameRGB / WriteMono / WriteRGB on a real camera frame; orientation; planes, plain and tile-compressed; HDU selection (primary if NAXIS >= 2, else HDU 1, incl. a 1-D primary); NAXIS3 = 0 cube and zero WriteMono / WriteRGB dimensions refused; WriteMono/WriteRGB keep float32 bits verbatim (signaling NaNs, payloads, -0, ±Inf, denormals) | `TestReadFrame_RealCamera`, `TestReadFrame_RowOrientation`, `TestReadFrameRGB_Planes`, `TestReadFrame_PrimaryNAXIS1UsesHDU1`, `TestFrame_ReadFrame_ErrorsOnZeroPlaneCube`, `TestFrame_WriteMono_ZeroDims`, `TestFrame_WriteRGB_ZeroDims`, `TestWriteMono_AstropyReadsVerbatim`, `TestWriteRGB_AstropyPlanes`, `TestWriteMono_TruncatesExisting` | Siril harness; astropy 8.0.2; cfitsio null flags (`cfnull`) for the `.fz`; fitsverify; documented contract (refusals) | `testdata/imaging/asi533_u16_crop.fits`, `testdata/imaging/u16_orient.fits`, `testdata/imaging/u16_rgb.fits`, `testdata/imaging/comp_u16_rgb.fits.fz`, `testdata/imaging/frame_naxis1_primary.fits` | bit-exact | 1 + 2 live |
| AsImage / image.Decode lossless or `ErrNotImageCompatible` (incl. float BZERO without BSCALE, exact and inexact); FloatColor.RGBA clamp | `TestAsImage_Matrix`, `TestAsImage_RowOrientation`, `TestImageDecode_Registered`, `TestFloatColor_RGBA` | astropy-read values; documented mapping | `testdata/imaging/stdimage_expected.json`, `testdata/imaging/f32_bzero_inexact_2d.fits` | exact | 1 |

### Statistics (`fits/stats`)

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| `MinMax` and `Mean`: NaN skipped, (0, 0) / 0 for empty or all-NaN input, ±Inf by IEEE rules; `Mean` within its documented recursive-summation bound n·ε·mean(\|x\|) of the exact mean, on every golden input (all element types); comparisons NaN-aware (a NaN result fails) | `TestMinMax`, `TestMean`, `TestAllNaN`, `TestEmptySlice` | numpy 2.4.6 min/max/mean; exact rationals | `stats/testdata/golden/stats_golden.json` | min/max exact; Mean γ_{n+1}·Σ\|x\|/n (Higham) | 1 |
| `MeanStdev`: float64 one-pass sample stdev (N−1), agrees with numpy `std(ddof=1)` within the one-pass bound; loses precision on large offsets (1e9 + {0,1} → 95.834358265834 on amd64, 96.153395548339716 on arm64, vs true 0.50025) | `TestMeanStdev_Numpy`, `TestMeanStdev_OnePassFormula`, `TestMeanStdev_OffsetCancellation` | numpy 2.4.6 `std(ddof=1)`; Python reproduction of the formula per FMA-fusion model (none / fused accumulation / fused accumulation and subtraction, exact-rational `fma`); exact rationals | `stats/testdata/golden/stats_golden.json`, `stats/testdata/gen_golden.py` | gamma_n rounding bound (NaN-aware); bit-exact vs the build's fusion model (amd64, GOAMD64=v3, arm64 measured; other architectures skip the bit pin loudly, but the non-finite cases' mean bits and stdev class (NaN / ±Inf, the same in all three fusion models) are checked on every build) | 1 |
| `MeanStdevTwoPass`: accurate generic sample mean/stdev (corrected two-pass with Neumaier sums, incl. a value larger than the running sum after many small ones); int8/uint32/uint64/int64 and named-float32 inputs; int64/uint64 beyond 2^53 rounded to float64 first (data clustered at 2^60 and 2^63); the mean is correctly rounded (20000 regenerated sets of integers ±(2^40..2^62 + k), as int64 and float64; a remainder from a rounded product, a missing remainder correction or a dropped low word each fail 2948-3328 of them) | `TestMeanStdevTwoPass_Exact`, `TestMeanStdevTwoPass_RoundsLargeIntegers`, `TestMeanStdevTwoPass_MeanCorrectlyRounded` | Python `fractions` / exact integer division, sqrt via decimal at 80 digits | `stats/testdata/golden/stats_golden.json`, `stats/testdata/golden/mean_rounding.json` (sha256 of the 20000 means + 35 samples) | mean bit-exact (correctly rounded); stdev ≤ 2 ULP (measured 1) | 1 |
| `MeanStdevSiril`, `Percentile`/`Median`/`MAD` (float32), `FilterNonZero` (every committed or generated input, incl. 2e7 elements, as a list or count + sha256), QuickSelect vs `quickmedian_float`; the 424 golden inputs (incl. 2^24+3 and 2e7 elements) required by name and size; the MeanStdevSiril doc's 43/50-of-424 rounding counts | `TestSirilOracle_Float32`, `TestMeanStdevSiril_Edges`, `TestMeanStdevSiril_DocCounts`, `TestQuickSelect_AgainstSort` | Siril 1.4.0 `siril_stats_float_sd`, `findMinMaxPercentile`, `siril_stats_float_mad`, `reassign_to_non_null_data_float` compiled verbatim | `stats/testdata/sref/golden` | bit-exact (n up to 2e7) | 1 |
| Generic (non-float32) Percentile/Median/MAD path for every element type (int8, int16, int32, int64 to 2^61, uint8, uint16, uint32 above 2^31, uint64 above 2^63, float64, a named float32 type), incl. float64 and uint16 inputs of 2e7 elements (CDF walk above 2^24) | `TestSirilOracle_Generic` | the same Siril source with `float` → `double` | `stats/testdata/sref/golden` | bit-exact | 1 |
| `Percentile` is within one histogram bin of the order statistic x_(ceil(p·n)); p = 0 / 1 give min / max exactly (finite ranges) | `TestPercentile_WithinOneBin` | numpy 2.4.6 `np.sort` order statistics | `stats/testdata/golden/stats_golden.json` | one bin (max−min)/(min(n, 65536)−1) plus 4 ULP of the range at the working precision (NaN-aware) | 1 |
| `SigmaClip` follows `astropy.stats.SigmaClip`: symmetric and asymmetric (sigma_lower ≠ sigma_upper) bounds, maxiters 1 / 10 / None (`maxIter <= 0`, no cap: a fixture needs 101-183 iterations), iteration count, mean and median centre, all element types; an iteration that rejects every value (k = 0.5 on {0 x5, 1 x5, 100}, float32 and float64, maxiters None / 10 / 2: the zero result with astropy's NIter 3 / 3 / 2) | `TestSigmaClip_Astropy`, `TestSigmaClip_Degenerate` | astropy 8.0.2 `SigmaClip` with the exactly rounded mean and population std (exact integer arithmetic) or Siril's histogram median as cenfunc/stdfunc | `stats/testdata/golden/stats_golden.json` | NGood and NIter exact, median bit-exact, survivors' mean bit-exact (correctly rounded), stdev ≤ 2 ULP of exact (measured 1) | 1 |
| `BuildHistogram` matches `numpy.histogram`, including numpy's one-off index corrections downwards and upwards (values on its edges and one float64 either side, 26 ranges: 106 and 238 corrections); where numpy raises (single values at 2^53 and at 2^52 with 2 bins, max−min overflowing, subnormal and 2^60 ranges too narrow for nbins edges, integers clustered at 2^60 / 2^63) it returns `ErrHistogramRange` and a zero Histogram | `TestBuildHistogram_Numpy`, `TestBuildHistogram_NumpyEdgeCorrections`, `TestBuildHistogram_NumpyRanges`, `TestBuildHistogram_Inf` | numpy 2.4.6 (bins, or its ValueError / IndexError; its index computation replayed to count the corrections) | `stats/testdata/golden/stats_golden.json`, `stats/testdata/golden/histogram_ranges.json`, `stats/testdata/golden/histogram_edges.json` | counts/edges/CDF bit-exact; error exactly where numpy raises | 1 |
| `Histogram.Percentile`: the left-continuous inverse of the piecewise-linear CDF (p = 0: lower edge of the first non-empty bin; flat CDF: left end) | `TestBuildHistogram_Numpy`, `TestHistogramPercentile_Scipy` | exact rationals for the documented definition (all p, incl. clamped and flat points); scipy 1.17.1 `stats.rv_histogram(...).ppf` (215 points, 0 < p < 1 off flat parts) | `stats/testdata/golden/stats_golden.json` | definition: 4 ULP of the edges; scipy: (nbins+4)·ε·Total/count bin widths (scipy's float CDF) + 4 ULP; both NaN-aware | 1 |
| Non-finite contract (Percentile/Median/MAD → NaN for ±Inf data, a NaN p, or max−min overflowing the working precision, e.g. float32 {−3e38, 3e38, 1}, which used to give a median of 3e38) | `TestPercentile_NonFinite`, `TestPercentile_RangeOverflow` | documented contract only (Siril undefined there) | — | exact | 1 |
| Goldens reproduce (all five files; each golden directory holds exactly the files its generator writes); Siril excerpt verbatim; large inputs, the 400 small Siril cases and the Siril `gen` inputs are xorshift32 formulas / portable generators run identically by the C harness, `gen_golden.py` and the Go tests (not committed; every generated input over 64 values carries the sha256 of its bits, checked first); large outputs (histogram edges/CDF, FilterNonZero lists of generated inputs) are sha256 plus samples; no `std::*_distribution` (clang++ 22.1.8 reproduces the Siril file byte for byte) | `TestGoldensReproduce` | Siril 1.4.0 tarball sha256 | `stats/testdata/sref/extract_siril.py`, `stats/testdata/sref/sref.cc` | byte-exact | 2 |

### Stretches (`fits/stretch`)

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| Oracle harness reproduces the 18 siril-cli goldens and every oracle golden; goldens are exact digests; a real siril-cli 1.4.0 reproduces the 18 presets and all 146 oracle goldens that are siril commands, the 7 with ±Inf input included (given to siril-cli as markers that pixel math turns into ±Inf, since its FITS loader rescales data whose maximum exceeds 10) | `TestSirilGoldensReproduce`, `TestSirilCLIDigests`, `TestSirilCLIDigestsFromHistory`, `TestSrefManifest`, `TestSirilCLIReproduce` | Siril 1.4.0 source compiled verbatim (`-ffp-contract=off`); the original siril-cli pixel files (git `c17580b`); siril-cli 1.4.0 built from the release tarball (`stretch/testdata/build_siril_cli.sh`) | `stretch/testdata/golden/sref`, `stretch/testdata/golden/stretch_*.json`, `stretch/testdata/golden/siril_cli_sha256.json` | byte-exact; each golden is a sha256 over every pixel (plus a sample), so one differing bit fails | 1; 2: `TestSirilGoldensReproduce`, `TestSirilCLIDigestsFromHistory`, `TestSirilCLIReproduce` |
| Asinh (mono, RGB, `-human`, offset incl. negative and 0.99, beta=0, beta in (0, 1), beta = 3e38, input > 1, tiny and subnormal input incl. Siril's 0/0 NaN at beta = 0.5; NaN/infinite beta or offset, beta < 0 and offset >= 1 are `ErrAsinhParams`, compared with siril-cli's decision on 15 commands) | `TestAsinhTiny_SirilOracle`, `TestAsinh_Siril`, `TestAsinh_SirilOracle`, `TestAsinhEdges_SirilOracle`, `TestAsinhEdgeParams_SirilOracle`, `TestAsinhMTF_ParamErrors`, `TestSirilCLIArgs` | siril-cli goldens; Siril 1.4.0 oracle | `stretch/testdata/golden/stretch_asinh_strong.json` | 0 ULP | 1 |
| Autostretch (mono, unlinked, linked, inverted branch, no-good-pixel fallback; negative pixels count in the statistics; a median of exactly 0.5 takes the normal branch; linked c0/n, c1/n in float32 on 96 random images; ±Inf excluded from the statistics, equal to Siril with those pixels set to 0; targetBG outside [0, 1] is `ErrAutostretchParams` exactly where siril-cli refuses it, NaN arguments refused, ±Inf shadows clip gives siril-cli's pixels) | `TestSirilCLIArgs`, `TestAutostretchParams_PanicOnBadArgs`, `TestAutostretch_Siril`, `TestAutostretchLinked_Siril`, `TestAutostretch_SirilOracle`, `TestAutostretchParams_SirilOracle`, `TestAutostretchInvertedBoundary_SirilOracle`, `TestNegative_SirilOracle`, `TestAutostretchMedianHalf_SirilOracle`, `TestAutostretchSweep_SirilOracle`, `TestInf_SirilOracle` | siril-cli goldens; Siril oracle | `stretch/testdata/golden/stretch_auto_linked.json`, `stretch/testdata/golden/sweep/autostretch.json` | 0 ULP | 1 |
| MTF, Linear (incl. negative input, -0 -> +0 at BP = 0, BP = 1; BP outside [0, 1] or NaN is `ErrLinearParams`, exactly siril-cli's refusals; MTF with midtones 0 and 1, shadows >= 1 or < 0, highlights <= 0 or > 1, shadows >= highlights; NaN/infinite MTF parameters or midtones outside [0, 1] are `ErrMTFParams`) | `TestLinearBP1_SirilOracle`, `TestSirilCLIArgs`, `TestMTF_SirilOracle`, `TestMTFEdgeParams_SirilOracle`, `TestAsinhMTF_ParamErrors`, `TestMTFValue`, `TestLinear_Siril`, `TestLinear_SirilOracle`, `TestLinear_Basics`, `TestLinear_NegativeZero` | Siril oracle / siril-cli | `stretch/testdata/golden/sref` | 0 ULP, sign of zero included | 1 |
| GHT (output not clamped, small D, B threshold incl. B = ±0.001 and B = ±(double)1e-3f = ±0.0010000000474974513 exactly, B < 0 rounding; HP = 0 with LP = SP = 0 via `GHTParams.HPZero`; accept/refuse equal to siril-cli on 28 edge commands except the listed GUI extensions and NaN; GHTsetup coefficients bit for bit on 1024 random parameter sets in four B families; pixels exactly at LP/SP/HP and either side; 0 and -0 stay +0; D = 0 (the GUI's clamp alone) in every B family; out-of-range parameters are an error; the float32 conversions that stop FMA fusion, under GOAMD64=v3 and arm64) | `TestGHTD0_SirilOracle`, `TestFMAGuards_GOAMD64v3`, `TestCrossArch_arm64`, `TestGHT_Siril`, `TestGHT_SirilOracle`, `TestGHTEdges_SirilOracle`, `TestGHTOutputAbove1_SirilOracle`, `TestGHT_Basics`, `TestGHTSweep_SirilOracle`, `TestGHTBThreshold_SirilOracle`, `TestGHTBThresholdF32_SirilOracle`, `TestGHTHPZero_SirilOracle`, `TestGHT_RejectsOutOfRange`, `TestSirilCLIArgs` | Siril oracle (`sref ghtcoef` prints the verbatim GHTsetup struct); siril-cli's range checks (`core/command.c:2363-2387`); glibc libm probe | `stretch/testdata/golden/stretch_ght_strong.json`, `stretch/testdata/golden/sweep/ght_*.json.gz` | bit-exact except at the pixels each golden lists as libm-sensitive (where the verbatim oracle and the same code with double-rounded logf/powf/expf disagree: 0 to 58 per golden); there a value is accepted only if a 1-ULP change of one glibc logf/powf/expf result reproduces Siril (glibc is not correctly rounded) | 1; 2: `TestFMAGuards_GOAMD64v3`, `TestCrossArch_arm64` |
| NaN, ±Inf and negative input of every operation | `TestNaN_SirilOracle`, `TestNaN_AutostretchParams_SirilOracle`, `TestLabNaN_cv2`, `TestCvRound_cv2`, `TestInf_SirilOracle`, `TestNegative_SirilOracle` | Siril oracle on raw NaN / ±Inf / negative input; cv2 5.0.0 | `stretch/testdata/golden/nan_cv2.json` | 0 ULP | 1 |
| CLAHE equals Siril's OpenCV pipeline (incl. images smaller than the tile grid and 1 pixel wide; the clip limit converted as x86-64's `static_cast<int>`, so clip*tilePixels/65536 >= 2^31 or +Inf is limit 1, NaN no limit, on amd64 and arm64; clip and tile accept/refuse compared with siril-cli) | `TestCLAHEClipConversion_SirilOracle`, `TestSirilCLIArgs`, `TestCLAHETiny_SirilOracle`, `TestCLAHE_SirilOracle`, `TestCLAHERGB_cv2Rebuild`, `TestCLAHECore8bit_cv2`, `TestCLAHECore16bit_cv2`, `TestCLAHECoreClip_cv2`, `TestBGR2Lab_cv2`, `TestLab2BGR_cv2`, `TestCLAHE_BadArgsPanic` | Siril `cvClahe_float` + OpenCV 5.0.0; opencv-python 5.0.0.93 | `stretch/testdata/golden/lab_cv2.json`, `stretch/testdata/golden/clahe_clip_cv2.json` | bit-exact | 1 |
| Presets reproduce their Siril commands; JPEG output and its quality argument (quantization tables equal libjpeg-turbo's at each quality, out-of-range clamped); quantisation | `TestDefaultPresets_Siril`, `TestDefaultPresets_SirilOracle`, `TestFindPreset`, `TestToRGBImage`, `TestEncodeJPEG`, `TestEncodeJPEG_QualityTables`, `TestNormalizeChannels` (incl. both sides of the [0, 1] skip boundary; ±Inf left out of the range and mapped to 1 / 0; a finite range wider than MaxFloat32 scaled in float64) | siril-cli goldens; Siril oracle; libjpeg-turbo 3.1.2 via cv2 (decode and `imencode` tables); exact rationals | `stretch/testdata/golden/inputs/synth.f32`, `stretch/testdata/golden/jpeg_dqt_cv2.json` | 0 ULP / exact | 1 + 2 live |

### Test infrastructure

| Capability | Proving tests | Truth source | Samples / goldens | Tolerance (why) | Tier |
|---|---|---|---|---|---|
| Tier 2 gating: a missing tool, harness or variable (`FITS_ASTROPY_PYTHON`, `FITS_SIRIL_SRC`, `FITS_HEALPIX_CXX`, `FITS_SIRIL_CLI`) is a failure naming it; cfitsio goldens reproduce and cover every sample; fitsverify's verdicts on two cfitsio samples | `TestMissingToolFails`, `TestAstropyPythonBrokenFails`, `TestAstropyPythonLive`, `TestRequireToolLive`, `TestGoldensCoverEveryTestdataSample`, `TestCfrefDumpReproducesGoldens`, `TestCfrefDumpRejectsNonFITS`, `TestFitsverifySamples` | the harnesses themselves | `testdata/cfref/golden`, `internal/testutil/interop.go` | byte-exact | 1 + 2 live; 2: `TestMissingToolFails`, `TestAstropyPythonBrokenFails`, `TestAstropyPythonLive`, `TestRequireToolLive`, `TestCfrefDumpReproducesGoldens` |
| Tier 1 starts no program: no Go file a plain `go test` compiles (linux, darwin, windows; amd64, arm64) uses `os/exec` other than `exec.Command(os.Args[0], ...)`, `os.StartProcess` or `syscall.Exec`; every `*_interop_test.go` carries `//go:build interop`; the scanner catches each form | `TestTier1StartsNoPrograms`, `TestTier1ScannerCatches` | the Go sources (go/parser, go/build constraints) | `internal/testutil/hermetic_test.go` | exact | 1 |
| Write goldens: a reference tool's reading of a file is keyed by the tool, its committed harness source, the arguments and the sha256 of every input; tier 1 replays the recorded reading of exactly the bytes the library writes (so it also requires those bytes), tier 2 re-runs the tool and compares; large pixel readings as an exact digest resolved against a candidate | `TestWriteImageCfitsioReadsBack`, `TestCompressedWriterInterop`, `TestFlushEditMatrixMatchesCfitsio`, `TestEditFileCopyMatchesCfitsioCopyHDU`, `TestWriteTableInteropCfitsio`, `TestWriteChecksumCfitsioVerifies` | cfitsio, fitsverify, funpack, astropy, wcslib readings recorded by `FITS_RECORD=1` tier 2 runs | `internal/testutil/record.go`, `internal/testutil/digestfloats.go`, `testdata/recorded` | exact (a changed byte is a missing recording, a failure) | 1 + 2 live |
| Large cfitsio golden arrays (more than 1024 values: image pixels, flattened table columns) are stored as an exact sha256 digest over a canonical big-endian encoding plus explicit samples (corners, first min/max, every null as runs, 64 strided); a result compares equal only if every element, null flag, sign of zero and NaN-ness matches, and a mismatch names the differing sampled indices; every committed digest is self-consistent and no array over 1024 values is left explicit | `TestCfrefDigestCatchesOneChange`, `TestCfrefDigestsConsistent`, `TestDigestCanonicalEncoding`, `TestCfrefDumpReproducesGoldens` | libcfitsio 4.7.0 via `testdata/cfref/cfref.c` (self-contained SHA-256, agreeing with Go `crypto/sha256` on every resolved digest) | `testdata/cfref/golden`, `testdata/cfref/PROVENANCE.md`, `internal/testutil/cfref_values.go` | exact (a hash mismatch is a failure; one ULP on one of 90000 pixels is caught) | 1; 2: `TestCfrefDumpReproducesGoldens` |
| Golden sets and samples reproduce from their generators (rule 2); the cfref pool (mid-size arrays stored once, kind/count/sha256 checked on load) | `TestCrefGoldensReproduce`, `TestCfrefDumpReproducesGoldens`, `TestCfrefPoolComplete`, `TestCfrefPoolRefusesBadFiles`, `TestCompressedCorpusReproduce`, `TestImagesCorpusReproduce`, `TestTableSamplesReproduce`, `TestHeaderSamplesReproduce`, `TestImagingSamplesReproduce`, `TestHealpixMapSamplesReproduce`, `TestTNXSyntheticHeadersReproduce`, `TestHighOrderCaseGeneratorReproduce`, `TestCfnullGoldensReproduce`, `TestSirilrefGoldensReproduce`, `TestImagingExpectedReproduce`, `TestHealpixGoldensReproduce`, `TestSirilGoldensReproduce`, `TestCLAHECv2GoldensReproduce` | libcfitsio 4.7.0, Siril 1.4.0 source, astropy 8.0.2, healpy 1.20.1 + healpix_cxx 8.3.0, opencv-python 5.0.0.93 | `internal/testutil/reproduce.go` | byte-exact | 1; 2: `TestCrefGoldensReproduce`, `TestCfrefDumpReproducesGoldens`, `TestCompressedCorpusReproduce`, `TestImagesCorpusReproduce`, `TestTableSamplesReproduce`, `TestHeaderSamplesReproduce`, `TestImagingSamplesReproduce`, `TestHealpixMapSamplesReproduce`, `TestTNXSyntheticHeadersReproduce`, `TestHighOrderCaseGeneratorReproduce`, `TestCfnullGoldensReproduce`, `TestSirilrefGoldensReproduce`, `TestImagingExpectedReproduce`, `TestHealpixGoldensReproduce`, `TestSirilGoldensReproduce`, `TestCLAHECv2GoldensReproduce` |
| HEALPix and imaging goldens: float arrays as base64 big-endian float64 (exact), integer and exactly-checked float result arrays over 256 values as a sha256 digest of their big-endian int64 / IEEE bits plus explicit samples; the decoders, the digest and its sample checks are pinned on hand-made vectors | `TestGoldenCodec`, `TestImagingGoldenCodec`, `TestHealpixGoldensReproduce`, `TestImagingExpectedReproduce` | published sha256 vector; Python `struct`/`hashlib`/`base64` (the generators) | `healpix/testdata/golden`, `testdata/imaging/readfloat32_expected.json` | exact (a digest mismatch is a failure) | 1; 2: `TestHealpixGoldensReproduce`, `TestImagingExpectedReproduce` |
| Every tile-compressed imaging sample has its cfitsio null/physical golden and every golden is checked against the reader | `TestCfnullGoldensCoverImagingSamples` | cfitsio `fits_read_pixnull` TDOUBLE | `testdata/imaging/cfnull/golden` | exact | 1 |
| This ledger names only existing tests and samples | `TestValidationLedger` | the repository | `docs/VALIDATION.md` | exact | 1 |
| Size budget of the committed reference data: no tracked file over 1 MB (except `reference/fits_standard_v4.0.pdf`), at most 20 MB tracked under testdata directories, at most 25 MB tracked in total; large text goldens are committed gzipped (`<name>.gz`, gzip level 9 without name or time) and read through `testutil.ReadGolden`, and every reproduce test and `make verify` compares the decompressed bytes with the harness output | `TestRepositorySizeBudget`, `TestHdrGoldensReproduce`, `TestCfrefDumpReproducesGoldens`, `TestHealpixGoldensReproduce`, `TestImagesCorpusReproduce`, `TestSirilGoldensReproduce`, `TestGoldensReproduce` | the git index, read in Go (`testutil.TrackedFiles`; tier 2 checks it against `git ls-files`, `TestTrackedFilesMatchesGit`) | `internal/testutil/golden.go`, `internal/testutil/gitindex.go`, `testdata/images/cmp_golden.py` | exact (bytes) | 1; 2: `TestHdrGoldensReproduce`, `TestCfrefDumpReproducesGoldens`, `TestHealpixGoldensReproduce`, `TestImagesCorpusReproduce`, `TestSirilGoldensReproduce`, `TestGoldensReproduce`, `TestTrackedFilesMatchesGit` |
| Bit-exact claims hold on builds whose compiler fuses multiply-adds: no fusable float product in library code (only the frozen `stats.MeanStdev`, pinned per build); tier 1 run on GOAMD64=v3 and arm64 (`TestFMAGuards_GOAMD64v3`, `TestCrossArch_arm64`, tier 2) | `TestNoFusableFloatProducts`, `TestMeanStdev_OnePassFormula` | Go spec (explicit conversions round); the same C/Python goldens on every build | `fma_guard_test.go` | exact | 1; 2: `TestFMAGuards_GOAMD64v3`, `TestCrossArch_arm64` |

## Deliberate differences from the reference

Each is documented in the function doc and asserted by a test.

- **Images**: a value that cannot be stored is `ErrNumOverflow` with
  nothing written, where cfitsio stores the range limit and returns
  NUM_OVERFLOW (until round 2 the library clamped silently, the PR #2
  behaviour, 0258acd); also refused where cfitsio stores a wrong value
  without a status (BITPIX 64 scaled value exactly 2^63, which cfitsio
  stores as INT64_MIN via an undefined cast; -2^63 with BZERO 1000; a
  finite value scaled to ±Inf on float BITPIX); a header
  with a structural keyword repeated with a different value is refused
  (cfitsio opens it; the same value repeated opens);
  `WriteImage` refuses scaling cards cfitsio would read differently from
  the card (cfitsio's writer refuses only BSCALE 0); float32 reads beyond range give ±Inf (cfitsio ±FLT_MAX);
  `ReadPixelsMasked` treats ±Inf as valid; truncated XTENSION header is an
  error; integer targets refused for BITPIX 64 with a nonzero integral BZERO
  other than 2^63; `WriteImage`/`OverwritePixels` with BITPIX 64 and an
  integral BZERO other than 0 and 2^63 store v - BZERO exactly where
  cfitsio's `fits_write_img` rounds (v - BZERO) through a double with
  status 0 (v = 2^62 + 1, BZERO 1000: 4611686018427386905 here,
  4611686018427386880 in cfitsio; values within 2^53 of BZERO agree; the
  committed write matrix holds only those). ReadFloat32 / AsImage refuse a BSCALE or BZERO card that
  is not a finite number (cfitsio `fits_read_img` uses the default), and a
  last BSCALE / BZERO / BLANK card that is not a valid value card (cfitsio
  and astropy read those differently from each other, or drop them).
- **Framing**: a data unit cut short is `ErrTruncatedData` at Open (cfitsio
  opens and fails the read); an extension missing BITPIX/NAXISn is
  `ErrMissingRequiredKeyword` (cfitsio silently reports one HDU fewer;
  the same for a table extension without PCOUNT/GCOUNT/TFIELDS); a
  non-string XTENSION is refused (cfitsio reads it as an IMAGE).
  `ColumnByName` and `HDUByName` are case-sensitive (cfitsio's
  `fits_get_colnum` and `fits_movnam_hdu` are not by default). An extension
  of an unknown type is refused with `ErrUnsupportedExtension` (cfitsio
  reads it as an IMAGE_HDU).
- **Headers**: TRUE/Tx, hex tokens and complex with an empty part stay
  `TypeUnparsed` (cfitsio reads them); HIERARCH blank runs collapsed on read,
  refused on write; joined CONTINUE comment kept whole past 72 characters;
  `Int` does not truncate float cards (`ErrWrongType`). astropy joins CONTINUE
  comments differently and is not used as truth there.
- **Tables**: Inf/-0/subnormals kept in E/D/C/M; ASCII integers exact and long
  fields correctly rounded where cfitsio goes through a double; ASCII TNULL
  padded per §7.2.2; invalid logical bytes are an error (cfitsio reads
  0x01 as true); stored 0 with
  TSCAL -1 is +0 (cfitsio; astropy gives -0). TSCAL/TZERO literals are
  exact decimals: TZERO = 9223372036854775809 or 9223372036854775807 is
  that integer (cfitsio and astropy compare the double and apply the
  uint64 convention), integer targets on the float64 scaling path get the
  exact value, and a float TZERO literal near but not at 2^63 refuses
  integer targets. ASCII "-0" under TSCAL/TZERO is exactly TZERO
  (cfitsio adds the double of TZERO). E/D values and ASCII fields under
  TSCAL 1 and an integer TZERO no double holds are stored + TZERO
  correctly rounded (cfitsio and astropy add TZERO's double); float32
  targets on those exact paths get the float32 nearest the exact value
  (cfitsio's TFLOAT rounds through its double when it scales). String-array layouts cfitsio and the standard read
  differently (a TDIM shorter than the cell, rAw with r not a multiple of
  w) are errors (cfitsio splits them into repeat/width strings). TSCAL/TZERO
  on C/M columns is an error (cfitsio adds TZERO to the imaginary part,
  astropy drops it; the standard adds it to the real part only). A
  descriptor past the heap after a THEAP gap is an error (cfitsio reads the
  padding; fitsverify does not check it).
- **Compression**: fallback tiles of scaled/BLANK images follow funpack and
  astropy, not `fits_read_img`; lossy HCOMPRESS_1 pixels of a 16-bit image
  outside the type range are clipped as funpack and cfitsio's TSHORT read
  do, not returned unclipped as its TDOUBLE read does nor refused as its
  TUSHORT read of a uint16 image does (status 412); an
  HCOMPRESS_1 reconstruction outside int32 is refused (cfitsio wraps it); writer drops caller Z* keywords and
  duplicates (cfitsio copies them); writer refuses HCOMPRESS edge tiles of
  1..3 px (cfitsio enlarges); a tile on which cfitsio's quantizer overflows
  (int64 zero point, int32 pixel, a pixel on NULL_VALUE) is stored losslessly
  (cfitsio writes wrong pixels with status 0); a real float pixel equal to FLOATNULLVALUE (float64:
  DOUBLENULLVALUE) stays a value (`fits_img_compress` nulls it; funpack returns
  NaN for the float32 one); a caller BLANK on
  float input and a non-logical caller SIMPLE are refused (cfitsio copies
  them); no default `EXTNAME = 'COMPRESSED_IMAGE'`.
- **WCS**: the refusals in `TestHdrStricterThanWcslib` and the `TestNewRefuses*`
  tests; TPV missing PVi_1 defaults to 1 (registry; wcslib uses 0); a
  `RA---TPV` header (no `-SIP`) with SIP keywords ignores them, as astropy
  8.0.2 does, where wcslib 8.10 applies both distortions: on the real
  ATLAS header of PR #4 (`hdrcases_differ/real_atlas_k2_tpv_sip`) the SIP
  set is a copy of the TPV correction (wcslib on each part alone: 2.31e-7
  arcsec apart, (-12.28, 6.06) px at the stamp corner) and wcslib's
  reading of the header as written is 24.45-26.84 arcsec (13.28-14.59 px)
  off at every point, while astropy's (and ours, to 3.7e-10 arcsec) is the
  TPV one (`TestHdrDeliberateDifferences`, `TestATLASTPVPlusRedundantSIP`,
  `TestATLASSIPCopiesTPV`; `RA---TPV-SIP` still applies both); SIP keywords on a TNX header are
  refused (wcslib applies SIP + TAN and ignores the WAT correction when the
  SIP cards come last); sky->pixel on a DSS header far outside the cutout
  returns the converged pixel where wcslib's disx2p can fail to converge
  (6000 px out on a degree-5 plate solution: `dss_far_outside_plate`, where
  the pixel is checked against wcslib's own pixel->sky); accepted although
  wcslib refuses them (`TestHdrLenientThanWcslib`): SIP without AP/BP, or
  with AP or BP lacking any coefficient (wcslib: "Degree must not exceed
  9, got -1"; its inverse only seeds an iteration on A/B, ours too), and a
  DSS header whose PLATEID card follows the plate keywords (wcslib: "No
  distortion type, yet DQia keyvalues are present");
  `FK5` precession is IAU 1976 (astropy FK5→FK5 uses IAU 2006, ~0.3 arcsec
  per century apart).
- **HEALPix**: Ang2Pix/Vec2Pix use Go's math.Cos/Sin/Atan2, so a
  direction within about an ULP of a pixel boundary (boundary-ambiguous:
  its pixel changes under a few-ULP nudge) can be given the pixel on the
  other side than healpy's; such directions have no meaningful correct
  side, and the tests count and cap them. Vec2Pix rescales every
  vector, so it differs from healpy where healpy's squared components
  overflow or lose bits (e.g. near-pole vectors scaled by 2^-490).
- **Stats/stretch**: `SigmaClip(CenterMedian)` uses the std about the mean
  (astropy semantics) and drops ±Inf. SigmaClip's mean is correctly rounded
  and its std within 2 ULP, so it is checked against astropy given the exact
  mean and std; astropy with numpy's own `np.mean`/`np.std` can keep a
  different set on data with a large offset (numpy's std of integers at
  2^63 + [0, 10000) is 3119.4, the exact value 2981.4).

## Not supported / known gaps

### Not supported (returns an error, or not implemented)

- Headers with a duplicate structural keyword (`ErrDuplicateKeyword`),
  `WriteImage` with scaling cards cfitsio would read differently and
  `OverwritePixels` on a BSCALE = 0 HDU (`ErrInvalidScaling`); values
  outside the stored range of `WriteImage` / `OverwritePixels`
  (`ErrNumOverflow`; there is no clamping write).
- Random-groups HDUs (`ErrRandomGroups`), tiled table compression
  (`ErrTiledTableCompression`), time
  coordinates (§9), iterator/row-filter/grouping/region/extended-filename
  features of cfitsio, IRAF `.imh`.
- Tables: string-array layouts the references disagree on (a TDIM that
  covers fewer characters than the repeat count, rAw with r not a multiple
  of w), TDIM string arrays in VLA columns, and writing string arrays (all
  errors); TSCAL/TZERO on complex columns (an error: cfitsio and astropy
  each disagree with the standard, `TestComplexScalingReferencesDisagree`);
  integer targets under a float TZERO literal near but not at 2^63 (an
  error); a TScal/TZero whose exact literal does not fit a card (an
  error); general TSCAL inverse on write; TSCAL/TZERO/TDISP in the ASCII
  writer (hdr keys refused). 0P/0Q columns are proven only against the
  standard (no reference tool reads them).
- Compression read: RICE_1 tiles of an unquantized float image other than
  4-byte tiles (cfitsio reads 1-/2-byte ones past the tile end);
  HCOMPRESS_1/PLIO_1 tiles of an unquantized float image (no writer, no
  reference); quantized GZIP tiles of 8 bytes (cfitsio returns unwritten
  memory, funpack and astropy zeros); tile values outside the ZBITPIX
  range, including PLIO_1 ones and HCOMPRESS_1 ones of 8-bit images
  (and any HCOMPRESS_1 reconstruction outside int32) (refused for every target type; cfitsio fails only a read as the
  image type, a TDOUBLE read returns them; lossy HCOMPRESS_1 pixels of a
  16-bit image are clipped instead, see Deliberate differences), and a
  null tile value outside it in an unmasked read
  (the masked reads return the image); a BLANK/ZBLANK keyword or ZBLANK
  column value outside int32 (cfitsio reads them as a C int: it ignores
  such a keyword, falling back from ZBLANK to BLANK, and fails on such a
  column value with status 412; astropy 8.0.2 refuses the file); PLIO_1 with BZERO 32768 at ZBITPIX 8
  or 32 (cfitsio's writer and reader disagree there); a ZCMPTYPE cfitsio
  does not recognise (case-sensitive); ZNAXIS < 1 or more than 6 axes
  (cfitsio BAD_NAXIS); a ZNAXISn below 1 (cfitsio: 414); a ZBITPIX that is
  no BITPIX (cfitsio decodes the tiles anyway); an HCOMPRESS_1
  reconstruction outside the int32 range (cfitsio and funpack wrap it).
  NOCOMPRESS has no cfitsio truth (astropy only: cfitsio 4.7.0 cannot read
  it).
- Extensions other than IMAGE, TABLE, BINTABLE and their legacy aliases
  (e.g. FOREIGN, DUMP): framed, then refused with `ErrUnsupportedExtension`.
- Compression write: more than 6 axes, int64/uint64 (cfitsio cannot write
  them), float HCOMPRESS/PLIO (float32 and float64), a QLevel of NaN or
  -Inf,
  noise-relative HCOMPRESS scale, PLIO uint16, NOCOMPRESS, integer
  BLANK/ZBLANK nulls (a caller ZBLANK is dropped), a caller BLANK on float
  input, a caller SIMPLE that is not logical and a null value of the other
  width (NullValue to the float64 writer, NullValue64 to the float32 one)
  (refused).
- WCS: lookup-table distortions (CPDIS/CQDIS/D2IM), SIP on alternate
  descriptions, the other refusals listed in README; RADESYS/EQUINOX/MJD-OBS
  are parsed but not used by `transform`; FK4 at other equinoxes; proper
  motion.
- Editing (`OpenForEdit`, `EditFile`) off Linux and macOS (Windows, the
  BSDs, ...): refused with `ErrEditUnsupported` (the edit lock, the
  reader/writer lock and the xattr carry-over exist for Linux and macOS
  only). That code is only compiled (`GOOS=windows go vet`).
- Editing on macOS: supported, with the Linux guarantees (lock_darwin.go:
  OFD record locks for both locks, F_FULLFSYNC for every sync, xattrs and
  resource fork carried over), proven by the same suite run on a real Mac
  (**pending Mac run**: until it has run, the macOS rows above are
  compiled with `GOOS=darwin go vet` only). Differences from Linux: an ACL
  is not carried over by the temp-file replacement (the new file has the
  ACL entries it inherits from its directory; `TestReplaceDropsACLDarwin`);
  EditFile needs write permission on the file (the edit lock is a record
  lock for writing); on a file system without record locks editing is
  refused outright (on Linux the flock edit lock still works there and
  only in-place writes are refused).
- Headers with a card named END that is not the END card, `HIERARCH END
  = ...` included (`ErrMalformedEND`); cfitsio opens some of them (the
  card ends the header there).
- Writing a HIERARCH card whose path is a standard keyword name
  (`HIERARCH BSCALE`, `HIERARCH END`, ...): Flush, `CopyHDUWithHeader`
  and every writer (WriteImage, the Append* writers, on Create files and
  ModeEdit appends) refuse it (cfitsio reads it as that keyword).
- Reading an image whose BSCALE/BZERO/BLANK cards cfitsio and astropy read
  differently (a value that is not a FITS number, a strtod range error, a
  misplaced `=` or lowercase keyword, `HIERARCH BSCALE`/`BZERO`/`BLANK`, a
  BLANK outside the BITPIX range or equal to 1234554321): every read
  returns `ErrInvalidScaling` (`TestScalingKeywordsRule`).
- Stats: `BuildHistogram` on a range numpy cannot bin returns
  `ErrHistogramRange` (numpy raises there); `Percentile`/`Median`/`MAD`
  return NaN when the data hold ±Inf or max−min overflows (Siril's code is
  undefined there).
- Stretch: GHT inverse/modified-asinh/saturation modes, 16-bit Siril paths;
  GHT parameters outside D in [0, 10], B in [-5, 15], 0 <= LP <= SP <= HP <= 1
  or NaN (`ErrGHTParams`); linstretch BP outside [0, 1] or NaN
  (`ErrLinearParams`); autostretch targetBG outside [0, 1] or NaN, or a NaN
  shadows clip (`ErrAutostretchParams`); asinh beta or offset NaN or
  infinite, beta < 0, offset >= 1 (`ErrAsinhParams`); MTF parameters NaN or
  infinite, or midtones outside [0, 1] (`ErrMTFParams`). The accepted
  values outside siril-cli's ranges (asinh beta in [0, 1) and negative
  offsets, MTF midtones 0 and 1 and shadows/highlights outside [0, 1],
  CLAHE clip <= 0 = OpenCV's no limit) are pinned by oracle goldens
  (`TestAsinhEdgeParams_SirilOracle`, `TestMTFEdgeParams_SirilOracle`,
  `TestCLAHECore16bit_cv2`).
- Platforms: tier 1 and tier 2 pass on Linux amd64 (GOAMD64=v1 and v3);
  tier 1 passes on arm64 (`GOARCH=arm64 go test -exec
  qemu-aarch64-static`); tier 1 is meant to pass unchanged on macOS (amd64,
  arm64), where the recorded readings replace the tools. Other architectures are
  not run. Library code wraps every float product in an explicit
  conversion so the compiler cannot fuse it into an FMA (Go spec); the
  frozen `stats.MeanStdev` is the exception and is pinned per build. For
  `stretch` the tier 2 run re-runs the package's tier 1 tests with GOAMD64=v3
  (`TestFMAGuards_GOAMD64v3`, where the compiler fuses) and built for arm64
  under qemu (`TestCrossArch_arm64`), so a dropped FMA guard or a reliance
  on amd64's float-to-int conversions fails on an amd64 host. The
  C/Python truth sources ran on x86-64 only (CLAHE goldens depend on
  OpenCV's SSE2 rounding).

### Known gaps

What is still open after the last verifier round. No item in (a) or (b)
returns silently wrong data for any input the tests know of: each returns
an error, matches the reference, or follows a documented rule. (c) lists
what the library does not try to handle.

**(a) Real-world gaps, each with its current behaviour**

- TNX Legendre (type 2) surfaces: accepted and computed as WCSTools 3.9.7
  computes them, but proven on synthetic headers only (`syn_leg_full`,
  `syn_leg_lonpole_card`, `syn_leg_wat_1999`; the basis against its closed
  forms, `TestTNXBasisDerivatives`). No public real Legendre header was
  found (searched: GitHub code search for TNX WAT strings and `lngcor`,
  NOIRLab Astro Data Archive header scans of 746 headers including a
  by-epoch raw scan 1997-2024 and every remaining instrument name,
  IRAF/MSCRED/PyRAF/astropy/astropy-data/photutils-datasets/ccdproc/
  stsci.tools/starlink-pyast source and data trees, 15 NDWFS / Local Group
  Survey / ChaMPlane HLSP headers, web search; see
  `wcs/transform/testdata/tnx/PROVENANCE.md`).
- TNX with SIP keywords (with or without `-SIP`): refused with an error by
  `wcs.Parse`, `New` and `NewTNX`. This is stricter than wcslib on purpose:
  wcslib 8.10 applies SIP + TAN and ignores the WAT correction when the
  SIP cards come last (and fails when they come first), while WCSTools
  applies TNX; the two references disagree, so no answer is returned
  (`TestHdrStricterThanWcslib`, `TestNewTNXRefusesSIP`).
- No real modern camera or mission `.fz` file (DECam, HSC, ...) is in the
  compressed read corpus. The real compressed files are astropy's M13 (a
  2006 SkyView 300x300 int16 image) in GZIP_1 (4-byte tiles), RICE_1
  (BYTEPIX 4) and PLIO_1 (cut to their first 40 row tiles byte for byte)
  and HCOMPRESS_1 (verbatim), `compressed_with_nan.fits` (fpack 2022) and
  `compressed_image.fits`. A row-tiled modern `.fz` can be cut by whole
  tiles the same way; that has not been done yet.
- HPX and XPH next to their poles keep wcslib 8.10's double-precision
  limit, so the port stays identical to wcslib: up to 2.98e-7 deg forward /
  4.67e-7 deg inverse for HPX within ~1e-6 deg of a pole, 4.08e-7 deg for
  XPH forward at theta = -89.9999997, measured against the exact equations
  (`TestHPXPaperEquations`, `TestXPHPolarPrecision` fail if they grow).
- `HDUByName` matches EXTNAME exactly (case-sensitive, no HDUNAME
  fallback, a non-string EXTNAME is no name), so it can return a later HDU
  than cfitsio's `fits_movnam_hdu`, or none; the HDU it returns always has
  exactly the name asked for (`TestHDUByNameDifferencesFromCfitsio`).

**(b) Not pinned by a dedicated test (code believed correct)**

Mutation testing found these rules where a known-wrong change survives
every committed test. The code follows the reference (or is shown
equivalent on the code) but no sample sits on the deciding input:

- Image writes: cfitsio's 0.49 overflow margin in `roundStored` (moving it
  to 0.48 or 0.5 on either side survives; `write_overflow.json` has no
  value in (hi+0.48, hi+0.5) and the lower side only at lo-0.5), and the
  `(LONGLONG)(d ± 0.5)` rounding at its double-rounding edges
  (d = 0.49999999999999994, odd integers in [2^52, 2^53); `math.Round`
  survives).
- Tables, float32 targets on the exact paths (the float32 nearest the
  exact value, no double rounding): no sample has a uint64 value (K with
  TZERO 2^63, integer-TZERO sums above 2^63), an unscaled ASCII fraction
  or an ASCII integer of more than 18 digits on the float64 path whose
  double lands on a float32 tie.
- ASCII tables on the float64 scaling path: the integer-target refusal of
  a fraction with a tiny exponent (`1.5E-900` under TSCAL 2), and the
  documented float64 arithmetic of float targets (exact rounding of
  decimal*TSCAL+TZERO survives).
- Scaled E/D variable-length arrays (1PE/1PD/1QE/1QD with TSCAL/TZERO):
  applied as for fixed columns, but no committed sample has one.
- Header encoder: the imaginary-part half of the integer-complex guard
  (`integerComplex`; every test decides on the real part).
- HIERARCH aliases of standard keywords that are 8 characters long or
  contain `-` or `_` (`HIERARCH DATE-OBS`, `HIERARCH ZDITHER0`): refused,
  but every test alias is shorter (`hierarchAlias`).
- `Create` on a file system without record locks (truncates and proceeds;
  only the read and in-place write sides are tested with the injected
  `ofdLock`), and a read-only `Open` not being blocked while `Create`'s
  `*File` is open.
- `NewTNX`: the latitude-axis half of the CUNIT refusal (the refused
  sample sets both axes) and the Y half of the empty-surface-domain check.
- Projections: QSC's forward `zeco < 1e-8` switch (`<=` survives; no input
  ties it, `TestQSCZecoSwitchTie`), and the equivalent wcslib guards listed
  in `wcs/transform/projection.go` ("Guards that no input reaches"): the
  SIN/SZP root clamp, `<` against `<=` in SIN's and AIR's forward and
  XPH's inverse (`TestSmallAngleSwitchTies`), prjbchk on outputs that never
  leave the range (`TestPrjbchkRangeByConstruction`), the asind/acosd
  tolerance (`TestWcstrigGuards`), the TSC/QSC face overshoot
  (`TestQuadcubeFaceBound`), QSC's inverse zeta tolerance
  (`TestQSCInverseZetaRange`) and CSC's float32-underflow guards
  (`TestCSCUnderflowGuardPremises`). Measured only: CSC's 1e-7 forward
  overshoot tolerance (wcslib rejects none of 6.5e8 grid points), and
  AIR's `k == 30` return below theta_b ~ 20 deg (an exhaustive scan of
  2.2e6 doubles r found no difference).
- WCS pipeline branches no golden can observe (documented on the code):
  the DSS Newton start point (`dssPoly.Inverse`), in `sphS2X` the exact
  dlng mod 180 = 0 latitude form, the 1e-5 threshold of the rearranged
  longitude formula and the eul[1] >= 90 longitude origin at the native
  pole, in `celPole` the 360 deg sign normalisation and the "native pole
  undetermined" guard.
- Compressed write quantizer: the int64 range check on the zero point
  (an overflowed iqfactor is refused by the per-pixel int32 check anyway),
  the per-pixel guard in the null-holding branches (an overflowed minimum
  decodes back to itself on amd64/arm64), and the ZERO_VALUE half of the
  sentinel check (SUBTRACTIVE_DITHER_2; no reaching input found in a 3
  million sample search).
- Editing: power loss is not simulated (it needs a block-device fault
  injector and root); only the fsync order is tested (`TestFsyncOrder`).
  Keeping another user's ownership on the resize path is proven only with
  an injected EPERM.

**(c) Out of scope**

- AIR inverse at |x|, |y| ≳ 1e6 deg (thousands of turns; no real image
  reaches it): wcslib's regula falsi stops at a 1e-12 tolerance that is
  below one ULP of r, so whether a point is accepted depends on the last
  bit of the log and cosd results, and with Go's libm the decision can
  differ from wcslib's (about 9% of the audit's probe points at the two
  committed theta_b). The committed air_* points at r = 1e6..1e15 deg are
  checked as every other wref row.
- Advisory locking against other programs: the edit and reader/writer
  locks only stop this package. Another program writing the file, or a
  reader of this package that has no file to lock (`OpenReader` of a
  reader that is not an `*os.File`, `OpenReadAll`), can see a torn header
  or data unit, returned without an error; on NFS flock is emulated with
  record locks, where the two locks may interact (not tested).
- `OverwritePixels` of a data unit whose last-block padding is not zero:
  the write replaces the padding with zeros (FITS 4.0 §3.3.2), and its
  crash recovery journals the padding too, but the crash matrix
  (`TestInPlaceWriteCrashAtEveryStep`) uses only files whose padding is
  already zero.

**Upstream and environment observations (not gaps of this library)**

- cfitsio 4.7.0's TDOUBLE `fits_read_img` of a compressed integer image
  with 8-byte GZIP tiles returns its tile buffer unwritten, i.e.
  uninitialized memory with status 0 (zeros or heap addresses depending on
  the heap; the cfref harness runs with glibc's malloc fill so the goldens
  show the 0xA5 pattern; its TLONG read is right; asserted in
  `TestCompressedCorpusRead` and `TestCompressedInt64Read`); funpack
  returns zeros for 8-byte tiles on an 8-bit image. The imgref harness of
  `testdata/images` reads the same files through the same path, so its
  goldens of them hold whatever that process's heap held. Its uint64
  (BZERO 2^63) TLONG read overflows (status 412), and `fits_read_img`
  skips the null check of a tile in an image without ZBLANK.
- wcslib 8.10 crashes (segfault in `dssset`'s error message, whose format
  arguments are out of order) on a DSS header with a non-zero AMDX14..; we
  refuse it (`TestParseDSSRefusals`), so it has no `hdrcases_reject` case.
- cfitsio 4.7.0 bugs worked around in the harnesses: `fits_read_key(TULONGLONG)`
  aborts outside uint64; TLONG `fits_read_imgnull` segfaults on HCOMPRESS
  with BZERO (see `testdata/cfref/PROVENANCE.md`); funpack mis-decodes
  quantized floats with nulls in tiles narrower than NAXIS1; HPXcvt `-xn`
  writes its XPH image rotated 180 deg against its own header.
- The siril-cli build that made the original 18 stretch goldens (git
  `c17580b`) was not recorded and cannot be identified now. The goldens no
  longer rest on it: siril-cli 1.4.0 built from the release tarball
  against a locally unpacked OpenCV 4.13 (the system has OpenCV 5;
  `stretch/testdata/build_siril_cli.sh`) reproduces all 18 digests and all
  146 oracle goldens that are siril commands (`TestSirilCLIReproduce`).
  Its FITS loader rescales float data whose maximum exceeds 10, so ±Inf
  input is given as the markers ±7 and turned into ±Inf by siril-cli's
  pixel math before the command (`gen_golden.py`, checked bit for bit).
- realhdr/make_cases.py, extract_kpno.py and testdata/tables/extract_anet.py
  reproduce their outputs only with the downloaded sources (checked by hand;
  the outputs are sha256-pinned).

### Equivalent mutants: no input can tell them apart (closed by argument)

Mutation testing leaves these stretch mutants alive under every committed
test. Each is shown equivalent here; none is an open gap.

- `rgbblend` (asinh RGB clip mode): `tfmax + d > 1` versus `>= 1` and
  `min(m_CB, candidate)` versus `candidate` give the same k for every input
  with m_CB = 1, the only value Siril uses (the argument is in the
  `rgbblend` comment in `stretch/stretch.go`; a 1.8e7-pixel search found no
  distinguishing input). With m_CB = 1 the blend factor k is always 1
  (`maskCond` is never true and `limited` is always finite), so the output
  is `0*tf + 1*sf` and tf only matters when it is NaN or Inf; so tf = 1
  instead of 0 for a zero channel, tf = 0 instead of 1 for beta = 0 (both
  finite), `d > 0` for `d != 0` and tfmax taken over sf (both only reach
  `limited`, multiplied by 0) are equivalent too.
- `quantize8` (`ToRGBImage`): `v <= 0` for `!(v > 0)` differs only for NaN,
  which then reaches `uint8(NaN)`; the Go spec leaves that conversion
  implementation-specific and it gives 0 on amd64 and arm64, the platforms
  run (`TestCrossArch_arm64`). The code itself does not depend on it.
  `v > 1` for `v >= 1`: 1*255 + 0.5 truncates to 255.
- MAD: `float32(|float64(v) - float64(m)|)` for the float32 `|v - m|`
  (rounding to 53 and then 24 bits equals rounding to 24 once, since
  53 >= 2*24 + 2). The inverted-branch test `float32(median) > 0.5` for
  `median > 0.5` (the float32 histogram median converts exactly).
- asinh mono: k = 1 instead of 0 for x' = 0 (x'*k is 0 or -0 either way).
- CLAHE: the clip test `>=` for `>` (an equal bin is "clipped" to itself
  with no excess); BORDER_REFLECT_101 for a negative index, `-p - 1` for
  `-p` (a negative index needs padding wider than the image, which only
  happens with a tile height or width of 1, and then only in tiles beyond
  the last one any real pixel interpolates from;
  `TestCLAHETiny_SirilOracle` pins the reachable reflections); `p % 1` for
  0 when n = 1. `cvRound` and `cvttsd2si` test NaN and the int32 range in
  one comparison, so there is no separate NaN test to remove; their range
  bounds (`< 2^31` versus `<= 2^31`) are equivalent on amd64, where Go's own
  conversion of 2^31 is MinInt32, and are killed by `TestCrossArch_arm64`.
- Syntactic no-ops: `+ 0`, and `f64(D*HP)` for `f64(float32(D*HP))` with
  float32 operands (the product is already float32).

## Regenerating the goldens

Every generator is deterministic; each harness's `verify` target (or
`*Reproduce*` test) regenerates into a scratch directory and compares with
the committed bytes. Use the venv interpreter for Python
(`venv/bin/python`, Python 3.11.11 with the versions above). C harness
binaries are gitignored build artifacts.

Large text goldens are committed gzipped as `<name>.gz` (gzip level 9,
no name, mtime 0): the wcslib header-pipeline goldens
(`wcs/transform/testdata/wref/golden/hdr`), the cfref pool
(`testdata/cfref/golden/pool`), healpix and testdata/images goldens over
16 KB, the stretch sweeps and `stats/testdata/golden/stats_golden.json`.
The generators write them that way; the tests read them through
`testutil.ReadGolden`, and the reproduce checks compare the decompressed
bytes with the harness output (the compressed bytes depend on the zlib
build and are not compared).

### testdata/cfref (cfitsio dump of every top-level sample)

```sh
make -C testdata/cfref            # build ./cfref
make -C testdata/cfref goldens    # golden/<sample>.json for every testdata/*.fits|fit|std
make -C testdata/cfref verify     # regenerate into verify.out/ and cmp
```

Checked by `TestCfrefDumpReproducesGoldens`. Arrays of more than 1024 values
are written as exact digests (format in `testdata/cfref/PROVENANCE.md`); the
same harness output is used by `testdata/images`, `testdata/tables`,
`testdata/comp_bad` and `header/testdata`, so after a harness change
regenerate those too. Details: `testdata/cfref/PROVENANCE.md`.

### Compressed read corpus (testdata/comp_*)

```sh
FITS_ASTROPY_PYTHON=$PWD/venv/bin/python sh testdata/gen_fpack.sh   # comp_fz_*, comp_fzsrc_*, comp_bad/*, comp_fzun_same.txt
venv/bin/python testdata/generate_compressed.py   # comp_*, comp_raw_*
FITS_ASTROPY_PYTHON=$PWD/venv/bin/python make -C testdata/cfref goldens   # cfitsio truth for every sample, the twins (gen_twins.sh) and comp_bad
make -C testdata/images all PYTHON=$PWD/venv/bin/python   # imgref pix/chksum goldens of the new samples
make -C header/testdata goldens                   # header goldens of the new samples
```

`TestCompressedCorpusReproduce` re-runs both generators in a scratch copy
and requires every committed `comp_*.fits` and `comp_bad/*.fits` byte for
byte; `TestCompressedGoldensReproduce` regenerates the uncommitted twins
and requires every golden (samples, twins, comp_bad) and pool file.

Details: `testdata/COMP_FZ_PROVENANCE.md`.

### testdata/images (image read/write truth corpus)

```sh
make -C testdata/images                                     # build imgref and cfref
make -C testdata/images all PYTHON=$PWD/venv/bin/python     # every sample + golden
make -C testdata/images verify PYTHON=$PWD/venv/bin/python  # regenerate, diff -r
```

Details: `testdata/images/PROVENANCE.md`.

### testdata/imaging (ReadFloat32 / frame / AsImage)

```sh
cd testdata/imaging
../../venv/bin/python gen_samples.py      # every *.fits except the camera crop
./make_compressed.sh                      # comp_*.fits.fz (fpack -C + zblank_variants.py)
make -C cfnull goldens                    # cfitsio null flags / physical values
../../venv/bin/python gen_expected.py     # readfloat32_expected.json, stdimage_expected.json
make -C sirilref goldens                  # Siril convention goldens
make -C cfnull verify; make -C sirilref verify; make -C sirilref verbatim
```

The camera crop needs the source frame (`crop_camera.py <Light_001.fits>`,
not committed). Details: `testdata/imaging/PROVENANCE.md`.

### testdata/tables

```sh
make -C testdata/cfref              # cfitsio harness
make -C testdata/tables samples     # generated samples (gen_tables.c, gen_theap.py, gen_gaia.py)
python3 testdata/tables/extract_anet.py --src <downloads>   # anet_*.fits from the real astrometry.net index files
make -C testdata/tables goldens     # golden/*.json (compact), refused/*.json
make -C testdata/tables verify
```

Checked by `TestTableGoldensReproduce` and `TestTableSamplesReproduce`.
The cfitsio goldens are stored without whitespace between tokens
(`testdata/tables/compact_json.py`, every token byte for byte, which
`json.Compact` reproduces). Details: `testdata/tables/PROVENANCE.md`.

### testdata/cfcopy (fits_copy_hdu goldens for EditFile)

```sh
make -C testdata/cfcopy           # build cfcopy and mkcomments
make -C testdata/cfcopy goldens   # reserve.fits, comments.fits, golden/*.fits
make -C testdata/cfcopy verify
```

Checked by `TestCfcopyReproducesGoldens`. Details: `testdata/cfcopy/PROVENANCE.md`.

### testdata/cfmovnam (fits_movnam_hdu golden for HDUByName)

```sh
make -C testdata/cfmovnam           # build cfmovnam, mkdupext and mknames
make -C testdata/cfmovnam goldens   # dupext.fits, names.fits, golden/{dupext,names}.json
make -C testdata/cfmovnam verify
```

Checked by `TestCfmovnamReproducesGoldens`. Details: `testdata/cfmovnam/PROVENANCE.md`.

### header/testdata (hdrref + cfref)

```sh
make -C header/testdata hdrref/hdrref
make -C header/testdata samples
make -C header/testdata goldens
make -C header/testdata verify PYTHON=venv/bin/python
```

Checked by `TestGoldensReproduce` (header). Details: `header/testdata/PROVENANCE.md`.

### compress/testdata/cref (quantizer, encoders, whole-image writer goldens)

```sh
make -C compress/testdata/cref           # build ./cref
make -C compress/testdata/cref goldens
make -C compress/testdata/cref verify
```

Details: `compress/testdata/cref/PROVENANCE.md`.

### wcs/transform/testdata/wref (wcslib 8.10)

```sh
make -C wcs/transform/testdata/wref                   # wref and hdr2sky
venv/bin/python wcs/transform/testdata/wref/gen_prj_cases.py     # hdrcases/prj_*, pvfid_*, prjswap_*
venv/bin/python wcs/transform/testdata/wref/gen_sample_cases.py  # sample_* cases + samples/ (needs HPXcvt, healpy)
python3 wcs/transform/testdata/wref/realhdr/make_cases.py --src <downloads>  # real_* cases
make -C wcs/transform/testdata/wref goldens           # golden/*.json + *.bin and golden/hdr/*.json.gz
make -C wcs/transform/testdata/wref verify            # regenerate all of them into verify.out/ and cmp
venv/bin/python wcs/transform/testdata/wref/gen_libm.py       # libm/trig.json, libm/atan2.json + .bin (mpmath)
venv/bin/python wcs/transform/testdata/wref/gen_xph_paper.py  # paper/xph_paper.json + .bin
venv/bin/python wcs/transform/testdata/wref/gen_hpx_paper.py  # paper/hpx_paper.json + .bin
```

Checked by `TestWrefGoldensReproduce`, `TestHdrGoldensReproduce`,
`TestPrjCaseGeneratorReproduce`, `TestSampleGeneratorsReproduce`,
`TestLibmGeneratorReproduces`, `TestXPHPaperGeneratorReproduces`. The
header goldens (format `hdr2sky-2`, 1.94 MB for 208 cases, largest 48 kB) store wcslib's
per-point status and output coordinates as exact binary64 records and
sha256 digests of the case text, cards and input points, which the test
rebuilds from the case file. Details:
`wcs/transform/testdata/wref/PROVENANCE.md`.

### WCS astropy crosscheck

```sh
venv/bin/python testdata/generate_wcs_crosscheck.py > testdata/wcs_crosscheck.json
```

Checked by `TestCrossCheckFixturesReproduce`.

The deliberate-difference cases (`wcs/transform/testdata/wref/hdrcases_differ`)
have no golden of their own: astropy's reading (`astropy_hdr.py`) and
wcslib's (`hdr2sky`) of each header as written are recorded tool readings,
re-recorded with

```sh
FITS_RECORD=1 FITS_ASTROPY_PYTHON=$PWD/venv/bin/python go test -tags interop -count=1 \
  -run '^(TestHdrDeliberateDifferences|TestATLASTPVPlusRedundantSIP)$' ./wcs/transform
```

### wcs/transform/testdata/frames (astropy / ERFA)

```sh
venv/bin/python wcs/transform/testdata/frames/gen_frames.py > wcs/transform/testdata/frames/frames_golden.json
```

Checked by `TestFramesGoldenReproduce`. Details: `wcs/transform/testdata/frames/PROVENANCE.md`.

### wcs/transform/testdata/tnx (WCSTools 3.9.7)

```sh
cd wcs/transform/testdata/tnx
python3 make_synthetic.py           # headers/syn_*.hdr, refused/*.hdr
python3 gen_golden.py               # golden/*.json, refused/*.json via xy2sky / sky2xy
python3 gen_golden.py --check
python3 extract_kpno.py <downloaded .fits.fz>   # headers/kpno_kp4m_mosaic1_ccd1.hdr (needs funpack)
```

Checked by `TestTNXGoldensReproduce`. Details: `wcs/transform/testdata/tnx/PROVENANCE.md`.

### healpix/testdata (healpy + healpix_cxx)

```sh
cd healpix/testdata
make -C hpref                                   # healpix_cxx 8.3.0 harness
../../venv/bin/python gen_golden.py golden      # golden/*.json
cd maps && ../../../venv/bin/python gen_map.py  # healpy_map_nside16_*
```

Details: `healpix/testdata/PROVENANCE.md`.

### stats/testdata (numpy/astropy + Siril)

```sh
make -C stats/testdata/sref goldens             # Siril oracle golden
venv/bin/python stats/testdata/gen_golden.py    # golden/stats_golden.json.gz (and the histogram/mean goldens)
make -C stats/testdata/sref check-excerpt       # Siril excerpt is verbatim
```

Checked by `TestGoldensReproduce` (stats). Details: `stats/testdata/PROVENANCE.md`,
`stats/testdata/sref/PROVENANCE.md`.

### stretch/testdata (Siril 1.4.0 + OpenCV)

```sh
make -C stretch/testdata/sref                    # sref, sref_libm and sref_clahe
make -C stretch/testdata/sref goldens            # inputs/, stretch_*.json, sref/* via gen_siril_goldens.py
make -C stretch/testdata/sref verify             # regenerate into a scratch dir, byte-compare
make -C stretch/testdata/sref verify-excerpt
venv/bin/python stretch/testdata/gen_clahe_cv2.py   # cv2 goldens
venv/bin/python stretch/testdata/siril_cli_digests.py --check   # siril-cli digests from git history
stretch/testdata/build_siril_cli.sh                            # siril-cli 1.4.0 under /var/lib/aurbuild/pkgbuild/siril
venv/bin/python stretch/testdata/gen_golden.py <siril-cli>      # check presets, sref goldens, argument table
venv/bin/python stretch/testdata/gen_golden.py <siril-cli> --write-args   # rewrite golden/siril_cli_args.json
```

The 18 preset goldens (`golden/stretch_*.json`) are written by
`gen_siril_goldens.py` from the verbatim Siril 1.4.0 oracle, which refuses
to write them unless each output digest equals the digest of the original
siril-cli pixel file (`golden/siril_cli_sha256.json`, recomputed from git
history by `siril_cli_digests.py`). The siril-cli build that made those
files was not recorded; `gen_golden.py` checks a siril-cli 1.4.0 built by
`build_siril_cli.sh` against the digests, the oracle goldens and
`golden/siril_cli_args.json` (`TestSirilCLIReproduce`). Details:
`stretch/testdata/PROVENANCE.md`.
