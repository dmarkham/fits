# fits

Pure-Go FITS library for reading, writing, and processing astronomical
data files.

## Philosophy

This library is not trying to replace cfitsio, wcslib, healpy, or any
other established C/C++/Python implementation. Those tools have decades
of production mileage and are the right choice for many workflows. Our
goal is to bring **Go onto level footing** — giving Go programs native
access to the same FITS capabilities without cgo, while staying as
compatible as possible with the files and conventions those reference
implementations define.

Every algorithm is ported from the reference source (cfitsio, wcslib,
astrometry.net, Siril) and cross-validated against it. Where our output
differs from the reference, we treat that as a bug in our code, not a
difference of opinion.

## Guiding principle

This library is 100% scientific fidelity. It returns the bytes, values, and
types exactly as the file stores them. It does not auto-scale, auto-window,
auto-percentile, auto-display, silently coerce types, or silently map
NaN/BLANK to zero. Any interpretive helper that rescales, reshapes, or
transforms data for display belongs in a sibling package or a separate
module — never in the root `fits` package.

## Goals

- Pure Go, no cgo, standard library only.
- Idiomatic Go: `io.ReadSeeker`, `error` returns, method-based API, generics
  for typed pixel and column access.
- `image.Image` compatibility for the narrow cases where it is lossless.
- Context support on blocking I/O.

## Capabilities

- **Read/write** primary images and IMAGE extensions with every standard
  BITPIX (±8/16/32/64 int, -32/-64 float). BSCALE/BZERO scaling applied
  transparently.
- **Binary tables** with every scalar TFORM type, vector (repeat>1)
  columns, and variable-length arrays (P/Q descriptors + heap).
- **ASCII tables** (read).
- **Header** with full FITS v4 keyword compliance: fixed/free format,
  CONTINUE long strings, HIERARCH non-standard, all value types.
- **In-place edit surface**: header mutation with journaled tail-shift
  crash safety, same-shape pixel overwrite, HDU append at EOF.
- **Streaming rebuild edit surface** (`EditFile`): structural changes via
  temp file + atomic rename.
- **image.Image adapter** for 2D integer BITPIX + `FloatGray` for float
  BITPIX (no auto-scaling — see guiding principle).
- **World Coordinate System** (`fits/wcs` + `fits/wcs/transform`): all 27
  projections from Paper II + HEALPix, SIP/TPV/TNX distortion, sky-frame
  conversions (ICRS/FK5/FK4/galactic/ecliptic/supergalactic),
  cross-validated against wcslib/astropy.
- **Tile compression** (`fits/compress`): both read AND write for every
  algorithm in the Pence et al. 2010 FITS tile compression convention —
  RICE_1, GZIP_1, GZIP_2, HCOMPRESS_1, PLIO_1, and NOCOMPRESS. Read-side
  cross-validated byte-exact against astropy-generated fixtures; write-
  side cross-validated by having astropy read back Go-written files
  (all 6 algorithms pass).
- **Float tile compression** with per-tile quantization (RICE_1, GZIP_1,
  GZIP_2) via the cfitsio `fits_quantize_float` algorithm ported
  bit-exactly to Go and validated against live libcfitsio on 11 golden
  fixtures. Supports NO_DITHER, SUBTRACTIVE_DITHER_1, and
  SUBTRACTIVE_DITHER_2 via the IRAF Park-Miller PRNG. NaN handling
  (ZBLANK) and GZIP_COMPRESSED_DATA fallback for constant tiles both
  cross-validated against astropy. Covers JWST, HST, Gaia, Planck,
  WMAP, and every other mission using tile-compressed FITS.

## Requirements

- Go 1.23 or later (generics + `iter.Seq2`).
- Linux, macOS, and Windows. WASM is not supported (the edit surface
  requires a real filesystem; read-only via `OpenReader` / `OpenReadAll`
  may work but is untested).

## Concurrency

`*File` is **not safe for concurrent use by multiple goroutines.** Callers
that want parallelism must coordinate externally (one `*File` per goroutine,
or an external mutex).

## Cross-validation

The library is validated against two independent reference implementations:

- **cfitsio** (C, HEASARC). Byte-for-byte round-trip against every fixture
  in cfitsio's own test suite. Float quantization ported with a C reference
  harness (`compress/testdata/cref/`) that links libcfitsio and diffs
  results against the Go port; 11/11 fixtures agree bit-exactly on noise
  estimation, quantized int32 output, and bscale/bzero.
- **astropy / wcslib** (Python). Every one of the 27 WCS projections and
  all 6 sky frames round-trip within numerical tolerance against live
  astropy.wcs. All 6 tile compression algorithms survive a Go → astropy
  round trip, including float quantization with SUBTRACTIVE_DITHER_2
  exact-zero preservation and the GZIP_COMPRESSED_DATA fallback column.

See `plan.md` for the full design and scope.

## License

MIT. See `LICENSE`.
