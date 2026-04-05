# fits

Pure-Go port of the CFITSIO library for reading and writing FITS files.

**Status:** in development. See `plan.md` for the design and scope.

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
- **Tile compression** (`fits/compress`): read-side decoders for RICE_1,
  GZIP_1, GZIP_2, and NOCOMPRESS algorithms — covers JWST, HST, and most
  modern mission data products. HCOMPRESS_1 and PLIO_1 are recognized
  but return `ErrUnsupportedCompression` until their decoders are wired
  in a follow-up pass.

## Requirements

- Go 1.23 or later (generics + `iter.Seq2`).
- Linux. macOS, Windows, and WASM are not supported in v1.

## Concurrency

`*File` is **not safe for concurrent use by multiple goroutines.** Callers
that want parallelism must coordinate externally (one `*File` per goroutine,
or an external mutex).

## License

MIT. See `LICENSE`.
