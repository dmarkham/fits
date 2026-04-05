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

## Requirements

- Go 1.23 or later (generics + `iter.Seq2`).
- Linux. macOS, Windows, and WASM are not supported in v1.

## Concurrency

`*File` is **not safe for concurrent use by multiple goroutines.** Callers
that want parallelism must coordinate externally (one `*File` per goroutine,
or an external mutex).

## License

MIT. See `LICENSE`.
