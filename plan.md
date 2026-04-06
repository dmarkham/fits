# FITS Go Library — Plan

## Context

Pure-Go FITS library. Target module path: `github.com/dmarkham/fits`. Reference materials in `reference/`: FITS standard v3.0 and v4.0 PDFs plus a full clone of HEASARC/cfitsio.

## Philosophy — compatibility, not replacement

This library is not trying to replace cfitsio, wcslib, healpy, Siril, or any other established C/C++/Python tool. Those projects have decades of production mileage across thousands of missions and pipelines. They are the right choice for many workflows and will remain so.

Our goal is to bring **Go onto level footing** with the existing ecosystem — giving Go programs the same FITS read/write, WCS, compression, HEALPix, and image-processing capabilities without requiring cgo, while staying as compatible as possible with the files and conventions the reference implementations define.

Every algorithm is ported from a specific reference source (cfitsio for I/O and compression, wcslib for projections, astrometry.net for HEALPix, Siril/RawTherapee for image stretching) and cross-validated against it with golden fixtures. Where our output differs from the reference, we treat that as a bug in our code, not a difference of opinion.

## Guiding principle — no auto-anything

This library is 100% scientific fidelity. It returns the bytes, values, and types exactly as the file stores them. It does not auto-scale, auto-window, auto-percentile, auto-display, silently coerce types, or silently map NaN/BLANK to zero. Any interpretive helper that rescales, reshapes, or transforms data for display belongs in a sibling package (e.g. `fits/display`) or a separate module — never in the root `fits` package. When in doubt, fail loudly at the boundary rather than guess.

## Goals

1. **Pure Go.** No cgo. Standard library only where possible; third-party dependencies only if absolutely necessary (none currently planned).
2. **Idiomatic Go.** `io.ReadSeeker` not file pointers, `error` not status codes, method-based not `ffxxxx(fptr, ...)`-style, generics for typed pixel/column access.
3. **`image.Image` compatibility** for the narrow cases where it's lossless (2D integer BITPIX, including the TZERO unsigned trick). A native-fidelity `FloatGray` type for float BITPIX. Multi-dim and anything else returns `ErrNotImageCompatible`.
4. **Context support** where blocking I/O happens.

## Project setup

- **Module path**: `github.com/dmarkham/fits`
- **Repository visibility**: private
- **License**: MIT
- **Default branch**: `master` (per CLAUDE.md convention)
- **Minimum Go version**: 1.23 (generics + `iter.Seq2`)
- **Platform support**: **Linux only for v1.** macOS / Windows / WASM are not tested or supported in v1. The `EditFile` atomic-rename and fsync-on-directory paths assume POSIX semantics; we do not ship platform-specific code paths in v1. Non-Linux builds may work accidentally for read-only paths but are not a commitment.
- **Concurrency model**: `*File` is **not safe for concurrent use by multiple goroutines.** Callers that want parallelism must coordinate externally (one `*File` per goroutine, or an external mutex). The library does not take internal locks in v1. Adding goroutine safety later is a non-breaking enhancement, so we keep the option open but do not commit to it. Document this prominently in the package godoc.
- **Interop cross-check — MANDATORY in CI, not optional.** Every file the library writes must round-trip through cfitsio's `fitsverify` with zero errors. Every test fixture must parse-equivalent through `astropy.io.fits`. We do not ship v1 without both passing.
  - Both toolchains are already installed and in active use by the neighboring `/home/dmarkham/git/dmarkham/fits-manager` project (Python + `astropy==7.2.0` + `fitsio==1.3.0` + system `cfitsio`). Our CI harness uses the same setup — either reuse the fits-manager venv or mirror its `requirements.txt` pins.
- **Prior art (surveyed April 2026)**:
  - [`github.com/astrogo/fitsio`](https://github.com/astrogo/fitsio) — was the dominant Go FITS library. **Archived on GitHub March 2025**, moved to codeberg. Last stable release v0.1.0 November 2019. BSD-3. Pure Go, read/write images + tables + VLAs. No generics (predates Go 1.18), no compression, no WCS, no `image.Image`, no edit surface. API pitfall to avoid: `Rows.Scan(args ...interface{})` — type-unsafe, defeats compile-time checking; our generics replace this directly.
  - [`github.com/ATTron/fits`](https://pkg.go.dev/github.com/ATTron/fits) — active exploratory project framed as "assess suitability of golang for scientific applications." Reads images + tables, not aiming at production.
  - [`github.com/observerly/iris`](https://github.com/observerly/iris) — specialized for FITS + ASCOM Alpaca observatory automation, not a general FITS library.
  - **Conclusion**: the Go FITS ecosystem has no actively-maintained comprehensive library as of April 2026. astrogo was the dominant player and has been archived for over a year. Our port fills a genuine vacuum. Gap-fillers over the existing set: generics-first typed access, full header compliance (CONTINUE, HIERARCH), `image.Image` interop, in-place + rebuild edit surface, committed post-v1 WCS math, and active maintenance. Clean-room implementation from the FITS v4 spec — we do not copy code from any existing library.

## Scope

### In scope for v1

- File I/O over `io.ReadSeeker` / `io.ReadWriteSeeker`; `os.File` and in-memory drivers
- 2880-byte block framing (§3.3.1, §7.3.3)
- HDU discovery and navigation; eager header parsing, lazy data decode
- Primary HDU and IMAGE extensions (§7.1), read and write, all mandatory BITPIX values (±8/16/32/64 int, -32/-64 float)
- Full header / keyword support per §4 and Appendix A:
  - Fixed- and free-format cards, all value types (string, logical, int, float, complex)
  - `CONTINUE` long-string convention (§4.2.1.2)
  - `COMMENT`, `HISTORY`, blank, `END`
  - `HIERARCH` non-standard (widely deployed, supported by cfitsio)
  - Mandatory keyword ordering and validation (Tables 7, 10, 14, 17)
- Big-endian byte handling (§5.2)
- BSCALE/BZERO/BLANK read-side application and inverse on write (including the TZERO unsigned-int trick in columns, §7.3.2 / Table 19)
- ASCII tables (§7.2) read and write: `A`, `I`, `F`, `E`, `D` TFORM codes, TBCOL offsets, TNULL
- Binary tables (§7.3) read and write:
  - All scalar TFORM types (`L`, `X`, `B`, `I`, `J`, `K`, `A`, `E`, `D`, `C`, `M`)
  - Repeat counts (vector columns)
  - TDIMn multi-dim cell shape
  - Unsigned integer convention via TZERO
  - **Variable-length arrays** (`P`/`Q` descriptors + heap, §7.3.5)
- CHECKSUM and DATASUM (Appendix J) compute and verify
- `image.Image` adapter + `image.RegisterFormat` for the cases it can handle losslessly
- Round-trip fidelity: read → write must be byte-identical for unmodified files (modulo recomputed checksums)

### Deferred (not v1, explicitly)

- Tiled image/table compression (§10) — the `imcompress.c` / `quantize.c` / `ricecomp.c` / `fits_hcompress.c` bundle is ~15 k LOC. v1 detects compressed HDUs and returns a clear "not supported" error. Ships later as `fits/compress`.
- Random groups (§6) — deprecated, radio-interferometry only.
- **WCS — World Coordinate System (§8).** Not in v1. **Committed as the next major development target after v1 ships.** Two-package split for architectural cleanliness, shipping together as one "full WCS" deliverable:
  - `fits/wcs` — metadata parsing. Parses `CTYPEn`, `CRVALn`, `CRPIXn`, `CDELTn`, `CUNITn`, `PCi_j` / `CDi_j`, `PVi_m` / `PSi_m`, `LONPOLE`, `LATPOLE`, `RADESYS`, `EQUINOX`, etc. into a typed `*WCSHeader` struct. ~200 LOC, no math. Consumes `*header.Header` as input. Validates keyword combinations (e.g. mutually exclusive CDELT+PC vs CD matrix forms).
  - `fits/wcs/transform` — full projection math. **The bulk of the post-v1 work.** Implements all 27 spherical projections (TAN, SIN, ARC, ZEA, STG, SZP, ZPN, AZP, CYP, CEA, CAR, MER, SFL, PAR, MOL, AIT, COP, COE, COD, COO, BON, PCO, TSC, CSC, QSC, HPX, XPH) — forward and inverse — plus SIP/TPV/TNX distortion conventions and sky frame transforms (equatorial ↔ galactic ↔ ecliptic ↔ supergalactic with epoch/equinox handling). ~2000–3000 LOC of dense spherical trigonometry. Consumes `*WCSHeader`. Tests cross-validated against `wcslib` (the C reference implementation) on a fixture set of real observatory data (HST, JWST, Chandra, Gaia).
  - Rationale for split (not for deferral): parsing is a small separable layer with a stable keyword vocabulary; math is a large body of projection code with its own correctness story and test matrix. Keeping them in separate packages matches cfitsio→wcslib, astropy.io.fits→astropy.wcs, FITSIO.jl→WCS.jl and lets the math package evolve on its own release cadence without forcing core library bumps.
  - **Explicit non-goal of the post-v1 WCS work**: half-implemented projection coverage. The full 27-projection set ships or nothing does. Shipping a subset (e.g. "TAN and SIN only") is a correctness trap — users get silently wrong coordinates for images using unsupported projections and don't find out until their pipeline is wrong.
- **Time coordinates (§9).** Not in v1. Planned as two separate future packages:
  - `fits/time` — metadata parsing only. Parses `TIMESYS`, `DATE-OBS`, `MJD-OBS`, `JDREF`, `MJDREF`, `DATEREF`, `TREFPOS`, `TREFDIR`, `TIMEUNIT`, `TSTART`, `TSTOP`, `TELAPSE`, `TIMEDEL`, `TIMEPIXR`, `TIMEOFFS` into a typed `*TimeHeader` struct. ~100 LOC, no time math. Convenience `ObsTime() (time.Time, error)` method documented with its limitations (ignores leap seconds, assumes UTC when TIMESYS is empty).
  - `fits/time/convert` — rigorous timescale conversion. UTC↔TAI (leap-second table, requires periodic updates against IERS bulletins), TAI↔TT (fixed offset), TT↔TDB (1658-term relativistic series), TCG/TCB. May instead delegate to an existing Go time library (`github.com/soniakeys/meeus`, `github.com/hablullah/go-juliandays`) rather than re-implement. Decision deferred until the package is scoped.
  - Rationale for split: same as WCS. Leap-second tables and relativistic series expansions are their own specialty with their own maintenance cycle. Bundling them into the I/O library couples release schedules and pulls every user into the math dependency whether they need it or not.
- Row-filter / column expression parser (`eval.l` / `eval.y`, ~15 k LOC)
- Iterator API (`iter_*.c`) — Go users will use `range` / channels
- Hierarchical grouping (`group.c`, `grparser.c`)
- Region parsing (`region.c`)
- Extended filename syntax (`file.fits[1:10,*]`)
- Network drivers (HTTP / FTP / gsiftp) — users layer on via `io.ReadSeeker`
- **IRAF format** (`iraffits.c`, ~1800 LOC in cfitsio). IRAF's `.imh`/`.pix` pair is a non-FITS legacy format from NOAO. cfitsio reads it via a synthetic-FITS-cards translator for backward compatibility with pre-2010 astronomy pipelines. IRAF itself was retired by NOIRLab in 2013. **No future package planned** — if a user ever needs it, they can write a standalone `iraffits → fits` byte-stream translator and feed it to `OpenReader`; no core changes needed.
- **Histogram / binning** (`histo.c`, ~2500 LOC in cfitsio). cfitsio's server-side histogramming turns binary table event lists into images via a filter-expression DSL coupled to `eval.l`/`eval.y`. **No future package planned in the core module** — this is an analysis concern, not an I/O concern. Users write the binning loop directly against `ReadColumn` output (typically ~20 lines). If the pattern becomes common, it lives in a downstream analysis library, not here.

### Explicit non-goals

- Bug-for-bug compatibility with cfitsio quirks. Behavior-for-behavior on the spec: yes. Reproducing cfitsio internal state semantics: no.
- Auto-scaling of any kind in the main library (see guiding principle).
- **Fortran bindings.** cfitsio ships ~8000 LOC of `f77_wrap*.c` + `cfortran.h` glue to expose its C API to Fortran 77 callers (name mangling, string convention translation, 1-based ↔ 0-based indexing, `LOGICAL*1` ↔ `int`). Legacy interop exists in cfitsio because NASA HEASARC has major Fortran components (HEASOFT, AIPS, Miriad, XSPEC legacy parts). **Never in this library.** We are Go-native. Fortran interop would require cgo at minimum, violating the "no cgo" goal, and the target audience for our library is Go programs — Fortran users already have cfitsio and don't need us. Mixed codebases can use both libraries independently on the same files on disk. This is not "not yet" — this is "not ever, by design."

## Design decisions

### Approved

**5.1 Generics-only API.** One `ReadPixels[T Numeric]`, one `ReadColumn[T]`, one `ReadVarColumn[T]`. Deletes cfitsio's 12-file combinatorial explosion from `getcol*.c` / `putcol*.c`. Compiler monomorphizes per `T`, zero runtime cost. Honest about the fact that the type parameter is the *output* type, not an assertion about file BITPIX. Runtime mismatch surfaces as `ErrTypeMismatch{Requested, BITPIX, Lossy}` with an optional `CanConvert[T](hdu)` pre-flight helper.

**5.4 Strict `image.Image` mapping.** `AsImage()` returns a native `image.Image` **only** for 2D integer BITPIX (including the TZERO unsigned trick → `image.Gray16`). Float BITPIX returns a `FloatGray` that implements `image.Image` but holds raw `float32` — faithful to the values, not rescaled. Multi-dim returns `ErrNotImageCompatible`; users call `hdu.Slice(...)` explicitly to pin higher axes. `image.Decode` registered via `image.RegisterFormat("fits", "SIMPLE  =", ...)` but returns a clean error on anything it can't losslessly represent. **No `AsImageScaled`, no auto-zscale, no auto-percentile.** Any display helpers live in a sibling package.

**5.5 `VarColumn[T]` struct with flat storage + zero-copy sub-slice views.** Variable-length columns read into one contiguous `[]T` slab plus a `[]int64` offsets array (len = nrows+1). `At(row) []T` returns a zero-copy sub-slice. `Rows()` returns an `iter.Seq2[int, []T]` for range iteration. `Raw()` exposes `(values, offsets)` for GPU / SIMD / bulk consumers. Two allocations for 10 M rows instead of 10 M+1. Mirrors the on-disk FITS heap layout faithfully.

**5.2 Minimal eager, lazy full parse.** On `Open`, walk the file and for each HDU parse *only* the handful of structural keywords needed to compute data size and answer hot questions: `XTENSION`/`SIMPLE`, `BITPIX`, `NAXIS`, `NAXISn`, `PCOUNT`, `GCOUNT`, `EXTNAME`, `EXTVER`. Retain the raw header bytes `[]byte` per HDU. Parse the full `[]Card` lazily on first `Header()` touch, cache the result behind `sync.Once`. TFORM parsing for binary tables also happens on first `Columns()` / `ReadColumn` access.
- Rationale: open time drops ~15× for large MEF files (JWST/HST with 100+ HDUs), memory per open file drops ~3× vs full-eager parsed `[]Card`, and `NumHDU`/`HDUByName`/`Type`/`Shape`/`BITPIX` stay O(1) because the hot metadata was extracted during the Open scan. Full-lazy (option C) was rejected because `HDUByName` forces a per-HDU keyword parse anyway, collapsing it into this approach with worse ergonomics.
- Implementation: `hduRecord` struct holds `{offsets, bitpix, naxis, shape, pcount, gcount, extname, extver, rawHeader []byte, parsed *header.Header, once sync.Once}`. Public API is unchanged.

**5.3 `io.ReadSeeker` baseline, `OpenReadAll` slurp helper for streams, no streaming mode in v1.** Every internal read path assumes seeking is free. Users with an `io.Reader` call `OpenReadAll(r)` which buffers the stream into memory and hands back a regular `*File`. No auto-buffering, no temp-file spill, no dual seek/stream code paths — the library has one I/O model.
- Rationale: doubling the internal implementation to support a streaming mode buys nothing at v1 — real large FITS files live on POSIX/Lustre filesystems (ALMA, JWST, radio interferometer data) and come with file paths. Stream-sourced FITS is always either small enough to slurp or network-sourced from archives with HTTP range-request support (wrap as seeker externally). The "50 GB from stdin" case has no real-world example.
- Reserved for later (not v1): `NewStreamReader(r io.Reader) *StreamFile` for a reduced-API streaming mode (forward-sequential only, no VLA/subset/random access, feature gaps enforced at the type level). Purely additive, won't break v1 API if it ships in v2.
- Not provided by this library: HTTP range-request seeker adapters. Users bring their own or use third-party helpers. Revisit if the pattern is common enough.

**5.6 `internal/block` stays internal for v1.** Raw 2880-byte blocks are not part of the public API. Our own `cmd/fitsdump`, `cmd/fitscopy`, and internal tests use `internal/block` freely because they sit in the same module; external consumers cannot import it.
- Rationale: the block implementation is going to evolve (read cache, possible `io.ReaderAt` switch for concurrent reads, possible memory-mapping on Linux, prefetch strategy changes). Every one of those becomes a breaking change if block is public. No concrete external consumer exists at v1, so there is no use case to design the exposure shape around.
- Escape hatches that cover the "I need the raw bytes" story without exposing blocks: `Card.Raw [80]byte` preserves original card bytes for round-tripping; `(*BinaryTableHDU).RowBytes(row int64) ([]byte, error)` exposes raw row bytes; the minimal-eager `hduRecord.rawHeader []byte` is already retained internally and can be surfaced via a method trivially if ever needed. Document these prominently so users see the higher-level answer first.
- Promotion path if demand appears in v2: add curated methods on `*File` (e.g. `(*File).Blocks() iter.Seq2[int64, [2880]byte]`, `(*File).ReadBlock(n int64) ([2880]byte, error)`). **Do not promote the package itself to `fits/block`** — method-on-File keeps the implementation shape free to change. Purely additive, no v1 breakage.

**5.8 Header-in-place + same-shape pixel overwrite + `EditFile` streaming rebuild for structural edits.** Two distinct edit surfaces:

**In-place edit surface (`OpenForEdit`)** — fast, for operations where in-place is meaningfully cheaper than a rebuild:
- Keyword update / delete / add when the header still fits in its current 2880-byte block allocation: rewrite the affected block(s) in place, fsync. O(KB).
- Keyword add that pushes the header past its current block allocation: shift the tail of the file forward by N × 2880 bytes using a journal file in the same directory (back-to-front tail copy, journal records the shift, recovery on next open replays or rolls back). O(file size) in I/O, but single-pass and in-place — no doubled disk footprint. Crash-safe by journal.
- Same-shape pixel overwrite (`(*ImageHDU).OverwritePixels[T]`): seek to `dataStart`, write bytes, apply inverse BSCALE/BZERO. Errors if `T` / BITPIX / NAXIS mismatch. O(data size).
- HDU append at EOF (`(*File).AppendImage` / `AppendBinaryTable` / `AppendASCIITable` in edit mode): fast, no shift needed.

**Structural edit surface (`EditFile`)** — for operations that cfitsio also does in O(file size) anyway:
- Image resize (change BITPIX or NAXISn), row insert/delete, column insert/delete, HDU insert in middle, HDU delete, VLA heap reorganization.
- `fits.EditFile(path, fn)` opens source read-only, creates a temp file in the same directory, passes the user a `*File` + `*Writer` to stream HDUs through (with `CopyHDU` / `CopyHDUWithHeader` / `SkipHDU` / `AppendImage` / ...), fsync + atomic `os.Rename` over source on success.
- Crash-safe by construction (temp file + rename). Callback can mix copy and append calls freely.

**HDU data model stays immutable.** Headers returned from read-mode `Open()` can be mutated locally (they're Go in-memory objects) but `Flush()` errors because the file is read-only. Headers returned from `OpenForEdit()` mutate locally and persist on `Flush()` / `Close()`.

- Rationale: option (B) from the design discussion. Covers the ~90% case (keyword edits, checksum recompute, pixel overwrite, HDU append) at cfitsio speeds without the ~5800 LOC / bug-prone `edithdu.c` + `editcol.c` buffer-shuffling nightmare. The operations punted to `EditFile` rebuild are the ones where cfitsio's in-place implementation is only marginally faster than a rebuild anyway (all involve O(file size) shifts).
- Implementation cost: ~800 LOC for in-place paths (header-in-block rewrite, boundary-shift with journal, same-shape pixel overwrite, append) + ~400 LOC for `EditFile` / `Writer`. ~1200 LOC total for the edit surface.
- What we do NOT replicate from cfitsio: `ffirow`/`ffdrow`/`fficol`/`ffdcol`/`ffrsim`/`ffiimg`-style in-place structural mutation. Those live behind `EditFile` rebuild in v1. If a concrete v2 use case demands faster structural edits, we add them incrementally — the `EditFile` API doesn't break.
- Crash-safety story: cfitsio has none — a crash mid-edit may corrupt your file. v1 ships with journaled tail shifts + fsync barriers + atomic rename for `EditFile`. If the process dies during an in-place edit, the next `Open` detects the journal and replays or rolls back.

**5.7 Hybrid error strategy.** Sentinels for common checks (`ErrKeyNotFound`, `ErrNotImageHDU`, `ErrNotImageCompatible`, `ErrTypeMismatch`, `ErrReadOnly`, `ErrShapeMismatch`) via `errors.Is`. Typed errors with context (`*KeyNotFoundError{Key, HDU}`, `*TypeMismatchError{Requested, BITPIX, Lossy}`) via `errors.As` when callers need specifics. Standard Go practice.

**5.9 Context on blocking I/O only.** `OpenContext(ctx, path)`, `OpenForEditContext(ctx, path)`, `ReadPixelsContext[T](ctx, h)`, `EditFileContext(ctx, path, fn)` as canonical forms; `Open` / `OpenForEdit` / `ReadPixels` / `EditFile` are thin wrappers passing `context.Background()`. Precedent: `database/sql`. Pure in-memory ops (`Header.Get`, `Column()` metadata, `VarColumn.At`) don't take context.

**5.10 Case-insensitive keyword lookup, uppercase storage.** Headers store cards in uppercase per FITS §4.1.2.1, but `Header.Get` / `String` / `Int` / `Float` / `Bool` / `Set` / `Add` / `Delete` all normalize the requested keyword to uppercase before lookup. User code `h.Get("object")` and `h.Get("OBJECT")` both work. Internal storage is a `map[string]*Card` keyed on the uppercase form; `Cards()` returns entries in insertion order. Matches cfitsio and astropy conventions. Cost: one `strings.ToUpper` per lookup, negligible.

**5.11 Strict-fail on malformed files, no partial recovery in v1.** `Open` / `OpenReader` / `OpenForEdit` return a clear typed error and `nil *File` on any malformed input. No partial-recovery `*File` returned alongside an error — the pattern is unusual for Go, the forensic use case is rare, and keeping the happy path clean is worth more than the corner-case convenience.
- `ErrEmptyFile` — zero bytes.
- `ErrNotFITS` — first block does not begin `SIMPLE  =`.
- `*ErrUnterminatedHeader{HDU int, BytesScanned int64}` — no `END` card found within a generous block limit.
- `*ErrTruncatedData{HDU int, Expected int64, Got int64}` — data section ends before the declared `|BITPIX| × GCOUNT × (PCOUNT + Π NAXISi) / 8` size.
- `*ErrMissingRequiredKeyword{HDU int, Keyword string}` — any of `SIMPLE`/`BITPIX`/`NAXIS`/`NAXISn`/`XTENSION`/`PCOUNT`/`GCOUNT`/`TFIELDS` missing when required.
- Trailing garbage after the last padded data block is **silently ignored** (common from aborted writes; not worth failing over).
- Error messages include the HDU index and byte offset where the problem was detected. No stack traces, no hints about recovery.
- Forensic / recovery tooling is a **deferred v2 concern**. If users need to salvage data from corrupted files, we consider it when a concrete use case shows up. The promotion path from decision 5.6 (`(*File).Blocks()` iterator) covers this cleanly if it becomes needed.

## Done criteria for v1

Before tagging `v1.0.0` we must pass **every** item:

1. All 37 implementation steps from §"Implementation steps" complete.
2. Unit + integration tests green across every package.
3. `fitsverify` (from cfitsio) accepts every file the library writes across the full fixture set with zero errors and zero warnings.
4. `astropy.io.fits` parses every file the library writes with semantic match (headers, data, metadata all equal).
5. Byte-exact round-trip on all cfitsio fixtures (`testprog.std`, `iter_image.fit`, `iter_a.fit`, `iter_b.fit`, `iter_c.fit`, `vari.fits`, `testf77.std`) via `cmd/fitscopy`.
6. Byte-exact round-trip on at least one real ASI533MC Pro light frame from `~/Astrophotography/raw/` (including WCS/SIP keywords preserved unchanged).
7. Fuzz tests run ≥ 10 M iterations each (`FuzzParseHeader`, `FuzzTForm`, `FuzzBlockFraming`, `FuzzCardDecode`) without panic, data-race, or memory-safety finding.
8. `BenchmarkReadPixelsFloat32_1k²` within 2× of cfitsio throughput on the same hardware (measured on Linux, via the `fits-manager` venv's `fitsio` Python binding as the reference).
9. Public API reviewed against Go idioms (naming, receiver choices, method granularity, error wrapping) and against the guiding principle (no auto-anything).

All nine are required; we do not ship a partial v1.

## Package layout

```
fits/
├── fits.go              package fits — root types: File, HDU, Open, Create
├── hdu.go                            HDU interface, HDUType
├── image.go                          ImageHDU + generic pixel read/write
├── ascii_table.go                    ASCIITableHDU
├── bin_table.go                      BinaryTableHDU, Column, VarColumn
├── errors.go                         sentinel errors, typed errors
├── stdimage.go                       image.Image adapter, FloatGray, RegisterFormat
├── header/              package header — Header, Card, parser, encoder, keyword constants
│   ├── card.go
│   ├── parser.go
│   ├── encoder.go
│   ├── header.go
│   └── keywords.go
├── internal/block/      package block      — 2880-byte framed reader/writer
├── internal/bitpix/     package bitpix     — BITPIX enum, size/type mapping, scale/zero/blank
├── internal/bigendian/  package bigendian  — typed big-endian read/write, bulk byteswap
├── internal/tform/      package tform      — TFORMn and TDIMn parser
├── internal/checksum/   package checksum   — Appendix J checksum and ASCII encoding
├── cmd/
│   ├── fitsdump/                    human-readable HDU/header dump
│   └── fitscopy/                    round-trip copy tool for regression testing
├── testdata/                        cfitsio fixtures + synthetic test files
└── reference/                       spec PDFs + cfitsio C source (not shipped)
```

**Layout notes.**
- Root package owns user vocabulary (`fits.File`, `fits.HDU`, `fits.ImageHDU`, ...). Matches Go stdlib conventions (`image.Image` in `image`, not `image/core`).
- `header` is its own package because card parsing is self-contained and heavily tested in isolation. May alias `type Header = header.Header` into root for ergonomics.
- Everything else is `internal/` so users can't freeze implementation details into the public API before we know their final shape.
- `cmd/fitscopy` is the byte-for-byte regression tool: `fitscopy in.fits out.fits && cmp in.fits out.fits`.
- No `fits/compress` subpackage in v1.

## Core API sketch (reflects approved decisions)

```go
package fits

// --- file ---

type File struct { /* unexported */ }

type Mode int
const (
    ModeRead   Mode = iota // Open, OpenReader, OpenReadAll
    ModeEdit               // OpenForEdit — in-place mutations + appends
    ModeCreate             // Create — new file
)

func Open(name string) (*File, error)                       // read-only os.File
func OpenContext(ctx context.Context, name string) (*File, error)
func OpenReader(r io.ReadSeeker) (*File, error)             // read-only, any seeker
func OpenReadAll(r io.Reader) (*File, error)                // slurp to memory helper
func OpenForEdit(name string) (*File, error)                // read-write in-place
func OpenForEditContext(ctx context.Context, name string) (*File, error)
func Create(name string) (*File, error)                     // new file, write mode
func (f *File) Mode() Mode
func (f *File) Close() error
func (f *File) NumHDU() int
func (f *File) HDU(i int) (HDU, error)                      // 0-indexed; 0 = primary
func (f *File) HDUByName(name string) (HDU, error)          // EXTNAME lookup
func (f *File) Primary() (*ImageHDU, error)

// Flush persists any pending in-memory header mutations. ModeEdit / ModeCreate
// only. Returns ErrReadOnly for ModeRead files. If a header grew past its
// current block allocation and a tail shift is required, Flush uses the
// journal protocol for crash-safety.
func (f *File) Flush() error

// Append writes a new HDU at EOF. ModeEdit / ModeCreate only. Generic
// convenience wrappers for each HDU kind:
func AppendImage[T Numeric](f *File, hdr *header.Header, shape []int64, data []T) (*ImageHDU, error)
func AppendBinaryTable(f *File, hdr *header.Header, cols []Column, rows RowSource) (*BinaryTableHDU, error)
func AppendASCIITable(f *File, hdr *header.Header, cols []Column, rows RowSource) (*ASCIITableHDU, error)

// EditFile is the streaming-rebuild helper for structural edits (image
// resize, row/column insert/delete, HDU insert-in-middle, HDU delete, VLA
// heap reorganization). Opens src read-only, creates a temp file in the same
// directory, passes (in, out) to the callback, fsyncs, and atomically renames
// over src on success. On any error the temp file is removed and src is
// untouched. Crash-safe by construction.
func EditFile(path string, fn func(in *File, out *Writer) error) error
func EditFileContext(ctx context.Context, path string, fn func(in *File, out *Writer) error) error

// Writer is the streaming output side of EditFile. Users emit HDUs to it
// using CopyHDU (pass-through), CopyHDUWithHeader (new header, same data
// bytes — no parse round-trip), SkipHDU (omit from output), or the
// Append*/WriteImage-family functions to add new HDUs.
type Writer struct { /* unexported */ }
func (w *Writer) CopyHDU(h HDU) error                               // bit-for-bit passthrough
func (w *Writer) CopyHDUWithHeader(h HDU, hdr *header.Header) error // new header, same data
func (w *Writer) SkipHDU(h HDU) error                               // omit (cannot skip primary)
// Append*/WriteImage-family functions also accept a *Writer for emitting new HDUs.

// --- HDU interface ---

type HDUType int
const (
    TypeImage HDUType = iota
    TypeASCIITable
    TypeBinaryTable
)

type HDU interface {
    Type() HDUType
    Header() *header.Header
    // narrow via type assertion to *ImageHDU / *BinaryTableHDU / *ASCIITableHDU
}

// --- header (aliased from header package) ---

type Header = header.Header
type Card   = header.Card

// header package methods:
func (h *Header) Get(name string) (Card, bool)
func (h *Header) String(name string) (string, error)
func (h *Header) Int(name string) (int64, error)
func (h *Header) Float(name string) (float64, error)
func (h *Header) Bool(name string) (bool, error)
func (h *Header) Set(name string, value any, comment string) error   // update or create
func (h *Header) Add(name string, value any, comment string) error   // append; allows duplicates (HISTORY, COMMENT)
func (h *Header) Delete(name string) error
func (h *Header) Clone() *Header                                     // deep copy, safe to mutate
func (h *Header) Cards() []Card                                      // in order
func (h *Header) History() []string
func (h *Header) Comments() []string

// Mutation persistence semantics:
//   - Headers from ModeRead files: Set/Add/Delete succeed in memory, but
//     (*File).Flush() returns ErrReadOnly. Mutations do not reach disk.
//   - Headers from ModeEdit files: Set/Add/Delete are tracked; (*File).Flush()
//     writes them back in place, using the journaled tail-shift protocol if
//     the header grew past its current block allocation.
//   - Headers created with header.New() and passed to Append*/Writer methods:
//     serialized as part of the new HDU.

type Card struct {
    Key     string     // "NAXIS1", "HIERARCH ESO DET NAME", ...
    Value   any        // string | int64 | float64 | bool | complex128 | nil
    Comment string
    Raw     [80]byte   // original bytes, for round-tripping
}

// --- image HDU ---

type Numeric interface {
    ~uint8  | ~int8  | ~int16 | ~uint16 | ~int32 | ~uint32 |
    ~int64  | ~uint64 | ~float32 | ~float64
}

type ImageHDU struct { /* unexported */ }

func (h *ImageHDU) Type() HDUType
func (h *ImageHDU) Header() *header.Header
func (h *ImageHDU) BITPIX() int                              // 8, 16, 32, 64, -32, -64
func (h *ImageHDU) NAXIS()  int
func (h *ImageHDU) Shape()  []int64                          // [NAXIS1, NAXIS2, ...]
func (h *ImageHDU) BSCALE() float64
func (h *ImageHDU) BZERO()  float64

// Read all pixels, generic. Applies BSCALE/BZERO. Returns ErrTypeMismatch
// if T cannot losslessly hold the file's data type.
func ReadPixels[T Numeric](h *ImageHDU) ([]T, error)
func ReadPixelsContext[T Numeric](ctx context.Context, h *ImageHDU) ([]T, error)

// Like ReadPixels, but also returns a validity mask (true = valid, false = NaN/BLANK).
func ReadPixelsMasked[T Numeric](h *ImageHDU) (data []T, mask []bool, err error)

// N-dimensional hyperslab: [lower..upper) with stride.
func ReadSubset[T Numeric](h *ImageHDU, lower, upper, stride []int64) ([]T, error)

// Pre-flight check: can the file's data be represented as T without loss?
func CanConvert[T Numeric](h *ImageHDU) bool

// Write side — for use with Create / ModeEdit Files or inside an EditFile callback Writer:
func WriteImage[T Numeric](dst WriteTarget, h *header.Header, shape []int64, data []T) (*ImageHDU, error)

// WriteTarget abstracts "a *File in write/edit mode" and "a *Writer from EditFile"
// so Append/WriteImage-family functions work in both contexts.
type WriteTarget interface { /* unexported methods */ }

// OverwritePixels writes new pixel values over the existing data allocation
// of an ImageHDU. Requires ModeEdit. Returns ErrShapeMismatch if len(data)
// does not match the HDU's pixel count, or ErrTypeMismatch if T cannot be
// losslessly encoded into the HDU's BITPIX. Applies inverse BSCALE/BZERO.
// Does NOT support resize — for resize, use EditFile.
func OverwritePixels[T Numeric](h *ImageHDU, data []T) error

// Pin one or more axes of a higher-dim cube to produce a lower-dim view.
// Used to extract a 2D slice before calling AsImage().
func (h *ImageHDU) Slice(axis int, index int64) (*ImageHDU, error)

// --- image.Image adapter ---

// AsImage returns a native image.Image for 2D integer-BITPIX HDUs. For float
// BITPIX it returns a *FloatGray (which implements image.Image but holds raw
// float32 values — faithful, not scaled). Returns ErrNotImageCompatible for
// multi-dim, complex, or anything else that cannot be represented without loss.
func (h *ImageHDU) AsImage() (image.Image, error)

// FloatGray is image.Image over raw float32 pixels. Does NOT rescale.
// Rendering through png.Encode without a caller-provided transform produces
// garbage — that is intentional. The purpose of this type is to let float FITS
// data participate in the image.Image ecosystem with lossless fidelity.
type FloatGray struct {
    Rect   image.Rectangle
    Stride int
    Pix    []float32
}

func (f *FloatGray) ColorModel() color.Model
func (f *FloatGray) Bounds()     image.Rectangle
func (f *FloatGray) At(x, y int) color.Color   // returns the raw float32, boxed

// image.RegisterFormat registration happens in init():
func Decode(r io.Reader)       (image.Image,  error)
func DecodeConfig(r io.Reader) (image.Config, error)

// --- table HDUs ---

type ColumnType int
const (
    ColByte ColumnType = iota
    ColInt16
    ColInt32
    ColInt64
    ColFloat32
    ColFloat64
    ColString
    ColLogical
    ColBit
    ColComplex64
    ColComplex128
    ColVarArray
)

type Column struct {
    Index   int        // 1-based per FITS convention
    Name    string     // TTYPEn
    Unit    string     // TUNITn
    TForm   string     // raw TFORMn
    Repeat  int64
    Type    ColumnType
    Scale   float64    // TSCALn, default 1
    Zero    float64    // TZEROn, default 0
    Null    any        // TNULLn if set
    Display string     // TDISPn
    Dim     []int64    // TDIMn if set
}

type BinaryTableHDU struct { /* unexported */ }
type ASCIITableHDU  struct { /* unexported */ }

func (t *BinaryTableHDU) NumRows() int64
func (t *BinaryTableHDU) Columns() []Column
func (t *BinaryTableHDU) ColumnByName(name string) (Column, bool)
func (t *BinaryTableHDU) ColumnIndex(name string) int
func (t *BinaryTableHDU) RowBytes(row int64) ([]byte, error)   // escape hatch

// Scalar column (Repeat == 1).
func ReadColumn[T Numeric | ~string | ~bool](t *BinaryTableHDU, col int) ([]T, error)

// Fixed-width vector column (Repeat > 1) returns rowcount × repeat flattened.
// Reshape with Repeat and TDIM info.
func ReadVectorColumn[T Numeric](t *BinaryTableHDU, col int) ([]T, error)

// Variable-length array column (P/Q TFORM). See VarColumn.
func ReadVarColumn[T Numeric](t *BinaryTableHDU, col int) (*VarColumn[T], error)

// VarColumn stores all rows' variable-length data in one contiguous slab.
// At(row) returns a zero-copy sub-slice into the backing Values buffer.
type VarColumn[T Numeric] struct {
    // unexported: values []T, offsets []int64
}
func (v *VarColumn[T]) Len() int                                 // number of rows
func (v *VarColumn[T]) At(row int) []T                           // zero-copy view
func (v *VarColumn[T]) Rows() iter.Seq2[int, []T]                // range iterator
func (v *VarColumn[T]) Raw() (values []T, offsets []int64)       // escape hatch
```

### User-code examples

```go
// Read primary image
f, _ := fits.Open("m51.fits")
defer f.Close()
img, _ := f.Primary()
obj, _ := img.Header().String("OBJECT")
pix, _ := fits.ReadPixels[float32](img)
fmt.Println(obj, len(pix), img.Shape())

// Iterate HDUs
for i := 0; i < f.NumHDU(); i++ {
    h, _ := f.HDU(i)
    switch v := h.(type) {
    case *fits.ImageHDU:       fmt.Println(i, "image", v.Shape())
    case *fits.BinaryTableHDU: fmt.Println(i, "bintable", v.NumRows())
    case *fits.ASCIITableHDU:  fmt.Println(i, "asciitable")
    }
}

// Read a table column
tbl := h.(*fits.BinaryTableHDU)
flux, _ := fits.ReadColumn[float64](tbl, tbl.ColumnIndex("FLUX"))

// Variable-length array (e.g. per-event photon PHA lists)
col, _ := fits.ReadVarColumn[float32](tbl, tbl.ColumnIndex("PHA"))
for i, row := range col.Rows() {
    process(i, row)                     // row is a zero-copy []float32
}

// image.Image — only lossless for 2D integer BITPIX
raw, _ := os.ReadFile("star_u8.fits")
m, name, _ := image.Decode(bytes.NewReader(raw))   // name == "fits"
png.Encode(out, m)

// Write a new image
f2, _ := fits.Create("out.fits")
h := header.New()
h.Set("OBJECT", "synthetic", "")
fits.WriteImage(f2, h, []int64{1024, 1024}, mkData())
f2.Close()

// In-place header edit on an existing file (ModeEdit)
fe, _ := fits.OpenForEdit("obs.fits")
p, _ := fe.Primary()
p.Header().Add("HISTORY", "Calibrated 2026-04-05", "")
p.Header().Set("DATE-OBS", "2026-04-05T12:34:56", "")
fe.Flush()   // writes modified header blocks in place (journaled if tail shift required)
fe.Close()

// Overwrite pixel data in place (same shape)
fe, _ = fits.OpenForEdit("mask.fits")
img, _ := fe.Primary()
fits.OverwritePixels[uint8](img, newMask)
fe.Close()

// Append a new HDU at EOF in edit mode
fe, _ = fits.OpenForEdit("obs.fits")
h2 := header.New()
h2.Set("EXTNAME", "CALIBRATED", "")
fits.AppendImage(fe, h2, []int64{1024, 1024}, calibrated)
fe.Close()

// Structural edit: drop an HDU by name via EditFile streaming rebuild
err := fits.EditFile("obs.fits", func(in *fits.File, out *fits.Writer) error {
    for i := 0; i < in.NumHDU(); i++ {
        h, _ := in.HDU(i)
        if name, _ := h.Header().String("EXTNAME"); name == "JUNK" {
            out.SkipHDU(h)    // omit from output
            continue
        }
        if err := out.CopyHDU(h); err != nil {
            return err
        }
    }
    return nil
})
```

## Implementation steps

Sequential. No phases.

1. **Module init.** `go mod init github.com/dmarkham/fits`. `.gitignore`, `LICENSE`, `README.md` stub. Minimum Go 1.23 (generics + `iter.Seq2`). Create the directory skeleton from the layout above.
2. **`internal/block`: 2880-byte block reader/writer.** Over `io.ReadSeeker` / `io.WriteSeeker`. `ReadBlock`, `WriteBlock`, `BlockCount`. Tiny 1-block read cache on top of the OS page cache. Mirrors cfitsio `buffers.c` (`ffmbyt` / `ffgbyt` / `ffpbyt`) without the 40-buffer LRU. Spec: §3.
3. **`internal/bigendian`: typed big-endian helpers.** `Int16/32/64`, `Float32/64`, `Put*` counterparts, `BulkSwap16/32/64` for converting `[]byte` slabs in place to native-ordered typed slices. Built on `encoding/binary.BigEndian`. Spec: §5.2–5.3.
4. **`header/card.go`: Card type + encode/decode.** One struct, `any` value. Decoder handles strings (escaped `''`), logicals, ints, floats with `D`-exponent, complex. Reference: `fitscore.c` `ffpsvc` / `ffdtyp`. Spec: §4.2.
5. **`header/parser.go`: card stream parser.** Input: header byte slice (multiple of 2880). Output: `[]Card`. Handles fixed / free format, `CONTINUE`, `HIERARCH`, `COMMENT` / `HISTORY` / blank, locates `END`. Reference: `getkey.c` `ffgthd` / `ffgknm`. Spec: §4.1, Appendix A.
6. **`header/encoder.go`: card serializer.** Symmetric to parser. Fixed-format for mandatory keys. CONTINUE splitting for long strings. Reference: `putkey.c` `ffmkky` / `ffpkys`.
7. **`header/keywords.go` + `header/header.go`: Header container.** Ordered cards + name→index map (last-wins for duplicates, preserve insertion order on serialization). Typed getters / setters. Constants for every mandatory / reserved keyword. Spec: §4.4, Appendix C.
8. **`internal/bitpix`: BITPIX enum + Go type mapping.** `type BITPIX int8`, `Size()`, `Signed()`, mapping helpers. Spec: §4.4.1.1, Table 8.
9. **`internal/checksum`: Appendix J checksum.** Port `checksum.c` `ffcsum` (16-bit-folded 32-bit 1's complement over big-endian uint16), `ffesum` (base-64-ish ASCII encoding with the `[\]^_{|}~` exclusions), `ffdsum`. Spec: Appendix J.
10. **`fits.File` + `Open` / `OpenReader` / `OpenReadAll` / `Create`.** Holds the seeker, mode, and a slice of `hduOffset` records `{headerStart, dataStart, dataEnd, kind}`.
11. **HDU-scan pass (minimal eager, per decision 5.2).** Walk file: read header blocks until `END`, run a *restricted* keyword parse that extracts only `XTENSION`/`SIMPLE`, `BITPIX`, `NAXIS`, `NAXISn`, `PCOUNT`, `GCOUNT`, `EXTNAME`, `EXTVER`. Retain the raw header bytes as `[]byte` in the `hduRecord`. Compute data size (Eq. 2: `|BITPIX| × GCOUNT × (PCOUNT + NAXIS1 × ... × NAXISm) / 8`), round up to block boundary, seek past data, record offsets, repeat to EOF. Reference: `cfileio.c` `ffrhdu` / `ffphdu`. Spec: §3.1, §3.3.2, §4.4.1.
12. **`HDU` interface + factory + lazy full-header parse.** Inspect the extracted `XTENSION` / `SIMPLE` to pick `*ImageHDU` / `*BinaryTableHDU` / `*ASCIITableHDU`. Handle SIF (primary only) and MEF (primary + extensions). `(*HDU).Header()` uses `sync.Once` to parse the retained `rawHeader` bytes into a full `header.Header` on first access and cache the result. TFORM parsing for binary tables hangs off `BinaryTableHDU.Columns()` / `ReadColumn` via the same pattern.
13. **`ImageHDU.ReadPixels[T]` — read path.** Seek to `dataStart`, read bytes, big-endian bulk-swap, apply BSCALE / BZERO / BLANK / NaN, return `[]T`. Single generic implementation replaces cfitsio's `getcolb.c`..`getcold.c`. Spec: §5, §7.1, Table 11.
14. **`ImageHDU.ReadSubset[T]`.** Hyperslab reader. Computes byte ranges per-axis, issues seeks when gaps exceed a threshold (~64 KB). Reference: `ffgsv` family. Spec: §3.3.2.
15. **Image write path.** `WriteImage[T]` emits `SIMPLE` / `XTENSION=IMAGE`, `BITPIX`, `NAXIS` / `NAXISn`, `END`, block-pads, writes big-endian data with inverse BSCALE / BZERO, zero-pads final block. Reference: `putcol*.c`, `ffiimg`.
16. **`internal/tform`: TFORM and TDIM parser.** ASCII table formats (`Aw`, `Iw`, `Fw.d`, `Ew.d`, `Dw.d`; Table 15) and binary table formats (`rL rX rB rI rJ rK rA rE rD rC rM rP(T,max) rQ(T,max)`; Table 18). Reference: `fitscore.c` `ffbnfm` / `ffbnfmll` / `ffasfm`.
17. **`BinaryTableHDU` column metadata.** Parse `TFIELDS` + all per-column keywords; compute row stride + per-column byte offsets. Spec: Table 17, §7.3.2.
18. **`BinaryTableHDU.ReadColumn[T]` — scalar.** One column across all rows, row-strided read, bulk byteswap + TSCAL / TZERO application. Reference: `getcol*.c`.
19. **`BinaryTableHDU.ReadVectorColumn[T]`** for repeat > 1.
20. **Variable-length arrays (P / Q).** Parse row descriptors (nelem, heap offset), seek to `dataStart + THEAP + offset`, read the whole heap region in one pass, construct `VarColumn[T]` with the flat slab + offsets table. Reference: `ffgdes` / `ffgdesll`. Spec: §7.3.5.
21. **`ASCIITableHDU`.** Row = NAXIS1 chars, split by TBCOLn windows, parse per TFORMn with `strconv`. Spec: §7.2.
22. **Table write paths.** Binary: write header, allocate heap for P / Q writes, write rows, write heap, pad. ASCII: Fortran-rule formatting, right-justified. Reference: `putcol*.c`, `ffphbn` / `ffphtb`.
23. **Checksum integration.** On `Close` in write / edit mode, recompute `DATASUM` per HDU that changed, then `CHECKSUM` over header + data using the placeholder-then-rewrite trick from §J.1. `HDU.VerifyChecksum()` on read.
24. **`image.Image` adapter.** `AsImage()` returning `*image.Gray` for BITPIX=8, `*image.Gray16` for BITPIX=16 (with TZERO=32768 unsigned handling), `*FloatGray` for -32 / -64. Errors on multi-dim, complex, NAXIS < 2. Implement `Decode` / `DecodeConfig`. Register in `init()` with magic bytes `"SIMPLE  ="`.
25. **`OpenForEdit` + dirty-header tracking (decision 5.8).** Open file read-write, load HDU index + minimal eager metadata as for `Open`. Each `hduRecord` gains a `dirty bool` flag flipped when the user mutates the associated `Header`. `Flush` walks dirty records, re-serializes each header, compares byte length against current block allocation.
26. **`File.Flush` — in-place header rewrite within block allocation.** For each dirty HDU whose re-serialized header still fits in its current allocation: seek to `headerStart`, write the new header bytes (padded with spaces to the block boundary), fsync. Fast path. Reference: `putkey.c` `ffpkys` / `ffpky` / `ffmkey` — but we only do the rewrite, not the cfitsio edit semantics.
27. **`File.Flush` — journaled tail-shift for header growth across block boundary.** If any dirty HDU's new header is larger than its current allocation: write a journal file `<path>.fits.journal` describing `{hdu_index, old_block_count, new_block_count, shift_offset}`, fsync journal, shift the file tail forward in back-to-front chunks, write the new header, fsync data, fsync directory, delete journal. On next `Open` / `OpenForEdit`, detect a stale journal and replay or roll back. New implementation — cfitsio has no crash-safety story here.
28. **`OverwritePixels[T]`.** Edit-mode only. Validate shape and BITPIX compatibility, seek to `dataStart`, apply inverse BSCALE/BZERO, write bytes in big-endian, fsync. Returns `ErrShapeMismatch` / `ErrTypeMismatch` on invalid calls. Rejects resize requests (those go through `EditFile`).
29. **`(*File).AppendImage` / `AppendBinaryTable` / `AppendASCIITable` in edit mode.** Seek to EOF (past the last HDU's padded data), emit new HDU using the same serialization path as `Create` mode, update HDU index. Fast — no shift needed because append is at EOF.
30. **`EditFile` + `Writer` streaming rebuild.** Open source read-only, compute temp path `<path>.fits.tmp-<pid>` in the same directory (same filesystem → atomic rename). Create temp file in `ModeCreate`. Wrap temp file in `*Writer`. Invoke callback with `(in, writer)`. On clean return: fsync temp, `os.Rename(temp, path)`, fsync directory. On error or panic: remove temp, leave source untouched. `Writer.CopyHDU` passes header bytes and data bytes through verbatim (no parse round-trip). `Writer.CopyHDUWithHeader` re-emits the header but passes data through as raw bytes. `Writer.SkipHDU` marks the HDU as omitted; rejects skipping the primary.
31. **`cmd/fitsdump`.** Dumps headers and per-HDU summary. Dogfooding + debug output for tests.
32. **`cmd/fitscopy`.** Primary regression tool: open → iterate → write new → byte-compare.
33. **Test fixtures.** Copy `iter_image.fit`, `iter_a.fit`, `iter_b.fit`, `iter_c.fit`, `vari.fits`, `testprog.std`, `testf77.std` from `reference/cfitsio/` into `testdata/`. Add a Python script using `astropy.io.fits` to generate a small synthetic set with known values (one-time, committed).
34. **Unit tests per package.** `header/parser_test.go` covering every card edge case. `internal/tform_test.go` covering every TFORM. `internal/checksum_test.go` against §J.3 worked examples. Edit-surface tests: in-place header rewrite within block, journaled tail-shift, `OverwritePixels` validation, `EditFile` round-trip, `Writer.SkipHDU` + `CopyHDU` mixing.
35. **Integration round-trip tests.** `fitscopy` on every fixture, `bytes.Equal`. `vari.fits` is the critical variable-length-array case. `testprog.std` covers nearly every cfitsio feature. Edit round-trip: open for edit, mutate a header, flush, reopen, verify persistence + `VerifyChecksum()`.
36. **Crash-safety tests for the edit surface.** Simulate process death at each fsync point during a journaled tail-shift; verify recovery on next open produces either the pre-edit or post-edit state, never a corrupted intermediate. Simulate rename failure during `EditFile`; verify source is untouched.
37. **Fuzz + benchmarks.** `FuzzParseHeader`, `FuzzTForm`, `FuzzBlockFraming` — all seeded with real fixtures. Benchmarks: `BenchmarkReadPixelsFloat32_1k²`, `BenchmarkReadColumnInt32_1Mrow`, `BenchmarkWriteImage`, `BenchmarkEditHeaderInPlace`, `BenchmarkEditFileRebuild`. Target: within 2× of cfitsio throughput for reads; in-place header edit should complete in milliseconds regardless of file size.

## Testing strategy

1. **Unit tests per internal package.** Hand-crafted card / TFORM / checksum inputs for every spec edge case.
2. **Golden cfitsio fixtures in `testdata/`.** `iter_image.fit`, `iter_a.fit`, `iter_b.fit`, `iter_c.fit`, `vari.fits`, `testprog.std`, `testf77.std`. `testprog.std` is cfitsio's own correctness oracle — 69 KB exercising nearly every feature.
3. **Byte-for-byte round-trip tests.** Open → iterate all HDUs → rewrite via `fits.Create` → `bytes.Equal`. `vari.fits` is the must-pass for variable-length arrays + heap.
4. **Semantic cross-check with `astropy.io.fits`.** One-time Python dump of each fixture as JSON / NumPy, committed as `.golden`. Tests assert Go's parsed view matches.
5. **Checksum verification.** Every fixture that sets `CHECKSUM` / `DATASUM` must verify.
6. **Optional CI: live cfitsio cross-check.** Build cfitsio in CI, run `fitsverify` on our output files. Catches "valid enough but cfitsio would reject" cases.
7. **Fuzz tests.** `FuzzParseHeader` (seeded from real headers, must not panic, must round-trip any successful parse). `FuzzTForm`. `FuzzBlockFraming`.
8. **Benchmarks.** Single-HDU pixel reads, column reads, image writes. 2× cfitsio throughput as a target.
9. **`testdata/README.md`** documenting the provenance and purpose of every fixture.

## Reference files (start here when implementing each step)

- `reference/fits_standard_v4.0.pdf` — §3 (file organization), §4 (header cards), §7 (extensions), §J (checksum), Appendix A (keyword syntax), Appendix C (keyword dictionary)
- `reference/cfitsio/fitsio.h` — public API surface the port must cover (names we are replacing, not copying)
- `reference/cfitsio/cfileio.c` — file open, HDU discovery, block I/O (mirrors `fits.File` + HDU scan)
- `reference/cfitsio/fitscore.c` — keyword parsing, header validation, TFORM parsing
- `reference/cfitsio/buffers.c` — 2880-byte block I/O
- `reference/cfitsio/checksum.c` — Appendix J implementation
- `reference/cfitsio/putkey.c` / `getkey.c` / `modkey.c` — card encode / decode / edit
- `reference/cfitsio/getcol*.c` / `putcol*.c` — per-type pixel / column access (collapsed into one generic function in our port)
- `reference/cfitsio/testprog.c` + `testprog.std` — end-to-end correctness oracle
