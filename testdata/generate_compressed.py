#!/usr/bin/env python3
"""Generate tile-compressed FITS test fixtures using astropy.

Usage (the repo venv; versions pinned in testdata/COMP_FZ_PROVENANCE.md):

    venv/bin/python testdata/generate_compressed.py
    make -C testdata/cfref goldens     # cfitsio truth for every fixture

This script produces a matched set of compressed and raw-uncompressed
FITS files in testdata/. Every compressed file has a corresponding
comp_raw_<name>.fits file containing the same deterministic data through
a plain ImageHDU. For lossless fixtures the raw file equals the decoded
pixels; for quantized (lossy) float fixtures it holds the ORIGINAL floats,
so the Go tests compare those fixtures against the cfitsio decode of the
compressed file (testdata/cfref/golden/comp_<name>.fits.json), not against
the raw twin.

The fixture set covers:

  - Every common tile-compression algorithm (RICE_1, GZIP_1, GZIP_2,
    HCOMPRESS_1, PLIO_1, NOCOMPRESS)
  - Both integer and float pixel types
  - Single-row tiles (the astropy default) and non-trivial 2D tiles
  - Partial edge tiles (image size not divisible by tile size)
  - NO_DITHER, SUBTRACTIVE_DITHER_1 and SUBTRACTIVE_DITHER_2 quantization
    for floats, each selected explicitly with quantize_method and a fixed
    dither_seed (astropy's default dither_seed=0 seeds from the clock, and
    dither_seed alone does not select the method)

All data is generated with fixed seeds and every dither seed is explicit,
so re-running reproduces every byte (checked: run twice, cmp).
"""

import os
import numpy as np
import astropy.io.fits as fits

import functools
import gzip
import astropy.io.fits.hdu.compressed._codecs as _codecs

# astropy's GZIP_1/GZIP_2 codecs call gzip.compress(), which stamps the
# current time into every tile's gzip header; pin it so the bytes are
# reproducible (the timestamp is not part of the pixel data).
_codecs.gzip_compress = functools.partial(gzip.compress, mtime=0)

# astropy.io.fits quantize_method values (= cfitsio's).
NO_DITHER = -1
SUBTRACTIVE_DITHER_1 = 1
SUBTRACTIVE_DITHER_2 = 2


TESTDATA = os.path.dirname(os.path.abspath(__file__))
RNG = np.random.default_rng(42)


def write_pair(name, data, cls_kwargs, ext_name=None, raw_of=None, raw_decoded=False):
    """Write a compressed fixture and its raw-uncompressed counterpart.

    raw_of names an earlier fixture whose comp_raw_ file already holds the
    same plain image: no copy is written (the twins would be
    byte-identical; compressed_test.go maps the fixture to that file).
    raw_decoded writes astropy's own decode of the compressed file as the
    raw twin instead of the input (for quantized NOCOMPRESS floats, which
    cfitsio cannot read: astropy is their truth)."""
    comp_path = os.path.join(TESTDATA, f"comp_{name}.fits")
    raw_path = os.path.join(TESTDATA, f"comp_raw_{raw_of or name}.fits")

    comp_hdu = fits.CompImageHDU(data=data, name=ext_name or "IMG", **cls_kwargs)
    fits.HDUList([fits.PrimaryHDU(), comp_hdu]).writeto(comp_path, overwrite=True, checksum=False)

    if raw_of is None:
        raw = data
        if raw_decoded:
            with fits.open(comp_path) as hl:
                raw = hl[1].data.copy()
        raw_hdu = fits.ImageHDU(data=raw, name=ext_name or "IMG")
        fits.HDUList([fits.PrimaryHDU(), raw_hdu]).writeto(raw_path, overwrite=True)
    else:
        with fits.open(raw_path) as hl:
            assert np.array_equal(hl[1].data, data), f"{name}: raw_of={raw_of} holds other pixels"

    print(f"  {name:30s} comp={os.path.getsize(comp_path):8d} raw={os.path.basename(raw_path)}")


def main():
    import astropy, numpy
    print(f"Generating compressed FITS fixtures in testdata/ (astropy {astropy.__version__}, numpy {numpy.__version__}):")

    # ---------- Integer image, every algorithm ----------
    # 40 x 37 images (1480 pixels): the fixtures prove algorithms and tile
    # geometry, not size (they were 100 x 100 before the round-2 size
    # budget; every property below survives the smaller image).
    i16_small = (np.arange(37 * 40, dtype=np.int32) * 23 - 9000).astype(np.int16).reshape((37, 40))
    i32_small = (np.arange(37 * 40, dtype=np.int32) * 40009 - 30000000).reshape((37, 40))

    write_pair("rice_i16", i16_small, {"compression_type": "RICE_1"})
    write_pair("rice_i32", i32_small, {"compression_type": "RICE_1"})
    write_pair("gzip1_i16", i16_small, {"compression_type": "GZIP_1"}, raw_of="rice_i16")
    write_pair("gzip1_i32", i32_small, {"compression_type": "GZIP_1"}, raw_of="rice_i32")
    write_pair("gzip2_i16", i16_small, {"compression_type": "GZIP_2"}, raw_of="rice_i16")
    write_pair("gzip2_i32", i32_small, {"compression_type": "GZIP_2"}, raw_of="rice_i32")
    write_pair("hcompress_i16", i16_small, {"compression_type": "HCOMPRESS_1"}, raw_of="rice_i16")
    write_pair("nocompress_i16", i16_small, {"compression_type": "NOCOMPRESS"}, raw_of="rice_i16")

    # ---------- Float image with quantization ----------
    # Smooth gradient + gaussian noise, reproducible via fixed RNG.
    f32_smooth = (
        np.linspace(0, 1000, 37 * 40).reshape((37, 40))
        + RNG.normal(0, 1, (37, 40))
    ).astype(np.float32)

    # RICE_1 with NO_DITHER.
    write_pair(
        "rice_f32_nodither",
        f32_smooth,
        {"compression_type": "RICE_1", "quantize_level": 16,
         "quantize_method": NO_DITHER},
    )
    # RICE_1 with SUBTRACTIVE_DITHER_1, ZDITHER0 = 1 (exercises the
    # Park-Miller dither sequence; first tile uses random index 0).
    write_pair(
        "rice_f32_dither1",
        f32_smooth,
        {"compression_type": "RICE_1", "quantize_level": 16,
         "quantize_method": SUBTRACTIVE_DITHER_1, "dither_seed": 1},
        raw_of="rice_f32_nodither",
    )
    # RICE_1 with SUBTRACTIVE_DITHER_2 (zero pixels preserved exactly).
    f32_zeros = f32_smooth.copy()
    f32_zeros[::7, ::5] = 0.0
    write_pair(
        "rice_f32_dither2",
        f32_zeros,
        {"compression_type": "RICE_1", "quantize_level": 16,
         "quantize_method": SUBTRACTIVE_DITHER_2, "dither_seed": 77},
    )
    # GZIP_1 float (quantized, NO_DITHER).
    write_pair(
        "gzip1_f32",
        f32_smooth,
        {"compression_type": "GZIP_1", "quantize_level": 16,
         "quantize_method": NO_DITHER},
        raw_of="rice_f32_nodither",
    )
    # GZIP_2 float (quantized, SUBTRACTIVE_DITHER_1).
    write_pair(
        "gzip2_f32",
        f32_smooth,
        {"compression_type": "GZIP_2", "quantize_level": 16,
         "quantize_method": SUBTRACTIVE_DITHER_1, "dither_seed": 5},
        raw_of="rice_f32_nodither",
    )
    # HCOMPRESS float (quantized, NO_DITHER).
    write_pair(
        "hcompress_f32",
        f32_smooth,
        {"compression_type": "HCOMPRESS_1", "quantize_level": 16,
         "quantize_method": NO_DITHER},
        raw_of="rice_f32_nodither",
    )

    # ---------- PLIO (sparse integer mask) ----------
    # PLIO needs non-negative integer mask-like data. Make a small mask.
    mask_i16 = np.zeros((50, 50), dtype=np.int16)
    mask_i16[10:20, 15:35] = 1
    mask_i16[30:40, 5:45] = 2
    write_pair("plio_mask", mask_i16, {"compression_type": "PLIO_1"})

    # ---------- Tile geometry edge cases ----------
    # Non-trivial 2D tiles (16x16) over a 40x37 image → edge tiles are
    # partial (8 pixels wide in the last column, 5 pixels tall in the last
    # row).
    write_pair(
        "rice_multitile",
        i16_small,
        {"compression_type": "RICE_1", "tile_shape": (16, 16)},
        raw_of="rice_i16",
    )
    # Whole-image-as-one-tile.
    write_pair(
        "rice_onetile",
        i16_small,
        {"compression_type": "RICE_1", "tile_shape": (37, 40)},
        raw_of="rice_i16",
    )
    # Odd-sized image (23x19) with default tiles — exercises partial edge
    # behavior at a non-power-of-two boundary.
    i16_odd = (np.arange(23 * 19, dtype=np.int32) * 137 - 30000).astype(np.int16).reshape((19, 23))
    write_pair(
        "rice_odd",
        i16_odd,
        {"compression_type": "RICE_1"},
    )

    # ---------- NOCOMPRESS at the other ZBITPIX ----------
    # cfitsio 4.7.0 cannot read NOCOMPRESS, so astropy is the truth. For
    # integers the raw twin is the input; astropy writes float NOCOMPRESS
    # quantized (ZSCALE/ZZERO columns, SUBTRACTIVE_DITHER_1 with the seed
    # given here), so the raw twin of a float fixture is astropy's own
    # decode of the written file.
    rng = np.random.default_rng(7)
    nc = rng.normal(100, 5, (12, 16))
    write_pair("nocompress_u8", (nc + 28).astype(np.uint8), {"compression_type": "NOCOMPRESS"})
    write_pair("nocompress_i32", (nc * 1e5).astype(np.int32), {"compression_type": "NOCOMPRESS"})
    nc64 = (nc * 1e5).astype(np.int64) * (2**28) + 3   # past 2^53: not exact in float64
    write_pair("nocompress_i64", nc64, {"compression_type": "NOCOMPRESS"})
    for bp, dt in (("f32", np.float32), ("f64", np.float64)):
        write_pair(f"nocompress_{bp}", nc.astype(dt),
                   {"compression_type": "NOCOMPRESS", "quantize_level": 16,
                    "quantize_method": SUBTRACTIVE_DITHER_1, "dither_seed": 9},
                   raw_decoded=True)

    # ---------- Refused: a 7-axis image (comp_bad/) ----------
    # astropy writes ZNAXIS = 7; cfitsio's MAX_COMPRESS_DIM is 6 and it
    # refuses the HDU (BAD_NAXIS, imcompress.c:5370-5379), so the library
    # refuses it too. No raw twin: it is a refusal sample.
    d7 = np.arange(20, dtype=np.int16).reshape((1, 1, 1, 1, 1, 4, 5))
    p7 = os.path.join(TESTDATA, "comp_bad", "gzip1_7d.fits")
    fits.HDUList([fits.PrimaryHDU(), fits.CompImageHDU(data=d7, compression_type="GZIP_1")]).writeto(p7, overwrite=True, checksum=False)
    print(f"  {'comp_bad/gzip1_7d':30s} comp={os.path.getsize(p7):8d}")


if __name__ == "__main__":
    main()
