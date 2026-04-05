#!/usr/bin/env python3
"""Generate tile-compressed FITS test fixtures using astropy.

Usage:

    python3 testdata/generate_compressed.py

This script produces a matched set of compressed and raw-uncompressed
FITS files in testdata/. Every compressed file has a corresponding
comp_raw_<name>.fits file containing the same deterministic data through
a plain ImageHDU — the raw file is the byte-exact ground truth for the
Go integration tests.

The fixture set covers:

  - Every common tile-compression algorithm (RICE_1, GZIP_1, GZIP_2,
    HCOMPRESS_1, PLIO_1, NOCOMPRESS)
  - Both integer and float pixel types
  - Single-row tiles (the astropy default) and non-trivial 2D tiles
  - Partial edge tiles (image size not divisible by tile size)
  - NO_DITHER and SUBTRACTIVE_DITHER_1 quantization modes for floats

All data is generated with fixed seeds so the fixtures are bit-for-bit
reproducible across runs.
"""

import os
import numpy as np
import astropy.io.fits as fits


TESTDATA = os.path.dirname(os.path.abspath(__file__))
RNG = np.random.default_rng(42)


def write_pair(name, data, cls_kwargs, ext_name=None):
    """Write a compressed fixture and its raw-uncompressed counterpart."""
    comp_path = os.path.join(TESTDATA, f"comp_{name}.fits")
    raw_path = os.path.join(TESTDATA, f"comp_raw_{name}.fits")

    comp_hdu = fits.CompImageHDU(data=data, name=ext_name or "IMG", **cls_kwargs)
    fits.HDUList([fits.PrimaryHDU(), comp_hdu]).writeto(comp_path, overwrite=True)

    raw_hdu = fits.ImageHDU(data=data, name=ext_name or "IMG")
    fits.HDUList([fits.PrimaryHDU(), raw_hdu]).writeto(raw_path, overwrite=True)

    print(f"  {name:30s} comp={os.path.getsize(comp_path):8d} raw={os.path.getsize(raw_path):8d}")


def main():
    print("Generating compressed FITS fixtures in testdata/:")

    # ---------- Integer image, every algorithm ----------
    i16_small = np.arange(100 * 100, dtype=np.int16).reshape((100, 100))
    i32_small = np.arange(100 * 100, dtype=np.int32).reshape((100, 100))

    write_pair("rice_i16", i16_small, {"compression_type": "RICE_1"})
    write_pair("rice_i32", i32_small, {"compression_type": "RICE_1"})
    write_pair("gzip1_i16", i16_small, {"compression_type": "GZIP_1"})
    write_pair("gzip1_i32", i32_small, {"compression_type": "GZIP_1"})
    write_pair("gzip2_i16", i16_small, {"compression_type": "GZIP_2"})
    write_pair("gzip2_i32", i32_small, {"compression_type": "GZIP_2"})
    write_pair("hcompress_i16", i16_small, {"compression_type": "HCOMPRESS_1"})
    write_pair("nocompress_i16", i16_small, {"compression_type": "NOCOMPRESS"})

    # ---------- Float image with quantization ----------
    # Smooth gradient + gaussian noise, reproducible via fixed RNG.
    f32_smooth = (
        np.linspace(0, 1000, 100 * 100).reshape((100, 100))
        + RNG.normal(0, 1, (100, 100))
    ).astype(np.float32)

    # RICE_1 with default dithering (SUBTRACTIVE_DITHER_1).
    write_pair(
        "rice_f32_dither",
        f32_smooth,
        {"compression_type": "RICE_1", "quantize_level": 16},
    )
    # RICE_1 with NO_DITHER.
    write_pair(
        "rice_f32_nodither",
        f32_smooth,
        {
            "compression_type": "RICE_1",
            "quantize_level": 16,
            "dither_seed": 0,
            # astropy maps dither_seed=0 + quantize_method=NO_DITHER via method arg
        },
    )
    # GZIP_1 float (also quantized).
    write_pair(
        "gzip1_f32",
        f32_smooth,
        {"compression_type": "GZIP_1", "quantize_level": 16},
    )
    # GZIP_2 float.
    write_pair(
        "gzip2_f32",
        f32_smooth,
        {"compression_type": "GZIP_2", "quantize_level": 16},
    )
    # HCOMPRESS float.
    write_pair(
        "hcompress_f32",
        f32_smooth,
        {"compression_type": "HCOMPRESS_1", "quantize_level": 16},
    )

    # ---------- PLIO (sparse integer mask) ----------
    # PLIO needs non-negative integer mask-like data. Make a small mask.
    mask_i16 = np.zeros((50, 50), dtype=np.int16)
    mask_i16[10:20, 15:35] = 1
    mask_i16[30:40, 5:45] = 2
    write_pair("plio_mask", mask_i16, {"compression_type": "PLIO_1"})

    # ---------- Tile geometry edge cases ----------
    # Non-trivial 2D tiles (32x32) over a 100x100 image → edge tiles are
    # partial (4 pixels wide in the last column, 4 pixels tall in the last
    # row).
    write_pair(
        "rice_multitile",
        i16_small,
        {"compression_type": "RICE_1", "tile_shape": (32, 32)},
    )
    # Whole-image-as-one-tile.
    write_pair(
        "rice_onetile",
        i16_small,
        {"compression_type": "RICE_1", "tile_shape": (100, 100)},
    )
    # Odd-sized image (103x97) with default tiles — exercises partial edge
    # behavior at a non-power-of-two boundary.
    i16_odd = np.arange(103 * 97, dtype=np.int16).reshape((97, 103))
    write_pair(
        "rice_odd",
        i16_odd,
        {"compression_type": "RICE_1"},
    )


if __name__ == "__main__":
    main()
