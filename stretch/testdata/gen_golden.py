#!/usr/bin/env python3
"""Generate golden stretch fixtures using siril-cli as the oracle.

Creates a synthetic FITS image, runs each of the 18 fitview stretch
presets through siril-cli, and captures the stretched pixel values as
golden reference JSON files. The Go test loads these and diffs against
the pure-Go stretch implementation.

Usage:
    python gen_golden.py <siril_cli_path> <output_dir>
"""

import json
import math
import os
import subprocess
import sys
import tempfile

import numpy as np
from astropy.io import fits as pyfits

SIRIL_CLI = sys.argv[1] if len(sys.argv) > 1 else "siril-cli"
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else "golden"
os.makedirs(OUTDIR, exist_ok=True)

# The 18 presets from fitview/siril.go, exactly.
PRESETS = [
    ("kstars1",        "autostretch -2.8 0.25"),
    ("kstars2",        "autostretch -1.5 0.25"),
    ("kstars3",        "autostretch -1.0 0.125"),
    ("kstars4",        "autostretch -0.5 0.125"),
    ("kstars5",        "autostretch -2.0 0.125"),
    ("kstars6",        "autostretch -4.0 0.125"),
    ("kstars7",        "autostretch -5.5 0.25"),
    ("custom_dark",    "autostretch -6.5 0.03"),
    ("custom_extreme", "autostretch -7.5 0.02"),
    ("custom_ultra",   "autostretch -9.0 0.01"),
    ("custom_mask",    "autostretch -12.0 0.005"),
    ("auto_linked",    "autostretch -linked"),
    ("asinh_light",    "asinh 1 0"),
    ("asinh_medium",   "asinh 50 0"),
    ("asinh_strong",   "asinh 500 0"),
    ("ght_mild",       "ght -D=2.0 -B=0.0 -SP=0.5"),
    ("ght_strong",     "ght -D=5.0 -B=2.0 -SP=0.3"),
    ("linear",         "linstretch -BP=0.05"),
    # CLAHE omitted for now — requires OpenCV-specific behavior
]


def make_synthetic_fits(path, nx=64, ny=48):
    """Create a synthetic 3-channel float32 FITS image.

    The image has:
    - A smooth gradient in x (to test stretching across a range)
    - Gaussian noise (to test statistics computation)
    - A few exact zeros (to test black point handling)
    - Values in [0, 1] (normalized, like a typical calibrated astro image)
    """
    rng = np.random.RandomState(42)
    nch = 3
    data = np.zeros((nch, ny, nx), dtype=np.float32)

    for ch in range(nch):
        # Gradient: 0 at left, ~0.3 at right (typical faint nebula range)
        gradient = np.linspace(0.001, 0.3, nx).reshape(1, nx)
        gradient = np.broadcast_to(gradient, (ny, nx)).copy()

        # Per-channel offset so R, G, B have slightly different medians.
        gradient += 0.02 * (ch - 1)

        # Add Gaussian noise (sigma ~ 0.01, typical CCD noise)
        noise = rng.normal(0, 0.01, (ny, nx)).astype(np.float32)
        data[ch] = np.clip(gradient + noise, 0, 1)

    # Scatter some exact zeros.
    data[:, 0, 0] = 0
    data[:, ny-1, nx-1] = 0

    hdu = pyfits.PrimaryHDU(data)
    hdu.header['BITPIX'] = -32
    hdu.writeto(path, overwrite=True)
    return data


def run_siril_stretch(siril_cli, fits_path, command, out_fits_path, tmpdir):
    """Run a Siril stretch command and save the result as FITS."""
    fits_dir = os.path.dirname(os.path.abspath(fits_path))
    fits_base = os.path.splitext(os.path.basename(fits_path))[0]
    out_dir = os.path.dirname(os.path.abspath(out_fits_path))
    out_base = os.path.splitext(os.path.basename(out_fits_path))[0]

    script = f"""requires 1.3.6
cd {fits_dir}
load {fits_base}
{command}
cd {out_dir}
save {out_base}
close
"""
    script_path = os.path.join(tmpdir, "stretch.ssf")
    with open(script_path, "w") as f:
        f.write(script)

    result = subprocess.run(
        [siril_cli, "-s", script_path],
        capture_output=True, text=True, timeout=60
    )
    if result.returncode != 0:
        print(f"  WARNING: siril-cli failed for '{command}':")
        print(result.stderr[-500:] if result.stderr else "(no stderr)")
        return False
    return True


def read_fits_pixels(path):
    """Read a FITS file and return pixel data as a numpy array."""
    with pyfits.open(path) as hdul:
        return hdul[0].data.astype(np.float64)


def main():
    print(f"Generating stretch golden fixtures with siril-cli: {SIRIL_CLI}")
    print(f"Output: {os.path.abspath(OUTDIR)}/")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create synthetic input.
        input_path = os.path.join(tmpdir, "synthetic.fits")
        input_data = make_synthetic_fits(input_path)

        # Save the input pixels as golden reference.
        input_golden = {
            "shape": list(input_data.shape),
            "pixels": input_data.flatten().tolist(),
        }
        with open(os.path.join(OUTDIR, "input.json"), "w") as f:
            json.dump(input_golden, f, indent=2)
        print(f"  input: shape={input_data.shape}")

        # Run each preset.
        for name, command in PRESETS:
            out_path = os.path.join(tmpdir, f"out_{name}.fits")
            ok = run_siril_stretch(SIRIL_CLI, input_path, command, out_path, tmpdir)
            if not ok:
                print(f"  {name}: FAILED (skipped)")
                continue

            # Siril saves as .fit by default when using save command
            # Try both extensions.
            actual_path = out_path
            if not os.path.exists(actual_path):
                actual_path = out_path.replace(".fits", ".fit")
            if not os.path.exists(actual_path):
                print(f"  {name}: output not found at {out_path} or {actual_path}")
                continue

            out_data = read_fits_pixels(actual_path)
            golden = {
                "name": name,
                "command": command,
                "shape": list(out_data.shape),
                "pixels": out_data.flatten().tolist(),
            }
            golden_path = os.path.join(OUTDIR, f"stretch_{name}.json")
            with open(golden_path, "w") as f:
                json.dump(golden, f)  # no indent — large files
            print(f"  {name}: shape={out_data.shape} "
                  f"range=[{out_data.min():.6f}, {out_data.max():.6f}]")

    nfiles = len([f for f in os.listdir(OUTDIR) if f.endswith(".json")])
    print(f"Done: {nfiles} golden files written.")


if __name__ == "__main__":
    main()
