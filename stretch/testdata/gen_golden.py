#!/usr/bin/env python3
"""Check a real siril-cli against the stretch goldens; write its argument table.

    python gen_golden.py <siril_cli>               # check everything below; exit 1 on a difference
    python gen_golden.py <siril_cli> --write-args  # rewrite golden/siril_cli_args.json

1. The 18 fitview presets on the synthetic image: the canonical sha256 of
   each siril-cli output (goldenfmt.py: big-endian float32 bits, NaN as
   0x7fc00000) must equal golden/siril_cli_sha256.json.
2. Every oracle golden golden/sref/<case>.json that names a siril command
   ("command"): its input is rebuilt, written as a float32 FITS file, run
   through siril-cli, and the output digest must equal the golden's.
   siril-cli's FITS loader divides float data whose maximum exceeds 10 by
   that maximum (io/image_format_fits.c:705, :794), so a +-Inf pixel
   cannot be loaded as is. For an input with +-Inf (the inf_* cases) the
   file holds markers +-7 at those pixels (values the input does not
   have) and the script runs `pm "iif($T==7, 1e39, iif($T==-7, -1e39,
   $T))"` before the command: pixel math evaluates in double and stores
   (float)1e39 = +Inf; every other pixel passes through unchanged, which
   the script checks by saving the image after pm alone.
3. golden/siril_cli_args.json: for each command in ARG_CASES, whether
   siril-cli accepts it (exit status 0 and an output file), Siril's error
   message if not, and the output pixels if so (or "out_varies" when two
   runs give different pixels), on a small mono image.
   TestSirilCLIArgs compares the Go functions' errors and pixels with it.

History: this script, run with a siril-cli whose version was not recorded
("requires 1.3.6" below is only Siril's minimum-version directive), wrote
the original pixel files golden/stretch_<preset>.json (commit c17580b).
Those files are now digests written by gen_siril_goldens.py from the
verbatim Siril 1.4.0 oracle (siril_cli_digests.py keeps the siril-cli
digests from git history). siril-cli 1.4.0 built from the release tarball
(build_siril_cli.sh) reproduces all of them; see PROVENANCE.md.
"""

import json
import os
import subprocess
import sys
import tempfile

import numpy as np
from astropy.io import fits as pyfits

SIRIL_CLI = sys.argv[1] if len(sys.argv) > 1 else "siril-cli"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from goldenfmt import arr32, canon32, derive, dump, sha32, unhex32  # noqa: E402

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
    # CLAHE goldens: gen_siril_goldens.py (verbatim Siril + OpenCV) and
    # gen_clahe_cv2.py (opencv-python)
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


def run_siril_stretch(siril_cli, fits_path, command, out_fits_path, tmpdir, pre=""):
    """Run a Siril stretch command and save the result as FITS. Returns
    (ok, log): ok is False when siril-cli refuses the command (non-zero
    exit or no output file)."""
    fits_dir = os.path.dirname(os.path.abspath(fits_path))
    fits_base = os.path.splitext(os.path.basename(fits_path))[0]
    out_dir = os.path.dirname(os.path.abspath(out_fits_path))
    out_base = os.path.splitext(os.path.basename(out_fits_path))[0]
    for ext in (".fit", ".fits"):
        if os.path.exists(os.path.join(out_dir, out_base + ext)):
            os.remove(os.path.join(out_dir, out_base + ext))

    script = f"""requires 1.3.6
cd {fits_dir}
load {fits_base}
{pre}{command}
cd {out_dir}
save {out_base}
close
"""
    script_path = os.path.join(os.path.abspath(tmpdir), "stretch.ssf")
    with open(script_path, "w") as f:
        f.write(script)

    result = subprocess.run(
        [siril_cli, "-s", script_path],
        capture_output=True, text=True, timeout=60
    )
    out = result.stdout + result.stderr
    produced = any(os.path.exists(os.path.join(out_dir, out_base + e)) for e in (".fit", ".fits"))
    return result.returncode == 0 and produced, out


def read_fits_pixels(path):
    """Read a FITS file and return pixel data as a numpy array."""
    if not os.path.exists(path):
        path = path.replace(".fits", ".fit")
    with pyfits.open(path) as hdul:
        return hdul[0].data.astype(np.float64)


def write_fits(path, data):
    """data: (nch, ny, nx) float32; a mono image is written as 2-D."""
    d = data if data.shape[0] == 3 else data[0]
    pyfits.PrimaryHDU(np.ascontiguousarray(d, dtype=">f4")).writeto(path, overwrite=True)


INF_MARK = 7.0
INF_PM = 'pm "iif($T==7, 1e39, iif($T==-7, -1e39, $T))"\n'


def run_cli(data, command, tmpdir):
    """(ok, output (nch, ny, nx) or None, log). +-Inf pixels go in as the
    markers +-7 and are restored by pixel math inside siril-cli (module
    doc); the restore is checked bit for bit on its own first."""
    inp = os.path.join(tmpdir, "in.fits")
    out = os.path.join(tmpdir, "out.fits")
    pre = ""
    if np.isinf(data).any():
        if (np.abs(data[np.isfinite(data)]) == INF_MARK).any() or np.nanmax(np.where(np.isinf(data), 0, data)) > 10:
            sys.exit(f"{command}: the input has the Inf marker or values above 10")
        marked = np.where(np.isposinf(data), INF_MARK, np.where(np.isneginf(data), -INF_MARK, data)).astype(np.float32)
        write_fits(inp, marked)
        ok, log = run_siril_stretch(SIRIL_CLI, inp, "", out, tmpdir, pre=INF_PM)
        if not ok or sha32(read_fits_pixels(out).astype(np.float32).reshape(data.shape)) != sha32(data):
            sys.exit(f"{command}: pixel math did not restore the +-Inf input exactly\n{log[-500:]}")
        pre = INF_PM
    else:
        write_fits(inp, data)
    ok, log = run_siril_stretch(SIRIL_CLI, inp, command, out, tmpdir, pre=pre)
    if not ok:
        return False, None, log
    return True, read_fits_pixels(out).astype(np.float32).reshape(data.shape), log


def rebuild_input(rec):
    """An sref golden's input array (goldenfmt.py input records)."""
    if "const" in rec:
        return np.full(rec["shape"], unhex32(rec["const"])[0], np.float32)
    if "file" in rec:
        a = np.fromfile(os.path.join(HERE, "golden", rec["file"]), "<f4").reshape(rec["file_shape"])
        return derive(a, rec["ops"]) if rec.get("ops") else a
    return unhex32(rec["bits"]).reshape(rec["shape"])


def check_presets(tmpdir):
    want = json.load(open(os.path.join(HERE, "golden", "siril_cli_sha256.json")))
    bad = 0
    input_path = os.path.join(tmpdir, "synthetic.fits")
    input_data = make_synthetic_fits(input_path)
    if sha32(input_data) != want["input"]["sha256"]:
        sys.exit("synthetic input does not match the committed input digest")
    for name, command in PRESETS:
        out_path = os.path.join(tmpdir, f"out_{name}.fits")
        ok, log = run_siril_stretch(SIRIL_CLI, input_path, command, out_path, tmpdir)
        if not ok:
            print(f"  {name}: siril-cli FAILED\n{log[-500:]}")
            bad += 1
            continue
        got = read_fits_pixels(out_path).astype(np.float32)
        same = sha32(got) == want["presets"][name]["sha256"]
        print(f"  preset {name}: {'identical' if same else 'DIFFERENT'}")
        bad += not same
    return bad


def check_sref(tmpdir):
    names = open(os.path.join(HERE, "golden", "sref", "MANIFEST")).read().split()
    bad = n = ninf = 0
    for name in names:
        c = json.load(open(os.path.join(HERE, "golden", "sref", name + ".json")))
        if not c.get("command"):
            continue
        n += 1
        data = rebuild_input(c["input"])
        if sha32(data) != c["input"]["sha256"]:
            sys.exit(f"{name}: rebuilt input does not match its digest")
        ok, got, log = run_cli(data, c["command"], tmpdir)
        if not ok:
            print(f"  {name}: siril-cli refused {c['command']!r}\n{log[-500:]}")
            bad += 1
            continue
        same = sha32(got) == c["out"]["sha256"]
        inf = " (+-Inf input through pixel math)" if np.isinf(data).any() else ""
        ninf += bool(inf)
        print(f"  {name}: {'identical' if same else 'DIFFERENT'} ({c['command']}){inf}")
        bad += not same
    print(f"  {n} oracle goldens with a siril command, all compared bit for bit ({ninf} with +-Inf input)")
    return bad


# Commands whose acceptance (and output) siril-cli decides; on ARG_INPUT.
THR = repr(float(np.float32(1e-3)))
ARG_CASES = [
    "linstretch -BP=-0.5", "linstretch -BP=-0", "linstretch -BP=0", "linstretch -BP=0.05",
    "linstretch -BP=0.999", "linstretch -BP=1", "linstretch -BP=1.0000001", "linstretch -BP=1.5",
    "linstretch -BP=nan", "linstretch -BP=inf",
    "autostretch", "autostretch -2.8 0", "autostretch -2.8 -0", "autostretch -2.8 0.25", "autostretch -2.8 1",
    "autostretch -2.8 1.00000001", "autostretch -2.8 1.0000001", "autostretch -2.8 -0.01", "autostretch -2.8 1.5",
    "autostretch -2.8 nan", "autostretch -2.8 inf", "autostretch nan 0.25", "autostretch inf 0.25",
    "autostretch -inf 0.25", "autostretch 5 0.25", "autostretch -linked", "autostretch -linked -2.8 0.25",
    "autostretch -linked -2.8 1.5", "autostretch -linked -2.8 -0.5", "autostretch -linked nan 0.25",
    "ght -D=2 -SP=0.5", "ght -D=0 -SP=0.5", "ght -D=-1 -SP=0.5", "ght -D=10 -SP=0.5", "ght -D=10.000001 -SP=0.5",
    "ght -D=nan -SP=0.5", "ght -D=inf -SP=0.5", "ght -D=2 -B=-0.5 -SP=0.5", "ght -D=2 -B=15 -SP=0.5",
    "ght -D=2 -B=15.01 -SP=0.5", "ght -D=2 -B=nan -SP=0.5", "ght -D=2 -B=0.0005 -SP=0.5",
    "ght -D=2 -B=0.001 -SP=0.5", f"ght -D=2 -B={THR} -SP=0.5",
    "ght -D=2 -SP=1.01", "ght -D=2 -SP=-0.1", "ght -D=2 -SP=nan", "ght -D=2 -LP=0.3 -SP=0.2",
    "ght -D=2 -LP=0.2 -SP=0.2", "ght -D=2 -LP=-0.1 -SP=0.2", "ght -D=2 -LP=nan -SP=0.2",
    "ght -D=2 -SP=0.5 -HP=0.4", "ght -D=2 -SP=0.5 -HP=0.5", "ght -D=2 -SP=0.5 -HP=1.01", "ght -D=2 -SP=0.5 -HP=nan",
    "ght -D=1 -SP=0 -HP=0", "ght -D=3 -B=2 -SP=0 -HP=0", "ght -D=1 -SP=0.1 -HP=0",
    # asinh: siril-cli refuses only beta < 1 (command.c:2655); any offset.
    "asinh 1", "asinh 10 0.1", "asinh 0.999", "asinh 0", "asinh -1", "asinh nan", "asinh inf",
    "asinh 1e30", "asinh 10 -0.1", "asinh 10 0.99", "asinh 10 1", "asinh 10 1.5", "asinh 10 nan",
    "asinh 10 inf", "asinh 10 -inf",
    # mtf: shadows in [0, 1), midtones in (0, 1), highlights in (0, 1] (command.c:3146).
    "mtf 0 0.5 1", "mtf 0.01 0.2 0.8", "mtf 0.5 0.3 0.4", "mtf 0.3 0.6 0.3", "mtf 0 0 1", "mtf 0 1 1",
    "mtf 0 -0.5 1", "mtf 0 1.5 1", "mtf -0.1 0.5 1", "mtf 1 0.5 1", "mtf 0 0.5 0", "mtf 0 0.5 1.1",
    "mtf nan 0.5 1", "mtf 0 nan 1", "mtf 0 0.5 nan", "mtf 0 inf 1", "mtf -inf 0.5 1", "mtf 0 0.5 inf",
    # clahe: clip <= 0 and tile < 1 refused (command.c:2702-2711). On this
    # 8x4 image with tile 1, tileSizeTotal is 32: clip 2^42 gives
    # clip*32/65536 = 2^31, which OpenCV's int conversion turns into
    # INT_MIN, i.e. limit 1; 2^42 - 2048 gives 2^31 - 1.
    "clahe 2 1", "clahe 2 2", "clahe 0 1", "clahe -1 1", "clahe 2 0", "clahe nan 1", "clahe inf 1",
    "clahe 1e30 1", "clahe 4398046511104 1", "clahe 4398046509056 1", "clahe 1e5 1",
]


def arg_input():
    """A small mono image: 0, 1, -0, four equal pixels (so that CLAHE's
    clip limit 1 changes the output) and a spread of values in between."""
    v = np.concatenate([[0.0, 1.0, -0.0, 0.5, 0.25, 0.25, 0.25, 0.25], np.linspace(0.003, 0.97, 24)]).astype(np.float32)
    return v.reshape(1, 4, 8)


def siril_error(log):
    lines = [l.split("log: ", 1)[1] for l in log.splitlines() if l.startswith("log: ") and
             ("rror" in l or "must" in l or "nvalid" in l or "range" in l)]
    return lines[0].strip() if lines else ""


def arg_table(tmpdir):
    data = arg_input()
    rows = []
    for command in ARG_CASES:
        ok, got, log = run_cli(data, command, tmpdir)
        row = {"command": command, "accepted": ok}
        if ok:
            # Run twice: with a NaN parameter Siril can read uninitialised
            # coefficients (GHTsetup), and the pixels then vary run to run.
            again = run_cli(data, command, tmpdir)[1]
            if sha32(again) == sha32(got):
                row["out"] = arr32(got)
            else:
                row["out_varies"] = True
        else:
            row["error"] = siril_error(log)
        rows.append(row)
        print(f"  {command}: {('accepted' + (' (output varies run to run)' if row.get('out_varies') else '')) if ok else 'refused: ' + row['error']}")
    # The clahe rows must tell OpenCV's clip conversion apart: 2^42 acts
    # as limit 1 (as clip 2), 2^42 - 2048 as no limit (as clip 1e5).
    out = {r["command"]: r.get("out", {}).get("sha256") for r in rows}
    if not (out["clahe 4398046511104 1"] == out["clahe 2 1"] == out["clahe inf 1"] == out["clahe 1e30 1"]
            != out["clahe 4398046509056 1"] == out["clahe 1e5 1"] == out["clahe nan 1"]):
        sys.exit("clahe argument rows do not separate limit 1 from no limit")
    return {"siril_cli": cli_version(), "input": arr32(data), "cases": rows}


def cli_version():
    r = subprocess.run([SIRIL_CLI, "--version"], capture_output=True, text=True)
    return [l for l in (r.stdout + r.stderr).splitlines() if l.startswith("siril ")][0].strip()


def main():
    print(f"Checking siril-cli {SIRIL_CLI} ({cli_version()})")
    path = os.path.join(HERE, "golden", "siril_cli_args.json")
    with tempfile.TemporaryDirectory() as tmpdir:
        if "--write-args" in sys.argv:
            dump(path, arg_table(tmpdir))
            print(f"wrote {path}")
            return
        bad = check_presets(tmpdir)
        bad += check_sref(tmpdir)
        table = arg_table(tmpdir)
        want = json.load(open(path))
        if table != want:
            print("  golden/siril_cli_args.json differs from this siril-cli")
            bad += 1
    if bad:
        sys.exit(f"{bad} differences from the committed goldens")
    print("siril-cli agrees with every committed golden.")


if __name__ == "__main__":
    main()
