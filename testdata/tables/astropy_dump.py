"""Dump the numeric table cells of a FITS file as astropy reads them.

Usage:
  python astropy_dump.py <file.fits>...      JSON list, one column list per file
  python astropy_dump.py --golden <file.fits> one golden object (golden_astropy/)

The Go suite uses astropy as a second, independent reader next to
cfitsio. The committed goldens (golden_astropy/<sample>.json, written by
`make astropy-goldens`, checked byte for byte by
TestAstropyGoldensReproduce) are what TestTablesMatchAstropy compares
against, so the default `go test` needs no Python; the list mode is used at
test time for files the Go tests write themselves.

Each column entry is
  {"hdu": i, "col": n, "name": ..., "kind": "int"|"uint"|"float"|"complex"|"str",
   "tnull": TNULLn or null, "rows": [[v, ...], ...]}
for every numeric column (B I J K E D C M and their P/Q forms, ASCII
I F E D) and every fixed character column (binary rA, ASCII Aw; "str":
the strings astropy returns, one per cell or one per element of a TDIM
string array). Floats are written with float.hex() so the comparison is
bit-exact ("nan" for NaN, "inf"/"-inf"); integers as JSON integers.
astropy applies TSCAL/TZERO; it does not mask TNULL values.

Golden mode stores a column of more than DIGEST_MIN elements as an exact
digest plus samples (see digest below); it also records astropy's HDU
classes, and instead of failing on
the first problem it records it: a column astropy cannot read becomes
{"hdu", "col", "name", "error": "<ExceptionType>: <message>"}, and a file
astropy cannot open at all {"error": ...}. A TFORM astropy rejects outright
(1PX / 1QX: "Invalid column format") is patched in memory to the
descriptor of the same size (1PB / 1QB), that column is left out and the
patch is listed under "patched", so every other column of the file is
still dumped.
"""

import hashlib
import io
import json
import re
import struct
import sys
import warnings

import astropy
import numpy as np
from astropy.io import fits
from astropy.io.fits.verify import VerifyError


def fhex(x):
    x = float(x)
    if x != x:
        return "nan"
    return x.hex()


def cell(v, kind):
    a = np.asarray(v).ravel()
    if kind == "complex":
        return [[fhex(z.real), fhex(z.imag)] for z in a]
    if kind == "float":
        return [fhex(z) for z in a]
    return [int(z) for z in a]


def kind_of(dt):
    if dt.kind == "c":
        return "complex"
    if dt.kind == "f":
        return "float"
    if dt.kind == "u":
        return "uint"
    if dt.kind == "i":
        return "int"
    return None


def tnull_of(col):
    v = col.null
    if v is None or v == "":
        return None
    if isinstance(v, (bool, np.bool_)):
        return None
    if isinstance(v, (int, np.integer)):
        return int(v)
    return str(v)


def column(i, n, col, data):
    fmt = str(col.format)
    if "X" in fmt or "L" in fmt:
        return None
    if "A" in fmt and ("P" in fmt or "Q" in fmt):
        return None
    values = data[col.name] if col.name else data.field(n - 1)
    base = {"hdu": i, "col": n, "name": col.name, "tnull": tnull_of(col)}
    if "A" in fmt:
        # Fixed character column: each row's strings as astropy returns
        # them (one per cell, or the elements of a TDIM string array).
        return dict(base, kind="str", rows=[[str(z) for z in np.asarray(r).ravel()] for r in values])
    if "P" in fmt or "Q" in fmt:
        rows = list(values)
        if not rows:
            return None
        k = None
        for r in rows:
            if len(r):
                k = kind_of(np.asarray(r).dtype)
                break
        k = k or "float"  # all rows empty: nothing to compare
        return dict(base, kind=k, rows=[cell(r, k) if len(r) else [] for r in rows])
    k = kind_of(values.dtype)
    if k is None:
        return None
    return dict(base, kind=k, rows=[cell(r, k) for r in values])


# A golden column with more than DIGEST_MIN elements in all its rows keeps
# its values as a digest instead of "rows": "nelem" (total elements),
# "sha256" over the canonical big-endian encoding below, and "samples",
# [flat element index, value] pairs for the first and last 16 elements and
# 64 strided ones, so a mismatch is still debuggable. Canonical encoding,
# row after row: the row's element count as int64, then per element the
# float64 bits (float; NaN as 0x7ff8000000000000), two float64 bit
# patterns (complex), int64 (int) or uint64 (uint).
DIGEST_MIN = 1024


def canon(kind, v):
    if kind == "float":
        return b"\x7f\xf8" + b"\0" * 6 if v != v else struct.pack(">d", v)
    if kind == "complex":
        return canon("float", v.real) + canon("float", v.imag)
    if kind == "uint":
        return struct.pack(">Q", int(v))
    return struct.pack(">q", int(v))


def digest(c, rawrows):
    if c["kind"] == "str":
        return c
    flat = [z for r in rawrows for z in np.asarray(r).ravel()]
    if len(flat) <= DIGEST_MIN:
        return c
    h = hashlib.sha256()
    for r in rawrows:
        a = np.asarray(r).ravel()
        h.update(struct.pack(">q", len(a)))
        for z in a:
            h.update(canon(c["kind"], complex(z) if c["kind"] == "complex" else (float(z) if c["kind"] == "float" else int(z))))
    n = len(flat)
    idx = sorted(set(list(range(16)) + list(range(n - 16, n)) + [i * n // 64 for i in range(64)]))
    del c["rows"]
    c["nelem"] = n
    c["sha256"] = h.hexdigest()
    c["samples"] = [[i, cell(flat[i], c["kind"])[0]] for i in idx]
    return c


def errtext(e):
    return "%s: %s" % (type(e).__name__, str(e).splitlines()[0] if str(e) else "")


def dump_hdul(hl, golden, skip):
    out = []
    for i, h in enumerate(hl):
        if not isinstance(h, (fits.BinTableHDU, fits.TableHDU)):
            continue
        data = h.data
        for n, col in enumerate(h.columns, start=1):
            if (i, n) in skip:
                continue
            if not golden:
                c = column(i, n, col, data)
            else:
                try:
                    c = column(i, n, col, data)
                    if c is not None:
                        c = digest(c, list(data[col.name] if col.name else data.field(n - 1)))
                except Exception as e:  # recorded, not fatal
                    c = {"hdu": i, "col": n, "name": col.name, "error": errtext(e)}
            if c is not None:
                out.append(c)
    return out


def dump(path):
    with fits.open(path) as hl:
        return dump_hdul(hl, False, set())


TFORM_RE = re.compile(rb"TFORM(\d+) *= '( *1)([PQ])X")


def golden(path):
    g = {"astropy": astropy.__version__, "numpy": np.__version__}
    raw = open(path, "rb").read()
    patched, skip = [], set()
    for _ in range(8):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                with fits.open(io.BytesIO(raw)) as hl:
                    g["hdus"] = [type(h).__name__ for h in hl]
                    for h in hl:
                        if isinstance(h, (fits.BinTableHDU, fits.TableHDU)):
                            h.data  # noqa: B018 -- builds the column definitions
                    g["columns"] = dump_hdul(hl, True, skip)
            break
        except VerifyError as e:
            cm = TFORM_RE.search(raw) if "Invalid column format: 1" in str(e) else None
            if cm is None:
                g = dict(g, error=errtext(e))
                break
            # Patch the rejected bit-array VLA to a byte VLA with the same
            # descriptor size, and leave that column out.
            off = cm.start(3) + 1
            raw = raw[:off] + b"B" + raw[off + 1:]
            hdu = raw.count(b"XTENSION= ", 0, cm.start())
            skip.add((hdu, int(cm.group(1))))
            d = cm.group(3).decode()
            patched.append("HDU %d TFORM%s 1%sX -> 1%sB (column not dumped)" % (hdu, cm.group(1).decode(), d, d))
        except Exception as e:
            g = dict(g, error=errtext(e))
            break
    if patched:
        g["patched"] = patched
    return g


if __name__ == "__main__":
    if sys.argv[1:2] == ["--golden"]:
        if len(sys.argv) != 3:
            sys.exit("usage: astropy_dump.py --golden <file.fits>")
        sys.stdout.write(json.dumps(golden(sys.argv[2]), sort_keys=True, separators=(",", ":")) + "\n")
    else:
        json.dump([dump(p) for p in sys.argv[1:]], sys.stdout)
