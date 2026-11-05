"""Compact, exact golden encoding shared by the stretch golden generators.

Mirrored by stretch/golden_test.go; see PROVENANCE.md, "Golden format".

An array (an image, a vector set, a lookup table) is stored as

    {"shape": [...], "sha256": "<hex>", "bits": "<hex>"}                 (<= 1024 values)
    {"shape": [...], "sha256": "<hex>", "idx": [...], "bits": "<hex>"}  (larger)

sha256 is taken over the canonical encoding of EVERY element in row-major
(channel-major for images) order: float32 as its big-endian IEEE bits with
every NaN written as 0x7fc00000 (so -0 and +0 differ, NaN payloads do not);
uint16 as big-endian. "bits" holds the same canonical words in hex (8 hex
digits per float32, 4 per uint16) for every element (small arrays) or for
the elements at the sorted indices "idx": the first and last 8, the four
corners of every plane (3-D shapes [1 or 3, ny, nx]), the first minimum and maximum (NaN ignored), the
first 8 NaNs and 64 evenly strided elements. A test fails on any hash
mismatch; the sample only makes the failure message point at pixels.

Inputs: an input is either inline ({"shape", "sha256", "bits"} with every
value), a committed raw file ({"file": "inputs/<name>.f32"}, little-endian,
channel-major) optionally transformed by "ops" (see derive()), or a constant
({"const": "<hex>"}). Every input record carries the sha256 of the final
array so the Go test proves it rebuilt the same input.
"""

import gzip
import hashlib
import json
import os
import struct

import numpy as np

NAN32 = 0x7FC00000
INLINE_MAX = 1024


def canon32(a):
    f = np.ascontiguousarray(a, dtype=np.float32).reshape(-1)
    u = f.view(np.uint32).copy()
    u[np.isnan(f)] = NAN32
    return u


def sha32(a):
    return hashlib.sha256(canon32(a).astype(">u4").tobytes()).hexdigest()


def sha16(a):
    return hashlib.sha256(np.ascontiguousarray(a, dtype=np.uint16).reshape(-1).astype(">u2").tobytes()).hexdigest()


def hx32(u):
    return "".join(f"{int(v):08x}" for v in u)


def hx16(u):
    return "".join(f"{int(v):04x}" for v in u)


def f32hex(v):
    """One float32 as its 8 hex digits (canonical NaN)."""
    return hx32(canon32(np.array([v], np.float32)))


def f64hex(v):
    """float64 -> big-endian IEEE bits in hex (16 digits)."""
    return struct.pack(">d", v).hex()


def unhex32(s):
    return np.array([int(s[i:i + 8], 16) for i in range(0, len(s), 8)], np.uint32).view(np.float32)


def sample_indices(a, shape):
    """Deterministic sample (module doc). a: flat float array (or ints)."""
    n = a.size
    idx = set(range(min(8, n))) | set(range(max(0, n - 8), n))
    if len(shape) == 3 and shape[0] in (1, 3):  # image planes [nch, ny, nx]
        ny, nx = shape[1], shape[2]
        plane = ny * nx
        for c in range(n // plane):
            b = c * plane
            idx |= {b, b + nx - 1, b + (ny - 1) * nx, b + plane - 1}
    fa = a.astype(np.float64)
    finite = ~np.isnan(fa)
    if finite.any():
        w = np.where(finite)[0]
        idx.add(int(w[np.argmin(fa[w])]))
        idx.add(int(w[np.argmax(fa[w])]))
    idx |= set(int(i) for i in np.where(~finite)[0][:8])
    idx |= {k * n // 64 for k in range(64)}
    return sorted(i for i in idx if 0 <= i < n)


def arr32(a, shape=None):
    a = np.asarray(a, dtype=np.float32)
    shape = list(shape if shape is not None else a.shape)
    flat = a.reshape(-1)
    u = canon32(flat)
    rec = {"shape": shape, "sha256": sha32(flat)}
    if flat.size <= INLINE_MAX:
        rec["bits"] = hx32(u)
    else:
        idx = sample_indices(flat, shape)
        rec["idx"] = idx
        rec["bits"] = hx32(u[idx])
    return rec


def arr16(a, shape=None):
    a = np.asarray(a, dtype=np.uint16)
    shape = list(shape if shape is not None else a.shape)
    flat = a.reshape(-1)
    rec = {"shape": shape, "sha256": sha16(flat)}
    if flat.size <= INLINE_MAX:
        rec["bits"] = hx16(flat)
    else:
        idx = sample_indices(flat, shape)
        rec["idx"] = idx
        rec["bits"] = hx16(flat[idx])
    return rec


def inline_input(a):
    """An input stored with every value (any size the caller allows)."""
    a = np.asarray(a, dtype=np.float32)
    return {"shape": list(a.shape), "sha256": sha32(a), "bits": hx32(canon32(a))}


class Inputs:
    """Committed raw input files under <out>/inputs/ (little-endian)."""

    def __init__(self, out):
        self.dir = os.path.join(out, "inputs")
        os.makedirs(self.dir, exist_ok=True)
        self.arrays = {}

    def add(self, name, a):
        a = np.ascontiguousarray(a)
        path = os.path.join(self.dir, name)
        dt = "<f4" if a.dtype == np.float32 else "<u2"
        a.astype(dt).tofile(path)
        self.arrays[name] = a
        return a

    def ref(self, name, ops=None):
        """Input record: committed file `name`, then `ops` (derive())."""
        base = self.arrays[name]
        rec = {"file": "inputs/" + name, "file_shape": list(base.shape)}
        a = base
        if ops:
            rec["ops"] = ops
            a = derive(base, ops)
        rec["shape"] = list(a.shape)
        rec["sha256"] = sha32(a) if a.dtype == np.float32 else sha16(a)
        return rec, a


def const_input(value, shape):
    a = np.full(shape, np.float32(value), np.float32)
    return {"const": f32hex(value), "shape": list(shape), "sha256": sha32(a)}, a


def derive(a, ops):
    """Apply input transformations (mirrored by golden_test.go derive):

    ["channels", [c...]]          keep these channels, in this order
    ["fill", c, hex]              set every pixel of channel c
    ["nanpattern"]                NaN where (i + 4c) % 11 == 0, and at pixel 100 of every channel
    ["crop", y0, y1, x0, x1]      rows y0..y1-1, columns x0..x1-1
    ["mul", hex] / ["add", hex]   every pixel times / plus a float32 (one float32 rounding)
    ["put", c, [i...], [hex...]]  pixel i (flat index in the plane) of channel c; c = -1: every channel
    """
    a = np.array(a, dtype=np.float32, copy=True)
    for op in ops:
        k = op[0]
        if k == "channels":
            a = a[op[1]].copy()
        elif k == "fill":
            a[op[1]] = unhex32(op[2])[0]
        elif k == "nanpattern":
            nch = a.shape[0]
            n = a.shape[1] * a.shape[2]
            idx = np.arange(n)
            for c in range(nch):
                flat = a[c].reshape(-1)
                flat[(idx + 4 * c) % 11 == 0] = np.nan
                flat[100] = np.nan
        elif k == "crop":
            y0, y1, x0, x1 = op[1:]
            a = a[:, y0:y1, x0:x1].copy()
        elif k == "mul":
            a = (a * unhex32(op[1])[0]).astype(np.float32)
        elif k == "add":
            a = (a + unhex32(op[1])[0]).astype(np.float32)
        elif k == "put":
            c, ii, vv = op[1], op[2], unhex32("".join(op[3]))
            chans = range(a.shape[0]) if c < 0 else [c]
            for ch in chans:
                flat = a[ch].reshape(-1)
                flat[ii] = vv
        else:
            raise ValueError(op)
    return a


def dump(path, obj, gz=False):
    """Compact JSON, one trailing newline (byte-for-byte reproducible).
    With gz the JSON is written to path + ".gz" (gzip level 9, mtime 0):
    the Go tests read it through testutil.ReadGolden and the reproduce
    test compares the decompressed bytes."""
    text = (json.dumps(obj, separators=(",", ":")) + "\n").encode()
    if gz:
        with open(path + ".gz", "wb") as f:
            f.write(gzip.compress(text, compresslevel=9, mtime=0))
        return
    with open(path, "wb") as f:
        f.write(text)
