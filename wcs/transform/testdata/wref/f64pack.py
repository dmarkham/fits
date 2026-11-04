"""Compact storage for the mpmath reference tables (paper/, libm/).

A table is two files:

  <name>.json  small metadata: whatever the generator records (source,
               versions, row counts and column layout) plus
                 "payload":      "<name>.bin",
                 "values":       number of doubles in the payload,
                 "sha256":       sha256 of the UNCOMPRESSED payload,
                 "zlib_version": the zlib that compressed it;
  <name>.bin   zlib (RFC 1950, level 9) of every value as its exact
               IEEE-754 bits, big-endian float64, in the order the
               metadata describes.

Nothing is rounded: the Go reader (f64pack_test.go) recovers every double
bit for bit, and checks the length and the sha256. Until round 2 of the
proof audit the same doubles were written as decimal JSON rows (545 KB for
the three tables, now 108 KB). zlib's output for a given payload is fixed
for a given zlib version and level (Python's zlib.compress at level 9 and
wref.c's compress2 at level 9 make the same stream); the reproduce tests
compare the decoded payload and, when the local zlib is the recorded one,
the compressed bytes too.
"""
import hashlib
import json
import os
import struct
import zlib


def write(dirpath, name, meta, values):
    """Write dirpath/name.json and dirpath/name.bin; returns the json path."""
    values = [float(v) for v in values]
    payload = struct.pack(">%dd" % len(values), *values)
    m = dict(meta)
    m["payload"] = name + ".bin"
    m["values"] = len(values)
    m["sha256"] = hashlib.sha256(payload).hexdigest()
    m["zlib_version"] = zlib.ZLIB_RUNTIME_VERSION
    os.makedirs(dirpath, exist_ok=True)
    with open(os.path.join(dirpath, name + ".bin"), "wb") as f:
        f.write(zlib.compress(payload, 9))
    path = os.path.join(dirpath, name + ".json")
    with open(path, "w") as f:
        f.write(json.dumps(m, indent=1) + "\n")
    return path
