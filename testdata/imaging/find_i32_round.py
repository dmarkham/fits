#!/usr/bin/env python3
"""Search every BITPIX=32 value for the ones where float32(float64(u)/(2^32-1))
(double rounding) is not the float32 nearest the exact u/(2^32-1), with
u = raw + 2^31 (the full-range offset ReadFloat32 uses).

    <venv>/bin/python find_i32_round.py

numpy filters the u whose float64 quotient lies within 4 float64 ULP of a
float32 rounding midpoint; fractions.Fraction then decides each candidate
exactly. Prints the u values; gen_samples.py commits them as U32_ROUND.
"""
from fractions import Fraction

import numpy as np

SPAN = 2**32 - 1


def nearest_f32(fr):
    """float32 nearest the Fraction fr (ties to even), via exact comparison."""
    f = np.float32(float(fr))  # candidate within 1 float32 ULP
    best = None
    for c in (np.nextafter(f, np.float32(-np.inf)), f, np.nextafter(f, np.float32(np.inf))):
        d = abs(Fraction(float(c)) - fr)
        key = (d, int(np.frombuffer(np.float32(c).tobytes(), np.uint32)[0]) & 1)
        if best is None or key < best[0]:
            best = (key, c)
    return best[1]


found = []
CH = 1 << 26
for start in range(0, 1 << 32, CH):
    u = np.arange(start, start + CH, dtype=np.uint64)
    q = u.astype(np.float64) / float(SPAN)
    f = q.astype(np.float32)
    f64 = f.astype(np.float64)
    up = np.nextafter(f, np.float32(np.inf)).astype(np.float64)
    dn = np.nextafter(f, np.float32(-np.inf)).astype(np.float64)
    mid = np.where(q >= f64, (f64 + up) / 2, (f64 + dn) / 2)
    near = np.abs(q - mid) <= 4 * np.spacing(q)
    for uu in u[near]:
        uu = int(uu)
        naive = np.float32(float(uu) / float(SPAN))
        exact = nearest_f32(Fraction(uu, SPAN))
        if naive != exact:
            found.append(uu)
print(found)
