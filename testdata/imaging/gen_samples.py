#!/usr/bin/env python3
"""Write the ReadFloat32 / ReadFrame / AsImage sample matrix.

    <venv>/bin/python gen_samples.py           # in testdata/imaging

Files are written byte by byte (header cards + big-endian data) so every
BITPIX x BSCALE/BZERO x BLANK combination is exactly what the case asks
for; astropy is not used to write them (it would rescale). The real camera
crop is made by crop_camera.py; expected values by gen_expected.py.
Deterministic: re-running reproduces every byte.
"""
import os
import struct

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def card(key, value=None, comment=""):
    if key in ("COMMENT", "HISTORY", "END"):
        return f"{key:<8}{value or ''}".ljust(80)[:80]
    if isinstance(value, bool):
        v = "T" if value else "F"
        s = f"{key:<8}= {v:>20}"
    elif isinstance(value, str):
        s = f"{key:<8}= '{value:<8}'"
    elif isinstance(value, int):
        s = f"{key:<8}= {value:>20d}"
    else:
        s = f"{key:<8}= {repr(float(value)).upper():>20}"
    if comment:
        s += " / " + comment
    assert len(s) <= 80, s
    return s.ljust(80)


DTYPE = {8: ">u1", 16: ">i2", 32: ">i4", 64: ">i8", -32: ">f4", -64: ">f8"}


def hdu_bytes(bitpix, raw, shape, bscale=None, bzero=None, blank=None,
              extra=(), comment=None, raw_cards=(), xtension=False, extend=False,
              verbatim=()):
    """One HDU (header + data, each padded to 2880 bytes). raw: on-disk
    values (already offset; e.g. int16 for BITPIX 16). raw_cards: cards
    given as (key, value text) written verbatim (for deliberately invalid
    values). verbatim: whole card images (columns 1-80) written as given,
    after raw_cards (for malformed value indicators and duplicates)."""
    raw = np.asarray(raw).astype(DTYPE[bitpix])
    assert raw.size == int(np.prod(shape)) or (raw.size == 0 and 0 in shape)
    first = card("XTENSION", "IMAGE") if xtension else card("SIMPLE", True)
    cards = [first, card("BITPIX", bitpix), card("NAXIS", len(shape))]
    for i, n in enumerate(shape):
        cards.append(card(f"NAXIS{i+1}", int(n)))
    if xtension:
        cards += [card("PCOUNT", 0), card("GCOUNT", 1)]
    if extend:
        cards.append(card("EXTEND", True))
    for k, v in raw_cards:
        cards.append(f"{k:<8}= {v:>20}".ljust(80))
    for v in verbatim:
        assert len(v) <= 80, v
        cards.append(v.ljust(80))
    if bscale is not None:
        cards.append(card("BSCALE", bscale))
    if bzero is not None:
        cards.append(card("BZERO", bzero))
    if blank is not None:
        cards.append(card("BLANK", int(blank)))
    for k, v in extra:
        cards.append(card(k, v))
    if comment:
        cards.append(card("COMMENT", " " + comment))
    cards.append(card("END"))
    hdr = "".join(cards).encode("ascii")
    hdr += b" " * (-len(hdr) % 2880)
    data = raw.tobytes()
    data += b"\0" * (-len(data) % 2880)
    return hdr + data


def write(name, bitpix, raw, shape, **kw):
    with open(os.path.join(HERE, name), "wb") as f:
        f.write(hdu_bytes(bitpix, raw, shape, **kw))


rng = np.random.default_rng(20261101)

# --- BITPIX 8 -------------------------------------------------------------
write("u8_all.fits", 8, np.arange(256), (16, 16), comment="every BITPIX=8 value 0..255")
write("u8_blank.fits", 8, [0, 1, 7, 128, 254, 255, 255, 3], (8,), blank=255,
      comment="BLANK=255 at pixels 5 and 6")
write("i8_signed.fits", 8, np.arange(256), (16, 16), bzero=-128,
      comment="BZERO=-128 signed-byte convention, physical -128..127")
write("u8_bscale.fits", 8, [0, 1, 100, 255], (4,), bscale=2.0, bzero=10.0)

# --- BITPIX 16 ------------------------------------------------------------
allv = np.arange(-32768, 32768, dtype=np.int64)
write("i16_all.fits", 16, allv, (256, 256), bscale=1.0, bzero=0.0,
      comment="every int16 value, BZERO=0 present")
write("u16_all.fits", 16, allv, (256, 256), bscale=1.0, bzero=32768.0,
      comment="every uint16 value 0..65535 (BZERO=32768)")
small16 = [-32768, -32767, -1000, -1, 0, 1, 1000, 32766, 32767]
write("i16_nobzero.fits", 16, small16, (9,), comment="no BSCALE/BZERO cards")
write("i16_bscale_neg.fits", 16, small16, (9,), bscale=-1.0, bzero=0.0)
write("i16_bscale_frac.fits", 16, small16, (9,), bscale=0.5, bzero=-3.0)
write("u16_blank.fits", 16, [-32768, -32767, -1, 0, 1, 32767, -32768, 12345], (8,),
      bzero=32768.0, blank=-32768, comment="BLANK=-32768 (physical 0) at pixels 0 and 6")
write("i16_bscale0.fits", 16, [1, 2, 3], (3,), bscale=0.0, comment="BSCALE=0 is invalid")

# --- BITPIX 32 ------------------------------------------------------------
i32 = np.concatenate([[-2**31, -2**31 + 1, -65536, -1, 0, 1, 255, 65535, 65536, 2**31 - 2, 2**31 - 1],
                      rng.integers(-2**31, 2**31, size=53)])
write("i32.fits", 32, i32, (8, 8))
write("u32.fits", 32, i32, (8, 8), bzero=2147483648.0, comment="BZERO=2^31 unsigned 32-bit")
i32b = i32.copy(); i32b[[3, 40]] = -999
write("i32_blank_bscale.fits", 32, i32b, (8, 8), bscale=2.5, bzero=1.0e6, blank=-999,
      comment="BLANK=-999 at pixels 3 and 40")
write("i32_bscale_neg.fits", 32, i32, (8, 8), bscale=-3.0, bzero=7.0)

# --- BITPIX 64 ------------------------------------------------------------
i64 = np.concatenate([np.array([-2**63, -2**63 + 1, -2**53 - 1, -1, 0, 1, 2**53 + 1, 2**62 + 1,
                                2**63 - 2, 2**63 - 1], dtype=np.int64),
                      rng.integers(-2**63, 2**63 - 1, size=54, dtype=np.int64)])
write("i64.fits", 64, i64, (8, 8))
write("u64.fits", 64, i64, (8, 8), bzero=9.223372036854775808e18,
      comment="BZERO=2^63 unsigned 64-bit (real-valued card)")
i64b = i64.copy(); i64b[[2, 9]] = -7
write("i64_blank.fits", 64, i64b, (8, 8), blank=-7, comment="BLANK=-7 at pixels 2 and 9")

# --- float ----------------------------------------------------------------
f32u = np.array([0.0, 0.25, 0.5, 1.0, -0.0, np.nan, 1e-40, 0.75, 0.1, 0.9999999, 2.0**-24, np.nan],
                dtype=np.float32)
write("f32_unit.fits", -32, f32u, (12,), comment="values in [0,1] plus NaN and a denormal")
f32w = np.array([0.0, 1.0, 100.0, 65535.0, -5.0, np.inf, -np.inf, np.nan, 3.4e38], dtype=np.float32)
write("f32_wide.fits", -32, f32w, (9,), comment="values outside [0,1], Inf, NaN")
write("f32_bscale.fits", -32, np.array([0.0, 0.25, -1.0, np.nan], dtype=np.float32), (4,),
      bscale=2.0, bzero=1.0)
f64 = np.array([0.0, 0.5, 1.0, 1e-50, np.nan, 0.1, 1.0 / 3.0, 0.999999999, 2.0**-149, 2.0**-150,
                2.0**-150 * 3, -0.0], dtype=np.float64)
write("f64.fits", -64, f64, (12,), comment="underflow, NaN, rounding to float32, values in [0,1]")
write("f64_wide.fits", -64, np.array([0.0, 65535.0, -1e30, 3.4028234663852886e38, np.inf, np.nan]), (6,),
      comment="values outside [0,1] that fit float32")
write("f64_overflow.fits", -64, np.array([0.5, 1e300, 0.25]), (3,),
      comment="1e300 does not fit float32: ReadFloat32 must fail")
write("f32_overflow.fits", -32, np.array([0.5, 3e38, 0.25], dtype=np.float32), (3,), bscale=2.0,
      comment="physical 6e38 does not fit float32: ReadFloat32 must fail")

# --- frame helpers ----------------------------------------------------------
rgb = rng.integers(-32768, 32768, size=3 * 5 * 4)
write("u16_rgb.fits", 16, rgb, (5, 4, 3), bzero=32768.0, comment="planar RGB, NAXIS3=3")
orient = np.arange(4 * 3) * 1000 - 32768  # pixel (x, y) stores (y*4+x)*1000 physical
write("u16_orient.fits", 16, orient, (4, 3), bzero=32768.0,
      comment="pixel (x,y) holds physical 1000*(4y+x)")

# --- 2D variants for AsImage / image.Decode --------------------------------
write("u8_bscale_2d.fits", 8, np.arange(12), (4, 3), bscale=2.0, bzero=10.0)
write("u16_blank_2d.fits", 16, [-32768, 0, 5, -32768, 7, 9], (3, 2), bzero=32768.0, blank=-32768,
      comment="BLANK at pixels 0 and 3")
write("u16_noblankhit_2d.fits", 16, [0, 1, 2, 3, 4, 5], (3, 2), bzero=32768.0, blank=-32768,
      comment="BLANK card present, no pixel equals it")
write("u16_bscale2_2d.fits", 16, [0, 1, 2, 3, 4, 5], (3, 2), bscale=2.0, bzero=32768.0)
write("f32_2d.fits", -32, np.array([0.0, 0.5, np.nan, -1.0, 1e30, np.inf], dtype=np.float32), (3, 2))
write("f32_bscale_2d.fits", -32, np.array([0.0, 0.5, 1.0, 2.0], dtype=np.float32), (2, 2), bscale=3.0)
write("f64_exact_2d.fits", -64, np.array([0.5, 0.25, -0.0, np.nan, 7.0, 2.0**-149, 3.0, 65535.0]), (4, 2),
      comment="every value is exactly representable in float32")
write("f64_2d.fits", -64, np.array([0.1, 1.0 / 3.0, 2.0**-150, np.nan, -0.0, 7.0]), (2, 3))

# --- added by the second proof pass (no rng use: earlier bytes unchanged) --
# BITPIX 64 values u = raw + 2^63 where float32(float64(u) / (2^64-1))
# (double rounding) differs from the float32 nearest u / (2^64-1); found by
# a numpy/fractions search near float32 rounding midpoints. The last one is
# 2^39+1 (raw), u/(2^64-1) just above 1/2.
U_ROUND = [8784523776, 181791301632, 3551877791744, 75104839008256, 1500332270026752,
           11940132026318848, 28795950759673855, 56535610147471360, 119150402905047041,
           254187340082184192, 493261824121110529, 996722511943041022, 2097847776112017407,
           4425969021208231934, 2**63 + 2**39 + 1]
write("i64_round.fits", 64, [u - 2**63 for u in U_ROUND] + [-2**63, 2**63 - 1], (17,),
      comment="u=raw+2^63 where double rounding of u/(2^64-1) is wrong")
write("i64_round_neg.fits", 64, [2**63 - 1 - u for u in U_ROUND] + [-2**63, 2**63 - 1], (17,),
      bscale=-1.0, bzero=0.0, comment="BSCALE=-1 mirror: (2^64-1)-u is each U_ROUND value")
write("f32_bscale_inexact_2d.fits", -32, np.array([1.0, 2.0, 3.0, 0.5], dtype=np.float32), (2, 2),
      bscale=0.1, comment="physical 0.1*raw is not exactly representable in float32")
write("f32_bscale0_2d.fits", -32, np.array([1.0, 2.0, 3.0, 0.5], dtype=np.float32), (2, 2),
      bscale=0.0, comment="BSCALE=0 is invalid")

# --- added by the third proof pass (no rng use: earlier bytes unchanged) ---
# The only two BITPIX 32 values u = raw + 2^31 where float32(float64(u) /
# (2^32-1)) differs from the float32 nearest u / (2^32-1): an exhaustive
# search over all 2^32 values (find_i32_round.py, exact Fraction check).
U32_ROUND = [4294966655, 4294967167]
write("i32_round.fits", 32, [u - 2**31 for u in U32_ROUND] + [-2**31, 2**31 - 1], (4,),
      comment="u=raw+2^31 where double rounding of u/(2^32-1) is wrong")
write("i32_round_neg.fits", 32, [2**31 - 1 - u for u in U32_ROUND] + [-2**31, 2**31 - 1], (4,),
      bscale=-1.0, bzero=0.0, comment="BSCALE=-1 mirror: (2^32-1)-u is each U32_ROUND value")

print("samples written")

# --- added by the fourth proof pass (no rng use: earlier bytes unchanged) --
# uint64 with BZERO = 2^63 written as the integer literal 9223372036854775808,
# the form cfitsio and astropy write (u64.fits uses the real-valued spelling).
write("u64_intcard.fits", 64, i64, (8, 8), bzero=2**63,
      comment="BZERO=2^63 unsigned 64-bit (integer-literal card)")

# --- added by the round-2 proof pass (no rng use: earlier bytes unchanged) --
# Float BITPIX with BZERO and no BSCALE card (BSCALE defaults to 1): the
# offset must still be applied. Every physical value is exact in float32.
write("f32_bzero_2d.fits", -32, np.array([0.5, 1.5, -2.25, np.nan, 100.0, -0.0], dtype=np.float32), (3, 2),
      bzero=10.0, comment="BZERO=10 without BSCALE: physical raw+10")
write("f64_bzero_2d.fits", -64, np.array([0.5, 2.0, np.nan, 1000.0]), (2, 2), bzero=-3.0,
      comment="BZERO=-3 without BSCALE: physical raw-3")
write("f32_bzero_inexact_2d.fits", -32, np.array([0.001, 0.5, 2.0, 1e-7], dtype=np.float32), (2, 2),
      bzero=10.0, comment="BZERO=10: physical raw+10 not exact in float32")
# BSCALE / BZERO cards that are present but not numbers (fitsverify
# rejects them; cfitsio 4.7.0 fits_read_img silently uses the default):
# ReadFloat32 and AsImage must refuse them, not return unscaled data.
write("i16_bscale_str_2d.fits", 16, [0, 1, 100, -5, 7, 9], (3, 2), bzero=32768.0,
      raw_cards=[("BSCALE", "'abc'")], comment="BSCALE is a string: invalid")
write("i16_bscale_bad_2d.fits", 16, [0, 1, 100, -5, 7, 9], (3, 2), bzero=32768.0,
      raw_cards=[("BSCALE", "2.0.0")], comment="BSCALE is not a number: invalid")
write("f32_bzero_huge_2d.fits", -32, np.array([0.5, 1.5, 2.5, 3.5], dtype=np.float32), (2, 2),
      raw_cards=[("BZERO", "1E400")], comment="BZERO overflows a double: invalid")
# ReadFrame HDU selection: a 1-D primary (NAXIS=1, so not a frame) in
# front of a 2-D IMAGE extension holding the u16_orient.fits pixels.
# ReadFrame must return HDU 1's image.
with open(os.path.join(HERE, "frame_naxis1_primary.fits"), "wb") as f:
    f.write(hdu_bytes(16, [10, 20, 30, 40, 50], (5,), bzero=32768.0, extend=True,
                      comment="1-D primary; the frame is in HDU 1"))
    f.write(hdu_bytes(16, orient, (4, 3), bzero=32768.0, xtension=True,
                      comment="pixel (x,y) holds physical 1000*(4y+x)"))

# --- added by the round-2 second pass (appended; no rng use) ---------------
# Scaling cards whose value is a quoted numeric string. The FITS value is
# a character string, not a number (FITS 4.0 4.2.1), so the cards are
# invalid; cfitsio 4.7.0 fits_read_img silently uses the default (BZERO 0,
# BSCALE 1, status 0). ReadFloat32 / AsImage must refuse them on a fresh
# Open and after the header is parsed alike.
write("f32_bzero_qnum_2d.fits", -32, np.array([0.5, 1.5, -2.25, np.nan, 100.0, -0.0], dtype=np.float32),
      (3, 2), raw_cards=[("BZERO", "'10.0'")], comment="BZERO is the string '10.0': invalid")
write("i16_bscale_qnum_2d.fits", 16, [0, 1, 100, -5, 7, 9], (3, 2), bzero=32768.0,
      raw_cards=[("BSCALE", "'0.5'")], comment="BSCALE is the string '0.5': invalid")
# A C99 hex float: strtod (so cfitsio) reads 8, but it is not a FITS
# integer or real (4.2.3 / 4.2.4), so it is invalid.
write("f32_bzero_hex_2d.fits", -32, np.array([0.5, 1.5, -2.25, np.nan, 100.0, -0.0], dtype=np.float32),
      (3, 2), raw_cards=[("BZERO", "0x1p3")], comment="BZERO = 0x1p3 is not a FITS number")
# D exponents are valid FITS reals (4.2.4): cfitsio and astropy read
# BSCALE = 5.0D-1 and BZERO = -3.0D0 as 0.5 and -3, and BZERO = 1.0D1 as 10.
write("i16_dexp_2d.fits", 16, small16[:6], (3, 2),
      raw_cards=[("BSCALE", "5.0D-1"), ("BZERO", "-3.0D0")], comment="BSCALE/BZERO with D exponents")
write("f32_bzero_dexp_2d.fits", -32, np.array([0.5, 1.5, -2.25, np.nan, 100.0, -0.0], dtype=np.float32),
      (3, 2), raw_cards=[("BZERO", "1.0D1")], comment="BZERO = 1.0D1 (D exponent): physical raw+10")

# --- added by round-2 proof2 (appended; no rng use) -------------------------
# card_*: BSCALE / BZERO / BLANK card variants on which cfitsio 4.7.0
# fits_read_img and astropy 8.0.2 agree or disagree (PROVENANCE.md has the
# reading of each). The rule they prove: the LAST card of each keyword is
# the one applied (cfitsio and astropy both use it), and it must be a valid
# value card (keyword upper case in columns 1-8, "= " in columns 9-10, a
# FITS number; an integer for BLANK) or the read is an error. cfnull
# goldens hold cfitsio's reading of every card_* file.
F32 = np.array([0.5, 1.5, -2.25, np.nan, 100.0, -0.0], dtype=np.float32)
I16 = [0, 1, 100, -5, 7, 9]
for name, cards in [
    ("dup", ["BZERO   =                 10.0", "BZERO   =                 20.0"]),
    ("dupbadfirst", ["BZERO   =                'abc'", "BZERO   =                 20.0"]),
    ("dupbadlast", ["BZERO   =                 20.0", "BZERO   =                'abc'"]),
    ("noeqsp", ["BZERO   =10.0"]),
    ("col10eq", ["BZERO    =                10.0"]),
    ("col11eq", ["BZERO     =  10.0"]),
    ("spacein", ["BZERO   =               1 0.0"]),
    ("trail", ["BZERO   =               10.0 x"]),
    ("lowerd", ["BZERO   =                1.0d1"]),
    ("lowere", ["BZERO   =                1.0e1"]),
    ("lowkey", ["bzero   =                 10.0"]),
    ("empty", ["BZERO   =                      / no value"]),
    ("int", ["BZERO   =                   10"]),
    ("plus", ["BZERO   =               +1.0E1"]),
    ("dotlead", ["BZERO   =                 .5E1"]),
    ("dottrail", ["BZERO   =                  10."]),
]:
    write(f"card_f32_bzero_{name}_2d.fits", -32, F32, (3, 2), verbatim=cards,
          comment=f"BZERO card variant: {name}")
write("card_i16_bscale_dup_2d.fits", 16, I16, (3, 2),
      verbatim=["BSCALE  =                  2.0", "BZERO   =              32768.0", "BSCALE  =                 -1.0"],
      comment="BSCALE twice: 2.0 then -1.0 (the last is applied)")
for name, cards in [
    ("dup", ["BLANK   =                    1", "BLANK   =                  100"]),
    ("float", ["BLANK   =                100.0"]),
    ("str", ["BLANK   =                'abc'"]),
    ("noeq", ["BLANK   =100"]),
]:
    write(f"card_i16_blank_{name}_2d.fits", 16, I16, (3, 2), bscale=2.0, verbatim=cards,
          comment=f"BLANK card variant: {name}")
# "BZERO=10.0" ('=' in column 6): cfitsio reads 10, astropy 0 -> error.
write("card_f32_bzero_eqcol6_2d.fits", -32, F32, (3, 2), verbatim=["BZERO=10.0"],
      comment="BZERO card variant: eqcol6")
# A BZERO card after END (in the header block's padding): neither cfitsio
# nor astropy reads it (both +10).
b = bytearray(hdu_bytes(-32, F32, (3, 2), verbatim=["BZERO   =                 10.0"],
                        comment="BZERO card variant: afterend"))
end = next(i for i in range(0, 2880, 80) if b[i:i + 8] == b"END     ")
b[end + 80:end + 160] = "BZERO   =                 20.0".ljust(80).encode("ascii")
with open(os.path.join(HERE, "card_f32_bzero_afterend_2d.fits"), "wb") as f:
    f.write(bytes(b))
# BLANK on float BITPIX is undefined (FITS 4.0 4.4.2.4) and not read, even
# when its value is invalid: cfitsio and astropy both ignore it.
write("card_f32_blank_str_2d.fits", -32, F32, (3, 2), verbatim=["BLANK   =                'abc'"],
      comment="BLANK on float data: ignored")
