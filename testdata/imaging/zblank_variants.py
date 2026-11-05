#!/usr/bin/env python3
"""Derive the ZBLANK samples from fpack output (stdlib only; run by
make_compressed.sh). fpack copies an integer image's BLANK card into the
compressed header; these variants use the other null conventions of the
tiled-image compression convention, as cfitsio reads them
(imcomp_get_compressed_image_par: a ZBLANK column, else the ZBLANK
keyword, else BLANK):

  comp_u16_zblank.fits.fz         BLANK card renamed ZBLANK
  comp_i32_zblank_bscale.fits.fz  same, BITPIX 32 with BSCALE/BZERO
  comp_u16_zblank_blank.fits.fz   ZBLANK = -32768 and BLANK = -32767
                                  (cfitsio uses ZBLANK; the stored -32767 at
                                  pixel 1 is a valid pixel)
  comp_u16_zblank_col.fits.fz     a per-tile ZBLANK column (1J, -32768)
                                  instead of the keyword
  comp_u8_zbscale.fits.fz         comp_u8_blank plus ZBSCALE = 2 and
                                  ZBZERO = 5 cards, which are not part of
                                  the convention and which cfitsio ignores

and two scaling-card variants of comp_i16_bscale_neg (BSCALE -1.0) for the
scaling-card rule on a tile-compressed header:

  comp_card_i16_bscale_str.fits.fz      BSCALE  = 'abc' (not a number)
  comp_card_i16_bscale_hierarch.fits.fz BSCALE card replaced by
                                        HIERARCH BSCALE = -1.0 (cfitsio
                                        applies it, astropy ignores it)

Only header cards and the main table are rewritten; the compressed tile
bytes (heap) are copied verbatim. Deterministic.
"""
import os
import struct

HERE = os.path.dirname(os.path.abspath(__file__))
B = 2880


def blocks(n):
    return (n + B - 1) // B * B


def split(path):
    d = open(os.path.join(HERE, path), "rb").read()
    # primary header + data (NAXIS=0: header only)
    i = 0
    while True:
        card = d[i:i + 80]
        i += 80
        if card.startswith(b"END     "):
            break
    prim = d[:blocks(i)]
    rest = d[len(prim):]
    cards = []
    j = 0
    while True:
        c = rest[j:j + 80].decode("ascii")
        j += 80
        if c.startswith("END     "):
            break
        cards.append(c)
    data = rest[blocks(j):]
    return prim, cards, data


def card(key, value, comment=""):
    if isinstance(value, str):
        s = f"{key:<8}= '{value:<8}'"
    else:
        s = f"{key:<8}= {value:>20}"
    if comment:
        s += " / " + comment
    return s.ljust(80)[:80]


def get(cards, key):
    for c in cards:
        if c[:8].rstrip() == key:
            return c[10:].split("/")[0].strip().strip("'").strip()
    raise KeyError(key)


def setcard(cards, key, value, comment=""):
    for n, c in enumerate(cards):
        if c[:8].rstrip() == key:
            cards[n] = card(key, value, comment)
            return
    raise KeyError(key)


def write(name, prim, cards, data):
    hdr = "".join(cards + ["END".ljust(80)]).encode("ascii")
    hdr += b" " * (blocks(len(hdr)) - len(hdr))
    with open(os.path.join(HERE, name), "wb") as f:
        f.write(prim + hdr + data)


def rename_blank(src, dst, extra_blank=None):
    prim, cards, data = split(src)
    n = [k for k, c in enumerate(cards) if c[:8] == "BLANK   "]
    assert len(n) == 1, src
    v = int(get(cards, "BLANK"))
    cards[n[0]] = card("ZBLANK", v, "null value of the uncompressed integer image")
    if extra_blank is not None:
        cards.insert(n[0] + 1, card("BLANK", extra_blank, "not the null value: ZBLANK wins (cfitsio)"))
    write(dst, prim, cards, data)


def zblank_column(src, dst):
    prim, cards, data = split(src)
    assert get(cards, "TFIELDS") == "1"
    width, rows = int(get(cards, "NAXIS1")), int(get(cards, "NAXIS2"))
    pcount = int(get(cards, "PCOUNT"))
    blank = int(get(cards, "BLANK"))
    main, heap = data[:width * rows], data[width * rows:width * rows + pcount]
    newmain = b"".join(main[r * width:(r + 1) * width] + struct.pack(">i", blank) for r in range(rows))
    setcard(cards, "NAXIS1", width + 4, "width of table in bytes")
    setcard(cards, "TFIELDS", 2, "number of fields in each row")
    k = [i for i, c in enumerate(cards) if c[:8] == "TFORM1  "][0]
    cards.insert(k + 1, card("TTYPE2", "ZBLANK", "null value of each tile"))
    cards.insert(k + 2, card("TFORM2", "1J", "data format of field: 4-byte INTEGER"))
    cards = [c for c in cards if c[:8] != "BLANK   "]
    body = newmain + heap
    write(dst, prim, cards, body + b"\0" * (blocks(len(body)) - len(body)))


def add_zbscale(src, dst):
    prim, cards, data = split(src)
    k = [i for i, c in enumerate(cards) if c[:8] == "ZNAXIS1 "][0]
    cards.insert(k + 1, card("ZBSCALE", 2.0, "non-standard: cfitsio ignores it"))
    cards.insert(k + 2, card("ZBZERO", 5.0, "non-standard: cfitsio ignores it"))
    write(dst, prim, cards, data)


def replace_card(src, dst, key, text):
    prim, cards, data = split(src)
    n = [k for k, c in enumerate(cards) if c[:8] == f"{key:<8}"]
    assert len(n) == 1, (src, key)
    cards[n[0]] = text.ljust(80)[:80]
    write(dst, prim, cards, data)


rename_blank("comp_u16_blank.fits.fz", "comp_u16_zblank.fits.fz")
add_zbscale("comp_u8_blank.fits.fz", "comp_u8_zbscale.fits.fz")
rename_blank("comp_i32_blank_bscale.fits.fz", "comp_i32_zblank_bscale.fits.fz")
rename_blank("comp_u16_blank.fits.fz", "comp_u16_zblank_blank.fits.fz", extra_blank=-32767)
zblank_column("comp_u16_blank.fits.fz", "comp_u16_zblank_col.fits.fz")
replace_card("comp_i16_bscale_neg.fits.fz", "comp_card_i16_bscale_str.fits.fz", "BSCALE", "BSCALE  = 'abc'")
replace_card("comp_i16_bscale_neg.fits.fz", "comp_card_i16_bscale_hierarch.fits.fz", "BSCALE", "HIERARCH BSCALE = -1.0")
print("zblank variants written")
