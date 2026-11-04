#!/usr/bin/env python3
"""Write headers/kpno_kp4m_mosaic1_ccd1.hdr from the NOIRLab archive file.

Source: NOIRLab Astro Data Archive, public (no login), file md5
42f68c5bef0524074d565f93b687c873 (KPNO 4 m Mosaic-1 InstCal, NGC4594,
DATE-OBS 2006-02-26, FILENAME K4M06A_20060223-R-kp4m20060226T094108):

    curl -o kpno.fits.fz \\
      https://astroarchive.noirlab.edu/api/retrieve/42f68c5bef0524074d565f93b687c873/
    python3 extract_kpno.py kpno.fits.fz      # writes the .hdr
    python3 extract_kpno.py --check kpno.fits.fz

The .fits.fz (41 MB) is not committed; only the header cards of the
'ccd1' extension are. Steps (stdlib + cfitsio funpack only):

  1. check the md5 of the download;
  2. funpack -C -S (cfitsio) restores the uncompressed image header;
  3. take the extension with EXTNAME = 'ccd1' card for card, in order,
     except: XTENSION -> SIMPLE = T (a standalone header file), and
     PCOUNT, GCOUNT, CHECKSUM, DATASUM dropped (extension-only or
     rewritten by funpack). Trailing blanks of each line are stripped.

Every WCS card (CTYPE/CRVAL/CRPIX/CD/EQUINOX/WATn_nnn) is verbatim, in the
original order: the WAT1_006/007 and WAT2_006/007 cards come after
WAT2_005 and other keywords (interleaved WAT cards).
"""
import hashlib
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MD5 = "42f68c5bef0524074d565f93b687c873"
OUT = os.path.join(HERE, "headers", "kpno_kp4m_mosaic1_ccd1.hdr")
DROP = {"PCOUNT", "GCOUNT", "CHECKSUM", "DATASUM"}


def hdus(blob):
    """Yield the header cards of each HDU of an uncompressed FITS stream."""
    pos = 0
    while pos < len(blob):
        cards = []
        while True:
            block = blob[pos:pos + 2880]
            if len(block) < 2880:
                sys.exit("truncated FITS stream")
            pos += 2880
            done = False
            for i in range(0, 2880, 80):
                c = block[i:i + 80].decode("ascii")
                cards.append(c)
                if c[:8] == "END     ":
                    done = True
                    break
            if done:
                break
        kv = {c[:8].rstrip(): c[10:].split("/")[0].strip() for c in cards}
        naxis = int(kv.get("NAXIS", "0"))
        n = 1 if naxis else 0
        for i in range(1, naxis + 1):
            n *= int(kv["NAXIS%d" % i])
        nbytes = abs(int(kv["BITPIX"])) // 8 * n
        pos += nbytes + (-nbytes) % 2880
        yield cards


def main():
    args = [a for a in sys.argv[1:] if a != "--check"]
    check = len(args) != len(sys.argv) - 1
    if len(args) != 1:
        sys.exit(__doc__)
    with open(args[0], "rb") as f:
        if hashlib.md5(f.read()).hexdigest() != MD5:
            sys.exit("%s: md5 is not %s" % (args[0], MD5))
    blob = subprocess.run(["funpack", "-C", "-S", args[0]], check=True,
                          capture_output=True).stdout
    for cards in hdus(blob):
        if any(c.startswith("EXTNAME = 'ccd1") for c in cards):
            break
    else:
        sys.exit("no ccd1 extension")
    lines = []
    for c in cards:
        key = c[:8].rstrip()
        if key in DROP:
            continue
        if key == "XTENSION":
            c = "SIMPLE  =                    T / file does conform to FITS standard"
        lines.append(c.rstrip())
    text = "\n".join(lines) + "\n"
    if check:
        if open(OUT).read() != text:
            sys.exit("%s differs" % OUT)
        return
    with open(OUT, "w") as f:
        f.write(text)


if __name__ == "__main__":
    main()
