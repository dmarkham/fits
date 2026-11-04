#!/usr/bin/env python3
"""Write the hand-built header samples in this directory.

These samples are INPUT bytes, not truth: each is a primary HDU with
NAXIS = 0 whose header holds deliberately malformed or edge-case cards,
written byte for byte with plain Python (no FITS library, no Go code).
What the cards MEAN is defined by cfitsio reading the same bytes: the
goldens in golden/ are produced by testdata/cfref (cfref dump) and
hdrref (hdrref longstr), see Makefile and PROVENANCE.md.

Re-running this script reproduces the committed samples byte for byte.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def card(s):
    if len(s) > 80:
        raise SystemExit("card longer than 80: %r (%d)" % (s, len(s)))
    return s.ljust(80).encode("ascii")


def fits(cards):
    body = b"".join(card(c) for c in cards) + card("END")
    pad = (-len(body)) % 2880
    return body + b" " * pad


PRIMARY = [
    "SIMPLE  =                    T / file does conform to FITS standard",
    "BITPIX  =                    8 / bits per pixel",
    "NAXIS   =                    0 / no data",
]

# Value cards that do not follow the FITS standard. cfitsio opens the file
# and reads every other card; each of these is judged on its own.
MALFORMED = PRIMARY + [
    "DATE-OBS= 2020-01-01T00:00:00 / unquoted date string",
    "X       = TRUE / not a FITS logical",
    "TLONG   = Tx / starts with T",
    "BIGINT  = 123456789012345678901234 / 24-digit integer",
    "BZEROU  =  9223372036854775808 / 2**63, as cfitsio writes uint64 BZERO",
    "U64MAX  = 18446744073709551615 / 2**64-1",
    "NEGBIG  = -9223372036854775809 / below int64 min",
    "I64MIN  = -9223372036854775808 / int64 min",
    "UNTERM  = 'no closing quote",
    "FLTOVF  = 1.0E400 / overflows float64",
    "NANTOK  = NaN / not a FITS number",
    "INFTOK  = Inf / not a FITS number",
    "FUNDER  = 1.0E-400 / underflows to zero",
    "FDENORM = 4.9E-324 / subnormal",
    "FZERO   = 0.0E-400 / zero mantissa",
    "HEXV    = 0x10 / hexadecimal, C strtod reads it",
    "LOWD    = 1.5d3 / lowercase d exponent",
    "UPPD    = 1.5D3 / uppercase D exponent",
    "LOWE    = 1e5 / lowercase e exponent, no point",
    "TWOTOK  = 1 2 / two tokens",
    "CPXBAD  = (1.0, 2.0 / no closing paren",
    "CPXOK   = (1.5, -2.5) / complex",
    "CPXINT  = (3, -4) / integer complex",
    "NOSPACE =1 / value indicator without the space",
    "EQ10     = 5 / '=' in column 10",
    "SLASH   = / undefined value, comment only",
    "UNDEF   =",
    "CMTSP   =                    5 /   three spaces",
    "CMTNOSP =                    5 /no space",
    "NOSLASH =                    5 trailing junk without slash",
    "EMPTYSTR= ''",
    "SPACESTR= '   '",
    "LEADSP  = '  lead' / leading blanks are significant",
    "TRAILSP = 'trail   ' / trailing blanks are not",
    "QUOTES  = 'it''s' / doubled quote",
    "FREEFMT =      3.25 / free-format float",
    "FPOINT  = 120. / trailing point",
    "FNOINT  = .5 / leading point",
    "lowkey  =                    3 / lower-case keyword",
    "HIERARCH NO EQUALS SIGN here",
    "HIERARCH  ESO   DBL   SPACE = 7 / runs of blanks in the path",
    "HIERARCH= 5 / HIERARCH without the following blank",
    "FINAL   = 'still read' / proves later cards decode",
]

# CONTINUE edge cases.
FULL = "FULL    = 'abc&' / "
FULL = FULL + "c" * (80 - len(FULL))  # first card is exactly 80 non-blank wide
CONTINUE = PRIMARY + [
    "LONGSTR = 'first part &' / c1",
    "CONTINUE  'second part &' / c2",
    "CONTINUE  'third' / c3",
    "LONELY  = 'abc' / no ampersand",
    "CONTINUE  'def' / lone continue",
    "NUMPRE  =                    5",
    "CONTINUE  'after int'",
    "AMPEND  = 'dangling &' / & but no CONTINUE follows",
    "NEXT    =                    1",
    "AMPEMPTY= 'x&'",
    "CONTINUE  ''",
    "AMPCMT  = 'y&' / first",
    "CONTINUE  '' / only a comment",
    "CONTEQ  = 'z&'",
    "CONTINUE= 'value indicator: not a continuation'",
    "UNQ     = 'u&'",
    "CONTINUE  unquoted",
    "HIERARCH ESO OBS LONGNAME = 'hier part one &' / hc1",
    "CONTINUE  'hier part two' / hc2",
    FULL,
    "CONTINUE  'def' / more",
    "QSPLIT  = 'it''&'",
    "CONTINUE  '''s'",
    "FINAL   = 'end'",
]

# Further edge cases, most of them first written by an adversarial review
# of the decoder: keyword names cut by a blank or '=' inside columns 1-8
# or running past column 8 (cfitsio's ffgknm), number and complex forms
# outside FITS, and cards that reach cfitsio's buffer limits (a value
# token over 70 characters, a comment over 72) so the truncation is pinned.
ADVERSARIAL = PRIMARY + [
    "A1      =     'abc' 'def'",
    "A2      = 'a''",
    "A3      = 1.0E+3a",
    "A4      = +.5",
    "A5      = 1.",
    "A6      = .e5",
    "A7      = 1E",
    "A8      = - 5",
    "A9      = 00012",
    "B1      = -0",
    "B2      = (1,2) extra",
    "B3      = ( 1 , 2 )",
    "B4      = (1.5D2, 3)",
    "B5      = (1,2,3)",
    "B6      = 'str'/nocomment-space",
    "B7      = T/c",
    "B8      =T",
    "B9  =    5",
    "        = 5",
    "C1      = 1.5e3",
    "C2      = 1.5E3D",
    "C3      = 1D5",
    "C4      = 0x1p3",
    "C5      = 1_000",
    "C6      = ''''",
    "C7      = '   '",
    "C8      = 'a' / ",
    "C9      = 12345678901234567890",
    "D1      = -9223372036854775809",
    "D2      = +5",
    "D3      = 2.2250738585072011e-308",
    "D4      = 0.0E-400",
    "D5      = 1.7976931348623159E308",
    "D6      = 1.7976931348623158E308",
    "D7      = 1\t/ tab",
    "D8      = F00",
    "COMMENT = 'x'",
    "HISTORY= x",
    "key     = 5",
    "D9      = 'abc" + " " * 60,
    "E1      = 'a&'",
    "CONTINUE  'b&'",
    "CONTINUE  '' / only comment",
    "CONTINUE  'c'",
    "E2      = 1E+",
    "E3      = 1.5E+3.0",
    "E4      = ++5",
    "E5      = 5-",
    "E6      = (1, )",
    "E7      = (,1)",
    "E8      = ( )",
    "E9      = (1 2)",
    "F1      = 'x' junk / c",
    "F2      =    / just slash",
    "F3      = //double",
    "F4      = '' ",
    "F5      = 1.0E-320",
    "F6      = 4.9E-324",
    "F7      = 0.000000000000000000000000000000E5",
    "F8      = 1e5",
    "F9      = 1.0d0",
    "G1      = 1.0E5 / c",
    "G2      = TRUE",
    "G3      = F",
    "G4      = t",
    "G5      = 'a' 'b'",
    "XTENSION= 'IMAGE   '",
    "G6      =                     1234567890123456789012",
    # keyword names (cfitsio's ffgknm ends the name at the first blank or '=')
    "K1 K2   = 5 / blank inside columns 1-8",
    "KEYLONG10= 5 / no blank in columns 1-9",
    "NOEQUALSIGN" + "x" * 69,
    "K3=4 / '=' in column 3",
    # cfitsio buffer limits: value token > 70, quoted string > 70,
    # complex > 70, comment > 72
    "L1=" + "1" * 77,
    "L2=" + "x" * 77,
    "L3='" + "q" * 75 + "'",
    "L4=(1" + "0" * 71 + ",2)",
    "L5=1/" + "c" * 75,
    "L6=1 / " + "d" * 73,
    "FINAL   = 'end'",
]

SAMPLES = {
    "malformed.fits": MALFORMED,
    "continue.fits": CONTINUE,
    "adversarial.fits": ADVERSARIAL,
}


def main():
    for name, cards in SAMPLES.items():
        with open(os.path.join(HERE, name), "wb") as f:
            f.write(fits(cards))
        print("wrote", name)


if __name__ == "__main__":
    sys.exit(main())
