/*
 * gen_tables -- writes the synthetic table samples in testdata/tables with
 * the system libcfitsio (4.7.0 for the committed files; see PROVENANCE.md).
 *
 *   gen_tables <outdir>
 *
 * Every file is written by cfitsio itself: cfitsio picks the byte layout,
 * the heap, the TZERO keywords for its unsigned/signed-byte TFORM codes
 * (U, V, W, S) and the inverse scaling of physical values. The Go tests
 * then compare the Go reader against cfitsio's own reading of these files
 * (testdata/cfref "dump" goldens), so nothing here depends on the Go code.
 *
 * Files:
 *   cf_bintable_fixed.fits  every fixed-width TFORM (L X B I J K A E D C M),
 *                           scalar and vector, TDIM, TNULL, TSCAL/TZERO,
 *                           the unsigned TZERO convention for B(S)/I(U)/J(V)
 *   cf_bintable_uint64.fits K with TZERO = 2^63 (cfitsio TFORM 1W)
 *   cf_vla_p.fits           1P VLAs of every element type
 *                           (L X B I J K A E D C M), TNULL, TSCAL/TZERO
 *   cf_vla_q.fits           the same with 1Q descriptors
 *   cf_ascii.fits           ASCII table: Iw, Fw.d, Ew.d, Dw.d, Aw, TNULL
 *                           strings, blank fields, D exponents, TSCAL/TZERO
 *   cf_ascii_implied.fits   the same columns with implied decimal points
 *                           and blanks embedded in the digits
 *   cf_bintable_kzero.fits  K columns with integer TZERO other than 0 and
 *                           2^63: physical values above 2^63-1 and below
 *                           -2^63, and 2^62+1001 (not a double)
 *   cf_bintable_substr.fits rAw substring array (20A5) next to a plain 12A
 *   cf_ascii_scaled.fits    ASCII TSCAL/TZERO on F, E and D fields, -0,
 *                           TNULL '-99' against a right-justified -99,
 *                           exponents with leading zeros
 *   cf_bintable_ranges.fits physical values on and one step past the
 *                           range of every Go integer type (exact K,
 *                           double D, float-scaled K), fractions, Inf,
 *                           NaN, VLAs
 *   cf_bintable_tzero_big.fits integer TZERO literals no double holds
 *                           (2^63-1, 2^62+1, -(2^63-1)) on K and J, and a
 *                           float literal 4611686018427387905.0
 *   cf_bintable_a3dtable.fits a binary table whose XTENSION is the legacy
 *                           'A3DTABLE' (cfitsio and astropy read it as a
 *                           BINTABLE; cfitsio will not write it, so the
 *                           card is byte-patched after the write)
 *   cf_bintable_3dtable.fits the same table with XTENSION '3DTABLE', the
 *                           other legacy name cfitsio reads as BINARY_TBL
 *   cf_bintable_badlogical.fits 1L and 3L columns with one byte each
 *                           patched to 0x01 (an invalid logical byte)
 *   cf_ascii_exact.fits     ASCII integers beyond int64/uint64 and at
 *                           their edges, "N.000" fields, I20 with TZERO
 *                           1000, digitless fields ('.', 'E5', '+')
 *   cf_bintable_scaled_exact.fits fractional TZERO, TSCAL 2 / 0.1 and
 *                           integer TZERO beyond int64 (2^63 on J,
 *                           2^63+1, 2^64): the exact physical value is
 *                           not the float64 arithmetic result
 *   cf_bintable_strarray.fits character columns with a multi-axis TDIM
 *                           (arrays of strings), and TDIM on 8A2 columns
 *   cf_ascii_scaled_exact.fits the ASCII counterparts, fractions under an
 *                           integer TZERO, fields beyond the double range,
 *                           a TNULL longer than its field
 *   cf_ascii_scaled_edge.fits negative zeros under integer / fractional /
 *                           big TZERO and TSCAL 3; integers beyond uint64
 *                           and fractions under TZERO beyond 2^64
 *   cf_bintable_scaled_special.fits D-exponent TSCAL/TZERO literals and
 *                           +-Inf / NaN in scaled E and D columns
 *   cf_ascii_tzero_round.fits / cf_bintable_tzero_round.fits values under
 *                           an integer TZERO whose correctly rounded
 *                           float64 / float32 is not what double
 *                           arithmetic gives
 *   refused/cf_bintable_cplx_scaled.fits TSCAL/TZERO on C and M columns,
 *                           which the library refuses (outdir/refused must
 *                           exist)
 *
 * Re-running the generator must reproduce the committed bytes exactly
 * (make verify). cfitsio writes a DATE card only when asked, and this
 * program never asks, so the output is deterministic.
 */
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

static void check(int st, const char *what) {
  if (st) {
    fprintf(stderr, "gen_tables: %s: ", what);
    fits_report_error(stderr, st);
    exit(1);
  }
}

static fitsfile *create(const char *dir, const char *name) {
  char path[4096];
  int st = 0;
  fitsfile *f;
  snprintf(path, sizeof path, "!%s/%s", dir, name);
  fits_create_file(&f, path, &st);
  check(st, path);
  fits_create_img(f, 8, 0, NULL, &st);
  check(st, "primary");
  return f;
}

static void closef(fitsfile *f) {
  int st = 0;
  fits_close_file(f, &st);
  check(st, "close");
}

/* ---------------------------------------------------------------- fixed */
#define NR 6
static void fixed(const char *dir) {
  fitsfile *f = create(dir, "cf_bintable_fixed.fits");
  int st = 0;
  char *ttype[] = {"L1",  "L3",  "X1",   "X13",  "B1",   "B4",  "I1",  "I2",
                   "J1",  "K1",  "A12",  "E1",   "E3",   "D1",  "D6",  "C1",
                   "C2",  "M1",  "M2",   "S1",   "U1",   "V1",  "JSC", "ISC",
                   "ESC", "BNUL", "INUL", "JNUL", "KNUL"};
  char *tform[] = {"1L",  "3L",  "1X",  "13X", "1B",  "4B",  "1I",  "2I",
                   "1J",  "1K",  "12A", "1E",  "3E",  "1D",  "6D",  "1C",
                   "2C",  "1M",  "2M",  "1S",  "1U",  "1V",  "1J",  "1I",
                   "1E",  "1B",  "1I",  "1J",  "1K"};
  enum { NC = 29 };
  fits_create_tbl(f, BINARY_TBL, NR, NC, ttype, tform, NULL, "FIXED", &st);
  check(st, "create fixed");

  /* TDIM on the vector columns (B4 -> 2x2, D6 -> 3x2, E3 -> 3). */
  long d22[] = {2, 2}, d32[] = {3, 2}, d3[] = {3};
  fits_write_tdim(f, 6, 2, d22, &st);
  fits_write_tdim(f, 15, 2, d32, &st);
  fits_write_tdim(f, 13, 1, d3, &st);
  check(st, "tdim");

  /* Scaled columns: JSC TSCAL=0.5 TZERO=10, ISC TSCAL=1.5 TZERO=-3 with
   * TNULL=-32768, ESC (float) TSCAL=2 TZERO=1. */
  double v;
  v = 0.5; fits_write_key(f, TDOUBLE, "TSCAL23", &v, NULL, &st);
  v = 10;  fits_write_key(f, TDOUBLE, "TZERO23", &v, NULL, &st);
  v = 1.5; fits_write_key(f, TDOUBLE, "TSCAL24", &v, NULL, &st);
  v = -3;  fits_write_key(f, TDOUBLE, "TZERO24", &v, NULL, &st);
  long ln = -32768; fits_write_key(f, TLONG, "TNULL24", &ln, NULL, &st);
  v = 2;   fits_write_key(f, TDOUBLE, "TSCAL25", &v, NULL, &st);
  v = 1;   fits_write_key(f, TDOUBLE, "TZERO25", &v, NULL, &st);
  /* Integer TNULLs. */
  ln = 255;          fits_write_key(f, TLONG, "TNULL26", &ln, NULL, &st);
  ln = -32768;       fits_write_key(f, TLONG, "TNULL27", &ln, NULL, &st);
  ln = 2147483647;   fits_write_key(f, TLONG, "TNULL28", &ln, NULL, &st);
  LONGLONG kn = LLONG_MIN; fits_write_key(f, TLONGLONG, "TNULL29", &kn, NULL, &st);
  check(st, "scale keys");
  fits_set_hdustruc(f, &st); /* re-read column parameters */
  check(st, "hdustruc");

  /* L1: T F null T F null ; L3: patterns with nulls */
  char l1[NR] = {1, 0, 2, 1, 0, 2};
  fits_write_colnull(f, TLOGICAL, 1, 1, 1, NR, l1, &(char){2}, &st);
  char l3[NR * 3] = {1, 0, 1, 0, 0, 0, 1, 1, 1, 2, 1, 0, 0, 2, 2, 1, 0, 2};
  fits_write_colnull(f, TLOGICAL, 2, 1, 1, NR * 3, l3, &(char){2}, &st);
  check(st, "L");

  /* X1, X13 bits */
  char x1[NR] = {1, 0, 1, 1, 0, 0};
  /* cfitsio's TBIT writer does not wrap into the next row, so write the
   * bits of each row separately. */
  char x13[NR * 13];
  for (int i = 0; i < NR * 13; i++) x13[i] = ((i * 7 + i / 5) % 3) == 0;
  for (int r = 0; r < NR; r++) {
    fits_write_col(f, TBIT, 3, r + 1, 1, 1, &x1[r], &st);
    fits_write_col(f, TBIT, 4, r + 1, 1, 13, &x13[r * 13], &st);
  }
  check(st, "X");

  unsigned char b1[NR] = {0, 1, 127, 128, 254, 255};
  fits_write_col(f, TBYTE, 5, 1, 1, NR, b1, &st);
  unsigned char b4[NR * 4];
  for (int i = 0; i < NR * 4; i++) b4[i] = (unsigned char)(i * 11);
  fits_write_col(f, TBYTE, 6, 1, 1, NR * 4, b4, &st);

  short i1[NR] = {-32768, -1, 0, 1, 12345, 32767};
  fits_write_col(f, TSHORT, 7, 1, 1, NR, i1, &st);
  short i2[NR * 2] = {1, -1, 2, -2, 300, -300, 32767, -32768, 0, 0, 7, 8};
  fits_write_col(f, TSHORT, 8, 1, 1, NR * 2, i2, &st);

  int j1[NR] = {INT32_MIN, -1, 0, 1, 123456789, INT32_MAX};
  fits_write_col(f, TINT, 9, 1, 1, NR, j1, &st);

  /* K: extremes, a value not representable in float64 (2^62+1) and real
   * Gaia DR3 source_id style 64-bit identifiers. */
  LONGLONG k1[NR] = {LLONG_MIN, -1, (1LL << 62) + 1, 4295806720896741504LL,
                     6917528997577384320LL, LLONG_MAX};
  fits_write_col(f, TLONGLONG, 10, 1, 1, NR, k1, &st);
  check(st, "B I J K");

  char *a12[NR] = {"hello", "", "trailing  ", "exactly12chr", "a b c", " lead"};
  fits_write_col(f, TSTRING, 11, 1, 1, NR, a12, &st);
  check(st, "A");

  float e1[NR] = {1.5f, -0.0f, 3.4028235e38f, 1.0e-45f, NAN, INFINITY};
  fits_write_col(f, TFLOAT, 12, 1, 1, NR, e1, &st);
  float e3[NR * 3];
  for (int i = 0; i < NR * 3; i++) e3[i] = (float)i / 3.0f - 2.0f;
  e3[4] = NAN;
  fits_write_col(f, TFLOAT, 13, 1, 1, NR * 3, e3, &st);
  double d1[NR] = {M_PI, -1e300, 5e-324, 0.1, NAN, -INFINITY};
  fits_write_col(f, TDOUBLE, 14, 1, 1, NR, d1, &st);
  double d6[NR * 6];
  for (int i = 0; i < NR * 6; i++) d6[i] = i * 1.25 - 7.0;
  fits_write_col(f, TDOUBLE, 15, 1, 1, NR * 6, d6, &st);
  check(st, "E D");

  /* No NaN in complex cells: the shared cfref harness prints complex values
   * without their null flag, so a NaN (null) complex would dump undefined
   * memory and the golden would not reproduce. */
  float c1[NR * 2] = {1, -2, 3, -4, 0.5f, 0.25f, -1e30f, 1e-30f, 0, 0, -7, 1};
  fits_write_col(f, TCOMPLEX, 16, 1, 1, NR, c1, &st);
  float c2[NR * 4];
  for (int i = 0; i < NR * 4; i++) c2[i] = (float)(i - 10) * 0.75f;
  fits_write_col(f, TCOMPLEX, 17, 1, 1, NR * 2, c2, &st);
  double m1[NR * 2] = {M_PI, -M_E, 1e300, -1e-300, 0.1, 0.2, -0.0, 0.0, 7, 8, 9, -10};
  fits_write_col(f, TDBLCOMPLEX, 18, 1, 1, NR, m1, &st);
  double m2[NR * 4];
  for (int i = 0; i < NR * 4; i++) m2[i] = i * 1.0 / 7.0;
  fits_write_col(f, TDBLCOMPLEX, 19, 1, 1, NR * 2, m2, &st);
  check(st, "C M");

  /* S = signed byte via B + TZERO=-128; U = uint16 via I + TZERO=32768;
   * V = uint32 via J + TZERO=2147483648. cfitsio writes the TZERO card. */
  signed char s1[NR] = {-128, -1, 0, 1, 100, 127};
  fits_write_col(f, TSBYTE, 20, 1, 1, NR, s1, &st);
  unsigned short u1[NR] = {0, 1, 32767, 32768, 40000, 65535};
  fits_write_col(f, TUSHORT, 21, 1, 1, NR, u1, &st);
  unsigned int v1[NR] = {0, 1, 2147483647u, 2147483648u, 3000000000u, 4294967295u};
  fits_write_col(f, TUINT, 22, 1, 1, NR, v1, &st);
  check(st, "S U V");

  /* Scaled physical values; cfitsio applies the inverse on write. */
  double jsc[NR] = {10, 10.5, -5, 1e6, 0, 11};
  fits_write_col(f, TDOUBLE, 23, 1, 1, NR, jsc, &st);
  double isc[NR] = {-3, 0, 1.5, 4.5, -1.5, -3};
  fits_write_col(f, TDOUBLE, 24, 1, 1, NR, isc, &st);
  fits_write_col_null(f, 24, 6, 1, 1, &st); /* row 6 null -> TNULL raw */
  double esc[NR] = {1, 3, -1, 2.5, 1e10, 0};
  fits_write_col(f, TDOUBLE, 25, 1, 1, NR, esc, &st);
  check(st, "scaled");

  /* TNULL columns: rows 2 and 5 null. */
  unsigned char bn[NR] = {1, 0, 2, 3, 0, 254};
  fits_write_col(f, TBYTE, 26, 1, 1, NR, bn, &st);
  short in[NR] = {-1, 0, 2, 3, 0, 32767};
  fits_write_col(f, TSHORT, 27, 1, 1, NR, in, &st);
  int jn[NR] = {-1, 0, 2, 3, 0, INT32_MIN};
  fits_write_col(f, TINT, 28, 1, 1, NR, jn, &st);
  LONGLONG knn[NR] = {-1, 0, 2, 3, 0, LLONG_MAX};
  fits_write_col(f, TLONGLONG, 29, 1, 1, NR, knn, &st);
  for (int c = 26; c <= 29; c++) {
    fits_write_col_null(f, c, 2, 1, 1, &st);
    fits_write_col_null(f, c, 5, 1, 1, &st);
  }
  check(st, "tnull cols");
  closef(f);
}

/* --------------------------------------------------------------- uint64 */
static void uint64file(const char *dir) {
  fitsfile *f = create(dir, "cf_bintable_uint64.fits");
  int st = 0;
  char *ttype[] = {"W1"}, *tform[] = {"1W"};
  fits_create_tbl(f, BINARY_TBL, 5, 1, ttype, tform, NULL, "UINT64", &st);
  ULONGLONG w[5] = {0, 1, 9223372036854775807ULL, 9223372036854775808ULL, 18446744073709551615ULL};
  fits_write_col(f, TULONGLONG, 1, 1, 1, 5, w, &st);
  check(st, "uint64");
  closef(f);
}

/* ---------------------------------------------------- K, integer TZERO */
/* Raw (stored) values are written with scaling switched off, so the file
 * holds exactly these integers; the physical value is raw + TZERO, which
 * for these rows leaves the int64 range (KZP above 2^63-1, KZN below
 * -2^63) or is an int64 that no double can hold (KZI). */
static void kzero(const char *dir) {
  fitsfile *f = create(dir, "cf_bintable_kzero.fits");
  int st = 0;
  char *ttype[] = {"KZP", "KZN", "KZI"}, *tform[] = {"1K", "1K", "1K"};
  fits_create_tbl(f, BINARY_TBL, NR, 3, ttype, tform, NULL, "KZERO", &st);
  LONGLONG z = 1000;
  fits_write_key(f, TLONGLONG, "TZERO1", &z, NULL, &st);
  z = -1000;
  fits_write_key(f, TLONGLONG, "TZERO2", &z, NULL, &st);
  z = 1000;
  fits_write_key(f, TLONGLONG, "TZERO3", &z, NULL, &st);
  check(st, "kzero keys");
  fits_set_hdustruc(f, &st);
  for (int c = 1; c <= 3; c++) fits_set_tscale(f, c, 1.0, 0.0, &st);
  LONGLONG kzp[NR] = {0, -5, 4611686018427387905LL, LLONG_MAX - 1000, LLONG_MAX - 999, LLONG_MAX};
  LONGLONG kzn[NR] = {0, 5, -4611686018427387905LL, LLONG_MIN + 1000, LLONG_MIN + 999, LLONG_MIN};
  LONGLONG kzi[NR] = {-4611686018427387907LL, 4611686018427387905LL, -1000, 0, 4611686018427387903LL, LLONG_MIN + 1000};
  fits_write_col(f, TLONGLONG, 1, 1, 1, NR, kzp, &st);
  fits_write_col(f, TLONGLONG, 2, 1, 1, NR, kzn, &st);
  fits_write_col(f, TLONGLONG, 3, 1, 1, NR, kzi, &st);
  check(st, "kzero");
  closef(f);
}

/* ------------------------------------------------------ substring array */
static void substr(const char *dir) {
  fitsfile *f = create(dir, "cf_bintable_substr.fits");
  int st = 0;
  char *ttype[] = {"SUB", "PLAIN"}, *tform[] = {"20A5", "12A"};
  fits_create_tbl(f, BINARY_TBL, 2, 2, ttype, tform, NULL, "SUBSTR", &st);
  check(st, "create substr");
  char *sub[] = {"one", "two", "three", "four", "a", "bb", "ccc", "dddd"};
  fits_write_col(f, TSTRING, 1, 1, 1, 8, sub, &st);
  char *plain[] = {"hello", "world"};
  fits_write_col(f, TSTRING, 2, 1, 1, 2, plain, &st);
  check(st, "substr");
  closef(f);
}

/* ------------------------------------------------------------------ VLA */
#define VR 5
static void vla(const char *dir, const char *name, char d) {
  fitsfile *f = create(dir, name);
  int st = 0;
  const char *codes = "LXBIJKAEDCM";
  enum { NC = 13 };
  char *ttype[NC], *tform[NC];
  char tt[NC][16], tf[NC][24];
  for (int c = 0; c < 11; c++) {
    snprintf(tt[c], 16, "V%c", codes[c]);
    snprintf(tf[c], 24, "1%c%c(%d)", d, codes[c], codes[c] == 'X' ? 20 : 12);
  }
  /* Extra: scaled VLA (1PJ with TSCAL/TZERO) and unsigned VLA (1PU). */
  snprintf(tt[11], 16, "VJSC"); snprintf(tf[11], 24, "1%cJ(4)", d);
  snprintf(tt[12], 16, "VU");   snprintf(tf[12], 24, "1%cU(4)", d);
  for (int c = 0; c < NC; c++) { ttype[c] = tt[c]; tform[c] = tf[c]; }
  fits_create_tbl(f, BINARY_TBL, VR, NC, ttype, tform, NULL, d == 'P' ? "VLA_P" : "VLA_Q", &st);
  check(st, "create vla");
  /* TNULL on the B, I, J, K VLAs. */
  long ln;
  ln = 77; fits_write_key(f, TLONG, "TNULL3", &ln, NULL, &st);
  ln = 77; fits_write_key(f, TLONG, "TNULL4", &ln, NULL, &st);
  ln = 77; fits_write_key(f, TLONG, "TNULL5", &ln, NULL, &st);
  ln = 77; fits_write_key(f, TLONG, "TNULL6", &ln, NULL, &st);
  double v = 0.25; fits_write_key(f, TDOUBLE, "TSCAL12", &v, NULL, &st);
  v = 100; fits_write_key(f, TDOUBLE, "TZERO12", &v, NULL, &st);
  check(st, "vla keys");
  fits_set_hdustruc(f, &st);
  check(st, "hdustruc");

  long lens[VR] = {3, 0, 1, 7, 12};
  for (int r = 1; r <= VR; r++) {
    long n = lens[r - 1];
    /* L with a null in row 4 */
    char lv[12];
    for (int k = 0; k < n; k++) lv[k] = (char)((k + r) % 2);
    if (r == 4) lv[2] = 2;
    if (n) fits_write_colnull(f, TLOGICAL, 1, r, 1, n, lv, &(char){2}, &st);
    else fits_write_descript(f, 1, r, 0, 0, &st);
    check(st, "VL");
    /* X: n+r bits (lengths 4,2,4,11,17 bits) */
    long nb = n + r;
    char xb[20];
    for (int k = 0; k < nb; k++) xb[k] = ((k * 5 + r) % 3) == 0;
    fits_write_col(f, TBIT, 2, r, 1, nb, xb, &st);
    check(st, "VX");
    unsigned char bv[12];
    short iv[12];
    int jv[12];
    LONGLONG kv[12];
    float ev[12];
    double dv[12];
    float cv[24];
    double mv[24];
    for (int k = 0; k < n; k++) {
      bv[k] = (unsigned char)(k * 20 + r);
      iv[k] = (short)(k * -1000 + r);
      jv[k] = k * 100000 - r;
      kv[k] = (LONGLONG)k * 1000000000000LL + r;
      ev[k] = (float)k * 0.5f - r;
      dv[k] = k / 3.0 + r;
      cv[2 * k] = (float)k; cv[2 * k + 1] = -(float)r;
      mv[2 * k] = k * 0.1; mv[2 * k + 1] = r * -0.01;
    }
    if (n > 2) { bv[1] = 77; iv[1] = 77; jv[1] = 77; kv[1] = 77; } /* -> null */
    if (n > 3) kv[3] = (1LL << 62) + 3; /* not representable in float64 */
    if (n) {
      fits_write_col(f, TBYTE, 3, r, 1, n, bv, &st);
      fits_write_col(f, TSHORT, 4, r, 1, n, iv, &st);
      fits_write_col(f, TINT, 5, r, 1, n, jv, &st);
      fits_write_col(f, TLONGLONG, 6, r, 1, n, kv, &st);
      fits_write_col(f, TFLOAT, 8, r, 1, n, ev, &st);
      fits_write_col(f, TDOUBLE, 9, r, 1, n, dv, &st);
      fits_write_col(f, TCOMPLEX, 10, r, 1, n, cv, &st);
      fits_write_col(f, TDBLCOMPLEX, 11, r, 1, n, mv, &st);
    } else {
      for (int c = 3; c <= 11; c++) if (c != 7) fits_write_descript(f, c, r, 0, 0, &st);
    }
    check(st, "VB..VM");
    /* A: strings of length n */
    char s[16];
    for (int k = 0; k < n; k++) s[k] = (char)('a' + (k + r) % 26);
    s[n] = 0;
    if (r == 3) strcpy(s, "x");
    char *sp = s;
    if (n || r == 3) fits_write_col(f, TSTRING, 7, r, 1, 1, &sp, &st);
    else fits_write_descript(f, 7, r, 0, 0, &st);
    check(st, "VA");
    /* scaled J VLA: physical values, cfitsio stores (v-100)/0.25 */
    double ps[4] = {100, 100.25, 99.5, 1000};
    long ns = r % 4 + 1;
    fits_write_col(f, TDOUBLE, 12, r, 1, ns, ps, &st);
    unsigned short us[4] = {0, 32768, 65535, (unsigned short)(r * 1000)};
    fits_write_col(f, TUSHORT, 13, r, 1, ns, us, &st);
    check(st, "VJSC VU");
  }
  closef(f);
}

/* Columns: NAME A8, IVAL I6, FIMP F8.2, FDOT F8.3, EVAL E11.3, DVAL D14.6,
 * ISC I5 (TSCAL 2.5, TZERO -1). TBCOLs are explicit so that the raw row
 * text written below lines up with them. */
#define ANC 7
static const long a_tbcol[ANC] = {1, 10, 17, 26, 35, 47, 62};
static const int a_width[ANC] = {8, 6, 8, 8, 11, 14, 5};

/* Field text per row. Numeric fields are right-justified in their window;
 * a leading '<' means "left-justify" (cfitsio matches TNULL as a prefix of
 * the untrimmed field, so a TNULL marker must start the field). */
static void ascii(const char *dir, const char *name, const char *rows[][ANC], int nrows) {
  fitsfile *f = create(dir, name);
  int st = 0;
  char *ttype[] = {"NAME", "IVAL", "FIMP", "FDOT", "EVAL", "DVAL", "ISC"};
  char *tform[] = {"A8", "I6", "F8.2", "F8.3", "E11.3", "D14.6", "I5"};
  fits_create_tbl(f, ASCII_TBL, 0, ANC, ttype, tform, NULL, "ASCII", &st);
  check(st, "create ascii");
  char key[16];
  for (int c = 0; c < ANC; c++) {
    snprintf(key, sizeof key, "TBCOL%d", c + 1);
    fits_update_key(f, TLONG, key, (void *)&a_tbcol[c], NULL, &st);
  }
  long naxis1 = 66;
  fits_update_key(f, TLONG, "NAXIS1", &naxis1, NULL, &st);
  fits_write_key(f, TSTRING, "TNULL1", "NONE", NULL, &st);
  fits_write_key(f, TSTRING, "TNULL2", "*", NULL, &st);
  fits_write_key(f, TSTRING, "TNULL3", "null", NULL, &st);
  fits_write_key(f, TSTRING, "TNULL5", "-", NULL, &st);
  double v = 2.5; fits_write_key(f, TDOUBLE, "TSCAL7", &v, NULL, &st);
  v = -1; fits_write_key(f, TDOUBLE, "TZERO7", &v, NULL, &st);
  check(st, "ascii keys");
  fits_set_hdustruc(f, &st);
  check(st, "hdustruc");
  for (int r = 0; r < nrows; r++) {
    char buf[67];
    memset(buf, ' ', 66);
    buf[66] = 0;
    for (int c = 0; c < ANC; c++) {
      const char *s = rows[r][c];
      int left = s[0] == '<';
      if (left) s++;
      int n = (int)strlen(s);
      if (n > a_width[c]) { fprintf(stderr, "%s: field too wide row %d col %d\n", name, r, c); exit(1); }
      int at = (int)a_tbcol[c] - 1 + (left ? 0 : a_width[c] - n);
      memcpy(buf + at, s, (size_t)n);
    }
    fits_insert_rows(f, r, 1, &st);
    fits_write_tblbytes(f, r + 1, 1, 66, (unsigned char *)buf, &st);
    check(st, "ascii row");
  }
  closef(f);
}

/* Standard-conforming fields (fitsverify clean). */
static const char *ascii_std[][ANC] = {
    {"<alpha", "123", "123.45", "1.250", "1.250E+02", "1.500000D+00", "4"},
    /* all-blank fields: cfitsio reads numeric blanks as 0, not null */
    {"", "", "", "", "", "", ""},
    /* TNULL markers at the start of the field */
    {"<NONE", "<*", "<null", "7.000", "<-", "0.0", "0"},
    /* explicit signs, exponent forms */
    {"<gamma", "+12", "-1.50", "+2.500", "-1.5E+3", "-2.5D-03", "1"},
    {"<delta", "0", "1.", ".5", "5.E0", "1.0000000E+300", "99"},
    /* a string that looks like a number; full-width fields */
    {"<12345678", "-99999", "-9999.99", "9999.999", "-1.234E-100", "-1.234567D+100", "-9999"},
};

/* Fields without a decimal point: readers apply the implied decimal of
 * TFORM Fw.d / Ew.d / Dw.d (cfitsio does; fitsverify reports these as
 * errors, which is why they live in a file of their own).
 *   F8.2 "12345" = 123.45   F8.3 "-500" = -0.5   E11.3 "12500E-2" = 0.125
 *   D14.6 "1234567D2" = 123.4567 */
static const char *ascii_implied[][ANC] = {
    {"<beta", "-456", "12345", "-500", "12500E-2", "1234567D2", "-3"},
    {"<one", "1", "1", "1", "1", "1", "1"},
    {"<neg", "-7", "-7", "-7000", "-7E3", "-7D-3", "-7"},
    /* blanks between sign, digits, point and exponent are ignored:
     *   I6 "1 2" = 12   F8.2 "1 2.5" = 12.5   F8.3 "- 5" = -0.005
     *   E11.3 "1.5 E 2" = 150   D14.6 "2.5 D 1" = 25   I5 "- 3" */
    {"<blank", "1 2", "1 2.5", "- 5", "1.5 E 2", "2.5 D 1", "- 3"},
};

/* ASCII TSCAL/TZERO on float fields, negative zero, TNULL padding,
 * exponents with leading zeros. Layout: FSC F10.3 (TSCAL 2.5, TZERO -1),
 * ESC E12.4 (TSCAL 0.5, TZERO 100), DSC D14.6 (TSCAL -2, TZERO 0.25),
 * FNZ F8.3, INL I6 (TNULL '-99'), EXP E14.3. */
#define SNC 6
static void ascii_scaled(const char *dir) {
  static const long tbcol[SNC] = {1, 12, 25, 40, 49, 56};
  static const int width[SNC] = {10, 12, 14, 8, 6, 14};
  static const char *rows[][SNC] = {
      {"4.000", "1.0000E+02", "1.500000D+00", "-0.000", "-99", "1.000E+0000002"},
      {"-0.400", "-2.0000E+00", "-1.250000D-01", "-0.0", "<-99", "2.500E-0000001"},
      {"", "", "", "0.000", "12", "-7.000E+00"},
      {"1234.567", "3.1416E-03", "1.000000D+10", "-0.001", "-12", "1.000E+300"},
      {"-0.000", "-0.0000E+00", "-0.000000D+00", "-0.", "0", "-0.000E+00"},
  };
  enum { NROW = 5, W = 69 };
  fitsfile *f = create(dir, "cf_ascii_scaled.fits");
  int st = 0;
  char *ttype[] = {"FSC", "ESC", "DSC", "FNZ", "INL", "EXP"};
  char *tform[] = {"F10.3", "E12.4", "D14.6", "F8.3", "I6", "E14.3"};
  fits_create_tbl(f, ASCII_TBL, 0, SNC, ttype, tform, NULL, "SCALED", &st);
  check(st, "create ascii scaled");
  char key[16];
  for (int c = 0; c < SNC; c++) {
    snprintf(key, sizeof key, "TBCOL%d", c + 1);
    fits_update_key(f, TLONG, key, (void *)&tbcol[c], NULL, &st);
  }
  long naxis1 = W;
  fits_update_key(f, TLONG, "NAXIS1", &naxis1, NULL, &st);
  double v;
  v = 2.5;  fits_write_key(f, TDOUBLE, "TSCAL1", &v, NULL, &st);
  v = -1;   fits_write_key(f, TDOUBLE, "TZERO1", &v, NULL, &st);
  v = 0.5;  fits_write_key(f, TDOUBLE, "TSCAL2", &v, NULL, &st);
  v = 100;  fits_write_key(f, TDOUBLE, "TZERO2", &v, NULL, &st);
  v = -2;   fits_write_key(f, TDOUBLE, "TSCAL3", &v, NULL, &st);
  v = 0.25; fits_write_key(f, TDOUBLE, "TZERO3", &v, NULL, &st);
  fits_write_key(f, TSTRING, "TNULL5", "-99", NULL, &st);
  check(st, "ascii scaled keys");
  fits_set_hdustruc(f, &st);
  for (int r = 0; r < NROW; r++) {
    char buf[W + 1];
    memset(buf, ' ', W);
    buf[W] = 0;
    for (int c = 0; c < SNC; c++) {
      const char *s = rows[r][c];
      int left = s[0] == '<';
      if (left) s++;
      int n = (int)strlen(s);
      if (n > width[c]) { fprintf(stderr, "scaled: field too wide row %d col %d\n", r, c); exit(1); }
      int at = (int)tbcol[c] - 1 + (left ? 0 : width[c] - n);
      memcpy(buf + at, s, (size_t)n);
    }
    fits_insert_rows(f, r, 1, &st);
    fits_write_tblbytes(f, r + 1, 1, W, (unsigned char *)buf, &st);
    check(st, "ascii scaled row");
  }
  closef(f);
}

/* -------------------------------------------------------- target ranges */
/* Two-row columns whose physical values sit exactly on, or one step past,
 * the range of every Go integer type, in three arithmetic paths:
 *   K_*   1K, no scaling (exact integers)
 *   D_*   1D (the double itself)
 *   S_*   1K with TSCAL -1 (float arithmetic, stored = -value)
 * For each of I8 I16 I32 U8 U16 U32: _IN {lo, hi}, _LO {lo-1, lo},
 * _HI {hi, hi+1}. The 64-bit ranges use doubles that exist: D_I64 and
 * S2_I64 (1K, TSCAL 2) hold -2^63 and 2^63-1024 (in), -2^63-2048 and 2^63
 * (out); D_U64 holds 2^64-2048 (in) and 2^64 (out), S2_U64 2^64-2048
 * and 2*(2^63-1) = 2^64-2 (both in range exactly; the float64 arithmetic
 * of the latter is 2^64). Plus
 * fractions, +-Inf, NaN, full-range B/I/J and integral E columns, and
 * VLAs on and past the int16 range. */
#define RMAX 96
static int rn;
static char rname[RMAX][24], rform[RMAX][16];
static double rscal[RMAX];
static int rkind[RMAX]; /* 0 = stored LONGLONG, 1 = double, 2 = VLA/0P (written below) */
static LONGLONG rk[RMAX][2];
static double rd[RMAX][2];

static void radd(const char *name, const char *form, double scal, int kind) {
  snprintf(rname[rn], 24, "%s", name);
  snprintf(rform[rn], 16, "%s", form);
  rscal[rn] = scal;
  rkind[rn] = kind;
}
static void radd_k(const char *name, const char *form, double scal, LONGLONG a, LONGLONG b) {
  radd(name, form, scal, 0);
  rk[rn][0] = a;
  rk[rn][1] = b;
  rn++;
}
static void radd_d(const char *name, const char *form, double a, double b) {
  radd(name, form, 1, 1);
  rd[rn][0] = a;
  rd[rn][1] = b;
  rn++;
}
static void radd_v(const char *name, const char *form) {
  radd(name, form, 1, 2);
  rn++;
}
static int rcol(const char *name) {
  for (int c = 0; c < rn; c++)
    if (!strcmp(rname[c], name)) return c + 1;
  fprintf(stderr, "gen_tables: no column %s\n", name);
  exit(1);
}

static void ranges(const char *dir) {
  static const struct { const char *n; LONGLONG lo, hi; } b[] = {
      {"I8", -128, 127}, {"I16", -32768, 32767}, {"I32", -2147483648LL, 2147483647LL},
      {"U8", 0, 255},    {"U16", 0, 65535},      {"U32", 0, 4294967295LL},
  };
  char nm[24];
  rn = 0;
  for (int i = 0; i < 6; i++) {
    LONGLONG lo = b[i].lo, hi = b[i].hi;
    snprintf(nm, 24, "K_%s_IN", b[i].n); radd_k(nm, "1K", 1, lo, hi);
    snprintf(nm, 24, "K_%s_LO", b[i].n); radd_k(nm, "1K", 1, lo - 1, lo);
    snprintf(nm, 24, "K_%s_HI", b[i].n); radd_k(nm, "1K", 1, hi, hi + 1);
    snprintf(nm, 24, "D_%s_IN", b[i].n); radd_d(nm, "1D", (double)lo, (double)hi);
    snprintf(nm, 24, "D_%s_LO", b[i].n); radd_d(nm, "1D", (double)(lo - 1), (double)lo);
    snprintf(nm, 24, "D_%s_HI", b[i].n); radd_d(nm, "1D", (double)hi, (double)(hi + 1));
    snprintf(nm, 24, "S_%s_IN", b[i].n); radd_k(nm, "1K", -1, -lo, -hi);
    snprintf(nm, 24, "S_%s_LO", b[i].n); radd_k(nm, "1K", -1, -(lo - 1), -lo);
    snprintf(nm, 24, "S_%s_HI", b[i].n); radd_k(nm, "1K", -1, -hi, -(hi + 1));
  }
  const double p63 = 9223372036854775808.0, p64 = 18446744073709551616.0;
  radd_d("D_I64_IN", "1D", -p63, p63 - 1024);
  radd_d("D_I64_LO", "1D", -p63 - 2048, 0);
  radd_d("D_I64_HI", "1D", 0, p63);
  radd_d("D_U64_IN", "1D", 0, p64 - 2048);
  radd_d("D_U64_LO", "1D", -1, 0);
  radd_d("D_U64_HI", "1D", 0, p64);
  radd_k("S2_I64_IN", "1K", 2, -4611686018427387904LL, 4611686018427387392LL);
  radd_k("S2_I64_LO", "1K", 2, -4611686018427388928LL, 0);
  radd_k("S2_I64_HI", "1K", 2, 0, 4611686018427387904LL);
  radd_k("S2_U64_IN", "1K", 2, 0, 9223372036854774784LL);
  radd_k("S2_U64_LO", "1K", 2, -1, 0);
  radd_k("S2_U64_HI", "1K", 2, 0, LLONG_MAX);
  radd_d("D_FRAC", "1D", 0.5, 1);
  radd_d("D_INF", "1D", INFINITY, 1);
  radd_d("D_NINF", "1D", -INFINITY, 1);
  radd_d("D_NAN", "1D", NAN, 1);
  radd_k("S_FRAC", "1K", 0.5, 1, 2);
  radd_k("B_FULL", "1B", 1, 0, 255);
  radd_k("I_FULL", "1I", 1, -32768, 32767);
  radd_k("J_FULL", "1J", 1, -2147483648LL, 2147483647LL);
  radd_d("E_INT", "1E", 16777216.0, -3);
  radd_v("V_I16_IN", "1PJ(2)");
  radd_v("V_I16_HI", "1QJ(1)");

  fitsfile *f = create(dir, "cf_bintable_ranges.fits");
  int st = 0;
  char *ttype[RMAX], *tform[RMAX];
  for (int c = 0; c < rn; c++) {
    ttype[c] = rname[c];
    tform[c] = rform[c];
  }
  fits_create_tbl(f, BINARY_TBL, 2, rn, ttype, tform, NULL, "RANGES", &st);
  check(st, "create ranges");
  char key[16];
  for (int c = 0; c < rn; c++) {
    if (rscal[c] != 1) {
      snprintf(key, sizeof key, "TSCAL%d", c + 1);
      fits_write_key(f, TDOUBLE, key, &rscal[c], NULL, &st);
    }
  }
  check(st, "ranges keys");
  fits_set_hdustruc(f, &st);
  for (int c = 0; c < rn; c++) {
    if (rkind[c] == 0) {
      /* stored integers, written with scaling switched off */
      fits_set_tscale(f, c + 1, 1.0, 0.0, &st);
      fits_write_col(f, TLONGLONG, c + 1, 1, 1, 2, rk[c], &st);
    } else if (rkind[c] == 1) {
      fits_write_col(f, TDOUBLE, c + 1, 1, 1, 2, rd[c], &st);
    }
    check(st, rname[c]);
  }
  LONGLONG vin[2] = {-32768, 32767}, vzero[1] = {0}, vhi[1] = {32768};
  fits_write_col(f, TLONGLONG, rcol("V_I16_IN"), 1, 1, 2, vin, &st);
  fits_write_col(f, TLONGLONG, rcol("V_I16_IN"), 2, 1, 1, vzero, &st);
  fits_write_col(f, TLONGLONG, rcol("V_I16_HI"), 1, 1, 1, vhi, &st);
  fits_write_descript(f, rcol("V_I16_HI"), 2, 0, 0, &st);
  check(st, "ranges VLAs");
  closef(f);
}

/* ------------------------------------------- integer TZERO beyond 2^53 */
/* Integer TZERO literals that no double holds (2^63-1, 2^62+1, -(2^63-1))
 * on K and J columns, and one written as a float literal
 * ("4611686018427387905.0"). Stored values are written with scaling off;
 * the physical value is stored + TZERO exactly. */
static void tzerobig(const char *dir) {
  enum { NC = 7, NRW = 3 };
  fitsfile *f = create(dir, "cf_bintable_tzero_big.fits");
  int st = 0;
  char *ttype[NC] = {"KMAX", "KMAXU", "KHALF", "KNEG", "KNEGB", "JHALF", "KFLT"};
  char *tform[NC] = {"1K", "1K", "1K", "1K", "1K", "1J", "1K"};
  fits_create_tbl(f, BINARY_TBL, NRW, NC, ttype, tform, NULL, "TZEROBIG", &st);
  check(st, "create tzero_big");
  LONGLONG z[NC - 1] = {LLONG_MAX, LLONG_MAX, 4611686018427387905LL, -LLONG_MAX, -LLONG_MAX, 4611686018427387905LL};
  char key[16];
  for (int c = 0; c < NC - 1; c++) {
    snprintf(key, sizeof key, "TZERO%d", c + 1);
    fits_write_key(f, TLONGLONG, key, &z[c], NULL, &st);
  }
  fits_write_record(f, "TZERO7  =  4611686018427387905.0 / float literal of 2^62+1", &st);
  check(st, "tzero_big keys");
  fits_set_hdustruc(f, &st);
  for (int c = 1; c <= NC; c++) fits_set_tscale(f, c, 1.0, 0.0, &st);
  LONGLONG v[NC][NRW] = {
      {0, -1, LLONG_MIN},                  /* 2^63-1, 2^63-2, -1 */
      {1, 0, LLONG_MAX},                   /* 2^63, 2^63-1, 2^64-2 */
      {0, 1, -4611686018427387905LL},      /* 2^62+1, 2^62+2, 0 */
      {0, -1, 1},                          /* -2^63+1, -2^63, -2^63+2 */
      {-2, 0, LLONG_MIN},                  /* below -2^63 */
      {1, -1, 2147483647LL},               /* 2^62+2, 2^62, 2^62+2^31+2^31-... */
      {0, 1, 2},
  };
  for (int c = 0; c < NC; c++) fits_write_col(f, TLONGLONG, c + 1, 1, 1, NRW, v[c], &st);
  check(st, "tzero_big");
  closef(f);
}

/* ---------------------------------------------- ASCII exact integers */
/* Integer-valued ASCII fields outside the int64 range, at its edges, with
 * a fractional part that is all zeros, and with a fraction next to an
 * integer-valued double; an I20 with integer TZERO 1000; and fields with
 * a sign, point or exponent but no digits ('.', 'E5', '+'), which cfitsio
 * reads as 0. */
#define XNC 9
static void ascii_exact(const char *dir) {
  static const char *ttype[XNC] = {"IU", "IN", "IHI", "ILO", "FU", "FN", "FX", "IZ", "DIG"};
  static const char *tforms[XNC] = {"I20", "I20", "I20", "I20", "F24.3", "F24.3", "F24.3", "I20", "F8.2"};
  static const int width[XNC] = {20, 20, 20, 20, 24, 24, 24, 20, 8};
  static const char *rows[][XNC] = {
      {"18446744073709551615", "-9223372036854775808", "18446744073709551616", "-9223372036854775809",
       "18446744073709551615.000", "-9223372036854775808.000", "9223372036854775807.500", "9223372036854775807", "."},
      {"9223372036854775808", "9223372036854775807", "1", "1",
       "9223372036854775809.000", "9223372036854775807.000", "1.000", "-5", "E5"},
      {"0", "-1", "2", "2", "0.000", "12.000", "2.000", "0", "+"},
  };
  enum { NROW = 3 };
  long tbcol[XNC];
  int w = 0;
  for (int c = 0; c < XNC; c++) {
    tbcol[c] = w + 1;
    w += width[c] + 1;
  }
  w -= 1;
  fitsfile *f = create(dir, "cf_ascii_exact.fits");
  int st = 0;
  fits_create_tbl(f, ASCII_TBL, 0, XNC, (char **)ttype, (char **)tforms, NULL, "EXACT", &st);
  check(st, "create ascii exact");
  char key[16];
  for (int c = 0; c < XNC; c++) {
    snprintf(key, sizeof key, "TBCOL%d", c + 1);
    fits_update_key(f, TLONG, key, &tbcol[c], NULL, &st);
  }
  long naxis1 = w;
  fits_update_key(f, TLONG, "NAXIS1", &naxis1, NULL, &st);
  LONGLONG z = 1000;
  fits_write_key(f, TLONGLONG, "TZERO8", &z, NULL, &st);
  check(st, "ascii exact keys");
  fits_set_hdustruc(f, &st);
  for (int r = 0; r < NROW; r++) {
    char buf[512];
    memset(buf, ' ', (size_t)w);
    for (int c = 0; c < XNC; c++) {
      const char *s = rows[r][c];
      int n = (int)strlen(s);
      if (n > width[c]) { fprintf(stderr, "exact: field too wide row %d col %d\n", r, c); exit(1); }
      memcpy(buf + tbcol[c] - 1 + (width[c] - n), s, (size_t)n);
    }
    fits_insert_rows(f, r, 1, &st);
    fits_write_tblbytes(f, r + 1, 1, w, (unsigned char *)buf, &st);
    check(st, "ascii exact row");
  }
  closef(f);
}

/* ------------------------------------------ exact scaling, binary */
/* Columns whose physical value (stored*TSCAL + TZERO, with the literals
 * read as exact decimals) is not what float64 arithmetic gives: a
 * fractional TZERO whose sum rounds to an integral double (KHALF, DHALF),
 * TSCAL 2 on stored values beyond 2^53 (KS2), integer TZERO literals
 * outside int64 (J63 = 2^63 on a J column, K63P1 = 2^63+1, which rounds
 * to 2^63 as a double; stored 1024 makes 2^63+1025, whose correctly
 * rounded double is 2^63+2048 while 1024.0 + 2^63 rounds to 2^63; K264 =
 * 2^64), TSCAL 0.1 on multiples of 10
 * (exactly integral, not in float64), an exact float column (ESC3) and
 * TSCAL 4 reaching exactly 2^64, one past the uint64 range (KS4); and
 * TNULL = 0 on a J column holding 0, 7, 0 (JNUL0).
 * Every TSCAL/TZERO card is written as a raw record so the literal is
 * exactly the text here; stored values are written with scaling off. */
static void scaledexact(const char *dir) {
  enum { NC = 10, NRW = 3 };
  fitsfile *f = create(dir, "cf_bintable_scaled_exact.fits");
  int st = 0;
  char *ttype[NC] = {"KHALF", "KS2", "J63", "K63P1", "K264", "KTENTH", "DHALF", "ESC3", "KS4", "JNUL0"};
  char *tform[NC] = {"1K", "1K", "1J", "1K", "1K", "1K", "1D", "1E", "1K", "1J"};
  fits_create_tbl(f, BINARY_TBL, NRW, NC, ttype, tform, NULL, "SCALEXACT", &st);
  check(st, "create scaled_exact");
  const char *cards[] = {
      "TZERO1  =                  0.5 / fractional TZERO",
      "TSCAL2  =                    2",
      "TZERO3  =  9223372036854775808 / 2^63 on a J column",
      "TZERO4  =  9223372036854775809 / 2^63+1 (2^63 as a double)",
      "TZERO5  = 18446744073709551616 / 2^64",
      "TSCAL6  =                  0.1",
      "TZERO7  =                  0.5",
      "TSCAL8  =                    3",
      "TZERO8  =                  0.5",
      "TSCAL9  =                    4",
      "TNULL10 =                    0 / a stored 0 is null",
  };
  for (size_t i = 0; i < sizeof cards / sizeof cards[0]; i++) fits_write_record(f, cards[i], &st);
  check(st, "scaled_exact keys");
  fits_set_hdustruc(f, &st);
  for (int c = 1; c <= NC; c++) fits_set_tscale(f, c, 1.0, 0.0, &st);
  LONGLONG k[6][NRW] = {
      {1152921504606846976LL, 1152921504606846978LL, 2305843009213693952LL}, /* 2^60(+2), 2^61: +0.5 */
      {2305843009213693953LL, -2305843009213693953LL, 3},                    /* *2: +-(2^62+2), 6 */
      {5, -5, 2147483647LL},                                                 /* 2^63+5, 2^63-5, 2^63+2^31-1 */
      {0, -1, 1024},                                                         /* 2^63+1, 2^63, 2^63+1025 */
      {-1, LLONG_MIN, -2},                                                   /* 2^64-1, 2^63, 2^64-2 */
      {10, 20, 30},                                                          /* 1, 2, 3 */
  };
  for (int c = 0; c < 6; c++) fits_write_col(f, TLONGLONG, c + 1, 1, 1, NRW, k[c], &st);
  double d[NRW] = {9007199254740992.0, 9007199254740994.0, 18014398509481984.0}; /* 2^53, 2^53+2, 2^54: +0.5 */
  fits_write_col(f, TDOUBLE, 7, 1, 1, NRW, d, &st);
  float e[NRW] = {0.5f, 1.5f, 2.5f}; /* *3 + 0.5: 2, 5, 8 */
  fits_write_col(f, TFLOAT, 8, 1, 1, NRW, e, &st);
  LONGLONG k4[NRW] = {4611686018427387904LL, 1, -1}; /* *4: 2^64 (one past uint64), 4, -4 */
  fits_write_col(f, TLONGLONG, 9, 1, 1, NRW, k4, &st);
  int jn[NRW] = {0, 7, 0};
  fits_write_col(f, TINT, 10, 1, 1, NRW, jn, &st);
  check(st, "scaled_exact");
  closef(f);
}

/* ------------------------------------------- exact scaling, ASCII */
/* The ASCII counterparts: I20 with TSCAL 2 beyond 2^53 (IS2), with a
 * fractional TZERO (IH) and with TZERO 2^63 (I63); F30.3 fractions under
 * an integer TZERO (FZ: exact decimal + TZERO, correctly rounded);
 * F24.20 fractions of 21 significant digits under TZERO 2^53 (FZ53: the
 * exact sum rounds differently from the float64 of the field plus
 * TZERO); E12.4 fields beyond the double range (EU: 1.5E-900 is not 0
 * for an integer target); I4 with a TNULL longer than the field (INW:
 * the standard never matches it); fractions under TSCAL 10 that are
 * exact integers but not in float64 (FS10: 0.3*10 = 3.0000000000000004)
 * and 21-digit fractions under TSCAL 2 whose float64 is integral but whose
 * exact value is not (FS2). */
#define ENC 9
static void ascii_scaled_exact(const char *dir) {
  static const char *ttype[ENC] = {"IS2", "IH", "I63", "FZ", "FZ53", "EU", "INW", "FS10", "FS2"};
  static const char *tforms[ENC] = {"I20", "I20", "I20", "F30.3", "F24.20", "E12.4", "I4", "F6.1", "F24.20"};
  static const int width[ENC] = {20, 20, 20, 30, 24, 12, 4, 6, 24};
  static const char *rows[][ENC] = {
      {"4611686018427387905", "1152921504606846976", "5", "0.500", "1.00000000000000000001", "1.5E-900", "-999", "0.3", "1.00000000000000000001"},
      {"0", "1152921504606846978", "-5", "-1000.250", "3.00000000000000000001", "2.0E+00", "12", "0.7", "3.00000000000000000001"},
      {"3", "2305843009213693952", "0", "9223372036854774807.500", "0.00000000000000000001", "-2.5E-850", "-99", "1.1", "2.5"},
  };
  enum { NROW = 3 };
  long tbcol[ENC];
  int w = 0;
  for (int c = 0; c < ENC; c++) {
    tbcol[c] = w + 1;
    w += width[c] + 1;
  }
  w -= 1;
  fitsfile *f = create(dir, "cf_ascii_scaled_exact.fits");
  int st = 0;
  fits_create_tbl(f, ASCII_TBL, 0, ENC, (char **)ttype, (char **)tforms, NULL, "ASCALEXACT", &st);
  check(st, "create ascii scaled exact");
  char key[16];
  for (int c = 0; c < ENC; c++) {
    snprintf(key, sizeof key, "TBCOL%d", c + 1);
    fits_update_key(f, TLONG, key, &tbcol[c], NULL, &st);
  }
  long naxis1 = w;
  fits_update_key(f, TLONG, "NAXIS1", &naxis1, NULL, &st);
  const char *cards[] = {
      "TSCAL1  =                    2",
      "TZERO2  =                  0.5",
      "TZERO3  =  9223372036854775808",
      "TZERO4  =                 1000",
      "TZERO5  =     9007199254740992 / 2^53",
      "TNULL7  = '-99999'           / longer than the I4 field",
      "TSCAL8  =                   10",
      "TSCAL9  =                    2",
  };
  for (size_t i = 0; i < sizeof cards / sizeof cards[0]; i++) fits_write_record(f, cards[i], &st);
  check(st, "ascii scaled exact keys");
  fits_set_hdustruc(f, &st);
  for (int r = 0; r < NROW; r++) {
    char buf[512];
    memset(buf, ' ', (size_t)w);
    for (int c = 0; c < ENC; c++) {
      const char *s = rows[r][c];
      int n = (int)strlen(s);
      if (n > width[c]) { fprintf(stderr, "scaled exact: field too wide row %d col %d\n", r, c); exit(1); }
      memcpy(buf + tbcol[c] - 1 + (width[c] - n), s, (size_t)n);
    }
    fits_insert_rows(f, r, 1, &st);
    fits_write_tblbytes(f, r + 1, 1, w, (unsigned char *)buf, &st);
    check(st, "ascii scaled exact row");
  }
  closef(f);
}

/* ---------------------------------------- scaling edge cases, ASCII */
/* Negative zeros and big TZERO literals: a field with a minus sign and no
 * non-zero digit ("-0", "-0.0", "-", "-.") under an integer TZERO no
 * double holds (NZ53 = 2^53+1, NZ63 = 2^63+1, NZK = 1000), a fractional
 * TZERO (NZF) and TSCAL 3 (NZS): the physical value is exactly TZERO
 * (0*TSCAL + TZERO), not the double of TZERO. I30 integers beyond uint64
 * under TZERO = 2^64+2051 (BZ: the exact sum 2^65+4098 rounds to
 * 2^65+8192, the sum of the two doubles to 2^65; -2^64 + TZERO = 2051)
 * and F30.1 fractions under TZERO = 2^64+2048, a double rounding tie
 * (BZF: 0.5 + TZERO is 2^64+4096 correctly rounded, 2^64 through
 * doubles). */
#define GNC 7
static void ascii_scaled_edge(const char *dir) {
  static const char *ttype[GNC] = {"NZ53", "NZ63", "NZK", "NZF", "NZS", "BZ", "BZF"};
  static const char *tforms[GNC] = {"I20", "I20", "I20", "F20.3", "F20.3", "I30", "F30.1"};
  static const int width[GNC] = {20, 20, 20, 20, 20, 30, 30};
  static const char *rows[][GNC] = {
      {"-0", "-0", "-0", "-0.0", "-0.0", "18446744073709553663", "0.5"},
      {"0", "0", "-", "-.", "0.0", "-18446744073709551616", "-0.5"},
      {"-", "-0", "0", "-0.000", "-0.000", "-18446744073709553667", "1.5"},
  };
  enum { NROW = 3 };
  long tbcol[GNC];
  int w = 0;
  for (int c = 0; c < GNC; c++) {
    tbcol[c] = w + 1;
    w += width[c] + 1;
  }
  w -= 1;
  fitsfile *f = create(dir, "cf_ascii_scaled_edge.fits");
  int st = 0;
  fits_create_tbl(f, ASCII_TBL, 0, GNC, (char **)ttype, (char **)tforms, NULL, "ASCALEDGE", &st);
  check(st, "create ascii scaled edge");
  char key[16];
  for (int c = 0; c < GNC; c++) {
    snprintf(key, sizeof key, "TBCOL%d", c + 1);
    fits_update_key(f, TLONG, key, &tbcol[c], NULL, &st);
  }
  long naxis1 = w;
  fits_update_key(f, TLONG, "NAXIS1", &naxis1, NULL, &st);
  const char *cards[] = {
      "TZERO1  =     9007199254740993 / 2^53+1",
      "TZERO2  =  9223372036854775809 / 2^63+1",
      "TZERO3  =                 1000",
      "TZERO4  = 1.00000000000000000001",
      "TSCAL5  =                    3",
      "TZERO5  =     9007199254740993 / 2^53+1",
      "TZERO6  = 18446744073709553667 / 2^64+2051",
      "TZERO7  = 18446744073709553664 / 2^64+2048",
  };
  for (size_t i = 0; i < sizeof cards / sizeof cards[0]; i++) fits_write_record(f, cards[i], &st);
  check(st, "ascii scaled edge keys");
  fits_set_hdustruc(f, &st);
  for (int r = 0; r < NROW; r++) {
    char buf[512];
    memset(buf, ' ', (size_t)w);
    for (int c = 0; c < GNC; c++) {
      const char *s = rows[r][c];
      int n = (int)strlen(s);
      if (n > width[c]) { fprintf(stderr, "scaled edge: field too wide row %d col %d\n", r, c); exit(1); }
      memcpy(buf + tbcol[c] - 1 + (width[c] - n), s, (size_t)n);
    }
    fits_insert_rows(f, r, 1, &st);
    fits_write_tblbytes(f, r + 1, 1, w, (unsigned char *)buf, &st);
    check(st, "ascii scaled edge row");
  }
  closef(f);
}

/* ----------------------------------- scaling edge cases, binary */
/* TSCAL/TZERO literals with a D exponent (KD18: TZERO 4.611686018427387905D18
 * = 2^62+1, which no double holds; KSD: TSCAL 1.0D-1 on 10, 20, 30, exactly
 * 1, 2, 3; JD: TSCAL 5.0D-1, TZERO -2.5D0) and +-Inf / NaN stored in scaled
 * E and D columns (EINF: TSCAL 2, TZERO 0.5; DINF: TSCAL -3, TZERO 1;
 * DINFZ: TZERO 0.5 only; EINFV: a 2E vector with TSCAL 0.5). The literal
 * cards are written as raw records; stored values with scaling off. */
static void scaledspecial(const char *dir) {
  enum { NC = 7, NRW = 4 };
  fitsfile *f = create(dir, "cf_bintable_scaled_special.fits");
  int st = 0;
  char *ttype[NC] = {"KD18", "KSD", "JD", "EINF", "DINF", "DINFZ", "EINFV"};
  char *tform[NC] = {"1K", "1K", "1J", "1E", "1D", "1D", "2E"};
  fits_create_tbl(f, BINARY_TBL, NRW, NC, ttype, tform, NULL, "SCALSPECIAL", &st);
  check(st, "create scaled_special");
  const char *cards[] = {
      "TZERO1  = 4.611686018427387905D18 / 2^62+1 with a D exponent",
      "TSCAL2  =               1.0D-1",
      "TSCAL3  =               5.0D-1",
      "TZERO3  =               -2.5D0",
      "TSCAL4  =                    2",
      "TZERO4  =                  0.5",
      "TSCAL5  =                   -3",
      "TZERO5  =                    1",
      "TZERO6  =                  0.5",
      "TSCAL7  =                  0.5",
  };
  for (size_t i = 0; i < sizeof cards / sizeof cards[0]; i++) fits_write_record(f, cards[i], &st);
  check(st, "scaled_special keys");
  fits_set_hdustruc(f, &st);
  for (int c = 1; c <= NC; c++) fits_set_tscale(f, c, 1.0, 0.0, &st);
  LONGLONG kd[NRW] = {0, 1, -4611686018427387905LL, -1};   /* 2^62+1, 2^62+2, 0, 2^62 */
  LONGLONG ks[NRW] = {10, 20, 30, -70};                     /* 1, 2, 3, -7 */
  int jd[NRW] = {5, 7, 2147483647, -2147483647};            /* 0, 1, 1073741821, -1073741826 */
  float einf[NRW] = {INFINITY, -INFINITY, 1.25f, NAN};      /* Inf, -Inf, 3, null */
  double dinf[NRW] = {-INFINITY, INFINITY, 2.0, NAN};       /* Inf, -Inf, -5, null */
  double dinfz[NRW] = {INFINITY, 1.5, -INFINITY, NAN};      /* Inf, 2, -Inf, null */
  float einfv[2 * NRW] = {INFINITY, 4.0f, 2.0f, -INFINITY, 6.0f, 8.0f, NAN, 1.0f};
  fits_write_col(f, TLONGLONG, 1, 1, 1, NRW, kd, &st);
  fits_write_col(f, TLONGLONG, 2, 1, 1, NRW, ks, &st);
  fits_write_col(f, TINT, 3, 1, 1, NRW, jd, &st);
  fits_write_col(f, TFLOAT, 4, 1, 1, NRW, einf, &st);
  fits_write_col(f, TDOUBLE, 5, 1, 1, NRW, dinf, &st);
  fits_write_col(f, TDOUBLE, 6, 1, 1, NRW, dinfz, &st);
  fits_write_col(f, TFLOAT, 7, 1, 1, 2 * NRW, einfv, &st);
  check(st, "scaled_special");
  closef(f);
}

/* ------------------------------------- rounding under an integer TZERO */
/* ASCII fields whose physical value (field + TZERO, exactly) rounds to a
 * different double than the field's double plus TZERO's double:
 * BIGI: I22 integers beyond int64 under TZERO -2 (18446744073709553665 =
 * 2^64+2049 is 2^64+2047 exactly, 2^64 correctly rounded, 2^64+4096
 * through doubles; 18446745173221179395 is 2^64+2^40+1, whose double
 * 2^64+2^40 is a float32 tie: float32 2^64 through the double, 2^64+2^41
 * nearest); TINY: E12.4 exponents beyond -800 under TZERO =
 * 2^53+1, which no double holds ("1.5E-900" is just above the rounding
 * tie 2^53+1: 2^53+2; "-1.5E-900" just below: 2^53); WIDE: F820.812
 * fields of 802 significant digits (beyond the 800-digit exact path) a
 * hair above / below 0 under the same TZERO; FRAC: F30.3 fractions whose
 * correctly rounded float64 is a float32 rounding tie (536870911.5 +
 * 2^53+1 = 2^53+2^29+0.5: float64 2^53+2^29, then float32 2^53; the
 * float32 nearest the exact value is 2^53+2^30); EUS: E12.4 exponents
 * beyond +-800 on the float64 path (TSCAL 2): 0, -Inf, +Inf, 0. */
#define RNC 5
static void ascii_round(const char *dir) {
  static const char *ttype[RNC] = {"BIGI", "TINY", "WIDE", "FRAC", "EUS"};
  static const char *tforms[RNC] = {"I22", "E12.4", "F820.812", "F30.3", "E12.4"};
  static const int width[RNC] = {22, 12, 820, 30, 12};
  static char wide[3][821];
  /* 0.(10 zeros)1(800 zeros)1: 802 significant digits */
  for (int k = 0; k < 2; k++) {
    char *p = wide[k];
    if (k == 1) *p++ = '-';
    *p++ = '0';
    *p++ = '.';
    for (int i = 0; i < 10; i++) *p++ = '0';
    *p++ = '1';
    for (int i = 0; i < 800; i++) *p++ = '0';
    *p++ = '1';
    *p = 0;
  }
  strcpy(wide[2], "-0.25");
  const char *rows[][RNC] = {
      {"18446744073709553665", "1.5E-900", wide[0], "536870911.500", "1.5E-900"},
      {"18446744073709551617", "-1.5E-900", wide[1], "-536870912.500", "-1.5E+900"},
      {"-9223372036854775807", "0.5", wide[2], "0.250", "1.5E+900"},
      {"18446745173221179395", "-0.5", "0.75", "-0.250", "-1.5E-900"},
  };
  enum { NROW = 4 };
  long tbcol[RNC];
  int w = 0;
  for (int c = 0; c < RNC; c++) {
    tbcol[c] = w + 1;
    w += width[c] + 1;
  }
  w -= 1;
  fitsfile *f = create(dir, "cf_ascii_tzero_round.fits");
  int st = 0;
  fits_create_tbl(f, ASCII_TBL, 0, RNC, (char **)ttype, (char **)tforms, NULL, "ASCROUND", &st);
  check(st, "create ascii round");
  char key[16];
  for (int c = 0; c < RNC; c++) {
    snprintf(key, sizeof key, "TBCOL%d", c + 1);
    fits_update_key(f, TLONG, key, &tbcol[c], NULL, &st);
  }
  long naxis1 = w;
  fits_update_key(f, TLONG, "NAXIS1", &naxis1, NULL, &st);
  const char *cards[] = {
      "TZERO1  =                   -2",
      "TZERO2  =     9007199254740993 / 2^53+1",
      "TZERO3  =     9007199254740993 / 2^53+1",
      "TZERO4  =     9007199254740993 / 2^53+1",
      "TSCAL5  =                    2",
  };
  for (size_t i = 0; i < sizeof cards / sizeof cards[0]; i++) fits_write_record(f, cards[i], &st);
  check(st, "ascii round keys");
  fits_set_hdustruc(f, &st);
  for (int r = 0; r < NROW; r++) {
    char buf[2048];
    memset(buf, ' ', (size_t)w);
    for (int c = 0; c < RNC; c++) {
      const char *s = rows[r][c];
      int n = (int)strlen(s);
      if (n > width[c]) { fprintf(stderr, "ascii round: field too wide row %d col %d\n", r, c); exit(1); }
      memcpy(buf + tbcol[c] - 1 + (width[c] - n), s, (size_t)n);
    }
    fits_insert_rows(f, r, 1, &st);
    fits_write_tblbytes(f, r + 1, 1, w, (unsigned char *)buf, &st);
    check(st, "ascii round row");
  }
  closef(f);
}

/* Binary columns where the correctly rounded physical value is not what
 * double arithmetic (or a float32 of a double) gives. KF32: unscaled 1K
 * whose float32 nearest value differs from the float32 of its double
 * (2^60+2^36+1: double 2^60+2^36, a float32 tie that rounds to 2^60;
 * the float32 nearest 2^60+2^36+1 is 2^60+2^37). KZ: the same values
 * under TZERO 1000 (stored with scaling off). DZ53 / EZ53: 1D / 1E under
 * TZERO = 2^53+1, which no double holds (1.0 + TZERO is 2^53+2; the
 * double of TZERO is 2^53, plus 1.0 rounds to 2^53); 2^29 + TZERO =
 * 2^53+2^29+1 rounds to the float32 tie 2^53+2^29 as a double. DZ64: 1D
 * under TZERO = 2^64+2049 (-2.0 + TZERO is 2^64+2047: 2^64; the double of
 * TZERO is 2^64+4096), with -Inf and NaN (null). JF: unscaled 1J whose
 * float32 rounds. */
static void bin_round(const char *dir) {
  enum { NC = 6, NRW = 4 };
  fitsfile *f = create(dir, "cf_bintable_tzero_round.fits");
  int st = 0;
  char *ttype[NC] = {"KF32", "KZ", "DZ53", "EZ53", "DZ64", "JF"};
  char *tform[NC] = {"1K", "1K", "1D", "1E", "1D", "1J"};
  fits_create_tbl(f, BINARY_TBL, NRW, NC, ttype, tform, NULL, "TZROUND", &st);
  check(st, "create tzero round");
  const char *cards[] = {
      "TZERO2  =                 1000",
      "TZERO3  =     9007199254740993 / 2^53+1",
      "TZERO4  =     9007199254740993 / 2^53+1",
      "TZERO5  = 18446744073709553665 / 2^64+2049",
  };
  for (size_t i = 0; i < sizeof cards / sizeof cards[0]; i++) fits_write_record(f, cards[i], &st);
  check(st, "tzero round keys");
  fits_set_hdustruc(f, &st);
  for (int c = 1; c <= NC; c++) fits_set_tscale(f, c, 1.0, 0.0, &st);
  const LONGLONG a = (1LL << 60) + (1LL << 36);
  LONGLONG kf[NRW] = {a + 1, a, -(a + 1), 16777217};
  LONGLONG kz[NRW] = {a + 1 - 1000, (1LL << 62) + 1 - 1000, -(a + 1) - 1000, 9007199254740993LL - 1000};
  double dz[NRW] = {1.0, 0.5, 3.0, 536870912.0};
  float ez[NRW] = {-1.0f, 1.0f, 2.0f, 536870912.0f};
  double dz64[NRW] = {-2.0, 0.0, -INFINITY, NAN};
  int jf[NRW] = {16777217, -16777219, 2147483647, 0};
  fits_write_col(f, TLONGLONG, 1, 1, 1, NRW, kf, &st);
  fits_write_col(f, TLONGLONG, 2, 1, 1, NRW, kz, &st);
  fits_write_col(f, TDOUBLE, 3, 1, 1, NRW, dz, &st);
  fits_write_col(f, TFLOAT, 4, 1, 1, NRW, ez, &st);
  fits_write_col(f, TDOUBLE, 5, 1, 1, NRW, dz64, &st);
  fits_write_col(f, TINT, 6, 1, 1, NRW, jf, &st);
  check(st, "tzero round");
  closef(f);
}

/* ------------------------------------------- scaled complex (refused) */
/* TSCAL/TZERO on 1C and 1M columns (refused/cf_bintable_cplx_scaled.fits).
 * FITS 4.0 7.3.2 scales the real part by TSCAL and adds TZERO, the
 * imaginary part only by TSCAL: (1.5, 2.5) with TSCAL 2 / TZERO 0.5 is
 * (3.5, 5.0). The goldens record that cfitsio 4.7.0 adds TZERO to both
 * parts and astropy 8.0.2 drops the imaginary part, so the library
 * refuses these columns. Stored values written with scaling off. */
static void cplxscaled(const char *dir) {
  enum { NC = 2, NRW = 2 };
  fitsfile *f = create(dir, "refused/cf_bintable_cplx_scaled.fits");
  int st = 0;
  char *ttype[NC] = {"C", "M"};
  char *tform[NC] = {"1C", "1M"};
  fits_create_tbl(f, BINARY_TBL, NRW, NC, ttype, tform, NULL, "CPLXSCAL", &st);
  const char *cards[] = {
      "TSCAL1  =                    2",
      "TZERO1  =                  0.5",
      "TSCAL2  =                   -3",
      "TZERO2  =                 0.25",
  };
  for (size_t i = 0; i < sizeof cards / sizeof cards[0]; i++) fits_write_record(f, cards[i], &st);
  check(st, "cplx_scaled keys");
  fits_set_hdustruc(f, &st);
  for (int c = 1; c <= NC; c++) fits_set_tscale(f, c, 1.0, 0.0, &st);
  float cv[2 * NRW] = {1.5f, 2.5f, -1.0f, 4.0f};
  double mv[2 * NRW] = {1.5, 2.5, 0.0, -2.0};
  fits_write_col(f, TCOMPLEX, 1, 1, 1, NRW, cv, &st);
  fits_write_col(f, TDBLCOMPLEX, 2, 1, 1, NRW, mv, &st);
  check(st, "cplx_scaled");
  closef(f);
}

/* ------------------------------------------------------ string arrays */
/* Character columns with a multi-axis TDIM (FITS 4.0 7.3.2: the first
 * axis is the string length): S52 10A (5,2), S322 12A (3,2,2) and S8 8A
 * (8) (one string). TDIM is written before the data so cfitsio itself
 * splits the strings it is given. */
static void strarray(const char *dir) {
  enum { NC = 5, NRW = 3 };
  fitsfile *f = create(dir, "cf_bintable_strarray.fits");
  int st = 0;
  char *ttype[NC] = {"S52", "S322", "S8", "W8", "W42"};
  char *tform[NC] = {"10A", "12A", "8A", "8A2", "8A2"};
  fits_create_tbl(f, BINARY_TBL, NRW, NC, ttype, tform, NULL, "STRARRAY", &st);
  long d52[2] = {5, 2}, d322[3] = {3, 2, 2}, d8[1] = {8}, d42[2] = {4, 2};
  fits_write_tdim(f, 1, 2, d52, &st);
  fits_write_tdim(f, 2, 3, d322, &st);
  fits_write_tdim(f, 3, 1, d8, &st);
  fits_write_tdim(f, 4, 1, d8, &st);  /* 8A2 with TDIM (8): one 8-character string */
  fits_write_tdim(f, 5, 2, d42, &st); /* 8A2 with TDIM (4,2): two strings of 4 */
  check(st, "strarray tdim");
  fits_set_hdustruc(f, &st);
  char *s52[] = {"ab", "cdefg", "hello", "world", "", "x"};
  char *s322[] = {"a", "bb", "ccc", "d", "e", "f", "ggg", "hh", "i", "", "", "z"};
  char *s8[] = {"eight ch", "", "x y"};
  fits_write_col(f, TSTRING, 1, 1, 1, 6, s52, &st);
  fits_write_col(f, TSTRING, 2, 1, 1, 12, s322, &st);
  fits_write_col(f, TSTRING, 3, 1, 1, 3, s8, &st);
  char *w8[] = {"abcdefgh", "x", ""};
  char *w42[] = {"abcd", "efgh", "ab", "c", "", "wxyz"};
  fits_write_col(f, TSTRING, 4, 1, 1, 3, w8, &st);
  fits_write_col(f, TSTRING, 5, 1, 1, 6, w42, &st);
  check(st, "strarray");
  closef(f);
}

/* --------------------------------------------------------- A3DTABLE */
static void a3dtable(const char *dir, const char *name, const char *xt) {
  fitsfile *f = create(dir, name);
  int st = 0;
  char *ttype[] = {"J1", "E1", "A8"}, *tform[] = {"1J", "1E", "8A"};
  int j[3] = {-7, 0, 2147483647};
  float e[3] = {1.5f, -0.25f, 3e38f};
  char *a[3] = {"alpha", "", "eightchr"};
  fits_create_tbl(f, BINARY_TBL, 3, 3, ttype, tform, NULL, "LEGACY", &st);
  fits_write_col(f, TINT, 1, 1, 1, 3, j, &st);
  fits_write_col(f, TFLOAT, 2, 1, 1, 3, e, &st);
  fits_write_col(f, TSTRING, 3, 1, 1, 3, a, &st);
  check(st, "a3dtable");
  closef(f);
  /* Replace the value of the extension's first card in place. */
  char path[4096];
  snprintf(path, sizeof path, "%s/%s", dir, name);
  FILE *fp = fopen(path, "r+b");
  char card[21] = {0};
  if (!fp || fseek(fp, 2880, SEEK_SET) || fread(card, 1, 20, fp) != 20 ||
      memcmp(card, "XTENSION= 'BINTABLE'", 20) || fseek(fp, 2880 + 11, SEEK_SET) ||
      fwrite(xt, 1, 8, fp) != 8 || fclose(fp)) {
    fprintf(stderr, "gen_tables: patching %s failed\n", path);
    exit(1);
  }
}

/* ------------------------------------------------------- bad logicals */
/* A 1L and a 3L column written by cfitsio, then two logical bytes
 * byte-patched to 0x01 (neither 'T', 'F' nor 0): row 1 of L1 and element 2
 * of row 2 of L3. The library refuses them ("invalid logical byte"); the
 * golden records what cfitsio makes of them. */
static void badlogical(const char *dir) {
  const char *name = "cf_bintable_badlogical.fits";
  fitsfile *f = create(dir, name);
  int st = 0;
  char *ttype[] = {"L1", "L3", "J1"}, *tform[] = {"1L", "3L", "1J"};
  char l1[3] = {1, 0, 1}, l3[9] = {1, 1, 0, 0, 1, 0, 1, 1, 1};
  int j[3] = {1, 2, 3};
  LONGLONG hstart, dstart, dend;
  fits_create_tbl(f, BINARY_TBL, 3, 3, ttype, tform, NULL, "BADLOG", &st);
  fits_write_col(f, TLOGICAL, 1, 1, 1, 3, l1, &st);
  fits_write_col(f, TLOGICAL, 2, 1, 1, 9, l3, &st);
  fits_write_col(f, TINT, 3, 1, 1, 3, j, &st);
  fits_get_hduaddrll(f, &hstart, &dstart, &dend, &st);
  check(st, "badlogical");
  closef(f);
  char path[4096];
  snprintf(path, sizeof path, "%s/%s", dir, name);
  FILE *fp = fopen(path, "r+b");
  /* Row width 1 + 3 + 4 = 8 bytes. */
  if (!fp || fseek(fp, (long)dstart + 0 * 8 + 0, SEEK_SET) || fputc(1, fp) == EOF ||
      fseek(fp, (long)dstart + 1 * 8 + 1 + 1, SEEK_SET) || fputc(1, fp) == EOF || fclose(fp)) {
    fprintf(stderr, "gen_tables: patching %s failed\n", path);
    exit(1);
  }
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: gen_tables <outdir>\n");
    return 2;
  }
  fixed(argv[1]);
  uint64file(argv[1]);
  vla(argv[1], "cf_vla_p.fits", 'P');
  vla(argv[1], "cf_vla_q.fits", 'Q');
  ascii(argv[1], "cf_ascii.fits", ascii_std, (int)(sizeof ascii_std / sizeof ascii_std[0]));
  ascii(argv[1], "cf_ascii_implied.fits", ascii_implied, (int)(sizeof ascii_implied / sizeof ascii_implied[0]));
  kzero(argv[1]);
  substr(argv[1]);
  ascii_scaled(argv[1]);
  ranges(argv[1]);
  tzerobig(argv[1]);
  ascii_exact(argv[1]);
  a3dtable(argv[1], "cf_bintable_a3dtable.fits", "A3DTABLE");
  a3dtable(argv[1], "cf_bintable_3dtable.fits", "3DTABLE ");
  badlogical(argv[1]);
  scaledexact(argv[1]);
  ascii_scaled_exact(argv[1]);
  strarray(argv[1]);
  ascii_scaled_edge(argv[1]);
  scaledspecial(argv[1]);
  cplxscaled(argv[1]);
  ascii_round(argv[1]);
  bin_round(argv[1]);
  return 0;
}
