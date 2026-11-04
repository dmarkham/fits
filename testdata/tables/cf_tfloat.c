/*
 * cf_tfloat -- dumps cfitsio's float32 (TFLOAT) read of every numeric table
 * column in a FITS file, for TestTablesFloat32MatchCfitsio.
 *
 *   cf_tfloat <file.fits> > golden_tfloat/<file.fits>.json
 *
 * For every BINARY_TBL / ASCII_TBL HDU and every B, I, J, K, E, D column
 * (fixed or P/Q; ASCII I, F, E, D), each row is read with
 * fits_read_colnull(TFLOAT) under the column's own TSCAL/TZERO, and each
 * element is written as the 8 hex digits of its float32 bits, or null when
 * cfitsio flags it null. Output is compact JSON (no whitespace):
 *
 *   {"file":"x.fits","hdus":[{"index":1,"columns":[{"n":1,"ttype":"K",
 *    "rows":[["5d800001"],...]}]}]}
 *
 * "index" is the 0-based HDU number. A read that fails is recorded as
 * "error":<cfitsio status> for the column (no rows), except NUM_OVERFLOW
 * (412), after which cfitsio has stored +-FLT_MAX for the values beyond
 * the float32 range: the column gets "overflow":true and its rows.
 * Nothing here depends on the Go code. Built by `make cf_tfloat.out`
 * (system libcfitsio 4.7.0).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

static void jstr(const char *s) {
  putchar('"');
  for (; *s; s++) {
    if (*s == '"' || *s == '\\') putchar('\\');
    putchar(*s);
  }
  putchar('"');
}

static int numeric(int tc) {
  switch (tc < 0 ? -tc : tc) {
  case TBYTE: case TSBYTE: case TSHORT: case TUSHORT: case TINT: case TUINT:
  case TLONG: case TULONG: case TLONGLONG: case TULONGLONG: case TFLOAT: case TDOUBLE:
    return 1;
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: cf_tfloat <file.fits>\n");
    return 2;
  }
  fitsfile *f;
  int st = 0, nhdu = 0;
  if (fits_open_file(&f, argv[1], READONLY, &st) || fits_get_num_hdus(f, &nhdu, &st)) {
    fits_report_error(stderr, st);
    return 1;
  }
  const char *base = strrchr(argv[1], '/');
  printf("{\"file\":");
  jstr(base ? base + 1 : argv[1]);
  printf(",\"hdus\":[");
  int firsth = 1;
  for (int h = 1; h <= nhdu; h++) {
    int ht = 0;
    st = 0;
    if (fits_movabs_hdu(f, h, &ht, &st)) {
      fits_report_error(stderr, st);
      return 1;
    }
    if (ht != BINARY_TBL && ht != ASCII_TBL) continue;
    LONGLONG nrows = 0;
    int ncols = 0;
    fits_get_num_rowsll(f, &nrows, &st);
    fits_get_num_cols(f, &ncols, &st);
    if (st) {
      fits_report_error(stderr, st);
      return 1;
    }
    printf("%s{\"index\":%d,\"columns\":[", firsth ? "" : ",", h - 1);
    firsth = 0;
    int firstc = 1;
    for (int c = 1; c <= ncols; c++) {
      int tc = 0;
      LONGLONG rep = 0, width = 0;
      st = 0;
      if (fits_get_coltypell(f, c, &tc, &rep, &width, &st) || !numeric(tc)) continue;
      char key[FLEN_KEYWORD], ttype[FLEN_VALUE] = "";
      snprintf(key, sizeof key, "TTYPE%d", c);
      int s = 0;
      fits_read_key(f, TSTRING, key, ttype, NULL, &s);
      printf("%s{\"n\":%d,\"ttype\":", firstc ? "" : ",", c);
      firstc = 0;
      jstr(ttype);
      /* every row is read before anything is printed, so a failed read
       * leaves only the error */
      float **vals = calloc((size_t)nrows + 1, sizeof *vals);
      char **nuls = calloc((size_t)nrows + 1, sizeof *nuls);
      LONGLONG *ns = calloc((size_t)nrows + 1, sizeof *ns);
      int err = 0, over = 0;
      for (LONGLONG r = 1; r <= nrows && !err; r++) {
        LONGLONG n = rep;
        if (tc < 0) {
          LONGLONG off = 0;
          if (fits_read_descriptll(f, c, r, &n, &off, &st)) { err = st; break; }
        }
        ns[r] = n;
        vals[r] = malloc((size_t)(n > 0 ? n : 1) * sizeof(float));
        nuls[r] = malloc((size_t)(n > 0 ? n : 1));
        int anyn = 0;
        if (n > 0 && fits_read_colnull(f, TFLOAT, c, r, 1, n, vals[r], nuls[r], &anyn, &st)) {
          if (st == NUM_OVERFLOW) { /* values past FLT_MAX are set to +-FLT_MAX */
            over = 1;
            st = 0;
          } else {
            err = st;
          }
        }
      }
      if (err) {
        printf(",\"error\":%d}", err);
      } else {
        if (over) printf(",\"overflow\":true");
        printf(",\"rows\":[");
        for (LONGLONG r = 1; r <= nrows; r++) {
          printf("%s[", r > 1 ? "," : "");
          for (LONGLONG k = 0; k < ns[r]; k++) {
            if (k) putchar(',');
            if (nuls[r][k]) {
              printf("null");
            } else {
              unsigned int b;
              memcpy(&b, &vals[r][k], 4);
              printf("\"%08x\"", b);
            }
          }
          putchar(']');
        }
        printf("]}");
      }
      for (LONGLONG r = 1; r <= nrows; r++) { free(vals[r]); free(nuls[r]); }
      free(vals); free(nuls); free(ns);
    }
    printf("]}");
  }
  printf("]}\n");
  st = 0;
  fits_close_file(f, &st);
  return 0;
}
