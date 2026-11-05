/*
 * hdrref -- cfitsio truth harness for the Go header package tests.
 *
 *   hdrref write <out.fits>     cfitsio WRITES a reference header sample
 *                               (long strings, HIERARCH long strings, long
 *                               COMMENT/HISTORY, numbers) through its public
 *                               fits_write_* API.
 *   hdrref longstr <file.fits>  JSON, per HDU, for every card whose value
 *                               cfitsio types as a string ('C'): the value
 *                               and the COMMENT returned by
 *                               fits_read_key_longstr (CONTINUE resolved,
 *                               continuation comments joined by cfitsio);
 *                               for every integer-typed card ('I') whose
 *                               value is an unsigned decimal in uint64
 *                               range: the value read with
 *                               fits_read_key(TULONGLONG); for every
 *                               complex card ('X'): fits_read_key
 *                               (TDBLCOMPLEX) as [re, im].
 *   hdrref testkey <name>...    JSON: fits_test_keyword status per name
 *                               (0 = legal 8-char FITS keyword).
 *   hdrref readkey <file.fits> <name>...
 *                               JSON: for each name, as given (any case),
 *                               fits_read_keyword on the primary header with
 *                               the position reset to its start (so the
 *                               search finds the FIRST matching card): the
 *                               value text, comment and status (202 =
 *                               KEY_NO_EXIST). Also the number of cards
 *                               whose keyword name equals the upper-cased
 *                               name.
 *   hdrref delkey <in.fits> <out.fits> <name>
 *                               copy in to out and fits_delete_key(name) on
 *                               the primary header (position reset first).
 *   hdrref writekeys <out.fits> cfitsio WRITES integer and logical keywords
 *                               of every printed width through
 *                               fits_write_key: Pkk = 10^kk - 1 and
 *                               Nkk = -(10^kk - 1) for kk = 01..18, IMAX /
 *                               IMIN = int64 max / min, NBIG = -10^18,
 *                               UMAX = uint64 max, LT / LF = T / F. Keys
 *                               with an odd kk (and IMAX, UMAX, LT) carry
 *                               the comment "c"; the others none.
 *                               Last, "lowkey" = 5 / c: a lower-case name
 *                               given to fits_write_key.
 *
 * Everything here comes from the system libcfitsio (4.7.0 for the
 * committed goldens; see ../PROVENANCE.md). Nothing uses the Go library.
 * Card-level detail (raw card text, per-card value/comment, typed reads)
 * comes from the shared harness testdata/cfref ("cfref dump"); this
 * program only adds what cfref does not record.
 *
 * Strings are JSON-escaped; bytes >= 0x7F are emitted as \u00XX meaning
 * "the byte 0xXX" (same convention as cfref).
 */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

static void jstr(const char *s) {
  putchar('"');
  for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
    if (*p == '"') fputs("\\\"", stdout);
    else if (*p == '\\') fputs("\\\\", stdout);
    else if (*p < 0x20 || *p >= 0x7f) printf("\\u%04x", *p);
    else putchar(*p);
  }
  putchar('"');
}

static void jerr(int st) {
  char txt[FLEN_STATUS];
  fits_get_errstatus(st, txt);
  printf("\"status %d: %s\"", st, txt);
  fits_clear_errmsg();
}

static void die(const char *what, int st) {
  char txt[FLEN_STATUS];
  fits_get_errstatus(st, txt);
  fprintf(stderr, "hdrref: %s: status %d: %s\n", what, st, txt);
  exit(1);
}

/* Repeat pattern p until n characters. */
static char *rep(const char *p, int n) {
  char *s = malloc(n + 1);
  int k = strlen(p);
  for (int i = 0; i < n; i++) s[i] = p[i % k];
  s[n] = 0;
  return s;
}

static int do_write(const char *out) {
  fitsfile *f;
  int st = 0;
  char path[4096];
  snprintf(path, sizeof path, "!%s", out); /* '!' = overwrite */
  if (fits_create_file(&f, path, &st)) die("create", st);
  if (fits_create_img(f, BYTE_IMG, 0, NULL, &st)) die("create_img", st);

  /* Long strings (CONTINUE convention). The quote characters exercise
   * escaping across card boundaries. */
  char *s200 = rep("0123456789abcdef'ghij", 200);
  char *c120 = rep("comment text ", 120);
  if (fits_write_key_longstr(f, "LONGSTR", s200, c120, &st)) die("LONGSTR", st);
  char *s68 = rep("x", 68), *s69 = rep("y", 69);
  if (fits_write_key_longstr(f, "STR68", s68, "68 chars", &st)) die("STR68", st);
  if (fits_write_key_longstr(f, "STR69", s69, "69 chars", &st)) die("STR69", st);
  if (fits_write_key_longstr(f, "LEADSP", "   leading blanks", "kept", &st)) die("LEADSP", st);
  /* HIERARCH long string: cfitsio prepends HIERARCH to a long name. */
  char *h150 = rep("hierarch-value ", 150);
  if (fits_write_key_longstr(f, "ESO OBS LONGNAME", h150, "hierarch comment", &st)) die("HIERARCH long", st);
  if (fits_write_key_str(f, "ESO DET NAME", "short", "hierarch short", &st)) die("HIERARCH short", st);

  /* Long commentary text: cfitsio splits into 72-character cards. */
  char *hist = rep("history words ", 200);
  if (fits_write_history(f, hist, &st)) die("history", st);
  char *com = rep("comment words ", 150);
  if (fits_write_comment(f, com, &st)) die("comment", st);

  /* Numbers. */
  double d;
  d = 120.0;   if (fits_write_key(f, TDOUBLE, "FLT120", &d, "integral double", &st)) die("FLT120", st);
  d = 1e-7;    if (fits_write_key(f, TDOUBLE, "FSMALL", &d, "small", &st)) die("FSMALL", st);
  d = 1.5e300; if (fits_write_key(f, TDOUBLE, "FBIG", &d, "big", &st)) die("FBIG", st);
  d = 0.1;     if (fits_write_key(f, TDOUBLE, "FTENTH", &d, "0.1", &st)) die("FTENTH", st);
  ULONGLONG u = 18446744073709551615ULL;
  if (fits_write_key(f, TULONGLONG, "U64MAX", &u, "2**64-1", &st)) die("U64MAX", st);
  LONGLONG ll = -9223372036854775807LL - 1;
  if (fits_write_key(f, TLONGLONG, "I64MIN", &ll, "int64 min", &st)) die("I64MIN", st);
  int lg = 1;
  if (fits_write_key(f, TLOGICAL, "LOGT", &lg, "true", &st)) die("LOGT", st);
  double cx[2] = {1.5, -2.25};
  if (fits_write_key_dblcmp(f, "CPX", cx, -15, "complex", &st)) die("CPX", st);

  if (fits_close_file(f, &st)) die("close", st);
  free(s200); free(c120); free(s68); free(s69); free(h150); free(hist); free(com);
  return 0;
}

/* Only unsigned decimals in [0, 2^64-1] are read as TULONGLONG: cfitsio
 * 4.7.0's ffc2ujj overflows a message buffer (strncat, glibc fortify
 * abort) when reporting an out-of-range value such as a 24-digit integer,
 * and it wraps negative values. Those cards are covered by cfref's
 * TLONGLONG/TDOUBLE reads instead. */
static int u64_in_range(const char *v) {
  const char *max = "18446744073709551615";
  if (*v == '+') v++;
  size_t n = strlen(v);
  if (n == 0 || n > 20) return 0;
  for (size_t i = 0; i < n; i++) if (v[i] < '0' || v[i] > '9') return 0;
  while (n > 1 && *v == '0') { v++; n--; }
  if (n < 20) return 1;
  return strcmp(v, max) <= 0;
}

static int do_longstr(const char *in) {
  fitsfile *f;
  int st = 0, nhdu = 0;
  if (fits_open_file(&f, in, READONLY, &st)) die("open", st);
  fits_get_num_hdus(f, &nhdu, &st);
  const char *base = strrchr(in, '/');
  base = base ? base + 1 : in;
  printf("{\n \"file\": "); jstr(base);
  printf(",\n \"cfitsio\": \"%d.%d.%d\",\n \"hdus\": [", CFITSIO_MAJOR, CFITSIO_MINOR, CFITSIO_MICRO);
  for (int h = 1; h <= nhdu; h++) {
    int ht = 0;
    st = 0;
    fits_movabs_hdu(f, h, &ht, &st);
    int nkeys = 0, more = 0;
    fits_get_hdrspace(f, &nkeys, &more, &st);
    printf("%s\n  {\"index\": %d, \"keys\": [", h > 1 ? "," : "", h - 1);
    int first = 1;
    for (int i = 1; i <= nkeys; i++) {
      char card[FLEN_CARD], name[FLEN_KEYWORD], value[FLEN_VALUE], comment[FLEN_COMMENT];
      int len = 0, s = 0;
      char dt = 0;
      if (fits_read_record(f, i, card, &s)) continue;
      s = 0; fits_get_keyname(card, name, &len, &s);
      s = 0;
      if (fits_parse_value(card, value, comment, &s) || !value[0]) continue;
      s = 0;
      if (fits_get_keytype(value, &dt, &s) || (dt != 'C' && dt != 'I' && dt != 'X')) continue;
      /* Position before card i so a by-name read finds card i itself. */
      s = 0; fits_read_record(f, i - 1, card, &s);
      printf("%s\n   {\"n\": %d, \"key\": ", first ? "" : ",", i);
      jstr(name);
      first = 0;
      s = 0;
      if (dt == 'C') {
        char *ls = NULL, cm[FLEN_COMMENT] = "";
        if (fits_read_key_longstr(f, name, &ls, cm, &s) == 0) {
          printf(", \"str\": "); jstr(ls);
          printf(", \"comment\": "); jstr(cm);
        } else { printf(", \"error\": "); jerr(s); }
        if (ls) { int s2 = 0; fits_free_memory(ls, &s2); }
      } else if (dt == 'X') {
        double cx[2];
        if (fits_read_key(f, TDBLCOMPLEX, name, cx, NULL, &s) == 0) printf(", \"cplx\": [%.17g, %.17g]", cx[0], cx[1]);
        else { printf(", \"error\": "); jerr(s); }
      } else if (!u64_in_range(value)) {
        printf(", \"skipped\": \"not an unsigned decimal in uint64 range\"");
      } else {
        ULONGLONG u = 0;
        if (fits_read_key(f, TULONGLONG, name, &u, NULL, &s) == 0) printf(", \"u64\": %llu", u);
        else { printf(", \"error\": "); jerr(s); }
      }
      printf("}");
    }
    printf("\n  ]}");
  }
  printf("\n ]\n}\n");
  st = 0;
  fits_close_file(f, &st);
  return 0;
}

static int do_testkey(int n, char **names) {
  printf("[");
  for (int i = 0; i < n; i++) {
    int s = 0;
    fits_test_keyword(names[i], &s);
    printf("%s\n {\"key\": ", i ? "," : "");
    jstr(names[i]);
    printf(", \"status\": %d}", s);
    fits_clear_errmsg();
  }
  printf("\n]\n");
  return 0;
}

static int do_readkey(const char *path, int n, char **names) {
  fitsfile *f;
  int st = 0, nkeys = 0;
  if (fits_open_file(&f, path, READONLY, &st)) { fits_report_error(stderr, st); return 1; }
  fits_get_hdrspace(f, &nkeys, NULL, &st);
  printf("[");
  for (int i = 0; i < n; i++) {
    char value[FLEN_VALUE] = "", comment[FLEN_COMMENT] = "", card[FLEN_CARD], up[FLEN_KEYWORD];
    int s = 0, count = 0;
    size_t k;
    for (k = 0; names[i][k] && k < FLEN_KEYWORD - 1; k++) up[k] = (char)toupper((unsigned char)names[i][k]);
    up[k] = 0;
    for (int j = 1; j <= nkeys; j++) {
      char nm[FLEN_KEYWORD];
      int len = 0, s2 = 0;
      if (fits_read_record(f, j, card, &s2)) continue;
      fits_get_keyname(card, nm, &len, &s2);
      if (!strcmp(nm, up)) count++;
    }
    s = 0; fits_read_record(f, 0, card, &s); /* back to the start of the header */
    s = 0;
    fits_read_keyword(f, names[i], value, comment, &s);
    printf("%s\n {\"key\": ", i ? "," : "");
    jstr(names[i]);
    printf(", \"status\": %d, \"count\": %d", s, count);
    if (!s) { printf(", \"value\": "); jstr(value); printf(", \"comment\": "); jstr(comment); }
    printf("}");
    fits_clear_errmsg();
  }
  printf("\n]\n");
  st = 0;
  fits_close_file(f, &st);
  return 0;
}

static int do_delkey(const char *in, const char *out, const char *name) {
  fitsfile *f;
  int st = 0;
  char cmd[4096], card[FLEN_CARD];
  snprintf(cmd, sizeof cmd, "cp '%s' '%s'", in, out);
  if (system(cmd) != 0) return 1;
  if (fits_open_file(&f, out, READWRITE, &st)) { fits_report_error(stderr, st); return 1; }
  fits_read_record(f, 0, card, &st);
  fits_delete_key(f, name, &st);
  fits_close_file(f, &st);
  if (st) { fits_report_error(stderr, st); return 1; }
  return 0;
}

static void wkey(fitsfile *f, int type, const char *name, void *v, int withComment, int *st) {
  if (fits_write_key(f, type, name, v, withComment ? "c" : NULL, st)) die(name, *st);
}

static int do_writekeys(const char *out) {
  fitsfile *f;
  int st = 0;
  char path[4096], name[FLEN_KEYWORD];
  snprintf(path, sizeof path, "!%s", out);
  if (fits_create_file(&f, path, &st)) die("create", st);
  if (fits_create_img(f, BYTE_IMG, 0, NULL, &st)) die("create_img", st);
  LONGLONG p = 0;
  for (int k = 1; k <= 18; k++) {
    p = p * 10 + 9;
    LONGLONG n = -p;
    snprintf(name, sizeof name, "P%02d", k); wkey(f, TLONGLONG, name, &p, k % 2, &st);
    snprintf(name, sizeof name, "N%02d", k); wkey(f, TLONGLONG, name, &n, k % 2, &st);
  }
  LONGLONG imax = 9223372036854775807LL, imin = -9223372036854775807LL - 1, nbig = -1000000000000000000LL;
  ULONGLONG umax = 18446744073709551615ULL;
  int lt = 1, lf = 0;
  wkey(f, TLONGLONG, "IMAX", &imax, 1, &st);
  wkey(f, TLONGLONG, "IMIN", &imin, 0, &st);
  wkey(f, TLONGLONG, "NBIG", &nbig, 0, &st);
  wkey(f, TULONGLONG, "UMAX", &umax, 1, &st);
  wkey(f, TLOGICAL, "LT", &lt, 1, &st);
  wkey(f, TLOGICAL, "LF", &lf, 0, &st);
  /* A lower-case keyword name: cfitsio writes it upper-cased. */
  LONGLONG five = 5;
  wkey(f, TLONGLONG, "lowkey", &five, 1, &st);
  if (fits_close_file(f, &st)) die("close", st);
  return 0;
}

int main(int argc, char **argv) {
  if (argc == 3 && strcmp(argv[1], "writekeys") == 0) return do_writekeys(argv[2]);
  if (argc >= 3 && strcmp(argv[1], "write") == 0) return do_write(argv[2]);
  if (argc == 3 && strcmp(argv[1], "longstr") == 0) return do_longstr(argv[2]);
  if (argc >= 2 && strcmp(argv[1], "testkey") == 0) return do_testkey(argc - 2, argv + 2);
  if (argc >= 3 && strcmp(argv[1], "readkey") == 0) return do_readkey(argv[2], argc - 3, argv + 3);
  if (argc == 5 && strcmp(argv[1], "delkey") == 0) return do_delkey(argv[2], argv[3], argv[4]);
  fprintf(stderr, "usage: hdrref write <out.fits> | longstr <file.fits> | testkey <name>... | readkey <file.fits> <name>... | delkey <in> <out> <name> | writekeys <out.fits>\n");
  return 2;
}
