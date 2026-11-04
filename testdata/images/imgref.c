/*
 * imgref -- cfitsio truth harness for the image / checksum tests of the Go
 * fits library. It is linked only against the system libcfitsio; nothing
 * here uses the Go library.
 *
 *   imgref write <out.fits> <imgtype> <datatype> <bscale|-> <bzero|-> <blank|-> <naxes> <v1> <v2> ...
 *       Create a one-HDU image with fits_create_img(<imgtype>) (8, 10 =
 *       SBYTE_IMG, 16, 20 = USHORT_IMG, 32, 40 = ULONG_IMG, 64, 80 =
 *       ULONGLONG_IMG, -32, -64; cfitsio itself writes the BZERO/BSCALE
 *       cards of the unsigned / signed-byte conventions). When bscale/bzero
 *       are given they are written as BSCALE/BZERO cards and activated with
 *       fits_set_bscale; blank is written as BLANK. The values (decimal
 *       text, "nan"/"inf"/"-inf" for floats, or "bits:0x<hex>" for the exact
 *       IEEE bit pattern of a TFLOAT/TDOUBLE value, e.g. a signaling NaN
 *       that strtof cannot produce) are converted with strtoll /
 *       strtoull / strtod into the C type of <datatype> (a cfitsio type
 *       code: TBYTE 11, TSBYTE 12, TUSHORT 20, TSHORT 21, TUINT 30, TINT 31,
 *       TLONGLONG 81, TULONGLONG 80, TFLOAT 42, TDOUBLE 82) and written with
 *       fits_write_img. <naxes> is comma separated, NAXIS1 first.
 *       Prints {"status": N, "overflow": true|false} on stdout. cfitsio's
 *       NUM_OVERFLOW (412) means it clamped at least one value; the file is
 *       still written and the status is reported, not hidden. When
 *       fits_set_bscale refuses the scaling (ZERO_SCALE 322 for BSCALE 0)
 *       it prints {"status": N, "refused": "set_bscale"} instead, exits 0
 *       and leaves no file.
 *
 *   imgref overwrite <in.fits> <out.fits> <datatype> <v1> <v2> ...
 *       Copy <in> to <out>, open the copy READWRITE and write the values
 *       over the primary image with fits_write_img(<datatype>) (codes as
 *       for write): cfitsio scales them with the BSCALE/BZERO/BLANK it read
 *       from the existing header (ffgphd), the truth for OverwritePixels.
 *       Prints {"status": N} (412 NUM_OVERFLOW: clamped; the file is still
 *       written).
 *
 *   imgref chksum <file>
 *       JSON per HDU: fits_get_chksum datasum / hdusum (computed by cfitsio
 *       from the bytes on disk, independent of any CHECKSUM card) and
 *       fits_verify_chksum status (1 ok, 0 absent, -1 wrong). If cfitsio
 *       cannot compute the sums the HDU entry carries "error" instead; if
 *       fits_verify_chksum itself fails (a keyword with no value: 204) it
 *       carries the sums and "verify_error".
 *
 *   imgref addchk <in.fits> <out.fits>
 *       Copy <in> to <out> and add DATASUM/CHECKSUM to every HDU exactly as
 *       cfitsio's fits_write_chksum computes them, but with fixed comments
 *       (fits_write_chksum puts the current time in the CHECKSUM comment,
 *       which would make the committed sample non-reproducible).
 *
 *   imgref esum <u32> <complement 0|1>
 *       fits_encode_chksum of the value: the 16-character ASCII form.
 *
 *   imgref pix <file>
 *       Pixel reads the shared cfref dump does not make, for every image HDU:
 *       "f32": fits_read_imgnull(TFLOAT), the truth for float32 targets
 *       (cfitsio converts integers to float in one rounding, which differs
 *       from rounding a double for BITPIX 64); each float32 is written as
 *       its exact double ("%.17g"), null for cfitsio's null flag,
 *       "NaN"/"Inf"/"-Inf".
 *       "u64": fits_read_imgnull(TULONGLONG) for BITPIX 64 images (cfitsio
 *       4.7.0's fits_get_img_equivtype reports LONGLONG_IMG for the uint64
 *       BZERO = 2^63 convention, so cfref dumps only an overflowing TLONGLONG
 *       read for those).
 *       "f64nn": for integer BITPIX images that have null pixels, the
 *       TDOUBLE read without null checking (fits_read_img with nulval 0):
 *       the scaled stored value of every pixel, BLANK pixels included.
 *       Arrays are null for non-image or zero-pixel HDUs;
 *       "f32_error" / "u64_error" / "f64nn_error" carry cfitsio's status
 *       when a read fails; "f32_overflow": true marks a TFLOAT read where
 *       cfitsio reported NUM_OVERFLOW and clamped the out-of-range values
 *       to +-FLT_MAX. An array of more than DIGEST_MIN (1024) values is
 *       written as an exact sha256 digest plus samples, in the encoding of
 *       ../cfref/cfref.c (jdigest; the code below is a verbatim copy).
 *
 *   imgref open <file>
 *       What cfitsio makes of a file the Go library refuses: the status of
 *       fits_open_diskfile ("open_status", 0 = opened), fits_get_num_hdus
 *       ("nhdu"), for each HDU the status of a whole-image TDOUBLE
 *       fits_read_img ("read_status", 0 = read; -1 = not an image), and
 *       the status of fits_movabs_hdu to HDU nhdu+1 ("next_status": 107
 *       END_OF_FILE means cfitsio sees no further HDU). Always exits 0.
 *
 *   imgref strtod
 *       For each line of stdin (the line without its newline): the C
 *       library's strtod as cfitsio's ffc2dd uses it ("<hex bits of the
 *       double, every NaN as 7ff8000000000000> <bytes consumed> <1 if errno
 *       is ERANGE, else 0>") and strtoll base 10 as ffc2jj uses it ("<value>
 *       <bytes consumed> <ERANGE 0|1>"), tab separated, one line out per
 *       line in. The truth for the library's cStrtod / cStrtoll.
 *
 *   imgref c2
 *       For each line of stdin: cfitsio's own value converters, called
 *       directly (fitsio2.h): ffc2dd ("<status> <hex bits of the double>",
 *       bits 0 when status != 0) and ffc2jj ("<status> <value>", value 0
 *       when status != 0), tab separated. The truth for cfC2DD / cfC2JJ.
 *
 *   imgref subset <file> <hdu0> <i64|u64|f64|f32> <fpix> <lpix> <inc>
 *       fits_read_subset of HDU <hdu0> (0-based) with 1-based inclusive
 *       fpixel/lpixel and inc (comma separated, NAXIS1 first), read as
 *       TLONGLONG / TULONGLONG / TDOUBLE / TFLOAT (%.9g). JSON array on stdout; non-finite
 *       doubles as "NaN"/"Inf"/"-Inf".
 */
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>
#include <fitsio2.h> /* ffc2dd, ffc2jj */

static int die(const char *what, int st) {
  char txt[FLEN_STATUS];
  fits_get_errstatus(st, txt);
  fprintf(stderr, "imgref: %s: status %d: %s\n", what, st, txt);
  fits_report_error(stderr, st);
  return 1;
}

static void jdbl(double d) {
  if (isnan(d)) fputs("\"NaN\"", stdout);
  else if (isinf(d)) fputs(d > 0 ? "\"Inf\"" : "\"-Inf\"", stdout);
  else printf("%.17g", d);
}

static int parse_axes(const char *s, long long *ax, int max) {
  int n = 0;
  const char *p = s;
  while (*p && n < max) {
    char *end;
    ax[n++] = strtoll(p, &end, 10);
    if (end == p) return -1;
    p = end;
    if (*p == ',') p++;
  }
  return n;
}

static size_t typesize(int dt) {
  switch (dt) {
  case TBYTE: case TSBYTE: return 1;
  case TUSHORT: case TSHORT: return 2;
  case TUINT: case TINT: return sizeof(int);
  case TLONGLONG: case TULONGLONG: return 8;
  case TFLOAT: return 4;
  case TDOUBLE: return 8;
  }
  return 0;
}

static int store(int dt, void *buf, long long i, const char *txt) {
  char *end;
  errno = 0;
  if (!strncmp(txt, "bits:", 5) && (dt == TFLOAT || dt == TDOUBLE)) {
    unsigned long long u = strtoull(txt + 5, &end, 16);
    if (*end != 0 || errno == ERANGE) return -1;
    if (dt == TFLOAT) {
      unsigned int u32 = (unsigned int)u;
      if (u >> 32) return -1;
      memcpy((float *)buf + i, &u32, 4);
    } else {
      memcpy((double *)buf + i, &u, 8);
    }
    return 0;
  }
  switch (dt) {
  case TBYTE: ((unsigned char *)buf)[i] = (unsigned char)strtoull(txt, &end, 10); break;
  case TSBYTE: ((signed char *)buf)[i] = (signed char)strtoll(txt, &end, 10); break;
  case TUSHORT: ((unsigned short *)buf)[i] = (unsigned short)strtoull(txt, &end, 10); break;
  case TSHORT: ((short *)buf)[i] = (short)strtoll(txt, &end, 10); break;
  case TUINT: ((unsigned int *)buf)[i] = (unsigned int)strtoull(txt, &end, 10); break;
  case TINT: ((int *)buf)[i] = (int)strtoll(txt, &end, 10); break;
  case TLONGLONG: ((long long *)buf)[i] = strtoll(txt, &end, 10); break;
  case TULONGLONG: ((unsigned long long *)buf)[i] = strtoull(txt, &end, 10); break;
  case TFLOAT: ((float *)buf)[i] = strtof(txt, &end); break;
  case TDOUBLE: ((double *)buf)[i] = strtod(txt, &end); break;
  default: return -1;
  }
  if (*end != 0 || errno == ERANGE) {
    /* strtof/strtod set ERANGE for denormals; accept those, reject the rest */
    if (*end != 0) return -1;
    if (dt != TFLOAT && dt != TDOUBLE) return -1;
  }
  return 0;
}

static int cmd_write(int argc, char **argv) {
  /* argv: write out imgtype datatype bscale bzero blank naxes values... */
  if (argc < 9) { fprintf(stderr, "write: too few arguments\n"); return 2; }
  const char *out = argv[2];
  int imgtype = atoi(argv[3]);
  int dt = atoi(argv[4]);
  const char *sbscale = argv[5], *sbzero = argv[6], *sblank = argv[7];
  long long ax[9];
  int naxis = parse_axes(argv[8], ax, 9);
  if (naxis < 1) { fprintf(stderr, "write: bad naxes %s\n", argv[8]); return 2; }
  long long n = 1;
  for (int i = 0; i < naxis; i++) n *= ax[i];
  if (argc - 9 != n) { fprintf(stderr, "write: %d values for %lld pixels\n", argc - 9, n); return 2; }
  size_t sz = typesize(dt);
  if (!sz) { fprintf(stderr, "write: bad datatype %d\n", dt); return 2; }
  void *buf = calloc(n, sz);
  for (long long i = 0; i < n; i++)
    if (store(dt, buf, i, argv[9 + i])) { fprintf(stderr, "write: bad value %s\n", argv[9 + i]); return 2; }

  fitsfile *f;
  int st = 0;
  remove(out);
  if (fits_create_diskfile(&f, out, &st)) return die("create", st);
  if (fits_create_imgll(f, imgtype, naxis, ax, &st)) return die("create_img", st);
  double bscale = 1.0, bzero = 0.0;
  int scaled = 0;
  if (strcmp(sbscale, "-")) { bscale = strtod(sbscale, NULL); scaled = 1; }
  if (strcmp(sbzero, "-")) { bzero = strtod(sbzero, NULL); scaled = 1; }
  if (scaled) {
    if (fits_update_key_dbl(f, "BSCALE", bscale, -17, "physical = BZERO + BSCALE * raw", &st)) return die("BSCALE", st);
    if (fits_update_key_dbl(f, "BZERO", bzero, -17, "physical = BZERO + BSCALE * raw", &st)) return die("BZERO", st);
    if (fits_set_bscale(f, bscale, bzero, &st)) {
      int rst = st;
      st = 0;
      fits_delete_file(f, &st);
      fits_clear_errmsg();
      printf("{\"status\": %d, \"refused\": \"set_bscale\"}\n", rst);
      free(buf);
      return 0;
    }
  }
  if (strcmp(sblank, "-")) {
    long long blank = strtoll(sblank, NULL, 10);
    if (fits_update_key(f, TLONGLONG, "BLANK", &blank, "raw value of undefined pixels", &st)) return die("BLANK", st);
  }
  int wst = 0;
  fits_write_img(f, dt, 1, n, buf, &wst);
  int overflow = wst == NUM_OVERFLOW;
  if (wst && !overflow) return die("write_img", wst);
  fits_clear_errmsg();
  st = 0;
  if (fits_close_file(f, &st)) return die("close", st);
  printf("{\"status\": %d, \"overflow\": %s}\n", wst, overflow ? "true" : "false");
  free(buf);
  return 0;
}

static int cmd_overwrite(int argc, char **argv) {
  /* argv: overwrite in out datatype values... */
  if (argc < 5) { fprintf(stderr, "overwrite: too few arguments\n"); return 2; }
  const char *in = argv[2], *out = argv[3];
  int dt = atoi(argv[4]);
  long long n = argc - 5;
  size_t sz = typesize(dt);
  if (!sz) { fprintf(stderr, "overwrite: bad datatype %d\n", dt); return 2; }
  void *buf = calloc(n, sz);
  for (long long i = 0; i < n; i++)
    if (store(dt, buf, i, argv[5 + i])) { fprintf(stderr, "overwrite: bad value %s\n", argv[5 + i]); return 2; }
  FILE *fi = fopen(in, "rb"), *fo = fopen(out, "wb");
  if (!fi || !fo) { fprintf(stderr, "overwrite: cannot copy %s to %s\n", in, out); return 2; }
  char cb[2880];
  size_t got;
  while ((got = fread(cb, 1, sizeof cb, fi)) > 0)
    if (fwrite(cb, 1, got, fo) != got) { fprintf(stderr, "overwrite: write %s\n", out); return 2; }
  fclose(fi);
  if (fclose(fo)) { fprintf(stderr, "overwrite: close %s\n", out); return 2; }
  fitsfile *f;
  int st = 0, naxis = 0;
  long long ax[9];
  if (fits_open_diskfile(&f, out, READWRITE, &st)) return die("open", st);
  if (fits_get_img_dim(f, &naxis, &st) || fits_get_img_sizell(f, 9, ax, &st)) return die("img size", st);
  long long npix = 1;
  for (int i = 0; i < naxis; i++) npix *= ax[i];
  if (npix != n) { fprintf(stderr, "overwrite: %lld values for %lld pixels\n", n, npix); return 2; }
  int wst = 0;
  fits_write_img(f, dt, 1, n, buf, &wst);
  if (wst && wst != NUM_OVERFLOW) return die("write_img", wst);
  fits_clear_errmsg();
  st = 0;
  if (fits_close_file(f, &st)) return die("close", st);
  printf("{\"status\": %d}\n", wst);
  free(buf);
  return 0;
}

static int cmd_chksum(const char *path) {
  fitsfile *f;
  int st = 0, nhdu = 0, type = 0;
  if (fits_open_diskfile(&f, path, READONLY, &st)) return die("open", st);
  fits_get_num_hdus(f, &nhdu, &st);
  const char *base = strrchr(path, '/');
  base = base ? base + 1 : path;
  printf("{\"file\": \"%s\", \"cfitsio\": \"%d.%d.%d\", \"hdus\": [", base, CFITSIO_MAJOR, CFITSIO_MINOR, CFITSIO_MICRO);
  for (int h = 1; h <= nhdu; h++) {
    st = 0;
    if (fits_movabs_hdu(f, h, &type, &st)) return die("movabs", st);
    unsigned long ds = 0, hs = 0;
    printf("%s\n {\"index\": %d", h > 1 ? "," : "", h - 1);
    if (fits_get_chksum(f, &ds, &hs, &st)) {
      /* recorded, not fatal: e.g. a file whose last data block is short */
      char txt[FLEN_STATUS];
      fits_get_errstatus(st, txt);
      fits_clear_errmsg();
      printf(", \"error\": \"get_chksum status %d: %s\"}", st, txt);
      continue;
    }
    int dok = 0, hok = 0;
    if (fits_verify_chksum(f, &dok, &hok, &st)) {
      /* recorded, not fatal: e.g. a DATASUM card with no value (204) */
      char txt[FLEN_STATUS];
      fits_get_errstatus(st, txt);
      fits_clear_errmsg();
      printf(", \"datasum\": %lu, \"hdusum\": %lu, \"verify_error\": \"verify_chksum status %d: %s\"}", ds, hs, st, txt);
      st = 0;
      continue;
    }
    printf(", \"datasum\": %lu, \"hdusum\": %lu, \"verify_datasum\": %d, \"verify_hdusum\": %d}",
           ds, hs, dok, hok);
  }
  printf("\n]}\n");
  st = 0;
  fits_close_file(f, &st);
  return 0;
}

static int cmd_addchk(const char *in, const char *out) {
  fitsfile *fi, *fo;
  int st = 0, nhdu = 0, type = 0;
  remove(out);
  if (fits_open_diskfile(&fi, in, READONLY, &st)) return die("open in", st);
  if (fits_create_diskfile(&fo, out, &st)) return die("create out", st);
  if (fits_copy_file(fi, fo, 1, 1, 1, &st)) return die("copy", st);
  fits_close_file(fi, &st);
  fits_get_num_hdus(fo, &nhdu, &st);
  for (int h = 1; h <= nhdu; h++) {
    if (fits_movabs_hdu(fo, h, &type, &st)) return die("movabs", st);
    unsigned long ds = 0, hs = 0;
    /* Same steps as cfitsio ffpcks: DATASUM first (data only), then a
     * zero CHECKSUM placeholder, then the complement of the HDU sum. */
    if (fits_get_chksum(fo, &ds, &hs, &st)) return die("get_chksum", st);
    char dsstr[32];
    snprintf(dsstr, sizeof dsstr, "%lu", ds);
    if (fits_update_key_str(fo, "DATASUM", dsstr, "data unit checksum", &st)) return die("DATASUM", st);
    if (fits_update_key_str(fo, "CHECKSUM", "0000000000000000", "HDU checksum", &st)) return die("CHECKSUM", st);
    /* re-scan the header (writes the moved END card) before summing, as
     * ffpcks does with ffrdef */
    if (fits_set_hdustruc(fo, &st)) return die("set_hdustruc", st);
    if (fits_get_chksum(fo, &ds, &hs, &st)) return die("get_chksum 2", st);
    char enc[17];
    fits_encode_chksum(hs, 1, enc);
    if (fits_update_key_str(fo, "CHECKSUM", enc, "HDU checksum", &st)) return die("CHECKSUM 2", st);
    if (fits_set_hdustruc(fo, &st)) return die("set_hdustruc 2", st);
    int dok = 0, hok = 0;
    if (fits_verify_chksum(fo, &dok, &hok, &st)) return die("verify", st);
    if (dok != 1 || hok != 1) { fprintf(stderr, "addchk: HDU %d verify %d/%d\n", h, dok, hok); return 1; }
  }
  if (fits_close_file(fo, &st)) return die("close", st);
  return 0;
}

/* ---- copied verbatim from ../cfref/cfref.c (jstr, SHA-256, jdigest, jarr) */
static void jstr_n(const char *s, size_t n) {
  putchar('"');
  for (size_t i = 0; i < n; i++) {
    unsigned char c = (unsigned char)s[i];
    if (c == '"') fputs("\\\"", stdout);
    else if (c == '\\') fputs("\\\\", stdout);
    else if (c < 0x20 || c >= 0x7f) printf("\\u%04x", c);
    else putchar(c);
  }
  putchar('"');
}
static void jstr(const char *s) { jstr_n(s, strlen(s)); }

/* ---- SHA-256 (FIPS 180-4), self-contained so the harness links only
 * against libcfitsio. Checked against Go's crypto/sha256 by every digest
 * comparison and against known vectors by TestCfrefSha256Selftest. */
typedef struct { unsigned int h[8]; unsigned char buf[64]; unsigned long long len; int nbuf; } sha256_ctx;
static const unsigned int sha_k[64] = {
  0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
  0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
  0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
  0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
  0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
  0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
  0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
  0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
#define ROR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))
static void sha_block(sha256_ctx *c, const unsigned char *p) {
  unsigned int w[64], a, b, cc, d, e, f, g, h;
  for (int i = 0; i < 16; i++)
    w[i] = (unsigned int)p[4 * i] << 24 | (unsigned int)p[4 * i + 1] << 16 | (unsigned int)p[4 * i + 2] << 8 | p[4 * i + 3];
  for (int i = 16; i < 64; i++) {
    unsigned int s0 = ROR(w[i - 15], 7) ^ ROR(w[i - 15], 18) ^ (w[i - 15] >> 3);
    unsigned int s1 = ROR(w[i - 2], 17) ^ ROR(w[i - 2], 19) ^ (w[i - 2] >> 10);
    w[i] = w[i - 16] + s0 + w[i - 7] + s1;
  }
  a = c->h[0]; b = c->h[1]; cc = c->h[2]; d = c->h[3]; e = c->h[4]; f = c->h[5]; g = c->h[6]; h = c->h[7];
  for (int i = 0; i < 64; i++) {
    unsigned int t1 = h + (ROR(e, 6) ^ ROR(e, 11) ^ ROR(e, 25)) + ((e & f) ^ (~e & g)) + sha_k[i] + w[i];
    unsigned int t2 = (ROR(a, 2) ^ ROR(a, 13) ^ ROR(a, 22)) + ((a & b) ^ (a & cc) ^ (b & cc));
    h = g; g = f; f = e; e = d + t1; d = cc; cc = b; b = a; a = t1 + t2;
  }
  c->h[0] += a; c->h[1] += b; c->h[2] += cc; c->h[3] += d; c->h[4] += e; c->h[5] += f; c->h[6] += g; c->h[7] += h;
}
static void sha_init(sha256_ctx *c) {
  static const unsigned int iv[8] = {0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
  memcpy(c->h, iv, sizeof iv); c->len = 0; c->nbuf = 0;
}
static void sha_update(sha256_ctx *c, const unsigned char *p, size_t n) {
  c->len += n;
  while (n--) { c->buf[c->nbuf++] = *p++; if (c->nbuf == 64) { sha_block(c, c->buf); c->nbuf = 0; } }
}
static void sha_hex(sha256_ctx *c, char out[65]) {
  unsigned long long bits = c->len * 8;
  unsigned char pad = 0x80, z = 0, l[8];
  sha_update(c, &pad, 1);
  while (c->nbuf != 56) sha_update(c, &z, 1);
  for (int i = 0; i < 8; i++) l[i] = (unsigned char)(bits >> (56 - 8 * i));
  sha_update(c, l, 8);
  for (int i = 0; i < 8; i++) sprintf(out + 8 * i, "%08x", c->h[i]);
}
static void sha_u64(sha256_ctx *c, unsigned long long v) {
  unsigned char b[8];
  for (int i = 0; i < 8; i++) b[i] = (unsigned char)(v >> (56 - 8 * i));
  sha_update(c, b, 8);
}
static unsigned long long dbl_bits(double d) {
  unsigned long long u;
  if (isnan(d)) return 0x7ff8000000000000ULL; /* one canonical NaN */
  memcpy(&u, &d, 8);
  return u;
}

/* ---- Value arrays: explicit, or a digest when large.
 *
 * An array of more than DIGEST_MIN elements is written as a digest object
 * in place of the JSON array:
 *   {"n": N, "sha256": "<hex>", "nulls": [[start,count],...],
 *    "sample": [[index,value],...]}
 * sha256 is taken over the canonical big-endian encoding of all N
 * elements in order; per element (null = cfitsio's null flag set):
 *   f64        flag byte (1 valid, 0 null) + 8-byte IEEE bits (every NaN
 *              as 0x7ff8000000000000; a null -> 8 zero bytes)
 *   i64 / u64  flag byte + 8-byte two's complement / unsigned
 *   str        flag byte + 4-byte length + the bytes (null -> length 0)
 *   log        flag byte + 1 byte 0/1 (null -> 0)
 *   bits       1 byte 0/1
 *   cplx       8-byte bits of re + 8-byte bits of im (NaN as above)
 * "nulls" lists every null as runs. "sample" holds explicit values (the
 * same JSON form as an explicit array element) at a deterministic index
 * set: the first and last 8, the image corners of the first and last
 * plane, the first minimum and first maximum (numeric kinds, NaN and
 * nulls skipped) and 64 strided indices k*(N-1)/63; ascending, each once.
 * internal/testutil verifies the hash AND every sample.
 *
 * DIGEST_MIN is 1024 so the 960-pixel synthetic compression corpus
 * (comp_fz_* / comp_fzun_* / comp_fzsrc_*) stays explicit: its tests
 * relate cfitsio's decode to the funpack twin and the source element by
 * element (fallback rows read unscaled, twin within 1 float32 ULP,
 * quantized floats whose TDOUBLE values the library rounds to float32),
 * which an exact digest cannot express. Larger arrays (real images, the
 * astropy fixtures, the iter_* tables) are checked through
 * testutil.Values (Check / Resolve). */
#define DIGEST_MIN 1024
enum { K_F64, K_I64, K_U64, K_STR, K_LOG, K_BIT, K_CPLX };

static int has_null(int kind) { return kind != K_BIT && kind != K_CPLX; }

static void jelem(int kind, const void *v, const char *nul, LONGLONG k) {
  if (nul && nul[k] && has_null(kind)) { fputs("null", stdout); return; }
  switch (kind) {
  case K_F64: jdbl(((const double *)v)[k]); break;
  case K_I64: printf("%lld", ((const LONGLONG *)v)[k]); break;
  case K_U64: printf("%llu", ((const ULONGLONG *)v)[k]); break;
  case K_STR: jstr(((char *const *)v)[k]); break;
  case K_LOG: case K_BIT: printf("%d", ((const char *)v)[k] ? 1 : 0); break;
  case K_CPLX:
    putchar('['); jdbl(((const double *)v)[2 * k]); putchar(','); jdbl(((const double *)v)[2 * k + 1]); putchar(']');
    break;
  }
}

static int cmp_ll(const void *a, const void *b) {
  LONGLONG x = *(const LONGLONG *)a, y = *(const LONGLONG *)b;
  return x < y ? -1 : x > y;
}

/* naxes: an image's axes (adds the plane corners to the sample), or NULL. */
static void jdigest(int kind, const void *v, const char *nul, LONGLONG n, int naxis, const LONGLONG *naxes) {
  sha256_ctx c;
  sha_init(&c);
  for (LONGLONG k = 0; k < n; k++) {
    int isnull = nul && nul[k] && has_null(kind);
    unsigned char flag = isnull ? 0 : 1;
    switch (kind) {
    case K_F64: sha_update(&c, &flag, 1); sha_u64(&c, isnull ? 0 : dbl_bits(((const double *)v)[k])); break;
    case K_I64: sha_update(&c, &flag, 1); sha_u64(&c, isnull ? 0 : (unsigned long long)((const LONGLONG *)v)[k]); break;
    case K_U64: sha_update(&c, &flag, 1); sha_u64(&c, isnull ? 0 : ((const ULONGLONG *)v)[k]); break;
    case K_STR: {
      const char *s = isnull ? "" : ((char *const *)v)[k];
      unsigned int len = (unsigned int)strlen(s);
      unsigned char lb[4] = {(unsigned char)(len >> 24), (unsigned char)(len >> 16), (unsigned char)(len >> 8), (unsigned char)len};
      sha_update(&c, &flag, 1); sha_update(&c, lb, 4); sha_update(&c, (const unsigned char *)s, len);
      break;
    }
    case K_LOG: {
      unsigned char b = isnull ? 0 : (((const char *)v)[k] ? 1 : 0);
      sha_update(&c, &flag, 1); sha_update(&c, &b, 1);
      break;
    }
    case K_BIT: { unsigned char b = ((const char *)v)[k] ? 1 : 0; sha_update(&c, &b, 1); break; }
    case K_CPLX:
      sha_u64(&c, dbl_bits(((const double *)v)[2 * k]));
      sha_u64(&c, dbl_bits(((const double *)v)[2 * k + 1]));
      break;
    }
  }
  char hex[65];
  sha_hex(&c, hex);
  printf("{\"n\": %lld, \"sha256\": \"%s\", \"nulls\": [", n, hex);
  int first = 1;
  if (nul && has_null(kind)) {
    for (LONGLONG k = 0; k < n;) {
      if (!nul[k]) { k++; continue; }
      LONGLONG s = k;
      while (k < n && nul[k]) k++;
      printf("%s[%lld,%lld]", first ? "" : ",", s, k - s);
      first = 0;
    }
  }
  printf("], \"sample\": [");
  LONGLONG idx[8 + 8 + 8 + 2 + 64];
  int ni = 0;
  for (LONGLONG k = 0; k < 8 && k < n; k++) idx[ni++] = k;
  for (LONGLONG k = n - 8; k < n; k++) if (k >= 0) idx[ni++] = k;
  if (naxes && naxis >= 2 && naxes[0] > 0 && naxes[1] > 0) {
    LONGLONG nx = naxes[0], ny = naxes[1], plane = nx * ny;
    LONGLONG bases[2] = {0, (n / plane - 1) * plane};
    for (int b = 0; b < 2; b++) {
      idx[ni++] = bases[b]; idx[ni++] = bases[b] + nx - 1;
      idx[ni++] = bases[b] + (ny - 1) * nx; idx[ni++] = bases[b] + plane - 1;
    }
  }
  if (kind == K_F64 || kind == K_I64 || kind == K_U64) {
    LONGLONG kmin = -1, kmax = -1;
    for (LONGLONG k = 0; k < n; k++) {
      if (nul && nul[k]) continue;
      if (kind == K_F64) {
        const double *x = v;
        if (isnan(x[k])) continue;
        if (kmin < 0 || x[k] < x[kmin]) kmin = k;
        if (kmax < 0 || x[k] > x[kmax]) kmax = k;
      } else if (kind == K_I64) {
        const LONGLONG *x = v;
        if (kmin < 0 || x[k] < x[kmin]) kmin = k;
        if (kmax < 0 || x[k] > x[kmax]) kmax = k;
      } else {
        const ULONGLONG *x = v;
        if (kmin < 0 || x[k] < x[kmin]) kmin = k;
        if (kmax < 0 || x[k] > x[kmax]) kmax = k;
      }
    }
    if (kmin >= 0) { idx[ni++] = kmin; idx[ni++] = kmax; }
  }
  for (LONGLONG k = 0; k < 64 && n > 0; k++) idx[ni++] = k * (n - 1) / 63;
  qsort(idx, ni, sizeof idx[0], cmp_ll);
  first = 1;
  for (int i = 0; i < ni; i++) {
    if (i && idx[i] == idx[i - 1]) continue;
    printf("%s[%lld,", first ? "" : ",", idx[i]);
    jelem(kind, v, nul, idx[i]);
    putchar(']');
    first = 0;
  }
  printf("]}");
}

/* Print n values: an explicit JSON array, or a digest past DIGEST_MIN. */
static int digest_ok(LONGLONG n) { return n > DIGEST_MIN; }
static void jarr(int kind, const void *v, const char *nul, LONGLONG n, int naxis, const LONGLONG *naxes) {
  if (digest_ok(n)) { jdigest(kind, v, nul, n, naxis, naxes); return; }
  putchar('[');
  for (LONGLONG k = 0; k < n; k++) {
    if (k) putchar(',');
    jelem(kind, v, nul, k);
  }
  putchar(']');
}
/* ---- end of the copy */

static void pix_error(const char *field, int st) {
  char txt[FLEN_STATUS];
  fits_get_errstatus(st, txt);
  fits_clear_errmsg();
  printf("null, \"%s_error\": \"status %d: %s\"", field, st, txt);
}

static int cmd_pix(const char *path) {
  fitsfile *f;
  int st = 0, nhdu = 0, type = 0;
  if (fits_open_diskfile(&f, path, READONLY, &st)) return die("open", st);
  fits_get_num_hdus(f, &nhdu, &st);
  printf("{\"cfitsio\": \"%d.%d.%d\", \"hdus\": [", CFITSIO_MAJOR, CFITSIO_MINOR, CFITSIO_MICRO);
  for (int h = 1; h <= nhdu; h++) {
    st = 0;
    if (fits_movabs_hdu(f, h, &type, &st)) return die("movabs", st);
    printf("%s\n {\"index\": %d, \"f32\": ", h > 1 ? "," : "", h - 1);
    int bitpix = 0, naxis = 0;
    LONGLONG ax[999];
    LONGLONG n = 0;
    if (type == IMAGE_HDU && fits_get_img_paramll(f, 999, &bitpix, &naxis, ax, &st) == 0 && naxis > 0) {
      n = 1;
      for (int i = 0; i < naxis; i++) n *= ax[i];
    }
    st = 0;
    if (n == 0) { printf("null, \"u64\": null}"); continue; }
    char *nul = malloc(n);
    int anyn = 0, anynull = 0;
    float *v = malloc(n * sizeof *v);
    int ovf = 0;
    if (fits_read_imgnull(f, TFLOAT, 1, n, v, nul, &anyn, &st) == NUM_OVERFLOW) {
      /* cfitsio clamps out-of-range values to +-FLT_MAX and reports 412;
       * keep the values and say so */
      ovf = 1;
      st = 0;
      fits_clear_errmsg();
    }
    if (st) pix_error("f32", st);
    else {
      double *d = malloc(n * sizeof *d);
      for (LONGLONG k = 0; k < n; k++) d[k] = (double)v[k];
      jarr(K_F64, d, nul, n, naxis, ax);
      free(d);
      for (LONGLONG k = 0; k < n; k++) anynull |= nul[k];
      if (ovf) printf(", \"f32_overflow\": true");
    }
    free(v);
    printf(", \"u64\": ");
    st = 0;
    if (bitpix != 64) printf("null");
    else {
      unsigned long long *u = malloc(n * sizeof *u);
      if (fits_read_imgnull(f, TULONGLONG, 1, n, u, nul, &anyn, &st)) pix_error("u64", st);
      else jarr(K_U64, u, nul, n, naxis, ax);
      free(u);
    }
    if (bitpix > 0 && anynull) {
      /* nulval 0 turns cfitsio's null checking off (getcold.c: nultyp 1,
       * nulval 0 -> nulcheck 0): every pixel's scaled stored value */
      double *d = malloc(n * sizeof *d), zero = 0;
      int any2 = 0;
      st = 0;
      printf(", \"f64nn\": ");
      if (fits_read_img(f, TDOUBLE, 1, n, &zero, d, &any2, &st)) pix_error("f64nn", st);
      else jarr(K_F64, d, NULL, n, naxis, ax);
      free(d);
    }
    printf("}");
    free(nul);
  }
  printf("\n]}\n");
  st = 0;
  fits_close_file(f, &st);
  return 0;
}

static int cmd_open(const char *path) {
  fitsfile *f;
  int st = 0, nhdu = 0, type = 0;
  printf("{\"cfitsio\": \"%d.%d.%d\", ", CFITSIO_MAJOR, CFITSIO_MINOR, CFITSIO_MICRO);
  if (fits_open_diskfile(&f, path, READONLY, &st)) {
    printf("\"open_status\": %d, \"nhdu\": 0, \"hdus\": [], \"next_status\": 0}\n", st);
    fits_clear_errmsg();
    return 0;
  }
  fits_get_num_hdus(f, &nhdu, &st);
  printf("\"open_status\": 0, \"nhdu\": %d, \"hdus\": [", nhdu);
  for (int h = 1; h <= nhdu; h++) {
    int rs = -1;
    st = 0;
    if (fits_movabs_hdu(f, h, &type, &st) == 0 && type == IMAGE_HDU) {
      int bitpix = 0, naxis = 0;
      LONGLONG ax[999], n = 1;
      if (fits_get_img_paramll(f, 999, &bitpix, &naxis, ax, &st) == 0) {
        for (int i = 0; i < naxis; i++) n *= ax[i];
        if (naxis == 0) n = 0;
        double *v = malloc((n ? n : 1) * sizeof *v);
        int anyn = 0;
        if (n > 0) fits_read_img(f, TDOUBLE, 1, n, NULL, v, &anyn, &st);
        free(v);
      }
      rs = st;
    } else if (st) {
      rs = st;
    }
    printf("%s{\"index\": %d, \"read_status\": %d}", h > 1 ? ", " : "", h - 1, rs);
    fits_clear_errmsg();
  }
  st = 0;
  fits_movabs_hdu(f, nhdu + 1, &type, &st);
  printf("], \"next_status\": %d}\n", st);
  fits_clear_errmsg();
  st = 0;
  fits_close_file(f, &st);
  return 0;
}

static int cmd_subset(int argc, char **argv) {
  if (argc != 8) { fprintf(stderr, "subset: wrong argument count\n"); return 2; }
  const char *path = argv[2];
  int hdu = atoi(argv[3]);
  const char *ty = argv[4];
  long long fp[9], lp[9], inc[9];
  int n1 = parse_axes(argv[5], fp, 9), n2 = parse_axes(argv[6], lp, 9), n3 = parse_axes(argv[7], inc, 9);
  if (n1 < 1 || n1 != n2 || n1 != n3) { fprintf(stderr, "subset: bad pixel lists\n"); return 2; }
  long fpl[9], lpl[9], incl[9];
  long long n = 1;
  for (int i = 0; i < n1; i++) {
    fpl[i] = fp[i]; lpl[i] = lp[i]; incl[i] = inc[i];
    n *= (lp[i] - fp[i]) / inc[i] + 1;
  }
  fitsfile *f;
  int st = 0, type = 0, anyn = 0;
  if (fits_open_diskfile(&f, path, READONLY, &st)) return die("open", st);
  if (fits_movabs_hdu(f, hdu + 1, &type, &st)) return die("movabs", st);
  putchar('[');
  if (!strcmp(ty, "i64")) {
    long long *v = malloc(n * sizeof *v);
    if (fits_read_subset(f, TLONGLONG, fpl, lpl, incl, NULL, v, &anyn, &st)) return die("read_subset", st);
    for (long long k = 0; k < n; k++) printf("%s%lld", k ? "," : "", v[k]);
    free(v);
  } else if (!strcmp(ty, "u64")) {
    unsigned long long *v = malloc(n * sizeof *v);
    if (fits_read_subset(f, TULONGLONG, fpl, lpl, incl, NULL, v, &anyn, &st)) return die("read_subset", st);
    for (long long k = 0; k < n; k++) printf("%s%llu", k ? "," : "", v[k]);
    free(v);
  } else if (!strcmp(ty, "f64")) {
    double *v = malloc(n * sizeof *v);
    if (fits_read_subset(f, TDOUBLE, fpl, lpl, incl, NULL, v, &anyn, &st)) return die("read_subset", st);
    for (long long k = 0; k < n; k++) { if (k) putchar(','); jdbl(v[k]); }
    free(v);
  } else if (!strcmp(ty, "f32")) {
    float *v = malloc(n * sizeof *v);
    if (fits_read_subset(f, TFLOAT, fpl, lpl, incl, NULL, v, &anyn, &st)) return die("read_subset", st);
    for (long long k = 0; k < n; k++) { if (k) putchar(','); if (isnan(v[k]) || isinf(v[k])) jdbl(v[k]); else printf("%.9g", (double)v[k]); }
    free(v);
  } else { fprintf(stderr, "subset: bad type %s\n", ty); return 2; }
  printf("]\n");
  fits_close_file(f, &st);
  return 0;
}

static int cmd_strtod(void) {
  char line[1024];
  while (fgets(line, sizeof line, stdin)) {
    line[strcspn(line, "\n")] = 0;
    char *e;
    errno = 0;
    double d = strtod(line, &e);
    int er = errno == ERANGE;
    unsigned long long u;
    memcpy(&u, &d, 8);
    if (isnan(d)) u = 0x7ff8000000000000ULL;
    printf("%016llx %ld %d\t", u, (long)(e - line), er);
    errno = 0;
    long long ll = strtoll(line, &e, 10);
    printf("%lld %ld %d\n", ll, (long)(e - line), errno == ERANGE);
  }
  return 0;
}

static int cmd_c2(void) {
  char line[1024];
  while (fgets(line, sizeof line, stdin)) {
    line[strcspn(line, "\n")] = 0;
    int st = 0;
    double d = 0;
    ffc2dd(line, &d, &st);
    unsigned long long u = 0;
    if (st == 0) memcpy(&u, &d, 8);
    printf("%d %016llx\t", st, u);
    st = 0;
    LONGLONG j = 0;
    ffc2jj(line, &j, &st);
    printf("%d %lld\n", st, st == 0 ? (long long)j : 0LL);
    fits_clear_errmsg();
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc >= 2 && !strcmp(argv[1], "write")) return cmd_write(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "overwrite")) return cmd_overwrite(argc, argv);
  if (argc == 3 && !strcmp(argv[1], "chksum")) return cmd_chksum(argv[2]);
  if (argc == 4 && !strcmp(argv[1], "addchk")) return cmd_addchk(argv[2], argv[3]);
  if (argc == 4 && !strcmp(argv[1], "esum")) {
    char enc[17];
    fits_encode_chksum(strtoul(argv[2], NULL, 10), atoi(argv[3]), enc);
    printf("%s\n", enc);
    return 0;
  }
  if (argc == 3 && !strcmp(argv[1], "pix")) return cmd_pix(argv[2]);
  if (argc >= 2 && !strcmp(argv[1], "subset")) return cmd_subset(argc, argv);
  if (argc == 3 && !strcmp(argv[1], "open")) return cmd_open(argv[2]);
  if (argc == 2 && !strcmp(argv[1], "strtod")) return cmd_strtod();
  if (argc == 2 && !strcmp(argv[1], "c2")) return cmd_c2();
  fprintf(stderr, "usage: imgref write|overwrite|chksum|addchk|esum|pix|subset|open|strtod|c2 ...\n");
  return 2;
}
