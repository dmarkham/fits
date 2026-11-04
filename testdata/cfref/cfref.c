/*
 * cfref -- cfitsio truth harness for the Go fits library test suite.
 *
 *   cfref dump <file.fits>      JSON description of every HDU, on stdout
 *   cfref dump -pool <dir> <file.fits>
 *                               the same, with mid-size value arrays stored
 *                               once in <dir> (see jarr: "Pool")
 *
 * The JSON is the independent truth that the Go tests assert against: it is
 * produced only by the system libcfitsio (4.7.0 when the committed goldens
 * were generated; see PROVENANCE.md). Nothing here uses the Go library.
 *
 * What is dumped, per HDU (cfitsio numbering is 1-based, "index" here is
 * 0-based to match fits.File.HDU(i)):
 *   - fits_get_hduaddrll: header start, data start, data end (cfitsio's
 *     "dataend" is the byte offset of the next HDU, i.e. data + block pad).
 *   - fits_verify_chksum: DATASUM / CHECKSUM status (1 ok, 0 absent, -1 bad).
 *   - every header card: raw 80-byte text, keyword name (fits_get_keyname),
 *     value string (fits_parse_value), comment, key class, key type, and the
 *     decoded value: strings via fits_read_key_longstr (CONTINUE resolved),
 *     integers via fits_read_key_lnglng, floats via fits_read_key_dbl,
 *     logicals via fits_read_key_log.
 *   - images (including tile-compressed images, which cfitsio presents as
 *     images): BITPIX, equivalent type, axes, and all pixels as
 *       "phys"  TDOUBLE, scaled (BSCALE/BZERO applied), null -> JSON null
 *       "i64"   TLONGLONG, scaled, for integer equivalent types (TLONG
 *               for tile-compressed images: cfitsio's decompressor rejects
 *               TLONGLONG; long is 64-bit on LP64)
 *       "u64"   TULONGLONG instead of i64 when the equivalent type is uint64
 *       "raw"   TLONGLONG with scaling disabled (integer BITPIX with a
 *               non-default BSCALE/BZERO only)
 *       "phys_f32"  "phys" rounded to float32 by a C cast; only for a
 *               BITPIX -32 image whose "phys" is a digest (see jdigest)
 *     A JSON null means cfitsio's null flag was set (BLANK match for
 *     integers, NaN for floats).
 *   - tables (binary and ASCII): per column TTYPE/TFORM/TUNIT/TDISP/TSCAL/
 *     TZERO/TNULL, cfitsio type code, equivalent type, repeat, width, TDIM
 *     (fits_read_tdimll), ASCII TBCOL; and per row the values:
 *       "str"   TSTRING (one entry per string in the cell), null -> null
 *       "log"   TLOGICAL as 1/0, undefined (0x00) -> null
 *       "bits"  TBIT, one 0/1 per bit
 *       "cplx"  TDBLCOMPLEX as [re, im] pairs
 *       "phys"  TDOUBLE scaled, null -> null
 *       "i64" / "u64" / "raw" as for images
 *     Variable-length (P/Q) cells also carry "desc": [length, heap offset]
 *     from fits_read_descriptll, and their values cover exactly length
 *     elements. ASCII tables additionally carry each row's raw text.
 *
 * Size: an array of more than DIGEST_MIN values (an image's pixels, or a
 * fixed-width column's values over all rows, then written as "cells"
 * instead of "rows") is stored as an exact sha256 digest plus explicit
 * samples; see jdigest. Per-row table cells are always explicit. With
 * -pool, an image array of 65..DIGEST_MIN values is a reference to a pool
 * file holding it explicitly (jarr).
 *
 * Non-finite doubles are emitted as the strings "NaN", "Inf", "-Inf" (JSON
 * has no literal for them). Strings are JSON-escaped; bytes >= 0x7F are
 * emitted as \u00XX meaning "the byte 0xXX" (not a Unicode code point).
 * Any cfitsio error is recorded in an "errors" array next to where it
 * happened and the dump continues; the harness never silently drops data.
 */
#define _POSIX_C_SOURCE 200809L
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#include <stdint.h>
#endif
#include <fitsio.h>

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

static void jdbl(double d) {
  if (isnan(d)) fputs("\"NaN\"", stdout);
  else if (isinf(d)) fputs(d > 0 ? "\"Inf\"" : "\"-Inf\"", stdout);
  else printf("%.17g", d);
}

static void errmsg(char *buf, int st) {
  char txt[FLEN_STATUS];
  fits_get_errstatus(st, txt);
  snprintf(buf, 160, "status %d: %s", st, txt);
  fits_clear_errmsg();
}

/* Per-HDU error collector. */
static char errs[256][256];
static int nerrs;
static void adderr(const char *where, int st) {
  if (nerrs >= 256) return;
  char m[160];
  errmsg(m, st);
  snprintf(errs[nerrs++], 256, "%s: %s", where, m);
}

static const char *hdutype(int t) {
  switch (t) {
  case IMAGE_HDU: return "IMAGE";
  case ASCII_TBL: return "ASCII_TBL";
  case BINARY_TBL: return "BINARY_TBL";
  }
  return "UNKNOWN";
}

static int is_int_type(int t) {
  t = abs(t);
  return t == TBYTE || t == TSBYTE || t == TSHORT || t == TUSHORT || t == TINT ||
         t == TUINT || t == TLONG || t == TULONG || t == TLONGLONG || t == TULONGLONG;
}

/* Position the header pointer just before card i (1-based) so a by-name read
 * finds card i itself, not an earlier duplicate. */
static void seek_card(fitsfile *f, int i) {
  char tmp[FLEN_CARD];
  int st = 0;
  fits_read_record(f, i - 1, tmp, &st);
}

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

/* The sha256 of the canonical encoding of n values, as hex. */
static void jhash(int kind, const void *v, const char *nul, LONGLONG n, char hex[65]) {
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
  sha_hex(&c, hex);
}

/* naxes: an image's axes (adds the plane corners to the sample), or NULL. */
static void jdigest(int kind, const void *v, const char *nul, LONGLONG n, int naxis, const LONGLONG *naxes) {
  char hex[65];
  jhash(kind, v, nul, n, hex);
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

/* Print n values: an explicit JSON array, or a digest past DIGEST_MIN, or
 * (with "dump -pool DIR") a pool reference.
 *
 * Pool: an array of more than POOL_MIN and at most DIGEST_MIN values (the
 * 960-pixel compression corpus images) keeps every value explicit, but
 * once per distinct array: its explicit JSON array (exactly what jarr
 * would print, plus a newline) goes to DIR/<kind>-<sha256>.json, written
 * once, and the golden holds {"n": N, "sha256": "<hex>", "pool": "<kind>"}
 * (kind f64, i64, u64, str, log, bit, cplx; sha256 as for a digest). The
 * funpack twins and the sources of the corpus hold the same pixels as the
 * compressed files, and an integer image's phys, i64 and raw are often
 * equal, so most arrays are stored once. internal/testutil loads the pool
 * file and checks its count and hash before a test sees the values. */
#define POOL_MIN 64
static const char *pool_dir;
static const char *const kind_name[] = {"f64", "i64", "u64", "str", "log", "bit", "cplx"};
static int cell_mode; /* set while dumping one table cell: always explicit */
static int digest_ok(LONGLONG n) { return n > DIGEST_MIN && !cell_mode; }
static void jarr_explicit(int kind, const void *v, const char *nul, LONGLONG n) {
  putchar('[');
  for (LONGLONG k = 0; k < n; k++) {
    if (k) putchar(',');
    jelem(kind, v, nul, k);
  }
  putchar(']');
}
static void jarr(int kind, const void *v, const char *nul, LONGLONG n, int naxis, const LONGLONG *naxes) {
  if (digest_ok(n)) { jdigest(kind, v, nul, n, naxis, naxes); return; }
  if (pool_dir && !cell_mode && n > POOL_MIN) {
    char hex[65], path[4096];
    jhash(kind, v, nul, n, hex);
    snprintf(path, sizeof path, "%s/%s-%s.json", pool_dir, kind_name[kind], hex);
    if (access(path, F_OK) != 0) {
      /* stdout to the pool file while the explicit array is printed */
      int fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0644), saved;
      if (fd < 0) { perror(path); exit(1); }
      fflush(stdout);
      saved = dup(1);
      dup2(fd, 1);
      close(fd);
      jarr_explicit(kind, v, nul, n);
      putchar('\n');
      fflush(stdout);
      dup2(saved, 1);
      close(saved);
    }
    printf("{\"n\": %lld, \"sha256\": \"%s\", \"pool\": \"%s\"}", n, hex, kind_name[kind]);
    return;
  }
  jarr_explicit(kind, v, nul, n);
}
/* card_body prints card i's cfitsio parse: the JSON object of a "cards"
 * entry without "n" and "raw" (key, class, value, comment, type and the
 * typed value). */
static void card_body(fitsfile *f, int i, char *card) {
  char name[FLEN_KEYWORD] = "", value[FLEN_VALUE] = "", comment[FLEN_COMMENT] = "";
  int len = 0, s = 0;
  fits_get_keyname(card, name, &len, &s);
  printf("{\"key\": "); jstr(name);
  printf(", \"class\": %d", fits_get_keyclass(card));
  s = 0;
  if (fits_parse_value(card, value, comment, &s)) {
    char m[160]; errmsg(m, s);
    printf(", \"parse_error\": "); jstr(m);
  }
  printf(", \"value\": "); jstr(value);
  printf(", \"comment\": "); jstr(comment);
  if (value[0]) {
    char dt = 0; s = 0;
    if (fits_get_keytype(value, &dt, &s) == 0) {
      printf(", \"type\": \"%c\"", dt);
      s = 0; seek_card(f, i);
      if (dt == 'C') {
        char *ls = NULL;
        if (fits_read_key_longstr(f, name, &ls, NULL, &s) == 0) {
          printf(", \"str\": "); jstr(ls);
        } else { char m[160]; errmsg(m, s); printf(", \"read_error\": "); jstr(m); }
        if (ls) { int s2 = 0; fits_free_memory(ls, &s2); }
      } else if (dt == 'I') {
        LONGLONG v;
        if (fits_read_key_lnglng(f, name, &v, NULL, &s) == 0) printf(", \"int\": %lld", v);
        else { char m[160]; errmsg(m, s); printf(", \"read_error\": "); jstr(m); }
        s = 0; seek_card(f, i);
        double d;
        if (fits_read_key_dbl(f, name, &d, NULL, &s) == 0) { printf(", \"float\": "); jdbl(d); }
      } else if (dt == 'F') {
        double d;
        if (fits_read_key_dbl(f, name, &d, NULL, &s) == 0) { printf(", \"float\": "); jdbl(d); }
        else { char m[160]; errmsg(m, s); printf(", \"read_error\": "); jstr(m); }
      } else if (dt == 'L') {
        int v;
        if (fits_read_key_log(f, name, &v, NULL, &s) == 0) printf(", \"logical\": %s", v ? "true" : "false");
        else { char m[160]; errmsg(m, s); printf(", \"read_error\": "); jstr(m); }
      }
    } else {
      char m[160]; errmsg(m, s);
      printf(", \"type_error\": "); jstr(m);
    }
  }
  printf("}");
}

/* card_poolable: whether card i's parse is a function of its 80 bytes
 * alone. Every field is, except the "str" of a string card whose value
 * ends in '&': fits_read_key_longstr appends the CONTINUE cards after it. */
static int card_poolable(char *card) {
  char value[FLEN_VALUE] = "", comment[FLEN_COMMENT] = "";
  int s = 0;
  if (fits_parse_value(card, value, comment, &s) || value[0] != '\'') return 1;
  size_t n = strlen(value);
  while (n > 1 && (value[n - 1] == '\'' || value[n - 1] == ' ')) n--;
  return !(n > 1 && value[n - 1] == '&');
}

/* Cards. With "dump -pool DIR" a card whose parse depends on its bytes
 * alone (card_poolable) is stored once: its parse (card_body plus a
 * newline) goes to DIR/card-<sha256 of the card text>.json, and the golden
 * holds just the card text as a JSON string, its "n" being its position
 * (internal/testutil reads the parse back). The same card text in another
 * file must give the same parse, or the dump fails. Any other card, and
 * every card without -pool, is printed in full:
 * {"n": i, "raw": ..., <card_body fields>}. */
static void dump_cards(fitsfile *f) {
  int nkeys = 0, more = 0, st = 0;
  fits_get_hdrspace(f, &nkeys, &more, &st);
  printf("   \"nkeys\": %d,\n   \"cards\": [", nkeys);
  for (int i = 1; i <= nkeys; i++) {
    char card[FLEN_CARD] = "";
    int s = 0;
    printf("%s\n    ", i > 1 ? "," : "");
    if (fits_read_record(f, i, card, &s)) {
      char m[160]; errmsg(m, s);
      printf("{\"n\": %d, \"error\": ", i); jstr(m); printf("}");
      continue;
    }
    if (pool_dir && card_poolable(card)) {
      char hex[65], path[4096], tmp[4200];
      sha256_ctx c;
      sha_init(&c);
      sha_update(&c, (const unsigned char *)card, strlen(card));
      sha_hex(&c, hex);
      snprintf(path, sizeof path, "%s/card-%s.json", pool_dir, hex);
      snprintf(tmp, sizeof tmp, "%s.tmp", path);
      int fd = open(tmp, O_WRONLY | O_CREAT | O_TRUNC, 0644), saved;
      if (fd < 0) { perror(tmp); exit(1); }
      fflush(stdout);
      saved = dup(1);
      dup2(fd, 1);
      close(fd);
      card_body(f, i, card);
      putchar('\n');
      fflush(stdout);
      dup2(saved, 1);
      close(saved);
      if (access(path, F_OK) == 0) {
        /* already pooled by another card: the parse must be the same */
        FILE *x = fopen(path, "rb"), *y = fopen(tmp, "rb");
        int same = x && y, cx, cy;
        while (same && (cx = fgetc(x)) == (cy = fgetc(y)) && cx != EOF) {}
        if (same) same = cx == cy;
        if (x) fclose(x);
        if (y) fclose(y);
        unlink(tmp);
        if (!same) { fprintf(stderr, "cfref: card pool conflict for %s (card %d)\n", path, i); exit(1); }
      } else if (rename(tmp, path) != 0) { perror(path); exit(1); }
      jstr(card);
      continue;
    }
    printf("{\"n\": %d, \"raw\": ", i); jstr(card);
    printf(", ");
    /* card_body's object without its opening brace */
    {
      fflush(stdout);
      char *buf = NULL; size_t blen = 0;
      FILE *ms = open_memstream(&buf, &blen);
      FILE *saved = stdout;
      stdout = ms;
      card_body(f, i, card);
      fflush(ms);
      stdout = saved;
      fclose(ms);
      fputs(buf + 1, stdout);
      free(buf);
    }
  }
  printf("\n   ],\n");
}

static void jarr_dbl(const double *d, const char *nul, LONGLONG n) { jarr(K_F64, d, nul, n, 0, NULL); }
static void jarr_ll(const LONGLONG *v, const char *nul, LONGLONG n) { jarr(K_I64, v, nul, n, 0, NULL); }
static void jarr_ull(const ULONGLONG *v, const char *nul, LONGLONG n) { jarr(K_U64, v, nul, n, 0, NULL); }

/* readimg: fits_read_imgnull, except for a tile-compressed image of more
 * than 3 dimensions: cfitsio 4.7.0's fits_read_img path for compressed
 * images (fits_read_compressed_pixels) refuses those ("only 1D, 2D, or 3D
 * images are currently supported", status 414), while its subset reader
 * (fits_read_compressed_img) takes up to MAX_COMPRESS_DIM = 6. Read the
 * whole image through the subset reader with null flags there. (Probed:
 * right for 4-D and 5-D; for 6-D it returns shifted pixels and zeros,
 * see the cube6d16 case in compressed_corpus_test.go.) */
static int readimg(fitsfile *f, int type, LONGLONG n, void *v, char *nul, int *anyn, int *st,
                   int compressed, int naxis, const LONGLONG *naxes) {
  if (!compressed || naxis <= 3)
    return fits_read_imgnull(f, type, 1, n, v, nul, anyn, st);
  long fp[MAX_COMPRESS_DIM], lp[MAX_COMPRESS_DIM], inc[MAX_COMPRESS_DIM];
  if (naxis > MAX_COMPRESS_DIM) return *st = BAD_NAXIS;
  for (int i = 0; i < naxis; i++) { fp[i] = 1; lp[i] = (long)naxes[i]; inc[i] = 1; }
  switch (type) {
  case TDOUBLE: return fits_read_subsetnull_dbl(f, 1, naxis, (long *)lp, fp, lp, inc, v, nul, anyn, st);
  case TLONG: return fits_read_subsetnull_lng(f, 1, naxis, (long *)lp, fp, lp, inc, v, nul, anyn, st);
  default: return *st = BAD_DATATYPE;
  }
}

static void dump_image(fitsfile *f) {
  int st = 0, bitpix = 0, naxis = 0, eq = 0;
  LONGLONG naxes[999];
  memset(naxes, 0, sizeof naxes);
  if (fits_get_img_paramll(f, 999, &bitpix, &naxis, naxes, &st)) { adderr("fits_get_img_paramll", st); return; }
  st = 0;
  if (fits_get_img_equivtype(f, &eq, &st)) { adderr("fits_get_img_equivtype", st); st = 0; }
  int compressed = fits_is_compressed_image(f, &st); st = 0;
  LONGLONG n = naxis ? 1 : 0;
  for (int i = 0; i < naxis; i++) n *= naxes[i];
  printf("   \"image\": {\"bitpix\": %d, \"equivtype\": %d, \"compressed\": %s, \"naxis\": %d, \"naxes\": [",
         bitpix, eq, compressed ? "true" : "false", naxis);
  for (int i = 0; i < naxis; i++) printf("%s%lld", i ? "," : "", naxes[i]);
  printf("], \"npix\": %lld", n);

  /* Tile-compressed images keep the original BSCALE/BZERO under the same
   * names in the table header; cfitsio applies them on read. */
  double bscale = 1.0, bzero = 0.0;
  int s = 0;
  if (fits_read_key_dbl(f, "BSCALE", &bscale, NULL, &s)) bscale = 1.0;
  s = 0;
  if (fits_read_key_dbl(f, "BZERO", &bzero, NULL, &s)) bzero = 0.0;
  s = 0;
  LONGLONG blank = 0;
  int hasblank = fits_read_key_lnglng(f, "BLANK", &blank, NULL, &s) == 0;
  printf(", \"bscale\": "); jdbl(bscale);
  printf(", \"bzero\": "); jdbl(bzero);
  if (hasblank) printf(", \"blank\": %lld", blank);

  if (n > 0) {
    double *d = calloc(1, n * sizeof(double));
    char *nul = calloc(1, n);
    int anyn = 0;
    st = 0;
    printf(",\n    \"phys\": ");
    if (readimg(f, TDOUBLE, n, d, nul, &anyn, &st, compressed, naxis, naxes)) { printf("null"); adderr("read_img TDOUBLE", st); }
    else {
      jarr(K_F64, d, nul, n, naxis, naxes);
      /* A digested BITPIX -32 image also gets "phys_f32": the same TDOUBLE
       * values rounded to float32 by a C cast, so a float32 reading can
       * be checked against the digest exactly (an explicit "phys" serves
       * both directly). */
      if (bitpix == FLOAT_IMG && digest_ok(n)) {
        for (LONGLONG k = 0; k < n; k++) if (!nul[k]) d[k] = (double)(float)d[k];
        printf(",\n    \"phys_f32\": ");
        jdigest(K_F64, d, nul, n, naxis, naxes);
      }
    }
    free(d);
    int intimg = bitpix > 0;
    if (intimg && eq != ULONGLONG_IMG && eq != FLOAT_IMG && eq != DOUBLE_IMG) {
      LONGLONG *v = calloc(1, n * sizeof(LONGLONG));
      st = 0; anyn = 0;
      printf(",\n    \"i64\": ");
      /* cfitsio's tile decompressor rejects TLONGLONG (status 410) but
       * accepts TLONG, which is 64-bit on the LP64 hosts we support. */
      int it = (compressed && sizeof(long) == 8) ? TLONG : TLONGLONG;
      if (readimg(f, it, n, v, nul, &anyn, &st, compressed, naxis, naxes)) { printf("null"); adderr(it == TLONG ? "read_img TLONG" : "read_img TLONGLONG", st); }
      else jarr(K_I64, v, nul, n, naxis, naxes);
      free(v);
    }
    if (intimg && eq == ULONGLONG_IMG) {
      ULONGLONG *v = calloc(1, n * sizeof(ULONGLONG));
      st = 0; anyn = 0;
      printf(",\n    \"u64\": ");
      if (fits_read_imgnull(f, TULONGLONG, 1, n, v, nul, &anyn, &st)) { printf("null"); adderr("read_img TULONGLONG", st); }
      else jarr(K_U64, v, nul, n, naxis, naxes);
      free(v);
    }
    if (intimg && (bscale != 1.0 || bzero != 0.0)) {
      LONGLONG *v = calloc(1, n * sizeof(LONGLONG));
      /* A tile-compressed image is read raw through a second handle on
       * the same file: libcfitsio 4.7.0 segfaults in
       * imcomp_decompress_tile when an HCOMPRESS_1 image with BZERO is
       * re-read through the same handle after fits_set_bscale (its tile
       * cache; e.g. comp_fz_u16sat_hcomp_s2). A fresh handle reads the
       * same stored values. */
      fitsfile *rf = f;
      st = 0; anyn = 0;
      if (compressed) {
        char fname[FLEN_FILENAME];
        int hdunum = 0;
        fits_file_name(f, fname, &st);
        fits_get_hdu_num(f, &hdunum);
        if (!st) fits_open_file(&rf, fname, READONLY, &st);
        if (!st) fits_movabs_hdu(rf, hdunum, NULL, &st);
      }
      if (!st) fits_set_bscale(rf, 1.0, 0.0, &st);
      printf(",\n    \"raw\": ");
      int it = (compressed && sizeof(long) == 8) ? TLONG : TLONGLONG;
      if (st || readimg(rf, it, n, v, nul, &anyn, &st, compressed, naxis, naxes)) { printf("null"); adderr("read_img raw", st); }
      else jarr(K_I64, v, nul, n, naxis, naxes);
      st = 0;
      if (rf != f) fits_close_file(rf, &st);
      else fits_set_bscale(f, bscale, bzero, &st);
      free(v);
    }
    free(nul);
  }
  printf("\n   },\n");
}

static int keystr(fitsfile *f, const char *root, int col, char *out) {
  char key[FLEN_KEYWORD];
  int s = 0;
  snprintf(key, sizeof key, "%s%d", root, col);
  out[0] = 0;
  return fits_read_key(f, TSTRING, key, out, NULL, &s) == 0;
}
static int keydbl(fitsfile *f, const char *root, int col, double *out) {
  char key[FLEN_KEYWORD];
  int s = 0;
  snprintf(key, sizeof key, "%s%d", root, col);
  return fits_read_key_dbl(f, key, out, NULL, &s) == 0;
}

/* Dump nelem values of column c, row r, according to the column's types. */
static void dump_cell(fitsfile *f, int hdut, int c, LONGLONG r, int tc, int eq,
                      LONGLONG rep, LONGLONG width, LONGLONG nelem, int scaled,
                      double tscal, double tzero) {
  int st = 0, anyn = 0;
  int at = abs(tc), aeq = abs(eq);
  char where[64];
  if (at == TSTRING) {
    /* Fixed binary rAw: rep/width strings of width chars; VLA PA: one string
     * of nelem chars; ASCII: one string of width chars. */
    LONGLONG ns, w;
    if (tc < 0) { ns = 1; w = nelem; }
    else if (hdut == ASCII_TBL) { ns = 1; w = width; }
    else { w = width > 0 ? width : 1; ns = rep / w; if (ns < 1 && rep > 0) ns = 1; }
    char **p = calloc(1, (ns ? ns : 1) * sizeof(char *));
    char *buf = calloc((ns ? ns : 1), w + 1);
    char *nul = calloc(ns ? ns : 1, 1);
    for (LONGLONG k = 0; k < ns; k++) p[k] = buf + k * (w + 1);
    printf("\"str\": ");
    if (ns == 0 || (tc < 0 && nelem == 0)) {
      /* zero-length VLA string: cfitsio returns an empty string */
      if (tc < 0) {
        if (fits_read_colnull(f, TSTRING, c, r, 1, 1, p, nul, &anyn, &st)) {
          printf("null"); snprintf(where, sizeof where, "col %d row %lld TSTRING", c, r); adderr(where, st);
        } else { printf("["); if (nul[0]) printf("null"); else jstr(p[0]); printf("]"); }
      } else printf("[]");
    } else if (fits_read_colnull(f, TSTRING, c, r, 1, ns, p, nul, &anyn, &st)) {
      printf("null"); snprintf(where, sizeof where, "col %d row %lld TSTRING", c, r); adderr(where, st);
    } else {
      putchar('[');
      for (LONGLONG k = 0; k < ns; k++) { if (k) putchar(','); if (nul[k]) printf("null"); else jstr(p[k]); }
      putchar(']');
    }
    free(p); free(buf); free(nul);
    return;
  }
  if (nelem == 0) {
    printf("\"n\": 0");
    return;
  }
  if (at == TLOGICAL) {
    char *v = calloc(1, nelem), *nul = calloc(1, nelem);
    printf("\"log\": ");
    if (fits_read_colnull(f, TLOGICAL, c, r, 1, nelem, v, nul, &anyn, &st)) {
      printf("null"); snprintf(where, sizeof where, "col %d row %lld TLOGICAL", c, r); adderr(where, st);
    } else {
      putchar('[');
      for (LONGLONG k = 0; k < nelem; k++) { if (k) putchar(','); if (nul[k]) printf("null"); else printf("%d", v[k] ? 1 : 0); }
      putchar(']');
    }
    free(v); free(nul);
    return;
  }
  if (at == TBIT) {
    char *v = calloc(1, nelem);
    printf("\"bits\": ");
    if (fits_read_col(f, TBIT, c, r, 1, nelem, NULL, v, &anyn, &st)) {
      printf("null"); snprintf(where, sizeof where, "col %d row %lld TBIT", c, r); adderr(where, st);
    } else {
      putchar('[');
      for (LONGLONG k = 0; k < nelem; k++) printf("%s%d", k ? "," : "", v[k] ? 1 : 0);
      putchar(']');
    }
    free(v);
    return;
  }
  if (at == TCOMPLEX || at == TDBLCOMPLEX) {
    double *v = calloc(1, 2 * nelem * sizeof(double));
    char *nul = calloc(2 * nelem, 1);
    printf("\"cplx\": ");
    if (fits_read_colnull(f, TDBLCOMPLEX, c, r, 1, nelem, v, nul, &anyn, &st)) {
      printf("null"); snprintf(where, sizeof where, "col %d row %lld TDBLCOMPLEX", c, r); adderr(where, st);
    } else {
      putchar('[');
      for (LONGLONG k = 0; k < nelem; k++) {
        if (k) putchar(',');
        putchar('['); jdbl(v[2 * k]); putchar(','); jdbl(v[2 * k + 1]); putchar(']');
      }
      putchar(']');
    }
    free(v); free(nul);
    return;
  }
  /* numeric */
  double *d = calloc(1, nelem * sizeof(double));
  char *nul = calloc(1, nelem);
  printf("\"phys\": ");
  st = 0;
  if (fits_read_colnull(f, TDOUBLE, c, r, 1, nelem, d, nul, &anyn, &st)) {
    printf("null"); snprintf(where, sizeof where, "col %d row %lld TDOUBLE", c, r); adderr(where, st);
  } else jarr_dbl(d, nul, nelem);
  free(d);
  if (is_int_type(at) && is_int_type(aeq)) {
    st = 0;
    if (aeq == TULONGLONG) {
      ULONGLONG *v = calloc(1, nelem * sizeof(ULONGLONG));
      printf(", \"u64\": ");
      if (fits_read_colnull(f, TULONGLONG, c, r, 1, nelem, v, nul, &anyn, &st)) {
        printf("null"); snprintf(where, sizeof where, "col %d row %lld TULONGLONG", c, r); adderr(where, st);
      } else jarr_ull(v, nul, nelem);
      free(v);
    } else {
      LONGLONG *v = calloc(1, nelem * sizeof(LONGLONG));
      printf(", \"i64\": ");
      if (fits_read_colnull(f, TLONGLONG, c, r, 1, nelem, v, nul, &anyn, &st)) {
        printf("null"); snprintf(where, sizeof where, "col %d row %lld TLONGLONG", c, r); adderr(where, st);
      } else jarr_ll(v, nul, nelem);
      free(v);
    }
  }
  if (is_int_type(at) && scaled) {
    LONGLONG *v = calloc(1, nelem * sizeof(LONGLONG));
    st = 0;
    fits_set_tscale(f, c, 1.0, 0.0, &st);
    printf(", \"raw\": ");
    if (fits_read_colnull(f, TLONGLONG, c, r, 1, nelem, v, nul, &anyn, &st)) {
      printf("null"); snprintf(where, sizeof where, "col %d row %lld raw TLONGLONG", c, r); adderr(where, st);
    } else jarr_ll(v, nul, nelem);
    st = 0; fits_set_tscale(f, c, tscal, tzero, &st);
    free(v);
  }
  free(nul);
}

/* A fixed-width column (not P/Q) whose rows hold more than DIGEST_MIN
 * elements in total is written as one "cells" object instead of "rows":
 *   "cells": {"per_row": k, "<kind>": <digest>, ...}
 * where each digest covers the column's values for all rows in row order
 * (row r, element j is index r*k + j), read row by row with exactly the
 * calls dump_cell uses, and k is the strings per cell for A columns and the
 * repeat count otherwise. If any read fails, nothing is printed and the
 * caller falls back to the per-row dump (which records the error). */
static int dump_flat_column(fitsfile *f, int hdut, int c, LONGLONG nrows, int tc, int eq,
                            LONGLONG rep, LONGLONG width, int scaled, double tscal, double tzero) {
  int at = abs(tc), aeq = abs(eq), st = 0, anyn = 0, ok = 1;
  if (tc < 0 || rep <= 0 || nrows <= 0) return 0;
  if (at == TSTRING) {
    LONGLONG ns, w;
    if (hdut == ASCII_TBL) { ns = 1; w = width; }
    else { w = width > 0 ? width : 1; ns = rep / w; if (ns < 1 && rep > 0) ns = 1; }
    LONGLONG n = ns * nrows;
    if (ns <= 0 || !digest_ok(n)) return 0;
    char **p = calloc(1, n * sizeof(char *));
    char *buf = calloc(n, w + 1);
    char *nul = calloc(n, 1);
    for (LONGLONG k = 0; k < n; k++) p[k] = buf + k * (w + 1);
    for (LONGLONG r = 1; r <= nrows && ok; r++) {
      st = 0;
      if (fits_read_colnull(f, TSTRING, c, r, 1, ns, p + (r - 1) * ns, nul + (r - 1) * ns, &anyn, &st)) ok = 0;
    }
    if (ok) {
      printf("\"cells\": {\"per_row\": %lld, \"str\": ", ns);
      jdigest(K_STR, p, nul, n, 0, NULL);
      printf("}");
    }
    fits_clear_errmsg();
    free(p); free(buf); free(nul);
    return ok;
  }
  LONGLONG n = rep * nrows;
  if (!digest_ok(n)) return 0;
  if (at == TLOGICAL || at == TBIT) {
    char *v = calloc(1, n), *nul = calloc(n, 1);
    for (LONGLONG r = 1; r <= nrows && ok; r++) {
      st = 0;
      if (at == TLOGICAL ? fits_read_colnull(f, TLOGICAL, c, r, 1, rep, v + (r - 1) * rep, nul + (r - 1) * rep, &anyn, &st)
                         : fits_read_col(f, TBIT, c, r, 1, rep, NULL, v + (r - 1) * rep, &anyn, &st)) ok = 0;
    }
    if (ok) {
      printf("\"cells\": {\"per_row\": %lld, \"%s\": ", rep, at == TLOGICAL ? "log" : "bits");
      jdigest(at == TLOGICAL ? K_LOG : K_BIT, v, nul, n, 0, NULL);
      printf("}");
    }
    fits_clear_errmsg();
    free(v); free(nul);
    return ok;
  }
  if (at == TCOMPLEX || at == TDBLCOMPLEX) {
    double *v = calloc(1, 2 * n * sizeof(double));
    char *nul = calloc(2 * n, 1);
    for (LONGLONG r = 1; r <= nrows && ok; r++) {
      st = 0;
      if (fits_read_colnull(f, TDBLCOMPLEX, c, r, 1, rep, v + 2 * (r - 1) * rep, nul + 2 * (r - 1) * rep, &anyn, &st)) ok = 0;
    }
    if (ok) {
      printf("\"cells\": {\"per_row\": %lld, \"cplx\": ", rep);
      jdigest(K_CPLX, v, NULL, n, 0, NULL);
      printf("}");
    }
    fits_clear_errmsg();
    free(v); free(nul);
    return ok;
  }
  /* numeric: read every kind first, print only if all reads succeed */
  int isint = is_int_type(at) && is_int_type(aeq), isu = isint && aeq == TULONGLONG, israw = is_int_type(at) && scaled;
  double *d = calloc(1, n * sizeof(double));
  char *nd = calloc(n, 1), *ni = calloc(n, 1), *nr = calloc(n, 1);
  LONGLONG *vi = isint && !isu ? calloc(1, n * sizeof(LONGLONG)) : NULL;
  ULONGLONG *vu = isu ? calloc(1, n * sizeof(ULONGLONG)) : NULL;
  LONGLONG *vr = israw ? calloc(1, n * sizeof(LONGLONG)) : NULL;
  for (LONGLONG r = 1; r <= nrows && ok; r++) {
    LONGLONG o = (r - 1) * rep;
    st = 0;
    if (fits_read_colnull(f, TDOUBLE, c, r, 1, rep, d + o, nd + o, &anyn, &st)) { ok = 0; break; }
    st = 0;
    if (vu && fits_read_colnull(f, TULONGLONG, c, r, 1, rep, vu + o, ni + o, &anyn, &st)) { ok = 0; break; }
    if (vi && fits_read_colnull(f, TLONGLONG, c, r, 1, rep, vi + o, ni + o, &anyn, &st)) { ok = 0; break; }
    if (vr) {
      st = 0;
      fits_set_tscale(f, c, 1.0, 0.0, &st);
      if (fits_read_colnull(f, TLONGLONG, c, r, 1, rep, vr + o, nr + o, &anyn, &st)) ok = 0;
      st = 0; fits_set_tscale(f, c, tscal, tzero, &st);
    }
  }
  if (ok) {
    printf("\"cells\": {\"per_row\": %lld, \"phys\": ", rep);
    jdigest(K_F64, d, nd, n, 0, NULL);
    if (vu) { printf(", \"u64\": "); jdigest(K_U64, vu, ni, n, 0, NULL); }
    if (vi) { printf(", \"i64\": "); jdigest(K_I64, vi, ni, n, 0, NULL); }
    if (vr) { printf(", \"raw\": "); jdigest(K_I64, vr, nr, n, 0, NULL); }
    printf("}");
  }
  fits_clear_errmsg();
  free(d); free(nd); free(ni); free(nr); free(vi); free(vu); free(vr);
  return ok;
}

static void dump_table(fitsfile *f, int hdut) {
  int st = 0, ncols = 0;
  LONGLONG nrows = 0, rowlen = 0;
  long rl = 0;
  fits_get_num_rowsll(f, &nrows, &st);
  fits_get_num_cols(f, &ncols, &st);
  if (fits_read_key_lng(f, "NAXIS1", &rl, NULL, &st) == 0) rowlen = rl;
  if (st) { adderr("table params", st); return; }
  printf("   \"table\": {\"nrows\": %lld, \"ncols\": %d, \"rowlen\": %lld, \"columns\": [", nrows, ncols, rowlen);
  for (int c = 1; c <= ncols; c++) {
    int tc = 0, eq = 0;
    LONGLONG rep = 0, width = 0, erep = 0, ewidth = 0;
    char ttype[FLEN_VALUE], tform[FLEN_VALUE], tunit[FLEN_VALUE], tdisp[FLEN_VALUE], tnull[FLEN_VALUE];
    st = 0;
    if (fits_get_coltypell(f, c, &tc, &rep, &width, &st)) { adderr("fits_get_coltypell", st); st = 0; }
    if (fits_get_eqcoltypell(f, c, &eq, &erep, &ewidth, &st)) { adderr("fits_get_eqcoltypell", st); st = 0; }
    int hasname = keystr(f, "TTYPE", c, ttype);
    keystr(f, "TFORM", c, tform);
    int hasunit = keystr(f, "TUNIT", c, tunit);
    int hasdisp = keystr(f, "TDISP", c, tdisp);
    int hasnull = keystr(f, "TNULL", c, tnull);
    double tscal = 1.0, tzero = 0.0;
    int hasscal = keydbl(f, "TSCAL", c, &tscal);
    int haszero = keydbl(f, "TZERO", c, &tzero);
    if (!hasscal) tscal = 1.0;
    if (!haszero) tzero = 0.0;
    int scaled = tscal != 1.0 || tzero != 0.0;
    printf("%s\n    {\"n\": %d, \"ttype\": ", c > 1 ? "," : "", c);
    if (hasname) jstr(ttype); else printf("null");
    printf(", \"tform\": "); jstr(tform);
    printf(", \"tunit\": "); if (hasunit) jstr(tunit); else printf("null");
    printf(", \"tdisp\": "); if (hasdisp) jstr(tdisp); else printf("null");
    printf(", \"tnull\": "); if (hasnull) jstr(tnull); else printf("null");
    printf(", \"tscal\": "); jdbl(tscal);
    printf(", \"tzero\": "); jdbl(tzero);
    printf(", \"typecode\": %d, \"eqtype\": %d, \"repeat\": %lld, \"width\": %lld", tc, eq, rep, width);
    if (hdut == ASCII_TBL) {
      long tbcol = 0; int s = 0;
      char key[FLEN_KEYWORD]; snprintf(key, sizeof key, "TBCOL%d", c);
      if (fits_read_key_lng(f, key, &tbcol, NULL, &s) == 0) printf(", \"tbcol\": %ld", tbcol);
    } else {
      int ndim = 0; LONGLONG dims[99]; int s = 0;
      if (fits_read_tdimll(f, c, 99, &ndim, dims, &s) == 0) {
        printf(", \"tdim\": [");
        for (int k = 0; k < ndim; k++) printf("%s%lld", k ? "," : "", dims[k]);
        printf("]");
      } else { adderr("fits_read_tdimll", s); }
    }
    printf(", ");
    if (dump_flat_column(f, hdut, c, nrows, tc, eq, rep, width, scaled, tscal, tzero)) {
      printf("}");
      continue;
    }
    printf("\"rows\": [");
    for (LONGLONG r = 1; r <= nrows; r++) {
      printf("%s\n     {", r > 1 ? "," : "");
      LONGLONG nelem = rep;
      if (tc < 0) {
        LONGLONG len = 0, off = 0; int s = 0;
        if (fits_read_descriptll(f, c, r, &len, &off, &s)) {
          char where[64]; snprintf(where, sizeof where, "col %d row %lld descriptor", c, r); adderr(where, s);
          printf("\"desc\": null}");
          continue;
        }
        printf("\"desc\": [%lld,%lld], ", len, off);
        nelem = len;
      } else if (abs(tc) == TSTRING) {
        nelem = rep;
      }
      cell_mode = 1;
      dump_cell(f, hdut, c, r, tc, eq, rep, width, nelem, scaled, tscal, tzero);
      cell_mode = 0;
      printf("}");
    }
    printf("\n    ]}");
  }
  printf("\n   ]");
  if (hdut == ASCII_TBL && rowlen > 0) {
    unsigned char *buf = calloc(1, rowlen);
    printf(",\n   \"raw_rows\": [");
    for (LONGLONG r = 1; r <= nrows; r++) {
      int s = 0;
      printf("%s\n    ", r > 1 ? "," : "");
      if (fits_read_tblbytes(f, r, 1, rowlen, buf, &s)) { printf("null"); adderr("fits_read_tblbytes", s); }
      else jstr_n((const char *)buf, rowlen);
    }
    printf("\n   ]");
    free(buf);
  }
  printf("},\n");
}

static int dump(const char *path) {
  fitsfile *f;
  int st = 0, nhdu = 0, type = 0;
  const char *base = strrchr(path, '/');
  base = base ? base + 1 : path;
  if (fits_open_diskfile(&f, path, READONLY, &st)) {
    fits_report_error(stderr, st);
    return 1;
  }
  fits_get_num_hdus(f, &nhdu, &st);
  printf("{\n \"file\": "); jstr(base);
  printf(",\n \"cfitsio\": \"%d.%d.%d\",\n \"nhdu\": %d,\n \"hdus\": [", CFITSIO_MAJOR, CFITSIO_MINOR, CFITSIO_MICRO, nhdu);
  for (int h = 1; h <= nhdu; h++) {
    nerrs = 0;
    st = 0;
    printf("%s\n  {\n   \"index\": %d", h > 1 ? "," : "", h - 1);
    if (fits_movabs_hdu(f, h, &type, &st)) {
      adderr("fits_movabs_hdu", st);
      printf(",\n   \"type\": null,\n");
    } else {
      LONGLONG hs = 0, ds = 0, de = 0;
      printf(",\n   \"type\": \"%s\"", hdutype(type));
      st = 0;
      if (fits_get_hduaddrll(f, &hs, &ds, &de, &st)) adderr("fits_get_hduaddrll", st);
      printf(",\n   \"header_start\": %lld,\n   \"data_start\": %lld,\n   \"data_end\": %lld", hs, ds, de);
      int dok = 0, hok = 0;
      st = 0;
      if (fits_verify_chksum(f, &dok, &hok, &st)) adderr("fits_verify_chksum", st);
      printf(",\n   \"checksum\": {\"datasum\": %d, \"hdusum\": %d},\n", dok, hok);
      dump_cards(f);
      if (type == IMAGE_HDU) dump_image(f);
      else dump_table(f, type);
    }
    printf("   \"errors\": [");
    for (int e = 0; e < nerrs; e++) { printf("%s", e ? ", " : ""); jstr(errs[e]); }
    printf("]\n  }");
  }
  printf("\n ]\n}\n");
  st = 0;
  fits_close_file(f, &st);
  return 0;
}

int main(int argc, char **argv) {
  /* glibc: run with the tcache off and malloc perturbation on, so that
   * every malloc'd block is filled with 0xA5 bytes (the complement of 90)
   * and every freed one with 0x5A. libcfitsio 4.7.0 returns uninitialized
   * memory in one place: a TDOUBLE fits_read_img of an integer image with
   * 8-byte GZIP tiles (imcomp_decompress_tile has no branch for that tile
   * type and copies its malloc'd tile buffer out unwritten; also a
   * quantized float image with such tiles). With the fill those pixels read
   * as the double whose 8 bytes are 0xA5 on every run, instead of whatever
   * the heap held (zeros, or heap addresses). The tunables are read at
   * startup, so the harness re-executes itself once with them set. */
  if (!getenv("CFREF_MALLOC_FILL")) {
    setenv("CFREF_MALLOC_FILL", "1", 1);
#ifdef __APPLE__
    /* macOS libmalloc: MallocPreScribble fills every allocated block with
     * 0xAA (there is no choice of byte), so the unwritten pixels read as
     * the double whose 8 bytes are 0xAA here, not 0xA5 as on glibc. The
     * goldens are recorded on glibc; on macOS those few readings differ by
     * design and tier 2 only reports them (see testutil.RequireGoldenPlatform). */
    setenv("MallocPreScribble", "1", 1);
    {
      char self[4096];
      uint32_t n = sizeof self;
      if (_NSGetExecutablePath(self, &n) != 0) {
        fprintf(stderr, "cfref: executable path too long\n");
        return 1;
      }
      execv(self, argv);
    }
#else
    setenv("GLIBC_TUNABLES", "glibc.malloc.tcache_count=0:glibc.malloc.perturb=90", 1);
    execv("/proc/self/exe", argv);
#endif
    perror("cfref: re-exec with the malloc fill");
    return 1;
  }
  if (argc == 3 && strcmp(argv[1], "dump") == 0) return dump(argv[2]);
  if (argc == 5 && strcmp(argv[1], "dump") == 0 && strcmp(argv[2], "-pool") == 0) {
    pool_dir = argv[3];
    return dump(argv[4]);
  }
  fprintf(stderr, "usage: cfref dump [-pool <dir>] <file.fits>\n");
  return 2;
}
