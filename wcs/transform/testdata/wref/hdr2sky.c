/*
 * hdr2sky: full-pipeline wcslib truth harness.
 *
 * Where wref.c exercises bare projections (prjset/prjx2s/prjs2x), this
 * program runs the SAME raw FITS header text the Go package parses through
 * wcslib's complete header pipeline:
 *
 *     wcspih  ->  wcsset  ->  wcsp2s (pixel -> sky)
 *                         ->  wcss2p (sky -> pixel)
 *
 * so the Go side (header.ParseCards -> wcs.ParseAlt -> transform.New ->
 * PixelToSky/SkyToPixel) can be checked end to end: keyword parsing,
 * defaults (LONPOLE/LATPOLE/PV), linear matrices, distortions, projection
 * and spherical rotation, AND which points are rejected.
 *
 * Usage:
 *     hdr2sky <case-file>        > golden/hdr/<name>.json
 *
 * Case file format (plain text, one directive per line, '#' comments):
 *
 *     GRID x0 x1 nx y0 y1 ny      pixel grid, 1-based FITS pixel coords,
 *                                 nx (ny) evenly spaced values x0..x1
 *     PIX x y                     one extra pixel point (1-based)
 *     SKYGRID a0 a1 na d0 d1 nd   sky grid, degrees
 *     SKY a d                     one extra sky point, degrees
 *     ALT c                       alternate WCS letter (default primary)
 *     TOL_ARCSEC t                optional p2s tolerance override, copied to
 *                                 the JSON for the Go test
 *     TOL_PX t                    optional s2p tolerance override (pixels)
 *     NOTE text                   free text copied into the JSON
 *     NOINJECT                    hand the cards to wcspih exactly as written
 *                                 (no SIP AP/BP injection, see CAVEAT); used
 *                                 by hdrcases_lenient to show wcslib's
 *                                 refusal
 *     HEADER                      everything after this line is header card
 *                                 text: one card per line, <= 80 chars,
 *                                 padded to 80 with blanks. Blank lines are
 *                                 skipped. END is appended.
 *
 * Sky -> pixel is run on every explicit SKY/SKYGRID point AND on the world
 * coordinate of every pixel point wcslib accepted ("from":"p2s"), so each
 * case carries a round trip through valid territory as well as whatever
 * off-projection sky points the case author adds.
 *
 * Output JSON (stdout), format "hdr2sky-2": the wcslib version, sha256
 * digests of the case file text, of the 80-column cards handed to wcspih
 * (without the harness-injected ones, see below) and of the INPUT
 * coordinates (every pixel point x, y then every explicit sky point lng,
 * lat, as big-endian IEEE-754 binary64), the note / ALT / tolerance
 * directives, a summary of what wcsset derived (ctype, lonpole, latpole,
 * crval, cel.ref), a readable sample of each direction (see pickSample:
 * the first and last point, three strided interior points, the first two
 * rejected points and, among the accepted ones, the points holding the
 * smallest and largest value of each output coordinate; %.17g, non-finite
 * values as null), and the per-point results
 * as base64 packed records: status byte, stat byte and, only when status
 * == 0, the two output coordinates as big-endian binary64 (p2s: lng, lat;
 * s2p: x, y). The input coordinates are not stored: the Go test rebuilds
 * them from the case directives and the inputs digest proves the rebuild
 * is bit-identical. The s2p inputs are the explicit sky points followed by
 * the world coordinates of the accepted p2s points.
 *
 * --------------------------------------------------------------------------
 * CAVEAT: SIP without AP/BP.
 *
 * wcslib 8.10's SIP reader (dis.c, sipset) rejects a header that has
 * A_ORDER/B_ORDER but no AP_ORDER/BP_ORDER ("Degree must not exceed 9,
 * got -1"), i.e. it cannot express "forward polynomial only, invert
 * numerically". wcslib's disx2p always inverts the forward polynomial by
 * iteration and uses the AP/BP polynomial only as the starting guess. So
 * for such headers this harness INJECTS
 *
 *     AP_ORDER = 1, AP_0_0 = 0.0, BP_ORDER = 1, BP_0_0 = 0.0
 *
 * (an identity starting guess) into the copy handed to wcspih. The
 * injected cards are listed under "harness_injected" in the JSON and are
 * NOT part of "cards": the Go side parses the original header without
 * AP/BP. The resulting truth is wcslib's iterated inverse of A/B started
 * from the identity -- the exact inverse to wcslib's convergence
 * tolerance, which is what a SIP header without AP/BP means. If wcslib's
 * iteration fails to converge from the identity start it reports a
 * non-zero status for that point and the golden records the rejection.
 * --------------------------------------------------------------------------
 *
 * Requires: libwcs (Arch: extra/wcslib). Built by the Makefile in this
 * directory; the binary is a build artifact, the .c and goldens are
 * committed. See PROVENANCE.md.
 */

#define _DEFAULT_SOURCE /* strsep, strdup under -std=c99 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wcslib/wcslib.h>

#define MAXCARDS 4000
#define MAXPTS 200000

static double *pixx, *pixy, *skya, *skyd;
static int npix, nsky;

static void die(const char *msg, const char *arg) {
  fprintf(stderr, "hdr2sky: %s%s%s\n", msg, arg ? ": " : "", arg ? arg : "");
  exit(1);
}

static void addpix(double x, double y) {
  if (npix >= MAXPTS) die("too many pixel points", NULL);
  pixx[npix] = x; pixy[npix] = y; npix++;
}

static void addsky(double a, double d) {
  if (nsky >= MAXPTS) die("too many sky points", NULL);
  skya[nsky] = a; skyd[nsky] = d; nsky++;
}

/* grid helper: n evenly spaced values from v0 to v1 inclusive */
static double lin(double v0, double v1, int n, int i) {
  return n == 1 ? v0 : v0 + (v1 - v0) * (double)i / (double)(n - 1);
}

static void jstr(const char *s) {
  putchar('"');
  for (; *s; s++) {
    unsigned char c = (unsigned char)*s;
    if (c == '"' || c == '\\') { putchar('\\'); putchar(c); }
    else if (c == '\n') fputs("\\n", stdout);
    else if (c == '\t') fputs("\\t", stdout);
    else if (c < 0x20) printf("\\u%04x", c);
    else putchar(c);
  }
  putchar('"');
}

static void jnum(double v) {
  if (isfinite(v)) printf("%.17g", v);
  else fputs("null", stdout);
}

/* card keyword (first 8 chars, trailing blanks trimmed) equals kw */
static int haskw(const char *cards, int n, const char *kw) {
  size_t L = strlen(kw);
  for (int i = 0; i < n; i++) {
    const char *c = cards + 80 * i;
    if (strncmp(c, kw, L) == 0 && (L == 8 || c[L] == ' ' || c[L] == '='))
      return 1;
  }
  return 0;
}

/* big-endian IEEE-754 binary64 */
static void putbe(unsigned char *q, double v) {
  unsigned long long u; memcpy(&u, &v, 8);
  for (int k = 7; k >= 0; k--) { q[k] = (unsigned char)(u & 0xff); u >>= 8; }
}

static void b64out(const unsigned char *b, size_t n) {
  static const char tb[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  for (size_t i = 0; i < n; i += 3) {
    unsigned v = (unsigned)b[i] << 16;
    if (i + 1 < n) v |= (unsigned)b[i + 1] << 8;
    if (i + 2 < n) v |= b[i + 2];
    putchar(tb[(v >> 18) & 63]); putchar(tb[(v >> 12) & 63]);
    putchar(i + 1 < n ? tb[(v >> 6) & 63] : '=');
    putchar(i + 2 < n ? tb[v & 63] : '=');
  }
}

/* SHA-256 (FIPS 180-4), for the digests in the golden. */
static const unsigned int shaK[64] = {
  0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
  0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
  0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
  0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
  0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
  0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
  0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
  0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
#define ROR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))
static void sha256blk(unsigned int *h, const unsigned char *p) {
  unsigned int w[64], a, b, c, d, e, f, g, hh;
  for (int i = 0; i < 16; i++)
    w[i] = (unsigned)p[4*i] << 24 | (unsigned)p[4*i+1] << 16 | (unsigned)p[4*i+2] << 8 | p[4*i+3];
  for (int i = 16; i < 64; i++) {
    unsigned s0 = ROR(w[i-15], 7) ^ ROR(w[i-15], 18) ^ (w[i-15] >> 3);
    unsigned s1 = ROR(w[i-2], 17) ^ ROR(w[i-2], 19) ^ (w[i-2] >> 10);
    w[i] = w[i-16] + s0 + w[i-7] + s1;
  }
  a = h[0]; b = h[1]; c = h[2]; d = h[3]; e = h[4]; f = h[5]; g = h[6]; hh = h[7];
  for (int i = 0; i < 64; i++) {
    unsigned t1 = hh + (ROR(e, 6) ^ ROR(e, 11) ^ ROR(e, 25)) + ((e & f) ^ (~e & g)) + shaK[i] + w[i];
    unsigned t2 = (ROR(a, 2) ^ ROR(a, 13) ^ ROR(a, 22)) + ((a & b) ^ (a & c) ^ (b & c));
    hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
  }
  h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
}
static void sha256hex(const unsigned char *p, size_t n, char out[65]) {
  unsigned int h[8] = {0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
  unsigned char blk[128];
  size_t i = 0;
  for (; i + 64 <= n; i += 64) sha256blk(h, p + i);
  size_t r = n - i;
  memset(blk, 0, sizeof blk);
  memcpy(blk, p + i, r);
  blk[r] = 0x80;
  size_t tot = r + 1 + 8 <= 64 ? 64 : 128;
  unsigned long long bits = (unsigned long long)n * 8;
  for (int k = 0; k < 8; k++) blk[tot - 1 - k] = (unsigned char)(bits >> (8 * k));
  for (size_t k = 0; k < tot; k += 64) sha256blk(h, blk + k);
  for (int k = 0; k < 8; k++) sprintf(out + 8 * k, "%08x", h[k]);
}

/* pickSample marks the points of one direction that go into the readable
 * sample: the first and last, three strided interior points (n/4, n/2,
 * 3n/4), the first two points wcslib rejected (status != 0) and, among
 * the accepted points, the ones holding the minimum and maximum of each
 * output coordinate (first such index; NaN outputs never win). */
static void pickSample(int n, const unsigned char *st, const double *v, unsigned char *mark) {
  memset(mark, 0, (size_t)n + 1);
  if (n == 0) return;
  mark[0] = mark[n - 1] = 1;
  for (int k = 1; k <= 3; k++) mark[(int)((long long)n * k / 4)] = 1;
  for (int i = 0, r = 0; i < n && r < 2; i++) if (st[2 * i]) { mark[i] = 1; r++; }
  for (int c = 0; c < 2; c++) {
    int lo = -1, hi = -1;
    for (int i = 0; i < n; i++) {
      if (st[2 * i] || !isfinite(v[2 * i + c])) continue;
      if (lo < 0 || v[2 * i + c] < v[2 * lo + c]) lo = i;
      if (hi < 0 || v[2 * i + c] > v[2 * hi + c]) hi = i;
    }
    if (lo >= 0) mark[lo] = 1;
    if (hi >= 0) mark[hi] = 1;
  }
}

int main(int argc, char **argv) {
  if (argc != 2) { fprintf(stderr, "usage: hdr2sky case-file > golden.json\n"); return 2; }
  const char *path = argv[1];

  /* read whole case file (kept verbatim for the JSON) */
  FILE *f = fopen(path, "rb");
  if (!f) die("open", path);
  fseek(f, 0, SEEK_END); long fsz = ftell(f); fseek(f, 0, SEEK_SET);
  char *text = malloc(fsz + 1);
  if (fread(text, 1, fsz, f) != (size_t)fsz) die("read", path);
  text[fsz] = 0; fclose(f);

  pixx = malloc(sizeof(double) * MAXPTS); pixy = malloc(sizeof(double) * MAXPTS);
  skya = malloc(sizeof(double) * 2 * MAXPTS); skyd = malloc(sizeof(double) * 2 * MAXPTS);

  char *cards = calloc(1, 80 * (MAXCARDS + 8));
  int ncards = 0, inhdr = 0;
  char alt = ' ';
  double tolas = NAN, tolpx = NAN;
  char note[1024] = "";
  int noinject = 0;

  /* line-by-line parse of a copy (strtok-like on '\n') */
  char *copy = strdup(text), *save = copy, *line;
  int lineno = 0;
  while ((line = strsep(&save, "\n")) != NULL) {
    lineno++;
    size_t L = strcspn(line, "\r");
    line[L] = 0;
    if (inhdr) {
      if (L == 0) continue;
      if (L > 80) { char m[64]; snprintf(m, sizeof m, "line %d", lineno); die("header card longer than 80 chars", m); }
      if (ncards >= MAXCARDS) die("too many cards", NULL);
      char *c = cards + 80 * ncards;
      memset(c, ' ', 80); memcpy(c, line, L);
      if (strncmp(c, "END     ", 8) == 0) break;
      ncards++;
      continue;
    }
    char *p = line; while (*p == ' ' || *p == '\t') p++;
    if (*p == 0 || *p == '#') continue;
    double a0, a1, b0, b1; int na, nb; char cbuf[8];
    if (strncmp(p, "HEADER", 6) == 0) { inhdr = 1; }
    else if (sscanf(p, "GRID %lf %lf %d %lf %lf %d", &a0, &a1, &na, &b0, &b1, &nb) == 6) {
      if (na < 1 || nb < 1) die("bad GRID", p);
      for (int j = 0; j < nb; j++) for (int i = 0; i < na; i++) addpix(lin(a0, a1, na, i), lin(b0, b1, nb, j));
    } else if (sscanf(p, "PIX %lf %lf", &a0, &b0) == 2) { addpix(a0, b0); }
    else if (sscanf(p, "SKYGRID %lf %lf %d %lf %lf %d", &a0, &a1, &na, &b0, &b1, &nb) == 6) {
      if (na < 1 || nb < 1) die("bad SKYGRID", p);
      for (int j = 0; j < nb; j++) for (int i = 0; i < na; i++) addsky(lin(a0, a1, na, i), lin(b0, b1, nb, j));
    } else if (sscanf(p, "SKY %lf %lf", &a0, &b0) == 2) { addsky(a0, b0); }
    else if (sscanf(p, "ALT %1s", cbuf) == 1) { alt = cbuf[0]; }
    else if (sscanf(p, "TOL_ARCSEC %lf", &a0) == 1) { tolas = a0; }
    else if (sscanf(p, "TOL_PX %lf", &a0) == 1) { tolpx = a0; }
    else if (strcmp(p, "NOINJECT") == 0) { noinject = 1; }
    else if (strncmp(p, "NOTE ", 5) == 0) { snprintf(note, sizeof note, "%s", p + 5); }
    else { char m[64]; snprintf(m, sizeof m, "line %d", lineno); die("unknown directive", m); }
  }
  if (!inhdr || ncards == 0) die("case file has no HEADER section", path);
  if (npix == 0 && nsky == 0) die("case file has no points", path);

  /* Copy for wcspih, with the SIP-without-AP/BP workaround (see top). */
  char *wcards = calloc(1, 80 * (MAXCARDS + 8));
  memcpy(wcards, cards, 80 * ncards);
  int nw = ncards;
  const char *injected[8]; int ninj = 0;
  if (alt == ' ' && !noinject) {
    static const char *apfix[] = {"AP_ORDER= 1", "AP_0_0  = 0.0"};
    static const char *bpfix[] = {"BP_ORDER= 1", "BP_0_0  = 0.0"};
    if (haskw(cards, ncards, "A_ORDER") && !haskw(cards, ncards, "AP_ORDER"))
      for (int k = 0; k < 2; k++) injected[ninj++] = apfix[k];
    if (haskw(cards, ncards, "B_ORDER") && !haskw(cards, ncards, "BP_ORDER"))
      for (int k = 0; k < 2; k++) injected[ninj++] = bpfix[k];
    for (int k = 0; k < ninj; k++) {
      char *c = wcards + 80 * nw++;
      memset(c, ' ', 80); memcpy(c, injected[k], strlen(injected[k]));
    }
  }
  { char *c = wcards + 80 * nw++; memset(c, ' ', 80); memcpy(c, "END", 3); }

  int nreject = 0, nwcs = 0; struct wcsprm *wcs = NULL;
  int st = wcspih(wcards, nw, WCSHDR_all, 0, &nreject, &nwcs, &wcs);
  if (st) { char m[32]; snprintf(m, sizeof m, "status %d", st); die("wcspih failed", m); }
  struct wcsprm *w = NULL;
  for (int i = 0; i < nwcs; i++) if (wcs[i].alt[0] == alt) w = &wcs[i];
  if (!w) die("requested alternate WCS not found in header", NULL);
  wcserr_enable(1);
  if ((st = wcsset(w))) { wcsperr(w, "hdr2sky: "); die("wcsset failed", NULL); }
  if (w->lng < 0 || w->lat < 0) die("header has no celestial axis pair", NULL);

  int naxis = w->naxis;
  double *pix = calloc(naxis, sizeof(double)), *img = calloc(naxis, sizeof(double));
  double *wld = calloc(naxis, sizeof(double)), phi, theta; int stat;

  /* ---- compute every point first, then write the compact JSON ---- */
  int nexp = nsky; /* explicit sky points; accepted p2s results follow */
  unsigned char *p2st = calloc(npix + 1, 2), *s2st;
  double *p2v = calloc(2 * (size_t)npix + 2, sizeof(double)), *s2v;
  for (int i = 0; i < npix; i++) {
    /* The case's (x, y) drive the pixel axes that carry the celestial
     * world axes (pixel axis i <-> world axis i); any other axis sits at
     * pixel 1. */
    for (int k = 0; k < naxis; k++) pix[k] = 1.0;
    pix[w->lng] = pixx[i]; pix[w->lat] = pixy[i];
    stat = 0;
    int s = wcsp2s(w, 1, naxis, pix, img, &phi, &theta, wld, &stat);
    if (s < 0 || s > 255 || stat < 0 || stat > 255) die("status out of byte range", NULL);
    p2st[2 * i] = (unsigned char)s; p2st[2 * i + 1] = (unsigned char)stat;
    if (s == 0) {
      p2v[2 * i] = wld[w->lng]; p2v[2 * i + 1] = wld[w->lat];
      addsky(wld[w->lng], wld[w->lat]);
    }
  }
  s2st = calloc(nsky + 1, 2);
  s2v = calloc(2 * (size_t)nsky + 2, sizeof(double));
  for (int i = 0; i < nsky; i++) {
    for (int k = 0; k < naxis; k++) wld[k] = w->crval[k];
    wld[w->lng] = skya[i]; wld[w->lat] = skyd[i];
    stat = 0;
    int s = wcss2p(w, 1, naxis, wld, &phi, &theta, img, pix, &stat);
    if (s < 0 || s > 255 || stat < 0 || stat > 255) die("status out of byte range", NULL);
    s2st[2 * i] = (unsigned char)s; s2st[2 * i + 1] = (unsigned char)stat;
    if (s == 0) { s2v[2 * i] = pix[w->lng]; s2v[2 * i + 1] = pix[w->lat]; }
  }

  /* Digests. inputs_sha256 covers the INPUT coordinates the Go side
   * regenerates from the case directives: every pixel point (x, y), then
   * every explicit sky point (lng, lat), as big-endian IEEE-754 binary64. */
  char hcase[65], hcards[65], hin[65];
  sha256hex((const unsigned char *)text, (size_t)fsz, hcase);
  sha256hex((const unsigned char *)cards, 80 * (size_t)ncards, hcards);
  {
    unsigned char *b = malloc(16 * ((size_t)npix + (size_t)nexp) + 1), *q = b;
    for (int i = 0; i < npix; i++) { putbe(q, pixx[i]); putbe(q + 8, pixy[i]); q += 16; }
    for (int i = 0; i < nexp; i++) { putbe(q, skya[i]); putbe(q + 8, skyd[i]); q += 16; }
    sha256hex(b, (size_t)(q - b), hin);
    free(b);
  }

  int vers[3];
  printf("{\"format\": \"hdr2sky-2\",\n \"wcslib_version\": ");
  jstr(wcslib_version(vers));
  printf(",\n \"case_sha256\": \"%s\", \"cards_sha256\": \"%s\", \"ncards\": %d,\n", hcase, hcards, ncards);
  printf(" \"note\": "); jstr(note);
  { char a[2] = {alt == ' ' ? 0 : alt, 0}; printf(",\n \"alt\": "); jstr(a); }
  printf(", \"tol_arcsec\": "); jnum(tolas);
  printf(", \"tol_px\": "); jnum(tolpx);
  printf(",\n \"harness_injected\": [");
  for (int k = 0; k < ninj; k++) { printf(k ? ", " : ""); jstr(injected[k]); }
  printf("],\n \"wcs\": {\"naxis\": %d, \"lng\": %d, \"lat\": %d, \"ctype\": [", naxis, w->lng, w->lat);
  for (int i = 0; i < naxis; i++) { printf(i ? ", " : ""); jstr(w->ctype[i]); }
  printf("], \"crval\": ["); jnum(w->crval[w->lng]); printf(", "); jnum(w->crval[w->lat]);
  printf("], \"lonpole\": "); jnum(w->lonpole);
  printf(", \"latpole\": "); jnum(w->latpole);
  printf(", \"cel_ref\": ["); for (int i = 0; i < 4; i++) { printf(i ? ", " : ""); jnum(w->cel.ref[i]); }
  printf("]},\n");
  printf(" \"npix\": %d, \"nsky_case\": %d, \"inputs_sha256\": \"%s\",\n", npix, nexp, hin);

  /* Readable sample (pickSample, %.17g); the Go test checks it against
   * the packed records. */
  unsigned char *mark = malloc((size_t)(npix > nsky ? npix : nsky) + 1);
  pickSample(npix, p2st, p2v, mark);
  printf(" \"sample\": {\"p2s\": [");
  for (int n = 0, i = 0; i < npix; i++) {
    if (!mark[i]) continue;
    printf(n++ ? ", " : ""); printf("{\"i\": %d, \"x\": ", i); jnum(pixx[i]); printf(", \"y\": "); jnum(pixy[i]);
    printf(", \"status\": %d, \"lng\": ", p2st[2 * i]);
    if (p2st[2 * i] == 0) { jnum(p2v[2 * i]); printf(", \"lat\": "); jnum(p2v[2 * i + 1]); } else printf("null, \"lat\": null");
    printf("}");
  }
  printf("], \"s2p\": [");
  pickSample(nsky, s2st, s2v, mark);
  for (int n = 0, i = 0; i < nsky; i++) {
    if (!mark[i]) continue;
    printf(n++ ? ", " : ""); printf("{\"i\": %d, \"lng\": ", i); jnum(skya[i]); printf(", \"lat\": "); jnum(skyd[i]);
    printf(", \"status\": %d, \"x\": ", s2st[2 * i]);
    if (s2st[2 * i] == 0) { jnum(s2v[2 * i]); printf(", \"y\": "); jnum(s2v[2 * i + 1]); } else printf("null, \"y\": null");
    printf("}");
  }
  printf("]},\n");
  free(mark);

  /* Packed records, base64: per point status (1 byte), stat (1 byte),
   * then, only when status == 0, the two output coordinates as big-endian
   * binary64 (p2s: lng, lat in degrees; s2p: x, y 1-based pixels). */
  for (int dir = 0; dir < 2; dir++) {
    int n = dir ? nsky : npix;
    unsigned char *st = dir ? s2st : p2st;
    double *v = dir ? s2v : p2v;
    unsigned char *b = malloc(18 * (size_t)n + 1), *q = b;
    for (int i = 0; i < n; i++) {
      *q++ = st[2 * i]; *q++ = st[2 * i + 1];
      if (st[2 * i] == 0) { putbe(q, v[2 * i]); putbe(q + 8, v[2 * i + 1]); q += 16; }
    }
    printf(" \"%s\": \"", dir ? "s2p" : "p2s");
    b64out(b, (size_t)(q - b));
    printf("\"%s\n", dir ? "" : ",");
    free(b);
  }
  printf("}\n");

  wcsvfree(&nwcs, &wcs);
  return 0;
}
