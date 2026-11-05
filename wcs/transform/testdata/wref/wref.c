/*
 * wref.c -- wcslib bare-projection truth harness.
 *
 * Build:   make wref            (links the system libwcs, see PROVENANCE.md)
 * Use:     ./wref gen golden    (golden/<case>.json + .bin per case below)
 *          ./wref trig libm     (libm/wcstrig_guards.json, see write_trig)
 *          ./wref list          (case names, one per line)
 *
 * For every case this sets up a bare wcslib projection the way a FITS
 * header would (prjini, the code, then ONLY the PVi_m values the case
 * names; every other parameter keeps wcslib's default) and runs prjset,
 * prjs2x and prjx2s on fixed grids. The Go test (wref_test.go) builds the
 * same projection through transform.Select with the same PV keys (PV2_m on
 * the latitude axis) and requires, point by point:
 *
 *   - the same accept/reject decision (wcslib stat != 0 <=> Go ok=false);
 *   - the same coordinates, within a stated budget, where both accept.
 *
 * Angles are degrees on the wcslib side (its native unit, r0 = 180/pi);
 * the Go side works in radians and the test converts.
 *
 * Grids (identical for every case, so nothing is tailored to avoid a
 * boundary):
 *
 *   forward  (phi, theta): phi = -180..180 step 15 plus phi = 270 and
 *            phi = -270 (wcslib does not reduce phi; neither may we),
 *            theta = -90..90 step 15. Includes both poles, the equator,
 *            phi = +-180 and phi = +-90.
 *   inverse  (x, y): a coarse box [-400, 400]^2 step 50 deg, a fine box
 *            [-180, 180]^2 step 15 deg (lands exactly on many projection
 *            boundaries, e.g. x = +-180, y = +-90, x = +-45), and the (x, y)
 *            that wcslib's own forward produced for every accepted forward
 *            point (a round trip through valid territory).
 *
 * Output, per case (golden/<case>.json + golden/<case>.bin):
 *
 *   <case>.json, one small metadata object:
 *     { "name", "code", "wcslib_version", "zlib_version",
 *       "pv": [[m, value], ...]   the PVi_m given (all others default),
 *       "prjset": {phi0, theta0, x0, y0, r0, pv_after: prj.pv[0..3] after
 *                  prjset filled in defaults (null = UNDEFINED)},
 *       "forward": nf, "inverse": ni_own, "roundtrip": nrt,
 *       "payload": "<case>.bin", "payload_bytes": <uncompressed size> }
 *     Doubles here are the shortest of %.15g/%.16g/%.17g that reads back
 *     exactly.
 *
 *   <case>.bin, the points: zlib (RFC 1950, level 9) of a big-endian
 *   stream, every double as its exact IEEE-754 bits:
 *     u32 nf, u32 ni_own, u32 nrt
 *     nf x (f64 phi, f64 theta)          forward inputs
 *     nf x u8 stat                       wcslib prjs2x stat per row
 *     (f64 x, f64 y) for every row with stat == 0 (raw: a NaN, an
 *                                        infinity or UNDEFINED is kept)
 *     ni_own x (f64 x, f64 y)            inverse inputs (grids, probes)
 *     (ni_own + nrt) x u8 stat           wcslib prjx2s stat per row
 *     (f64 phi, f64 theta) for every inverse row with stat == 0
 *   The last nrt inverse rows are the round trips: their inputs are the
 *   (x, y) of the accepted forward rows with finite outputs, in forward
 *   order, and are not stored twice (the reader rebuilds them). Nothing is
 *   rounded: the reader recovers every input and wcslib output bit for
 *   bit. The previous one-row-per-line JSON (5.4 MB for 62 cases) held
 *   the same numbers in decimal; this is ~8 times smaller.
 *
 * History: the previous version of this file (commit 9d21799) put the
 * AZP/SZP parameters one index low (prj.pv[0] is not PV2_1), so
 * "azp_tilted" was really mu = 30, gamma = 0 and every SZP case had
 * theta_c = 90; its TSC case read 14 inverse points from a 13-point grid
 * (out of bounds); its forward grid never reached phi > 120; and it covered
 * 8 of the 28 projection codes. All of that is replaced here.
 */

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wcslib/prj.h>
#include <wcslib/wcsmath.h>
#include <wcslib/wcstrig.h>
#include <wcslib/wcs.h>
#include <zlib.h>

typedef struct {
  int m;
  double v;
} pvset_t;

typedef struct {
  const char *name;
  const char *code;
  int npv;
  pvset_t pv[6];
} case_t;

/* Every projection code wcslib 8.10 implements (Paper II + HEALPix, 28
 * codes) appears at least once with defaults; non-default parameter sets
 * follow. */
static const case_t cases[] = {
    /* ---- zenithal ---- */
    {"azp_default", "AZP", 0, {{0}}},
    {"azp_mu2", "AZP", 1, {{1, 2.0}}},
    {"azp_mu2_gamma30", "AZP", 2, {{1, 2.0}, {2, 30.0}}},
    {"azp_mu1p5_gammam25", "AZP", 2, {{1, 1.5}, {2, -25.0}}},
    {"azp_mum3_gamma20", "AZP", 2, {{1, -3.0}, {2, 20.0}}},
    /* mu = 1 with a tilt: azpx2s's theta = max(a, b) lands a few 1e-13 ..
     * 1e-9 deg below -90 next to the horizon, where the final prjbchk(1e-13)
     * decides (round 2 of the proof audit: a 1e-6 tolerance survived every
     * other case; a wcslib search found such inputs only for mu = 1). */
    {"azp_mu1_gamma45", "AZP", 2, {{1, 1.0}, {2, 45.0}}},
    {"szp_default", "SZP", 0, {{0}}},
    {"szp_mu2_phic180_thetac60", "SZP", 3, {{1, 2.0}, {2, 180.0}, {3, 60.0}}},
    {"szp_mu1p5_phic90_thetac45", "SZP", 3, {{1, 1.5}, {2, 90.0}, {3, 45.0}}},
    {"szp_mu0_phic30_thetac70", "SZP", 3, {{1, 0.0}, {2, 30.0}, {3, 70.0}}},
    {"tan_default", "TAN", 0, {{0}}},
    {"stg_default", "STG", 0, {{0}}},
    {"sin_default", "SIN", 0, {{0}}},
    {"sin_slant", "SIN", 2, {{1, 0.2}, {2, -0.3}}},
    {"sin_ncp60", "SIN", 2, {{1, 0.0}, {2, 0.57735026918962573}}},
    {"arc_default", "ARC", 0, {{0}}},
    {"zpn_linear", "ZPN", 2, {{0, 0.0}, {1, 1.0}}},
    {"zpn_quadratic", "ZPN", 3, {{0, 0.0}, {1, 1.0}, {2, -0.05}}},
    {"zpn_cubic", "ZPN", 4, {{0, 0.0}, {1, 1.0}, {2, 0.0}, {3, -0.02}}},
    {"zpn_cubic_inflect", "ZPN", 3, {{0, 0.0}, {1, 1.0}, {3, -0.2}}},
    {"zpn_quintic", "ZPN", 5, {{1, 0.9}, {2, 0.0}, {3, 0.1}, {4, 0.0}, {5, -0.02}}},
    /* The highest coefficient wcslib (and the Go code) reads: PV2_29. */
    {"zpn_pv29", "ZPN", 2, {{1, 1.0}, {29, -1e-12}}},
    {"zea_default", "ZEA", 0, {{0}}},
    {"air_default", "AIR", 0, {{0}}},
    {"air_thetab45", "AIR", 1, {{1, 45.0}}},
    /* ---- cylindrical ---- */
    {"cyp_default", "CYP", 0, {{0}}},
    {"cyp_mu0_lambda1", "CYP", 2, {{1, 0.0}, {2, 1.0}}},
    {"cyp_gall", "CYP", 2, {{1, 1.0}, {2, 0.70710678118654752}}},
    {"cyp_mu2_lambda0p5", "CYP", 2, {{1, 2.0}, {2, 0.5}}},
    {"cea_default", "CEA", 0, {{0}}},
    {"cea_lambda0p5", "CEA", 1, {{1, 0.5}}},
    {"car_default", "CAR", 0, {{0}}},
    {"mer_default", "MER", 0, {{0}}},
    /* ---- pseudo-cylindrical ---- */
    {"sfl_default", "SFL", 0, {{0}}},
    {"par_default", "PAR", 0, {{0}}},
    {"mol_default", "MOL", 0, {{0}}},
    {"ait_default", "AIT", 0, {{0}}},
    /* ---- conic ---- */
    {"cop_ta45", "COP", 1, {{1, 45.0}}},
    {"cop_ta45_eta10", "COP", 2, {{1, 45.0}, {2, 10.0}}},
    {"cop_tam30_eta15", "COP", 2, {{1, -30.0}, {2, 15.0}}},
    {"coe_ta45", "COE", 1, {{1, 45.0}}},
    {"coe_ta45_eta10", "COE", 2, {{1, 45.0}, {2, 10.0}}},
    {"coe_tam30_eta15", "COE", 2, {{1, -30.0}, {2, 15.0}}},
    {"cod_ta45", "COD", 1, {{1, 45.0}}},
    {"cod_ta45_eta10", "COD", 2, {{1, 45.0}, {2, 10.0}}},
    {"cod_tam30_eta15", "COD", 2, {{1, -30.0}, {2, 15.0}}},
    {"coo_ta45", "COO", 1, {{1, 45.0}}},
    {"coo_ta45_eta10", "COO", 2, {{1, 45.0}, {2, 10.0}}},
    {"coo_tam30_eta15", "COO", 2, {{1, -30.0}, {2, 15.0}}},
    /* ---- polyconic / pseudoconic ---- */
    {"bon_t30", "BON", 1, {{1, 30.0}}},
    {"bon_tm45", "BON", 1, {{1, -45.0}}},
    {"bon_t0", "BON", 1, {{1, 0.0}}},
    {"pco_default", "PCO", 0, {{0}}},
    /* ---- quad-cube ---- */
    {"tsc_default", "TSC", 0, {{0}}},
    {"csc_default", "CSC", 0, {{0}}},
    {"qsc_default", "QSC", 0, {{0}}},
    /* ---- HEALPix ---- */
    {"hpx_default", "HPX", 0, {{0}}},
    {"hpx_h4_k3", "HPX", 2, {{1, 4.0}, {2, 3.0}}},
    {"hpx_h3_k2", "HPX", 2, {{1, 3.0}, {2, 2.0}}},
    {"hpx_h6_k1", "HPX", 2, {{1, 6.0}, {2, 1.0}}},
    {"hpx_h5_k5", "HPX", 2, {{1, 5.0}, {2, 5.0}}},
    {"xph_default", "XPH", 0, {{0}}},
};
static const int ncases = (int)(sizeof(cases) / sizeof(cases[0]));

/* Shortest %.Ng (N = 15..17) that reads back as exactly v. */
static void jd(FILE *f, double v) {
  if (isnan(v) || isinf(v) || v == UNDEFINED) {
    fprintf(f, "null");
    return;
  }
  char buf[64];
  for (int prec = 15; prec <= 17; prec++) {
    snprintf(buf, sizeof(buf), "%.*g", prec, v);
    if (strtod(buf, NULL) == v) break;
  }
  fprintf(f, "%s", buf);
}

/* Big-endian encoders for the binary payload. */
static void put_u32(unsigned char *b, size_t *n, uint32_t v) {
  for (int i = 3; i >= 0; i--) b[(*n)++] = (unsigned char)(v >> (8 * i));
}

static void put_f64(unsigned char *b, size_t *n, double v) {
  uint64_t u;
  memcpy(&u, &v, sizeof u);
  for (int i = 7; i >= 0; i--) b[(*n)++] = (unsigned char)(u >> (8 * i));
}

/* Write dir/name = zlib (RFC 1950) compression of buf at level 9, via a
 * temp file renamed on success. zlib's output is deterministic for a given
 * zlib version and level; the version is recorded in the metadata. */
static int write_z(const char *dir, const char *name, const unsigned char *buf, size_t n) {
  uLongf zn = compressBound((uLong)n);
  unsigned char *z = malloc(zn);
  if (!z) return -1;
  if (compress2(z, &zn, buf, (uLong)n, 9) != Z_OK) {
    fprintf(stderr, "wref: zlib compress2 failed for %s\n", name);
    free(z);
    return -1;
  }
  char path[1024], tmp[1100];
  snprintf(path, sizeof(path), "%s/%s", dir, name);
  snprintf(tmp, sizeof(tmp), "%s.tmp", path);
  FILE *f = fopen(tmp, "wb");
  if (!f) {
    fprintf(stderr, "wref: open %s: %s\n", tmp, strerror(errno));
    free(z);
    return -1;
  }
  size_t w = fwrite(z, 1, zn, f);
  free(z);
  if (w != zn || fclose(f) != 0 || rename(tmp, path) != 0) {
    fprintf(stderr, "wref: write %s: %s\n", path, strerror(errno));
    return -1;
  }
  return 0;
}

/* Parameter sets every one of which wcslib's prjset REJECTS (or, for the
 * last few, accepts as a boundary case); written with their prjset status
 * to golden/prjset_status.json so the Go side can be required to refuse
 * exactly what wcslib refuses. */
static const case_t param_cases[] = {
    {"azp_mu_m1", "AZP", 1, {{1, -1.0}}},
    {"azp_gamma90", "AZP", 2, {{1, 2.0}, {2, 90.0}}},
    {"szp_w3_zero", "SZP", 3, {{1, -1.0}, {2, 0.0}, {3, 90.0}}},
    {"zpn_all_zero", "ZPN", 0, {{0}}},
    {"zpn_deg3_pv1_zero", "ZPN", 2, {{1, 0.0}, {3, 0.1}}},
    {"air_thetab_m90", "AIR", 1, {{1, -90.0}}},
    {"cyp_lambda0", "CYP", 2, {{1, 1.0}, {2, 0.0}}},
    {"cyp_mu_plus_lambda0", "CYP", 2, {{1, -1.0}, {2, 1.0}}},
    {"cea_lambda0", "CEA", 1, {{1, 0.0}}},
    {"cea_lambda1p5", "CEA", 1, {{1, 1.5}}},
    {"cop_no_pv1", "COP", 0, {{0}}},
    {"cop_ta0", "COP", 1, {{1, 0.0}}},
    {"coe_no_pv1", "COE", 0, {{0}}},
    {"cod_no_pv1", "COD", 0, {{0}}},
    {"coo_no_pv1", "COO", 0, {{0}}},
    {"bon_no_pv1", "BON", 0, {{0}}},
    {"hpx_h0", "HPX", 2, {{1, 0.0}, {2, 3.0}}},
    {"hpx_k_m1", "HPX", 2, {{1, 4.0}, {2, -1.0}}},
    /* accepted */
    {"bon_t0_ok", "BON", 1, {{1, 0.0}}},
    {"cea_lambda1_ok", "CEA", 1, {{1, 1.0}}},
    {"zpn_pv0_only_ok", "ZPN", 1, {{0, 0.5}}},
};
static const int nparam = (int)(sizeof(param_cases) / sizeof(param_cases[0]));

static int write_param_status(const char *outdir, const char *vers) {
  char path[1024], tmp[1100];
  snprintf(path, sizeof(path), "%s/prjset_status.json", outdir);
  snprintf(tmp, sizeof(tmp), "%s.tmp", path);
  FILE *f = fopen(tmp, "w");
  if (!f) return -1;
  fprintf(f, "{\n\"wcslib_version\": \"%s\",\n\"cases\": [\n", vers);
  for (int i = 0; i < nparam; i++) {
    const case_t *c = &param_cases[i];
    struct prjprm prj;
    prjini(&prj);
    strncpy(prj.code, c->code, 3);
    prj.code[3] = '\0';
    for (int k = 0; k < c->npv; k++) prj.pv[c->pv[k].m] = c->pv[k].v;
    int status = prjset(&prj);
    fprintf(f, "{\"name\": \"%s\", \"code\": \"%s\", \"pv\": [", c->name, c->code);
    for (int k = 0; k < c->npv; k++) {
      fprintf(f, "%s[%d, ", k ? ", " : "", c->pv[k].m);
      jd(f, c->pv[k].v);
      fprintf(f, "]");
    }
    fprintf(f, "], \"prjset_status\": %d}%s\n", status, i + 1 < nparam ? "," : "");
    prjfree(&prj);
  }
  fprintf(f, "]\n}\n");
  if (fclose(f) != 0 || rename(tmp, path) != 0) return -1;
  printf("  prjset_status.json: %d parameter sets\n", nparam);
  return 0;
}

#define MAXPT 16384

/* Extra probe points for branches the regular 15-degree / 50-unit grids
 * never reach:
 *  - QSC forward: wcslib qscs2x switches to a small-angle formula for
 *    1 - zeta < 1e-8, i.e. within ~0.008 deg of a cube-face centre. Points
 *    0.001-0.005 deg from all six face centres (faces 1-4 on the equator,
 *    0 and 5 at the poles), and the matching (x, y) near each face centre
 *    of the unfolded cube for the inverse. Also used for TSC and CSC.
 *  - MOL inverse: |y| = sqrt(2) r0 is the pole line; wcslib accepts it
 *    only at |x| <= tol and rejects any other x there (istat < 0 &&
 *    |x| > tol). Points on that line at x = 0, 0.5, 5 and just inside and
 *    outside it. */
static int extra_fwd(const char *code, double *phi, double *the) {
  static const double qsc[][2] = {
      {0.001, 0.002},   {90.003, -0.004},  {180.002, 0.001},
      {-89.996, 0.003}, {0.0, 89.997},     {10.0, -89.995},
      {0.0003, 0.0001}, {-0.004, -0.0005}, {135.0, 89.999},
  };
  int n = 0;
  if (!strcmp(code, "QSC") || !strcmp(code, "TSC") || !strcmp(code, "CSC")) {
    for (size_t k = 0; k < sizeof(qsc) / sizeof(qsc[0]); k++, n++) {
      phi[n] = qsc[k][0];
      the[n] = qsc[k][1];
    }
  }
  return n;
}

/* ---- Face-centre probes (TSC, CSC, QSC) ----
 *
 * qscs2x replaces zeco = 1 - zeta by a small-angle formula when
 * zeco < 1e-8, i.e. within ~0.0081 deg of a cube-face centre. The two
 * formulas agree there to ~5e-11 deg, so a moved switch is visible only
 * inside that ring and only at the last-bit level (round 2 of the proof
 * audit: the switch moved to 0.5e-8, and written <=, survived every test;
 * the face-centre points of extra_fwd sit at zeco ~ 1e-9 or below).
 * These points put zeco, as qscs2x itself computes it (sincosd, the same
 * products, the same face selection: qsc_zeco below), at
 * 1e-8 (1 +- {5e-4, 1e-3, 1e-2, 1e-1}), 0.5e-8 and 2e-8 on all six faces,
 * in 8 directions per face (20 + 45 k deg, so both |xi| > |eta| and
 * |eta| > |xi| and every sign occur). The distance from the face centre
 * is set by a few Newton-like rescalings of the angular offset (zeco
 * grows as its square); only the probe POSITION comes from this search,
 * the truth is wcslib's output there. A point whose computed zeco lands
 * on the wrong side of 1e-8 (impossible at these margins: zeco is exact,
 * a multiple of 2^-53) aborts the run. TSC and CSC get the same points
 * (no switch there, but more coverage next to the face centres). */
static double qsc_zeco(double phi, double the) {
  double sinphi, cosphi, sinthe, costhe;
  sincosd(phi, &sinphi, &cosphi);
  sincosd(the, &sinthe, &costhe);
  double l = costhe * cosphi, m = costhe * sinphi, n = sinthe;
  double zeta = n;
  if (l > zeta) zeta = l;
  if (m > zeta) zeta = m;
  if (-l > zeta) zeta = -l;
  if (-m > zeta) zeta = -m;
  if (-n > zeta) zeta = -n;
  return 1.0 - zeta;
}

static const double centre_z[] = {0.5,    0.9,    0.99,   0.999,  0.9995,
                                  1.0005, 1.001,  1.01,   1.1,    2.0};

static int centre_fwd(const char *code, double *phi, double *the) {
  if (strcmp(code, "QSC") && strcmp(code, "TSC") && strcmp(code, "CSC")) return 0;
  const double r2d = 180.0 / 3.14159265358979323846;
  const double d2r = 3.14159265358979323846 / 180.0;
  /* face centres (phi, theta), faces 0..5 */
  static const double fc[6][2] = {{0.0, 90.0},  {0.0, 0.0},   {90.0, 0.0},
                                  {180.0, 0.0}, {-90.0, 0.0}, {0.0, -90.0}};
  int n = 0;
  for (int f = 0; f < 6; f++)
    for (int k = 0; k < 8; k++) {
      double a = (20.0 + 45.0 * k) * d2r;
      for (size_t iz = 0; iz < sizeof(centre_z) / sizeof(centre_z[0]); iz++) {
        double z = 1e-8 * centre_z[iz];
        double d = 2.0 * asin(sqrt(z / 2.0)) * r2d; /* 1 - cos d = z */
        double p = 0.0, t = 0.0;
        for (int it = 0; it < 6; it++) {
          if (f == 0 || f == 5) {
            p = 20.0 + 45.0 * k;
            if (p > 180.0) p -= 360.0;
            t = fc[f][1] > 0 ? 90.0 - d : -90.0 + d;
          } else {
            t = d * sin(a);
            p = fc[f][0] + d * cos(a) / cos(t * d2r);
            if (p > 180.0) p -= 360.0;
          }
          double zc = qsc_zeco(p, t);
          if (zc <= 0.0) break;
          d *= sqrt(z / zc);
        }
        double zc = qsc_zeco(p, t);
        if ((zc < 1e-8) != (centre_z[iz] < 1.0) || fabs(zc / z - 1.0) > 1e-4) {
          fprintf(stderr, "wref: %s: face-centre probe f=%d k=%d z=%g landed at zeco=%.17g\n",
                  code, f, k, z, zc);
          exit(1);
        }
        phi[n] = p;
        the[n] = t;
        n++;
      }
    }
  return n;
}

static int extra_inv(const char *code, double *x, double *y) {
  static const double cube[][2] = {
      {0.001, 0.002},  {90.003, -0.004}, {180.002, 0.001}, {270.004, 0.003},
      {0.0, 89.997},   {0.001, -89.995}, {-0.0005, 0.0003},
  };
  /* sqrt(2) * 180/pi, the MOL pole line */
  const double yp = 81.028468454139556;
  const double mol[][2] = {
      {0.0, yp},  {0.5, yp},   {5.0, yp},  {-3.0, yp},  {0.0, -yp},
      {2.0, -yp}, {1e-13, yp}, {0.5, yp - 1e-9}, {0.5, yp + 1e-9},
  };
  /* AIR far out: airx2s brackets the root by halving xi's cosine 30
   * times, which reaches r ~ 5e8 deg; past that it rejects (k == 30), and
   * its regula falsi would give up too (k == 100). Radii either side, on
   * two directions, so both decisions are pinned against wcslib. */
  static const double airr[] = {1e6, 1e8, 3e8, 4e8, 5e8, 6e8, 1e9, 1e10, 1e12, 1e15};
  int n = 0;
  if (!strcmp(code, "AIR")) {
    for (size_t k = 0; k < sizeof(airr) / sizeof(airr[0]); k++) {
      x[n] = airr[k];
      y[n] = 0.0;
      n++;
      x[n] = -0.6 * airr[k];
      y[n] = 0.8 * airr[k];
      n++;
    }
  } else if (!strcmp(code, "QSC") || !strcmp(code, "TSC") || !strcmp(code, "CSC")) {
    for (size_t k = 0; k < sizeof(cube) / sizeof(cube[0]); k++, n++) {
      x[n] = cube[k][0];
      y[n] = cube[k][1];
    }
  } else if (!strcmp(code, "MOL")) {
    for (size_t k = 0; k < sizeof(mol) / sizeof(mol[0]); k++, n++) {
      x[n] = mol[k][0];
      y[n] = mol[k][1];
    }
  }
  return n;
}

/* ---- Small-angle band probes (generic, the same for every case) ----
 *
 * Several wcslib projections switch to a small-angle approximation below
 * a fixed threshold (PCO: |y| < 1e-4 deg in pcox2s and |theta| < 1e-4
 * deg in pcos2x; XPH: sigma < 1e-4). A port that moves such a switch
 * changes results only inside a narrow band next to the equator (or
 * pole), which the 15-deg grid never visits, so a wrong switch would go
 * unseen (the verifier widened PCO's two switches to 1e-1 and 1e-3 deg
 * and both survived). These points straddle the band on both sides of
 * every such threshold, off the central meridian (where the
 * approximation is exact):
 *   inverse: y = +-{5e-5, 9.9e-5, 9.95e-5, 9.99e-5, 1e-4, 1.005e-4,
 *            1.01e-4, 2e-4, 1e-3, 5e-3, 2e-2, 5e-2, 9e-2} deg at
 *            x = -170, -45, 10, 170 deg;
 *   forward: theta = +-{5e-5, 9.9e-5, 9.95e-5, 9.99e-5, 1e-4, 1.005e-4,
 *            1.01e-4, 2e-4, 5e-4, 9e-4, 3e-3, 9e-3, 3e-2} deg at
 *            phi = -170, -60, 30, 170 deg.
 * The pairs 9.9e-5 / 1.01e-4 pin a 1e-4 switch to 1%; 1e-4 itself and
 * 1.005e-4 pin the comparison (< against <=) and the switch to 0.5% from
 * above; 9.95e-5 and 9.99e-5 (added in round 2 of the proof audit, when
 * the verifier moved both PCO switches 0.5% and 0.9% DOWN and they
 * survived) pin it to 0.5% and 0.1% from below. */
static const double band_y[] = {5e-5, 9.9e-5, 9.95e-5, 9.99e-5, 1e-4, 1.005e-4, 1.01e-4, 2e-4, 1e-3, 5e-3, 2e-2, 5e-2, 9e-2};
static const double band_x[] = {-170.0, -45.0, 10.0, 170.0};
static const double band_t[] = {5e-5, 9.9e-5, 9.95e-5, 9.99e-5, 1e-4, 1.005e-4, 1.01e-4, 2e-4, 5e-4, 9e-4, 3e-3, 9e-3, 3e-2};
static const double band_p[] = {-170.0, -60.0, 30.0, 170.0};

static int band_inv(double *x, double *y) {
  int n = 0;
  for (size_t j = 0; j < sizeof(band_y) / sizeof(band_y[0]); j++)
    for (int s = -1; s <= 1; s += 2)
      for (size_t i = 0; i < sizeof(band_x) / sizeof(band_x[0]); i++, n++) {
        x[n] = band_x[i];
        y[n] = s * band_y[j];
      }
  return n;
}

static int band_fwd(double *phi, double *the) {
  int n = 0;
  for (size_t j = 0; j < sizeof(band_t) / sizeof(band_t[0]); j++)
    for (int s = -1; s <= 1; s += 2)
      for (size_t i = 0; i < sizeof(band_p) / sizeof(band_p[0]); i++, n++) {
        phi[n] = band_p[i];
        the[n] = s * band_t[j];
      }
  return n;
}

/* ---- Pole band probes (every case whose fiducial point is the native
 * pole, prjset theta0 = 90: the zenithal codes and XPH) ----
 *
 * The zenithal projections and XPH switch to a small-angle formula next
 * to the pole, at thresholds the regular grid and the 90 - 10^-k probes
 * (1e-4, 1e-3 deg) do not straddle:
 *   SIN forward: 90 - |theta| < 1e-5 rad (5.73e-4 deg);
 *   SIN (slant) and SZP inverse: x^2 + y^2 < 1e-10 rad^2 (r < 5.73e-4 deg);
 *   AIR forward xi < 1e-4 rad and inverse r < w2 1e-4, the same
 *     xi = (90 - theta)/2, i.e. 90 - theta < 1.146e-2 deg;
 *   XPH forward theta < w[5] and inverse sigma < 1e-4, the same sigma,
 *     i.e. 90 - theta < 4.68e-3 deg.
 * A port that moves one of them (the verifier moved the SIN/SZP switch
 * to 1e-11 and 0.5e-10, AIR's to 2e-4 and XPH's to 1.2e-4 .. 2e-4) gives
 * results 1e-10 .. 5e-9 deg from wcslib only inside that annulus. These
 * points cover it, at both poles, off every axis and facet edge:
 *   theta = +-(90 - d) at phi = -150.3, -41.7, 22.9, 113.1 deg, with
 *   d = 1e-5 10^(k/10) deg, k = 0..40 (1e-5 .. 0.1 deg, 26% apart), and
 *   d = T (1 +- 5e-4), T (1 +- 1e-3), T (1 +- 1e-2) for T = 1e-5 rad,
 *   2e-4 rad and 1e-4 sqrt(2/3) rad (in degrees), plus theta = w[5]
 *   exactly for XPH. (The 5e-4 factors were added in round 2 of the proof
 *   audit: at T (1 + 1e-3) SIN's computed t is not below 1.001e-5, so a
 *   switch moved exactly 0.1% up survived.)
 * Their forward images join the inverse set as round trips, so the
 * inverse switches (at the same xi, sigma or r) are straddled too. */
static const double pole_p[] = {-150.3, -41.7, 22.9, 113.1};

static int pole_fwd(const struct prjprm *prj, double *phi, double *the) {
  const double r2d = 180.0 / 3.14159265358979323846;
  const double thr[] = {1e-5 * r2d, 2e-4 * r2d, 1e-4 * sqrt(2.0 / 3.0) * r2d};
  const double fac[] = {1.0 - 1e-3, 1.0 + 1e-3, 1.0 - 1e-2, 1.0 + 1e-2, 1.0 - 5e-4, 1.0 + 5e-4};
  double d[64];
  int nd = 0;
  for (int k = 0; k <= 40; k++) d[nd++] = 1e-5 * pow(10.0, k / 10.0);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 6; j++) d[nd++] = thr[i] * fac[j];
  int n = 0;
  for (int k = 0; k < nd; k++)
    for (int s = -1; s <= 1; s += 2)
      for (size_t i = 0; i < sizeof(pole_p) / sizeof(pole_p[0]); i++, n++) {
        phi[n] = pole_p[i];
        the[n] = s * (90.0 - d[k]);
      }
  if (!strcmp(prj->code, "XPH"))
    for (size_t i = 0; i < sizeof(pole_p) / sizeof(pole_p[0]); i++, n++) {
      phi[n] = pole_p[i];
      the[n] = prj->w[5];
    }
  return n;
}

/* ---- Tie probes: inputs ON an inverse small-angle threshold ----
 *
 * sinx2s and szpx2s switch on r2 < 1e-10 with r2 = x0 x0 + y0 y0,
 * x0 = x w[0], y0 = y w[0]; airx2s on r < w[5] with
 * r = sqrt(x x + y y) / w[0]. A port that writes <= for < differs only
 * where the computed value equals the threshold exactly. These inputs make
 * it so: y is fixed (1e-4 and 2e-4 deg for SIN/SZP, 0 and 4e-3 deg for
 * AIR; at y = 0 no x makes fl(x0 x0) exactly 1e-10) and x
 * is stepped one ULP at a time from the threshold's nominal x until
 * wcslib's own expression for r2 (or r), evaluated here in the same order
 * with wcslib's prjprm w[], equals the threshold; the tie point and the
 * doubles either side of it in x are probes. Only the probe position
 * comes from this search; the truth is wcslib's output there. (The
 * forward switches, SIN's t < 1e-5 and AIR's xi < w[4], and XPH's inverse
 * sigma < 1e-4 cannot be tied at all: TestSmallAngleSwitchTies.) */
static int tie_inv(const struct prjprm *prj, double *x, double *y) {
  int air = !strcmp(prj->code, "AIR");
  if (!air && strcmp(prj->code, "SIN") && strcmp(prj->code, "SZP")) return 0;
  const double ys[2] = {air ? 0.0 : 1e-4, air ? 4e-3 : 2e-4};
  int n = 0;
  for (int iy = 0; iy < 2; iy++) {
    volatile double yv = ys[iy];
    double x0;
    if (air) {
      /* r = w[5]: x = sqrt((w[5] w[0])^2 - y^2) */
      double rr = prj->w[5] * prj->w[0];
      x0 = sqrt(rr * rr - yv * yv);
    } else {
      x0 = sqrt(1e-10 / (prj->w[0] * prj->w[0]) - yv * yv);
    }
    int found = 0;
    volatile double xu = x0, xd = x0;
    for (int it = 0; it < 100000 && !found; it++) {
      for (int k = 0; k < 2 && !found; k++) {
        volatile double xv = k ? xd : xu;
        int tie;
        if (air) {
          volatile double yj2 = yv * yv;
          volatile double s = xv * xv + yj2;
          volatile double r = sqrt(s) / prj->w[0];
          tie = r == prj->w[5];
        } else {
          volatile double a = xv * prj->w[0];
          volatile double b = yv * prj->w[0];
          volatile double b2 = b * b;
          volatile double r2 = a * a + b2;
          tie = r2 == 1.0e-10;
        }
        if (tie) {
          x[n] = nextafter(xv, 0.0);
          y[n++] = yv;
          x[n] = xv;
          y[n++] = yv;
          x[n] = nextafter(xv, 1.0);
          y[n++] = yv;
          found = 1;
        }
      }
      xu = nextafter(xu, 1.0);
      xd = nextafter(xd, 0.0);
    }
    if (!found) fprintf(stderr, "wref: %s: no tie probe found at y = %g\n", prj->code, yv);
  }
  return n;
}

/* ---- XPH seam probes: inputs whose native longitude leaves [-180, 180]
 * by rounding ----
 *
 * Next to XPH's north pole on the half-line x = 0, y > 0 (the phi = +-180
 * seam), xphx2s forms phi = +-180 + (45 + xi1 / sigma) - 45 with sigma
 * = (90 - |eta|) / 45 small, so phi lands a few 1e-14 .. 1e-7 deg past
 * +-180, and the final prjbchk(1e-12) decides: within 1e-12 it clamps to
 * +-180 and accepts, further out it rejects. No grid or boundary probe
 * reached that window (round 2 of the proof audit: a prjbchk tolerance of
 * 1e-6 survived every test; a 2e7-point random search with wcslib found
 * 47806 such inputs, all within ~1e-5 deg of the pole). These probes are
 * found the same way, with wcslib's own xphx2s with the bounds check
 * switched off: on y = u for u = 1e-6 .. 0.1 and x = +-k 1e-17 u
 * (k = 1, 2, ...), the first x whose raw phi is past +-180 by at most
 * 1e-12 (clamped) and the first by 1e-12 .. 1e-6 (rejected). Only the
 * probe POSITION comes from this search; the truth is wcslib's output
 * there, bounds check on. */
static int xph_seam_inv(const struct prjprm *prj, double *x, double *y) {
  if (strcmp(prj->code, "XPH")) return 0;
  static const double us[] = {1e-6, 3e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1};
  struct prjprm q = *prj;
  q.bounds &= ~4;
  int n = 0;
  for (size_t i = 0; i < sizeof(us) / sizeof(us[0]); i++)
    for (int s = -1; s <= 1; s += 2) {
      int fc = 0, fw = 0;
      for (int k = 1; k <= 200000 && !(fc && fw); k++) {
        double xv = s * k * 1e-17 * us[i], yv = us[i], ph, th;
        int st = 0;
        q.prjx2s(&q, 1, 1, 1, 1, &xv, &yv, &ph, &th, &st);
        if (st) continue;
        double e = fabs(ph) - 180.0;
        if (!fc && e > 0.0 && e <= 1e-12) {
          fc = 1;
          x[n] = xv;
          y[n++] = yv;
        }
        if (!fw && e > 1e-12 && e <= 1e-6) {
          fw = 1;
          x[n] = xv;
          y[n++] = yv;
        }
      }
    }
  return n;
}

/* ---- Boundary probes (generic, the same for every case) ----
 *
 * The regular grids are 15-50 deg apart, so on their own they prove only
 * that a point well off the projection is refused; a boundary placed a
 * few percent too far out would pass. These probes pin every boundary to
 * within BDELTA:
 *
 *  1. Scan lines. Inverse: 28 horizontal and 28 vertical lines of the
 *     (x, y) plane (y, x = -397.3 + 30 k, off the regular grid) and 24 rays
 *     from the origin (angle 1.3 + 15 k deg), plus the 9 + 9 exact lines
 *     x, y = 45 k (k = -4..4: the axes, facet edges, PAR's pole lines),
 *     each sampled every 0.5 deg over [-420, 420]. Forward: 24 parallels (theta = -87.9 + 7.5 k) over
 *     phi in [-180, 180] and 24 meridians (phi = -175.7 + 15 k) over theta
 *     in [-90, 90], sampled every 0.25 deg. A point counts as accepted when
 *     wcslib's stat is 0 and its outputs are finite (the test's rule).
 *     Wherever two neighbouring samples differ, the crossing is bisected
 *     (80 halvings) and the two points BDELTA either side of it are probes,
 *     if the crossing is transversal (see transversal()).
 *  2. Pole images. For phi in {0, 45, 90, 180, -90} and theta = +-90,
 *     every finite forward image (x_p, y_p) of the pole, and the points
 *     d (cos a, sin a) from it for d in {1e-9, 1e-7, 1e-5, 1e-3} deg and
 *     a = 0, 45, ..., 315 deg: degenerate boundaries (a pole that is a
 *     single point of a line, PAR/MOL/SFL...) and near-pole small-angle
 *     branches (XPH sigma < 1e-4, QSC) are reached here.
 *  3. Forward near-pole and seam points: theta = +-(90 - 10^-k), k = 1..7,
 *     at phi in {0, 30, 90, 135, -150, 179.9}; phi = +-(180 -+ BDELTA) at
 *     theta in {0, 30, -60, 80}.
 *
 * Every accepted forward probe's (x, y) also joins the inverse set (the
 * round trip, as for the regular grid).
 *
 * BDELTA = 1e-7 deg: about 1e6 times the ULP of a coordinate <= 420, so a
 * correct port is never on the other side by rounding, and far inside the
 * smallest seeded boundary error the verifier used (a CEA edge 1% out,
 * PAR's pole-line tolerance at 1 deg, XPH bounds 1 deg out) and below
 * wcslib's coarsest tolerance (CSC 1e-7 in face units is 4.5e-6 deg). */
#define BDELTA 1e-7

static int accept_inv(struct prjprm *prj, double x, double y) {
  double p, t;
  int st = 0;
  prj->prjx2s(prj, 1, 1, 1, 1, &x, &y, &p, &t, &st);
  return st == 0 && isfinite(p) && isfinite(t);
}

static int accept_fwd(struct prjprm *prj, double p, double t) {
  double x, y;
  int st = 0;
  prj->prjs2x(prj, 1, 1, 1, 1, &p, &t, &x, &y, &st);
  return st == 0 && isfinite(x) && isfinite(y);
}

/* A crossing at o + sc u is used only if it is transversal: shifting the
 * line sideways by BSHIFT (1e-9 deg, either way) must leave a crossing
 * within BWIN (1e-6 deg) of sc, i.e. the boundary meets the line at more
 * than about 0.06 deg. A boundary that only grazes the line (PCO's
 * phi = 180 + tol edge along x = 180, where phi - 180 grows as y^2) is a
 * knife edge: there points 1e-7 deg along the line differ in the deciding
 * quantity by 1e-14, below the last-bit noise of any port. Grazing
 * crossings are skipped; the transversal ones pin every boundary. The
 * exception is an accepted piece that is itself only a line (see below). */
#define BSHIFT 1e-9
#define BWIN 1e-6
static int transversal(struct prjprm *prj, int (*acc)(struct prjprm *, double, double),
                       double ox, double oy, double ux, double uy, double sc) {
  int nsame = 0, nrej = 0;
  for (int sgn = -1; sgn <= 1; sgn += 2) {
    double px = ox - sgn * BSHIFT * uy, py = oy + sgn * BSHIFT * ux;
    int a = acc(prj, px + (sc - BWIN) * ux, py + (sc - BWIN) * uy);
    int b = acc(prj, px + (sc + BWIN) * ux, py + (sc + BWIN) * uy);
    if (a == b) {
      nsame++;
      nrej += !a;
    }
  }
  /* Both shifted lines entirely off the projection: the accepted part of
   * this line is a lower-dimensional piece of it (PAR's pole line
   * y = +-90, accepted only at x = 0), decided by exact arithmetic on the
   * exact line, so its crossings are kept. */
  return nsame == 0 || nrej == 2;
}

/* Scan the line o + s u, s in [s0, s1] step ds; append the BDELTA probes
 * either side of every accept/reject crossing to (a[], b[]). */
static int scan_line(struct prjprm *prj, int inverse, double ox, double oy,
                     double ux, double uy, double s0, double s1, double ds,
                     double *a, double *b, int n, int nmax) {
  int (*acc)(struct prjprm *, double, double) = inverse ? accept_inv : accept_fwd;
  int prev = acc(prj, ox + s0 * ux, oy + s0 * uy);
  double sp = s0;
  int ns = (int)((s1 - s0) / ds + 0.5);
  for (int k = 1; k <= ns; k++) {
    double s = s0 + ds * k;
    int cur = acc(prj, ox + s * ux, oy + s * uy);
    if (cur != prev) {
      double lo = sp, hi = s;
      for (int it = 0; it < 80; it++) {
        double mid = 0.5 * (lo + hi);
        if (acc(prj, ox + mid * ux, oy + mid * uy) == prev)
          lo = mid;
        else
          hi = mid;
      }
      double sc = 0.5 * (lo + hi);
      if (!transversal(prj, acc, ox, oy, ux, uy, sc)) {
        prev = cur;
        sp = s;
        continue;
      }
      if (n + 2 > nmax) return -1;
      a[n] = ox + (sc - BDELTA) * ux;
      b[n] = oy + (sc - BDELTA) * uy;
      n++;
      a[n] = ox + (sc + BDELTA) * ux;
      b[n] = oy + (sc + BDELTA) * uy;
      n++;
    }
    prev = cur;
    sp = s;
  }
  return n;
}

static int boundary_inv(struct prjprm *prj, double *x, double *y, int nmax) {
  const double d2r = 3.14159265358979323846 / 180.0;
  int n = 0;
  for (int k = 0; k < 28 && n >= 0; k++) {
    double c = -397.3 + 30.0 * k;
    n = scan_line(prj, 1, -420.0, c, 1.0, 0.0, 0.0, 840.0, 0.5, x, y, n, nmax);
    if (n < 0) break;
    n = scan_line(prj, 1, c, -420.0, 0.0, 1.0, 0.0, 840.0, 0.5, x, y, n, nmax);
  }
  /* Exact lines at multiples of 45 deg (the axes, HPX/XPH facet edges, a
   * pole line such as PAR's y = +-90, on which wcslib accepts only x = 0). */
  for (int k = -4; k <= 4 && n >= 0; k++) {
    n = scan_line(prj, 1, -420.0, 45.0 * k, 1.0, 0.0, 0.0, 840.0, 0.5, x, y, n, nmax);
    if (n < 0) break;
    n = scan_line(prj, 1, 45.0 * k, -420.0, 0.0, 1.0, 0.0, 840.0, 0.5, x, y, n, nmax);
  }
  for (int k = 0; k < 24 && n >= 0; k++) {
    double a = (1.3 + 15.0 * k) * d2r;
    n = scan_line(prj, 1, 0.0, 0.0, cos(a), sin(a), 0.0, 420.0, 0.5, x, y, n, nmax);
  }
  if (n < 0) return -1;
  /* Pole images and the points around them. */
  static const double pphi[] = {0.0, 45.0, 90.0, 180.0, -90.0};
  static const double dd[] = {1e-9, 1e-7, 1e-5, 1e-3};
  double seen[20][2];
  int nseen = 0;
  for (int ip = 0; ip < 5; ip++)
    for (int pole = -1; pole <= 1; pole += 2) {
      double p = pphi[ip], t = 90.0 * pole, xp, yp;
      int st = 0;
      prj->prjs2x(prj, 1, 1, 1, 1, &p, &t, &xp, &yp, &st);
      if (st != 0 || !isfinite(xp) || !isfinite(yp)) continue;
      int dup = 0;
      for (int s = 0; s < nseen; s++)
        if (seen[s][0] == xp && seen[s][1] == yp) dup = 1;
      if (dup) continue;
      seen[nseen][0] = xp;
      seen[nseen][1] = yp;
      nseen++;
      if (n + 1 + 4 * 8 > nmax) return -1;
      x[n] = xp;
      y[n] = yp;
      n++;
      for (int id = 0; id < 4; id++)
        for (int ia = 0; ia < 8; ia++) {
          double a = 45.0 * ia * d2r;
          x[n] = xp + dd[id] * cos(a);
          y[n] = yp + dd[id] * sin(a);
          n++;
        }
    }
  return n;
}

static int boundary_fwd(struct prjprm *prj, double *phi, double *the, int nmax) {
  int n = 0;
  for (int k = 0; k < 24 && n >= 0; k++) {
    double t = -87.9 + 7.5 * k;
    n = scan_line(prj, 0, -180.0, t, 1.0, 0.0, 0.0, 360.0, 0.25, phi, the, n, nmax);
    if (n < 0) break;
    double p = -175.7 + 15.0 * k;
    n = scan_line(prj, 0, p, -90.0, 0.0, 1.0, 0.0, 180.0, 0.25, phi, the, n, nmax);
  }
  if (n < 0) return -1;
  /* Not phi = +-180 itself: next to a pole, the inverse of a point exactly
   * on the phi = +-180 seam is decided by the last digits of the
   * projection constants (COO: four ULP of w[2] move phi by 1.8e-9 deg
   * across 180 at theta = 89.99999), a knife edge no port can be held
   * to. The seam is probed BDELTA either side below, and on it at the
   * regular grid's theta. */
  static const double pp[] = {0.0, 30.0, 90.0, 135.0, -150.0, 179.9};
  for (int k = 1; k <= 7; k++)
    for (int ip = 0; ip < 6; ip++)
      for (int pole = -1; pole <= 1; pole += 2) {
        if (n >= nmax) return -1;
        phi[n] = pp[ip];
        the[n] = pole * (90.0 - pow(10.0, -k));
        n++;
      }
  static const double st[] = {0.0, 30.0, -60.0, 80.0};
  for (int it = 0; it < 4; it++) {
    const double ps[] = {180.0 - BDELTA, 180.0 + BDELTA, -180.0 + BDELTA,
                         -180.0 - BDELTA};
    for (int ip = 0; ip < 4; ip++) {
      if (n >= nmax) return -1;
      phi[n] = ps[ip];
      the[n] = st[it];
      n++;
    }
  }
  return n;
}

static double fphi[MAXPT], fthe[MAXPT], fx[MAXPT], fy[MAXPT];
static int fstat[MAXPT];
static double ix[MAXPT], iy[MAXPT], iphi[MAXPT], ithe[MAXPT];
static int istat[MAXPT];

static int run_case(const case_t *c, const char *outdir, const char *vers) {
  struct prjprm prj;
  prjini(&prj);
  strncpy(prj.code, c->code, 3);
  prj.code[3] = '\0';
  for (int i = 0; i < c->npv; i++) prj.pv[c->pv[i].m] = c->pv[i].v;

  int status = prjset(&prj);
  if (status != 0) {
    fprintf(stderr, "wref: %s: prjset status=%d (%s)\n", c->name, status,
            prj_errmsg[status]);
    return -1;
  }

  /* Forward grid. */
  int nf = 0;
  for (int it = 0; it <= 12; it++) {
    double the = -90.0 + 15.0 * it;
    for (int ip = 0; ip <= 26; ip++) {
      double phi;
      if (ip <= 24)
        phi = -180.0 + 15.0 * ip;
      else if (ip == 25)
        phi = 270.0;
      else
        phi = -270.0;
      fphi[nf] = phi;
      fthe[nf] = the;
      fstat[nf] = 0;
      nf++;
    }
  }
  /* Code-specific forward probes for branches the regular grid cannot
   * reach (see extra_fwd). */
  nf += extra_fwd(c->code, fphi + nf, fthe + nf);
  nf += centre_fwd(c->code, fphi + nf, fthe + nf);
  nf += band_fwd(fphi + nf, fthe + nf);
  if (prj.theta0 == 90.0) nf += pole_fwd(&prj, fphi + nf, fthe + nf);
  {
    int nb = boundary_fwd(&prj, fphi + nf, fthe + nf, MAXPT / 2 - nf);
    if (nb < 0) {
      fprintf(stderr, "wref: %s: too many forward boundary probes\n", c->name);
      return -1;
    }
    nf += nb;
  }
  for (int k = 0; k < nf; k++) fstat[k] = 0;
  prj.prjs2x(&prj, nf, 0, 1, 1, fphi, fthe, fx, fy, fstat);

  /* Inverse grid. */
  int ni = 0;
  for (int j = 0; j <= 16; j++)
    for (int i = 0; i <= 16; i++) {
      ix[ni] = -400.0 + 50.0 * i;
      iy[ni] = -400.0 + 50.0 * j;
      ni++;
    }
  for (int j = 0; j <= 24; j++)
    for (int i = 0; i <= 24; i++) {
      ix[ni] = -180.0 + 15.0 * i;
      iy[ni] = -180.0 + 15.0 * j;
      ni++;
    }
  ni += extra_inv(c->code, ix + ni, iy + ni);
  ni += band_inv(ix + ni, iy + ni);
  ni += tie_inv(&prj, ix + ni, iy + ni);
  ni += xph_seam_inv(&prj, ix + ni, iy + ni);
  {
    int nb = boundary_inv(&prj, ix + ni, iy + ni, MAXPT / 2 - ni);
    if (nb < 0) {
      fprintf(stderr, "wref: %s: too many inverse boundary probes\n", c->name);
      return -1;
    }
    ni += nb;
  }
  /* Everything after ni_own is a round trip: the (x, y) of an accepted,
   * finite forward row, in forward order. The golden does not store these
   * inputs again; the Go test rebuilds them from the forward rows. */
  const int ni_own = ni;
  for (int k = 0; k < nf; k++) {
    if (fstat[k] != 0 || !isfinite(fx[k]) || !isfinite(fy[k])) continue;
    if (ni >= MAXPT) {
      fprintf(stderr, "wref: %s: too many points\n", c->name);
      return -1;
    }
    ix[ni] = fx[k];
    iy[ni] = fy[k];
    ni++;
  }
  for (int k = 0; k < ni; k++) istat[k] = 0;
  prj.prjx2s(&prj, ni, 0, 1, 1, ix, iy, iphi, ithe, istat);

  /* Payload (see the header comment), then zlib. */
  size_t cap = 12 + (size_t)nf * 33 + (size_t)ni * 33;
  unsigned char *buf = malloc(cap);
  if (!buf) return -1;
  size_t n = 0;
  put_u32(buf, &n, (uint32_t)nf);
  put_u32(buf, &n, (uint32_t)ni_own);
  put_u32(buf, &n, (uint32_t)(ni - ni_own));
  for (int k = 0; k < nf; k++) {
    put_f64(buf, &n, fphi[k]);
    put_f64(buf, &n, fthe[k]);
  }
  for (int k = 0; k < nf; k++) {
    if (fstat[k] < 0 || fstat[k] > 255) {
      fprintf(stderr, "wref: %s: forward stat %d does not fit a byte\n", c->name, fstat[k]);
      return -1;
    }
    buf[n++] = (unsigned char)fstat[k];
  }
  for (int k = 0; k < nf; k++) {
    if (fstat[k]) continue;
    put_f64(buf, &n, fx[k]);
    put_f64(buf, &n, fy[k]);
  }
  for (int k = 0; k < ni_own; k++) {
    put_f64(buf, &n, ix[k]);
    put_f64(buf, &n, iy[k]);
  }
  for (int k = 0; k < ni; k++) {
    if (istat[k] < 0 || istat[k] > 255) {
      fprintf(stderr, "wref: %s: inverse stat %d does not fit a byte\n", c->name, istat[k]);
      return -1;
    }
    buf[n++] = (unsigned char)istat[k];
  }
  for (int k = 0; k < ni; k++) {
    if (istat[k]) continue;
    put_f64(buf, &n, iphi[k]);
    put_f64(buf, &n, ithe[k]);
  }
  char binname[256];
  snprintf(binname, sizeof(binname), "%s.bin", c->name);
  if (write_z(outdir, binname, buf, n) != 0) {
    free(buf);
    return -1;
  }
  free(buf);

  char path[1024], tmp[1100];
  snprintf(path, sizeof(path), "%s/%s.json", outdir, c->name);
  snprintf(tmp, sizeof(tmp), "%s.tmp", path);
  FILE *f = fopen(tmp, "w");
  if (!f) {
    fprintf(stderr, "wref: open %s: %s\n", tmp, strerror(errno));
    return -1;
  }
  fprintf(f, "{\"name\": \"%s\", \"code\": \"%s\", \"wcslib_version\": \"%s\", \"zlib_version\": \"%s\",\n",
          c->name, c->code, vers, zlibVersion());
  fprintf(f, "\"pv\": [");
  for (int i = 0; i < c->npv; i++) {
    fprintf(f, "%s[%d, ", i ? ", " : "", c->pv[i].m);
    jd(f, c->pv[i].v);
    fprintf(f, "]");
  }
  fprintf(f, "],\n\"prjset\": {\"phi0\": ");
  jd(f, prj.phi0);
  fprintf(f, ", \"theta0\": ");
  jd(f, prj.theta0);
  fprintf(f, ", \"x0\": ");
  jd(f, prj.x0);
  fprintf(f, ", \"y0\": ");
  jd(f, prj.y0);
  fprintf(f, ", \"r0\": ");
  jd(f, prj.r0);
  fprintf(f, ", \"pv_after\": [");
  for (int i = 0; i < 4; i++) {
    if (i) fprintf(f, ", ");
    jd(f, prj.pv[i]);
  }
  fprintf(f, "]},\n\"forward\": %d, \"inverse\": %d, \"roundtrip\": %d, \"payload\": \"%s\", \"payload_bytes\": %zu}\n",
          nf, ni_own, ni - ni_own, binname, n);
  if (fclose(f) != 0 || rename(tmp, path) != 0) {
    fprintf(stderr, "wref: write %s: %s\n", path, strerror(errno));
    return -1;
  }
  int nfb = 0, nib = 0;
  for (int k = 0; k < nf; k++) nfb += fstat[k] != 0;
  for (int k = 0; k < ni; k++) nib += istat[k] != 0;
  printf("  %-28s %s fwd %4d (%3d rejected)  inv %4d (%4d rejected)\n",
         c->name, c->code, nf, nfb, ni, nib);
  prjfree(&prj);
  return 0;
}

/* wcstrig guards: asind/acosd at arguments just inside, on and outside
 * +-(1 + WCSTRIG_TOL) (wcstrig.h: 1e-10). The projections cannot reach
 * these arguments (a dense scan found none: every caller's argument is
 * within a few ULP of [-1, 1] or far outside it), so the tolerance is
 * pinned here, function by function, against wcslib's own asind/acosd.
 * Written to <outdir>/wcstrig_guards.json: [v, asind(v), acosd(v)] with
 * null for NaN. */
static int write_trig(const char *outdir, const char *vers) {
  static const double d[] = {0.0,    1e-16, 2.3e-16, 1e-12, 5e-11, 9.9e-11, 9.99e-11, 1e-10,
                             1.001e-10, 1.01e-10, 2e-10, 1e-9, 1e-6, 1e-3, 0.5};
  char path[1024], tmp[1100];
  snprintf(path, sizeof(path), "%s/wcstrig_guards.json", outdir);
  snprintf(tmp, sizeof(tmp), "%s.tmp", path);
  FILE *f = fopen(tmp, "w");
  if (!f) return -1;
  fprintf(f, "{\"wcslib_version\": \"%s\", \"wcstrig_tol\": ", vers);
  jd(f, WCSTRIG_TOL);
  fprintf(f, ",\n\"rows\": [\n");
  int first = 1, nrow = 0;
  for (int base = -1; base <= 1; base += 2)
    for (int dir = -1; dir <= 1; dir += 2)
      for (size_t k = 0; k < sizeof(d) / sizeof(d[0]); k++) {
        if (dir < 0 && d[k] == 0.0) continue;
        /* outward (dir = +1) moves |v| above 1, inward below it */
        double v = base * (1.0 + dir * d[k]);
        fprintf(f, "%s[", first ? "" : ",\n");
        first = 0;
        jd(f, v);
        fprintf(f, ",");
        jd(f, asind(v));
        fprintf(f, ",");
        jd(f, acosd(v));
        fprintf(f, "]");
        nrow++;
      }
  fprintf(f, "\n]}\n");
  if (fclose(f) != 0 || rename(tmp, path) != 0) return -1;
  printf("  wcstrig_guards.json: %d arguments\n", nrow);
  return 0;
}

int main(int argc, char **argv) {
  if (argc >= 2 && strcmp(argv[1], "list") == 0) {
    for (int i = 0; i < ncases; i++) printf("%s\n", cases[i].name);
    return 0;
  }
  if (argc >= 3 && strcmp(argv[1], "trig") == 0) {
    int iv[3];
    return write_trig(argv[2], wcslib_version(iv)) != 0;
  }
  if (argc < 3 || strcmp(argv[1], "gen") != 0) {
    fprintf(stderr, "usage: %s gen <output_dir> | %s trig <output_dir> | %s list\n", argv[0], argv[0], argv[0]);
    return 2;
  }
  int iv[3];
  const char *vers = wcslib_version(iv);
  printf("wref (wcslib %s): %d cases -> %s\n", vers, ncases, argv[2]);
  for (int i = 0; i < ncases; i++)
    if (run_case(&cases[i], argv[2], vers) != 0) return 1;
  if (write_param_status(argv[2], vers) != 0) return 1;
  return 0;
}
