/*
 * wref.c — wcslib reference harness for validating the Go port of the
 * WCS projection code against the real wcslib implementation.
 *
 * Build:   make
 * Use:     ./wref gen golden/
 *
 * For each built-in test case this emits one JSON file in the output
 * directory, e.g. golden/szp_standard.json. The file contains:
 *
 *   code            3-letter projection code
 *   pv[]            projection parameters (PV2_0 .. PV2_20 or similar)
 *   forward:        array of {phi, theta, x, y, status}
 *   inverse:        array of {x, y, phi, theta, status}
 *
 * Angles in degrees (wcslib's native unit). The Go test loads these
 * files, converts inputs to radians, runs its pure-Go port, and diffs
 * the outputs at <1e-12 absolute tolerance (bit-exact within
 * double-precision rounding of the deg↔rad conversion).
 *
 * Covered projections:
 *   SZP, CSC, ZPN, QSC, MOL (iterative forward), PCO (regression),
 *   AZP (closed-form regression baseline)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <wcslib/prj.h>

/* ---------- JSON output helpers ---------- */

static void jprintDouble(FILE *f, double v) {
    /* Emit JSON null for NaN / inf so the Go parser accepts it; tests
     * use the "stat" field to decide whether a point is valid, not
     * the numeric output. */
    if (isnan(v) || isinf(v)) fprintf(f, "null");
    else                      fprintf(f, "%.17g", v);
}

/* ---------- one projection run ---------- */

typedef struct {
    const char *name;          /* output filename stem */
    const char *code;          /* 3-letter projection code, uppercase */
    int  npv;                  /* number of PV parameters */
    double pv[30];             /* pv[0..npv-1] copied into prj->pv[0..npv-1] */
    int  nphi;                 /* number of forward (phi, theta) test pairs */
    double phitheta[64][2];    /* (phi, theta) in degrees */
    int  nxy;                  /* number of inverse (x, y) test pairs */
    double xy[64][2];          /* (x, y) in degrees on the projection plane */
} case_t;

static int run_case(const case_t *c, const char *outdir) {
    struct prjprm prj;
    memset(&prj, 0, sizeof(prj));
    prjini(&prj);
    strncpy(prj.code, c->code, 3);
    prj.code[3] = '\0';
    for (int i = 0; i < c->npv; i++) prj.pv[i] = c->pv[i];

    int status = prjset(&prj);
    if (status != 0) {
        fprintf(stderr, "wref: %s: prjset status=%d (%s)\n", c->name, status, prj_errmsg[status]);
        prjfree(&prj);
        return -1;
    }

    /* Forward: (phi, theta) -> (x, y) */
    double fphi[64], ftheta[64], fx[64], fy[64];
    int fstat[64];
    for (int i = 0; i < c->nphi; i++) {
        fphi[i]   = c->phitheta[i][0];
        ftheta[i] = c->phitheta[i][1];
    }
    int fwdStatus = 0;
    if (c->nphi > 0) {
        fwdStatus = prj.prjs2x(&prj, c->nphi, 0, 1, 1, fphi, ftheta, fx, fy, fstat);
    }

    /* Inverse: (x, y) -> (phi, theta) */
    double ix[64], iy[64], iphi[64], itheta[64];
    int istat[64];
    for (int i = 0; i < c->nxy; i++) {
        ix[i] = c->xy[i][0];
        iy[i] = c->xy[i][1];
    }
    int invStatus = 0;
    if (c->nxy > 0) {
        invStatus = prj.prjx2s(&prj, c->nxy, 0, 1, 1, ix, iy, iphi, itheta, istat);
    }

    /* Write JSON. */
    char path[1024];
    snprintf(path, sizeof(path), "%s/%s.json", outdir, c->name);
    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "wref: open %s: %s\n", path, strerror(errno));
        prjfree(&prj);
        return -1;
    }
    fprintf(f, "{\n");
    fprintf(f, "  \"name\": \"%s\",\n", c->name);
    fprintf(f, "  \"code\": \"%s\",\n", c->code);
    fprintf(f, "  \"pv\": [");
    for (int i = 0; i < c->npv; i++) {
        if (i > 0) fprintf(f, ", ");
        jprintDouble(f, c->pv[i]);
    }
    fprintf(f, "],\n");
    fprintf(f, "  \"npv\": %d,\n", c->npv);
    fprintf(f, "  \"forward\": [\n");
    for (int i = 0; i < c->nphi; i++) {
        fprintf(f, "    {\"phi\": ");   jprintDouble(f, fphi[i]);
        fprintf(f, ", \"theta\": ");    jprintDouble(f, ftheta[i]);
        fprintf(f, ", \"x\": ");        jprintDouble(f, fx[i]);
        fprintf(f, ", \"y\": ");        jprintDouble(f, fy[i]);
        fprintf(f, ", \"stat\": %d}", fstat[i]);
        if (i+1 < c->nphi) fprintf(f, ",");
        fprintf(f, "\n");
    }
    fprintf(f, "  ],\n");
    fprintf(f, "  \"inverse\": [\n");
    for (int i = 0; i < c->nxy; i++) {
        fprintf(f, "    {\"x\": ");     jprintDouble(f, ix[i]);
        fprintf(f, ", \"y\": ");        jprintDouble(f, iy[i]);
        fprintf(f, ", \"phi\": ");      jprintDouble(f, iphi[i]);
        fprintf(f, ", \"theta\": ");    jprintDouble(f, itheta[i]);
        fprintf(f, ", \"stat\": %d}", istat[i]);
        if (i+1 < c->nxy) fprintf(f, ",");
        fprintf(f, "\n");
    }
    fprintf(f, "  ],\n");
    fprintf(f, "  \"forward_status\": %d,\n", fwdStatus);
    fprintf(f, "  \"inverse_status\": %d\n", invStatus);
    fprintf(f, "}\n");
    fclose(f);

    printf("  %s: %s npv=%d fwd=%d inv=%d\n", c->name, c->code, c->npv, c->nphi, c->nxy);

    prjfree(&prj);
    return 0;
}

/* ---------- test case definitions ----------
 *
 * For each projection we need a representative set of inputs. For the
 * forward grid we pick (phi, theta) pairs spanning the natural domain:
 * a few equatorial points, a few mid-latitude, some near-pole, and
 * both hemispheres. For the inverse grid we use the projected (x, y)
 * values that wcslib itself produces from the forward grid — the Go
 * test loads both halves and diffs against these canonical values.
 *
 * Since we don't want wref to do its own forward pass first, we just
 * hardcode plausible (x, y) inputs in the projection plane and let the
 * inverse do its thing. The Go port will be validated against these
 * same (x, y) inputs, so any disagreement is real.
 */

#define STDGRID_PT 25
static const double stdgrid[STDGRID_PT][2] = {
    /* phi (deg), theta (deg) — 5x5 grid skipping the exact pole */
    {  0,  0}, { 30,  0}, { 60,  0}, { 90,  0}, {120,  0},
    {  0, 30}, { 30, 30}, { 60, 30}, { 90, 30}, {120, 30},
    {  0, 60}, { 30, 60}, { 60, 60}, { 90, 60}, {120, 60},
    {  0, 85}, { 30, 85}, { 60, 85}, { 90, 85}, {120, 85},
    {  0,-30}, { 30,-30}, { 60,-60}, { 90,-60}, {120,-85},
};

/* Inverse grid: pick (x, y) on the projection plane that all the
 * projections can reasonably evaluate. For zenithal projections
 * these are radial distances in degrees from the origin. Quadcubes
 * need points within the 6-face cross. These values are well-behaved
 * for SZP, AZP, ZPN, PCO; CSC/QSC need a separate in-face grid. */
#define ZENITHAL_INV_PT 13
static const double zenithal_inv[ZENITHAL_INV_PT][2] = {
    {  0,   0},
    { 10,   0}, {-10,   0}, {  0,  10}, {  0, -10},
    { 20,  20}, {-20,  20}, { 20, -20}, {-20, -20},
    { 30,  45}, {-30,  45}, { 30, -45}, {-30, -45},
};

/* Cube-face grid: each face is π/2 × π/2 (90° × 90°). Face layout in the
 * cross pattern (wcslib puts face 0 at origin in native coords). */
#define CUBE_INV_PT 14
static const double cube_inv[CUBE_INV_PT][2] = {
    /* face 0 (equator, front) */
    {  0,  0}, { 20, 15}, {-20, -15}, { 30,  0}, {  0,  25},
    /* faces 1, 2, 3 (equator, right/back/left) */
    { 90,  0}, {135, 10}, {180,  0}, {-90,  0}, {-135, -10},
    /* faces 4, 5 (top, bottom) */
    {  0,  90}, { 15,  110}, {  0, -90}, {-15, -110},
};

/* MOL inverse grid — equal-area sphere maps onto an ellipse of semi-axes
 * (2*sqrt(2), sqrt(2)) × r0 in the default r0=180/π case, i.e. about
 * (162°, 81°). Stay inside. */
#define MOL_INV_PT 12
static const double mol_inv[MOL_INV_PT][2] = {
    {   0,   0}, {  60,   0}, { -60,   0}, { 120,   0}, {-120,   0},
    {   0,  40}, {  60,  40}, {-100,  20}, {  80, -30}, { -80, -40},
    {   0,  70}, {   0, -70},
};

static case_t cases[] = {
    /* ----- SZP: slant zenithal perspective ----- */
    {
        "szp_standard", "SZP", 3, {2.0, 180.0, 60.0}, /* mu=2, phi_c=180, theta_c=60 */
        STDGRID_PT, {{0}}, /* filled below */
        ZENITHAL_INV_PT, {{0}},
    },
    {
        "szp_infinity", "SZP", 3, {0.0, 0.0, 90.0}, /* mu=0 → reduces to TAN */
        STDGRID_PT, {{0}},
        ZENITHAL_INV_PT, {{0}},
    },
    {
        "szp_tilted", "SZP", 3, {1.5, 90.0, 45.0},
        STDGRID_PT, {{0}},
        ZENITHAL_INV_PT, {{0}},
    },

    /* ----- AZP regression ----- */
    {
        "azp_standard", "AZP", 2, {2.0, 0.0},
        STDGRID_PT, {{0}},
        ZENITHAL_INV_PT, {{0}},
    },
    {
        "azp_tilted", "AZP", 2, {2.0, 30.0},
        STDGRID_PT, {{0}},
        ZENITHAL_INV_PT, {{0}},
    },

    /* ----- ZPN: polynomial zenithal ----- */
    {
        "zpn_linear", "ZPN", 2, {0.0, 1.0}, /* R = zeta */
        STDGRID_PT, {{0}},
        ZENITHAL_INV_PT, {{0}},
    },
    {
        "zpn_quadratic", "ZPN", 3, {0.0, 1.0, -0.05},
        STDGRID_PT, {{0}},
        ZENITHAL_INV_PT, {{0}},
    },
    {
        "zpn_cubic", "ZPN", 4, {0.0, 1.0, 0.0, -0.02},
        STDGRID_PT, {{0}},
        ZENITHAL_INV_PT, {{0}},
    },
    {
        "zpn_quintic", "ZPN", 6, {0.0, 0.9, 0.0, 0.1, 0.0, -0.02},
        STDGRID_PT, {{0}},
        ZENITHAL_INV_PT, {{0}},
    },

    /* ----- TSC: tangential spherical cube (regression) ----- */
    {
        "tsc_allfaces", "TSC", 0, {0.0},
        STDGRID_PT, {{0}},
        CUBE_INV_PT, {{0}},
    },

    /* ----- CSC: COBE quadrilateralized ----- */
    {
        "csc_allfaces", "CSC", 0, {0.0},
        STDGRID_PT, {{0}},
        CUBE_INV_PT, {{0}},
    },

    /* ----- QSC: quadrilateralized ----- */
    {
        "qsc_allfaces", "QSC", 0, {0.0},
        STDGRID_PT, {{0}},
        CUBE_INV_PT, {{0}},
    },

    /* ----- MOL regression ----- */
    {
        "mol_standard", "MOL", 0, {0.0},
        STDGRID_PT, {{0}},
        MOL_INV_PT, {{0}},
    },

    /* ----- PCO regression (after our earlier port) ----- */
    {
        "pco_standard", "PCO", 0, {0.0},
        STDGRID_PT, {{0}},
        MOL_INV_PT, {{0}}, /* PCO and MOL share a similar inverse domain shape */
    },
};
static const int ncases = (int)(sizeof(cases) / sizeof(cases[0]));

/* Copy stdgrid / zenithal_inv / cube_inv / mol_inv into the per-case
 * arrays since C doesn't let us designate a static initializer to
 * another static array. */
static void fill_cases(void) {
    for (int i = 0; i < ncases; i++) {
        case_t *c = &cases[i];
        /* Forward grid is always stdgrid. */
        for (int j = 0; j < c->nphi; j++) {
            c->phitheta[j][0] = stdgrid[j][0];
            c->phitheta[j][1] = stdgrid[j][1];
        }
        /* Inverse grid depends on the family. */
        const double (*src)[2] = NULL;
        if (strcmp(c->code, "CSC") == 0 || strcmp(c->code, "QSC") == 0) {
            src = cube_inv;
        } else if (strcmp(c->code, "MOL") == 0 || strcmp(c->code, "PCO") == 0) {
            src = mol_inv;
        } else {
            src = zenithal_inv;
        }
        for (int j = 0; j < c->nxy; j++) {
            c->xy[j][0] = src[j][0];
            c->xy[j][1] = src[j][1];
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s gen <output_dir> | %s list\n", argv[0], argv[0]);
        return 2;
    }
    fill_cases();
    if (strcmp(argv[1], "list") == 0) {
        for (int i = 0; i < ncases; i++) printf("%s\n", cases[i].name);
        return 0;
    }
    if (strcmp(argv[1], "gen") != 0 || argc < 3) {
        fprintf(stderr, "usage: %s gen <output_dir>\n", argv[0]);
        return 2;
    }
    const char *outdir = argv[2];
    printf("wref: %d cases -> %s\n", ncases, outdir);
    for (int i = 0; i < ncases; i++) {
        if (run_case(&cases[i], outdir) != 0) return 1;
    }
    return 0;
}
