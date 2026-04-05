/*
 * cref.c — C reference harness for validating the Go port of cfitsio's
 * float quantization pipeline (fits_img_stats_float + fits_quantize_float).
 *
 * Build:   make
 * Use:     ./cref gen <output_dir>
 *
 * For each built-in test case this emits three files in <output_dir>:
 *
 *   <case>.in.bin        float32 input array, little-endian
 *   <case>.out.bin       int32  quantized output, little-endian
 *   <case>.out.json      scalar stats (noise1..noise5, mean, sigma, min,
 *                        max, ngood, bscale, bzero, iminval, imaxval,
 *                        quantize_return)
 *
 * The Go side reads the same .in.bin, runs the pure-Go port, and diffs
 * the result against .out.bin + .out.json. Any drift from cfitsio's
 * exact bytes is a bug in the port.
 *
 * Only the public fits_img_stats_float is declared in fitsio.h;
 * fits_quantize_float is internal but exported — we forward-declare it
 * here so we can call the real cfitsio implementation as the oracle.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <errno.h>
#include <fitsio.h>

/* Forward-declaration for the internal function. Signature mirrors
 * reference/cfitsio/quantize.c exactly. */
int fits_quantize_float(long row, float fdata[], long nxpix, long nypix,
    int nullcheck, float in_null_value, float qlevel, int dither_method,
    int idata[], double *bscale, double *bzero, int *iminval, int *imaxval);

/* Simple deterministic LCG for generating repeatable test data without
 * pulling in stdlib rand() (whose sequence differs across platforms). */
typedef struct {
    uint64_t state;
} lcg_t;

static void lcg_init(lcg_t *r, uint64_t seed) { r->state = seed * 6364136223846793005ULL + 1442695040888963407ULL; }

static uint32_t lcg_u32(lcg_t *r) {
    r->state = r->state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (uint32_t)(r->state >> 33);
}

/* Uniform [0,1) as double. */
static double lcg_uniform(lcg_t *r) {
    return (double)lcg_u32(r) / 4294967296.0;
}

/* Box-Muller, returns one N(0,1) sample per call (throws away the pair
 * partner — simpler + deterministic). */
static double lcg_gauss(lcg_t *r) {
    double u1 = lcg_uniform(r);
    double u2 = lcg_uniform(r);
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
}

/* ---------- test case definitions ---------- */

typedef struct {
    const char *name;
    long nx, ny;
    int nullcheck;
    float nullvalue;
    float qlevel;
    int dither_method;   /* 0=no dither, 1=SUB_DITHER_1, 2=SUB_DITHER_2 */
    long row;            /* if >0 and dither_method>0, dithered encode */
    /* generator */
    void (*fill)(float *out, long n);
} case_t;

static void fill_gauss_100_5(float *out, long n) {
    lcg_t r; lcg_init(&r, 42);
    for (long i = 0; i < n; i++) out[i] = (float)(100.0 + 5.0 * lcg_gauss(&r));
}

/* NULLF is the float sentinel used to mark null pixels for test cases
 * that exercise nullcheck=1. cfitsio's float tile writer also works by
 * substituting NaN with a sentinel before calling fits_quantize_float,
 * so this mirrors the real pipeline. NaN-as-sentinel doesn't work
 * because C's NaN != NaN breaks the equality check inside quantize.c. */
#define NULLF (-1.0e30f)

static void fill_gauss_with_nulls(float *out, long n) {
    lcg_t r; lcg_init(&r, 43);
    for (long i = 0; i < n; i++) {
        double v = 100.0 + 5.0 * lcg_gauss(&r);
        /* Every 37th pixel becomes a null (sentinel). */
        if (i % 37 == 0) out[i] = NULLF;
        else out[i] = (float)v;
    }
}

/* Generator for the full NaN-handling pipeline: emits NaN into the
 * output array, and the run_case driver notices nullcheck=2 to do
 * the NaN->sentinel substitution before calling fits_quantize_float.
 * This matches the cfitsio writer's two-phase NaN handling. */
static void fill_gauss_with_real_nans(float *out, long n) {
    lcg_t r; lcg_init(&r, 46);
    for (long i = 0; i < n; i++) {
        double v = 100.0 + 5.0 * lcg_gauss(&r);
        if (i % 41 == 0) out[i] = NAN;
        else out[i] = (float)v;
    }
}

static void fill_gauss_with_zeros(float *out, long n) {
    lcg_t r; lcg_init(&r, 44);
    for (long i = 0; i < n; i++) {
        double v = 100.0 + 5.0 * lcg_gauss(&r);
        if (i % 31 == 0) out[i] = 0.0f;
        else out[i] = (float)v;
    }
}

static void fill_gradient_plus_noise(float *out, long n) {
    lcg_t r; lcg_init(&r, 45);
    long nx = 64;
    for (long i = 0; i < n; i++) {
        double x = (double)(i % nx);
        double y = (double)(i / nx);
        double gradient = 50.0 + 0.5 * x + 0.25 * y;
        out[i] = (float)(gradient + 2.0 * lcg_gauss(&r));
    }
}

static void fill_constant(float *out, long n) {
    for (long i = 0; i < n; i++) out[i] = 42.5f;
}

static void fill_tiny(float *out, long n) {
    /* Exercises the nx<9 fallback path in FnNoise5. */
    const float vals[6] = {1.0f, 2.0f, 1.5f, 2.5f, 1.25f, 2.75f};
    for (long i = 0; i < n; i++) out[i] = vals[i % 6];
}

static case_t cases[] = {
    {"gauss_64x64_q4_nodither",        64, 64, 0, 0.0f, 4.0f,  0, 0,  fill_gauss_100_5},
    {"gauss_64x64_q4_dither1",         64, 64, 0, 0.0f, 4.0f,  1, 1,  fill_gauss_100_5},
    {"gauss_64x64_q4_dither1_row7",    64, 64, 0, 0.0f, 4.0f,  1, 7,  fill_gauss_100_5},
    {"gauss_64x64_q4_dither2_zeros",   64, 64, 0, 0.0f, 4.0f,  2, 1,  fill_gauss_with_zeros},
    {"gauss_64x64_q4_nulls",           64, 64, 1, NULLF, 4.0f, 1, 1,  fill_gauss_with_nulls},
    {"gauss_64x64_q4_real_nans",       64, 64, 2, NULLF, 4.0f, 1, 1,  fill_gauss_with_real_nans},
    {"gauss_64x64_q16_dither1",        64, 64, 0, 0.0f, 16.0f, 1, 1,  fill_gauss_100_5},
    {"gauss_64x64_abs_step",           64, 64, 0, 0.0f, -0.5f, 0, 0,  fill_gauss_100_5},
    {"gradient_64x64_q4_dither1",      64, 64, 0, 0.0f, 4.0f,  1, 1,  fill_gradient_plus_noise},
    {"constant_16x16",                 16, 16, 0, 0.0f, 4.0f,  0, 0,  fill_constant},
    {"tiny_6x1",                        6,  1, 0, 0.0f, 4.0f,  0, 0,  fill_tiny},
};
static const int ncases = (int)(sizeof(cases) / sizeof(cases[0]));

/* ---------- output helpers ---------- */

static int write_file(const char *path, const void *data, size_t n) {
    FILE *f = fopen(path, "wb");
    if (!f) { fprintf(stderr, "cref: open %s: %s\n", path, strerror(errno)); return -1; }
    size_t w = fwrite(data, 1, n, f);
    fclose(f);
    if (w != n) { fprintf(stderr, "cref: short write %s\n", path); return -1; }
    return 0;
}

/* Print a double so that round-trip parsing is lossless. printf %.17g
 * guarantees this for IEEE doubles. */
static void fprint_double(FILE *f, double v) {
    if (isnan(v))      fprintf(f, "\"nan\"");
    else if (isinf(v)) fprintf(f, v > 0 ? "\"inf\"" : "\"-inf\"");
    else               fprintf(f, "%.17g", v);
}

static int run_case(const case_t *c, const char *outdir) {
    long n = c->nx * c->ny;
    float *fdata = (float *)malloc((size_t)n * sizeof(float));
    float *fcopy = (float *)malloc((size_t)n * sizeof(float));
    int *idata   = (int   *)malloc((size_t)n * sizeof(int));
    if (!fdata || !fcopy || !idata) { fprintf(stderr, "cref: OOM\n"); return -1; }

    c->fill(fdata, n);
    /* Full-pipeline mode (nullcheck==2): substitute NaN with sentinel
     * before handing off to quantize. This mirrors what the cfitsio
     * writer does for float tiles with real NaN inputs. The input
     * .in.bin is written BEFORE substitution so the Go test sees the
     * same raw bytes the caller would. */
    int quant_nullcheck = c->nullcheck;
    if (c->nullcheck == 2) {
        for (long i = 0; i < n; i++) {
            if (isnan(fdata[i])) fdata[i] = c->nullvalue;  /* pre-subst for quantize */
        }
        quant_nullcheck = 1;
    }
    memcpy(fcopy, fdata, (size_t)n * sizeof(float));  /* quantize may mutate (it doesn't, but defensive) */

    /* --- fits_img_stats_float --- */
    long ngood = 0;
    float minv = 0, maxv = 0;
    double mean = 0, sigma = 0, noise1 = 0, noise2 = 0, noise3 = 0, noise5 = 0;
    int status = 0;
    fits_img_stats_float(fdata, c->nx, c->ny, quant_nullcheck, c->nullvalue,
        &ngood, &minv, &maxv, &mean, &sigma, &noise1, &noise2, &noise3, &noise5, &status);
    if (status != 0) { fprintf(stderr, "cref: %s: fits_img_stats_float status=%d\n", c->name, status); return -1; }

    /* --- fits_quantize_float --- */
    double bscale = 0, bzero = 0;
    int iminval = 0, imaxval = 0;
    int qret = fits_quantize_float(c->row, fcopy, c->nx, c->ny, quant_nullcheck,
        c->nullvalue, c->qlevel, c->dither_method, idata, &bscale, &bzero, &iminval, &imaxval);

    /* --- write outputs --- */
    char path[1024];
    snprintf(path, sizeof(path), "%s/%s.in.bin", outdir, c->name);
    if (write_file(path, fdata, (size_t)n * sizeof(float)) != 0) return -1;
    snprintf(path, sizeof(path), "%s/%s.out.bin", outdir, c->name);
    if (qret == 1) {
        if (write_file(path, idata, (size_t)n * sizeof(int32_t)) != 0) return -1;
    } else {
        /* Write empty file as marker that quantize returned 0 (not quantized). */
        if (write_file(path, NULL, 0) != 0) return -1;
    }

    snprintf(path, sizeof(path), "%s/%s.out.json", outdir, c->name);
    FILE *jf = fopen(path, "w");
    if (!jf) { fprintf(stderr, "cref: open %s: %s\n", path, strerror(errno)); return -1; }
    fprintf(jf, "{\n");
    fprintf(jf, "  \"name\": \"%s\",\n", c->name);
    fprintf(jf, "  \"nx\": %ld,\n", c->nx);
    fprintf(jf, "  \"ny\": %ld,\n", c->ny);
    fprintf(jf, "  \"nullcheck\": %d,\n", c->nullcheck);
    fprintf(jf, "  \"nullvalue\": ");  fprint_double(jf, (double)c->nullvalue); fprintf(jf, ",\n");
    fprintf(jf, "  \"qlevel\": ");  fprint_double(jf, (double)c->qlevel); fprintf(jf, ",\n");
    fprintf(jf, "  \"dither_method\": %d,\n", c->dither_method);
    fprintf(jf, "  \"row\": %ld,\n", c->row);
    fprintf(jf, "  \"stats\": {\n");
    fprintf(jf, "    \"ngood\": %ld,\n", ngood);
    fprintf(jf, "    \"min\": ");    fprint_double(jf, (double)minv);   fprintf(jf, ",\n");
    fprintf(jf, "    \"max\": ");    fprint_double(jf, (double)maxv);   fprintf(jf, ",\n");
    fprintf(jf, "    \"mean\": ");   fprint_double(jf, mean);           fprintf(jf, ",\n");
    fprintf(jf, "    \"sigma\": ");  fprint_double(jf, sigma);          fprintf(jf, ",\n");
    fprintf(jf, "    \"noise1\": "); fprint_double(jf, noise1);         fprintf(jf, ",\n");
    fprintf(jf, "    \"noise2\": "); fprint_double(jf, noise2);         fprintf(jf, ",\n");
    fprintf(jf, "    \"noise3\": "); fprint_double(jf, noise3);         fprintf(jf, ",\n");
    fprintf(jf, "    \"noise5\": "); fprint_double(jf, noise5);         fprintf(jf, "\n");
    fprintf(jf, "  },\n");
    fprintf(jf, "  \"quantize\": {\n");
    fprintf(jf, "    \"return\": %d,\n", qret);
    fprintf(jf, "    \"bscale\": ");  fprint_double(jf, bscale);  fprintf(jf, ",\n");
    fprintf(jf, "    \"bzero\": ");   fprint_double(jf, bzero);   fprintf(jf, ",\n");
    fprintf(jf, "    \"iminval\": %d,\n", iminval);
    fprintf(jf, "    \"imaxval\": %d\n",  imaxval);
    fprintf(jf, "  }\n");
    fprintf(jf, "}\n");
    fclose(jf);

    printf("  %s: ngood=%ld min=%g max=%g noise3=%g qret=%d bscale=%g bzero=%g\n",
        c->name, ngood, (double)minv, (double)maxv, noise3, qret, bscale, bzero);

    free(fdata); free(fcopy); free(idata);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s gen <output_dir> | %s list\n", argv[0], argv[0]);
        return 2;
    }
    if (strcmp(argv[1], "list") == 0) {
        for (int i = 0; i < ncases; i++) printf("%s\n", cases[i].name);
        return 0;
    }
    if (strcmp(argv[1], "gen") != 0 || argc < 3) {
        fprintf(stderr, "usage: %s gen <output_dir>\n", argv[0]);
        return 2;
    }
    const char *outdir = argv[2];
    printf("cref: %d test cases -> %s\n", ncases, outdir);
    for (int i = 0; i < ncases; i++) {
        if (run_case(&cases[i], outdir) != 0) return 1;
    }
    return 0;
}
