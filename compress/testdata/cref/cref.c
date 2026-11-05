/*
 * cref.c — C reference harness linked against the system libcfitsio. The
 * Go compress package and the root package's compressed-image writer are
 * diffed against what the real cfitsio produces.
 *
 * Build:   make
 * Modes:   ./cref gen <output_dir>     quantizer / noise-statistics goldens
 *                                      (described below) + the dither table
 *          ./cref encode <output_dir>  writer-parity goldens (RICE, HCOMPRESS,
 *                                      PLIO codec streams and whole images
 *                                      compressed by fits_img_compress; see
 *                                      the "encode" section further down)
 *          ./cref decode <file>        print cfitsio's decode of a file's
 *                                      first image as JSON (used live by the
 *                                      Go interop tests)
 *          ./cref list                 names of the gen cases
 *
 * gen: for each built-in test case this emits <output_dir>/<case>.out.json:
 * the generator ("fill") and the sha256 of the float32 little-endian
 * input it produced (raw NaNs included: hashed before the NaN -> sentinel
 * substitution of the nullcheck=2 case), the scalar stats (noise1..noise5,
 * mean, sigma, min, max, ngood) and quantize results (return, bscale,
 * bzero, iminval, imaxval), and "out": the int32 quantized output as an
 * exact digest (see jvalues; null when fits_quantize_float returned 0).
 *
 * Inputs are not committed: the Go side regenerates them with
 * internal/crefgen (a port of the generators below), checks the sha256,
 * runs the pure-Go port and diffs the result against the JSON. Any drift
 * from cfitsio's exact values is a bug in the port.
 *
 * gen also writes golden/fits_rand_value.json: libcfitsio's dither table
 * (fits_rand_value[0..N_RANDOM-1] after fits_init_randoms) as a digest.
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
#include <float.h>
#include <errno.h>
#include <fitsio.h>
#include <fitsio2.h>  /* fits_rcomp*, pl_p2li (exported by libcfitsio) */

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
    /* generator (fill_name: its name, recorded in the JSON so the Go
     * port in internal/crefgen regenerates the same input) */
    void (*fill)(float *out, long n);
    const char *fill_name;
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
/* a denormal null value (imgw/rice_f32_nullval_denorm_t40x5) */
#define DENORM_NULL 1.04488e-40f

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

/* Runs of exactly 3 equal pixels (and some of 5) inside noise. A 3-run
 * makes v5 == v6 == v7 without v3..v7 all equal, so FnNoise5 adds a 3rd/
 * 5th-order difference but no 2nd-order one (nvals2 < nvals) and then
 * takes quick_select_float(differences2, nvals) -- over nvals entries,
 * including ones left over (in quick_select's permuted order) from the
 * previous row. Pins that quirk in noise2 and in the quantizer. */
static void fill_gauss_with_runs3(float *out, long n) {
    lcg_t r; lcg_init(&r, 47);
    for (long i = 0; i < n; i++) {
        double v = 100.0 + 5.0 * lcg_gauss(&r);
        long k = i % 23;
        if (k >= 4 && k <= 6) out[i] = 97.25f;
        else if (k >= 14 && k <= 18) out[i] = 103.5f;
        else out[i] = (float)v;
    }
}
static void fill_runs3_with_nulls(float *out, long n) {
    fill_gauss_with_runs3(out, n);
    for (long i = 0; i < n; i++) if (i % 29 == 11) out[i] = NULLF;
}

/* Two 200-pixel rows: row 0 plain noise (sigma 5), row 1 noisier
 * (sigma 10) with a 3-run every 8 pixels. Row 1's 2nd-order median then
 * includes the entries row 0's quick_select_float left at the tail of
 * differences2, so the result depends on the exact permutation
 * quick_select leaves (a sort would leave the largest values there). */
static void fill_runs3_stale(float *out, long n) {
    lcg_t r; lcg_init(&r, 48);
    for (long i = 0; i < n; i++) {
        long row = i / 200, col = i % 200;
        double v = 100.0 + (row ? 10.0 : 5.0) * lcg_gauss(&r);
        if (row == 1 && col % 8 >= 2 && col % 8 <= 4) v = 100.0 + (double)(col / 8);
        out[i] = (float)v;
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

/* Quantizer range limits (fits_quantize_float, quantize.c). With
 * float minval/maxval and double delta:
 *   (maxval - minval) / delta > 2*2147483647 - N_RESERVED_VALUES -> not quantized
 *   (maxval - minval) / delta < 2147483647 - N_RESERVED_VALUES   -> zeropt = minval
 *   otherwise                                  -> zeropt = (minval + maxval) / 2.
 * maxval - minval and minval + maxval are FLOAT operations in C (both
 * operands float), rounded to float before the division. With an
 * absolute qlevel (delta = -qlevel, a float) the ratio of two floats can
 * never fall inside either 10-value margin, so the margins are reached
 * with a positive qlevel (delta = noise / qlevel, a double). */
static void fill_ramp(float *out, long n, double lo, double hi) {
    for (long i = 0; i < n; i++) out[i] = (float)(lo + (double)i * (hi - lo) / (double)(n - 1));
    out[0] = (float)lo;
    out[n - 1] = (float)hi;
}
/* exact span 2^32 - 13 steps of 2^-32 (quantizable), but maxval - minval
 * rounds to 1.0f in float: 2^32 steps, not quantized */
#define QLEVEL_2M32 (-2.3283064365386963e-10f) /* -2^-32 */
static void fill_range_fsub(float *out, long n) {
    fill_ramp(out, n, -(ldexp(1, -24) - 13 * ldexp(1, -32)), 1 - ldexp(1, -24));
}
/* 3e9 steps of 1: the centered branch (zeropt = (min+max)/2) */
static void fill_wide_center(float *out, long n) { fill_ramp(out, n, 0, 3.0e9); }
/* centered, and minval + maxval = 2147483400 rounds to 2147483392 in
 * float: bzero 1073741696, not 1073741700 */
static void fill_center_fsum(float *out, long n) { fill_ramp(out, n, -120, 2147483520.0); }
/* A 0,0,1,1 pattern (noise2 = 1.0483579 is the noise estimate) with one
 * extreme pixel M. With qlevel Q: (M - 0) / (1.0483579 / Q) =
 *   MARGIN: 4294967292.80, inside (2^32 - 12, 2^32 - 2]: not quantized
 *   HALF:   2147483638.39, inside [2^31 - 10, 2^31 - 1): centered
 * (found by searching Q and M with fits_img_stats_float). */
#define Q_MARGIN 0x1.ffffaap-1f
#define M_MARGIN 0x1.0c615cp+32f
#define Q_HALF   0x1.fffffap-1f
#define M_HALF   0x1.0c6132p+31f
static void fill_pattern_extreme(float *out, long n, float m) {
    for (long i = 0; i < n; i++) out[i] = (float)((i % 4) >= 2);
    out[33] = m;
}
static void fill_margin_pos(float *out, long n) { fill_pattern_extreme(out, n, M_MARGIN); }
/* (M - 0) / (1.0483579 / Q) = 4294967286.75: also inside the margin, but
 * the centered extremes quantize to about -+2147483643, clear of the
 * NULL_VALUE / ZERO_VALUE sentinels (at 4294967292.8 the minimum lands
 * on them, so a writer that refuses sentinel values would hide a missing
 * margin) */
#define Q_MARGIN_FAR 0x1.0000b8p+0f
#define M_MARGIN_FAR 0x1.0c606ep+32f
static void fill_margin_far(float *out, long n) { fill_pattern_extreme(out, n, M_MARGIN_FAR); }
static void fill_halfmargin_pos(float *out, long n) { fill_pattern_extreme(out, n, M_HALF); }
/* 4 x 5: FnNoise3's nx < 5 flatten; the extremes are not in the first row */
static void fill_tiny_4x5(float *out, long n) {
    for (long i = 0; i < n; i++) out[i] = 10.0f + (float)(i % 3);
    out[13] = 50.0f;
    out[17] = -7.0f;
}

/* Round 2: literal inputs that pin single branches of the noise
 * estimators and the quantizer (each one changes a statistic when the
 * branch is mutated). The same literals are in internal/crefgen. */
static void fill_lit(float *out, long n, const float *lit, long nlit) {
    for (long i = 0; i < n; i++) out[i] = lit[i % nlit];
}
/* FnNoise5, 9 x 3: row 3 has nvals == 1 and nvals2 == 0 (v5 == v6 == v7
 * at its one step, v3..v7 not all equal): noise2 takes no entry from it */
static const float lit_n5_one[27] = {3, 1, 4, 1, 5, 9, 2, 6, 5,
                                     2, 7, 1, 8, 2, 8, 1, 8, 2,
                                     1, 7, 2, 3, 5, 5, 5, 8, 4};
static void fill_n5_one(float *o, long n) { fill_lit(o, n, lit_n5_one, 27); }
/* FnNoise5, 10 x 3: row 3 has nvals == 2 and nvals2 == 1 */
static const float lit_n5_two[30] = {3, 1, 4, 1, 5, 9, 2, 6, 5, 3,
                                     2, 7, 1, 8, 2, 8, 1, 8, 2, 8,
                                     3, 8, 1, 2, 5, 5, 5, 9, 4, 6};
static void fill_n5_two(float *o, long n) { fill_lit(o, n, lit_n5_two, 30); }
/* FnNoise5 with nulls (-1), 9 x 3: row 3 has exactly 8 valid pixels,
 * so nvals == 0 and the row adds nothing to noise2/3/5 */
static const float lit_n5_zero[27] = {3, 1, 4, 1, 5, 9, 2, 6, 5,
                                      2, 7, 1, 8, 2, 8, 1, 8, 2,
                                      1, 7, 2, -1, 3, 5, 6, 8, 4};
static void fill_n5_zero(float *o, long n) { fill_lit(o, n, lit_n5_zero, 27); }
/* FnNoise1 with nulls (-1), 5 x 3: row 3 has 2 valid pixels (one
 * difference: nvals < 2, skipped) */
static const float lit_n1_short[15] = {1, 4, 2, 8, 3,
                                       5, 1, 9, 2, 6,
                                       -1, 7, -1, -1, 2};
static void fill_n1_short(float *o, long n) { fill_lit(o, n, lit_n1_short, 15); }
/* FnMeanSigma with one valid pixel (nulls -1), 3 x 3: sigma = 0 */
static const float lit_one_good[9] = {-1, -1, -1, -1, 7, -1, -1, -1, -1};
static void fill_one_good(float *o, long n) { fill_lit(o, n, lit_one_good, 9); }
/* FnNoise1's 5-sigma clip on an exact tie: 26 first differences, 25
 * zeros and one 26 (mean 1, stdev 5, |26 - 1| == 5 * 5 exactly). cfitsio
 * keeps a difference only when |d - mean| < 5 stdev, so the 26 is
 * clipped and the row's stdev is 0. */
static void fill_clip_tie(float *o, long n) {
    for (long i = 0; i < n; i++) o[i] = (i % 27) < 13 ? 26.0f : 0.0f;
}
/* Absolute qlevel with nulls (-1e30), 10 x 4: the minimum (3) is the
 * first valid pixel of row 2, read by FnNoise3's first-4-pixel loop */
static const float lit_abs_nulls[40] = {12, 15, 11, 18, 14, 16, 13, 17, 19, 12,
                                        NULLF, 3, 14, 16, 12, 15, 18, 11, 13, 17,
                                        16, 12, 15, 11, NULLF, 14, 18, 13, 12, 16,
                                        14, 17, 12, 15, 13, 11, 16, 18, 12, 14};
static void fill_abs_nulls(float *o, long n) { fill_lit(o, n, lit_abs_nulls, 40); }
/* every pixel null (-1e30): fits_quantize_float's ngood == 0 branch
 * (dummy minval 0, maxval 1, stdev 1) */
static void fill_all_null(float *o, long n) { for (long i = 0; i < n; i++) o[i] = NULLF; }
/* every pixel NaN (the writer's substitution path, nullcheck 2) */
static void fill_all_nan(float *o, long n) { for (long i = 0; i < n; i++) o[i] = NAN; }

#define F(fn) fn, #fn
static case_t cases[] = {
    {"gauss_64x64_q4_nodither",        64, 64, 0, 0.0f, 4.0f,  0, 0, F(fill_gauss_100_5)},
    {"gauss_64x64_q4_dither1",         64, 64, 0, 0.0f, 4.0f,  1, 1, F(fill_gauss_100_5)},
    {"gauss_64x64_q4_dither1_row7",    64, 64, 0, 0.0f, 4.0f,  1, 7, F(fill_gauss_100_5)},
    {"gauss_64x64_q4_dither2_zeros",   64, 64, 0, 0.0f, 4.0f,  2, 1, F(fill_gauss_with_zeros)},
    {"gauss_64x64_q4_nulls",           64, 64, 1, NULLF, 4.0f, 1, 1, F(fill_gauss_with_nulls)},
    {"gauss_64x64_q4_real_nans",       64, 64, 2, NULLF, 4.0f, 1, 1, F(fill_gauss_with_real_nans)},
    {"gauss_64x64_q16_dither1",        64, 64, 0, 0.0f, 16.0f, 1, 1, F(fill_gauss_100_5)},
    {"gauss_64x64_abs_step",           64, 64, 0, 0.0f, -0.5f, 0, 0, F(fill_gauss_100_5)},
    {"gradient_64x64_q4_dither1",      64, 64, 0, 0.0f, 4.0f,  1, 1, F(fill_gradient_plus_noise)},
    {"constant_16x16",                 16, 16, 0, 0.0f, 4.0f,  0, 0, F(fill_constant)},
    {"tiny_6x1",                        6,  1, 0, 0.0f, 4.0f,  0, 0, F(fill_tiny)},
    {"runs3_64x64_q4_dither1",         64, 64, 0, 0.0f, 4.0f,  1, 3, F(fill_gauss_with_runs3)},
    {"runs3_8x40_q4_nulls_dither2",     8, 40, 1, NULLF, 4.0f, 2, 5, F(fill_runs3_with_nulls)},
    {"runs3_stale_200x2_q4_dither1",  200,  2, 0, 0.0f, 4.0f,  1, 1, F(fill_runs3_stale)},
    /* appended: the quantizer's range limits and the branches the writer
     * never reaches (qlevel 0, FnNoise1's nx < 3, FnNoise3's nx < 5) */
    {"gauss_64x64_q0_dither1",        64, 64, 0, 0.0f, 0.0f,  1, 1,  F(fill_gauss_100_5)},
    {"tiny_2x8",                        2,  8, 0, 0.0f, 4.0f,  0, 0,  F(fill_gauss_100_5)},
    {"tiny_4x5_abs",                    4,  5, 0, 0.0f, -0.5f, 0, 0,  F(fill_tiny_4x5)},
    {"tiny_2x2_abs",                    2,  2, 0, 0.0f, -0.5f, 0, 0,  F(fill_tiny)},
    {"range_fsub_256x1",              256,  1, 0, 0.0f, QLEVEL_2M32, 1, 1, F(fill_range_fsub)},
    {"wide_center_64x1",               64,  1, 0, 0.0f, -1.0f, 1, 1,  F(fill_wide_center)},
    {"center_fsum_64x1",               64,  1, 0, 0.0f, -1.0f, 1, 1,  F(fill_center_fsum)},
    {"margin_64x1",                    64,  1, 0, 0.0f, Q_MARGIN, 1, 1, F(fill_margin_pos)},
    {"halfmargin_64x1",                64,  1, 0, 0.0f, Q_HALF,   1, 1, F(fill_halfmargin_pos)},
    /* round 2 (appended): single branches of the noise estimators and
     * of the quantizer's null handling */
    {"n5_one_9x3",                      9,  3, 0, 0.0f, 4.0f,  0, 0,  F(fill_n5_one)},
    {"n5_two_10x3",                    10,  3, 0, 0.0f, 4.0f,  0, 0,  F(fill_n5_two)},
    {"n5_zero_nulls_9x3",               9,  3, 1, -1.0f, 4.0f, 0, 0,  F(fill_n5_zero)},
    {"n1_short_nulls_5x3",              5,  3, 1, -1.0f, 4.0f, 0, 0,  F(fill_n1_short)},
    {"one_good_3x3",                    3,  3, 1, -1.0f, 4.0f, 0, 0,  F(fill_one_good)},
    {"clip_tie_27x1",                  27,  1, 0, 0.0f, 4.0f,  0, 0,  F(fill_clip_tie)},
    {"abs_nulls_10x4",                 10,  4, 1, NULLF, -0.5f, 0, 0, F(fill_abs_nulls)},
    {"abs_nulls_10x4_dither1",         10,  4, 1, NULLF, -0.5f, 1, 3, F(fill_abs_nulls)},
    {"all_null_4x4_q4_dither1",         4,  4, 1, NULLF, 4.0f, 1, 1,  F(fill_all_null)},
    {"all_null_12x3_q4_nodither",      12,  3, 1, NULLF, 4.0f, 0, 0,  F(fill_all_null)},
    {"all_null_4x4_abs_dither2",        4,  4, 1, NULLF, -0.5f, 2, 2, F(fill_all_null)},
    {"all_nan_10x10_q16_dither2",      10, 10, 2, NULLF, 16.0f, 2, 7, F(fill_all_nan)},
    /* round 3 (appended) */
    {"margin_far_64x1",                64,  1, 0, 0.0f, Q_MARGIN_FAR, 0, 0, F(fill_margin_far)},
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

/* ---- SHA-256 (FIPS 180-4), self-contained so the harness links only
 * against libcfitsio (the same code as testdata/cfref/cfref.c; checked
 * against Go's crypto/sha256 by every digest comparison). */
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

/* sha256 (hex) of a raw byte buffer: the generated input of a case. The
 * inputs are not committed; the Go tests regenerate them with a port of
 * these generators (internal/crefgen) and must reproduce this hash
 * before using them. */
static void sha_bytes_hex(const void *p, size_t n, char out[65]) {
  sha256_ctx c;
  sha_init(&c);
  sha_update(&c, p, n);
  sha_hex(&c, out);
}

static void jdbl_f(FILE *o, double d) {
  if (isnan(d)) fputs("\"NaN\"", o);
  else if (isinf(d)) fputs(d > 0 ? "\"Inf\"" : "\"-Inf\"", o);
  else fprintf(o, "%.17g", d);
}

/* ---- Value arrays: explicit when short, else an exact digest.
 *
 * The format of testdata/cfref/cfref.c (jdigest), decoded by
 * internal/testutil.Values: an array of more than DIGEST_MIN elements is
 *   {"n": N, "sha256": "<hex>", "nulls": [[start,count],...],
 *    "sample": [[index,value],...]}
 * sha256 over the canonical big-endian encoding of all N elements:
 *   f64  flag byte (1 valid, 0 null) + 8-byte IEEE bits (null -> 0)
 *   i64  flag byte + 8-byte two's complement
 * "nulls" lists every null as runs; "sample" holds explicit values at
 * the first and last 8 indices, the corners of the first and last 2-D
 * plane, the first minimum and maximum and 64 strided indices
 * k*(N-1)/63. A hash mismatch is a failure; the samples say where. */
#define DIGEST_MIN 64
enum { K_F64, K_I64 };

static int cmp_long(const void *a, const void *b) {
  long x = *(const long *)a, y = *(const long *)b;
  return x < y ? -1 : x > y;
}

static void jelem_f(FILE *o, int kind, const void *v, const char *nul, long k) {
  if (nul && nul[k]) { fputs("null", o); return; }
  if (kind == K_F64) jdbl_f(o, ((const double *)v)[k]);
  else fprintf(o, "%lld", ((const long long *)v)[k]);
}

/* naxes may be NULL (no plane corners). */
static void jvalues(FILE *o, int kind, const void *v, const char *nul, long n, int naxis, const long *naxes) {
  if (n <= DIGEST_MIN) {
    fputc('[', o);
    for (long k = 0; k < n; k++) { if (k) fputc(',', o); jelem_f(o, kind, v, nul, k); }
    fputc(']', o);
    return;
  }
  sha256_ctx c;
  sha_init(&c);
  for (long k = 0; k < n; k++) {
    int isnull = nul && nul[k];
    unsigned char flag = isnull ? 0 : 1;
    unsigned long long u = 0;
    if (!isnull) {
      if (kind == K_F64) {
        double d = ((const double *)v)[k];
        if (isnan(d)) u = 0x7ff8000000000000ULL; /* one canonical NaN */
        else memcpy(&u, &d, 8);
      } else u = (unsigned long long)((const long long *)v)[k];
    }
    sha_update(&c, &flag, 1);
    sha_u64(&c, u);
  }
  char hex[65];
  sha_hex(&c, hex);
  fprintf(o, "{\"n\": %ld, \"sha256\": \"%s\", \"nulls\": [", n, hex);
  int first = 1;
  if (nul) {
    for (long k = 0; k < n;) {
      if (!nul[k]) { k++; continue; }
      long s = k;
      while (k < n && nul[k]) k++;
      fprintf(o, "%s[%ld,%ld]", first ? "" : ",", s, k - s);
      first = 0;
    }
  }
  fprintf(o, "], \"sample\": [");
  long idx[8 + 8 + 8 + 2 + 64];
  int ni = 0;
  for (long k = 0; k < 8 && k < n; k++) idx[ni++] = k;
  for (long k = n - 8; k < n; k++) if (k >= 0) idx[ni++] = k;
  if (naxes && naxis >= 2) {
    long nx = naxes[0], ny = naxes[1], plane = nx * ny;
    long bases[2] = {0, (n / plane - 1) * plane};
    for (int b = 0; b < 2; b++) {
      idx[ni++] = bases[b]; idx[ni++] = bases[b] + nx - 1;
      idx[ni++] = bases[b] + (ny - 1) * nx; idx[ni++] = bases[b] + plane - 1;
    }
  }
  long kmin = -1, kmax = -1;
  for (long k = 0; k < n; k++) {
    if (nul && nul[k]) continue;
    if (kind == K_F64) {
      const double *x = v;
      if (isnan(x[k])) continue;
      if (kmin < 0 || x[k] < x[kmin]) kmin = k;
      if (kmax < 0 || x[k] > x[kmax]) kmax = k;
    } else {
      const long long *x = v;
      if (kmin < 0 || x[k] < x[kmin]) kmin = k;
      if (kmax < 0 || x[k] > x[kmax]) kmax = k;
    }
  }
  if (kmin >= 0) { idx[ni++] = kmin; idx[ni++] = kmax; }
  for (long k = 0; k < 64; k++) idx[ni++] = k * (n - 1) / 63;
  qsort(idx, (size_t)ni, sizeof idx[0], cmp_long);
  first = 1;
  for (int i = 0; i < ni; i++) {
    if (i && idx[i] == idx[i - 1]) continue;
    fprintf(o, "%s[%ld,", first ? "" : ",", idx[i]);
    jelem_f(o, kind, v, nul, idx[i]);
    fputc(']', o);
    first = 0;
  }
  fprintf(o, "]}");
}

static int run_case(const case_t *c, const char *outdir) {
    long n = c->nx * c->ny;
    float *fdata = (float *)malloc((size_t)n * sizeof(float));
    float *fcopy = (float *)malloc((size_t)n * sizeof(float));
    int *idata   = (int   *)malloc((size_t)n * sizeof(int));
    if (!fdata || !fcopy || !idata) { fprintf(stderr, "cref: OOM\n"); return -1; }

    c->fill(fdata, n);
    /* The input hash is taken BEFORE any NaN substitution: the Go test
     * regenerates the same raw values (including the NaNs). */
    char in_sha[65];
    sha_bytes_hex(fdata, (size_t)n * sizeof(float), in_sha);
    char path[1024];

    /* Full-pipeline mode (nullcheck==2): substitute NaN with sentinel
     * before handing off to quantize. This mirrors what the cfitsio
     * writer does for float tiles with real NaN inputs. */
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

    /* --- write outputs: one JSON file; the quantized integers as an
     * exact digest ("out": null when fits_quantize_float returned 0) --- */
    snprintf(path, sizeof(path), "%s/%s.out.json", outdir, c->name);
    FILE *jf = fopen(path, "w");
    if (!jf) { fprintf(stderr, "cref: open %s: %s\n", path, strerror(errno)); return -1; }
    fprintf(jf, "{\"name\": \"%s\", \"nx\": %ld, \"ny\": %ld, \"nullcheck\": %d, \"nullvalue\": ",
        c->name, c->nx, c->ny, c->nullcheck);
    fprint_double(jf, (double)c->nullvalue);
    fprintf(jf, ", \"qlevel\": ");  fprint_double(jf, (double)c->qlevel);
    fprintf(jf, ", \"dither_method\": %d, \"row\": %ld, \"fill\": \"%s\", \"in_sha256\": \"%s\",\n",
        c->dither_method, c->row, c->fill_name, in_sha);
    fprintf(jf, " \"stats\": {\"ngood\": %ld, \"min\": ", ngood);
    fprint_double(jf, (double)minv);
    fprintf(jf, ", \"max\": ");    fprint_double(jf, (double)maxv);
    fprintf(jf, ", \"mean\": ");   fprint_double(jf, mean);
    fprintf(jf, ", \"sigma\": ");  fprint_double(jf, sigma);
    fprintf(jf, ", \"noise1\": "); fprint_double(jf, noise1);
    fprintf(jf, ", \"noise2\": "); fprint_double(jf, noise2);
    fprintf(jf, ", \"noise3\": "); fprint_double(jf, noise3);
    fprintf(jf, ", \"noise5\": "); fprint_double(jf, noise5);
    fprintf(jf, "},\n \"quantize\": {\"return\": %d, \"bscale\": ", qret);
    fprint_double(jf, bscale);
    fprintf(jf, ", \"bzero\": ");   fprint_double(jf, bzero);
    fprintf(jf, ", \"iminval\": %d, \"imaxval\": %d},\n \"out\": ", iminval, imaxval);
    if (qret == 1) {
        long long *iv = malloc((size_t)n * sizeof(long long));
        for (long i = 0; i < n; i++) iv[i] = idata[i];
        long naxes[2] = {c->nx, c->ny};
        jvalues(jf, K_I64, iv, NULL, n, 2, naxes);
        free(iv);
    } else {
        fputs("null", jf);
    }
    fprintf(jf, "}\n");
    fclose(jf);

    printf("  %s: ngood=%ld min=%g max=%g noise3=%g qret=%d bscale=%g bzero=%g\n",
        c->name, ngood, (double)minv, (double)maxv, noise3, qret, bscale, bzero);

    free(fdata); free(fcopy); free(idata);
    return 0;
}

/* ======================================================================
 * Writer-parity goldens ("encode" mode).
 *
 * These cases feed adversarial inputs to the real libcfitsio ENCODERS and
 * record their exact output, so the Go writer can be diffed byte-for-byte:
 *
 * Every <case>.json records the generator ("kind", "gen_seed") and
 * "in_sha256", the sha256 of the input bytes below; the inputs are not
 * committed (internal/crefgen regenerates them and checks the hash).
 *
 *   golden/rice/<case>.rice     fits_rcomp / fits_rcomp_short /
 *                               fits_rcomp_byte output bytes
 *   golden/rice/<case>.json     bytepix, n, blocksize, nbytes; input hashed
 *                               little-endian, BYTEPIX wide (BITPIX 8 data
 *                               is the raw unsigned byte)
 *
 *   golden/hcomp/<case>.hc      fits_hcompress (BITPIX 8/16) or
 *                               fits_hcompress64 (BITPIX 32) output bytes,
 *                               called exactly as imcompress.c does:
 *                               fits_hcompress(a, NAXIS1, NAXIS2, scale, ..)
 *   golden/hcomp/<case>.json    bitpix, naxis1, naxis2, scale, nbytes; input
 *                               hashed as int32 little-endian
 *
 *   golden/plio/<case>.plio     pl_p2li line list, big-endian int16 words
 *                               (the bytes of the 1PI COMPRESSED_DATA cell)
 *   golden/plio/<case>.json     n, nwords; input hashed as int32 LE
 *
 *   golden/img/<case>.fits      a whole image compressed by cfitsio's own
 *                               writer (fits_img_compress, the fpack path)
 *   golden/img/<case>.json      every compression parameter; input hashed
 *                               little-endian in the image's native type
 *                               (u8 s8 i16 u16 i32 u32 f32); for lossy
 *                               cases (float32, HCOMPRESS scale != 0)
 *                               "dec": cfitsio's decode of <case>.fits
 *                               (fits_read_img TDOUBLE, NaN -> null flag,
 *                               float32 images rounded to float32) as an
 *                               exact digest. For lossless integer cases
 *                               the harness instead checks that cfitsio's
 *                               decode equals the input and records
 *                               "lossless_verified". "status" is the
 *                               cfitsio status of the compression (non-zero
 *                               means cfitsio REFUSED the input and there is
 *                               no .fits: the Go writer must refuse it too)
 *
 * Note: cfitsio sizes its HCOMPRESS output buffer as 2.2 (8/16-bit) or
 * 4.4 (32-bit) bytes per pixel + 26 (imcomp_calc_max_elem) and does not
 * check it; full-range random tiles exceed that and corrupt the heap
 * (observed with libcfitsio 4.7.0). The HCOMPRESS image cases therefore
 * use noise-like content; full-range data is covered by golden/hcomp,
 * which calls fits_hcompress* directly with a large enough buffer.
 *
 *   golden/MANIFEST             "<kind> <case>" for every case above, so the
 *                               Go tests can fail when a golden goes missing.
 * ====================================================================== */

#define BP_MIN(bp) ((bp) == 1 ? 0LL : (bp) == 2 ? -32768LL : (long long)INT32_MIN)
#define BP_MAX(bp) ((bp) == 1 ? 255LL : (bp) == 2 ? 32767LL : (long long)INT32_MAX)

static long long clampll(long long v, long long lo, long long hi) {
    return v < lo ? lo : v > hi ? hi : v;
}

/* Uniform integer in [lo, hi] (inclusive) from two LCG draws. */
static long long lcg_range(lcg_t *r, long long lo, long long hi) {
    unsigned long long span = (unsigned long long)(hi - lo) + 1ULL;
    unsigned long long u = ((unsigned long long)lcg_u32(r) << 32) | lcg_u32(r);
    return lo + (long long)(u % span);
}

/* ---- integer pattern generators (shared by rice / hcomp / img) ---- */

/* The four extreme bit patterns of the range, cycled: 0, max, min, -1 (for
 * BITPIX 8: 0x00 0x7f 0x80 0xff). Every successive difference overflows
 * the signed pixel width, which is the wrap case in ricecomp.c. */
static void gen_extremes(long long *o, long n, int bp) {
    long long pat1[4] = {0, 127, 128, 255};
    long long patn[4] = {0, BP_MAX(bp), BP_MIN(bp), -1};
    for (long i = 0; i < n; i++) o[i] = (bp == 1) ? pat1[i % 4] : patn[i % 4];
}
static void gen_jumps(long long *o, long n, int bp, uint64_t seed) {
    lcg_t r; lcg_init(&r, seed);
    for (long i = 0; i < n; i++) o[i] = lcg_range(&r, BP_MIN(bp), BP_MAX(bp));
}
static void gen_smooth(long long *o, long n, int bp, uint64_t seed) {
    lcg_t r; lcg_init(&r, seed);
    long long mid = (bp == 1) ? 100 : (bp == 2) ? 1000 : 100000;
    for (long i = 0; i < n; i++)
        o[i] = clampll(mid + (i % 97) + (long long)llround(3.0 * lcg_gauss(&r)), BP_MIN(bp), BP_MAX(bp));
}
/* Per 32-pixel block: constant, low-noise, or full-range, so one stream
 * exercises the fs==0, normal and fs>=fsmax (literal) block encodings. */
static void gen_mixed(long long *o, long n, int bp, uint64_t seed) {
    lcg_t r; lcg_init(&r, seed);
    for (long i = 0; i < n; i++) {
        int kind = (int)((i / 32) % 3);
        if (kind == 0) o[i] = (bp == 1) ? 7 : -5;
        else if (kind == 1) o[i] = clampll(((bp == 1) ? 50 : -20) + (long long)llround(2.0 * lcg_gauss(&r)), BP_MIN(bp), BP_MAX(bp));
        else o[i] = lcg_range(&r, BP_MIN(bp), BP_MAX(bp));
    }
}

typedef struct {
    const char *name;
    int bytepix;
    long n;
    int blocksize;
    int kind; /* 0 extremes 1 jumps 2 smooth 3 mixed 4 zeros */
} rice_case_t;

static const rice_case_t rice_cases[] = {
#define RICE_FOR(bp, tag) \
    {"rice_" tag "_extremes", bp, 1000, 32, 0}, \
    {"rice_" tag "_jumps",    bp, 1000, 32, 1}, \
    {"rice_" tag "_smooth",   bp, 1000, 32, 2}, \
    {"rice_" tag "_mixed",    bp, 1024, 32, 3}, \
    {"rice_" tag "_zeros",    bp,  100, 32, 4}, \
    {"rice_" tag "_tail",     bp, 1001, 32, 3}, \
    {"rice_" tag "_bs16",     bp, 1000, 16, 3}, \
    {"rice_" tag "_single",   bp,    1, 32, 1}
    RICE_FOR(1, "b"),
    RICE_FOR(2, "s"),
    RICE_FOR(4, "i"),
#undef RICE_FOR
};
static const int n_rice_cases = (int)(sizeof(rice_cases) / sizeof(rice_cases[0]));

static void fill_int_kind(long long *o, long n, int bp, int kind, uint64_t seed) {
    switch (kind) {
    case 0: gen_extremes(o, n, bp); break;
    case 1: gen_jumps(o, n, bp, seed); break;
    case 2: gen_smooth(o, n, bp, seed); break;
    case 3: gen_mixed(o, n, bp, seed); break;
    default: for (long i = 0; i < n; i++) o[i] = 0; break;
    }
}

static FILE *manifest;
static const char *only; /* optional single-case filter (debugging) */

static int put_le(unsigned char *dst, long long v, int width) {
    for (int b = 0; b < width; b++) dst[b] = (unsigned char)(((unsigned long long)v >> (8 * b)) & 0xff);
    return width;
}

static int run_rice_case(const rice_case_t *c, const char *dir) {
    long n = c->n;
    long long *v = calloc((size_t)n, sizeof(long long));
    unsigned char *raw = calloc((size_t)n, 8);
    size_t clen = (size_t)n * 8 + 1024;
    unsigned char *out = calloc(clen, 1);
    if (!v || !raw || !out) { fprintf(stderr, "cref: OOM\n"); return -1; }
    uint64_t seed = 1000 + (uint64_t)c->bytepix * 10 + (uint64_t)c->kind;
    fill_int_kind(v, n, c->bytepix, c->kind, seed);
    int nbytes;
    if (c->bytepix == 1) {
        signed char *a = malloc((size_t)n);
        for (long i = 0; i < n; i++) { a[i] = (signed char)(unsigned char)v[i]; raw[i] = (unsigned char)v[i]; }
        nbytes = fits_rcomp_byte(a, (int)n, out, (int)clen, c->blocksize);
        free(a);
    } else if (c->bytepix == 2) {
        short *a = malloc((size_t)n * sizeof(short));
        for (long i = 0; i < n; i++) { a[i] = (short)v[i]; put_le(raw + 2 * i, v[i], 2); }
        nbytes = fits_rcomp_short(a, (int)n, out, (int)clen, c->blocksize);
        free(a);
    } else {
        int *a = malloc((size_t)n * sizeof(int));
        for (long i = 0; i < n; i++) { a[i] = (int)v[i]; put_le(raw + 4 * i, v[i], 4); }
        nbytes = fits_rcomp(a, (int)n, out, (int)clen, c->blocksize);
        free(a);
    }
    if (nbytes < 0) { fprintf(stderr, "cref: %s: rice compress failed %d\n", c->name, nbytes); return -1; }
    char path[1024], in_sha[65];
    sha_bytes_hex(raw, (size_t)n * (size_t)c->bytepix, in_sha);
    snprintf(path, sizeof path, "%s/%s.rice", dir, c->name);
    if (write_file(path, out, (size_t)nbytes) != 0) return -1;
    snprintf(path, sizeof path, "%s/%s.json", dir, c->name);
    FILE *jf = fopen(path, "w");
    if (!jf) return -1;
    fprintf(jf, "{\"name\": \"%s\", \"bytepix\": %d, \"n\": %ld, \"blocksize\": %d, \"nbytes\": %d, \"kind\": %d, \"gen_seed\": %llu, \"in_sha256\": \"%s\"}\n",
        c->name, c->bytepix, n, c->blocksize, nbytes, c->kind, (unsigned long long)seed, in_sha);
    fclose(jf);
    fprintf(manifest, "rice %s\n", c->name);
    free(v); free(raw); free(out);
    return 0;
}

/* ---- HCOMPRESS ---- */

typedef struct {
    const char *name;
    int bitpix;          /* 8, 16 or 32: picks fits_hcompress vs fits_hcompress64 */
    int naxis1, naxis2;  /* NAXIS1 is the fastest axis */
    int scale;
    int kind;            /* 0 extremes 1 jumps 2 smooth 3 mixed 5 noise-image */
} hc_case_t;

static const hc_case_t hc_cases[] = {
    {"hc_i16_smooth_40x25_s0",   16, 40, 25,   0, 2},
    {"hc_i16_noise_40x25_s0",    16, 40, 25,   0, 5},
    {"hc_i16_noise_40x25_s4",    16, 40, 25,   4, 5},
    {"hc_i16_jumps_16x16_s0",    16, 16, 16,   0, 1},
    {"hc_i16_extremes_12x8_s0",  16, 12,  8,   0, 0},
    {"hc_i16_odd_5x7_s0",        16,  5,  7,   0, 3},
    {"hc_i16_min_4x4_s0",        16,  4,  4,   0, 1},
    {"hc_u8_jumps_37x23_s0",      8, 37, 23,   0, 1},
    {"hc_u8_noise_37x23_s3",      8, 37, 23,   3, 5},
    {"hc_i32_noise_40x25_s0",    32, 40, 25,   0, 5},
    {"hc_i32_noise_40x25_s100",  32, 40, 25, 100, 5},
    {"hc_i32_jumps_33x17_s0",    32, 33, 17,   0, 1},
    {"hc_i32_extremes_16x12_s0", 32, 16, 12,   0, 0},
    /* round 2 (appended): quadrants with one side a power of two and the
     * other one larger (nqmax and log2n in qtree_encode, e.g. nqx 4,
     * nqy 5), in both orientations */
    {"hc_i16_jumps_8x10_s0",     16,  8, 10,   0, 1},
    {"hc_i16_jumps_10x8_s0",     16, 10,  8,   0, 1},
    {"hc_i16_smooth_16x18_s0",   16, 16, 18,   0, 2},
    {"hc_i16_smooth_18x16_s0",   16, 18, 16,   0, 2},
    {"hc_i16_jumps_8x9_s0",      16,  8,  9,   0, 1},
    {"hc_i16_jumps_9x8_s0",      16,  9,  8,   0, 1},
    {"hc_i32_jumps_16x18_s0",    32, 16, 18,   0, 1},
};
static const int n_hc_cases = (int)(sizeof(hc_cases) / sizeof(hc_cases[0]));

static void gen_noise_image(long long *o, long n, int bp, uint64_t seed) {
    lcg_t r; lcg_init(&r, seed);
    double mid = (bp == 1) ? 120 : (bp == 2) ? 1000 : 3.0e6;
    double sig = (bp == 1) ? 10 : (bp == 2) ? 50 : 4.0e4;
    for (long i = 0; i < n; i++)
        o[i] = clampll((long long)llround(mid + sig * lcg_gauss(&r) + (double)(i % 40)), BP_MIN(bp), BP_MAX(bp));
}

static int run_hc_case(const hc_case_t *c, const char *dir) {
    long n = (long)c->naxis1 * c->naxis2;
    int bp = c->bitpix / 8;
    long long *v = calloc((size_t)n, sizeof(long long));
    unsigned char *raw = calloc((size_t)n, 4);
    long nbytes = n * 16 + 1024;
    char *out = calloc((size_t)nbytes, 1);
    int status = 0;
    uint64_t seed = (c->kind == 5) ? 2000 + (uint64_t)c->bitpix : 3000 + (uint64_t)c->bitpix * 10 + (uint64_t)c->kind;
    if (c->kind == 5) gen_noise_image(v, n, bp, seed);
    else fill_int_kind(v, n, bp, c->kind, seed);
    for (long i = 0; i < n; i++) put_le(raw + 4 * i, v[i], 4);
    if (c->bitpix == 32) {
        LONGLONG *a = malloc((size_t)n * sizeof(LONGLONG));
        for (long i = 0; i < n; i++) a[i] = v[i];
        fits_hcompress64(a, c->naxis1, c->naxis2, c->scale, out, &nbytes, &status);
        free(a);
    } else {
        int *a = malloc((size_t)n * sizeof(int));
        for (long i = 0; i < n; i++) a[i] = (int)v[i];
        fits_hcompress(a, c->naxis1, c->naxis2, c->scale, out, &nbytes, &status);
        free(a);
    }
    if (status) { fprintf(stderr, "cref: %s: hcompress status %d\n", c->name, status); return -1; }
    char path[1024], in_sha[65];
    sha_bytes_hex(raw, (size_t)n * 4, in_sha);
    snprintf(path, sizeof path, "%s/%s.hc", dir, c->name);
    if (write_file(path, out, (size_t)nbytes) != 0) return -1;
    snprintf(path, sizeof path, "%s/%s.json", dir, c->name);
    FILE *jf = fopen(path, "w");
    if (!jf) return -1;
    fprintf(jf, "{\"name\": \"%s\", \"bitpix\": %d, \"naxis1\": %d, \"naxis2\": %d, \"scale\": %d, \"nbytes\": %ld, \"kind\": %d, \"gen_seed\": %llu, \"in_sha256\": \"%s\"}\n",
        c->name, c->bitpix, c->naxis1, c->naxis2, c->scale, nbytes, c->kind, (unsigned long long)seed, in_sha);
    fclose(jf);
    fprintf(manifest, "hcomp %s\n", c->name);
    free(v); free(raw); free(out);
    return 0;
}

/* ---- PLIO ---- */

typedef struct {
    const char *name;
    long n;
    int kind; /* 0 mask runs, 1 long runs, 2 big values, 3 isolated pixels, 4 zeros, 5 random 0..2^24-1, 6 +-4095/4096/4097 steps */
} plio_case_t;

static const plio_case_t plio_cases[] = {
    {"plio_mask",       1000, 0},
    {"plio_long_runs", 12000, 1},
    {"plio_big_values", 1000, 2},
    {"plio_isolated",   1000, 3},
    {"plio_zeros",       500, 4},
    {"plio_random24",    700, 5},
    {"plio_single",        1, 2},
    /* > 32767 words: exercises line-list header word 5 ((op-1)/32768) */
    {"plio_long_list", 12000, 5},
    /* runs whose value changes by exactly +-4095, +-4096 and +-4097: the
     * |dv| > 4095 boundary between the 1-word (12-bit) and 2-word
     * (24-bit) set-high-value instructions, for isolated pixels and runs */
    {"plio_dv_boundary", 600, 6},
};
static const int n_plio_cases = (int)(sizeof(plio_cases) / sizeof(plio_cases[0]));

static void gen_plio(long long *o, long n, int kind, uint64_t seed) {
    lcg_t r; lcg_init(&r, seed);
    for (long i = 0; i < n; i++) {
        switch (kind) {
        case 0: o[i] = ((i / 13) % 5 == 0) ? 0 : ((i / 13) % 5 == 1) ? 1 : ((i / 13) % 5 == 2) ? 2 : ((i / 13) % 5 == 3) ? 0 : 7; break;
        case 1: o[i] = (i < 5000) ? 0 : (i < 9500) ? 3 : (i < 9600) ? 4100 : 1; break;
        case 2: { static const long long big[6] = {16777215, 0, 4096, 1, 16777215, 8000000}; o[i] = big[(i / 9) % 6]; } break;
        case 3: o[i] = (i % 7 == 3) ? (long long)(1 + (i % 11)) : 0; break;
        case 4: o[i] = 0; break;
        case 6: {
            /* successive distinct values: 1 -> 4096 (+4095) -> 8192 (+4096)
             * -> 4096 (-4096) -> 8193 (+4097) -> 4096 (-4097) -> 8191
             * (+4095) -> 4096 (-4095) -> ... ; run length 1, 2 or 3 and a
             * zero gap every few runs */
            static const long long dv[14] = {1, 4096, 8192, 4096, 8193, 4096, 8191, 4096,
                                             1, 4097, 1, 4098, 4097, 1};
            long seg = i / 3, len = 1 + (seg % 3), pos = i % 3;
            o[i] = (seg % 5 == 4) ? 0 : (pos < len) ? dv[seg % 14] : dv[(seg + 1) % 14];
        } break;
        default: o[i] = lcg_range(&r, 0, 16777215); break;
        }
    }
}

static int run_plio_case(const plio_case_t *c, const char *dir) {
    long n = c->n;
    long long *v = calloc((size_t)n, sizeof(long long));
    int *a = calloc((size_t)n, sizeof(int));
    unsigned char *raw = calloc((size_t)n, 4);
    short *ll = calloc((size_t)n * 4 + 64, sizeof(short));
    uint64_t seed = 4000 + (uint64_t)c->kind;
    gen_plio(v, n, c->kind, seed);
    for (long i = 0; i < n; i++) { a[i] = (int)v[i]; put_le(raw + 4 * i, v[i], 4); }
    int nw = pl_p2li(a, 1, ll, (int)n);
    if (nw <= 0) { fprintf(stderr, "cref: %s: pl_p2li returned %d\n", c->name, nw); return -1; }
    unsigned char *be = calloc((size_t)nw, 2);
    for (int i = 0; i < nw; i++) { be[2 * i] = (unsigned char)(((unsigned short)ll[i]) >> 8); be[2 * i + 1] = (unsigned char)(((unsigned short)ll[i]) & 0xff); }
    char path[1024], in_sha[65];
    sha_bytes_hex(raw, (size_t)n * 4, in_sha);
    snprintf(path, sizeof path, "%s/%s.plio", dir, c->name);
    if (write_file(path, be, (size_t)nw * 2) != 0) return -1;
    snprintf(path, sizeof path, "%s/%s.json", dir, c->name);
    FILE *jf = fopen(path, "w");
    if (!jf) return -1;
    fprintf(jf, "{\"name\": \"%s\", \"n\": %ld, \"nwords\": %d, \"kind\": %d, \"gen_seed\": %llu, \"in_sha256\": \"%s\"}\n",
        c->name, n, nw, c->kind, (unsigned long long)seed, in_sha);
    fclose(jf);
    fprintf(manifest, "plio %s\n", c->name);
    free(v); free(a); free(raw); free(ll); free(be);
    return 0;
}

/* ---- whole images through cfitsio's writer (fits_img_compress) ---- */

typedef struct {
    const char *name;
    const char *type;   /* u8 s8 i16 u16 i32 u32 f32 */
    int naxis;
    long naxes[7];
    int algo;           /* RICE_1 GZIP_1 GZIP_2 HCOMPRESS_1 PLIO_1 */
    long tile[7];       /* all zero = cfitsio default tiling */
    float qlevel;       /* 0 = cfitsio default */
    int qmethod;        /* 0 = default, NO_DITHER (-1), 1, 2 */
    int seed;           /* ZDITHER0 (float only) */
    float hscale;       /* HCOMPRESS scale as passed to fits_set_hcomp_scale */
    int kind;           /* int: 0 adversarial 1 mask 2 noise | float: 0 noise 1 zeros 2 nan 3 const 4 huge 5 IEEE specials 6 huge + specials
                         * 7 NaN + real -1e30 8 NaN + real FLOATNULLVALUE + -1e30 9 range margin 10 centered range 11 half margin
                         * 12 NaN every 37th and real -1e30 every 41st pixel
                         * 13 per 200 pixels: a NaN (3) and a real FLOATNULLVALUE (10) */
    int hdr;            /* 1: the source image carries extra keywords, 2: a BLANK keyword (float image);
                         * the source is kept as <case>.src.fits.
                         * 3: no source image: the pixels are written straight into the compressed
                         * image with fits_write_imgnull and null value -1e30 (the generator's NaNs
                         * passed as -1e30), so cfitsio quantizes with that null value instead of
                         * FLOATNULLVALUE, as the Go writer does; every tile must hold a null.
                         * 4: the source image's SIMPLE card is patched to F (cfitsio copies it
                         * to ZSIMPLE = F) */
    double nullv;       /* hdr 3: the null value given to fits_write_imgnull (NaN pixels are
                         * passed as it); 0 = -1e30. Recorded as "nullvalue" in the JSON. */
} img_case_t;

#define D1(a) 1, {a, 0, 0}
#define D2(a, b) 2, {a, b, 0}
#define D3(a, b, c) 3, {a, b, c}
#define D4(a, b, c, d) 4, {a, b, c, d}
#define D6(a, b, c, d, e, f) 6, {a, b, c, d, e, f}
#define D7(a, b, c, d, e, f, g) 7, {a, b, c, d, e, f, g}
#define T0 {0, 0, 0}
#define T2(a, b) {a, b, 0}
#define T4(a, b, c, d) {a, b, c, d}

static const img_case_t img_cases[] = {
#define INT_SET(tag) \
    {"rice_" tag,            tag, D2(40, 25), RICE_1,      T0,        0, 0, 0, 0, 0, 0}, \
    {"rice_" tag "_t16x8",   tag, D2(40, 25), RICE_1,      T2(16, 8), 0, 0, 0, 0, 0, 0}, \
    {"gzip1_" tag,           tag, D2(40, 25), GZIP_1,      T0,        0, 0, 0, 0, 0, 0}, \
    {"gzip2_" tag "_t16x8",  tag, D2(40, 25), GZIP_2,      T2(16, 8), 0, 0, 0, 0, 0, 0}, \
    {"hcomp_" tag,           tag, D2(40, 25), HCOMPRESS_1, T0,        0, 0, 0, 0, 2, 0}, \
    {"plio_" tag,            tag, D2(40, 25), PLIO_1,      T0,        0, 0, 0, 0, 1, 0}
    INT_SET("u8"),
    INT_SET("s8"),
    INT_SET("i16"),
    INT_SET("u16"),
    INT_SET("i32"),
    INT_SET("u32"),
#undef INT_SET
    {"hcomp_i16_noise_s4",     "i16", D2(40, 25), HCOMPRESS_1, T0,         0, 0, 0, -4,   2, 0},
    {"hcomp_i32_noise_s25",    "i32", D2(40, 25), HCOMPRESS_1, T0,         0, 0, 0, -25,  2, 0},
    {"hcomp_u16_noise_s2",     "u16", D2(40, 25), HCOMPRESS_1, T0,         0, 0, 0, -2,   2, 0},
    {"hcomp_i16_40x45",        "i16", D2(40, 45), HCOMPRESS_1, T0,         0, 0, 0, 0,    2, 0},
    {"hcomp_i32_t20x10",       "i32", D2(40, 25), HCOMPRESS_1, T2(20, 10), 0, 0, 0, 0,    2, 0},
    {"hcomp_u8_adversarial",   "u8",  D2(40, 25), HCOMPRESS_1, T0,         0, 0, 0, 0,    0, 0},
    {"hcomp_s8_adversarial",   "s8",  D2(40, 25), HCOMPRESS_1, T0,         0, 0, 0, 0,    0, 0},
    {"rice_i16_cube_t5x4x3",   "i16", D3(12, 9, 4), RICE_1,    {5, 4, 3},  0, 0, 0, 0,    0, 0},
    {"gzip2_i32_cube_t12x9x1", "i32", D3(12, 9, 4), GZIP_2,    {12, 9, 1},        0, 0, 0, 0,    0, 0},
    {"plio_i16_negative",      "i16", D2(40, 25), PLIO_1,      T0,         0, 0, 0, 0,    0, 0},
    {"plio_i32_adversarial",   "i32", D2(40, 25), PLIO_1,      T0,         0, 0, 0, 0,    0, 0},
    /* float32: quantized through fits_quantize_float, then the integer codec */
    {"rice_f32_sd1_z1",        "f32", D2(40, 25), RICE_1, T0,        4,    1,  1,     0, 0, 0},
    {"rice_f32_sd1_z7_t16x8",  "f32", D2(40, 25), RICE_1, T2(16, 8), 4,    1,  7,     0, 0, 0},
    {"rice_f32_sd1_z10000",    "f32", D2(40, 25), RICE_1, T0,        4,    1,  10000, 0, 0, 0},
    {"rice_f32_sd2_z1_zeros",  "f32", D2(40, 25), RICE_1, T0,        4,    2,  1,     0, 1, 0},
    {"rice_f32_nodither_q16",  "f32", D2(40, 25), RICE_1, T0,        16,   -1, 0,     0, 0, 0},
    {"rice_f32_nodither_qdef", "f32", D2(40, 25), RICE_1, T0,        0,    -1, 0,     0, 0, 0},
    {"rice_f32_defaults_z1",   "f32", D2(40, 25), RICE_1, T0,        0,    0,  1,     0, 0, 0},
    {"gzip1_f32_sd1_z3",       "f32", D2(40, 25), GZIP_1, T0,        4,    1,  3,     0, 0, 0},
    {"gzip2_f32_sd2_z5_t16x8", "f32", D2(40, 25), GZIP_2, T2(16, 8), 4,    2,  5,     0, 1, 0},
    {"rice_f32_abs_z2",        "f32", D2(40, 25), RICE_1, T0,        -0.5f, 1, 2,     0, 0, 0},
    {"rice_f32_nan_z1_t16x8",  "f32", D2(40, 25), RICE_1, T2(16, 8), 4,    1,  1,     0, 2, 0},
    {"gzip2_f32_nan_nodither", "f32", D2(40, 25), GZIP_2, T0,        16,   -1, 0,     0, 2, 0},
    {"rice_f32_const_t40x5",   "f32", D2(40, 25), RICE_1, T2(40, 5), 4,    1,  1,     0, 3, 0},
    {"rice_f32_huge",          "f32", D2(40, 25), RICE_1, T2(40, 5), 4,    1,  1,     0, 4, 0},
    {"rice_f32_wrap_100x100",  "f32", D2(100, 100), RICE_1, T2(100, 100), 4, 1, 9000, 0, 0, 0},
    {"rice_f32_cube",          "f32", D3(12, 9, 4), RICE_1, {12, 9, 1}, 4, 1, 4,     0, 0, 0},
    /* Appended cases (appending keeps the seeds, 5000 + index, of the
     * cases above). HCOMPRESS_1 default tile height (imcomp_init_table):
     * whole image up to 30 rows, then 16 (31), 24 (33), 20 (49), 28 (241). */
    {"hcomp_i16_8x30",         "i16", D2(8, 30),  HCOMPRESS_1, T0, 0, 0, 0, 0, 2, 0},
    {"hcomp_i16_8x31",         "i16", D2(8, 31),  HCOMPRESS_1, T0, 0, 0, 0, 0, 2, 0},
    {"hcomp_i16_8x33",         "i16", D2(8, 33),  HCOMPRESS_1, T0, 0, 0, 0, 0, 2, 0},
    {"hcomp_i16_8x49",         "i16", D2(8, 49),  HCOMPRESS_1, T0, 0, 0, 0, 0, 2, 0},
    {"hcomp_i16_4x241",        "i16", D2(4, 241), HCOMPRESS_1, T0, 0, 0, 0, 0, 2, 0},
    /* 1-D images and 3-D HCOMPRESS/PLIO (default tiling) */
    {"rice_i16_1d",            "i16", D1(1000),     RICE_1,      T0, 0, 0, 0, 0, 0, 0},
    {"gzip1_f32_1d_sd1_z2",    "f32", D1(1000),     GZIP_1,      T0, 4, 1, 2, 0, 2, 0},
    {"hcomp_i16_cube",         "i16", D3(12, 9, 3), HCOMPRESS_1, T0, 0, 0, 0, 0, 2, 0},
    {"plio_i16_cube",          "i16", D3(12, 9, 3), PLIO_1,      T0, 0, 0, 0, 0, 1, 0},
    /* caller keywords: SIMPLE/EXTEND/BLOCKED/CHECKSUM/DATASUM are renamed
     * Z*, the rest copied (imcomp_copy_img2comp) */
    {"rice_i16_hdr",           "i16", D2(40, 25), RICE_1, T0, 0, 0, 0, 0, 0, 1},
    /* ZDITHER0 = 10000: the single 10000-pixel tile starts at random-table
     * row 9999 and walks off its end, wrapping iseed 9999 -> 0 inside the
     * tile (both the no-null and the null branch of fits_quantize_float) */
    {"rice_f32_wrapin_z10000",     "f32", D2(100, 100), RICE_1, T2(100, 100), 4, 1, 10000, 0, 0, 0},
    {"rice_f32_wrapin_nan_z10000", "f32", D2(100, 100), RICE_1, T2(100, 100), 4, 1, 10000, 0, 2, 0},
    /* float IEEE specials (kind 5: +-Inf, NaN, denormals, -0, +-FLT_MIN;
     * kind 6: the same over the rice_f32_huge content, so the first tile
     * cannot be quantized and goes to the GZIP_COMPRESSED_DATA fallback) */
    {"rice_f32_inf_sd1_z1",          "f32", D2(40, 25), RICE_1, T0,        4, 1,  1, 0, 5, 0},
    {"rice_f32_inf_sd2_z3_t16x8",    "f32", D2(40, 25), RICE_1, T2(16, 8), 4, 2,  3, 0, 5, 0},
    {"rice_f32_inf_nodither",        "f32", D2(40, 25), RICE_1, T0,        0, -1, 0, 0, 5, 0},
    {"gzip2_f32_inf_fallback_t40x5", "f32", D2(40, 25), GZIP_2, T2(40, 5), 4, 1,  1, 0, 6, 0},
    /* HCOMPRESS_1 tile limits (imcomp_init_table): a 3-pixel tile is
     * refused; an 8-pixel tile on 11 columns (3-pixel edge) is silently
     * enlarged to 11 by cfitsio (our writer refuses it); a 4-pixel tile
     * and 4-pixel edge tiles are accepted. (A 4x4 tile is not used: 16
     * pixels overflow cfitsio's undersized HCOMPRESS buffer, see above;
     * golden/hcomp/hc_i16_min_4x4_s0 covers it with a large buffer.) */
    {"hcomp_i16_t3x8",               "i16", D2(40, 25), HCOMPRESS_1, T2(3, 8), 0, 0, 0, 0, 2, 0},
    {"hcomp_i16_11x8_t8x8",          "i16", D2(11, 8),  HCOMPRESS_1, T2(8, 8), 0, 0, 0, 0, 2, 0},
    {"hcomp_i16_12x9_t8x5",          "i16", D2(12, 9),  HCOMPRESS_1, T2(8, 5), 0, 0, 0, 0, 2, 0},
    {"hcomp_i16_12x24_t4x12",        "i16", D2(12, 24), HCOMPRESS_1, T2(4, 12), 0, 0, 0, 0, 2, 0},
    /* a requested tile longer than the image axis is written as requested
     * (ZTILE1 = 100 on an 11-pixel axis), not clamped */
    {"rice_i16_11x8_t100x3",         "i16", D2(11, 8),  RICE_1, T2(100, 3), 0, 0, 0, 0, 0, 0},
    {"gzip1_f32_11x8_t100x3",        "f32", D2(11, 8),  GZIP_1, T2(100, 3), 4, 1, 1, 0, 0, 0},
    /* more than 3 dimensions: cfitsio compresses up to MAX_COMPRESS_DIM = 6
     * axes and refuses 7 (fits_img_compress "too many dimensions") */
    {"rice_i16_4d_t10x6x2x1",        "i16", D4(10, 6, 4, 2), RICE_1, T4(10, 6, 2, 1), 0, 0, 0, 0, 2, 0},
    {"gzip2_i16_6d",                 "i16", D6(5, 4, 3, 2, 2, 2), GZIP_2, T0, 0, 0, 0, 0, 2, 0},
    {"rice_i16_7d",                  "i16", D7(4, 3, 2, 2, 2, 2, 2), RICE_1, T0, 0, 0, 0, 0, 2, 0},
    /* Round 2 (appended). Null sentinel: tile 1 holds a NaN and a real
     * -1e30 (the writer's default NullValue, Nullcheck off): the real pixel
     * stays a value and the tile goes to GZIP_COMPRESSED_DATA. */
    {"rice_f32_nan_realnull_t40x5",  "f32", D2(40, 25), RICE_1, T2(40, 5), 4, 1, 1, 0, 7, 0},
    /* FLOATNULLVALUE (-9.11912e-36) as a real pixel in a tile with a NaN:
     * cfitsio's writer turns it into a null (imcomp_float2nan); the Go
     * writer keeps it (a documented difference). */
    {"rice_f32_floatnullvalue_t40x5", "f32", D2(40, 25), RICE_1, T2(40, 5), 4, 1, 1, 0, 8, 0},
    /* quantizer range limits through the whole writer: row 2 spans
     * 4294967292.8 steps (inside the N_RESERVED_VALUES margin:
     * GZIP_COMPRESSED_DATA), row 1 has a NaN (so ZBLANK is written); a
     * centered tile of 3e9 steps and one of 2147483638.4 steps */
    {"rice_f32_margin_64x2",         "f32", D2(64, 2), RICE_1, T0, Q_MARGIN, 1, 1, 0, 9, 0},
    {"rice_f32_center_64",           "f32", D1(64),    RICE_1, T0, -1, 1, 1, 0, 10, 0},
    {"rice_f32_halfmargin_64",       "f32", D1(64),    RICE_1, T0, Q_HALF, 1, 1, 0, 11, 0},
    /* dither-table wrap at ZDITHER0 9000 inside a tile that has nulls
     * (the null branch's re-seed lands on a non-zero offset) */
    {"rice_f32_wrap_nan_z9000",      "f32", D2(100, 100), RICE_1, T2(100, 100), 4, 1, 9000, 0, 2, 0},
    /* a float source image with BLANK = 5: cfitsio copies BLANK into the
     * compressed header and every reader then nulls the real pixels whose
     * quantized integer is 5; the Go writer refuses BLANK for float input */
    {"rice_f32_blank_hdr",           "f32", D2(40, 25), RICE_1, T0, 4, 1, 1, 0, 0, 2},
    /* whole-row tiles, most with a NaN and a real -1e30: many
     * GZIP_COMPRESSED_DATA tiles with nulls (cfitsio's fits_read_imgnull
     * leaves stale null flags across them; see "decodenull") */
    {"rice_f32_nan_realnull_rows",   "f32", D2(40, 25), RICE_1, T0, 4, 1, 1, 0, 12, 0},
};
static const int n_img_cases = (int)(sizeof(img_cases) / sizeof(img_cases[0]));

/* golden/imgw: images cfitsio compresses while they are WRITTEN
 * (fits_write_imgnull into a compressing HDU, hdr 3) instead of through
 * fits_img_compress. Same JSON as img/; input seed 7000 + index; MANIFEST
 * kind "imgw". A real FLOATNULLVALUE beside a NaN in QUANTIZED tiles (img/
 * only reaches it in a GZIP_COMPRESSED_DATA fallback tile): with null value
 * -1e30 cfitsio quantizes the real pixel as a value, as the Go writer
 * does. */
static const img_case_t imgw_cases[] = {
    {"rice_f32_fnv_quant_t40x5",     "f32", D2(40, 10), RICE_1, T2(40, 5), 4,  1, 1, 0, 13, 3},
    {"gzip2_f32_fnv_quant_nodither", "f32", D2(40, 10), GZIP_2, T2(40, 5), 16, -1, 0, 0, 13, 3},
    /* Round 3 (appended): a caller null value. fits_write_imgnull nulls
     * the pixels equal to it (compared raw, before any denormal is read as
     * 0): -999 (float32 and float64) and a denormal 1.04488e-40, beside
     * another denormal (-3e-42) that stays a value. The Go writer gets the
     * same pixels with Nullcheck and NullValue / NullValue64. */
    {"rice_f32_nullval_m999_t40x5",  "f32", D2(40, 10), RICE_1, T2(40, 5), 4,  1, 1, 0, 17, 3, -999.0},
    {"rice_f32_nullval_denorm_t40x5","f32", D2(40, 10), RICE_1, T2(40, 5), 4,  1, 1, 0, 18, 3, (double)DENORM_NULL},
    {"rice_f64_nullval_m999_t40x5",  "f64", D2(40, 10), RICE_1, T2(40, 5), 4,  1, 1, 0, 17, 3, -999.0},
};
static const int n_imgw_cases = (int)(sizeof(imgw_cases) / sizeof(imgw_cases[0]));

/* golden/img64: BITPIX -64 images through fits_img_compress (each tile
 * quantized by fits_quantize_double). Same JSON as img/ (input hashed as
 * float64 little-endian; "dec" is cfitsio's TDOUBLE decode, not rounded);
 * input seed 8000 + index; MANIFEST kind "img64". */
static const img_case_t img64_cases[] = {
    {"rice_f64_sd1_z1",              "f64", D2(40, 25),   RICE_1, T0,           4,     1,  1,    0, 0, 0},
    {"rice_f64_sd2_z3_zeros_t16x8",  "f64", D2(40, 25),   RICE_1, T2(16, 8),    4,     2,  3,    0, 1, 0},
    {"gzip1_f64_nodither_q16",       "f64", D2(40, 25),   GZIP_1, T0,           16,    -1, 0,    0, 0, 0},
    {"gzip2_f64_sd1_z5_nan_t16x8",   "f64", D2(40, 25),   GZIP_2, T2(16, 8),    4,     1,  5,    0, 2, 0},
    {"rice_f64_abs_z2",              "f64", D2(40, 25),   RICE_1, T0,           -0.5f, 1,  2,    0, 0, 0},
    {"rice_f64_const_t40x5",         "f64", D2(40, 25),   RICE_1, T2(40, 5),    4,     1,  1,    0, 3, 0},
    {"rice_f64_huge",                "f64", D2(40, 25),   RICE_1, T2(40, 5),    4,     1,  1,    0, 4, 0},
    {"rice_f64_inf_sd1_z1",          "f64", D2(40, 25),   RICE_1, T0,           4,     1,  1,    0, 5, 0},
    {"gzip2_f64_inf_fallback_t40x5", "f64", D2(40, 25),   GZIP_2, T2(40, 5),    4,     1,  1,    0, 6, 0},
    {"rice_f64_nan_realnull_t40x5",  "f64", D2(40, 25),   RICE_1, T2(40, 5),    4,     1,  1,    0, 7, 0},
    {"rice_f64_dnv_t40x5",           "f64", D2(40, 25),   RICE_1, T2(40, 5),    4,     1,  1,    0, 8, 0},
    {"rice_f64_precise_1e8",         "f64", D2(40, 25),   RICE_1, T0,           4,     1,  1,    0, 14, 0},
    {"rice_f64_cube",                "f64", D3(12, 9, 4), RICE_1, {12, 9, 1},   4,     1,  4,    0, 0, 0},
    {"rice_f64_wrap_nan_z9000",      "f64", D2(100, 100), RICE_1, T2(100, 100), 4,     1,  9000, 0, 2, 0},
    {"rice_f64_defaults",            "f64", D2(40, 25),   RICE_1, T0,           0,     0,  1,    0, 0, 0},
    {"rice_f64_nodither_qdef",       "f64", D2(40, 25),   RICE_1, T0,           0,     -1, 0,    0, 0, 0},
    {"rice_f64_beyond_flt_max",      "f64", D2(40, 25),   RICE_1, T0,           4,     1,  1,    0, 15, 0},
    {"rice_f64_beyond_flt_max_abs",  "f64", D2(40, 25),   RICE_1, T2(40, 5),    -1.0e32f, 1, 1,  0, 15, 0},
    /* Round 3 (appended). SUBTRACTIVE_DITHER_2 in tiles that hold both
     * NaNs and exact zeros (the zero rule of the quantizer's null branch). */
    {"rice_f64_sd2_nan_zeros_t16x8", "f64", D2(40, 25),   RICE_1, T2(16, 8),    4,     2,  3,    0, 16, 0},
    /* An absolute qlevel so small that zeropt / delta (1e22) overflows
     * fits_quantize_double's int64 iqfactor: cfitsio writes the tile with
     * status 0 and it decodes to about -9.2e6 instead of 1e10 (a cfitsio
     * defect; the Go writer stores the tile losslessly instead). */
    {"rice_f64_zpt_overflow_abs",    "f64", D1(64),       RICE_1, T0,           -1.0e-12f, -1, 0, 0, 19, 0},
};
static const int n_img64_cases = (int)(sizeof(img64_cases) / sizeof(img64_cases[0]));

/* golden/img2: more fits_img_compress cases (round 3), kept apart from
 * golden/img so the img set the reader tests pin stays as it is. Same JSON
 * as img/; input seed 9000 + index; MANIFEST kind "img2". */
static const img_case_t img2_cases[] = {
    /* SUBTRACTIVE_DITHER_2 in tiles with NaNs and exact zeros */
    {"rice_f32_sd2_nan_zeros_t16x8", "f32", D2(40, 25), RICE_1, T2(16, 8), 4, 2, 3, 0, 16, 0},
    /* a source image with SIMPLE = F (and the hdr 1 keywords): cfitsio
     * writes ZSIMPLE = F */
    {"rice_i16_simple_f",            "i16", D2(40, 25), RICE_1, T0,        0, 0, 0, 0, 0,  4},
    /* cfitsio defects, status 0 and wrong pixels: a constant tile of 100
     * with qlevel -1e-20 overflows the int64 iqfactor (zeropt / delta =
     * 1e22); a 4294966272-step ramp from 1e10 takes the centered branch
     * and its float sum (minval + maxval) rounds 1024 low, so the maximum
     * quantizes to 2^31 and wraps. The Go writer stores both losslessly. */
    {"rice_f32_zpt_overflow_abs",    "f32", D1(64),     RICE_1, T0,        -1.0e-20f, -1, 0, 0, 19, 0},
    {"rice_f32_center_overflow",     "f32", D1(64),     RICE_1, T0,        -1.0f,     -1, 0, 0, 20, 0},
    /* a cfitsio defect, status 0: a real pixel quantized onto a sentinel.
     * Row 1's centered minimum becomes NULL_VALUE (the image has ZBLANK
     * from row 2's NaN): every reader returns it as a null. The Go writer
     * stores row 1 losslessly. */
    {"rice_f32_null_sentinel_hit",   "f32", D2(64, 2),  RICE_1, T0,        -0x1.528c8cp+0f, -1, 0, 0, 21, 0},
};
static const int n_img2_cases = (int)(sizeof(img2_cases) / sizeof(img2_cases[0]));

static int type_bitpix(const char *t) {
    if (!strcmp(t, "f64")) return DOUBLE_IMG;
    if (!strcmp(t, "u8")) return BYTE_IMG;
    if (!strcmp(t, "s8")) return SBYTE_IMG;
    if (!strcmp(t, "i16")) return SHORT_IMG;
    if (!strcmp(t, "u16")) return USHORT_IMG;
    if (!strcmp(t, "i32")) return LONG_IMG;
    if (!strcmp(t, "u32")) return ULONG_IMG;
    return FLOAT_IMG;
}
static int type_datatype(const char *t) {
    if (!strcmp(t, "f64")) return TDOUBLE;
    if (!strcmp(t, "u8")) return TBYTE;
    if (!strcmp(t, "s8")) return TSBYTE;
    if (!strcmp(t, "i16")) return TSHORT;
    if (!strcmp(t, "u16")) return TUSHORT;
    if (!strcmp(t, "i32")) return TINT;
    if (!strcmp(t, "u32")) return TUINT;
    return TFLOAT;
}
static int type_width(const char *t) {
    if (!strcmp(t, "f64")) return 8;
    if (!strcmp(t, "u8") || !strcmp(t, "s8")) return 1;
    if (!strcmp(t, "i16") || !strcmp(t, "u16")) return 2;
    return 4;
}
static void type_range(const char *t, long long *lo, long long *hi) {
    if (!strcmp(t, "u8")) { *lo = 0; *hi = 255; }
    else if (!strcmp(t, "s8")) { *lo = -128; *hi = 127; }
    else if (!strcmp(t, "i16")) { *lo = -32768; *hi = 32767; }
    else if (!strcmp(t, "u16")) { *lo = 0; *hi = 65535; }
    else if (!strcmp(t, "i32")) { *lo = INT32_MIN; *hi = INT32_MAX; }
    else { *lo = 0; *hi = 4294967295LL; }
}

/* Integer image content. kind 0 (adversarial): 5 bands of rows -
 * smooth ramp, full-range random, the 4 range extremes cycled, low noise
 * around the middle, constant. kind 1 (mask): PLIO-style runs inside
 * [0, min(max, 2^24-1)], reaching above 32767 where the type allows it.
 * kind 2 (noise): smooth noisy sky for lossy HCOMPRESS. */
static void gen_img_int(long long *o, long n, long nx, const char *t, int kind, uint64_t seed) {
    lcg_t r; lcg_init(&r, seed);
    long long lo, hi; type_range(t, &lo, &hi);
    long long mid = lo + (hi - lo) / 2;
    for (long i = 0; i < n; i++) {
        long row = i / nx, col = i % nx;
        if (kind == 1) {
            long long top = hi < 16777215 ? hi : 16777215;
            long long vals[6] = {0, 1, 2, top, top / 2 + 1, 0};
            o[i] = ((row / 3) % 2 == 0 && col > 5 && col < 30) ? vals[(col / 4) % 6] : 0;
            continue;
        }
        if (kind == 2) {
            double sig = (hi - lo) > 70000 ? 4.0e4 : (hi - lo) > 300 ? 50 : 8;
            double base = (hi - lo) > 70000 ? (double)mid / 1000.0 : (double)mid;
            if (lo == 0 && (hi - lo) > 70000) base = 3.0e6;
            o[i] = clampll((long long)llround(base + sig * lcg_gauss(&r) + (double)col), lo, hi);
            continue;
        }
        switch ((row / 5) % 5) {
        case 0: o[i] = clampll(lo + (hi - lo) / 3 + col * 3 + row, lo, hi); break;
        case 1: o[i] = lcg_range(&r, lo, hi); break;
        case 2: { long long e[4] = {lo, hi, mid, lo + 1}; o[i] = e[i % 4]; } break;
        case 3: o[i] = clampll(mid + (long long)llround(4.0 * lcg_gauss(&r)), lo, hi); break;
        default: o[i] = mid; break;
        }
    }
}

static void gen_img_float(float *o, long n, long nx, int kind, uint64_t seed) {
    lcg_t r; lcg_init(&r, seed);
    if (kind == 9) {
        /* row 1: the 0,0,1,1 pattern and a NaN; row 2: margin_64x1 */
        fill_pattern_extreme(o, nx, 0.0f);
        o[5] = NAN;
        fill_margin_pos(o + nx, n - nx);
        return;
    }
    if (kind == 10) { fill_wide_center(o, n); return; }
    if (kind == 11) { fill_halfmargin_pos(o, n); return; }
    if (kind == 19) { for (long i = 0; i < n; i++) o[i] = 100.0f; return; }
    if (kind == 20) { fill_ramp(o, n, 1.0e10, 1.0e10 + 4294966784.0); return; }
    if (kind == 21) {
        /* row 1: a ramp whose centered minimum quantizes to NULL_VALUE
         * (qlevel -0x1.528c8cp+0, NO_DITHER); row 2: 0..63 with a NaN, so
         * the image has ZBLANK */
        fill_ramp(o, nx, -0x1.b0a9fep+23, 0x1.51b436p+32);
        for (long i = nx; i < n; i++) o[i] = (float)(i - nx);
        o[nx + 5] = NAN;
        return;
    }

    for (long i = 0; i < n; i++) {
        long row = i / nx, col = i % nx;
        double v = 100.0 + 0.3 * col + 0.2 * row + 5.0 * lcg_gauss(&r);
        switch (kind) {
        case 1: o[i] = (i % 23 == 0) ? 0.0f : (float)v; break;
        case 2: o[i] = (i % 31 == 0) ? NAN : (float)v; break;
        case 3: o[i] = (row < 10) ? 42.5f : (float)v; break;
        case 4: o[i] = (row < 5) ? (float)(((i % 2) ? 1.0e20 : 0.0) + v) : (float)v; break;
        case 6: o[i] = (row < 5) ? (float)(((i % 2) ? 1.0e20 : 0.0) + v) : (float)v; break; /* huge rows, then specials */
        default: o[i] = (float)v; break;
        }
        if (kind == 12) {
            if (i % 37 == 5) o[i] = NAN;
            if (i % 41 == 7) o[i] = -1.0e30f;
        }
        if (kind == 13) {
            if (i % 200 == 3) o[i] = NAN;
            if (i % 200 == 10) o[i] = FLOATNULLVALUE;
        }
        if (kind == 16) {
            if (i % 23 == 0) o[i] = 0.0f;
            if (i % 31 == 0) o[i] = NAN;
        }
        if (kind == 17 || kind == 18) {
            if (i % 200 == 3) o[i] = NAN;
            if (i % 50 == 7) o[i] = kind == 17 ? -999.0f : DENORM_NULL;
            if (kind == 18 && i % 50 == 9) o[i] = -3.0e-42f;
        }
        if (kind == 7 || kind == 8) {
            if (i == 3) o[i] = NAN;
            if (i == 20) o[i] = -1.0e30f;
            if (kind == 8 && i == 10) o[i] = FLOATNULLVALUE;
        }
        if (kind == 5 || kind == 6) {
            /* IEEE specials. cfitsio reads the tile with fits_read_subset
             * and a null value (imcomp_compress_image), so fnan(): an all-
             * ones exponent (NaN, +Inf, -Inf) is a null, an all-zero
             * exponent (+-0, denormals) is read as 0. FLT_MIN (smallest
             * normal) is a real value. */
            if (i % 31 == 0) o[i] = ((i / 31) % 2) ? INFINITY : -INFINITY;
            else if (i % 37 == 0) o[i] = NAN;
            else if (i % 17 == 0) o[i] = ((i / 17) % 2) ? 1.0e-40f : -3.0e-42f;
            else if (i % 13 == 0) o[i] = -0.0f;
            else if (i % 19 == 0) o[i] = ((i / 19) % 2) ? FLT_MIN : -FLT_MIN;
        }
    }
}

/* float64 image content (img64/): kind 0 noise, 1 zeros, 2 NaN, 3 const
 * rows, 4 huge rows (unquantizable tiles), 5 IEEE specials, 6 huge +
 * specials, 7 NaN + real -1e30, 8 NaN + real DOUBLENULLVALUE + -1e30,
 * 14 a 1e8 background with 0.01 noise (needs double precision: float32
 * cannot hold it), 15 every pixel above FLT_MAX (1e39 + 1e33 noise: the
 * min/max scans must start from DBL_MAX). Computed in double, never rounded to float. The noise
 * is a sum of 4 uniforms (unit variance), not lcg_gauss: log/cos differ
 * in the last bit between libm and Go, which float32 rounding hides but a
 * float64 input would not; + - * are correctly rounded everywhere. */
static double lcg_noise4(lcg_t *r) {
    double s = lcg_uniform(r) + lcg_uniform(r) + lcg_uniform(r) + lcg_uniform(r);
    return (s - 2.0) * 1.7320508075688772;
}
static void gen_img_double(double *o, long n, long nx, int kind, uint64_t seed) {
    lcg_t r; lcg_init(&r, seed);
    for (long i = 0; i < n; i++) {
        long row = i / nx, col = i % nx;
        double v = 100.0 + 0.3 * col + 0.2 * row + 5.0 * lcg_noise4(&r);
        switch (kind) {
        case 1: o[i] = (i % 23 == 0) ? 0.0 : v; break;
        case 2: o[i] = (i % 31 == 0) ? NAN : v; break;
        case 3: o[i] = (row < 10) ? 42.5 : v; break;
        case 4: case 6: o[i] = (row < 5) ? ((i % 2) ? 1.0e20 : 0.0) + v : v; break;
        case 14: o[i] = 1.0e8 + 1.0e-4 * col + 0.01 * (v - 100.0 - 0.3 * col - 0.2 * row); break;
        case 15: o[i] = 1.0e39 + 1.0e33 * (v - 100.0); break;
        case 19: o[i] = 1.0e10 + (double)i * 1.0e-6; break;
        default: o[i] = v; break;
        }
        if (kind == 16) {
            if (i % 23 == 0) o[i] = 0.0;
            if (i % 31 == 0) o[i] = NAN;
        }
        if (kind == 17) {
            if (i % 200 == 3) o[i] = NAN;
            if (i % 50 == 7) o[i] = -999.0;
        }
        if (kind == 7 || kind == 8) {
            if (i == 3) o[i] = NAN;
            if (i == 20) o[i] = -1.0e30;
            if (kind == 8 && i == 10) o[i] = DOUBLENULLVALUE;
        }
        if (kind == 5 || kind == 6) {
            if (i % 31 == 0) o[i] = ((i / 31) % 2) ? INFINITY : -INFINITY;
            else if (i % 37 == 0) o[i] = NAN;
            else if (i % 17 == 0) o[i] = ((i / 17) % 2) ? 1.0e-310 : -4.0e-320;
            else if (i % 13 == 0) o[i] = -0.0;
            else if (i % 19 == 0) o[i] = ((i / 19) % 2) ? DBL_MIN : -DBL_MIN;
        }
    }
}

static const char *algo_name(int a) {
    switch (a) {
    case RICE_1: return "RICE_1";
    case GZIP_1: return "GZIP_1";
    case GZIP_2: return "GZIP_2";
    case HCOMPRESS_1: return "HCOMPRESS_1";
    case PLIO_1: return "PLIO_1";
    }
    return "?";
}

static int run_img_case(const img_case_t *c, const char *dir, int idx, const char *tag, uint64_t seedbase) {
    long n = 1;
    for (int d = 0; d < c->naxis; d++) n *= c->naxes[d];
    int isfloat = !strcmp(c->type, "f32");
    int isdouble = !strcmp(c->type, "f64");
    int width = type_width(c->type);
    unsigned char *raw = calloc((size_t)n, 8);
    double *dv = calloc((size_t)n, sizeof(double));
    long long *iv = calloc((size_t)n, sizeof(long long));
    float *fv = calloc((size_t)n, sizeof(float));
    void *data = NULL;
    uint64_t seed = seedbase + (uint64_t)idx;
    if (isfloat) {
        gen_img_float(fv, n, c->naxes[0], c->kind, seed);
        memcpy(raw, fv, (size_t)n * 4);
        data = fv;
    } else if (isdouble) {
        gen_img_double(dv, n, c->naxes[0], c->kind, seed);
        memcpy(raw, dv, (size_t)n * 8);
        data = dv;
    } else {
        gen_img_int(iv, n, c->naxes[0], c->type, c->kind, seed);
        data = malloc((size_t)n * 4);
        for (long i = 0; i < n; i++) {
            put_le(raw + (size_t)i * width, iv[i], width);
            switch (type_datatype(c->type)) {
            case TBYTE: ((unsigned char *)data)[i] = (unsigned char)iv[i]; break;
            case TSBYTE: ((signed char *)data)[i] = (signed char)iv[i]; break;
            case TSHORT: ((short *)data)[i] = (short)iv[i]; break;
            case TUSHORT: ((unsigned short *)data)[i] = (unsigned short)iv[i]; break;
            case TINT: ((int *)data)[i] = (int)iv[i]; break;
            default: ((unsigned int *)data)[i] = (unsigned int)iv[i]; break;
            }
        }
    }
    char path[1024], tmp[1100], cpath[1100], in_sha[65];
    sha_bytes_hex(raw, (size_t)n * (size_t)width, in_sha);

    int status = 0;
    fitsfile *in = NULL, *out = NULL;
    int cstatus = 0, lossless_verified = 0;
    double *decd = NULL; /* cfitsio's decode, recorded as a digest (lossy cases) */
    char *decnul = NULL;
    char errmsg[FLEN_ERRMSG] = "";
    if (c->hdr == 3) {
        /* straight into a compressed image, null value c->nullv (0: -1e30) */
        double nulv = c->nullv != 0 ? c->nullv : (double)NULLF;
        float fnulv = (float)nulv;
        float *wv = malloc((size_t)n * sizeof(float));
        double *wd = malloc((size_t)n * sizeof(double));
        for (long i = 0; i < n; i++) {
            if (isdouble) wd[i] = isnan(dv[i]) ? nulv : dv[i];
            else wv[i] = isnan(fv[i]) ? fnulv : fv[i];
        }
        snprintf(cpath, sizeof cpath, "!%s/%s.fits", dir, c->name);
        fits_create_file(&out, cpath, &status);
        fits_create_img(out, BYTE_IMG, 0, NULL, &status);
        fits_set_compression_type(out, c->algo, &status);
        fits_set_tile_dim(out, c->naxis, (long *)c->tile, &status);
        if (c->qlevel != 0) fits_set_quantize_level(out, c->qlevel, &status);
        if (c->qmethod != 0) fits_set_quantize_method(out, c->qmethod, &status);
        if (c->seed != 0) fits_set_dither_seed(out, c->seed, &status);
        if (isdouble) {
            fits_create_img(out, DOUBLE_IMG, c->naxis, (long *)c->naxes, &status);
            fits_write_imgnull(out, TDOUBLE, 1, n, wd, &nulv, &status);
        } else {
            fits_create_img(out, FLOAT_IMG, c->naxis, (long *)c->naxes, &status);
            fits_write_imgnull(out, TFLOAT, 1, n, wv, &fnulv, &status);
        }
        fits_close_file(out, &status);
        free(wv); free(wd);
        if (status) { fits_report_error(stderr, status); return -1; }
        goto decode;
    }

    /* 1. uncompressed source image, written by cfitsio */
    snprintf(tmp, sizeof tmp, "!%s/%s.src.%s", dir, c->name, c->hdr ? "fits" : "tmp");
    fits_create_file(&in, tmp, &status);
    fits_create_img(in, type_bitpix(c->type), c->naxis, (long *)c->naxes, &status);
    if (c->hdr == 2) {
        long blank = 5;
        fits_write_key(in, TLONG, "BLANK", &blank, "null value (illegal for BITPIX -32)", &status);
    } else if (c->hdr) {
        /* Fixed strings, not real checksums (fits_write_chksum stamps the
         * date into the comment, which would make the golden
         * unreproducible); cfitsio copies them opaquely. */
        int t = 1;
        fits_write_key(in, TLOGICAL, "BLOCKED", &t, "deprecated, renamed ZBLOCKED", &status);
        fits_write_key_str(in, "OBJECT", "NGC 1234", "target", &status);
        fits_write_key_dbl(in, "EXPTIME", 300.5, -6, "exposure (s)", &status);
        fits_write_key_str(in, "CHECKSUM", "AAAAAAAAAAAAAAAA", "fixed test value", &status);
        fits_write_key_str(in, "DATASUM", "1234567890", "fixed test value", &status);
        fits_write_history(in, "source image written by cref", &status);
    }
    fits_write_img(in, type_datatype(c->type), 1, n, data, &status);
    if (status) { fits_report_error(stderr, status); return -1; }
    if (c->hdr == 4) {
        /* SIMPLE = F: patch column 30 of the first card, then reopen */
        fits_close_file(in, &status);
        snprintf(tmp, sizeof tmp, "%s/%s.src.fits", dir, c->name);
        FILE *pf = fopen(tmp, "r+b");
        if (!pf || fseek(pf, 29, SEEK_SET) || fputc('F', pf) == EOF) { fprintf(stderr, "cref: patch %s\n", tmp); return -1; }
        fclose(pf);
        fits_open_file(&in, tmp, READONLY, &status);
        if (status) { fits_report_error(stderr, status); return -1; }
    }

    /* 2. compress it with cfitsio's own writer, as fpack does */
    snprintf(cpath, sizeof cpath, "!%s/%s.fits", dir, c->name);
    fits_create_file(&out, cpath, &status);
    fits_create_img(out, BYTE_IMG, 0, NULL, &status);
    fits_set_compression_type(out, c->algo, &status);
    if (c->tile[0] != 0) fits_set_tile_dim(out, c->naxis, (long *)c->tile, &status);
    if (c->qlevel != 0) fits_set_quantize_level(out, c->qlevel, &status);
    if (c->qmethod != 0) fits_set_quantize_method(out, c->qmethod, &status);
    if (c->seed != 0) fits_set_dither_seed(out, c->seed, &status);
    if (c->hscale != 0) fits_set_hcomp_scale(out, c->hscale, &status);
    if (status) { fits_report_error(stderr, status); return -1; }
    fits_img_compress(in, out, &cstatus);
    if (cstatus) {
        char m[FLEN_ERRMSG];
        /* keep the most specific (first) message cfitsio queued */
        if (fits_read_errmsg(m)) snprintf(errmsg, sizeof errmsg, "%s", m);
        while (fits_read_errmsg(m)) {}
    }
    int st2 = 0;
    fits_close_file(out, &st2);
    st2 = 0;
    if (c->hdr) fits_close_file(in, &st2);
    else fits_delete_file(in, &st2);
decode:
    if (cstatus) {
        /* cfitsio refused: leave no .fits behind */
        snprintf(cpath, sizeof cpath, "%s/%s.fits", dir, c->name);
        remove(cpath);
    } else {
        /* 3. cfitsio's decode of its own file */
        fitsfile *rd = NULL;
        snprintf(cpath, sizeof cpath, "%s/%s.fits", dir, c->name);
        fits_open_file(&rd, cpath, READONLY, &status);
        fits_movabs_hdu(rd, 2, NULL, &status);
        double *dec = calloc((size_t)n, sizeof(double));
        double nanv = NAN;
        int anynul = 0;
        if (c->naxis <= 3) {
            fits_read_img(rd, TDOUBLE, 1, n, &nanv, dec, &anynul, &status);
        } else {
            /* fits_read_img on a compressed image takes 1-3 dimensions
             * only (status 414); decompress into memory as funpack does
             * (fits_img_decompress) and read the plain image. */
            fitsfile *mem = NULL;
            fits_create_file(&mem, "mem://", &status);
            fits_img_decompress(rd, mem, &status);
            fits_read_img(mem, TDOUBLE, 1, n, &nanv, dec, &anynul, &status);
            int st3 = 0;
            fits_close_file(mem, &st3);
        }
        fits_close_file(rd, &status);
        if (status) { fits_report_error(stderr, status); return -1; }
        if (isfloat || isdouble || c->hscale != 0) {
            /* nulls (NaN) as null flags; a float32 image's values rounded
             * to float32, its logical pixel type (what the Go reader
             * returns; the TDOUBLE values are float32 arithmetic results
             * widened by cfitsio, so this loses nothing the Go side reads) */
            decd = malloc((size_t)n * sizeof(double));
            decnul = calloc((size_t)n, 1);
            for (long i = 0; i < n; i++) {
                if (isnan(dec[i])) { decnul[i] = 1; decd[i] = 0; continue; }
                decd[i] = isfloat ? (double)(float)dec[i] : dec[i];
            }
        } else {
            /* lossless integer case: cfitsio must give back the input */
            for (long i = 0; i < n; i++) {
                if (dec[i] != (double)iv[i]) {
                    fprintf(stderr, "cref: %s: cfitsio decode differs from input at %ld (%.17g vs %lld)\n",
                        c->name, i, dec[i], iv[i]);
                    return -1;
                }
            }
            lossless_verified = 1;
        }
        free(dec);
    }

    snprintf(path, sizeof path, "%s/%s.json", dir, c->name);
    FILE *jf = fopen(path, "w");
    if (!jf) return -1;
    fprintf(jf, "{\"name\": \"%s\", \"type\": \"%s\", \"naxes\": [", c->name, c->type);
    for (int d = 0; d < c->naxis; d++) fprintf(jf, "%s%ld", d ? ", " : "", c->naxes[d]);
    fprintf(jf, "], \"algo\": \"%s\", \"tile\": [", algo_name(c->algo));
    if (c->tile[0] != 0) for (int d = 0; d < c->naxis; d++) fprintf(jf, "%s%ld", d ? ", " : "", c->tile[d]);
    fprintf(jf, "], \"qlevel\": ");
    fprint_double(jf, (double)c->qlevel);
    fprintf(jf, ", \"qmethod\": %d, \"seed\": %d, \"hscale\": ", c->qmethod, c->seed);
    fprint_double(jf, (double)c->hscale);
    fprintf(jf, ", \"lossless_verified\": %s", lossless_verified ? "true" : "false");
    if (c->hdr == 3) {
        fprintf(jf, ", \"write_imgnull\": true");
        if (c->nullv != 0) { fprintf(jf, ", \"nullvalue\": "); fprint_double(jf, c->nullv); }
    }
    else if (c->hdr) fprintf(jf, ", \"src_header\": true");
    fprintf(jf, ", \"status\": %d, \"error\": \"", cstatus);
    for (const char *p = errmsg; *p; p++) if (*p != '"' && *p != '\\') fputc(*p, jf);
    fprintf(jf, "\", \"kind\": %d, \"gen_seed\": %llu, \"in_sha256\": \"%s\"", c->kind, (unsigned long long)seed, in_sha);
    if (decd) {
        fprintf(jf, ",\n \"dec\": ");
        jvalues(jf, K_F64, decd, decnul, n, c->naxis, c->naxes);
        free(decd); free(decnul);
    }
    fprintf(jf, "}\n");
    fclose(jf);
    fprintf(manifest, "%s %s\n", tag, c->name);
    free(raw); free(iv); free(fv); free(dv);
    if (!isfloat && !isdouble) free(data);
    return 0;
}

/* cfitsio's default HCOMPRESS_1 tile shape for a range of image heights:
 * fits_create_img on a compressing mem:// file runs imcomp_init_table,
 * which writes ZTILE1/ZTILE2; no pixels are needed. One line per height:
 * "naxis1 naxis2 ztile1 ztile2". */
static int run_hcomp_default_tiles(const char *path) {
    FILE *o = fopen(path, "w");
    if (!o) { fprintf(stderr, "cref: open %s: %s\n", path, strerror(errno)); return -1; }
    static const long extra[] = {1681, 8401, 13443, 166323};
    long nh = 2000 - 4 + 1 + (long)(sizeof extra / sizeof extra[0]);
    for (long k = 0; k < nh; k++) {
        long naxes[2] = {8, k < 1997 ? 4 + k : extra[k - 1997]};
        int status = 0;
        fitsfile *f = NULL;
        long zt1 = 0, zt2 = 0;
        fits_create_file(&f, "mem://", &status);
        fits_create_img(f, BYTE_IMG, 0, NULL, &status);
        fits_set_compression_type(f, HCOMPRESS_1, &status);
        fits_create_img(f, SHORT_IMG, 2, naxes, &status);
        fits_read_key(f, TLONG, "ZTILE1", &zt1, NULL, &status);
        fits_read_key(f, TLONG, "ZTILE2", &zt2, NULL, &status);
        if (status) { fits_report_error(stderr, status); fclose(o); return -1; }
        int st2 = 0;
        fits_delete_file(f, &st2);
        fprintf(o, "%ld %ld %ld %ld\n", naxes[0], naxes[1], zt1, zt2);
    }
    fclose(o);
    return 0;
}

static int run_encode(const char *outdir) {
    char dir[256];
    snprintf(dir, sizeof dir, "%s/MANIFEST", outdir);
    manifest = fopen(dir, "w");
    if (!manifest) { fprintf(stderr, "cref: open %s: %s\n", dir, strerror(errno)); return 1; }
    snprintf(dir, sizeof dir, "%s/rice", outdir);
    for (int i = 0; i < n_rice_cases; i++) if ((!only || !strcmp(only, rice_cases[i].name)) && run_rice_case(&rice_cases[i], dir) != 0) return 1;
    snprintf(dir, sizeof dir, "%s/hcomp", outdir);
    for (int i = 0; i < n_hc_cases; i++) if ((!only || !strcmp(only, hc_cases[i].name)) && run_hc_case(&hc_cases[i], dir) != 0) return 1;
    snprintf(dir, sizeof dir, "%s/plio", outdir);
    for (int i = 0; i < n_plio_cases; i++) if ((!only || !strcmp(only, plio_cases[i].name)) && run_plio_case(&plio_cases[i], dir) != 0) return 1;
    snprintf(dir, sizeof dir, "%s/img", outdir);
    for (int i = 0; i < n_img_cases; i++) if ((!only || !strcmp(only, img_cases[i].name)) && run_img_case(&img_cases[i], dir, i, "img", 5000) != 0) return 1;
    snprintf(dir, sizeof dir, "%s/imgw", outdir);
    for (int i = 0; i < n_imgw_cases; i++) if ((!only || !strcmp(only, imgw_cases[i].name)) && run_img_case(&imgw_cases[i], dir, i, "imgw", 7000) != 0) return 1;
    snprintf(dir, sizeof dir, "%s/img64", outdir);
    for (int i = 0; i < n_img64_cases; i++) if ((!only || !strcmp(only, img64_cases[i].name)) && run_img_case(&img64_cases[i], dir, i, "img64", 8000) != 0) return 1;
    snprintf(dir, sizeof dir, "%s/img2", outdir);
    for (int i = 0; i < n_img2_cases; i++) if ((!only || !strcmp(only, img2_cases[i].name)) && run_img_case(&img2_cases[i], dir, i, "img2", 9000) != 0) return 1;
    fclose(manifest);
    if (!only) {
        snprintf(dir, sizeof dir, "%s/hcomp_default_tiles.txt", outdir);
        if (run_hcomp_default_tiles(dir) != 0) return 1;
    }
    printf("cref encode: %d rice, %d hcomp, %d plio, %d img, %d imgw, %d img64, %d img2 cases -> %s\n",
        n_rice_cases, n_hc_cases, n_plio_cases, n_img_cases, n_imgw_cases, n_img64_cases, n_img2_cases, outdir);
    return 0;
}

/* ======================================================================
 * "decode" mode: cref decode <file> prints cfitsio's reading of the first
 * HDU holding an image with NAXIS > 0 as JSON:
 *   {"cfitsio": "4.7.0", "naxes": [...], "pixels": [v, ...]}
 * Each pixel is fits_read_img(TDOUBLE, nulval = NaN) with BSCALE/BZERO
 * applied -- the read the img goldens' "dec" digests use; a null pixel
 * (NaN) is printed as null, +-Inf as +-DBL_MAX (cfitsio itself reads Inf
 * as null, so this does not occur in practice). Values are printed with
 * %.17g (round-trip exact). Used by the Go interop tests to judge files
 * the Go writer produced. (testdata/cfref "dump" also reads TLONG, which
 * makes libcfitsio 4.7.0 crash on HCOMPRESS_1 images with BZERO; this
 * mode reads TDOUBLE only.)
 *
 * "decodenull" mode reads with fits_read_imgnull (null flags) instead.
 * libcfitsio 4.7.0 gets those flags wrong for float images with more than
 * one GZIP_COMPRESSED_DATA fallback tile holding NaNs: a tile keeps the
 * null flags of the earlier fallback tiles at the same offsets, so real
 * pixels are flagged null (golden img/rice_f32_nan_realnull_rows, also
 * on cfitsio's own file; TestCfitsioImgnullStaleFallbackFlags). fits_read_img
 * with a null value is right there, which is why "decode" uses it.
 * ====================================================================== */
static int run_decode(const char *path, int imgnull) {
    fitsfile *f = NULL;
    int st = 0, naxis = 0, anyn = 0;
    long naxes[9] = {0};
    fits_open_file(&f, path, READONLY, &st);
    /* first HDU with an image of NAXIS > 0 (a funpacked ZSIMPLE image is
     * the primary array; otherwise HDU 2) */
    for (int h = 1; !st; h++) {
        int ht = 0;
        if (fits_movabs_hdu(f, h, &ht, &st)) break;
        if (ht == IMAGE_HDU && !fits_get_img_dim(f, &naxis, &st) && naxis > 0) break;
    }
    fits_get_img_dim(f, &naxis, &st);
    if (naxis > 9) naxis = 9;
    fits_get_img_size(f, naxis, naxes, &st);
    if (st) { fits_report_error(stderr, st); return 1; }
    long n = 1;
    for (int i = 0; i < naxis; i++) n *= naxes[i];
    double *d = malloc((size_t)n * sizeof(double));
    char *nul = calloc((size_t)n, 1);
    if (!d || !nul) { fprintf(stderr, "cref: OOM\n"); return 1; }
    if (imgnull) {
        fits_read_imgnull(f, TDOUBLE, 1, n, d, nul, &anyn, &st);
    } else {
        double nanv = NAN;
        fits_read_img(f, TDOUBLE, 1, n, &nanv, d, &anyn, &st);
    }
    fits_close_file(f, &st);
    if (st) { fits_report_error(stderr, st); return 1; }
    float ver = 0;
    fits_get_version(&ver);
    printf("{\"cfitsio\": \"%.3f\", \"naxes\": [", ver);
    for (int i = 0; i < naxis; i++) printf("%s%ld", i ? ", " : "", naxes[i]);
    printf("], \"pixels\": [");
    for (long i = 0; i < n; i++) {
        if (i) putchar(',');
        if (nul[i] || isnan(d[i])) printf("null");
        else if (isinf(d[i])) printf("%s1.7976931348623157e308", d[i] < 0 ? "-" : ""); /* +-DBL_MAX: JSON has no Inf */
        else printf("%.17g", d[i]);
    }
    printf("]}\n");
    free(d); free(nul);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s gen <output_dir> | %s encode <output_dir> | %s list\n", argv[0], argv[0], argv[0]);
        return 2;
    }
    if (strcmp(argv[1], "list") == 0) {
        for (int i = 0; i < ncases; i++) printf("%s\n", cases[i].name);
        return 0;
    }
    if (strcmp(argv[1], "decode") == 0 && argc >= 3) return run_decode(argv[2], 0);
    if (strcmp(argv[1], "decodenull") == 0 && argc >= 3) return run_decode(argv[2], 1);
    if (strcmp(argv[1], "encode") == 0 && argc >= 3) {
        if (argc >= 4) only = argv[3];
        return run_encode(argv[2]);
    }
    if (strcmp(argv[1], "gen") != 0 || argc < 3) {
        fprintf(stderr, "usage: %s gen <output_dir> | %s encode <output_dir>\n", argv[0], argv[0]);
        return 2;
    }
    const char *outdir = argv[2];
    printf("cref: %d test cases -> %s\n", ncases, outdir);
    for (int i = 0; i < ncases; i++) {
        if (run_case(&cases[i], outdir) != 0) return 1;
    }
    /* The dither random table itself: fits_rand_value[0..N_RANDOM-1] as
     * built by libcfitsio's fits_init_randoms, each float32 widened to
     * double, as an exact digest (fits_rand_value.json). */
    extern float *fits_rand_value;
    if (fits_init_randoms() != 0) { fprintf(stderr, "cref: fits_init_randoms failed\n"); return 1; }
    char rpath[1024];
    snprintf(rpath, sizeof rpath, "%s/fits_rand_value.json", outdir);
    FILE *rf = fopen(rpath, "w");
    if (!rf) { fprintf(stderr, "cref: open %s: %s\n", rpath, strerror(errno)); return 1; }
    double *rd = malloc((size_t)N_RANDOM * sizeof(double));
    for (int i = 0; i < N_RANDOM; i++) rd[i] = (double)fits_rand_value[i];
    fprintf(rf, "{\"n_random\": %d, \"values\": ", N_RANDOM);
    jvalues(rf, K_F64, rd, NULL, N_RANDOM, 0, NULL);
    fprintf(rf, "}\n");
    fclose(rf);
    free(rd);
    printf("  fits_rand_value: %d floats\n", N_RANDOM);
    return 0;
}
