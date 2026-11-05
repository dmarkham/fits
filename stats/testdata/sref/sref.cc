/*
 * sref.cc -- Siril reference harness for the Go stats package.
 *
 * Everything that produces an expected value is Siril's own C/C++ code,
 * copied verbatim (see extract_siril.py, PROVENANCE.md):
 *
 *   findMinMaxPercentile      rt_algo_excerpt.cc  (Siril src/rt/rt_algo.cc)
 *   histogram_median_float    siril_excerpt.c     (src/algos/sorting.c)
 *   quickmedian_float         siril_excerpt.c     (src/algos/sorting.c)
 *   siril_stats_float_sd      siril_excerpt.c     (src/algos/statistics_float.c)
 *   siril_stats_float_mad     siril_excerpt.c     (src/algos/statistics_float.c)
 *   reassign_to_non_null_data_float              (src/algos/statistics_float.c)
 *
 * The Go generic (non-float32) percentile path runs the same algorithm in
 * float64. Its oracle is the same verbatim findMinMaxPercentile source
 * compiled a second time with float replaced by double (the F64 include
 * below); nothing else about the code changes.
 *
 * This file only generates inputs, calls those functions and prints the
 * results as JSON. Every float is printed as its IEEE bit pattern in hex
 * so the Go test can compare bit-for-bit.
 *
 * Build:   make
 * Use:     ./sref gen <output_dir>      write golden/siril_stats_golden.json
 *          ./sref median-f32 < hex      Siril median of float32 bit patterns
 *          ./sref median-f64 < hex      float64 findMinMaxPercentile median
 * The median modes are used by ../gen_golden.py for the sigma-clip centre.
 */

#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>

namespace rtengine {
    template<typename T>
    constexpr const T& LIM(const T& val, const T& low, const T& high) {
        return val < low ? low : (val > high ? high : val);
    }
}

// float32 findMinMaxPercentile, verbatim (C linkage, used by Siril's C code).
#include "rt_algo_excerpt.cc"

// The same source compiled at double precision for the Go generic path.
#define float double
#define findMinMaxPercentile findMinMaxPercentile_f64
#include "rt_algo_excerpt.cc"
#undef findMinMaxPercentile
#undef float

#include "siril_stubs.h"
#include "sha256.h"

// ----------------------------------------------------------------
// Bits helpers
// ----------------------------------------------------------------

static uint32_t f32bits(float f) { uint32_t u; memcpy(&u, &f, 4); return u; }
static uint64_t f64bits(double d) { uint64_t u; memcpy(&u, &d, 8); return u; }
static float bitsf32(uint32_t u) { float f; memcpy(&f, &u, 4); return f; }
static double bitsf64(uint64_t u) { double d; memcpy(&d, &u, 8); return d; }

// ----------------------------------------------------------------
// Deterministic formula inputs. The Go test (siril_test.go, formulaData)
// implements the identical integer recurrence, so these large inputs
// need not be committed. Only exactly-rounded IEEE operations are used:
// an integer times 2^-24, and one float division of two integers.
// ----------------------------------------------------------------

static uint32_t xorshift32(uint32_t &x) {
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return x;
}

static std::vector<float> formula(const std::string &kind, uint32_t seed, size_t n) {
    std::vector<float> out(n);
    uint32_t x = seed;
    for (size_t i = 0; i < n; i++) {
        uint32_t b = xorshift32(x);
        if (kind == "uniform24") {
            out[i] = (float)(b >> 8) * 5.9604644775390625e-8f; // exact: 24-bit int * 2^-24
        } else if (kind == "astro_u16") {
            uint32_t adu = 2900 + (b & 0xff) + ((b >> 8) & 0xff) + ((b >> 16) & 0xff) + (b >> 24);
            if (b % 509 == 0) adu += (b >> 7) % 60000;
            if (adu > 65535) adu = 65535;
            out[i] = (float)adu / 65535.0f;
        } else {
            fprintf(stderr, "unknown formula %s\n", kind.c_str());
            exit(1);
        }
    }
    return out;
}

// Integer-valued and float64 formula inputs for the generic (non-float32)
// path, as doubles. Every value is an integer or an integer times a power
// of two that a double holds exactly, so the Go test converts it to the
// target element type without rounding (convertSlice checks this).
//
//   uniform24   float64 of the float32 uniform24 value
//   astro_adu   the astro_u16 ADU integer (before the division)
//   u16         b >> 16                      (uint16)
//   i8          the low byte as a signed int  (int8, -128..127)
//   u32         b                             (uint32, up to 2^32-1)
//   u64         odd b: 2^63 + b*2^11, even b: b*2^20   (uint64, above 2^63)
//   i64         (int32)b * 2^30               (int64, |v| up to 2^61)
//   uniform53   53 bits from two draws times 2^-53     (float64 in [0,1))
static std::vector<double> formula64(const std::string &kind, uint32_t seed, size_t n) {
    if (kind == "uniform24") {
        auto d32 = formula(kind, seed, n);
        return std::vector<double>(d32.begin(), d32.end());
    }
    std::vector<double> out(n);
    uint32_t x = seed;
    for (size_t i = 0; i < n; i++) {
        uint32_t b = xorshift32(x);
        if (kind == "astro_adu") {
            uint32_t adu = 2900 + (b & 0xff) + ((b >> 8) & 0xff) + ((b >> 16) & 0xff) + (b >> 24);
            if (b % 509 == 0) adu += (b >> 7) % 60000;
            if (adu > 65535) adu = 65535;
            out[i] = (double)adu;
        } else if (kind == "u16") {
            out[i] = (double)(b >> 16);
        } else if (kind == "i8") {
            int v = (int)(b & 0xff);
            if (v >= 128) v -= 256;
            out[i] = (double)v;
        } else if (kind == "u32") {
            out[i] = (double)b;
        } else if (kind == "u64") {
            uint64_t v = (b & 1) ? (UINT64_C(1) << 63) + ((uint64_t)b << 11) : (uint64_t)b << 20;
            out[i] = (double)v;
        } else if (kind == "i64") {
            int64_t v = (b >= UINT32_C(0x80000000)) ? (int64_t)b - (INT64_C(1) << 32) : (int64_t)b;
            out[i] = (double)(v * (INT64_C(1) << 30));
        } else if (kind == "uniform53") {
            uint32_t b2 = xorshift32(x);
            uint64_t m = (((uint64_t)b << 32) | b2) >> 11;
            out[i] = (double)m * 0x1p-53;
        } else {
            fprintf(stderr, "unknown formula64 %s\n", kind.c_str());
            exit(1);
        }
    }
    return out;
}

// ----------------------------------------------------------------
// Input generators for the "gen" cases (and the small committed ones).
// Portable: xorshift32 (not std::*_distribution, whose output is
// implementation-defined) and exactly specified IEEE arithmetic, so any
// C++ compiler and standard library reproduces the file, and the Go test
// (genValues32 / genValues64 in siril_test.go) regenerates the same bits;
// in_sha256 proves it.
// ----------------------------------------------------------------

// A generator state from a small seed: the seed is spread by a
// multiplicative hash and the first draws, which are poorly mixed for
// small states, are discarded.
static uint32_t mixed_state(uint32_t seed) {
    uint32_t x = seed * UINT32_C(2654435761);
    if (x == 0) x = 1;
    for (int i = 0; i < 8; i++) xorshift32(x);
    return x;
}

// Uniform on [lo, hi): lo + (hi - lo) * u with u a multiple of 2^-24 in
// [0, 1), each operation rounded once in float.
static std::vector<float> gen_uniform(uint32_t seed, size_t n, float lo, float hi) {
    uint32_t x = mixed_state(seed);
    std::vector<float> out(n);
    for (auto &v : out) {
        float u = (float)(xorshift32(x) >> 8) * 5.9604644775390625e-8f;
        v = lo + (hi - lo) * u;
    }
    return out;
}

static std::vector<float> gen_gradient(int n, float lo, float hi) {
    std::vector<float> out(n);
    for (int i = 0; i < n; i++)
        out[i] = lo + (hi - lo) * float(i) / float(n - 1);
    return out;
}

// Gen specs: the JSON object the Go test reads back (floats as IEEE bit
// patterns in hex, so nothing is lost in decimal).
static std::string spec32(const char *kind, uint32_t seed, float lo, float hi) {
    char b[160];
    snprintf(b, sizeof(b), "{\"kind\":\"%s\",\"seed\":%" PRIu32 ",\"lo\":\"%08" PRIx32 "\",\"hi\":\"%08" PRIx32 "\"}",
             kind, seed, f32bits(lo), f32bits(hi));
    return b;
}

// Ties: 0.1f, 0.2f, 0.3f repeating with stride 7.
static std::vector<float> gen_ties3(size_t n) {
    static const float three[] = {0.1f, 0.2f, 0.3f};
    std::vector<float> d(n);
    for (size_t i = 0; i < n; i++) d[i] = three[(i * 7) % 3];
    return d;
}

// 30000 + k/4, k = 0..40: ADU levels in quarter steps (exact).
static std::vector<float> gen_adu_quarter(uint32_t seed, size_t n) {
    uint32_t x = mixed_state(seed);
    std::vector<float> out(n);
    for (auto &v : out) v = 30000.0f + (float)(xorshift32(x) % 41) * 0.25f;
    return out;
}

// ----------------------------------------------------------------
// JSON output
// ----------------------------------------------------------------

static const float PS32[] = {0.0f, 0.01f, 0.1f, 0.25f, 0.5f, 0.75f, 0.9f, 0.99f, 1.0f};
static const double PS64[] = {0.0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0};
static const char *PNAMES[] = {"0", "0.01", "0.1", "0.25", "0.5", "0.75", "0.9", "0.99", "1"};
static const int NPS = 9;

// A case's input is committed ("input") unless it has a formula (the
// xorshift32 kinds above) or a gen spec (the committed-input generators
// below, mirrored in Go by genValues32/genValues64): then only the spec is
// written. Inputs and non-null lists of more than DIGEST_MIN values that
// are not committed are identified by sha256 over their big-endian bit
// patterns ("in_sha256", "nonnull_sha256"), so the Go test proves it
// regenerated the same input and that FilterNonZero kept the same values.
static const size_t DIGEST_MIN = 64;

struct F32Case {
    std::string name;
    std::vector<float> data;
    std::string formula; // formula kind, or empty
    uint32_t seed;
    std::string gen;     // gen spec (a JSON object), or empty
};

static std::string sha32(const std::vector<float> &d) {
    Sha256 h;
    for (float v : d) h.be32(f32bits(v));
    return h.hex();
}

static std::string sha64(const std::vector<double> &d) {
    Sha256 h;
    for (double v : d) h.be64(f64bits(v));
    return h.hex();
}

static void hexlist32(FILE *f, const std::vector<float> &d) {
    fputc('"', f);
    for (size_t i = 0; i < d.size(); i++)
        fprintf(f, "%s%08" PRIx32, i ? " " : "", f32bits(d[i]));
    fputc('"', f);
}

static void write_f32_case(FILE *f, const F32Case &c, bool last) {
    const std::vector<float> &raw = c.data;
    const size_t n = raw.size();

    // Go skips NaN in Percentile/Median/MAD; Siril's callers strip them
    // first. Evaluate Siril on the NaN-stripped data.
    std::vector<float> clean;
    for (float v : raw) if (!std::isnan(v)) clean.push_back(v);
    const size_t nc = clean.size();

    const bool committed = c.formula.empty() && c.gen.empty();
    const bool digest = !committed && n > DIGEST_MIN;
    fprintf(f, "{\"name\":\"%s\",\"n\":%zu,", c.name.c_str(), n);
    if (committed) {
        fprintf(f, "\"input\":");
        hexlist32(f, raw);
        fprintf(f, ",");
    } else if (!c.formula.empty()) {
        fprintf(f, "\"formula\":{\"kind\":\"%s\",\"seed\":%" PRIu32 "},", c.formula.c_str(), c.seed);
    } else {
        fprintf(f, "\"gen\":%s,", c.gen.c_str());
    }
    if (digest) fprintf(f, "\"in_sha256\":\"%s\",", sha32(raw).c_str());

    // siril_stats_float_sd on the raw data (no NaN filtering, as Siril).
    {
        float mean = 0;
        float sd = siril_stats_float_sd(raw.data(), (int)n, &mean);
        fprintf(f, "\"sd\":{\"mean\":\"%08" PRIx32 "\",\"stdev\":\"%08" PRIx32 "\"},", f32bits(mean), f32bits(sd));
    }

    // reassign_to_non_null_data_float: the Siril non-zero/non-NaN filter,
    // as a list for committed inputs, as a count and digest for the other
    // inputs above DIGEST_MIN values.
    if (committed || digest) {
        size_t good = 0;
        for (float v : raw) if (v != 0.f && !std::isnan(v)) good++;
        std::vector<float> nn;
        if (good > 0) {
            std::vector<float> copy = raw;
            float *out = sref_reassign_to_non_null_data_float(copy.data(), n, good);
            nn.assign(out, out + good);
            free(out);
        }
        if (committed) {
            fprintf(f, "\"nonnull\":");
            hexlist32(f, nn);
            fprintf(f, ",");
        } else {
            fprintf(f, "\"nonnull_n\":%zu,\"nonnull_sha256\":\"%s\",", nn.size(), sha32(nn).c_str());
        }
    }

    if (nc > 0) {
        fprintf(f, "\"pct\":{");
        for (int i = 0; i < NPS; i++) {
            float v;
            findMinMaxPercentile(clean.data(), nc, PS32[i], &v, PS32[i], &v, 1);
            fprintf(f, "%s\"%s\":\"%08" PRIx32 "\"", i ? "," : "", PNAMES[i], f32bits(v));
        }
        fprintf(f, "},");
        std::vector<float> work = clean;
        double med = histogram_median_float(work.data(), nc, SINGLE_THREADED);
        double mad = siril_stats_float_mad(clean.data(), nc, med, SINGLE_THREADED, NULL);
        work = clean;
        double qmed = quickmedian_float(work.data(), nc);
        fprintf(f, "\"median\":\"%08" PRIx32 "\",", f32bits((float)med));
        fprintf(f, "\"mad\":\"%08" PRIx32 "\",", f32bits((float)mad));
        fprintf(f, "\"quickmedian\":\"%016" PRIx64 "\",", f64bits(qmed));
    }
    fprintf(f, "\"nclean\":%zu}%s\n", nc, last ? "" : ",");
}

// ---- generic (float64) path ----

struct F64Case {
    std::string name;
    std::string dtype; // Go element type the test converts to
    std::vector<double> data;
    std::string formula; // formula64(kind, seed, n) kind, or empty
    uint32_t seed;
    std::string gen;     // gen spec (a JSON object), or empty
};

static double median_f64(std::vector<double> d) {
    double med;
    findMinMaxPercentile_f64(d.data(), d.size(), 0.5, &med, 0.5, &med, 1);
    return med;
}

static void write_f64_case(FILE *f, const F64Case &c, bool last) {
    const auto &d = c.data;
    const bool committed = c.formula.empty() && c.gen.empty();
    fprintf(f, "{\"name\":\"%s\",\"dtype\":\"%s\",\"n\":%zu,", c.name.c_str(), c.dtype.c_str(), d.size());
    if (committed) {
        fprintf(f, "\"input\":\"");
        for (size_t i = 0; i < d.size(); i++)
            fprintf(f, "%s%016" PRIx64, i ? " " : "", f64bits(d[i]));
        fprintf(f, "\",");
    } else if (!c.formula.empty()) {
        fprintf(f, "\"formula\":{\"kind\":\"%s\",\"seed\":%" PRIu32 "},", c.formula.c_str(), c.seed);
    } else {
        fprintf(f, "\"gen\":%s,", c.gen.c_str());
    }
    if (!committed && d.size() > DIGEST_MIN) fprintf(f, "\"in_sha256\":\"%s\",", sha64(d).c_str());
    fprintf(f, "\"pct\":{");
    for (int i = 0; i < NPS; i++) {
        double v;
        findMinMaxPercentile_f64(d.data(), d.size(), PS64[i], &v, PS64[i], &v, 1);
        fprintf(f, "%s\"%s\":\"%016" PRIx64 "\"", i ? "," : "", PNAMES[i], f64bits(v));
    }
    double med = median_f64(d);
    // MAD in float64: the siril_stats_float_mad expression, in double.
    std::vector<double> tmp(d.size());
    for (size_t i = 0; i < d.size(); i++) tmp[i] = std::fabs(d[i] - med);
    double mad = median_f64(tmp);
    fprintf(f, "},\"median\":\"%016" PRIx64 "\",\"mad\":\"%016" PRIx64 "\"}%s\n",
            f64bits(med), f64bits(mad), last ? "" : ",");
}

// Integers lo..hi (b mod the range size; the slight modulo bias is
// irrelevant here), converted to T.
template <typename T>
static std::vector<double> gen_ints(uint32_t seed, size_t n, long lo, long hi) {
    uint32_t x = mixed_state(seed);
    const uint64_t span = (uint64_t)(hi - lo) + 1;
    std::vector<double> out(n);
    for (auto &v : out) v = (double)(T)(lo + (long)(xorshift32(x) % span));
    return out;
}

// Uniform on [lo, hi): lo + (hi - lo) * u with u = 53 random bits * 2^-53.
static std::vector<double> gen_f64(uint32_t seed, size_t n, double lo, double hi) {
    uint32_t x = mixed_state(seed);
    std::vector<double> out(n);
    for (auto &v : out) {
        uint32_t a = xorshift32(x), b = xorshift32(x);
        double u = (double)((((uint64_t)a << 32) | b) >> 11) * 0x1p-53;
        v = lo + (hi - lo) * u;
    }
    return out;
}

static int gen(const char *outdir) {
    if (!sha256_selftest()) { fprintf(stderr, "sha256 self-test failed\n"); return 1; }
    std::vector<F32Case> cs;
    auto add = [&](const std::string &name, std::vector<float> d) { cs.push_back({name, std::move(d), "", 0, ""}); };
    auto addgen = [&](const std::string &name, std::vector<float> d, std::string spec) {
        cs.push_back({name, std::move(d), "", 0, std::move(spec)});
    };

    addgen("gradient_200", gen_gradient(200, 0.01f, 0.30f), spec32("gradient", 0, 0.01f, 0.30f));
    addgen("gradient_1000", gen_gradient(1000, 0.0f, 1.0f), spec32("gradient", 0, 0.0f, 1.0f));
    addgen("random_500", gen_uniform(42, 500, 0.0f, 1.0f), spec32("uniform", 42, 0.0f, 1.0f));
    // random_10000 and astro_u16_5000 are formula inputs (below).
    {
        // A gradient whose first 60 values are zeroed.
        auto d = gen_gradient(300, 0.01f, 0.25f);
        for (int i = 0; i < 60; i++) d[i] = 0.0f;
        addgen("with_zeros", d, spec32("gradient_zero60", 0, 0.01f, 0.25f));
    }
    add("all_same", std::vector<float>(100, 0.42f));
    add("single", {3.14f});
    add("empty", {});
    add("n2", gen_uniform(2, 2, 0.0f, 1.0f));
    add("n3", gen_uniform(3, 3, 0.0f, 1.0f));
    add("n8", gen_uniform(8, 8, 0.0f, 1.0f));
    add("n9", gen_uniform(9, 9, 0.0f, 1.0f));
    add("n10", gen_uniform(10, 10, 0.0f, 1.0f));
    addgen("negative_301", gen_uniform(301, 301, -1.0f, 1.0f), spec32("uniform", 301, -1.0f, 1.0f));
    addgen("ties_1001", gen_ties3(1001), "{\"kind\":\"ties3\"}");
    {
        // Uniform with every 9th value (from the first) NaN.
        auto d = gen_uniform(77, 200, 0.0f, 1.0f);
        for (size_t i = 0; i < d.size(); i += 9) d[i] = NAN;
        addgen("with_nan_200", d, spec32("uniform_nan9", 77, 0.0f, 1.0f));
    }
    add("mixed_zero_nan", {0.0f, -0.0f, 0.5f, NAN, 1e-40f, -0.25f, 0.0f, 3.0f, NAN, -0.0f, 0.125f});
    add("all_zero", std::vector<float>(7, 0.0f));
    addgen("adu_offset_1000", gen_adu_quarter(1000, 1000), "{\"kind\":\"adu_quarter\",\"seed\":1000}");
    // Many small random arrays, where siril_stats_float_sd's float32
    // roundings matter (siril.go's MeanStdevSiril doc counts them). They
    // are uniform24 formula inputs (not committed): one xorshift32 stream
    // gives each case its length (2..33) and its formula seed.
    {
        uint32_t x = mixed_state(2024);
        for (int k = 0; k < 400; k++) {
            int n = 2 + (int)(xorshift32(x) % 32);
            uint32_t seed = xorshift32(x);
            char name[32];
            snprintf(name, sizeof(name), "small_%03d", k);
            cs.push_back({name, formula("uniform24", seed, n), "uniform24", seed, ""});
        }
    }
    // Formula inputs (not committed; regenerated in Go).
    cs.push_back({"random_10000", formula("uniform24", 99, 10000), "uniform24", 99, ""});
    cs.push_back({"astro_u16_5000", formula("astro_u16", 5, 5000), "astro_u16", 5, ""});
    cs.push_back({"f_astro_65537", formula("astro_u16", 11, 65537), "astro_u16", 11, ""});
    cs.push_back({"f_uniform_16777219", formula("uniform24", 12, 16777219), "uniform24", 12, ""});
    cs.push_back({"f_uniform_20000000", formula("uniform24", 13, 20000000), "uniform24", 13, ""});
    cs.push_back({"f_astro_20000000", formula("astro_u16", 14, 20000000), "astro_u16", 14, ""});

    std::vector<F64Case> gs;
    auto ints = [](uint32_t seed, long lo, long hi) {
        char b[96];
        snprintf(b, sizeof(b), "{\"kind\":\"ints\",\"seed\":%" PRIu32 ",\"lo\":%ld,\"hi\":%ld}", seed, lo, hi);
        return std::string(b);
    };
    gs.push_back({"int16_random_500", "int16", gen_ints<int16_t>(1, 500, 0, 30000), "", 0, ints(1, 0, 30000)});
    gs.push_back({"uint8_random_500", "uint8", gen_ints<uint8_t>(2, 500, 0, 255), "", 0, ints(2, 0, 255)});
    gs.push_back({"int32_negative_200", "int32", gen_ints<int32_t>(3, 200, -1000, 1000), "", 0, ints(3, -1000, 1000)});
    auto f64 = [&](const std::string &name, const std::string &dtype, const std::string &kind, uint32_t seed, size_t n) {
        gs.push_back({name, dtype, formula64(kind, seed, n), kind, seed, ""});
    };
    auto put = [&](const std::string &name, const std::string &dtype, std::vector<double> d) {
        gs.push_back({name, dtype, std::move(d), "", 0, ""});
    };
    f64("uint16_5000", "uint16", "u16", 4, 5000);
    gs.push_back({"int64_100", "int64", gen_ints<int64_t>(5, 100, -100000, 100000), "", 0, ints(5, -100000, 100000)});
    f64("float64_random_1000", "float64", "uniform53", 6, 1000);
    {
        char b[160];
        snprintf(b, sizeof(b), "{\"kind\":\"uniform\",\"seed\":7,\"lo\":\"%016" PRIx64 "\",\"hi\":\"%016" PRIx64 "\"}",
                 f64bits(1e6), f64bits(1e6 + 1.0));
        gs.push_back({"float64_offset_777", "float64", gen_f64(7, 777, 1e6, 1e6 + 1.0), "", 0, b});
    }
    put("float64_n2", "float64", gen_f64(8, 2, 0.0, 1.0));
    put("float64_n3", "float64", gen_f64(9, 3, 0.0, 1.0));
    put("float64_single", "float64", {2.5});
    put("float64_all_same", "float64", std::vector<double>(50, 0.7));
    f64("f_float64_uniform_70001", "float64", "uniform24", 15, 70001);
    // The generic CDF walk above 2^24 values: counts near p*n (>= 1.8e7)
    // are not float32 values, so a float32 count or float32 c0/c1 changes
    // the result.
    f64("f_float64_uniform_20000000", "float64", "uniform24", 16, 20000000);
    f64("f_uint16_adu_20000000", "uint16", "astro_adu", 17, 20000000);
    // Element types whose float64 conversion has its own failure mode: the
    // sign of int8, uint32 above 2^31, uint64 above 2^63, int64 beyond
    // 2^53, and a named type whose underlying type is float32 (it takes
    // the float64 path, not the []float32 Siril path).
    f64("int8_3000", "int8", "i8", 18, 3000);
    f64("uint32_3000", "uint32", "u32", 19, 3000);
    f64("uint64_3000", "uint64", "u64", 20, 3000);
    f64("int64_wide_3000", "int64", "i64", 21, 3000);
    f64("namedf32_5000", "namedf32", "uniform24", 22, 5000);

    char path[1024];
    snprintf(path, sizeof(path), "%s/siril_stats_golden.json", outdir);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return 1; }
    // One case per line: compact, and a diff stays per case.
    fprintf(f, "{\"generator\":\"stats/testdata/sref (make goldens)\",\n");
    fprintf(f, "\"percentiles\":[");
    for (int i = 0; i < NPS; i++) fprintf(f, "%s\"%s\"", i ? "," : "", PNAMES[i]);
    fprintf(f, "],\n\"f32\":[\n");
    for (size_t i = 0; i < cs.size(); i++) write_f32_case(f, cs[i], i + 1 == cs.size());
    fprintf(f, "],\n\"f64\":[\n");
    for (size_t i = 0; i < gs.size(); i++) write_f64_case(f, gs[i], i + 1 == gs.size());
    fprintf(f, "]}\n");
    fclose(f);
    printf("sref: %zu float32 cases, %zu float64 cases -> %s\n", cs.size(), gs.size(), path);
    return 0;
}

// median-f32 / median-f64: read whitespace-separated hex bit patterns on
// stdin, print the median's bit pattern. Exit 2 on empty input.
static int median_mode(bool f64) {
    std::vector<float> d32;
    std::vector<double> d64;
    char buf[64];
    while (scanf("%63s", buf) == 1) {
        unsigned long long u = strtoull(buf, NULL, 16);
        if (f64) d64.push_back(bitsf64(u)); else d32.push_back(bitsf32((uint32_t)u));
    }
    if (f64) {
        if (d64.empty()) return 2;
        printf("%016" PRIx64 "\n", f64bits(median_f64(d64)));
    } else {
        if (d32.empty()) return 2;
        double med = histogram_median_float(d32.data(), d32.size(), SINGLE_THREADED);
        printf("%08" PRIx32 "\n", f32bits((float)med));
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 3 && strcmp(argv[1], "gen") == 0) return gen(argv[2]);
    if (argc == 2 && strcmp(argv[1], "median-f32") == 0) return median_mode(false);
    if (argc == 2 && strcmp(argv[1], "median-f64") == 0) return median_mode(true);
    fprintf(stderr, "usage: %s gen <output_dir> | median-f32 | median-f64\n", argv[0]);
    return 2;
}
