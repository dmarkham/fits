/*
 * sref.cc — Siril/RawTherapee reference harness for validating the Go
 * stats package against the real findMinMaxPercentile implementation.
 *
 * Build:   make
 * Use:     ./sref gen <output_dir>
 *
 * Compiles the actual RawTherapee findMinMaxPercentile from Siril's
 * source tree and calls it on deterministic test inputs, dumping the
 * results as JSON. The Go test diffs against these goldens.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <random>

// Minimal stubs so we can include rt_algo.cc without pulling all of Siril.
namespace rtengine {
    template<typename T>
    constexpr const T& LIM(const T& val, const T& low, const T& high) {
        return val < low ? low : (val > high ? high : val);
    }
}

// Include the actual findMinMaxPercentile implementation.
// We skip the gauss.h include (not needed for percentile) and the
// buildBlendMask function (not needed).
#define GAUSS_H_  // skip gauss.h include
#include "rt_algo_excerpt.cc"

// ----------------------------------------------------------------
// Test data generators (must match gen_golden.py exactly)
// ----------------------------------------------------------------

static std::vector<float> gen_float32_random(int seed, int n) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    std::vector<float> out(n);
    for (int i = 0; i < n; i++) out[i] = dist(rng);
    return out;
}

static std::vector<float> gen_gradient(int n, float lo, float hi) {
    std::vector<float> out(n);
    for (int i = 0; i < n; i++)
        out[i] = lo + (hi - lo) * float(i) / float(n - 1);
    return out;
}

// ----------------------------------------------------------------
// Siril-compatible stats using findMinMaxPercentile
// ----------------------------------------------------------------

static float siril_median(float* data, size_t n) {
    float med;
    findMinMaxPercentile(data, n, 0.5f, &med, 0.5f, &med, 1);
    return med;
}

static float siril_percentile(float* data, size_t n, float p) {
    float val;
    findMinMaxPercentile(data, n, p, &val, p, &val, 1);
    return val;
}

static float siril_mad(float* data, size_t n, float median) {
    std::vector<float> absDev(n);
    for (size_t i = 0; i < n; i++)
        absDev[i] = std::fabs(data[i] - median);
    float mad;
    findMinMaxPercentile(absDev.data(), n, 0.5f, &mad, 0.5f, &mad, 1);
    return mad;
}

// Sigma clipping matching Siril's iterative approach.
struct SigmaClipResult {
    double mean, stdev, median;
    int ngood;
};

static SigmaClipResult siril_sigma_clip(float* data, size_t n,
    double klow, double khigh, int maxiter, bool use_median) {
    // Copy data.
    std::vector<float> work(data, data + n);

    for (int iter = 0; iter < maxiter; iter++) {
        size_t m = work.size();
        if (m == 0) break;

        double center;
        if (use_median) {
            center = siril_median(work.data(), m);
        } else {
            double sum = 0;
            for (size_t i = 0; i < m; i++) sum += work[i];
            center = sum / m;
        }

        double sum2 = 0;
        for (size_t i = 0; i < m; i++) {
            double d = work[i] - center;
            sum2 += d * d;
        }
        double sigma = std::sqrt(sum2 / m);
        if (sigma == 0) break;

        double lo = center - klow * sigma;
        double hi = center + khigh * sigma;
        std::vector<float> kept;
        for (size_t i = 0; i < m; i++) {
            if (work[i] >= lo && work[i] <= hi)
                kept.push_back(work[i]);
        }
        if (kept.size() == work.size()) break;
        work = kept;
    }

    SigmaClipResult r = {0, 0, 0, (int)work.size()};
    if (work.empty()) return r;

    double sum = 0;
    for (auto v : work) sum += v;
    r.mean = sum / work.size();

    double sum2 = 0;
    for (auto v : work) { double d = v - r.mean; sum2 += d * d; }
    r.stdev = std::sqrt(sum2 / work.size());

    r.median = siril_median(work.data(), work.size());
    return r;
}

// ----------------------------------------------------------------
// Test cases
// ----------------------------------------------------------------

struct TestCase {
    const char* name;
    std::vector<float> data;
};

static void write_json(const char* outdir, const std::vector<TestCase>& cases) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/siril_stats_golden.json", outdir);
    FILE* f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }

    fprintf(f, "{\n");
    for (size_t ci = 0; ci < cases.size(); ci++) {
        const auto& c = cases[ci];
        std::vector<float> data = c.data; // copy for mutation
        size_t n = data.size();

        fprintf(f, "  \"%s\": {\n", c.name);
        fprintf(f, "    \"len\": %zu,\n", n);

        if (n == 0) {
            fprintf(f, "    \"median\": null,\n");
            fprintf(f, "    \"mad\": null,\n");
            fprintf(f, "    \"percentiles\": {},\n");
            fprintf(f, "    \"sigma_clip\": {}\n");
        } else {
            float med = siril_median(data.data(), n);
            float mad = siril_mad(data.data(), n, med);

            fprintf(f, "    \"median\": %.17g,\n", (double)med);
            fprintf(f, "    \"mad\": %.17g,\n", (double)mad);

            // Percentiles.
            fprintf(f, "    \"percentiles\": {\n");
            float pvals[] = {0.01f, 0.1f, 0.25f, 0.5f, 0.75f, 0.9f, 0.99f};
            for (int pi = 0; pi < 7; pi++) {
                float pv = siril_percentile(data.data(), n, pvals[pi]);
                fprintf(f, "      \"%.2g\": %.17g%s\n", pvals[pi], (double)pv,
                    pi < 6 ? "," : "");
            }
            fprintf(f, "    },\n");

            // Sigma clipping.
            fprintf(f, "    \"sigma_clip\": {\n");
            double ks[] = {2.0, 3.0, 5.0};
            const char* centers[] = {"mean", "median"};
            int first = 1;
            for (int ki = 0; ki < 3; ki++) {
                for (int ci2 = 0; ci2 < 2; ci2++) {
                    if (!first) fprintf(f, ",\n");
                    first = 0;
                    auto r = siril_sigma_clip(data.data(), n, ks[ki], ks[ki],
                        10, ci2 == 1);
                    fprintf(f, "      \"k%.1f_%s\": {\"mean\": %.17g, \"stdev\": %.17g, "
                        "\"median\": %.17g, \"ngood\": %d}",
                        ks[ki], centers[ci2], r.mean, r.stdev, r.median, r.ngood);
                }
            }
            fprintf(f, "\n    }\n");
        }

        fprintf(f, "  }%s\n", ci + 1 < cases.size() ? "," : "");
    }
    fprintf(f, "}\n");
    fclose(f);
    printf("Written to %s\n", path);
}

int main(int argc, char** argv) {
    if (argc < 3 || strcmp(argv[1], "gen") != 0) {
        fprintf(stderr, "usage: %s gen <output_dir>\n", argv[0]);
        return 2;
    }

    // Build test cases. Use deterministic data that we can reproduce
    // identically in the Go test.
    std::vector<TestCase> cases;

    // Gradient 0.01 to 0.30, 200 points.
    cases.push_back({"gradient_200", gen_gradient(200, 0.01f, 0.30f)});

    // Gradient 1000 points (enough bins for tight convergence).
    cases.push_back({"gradient_1000", gen_gradient(1000, 0.0f, 1.0f)});

    // Uniform random, 500 points.
    cases.push_back({"random_500", gen_float32_random(42, 500)});

    // Large random, 10000 points.
    cases.push_back({"random_10000", gen_float32_random(99, 10000)});

    // With zeros (simulating masked astro image).
    {
        auto d = gen_gradient(300, 0.01f, 0.25f);
        for (int i = 0; i < 60; i++) d[i] = 0.0f; // 20% zeros
        cases.push_back({"with_zeros", d});
    }

    // All same value.
    cases.push_back({"all_same", std::vector<float>(100, 0.42f)});

    // Single element.
    cases.push_back({"single", {3.14f}});

    printf("sref: %zu test cases\n", cases.size());
    write_json(argv[2], cases);
    return 0;
}
