/*
 * siril_stubs.h -- the minimum of Siril's headers that siril_excerpt.c
 * needs, so the verbatim functions compile outside the Siril tree.
 * Nothing here takes part in the arithmetic: threading is forced to a
 * single thread (findMinMaxPercentile's per-thread histograms sum to
 * the same integer counts anyway) and allocation errors abort.
 */
#ifndef SREF_SIRIL_STUBS_H
#define SREF_SIRIL_STUBS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/* src/core/processing.h */
typedef int threading_type;
#define SINGLE_THREADED 1
static inline int check_threading(threading_type *threads) { (void)threads; return 1; }

/* src/core/siril.h */
#define PRINT_ALLOC_ERR fprintf(stderr, "Out of memory in %s (%s:%d) - aborting\n", __func__, __FILE__, __LINE__)

/* src/rt/rt_algo.h (C linkage, defined in rt_algo_excerpt.cc) */
void findMinMaxPercentile(const float* data, size_t size, float minPrct, float* minOut, float maxPrct, float* maxOut, int threads);

/* Exported wrappers for the file-static Siril helpers (siril_excerpt.c). */
float siril_stats_float_sd(const float data[], const int N, float *m);
double siril_stats_float_mad(const float *data, const size_t n, const double m, threading_type threads, float *buffer);
double sortnet_median_float(float *a, size_t n);
double quickmedian_float(float *a, size_t n);
double histogram_median_float(float *a, size_t n, threading_type threading);
float *sref_reassign_to_non_null_data_float(float *data, size_t inputlen, size_t outputlen);

#ifdef __cplusplus
}
#endif

#endif
