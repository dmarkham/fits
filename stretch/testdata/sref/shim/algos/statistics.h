/* Stand-in for Siril algos/statistics.h. compute_all_channels_statistics_
 * single_image is implemented in ../glue.c on top of the verbatim Siril
 * statistics functions in ../siril/stats_excerpt.c. */
#ifndef SREF_SHIM_STATISTICS_H
#define SREF_SHIM_STATISTICS_H
#include "core/siril.h"
#ifdef __cplusplus
extern "C" {
#endif
/* algos/statistics.h STATS_* flags (values are irrelevant to the oracle;
 * glue.c always computes the median and MAD it is asked for). */
#define STATS_BASIC (1 << 1)
#define STATS_MAD (1 << 3)
int compute_all_channels_statistics_single_image(fits *fit, int option,
		threading_type threading, imstats **stats);
void free_stats(imstats *stat);
double histogram_median_float(float *a, size_t n, threading_type threading);
double siril_stats_float_mad(const float *data, const size_t n, const double m, threading_type threads, float *buffer);
#ifdef __cplusplus
}
#endif
#endif
