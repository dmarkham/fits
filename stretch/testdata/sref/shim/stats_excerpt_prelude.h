/* Includes needed by ../siril/stats_excerpt.c (verbatim Siril functions). */
#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
void findMinMaxPercentile(const float* data, size_t size, float minPrct, float* minOut, float maxPrct, float* maxOut, int threads);
