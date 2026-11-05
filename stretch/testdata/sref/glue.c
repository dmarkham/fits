/* HAND-WRITTEN glue between the verbatim Siril 1.4.0 code in siril/ and the
 * oracle driver. It replaces Siril's stats cache / threading plumbing; the
 * arithmetic it performs is limited to calling the verbatim functions in the
 * order Siril's statistics_internal_float() (algos/statistics_float.c:295)
 * calls them for STATS_BASIC|STATS_MAD on a float image:
 *   ngoodpix  <- FnMeanSigma_float (non-zero, non-NaN count)       :403
 *   data      <- reassign_to_non_null_data_float if total != ngood  :427
 *   median    <- histogram_median_float(data, ngoodpix)             :442
 *   mad       <- siril_stats_float_mad(data, ngoodpix, median)      :462
 * and returns an error (Siril: stats NULL) when ngoodpix == 0 (:418). */
#include "core/siril.h"
#include "core/proto.h"
#include "algos/statistics.h"
#include "algos/colors.h"

struct cominfo com = { 1 };

/* core/utils.c:825, float branch (the oracle never builds ushort fits). */
double get_normalized_value(fits *fit) {
	if (fit->type == DATA_USHORT)
		return USHRT_MAX_DOUBLE;
	return 1.0;
}

void invalidate_stats_from_fit(fits *fit) { (void)fit; }

/* core/processing.c:1417, verbatim. */
int check_threading(threading_type *threads) {
	if (*threads == MULTI_THREADED)
		*threads = com.max_thread;
	return *threads;
}

void rgb_to_hslf(float r, float g, float b, float *h, float *s, float *l) { abort(); }
void hsl_to_rgbf(float h, float s, float l, float *r, float *g, float *b) { abort(); }
void rgbw_to_hslw(WORD r, WORD g, WORD b, WORD *h, WORD *s, WORD *l) { abort(); }
void hslw_to_rgbw(WORD h, WORD s, WORD l, WORD *r, WORD *g, WORD *b) { abort(); }

void free_stats(imstats *stat) { free(stat); }

/* The verbatim statistics functions (two of them are static in Siril, so
 * the excerpt is compiled as part of this translation unit). */
#include "siril/stats_excerpt.c"

int compute_all_channels_statistics_single_image(fits *fit, int option,
		threading_type threading, imstats **stats) {
	(void)option;
	int retval = 0;
	for (int c = 0; c < (int)fit->naxes[2]; c++) {
		imstats *st = calloc(1, sizeof(imstats));
		stats[c] = st;
		float *data = fit->fpdata[c];
		st->total = (long)fit->rx * fit->ry;
		st->normValue = 1.0; /* FLOAT_IMG */
		double mean, sigma;
		int status = 0;
		FnMeanSigma_float(data, st->total, &st->ngoodpix, &mean, &sigma, &status);
		if (st->ngoodpix == 0) {
			free(st);
			stats[c] = NULL;
			retval = -1;
			continue;
		}
		int free_data = 0;
		if (st->total != st->ngoodpix) {
			data = reassign_to_non_null_data_float(data, st->total, st->ngoodpix, 0);
			free_data = 1;
		}
		st->median = histogram_median_float(data, st->ngoodpix, threading) * st->normValue;
		st->mad = siril_stats_float_mad(data, st->ngoodpix, st->median / st->normValue, threading, NULL) * st->normValue;
		if (free_data)
			free(data);
	}
	return retval;
}
