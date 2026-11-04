/*
 * siril_wrappers.c -- compiles the verbatim Siril functions and exports
 * the one that Siril declares static. Ours, not Siril's.
 */
#include "siril_excerpt.c"

/* reassign_to_non_null_data_float is static in statistics_float.c. Never
 * frees the input (free_input = 0). */
float *sref_reassign_to_non_null_data_float(float *data, size_t inputlen, size_t outputlen) {
	return reassign_to_non_null_data_float(data, inputlen, outputlen, 0);
}
