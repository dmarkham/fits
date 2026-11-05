/* Stand-in for Siril core/proto.h: roundf_to_WORD is copied verbatim from
 * Siril 1.4.0 core/proto.h:336; the rest are declarations implemented in
 * ../glue.c. */
#ifndef SREF_SHIM_PROTO_H
#define SREF_SHIM_PROTO_H
#include "core/siril.h"
#ifdef __cplusplus
extern "C" {
#endif
static inline WORD roundf_to_WORD(float f) {
	f = f + 0.5f;
	f = (f > (float)USHRT_MAX) ? (float)USHRT_MAX : f;
	f = (f < 0.0f) ? 0.0f : f;
	return (WORD)f;
}
double get_normalized_value(fits *fit);
void invalidate_stats_from_fit(fits *fit);
int check_threading(threading_type *threads);
#ifdef __cplusplus
}
#endif
#endif
