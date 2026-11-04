/* HAND-WRITTEN, used only to build sref_libm (see Makefile): the verbatim
 * filters/ght.c compiled with logf/powf/expf replaced by the double
 * function rounded to float, which is what the Go port computes. It is NOT
 * a truth source: gen_siril_goldens.py runs it next to the real oracle only
 * to list the pixels where glibc's float logf/powf/expf (not correctly
 * rounded) change Siril's GHT output, so the golden can store those pixels
 * explicitly (the "libm" record of a GHT golden). Expected values always
 * come from the verbatim build. */
#ifndef SREF_LIBM_DOUBLE_H
#define SREF_LIBM_DOUBLE_H
#include <math.h>
static inline float sref_logf(float x) { return (float)log((double)x); }
static inline float sref_powf(float x, float y) { return (float)pow((double)x, (double)y); }
static inline float sref_expf(float x) { return (float)exp((double)x); }
#define logf sref_logf
#define powf sref_powf
#define expf sref_expf
#endif
