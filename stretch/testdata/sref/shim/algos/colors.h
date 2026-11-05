/* Stand-in for Siril algos/colors.h: only the saturation GHT (not used by
 * the oracle) calls these; glue.c defines them to abort. */
#ifndef SREF_SHIM_COLORS_H
#define SREF_SHIM_COLORS_H
#include "core/siril.h"
void rgb_to_hslf(float r, float g, float b, float *h, float *s, float *l);
void hsl_to_rgbf(float h, float s, float l, float *r, float *g, float *b);
void rgbw_to_hslw(WORD r, WORD g, WORD b, WORD *h, WORD *s, WORD *l);
void hslw_to_rgbw(WORD h, WORD s, WORD l, WORD *r, WORD *g, WORD *b);
#endif
