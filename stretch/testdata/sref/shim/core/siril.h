/* Minimal stand-in for Siril's core/siril.h so the verbatim Siril 1.4.0
 * sources in ../siril compile outside the Siril tree. HAND-WRITTEN shim:
 * every macro/enum/constant that affects arithmetic is copied from Siril
 * 1.4.0 core/siril.h (line numbers given); the fits/imstats structs keep only
 * the members the extracted code touches, with Siril's member types. */
#ifndef SREF_SHIM_SIRIL_H
#define SREF_SHIM_SIRIL_H

#include <glib.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#define _(String) (String)
#define N_(String) (String)

/* core/siril.h:44-63 */
#ifndef __cplusplus
#undef max
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#undef min
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define SQR(x) ((x)*(x))
#endif

/* core/siril.h:70-86 */
#define USHRT_MAX_DOUBLE ((double)USHRT_MAX)
#define USHRT_MAX_SINGLE ((float)USHRT_MAX)
#define UCHAR_MAX_DOUBLE ((double)UCHAR_MAX)
typedef unsigned char BYTE;
typedef unsigned short WORD;

/* core/siril.h:104 */
#define MAD_NORM 1.4826

/* core/siril.h:148-150 */
#define RLAYER 0
#define GLAYER 1
#define BLAYER 2

/* cfitsio fitsio.h */
#define BYTE_IMG 8
#define USHORT_IMG 20
#define FLOAT_IMG -32

/* core/siril.h:245-250 */
typedef enum {
	CLIP,
	RESCALE,
	RGBBLEND,
	RESCALEGLOBAL
} clip_mode_t;

/* core/siril.h:284-288 */
typedef enum {
	MULTI_THREADED = 0,
	SINGLE_THREADED = 1
} threading_type;

/* core/siril.h:463 */
typedef enum { DATA_USHORT, DATA_FLOAT, DATA_UNSUPPORTED } data_type;

typedef struct {
	double median, mad, normValue;
	long total, ngoodpix;
} imstats;

typedef struct {
	int bitpix;
	data_type type;
	unsigned int rx, ry;
	long naxes[3];
	WORD *data, *pdata[3];
	float *fdata, *fpdata[3];
	double maxi;
} fits;

struct cominfo { int max_thread; };
extern struct cominfo com;

#define PRINT_ALLOC_ERR fprintf(stderr, "Out of memory in %s (%s:%d)\n", __func__, __FILE__, __LINE__)

#ifdef __cplusplus
}
#endif
#endif
