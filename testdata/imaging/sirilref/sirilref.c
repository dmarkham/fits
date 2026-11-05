/*
 * sirilref -- Siril's FITS-to-float conversion, compiled from Siril's own
 * source, as the truth for fits.ReadFloat32's "Siril convention" claims.
 *
 *   sirilref <file.fits>   ->  JSON on stdout:
 *     {"file":..., "bitpix": <after manage_bitpix>, "orig_bitpix":...,
 *      "data_max":..., "normalized": true|false, "status": <cfitsio>,
 *      "unsupported": true|false, "bits": [<float32 bit patterns>]}
 *   "bits" holds more than DIGEST_MIN (1024) values as an exact digest
 *   instead: {"n": N, "sha256": <hex of sha256 over the big-endian int64
 *   of every bit pattern>, "sample": [[index, bits], ...] at the first and
 *   last 4 indices, 32 evenly spaced ones and the first minimum and
 *   maximum} (decoded by imaging_codec_test.go).
 *
 * This reproduces readfits(..., force_float=TRUE) — Siril's path that
 * yields 32-bit float pixels in [0, 1] — for the first image HDU:
 *   read_fits_metadata: fits_get_img_param + manage_bitpix, LONGLONG_IMG
 *     rejected;
 *   read_fits_header: keywords.data_max = DATAMAX keyword (0 if absent),
 *     overridden by fit_stats' maximum for FLOAT_IMG files whose PROGRAM
 *     is not Siril and for DOUBLE_IMG (fits_keywords.c:1135-1142);
 *   read_fits_with_convert, FLOAT_IMG case: fits_read_img(TFLOAT, nulval
 *     "zero" = no null checking) then convert_floats() when the original
 *     bitpix is BYTE/SHORT/USHORT or data_max > 10.
 *
 * The functions marked VERBATIM below are copied unchanged from Siril
 * src/io/image_format_fits.c (git 55111156048b, 1.4.0-rc1-235, 2026-01-10;
 * identical in the siril-1.4.0 release tarball, checked by `make verbatim`
 * with SIRIL_SRC set). The constants are Siril's src/core/siril.h:70-77.
 * Siril is GPL-3.0-or-later; this file is a test harness only and is not
 * linked into the Go library.
 */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

/* ---- SHA-256 (FIPS 180-4), the self-contained implementation of
 * testdata/cfref/cfref.c, so the harness links only against libcfitsio. */
typedef struct { unsigned int h[8]; unsigned char buf[64]; unsigned long long len; int nbuf; } sha256_ctx;
static const unsigned int sha_k[64] = {
  0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
  0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
  0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
  0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
  0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
  0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
  0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
  0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
#define ROR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))
static void sha_block(sha256_ctx *c, const unsigned char *p) {
  unsigned int w[64], a, b, cc, d, e, f, g, h;
  for (int i = 0; i < 16; i++)
    w[i] = (unsigned int)p[4 * i] << 24 | (unsigned int)p[4 * i + 1] << 16 | (unsigned int)p[4 * i + 2] << 8 | p[4 * i + 3];
  for (int i = 16; i < 64; i++) {
    unsigned int s0 = ROR(w[i - 15], 7) ^ ROR(w[i - 15], 18) ^ (w[i - 15] >> 3);
    unsigned int s1 = ROR(w[i - 2], 17) ^ ROR(w[i - 2], 19) ^ (w[i - 2] >> 10);
    w[i] = w[i - 16] + s0 + w[i - 7] + s1;
  }
  a = c->h[0]; b = c->h[1]; cc = c->h[2]; d = c->h[3]; e = c->h[4]; f = c->h[5]; g = c->h[6]; h = c->h[7];
  for (int i = 0; i < 64; i++) {
    unsigned int t1 = h + (ROR(e, 6) ^ ROR(e, 11) ^ ROR(e, 25)) + ((e & f) ^ (~e & g)) + sha_k[i] + w[i];
    unsigned int t2 = (ROR(a, 2) ^ ROR(a, 13) ^ ROR(a, 22)) + ((a & b) ^ (a & cc) ^ (b & cc));
    h = g; g = f; f = e; e = d + t1; d = cc; cc = b; b = a; a = t1 + t2;
  }
  c->h[0] += a; c->h[1] += b; c->h[2] += cc; c->h[3] += d; c->h[4] += e; c->h[5] += f; c->h[6] += g; c->h[7] += h;
}
static void sha_init(sha256_ctx *c) {
  static const unsigned int iv[8] = {0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
  memcpy(c->h, iv, sizeof iv); c->len = 0; c->nbuf = 0;
}
static void sha_update(sha256_ctx *c, const unsigned char *p, size_t n) {
  c->len += n;
  while (n--) { c->buf[c->nbuf++] = *p++; if (c->nbuf == 64) { sha_block(c, c->buf); c->nbuf = 0; } }
}
static void sha_hex(sha256_ctx *c, char out[65]) {
  unsigned long long bits = c->len * 8;
  unsigned char pad = 0x80, z = 0, l[8];
  sha_update(c, &pad, 1);
  while (c->nbuf != 56) sha_update(c, &z, 1);
  for (int i = 0; i < 8; i++) l[i] = (unsigned char)(bits >> (56 - 8 * i));
  sha_update(c, l, 8);
  for (int i = 0; i < 8; i++) sprintf(out + 8 * i, "%08x", c->h[i]);
}
static void sha_u64(sha256_ctx *c, unsigned long long v) {
  unsigned char b[8];
  for (int i = 0; i < 8; i++) b[i] = (unsigned char)(v >> (56 - 8 * i));
  sha_update(c, b, 8);
}

#define DIGEST_MIN 1024

/* print_bits writes the "bits" value: an explicit list, or the digest. */
static void print_bits(const float *fdata, size_t n) {
	if (n <= DIGEST_MIN) {
		printf("[");
		for (size_t i = 0; i < n; i++) {
			uint32_t b;
			memcpy(&b, &fdata[i], 4);
			printf(i ? ",%u" : "%u", b);
		}
		printf("]");
		return;
	}
	sha256_ctx c;
	sha_init(&c);
	size_t imin = 0, imax = 0;
	uint32_t vmin = UINT32_MAX, vmax = 0;
	for (size_t i = 0; i < n; i++) {
		uint32_t b;
		memcpy(&b, &fdata[i], 4);
		sha_u64(&c, (unsigned long long)b);
		if (b < vmin) { vmin = b; imin = i; }
		if (b > vmax) { vmax = b; imax = i; }
	}
	char hex[65];
	sha_hex(&c, hex);
	printf("{\"n\":%zu,\"sha256\":\"%s\",\"sample\":[", n, hex);
	int first = 1;
	for (size_t i = 0; i < n; i++) {
		int take = i < 4 || i >= n - 4 || i == imin || i == imax;
		for (size_t k = 0; k < 32 && !take; k++)
			take = (k * (n - 1) / 31 == i);
		if (!take) continue;
		uint32_t b;
		memcpy(&b, &fdata[i], 4);
		printf(first ? "[%zu,%u]" : ",[%zu,%u]", i, b);
		first = 0;
	}
	printf("]}");
}

typedef int gboolean;
#define TRUE 1
#define FALSE 0
#define _(s) (s)
#define siril_log_message(...) ((void)0)
#define siril_debug_print(...) ((void)0)
#define report_fits_error(s) ((void)0)
#define PRINT_ALLOC_ERR ((void)0)

/* Siril src/core/siril.h:70-77 */
#define USHRT_MAX_DOUBLE ((double)USHRT_MAX)
#define USHRT_MAX_SINGLE ((float)USHRT_MAX)
#define UCHAR_MAX_SINGLE ((float)UCHAR_MAX)
#define INV_USHRT_MAX_SINGLE .000015259022f	 // 1/65535
#define INV_UCHAR_MAX_SINGLE .0039215686f	 // 1/255

/* VERBATIM: image_format_fits.c fit_stats */
int fit_stats(fitsfile *fptr, float *mini, float *maxi) {
	int status = 0;
	int ii;
	long npixels = 1L;
	long anaxes[3] = { 1L, 1L, 1L }, firstpix[3] = { 1L, 1L, 1L };
	float *pix;
	float minval = 1.E33, maxval = -1.E33;

	/* initialize value in case where it does not work */
	*mini = 0;
	*maxi = 0;

	fits_get_img_size(fptr, 3, anaxes, &status);

	if (status) {
		report_fits_error(status); /* print error message */
		return(status);
	}

	npixels = anaxes[0];  /* no. of pixels to read in each row */
	pix = malloc(npixels * sizeof(float)); /* memory for 1 row */
	if (pix == NULL) {
		PRINT_ALLOC_ERR;
		return (1);
	}

	/* loop over all planes of the cube (2D images have 1 plane) */
	for (firstpix[2] = 1; firstpix[2] <= anaxes[2]; firstpix[2]++) {
		/* loop over all rows of the plane */
		for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
			/* give starting pixel coordinate and number of pixels to read */
			if (fits_read_pix(fptr, TFLOAT, firstpix, npixels, NULL, pix,
						NULL, &status))
				break; /* jump out of loop on error */

			for (ii = 0; ii < npixels; ii++) {
				if (pix[ii] < minval)
					minval = pix[ii]; /* find min and  */
				if (pix[ii] > maxval)
					maxval = pix[ii]; /* max values    */
			}
		}
	}    /* end of loop over planes */
	free(pix);

	if (status) {
		report_fits_error(status); /* print any error message */
	} else {
		siril_debug_print("  minimum value = %f\n", minval);
		siril_debug_print("  maximum value = %f\n", maxval);
		*maxi = maxval;
		*mini = minval;
	}
	return status;
}

/* VERBATIM: image_format_fits.c convert_floats */
static void convert_floats(int bitpix, float *data, size_t nbdata, gboolean verbose) {
	size_t i;
	switch (bitpix) {
		case BYTE_IMG:
			for (i = 0; i < nbdata; i++)
				data[i] = data[i] * INV_UCHAR_MAX_SINGLE;
			break;
		case FLOAT_IMG:
			if (verbose)
				siril_log_message(_("Normalizing input data to our float range [0, 1]\n")); // @suppress("No break at end of case")
				// Fallthrough is deliberate, this handles FITS with floating point data
				// where MAX is 65565 (e.g. JWST)
		default:
		case USHORT_IMG:	// siril 0.9 native
			for (i = 0; i < nbdata; i++)
				data[i] = data[i] * INV_USHRT_MAX_SINGLE;
			break;
		case SHORT_IMG:
			// add 2^15 to the read data to obtain unsigned
			for (i = 0; i < nbdata; i++) {
				data[i] = (32768.f + data[i]) * INV_USHRT_MAX_SINGLE;
			}
			break;
	}
}

/* VERBATIM: image_format_fits.c manage_bitpix */
void manage_bitpix(fitsfile *fptr, int *bitpix, int *orig_bitpix) {
	double offset;
	int status = 0;
	fits_read_key(fptr, TDOUBLE, "BZERO", &offset, NULL, &status);
	if (!status) {
		if (*bitpix == SHORT_IMG && offset != 0.0) {
			*bitpix = USHORT_IMG;
		}
		else if (*bitpix == LONG_IMG && offset != 0.0)
			*bitpix = ULONG_IMG;
	} else {
		/* but some software just put unsigned 16-bit data in the file
		 * and don't set the BZERO keyword... */
		if (status == KEY_NO_EXIST && *bitpix == SHORT_IMG)
			*bitpix = USHORT_IMG;
	}
	// and we store the original bitpix to reuse it during later partial
	// reads when we have no access to the header or a fits * struct.
	*orig_bitpix = *bitpix;
}

int main(int argc, char **argv) {
	if (argc != 2) {
		fprintf(stderr, "usage: sirilref <file.fits>\n");
		return 2;
	}
	const char *base = strrchr(argv[1], '/');
	base = base ? base + 1 : argv[1];
	fitsfile *fptr;
	int status = 0, bitpix = 0, orig_bitpix = 0, naxis = 0;
	long naxes[3] = { 1, 1, 1 };
	/* siril_fits_move_first_image: first image HDU with NAXIS > 0. */
	fits_open_image(&fptr, argv[1], READONLY, &status);
	if (status) {
		printf("{\"file\":\"%s\",\"status\":%d}\n", base, status);
		return 0;
	}
	/* read_fits_metadata */
	fits_get_img_param(fptr, 3, &bitpix, &naxis, naxes, &status);
	manage_bitpix(fptr, &bitpix, &orig_bitpix);
	if (naxis == 2) naxes[2] = 1;
	if (bitpix == LONGLONG_IMG) {
		printf("{\"file\":\"%s\",\"bitpix\":%d,\"unsupported\":true,\"status\":0}\n", base, bitpix);
		fits_close_file(fptr, &status);
		return 0;
	}
	/* read_fits_header: keywords.data_max */
	double data_max = 0.0;
	int st2 = 0;
	fits_read_key(fptr, TDOUBLE, "DATAMAX", &data_max, NULL, &st2);
	if (st2) data_max = 0.0;
	char program[FLEN_VALUE] = "";
	st2 = 0;
	fits_read_key(fptr, TSTRING, "PROGRAM", program, NULL, &st2);
	gboolean not_from_siril = (strstr(program, "Siril") == NULL);
	if ((bitpix == FLOAT_IMG && not_from_siril) || bitpix == DOUBLE_IMG) {
		float mini, maxi;
		fit_stats(fptr, &mini, &maxi);
		data_max = (double) maxi;
	}
	/* read_fits_with_convert(force_float=TRUE): FLOAT_IMG branch */
	size_t nbdata = (size_t)naxes[0] * naxes[1] * naxes[2];
	float *fdata = malloc(nbdata * sizeof(float));
	int zero = 0;
	fits_read_img(fptr, TFLOAT, 1, nbdata, &zero, fdata, &zero, &status);
	gboolean normalized = FALSE;
	if ((bitpix == USHORT_IMG || bitpix == SHORT_IMG
			|| bitpix == BYTE_IMG) || data_max > 10.0) {
		convert_floats(bitpix, fdata, nbdata, TRUE);
		normalized = TRUE;
	}
	char dmax[40];
	if (isnan(data_max)) snprintf(dmax, sizeof dmax, "\"NaN\"");
	else if (isinf(data_max)) snprintf(dmax, sizeof dmax, data_max > 0 ? "\"Inf\"" : "\"-Inf\"");
	else snprintf(dmax, sizeof dmax, "%.17g", data_max);
	printf("{\"file\":\"%s\",\"bitpix\":%d,\"orig_bitpix\":%d,\"data_max\":%s,\"normalized\":%s,\"status\":%d,\"unsupported\":false,\"bits\":",
			base, bitpix, orig_bitpix, dmax, normalized ? "true" : "false", status);
	print_bits(fdata, status ? 0 : nbdata);
	printf("}\n");
	free(fdata);
	st2 = 0;
	fits_close_file(fptr, &st2);
	return 0;
}
