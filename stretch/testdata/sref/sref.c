/* sref: Siril 1.4.0 stretch oracle driver.
 *
 * Runs the VERBATIM Siril 1.4.0 processing functions in siril/ on a raw
 * little-endian float32 image (channel-major, Siril's fdata layout) and
 * writes the result the same way. The per-command argument handling below is
 * HAND-WRITTEN and mirrors core/command.c of Siril 1.4.0 line by line where
 * it affects the numbers (cited per op); it skips Siril's range validation so
 * function-level cases the CLI refuses (e.g. GHT with B < 0, which the GUI
 * allows) can be generated -- each such case is labelled in the goldens.
 *
 *   sref asinh       NCH NX NY IN OUT beta offset human(0|1)
 *   sref autostretch NCH NX NY IN OUT linked(0|1) shadows_clip target_bg
 *   sref asparams    NCH NX NY IN OUT linked(0|1) shadows_clip target_bg
 *                    (OUT gets "retval s m h" lines, one per parameter set)
 *   sref mtf         NCH NX NY IN OUT shadows midtones highlights
 *   sref ght         NCH NX NY IN OUT D B LP SP HP
 *   sref ghtcoef     NCH NX NY IN OUT D B LP SP HP
 *                    (OUT gets the GHTsetup coefficients b1 a2 b2 c2 d2 e2
 *                    a3 b3 c3 d3 e3 a4 b4 as float32 hex bits; IN is read
 *                    but unused)
 *   sref linstretch  NCH NX NY IN OUT BP
 */
#include "core/siril.h"
#include "filters/mtf.h"
#include "filters/ght.h"
#include <string.h>

int command_asinh(fits *fit, float beta, float offset, gboolean human_luminance, clip_mode_t clipmode);

static double arg(const char *s) {
	gchar *end;
	double v = g_ascii_strtod(s, &end);
	if (end == s) {
		fprintf(stderr, "bad number %s\n", s);
		exit(2);
	}
	return v;
}

int main(int argc, char **argv) {
	if (argc < 7) {
		fprintf(stderr, "usage: sref OP NCH NX NY IN OUT args...\n");
		return 2;
	}
	const char *op = argv[1];
	int nch = atoi(argv[2]), nx = atoi(argv[3]), ny = atoi(argv[4]);
	const char *in = argv[5], *out = argv[6];
	char **a = argv + 7;
	int na = argc - 7;
	size_t layer = (size_t)nx * ny, n = layer * nch;

	fits fit = { 0 };
	fit.bitpix = FLOAT_IMG;
	fit.type = DATA_FLOAT;
	fit.rx = nx;
	fit.ry = ny;
	fit.naxes[0] = nx;
	fit.naxes[1] = ny;
	fit.naxes[2] = nch;
	fit.fdata = malloc(n * sizeof(float));
	for (int c = 0; c < 3; c++)
		fit.fpdata[c] = fit.fdata + (nch == 3 ? c * layer : 0);

	FILE *f = fopen(in, "rb");
	if (!f || fread(fit.fdata, sizeof(float), n, f) != n) {
		fprintf(stderr, "cannot read %zu floats from %s\n", n, in);
		return 1;
	}
	fclose(f);

	if (!strcmp(op, "asinh") && na == 3) {
		/* command.c:2643 process_asinh: clip mode defaults to RGBBLEND */
		command_asinh(&fit, arg(a[0]), arg(a[1]), atoi(a[2]) != 0, RGBBLEND);
	} else if ((!strcmp(op, "autostretch") || !strcmp(op, "asparams")) && na == 3) {
		/* command.c:3317 process_autostretch */
		gboolean linked = atoi(a[0]) != 0;
		float shadows_clipping = arg(a[1]);
		float target_bg = arg(a[2]);
		FILE *pf = NULL;
		if (!strcmp(op, "asparams"))
			pf = fopen(out, "w");
		if (linked) {
			struct mtf_params params;
			int r = find_linked_midtones_balance(&fit, shadows_clipping, target_bg, &params);
			params.do_red = TRUE;
			params.do_green = TRUE;
			params.do_blue = TRUE;
			if (pf)
				fprintf(pf, "%d %.9g %.9g %.9g\n", r, params.shadows, params.midtones, params.highlights);
			else
				apply_linked_mtf_to_fits(&fit, &fit, params, TRUE);
		} else {
			struct mtf_params params[3];
			int r = find_unlinked_midtones_balance(&fit, shadows_clipping, target_bg, params);
			if (pf) {
				for (int c = 0; c < nch; c++)
					fprintf(pf, "%d %.9g %.9g %.9g\n", r, params[c].shadows, params[c].midtones, params[c].highlights);
			} else
				apply_unlinked_mtf_to_fits(&fit, &fit, params);
		}
		if (pf) {
			fclose(pf);
			return 0;
		}
	} else if (!strcmp(op, "mtf") && na == 3) {
		struct mtf_params params = { .shadows = arg(a[0]), .midtones = arg(a[1]), .highlights = arg(a[2]),
			.do_red = TRUE, .do_green = TRUE, .do_blue = TRUE };
		apply_linked_mtf_to_fits(&fit, &fit, params, TRUE);
	} else if (!strcmp(op, "ght") && na == 5) {
		/* command.c:2235 process_ght_args, STRETCH_PAYNE_NORMAL, default
		 * colour model COL_INDEP and clip mode RGBBLEND */
		double D = arg(a[0]), B = arg(a[1]), LP = arg(a[2]), SP = arg(a[3]), HP = arg(a[4]);
		if (fabs(B) < 1.e-3f)
			B = 0.f;
		ght_params params = { 0 };
		params.B = (float) B;
		params.D = (float) expm1(D);
		params.LP = (float) LP;
		params.SP = (float) SP;
		params.HP = (float) HP;
		params.BP = (float) 0.0;
		params.stretchtype = STRETCH_PAYNE_NORMAL;
		params.payne_colourstretchmodel = COL_INDEP;
		params.do_red = params.do_green = params.do_blue = TRUE;
		params.clip_mode = RGBBLEND;
		apply_linked_ght_to_fits(&fit, &fit, &params, TRUE);
	} else if (!strcmp(op, "ghtcoef") && na == 5) {
		/* the parameter conversion of the "ght" op above, then the verbatim
		 * GHTsetup (filters/ght.c:34); the struct is zeroed first because
		 * GHTsetup leaves e2/e3 unset for B == -1 and B == 0 */
		double D = arg(a[0]), B = arg(a[1]), LP = arg(a[2]), SP = arg(a[3]), HP = arg(a[4]);
		if (fabs(B) < 1.e-3f)
			B = 0.f;
		ght_compute_params c;
		memset(&c, 0, sizeof c);
		GHTsetup(&c, (float) B, (float) expm1(D), (float) LP, (float) SP, (float) HP, STRETCH_PAYNE_NORMAL);
		const float v[] = { c.b1, c.a2, c.b2, c.c2, c.d2, c.e2, c.a3, c.b3, c.c3, c.d3, c.e3, c.a4, c.b4 };
		FILE *pf = fopen(out, "w");
		for (size_t k = 0; k < sizeof v / sizeof v[0]; k++) {
			uint32_t u;
			memcpy(&u, &v[k], 4);
			fprintf(pf, "%08x\n", u);
		}
		fclose(pf);
		return 0;
	} else if (!strcmp(op, "linstretch") && na == 1) {
		/* command.c:2547 process_linstretch -> process_ght_args with
		 * STRETCH_LINEAR: D = B = LP = SP = HP = 0 */
		ght_params params = { 0 };
		params.BP = (float) arg(a[0]);
		params.stretchtype = STRETCH_LINEAR;
		params.payne_colourstretchmodel = COL_INDEP;
		params.do_red = params.do_green = params.do_blue = TRUE;
		params.clip_mode = RGBBLEND;
		apply_linked_ght_to_fits(&fit, &fit, &params, TRUE);
	} else {
		fprintf(stderr, "unknown op or wrong argument count: %s\n", op);
		return 2;
	}

	f = fopen(out, "wb");
	if (!f || fwrite(fit.fdata, sizeof(float), n, f) != n) {
		fprintf(stderr, "cannot write %s\n", out);
		return 1;
	}
	fclose(f);
	return 0;
}
