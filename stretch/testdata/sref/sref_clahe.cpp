/* sref_clahe: runs Siril 1.4.0's cvClahe_float (opencv/opencv.cpp:806,
 * VERBATIM in siril/opencv_excerpt.cpp) against the system OpenCV, exactly
 * as the `clahe cliplimit tilesize` command does on a 32-bit float image
 * (command.c:2698 process_clahe -> filters/clahe.c:94 cvClahe).
 *
 *   sref_clahe NCH NX NY IN OUT cliplimit tilesize
 *
 * IN/OUT are raw little-endian float32, channel-major (Siril fdata layout).
 * Also: sref_clahe --version prints the OpenCV version linked. */
#include "siril/opencv_excerpt.cpp"

extern "C" void invalidate_stats_from_fit(fits *fit) { (void)fit; }

int main(int argc, char **argv) {
	if (argc == 2 && !strcmp(argv[1], "--version")) {
		printf("%s\n", CV_VERSION);
		return 0;
	}
	if (argc != 8) {
		fprintf(stderr, "usage: sref_clahe NCH NX NY IN OUT cliplimit tilesize\n");
		return 2;
	}
	int nch = atoi(argv[1]), nx = atoi(argv[2]), ny = atoi(argv[3]);
	size_t layer = (size_t)nx * ny, n = layer * nch;
	fits fit = {};
	fit.bitpix = FLOAT_IMG;
	fit.type = DATA_FLOAT;
	fit.rx = nx;
	fit.ry = ny;
	fit.naxes[0] = nx;
	fit.naxes[1] = ny;
	fit.naxes[2] = nch;
	fit.fdata = (float *)malloc(n * sizeof(float));
	for (int c = 0; c < 3; c++)
		fit.fpdata[c] = fit.fdata + (nch == 3 ? c * layer : 0);
	FILE *f = fopen(argv[4], "rb");
	if (!f || fread(fit.fdata, sizeof(float), n, f) != n) {
		fprintf(stderr, "cannot read %s\n", argv[4]);
		return 1;
	}
	fclose(f);
	/* command.c:2700-2712 */
	double clip_limit = g_ascii_strtod(argv[6], NULL);
	int size = g_ascii_strtoull(argv[7], NULL, 10);
	if (cvClahe_float(&fit, clip_limit, size))
		return 1;
	f = fopen(argv[5], "wb");
	if (!f || fwrite(fit.fdata, sizeof(float), n, f) != n) {
		fprintf(stderr, "cannot write %s\n", argv[5]);
		return 1;
	}
	fclose(f);
	return 0;
}
