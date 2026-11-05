/* cfnull: what cfitsio reads from each image of a (tile-compressed) FITS
 * file: the scaling it applies, and every pixel as a double with cfitsio's
 * own null flag (fits_read_pixnull, TDOUBLE), or the cfitsio error.
 *
 *   cfnull FILE   -> JSON on stdout
 *
 * JSON: {"file": base name, "cfitsio": version, "image_hdu": 1-based HDU,
 *        "status": cfitsio status of the pixel read (0 = ok),
 *        "error": cfitsio error text or "", "npix": n,
 *        "phys": [value as %.17g or null where cfitsio flags a null]}
 * The first HDU with NAXIS > 0 is used (HDU 2 for a tile-compressed file).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

int main(int argc, char **argv)
{
    fitsfile *f;
    int status = 0, nhdu = 0, naxis = 0, hdu = 0, anynul = 0;
    long naxes[9];
    float ver;
    if (argc != 2) { fprintf(stderr, "usage: cfnull FILE\n"); return 2; }
    if (fits_open_file(&f, argv[1], READONLY, &status)) { fits_report_error(stderr, status); return 1; }
    fits_get_num_hdus(f, &nhdu, &status);
    for (int i = 1; i <= nhdu && !status; i++) {
        int type;
        fits_movabs_hdu(f, i, &type, &status);
        if (type == IMAGE_HDU || fits_is_compressed_image(f, &status)) {
            fits_get_img_dim(f, &naxis, &status);
            if (naxis > 0) { hdu = i; break; }
        }
    }
    if (status || !hdu) { fits_report_error(stderr, status); return 1; }
    fits_get_img_size(f, naxis, naxes, &status);
    long long n = 1;
    for (int i = 0; i < naxis; i++) n *= naxes[i];
    double *v = calloc((size_t)n, sizeof *v);
    char *nul = calloc((size_t)n, 1);
    long fpix[9];
    for (int i = 0; i < 9; i++) fpix[i] = 1;
    int rs = 0;
    fits_read_pixnull(f, TDOUBLE, fpix, n, v, nul, &anynul, &rs);
    char err[FLEN_ERRMSG] = "";
    if (rs) fits_get_errstatus(rs, err);
    const char *base = strrchr(argv[1], '/');
    base = base ? base + 1 : argv[1];
    printf("{\"file\":\"%s\",\"cfitsio\":\"%.3f\",\"image_hdu\":%d,\"status\":%d,\"error\":\"%s\",\"npix\":%lld,\"phys\":",
           base, fits_get_version(&ver), hdu, rs, err, n);
    if (rs) {
        printf("null}\n");
    } else {
        printf("[");
        for (long long i = 0; i < n; i++)
            nul[i] ? printf("%snull", i ? "," : "") : printf("%s%.17g", i ? "," : "", v[i]);
        printf("]}\n");
    }
    fits_close_file(f, &status);
    return 0;
}
