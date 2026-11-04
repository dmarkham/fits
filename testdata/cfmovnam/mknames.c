/*
 * mknames: writes names.fits, the cfitsio-made input of the scenario where
 * HDUByName deliberately differs from fits_movnam_hdu (see the Makefile
 * and PROVENANCE.md).
 *
 *   mknames OUT
 *
 * Every HDU is a NAXIS = 0 IMAGE (headers only).
 *
 *   HDU 0  primary, no EXTNAME
 *   HDU 1  IMAGE  EXTNAME = 'sci'
 *   HDU 2  IMAGE  EXTNAME = 'SCI'
 *   HDU 3  IMAGE  HDUNAME = 'ERR', no EXTNAME
 *   HDU 4  IMAGE  EXTNAME = 'ERR'
 *   HDU 5  IMAGE  EXTNAME = 5 (an integer), then EXTNAME = 'NUM'
 *   HDU 6  IMAGE  EXTNAME = 'NUM'
 *
 * Exit status 0 on success; on a cfitsio error the status and message go
 * to stderr and the exit status is 1.
 */
#include <stdio.h>
#include <fitsio.h>

static int fail(int status)
{
    char msg[FLEN_ERRMSG];
    fits_get_errstatus(status, msg);
    fprintf(stderr, "mknames: status %d: %s\n", status, msg);
    while (fits_read_errmsg(msg))
        fprintf(stderr, "mknames:   %s\n", msg);
    return 1;
}

/* Every HDU is a NAXIS = 0 image: the lookup needs headers only, and
 * without data units the file is one block per HDU. */
static int image(fitsfile *f, short base, int *status)
{
    (void)base;
    return fits_create_img(f, SHORT_IMG, 0, NULL, status);
}

static int key(fitsfile *f, const char *name, const char *value, int *status)
{
    return fits_write_key(f, TSTRING, name, (void *)value, NULL, status);
}

int main(int argc, char **argv)
{
    fitsfile *f = NULL;
    int status = 0, five = 5;
    char name[FLEN_FILENAME + 1];

    if (argc != 2) {
        fprintf(stderr, "usage: mknames OUT\n");
        return 2;
    }
    snprintf(name, sizeof name, "!%s", argv[1]);
    if (fits_create_file(&f, name, &status))
        return fail(status);
    if (image(f, 0, &status))
        return fail(status);
    if (image(f, 10, &status) || key(f, "EXTNAME", "sci", &status))
        return fail(status);
    if (image(f, 20, &status) || key(f, "EXTNAME", "SCI", &status))
        return fail(status);
    if (image(f, 30, &status) || key(f, "HDUNAME", "ERR", &status))
        return fail(status);
    if (image(f, 40, &status) || key(f, "EXTNAME", "ERR", &status))
        return fail(status);
    if (image(f, 50, &status) ||
        fits_write_key(f, TINT, "EXTNAME", &five, "an integer EXTNAME", &status) ||
        key(f, "EXTNAME", "NUM", &status))
        return fail(status);
    if (image(f, 60, &status) || key(f, "EXTNAME", "NUM", &status))
        return fail(status);
    if (fits_close_file(f, &status))
        return fail(status);
    return 0;
}
