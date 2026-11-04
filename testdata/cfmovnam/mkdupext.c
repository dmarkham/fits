/*
 * mkdupext: writes dupext.fits, the cfitsio-made input of the duplicate
 * EXTNAME lookup scenario (see the Makefile and PROVENANCE.md).
 *
 *   mkdupext OUT
 *
 * fits_write_key always appends a new card, so writing EXTNAME twice gives
 * a header with two EXTNAME cards:
 *
 *   HDU 0  primary, no EXTNAME
 *   HDU 1  IMAGE  EXTNAME = 'SCI', then EXTNAME = 'ERR'
 *   HDU 2  IMAGE  EXTNAME = 'ERR'
 *   HDU 3  IMAGE  EXTNAME = 'DQ', then EXTNAME = 'SCI'
 *   HDU 4  BINTABLE  EXTNAME = 'TAB' twice
 *   HDU 5  IMAGE  EXTNAME = 'LAST', then EXTNAME = 'DQ'
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
    fprintf(stderr, "mkdupext: status %d: %s\n", status, msg);
    while (fits_read_errmsg(msg))
        fprintf(stderr, "mkdupext:   %s\n", msg);
    return 1;
}

static int image(fitsfile *f, const char *n1, const char *n2, short base, int *status)
{
    long naxes[1] = {3};
    short pix[3];
    int i;
    for (i = 0; i < 3; i++)
        pix[i] = (short)(base + i);
    fits_create_img(f, SHORT_IMG, 1, naxes, status);
    if (n1)
        fits_write_key(f, TSTRING, "EXTNAME", (void *)n1, "first EXTNAME card", status);
    if (n2)
        fits_write_key(f, TSTRING, "EXTNAME", (void *)n2, "second EXTNAME card", status);
    fits_write_img(f, TSHORT, 1, 3, pix, status);
    return *status;
}

int main(int argc, char **argv)
{
    fitsfile *f = NULL;
    int status = 0;
    char name[FLEN_FILENAME + 1];
    char *ttype[] = {"N"}, *tform[] = {"1J"}, *tunit[] = {""};
    int col[2] = {7, 8};

    if (argc != 2) {
        fprintf(stderr, "usage: mkdupext OUT\n");
        return 2;
    }
    snprintf(name, sizeof name, "!%s", argv[1]);
    if (fits_create_file(&f, name, &status))
        return fail(status);
    if (image(f, NULL, NULL, 0, &status))
        return fail(status);
    if (image(f, "SCI", "ERR", 10, &status))
        return fail(status);
    if (image(f, "ERR", NULL, 20, &status))
        return fail(status);
    if (image(f, "DQ", "SCI", 30, &status))
        return fail(status);
    if (fits_create_tbl(f, BINARY_TBL, 2, 1, ttype, tform, tunit, NULL, &status))
        return fail(status);
    if (fits_write_key(f, TSTRING, "EXTNAME", "TAB", "first EXTNAME card", &status) ||
        fits_write_key(f, TSTRING, "EXTNAME", "TAB", "second EXTNAME card", &status) ||
        fits_write_col(f, TINT, 1, 1, 1, 2, col, &status))
        return fail(status);
    if (image(f, "LAST", "DQ", 50, &status))
        return fail(status);
    if (fits_close_file(f, &status))
        return fail(status);
    return 0;
}
