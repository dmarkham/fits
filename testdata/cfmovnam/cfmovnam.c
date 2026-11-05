/*
 * cfmovnam: cfitsio truth for HDUByName (fits_movnam_hdu).
 *
 *   cfmovnam IN NAME [NAME ...]
 *
 * For each NAME, opens IN and calls fits_movnam_hdu(ANY_HDU, NAME, 0),
 * which scans the HDUs from the primary and reads each one's EXTNAME with
 * fits_read_key_str (the first EXTNAME card of the header). It prints one
 * JSON object: {"file": IN's base name, "hdu": {NAME: 0-based HDU index,
 * or -1 when cfitsio reports BAD_HDU_NUM (no match)}} in argument order.
 *
 * Exit status 0 on success; on any other cfitsio error the status and
 * message go to stderr and the exit status is 1.
 */
#include <stdio.h>
#include <string.h>
#include <fitsio.h>

static int fail(int status, const char *what)
{
    char msg[FLEN_ERRMSG];
    fits_get_errstatus(status, msg);
    fprintf(stderr, "cfmovnam: %s: status %d: %s\n", what, status, msg);
    while (fits_read_errmsg(msg))
        fprintf(stderr, "cfmovnam:   %s\n", msg);
    return 1;
}

int main(int argc, char **argv)
{
    int i;
    const char *base;

    if (argc < 3) {
        fprintf(stderr, "usage: cfmovnam IN NAME [NAME ...]\n");
        return 2;
    }
    base = strrchr(argv[1], '/');
    base = base ? base + 1 : argv[1];
    printf("{\"file\":\"%s\",\"hdu\":{", base);
    for (i = 2; i < argc; i++) {
        fitsfile *f = NULL;
        int status = 0, hdu = -1;
        if (fits_open_file(&f, argv[1], READONLY, &status))
            return fail(status, "open");
        if (fits_movnam_hdu(f, ANY_HDU, argv[i], 0, &status) == 0) {
            fits_get_hdu_num(f, &hdu);
            hdu--;
        } else if (status == BAD_HDU_NUM) {
            status = 0;
            fits_clear_errmsg();
        } else {
            return fail(status, "movnam");
        }
        if (fits_close_file(f, &status))
            return fail(status, "close");
        printf("%s\"%s\":%d", i > 2 ? "," : "", argv[i], hdu);
    }
    printf("}}\n");
    return 0;
}
