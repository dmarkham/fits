/*
 * mkcomments: writes comments.fits, the cfitsio-made input of the
 * comments_1_0 scenario (see the Makefile and PROVENANCE.md).
 *
 *   mkcomments OUT
 *
 * When fits_copy_hdu turns a primary array into an IMAGE extension it
 * drops EXTEND and the two standard "FITS (Flexible Image Transport
 * System) format ..." COMMENT cards, recognising them by their first 8,
 * 58 and 47 characters (edithdu.c ffcphd). The primary written here
 * carries, besides the standard cards fits_create_img writes, COMMENT
 * cards that agree with those prefixes but differ later (cfitsio drops
 * them) or differ inside the prefix (cfitsio keeps them), and a keyword
 * EXTENDED that starts with EXTEND (kept). The IMAGE extension after it
 * carries the two standard COMMENT cards, which the reverse conversion
 * keeps after the ones it inserts (cfitsio does not look for duplicates).
 *
 * Exit status 0 on success; on a cfitsio error the status and message go
 * to stderr and the exit status is 1.
 */
#include <stdio.h>
#include <fitsio.h>

static const char *primaryCards[] = {
    /* same 58-character prefix as standard comment 1, different text */
    "COMMENT   FITS (Flexible Image Transport System) format is also an archive form",
    /* differs from standard comment 1 at character 58 */
    "COMMENT   FITS (Flexible Image Transport System) format iX",
    /* same 47-character prefix as standard comment 2, different text */
    "COMMENT   and Astrophysics', volume 376, page 3xx of another journal",
    /* differs from standard comment 2 at character 47 */
    "COMMENT   and Astrophysics', volume 376, page 4",
    "EXTENDED= 'not EXTEND'         / starts with EXTEND",
    NULL,
};

static const char *extCards[] = {
    "COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy",
    "COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H",
    NULL,
};

static int fail(int status)
{
    char msg[FLEN_ERRMSG];
    fits_get_errstatus(status, msg);
    fprintf(stderr, "mkcomments: status %d: %s\n", status, msg);
    while (fits_read_errmsg(msg))
        fprintf(stderr, "mkcomments:   %s\n", msg);
    return 1;
}

int main(int argc, char **argv)
{
    fitsfile *f = NULL;
    int status = 0, i;
    char name[FLEN_FILENAME + 1];
    long naxes[1] = {3};
    unsigned char pix[3] = {1, 2, 3};
    short pix2[3] = {-1, 0, 1};

    if (argc != 2) {
        fprintf(stderr, "usage: mkcomments OUT\n");
        return 2;
    }
    snprintf(name, sizeof name, "!%s", argv[1]);
    fits_create_file(&f, name, &status);
    fits_create_img(f, BYTE_IMG, 1, naxes, &status);
    for (i = 0; primaryCards[i]; i++)
        fits_write_record(f, primaryCards[i], &status);
    fits_write_img(f, TBYTE, 1, 3, pix, &status);
    fits_create_img(f, SHORT_IMG, 1, naxes, &status);
    for (i = 0; extCards[i]; i++)
        fits_write_record(f, extCards[i], &status);
    fits_write_img(f, TSHORT, 1, 3, pix2, &status);
    fits_close_file(f, &status);
    if (status)
        return fail(status);
    return 0;
}
