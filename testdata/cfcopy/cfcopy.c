/*
 * cfcopy: cfitsio truth harness for EditFile's HDU copy / reorder / skip.
 *
 *   cfcopy [-r MOREKEYS] OUT IN HDU [HDU ...]
 *
 * Creates OUT and copies the listed HDUs of IN (0-based, in the given
 * order, repeats allowed) into it with fits_copy_hdu. cfitsio converts a
 * primary array that lands after position 0 into an IMAGE extension, an
 * IMAGE extension that lands at position 0 into a primary array, and puts
 * a dummy (NAXIS = 0) primary array in front of a table copied to position
 * 0 (edithdu.c ffcphd). The Go test compares its own EditFile output with
 * the cfitsio dump (testdata/cfref) of the file this program writes.
 *
 * -r MOREKEYS is passed to fits_copy_hdu as its morekeys argument: cfitsio
 * reserves room for that many more keywords in each copied header, which
 * shows up in the file as blank cards before END. The Makefile uses it to
 * make reserve.fits, a source whose headers carry reserved space.
 *
 * Exit status 0 on success; on a cfitsio error the status and message go
 * to stderr and the exit status is 1.
 */
#include <stdio.h>
#include <stdlib.h>
#include <fitsio.h>

static int fail(int status, const char *what)
{
    char msg[FLEN_ERRMSG];
    fits_get_errstatus(status, msg);
    fprintf(stderr, "cfcopy: %s: status %d: %s\n", what, status, msg);
    while (fits_read_errmsg(msg))
        fprintf(stderr, "cfcopy:   %s\n", msg);
    return 1;
}

int main(int argc, char **argv)
{
    fitsfile *in = NULL, *out = NULL;
    int status = 0, i;
    char outname[FLEN_FILENAME + 1];

    int morekeys = 0;

    if (argc > 2 && argv[1][0] == '-' && argv[1][1] == 'r' && argv[1][2] == 0) {
        morekeys = atoi(argv[2]);
        argv += 2;
        argc -= 2;
    }
    if (argc < 4) {
        fprintf(stderr, "usage: cfcopy [-r MOREKEYS] OUT IN HDU [HDU ...]\n");
        return 2;
    }
    /* "!" makes cfitsio overwrite an existing OUT. */
    snprintf(outname, sizeof outname, "!%s", argv[1]);
    if (fits_open_file(&in, argv[2], READONLY, &status))
        return fail(status, "open input");
    if (fits_create_file(&out, outname, &status))
        return fail(status, "create output");
    for (i = 3; i < argc; i++) {
        int hdu = atoi(argv[i]);
        if (fits_movabs_hdu(in, hdu + 1, NULL, &status))
            return fail(status, "move to input HDU");
        if (fits_copy_hdu(in, out, morekeys, &status))
            return fail(status, "copy HDU");
    }
    if (fits_close_file(out, &status))
        return fail(status, "close output");
    if (fits_close_file(in, &status))
        return fail(status, "close input");
    return 0;
}
