/*
 * gen_fpack.c -- deterministic source images and the cfitsio-only
 * compressed variants for the tile-compression READ corpus
 * (testdata/comp_fz_*.fits). Driven by testdata/gen_fpack.sh; see
 * testdata/COMP_FZ_PROVENANCE.md.
 *
 *   gen_fpack src <dir>
 *       write every source image <dir>/<name>.fits with libcfitsio
 *       (fits_create_img + fits_write_img, plus BLANK/BSCALE/BZERO keywords
 *       where the case needs them). Pixel values come from a fixed 64-bit
 *       LCG, so the bytes are identical on every run.
 *
 *   gen_fpack hcomp <in.fits> <out.fits> <scale> <smooth> [qlevel seed]
 *       HCOMPRESS_1 the primary image of <in.fits> with fits_img_compress,
 *       setting fits_set_hcomp_scale(scale) and fits_set_hcomp_smooth(smooth).
 *       fpack has no command-line flag for SMOOTH, which is why this exists.
 *       For float input, qlevel/seed set the quantization level and the
 *       SUBTRACTIVE_DITHER_1 seed (fits_set_dither_seed).
 *
 *   gen_fpack zblankcol <in.fz> <out.fits> [1K <value>]
 *       copy a cfitsio tile-compressed file and move its null value from the
 *       ZBLANK (or BLANK) header keyword into a per-tile ZBLANK table column
 *       (1J), as the tiled-image convention allows. Row r gets the value
 *       blank + (r - 1) * step where step comes from the "ZBLKSTEP" keyword
 *       (written by the source generator, 0 for float images), so each tile
 *       really has its own null value. With "1K <value>" after the file
 *       names the column is 1K and tile 1's entry is <value> (a null value
 *       outside the C int that cfitsio reads it as).
 *
 *   gen_fpack editkey <in.fits> <out.fits> <KEY> [value]
 *       byte copy of <in.fits> with keyword KEY of HDU 2 deleted (no value)
 *       or set to the floating-point value. Used to make a
 *       SUBTRACTIVE_DITHER_1 file without ZDITHER0 (as written before
 *       October 2009; cfitsio then starts the dither at row 1) and a
 *       quantized float image carrying BSCALE/BZERO (fpack refuses to write
 *       one; cfitsio reads it by combining the scalings).
 *
 *   gen_fpack setkey <in.fits> <out.fits> <KEY> int|str <value>
 *       byte copy of <in.fits> with keyword KEY of HDU 2 set to an integer
 *       or a string value (ZBLANK next to BLANK, a ZQUANTIZ cfitsio does
 *       not know, swapped ZNAMEn).
 *
 *   gen_fpack uncompcol <in.fz> <out.fits> [row ...]
 *       write tiles into the pre-2011 UNCOMPRESSED_DATA column, holding the
 *       raw stored values (read back through cfitsio with BSCALE/BZERO
 *       turned off) in the column type of ZBITPIX (1PB/1PI/1PJ/1PE/1PD).
 *       Without rows: the GZIP_COMPRESSED_DATA fallback tiles of a float
 *       image are moved and that column is deleted. With rows (1-based):
 *       those tiles are moved out of COMPRESSED_DATA (its descriptor is
 *       set to length 0), for integer images too.
 *
 *   gen_fpack rewidth <in.fz> <out.fits> <W> [extra] [bump]
 *       copy an integer GZIP_1, GZIP_2 or RICE_1 file with row tiles and
 *       re-encode every tile with W-byte tile integers (W = 1 unsigned, 2,
 *       4 or 8 signed, big-endian), using libcfitsio's own coders
 *       (compress2mem_from_mem, fits_rcomp / fits_rcomp_short /
 *       fits_rcomp_byte). The tile values are the stored values read back
 *       through cfitsio with BSCALE/BZERO off. This is the layout old
 *       cfitsio versions wrote (e.g. 16-bit images in 4-byte GZIP tiles,
 *       or RICE_1 with BYTEPIX 4): cfitsio's reader picks the tile type
 *       from the inflated length (GZIP) or BYTEPIX (RICE, ZVAL2 is set to
 *       W). <extra> appends that many zero bytes to every GZIP tile before
 *       deflating (cfitsio refuses such a tile); <bump> is added to the
 *       first value of tile 1 (to make a value the image type cannot hold);
 *       <i:v,...> sets value i of tile 1 to v (the edges of the ZBITPIX
 *       range, and one past them); R@i:v sets value i of tile R. A
 *       ZBITPIX 64 input (written by astropy) is read with TLONG.
 *
 *   gen_fpack rewidthq <in.fz> <out.fits> <W>
 *       the same for the tiles of a quantized float image: RICE_1 with W =
 *       1 or 2 byte Rice coding of the stored quantized integers, GZIP_1 /
 *       GZIP_2 with the integers deflated (shuffled for GZIP_2) as W = 1,
 *       2 or 8 byte values instead of 4.
 *
 *   gen_fpack emptytile <in> <out> <row>
 *       COMPRESSED_DATA of one tile set to length 0, with no fallback column.
 *
 *   gen_fpack gzipcol <in.fz> <out> <row>
 *       one integer tile moved into a GZIP_COMPRESSED_DATA column (which
 *       cfitsio only writes, and only reads, for float images).
 *
 *   gen_fpack croprows <in> <out> <N>
 *       the first N tiles of a row-tiled compressed image, byte for byte
 *       (the real M13 files, cut to 40 of their 300 rows).
 *
 *   gen_fpack mk64 <out.fits> RICE_1|HCOMPRESS_1|PLIO_1 [V]
 *       a ZBITPIX 64 tile-compressed image built from libcfitsio's own
 *       coders (fpack and fits_img_compress refuse 64-bit images).
 *
 *   uncompcol also takes a negative row -r: tile r moved with its last
 *   value left out (an UNCOMPRESSED_DATA tile of the wrong length).
 *
 * Nothing here uses the Go library.
 */
#define _DEFAULT_SOURCE
#include <fitsio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int compress2mem_from_mem(char *inmemptr, size_t inmemsize, char **buffptr,
                          size_t *buffsize,
                          void *(*mem_realloc)(void *p, size_t newsize),
                          size_t *filesize, int *status);
int uncompress2mem_from_mem(char *inmemptr, size_t inmemsize, char **buffptr,
                            size_t *buffsize,
                            void *(*mem_realloc)(void *p, size_t newsize),
                            size_t *filesize, int *status);
int fits_rcomp(int a[], int nx, unsigned char *c, int clen, int nblock);
int fits_rcomp_short(short a[], int nx, unsigned char *c, int clen, int nblock);
int fits_rcomp_byte(signed char a[], int nx, unsigned char *c, int clen,
                    int nblock);
int fits_rdecomp(unsigned char *c, int clen, unsigned int array[], int nx,
                 int nblock);
int pl_p2li(int *pxsrc, int xs, short *lldst, int npix);

#define NX 40
#define NY 24
#define N (NX * NY)

static unsigned long long st = 12345;
static double U(void) {
  st = st * 6364136223846793005ULL + 1442695040888963407ULL;
  return (double)(st >> 33) / 4294967296.0;
}
static double G(void) {
  double a = U(), b = U();
  if (a < 1e-12) a = 1e-12;
  return sqrt(-2 * log(a)) * cos(2 * M_PI * b);
}

static void chk(int s, const char *w) {
  if (s) {
    fprintf(stderr, "%s: ", w);
    fits_report_error(stderr, s);
    exit(1);
  }
}

typedef void (*kwfn)(fitsfile *);

static void img(const char *dir, const char *name, int bitpix, int dtype,
                void *data, long n, long *ax, int nax, kwfn kw) {
  fitsfile *f;
  int s = 0;
  char p[1024];
  snprintf(p, sizeof p, "!%s/%s.fits", dir, name);
  fits_create_file(&f, p, &s);
  fits_create_img(f, bitpix, nax, ax, &s);
  if (kw) kw(f);
  fits_write_img(f, dtype, 1, n, data, &s);
  fits_close_file(f, &s);
  chk(s, name);
}

static void kw_blank16(fitsfile *f) {
  int s = 0;
  long b = -32768;
  fits_write_key(f, TLONG, "BLANK", &b, "null value", &s);
  chk(s, "BLANK");
}
static void kw_blank65541(fitsfile *f) {
  int s = 0;
  long b = 65541;
  fits_write_key(f, TLONG, "BLANK", &b, "null value", &s);
  chk(s, "BLANK");
}
static void kw_blank_u16plio(fitsfile *f) {
  int s = 0;
  long b = -32743;
  fits_write_key(f, TLONG, "BLANK", &b, "null value (stored)", &s);
  chk(s, "BLANK");
}
/* per-tile blank: tile (row) r uses -32768 + r*7 (r 0-based). */
static void kw_blankstep(fitsfile *f) {
  int s = 0;
  long b = -32768, step = 7;
  fits_write_key(f, TLONG, "BLANK", &b, "null value of row 1", &s);
  fits_write_key(f, TLONG, "ZBLKSTEP", &step, "gen_fpack: per-row null step", &s);
  chk(s, "BLANK/ZBLKSTEP");
}
static void kw_scaled(fitsfile *f) {
  int s = 0;
  double bs = 0.5, bz = 100;
  fits_write_key(f, TDOUBLE, "BSCALE", &bs, "", &s);
  fits_write_key(f, TDOUBLE, "BZERO", &bz, "", &s);
  fits_set_bscale(f, 1.0, 0.0, &s); /* write the stored integers as given */
  chk(s, "BSCALE/BZERO");
}

static int do_src(const char *dir) {
  long ax[2] = {NX, NY};
  long i;
  static unsigned char u8[N];
  static signed char s8[N];
  static short i16[N], i16b[N], i16s[N], i16z[N], i16m[N];
  static unsigned short u16[N];
  static int i32[N], i32m[N];
  static unsigned int u32[N];
  static float f32[N], f32n[N], f32z[N], f32c[N];
  static double f64[N], f64n[N], f64c[N], f64z[N];
  for (i = 0; i < N; i++) {
    long x = i % NX, y = i / NX;
    double g = G();
    double smooth = 30 * sin(x / 6.0) * cos(y / 5.0); /* structure for HCOMPRESS */
    u8[i] = (unsigned char)((int)(128 + 40 * g + (i % 50)) & 255);
    s8[i] = (signed char)(int)fmax(-128, fmin(127, 30 * g));
    i16[i] = (short)(int)(200 * g - 100 + (i % 300) + 10 * smooth);
    u16[i] = (unsigned short)(int)fmin(65535, fmax(0, 40000 + 3000 * g + (i % 5000)));
    i32[i] = (int)(1e5 * g) - 123456 + (int)(1000 * smooth);
    u32[i] = (unsigned int)(3000000000.0 + 1e6 * g);
    f32[i] = (float)(100 + 5 * g + smooth);
    f32n[i] = (i % 31 == 0) ? NAN : f32[i];
    f32z[i] = (i % 23 == 0) ? 0.0f : f32[i];
    /* row 5 constant: zero noise -> cfitsio cannot quantize that tile and
       stores it in the GZIP_COMPRESSED_DATA fallback column */
    f32c[i] = (y == 5) ? 7.25f : f32[i];
    f64[i] = 1000 + 20 * g + smooth;
    f64n[i] = (i % 37 == 0) ? NAN : f64[i];
    f64c[i] = (y == 5) ? 7.25 : f64[i]; /* as f32flat, for BITPIX -64 */
    f64z[i] = (i % 23 == 0) ? 0.0 : f64[i]; /* as f32zero, for BITPIX -64 SD2 */
    /* PLIO at ZBITPIX 32: runs of 5 pixels with values across 0 .. 2^24-1
       (PLIO's range), including the 12-bit set-value boundary and the
       maximum. (Pixel-to-pixel noise makes PLIO expand the data, and
       fpack 1.7.0 (CFITSIO 4.7.0) then overruns its buffer; see COMP_FZ_PROVENANCE.md.) */
    i32m[i] = (y < 4) ? 0 : (y < 8) ? (int)(((x / 5) * 1048573L + y * 4099L) % 16777216L)
            : (y < 12) ? 16777215 : (x < 20) ? 4095 : 4096;
    i16b[i] = (i % 29 == 0) ? -32768 : i16[i];
    i16s[i] = (short)(int)(50 * g);
    /* PLIO needs small non-negative values */
    i16m[i] = (short)((x > 5 && x < 20 && y > 3 && y < 12) ? 3 : (x + y) % 9 == 0 ? 1 : 0);
    /* per-row null: -32768 + 7*y every 13th pixel */
    i16z[i] = (i % 13 == 0) ? (short)(-32768 + 7 * y) : i16[i];
  }
  img(dir, "u8", BYTE_IMG, TBYTE, u8, N, ax, 2, 0);
  img(dir, "s8", SBYTE_IMG, TSBYTE, s8, N, ax, 2, 0);
  img(dir, "i16", SHORT_IMG, TSHORT, i16, N, ax, 2, 0);
  img(dir, "u16", USHORT_IMG, TUSHORT, u16, N, ax, 2, 0);
  img(dir, "i16blank", SHORT_IMG, TSHORT, i16b, N, ax, 2, kw_blank16);
  img(dir, "i16rowblank", SHORT_IMG, TSHORT, i16z, N, ax, 2, kw_blankstep);
  img(dir, "i16scaled", SHORT_IMG, TSHORT, i16s, N, ax, 2, kw_scaled);
  img(dir, "i16mask", SHORT_IMG, TSHORT, i16m, N, ax, 2, 0);
  img(dir, "i32", LONG_IMG, TINT, i32, N, ax, 2, 0);
  img(dir, "i32mask", LONG_IMG, TINT, i32m, N, ax, 2, 0);
  img(dir, "u32", ULONG_IMG, TUINT, u32, N, ax, 2, 0);
  img(dir, "f32", FLOAT_IMG, TFLOAT, f32, N, ax, 2, 0);
  img(dir, "f32nan", FLOAT_IMG, TFLOAT, f32n, N, ax, 2, 0);
  img(dir, "f32zero", FLOAT_IMG, TFLOAT, f32z, N, ax, 2, 0);
  img(dir, "f32flat", FLOAT_IMG, TFLOAT, f32c, N, ax, 2, 0);
  img(dir, "f64", DOUBLE_IMG, TDOUBLE, f64, N, ax, 2, 0);
  img(dir, "f64nan", DOUBLE_IMG, TDOUBLE, f64n, N, ax, 2, 0);
  img(dir, "f64flat", DOUBLE_IMG, TDOUBLE, f64c, N, ax, 2, 0);
  img(dir, "f64zero", DOUBLE_IMG, TDOUBLE, f64z, N, ax, 2, 0);
  {
    long a3[3] = {20, 12, 4};
    static short c[20 * 12 * 4];
    static float cf[20 * 12 * 4];
    for (i = 0; i < 20 * 12 * 4; i++) {
      c[i] = (short)(i % 997 - 300);
      cf[i] = (float)(50 + 3 * G() + (i / 240) * 10);
    }
    img(dir, "cube16", SHORT_IMG, TSHORT, c, 20 * 12 * 4, a3, 3, 0);
    img(dir, "cubef32", FLOAT_IMG, TFLOAT, cf, 20 * 12 * 4, a3, 3, 0);
  }
  /* Added by the fourth proof pass, after every G() call above, so the
     earlier sources are unchanged. */
  {
    /* uint16 PLIO with a BLANK keyword: runs of small physical values
       0..39 (PLIO-friendly), BLANK = -32743 is the STORED value of the
       physical 25 that 24 pixels hold. cfitsio compares the null value
       with the PLIO-decoded value (stored + 32768), not the stored one. */
    static unsigned short u16p[N];
    for (i = 0; i < N; i++) {
      long x = i % NX, y = i / NX;
      u16p[i] = (unsigned short)((x / 5 + y) % 40);
    }
    img(dir, "u16pliob", USHORT_IMG, TUSHORT, u16p, N, ax, 2, kw_blank_u16plio);
    /* 4-D and 6-D int16 images (cfitsio compresses up to 6 dimensions,
       MAX_COMPRESS_DIM). */
    long a4[4] = {10, 6, 4, 2}, a6[6] = {5, 4, 3, 2, 2, 2};
    static short c4[10 * 6 * 4 * 2], c6[5 * 4 * 3 * 2 * 2 * 2];
    for (i = 0; i < 10 * 6 * 4 * 2; i++) c4[i] = (short)((i * 37) % 1001 - 400);
    for (i = 0; i < 5 * 4 * 3 * 2 * 2 * 2; i++) c6[i] = (short)((i * 53) % 777 - 300);
    img(dir, "cube4d16", SHORT_IMG, TSHORT, c4, 10 * 6 * 4 * 2, a4, 4, 0);
    img(dir, "cube6d16", SHORT_IMG, TSHORT, c6, 5 * 4 * 3 * 2 * 2 * 2, a6, 6, 0);
  }
  {
    /* A 500-row binary table for fpack -table (tiled TABLE compression,
       ZTABLE = T), which the library refuses; fpack leaves smaller tables
       uncompressed. */
    enum { NR = 500 };
    static float flux[NR];
    static int id[NR];
    static char names[NR][9], *np[NR];
    for (i = 0; i < NR; i++) {
      flux[i] = (float)(10 + (i % 37) * 0.25);
      id[i] = (int)(1000 + 3 * i);
      snprintf(names[i], sizeof names[i], "S%05ld", i);
      np[i] = names[i];
    }
    char path[4096];
    fitsfile *f;
    int s = 0;
    char *tt[] = {"FLUX", "ID", "NAME"}, *tf[] = {"1E", "1J", "8A"};
    snprintf(path, sizeof path, "!%s/tbl500.fits", dir);
    fits_create_file(&f, path, &s);
    fits_create_tbl(f, BINARY_TBL, NR, 3, tt, tf, NULL, "SOURCES", &s);
    fits_write_col(f, TFLOAT, 1, 1, 1, NR, flux, &s);
    fits_write_col(f, TINT, 2, 1, 1, NR, id, &s);
    fits_write_col(f, TSTRING, 3, 1, 1, NR, np, &s);
    fits_close_file(f, &s);
    chk(s, "tbl500");
  }
  /* Added by the round-2 proof pass, after every earlier G() call, so the
     earlier sources are unchanged. */
  {
    /* f64wrap: one 102x100 tile of 10200 pixels (more than the 10000
       random values of the dither sequence), so SUBTRACTIVE_DITHER_1
       dequantization must wrap the sequence inside the tile for
       BITPIX -64. */
    enum { WX = 102, WY = 100 };
    static double fw[WX * WY];
    long aw[2] = {WX, WY};
    for (i = 0; i < WX * WY; i++)
      fw[i] = 1000 + 20 * G() + 30 * sin((i % WX) / 6.0) * cos((i / WX) / 5.0);
    img(dir, "f64wrap", DOUBLE_IMG, TDOUBLE, fw, WX * WY, aw, 2, 0);
    /* f32flatnan / f64flatnan: f32flat / f64flat with NaN at x = 3 and 17
       of the constant row 5 (cfitsio cannot quantize that row and stores
       it, NaN included, in the GZIP_COMPRESSED_DATA fallback) and at x =
       8 of row 2 (a quantized row: ZBLANK). */
    static float f32fn[N];
    static double f64fn[N];
    for (i = 0; i < N; i++) {
      long x = i % NX, y = i / NX;
      int nan = (y == 5 && (x == 3 || x == 17)) || (y == 2 && x == 8);
      f32fn[i] = nan ? NAN : f32c[i];
      f64fn[i] = nan ? NAN : f64c[i];
    }
    img(dir, "f32flatnan", FLOAT_IMG, TFLOAT, f32fn, N, ax, 2, 0);
    img(dir, "f64flatnan", DOUBLE_IMG, TDOUBLE, f64fn, N, ax, 2, 0);
    /* i16small / i32small (40x6): stored values 0..255, so they fit 1-,
       2- and 4-byte tile integers (gen_fpack rewidth), with zeros and
       fives. BLANK -32768 on i16small and 65541 on i32small match no
       stored value, but cfitsio compares the null value with the tile
       integers after casting it to the tile type ((unsigned char) -32768
       = 0, (unsigned char) and (short) 65541 = 5; imcompress.c fffi1i2 /
       fffi2int calls), so narrow tiles have nulls. */
    enum { SY = 6 };
    static short s16[NX * SY];
    static int s32[NX * SY];
    long as[2] = {NX, SY};
    for (i = 0; i < NX * SY; i++) {
      s16[i] = (short)((i * 37 + (i / NX) * 11) % 256);
      s32[i] = (int)((i * 53 + 5) % 256);
      if (i % 17 == 0) s16[i] = 0, s32[i] = 5;
    }
    img(dir, "i16small", SHORT_IMG, TSHORT, s16, NX * SY, as, 2, kw_blank16);
    img(dir, "i32small", LONG_IMG, TINT, s32, NX * SY, as, 2, kw_blank65541);
  }
  /* Added by the round-2 fix pass (its own LCG, so no earlier source
     changes). u8sat / i16sat / u16sat / i32sat (40x6): noise of +-amp just
     below the type maximum and a 16x4 block AT the maximum. Lossy
     HCOMPRESS (fpack -h -s 2) reconstructs pixels of that block above the
     maximum: cfitsio clips them for 16-bit images (imcompress.c resets the
     overflow for HCOMPRESS_1), funpack fails for 8-bit (status 412), and
     the 32-bit reconstruction wraps in fits_hdecompress64's (int) cast. */
  {
    enum { SY = 6 };
    static unsigned char a8[NX * SY];
    static short a16[NX * SY];
    static unsigned short au16[NX * SY];
    static int a32[NX * SY];
    long as[2] = {NX, SY};
    unsigned long lcg = 12345;
    for (i = 0; i < NX * SY; i++) {
      long x = i % NX, y = i / NX, r;
      int blk = x >= 10 && x < 26 && y >= 1 && y < 5;
      lcg = (lcg * 1103515245UL + 12345UL) % 2147483648UL;
      r = (long)lcg;
      a8[i] = (unsigned char)(blk ? 255 : 255 - 21 + (r % 41 - 20));
      a16[i] = (short)(blk ? 32767 : 32767 - 301 + (r % 601 - 300));
      au16[i] = (unsigned short)(blk ? 65535 : 65535 - 301 + (r % 601 - 300));
      a32[i] = (int)(blk ? 2147483647L : 2147483647L - 301 + (r % 601 - 300));
    }
    img(dir, "u8sat", BYTE_IMG, TBYTE, a8, NX * SY, as, 2, 0);
    img(dir, "i16sat", SHORT_IMG, TSHORT, a16, NX * SY, as, 2, 0);
    img(dir, "u16sat", USHORT_IMG, TUSHORT, au16, NX * SY, as, 2, 0);
    img(dir, "i32sat", LONG_IMG, TINT, a32, NX * SY, as, 2, 0);
  }
  return 0;
}

static int do_hcomp(int argc, char **argv) {
  fitsfile *in, *out;
  int s = 0;
  char p[1024];
  if (argc < 6) {
    fprintf(stderr, "usage: gen_fpack hcomp in out scale smooth [qlevel seed]\n");
    return 2;
  }
  fits_open_file(&in, argv[2], READONLY, &s);
  snprintf(p, sizeof p, "!%s", argv[3]);
  fits_create_file(&out, p, &s);
  fits_create_img(out, 8, 0, NULL, &s);
  fits_set_compression_type(out, HCOMPRESS_1, &s);
  {
    long t[2] = {NX, 8}; /* 3 tiles of 40x8 */
    fits_set_tile_dim(out, 2, t, &s);
  }
  fits_set_hcomp_scale(out, (float)atof(argv[4]), &s);
  fits_set_hcomp_smooth(out, atoi(argv[5]), &s);
  if (argc >= 8) {
    fits_set_quantize_level(out, (float)atof(argv[6]), &s);
    fits_set_quantize_method(out, SUBTRACTIVE_DITHER_1, &s);
    fits_set_dither_seed(out, atoi(argv[7]), &s);
  }
  fits_img_compress(in, out, &s);
  fits_close_file(out, &s);
  fits_close_file(in, &s);
  chk(s, "hcomp");
  return 0;
}

static int do_zblankcol(int argc, char **argv) {
  fitsfile *f;
  int s = 0, col = 0, wide;
  long blank = 0, step = 0, nrows = 0, r;
  char cmd[2048];
  if (argc < 4) {
    fprintf(stderr, "usage: gen_fpack zblankcol in.fz out.fits\n");
    return 2;
  }
  /* work on a byte copy, opened as a plain binary table ([compress] off) */
  snprintf(cmd, sizeof cmd, "cp '%s' '%s'", argv[2], argv[3]);
  if (system(cmd) != 0) return 1;
  /* fits_open_data with the HDU opened explicitly; cfitsio exposes the
     underlying BINTABLE to the table routines. */
  fits_open_file(&f, argv[3], READWRITE, &s);
  fits_movabs_hdu(f, 2, NULL, &s);
  chk(s, "open");
  {
    int ts = 0;
    if (fits_read_key(f, TLONG, "ZBLANK", &blank, NULL, &ts) == 0) {
      fits_delete_key(f, "ZBLANK", &s);
    } else {
      ts = 0;
      fits_read_key(f, TLONG, "BLANK", &blank, NULL, &ts);
      chk(ts, "no ZBLANK/BLANK keyword");
      fits_delete_key(f, "BLANK", &s);
    }
    ts = 0;
    if (fits_read_key(f, TLONG, "ZBLKSTEP", &step, NULL, &ts) != 0) step = 0;
  }
  fits_get_num_rows(f, &nrows, &s);
  fits_get_num_cols(f, &col, &s);
  wide = argc >= 6 && !strcmp(argv[4], "1K");
  fits_insert_col(f, col + 1, "ZBLANK", wide ? "1K" : "1J", &s);
  chk(s, "insert ZBLANK column");
  for (r = 1; r <= nrows; r++) {
    if (wide) { /* 1K column, tile 1's null value given (outside int32) */
      LONGLONG v = r == 1 ? atoll(argv[5]) : blank + (r - 1) * step;
      fits_write_col(f, TLONGLONG, col + 1, r, 1, 1, &v, &s);
    } else {
      int v = (int)(blank + (r - 1) * step);
      fits_write_col(f, TINT, col + 1, r, 1, 1, &v, &s);
    }
  }
  fits_close_file(f, &s);
  chk(s, "zblankcol");
  return 0;
}

static int do_editkey(int argc, char **argv) {
  fitsfile *f;
  int s = 0;
  char cmd[2048];
  if (argc < 5) {
    fprintf(stderr, "usage: gen_fpack editkey in out KEY [value]\n");
    return 2;
  }
  snprintf(cmd, sizeof cmd, "cp '%s' '%s'", argv[2], argv[3]);
  if (system(cmd) != 0) return 1;
  fits_open_file(&f, argv[3], READWRITE, &s);
  fits_movabs_hdu(f, 2, NULL, &s);
  if (argc >= 6) {
    double v = atof(argv[5]);
    fits_update_key(f, TDOUBLE, argv[4], &v, "gen_fpack editkey", &s);
  } else {
    fits_delete_key(f, argv[4], &s);
  }
  fits_close_file(f, &s);
  chk(s, "editkey");
  return 0;
}

static int do_setkey(int argc, char **argv) {
  fitsfile *f;
  int s = 0;
  char cmd[2048];
  if (argc < 7) {
    fprintf(stderr, "usage: gen_fpack setkey in out KEY int|str value\n");
    return 2;
  }
  snprintf(cmd, sizeof cmd, "cp '%s' '%s'", argv[2], argv[3]);
  if (system(cmd) != 0) return 1;
  fits_open_file(&f, argv[3], READWRITE, &s);
  fits_movabs_hdu(f, 2, NULL, &s);
  if (!strcmp(argv[5], "int")) {
    long v = atol(argv[6]);
    fits_update_key(f, TLONG, argv[4], &v, "gen_fpack setkey", &s);
  } else {
    fits_update_key(f, TSTRING, argv[4], argv[6], "gen_fpack setkey", &s);
  }
  fits_close_file(f, &s);
  chk(s, "setkey");
  return 0;
}

/* patch in out OLD NEW: copy in to out, replacing the one occurrence of
   the byte string OLD (on an 80-byte card boundary) with NEW of the same
   length. For header values cfitsio itself refuses to write. */
static int do_patch(int argc, char **argv) {
  FILE *fp;
  static char buf[1 << 22];
  size_t n, len, i, hit = 0, at = 0;
  if (argc < 6 || strlen(argv[4]) != strlen(argv[5])) {
    fprintf(stderr, "usage: gen_fpack patch in out OLD NEW (same length)\n");
    return 2;
  }
  fp = fopen(argv[2], "rb");
  if (!fp) return 1;
  n = fread(buf, 1, sizeof buf, fp);
  fclose(fp);
  len = strlen(argv[4]);
  for (i = 0; i + len <= n; i += 80)
    if (!memcmp(buf + i, argv[4], len)) hit++, at = i;
  if (hit != 1) {
    fprintf(stderr, "gen_fpack patch: '%s' found %zu times\n", argv[4], hit);
    return 1;
  }
  memcpy(buf + at, argv[5], len);
  fp = fopen(argv[3], "wb");
  if (!fp || fwrite(buf, 1, n, fp) != n) return 1;
  return fclose(fp) != 0;
}

static int do_uncompcol(int argc, char **argv) {
  fitsfile *in, *f;
  int s = 0, gcol = 0, ccol = 0, ncol = 0, anyn = 0, zbitpix = 0, dtype = 0;
  long nrows = 0, r, nx = 0, ztile2 = 0;
  const char *tform;
  char cmd[2048];
  if (argc < 4) {
    fprintf(stderr, "usage: gen_fpack uncompcol in.fz out.fits [row ...]\n");
    return 2;
  }
  snprintf(cmd, sizeof cmd, "cp '%s' '%s'", argv[2], argv[3]);
  if (system(cmd) != 0) return 1;
  fits_open_file(&in, argv[2], READONLY, &s);
  fits_movabs_hdu(in, 2, NULL, &s);
  fits_set_bscale(in, 1.0, 0.0, &s); /* raw stored values */
  fits_open_file(&f, argv[3], READWRITE, &s);
  fits_movabs_hdu(f, 2, NULL, &s);
  fits_read_key(f, TLONG, "ZNAXIS1", &nx, NULL, &s);
  fits_read_key(f, TLONG, "ZTILE2", &ztile2, NULL, &s);
  fits_read_key(f, TINT, "ZBITPIX", &zbitpix, NULL, &s);
  chk(s, "open");
  if (ztile2 != 1) {
    fprintf(stderr, "uncompcol: expects row tiles (ZTILE2 = 1)\n");
    return 1;
  }
  switch (zbitpix) {
  case 8: tform = "1PB"; dtype = TBYTE; break;
  case 16: tform = "1PI"; dtype = TSHORT; break;
  case 32: tform = "1PJ"; dtype = TINT; break;
  case -32: tform = "1PE"; dtype = TFLOAT; break;
  case -64: tform = "1PD"; dtype = TDOUBLE; break;
  default: fprintf(stderr, "uncompcol: ZBITPIX %d\n", zbitpix); return 1;
  }
  fits_get_num_rows(f, &nrows, &s);
  {
    char *move = calloc(nrows, 1);
    double *row = calloc(nx, sizeof(double)); /* big enough for any type */
    int i;
    if (argc == 4) {
      fits_get_colnum(f, CASEINSEN, "GZIP_COMPRESSED_DATA", &gcol, &s);
      chk(s, "no GZIP_COMPRESSED_DATA column");
      for (r = 1; r <= nrows; r++) {
        long len = 0, off = 0;
        fits_read_descript(f, gcol, r, &len, &off, &s);
        move[r - 1] = len != 0;
      }
      fits_delete_col(f, gcol, &s);
    } else {
      fits_get_colnum(f, CASEINSEN, "COMPRESSED_DATA", &ccol, &s);
      chk(s, "no COMPRESSED_DATA column");
      for (i = 4; i < argc; i++) {
        r = atol(argv[i]);
        /* a negative row -r moves row r with its last value left out
           (a tile of the wrong length) */
        int shortrow = r < 0;
        if (shortrow) r = -r;
        if (r < 1 || r > nrows) {
          fprintf(stderr, "uncompcol: row %ld out of range\n", r);
          return 1;
        }
        move[r - 1] = shortrow ? 2 : 1;
      }
    }
    fits_get_num_cols(f, &ncol, &s);
    fits_insert_col(f, ncol + 1, "UNCOMPRESSED_DATA", (char *)tform, &s);
    chk(s, "insert UNCOMPRESSED_DATA");
    for (r = 1; r <= nrows; r++) {
      if (!move[r - 1]) continue;
      fits_read_img(in, dtype, (r - 1) * nx + 1, nx, NULL, row, &anyn, &s);
      fits_write_col(f, dtype, ncol + 1, r, 1, move[r - 1] == 2 ? nx - 1 : nx, row, &s);
      if (ccol) fits_write_descript(f, ccol, r, 0, 0, &s);
    }
    free(move);
    free(row);
  }
  fits_close_file(f, &s);
  fits_close_file(in, &s);
  chk(s, "uncompcol");
  return 0;
}

static int do_rewidth(int argc, char **argv) {
  fitsfile *in, *f;
  int s = 0, ccol = 0, anyn = 0, zbitpix = 0, w, extra = 0, gz2 = 0, rice = 0;
  long nrows = 0, r, nx = 0, ztile2 = 0, i, bump = 0, blocksize = 32;
  const char *setlist = NULL;
  char cmd[2048], ctype[FLEN_VALUE];
  if (argc < 5) {
    fprintf(stderr, "usage: gen_fpack rewidth in.fz out.fits W [extra] [bump] [i:v,...]\n");
    return 2;
  }
  w = atoi(argv[4]);
  if (argc >= 6) extra = atoi(argv[5]);
  if (argc >= 7) bump = atol(argv[6]);
  if (argc >= 8) setlist = argv[7];
  if (w != 1 && w != 2 && w != 4 && w != 8) {
    fprintf(stderr, "rewidth: W must be 1, 2, 4 or 8\n");
    return 2;
  }
  snprintf(cmd, sizeof cmd, "cp '%s' '%s'", argv[2], argv[3]);
  if (system(cmd) != 0) return 1;
  fits_open_file(&in, argv[2], READONLY, &s);
  fits_movabs_hdu(in, 2, NULL, &s);
  fits_set_bscale(in, 1.0, 0.0, &s); /* raw stored values */
  fits_open_file(&f, argv[3], READWRITE, &s);
  fits_movabs_hdu(f, 2, NULL, &s);
  fits_read_key(f, TLONG, "ZNAXIS1", &nx, NULL, &s);
  fits_read_key(f, TLONG, "ZTILE2", &ztile2, NULL, &s);
  fits_read_key(f, TINT, "ZBITPIX", &zbitpix, NULL, &s);
  fits_read_key(f, TSTRING, "ZCMPTYPE", ctype, NULL, &s);
  fits_get_colnum(f, CASEINSEN, "COMPRESSED_DATA", &ccol, &s);
  chk(s, "open");
  if (ztile2 != 1 || zbitpix < 0) {
    fprintf(stderr, "rewidth: expects an integer image with row tiles\n");
    return 1;
  }
  gz2 = !strcmp(ctype, "GZIP_2");
  rice = !strcmp(ctype, "RICE_1");
  if (!rice && !gz2 && strcmp(ctype, "GZIP_1")) {
    fprintf(stderr, "rewidth: ZCMPTYPE %s\n", ctype);
    return 1;
  }
  if (rice) {
    int ts = 0;
    if (w == 8 || extra) {
      fprintf(stderr, "rewidth: RICE_1 takes W 1, 2 or 4 and no extra\n");
      return 1;
    }
    fits_read_key(f, TLONG, "ZVAL1", &blocksize, NULL, &ts);
    chk(ts, "ZVAL1");
  }
  fits_get_num_rows(f, &nrows, &s);
  {
    long long *v = calloc(nx, sizeof(long long));
    double *dv = calloc(nx, sizeof(double)); /* exact for 32-bit stored values */
    size_t rawlen = (size_t)nx * w + extra;
    unsigned char *raw = calloc(rawlen, 1), *sh = calloc(rawlen, 1);
    unsigned char *rc = calloc((size_t)nx * 16 + 64, 1);
    for (r = 1; r <= nrows; r++) {
      unsigned char *src = raw, *out = NULL;
      size_t outlen = 0;
      if (zbitpix == 64) {
        /* 64-bit stored values: TLONG (a 64-bit long on LP64; cfitsio's
           TDOUBLE read of 8-byte tiles returns zeros) */
        long *lv = (long *)v;
        fits_read_img(in, TLONG, (r - 1) * nx + 1, nx, NULL, lv, &anyn, &s);
      } else {
        fits_read_img(in, TDOUBLE, (r - 1) * nx + 1, nx, NULL, dv, &anyn, &s);
        for (i = 0; i < nx; i++) v[i] = (long long)dv[i];
      }
      chk(s, "read row");
      if (r == 1) v[0] += bump;
      if (setlist) { /* "i:v,R@i:v": value i of tile 1 (or of tile R) set to v */
        const char *p = setlist;
        while (*p) {
          char *e;
          long row = 1, k = strtol(p, &e, 10);
          if (*e == '@') {
            row = k;
            k = strtol(e + 1, &e, 10);
          }
          if (*e != ':' || k < 0 || k >= nx || row < 1 || row > nrows) {
            fprintf(stderr, "rewidth: bad set list %s\n", setlist);
            return 1;
          }
          long long val = strtoll(e + 1, &e, 10);
          if (row == r) v[k] = val;
          p = *e == ',' ? e + 1 : e;
        }
      }
      for (i = 0; i < nx; i++) {
        long long x = v[i];
        int b;
        if ((w == 1 && (x < 0 || x > 255)) || (w == 2 && (x < -32768 || x > 32767)) ||
            (w == 4 && (x < -2147483648LL || x > 2147483647LL))) {
          fprintf(stderr, "rewidth: value %lld does not fit %d bytes\n", x, w);
          return 1;
        }
        for (b = 0; b < w; b++) raw[i * w + b] = (unsigned char)(x >> (8 * (w - 1 - b)));
      }
      if (rice) {
        int n;
        if (w == 4) {
          int *a = malloc(nx * sizeof(int));
          for (i = 0; i < nx; i++) a[i] = (int)v[i];
          n = fits_rcomp(a, (int)nx, rc, (int)(nx * 16 + 64), (int)blocksize);
          free(a);
        } else if (w == 2) {
          short *a = malloc(nx * sizeof(short));
          for (i = 0; i < nx; i++) a[i] = (short)v[i];
          n = fits_rcomp_short(a, (int)nx, rc, (int)(nx * 16 + 64), (int)blocksize);
          free(a);
        } else {
          signed char *a = malloc(nx);
          for (i = 0; i < nx; i++) a[i] = (signed char)(unsigned char)v[i];
          n = fits_rcomp_byte(a, (int)nx, rc, (int)(nx * 16 + 64), (int)blocksize);
          free(a);
        }
        if (n <= 0) {
          fprintf(stderr, "rewidth: fits_rcomp failed\n");
          return 1;
        }
        out = rc;
        outlen = (size_t)n;
      } else {
        size_t bs = rawlen + 1024, fs = 0;
        if (gz2 && w > 1) { /* cfitsio fits_shuffle_Nbytes: byte b of every value, then b+1 */
          int b;
          for (i = 0; i < nx; i++)
            for (b = 0; b < w; b++) sh[(size_t)b * nx + i] = raw[i * w + b];
          memcpy(sh + (size_t)nx * w, raw + (size_t)nx * w, extra);
          src = sh;
        }
        out = malloc(bs);
        compress2mem_from_mem((char *)src, rawlen, (char **)&out, &bs, realloc, &fs, &s);
        chk(s, "compress2mem_from_mem");
        outlen = fs;
      }
      fits_write_descript(f, ccol, r, 0, 0, &s);
      fits_write_col(f, TBYTE, ccol, r, 1, (long)outlen, out, &s);
      chk(s, "write tile");
      if (out != rc) free(out);
    }
    free(v);
    free(dv);
    free(raw);
    free(sh);
    free(rc);
  }
  if (rice) {
    long lw = w;
    fits_update_key(f, TLONG, "ZVAL2", &lw, "gen_fpack rewidth: BYTEPIX", &s);
  }
  fits_compress_heap(f, &s);
  fits_close_file(f, &s);
  fits_close_file(in, &s);
  chk(s, "rewidth");
  return 0;
}

/* emptytile in out row: COMPRESSED_DATA of tile <row> (1-based) set to
   length 0, with no fallback column to take its place. */
static int do_emptytile(int argc, char **argv) {
  fitsfile *f;
  int s = 0, ccol = 0;
  char cmd[2048];
  if (argc < 5) {
    fprintf(stderr, "usage: gen_fpack emptytile in out row\n");
    return 2;
  }
  snprintf(cmd, sizeof cmd, "cp '%s' '%s'", argv[2], argv[3]);
  if (system(cmd) != 0) return 1;
  fits_open_file(&f, argv[3], READWRITE, &s);
  fits_movabs_hdu(f, 2, NULL, &s);
  fits_get_colnum(f, CASEINSEN, "COMPRESSED_DATA", &ccol, &s);
  fits_write_descript(f, ccol, atol(argv[4]), 0, 0, &s);
  fits_compress_heap(f, &s);
  fits_close_file(f, &s);
  chk(s, "emptytile");
  return 0;
}

/* gzipcol in.fz out row: move integer tile <row> (1-based) of a row-tiled
   image out of COMPRESSED_DATA into a GZIP_COMPRESSED_DATA column (1PB)
   holding its stored values, big-endian in the ZBITPIX width, gzipped by
   libcfitsio (compress2mem_from_mem). cfitsio writes that column only for
   float tiles and refuses it on an integer image. */
static int do_gzipcol(int argc, char **argv) {
  fitsfile *in, *f;
  int s = 0, ccol = 0, ncol = 0, anyn = 0, zbitpix = 0, w, b;
  long nx = 0, r, i;
  char cmd[2048];
  if (argc < 5) {
    fprintf(stderr, "usage: gen_fpack gzipcol in.fz out row\n");
    return 2;
  }
  r = atol(argv[4]);
  snprintf(cmd, sizeof cmd, "cp '%s' '%s'", argv[2], argv[3]);
  if (system(cmd) != 0) return 1;
  fits_open_file(&in, argv[2], READONLY, &s);
  fits_movabs_hdu(in, 2, NULL, &s);
  fits_set_bscale(in, 1.0, 0.0, &s);
  fits_open_file(&f, argv[3], READWRITE, &s);
  fits_movabs_hdu(f, 2, NULL, &s);
  fits_read_key(f, TLONG, "ZNAXIS1", &nx, NULL, &s);
  fits_read_key(f, TINT, "ZBITPIX", &zbitpix, NULL, &s);
  fits_get_colnum(f, CASEINSEN, "COMPRESSED_DATA", &ccol, &s);
  chk(s, "open");
  if (zbitpix <= 0) {
    fprintf(stderr, "gzipcol: expects an integer image\n");
    return 1;
  }
  w = zbitpix / 8;
  {
    double *dv = calloc(nx, sizeof(double));
    unsigned char *raw = calloc((size_t)nx * w, 1);
    char *out = malloc((size_t)nx * w + 1024);
    size_t bs = (size_t)nx * w + 1024, fs = 0;
    fits_read_img(in, TDOUBLE, (r - 1) * nx + 1, nx, NULL, dv, &anyn, &s);
    for (i = 0; i < nx; i++) {
      long long x = (long long)dv[i];
      for (b = 0; b < w; b++) raw[i * w + b] = (unsigned char)(x >> (8 * (w - 1 - b)));
    }
    compress2mem_from_mem((char *)raw, (size_t)nx * w, &out, &bs, realloc, &fs, &s);
    fits_get_num_cols(f, &ncol, &s);
    fits_insert_col(f, ncol + 1, "GZIP_COMPRESSED_DATA", "1PB", &s);
    fits_write_col(f, TBYTE, ncol + 1, r, 1, (long)fs, out, &s);
    fits_write_descript(f, ccol, r, 0, 0, &s);
    free(dv);
    free(raw);
    free(out);
  }
  fits_compress_heap(f, &s);
  fits_close_file(f, &s);
  fits_close_file(in, &s);
  chk(s, "gzipcol");
  return 0;
}

/* rewidthq in.fz out W: re-encode every RICE_1 tile of a quantized float
   image (BYTEPIX 4, row tiles) with W-byte (1 or 2) Rice coding, using
   libcfitsio's own coders: fits_rdecomp gives the stored quantized
   integers, fits_rcomp_byte / fits_rcomp_short write them again, and
   ZVAL2 becomes W. Tiles in a fallback column are left as they are.
   cfitsio reads such tiles with unquantize_i1r4 / unquantize_i2r4. */
static int do_rewidthq(int argc, char **argv) {
  fitsfile *f;
  int s = 0, ccol = 0, w, zbitpix = 0, gz = 0, gz2 = 0;
  long nrows = 0, r, i, nx = 0, blocksize = 32, ztile2 = 0;
  char cmd[2048], ctype[FLEN_VALUE];
  if (argc < 5) {
    fprintf(stderr, "usage: gen_fpack rewidthq in.fz out W\n");
    return 2;
  }
  w = atoi(argv[4]);
  snprintf(cmd, sizeof cmd, "cp '%s' '%s'", argv[2], argv[3]);
  if (system(cmd) != 0) return 1;
  fits_open_file(&f, argv[3], READWRITE, &s);
  fits_movabs_hdu(f, 2, NULL, &s);
  fits_read_key(f, TLONG, "ZNAXIS1", &nx, NULL, &s);
  fits_read_key(f, TLONG, "ZTILE2", &ztile2, NULL, &s);
  fits_read_key(f, TINT, "ZBITPIX", &zbitpix, NULL, &s);
  fits_read_key(f, TSTRING, "ZCMPTYPE", ctype, NULL, &s);
  fits_get_colnum(f, CASEINSEN, "COMPRESSED_DATA", &ccol, &s);
  fits_get_num_rows(f, &nrows, &s);
  chk(s, "open");
  gz2 = !strcmp(ctype, "GZIP_2");
  gz = gz2 || !strcmp(ctype, "GZIP_1");
  if (!gz) {
    int ts = 0;
    fits_read_key(f, TLONG, "ZVAL1", &blocksize, NULL, &ts);
    chk(ts, "ZVAL1");
  }
  if (zbitpix > 0 || ztile2 != 1 || (!gz && strcmp(ctype, "RICE_1"))) {
    fprintf(stderr, "rewidthq: expects a RICE_1/GZIP_1/GZIP_2 float image with row tiles\n");
    return 1;
  }
  if (w != 1 && w != 2 && !(gz && w == 8)) {
    fprintf(stderr, "rewidthq: W must be 1 or 2 (or 8 for GZIP)\n");
    return 2;
  }
  {
    unsigned char *cb = malloc((size_t)nx * 16 + 64), *rc = malloc((size_t)nx * 16 + 64);
    unsigned char *raw = malloc((size_t)nx * 8), *sh = malloc((size_t)nx * 8);
    unsigned int *q = malloc(nx * sizeof(int));
    short *q2 = malloc(nx * sizeof(short));
    signed char *q1 = malloc(nx);
    for (r = 1; r <= nrows; r++) {
      long len = 0, off = 0;
      int n = 0, b;
      unsigned char *out = rc;
      fits_read_descript(f, ccol, r, &len, &off, &s);
      chk(s, "descript");
      if (len == 0) continue;
      fits_read_col(f, TBYTE, ccol, r, 1, len, NULL, cb, NULL, &s);
      chk(s, "read tile");
      if (gz) {
        /* inflate (cfitsio's uncompress2mem_from_mem), un-shuffle GZIP_2,
           big-endian 4-byte quantized integers */
        size_t isz = (size_t)nx * 4, fs = 0;
        unsigned char *in = malloc(isz), *u = in;
        uncompress2mem_from_mem((char *)cb, (size_t)len, (char **)&in, &isz, realloc, &fs, &s);
        chk(s, "uncompress2mem_from_mem");
        if (fs != (size_t)nx * 4) {
          fprintf(stderr, "rewidthq: row %ld inflates to %zu bytes, not 4 per pixel\n", r, fs);
          return 1;
        }
        u = in;
        if (gz2) {
          for (i = 0; i < nx; i++)
            for (b = 0; b < 4; b++) sh[i * 4 + b] = in[(size_t)b * nx + i];
          u = sh;
        }
        for (i = 0; i < nx; i++)
          q[i] = (unsigned int)u[i * 4] << 24 | (unsigned int)u[i * 4 + 1] << 16 | (unsigned int)u[i * 4 + 2] << 8 | u[i * 4 + 3];
        free(in);
      } else if (fits_rdecomp(cb, (int)len, q, (int)nx, (int)blocksize)) {
        fprintf(stderr, "rewidthq: fits_rdecomp failed\n");
        return 1;
      }
      for (i = 0; i < nx; i++) {
        int v = (int)q[i];
        if ((w == 1 && (v < 0 || v > 255)) || (w == 2 && (v < -32768 || v > 32767))) {
          fprintf(stderr, "rewidthq: row %ld value %d does not fit %d bytes\n", r, v, w);
          return 1;
        }
        q2[i] = (short)v;
        q1[i] = (signed char)(unsigned char)v;
        for (b = 0; b < w; b++) raw[i * w + b] = (unsigned char)((long long)v >> (8 * (w - 1 - b)));
      }
      if (gz) {
        size_t bs = (size_t)nx * 8 + 1024, fs = 0;
        unsigned char *src = raw;
        if (gz2 && w > 1) { /* fits_shuffle_Nbytes */
          for (i = 0; i < nx; i++)
            for (b = 0; b < w; b++) sh[(size_t)b * nx + i] = raw[i * w + b];
          src = sh;
        }
        out = malloc(bs);
        compress2mem_from_mem((char *)src, (size_t)nx * w, (char **)&out, &bs, realloc, &fs, &s);
        chk(s, "compress2mem_from_mem");
        n = (int)fs;
      } else {
        n = w == 2 ? fits_rcomp_short(q2, (int)nx, rc, (int)(nx * 16 + 64), (int)blocksize)
                   : fits_rcomp_byte(q1, (int)nx, rc, (int)(nx * 16 + 64), (int)blocksize);
        if (n <= 0) {
          fprintf(stderr, "rewidthq: fits_rcomp failed\n");
          return 1;
        }
      }
      fits_write_descript(f, ccol, r, 0, 0, &s);
      fits_write_col(f, TBYTE, ccol, r, 1, n, out, &s);
      chk(s, "write tile");
      if (out != rc) free(out);
    }
    free(cb);
    free(rc);
    free(raw);
    free(sh);
    free(q);
    free(q2);
    free(q1);
  }
  if (!gz) {
    long lw = w;
    fits_update_key(f, TLONG, "ZVAL2", &lw, "gen_fpack rewidthq: BYTEPIX", &s);
  }
  fits_compress_heap(f, &s);
  fits_close_file(f, &s);
  chk(s, "rewidthq");
  return 0;
}

/* mk64 out.fits ALGO [V]: a ZBITPIX 64 tile-compressed image (12x8, two
   12x4 tiles) built by hand from libcfitsio's own coders, which fpack and
   fits_img_compress will not do for 64-bit images: RICE_1 (fits_rcomp,
   BYTEPIX 4, BLOCKSIZE 32), HCOMPRESS_1 (fits_hcompress64, SCALE 0,
   lossless) or PLIO_1 (pl_p2li, a 1PI column). cfitsio's reader decodes
   these into C ints (imcompress.c imcomp_decompress_tile: RICE BYTEPIX 4,
   fits_hdecompress64 then an (int) cast, pl_l2pi) and widens them to the
   64-bit image, so the values stay in the int32 range: both int32 extremes
   (RICE_1, HCOMPRESS_1) or 0 .. 2^24 - 1 (PLIO_1). With V pixel 5 is V
   (HCOMPRESS_1 with one past either int32 end, which cfitsio's (int) cast
   wraps). */
static int do_mk64(int argc, char **argv) {
  enum { MX = 12, MY = 8, TY = 4 };
  fitsfile *f;
  int s = 0, plio, hc, rice, big;
  long long bigv = 0;
  long t, i;
  char path[4096];
  long long v[MX * MY];
  if (argc < 4) {
    fprintf(stderr, "usage: gen_fpack mk64 out.fits RICE_1|HCOMPRESS_1|PLIO_1 [V]\n");
    return 2;
  }
  rice = !strcmp(argv[3], "RICE_1");
  hc = !strcmp(argv[3], "HCOMPRESS_1");
  plio = !strcmp(argv[3], "PLIO_1");
  big = argc >= 5;
  if (big) bigv = strtoll(argv[4], NULL, 10);
  if (!rice && !hc && !plio) {
    fprintf(stderr, "mk64: ZCMPTYPE %s\n", argv[3]);
    return 2;
  }
  for (i = 0; i < MX * MY; i++) {
    if (plio) {
      v[i] = (i / 3 % 4 == 0) ? 0 : (long long)((i / 3) * 524287L % 16777216L);
    } else {
      v[i] = (long long)((i * 7919L) % 20001L) - 10000 + (i / MX) * 100000LL;
    }
  }
  if (!plio) {
    v[1] = -2147483648LL;
    v[2] = 2147483647LL;
    v[MX * TY + 3] = -2147483647LL;
  }
  if (big) v[5] = bigv;
  snprintf(path, sizeof path, "!%s", argv[2]);
  fits_create_file(&f, path, &s);
  fits_create_img(f, 8, 0, NULL, &s);
  {
    char *tt[] = {"COMPRESSED_DATA"}, *tf[] = {plio ? "1PI" : "1PB"};
    fits_create_tbl(f, BINARY_TBL, MY / TY, 1, tt, tf, NULL, "COMPRESSED_IMAGE", &s);
  }
  {
    int yes = 1, zbp = 64, zn = 2;
    long n1 = MX, n2 = MY, t1 = MX, t2 = TY;
    fits_write_key(f, TLOGICAL, "ZIMAGE", &yes, "extension contains compressed image", &s);
    fits_write_key(f, TLONG, "ZTILE1", &t1, "size of tiles to be compressed", &s);
    fits_write_key(f, TLONG, "ZTILE2", &t2, "size of tiles to be compressed", &s);
    fits_write_key(f, TSTRING, "ZCMPTYPE", argv[3], "compression algorithm", &s);
    if (rice) {
      long bs = 32, bp = 4;
      fits_write_key(f, TSTRING, "ZNAME1", "BLOCKSIZE", "compression block size", &s);
      fits_write_key(f, TLONG, "ZVAL1", &bs, "pixels per block", &s);
      fits_write_key(f, TSTRING, "ZNAME2", "BYTEPIX", "bytes per pixel (1, 2, 4, or 8)", &s);
      fits_write_key(f, TLONG, "ZVAL2", &bp, "bytes per pixel (1, 2, 4, or 8)", &s);
    } else if (hc) {
      long sc = 0, sm = 0;
      fits_write_key(f, TSTRING, "ZNAME1", "SCALE", "HCOMPRESS scale factor", &s);
      fits_write_key(f, TLONG, "ZVAL1", &sc, "HCOMPRESS scale factor", &s);
      fits_write_key(f, TSTRING, "ZNAME2", "SMOOTH", "HCOMPRESS smooth option", &s);
      fits_write_key(f, TLONG, "ZVAL2", &sm, "HCOMPRESS smooth option", &s);
    }
    fits_write_key(f, TINT, "ZBITPIX", &zbp, "number of bits per data pixel", &s);
    fits_write_key(f, TINT, "ZNAXIS", &zn, "number of data axes", &s);
    fits_write_key(f, TLONG, "ZNAXIS1", &n1, "length of data axis 1", &s);
    fits_write_key(f, TLONG, "ZNAXIS2", &n2, "length of data axis 2", &s);
  }
  chk(s, "mk64 header");
  for (t = 0; t < MY / TY; t++) {
    long long *tv = v + t * MX * TY;
    long n = MX * TY;
    if (rice) {
      int a[MX * TY];
      unsigned char c[MX * TY * 16 + 64];
      int len;
      for (i = 0; i < n; i++) a[i] = (int)tv[i];
      len = fits_rcomp(a, (int)n, c, (int)sizeof c, 32);
      if (len <= 0) {
        fprintf(stderr, "mk64: fits_rcomp failed\n");
        return 1;
      }
      fits_write_col(f, TBYTE, 1, t + 1, 1, len, c, &s);
    } else if (hc) {
      LONGLONG a[MX * TY];
      char c[MX * TY * 16 + 1024];
      long len = sizeof c;
      for (i = 0; i < n; i++) a[i] = tv[i];
      /* nx, ny as imcomp_compress_tile passes them: tile axis 1, axis 2 */
      fits_hcompress64(a, TY, MX, 0, c, &len, &s);
      chk(s, "fits_hcompress64");
      fits_write_col(f, TBYTE, 1, t + 1, 1, len, c, &s);
    } else {
      int a[MX * TY];
      short c[MX * TY * 4 + 64];
      int len;
      for (i = 0; i < n; i++) a[i] = (int)tv[i];
      len = pl_p2li(a, 1, c, (int)n);
      fits_write_col(f, TSHORT, 1, t + 1, 1, len, c, &s);
    }
    chk(s, "mk64 tile");
  }
  fits_close_file(f, &s);
  chk(s, "mk64");
  return 0;
}

/* croprows in out N: keep the first N tiles (table rows) of a row-tiled
   compressed image (ZTILE2 = 1) byte for byte: the later rows are deleted
   (fits_delete_rows), ZNAXIS2 becomes N, the heap is compacted
   (fits_compress_heap) and the checksum keywords, which no longer match,
   are removed (CHECKSUM/DATASUM of the HDU, ZDATASUM of the uncompressed
   image). The kept tiles' compressed bytes are unchanged. */
static int do_croprows(int argc, char **argv) {
  fitsfile *f;
  int s = 0, ts;
  long nrows = 0, ztile2 = 0, n;
  char cmd[2048];
  const char *drop[] = {"CHECKSUM", "DATASUM", "ZDATASUM"};
  size_t k;
  if (argc < 5) {
    fprintf(stderr, "usage: gen_fpack croprows in out N\n");
    return 2;
  }
  n = atol(argv[4]);
  snprintf(cmd, sizeof cmd, "cp '%s' '%s'", argv[2], argv[3]);
  if (system(cmd) != 0) return 1;
  fits_open_file(&f, argv[3], READWRITE, &s);
  fits_movabs_hdu(f, 2, NULL, &s);
  fits_read_key(f, TLONG, "ZTILE2", &ztile2, NULL, &s);
  fits_get_num_rows(f, &nrows, &s);
  chk(s, "open");
  if (ztile2 != 1 || n < 1 || n >= nrows) {
    fprintf(stderr, "croprows: expects row tiles and 1 <= N < %ld\n", nrows);
    return 1;
  }
  fits_update_key(f, TLONG, "ZNAXIS2", &n, NULL, &s);
  chk(s, "ZNAXIS2");
  fits_delete_rows(f, n + 1, nrows - n, &s);
  chk(s, "delete rows");
  fits_compress_heap(f, &s);
  chk(s, "compress heap");
  for (k = 0; k < sizeof drop / sizeof drop[0]; k++) {
    ts = 0;
    fits_delete_key(f, drop[k], &ts);
  }
  chk(s, "croprows");
  fits_movabs_hdu(f, 1, NULL, &s);
  for (k = 0; k < 2; k++) {
    ts = 0;
    fits_delete_key(f, drop[k], &ts);
  }
  fits_close_file(f, &s);
  chk(s, "croprows close");
  return 0;
}

int main(int argc, char **argv) {
  if (argc >= 3 && !strcmp(argv[1], "src")) return do_src(argv[2]);
  if (argc >= 2 && !strcmp(argv[1], "hcomp")) return do_hcomp(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "zblankcol")) return do_zblankcol(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "editkey")) return do_editkey(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "setkey")) return do_setkey(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "uncompcol")) return do_uncompcol(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "patch")) return do_patch(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "rewidth")) return do_rewidth(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "rewidthq")) return do_rewidthq(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "emptytile")) return do_emptytile(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "gzipcol")) return do_gzipcol(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "croprows")) return do_croprows(argc, argv);
  if (argc >= 2 && !strcmp(argv[1], "mk64")) return do_mk64(argc, argv);
  fprintf(stderr, "usage: gen_fpack src <dir> | hcomp ... | zblankcol ...\n");
  return 2;
}
