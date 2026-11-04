package fits

import (
	"bytes"
	"compress/gzip"
	"compress/zlib"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strings"
	"sync"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// The tile-compression READ corpus: testdata/comp_fz_<case>.fits were
// written by cfitsio (fpack, or fits_img_compress for what fpack cannot
// make), testdata/comp_fzun_<case>.fits are their funpack twins (a twin
// byte-identical to an earlier one is committed once: testdata/
// comp_fzun_same.txt maps it, see twinFile), and
// testdata/cfref/golden/<file>.json is what libcfitsio itself decodes from
// each file. testdata/gen_fpack.sh regenerates all of it byte-for-byte;
// tool versions are in testdata/COMP_FZ_PROVENANCE.md.
//
// Every case is compared BIT-EXACT against the cfitsio decode of the
// compressed file. For BITPIX=-32 the logical pixel is a float32 (what
// funpack writes and what an uncompressed image would hold), so the truth
// is float32(cfitsio TDOUBLE value); the test also checks that this equals
// the funpack twin, except where funpack 1.7.0 (CFITSIO 4.7.0) is known to be wrong
// (twinWrong), in which case the source image is the arbiter.

type fzTwin int

const (
	twinExact     fzTwin = iota // funpack twin == cfitsio decode of the .fz
	twinWrong                   // funpack 1.7.0 (CFITSIO 4.7.0) output is wrong; see comment on the case
	twinNoNull                  // twin cannot carry per-tile nulls (single BLANK); values match off-null
	twinRescaled                // twin stores float32 before BSCALE/BZERO: within 1 float32 ULP
	twinBothValid               // twin carries BLANK, not the ZBLANK that wins: compare where both are valid
)

type fzCase struct {
	name    string
	cmp     string            // ZCMPTYPE
	zbitpix int               // ZBITPIX
	cards   map[string]string // header facts the case exists to exercise ("" = present, any value)
	cols    []string          // binary-table columns that must exist
	twin    fzTwin
	src     string // comp_fzsrc_<src>.fits: NaN positions must match the source
	exact   bool   // with src: decoded values must equal the source bit-for-bit (lossless)
	// with src: these image rows (NAXIS2 index, 0-based) must equal the
	// source bit-for-bit (tiles stored losslessly in GZIP_COMPRESSED_DATA).
	exactRows []int
	nullAt    []int // the only null pixel indices (known from the source)
	// twinTruth: libcfitsio's direct decode of the .fz (cfref golden) is
	// wrong for this file, so the funpack twin (fits_img_decompress, an
	// independent cfitsio path, which astropy 8.0.2 also agrees with) is
	// the truth; the test fails if the golden stops disagreeing.
	twinTruth bool
	// rawRows are image rows (0-based) held in a GZIP_COMPRESSED_DATA or
	// UNCOMPRESSED_DATA fallback tile of an image with BSCALE/BZERO or
	// BLANK. libcfitsio's fits_read_img returns those tiles' stored values
	// unscaled and without null checking (imcomp_decompress_tile reads them
	// with fffr4r4(..., 1., 0., ...) or fits_read_col), while it scales the
	// image's other tiles, so one image mixes two units. funpack writes
	// every tile's stored values under the image's BSCALE/BZERO/BLANK, and
	// astropy decodes the same way; the library follows them. In these rows
	// the funpack twin (read by cfitsio as a plain image) is the truth, and
	// the test requires cfitsio to still differ there so the note stays
	// true. The same holds for a NaN in a regular tile of a float image
	// without a ZBLANK: fits_read_img turns null checking off for such a
	// tile (imcomp_decompress_tile: cn_zblank == 0 => nullcheck = 0) and
	// returns the NaN as a valid value, while the funpack twin and the
	// plain-image rule flag it (f32_rice_zqnone row 10).
	rawRows []int
	// scaledF64: a BITPIX -32 image whose stored float32 values are scaled
	// after decoding (unquantized, with ZSCALE/ZZERO keywords). Read as
	// float64 the scaling is done in double precision, exactly as cfitsio's
	// TDOUBLE read and ReadPixels on an uncompressed -32 image with
	// BSCALE; that read is compared with the unrounded cfitsio values, the
	// float32 read with them rounded to float32.
	scaledF64 bool
	// noCards: header cards that must be ABSENT (e.g. ZVAL2, whose absence
	// means RICE_1 BYTEPIX 4).
	noCards []string
	// blankRefused: the BLANK card lies outside the ZBITPIX range, which the
	// scaling-card rule refuses (scalingCards): cfitsio compares it cast to
	// the tile type, astropy refuses the file (RICE_1: "out of bounds for
	// int16", BLANK 2^32: "too large") or compares the integer (HCOMPRESS_1,
	// PLIO_1: no nulls). Every pixel read and BSCALE/BZERO return
	// ErrInvalidScaling (checkCorpusScalingRefused); the structure stays
	// readable. The cfitsio goldens record what cfitsio does.
	blankRefused bool
	// tileWidth: GZIP_1/GZIP_2 tile 1 inflates (Go's compress/gzip or
	// compress/zlib, not the library decoder) to tileWidth bytes per pixel
	// of the tile; 0 = the ZBITPIX width, not checked. For 8-byte tiles the
	// test also pins libcfitsio's inconsistent reads: its TDOUBLE
	// fits_read_img (the golden's phys) returns memory it never wrote
	// (unwrittenDouble), its TLONG read (the golden's i64) the right values.
	tileWidth int
	// i64Truth: the golden's i64 (libcfitsio's TLONG fits_read_img) is the
	// truth: for 8-byte tiles on an 8-bit image cfitsio's TDOUBLE read
	// returns unwritten memory and funpack zeros (imcomp_decompress_tile
	// has no TDOUBLE-tile branch for a TBYTE or TDOUBLE output).
	i64Truth bool
	// sameAs: the truth must also equal the funpack twin of that case (the
	// same source image), an independent check of the chosen truth.
	sameAs string
	// at: truth values the case exists for (pixel index -> physical
	// value), checked so a generator change cannot drop them.
	at map[int]float64
	// zeros: the truth has valid zero pixels (what a null value cast to
	// the 16-bit image type would flag).
	zeros bool
	// clipTruth: lossy HCOMPRESS_1 of a 16-bit image reconstructs some
	// pixels outside the ZBITPIX range. cfitsio's TSHORT read clips them
	// (imcomp_decompress_tile resets the overflow status for HCOMPRESS_1
	// after fffi4i2; its TUSHORT read of a uint16 image fails with 412)
	// and so does funpack, whose twin is the truth; cfitsio's TDOUBLE read
	// (the golden) returns them unclipped.
	// The test requires the golden to differ from the twin exactly at
	// pixels whose TDOUBLE value is outside the physical range of the
	// image type, with the twin holding that value clipped to the range.
	clipTruth bool
	// minNulls: the truth has at least this many null pixels.
	minNulls int
}

// twinFile is the funpack twin of a corpus case: comp_fzun_<case>.fits,
// or, when that twin was byte-identical to an earlier one, the file
// testdata/comp_fzun_same.txt (written by gen_fpack.sh) names for it. The
// twins are not committed (testdata/gen_twins.sh regenerates them); their
// cfitsio goldens are, in testdata/cfref/golden/twins.
func twinFile(t *testing.T, name string) string {
	t.Helper()
	twinSameOnce.Do(func() {
		b, err := os.ReadFile(filepath.Join("testdata", "comp_fzun_same.txt"))
		if err != nil {
			twinSameErr = err
			return
		}
		twinSame = map[string]string{}
		for _, ln := range strings.Split(strings.TrimSpace(string(b)), "\n") {
			f := strings.Fields(ln)
			if len(f) != 2 {
				twinSameErr = fmt.Errorf("comp_fzun_same.txt: bad line %q", ln)
				return
			}
			twinSame[f[0]] = f[1]
		}
	})
	if twinSameErr != nil {
		t.Fatalf("committed twin map missing (failure, not a skip): %v", twinSameErr)
	}
	f := "comp_fzun_" + name + ".fits"
	if k, ok := twinSame[f]; ok {
		return k
	}
	return f
}

var (
	twinSameOnce sync.Once
	twinSame     map[string]string
	twinSameErr  error
)

// unwrittenDouble is the float64 whose bytes are all 0xA5: what the cfref
// harness reads where libcfitsio returns memory it never wrote (cfref.c
// runs with glibc's malloc perturbation, which fills every malloc'd block
// with 0xA5). libcfitsio 4.7.0's TDOUBLE fits_read_img of an integer image
// with 8-byte GZIP tiles does that: imcomp_decompress_tile has no branch
// for such a tile and copies its tile buffer out unwritten, with status 0
// (without the fill the values are whatever the heap held: zeros or heap
// addresses).
const unwrittenDouble = 0xa5a5a5a5a5a5a5a5

// cfrefI64 is the golden's i64 view of the image HDU (libcfitsio's TLONG
// fits_read_img), as float64 values and null flags.
func cfrefI64(t *testing.T, name string) ([]float64, []bool) {
	t.Helper()
	g := testutil.LoadCfrefGolden(t, name)
	for i := range g.HDUs {
		h := &g.HDUs[i]
		if h.Image == nil || h.Image.Npix == 0 {
			continue
		}
		v := h.Image.I64.Values(t)
		vals, null := make([]float64, len(v)), make([]bool, len(v))
		for k, e := range v {
			vals[k], null[k] = float64(e.V), !e.Valid
			if !e.Valid {
				vals[k] = math.NaN()
			}
		}
		return vals, null
	}
	t.Fatalf("%s: golden has no image pixels", name)
	return nil, nil
}

// caseTruth is the pixel truth of corpus case c as the other corpus tests
// use it: the cfitsio decode of the .fz, or the funpack twin (twinTruth),
// or cfitsio's TLONG read (i64Truth).
func caseTruth(t *testing.T, c fzCase, f32 bool) ([]float64, []bool) {
	t.Helper()
	switch {
	case c.twinTruth, c.clipTruth:
		_, v, n := cfrefPixels(t, twinFile(t, c.name), f32)
		return v, n
	case c.i64Truth:
		return cfrefI64(t, "comp_fz_"+c.name+".fits")
	}
	_, v, n := cfrefPixels(t, "comp_fz_"+c.name+".fits", f32)
	return v, n
}

var fzCorpus = []fzCase{
	// Integer BITPIX; BZERO/BSCALE live under their plain names in the
	// compressed header (cfitsio imcompress.c fits_img_compress / imcomp_get_compressed_image_par).
	{name: "u8_rice", cmp: "RICE_1", zbitpix: 8, cards: map[string]string{"ZVAL2": "1"}},
	{name: "u8_gzip1", cmp: "GZIP_1", zbitpix: 8},
	{name: "u8_hcomp", cmp: "HCOMPRESS_1", zbitpix: 8},
	{name: "s8_rice", cmp: "RICE_1", zbitpix: 8, cards: map[string]string{"BZERO": "-128"}},
	{name: "s8_gzip2", cmp: "GZIP_2", zbitpix: 8, cards: map[string]string{"BZERO": "-128"}},
	{name: "u16_rice", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"BZERO": "32768"}},
	{name: "u16_gzip2", cmp: "GZIP_2", zbitpix: 16, cards: map[string]string{"BZERO": "32768"}},
	{name: "u16_hcomp", cmp: "HCOMPRESS_1", zbitpix: 16, cards: map[string]string{"BZERO": "32768"}},
	{name: "u32_rice", cmp: "RICE_1", zbitpix: 32, cards: map[string]string{"BZERO": "2147483648"}},
	{name: "u32_hcomp", cmp: "HCOMPRESS_1", zbitpix: 32, cards: map[string]string{"BZERO": "2147483648"}},
	// Lossy HCOMPRESS on full-range uint32: undigitize (coefficient times
	// scale) overflows 32 bits, so it must run in 64 bits as cfitsio's
	// fits_hdecompress64 does.
	{name: "u32_hcomp_s4", cmp: "HCOMPRESS_1", zbitpix: 32, cards: map[string]string{"BZERO": "2147483648", "ZVAL1": "4.0000000E+00"}},
	{name: "u32_hcomp_s1024", cmp: "HCOMPRESS_1", zbitpix: 32, cards: map[string]string{"BZERO": "2147483648", "ZVAL1": "1.0240000E+03"}},
	{name: "i32_hcomp", cmp: "HCOMPRESS_1", zbitpix: 32},
	{name: "i32_gzip2", cmp: "GZIP_2", zbitpix: 32},
	{name: "i16scaled_rice", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"BSCALE": "0.5", "BZERO": "100."}},
	{name: "i16blank_rice", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"BLANK": "-32768"}},
	{name: "i16blank_gzip1", cmp: "GZIP_1", zbitpix: 16, cards: map[string]string{"BLANK": "-32768"}},
	{name: "i16mask_plio", cmp: "PLIO_1", zbitpix: 16},
	{name: "i16_hcomp_s4", cmp: "HCOMPRESS_1", zbitpix: 16, cards: map[string]string{"ZVAL1": "4.0000000E+00"}},
	{name: "cube16_rice_t3", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"ZNAXIS": "3", "ZTILE3": "2"}},
	{name: "cube16_gzip2_t3", cmp: "GZIP_2", zbitpix: 16, cards: map[string]string{"ZNAXIS": "3", "ZTILE3": "4"}},
	// HCOMPRESS SMOOTH (ZVAL2 = 1), lossy scale so smoothing changes pixels.
	{name: "i16_hcomp_smooth", cmp: "HCOMPRESS_1", zbitpix: 16, cards: map[string]string{"ZNAME2": "SMOOTH", "ZVAL2": "1"}},
	{name: "i32_hcomp_smooth", cmp: "HCOMPRESS_1", zbitpix: 32, cards: map[string]string{"ZNAME2": "SMOOTH", "ZVAL2": "1"}},
	{name: "f32_hcomp_smooth", cmp: "HCOMPRESS_1", zbitpix: -32, cards: map[string]string{"ZVAL2": "1", "ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}},
	// Quantized float.
	{name: "f32_rice_sd1", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1", "ZDITHER0": "17"}, cols: []string{"ZSCALE", "ZZERO"}},
	{name: "f32_rice_sd2", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_2", "ZDITHER0": "17"}},
	{name: "f32_rice_nodith", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "NO_DITHER"}},
	{name: "f32_rice_abs", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1", "ZDITHER0": "23"}},
	{name: "f32_gzip1_sd1", cmp: "GZIP_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}},
	{name: "f32_gzip2_sd1", cmp: "GZIP_2", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}},
	{name: "f32nan_rice_sd1", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZBLANK": "-2147483647"}, src: "f32nan"},
	// funpack 1.7.0 (CFITSIO 4.7.0) writes wrong pixels for this 20x8-tiled NaN image (840
	// of 960 differ from the library read, NaN positions differ from the
	// source, max error 72.9 vs 0.99 for the library read); the library
	// decode matches the source, so the source is the arbiter.
	{name: "f32nan_rice_t", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZTILE1": "20", "ZTILE2": "8"}, twin: twinWrong, src: "f32nan"},
	{name: "f32_hcomp_q", cmp: "HCOMPRESS_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}},
	{name: "f32nan_hcomp_q", cmp: "HCOMPRESS_1", zbitpix: -32, cards: map[string]string{"ZBLANK": "-2147483647"}, src: "f32nan"},
	// Row 5 is constant: cfitsio cannot quantize it (zero noise) and stores
	// it gzipped, unquantized, in GZIP_COMPRESSED_DATA.
	{name: "f32flat_rice", cmp: "RICE_1", zbitpix: -32, cols: []string{"GZIP_COMPRESSED_DATA"}, src: "f32flat", exactRows: []int{5}},
	{name: "cubef32_rice_t3", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZNAXIS": "3", "ZTILE3": "2"}},
	{name: "f64_rice_sd1", cmp: "RICE_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1", "ZVAL2": "4"}},
	{name: "f64_gzip1_sd1", cmp: "GZIP_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}},
	{name: "f64_gzip2_sd1", cmp: "GZIP_2", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}},
	{name: "f64_hcomp_q", cmp: "HCOMPRESS_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}},
	{name: "f64nan_rice_sd1", cmp: "RICE_1", zbitpix: -64, cards: map[string]string{"ZBLANK": "-2147483647"}, src: "f64nan"},
	// Lossless float: ZQUANTIZ = 'NONE', no ZSCALE/ZZERO columns.
	{name: "f32_gzip1_lossless", cmp: "GZIP_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "NONE"}},
	{name: "f32_gzip2_lossless", cmp: "GZIP_2", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "NONE"}},
	{name: "f32nan_gzip2_lossless", cmp: "GZIP_2", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "NONE"}, src: "f32nan", exact: true},
	{name: "f64_gzip1_lossless", cmp: "GZIP_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "NONE"}},
	{name: "f64_gzip2_lossless", cmp: "GZIP_2", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "NONE"}},
	// Per-tile ZBLANK column (tile r's null value is -32768 + 7*(r-1)). The
	// funpack twin has no way to express per-tile nulls (it drops them), so
	// only its non-null pixels are compared.
	{name: "i16_zblankcol", cmp: "RICE_1", zbitpix: 16, cols: []string{"ZBLANK"}, twin: twinNoNull},
	{name: "f32nan_zblankcol", cmp: "RICE_1", zbitpix: -32, cols: []string{"ZBLANK"}, src: "f32nan"},
	// SUBTRACTIVE_DITHER_1 written with ZDITHER0 = 1 and the keyword then
	// removed (files before Oct 2009 have none): cfitsio defaults it to 1.
	{name: "f32_rice_nozd0", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}},
	// Quantized float with BSCALE 2 / BZERO 10 added to the header; cfitsio
	// combines them with ZSCALE/ZZERO. funpack writes the float32 values
	// before scaling, so its twin rounds differently (<= 1 float32 ULP).
	{name: "f32scaled_rice", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"BSCALE": "2.", "BZERO": "10."}, twin: twinRescaled},
	// u16_rice with ZTILE1/ZTILE2 removed: the convention's defaults
	// (ZTILE1 = ZNAXIS1, others 1) give the same row tiles.
	{name: "u16_rice_noztile", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"BZERO": "32768"}},
	// The GZIP_COMPRESSED_DATA tile of f32flat_rice moved to a pre-2011
	// UNCOMPRESSED_DATA (1PE) column.
	{name: "f32flat_uncomp", cmp: "RICE_1", zbitpix: -32, cols: []string{"UNCOMPRESSED_DATA"}, src: "f32flat", exactRows: []int{5}},
	// BITPIX -64 constant row: GZIP_COMPRESSED_DATA holds raw doubles.
	{name: "f64flat_rice", cmp: "RICE_1", zbitpix: -64, cols: []string{"GZIP_COMPRESSED_DATA"}, src: "f64flat", exactRows: []int{5}},
	// The same tile as pre-2011 UNCOMPRESSED_DATA, 1PD.
	{name: "f64flat_uncomp", cmp: "RICE_1", zbitpix: -64, cols: []string{"UNCOMPRESSED_DATA"}, cards: map[string]string{"TFORM4": "1PD(40)"}, src: "f64flat", exactRows: []int{5}},
	// Fallback tiles in a quantized float image with BSCALE 2 / BZERO 10:
	// row 5 is 7.25 stored, 24.5 physical (see rawRows).
	{name: "f32flat_rice_scaled", cmp: "RICE_1", zbitpix: -32, cols: []string{"GZIP_COMPRESSED_DATA"}, cards: map[string]string{"BSCALE": "2.", "BZERO": "10."}, twin: twinRescaled, rawRows: []int{5}},
	{name: "f32flat_uncomp_scaled", cmp: "RICE_1", zbitpix: -32, cols: []string{"UNCOMPRESSED_DATA"}, cards: map[string]string{"BSCALE": "2.", "BZERO": "10."}, twin: twinRescaled, rawRows: []int{5}},
	// Integer tiles moved to UNCOMPRESSED_DATA (rows 4 and 11, 3, 4, 1 and
	// 4, 1-based): 1PJ, 1PB, 1PI with BZERO 32768, 1PI with BLANK.
	{name: "i32_uncomp", cmp: "RICE_1", zbitpix: 32, cols: []string{"UNCOMPRESSED_DATA"}, cards: map[string]string{"TFORM2": "1PJ(40)"}},
	{name: "u8_uncomp", cmp: "RICE_1", zbitpix: 8, cols: []string{"UNCOMPRESSED_DATA"}, cards: map[string]string{"TFORM2": "1PB(40)"}},
	{name: "u16_uncomp", cmp: "RICE_1", zbitpix: 16, cols: []string{"UNCOMPRESSED_DATA"}, cards: map[string]string{"BZERO": "32768"}, rawRows: []int{3}},
	{name: "i16blank_uncomp", cmp: "RICE_1", zbitpix: 16, cols: []string{"UNCOMPRESSED_DATA"}, cards: map[string]string{"BLANK": "-32768"}, rawRows: []int{0, 3}},
	// i16blank_rice plus ZBLANK = 26: cfitsio takes ZBLANK over BLANK, so
	// the 6 pixels equal to 26 are null and the -32768 pixels are values.
	{name: "i16blank_zblankkw", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"ZBLANK": "26", "BLANK": "-32768"}, twin: twinBothValid},
	// Lossless float (no ZSCALE/ZZERO columns) with ZSCALE 2 / ZZERO
	// keywords: cfitsio applies them to the stored floats.
	{name: "f32_lossless_zscalekw", cmp: "GZIP_1", zbitpix: -32, cards: map[string]string{"ZSCALE": "2.", "ZZERO": "10."}, scaledF64: true},
	// Lossless float with a ZQUANTIZ cfitsio does not know: without
	// ZSCALE/ZZERO columns cfitsio never reads ZQUANTIZ.
	{name: "f32_lossless_futurezq", cmp: "GZIP_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "FUTURE_METHOD"}},
	// ZNAME1/ZNAME2 swapped: cfitsio takes RICE BLOCKSIZE/BYTEPIX and
	// HCOMPRESS SCALE/SMOOTH from ZVAL1/ZVAL2 by position.
	{name: "u16_rice_znameswap", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"ZNAME1": "BYTEPIX", "ZNAME2": "BLOCKSIZE", "ZVAL1": "32", "ZVAL2": "2"}},
	{name: "i16_hcomp_smooth_znameswap", cmp: "HCOMPRESS_1", zbitpix: 16, cards: map[string]string{"ZNAME1": "SMOOTH", "ZNAME2": "SCALE", "ZVAL2": "1"}},
	// Pre-2009 RICE header: ZNAME2 = 'NOISEBIT' means ZVAL2 is not
	// BYTEPIX, and cfitsio uses BYTEPIX 4 (ZVAL2 16 here would be invalid).
	{name: "f32_rice_noisebit", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZNAME2": "NOISEBIT", "ZVAL2": "1", "ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}, sameAs: "f32_rice_sd1"},
	// ZVAL1 = 2 / ZVAL2 = 32 (values in the wrong slots): BLOCKSIZE < 16
	// and BYTEPIX > 8, so cfitsio swaps them back to 32 / 2.
	{name: "u16_rice_zvalswap", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"ZVAL1": "2", "ZVAL2": "32", "BZERO": "32768"}},
	// Real files from astropy's test data (BSD-3-Clause), verbatim.
	// compressed_with_nan.fits: fpack'd 2022 by CFITSIO, BITPIX -64,
	// SUBTRACTIVE_DITHER_1 with ZDITHER0 2961, 6x6 tiles, source
	// arange(144) with NaN at [1,1]. funpack 1.7.0 (CFITSIO 4.7.0) decodes it wrongly (30
	// pixels, NaN lost) like the other tiled NaN case.
	{name: "real_with_nan", cmp: "RICE_1", zbitpix: -64, cards: map[string]string{"ZDITHER0": "2961", "ZTILE1": "6"}, twin: twinWrong, nullAt: []int{13}},
	// compressed_image.fits: 10x10 int16 RICE_1 with ZTENSION.
	{name: "real_image", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"ZTENSION": "IMAGE"}},

	// The remaining cells of the claimed matrix (corpusMatrix below).
	{name: "s8_gzip1", cmp: "GZIP_1", zbitpix: 8, cards: map[string]string{"BZERO": "-128"}},
	{name: "u16_gzip1", cmp: "GZIP_1", zbitpix: 16, cards: map[string]string{"BZERO": "32768"}},
	{name: "i32_gzip1", cmp: "GZIP_1", zbitpix: 32},
	{name: "u32_gzip1", cmp: "GZIP_1", zbitpix: 32, cards: map[string]string{"BZERO": "2147483648"}},
	{name: "u8_gzip2", cmp: "GZIP_2", zbitpix: 8},
	{name: "u32_gzip2", cmp: "GZIP_2", zbitpix: 32, cards: map[string]string{"BZERO": "2147483648"}},
	{name: "s8_hcomp", cmp: "HCOMPRESS_1", zbitpix: 8, cards: map[string]string{"BZERO": "-128"}},
	// PLIO_1 at ZBITPIX 8 (plain and signed-byte BZERO -128) and 32
	// (values up to 2^24-1), and the uint16 convention: cfitsio stores
	// value + 32768 for PLIO and subtracts it on read (imcompress.c:2343,
	// 6477).
	{name: "u8_plio", cmp: "PLIO_1", zbitpix: 8},
	{name: "s8_plio", cmp: "PLIO_1", zbitpix: 8, cards: map[string]string{"BZERO": "-128"}},
	{name: "i32_plio", cmp: "PLIO_1", zbitpix: 32},
	{name: "u16_plio", cmp: "PLIO_1", zbitpix: 16, cards: map[string]string{"BZERO": "32768", "ZTILE2": "6"}},
	// uint16 PLIO with BLANK = -32743 (the stored value of the physical 25
	// that 30 pixels hold): cfitsio compares the null value with the PLIO
	// value (stored + 32768) and flags no pixel. Comparing the stored value
	// instead would null those 30 pixels. funpack writes the stored values
	// under the same BLANK into a plain image, where the plain-image rule
	// (stored value == BLANK) nulls them: cfitsio's two paths disagree, and
	// the library follows fits_read_img (no value is dropped). The twin is
	// compared where both are valid.
	{name: "u16blank_plio", cmp: "PLIO_1", zbitpix: 16, cards: map[string]string{"BZERO": "32768", "BLANK": "-32743", "ZTILE2": "6"}, nullAt: []int{}, twin: twinBothValid},
	// 4-D and 6-D: cfitsio's fits_read_img refuses compressed images of
	// more than 3 dimensions (status 414), so their golden is cfitsio's
	// subset reader (fits_read_compressed_img; cfref.c readimg). That
	// reader is right for 4-D (and 5-D, probed) but misreads 6-D (pixels
	// shifted, the second half 0); there the funpack twin, equal to
	// astropy 8.0.2's decode and to the generated source, is the truth.
	{name: "cube4d16_rice", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"ZNAXIS": "4", "ZNAXIS4": "2", "ZTILE3": "2"}},
	{name: "cube6d16_gzip2", cmp: "GZIP_2", zbitpix: 16, cards: map[string]string{"ZNAXIS": "6", "ZNAXIS6": "2", "ZTILE4": "2"}, twinTruth: true},
	// BITPIX -64 with SUBTRACTIVE_DITHER_2 (exact zeros kept) and
	// NO_DITHER, and every algorithm x dither method for floats.
	{name: "f64_rice_sd2", cmp: "RICE_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_2", "ZDITHER0": "17"}},
	{name: "f64_rice_nodith", cmp: "RICE_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "NO_DITHER"}},
	{name: "f64_gzip2_sd2", cmp: "GZIP_2", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_2"}},
	{name: "f64_hcomp_nodith", cmp: "HCOMPRESS_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "NO_DITHER"}},
	{name: "f32_gzip1_sd2", cmp: "GZIP_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_2"}},
	{name: "f32_gzip1_nodith", cmp: "GZIP_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "NO_DITHER"}},
	{name: "f64_gzip1_sd2", cmp: "GZIP_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_2"}},
	{name: "f64_gzip1_nodith", cmp: "GZIP_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "NO_DITHER"}},
	{name: "f32_gzip2_sd2", cmp: "GZIP_2", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_2"}},
	{name: "f32_gzip2_nodith", cmp: "GZIP_2", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "NO_DITHER"}},
	{name: "f64_gzip2_nodith", cmp: "GZIP_2", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "NO_DITHER"}},
	{name: "f32_hcomp_nodith", cmp: "HCOMPRESS_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "NO_DITHER"}},
	// Real M13 (a 2006 SkyView 300x300 int16 image) as astropy 8.0.2 ships
	// it (BSD-3-Clause, verbatim). m13_gzip holds the 16-bit pixels in
	// 4-byte GZIP_1 tiles (1200 inflated bytes per 300-pixel row) and
	// m13_rice is RICE_1 without ZVAL2, so BYTEPIX 4: the layouts old
	// cfitsio wrote. cfitsio picks the tile type from the inflated length
	// (GZIP, imcompress.c:6372-6417) or from BYTEPIX (RICE, :6335-6346) and
	// converts to ZBITPIX. funpack and astropy 8.0.2 agree with it.
	{name: "real_m13_gzip", cmp: "GZIP_1", zbitpix: 16, cards: map[string]string{"ZNAXIS1": "300", "ZTILE2": "1"}, tileWidth: 4},
	{name: "real_m13_rice", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"ZNAXIS1": "300", "ZVAL1": "32"}, noCards: []string{"ZVAL2"}},
	{name: "real_m13_hcomp", cmp: "HCOMPRESS_1", zbitpix: 16, cards: map[string]string{"ZNAXIS1": "300", "ZTILE2": "300"}},
	{name: "real_m13_plio", cmp: "PLIO_1", zbitpix: 16, cards: map[string]string{"ZNAXIS1": "300"}},
	// The same tile-width rule on corpus images re-encoded by libcfitsio's
	// own coders (gen_fpack rewidth): every width cfitsio accepts.
	{name: "u8_gzip1_w4", cmp: "GZIP_1", zbitpix: 8, tileWidth: 4},
	{name: "u8_gzip2_w2", cmp: "GZIP_2", zbitpix: 8, tileWidth: 2},
	{name: "u16_gzip2_w4", cmp: "GZIP_2", zbitpix: 16, cards: map[string]string{"BZERO": "32768"}, tileWidth: 4},
	// libcfitsio's TDOUBLE read (the golden's phys) of 8-byte tiles on an
	// integer image returns memory it never wrote (imcomp_decompress_tile
	// has no branch for the tile type; unwrittenDouble); its TLONG read (the golden's i64),
	// funpack and astropy 8.0.2 give the source values, so the twin is the
	// truth.
	{name: "i16blank_gzip1_w8", cmp: "GZIP_1", zbitpix: 16, cards: map[string]string{"BLANK": "-32768"}, tileWidth: 8, twinTruth: true},
	{name: "i32_gzip2_w8", cmp: "GZIP_2", zbitpix: 32, tileWidth: 8, twinTruth: true},
	{name: "u8_rice_bp4", cmp: "RICE_1", zbitpix: 8, cards: map[string]string{"ZVAL2": "4"}},
	// BITPIX -64 twins of the scaled -32 cases: quantized + BSCALE/BZERO,
	// GZIP_COMPRESSED_DATA and UNCOMPRESSED_DATA fallback tiles + BSCALE/
	// BZERO (rawRows: see f32flat_rice_scaled), lossless + ZSCALE/ZZERO
	// keywords.
	{name: "f64scaled_rice", cmp: "RICE_1", zbitpix: -64, cards: map[string]string{"BSCALE": "2.", "BZERO": "10."}, twin: twinRescaled},
	{name: "f64flat_rice_scaled", cmp: "RICE_1", zbitpix: -64, cols: []string{"GZIP_COMPRESSED_DATA"}, cards: map[string]string{"BSCALE": "2.", "BZERO": "10."}, twin: twinRescaled, rawRows: []int{5}},
	{name: "f64flat_uncomp_scaled", cmp: "RICE_1", zbitpix: -64, cols: []string{"UNCOMPRESSED_DATA"}, cards: map[string]string{"BSCALE": "2.", "BZERO": "10."}, twin: twinRescaled, rawRows: []int{5}},
	{name: "f64_lossless_zscalekw", cmp: "GZIP_1", zbitpix: -64, cards: map[string]string{"ZSCALE": "2.", "ZZERO": "10."}},
	// ZCMPTYPE = 'RICE_ONE': cfitsio's alias for RICE_1 (imcompress.c:5279).
	{name: "u16_rice_riceone", cmp: "RICE_ONE", zbitpix: 16, cards: map[string]string{"BZERO": "32768"}},

	// Round 2. BITPIX -64 SUBTRACTIVE_DITHER_1 in one 102x100 tile (10200
	// pixels): the 10000-value dither sequence wraps inside the tile.
	{name: "f64_rice_wrap", cmp: "RICE_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1", "ZDITHER0": "17", "ZTILE1": "102", "ZTILE2": "100"}},
	// NaN inside float fallback tiles: row 5 is constant except two NaN,
	// so cfitsio stores it raw (NaN included) in GZIP_COMPRESSED_DATA,
	// moved to UNCOMPRESSED_DATA in the _uncomp cases; pixel 88 is a NaN
	// in a quantized row (ZBLANK). cfitsio flags all three null.
	{name: "f32flatnan_rice", cmp: "RICE_1", zbitpix: -32, cols: []string{"GZIP_COMPRESSED_DATA"}, nullAt: []int{88, 203, 217}},
	{name: "f64flatnan_rice", cmp: "RICE_1", zbitpix: -64, cols: []string{"GZIP_COMPRESSED_DATA"}, nullAt: []int{88, 203, 217}},
	{name: "f32flatnan_uncomp", cmp: "RICE_1", zbitpix: -32, cols: []string{"UNCOMPRESSED_DATA"}, cards: map[string]string{"TFORM4": "1PE(40)"}, nullAt: []int{88, 203, 217}},
	{name: "f64flatnan_uncomp", cmp: "RICE_1", zbitpix: -64, cols: []string{"UNCOMPRESSED_DATA"}, cards: map[string]string{"TFORM4": "1PD(40)"}, nullAt: []int{88, 203, 217}},
	// ZQUANTIZ = 'NONE' on f32_rice_sd1 (ZSCALE/ZZERO columns kept):
	// cfitsio NO_QUANTIZE takes the 4-byte tile integers as IEEE floats and
	// still scales them by the tile's ZSCALE/ZZERO (imcompress.c:5330,
	// fffr4r4); funpack agrees. Treating the file as quantized gives the
	// sd1 pixels, which is what astropy 8.0.2 returns for this file: the
	// library follows cfitsio.
	// Tile integer -1 of row 10 is a NaN bit pattern (see rawRows).
	{name: "f32_rice_zqnone", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "NONE"}, cols: []string{"ZSCALE", "ZZERO"}, rawRows: []int{10}},
	// The same with BSCALE 2 / BZERO 10, which cfitsio folds into the tile
	// ZSCALE/ZZERO (funpack writes the floats before that scaling).
	{name: "f32_rice_zqnone_scaled", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "NONE", "BSCALE": "2.", "BZERO": "10."}, cols: []string{"ZSCALE", "ZZERO"}, rawRows: []int{10}, twin: twinRescaled},
	// RICE_1 BYTEPIX and GZIP tile widths other than ZBITPIX's, re-encoded
	// by libcfitsio's own coders (gen_fpack rewidth / rewidthq) on sources
	// with stored values 0..255. cfitsio compares the null value cast to
	// the tile type: BLANK -32768 is 0 in 1-byte tiles, BLANK 65541 is 5
	// in 1- and 2-byte tiles, so those pixels are null (the funpack twin,
	// which does not null-check integer tiles, has them as values).
	{name: "i16small_rice_bp1", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"ZVAL2": "1", "BLANK": "-32768"}, twin: twinBothValid},
	{name: "i16small_gzip1_w1", cmp: "GZIP_1", zbitpix: 16, cards: map[string]string{"BLANK": "-32768"}, tileWidth: 1, twin: twinBothValid},
	{name: "i32small_rice_bp1", cmp: "RICE_1", zbitpix: 32, cards: map[string]string{"ZVAL2": "1", "BLANK": "65541"}, twin: twinBothValid},
	{name: "i32small_rice_bp2", cmp: "RICE_1", zbitpix: 32, cards: map[string]string{"ZVAL2": "2", "BLANK": "65541"}, twin: twinBothValid},
	{name: "i32small_gzip2_w1", cmp: "GZIP_2", zbitpix: 32, cards: map[string]string{"BLANK": "65541"}, tileWidth: 1, twin: twinBothValid},
	{name: "i32small_gzip2_w2", cmp: "GZIP_2", zbitpix: 32, cards: map[string]string{"BLANK": "65541"}, tileWidth: 2, twin: twinBothValid},
	// The same cast on native 2-byte tiles: BLANK 65536 is (short) 0.
	{name: "i16small_rice_blankcast", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"ZVAL2": "2", "BLANK": "65536"}, twin: twinBothValid, blankRefused: true},
	{name: "u8_rice_bp2", cmp: "RICE_1", zbitpix: 8, cards: map[string]string{"ZVAL2": "2"}, sameAs: "u8_rice"},
	// 8-byte GZIP_2 tiles on an 8-bit image: libcfitsio's TDOUBLE read
	// returns unwritten memory and funpack (fits_img_decompress, TBYTE
	// output) zeros;
	// cfitsio's TLONG read gives the source values, equal to astropy 8.0.2
	// and to the u8 twin, and is the truth.
	{name: "u8_gzip2_w8", cmp: "GZIP_2", zbitpix: 8, tileWidth: 8, twin: twinWrong, i64Truth: true, sameAs: "u8_gzip2"},
	// RICE_1 ZVAL2 = 8: cfitsio uses the 4-byte coder for every BYTEPIX
	// other than 1 and 2, so this reads as f64_rice_sd1.
	{name: "f64_rice_bp8", cmp: "RICE_1", zbitpix: -64, cards: map[string]string{"ZVAL2": "8", "ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}, sameAs: "f64_rice_sd1"},
	// Quantized float tiles in 2- and 1-byte Rice (unquantize_i2r4,
	// _i2r8; fffi1r4 for NO_DITHER).
	{name: "f32_rice_bp2", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZVAL2": "2", "ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}, sameAs: "f32_rice_sd1"},
	{name: "f64_rice_bp2", cmp: "RICE_1", zbitpix: -64, cards: map[string]string{"ZVAL2": "2", "ZQUANTIZ": "NO_DITHER"}, sameAs: "f64_rice_nodith"},
	{name: "f32_rice_bp1", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZVAL2": "1", "ZQUANTIZ": "NO_DITHER"}},

	// Round 3. The edges of the ZBITPIX range in wider tiles: tile 1
	// values 0 and 1 are the lowest and highest stored value the image
	// type holds (one past them: comp_bad/*_lo, *_hi, refused).
	{name: "u8_gzip2_w2_edge", cmp: "GZIP_2", zbitpix: 8, tileWidth: 2, at: map[int]float64{0: 0, 1: 255}},
	{name: "u16_gzip1_w4_edge", cmp: "GZIP_1", zbitpix: 16, cards: map[string]string{"BZERO": "32768"}, tileWidth: 4, at: map[int]float64{0: 0, 1: 65535}},
	{name: "i32_gzip1_w8_edge", cmp: "GZIP_1", zbitpix: 32, tileWidth: 8, twinTruth: true, at: map[int]float64{0: math.MinInt32, 1: math.MaxInt32}},
	// ZBLANK 406 on 1-byte RICE tiles of a quantized float: cfitsio
	// compares (unsigned char) 406 = 150 with the tile bytes (signed it
	// would be -106, and the uncast 406 matches no byte).
	{name: "f32_rice_bp1_zblank406", cmp: "RICE_1", zbitpix: -32, cards: map[string]string{"ZVAL2": "1", "ZBLANK": "406"}, nullAt: []int{622, 628, 668, 751, 753, 792}},
	// HCOMPRESS_1 and PLIO_1 tiles are 4-byte ints compared with the int
	// null value: BLANK 65536 flags none of the zeros of a 16-bit image.
	{name: "i16small_hcomp_blank65536", cmp: "HCOMPRESS_1", zbitpix: 16, cards: map[string]string{"BLANK": "65536"}, nullAt: []int{}, zeros: true, blankRefused: true},
	{name: "i16small_plio_blank65536", cmp: "PLIO_1", zbitpix: 16, cards: map[string]string{"BLANK": "65536"}, nullAt: []int{}, zeros: true, blankRefused: true},
	// Quantized float GZIP tiles of 1 and 2 bytes (cfitsio takes the tile
	// type from the inflated length: unquantize_i1r4, _i2r4, _i2r8).
	// astropy 8.0.2 agrees.
	{name: "f32_gzip1_q_w1", cmp: "GZIP_1", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "NO_DITHER"}, tileWidth: 1},
	{name: "f32_gzip2_q_w2", cmp: "GZIP_2", zbitpix: -32, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}, tileWidth: 2, sameAs: "f32_gzip2_sd1"},
	{name: "f64_gzip1_q_w2", cmp: "GZIP_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "SUBTRACTIVE_DITHER_1"}, tileWidth: 2, sameAs: "f64_gzip1_sd1"},
	// ZQUANTIZ 'NONE' on a -64 RICE_1 file with ZSCALE/ZZERO columns:
	// cfitsio reads the 4-byte tile integers as IEEE float32 (fffr4r8) and
	// scales them; two of them are NaN patterns (rows 5 and 15, see
	// rawRows). astropy 8.0.2 dequantizes instead.
	{name: "f64_rice_zqnone", cmp: "RICE_1", zbitpix: -64, cards: map[string]string{"ZQUANTIZ": "NONE", "ZVAL2": "4"}, cols: []string{"ZSCALE", "ZZERO"}, rawRows: []int{5, 15}},

	// Round 2 fix pass. Lossy HCOMPRESS of 16-bit images saturated at the
	// type maximum: pixels reconstructed above it are clipped by cfitsio's
	// TSHORT read and by funpack (the twin, the truth; cfitsio's TUSHORT
	// read of the uint16 file fails with 412);
	// cfitsio's TDOUBLE read returns them unclipped (clipTruth pins that).
	{name: "u16sat_hcomp_s2", cmp: "HCOMPRESS_1", zbitpix: 16, cards: map[string]string{"BZERO": "32768", "ZVAL1": "2.0000000E+00"}, clipTruth: true},
	{name: "i16sat_hcomp_s2", cmp: "HCOMPRESS_1", zbitpix: 16, cards: map[string]string{"ZVAL1": "2.0000000E+00"}, clipTruth: true},
	// The null value at both ends of the int32 range (4-byte GZIP_1 tiles,
	// pixel 0 = -2^31, pixel 1 = 2^31 - 1): BLANK -2^31 and 2^31 - 1, and a
	// 1K ZBLANK column of -2^31 (one past either end: comp_bad, refused).
	{name: "i32_gzip1_blankmin", cmp: "GZIP_1", zbitpix: 32, cards: map[string]string{"BLANK": "-2147483648"}, tileWidth: 4, nullAt: []int{0}, at: map[int]float64{1: math.MaxInt32}},
	{name: "i32_gzip1_blankmax", cmp: "GZIP_1", zbitpix: 32, cards: map[string]string{"BLANK": "2147483647"}, tileWidth: 4, nullAt: []int{1}, at: map[int]float64{0: math.MinInt32}},
	{name: "i32_gzip1_zblankcolmin", cmp: "GZIP_1", zbitpix: 32, cols: []string{"ZBLANK"}, noCards: []string{"BLANK"}, tileWidth: 4, nullAt: []int{0}, at: map[int]float64{1: math.MaxInt32}, twin: twinNoNull},
	// A ZBLANK column next to a BLANK keyword outside int32: cfitsio uses
	// the column and never reads the keyword as the null value.
	{name: "i16_zblankcol_blank2p32", cmp: "RICE_1", zbitpix: 16, cards: map[string]string{"BLANK": "4294967296"}, cols: []string{"ZBLANK"}, twin: twinNoNull, blankRefused: true},
	// A BLANK the decoded HCOMPRESS_1 / PLIO_1 ints hold: flagged.
	{name: "i16small_hcomp_blank0", cmp: "HCOMPRESS_1", zbitpix: 16, cards: map[string]string{"BLANK": "0"}, minNulls: 10},
	{name: "i16mask_plio_blank3", cmp: "PLIO_1", zbitpix: 16, cards: map[string]string{"BLANK": "3"}, minNulls: 100},
}

// corpusMatrix is every cell the README and the ledger claim for the
// compressed reader: each algorithm x each ZBITPIX x each BZERO
// convention (integers), and each algorithm x float ZBITPIX x each
// ZQUANTIZ (floats). A cell is "<ZCMPTYPE> <ZBITPIX> <BZERO or ZQUANTIZ>".
// Cells cfitsio itself cannot write are not claimed and are listed in
// corpusExcluded with the reason. TestCompressedCorpusComplete reads the
// cells of every corpus file from its cfitsio golden and fails on a
// missing cell.
func corpusMatrix() []string {
	var cells []string
	ints := []string{"8 0", "8 -128", "16 0", "16 32768", "32 0", "32 2147483648"}
	for _, a := range []string{"RICE_1", "GZIP_1", "GZIP_2", "HCOMPRESS_1", "PLIO_1"} {
		for _, c := range ints {
			cells = append(cells, a+" "+c)
		}
	}
	for _, a := range []string{"RICE_1", "GZIP_1", "GZIP_2", "HCOMPRESS_1"} {
		for _, bp := range []string{"-32", "-64"} {
			for _, q := range []string{"SUBTRACTIVE_DITHER_1", "SUBTRACTIVE_DITHER_2", "NO_DITHER", "NONE"} {
				cells = append(cells, a+" "+bp+" "+q)
			}
		}
	}
	cells = append(cells, "RICE_ONE 16 32768")
	return cells
}

var corpusExcluded = map[string]string{
	"PLIO_1 32 2147483648": "uint32 values exceed PLIO's 0..2^24-1 (cfitsio status 413, cref plio_u32)",
	// RICE_1 -32/-64 NONE: fpack will not write them; f32_rice_zqnone and
	// f64_rice_zqnone are the cfitsio NO_QUANTIZE readings of RICE_1 files
	// with ZQUANTIZ 'NONE'.
	"HCOMPRESS_1 -32 NONE":                 "fpack: '-q 0' only allowed with GZIP",
	"HCOMPRESS_1 -64 NONE":                 "fpack: '-q 0' only allowed with GZIP",
	"HCOMPRESS_1 -32 SUBTRACTIVE_DITHER_2": "cfitsio does not allow SUBTRACTIVE_DITHER_2 with HCOMPRESS (fpack falls back to _1)",
	"HCOMPRESS_1 -64 SUBTRACTIVE_DITHER_2": "cfitsio does not allow SUBTRACTIVE_DITHER_2 with HCOMPRESS (fpack falls back to _1)",
}

// TestCompressedCorpusComplete pins the corpus: every comp_fz_ file has a
// case (and twin), every case has its files. Adding or losing a sample
// without updating the test fails here.
func TestCompressedCorpusComplete(t *testing.T) {
	const wantCases = 148
	have := map[string]bool{}
	for _, c := range fzCorpus {
		g := testutil.LoadCfrefGolden(t, "comp_fz_"+c.name+".fits")
		for i := range g.HDUs {
			h := &g.HDUs[i]
			if h.Image == nil || !h.Image.Compressed {
				continue
			}
			cmp, _ := goldenCardValue(h, "ZCMPTYPE")
			var cell string
			if h.Image.Bitpix > 0 {
				cell = fmt.Sprintf("%s %d %.0f", cmp, h.Image.Bitpix, float64(h.Image.Bzero))
			} else {
				// cfitsio: a float image is quantized unless it has no
				// ZSCALE column (then lossless, whatever ZQUANTIZ says) or
				// ZQUANTIZ is 'NONE'; an absent ZQUANTIZ is NO_DITHER
				// (imcompress.c:5322-5330; compressed.go parses it so).
				q, ok := goldenCardValue(h, "ZQUANTIZ")
				quant := false
				for _, cd := range h.Cards {
					if strings.HasPrefix(string(cd.Key), "TTYPE") && cardStr(cd) == "ZSCALE" {
						quant = true
					}
				}
				switch {
				case !quant:
					q = "NONE"
				case !ok:
					q = "NO_DITHER"
				}
				cell = fmt.Sprintf("%s %d %s", cmp, h.Image.Bitpix, q)
			}
			have[cell] = true
		}
	}
	for _, cell := range corpusMatrix() {
		why, excluded := corpusExcluded[cell]
		switch {
		case excluded && have[cell]:
			t.Errorf("cell %q is excluded (%s) but the corpus has a sample: claim it", cell, why)
		case !excluded && !have[cell]:
			t.Errorf("claimed cell %q has no sample in the compressed read corpus", cell)
		}
	}
	for cell := range corpusExcluded {
		if !slices.Contains(corpusMatrix(), cell) {
			t.Errorf("excluded cell %q is not a matrix cell", cell)
		}
	}
	// The scaled float layouts (image BSCALE/BZERO on quantized and on
	// fallback tiles, ZSCALE/ZZERO keywords on a lossless image) for BOTH
	// float ZBITPIX: they are not matrix cells but the reader claims them.
	for _, bp := range []string{"32", "64"} {
		for _, n := range []string{"f%sscaled_rice", "f%sflat_rice_scaled", "f%sflat_uncomp_scaled", "f%s_lossless_zscalekw"} {
			name := fmt.Sprintf(n, bp)
			if !slices.ContainsFunc(fzCorpus, func(c fzCase) bool { return c.name == name }) {
				t.Errorf("scaled float layout %s has no corpus case", name)
			}
		}
	}
	if len(fzCorpus) != wantCases {
		t.Fatalf("fzCorpus has %d cases, want %d", len(fzCorpus), wantCases)
	}
	glob := func(p string) map[string]bool {
		m, err := filepath.Glob(testutil.TestdataPath(t, p))
		if err != nil {
			t.Fatal(err)
		}
		out := map[string]bool{}
		for _, f := range m {
			out[filepath.Base(f)] = true
		}
		return out
	}
	fz := glob("comp_fz_*.fits")
	un := map[string]bool{}
	for g := range glob(filepath.Join("cfref", "golden", "twins", "comp_fzun_*.fits.json")) {
		un[strings.TrimSuffix(g, ".json")] = true
	}
	if committed := glob("comp_fzun_*.fits"); len(committed) != 0 {
		t.Errorf("%d funpack twins are committed in testdata/; they belong in testdata/cfref/twins.out (gen_twins.sh)", len(committed))
	}
	seen := map[string]bool{}
	usedTwin := map[string]bool{}
	names := make([]string, 0, len(fzCorpus)+len(fz64Corpus))
	for _, c := range fzCorpus {
		names = append(names, c.name)
	}
	for _, c := range fz64Corpus {
		names = append(names, c.name) // ZBITPIX 64: TestCompressedInt64Read
	}
	for _, name := range names {
		c := struct{ name string }{name}
		if seen[c.name] {
			t.Errorf("duplicate case %s", c.name)
		}
		seen[c.name] = true
		if !fz["comp_fz_"+c.name+".fits"] {
			t.Errorf("case %s: comp_fz_%s.fits missing", c.name, c.name)
		}
		tw := twinFile(t, c.name)
		if !un[tw] {
			t.Errorf("case %s: twin golden cfref/golden/twins/%s.json missing", c.name, tw)
		}
		if own := "comp_fzun_" + c.name + ".fits"; tw != own && un[own] {
			t.Errorf("case %s: %s is committed and also mapped to %s in comp_fzun_same.txt", c.name, own, tw)
		}
		usedTwin[tw] = true
	}
	for f := range fz {
		if !seen[strings.TrimSuffix(strings.TrimPrefix(f, "comp_fz_"), ".fits")] {
			t.Errorf("%s has no test case", f)
		}
	}
	for f := range un {
		if !usedTwin[f] {
			t.Errorf("%s is no case's twin", f)
		}
	}
	if total := wantCases + len(fz64Corpus); len(fz) != total || len(un) != len(usedTwin) || len(twinSame)+len(un) != total {
		t.Errorf("found %d comp_fz_ files (want %d), %d comp_fzun_ twins + %d mapped (want %d)", len(fz), total, len(un), len(twinSame), total)
	}
}

// inflatedTileWidth inflates tile 1 of a GZIP_1/GZIP_2 corpus file with
// the standard library (gzip, or zlib as astropy writes) and returns its
// length divided by the tile's pixel count: the width of the tile values
// the writer stored, independent of the library's decoder.
func inflatedTileWidth(t *testing.T, name string) int {
	t.Helper()
	f, err := Open(testutil.TestdataPath(t, name))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	hdu, err := f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	h := hdu.(*CompressedImageHDU)
	meta, err := h.metadata()
	if err != nil {
		t.Fatal(err)
	}
	cols, err := h.tbl.Columns()
	if err != nil {
		t.Fatal(err)
	}
	p, err := readVarBytes(h.tbl, findCompColumn(cols, "COMPRESSED_DATA"), h.tbl.NumRows())
	if err != nil {
		t.Fatal(err)
	}
	var r io.Reader
	if gz, err := gzip.NewReader(bytes.NewReader(p[0])); err == nil {
		r = gz
	} else if zr, err := zlib.NewReader(bytes.NewReader(p[0])); err == nil {
		r = zr
	} else {
		t.Fatalf("%s tile 1 is neither gzip nor zlib", name)
	}
	b, err := io.ReadAll(r)
	if err != nil {
		t.Fatal(err)
	}
	box := tileBox(meta.Shape, meta.Tile, make([]int64, meta.Naxis))
	n := int64(1)
	for _, d := range box.size {
		n *= d
	}
	if int64(len(b))%n != 0 {
		t.Fatalf("%s tile 1 inflates to %d bytes, not a multiple of its %d pixels", name, len(b), n)
	}
	return int(int64(len(b)) / n)
}

// cfrefPixels returns the image HDU's pixel truth from a cfitsio golden:
// values (float32-rounded when f32) and the null flags.
func cfrefPixels(t *testing.T, name string, f32 bool) (*testutil.CfrefHDU, []float64, []bool) {
	t.Helper()
	g := testutil.LoadCfrefGolden(t, name)
	for i := range g.HDUs {
		h := &g.HDUs[i]
		if h.Image == nil || h.Image.Npix == 0 {
			continue
		}
		if len(h.Errors) != 0 {
			t.Fatalf("%s: cfitsio errors in golden: %v", name, h.Errors)
		}
		vals := make([]float64, h.Image.Phys.Len())
		null := make([]bool, h.Image.Phys.Len())
		var phys []testutil.NullFloat
		if f32 {
			phys = h.Image.PhysAsFloat32(t, name+" cfitsio TDOUBLE as float32", func() []testutil.NullFloat {
				return roundF32(ourPhysView(t, name, i))
			})
		} else {
			phys = h.Image.Phys.Resolve(t, name+" cfitsio TDOUBLE", func() []testutil.NullFloat { return ourPhysView(t, name, i) })
		}
		for k, p := range phys {
			if !p.Valid {
				null[k] = true
				vals[k] = math.NaN()
				continue
			}
			vals[k] = p.V
		}
		return h, vals, null
	}
	t.Fatalf("%s: golden has no image pixels", name)
	return nil, nil, nil
}

// ourPhysView is the library's masked float64 reading of HDU i of the
// testdata sample name, in cfitsio's TDOUBLE terms (see cfitsioPhysView):
// the candidate that testutil.Values.Resolve verifies against a digested
// golden before a test uses the golden's pixels one by one.
func ourPhysView(t *testing.T, name string, i int) []testutil.NullFloat {
	t.Helper()
	if strings.HasPrefix(name, "comp_fzun_") {
		// A twin is not committed (only its golden): the candidate is the
		// library's read of the compressed sample it was made from, which
		// a twinExact case must equal anyway (comp_fzun_<case> is the
		// twin of comp_fz_<case>; a mapped twin is kept under the name of
		// the first case that has it).
		name, i = "comp_fz_"+strings.TrimPrefix(name, "comp_fzun_"), 1
	}
	f, err := Open(testutil.TestdataPath(t, name))
	if err != nil {
		t.Fatalf("%s: %v", name, err)
	}
	defer f.Close()
	hdu, err := f.HDU(i)
	if err != nil {
		t.Fatalf("%s HDU %d: %v", name, i, err)
	}
	switch h := hdu.(type) {
	case *ImageHDU:
		return cfitsioPhysView(t, h)
	case *CompressedImageHDU:
		v, ok, err := ReadPixelsCompressedMasked[float64](h)
		if err != nil {
			return nil
		}
		return testutil.NullFloats(v, ok)
	}
	return nil
}

// roundF32 rounds every valid value to float32 (nil stays nil).
func roundF32(v []testutil.NullFloat) []testutil.NullFloat {
	if v == nil {
		return nil
	}
	out := append([]testutil.NullFloat(nil), v...)
	for k := range out {
		if out[k].Valid {
			out[k].V = float64(float32(out[k].V))
		}
	}
	return out
}

func goldenCardValue(h *testutil.CfrefHDU, key string) (string, bool) {
	c := h.Card(key)
	if c == nil {
		return "", false
	}
	if c.Str != nil {
		return strings.TrimSpace(string(*c.Str)), true
	}
	return strings.TrimSpace(string(c.Value)), true
}

// sameFloat (images_truth_test.go): bit-exact, every NaN equal to every NaN.

type pixDiff struct {
	n     int
	first int
	got   float64
	want  float64
	maxAb float64
}

func (d pixDiff) String() string {
	return fmt.Sprintf("%d pixels differ; first at %d: got %v want %v; max |diff| %g", d.n, d.first, d.got, d.want, d.maxAb)
}

func diffPix(got, want []float64, skip []bool) pixDiff {
	d := pixDiff{first: -1}
	for i := range want {
		if skip != nil && skip[i] {
			continue
		}
		if !sameFloat(got[i], want[i]) {
			if d.first < 0 {
				d.first, d.got, d.want = i, got[i], want[i]
			}
			d.n++
			if a := math.Abs(got[i] - want[i]); a > d.maxAb || math.IsNaN(a) {
				d.maxAb = a
			}
		}
	}
	return d
}

func toF64[T Numeric](v []T) []float64 {
	out := make([]float64, len(v))
	for i, x := range v {
		out[i] = float64(x)
	}
	return out
}

// readNative reads h in the Go type that holds its logical pixels exactly
// (uint16 for BZERO=32768, int8 for BZERO=-128, ...), widened to float64.
func readNative(h *CompressedImageHDU, zbitpix int, bscale, bzero float64) (string, []float64, error) {
	switch {
	case zbitpix == 8 && bscale == 1 && bzero == -128:
		v, err := ReadPixelsCompressed[int8](h)
		return "int8", toF64(v), err
	case zbitpix == 8 && bscale == 1 && bzero == 0:
		v, err := ReadPixelsCompressed[uint8](h)
		return "uint8", toF64(v), err
	case zbitpix == 16 && bscale == 1 && bzero == 32768:
		v, err := ReadPixelsCompressed[uint16](h)
		return "uint16", toF64(v), err
	case zbitpix == 16 && bscale == 1 && bzero == 0:
		v, err := ReadPixelsCompressed[int16](h)
		return "int16", toF64(v), err
	case zbitpix == 32 && bscale == 1 && bzero == 2147483648:
		v, err := ReadPixelsCompressed[uint32](h)
		return "uint32", toF64(v), err
	case zbitpix == 32 && bscale == 1 && bzero == 0:
		v, err := ReadPixelsCompressed[int32](h)
		return "int32", toF64(v), err
	case zbitpix == -32:
		v, err := ReadPixelsCompressed[float32](h)
		return "float32", toF64(v), err
	}
	v, err := ReadPixelsCompressed[float64](h)
	return "float64", v, err
}

// checkCorpusScalingRefused: for a blankRefused case, every pixel read of
// the compressed HDU and BSCALE / BZERO fail with ErrInvalidScaling while
// BITPIX and Shape stay readable. It reports whether c is such a case.
func checkCorpusScalingRefused(t *testing.T, c fzCase) bool {
	t.Helper()
	if !c.blankRefused {
		return false
	}
	f, h := openCompressed(t, filepathTestdata(t, "comp_fz_"+c.name+".fits"))
	defer f.Close()
	if h.BITPIX() != c.zbitpix || len(h.Shape()) == 0 {
		t.Errorf("BITPIX %d, Shape %v: the structure must stay readable", h.BITPIX(), h.Shape())
	}
	is := func(what string, err error) {
		t.Helper()
		if !errors.Is(err, ErrInvalidScaling) || !strings.Contains(err.Error(), "BLANK") {
			t.Errorf("%s: %v, want ErrInvalidScaling naming BLANK", what, err)
		}
	}
	_, err := h.BSCALE()
	is("BSCALE", err)
	_, err = h.BZERO()
	is("BZERO", err)
	if err := h.Validate(); err != nil {
		t.Errorf("Validate: %v; the refused cards are not a structural error", err)
	}
	_, err = ReadPixelsCompressed[float64](h)
	is("ReadPixelsCompressed[float64]", err)
	_, err = ReadPixelsCompressed[int16](h)
	is("ReadPixelsCompressed[int16]", err)
	_, _, err = ReadPixelsCompressedMasked[float64](h)
	is("ReadPixelsCompressedMasked[float64]", err)
	_, err = ReadFloat32(h)
	is("ReadFloat32", err)
	return true
}

func TestCompressedCorpusRead(t *testing.T) {
	for _, c := range fzCorpus {
		t.Run(c.name, func(t *testing.T) {
			if checkCorpusScalingRefused(t, c) {
				return
			}
			fzName := "comp_fz_" + c.name + ".fits"
			unName := twinFile(t, c.name)
			f32 := c.zbitpix == -32
			gh, want, null := cfrefPixels(t, fzName, f32)
			im := gh.Image

			// The case really exercises what it claims (guards the corpus
			// against a generator that silently changes a variant).
			if v, _ := goldenCardValue(gh, "ZCMPTYPE"); v != c.cmp {
				t.Fatalf("golden ZCMPTYPE %q, case says %q", v, c.cmp)
			}
			if im.Bitpix != c.zbitpix || !im.Compressed {
				t.Fatalf("golden bitpix %d compressed %v, case says ZBITPIX %d", im.Bitpix, im.Compressed, c.zbitpix)
			}
			keys := make([]string, 0, len(c.cards))
			for k := range c.cards {
				keys = append(keys, k)
			}
			sort.Strings(keys)
			for _, k := range keys {
				v, ok := goldenCardValue(gh, k)
				if !ok || (c.cards[k] != "" && v != c.cards[k]) {
					t.Fatalf("golden card %s = %q (present %v), case needs %q", k, v, ok, c.cards[k])
				}
			}
			for _, k := range c.noCards {
				if _, ok := goldenCardValue(gh, k); ok {
					t.Fatalf("golden has card %s, case needs it absent", k)
				}
			}
			if c.tileWidth != 0 {
				w := inflatedTileWidth(t, fzName)
				if w != c.tileWidth {
					t.Fatalf("tile 1 inflates to %d bytes per pixel, case says %d", w, c.tileWidth)
				}
			}
			for _, col := range c.cols {
				found := false
				for _, cd := range gh.Cards {
					if strings.HasPrefix(string(cd.Key), "TTYPE") && cardStr(cd) == col {
						found = true
					}
				}
				if !found {
					t.Fatalf("golden has no %s column", col)
				}
			}
			_, twin, twinNull := cfrefPixels(t, unName, f32)
			if c.tileWidth == 8 {
				// libcfitsio's TDOUBLE read of 8-byte integer tiles copies
				// out its tile buffer unwritten (uninitialized memory; the
				// cfref harness fills malloc'd memory with 0xA5 bytes, so
				// every pixel is unwrittenDouble); its TLONG read is the
				// value the truth below holds.
				for k, p := range gh.Image.Phys.Values(t) {
					if !p.Valid || math.Float64bits(p.V) != unwrittenDouble {
						t.Fatalf("cfitsio TDOUBLE pixel %d = %v (valid %v); the 8-byte-tile note says cfitsio leaves it unwritten (%#x): update it", k, p.V, p.Valid, uint64(unwrittenDouble))
					}
				}
				iv, in := cfrefI64(t, fzName)
				truth, truthNull := twin, twinNull
				if !c.twinTruth {
					truth, truthNull = iv, in
				}
				if d := diffPix(iv, truth, nil); d.n > 0 || !slices.Equal(in, truthNull) {
					t.Fatalf("cfitsio TLONG read of 8-byte tiles differs from the truth: %s", d)
				}
			}
			if c.i64Truth {
				want, null = cfrefI64(t, fzName)
				if f32 {
					t.Fatalf("i64Truth on a float image")
				}
			}
			if c.sameAs != "" {
				_, sv, sn := cfrefPixels(t, twinFile(t, c.sameAs), f32)
				if d := diffPix(want, sv, nil); d.n > 0 || !slices.Equal(null, sn) {
					t.Fatalf("truth differs from the twin of %s (same source): %s", c.sameAs, d)
				}
			}
			if c.twinTruth {
				if d := diffPix(want, twin, nil); d.n == 0 {
					t.Fatalf("libcfitsio's decode now equals the funpack twin; drop twinTruth")
				} else {
					t.Logf("libcfitsio subset decode differs from the funpack twin (%s); twin used as truth", d)
				}
				want, null = twin, twinNull
			}
			if c.clipTruth {
				lo, hi := -32768.0, 32767.0
				if c.zbitpix != 16 {
					t.Fatalf("clipTruth on ZBITPIX %d", c.zbitpix)
				}
				lo, hi = lo*float64(im.Bscale)+float64(im.Bzero), hi*float64(im.Bscale)+float64(im.Bzero)
				clipped := 0
				for i := range want {
					v := want[i]
					if v < lo || v > hi {
						clipped++
						if twin[i] != min(max(v, lo), hi) || twinNull[i] {
							t.Fatalf("pixel %d: cfitsio TDOUBLE %v is outside [%v, %v]; the twin holds %v, not the clipped value", i, v, lo, hi, twin[i])
						}
					} else if !sameFloat(twin[i], v) || twinNull[i] != null[i] {
						t.Fatalf("pixel %d: cfitsio TDOUBLE %v (in range) differs from the twin %v", i, v, twin[i])
					}
				}
				if clipped == 0 {
					t.Fatalf("no reconstructed pixel is outside the ZBITPIX range: the case proves nothing")
				}
				t.Logf("cfitsio's TDOUBLE read returns %d pixels outside [%v, %v]; funpack and cfitsio's TSHORT read clip them (twin used as truth)", clipped, lo, hi)
				want, null = twin, twinNull
			}
			if c.rawRows != nil {
				nx := int(im.Naxes[0])
				differ := 0
				for _, r := range c.rawRows {
					for i := r * nx; i < (r+1)*nx; i++ {
						if null[i] != twinNull[i] || (!null[i] && !sameFloat(want[i], twin[i])) {
							differ++
						}
						want[i], null[i] = twin[i], twinNull[i]
					}
				}
				if differ == 0 {
					t.Errorf("libcfitsio now agrees with funpack on fallback rows %v; drop rawRows and the note", c.rawRows)
				} else {
					t.Logf("libcfitsio fits_read_img leaves %d pixels of fallback rows %v unscaled / unnulled; funpack twin used there", differ, c.rawRows)
				}
			}
			for k, v := range c.at {
				if null[k] || want[k] != v {
					t.Fatalf("truth pixel %d = %v (null %v), case needs %v", k, want[k], null[k], v)
				}
			}
			if c.zeros && !slices.ContainsFunc(want, func(v float64) bool { return v == 0 }) {
				t.Fatalf("truth has no valid zero pixel")
			}
			nulls := 0
			for _, n := range null {
				if n {
					nulls++
				}
			}
			if nulls < c.minNulls {
				t.Fatalf("truth has %d null pixels, the case needs at least %d", nulls, c.minNulls)
			}

			f, err := Open(testutil.TestdataPath(t, fzName))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			hdu, err := f.HDU(1)
			if err != nil {
				t.Fatal(err)
			}
			h, ok := hdu.(*CompressedImageHDU)
			if !ok {
				t.Fatalf("HDU 1 is %T", hdu)
			}
			wantCmp := c.cmp
			if wantCmp == "RICE_ONE" {
				wantCmp = "RICE_1" // CompressionType names the algorithm
			}
			if h.CompressionType() != wantCmp || h.BITPIX() != c.zbitpix {
				t.Errorf("CompressionType/BITPIX = %s/%d, want %s/%d", h.CompressionType(), h.BITPIX(), wantCmp, c.zbitpix)
			}
			if fmt.Sprint(h.Shape()) != fmt.Sprint(im.Naxes) {
				t.Errorf("Shape %v, cfitsio %v", h.Shape(), im.Naxes)
			}
			if bscaleOf(t, h) != float64(im.Bscale) || bzeroOf(t, h) != float64(im.Bzero) {
				t.Errorf("BSCALE/BZERO = %v/%v, cfitsio %v/%v", bscaleOf(t, h), bzeroOf(t, h), float64(im.Bscale), float64(im.Bzero))
			}

			// Integer null pixels have no defined value in the unmasked
			// read (they are the stored BLANK, scaled); float nulls are NaN.
			var skip []bool
			if c.zbitpix > 0 {
				skip = null
			}
			got64, err := ReadPixelsCompressed[float64](h)
			if err != nil {
				t.Fatalf("ReadPixelsCompressed[float64]: %v", err)
			}
			want64 := want
			if c.scaledF64 {
				_, want64, _ = cfrefPixels(t, fzName, false)
				if d := diffPix(want64, want, nil); d.n == 0 {
					t.Fatalf("scaledF64 case has no pixel where double and float32 scaling differ")
				}
			}
			if d := diffPix(got64, want64, skip); d.n > 0 {
				t.Errorf("ReadPixelsCompressed[float64] vs cfitsio: %s", d)
			}
			typ, gotN, err := readNative(h, c.zbitpix, float64(im.Bscale), float64(im.Bzero))
			if err != nil {
				t.Fatalf("ReadPixelsCompressed[%s]: %v", typ, err)
			}
			if d := diffPix(gotN, want, skip); d.n > 0 {
				t.Errorf("ReadPixelsCompressed[%s] vs cfitsio: %s", typ, d)
			}

			vals, mask, err := ReadPixelsCompressedMasked[float64](h)
			if err != nil {
				t.Fatalf("ReadPixelsCompressedMasked: %v", err)
			}
			badMask, badVal := 0, 0
			for i := range want {
				if mask[i] == null[i] {
					badMask++
					continue
				}
				if null[i] {
					if vals[i] != 0 {
						badVal++
					}
				} else if !sameFloat(vals[i], want64[i]) {
					badVal++
				}
			}
			if badMask > 0 || badVal > 0 {
				t.Errorf("ReadPixelsCompressedMasked: %d mask flags and %d values differ from cfitsio (%d nulls)", badMask, badVal, nulls)
			}

			// The funpack twin (an independent cfitsio code path writing a
			// plain image) must agree with the library decode.
			switch c.twin {
			case twinExact:
				if d := diffPix(twin, want, nil); d.n > 0 {
					t.Errorf("funpack twin disagrees with cfitsio decode: %s", d)
				}
			case twinNoNull:
				if d := diffPix(twin, want, null); d.n > 0 {
					t.Errorf("funpack twin (off-null) disagrees with cfitsio decode: %s", d)
				}
			case twinRescaled:
				for i := range want {
					a, b := float32(want[i]), float32(twin[i])
					if math.IsNaN(want[i]) && math.IsNaN(twin[i]) {
						continue
					}
					if math.IsNaN(want[i]) != math.IsNaN(twin[i]) || (a != b && math.Nextafter32(a, b) != b) {
						t.Fatalf("pixel %d: twin %v not within 1 float32 ULP of cfitsio decode %v", i, twin[i], want[i])
					}
				}
			case twinBothValid:
				skip := make([]bool, len(null))
				for i := range skip {
					skip[i] = null[i] || twinNull[i]
				}
				if d := diffPix(twin, want, skip); d.n > 0 {
					t.Errorf("funpack twin (where both valid) disagrees with cfitsio decode: %s", d)
				}
			case twinWrong:
				if d := diffPix(twin, want, nil); d.n == 0 {
					t.Errorf("funpack twin now agrees with the library decode; drop twinWrong for %s", c.name)
				} else {
					t.Logf("known funpack 1.7.0 (CFITSIO 4.7.0) error on this file: %s", d)
				}
			}

			if c.nullAt != nil {
				var at []int
				for i, n := range null {
					if n {
						at = append(at, i)
					}
				}
				if fmt.Sprint(at) != fmt.Sprint(c.nullAt) {
					t.Errorf("cfitsio null pixels %v, source has %v", at, c.nullAt)
				}
			}
			if c.src != "" {
				_, srcV, srcNull := cfrefPixels(t, "comp_fzsrc_"+c.src+".fits", f32)
				for i := range want {
					if srcNull[i] != null[i] {
						t.Fatalf("pixel %d: cfitsio decode null=%v, source null=%v", i, null[i], srcNull[i])
					}
				}
				if c.exact {
					if d := diffPix(want, srcV, nil); d.n > 0 {
						t.Errorf("lossless case differs from source: %s", d)
					}
				}
				nx := int(im.Naxes[0])
				for _, r := range c.exactRows {
					for i := r * nx; i < (r+1)*nx; i++ {
						if !sameFloat(want[i], srcV[i]) || !sameFloat(got64[i], srcV[i]) {
							t.Fatalf("row %d pixel %d: cfitsio %v, ours %v, source %v (row must be lossless)", r, i, want[i], got64[i], srcV[i])
						}
					}
				}
			}
		})
	}
}

func cardStr(c testutil.CfrefCard) string {
	if c.Str != nil {
		return strings.TrimSpace(string(*c.Str))
	}
	return strings.Trim(strings.TrimSpace(string(c.Value)), "'")
}

// TestReadFloat32NullsMatchCfitsioCorpus: for every integer case of the
// compressed read corpus, the NaN pixels of ReadFloat32 are exactly the
// pixels cfitsio flags null (cfref golden; twin for twinTruth cases). This
// pins the null rule on the ReadFloat32 path too, e.g. PLIO_1 uint16 with
// BLANK (u16blank_plio: cfitsio flags none, so no NaN).
func TestReadFloat32NullsMatchCfitsioCorpus(t *testing.T) {
	n := 0
	for _, c := range fzCorpus {
		if c.zbitpix < 0 {
			continue
		}
		t.Run(c.name, func(t *testing.T) {
			if checkCorpusScalingRefused(t, c) {
				return
			}
			_, null := caseTruth(t, c, false)
			if c.rawRows != nil {
				// fits_read_img does not null-check fallback rows; the
				// funpack twin is the truth there (see fzCase.rawRows).
				_, _, twinNull := cfrefPixels(t, twinFile(t, c.name), false)
				gh, _, _ := cfrefPixels(t, "comp_fz_"+c.name+".fits", false)
				nx := int(gh.Image.Naxes[0])
				for _, r := range c.rawRows {
					copy(null[r*nx:(r+1)*nx], twinNull[r*nx:(r+1)*nx])
				}
			}
			f, h := openCompressed(t, filepathTestdata(t, "comp_fz_"+c.name+".fits"))
			defer f.Close()
			got, err := ReadFloat32(h)
			if err != nil {
				t.Fatalf("ReadFloat32: %v", err)
			}
			if len(got) != len(null) {
				t.Fatalf("%d pixels, cfitsio %d", len(got), len(null))
			}
			for i := range null {
				if math.IsNaN(float64(got[i])) != null[i] {
					t.Fatalf("pixel %d: ReadFloat32 NaN=%v, cfitsio null=%v", i, math.IsNaN(float64(got[i])), null[i])
				}
			}
		})
		n++
	}
	if n < 30 {
		t.Fatalf("only %d integer corpus cases", n)
	}
}
