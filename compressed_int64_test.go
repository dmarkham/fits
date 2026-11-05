package fits

import (
	"errors"
	"math"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// fz64Corpus: the ZBITPIX 64 cases of the compressed read corpus
// (gen_fpack.sh, written by astropy 8.0.2 since fpack and
// fits_img_compress refuse 64-bit images; the twin is astropy's decode as
// a plain image, as funpack cannot decompress them). The 1- and 4-byte
// tile cases are re-encoded by libcfitsio's own GZIP coder (gen_fpack
// rewidth).
var fz64Corpus = []struct {
	name, cmp string
	unsigned  bool // uint64: BZERO = 2^63
	tileWidth int
}{
	{"i64_gzip1", "GZIP_1", false, 8},
	{"i64_gzip2", "GZIP_2", false, 8},
	{"u64_gzip1", "GZIP_1", true, 8},
	{"u64_gzip2", "GZIP_2", true, 8},
	{"i64small_gzip1_w4", "GZIP_1", false, 4},
	{"i64small_gzip1_w1", "GZIP_1", false, 1},
	// BLANK = 37 (pixel 1 holds it) in 8-, 2- and 1-byte tiles.
	{"i64small_gzip1_blank", "GZIP_1", false, 8},
	{"i64small_gzip1_w2_blank", "GZIP_1", false, 2},
	{"i64small_gzip1_w1_blank", "GZIP_1", false, 1},
	// RICE_1 (BYTEPIX 4), HCOMPRESS_1 (lossless) and PLIO_1 at ZBITPIX 64,
	// built from libcfitsio's own coders (gen_fpack mk64): cfitsio decodes
	// them into C ints and widens them, so the values (both int32
	// extremes; PLIO 0 .. 2^24 - 1) stay in the int32 range. tileWidth 0:
	// not a GZIP tile.
	{"i64_rice", "RICE_1", false, 0},
	{"i64_hcomp", "HCOMPRESS_1", false, 0},
	{"i64_plio", "PLIO_1", false, 0},
}

// TestCompressedInt64Read: ZBITPIX 64 tile-compressed images are read
// exactly. The truth is libcfitsio's: for int64 its TLONG fits_read_img
// (the golden's i64: a 64-bit long on LP64), for uint64, whose TLONG read
// overflows (status 412, pinned), its stored values read with scaling
// off (the golden's raw) plus BZERO 2^63. Both must equal the same view of
// the twin, astropy 8.0.2's decode written as a plain image and read by
// cfitsio, so two independent decoders agree on every value. For 8-byte
// tiles cfitsio's TDOUBLE read returns memory it never wrote (pinned as
// unwrittenDouble, as for the 8-byte tiles of the 8/16/32-bit corpus). Values reach past 2^53 and both
// int64/uint64 extremes, so the check is on the integers, not on float64.
func TestCompressedInt64Read(t *testing.T) {
	cells := map[string]bool{}
	for _, c := range fz64Corpus {
		t.Run(c.name, func(t *testing.T) {
			fzName := "comp_fz_" + c.name + ".fits"
			g := testutil.LoadCfrefGolden(t, fzName)
			if g.NHDU != 2 || g.HDUs[1].Image == nil {
				t.Fatalf("cfitsio golden: %d HDUs", g.NHDU)
			}
			gh := &g.HDUs[1]
			im := gh.Image
			if cmp, _ := goldenCardValue(gh, "ZCMPTYPE"); cmp != c.cmp || im.Bitpix != 64 || !im.Compressed {
				t.Fatalf("golden ZCMPTYPE %q BITPIX %d compressed %v; case says %s ZBITPIX 64", cmp, im.Bitpix, im.Compressed, c.cmp)
			}
			if (im.Bzero == 9223372036854775808) != c.unsigned {
				t.Fatalf("golden BZERO %v; case unsigned %v", float64(im.Bzero), c.unsigned)
			}
			if c.tileWidth != 0 {
				if w := inflatedTileWidth(t, fzName); w != c.tileWidth {
					t.Fatalf("tile 1 inflates to %d bytes per pixel, case says %d", w, c.tileWidth)
				}
			}
			cells[c.cmp+map[bool]string{false: " 64 0", true: " 64 9223372036854775808"}[c.unsigned]] = true
			tg := testutil.LoadCfrefGolden(t, twinFile(t, c.name))
			twin := tg.HDUs[0].Image
			if twin == nil || (twin.Bitpix != 64 && !strings.HasSuffix(c.name, "_blank")) || twin.Compressed {
				t.Fatalf("twin golden is not a plain BITPIX 64 image")
			}
			n := int(im.Npix)
			var truthI []int64
			var truthU []uint64
			blank, hasBlank := goldenCardValue(gh, "BLANK")
			if hasBlank != strings.HasSuffix(c.name, "_blank") {
				t.Fatalf("golden BLANK present %v; case name %s", hasBlank, c.name)
			}
			var truthNull []bool
			if !c.unsigned {
				fz := im.I64.Values(t)
				if len(fz) != n {
					t.Fatalf("golden i64 %d values, npix %d", len(fz), n)
				}
				if hasBlank {
					// astropy decodes an integer image with BLANK to
					// float64 with NaN for the nulls: its twin is BITPIX
					// -64, equal where cfitsio has a value (all below
					// 2^53) and NaN exactly where cfitsio flags a null.
					tw := tg.HDUs[0].Image
					if tw.Bitpix != -64 {
						t.Fatalf("twin BITPIX %d, want astropy's float64 decode", tw.Bitpix)
					}
					tv := tw.Phys.Values(t)
					nulls := 0
					for k := range fz {
						tnan := !tv[k].Valid || math.IsNaN(tv[k].V)
						if !fz[k].Valid {
							nulls++
						}
						if tnan != !fz[k].Valid || (fz[k].Valid && tv[k].V != float64(fz[k].V)) {
							t.Fatalf("pixel %d: cfitsio TLONG %+v, astropy decode %+v", k, fz[k], tv[k])
						}
						truthI = append(truthI, fz[k].V)
						truthNull = append(truthNull, !fz[k].Valid)
					}
					if nulls == 0 {
						t.Fatalf("BLANK %s flags no pixel: the case proves nothing", blank)
					}
				} else {
					tw := twin.I64.Values(t)
					for k := range fz {
						if !fz[k].Valid || fz[k] != tw[k] {
							t.Fatalf("pixel %d: cfitsio TLONG %+v, astropy decode %+v", k, fz[k], tw[k])
						}
						truthI = append(truthI, fz[k].V)
						truthNull = append(truthNull, false)
					}
				}
			} else {
				if len(gh.Errors) == 0 || !strings.Contains(gh.Errors[0], "read_img TLONG: status 412") {
					t.Fatalf("golden errors %v: cfitsio's uint64 TLONG read no longer overflows; use it", gh.Errors)
				}
				fz, tw := im.Raw.Values(t), twin.Raw.Values(t)
				if len(fz) != n || len(tw) != n {
					t.Fatalf("golden raw %d values, twin %d, npix %d", len(fz), len(tw), n)
				}
				for k := range fz {
					if !fz[k].Valid || fz[k] != tw[k] {
						t.Fatalf("pixel %d: cfitsio stored %+v, astropy decode stored %+v", k, fz[k], tw[k])
					}
					truthU = append(truthU, uint64(fz[k].V)^(1<<63)) // + 2^63 (mod 2^64)
					truthNull = append(truthNull, false)
				}
			}
			if c.tileWidth == 8 {
				for k, p := range im.Phys.Values(t) {
					if !p.Valid || math.Float64bits(p.V) != unwrittenDouble {
						t.Fatalf("cfitsio TDOUBLE pixel %d = %+v; the 8-byte-tile note says cfitsio leaves it unwritten: update it", k, p)
					}
				}
			}

			f, h := openCompressed(t, filepathTestdata(t, fzName))
			defer f.Close()
			if h.BITPIX() != 64 || h.CompressionType() != c.cmp {
				t.Fatalf("BITPIX %d CompressionType %s", h.BITPIX(), h.CompressionType())
			}
			if !c.unsigned {
				got, err := ReadPixelsCompressed[int64](h)
				if err != nil {
					t.Fatalf("ReadPixelsCompressed[int64]: %v", err)
				}
				var bv int64
				if hasBlank {
					if bv, err = strconv.ParseInt(blank, 10, 64); err != nil {
						t.Fatal(err)
					}
				}
				for k := range truthI {
					// A null pixel's unmasked value is its stored value:
					// the BLANK itself (it fits every tile width here).
					if (truthNull[k] && got[k] != bv) || (!truthNull[k] && got[k] != truthI[k]) {
						t.Fatalf("pixel %d: %d, cfitsio %d (null %v, BLANK %s)", k, got[k], truthI[k], truthNull[k], blank)
					}
				}
				mv, mm, err := ReadPixelsCompressedMasked[int64](h)
				if err != nil {
					t.Fatalf("ReadPixelsCompressedMasked[int64]: %v", err)
				}
				for k := range truthI {
					if mm[k] == truthNull[k] || (mm[k] && mv[k] != truthI[k]) {
						t.Fatalf("masked int64 pixel %d: %d valid %v; cfitsio %d null %v", k, mv[k], mm[k], truthI[k], truthNull[k])
					}
				}
				if _, err := ReadPixelsCompressed[uint64](h); !errors.Is(err, ErrTypeMismatch) {
					t.Errorf("ReadPixelsCompressed[uint64] on int64: err %v, want ErrTypeMismatch", err)
				}
			} else {
				got, err := ReadPixelsCompressed[uint64](h)
				if err != nil {
					t.Fatalf("ReadPixelsCompressed[uint64]: %v", err)
				}
				for k := range truthU {
					if got[k] != truthU[k] {
						t.Fatalf("pixel %d: %d, cfitsio stored + 2^63 %d", k, got[k], truthU[k])
					}
				}
				if _, err := ReadPixelsCompressed[int64](h); !errors.Is(err, ErrTypeMismatch) {
					t.Errorf("ReadPixelsCompressed[int64] on uint64: err %v, want ErrTypeMismatch", err)
				}
			}
			if _, err := ReadPixelsCompressed[int32](h); !errors.Is(err, ErrTypeMismatch) {
				t.Errorf("ReadPixelsCompressed[int32]: err %v, want ErrTypeMismatch", err)
			}
			f64, mask, err := ReadPixelsCompressedMasked[float64](h)
			if err != nil {
				t.Fatalf("ReadPixelsCompressedMasked[float64]: %v", err)
			}
			for k := 0; k < n; k++ {
				var want float64
				if c.unsigned {
					want = float64(truthU[k])
				} else {
					want = float64(truthI[k])
				}
				if truthNull[k] {
					if mask[k] || f64[k] != 0 {
						t.Fatalf("masked float64 pixel %d: %v valid %v; cfitsio flags it null", k, f64[k], mask[k])
					}
					continue
				}
				if !mask[k] || f64[k] != want || math.IsNaN(f64[k]) {
					t.Fatalf("masked float64 pixel %d: %v valid %v, want %v", k, f64[k], mask[k], want)
				}
			}
		})
	}
	for _, cell := range []string{"GZIP_1 64 0", "GZIP_2 64 0", "GZIP_1 64 9223372036854775808", "GZIP_2 64 9223372036854775808", "RICE_1 64 0", "HCOMPRESS_1 64 0", "PLIO_1 64 0"} {
		if !cells[cell] {
			t.Errorf("ZBITPIX 64 cell %q has no sample", cell)
		}
	}
}
