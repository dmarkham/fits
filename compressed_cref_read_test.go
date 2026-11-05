package fits

import (
	"encoding/binary"
	"encoding/json"
	"errors"
	"math"
	"math/big"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/compress"
	"github.com/dmarkham/fits/internal/testutil"
)

// The writer-parity golden set (compress/testdata/cref/golden/img) holds,
// for every case, a file written by libcfitsio 4.7.0's fits_img_compress
// (<case>.fits), its input pixels (<case>.in.bin) and, for lossy cases,
// cfitsio's own decode of that file (<case>.dec.bin: fits_read_img
// TDOUBLE, nulls as NaN; cref.c). For lossless integer cases cref.c checked
// that cfitsio's decode equals the input ("lossless_verified"), so the
// input is cfitsio's decode.
//
// TestReaderDecodesCfitsioWrittenGoldens ties OUR READER to those files:
// every cfitsio-written image is read with ReadPixelsCompressed (float64
// and the image's native Go type) and ReadPixelsCompressedMasked and must
// equal cfitsio's decode exactly. For ZBITPIX -32 the logical pixel is a
// float32, so the truth is float32(cfitsio TDOUBLE) and the float64 read
// is that float32 widened (the documented float64-read rule).
//
// Among the cases: PLIO_1 at ZBITPIX 8/16/32 and the uint16 convention
// (plio_u16: cfitsio adds 32768 before PLIO and subtracts it on read), the
// dither random-sequence wrap inside one tile (rice_f32_wrap_100x100,
// rice_f32_wrapin_z10000, rice_f32_wrapin_nan_z10000: 10000-pixel tiles),
// ZDITHER0 9000/10000, IEEE specials, NaN nulls, HCOMPRESS lossy integer
// scales, 1-D and 3-D images.
func TestReaderDecodesCfitsioWrittenGoldens(t *testing.T) {
	ran := 0
	for _, g := range imgGoldenCases(t) {
		if g.Status != 0 {
			continue // cfitsio refused to write it: no file
		}
		ran++
		t.Run(g.Name, func(t *testing.T) {
			f, h := openCompressed(t, filepath.Join(imgGoldenDir, g.Name+".fits"))
			defer f.Close()

			got, err := ReadPixelsCompressed[float64](h)
			if err != nil {
				t.Fatalf("ReadPixelsCompressed[float64]: %v", err)
			}
			want := crefImgTruth(t, g, got)
			if d := diffPix(got, want, nil); d.n > 0 {
				t.Errorf("ReadPixelsCompressed[float64] vs cfitsio: %v", d)
			}

			name, nat, err := readCrefNative(h, g.Type)
			if err != nil {
				t.Fatalf("ReadPixelsCompressed[%s]: %v", name, err)
			}
			if d := diffPix(nat, want, nil); d.n > 0 {
				t.Errorf("ReadPixelsCompressed[%s] vs cfitsio: %v", name, d)
			}

			mv, mask, err := ReadPixelsCompressedMasked[float64](h)
			if err != nil {
				t.Fatalf("ReadPixelsCompressedMasked: %v", err)
			}
			bad := 0
			for i := range want {
				null := math.IsNaN(want[i])
				if mask[i] == null || (!null && !sameFloat(mv[i], want[i])) {
					if bad < 3 {
						t.Errorf("masked pixel %d: value %v valid %v, cfitsio %v", i, mv[i], mask[i], want[i])
					}
					bad++
				}
			}
			if bad > 0 {
				t.Errorf("ReadPixelsCompressedMasked: %d/%d pixels differ from cfitsio", bad, len(want))
			}

			// ReadFloat32 on the camera formats is Siril's
			// float32(u) * c with c the Siril literal (src/core/siril.h);
			// u is cfitsio's decoded value.
			if c, ok := map[string]float32{"u8": .0039215686, "u16": .000015259022}[g.Type]; ok {
				rf, err := ReadFloat32(h)
				if err != nil {
					t.Fatalf("ReadFloat32: %v", err)
				}
				for i := range want {
					if w := float32(want[i]) * c; rf[i] != w {
						t.Errorf("ReadFloat32 pixel %d = %v, want %v (cfitsio value %v)", i, rf[i], w, want[i])
						break
					}
				}
			}
		})
	}
	// 96 MANIFEST cases minus the 5 cfitsio refusals.
	if ran != 91 {
		t.Fatalf("decoded %d cfitsio-written goldens, want 91", ran)
	}
}

// crefImgTruth is cfitsio's decode of golden img/<case>.fits: for lossy
// cases the golden's exact digest of it (float32-rounded for f32 images),
// against which our float64 reading got is verified exactly before it is
// returned as the truth (decTruth); for lossless integer cases the
// regenerated, cfitsio-verified input.
func crefImgTruth(t *testing.T, g imgGolden, got []float64) []float64 {
	t.Helper()
	n := g.npix()
	if g.Type == "f32" || g.HScale != 0 {
		return g.decTruth(t, "ReadPixelsCompressed[float64]", got)
	}
	if !g.LosslessVerified {
		t.Fatalf("%s: lossless case without cfitsio's lossless_verified", g.Name)
	}
	in := g.readInput(t)
	out := make([]float64, n)
	for i := range out {
		switch g.Type {
		case "u8":
			out[i] = float64(in[i])
		case "s8":
			out[i] = float64(int8(in[i]))
		case "i16":
			out[i] = float64(int16(binary.LittleEndian.Uint16(in[2*i:])))
		case "u16":
			out[i] = float64(binary.LittleEndian.Uint16(in[2*i:]))
		case "i32":
			out[i] = float64(int32(binary.LittleEndian.Uint32(in[4*i:])))
		case "u32":
			out[i] = float64(binary.LittleEndian.Uint32(in[4*i:]))
		default:
			t.Fatalf("%s: unknown type %q", g.Name, g.Type)
		}
	}
	return out
}

// readCrefNative reads h in the Go type of the golden's input pixels.
func readCrefNative(h *CompressedImageHDU, typ string) (string, []float64, error) {
	switch typ {
	case "u8":
		v, err := ReadPixelsCompressed[uint8](h)
		return "uint8", toF64(v), err
	case "s8":
		v, err := ReadPixelsCompressed[int8](h)
		return "int8", toF64(v), err
	case "i16":
		v, err := ReadPixelsCompressed[int16](h)
		return "int16", toF64(v), err
	case "u16":
		v, err := ReadPixelsCompressed[uint16](h)
		return "uint16", toF64(v), err
	case "i32":
		v, err := ReadPixelsCompressed[int32](h)
		return "int32", toF64(v), err
	case "u32":
		v, err := ReadPixelsCompressed[uint32](h)
		return "uint32", toF64(v), err
	}
	v, err := ReadPixelsCompressed[float32](h)
	return "float32", toF64(v), err
}

// compressedIntTarget reads a compressed image as one Go integer type.
type compressedIntTarget struct {
	name   string
	lo, hi *big.Float
	read   func(*CompressedImageHDU) ([]float64, []bool, error)
}

func cTarget[T int8 | uint8 | int16 | uint16 | int32 | uint32 | int64 | uint64](name string, lo, hi *big.Float) compressedIntTarget {
	return compressedIntTarget{name: name, lo: lo, hi: hi,
		read: func(h *CompressedImageHDU) ([]float64, []bool, error) {
			v, err := ReadPixelsCompressed[T](h)
			if err != nil {
				return nil, nil, err
			}
			_, mask, err := ReadPixelsCompressedMasked[T](h)
			if err != nil {
				return nil, nil, err
			}
			return toF64(v), mask, nil
		}}
}

var compressedIntTargets = []compressedIntTarget{
	cTarget[int8]("int8", big.NewFloat(math.MinInt8), big.NewFloat(math.MaxInt8)),
	cTarget[uint8]("uint8", big.NewFloat(0), big.NewFloat(math.MaxUint8)),
	cTarget[int16]("int16", big.NewFloat(math.MinInt16), big.NewFloat(math.MaxInt16)),
	cTarget[uint16]("uint16", big.NewFloat(0), big.NewFloat(math.MaxUint16)),
	cTarget[int32]("int32", big.NewFloat(math.MinInt32), big.NewFloat(math.MaxInt32)),
	cTarget[uint32]("uint32", big.NewFloat(0), big.NewFloat(math.MaxUint32)),
	cTarget[int64]("int64", new(big.Float).SetInt64(math.MinInt64), new(big.Float).SetInt64(math.MaxInt64)),
	cTarget[uint64]("uint64", big.NewFloat(0), new(big.Float).SetUint64(math.MaxUint64)),
}

// TestReadPixelsCompressedIntegerTargets applies the plain-image rule to
// tile-compressed images: an integer T is accepted only when ZBITPIX is an
// integer type, BSCALE and BZERO are integral, and T holds the physical
// value of every possible stored value; otherwise ReadPixelsCompressed
// (and the masked read) return a *TypeMismatchError wrapping
// ErrTypeMismatch. Accepted reads must equal libcfitsio's decode of the
// file (testdata/cfref golden) on every non-null pixel. Every case of the
// compressed read corpus is read as every Go integer type.
func TestReadPixelsCompressedIntegerTargets(t *testing.T) {
	for _, c := range fzCorpus {
		t.Run(c.name, func(t *testing.T) {
			if checkCorpusScalingRefused(t, c) {
				return
			}
			fzName := "comp_fz_" + c.name + ".fits"
			want, null := caseTruth(t, c, false)
			f, h := openCompressed(t, filepathTestdata(t, fzName))
			defer f.Close()
			skipRow := map[int]bool{}
			for _, r := range c.rawRows {
				skipRow[r] = true // cfitsio's fits_read_img leaves these unscaled (see fzCase.rawRows)
			}
			nx := int(h.Shape()[0])
			bscale, bzero := bscaleOf(t, h), bzeroOf(t, h)
			var lo, hi *big.Float
			integral := false
			if c.zbitpix > 0 {
				lo, hi, integral = physRange(c.zbitpix, bscale, bzero)
			}
			accepted := 0
			for _, tg := range compressedIntTargets {
				allowed := integral && lo.Cmp(tg.lo) >= 0 && hi.Cmp(tg.hi) <= 0
				got, mask, err := tg.read(h)
				if !allowed {
					if !errors.Is(err, ErrTypeMismatch) {
						t.Errorf("ReadPixelsCompressed[%s] (ZBITPIX %d, BSCALE %g, BZERO %g): err %v, want ErrTypeMismatch",
							tg.name, c.zbitpix, bscale, bzero, err)
						continue
					}
					var tm *TypeMismatchError
					if !errors.As(err, &tm) || tm.Requested != tg.name || tm.BITPIX != c.zbitpix || !tm.Lossy || tm.BSCALE != bscale || tm.BZERO != bzero {
						t.Errorf("ReadPixelsCompressed[%s]: TypeMismatchError %+v, want Requested %s BITPIX %d Lossy BSCALE %g BZERO %g",
							tg.name, tm, tg.name, c.zbitpix, bscale, bzero)
					}
					continue
				}
				accepted++
				if err != nil {
					t.Errorf("ReadPixelsCompressed[%s]: %v (T holds every physical value)", tg.name, err)
					continue
				}
				bad := 0
				for k := range want {
					if skipRow[k/nx] {
						continue
					}
					if mask[k] == null[k] {
						t.Errorf("ReadPixelsCompressedMasked[%s] pixel %d: valid %v, cfitsio null %v", tg.name, k, mask[k], null[k])
						break
					}
					if !null[k] && got[k] != want[k] {
						if bad == 0 {
							t.Errorf("ReadPixelsCompressed[%s] pixel %d = %v, cfitsio %v", tg.name, k, got[k], want[k])
						}
						bad++
					}
				}
				if bad > 0 {
					t.Errorf("ReadPixelsCompressed[%s]: %d/%d pixels differ from cfitsio", tg.name, bad, len(want))
				}
			}
			if c.zbitpix > 0 && integral && accepted == 0 {
				t.Errorf("no integer type accepted for integral scaling")
			}
		})
	}
}

func filepathTestdata(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join("testdata", name)
}

// TestZCMPTYPEComparedLikeCfitsio: cfitsio recognises ZCMPTYPE with a
// case-sensitive comparison (imcompress.c:5279-5292, FSTRCMP), so
// 'RICE_ONE' is RICE_1 (comp_fz_u16_rice_riceone, in the read corpus) and
// 'rice_1' is refused: libcfitsio cannot move to that HDU (the committed
// cfref golden of the sample stops at the primary HDU) and funpack 1.7.0 (CFITSIO 4.7.0)
// fails with status 414 "Unknown image compression type". The library
// refuses it too, with ErrUnknownCompression, instead of decoding it.
func TestZCMPTYPEComparedLikeCfitsio(t *testing.T) {
	const sample = "testdata/comp_bad/u16_rice_lowercase.fits"
	g := loadGoldenFile(t, "testdata/comp_bad/golden/u16_rice_lowercase.fits.json")
	if g.NHDU != 1 {
		t.Fatalf("cfitsio golden has %d HDUs; it should refuse the 'rice_1' HDU", g.NHDU)
	}
	f, err := Open(sample)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	hdu, err := f.HDU(1)
	if err != nil {
		t.Fatalf("HDU(1): %v (the HDU is framed correctly; only its ZCMPTYPE is refused)", err)
	}
	h, ok := hdu.(*CompressedImageHDU)
	if !ok {
		t.Fatalf("HDU 1 is %T", hdu)
	}
	if _, err := ReadPixelsCompressed[float64](h); !errors.Is(err, compress.ErrUnknownCompression) {
		t.Errorf("ReadPixelsCompressed on ZCMPTYPE 'rice_1': err %v, want ErrUnknownCompression (cfitsio refuses it)", err)
	}
	if _, err := ReadFloat32(h); !errors.Is(err, compress.ErrUnknownCompression) {
		t.Errorf("ReadFloat32 on ZCMPTYPE 'rice_1': err %v, want ErrUnknownCompression", err)
	}
}

// TestCfnullGoldensCoverImagingSamples: every tile-compressed imaging
// sample (testdata/imaging/*.fz) has its cfitsio golden
// (testdata/imaging/cfnull/golden: fits_read_pixnull TDOUBLE, null flags)
// and every golden has its sample, and every golden is checked here
// directly: ReadPixelsCompressedMasked[float64] equals cfitsio's physical
// values and null flags exactly, including the files whose ReadFloat32
// rule is "error" (BSCALE = 0: cfitsio silently returns zeros, which the
// generic reader matches and ReadFloat32 refuses, pinned by
// TestReadFloat32_ExpectedMatrix). The samples whose scaling cards the
// scaling-card rule refuses (readfloat32_expected.json invalid_card:
// comp_card_i16_bscale_str, comp_card_i16_bscale_hierarch) instead fail
// every pixel read with ErrInvalidScaling (the golden records cfitsio's
// silent reading).
func TestCfnullGoldensCoverImagingSamples(t *testing.T) {
	exp := loadRF32Expected(t)
	refused := 0
	samples, _ := filepath.Glob(imagingPath("*.fz"))
	goldens, _ := filepath.Glob(imagingPath(filepath.Join("cfnull", "golden", "*.json")))
	if len(samples) < 12 || len(goldens) != len(samples) {
		t.Fatalf("%d .fz samples, %d cfnull goldens", len(samples), len(goldens))
	}
	for _, s := range samples {
		name := filepath.Base(s)
		t.Run(name, func(t *testing.T) {
			b, err := os.ReadFile(imagingPath(filepath.Join("cfnull", "golden", name+".json")))
			if err != nil {
				t.Fatalf("committed golden missing (failure, not a skip): %v", err)
			}
			var g struct {
				File     string     `json:"file"`
				Cfitsio  string     `json:"cfitsio"`
				ImageHDU int        `json:"image_hdu"`
				Status   int        `json:"status"`
				Npix     int        `json:"npix"`
				Phys     []*float64 `json:"phys"`
			}
			if err := json.Unmarshal(b, &g); err != nil {
				t.Fatal(err)
			}
			if g.File != name || g.Cfitsio != "4.070" || g.Status != 0 || g.ImageHDU != 2 || len(g.Phys) != g.Npix {
				t.Fatalf("golden header %+v", g)
			}
			f, h := openCompressed(t, s)
			defer f.Close()
			if card := exp[name].InvalidCard; card != "" {
				refused++
				_, _, err := ReadPixelsCompressedMasked[float64](h)
				_, err2 := ReadPixelsCompressed[float64](h)
				_, err3 := h.BSCALE()
				if err := h.Validate(); err != nil {
					t.Errorf("Validate: %v; the refused card is not a structural error", err)
				}
				for _, e := range []error{err, err2, err3} {
					if !errors.Is(e, ErrInvalidScaling) || !strings.Contains(e.Error(), card) {
						t.Errorf("%v, want ErrInvalidScaling naming %s", e, card)
					}
				}
				return
			}
			got, mask, err := ReadPixelsCompressedMasked[float64](h)
			if err != nil {
				t.Fatalf("ReadPixelsCompressedMasked: %v (cfitsio reads it)", err)
			}
			if len(got) != g.Npix {
				t.Fatalf("%d pixels, cfitsio %d", len(got), g.Npix)
			}
			for i, p := range g.Phys {
				if (p == nil) == mask[i] || (p != nil && got[i] != *p) {
					t.Errorf("pixel %d: %v valid %v, cfitsio %v", i, got[i], mask[i], p)
				}
			}
		})
	}
	for _, gp := range goldens {
		if _, err := os.Stat(imagingPath(strings.TrimSuffix(filepath.Base(gp), ".json"))); err != nil {
			t.Errorf("golden %s has no sample", filepath.Base(gp))
		}
	}
	if refused != 2 {
		t.Errorf("%d samples refused by the scaling-card rule, want 2 (comp_card_i16_bscale_str, _hierarch)", refused)
	}
}

// TestCompressedBadSamplesRefused: the byte-patched / key-edited files in
// testdata/comp_bad (gen_fpack.sh) hold layouts the library refuses with
// an error rather than decode. Each guard was removable without any test
// failing before these samples existed. The committed cfitsio goldens
// record what libcfitsio 4.7.0 does with each:
//
//   - ZTILE1 = 0 and a missing COMPRESSED_DATA column (TTYPE1 renamed):
//     cfitsio will not move to the HDU (the golden stops at the primary).
//   - PLIO_1 with BZERO = 32768 at ZBITPIX 32 and 8: cfitsio's reader
//     applies its uint16 shift whatever ZBITPIX is while its writer adds
//     it only for 16-bit images, so there is no consistent reading; the
//     golden shows cfitsio returning values, the library refuses.
//
// cfitsio is "hdu" (cfitsio will not move to the HDU: the golden stops at
// the primary), "read" (cfitsio reads the HDU's pixels without an error)
// or a status ("415": cfitsio moves to the HDU and every pixel read fails
// with that status).
func TestCompressedBadSamplesRefused(t *testing.T) {
	cases := []struct {
		file, want string
		cfitsio    string
		noNull     bool // the cfitsio golden flags no pixel null
	}{
		{"u16_rice_ztile0.fits", "ZTILE1 = 0", "hdu", false},
		{"u16_rice_nocompdata.fits", "missing COMPRESSED_DATA column", "hdu", false},
		{"i32_plio_bzero32768.fits", "PLIO_1 with BZERO = 32768 on ZBITPIX 32", "read", false},
		{"u8_plio_bzero32768.fits", "PLIO_1 with BZERO = 32768 on ZBITPIX 8", "read", false},
		// Tile 3 (row 3) with an empty COMPRESSED_DATA and no fallback
		// column: cfitsio fails every read (status 415, NO_COMPRESSED_TILE);
		// read as zeros it would drop the tile silently.
		{"u16_rice_emptytile.fits", "tile 2 has no payload", "415", false},
		// An integer tile in GZIP_COMPRESSED_DATA: cfitsio writes that
		// column only for float tiles and fails the read (status 414,
		// "incompatible data type in gzipped floating-point tile").
		{"u16_rice_gzipcol.fits", "GZIP_COMPRESSED_DATA on an integer image (ZBITPIX 16)", "414", false},
		// ZNAXIS2 patched 24 -> 23: 24 table rows for 23 tiles. cfitsio
		// will not move to the HDU.
		{"u16_rice_extrarow.fits", "tile count mismatch: 24 rows vs 23 tiles", "hdu", false},
		// UNCOMPRESSED_DATA tile 4 one value short: cfitsio reads 39 values
		// and leaves the 40th pixel of the row undefined, without an error;
		// the library refuses.
		{"u16_uncomp_short.fits", "tile 3 UNCOMPRESSED_DATA has 78 bytes, want 80", "read", false},
		// A GZIP_1 tile one byte longer than 2 x tilelen: cfitsio reads
		// the HDU but every pixel read fails (status 414, "uncompressed
		// tile has wrong size"); funpack refuses it
		// (TestCompressedBadTilesFunpackRefuses).
		{"u16_gzip1_extrabyte.fits", "inflates to", "414", false},
		// 4-byte tiles on a 16-bit image holding 70000 - 32768: cfitsio's
		// TDOUBLE read (the golden) returns it, a read as the image type
		// fails with status 412 and funpack refuses the file; refused.
		{"i16_gzip1_w4_overflow.fits", "does not fit ZBITPIX 16", "read", false},
		// ZNAXIS outside 1..6: cfitsio refuses to move to the HDU
		// (BAD_NAXIS, imcompress.c:5370-5379). 7 axes written by astropy
		// 8.0.2 (generate_compressed.py).
		{"u16_rice_znaxis0.fits", "ZNAXIS 0", "hdu", false},
		{"u16_rice_znaxism1.fits", "ZNAXIS -1", "hdu", false},
		{"gzip1_7d.fits", "ZNAXIS 7", "hdu", false},
		// One past each edge of the ZBITPIX range in wider tiles (the
		// edges themselves: comp_fz_*_edge). cfitsio's TDOUBLE/TLONG
		// reads return the value, a read as the image type fails (412),
		// funpack refuses (TestCompressedBadTilesFunpackRefuses).
		{"u8_gzip2_w2_lo.fits", "tile value -1 does not fit ZBITPIX 8", "read", false},
		{"u8_gzip2_w2_hi.fits", "tile value 256 does not fit ZBITPIX 8", "read", false},
		{"u16_gzip1_w4_lo.fits", "tile value -32769 does not fit ZBITPIX 16", "read", false},
		{"u16_gzip1_w4_hi.fits", "tile value 32768 does not fit ZBITPIX 16", "read", false},
		{"i32_gzip1_w8_lo.fits", "tile value -2147483649 does not fit ZBITPIX 32", "read", false},
		{"i32_gzip1_w8_hi.fits", "tile value 2147483648 does not fit ZBITPIX 32", "read", false},
		// BLANK outside the C int cfitsio reads it as: cfitsio's ffgky
		// TINT fails, so it sets no null value and flags nothing (golden);
		// astropy 8.0.2 refuses the file (TestCompressedBlankRangeAstropy).
		// A truncated BLANK would flag the zeros (2^32) or 7 (2^40 + 7).
		{"i32_rice_blank2p32.fits", "BLANK = 4294967296 is outside the BITPIX 32 range", "read", true},
		{"i64_gzip1_blank2p40.fits", "BLANK 1099511627783 is outside the int32 range", "read", true},
		// ZBLANK 2^33 next to BLANK -32768: cfitsio's ZBLANK read fails
		// and it uses BLANK instead (golden: the -32768 pixels null).
		{"i16_rice_zblank2p33.fits", "ZBLANK 8589934592 is outside the int32 range", "read", false},
		// A per-tile 1K ZBLANK column holding 2^32: cfitsio fails every
		// read reading the tile's null value (status 412).
		{"i16_rice_zblankcol1k.fits", "tile 0: ZBLANK column value 4294967296 is outside the int32 range", "412", false},
		// Quantized float GZIP tiles of 8 bytes: cfitsio has no 8-byte
		// branch for quantized tiles; cfitsio returns memory it never
		// wrote (golden: unwrittenDouble), funpack and astropy 8.0.2 zeros,
		// for an image that is not zero.
		{"f32_gzip1_q_w8.fits", "tile 0: quantized GZIP_1 tile inflates to 8 bytes per pixel", "read", false},
		// Round 2 fix pass. Lossy HCOMPRESS_1 of an 8-bit image saturated at
		// 255 reconstructs pixels above it: cfitsio's TDOUBLE read returns
		// them (golden), funpack fails (412), its TBYTE read clips them
		// with NUM_OVERFLOW; refused as any other tile value the ZBITPIX
		// type cannot hold.
		{"u8sat_hcomp_s2.fits", "does not fit ZBITPIX 8", "read", false},
		// The 32-bit one: fits_hdecompress64 casts the reconstruction to a
		// C int, so cfitsio (golden) and funpack return the saturated
		// pixels wrapped to about -2^31; refused.
		{"i32sat_hcomp_s2.fits", "is outside the int32 range of HCOMPRESS_1 tile integers", "read", false},
		// ZBITPIX 64 HCOMPRESS_1 with pixel 5 one past either end of the
		// int32 range: cfitsio's (int) cast reads them as -2^31 and 2^31 - 1
		// (golden pixel 5); refused.
		{"i64_hcomp_2p31.fits", "HCOMPRESS_1 tile value 2147483648 is outside the int32 range", "read", false},
		{"i64_hcomp_m2p31m1.fits", "HCOMPRESS_1 tile value -2147483649 is outside the int32 range", "read", false},
		// PLIO_1 values the ZBITPIX type cannot hold: cfitsio's TDOUBLE
		// read returns them, its read as the image type and funpack fail
		// (412). uint16 PLIO with BZERO 0 (values to 65535), the 0..2^24-1
		// values as ZBITPIX 8, and as ZBITPIX 16 with BZERO 32768 (value -
		// 32768 above 32767).
		{"u16_plio_bzero0.fits", "does not fit ZBITPIX 16", "read", false},
		{"i32_plio_zbitpix8.fits", "does not fit ZBITPIX 8", "read", false},
		{"i32_plio_zbitpix16_bzero32768.fits", "does not fit ZBITPIX 16", "read", false},
		// One past each end of the int32 range for BLANK (the ends
		// themselves: comp_fz_i32_gzip1_blankmin/_blankmax): cfitsio's C
		// int read of the keyword fails and it flags nothing (golden).
		{"i32_gzip1_blank2p31.fits", "BLANK = 2147483648 is outside the BITPIX 32 range", "read", true},
		{"i32_gzip1_blankm2p31m1.fits", "BLANK = -2147483649 is outside the BITPIX 32 range", "read", true},
		// ... and for a 1K ZBLANK column (-2^31 itself:
		// comp_fz_i32_gzip1_zblankcolmin): cfitsio fails every read (412).
		{"i32_gzip1_zblankcol2p31.fits", "tile 0: ZBLANK column value 2147483648 is outside the int32 range", "412", false},
		{"i32_gzip1_zblankcolm2p31m1.fits", "tile 0: ZBLANK column value -2147483649 is outside the int32 range", "412", false},
		// A ZBITPIX that is no FITS BITPIX: cfitsio reads the RICE tiles
		// anyway (golden); refused (it used to panic).
		{"u16_rice_zbitpix12.fits", "ZBITPIX 12 is not a valid BITPIX", "read", false},
		{"u16_rice_zbitpix0.fits", "ZBITPIX 0 is not a valid BITPIX", "read", false},
		// ZNAXIS1 = -1, ZTILE1 = 1, ZNAXIS2 = -24 (a tile count that
		// matches the 24 rows): cfitsio refuses the HDU (status 414); it
		// used to panic here.
		{"u16_rice_znaxisneg.fits", "ZNAXIS1 = -1", "hdu", false},
	}
	files, _ := filepath.Glob(filepath.Join("testdata", "comp_bad", "*.fits"))
	// + u16_rice_lowercase.fits (TestZCMPTYPEComparedLikeCfitsio),
	// ztable_tbl500.fits (TestTiledTableCompressionRefused) and
	// i16_gzip1_w4_nullout.fits, i16_gzip1_w4_nullout2.fits
	// (TestCompressedNullOutsideZBITPIX).
	if len(files) != len(cases)+4 {
		t.Fatalf("testdata/comp_bad has %d files, the tests name %d", len(files), len(cases)+4)
	}
	for _, c := range cases {
		t.Run(c.file, func(t *testing.T) {
			g := loadGoldenFile(t, "testdata/comp_bad/golden/"+c.file+".json")
			if g.Cfitsio != "4.7.0" {
				t.Fatalf("golden made with cfitsio %s", g.Cfitsio)
			}
			switch c.cfitsio {
			case "hdu":
				if g.NHDU != 1 {
					t.Fatalf("cfitsio golden has %d HDUs; cfitsio should refuse HDU 1", g.NHDU)
				}
			case "read":
				if g.NHDU != 2 || g.HDUs[1].Image == nil || !g.HDUs[1].Image.Phys.Present() || len(g.HDUs[1].Errors) != 0 {
					t.Fatalf("cfitsio golden: %d HDUs, errors %v; want cfitsio reading HDU 1's pixels", g.NHDU, g.HDUs[1].Errors)
				}
				if c.noNull && (g.HDUs[1].Image.Phys.NullCount() != 0 || g.HDUs[1].Image.I64.NullCount() != 0) {
					t.Fatalf("cfitsio golden flags %d pixels null; the case says cfitsio sets no null value", g.HDUs[1].Image.Phys.NullCount())
				}
			default:
				if g.NHDU != 2 || g.HDUs[1].Image == nil || g.HDUs[1].Image.Phys.Present() || len(g.HDUs[1].Errors) == 0 || !strings.Contains(g.HDUs[1].Errors[0], "status "+c.cfitsio) {
					t.Fatalf("cfitsio golden: %d HDUs, errors %v; want every pixel read to fail with status %s", g.NHDU, g.HDUs[1].Errors, c.cfitsio)
				}
			}
			f, err := Open(filepath.Join("testdata", "comp_bad", c.file))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			hdu, err := f.HDU(1)
			if err != nil {
				if !strings.Contains(err.Error(), c.want) {
					t.Fatalf("HDU(1): %v, want an error mentioning %q", err, c.want)
				}
				return
			}
			h, ok := hdu.(*CompressedImageHDU)
			if !ok {
				t.Fatalf("HDU 1 is %T", hdu)
			}
			_ = h.Shape() // must not panic on a bad ZNAXIS
			if _, err := ReadPixelsCompressed[float64](h); err == nil || !strings.Contains(err.Error(), c.want) {
				t.Errorf("ReadPixelsCompressed: err %v, want one mentioning %q", err, c.want)
			}
			f32want := c.want
			if h.BITPIX() == 64 {
				// ReadFloat32 refuses every ZBITPIX 64 image first.
				f32want = "ZBITPIX 64 tile-compressed images are not supported"
			}
			if _, err := ReadFloat32(h); err == nil || !strings.Contains(err.Error(), f32want) {
				t.Errorf("ReadFloat32: err %v, want one mentioning %q", err, f32want)
			}
		})
	}
}

// TestCompressedBadTilesFunpackRefuses: funpack 1.7.0 (CFITSIO 4.7.0)
// refuses the bad-tile samples the library refuses, with the cfitsio
// status the library's error names (412 overflow, 414 wrong tile size),
// and the null value outside ZBITPIX that the library's unmasked read
// refuses.
func TestCompressedBadTilesFunpackRefuses(t *testing.T) {
	for file, status := range map[string]string{
		"i16_gzip1_w4_overflow.fits":         "FITSIO status = 412",
		"u16_gzip1_extrabyte.fits":           "FITSIO status = 414",
		"u8_gzip2_w2_lo.fits":                "FITSIO status = 412",
		"u8_gzip2_w2_hi.fits":                "FITSIO status = 412",
		"u16_gzip1_w4_lo.fits":               "FITSIO status = 412",
		"u16_gzip1_w4_hi.fits":               "FITSIO status = 412",
		"i32_gzip1_w8_lo.fits":               "FITSIO status = 412",
		"i32_gzip1_w8_hi.fits":               "FITSIO status = 412",
		"i16_gzip1_w4_nullout.fits":          "FITSIO status = 412",
		"i16_rice_zblankcol1k.fits":          "FITSIO status = 412",
		"i16_gzip1_w4_nullout2.fits":         "FITSIO status = 412",
		"u8sat_hcomp_s2.fits":                "FITSIO status = 412",
		"u16_plio_bzero0.fits":               "FITSIO status = 412",
		"i32_plio_zbitpix8.fits":             "FITSIO status = 412",
		"i32_plio_zbitpix16_bzero32768.fits": "FITSIO status = 412",
		"i32_gzip1_zblankcol2p31.fits":       "FITSIO status = 412",
		"i32_gzip1_zblankcolm2p31m1.fits":    "FITSIO status = 412",
	} {
		p := filepath.Join("testdata", "comp_bad", file)
		r := testutil.Run(t, testutil.ToolFunpack, testutil.Call{Args: []string{"-S", p}, Inputs: []string{p}})
		out := append(r.Stdout, r.Stderr...)
		if r.Exit == 0 || !strings.Contains(string(out), status) {
			t.Errorf("funpack -S %s: exit %d, output %q; want a failure with %q", file, r.Exit, out, status)
		}
	}
}

// TestTiledTableCompressionRefused: tiled TABLE compression (FITS 4.0
// §10.3, ZTABLE = T; testdata/comp_bad/ztable_tbl500.fits is fpack -table
// of a 500-row cfitsio-written table) is not supported. Read as a plain
// binary table it is wrong without any error: the cfitsio golden shows
// what a plain reader sees (NAXIS2 = 1 row of gzip/RICE bytes for each
// column, ZNAXIS2 = 500 the real row count; cfitsio's plain reader is
// no better). HDU(1) must refuse it with ErrTiledTableCompression.
func TestTiledTableCompressionRefused(t *testing.T) {
	const name = "ztable_tbl500.fits"
	g := loadGoldenFile(t, "testdata/comp_bad/golden/"+name+".json")
	if g.NHDU != 2 || g.HDUs[1].Table == nil || g.HDUs[1].Table.NRows != 1 {
		t.Fatalf("cfitsio golden: %d HDUs, want HDU 1 a 1-row binary table", g.NHDU)
	}
	if v, ok := goldenCardValue(&g.HDUs[1], "ZTABLE"); !ok || v != "T" {
		t.Fatalf("golden ZTABLE = %q (present %v)", v, ok)
	}
	if v, _ := goldenCardValue(&g.HDUs[1], "ZNAXIS2"); v != "500" {
		t.Fatalf("golden ZNAXIS2 = %q, want 500", v)
	}
	f, err := Open(filepath.Join("testdata", "comp_bad", name))
	if err != nil {
		t.Fatalf("Open: %v (the file is framed correctly)", err)
	}
	defer f.Close()
	h, err := f.HDU(1)
	if !errors.Is(err, ErrTiledTableCompression) {
		t.Fatalf("HDU(1) = %T, err %v; want ErrTiledTableCompression", h, err)
	}
	if _, err := f.HDU(0); err != nil {
		t.Errorf("HDU(0): %v", err)
	}
}
