package fits

import (
	"encoding/json"
	"errors"
	"fmt"
	"math"
	"math/big"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// physRange returns the exact physical range [lo, hi] of an integer BITPIX
// image with the given BSCALE/BZERO (as cfitsio reports them), and whether
// every physical value is an integer (BSCALE and BZERO integral). It is
// computed with math/big, independently of the library's int128 code.
func physRange(bitpix int, bscale, bzero float64) (lo, hi *big.Float, integral bool) {
	var rmin, rmax *big.Float
	switch bitpix {
	case 8:
		rmin, rmax = big.NewFloat(0), big.NewFloat(255)
	case 16:
		rmin, rmax = big.NewFloat(math.MinInt16), big.NewFloat(math.MaxInt16)
	case 32:
		rmin, rmax = big.NewFloat(math.MinInt32), big.NewFloat(math.MaxInt32)
	case 64:
		rmin, rmax = new(big.Float).SetInt64(math.MinInt64), new(big.Float).SetInt64(math.MaxInt64)
	}
	f := func(r *big.Float) *big.Float {
		v := new(big.Float).SetPrec(256).Mul(new(big.Float).SetPrec(256).Set(r), big.NewFloat(bscale))
		return v.Add(v, big.NewFloat(bzero))
	}
	a, b := f(rmin), f(rmax)
	if a.Cmp(b) > 0 {
		a, b = b, a
	}
	return a, b, bscale == math.Trunc(bscale) && bzero == math.Trunc(bzero)
}

// intTarget describes one Go integer target type.
type intTarget struct {
	name   string
	lo, hi *big.Float
	read   func(*ImageHDU) ([]*big.Int, error)
	can    func(*ImageHDU) bool
}

func target[T int8 | uint8 | int16 | uint16 | int32 | uint32 | int64 | uint64](name string, lo, hi *big.Float) intTarget {
	return intTarget{name: name, lo: lo, hi: hi,
		read: func(h *ImageHDU) ([]*big.Int, error) {
			v, err := ReadPixels[T](h)
			if err != nil {
				return nil, err
			}
			out := make([]*big.Int, len(v))
			var zero T
			for i := range v {
				if zero-1 > 0 { // unsigned
					out[i] = new(big.Int).SetUint64(uint64(v[i]))
				} else {
					out[i] = big.NewInt(int64(v[i]))
				}
			}
			return out, nil
		},
		can: func(h *ImageHDU) bool { return CanConvert[T](h) },
	}
}

var intTargets = []intTarget{
	target[int8]("int8", big.NewFloat(math.MinInt8), big.NewFloat(math.MaxInt8)),
	target[uint8]("uint8", big.NewFloat(0), big.NewFloat(math.MaxUint8)),
	target[int16]("int16", big.NewFloat(math.MinInt16), big.NewFloat(math.MaxInt16)),
	target[uint16]("uint16", big.NewFloat(0), big.NewFloat(math.MaxUint16)),
	target[int32]("int32", big.NewFloat(math.MinInt32), big.NewFloat(math.MaxInt32)),
	target[uint32]("uint32", big.NewFloat(0), big.NewFloat(math.MaxUint32)),
	target[int64]("int64", new(big.Float).SetInt64(math.MinInt64), new(big.Float).SetInt64(math.MaxInt64)),
	target[uint64]("uint64", big.NewFloat(0), new(big.Float).SetUint64(math.MaxUint64)),
}

// cfitsio image equivalent types (fitsio.h) and the Go type each names.
var equivGoType = map[int]string{8: "uint8", 10: "int8", 16: "int16", 20: "uint16", 32: "int32", 40: "uint32", 64: "int64", 80: "uint64"}

// parseHeaderFirst, when set (TestReadsAfterHeaderParseMatchCfitsio), is
// called on every HDU before the read checks, to load the full parsed
// header so BSCALE/BZERO/BLANK come from it instead of the raw cards.
var parseHeaderFirst func(t *testing.T, s imageSample, i int, img *ImageHDU)

// imageHDUsUnderTest yields every plain image HDU with pixels in the corpus
// that cfitsio could read, with its golden.
func forEachImage(t *testing.T, fn func(t *testing.T, s imageSample, i int, h *ImageHDU, g *testutil.CfrefHDU)) {
	for _, s := range allImageSamples(t) {
		f, err := Open(s.path)
		if err != nil {
			t.Errorf("%s: Open: %v (cfitsio opens it)", s.name, err)
			continue
		}
		for i := range s.golden.HDUs {
			g := &s.golden.HDUs[i]
			if g.Image == nil || g.Image.Compressed || g.Image.Npix == 0 || !g.Image.Phys.Present() {
				continue
			}
			hdu, err := f.HDU(i)
			if errors.Is(err, ErrUnsupportedExtension) {
				// framing/xtension_foreign.fits: cfitsio reads an unknown
				// extension type as an image; the library refuses it
				// (TestUnsupportedExtensionFramedAndRefused).
				continue
			}
			if err != nil {
				t.Errorf("%s HDU %d: %v", s.name, i, err)
				continue
			}
			img, ok := hdu.(*ImageHDU)
			if !ok {
				t.Errorf("%s HDU %d: %T, cfitsio says plain image", s.name, i, hdu)
				continue
			}
			if parseHeaderFirst != nil {
				parseHeaderFirst(t, s, i, img)
			}
			t.Run(fmt.Sprintf("%s/%d", s.name, i), func(t *testing.T) { fn(t, s, i, img, g) })
		}
		f.Close()
	}
}

// TestReadPixelsFloat64MatchesCfitsio: ReadPixels[float64] equals cfitsio's
// TDOUBLE reading (BSCALE/BZERO applied) bit for bit on every image of the
// corpus. A BLANK pixel of an integer image (cfitsio null) must read as
// its scaled stored value, which is cfitsio's TDOUBLE reading with null
// checking off (imgref pix "f64nn"), as the ReadPixels doc says.
func TestReadPixelsFloat64MatchesCfitsio(t *testing.T) {
	blanksChecked := 0
	forEachImage(t, func(t *testing.T, s imageSample, i int, img *ImageHDU, g *testutil.CfrefHDU) {
		got, err := ReadPixels[float64](img)
		if err != nil {
			t.Fatal(err)
		}
		if len(got) != g.Image.Phys.Len() {
			t.Fatalf("%d pixels, cfitsio %d", len(got), g.Image.Phys.Len())
		}
		var nn []testutil.NullFloat
		if g.Image.Bitpix > 0 && g.Image.Phys.NullCount() > 0 {
			p := s.pix.HDUs[i]
			if p.F64NNError != "" || !p.F64NN.Present() {
				t.Fatalf("image has BLANK pixels but the golden has no TDOUBLE reading without null checks (%s)", p.F64NNError)
			}
			nn = p.F64NN.Resolve(t, "cfitsio TDOUBLE without null checks", func() []testutil.NullFloat {
				v, err := ReadPixels[float64](img)
				if err != nil {
					return nil
				}
				return testutil.NullFloats(v, nil)
			})
		}
		for k, w := range g.Image.Phys.Resolve(t, "cfitsio TDOUBLE", func() []testutil.NullFloat { return cfitsioPhysView(t, img) }) {
			if p, ok := cfitsioFloatQuirk(t, img, k, w); ok {
				if !sameFloat(got[k], p) {
					t.Fatalf("pixel %d: %v, stored %v (cfitsio reports %v)", k, got[k], p, w)
				}
				continue
			}
			if !w.Valid {
				if g.Image.Bitpix < 0 {
					if !math.IsNaN(got[k]) {
						t.Fatalf("pixel %d: %v, cfitsio null (NaN)", k, got[k])
					}
					continue
				}
				if !sameFloat(got[k], nn[k].V) {
					t.Fatalf("BLANK pixel %d: %v, cfitsio's scaled stored value %v", k, got[k], nn[k].V)
				}
				blanksChecked++
				continue
			}
			if !sameFloat(got[k], w.V) {
				t.Fatalf("pixel %d: %v, cfitsio %v", k, got[k], w.V)
			}
		}
	})
	if blanksChecked == 0 {
		t.Error("no BLANK pixel of an integer image was compared with its scaled stored value")
	}
}

// TestReadPixelsFloat32MatchesCfitsio: ReadPixels[float32] equals cfitsio's
// TFLOAT reading bit for bit (for BITPIX 64 cfitsio rounds the integer to
// float once; rounding through float64 first differs). BLANK pixels of
// integer images read as their scaled stored value rounded to float32.
func TestReadPixelsFloat32MatchesCfitsio(t *testing.T) {
	forEachImage(t, func(t *testing.T, s imageSample, i int, img *ImageHDU, g *testutil.CfrefHDU) {
		p := s.pix.HDUs[i]
		if p.F32Error != "" || !p.F32.Present() {
			t.Fatalf("golden has no TFLOAT reading: %s", p.F32Error)
		}
		got, err := ReadPixels[float32](img)
		if err != nil {
			t.Fatal(err)
		}
		if len(got) != p.F32.Len() {
			t.Fatalf("%d pixels, cfitsio %d", len(got), p.F32.Len())
		}
		var nn []testutil.NullFloat
		if g.Image.Bitpix > 0 && p.F32.NullCount() > 0 {
			nn = p.F64NN.Resolve(t, "cfitsio TDOUBLE without null checks", func() []testutil.NullFloat {
				v, err := ReadPixels[float64](img)
				if err != nil {
					return nil
				}
				return testutil.NullFloats(v, nil)
			})
		}
		for k, w := range p.F32.Resolve(t, "cfitsio TFLOAT", func() []testutil.NullFloat { return cfitsioF32View(t, img, p.F32Ovf) }) {
			if pv, quirk := cfitsioFloatQuirk(t, img, k, w); quirk {
				if !sameFloat(float64(got[k]), float64(float32(pv))) {
					t.Fatalf("pixel %d: %v, stored %v (cfitsio TFLOAT reports %v)", k, got[k], float32(pv), w)
				}
				continue
			}
			if !w.Valid {
				if g.Image.Bitpix < 0 {
					if !math.IsNaN(float64(got[k])) {
						t.Fatalf("pixel %d: %v, cfitsio null (NaN)", k, got[k])
					}
					continue
				}
				if !sameFloat(float64(got[k]), float64(float32(nn[k].V))) {
					t.Fatalf("BLANK pixel %d: %v, cfitsio's scaled stored value %v", k, got[k], float32(nn[k].V))
				}
				continue
			}
			// cfitsio clamps values beyond float32 range to +-FLT_MAX
			// (status 412); Go's float32 conversion gives +-Inf.
			if p.F32Ovf && math.Abs(w.V) == math.MaxFloat32 && math.IsInf(float64(got[k]), int(math.Copysign(1, w.V))) {
				continue
			}
			if !sameFloat(float64(got[k]), w.V) {
				t.Fatalf("pixel %d: %v, cfitsio TFLOAT %v", k, got[k], w.V)
			}
		}
	})
}

// TestReadPixelsIntegerTargets: for every integer Go type T, CanConvert[T]
// and ReadPixels[T] accept the image exactly when every possible physical
// value (raw range x BSCALE + BZERO, cfitsio's BSCALE/BZERO) is an integer
// inside T's range; accepted reads equal cfitsio's exact TLONGLONG /
// TULONGLONG reading, refused reads return ErrTypeMismatch (never a
// silently wrapped value). When cfitsio's own equivalent type is an
// integer type, that Go type must be accepted.
func TestReadPixelsIntegerTargets(t *testing.T) {
	forEachImage(t, func(t *testing.T, s imageSample, i int, img *ImageHDU, g *testutil.CfrefHDU) {
		im := g.Image
		if im.Bitpix < 0 {
			for _, tg := range intTargets {
				if tg.can(img) {
					t.Errorf("CanConvert[%s] true for float BITPIX %d", tg.name, im.Bitpix)
				}
				if _, err := tg.read(img); !errors.Is(err, ErrTypeMismatch) {
					t.Errorf("ReadPixels[%s] on float BITPIX: err %v, want ErrTypeMismatch", tg.name, err)
				}
			}
			return
		}
		lo, hi, integral := physRange(im.Bitpix, float64(im.Bscale), float64(im.Bzero))
		// Exact expected values: cfitsio TULONGLONG for BITPIX 64 when it
		// succeeds without overflow, else TLONGLONG.
		var u64 []testutil.NullUint64
		if s.pix != nil && s.pix.HDUs[i].U64Error == "" && s.pix.HDUs[i].U64.Present() {
			u64 = s.pix.HDUs[i].U64.Resolve(t, "cfitsio TULONGLONG", func() []testutil.NullUint64 {
				v, ok, err := ReadPixelsMasked[uint64](img)
				if err != nil {
					return nil
				}
				return testutil.NullUint64s(v, ok)
			})
		}
		var want []*big.Int
		null := make([]bool, im.Npix)
		switch {
		case im.I64.Len() == int(im.Npix):
			for k, v := range im.I64.Resolve(t, "cfitsio TLONGLONG", func() []testutil.NullInt64 { return int64View(img) }) {
				null[k] = !v.Valid
				want = append(want, big.NewInt(v.V))
			}
		case len(u64) == int(im.Npix):
			for k, v := range u64 {
				null[k] = !v.Valid
				want = append(want, new(big.Int).SetUint64(v.V))
			}
		case integral && im.Phys.Len() == int(im.Npix):
			// cfref dumps TLONGLONG only for integer equivalent types; for
			// e.g. BITPIX 32 with BZERO 100 (cfitsio: DOUBLE_IMG) the
			// TDOUBLE values are exact integers below 2^53.
			for k, v := range im.Phys.Resolve(t, "cfitsio TDOUBLE", func() []testutil.NullFloat { return cfitsioPhysView(t, img) }) {
				null[k] = !v.Valid
				if v.Valid && (v.V != math.Trunc(v.V) || math.Abs(v.V) >= 1<<53) {
					want = nil
					break
				}
				bi, _ := new(big.Float).SetFloat64(v.V).Int(nil)
				want = append(want, bi)
			}
		}
		allowedAny := false
		for _, tg := range intTargets {
			allowed := integral && lo.Cmp(tg.lo) >= 0 && hi.Cmp(tg.hi) <= 0
			// cfitsio 4.7.0 does not range-check 64-bit equivalent types:
			// it reports LONGLONG_IMG for BITPIX 64 with any integral
			// BZERO (even the uint64 convention, which it then reads as
			// TULONGLONG), so for BITPIX 64 the equivalent type is not a
			// range statement.
			if gt, ok := equivGoType[im.Equivtype]; ok && gt == tg.name && !allowed && im.Bitpix != 64 {
				t.Errorf("cfitsio equivalent type is %s but the exact physical range [%v, %v] does not fit it", gt, lo, hi)
			}
			if got := tg.can(img); got != allowed {
				t.Errorf("CanConvert[%s] = %v, want %v (physical range [%.0f, %.0f], integral %v)", tg.name, got, allowed, lo, hi, integral)
			}
			got, err := tg.read(img)
			if !allowed {
				if !errors.Is(err, ErrTypeMismatch) {
					t.Errorf("ReadPixels[%s]: err %v, want ErrTypeMismatch", tg.name, err)
				}
				continue
			}
			allowedAny = true
			if err != nil {
				t.Errorf("ReadPixels[%s]: %v", tg.name, err)
				continue
			}
			if want == nil {
				t.Fatalf("no exact cfitsio integer reading in the golden (errors %v)", g.Errors)
			}
			for k := range got {
				if !null[k] && got[k].Cmp(want[k]) != 0 {
					t.Errorf("ReadPixels[%s] pixel %d = %v, cfitsio %v", tg.name, k, got[k], want[k])
					break
				}
			}
		}
		if integral && !allowedAny && !(im.Bitpix == 64 && im.Bzero != 0) {
			t.Errorf("no integer type accepted for integral scaling")
		}
	})
}

// TestReadPixelsMaskedMatchesCfitsio: the validity mask equals cfitsio's
// null flags (BLANK for integer BITPIX, NaN for float BITPIX) and valid
// pixels equal cfitsio's scaled TDOUBLE values.
func TestReadPixelsMaskedMatchesCfitsio(t *testing.T) {
	sawBlank := map[string]bool{}
	forEachImage(t, func(t *testing.T, s imageSample, i int, img *ImageHDU, g *testutil.CfrefHDU) {
		got, mask, err := ReadPixelsMasked[float64](img)
		if err != nil {
			t.Fatal(err)
		}
		if len(mask) != g.Image.Phys.Len() {
			t.Fatalf("%d mask entries, cfitsio %d pixels", len(mask), g.Image.Phys.Len())
		}
		for k, w := range g.Image.Phys.Resolve(t, "cfitsio TDOUBLE", func() []testutil.NullFloat { return cfitsioPhysView(t, img) }) {
			if p, ok := cfitsioFloatQuirk(t, img, k, w); ok {
				// stored +-Inf / -0 / subnormal: a valid value
				if !mask[k] || !sameFloat(got[k], p) {
					t.Fatalf("pixel %d: mask %v value %v, stored %v", k, mask[k], got[k], p)
				}
				continue
			}
			if mask[k] != w.Valid {
				t.Fatalf("pixel %d: mask %v, cfitsio null flag %v", k, mask[k], !w.Valid)
			}
			if !w.Valid {
				sawBlank[s.name] = true
				continue
			}
			if !sameFloat(got[k], w.V) {
				t.Fatalf("pixel %d: %v, cfitsio %v", k, got[k], w.V)
			}
		}
	})
	// The corpus must actually contain nulls of both kinds.
	for _, n := range []string{"cfwrite/i16_blank.fits", "cfwrite/u8_blank.fits", "cfwrite/i32_blank_bzero.fits", "astropy_data/blank.fits", "testprog.std", "cfwrite/f32_plain.fits"} {
		if !sawBlank[n] {
			t.Errorf("%s: no null pixels seen; the masked-read check did not exercise BLANK/NaN", n)
		}
	}
}

// TestReadSubsetMatchesCfitsio: ReadSubset equals cfitsio's
// fits_read_subset for >2-D cubes, BZERO unsigned data, BITPIX 64,
// non-integer BSCALE, strided reads and float32 targets
// (testdata/images/golden/subsets.json).
func TestReadSubsetMatchesCfitsio(t *testing.T) {
	var doc struct {
		Cfitsio string `json:"cfitsio"`
		Subsets []struct {
			File   string            `json:"file"`
			HDU    int               `json:"hdu"`
			Type   string            `json:"type"`
			Fpixel []int64           `json:"fpixel"`
			Lpixel []int64           `json:"lpixel"`
			Inc    []int64           `json:"inc"`
			Values []json.RawMessage `json:"values"`
		} `json:"subsets"`
	}
	loadJSON(t, imagesPath(t, "golden", "subsets.json"), &doc)
	if doc.Cfitsio != "4.7.0" || len(doc.Subsets) == 0 {
		t.Fatalf("subsets golden: cfitsio %q, %d subsets", doc.Cfitsio, len(doc.Subsets))
	}
	for _, sc := range doc.Subsets {
		name := fmt.Sprintf("%s[%d] %s %v-%v/%v", sc.File, sc.HDU, sc.Type, sc.Fpixel, sc.Lpixel, sc.Inc)
		t.Run(name, func(t *testing.T) {
			f, err := Open(imagesPath(t, strings.Split(sc.File, "/")...))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			h, err := f.HDU(sc.HDU)
			if err != nil {
				t.Fatal(err)
			}
			img := h.(*ImageHDU)
			lower := make([]int64, len(sc.Fpixel))
			upper := make([]int64, len(sc.Fpixel))
			for i := range lower {
				lower[i] = sc.Fpixel[i] - 1 // 1-based inclusive -> 0-based
				upper[i] = sc.Lpixel[i]     // inclusive -> exclusive
			}
			var got []string
			var gotF []float64
			switch sc.Type {
			case "i64":
				v, err := ReadSubset[int64](img, lower, upper, sc.Inc)
				if err != nil {
					t.Fatal(err)
				}
				for _, x := range v {
					got = append(got, fmt.Sprint(x))
				}
			case "u64":
				v, err := ReadSubset[uint64](img, lower, upper, sc.Inc)
				if err != nil {
					t.Fatal(err)
				}
				for _, x := range v {
					got = append(got, fmt.Sprint(x))
				}
			case "f32":
				v, err := ReadSubset[float32](img, lower, upper, sc.Inc)
				if err != nil {
					t.Fatal(err)
				}
				for _, x := range v {
					gotF = append(gotF, float64(x))
					got = append(got, fmt.Sprint(x))
				}
			case "f64":
				v, err := ReadSubset[float64](img, lower, upper, sc.Inc)
				if err != nil {
					t.Fatal(err)
				}
				gotF = v
				for _, x := range v {
					got = append(got, fmt.Sprint(x))
				}
			}
			if len(got) != len(sc.Values) {
				t.Fatalf("%d values, cfitsio %d", len(got), len(sc.Values))
			}
			for k, w := range sc.Values {
				if sc.Type == "f64" || sc.Type == "f32" {
					var wf testutil.NullFloat
					if err := wf.UnmarshalJSON(w); err != nil || !wf.Valid {
						t.Fatalf("bad golden value %s", w)
					}
					if sc.Type == "f32" {
						// printed with %.9g, which round-trips float32 exactly
						wf.V = float64(float32(wf.V))
					}
					if !sameFloat(gotF[k], wf.V) {
						t.Fatalf("element %d = %v, cfitsio %v", k, gotF[k], wf.V)
					}
					continue
				}
				if got[k] != string(w) {
					t.Fatalf("element %d = %s, cfitsio %s", k, got[k], w)
				}
			}
		})
	}
}
