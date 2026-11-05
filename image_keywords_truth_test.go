package fits

// Image scaling keywords (BSCALE, BZERO, BLANK) by the scaling-card rule
// (scalingCards), and the structural-keyword and write-side refusals.
// Truth: testdata/images/keywords (byte-assembled headers, gen.py
// KEYWORD_CASES) with cfitsio's reading of each (cfref dump: TDOUBLE with
// null flags; imgref pix: TFLOAT, and TDOUBLE with null checks off),
// testdata/images/duplicate (cfref dump), and
// testdata/images/golden/write_refused.json (cfitsio fits_set_bscale
// status). See testdata/images/PROVENANCE.md.

import (
	"bytes"
	"encoding/hex"
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// keywordSamples lists testdata/images/keywords; the cases every test run
// must see are named so deleting one fails.
func keywordSamples(t *testing.T) []string {
	t.Helper()
	entries, err := os.ReadDir(imagesPath(t, "keywords"))
	if err != nil {
		t.Fatalf("committed samples missing (failure, not a skip): %v", err)
	}
	var names []string
	for _, e := range entries {
		names = append(names, e.Name())
	}
	sort.Strings(names)
	have := map[string]bool{}
	for _, n := range names {
		have[n] = true
	}
	for _, n := range []string{
		"dexp_upper.fits", "dexp_lower.fits", "dexp_hex.fits", "dexp_hex_two_d.fits", "bscale_eq_col9_after_eq.fits",
		"dup_bscale.fits", "dup_bzero.fits", "dup_blank.fits", "dup_bscale_bad_last.fits", "dup_bscale_bad_first.fits", "dup_blank_bad_last.fits",
		"bscale_string.fits", "bzero_string.fits", "bscale_inf.fits", "bscale_nan.fits", "bscale_underflow.fits", "bscale_subnormal.fits",
		"bscale_subnormal_exact.fits", "bscale_overflow.fits", "bscale_noeq_space.fits", "bscale_key_eq.fits", "bscale_hex.fits",
		"bscale_hex_noexp.fits", "bscale_empty.fits", "bscale_slash.fits", "bscale_no_equals.fits", "bscale_comment_equals.fits",
		"bscale_lower.fits", "bscale_hierarch.fits", "bscale_zero.fits", "bscale_junk.fits", "bscale_trailing_e.fits",
		"blank_plus.fits", "blank_float.fits", "blank_empty.fits", "blank_overflow.fits", "blank_string.fits", "blank_hex.fits",
		"blank_wrap16.fits", "blank_wrap32.fits", "blank_magic.fits", "bscale_min_normal_rounded.fits", "bscale_below_min_normal.fits",
	} {
		if !have[n] {
			t.Errorf("testdata/images/keywords/%s missing", n)
		}
	}
	if len(names) < 51 {
		t.Fatalf("testdata/images/keywords has %d samples, want >= 51", len(names))
	}
	return names
}

// checkKeywordReads compares every read of img with cfitsio's reading of
// the same file. The samples store 0 and 1 in their first two pixels, so
// the BZERO and BSCALE cfitsio applied are its values there (read with
// null checks off).
func checkKeywordReads(t *testing.T, img *ImageHDU, g *testutil.CfrefHDU, p pixHDU) {
	t.Helper()
	phys := g.Image.Phys.Values(t)
	f32 := p.F32.Values(t)
	var nn []testutil.NullFloat
	if p.F64NN.Present() {
		nn = p.F64NN.Values(t)
	}
	if g.Image.Phys.NullCount() > 0 && nn == nil {
		t.Fatalf("golden has null pixels but no TDOUBLE reading without null checks")
	}
	base := phys
	if nn != nil {
		base = nn
	}
	wantZero, wantScale := base[0].V, base[1].V-base[0].V
	if !sameFloat(bzeroOf(t, img), wantZero) || !sameFloat(bscaleOf(t, img), wantScale) {
		t.Errorf("BSCALE/BZERO = %v/%v, cfitsio applies %v/%v", bscaleOf(t, img), bzeroOf(t, img), wantScale, wantZero)
	}
	got, err := ReadPixels[float64](img)
	if err != nil {
		t.Fatal(err)
	}
	m, mask, err := ReadPixelsMasked[float64](img)
	if err != nil {
		t.Fatal(err)
	}
	g32, err := ReadPixels[float32](img)
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != len(phys) || len(f32) != len(phys) {
		t.Fatalf("%d pixels, cfitsio %d / %d", len(got), len(phys), len(f32))
	}
	for k, w := range phys {
		if mask[k] != w.Valid || f32[k].Valid != w.Valid {
			t.Errorf("pixel %d: mask %v, cfitsio null flag TDOUBLE %v TFLOAT %v", k, mask[k], !w.Valid, !f32[k].Valid)
			continue
		}
		want, want32 := w.V, f32[k].V
		if !w.Valid {
			want, want32 = nn[k].V, float64(float32(nn[k].V))
		} else if !sameFloat(m[k], w.V) {
			t.Errorf("ReadPixelsMasked pixel %d = %v, cfitsio %v", k, m[k], w.V)
		}
		if !sameFloat(got[k], want) {
			t.Errorf("ReadPixels[float64] pixel %d = %v, cfitsio %v", k, got[k], want)
		}
		if !sameFloat(float64(g32[k]), want32) {
			t.Errorf("ReadPixels[float32] pixel %d = %v, cfitsio TFLOAT %v", k, g32[k], want32)
		}
	}
}

// keywordRule is one entry of testdata/images/golden/keywords_rule.json
// (gen.py gen_keyword_rule): the scaling-card rule's outcome for a
// keywords/ sample and astropy 8.0.2's reading of it.
type keywordRule struct {
	Rule      string `json:"rule"`
	RefsAgree bool   `json:"refs_agree"`
	Astropy   struct {
		Dtype  string    `json:"dtype"`
		Values []*string `json:"values"` // Python repr of each value, null = NaN
		Error  string    `json:"error"`
	} `json:"astropy"`
}

func loadKeywordRules(t *testing.T) map[string]keywordRule {
	t.Helper()
	var m map[string]keywordRule
	loadJSON(t, imagesPath(t, "golden", "keywords_rule.json"), &m)
	return m
}

// checkKeywordsAstropy compares the masked float64 read with astropy's
// scaled getdata of the sample, in astropy's dtype (it scales BITPIX 16
// into float32), BLANK pixels as NaN.
func checkKeywordsAstropy(t *testing.T, img *ImageHDU, r keywordRule) {
	t.Helper()
	if r.Astropy.Error != "" {
		t.Fatalf("astropy refuses an accepted sample (%s)", r.Astropy.Error)
	}
	v, ok, err := ReadPixelsMasked[float64](img)
	if err != nil {
		t.Fatal(err)
	}
	if len(v) != len(r.Astropy.Values) {
		t.Fatalf("%d pixels, astropy %d", len(v), len(r.Astropy.Values))
	}
	for k, a := range r.Astropy.Values {
		got := v[k]
		if r.Astropy.Dtype == "float32" {
			got = float64(float32(got))
		}
		if a == nil {
			if ok[k] {
				t.Errorf("pixel %d = %v, astropy NaN (null)", k, got)
			}
			continue
		}
		want, err := strconv.ParseFloat(*a, 64)
		if err != nil {
			t.Fatal(err)
		}
		if !ok[k] || !sameFloat(got, want) {
			t.Errorf("pixel %d = %v (valid %v), astropy %v", k, got, ok[k], want)
		}
	}
}

// checkKeywordRefused: every read of img that applies the scaling cards
// fails with ErrInvalidScaling, and so does OverwritePixels in ModeEdit
// (the file is left as it was).
func checkKeywordRefused(t *testing.T, img *ImageHDU) {
	t.Helper()
	is := func(what string, err error) {
		t.Helper()
		if !errors.Is(err, ErrInvalidScaling) {
			t.Errorf("%s: %v, want ErrInvalidScaling", what, err)
		}
	}
	_, err := img.BSCALE()
	is("BSCALE", err)
	_, err = img.BZERO()
	is("BZERO", err)
	_, _, err = img.BLANK()
	is("BLANK", err)
	_, err = ReadPixels[float64](img)
	is("ReadPixels[float64]", err)
	_, err = ReadPixels[float32](img)
	is("ReadPixels[float32]", err)
	_, err = ReadPixels[int64](img)
	is("ReadPixels[int64]", err)
	_, _, err = ReadPixelsMasked[float64](img)
	is("ReadPixelsMasked", err)
	_, err = ReadFloat32(img)
	is("ReadFloat32", err)
	if CanConvert[float64](img) {
		t.Error("CanConvert[float64] = true")
	}
	if img.rec.file.mode == ModeEdit {
		n := img.pixelCount()
		var err error
		if img.BITPIX() == 16 {
			err = OverwritePixels(img, make([]int16, n))
		} else {
			err = OverwritePixels(img, make([]int32, n))
		}
		is("OverwritePixels", err)
	}
}

// TestScalingKeywordsRule: the scaling-card rule (scalingCards) on every
// keywords/ sample (D and d exponents, duplicate BSCALE/BZERO/BLANK incl.
// a bad last or first card, quoted, logical, complex, Inf/NaN,
// overflowing, underflowing, subnormal and hex values, no '=' in column 9
// or none at all, an empty value, lowercase and HIERARCH names, BLANK with
// '+', a fraction, hex, overflow, outside the stored type, and cfitsio's
// NULL_UNDEFINED 1234554321), the same three ways: from the raw cards,
// after the full header is parsed, and after an unrelated in-memory edit
// in ModeEdit (which is then flushed, and cfitsio must still read the file
// as before). The expected outcome is keywords_rule.json (gen.py): an
// accepted sample reads exactly as cfitsio reads it (BSCALE/BZERO/BLANK,
// ReadPixels float64 and float32, ReadPixelsMasked: cfref dump and imgref
// pix) and as astropy reads it (its scaled getdata); for every other
// sample, where the references disagree, refuse the file, or drop the
// card, every read and OverwritePixels fail with ErrInvalidScaling. The
// generator asserts that no accepted sample is read differently by the two
// references. Before this rule the library read every sample as cfitsio
// does (a bad card reset to the default, HIERARCH BSCALE applied), so 39
// samples read differently from astropy, or after Header() than before
// on the ReadFloat32 path.
func TestScalingKeywordsRule(t *testing.T) {
	dir := t.TempDir()
	rules := loadKeywordRules(t)
	parsedRuns, accepted, refused := 0, 0, 0
	for _, name := range keywordSamples(t) {
		r, ok := rules[name]
		if !ok {
			t.Fatalf("keywords_rule.json has no entry for %s", name)
		}
		if r.Rule == "accept" {
			accepted++
		} else {
			refused++
		}
		t.Run(name, func(t *testing.T) {
			sample := "keywords/" + name
			path := imagesPath(t, "keywords", name)
			g := loadImagesCfref(t, sample)
			pix := loadPixGolden(t, sample)
			if len(g.HDUs) != 1 || g.HDUs[0].Image == nil || len(g.HDUs[0].Errors) != 0 {
				t.Fatalf("golden: want one image cfitsio read without errors, got %+v", g.HDUs)
			}
			gh, ph := &g.HDUs[0], pix.HDUs[0]
			check := func(t *testing.T, img *ImageHDU) {
				if r.Rule != "accept" {
					checkKeywordRefused(t, img)
					return
				}
				checkKeywordReads(t, img, gh, ph)
				checkKeywordsAstropy(t, img, r)
			}

			f, err := Open(path)
			if err != nil {
				t.Fatalf("Open: %v (cfitsio opens it)", err)
			}
			defer f.Close()
			img, err := f.Primary()
			if err != nil {
				t.Fatal(err)
			}
			t.Run("raw", func(t *testing.T) {
				if img.rec.parsed != nil {
					t.Fatal("header already parsed")
				}
				check(t, img)
			})
			_, herr := HeaderOf(img)
			if herr == nil {
				parsedRuns++
				t.Run("parsed", func(t *testing.T) {
					if img.rec.parsed == nil {
						t.Fatal("HeaderOf succeeded but the parsed header is not cached")
					}
					check(t, img)
				})
			} else {
				t.Logf("header package cannot parse the header (%v); the raw-card reads above are the reads", herr)
			}

			cp := filepath.Join(dir, name)
			copyFile(t, path, cp)
			orig := mustRead(t, cp)
			e, err := OpenForEdit(cp)
			if err != nil {
				t.Fatal(err)
			}
			eimg, err := e.Primary()
			if err != nil {
				t.Fatal(err)
			}
			h, err := HeaderOf(eimg)
			if err != nil {
				e.Close()
				if herr == nil {
					t.Fatalf("HeaderOf in ModeEdit: %v (it parsed in ModeRead)", err)
				}
				return
			}
			if err := h.Set("OBJECT", "edited", ""); err != nil {
				t.Fatal(err)
			}
			t.Run("edited", func(t *testing.T) { check(t, eimg) })
			if r.Rule != "accept" && !bytes.Equal(mustRead(t, cp), orig) {
				t.Fatal("a refused OverwritePixels changed the file")
			}
			if err := e.Flush(); err != nil {
				// A HIERARCH card named like a standard keyword
				// (bscale_hierarch) is refused by Flush too (edit.go
				// hierarchAlias), which leaves the file as it was.
				if !strings.HasPrefix(name, "bscale_hierarch") || !errors.Is(err, header.ErrInvalidKeyword) {
					t.Fatalf("Flush: %v", err)
				}
				e.Close()
				if !bytes.Equal(mustRead(t, cp), orig) {
					t.Errorf("Flush refused (%v) but changed the file", err)
				}
				return
			}
			if err := e.Close(); err != nil {
				t.Fatal(err)
			}
			live := testutil.CfrefDump(t, cp)
			if d := live.HDUs[0].Image.Phys.Same(gh.Image.Phys); d != "" {
				t.Errorf("cfitsio reads the flushed file differently: %s", d)
			}
		})
	}
	if parsedRuns < 30 || accepted < 10 || refused < 30 {
		t.Errorf("%d samples exercised the parsed-header path, %d accepted, %d refused", parsedRuns, accepted, refused)
	}
}

// TestDuplicateStructuralKeywordsRefused: a header with a second BITPIX,
// NAXIS, NAXISn (n <= NAXIS), SIMPLE, XTENSION, PCOUNT, GCOUNT or TFIELDS
// is refused by every open path with *DuplicateKeywordError naming the HDU
// and keyword. cfitsio opens each of these files (golden), reading
// BITPIX/NAXIS/NAXISn by position (the first) and an image's PCOUNT/GCOUNT
// in its keyword loop (the last), so a stray second BITPIX would change the
// element size (dup_bitpix: cfitsio 16, the old last-wins scan 8, data
// misread) and a second NAXIS1 the shape: no single reading exists.
func TestDuplicateStructuralKeywordsRefused(t *testing.T) {
	cases := map[string]DuplicateKeywordError{
		"dup_bitpix.fits":         {HDU: 0, Keyword: "BITPIX"},
		"dup_naxis.fits":          {HDU: 0, Keyword: "NAXIS"},
		"dup_naxis1_smaller.fits": {HDU: 0, Keyword: "NAXIS1"},
		"dup_naxis1_larger.fits":  {HDU: 0, Keyword: "NAXIS1"},
		// a duplicated stray NAXIS3 (n > NAXIS, ignored) after it must not
		// hide the real NAXIS1 duplicate
		"dup_naxis1_and_stray_naxis3.fits": {HDU: 0, Keyword: "NAXIS1"},
		"dup_simple.fits":                  {HDU: 0, Keyword: "SIMPLE"},
		"dup_pcount_ext.fits":              {HDU: 1, Keyword: "PCOUNT"},
		"dup_gcount_ext.fits":              {HDU: 1, Keyword: "GCOUNT"},
		"dup_xtension_ext.fits":            {HDU: 1, Keyword: "XTENSION"},
		"dup_tfields_ext.fits":             {HDU: 1, Keyword: "TFIELDS"},
	}
	entries, err := os.ReadDir(imagesPath(t, "duplicate"))
	if err != nil {
		t.Fatalf("committed samples missing (failure, not a skip): %v", err)
	}
	if len(entries) != len(cases) {
		t.Fatalf("testdata/images/duplicate has %d files, the test names %d", len(entries), len(cases))
	}
	for name, want := range cases {
		t.Run(name, func(t *testing.T) {
			g := loadImagesCfref(t, "duplicate/"+name)
			if g.NHDU < want.HDU+1 {
				t.Fatalf("cfitsio golden has %d HDUs; the documented difference is that cfitsio opens HDU %d", g.NHDU, want.HDU)
			}
			for _, op := range openPaths {
				_, err := op.open(t, imagesPath(t, "duplicate", name))
				var e *DuplicateKeywordError
				if !errors.As(err, &e) || *e != want || !errors.Is(err, ErrDuplicateKeyword) {
					t.Errorf("%s: err %v, want %+v", op.name, err, want)
				}
			}
		})
	}
}

// writeRefusedGolden is golden/write_refused.json: cfitsio's status when
// fits_set_bscale is given each scaling.
type writeRefusedGolden map[string]struct {
	Imgtype  int    `json:"imgtype"`
	Datatype int    `json:"datatype"`
	Bscale   string `json:"bscale"`
	Bzero    string `json:"bzero"`
	Status   int    `json:"status"`
	Refused  string `json:"refused"`
}

// writeScaled writes 1,2,3,4 as T with hdr into a fresh file and returns
// WriteImage's error and the number of HDUs the file then holds.
func writeScaled[T Numeric](t *testing.T, path string, hdr *header.Header) (error, int) {
	t.Helper()
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	_, werr := WriteImage(f, hdr, []int64{4}, []T{1, 2, 3, 4})
	n := f.NumHDU()
	if err := f.Close(); err != nil && werr == nil {
		t.Fatal(err)
	}
	return werr, n
}

// writeRefusedKeywordSamples are the keywords/ samples whose header
// WriteImage refuses as a caller header: every sample the scaling-card
// rule refuses (keywords_rule.json), a keyword given twice (the readers
// apply the last; a caller who wrote two meant one of them) and BSCALE 0.
func writeRefusedKeywordSamples(t *testing.T) map[string]bool {
	refuse := map[string]bool{"bscale_zero.fits": true}
	for n, r := range loadKeywordRules(t) {
		if r.Rule != "accept" || strings.HasPrefix(n, "dup_") {
			refuse[n] = true
		}
	}
	return refuse
}

// TestWriteImageRefusesInvalidScaling: WriteImage returns ErrInvalidScaling
// and writes no HDU for (a) BSCALE = 0, as a real or an integer card, for
// every BITPIX cfitsio refuses it on (golden: fits_set_bscale status 322
// ZERO_SCALE); every value would otherwise be stored as BZERO's inverse
// and read back as BZERO; (b) caller headers cfitsio would read
// differently from the values WriteImage stores with: a second BSCALE,
// BZERO or BLANK card (cfitsio applies the last: keywords/dup_*.fits),
// a quoted or otherwise non-numeric BSCALE/BZERO and a non-integer BLANK
// (cfitsio ignores them: keywords/bscale_string.fits, blank_float.fits),
// and non-finite or out-of-range BSCALE/BZERO (cfitsio's ffc2dd rejects
// them on read; the header package itself refuses to build NaN/Inf cards,
// so those reach WriteImage only as card bytes taken from a file). A plain
// BSCALE 2 still writes.
func TestWriteImageRefusesInvalidScaling(t *testing.T) {
	dir := t.TempDir()
	var g writeRefusedGolden
	loadJSON(t, imagesPath(t, "golden", "write_refused.json"), &g)
	if len(g) < 5 {
		t.Fatalf("write_refused.json has %d cases", len(g))
	}
	try := func(t *testing.T, datatype int, hdr *header.Header) (error, int) {
		p := filepath.Join(dir, strings.ReplaceAll(t.Name(), "/", "_")+".fits")
		switch datatype {
		case 11:
			return writeScaled[uint8](t, p, hdr)
		case 21:
			return writeScaled[int16](t, p, hdr)
		case 81:
			return writeScaled[int64](t, p, hdr)
		case 42:
			return writeScaled[float32](t, p, hdr)
		}
		t.Fatalf("datatype %d not mapped", datatype)
		return nil, 0
	}
	refused := func(t *testing.T, err error, n int) {
		t.Helper()
		want := ErrInvalidScaling
		if strings.HasPrefix(t.Name(), "TestWriteImageRefusesInvalidScaling/header_of_keywords/bscale_lower.fits") {
			want = header.ErrInvalidKeyword
		}
		if !errors.Is(err, want) || n != 0 {
			t.Errorf("WriteImage: err %v, %d HDUs written; want %v and none", err, n, want)
		}
	}
	for name, c := range g {
		t.Run(name, func(t *testing.T) {
			if c.Status != 322 || c.Refused != "set_bscale" || c.Bscale != "0" {
				t.Fatalf("golden %+v: want cfitsio fits_set_bscale ZERO_SCALE (322) for BSCALE 0", c)
			}
			for _, card := range []any{0.0, int64(0)} {
				h := header.New()
				if err := h.Set("BSCALE", card, ""); err != nil {
					t.Fatal(err)
				}
				if c.Bzero != "0" {
					if err := h.Set("BZERO", 5.0, ""); err != nil {
						t.Fatal(err)
					}
				}
				err, n := try(t, c.Datatype, h)
				refused(t, err, n)
			}
		})
	}
	// Caller headers taken from the keywords/ samples (HeaderOf keeps each
	// card's bytes, so WriteImage emits them as cfitsio read them there).
	// WriteImage must refuse exactly the headers whose scaling cards
	// cfitsio would read differently from the caller's card (a second
	// BSCALE/BZERO/BLANK: cfitsio applies the last; a value cfitsio cannot
	// convert and replaces by the default: D exponent in lowercase, quoted,
	// logical, complex, Inf/NaN, over/underflow, subnormal, trailing junk;
	// BSCALE 0, also from an empty value; BLANK = 1234554321, which cfitsio
	// takes for "no BLANK"). Every other header writes, and cfitsio's live
	// reading of the written file equals the library's.
	refuse := writeRefusedKeywordSamples(t)
	// BSCALE 2^-1022 is a valid scaling card, but every nonzero value
	// overflows under it (cfitsio: NUM_OVERFLOW,
	// golden/keywords_overwrite.json; TestWriteImageKeywordsMatchCfitsio).
	overflow := map[string]bool{"bscale_min_normal_rounded.fits": true}
	written, refusedN, overflowN := 0, 0, 0
	for _, name := range keywordSamples(t) {
		t.Run("header of keywords/"+name, func(t *testing.T) {
			src := openImage(t, imagesPath(t, "keywords", name))
			hdr, err := HeaderOf(src)
			if err != nil {
				if refuse[name] {
					t.Fatalf("expected a refusal but the header does not parse: %v", err)
				}
				t.Logf("header does not parse (%v): nothing to write with", err)
				return
			}
			p := filepath.Join(dir, "w_"+name)
			err, n := writeScaled[int16](t, p, hdr.Clone())
			if refuse[name] {
				refusedN++
				refused(t, err, n)
				return
			}
			if overflow[name] {
				overflowN++
				if !errors.Is(err, ErrNumOverflow) || n != 0 {
					t.Errorf("WriteImage: err %v, %d HDUs written; want ErrNumOverflow and none", err, n)
				}
				return
			}
			if err != nil || n != 1 {
				t.Fatalf("WriteImage: err %v, %d HDUs", err, n)
			}
			written++
			img := openImage(t, p)
			v, ok, err := ReadPixelsMasked[float64](img)
			if err != nil {
				t.Fatal(err)
			}
			live := testutil.CfrefDump(t, p)
			if d := live.HDUs[0].Image.Phys.Diff(testutil.NullFloats(v, ok)); d != "" {
				t.Errorf("cfitsio reads the written file differently from the library: %s", d)
			}
		})
	}
	if refusedN != len(refuse) || overflowN != len(overflow) || written < 6 {
		t.Errorf("%d of %d expected refusals and %d of %d overflows ran, %d headers written",
			refusedN, len(refuse), overflowN, len(overflow), written)
	}
	// Duplicates built in memory (Header.Add allows them).
	for _, c := range []struct {
		name string
		add  func(h *header.Header) error
	}{
		{"Add BSCALE twice", func(h *header.Header) error {
			if err := h.Set("BSCALE", 2.0, ""); err != nil {
				return err
			}
			return h.Add("BSCALE", 3.0, "")
		}},
		{"Add BZERO twice", func(h *header.Header) error {
			if err := h.Set("BZERO", 10.0, ""); err != nil {
				return err
			}
			return h.Add("BZERO", 20.0, "")
		}},
		{"Add BLANK twice", func(h *header.Header) error {
			if err := h.Set("BLANK", int64(1), ""); err != nil {
				return err
			}
			return h.Add("BLANK", int64(2), "")
		}},
		{"string BSCALE", func(h *header.Header) error { return h.Set("BSCALE", "2.0", "") }},
		{"real BLANK", func(h *header.Header) error { return h.Set("BLANK", 2.0, "") }},
	} {
		t.Run(c.name, func(t *testing.T) {
			h := header.New()
			if err := c.add(h); err != nil {
				t.Fatal(err)
			}
			err, n := try(t, 21, h)
			refused(t, err, n)
		})
	}
	// A scaling card that cannot be encoded never reaches the disk (a
	// value edited to NaN in the Cards() slice of an in-memory header):
	// it must still be refused as invalid scaling, not dropped (which
	// would store the data unscaled).
	for _, key := range []string{"BSCALE", "BZERO"} {
		t.Run(key+" edited to NaN in Cards()", func(t *testing.T) {
			h := header.New()
			if err := h.Set(key, 2.0, ""); err != nil {
				t.Fatal(err)
			}
			h.Cards()[0].Value = math.NaN()
			err, n := try(t, 21, h)
			refused(t, err, n)
		})
	}
	t.Run("BSCALE 2 writes", func(t *testing.T) {
		h := header.New()
		if err := h.Set("BSCALE", 2.0, ""); err != nil {
			t.Fatal(err)
		}
		err, n := try(t, 21, h)
		if err != nil || n != 1 {
			t.Fatalf("err %v, %d HDUs", err, n)
		}
	})
}

// TestOverwritePixelsRefusesZeroBSCALE: OverwritePixels on an HDU whose
// BSCALE is 0 (keywords/bscale_zero.fits, and a BSCALE 0 set in memory on
// a parsed header) returns ErrInvalidScaling and leaves the file's bytes
// unchanged, instead of storing values that read back as BZERO. (An empty
// BSCALE value, which cfitsio's strtod reads as 0, is refused by the
// scaling-card rule itself: TestScalingKeywordsRule.)
func TestOverwritePixelsRefusesZeroBSCALE(t *testing.T) {
	dir := t.TempDir()
	for _, c := range []struct {
		name, src string
		edit      bool
	}{
		{"bscale_zero", "keywords/bscale_zero.fits", false},
		{"in-memory BSCALE 0", "cfwrite/i16_plain.fits", true},
	} {
		t.Run(c.name, func(t *testing.T) {
			p := filepath.Join(dir, strings.ReplaceAll(c.name, " ", "_")+".fits")
			copyFile(t, imagesPath(t, filepath.FromSlash(c.src)), p)
			before, err := os.ReadFile(p)
			if err != nil {
				t.Fatal(err)
			}
			e, err := OpenForEdit(p)
			if err != nil {
				t.Fatal(err)
			}
			img, err := e.Primary()
			if err != nil {
				t.Fatal(err)
			}
			if c.edit {
				h, err := HeaderOf(img)
				if err != nil {
					t.Fatal(err)
				}
				if err := h.Set("BSCALE", 0.0, ""); err != nil {
					t.Fatal(err)
				}
			}
			if bscaleOf(t, img) != 0 {
				t.Fatalf("BSCALE() = %v, want 0", bscaleOf(t, img))
			}
			vals := make([]int16, img.pixelCount())
			for i := range vals {
				vals[i] = int16(i + 1)
			}
			if err := OverwritePixels(img, vals); !errors.Is(err, ErrInvalidScaling) {
				t.Errorf("OverwritePixels: err %v, want ErrInvalidScaling", err)
			}
			e.Close()
			after, err := os.ReadFile(p)
			if err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(before, after) {
				t.Error("file changed")
			}
		})
	}
}

// TestOverwritePixelsRefusesUnencodableScaling: a BSCALE/BZERO/BLANK card
// of a parsed header whose value was edited in the Cards() slice to
// something that cannot be written (NaN) never reaches the disk, so the
// data cannot be stored with it: OverwritePixels returns ErrInvalidScaling
// and leaves the file unchanged (instead of storing the values with the
// default scaling while the file keeps its old card, e.g. BZERO 1000 on
// cfwrite/i16_bzero1000.fits, which cfitsio would then read 1000 off).
func TestOverwritePixelsRefusesUnencodableScaling(t *testing.T) {
	dir := t.TempDir()
	for _, key := range []string{"BSCALE", "BZERO"} {
		t.Run(key, func(t *testing.T) {
			p := filepath.Join(dir, key+".fits")
			copyFile(t, imagesPath(t, "cfwrite", "i16_bzero1000.fits"), p)
			before, err := os.ReadFile(p)
			if err != nil {
				t.Fatal(err)
			}
			e, err := OpenForEdit(p)
			if err != nil {
				t.Fatal(err)
			}
			img, err := e.Primary()
			if err != nil {
				t.Fatal(err)
			}
			h, err := HeaderOf(img)
			if err != nil {
				t.Fatal(err)
			}
			found := false
			for i, c := range h.Cards() {
				if c.Key == key {
					h.Cards()[i].Value, h.Cards()[i].Raw = math.NaN(), [80]byte{}
					found = true
				}
			}
			if !found {
				t.Fatalf("no %s card in the sample", key)
			}
			vals := make([]int16, img.pixelCount())
			if err := OverwritePixels(img, vals); !errors.Is(err, ErrInvalidScaling) {
				t.Errorf("OverwritePixels: err %v, want ErrInvalidScaling", err)
			}
			e.Close()
			after, err := os.ReadFile(p)
			if err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(before, after) {
				t.Error("file changed")
			}
		})
	}
}

// TestUint64HeaderEditFlushCfitsioReads: the uint64 convention's integer
// card BZERO = 9223372036854775808 parses (header TypeBigInt), so a uint64
// image can be opened for editing, its header edited and checksummed, and
// flushed; cfitsio then reads every value exactly as before (TULONGLONG),
// verifies the checksums, and fitsverify reports nothing. Covers the files
// cfitsio and astropy write.
func TestUint64HeaderEditFlushCfitsioReads(t *testing.T) {
	dir := t.TempDir()
	for _, s := range []string{"cfwrite/u64_conv.fits", "astropy_written/u64.fits", "astropy_written/u64_cube.fits"} {
		t.Run(s, func(t *testing.T) {
			want := loadPixGolden(t, s).HDUs[0]
			if want.U64Error != "" || want.U64.Len() == 0 {
				t.Fatalf("golden has no TULONGLONG reading: %s", want.U64Error)
			}
			p := filepath.Join(dir, filepath.Base(s))
			copyFile(t, imagesPath(t, filepath.FromSlash(s)), p)
			e, err := OpenForEdit(p)
			if err != nil {
				t.Fatal(err)
			}
			img, err := e.Primary()
			if err != nil {
				t.Fatal(err)
			}
			h, err := HeaderOf(img)
			if err != nil {
				t.Fatalf("HeaderOf: %v", err)
			}
			if c, ok := h.Get("BZERO"); !ok || c.Type != header.TypeBigInt {
				t.Fatalf("BZERO card %+v, want TypeBigInt", c)
			}
			if err := h.Set("OBJECT", "edited", ""); err != nil {
				t.Fatal(err)
			}
			if _, _, err := WriteChecksum(img); err != nil {
				t.Fatalf("WriteChecksum: %v", err)
			}
			if err := e.Flush(); err != nil {
				t.Fatalf("Flush: %v", err)
			}
			if err := e.Close(); err != nil {
				t.Fatal(err)
			}
			if d := imgrefPix(t, p).HDUs[0].U64.Same(want.U64); d != "" {
				t.Errorf("cfitsio TULONGLONG after the edit: %s", d)
			}
			for i, c := range imgrefChksum(t, p).HDUs {
				if c.VerifyDatasum != 1 || c.VerifyHDUSum != 1 {
					t.Errorf("HDU %d: cfitsio fits_verify_chksum %d/%d, want 1/1", i, c.VerifyDatasum, c.VerifyHDUSum)
				}
			}
			if r := testutil.Fitsverify(t, p); !r.OK() {
				t.Errorf("fitsverify %d warnings, %d errors:\n%s", r.Warnings, r.Errors, r.Report)
			}
			got, err := ReadPixels[uint64](openImage(t, p))
			if err != nil {
				t.Fatal(err)
			}
			want.U64.Check(t, "ReadPixels[uint64] after the edit", testutil.NullUint64s(got, nil))
		})
	}
}

// TestCStrtodMatchesGlibc: the library's port of C strtod and strtoll
// (cStrtod, cStrtoll: what cfitsio's ffc2dd / ffc2jj call to convert
// BSCALE, BZERO and BLANK, and atof for DATASUM) returns the same value
// bits, bytes consumed and ERANGE flag as the C library cfitsio runs on
// (glibc), for 3084 inputs: boundary cases (hex floats with and without
// exponent, inf/nan forms, "nan(...)", over/underflow, subnormals exact
// and inexact, values that round up to the smallest normal, int64
// limits) and seeded random tokens and numbers (golden/strtod.txt from
// imgref strtod).
func TestCStrtodMatchesGlibc(t *testing.T) {
	p := imagesPath(t, "golden", "strtod.txt")
	b := testutil.MustReadGolden(t, p)
	lines := strings.Split(strings.TrimSuffix(string(b), "\n"), "\n")
	if len(lines) < 3000 {
		t.Fatalf("strtod.txt has %d lines", len(lines))
	}
	sawErange, sawSubnormal, sawHex := 0, 0, 0
	for _, ln := range lines {
		parts := strings.Split(ln, "\t")
		if len(parts) != 3 {
			t.Fatalf("bad line %q", ln)
		}
		in := parts[0]
		v, n, er := cStrtod(in)
		bits := math.Float64bits(v)
		if math.IsNaN(v) {
			bits = 0x7ff8000000000000
		}
		e := 0
		if er {
			e, sawErange = 1, sawErange+1
		}
		if v != 0 && math.Abs(v) < 0x1p-1022 {
			sawSubnormal++
		}
		if strings.Contains(strings.ToLower(in), "0x") {
			sawHex++
		}
		if got := fmt.Sprintf("%016x %d %d", bits, n, e); got != parts[1] {
			t.Errorf("strtod(%q) = %s, glibc %s", in, got, parts[1])
		}
		x, n2, er2 := cStrtoll(in)
		e2 := 0
		if er2 {
			e2 = 1
		}
		if got := fmt.Sprintf("%d %d %d", x, n2, e2); got != parts[2] {
			t.Errorf("strtoll(%q) = %s, glibc %s", in, got, parts[2])
		}
	}
	if sawErange < 20 || sawSubnormal < 20 || sawHex < 500 {
		t.Errorf("weak corpus: %d ERANGE, %d subnormal, %d hex inputs", sawErange, sawSubnormal, sawHex)
	}
}

// TestC2DDAndC2JJMatchCfitsio: cfC2DD and cfC2JJ, the library's ports of
// cfitsio's ffc2dd / ffc2jj (BSCALE/BZERO and BLANK conversion), accept
// exactly the value strings cfitsio's own functions accept, with the same
// value bits, on the strtod boundary inputs, the random tokens and the
// wrapper cases (trailing blank and blank-separated junk accepted, only
// the first 'D' replaced, a 'D' string over 72 characters refused,
// Inf/NaN/ERANGE refused): golden/c2.txt, imgref c2 calling ffc2dd and
// ffc2jj from libcfitsio 4.7.0.
func TestC2DDAndC2JJMatchCfitsio(t *testing.T) {
	p := imagesPath(t, "golden", "c2.txt")
	b := testutil.MustReadGolden(t, p)
	lines := strings.Split(strings.TrimSuffix(string(b), "\n"), "\n")
	if len(lines) < 1500 {
		t.Fatalf("c2.txt has %d lines", len(lines))
	}
	okDD, okJJ := 0, 0
	have := map[string]bool{}
	for _, ln := range lines {
		parts := strings.Split(ln, "\t")
		if len(parts) != 3 {
			t.Fatalf("bad line %q", ln)
		}
		in := parts[0]
		have[in] = true
		var ddSt, jjSt int
		var ddBits uint64
		var jj int64
		if _, err := fmt.Sscanf(parts[1], "%d %x", &ddSt, &ddBits); err != nil {
			t.Fatalf("line %q: %v", ln, err)
		}
		if _, err := fmt.Sscanf(parts[2], "%d %d", &jjSt, &jj); err != nil {
			t.Fatalf("line %q: %v", ln, err)
		}
		f, ok := cfC2DD(in)
		if ok != (ddSt == 0) || ok && math.Float64bits(f) != ddBits {
			t.Errorf("cfC2DD(%q) = %v (%016x), %v; cfitsio ffc2dd status %d bits %016x", in, f, math.Float64bits(f), ok, ddSt, ddBits)
		}
		x, ok2 := cfC2JJ(in)
		if ok2 != (jjSt == 0) || ok2 && x != jj {
			t.Errorf("cfC2JJ(%q) = %d, %v; cfitsio ffc2jj status %d value %d", in, x, ok2, jjSt, jj)
		}
		if ok {
			okDD++
		}
		if ok2 {
			okJJ++
		}
	}
	for _, in := range []string{"2 ", "2 x", "0xDD", "1DD", "1.5D0 ", strings.Repeat("1", 71) + "D0"} {
		if !have[in] {
			t.Errorf("c2.txt lacks the input %q", in)
		}
	}
	if okDD < 100 || okJJ < 50 {
		t.Errorf("weak corpus: %d ffc2dd and %d ffc2jj conversions succeed", okDD, okJJ)
	}
}

// TestOverwritePixelsKeywordsMatchCfitsio: OverwritePixels over the image of
// every keywords/ sample stores exactly the bytes cfitsio's fits_write_img
// stores over the same HDU (golden/keywords_overwrite.json, imgref
// overwrite) when the scaling-card rule accepts the header (the last of
// duplicate cards, as cfitsio applies it). The refusals (file unchanged)
// are every sample the rule refuses (keywords_rule.json: ErrInvalidScaling,
// where cfitsio would store with its own reading of the card), BSCALE 0
// (ErrInvalidScaling, where cfitsio stores clamped garbage, status 412),
// and a BSCALE under which the values overflow the stored range (2^-1022:
// ErrNumOverflow, where cfitsio reports NUM_OVERFLOW, 412, after storing
// the range limits).
func TestOverwritePixelsKeywordsMatchCfitsio(t *testing.T) {
	var g map[string]struct {
		Datatype int      `json:"datatype"`
		Values   []string `json:"values"`
		Status   int      `json:"status"`
		Data     string   `json:"data"`
	}
	loadJSON(t, imagesPath(t, "golden", "keywords_overwrite.json"), &g)
	samples := keywordSamples(t)
	if len(g) != len(samples) {
		t.Fatalf("keywords_overwrite.json has %d entries, keywords/ has %d samples", len(g), len(samples))
	}
	dir := t.TempDir()
	rules := loadKeywordRules(t)
	refusedN, ruleN, overflowN := 0, 0, 0
	for _, name := range samples {
		c, ok := g[name]
		if !ok {
			t.Fatalf("%s: no overwrite golden", name)
		}
		t.Run(name, func(t *testing.T) {
			p := filepath.Join(dir, name)
			copyFile(t, imagesPath(t, "keywords", name), p)
			before, err := os.ReadFile(p)
			if err != nil {
				t.Fatal(err)
			}
			e, err := OpenForEdit(p)
			if err != nil {
				t.Fatal(err)
			}
			img, err := e.Primary()
			if err != nil {
				t.Fatal(err)
			}
			var werr error
			switch c.Datatype {
			case 21:
				werr = OverwritePixels(img, parseIntVals[int16](t, c.Values, 16))
			case 31:
				werr = OverwritePixels(img, parseIntVals[int32](t, c.Values, 32))
			default:
				t.Fatalf("datatype %d", c.Datatype)
			}
			refusedByRule := rules[name].Rule != "accept"
			zero := !refusedByRule && bscaleOf(t, img) == 0
			if err := e.Close(); err != nil {
				t.Fatal(err)
			}
			after, err := os.ReadFile(p)
			if err != nil {
				t.Fatal(err)
			}
			if refusedByRule {
				ruleN++
				if !errors.Is(werr, ErrInvalidScaling) || !bytes.Equal(before, after) {
					t.Errorf("scaling-card rule %q: err %v, file changed %v; want ErrInvalidScaling, unchanged", rules[name].Rule, werr, !bytes.Equal(before, after))
				}
				return
			}
			if zero {
				refusedN++
				if !errors.Is(werr, ErrInvalidScaling) || !bytes.Equal(before, after) {
					t.Errorf("BSCALE 0 (cfitsio status %d): err %v, file changed %v; want ErrInvalidScaling, unchanged", c.Status, werr, !bytes.Equal(before, after))
				}
				return
			}
			if c.Status == 412 {
				overflowN++
				if !errors.Is(werr, ErrNumOverflow) || !bytes.Equal(before, after) {
					t.Errorf("cfitsio NUM_OVERFLOW (412): err %v, file changed %v; want ErrNumOverflow, unchanged", werr, !bytes.Equal(before, after))
				}
				return
			}
			if werr != nil || c.Status != 0 {
				t.Fatalf("OverwritePixels: %v (cfitsio writes, status %d)", werr, c.Status)
			}
			_, data := dataUnit(t, p)
			want, err := hex.DecodeString(c.Data)
			if err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(data[:len(want)], want) {
				t.Errorf("data % x, cfitsio fits_write_img % x", data[:len(want)], want)
			}
		})
	}
	if refusedN != 1 || ruleN < 30 {
		t.Errorf("%d BSCALE-0 refusals (want 1, bscale_zero), %d scaling-card rule refusals (want >= 30)", refusedN, ruleN)
	}
	// BSCALE 2^-1022 (bscale_min_normal_rounded): every nonzero value
	// overflows
	if overflowN != 1 {
		t.Errorf("%d overflow refusals, want 1", overflowN)
	}
}
