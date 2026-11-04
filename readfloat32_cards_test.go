package fits

import (
	"bytes"
	"encoding/json"
	"errors"
	"fmt"
	"math"
	"math/big"
	"sort"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
)

// cfnullGolden is cfitsio 4.7.0's reading of a sample (testdata/imaging/
// cfnull/cards: fits_read_pixnull TDOUBLE, null = cfitsio's null flag).
type cfnullGolden struct {
	Status int        `json:"status"`
	Npix   int        `json:"npix"`
	Phys   []*float64 `json:"phys"`
}

func loadCfnull(t *testing.T, name string) cfnullGolden {
	t.Helper()
	var g cfnullGolden
	if err := json.Unmarshal(mustReadFile(t, imagingPath("cfnull/cards/"+name+".json")), &g); err != nil {
		t.Fatal(err)
	}
	return g
}

// cardVariantNames lists the card_* samples of readfloat32_expected.json.
func cardVariantNames(t *testing.T, exp map[string]rf32Expected) []string {
	t.Helper()
	var names []string
	for n := range exp {
		if strings.HasPrefix(n, "card_") {
			names = append(names, n)
		}
	}
	sort.Strings(names)
	if len(names) < 24 {
		t.Fatalf("readfloat32_expected.json lists %d card_* samples, want >= 24", len(names))
	}
	return names
}

// outcome is what one read returned, for comparing the two header paths.
func outcome(pix []float32, err error) string {
	if err != nil {
		return "error: " + err.Error()
	}
	return fmt.Sprint(f32Canon(pix))
}

// TestReadFloat32_ScalingCardVariants: the BSCALE / BZERO / BLANK card
// rule of ReadFloat32 and AsImage (testdata/imaging card_* samples) on a
// fresh Open (stored cards) and after Header() (parsed header) alike, with
// identical results on both paths:
//
//   - the LAST card of each keyword is applied: cfitsio 4.7.0
//     fits_read_img and astropy 8.0.2 both apply it (BZERO 10 then 20
//     reads 20; 'abc' then 20 reads 20; BSCALE 2 then -1 reads -1; BLANK 1
//     then 100 flags 100);
//   - that card must be a valid value card: keyword upper case in columns
//     1-8, "= " in columns 9-10, a FITS number (BLANK: an integer). Every
//     other form is an error, because the references disagree or refuse:
//     "BZERO   =10.0", "BZERO    =  10.0", "BZERO     =  10.0" (cfitsio
//     10, astropy 0), "bzero = 10.0" (cfitsio 0, astropy 10), "1.0d1"
//     (cfitsio 0, astropy 10), "1 0.0" and "10.0 x" (cfitsio's strtod
//     prefix 1 / 10, astropy refuses the file), an empty value (cfitsio 0,
//     astropy fails), a last card 'abc' (cfitsio 0, astropy fails),
//     BLANK = 100.0 / 'abc' (both silently ignore the null marker) and
//     "BLANK   =100" (cfitsio applies it, astropy does not);
//     "BZERO=10.0" (cfitsio 10, astropy 0) is an error too;
//   - valid spellings read as both references read them: 10, +1.0E1,
//     .5E1, 10., and 1.0e1 (lower-case e: both read 10);
//   - a BZERO card after END is not read (both: +10), and BLANK on float
//     BITPIX is not read even when its value is invalid (both ignore it).
//
// Truth: gen_expected.py (rule, bits; it asserts that cfitsio's and
// astropy's physical values of every accepted sample equal the documented
// ones, and records both readings under "refs"), and, checked here a
// second time, the cfnull/cards golden (cfitsio's own pixels).
func TestReadFloat32_ScalingCardVariants(t *testing.T) {
	exp := loadRF32Expected(t)
	std := loadStdimageExpected(t)
	rules := map[string]int{}
	for _, name := range cardVariantNames(t, exp) {
		e := exp[name]
		rules[e.Rule]++
		t.Run(name, func(t *testing.T) {
			var outcomes [2]string
			for i, parsed := range []bool{false, true} {
				f, err := Open(imagingPath(name))
				if err != nil {
					t.Fatal(err)
				}
				p, err := f.Primary()
				if err != nil {
					t.Fatal(err)
				}
				if parsed {
					p.Header()
				}
				got, err := ReadFloat32(p)
				outcomes[i] = outcome(got, err)
				img, ierr := p.AsImage()
				f.Close()
				if e.Rule == "error" {
					if err == nil || !strings.Contains(err.Error(), e.InvalidCard) {
						t.Errorf("parsed=%v: ReadFloat32 = %v, %v; want an error naming %s", parsed, got, err, e.InvalidCard)
					}
					if !errors.Is(ierr, ErrNotImageCompatible) || !strings.Contains(ierr.Error(), e.InvalidCard) {
						t.Errorf("parsed=%v: AsImage = %T, %v; want ErrNotImageCompatible naming %s", parsed, img, ierr, e.InvalidCard)
					}
					continue
				}
				if err != nil {
					t.Errorf("parsed=%v: ReadFloat32: %v", parsed, err)
					continue
				}
				e.Bits.Check(t, fmt.Sprintf("parsed=%v ReadFloat32 (rule %s)", parsed, e.Rule), f32Canon(got))
				if e.Rule == "float" {
					// Second truth: float32 of cfitsio's own physical values.
					g := loadCfnull(t, name)
					if g.Status != 0 || len(g.Phys) != len(got) {
						t.Fatalf("cfnull golden status %d, %d pixels; ReadFloat32 %d", g.Status, len(g.Phys), len(got))
					}
					for k, p := range g.Phys {
						if p == nil {
							if !math.IsNaN(float64(got[k])) {
								t.Errorf("parsed=%v pixel %d: %v, cfitsio null", parsed, k, got[k])
							}
						} else if float32(*p) != got[k] {
							t.Errorf("parsed=%v pixel %d: %v, cfitsio %v", parsed, k, got[k], *p)
						}
					}
				}
				s := std[name]
				if s.Kind == "reject" {
					if !errors.Is(ierr, ErrNotImageCompatible) {
						t.Errorf("parsed=%v: AsImage = %T, %v; want ErrNotImageCompatible", parsed, img, ierr)
					}
				} else if ierr != nil {
					t.Errorf("parsed=%v: AsImage: %v", parsed, ierr)
				} else {
					checkImage(t, img, s)
				}
			}
			if outcomes[0] != outcomes[1] {
				t.Errorf("fresh Open and after Header() differ:\n  fresh:  %s\n  parsed: %s", outcomes[0], outcomes[1])
			}
			cfg, cerr := DecodeConfig(bytes.NewReader(mustReadFile(t, imagingPath(name))))
			if e.Rule == "error" && cerr == nil {
				t.Errorf("DecodeConfig = %+v, nil; want an error naming %s", cfg, e.InvalidCard)
			}
			if _, _, _, ferr := ReadFrame(imagingPath(name)); (ferr != nil) != (e.Rule == "error") {
				t.Errorf("ReadFrame err = %v, want error=%v", ferr, e.Rule == "error")
			}
		})
	}
	for _, r := range []string{"float", "exact", "error"} {
		if rules[r] == 0 {
			t.Errorf("no card_* sample of rule %q", r)
		}
	}
}

// TestReadFloat32_InMemoryScalingCard: after Header().Set the in-memory
// card is the one applied (a card built in memory is judged by its typed
// value). Setting BZERO = 20 on f32_bzero_2d.fits (BZERO 10) must give
// the pixels of card_f32_bzero_dup_2d.fits (same stored floats, BZERO 10
// then 20), whose truth is cfitsio's and astropy's reading; setting a
// string is an error.
func TestReadFloat32_InMemoryScalingCard(t *testing.T) {
	want := loadRF32Expected(t)["card_f32_bzero_dup_2d.fits"]
	if want.Rule != "float" || !want.Bits.Present() {
		t.Fatalf("expected entry %+v", want)
	}
	for _, c := range []struct {
		val     any
		wantErr bool
	}{{20.0, false}, {int64(20), false}, {"20", true}} {
		f, err := Open(imagingPath("f32_bzero_2d.fits"))
		if err != nil {
			t.Fatal(err)
		}
		p, err := f.Primary()
		if err != nil {
			t.Fatal(err)
		}
		if err := p.Header().Set("BZERO", c.val, ""); err != nil {
			t.Fatal(err)
		}
		got, err := ReadFloat32(p)
		f.Close()
		if c.wantErr {
			if err == nil || !strings.Contains(err.Error(), "BZERO") {
				t.Errorf("BZERO = %#v: ReadFloat32 = %v, %v; want an error naming BZERO", c.val, got, err)
			}
			continue
		}
		if err != nil {
			t.Fatalf("BZERO = %#v: %v", c.val, err)
		}
		want.Bits.Check(t, fmt.Sprintf("BZERO = %#v", c.val), f32Canon(got))
	}
}

// TestReadFloat32_InMemoryBlankCard: an in-memory BLANK card is applied
// when it is an integer and refused otherwise. Setting BLANK = 100 on
// card_i16_blank_float_2d.fits (BSCALE 2, BLANK = 100.0, refused) must
// give the pixels of card_i16_blank_dup_2d.fits (same stored values,
// BSCALE 2, last BLANK 100), whose truth is cfitsio's and astropy's
// reading; setting BLANK = 100.0 or 2^70 (not an int64) is an error.
func TestReadFloat32_InMemoryBlankCard(t *testing.T) {
	want := loadRF32Expected(t)["card_i16_blank_dup_2d.fits"]
	if want.Rule != "exact" || !want.Bits.Present() {
		t.Fatalf("expected entry %+v", want)
	}
	for _, c := range []struct {
		val     any
		wantErr bool
	}{{int64(100), false}, {100.0, true}, {new(big.Int).Lsh(big.NewInt(1), 70), true}} {
		f, err := Open(imagingPath("card_i16_blank_float_2d.fits"))
		if err != nil {
			t.Fatal(err)
		}
		p, err := f.Primary()
		if err != nil {
			t.Fatal(err)
		}
		if err := p.Header().Set("BLANK", c.val, ""); err != nil {
			t.Fatal(err)
		}
		got, err := ReadFloat32(p)
		f.Close()
		if c.wantErr {
			if err == nil || !strings.Contains(err.Error(), "BLANK") {
				t.Errorf("BLANK = %#v: ReadFloat32 = %v, %v; want an error naming BLANK", c.val, got, err)
			}
			continue
		}
		if err != nil {
			t.Fatalf("BLANK = %#v: %v", c.val, err)
		}
		want.Bits.Check(t, fmt.Sprintf("BLANK = %#v", c.val), f32Canon(got))
	}
}

// TestReadFloat32_InMemoryBigIntScaling: an in-memory BZERO beyond int64
// (*big.Int 2^63, the uint64 convention) is applied: setting it on
// u64.fits (BZERO 2^63 as a real-valued card) gives the pixels of
// u64_intcard.fits (the same stored values, BZERO 2^63 as an integer
// literal), whose truth is astropy's stored values mapped exactly
// (gen_expected.py).
func TestReadFloat32_InMemoryBigIntScaling(t *testing.T) {
	exp := loadRF32Expected(t)
	want := exp["u64_intcard.fits"]
	if want.Rule != "exact" || !want.Bits.Present() {
		t.Fatalf("expected entry %+v", want)
	}
	f, err := Open(imagingPath("u64.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	p, err := f.Primary()
	if err != nil {
		t.Fatal(err)
	}
	if err := p.Header().Set("BZERO", new(big.Int).Lsh(big.NewInt(1), 63), ""); err != nil {
		t.Fatal(err)
	}
	if c, ok := p.Header().Get("BZERO"); !ok || c.Type != header.TypeBigInt {
		t.Fatalf("BZERO card %+v is not TypeBigInt", c)
	}
	got, err := ReadFloat32(p)
	if err != nil {
		t.Fatal(err)
	}
	want.Bits.Check(t, "BZERO = 2^63 (*big.Int)", f32Canon(got))

	// The integer mapping does not depend on the size of BZERO, so the
	// value itself is checked on float data: BZERO = 2^64 on the stored
	// floats 0.5 1.5 -2.25 NaN 100 -0 of f32_bzero_2d.fits gives
	// 2^64 + v, which is 2^64 in float64 for every |v| <= 100 (exact
	// arithmetic: the float64 spacing at 2^64 is 4096), NaN stays NaN.
	g, err := Open(imagingPath("f32_bzero_2d.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer g.Close()
	q, err := g.Primary()
	if err != nil {
		t.Fatal(err)
	}
	if err := q.Header().Set("BZERO", new(big.Int).Lsh(big.NewInt(1), 64), ""); err != nil {
		t.Fatal(err)
	}
	px, err := ReadFloat32(q)
	if err != nil {
		t.Fatal(err)
	}
	for i, v := range px {
		if i == 3 {
			if !math.IsNaN(float64(v)) {
				t.Errorf("pixel 3 = %v, want NaN", v)
			}
		} else if v != 0x1p64 {
			t.Errorf("pixel %d = %v, want 2^64", i, v)
		}
	}
}

// TestStrictFITSNumber_Grammar: strictFITSNumber accepts exactly the
// FITS 4.0 §4.2.3/§4.2.4 integer and real syntax (optional sign, digits
// with at most one decimal point and at least one digit, optional E or D
// exponent with an optionally signed digit string), plus a lower-case e
// exponent (cfitsio and astropy both read it), as long as the value is a
// finite double. The grammar is gen_expected.py's NUMBER regular
// expression; this pins the tokens on which a looser parser (ParseFloat
// alone, strtod) would differ.
func TestStrictFITSNumber_Grammar(t *testing.T) {
	ok := map[string]float64{"10": 10, "+1.0E1": 10, ".5E1": 5, "10.": 10, "1.0D1": 10, "1.0e1": 10,
		"-0": 0, "1E+05": 1e5, "1D-2": 0.01, "007": 7, "-.25": -0.25, "1E0400": 0, "1e-400": 0}
	delete(ok, "1E0400") // 1E400: finite FITS syntax, but not a finite double (below)
	for tok, want := range ok {
		if got, valid := strictFITSNumber(tok); !valid || got != want {
			t.Errorf("strictFITSNumber(%q) = %v, %v; want %v, true", tok, got, valid, want)
		}
	}
	for _, tok := range []string{"", "+", "-", ".", "E5", "1E", "1E+", "1E-", "1E1_0", "1_0", "0x1p3", "0x10",
		"1.0d1", "1d1", "inf", "Inf", "+Inf", "NaN", "nan", "1E400", "1.0.0", "1E1.5", "1EE1", "1E1D2",
		"++1", "+-1", "1E+-1", "1 0", "10 x", "'10'", "1,5", "1E400", "1e1e1"} {
		if got, valid := strictFITSNumber(tok); valid {
			t.Errorf("strictFITSNumber(%q) = %v, true; want refused", tok, got)
		}
	}
}
