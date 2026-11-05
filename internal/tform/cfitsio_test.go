package tform

import (
	"os"
	"path/filepath"
	"strconv"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// cfitsio type codes (fitsio.h) for the binary TFORM letters.
var cfitsioBinCode = map[BinaryType]int{
	BinBit: 1, BinUint8: 11, BinLogical: 14, BinChar: 16, BinInt16: 21,
	BinInt32: 41, BinFloat32: 42, BinInt64: 81, BinFloat64: 82,
	BinComplex64: 83, BinComplex128: 163,
}

// cfitsio type codes for the ASCII TFORM letters (I -> TLONG, F/D ->
// TDOUBLE, E -> TFLOAT).
var cfitsioASCIICode = map[ASCIIType]int{
	AsciiChar: 16, AsciiInt: 41, AsciiFloatF: 82, AsciiFloatE: 42, AsciiFloatD: 82,
}

// TestParseMatchesCfitsioColumns parses the TFORM (and TDIM) of every
// table column in every committed cfitsio golden and requires the type,
// repeat, width and dimensions cfitsio reported for that column.
func TestParseMatchesCfitsioColumns(t *testing.T) {
	var goldens []string
	for _, pat := range []string{
		testutil.TestdataPath(t, "cfref", "golden", "*.json"),
		testutil.TestdataPath(t, "tables", "golden", "*.json"),
	} {
		g, _ := filepath.Glob(pat)
		goldens = append(goldens, g...)
	}
	if len(goldens) < 40 {
		t.Fatalf("found %d goldens, expected the committed cfref and tables goldens", len(goldens))
	}
	ncols := 0
	for _, g := range goldens {
		b, err := os.ReadFile(g)
		if err != nil {
			t.Fatal(err)
		}
		cf, err := testutil.DecodeCfref(b)
		if err != nil {
			t.Fatalf("%s: %v", g, err)
		}
		for _, h := range cf.HDUs {
			if h.Table == nil {
				continue
			}
			// Tile-compressed images are dumped as images, so only real
			// tables reach here.
			for _, c := range h.Table.Columns {
				ncols++
				where := filepath.Base(g) + " HDU " + itoa(h.Index) + " " + string(c.TForm)
				if h.Type == "ASCII_TBL" {
					f, err := ParseASCII(string(c.TForm))
					if err != nil {
						t.Errorf("%s: %v", where, err)
						continue
					}
					code := cfitsioASCIICode[f.Type]
					if f.Type == AsciiInt && (c.Typecode == 21 || c.Typecode == 81) {
						code = c.Typecode // cfitsio picks TSHORT/TLONG/TLONGLONG by width
					}
					// cfitsio ffasfm: Fw.d is TFLOAT up to 7 characters,
					// Ew.d TDOUBLE beyond 6 decimals.
					if f.Type == AsciiFloatF && f.Width <= 7 {
						code = 42
					}
					if f.Type == AsciiFloatE && f.Fraction > 6 {
						code = 82
					}
					if code != c.Typecode || int64(f.Width) != c.Width {
						t.Errorf("%s: type %c width %d, cfitsio typecode %d width %d", where, f.Type, f.Width, c.Typecode, c.Width)
					}
					continue
				}
				f, err := ParseBinary(string(c.TForm))
				if err != nil {
					t.Errorf("%s: %v", where, err)
					continue
				}
				elem := f.Type
				repeat := f.Repeat
				if f.Type == BinPArrayDesc32 || f.Type == BinQArrayDesc64 {
					elem = f.VarType
				}
				code := cfitsioBinCode[elem]
				if f.Type == BinPArrayDesc32 || f.Type == BinQArrayDesc64 {
					code = -code
				}
				if code != c.Typecode {
					t.Errorf("%s: typecode %d, cfitsio %d", where, code, c.Typecode)
				}
				if f.Type != BinPArrayDesc32 && f.Type != BinQArrayDesc64 && repeat != c.Repeat {
					t.Errorf("%s: repeat %d, cfitsio %d", where, repeat, c.Repeat)
				}
				// cfitsio's width: element bytes, or the substring width
				// for A (defaulting to the repeat).
				wantW := int64(elem.ElementSize())
				if elem == BinBit {
					wantW = 1
				}
				if elem == BinChar {
					wantW = f.Width
					if wantW == 0 {
						wantW = repeat
					}
					if f.Type != BinChar {
						wantW = f.VarMaxLen
					}
				}
				var tdim string
				for _, card := range h.Cards {
					if string(card.Key) == "TDIM"+itoa(c.N) && card.Str != nil {
						tdim = string(*card.Str)
					}
				}
				if f.Type == BinChar && tdim != "" {
					// cfitsio takes the string width from the first TDIM
					// axis when it is not longer than the repeat (string
					// arrays; fitscore.c ffgtbp TDIMn), also over an rAw
					// width (8A2 with TDIM (8) is one string of 8).
					if d, err := ParseDim(tdim); err == nil && len(d) > 0 && d[0] > 0 && d[0] <= repeat {
						wantW = d[0]
					}
				}
				if wantW != c.Width {
					t.Errorf("%s: width %d, cfitsio %d", where, wantW, c.Width)
				}
				if tdim != "" {
					d, err := ParseDim(tdim)
					if err != nil {
						t.Errorf("%s: TDIM %q: %v", where, tdim, err)
					} else if !equalDims(d, c.TDim) {
						t.Errorf("%s: TDIM %v, cfitsio %v", where, d, c.TDim)
					}
				}
			}
		}
	}
	t.Logf("checked %d columns in %d goldens", ncols, len(goldens))
}

func equalDims(a, b []int64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func itoa(i int) string { return strconv.Itoa(i) }
