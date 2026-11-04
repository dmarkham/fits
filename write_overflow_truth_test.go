package fits

import (
	"bytes"
	"encoding/hex"
	"errors"
	"fmt"
	"math"
	"math/big"
	"os"
	"path/filepath"
	"strconv"
	"testing"

	"github.com/dmarkham/fits/internal/bigendian"
)

// writeOverflowGolden is testdata/images/golden/write_overflow.json: per
// case, each value written alone by cfitsio's fits_write_img with the
// case's BSCALE/BZERO (gen.py WRITE_OVERFLOW), with cfitsio's status
// (412 NUM_OVERFLOW or 0) and the stored element (hex).
type writeOverflowGolden map[string]struct {
	GoType   string   `json:"gotype"`
	Imgtype  int      `json:"imgtype"`
	Datatype int      `json:"datatype"`
	Bscale   *string  `json:"bscale"`
	Bzero    *string  `json:"bzero"`
	Values   []string `json:"values"`
	Status   []int    `json:"status"`
	Data     []string `json:"data"`
}

// silentOverflow lists the integer values cfitsio stores wrongly WITHOUT
// reporting NUM_OVERFLOW, which the library refuses as well: (case, value
// index) -> what cfitsio stores.
var silentOverflow = map[string]map[int]int64{
	// v = -2^63 with BZERO 1000: the stored value would be -2^63-1000;
	// cfitsio's double arithmetic absorbs the 1000 and stores -2^63,
	// which reads back as -2^63+1000.
	"i64_bzero1000": {0: math.MinInt64},
	// v = 2^62 and 2^62-1 with BSCALE 0.5, BZERO 0.25: the scaled value
	// is 2^63 as a double; cfitsio's range test (dvalue > DLONGLONG_MAX,
	// which is 2^63 as a double) lets it through to an out-of-range
	// (LONGLONG) cast that stores INT64_MIN: read back as -2^62+0.25.
	"i64_bscale_half": {2: math.MinInt64, 3: math.MinInt64},
}

// overflowCase turns a golden entry into a write-matrix case holding the
// values idx.
func overflowCase(t *testing.T, name string, g writeOverflowGolden, idx []int) imgCase {
	t.Helper()
	e := g[name]
	c := imgCase{Name: name, GoType: e.GoType, CfImgtype: e.Imgtype, CfDatatype: e.Datatype, Naxes: []int64{int64(len(idx))}}
	parse := func(s *string) *float64 {
		if s == nil {
			return nil
		}
		v, err := strconv.ParseFloat(*s, 64)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		return &v
	}
	c.Bscale, c.Bzero = parse(e.Bscale), parse(e.Bzero)
	for _, i := range idx {
		c.Values = append(c.Values, e.Values[i])
	}
	return c
}

// storedMiss is |physical - v|, exact, for a stored element (integer
// BITPIX; physical = stored*BSCALE + BZERO) and the value v written.
func storedMiss(t *testing.T, bitpix int, raw []byte, bscale, bzero float64, v string) *big.Float {
	t.Helper()
	var s int64
	switch bitpix {
	case 8:
		s = int64(raw[0])
	case 16:
		s = int64(bigendian.Int16(raw))
	case 32:
		s = int64(bigendian.Int32(raw))
	case 64:
		s = bigendian.Int64(raw)
	}
	const prec = 256
	phys := new(big.Float).SetPrec(prec).SetInt64(s)
	phys.Mul(phys, new(big.Float).SetPrec(prec).SetFloat64(bscale))
	phys.Add(phys, new(big.Float).SetPrec(prec).SetFloat64(bzero))
	want, ok := new(big.Float).SetPrec(prec).SetString(v)
	if !ok {
		t.Fatalf("value %q", v)
	}
	return phys.Sub(phys, want).Abs(phys)
}

// overflowRefused reports whether the library must refuse value i of case
// name, from cfitsio's behaviour only: cfitsio reported NUM_OVERFLOW; or
// it stored, without reporting, an integer the documented silentOverflow
// entry describes (asserted against the golden, and shown to be wrong:
// it reads back further from the value than one BSCALE step); or, for
// float BITPIX, it stored an infinity for a finite value.
func overflowRefused(t *testing.T, name string, g writeOverflowGolden, i int) bool {
	t.Helper()
	e := g[name]
	raw, err := hex.DecodeString(e.Data[i])
	if err != nil {
		t.Fatal(err)
	}
	switch e.Status[i] {
	case 412:
		return true
	case 0:
	default:
		t.Fatalf("%s[%d]: cfitsio status %d", name, i, e.Status[i])
	}
	if e.Imgtype < 0 {
		in, err := strconv.ParseFloat(e.Values[i], 64)
		if err != nil {
			t.Fatal(err)
		}
		var stored float64
		if e.Imgtype == -32 {
			stored = float64(bigendian.Float32(raw))
		} else {
			stored = bigendian.Float64(raw)
		}
		return !math.IsInf(in, 0) && !math.IsNaN(in) && math.IsInf(stored, 0)
	}
	want, silent := silentOverflow[name][i]
	if !silent {
		return false
	}
	var bscale, bzero float64 = 1, 0
	if e.Bscale != nil {
		bscale, _ = strconv.ParseFloat(*e.Bscale, 64)
	}
	if e.Bzero != nil {
		bzero, _ = strconv.ParseFloat(*e.Bzero, 64)
	}
	if got := bigendian.Int64(raw); e.Imgtype != 64 || got != want {
		t.Fatalf("%s[%d]: cfitsio stored %d (BITPIX %d), silentOverflow records %d", name, i, got, e.Imgtype, want)
	}
	if storedMiss(t, e.Imgtype, raw, bscale, bzero, e.Values[i]).Cmp(big.NewFloat(math.Abs(bscale))) <= 0 {
		t.Fatalf("%s[%d]: cfitsio's stored value reads back within one BSCALE step of %s: not a silent overflow", name, i, e.Values[i])
	}
	return true
}

// TestWriteImageOverflowMatchesCfitsio: WriteImage returns ErrNumOverflow
// and writes nothing for a value that cannot be stored under the image's
// BITPIX and BSCALE/BZERO, instead of clamping it silently, and for every
// other value stores exactly cfitsio's bytes. Truth: each value written
// alone by cfitsio (golden/write_overflow.json): refused exactly where
// cfitsio reports NUM_OVERFLOW (412) on the lower and upper boundary of
// BITPIX 8/16/32/64 (incl. a value of exactly lo-0.5, unsigned and int8
// conventions with negative and non-integer BSCALE), plus the values
// cfitsio stores wrongly without reporting it (silentOverflow; float
// BITPIX scaled to infinity). Writing a whole case at once is refused when
// any value is.
func TestWriteImageOverflowMatchesCfitsio(t *testing.T) {
	var g writeOverflowGolden
	loadJSON(t, imagesPath(t, "golden", "write_overflow.json"), &g)
	if len(g) < 18 {
		t.Fatalf("write_overflow.json has %d cases", len(g))
	}
	dir := t.TempDir()
	refusedN, storedN := 0, 0
	for name, e := range g {
		if len(e.Status) != len(e.Values) || len(e.Data) != len(e.Values) {
			t.Fatalf("%s: malformed golden", name)
		}
		t.Run(name, func(t *testing.T) {
			anyRefused := false
			for i := range e.Values {
				refuse := overflowRefused(t, name, g, i)
				anyRefused = anyRefused || refuse
				c := overflowCase(t, name, g, []int{i})
				p := filepath.Join(dir, fmt.Sprintf("%s_%d.fits", name, i))
				err, n := writeOverflowCase(t, p, c)
				if refuse {
					refusedN++
					if !errors.Is(err, ErrNumOverflow) || n != 0 {
						t.Errorf("value %s: err %v, %d HDUs; want ErrNumOverflow and none (cfitsio status %d, stored %s)",
							e.Values[i], err, n, e.Status[i], e.Data[i])
					}
					continue
				}
				storedN++
				if err != nil || n != 1 {
					t.Errorf("value %s: err %v, %d HDUs (cfitsio stores %s)", e.Values[i], err, n, e.Data[i])
					continue
				}
				_, data := dataUnit(t, p)
				if got := hex.EncodeToString(data[:len(e.Data[i])/2]); got != e.Data[i] {
					t.Errorf("value %s: stored %s, cfitsio %s", e.Values[i], got, e.Data[i])
				}
			}
			all := make([]int, len(e.Values))
			for i := range all {
				all[i] = i
			}
			err, n := writeOverflowCase(t, filepath.Join(dir, name+"_all.fits"), overflowCase(t, name, g, all))
			if anyRefused != (err != nil) || anyRefused && (!errors.Is(err, ErrNumOverflow) || n != 0) {
				t.Errorf("all values at once: err %v, %d HDUs (a value overflows: %v)", err, n, anyRefused)
			}
		})
	}
	if refusedN < 25 || storedN < 60 {
		t.Errorf("%d refusals, %d stored values", refusedN, storedN)
	}
}

// TestOverwritePixelsOverflowMatchesCfitsio: OverwritePixels over an HDU
// with the case's scaling returns ErrNumOverflow and leaves the file's
// bytes unchanged for exactly the values WriteImage refuses (cfitsio's
// NUM_OVERFLOW and the silent overflows, golden/write_overflow.json), and
// stores cfitsio's bytes for every other value.
func TestOverwritePixelsOverflowMatchesCfitsio(t *testing.T) {
	var g writeOverflowGolden
	loadJSON(t, imagesPath(t, "golden", "write_overflow.json"), &g)
	dir := t.TempDir()
	for name, e := range g {
		t.Run(name, func(t *testing.T) {
			first := -1
			for i := range e.Values {
				if !overflowRefused(t, name, g, i) {
					first = i
					break
				}
			}
			if first < 0 {
				t.Fatal("no storable value to create the HDU with")
			}
			for i := range e.Values {
				refuse := overflowRefused(t, name, g, i)
				p := filepath.Join(dir, fmt.Sprintf("%s_%d.fits", name, i))
				if err, n := writeOverflowCase(t, p, overflowCase(t, name, g, []int{first})); err != nil || n != 1 {
					t.Fatalf("WriteImage: err %v, %d HDUs", err, n)
				}
				before, err := os.ReadFile(p)
				if err != nil {
					t.Fatal(err)
				}
				f, err := OpenForEdit(p)
				if err != nil {
					t.Fatal(err)
				}
				img, err := f.Primary()
				if err != nil {
					t.Fatal(err)
				}
				_, over := caseData(t, overflowCase(t, name, g, []int{i}), false)
				werr := over(img)
				if err := f.Close(); err != nil {
					t.Fatal(err)
				}
				after, err := os.ReadFile(p)
				if err != nil {
					t.Fatal(err)
				}
				if refuse {
					if !errors.Is(werr, ErrNumOverflow) || !bytes.Equal(before, after) {
						t.Errorf("value %s: err %v, file changed %v; want ErrNumOverflow, unchanged", e.Values[i], werr, !bytes.Equal(before, after))
					}
					continue
				}
				if werr != nil {
					t.Errorf("value %s: %v", e.Values[i], werr)
					continue
				}
				_, data := dataUnit(t, p)
				if got := hex.EncodeToString(data[:len(e.Data[i])/2]); got != e.Data[i] {
					t.Errorf("value %s: stored %s, cfitsio %s", e.Values[i], got, e.Data[i])
				}
			}
		})
	}
}

// writeOverflowCase writes case c with WriteImage into a new file at p and
// returns WriteImage's error and the number of HDUs written.
func writeOverflowCase(t *testing.T, p string, c imgCase) (error, int) {
	t.Helper()
	w, _ := caseData(t, c, false)
	f, err := Create(p)
	if err != nil {
		t.Fatal(err)
	}
	_, werr := w(f, caseHeader(t, c), c.Naxes)
	n := f.NumHDU()
	if err := f.Close(); err != nil && werr == nil {
		t.Fatal(err)
	}
	return werr, n
}
