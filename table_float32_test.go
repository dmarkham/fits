package fits

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// float32 targets against cfitsio's own float32 read (TFLOAT).
//
// testdata/tables/golden_tfloat/<sample>.json is cfitsio 4.7.0's
// fits_read_colnull(TFLOAT) of every numeric column of every table sample
// (cf_tfloat.c). cfitsio converts an unscaled integer or float straight to
// float32 (one rounding, the correctly rounded value) and everything else
// through its double result ((float)(stored*scale + zero)). Our float32
// targets get the float32 nearest the exact value on the exact paths and
// the float32 of the float64 arithmetic on the float64 path (TSCAL != 1,
// fractional TZERO), which is what cfitsio computes there.
//
// So every cell must be bit-identical to cfitsio's TFLOAT, except:
//   - cfitsio's null-checking read flags ±Inf null and flushes -0 and
//     subnormals to +0 (cfitsioFlushed; those cells are proven against
//     astropy), and its ASCII TNULL prefix match (asciiTNULLPrefixDivergences);
//   - beyond the float32 range cfitsio stores ±FLT_MAX with NUM_OVERFLOW,
//     we give ±Inf (the documented difference, as for images);
//   - cells where cfitsio's TFLOAT is the float32 of its own TDOUBLE read
//     (cfref golden) and ours is not: cfitsio rounded through a double
//     (a sum with a TZERO no double holds, a long ASCII field; its NaN
//     for a field of more than ~308 digits counts as such);
//     the value we return there is checked against the exact model in
//     TestTablesMatchCfitsio (checkTarget, float32 = the float32 nearest
//     the exact value). Samples with such cells are pinned below.
//
// Anything else is a failure: a float32 read that rounds twice where
// cfitsio rounds once (the old float32(float64) rule on 1K 2^60+2^36+1:
// ours 0x5d800000, cfitsio 0x5d800001) is caught here.
var tfloatThroughDouble = map[string]string{
	"tables/cf_bintable_tzero_big.fits":   "TZERO held as a double",
	"tables/cf_bintable_tzero_round.fits": "TZERO 2^53+1 / 2^64+2049 held as a double; K + 1000 in double precision",
	"tables/cf_ascii_scaled_edge.fits":    "ASCII fields and TZERO beyond 2^53 read through doubles",
	"tables/cf_ascii_tzero_round.fits":    "ASCII fields and TZERO 2^53+1 read through doubles",
}

// tfloatOverflow lists the columns where cfitsio's TFLOAT read returns
// NUM_OVERFLOW (values beyond FLT_MAX, which it sets to ±FLT_MAX).
var tfloatOverflow = map[string]bool{"tables/cf_bintable_fixed.fits D1": true}

type tfloatGolden struct {
	File string `json:"file"`
	HDUs []struct {
		Index   int `json:"index"`
		Columns []struct {
			N        int         `json:"n"`
			TType    string      `json:"ttype"`
			Error    int         `json:"error"`
			Overflow bool        `json:"overflow"`
			Rows     [][]*string `json:"rows"`
		} `json:"columns"`
	} `json:"hdus"`
}

func TestTablesFloat32MatchCfitsio(t *testing.T) {
	dir := testutil.TestdataPath(t, "tables")
	n := 0
	for _, s := range tableSamples(t) {
		if filepath.Dir(s.path) != dir {
			continue // the top-level cfitsio test files have no TFLOAT golden
		}
		n++
		t.Run(s.name, func(t *testing.T) {
			gp := filepath.Join(dir, "golden_tfloat", filepath.Base(s.path)+".json")
			testutil.RequireFile(t, gp)
			b, err := os.ReadFile(gp)
			if err != nil {
				t.Fatal(err)
			}
			var g tfloatGolden
			if err := json.Unmarshal(b, &g); err != nil {
				t.Fatal(err)
			}
			cf := loadGoldenFile(t, s.golden)
			f, err := Open(s.path)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			through := 0
			for _, gh := range g.HDUs {
				var ch *testutil.CfrefHDU
				for i := range cf.HDUs {
					if cf.HDUs[i].Index == gh.Index {
						ch = &cf.HDUs[i]
					}
				}
				if ch == nil || ch.Table == nil {
					t.Fatalf("HDU %d: not a table in the cfref golden", gh.Index)
				}
				h, err := f.HDU(gh.Index)
				if err != nil {
					t.Fatal(err)
				}
				for _, gc := range gh.Columns {
					where := fmt.Sprintf("%s HDU %d col %d %s", s.name, gh.Index, gc.N, gc.TType)
					if gc.Error != 0 {
						t.Errorf("%s: cfitsio TFLOAT read failed with status %d", where, gc.Error)
						continue
					}
					if gc.Overflow != tfloatOverflow[s.name+" "+gc.TType] {
						t.Errorf("%s: cfitsio NUM_OVERFLOW = %v, tfloatOverflow says %v", where, gc.Overflow, !gc.Overflow)
					}
					var rc testutil.CfrefColumn
					for _, c := range ch.Table.Columns {
						if c.N == gc.N {
							rc = c
						}
					}
					var td []testutil.NullFloat
					for _, cell := range rc.RowCells(t) {
						td = append(td, cell.Phys...)
					}
					var rd targetReaders[float32]
					var ascii bool
					switch v := h.(type) {
					case *BinaryTableHDU:
						if rc.Typecode < 0 {
							rd = varReaders[float32](v, gc.N)
						} else {
							rd = fixedReaders[float32](v, gc.N, false)
						}
					case *ASCIITableHDU:
						rd, ascii = asciiReaders[float32](v, gc.N), true
					}
					got, valid, err := rd.masked()
					if err != nil {
						t.Errorf("%s: %v", where, err)
						continue
					}
					var want []*string
					for _, r := range gc.Rows {
						want = append(want, r...)
					}
					if len(want) != len(got) || len(td) != len(got) {
						t.Errorf("%s: %d values, cfitsio TFLOAT %d, TDOUBLE %d", where, len(got), len(want), len(td))
						continue
					}
					for k, w := range want {
						gb := math.Float32bits(got[k])
						switch {
						case w == nil && !valid[k]:
						case w == nil:
							if cfitsioFlushed(float64(got[k]), true, testutil.NullFloat{}) {
								continue // ±Inf, which cfitsio flags null
							}
							if ascii && asciiTNULLPrefix(s.name, gh.Index, gc.TType, k) {
								continue
							}
							t.Errorf("%s elem %d: ours %v, cfitsio TFLOAT null", where, k, got[k])
						case !valid[k]:
							t.Errorf("%s elem %d: ours null, cfitsio TFLOAT %s", where, k, *w)
						default:
							wb64, perr := strconv.ParseUint(*w, 16, 32)
							if perr != nil {
								t.Fatalf("%s elem %d: bad golden bits %q", where, k, *w)
							}
							wb := uint32(wb64)
							wf := math.Float32frombits(wb)
							switch {
							case wb == gb:
							case gc.Overflow && math.IsInf(float64(got[k]), 0) && math.Abs(float64(wf)) == math.MaxFloat32 && math.Signbit(float64(wf)) == math.Signbit(float64(got[k])):
							case wb == 0 && cfitsioFlushed(float64(got[k]), true, testutil.NullFloat{V: 0, Valid: true}):
							case td[k].Valid && (wb == math.Float32bits(float32(td[k].V)) || (math.IsNaN(td[k].V) && math.IsNaN(float64(wf)))):
								// NaN: cfitsio's digit accumulator overflows on
								// fields of more than ~308 digits (Inf/Inf).
								through++
							default:
								t.Errorf("%s elem %d: ours %08x (%v), cfitsio TFLOAT %08x (%v), cfitsio TDOUBLE %v: not explained", where, k, gb, got[k], wb, wf, td[k])
							}
						}
					}
				}
			}
			why, listed := tfloatThroughDouble[s.name]
			switch {
			case listed && through == 0:
				t.Errorf("cfitsio's TFLOAT now agrees on every cell (%q no longer holds); drop the entry", why)
			case !listed && through > 0:
				t.Errorf("%d cells where cfitsio's TFLOAT rounds through its double and ours does not; pin the sample in tfloatThroughDouble", through)
			case listed:
				t.Logf("pinned: %d cells where cfitsio rounds through a double (%s)", through, why)
			}
		})
	}
	if n == 0 {
		t.Fatal("no table samples")
	}
}

func asciiTNULLPrefix(sample string, hdu int, col string, row int) bool {
	for _, d := range asciiTNULLPrefixDivergences {
		if d.sample == sample && d.hdu == hdu && d.col == col && d.row == row {
			return true
		}
	}
	return false
}
