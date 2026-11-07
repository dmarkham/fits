package fits

import (
	"bytes"
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// astropy is the second independent reader for table samples. Every
// numeric cell astropy reads (B I J K E D C M, their VLAs, ASCII I/F/E/D)
// must equal ours bit-for-bit. The truth is committed
// (testdata/tables/golden_astropy/<sample>.json, written by
// `astropy_dump.py --golden`, reproduced byte for byte by
// TestAstropyGoldensReproduce), so this runs in the default `go test`.
//
// A cell we report as null must be null in astropy too: NaN for a float
// cell (astropy keeps ±Inf, so a reader that turned ±Inf into null fails
// here), the column's TNULL value for an integer cell (astropy does not
// mask TNULL). ASCII null and blank fields are not compared (astropy reads
// blank numeric fields as NaN or 0; FITS 7.2.2 is pinned against the
// field text in TestTablesMatchCfitsio).
//
// Where cfitsio and astropy disagree this test is what pins our choice to
// the standard: denormals, -0 and ±Inf in E/D/C/M cells are kept as IEEE
// values (astropy) where cfitsio's null-checking read flushes them to 0 or
// flags them null; ASCII I20 "9223372036854775807" is exact (astropy)
// where cfitsio overflows through a double.

type astropyCol struct {
	HDU   int               `json:"hdu"`
	Col   int               `json:"col"`
	Name  string            `json:"name"`
	Kind  string            `json:"kind"`
	TNull json.RawMessage   `json:"tnull"`
	Rows  []json.RawMessage `json:"rows"`
	// Golden mode: a column astropy could not read, or a large column
	// stored as a digest (see astropy_dump.py).
	Error   string              `json:"error"`
	NElem   int                 `json:"nelem"`
	SHA256  string              `json:"sha256"`
	Samples [][]json.RawMessage `json:"samples"`
}

// astropyGolden is one committed golden_astropy/<sample>.json.
type astropyGolden struct {
	Astropy string       `json:"astropy"`
	Numpy   string       `json:"numpy"`
	HDUs    []string     `json:"hdus"`
	Error   string       `json:"error"`
	Patched []string     `json:"patched"`
	Columns []astropyCol `json:"columns"`
}

func astropyGoldenPath(t *testing.T, s tableSample) string {
	return testutil.TestdataPath(t, "tables", "golden_astropy", filepath.Base(s.path)+".json")
}

func loadAstropyGolden(t *testing.T, path string) astropyGolden {
	t.Helper()
	testutil.RequireFile(t, path)
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	var g astropyGolden
	if err := json.Unmarshal(b, &g); err != nil {
		t.Fatalf("%s: %v", path, err)
	}
	if g.Astropy != "8.0.2" || g.Numpy != "2.4.6" {
		t.Fatalf("%s: made with astropy %s / numpy %s, PROVENANCE.md pins 8.0.2 / 2.4.6", path, g.Astropy, g.Numpy)
	}
	return g
}

// What astropy 8.0.2 cannot read, per sample: the committed golden must
// record exactly these failures (each a substring of astropy's own error
// text, so a stale reason fails), and every other column is compared.
type astropyExpectation struct {
	colErrors map[string]string // column name -> astropy error substring
	hduClass  map[int]string    // HDU index -> astropy class (a table astropy does not see as one)
	patched   []string          // TFORM patches astropy_dump.py had to make
	// unscaledVLA: scaled VLA columns astropy mis-reads (see
	// checkAstropyUnscaled); not compared with ours.
	unscaledVLA map[string]bool
	// wholeCell: rAw columns astropy reads as one string per cell (it does
	// not split substring arrays); checked against the raw cell.
	wholeCell map[string]bool
	// rawBytes: rA columns that hold binary data (astrometry.net index
	// files pack integers and doubles into A columns). astropy returns
	// such a column as bytes when a cell is not ASCII, each cell with its
	// trailing NULs dropped (numpy) and trailing ASCII whitespace stripped
	// (chararray), dumped as the Python repr, and as str when every cell
	// decodes; checked against the raw cell, and so is our read.
	rawBytes map[string]bool
	why      string
}

// anetRawBytes are the astrometry.net index columns that hold packed
// binary data in A columns (the quads and every kdtree array).
var anetRawBytes = map[string]bool{
	"quads": true, "kdtree_lr_codes": true, "kdtree_split_codes": true, "kdtree_range_codes": true,
	"kdtree_data_codes": true, "kdtree_lr_stars": true, "kdtree_split_stars": true,
	"kdtree_range_stars": true, "kdtree_data_stars": true,
}

const anetWhy = "astrometry.net stores packed integers and doubles in A columns; astropy returns the non-ASCII ones as bytes, compared against the raw cells (checkAstropyRawBytes)"

var astropyExpect = map[string]astropyExpectation{
	"tables/anet_5205-00.fits": {rawBytes: anetRawBytes, why: anetWhy},
	"tables/anet_5206-38.fits": {rawBytes: anetRawBytes, why: anetWhy},
	"tables/anet_4107.fits":    {rawBytes: anetRawBytes, why: anetWhy},
	"tables/cf_vla_p.fits": {patched: []string{"TFORM2 1PX -> 1PB"}, unscaledVLA: map[string]bool{"VJSC": true, "VU": true},
		why: "astropy 8.0.2 rejects TFORM 1PX (VerifyError: Invalid column format) and mis-scales P columns with TSCAL/TZERO (first row only, wrapped); every other column is compared"},
	"tables/cf_vla_q.fits": {patched: []string{"TFORM2 1QX -> 1QB"}, unscaledVLA: map[string]bool{"VJSC": true, "VU": true},
		why: "astropy 8.0.2 rejects TFORM 1QX (VerifyError: Invalid column format) and mis-scales Q columns with TSCAL/TZERO (first row only, wrapped); every other column is compared"},
	"tables/cf_ascii_implied.fits": {colErrors: map[string]string{
		"IVAL": "invalid literal for int() with base 10: b'   1 2'", "FIMP": "could not convert string to float: np.bytes_(b'   1 2.5')",
		"FDOT": "could not convert string to float: np.bytes_(b'     - 5')", "EVAL": "could not convert string to float: np.bytes_(b'    1.5 E 2')",
		"DVAL": "could not convert string to float", "ISC": "invalid literal for int() with base 10: b'  - 3'"},
		why: "astropy 8.0.2 refuses blanks embedded in a field (cfitsio and FITS 7.2.5 ignore them)"},
	"tables/cf_ascii_exact.fits": {colErrors: map[string]string{
		"IU": "OverflowError: Python int too large to convert to C long", "IHI": "OverflowError", "ILO": "OverflowError",
		"DIG": "could not convert string to float: np.bytes_(b'       .')"},
		why: "astropy 8.0.2 fails on ASCII integers beyond int64 and on digitless fields; checked against the field text read as an exact decimal"},
	"tables/cf_bintable_tzero_big.fits": {colErrors: map[string]string{
		"KMAX": "bzero64", "KMAXU": "bzero64", "KHALF": "bzero64", "KNEG": "bzero64", "KNEGB": "bzero64", "KFLT": "bzero64"},
		why: "astropy 8.0.2 crashes reading a K column with an integer TZERO other than 2^63 (UnboundLocalError: bzero64); checked against stored+TZERO evaluated exactly"},
	"tables/cf_bintable_kzero.fits": {colErrors: map[string]string{"KZP": "bzero64", "KZN": "bzero64", "KZI": "bzero64"},
		why: "astropy 8.0.2 crashes reading a K column with integer TZERO 1000 (UnboundLocalError: bzero64 in fitsrec._convert_other); checked against stored+TZERO evaluated exactly"},
	"tables/cf_bintable_scaled_exact.fits": {colErrors: map[string]string{"KHALF": "bzero64", "K63P1": "bzero64", "K264": "bzero64"},
		why: "astropy 8.0.2 crashes reading a K column with a TZERO other than 0 and 2^63 (UnboundLocalError: bzero64); checked against stored*TSCAL+TZERO evaluated exactly"},
	"tables/cf_ascii_scaled_edge.fits": {colErrors: map[string]string{
		"NZ53": "invalid literal for int() with base 10: b'                   -'", "NZK": "invalid literal for int() with base 10: b'                   -'",
		"NZF": "could not convert string to float: np.bytes_(b'                  -.')", "BZ": "OverflowError: Python int too large to convert to C long"},
		why: "astropy 8.0.2 fails on digitless fields (\"-\", \"-.\") and on ASCII integers beyond int64; checked against the field text read as an exact decimal"},
	"tables/cf_bintable_scaled_special.fits": {colErrors: map[string]string{"KD18": "bzero64"},
		why: "astropy 8.0.2 crashes reading a K column with an integer TZERO other than 0 and 2^63 (UnboundLocalError: bzero64); checked against stored+TZERO evaluated exactly"},
	"tables/cf_ascii_tzero_round.fits": {colErrors: map[string]string{"BIGI": "OverflowError: Python int too large to convert to C long"},
		why: "astropy 8.0.2 fails on ASCII integers beyond int64; checked against the field text read as an exact decimal"},
	"tables/cf_bintable_tzero_round.fits": {colErrors: map[string]string{"KZ": "bzero64"},
		why: "astropy 8.0.2 crashes reading a K column with integer TZERO 1000 (UnboundLocalError: bzero64); checked against stored+TZERO evaluated exactly"},
	"tables/cf_bintable_substr.fits": {wholeCell: map[string]bool{"SUB": true},
		why: "astropy 8.0.2 reads the 20A5 substring array as one 20-character string (cfitsio splits it into 4)"},
	"tables/cf_bintable_3dtable.fits": {hduClass: map[int]string{1: "NonstandardExtHDU"},
		why: "astropy 8.0.2 does not know XTENSION '3DTABLE' (NonstandardExtHDU); cfitsio 4.7.0 reads it as BINARY_TBL (golden)"},
}

// Samples with only L/X and VLA string columns (nothing for
// astropy_dump.py to dump).
var noNumericColumns = map[string]bool{
	"tables/astropy_logical_null.fits":         true,
	"tables/astropy_vla_logical_all_zero.fits": true,
	"tables/astropy_vla_logical_null.fits":     true,
	"tables/astropy_vla_logical_pre_fix.fits":  true,
	"tables/cf_bintable_3dtable.fits":          true,
}

// astropyDivergences are single cells where astropy departs from the FITS
// standard and from cfitsio; the test checks astropy still says this and
// does not compare the cell. Key: sample, column name, row.
type astropyCell struct {
	sample, col string
	row         int
}

var astropyDivergences = map[astropyCell]struct{ astropy, why string }{
	// TNULL5 = '-99' and the field is "   -99" (right-justified). §7.2.2
	// pads TNULL with blanks to the field width ("-99   "), so this is
	// the number -99; cfitsio (prefix match) agrees. astropy strips the
	// field before matching and replaces it with 0.
	{"tables/cf_ascii_scaled.fits", "INL", 0}: {"0", "astropy matches TNULL after stripping blanks"},
	// F24.20 "1.00000000000000000001" with TZERO 2^53: the exact sum
	// 2^53 + 1.00000000000000000001 rounds to 2^53+2; astropy parses the
	// field to the double 1.0 first and then adds TZERO in float64, a tie
	// that rounds to even (2^53).
	{"tables/cf_ascii_scaled_exact.fits", "FZ53", 0}: {"4340000000000000", "astropy rounds the field to a double before adding TZERO (rounds twice)"},
	// F30.1 "0.5" and "1.5" with TZERO 2^64+2048: the exact sums round to
	// 2^64+4096; astropy adds the field's double to TZERO's double (2^64,
	// a tie rounded to even), rounding twice.
	{"tables/cf_ascii_scaled_edge.fits", "BZF", 0}: {"43f0000000000000", "astropy rounds the field and TZERO to doubles before adding (rounds twice)"},
	{"tables/cf_ascii_scaled_edge.fits", "BZF", 2}: {"43f0000000000000", "astropy rounds the field and TZERO to doubles before adding (rounds twice)"},
	// TZERO = 2^53+1 (no double holds it) under fields and stored
	// floats: the exact sums round to the values in the comments; astropy
	// adds TZERO's double 2^53 (rounding twice).
	{"tables/cf_ascii_tzero_round.fits", "TINY", 0}:    {"4340000000000000", "astropy reads 1.5E-900 as 0 and adds TZERO's double 2^53 (exact: just above the tie 2^53+1, so 2^53+2)"},
	{"tables/cf_ascii_tzero_round.fits", "TINY", 2}:    {"4340000000000000", "astropy adds 0.5 to TZERO's double 2^53 (exact 2^53+1.5: 2^53+2)"},
	{"tables/cf_ascii_tzero_round.fits", "WIDE", 0}:    {"4340000000000000", "astropy reads the field's double and adds TZERO's double 2^53 (exact: just above the tie 2^53+1, so 2^53+2)"},
	{"tables/cf_ascii_tzero_round.fits", "WIDE", 3}:    {"4340000000000000", "astropy adds 0.75 to TZERO's double 2^53 (exact 2^53+1.75: 2^53+2)"},
	{"tables/cf_ascii_tzero_round.fits", "FRAC", 2}:    {"4340000000000000", "astropy adds 0.25 to TZERO's double 2^53 (exact 2^53+1.25: 2^53+2)"},
	{"tables/cf_bintable_tzero_round.fits", "DZ53", 0}: {"4340000000000000", "astropy adds TZERO's double 2^53 to 1.0 (exact 2^53+2)"},
	{"tables/cf_bintable_tzero_round.fits", "DZ53", 1}: {"4340000000000000", "astropy adds TZERO's double 2^53 to 0.5 (exact 2^53+1.5: 2^53+2)"},
	{"tables/cf_bintable_tzero_round.fits", "EZ53", 0}: {"433fffffffffffff", "astropy adds TZERO's double 2^53 to -1.0 (exact 2^53)"},
	{"tables/cf_bintable_tzero_round.fits", "EZ53", 1}: {"4340000000000000", "astropy adds TZERO's double 2^53 to 1.0 (exact 2^53+2)"},
	{"tables/cf_bintable_tzero_round.fits", "EZ53", 2}: {"4340000000000001", "astropy adds TZERO's double 2^53 to 2.0 (exact 2^53+3, a tie: 2^53+4)"},
	// TZERO = 2^64+2049: -2.0 + TZERO = 2^64+2047 rounds to 2^64; astropy
	// adds TZERO's double 2^64+4096.
	{"tables/cf_bintable_tzero_round.fits", "DZ64", 0}: {"43f0000000000001", "astropy adds TZERO's double 2^64+4096 to -2.0 (exact 2^64+2047: 2^64)"},
	// Stored 0 with TSCAL -1 and no TZERO: astropy computes stored*TSCAL
	// only, giving -0; the FITS formula TZERO + TSCAL*stored is
	// 0 + (-0) = +0 in IEEE arithmetic, which cfitsio and we return.
	{"tables/cf_ascii_tzero_round.fits", "EUS", 3}:    {"8000000000000000", "astropy skips adding TZERO = 0 (-1.5E-900 * 2 is -0; + 0 gives +0)"},
	{"tables/cf_bintable_ranges.fits", "S_U8_IN", 0}:  {"8000000000000000", "astropy skips adding TZERO = 0 (-0 instead of +0)"},
	{"tables/cf_bintable_ranges.fits", "S_U8_LO", 1}:  {"8000000000000000", "astropy skips adding TZERO = 0 (-0 instead of +0)"},
	{"tables/cf_bintable_ranges.fits", "S_U16_IN", 0}: {"8000000000000000", "astropy skips adding TZERO = 0 (-0 instead of +0)"},
	{"tables/cf_bintable_ranges.fits", "S_U16_LO", 1}: {"8000000000000000", "astropy skips adding TZERO = 0 (-0 instead of +0)"},
	{"tables/cf_bintable_ranges.fits", "S_U32_IN", 0}: {"8000000000000000", "astropy skips adding TZERO = 0 (-0 instead of +0)"},
	{"tables/cf_bintable_ranges.fits", "S_U32_LO", 1}: {"8000000000000000", "astropy skips adding TZERO = 0 (-0 instead of +0)"},
}

func TestTablesMatchAstropy(t *testing.T) {
	used := map[astropyCell]bool{}
	for _, s := range tableSamples(t) {
		t.Run(s.name, func(t *testing.T) {
			g := loadAstropyGolden(t, astropyGoldenPath(t, s))
			exp := astropyExpect[s.name]
			if g.Error != "" {
				t.Fatalf("astropy could not open the sample: %s", g.Error)
			}
			for i, cls := range g.HDUs {
				if want, ok := exp.hduClass[i]; ok && cls != want {
					t.Errorf("astropy reads HDU %d as %s, the expectation says %s (%s)", i, cls, want, exp.why)
				}
			}
			for i, want := range exp.hduClass {
				if i >= len(g.HDUs) || g.HDUs[i] != want {
					t.Errorf("astropy HDU classes %v: HDU %d is not %s; drop the expectation (%s)", g.HDUs, i, want, exp.why)
				}
			}
			for _, p := range exp.patched {
				if !strings.Contains(strings.Join(g.Patched, "\n"), p) {
					t.Errorf("golden patched %q, want %q (%s)", g.Patched, p, exp.why)
				}
			}
			if len(g.Patched) != len(exp.patched) {
				t.Errorf("golden patched %q, expectation %q", g.Patched, exp.patched)
			}
			seenErr := map[string]bool{}
			numeric := 0
			f, err := Open(s.path)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			for _, ac := range g.Columns {
				if ac.Error != "" {
					want, ok := exp.colErrors[ac.Name]
					if !ok || !strings.Contains(ac.Error, want) {
						t.Errorf("column %s: astropy fails with %q, not the expected %q", ac.Name, ac.Error, want)
					}
					seenErr[ac.Name] = true
					continue
				}
				numeric++
				if exp.wholeCell[ac.Name] {
					seenErr["whole:"+ac.Name] = true
					checkAstropyWholeCell(t, f, ac)
					continue
				}
				if exp.unscaledVLA[ac.Name] {
					seenErr["unscaled:"+ac.Name] = true
					checkAstropyUnscaled(t, s, ac)
					continue
				}
				if exp.rawBytes[ac.Name] {
					seenErr["raw:"+ac.Name] = true
					checkAstropyRawBytes(t, f, ac)
					continue
				}
				checkAstropyColumn(t, s.name, f, ac, used)
			}
			for name := range exp.wholeCell {
				if !seenErr["whole:"+name] {
					t.Errorf("column %s: not in the astropy golden; drop the wholeCell expectation", name)
				}
			}
			for name := range exp.unscaledVLA {
				if !seenErr["unscaled:"+name] {
					t.Errorf("column %s: not in the astropy golden; drop the unscaledVLA expectation", name)
				}
			}
			for name := range exp.rawBytes {
				if !seenErr["raw:"+name] {
					t.Errorf("column %s: not in the astropy golden; drop the rawBytes expectation", name)
				}
			}
			for name := range exp.colErrors {
				if !seenErr[name] {
					t.Errorf("column %s: astropy now reads it (no error in the golden); drop the expectation (%s)", name, exp.why)
				}
			}
			if len(g.Columns) == 0 && !noNumericColumns[s.name] {
				t.Fatalf("astropy reported no numeric columns")
			}
			if len(g.Columns) > 0 && noNumericColumns[s.name] {
				t.Errorf("astropy reports %d numeric columns; drop the noNumericColumns entry", len(g.Columns))
			}
			_ = numeric
		})
	}
	for k := range astropyDivergences {
		if !used[k] {
			t.Errorf("astropyDivergences %v names no cell any golden has", k)
		}
	}
}

// checkAstropyWholeCell pins a wholeCell expectation: astropy's one string
// per row is the whole cell (cut at NUL, trailing blanks removed).
func checkAstropyWholeCell(t *testing.T, f *File, ac astropyCol) {
	t.Helper()
	bt := mustBinTable(t, f, ac.HDU)
	c, err := bt.column(ac.Col)
	if err != nil {
		t.Fatal(err)
	}
	for r, raw := range ac.Rows {
		got := astropyRow(t, raw, "str")
		row, err := bt.RowBytes(int64(r))
		if err != nil {
			t.Fatal(err)
		}
		want := cellString(row[c.byteOffset : c.byteOffset+c.byteSize])
		if len(got) != 1 || got[0] != want {
			t.Errorf("%s row %d: astropy %q, the whole cell is %q: astropy now splits rAw; drop the wholeCell expectation", ac.Name, r, got, want)
		}
	}
}

// checkAstropyRawBytes pins a rawBytes expectation against the raw cell
// bytes: astropy's cell is the raw cell without trailing NULs and ASCII
// whitespace, as the repr of a bytes object when the column did not
// decode as ASCII and as the string otherwise, whose rstrip also takes
// the separators \x1c..\x1f (astropy 8.0.2 / numpy 2.4.6, checked on
// every cell of the fixtures); ours is the FITS string
// of the cell (cut at the first NUL, trailing blanks removed, §7.3.3.1).
func checkAstropyRawBytes(t *testing.T, f *File, ac astropyCol) {
	t.Helper()
	bt := mustBinTable(t, f, ac.HDU)
	c, err := bt.column(ac.Col)
	if err != nil {
		t.Fatal(err)
	}
	arr, err := ReadStringArrayColumn(bt, c.Index)
	if err != nil {
		t.Fatal(err)
	}
	if len(arr) != len(ac.Rows) {
		t.Fatalf("%s: %d rows, astropy %d", ac.Name, len(arr), len(ac.Rows))
	}
	for r, raw := range ac.Rows {
		got := astropyRow(t, raw, "str")
		row, err := bt.RowBytes(int64(r))
		if err != nil {
			t.Fatal(err)
		}
		cell := row[c.byteOffset : c.byteOffset+c.byteSize]
		// bytes.rstrip strips ASCII whitespace; str.rstrip also the
		// separators \x1c..\x1f (Unicode whitespace).
		noNUL := bytes.TrimRight(cell, "\x00")
		want := string(bytes.TrimRight(noNUL, " \t\n\r\x0b\x0c\x1c\x1d\x1e\x1f"))
		if len(got) == 1 && (strings.HasPrefix(got[0], "b'") || strings.HasPrefix(got[0], "b\"")) {
			want = pyBytesRepr(bytes.TrimRight(noNUL, " \t\n\r\x0b\x0c"))
		}
		if len(got) != 1 || got[0] != want {
			t.Errorf("%s row %d: astropy %q, the raw cell %q gives %q: astropy changed; review the rawBytes expectation", ac.Name, r, got, cell, want)
		}
		if len(arr[r]) != 1 || arr[r][0] != cellString(cell) {
			t.Errorf("%s row %d: we read %q, the raw cell %q is the string %q", ac.Name, r, arr[r], cell, cellString(cell))
		}
	}
}

// pyBytesRepr is Python's repr of a bytes object: single quotes unless
// the bytes hold a single quote and no double quote; backslash, the
// quote, \t \n \r escaped by name, every other byte outside 0x20..0x7e as
// \xNN.
func pyBytesRepr(b []byte) string {
	q := byte('\'')
	if bytes.IndexByte(b, '\'') >= 0 && bytes.IndexByte(b, '"') < 0 {
		q = '"'
	}
	var sb strings.Builder
	sb.WriteByte('b')
	sb.WriteByte(q)
	for _, c := range b {
		switch {
		case c == q || c == '\\':
			sb.WriteByte('\\')
			sb.WriteByte(c)
		case c == '\t':
			sb.WriteString("\\t")
		case c == '\n':
			sb.WriteString("\\n")
		case c == '\r':
			sb.WriteString("\\r")
		case c < 0x20 || c >= 0x7f:
			fmt.Fprintf(&sb, "\\x%02x", c)
		default:
			sb.WriteByte(c)
		}
	}
	sb.WriteByte(q)
	return sb.String()
}

// checkAstropyUnscaled pins the astropy bug behind an unscaledVLA
// expectation: for a P/Q column with TSCAL/TZERO astropy 8.0.2 scales the
// FIRST row only, truncated and wrapped into the stored integer type, and
// returns every other row's stored values unscaled. Both halves are
// checked against cfitsio (its physical read and its read with scaling
// off), so an astropy change is noticed.
func checkAstropyUnscaled(t *testing.T, s tableSample, ac astropyCol) {
	t.Helper()
	gc := loadGoldenFile(t, s.golden).HDUs[ac.HDU].Table.Column(ac.Name)
	if gc == nil {
		t.Fatalf("%s: no cfitsio golden column", ac.Name)
	}
	wrap := func(f float64) int64 {
		switch -gc.Typecode {
		case cfShort:
			return int64(int16(int64(f)))
		case cfLong:
			return int64(int32(int64(f)))
		}
		t.Fatalf("%s: typecode %d not modelled", ac.Name, gc.Typecode)
		return 0
	}
	rows := gc.RowCells(t)
	if len(rows) != len(ac.Rows) {
		t.Fatalf("%s: %d rows, cfitsio %d", ac.Name, len(ac.Rows), len(rows))
	}
	for r, raw := range ac.Rows {
		got := astropyRow(t, raw, ac.Kind)
		if len(got) != len(rows[r].Raw) || len(got) != len(rows[r].Phys) {
			t.Errorf("%s row %d: astropy %d elements, cfitsio %d", ac.Name, r, len(got), len(rows[r].Raw))
			continue
		}
		for k := range got {
			want := rows[r].Raw[k].V
			if r == 0 {
				want = wrap(float64(rows[r].Phys[k].V))
			}
			if got[k] != strconv.FormatInt(want, 10) {
				t.Errorf("%s row %d elem %d: astropy %s, want %d (first row scaled and wrapped, others stored): astropy changed; review the unscaledVLA expectation", ac.Name, r, k, got[k], want)
			}
		}
	}
}

func checkAstropyColumn(t *testing.T, sample string, f *File, ac astropyCol, used map[astropyCell]bool) {
	t.Helper()
	where := sample + " HDU " + strconv.Itoa(ac.HDU) + " col " + strconv.Itoa(ac.Col) + " " + ac.Name
	h, err := f.HDU(ac.HDU)
	if err != nil {
		t.Fatal(err)
	}
	// Our values per row, as strings in the same encoding as the dump.
	var ours [][]string
	var valid [][]bool
	isASCII := false
	var col Column
	switch v := h.(type) {
	case *BinaryTableHDU:
		cols, err := v.Columns()
		if err != nil {
			t.Fatal(err)
		}
		col = cols[ac.Col-1]
		ours, valid, err = goCells(v, col, ac.Kind)
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
	case *ASCIITableHDU:
		isASCII = true
		ours, valid, err = goASCIICells(v, ac.Col, ac.Kind)
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
	default:
		t.Fatalf("%s: HDU is %T", where, h)
	}
	// A cell we call null must be null in astropy: NaN for floats, the
	// TNULL value for integers (astropy does not mask it, and applies
	// TSCAL/TZERO to it: TNULL*TSCAL + TZERO, an integer for the unsigned
	// conventions, a float64 otherwise).
	nullAs := func() string {
		if ac.Kind == "complex" {
			return "<complex has no null>"
		}
		var tn int64
		if err := json.Unmarshal(ac.TNull, &tn); err != nil || len(ac.TNull) == 0 || string(ac.TNull) == "null" {
			if ac.Kind == "float" {
				return "NaN"
			}
			return "<no TNULL in astropy>"
		}
		if ac.Kind == "float" {
			return fbits(float64(tn)*col.Scale + col.Zero)
		}
		return strconv.FormatInt(tn+int64(col.Zero), 10)
	}
	nullOK := func(a string) bool { return a == nullAs() || (ac.Kind == "float" && a == "NaN") }
	if ac.SHA256 != "" {
		checkAstropyDigest(t, where, ac, ours, valid, nullAs(), isASCII)
		return
	}
	if len(ours) != len(ac.Rows) {
		t.Errorf("%s: %d rows, astropy %d", where, len(ours), len(ac.Rows))
		return
	}
	compared := 0
	for r, raw := range ac.Rows {
		want := astropyRow(t, raw, ac.Kind)
		if len(want) != len(ours[r]) {
			t.Errorf("%s row %d: %d elements, astropy %d", where, r, len(ours[r]), len(want))
			continue
		}
		for k := range want {
			key := astropyCell{sample, ac.Name, r}
			if d, ok := astropyDivergences[key]; ok {
				if used != nil {
					used[key] = true
				}
				if want[k] != d.astropy {
					t.Errorf("%s row %d: astropy now reads %s, not %s (%s); drop the divergence entry", where, r, want[k], d.astropy, d.why)
				}
				continue
			}
			if isASCII {
				if !valid[r][k] {
					continue // null or blank field: see the top of the file
				}
			} else if !valid[r][k] {
				if ok := nullOK(want[k]); !ok {
					t.Errorf("%s row %d elem %d: we read null, astropy %s (want %s)", where, r, k, want[k], nullAs())
				}
				continue
			}
			compared++
			if ac.Kind == "str" {
				// Trailing blanks are not significant (§7.3.3.1); astropy
				// keeps them on the elements of a string array.
				want[k] = strings.TrimRight(want[k], " ")
			}
			if ours[r][k] != want[k] {
				t.Errorf("%s row %d elem %d: %s, astropy %s", where, r, k, ours[r][k], want[k])
			}
		}
	}
	if !isASCII && compared == 0 && len(ac.Rows) > 0 && rowsHaveElems(ac.Rows) {
		t.Errorf("%s: no cell compared", where)
	}
}

// checkAstropyDigest compares a digested astropy column (astropy_dump.py
// digest) with ours: the sha256 of our values in the same canonical
// encoding (a null cell as astropy holds it), plus every sampled element.
func checkAstropyDigest(t *testing.T, where string, ac astropyCol, ours [][]string, valid [][]bool, nullAs string, isASCII bool) {
	t.Helper()
	if isASCII {
		t.Fatalf("%s: digested ASCII column: blank-field handling is not modelled", where)
	}
	var flat []string
	h := sha256.New()
	var b [8]byte
	put := func(u uint64) {
		binary.BigEndian.PutUint64(b[:], u)
		h.Write(b[:])
	}
	enc := func(s string) {
		switch ac.Kind {
		case "float":
			put(canonBits(t, s))
		case "complex":
			re, im, _ := strings.Cut(s, ",")
			put(canonBits(t, re))
			put(canonBits(t, im))
		case "uint":
			u, err := strconv.ParseUint(s, 10, 64)
			if err != nil {
				t.Fatalf("%s: %q: %v", where, s, err)
			}
			put(u)
		default:
			i, err := strconv.ParseInt(s, 10, 64)
			if err != nil {
				t.Fatalf("%s: %q: %v", where, s, err)
			}
			put(uint64(i))
		}
	}
	for r := range ours {
		put(uint64(len(ours[r])))
		for k, v := range ours[r] {
			if !valid[r][k] {
				v = nullAs
			}
			flat = append(flat, v)
			enc(v)
		}
	}
	if len(flat) != ac.NElem {
		t.Errorf("%s: %d elements, astropy %d", where, len(flat), ac.NElem)
		return
	}
	bad := 0
	for _, smp := range ac.Samples {
		var i int
		if err := json.Unmarshal(smp[0], &i); err != nil {
			t.Fatal(err)
		}
		want := astropyRow(t, json.RawMessage("["+string(smp[1])+"]"), ac.Kind)[0]
		if flat[i] != want {
			bad++
			if bad <= 5 {
				t.Errorf("%s elem %d: %s, astropy %s", where, i, flat[i], want)
			}
		}
	}
	if got := hex.EncodeToString(h.Sum(nil)); got != ac.SHA256 {
		t.Errorf("%s: sha256 of our %d values %s, astropy %s (%d of %d sampled elements differ)", where, len(flat), got, ac.SHA256, bad, len(ac.Samples))
	}
}

// canonBits is the canonical float64 bit pattern of a canonFloat string
// (NaN as 0x7ff8000000000000).
func canonBits(t *testing.T, s string) uint64 {
	t.Helper()
	if s == "NaN" {
		return 0x7ff8000000000000
	}
	u, err := strconv.ParseUint(s, 16, 64)
	if err != nil {
		t.Fatalf("bad float bits %q: %v", s, err)
	}
	return u
}

func rowsHaveElems(rows []json.RawMessage) bool {
	for _, r := range rows {
		if string(r) != "[]" {
			return true
		}
	}
	return false
}

// astropyRow decodes one dumped row into canonical strings.
func astropyRow(t *testing.T, raw json.RawMessage, kind string) []string {
	t.Helper()
	switch kind {
	case "str":
		var xs []string
		if err := json.Unmarshal(raw, &xs); err != nil {
			t.Fatal(err)
		}
		return xs
	case "float":
		var xs []string
		if err := json.Unmarshal(raw, &xs); err != nil {
			t.Fatal(err)
		}
		for i, x := range xs {
			xs[i] = canonFloat(t, x)
		}
		return xs
	case "complex":
		var xs [][2]string
		if err := json.Unmarshal(raw, &xs); err != nil {
			t.Fatal(err)
		}
		out := make([]string, len(xs))
		for i, x := range xs {
			out[i] = canonFloat(t, x[0]) + "," + canonFloat(t, x[1])
		}
		return out
	}
	var xs []json.Number
	if err := json.Unmarshal(raw, &xs); err != nil {
		t.Fatal(err)
	}
	out := make([]string, len(xs))
	for i, x := range xs {
		out[i] = x.String()
	}
	return out
}

func canonFloat(t *testing.T, s string) string {
	t.Helper()
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		t.Fatalf("bad astropy float %q: %v", s, err)
	}
	return fbits(f)
}

func fbits(f float64) string {
	if math.IsNaN(f) {
		return "NaN"
	}
	return strconv.FormatUint(math.Float64bits(f), 16)
}

func goCells(h *BinaryTableHDU, c Column, kind string) ([][]string, [][]bool, error) {
	nrows := int(h.NumRows())
	ours := make([][]string, nrows)
	valid := make([][]bool, nrows)
	isVar := c.Type == ColVarArray
	if isVar && (c.bin.VarType == 'C' || c.bin.VarType == 'M') {
		kind = "complex" // astropy cannot tell for all-empty rows
	}
	switch kind {
	case "str":
		arr, err := ReadStringArrayColumn(h, c.Index)
		if err != nil {
			return nil, nil, err
		}
		for r := 0; r < nrows; r++ {
			ours[r] = arr[r]
			valid[r] = make([]bool, len(arr[r]))
			for k := range valid[r] {
				valid[r][k] = true
			}
		}
		return ours, valid, nil
	case "complex":
		if isVar {
			v, err := ReadVarComplexColumn[complex128](h, c.Index)
			if err != nil {
				return nil, nil, err
			}
			for r := 0; r < nrows; r++ {
				for _, z := range v.At(r) {
					ours[r] = append(ours[r], fbits(real(z))+","+fbits(imag(z)))
					valid[r] = append(valid[r], true)
				}
			}
			return ours, valid, nil
		}
		v, err := ReadComplexColumn[complex128](h, c.Index)
		if err != nil {
			return nil, nil, err
		}
		rep := int(c.Repeat)
		for r := 0; r < nrows; r++ {
			for _, z := range v[r*rep : (r+1)*rep] {
				ours[r] = append(ours[r], fbits(real(z))+","+fbits(imag(z)))
				valid[r] = append(valid[r], true)
			}
		}
		return ours, valid, nil
	case "float":
		return goNumCells[float64](h, c, isVar, func(x float64) string { return fbits(x) })
	case "uint":
		return goNumCells[uint64](h, c, isVar, func(x uint64) string { return strconv.FormatUint(x, 10) })
	}
	return goNumCells[int64](h, c, isVar, func(x int64) string { return strconv.FormatInt(x, 10) })
}

func goNumCells[T Numeric](h *BinaryTableHDU, c Column, isVar bool, str func(T) string) ([][]string, [][]bool, error) {
	nrows := int(h.NumRows())
	ours := make([][]string, nrows)
	valid := make([][]bool, nrows)
	if isVar {
		v, ok, err := ReadVarColumnMasked[T](h, c.Index)
		if err != nil {
			return nil, nil, err
		}
		for r := 0; r < nrows; r++ {
			for i, x := range v.At(r) {
				ours[r] = append(ours[r], str(x))
				valid[r] = append(valid[r], ok.At(r)[i])
			}
		}
		return ours, valid, nil
	}
	v, ok, err := ReadVectorColumnMasked[T](h, c.Index)
	if err != nil {
		return nil, nil, err
	}
	rep := int(c.Repeat)
	for r := 0; r < nrows; r++ {
		for k := 0; k < rep; k++ {
			ours[r] = append(ours[r], str(v[r*rep+k]))
			valid[r] = append(valid[r], ok[r*rep+k])
		}
	}
	return ours, valid, nil
}

// goASCIICells reads an ASCII column. All-blank fields are marked invalid
// for the comparison: cfitsio (and we) read them as 0, astropy as NaN/0.
// It also returns each row's field text.
func goASCIICells(h *ASCIITableHDU, col int, kind string) ([][]string, [][]bool, error) {
	cols, err := h.Columns()
	if err != nil {
		return nil, nil, err
	}
	c := cols[col-1]
	rowW := int(h.rec.shape[0])
	nrows := int(h.NumRows())
	buf := make([]byte, rowW*nrows)
	if err := h.rec.file.br.ReadRange(h.rec.dataStart, buf); err != nil {
		return nil, nil, err
	}
	fields := make([]string, nrows)
	for r := range fields {
		fields[r] = string(buf[r*rowW+c.byteOffset : r*rowW+c.byteOffset+c.byteSize])
	}
	blank := func(r int) bool { return strings.TrimSpace(fields[r]) == "" }
	ours := make([][]string, nrows)
	valid := make([][]bool, nrows)
	if kind == "str" {
		v, ok, err := ReadASCIIColumnMasked[string](h, col)
		if err != nil {
			return nil, nil, err
		}
		for r := range v {
			ours[r] = []string{v[r]}
			valid[r] = []bool{ok[r]}
		}
		return ours, valid, nil
	}
	if kind == "float" {
		v, ok, err := ReadASCIIColumnMasked[float64](h, col)
		if err != nil {
			return nil, nil, err
		}
		for r := range v {
			ours[r] = []string{fbits(v[r])}
			valid[r] = []bool{ok[r] && !blank(r)}
		}
		return ours, valid, nil
	}
	v, ok, err := ReadASCIIColumnMasked[int64](h, col)
	if err != nil {
		return nil, nil, err
	}
	for r := range v {
		ours[r] = []string{strconv.FormatInt(v[r], 10)}
		valid[r] = []bool{ok[r] && !blank(r)}
	}
	return ours, valid, nil
}

// TestSpecialValuesSample proves the IEEE special values (±Inf, -0, the
// smallest subnormals, NaN, the largest finite values) in E, D, C and M
// cells, fixed and VLA, against a file written by astropy itself
// (testdata/tables/astropy_only/gen_special.py) and astropy's committed
// read of it. ±Inf, -0 and subnormals must be valid IEEE values, NaN the
// null of E/D; a reader that turned ±Inf into null, or flushed -0 or a
// subnormal to +0, fails. The golden itself must hold the generator's
// bit patterns, so the sample cannot lose them.
func TestSpecialValuesSample(t *testing.T) {
	path := testutil.TestdataPath(t, "tables", "astropy_only", "astropy_special.fits")
	testutil.RequireFile(t, path)
	g := loadAstropyGolden(t, testutil.TestdataPath(t, "tables", "golden_astropy", "astropy_special.fits.json"))
	inf := math.Inf(1)
	f32 := []float64{inf, -inf, math.Copysign(0, -1), float64(math.SmallestNonzeroFloat32), -float64(math.SmallestNonzeroFloat32), math.NaN(), 1.5, math.MaxFloat32}
	f64 := []float64{inf, -inf, math.Copysign(0, -1), math.SmallestNonzeroFloat64, -math.SmallestNonzeroFloat64, math.NaN(), 1.5, math.MaxFloat64}
	n := len(f32)
	enc := func(xs ...float64) string {
		var b []string
		for _, x := range xs {
			b = append(b, fbits(x))
		}
		return strings.Join(b, ",")
	}
	want := map[string][][]string{}
	for r := 0; r < n; r++ {
		want["E1"] = append(want["E1"], []string{enc(f32[r])})
		want["D1"] = append(want["D1"], []string{enc(f64[r])})
		want["E2"] = append(want["E2"], []string{enc(f32[r]), enc(f32[(r+3)%n])})
		want["C1"] = append(want["C1"], []string{enc(f32[r], f32[(r+2)%n])})
		want["M1"] = append(want["M1"], []string{enc(f64[r], f64[(r+2)%n])})
		var ve, vd, vc, vm []string
		for k := 0; k < r%4; k++ {
			ve = append(ve, enc(f32[k]))
		}
		for k := r; k < n; k++ {
			vd = append(vd, enc(f64[k]))
		}
		for k := 0; k < r%3; k++ {
			vc = append(vc, enc(f32[k], f32[(k+5)%n]))
		}
		for k := r; k < min(n, r+2); k++ {
			vm = append(vm, enc(f64[k], f64[(k+1)%n]))
		}
		want["VE"] = append(want["VE"], ve)
		want["VD"] = append(want["VD"], vd)
		want["VC"] = append(want["VC"], vc)
		want["VM"] = append(want["VM"], vm)
	}
	f, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	if len(g.Columns) != len(want) {
		t.Fatalf("golden has %d columns, generator writes %d", len(g.Columns), len(want))
	}
	for _, ac := range g.Columns {
		w := want[ac.Name]
		if ac.Error != "" || len(ac.Rows) != len(w) {
			t.Fatalf("%s: golden error %q, %d rows", ac.Name, ac.Error, len(ac.Rows))
		}
		for r, raw := range ac.Rows {
			if got := strings.Join(astropyRow(t, raw, ac.Kind), ";"); got != strings.Join(w[r], ";") {
				t.Errorf("%s row %d: astropy golden %s, generator wrote %s", ac.Name, r, got, strings.Join(w[r], ";"))
			}
		}
		checkAstropyColumn(t, "astropy_special.fits", f, ac, nil)
	}
	// Only the NaN cells of E/D are null; ±Inf, -0 and subnormals are values.
	bt := mustBinTable(t, f, 1)
	for _, c := range []int{1, 2} {
		_, ok, err := ReadColumnMasked[float64](bt, c)
		if err != nil {
			t.Fatal(err)
		}
		for r, v := range ok {
			if v != (r != 5) {
				t.Errorf("column %d row %d: valid = %v (only the NaN row 5 is null)", c, r, v)
			}
		}
	}
}
