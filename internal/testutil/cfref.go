package testutil

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// This file decodes the JSON written by testdata/cfref/cfref.c ("cfref dump").
// The field meanings are documented at the top of cfref.c; in short, every
// value here is what the system libcfitsio returned, and a JSON null in a
// value array means cfitsio set its null flag for that element. Arrays of
// more than 1024 values are stored as exact digests, and mid-size arrays
// of the testdata/cfref goldens once each in golden/pool (cfref_values.go):
// image pixels and flattened table columns are Values, compared with
// Check / Diff or recovered with Resolve, never indexed blindly.

// CfrefFile is one "cfref dump" result.
type CfrefFile struct {
	File    string     `json:"file"`    // base name of the dumped file
	Cfitsio string     `json:"cfitsio"` // libcfitsio version, e.g. "4.7.0"
	NHDU    int        `json:"nhdu"`
	HDUs    []CfrefHDU `json:"hdus"`
}

// CfrefHDU describes one HDU. Index is 0-based (fits.File.HDU(i)).
type CfrefHDU struct {
	Index int    `json:"index"`
	Type  string `json:"type"` // "IMAGE", "ASCII_TBL", "BINARY_TBL"
	// From fits_get_hduaddrll. DataEnd is cfitsio's "dataend": the offset of
	// the next HDU, i.e. the end of the data unit including block padding.
	HeaderStart int64 `json:"header_start"`
	DataStart   int64 `json:"data_start"`
	DataEnd     int64 `json:"data_end"`
	// From fits_verify_chksum: 1 = correct, 0 = keyword absent, -1 = wrong.
	Checksum struct {
		Datasum int `json:"datasum"`
		HDUSum  int `json:"hdusum"`
	} `json:"checksum"`
	NKeys  int         `json:"nkeys"`
	Cards  []CfrefCard `json:"cards"`
	Image  *CfrefImage `json:"image"`
	Table  *CfrefTable `json:"table"`
	Errors []string    `json:"errors"` // cfitsio errors met while dumping this HDU
}

// CfrefCard is one header card (END excluded, as cfitsio counts keys).
type CfrefCard struct {
	N int `json:"n"` // 1-based card number
	// Raw is the card text from fits_read_record; cfitsio strips trailing
	// blanks, so compare against strings.TrimRight(card, " ").
	Raw     Latin1 `json:"raw"`
	Key     Latin1 `json:"key"` // fits_get_keyname (HIERARCH prefix removed)
	Class   int    `json:"class"`
	Value   Latin1 `json:"value"` // fits_parse_value value string (quotes kept)
	Comment Latin1 `json:"comment"`
	Type    string `json:"type"` // fits_get_keytype: C, L, I, F, X; "" if no value
	// Decoded values (only the one matching Type is set; I sets Int and Float).
	Str     *Latin1 `json:"str"` // fits_read_key_longstr: CONTINUE resolved
	Int     *int64  `json:"int"`
	Float   *Float  `json:"float"`
	Logical *bool   `json:"logical"`
	// Non-empty when cfitsio reported a problem with this card.
	Error      string `json:"error"`
	ParseError string `json:"parse_error"`
	TypeError  string `json:"type_error"`
	ReadError  string `json:"read_error"`
}

// cardFields is CfrefCard as a pooled card file holds it (no n, no raw).
type cardFields struct {
	Key        Latin1  `json:"key"`
	Class      int     `json:"class"`
	Value      Latin1  `json:"value"`
	Comment    Latin1  `json:"comment"`
	Type       string  `json:"type"`
	Str        *Latin1 `json:"str"`
	Int        *int64  `json:"int"`
	Float      *Float  `json:"float"`
	Logical    *bool   `json:"logical"`
	ParseError string  `json:"parse_error"`
	TypeError  string  `json:"type_error"`
	ReadError  string  `json:"read_error"`
}

// UnmarshalJSON accepts a full card object, or (cfref "dump -pool") a
// card stored once in testdata/cfref/golden/pool: the golden then holds
// just the card text, and pool/card-<sha256 of that text>.json holds
// cfitsio's parse of it (every field but n and raw). N is the card's
// position in its HDU (set by DecodeCfref). cfref.c pools a card only when
// its parse depends on its 80 bytes alone, and refuses to write a second
// parse for the same text.
func (c *CfrefCard) UnmarshalJSON(b []byte) error {
	b = bytes.TrimSpace(b)
	if len(b) > 0 && b[0] == '"' {
		var raw Latin1
		if err := raw.UnmarshalJSON(b); err != nil {
			return err
		}
		sum := sha256.Sum256([]byte(raw))
		path := filepath.Join(cfrefPoolDir(), "card-"+hex.EncodeToString(sum[:])+".json")
		pb, err := ReadGolden(path)
		if err != nil {
			return fmt.Errorf("testutil: pooled card %q: %w", string(raw), err)
		}
		dec := json.NewDecoder(bytes.NewReader(pb))
		dec.DisallowUnknownFields()
		var f cardFields
		if err := dec.Decode(&f); err != nil {
			return fmt.Errorf("testutil: pooled card file %s: %w", path, err)
		}
		*c = CfrefCard{Raw: raw, Key: f.Key, Class: f.Class, Value: f.Value, Comment: f.Comment,
			Type: f.Type, Str: f.Str, Int: f.Int, Float: f.Float, Logical: f.Logical,
			ParseError: f.ParseError, TypeError: f.TypeError, ReadError: f.ReadError}
		return nil
	}
	type plain CfrefCard
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	var p plain
	if err := dec.Decode(&p); err != nil {
		return err
	}
	if p.N < 1 {
		return fmt.Errorf("testutil: harness card without n: %s", b)
	}
	*c = CfrefCard(p)
	return nil
}

// CfrefImage holds an image HDU's parameters and every pixel.
type CfrefImage struct {
	Bitpix     int     `json:"bitpix"`
	Equivtype  int     `json:"equivtype"` // fits_get_img_equivtype
	Compressed bool    `json:"compressed"`
	Naxis      int     `json:"naxis"`
	Naxes      []int64 `json:"naxes"`
	Npix       int64   `json:"npix"`
	Bscale     Float   `json:"bscale"`
	Bzero      Float   `json:"bzero"`
	Blank      *int64  `json:"blank"`
	// Pixels in FITS order (NAXIS1 fastest). Phys is TDOUBLE with scaling;
	// I64 is TLONGLONG (integer equivalent types except uint64); U64 is
	// TULONGLONG (uint64 equivalent type); Raw is the stored integer with
	// scaling disabled (only when BSCALE/BZERO are non-default). Past
	// DIGEST_MIN pixels each is stored as a digest: compare with Check /
	// CheckFunc (cfref_values.go), which stay exact.
	Phys Values[NullFloat]  `json:"phys"`
	I64  Values[NullInt64]  `json:"i64"`
	U64  Values[NullUint64] `json:"u64"`
	Raw  Values[NullInt64]  `json:"raw"`
	// PhysF32 is Phys rounded to float32 (a C cast in the harness), written
	// only for a BITPIX -32 image whose Phys is a digest, so a float32
	// reading can be resolved against it (PhysAsFloat32).
	PhysF32 Values[NullFloat] `json:"phys_f32"`
}

// PhysAsFloat32 returns the TDOUBLE pixels rounded to float32 (as
// float64), element by element: from an explicit Phys, or by resolving
// cand (a float32 reading widened to float64) against PhysF32.
func (im *CfrefImage) PhysAsFloat32(t testing.TB, where string, cand func() []NullFloat) []NullFloat {
	t.Helper()
	if !im.Phys.Digested() {
		out := append([]NullFloat(nil), im.Phys.Values(t)...)
		for k := range out {
			if out[k].Valid {
				out[k].V = float64(float32(out[k].V))
			}
		}
		return out
	}
	if !im.PhysF32.Present() {
		t.Fatalf("%s: digested pixels without the float32-rounded digest (BITPIX %d)", where, im.Bitpix)
	}
	return im.PhysF32.Resolve(t, where, cand)
}

// CfrefTable holds a binary or ASCII table.
type CfrefTable struct {
	NRows   int64         `json:"nrows"`
	NCols   int           `json:"ncols"`
	RowLen  int64         `json:"rowlen"` // NAXIS1
	Columns []CfrefColumn `json:"columns"`
	// ASCII tables only: each row's raw bytes (fits_read_tblbytes).
	RawRows []Latin1 `json:"raw_rows"`
}

// CfrefColumn is one table column with every row's values.
type CfrefColumn struct {
	N        int     `json:"n"`
	TType    *Latin1 `json:"ttype"`
	TForm    Latin1  `json:"tform"`
	TUnit    *Latin1 `json:"tunit"`
	TDisp    *Latin1 `json:"tdisp"`
	TNull    *Latin1 `json:"tnull"`
	TScal    Float   `json:"tscal"`
	TZero    Float   `json:"tzero"`
	Typecode int     `json:"typecode"` // fits_get_coltypell (negative = VLA)
	Eqtype   int     `json:"eqtype"`   // fits_get_eqcoltypell
	Repeat   int64   `json:"repeat"`
	Width    int64   `json:"width"`
	TDim     []int64 `json:"tdim"`  // binary tables: fits_read_tdimll
	TBCol    *int64  `json:"tbcol"` // ASCII tables
	// The harness writes every row explicitly ("rows"), except that a
	// fixed-width column with more than DIGEST_MIN elements in all its
	// rows is stored flattened as digests ("cells"). Use RowCells for
	// per-row access (it fails the test on a digested column) and Flat for
	// whole-column checks in either form.
	rows  []CfrefCell
	cells *CfrefCells
}

// UnmarshalJSON decodes a column strictly (unknown fields are errors).
func (c *CfrefColumn) UnmarshalJSON(b []byte) error {
	type plain CfrefColumn
	var v struct {
		plain
		Rows  []CfrefCell `json:"rows"`
		Cells *CfrefCells `json:"cells"`
	}
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	if err := dec.Decode(&v); err != nil {
		return err
	}
	if (v.Rows == nil) == (v.Cells == nil) {
		return fmt.Errorf("testutil: column %d must have exactly one of rows / cells", v.N)
	}
	*c = CfrefColumn(v.plain)
	c.rows, c.cells = v.Rows, v.Cells
	return nil
}

// Digested reports whether the column is stored as flattened digests.
func (c *CfrefColumn) Digested() bool { return c.cells != nil }

// CfrefCells is a column's values for all rows, flattened in row order:
// row r, element j is index r*PerRow+j. PerRow is the strings per cell for
// an A column and the repeat count otherwise; -1 when rows differ in
// length (a variable-length column given explicitly).
type CfrefCells struct {
	PerRow int                `json:"per_row"`
	Str    Values[NullString] `json:"str"`
	Log    Values[NullBool]   `json:"log"`
	Bits   Values[Bit]        `json:"bits"`
	Cplx   Values[Cplx]       `json:"cplx"`
	Phys   Values[NullFloat]  `json:"phys"`
	I64    Values[NullInt64]  `json:"i64"`
	U64    Values[NullUint64] `json:"u64"`
	Raw    Values[NullInt64]  `json:"raw"`
}

// RowCells returns the explicit per-row cells; a column stored as a
// digest (Cells) fails the test: compare it through Flat instead.
func (c *CfrefColumn) RowCells(t testing.TB) []CfrefCell {
	t.Helper()
	if c.cells != nil {
		t.Fatalf("testutil: column %d (%d values) is stored as a digest; compare it through Flat()", c.N, c.Flat().Len())
	}
	return c.rows
}

// Flat returns the column's values for all rows, flattened, in either
// form. For explicit rows a kind is Present only if some row carries it.
func (c *CfrefColumn) Flat() *CfrefCells {
	if c.cells != nil {
		return c.cells
	}
	var f CfrefCells
	var str []NullString
	var log []NullBool
	var bits []Bit
	var cplx []Cplx
	var phys []NullFloat
	var i64, raw []NullInt64
	var u64 []NullUint64
	f.PerRow = -2
	for _, r := range c.rows {
		n := len(r.Str) + len(r.Log) + len(r.Bits) + len(r.Cplx) + len(r.Phys)
		switch {
		case f.PerRow == -2:
			f.PerRow = n
		case f.PerRow != n:
			f.PerRow = -1
		}
		str = append(str, r.Str...)
		log = append(log, r.Log...)
		for _, b := range r.Bits {
			bits = append(bits, Bit(b))
		}
		for _, x := range r.Cplx {
			cplx = append(cplx, Cplx(x))
		}
		phys = append(phys, r.Phys...)
		i64 = append(i64, r.I64...)
		u64 = append(u64, r.U64...)
		raw = append(raw, r.Raw...)
	}
	if f.PerRow == -2 {
		f.PerRow = 0
	}
	has := func(p func(r *CfrefCell) bool) bool {
		for i := range c.rows {
			if p(&c.rows[i]) {
				return true
			}
		}
		return false
	}
	if has(func(r *CfrefCell) bool { return r.Str != nil }) {
		f.Str = ExplicitValues(str)
	}
	if has(func(r *CfrefCell) bool { return r.Log != nil }) {
		f.Log = ExplicitValues(log)
	}
	if has(func(r *CfrefCell) bool { return r.Bits != nil }) {
		f.Bits = ExplicitValues(bits)
	}
	if has(func(r *CfrefCell) bool { return r.Cplx != nil }) {
		f.Cplx = ExplicitValues(cplx)
	}
	if has(func(r *CfrefCell) bool { return r.Phys != nil }) {
		f.Phys = ExplicitValues(phys)
	}
	if has(func(r *CfrefCell) bool { return r.I64 != nil }) {
		f.I64 = ExplicitValues(i64)
	}
	if has(func(r *CfrefCell) bool { return r.U64 != nil }) {
		f.U64 = ExplicitValues(u64)
	}
	if has(func(r *CfrefCell) bool { return r.Raw != nil }) {
		f.Raw = ExplicitValues(raw)
	}
	return &f
}

// Len is the number of flattened values of the column's value kind.
func (f *CfrefCells) Len() int {
	for _, n := range []int{f.Str.Len(), f.Log.Len(), f.Bits.Len(), f.Cplx.Len(), f.Phys.Len()} {
		if n > 0 {
			return n
		}
	}
	return 0
}

// CfrefCell is one row of one column. Exactly one of the value groups is
// set, chosen by the column type: Str (A), Log (L), Bits (X), Cplx (C/M),
// or Phys plus I64/U64/Raw (numeric). For a VLA column Desc is
// [length, heap offset] and the values cover exactly length elements; an
// empty VLA cell has N == 0 and no values.
type CfrefCell struct {
	Desc *[2]int64    `json:"desc"`
	N    *int         `json:"n"`
	Str  []NullString `json:"str"`
	Log  []NullBool   `json:"log"`
	Bits []int        `json:"bits"`
	Cplx [][2]Float   `json:"cplx"`
	Phys []NullFloat  `json:"phys"`
	I64  []NullInt64  `json:"i64"`
	U64  []NullUint64 `json:"u64"`
	Raw  []NullInt64  `json:"raw"`
}

// Latin1 is a harness string. cfref.c escapes every byte >= 0x7F as
// \u00XX meaning that byte, so Latin1 maps each decoded rune back to a
// single byte: string(l) is the exact byte sequence cfitsio returned.
type Latin1 string

func (l *Latin1) UnmarshalJSON(b []byte) error {
	var s string
	if err := json.Unmarshal(b, &s); err != nil {
		return err
	}
	out := make([]byte, 0, len(s))
	for _, r := range s {
		if r > 0xFF {
			return fmt.Errorf("testutil: harness string has rune %U > U+00FF", r)
		}
		out = append(out, byte(r))
	}
	*l = Latin1(out)
	return nil
}

// Float is a float64 that also accepts the harness's "NaN", "Inf", "-Inf".
type Float float64

func (f *Float) UnmarshalJSON(b []byte) error {
	v, null, err := parseFloat(b)
	if err != nil {
		return err
	}
	if null {
		return fmt.Errorf("testutil: unexpected null float")
	}
	*f = Float(v)
	return nil
}

func parseFloat(b []byte) (v float64, null bool, err error) {
	b = bytes.TrimSpace(b)
	switch string(b) {
	case "null":
		return 0, true, nil
	case `"NaN"`:
		return math.NaN(), false, nil
	case `"Inf"`:
		return math.Inf(1), false, nil
	case `"-Inf"`:
		return math.Inf(-1), false, nil
	}
	v, err = strconv.ParseFloat(string(b), 64)
	if err != nil {
		return 0, false, fmt.Errorf("testutil: bad float %s: %w", b, err)
	}
	return v, false, nil
}

// NullFloat is a value with cfitsio's null flag (Valid false = null).
type NullFloat struct {
	V     float64
	Valid bool
}

func (n *NullFloat) UnmarshalJSON(b []byte) error {
	v, null, err := parseFloat(b)
	if err != nil {
		return err
	}
	*n = NullFloat{V: v, Valid: !null}
	return nil
}

// NullInt64 is an int64 with cfitsio's null flag.
type NullInt64 struct {
	V     int64
	Valid bool
}

func (n *NullInt64) UnmarshalJSON(b []byte) error {
	s := string(bytes.TrimSpace(b))
	if s == "null" {
		*n = NullInt64{}
		return nil
	}
	v, err := strconv.ParseInt(s, 10, 64)
	if err != nil {
		return fmt.Errorf("testutil: bad int64 %s: %w", s, err)
	}
	*n = NullInt64{V: v, Valid: true}
	return nil
}

// NullUint64 is a uint64 with cfitsio's null flag.
type NullUint64 struct {
	V     uint64
	Valid bool
}

func (n *NullUint64) UnmarshalJSON(b []byte) error {
	s := string(bytes.TrimSpace(b))
	if s == "null" {
		*n = NullUint64{}
		return nil
	}
	v, err := strconv.ParseUint(s, 10, 64)
	if err != nil {
		return fmt.Errorf("testutil: bad uint64 %s: %w", s, err)
	}
	*n = NullUint64{V: v, Valid: true}
	return nil
}

// NullBool is a FITS logical: Valid false means undefined (byte 0x00).
type NullBool struct {
	V     bool
	Valid bool
}

func (n *NullBool) UnmarshalJSON(b []byte) error {
	switch string(bytes.TrimSpace(b)) {
	case "null":
		*n = NullBool{}
	case "1":
		*n = NullBool{V: true, Valid: true}
	case "0":
		*n = NullBool{Valid: true}
	default:
		return fmt.Errorf("testutil: bad logical %s", b)
	}
	return nil
}

// NullString is a string with cfitsio's null flag.
type NullString struct {
	V     string
	Valid bool
}

func (n *NullString) UnmarshalJSON(b []byte) error {
	if string(bytes.TrimSpace(b)) == "null" {
		*n = NullString{}
		return nil
	}
	var l Latin1
	if err := l.UnmarshalJSON(b); err != nil {
		return err
	}
	*n = NullString{V: string(l), Valid: true}
	return nil
}

// Card returns the first card with the given keyword name, or nil.
func (h *CfrefHDU) Card(key string) *CfrefCard {
	for i := range h.Cards {
		if string(h.Cards[i].Key) == key {
			return &h.Cards[i]
		}
	}
	return nil
}

// Column returns the column whose TTYPE is name, or nil.
func (tb *CfrefTable) Column(name string) *CfrefColumn {
	for i := range tb.Columns {
		if c := tb.Columns[i].TType; c != nil && string(*c) == name {
			return &tb.Columns[i]
		}
	}
	return nil
}

// CfrefDir is <repo>/testdata/cfref.
func CfrefDir(t testing.TB) string {
	t.Helper()
	return TestdataPath(t, "cfref")
}

// CfrefGoldenPath is the committed golden for a testdata sample, e.g.
// CfrefGoldenPath(t, "testprog.std") -> testdata/cfref/golden/testprog.std.json.
// The decompressed twins of the compressed read corpus (comp_fzun_<case>
// .fits, not committed: testdata/gen_twins.sh regenerates them) have their
// goldens in testdata/cfref/golden/twins.
func CfrefGoldenPath(t testing.TB, name string) string {
	t.Helper()
	name = strings.TrimSuffix(filepath.Base(name), ".json")
	if strings.HasPrefix(name, "comp_fzun_") {
		return filepath.Join(CfrefDir(t), "golden", "twins", name+".json")
	}
	return filepath.Join(CfrefDir(t), "golden", name+".json")
}

// LoadCfrefGolden decodes the committed cfitsio golden for the testdata
// sample called name (base name, with or without the .json suffix). A
// missing golden fails the test.
func LoadCfrefGolden(t testing.TB, name string) *CfrefFile {
	t.Helper()
	p := CfrefGoldenPath(t, name)
	RequireFile(t, p)
	b, err := os.ReadFile(p)
	if err != nil {
		t.Fatalf("read golden %s: %v", p, err)
	}
	cf, err := DecodeCfref(b)
	if err != nil {
		t.Fatalf("decode golden %s: %v", p, err)
	}
	return cf
}

// DecodeCfref parses "cfref dump" JSON strictly (unknown fields are errors,
// so a harness change cannot be silently ignored by the decoder).
func DecodeCfref(b []byte) (*CfrefFile, error) {
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	var cf CfrefFile
	if err := dec.Decode(&cf); err != nil {
		return nil, err
	}
	if len(cf.HDUs) != cf.NHDU {
		return nil, fmt.Errorf("nhdu %d but %d HDUs dumped", cf.NHDU, len(cf.HDUs))
	}
	for h := range cf.HDUs {
		for i := range cf.HDUs[h].Cards {
			c := &cf.HDUs[h].Cards[i]
			if c.N == 0 {
				c.N = i + 1 // a pooled card: its position
			} else if c.N != i+1 {
				return nil, fmt.Errorf("HDU %d card %d has n %d", h, i+1, c.N)
			}
		}
	}
	return &cf, nil
}

// CfrefDump is the cfitsio harness's reading ("cfref dump") of path,
// typically a file the Go code under test just wrote: the recorded
// reading of these exact bytes in tier 1, the harness run live in tier 2
// (Run). If cfitsio cannot open the file at all the test fails with
// cfitsio's message; use CfrefDumpErr to test for rejection.
func CfrefDump(t testing.TB, path string) *CfrefFile {
	t.Helper()
	cf, err := CfrefDumpErr(t, path)
	if err != nil {
		t.Fatalf("cfitsio harness on %s: %v", path, err)
	}
	return cf
}

// committedSampleGolden is the golden name of path when path is a sample
// committed directly in testdata/ that has a cfref golden, "" otherwise.
func committedSampleGolden(t testing.TB, path string) string {
	t.Helper()
	abs, err := filepath.Abs(path)
	if err != nil || filepath.Dir(abs) != TestdataPath(t) {
		return ""
	}
	if _, err := os.Stat(CfrefGoldenPath(t, abs)); err != nil {
		return ""
	}
	return filepath.Base(abs)
}

// CfrefDumpErr is CfrefDump returning cfitsio's refusal as an error.
func CfrefDumpErr(t testing.TB, path string) (*CfrefFile, error) {
	t.Helper()
	if g := committedSampleGolden(t, path); g != "" {
		// A committed sample: its committed golden is the harness's dump
		// of it (TestCfrefDumpReproducesGoldens, tier 2).
		return LoadCfrefGolden(t, g), nil
	}
	r := Run(t, ToolCfref, Call{Args: []string{"dump", path}, Inputs: []string{path}})
	if r.Exit != 0 {
		return nil, fmt.Errorf("cfref dump failed (exit status %d): %s", r.Exit, strings.TrimSpace(string(r.Stderr)))
	}
	cf, err := DecodeCfref(r.Stdout)
	if err != nil {
		return nil, fmt.Errorf("decode cfref output: %w", err)
	}
	return cf, nil
}
