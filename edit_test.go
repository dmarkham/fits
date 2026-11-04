package fits

import (
	"bytes"
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// Tests for OpenForEdit + Flush. Each result is read back by cfitsio (the
// cfref harness), not only by this package.

func TestEditHeaderInPlace(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "edit.fits")
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	data := make([]float32, 16)
	for i := range data {
		data[i] = float32(i) * 1.5
	}
	if _, err := WriteImage(f, nil, []int64{4, 4}, data); err != nil {
		t.Fatal(err)
	}
	f.Close()
	before := mustRead(t, path)

	fe, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	p, _ := fe.Primary()
	if err := p.Header().Set("OBJECT", "test object", "added by test"); err != nil {
		t.Fatal(err)
	}
	if err := fe.CloseAndFlush(); err != nil {
		t.Fatal(err)
	}
	after := mustRead(t, path)
	if len(after) != len(before) {
		t.Fatalf("in-place edit changed the file size %d -> %d", len(before), len(after))
	}
	if !bytes.Equal(after[2880:], before[2880:]) {
		t.Fatal("in-place header edit changed bytes outside the header")
	}

	cf := testutil.CfrefDump(t, path)
	c := cf.HDUs[0].Card("OBJECT")
	if c == nil || c.Str == nil || string(*c.Str) != "test object" || strings.TrimSpace(string(c.Comment)) != "added by test" {
		t.Fatalf("cfitsio OBJECT card = %+v", c)
	}
	for i, v := range data {
		if px := cf.HDUs[0].Image.Phys.Values(t)[i]; !px.Valid || px.V != float64(v) {
			t.Fatalf("cfitsio pixel %d = %+v, want %v", i, px, v)
		}
	}
}

// Close without Flush discards header changes (documented on OpenForEdit).
func TestCloseDoesNotFlush(t *testing.T) {
	path, orig := copySample(t, "iter_a.fit")
	f, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	p, _ := f.Primary()
	p.Header().Set("DISCARD", int64(1), "")
	f.Close()
	if !bytes.Equal(mustRead(t, path), orig) {
		t.Fatal("Close wrote header changes")
	}
}

// Opening every committed sample for edit, loading every header and
// flushing with no change must leave the file byte-identical. (Before,
// Flush re-encoded any loaded header whose encoding differed from the
// disk bytes, and the encoder drops CONTINUE cards: testprog.std lost the
// tail of its KEY_PKLS long string.)
func TestFlushWithoutChangesIsByteIdentical(t *testing.T) {
	var samples []string
	for _, pat := range []string{"*.fits", "*.fit", "*.std"} {
		m, _ := filepath.Glob(testutil.TestdataPath(t, pat))
		samples = append(samples, m...)
	}
	if len(samples) < 40 {
		t.Fatalf("only %d samples found in testdata", len(samples))
	}
	for _, s := range samples {
		name := filepath.Base(s)
		t.Run(name, func(t *testing.T) {
			path, orig := copySample(t, name)
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			for i := range f.NumHDU() {
				h, err := f.HDU(i)
				if err != nil {
					t.Fatal(err)
				}
				_ = h.Header()
			}
			if err := f.CloseAndFlush(); err != nil {
				t.Fatalf("Flush: %v", err)
			}
			if got := mustRead(t, path); !bytes.Equal(got, orig) {
				t.Fatalf("no-op Flush changed the file (first diff at byte %d)", firstDiff(got, orig))
			}
		})
	}
}

// lossyHeaderCases are header edits that header.Encode may not represent
// exactly. testprog.std's primary has a CONTINUE long string (KEY_PKLS);
// adding a card makes the writer re-encode that header. Each case states
// the outcome the writer must produce (outcome), and the test fails on any
// other:
//
//   - written: the edit is written and cfitsio reads back what was set
//     (check): an existing or newly set CONTINUE long string, and a float
//     with an integer value (written with a decimal point, type F).
//   - refusedBySet: Set/Add itself refuses the edit (NaN, a keyword named
//     END) and the file stays byte-identical.
//   - refusedByWriter: Flush / CopyHDUWithHeader refuses (commentary text
//     longer than one card changes the card count) and the file stays
//     byte-identical.
const (
	written = iota
	refusedBySet
	refusedByWriter
)

var lossyHeaderCases = []struct {
	name    string
	hdu     int
	outcome int
	mutate  func(h *header.Header) error
	check   func(t *testing.T, path string) // when the write succeeds
}{
	{"continue-survives", 0, written, func(h *header.Header) error {
		return h.Set("NEWKEY", int64(7), "forces re-encode")
	}, func(t *testing.T, path string) {
		want := testutil.LoadCfrefGolden(t, "testprog.std").HDUs[0].Card("KEY_PKLS")
		if want == nil || want.Str == nil {
			t.Fatal("golden has no KEY_PKLS long string")
		}
		c := testutil.CfrefDump(t, path).HDUs[0].Card("KEY_PKLS")
		if c == nil || c.Str == nil || string(*c.Str) != string(*want.Str) {
			t.Fatalf("cfitsio reads KEY_PKLS as %+v, want %q", c, *want.Str)
		}
	}},
	{"continue-new-longstr", 2, written, func(h *header.Header) error {
		return h.Set("LONGSTR", strings.Repeat("abcdefghij", 12), "a 120-character string")
	}, func(t *testing.T, path string) {
		c := testutil.CfrefDump(t, path).HDUs[2].Card("LONGSTR")
		if c == nil || c.Str == nil || string(*c.Str) != strings.Repeat("abcdefghij", 12) {
			t.Fatalf("cfitsio reads LONGSTR as %+v, want the 120-character string", c)
		}
	}},
	{"long-history", 2, refusedByWriter, func(h *header.Header) error {
		return h.Add("HISTORY", nil, strings.Repeat("h", 100))
	}, func(t *testing.T, path string) {
		var text string
		for _, c := range testutil.CfrefDump(t, path).HDUs[2].Cards {
			if c.Key == "HISTORY" {
				text += strings.TrimPrefix(string(c.Raw), "HISTORY ")
			}
		}
		if strings.Count(text, "h") != 100 {
			t.Fatalf("HISTORY text lost characters: %q", text)
		}
	}},
	{"nan", 2, refusedBySet, func(h *header.Header) error {
		return h.Set("NANKEY", math.NaN(), "")
	}, func(t *testing.T, path string) {
		if v := testutil.Fitsverify(t, path); v.Errors != 0 {
			t.Fatalf("wrote invalid FITS:\n%s", v.Report)
		}
	}},
	{"integer-valued-float", 2, written, func(h *header.Header) error {
		return h.Set("FLTKEY", float64(2), "a float")
	}, func(t *testing.T, path string) {
		c := testutil.CfrefDump(t, path).HDUs[2].Card("FLTKEY")
		if c == nil || c.Type != "F" || c.Float == nil || float64(*c.Float) != 2 {
			t.Fatalf("cfitsio reads FLTKEY as %+v, want a float 2 (type F)", c)
		}
	}},
	{"end-keyword", 2, refusedBySet, func(h *header.Header) error {
		return h.Set("END", int64(1), "")
	}, func(t *testing.T, path string) {
		if c := testutil.CfrefDump(t, path).HDUs[2].Card("END"); c == nil {
			t.Fatal("the END keyword card was lost")
		}
	}},
}

// Flush must never write a header that reads back differently from what
// was set (see lossyHeaderCases).
func TestFlushRefusesLossyEncoding(t *testing.T) {
	for _, tc := range lossyHeaderCases {
		t.Run(tc.name, func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			if merr := tc.mutate(hduHeader(t, f, tc.hdu)); merr != nil {
				t.Logf("header refused the edit: %v", merr)
				if tc.outcome != refusedBySet {
					t.Errorf("Set/Add refused the edit (%v); outcome %d expected", merr, tc.outcome)
				}
				if err := f.Flush(); err != nil {
					t.Fatalf("Flush after a refused edit: %v", err)
				}
				f.Close()
				if !bytes.Equal(mustRead(t, path), orig) {
					t.Fatal("a refused header edit changed the file")
				}
				return
			}
			if tc.outcome == refusedBySet {
				t.Errorf("Set/Add accepted the edit; it must refuse it")
			}
			ferr := f.Flush()
			f.Close()
			if (ferr != nil) != (tc.outcome == refusedByWriter) {
				t.Errorf("Flush err = %v, outcome %d expected", ferr, tc.outcome)
			}
			if ferr != nil {
				t.Logf("Flush refused: %v", ferr)
				if !bytes.Equal(mustRead(t, path), orig) {
					t.Fatal("Flush failed but changed the file")
				}
				return
			}
			tc.check(t, path)
		})
	}
}

// Header edits that change the data layout or the HDU type cannot be
// applied in place: Flush refuses with the reason (every structural
// keyword the OpenForEdit doc lists) and the file is unchanged.
func TestFlushRefusesStructuralChange(t *testing.T) {
	bump := func(key string, by int64) func(*testing.T, *header.Header) {
		return func(t *testing.T, h *header.Header) {
			v, err := h.Int(key)
			if err != nil {
				t.Fatalf("%s: %v", key, err)
			}
			if err := h.Set(key, v+by, ""); err != nil {
				t.Fatal(err)
			}
		}
	}
	set := func(key string, v any) func(*testing.T, *header.Header) {
		return func(t *testing.T, h *header.Header) {
			if err := h.Set(key, v, ""); err != nil {
				t.Fatal(err)
			}
		}
	}
	cases := []struct {
		name, sample string
		hdu          int
		mutate       func(*testing.T, *header.Header)
		want         string
	}{
		{"NAXIS1", "testprog.std", 1, bump("NAXIS1", 8), "NAXIS/NAXISn changed"},
		{"NAXIS", "iter_image.fit", 0, bump("NAXIS", -1), "NAXIS/NAXISn changed"},
		{"BITPIX", "testprog.std", 1, bump("BITPIX", 8), "BITPIX changed"},
		{"PCOUNT", "testprog.std", 1, bump("PCOUNT", 8), "PCOUNT/GCOUNT changed"},
		{"GCOUNT", "testprog.std", 1, bump("GCOUNT", 8), "PCOUNT/GCOUNT changed"},
		{"TFIELDS", "testprog.std", 1, bump("TFIELDS", 8), "TFIELDS changed"},
		{"SIMPLE", "iter_image.fit", 0, set("SIMPLE", false), "SIMPLE/XTENSION changed"},
		{"XTENSION", "testprog.std", 1, set("XTENSION", "A3DTABLE"), "SIMPLE/XTENSION changed"},
		{"ZIMAGE", "comp_rice_i16.fits", 1, set("ZIMAGE", false), "ZIMAGE changed"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			path, orig := copySample(t, c.sample)
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			c.mutate(t, hduHeader(t, f, c.hdu))
			err = f.Flush()
			f.Close()
			if err == nil {
				t.Fatalf("Flush accepted a change to %s", c.name)
			}
			if !strings.Contains(err.Error(), c.want) {
				t.Fatalf("Flush error = %q, want it to say %q", err, c.want)
			}
			if !bytes.Equal(mustRead(t, path), orig) {
				t.Fatal("file changed")
			}
		})
	}
}

// After a resizing Flush (which replaces the file) the same *File keeps
// working: a second in-place edit and an appended image land in the new
// file, and cfitsio reads all of it.
func TestFlushResizeThenKeepEditing(t *testing.T) {
	path, _ := copySample(t, "testprog.std")
	f, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	growHDU1(t, f)
	if err := f.Flush(); err != nil {
		t.Fatal(err)
	}
	if err := hduHeader(t, f, 6).Set("AFTER", "second edit", "after resize"); err != nil {
		t.Fatal(err)
	}
	pix := []int32{7, -8, 9}
	app, err := AppendImage(f, nil, []int64{3}, pix)
	if err != nil {
		t.Fatal(err)
	}
	vla, err := AppendBinaryTable(f, nil, []ColumnData{{Name: "V", DataVarInt32: [][]int32{{1, 2}, {3}}}})
	if err != nil {
		t.Fatal(err)
	}
	if err := f.Flush(); err != nil {
		t.Fatal(err)
	}
	// Edit the headers of the HDUs appended in this session; the image
	// header grows past one block (another resize).
	for i := 1; i <= 40; i++ {
		if err := app.Header().Set(fmt.Sprintf("APP%05d", i), int64(i), ""); err != nil {
			t.Fatal(err)
		}
	}
	if err := vla.Header().Set("EXTNAME", "VLA", ""); err != nil {
		t.Fatal(err)
	}
	if err := f.CloseAndFlush(); err != nil {
		t.Fatal(err)
	}

	want := testutil.LoadCfrefGolden(t, "testprog.std")
	got := testutil.CfrefDump(t, path)
	if got.NHDU != want.NHDU+2 {
		t.Fatalf("cfitsio reads %d HDUs, want %d", got.NHDU, want.NHDU+2)
	}
	trimmed := *got
	trimmed.HDUs = got.HDUs[:want.NHDU]
	trimmed.NHDU = want.NHDU
	compareEditedDump(t, want, &trimmed, shiftScenario{
		expect: map[int]func([]testutil.CfrefCard) []wantCard{
			1: expectGrowHDU1,
			6: func(cards []testutil.CfrefCard) []wantCard {
				return append(unchanged(cards), wantCard{key: "AFTER", value: "'second edit'", comment: "after resize"})
			},
		},
		blocks: map[int]int64{1: 3},
	})
	img := got.HDUs[want.NHDU].Image
	if img == nil || img.Bitpix != 32 {
		t.Fatalf("appended HDU: %+v", img)
	}
	for i, v := range pix {
		if !img.I64.Values(t)[i].Valid || img.I64.Values(t)[i].V != int64(v) {
			t.Errorf("appended pixel %d = %+v, want %d", i, img.I64.Values(t)[i], v)
		}
	}
	ah := got.HDUs[want.NHDU]
	if ah.DataStart-ah.HeaderStart != 2*2880 || ah.Card("APP00040") == nil || ah.Card("APP00040").Int == nil || *ah.Card("APP00040").Int != 40 {
		t.Errorf("appended image header not grown/edited as set (header %d bytes)", ah.DataStart-ah.HeaderStart)
	}
	vh := got.HDUs[want.NHDU+1]
	if c := vh.Card("EXTNAME"); c == nil || c.Str == nil || string(*c.Str) != "VLA" {
		t.Errorf("appended table EXTNAME = %+v", c)
	}
	if col := vh.Table.Column("V"); col == nil || len(col.RowCells(t)) != 2 || len(col.RowCells(t)[0].I64) != 2 || col.RowCells(t)[0].I64[1].V != 2 || col.RowCells(t)[1].I64[0].V != 3 {
		t.Errorf("appended VLA column as read by cfitsio: %+v", col)
	}
	if v := testutil.Fitsverify(t, path); v.Errors != 0 {
		t.Errorf("fitsverify: %d errors\n%s", v.Errors, v.Report)
	}
}

func TestOverwritePixels(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "ow.fits")
	f, _ := Create(path)
	if _, err := WriteImage(f, nil, []int64{8, 2}, make([]int16, 16)); err != nil {
		t.Fatal(err)
	}
	f.Close()

	fe, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	p, _ := fe.Primary()
	newData := []int16{100, 101, 102, 103, 104, 105, 106, 107,
		-200, -201, -202, -203, -204, -205, -206, 32767}
	if err := OverwritePixels(p, newData); err != nil {
		t.Fatal(err)
	}
	fe.Close()

	img := testutil.CfrefDump(t, path).HDUs[0].Image
	for i, v := range newData {
		if !img.I64.Values(t)[i].Valid || img.I64.Values(t)[i].V != int64(v) {
			t.Fatalf("cfitsio pixel %d = %+v, want %d", i, img.I64.Values(t)[i], v)
		}
	}
}

func TestOverwritePixelsShapeMismatch(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "mis.fits")
	f, _ := Create(path)
	WriteImage(f, nil, []int64{4, 4}, make([]int32, 16))
	f.Close()
	orig := mustRead(t, path)

	fe, _ := OpenForEdit(path)
	p, _ := fe.Primary()
	if err := OverwritePixels(p, []int32{1, 2, 3}); !errors.Is(err, ErrShapeMismatch) {
		t.Fatalf("short data: err %v, want ErrShapeMismatch", err)
	}
	if err := OverwritePixels(p, make([]int32, 17)); !errors.Is(err, ErrShapeMismatch) {
		t.Fatalf("long data: err %v, want ErrShapeMismatch", err)
	}
	fe.Close()
	if !bytes.Equal(mustRead(t, path), orig) {
		t.Fatal("rejected OverwritePixels changed the file")
	}
}

// AppendImage and AppendBinaryTable on a file opened for edit add HDUs at
// EOF; cfitsio reads the original HDUs unchanged and the new ones with the
// written values.
func TestAppendInEdit(t *testing.T) {
	path, orig := copySample(t, "iter_a.fit")
	fe, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	data := []float32{1, 2.5, -3, 4}
	if _, err := AppendImage(fe, nil, []int64{2, 2}, data); err != nil {
		t.Fatal(err)
	}
	if _, err := AppendBinaryTable(fe, nil, []ColumnData{
		{Name: "J", DataInt32: []int32{1, -2, 3}},
		{Name: "D", DataFloat64: []float64{0.5, 1e300, -2}},
		{Name: "S", DataString: []string{"abc", "def", "ghi"}},
	}); err != nil {
		t.Fatal(err)
	}
	fe.Close()

	got := mustRead(t, path)
	if !bytes.Equal(got[:len(orig)], orig) {
		t.Fatal("append modified the existing HDUs")
	}
	cf := testutil.CfrefDump(t, path)
	if cf.NHDU != 4 {
		t.Fatalf("cfitsio reads %d HDUs, want 4", cf.NHDU)
	}
	img := cf.HDUs[2].Image
	for i, v := range data {
		if !img.Phys.Values(t)[i].Valid || img.Phys.Values(t)[i].V != float64(v) {
			t.Errorf("image pixel %d = %+v, want %v", i, img.Phys.Values(t)[i], v)
		}
	}
	tb := cf.HDUs[3].Table
	if tb == nil || tb.NRows != 3 {
		t.Fatalf("table: %+v", tb)
	}
	j, d, s := tb.Column("J"), tb.Column("D"), tb.Column("S")
	if j == nil || d == nil || s == nil {
		t.Fatal("columns J, D, S not found by cfitsio")
	}
	for r := 0; r < 3; r++ {
		if v := j.RowCells(t)[r].I64[0]; !v.Valid || v.V != int64([]int32{1, -2, 3}[r]) {
			t.Errorf("J[%d] = %+v", r, v)
		}
		if v := d.RowCells(t)[r].Phys[0]; !v.Valid || v.V != []float64{0.5, 1e300, -2}[r] {
			t.Errorf("D[%d] = %+v", r, v)
		}
		if v := s.RowCells(t)[r].Str[0]; !v.Valid || v.V != []string{"abc", "def", "ghi"}[r] {
			t.Errorf("S[%d] = %+v", r, v)
		}
	}
	if v := testutil.Fitsverify(t, path); v.Errors != 0 {
		t.Errorf("fitsverify: %d errors\n%s", v.Errors, v.Report)
	}
	_ = os.Remove
}

// TestCloseAndFlushLeavesFileOpenOnFlushError pins the CloseAndFlush doc:
// when Flush fails, its error is returned and the *File stays open, so the
// caller can fix the header and flush again (or Close). The failure is a
// structural keyword change (NAXIS1), which Flush refuses; the header is
// then put back, a new card added, and the second CloseAndFlush must
// succeed and leave a file cfitsio reads with that card and the original
// pixels.
func TestCloseAndFlushLeavesFileOpenOnFlushError(t *testing.T) {
	path := filepath.Join(t.TempDir(), "cf.fits")
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	data := []int16{1, 2, 3, 4, 5, 6}
	if _, err := WriteImage(f, nil, []int64{3, 2}, data); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}
	before := mustRead(t, path)

	fe, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	p, _ := fe.Primary()
	if err := p.Header().Set("NAXIS1", int64(2), ""); err != nil {
		t.Fatal(err)
	}
	if err := fe.CloseAndFlush(); err == nil || !strings.Contains(err.Error(), "NAXIS") {
		t.Fatalf("CloseAndFlush with a changed NAXIS1: err %v, want the Flush refusal", err)
	}
	if !bytes.Equal(mustRead(t, path), before) {
		t.Fatal("refused Flush changed the file")
	}
	// Still open: the HDU and its header are usable and a retry works.
	p, err = fe.Primary()
	if err != nil {
		t.Fatalf("Primary after the failed CloseAndFlush: %v (file was closed)", err)
	}
	if err := p.Header().Set("NAXIS1", int64(3), ""); err != nil {
		t.Fatal(err)
	}
	if err := p.Header().Set("OBSERVER", "retry", ""); err != nil {
		t.Fatal(err)
	}
	if err := fe.CloseAndFlush(); err != nil {
		t.Fatalf("second CloseAndFlush: %v", err)
	}
	cf := testutil.CfrefDump(t, path)
	c := cf.HDUs[0].Card("OBSERVER")
	if c == nil || c.Str == nil || string(*c.Str) != "retry" {
		t.Fatalf("cfitsio OBSERVER card = %+v", c)
	}
	if im := cf.HDUs[0].Image; im == nil || im.Phys.Len() != len(data) {
		t.Fatalf("cfitsio image %+v", im)
	} else {
		for i, v := range im.Phys.Values(t) {
			if !v.Valid || v.V != float64(data[i]) {
				t.Fatalf("cfitsio pixel %d = %+v, want %d", i, v, data[i])
			}
		}
	}
}
