package fits

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"sync/atomic"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Tests for the documented error contract: every typed error and sentinel
// a doc promises is asserted by name (errors.Is / errors.As), and every
// refusal leaves the file untouched.

// cfitsioOpen is the `imgref open` golden for a file the library refuses.
type cfitsioOpen struct {
	Cfitsio    string `json:"cfitsio"`
	OpenStatus int    `json:"open_status"`
	NHDU       int    `json:"nhdu"`
	HDUs       []struct {
		Index      int `json:"index"`
		ReadStatus int `json:"read_status"`
	} `json:"hdus"`
	NextStatus int `json:"next_status"`
}

// TestRefusedFramingMatchesCfitsio: the byte-assembled files in
// testdata/images/refused (gen.py gen_refused) are refused by Open with
// the documented typed error, and cfitsio 4.7.0 (imgref open, golden
// refused__<name>.open.json) is shown to fail on the same HDU:
//
//   - a data unit whose declared bytes run past EOF: *ErrTruncatedData
//     with the HDU, the declared and the present byte counts. cfitsio
//     opens the file and fails the read (status 107/108); we fail at
//     Open, deliberately stricter (scan.go).
//   - a missing mandatory keyword: *ErrMissingRequiredKeyword with the
//     HDU and keyword. cfitsio refuses to open a primary without BITPIX
//     (222), NAXIS (223) or NAXISn (224), and will not move to an
//     extension without them (next_status 222/224; it reports one HDU
//     fewer instead of failing, where we fail loudly).
//   - a BINTABLE/TABLE extension without PCOUNT, GCOUNT or TFIELDS:
//     *ErrMissingRequiredKeyword. cfitsio will not move to it (228, 229,
//     230). An IMAGE extension may omit PCOUNT/GCOUNT (cfitsio accepts
//     it, so do we).
//   - XTENSION with a non-string value: *ErrMissingRequiredKeyword
//     XTENSION. cfitsio reads that HDU as an IMAGE; we refuse it
//     (deliberate difference, docs/VALIDATION.md).
func TestRefusedFramingMatchesCfitsio(t *testing.T) {
	type want struct {
		trunc   *ErrTruncatedData
		missing *ErrMissingRequiredKeyword
		check   func(t *testing.T, g cfitsioOpen)
	}
	readFails := func(hdu int) func(*testing.T, cfitsioOpen) {
		return func(t *testing.T, g cfitsioOpen) {
			if g.OpenStatus != 0 || hdu >= len(g.HDUs) || g.HDUs[hdu].ReadStatus == 0 {
				t.Errorf("cfitsio golden %+v: want the file opened and HDU %d's read failing", g, hdu)
			}
		}
	}
	openFails := func(status int) func(*testing.T, cfitsioOpen) {
		return func(t *testing.T, g cfitsioOpen) {
			if g.OpenStatus != status {
				t.Errorf("cfitsio open status %d, want %d", g.OpenStatus, status)
			}
		}
	}
	extRefused := func(status int) func(*testing.T, cfitsioOpen) {
		return func(t *testing.T, g cfitsioOpen) {
			if g.OpenStatus != 0 || g.NHDU != 1 || g.NextStatus != status {
				t.Errorf("cfitsio golden %+v: want 1 HDU and status %d moving to HDU 2", g, status)
			}
		}
	}
	cases := map[string]want{
		"truncated_data.fits":          {trunc: &ErrTruncatedData{HDU: 0, Expected: 4000, Got: 2880}, check: readFails(0)},
		"truncated_data_ext.fits":      {trunc: &ErrTruncatedData{HDU: 1, Expected: 6, Got: 2}, check: readFails(1)},
		"missing_naxis2.fits":          {missing: &ErrMissingRequiredKeyword{HDU: 0, Keyword: "NAXIS2"}, check: openFails(224)},
		"missing_bitpix.fits":          {missing: &ErrMissingRequiredKeyword{HDU: 0, Keyword: "BITPIX"}, check: openFails(222)},
		"missing_naxis.fits":           {missing: &ErrMissingRequiredKeyword{HDU: 0, Keyword: "NAXIS"}, check: openFails(223)},
		"ext_missing_naxis1.fits":      {missing: &ErrMissingRequiredKeyword{HDU: 1, Keyword: "NAXIS1"}, check: extRefused(224)},
		"ext_missing_bitpix.fits":      {missing: &ErrMissingRequiredKeyword{HDU: 1, Keyword: "BITPIX"}, check: extRefused(222)},
		"ext_xtension_not_string.fits": {missing: &ErrMissingRequiredKeyword{HDU: 1, Keyword: "XTENSION"}},
		// Table extensions: cfitsio will not move to them without PCOUNT
		// (228), GCOUNT (229) or TFIELDS (230).
		"tbl_no_pcount.fits":    {missing: &ErrMissingRequiredKeyword{HDU: 1, Keyword: "PCOUNT"}, check: extRefused(228)},
		"tbl_no_gcount.fits":    {missing: &ErrMissingRequiredKeyword{HDU: 1, Keyword: "GCOUNT"}, check: extRefused(229)},
		"tbl_no_tfields.fits":   {missing: &ErrMissingRequiredKeyword{HDU: 1, Keyword: "TFIELDS"}, check: extRefused(230)},
		"ascii_no_tfields.fits": {missing: &ErrMissingRequiredKeyword{HDU: 1, Keyword: "TFIELDS"}, check: extRefused(230)},
	}
	entries, err := os.ReadDir(imagesPath(t, "refused"))
	if err != nil {
		t.Fatalf("committed samples missing (failure, not a skip): %v", err)
	}
	if len(entries) != len(cases) {
		t.Fatalf("testdata/images/refused has %d files, the test names %d", len(entries), len(cases))
	}
	for name, w := range cases {
		t.Run(name, func(t *testing.T) {
			var g cfitsioOpen
			loadJSON(t, imagesPath(t, "golden", "refused__"+name+".open.json"), &g)
			if g.Cfitsio != "4.7.0" {
				t.Fatalf("golden made with cfitsio %s", g.Cfitsio)
			}
			if w.check != nil {
				w.check(t, g)
			} else if g.OpenStatus != 0 || g.NHDU != 2 {
				t.Errorf("cfitsio golden %+v: the documented difference is that cfitsio reads this HDU", g)
			}
			path := imagesPath(t, "refused", name)
			for _, op := range openPaths {
				_, err := op.open(t, path)
				switch {
				case w.trunc != nil:
					var e *ErrTruncatedData
					if !errors.As(err, &e) || *e != *w.trunc {
						t.Errorf("%s: err %v (%#v), want %+v", op.name, err, e, *w.trunc)
					}
				case w.missing != nil:
					var e *ErrMissingRequiredKeyword
					if !errors.As(err, &e) || *e != *w.missing {
						t.Errorf("%s: err %v, want %+v", op.name, err, *w.missing)
					}
				}
			}
		})
	}
}

// topSamples are the top-level committed samples (cfitsio's own test
// files and the compression corpus).
func topSamples(t *testing.T) []string {
	t.Helper()
	out := allSamples(t)
	if len(out) < 100 {
		t.Fatalf("only %d top-level samples", len(out))
	}
	return out
}

// copySampleBytes copies a committed sample into a fresh temp dir.
func copySampleBytes(t *testing.T, src string) (string, []byte) {
	t.Helper()
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatalf("committed sample unreadable: %v", err)
	}
	p := filepath.Join(t.TempDir(), filepath.Base(src))
	if err := os.WriteFile(p, b, 0o644); err != nil {
		t.Fatal(err)
	}
	return p, b
}

// TestReadOnlyRefusesEveryWrite: a file opened read-only (Open, OpenReader,
// OpenReadAll) refuses every write with ErrReadOnly (Flush after a header
// edit, OverwritePixels, WriteImage, AppendImage, AppendBinaryTable,
// AppendCompressedImage, WriteChecksum), for every top-level sample, and
// the file on disk is byte-identical afterwards. OverwritePixels also
// refuses ModeCreate (it needs ModeEdit).
func TestReadOnlyRefusesEveryWrite(t *testing.T) {
	for _, src := range topSamples(t) {
		t.Run(filepath.Base(src), func(t *testing.T) {
			path, orig := copySampleBytes(t, src)
			for _, op := range openPaths {
				f, err := op.open(t, path)
				if err != nil {
					t.Fatalf("%s: %v", op.name, err)
				}
				for i := 0; i < f.NumHDU(); i++ {
					h, err := f.HDU(i)
					if err != nil {
						continue // random groups / undecodable: nothing to write through
					}
					if hdr, err := HeaderOf(h); err == nil {
						if err := hdr.Set("XRDONLY", int64(1), "edit that must not persist"); err != nil {
							t.Fatal(err)
						}
					}
					if img, ok := h.(*ImageHDU); ok {
						if err := OverwritePixels(img, make([]float64, img.pixelCount())); !errors.Is(err, ErrReadOnly) {
							t.Errorf("%s HDU %d: OverwritePixels err %v, want ErrReadOnly", op.name, i, err)
						}
					}
					if _, _, err := WriteChecksum(h); !errors.Is(err, ErrReadOnly) {
						t.Errorf("%s HDU %d: WriteChecksum err %v, want ErrReadOnly", op.name, i, err)
					}
				}
				if err := f.Flush(); !errors.Is(err, ErrReadOnly) {
					t.Errorf("%s: Flush after a header edit: err %v, want ErrReadOnly", op.name, err)
				}
				if _, err := WriteImage(f, nil, []int64{2}, []int16{1, 2}); !errors.Is(err, ErrReadOnly) {
					t.Errorf("%s: WriteImage err %v, want ErrReadOnly", op.name, err)
				}
				if _, err := AppendImage(f, nil, []int64{2}, []int16{1, 2}); !errors.Is(err, ErrReadOnly) {
					t.Errorf("%s: AppendImage err %v, want ErrReadOnly", op.name, err)
				}
				if _, err := AppendBinaryTable(f, nil, []ColumnData{{Name: "A", DataInt32: []int32{1}}}); !errors.Is(err, ErrReadOnly) {
					t.Errorf("%s: AppendBinaryTable err %v, want ErrReadOnly", op.name, err)
				}
				if _, err := AppendCompressedImage(f, nil, []int64{2}, []int16{1, 2}, CompressOptions{}); !errors.Is(err, ErrReadOnly) {
					t.Errorf("%s: AppendCompressedImage err %v, want ErrReadOnly", op.name, err)
				}
				f.Close()
				if b, _ := os.ReadFile(path); !bytes.Equal(b, orig) {
					t.Fatalf("%s: the refused writes changed the file", op.name)
				}
			}
		})
	}
	// ModeCreate: OverwritePixels needs ModeEdit.
	path := filepath.Join(t.TempDir(), "create.fits")
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	img, err := WriteImage(f, nil, []int64{3}, []int16{1, 2, 3})
	if err != nil {
		t.Fatal(err)
	}
	before, _ := os.ReadFile(path)
	if err := OverwritePixels(img, []int16{4, 5, 6}); !errors.Is(err, ErrReadOnly) {
		t.Errorf("OverwritePixels in ModeCreate: err %v, want ErrReadOnly", err)
	}
	f.Close()
	if after, _ := os.ReadFile(path); !bytes.Equal(before, after) {
		t.Error("refused OverwritePixels changed the file")
	}
}

// overwriteTarget runs OverwritePixels[T] with n zero values.
type overwriteTarget struct {
	name string
	bp   int // the BITPIX WriteImage uses for T
	run  func(h *ImageHDU, n int64) error
}

func ow[T Numeric](name string, bp int) overwriteTarget {
	return overwriteTarget{name, bp, func(h *ImageHDU, n int64) error { return OverwritePixels(h, make([]T, n)) }}
}

var overwriteTargets = []overwriteTarget{
	ow[uint8]("uint8", 8), ow[int8]("int8", 8), ow[int16]("int16", 16), ow[uint16]("uint16", 16),
	ow[int32]("int32", 32), ow[uint32]("uint32", 32), ow[int64]("int64", 64), ow[uint64]("uint64", 64),
	ow[float32]("float32", -32), ow[float64]("float64", -64),
}

// TestOverwritePixelsTypeAndShapeMatrix: for every cfitsio-written image
// in testdata/images/cfwrite and every Go type, OverwritePixels in
// ModeEdit refuses a T that is not the type WriteImage would use for the
// HDU's BITPIX (a *TypeMismatchError wrapping ErrTypeMismatch, Requested
// and BITPIX set), refuses a short or long slice of an accepted T with
// ErrShapeMismatch, and leaves the file byte-identical every time.
func TestOverwritePixelsTypeAndShapeMatrix(t *testing.T) {
	files, _ := filepath.Glob(imagesPath(t, "cfwrite", "*.fits"))
	if len(files) < 30 {
		t.Fatalf("only %d cfwrite samples", len(files))
	}
	bitpixSeen := map[int]bool{}
	for _, src := range files {
		t.Run(filepath.Base(src), func(t *testing.T) {
			path, orig := copySampleBytes(t, src)
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			p, err := f.Primary()
			if err != nil {
				t.Fatal(err)
			}
			n := p.pixelCount()
			if n == 0 {
				return
			}
			bp := p.BITPIX()
			bitpixSeen[bp] = true
			accepted := 0
			for _, tg := range overwriteTargets {
				if tg.bp != bp {
					err := tg.run(p, n)
					var tm *TypeMismatchError
					if !errors.Is(err, ErrTypeMismatch) || !errors.As(err, &tm) || tm.Requested != tg.name || tm.BITPIX != bp {
						t.Errorf("OverwritePixels[%s] on BITPIX %d: err %v, want a TypeMismatchError{Requested %s, BITPIX %d}", tg.name, bp, err, tg.name, bp)
					}
					continue
				}
				accepted++
				for _, m := range []int64{n - 1, n + 1} {
					if err := tg.run(p, m); !errors.Is(err, ErrShapeMismatch) {
						t.Errorf("OverwritePixels[%s] with %d of %d values: err %v, want ErrShapeMismatch", tg.name, m, n, err)
					}
				}
			}
			if accepted == 0 {
				t.Errorf("no Go type accepted for BITPIX %d", bp)
			}
			if b, _ := os.ReadFile(path); !bytes.Equal(b, orig) {
				t.Fatal("a refused OverwritePixels changed the file")
			}
		})
	}
	for _, bp := range []int{8, 16, 32, 64, -32, -64} {
		if !bitpixSeen[bp] {
			t.Errorf("no cfwrite sample with BITPIX %d", bp)
		}
	}
}

// TestValidationBranchesRefuse: every input check of ReadSubset,
// WriteImage/AppendImage and File.HDU returns an error (never a panic,
// never padding read as data, never a partial write).
func TestValidationBranchesRefuse(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "v.fits")
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, nil, []int64{4, 2}, []int16{1, 2, 3, 4, 5, 6, 7, 8}); err != nil {
		t.Fatal(err)
	}
	before, _ := os.ReadFile(path)
	for _, c := range []struct {
		shape []int64
		n     int
	}{{[]int64{4, 2}, 7}, {[]int64{4, 2}, 9}, {[]int64{3}, 0}, {[]int64{}, 1}} {
		if _, err := AppendImage(f, nil, c.shape, make([]int16, c.n)); !errors.Is(err, ErrShapeMismatch) {
			t.Errorf("AppendImage shape %v with %d values: err %v, want ErrShapeMismatch", c.shape, c.n, err)
		}
		if _, err := WriteImage(f, nil, c.shape, make([]float32, c.n)); !errors.Is(err, ErrShapeMismatch) {
			t.Errorf("WriteImage shape %v with %d values: err %v, want ErrShapeMismatch", c.shape, c.n, err)
		}
	}
	for _, shape := range [][]int64{{0}, {-1}, {4, 0}} {
		if _, err := WriteImage(f, nil, shape, []int16{}); err == nil {
			t.Errorf("WriteImage accepted the non-positive shape %v", shape)
		}
	}
	f.Close()
	if after, _ := os.ReadFile(path); !bytes.Equal(before, after) {
		t.Fatal("refused WriteImage/AppendImage calls wrote to the file")
	}

	r, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	for _, i := range []int{-1, r.NumHDU(), r.NumHDU() + 5} {
		if h, err := r.HDU(i); err == nil || h != nil {
			t.Errorf("HDU(%d) of %d: (%v, %v), want an error", i, r.NumHDU(), h, err)
		}
	}
	if _, err := r.HDUByName("NOSUCH"); err == nil {
		t.Error("HDUByName of a missing EXTNAME succeeded")
	}
	p, err := r.Primary()
	if err != nil {
		t.Fatal(err)
	}
	if got, err := ReadSubset[int16](p, []int64{0, 0}, []int64{4, 2}, []int64{1, 1}); err != nil || fmt.Sprint(got) != "[1 2 3 4 5 6 7 8]" {
		t.Fatalf("the valid full subset: %v %v", got, err)
	}
	bad := []struct{ lo, up, st []int64 }{
		{[]int64{0, 0}, []int64{4, 3}, []int64{1, 1}},  // upper past the shape (would read padding)
		{[]int64{0, 0}, []int64{5, 2}, []int64{1, 1}},  // upper past axis 1
		{[]int64{-1, 0}, []int64{4, 2}, []int64{1, 1}}, // negative lower
		{[]int64{2, 0}, []int64{2, 2}, []int64{1, 1}},  // empty range
		{[]int64{3, 0}, []int64{2, 2}, []int64{1, 1}},  // lower > upper
		{[]int64{0, 0}, []int64{4, 2}, []int64{0, 1}},  // stride 0
		{[]int64{0, 0}, []int64{4, 2}, []int64{1, -2}}, // negative stride
		{[]int64{0}, []int64{4, 2}, []int64{1, 1}},     // wrong lengths
		{[]int64{0, 0}, []int64{4}, []int64{1, 1}},
		{[]int64{0, 0}, []int64{4, 2}, []int64{1}},
		{[]int64{0, 0, 0}, []int64{4, 2, 1}, []int64{1, 1, 1}},
	}
	for _, b := range bad {
		if got, err := ReadSubset[int16](p, b.lo, b.up, b.st); err == nil {
			t.Errorf("ReadSubset(%v, %v, %v) = %v, want an error", b.lo, b.up, b.st, got)
		}
	}
}

// cancelAfter is a context whose Err turns to context.Canceled after n
// calls, to cancel in the middle of an operation that polls it.
type cancelAfter struct {
	context.Context
	n atomic.Int64
}

func newCancelAfter(n int64) *cancelAfter {
	c := &cancelAfter{Context: context.Background()}
	c.n.Store(n)
	return c
}

func (c *cancelAfter) Err() error {
	if c.n.Add(-1) < 0 {
		return context.Canceled
	}
	return nil
}

// TestContextFunctions: each *Context function returns context.Canceled
// for a cancelled context without touching anything, and the opens also
// stop when the context is cancelled during the HDU scan (OpenContext,
// OpenForEditContext and EditFileContext consult it once per HDU).
func TestContextFunctions(t *testing.T) {
	done, cancel := context.WithCancel(context.Background())
	cancel()
	src := testutil.TestdataPath(t, "testprog.std") // 7 HDUs

	if _, err := OpenContext(done, src); !errors.Is(err, context.Canceled) {
		t.Errorf("OpenContext(cancelled): %v", err)
	}
	// Cancelled after the entry check and two HDUs of the scan.
	if _, err := OpenContext(newCancelAfter(3), src); !errors.Is(err, context.Canceled) {
		t.Errorf("OpenContext cancelled during the scan: %v", err)
	}
	if f, err := OpenContext(newCancelAfter(1000), src); err != nil {
		t.Errorf("OpenContext with a live context: %v", err)
	} else {
		f.Close()
	}

	path, orig := copySampleBytes(t, src)
	if _, err := OpenForEditContext(done, path); !errors.Is(err, context.Canceled) {
		t.Errorf("OpenForEditContext(cancelled): %v", err)
	}
	if _, err := OpenForEditContext(newCancelAfter(3), path); !errors.Is(err, context.Canceled) {
		t.Errorf("OpenForEditContext cancelled during the scan: %v", err)
	}

	// EditFileContext: nothing happens — no recovery pass (a stale temp
	// file of an interrupted edit stays), no new temp file, fn not called.
	stale := path + tempSuffix
	if err := os.WriteFile(stale, []byte("stale"), 0o644); err != nil {
		t.Fatal(err)
	}
	called := false
	fn := func(in *File, out *Writer) error { called = true; return nil }
	if err := EditFileContext(done, path, fn); !errors.Is(err, context.Canceled) {
		t.Errorf("EditFileContext(cancelled): %v", err)
	}
	if b, err := os.ReadFile(stale); err != nil || string(b) != "stale" {
		t.Errorf("EditFileContext(cancelled) touched the stale temp file: %q %v", b, err)
	}
	os.Remove(stale)
	if err := EditFileContext(newCancelAfter(3), path, fn); !errors.Is(err, context.Canceled) {
		t.Errorf("EditFileContext cancelled during the source scan: %v", err)
	}
	if called {
		t.Error("fn was called under a cancelled context")
	}
	if _, err := os.Stat(stale); !os.IsNotExist(err) {
		t.Errorf("EditFileContext left a temp file: %v", err)
	}
	if b, _ := os.ReadFile(path); !bytes.Equal(b, orig) {
		t.Error("a cancelled EditFileContext changed the source")
	}

	newPath := filepath.Join(t.TempDir(), "new.fits")
	if _, err := CreateContext(done, newPath); !errors.Is(err, context.Canceled) {
		t.Errorf("CreateContext(cancelled): %v", err)
	}
	if _, err := os.Stat(newPath); !os.IsNotExist(err) {
		t.Errorf("CreateContext(cancelled) created the file: %v", err)
	}

	f, err := Open(src)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	p, err := f.Primary()
	if err != nil {
		t.Fatal(err)
	}
	if v, err := ReadPixelsContext[float64](done, p); !errors.Is(err, context.Canceled) || v != nil {
		t.Errorf("ReadPixelsContext(cancelled): %v %v", len(v), err)
	}
	if _, err := ReadPixelsContext[float64](context.Background(), p); err != nil {
		t.Errorf("ReadPixelsContext: %v", err)
	}
}

// TestVarColumnViewsAgree: for every VLA column of the cfitsio VLA
// samples (whose At() values TestTablesMatchCfitsio pins to cfitsio),
// Rows() yields every row in order with the same zero-copy slice, Raw()
// returns the backing slab and offsets with values[off[r]:off[r+1]] equal
// to At(r), and early break from Rows() stops the iteration.
func TestVarColumnViewsAgree(t *testing.T) {
	checked := 0
	for _, name := range []string{"cf_vla_p.fits", "cf_vla_q.fits", "astropy_variable_length_table.fits"} {
		f, err := Open(testutil.TestdataPath(t, filepath.Join("tables", name)))
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		for i := 1; i < f.NumHDU(); i++ {
			h, err := f.HDU(i)
			if err != nil {
				t.Fatal(err)
			}
			bt, ok := h.(*BinaryTableHDU)
			if !ok {
				continue
			}
			cols, err := bt.Columns()
			if err != nil {
				t.Fatal(err)
			}
			for _, c := range cols {
				v, err := ReadVarColumn[float64](bt, c.Index)
				if err != nil {
					continue // not a numeric VLA column
				}
				checked++
				where := fmt.Sprintf("%s HDU %d %s", name, i, c.Name)
				vals, off := v.Raw()
				if len(off) != v.Len()+1 || off[0] != 0 || off[len(off)-1] != int64(len(vals)) {
					t.Errorf("%s: Raw offsets %v for %d rows and %d values", where, off, v.Len(), len(vals))
					continue
				}
				next := 0
				for r, row := range v.Rows() {
					if r != next {
						t.Errorf("%s: Rows yielded row %d, want %d", where, r, next)
					}
					next++
					at := v.At(r)
					if len(row) != len(at) || (len(at) > 0 && &row[0] != &at[0]) {
						t.Errorf("%s row %d: Rows slice is not At's zero-copy view", where, r)
					}
					raw := vals[off[r]:off[r+1]]
					if len(raw) != len(at) || (len(at) > 0 && &raw[0] != &at[0]) {
						t.Errorf("%s row %d: Raw slab does not back At", where, r)
					}
				}
				if next != v.Len() {
					t.Errorf("%s: Rows yielded %d rows, want %d", where, next, v.Len())
				}
				if v.Len() > 1 {
					n := 0
					for range v.Rows() {
						n++
						break
					}
					if n != 1 {
						t.Errorf("%s: break from Rows did not stop it", where)
					}
				}
			}
		}
		f.Close()
	}
	if checked < 10 {
		t.Fatalf("only %d VLA columns checked", checked)
	}
}

// TestColumnByNameIsCaseSensitive: on every table of the committed table
// samples, ColumnByName finds each column under its exact TTYPE as
// cfitsio reads it (golden), with that column's index, and does not find
// it under a different case (the documented case-sensitive rule; cfitsio's
// fits_get_colnum is case-insensitive by default, a deliberate
// difference).
func TestColumnByNameIsCaseSensitive(t *testing.T) {
	exact, folded := 0, 0
	for _, s := range tableSamples(t) {
		g := loadGoldenFile(t, s.golden)
		f, err := Open(s.path)
		if err != nil {
			t.Fatalf("%s: %v", s.name, err)
		}
		for _, gh := range g.HDUs {
			if gh.Type != "BINARY_TBL" && gh.Type != "ASCII_TBL" {
				continue
			}
			h, err := f.HDU(gh.Index)
			if err != nil {
				continue // e.g. tile-compressed
			}
			var byName func(string) (Column, bool)
			switch tb := h.(type) {
			case *BinaryTableHDU:
				byName = tb.ColumnByName
			case *ASCIITableHDU:
				byName = tb.ColumnByName
			default:
				continue
			}
			names := map[string]int{}
			for _, cd := range gh.Cards {
				k := string(cd.Key)
				if strings.HasPrefix(k, "TTYPE") && cd.Str != nil {
					var idx int
					fmt.Sscanf(k[5:], "%d", &idx)
					names[strings.TrimRight(string(*cd.Str), " ")] = idx
				}
			}
			for name, idx := range names {
				c, ok := byName(name)
				if !ok || c.Index != idx {
					// a duplicate TTYPE resolves to the first column
					if !ok {
						t.Errorf("%s HDU %d: ColumnByName(%q) not found", s.name, gh.Index, name)
					}
					continue
				}
				exact++
				for _, other := range []string{strings.ToLower(name), strings.ToUpper(name)} {
					if other == name {
						continue
					}
					if _, dup := names[other]; dup {
						continue
					}
					if _, ok := byName(other); ok {
						t.Errorf("%s HDU %d: ColumnByName(%q) found column %q (case-insensitive)", s.name, gh.Index, other, name)
					}
					folded++
				}
			}
			if _, ok := byName("NO_SUCH_COLUMN_X"); ok {
				t.Errorf("%s: ColumnByName of a missing name succeeded", s.name)
			}
		}
		f.Close()
	}
	if exact < 50 || folded < 50 {
		t.Fatalf("only %d exact and %d case-changed lookups checked", exact, folded)
	}
}

// TestA3DTABLEIsABinaryTable: XTENSION = 'A3DTABLE', the pre-standard name
// cfitsio (getkey.c) and astropy still read as a BINTABLE, opens as a
// *BinaryTableHDU. The sample tables/cf_bintable_a3dtable.fits is
// cfitsio-written and byte-patched (gen_tables.c); its column values are
// compared with cfitsio's reading by TestTablesMatchCfitsio and with
// astropy's by TestTablesMatchAstropy.
//
// The same holds for XTENSION = '3DTABLE' (tables/cf_bintable_3dtable.fits),
// the other legacy name cfitsio reads as BINARY_TBL; astropy 8.0.2 does
// not know it. The IUEIMAGE name cfitsio reads as IMAGE_HDU is pinned by
// images/framing/xtension_iueimage.fits (TestHDUAddressesMatchCfitsio).
func TestA3DTABLEIsABinaryTable(t *testing.T) {
	for _, c := range []struct{ file, xt string }{
		{"cf_bintable_a3dtable.fits", "XTENSION= 'A3DTABLE'"},
		{"cf_bintable_3dtable.fits", "XTENSION= '3DTABLE '"},
	} {
		t.Run(c.file, func(t *testing.T) {
			g := loadGoldenFile(t, testutil.TestdataPath(t, "tables/golden/"+c.file+".json"))
			if g.NHDU != 2 || g.HDUs[1].Type != "BINARY_TBL" {
				t.Fatalf("cfitsio golden: %d HDUs, HDU 1 %q", g.NHDU, g.HDUs[1].Type)
			}
			if raw := string(g.HDUs[1].Cards[0].Raw); !strings.HasPrefix(raw, c.xt) {
				t.Fatalf("sample's first extension card is %q", raw)
			}
			f, err := Open(testutil.TestdataPath(t, "tables/"+c.file))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			h, err := f.HDU(1)
			if err != nil {
				t.Fatal(err)
			}
			bt, ok := h.(*BinaryTableHDU)
			if !ok {
				t.Fatalf("HDU 1 is %T, want *BinaryTableHDU", h)
			}
			if v, err := ReadColumn[int32](bt, 1); err != nil || fmt.Sprint(v) != "[-7 0 2147483647]" {
				t.Errorf("J1 = %v %v", v, err)
			}
		})
	}
}
