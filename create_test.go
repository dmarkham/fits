package fits

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Create truncates an existing file, and each WriteImage call leaves a
// complete, block-padded HDU on disk before Close: cfitsio reads the file
// while it is still open for writing.
func TestCreateTruncatesAndWritesCompleteHDUs(t *testing.T) {
	path := filepath.Join(t.TempDir(), "c.fits")
	if err := os.WriteFile(path, bytes.Repeat([]byte("J"), 5*2880+17), 0o644); err != nil {
		t.Fatal(err)
	}
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	pix := []int16{3, -1, 4, -1, 5, -9}
	if _, err := WriteImage(f, nil, []int64{3, 2}, pix); err != nil {
		t.Fatal(err)
	}
	st, _ := os.Stat(path)
	if st.Size()%2880 != 0 || st.Size() != 2*2880 {
		t.Fatalf("file size after one HDU = %d, want 5760 (old contents truncated, data padded)", st.Size())
	}
	img := testutil.CfrefDump(t, path).HDUs[0].Image
	for i, v := range pix {
		if !img.I64.Values(t)[i].Valid || img.I64.Values(t)[i].V != int64(v) {
			t.Fatalf("cfitsio pixel %d = %+v, want %d", i, img.I64.Values(t)[i], v)
		}
	}
}

// CreateWriter on an empty stream writes a file cfitsio reads. It does not
// truncate and it appends at the stream's end (documented), so a non-empty
// stream keeps its old bytes in front.
func TestCreateWriter(t *testing.T) {
	path := filepath.Join(t.TempDir(), "w.fits")
	fh, err := os.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	f, err := CreateWriter(fh)
	if err != nil {
		t.Fatal(err)
	}
	pix := []float64{0.25, -1e10, 3}
	if _, err := WriteImage(f, nil, []int64{3}, pix); err != nil {
		t.Fatal(err)
	}
	f.Close() // must not close fh
	if _, err := fh.Write(nil); err != nil {
		t.Fatalf("File.Close closed the caller's stream: %v", err)
	}
	fh.Close()
	img := testutil.CfrefDump(t, path).HDUs[0].Image
	for i, v := range pix {
		if !img.Phys.Values(t)[i].Valid || img.Phys.Values(t)[i].V != v {
			t.Fatalf("cfitsio pixel %d = %+v, want %v", i, img.Phys.Values(t)[i], v)
		}
	}

	// Non-empty stream: HDUs written after the existing bytes would make
	// invalid FITS (rule 7), so CreateWriter refuses and writes nothing.
	// The stream's position is not the test: one left at 0 or at the end
	// is refused alike.
	before := mustRead(t, path)
	for _, whence := range []int{io.SeekStart, io.SeekEnd} {
		fh, err = os.OpenFile(path, os.O_RDWR, 0)
		if err != nil {
			t.Fatal(err)
		}
		if _, err := fh.Seek(0, whence); err != nil {
			t.Fatal(err)
		}
		if f, err := CreateWriter(fh); err == nil {
			WriteImage(f, nil, []int64{3}, pix)
			fh.Close()
			t.Fatalf("CreateWriter accepted a %d-byte stream (now %d bytes)", len(before), len(mustRead(t, path)))
		}
		fh.Close()
		if !bytes.Equal(mustRead(t, path), before) {
			t.Fatal("refused CreateWriter changed the stream")
		}
	}
}

// TestCreateExtensionFirstGetsPrimary: on a freshly Created file, the
// writers whose HDU is an extension (AppendBinaryTable, AppendASCIITable,
// AppendCompressedImage, AppendCompressedFloat32Image) used to write that
// XTENSION as the first HDU with err = nil: fitsverify failed the file
// (8-14 errors) and Open returned "not a FITS file". They now write an
// empty primary first (BITPIX 16, NAXIS 0, as cfitsio's ffcrtb does for a
// table created in an empty file). Checked with fitsverify, the cfitsio
// dump of the file (2 HDUs, an empty IMAGE then the extension) and Open.
func TestCreateExtensionFirstGetsPrimary(t *testing.T) {
	appenders := []struct {
		name, cfType string
		run          func(f *File) error
	}{
		{"AppendBinaryTable", "BINARY_TBL", func(f *File) error {
			_, err := AppendBinaryTable(f, nil, []ColumnData{{Name: "A", DataInt32: []int32{1, 2, 3}}})
			return err
		}},
		{"AppendASCIITable", "ASCII_TBL", func(f *File) error {
			_, err := AppendASCIITable(f, nil, []ASCIIColumnData{{Name: "A", TForm: "I4", DataInt64: []int64{1, 2, 3}}})
			return err
		}},
		{"AppendCompressedImage", "IMAGE", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{4, 3}, []int16{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, CompressOptions{})
			return err
		}},
		{"AppendCompressedFloat32Image", "IMAGE", func(f *File) error {
			_, err := AppendCompressedFloat32Image(f, nil, []int64{4, 3}, make([]float32, 12), CompressFloatOptions{})
			return err
		}},
	}
	for _, a := range appenders {
		t.Run(a.name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "x.fits")
			f, err := Create(path)
			if err != nil {
				t.Fatal(err)
			}
			if err := a.run(f); err != nil {
				t.Fatal(err)
			}
			if f.NumHDU() != 2 {
				t.Errorf("NumHDU after the append = %d, want 2 (primary + extension)", f.NumHDU())
			}
			if err := f.Close(); err != nil {
				t.Fatal(err)
			}
			if v := testutil.Fitsverify(t, path); v.Errors != 0 || v.Warnings != 0 {
				t.Fatalf("fitsverify: %d errors, %d warnings\n%s", v.Errors, v.Warnings, v.Report)
			}
			g := testutil.CfrefDump(t, path)
			if g.NHDU != 2 || g.HDUs[0].Type != "IMAGE" || g.HDUs[0].Image == nil || g.HDUs[0].Image.Naxis != 0 ||
				g.HDUs[0].Image.Bitpix != 16 || g.HDUs[1].Type != a.cfType {
				t.Fatalf("cfitsio reads %d HDUs: %+v", g.NHDU, g.HDUs)
			}
			r, err := Open(path)
			if err != nil {
				t.Fatalf("Open: %v", err)
			}
			defer r.Close()
			if r.NumHDU() != 2 {
				t.Fatalf("Open: %d HDUs", r.NumHDU())
			}
		})
	}
}
