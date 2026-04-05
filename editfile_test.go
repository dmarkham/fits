package fits

import (
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/header"
)

func TestEditFileDropHDU(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "multi.fits")

	// Create a file with primary + 2 extensions (EXTNAME = KEEP, JUNK).
	f, _ := Create(path)
	WriteImage(f, nil, []int64{}, []float32{}) // primary
	h1 := header.New()
	h1.Set("EXTNAME", "KEEP", "")
	AppendImage(f, h1, []int64{4}, []float32{1, 2, 3, 4})
	h2 := header.New()
	h2.Set("EXTNAME", "JUNK", "")
	AppendImage(f, h2, []int64{2}, []float32{9, 9})
	f.Close()

	// Drop the JUNK HDU through EditFile.
	err := EditFile(path, func(in *File, out *Writer) error {
		for i := range in.NumHDU() {
			h, err := in.HDU(i)
			if err != nil {
				return err
			}
			if name, _ := h.Header().String("EXTNAME"); name == "JUNK" {
				if err := out.SkipHDU(h); err != nil {
					return err
				}
				continue
			}
			if err := out.CopyHDU(h); err != nil {
				return err
			}
		}
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}

	// Reopen and verify.
	fr, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fr.Close()
	if fr.NumHDU() != 2 {
		t.Fatalf("NumHDU after drop = %d", fr.NumHDU())
	}
	h, _ := fr.HDU(1)
	if name, _ := h.Header().String("EXTNAME"); name != "KEEP" {
		t.Fatalf("HDU 1 EXTNAME = %q", name)
	}
}

func TestEditFileFailureLeavesSourceIntact(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "intact.fits")
	f, _ := Create(path)
	WriteImage(f, nil, []int64{4}, []float32{1, 2, 3, 4})
	f.Close()

	// Read original bytes.
	orig, err := readFileContents(path)
	if err != nil {
		t.Fatal(err)
	}

	// Callback returns an error.
	EditFile(path, func(in *File, out *Writer) error {
		out.CopyHDU(&ImageHDU{rec: &hduRecord{}}) // bogus, expect error
		return errTestForce
	})

	// Source untouched.
	after, _ := readFileContents(path)
	if len(after) != len(orig) {
		t.Fatalf("source size changed: %d → %d", len(orig), len(after))
	}
}

var errTestForce = &testError{"forced test failure"}

type testError struct{ msg string }

func (e *testError) Error() string { return e.msg }

// readFileContents is a minimal local helper to avoid pulling in os here.
func readFileContents(path string) ([]byte, error) {
	f, err := Open(path)
	if err != nil {
		return nil, err
	}
	// We really just want raw bytes, but Open parses. Go through the OS
	// directly via the lower level.
	defer f.Close()
	// Collect the bytes by reading the whole file via the block reader.
	totalBlocks := f.br.BlockCount()
	out := make([]byte, 0, totalBlocks*2880)
	for i := int64(0); i < totalBlocks; i++ {
		b, err := f.br.ReadBlock(i)
		if err != nil {
			return nil, err
		}
		out = append(out, b[:]...)
	}
	return out, nil
}
