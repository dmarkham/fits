package fits

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// After a Flush, the *File keeps describing the file on disk: CopyTo
// emits the new header bytes, an EditFile that copies HDUs out of the
// flushed *File copies the new headers, and HDUByName finds HDUs by their
// new EXTNAME. testprog.std HDU 1 ('Test-BINTABLE', a name HDUs 4 and 5
// share) is renamed RENAMED, on the in-place path and on the resize path
// (30 more cards: the header grows by a block).
//
// Truth: the bytes on disk after the Flush (CopyTo and the EditFile copy
// must equal them) and cfitsio's reading of the copy (EXTNAME RENAMED on
// HDU 1, 'Test-BINTABLE' still on HDUs 4 and 5).
func TestFlushedFileCopiesAndIndexesNewHeaders(t *testing.T) {
	for _, path := range []struct {
		name   string
		resize bool
	}{{"inplace", false}, {"resize", true}} {
		t.Run(path.name, func(t *testing.T) {
			p, orig := copySample(t, "testprog.std")
			f, err := OpenForEdit(p)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			if err := hduHeader(t, f, 1).Set("EXTNAME", "RENAMED", "renamed by test"); err != nil {
				t.Fatal(err)
			}
			if path.resize {
				growHDU1(t, f)
			}
			if err := f.Flush(); err != nil {
				t.Fatal(err)
			}
			disk := mustRead(t, p)
			if grew := len(disk) != len(orig); grew != path.resize {
				t.Fatalf("file size %d -> %d: wrong Flush path for %s", len(orig), len(disk), path.name)
			}

			// CopyTo emits exactly the file on disk.
			var buf bytes.Buffer
			if err := f.CopyTo(&buf); err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(buf.Bytes(), disk) {
				t.Fatalf("CopyTo after Flush differs from the file on disk at byte %d (len %d, disk %d)", firstDiff(buf.Bytes(), disk), buf.Len(), len(disk))
			}
			cp := filepath.Join(t.TempDir(), "copy.fits")
			if err := os.WriteFile(cp, buf.Bytes(), 0o644); err != nil {
				t.Fatal(err)
			}
			cf := testutil.CfrefDump(t, cp)
			for i, want := range map[int]string{1: "RENAMED", 3: "Test-ASCII", 4: "Test-BINTABLE", 5: "Test-BINTABLE"} {
				c := cf.HDUs[i].Card("EXTNAME")
				if c == nil || c.Str == nil || string(*c.Str) != want {
					t.Errorf("cfitsio reads HDU %d EXTNAME of the copy as %+v, want %q", i, c, want)
				}
			}

			// EditFile copying HDUs out of the flushed *File.
			dst, _ := copySample(t, "iter_a.fit")
			if err := EditFile(dst, func(_ *File, out *Writer) error {
				for i := range f.NumHDU() {
					h, err := f.HDU(i)
					if err != nil {
						return err
					}
					if err := out.CopyHDU(h); err != nil {
						return err
					}
				}
				return nil
			}); err != nil {
				t.Fatal(err)
			}
			if got := mustRead(t, dst); !bytes.Equal(got, disk) {
				t.Fatalf("EditFile copy of the flushed HDUs differs from the file on disk at byte %d", firstDiff(got, disk))
			}

			// HDUByName follows the rename.
			h, err := f.HDUByName("RENAMED")
			if err != nil || h.Index() != 1 {
				t.Fatalf("HDUByName(RENAMED) after Flush = %v, %v; want HDU 1", h, err)
			}
			if h, err := f.HDUByName("Test-BINTABLE"); err != nil || h.Index() != 4 {
				t.Fatalf("HDUByName(Test-BINTABLE) after Flush = %v, %v; want HDU 4 (HDU 1 was renamed)", h, err)
			}
			// A fresh Open agrees.
			r, err := Open(p)
			if err != nil {
				t.Fatal(err)
			}
			defer r.Close()
			if h, err := r.HDUByName("RENAMED"); err != nil || h.Index() != 1 {
				t.Fatalf("reopened HDUByName(RENAMED) = %v, %v", h, err)
			}
		})
	}
}

// HDUByName matches the EXTNAME the HDU's Header() holds now, flushed or
// not, as cfitsio's fits_movnam_hdu reads the EXTNAME keyword of each
// header it visits (an edit made with fits_update_key is in that header
// before it is written out). Before the edit is flushed the HDU is found
// by its new name and no longer by its old one; Close without Flush
// leaves the file, and so a reopened file's names, unchanged.
func TestHDUByNameSeesUnflushedRename(t *testing.T) {
	p, orig := copySample(t, "testprog.std")
	f, err := OpenForEdit(p)
	if err != nil {
		t.Fatal(err)
	}
	if err := hduHeader(t, f, 1).Set("EXTNAME", "RENAMED", ""); err != nil {
		t.Fatal(err)
	}
	if h, err := f.HDUByName("RENAMED"); err != nil || h.Index() != 1 {
		t.Fatalf("HDUByName(RENAMED) before Flush = %v, %v; want HDU 1", h, err)
	}
	if h, err := f.HDUByName("Test-BINTABLE"); err != nil || h.Index() != 4 {
		t.Fatalf("HDUByName(Test-BINTABLE) before Flush = %v, %v; want HDU 4, whose header still says so", h, err)
	}
	// Removing EXTNAME takes the HDU out of the name index too.
	if err := hduHeader(t, f, 3).Delete("EXTNAME"); err != nil {
		t.Fatal(err)
	}
	if h, err := f.HDUByName("Test-ASCII"); err == nil {
		t.Fatalf("HDUByName(Test-ASCII) after deleting its EXTNAME = HDU %d", h.Index())
	}
	f.Close()
	if !bytes.Equal(mustRead(t, p), orig) {
		t.Fatal("Close wrote the edit")
	}
	r, err := Open(p)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	if h, err := r.HDUByName("Test-BINTABLE"); err != nil || h.Index() != 1 {
		t.Fatalf("reopened HDUByName(Test-BINTABLE) = %v, %v; want HDU 1", h, err)
	}
}
