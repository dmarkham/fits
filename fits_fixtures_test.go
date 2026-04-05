package fits

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// fixtures lists the cfitsio reference files that the library's regression
// suite must round-trip byte-for-byte (plan §"Done criteria" item 5).
var fixtures = []string{
	"iter_image.fit",
	"iter_a.fit",
	"iter_b.fit",
	"iter_c.fit",
	"vari.fits",
	"testprog.std",
	"testf77.std",
}

// TestRoundTripFixtures opens each fixture, walks every HDU through
// (*File).CopyTo, and asserts the output equals the input.
func TestRoundTripFixtures(t *testing.T) {
	for _, name := range fixtures {
		name := name
		t.Run(name, func(t *testing.T) {
			path := filepath.Join("testdata", name)
			orig, err := os.ReadFile(path)
			if err != nil {
				t.Skipf("fixture missing: %v", err)
				return
			}

			f, err := Open(path)
			if err != nil {
				t.Fatalf("Open: %v", err)
			}
			defer f.Close()

			var buf bytes.Buffer
			if err := f.CopyTo(&buf); err != nil {
				t.Fatalf("CopyTo: %v", err)
			}
			out := buf.Bytes()
			if !bytes.Equal(out, orig) {
				// Point out where the first difference is so a human can dig in.
				n := len(orig)
				if len(out) < n {
					n = len(out)
				}
				var diff int64 = -1
				for i := 0; i < n; i++ {
					if out[i] != orig[i] {
						diff = int64(i)
						break
					}
				}
				t.Fatalf("byte mismatch: lens %d (out) vs %d (orig), first diff at offset %d", len(out), len(orig), diff)
			}
		})
	}
}
