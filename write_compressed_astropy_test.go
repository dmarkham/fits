package fits

import (
	"encoding/json"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/compress"
)

// TestWriteCompressedReadByAstropy is the interop acid test: we write a
// tile-compressed FITS file via AppendCompressedImage, then invoke
// astropy.io.fits (via the fits-manager venv) to decompress it, and
// compare the decoded pixel array against the original.
//
// This closes the loop — if astropy can read our compressed output,
// any other FITS-aware tool can too.
func TestWriteCompressedReadByAstropy(t *testing.T) {
	venvPython := "/home/dmarkham/git/dmarkham/fits-manager/venv/bin/python"
	if _, err := os.Stat(venvPython); err != nil {
		t.Skip("fits-manager venv not available")
	}

	cases := []struct {
		name string
		algo compress.Algorithm
	}{
		{"RICE_1", compress.RICE1},
		{"GZIP_1", compress.GZIP1},
		{"GZIP_2", compress.GZIP2},
		{"NOCOMPRESS", compress.NoCompress},
		{"HCOMPRESS_1", compress.HCOMPRESS1},
		{"PLIO_1", compress.PLIO1},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dir := t.TempDir()
			path := filepath.Join(dir, "go_"+c.name+".fits")

			// Write a deterministic int16 image via the Go library.
			// Each algorithm has domain constraints; pick test data
			// that matches the algorithm's intended use case.
			shape := []int64{64, 48}
			pixels := make([]int16, 64*48)
			switch c.algo {
			case compress.PLIO1:
				// PLIO is a sparse mask codec. Sparse-run pattern.
				// NAXIS1=64 is the fast axis → 64 columns per row.
				const ncols = 64
				for row := 10; row < 20; row++ {
					for col := 5; col < 40; col++ {
						pixels[row*ncols+col] = 1
					}
				}
				for row := 30; row < 45; row++ {
					for col := 10; col < 30; col++ {
						pixels[row*ncols+col] = 2
					}
				}
			default:
				for i := range pixels {
					pixels[i] = int16(i % 1000) // non-negative, max 999
				}
			}
			opts := CompressOptions{Algorithm: c.algo}
			// HCOMPRESS needs 2D tiles; use the whole image as one tile
			// so astropy's "expects two dimensional tiles" check passes.
			if c.algo == compress.HCOMPRESS1 {
				opts.TileShape = []int64{64, 48}
			}

			f, err := Create(path)
			if err != nil {
				t.Fatal(err)
			}
			if _, err := WriteImage(f, nil, []int64{}, []float32{}); err != nil {
				t.Fatal(err)
			}
			if _, err := AppendCompressedImage(f, nil, shape, pixels, opts); err != nil {
				t.Fatal(err)
			}
			if err := f.Close(); err != nil {
				t.Fatal(err)
			}

			// Read via astropy, dump pixels as JSON.
			script := `
import sys, json
import astropy.io.fits as fits
with fits.open(sys.argv[1]) as h:
    d = h[1].data
    print(json.dumps({"shape": list(d.shape), "pixels": d.flatten().tolist()}))
`
			cmd := exec.Command(venvPython, "-c", script, path)
			out, err := cmd.CombinedOutput()
			if err != nil {
				t.Fatalf("astropy read failed: %v\n%s", err, out)
			}
			var result struct {
				Shape  []int   `json:"shape"`
				Pixels []int64 `json:"pixels"`
			}
			if err := json.Unmarshal(out, &result); err != nil {
				t.Fatalf("parse astropy output: %v\n%s", err, out)
			}
			if len(result.Pixels) != len(pixels) {
				t.Fatalf("astropy returned %d pixels, expected %d", len(result.Pixels), len(pixels))
			}
			for i := range pixels {
				if int64(pixels[i]) != result.Pixels[i] {
					t.Fatalf("pixel %d: astropy got %d, we wrote %d", i, result.Pixels[i], pixels[i])
				}
			}
		})
	}
}
