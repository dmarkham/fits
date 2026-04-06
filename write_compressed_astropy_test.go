package fits

import (
	"encoding/json"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/compress"
)

// astropyPython returns the path to a Python interpreter with astropy
// installed. Set FITS_ASTROPY_PYTHON to override; the default looks for
// a venv in the sibling fits-manager project.
func astropyPython() string {
	if p := os.Getenv("FITS_ASTROPY_PYTHON"); p != "" {
		return p
	}
	return "python3" // fallback to system python
}

// TestWriteCompressedReadByAstropy is the interop acid test: we write a
// tile-compressed FITS file via AppendCompressedImage, then invoke
// astropy.io.fits (via the fits-manager venv) to decompress it, and
// compare the decoded pixel array against the original.
//
// This closes the loop — if astropy can read our compressed output,
// any other FITS-aware tool can too.
func TestWriteCompressedReadByAstropy(t *testing.T) {
	venvPython := astropyPython()
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

// TestWriteCompressedFloat32_FallbackReadByAstropy verifies that
// constant-tile files written via the GZIP_COMPRESSED_DATA fallback
// path are readable by astropy. The constant tile must round-trip
// bit-exactly (no quantization loss).
func TestWriteCompressedFloat32_FallbackReadByAstropy(t *testing.T) {
	venvPython := astropyPython()
	if _, err := os.Stat(venvPython); err != nil {
		t.Skip("fits-manager venv not available")
	}
	dir := t.TempDir()
	path := filepath.Join(dir, "float_fallback_astropy.fits")

	shape := []int64{64, 48}
	n := shape[0] * shape[1]
	pixels := make([]float32, n)
	for i := int64(0); i < n; i++ {
		if i/shape[0] < 32 {
			pixels[i] = 42.0
		} else {
			pixels[i] = 100.0 + 5.0*float32((i%7)-3)
		}
	}

	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, nil, []int64{}, []float32{}); err != nil {
		t.Fatal(err)
	}
	if _, err := AppendCompressedFloat32Image(f, nil, shape, pixels, CompressFloatOptions{
		Algorithm: compress.RICE1,
		TileShape: []int64{64, 16},
		QLevel:    4,
	}); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}

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
		Shape  []int     `json:"shape"`
		Pixels []float64 `json:"pixels"`
	}
	if err := json.Unmarshal(out, &result); err != nil {
		t.Fatalf("parse astropy output: %v\n%s", err, out)
	}
	if int64(len(result.Pixels)) != n {
		t.Fatalf("astropy returned %d pixels, expected %d", len(result.Pixels), n)
	}
	// Constant region must round-trip exactly (fallback path).
	for i := int64(0); i < 32*shape[0]; i++ {
		if result.Pixels[i] != 42.0 {
			t.Fatalf("constant region pixel %d: astropy got %g, we wrote 42", i, result.Pixels[i])
		}
	}
	// Noisy region within tolerance.
	worst := 0.0
	for i := 32 * shape[0]; i < n; i++ {
		diff := result.Pixels[i] - float64(pixels[i])
		if diff < 0 {
			diff = -diff
		}
		if diff > worst {
			worst = diff
		}
	}
	t.Logf("noisy region worst error (astropy-read): %g", worst)
	if worst > 10 {
		t.Errorf("noisy region too lossy: %g", worst)
	}
}

// TestWriteCompressedFloat32_ReadByAstropy closes the interop loop for
// the float quantization write path: we write a float32 tile-compressed
// image via AppendCompressedFloat32Image, then have astropy decompress
// it and compare the reconstructed pixel values against the originals
// within the expected quantization precision.
func TestWriteCompressedFloat32_ReadByAstropy(t *testing.T) {
	venvPython := astropyPython()
	if _, err := os.Stat(venvPython); err != nil {
		t.Skip("fits-manager venv not available")
	}

	type variant struct {
		name   string
		algo   compress.Algorithm
		dither compress.DitherMethod
	}
	cases := []variant{
		{"RICE_1_dither1", compress.RICE1, compress.SubtractiveDither1},
		{"RICE_1_dither2", compress.RICE1, compress.SubtractiveDither2},
		{"GZIP_1_dither1", compress.GZIP1, compress.SubtractiveDither1},
		{"GZIP_2_dither1", compress.GZIP2, compress.SubtractiveDither1},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			dir := t.TempDir()
			path := filepath.Join(dir, "go_"+c.name+".fits")

			// Deterministic Gaussian noise + a few exact zeros so the
			// dither_method=2 case has something to preserve.
			shape := []int64{64, 48}
			n := shape[0] * shape[1]
			pixels := make([]float32, n)
			seed := uint64(42)
			for i := int64(0); i < n; i++ {
				seed = seed*6364136223846793005 + 1442695040888963407
				u1 := float64((seed>>33)&0xFFFFFFFF) / 4294967296.0
				seed = seed*6364136223846793005 + 1442695040888963407
				u2 := float64((seed>>33)&0xFFFFFFFF) / 4294967296.0
				if u1 < 1e-12 {
					u1 = 1e-12
				}
				// Box-Muller
				g := 5.0 * math.Sqrt(-2*math.Log(u1)) * math.Cos(2*math.Pi*u2)
				pixels[i] = float32(100.0 + g)
			}
			// Scatter some exact zeros (to test dither method 2).
			if c.dither == compress.SubtractiveDither2 {
				for i := int64(0); i < n; i += 57 {
					pixels[i] = 0
				}
			}

			f, err := Create(path)
			if err != nil {
				t.Fatal(err)
			}
			if _, err := WriteImage(f, nil, []int64{}, []float32{}); err != nil {
				t.Fatal(err)
			}
			if _, err := AppendCompressedFloat32Image(f, nil, shape, pixels, CompressFloatOptions{
				Algorithm:    c.algo,
				TileShape:    []int64{64, 16},
				QLevel:       4,
				DitherMethod: c.dither,
				ZDither0:     1,
			}); err != nil {
				t.Fatal(err)
			}
			if err := f.Close(); err != nil {
				t.Fatal(err)
			}

			// astropy reads, dumps pixel values as JSON.
			script := `
import sys, json
import astropy.io.fits as fits
import numpy as np
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
				Shape  []int     `json:"shape"`
				Pixels []float64 `json:"pixels"`
			}
			if err := json.Unmarshal(out, &result); err != nil {
				t.Fatalf("parse astropy output: %v\n%s", err, out)
			}
			if int64(len(result.Pixels)) != n {
				t.Fatalf("astropy returned %d pixels, expected %d", len(result.Pixels), n)
			}

			// Tolerance: delta = sigma/qlevel ~ 5/4 = 1.25; with
			// dither the reconstruction error is bounded by delta/2
			// on average and occasionally larger. 2*delta is a safe
			// upper bound.
			tolerance := 2.5
			worst := 0.0
			var nbad int
			for i := int64(0); i < n; i++ {
				got := result.Pixels[i]
				want := float64(pixels[i])
				diff := got - want
				if diff < 0 {
					diff = -diff
				}
				if diff > worst {
					worst = diff
				}
				if diff > tolerance {
					nbad++
				}
			}
			if nbad > 0 {
				t.Errorf("%s: %d/%d pixels outside tolerance %g (worst %g)",
					c.name, nbad, n, tolerance, worst)
			}
			// For dither2, zero pixels must be exactly preserved.
			if c.dither == compress.SubtractiveDither2 {
				for i := int64(0); i < n; i += 57 {
					if result.Pixels[i] != 0 {
						t.Errorf("dither2: zero at pixel %d round-tripped to %g", i, result.Pixels[i])
						break
					}
				}
			}
			t.Logf("%s: worst reconstruction error %g (tolerance %g)", c.name, worst, tolerance)
		})
	}
}

