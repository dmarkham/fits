package fits

import (
	"math"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/compress"
)

// TestWriteCompressedFloat32_GzipFallback exercises the
// GZIP_COMPRESSED_DATA fallback column: when a tile cannot be quantized
// (constant value / zero noise), the writer should store the raw float
// bytes gzipped in a secondary column. The reader then routes that
// column back through the gzip decoder and skips dequantization. Pixels
// reconstructed from a fallback tile must equal the original bit-for-bit.
func TestWriteCompressedFloat32_GzipFallback(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "float_fallback.fits")

	// 64x48 image with a constant region (tiles 0-1) and a noisy
	// region (tiles 2). 16-row tiles -> 3 tiles total.
	shape := []int64{64, 48}
	n := shape[0] * shape[1]
	pixels := make([]float32, n)
	// First two tiles (rows 0..31) = constant 42.0. Third tile (rows
	// 32..47) = Gaussian noise around 100.
	for i := int64(0); i < n; i++ {
		row := i / shape[0]
		if row < 32 {
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
		Algorithm:    compress.RICE1,
		TileShape:    []int64{64, 16}, // 3 tiles
		QLevel:       4,
		DitherMethod: compress.SubtractiveDither1,
	}); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}

	fr, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fr.Close()
	h, _ := fr.HDU(1)
	cimg := h.(*CompressedImageHDU)
	got, err := ReadPixelsCompressed[float32](cimg)
	if err != nil {
		t.Fatalf("decompress: %v", err)
	}
	// Constant-region tiles should round-trip exactly because they
	// went through the raw-float gzip fallback (no quantization loss).
	for i := int64(0); i < 32*shape[0]; i++ {
		if got[i] != pixels[i] {
			t.Fatalf("constant region pixel %d: got %g want %g (fallback not exact)",
				i, got[i], pixels[i])
		}
	}
	// Noisy region should round-trip within quantization tolerance.
	worst := float32(0)
	for i := 32 * shape[0]; i < n; i++ {
		diff := got[i] - pixels[i]
		if diff < 0 {
			diff = -diff
		}
		if diff > worst {
			worst = diff
		}
	}
	t.Logf("noisy region worst reconstruction error: %g", worst)
	if worst > 10 {
		t.Errorf("noisy region reconstruction too lossy: %g", worst)
	}
}

// TestWriteCompressedFloat32_Roundtrip writes a float32 image with
// tile compression via AppendCompressedFloat32Image, then reads it
// back with Open + ReadPixelsCompressed[float32] and verifies pixel
// values match within the expected quantization precision
// (~sigma/qlevel ~ 5/4 ~ 1.25).
func TestWriteCompressedFloat32_Roundtrip(t *testing.T) {
	algos := []struct {
		name string
		algo compress.Algorithm
	}{
		{"RICE_1", compress.RICE1},
		{"GZIP_1", compress.GZIP1},
		{"GZIP_2", compress.GZIP2},
	}
	for _, tc := range algos {
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			path := filepath.Join(dir, "float_"+tc.name+".fits")

			// 64x48 deterministic Gaussian-ish image with a mean of 100
			// and sigma ~5. Use a simple LCG so results are reproducible.
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
				g := math.Sqrt(-2*math.Log(u1)) * math.Cos(2*math.Pi*u2)
				pixels[i] = float32(100 + 5*g)
			}

			f, err := Create(path)
			if err != nil {
				t.Fatal(err)
			}
			if _, err := WriteImage(f, nil, []int64{}, []float32{}); err != nil {
				t.Fatal(err)
			}
			// Use a 2D tile shape explicitly — whole-row tiles are 1D
			// after the first axis and some read-side paths prefer a
			// proper 2D tile.
			if _, err := AppendCompressedFloat32Image(f, nil, shape, pixels, CompressFloatOptions{
				Algorithm:    tc.algo,
				TileShape:    []int64{shape[0], 16},
				QLevel:       4,
				DitherMethod: compress.SubtractiveDither1,
				ZDither0:     1,
			}); err != nil {
				t.Fatal(err)
			}
			if err := f.Close(); err != nil {
				t.Fatal(err)
			}

			// Read back.
			fr, err := Open(path)
			if err != nil {
				t.Fatal(err)
			}
			defer fr.Close()
			h, _ := fr.HDU(1)
			cimg, ok := h.(*CompressedImageHDU)
			if !ok {
				t.Fatalf("HDU(1) is %T, want *CompressedImageHDU", h)
			}
			if cimg.CompressionType() != tc.name {
				t.Fatalf("algorithm %q, want %s", cimg.CompressionType(), tc.name)
			}
			if cimg.BITPIX() != -32 {
				t.Fatalf("BITPIX %d, want -32 (FLOAT_IMG)", cimg.BITPIX())
			}
			got, err := ReadPixelsCompressed[float32](cimg)
			if err != nil {
				t.Fatalf("decompress: %v", err)
			}
			if int64(len(got)) != n {
				t.Fatalf("len %d want %d", len(got), n)
			}

			// The quantization step is delta = sigma/qlevel ~= 1.25.
			// With dither, reconstruction error per pixel is bounded
			// by ~0.5*delta ~ 0.625. Allow 2*delta as a safety margin.
			tolerance := float32(2.5)
			worst := float32(0)
			nbad := 0
			for i := int64(0); i < n; i++ {
				diff := got[i] - pixels[i]
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
				t.Errorf("%s: %d pixels outside tolerance %g (worst %g)",
					tc.name, nbad, tolerance, worst)
			}
			t.Logf("%s: worst reconstruction error %g (tolerance %g)",
				tc.name, worst, tolerance)
		})
	}
}
