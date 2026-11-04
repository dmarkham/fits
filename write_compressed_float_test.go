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
// reconstructed from a fallback tile must equal the original bit-for-bit;
// the quantized tile must be within half a quantization step.
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
	// Noisy region (tile 3): within half a quantization step of the
	// input, the bound fits_quantize_float guarantees (ZSCALE of tile 3).
	tc := readTableCols(t, cimg)
	if len(tc.f64["ZSCALE"]) != 3 || tc.f64["ZSCALE"][0] != 0 || tc.f64["ZSCALE"][1] != 0 || tc.f64["ZSCALE"][2] <= 0 {
		t.Fatalf("ZSCALE %v: want fallback rows 1-2 = 0 (as cfitsio leaves them) and a positive step for row 3", tc.f64["ZSCALE"])
	}
	if d := tc.data["GZIP_COMPRESSED_DATA"]; len(d) != 3 || len(d[0]) == 0 || len(d[1]) == 0 || len(d[2]) != 0 {
		t.Fatalf("GZIP_COMPRESSED_DATA: want cells for tiles 1-2 only")
	}
	bound := 0.5*tc.f64["ZSCALE"][2]*(1+1e-12) + 1e-5
	for i := 32 * shape[0]; i < n; i++ {
		if diff := math.Abs(float64(got[i]) - float64(pixels[i])); diff > bound {
			t.Fatalf("noisy pixel %d: |%g - %g| = %g > half step %g", i, got[i], pixels[i], diff, bound)
		}
	}
}

// TestWriteCompressedFloat32_Roundtrip writes a float32 image with
// tile compression via AppendCompressedFloat32Image, then reads it
// back with Open + ReadPixelsCompressed[float32] and verifies every pixel
// is within half a quantization step (0.5*ZSCALE of its tile) of the
// input. External readers of the same files are checked in
// TestCompressedWriterInterop.
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

			// Every pixel must be within half a quantization step of the
			// input: the bound of fits_quantize_float with dithering
			// (ZSCALE is that step, per 16-row tile).
			tcols := readTableCols(t, cimg)
			zs := tcols.f64["ZSCALE"]
			if len(zs) != 3 {
				t.Fatalf("ZSCALE has %d rows, want 3", len(zs))
			}
			worst := 0.0
			for i := int64(0); i < n; i++ {
				bound := 0.5*zs[i/(shape[0]*16)]*(1+1e-12) + 1e-5
				diff := math.Abs(float64(got[i]) - float64(pixels[i]))
				if diff > worst {
					worst = diff
				}
				if diff > bound {
					t.Fatalf("%s: pixel %d: |%g - %g| = %g > half step %g", tc.name, i, got[i], pixels[i], diff, bound)
				}
			}
			t.Logf("%s: worst reconstruction error %g; half steps per tile %g %g %g", tc.name, worst, zs[0]/2, zs[1]/2, zs[2]/2)
		})
	}
}
