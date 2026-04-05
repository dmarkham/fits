package fits

import (
	"path/filepath"
	"testing"
)

// Decompression benchmarks — measure end-to-end file-open → decompress
// throughput for each algorithm on the astropy test fixtures. Baseline
// target: within 2× of cfitsio/astropy throughput on the same hardware.

func BenchmarkDecompressRICE1_i16(b *testing.B) {
	benchDecompressI16(b, "comp_rice_i16.fits")
}

func BenchmarkDecompressGZIP1_i16(b *testing.B) {
	benchDecompressI16(b, "comp_gzip1_i16.fits")
}

func BenchmarkDecompressGZIP2_i16(b *testing.B) {
	benchDecompressI16(b, "comp_gzip2_i16.fits")
}

func BenchmarkDecompressNoCompress_i16(b *testing.B) {
	benchDecompressI16(b, "comp_nocompress_i16.fits")
}

func BenchmarkDecompressRICE1_f32(b *testing.B) {
	benchDecompressF32(b, "comp_rice_f32_nodither.fits")
}

func BenchmarkDecompressGZIP2_f32(b *testing.B) {
	benchDecompressF32(b, "comp_gzip2_f32.fits")
}

func benchDecompressI16(b *testing.B, name string) {
	path := filepath.Join("testdata", name)
	b.ResetTimer()
	b.SetBytes(100 * 100 * 2)
	for b.Loop() {
		f, err := Open(path)
		if err != nil {
			b.Fatal(err)
		}
		h, _ := f.HDU(1)
		cimg := h.(*CompressedImageHDU)
		if _, err := ReadPixelsCompressed[int16](cimg); err != nil {
			b.Fatal(err)
		}
		f.Close()
	}
}

func benchDecompressF32(b *testing.B, name string) {
	path := filepath.Join("testdata", name)
	b.ResetTimer()
	b.SetBytes(100 * 100 * 4)
	for b.Loop() {
		f, err := Open(path)
		if err != nil {
			b.Fatal(err)
		}
		h, _ := f.HDU(1)
		cimg := h.(*CompressedImageHDU)
		if _, err := ReadPixelsCompressed[float32](cimg); err != nil {
			b.Fatal(err)
		}
		f.Close()
	}
}
