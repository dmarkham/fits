package fits

import (
	"path/filepath"
	"testing"
)

// BenchmarkReadPixelsFloat32_1k measures the dominant hot path: reading a
// 1024×1024 float32 image via ReadPixels[float32].
//
// Plan target: within 2× of cfitsio throughput on the same hardware.
func BenchmarkReadPixelsFloat32_1k(b *testing.B) {
	dir := b.TempDir()
	path := filepath.Join(dir, "bench.fits")
	data := make([]float32, 1024*1024)
	for i := range data {
		data[i] = float32(i)
	}
	f, _ := Create(path)
	WriteImage(f, nil, []int64{1024, 1024}, data)
	f.Close()

	b.ResetTimer()
	b.SetBytes(int64(len(data) * 4))
	for b.Loop() {
		fr, _ := Open(path)
		img, _ := fr.Primary()
		_, err := ReadPixels[float32](img)
		if err != nil {
			b.Fatal(err)
		}
		fr.Close()
	}
}

// BenchmarkWriteImageFloat32_1k measures the write path.
func BenchmarkWriteImageFloat32_1k(b *testing.B) {
	dir := b.TempDir()
	data := make([]float32, 1024*1024)
	b.ResetTimer()
	b.SetBytes(int64(len(data) * 4))
	for b.Loop() {
		path := filepath.Join(dir, "w.fits")
		f, _ := Create(path)
		WriteImage(f, nil, []int64{1024, 1024}, data)
		f.Close()
	}
}

// BenchmarkReadColumnInt32_1Mrow measures scalar column read throughput.
func BenchmarkReadColumnInt32_1Mrow(b *testing.B) {
	dir := b.TempDir()
	path := filepath.Join(dir, "col.fits")
	const n = 1_000_000
	ids := make([]int32, n)
	for i := range ids {
		ids[i] = int32(i)
	}
	f, _ := Create(path)
	WriteImage(f, nil, []int64{}, []float32{}) // primary
	AppendBinaryTable(f, nil, []ColumnData{{Name: "ID", DataInt32: ids}})
	f.Close()

	b.ResetTimer()
	b.SetBytes(n * 4)
	for b.Loop() {
		fr, _ := Open(path)
		hdu, _ := fr.HDU(1)
		tbl := hdu.(*BinaryTableHDU)
		_, err := ReadColumn[int32](tbl, 1)
		if err != nil {
			b.Fatal(err)
		}
		fr.Close()
	}
}

// BenchmarkEditHeaderInPlace measures the header-only in-place edit path.
// Target: milliseconds regardless of file size.
func BenchmarkEditHeaderInPlace(b *testing.B) {
	dir := b.TempDir()
	path := filepath.Join(dir, "edit.fits")
	data := make([]float32, 1024*1024)
	f, _ := Create(path)
	WriteImage(f, nil, []int64{1024, 1024}, data)
	f.Close()

	b.ResetTimer()
	for b.Loop() {
		fe, _ := OpenForEdit(path)
		p, _ := fe.Primary()
		// Flip a keyword value without changing length.
		p.Header().Set("OBJECT", "benchmark", "")
		fe.Flush()
		fe.Close()
	}
}
