//go:build arm64

package stats

// fusionModel: see fusion_amd64_test.go. On arm64 Go fuses both the
// accumulation `sum2 += x*x` (FMADD) and `sum2 - float64(n)*mean*mean`
// (FMSUB) (measured under qemu-aarch64: TestMeanStdev_OnePassFormula
// passes bit for bit with onepass_stdev_fma_all on every golden case).
const fusionModel = "fma_all"
