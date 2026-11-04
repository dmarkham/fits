//go:build amd64.v3

package stats

// fusionModel: see fusion_amd64_test.go. GOAMD64=v3 has FMA; Go 1.27
// fuses the accumulation `sum2 += x*x` but not the final subtraction
// (measured: TestMeanStdev_OnePassFormula passes bit for bit with
// onepass_stdev_fma_acc on every golden case).
const fusionModel = "fma_acc"
