//go:build amd64 && !amd64.v3

package stats

// fusionModel names how the Go compiler treats the multiply-adds of the
// frozen MeanStdev (which has no explicit float64() conversions) on this
// build; gen_golden.py records the one-pass result for each model. On
// amd64 below GOAMD64=v3 there is no FMA instruction: nothing is fused.
const fusionModel = "none"
