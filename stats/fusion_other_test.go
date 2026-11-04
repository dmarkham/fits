//go:build !amd64 && !arm64

package stats

// fusionModel: see fusion_amd64_test.go. No measured model for this
// architecture: the bit-exact one-pass pin is skipped (loudly) and only
// the rounding-bound tests apply.
const fusionModel = ""
