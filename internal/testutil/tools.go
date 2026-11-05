package testutil

import (
	"regexp"
	"strconv"
	"testing"
)

// FitsverifyResult is the outcome of running fitsverify on one file.
type FitsverifyResult struct {
	Warnings int
	Errors   int
	Report   string // full fitsverify output, for t.Log on failure
}

// OK reports whether fitsverify found neither warnings nor errors.
func (r FitsverifyResult) OK() bool { return r.Warnings == 0 && r.Errors == 0 }

var fitsverifySummary = regexp.MustCompile(`\*\*\*\* Verification found (\d+) warning\(s\) and (\d+) error\(s\)\. \*\*\*\*`)

// Fitsverify is HEASARC fitsverify's verdict on path: the recorded
// reading of these exact bytes in tier 1, fitsverify run live in tier 2
// (Run). fitsverify exits with (errors + warnings), so a non-zero exit is
// expected for a file with findings. Output that has no summary line
// fails the test (fitsverify could not judge the file at all).
func Fitsverify(t testing.TB, path string) FitsverifyResult {
	t.Helper()
	r := Run(t, ToolFitsverify, Call{Args: []string{path}, Inputs: []string{path}})
	return parseFitsverify(t, path, string(r.Stdout)+string(r.Stderr))
}

func parseFitsverify(t testing.TB, path, report string) FitsverifyResult {
	t.Helper()
	m := fitsverifySummary.FindAllStringSubmatch(report, -1)
	if len(m) != 1 {
		t.Fatalf("fitsverify output for %s has %d summary lines, want 1:\n%s", path, len(m), report)
	}
	w, _ := strconv.Atoi(m[0][1])
	e, _ := strconv.Atoi(m[0][2])
	return FitsverifyResult{Warnings: w, Errors: e, Report: report}
}
