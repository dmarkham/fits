//go:build interop

package stretch

import (
	"os"
	"os/exec"
	"runtime"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestFMAGuards_GOAMD64v3: the explicit float32(...) conversions around
// products that feed an addition (ghtSetup, the Lab conversions) exist so
// the compiler cannot fuse them into an FMA, which would round once where
// Siril's x86-64 build rounds twice. The default amd64 build (GOAMD64=v1)
// has no FMA instruction, so there a missing conversion changes nothing and
// the oracle tests cannot see it. With GOAMD64=v3 the compiler does fuse
// x*y + z, so this test re-runs the package's tests in that mode: a
// dropped conversion then fails TestGHTSweep_SirilOracle and the other
// oracle tests (each of the mutants c.b4, q1, c2 and b1 of the B = -1
// setup without its conversion passes under v1 and fails here). On arm64
// the ordinary run already fuses, so it needs no second run.
//
// Tier 2 (it starts the go tool): it needs the go tool (a failure without
// it) and an x86-64-v3 CPU (AVX2, BMI1/2, FMA, F16C, LZCNT, MOVBE; a
// platform skip without one).
func TestFMAGuards_GOAMD64v3(t *testing.T) {
	if os.Getenv("FITS_FMA_CHILD") != "" {
		t.Skip("running inside the GOAMD64=v3 child")
	}
	if runtime.GOARCH != "amd64" {
		t.Skipf("GOARCH=%s: not amd64; arm64 fuses multiply-adds in the ordinary run", runtime.GOARCH)
	}
	if testing.Short() {
		t.Skip("-short: skips the GOAMD64=v3 rebuild")
	}
	goTool := testutil.RequireTool(t, "go")
	cpu, err := os.ReadFile("/proc/cpuinfo")
	if err != nil {
		t.Skipf("PLATFORM SKIP: cannot read /proc/cpuinfo to check for x86-64-v3: %v", err)
	}
	flags := map[string]bool{}
	for _, line := range strings.Split(string(cpu), "\n") {
		if k, v, ok := strings.Cut(line, ":"); ok && strings.TrimSpace(k) == "flags" {
			for _, f := range strings.Fields(v) {
				flags[f] = true
			}
			break
		}
	}
	for _, f := range []string{"avx2", "bmi1", "bmi2", "fma", "f16c", "abm", "movbe"} {
		if !flags[f] {
			t.Skipf("PLATFORM SKIP: CPU lacks %s, cannot run GOAMD64=v3 code", f)
		}
	}
	out, err := childTest(goTool, []string{"GOAMD64=v3"})
	if err != nil {
		t.Fatalf("package tests fail with GOAMD64=v3 (a multiply-add was fused): %v\n%s", err, out)
	}
	t.Logf("GOAMD64=v3: %s", strings.TrimSpace(string(out)))
}

// childTest runs this package's tier 1 tests in a child `go test` (no
// interop tag, so neither the tests in this file nor any other tool run)
// with extra environment settings.
func childTest(goTool string, set []string, args ...string) ([]byte, error) {
	cmd := exec.Command(goTool, append([]string{"test", "-count=1"}, append(args, ".")...)...)
	env := append([]string{"FITS_FMA_CHILD=1"}, set...)
	for _, e := range os.Environ() {
		k, _, _ := strings.Cut(e, "=")
		drop := false
		for _, s := range set {
			if sk, _, _ := strings.Cut(s, "="); sk == k {
				drop = true
			}
		}
		if !drop {
			env = append(env, e)
		}
	}
	cmd.Env = env
	return cmd.CombinedOutput()
}

// TestCrossArch_arm64: re-runs this package's tests built for arm64 under
// qemu-aarch64-static. Go's float-to-int conversion of an out-of-range
// value is implementation-specific and differs there (it saturates; amd64
// gives the minimum integer), and arm64 fuses multiply-adds by default, so
// code that relied on amd64's conversions or that lost an FMA guard
// passes on amd64 and fails here. Example: cvttsd2si's range test written
// v <= 2^31 instead of v < 2^31 passes on amd64 (Go's own int32(2^31) is
// MinInt32 there) and fails TestCLAHEClipConversion_SirilOracle here.
//
// Tier 2: needs the go tool and qemu-aarch64-static (a failure without
// them).
func TestCrossArch_arm64(t *testing.T) {
	if os.Getenv("FITS_FMA_CHILD") != "" {
		t.Skip("running inside a child build")
	}
	if runtime.GOARCH == "arm64" {
		t.Skip("already running on arm64")
	}
	if testing.Short() {
		t.Skip("-short: skips the arm64 rebuild")
	}
	goTool := testutil.RequireTool(t, "go")
	qemu := testutil.RequireTool(t, "qemu-aarch64-static")
	out, err := childTest(goTool, []string{"GOARCH=arm64", "CGO_ENABLED=0"}, "-exec", qemu)
	if err != nil {
		t.Fatalf("package tests fail built for arm64 (run under %s): %v\n%s", qemu, err, out)
	}
	t.Logf("GOARCH=arm64: %s", strings.TrimSpace(string(out)))
}
