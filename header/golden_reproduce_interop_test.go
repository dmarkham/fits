//go:build interop

package header

import (
	"bytes"
	"encoding/json"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestGoldensReproduce: the committed goldens are exactly what the live
// harnesses print for the committed samples (cfref dump, hdrref longstr,
// hdrref testkey), so a golden cannot drift from its generator.
func TestGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	cfref := testutil.CfrefBinary(t)
	check := func(golden string, got []byte) {
		t.Helper()
		testutil.RequireFile(t, golden)
		want, _ := os.ReadFile(golden)
		if !bytes.Equal(got, want) {
			t.Errorf("%s does not match the live harness output (regenerate with make -C header/testdata goldens)", golden)
		}
	}
	for _, n := range headerSamples {
		p := headerTestdata(t, n)
		out, err := exec.Command(cfref, "dump", p).Output()
		if err != nil {
			t.Fatalf("cfref dump %s: %v", n, err)
		}
		// The Makefile runs the harness on the relative name.
		out = bytes.ReplaceAll(out, []byte(`"file": "`+p+`"`), []byte(`"file": "`+n+`"`))
		check(headerTestdata(t, "golden", n+".json"), out)
		check(headerTestdata(t, "golden", n+".longstr.json"), liveHdrref(t, "longstr", p))
	}
	gl, _ := filepath.Glob(headerTestdata(t, "golden", "testdata", "*.longstr.json"))
	for _, g := range gl {
		name := strings.TrimSuffix(filepath.Base(g), ".longstr.json")
		check(g, liveHdrref(t, "longstr", testutil.TestdataPath(t, name)))
	}
	var keys []string
	b, _ := os.ReadFile(headerTestdata(t, "golden", "testkey.json"))
	var kk []struct {
		Key testutil.Latin1 `json:"key"`
	}
	if err := json.Unmarshal(b, &kk); err != nil {
		t.Fatal(err)
	}
	for _, k := range kk {
		keys = append(keys, string(k.Key))
	}
	check(headerTestdata(t, "golden", "testkey.json"), liveHdrref(t, append([]string{"testkey"}, keys...)...))
}

// TestHeaderSamplesReproduce: header/testdata `make verify` regenerates
// the samples (gen_samples.py: malformed, continue, adversarial; hdrref
// write: cfwritten.fits) and every golden (cfref dump, hdrref longstr for
// the samples and every repo testdata/ file, hdrref testkey) into a
// scratch directory and byte-compares them with the committed files.
// TestGoldensReproduce checks the goldens only.
func TestHeaderSamplesReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	testutil.MakeVerify(t, "testdata", []string{"hdrref/hdrref"}, []string{"verify"})
}

// liveHdrref runs the built hdrref harness (regeneration checks compare
// its live output with the committed goldens; nothing is recorded).
func liveHdrref(t testing.TB, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(testutil.Harness(t, toolHdrref.Harness), args...)
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("hdrref %v: %v: %s", args, err, strings.TrimSpace(stderr.String()))
	}
	return out
}
