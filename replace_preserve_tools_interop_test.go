//go:build interop && (linux || darwin)

package fits

import (
	"os/exec"
	"strings"
	"testing"
)

func toolOutput(t *testing.T, name string, args ...string) string {
	t.Helper()
	out, err := exec.Command(name, args...).CombinedOutput()
	if err != nil {
		t.Fatalf("%s %v: %v\n%s", name, args, err, out)
	}
	// Drop the "# file: <path>" line, which names the path.
	var keep []string
	for _, l := range strings.Split(string(out), "\n") {
		if !strings.HasPrefix(l, "# file:") {
			keep = append(keep, l)
		}
	}
	return strings.Join(keep, "\n")
}

// requireTools fails the test (tier 2) when one of the tools used as
// the truth here is missing.
func requireTools(t *testing.T, tools ...string) {
	t.Helper()
	for _, tool := range tools {
		if _, err := exec.LookPath(tool); err != nil {
			t.Fatalf("tier 2 (-tags interop) needs %s as the truth for this test: %v", tool, err)
		}
	}
}
