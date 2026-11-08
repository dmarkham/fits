//go:build interop

package fits

import (
	"bytes"
	"encoding/json"
	"os/exec"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// The committed goldens are exactly what the cfitsio harness says today
// (byte for byte; `make -C testdata/tables verify` also regenerates the
// samples themselves and compares them).
func TestTableGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	bin := testutil.CfrefBinary(t)
	for _, s := range tableSamples(t) {
		if !strings.HasPrefix(s.name, "tables/") {
			continue // covered by internal/testutil.TestCfrefDumpReproducesGoldens
		}
		live, err := exec.Command(bin, "dump", s.path).Output()
		if err != nil {
			t.Fatalf("%s: cfref dump: %v", s.name, err)
		}
		want, err := testutil.ReadGolden(s.golden)
		if err != nil {
			t.Fatal(err)
		}
		// The harness records the base name only, so the dump does not
		// depend on where the sample lives. The goldens are stored
		// compact (compact_json.py: the whitespace between tokens
		// dropped, every token byte for byte), which json.Compact
		// reproduces exactly.
		var compact bytes.Buffer
		if err := json.Compact(&compact, live); err != nil {
			t.Fatalf("%s: cfref dump is not JSON: %v", s.name, err)
		}
		if !bytes.Equal(compact.Bytes(), want) {
			t.Errorf("%s: live cfitsio dump differs from the committed golden (regenerate with make -C testdata/tables goldens and review)", s.name)
		}
	}
}
