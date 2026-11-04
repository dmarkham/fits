//go:build interop

package testutil

import (
	"bytes"
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
)

// Tier 2: run the tool live. On the reference platform (linux/amd64, where
// the committed readings are defined, see ReferencePlatform) the live
// reading is what the test asserts on, and it must equal the recording,
// so the committed readings tier 1 replays are the tools' real ones. On
// any other platform the tools still run live, every live reading that
// differs from the reference recording is reported (t.Logf, "REFERENCE
// TOOL DIFFERS"), and the test asserts on the recorded reference reading:
// those differences belong to the platform's build of the tool, and the
// library is checked against the reference. With FITS_RECORD=1 (reference
// platform only) the reading is (re)recorded once the test has passed.
func execute(t testing.TB, tool Tool, c *Call, r *record, p placeholders) *record {
	t.Helper()
	bin, pre := toolCommand(t, tool)
	cmd := exec.Command(bin, append(pre, c.Args...)...)
	cmd.Dir = c.Dir
	if c.Stdin != nil {
		cmd.Stdin = bytes.NewReader(c.Stdin)
	}
	var stdout, stderr bytes.Buffer
	cmd.Stdout, cmd.Stderr = &stdout, &stderr
	err := cmd.Run()
	var exitErr *exec.ExitError
	switch {
	case err == nil:
	case errors.As(err, &exitErr):
		r.Exit = exitErr.ExitCode()
	default:
		t.Fatalf("run %s %v: %v", tool.Name, c.Args, err)
	}
	if c.Compact != nil {
		r.liveStdout = []byte(p.hide(stdout.String()))
		r.setStdout([]byte(p.hide(string(c.Compact(stdout.Bytes())))))
	} else {
		r.setStdout([]byte(p.hide(stdout.String())))
	}
	r.Stderr = p.hide(stderr.String())
	for _, o := range c.Outputs {
		b, err := os.ReadFile(o)
		if err != nil {
			t.Fatalf("%s %v did not write %s: %v\nstderr: %s", tool.Name, c.Args, o, err, stderr.String())
		}
		if c.FixOutput != nil {
			if err := os.WriteFile(o, c.FixOutput(b), 0o644); err != nil {
				t.Fatal(err)
			}
		}
		r.Outputs = append(r.Outputs, fileDigest(t, o, !c.OutputsByDigest))
	}
	r.RecordedOn = runtime.GOOS + "/" + runtime.GOARCH
	r.Test = t.Name()
	r.Versions = toolVersions(t, tool)

	rec := readRecord(t, r)
	recording := os.Getenv(EnvRecord) == "1"
	if recording && !ReferencePlatform() {
		t.Fatalf("%s=1: readings are recorded on %s only (this is %s/%s)", EnvRecord, referencePlatform, runtime.GOOS, runtime.GOARCH)
	}
	switch {
	case rec == nil && !recording && !ReferencePlatform():
		t.Errorf("no recorded %s reading %v of %s for this input; on %s/%s that means an earlier step of the test produced different bytes than on %s (a library difference, not a tool one)", tool.Name, r.Args, inputNames(r), runtime.GOOS, runtime.GOARCH, referencePlatform)
	case rec == nil && !recording:
		t.Errorf("no recorded %s reading %v of %s (tier 1 replays it): record it with\n  %s", tool.Name, r.Args, inputNames(r), recordHint(t))
	case rec != nil:
		if d := sameReading(r, rec); d != "" {
			if !ReferencePlatform() {
				t.Logf("REFERENCE TOOL DIFFERS: live %s reading %v of %s on %s/%s differs from the %s recording (%s); asserting on the recording",
					tool.Name, r.Args, inputNames(r), runtime.GOOS, runtime.GOARCH, referencePlatform, d)
				// The reference reading's outputs replace the live ones (Run
				// writes them back, or the digest stub tier 1 uses), so a
				// later step reads what the reference tool wrote.
				for _, o := range c.Outputs {
					_ = os.Remove(o)
				}
				return replay(t, tool, c, r)
			} else if !recording {
				t.Errorf("live %s reading %v of %s differs from the committed recording %s: %s\n(re-record with %s)",
					tool.Name, r.Args, inputNames(r), bundlePath(t), d, recordHint(t))
			}
		}
	}
	if recording {
		t.Cleanup(func() {
			if t.Failed() {
				t.Logf("%s: test failed, %s reading not recorded", bundlePath(t), tool.Name)
				return
			}
			saveRecord(t, r)
		})
	}
	// Outputs are already on disk; Run writes them back unchanged.
	return r
}

func inputNames(r *record) string {
	var n []string
	for _, f := range r.Inputs {
		n = append(n, f.Name+" sha256 "+f.SHA256[:12])
	}
	return strings.Join(n, ", ")
}

var (
	saveMu   sync.Mutex
	reRecord = map[string]bool{} // tests whose old readings this process dropped
)

// saveRecord adds a reading to the current test's bundle, first dropping
// every reading the same (sub)test recorded before, so a re-recorded test
// leaves no stale readings behind. Readings are recorded only by a test
// that passed: the tools verified the bytes first.
func saveRecord(t testing.TB, r *record) {
	t.Helper()
	saveMu.Lock()
	defer saveMu.Unlock()
	path := bundlePath(t)
	var recs []*record
	if _, err := os.Stat(path); err == nil {
		old, err := readBundle(path)
		if err != nil {
			t.Fatal(err)
		}
		recs = old
	}
	oldPlain, _ := encodeBundle(append([]*record(nil), recs...))
	first := !reRecord[r.Test]
	reRecord[r.Test] = true
	kept := recs[:0]
	for _, o := range recs {
		if (first && o.Test == r.Test) || (o.Test == r.Test && o.Key == r.Key) {
			continue
		}
		kept = append(kept, o)
	}
	recs = append(kept, r)
	plain, gz := encodeBundle(recs)
	if bytes.Equal(plain, oldPlain) {
		return
	}
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		t.Fatal(err)
	}
	tmp := path + ".tmp"
	if err := os.WriteFile(tmp, gz, 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.Rename(tmp, path); err != nil {
		t.Fatal(err)
	}
}

// toolCommand resolves the program of a tool, failing (never skipping)
// when it is not available: the test asked for tier 2.
func toolCommand(t testing.TB, tool Tool) (string, []string) {
	t.Helper()
	switch {
	case tool.Harness != "":
		return Harness(t, tool.Harness), nil
	case tool.Program != "":
		return RequireTool(t, tool.Program), nil
	case tool.Python:
		return PythonWith(t, tool.Modules...), nil
	}
	t.Fatalf("tool %s names no program", tool.Name)
	return "", nil
}

// toolVersions is informational (not compared): the versions of the
// Python modules a Python tool uses.
func toolVersions(t testing.TB, tool Tool) string {
	if !tool.Python {
		return ""
	}
	return pythonVersions(t)
}
