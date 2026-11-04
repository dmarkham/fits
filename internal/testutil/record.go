package testutil

import (
	"bytes"
	"compress/gzip"
	"crypto/sha256"
	"encoding/base64"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strings"
	"sync"
	"testing"
	"unicode/utf8"
)

// Recorded tool readings (the write golden chain, docs/VALIDATION.md).
//
// A test that wants an external program's reading of a file (typically a
// file the library just wrote) calls Run. The reading is identified by a
// key over everything that determines it: the tool, the committed source
// of the harness (Tool.Sources), the arguments and standard input with
// the file paths replaced by placeholders, and the base name and sha256
// of every input file. The reading itself (exit status, stdout, stderr,
// the files the program wrote) is committed in the bundle of the test
// that recorded it, testdata/recorded/<package>/<TestName>.jsonl.gz.
//
//   - Tier 1 (no build tag) never starts the program: Run replays the
//     committed reading. A missing reading fails the test: the library
//     no longer writes the bytes the tools verified (or the case is new),
//     so tier 1 requires the writer to reproduce those bytes exactly. The
//     test's own assertions then run on the replayed reading, so every
//     proof still holds.
//   - Tier 2 (-tags interop, record_live.go) runs the program live, hands
//     the live reading to the same assertions and, on amd64, requires it
//     to equal the committed one. With FITS_RECORD=1 a test that passes
//     writes (or rewrites) its readings: the tools verify the bytes first,
//     the recording follows.

// EnvRecord, set to "1" in a tier 2 run, records the live readings of
// every test that passes into testdata/recorded.
const EnvRecord = "FITS_RECORD"

// Tool names one external program and how tier 2 finds it.
type Tool struct {
	// Name is the directory under testdata/recorded.
	Name string
	// Harness is the repository-relative path of a built reference
	// harness (built by `make -C <its dir>`); Program a program found on
	// PATH; Python the FITS_ASTROPY_PYTHON interpreter (with Modules
	// importable). Exactly one is set.
	Harness string
	Program string
	Python  bool
	Modules []string
	// Sources are repository-relative files (the harness source) whose
	// bytes are part of every key, so a changed harness invalidates its
	// readings.
	Sources []string
}

// The reference tools shared by several packages.
var (
	ToolCfref      = Tool{Name: "cfref", Harness: "testdata/cfref/cfref", Sources: []string{"testdata/cfref/cfref.c"}}
	ToolFitsverify = Tool{Name: "fitsverify", Program: "fitsverify"}
	ToolFunpack    = Tool{Name: "funpack", Program: "funpack"}
	ToolFpack      = Tool{Name: "fpack", Program: "fpack"}
	ToolAstropy    = Tool{Name: "astropy", Python: true}
)

// Call is one run of a Tool.
type Call struct {
	// Args are the program's arguments. Any occurrence of an Inputs,
	// Outputs or Dir path is replaced by a placeholder in the key and the
	// recording, and substituted back on replay.
	Args []string
	// Inputs are the files the program reads; their bytes key the reading.
	Inputs []string
	// Outputs are files the program writes; their bytes are recorded and
	// written back on replay.
	Outputs []string
	// OutputsByDigest records only the sha256 of each output: for an
	// output that is handed to further tool calls only (a funpacked
	// image read by the cfitsio harness), never read by the test itself.
	// On replay the file holds a stub naming that digest, which the next
	// call keys on as if it were the real bytes.
	OutputsByDigest bool
	// FixOutput, when set, rewrites each output's bytes (tier 2, right
	// after the program ran, before they are digested): for a program
	// that stamps the wall-clock time into its output, so the bytes, and
	// the readings keyed on them, are the same on every run.
	FixOutput func(b []byte) []byte
	// Compact, when set, is applied to the live standard output before it
	// is recorded (and compared with the recording): a reading too large
	// to commit is stored as an exact digest plus a few explicit values,
	// which the caller resolves on replay against a candidate it computes
	// (see the callers). The caller gets the full live output in tier 2
	// and the compact form in tier 1, and must accept both.
	Compact func(stdout []byte) []byte
	// Stdin is the program's standard input (paths replaced as in Args).
	Stdin []byte
	// Dir is the working directory (tier 2), "" for the test's.
	Dir string
}

// Result is a program's reading.
type Result struct {
	Exit   int // exit status; non-zero is for the caller to judge
	Stdout []byte
	Stderr []byte
}

// record is the committed form of one reading.
type record struct {
	Tool        string       `json:"tool"`
	Key         string       `json:"key"`
	Sources     []recordFile `json:"sources,omitempty"`
	Args        []string     `json:"args"`
	Inputs      []recordFile `json:"inputs,omitempty"`
	Stdin       string       `json:"stdin_sha256,omitempty"`
	Exit        int          `json:"exit"`
	Stdout      string       `json:"stdout,omitempty"`
	StdoutLines []int        `json:"stdout_lines,omitempty"` // in a bundle: stdout as indices into its line pool
	StdoutB64   string       `json:"stdout_base64,omitempty"`
	Stderr      string       `json:"stderr,omitempty"`
	StderrB64   string       `json:"stderr_base64,omitempty"`
	Outputs     []recordFile `json:"outputs,omitempty"`
	RecordedOn  string       `json:"recorded_on"`
	Versions    string       `json:"versions,omitempty"`
	Test        string       `json:"test,omitempty"` // the test that recorded it (informational)
	stdoutRaw   []byte       // decoded stdout (not serialized)
	liveStdout  []byte       // tier 2: the full live stdout when Compact shrank it
}

// recordFile is an input (name, digest) or an output (with its bytes).
type recordFile struct {
	Name   string `json:"name"`
	SHA256 string `json:"sha256"`
	Size   int    `json:"size"`
	Data   string `json:"data_base64,omitempty"`
}

// placeholders maps each path of a call to its placeholder, longest
// path first so a path is never replaced inside a longer one.
type placeholders [][2]string

func (c *Call) placeholders() placeholders {
	var p placeholders
	for i, f := range c.Inputs {
		p = append(p, [2]string{f, fmt.Sprintf("{in%d}", i)})
	}
	for i, f := range c.Outputs {
		p = append(p, [2]string{f, fmt.Sprintf("{out%d}", i)})
	}
	if c.Dir != "" {
		p = append(p, [2]string{c.Dir, "{dir}"})
	}
	sort.SliceStable(p, func(i, j int) bool { return len(p[i][0]) > len(p[j][0]) })
	return p
}

func (p placeholders) hide(s string) string {
	for _, r := range p {
		if r[0] != "" {
			s = strings.ReplaceAll(s, r[0], r[1])
		}
	}
	if strings.Contains(s, "\n") {
		// cfitsio wraps its error messages at a fixed width, line breaks
		// falling inside a path; such a path is a placeholder too.
		for _, r := range p {
			if len(r[0]) < 2 {
				continue
			}
			var re strings.Builder
			for i, c := range r[0] {
				if i > 0 {
					re.WriteString("\n?")
				}
				re.WriteString(regexp.QuoteMeta(string(c)))
			}
			s = regexp.MustCompile(re.String()).ReplaceAllLiteralString(s, r[1])
		}
	}
	return s
}

func (p placeholders) show(s string) string {
	for i := len(p) - 1; i >= 0; i-- {
		s = strings.ReplaceAll(s, p[i][1], p[i][0])
	}
	return s
}

func sha(b []byte) string {
	s := sha256.Sum256(b)
	return hex.EncodeToString(s[:])
}

// stubMagic starts a replayed OutputsByDigest file: the stub stands for
// bytes with the digest and size that follow it.
const stubMagic = "FITS-TESTUTIL-RECORDED-OUTPUT\x00"

// fileDigest reads a file and returns its record entry.
func fileDigest(t testing.TB, path string, keep bool) recordFile {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("recorded tool call: %v", err)
	}
	if rest, ok := strings.CutPrefix(string(b), stubMagic); ok {
		var f recordFile
		if err := json.Unmarshal([]byte(rest), &f); err != nil {
			t.Fatalf("recorded output stub %s: %v", path, err)
		}
		f.Name = filepath.Base(path)
		return f
	}
	f := recordFile{Name: filepath.Base(path), SHA256: sha(b), Size: len(b)}
	if keep {
		f.Data = base64.StdEncoding.EncodeToString(b)
	}
	return f
}

// prepare computes the key of a call and the parts of its record that
// do not depend on running it.
func prepare(t testing.TB, tool Tool, c *Call) (*record, placeholders) {
	t.Helper()
	p := c.placeholders()
	r := &record{Tool: tool.Name}
	root := RepoRoot(t)
	for _, s := range tool.Sources {
		f := fileDigest(t, filepath.Join(root, s), false)
		f.Name = s
		r.Sources = append(r.Sources, f)
	}
	tmp := os.TempDir()
	for _, a := range c.Args {
		h := p.hide(a)
		if strings.Contains(h, tmp) {
			// A temp path in the key makes it differ on every run.
			t.Fatalf("recorded tool call: argument %q names a temporary path that is not one of the call's Inputs, Outputs or Dir", a)
		}
		r.Args = append(r.Args, h)
	}
	for _, in := range c.Inputs {
		r.Inputs = append(r.Inputs, fileDigest(t, in, false))
	}
	if c.Stdin != nil {
		h := p.hide(string(c.Stdin))
		if strings.Contains(h, tmp) {
			t.Fatalf("recorded tool call: standard input names a temporary path that is not one of the call's Inputs, Outputs or Dir")
		}
		r.Stdin = sha([]byte(h))
	}
	type keyed struct {
		Tool    string       `json:"tool"`
		Sources []recordFile `json:"sources"`
		Args    []string     `json:"args"`
		Inputs  []recordFile `json:"inputs"`
		Stdin   string       `json:"stdin"`
		Outputs []string     `json:"outputs"`
	}
	k := keyed{Tool: r.Tool, Sources: r.Sources, Args: r.Args, Inputs: r.Inputs, Stdin: r.Stdin}
	for _, o := range c.Outputs {
		k.Outputs = append(k.Outputs, filepath.Base(o))
	}
	kb, _ := json.Marshal(k)
	r.Key = sha(kb)
	return r, p
}

// Readings are committed per package and top-level test, as one gzipped
// JSON-lines bundle: testdata/recorded/<package>/<TestName>.jsonl.gz
// (<package> is the package directory relative to the module root, "_"
// for the root package). One bundle per test keeps the many similar
// readings of a test matrix together, so they compress well.

// packageDir is the test's package directory relative to the module
// root, "_" for the root.
func packageDir(t testing.TB) string {
	t.Helper()
	wd, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	rel, err := filepath.Rel(RepoRoot(t), wd)
	if err != nil || rel == "." {
		return "_"
	}
	return filepath.ToSlash(rel)
}

// bundlePath is the bundle of the current top-level test.
func bundlePath(t testing.TB) string {
	t.Helper()
	top, _, _ := strings.Cut(t.Name(), "/")
	return TestdataPath(t, "recorded", filepath.FromSlash(packageDir(t)), top+".jsonl.gz")
}

var (
	bundleMu    sync.Mutex
	bundleIndex map[string]*record // key -> reading, every bundle of this package
	bundleErr   error
)

// readRecord returns the committed reading for r's key, nil if none.
func readRecord(t testing.TB, r *record) *record {
	t.Helper()
	bundleMu.Lock()
	defer bundleMu.Unlock()
	if bundleIndex == nil {
		bundleIndex, bundleErr = loadBundles(TestdataPath(t, "recorded", filepath.FromSlash(packageDir(t))))
	}
	if bundleErr != nil {
		t.Fatalf("recorded readings: %v", bundleErr)
	}
	got, ok := bundleIndex[r.Key]
	if !ok {
		return nil
	}
	if got.Tool != r.Tool {
		t.Fatalf("recorded reading %s is a %s reading, want %s", r.Key, got.Tool, r.Tool)
	}
	return got
}

// loadBundles indexes every bundle in dir by key.
func loadBundles(dir string) (map[string]*record, error) {
	idx := map[string]*record{}
	files, _ := filepath.Glob(filepath.Join(dir, "*.jsonl.gz"))
	for _, f := range files {
		recs, err := readBundle(f)
		if err != nil {
			return nil, err
		}
		for _, r := range recs {
			idx[r.Key] = r
		}
	}
	return idx, nil
}

// readBundle decodes one bundle (strictly: unknown fields are errors).
func readBundle(path string) ([]*record, error) {
	z, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	b, err := Gunzip(z)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	var pool bundlePool
	if err := dec.Decode(&pool); err != nil {
		return nil, fmt.Errorf("%s: line pool: %w", path, err)
	}
	var out []*record
	for dec.More() {
		var r record
		if err := dec.Decode(&r); err != nil {
			return nil, fmt.Errorf("%s: %w", path, err)
		}
		if r.StdoutLines != nil {
			lines := make([]string, len(r.StdoutLines))
			for i, k := range r.StdoutLines {
				if k < 0 || k >= len(pool.Lines) {
					return nil, fmt.Errorf("%s: %s: line %d not in the pool", path, r.Key, k)
				}
				lines[i] = pool.Lines[k]
			}
			r.Stdout, r.StdoutLines = strings.Join(lines, "\n"), nil
		}
		if r.StderrB64 != "" {
			raw, err := base64.StdEncoding.DecodeString(r.StderrB64)
			if err != nil {
				return nil, fmt.Errorf("%s: %s: stderr: %w", path, r.Key, err)
			}
			r.Stderr, r.StderrB64 = string(raw), ""
		}
		if r.StdoutB64 != "" {
			raw, err := base64.StdEncoding.DecodeString(r.StdoutB64)
			if err != nil {
				return nil, fmt.Errorf("%s: %s: stdout: %w", path, r.Key, err)
			}
			r.stdoutRaw = raw
		} else {
			r.stdoutRaw = []byte(r.Stdout)
		}
		out = append(out, &r)
	}
	return out, nil
}

// setStdout stores stdout as text when it is valid UTF-8 (readable
// diffs), base64 otherwise.
func (r *record) setStdout(b []byte) {
	r.stdoutRaw = b
	r.Stdout, r.StdoutB64 = "", ""
	if utf8.Valid(b) {
		r.Stdout = string(b)
	} else {
		r.StdoutB64 = base64.StdEncoding.EncodeToString(b)
	}
}

// bundlePool is the first JSON value of a bundle: every distinct line of
// the text standard outputs, each stored once (a test matrix dumps many
// nearly identical files, too large for gzip's window to share).
type bundlePool struct {
	Lines []string `json:"lines"`
}

// encodeBundle is the committed bytes of a bundle: the line pool, then
// one compact JSON record per line, sorted by test and key, gzipped.
func encodeBundle(recs []*record) (plain, gz []byte) {
	sort.Slice(recs, func(i, j int) bool {
		if recs[i].Test != recs[j].Test {
			return recs[i].Test < recs[j].Test
		}
		return recs[i].Key < recs[j].Key
	})
	var pool bundlePool
	ids := map[string]int{}
	pooled := make([]record, len(recs))
	for i, r := range recs {
		pooled[i] = *r
		if !utf8.ValidString(r.Stderr) {
			pooled[i].Stderr, pooled[i].StderrB64 = "", base64.StdEncoding.EncodeToString([]byte(r.Stderr))
		}
		if r.StdoutB64 != "" || r.Stdout == "" {
			continue
		}
		lines := strings.Split(r.Stdout, "\n")
		pooled[i].Stdout, pooled[i].StdoutLines = "", make([]int, len(lines))
		for j, l := range lines {
			k, ok := ids[l]
			if !ok {
				k = len(pool.Lines)
				ids[l] = k
				pool.Lines = append(pool.Lines, l)
			}
			pooled[i].StdoutLines[j] = k
		}
	}
	var buf bytes.Buffer
	enc := json.NewEncoder(&buf)
	enc.SetEscapeHTML(false)
	_ = enc.Encode(pool)
	for i := range pooled {
		_ = enc.Encode(&pooled[i])
	}
	var z bytes.Buffer
	w, _ := gzip.NewWriterLevel(&z, gzip.BestCompression)
	_, _ = w.Write(buf.Bytes())
	_ = w.Close()
	return buf.Bytes(), z.Bytes()
}

// sameReading reports how two readings of the same key differ, "" if not.
func sameReading(a, b *record) string {
	switch {
	case a.Exit != b.Exit:
		return fmt.Sprintf("exit %d, recorded %d", a.Exit, b.Exit)
	case !bytes.Equal(a.stdoutRaw, b.stdoutRaw):
		return fmt.Sprintf("stdout differs (%d vs %d bytes) at byte %d", len(a.stdoutRaw), len(b.stdoutRaw), firstByteDiff(a.stdoutRaw, b.stdoutRaw))
	case a.Stderr != b.Stderr:
		return fmt.Sprintf("stderr %q, recorded %q", a.Stderr, b.Stderr)
	case len(a.Outputs) != len(b.Outputs):
		return fmt.Sprintf("%d output files, recorded %d", len(a.Outputs), len(b.Outputs))
	}
	for i := range a.Outputs {
		if a.Outputs[i].SHA256 != b.Outputs[i].SHA256 {
			return fmt.Sprintf("output %s differs", a.Outputs[i].Name)
		}
	}
	return ""
}

func firstByteDiff(a, b []byte) int {
	for i := 0; i < len(a) && i < len(b); i++ {
		if a[i] != b[i] {
			return i
		}
	}
	if len(a) < len(b) {
		return len(a)
	}
	return len(b)
}

// recordHint is how to (re)record the current test's readings.
func recordHint(t testing.TB) string {
	name := t.Name()
	if i := strings.IndexByte(name, '/'); i >= 0 {
		name = name[:i]
	}
	pkg := "."
	if wd, err := os.Getwd(); err == nil {
		if rel, err := filepath.Rel(RepoRoot(t), wd); err == nil && rel != "." {
			pkg = "./" + filepath.ToSlash(rel)
		}
	}
	return fmt.Sprintf("%s=1 go test -tags interop -count=1 -run '^%s$' %s", EnvRecord, regexp.QuoteMeta(name), pkg)
}

// Run returns tool's reading for call: replayed from testdata/recorded in
// tier 1, run live (and checked against the recording) in tier 2.
func Run(t testing.TB, tool Tool, c Call) Result {
	t.Helper()
	r, p := prepare(t, tool, &c)
	got := execute(t, tool, &c, r, p)
	for i, o := range got.Outputs {
		if o.Data == "" {
			// OutputsByDigest: live, the real file is on disk; replayed,
			// a stub stands for it.
			if _, err := os.Stat(c.Outputs[i]); err == nil {
				continue
			}
			stub, _ := json.Marshal(recordFile{SHA256: o.SHA256, Size: o.Size})
			if err := os.WriteFile(c.Outputs[i], append([]byte(stubMagic), stub...), 0o644); err != nil {
				t.Fatal(err)
			}
			continue
		}
		b, err := base64.StdEncoding.DecodeString(o.Data)
		if err != nil {
			t.Fatalf("recorded %s output %s: %v", tool.Name, o.Name, err)
		}
		if err := os.WriteFile(c.Outputs[i], b, 0o644); err != nil {
			t.Fatal(err)
		}
	}
	out := got.stdoutRaw
	if got.liveStdout != nil {
		out = got.liveStdout
	}
	return Result{Exit: got.Exit, Stdout: []byte(p.show(string(out))), Stderr: []byte(p.show(got.Stderr))}
}

// replay is the tier 1 reading: the committed record, or a failure that
// says the bytes changed.
func replay(t testing.TB, tool Tool, c *Call, r *record) *record {
	t.Helper()
	got := readRecord(t, r)
	if got == nil {
		var in []string
		for i, f := range r.Inputs {
			in = append(in, fmt.Sprintf("%s (%s: %d bytes, sha256 %s)", f.Name, c.Inputs[i], f.Size, f.SHA256))
		}
		t.Fatalf("no recorded %s reading %v for %s:\n"+
			"the file the library wrote is not the one the reference tools verified (tier 1 requires the exact bytes), "+
			"or the case is new. Verify it live and record the reading with:\n  %s",
			tool.Name, r.Args, strings.Join(in, ", "), recordHint(t))
	}
	if len(got.Outputs) != len(c.Outputs) {
		t.Fatalf("recorded %s reading %s has %d outputs, the call names %d", tool.Name, r.Key, len(got.Outputs), len(c.Outputs))
	}
	return got
}

// RunOK is Run for a call that must succeed: a non-zero exit fails the
// test with the program's stderr.
func RunOK(t testing.TB, tool Tool, c Call) Result {
	t.Helper()
	r := Run(t, tool, c)
	if r.Exit != 0 {
		t.Fatalf("%s %v: exit status %d\n%s%s", tool.Name, c.Args, r.Exit, r.Stderr, r.Stdout)
	}
	return r
}

// AstropyScript is the standard output of the committed Python script
// (run with FITS_ASTROPY_PYTHON in tier 2) on files, which are its
// arguments; the script's bytes key the reading too. A failing script
// fails the test.
func AstropyScript(t testing.TB, script string, files ...string) []byte {
	t.Helper()
	RequireFile(t, script)
	return RunOK(t, ToolAstropy, Call{Args: append([]string{script}, files...), Inputs: append([]string{script}, files...)}).Stdout
}

// AstropyCode is AstropyScript for inline Python code (python -c code
// files...).
func AstropyCode(t testing.TB, code string, files ...string) []byte {
	t.Helper()
	return RunOK(t, ToolAstropy, Call{Args: append([]string{"-c", code}, files...), Inputs: files}).Stdout
}
