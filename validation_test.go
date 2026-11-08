package fits

import (
	"errors"
	"io/fs"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// ledgerDocs are the documents whose references TestValidationLedger
// checks. docs/VALIDATION.md is the claim -> test -> truth ledger; the
// others are the top-level docs that cite tests and samples.
var ledgerDocs = []string{
	"docs/VALIDATION.md",
	"README.md",
	"plan.md",
	"testdata/README.md",
}

// ledgerPathPrefixes: a backticked token starting with one of these is a
// repository path (or glob) that must match at least one file.
var ledgerPathPrefixes = []string{
	"testdata/", "header/", "compress/", "wcs/", "healpix/", "stats/",
	"stretch/", "cmd/", "internal/", "docs/",
}

var (
	ledgerTokenRe = regexp.MustCompile("`([^`\n]+)`")
	ledgerTestRe  = regexp.MustCompile(`^((?:Test|Fuzz)[A-Z0-9_][A-Za-z0-9_]*)(/\S*)?$`)
	ledgerFuncRe  = regexp.MustCompile(`(?m)^func ((?:Test|Fuzz)[A-Za-z0-9_]*)\(`)
)

// repoTestFuncs returns the name of every Test/Fuzz function in the
// repository's _test.go files (venv, reference and .git excluded).
func repoTestFuncs(t *testing.T) map[string]string {
	t.Helper()
	funcs := map[string]string{}
	err := filepath.WalkDir(".", func(p string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if d.IsDir() {
			switch d.Name() {
			case ".git", "venv", "reference", "node_modules":
				return filepath.SkipDir
			}
			return nil
		}
		if !strings.HasSuffix(p, "_test.go") {
			return nil
		}
		b, err := os.ReadFile(p)
		if err != nil {
			return err
		}
		for _, m := range ledgerFuncRe.FindAllSubmatch(b, -1) {
			funcs[string(m[1])] = p
		}
		return nil
	})
	if err != nil {
		t.Fatalf("walking the repository: %v", err)
	}
	if len(funcs) < 400 {
		t.Fatalf("found only %d Test/Fuzz functions; the walk is broken", len(funcs))
	}
	return funcs
}

// TestValidationLedger keeps docs/VALIDATION.md (and the top-level docs)
// from rotting: every backticked test name must be a real Test/Fuzz
// function, every backticked repository path or glob must match a file,
// and every row of a ledger table must name at least one test.
func TestValidationLedger(t *testing.T) {
	funcs := repoTestFuncs(t)
	for _, doc := range ledgerDocs {
		b, err := os.ReadFile(doc)
		if err != nil {
			t.Fatalf("%s: %v", doc, err)
		}
		checkLedgerDoc(t, doc, string(b), funcs)
	}
}

func checkLedgerDoc(t *testing.T, doc, text string, funcs map[string]string) {
	t.Helper()
	isLedger := doc == "docs/VALIDATION.md"
	inFence, inLedgerSection := false, false
	rows, tests := 0, map[string]bool{}
	for i, line := range strings.Split(text, "\n") {
		lineNo := i + 1
		trim := strings.TrimSpace(line)
		if strings.HasPrefix(trim, "```") {
			inFence = !inFence
			continue
		}
		if inFence {
			continue
		}
		if strings.HasPrefix(line, "## ") {
			inLedgerSection = isLedger && strings.TrimSpace(line) == "## Ledger"
		}
		rowTests := 0
		for _, m := range ledgerTokenRe.FindAllStringSubmatch(line, -1) {
			tok := m[1]
			if tm := ledgerTestRe.FindStringSubmatch(tok); tm != nil {
				if _, ok := funcs[tm[1]]; !ok {
					t.Errorf("%s:%d: test %s does not exist", doc, lineNo, tm[1])
				}
				rowTests++
				tests[tm[1]] = true
				continue
			}
			if isLedgerPath(tok) {
				checkLedgerPath(t, doc, lineNo, tok)
			}
		}
		if inLedgerSection && strings.HasPrefix(trim, "|") && !isTableHeaderOrRule(trim) {
			rows++
			if rowTests == 0 {
				t.Errorf("%s:%d: ledger row names no proving test: %s", doc, lineNo, trim)
			}
		}
	}
	if inFence {
		t.Errorf("%s: unterminated code fence", doc)
	}
	if isLedger {
		// The ledger covers every area; a large drop means rows were lost.
		if rows < 100 || len(tests) < 300 {
			t.Errorf("%s: only %d ledger rows naming %d distinct tests (want >= 100 rows, >= 300 tests)", doc, rows, len(tests))
		}
	}
}

func isLedgerPath(tok string) bool {
	if strings.ContainsAny(tok, " <>") {
		return false
	}
	for _, p := range ledgerPathPrefixes {
		if strings.HasPrefix(tok, p) {
			return true
		}
	}
	return false
}

func checkLedgerPath(t *testing.T, doc string, lineNo int, tok string) {
	t.Helper()
	if strings.ContainsAny(tok, "*?[") {
		m, err := filepath.Glob(tok)
		if err != nil || len(m) == 0 {
			t.Errorf("%s:%d: glob %s matches no file", doc, lineNo, tok)
		}
		return
	}
	if !testutil.GoldenExists(strings.TrimSuffix(tok, "/")) {
		t.Errorf("%s:%d: path %s does not exist", doc, lineNo, tok)
	}
}

// isTableHeaderOrRule reports whether a markdown table line is the header
// row or the |---| separator of a ledger table.
func isTableHeaderOrRule(row string) bool {
	if strings.Trim(row, "|-: ") == "" {
		return true
	}
	cells := strings.Split(strings.Trim(row, "|"), "|")
	return len(cells) > 0 && strings.TrimSpace(cells[0]) == "Capability"
}

// Size budget for committed reference data (RULES, round 2): no tracked
// file over 1 MB except the FITS 4.0 standard, at most 20 MB of tracked
// files under testdata directories, at most 25 MB tracked in total.
// (1 MB = 1,000,000 bytes.) Large goldens are kept small by digests,
// compact encodings and gzip (testutil.ReadGolden), never by dropping
// checks.
const (
	sizeMaxFile     = 1_000_000
	sizeMaxTestdata = 20_000_000
	sizeMaxTotal    = 25_000_000
)

var sizeExempt = map[string]bool{}

// TestRepositorySizeBudget measures every file git tracks (its size in the
// working tree) against the budget above. The tracked list is read from
// the git index in pure Go (testutil.TrackedFiles; tier 2 checks it
// against `git ls-files`); outside a git checkout (a module download)
// every file of the tree is measured instead.
func TestRepositorySizeBudget(t *testing.T) {
	root := testutil.RepoRoot(t)
	files, ok, err := testutil.TrackedFiles(root)
	if err != nil {
		t.Fatalf("reading the git index of %s: %v", root, err)
	}
	if !ok {
		t.Logf("%s is not a git checkout: measuring every file of the tree", root)
		if files, err = testutil.WalkFiles(root); err != nil {
			t.Fatal(err)
		}
	}
	var total, testdata int64
	var n int
	for _, f := range files {
		st, err := os.Lstat(filepath.Join(root, f))
		if err != nil {
			if errors.Is(err, fs.ErrNotExist) {
				continue // deleted in the working tree, not yet committed
			}
			t.Fatal(err)
		}
		n++
		size := st.Size()
		total += size
		inTestdata := strings.HasPrefix(f, "testdata/") || strings.Contains(f, "/testdata/")
		if inTestdata {
			testdata += size
		}
		if size > sizeMaxFile && !sizeExempt[f] {
			t.Errorf("%s is %d bytes, over the %d-byte per-file budget (store it as a digest, compactly or gzipped: testutil.ReadGolden)", f, size, sizeMaxFile)
		}
	}
	if n < 1000 {
		t.Fatalf("only %d tracked files found", n)
	}
	if testdata > sizeMaxTestdata {
		t.Errorf("tracked testdata is %d bytes, over the %d-byte budget", testdata, sizeMaxTestdata)
	}
	if total > sizeMaxTotal {
		t.Errorf("tracked files total %d bytes, over the %d-byte budget", total, sizeMaxTotal)
	}
	t.Logf("%d tracked files: %d bytes in total, %d bytes under testdata", n, total, testdata)
}
