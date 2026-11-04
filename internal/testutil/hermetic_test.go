package testutil

import (
	"go/ast"
	"go/build"
	"go/parser"
	"go/token"
	"io/fs"
	"path/filepath"
	"strings"
	"testing"
)

// The tier rule (docs/VALIDATION.md, "Test tiers"): plain `go test ./...`
// needs nothing but Go and the committed files. No Go file that tier 1
// compiles -- on any of the platforms below, without the interop tag --
// may start another program, except the test binary itself (os.Args[0],
// for cross-process tests). Starting a reference tool belongs in a file
// with //go:build interop.

var tier1Platforms = [][2]string{
	{"linux", "amd64"}, {"linux", "arm64"},
	{"darwin", "amd64"}, {"darwin", "arm64"},
	{"windows", "amd64"}, {"windows", "arm64"},
}

// moduleGoFiles lists every .go file of the module that the go tool can
// compile (testdata, reference trees and the venv are not packages).
func moduleGoFiles(t *testing.T) []string {
	t.Helper()
	root := RepoRoot(t)
	var files []string
	err := filepath.WalkDir(root, func(p string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if d.IsDir() {
			n := d.Name()
			if p != root && (n == "testdata" || n == "reference" || n == "venv" || n == "node_modules" ||
				strings.HasPrefix(n, ".") || strings.HasPrefix(n, "_")) {
				return filepath.SkipDir
			}
			return nil
		}
		if strings.HasSuffix(p, ".go") {
			files = append(files, p)
		}
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(files) < 100 {
		t.Fatalf("found only %d Go files under %s", len(files), root)
	}
	return files
}

// inTier1 reports whether a plain `go test` builds the file on some
// tier 1 platform (build tags and _GOOS/_GOARCH file names; the interop
// tag is not set).
func inTier1(t *testing.T, path string) bool {
	t.Helper()
	dir, name := filepath.Split(path)
	for _, pl := range tier1Platforms {
		ctx := build.Default
		ctx.GOOS, ctx.GOARCH = pl[0], pl[1]
		ctx.BuildTags = nil
		ctx.CgoEnabled = false
		ok, err := ctx.MatchFile(dir, name)
		if err != nil {
			t.Fatalf("%s: %v", path, err)
		}
		if ok {
			return true
		}
	}
	return false
}

// selfExec reports whether call starts the test binary itself:
// exec.Command(os.Args[0], ...) or exec.CommandContext(ctx, os.Args[0], ...).
func selfExec(call *ast.CallExpr, fn string) bool {
	i := 0
	if fn == "CommandContext" {
		i = 1
	}
	if len(call.Args) <= i {
		return false
	}
	ix, ok := call.Args[i].(*ast.IndexExpr)
	if !ok {
		return false
	}
	sel, ok := ix.X.(*ast.SelectorExpr)
	if !ok || sel.Sel.Name != "Args" {
		return false
	}
	pkg, ok := sel.X.(*ast.Ident)
	lit, ok2 := ix.Index.(*ast.BasicLit)
	return ok && ok2 && pkg.Name == "os" && lit.Value == "0"
}

// tier1Violations lists the ways one file starts a program.
func tier1Violations(fset *token.FileSet, f *ast.File) []string {
	execNames := map[string]bool{}
	for _, im := range f.Imports {
		if im.Path.Value == `"os/exec"` {
			if im.Name != nil {
				execNames[im.Name.Name] = true
			} else {
				execNames["exec"] = true
			}
		}
	}
	var bad []string
	report := func(n ast.Node, what string) {
		bad = append(bad, fset.Position(n.Pos()).String()+": "+what)
	}
	allowed := map[*ast.SelectorExpr]bool{}
	ast.Inspect(f, func(n ast.Node) bool {
		if call, ok := n.(*ast.CallExpr); ok {
			if sel, ok := call.Fun.(*ast.SelectorExpr); ok {
				if id, ok := sel.X.(*ast.Ident); ok && execNames[id.Name] &&
					(sel.Sel.Name == "Command" || sel.Sel.Name == "CommandContext") && selfExec(call, sel.Sel.Name) {
					allowed[sel] = true
				}
			}
		}
		sel, ok := n.(*ast.SelectorExpr)
		if !ok {
			return true
		}
		id, ok := sel.X.(*ast.Ident)
		if !ok {
			return true
		}
		switch {
		case execNames[id.Name]:
			switch sel.Sel.Name {
			case "Cmd", "ExitError":
				// types used to drive the self-exec child
			default:
				if !allowed[sel] {
					report(sel, "os/exec."+sel.Sel.Name+" (only exec.Command(os.Args[0], ...) is allowed in tier 1)")
				}
			}
		case id.Name == "os" && sel.Sel.Name == "StartProcess":
			report(sel, "os.StartProcess")
		case id.Name == "syscall" && (sel.Sel.Name == "Exec" || sel.Sel.Name == "ForkExec" || sel.Sel.Name == "StartProcess"):
			report(sel, "syscall."+sel.Sel.Name)
		}
		return true
	})
	return bad
}

// TestTier1StartsNoPrograms scans every Go file tier 1 compiles.
func TestTier1StartsNoPrograms(t *testing.T) {
	fset := token.NewFileSet()
	tier1, tier2 := 0, 0
	for _, p := range moduleGoFiles(t) {
		f, err := parser.ParseFile(fset, p, nil, parser.SkipObjectResolution)
		if err != nil {
			t.Fatal(err)
		}
		rel, _ := filepath.Rel(RepoRoot(t), p)
		if !inTier1(t, p) {
			tier2++
			continue
		}
		tier1++
		if strings.HasSuffix(p, "_interop_test.go") {
			t.Errorf("%s: an *_interop_test.go file must carry //go:build interop", rel)
		}
		for _, v := range tier1Violations(fset, f) {
			t.Errorf("tier 1 (plain go test) starts a program: %s; move it to a //go:build interop file", v)
		}
	}
	t.Logf("%d Go files in tier 1, %d only in other builds (interop or other platforms)", tier1, tier2)
}

// The scanner itself: it must catch what it is for.
func TestTier1ScannerCatches(t *testing.T) {
	src := `package x
import (
	"os"
	"os/exec"
	run "os/exec"
	"syscall"
)
func f() {
	exec.Command(os.Args[0], "-test.run=X")
	var c *exec.Cmd
	_ = c
	exec.Command("fitsverify", "x")
	exec.LookPath("cc")
	run.Command("python3")
	os.StartProcess("/bin/sh", nil, nil)
	syscall.Exec("/bin/sh", nil, nil)
}`
	fset := token.NewFileSet()
	f, err := parser.ParseFile(fset, "x.go", src, 0)
	if err != nil {
		t.Fatal(err)
	}
	bad := tier1Violations(fset, f)
	want := []string{"exec.Command", "exec.LookPath", "exec.Command", "os.StartProcess", "syscall.Exec"}
	if len(bad) != len(want) {
		t.Fatalf("violations %d, want %d:\n%s", len(bad), len(want), strings.Join(bad, "\n"))
	}
	for i, w := range want {
		if !strings.Contains(bad[i], strings.SplitN(w, ".", 2)[1]) {
			t.Errorf("violation %d = %q, want %s", i, bad[i], w)
		}
	}
}
