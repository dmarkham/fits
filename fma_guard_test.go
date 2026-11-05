package fits

import (
	"go/ast"
	"go/build"
	"go/importer"
	"go/parser"
	"go/token"
	"go/types"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestNoFusableFloatProducts keeps the library's float arithmetic free of
// multiply-adds the Go compiler may fuse into an FMA instruction (the Go
// spec allows it; amd64 at GOAMD64=v3 and arm64 do it), which would round
// differently from the C and Python truth sources and break the bit-exact
// claims on those builds. Every non-constant float32/float64 product that
// does not feed another product or a division, and is not already the
// operand of a conversion, must be wrapped in an explicit conversion
// (float64(x*y) + z): an explicit conversion rounds and prevents fusion.
// The one exception is the frozen stats.MeanStdev, whose result is pinned
// per build instead (TestMeanStdev_OnePassFormula). The whole suite is
// also run on GOAMD64=v3 and arm64 (docs/VALIDATION.md).
func TestNoFusableFloatProducts(t *testing.T) {
	allowed := map[string]bool{"stats/stats.go:MeanStdev": true}
	root := testutil.RepoRoot(t)
	fset := token.NewFileSet()
	l := &fmaLoader{root: root, fset: fset, std: importer.ForCompiler(fset, "source", nil),
		pkgs: map[string]*types.Package{}, files: map[string][]*ast.File{}, info: map[string]*types.Info{}}
	checked := 0
	for _, rel := range []string{".", "compress", "header", "healpix", "stats", "stretch", "wcs", "wcs/transform",
		"internal/bigendian", "internal/bitpix", "internal/block", "internal/checksum", "internal/tform"} {
		path := fmaModPath
		if rel != "." {
			path += "/" + rel
		}
		if _, err := l.Import(path); err != nil {
			t.Fatalf("type-checking %s: %v", path, err)
		}
		info := l.info[path]
		for _, f := range l.files[path] {
			relname, _ := filepath.Rel(root, fset.Position(f.Pos()).Filename)
			for _, d := range f.Decls {
				if fd, ok := d.(*ast.FuncDecl); ok && allowed[relname+":"+fd.Name.Name] {
					continue
				}
				var stack []ast.Node
				ast.Inspect(d, func(n ast.Node) bool {
					if n == nil {
						stack = stack[:len(stack)-1]
						return true
					}
					defer func() { stack = append(stack, n) }()
					be, ok := n.(*ast.BinaryExpr)
					if !ok || be.Op != token.MUL {
						return true
					}
					tv, ok := info.Types[be]
					if !ok || tv.Value != nil {
						return true
					}
					b, ok := tv.Type.Underlying().(*types.Basic)
					if !ok || (b.Kind() != types.Float64 && b.Kind() != types.Float32) {
						return true
					}
					checked++
					var parent ast.Node
					for i := len(stack) - 1; i >= 0; i-- {
						if _, isParen := stack[i].(*ast.ParenExpr); !isParen {
							parent = stack[i]
							break
						}
					}
					switch p := parent.(type) {
					case *ast.BinaryExpr:
						if p.Op == token.MUL || p.Op == token.QUO {
							return true
						}
					case *ast.CallExpr:
						if ptv, ok := info.Types[p.Fun]; ok && ptv.IsType() {
							return true
						}
					}
					t.Errorf("%s: float product not wrapped in an explicit conversion (the compiler may fuse it into an FMA)", fset.Position(be.Pos()))
					return true
				})
			}
		}
	}
	if checked < 500 {
		t.Fatalf("only %d float products found; the walk is broken", checked)
	}
}

const fmaModPath = "github.com/dmarkham/fits"

type fmaLoader struct {
	root  string
	fset  *token.FileSet
	std   types.Importer
	pkgs  map[string]*types.Package
	files map[string][]*ast.File
	info  map[string]*types.Info
}

func (l *fmaLoader) Import(path string) (*types.Package, error) {
	if p, ok := l.pkgs[path]; ok {
		return p, nil
	}
	if path == fmaModPath || strings.HasPrefix(path, fmaModPath+"/") {
		rel := strings.TrimPrefix(strings.TrimPrefix(path, fmaModPath), "/")
		return l.load(path, filepath.Join(l.root, rel))
	}
	return l.std.Import(path)
}

func (l *fmaLoader) load(path, dir string) (*types.Package, error) {
	ents, err := os.ReadDir(dir)
	if err != nil {
		return nil, err
	}
	var files []*ast.File
	ctx := build.Default
	for _, e := range ents {
		n := e.Name()
		if !strings.HasSuffix(n, ".go") || strings.HasSuffix(n, "_test.go") {
			continue
		}
		if ok, _ := ctx.MatchFile(dir, n); !ok {
			continue
		}
		f, err := parser.ParseFile(l.fset, filepath.Join(dir, n), nil, parser.ParseComments)
		if err != nil {
			return nil, err
		}
		files = append(files, f)
	}
	info := &types.Info{Types: map[ast.Expr]types.TypeAndValue{}}
	conf := types.Config{Importer: l}
	p, err := conf.Check(path, l.fset, files, info)
	if err != nil {
		return nil, err
	}
	l.pkgs[path], l.files[path], l.info[path] = p, files, info
	return p, nil
}
