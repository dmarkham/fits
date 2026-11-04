//go:build interop

package transform

import (
	"bytes"
	"encoding/json"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// requireTNXTool is RequireTool with a hint that fits the TNX tools.
func requireTNXTool(t *testing.T, name string) string {
	t.Helper()
	p, err := exec.LookPath(name)
	if err != nil {
		t.Fatalf("%s not found on PATH (%v); the TNX generators need python3 (stdlib only) and WCSTools 3.9.7 xy2sky/sky2xy (Arch: wcstools package; see testdata/tnx/PROVENANCE.md)", name, err)
	}
	return p
}

// TestTNXGoldensReproduce re-runs gen_golden.py --check, which regenerates
// every golden with the installed WCSTools and fails if any byte differs.
// Needs python3 and xy2sky/sky2xy (tier 2: a failure without them).
func TestTNXGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := requireTNXTool(t, "python3")
	requireTNXTool(t, "xy2sky")
	requireTNXTool(t, "sky2xy")
	cmd := exec.Command(py, "gen_golden.py", "--check")
	cmd.Dir = filepath.Join("testdata", "tnx")
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("gen_golden.py --check: %v\n%s", err, out)
	}
}

// TestTNXSyntheticHeadersReproduce: make_synthetic.py (stdlib Python)
// must write exactly the syn_* cases of tnxRequiredCases and the
// refused/*.hdr of tnxRefusedCases, byte for byte.
// It runs `make_synthetic.py --check` on the committed directory (fails on
// a changed, missing or extra syn_*.hdr), then deletes headers/syn_*.hdr in
// a scratch copy (and refused/*.hdr) and requires the generator to write them back
// identically (gen_golden.py --check then proves the goldens from the
// headers).
func TestTNXSyntheticHeadersReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := requireTNXTool(t, "python3")
	src := filepath.Join("testdata", "tnx")
	testutil.RunIn(t, src, py, "make_synthetic.py", "--check")
	dir := t.TempDir()
	testutil.CopyTree(t, src, dir)
	syn, err := filepath.Glob(filepath.Join(dir, "headers", "syn_*.hdr"))
	if err != nil {
		t.Fatal(err)
	}
	got := map[string]bool{}
	for _, p := range syn {
		got[strings.TrimSuffix(filepath.Base(p), ".hdr")] = true
		if err := os.Remove(p); err != nil {
			t.Fatal(err)
		}
	}
	for n := range tnxRequiredCases {
		if strings.HasPrefix(n, "syn_") && !got[n] {
			t.Errorf("required synthetic header %s.hdr is not committed", n)
		}
	}
	refused, err := filepath.Glob(filepath.Join(dir, "refused", "*.hdr"))
	if err != nil {
		t.Fatal(err)
	}
	if len(refused) != len(tnxRefusedCases) {
		t.Errorf("%d refused/*.hdr committed, tnxRefusedCases has %d", len(refused), len(tnxRefusedCases))
	}
	for _, p := range refused {
		if err := os.Remove(p); err != nil {
			t.Fatal(err)
		}
	}
	testutil.RunIn(t, dir, py, "make_synthetic.py")
	testutil.SameFiles(t, dir, src, nil)
}

// TestCrossCheckFixturesReproduce re-runs the generator with the astropy
// interpreter (FITS_ASTROPY_PYTHON) and requires byte-identical output.
// Tier 2: a failure without astropy.
func TestCrossCheckFixturesReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.AstropyPython(t)
	script := testutil.TestdataPath(t, "generate_wcs_crosscheck.py")
	testutil.RequireFile(t, script)
	got, err := exec.Command(py, script).Output()
	if err != nil {
		t.Fatalf("generator failed: %v", err)
	}
	want, err := os.ReadFile(testutil.TestdataPath(t, "wcs_crosscheck.json"))
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(got, want) {
		var doc crossCheckDoc
		_ = json.Unmarshal(got, &doc)
		t.Fatalf("generator output differs from the committed testdata/wcs_crosscheck.json (this run: astropy %s, numpy %s); regenerate only after checking why", doc.Astropy, doc.Numpy)
	}
}

// TestFramesGoldenReproduce re-runs gen_frames.py and requires
// byte-identical output, so the committed golden is provably what the
// pinned astropy/ERFA produce. Needs FITS_ASTROPY_PYTHON (tier 2: a
// failure without it).
func TestFramesGoldenReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.AstropyPython(t)
	want, err := os.ReadFile(framesGolden)
	if err != nil {
		t.Fatalf("committed frames golden missing: %v", err)
	}
	cmd := exec.Command(py, filepath.Join("testdata", "frames", "gen_frames.py"))
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	got, err := cmd.Output()
	if err != nil {
		t.Fatalf("gen_frames.py: %v\n%s", err, stderr.String())
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("gen_frames.py output (%d bytes) differs from committed %s (%d bytes): wrong astropy/pyerfa version or a stale golden", len(got), framesGolden, len(want))
	}
}

// TestHdrGoldensReproduce rebuilds hdr2sky against the installed wcslib and
// requires it to regenerate every committed golden byte for byte. Needs a C
// compiler and libwcs (tier 2: a failure without them).
func TestHdrGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	unavailable := t.Fatalf
	if _, err := exec.LookPath("make"); err != nil {
		unavailable("make not found")
	}
	out, err := exec.Command("make", "-C", "testdata/wref", "hdr2sky").CombinedOutput()
	if err != nil {
		unavailable("cannot build hdr2sky (needs cc + wcslib %s): %v\n%s", hdrWcslibVersion, err, out)
	}
	bin := filepath.Join("testdata", "wref", "hdr2sky")
	for _, name := range hdrCaseNames(t) {
		got, err := exec.Command(bin, filepath.Join(hdrCaseDir, name+".case")).Output()
		if err != nil {
			t.Errorf("%s: hdr2sky failed: %v", name, err)
			continue
		}
		want, err := testutil.ReadGolden(filepath.Join(hdrGoldenDir, name+".json"))
		if err != nil {
			t.Errorf("%s: golden missing: %v", name, err)
			continue
		}
		if !bytes.Equal(got, want) {
			t.Errorf("%s: hdr2sky output differs from the committed golden (wrong wcslib version, or golden not regenerated)", name)
		}
	}
}

// TestHPXPaperGeneratorReproduces re-runs gen_hpx_paper.py and requires
// the committed table back (f64packSame) (tier 2: needs the venv's
// mpmath, a failure without it).
func TestHPXPaperGeneratorReproduces(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.AstropyPython(t)
	if b, err := exec.Command(py, "-c", "import mpmath").CombinedOutput(); err != nil {
		t.Fatalf("%s cannot import mpmath (%v: %s)", py, err, b)
	}
	out := t.TempDir()
	if b, err := exec.Command(py, filepath.Join(wrefDir, "gen_hpx_paper.py"), out).CombinedOutput(); err != nil {
		t.Fatalf("gen_hpx_paper.py: %v\n%s", err, b)
	}
	f64packSame(t, filepath.Join(out, "paper", "hpx_paper.json"), hpxPaperGolden)
}

// TestLibmGeneratorReproduces re-runs gen_libm.py and requires the
// committed tables back (f64packSame). Needs the venv interpreter with
// mpmath (tier 2: a failure without it).
func TestLibmGeneratorReproduces(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.AstropyPython(t)
	if b, err := exec.Command(py, "-c", "import mpmath").CombinedOutput(); err != nil {
		t.Fatalf("%s cannot import mpmath (%v: %s)", py, err, b)
	}
	out := t.TempDir()
	if b, err := exec.Command(py, filepath.Join(wrefDir, "gen_libm.py"), out).CombinedOutput(); err != nil {
		t.Fatalf("gen_libm.py: %v\n%s", err, b)
	}
	f64packSame(t, filepath.Join(out, "libm", "trig.json"), libmGolden)
	f64packSame(t, filepath.Join(out, "libm", "atan2.json"), libmAtan2Golden)
}

// TestPrjCaseGeneratorReproduce re-runs gen_prj_cases.py (stdlib Python)
// into a temp directory and requires exactly the committed prj_*/pvfid_*/prjswap_*
// case files back. Tier 2: a failure without python3.
func TestPrjCaseGeneratorReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py, err := exec.LookPath("python3")
	if err != nil {
		t.Fatalf("python3 not on PATH: %v", err)
	}
	out := t.TempDir()
	if b, err := exec.Command(py, filepath.Join(wrefDir, "gen_prj_cases.py"), out).CombinedOutput(); err != nil {
		t.Fatalf("gen_prj_cases.py: %v\n%s", err, b)
	}
	diffGenerated(t, filepath.Join(out, "hdrcases"), hdrCaseDir, "prj_*.case")
	diffGenerated(t, filepath.Join(out, "hdrcases"), hdrCaseDir, "pvfid_*.case")
	diffGenerated(t, filepath.Join(out, "hdrcases"), hdrCaseDir, "prjswap_*.case")
}

// TestHighOrderCaseGeneratorReproduce re-runs gen_highorder_cases.py
// (stdlib Python) into a temp directory and requires exactly the committed
// hiorder_*.case files back (SIP order 7 and 9, TPV with every PVi_0..39).
// Tier 2: a failure without python3.
func TestHighOrderCaseGeneratorReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py, err := exec.LookPath("python3")
	if err != nil {
		t.Fatalf("python3 not on PATH: %v", err)
	}
	out := t.TempDir()
	if b, err := exec.Command(py, filepath.Join(wrefDir, "gen_highorder_cases.py"), out).CombinedOutput(); err != nil {
		t.Fatalf("gen_highorder_cases.py: %v\n%s", err, b)
	}
	diffGenerated(t, filepath.Join(out, "hdrcases"), hdrCaseDir, "hiorder_*.case")
	if n := len(globNames(t, hdrCaseDir, "hiorder_*.case", ".case")); n != 4 {
		t.Fatalf("%d hiorder_*.case files, want 4", n)
	}
}

// TestSampleGeneratorsReproduce re-runs gen_sample_cases.py (needs the
// astropy 8.0.2 test-data headers, healpy and /usr/bin/HPXcvt) into a temp
// directory and requires exactly the committed sample_*.case files and the
// committed HPXcvt image + healpy truth back, byte for byte. Tier 2: a
// failure without those tools.
func TestSampleGeneratorsReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.AstropyPython(t)
	if b, err := exec.Command(py, "-c", "import healpy").CombinedOutput(); err != nil {
		t.Fatalf("%s cannot import healpy (%v: %s)", py, err, b)
	}
	testutil.RequireTool(t, "HPXcvt")
	out := t.TempDir()
	if b, err := exec.Command(py, filepath.Join(wrefDir, "gen_sample_cases.py"), out).CombinedOutput(); err != nil {
		t.Fatalf("gen_sample_cases.py: %v\n%s", err, b)
	}
	diffGenerated(t, filepath.Join(out, "hdrcases"), hdrCaseDir, "sample_*.case")
	diffGenerated(t, filepath.Join(out, "samples"), wrefSampleDir, "*")
}

// TestWrefGoldensReproduce rebuilds wref against the installed wcslib and
// requires it to regenerate every committed bare-projection golden (the
// metadata .json and the payload .bin) byte for byte, with no golden
// missing and none orphaned. Needs a C compiler and libwcs (tier 2: a
// failure without them).
func TestWrefGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	unavailable := t.Fatalf
	if _, err := exec.LookPath("make"); err != nil {
		unavailable("make not found")
	}
	out, err := exec.Command("make", "-C", "testdata/wref", "wref").CombinedOutput()
	if err != nil {
		unavailable("cannot build wref (needs cc + wcslib %s): %v\n%s", wrefWcslibVersion, err, out)
	}
	bin, _ := filepath.Abs(filepath.Join("testdata", "wref", "wref"))
	dir := t.TempDir()
	if out, err := exec.Command(bin, "gen", dir).CombinedOutput(); err != nil {
		t.Fatalf("wref gen: %v\n%s", err, out)
	}
	list, err := exec.Command(bin, "list").Output()
	if err != nil {
		t.Fatalf("wref list: %v", err)
	}
	want := strings.Fields(string(list))
	sort.Strings(want)
	have := wrefNames(t)
	if strings.Join(want, ",") != strings.Join(have, ",") {
		t.Errorf("committed goldens %v != wref cases %v", have, want)
	}
	// The metadata and the uncompressed payload must come back byte for
	// byte. The compressed .bin must too when the local zlib is the one
	// that wrote the committed golden (zlib_version in the metadata); with
	// a different zlib the deflate stream may legitimately differ, and
	// only the payload it decodes to is the golden's content.
	for _, name := range want {
		gm, graw := wrefPayload(t, dir, name)
		cm, craw := wrefPayload(t, wrefGoldenDir, name)
		gj, _ := os.ReadFile(filepath.Join(dir, name+".json"))
		cj, _ := os.ReadFile(filepath.Join(wrefGoldenDir, name+".json"))
		sameZlib := gm.ZlibVersion == cm.ZlibVersion
		if !sameZlib {
			t.Logf("%s: local zlib %s, golden written with %s: comparing the decoded payload only", name, gm.ZlibVersion, cm.ZlibVersion)
			gj = bytes.Replace(gj, []byte(`"zlib_version": "`+gm.ZlibVersion+`"`), []byte(`"zlib_version": "`+cm.ZlibVersion+`"`), 1)
		}
		if !bytes.Equal(gj, cj) {
			t.Errorf("%s.json: wref output differs from the committed golden (regenerate with make wref-goldens and review the diff)", name)
		}
		if !bytes.Equal(graw, craw) {
			t.Errorf("%s.bin: wref payload differs from the committed golden's (regenerate with make wref-goldens and review the diff)", name)
		}
		if sameZlib {
			gz, _ := os.ReadFile(filepath.Join(dir, name+".bin"))
			cz, _ := os.ReadFile(filepath.Join(wrefGoldenDir, name+".bin"))
			if !bytes.Equal(gz, cz) {
				t.Errorf("%s.bin: same zlib %s but the compressed bytes differ from the committed golden", name, gm.ZlibVersion)
			}
		}
	}
	if out, err := exec.Command(bin, "trig", dir).CombinedOutput(); err != nil {
		t.Fatalf("wref trig: %v\n%s", err, out)
	}
	if got, err := os.ReadFile(filepath.Join(dir, "wcstrig_guards.json")); err != nil {
		t.Error(err)
	} else if committed, err := os.ReadFile(wcstrigGuardGolden); err != nil || !bytes.Equal(got, committed) {
		t.Errorf("wref trig output differs from the committed %s (err %v)", wcstrigGuardGolden, err)
	}
	for _, name := range []string{wrefPrjsetStatus} {
		got, err := os.ReadFile(filepath.Join(dir, name+".json"))
		if err != nil {
			t.Errorf("%s: %v", name, err)
			continue
		}
		committed, err := os.ReadFile(filepath.Join(wrefGoldenDir, name+".json"))
		if err != nil {
			t.Errorf("%s: committed golden missing: %v", name, err)
			continue
		}
		if !bytes.Equal(got, committed) {
			t.Errorf("%s: wref output differs from the committed golden (regenerate with make wref-goldens and review the diff)", name)
		}
	}
}

// TestXPHPaperGeneratorReproduces re-runs gen_xph_paper.py and requires
// the committed golden back byte for byte (tier 2: needs the venv's
// mpmath, a failure without it).
func TestXPHPaperGeneratorReproduces(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.AstropyPython(t)
	if b, err := exec.Command(py, "-c", "import mpmath").CombinedOutput(); err != nil {
		t.Fatalf("%s cannot import mpmath (%v: %s)", py, err, b)
	}
	out := t.TempDir()
	if b, err := exec.Command(py, filepath.Join(wrefDir, "gen_xph_paper.py"), out).CombinedOutput(); err != nil {
		t.Fatalf("gen_xph_paper.py: %v\n%s", err, b)
	}
	f64packSame(t, filepath.Join(out, "paper", "xph_paper.json"), xphPaperGolden)
}
