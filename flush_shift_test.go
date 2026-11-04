package fits

import (
	"bytes"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// The tests in this file prove that a Flush that has to resize a header
// (move every later byte of the file) is correct and crash safe:
//
//   - Correct: the file a fault-free Flush produces is read by cfitsio (the
//     cfref harness) as exactly the input file (its committed cfitsio
//     golden) plus the requested card changes. Every other card, every
//     pixel and every table cell is unchanged, and fitsverify reports
//     nothing new.
//   - Crash safe: a crash is simulated after each step of the Flush (the
//     flushFaultHook stops Flush right there and the descriptor is closed
//     with no cleanup, exactly what a process death leaves on disk). After
//     OpenForEdit runs its recovery, the file must be byte-identical to
//     either the original file or the fault-free result, and once a step
//     has produced the result every later step must too (a completed edit
//     is never thrown away).

// shiftScenario is one Flush that needs a header resize.
type shiftScenario struct {
	name   string
	sample string // committed testdata sample the edit starts from
	// mutate applies the edit through the public Header() surface.
	mutate func(t *testing.T, f *File)
	// expect turns the cfitsio golden card list of HDU i into the card
	// list cfitsio must read after the edit (nil = HDU unchanged).
	expect map[int]func(cards []testutil.CfrefCard) []wantCard
	// blocks is the header block count each edited HDU must end with.
	blocks map[int]int64
}

// wantCard is one expected card: either an unchanged card (raw text as
// cfitsio reported it in the golden) or a new/changed card identified by
// cfitsio's parsed key, value and comment.
type wantCard struct {
	raw     string // set for unchanged cards
	key     string
	value   string // fits_parse_value text, quotes kept for strings
	comment string
}

func unchanged(cards []testutil.CfrefCard) []wantCard {
	out := make([]wantCard, len(cards))
	for i, c := range cards {
		out[i] = wantCard{raw: string(c.Raw)}
	}
	return out
}

func hduHeader(t *testing.T, f *File, i int) *header.Header {
	t.Helper()
	h, err := f.HDU(i)
	if err != nil {
		t.Fatalf("HDU(%d): %v", i, err)
	}
	return h.Header()
}

// growHDU1 adds 30 cards to testprog.std HDU 1 (a binary table whose 44
// cards + END fill 2 blocks): 75 cards need 3 blocks, so the 6 HDUs after
// it move forward by 2880 bytes.
func growHDU1(t *testing.T, f *File) {
	h := hduHeader(t, f, 1)
	for i := 1; i <= 30; i++ {
		if err := h.Set(fmt.Sprintf("GROW%04d", i), int64(i), fmt.Sprintf("grow card %d", i)); err != nil {
			t.Fatal(err)
		}
	}
}

func expectGrowHDU1(cards []testutil.CfrefCard) []wantCard {
	out := unchanged(cards)
	for i := 1; i <= 30; i++ {
		out = append(out, wantCard{key: fmt.Sprintf("GROW%04d", i), value: fmt.Sprint(i), comment: fmt.Sprintf("grow card %d", i)})
	}
	return out
}

// shrinkHDU4 edits testprog.std HDU 4, whose 38 keywords are followed by
// 34 blank cards (reserve space left by cfitsio) and END at card 73, so
// the header fills 3 blocks. Changing one value and deleting one blank
// card leaves 71 cards + END = 72 cards = 2 blocks, so the 3 HDUs after it
// move back by 2880 bytes. cfitsio does not count trailing blank cards as
// keywords, so it must still read the same 38 keywords.
func shrinkHDU4(t *testing.T, f *File) {
	h := hduHeader(t, f, 4)
	if err := h.Set("KEY_PREC", "changed by shrink test", "comment here"); err != nil {
		t.Fatal(err)
	}
	if err := h.Delete(""); err != nil {
		t.Fatal(err)
	}
}

func expectShrinkHDU4(cards []testutil.CfrefCard) []wantCard {
	out := unchanged(cards)
	for i, c := range cards {
		if string(c.Key) == "KEY_PREC" {
			out[i] = wantCard{key: "KEY_PREC", value: "'changed by shrink test'", comment: "comment here"}
		}
	}
	return out
}

// inPlaceHDU2 changes a value in HDU 2 (an image extension) without
// changing its block count.
func inPlaceHDU2(t *testing.T, f *File) {
	h := hduHeader(t, f, 2)
	if err := h.Set("INPLACE", "same block count", "in-place edit"); err != nil {
		t.Fatal(err)
	}
}

func expectInPlaceHDU2(cards []testutil.CfrefCard) []wantCard {
	return append(unchanged(cards), wantCard{key: "INPLACE", value: "'same block count'", comment: "in-place edit"})
}

var shiftScenarios = []shiftScenario{
	{
		// No block count changes: the undo-journal path.
		name: "inplace", sample: "testprog.std", mutate: inPlaceHDU2,
		expect: map[int]func([]testutil.CfrefCard) []wantCard{2: expectInPlaceHDU2},
		blocks: map[int]int64{2: 1},
	},
	{
		// Two headers rewritten in place by one Flush (HDU 2 image, HDU 3
		// ASCII table): both or neither.
		name: "inplace-two-hdus", sample: "testprog.std",
		mutate: func(t *testing.T, f *File) {
			inPlaceHDU2(t, f)
			if err := hduHeader(t, f, 3).Set("EXTNAME", "renamed-ascii", "new name"); err != nil {
				t.Fatal(err)
			}
		},
		expect: map[int]func([]testutil.CfrefCard) []wantCard{
			2: expectInPlaceHDU2,
			3: func(cards []testutil.CfrefCard) []wantCard {
				out := unchanged(cards)
				for i, c := range cards {
					if string(c.Key) == "EXTNAME" {
						out[i] = wantCard{key: "EXTNAME", value: "'renamed-ascii'", comment: "new name"}
					}
				}
				return out
			},
		},
		blocks: map[int]int64{2: 1, 3: 1},
	},
	{
		name: "grow", sample: "testprog.std", mutate: growHDU1,
		expect: map[int]func([]testutil.CfrefCard) []wantCard{1: expectGrowHDU1},
		blocks: map[int]int64{1: 3},
	},
	{
		name: "shrink", sample: "testprog.std", mutate: shrinkHDU4,
		expect: map[int]func([]testutil.CfrefCard) []wantCard{4: expectShrinkHDU4},
		blocks: map[int]int64{4: 2},
	},
	{
		// One Flush that grows HDU 1, rewrites HDU 2 in place and shrinks
		// HDU 4: recovery must not leave some of these applied and others
		// not.
		name: "grow+inplace+shrink", sample: "testprog.std",
		mutate: func(t *testing.T, f *File) { growHDU1(t, f); inPlaceHDU2(t, f); shrinkHDU4(t, f) },
		expect: map[int]func([]testutil.CfrefCard) []wantCard{1: expectGrowHDU1, 2: expectInPlaceHDU2, 4: expectShrinkHDU4},
		blocks: map[int]int64{1: 3, 2: 1, 4: 2},
	},
	{
		// The primary HDU of iter_a.fit (1 block, 20 cards) grows to 2
		// blocks; the binary table after it moves.
		name: "grow-primary", sample: "iter_a.fit",
		mutate: func(t *testing.T, f *File) {
			h := hduHeader(t, f, 0)
			for i := 1; i <= 20; i++ {
				if err := h.Set(fmt.Sprintf("PGROW%03d", i), float64(i)+0.5, "primary grow"); err != nil {
					t.Fatal(err)
				}
			}
		},
		expect: map[int]func([]testutil.CfrefCard) []wantCard{0: func(cards []testutil.CfrefCard) []wantCard {
			out := unchanged(cards)
			for i := 1; i <= 20; i++ {
				out = append(out, wantCard{key: fmt.Sprintf("PGROW%03d", i), value: fmt.Sprintf("%g", float64(i)+0.5), comment: "primary grow"})
			}
			return out
		}},
		blocks: map[int]int64{0: 2},
	},
}

func copyFileT(t *testing.T, src, dst string) {
	t.Helper()
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(dst, b, 0o644); err != nil {
		t.Fatal(err)
	}
}

// editAndFlush opens path for edit, applies mutate and flushes. It returns
// Flush's error. The file is closed without any cleanup either way.
func editAndFlush(t *testing.T, path string, mutate func(*testing.T, *File)) error {
	t.Helper()
	f, err := OpenForEdit(path)
	if err != nil {
		t.Fatalf("OpenForEdit: %v", err)
	}
	mutate(t, f)
	ferr := f.Flush()
	if cerr := f.Close(); cerr != nil && ferr == nil {
		ferr = cerr
	}
	return ferr
}

// sidecars lists every directory entry other than the FITS file itself.
func sidecars(t *testing.T, dir, keep string) []string {
	t.Helper()
	ents, err := os.ReadDir(dir)
	if err != nil {
		t.Fatal(err)
	}
	var out []string
	for _, e := range ents {
		if e.Name() != keep {
			out = append(out, e.Name())
		}
	}
	return out
}

var errInjectedCrash = errors.New("injected crash")

// flushSteps runs scenario sc fault-free on a copy of the sample, records
// the fault-injection steps Flush passed through, and returns them with
// the resulting bytes.
func flushSteps(t *testing.T, sc shiftScenario) (steps []string, post []byte) {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, sc.sample)
	copyFileT(t, testutil.TestdataPath(t, sc.sample), path)
	flushFaultHook = func(step string) error { steps = append(steps, step); return nil }
	defer func() { flushFaultHook = nil }()
	if err := editAndFlush(t, path, sc.mutate); err != nil {
		t.Fatalf("fault-free Flush: %v", err)
	}
	flushFaultHook = nil
	if s := sidecars(t, dir, sc.sample); len(s) != 0 {
		t.Fatalf("fault-free Flush left files behind: %v", s)
	}
	post, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	return steps, post
}

// TestFlushHeaderResizeMatchesCfitsio proves the fault-free result: cfitsio
// reads the edited file as the original (committed golden) plus exactly
// the requested card changes, with every pixel and table cell unchanged,
// and fitsverify finds nothing it did not find in the original.
func TestFlushHeaderResizeMatchesCfitsio(t *testing.T) {
	for _, sc := range shiftScenarios {
		t.Run(sc.name, func(t *testing.T) {
			src := testutil.TestdataPath(t, sc.sample)
			testutil.RequireFile(t, src)
			want := testutil.LoadCfrefGolden(t, sc.sample)
			_, post := flushSteps(t, sc)
			dir := t.TempDir()
			outPath := filepath.Join(dir, "post.fits")
			if err := os.WriteFile(outPath, post, 0o644); err != nil {
				t.Fatal(err)
			}

			got := testutil.CfrefDump(t, outPath)
			compareEditedDump(t, want, got, sc)

			// Our own reader must agree with cfitsio on the new layout.
			f, err := Open(outPath)
			if err != nil {
				t.Fatalf("Open edited file: %v", err)
			}
			defer f.Close()
			if f.NumHDU() != got.NHDU {
				t.Fatalf("NumHDU %d, cfitsio %d", f.NumHDU(), got.NHDU)
			}
			for i, h := range got.HDUs {
				r := f.hdus[i]
				if r.headerStart != h.HeaderStart || r.dataStart != h.DataStart || r.paddedEnd != h.DataEnd {
					t.Errorf("HDU %d layout %d/%d/%d, cfitsio %d/%d/%d", i, r.headerStart, r.dataStart, r.paddedEnd, h.HeaderStart, h.DataStart, h.DataEnd)
				}
			}

			before := testutil.Fitsverify(t, src)
			after := testutil.Fitsverify(t, outPath)
			if after.Errors != before.Errors || after.Warnings != before.Warnings {
				t.Errorf("fitsverify: original %d warnings / %d errors, edited %d / %d\n%s",
					before.Warnings, before.Errors, after.Warnings, after.Errors, after.Report)
			}
		})
	}
}

// compareEditedDump checks a cfitsio dump of an edited file against the
// golden of the original plus the scenario's expected card changes.
func compareEditedDump(t *testing.T, want, got *testutil.CfrefFile, sc shiftScenario) {
	t.Helper()
	if got.NHDU != want.NHDU {
		t.Fatalf("cfitsio sees %d HDUs, original has %d", got.NHDU, want.NHDU)
	}
	var shift int64
	for i := range want.HDUs {
		w, g := want.HDUs[i], got.HDUs[i]
		if g.HeaderStart != w.HeaderStart+shift {
			t.Errorf("HDU %d starts at %d, want %d (original %d shifted by %d)", i, g.HeaderStart, w.HeaderStart+shift, w.HeaderStart, shift)
		}
		hdrBlocks := (g.DataStart - g.HeaderStart) / 2880
		if nb, ok := sc.blocks[i]; ok && hdrBlocks != nb {
			t.Errorf("HDU %d header has %d blocks, want %d", i, hdrBlocks, nb)
		}
		shift += (g.DataStart - g.HeaderStart) - (w.DataStart - w.HeaderStart)
		if g.DataEnd-g.DataStart != w.DataEnd-w.DataStart {
			t.Errorf("HDU %d data length %d, original %d", i, g.DataEnd-g.DataStart, w.DataEnd-w.DataStart)
		}
		if !reflect.DeepEqual(g.Image, w.Image) {
			t.Errorf("HDU %d: cfitsio image values differ from the original", i)
		}
		if !reflect.DeepEqual(g.Table, w.Table) {
			t.Errorf("HDU %d: cfitsio table values differ from the original", i)
		}
		if !reflect.DeepEqual(g.Errors, w.Errors) {
			t.Errorf("HDU %d: cfitsio errors %v, original %v", i, g.Errors, w.Errors)
		}
		exp := unchanged(w.Cards)
		if f := sc.expect[i]; f != nil {
			exp = f(w.Cards)
		}
		if len(g.Cards) != len(exp) {
			t.Errorf("HDU %d: cfitsio reads %d cards, want %d", i, len(g.Cards), len(exp))
			continue
		}
		for j, e := range exp {
			c := g.Cards[j]
			if e.raw != "" {
				if string(c.Raw) != e.raw {
					t.Errorf("HDU %d card %d = %q, want unchanged %q", i, j+1, c.Raw, e.raw)
				}
				continue
			}
			if string(c.Key) != e.key || string(c.Value) != e.value || strings.TrimSpace(string(c.Comment)) != e.comment || c.Error != "" || c.ParseError != "" {
				t.Errorf("HDU %d card %d = %q (key %q value %q comment %q err %q/%q), want key %q value %q comment %q",
					i, j+1, c.Raw, c.Key, c.Value, c.Comment, c.Error, c.ParseError, e.key, e.value, e.comment)
			}
		}
	}
}

// crashRecoverers are the two documented ways a crashed edit is recovered:
// OpenForEdit, and EditFile (which recovers before it opens the file; here
// with a callback that copies every HDU unchanged, so a correctly recovered
// file comes out byte-identical).
var crashRecoverers = []struct {
	name string
	run  func(path string) error
}{
	{"OpenForEdit", func(path string) error {
		f, err := OpenForEdit(path)
		if err != nil {
			return err
		}
		return f.Close()
	}},
	{"EditFile", func(path string) error { return EditFile(path, copyAll) }},
}

// TestFlushHeaderResizeCrashAtEveryStep injects a crash after every step
// of a resizing Flush and checks what OpenForEdit, and separately EditFile,
// recovers.
func TestFlushHeaderResizeCrashAtEveryStep(t *testing.T) {
	for _, sc := range shiftScenarios {
		t.Run(sc.name, func(t *testing.T) {
			src := testutil.TestdataPath(t, sc.sample)
			testutil.RequireFile(t, src)
			pre, err := os.ReadFile(src)
			if err != nil {
				t.Fatal(err)
			}
			steps, post := flushSteps(t, sc)
			if bytes.Equal(pre, post) {
				t.Fatal("scenario does not change the file")
			}
			if len(steps) == 0 {
				t.Fatal("Flush passed no fault-injection steps")
			}
			t.Logf("steps: %v", steps)
			for _, rv := range crashRecoverers {
				t.Run(rv.name, func(t *testing.T) {
					reachedPost := false
					for k, step := range steps {
						dir := t.TempDir()
						path := filepath.Join(dir, sc.sample)
						copyFileT(t, src, path)
						n := 0
						flushFaultHook = func(s string) error {
							if n == k {
								return fmt.Errorf("%w at %s", errInjectedCrash, s)
							}
							n++
							return nil
						}
						ferr := editAndFlush(t, path, sc.mutate)
						flushFaultHook = nil
						if !errors.Is(ferr, errInjectedCrash) {
							t.Fatalf("step %d (%s): Flush returned %v, want the injected crash", k, step, ferr)
						}

						// Before recovery a read-only Open must not show a
						// half-applied Flush: it either refuses (a journal is
						// pending) or the bytes on disk are the original or
						// the edited file.
						checkReadOnlyOpenAfterCrash(t, step, path, pre, post)

						if err := rv.run(path); err != nil {
							t.Errorf("crash after %q: %s after crash: %v", step, rv.name, err)
							continue
						}
						got, err := os.ReadFile(path)
						if err != nil {
							t.Fatal(err)
						}
						isPre, isPost := bytes.Equal(got, pre), bytes.Equal(got, post)
						switch {
						case isPost:
							reachedPost = true
						case isPre:
							if reachedPost {
								t.Errorf("crash after %q: recovery returned the original file although an earlier step already had the edit committed (completed edit destroyed)", step)
							}
						default:
							t.Errorf("crash after %q: recovered file is neither the original nor the edited file (len %d, original %d, edited %d, first diff vs original at %d)",
								step, len(got), len(pre), len(post), firstDiff(got, pre))
						}
						if s := sidecars(t, dir, sc.sample); len(s) != 0 {
							t.Errorf("crash after %q: files left after recovery: %v", step, s)
						}
					}
					if !reachedPost {
						t.Errorf("no crash point yielded the edited file: a crash after the last step (%s) lost a completed edit", steps[len(steps)-1])
					}
				})
			}
		})
	}
}

// firstDiff is defined in image_write_readback_test.go.

// A resizing Flush replaces the file through a temp file and a rename. The
// new file keeps the old permission bits exactly (not widened, not
// narrowed by the umask), and a symbolic link to the file stays a link
// whose target holds the edited file.
func TestFlushHeaderResizeKeepsModeAndSymlink(t *testing.T) {
	for _, mode := range []os.FileMode{0o600, 0o640, 0o666, 0o606} {
		t.Run(mode.String(), func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			if err := os.Chmod(path, mode); err != nil {
				t.Fatal(err)
			}
			if err := editAndFlush(t, path, growHDU1); err != nil {
				t.Fatal(err)
			}
			st, err := os.Stat(path)
			if err != nil {
				t.Fatal(err)
			}
			if st.Mode().Perm() != mode {
				t.Fatalf("mode after resizing Flush = %v, want %v", st.Mode().Perm(), mode)
			}
			if len(mustRead(t, path)) != len(orig)+2880 {
				t.Fatal("Flush did not resize the header")
			}
		})
	}
	t.Run("symlink", func(t *testing.T) {
		path, orig := copySample(t, "testprog.std")
		link := filepath.Join(t.TempDir(), "link.fits")
		if err := os.Symlink(path, link); err != nil {
			t.Fatal(err)
		}
		if err := editAndFlush(t, link, growHDU1); err != nil {
			t.Fatal(err)
		}
		if st, err := os.Lstat(link); err != nil || st.Mode()&os.ModeSymlink == 0 {
			t.Fatalf("link replaced by a regular file (%v)", err)
		}
		plain, _ := copySample(t, "testprog.std")
		if err := editAndFlush(t, plain, growHDU1); err != nil {
			t.Fatal(err)
		}
		want := mustRead(t, plain)
		if got := mustRead(t, path); len(want) != len(orig)+2880 || !bytes.Equal(got, want) {
			t.Fatalf("link target (%d bytes) is not the file a Flush through the real path gives (%d bytes)", len(got), len(want))
		}
		if s := sidecars(t, filepath.Dir(path), filepath.Base(path)); len(s) != 0 {
			t.Fatalf("files left next to the target: %v", s)
		}
	})
}

// checkReadOnlyOpenAfterCrash opens path read-only after an injected crash.
// While a journal is pending Open must fail with ErrPendingJournal (it
// cannot roll back a file it only reads); otherwise the file must already
// be exactly pre or post.
func checkReadOnlyOpenAfterCrash(t *testing.T, step, path string, pre, post []byte) {
	t.Helper()
	raw := mustRead(t, path)
	consistent := bytes.Equal(raw, pre) || bytes.Equal(raw, post)
	_, jerr := os.Stat(journalPath(path))
	pending := jerr == nil
	rf, oerr := Open(path)
	if oerr == nil {
		rf.Close()
	}
	switch {
	case oerr == nil && !consistent:
		t.Errorf("crash after %q: read-only Open succeeded on a half-applied Flush (first diff vs original at %d)", step, firstDiff(raw, pre))
	case oerr == nil && pending:
		t.Errorf("crash after %q: read-only Open ignored the pending journal", step)
	case oerr != nil && !errors.Is(oerr, ErrPendingJournal):
		t.Errorf("crash after %q: read-only Open: %v, want ErrPendingJournal or success", step, oerr)
	case oerr != nil && !pending:
		t.Errorf("crash after %q: read-only Open refused with no journal pending: %v", step, oerr)
	}
}

// TestFlushKeepsBytesAfterLastHDU: the bytes after the last HDU (a
// trailing zero block, a junk block, a few junk bytes: framing/ samples,
// which cfitsio opens ignoring them) and a last data unit that lacks its
// final padding (missing_final_pad*) survive a Flush exactly as they were,
// whether the header is overwritten in place or has to grow (temp-file
// rewrite): everything after the edited header is byte-identical to what
// followed it before, so the two kinds of Flush differ only in the header
// itself. TrailingBytes reports the same count afterwards.
func TestFlushKeepsBytesAfterLastHDU(t *testing.T) {
	samples := []string{"trail_zero_block.fits", "trail_junk_block.fits", "trail_junk_bytes.fits",
		"missing_final_pad.fits", "missing_final_pad_ext.fits"}
	edits := map[string]func(t *testing.T, f *File){
		"in-place": func(t *testing.T, f *File) {
			if err := hduHeader(t, f, 0).Set("OBJECT", "edited", ""); err != nil {
				t.Fatal(err)
			}
		},
		"resize": func(t *testing.T, f *File) {
			h := hduHeader(t, f, 0)
			for i := 1; i <= 40; i++ {
				if err := h.Set(fmt.Sprintf("GROW%04d", i), int64(i), ""); err != nil {
					t.Fatal(err)
				}
			}
		},
	}
	for _, s := range samples {
		for name, edit := range edits {
			t.Run(s+"/"+name, func(t *testing.T) {
				src := imagesPath(t, "framing", s)
				orig := mustRead(t, src)
				rf, err := Open(src)
				if err != nil {
					t.Fatal(err)
				}
				trailing := rf.TrailingBytes()
				oldEnd := int64(len(rf.hdus[0].rawHeader))
				rf.Close()
				path := filepath.Join(t.TempDir(), s)
				copyFileT(t, src, path)
				if err := editAndFlush(t, path, edit); err != nil {
					t.Fatal(err)
				}
				got := mustRead(t, path)
				ef, err := Open(path)
				if err != nil {
					t.Fatal(err)
				}
				newEnd := int64(len(ef.hdus[0].rawHeader))
				gotTrailing := ef.TrailingBytes()
				ef.Close()
				if (name == "resize") != (newEnd != oldEnd) {
					t.Fatalf("header %d -> %d bytes: not a %s Flush", oldEnd, newEnd, name)
				}
				if !bytes.Equal(got[newEnd:], orig[oldEnd:]) {
					t.Errorf("bytes after the header differ (%d after, %d before; first diff at %d)",
						len(got)-int(newEnd), len(orig)-int(oldEnd), firstDiff(got[newEnd:], orig[oldEnd:]))
				}
				if gotTrailing != trailing {
					t.Errorf("TrailingBytes %d after Flush, %d before", gotTrailing, trailing)
				}
			})
		}
	}
}
