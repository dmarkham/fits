package fits

import (
	"bytes"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Crash and error behaviour of the writes that change a file opened with
// OpenForEdit in place: OverwritePixels and the appends (AppendImage /
// WriteImage, AppendBinaryTable, AppendCompressedImage), next to the
// in-place header Flush they share the undo journal with (journal.go).
//
// Each operation is modelled two ways:
//
//   - a process crash after each named step (flushFaultHook), as in
//     TestFlushHeaderResizeCrashAtEveryStep;
//   - a crash part way through the data write itself: the file holds the
//     original bytes with the first L bytes of the write replaced by the
//     new ones, for every L in steps of one card (80 bytes), with the
//     journal the call wrote before writing (if any).
//
// In every state a read-only Open must refuse (ErrPendingJournal) or see
// exactly the original or the finished file, and after OpenForEdit /
// EditFile recovery the file must be byte-identical to one of the two,
// and to the finished file only once the call had committed. The truth is
// byte equality with the original sample and the fault-free result.

// inPlaceOp is one in-place write on a copy of a committed sample.
type inPlaceOp struct {
	name   string
	sample string
	// junk is appended to the sample first (trailing bytes an append
	// overwrites).
	junk []byte
	// cut, when > 0, truncates the sample to that many bytes first (a
	// data unit that lacks its final padding, which the write restores).
	cut int
	run func(t *testing.T, f *File) error
	// recoverers to check (EditFile rewrites a file and drops trailing
	// bytes, so it is not a byte-identical recovery for junk).
	editFileRecovers bool
}

func iterImageOverwrite(t *testing.T, f *File) error {
	p, err := f.Primary()
	if err != nil {
		t.Fatal(err)
	}
	px := make([]int32, p.pixelCount())
	for i := range px {
		px[i] = int32(i*7919) - 1<<30
	}
	return OverwritePixels(p, px)
}

func appendInt16Image(t *testing.T, f *File) error {
	px := make([]int16, 100*50)
	for i := range px {
		px[i] = int16(i*31 - 20000)
	}
	_, err := AppendImage(f, nil, []int64{100, 50}, px)
	return err
}

func appendTable3(t *testing.T, f *File) error {
	_, err := AppendBinaryTable(f, nil, []ColumnData{
		{Name: "J", DataInt32: []int32{1, -2, 3}},
		{Name: "V", DataVarInt32: [][]int32{{1, 2}, {3}, {}}},
		{Name: "S", DataString: []string{"abc", "def", "ghi"}},
	})
	return err
}

func appendCompressed(t *testing.T, f *File) error {
	px := make([]int16, 64*40)
	for i := range px {
		px[i] = int16(i%251 - 100)
	}
	_, err := AppendCompressedImage(f, nil, []int64{64, 40}, px, CompressOptions{})
	return err
}

// f32PlainOverwrite overwrites the 8x8 float32 primary of
// images/cfwrite/f32_plain.fits.
func f32PlainOverwrite(t *testing.T, f *File) error {
	p, err := f.Primary()
	if err != nil {
		t.Fatal(err)
	}
	px := make([]float32, p.pixelCount())
	for i := range px {
		px[i] = float32(i)*0.25 - 3
	}
	return OverwritePixels(p, px)
}

func appendASCII(t *testing.T, f *File) error {
	_, err := AppendASCIITable(f, nil, []ASCIIColumnData{
		{Name: "N", TForm: "I6", DataInt64: []int64{1, -22, 333}},
		{Name: "X", TForm: "E14.6", DataFloat64: []float64{1.5, -2.25, 1e10}},
		{Name: "S", TForm: "A5", DataString: []string{"ab", "cdefg", ""}},
	})
	return err
}

func appendCompressedFloat32Op(t *testing.T, f *File) error {
	px := make([]float32, 48*30)
	for i := range px {
		px[i] = float32(i%97)*0.5 - 11
	}
	_, err := AppendCompressedFloat32Image(f, nil, []int64{48, 30}, px, CompressFloatOptions{})
	return err
}

var inPlaceOps = []inPlaceOp{
	{name: "OverwritePixels", sample: "iter_image.fit", run: iterImageOverwrite, editFileRecovers: true},
	// The data unit lacks its final padding (2924 of 5760 bytes): the
	// write grows the file, so rollback must also cut it back (the
	// journal's Truncate flag).
	{name: "OverwritePixelsUnpadded", sample: "images/cfwrite/f32_plain.fits", cut: 2924, run: f32PlainOverwrite},
	{name: "AppendImage", sample: "testprog.std", run: appendInt16Image, editFileRecovers: true},
	{name: "AppendBinaryTable", sample: "testprog.std", run: appendTable3, editFileRecovers: true},
	{name: "AppendASCIITable", sample: "testprog.std", run: appendASCII, editFileRecovers: true},
	{name: "AppendCompressedImage", sample: "testprog.std", run: appendCompressed, editFileRecovers: true},
	{name: "AppendCompressedFloat32Image", sample: "testprog.std", run: appendCompressedFloat32Op, editFileRecovers: true},
	// An append overwrites trailing bytes that are not part of any HDU;
	// recovery must put them back too.
	{name: "AppendImageOverJunk", sample: "testprog.std", junk: bytes.Repeat([]byte("junk"), 1500), run: appendInt16Image},
}

// setupOp writes the op's starting file into a new temp dir.
func setupOp(t *testing.T, op inPlaceOp) (dir, path string, pre []byte) {
	t.Helper()
	src := testutil.TestdataPath(t, op.sample)
	testutil.RequireFile(t, src)
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatal(err)
	}
	if op.cut > 0 {
		if op.cut >= len(b) {
			t.Fatalf("%s: cut %d >= sample size %d", op.sample, op.cut, len(b))
		}
		b = b[:op.cut]
	}
	pre = append(b, op.junk...)
	dir = t.TempDir()
	path = filepath.Join(dir, op.file())
	if err := os.WriteFile(path, pre, 0o644); err != nil {
		t.Fatal(err)
	}
	return dir, path, pre
}

// file is the name of the op's copy of its sample.
func (op inPlaceOp) file() string { return filepath.Base(op.sample) }

// runOp opens path for edit, runs the op and closes without cleanup.
func runOp(t *testing.T, op inPlaceOp, path string) error {
	t.Helper()
	f, err := OpenForEdit(path)
	if err != nil {
		t.Fatalf("OpenForEdit: %v", err)
	}
	rerr := op.run(t, f)
	if cerr := f.Close(); cerr != nil && rerr == nil {
		rerr = cerr
	}
	return rerr
}

func (op inPlaceOp) recoverers() []struct {
	name string
	run  func(path string) error
} {
	if op.editFileRecovers {
		return crashRecoverers
	}
	return crashRecoverers[:1]
}

// writeStart is where the op's data write begins in pre/post: the first
// differing byte.
func writeRange(pre, post []byte) (start, end int) {
	start = firstDiff(pre, post)
	end = max(len(pre), len(post))
	for end > start && end <= len(pre) && end <= len(post) && pre[end-1] == post[end-1] {
		end--
	}
	return start, end
}

func TestInPlaceWriteCrashAtEveryStep(t *testing.T) {
	for _, op := range inPlaceOps {
		t.Run(op.name, func(t *testing.T) {
			// Fault-free run: the steps and the finished file.
			dir, path, pre := setupOp(t, op)
			var steps []string
			flushFaultHook = func(s string) error { steps = append(steps, s); return nil }
			err := runOp(t, op, path)
			flushFaultHook = nil
			if err != nil {
				t.Fatalf("fault-free %s: %v", op.name, err)
			}
			post := mustRead(t, path)
			if bytes.Equal(pre, post) {
				t.Fatal("op does not change the file")
			}
			if s := sidecars(t, dir, op.file()); len(s) != 0 {
				t.Fatalf("fault-free run left %v", s)
			}
			t.Logf("steps: %v", steps)
			if len(steps) == 0 {
				t.Errorf("%s passes no crash-safety step (no journal)", op.name)
			}

			for _, rv := range op.recoverers() {
				t.Run("step/"+rv.name, func(t *testing.T) {
					reachedPost := false
					for k, step := range steps {
						dir, path, _ := setupOp(t, op)
						n := 0
						flushFaultHook = func(s string) error {
							if n == k {
								return fmt.Errorf("%w at %s", errInjectedCrash, s)
							}
							n++
							return nil
						}
						err := runOp(t, op, path)
						flushFaultHook = nil
						if !errors.Is(err, errInjectedCrash) {
							t.Fatalf("step %d (%s): %v, want the injected crash", k, step, err)
						}
						checkReadOnlyOpenAfterCrash(t, step, path, pre, post)
						if err := rv.run(path); err != nil {
							t.Errorf("crash after %q: %s: %v", step, rv.name, err)
							continue
						}
						got := mustRead(t, path)
						switch {
						case bytes.Equal(got, post):
							reachedPost = true
						case bytes.Equal(got, pre):
							if reachedPost {
								t.Errorf("crash after %q: committed write rolled back", step)
							}
						default:
							t.Errorf("crash after %q: recovered file is neither the original nor the result (first diff at %d)", step, firstDiff(got, pre))
						}
						if s := sidecars(t, dir, op.file()); len(s) != 0 {
							t.Errorf("crash after %q: left %v", step, s)
						}
					}
					if len(steps) > 0 && !reachedPost {
						t.Errorf("a crash after the last step (%s) lost the completed write", steps[len(steps)-1])
					}
				})
			}

			// A crash inside the data write. Take the journal the call
			// writes before touching the file (none if it writes none).
			var journal []byte
			{
				_, jpath, _ := setupOp(t, op)
				flushFaultHook = func(s string) error {
					if s == "journal-written" {
						return errInjectedCrash
					}
					return nil
				}
				err := runOp(t, op, jpath)
				flushFaultHook = nil
				if errors.Is(err, errInjectedCrash) {
					journal = mustRead(t, journalPath(jpath))
				}
			}
			start, end := writeRange(pre, post)
			torn := 0
			for L := 80; start+L < end; L += 80 {
				state := append([]byte(nil), pre...)
				if n := start + L; n > len(state) {
					state = append(state, make([]byte, n-len(state))...)
				}
				copy(state[start:start+L], post[start:start+L])
				if bytes.Equal(state, pre) || bytes.Equal(state, post) {
					continue
				}
				torn++
				for _, rv := range op.recoverers() {
					dir := t.TempDir()
					p := filepath.Join(dir, op.file())
					if err := os.WriteFile(p, state, 0o644); err != nil {
						t.Fatal(err)
					}
					if journal != nil {
						if err := os.WriteFile(journalPath(p), journal, 0o644); err != nil {
							t.Fatal(err)
						}
					}
					what := fmt.Sprintf("crash with %d of %d bytes written", L, end-start)
					if rf, err := Open(p); err == nil {
						rf.Close()
						t.Fatalf("%s: read-only Open succeeded on a half-written file", what)
					} else if journal != nil && !errors.Is(err, ErrPendingJournal) {
						t.Fatalf("%s: read-only Open: %v, want ErrPendingJournal", what, err)
					}
					if err := rv.run(p); err != nil {
						t.Fatalf("%s: %s: %v", what, rv.name, err)
					}
					if got := mustRead(t, p); !bytes.Equal(got, pre) {
						t.Fatalf("%s: %s recovered a file that is not the original (first diff at %d, len %d want %d)", what, rv.name, firstDiff(got, pre), len(got), len(pre))
					}
					if s := sidecars(t, dir, op.file()); len(s) != 0 {
						t.Fatalf("%s: left %v", what, s)
					}
				}
			}
			if torn == 0 {
				t.Fatal("no torn state modelled")
			}
			t.Logf("%d torn states (write of bytes %d..%d)", torn, start, end)
		})
	}
}

// A real I/O error (not a crash) during a journaled write is rolled back
// before the call returns it: the file is byte-identical to the original,
// no journal or temp file is left, and the *File stays usable (the in-place
// Flush can be retried and then writes the edit; a failed append adds no
// HDU).
func TestInPlaceWriteIOErrorRollsBack(t *testing.T) {
	flushOp := inPlaceOp{name: "Flush", sample: "testprog.std", editFileRecovers: true,
		run: func(t *testing.T, f *File) error {
			inPlaceHDU2(t, f)
			if err := hduHeader(t, f, 3).Set("EXTNAME", "renamed-ascii", ""); err != nil {
				t.Fatal(err)
			}
			return f.Flush()
		}}
	ops := append([]inPlaceOp{flushOp}, inPlaceOps...)
	for _, op := range ops {
		for _, at := range []string{"write", "fsync", "remove-journal"} {
			t.Run(op.name+"/"+at, func(t *testing.T) {
				// The finished file, for the retry check.
				_, ppath, _ := setupOp(t, op)
				if err := runOp(t, op, ppath); err != nil {
					t.Fatal(err)
				}
				post := mustRead(t, ppath)

				dir, path, pre := setupOp(t, op)
				f, err := OpenForEdit(path)
				if err != nil {
					t.Fatal(err)
				}
				defer f.Close()
				nhdu := f.NumHDU()
				injected := errors.New("injected I/O error")
				calls := 0
				ioFaultHook = func(op string) error {
					if op == at {
						calls++
						return injected
					}
					return nil
				}
				err = op.run(t, f)
				ioFaultHook = nil
				if calls == 0 {
					t.Fatalf("%s never reached I/O step %q", op.name, at)
				}
				if !errors.Is(err, injected) {
					t.Fatalf("%s with an I/O error at %q: err = %v, want the injected error", op.name, at, err)
				}
				if got := mustRead(t, path); !bytes.Equal(got, pre) {
					t.Fatalf("I/O error at %q was not rolled back: file differs from the original at byte %d (len %d, want %d)", at, firstDiff(got, pre), len(got), len(pre))
				}
				if s := sidecars(t, dir, op.file()); len(s) != 0 {
					t.Fatalf("I/O error at %q left %v", at, s)
				}
				if f.NumHDU() != nhdu {
					t.Fatalf("failed %s changed NumHDU %d -> %d", op.name, nhdu, f.NumHDU())
				}
				if op.name == "Flush" {
					// The edit is still pending; the retry writes it.
					if err := f.Flush(); err != nil {
						t.Fatalf("retry Flush: %v", err)
					}
					if got := mustRead(t, path); !bytes.Equal(got, post) {
						t.Fatalf("retried Flush wrote a different file (first diff at %d)", firstDiff(got, post))
					}
				}
			})
		}
	}
}

// OverwritePixels writes the zero padding after the data unit (FITS 4.0
// section 3.3.2) also when the file on disk lacks it or holds other bytes
// there. Truth: the cfitsio-written image of each testdata/images case;
// its padding is cut off (or overwritten with 0xAA), the case values are
// written back with OverwritePixels, and the file must equal cfitsio's
// byte for byte (with the documented deviations of expectedCaseData).
func TestOverwritePixelsRestoresPadding(t *testing.T) {
	ran := 0
	for _, c := range loadCases(t) {
		ref := imagesPath(t, "cfwrite", c.Name+".fits")
		want, err := os.ReadFile(ref)
		if err != nil {
			t.Fatalf("missing committed cfitsio file: %v", err)
		}
		r, err := Open(ref)
		if err != nil {
			t.Fatal(err)
		}
		rec := r.hdus[0]
		dataStart, dataEnd, padded := rec.dataStart, rec.dataEnd, rec.paddedEnd
		r.Close()
		if len(r.hdus) != 1 || padded != int64(len(want)) || dataEnd == padded {
			continue // the data fills its last block, or more HDUs follow
		}
		// Cases with documented deviations from cfitsio's rounding carry
		// them in the data (expectedCaseData); header and padding are
		// cfitsio's file.
		want = append([]byte(nil), want...)
		copy(want[dataStart:dataEnd], expectedCaseData(t, c))
		for _, damage := range []string{"cut", "nonzero"} {
			t.Run(c.Name+"/"+damage, func(t *testing.T) {
				ran++
				start := append([]byte(nil), want[:dataEnd]...)
				if damage == "nonzero" {
					start = append(start, bytes.Repeat([]byte{0xAA}, int(padded-dataEnd))...)
				}
				path := filepath.Join(t.TempDir(), c.Name+".fits")
				if err := os.WriteFile(path, start, 0o644); err != nil {
					t.Fatal(err)
				}
				f, err := OpenForEdit(path)
				if err != nil {
					t.Fatal(err)
				}
				img, err := f.Primary()
				if err != nil {
					t.Fatal(err)
				}
				_, over := caseData(t, c, false)
				if err := over(img); err != nil {
					t.Fatalf("OverwritePixels: %v", err)
				}
				if err := f.Close(); err != nil {
					t.Fatal(err)
				}
				if got := mustRead(t, path); !bytes.Equal(got, want) {
					t.Fatalf("file differs from cfitsio's at byte %d (len %d, want %d)", firstDiff(got, want), len(got), len(want))
				}
			})
		}
	}
	if ran == 0 {
		t.Fatal("no case with a partly filled last block")
	}
}
