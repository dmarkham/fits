package stretch

import (
	"errors"
	"testing"
)

// TestAsinhRGBChannelLengthMismatch: mismatched channel slices are an
// error and leave every pixel unchanged (previously the extra pixels of a
// longer g or b were silently left unstretched, and a shorter one
// panicked).
func TestAsinhRGBChannelLengthMismatch(t *testing.T) {
	for _, fn := range []struct {
		name string
		f    func(r, g, b []float32, beta, offset float32) error
	}{{"AsinhRGB", AsinhRGB}, {"AsinhRGBHuman", AsinhRGBHuman}} {
		for _, c := range []struct{ nr, ng, nb int }{{3, 4, 3}, {3, 3, 2}, {4, 3, 3}} {
			r, g, b := make([]float32, c.nr), make([]float32, c.ng), make([]float32, c.nb)
			for _, s := range [][]float32{r, g, b} {
				for i := range s {
					s[i] = 0.25
				}
			}
			err := fn.f(r, g, b, 10, 0)
			if !errors.Is(err, ErrAsinhParams) {
				t.Errorf("%s lengths %v: err = %v, want ErrAsinhParams", fn.name, c, err)
			}
			for _, s := range [][]float32{r, g, b} {
				for i, v := range s {
					if v != 0.25 {
						t.Fatalf("%s lengths %v: pixel %d changed to %v", fn.name, c, i, v)
					}
				}
			}
		}
	}
}
