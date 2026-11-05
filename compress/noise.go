package compress

import (
	"math"
	"sort"
)

// Image statistics for float32 tile quantization.
//
// Pure-Go port of cfitsio's fits_img_stats_float / FnMeanSigma_float /
// FnNoise1_float / FnNoise5_float from quantize.c. The function names
// and structure deliberately mirror the C so that future drift can be
// diffed line-by-line against reference/cfitsio/quantize.c.
//
// The noise estimator is the key input to fits_quantize_float: the
// quantization step delta is stdev/qlevel, where stdev is the minimum
// of the 2nd/3rd/5th-order MAD noise estimates. The 3rd-order estimate
// (FnNoise5's noise3 output) is the primary, with the 2nd and 5th
// providing sanity bounds in the presence of background structure.
//
// Reference:
//   reference/cfitsio/quantize.c:579  fits_img_stats_float
//   reference/cfitsio/quantize.c:1577 FnNoise5_float
//   reference/cfitsio/quantize.c:3373 FnNoise1_float
//   reference/cfitsio/quantize.c:3679 quick_select_float

// ImageStatsF32 carries the per-tile statistics that fits_img_stats_float
// fills in. Any field may be zero if the corresponding path in the C
// function wasn't taken (e.g. sigma-clipping iterations bailing out).
type ImageStatsF32 struct {
	NGood  int64   // number of non-null pixels
	Min    float32 // minimum non-null value
	Max    float32 // maximum non-null value
	Mean   float64 // mean of non-null pixels
	Sigma  float64 // RMS sigma of non-null pixels
	Noise1 float64 // 1st-order differences, sigma-clipped per row, median of rows
	Noise2 float64 // 2nd-order MAD of non-null pixels
	Noise3 float64 // 3rd-order MAD (primary noise estimator for quantization)
	Noise5 float64 // 5th-order MAD
}

// ImgStatsFloat32 is the Go counterpart of cfitsio's fits_img_stats_float.
// array is a 2D image of shape (nx, ny) laid out in row-major order,
// first axis varying fastest. When nullcheck is true, pixels equal to
// nullvalue are treated as null and excluded from every statistic.
//
// This function matches cfitsio's numerical output bit-for-bit on the
// fixtures in testdata/cref/golden/.
func ImgStatsFloat32(array []float32, nx, ny int64, nullcheck bool, nullvalue float32) ImageStatsF32 {
	var s ImageStatsF32

	// Mean + sigma on the flat pixel list (FnMeanSigma_float).
	s.NGood, s.Mean, s.Sigma = fnMeanSigmaFloat32(array, nx*ny, nullcheck, nullvalue)

	// 1st-order differences, sigma-clipped per row, median of rows
	// (FnNoise1_float).
	s.Noise1 = fnNoise1Float32(array, nx, ny, nullcheck, nullvalue)

	// 2nd/3rd/5th-order MAD + min/max (FnNoise5_float). This overrides
	// NGood with a possibly different count because FnNoise5 increments
	// ngoodpix as it walks the rows.
	ng, minv, maxv, n2, n3, n5 := fnNoise5Float32(array, nx, ny, nullcheck, nullvalue)
	s.NGood = ng
	s.Min = minv
	s.Max = maxv
	s.Noise2 = n2
	s.Noise3 = n3
	s.Noise5 = n5
	return s
}

// fnMeanSigmaFloat32 is the port of FnMeanSigma_float.
// Plain running sum and sum-of-squares. For bit-exact parity with cfitsio,
// the iteration order and the intermediate double accumulations must be
// identical.
func fnMeanSigmaFloat32(array []float32, npix int64, nullcheck bool, nullvalue float32) (ngood int64, mean, sigma float64) {
	var sum, sum2 float64
	if nullcheck {
		for i := int64(0); i < npix; i++ {
			if array[i] != nullvalue {
				ngood++
				x := float64(array[i])
				sum += x
				sum2 += float64(x * x)
			}
		}
	} else {
		ngood = npix
		for i := int64(0); i < npix; i++ {
			x := float64(array[i])
			sum += x
			sum2 += float64(x * x)
		}
	}
	if ngood > 1 {
		mean = sum / float64(ngood)
		sigma = math.Sqrt(sum2/float64(ngood) - float64(mean*mean))
	} else if ngood == 1 {
		mean = sum
		sigma = 0
	}
	return
}

// fnMeanSigmaFloat32Slice is fnMeanSigmaFloat32 without the nullcheck
// fast-path, used by the inner sigma-clipping loop of fnNoise1.
func fnMeanSigmaFloat32Slice(arr []float32) (mean, stdev float64) {
	if len(arr) == 0 {
		return 0, 0
	}
	var sum, sum2 float64
	for _, v := range arr {
		x := float64(v)
		sum += x
		sum2 += float64(x * x)
	}
	n := float64(len(arr))
	mean = sum / n
	// As FnMeanSigma_float: no clamp of a negative variance (NaN from
	// sqrt, as in C). For differences of float32 pixels it cannot occur:
	// the sums are exact and the rounding of mean*mean and sum2/n is
	// below the smallest non-zero variance.
	return mean, math.Sqrt(sum2/n - float64(mean*mean))
}

// fnNoise1Float32 ports FnNoise1_float.
//
// For each row: build first-difference array diff[i] = rowpix[i] -
// rowpix[i-1] over non-null pixels. Compute mean and stdev; iteratively
// sigma-clip (NITER=3, SIGMA_CLIP=5). The row's noise contribution is
// the resulting stdev. Across all rows, take the median (lower-middle
// for even row count), multiply by 1/sqrt(2) ~= 0.70710678.
func fnNoise1Float32(array []float32, nx, ny int64, nullcheck bool, nullvalue float32) float64 {
	if nx < 3 {
		return 0
	}
	const niter = 3
	const sigmaClip = 5.0
	differences := make([]float32, nx)
	diffs := make([]float64, 0, ny)
	for jj := int64(0); jj < ny; jj++ {
		rowpix := array[jj*nx : (jj+1)*nx]
		ii := int64(0)
		if nullcheck {
			for ii < nx && rowpix[ii] == nullvalue {
				ii++
			}
		}
		if ii == nx {
			continue
		}
		v1 := rowpix[ii]
		nvals := 0
		for ii++; ii < nx; ii++ {
			if nullcheck {
				for ii < nx && rowpix[ii] == nullvalue {
					ii++
				}
			}
			if ii == nx {
				break
			}
			differences[nvals] = v1 - rowpix[ii]
			nvals++
			v1 = rowpix[ii]
		}
		if nvals < 2 {
			continue
		}
		mean, stdev := fnMeanSigmaFloat32Slice(differences[:nvals])
		if stdev > 0 {
			for iter := 0; iter < niter; iter++ {
				kk := 0
				for i := 0; i < nvals; i++ {
					if math.Abs(float64(differences[i])-mean) < float64(sigmaClip*stdev) {
						if kk < i {
							differences[kk] = differences[i]
						}
						kk++
					}
				}
				if kk == nvals {
					break
				}
				nvals = kk
				mean, stdev = fnMeanSigmaFloat32Slice(differences[:nvals])
			}
		}
		diffs = append(diffs, stdev)
	}
	xnoise := medianDoubleLowerMiddle(diffs)
	return float64(0.70710678 * xnoise) // 1/sqrt(2): normalize first-difference variance to per-pixel sigma
}

// fnNoise5Float32 ports FnNoise5_float.
//
// Returns ngood, minval, maxval, noise2, noise3, noise5 in the same
// units / scaling as the C reference:
//
//	noise2 = 1.0483579 * median(|v5 - v7|)
//	noise3 = 0.6052697 * median(|2*v5 - v3 - v7|)
//	noise5 = 0.1772048 * median(|6*v5 - 4*v3 - 4*v7 + v1 + v9|)
//
// where the inner median is the lower-middle quickselect result per
// row, and the outer median across rows is the average of the two
// middles for even row counts (standard median).
func fnNoise5Float32(array []float32, nx, ny int64, nullcheck bool, nullvalue float32) (ngood int64, minval, maxval float32, noise2, noise3, noise5 float64) {
	return fnNoise5(array, nx, ny, nullcheck, nullvalue)
}

// fnNoise5 is FnNoise5_float (T float32) and FnNoise5_double (T float64):
// the two C functions are the same code with the pixel type (and the
// type of the difference arrays) changed, so the port is one generic
// function. All arithmetic on pixels is done in T, as in C.
func fnNoise5[T floatPix](array []T, nx, ny int64, nullcheck bool, nullvalue T) (ngood int64, minval, maxval T, noise2, noise3, noise5 float64) {
	// rows must have at least 9 pixels — otherwise flatten to a
	// single row and try again.
	if nx < 9 {
		nx = nx * ny
		ny = 1
	}

	// Fallback: even after flattening, too few pixels for any MAD.
	if nx < 9 {
		xminval, xmaxval := maxPix[T](), -maxPix[T]()
		var ngoodpix int64
		for ii := int64(0); ii < nx; ii++ {
			if nullcheck && array[ii] == nullvalue {
				continue
			}
			if array[ii] < xminval {
				xminval = array[ii]
			}
			if array[ii] > xmaxval {
				xmaxval = array[ii]
			}
			ngoodpix++
		}
		return ngoodpix, xminval, xmaxval, 0, 0, 0
	}

	xminval := maxPix[T]()
	xmaxval := -maxPix[T]()

	differences2 := make([]T, nx)
	differences3 := make([]T, nx)
	differences5 := make([]T, nx)
	diffs2 := make([]float64, 0, ny)
	diffs3 := make([]float64, 0, ny)
	diffs5 := make([]float64, 0, ny)
	var ngoodpix int64

	for jj := int64(0); jj < ny; jj++ {
		rowpix := array[jj*nx : (jj+1)*nx]

		// Find the first 8 valid pixels. Any row that runs out of
		// valid pixels before 8 is skipped entirely.
		var v1, v2, v3, v4, v5, v6, v7, v8, v9 T
		ii := int64(0)
		initial := [8]*T{&v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8}
		skipRow := false
		for slot := 0; slot < 8; slot++ {
			if nullcheck {
				for ii < nx && rowpix[ii] == nullvalue {
					ii++
				}
			}
			if ii == nx {
				skipRow = true
				break
			}
			*initial[slot] = rowpix[ii]
			ngoodpix++
			if rowpix[ii] < xminval {
				xminval = rowpix[ii]
			}
			if rowpix[ii] > xmaxval {
				xmaxval = rowpix[ii]
			}
			ii++
		}
		if skipRow {
			continue
		}

		// Walk the rest of the row, shifting the v1..v9 window by
		// one valid pixel at a time and emitting a difference each step.
		nvals := 0
		nvals2 := 0
		for ; ii < nx; ii++ {
			if nullcheck {
				for ii < nx && rowpix[ii] == nullvalue {
					ii++
				}
			}
			if ii == nx {
				break
			}
			v9 = rowpix[ii]
			if v9 < xminval {
				xminval = v9
			}
			if v9 > xmaxval {
				xmaxval = v9
			}

			// 2nd-order: skip identical-run regions because |v5-v7|=0
			// there.
			if !(v5 == v6 && v6 == v7) {
				differences2[nvals2] = absPix(v5 - v7)
				nvals2++
			}

			// 3rd/5th-order: skip regions where v3..v7 are all equal.
			if !(v3 == v4 && v4 == v5 && v5 == v6 && v6 == v7) {
				differences3[nvals] = absPix(T(2*v5) - v3 - v7)
				differences5[nvals] = absPix(T(T(T(T(6*v5)-T(4*v3))-T(4*v7))+v1) + v9)
				nvals++
			} else {
				// still counts as "good" for ngood bookkeeping
				ngoodpix++
			}

			v1, v2, v3, v4, v5, v6, v7, v8 = v2, v3, v4, v5, v6, v7, v8, v9
		}

		ngoodpix += int64(nvals)

		if nvals == 0 {
			continue
		}
		if nvals == 1 {
			if nvals2 == 1 {
				diffs2 = append(diffs2, float64(differences2[0]))
			}
			diffs3 = append(diffs3, float64(differences3[0]))
			diffs5 = append(diffs5, float64(differences5[0]))
		} else {
			if nvals2 > 1 {
				// cfitsio passes nvals, not nvals2: when a run of 3 equal
				// pixels skipped a 2nd-order difference, the median is
				// taken over nvals entries, including zeros or entries
				// left (in quick_select's order) by the previous row.
				// Kept for bit-exactness (goldens runs3_*).
				diffs2 = append(diffs2, float64(quickSelect(differences2[:nvals])))
			}
			diffs3 = append(diffs3, float64(quickSelect(differences3[:nvals])))
			diffs5 = append(diffs5, float64(quickSelect(differences5[:nvals])))
		}
	}

	var xnoise2, xnoise3, xnoise5 float64
	xnoise3 = medianDoublePair(diffs3)
	xnoise5 = medianDoublePair(diffs5)
	xnoise2 = medianDoublePair(diffs2)

	noise2 = float64(1.0483579 * xnoise2)
	noise3 = float64(0.6052697 * xnoise3)
	noise5 = float64(0.1772048 * xnoise5)
	return ngoodpix, xminval, xmaxval, noise2, noise3, noise5
}

// fnNoise3Float32MinMax ports FnNoise3_float in its noise=NULL shape,
// used by fits_quantize_float's negative-qlevel branch. We only need
// ngood, minval, maxval — the 3rd-order diff array is already covered
// by fnNoise5Float32 for the positive-qlevel path. The bookkeeping is
// preserved verbatim because fits_quantize_float compares `ngood == nx`
// to decide the zeropt calculation.
//
// Reference: reference/cfitsio/quantize.c:2703 (FnNoise3_float).
func fnNoise3Float32MinMax(array []float32, nx, ny int64, nullcheck bool, nullvalue float32) (ngood int64, minval, maxval float32) {
	return fnNoise3MinMax(array, nx, ny, nullcheck, nullvalue)
}

// fnNoise3MinMax is FnNoise3_float / FnNoise3_double (noise = NULL), generic
// over the pixel type like fnNoise5.
func fnNoise3MinMax[T floatPix](array []T, nx, ny int64, nullcheck bool, nullvalue T) (ngood int64, minval, maxval T) {
	if nx < 5 {
		nx = nx * ny
		ny = 1
	}
	if nx < 5 {
		xminval, xmaxval := maxPix[T](), -maxPix[T]()
		var ngoodpix int64
		for ii := int64(0); ii < nx; ii++ {
			if nullcheck && array[ii] == nullvalue {
				continue
			}
			if array[ii] < xminval {
				xminval = array[ii]
			}
			if array[ii] > xmaxval {
				xmaxval = array[ii]
			}
			ngoodpix++
		}
		return ngoodpix, xminval, xmaxval
	}

	xminval := maxPix[T]()
	xmaxval := -maxPix[T]()
	var ngoodpix int64

	for jj := int64(0); jj < ny; jj++ {
		rowpix := array[jj*nx : (jj+1)*nx]
		var v5 T
		ii := int64(0)
		// Need v1..v4 before the main loop over v5 can start. We
		// just read and min/max-update them; when noise=nil nothing
		// else uses them, but skipping short rows is still required
		// so that the ngoodpix bookkeeping matches cfitsio's exactly.
		skipRow := false
		for slot := 0; slot < 4; slot++ {
			if nullcheck {
				for ii < nx && rowpix[ii] == nullvalue {
					ii++
				}
			}
			if ii == nx {
				skipRow = true
				break
			}
			v := rowpix[ii]
			if v < xminval {
				xminval = v
			}
			if v > xmaxval {
				xmaxval = v
			}
			ii++
		}
		if skipRow {
			continue
		}
		for ; ii < nx; ii++ {
			if nullcheck {
				for ii < nx && rowpix[ii] == nullvalue {
					ii++
				}
			}
			if ii == nx {
				break
			}
			v5 = rowpix[ii]
			if v5 < xminval {
				xminval = v5
			}
			if v5 > xmaxval {
				xmaxval = v5
			}
			ngoodpix++ // noise==nil branch increments per v5
		}
		// Add 4 for v1..v4 (matches `ngoodpix += nvals + 4` with nvals=0).
		ngoodpix += 4
	}
	return ngoodpix, xminval, xmaxval
}

// floatPix is the pixel type of the float quantizer: float32
// (fits_quantize_float) or float64 (fits_quantize_double).
type floatPix interface{ ~float32 | ~float64 }

// maxPix is FLT_MAX for float32 and DBL_MAX for float64.
func maxPix[T floatPix]() T {
	var z T
	if any(z) == any(float32(0)) {
		return T(math.MaxFloat32)
	}
	m := math.MaxFloat64
	return T(m)
}

// absPix mirrors (float)fabs(x) / fabs(x): the absolute value of a T is
// exact, so the double detour of the float32 version changes nothing.
func absPix[T floatPix](x T) T {
	if x < 0 {
		return -x
	}
	return x
}

// quickSelect is cfitsio's quick_select_float / quick_select_double (quantize.c), a
// verbatim port: it returns the element at index (len(arr)-1)/2 of the
// sorted order and leaves arr in the same partially-partitioned order as
// the C code. The order matters: FnNoise5 reads entries a previous row's
// selection left behind (see the nvals note above), so a different
// permutation (e.g. a full sort) changes noise2.
func quickSelect[T floatPix](arr []T) T {
	n := len(arr)
	if n == 0 {
		return 0
	}
	low, high := 0, n-1
	median := (low + high) / 2
	for {
		if high <= low { // one element only
			return arr[median]
		}
		if high == low+1 { // two elements only
			if arr[low] > arr[high] {
				arr[low], arr[high] = arr[high], arr[low]
			}
			return arr[median]
		}
		// Find median of low, middle and high items; swap into position low.
		middle := (low + high) / 2
		if arr[middle] > arr[high] {
			arr[middle], arr[high] = arr[high], arr[middle]
		}
		if arr[low] > arr[high] {
			arr[low], arr[high] = arr[high], arr[low]
		}
		if arr[middle] > arr[low] {
			arr[middle], arr[low] = arr[low], arr[middle]
		}
		// Swap low item (now in position middle) into position low+1.
		arr[middle], arr[low+1] = arr[low+1], arr[middle]
		// Nibble from each end towards middle, swapping items when stuck.
		ll, hh := low+1, high
		for {
			ll++
			for arr[low] > arr[ll] {
				ll++
			}
			hh--
			for arr[hh] > arr[low] {
				hh--
			}
			if hh < ll {
				break
			}
			arr[ll], arr[hh] = arr[hh], arr[ll]
		}
		// Swap middle item (in position low) back into correct position.
		arr[low], arr[hh] = arr[hh], arr[low]
		// Re-set active partition.
		if hh <= median {
			low = ll
		}
		if hh >= median {
			high = hh - 1
		}
	}
}

// medianDoublePair ports the "compute median of row medians" logic
// used at the tail of FnNoise5 / FnNoise1:
//
//	if nrows == 0         -> 0
//	if nrows == 1         -> diffs[0]
//	else                  -> (diffs_sorted[(n-1)/2] + diffs_sorted[n/2]) / 2
func medianDoublePair(diffs []float64) float64 {
	if len(diffs) == 0 {
		return 0
	}
	if len(diffs) == 1 {
		return diffs[0]
	}
	sort.Float64s(diffs)
	return (diffs[(len(diffs)-1)/2] + diffs[len(diffs)/2]) / 2
}

// medianDoubleLowerMiddle is the noise1 variant — FnNoise1 also uses
// the paired-middle form after sorting diffs. Named separately to keep
// the reference mapping explicit, but delegates to medianDoublePair.
func medianDoubleLowerMiddle(diffs []float64) float64 {
	return medianDoublePair(diffs)
}
