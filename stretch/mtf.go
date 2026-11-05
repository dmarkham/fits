package stretch

import (
	"errors"
	"fmt"
	"math"

	"github.com/dmarkham/fits/stats"
)

// MTFParams is one set of Midtone Transfer Function parameters, as Siril's
// autostretch computes them (struct mtf_params in filters/mtf.h).
type MTFParams struct {
	Shadows    float32 // black clip point: pixels <= Shadows become 0
	Midtones   float32 // midtone balance m (0.5 = identity)
	Highlights float32 // white clip point: pixels >= Highlights become 1
	// StatsOK is false when a channel had no usable pixel (all zero, NaN
	// or ±Inf). Siril then falls back to (0, 0.2, 1) for every channel and
	// still applies it; so do the Autostretch functions.
	StatsOK bool
}

// MTF applies the Midtone Transfer Function to a single channel, exactly
// as Siril's MTF() (filters/mtf.c) does it in float32:
//
//	x <= shadows    -> 0
//	x >= highlights -> 1
//	otherwise       x' = (x - shadows) / (highlights - shadows)
//	                out = ((m - 1) * x') / ((2m - 1) * x' - m)
//
// midtones = 0.5 is the identity on (shadows, highlights) after
// rescaling. NaN pixels stay NaN.
//
// MTF takes every finite shadows and highlights and 0 <= midtones <= 1:
// these include values siril-cli's `mtf` command refuses (midtones 0 or 1,
// shadows outside [0, 1), highlights outside (0, 1]) but Siril's
// autostretch hands its MTF (midtones 0 or 1, shadows above 1, highlights
// below 0), and MTF gives Siril's MTF() result for them. A NaN or infinite
// parameter, or midtones outside [0, 1] (where the denominator
// (2m-1)x' - m can be 0 inside the range), is an error wrapping
// ErrMTFParams and the pixels are left unchanged: Siril would return NaN
// or infinite pixels there.
func MTF(pixels []float32, shadows, midtones, highlights float32) error {
	if !(midtones >= 0 && midtones <= 1) || !finite32(shadows) || !finite32(highlights) {
		return fmt.Errorf("%w: shadows %v, midtones %v, highlights %v; want finite shadows and highlights and 0 <= midtones <= 1",
			ErrMTFParams, shadows, midtones, highlights)
	}
	applyMTF(pixels, shadows, midtones, highlights)
	return nil
}

// ErrMTFParams is wrapped by the error MTF returns for parameters it
// refuses.
var ErrMTFParams = errors.New("stretch: MTF parameters out of range")

func finite32(v float32) bool { return v >= -math.MaxFloat32 && v <= math.MaxFloat32 }

// applyMTF is MTF without the parameter check: the Autostretch functions
// apply Siril's parameters whatever they are (a ±Inf shadows clip gives
// infinite or NaN shadows, as in siril-cli).
func applyMTF(pixels []float32, shadows, midtones, highlights float32) {
	for i, v := range pixels {
		pixels[i] = mtf(v, midtones, shadows, highlights)
	}
}

// mtf is Siril's float MTF(float x, float m, float lo, float hi).
func mtf(x, m, lo, hi float32) float32 {
	if x <= lo {
		return 0
	}
	if x >= hi {
		return 1
	}
	xp := (x - lo) / (hi - lo)
	return float32((m-1)*xp) / (float32((float32(2*m)-1)*xp) - m)
}

// MTFValue evaluates the MTF at a single point with shadows = 0 and
// highlights = 1: 0 for x <= 0, 1 for x >= 1, the formula otherwise.
func MTFValue(x, m float32) float32 { return mtf(x, m, 0, 1) }

// sirilChannelStats returns the median and MAD Siril's autostretch uses
// for one channel (statistics_internal_float with STATS_BASIC|STATS_MAD):
// both are taken over the pixels that are neither 0 nor NaN (negative
// pixels count); the median is the histogram median of
// findMinMaxPercentile (stats.Median on float32), and the MAD is that
// median of |x - (float)median| computed in float32
// (siril_stats_float_mad). ok is false when no pixel qualifies.
//
// ±Inf pixels are also left out. That is a deliberate difference from
// Siril, whose findMinMaxPercentile casts the infinite or NaN histogram
// position of an infinite pixel (or of every pixel, when the minimum or
// maximum is infinite) to uint16, which is undefined behaviour in C++; on
// x86-64 one Inf pixel turns that channel's parameters into NaN or
// meaningless values (0, 0, 1) or (0, 1, 1) (TestInf_SirilOracle records
// Siril's values). With the exclusion the parameters equal Siril's for the
// same image with the Inf pixels set to 0.
func sirilChannelStats(ch []float32) (median, mad float64, ok bool) {
	good := make([]float32, 0, len(ch))
	for _, v := range ch {
		if v != 0 && v == v && !math.IsInf(float64(v), 0) {
			good = append(good, v)
		}
	}
	if len(good) == 0 {
		return 0, 0, false
	}
	median = stats.Median(good)
	m := float32(median)
	for i, v := range good {
		good[i] = float32(math.Abs(float64(v - m)))
	}
	return median, stats.Median(good), true
}

// madNorm is Siril's MAD_NORM (core/siril.h:104), used as (float)MAD_NORM.
const madNorm = float32(1.4826)

type chanStat struct {
	median, mad float32 // as find_*_midtones_balance converts them
	inverted    bool    // median > 0.5 (compared in double)
}

// channelStats mirrors the statistics block shared by Siril's
// find_unlinked_midtones_balance and find_linked_midtones_balance
// (filters/mtf.c:259, :374) for a float image (normValue = 1). ok is false
// if any channel has no usable pixel (Siril: retval -1).
func channelStats(channels [][]float32) (st []chanStat, inverted int, ok bool) {
	ok = true
	for _, ch := range channels {
		med, mad, cok := sirilChannelStats(ch)
		if !cok {
			ok = false
		}
		const normValue = float32(1)
		s := chanStat{
			median:   float32(med) / normValue,
			mad:      float32(float32(float32(mad)/normValue) * madNorm),
			inverted: med/1.0 > 0.5,
		}
		// guard against the breakdown point
		if s.mad == 0 {
			s.mad = 0.001
		}
		if s.inverted {
			inverted++
		}
		st = append(st, s)
	}
	return st, inverted, ok
}

// AutostretchUnlinkedParams computes Siril's unlinked autostretch
// parameters for an image with the given channels (1 for mono, 3 for RGB),
// porting find_unlinked_midtones_balance (filters/mtf.c:374):
//
// If at least one channel has median <= 0.5 (normal case), each channel c
// gets shadows = max(0, median_c + shadowsClip*MADn_c), highlights = 1,
// midtones = MTF(median_c - shadows, targetBG). If EVERY channel has
// median > 0.5 (bright image), each channel instead gets shadows = 0,
// highlights = min(1, median_c - shadowsClip*MADn_c) and midtones =
// 1 - MTF(highlights - median_c, targetBG). MADn is the MAD times 1.4826
// (0.001 if the MAD is 0); median and MAD exclude pixels that are 0 or
// NaN, as in Siril, and also ±Inf, where Siril's result is undefined (see
// sirilChannelStats); negative pixels count. "median > 0.5" is strict: a
// median of exactly 0.5 takes the normal branch. If any channel has no
// usable pixel, every channel gets (0, 0.2, 1) with StatsOK = false.
//
// shadowsClip is normally negative (e.g. -2.8); targetBG is the target
// background level in [0, 1] (e.g. 0.25). It panics with an error wrapping
// ErrAutostretchParams if targetBG is outside [0, 1] or NaN, or shadowsClip
// is NaN (AutostretchUnlinked returns that error instead).
func AutostretchUnlinkedParams(channels [][]float32, shadowsClip, targetBG float32) []MTFParams {
	mustAutostretchArgs(shadowsClip, targetBG)
	st, inverted, ok := channelStats(channels)
	out := make([]MTFParams, len(channels))
	if !ok {
		for i := range out {
			out[i] = MTFParams{Shadows: 0, Midtones: 0.2, Highlights: 1}
		}
		return out
	}
	for i, s := range st {
		if inverted < len(channels) {
			c0 := s.median + float32(shadowsClip*s.mad)
			if c0 < 0 {
				c0 = 0
			}
			m2 := s.median - c0
			out[i] = MTFParams{Shadows: c0, Midtones: mtf(m2, targetBG, 0, 1), Highlights: 1, StatsOK: true}
		} else {
			c1 := s.median - float32(shadowsClip*s.mad)
			if c1 > 1 {
				c1 = 1
			}
			m2 := c1 - s.median
			out[i] = MTFParams{Shadows: 0, Midtones: 1 - mtf(m2, targetBG, 0, 1), Highlights: c1, StatsOK: true}
		}
	}
	return out
}

// AutostretchParams computes Siril's autostretch parameters for a mono
// image (AutostretchUnlinkedParams with one channel): the normal branch
// when the median of the usable pixels (not 0, NaN or ±Inf) is <= 0.5,
// the inverted branch (shadows 0, highlights < 1) when it is above 0.5.
// If the channel has no usable pixel the result is (0, 0.2, 1). It panics
// on the arguments Autostretch refuses (ErrAutostretchParams).
func AutostretchParams(pixels []float32, shadowsClip, targetBG float32) (shadows, midtones, highlights float32) {
	p := AutostretchUnlinkedParams([][]float32{pixels}, shadowsClip, targetBG)[0]
	return p.Shadows, p.Midtones, p.Highlights
}

// AutostretchLinkedParams computes one set of MTF parameters for all
// channels from their averaged statistics, porting Siril's
// find_linked_midtones_balance (filters/mtf.c:259): with n channels,
//
//	normal (some channel median <= 0.5):
//	  shadows = max(0, Σ(median_c + shadowsClip*MADn_c)/n), highlights = 1,
//	  midtones = MTF(Σmedian_c/n - shadows, targetBG)
//	inverted (every channel median > 0.5):
//	  shadows = 0, highlights = min(1, Σ(median_c - shadowsClip*MADn_c)/n),
//	  midtones = 1 - MTF(highlights - Σmedian_c/n, targetBG)
//
// with the sums accumulated in float32 in channel order and divided by n
// in float32, as Siril does. If a channel has no usable pixel (not 0, NaN
// or ±Inf) the result is (0, 0.2, 1). It panics on the arguments
// AutostretchLinked refuses (ErrAutostretchParams).
func AutostretchLinkedParams(channels [][]float32, shadowsClip, targetBG float32) (shadows, midtones, highlights float32) {
	mustAutostretchArgs(shadowsClip, targetBG)
	st, inverted, ok := channelStats(channels)
	if !ok {
		return 0, 0.2, 1
	}
	n := float32(len(channels))
	var c0, c1, m float32
	if inverted < len(channels) {
		for _, s := range st {
			c0 += s.median + float32(shadowsClip*s.mad)
			m += s.median
		}
		c0 /= n
		if c0 < 0 {
			c0 = 0
		}
		m2 := m/n - c0
		return c0, mtf(m2, targetBG, 0, 1), 1
	}
	for _, s := range st {
		m += s.median
		c1 += s.median - float32(shadowsClip*s.mad)
	}
	c1 /= n
	if c1 > 1 {
		c1 = 1
	}
	m2 := c1 - m/n
	return 0, 1 - mtf(m2, targetBG, 0, 1), c1
}

// ErrAutostretchParams is wrapped by the error the Autostretch functions
// return (and the value the *Params functions panic with) for arguments
// outside the accepted range: targetBG outside [0, 1] or NaN (siril-cli's
// process_autostretch refuses targetBG outside [0, 1], core/command.c:3338),
// or a NaN shadowsClip. siril-cli accepts a NaN for either and returns an
// image of NaN; that is refused here. ±Inf shadowsClip is accepted and
// gives Siril's result (shadows clamp to 0, highlights to 1).
var ErrAutostretchParams = errors.New("stretch: autostretch parameters out of range")

func autostretchArgs(shadowsClip, targetBG float32) error {
	if shadowsClip != shadowsClip {
		return fmt.Errorf("%w: shadowsClip is NaN", ErrAutostretchParams)
	}
	if !(targetBG >= 0 && targetBG <= 1) {
		return fmt.Errorf("%w: targetBG = %v, want 0 <= targetBG <= 1", ErrAutostretchParams, targetBG)
	}
	return nil
}

// mustAutostretchArgs: the *Params functions keep their (error-free)
// signatures and panic on arguments the Autostretch functions refuse.
func mustAutostretchArgs(shadowsClip, targetBG float32) {
	if err := autostretchArgs(shadowsClip, targetBG); err != nil {
		panic(err)
	}
}

// Autostretch applies Siril's autostretch to a mono image: compute the
// parameters with AutostretchParams, then apply MTF. For an RGB image use
// AutostretchUnlinked (Siril's `autostretch`), not Autostretch per
// channel: Siril only takes the inverted branch when every channel is
// bright. Out-of-range arguments are an error wrapping
// ErrAutostretchParams and leave the pixels unchanged.
func Autostretch(pixels []float32, shadowsClip, targetBG float32) error {
	if err := autostretchArgs(shadowsClip, targetBG); err != nil {
		return err
	}
	shadows, midtones, highlights := AutostretchParams(pixels, shadowsClip, targetBG)
	applyMTF(pixels, shadows, midtones, highlights)
	return nil
}

// AutostretchUnlinked applies Siril's `autostretch shadowsClip targetBG`
// to an image with the given channels: per-channel parameters from
// AutostretchUnlinkedParams, each applied to its own channel.
// Out-of-range arguments are an error wrapping ErrAutostretchParams and
// leave the pixels unchanged.
func AutostretchUnlinked(channels [][]float32, shadowsClip, targetBG float32) error {
	if err := autostretchArgs(shadowsClip, targetBG); err != nil {
		return err
	}
	for i, p := range AutostretchUnlinkedParams(channels, shadowsClip, targetBG) {
		applyMTF(channels[i], p.Shadows, p.Midtones, p.Highlights)
	}
	return nil
}

// AutostretchLinked applies Siril's `autostretch -linked shadowsClip
// targetBG`: one parameter set from AutostretchLinkedParams, applied to
// every channel. Out-of-range arguments are an error wrapping
// ErrAutostretchParams and leave the pixels unchanged.
func AutostretchLinked(channels [][]float32, shadowsClip, targetBG float32) error {
	if err := autostretchArgs(shadowsClip, targetBG); err != nil {
		return err
	}
	shadows, midtones, highlights := AutostretchLinkedParams(channels, shadowsClip, targetBG)
	for _, ch := range channels {
		applyMTF(ch, shadows, midtones, highlights)
	}
	return nil
}
