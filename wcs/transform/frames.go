package transform

import "math"

// Sky coordinate-frame conversions.
//
// This file implements the four frame conversions most commonly needed
// for astronomical image analysis:
//
//   - Equatorial ↔ Galactic (fixed rotation)
//   - Equatorial ↔ Ecliptic (rotation by the time-dependent obliquity
//     of the ecliptic)
//   - Equatorial ↔ Supergalactic (fixed rotation, via galactic)
//   - ICRS / FK5 / FK4 precession (between equatorial frames at
//     different equinoxes)
//
// All inputs and outputs are in DEGREES. Output longitudes are normalized
// to [0, 360); latitudes are in [-90, 90].
//
// The conversions use standard 3-vector rotations: convert (lon, lat) to
// a unit 3-vector, multiply by a fixed matrix, convert back. For
// time-dependent frames (precession, ecliptic obliquity) the matrix is
// computed from the given epoch.

// ------------------ Equatorial ↔ Galactic ------------------

// eqToGalMatrix rotates equatorial (RA, Dec) in ICRS → galactic (l, b).
// From Blaauw et al. 1960, expressed in the J2000 frame.
//
// The column layout means that when applied to the equatorial unit
// vector (1, 0, 0) (RA=0, Dec=0, the vernal equinox) the result is the
// galactic coordinates of the vernal equinox, approximately
// (l=96.3°, b=-60.2°). Equivalently, the rows of this matrix are the
// three galactic basis vectors expressed in equatorial coordinates.
var eqToGalMatrix = [3][3]float64{
	{-0.054875539396, -0.873437104728, -0.483834991775},
	{+0.494109453628, -0.444829594298, +0.746982248696},
	{-0.867666135683, -0.198076389613, +0.455983794521},
}

// EquatorialToGalactic converts (alpha, delta) in ICRS degrees to
// galactic (l, b) in degrees.
func EquatorialToGalactic(alpha, delta float64) (l, b float64) {
	return applyRotation(eqToGalMatrix, alpha, delta)
}

// GalacticToEquatorial is the inverse.
func GalacticToEquatorial(l, b float64) (alpha, delta float64) {
	return applyRotationInverse(eqToGalMatrix, l, b)
}

// ------------------ Equatorial ↔ Ecliptic ------------------

// obliquityIAU2006 returns the obliquity of the ecliptic (epsilon) at a
// given Julian Date in terrestrial time, in radians. From IAU 2006
// Resolution B1.1, eq. (21) in Capitaine et al. 2003.
func obliquityIAU2006(jd float64) float64 {
	// T = Julian centuries since J2000.0
	T := (jd - 2451545.0) / 36525.0
	// Polynomial in arcseconds of arc.
	eps := 84381.406 - 46.836769*T - 0.0001831*T*T + 0.00200340*T*T*T -
		5.76e-7*T*T*T*T - 4.34e-8*T*T*T*T*T
	return eps * math.Pi / (180 * 3600)
}

// EquatorialToEcliptic converts (alpha, delta) in ICRS degrees to
// ecliptic (lambda, beta) at the given Julian Date (typically derived
// from DATE-OBS or MJD-OBS).
func EquatorialToEcliptic(alpha, delta, jd float64) (lambda, beta float64) {
	eps := obliquityIAU2006(jd)
	sinEps, cosEps := math.Sincos(eps)
	a := degToRad(alpha)
	d := degToRad(delta)
	sinA, cosA := math.Sincos(a)
	sinD, cosD := math.Sincos(d)
	// Ecliptic β = asin(sin δ cos ε - cos δ sin α sin ε)
	// Ecliptic λ = atan2(sin α cos ε + tan δ sin ε, cos α)
	sinB := sinD*cosEps - cosD*sinA*sinEps
	beta = radToDeg(math.Asin(clampUnit(sinB)))
	y := sinA*cosEps + (sinD/cosD)*sinEps
	x := cosA
	lambda = normalizeDeg(radToDeg(math.Atan2(y, x)))
	return lambda, beta
}

// EclipticToEquatorial is the inverse of EquatorialToEcliptic.
func EclipticToEquatorial(lambda, beta, jd float64) (alpha, delta float64) {
	eps := obliquityIAU2006(jd)
	sinEps, cosEps := math.Sincos(eps)
	L := degToRad(lambda)
	B := degToRad(beta)
	sinL, cosL := math.Sincos(L)
	sinB, cosB := math.Sincos(B)
	// α = atan2(sin λ cos ε - tan β sin ε, cos λ)
	// δ = asin(sin β cos ε + cos β sin ε sin λ)
	y := sinL*cosEps - (sinB/cosB)*sinEps
	x := cosL
	alpha = normalizeDeg(radToDeg(math.Atan2(y, x)))
	sinD := sinB*cosEps + cosB*sinEps*sinL
	delta = radToDeg(math.Asin(clampUnit(sinD)))
	return alpha, delta
}

// ------------------ Equatorial ↔ Supergalactic ------------------

// supergalacticToGalMatrix rotates supergalactic (SGL, SGB) → galactic
// (l, b). de Vaucouleurs et al. 1976: supergalactic north pole is at
// galactic (l=47.37°, b=6.32°); supergalactic longitude zero is at
// galactic (l=137.37°, b=0°).
var supergalacticToGalMatrix = [3][3]float64{
	{-0.735742574804, +0.677261296414, 0.0},
	{-0.074553778365, -0.080991471307, +0.993922590400},
	{+0.673145302109, +0.731271165817, +0.110081262237},
}

// GalacticToSupergalactic converts galactic (l, b) to supergalactic
// (SGL, SGB).
func GalacticToSupergalactic(l, b float64) (sgl, sgb float64) {
	return applyRotationInverse(supergalacticToGalMatrix, l, b)
}

// SupergalacticToGalactic is the inverse.
func SupergalacticToGalactic(sgl, sgb float64) (l, b float64) {
	return applyRotation(supergalacticToGalMatrix, sgl, sgb)
}

// EquatorialToSupergalactic chains ICRS → galactic → supergalactic.
func EquatorialToSupergalactic(alpha, delta float64) (sgl, sgb float64) {
	l, b := EquatorialToGalactic(alpha, delta)
	return GalacticToSupergalactic(l, b)
}

// SupergalacticToEquatorial is the inverse.
func SupergalacticToEquatorial(sgl, sgb float64) (alpha, delta float64) {
	l, b := SupergalacticToGalactic(sgl, sgb)
	return GalacticToEquatorial(l, b)
}

// ------------------ ICRS / FK5 / FK4 precession ------------------

// PrecessFK5 precesses equatorial coordinates from equinox equinox1 to
// equinox2, both given as Besselian/Julian years (e.g. 2000.0 for J2000).
// Uses the IAU 1976 precession model via the standard p-matrix formulas.
//
// This is an approximation adequate for arcsecond-level precision over
// ~1 century; higher-precision work should use the IAU 2006 precession
// or delegate to an ephemeris library.
func PrecessFK5(alpha, delta, equinox1, equinox2 float64) (a2, d2 float64) {
	// T in Julian centuries from equinox1 to equinox2.
	T := (equinox2 - equinox1) / 100.0
	// IAU 1976 precession angles (zeta, z, theta) in arcseconds.
	T0 := (equinox1 - 2000.0) / 100.0
	zeta := (2306.2181+1.39656*T0-0.000139*T0*T0)*T +
		(0.30188-0.000344*T0)*T*T + 0.017998*T*T*T
	z := (2306.2181+1.39656*T0-0.000139*T0*T0)*T +
		(1.09468+0.000066*T0)*T*T + 0.018203*T*T*T
	theta := (2004.3109-0.85330*T0-0.000217*T0*T0)*T -
		(0.42665+0.000217*T0)*T*T - 0.041833*T*T*T
	// Convert arcseconds to radians.
	zeta *= math.Pi / (180 * 3600)
	z *= math.Pi / (180 * 3600)
	theta *= math.Pi / (180 * 3600)

	// Build the rotation matrix R = R_z(-z) R_y(theta) R_z(-zeta).
	m := precessionMatrix(zeta, z, theta)
	a2, d2 = applyRotation(m, alpha, delta)
	return a2, d2
}

// precessionMatrix builds the 3x3 precession rotation matrix from the
// three angles (all in radians).
func precessionMatrix(zeta, z, theta float64) [3][3]float64 {
	cZeta, sZeta := math.Cos(zeta), math.Sin(zeta)
	cZ, sZ := math.Cos(z), math.Sin(z)
	cTheta, sTheta := math.Cos(theta), math.Sin(theta)
	return [3][3]float64{
		{cZeta*cTheta*cZ - sZeta*sZ, -sZeta*cTheta*cZ - cZeta*sZ, -sTheta * cZ},
		{cZeta*cTheta*sZ + sZeta*cZ, -sZeta*cTheta*sZ + cZeta*cZ, -sTheta * sZ},
		{cZeta * sTheta, -sZeta * sTheta, cTheta},
	}
}

// ------------------ Helpers ------------------

// applyRotation multiplies the given 3x3 rotation matrix by the unit
// vector of (lon, lat) in degrees and returns the resulting (lon', lat').
func applyRotation(m [3][3]float64, lon, lat float64) (lon2, lat2 float64) {
	l := degToRad(lon)
	b := degToRad(lat)
	cosB := math.Cos(b)
	x := cosB * math.Cos(l)
	y := cosB * math.Sin(l)
	z := math.Sin(b)
	x2 := m[0][0]*x + m[0][1]*y + m[0][2]*z
	y2 := m[1][0]*x + m[1][1]*y + m[1][2]*z
	z2 := m[2][0]*x + m[2][1]*y + m[2][2]*z
	lon2 = normalizeDeg(radToDeg(math.Atan2(y2, x2)))
	lat2 = radToDeg(math.Asin(clampUnit(z2)))
	return lon2, lat2
}

// applyRotationInverse multiplies by the transpose of the given matrix.
// For an orthogonal (rotation) matrix this is the inverse rotation.
func applyRotationInverse(m [3][3]float64, lon, lat float64) (lon2, lat2 float64) {
	mT := [3][3]float64{
		{m[0][0], m[1][0], m[2][0]},
		{m[0][1], m[1][1], m[2][1]},
		{m[0][2], m[1][2], m[2][2]},
	}
	return applyRotation(mT, lon, lat)
}

// clampUnit clamps x to [-1, 1] — guards against 1 ULP overshoot before
// asin.
func clampUnit(x float64) float64 {
	if x > 1 {
		return 1
	}
	if x < -1 {
		return -1
	}
	return x
}
