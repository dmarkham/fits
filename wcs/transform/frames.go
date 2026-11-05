package transform

import "math"

// Sky coordinate-frame conversions.
//
// This file converts directions (no distance, no proper motion, no
// aberration) between these celestial frames:
//
//   - ICRS and FK5 J2000 ↔ Galactic (IAU 1958 system as realised in
//     astropy.coordinates)
//   - Galactic ↔ Supergalactic (de Vaucouleurs et al. 1976 pole and
//     origin, as astropy), and ICRS ↔ Supergalactic through Galactic
//   - ICRS ↔ mean ecliptic and equinox of date (IAU 2006 precession and
//     obliquity, no nutation; SOFA iauEcm06)
//   - FK5 equinox ↔ FK5 equinox (IAU 1976 precession, Lieske et al. 1977)
//   - FK4 B1950 (with E-terms of aberration) ↔ FK5 J2000 (Standish 1982
//     matrix plus the Murray 1989 FK4 rotation term, as astropy)
//
// All inputs and outputs are in DEGREES. Output longitudes are normalized
// to [0, 360); latitudes are in [-90, 90].
//
// Every function is proven point by point against astropy.coordinates 8.0.2
// (and ERFA's iauPrec76 for PrecessFK5) to 1 micro-arcsecond; see
// frames_golden_test.go and testdata/frames/PROVENANCE.md.
//
// The conversions use 3-vector rotations: convert (lon, lat) to a unit
// vector, multiply by a 3x3 rotation, convert back. The fixed matrices are
// built once from their defining angles. The helpers rot3X, rot3Y and
// rot3Z follow the SOFA (iauRx/iauRy/iauRz) and astropy rotation_matrix
// convention: they rotate the coordinate FRAME by the given angle, e.g.
// rot3Z(a) = [[cos a, sin a, 0], [-sin a, cos a, 0], [0, 0, 1]].

// ------------------ Equatorial ↔ Galactic ------------------

// Defining angles of the galactic frame in FK5 J2000 (astropy
// Galactic._ngp_J2000 / _lon0_J2000): the IAU 1958 FK4 B1950 definition
// (pole at 192.25, +27.4; galactic longitude of the celestial pole 123
// deg) carried into FK5 J2000.
const (
	galNGPRaJ2000  = 192.8594812065348
	galNGPDecJ2000 = 27.12825118085622
	galLon0J2000   = 122.9319185680026
)

// ICRS → FK5 J2000 frame bias (USNO Circular 179): offsets of the FK5
// J2000 pole and equinox from the ICRS, in milliarcseconds. Same constants
// and composition as astropy's _icrs_to_fk5_matrix.
const (
	biasEta0 = -19.9
	biasXi0  = 9.1
	biasDa0  = -22.9
)

var (
	// fk5ToGalMatrix rotates FK5 J2000 (alpha, delta) to galactic (l, b).
	fk5ToGalMatrix = mat3Mul(rot3Z(degToRad(180-galLon0J2000)),
		mat3Mul(rot3Y(degToRad(90-galNGPDecJ2000)), rot3Z(degToRad(galNGPRaJ2000))))

	// icrsToFK5Matrix is the frame bias matrix.
	icrsToFK5Matrix = mat3Mul(rot3X(-masToRad(biasEta0)),
		mat3Mul(rot3Y(masToRad(biasXi0)), rot3Z(masToRad(biasDa0))))

	// eqToGalMatrix rotates ICRS (alpha, delta) to galactic (l, b): frame
	// bias to FK5 J2000, then fk5ToGalMatrix (astropy's ICRS → FK5 →
	// Galactic path).
	eqToGalMatrix = mat3Mul(fk5ToGalMatrix, icrsToFK5Matrix)
)

// EquatorialToGalactic converts ICRS (alpha, delta) in degrees to galactic
// (l, b) in degrees, matching astropy ICRS → Galactic (the ICRS → FK5
// frame bias is applied, then the FK5 J2000 galactic rotation).
//
// Catalogues that compute (l, b) by applying the FK5 J2000 galactic
// matrix directly to ICRS positions (the Hipparcos and Gaia convention)
// differ from this by the frame bias, up to ~0.03 arcsec. Use
// FK5ToGalactic on the ICRS values to reproduce those.
func EquatorialToGalactic(alpha, delta float64) (l, b float64) {
	return applyRotation(eqToGalMatrix, alpha, delta)
}

// GalacticToEquatorial is the inverse of EquatorialToGalactic: galactic
// (l, b) to ICRS (alpha, delta), in degrees.
func GalacticToEquatorial(l, b float64) (alpha, delta float64) {
	return applyRotationInverse(eqToGalMatrix, l, b)
}

// FK5ToGalactic converts FK5 equinox J2000 (alpha, delta) in degrees to
// galactic (l, b) in degrees (astropy FK5(equinox=J2000) → Galactic).
func FK5ToGalactic(alpha, delta float64) (l, b float64) {
	return applyRotation(fk5ToGalMatrix, alpha, delta)
}

// GalacticToFK5 is the inverse of FK5ToGalactic.
func GalacticToFK5(l, b float64) (alpha, delta float64) {
	return applyRotationInverse(fk5ToGalMatrix, l, b)
}

// ------------------ Equatorial ↔ Ecliptic ------------------

// arcsecToRad converts arcseconds to radians.
func arcsecToRad(as float64) float64 { return as * math.Pi / (180 * 3600) }

// masToRad converts milliarcseconds to radians.
func masToRad(mas float64) float64 { return arcsecToRad(mas / 1000) }

// obliquityIAU2006 returns the mean obliquity of the ecliptic of date
// (epsilon_A) at a Julian Date in TT, in radians. IAU 2006 (Capitaine et
// al. 2003; SOFA iauObl06).
func obliquityIAU2006(jd float64) float64 {
	t := (jd - 2451545.0) / 36525.0
	return arcsecToRad(84381.406 + float64((-46.836769+float64((-0.0001831+float64((0.00200340+float64((-0.000000576+float64((-0.0000000434)*t))*t))*t))*t))*t))
}

// eclipticMatrix returns the rotation from ICRS to the mean ecliptic and
// equinox of date jd (TT): the IAU 2006 frame bias + precession matrix
// from the Fukushima-Williams angles (SOFA iauPfw06 + iauFw2m), followed
// by a rotation about the mean equinox by the mean obliquity of date. This
// is SOFA iauEcm06, which astropy's BarycentricMeanEcliptic uses.
func eclipticMatrix(jd float64) [3][3]float64 {
	t := (jd - 2451545.0) / 36525.0
	gamb := arcsecToRad(-0.052928 + float64((10.556378+float64((0.4932044+float64((-0.00031238+float64((-0.000002788+float64((0.0000000260)*t))*t))*t))*t))*t))
	phib := arcsecToRad(84381.412819 + float64((-46.811016+float64((0.0511268+float64((0.00053289+float64((-0.000000440+float64((-0.0000000176)*t))*t))*t))*t))*t))
	psib := arcsecToRad(-0.041775 + float64((5038.481484+float64((1.5584175+float64((-0.00018522+float64((-0.000026452+float64((-0.0000000148)*t))*t))*t))*t))*t))
	epsa := obliquityIAU2006(jd)
	// iauFw2m: R = Rx(-epsa) . Rz(-psib) . Rx(phib) . Rz(gamb).
	rbp := mat3Mul(rot3X(-epsa), mat3Mul(rot3Z(-psib), mat3Mul(rot3X(phib), rot3Z(gamb))))
	return mat3Mul(rot3X(epsa), rbp)
}

// EquatorialToEcliptic converts ICRS (alpha, delta) in degrees to mean
// ecliptic coordinates (lambda, beta) in degrees, referred to the mean
// ecliptic and equinox of the given Julian Date in TT (e.g. MJD-OBS +
// 2400000.5). Pass jd = 2451545.0 for the J2000 mean ecliptic.
//
// IAU 2006 precession and obliquity; nutation is not applied (the "mean"
// ecliptic). Matches astropy BarycentricMeanEcliptic(equinox=jd) for
// directions.
func EquatorialToEcliptic(alpha, delta, jd float64) (lambda, beta float64) {
	return applyRotation(eclipticMatrix(jd), alpha, delta)
}

// EclipticToEquatorial is the inverse of EquatorialToEcliptic: mean
// ecliptic of date (lambda, beta) to ICRS (alpha, delta), in degrees.
func EclipticToEquatorial(lambda, beta, jd float64) (alpha, delta float64) {
	return applyRotationInverse(eclipticMatrix(jd), lambda, beta)
}

// ------------------ Galactic ↔ Supergalactic ------------------

// Defining angles of the supergalactic frame (astropy
// Supergalactic._nsgp_gal): north supergalactic pole at galactic
// (l = 47.37, b = +6.32) deg; SGL = 0 lies at galactic (137.37, 0).
const (
	sgPoleL = 47.37
	sgPoleB = 6.32
)

// galToSGMatrix rotates galactic (l, b) to supergalactic (SGL, SGB):
// astropy rotation_matrix(90, z) @ rotation_matrix(90 - b_p, y) @
// rotation_matrix(l_p, z).
var galToSGMatrix = mat3Mul(rot3Z(degToRad(90)),
	mat3Mul(rot3Y(degToRad(90-sgPoleB)), rot3Z(degToRad(sgPoleL))))

// GalacticToSupergalactic converts galactic (l, b) to supergalactic
// (SGL, SGB), in degrees.
func GalacticToSupergalactic(l, b float64) (sgl, sgb float64) {
	return applyRotation(galToSGMatrix, l, b)
}

// SupergalacticToGalactic is the inverse.
func SupergalacticToGalactic(sgl, sgb float64) (l, b float64) {
	return applyRotationInverse(galToSGMatrix, sgl, sgb)
}

// EquatorialToSupergalactic chains ICRS → galactic → supergalactic.
func EquatorialToSupergalactic(alpha, delta float64) (sgl, sgb float64) {
	l, b := EquatorialToGalactic(alpha, delta)
	return GalacticToSupergalactic(l, b)
}

// SupergalacticToEquatorial is the inverse (supergalactic → ICRS).
func SupergalacticToEquatorial(sgl, sgb float64) (alpha, delta float64) {
	l, b := SupergalacticToGalactic(sgl, sgb)
	return GalacticToEquatorial(l, b)
}

// ------------------ FK5 precession ------------------

// PrecessFK5 precesses FK5 equatorial coordinates (degrees) from equinox
// equinox1 to equinox2, both JULIAN epochs (2000.0 = J2000.0; Besselian
// B1950.0 is Julian 1949.99979). Uses the IAU 1976 precession angles
// (zeta, z, theta) of Lieske et al. 1977 for arbitrary start and end
// epochs, as SOFA iauPrec76, and agrees with ERFA to micro-arcseconds
// (proven for spans up to 200 years).
//
// IAU 1976 is the model FK5 is defined with. It is not the IAU 2006 model
// that astropy's FK5 → FK5 transform uses; those differ by ~0.3 arcsec per
// century of span, a model difference rather than an error. Directions
// only: no proper motion is applied.
func PrecessFK5(alpha, delta, equinox1, equinox2 float64) (a2, d2 float64) {
	// T in Julian centuries from equinox1 to equinox2.
	T := (equinox2 - equinox1) / 100.0
	// IAU 1976 precession angles (zeta, z, theta) in arcseconds.
	T0 := (equinox1 - 2000.0) / 100.0
	zeta := float64((2306.2181+float64(1.39656*T0)-float64(0.000139*T0*T0))*T) +
		float64((0.30188-float64(0.000344*T0))*T*T) + float64(0.017998*T*T*T)
	z := float64((2306.2181+float64(1.39656*T0)-float64(0.000139*T0*T0))*T) +
		float64((1.09468+float64(0.000066*T0))*T*T) + float64(0.018203*T*T*T)
	theta := float64((2004.3109-float64(0.85330*T0)-float64(0.000217*T0*T0))*T) -
		float64((0.42665+float64(0.000217*T0))*T*T) - float64(0.041833*T*T*T)

	// R = Rz(-z) Ry(theta) Rz(-zeta)  (SOFA iauPmat76 composition).
	m := mat3Mul(rot3Z(-arcsecToRad(z)), mat3Mul(rot3Y(arcsecToRad(theta)), rot3Z(-arcsecToRad(zeta))))
	return applyRotation(m, alpha, delta)
}

// ------------------ FK4 ↔ FK5 ------------------

// fk4B1950ToJ2000 is the Standish (1982) rotation from FK4 B1950 (E-terms
// removed) to FK5 J2000 for epoch B1950 (astropy _B1950_TO_J2000_M).
var fk4B1950ToJ2000 = [3][3]float64{
	{0.9999256794956877, -0.0111814832204662, -0.0048590038153592},
	{0.0111814832391717, +0.9999374848933135, -0.0000271625947142},
	{0.0048590037723143, -0.0000271702937440, +0.9999881946023742},
}

// fk4Corr is the time-dependent correction to fk4B1950ToJ2000 per Julian
// century of observation epoch from 1950, because FK4 rotates relative to
// an inertial frame (Murray 1989 eq. 29; astropy _FK4_CORR, x 1e-6).
var fk4Corr = [3][3]float64{
	{-0.0026455262e-6, -1.1539918689e-6, +2.1111346190e-6},
	{+1.1540628161e-6, -0.0129042997e-6, +0.0236021478e-6},
	{-2.1112979048e-6, -0.0056024448e-6, +0.0102587734e-6},
}

// besselianToJD converts a Besselian epoch to a Julian Date (TT), SOFA
// iauEpb2jd.
func besselianToJD(epb float64) float64 {
	return 2400000.5 + 15019.81352 + float64((epb-1900.0)*365.242198781)
}

// fk4BMatrix returns the FK4 (no E-terms, B1950) → FK5 J2000 rotation for
// an observation epoch given as a Besselian year.
func fk4BMatrix(obstimeB float64) [3][3]float64 {
	jyear := 2000.0 + (besselianToJD(obstimeB)-2451545.0)/365.25
	t := (jyear - 1950.0) / 100.0
	var m [3][3]float64
	for i := range 3 {
		for j := range 3 {
			m[i][j] = fk4B1950ToJ2000[i][j] + float64(fk4Corr[i][j]*t)
		}
	}
	return m
}

// fk4ETerms returns the E-terms of aberration vector for the FK4 equinox
// at Julian Date jd (astropy fk4_e_terms; obliquity IAU 1980, SOFA
// iauObl80).
func fk4ETerms(jd float64) [3]float64 {
	const k = 0.0056932 // degrees
	T := (jd - besselianToJD(1950)) / 36525.0
	ek := float64(degToRad(k) * (float64((float64(-0.000000126*T)-0.00004193)*T) + 0.01673011))
	g := degToRad((float64((float64((float64(0.012*T)+1.65)*T)+6190.67)*T) + 1015489.951) / 3600.0)
	t80 := (jd - 2451545.0) / 36525.0
	o := arcsecToRad(84381.448 + float64((-46.8150+float64((-0.00059+float64(0.001813*t80))*t80))*t80))
	sg, cg := math.Sincos(g)
	so, co := math.Sincos(o)
	return [3]float64{float64(ek * sg), float64(-ek * cg * co), float64(-ek * cg * so)}
}

// FK4ToFK5 converts FK4 equinox B1950.0 (alpha, delta) in degrees,
// including the E-terms of aberration as FK4 catalogue positions do, to
// FK5 equinox J2000.0, for a source observed at Besselian epoch obstimeB
// (1950.0 for a plain B1950 catalogue position). Directions only: proper
// motion is not modelled; the obstime dependence is FK4's rotation
// relative to FK5, exactly as astropy FK4(equinox=B1950, obstime) →
// FK5(equinox=J2000).
//
// Only equinox B1950.0 is supported on the FK4 side.
func FK4ToFK5(alpha, delta, obstimeB float64) (a2, d2 float64) {
	v := lonLatToVec3(alpha, delta)
	e := fk4ETerms(besselianToJD(1950))
	// Remove the E-terms: v - A + (A.v) v (astropy fk4_to_fk4_no_e); the
	// length is not renormalised because vec3ToLonLat reads only the
	// direction.
	dot := float64(e[0]*v[0]) + float64(e[1]*v[1]) + float64(e[2]*v[2])
	for i := range 3 {
		v[i] = v[i] - e[i] + float64(dot*v[i])
	}
	return vec3ToLonLat(mat3Vec(fk4BMatrix(obstimeB), v))
}

// FK5ToFK4 is the inverse of FK4ToFK5: FK5 J2000.0 (alpha, delta) to FK4
// B1950.0 with E-terms, for a source observed at Besselian epoch obstimeB.
func FK5ToFK4(alpha, delta, obstimeB float64) (a2, d2 float64) {
	v := mat3Vec(mat3T(fk4BMatrix(obstimeB)), lonLatToVec3(alpha, delta))
	e := fk4ETerms(besselianToJD(1950))
	// Add the E-terms. astropy (fk4.py fk4_no_e_to_fk4) iterates
	// r = (v + A) / (1 + A.r) ten times; the scalar 1 + A.r only rescales
	// r, and vec3ToLonLat reads only the direction, so every pass gives
	// the same direction as v + A. That direction is used directly (an
	// iteration count or the denominator could not change any output, so
	// they are not kept as untestable code).
	r := [3]float64{v[0] + e[0], v[1] + e[1], v[2] + e[2]}
	return vec3ToLonLat(r)
}

// ------------------ Helpers ------------------

// rot3X, rot3Y, rot3Z are frame rotations by angle a (radians) about the
// x, y and z axes (SOFA iauRx / iauRy / iauRz, astropy rotation_matrix).
func rot3X(a float64) [3][3]float64 {
	s, c := math.Sincos(a)
	return [3][3]float64{{1, 0, 0}, {0, c, s}, {0, -s, c}}
}

func rot3Y(a float64) [3][3]float64 {
	s, c := math.Sincos(a)
	return [3][3]float64{{c, 0, -s}, {0, 1, 0}, {s, 0, c}}
}

func rot3Z(a float64) [3][3]float64 {
	s, c := math.Sincos(a)
	return [3][3]float64{{c, s, 0}, {-s, c, 0}, {0, 0, 1}}
}

// mat3Mul returns a·b.
func mat3Mul(a, b [3][3]float64) [3][3]float64 {
	var m [3][3]float64
	for i := range 3 {
		for j := range 3 {
			m[i][j] = float64(a[i][0]*b[0][j]) + float64(a[i][1]*b[1][j]) + float64(a[i][2]*b[2][j])
		}
	}
	return m
}

// mat3T returns the transpose (the inverse, for a rotation).
func mat3T(m [3][3]float64) [3][3]float64 {
	return [3][3]float64{
		{m[0][0], m[1][0], m[2][0]},
		{m[0][1], m[1][1], m[2][1]},
		{m[0][2], m[1][2], m[2][2]},
	}
}

func mat3Vec(m [3][3]float64, v [3]float64) [3]float64 {
	return [3]float64{
		float64(m[0][0]*v[0]) + float64(m[0][1]*v[1]) + float64(m[0][2]*v[2]),
		float64(m[1][0]*v[0]) + float64(m[1][1]*v[1]) + float64(m[1][2]*v[2]),
		float64(m[2][0]*v[0]) + float64(m[2][1]*v[1]) + float64(m[2][2]*v[2]),
	}
}

// lonLatToVec3 returns the unit vector of (lon, lat) in degrees.
func lonLatToVec3(lon, lat float64) [3]float64 {
	sl, cl := math.Sincos(degToRad(lon))
	sb, cb := math.Sincos(degToRad(lat))
	return [3]float64{float64(cb * cl), float64(cb * sl), sb}
}

// vec3ToLonLat returns (lon, lat) in degrees of a (not necessarily unit)
// vector, lon in [0, 360). Latitude uses atan2(z, hypot(x, y)), which
// keeps full precision near the poles where asin(z) loses half the digits.
func vec3ToLonLat(v [3]float64) (lon, lat float64) {
	lon = normalizeDeg(radToDeg(math.Atan2(v[1], v[0])))
	if lon >= 360 { // a longitude a hair below 0 can round up to 360
		lon = 0
	}
	lat = radToDeg(math.Atan2(v[2], math.Hypot(v[0], v[1])))
	return lon, lat
}

// applyRotation multiplies the given 3x3 rotation matrix by the unit
// vector of (lon, lat) in degrees and returns the resulting (lon', lat').
func applyRotation(m [3][3]float64, lon, lat float64) (lon2, lat2 float64) {
	return vec3ToLonLat(mat3Vec(m, lonLatToVec3(lon, lat)))
}

// applyRotationInverse multiplies by the transpose of the given matrix.
func applyRotationInverse(m [3][3]float64, lon, lat float64) (lon2, lat2 float64) {
	return applyRotation(mat3T(m), lon, lat)
}
