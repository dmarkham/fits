package stretch

import "math"

// GHTParams holds the parameters for the Generalized Hyperbolic
// Transformation (Strasser 2022). Port of Siril's ght command
// (ght.c).
type GHTParams struct {
	D  float64 // Stretch intensity (user-space; internally transformed to expm1(D))
	B  float64 // Localisation parameter (0 = exponential, -1 = log, >0 = power)
	LP float64 // Low protection point [0, SP]
	SP float64 // Symmetry point [0, 1]
	HP float64 // High protection point [SP, 1]
	BP float64 // Black point (only for linear type)
}

// ghtCoeffs holds the precomputed piecewise coefficients.
// These mirror Siril's ght_compute_params struct fields used for
// STRETCH_PAYNE_NORMAL.
type ghtCoeffs struct {
	b1             float64
	a2, b2, c2, d2 float64
	e2             float64
	a3, b3, c3, d3 float64
	e3             float64
	a4, b4         float64
}

// GHT applies the Generalized Hyperbolic Transformation to a pixel
// array. The stretch is applied independently per channel.
//
// The user-supplied D is transformed internally: actualD = expm1(D) = e^D - 1.
// This matches Siril's parameter convention (command.c:3316).
func GHT(pixels []float32, params GHTParams) {
	D := math.Expm1(params.D) // e^D - 1
	B := params.B
	if math.Abs(B) < 0.001 {
		B = 0
	}
	SP := params.SP
	LP := params.LP
	HP := params.HP
	if HP == 0 {
		HP = 1.0
	}

	if D == 0 {
		return // identity
	}

	// Precompute piecewise coefficients.
	c := ghtSetup(D, B, LP, SP, HP)

	for i, v := range pixels {
		pixels[i] = ghtEval(float64(v), D, B, LP, SP, HP, &c)
	}
}

// ghtSetup precomputes the piecewise polynomial coefficients for the
// GHT STRETCH_PAYNE_NORMAL mode. Direct port of Siril's GHTsetup
// (ght.c lines 34-115).
//
// Siril's approach: compute qlp, q0, qwp, q1 (the unnormalized
// function values at the boundaries), then q = 1/(q1-q0) as the
// global normalization factor. All coefficients are expressed in
// terms of q0 and q so the output spans [0, 1].
func ghtSetup(D, B, LP, SP, HP float64) ghtCoeffs {
	var c ghtCoeffs

	if B == -1.0 {
		// Logarithmic family (B == -1).
		qlp := -math.Log1p(D * (SP - LP))
		q0 := qlp - D*LP/(1.0+D*(SP-LP))
		qwp := math.Log1p(D * (HP - SP))
		q1 := qwp + D*(1.0-HP)/(1.0+D*(HP-SP))
		q := 1.0 / (q1 - q0)

		c.b1 = (1.0 + D*(SP-LP)) / (D * q)
		c.a2 = (-q0) * q
		c.b2 = -q
		c.c2 = 1.0 + D*SP
		c.d2 = -D
		c.a3 = (-q0) * q
		c.b3 = q
		c.c3 = 1.0 - D*SP
		c.d3 = D
		c.a4 = (qwp - q0 - D*HP/(1.0+D*(HP-SP))) * q
		c.b4 = q * D / (1.0 + D*(HP-SP))

	} else if B < 0.0 {
		// Generalized power, B < 0 (B != -1). Siril negates B.
		B = -B
		qlp := (1.0 - math.Pow(1.0+D*B*(SP-LP), (B-1.0)/B)) / (B - 1.0)
		q0 := qlp - D*LP*math.Pow(1.0+D*B*(SP-LP), -1.0/B)
		qwp := (math.Pow(1.0+D*B*(HP-SP), (B-1.0)/B) - 1.0) / (B - 1.0)
		q1 := qwp + D*(1.0-HP)*math.Pow(1.0+D*B*(HP-SP), -1.0/B)
		q := 1.0 / (q1 - q0)

		c.b1 = D * math.Pow(1.0+D*B*(SP-LP), -1.0/B) * q
		c.a2 = (1.0/(B-1.0) - q0) * q
		c.b2 = -q / (B - 1.0)
		c.c2 = 1.0 + D*B*SP
		c.d2 = -D * B
		c.e2 = (B - 1.0) / B
		c.a3 = (-1.0/(B-1.0) - q0) * q
		c.b3 = q / (B - 1.0)
		c.c3 = 1.0 - D*B*SP
		c.d3 = D * B
		c.e3 = (B - 1.0) / B
		c.a4 = (qwp - q0 - D*HP*math.Pow(1.0+D*B*(HP-SP), -1.0/B)) * q
		c.b4 = D * math.Pow(1.0+D*B*(HP-SP), -1.0/B) * q

	} else if B == 0.0 {
		// Exponential family (B == 0).
		qlp := math.Exp(-D * (SP - LP))
		q0 := qlp - D*LP*math.Exp(-D*(SP-LP))
		qwp := 2.0 - math.Exp(-D*(HP-SP))
		q1 := qwp + D*(1.0-HP)*math.Exp(-D*(HP-SP))
		q := 1.0 / (q1 - q0)

		c.b1 = D * math.Exp(-D*(SP-LP)) * q
		c.a2 = -q0 * q
		c.b2 = q
		c.c2 = -D * SP
		c.d2 = D
		c.a3 = (2.0 - q0) * q
		c.b3 = -q
		c.c3 = D * SP
		c.d3 = -D
		c.a4 = (qwp - q0 - D*HP*math.Exp(-D*(HP-SP))) * q
		c.b4 = D * math.Exp(-D*(HP-SP)) * q

	} else {
		// B > 0: power law.
		qlp := math.Pow(1.0+D*B*(SP-LP), -1.0/B)
		q0 := qlp - D*LP*math.Pow(1.0+D*B*(SP-LP), -(1.0+B)/B)
		qwp := 2.0 - math.Pow(1.0+D*B*(HP-SP), -1.0/B)
		q1 := qwp + D*(1.0-HP)*math.Pow(1.0+D*B*(HP-SP), -(1.0+B)/B)
		q := 1.0 / (q1 - q0)

		c.b1 = D * math.Pow(1.0+D*B*(SP-LP), -(1.0+B)/B) * q
		c.a2 = -q0 * q
		c.b2 = q
		c.c2 = 1.0 + D*B*SP
		c.d2 = -D * B
		c.e2 = -1.0 / B
		c.a3 = (2.0 - q0) * q
		c.b3 = -q
		c.c3 = 1.0 - D*B*SP
		c.d3 = D * B
		c.e3 = -1.0 / B
		c.a4 = (qwp - q0 - D*HP*math.Pow(1.0+D*B*(HP-SP), -(1.0+B)/B)) * q
		c.b4 = D * math.Pow(1.0+D*B*(HP-SP), -(1.0+B)/B) * q
	}

	return c
}

func ghtEval(in, D, B, LP, SP, HP float64, c *ghtCoeffs) float32 {
	in = math.Max(0, in)
	if D == 0 {
		return float32(in)
	}
	var res1, res2 float64
	if B == -1.0 {
		res1 = c.a2 + c.b2*math.Log(c.c2+c.d2*in)
		res2 = c.a3 + c.b3*math.Log(c.c3+c.d3*in)
	} else if B == 0.0 {
		res1 = c.a2 + c.b2*math.Exp(c.c2+c.d2*in)
		res2 = c.a3 + c.b3*math.Exp(c.c3+c.d3*in)
	} else {
		res1 = c.a2 + c.b2*math.Pow(c.c2+c.d2*in, c.e2)
		res2 = c.a3 + c.b3*math.Pow(c.c3+c.d3*in, c.e3)
	}
	var out float64
	switch {
	case in < LP:
		out = c.b1 * in
	case in < SP:
		out = res1
	case in < HP:
		out = res2
	default:
		out = c.a4 + c.b4*in
	}
	return float32(math.Max(0, math.Min(1, out)))
}
