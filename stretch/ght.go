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
// GHT. Port of Siril's GHTsetup (ght.c).
func ghtSetup(D, B, LP, SP, HP float64) ghtCoeffs {
	var c ghtCoeffs

	// The GHT is piecewise on [0,LP], [LP,SP], [SP,HP], [HP,1].
	// Each region uses a different functional form depending on B.
	// We compute the coefficients to ensure C1 continuity.

	if B == 0 {
		// Exponential family.
		c.c2 = -D * SP
		c.d2 = D
		c.c3 = D * SP
		c.d3 = -D
		// At x=SP: both sides meet. Region 2 value at SP: exp(-D*SP + D*SP) = exp(0) = 1.
		// Region 3 value at SP: exp(D*SP - D*SP) = 1.
		// b2 and b3 scale so that the function spans [0, 1] properly.
		// See Siril ght.c for the exact coefficient derivation.
		exp_c2_d2LP := math.Exp(c.c2 + c.d2*LP)
		exp_c2_d2SP := math.Exp(c.c2 + c.d2*SP)
		exp_c3_d3SP := math.Exp(c.c3 + c.d3*SP)
		exp_c3_d3HP := math.Exp(c.c3 + c.d3*HP)

		if exp_c2_d2SP == exp_c2_d2LP {
			c.b2 = 1
		} else {
			c.b2 = 1.0 / (exp_c2_d2SP - exp_c2_d2LP)
		}
		c.a2 = -c.b2 * exp_c2_d2LP

		if exp_c3_d3SP == exp_c3_d3HP {
			c.b3 = 1
		} else {
			c.b3 = -1.0 / (exp_c3_d3HP - exp_c3_d3SP)
		}
		c.a3 = 1.0 - c.b3*exp_c3_d3HP

	} else if B == -1 {
		// Logarithmic family.
		c.c2 = 1 + D*SP
		c.d2 = -D
		c.c3 = 1 - D*SP
		c.d3 = D
		log_c2_d2LP := math.Log(c.c2 + c.d2*LP)
		log_c2_d2SP := math.Log(c.c2 + c.d2*SP)
		log_c3_d3SP := math.Log(c.c3 + c.d3*SP)
		log_c3_d3HP := math.Log(c.c3 + c.d3*HP)

		if log_c2_d2SP == log_c2_d2LP {
			c.b2 = 1
		} else {
			c.b2 = 1.0 / (log_c2_d2SP - log_c2_d2LP)
		}
		c.a2 = -c.b2 * log_c2_d2LP

		if log_c3_d3SP == log_c3_d3HP {
			c.b3 = 1
		} else {
			c.b3 = -1.0 / (log_c3_d3HP - log_c3_d3SP)
		}
		c.a3 = 1.0 - c.b3*log_c3_d3HP

	} else if B < 0 {
		// Generalized power (B < 0, B != -1).
		c.e2 = (B - 1) / B
		c.c2 = 1 + D*SP
		c.d2 = -D
		c.c3 = 1 - D*SP
		c.d3 = D
		c.e3 = c.e2

		pow_LP := math.Pow(c.c2+c.d2*LP, c.e2)
		pow_SP2 := math.Pow(c.c2+c.d2*SP, c.e2)
		pow_SP3 := math.Pow(c.c3+c.d3*SP, c.e3)
		pow_HP := math.Pow(c.c3+c.d3*HP, c.e3)

		if pow_SP2 == pow_LP {
			c.b2 = 1
		} else {
			c.b2 = 1.0 / (pow_SP2 - pow_LP)
		}
		c.a2 = -c.b2 * pow_LP

		if pow_SP3 == pow_HP {
			c.b3 = 1
		} else {
			c.b3 = -1.0 / (pow_HP - pow_SP3)
		}
		c.a3 = 1.0 - c.b3*pow_HP

	} else {
		// B > 0: power law with negative exponent.
		c.e2 = -1.0 / B
		c.c2 = 1 + D*SP
		c.d2 = -D
		c.c3 = 1 - D*SP
		c.d3 = D
		c.e3 = c.e2

		pow_LP := math.Pow(c.c2+c.d2*LP, c.e2)
		pow_SP2 := math.Pow(c.c2+c.d2*SP, c.e2)
		pow_SP3 := math.Pow(c.c3+c.d3*SP, c.e3)
		pow_HP := math.Pow(c.c3+c.d3*HP, c.e3)

		if pow_SP2 == pow_LP {
			c.b2 = 1
		} else {
			c.b2 = 1.0 / (pow_SP2 - pow_LP)
		}
		c.a2 = -c.b2 * pow_LP

		if pow_SP3 == pow_HP {
			c.b3 = 1
		} else {
			c.b3 = -1.0 / (pow_HP - pow_SP3)
		}
		c.a3 = 1.0 - c.b3*pow_HP
	}

	// Linear region [0, LP]: slope matches region 2 at LP.
	if LP > 0 {
		valLP := ghtRegion2(LP, B, &c)
		c.b1 = valLP / LP
	}

	// Linear region [HP, 1]: slope matches region 3 at HP.
	valHP := ghtRegion3(HP, B, &c)
	if HP < 1 {
		c.a4 = 1.0
		c.b4 = (1.0 - valHP) / (1.0 - HP)
		c.a4 = valHP - c.b4*HP
	} else {
		c.a4 = 0
		c.b4 = 1
	}

	return c
}

func ghtRegion2(x, B float64, c *ghtCoeffs) float64 {
	if B == 0 {
		return c.a2 + c.b2*math.Exp(c.c2+c.d2*x)
	} else if B == -1 {
		return c.a2 + c.b2*math.Log(c.c2+c.d2*x)
	}
	return c.a2 + c.b2*math.Pow(c.c2+c.d2*x, c.e2)
}

func ghtRegion3(x, B float64, c *ghtCoeffs) float64 {
	if B == 0 {
		return c.a3 + c.b3*math.Exp(c.c3+c.d3*x)
	} else if B == -1 {
		return c.a3 + c.b3*math.Log(c.c3+c.d3*x)
	}
	return c.a3 + c.b3*math.Pow(c.c3+c.d3*x, c.e3)
}

func ghtEval(in, D, B, LP, SP, HP float64, c *ghtCoeffs) float32 {
	in = math.Max(0, in)
	if D == 0 {
		return float32(in)
	}
	var out float64
	switch {
	case in < LP:
		out = c.b1 * in
	case in < SP:
		out = ghtRegion2(in, B, c)
	case in < HP:
		out = ghtRegion3(in, B, c)
	default:
		out = c.a4 + c.b4*in
	}
	return float32(math.Max(0, math.Min(1, out)))
}
