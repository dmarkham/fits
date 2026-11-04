/* libm_probe: how often glibc's float logf/powf/expf/asinhf differ from the
 * correctly rounded float32 (approximated by rounding the double function,
 * which differs from correct rounding with negligible probability).
 * Evidence for the GHT test tolerance (stretch/ght_test.go). Prints counts
 * for 2e6 deterministic pseudo-random arguments in the ranges GHT uses.
 *
 *   make libm-probe */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
	srand(1);
	long np = 0, ne = 0, nl = 0, na = 0, N = 2000000;
	for (long i = 0; i < N; i++) {
		float x = (float)rand() / RAND_MAX * 2.0f + 1e-3f;
		float y = (float)rand() / RAND_MAX * 4.0f - 2.0f;
		float e = (float)rand() / RAND_MAX * 20.0f - 10.0f;
		float s = x * 100.0f;
		if (powf(x, y) != (float)pow((double)x, (double)y)) np++;
		if (expf(e) != (float)exp((double)e)) ne++;
		if (logf(x) != (float)log((double)x)) nl++;
		if (asinhf(s) != (float)asinh((double)s)) na++;
	}
	printf("of %ld arguments: powf %ld, expf %ld, logf %ld, asinhf %ld differ from the correctly rounded float\n",
		N, np, ne, nl, na);
	return 0;
}
