// hpref: healpix_cxx truth harness for the polar-boundary golden.
//
// It links the real healpix_cxx (the library healpy wraps; the source
// tree is NOT copied here, see Makefile HEALPIX_CXX and ../PROVENANCE.md)
// and answers one query per stdin line (floats in C99 hex):
//
//   a nside nest theta phi        healpix_cxx ang2pix (calls libm cos/sin)
//   v nside nest x y z            healpix_cxx vec2pix (calls libm atan2)
//   A nside nest theta phi c s    ang2pix's body with cos(theta)=c and
//                                 sin(theta)=s supplied by the caller
//   V nside nest x y z phi        vec2pix's body with atan2(y,x)=phi
//                                 supplied by the caller
//   C / S nside nest theta phi c s  loc2pix forced to the cos (C) or the
//                                 sin (S) formula, whatever theta is
//   Q / P nside nest x y z phi    vec2pix's body forced to the z-only (Q)
//                                 or the sth (P) formula, whatever |z| is
//
// C/S/Q/P are not truth: gen_golden.py uses them only to pick inputs
// where the branch choice changes the pixel (write_branch_switch).
//
// The A/V bodies are copied verbatim from healpix_base.h ang2pix/vec2pix,
// with only the libm call replaced by the supplied value, and call the
// real (protected) T_Healpix_Base<int64>::loc2pix. gen_golden.py feeds
// them 256-bit mpmath values rounded once to double (correctly rounded
// cos/sin/atan2) and checks a/v against healpy itself.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#define protected public
#include "healpix_base.h"
#undef protected

typedef T_Healpix_Base<int64> HB;

static int64 ang_body(const HB &hb, double theta, double phi, double c, double s)
  {
  return ((theta<0.01) || (theta > 3.14159-0.01)) ?
    hb.loc2pix(c,phi,s,true) :
    hb.loc2pix(c,phi,0.,false);
  }

static int64 vec_body(const HB &hb, const vec3 &vec, double phi)
  {
  double xl = 1./vec.Length();
  double nz = vec.z*xl;
  if (std::abs(nz)>0.99)
    return hb.loc2pix (nz,phi,sqrt(vec.x*vec.x+vec.y*vec.y)*xl,true);
  else
    return hb.loc2pix (nz,phi,0,false);
  }

static int64 vec_forced(const HB &hb, const vec3 &vec, double phi, bool sth)
  {
  double xl = 1./vec.Length();
  double nz = vec.z*xl;
  if (sth)
    return hb.loc2pix (nz,phi,sqrt(vec.x*vec.x+vec.y*vec.y)*xl,true);
  return hb.loc2pix (nz,phi,0,false);
  }

int main()
  {
  char line[512];
  while (fgets(line, sizeof line, stdin))
    {
    char op; long long nside; int nest;
    char *p = line;
    op = *p++;
    nside = strtoll(p, &p, 10);
    nest = (int)strtol(p, &p, 10);
    double f[6]; int n = 0;
    while (n < 6)
      {
      char *q;
      double v = strtod(p, &q);
      if (q == p) break;
      f[n++] = v; p = q;
      }
    HB hb(nside, nest ? NEST : RING, SET_NSIDE);
    int64 r;
    switch (op)
      {
      case 'a': r = hb.ang2pix(pointing(f[0], f[1])); break;
      case 'v': r = hb.vec2pix(vec3(f[0], f[1], f[2])); break;
      case 'A': r = ang_body(hb, f[0], f[1], f[2], f[3]); break;
      case 'V': r = vec_body(hb, vec3(f[0], f[1], f[2]), f[3]); break;
      case 'C': r = hb.loc2pix(f[2], f[1], 0., false); break;
      case 'S': r = hb.loc2pix(f[2], f[1], f[3], true); break;
      case 'Q': r = vec_forced(hb, vec3(f[0], f[1], f[2]), f[3], false); break;
      case 'P': r = vec_forced(hb, vec3(f[0], f[1], f[2]), f[3], true); break;
      default: fprintf(stderr, "bad op %c\n", op); return 2;
      }
    printf("%lld\n", (long long)r);
    }
  return 0;
  }
