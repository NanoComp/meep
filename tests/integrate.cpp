/* Copyright (C) 2004 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* Check of fields::integrate, by giving it random volumes in which to
   integrate purely linear functions of the coordinates--by
   construction, we should be able to integrate these exactly. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "meep.h"
using namespace meep;

const double size[3] = {3.0,3.0,2.6};

static double one(const vec &p) {
  (void) p;
  return 1.0;
}

typedef struct {
  double c, ax,ay,az, axy,ayz,axz, axyz;
  long double sum;
} linear_integrand_data;

/* integrand for integrating c + ax*x + ay*y + az*z. */
static void linear_integrand(fields_chunk *fc,
			     ivec is, ivec ie,
			     vec s0, vec s1, vec e0, vec e1,
			     double dV0, double dV1,
			     vec shift, complex<double> shift_phase, 
			     const symmetry &S, int sn,
			     void *data_)
{
  linear_integrand_data *data = (linear_integrand_data *) data_;
  
  /* we don't use the Bloch phase here, but it's important for
     integrating fields, obviously */
  (void) shift_phase;

  // get the integration weight along a given dimension
#define WEIGHT(i, n, dir) ((i > 1 && i < n - 2) ? 1.0 : (i == 0 ? s0.in_direction(direction(dir)) : (i == 1 ? s1.in_direction(direction(dir)) : i == n - 1 ? e0.in_direction(direction(dir)) : (i == n - 2 ? e1.in_direction(direction(dir)) : 1.0))))

  vec loc(fc->v.dim);
  LOOP_OVER_IVECS(fc->v, is, ie, idx) {
    // slightly evil use of loop vars defined internally by LOOP_OVER_IVECS...
    double w1 = WEIGHT(loop_i1, loop_n1, loop_d1);
    double dV = dV0 + dV1 * loop_i2;
    double w12 = w1 * WEIGHT(loop_i2, loop_n2, loop_d2) * dV;
    double w123 = w12 * WEIGHT(loop_i3, loop_n3, loop_d3);
    double inva = fc->v.inva;
    loc.set_direction(direction(loop_d1), loop_is1*0.5*inva + loop_i1*inva);
    loc.set_direction(direction(loop_d2), loop_is2*0.5*inva + loop_i2*inva);
    loc.set_direction(direction(loop_d3), loop_is3*0.5*inva + loop_i3*inva);

    vec locS(S.transform(loc, sn) + shift);
    if (0) master_printf("at (%g,%g,%g) with weight*dV = %g*%g\n",
			 locS.in_direction(X),
			 locS.in_direction(Y),
			 locS.in_direction(Z),
			 w123/dV, dV);
    data->sum +=
      w123 * (data->c
	      + data->ax * locS.in_direction(X)
	      + data->ay * locS.in_direction(Y)
	      + data->az * locS.in_direction(Z)
	      + data->axy * locS.in_direction(X) * locS.in_direction(Y)
	      + data->ayz * locS.in_direction(Z) * locS.in_direction(Y)
	      + data->axz * locS.in_direction(X) * locS.in_direction(Z)
	      + data->axyz * locS.in_direction(X) * locS.in_direction(Y)
	                   * locS.in_direction(Z)
	      );
  }

#undef WEIGHT
}

/* integrals of 1 and x, respectively, from a to b, or 1 and x if a==b: */
static double integral1(double a, double b) { return a==b ? 1 : b-a; }
static double integralx(double a, double b) { return a==b ? a : (b*b-a*a)*.5; }

static double correct_integral(const geometric_volume &gv,
			       const linear_integrand_data &data)
{
  double x1 = gv.in_direction_min(X);
  double x2 = gv.in_direction_max(X);
  double y1 = gv.in_direction_min(Y);
  double y2 = gv.in_direction_max(Y);
  double z1 = gv.in_direction_min(Z);
  double z2 = gv.in_direction_max(Z);

  return (data.c * integral1(x1,x2) * integral1(y1,y2) * integral1(z1,z2)
       + data.ax * integralx(x1,x2) * integral1(y1,y2) * integral1(z1,z2)
       + data.ay * integral1(x1,x2) * integralx(y1,y2) * integral1(z1,z2)
       + data.az * integral1(x1,x2) * integral1(y1,y2) * integralx(z1,z2)
      + data.axy * integralx(x1,x2) * integralx(y1,y2) * integral1(z1,z2)
      + data.ayz * integral1(x1,x2) * integralx(y1,y2) * integralx(z1,z2)
      + data.axz * integralx(x1,x2) * integral1(y1,y2) * integralx(z1,z2)
     + data.axyz * integralx(x1,x2) * integralx(y1,y2) * integralx(z1,z2)
	  );
}

// uniform pseudo-random number in [min,max]
static double urand(double min, double max)
{
  return (rand() * ((max - min) / RAND_MAX) + min);
}

static geometric_volume random_gv(void)
{
  geometric_volume gv(D3);

  gv.set_direction_min(X, urand(-100, 100));
  gv.set_direction_min(Y, urand(-100, 100));
  gv.set_direction_min(Z, urand(-100, 100));

  double s[3] = {0,0,0};
  switch (rand() % 4) { /* dimensionality */
  case 0:
    break;
  case 1: {
    int d = rand() % 3;
    s[d] = urand(0, size[d]);
    break;
  }
  case 2: {
    int d1 = rand() % 3;
    int d2 = (d1 + 1 + rand() % 2) % 3;
    s[d1] = urand(0, size[d1]);
    s[d2] = urand(0, size[d2]);
    break;
  }
  case 3:
    s[0] = urand(0, size[0]);
    s[1] = urand(0, size[1]);
    s[2] = urand(0, size[2]);
  }
  gv.set_direction_max(X, s[0] + gv.in_direction_min(X));
  gv.set_direction_max(Y, s[1] + gv.in_direction_min(Y));
  gv.set_direction_max(Z, s[2] + gv.in_direction_min(Z));
  return gv;
}

void check_integral(fields &f,
		    linear_integrand_data &d, const geometric_volume &gv)
{
  double x1 = gv.in_direction_min(X);
  double x2 = gv.in_direction_max(X);
  double y1 = gv.in_direction_min(Y);
  double y2 = gv.in_direction_max(Y);
  double z1 = gv.in_direction_min(Z);
  double z2 = gv.in_direction_max(Z);
  
  if (0)
    master_printf("Checking volume (%g,%g,%g) at (%g,%g,%g) with integral (%g, %g,%g,%g, %g,%g,%g, %g)...\n",
		  x2 - x1, y2 - y1, z2 - z1,
		  (x1+x2)/2, (y1+y2)/2, (z1+z2)/2,
		  d.c, d.ax,d.ay,d.az, d.axy,d.ayz,d.axz, d.axyz);
  else
    master_printf("Checking %d-dim. volume with %s integral...",
		  (x2 - x1 > 0) + (y2 - y1 > 0) + (z2 - z1 > 0), 
		  (d.c == 1.0 && !d.axy && !d.ax && !d.ay && !d.az
		   && !d.axy && !d.ayz && !d.axz) ? "unit" : "linear");
  d.sum = 0.0;
  f.integrate(linear_integrand, (void *) &d, gv);
  d.sum = sum_to_all(d.sum);
  if (fabs(d.sum - correct_integral(gv, d)) > 1e-9 * fabs(d.sum))
    abort("FAILED: %0.16g instead of %0.16g\n", 
	  (double) d.sum, correct_integral(gv, d));
  master_printf("...PASSED.\n");
}

void check_splitsym(const volume &v, int splitting, 
		    const symmetry &S, const char *Sname)
{
  const int num_random_trials = 100;
  structure s(v, one, splitting, S);
  fields f(&s);

  // periodic boundaries:
  f.use_bloch(X, 0);
  f.use_bloch(Y, 0);
  f.use_bloch(Z, 0);

  master_printf("\nCHECKS for splitting=%d, symmetry=%s\n...",
		splitting, Sname);
  for (int i = 0; i < num_random_trials; ++i) {
    geometric_volume gv(random_gv());
    linear_integrand_data d;

    // try integral of 1 first (easier to debug, I hope)
    d.c = 1.0;
    d.ax = d.ay = d.az = d.axy = d.ayz = d.axz = d.axyz = 0.0;
    check_integral(f, d, gv);

    d.c = urand(-1,1);
    d.ax = urand(-1,1); d.ay = urand(-1,1); d.az = urand(-1,1);
    d.axy = urand(-1,1); d.ayz = urand(-1,1); d.axz = urand(-1,1);
    d.axyz = urand(-1,1);
    check_integral(f, d, gv);
  }
}

int main(int argc, char **argv)
{
  const double a = 10.0;
  initialize mpi(argc, argv);
  const volume v = vol3d(size[0], size[1], size[2], a);

  srand(0); // use fixed random sequence

  for (int splitting = 0; splitting < 5; ++splitting) {
    check_splitsym(v, splitting, identity(), "identity");
    check_splitsym(v, splitting, mirror(X, v), "mirrorx");
    check_splitsym(v, splitting, mirror(Y, v), "mirrory");
    check_splitsym(v, splitting, mirror(X, v) + mirror(Y, v), "mirrorxy");
    check_splitsym(v, splitting, rotate4(Z, v), "rotate4");
  }
  return 0;
}
