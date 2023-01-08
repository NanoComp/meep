/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
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

#include <meep.hpp>
using namespace meep;
using std::complex;

double sz[3] = {3.0, 3.0, 2.6};

static double one(const vec &p) {
  (void)p;
  return 1.0;
}

typedef struct {
  direction dx, dy, dz;
  double c, ax, ay, az, axy, ayz, axz, axyz;
  long double sum;
} linear_integrand_data;

/* integrand for integrating c + ax*x + ay*y + az*z. */
static complex<double> linear_integrand(const complex<realnum> *fields, const vec &loc,
                                        void *data_) {
  linear_integrand_data *data = (linear_integrand_data *)data_;

  (void)fields; // unused

  // clean_vec is only necessary because we reference X/Y/Z for any gv.dim
  vec locS(clean_vec(loc));

  return (data->c + data->ax * locS.in_direction(data->dx) +
          data->ay * locS.in_direction(data->dy) + data->az * locS.in_direction(data->dz)

          + data->axy * locS.in_direction(data->dx) * locS.in_direction(data->dy)

          + data->ayz * locS.in_direction(data->dz) * locS.in_direction(data->dy)

          + data->axz * locS.in_direction(data->dx) * locS.in_direction(data->dz)

          + data->axyz * locS.in_direction(data->dx) * locS.in_direction(data->dy) *
                locS.in_direction(data->dz));
}

/* integrals of 1 and x, respectively, from a to b, or 1 and x if a==b: */
static double integral1(double a, double b, direction d) {
  if (d == R)
    return a == b ? 2 * pi * a : pi * (b * b - a * a);
  else
    return a == b ? 1 : b - a;
}

static double integralx(double a, double b, direction d) {
  if (d == R)
    return a == b ? 2 * pi * a * a : 2 * pi * (b * b * b - a * a * a) * .3333333333333333333333333;
  else
    return a == b ? a : (b * b - a * a) * .5;
}

static double correct_integral(const volume &v, const linear_integrand_data &data) {
  direction x = data.dx, y = data.dy, z = data.dz;
  double x1 = v.in_direction_min(x);
  double x2 = v.in_direction_max(x);
  double y1 = v.in_direction_min(y);
  double y2 = v.in_direction_max(y);
  double z1 = v.in_direction_min(z);
  double z2 = v.in_direction_max(z);

  return (data.c * integral1(x1, x2, x) * integral1(y1, y2, y) * integral1(z1, z2, z) +
          data.ax * integralx(x1, x2, x) * integral1(y1, y2, y) * integral1(z1, z2, z) +
          data.ay * integral1(x1, x2, x) * integralx(y1, y2, y) * integral1(z1, z2, z) +
          data.az * integral1(x1, x2, x) * integral1(y1, y2, y) * integralx(z1, z2, z) +
          data.axy * integralx(x1, x2, x) * integralx(y1, y2, y) * integral1(z1, z2, z) +
          data.ayz * integral1(x1, x2, x) * integralx(y1, y2, y) * integralx(z1, z2, z) +
          data.axz * integralx(x1, x2, x) * integral1(y1, y2, y) * integralx(z1, z2, z) +
          data.axyz * integralx(x1, x2, x) * integralx(y1, y2, y) * integralx(z1, z2, z));
}

// uniform pseudo-random number in [min,max]
static double urand(double min, double max) { return (rand() * ((max - min) / RAND_MAX) + min); }

static volume random_gv(ndim dim) {
  volume v(dim);

  double s[3] = {0, 0, 0};
  int idim = dim == Dcyl ? 1 : int(dim);
  switch (rand() % (idim + 2)) { /* dimensionality */
    case 0: break;
    case 1: {
      int d = rand() % (idim + 1);
      s[d] = urand(0, sz[d]);
      break;
    }
    case 2: {
      int d1 = rand() % (idim + 1);
      int d2 = (d1 + 1 + rand() % 2) % 3;
      s[d1] = urand(0, sz[d1]);
      s[d2] = urand(0, sz[d2]);
      break;
    }
    case 3:
      s[0] = urand(0, sz[0]);
      s[1] = urand(0, sz[1]);
      s[2] = urand(0, sz[2]);
  }

  switch (dim) {
    case D1:
      v.set_direction_min(X, 0);
      v.set_direction_max(X, 0);
      v.set_direction_min(Y, 0);
      v.set_direction_max(Y, 0);
      v.set_direction_min(Z, urand(-100, 100));
      v.set_direction_max(Z, s[0] + v.in_direction_min(Z));
      break;
    case D2:
      v.set_direction_min(X, urand(-100, 100));
      v.set_direction_min(Y, urand(-100, 100));
      v.set_direction_max(X, s[0] + v.in_direction_min(X));
      v.set_direction_max(Y, s[1] + v.in_direction_min(Y));
      v.set_direction_min(Z, 0);
      v.set_direction_max(Z, 0);
      break;
    case Dcyl:
      v.set_direction_min(X, 0);
      v.set_direction_max(X, 0);
      v.set_direction_min(Y, 0);
      v.set_direction_max(Y, 0);
      v.set_direction_min(R, 0.1 + urand(0, sz[0] - s[0]));
      v.set_direction_min(Z, urand(-100, 100));
      v.set_direction_max(R, s[0] + v.in_direction_min(R));
      v.set_direction_max(Z, s[1] + v.in_direction_min(Z));
      v.set_direction_min(P, 0);
      v.set_direction_max(P, 0);
      break;
    case D3:
      v.set_direction_min(X, urand(-100, 100));
      v.set_direction_min(Y, urand(-100, 100));
      v.set_direction_max(X, s[0] + v.in_direction_min(X));
      v.set_direction_max(Y, s[1] + v.in_direction_min(Y));
      v.set_direction_min(Z, urand(-100, 100));
      v.set_direction_max(Z, s[2] + v.in_direction_min(Z));
      break;
    default: meep::abort("unsupported dimensionality in integrate.cpp");
  }

  return v;
}

void check_integral(fields &f, linear_integrand_data &d, const volume &v, component cgrid) {
  double x1 = v.in_direction_min(d.dx);
  double x2 = v.in_direction_max(d.dx);
  double y1 = v.in_direction_min(d.dy);
  double y2 = v.in_direction_max(d.dy);
  double z1 = v.in_direction_min(d.dz);
  double z2 = v.in_direction_max(d.dz);

  master_printf("Check %d-dim. %s integral in %s cell with %s integrand...",
                (x2 - x1 > 0) + (y2 - y1 > 0) + (z2 - z1 > 0), component_name(cgrid),
                v.dim == D3 ? "3d" : (v.dim == D2 ? "2d" : (v.dim == Dcyl ? "cylindrical" : "1d")),
                (d.c == 1.0 && !d.axy && !d.ax && !d.ay && !d.az && !d.axy && !d.ayz && !d.axz)
                    ? "unit"
                    : "linear");
  if (0)
    master_printf("\n... grid_volume (%g,%g,%g) at (%g,%g,%g) with integral (%g, %g,%g,%g, "
                  "%g,%g,%g, %g)...\n",
                  x2 - x1, y2 - y1, z2 - z1, (x1 + x2) / 2, (y1 + y2) / 2, (z1 + z2) / 2, d.c, d.ax,
                  d.ay, d.az, d.axy, d.ayz, d.axz, d.axyz);

  double sum = real(f.integrate(0, 0, linear_integrand, (void *)&d, v));
  if (fabs(sum - correct_integral(v, d)) > 1e-9 * fabs(sum))
    meep::abort("FAILED: %0.16g instead of %0.16g\n", sum, correct_integral(v, d));
  master_printf("...PASSED.\n");
}

void check_splitsym(const grid_volume &gv, int splitting, const symmetry &S, const char *Sname) {
  const int num_random_trials = 100;
  structure s(gv, one, no_pml(), S, splitting);
  fields f(&s);

  // periodic boundaries:
  f.use_bloch(zero_vec(gv.dim));

  linear_integrand_data d;
  if (gv.dim == Dcyl) {
    d.dx = R;
    d.dy = P;
    d.dz = Z;
  }
  else {
    d.dx = X;
    d.dy = Y;
    d.dz = Z;
  }

  master_printf("\nCHECKS for splitting=%d, symmetry=%s\n...", splitting, Sname);
  for (int i = 0; i < num_random_trials; ++i) {
    volume v(random_gv(gv.dim));
    component cgrid;

    do {
      cgrid = component(rand() % (Dielectric + 1));
    } while (coordinate_mismatch(gv.dim, component_direction(cgrid)));

    // try integral of 1 first (easier to debug, I hope)
    d.c = 1.0;
    d.ax = d.ay = d.az = d.axy = d.ayz = d.axz = d.axyz = 0.0;
    check_integral(f, d, v, cgrid);

    d.c = urand(-1, 1);
    d.ax = urand(-1, 1);
    d.ay = urand(-1, 1);
    d.az = urand(-1, 1);
    d.axy = urand(-1, 1);
    d.ayz = urand(-1, 1);
    d.axz = urand(-1, 1);
    d.axyz = urand(-1, 1);
    if (gv.dim == Dcyl) // cyl. doesn't integrate linear functions of r exactly
      d.ax = d.axy = d.axz = d.axyz = 0;
    check_integral(f, d, v, cgrid);
  }
}

// check LOOP_OVER_VOL and LOOP_OVER_VOL_OWNED macros
void check_loop_vol(const grid_volume &gv, component c) {
  size_t count = 0, count_owned = 0;
  ptrdiff_t min_i = ptrdiff_t(gv.ntot()), max_i = 0;

  master_printf("Checking %s loops for %s grid_volume...\n", component_name(c),
                dimension_name(gv.dim));

  ivec vmin(gv.little_corner() + gv.iyee_shift(c));
  ivec vmax(gv.big_corner() + gv.iyee_shift(c));

  LOOP_OVER_VOL(gv, c, i) {
    IVEC_LOOP_ILOC(gv, ihere);
    IVEC_LOOP_LOC(gv, here);
    ivec ihere0(gv.iloc(c, i));
    vec here0(gv[ihere0]);
    if (ihere0 != ihere) meep::abort("FAILED: wrong LOOP_OVER_VOL iloc at i=%td\n", i);
    if (abs(here0 - here) > 1e-13)
      meep::abort("FAILED: wrong LOOP_OVER_VOL loc (err = %g) at i=%td\n", abs(here0 - here), i);
    ++count;
    if (i < min_i) min_i = i;
    if (i > max_i) max_i = i;
    if (gv.owns(ihere)) ++count_owned;
    if (ihere < vmin || ihere > vmax) meep::abort("FAILED: LOOP_OVER_VOL outside V at i=%td\n", i);
  }
  if (count != gv.ntot())
    meep::abort("FAILED: LOOP_OVER_VOL has %zd iterations instead of ntot=%zd\n", count, gv.ntot());
  if (count_owned != gv.nowned(c))
    meep::abort("FAILED: LOOP_OVER_VOL has %zd owned points instead of nowned=%zd\n", count_owned,
                gv.nowned(c));
  if (min_i != 0) meep::abort("FAILED: LOOP_OVER_VOL has minimum index %td instead of 0\n", min_i);
  if (size_t(max_i) != gv.ntot() - 1)
    meep::abort("FAILED: LOOP_OVER_VOL has max index %td instead of ntot-1\n", max_i);

  count = 0;
  LOOP_OVER_VOL_OWNED(gv, c, i) {
    IVEC_LOOP_ILOC(gv, ihere);
    IVEC_LOOP_LOC(gv, here);
    ivec ihere0(gv.iloc(c, i));
    vec here0(gv[ihere0]);
    if (ihere0 != ihere) meep::abort("FAILED: wrong LOOP_OVER_VOL_OWNED iloc at i=%td\n", i);
    if (abs(here0 - here) > 1e-13)
      meep::abort("FAILED: wrong LOOP_OVER_VOL_OWNED loc (err = %g) at i=%td\n", abs(here0 - here),
                  i);
    if (!gv.owns(ihere))
      meep::abort("FAILED: LOOP_OVER_VOL_OWNED includes non-owned at i=%td\n", i);
    ++count;
  }
  if (count != count_owned)
    meep::abort("FAILED: LOOP_OVER_VOL_OWNED has %zd iterations instead of %zd\n", count,
                count_owned);

  count = 0;
  LOOP_OVER_VOL_NOTOWNED(gv, c, i) {
    IVEC_LOOP_ILOC(gv, ihere);
    IVEC_LOOP_LOC(gv, here);
    ivec ihere0(gv.iloc(c, i));
    vec here0(gv[ihere0]);
    if (ihere0 != ihere) meep::abort("FAILED: wrong LOOP_OVER_VOL_NOTOWNED iloc at i=%td\n", i);
    if (abs(here0 - here) > 1e-13)
      meep::abort("FAILED: wrong LOOP_OVER_VOL_NOTOWNED loc (err = %g) at i=%td\n",
                  abs(here0 - here), i);
    if (gv.owns(ihere)) meep::abort("FAILED: LOOP_OVER_VOL_NOTOWNED includes owned at i=%td\n", i);
    if (ihere < vmin || ihere > vmax)
      meep::abort("FAILED: LOOP_OVER_VOL_NOTOWNED outside V at i=%td\n", i);
    ++count;
  }
  if (count != gv.ntot() - count_owned)
    meep::abort("FAILED: LOOP_OVER_VOL_NOTOWNED has %zd iterations instead of %zd\n", count,
                gv.ntot() - count_owned);

  master_printf("...PASSED.\n");
}

int main(int argc, char **argv) {
  const double a = 10.0;
  initialize mpi(argc, argv);
  verbosity = 0;
  const grid_volume v3d = vol3d(sz[0], sz[1], sz[2], a);
  const grid_volume v3d0 = vol3d(sz[0], sz[1], 0, a);
  const grid_volume v3d00 = vol3d(sz[0], 0, 0, a);
  const grid_volume v2d = vol2d(sz[0], sz[1], a);
  const grid_volume v1d = vol1d(sz[0], a);
  const grid_volume vcyl = volcyl(sz[0], sz[1], a);

  srand(0); // use fixed random sequence

  check_splitsym(v3d, 0, identity(), "identity");
  check_splitsym(v3d, 0, mirror(X, v3d), "mirrorx");
  return 0;

  for (int splitting = 0; splitting < 5; ++splitting) {
    check_splitsym(v3d, splitting, identity(), "identity");
    check_splitsym(v3d, splitting, mirror(X, v3d), "mirrorx");
    check_splitsym(v3d, splitting, mirror(Y, v3d), "mirrory");
    check_splitsym(v3d, splitting, mirror(X, v3d) + mirror(Y, v3d), "mirrorxy");
    check_splitsym(v3d, splitting, rotate4(Z, v3d), "rotate4");
  }

  for (int splitting = 0; splitting < 5; ++splitting) {
    check_splitsym(v2d, splitting, identity(), "identity");
    check_splitsym(v2d, splitting, mirror(X, v2d), "mirrorx");
    check_splitsym(v2d, splitting, mirror(Y, v2d), "mirrory");
    check_splitsym(v2d, splitting, mirror(X, v2d) + mirror(Y, v2d), "mirrorxy");
    check_splitsym(v2d, splitting, rotate4(Z, v2d), "rotate4");
  }

  const grid_volume vcyl_pad = volcyl(sz[0] + 0.2, sz[1], a);
  for (int splitting = 0; splitting < 5; ++splitting) {
    check_splitsym(vcyl_pad, splitting, identity(), "identity");
    check_splitsym(vcyl_pad, splitting, mirror(Z, vcyl), "mirrorz");
  }

  for (int splitting = 0; splitting < 5; ++splitting) {
    check_splitsym(v1d, splitting, identity(), "identity");
    check_splitsym(v1d, splitting, mirror(Z, v1d), "mirrorz");
  }

  return 0;
}
