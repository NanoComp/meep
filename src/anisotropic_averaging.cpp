#include <math.h>

#include "meep.h"
#include "threevec.h"

namespace meep {

static int count_quadrants(const geometric_volume &v) {
  return 1 << number_of_directions(v.dim);
}

static vec geo_center(const geometric_volume &v) {
  vec o = zero_vec(v.dim);
  LOOP_OVER_DIRECTIONS(v.dim,d)
    o.set_direction(d, 0.5*(v.in_direction_max(d) + v.in_direction_min(d)));
  return o;
}

static geometric_volume nth_quadrant(const geometric_volume &v, int n) {
  geometric_volume o = v;
  vec cent = geo_center(v);
  LOOP_OVER_DIRECTIONS(v.dim,d) {
    if (n & 1) o.set_direction_min(d, cent.in_direction(d));
    else o.set_direction_max(d, cent.in_direction(d));
    n = n >> 1;
  }
  return o;
}

#define USE_SPHERE_QUAD 0

#if USE_SPHERE_QUAD
#include "sphere-quad.h"

static int count_sphere_pts(const geometric_volume &v) {
     return num_sphere_quad[number_of_directions(v.dim)];
}

static inline double dmax(double a, double b) { return (a > b ? a : b); }

static double geo_diam(const geometric_volume &v) {
  double diam = 0.0;
  LOOP_OVER_DIRECTIONS(v.dim,d) {
    diam = dmax(diam, v.in_direction_max(d) - v.in_direction_min(d));
  }
  return diam;
}

static vec sphere_pt(const geometric_volume &v, int n, double &weight) {
     geometric_volume gv = v;
     vec cent = geo_center(v);
     vec pt;
     double R = geo_diam(gv) * 0.5;
     switch (v.dim) {
	 case D1:
	 {
	      weight = sphere_quad[0][n][3];
	      vec pt(sphere_quad[0][n][2]);
	      return cent + pt * R;
	 }
	 case D2:
	 {
	      weight = sphere_quad[1][n][3];
	      return cent + vec2d(sphere_quad[1][n][0], sphere_quad[1][n][1]) * R;
	 }
	 case D3:
	 {
	      weight = sphere_quad[2][n][3];
	      vec pt(sphere_quad[2][n][0], sphere_quad[2][n][1],
		     sphere_quad[2][n][2]);
	      return cent + pt * R;
	 }
	 case Dcyl:
	 {
	      weight = sphere_quad[1][n][3];
	      vec pt(sphere_quad[1][n][0], sphere_quad[1][n][1]);
	      return cent + pt * R;
	 }
     }
}
#endif /* USE_SPHERE_QUAD */

void prthv(threevec v) {
  master_printf("%10lg\t%10lg\t%10lg\n", v.val[0], v.val[1], v.val[2]);
}

void prtens(tensor t) {
  FOR3(i) prthv(t.row[i]);
}

void prgeo(geometric_volume v) {
  LOOP_OVER_DIRECTIONS(v.dim,d)
    master_printf("%g < %s < %g\n",
                  v.in_direction_min(d), direction_name(d), v.in_direction_max(d));
}

static double is_constant_inveps(material_function &eps, const geometric_volume &vol,
                                 double minvol) {
  double inveps = 1.0/eps.eps(geo_center(vol));
  if (vol.full_volume() <= minvol) return inveps;
  for (int i=0;i<count_quadrants(vol);i++) {
    geometric_volume here = nth_quadrant(vol, i);
    const double here_inveps = is_constant_inveps(eps, here, minvol);
    if (here_inveps != inveps) return -1.0;
  }
  return inveps;
}

static tensor doaverage_inveps(material_function &eps, const geometric_volume &vol,
                               double minvol) {
  if (vol.full_volume() <= minvol) return diagonal(1.0/eps.eps(geo_center(vol)));
  vec gradient = zero_vec(vol.dim);
  vec center = geo_center(vol);
  tensor mean = diagonal(0.0), meaninv = diagonal(0.0);
  for (int i=0;i<count_quadrants(vol);i++) {
    geometric_volume here = nth_quadrant(vol, i);
    tensor average_here = doaverage_inveps(eps, here, minvol);
    mean += (1.0/average_here);
    meaninv += average_here;
#if !USE_SPHERE_QUAD
    double invepshere = trace(average_here);
    gradient += (geo_center(here) - center)*invepshere;
#endif
  }
#if USE_SPHERE_QUAD
  for (int i = 0; i < count_sphere_pts(vol); ++i) {
    double weight;
    vec pt = sphere_pt(vol, i, weight);
    gradient += (pt - center) * (weight / eps.eps(pt));
  }
#endif
  mean = mean*(1.0/count_quadrants(vol));
  meaninv = meaninv*(1.0/count_quadrants(vol));
  threevec normdir;
  FOR3(i) normdir.val[i] = 0.0;
  LOOP_OVER_DIRECTIONS(vol.dim, d) normdir.val[d%3] = gradient.in_direction(d);
  if (abs(normdir)) normdir /= abs(normdir);
  else return meaninv;
  tensor project_norm(normdir);
  tensor project_parallel = diagonal(1.0) - project_norm;
  tensor invmean = 1.0/mean;
  return (project_parallel*invmean + invmean*project_parallel +
          project_norm*meaninv + meaninv*project_norm)*0.5;
}

double anisoaverage(component ec, direction d, material_function &eps,
                    const geometric_volume &vol, double minvol) {
  double const_inveps = is_constant_inveps(eps, vol, minvol);
  if (const_inveps >= 0.0) return (component_direction(ec) == d) ? const_inveps : 0.0;
  tensor avg = doaverage_inveps(eps, vol, minvol);
  int rownum = component_direction(ec) % 3;
  int colnum = d % 3;
  return avg.row[rownum].val[colnum];
}

} // namespace meep
