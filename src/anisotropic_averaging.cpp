#include "meep.h"
#include "threevec.h"

namespace meep {

static int count_quadrants(const geometric_volume &v) {
  return 1 << number_of_directions(v.dim);
}

static vec geo_center(geometric_volume v) {
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
}

static tensor doaverage_inveps(material_function &eps, const geometric_volume &vol,
                               double minvol) {
  if (vol.full_volume() <= minvol) return tensor(1.0/eps(geo_center(vol)));
  tensor averages[8]; // averages hold local averages of *inveps*
  vec gradient = zero_vec(vol.dim);
  tensor mean = diagonal(0.0), meaninv = diagonal(0.0);
  for (int i=0;i<count_quadrants(vol);i++) {
    geometric_volume here = nth_quadrant(vol, i);
    tensor average_here = doaverage_inveps(eps, here, minvol);
    mean += (1.0/average_here);
    meaninv += average_here;
    double invepshere = trace(average_here);
    gradient += (geo_center(here) - geo_center(vol))*invepshere;
  }
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
