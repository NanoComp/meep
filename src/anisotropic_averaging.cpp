#include "meep.h"
#include "threevec.h"

namespace meep {

static int count_quadrants(const geometric_volume &v) {
  return 2 << number_of_directions(v.dim);
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
  return o;
}

static tensor doaverage_inveps(material_function &eps, const geometric_volume &vol,
                               double minvol) {
  if (vol.full_volume() <= minvol) return tensor(1.0/eps.eps(geo_center(vol)));
  tensor averages[8]; // averages hold local averages of *inveps*
  vec gradient = zero_vec(vol.dim);
  for (int i=0;i<count_quadrants(vol);i++) {
    geometric_volume here = nth_quadrant(vol, i);
    averages[i] = doaverage_inveps(eps, here, minvol);
    double invepshere = trace(averages[i]);
    gradient += geo_center(here)*invepshere;
  }
  threevec normdir;
  FOR3(i) normdir.val[i] = 0.0;
  LOOP_OVER_DIRECTIONS(vol.dim, d) normdir.val[d%3] = gradient.in_direction(d);
  normdir /= abs(normdir);
  tensor project_norm(normdir);
  tensor project_parallel = tensor(1.0) - project_norm;
  tensor mean(0.0), meaninv(0.0);
  for (int i=0;i<count_quadrants(vol);i++) {
    mean += project_norm*(1.0/averages[i]);
    meaninv += project_parallel*averages[i];
  }
  return 1.0/mean + meaninv;
}

double anisoaverage(component ec, direction d, material_function &eps,
                    const geometric_volume &vol, double minvol) {
  tensor avg = doaverage_inveps(eps, vol, minvol);
  int rownum = component_direction(ec) % 3;
  int colnum = d % 3;
  return avg.row[rownum].val[colnum];
}

} // namespace meep
