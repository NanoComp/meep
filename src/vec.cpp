/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>

#include "meep_internals.hpp"
#include "meepgeom.hpp"

using namespace std;

namespace meep {

static bool isinteger(double value) { return value == floor(value); }

ivec grid_volume::round_vec(const vec &p) const {
  ivec result(dim);
  LOOP_OVER_DIRECTIONS(dim, d) { result.set_direction(d, my_round(p.in_direction(d) * 2 * a)); }
  return result;
}

void grid_volume::set_origin(const ivec &o) {
  io = o;
  origin = operator[](io); // adjust origin to match io
}

void grid_volume::set_origin(direction d, int o) {
  io.set_direction(d, o);
  origin = operator[](io); // adjust origin to match io
}

void grid_volume::set_origin(const vec &o) { set_origin(round_vec(o)); }

const char *dimension_name(ndim dim) {
  switch (dim) {
    case D1: return "1D";
    case D2: return "2D";
    case D3: return "3D";
    case Dcyl: return "Cylindrical";
  }
  return "Error in dimension_name";
}

const char *direction_name(direction d) {
  switch (d) {
    case X: return "x";
    case Y: return "y";
    case Z: return "z";
    case R: return "r";
    case P: return "phi";
    case NO_DIRECTION: return "no_direction";
  }
  return "Error in direction_name";
}

const char *component_name(component c) {
  if (is_derived(int(c))) return component_name(derived_component(c));
  switch (c) {
    case Er: return "er";
    case Ep: return "ep";
    case Ez: return "ez";
    case Hr: return "hr";
    case Hp: return "hp";
    case Hz: return "hz";
    case Ex: return "ex";
    case Ey: return "ey";
    case Hx: return "hx";
    case Hy: return "hy";
    case Dx: return "dx";
    case Dy: return "dy";
    case Dz: return "dz";
    case Dr: return "dr";
    case Dp: return "dp";
    case Bx: return "bx";
    case By: return "by";
    case Bz: return "bz";
    case Br: return "br";
    case Bp: return "bp";
    case Dielectric: return "eps";
    case Permeability: return "mu";
    case NO_COMPONENT: return "none";
  }
  return "Error in component_name";
}

const char *component_name(derived_component c) {
  if (!is_derived(int(c))) return component_name(component(c));
  switch (c) {
    case Sr: return "sr";
    case Sp: return "sp";
    case Sz: return "sz";
    case Sx: return "sx";
    case Sy: return "sy";
    case EnergyDensity: return "energy";
    case D_EnergyDensity: return "denergy";
    case H_EnergyDensity: return "henergy";
  }
  return "Error in component_name";
}

const char *component_name(int c) {
  return (is_derived(c) ? component_name(derived_component(c)) : component_name(component(c)));
}

component first_field_component(field_type ft) {
  switch (ft) {
    case E_stuff: return Ex;
    case H_stuff: return Hx;
    case D_stuff: return Dx;
    case B_stuff: return Bx;
    default: meep::abort("bug - only E/H/D/B stuff have components"); return NO_COMPONENT;
  }
}

vec min(const vec &vec1, const vec &vec2) {
  vec m(vec1.dim);
  LOOP_OVER_DIRECTIONS(vec1.dim, d) {
    m.set_direction(d, std::min(vec1.in_direction(d), vec2.in_direction(d)));
  }
  return m;
}

vec max(const vec &vec1, const vec &vec2) {
  vec m(vec1.dim);
  LOOP_OVER_DIRECTIONS(vec1.dim, d) {
    m.set_direction(d, std::max(vec1.in_direction(d), vec2.in_direction(d)));
  }
  return m;
}

ivec min(const ivec &ivec1, const ivec &ivec2) {
  ivec m(ivec1.dim);
  LOOP_OVER_DIRECTIONS(ivec1.dim, d) {
    m.set_direction(d, std::min(ivec1.in_direction(d), ivec2.in_direction(d)));
  }
  return m;
}

ivec max(const ivec &ivec1, const ivec &ivec2) {
  ivec m(ivec1.dim);
  LOOP_OVER_DIRECTIONS(ivec1.dim, d) {
    m.set_direction(d, std::max(ivec1.in_direction(d), ivec2.in_direction(d)));
  }
  return m;
}

volume::volume(const vec &vec1, const vec &vec2) {
  min_corner = min(vec1, vec2);
  max_corner = max(vec1, vec2);
  dim = vec1.dim;
}

volume::volume(const vec &pt) {
  dim = pt.dim;
  min_corner = pt;
  max_corner = pt;
}

double volume::computational_volume() const {
  double vol = 1.0;
  LOOP_OVER_DIRECTIONS(dim, d) { vol *= in_direction(d); }
  return vol;
}

double volume::integral_volume() const {
  double vol = 1.0;
  LOOP_OVER_DIRECTIONS(dim, d) {
    if (in_direction(d) != 0.0) vol *= in_direction(d);
  }
  if (dim == Dcyl) vol *= pi * (in_direction_max(R) + in_direction_min(R));
  return vol;
}

double volume::full_volume() const {
  double vol = computational_volume();
  if (dim == Dcyl) vol *= pi * (in_direction_max(R) + in_direction_min(R));
  return vol;
}

double volume::diameter() const {
  double diam = 0.0;
  LOOP_OVER_DIRECTIONS(dim, d) { diam = std::max(diam, in_direction(d)); }
  return diam;
}

volume volume::intersect_with(const volume &a) const {
  if (a.dim != dim) meep::abort("Can't intersect volumes of dissimilar dimensions.\n");
  volume result(dim);
  LOOP_OVER_DIRECTIONS(dim, d) {
    double minval = std::max(in_direction_min(d), a.in_direction_min(d));
    double maxval = std::min(in_direction_max(d), a.in_direction_max(d));
    if (minval > maxval) return volume(zero_vec(dim), zero_vec(dim));
    result.set_direction_min(d, minval);
    result.set_direction_max(d, maxval);
  }
  return result;
}

bool volume::intersects(const volume &a) const {
  if (a.dim != dim) meep::abort("Can't intersect volumes of dissimilar dimensions.\n");
  LOOP_OVER_DIRECTIONS(dim, d) {
    double minval = std::max(in_direction_min(d), a.in_direction_min(d));
    double maxval = std::min(in_direction_max(d), a.in_direction_max(d));
    if (minval > maxval) return false;
  }
  return true;
}

// Return normal direction to grid_volume, if the grid_volume is dim-1 dimensional;
// otherwise, return NO_DIRECTION.
direction volume::normal_direction() const {
  direction d = NO_DIRECTION;
  switch (dim) {
    case D1: d = Z; break;
    case D2:
      if (in_direction(X) == 0 && in_direction(Y) > 0)
        d = X;
      else if (in_direction(X) > 0 && in_direction(Y) == 0)
        d = Y;
      break;
    case Dcyl:
      if (in_direction(R) == 0 && in_direction(Z) > 0)
        d = R;
      else if (in_direction(R) > 0 && in_direction(Z) == 0)
        d = Z;
      break;
    case D3: {
      bool zx = in_direction(X) == 0;
      bool zy = in_direction(Y) == 0;
      bool zz = in_direction(Z) == 0;
      if (zx && !zy && !zz)
        d = X;
      else if (!zx && zy && !zz)
        d = Y;
      else if (!zx && !zy && zz)
        d = Z;
      break;
    }
  }
  return d;
}

/* Used for n=0,1,2 nested loops in macros.  We should arrange
   the ordering so that this gives most efficient traversal of
   a field array, where n=2 is the innermost loop. */
static direction yucky_dir(ndim dim, int n) {
  if (dim == Dcyl) switch (n) {
      case 0: return P;
      case 1: return R;
      case 2: return Z;
    }
  else if (dim == D2)
    return (direction)((n + 2) % 3); /* n = 0,1,2 gives Z, X, Y */
  return (direction)n;
}

int ivec::yucky_val(int n) const {
  if (has_direction(dim, yucky_dir(dim, n))) return in_direction(yucky_dir(dim, n));
  return 0;
}

int grid_volume::yucky_num(int n) const {
  if (has_direction(dim, yucky_dir(dim, n))) return num_direction(yucky_dir(dim, n));
  return 1;
}

direction grid_volume::yucky_direction(int n) const { return yucky_dir(dim, n); }

volume grid_volume::surroundings() const {
  return volume(operator[](little_corner()), operator[](big_corner()));
}

volume grid_volume::interior() const {
  return volume(operator[](little_corner()), operator[](big_corner() - one_ivec(dim) * 2));
}

void grid_volume::update_ntot() {
  the_ntot = 1;
  LOOP_OVER_DIRECTIONS(dim, d) { the_ntot *= (size_t)(num[d % 3] + 1); }
}

void grid_volume::set_num_direction(direction d, int value) {
  num[d % 3] = value;
  num_changed();
}

grid_volume::grid_volume(ndim td, double ta, int na, int nb, int nc) {
  dim = td;
  a = ta;
  inva = 1.0 / ta;
  num[0] = na;
  num[1] = nb;
  num[2] = nc;
  num_changed();
  set_origin(zero_vec(dim));
}

component grid_volume::eps_component() const {
  switch (dim) {
    case D1: return Hy;
    case D2: return Hz;
    case D3: return Dielectric;
    case Dcyl: return Hp;
  }
  meep::abort("Unsupported dimensionality eps.\n");
  return Ex;
}

vec grid_volume::yee_shift(component c) const { return operator[](iyee_shift(c)); }

/* Return array offsets to average with a given array location of c in
   order to get c on the "centered" grid.  Then, to get the
   centered grid point i, you should average c over the four
   locations: i, i+offset1, i+offset2, i+offset1+offset2.
   (offset2, and possibly offset1, may be zero if only 2 or 1
   locations need to be averaged). */
void grid_volume::yee2cent_offsets(component c, ptrdiff_t &offset1, ptrdiff_t &offset2) const {
  offset1 = offset2 = 0;
  LOOP_OVER_DIRECTIONS(dim, d) {
    if (!iyee_shift(c).in_direction(d)) {
      if (offset2) meep::abort("weird yee shift for component %s", component_name(c));
      if (offset1)
        offset2 = stride(d);
      else
        offset1 = stride(d);
    }
  }
}

/* Same as yee2cent_offsets, but averages centered grid to get c */
void grid_volume::cent2yee_offsets(component c, ptrdiff_t &offset1, ptrdiff_t &offset2) const {
  yee2cent_offsets(c, offset1, offset2);
  offset1 = -offset1;
  offset2 = -offset2;
}

bool volume::contains(const vec &p) const {
  LOOP_OVER_DIRECTIONS(dim, d) {
    if (p.in_direction(d) > in_direction_max(d) || p.in_direction(d) < in_direction_min(d))
      return false;
  }
  return true;
}

bool volume::contains(const volume &a) const {
  return contains(a.get_min_corner()) && contains(a.get_max_corner());
}

bool grid_volume::contains(const ivec &p) const {
  // containts returns true if the grid_volume has information about this grid
  // point.
  const ivec o = p - io;
  LOOP_OVER_DIRECTIONS(dim, d) {
    if (o.in_direction(d) < 0 || o.in_direction(d) >= (num_direction(d) + 1) * 2) return false;
  }
  return true;
}

bool grid_volume::contains(const vec &p) const {
  // containts returns true if the grid_volume has any information in it
  // relevant to the point p.  Basically has is like owns (see below)
  // except it is more lenient, in that more than one lattice may contain a
  // given point.
  const vec o = p - origin;
  LOOP_OVER_DIRECTIONS(dim, d) {
    if (o.in_direction(d) < -inva || o.in_direction(d) > num_direction(d) * inva + inva)
      return false;
  }
  return true;
}

/* Compute the corners (cs,ce) of the ib-th boundary for component c,
   returning true if ib is a valid index (ib = 0..#boundaries-1).  The
   boundaries are all the points that are in but not owned by the
   grid_volume, and are a set of *disjoint* regions.  The main purpose of
   this function is currently to support the LOOP_OVER_NOT_OWNED
   macro.  (In the future, it may be used for other
   boundary-element-type computations, too.) */
bool grid_volume::get_boundary_icorners(component c, int ib, ivec *cs, ivec *ce) const {
  ivec cl(little_corner() + iyee_shift(c));
  ivec cb(big_corner() + iyee_shift(c));
  ivec clo(little_owned_corner(c));
  ivec cbo(big_corner() - iyee_shift(c));
  *cs = cl;
  *ce = cb;
  bool ib_found = false;
  int jb = 0;
  LOOP_OVER_DIRECTIONS(dim, d) {
    if (cl.in_direction(d) < clo.in_direction(d)) {
      if (jb == ib) {
        ce->set_direction(d, cs->in_direction(d));
        ib_found = true;
        break;
      }
      cs->set_direction(d, clo.in_direction(d));
      jb++;
    }
    if (cb.in_direction(d) > cbo.in_direction(d)) {
      if (jb == ib) {
        cs->set_direction(d, ce->in_direction(d));
        ib_found = true;
        break;
      }
      ce->set_direction(d, cbo.in_direction(d));
      jb++;
    }
  }
  if (!ib_found) { // yucky interaction here with LOOP_OVER_VOL_NOTOWNED
    *cs = one_ivec(dim);
    *ce = -one_ivec(dim);
  }
  return ib_found;
}

// first "owned" point for c in grid_volume (see also grid_volume::owns)
ivec grid_volume::little_owned_corner(component c) const {
  ivec iloc(little_owned_corner0(c));
  if (dim == Dcyl && origin.r() == 0.0 && iloc.r() == 2) iloc.set_direction(R, 0);
  return iloc;
}

size_t grid_volume::nowned(component c) const {
  size_t n = 1;
  ivec pt = big_corner() - little_owned_corner(c);
  LOOP_OVER_DIRECTIONS(dim, d) { n *= pt.in_direction(d) / 2 + 1; }
  return n;
}

bool grid_volume::owns(const ivec &p) const {
  // owns returns true if the point "owned" by this grid_volume, meaning that it
  // is the grid_volume that would timestep the point.
  const ivec o = p - io;
  if (dim == Dcyl) {
    if (origin.r() == 0.0 && o.z() > 0 && o.z() <= nz() * 2 && o.r() == 0) return true;
    return o.r() > 0 && o.z() > 0 && o.r() <= nr() * 2 && o.z() <= nz() * 2;
  }
  else if (dim == D3) {
    return o.x() > 0 && o.x() <= nx() * 2 && o.y() > 0 && o.y() <= ny() * 2 && o.z() > 0 &&
           o.z() <= nz() * 2;
  }
  else if (dim == D2) { return o.x() > 0 && o.x() <= nx() * 2 && o.y() > 0 && o.y() <= ny() * 2; }
  else if (dim == D1) { return o.z() > 0 && o.z() <= nz() * 2; }
  else {
    meep::abort("Unsupported dimension in owns.\n");
    return false;
  }
}

int grid_volume::has_boundary(boundary_side b, direction d) const {
  switch (dim) {
    case Dcyl: return d == Z || (d == R && (b == High || get_origin().r() > 0));
    case D1: return d == Z;
    case D2: return d == X || d == Y;
    case D3: return d == X || d == Y || d == Z;
  }
  return 0; // This should never be reached.
}

ptrdiff_t grid_volume::index(component c, const ivec &p) const {
  const ivec offset = p - io - iyee_shift(c);
  ptrdiff_t idx = 0;
  LOOP_OVER_DIRECTIONS(dim, d) { idx += offset.in_direction(d) / 2 * stride(d); }
  return idx;
}

void grid_volume::set_strides() {
  FOR_DIRECTIONS(d) { the_stride[d] = 0; /* Yuck yuck yuck. */ }
  LOOP_OVER_DIRECTIONS(dim, d) {
    switch (d) {
      case Z: the_stride[d] = 1; break;
      case R: the_stride[d] = nz() + 1; break;
      case X: the_stride[d] = ptrdiff_t(nz() + 1) * (ny() + 1); break;
      case Y: the_stride[d] = nz() + 1; break;
      case P: break;            // There is no phi stride...
      case NO_DIRECTION: break; // no stride here, either
    }
  }
}

static inline void stupidsort(ptrdiff_t *ind, double *w, int l) {
  while (l) {
    if (fabs(w[0]) < 2e-15) {
      w[0] = w[l - 1];
      ind[0] = ind[l - 1];
      w[l - 1] = 0.0;
      ind[l - 1] = 0;
    }
    else {
      w += 1;
      ind += 1;
    }
    l -= 1;
  }
}

static inline void stupidsort(ivec *locs, double *w, int l) {
  while (l) {
    if (fabs(w[0]) < 2e-15) {
      w[0] = w[l - 1];
      locs[0] = locs[l - 1];
      w[l - 1] = 0.0;
      locs[l - 1] = 0;
    }
    else {
      w += 1;
      locs += 1;
    }
    l -= 1;
  }
}

void grid_volume::interpolate(component c, const vec &p, ptrdiff_t indices[8],
                              double weights[8]) const {
  ivec locs[8];
  interpolate(c, p, locs, weights);
  for (int i = 0; i < 8 && weights[i]; i++)
    if (!owns(locs[i])) weights[i] = 0.0;
  stupidsort(locs, weights, 8);
  for (int i = 0; i < 8 && weights[i]; i++)
    indices[i] = index(c, locs[i]);
  if (!contains(p) && weights[0]) {
    printf("Error at point %g %g\n", p.r(), p.z());
    printf("Interpolated to point %d %d\n", locs[0].r(), locs[0].z());
    printf("Or in other words... %g %g\n", operator[](locs[0]).r(), operator[](locs[0]).z());
    printf("I %s own the interpolated point.\n", owns(locs[0]) ? "actually" : "don't");
    print();
    meep::abort("Error made in interpolation of %s--fix this bug!!!\n", component_name(c));
  }
  // Throw out out of range indices:
  for (int i = 0; i < 8 && weights[i]; i++)
    if (indices[0] < 0 || size_t(indices[0]) >= ntot()) weights[i] = 0.0;
  // Stupid very crude code to compactify arrays:
  stupidsort(indices, weights, 8);
  if (!contains(p) && weights[0]) {
    printf("Error at point %g %g\n", p.r(), p.z());
    printf("Interpolated to point %d %d\n", locs[0].r(), locs[0].z());
    print();
    meep::abort("Error made in interpolation of %s--fix this bug!!!\n", component_name(c));
  }
}

void grid_volume::interpolate(component c, const vec &pc, ivec locs[8], double weights[8]) const {
  const double SMALL = 1e-13;
  const vec p = (pc - yee_shift(c)) * a;
  ivec middle(dim);
  LOOP_OVER_DIRECTIONS(dim, d) { middle.set_direction(d, ((int)floor(p.in_direction(d))) * 2 + 1); }
  middle += iyee_shift(c);
  const vec midv = operator[](middle);
  const vec dv = (pc - midv) * (2 * a);
  int already_have = 1;
  for (int i = 0; i < 8; i++) {
    locs[i] = round_vec(midv);
    weights[i] = 1.0;
  }
  LOOP_OVER_DIRECTIONS(dim, d) {
    for (int i = 0; i < already_have; i++) {
      locs[already_have + i] = locs[i];
      weights[already_have + i] = weights[i];
      locs[i].set_direction(d, middle.in_direction(d) - 1);
      weights[i] *= 0.5 * (1.0 - dv.in_direction(d));
      locs[already_have + i].set_direction(d, middle.in_direction(d) + 1);
      weights[already_have + i] *= 0.5 * (1.0 + dv.in_direction(d));
    }
    already_have *= 2;
  }
  for (int i = already_have; i < 8; i++)
    weights[i] = 0.0;
  double total_weight = 0.0;
  for (int i = 0; i < already_have; i++)
    total_weight += weights[i];
  for (int i = 0; i < already_have; i++)
    weights[i] += (1.0 - total_weight) * (1.0 / already_have);
  for (int i = 0; i < already_have; i++) {
    if (weights[i] < 0.0) {
      if (-weights[i] >= SMALL * 1e5)
        meep::abort("large negative interpolation weight[%d] = %e\n", i, weights[i]);
      weights[i] = 0.0;
    }
    else if (weights[i] < SMALL)
      weights[i] = 0.0;
  }
  stupidsort(locs, weights, already_have);
  // The rest of this code is a crude hack to get the weights right when we
  // are exactly between a few grid points.  i.e. to eliminate roundoff
  // error.
  bool all_same = true;
  for (int i = 0; i < 8 && weights[i]; i++)
    if (weights[i] != weights[0]) all_same = false;
  if (all_same) {
    int num_weights = 0;
    for (int i = 0; i < 8 && weights[i]; i++)
      num_weights++;
    for (int i = 0; i < 8 && weights[i]; i++)
      weights[i] = 1.0 / num_weights;
  }
}

volume empty_volume(ndim dim) {
  volume out(dim);
  LOOP_OVER_DIRECTIONS(dim, d) {
    out.set_direction_max(d, 0.0);
    out.set_direction_min(d, 0.0);
  }
  return out;
}

volume grid_volume::dV(const ivec &here, double diameter) const {
  const double hinva = 0.5 * inva * diameter;
  const grid_volume &gv = *this;
  const vec h = gv[here];
  volume out(dim);
  LOOP_OVER_DIRECTIONS(dim, d) {
    out.set_direction_max(d, h.in_direction(d) + hinva);
    out.set_direction_min(d, h.in_direction(d) - hinva);
  }
  if (dim == Dcyl && here.r() == 0) { out.set_direction_min(R, 0.0); }
  return out;
}

volume grid_volume::dV(component c, ptrdiff_t ind) const {
  if (!owns(iloc(c, ind))) return empty_volume(dim);
  return dV(iloc(c, ind));
}

double grid_volume::xmax() const {
  const double qinva = 0.25 * inva;
  return origin.x() + nx() * inva + qinva;
}

double grid_volume::xmin() const {
  const double qinva = 0.25 * inva;
  return origin.x() + qinva;
}

double grid_volume::ymax() const {
  const double qinva = 0.25 * inva;
  return origin.y() + ny() * inva + qinva;
}

double grid_volume::ymin() const {
  const double qinva = 0.25 * inva;
  return origin.y() + qinva;
}

double grid_volume::zmax() const {
  const double qinva = 0.25 * inva;
  return origin.z() + nz() * inva + qinva;
}

double grid_volume::zmin() const {
  const double qinva = 0.25 * inva;
  return origin.z() + qinva;
}

double grid_volume::rmax() const {
  const double qinva = 0.25 * inva;
  if (dim == Dcyl) return origin.r() + nr() * inva + qinva;
  meep::abort("No rmax in these dimensions.\n");
  return 0.0; // This is never reached.
}

double grid_volume::rmin() const {
  const double qinva = 0.25 * inva;
  if (dim == Dcyl) {
    if (origin.r() == 0.0) { return 0.0; }
    else { return origin.r() + qinva; }
  }
  meep::abort("No rmin in these dimensions.\n");
  return 0.0; // This is never reached.
}

double vec::project_to_boundary(direction d, double boundary_loc) {
  return fabs(boundary_loc - in_direction(d));
}

double grid_volume::boundary_location(boundary_side b, direction d) const {
  // Returns the location of metallic walls...
  if (b == High) switch (d) {
      case X: return loc(Ez, ntot() - 1).x();
      case Y: return loc(Ez, ntot() - 1).y();
      case R: return loc(Ep, ntot() - 1).r();
      case Z:
        if (dim == Dcyl)
          return loc(Ep, ntot() - 1).z();
        else
          return loc(Ex, ntot() - 1).z();
      case P: meep::abort("P has no boundary!\n");
      case NO_DIRECTION: meep::abort("NO_DIRECTION has no boundary!\n");
    }
  else
    switch (d) {
      case X: return loc(Ez, 0).x();
      case Y: return loc(Ez, 0).y();
      case R: return loc(Ep, 0).r();
      case Z:
        if (dim == Dcyl)
          return loc(Ep, 0).z();
        else
          return loc(Ex, 0).z();
      case P: meep::abort("P has no boundary!\n");
      case NO_DIRECTION: meep::abort("NO_DIRECTION has no boundary!\n");
    }
  return 0.0;
}

ivec grid_volume::big_corner() const {
  switch (dim) {
    case D1: return io + ivec(nz()) * 2;
    case D2: return io + ivec(nx(), ny()) * 2;
    case D3: return io + ivec(nx(), ny(), nz()) * 2;
    case Dcyl: return io + iveccyl(nr(), nz()) * 2;
  }
  return ivec(0); // This is never reached.
}

vec grid_volume::corner(boundary_side b) const {
  if (b == Low) return origin; // Low corner
  vec tmp = origin;
  LOOP_OVER_DIRECTIONS(dim, d) {
    tmp.set_direction(d, tmp.in_direction(d) + num_direction(d) * inva);
  }
  return tmp; // High corner
}

void grid_volume::print() const {
  LOOP_OVER_DIRECTIONS(dim, d) {
    printf("%s =%5g - %5g (%5g) \t", direction_name(d), origin.in_direction(d),
           origin.in_direction(d) + num_direction(d) / a, num_direction(d) / a);
  }
  printf("\n");
}

bool grid_volume::intersect_with(const grid_volume &vol_in, grid_volume *intersection,
                                 grid_volume *others, int *num_others) const {
  int temp_num[3] = {0, 0, 0};
  ivec new_io(dim);
  LOOP_OVER_DIRECTIONS(dim, d) {
    int minval = std::max(little_corner().in_direction(d), vol_in.little_corner().in_direction(d));
    int maxval = std::min(big_corner().in_direction(d), vol_in.big_corner().in_direction(d));
    if (minval >= maxval) return false;
    temp_num[d % 3] = (maxval - minval) / 2;
    new_io.set_direction(d, minval);
  }
  if (intersection != NULL) {
    *intersection = grid_volume(dim, a, temp_num[0], temp_num[1],
                                temp_num[2]); // fix me : ugly, need new constructor
    intersection->set_origin(new_io);
  }
  if (others != NULL) {
    int counter = 0;
    grid_volume vol_containing = *this;
    LOOP_OVER_DIRECTIONS(dim, d) {
      if (vol_containing.little_corner().in_direction(d) < vol_in.little_corner().in_direction(d)) {
        // shave off lower slice from vol_containing and add it to others
        grid_volume other = vol_containing;
        const int thick = (vol_in.little_corner().in_direction(d) -
                           vol_containing.little_corner().in_direction(d)) /
                          2;
        other.set_num_direction(d, thick);
        others[counter] = other;
        counter++;
        vol_containing.shift_origin(d, thick * 2);
        vol_containing.set_num_direction(d, vol_containing.num_direction(d) - thick);
        if (vol_containing.little_corner().in_direction(d) < vol_in.little_corner().in_direction(d))
          meep::abort("intersect_with: little corners differ by odd integer?");
      }
      if (vol_containing.big_corner().in_direction(d) > vol_in.big_corner().in_direction(d)) {
        // shave off upper slice from vol_containing and add it to others
        grid_volume other = vol_containing;
        const int thick =
            (vol_containing.big_corner().in_direction(d) - vol_in.big_corner().in_direction(d)) / 2;
        other.set_num_direction(d, thick);
        other.shift_origin(d, (vol_containing.num_direction(d) - thick) * 2);
        others[counter] = other;
        counter++;
        vol_containing.set_num_direction(d, vol_containing.num_direction(d) - thick);
        if (vol_containing.big_corner().in_direction(d) < vol_in.big_corner().in_direction(d))
          meep::abort("intersect_with: big corners differ by odd integer?");
      }
    }
    *num_others = counter;

    size_t initial_points = 1;
    LOOP_OVER_DIRECTIONS(dim, d) { initial_points *= num_direction(d); }
    size_t final_points, temp = 1;
    LOOP_OVER_DIRECTIONS(dim, d) { temp *= intersection->num_direction(d); }
    final_points = temp;
    for (int j = 0; j < *num_others; j++) {
      temp = 1;
      LOOP_OVER_DIRECTIONS(dim, d) { temp *= others[j].num_direction(d); }
      final_points += temp;
    }
    if (initial_points != final_points)
      meep::abort("intersect_with: initial_points != final_points,  %zd, %zd\n", initial_points,
                  final_points);
  }
  return true;
}

vec grid_volume::loc_at_resolution(ptrdiff_t index, double res) const {
  vec where = origin;
  for (int dd = X; dd <= R; dd++) {
    const direction d = (direction)dd;
    if (has_boundary(High, d)) {
      const double dist = boundary_location(High, d) - boundary_location(Low, d);
      const int nhere = std::max(1, (int)floor(dist * res + 0.5));
      where.set_direction(d, origin.in_direction(d) + ((index % nhere) + 0.5) * (1.0 / res));
      index /= nhere;
    }
  }
  return where;
}

size_t grid_volume::ntot_at_resolution(double res) const {
  size_t mytot = 1;
  for (int d = X; d <= R; d++)
    if (has_boundary(High, (direction)d)) {
      const double dist =
          boundary_location(High, (direction)d) - boundary_location(Low, (direction)d);
      mytot *= std::max(size_t(1), (size_t)(dist * res + 0.5));
    }
  return mytot;
}

vec grid_volume::loc(component c, ptrdiff_t ind) const { return operator[](iloc(c, ind)); }

ivec grid_volume::iloc(component c, ptrdiff_t ind) const {
  ivec out(dim);
  LOOP_OVER_DIRECTIONS(dim, d) {
    ptrdiff_t ind_over_stride = ind / stride(d);
    while (ind_over_stride < 0)
      ind_over_stride += num_direction(d) + 1;
    out.set_direction(d, 2 * (ind_over_stride % (num_direction(d) + 1)));
  }
  return out + iyee_shift(c) + io;
}

size_t grid_volume::surface_area() const {
  switch (dim) {
    case Dcyl: return 2 * (nr() + nz());
    case D1: return 2;
    case D2: return 2 * (nx() + ny());
    case D3: return 2 * (nx() * ny() + nx() * nz() + ny() * nz());
  }
  return 0; // This is never reached.
}

vec grid_volume::dr() const {
  switch (dim) {
    case Dcyl: return veccyl(inva, 0.0);
    case D1:
    case D2:
    case D3: meep::abort("Error in dr\n");
  }
  return vec(0); // This is never reached.
}

vec grid_volume::dx() const {
  switch (dim) {
    case D3: return vec(inva, 0, 0);
    case D2: return vec(inva, 0);
    case D1:
    case Dcyl: meep::abort("Error in dx.\n");
  }
  return vec(0); // This is never reached.
}

vec grid_volume::dy() const {
  switch (dim) {
    case D3: return vec(0, inva, 0);
    case D2: return vec(0, inva);
    case D1:
    case Dcyl: meep::abort("Error in dy.\n");
  }
  return vec(0); // This is never reached.
}

vec grid_volume::dz() const {
  switch (dim) {
    case Dcyl: return veccyl(0.0, inva);
    case D3: return vec(0, 0, inva);
    case D1: return vec(inva);
    case D2: meep::abort("dz doesn't exist in 2D\n");
  }
  return vec(0); // This is never reached.
}

grid_volume volone(double zsize, double a) {
  if (!isinteger(zsize * a))
    master_printf_stderr("Warning: grid volume is not an integer number of pixels; cell size will "
                         "be rounded to nearest pixel.\n");
  return grid_volume(D1, a, 0, 0, (int)(zsize * a + 0.5));
}

grid_volume voltwo(double xsize, double ysize, double a) {
  if (!isinteger(xsize * a) || !isinteger(ysize * a))
    master_printf_stderr("Warning: grid volume is not an integer number of pixels; cell size will "
                         "be rounded to nearest pixel.\n");
  return grid_volume(D2, a, (xsize == 0) ? 1 : (int)(xsize * a + 0.5),
                     (ysize == 0) ? 1 : (int)(ysize * a + 0.5), 0);
}

grid_volume vol1d(double zsize, double a) { return volone(zsize, a); }

grid_volume vol2d(double xsize, double ysize, double a) { return voltwo(xsize, ysize, a); }

grid_volume vol3d(double xsize, double ysize, double zsize, double a) {
  if (!isinteger(xsize * a) || !isinteger(ysize * a) || !isinteger(zsize * a))
    master_printf_stderr("Warning: grid volume is not an integer number of pixels; cell size will "
                         "be rounded to nearest pixel.\n");
  return grid_volume(D3, a, (xsize == 0) ? 1 : (int)(xsize * a + 0.5),
                     (ysize == 0) ? 1 : (int)(ysize * a + 0.5),
                     (zsize == 0) ? 1 : (int)(zsize * a + 0.5));
}

grid_volume volcyl(double rsize, double zsize, double a) {
  if (!isinteger(rsize) || !isinteger(zsize))
    master_printf_stderr("Warning: grid volume is not an integer number of pixels; cell size will "
                         "be rounded to nearest pixel.\n");

  if (zsize == 0.0)
    return grid_volume(Dcyl, a, (int)(rsize * a + 0.5), 0, 1);
  else
    return grid_volume(Dcyl, a, (int)(rsize * a + 0.5), 0, (int)(zsize * a + 0.5));
}

double grid_volume::get_cost() const {
  geom_box box = meep_geom::gv2box(surroundings());
  meep_geom::fragment_stats fstats(box);
  fstats.compute();
  return fstats.cost();
}

// return complex(left cost, right cost).  Should really be a tuple, but we don't want to require
// C++11? yet?
std::complex<double> grid_volume::get_split_costs(direction d, int split_point,
                                                  bool fragment_cost) const {
  double left_cost = 0, right_cost = 0;
  if (split_point > 0) {
    grid_volume v_left = *this;
    v_left.set_num_direction(d, split_point);
    left_cost = fragment_cost ? v_left.get_cost() : v_left.nowned_min();
  }
  if (split_point < num_direction(d)) {
    grid_volume v_right = *this;
    v_right.set_num_direction(d, num_direction(d) - split_point);
    v_right.shift_origin(d, split_point * 2);
    right_cost = fragment_cost ? v_right.get_cost() : v_right.nowned_min();
  }
  return std::complex<double>(left_cost, right_cost);
}

static double cost_diff(int desired_chunks, std::complex<double> costs) {
  double left_cost = real(costs) / (desired_chunks / 2);
  double right_cost = imag(costs) / (desired_chunks - desired_chunks / 2);
  return right_cost - left_cost;
}

void grid_volume::tile_split(int &best_split_point, direction &best_split_direction) const {
  const size_t ntot_thresh = 10;
  if (ntot() < ntot_thresh) {
    best_split_point = 0;
    best_split_direction = NO_DIRECTION;
  }
  else if (nx() > 1) {
    best_split_point = nx() / 2;
    best_split_direction = X;
  }
  else if (ny() > 1) {
    best_split_point = ny() / 2;
    best_split_direction = Y;
  }
  else {
    best_split_point = nz() / 2;
    best_split_direction = Z;
  }
}

void grid_volume::find_best_split(int desired_chunks, bool fragment_cost, int &best_split_point,
                                  direction &best_split_direction,
                                  double &left_effort_fraction) const {
  if (size_t(desired_chunks) > nowned_min()) {
    meep::abort("Cannot split %zd grid points into %d parts\n", nowned_min(), desired_chunks);
  }

  left_effort_fraction = 0;
  best_split_point = 0;
  best_split_direction = NO_DIRECTION;
  if (desired_chunks == 1) return;

  direction longest_axis = NO_DIRECTION;
  int num_in_longest_axis = 0;
  LOOP_OVER_DIRECTIONS(dim, d) {
    if (num_direction(d) > num_in_longest_axis) {
      longest_axis = d;
      num_in_longest_axis = num_direction(d);
    }
  }

  double best_split_measure = 1e20;
  LOOP_OVER_DIRECTIONS(dim, d) {
    int first = 0, last = num_direction(d);
    while (first < last) { // bisection search for balanced splitting
      int mid = (first + last) / 2;
      double mid_diff = cost_diff(desired_chunks, get_split_costs(d, mid, fragment_cost));
      if (mid_diff > 0) {
        if (first == mid) break;
        first = mid;
      }
      else if (mid_diff < 0)
        last = mid;
      else
        break;
    }
    int split_point = (first + last) / 2;
    std::complex<double> costs = get_split_costs(d, split_point, fragment_cost);
    double left_cost = real(costs), right_cost = imag(costs);
    double total_cost = left_cost + right_cost;
    double split_measure = std::max(left_cost / (desired_chunks / 2),
                                    right_cost / (desired_chunks - (desired_chunks / 2)));
    // Give a 30% preference to the longest axis, as a heuristic to prefer lower communication costs
    // when the split_measure is somewhat close.   TODO: use a data-driven communication cost
    // function.
    if (d == longest_axis) split_measure *= 0.7;
    if (split_measure < best_split_measure) {
      best_split_measure = split_measure;
      best_split_point = split_point;
      best_split_direction = d;
      left_effort_fraction = left_cost / total_cost;
    }
  }
}

grid_volume grid_volume::split_at_fraction(bool side_high, int split_pt, int split_dir) const {
  if (dim == Dcyl) split_dir %= 3;
  grid_volume retval(dim, a, 1, 1, 1);
  for (int i = 0; i < 3; i++)
    retval.num[i] = num[i];
  if (split_pt >= num[split_dir]) meep::abort("Aaack bad bug in split_at_fraction.\n");
  direction d = (direction)split_dir;
  if (dim == Dcyl && d == X) d = R;
  retval.set_origin(io);
  if (side_high) retval.shift_origin(d, split_pt * 2);

  if (side_high)
    retval.num[split_dir] -= split_pt;
  else
    retval.num[split_dir] = split_pt;
  retval.num_changed();
  return retval;
}

// Halve the grid_volume for symmetry exploitation...must contain icenter!
grid_volume grid_volume::halve(direction d) const {
  grid_volume retval(*this);
  retval.set_num_direction(d, 1 + (big_corner().in_direction(d) - icenter().in_direction(d)) / 2);
  retval.set_origin(d, icenter().in_direction(d) - 2);
  return retval;
}

grid_volume grid_volume::pad(direction d) const {
  grid_volume gv(*this);
  gv.pad_self(d);
  return gv;
}

void grid_volume::pad_self(direction d) {
  num[d % 3] += 2; // Pad in both directions by one grid point.
  num_changed();
  shift_origin(d, -2);
}

ivec grid_volume::icenter() const {
  /* Find the center of the user's cell.  This will be used as the
     symmetry point, and therefore icenter-io must be *even*
     in all components in order that rotations preserve the Yee lattice. */
  switch (dim) {
    case D1: return io + ivec(nz()).round_down_to_even();
    case D2: return io + ivec(nx(), ny()).round_down_to_even();
    case D3: return io + ivec(nx(), ny(), nz()).round_down_to_even();
    case Dcyl: return io + iveccyl(0, nz()).round_down_to_even();
  }
  meep::abort("Can't do symmetry with these dimensions.\n");
  return ivec(0); // This is never reached.
}

vec grid_volume::center() const { return operator[](icenter()); }

symmetry rotate4(direction axis, const grid_volume &gv) {
  symmetry s = identity();
  if (axis > 2) meep::abort("Can only rotate4 in 2D or 3D.\n");
  s.g = 4;
  FOR_DIRECTIONS(d) {
    s.S[d].d = d;
    s.S[d].flipped = false;
  }
  s.S[(axis + 1) % 3].d = (direction)((axis + 2) % 3);
  s.S[(axis + 1) % 3].flipped = true;
  s.S[(axis + 2) % 3].d = (direction)((axis + 1) % 3);
  s.symmetry_point = gv.center();
  s.i_symmetry_point = gv.icenter();
  return s;
}

symmetry rotate2(direction axis, const grid_volume &gv) {
  symmetry s = identity();
  if (axis > 2) meep::abort("Can only rotate2 in 2D or 3D.\n");
  s.g = 2;
  s.S[(axis + 1) % 3].flipped = true;
  s.S[(axis + 2) % 3].flipped = true;
  s.symmetry_point = gv.center();
  s.i_symmetry_point = gv.icenter();
  return s;
}

symmetry mirror(direction axis, const grid_volume &gv) {
  symmetry s = identity();
  s.g = 2;
  s.S[axis].flipped = true;
  s.symmetry_point = gv.center();
  s.i_symmetry_point = gv.icenter();
  return s;
}

symmetry r_to_minus_r_symmetry(double m) {
  symmetry s = identity();
  s.g = 2;
  s.S[R].flipped = true;
  s.S[P].flipped = true;
  s.symmetry_point = zero_vec(Dcyl);
  s.i_symmetry_point = zero_ivec(Dcyl);
  if (m == int(m)) // phase is purely real (+/- 1) when m an integer
    s.ph = (int(m) & 1) ? -1.0 : 1.0;
  else
    s.ph = polar(1.0, m * pi); // general case
  return s;
}

symmetry identity() { return symmetry(); }

symmetry::symmetry() {
  g = 1;
  ph = 1.0;
  FOR_DIRECTIONS(d) {
    S[d].d = d;
    S[d].flipped = false;
  }
  next = NULL;
}

symmetry::symmetry(const symmetry &s) {
  g = s.g;
  FOR_DIRECTIONS(d) {
    S[d].d = s.S[d].d;
    S[d].flipped = s.S[d].flipped;
  }
  ph = s.ph;
  symmetry_point = s.symmetry_point;
  i_symmetry_point = s.i_symmetry_point;
  if (s.next)
    next = new symmetry(*s.next);
  else
    next = NULL;
}

void symmetry::operator=(const symmetry &s) {
  g = s.g;
  FOR_DIRECTIONS(d) {
    S[d].d = s.S[d].d;
    S[d].flipped = s.S[d].flipped;
  }
  ph = s.ph;
  symmetry_point = s.symmetry_point;
  i_symmetry_point = s.i_symmetry_point;
  if (s.next)
    next = new symmetry(*s.next);
  else
    next = NULL;
}

bool symmetry::operator==(const symmetry &sym) const {
  int gtot = multiplicity();
  if (gtot != sym.multiplicity()) return false;
  for (int sn = 1; sn < gtot; ++sn)
    FOR_DIRECTIONS(d)
  if (transform(d, sn) != sym.transform(d, sn)) return false;
  return true;
}

symmetry::~symmetry() { delete next; }

int symmetry::multiplicity() const {
  if (next)
    return g * next->multiplicity();
  else
    return g;
}

int symmetry::multiplicity(ivec &x) const {
  int m = multiplicity();
  for (int n = 1; n < m; ++n) {
    if (transform(x, n) == x) return n;
  }
  return m;
}

symmetry symmetry::operator+(const symmetry &b) const {
  // The following optimization ignores identity when adding symmetries
  // together.  This is important because identity has an undefined
  // symmetry point.
  if (multiplicity() == 1)
    return b;
  else if (b.multiplicity() == 1)
    return *this;
  symmetry s = *this;
  symmetry *sn = &s;
  for (; sn->next; sn = sn->next)
    ;
  sn->next = new symmetry(b);
  return s;
}

symmetry symmetry::operator*(complex<double> p) const {
  symmetry s = *this;
  s.ph *= p;
  return s;
}

signed_direction signed_direction::operator*(complex<double> p) {
  signed_direction sd = *this;
  sd.phase *= p;
  return sd;
}

signed_direction symmetry::transform(direction d, int n) const {
  // Returns transformed direction + phase/flip; -n indicates inverse transform
  if (n == 0 || d == NO_DIRECTION) return signed_direction(d);
  int nme, nrest;
  if (n < 0) {
    nme = (g - (-n) % g) % g;
    nrest = -((-n) / g);
  }
  else {
    nme = n % g;
    nrest = n / g;
  }
  if (nme == 0) {
    if (nrest == 0)
      return signed_direction(d);
    else
      return next->transform(d, nrest);
  }
  else {
    signed_direction sd;
    if (nme == 1) sd = S[d];
    if (S[d].flipped)
      sd = flip(transform(S[d].d, nme - 1));
    else
      sd = transform(S[d].d, nme - 1);

    if (next && nrest) {
      if (sd.flipped)
        return flip(next->transform(sd.d, nrest)) * ph;
      else
        return next->transform(sd.d, nrest) * ph;
    }
    else { return sd * ph; }
  }
}

ivec symmetry::transform(const ivec &ov, int n) const {
  if (n == 0) return ov;
  ivec out = ov;
  LOOP_OVER_DIRECTIONS(ov.dim, d) {
    const signed_direction s = transform(d, n);
    const int sp_d = i_symmetry_point.in_direction(d);
    const int sp_sd = i_symmetry_point.in_direction(s.d);
    const int delta = ov.in_direction(d) - sp_d;
    if (s.flipped)
      out.set_direction(s.d, sp_sd - delta);
    else
      out.set_direction(s.d, sp_sd + delta);
  }
  return out;
}

ivec symmetry::transform_unshifted(const ivec &ov, int n) const {
  if (n == 0) return ov;
  ivec out(ov.dim);
  LOOP_OVER_DIRECTIONS(ov.dim, d) {
    const signed_direction s = transform(d, n);
    if (s.flipped)
      out.set_direction(s.d, -ov.in_direction(d));
    else
      out.set_direction(s.d, ov.in_direction(d));
  }
  return out;
}

vec symmetry::transform(const vec &ov, int n) const {
  if (n == 0) return ov;
  vec delta = ov;
  LOOP_OVER_DIRECTIONS(ov.dim, d) {
    const signed_direction s = transform(d, n);
    double deltad = ov.in_direction(d) - symmetry_point.in_direction(d);
    if (s.flipped)
      delta.set_direction(s.d, -deltad);
    else
      delta.set_direction(s.d, deltad);
  }
  return symmetry_point + delta;
}

volume symmetry::transform(const volume &v, int n) const {
  return volume(transform(v.get_min_corner(), n), transform(v.get_max_corner(), n));
}

component symmetry::transform(component c, int n) const {
  return direction_component(c, transform(component_direction(c), n).d);
}

derived_component symmetry::transform(derived_component c, int n) const {
  return direction_component(c, transform(component_direction(c), n).d);
}

int symmetry::transform(int c, int n) const {
  return (is_derived(c) ? int(transform(derived_component(c), n))
                        : int(transform(component(c), n)));
}

complex<double> symmetry::phase_shift(component c, int n) const {
  if (c == Dielectric || c == Permeability) return 1.0;
  complex<double> phase = transform(component_direction(c), n).phase;
  // flip tells us if we need to flip the sign.  For vectors (E), it is
  // just this simple:
  bool flip = transform(component_direction(c), n).flipped;
  if (is_magnetic(c) || is_B(c)) {
    // Because H is a pseudovector, here we have to figure out if the
    // transformation changes the handedness of the basis.
    bool have_one = false, have_two = false;
    FOR_DIRECTIONS(d) {
      if (transform(d, n).flipped) flip = !flip;
      int shift = (transform(d, n).d - d + 6) % 3;
      if (shift == 1) have_one = true;
      if (shift == 2) have_two = true;
    }
    if (have_one && have_two) flip = !flip;
  }
  if (flip)
    return -phase;
  else
    return phase;
}

complex<double> symmetry::phase_shift(derived_component c, int n) const {
  if (is_poynting(c)) {
    signed_direction ds = transform(component_direction(c), n);
    complex<double> ph = conj(ds.phase) * ds.phase; // E x H gets |phase|^2
    return (ds.flipped ? -ph : ph);
  }
  else /* energy density */
    return 1.0;
}

complex<double> symmetry::phase_shift(int c, int n) const {
  return (is_derived(c) ? phase_shift(derived_component(c), n) : phase_shift(component(c), n));
}

bool symmetry::is_primitive(const ivec &p) const {
  // This is only correct if p is somewhere on the yee lattice.
  if (multiplicity() == 1) return true;
  for (int i = 1; i < multiplicity(); i++) {
    const ivec pp = transform(p, i);
    switch (p.dim) {
      case D2:
        if (pp.x() + pp.y() < p.x() + p.y()) return false;
        if (pp.x() + pp.y() == p.x() + p.y() && p.y() > p.x() && pp.y() <= pp.x()) return false;
        break;
      case D3:
        if (pp.x() + pp.y() + pp.z() < p.x() + p.y() + p.z()) return false;
        if (pp.x() + pp.y() + pp.z() == p.x() + p.y() + p.z() &&
            pp.x() + pp.y() - pp.z() < p.x() + p.y() - p.z())
          return false;
        if (pp.x() + pp.y() + pp.z() == p.x() + p.y() + p.z() &&
            pp.x() + pp.y() - pp.z() == p.x() + p.y() - p.z() &&
            pp.x() - pp.y() - pp.z() < p.x() - p.y() - p.z())
          return false;
        break;
      case D1:
      case Dcyl:
        if (pp.z() < p.z()) return false;
        break;
    }
  }
  return true;
}

/* given a list of geometric volumes, produce a new list with appropriate
   weights that is minimized according to the symmetry.  */
volume_list *symmetry::reduce(const volume_list *gl) const {
  volume_list *glnew = 0;
  for (const volume_list *g = gl; g; g = g->next) {
    int sn;
    for (sn = 0; sn < multiplicity(); ++sn) {
      volume gS(transform(g->v, sn));
      int cS = transform(g->c, sn);
      volume_list *gn;
      for (gn = glnew; gn; gn = gn->next)
        if (gn->c == cS && gn->v.round_float() == gS.round_float()) break;
      if (gn) { // found a match
        gn->weight += g->weight * phase_shift(g->c, sn);
        break;
      }
    }
    if (sn == multiplicity() && g->weight != 0.0) { // no match, add to glnew
      volume_list *gn = new volume_list(g->v, g->c, g->weight, glnew);
      glnew = gn;
    }
  }

  // reduce v's redundant with themselves & delete elements with zero weight:
  volume_list *gprev = 0, *g = glnew;
  while (g) {
    // first, see if g->v is redundant with itself
    bool halve[5] = {false, false, false, false, false};
    complex<double> weight = g->weight;
    for (int sn = 1; sn < multiplicity(); ++sn)
      if (g->c == transform(g->c, sn) && g->v.round_float() == transform(g->v, sn).round_float()) {
        LOOP_OVER_DIRECTIONS(g->v.dim, d) {
          if (transform(d, sn).flipped) {
            halve[d] = true;
            break;
          }
        }
        g->weight += weight * phase_shift(g->c, sn);
      }
    LOOP_OVER_DIRECTIONS(g->v.dim, d) {
      if (halve[d])
        g->v.set_direction_max(d, g->v.in_direction_min(d) + 0.5 * g->v.in_direction(d));
    }
    // now, delete it if it has zero weight
    if (g->weight == 0.0) {
      if (gprev)
        gprev->next = g->next;
      else // g == glnew
        glnew = g->next;
      g->next = 0; // necessary so that g->next is not deleted recursively
      delete g;
      g = gprev ? gprev->next : glnew;
    }
    else
      g = (gprev = g)->next;
  }

  return glnew;
}

/***************************************************************************/

static double poynting_fun(const complex<realnum> *fields, const vec &loc, void *data_) {
  (void)loc;   // unused
  (void)data_; // unused
  return (real(conj(cdouble(fields[0])) * cdouble(fields[1])) -
          real(conj(cdouble(fields[2])) * cdouble(fields[3])));
}

static double energy_fun(const complex<realnum> *fields, const vec &loc, void *data_) {
  (void)loc; // unused
  int nfields = *((int *)data_) / 2;
  double sum = 0;
  for (int k = 0; k < nfields; ++k)
    sum += real(conj(cdouble(fields[2 * k])) * cdouble(fields[2 * k + 1]));
  return sum * 0.5;
}

field_rfunction derived_component_func(derived_component c, const grid_volume &gv, int &nfields,
                                       component cs[12]) {
  switch (c) {
    case Sx:
    case Sy:
    case Sz:
    case Sr:
    case Sp:
      switch (c) {
        case Sx:
          cs[0] = Ey;
          cs[1] = Hz;
          break;
        case Sy:
          cs[0] = Ez;
          cs[1] = Hx;
          break;
        case Sz:
          cs[0] = Ex;
          cs[1] = Hy;
          break;
        case Sr:
          cs[0] = Ep;
          cs[1] = Hz;
          break;
        case Sp:
          cs[0] = Ez;
          cs[1] = Hr;
          break;
        default: break; // never reached
      }
      nfields = 4;
      cs[2] = direction_component(Ex, component_direction(cs[1]));
      cs[3] = direction_component(Hx, component_direction(cs[0]));
      return poynting_fun;

    case EnergyDensity:
    case D_EnergyDensity:
    case H_EnergyDensity:
      nfields = 0;
      if (c != H_EnergyDensity) FOR_ELECTRIC_COMPONENTS(c0) {
          if (gv.has_field(c0)) {
            cs[nfields++] = c0;
            cs[nfields++] = direction_component(Dx, component_direction(c0));
          }
        }
      if (c != D_EnergyDensity) FOR_MAGNETIC_COMPONENTS(c0) {
          if (gv.has_field(c0)) {
            cs[nfields++] = c0;
            cs[nfields++] = direction_component(Bx, component_direction(c0));
          }
        }
      if (nfields > 12) meep::abort("too many field components");
      return energy_fun;

    default: meep::abort("unknown derived_component in derived_component_func");
  }
  return 0;
}

/***************************************************************************/
/* utility methods for pretty-printing. may be called with no arguments,   */
/* in which case static internal buffers are used; NUMBUFS defines the     */
/* number of str() calls with no arguments that may be appear              */
/* simultaneously as e.g. arguments to a single invocation of printf().    */
/***************************************************************************/
#define BUFLEN 100
#define NUMBUFS 10
const char *ivec::str(char *buffer, size_t buflen) {
  static char bufring[NUMBUFS][BUFLEN];
  static int nbuf = 0;
  if (buffer == 0) {
    buffer = bufring[nbuf];
    buflen = BUFLEN;
    nbuf = (nbuf + 1) % NUMBUFS;
  }
  if (dim == Dcyl)
    snprintf(buffer, buflen, "{%i,%i}", t[R], t[Z]);
  else
    snprintf(buffer, buflen, "{%i,%i,%i}", t[X], t[Y], t[Z]);
  return buffer;
}

const char *vec::str(char *buffer, size_t buflen) {
  static char bufring[NUMBUFS][BUFLEN];
  static int nbuf = 0;
  if (buffer == 0) {
    buffer = bufring[nbuf];
    buflen = BUFLEN;
    nbuf = (nbuf + 1) % NUMBUFS;
  }
  if (dim == Dcyl)
    snprintf(buffer, buflen, "{%f,%f}", t[R], t[Z]);
  else
    snprintf(buffer, buflen, "{%f,%f,%f}", t[X], t[Y], t[Z]);
  return buffer;
}

const char *volume::str(char *buffer, size_t buflen) {
  static char sbuf[1024]; // TODO: Use a bufring like the above?
  if (buffer == 0) {
    buffer = sbuf;
    buflen = sizeof(sbuf);
  }
  snprintf(buffer, buflen, "min_corner:%s, max_corner:%s", min_corner.str(), max_corner.str());
  return buffer;
}

const char *grid_volume::str(char *buffer, size_t buflen) {
  static char sbuf[1024]; // TODO: is this big enough?
  int written = 0;
  if (buffer == 0) {
    buffer = sbuf;
    buflen = sizeof(sbuf);
  }

  written += snprintf(buffer + written, buflen - written,
                      "grid_volume {\n  dim:%s, a:%f, inva:%f, num:{%d, %d, %d}\n",
                      dimension_name(dim), a, inva, num[0], num[1], num[2]);

  // Adapted from the print() method
  LOOP_OVER_DIRECTIONS(dim, d) {
    written += snprintf(buffer + written, buflen - written, "  %s =%5g - %5g (%5g) \t",
                        direction_name(d), origin.in_direction(d),
                        origin.in_direction(d) + num_direction(d) / a, num_direction(d) / a);
    if (buflen - written <= 0) break;
  }
  snprintf(buffer + written, buflen - written, "\n}");
  return buffer;
}

/********************************************************************/
/********************************************************************/
/********************************************************************/
grid_volume grid_volume::subvolume(ivec is, ivec ie, component c) {
  if (!(contains(is) && contains(ie))) meep::abort("invalid extents in subvolume");
  grid_volume sub;
  sub.dim = dim;
  sub.a = a;
  sub.inva = inva;
  sub.init_subvolume(is, ie, c);
  return sub;
}

void grid_volume::init_subvolume(ivec is, ivec ie, component c) {
  ivec origin(dim, 0);
  for (int i = 0; i < 3; i++)
    num[i] = 0;
  LOOP_OVER_DIRECTIONS(dim, d) {
    num[(int)d] = (ie - is).in_direction(d) / 2;
    origin.set_direction(d, is.in_direction(d) - iyee_shift(c).in_direction(d));
  }
  num_changed();
  // center_origin();
  set_origin(origin);
}

} // namespace meep
