/* Copyright (C) 2003 Massachusetts Institute of Technology
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

namespace meep {

static inline double yee_to_lattice(int n, double a, double inva=0.0) {
  if (inva == 0.0) inva = 1.0/a;
  return n*(0.5*inva);
}

static inline int lattice_to_yee(double x, double a) {
  return (int)floor(x*(2.0*a) + 0.5);
}

ivec volume::round_vec(const vec &p) const {
  ivec result(dim);
  LOOP_OVER_DIRECTIONS(dim, d)
    result.set_direction(d, lattice_to_yee(p.in_direction(d),a));
  return result;
}

void volume::update_io() {
  io = round_vec(origin);
}

static inline double rtl(double x, double a, double inva=0.0) {
  // Rounds to a value somewhere on the yee lattice.
  return ((int)floor(x*(2.0*a) + 0.5))*(0.5*inva);
}

static inline vec round_to_lattice(const vec &p, double a, double inva=0.0) {
  // Rounds to a value somewhere on the yee lattice.
  vec result(p.dim);
  LOOP_OVER_DIRECTIONS(p.dim, d) result.set_direction(d, rtl(p.in_direction(d),a,inva));
  return result;
}

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
  case Dielectric: return "eps";
  }
  return "Error in component_name";
}

const char *component_name(derived_component c) {
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
  return (is_derived(c) ? component_name(derived_component(c))
	  : component_name(component(c)));
}

void vec::print(file *f) const {
  if (dim == Dcyl) {
    i_fprintf(f, "%g %g 0", r(), z());
  } else if (dim == D3  ) {
    i_fprintf(f, "%g %g %g", x(), y(), z());
  } else if (dim == D2  ) {
    i_fprintf(f, "%g %g 0", x(), y());
  } else if (dim == D1  ) {
    i_fprintf(f, "0 0 %g", z());
  } else  {
    printf("I don't know how to print in this dimension!\n");
    i_fprintf(f, "I don't know how to print in this dimension!\n");
  }
}

void ivec::print(file *f) const {
  if (dim == Dcyl) {
    i_fprintf(f, "%d %d 0", r(), z());
  } else if (dim == D3  ) {
    i_fprintf(f, "%d %d %d", x(), y(), z());
  } else if (dim == D2  ) {
    i_fprintf(f, "%d %d 0", x(), y());
  } else if (dim == D1  ) {
    i_fprintf(f, "0 0 %d", z());
  } else  {
    printf("I don't know how to print in this dimension!\n");
    i_fprintf(f, "I don't know how to print in this dimension!\n");
  }
}

vec min(const vec &vec1, const vec &vec2) {
  vec m(vec1.dim);
  LOOP_OVER_DIRECTIONS(vec1.dim, d)
    m.set_direction(d, min(vec1.in_direction(d), vec2.in_direction(d)));
  return m;
}

vec max(const vec &vec1, const vec &vec2) {
  vec m(vec1.dim);
  LOOP_OVER_DIRECTIONS(vec1.dim, d)
    m.set_direction(d, max(vec1.in_direction(d), vec2.in_direction(d)));
  return m;
}

ivec min(const ivec &ivec1, const ivec &ivec2) {
  ivec m(ivec1.dim);
  LOOP_OVER_DIRECTIONS(ivec1.dim, d)
    m.set_direction(d, min(ivec1.in_direction(d), ivec2.in_direction(d)));
  return m;
}

ivec max(const ivec &ivec1, const ivec &ivec2) {
  ivec m(ivec1.dim);
  LOOP_OVER_DIRECTIONS(ivec1.dim, d)
    m.set_direction(d, max(ivec1.in_direction(d), ivec2.in_direction(d)));
  return m;
}

geometric_volume::geometric_volume(const vec &vec1, const vec &vec2) {
  min_corner = min(vec1, vec2);
  max_corner = max(vec1, vec2);
  dim = vec1.dim; 
}

geometric_volume::geometric_volume(const vec &pt) {
  dim = pt.dim; 
  min_corner = pt;
  max_corner = pt;
}

double geometric_volume::computational_volume() {
  double vol = 1.0; 
  LOOP_OVER_DIRECTIONS(dim,d) vol *= (in_direction_max(d) - in_direction_min(d));
  return vol;
}

double geometric_volume::full_volume() const {
  double vol = 1.0; 
  LOOP_OVER_DIRECTIONS(dim, d) vol *= (in_direction_max(d) - in_direction_min(d));
  if (dim == Dcyl) vol *= pi * (in_direction_max(R) + in_direction_min(R));
  return vol;
}

geometric_volume geometric_volume::intersect_with(const geometric_volume &a) const {
  if (a.dim != dim) abort("Can't intersect volumes of dissimilar dimensions.\n");
  geometric_volume result(dim);
  LOOP_OVER_DIRECTIONS(dim, d) {
    double minval = max(in_direction_min(d), a.in_direction_min(d));
    double maxval = min(in_direction_max(d), a.in_direction_max(d));
    if (minval > maxval)
      return geometric_volume(zero_vec(dim), zero_vec(dim));
    result.set_direction_min(d, minval);
    result.set_direction_max(d, maxval);
  }
  return result;
}

bool geometric_volume::intersects(const geometric_volume &a) const {
  if (a.dim != dim) abort("Can't intersect volumes of dissimilar dimensions.\n");
  LOOP_OVER_DIRECTIONS(dim, d) {
    double minval = max(in_direction_min(d), a.in_direction_min(d));
    double maxval = min(in_direction_max(d), a.in_direction_max(d));
    if (minval > maxval)
      return false;
  }
  return true;
}

// Return normal direction to volume, if the volume is dim-1 dimensional;
// otherwise, return NO_DIRECTION.
direction geometric_volume::normal_direction() const {
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
    if (zx && !zy && !zz) d = X;
    else if (!zx && zy && !zz) d = Y;
    else if (!zx && !zy && zz) d = Z;
    break;
  }
  }
  return d;
}

/* Used for n=0,1,2 nested loops in macros.  We should arrange
   the ordering so that this gives most efficient traversal of
   a field array, where n=2 is the innermost loop. */
static direction yucky_dir(ndim dim, int n) {
  if (dim == Dcyl)
    switch (n) {
    case 0: return P;
    case 1: return R;
    case 2: return Z;
    }
  else if (dim == D2)
    return (direction) ((n + 2) % 3); /* n = 0,1,2 gives Z, X, Y */
  return (direction) n ;
}

int ivec::yucky_val(int n) const {
  if (has_direction(dim, yucky_dir(dim, n)))
    return in_direction(yucky_dir(dim, n));
  return 0;
}

int volume::yucky_num(int n) const {
  if (has_direction(dim, yucky_dir(dim, n)))
    return num_direction(yucky_dir(dim, n));
  return 1;
}

direction volume::yucky_direction(int n) const {
  return yucky_dir(dim, n);
}

geometric_volume volume::surroundings() const {
  geometric_volume res(dim);
  LOOP_OVER_DIRECTIONS(dim, d) {
    res.set_direction_min(d, operator[](little_corner()).in_direction(d)+0.0/a);
    res.set_direction_max(d, operator[](big_corner()).in_direction(d)-0.0/a);
  }
  return res;
}

void volume::update_ntot() {
  the_ntot = 1;
  LOOP_OVER_DIRECTIONS(dim, d) the_ntot *= num[d%3] + 1;
}

void volume::set_num_direction(direction d, int value) {
  num[d%3] = value; num_changed();
}

volume::volume(ndim td, double ta, int na, int nb, int nc) {
  dim = td; a = ta; inva = 1.0 / ta;
  num[0] = na;
  num[1] = nb;
  num[2] = nc;
  num_changed();
  set_origin(zero_vec(dim));
}

component volume::eps_component() const {
  switch (dim) {
  case D1: return Hy;
  case D2: return Hz;
  case D3: return Dielectric;
  case Dcyl: return Hp;
  }
  abort("Unsupported dimensionality eps.\n");
  return Ex;
}

vec volume::yee_shift(component c) const {
  return operator[](iyee_shift(c));
}

/* Return array offsets to average with a given array location of c in
   order to get c on the "dielectric" grid.  Then, to get the
   dielectric grid point i, you should average c over the four
   locations: i, i+offset1, i+offset2, i+offset1+offset2. 
   (offset2, and possibly offset1, may be zero if only 2 or 1
   locations need to be averaged). */
void volume::yee2diel_offsets(component c, int &offset1, int &offset2) {
  offset1 = offset2 = 0;
  LOOP_OVER_DIRECTIONS(dim,d) {
    if (!iyee_shift(c).in_direction(d)) {
      if (offset2) 
	abort("weird yee shift for component %s", component_name(c));
      if (offset1) offset2 = stride(d);
      else offset1 = stride(d);
    }
  }
}

bool geometric_volume::contains(const vec &p) const {
  LOOP_OVER_DIRECTIONS(dim,d) {
    if (p.in_direction(d) > in_direction_max(d) ||
        p.in_direction(d) < in_direction_min(d)) return false;
  }
  return true;
}

bool volume::contains(const ivec &p) const {
  // containts returns true if the volume has information about this grid
  // point.
  const ivec o = p - io;
  LOOP_OVER_DIRECTIONS(dim, d)
    if (o.in_direction(d) < 0 || o.in_direction(d) >= (num_direction(d)+1)*2)
      return false;
  return true;
}

bool volume::contains(const vec &p) const {
  // containts returns true if the volume has any information in it
  // relevant to the point p.  Basically has is like owns (see below)
  // except it is more lenient, in that more than one lattice may contain a
  // given point.
  const vec o = p - origin;
  LOOP_OVER_DIRECTIONS(dim, d)
    if (o.in_direction(d) < -inva || o.in_direction(d) > num_direction(d)*inva+inva)
      return false;
  return true;
}

/* Compute the corners (cs,ce) of the ib-th boundary for component c,
   returning true if ib is a valid index (ib = 0..#boundaries-1).  The
   boundaries are all the points that are in but not owned by the
   volume, and are a set of *disjoint* regions.  The main purpose of
   this function is currently to support the LOOP_OVER_NOT_OWNED
   macro.  (In the future, it may be used for other
   boundary-element-type computations, too.) */
bool volume::get_boundary_icorners(component c, int ib,
				   ivec *cs, ivec *ce) const {
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

// first "owned" point for c in volume (see also volume::owns)
ivec volume::little_owned_corner(component c) const {
  ivec iloc(little_corner() + one_ivec(dim)*2 - iyee_shift(c));
  if (dim == Dcyl && origin.r() == 0.0 && iloc.r() == 2)
    iloc.set_direction(R, 0);
  return iloc;
}

int volume::nowned(component c) const {
  int n = 1;
  ivec v = big_corner() - little_owned_corner(c);
  LOOP_OVER_DIRECTIONS(dim, d) n *= v.in_direction(d) / 2 + 1;
  return n;
}

bool volume::owns(const ivec &p) const {
  // owns returns true if the point "owned" by this volume, meaning that it
  // is the volume that would timestep the point.
  const ivec o = p - io;
  if (dim == Dcyl) {
    if (origin.r() == 0.0 && o.z() > 0 && o.z() <= nz()*2 &&
        o.r() == 0) return true;
    return o.r() > 0 && o.z() > 0 &&
           o.r() <= nr()*2 && o.z() <= nz()*2;
  } else if (dim == D3) {
    return
      o.x() > 0 && o.x() <= nx()*2 &&
      o.y() > 0 && o.y() <= ny()*2 &&
      o.z() > 0 && o.z() <= nz()*2;
  } else if (dim == D2) {
    return
      o.x() > 0 && o.x() <= nx()*2 &&
      o.y() > 0 && o.y() <= ny()*2;
  } else if (dim == D1) {
    return o.z() > 0 && o.z() <= nz()*2;
  } else {
    abort("Unsupported dimension in owns.\n");
    return false;
  }
}

int volume::has_boundary(boundary_side b,direction d) const {
  switch (dim) {
  case Dcyl: return d == Z || (d == R && b == High);
  case D1: return d == Z;
  case D2: return d == X || d == Y;
  case D3: return d == X || d == Y || d == Z;
  }
  return 0; // This should never be reached.
}

int volume::index(component c, const ivec &p) const {
  const ivec offset = p - io - iyee_shift(c);
  int idx = 0;
  LOOP_OVER_DIRECTIONS(dim,d) idx += offset.in_direction(d)/2*stride(d);
  return idx;
}

void volume::set_strides() {
  for (int d=0;d<5;d++) the_stride0[d] = 0; // Yuck yuck yuck.
  LOOP_OVER_DIRECTIONS(dim,d)
    switch(d) {
    case Z: the_stride0[d] = 1; break;
    case R: the_stride0[d] = nz()+1; break;
    case X: the_stride0[d] = (nz()+1)*(ny() + 1); break;
    case Y: the_stride0[d] = nz() + 1; break;
    case P: break; // There is no phi stride...
    case NO_DIRECTION: break; // no stride here, either
    }
  for (int d=0;d<5;d++) the_stride[d] = the_stride0[d] ? the_stride0[d] : 1;
}

static inline void stupidsort(int *ind, double *w, int l) {
  while (l) {
    if (fabs(w[0]) < 2e-15) {
      w[0] = w[l-1];
      ind[0] = ind[l-1];
      w[l-1] = 0.0;
      ind[l-1] = 0;
    } else {
      w += 1;
      ind += 1;
    }
    l -= 1;
  }
}

static inline void stupidsort(ivec *locs, double *w, int l) {
  while (l) {
    if (fabs(w[0]) < 2e-15) {
      w[0] = w[l-1];
      locs[0] = locs[l-1];
      w[l-1] = 0.0;
      locs[l-1] = 0;
    } else {
      w += 1;
      locs += 1;
    }
    l -= 1;
  }
}

void volume::interpolate(component c, const vec &p,
                         int indices[8], double weights[8]) const {
  ivec locs[8];
  interpolate(c, p, locs, weights);
  for (int i=0;i<8&&weights[i];i++)
    if (!owns(locs[i])) weights[i] = 0.0;
  stupidsort(locs, weights, 8);
  for (int i=0;i<8&&weights[i];i++)
    indices[i] = index(c, locs[i]);
  if (!contains(p) && weights[0]) {
    printf("Error at point %g %g\n", p.r(), p.z());
    printf("Interpolated to point %d %d\n", locs[0].r(), locs[0].z());
    printf("Or in other words... %g %g\n",
           operator[](locs[0]).r(), operator[](locs[0]).z());
    printf("I %s own the interpolated point.\n",
           owns(locs[0])?"actually":"don't");
    print();
    abort("Error made in interpolation of %s--fix this bug!!!\n",
          component_name(c));
  }
  // Throw out out of range indices:
  for (int i=0;i<8&&weights[i];i++)
    if (indices[0] < 0 || indices[0] >= ntot()) weights[i] = 0.0;
  // Stupid very crude code to compactify arrays:
  stupidsort(indices, weights, 8);
  if (!contains(p) && weights[0]) {
    printf("Error at point %g %g\n", p.r(), p.z());
    printf("Interpolated to point %d %d\n", locs[0].r(), locs[0].z());
    print();
    abort("Error made in interpolation of %s--fix this bug!!!\n",
          component_name(c));
  }
}

void volume::interpolate(component c, const vec &pc,
                         ivec locs[8], double weights[8]) const {
  const double SMALL = 1e-13;
  const vec p = (pc - yee_shift(c))*a;
  ivec middle(dim);
  LOOP_OVER_DIRECTIONS(dim,d)
    middle.set_direction(d, ((int) floor(p.in_direction(d)))*2+1);
  middle += iyee_shift(c);
  const vec midv = operator[](middle);
  const vec dv = (pc - midv)*(2*a);
  int already_have = 1;
  for (int i=0;i<8;i++) {
    locs[i] = round_vec(midv);
    weights[i] = 1.0;
  }
  LOOP_OVER_DIRECTIONS(dim,d) {
    for (int i=0;i<already_have;i++) {
      locs[already_have+i] = locs[i];
      weights[already_have+i] = weights[i];
      locs[i].set_direction(d,middle.in_direction(d)-1);
      weights[i] *= 0.5*(1.0-dv.in_direction(d));
      locs[already_have+i].set_direction(d,middle.in_direction(d)+1);
      weights[already_have+i] *= 0.5*(1.0+dv.in_direction(d));
    }
    already_have *= 2;
  }
  for (int i=already_have;i<8;i++) weights[i] = 0.0;
  double total_weight = 0.0;
  for (int i=0;i<already_have;i++) total_weight += weights[i];
  for (int i=0;i<already_have;i++)
    weights[i] += (1.0 - total_weight)*(1.0/already_have);
  for (int i=0;i<already_have;i++) {
    if (weights[i] < 0.0) {
      if (-weights[i] >= SMALL)
        abort("large negative interpolation weight[%d] = %e\n", i, weights[i]);
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
  for (int i=0;i<8&&weights[i];i++)
    if (weights[i] != weights[0]) all_same = false;
  if (all_same) {
    int num_weights = 0;
    for (int i=0;i<8&&weights[i];i++) num_weights++;
    for (int i=0;i<8&&weights[i];i++) weights[i] = 1.0/num_weights;    
  }
}

geometric_volume empty_volume(ndim dim) {
  geometric_volume out(dim);
  LOOP_OVER_DIRECTIONS(dim,d) {
    out.set_direction_max(d,0.0);
    out.set_direction_min(d,0.0);
  }
  return out;
}

geometric_volume volume::dV(const ivec &here) const {
  const double hinva = 0.5*inva;
  const volume &v = *this;
  const vec h = v[here];
  geometric_volume out(dim);
  LOOP_OVER_DIRECTIONS(dim,d) {
    out.set_direction_max(d,h.in_direction(d)+hinva);
    out.set_direction_min(d,h.in_direction(d)-hinva);
  }
  if (dim == Dcyl && here.r() == 0) {
    out.set_direction_min(R,0.0);
  }
  return out;
}

geometric_volume volume::dV(component c, int ind) const {
  if (!owns(iloc(c, ind))) return empty_volume(dim);
  return dV(iloc(c,ind));
}

double volume::xmax() const {
  const double qinva = 0.25*inva;
  return origin.x() + nx()*inva + qinva;
}

double volume::xmin() const {
  const double qinva = 0.25*inva;
  return origin.x() + qinva;
}

double volume::ymax() const {
  const double qinva = 0.25*inva;
  return origin.y() + ny()*inva + qinva;
}

double volume::ymin() const {
  const double qinva = 0.25*inva;
  return origin.y() + qinva;
}

double volume::zmax() const {
  const double qinva = 0.25*inva;
  return origin.z() + nz()*inva + qinva;
}

double volume::zmin() const {
  const double qinva = 0.25*inva;
  return origin.z() + qinva;
}

double volume::rmax() const {
  const double qinva = 0.25*inva;
  if (dim == Dcyl) return origin.r() + nr()*inva + qinva;
  abort("No rmax in these dimensions.\n");
  return 0.0; // This is never reached.
}

double volume::rmin() const {
  const double qinva = 0.25*inva;
  if (dim == Dcyl) {
    if (origin.r() == 0.0) {
      return 0.0;
    } else {
      return origin.r() + qinva;
    }
  }
  abort("No rmin in these dimensions.\n");
  return 0.0; // This is never reached.
}

double vec::project_to_boundary(direction d, double boundary_loc) {
  return fabs(boundary_loc - in_direction(d));
}

double volume::boundary_location(boundary_side b, direction d) const {
  // Returns the location of metallic walls...
  if (b == High) switch (d) {
  case X: return loc(Ez,ntot()-1).x();
  case Y: return loc(Ez,ntot()-1).y();
  case R: return loc(Ep,ntot()-1).r();
  case Z: if (dim == Dcyl) return loc(Ep,ntot()-1).z();
          else return loc(Ex,ntot()-1).z();
  case P: abort("P has no boundary!\n");
  case NO_DIRECTION: abort("NO_DIRECTION has no boundary!\n");
  }
  else switch (d) {
  case X: return loc(Ez,0).x();
  case Y: return loc(Ez,0).y();
  case R: return loc(Ep,0).r();
  case Z: if (dim == Dcyl) return loc(Ep,0).z();
          else return loc(Ex,0).z();
  case P: abort("P has no boundary!\n");
  case NO_DIRECTION: abort("NO_DIRECTION has no boundary!\n");
  }
  return 0.0;
}

ivec volume::big_corner() const {
  switch (dim) {
  case D1: return io + ivec(nz())*2;
  case D2: return io + ivec(nx(),ny())*2;
  case D3: return io + ivec(nx(),ny(),nz())*2;
  case Dcyl: return io + iveccyl(nr(),nz())*2;
  }
  return ivec(0); // This is never reached.
}

vec volume::corner(boundary_side b) const { 
  if (b == Low) return origin; // Low corner
  vec tmp = origin;
  LOOP_OVER_DIRECTIONS(dim, d)
    tmp.set_direction(d, tmp.in_direction(d) + num_direction(d) * inva);
  return tmp; // High corner
}

void volume::print() const {
  LOOP_OVER_DIRECTIONS(dim, d)
    printf("%s =%5g - %5g (%5g) \t", 
      direction_name(d), origin.in_direction(d), 
      origin.in_direction(d)+num_direction(d)/a, num_direction(d)/a); 
  printf("\n");
}

bool volume::intersect_with(const volume &vol_in, volume *intersection, volume *others, int *num_others) const {
  int temp_num[3] = {0,0,0};
  vec origin_shift(dim);
  LOOP_OVER_DIRECTIONS(dim, d) {
    //printf("dim=%s vol_in.io.in_direction(%s) = %d\n", dimension_name(dim), direction_name(d), vol_in.io.in_direction(d));
    int minval = max(io.in_direction(d), vol_in.io.in_direction(d));
    int maxval = min(big_corner().in_direction(d), vol_in.big_corner().in_direction(d));
    if (minval >= maxval)
      return false;
    temp_num[d%3] = (maxval - minval)/2;
    origin_shift.set_direction(d, (minval - io.in_direction(d))/2/a);
  }
  if (intersection != NULL) {
    *intersection = volume(dim, a, temp_num[0], temp_num[1], temp_num[2]); // fix me : ugly, need new constructor
    intersection->set_origin(origin + origin_shift);
  }
  if (others != NULL) {
    int counter = 0;
    volume vol_containing = *this;
    LOOP_OVER_DIRECTIONS(dim, d) {
      if (vol_containing.io.in_direction(d) < vol_in.io.in_direction(d)) {
	// shave off lower slice from vol_containing and add it to others
	volume other = vol_containing;
	const int thick = (vol_in.io.in_direction(d) - vol_containing.io.in_direction(d))/2;
	other.set_num_direction(d, thick);
	others[counter] = other;
	counter++;
	vol_containing.origin_set_direction(d, vol_containing.origin.in_direction(d) + thick/a);
	vol_containing.set_num_direction(d, vol_containing.num_direction(d) - thick);
      }
      if (vol_containing.big_corner().in_direction(d) > vol_in.big_corner().in_direction(d)) {
	// shave off upper slice from vol_containing and add it to others
	volume other = vol_containing;
	const int thick = (vol_containing.big_corner().in_direction(d) - vol_in.big_corner().in_direction(d))/2;
	other.set_num_direction(d, thick);
	other.origin_set_direction(d, vol_containing.origin.in_direction(d) + (vol_containing.num_direction(d) - thick)/a); 
	others[counter] = other;
	counter++;
	vol_containing.set_num_direction(d, vol_containing.num_direction(d) - thick);
      }
    }
    *num_others = counter;
    
    int initial_points = 1;
    LOOP_OVER_DIRECTIONS(dim, d) initial_points *= num_direction(d);
    int final_points , temp = 1;
    LOOP_OVER_DIRECTIONS(dim, d) temp *= intersection->num_direction(d);    
    final_points = temp;
    for (int j=0; j<*num_others; j++) {
      temp = 1;
      LOOP_OVER_DIRECTIONS(dim, d) temp *= others[j].num_direction(d);
      final_points += temp;
    }
    if (initial_points != final_points)
      abort("intersect_with: initial_points != final_points,  %d, %d\n", 
	    initial_points, final_points);
  }
  return true;
}

vec volume::loc_at_resolution(int index, double res) const {
  vec where = origin;
  for (int dd=X;dd<=R;dd++) {
    const direction d = (direction) dd;
    if (has_boundary(High,d)) {
      const double dist = boundary_location(High,d)-boundary_location(Low,d);
      const int nhere = max(1,(int)floor(dist*res+0.5));
      where.set_direction(d,origin.in_direction(d) +
                          ((index % nhere)+0.5)*(1.0/res));
      index /= nhere;
    }
  }
  return where;
}

int volume::ntot_at_resolution(double res) const {
  int mytot = 1;
  for (int d=X;d<=R;d++)
    if (has_boundary(High,(direction)d)) {
      const double dist = boundary_location(High,(direction)d)
                        - boundary_location(Low,(direction)d);
      mytot *= max(1,(int)(dist*res+0.5));
    }
  return mytot;
}

vec volume::loc(component c, int ind) const {
  return operator[](iloc(c,ind));
}

ivec volume::iloc(component c, int ind) const {
  ivec out(dim);
  LOOP_OVER_DIRECTIONS(dim,d) {
    int ind_over_stride = ind/stride(d);
    while (ind_over_stride < 0) ind_over_stride += num_direction(d)+1;
    out.set_direction(d, 2*(ind_over_stride%(num_direction(d)+1)));
  }
  return out + iyee_shift(c) + io;
}

vec volume::dr() const {
  switch (dim) {
  case Dcyl: return veccyl(inva, 0.0);
  case D1: case D2: case D3: abort("Error in dr\n");
  }
  return vec(0); // This is never reached.
}

vec volume::dx() const {
  switch (dim) {
  case D3: return vec(inva,0,0);
  case D2: return vec(inva,0);
  case D1: case Dcyl: abort("Error in dx.\n");
  }
  return vec(0); // This is never reached.
}

vec volume::dy() const {
  switch (dim) {
  case D3: return vec(0,inva,0);
  case D2: return vec(0,inva);
  case D1: case Dcyl: abort("Error in dy.\n");
  }
  return vec(0); // This is never reached.
}

vec volume::dz() const {
  switch (dim) {
  case Dcyl: return veccyl(0.0,inva);
  case D3: return vec(0,0,inva);
  case D1: return vec(inva);
  case D2: abort("dz doesn't exist in 2D\n");
  }
  return vec(0); // This is never reached.
}

volume volone(double zsize, double a) {
  return volume(D1, a, 0, 0, (int) (zsize*a + 0.5));
}

volume voltwo(double xsize, double ysize, double a) {
  return volume(D2, a, (xsize==0)?1:(int) (xsize*a + 0.5),
                       (ysize==0)?1:(int) (ysize*a + 0.5),0);
}

volume vol1d(double zsize, double a) {
  return volone(zsize, a);
}

volume vol2d(double xsize, double ysize, double a) {
  return voltwo(xsize, ysize, a);
}

volume vol3d(double xsize, double ysize, double zsize, double a) {
  return volume(D3, a,(xsize==0)?1:(int) (xsize*a + 0.5),
                      (ysize==0)?1:(int) (ysize*a + 0.5),
                      (zsize==0)?1:(int) (zsize*a + 0.5));
}

volume volcyl(double rsize, double zsize, double a) {
  if (zsize == 0.0) return volume(Dcyl, a, (int) (rsize*a + 0.5), 0, 1);
  else return volume(Dcyl, a, (int) (rsize*a + 0.5), 0, (int) (zsize*a + 0.5));
}

int volume::can_split_evenly(int n) const {
  int bestd = -1, bestlen = 0;
  for (int i=0;i<3;i++)
    if (num[i] > bestlen && num[i] % n == 0 && num[i] > 1) {
      bestd = i;
      bestlen = num[i];
    }
  if (bestd == -1) return 0;
  else return 1 + bestd;
}

// static int greatest_prime_factor(int n) {
//   for (int i=2;i<=n;i++) {
//     if (n % i == 0) {
//       while (n % i == 0) n /= i;
//       if (n == 1) return i;
//     }
//   }
//   abort("Can't calculate gpf of %d!\n", n);
// }

volume volume::split(int n, int which) const {
  if (n > nowned_min())
    abort("Cannot split %d grid points into %d parts\n", nowned_min(), n);
  if (n == 1) return *this;

  // Try to get as close as we can...
  int biglen = 0;
  for (int i=0;i<3;i++) if (num[i] > biglen) biglen = num[i];
  const int split_point = (int)(biglen*(n/2)/(double)n + 0.5);
  const int num_low = (int)(split_point*n/(double)biglen + 0.5);
  if (which < num_low)
    return split_at_fraction(false, split_point).split(num_low,which);
  else
    return split_at_fraction(true, split_point).split(n-num_low,which-num_low);
}

volume volume::split_by_effort(int n, int which, int Ngv, const volume *gv, double *effort) const {
  const int grid_points_owned = nowned_min();
  if (n > grid_points_owned)
    abort("Cannot split %d grid points into %d parts\n", nowned_min(), n);
  if (n == 1) return *this;
  int biglen = 0;
  direction splitdir;
  LOOP_OVER_DIRECTIONS(dim, d) if (num_direction(d) > biglen) { biglen = num_direction(d); splitdir = d; } 
  double best_split_measure = 1e20, left_effort_fraction;
  int best_split_point;
  vec corner = zero_vec(dim);
  LOOP_OVER_DIRECTIONS(dim, d) corner.set_direction(d, origin.in_direction(d) + num_direction(d)/a); 

  for (int split_point = 1; split_point < biglen; split_point+=1) {
    volume v_left = *this;
    v_left.set_num_direction(splitdir, split_point);
    volume v_right = *this;
    v_right.set_num_direction(splitdir, num_direction(splitdir) - split_point);
    v_right.origin_set_direction(splitdir, v_right.origin.in_direction(splitdir)+split_point/a);

    double total_left_effort = 0, total_right_effort = 0;
    volume vol;
    if (Ngv == 0) {
      total_left_effort = v_left.ntot();
      total_right_effort = v_right.ntot();
    }
    else {
      for (int j = 0; j<Ngv; j++) {
	if (v_left.intersect_with(gv[j], &vol))
	  total_left_effort += effort[j] * vol.ntot();
	if (v_right.intersect_with(gv[j], &vol))
	  total_right_effort += effort[j] * vol.ntot();
      }
    }
    double split_measure = max(total_left_effort/(n/2), total_right_effort/(n-n/2));
    if (split_measure < best_split_measure) {
      best_split_measure = split_measure;
      best_split_point = split_point;
      left_effort_fraction = total_left_effort/(total_left_effort + total_right_effort);
    }
  }
  const int split_point = best_split_point;
    
  const int num_low = (int)(left_effort_fraction *n + 0.5);
  // Revert to split() when effort method gives less grid points than chunks
  if (num_low > best_split_point*(grid_points_owned/biglen) || 
      (n-num_low) > (grid_points_owned - best_split_point*(grid_points_owned/biglen)))
    return split(n, which);

  if (which < num_low)
    return split_at_fraction(false, split_point).split_by_effort(num_low,which, Ngv,gv,effort);
  else
    return split_at_fraction(true, split_point).split_by_effort(n-num_low,which-num_low, Ngv,gv,effort);
}

volume volume::split_once(int n, int which) const {
  if (n == 1) return *this;
  int cse = can_split_evenly(n);
  if (cse) {
    const int bestd = cse-1;
    return split_specifically(n, which, (direction) bestd);
  } else {
    abort("Can't split when dimensions don't work out right\n");
    return volume(); // This is never reached.
  }
}

volume volume::split_at_fraction(bool want_high, int numer) const {
  int bestd = -1, bestlen = 1;
  for (int i=0;i<3;i++)
    if (num[i] > bestlen) {
      bestd = i;
      bestlen = num[i];
    }
  if (bestd == -1) {
    for (int i=0;i<3;i++) printf("num[%d] = %d\n", i, num[i]);
    abort("Crazy weird splitting error.\n");
  }
  volume retval(dim, a, 1,1,1);
  for (int i=0;i<3;i++) retval.num[i] = num[i];
  if (numer >= num[bestd])
    abort("Aaack bad bug in split_at_fraction.\n");
  direction d = (direction) bestd;
  if (dim == Dcyl && d == X) d = R;
  retval.set_origin(origin);
  if (want_high)
    retval.origin_set_direction(d,origin.in_direction(d)+numer/a);

  if (want_high) retval.num[bestd] -= numer;
  else retval.num[bestd] = numer;
  retval.num_changed();
  return retval;
}

volume volume::split_specifically(int n, int which, direction d) const {
  volume retval(dim, a, 1,1,1);
  for (int i=0;i<3;i++) retval.num[i] = num[i];

  vec shift = zero_vec(dim);
  shift.set_direction(d, num_direction(d)/n*which/a);
  retval.set_origin(origin + shift);

  retval.num[d % 3] /= n;
  retval.num_changed();
  return retval;
}

volume volume::pad(direction d) const {
  volume v(*this);
  v.pad_self(d);
  return v;
}

void volume::pad_self(direction d) {
  num[d%3]+=2; // Pad in both directions by one grid point.
  ivec temp = io;
  temp.set_direction(d, io.in_direction(d) - 2);
  num_changed();
  set_origin(operator[](temp));
}

ivec volume::icenter() const {
  // Find the center of the user's cell (which must be the symmetry
  // point):
  switch (dim) {
  case D1: return io + ivec(nz());
  case D2: return io + ivec(nx(), ny());
  case D3: return io + ivec(nx(), ny(), nz());
  case Dcyl: return io + iveccyl(0, nz());
  }
  abort("Can't do symmetry with these dimensions.\n");
  return ivec(0); // This is never reached.
}

vec volume::center() const {
  return operator[](icenter());
}

symmetry rotate4(direction axis, const volume &v) {
  symmetry s = identity();
  if (axis > 2) abort("Can only rotate4 in 2D or 3D.\n");
  s.g = 4;
  for (int d=0;d<5;d++) {
    s.S[d].d = (direction)d;
    s.S[d].flipped = false;
  }
  s.S[(axis+1)%3].d = (direction)((axis+2)%3);
  s.S[(axis+1)%3].flipped = true;
  s.S[(axis+2)%3].d = (direction)((axis+1)%3);
  s.symmetry_point = v.center();
  s.i_symmetry_point = v.icenter();
  return s;
}

symmetry rotate2(direction axis, const volume &v) {
  symmetry s = identity();
  if (axis > 2) abort("Can only rotate2 in 2D or 3D.\n");
  s.g = 2;
  s.S[(axis+1)%3].flipped = true;
  s.S[(axis+2)%3].flipped = true;
  s.symmetry_point = v.center();
  s.i_symmetry_point = v.icenter();
  return s;
}

symmetry mirror(direction axis, const volume &v) {
  symmetry s = identity();
  s.g = 2;
  s.S[axis].flipped = true;
  s.symmetry_point = v.center();
  s.i_symmetry_point = v.icenter();
  return s;
}

symmetry r_to_minus_r_symmetry(int m) {
  symmetry s = identity();
  s.g = 2;
  s.S[R].flipped = true;
  s.S[P].flipped = true;
  s.symmetry_point = zero_vec(Dcyl);
  s.i_symmetry_point = zero_ivec(Dcyl);
  s.ph = (m & 1) ? -1.0 : 1.0;
  return s;
}

symmetry identity() {
  return symmetry();
}

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
  if (s.next) next = new symmetry(*s.next);
  else next = NULL;
}

void symmetry::operator=(const symmetry &s) {
  g = s.g;
  for (int d=0;d<5;d++) {
    S[d].d = s.S[d].d;
    S[d].flipped = s.S[d].flipped;
  }
  ph = s.ph;
  symmetry_point = s.symmetry_point;
  i_symmetry_point = s.i_symmetry_point;
  if (s.next) next = new symmetry(*s.next);
  else next = NULL;
}

bool symmetry::operator==(const symmetry &sym) const {
  int gtot = multiplicity();
  if (gtot != sym.multiplicity())
    return false;
  for (int sn = 1; sn < gtot; ++sn)
    FOR_DIRECTIONS(d)
      if (transform(d, sn) != sym.transform(d, sn))
	return false;
  return true;
}

symmetry::~symmetry() {
  delete next;
}

int symmetry::multiplicity() const {
  if (next) return g*next->multiplicity();
  else return g;
}

symmetry symmetry::operator+(const symmetry &b) const {
  // The following optimization ignores identity when adding symmetries
  // together.  This is important because identity has an undefined
  // symmetry point.
  if (multiplicity() == 1) return b;
  else if (b.multiplicity() == 1) return *this;
  symmetry s = *this;
  s.next = new symmetry(b);
  return s;
}

symmetry symmetry::operator*(double p) const {
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
  } else {
       nme = n % g;
       nrest = n / g;
  }
  if (nme == 0) {
    if (nrest == 0) return signed_direction(d);
    else return next->transform(d,nrest);
  } else {
    signed_direction sd;
    if (nme == 1) sd = S[d];
    if (S[d].flipped) sd = flip(transform(S[d].d, nme-1));
    else sd = transform(S[d].d, nme-1);

    if (next && nrest) {
      if (sd.flipped) return flip(next->transform(sd.d, nrest))*ph;
      else return next->transform(sd.d, nrest)*ph;
    } else {
      return sd*ph;
    }
  }
}

ivec symmetry::transform(const ivec &ov, int n) const {
  if (n == 0) return ov;
  ivec out = ov;
  LOOP_OVER_DIRECTIONS(ov.dim, d) {
    const signed_direction s = transform(d,n);
    const int sp_d  = i_symmetry_point.in_direction(d);
    const int sp_sd = i_symmetry_point.in_direction(s.d);
    const int delta = ov.in_direction(d) - sp_d;
    if (s.flipped) out.set_direction(s.d, sp_sd - delta);
    else out.set_direction(s.d, sp_sd + delta);
  }
  return out;
}

ivec symmetry::transform_unshifted(const ivec &ov, int n) const {
  if (n == 0) return ov;
  ivec out(ov.dim);
  LOOP_OVER_DIRECTIONS(ov.dim, d) {
    const signed_direction s = transform(d,n);
    if (s.flipped) out.set_direction(s.d, -ov.in_direction(d));
    else out.set_direction(s.d, ov.in_direction(d));
  }
  return out;
}

vec symmetry::transform(const vec &ov, int n) const {
  if (n == 0) return ov;
  vec delta = ov;
  LOOP_OVER_DIRECTIONS(ov.dim, d) {
    const signed_direction s = transform(d,n);
    double deltad = ov.in_direction(d) - symmetry_point.in_direction(d);
    if (s.flipped) delta.set_direction(s.d, -deltad);
    else delta.set_direction(s.d, deltad);
  }
  return symmetry_point + delta;
}

geometric_volume symmetry::transform(const geometric_volume &gv, int n) const {
  return geometric_volume(transform(gv.get_min_corner(),n),
                          transform(gv.get_max_corner(),n));
}

component symmetry::transform(component c, int n) const {
  return direction_component(c,transform(component_direction(c),n).d);
}

derived_component symmetry::transform(derived_component c, int n) const {
  return direction_component(c,transform(component_direction(c),n).d);
}

int symmetry::transform(int c, int n) const {
  return (is_derived(c) ? int(transform(derived_component(c), n))
	  : int(transform(component(c), n)));
}

complex<double> symmetry::phase_shift(component c, int n) const {
  if (c == Dielectric) return 1.0;
  complex<double> phase = transform(component_direction(c),n).phase;
  // flip tells us if we need to flip the sign.  For vectors (E), it is
  // just this simple:
  bool flip = transform(component_direction(c),n).flipped;
  if (is_magnetic(c)) {
    // Because H is a pseudovector, here we have to figure out if the
    // transformation changes the handedness of the basis.
    bool have_one = false, have_two = false;
    FOR_DIRECTIONS(d) {
      if (transform(d,n).flipped) flip = !flip;
      int shift = (transform(d,n).d - d + 6) % 3;
      if (shift == 1) have_one = true;
      if (shift == 2) have_two = true;
    }
    if (have_one && have_two) flip = !flip;
  }
  if (flip) return -phase;
  else return phase;
}

complex<double> symmetry::phase_shift(derived_component c, int n) const {
  if (is_poynting(c)) {
    signed_direction ds = transform(component_direction(c),n);
    complex<double> ph = conj(ds.phase) * ds.phase; // E x H gets |phase|^2
    return (ds.flipped ? -ph : ph);
  }
  else /* energy density */
    return 1.0;
}

complex<double> symmetry::phase_shift(int c, int n) const {
  return (is_derived(c) ? phase_shift(derived_component(c), n)
	  : phase_shift(component(c), n));
}

bool symmetry::is_primitive(const ivec &p) const {
  // This is only correct if p is somewhere on the yee lattice.
  if (multiplicity() == 1) return true;
  for (int i=1;i<multiplicity();i++) {
    const ivec pp = transform(p,i);
    switch (p.dim) {
    case D2:
      if (pp.x()+pp.y() < p.x()+p.y()) return false;
      if (pp.x()+pp.y() == p.x()+p.y() &&
          p.y() > p.x() && pp.y() <= pp.x()) return false;
      break;
    case D3:
      if (pp.x()+pp.y()+pp.z() <  p.x()+p.y()+p.z()) return false;
      if (pp.x()+pp.y()+pp.z() == p.x()+p.y()+p.z() &&
          pp.x()+pp.y()-pp.z() <  p.x()+p.y()-p.z()) return false;
      if (pp.x()+pp.y()+pp.z() == p.x()+p.y()+p.z() &&
          pp.x()+pp.y()-pp.z() == p.x()+p.y()-p.z() &&
          pp.x()-pp.y()-pp.z() <  p.x()-p.y()-p.z()) return false;
      break;
    case D1: case Dcyl:
      if (pp.z() < p.z()) return false;
      break;
    }
  }
  return true;
}

/* given a list of geometric volumes, produce a new list with appropriate
   weights that is minimized according to the symmetry.  */
geometric_volume_list *symmetry::reduce(const geometric_volume_list *gl) const {
  geometric_volume_list *glnew = 0;
  for (const geometric_volume_list *g = gl; g; g = g->next) {
    int sn;
    for (sn = 0; sn < multiplicity(); ++sn) {
      geometric_volume gS(transform(g->gv, sn));
      int cS = transform(g->c, sn);
      geometric_volume_list *gn;
      for (gn = glnew; gn; gn = gn->next)
	if (gn->c == cS && gn->gv.round_float() == gS.round_float())
	  break;
      if (gn) { // found a match
	gn->weight += g->weight * phase_shift(g->c, sn);
	break;
      }
    }
    if (sn == multiplicity() && g->weight != 0.0) { // no match, add to glnew
      geometric_volume_list *gn = 
	new geometric_volume_list(g->gv, g->c, g->weight, glnew);
      glnew = gn;
    }
  }

  // reduce gv's redundant with themselves & delete elements with zero weight:
  geometric_volume_list *gprev = 0, *g = glnew;
  while (g) {
    // first, see if g->gv is redundant with itself
    bool halve[5] = {false,false,false,false,false};
    complex<double> weight = g->weight;
    for (int sn = 1; sn < multiplicity(); ++sn)
      if (g->c == transform(g->c, sn) && 
	  g->gv.round_float() == transform(g->gv, sn).round_float()) {
	LOOP_OVER_DIRECTIONS(g->gv.dim, d)
	  if (transform(d,sn).flipped) {
	    halve[d] = true;
	    break;
	  }
	g->weight += weight * phase_shift(g->c, sn);
      }
    LOOP_OVER_DIRECTIONS(g->gv.dim, d)
      if (halve[d])
	g->gv.set_direction_max(d, g->gv.in_direction_min(d) +
				0.5 * g->gv.in_direction(d));
    
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

} // namespace meep
