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

#include "vec.h"
#include "dactyl.h"

static inline int min(int a, int b) { return (a<b)?a:b; };
static inline int max(int a, int b) { return (a>b)?a:b; };
static inline double min(double a, double b) { return (a<b)?a:b; };
static inline double max(double a, double b) { return (a>b)?a:b; };

static inline double int_to_lattice(int n, double a, double inva=0.0) {
  if (inva == 0.0) inva = 1.0/a;
  return (2*n)*(0.5*inva);
}

static inline int lattice_to_int(double x, double a, double inva=0.0) {
  if (inva == 0.0) inva = 1.0/a;
  return (int)(x*(2.0*a) + 0.5)/2;
}

static inline double yee_to_lattice(int n, double a, double inva=0.0) {
  if (inva == 0.0) inva = 1.0/a;
  return n*(0.5*inva);
}

static inline int lattice_to_yee(double x, double a, double inva=0.0) {
  return (int)(x*(2.0*a) + 0.5);
}

inline ivec volume::round_vec(const vec &p) const {
  ivec result(dim);
  LOOP_OVER_DIRECTIONS(dim, d)
    result.set_direction(d, lattice_to_yee(p.in_direction(d),a,inva));
  return result;
}

ivec volume::io() const {
  return round_vec(origin);
}

static inline double rtl(double x, double a, double inva=0.0) {
  // Rounds to a value somewhere on the yee lattice.
  return ((int)(x*(2.0*a) + 0.5))*(0.5*inva);
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
}

const char *direction_name(direction d) {
  switch (d) {
  case X: return "x";
  case Y: return "y";
  case Z: return "z";
  case R: return "r";
  case P: return "phi";
  }
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
  }
  abort("Unsupported case.\n");
}

#include "mympi.h"

void vec::print(file *f) const {
  if (dim == Dcyl) {
    i_fprintf(f, "%lg %lg 0", r(), z());
  } else if (dim == D3  ) {
    i_fprintf(f, "%lg %lg %lg", x(), y(), z());
  } else if (dim == D2  ) {
    i_fprintf(f, "%lg %lg 0", x(), y());
  } else if (dim == D1  ) {
    i_fprintf(f, "0 0 %lg", z());
  } else  {
    printf("I don't know how to print in this dimension!\n");
    i_fprintf(f, "I don't know how to print in this dimension!\n");
  }
}

geometric_volume::geometric_volume(const vec &vec1, const vec &vec2) {
  dim = vec1.dim; 
  LOOP_OVER_DIRECTIONS(dim, d) {
    set_direction_min(d, min(vec1.in_direction(d), vec2.in_direction(d)));
    set_direction_max(d, max(vec1.in_direction(d), vec2.in_direction(d)));
  }
};
double geometric_volume::computational_volume() {
  double vol = 1.0; 
  LOOP_OVER_DIRECTIONS(dim,d) vol *= (in_direction_max(d) - in_direction_min(d));
  return vol;
};

double geometric_volume::full_volume() {
  double vol = 1.0; 
  LOOP_OVER_DIRECTIONS(dim, d) vol *= (in_direction_max(d) - in_direction_min(d));
  if (dim == Dcyl) vol *= pi * (in_direction_max(R) + in_direction_min(R));
  return vol;
};

geometric_volume geometric_volume::intersect_with(const geometric_volume &a) const {
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
};

int volume::yucky_num(int n) const {
  if (has_direction(dim, yucky_direction(n)))
    return num_direction(yucky_direction(n));
  return 1;
}

direction volume::yucky_direction(int n) const {
  if (dim == Dcyl)
    switch (n) {
    case 0: return X;
    case 1: return R;
    case 2: return Z;
    }
  return (direction) n;
}

geometric_volume volume::surroundings() const {
  geometric_volume res(dim);
  LOOP_OVER_DIRECTIONS(dim, d) {
    res.set_direction_min(d, operator[](little_corner()).in_direction(d));
    res.set_direction_max(d, operator[](big_corner()).in_direction(d));
  }
  return res;
}

inline int right_ntot(ndim di, const int num[3]) {
  int result = 1;
  LOOP_OVER_DIRECTIONS(di, d) result *= num[d%3]+1;
  return result;
}

volume::volume(ndim d, double ta, int na, int nb, int nc) {
  dim = d;
  a = ta;
  inva = 1.0/a;
  num[0] = na;
  num[1] = nb;
  num[2] = nc;
  origin = zero_vec(dim);
  the_ntot = right_ntot(d, num);
  set_strides();
}

component volume::eps_component() const {
  switch (dim) {
  case D1: return Ex;
  case D2: return Ez;
  case D3: return Dielectric;
  case Dcyl: return Hp;
  }
  abort("Unsupported dimensionality eps.\n");
}

vec volume::yee_shift(component c) const {
  return operator[](iyee_shift(c));
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
  const ivec o = p - io();
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
  const double inva = 1.0/a;
  const vec o = p - origin;
  LOOP_OVER_DIRECTIONS(dim, d)
    if (o.in_direction(d) < -inva || o.in_direction(d) > num_direction(d)*inva+inva)
      return false;
  return true;
}

bool volume::owns(const ivec &p) const {
  // owns returns true if the point "owned" by this volume, meaning that it
  // is the volume that would timestep the point.
  const ivec o = p - io();
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
  }
}

int volume::has_boundary(boundary_side b,direction d) const {
  switch (dim) {
  case Dcyl: return d == Z || (d == R && b == High);
  case D1: return d == Z;
  case D2: return d == X || d == Y;
  case D3: return d == X || d == Y || d == Z;
  }
}

int volume::index(component c, const ivec &p) const {
  const ivec offset = p - io() - iyee_shift(c);
  int idx = 0;
  LOOP_OVER_DIRECTIONS(dim,d) idx += offset.in_direction(d)/2*stride(d);
  return idx;
}

void volume::set_strides() {
  for (int d=0;d<5;d++) the_stride[d] = 1; // Yuck yuck yuck.
  LOOP_OVER_DIRECTIONS(dim,d)
    switch(d) {
    case Z: the_stride[d] = 1; break;
    case R: the_stride[d] = nz()+1; break;
    case X: the_stride[d] = (nz()+1)*(ny() + 1); break;
    case Y: the_stride[d] = nz() + 1; break;
    }
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

void volume::interpolate(component c, const vec &p,
                         int indices[8], double weights[8]) const {
  if (dim == D1) {
    interpolate_one(c, p, indices, weights);
  } else if (dim == D2) {
    interpolate_two(c, p, indices, weights);
  } else if (dim == D3) {
    interpolate_fancy(c, p, indices, weights);
  } else if (dim == Dcyl) {
    interpolate_cyl(c, p, 0, indices, weights);
  } else {
    abort("Can't interpolate in these dimensions.\n");
  }
  // Throw out out of range indices:
  for (int i=0;i<8;i++)
    if (indices[0] < 0 || indices[0] >= ntot()) weights[i] = 0.0;
  // Stupid very crude code to compactify arrays:
  stupidsort(indices, weights, 8);
  if (!contains(p) && weights[0])
    abort("Error made in interpolation of %s--fix this bug!!!\n",
          component_name(c));
}

void volume::interpolate_one(component c, const vec &p,
                             int indices[8], double weights[8]) const {
  const double iz = (p.z() - yee_shift(c).z())*a;
  const int ioriginz = (int) (origin.z()*a + 0.5);
  const int izlo = (int)iz - ioriginz, izhi = izlo+1;
  const double dz = p.z()*a - (ioriginz + izlo + 0.5) - yee_shift(c).z()*a;
  indices[0] = izlo;
  const int indexshift = (c == Ex)? 1 : 0;
  for (int i=0;i<8;i++) {
    indices[i] = -1;
    weights[i] = 0;
  }
  if (izlo - indexshift < nz() && izlo - indexshift >= 0) {
    indices[0] = izlo;
    weights[0] = 0.5 - dz;
  }
  if (izhi - indexshift < nz() && izhi - indexshift >= 0) {
    indices[1] = izhi;
    weights[1] = 0.5 + dz;
  }
  stupidsort(indices, weights, 2);
  if (!contains(p) && weights[0]) {
    printf("Looks like a bug in one D interpolate!\n");
    printf("Point   %lg %lg\n", p.r(), p.z());
    abort("Looks like a bug in one D interpolate!\n");
  }
}

void volume::interpolate_cyl(component c, const vec &p, int m,
                             int indices[8], double weights[8]) const {
  const double ir = (p.r() - yee_shift(c).r())*a;
  const double iz = (p.z() - yee_shift(c).z())*a;
  // (0,0) must be on grid!
  const int ioriginr = (int) (origin.r()*a + 0.5);
  const int ioriginz = (int) (origin.z()*a + 0.5);
  int irlo = (int)ir - ioriginr, izlo = (int)iz - ioriginz,
    irhi = irlo+1, izhi = izlo+1;
  const double dr = p.r()*a - (ioriginr + irlo + 0.5) - yee_shift(c).r()*a;
  const double dz = p.z()*a - (ioriginz + izlo + 0.5) - yee_shift(c).z()*a;
  for (int i=0;i<8;i++) indices[i] = 0;
  for (int i=0;i<8;i++) weights[i] = 0;
  // Tabulate the actual weights:
  if (dr+dz > 0.0) {
    if (dr-dz > 0.0) { // North
      weights[0] = (1-2*dr)*0.25;
      weights[1] = (1-2*dr)*0.25;
      weights[2] = 2*dr*(0.5-dz) + (1-2*dr)*0.25;
      weights[3] = 2*dr*(0.5+dz) + (1-2*dr)*0.25;
    } else { // East
      weights[0] = (1-2*dz)*0.25;
      weights[2] = (1-2*dz)*0.25;
      weights[1] = 2*dz*(0.5-dr) + (1-2*dz)*0.25;
      weights[3] = 2*dz*(0.5+dr) + (1-2*dz)*0.25;
    }
  } else {
    if (dr-dz > 0.0) { // West
      weights[3] = (1+2*dz)*0.25;
      weights[1] = (1+2*dz)*0.25;
      weights[0] = -2*dz*(0.5-dr) + (1+2*dz)*0.25;
      weights[2] = -2*dz*(0.5+dr) + (1+2*dz)*0.25;
    } else { // South
      weights[2] = (1+2*dr)*0.25;
      weights[3] = (1+2*dr)*0.25;
      weights[0] = -2*dr*(0.5-dz) + (1+2*dr)*0.25;
      weights[1] = -2*dr*(0.5+dz) + (1+2*dr)*0.25;
    }
  }
  // These are the four nearest neighbor points:
  indices[0] = izlo + (1+nz())*irlo; // SW
  indices[1] = izhi + (1+nz())*irlo; // SE
  indices[2] = izlo + (1+nz())*irhi; // NW
  indices[3] = izhi + (1+nz())*irhi; // NE
  // Figure out which of these points is off the grid in z direction:
  switch ((component)c) {
  case Ep: case Er: case Hz:
    if (izlo <= 0 || izlo > nz()) {
      weights[0] = 0;
      weights[2] = 0;
    }
    if (izhi <= 0 || izhi > nz()) {
      weights[1] = 0;
      weights[3] = 0;
    }
    break;
  case Hp: case Hr: case Ez:
    if (izlo < 0 || izlo >= nz()) {
      weights[0] = 0;
      weights[2] = 0;
    }
    if (izhi < 0 || izhi >= nz()) {
      weights[1] = 0;
      weights[3] = 0;
    }
  }
  // Figure out which of the points is off the grid in r direction:
  const int have_r_zero = origin.r() == 0.0;
  const int odd_field_at_origin = (m == 1 && have_r_zero && (c == Er || c == Hz))
                               || (m == 0 && have_r_zero && (c == Hp));
  // Ep is also odd, but we don't need to interpolate it to zero.
  switch ((component)c) {
  case Ep: case Ez: case Hr:
    if (irhi > nr() || irhi < 0 || (irhi == 0 && !have_r_zero)) {
      weights[2] = 0;
      weights[3] = 0;
    }
    if (irlo > nr() || irlo < 0 || (irlo == 0 && !have_r_zero)) {
      weights[0] = 0;
      weights[1] = 0;
    }
    break;
  case Hp: case Hz: case Er:
    if (irhi >= nr() || irhi < 0) {
      weights[2] = 0;
      weights[3] = 0;
    } else if (irhi == -1 && have_r_zero) {
      indices[2] = izlo + (1+nz())*0; // NW
      indices[3] = izhi + (1+nz())*0; // NE
      if (c != Hp) { // These guys are odd!
        weights[2] = -weights[2];
        weights[3] = -weights[3];
      }
    }
    if (irlo >= nr() || irlo < 0) {
      weights[0] = 0;
      weights[1] = 0;
    } else if (irlo == -1 && have_r_zero) {
      if (c != Hp) { // These guys are odd!
        weights[2] -= weights[0];
        weights[3] -= weights[1];
      } else {
        weights[2] += weights[0];
        weights[3] += weights[1];
      }
      weights[0] = 0;
      weights[1] = 0;
    }
  }
  // Now I need to reorder them so the points with nonzero weights come
  // first.
  stupidsort(indices, weights, 4);
  for (int i=0;i<8&&weights[i];i++)
    if (!owns(iloc(c, indices[i])))
      weights[i] = 0.0;
  stupidsort(indices, weights, 4);
  if (!contains(p) && weights[0]) {
    printf("Error made in cyl interpolation--fix this bug!!!\n");
    printf("%s irlo %d irhi %d izlo %d izhi %d\n",
           component_name(c), irlo, irhi, izlo, izhi);
    printf("Point is at %lg %lg -- in real space this is %lg %lg\n",
           ir, iz, p.r(), p.z());
    printf("  dr %.20lg dz %.20lg\n", dr, dz);
    for (int i=0;i<8 &&weights[i];i++) {
      printf("  Point %lg %lg Weight %.25lg\n",
             loc(c, indices[i]).r(), loc(c, indices[i]).z(), weights[i]);
      if (!owns(iloc(c, indices[i]))) {
        printf("  ...we don't own this index!\n");
        weights[i] = 0.0;
      }
    }
    abort("aack");
  }
}

void volume::interpolate_two(component c, const vec &p,
                             int indices[8], double weights[8]) const {
  const double ix = (p.x() - yee_shift(c).x())*a;
  const double iy = (p.y() - yee_shift(c).y())*a;
  // (0,0) must be on grid!
  const int ioriginx = (int) (origin.x()*a + 0.5);
  const int ioriginy = (int) (origin.y()*a + 0.5);
  int ixlo = (int)ix - ioriginx, iylo = (int)iy - ioriginy,
    ixhi = ixlo+1, iyhi = iylo+1;
  const double dx = p.x()*a - (ioriginx + ixlo + 0.5) - yee_shift(c).x()*a;
  const double dy = p.y()*a - (ioriginy + iylo + 0.5) - yee_shift(c).y()*a;
  for (int i=0;i<8;i++) indices[i] = 0;
  for (int i=0;i<8;i++) weights[i] = 0;
  // Tabulate the actual weights:
  if (dx+dy > 0.0) {
    if (dx-dy > 0.0) { // North
      weights[0] = (1-2*dx)*0.25;
      weights[1] = (1-2*dx)*0.25;
      weights[2] = 2*dx*(0.5-dy) + (1-2*dx)*0.25;
      weights[3] = 2*dx*(0.5+dy) + (1-2*dx)*0.25;
    } else { // East
      weights[0] = (1-2*dy)*0.25;
      weights[2] = (1-2*dy)*0.25;
      weights[1] = 2*dy*(0.5-dx) + (1-2*dy)*0.25;
      weights[3] = 2*dy*(0.5+dx) + (1-2*dy)*0.25;
    }
  } else {
    if (dx-dy > 0.0) { // West
      weights[3] = (1+2*dy)*0.25;
      weights[1] = (1+2*dy)*0.25;
      weights[0] = -2*dy*(0.5-dx) + (1+2*dy)*0.25;
      weights[2] = -2*dy*(0.5+dx) + (1+2*dy)*0.25;
    } else { // South
      weights[2] = (1+2*dx)*0.25;
      weights[3] = (1+2*dx)*0.25;
      weights[0] = -2*dx*(0.5-dy) + (1+2*dx)*0.25;
      weights[1] = -2*dx*(0.5+dy) + (1+2*dx)*0.25;
    }
  }
  // These are the four nearest neighbor points:
  indices[0] = iylo + (1+ny())*ixlo; // SW
  indices[1] = iyhi + (1+ny())*ixlo; // SE
  indices[2] = iylo + (1+ny())*ixhi; // NW
  indices[3] = iyhi + (1+ny())*ixhi; // NE
  // Figure out which of these points is off the grid in y direction:
  switch ((component)c) {
  case Ez: case Ex: case Hy:
    if (iylo <= 0 || iylo > ny()) {
      weights[0] = 0;
      weights[2] = 0;
    }
    if (iyhi <= 0 || iyhi > ny()) {
      weights[1] = 0;
      weights[3] = 0;
    }
    break;
  case Hz: case Hx: case Ey:
    if (iylo < 0 || iylo >= ny()) {
      weights[0] = 0;
      weights[2] = 0;
    }
    if (iyhi < 0 || iyhi >= ny()) {
      weights[1] = 0;
      weights[3] = 0;
    }
  }
  // Figure out which of the points is off the grid in x direction:
  switch ((component)c) {
  case Ez: case Ey: case Hx:
    if (ixhi > nx() || ixhi < 0) {
      weights[2] = 0;
      weights[3] = 0;
    }
    if (ixlo > nx() || ixlo < 0) {
      weights[0] = 0;
      weights[1] = 0;
    }
    break;
  case Hz: case Hy: case Ex:
    if (ixhi >= nx() || ixhi < 0) {
      weights[2] = 0;
      weights[3] = 0;
    }
    if (ixlo >= nx() || ixlo < 0) {
      weights[0] = 0;
      weights[1] = 0;
    }
  }
  // Now I need to reorder them so the points with nonzero weights come
  // first.
  stupidsort(indices, weights, 4);
  for (int i=0;i<8&&weights[i];i++)
    if (!owns(iloc(c, indices[i])))
      weights[i] = 0.0;
  stupidsort(indices, weights, 4);
  for (int i=0;i<8&&weights[i];i++)
    if (!owns(iloc(c, indices[i])))
      abort("Aaack, I don't actually own this! (2D)\n");
  if (!contains(p) && weights[0]) {
    abort("Error made in 2D interpolation--fix this bug!!!\n");
  }
}

void volume::interpolate_fancy(component c, const vec &pc,
                               int indices[8], double weights[8]) const {
  const vec p = (pc - yee_shift(c))*a;
  ivec middle(dim);
  LOOP_OVER_DIRECTIONS(dim,d)
    middle.set_direction(d, (int) p.in_direction(d)*2+1);
  middle += iyee_shift(c);
  const vec midv = this->operator[](middle);
  const vec dv = (pc - midv)*(2*a);
  int already_have = 1;
  ivec locs[8];
  for (int i=0;i<8;i++) {
    locs[i] = round_vec(midv);
    weights[i] = 1.0;
  }
  LOOP_OVER_DIRECTIONS(dim,d) {
    for (int i=0;i<already_have;i++) {
      locs[already_have+i] = locs[i];
      weights[already_have+i] = weights[i];
      locs[i].set_direction(d,middle.in_direction(d)-1);
      if (dv.in_direction(d) < 0.0) weights[i] *= -dv.in_direction(d);
      else weights[i] = 0.0;
      locs[already_have+i].set_direction(d,middle.in_direction(d)+1);
      if (dv.in_direction(d) > 0.0) weights[already_have+i] *= dv.in_direction(d);
      else weights[already_have+i] = 0.0;
    }
    already_have *= 2;
  }
  for (int i=already_have;i<8;i++) weights[i] = 0.0;
  double total_weight = 0.0;
  for (int i=0;i<already_have;i++) total_weight += weights[i];
  for (int i=0;i<already_have;i++)
    weights[i] += (1.0 - total_weight)*(1.0/already_have);
  for (int i=0;i<already_have;i++) {
    if (!owns(locs[i])) weights[i] = 0.0;
    else indices[i] = index(c,locs[i]);
  }
  stupidsort(indices, weights, already_have);
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
  const double inva = 1.0/a, qinva = 0.25*inva;
  return origin.x() + nx()*inva + qinva;
}

double volume::xmin() const {
  const double inva = 1.0/a, qinva = 0.25*inva;
  return origin.x() + qinva;
}

double volume::ymax() const {
  const double inva = 1.0/a, qinva = 0.25*inva;
  return origin.y() + ny()*inva + qinva;
}

double volume::ymin() const {
  const double inva = 1.0/a, qinva = 0.25*inva;
  return origin.y() + qinva;
}

double volume::zmax() const {
  const double inva = 1.0/a, qinva = 0.25*inva;
  return origin.z() + nz()*inva + qinva;
}

double volume::zmin() const {
  const double inva = 1.0/a, qinva = 0.25*inva;
  return origin.z() + qinva;
}

double volume::rmax() const {
  const double inva = 1.0/a, qinva = 0.25*inva;
  if (dim == Dcyl) return origin.r() + nr()*inva + qinva;
  abort("No rmax in these dimensions.\n");
}

double volume::rmin() const {
  const double inva = 1.0/a, qinva = 0.25*inva;
  if (dim == Dcyl) {
    if (origin.r() == 0.0) {
      return 0.0;
    } else {
      return origin.r() + qinva;
    }
  }
  abort("No rmin in these dimensions.\n");
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
  }
  else switch (d) {
  case X: return loc(Ez,0).x();
  case Y: return loc(Ez,0).y();
  case R: return loc(Ep,0).r();
  case Z: if (dim == Dcyl) return loc(Ep,0).z();
          else return loc(Ex,0).z();
  }
}

ivec volume::big_corner() const {
  switch (dim) {
  case D1: return io() + ivec(nz())*2;
  case D2: return io() + ivec2d(nx(),ny())*2;
  case D3: return io() + ivec(nx(),ny(),nz())*2;
  case Dcyl: return io() + ivec(nr(),nz())*2;
  }
}

vec volume::loc_at_resolution(int index, double res) const {
  vec where = origin;
  for (int dd=X;dd<=R;dd++) {
    const direction d = (direction) dd;
    if (has_boundary(High,d)) {
      const double dist = boundary_location(High,d)-boundary_location(Low,d);
      const int nhere = max(1,(int)(dist*res+0.5));
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
  LOOP_OVER_DIRECTIONS(dim,d)
    out.set_direction(d, 2*((ind/stride(d))%(num_direction(d)+1)));
  return out + iyee_shift(c) + io();
}

vec volume::dr() const {
  switch (dim) {
  case Dcyl: return vec(inva, 0.0);
  }
}

vec volume::dx() const {
  switch (dim) {
  case D3: return vec(inva,0,0);
  case D2: return vec2d(inva,0);
  }
}

vec volume::dy() const {
  switch (dim) {
  case D3: return vec(0,inva,0);
  case D2: return vec2d(0,inva);
  }
}

vec volume::dz() const {
  switch (dim) {
  case Dcyl: return vec(0.0,inva);
  case D3: return vec(0,0,inva);
  case D1: return vec(inva);
  }
}

volume volone(double zsize, double a) {
  return volume(D1, a, 0, 0, (int) (zsize*a + 0.5));
}

volume voltwo(double xsize, double ysize, double a) {
  return volume(D2, a, (int) (xsize*a + 0.5), (int) (ysize*a + 0.5),0);
}

volume vol3d(double xsize, double ysize, double zsize, double a) {
  return volume(D3, a, (int) (xsize*a + 0.5), (int) (ysize*a + 0.5),
                (int) (zsize*a + 0.5));
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

static int greatest_prime_factor(int n) {
  for (int i=2;i<=n;i++) {
    if (n % i == 0) {
      while (n % i == 0) n /= i;
      if (n == 1) return i;
    }
  }
  abort("Can't calculate gpf of %d!\n", n);
}

volume volume::split(int n, int which) const {
  if (n == 1) return *this;
  int spl = greatest_prime_factor(n);
  if (false && can_split_evenly(spl)) { //disable prime factor
    // Deal with case where we can split evenly by this big prime factor.
    if (spl == n) {
      return split_once(n, which);
    } else {
      volume v = split_once(spl, which % spl);
      return v.split(n/spl, which/spl);
    }
  } else {
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
}

volume volume::split_once(int n, int which) const {
  if (n == 1) return *this;
  int cse = can_split_evenly(n);
  if (cse) {
    const int bestd = cse-1;
    return split_specifically(n, which, (direction) bestd);
  } else {
    abort("Can't split when dimensions don't work out right\n");
  }
}

volume volume::split_at_fraction(bool want_high, int numer) const {
  int bestd = -1, bestlen = 1;
  for (int i=0;i<3;i++)
    if (num[i] > bestlen) {
      bestd = i;
      bestlen = num[i];
    }
  if (bestd == -1) abort("Crazy weird splitting error.\n");
  volume retval(dim, a,1);
  for (int i=0;i<3;i++) retval.num[i] = num[i];
  if (numer >= num[bestd])
    abort("Aaack bad bug in split_at_fraction.\n");
  direction d = (direction) bestd;
  if (dim == Dcyl && d == X) d = R;
  retval.origin = origin;
  if (want_high)
    retval.origin.set_direction(d,origin.in_direction(d)+numer/a);

  if (want_high) retval.num[bestd] -= numer;
  else retval.num[bestd] = numer;
  retval.the_ntot = right_ntot(dim, retval.num);
  retval.set_strides();
  return retval;
}

volume volume::split_specifically(int n, int which, direction d) const {
  volume retval(dim, a,1);
  for (int i=0;i<3;i++) retval.num[i] = num[i];

  vec shift = zero_vec(dim);
  shift.set_direction(d, num_direction(d)/n*which/a);
  retval.origin = origin + shift;

  retval.num[d % 3] /= n;
  retval.the_ntot = right_ntot(dim, retval.num);
  retval.set_strides();
  return retval;
}

volume volume::pad(direction d) const {
  volume v = *this;
  v.num[d%3]++;
  v.the_ntot = right_ntot(dim, v.num);
  v.set_strides();
  return v;
}

vec volume::center() const {
  // Find the center of the user's cell (which must be the symmetry
  // point):
  const double inva = 1.0/a;
  vec almost_center;
  switch (dim) {
  case D2:
    almost_center = origin + vec2d(nx()/2*inva, ny()/2*inva);
    return operator[](round_vec(almost_center));
  case D3:
    almost_center = origin + vec(nx()/2*inva, ny()/2*inva, nz()/2*inva);
  }
  abort("Can't do symmetry with these dimensions.\n");
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
  s.a = v.a;
  s.inva = 1.0/v.a;
  return s;
}

symmetry rotate2(direction axis, const volume &v) {
  symmetry s = identity();
  if (axis > 2) abort("Can only rotate2 in 2D or 3D.\n");
  s.g = 2;
  s.S[(axis+1)%3].flipped = true;
  s.S[(axis+2)%3].flipped = true;
  s.symmetry_point = v.center();
  s.a = v.a;
  s.inva = 1.0/v.a;
  return s;
}

symmetry mirror(direction axis, const volume &v) {
  symmetry s = identity();
  s.g = 2;
  s.S[axis].flipped = true;
  s.symmetry_point = v.center();
  s.a = v.a;
  s.inva = 1.0/v.a;
  return s;
}

symmetry identity() {
  return symmetry();
}

symmetry::symmetry() {
  g = 1;
  for (int d=0;d<5;d++) {
    S[d].d = (direction)d;
    S[d].flipped = false;
  }
  // Set these to NaN just to make sure I don't use them.
  a = inva = sqrt(-1.0);
  next = NULL;
}

symmetry::symmetry(const symmetry &s) {
  g = s.g;
  for (int d=0;d<5;d++) {
    S[d].d = s.S[d].d;
    S[d].flipped = s.S[d].flipped;
  }
  symmetry_point = s.symmetry_point;
  if (s.next) next = new symmetry(*s.next);
  else next = NULL;
  a = s.a;
  inva = s.inva;
}

void symmetry::operator=(const symmetry &s) {
  g = s.g;
  for (int d=0;d<5;d++) {
    S[d].d = s.S[d].d;
    S[d].flipped = s.S[d].flipped;
  }
  symmetry_point = s.symmetry_point;
  if (s.next) next = new symmetry(*s.next);
  else next = NULL;
  a = s.a;
  inva = s.inva;
}

symmetry::~symmetry() {
  delete next;
}

int symmetry::multiplicity() const {
  if (next) return g*next->multiplicity();
  else return g;
}

symmetry symmetry::operator+(const symmetry &b) const {
  symmetry s = *this;
  s.next = new symmetry(b);
  return s;
}

signed_direction symmetry::transform(direction d, int n) const {
  // Returns direction or if opposite, 
  const int nme = n % g;
  const int nrest = n / g;
  if (nme == 0) {
    if (nrest == 0) return signed_direction(d);
    else return next->transform(d,nrest);
  } else {
    signed_direction sd;
    if (nme == 1) sd = S[d];
    if (S[d].flipped) sd = flip(transform(S[d].d, nme-1));
    else sd = transform(S[d].d, nme-1);

    if (next && nrest) {
      if (sd.flipped) return flip(next->transform(sd.d, nrest));
      else return next->transform(sd.d, nrest);
    } else {
      return sd;
    }
  }
}

ivec symmetry::transform(const ivec &ov, int n) const {
  if (n == 0) return ov;
  ivec out = ov;
  LOOP_OVER_DIRECTIONS(ov.dim, d) {
    const signed_direction s = transform(d,n);
    const int sp_d  = lattice_to_yee(symmetry_point.in_direction(d),a,inva);
    const int sp_sd = lattice_to_yee(symmetry_point.in_direction(s.d),a,inva);
    const int delta = ov.in_direction(d) - sp_d;
    if (s.flipped) out.set_direction(s.d, sp_sd - delta);
    else out.set_direction(s.d, sp_sd + delta);
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

component symmetry::transform(component c, int n) const {
  return direction_component(c,transform(component_direction(c),n).d);
}

complex<double> symmetry::phase_shift(component c, int n) const {
  if (is_magnetic(c)) {
    // Because H is a pseudovector, here we have to figure out if it is an
    // inversion... to do this we'll figure out what happens when we
    // transform the two other directions.
    const direction d0 = component_direction(c);
    const direction d1 = (direction) ((d0+1)%3);
    const direction d2 = (direction) ((d0+2)%3);
    // Note that according to the above definition, c1 and c2 are in
    // cyclical permutation order, i.e. (c2 - c1)%3 == 1.  If the
    // transformation of these two directions is no longer cyclical, then
    // we need to flip the sign of this component of H.
    bool flip = false;
    if (((3+transform(d2,n).d - transform(d1,n).d)%3) == 2) flip = !flip;
    if (transform(d1,n).flipped) flip = !flip;
    if (transform(d2,n).flipped) flip = !flip;
    if (flip) return -1.0;
    else return 1.0;
  } else {
    if (transform(component_direction(c),n).flipped) return -1.0;
    else return 1.0;
  }
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
      if (pp.x()+pp.y()+pp.z() < p.x()+p.y()+p.z()) return false;
      if (pp.x()+pp.y()+pp.z() == p.x()+p.y()+p.z() &&
          pp.y() > pp.x() && p.y() <= p.x()) return false;
      abort("there is a known bug in 3d is_primitive!\n");// FIXME!
      break;
    case D1: case Dcyl:
      if (pp.z() < p.z()) return false;
      break;
    }
  }
  return true;
}
