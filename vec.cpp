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

const char *dimension_name(ndim dim) {
  switch (dim) {
  case d1: return "1D";
  case dcyl: return "Cylindrical";
  }
  printf("Unsupported dimensionality name.\n");
  exit(1);
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
  printf("Unsupported case.\n");
  exit(1);
}

void vec::print(FILE *f) const {
  if (dim == dcyl) {
    fprintf(f, "%lg %lg 0", r(), z());
  } else if (dim == d1  ) {
    fprintf(f, "0 0 %lg", z());
  } else  {
    printf("I don't know how to print in this dimension!\n");
    fprintf(f, "I don't know how to print in this dimension!\n");
  }
}

volume::volume() {
  dim = d1;
  the_ntot = num[0] = num[1] = num[2] = 0;
  origin = vec(0);
}

inline int right_ntot(ndim d, const int num[3]) {
  switch (d) {
  case d1: return num[2] + 1;
  case dcyl: return (num[0]+1)*(num[2]+1);
  }
}

volume::volume(ndim d, double ta, int na, int nb, int nc) {
  dim = d;
  a = ta;
  inva = 1.0/a;
  num[0] = na;
  num[1] = nb;
  num[2] = nc;
  switch (d) {
  case dcyl: origin = vec(0,0);
    break;
  case d1: origin = vec(0);
    break;
    printf("Can't make volume with these dimensions!\n");
    exit(1);
  }
  the_ntot = right_ntot(d, num);
}

component volume::eps_component() const {
  switch (dim) {
  case d1: return Ex;
  case dcyl: return Hp;
  }
  printf("Unsupported dimensionality eps.\n");
  exit(1);
}

vec volume::yee_shift(component c) const {
  if (dim == dcyl) {
    vec offset(0,0);
    switch (c) {
    case Er:
    case Hz: offset = dr()*0.5; break;
    case Ez:
    case Hr: offset = dz()*0.5; break;
    case Hp: offset = (dz()+dr())*0.5; break;
    }
    return offset;
  } else if (dim == d1) {
    switch (c) {
    case Ex: return vec(0);
    case Hy: return dz()*0.5;
    }
  } else {
    printf("Invalid component!\n");
    exit(1);
  }
  printf("Unsupported dimension! yee shift of %s in %s\n",
         component_name(c), dimension_name(dim));
  exit(1);
}

int volume::contains(const vec &p) const {
  // containts returns true if the volume has any information in it
  // relevant to the point p.  Basically has is like owns (see below)
  // except it is more lenient, in that more than one lattice may contain a
  // given point.
  const double inva = 1.0/a;
  const double hinva = 0.5*inva;
  const vec o = p - origin;
  if (dim == dcyl) {
    return o.r() >= -hinva && o.z() >= -hinva &&
      o.r() <= nr()*inva - hinva && o.z() <= nz()*inva - hinva;
  } else if (dim == d1) {
    return o.z() >= -hinva && o.z() <= nz()*inva - hinva;
  } else {
    printf("Unsupported dimension.\n");
    exit(1);
  }
}

int volume::owns(const vec &p) const {
  // owns returns true if the point is closer to an active point on the
  // lattice than it is to any active point on a different lattice.  Thus
  // the "owns" borders of two adjacent volumes touch, while the "contains"
  // borders overlap.  "owns" is meant to indicate that when the point is a
  // lattice point, only one chunk actually *owns* that active lattice
  // point.
  const double inva = 1.0/a;
  const double qinva = 0.25*inva;
  const vec o = p - origin;
  if (dim == dcyl) {
    if (origin.r() == 0.0 && o.z() >= qinva && o.z() <= nz()*inva + qinva &&
        o.r() < qinva) return true;
    return o.r() >= qinva && o.z() >= qinva &&
      o.r() <= nr()*inva + qinva && o.z() <= nz()*inva + qinva;
  } else if (dim == d1) {
    return o.z() >= qinva && o.z() <= nz()*inva + qinva;
  } else {
    printf("Unsupported dimension.\n");
    exit(1);
  }
}

int volume::has_field(component c) const {
  const int is_r = c == Hr || c == Er;
  const int is_p = c == Hp || c == Ep;
  const int is_z = c == Hz || c == Ez;
  const int is_x = c == Hx || c == Ex;
  const int is_y = c == Hy || c == Ey;
  if (dim == dcyl) {
    return is_r || is_p || is_z;
  } else if (dim == d1) {
    return c == Ex || c == Hy;
  } else {
    printf("Aaack unsupported dim! (%d)\n", (int) dim);
    exit(1);
  }
}

int volume::index(component c, const vec &p) const {
  const vec offset = p - origin - yee_shift(c);
  int theindex = -1;
  if (dim == dcyl) {
    theindex = (int)((offset.z() + offset.r()*(nz()+1))*a + 0.5);
  } else if (dim == d1) {
    theindex = (int)(offset.z()*a + 0.5);
  } else {
    printf("Unsupported dimension.\n");
    exit(1);
  }
  return theindex;
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
  const vec offset = p - origin - yee_shift(c);
  if (dim == d1) {
    indices[0] = (int)(offset.z()*a);
    if (indices[0] < nz()) {
      indices[1] = indices[0] + 1;
      const double lowfrac = offset.z()*a - indices[0];
      if (lowfrac < 1e-15) {
        indices[0] = indices[1];
        weights[0] = 1.0;
        weights[1] = 0.0;
      } else {
        weights[0] = lowfrac;
        weights[1] = 1.0 - lowfrac;
      }
    } else {
      weights[0] = 1.0;
      weights[1] = 0.0;
      indices[1] = 0;
    }
    for (int i=2;i<8;i++) {
      indices[i] = 0;
      weights[i] = 0;
    }
  } else if (dim == dcyl) {
    interpolate_cyl(c, p, 0, indices, weights);
  } else {
    // by default use the 'index' routine.
    indices[0] = index(c,p);
    weights[0] = 1.0;
    for (int i=1;i<8;i++) {
      indices[i] = 0;
      weights[i] = 0;
    }
  }
  // Throw out out of range indices:
  for (int i=0;i<8;i++)
    if (indices[0] < 0 || indices[0] >= ntot()) weights[i] = 0.0;
  // Stupid very crude code to compactify arrays:
  stupidsort(indices, weights, 8);
}

void volume::interpolate_cyl(component c, const vec &p, int m,
                             int indices[8], double weights[8]) const {
  const vec indexpt = p - origin - yee_shift(c);
  const double ir = indexpt.r()*a, iz = indexpt.z()*a;
  int irlo = (int) floor(ir), izlo = (int) floor(iz),
    irhi = irlo+1, izhi = izlo+1;
  double dr = ir - (irlo + 0.5), dz = iz - (izlo + 0.5);
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
      weights[0] = (1+2*dz)*0.25;
      weights[2] = (1+2*dz)*0.25;
      weights[1] = -2*dz*(0.5-dr) + (1+2*dz)*0.25;
      weights[3] = -2*dz*(0.5+dr) + (1+2*dz)*0.25;
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
  const int odd_field_at_origin = m == 1 && have_r_zero && (c == Er || c == Hz);
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
    if (irhi > nr() || irhi < 0 || (irhi == 0 && !have_r_zero)) {
      weights[2] = 0;
      weights[3] = 0;
    } else if (irhi == -1) {
      indices[2] = izlo + (1+nz())*0; // NW
      indices[3] = izhi + (1+nz())*0; // NE
      if (c != Hp) { // These guys are odd!
        weights[2] = -weights[2];
        weights[3] = -weights[3];
      }
    }
    if (irlo > nr() || irlo < 0 || (irlo == 0 && !have_r_zero)) {
      weights[0] = 0;
      weights[1] = 0;
    } else if (irlo == -1) {
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
}

double volume::dv(component c, int ind) const {
  const double pi = 3.141592653589793238462643383276L;
  if (!owns(loc(c, ind))) return 0.0;
  switch (dim) {
  case dcyl: {
    const double r = loc(c,ind).r();
    const double Dr = inva*0.5;
    if (r != 0.0) return inva*pi*((r+Dr)*(r+Dr) - (r-Dr)*(r-Dr));
    else return inva*pi*(r+Dr)*(r+Dr);
  }
  case d1: return inva;
  }
}

vec volume::loc(component c, int ind) const {
  const vec offset = origin + yee_shift(c);
  switch (dim) {
  case dcyl: return offset + vec(inva*(ind/(nz()+1)), inva*(ind%(nz()+1)));
  case d1: return offset + vec(inva*ind);
  }
}

vec volume::dr() const {
  switch (dim) {
  case dcyl: return vec(inva, 0.0);
  }
}

vec volume::dz() const {
  switch (dim) {
  case dcyl: return vec(0.0,inva);
  case d1: return vec(inva);
  }
}

volume volone(double zsize, double a) {
  return volume(d1, a, 0, 0, (int) (zsize*a + 0.5));
}

volume volcyl(double rsize, double zsize, double a) {
  if (zsize == 0.0) return volume(dcyl, a, (int) (rsize*a + 0.5), 0, 1);
  else return volume(dcyl, a, (int) (rsize*a + 0.5), 0, (int) (zsize*a + 0.5));
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
  printf("Can't calculate gpf of %d!\n", n);
  exit(1);
}

volume volume::split(int n, int which) const {
  if (n == 1) return *this;
  int spl = greatest_prime_factor(n);
  if (spl == n) {
    return split_once(n, which);
  } else {
    volume v = split_once(spl, which % spl);
    return v.split(n/spl, which/spl);
  }
}

volume volume::split_once(int n, int which) const {
  if (n == 1) return *this;
  int cse = can_split_evenly(n);
  if (cse) {
    int bestd = cse-1;
    volume retval(dim, a,1);
    for (int i=0;i<3;i++) retval.num[i] = num[i];
    switch (dim) {
    case d1:
      retval.origin = origin + vec(nz()/n*which/a);
      break;
    case dcyl:
      if (bestd == 0) { // split on r
        retval.origin = origin + vec(nr()/n*which/a,0.0);
      } else { // split on z
        retval.origin = origin + vec(0.0,nz()/n*which/a);
      }
      break;
      printf("Can't split in this dimensionality!\n");
      exit(1);
    }
    retval.num[bestd] /= n;
    retval.the_ntot = right_ntot(dim, retval.num);
    return retval;
  } else {
    printf("Can't split when dimensions don't work out right\n");
    exit(1);
  }
}
