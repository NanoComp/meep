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

enum component { Ex=0, Ey, Er, Ep, Ez, Hx, Hy, Hr, Hp, Hz };
enum ndim { d1=0, d2, d3, dcyl };
enum field_type { E_stuff=0, H_stuff=1 };
enum boundary_side { High=0, Low };
enum direction { X=0,Y,Z,R,P };

inline int is_electric(component c) { return (int) c < 5; }
inline int is_magnetic(component c) { return (int) c >= 5; }
inline field_type type(component c) {
  if (is_electric(c)) return E_stuff;
  else return H_stuff;
}
const char *component_name(component c);
inline direction component_direction(component c) {
  switch (c) {
  case Ex: case Hx: return X;
  case Ey: case Hy: return Y;
  case Ez: case Hz: return Z;
  case Er: case Hr: return R;
  case Ep: case Hp: return P;
  }
}

class vec {
 public:
  vec() { dim = d2; t[X] = t[Y] = 0; };
  vec(double zz) { dim = d1; t[Z] = zz; };
  vec(double rr, double zz) { dim = dcyl; t[R] = rr; t[Z] = zz; };
  vec(double xx, double yy, double zz) {
    dim = d3; t[X] = xx; t[Y] = yy; t[Z] = zz; };
  friend vec vec2d(double xx, double yy);
  ~vec() {};

  vec operator+(const vec &a) const {
    switch (dim) {
    case dcyl: return vec(t[R]+a.t[R],t[Z]+a.t[Z]);
    case d3: return vec(t[X]+a.t[X],t[Y]+a.t[Y],t[Z]+a.t[Z]);
    case d2: return vec2d(t[X]+a.t[X],t[Y]+a.t[Y]);
    case d1: return vec(t[Z]+a.t[Z]);
    }
  };
  vec operator+=(const vec &a) {
    switch (dim) {
    case dcyl: t[R] += a.t[R]; t[Z] += a.t[Z]; return *this;
    case d3: t[X] += a.t[X]; t[Y] += a.t[Y]; t[Z] += a.t[Z]; return *this;
    case d2: t[X] += a.t[X]; t[Y] += a.t[Y]; return *this;
    case d1: t[Z] += a.t[Z]; return vec(t[Z]+a.t[Z]);
    }
  };
  vec operator-(const vec &a) const {
    switch (dim) {
    case dcyl: return vec(t[R]-a.t[R],t[Z]-a.t[Z]);
    case d3: return vec(t[X]-a.t[X],t[Y]-a.t[Y],t[Z]-a.t[Z]);
    case d2: return vec2d(t[X]-a.t[X],t[Y]-a.t[Y]);
    case d1: return vec(t[Z]-a.t[Z]);
    }
  };
  vec operator==(const vec &a) const;
  vec operator*(double s) const {
    switch (dim) {
    case dcyl: return vec(t[R]*s,t[Z]*s);
    case d3: return vec(t[X]*s,t[Y]*s,t[Z]*s);
    case d2: return vec2d(t[X]*s,t[Y]*s);
    case d1: return vec(t[Z]*s);
    }
  };
  ndim dim;

  double r() const { return t[R]; };
  double x() const { return t[X]; };
  double y() const { return t[Y]; };
  double z() const { return t[Z]; };
  double in_direction(direction d) const { return t[d]; };
  double set_direction(direction d, double val) { t[d] = val; };

  double project_to_boundary(direction, double boundary_loc);
  void print(FILE *) const;
  friend vec zero_vec(ndim);
 private:
  double t[5];
};

inline vec zero_vec(ndim di) {
  vec v; v.dim = di; for (int d=0;d<5;d++) v.set_direction((direction)d,0.0);
  return v;
}

inline vec vec2d(double xx, double yy) {
  vec v; v.t[X] = xx; v.t[Y] = yy; return v;
}

class plane {
 public:
  plane(const vec &center, const vec &half_diagonal) {
    c = center;
    d = half_diagonal;
  };

  vec center() const { return c; };
  vec diagonal() const { return d; };
 private:
  vec c, d;
};

class volume {
 public:
  volume();

  ndim dim;
  vec origin;
  double a, inva;

  int nz() const { return num[2]; }
  int ny() const { return num[1]; }
  int nx() const { return num[0]; }
  int nr() const { return num[0]; }
  int num_direction(direction d) const {
    switch (d) {
    case X: return nx();
    case Y: return ny();
    case Z: return nz();
    case R: return nr();
    }
  };

  int has_field(component) const;
  int has_boundary(boundary_side,direction) const;

  vec dr() const;
  vec dx() const;
  vec dy() const;
  vec dz() const;

  int ntot() const { return the_ntot; }
  int index(component, const vec &) const; // DEPRECATED!
  void interpolate(component, const vec &, int indices[8], double weights[8]) const;

  void interpolate_one(component c, const vec &p,
                       int indices[8], double weights[8]) const;
  void interpolate_two(component c, const vec &p,
                       int indices[8], double weights[8]) const;
  void interpolate_cyl(component c, const vec &p, int m,
                       int indices[8], double weights[8]) const;

  double dv(component c, int index) const;
  volume dV(component c, int index) const;
  double intersection(const volume &) const;
  double rmin() const;
  double rmax() const;
  double xmin() const;
  double xmax() const;
  double ymin() const;
  double ymax() const;
  double zmin() const;
  double zmax() const;
  vec loc(component, int index) const;
  vec yee_shift(component) const;
  component eps_component() const;

  double boundary_location(boundary_side, direction);

  int contains(const vec &) const;
  int owns(const vec &) const;

  friend volume volcyl(double rsize, double zsize, double a);
  friend volume volone(double zsize, double a);
  friend volume voltwo(double xsize, double ysize, double a);

  int can_split_evenly(int num) const;
  volume split(int num, int which) const;
  volume split_once(int num, int which) const;
  volume split_specifically(int num, int which, direction d) const;
 private:
  volume(ndim, double a, int na, int nb=1, int nc=1);
  int num[3];
  int the_ntot;
};
