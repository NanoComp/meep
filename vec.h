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

#ifndef VEC_H
#define VEC_H

#include <complex>

using namespace std;

enum component { Ex=0, Ey, Er, Ep, Ez, Hx, Hy, Hr, Hp, Hz };
enum ndim { D1=0, D2, D3, Dcyl };
enum field_type { E_stuff=0, H_stuff=1 };
enum boundary_side { High=0, Low };
enum direction { X=0,Y,Z,R,P };
struct signed_direction {
  signed_direction(direction dd=X,bool f=false) { d = dd; flipped = f; };
  direction d;
  bool flipped;
};

inline int number_of_directions(ndim dim) {
  return (int) (dim + 1 - 2 * (dim == Dcyl));
}

inline direction start_at_direction(ndim dim) {
  return (direction) (((dim == D1) || (dim == Dcyl)) ? 2 : 0);
}

inline direction stop_at_direction(ndim dim) {
  return (direction) (dim + 1 + 2 * (dim == D1));
}

#define LOOP_OVER_DIRECTIONS(dim, d) for (direction d = start_at_direction(dim), \
                                     loop_stop_directi = stop_at_direction(dim); \
                                     d < loop_stop_directi; d = (direction) (d+1))

inline signed_direction flip(signed_direction d) {
  signed_direction D2 = d;
  D2.flipped = !d.flipped;
  return D2;
}

inline bool has_direction(ndim dim, direction d) {
  switch (dim) {
  case D2: return d == X || d == Y;
  case D1: return d == Z;
  case D3: return d == X || d == Y || d == Z;
  case Dcyl: return d == Z || d == R;
  }
}

inline int is_electric(component c) { return (int) c < 5; }
inline int is_magnetic(component c) { return (int) c >= 5; }
inline field_type type(component c) {
  if (is_electric(c)) return E_stuff;
  else return H_stuff;
}
const char *component_name(component c);
const char *direction_name(direction);
const char *dimension_name(ndim);
inline direction component_direction(component c) {
  switch (c) {
  case Ex: case Hx: return X;
  case Ey: case Hy: return Y;
  case Ez: case Hz: return Z;
  case Er: case Hr: return R;
  case Ep: case Hp: return P;
  }
}
inline component direction_component(component c, direction d) {
  switch (d) {
  case X: if (is_electric(c)) return Ex; else return Hx;
  case Y: if (is_electric(c)) return Ey; else return Hy;
  case Z: if (is_electric(c)) return Ez; else return Hz;
  case R: if (is_electric(c)) return Er; else return Hr;
  case P: if (is_electric(c)) return Ep; else return Hp;
  }
}

class vec {
 public:
  vec() {};
  vec(double zz) { dim = D1; t[Z] = zz; };
  vec(double rr, double zz) { dim = Dcyl; t[R] = rr; t[Z] = zz; };
  vec(double xx, double yy, double zz) {
    dim = D3; t[X] = xx; t[Y] = yy; t[Z] = zz; };
  friend vec vec2d(double xx, double yy);
  ~vec() {};

  vec operator+(const vec &a) const {
    vec result = a;
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] += t[d];
    return result;
  };

  vec operator+=(const vec &a) {
    LOOP_OVER_DIRECTIONS(dim, d) t[d] += a.t[d];
  };

  vec operator-(const vec &a) const {
    vec result = a;
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] = t[d] - result.t[d];
    return result;
  };

  vec operator-=(const vec &a) {
    LOOP_OVER_DIRECTIONS(dim, d) t[d] -= a.t[d];
  };

  bool operator!=(const vec &a) const {
    LOOP_OVER_DIRECTIONS(dim, d) if (t[d] != a.t[d]) return true;
    return false;
  };

  bool operator==(const vec &a) const {
    LOOP_OVER_DIRECTIONS(dim, d) if (t[d] != a.t[d]) return false;
    return true;
  };

  vec operator*(double s) const {
    vec result = *this;
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] *= s;
    return result;
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
  vec v; v.dim = D2; v.t[X] = xx; v.t[Y] = yy; return v;
}

class ivec {
 public:
  ivec() { dim = D2; t[X] = t[Y] = 0; };
  ivec(int zz) { dim = D1; t[Z] = zz; };
  ivec(int rr, int zz) { dim = Dcyl; t[R] = rr; t[Z] = zz; };
  ivec(int xx, int yy, int zz) {
    dim = D3; t[X] = xx; t[Y] = yy; t[Z] = zz; };
  friend ivec ivec2d(int xx, int yy);
  ~ivec() {};

  ivec operator+(const ivec &a) const {
    switch (dim) {
    case Dcyl: return ivec(t[R]+a.t[R],t[Z]+a.t[Z]);
    case D3: return ivec(t[X]+a.t[X],t[Y]+a.t[Y],t[Z]+a.t[Z]);
    case D2: return ivec2d(t[X]+a.t[X],t[Y]+a.t[Y]);
    case D1: return ivec(t[Z]+a.t[Z]);
    }
  };
  ivec operator+=(const ivec &a) {
    switch (dim) {
    case Dcyl: t[R] += a.t[R]; t[Z] += a.t[Z]; return *this;
    case D3: t[X] += a.t[X]; t[Y] += a.t[Y]; t[Z] += a.t[Z]; return *this;
    case D2: t[X] += a.t[X]; t[Y] += a.t[Y]; return *this;
    case D1: t[Z] += a.t[Z]; return ivec(t[Z]+a.t[Z]);
    }
  };
  ivec operator-=(const ivec &a) {
    switch (dim) {
    case Dcyl: t[R] -= a.t[R]; t[Z] -= a.t[Z]; return *this;
    case D3: t[X] -= a.t[X]; t[Y] -= a.t[Y]; t[Z] -= a.t[Z]; return *this;
    case D2: t[X] -= a.t[X]; t[Y] -= a.t[Y]; return *this;
    case D1: t[Z] -= a.t[Z]; return ivec(t[Z]-a.t[Z]);
    }
  };
  ivec operator-(const ivec &a) const {
    switch (dim) {
    case Dcyl: return ivec(t[R]-a.t[R],t[Z]-a.t[Z]);
    case D3: return ivec(t[X]-a.t[X],t[Y]-a.t[Y],t[Z]-a.t[Z]);
    case D2: return ivec2d(t[X]-a.t[X],t[Y]-a.t[Y]);
    case D1: return ivec(t[Z]-a.t[Z]);
    }
  };
  bool operator!=(const ivec &a) const {
    switch (dim) {
    case Dcyl: return t[R]!=a.t[R] || t[Z]!=a.t[Z];
    case D3: return t[X]!=a.t[X] || t[Y]!=a.t[Y] || t[Z]!=a.t[Z];
    case D2: return t[X]!=a.t[X] || t[Y]!=a.t[Y];
    case D1: return t[Z]!=a.t[Z];
    }
  };
  bool operator==(const ivec &a) const {
    switch (dim) {
    case Dcyl: return t[R]==a.t[R] && t[Z]==a.t[Z];
    case D3: return t[X]==a.t[X] && t[Y]==a.t[Y] && t[Z]==a.t[Z];
    case D2: return t[X]==a.t[X] && t[Y]==a.t[Y];
    case D1: return t[Z]==a.t[Z];
    }
  };
  ivec operator*(int s) const {
    switch (dim) {
    case Dcyl: return ivec(t[R]*s,t[Z]*s);
    case D3: return ivec(t[X]*s,t[Y]*s,t[Z]*s);
    case D2: return ivec2d(t[X]*s,t[Y]*s);
    case D1: return ivec(t[Z]*s);
    }
  };
  vec operator*(double s) const {
    switch (dim) {
    case Dcyl: return vec(t[R]*s,t[Z]*s);
    case D3: return vec(t[X]*s,t[Y]*s,t[Z]*s);
    case D2: return vec2d(t[X]*s,t[Y]*s);
    case D1: return vec(t[Z]*s);
    }
  };
  ndim dim;

  int r() const { return t[R]; };
  int x() const { return t[X]; };
  int y() const { return t[Y]; };
  int z() const { return t[Z]; };
  int in_direction(direction d) const { return t[d]; };
  int set_direction(direction d, int val) { t[d] = val; };

  friend ivec zero_ivec(ndim);
 private:
  int t[5];
};

inline ivec zero_ivec(ndim di) {
  ivec v; v.dim = di; for (int d=0;d<5;d++) v.set_direction((direction)d,0);
  return v;
}

inline ivec ivec2d(int xx, int yy) {
  ivec v; v.t[X] = xx; v.t[Y] = yy; return v;
}

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
  vec operator[](const ivec &p) const { return p*(0.5*inva); };
  int index(component, const ivec &) const;
  ivec round_vec(const vec &) const;
  void interpolate(component, const vec &, int indices[8], double weights[8]) const;

  void interpolate_one(component c, const vec &p,
                       int indices[8], double weights[8]) const;
  void interpolate_two(component c, const vec &p,
                       int indices[8], double weights[8]) const;
  void interpolate_cyl(component c, const vec &p, int m,
                       int indices[8], double weights[8]) const;

  double dv(component c, int index) const;
  volume dV(component c, int index) const;
  volume dV(const ivec &) const;
  double intersection(const volume &) const;
  double rmin() const;
  double rmax() const;
  double xmin() const;
  double xmax() const;
  double ymin() const;
  double ymax() const;
  double zmin() const;
  double zmax() const;
  vec center() const;
  vec loc(component, int index) const;
  vec loc_at_resolution(int index, double res) const;
  int ntot_at_resolution(double res) const;
  ivec iloc(component, int index) const;
  vec yee_shift(component) const;
  component eps_component() const;

  double boundary_location(boundary_side, direction) const;
  ivec big_corner() const;
  ivec little_corner() const { return io(); };

  int contains(const vec &) const;
  int contains(const ivec &) const;
  int owns(const ivec &) const;

  friend volume volcyl(double rsize, double zsize, double a);
  friend volume volone(double zsize, double a);
  friend volume voltwo(double xsize, double ysize, double a);

  int can_split_evenly(int num) const;
  volume split(int num, int which) const;
  volume split_once(int num, int which) const;
  volume split_at_fraction(bool want_high, int numer) const;
  volume split_specifically(int num, int which, direction d) const;
  volume pad(direction d) const;
 private:
  volume(ndim, double a, int na, int nb=1, int nc=1);
  ivec io() const;
  ivec iyee_shift(component) const;
  int num[3];
  int the_ntot;
};

class symmetry {
 public:
  symmetry();
  symmetry(const symmetry &);
  ~symmetry();
  friend symmetry identity();
  friend symmetry rotate4(direction,const volume &);
  friend symmetry rotate2(direction,const volume &);
  friend symmetry mirror(direction,const volume &);

  signed_direction transform(direction d, int n) const;
  ivec transform(const ivec &, int n) const;
  vec transform(const vec &, int n) const;
  component transform(component, int n) const;
  complex<double> phase_shift(component, int n) const;
  int multiplicity() const;
  bool is_primitive(const ivec &) const;

  symmetry operator+(const symmetry &) const;
  void operator=(const symmetry &);
 private:
  signed_direction S[5];
  vec symmetry_point;
  double a, inva;
  int g; // g is the multiplicity of the symmetry.
  symmetry *next;
};

symmetry identity();
symmetry rotate4(direction,const volume &);
symmetry rotate2(direction,const volume &);
symmetry mirror(direction,const volume &);

#endif
