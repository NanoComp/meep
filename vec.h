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

#define FOR_ELECTRIC_COMPONENTS(c) for (component c = Ex; \
                                        c < Hx; c = (component) (c+1))
#define FOR_MAGNETIC_COMPONENTS(c) for (component c = Hz; \
                                        c > Ez; c = (component) (c-1))
#define FOR_COMPONENTS(c) for (component c = Ex,loop_stop_co=Ey; \
                               c != loop_stop_co; \
                               c = (component)((c+1)%10), \
                               loop_stop_co = Ex)
#define FOR_DIRECTIONS(d) for (direction d = X,loop_stop_di=Y; \
                               d != loop_stop_di; \
                               d = (direction)((d+1)%5), \
                               loop_stop_di = X)

#define LOOP_OVER_DIRECTIONS(dim, d) for (direction d = start_at_direction(dim), \
                                     loop_stop_directi = stop_at_direction(dim); \
                                     d < loop_stop_directi; d = (direction) (d+1))

inline signed_direction flip(signed_direction d) {
  signed_direction D2 = d;
  D2.flipped = !d.flipped;
  return D2;
}

inline bool has_direction(ndim dim, direction d) {
  LOOP_OVER_DIRECTIONS(dim, dd) if (dd == d) return true;
  return false;
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

class file;

class vec {
 public:
  vec() {};
  vec(ndim di) { dim = di; };
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
  void print(file *) const;
  friend vec zero_vec(ndim);
 private:
  double t[5];
};

inline vec zero_vec(ndim di) {
  vec v; v.dim = di; LOOP_OVER_DIRECTIONS(di, d) v.set_direction(d, 0.0);
  return v;
}

inline vec vec2d(double xx, double yy) {
  vec v; v.dim = D2; v.t[X] = xx; v.t[Y] = yy; return v;
}

class ivec {
 public:
  ivec() { dim = D2; t[X] = t[Y] = 0; };
  ivec(ndim di) { dim = di; };
  ivec(int zz) { dim = D1; t[Z] = zz; };
  ivec(int rr, int zz) { dim = Dcyl; t[R] = rr; t[Z] = zz; };
  ivec(int xx, int yy, int zz) {
    dim = D3; t[X] = xx; t[Y] = yy; t[Z] = zz; };
  friend ivec ivec2d(int xx, int yy);
  ~ivec() {};

  ivec operator+(const ivec &a) const {
    ivec result = a;
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] += t[d];
    return result;
  };

  ivec operator+=(const ivec &a) {
    LOOP_OVER_DIRECTIONS(dim, d) t[d] += a.t[d];
  };

  ivec operator-(const ivec &a) const {
    ivec result = a;
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] = t[d] - result.t[d];
    return result;
  };

  ivec operator-=(const ivec &a) {
    LOOP_OVER_DIRECTIONS(dim, d) t[d] -= a.t[d];
  };

  bool operator!=(const ivec &a) const {
    LOOP_OVER_DIRECTIONS(dim, d) if (t[d] != a.t[d]) return true;
    return false;
  };

  bool operator==(const ivec &a) const {
    LOOP_OVER_DIRECTIONS(dim, d) if (t[d] != a.t[d]) return false;
    return true;
  };

  ivec operator*(int s) const {
    ivec result = *this;
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] *= s;
    return result;
  };

  vec operator*(double s) const {
    vec result(dim);
    LOOP_OVER_DIRECTIONS(dim, d) result.set_direction(d, t[d] * s);
    return result;
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
  ivec v; v.dim = di; LOOP_OVER_DIRECTIONS(di, d) v.set_direction(d, 0);
  return v;
}

inline ivec ivec2d(int xx, int yy) {
  ivec v; v.t[X] = xx; v.t[Y] = yy; return v;
}

class geometric_volume {
 public:
  ndim dim;
  geometric_volume(ndim di) { dim = di; };
  geometric_volume(const vec &vec1, const vec &vec2);
  void set_direction_min(direction d, double val) { min_corner.set_direction(d, val); };
  void set_direction_max(direction d, double val) { max_corner.set_direction(d, val); };
  double in_direction_min(direction d) const { return min_corner.in_direction(d); };
  double in_direction_max(direction d) const { return max_corner.in_direction(d); };
  double computational_volume(); 
  double full_volume(); 
  geometric_volume intersect_with(const geometric_volume &a);
 private:
  vec min_corner, max_corner;
};

class volume {
 public:
  volume() {};

  ndim dim;
  vec origin;
  double a, inva;

  int num_direction(direction d) const {
    return num[((int) d) % 3];
  };
  int nr() const { return num_direction(R); }
  int nx() const { return num_direction(X); }
  int ny() const { return num_direction(Y); }
  int nz() const { return num_direction(Z); }

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

  bool contains(const vec &) const;
  bool contains(const ivec &) const;
  bool owns(const ivec &) const;

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
