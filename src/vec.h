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

#ifndef MEEP_VEC_H
#define MEEP_VEC_H

#include <complex>

using namespace std;

namespace meep {

const int NUM_FIELD_COMPONENTS = 15;
const int NUM_FIELD_TYPES = 4;

enum component { Ex=0, Ey, Er, Ep, Ez, Hx, Hy, Hr, Hp, Hz,
                 Dx, Dy, Dr, Dp, Dz, Dielectric };
enum ndim { D1=0, D2, D3, Dcyl };
enum field_type { E_stuff=0, H_stuff=1, D_stuff=2, P_stuff=3 };
enum boundary_side { High=0, Low };
enum direction { X=0,Y,Z,R,P, NO_DIRECTION };
struct signed_direction {
  signed_direction(direction dd=X,bool f=false, complex<double> ph=1.0) {
    d = dd; flipped = f; phase = ph;
  };
  signed_direction(const signed_direction &sd) {
    d = sd.d; flipped = sd.flipped; phase = sd.phase;
  }
  signed_direction operator*(complex<double> ph);
  direction d;
  bool flipped;
  complex<double> phase;
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

#define FOR_FIELD_TYPES(ft) for (field_type ft = E_stuff; \
                                 ft <= P_stuff; ft = (field_type) (ft+1))
#define FOR_ELECTRIC_COMPONENTS(c) for (component c = Ex; \
                                        c < Hx; c = (component) (c+1))
#define FOR_MAGNETIC_COMPONENTS(c) for (component c = Hz; \
                                        c > Ez; c = (component) (c-1))
#define FOR_D_COMPONENTS(c) for (component c = Dz; \
                                 c > Hz; c = (component) (c-1))
#define FOR_E_AND_D(e,d) for (component e = Ex, d = Dx; \
                              e <= Ez; e = (component) (e+1), d = (component) (d+1))
#define FOR_COMPONENTS(c) for (component c = Ex,loop_stop_co=Ey; \
                               c != loop_stop_co; \
                               c = (component)((c+1)%NUM_FIELD_COMPONENTS), \
                               loop_stop_co = Ex)
#define FOR_DIRECTIONS(d) for (direction d = X,loop_stop_di=Y; \
                               d != loop_stop_di; \
                               d = (direction)((d+1)%5), \
                               loop_stop_di = X)

#define LOOP_OVER_DIRECTIONS(dim, d) for (direction d = start_at_direction(dim), \
                                     loop_stop_directi = stop_at_direction(dim); \
                                     d < loop_stop_directi; d = (direction) (d+1))

// loop over indices idx from is to ie (inclusive) in v
#define LOOP_OVER_IVECS(v, is, ie, idx) \
  for (int loop_n1 = (ie.yucky_val(0) - is.yucky_val(0)) / 2 + 1, \
           loop_n2 = (ie.yucky_val(1) - is.yucky_val(1)) / 2 + 1, \
           loop_n3 = (ie.yucky_val(2) - is.yucky_val(2)) / 2 + 1, \
           loop_d1 = v.yucky_direction(0), \
           loop_d2 = v.yucky_direction(1), \
           loop_d3 = v.yucky_direction(2), \
           loop_is1 = is.yucky_val(0), \
           loop_is2 = is.yucky_val(1), \
           loop_is3 = is.yucky_val(2), \
           loop_s1 = v.stride((direction) loop_d1), \
           loop_s2 = v.stride((direction) loop_d2), \
           loop_s3 = v.stride((direction) loop_d3), \
           idx0 = (is - v.little_corner()).yucky_val(0) / 2 * loop_s1 \
                + (is - v.little_corner()).yucky_val(1) / 2 * loop_s2 \
                + (is - v.little_corner()).yucky_val(2) / 2 * loop_s3,\
           loop_i1 = 0; loop_i1 < loop_n1; loop_i1++) \
    for (int loop_i2 = 0; loop_i2 < loop_n2; loop_i2++) \
      for (int idx = idx0 + loop_i1*loop_s1 + loop_i2*loop_s2, \
           loop_i3 = 0; loop_i3 < loop_n3; loop_i3++, idx+=loop_s3)

// integration weight for using LOOP_OVER_IVECS with field::integrate
#define IVEC_LOOP_WEIGHTx(i, n, dir) ((i > 1 && i < n - 2) ? 1.0 : (i == 0 ? s0.in_direction(direction(dir)) : (i == 1 ? s1.in_direction(direction(dir)) : i == n - 1 ? e0.in_direction(direction(dir)) : (i == n - 2 ? e1.in_direction(direction(dir)) : 1.0))))
#define IVEC_LOOP_WEIGHT(k) IVEC_LOOP_WEIGHTx(loop_i##k, loop_n##k, loop_d##k)

#define LOOP_OVER_OWNED(v, idx) \
  for (int loop_n1 = v.yucky_num(0), \
           loop_n2 = v.yucky_num(1), \
           loop_n3 = v.yucky_num(2), \
           loop_s1 = v.stride(v.yucky_direction(0)), \
           loop_s2 = v.stride(v.yucky_direction(1)), \
           loop_s3 = v.stride(v.yucky_direction(2)), \
           loop_i1 = 0; loop_i1 < loop_n1; loop_i1++) \
      for (int loop_i2 = 0; loop_i2 < loop_n2; loop_i2++) \
        for (int idx = loop_i1*loop_s1 + loop_i2*loop_s2, \
                 loop_i3 = 0; loop_i3 < loop_n3; loop_i3++, idx+=loop_s3)

inline signed_direction flip(signed_direction d) {
  signed_direction d2 = d;
  d2.flipped = !d.flipped;
  return d2;
}

inline bool has_direction(ndim dim, direction d) {
  LOOP_OVER_DIRECTIONS(dim, dd) if (dd == d) return true;
  return false;
}

// true if d is polar while dim is cartesian, or vice versa 
inline bool coordinate_mismatch(ndim dim, direction d) {
  return ((dim >= D1 && dim <= D3 && (d == R || d == P)) ||
	  (dim == Dcyl && (d == X || d == Y)));
}

void abort(const char *, ...);

inline int is_electric(component c) { return c < Hx; }
inline int is_magnetic(component c) { return c >= Hx && c < Dx; }
inline int is_D(component c) { return c >= Dx && c < Dielectric; }
inline field_type type(component c) {
  if (is_electric(c)) return E_stuff;
  else if (is_magnetic(c)) return H_stuff;
  else if (is_D(c)) return D_stuff;
  abort("Invalid field in type.\n");
  return E_stuff; // This is never reached.
}
const char *component_name(component c);
const char *direction_name(direction);
const char *dimension_name(ndim);

inline direction component_direction(component c) {
  switch (c) {
  case Ex: case Hx: case Dx: return X;
  case Ey: case Hy: case Dy: return Y;
  case Ez: case Hz: case Dz: return Z;
  case Er: case Hr: case Dr: return R;
  case Ep: case Hp: case Dp: return P;
  case Dielectric: return NO_DIRECTION;
  }
  return X; // This code is never reached...
}
inline component direction_component(component c, direction d) {
  component start_point;
  if (is_electric(c)) start_point = Ex;
  else if (is_magnetic(c)) start_point = Hx;
  else if (is_D(c)) start_point = Dx;
  else if (c == Dielectric && d == NO_DIRECTION) return Dielectric;
  else abort("unknown field component %d", c);
  switch (d) {
  case X: return start_point;
  case Y: return (component) (start_point + 1);
  case Z: return (component) (start_point + 4);
  case R: return (component) (start_point + 2);
  case P: return (component) (start_point + 3);
  case NO_DIRECTION: abort("vector %d component in NO_DIRECTION", c);
  }
  return Ex; // This is never reached.
}

class file;

class vec {
 public:
  vec() {};
  vec(ndim di) { dim = di; };
  vec(ndim di, double val) { dim = di; t[0]=t[1]=t[2]=t[3]=t[4]=val; };
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
    return *this;
  };

  vec operator-(const vec &a) const {
    vec result = a;
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] = t[d] - result.t[d];
    return result;
  };

  vec operator-(void) const {
    vec result(dim);
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] = -t[d];
    return result;
  };

  vec operator-=(const vec &a) {
    LOOP_OVER_DIRECTIONS(dim, d) t[d] -= a.t[d];
    return *this;
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

  vec operator/(double s) const {
    vec result = *this;
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] *= (1.0/s);
    return result;
  };

  // I use & as a dot product.
  double operator&(const vec &a) const {
    double result = 0.0;
    LOOP_OVER_DIRECTIONS(dim, d) result += t[d] * a.t[d];
    return result;
  };
  ndim dim;

  double r() const { return t[R]; };
  double x() const { return t[X]; };
  double y() const { return t[Y]; };
  double z() const { return t[Z]; };
  double in_direction(direction d) const { return t[d]; };
  void set_direction(direction d, double val) { t[d] = val; };

  double project_to_boundary(direction, double boundary_loc);
  void print(file *) const;
  friend vec zero_vec(ndim);
 private:
  double t[5];
};

inline double abs(const vec &v) { return sqrt(v & v); }

inline vec zero_vec(ndim di) {
  vec v; v.dim = di; LOOP_OVER_DIRECTIONS(di, d) v.set_direction(d, 0.0);
  return v;
}

inline vec clean_vec(const vec &v, double val_unused = 0.0) {
  vec vc(v.dim, val_unused);
  LOOP_OVER_DIRECTIONS(v.dim, d) vc.set_direction(d, v.in_direction(d));
  return vc;
}

inline vec vec2d(double xx, double yy) {
  vec v; v.dim = D2; v.t[X] = xx; v.t[Y] = yy; return v;
}

class ivec {
 public:
  ivec() { dim = D2; t[X] = t[Y] = 0; };
  ivec(ndim di) { dim = di; };
  ivec(ndim di, int val) { dim = di; t[0]=t[1]=t[2]=t[3]=t[4]=val; };
  ivec(int zz) { dim = D1; t[Z] = zz; };
  ivec(int rr, int zz) { dim = Dcyl; t[R] = rr; t[Z] = zz; };
  ivec(int xx, int yy, int zz) {
    dim = D3; t[X] = xx; t[Y] = yy; t[Z] = zz; };
  friend ivec ivec2d(int xx, int yy);
  ~ivec() {};

  // Only an idiot (or a macro) would use a yucky function.  Don't be an
  // idiot.
  int yucky_val(int) const;

  ivec operator+(const ivec &a) const {
    ivec result = a;
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] += t[d];
    return result;
  };

  ivec operator+=(const ivec &a) {
    LOOP_OVER_DIRECTIONS(dim, d) t[d] += a.t[d];
    return *this;
  };

  ivec operator-(const ivec &a) const {
    ivec result = a;
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] = t[d] - result.t[d];
    return result;
  };

  ivec operator-(void) const {
    ivec result(dim);
    LOOP_OVER_DIRECTIONS(dim, d) result.t[d] = -t[d];
    return result;
  };

  ivec operator-=(const ivec &a) {
    LOOP_OVER_DIRECTIONS(dim, d) t[d] -= a.t[d];
    return *this;
  };

  bool operator!=(const ivec &a) const {
    LOOP_OVER_DIRECTIONS(dim, d) if (t[d] != a.t[d]) return true;
    return false;
  };

  bool operator==(const ivec &a) const {
    LOOP_OVER_DIRECTIONS(dim, d) if (t[d] != a.t[d]) return false;
    return true;
  };

  bool operator<=(const ivec &a) const {
    LOOP_OVER_DIRECTIONS(dim, d) if (t[d] > a.t[d]) return false;
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
  void set_direction(direction d, int val) { t[d] = val; };

  friend ivec zero_ivec(ndim);
  friend ivec one_ivec(ndim);
 private:
  int t[5];
};

inline ivec zero_ivec(ndim di) {
  ivec v; v.dim = di; LOOP_OVER_DIRECTIONS(di, d) v.set_direction(d, 0);
  return v;
}

inline ivec one_ivec(ndim di) {
  ivec v; v.dim = di; LOOP_OVER_DIRECTIONS(di, d) v.set_direction(d, 1);
  return v;
}

inline ivec unit_ivec(ndim di, direction d) {
  ivec v(zero_ivec(di));
  v.set_direction(d, 1);
  return v;
}

inline ivec ivec2d(int xx, int yy) {
  ivec v; v.t[X] = xx; v.t[Y] = yy; return v;
}

class geometric_volume {
 public:
  ndim dim;
  geometric_volume(ndim di) { dim = di; min_corner.dim = di; max_corner.dim = di; };
  geometric_volume(const vec &vec1, const vec &vec2);
  void set_direction_min(direction d, double val) { min_corner.set_direction(d, val); };
  void set_direction_max(direction d, double val) { max_corner.set_direction(d, val); };
  double in_direction_min(direction d) const { return min_corner.in_direction(d); };
  double in_direction_max(direction d) const { return max_corner.in_direction(d); };
  double computational_volume(); 
  double full_volume() const;
  bool contains(const vec &h) const;
  geometric_volume intersect_with(const geometric_volume &a) const;
  geometric_volume operator&(const geometric_volume &a) const {
    return intersect_with(a);
  };
  geometric_volume operator+(const vec &a) const {
    return geometric_volume(min_corner + a, max_corner + a);
  }
  geometric_volume operator+=(const vec &a) {
    min_corner += a; max_corner += a;
    return *this;
  }
  geometric_volume operator-(const vec &a) const {
    return geometric_volume(min_corner - a, max_corner - a);
  }
  geometric_volume operator-=(const vec &a) {
    min_corner -= a; max_corner -= a;
    return *this;
  }
  bool intersects(const geometric_volume &a) const;
  bool operator&&(const geometric_volume &a) const {
    return intersects(a);
  };
  vec get_min_corner() const { return min_corner; };
  vec get_max_corner() const { return max_corner; };
 private:
  vec min_corner, max_corner;
};

class volume {
 public:
  volume() {};

  ndim dim;
  vec origin;
  double a, inva;

  void print() const;
  int stride(direction d) const { return the_stride[d]; };
  int num_direction(direction d) const {
    return num[((int) d) % 3];
  };
  // Only an idiot (or a macro) would use a yucky function.  Don't be an
  // idiot.
  int yucky_num(int) const;
  direction yucky_direction(int) const;
  void set_num_direction(direction d, int value);
  int nr() const { return num_direction(R); }
  int nx() const { return num_direction(X); }
  int ny() const { return num_direction(Y); }
  int nz() const { return num_direction(Z); }

  bool has_field(component c) const {
    if (dim == D1) return c == Ex || c == Hy || c == Dx;
    return (dim == Dcyl)?component_direction(c)>Y:component_direction(c)<R;
  }
  int has_boundary(boundary_side,direction) const;

  vec dr() const;
  vec dx() const;
  vec dy() const;
  vec dz() const;

  int ntot() const { return the_ntot; }
  int nowned() const { int n = 1; LOOP_OVER_DIRECTIONS(dim,d) n *= num_direction(d); return n; }
  vec operator[](const ivec &p) const { return p*(0.5*inva); };
  int index(component, const ivec &) const;
  ivec round_vec(const vec &) const;
  void interpolate(component, const vec &, int indices[8], double weights[8]) const;
  void interpolate(component, const vec &, ivec locs[8], double weights[8]) const;

  geometric_volume dV(component c, int index) const;
  geometric_volume dV(const ivec &) const;
  bool intersect_with(const volume &vol_in, volume *intersection = NULL, volume *others = NULL, int *num_others = NULL) const;
  double rmin() const;
  double rmax() const;
  double xmin() const;
  double xmax() const;
  double ymin() const;
  double ymax() const;
  double zmin() const;
  double zmax() const;
  vec center() const;
  ivec icenter() const;
  vec loc(component, int index) const;
  vec loc_at_resolution(int index, double res) const;
  int ntot_at_resolution(double res) const;
  ivec iloc(component, int index) const;

  int yee_index(component c) const {
    int idx = 0;
    LOOP_OVER_DIRECTIONS(dim,d)
      idx += (1-iyee_shift(c).in_direction(d))*stride(d);
    return idx;
  }
  vec yee_shift(component) const;
  component eps_component() const;
  void yee2diel_offsets(component c, int &offset1, int &offset2);

  double boundary_location(boundary_side, direction) const;
  ivec big_corner() const;
  ivec little_corner() const { return io(); };
  vec corner(boundary_side b) const;

  bool contains(const vec &) const;
  bool contains(const ivec &) const;
  bool owns(const ivec &) const;
  geometric_volume surroundings() const;

  friend volume volcyl(double rsize, double zsize, double a);
  friend volume volone(double zsize, double a);
  friend volume vol1d(double zsize, double a);
  friend volume voltwo(double xsize, double ysize, double a);
  friend volume vol2d(double xsize, double ysize, double a);
  friend volume vol3d(double xsize, double ysize, double zsize, double a);

  int can_split_evenly(int num) const;
  volume split(int num, int which) const;
  volume split_by_effort(int num, int which, int Ngv = 0, const volume *gv = NULL, double *effort = NULL) const;
  volume split_once(int num, int which) const;
  volume split_at_fraction(bool want_high, int numer) const;
  volume split_specifically(int num, int which, direction d) const;
  volume pad(direction d) const;
  volume pad() const {
       volume v = *this;
       LOOP_OVER_DIRECTIONS(dim,d)
	    v = v.pad(d);
       return v;
  }
  ivec iyee_shift(component c) const {
    ivec out = zero_ivec(dim);
    LOOP_OVER_DIRECTIONS(dim,d)
      if (c == Dielectric ||
          ((is_electric(c) || is_D(c)) && d == component_direction(c)) ||
          (is_magnetic(c) && d != component_direction(c)))
        out.set_direction(d,1);
    return out;
  }
 private:
  volume(ndim, double a, int na, int nb=1, int nc=1);
  ivec io() const;
  void set_strides();
  int num[3];
  int the_stride[5];
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
  geometric_volume transform(const geometric_volume &, int n) const;
  component transform(component, int n) const;
  complex<double> phase_shift(component, int n) const;
  int multiplicity() const;
  bool is_primitive(const ivec &) const;

  symmetry operator+(const symmetry &) const;
  symmetry operator*(double) const;
  void operator=(const symmetry &);
 private:
  signed_direction S[5];
  complex<double> ph;
  vec symmetry_point;
  ivec i_symmetry_point;
  int g; // g is the multiplicity of the symmetry.
  symmetry *next;
  friend symmetry r_to_minus_r_symmetry(int m);
};

symmetry identity();
symmetry rotate4(direction,const volume &);
symmetry rotate2(direction,const volume &);
symmetry mirror(direction,const volume &);

} /* namespace meep */

#endif /* MEEP_VEC_H */
