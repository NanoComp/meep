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

inline int is_electric(component c) { return (int) c < 5; }
inline int is_magnetic(component c) { return (int) c >= 5; }
const char *component_name(component c);

class vec {
 public:
  vec() { dim = d2; tx = ty = 0; };
  vec(double zz) { dim = d1; tz = zz; };
  vec(double rr, double zz) { dim = dcyl; tr = rr; tz = zz; };
  ~vec() {};

  vec operator+(const vec &a) const {
    switch (dim) {
    case dcyl: return vec(tr+a.tr,tz+a.tz);
    case d1: return vec(tz+a.tz);
    }
  };
  vec operator-(const vec &a) const {
    switch (dim) {
    case dcyl: return vec(tr-a.tr,tz-a.tz);
    case d1: return vec(tz-a.tz);
    }
  };
  vec operator==(const vec &a) const;
  vec operator*(double s) const {
    switch (dim) {
    case dcyl: return vec(tr*s,tz*s);
    case d1: return vec(tz*s);
    }
  };
  ndim dim;

  double r() const { return tr; };
  double x() const { return tx; };
  double y() const { return ty; };
  double z() const { return tz; };

  void print(FILE *) const;
 private:
  double tx, ty, tz, tr;
};

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

  int has_field(component) const;

  vec dr() const;
  vec dz() const;

  int ntot() const { return the_ntot; }
  int index(component, const vec &) const; // DEPRECATED!
  void interpolate(component, const vec &, int indices[8], double weights[8]) const;

  void interpolate_cyl(component c, const vec &p, int m,
                       int indices[8], double weights[8]) const;

  double dv(component c, int index) const;
  vec loc(component, int index) const;
  vec yee_shift(component) const;
  component eps_component() const;

  int contains(const vec &) const;
  int owns(const vec &) const;

  friend volume volcyl(double rsize, double zsize, double a);
  friend volume volone(double zsize, double a);
 private:
  volume(ndim, double a, int na, int nb=1, int nc=1);

  int num[3];
  int the_ntot;
};
