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

#include "meep.h"

namespace meep {

#define DOCMP for (int cmp=0;cmp<2-is_real;cmp++)

inline double max(double a, double b) { return (a > b) ? a : b; }
inline double min(double a, double b) { return (a < b) ? a : b; }
inline int max(int a, int b) { return (a > b) ? a : b; }
inline int min(int a, int b) { return (a < b) ? a : b; }

inline ivec min(const ivec &a, const ivec &b) {
  ivec v(a.dim);
  LOOP_OVER_DIRECTIONS(a.dim, d) 
    v.set_direction(d, min(a.in_direction(d), b.in_direction(d)));
  return v;
}
inline ivec max(const ivec &a, const ivec &b) {
  ivec v(a.dim);
  LOOP_OVER_DIRECTIONS(a.dim, d) 
    v.set_direction(d, max(a.in_direction(d), b.in_direction(d)));
  return v;
}

inline int small_r_metal(int m) {
  return m-1;
}

inline int rmin_bulk(int m) {
  int r = 1 + small_r_metal(m);
  if (r < 1) r = 1;
  return r;
}

class polarizability {
 public:
  volume v;
  polarizability(const structure_chunk *, material_function &sig,
                 double om, double ga, double sigscale,
                 double energy_saturation = 0.0, bool mine = true);
  polarizability(const polarizability *);
  ~polarizability();
  double gamma, omeganot, *sigma, *s[NUM_FIELD_COMPONENTS];
  double energy_saturation, saturated_sigma;
  bool is_mine() { return is_it_mine; };
  bool is_it_mine;
  polarizability *next;

  polarizability_identifier get_identifier() const;
};

class polarization {
 public:
  polarization(const polarizability *the_pb, int is_real);
  ~polarization();
  double saturation_factor;
  double *(P[NUM_FIELD_COMPONENTS][2]), *(energy[NUM_FIELD_COMPONENTS]),
    *(s[NUM_FIELD_COMPONENTS]);
  int is_real;
  const polarizability *pb;
  polarization *next;

  complex<double> analytic_epsilon(double freq, const vec &) const;
  double local_energy(const ivec &);
  double total_energy(const geometric_volume &);
  static polarization *set_up_polarizations(const structure_chunk *s, int is_real);
  void use_real_fields();
  void initialize_energy(double energy(const vec &));
};

class src_vol {
 public:
  src_vol(component cc, src_time *st, int n, const int *ind, const complex<double> *amps);
  src_vol(const src_vol &sv);
  ~src_vol() { delete next; delete[] index; delete[] A;}

  src_time *t;
  int *index; // list of locations of sources in grid (indices)
  int npts; // number of points in list
  component c; // field component the source applies to
  complex<double> *A; // list of amplitudes

  complex<double> dipole(int j) { return A[j] * t->dipole(); }
  complex<double> dipole(int j, double time) { return A[j] * t->dipole(time); }
  complex<double> current(int j) { return A[j] * t->current(); }
  complex<double> current(int j, double T, double dt) { return A[j] * t->current(T,dt); }
  void update(double time, double dt) { t->update(time, dt); }

  bool operator==(const src_vol &sv) const {
    return sv.index[0]==index[0] && sv.index[sv.npts-1]==index[npts-1] && sv.c==c && sv.t==t;
  }

  src_vol *add_to(src_vol *others) const;
  src_vol *next;
};

const int num_bandpts = 32;

class bandsdata {
 public:
  bandsdata();
  ~bandsdata();

  complex<double> *f[num_bandpts][NUM_FIELD_COMPONENTS];
  // The following is the polarization at just one point, with Pz and Pp
  // added together (a crude compromize for speed, while still observing the
  // phonon bands).
  complex<double> *P;
  int tstart, tend, index[num_bandpts], maxbands, scale_factor;
  fields_chunk *chunk[num_bandpts];
  double a, inva, fmin, fmax, qmin, fpmin;
  int ntime;
  int verbosity;

  int get_freqs(complex<double> *data, int n,
                complex<double> *amps, double *freqs, double *decays);
  int look_for_more_bands(complex<double> *simple_data,
                          double *reff, double *refd,
                          complex<double> *refa,
                          complex<double> *refdata,
                          int numref);
};

class partial_flux_plane {
 public:
  partial_flux_plane(fields_chunk *, int);
  ~partial_flux_plane();

  component cE, cH;
  double *weights;
  int *indE, *indH, numpts;
  partial_flux_plane *next, *next_in_chunk;
  fields_chunk *f;

  double *oldE[2];

  void append(partial_flux_plane *);
  void append_in_chunk(partial_flux_plane *);
 private:
};

class weighted_flux_plane {
 public:
  int ymin, ymax, xconst;
  double dy_min, dy_max;
  int is_rflux;
  int verbosity;
  weighted_flux_plane() {};
  weighted_flux_plane(int ymin, int ymax, int xconst, 
	    double dy_min, double dy_max, int is_rflux);
  ~weighted_flux_plane() {};
  complex<double> flux(fields_chunk *f);  
};

symmetry r_to_minus_r_symmetry(int m);

} // namespace meep
