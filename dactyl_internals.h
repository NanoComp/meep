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

#define DOCMP for (int cmp=0;cmp<2-is_real;cmp++)

inline double max(double a, double b) { return (a > b) ? a : b; }
inline double min(double a, double b) { return (a < b) ? a : b; }
inline int max(int a, int b) { return (a > b) ? a : b; }
inline int min(int a, int b) { return (a < b) ? a : b; }

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
  polarizability(const mat *, double sig(const vec &),
                 double om, double ga, double sigscale,
                 double energy_saturation = 0.0);
  polarizability(const polarizability *);
  ~polarizability();
  double gamma, omeganot, *sigma, *s[10];
  double energy_saturation, saturated_sigma;
  polarizability *next;

  void use_pml();
};

class polarization {
 public:
  polarization(const polarizability *the_pb, int is_real);
  ~polarization();
  double saturation_factor;
  double *(P[10][2]), *(P_pml[10][2]), *(energy[10]), *(s[10]);
  int is_real;
  const polarizability *pb;
  polarization *next;

  double total_energy(const volume &);
  static polarization *set_up_polarizations(const mat *ma, int is_real);
  void use_real_fields();
};

class src {
 public:
  double freq, width, peaktime;
  complex<double> A[10], amp_shift;
  int i, cutoff;
  int is_continuous;
  src *next;
  int find_last_source(int guess=0);
  void use_real_sources();
  complex<double> get_amplitude_at_time(int t) const;
  double get_envelope_at_time(int t) const;
};

class bandsdata {
 public:
  bandsdata();
  ~bandsdata();

  complex<double> *f[10];
  // The following is the polarization at just one point, with Pz and Pp
  // added together (a crude compromize for speed, while still observing the
  // phonon bands).
  complex<double> *P;
  int tstart, tend, index, maxbands, scale_factor;
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
  complex<double> flux(fields *f);  
};
