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

#define MA(e,z) ((e)[(z)])
#define DOCMP for (int cmp=0;cmp<2;cmp++)

#define RE(f,z) ((f)[0][(z)])
#define IM(f,z) ((f)[1][(z)])

#define PMLZ(f) ((f)[cmp][(lr+1)/2][(iz)])
#define CM(f,z) ((f)[cmp][(z)])
#define IT(f,z) (cmp?((f)[0][(z)]):(-(f)[1][(z)]))

#define EIKZ(f,z) (cmp? (cosknz*IM(f,z)+sinknz*RE(f,z)) \
                      : (cosknz*RE(f,z)-sinknz*IM(f,z)))
#define EMIKZ(f,z) (cmp? (cosknz*IM(f,z)-sinknz*RE(f,z)) \
                       : (cosknz*RE(f,z)+sinknz*IM(f,z)))
#define IEIKZ(f,z) (cmp?   (cosknz*RE(f,z)+sinknz*IM(f,z)) \
                       : (-(cosknz*IM(f,z)-sinknz*RE(f,z))))

inline double max(double a, double b) { return (a > b) ? a : b; }
inline double min(double a, double b) { return (a < b) ? a : b; }
inline int max(int a, int b) { return (a > b) ? a : b; }
inline int min(int a, int b) { return (a < b) ? a : b; }

class polarizability_1d {
 public:
  int nr, nz, npmlr, npmlz;
  polarizability_1d(const mat_1d *, double sig(double),
                 double om, double ga, double sigscale);
  polarizability_1d(const polarizability_1d *);
  ~polarizability_1d();
  double gamma, omeganot, *sigma;
  polarizability_1d *next;

  void use_integer_pml(int npmlz);
};

class polarization_1d {
 public:
  polarization_1d(const polarizability_1d *the_pb);
  ~polarization_1d();
  double *(Px[2]);
  const polarizability_1d *pb;
  polarization_1d *next;

  static polarization_1d *set_up_polarizations(const mat_1d *ma);
};

class src_1d {
 public:
  double freq, width, peaktime;
  complex<double> amp, amp_shift;
  int z, cutoff;
  int is_real, is_continuous;
  src_1d *next;
  int find_last_source(int guess=0);
  void use_real_sources();
  complex<double> get_amplitude_at_time(int t) const;
  double get_envelope_at_time(int t) const;
};

class bandsdata_1d {
 public:
  bandsdata_1d();
  ~bandsdata_1d();

  complex<double> *hy, *hz, *ex, *ez;
  int tstart, tend, nz, maxbands, scale_factor;
  double a, inva, fmin, fmax, qmin, fpmin;
  int ntime;
  int verbosity;

  void get_fields(complex<double> *eigen, complex<double> *f_and_d,
                  int nbands, int n);
  int get_both_freqs(complex<double> *data1, complex<double> *data2, int n,
                     complex<double> *amps1, complex<double> *amps2, 
                     double *freqs, double *decays);
  int get_freqs(complex<double> *data, int n,
                complex<double> *amps, double *freqs, double *decays);
  int look_for_more_bands(complex<double> *simple_data,
                          double *reff, double *refd,
                          complex<double> *refa,
                          complex<double> *refdata,
                          int numref);
};
