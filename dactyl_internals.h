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

#define MA(e,r,z) ((e)[(z)+(r)*(nz+1)])
#define DOCMP for (int cmp=0;cmp<2;cmp++)

#define RE(f,r,z) ((f)[0][(z)+(r)*(nz+1)])
#define IM(f,r,z) ((f)[1][(z)+(r)*(nz+1)])

#define PMLR(f,r,z) ((f)[cmp][(z)+((r)-nr+npmlr)*(nz+1)])
#define PMLZ(f,r) ((f)[cmp][(lr+1)/2][(iz)+(r)*npmlz])
#define CM(f,r,z) ((f)[cmp][(z)+(r)*(nz+1)])
#define IT(f,r,z) (cmp?((f)[0][(z)+(r)*(nz+1)]):(-(f)[1][(z)+(r)*(nz+1)]))

#define EIKZ(f,r,z) (cmp? (cosknz*IM(f,r,z)+sinknz*RE(f,r,z)) \
                         : (cosknz*RE(f,r,z)-sinknz*IM(f,r,z)))
#define EMIKZ(f,r,z) (cmp? (cosknz*IM(f,r,z)-sinknz*RE(f,r,z)) \
                          : (cosknz*RE(f,r,z)+sinknz*IM(f,r,z)))
#define IEIKZ(f,r,z) (cmp? (cosknz*RE(f,r,z)+sinknz*IM(f,r,z)) \
                          : (-(cosknz*IM(f,r,z)-sinknz*RE(f,r,z))))

#define FIPHI(f,phase) (cmp ? (cos(phase)*imag(f)+sin(phase)*real(f)) \
                            : (cos(phase)*real(f)-sin(phase)*imag(f)))
#define FW(f,ipos) &((f)[0][ipos*ifreqmax]), &((f)[1][ipos*ifreqmax])
#define FPW(f,ipos,freq) (complex<double>((f)[0][(ipos)*ifreqmax+(freq)], \
                                  (f)[1][(ipos)*ifreqmax+(freq)]))
#define SWAP(a,b) {(a) += (b); (b) = (a)-(b); (a) -= (b); }

inline double max(double a, double b) { return (a > b) ? a : b; }
inline int max(int a, int b) { return (a > b) ? a : b; }

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
  int nr, nz, npmlr, npmlz;
  polarizability(const mat *, double sig(double,double),
                 double om, double ga, double sigscale);
  polarizability(const polarizability *);
  ~polarizability();
  double gamma, omeganot, *sigma, *sr, *sp, *sz;
  polarizability *next;

  void use_pml(int npmlr, int npmlz);
};

class polarization {
 public:
  polarization(const polarizability *the_pb);
  ~polarization();
  double *(Pr[2]), *(Pp[2]), *(Pz[2]);
  // The PML split polarization fields:
  double *(Prp[2]), *(Ppz[2]), *(Pzr[2]);
  double *(z_Prp[2][2]), *(z_Ppz[2][2]);
  const polarizability *pb;
  polarization *next;

  static polarization *set_up_polarizations(const mat *ma);
};

class src {
 public:
  double freq, width, peaktime;
  complex<double> ez, ep, er, amp_shift;
  int r, z, cutoff;
  int is_real;
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

  complex<double> *hr, *hp, *hz, *er, *ep, *ez;
  // The following is the polarization at just one point, with Pz and Pp
  // added together (a crude compromize for speed, while still observing the
  // phonon bands).
  complex<double> *P;
  int tstart, tend, z, nr, maxbands, scale_factor;
  double a, inva, fmin, fmax, qmin;
  int ntime;

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
