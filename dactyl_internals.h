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

class polarizability {
 public:
  int nr, nz;
  polarizability(const mat *, double sig(double,double),
                 double om, double ga, double sigscale);
  polarizability(const polarizability *);
  ~polarizability();
  double gamma, omeganot, *sigma, *sr, *sp, *sz;
  polarizability *next;
};

class polarization {
 public:
  polarization(const int nr, const int nz, const polarizability *the_pb);
  ~polarization();
  double *(Pr[2]), *(Pp[2]), *(Pz[2]);
  const polarizability *pb;
  polarization *next;

  static polarization *set_up_polarizations(const mat *ma);
  void step_polarization(const polarization *newer_pol);
};

class src {
 public:
  double freq, width, peaktime;
  complex<double> ez, ep, er;
  int r, z, cutoff;
  int is_real;
  src *next;
  int find_last_source(int guess=0);
  void use_real_sources();
};

class bandsdata {
 public:
  bandsdata();
  ~bandsdata();

  complex<double> *hr, *hp, *hz, *er, *ep, *ez;
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
};
