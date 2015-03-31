/* Copyright (C) 2005-2015 Massachusetts Institute of Technology
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

#include "meep.hpp"

namespace meep {

#define DOCMP for (int cmp=0;cmp<2-is_real;cmp++)
#define DOCMP2 for (int cmp=0;cmp<2;cmp++)

inline double max(double a, double b) { return (a > b) ? a : b; }
inline double min(double a, double b) { return (a < b) ? a : b; }
inline int max(int a, int b) { return (a > b) ? a : b; }
inline int min(int a, int b) { return (a < b) ? a : b; }
static inline int abs(int a) { return a < 0 ? -a : a; }
static inline double abs(double a) { return fabs(a); }

// note that C99 has a round() function, but I don't want to rely on it
static inline int my_round(double x) {
  return int(floor(fabs(x) + 0.5) * (x < 0 ? -1 : 1));
}

inline int small_r_metal(int m) {
  return m-1;
}

inline int rmin_bulk(int m) {
  int r = 1 + small_r_metal(m);
  if (r < 1) r = 1;
  return r;
}

class src_vol {
 public:
  src_vol(component cc, src_time *st, int n, int *ind, std::complex<double> *amps);
  src_vol(const src_vol &sv);
  ~src_vol() { delete next; delete[] index; delete[] A;}

  src_time *t;
  int *index; // list of locations of sources in grid (indices)
  int npts; // number of points in list
  component c; // field component the source applies to
  std::complex<double> *A; // list of amplitudes

  std::complex<double> dipole(int j) { return A[j] * t->dipole(); }
  std::complex<double> current(int j) { return A[j] * t->current(); }
  void update(double time, double dt) { t->update(time, dt); }

  bool operator==(const src_vol &sv) const {
    return sv.index[0]==index[0] && sv.index[sv.npts-1]==index[npts-1] && sv.c==c && sv.t==t;
  }

  src_vol *add_to(src_vol *others);
  src_vol *next;
};

const int num_bandpts = 32;

class bandsdata {
 public:
  bandsdata();
  ~bandsdata();

  std::complex<double> *f[num_bandpts][NUM_FIELD_COMPONENTS];
  // The following is the polarization at just one point, with Pz and Pp
  // added together (a crude compromize for speed, while still observing the
  // phonon bands).
  std::complex<double> *P;
  int tstart, tend, index[num_bandpts], maxbands, scale_factor;
  fields_chunk *chunk[num_bandpts];
  double dt, fmin, fmax, qmin, fpmin;
  int ntime;
  int verbosity;

  int get_freqs(std::complex<double> *data, int n,
                std::complex<double> *amps, double *freqs, double *decays);
  int look_for_more_bands(std::complex<double> *simple_data,
                          double *reff, double *refd,
                          std::complex<double> *refa,
                          std::complex<double> *refdata,
                          int numref);
};

symmetry r_to_minus_r_symmetry(int m);

#define MIN_OUTPUT_TIME 4.0 // output no more often than this many seconds


// functions in step_generic.cpp:

void step_curl(realnum *f, component c, const realnum *g1, const realnum *g2,
	       int s1, int s2, // strides for g1/g2 shift
	       const grid_volume &gv, double dtdx,
	       direction dsig, const double *sig, const double *kap, const double *siginv,
	       realnum *fu, direction dsigu, const double *sigu, const double *kapu, const double *siginvu,
	       double dt, const realnum *cnd, const realnum *cndinv,
	       realnum *fcnd);

void step_update_EDHB(realnum *f, component fc, const grid_volume &gv,
		      const realnum *g, const realnum *g1, const realnum *g2,
		      const realnum *u, const realnum *u1, const realnum *u2,
		      int s, int s1, int s2,
		      const realnum *chi2, const realnum *chi3,
		      realnum *fw, direction dsigw, const double *sigw, const double *kapw);

void step_beta(realnum *f, component c, const realnum *g,
	       const grid_volume &gv, double betadt,
	       direction dsig, const double *siginv,
	       realnum *fu, direction dsigu, const double *siginvu,
	       const realnum *cndinv, realnum *fcnd);

// functions in step_generic_stride1.cpp, generated from step_generic.cpp:

void step_curl_stride1(realnum *f, component c, const realnum *g1, const realnum *g2,
	       int s1, int s2, // strides for g1/g2 shift
	       const grid_volume &gv, double dtdx,
	       direction dsig, const double *sig, const double *kap, const double *siginv,
	       realnum *fu, direction dsigu, const double *sigu, const double *kapu, const double *siginvu,
	       double dt, const realnum *cnd, const realnum *cndinv,
               realnum *fcnd);

void step_update_EDHB_stride1(realnum *f, component fc, const grid_volume &gv,
		      const realnum *g, const realnum *g1, const realnum *g2,
		      const realnum *u, const realnum *u1, const realnum *u2,
		      int s, int s1, int s2,
		      const realnum *chi2, const realnum *chi3,
		      realnum *fw, direction dsigw, const double *sigw, const double *kapw);

void step_beta_stride1(realnum *f, component c, const realnum *g,
		       const grid_volume &gv, double betadt,
		       direction dsig, const double *siginv,
		       realnum *fu, direction dsigu, const double *siginvu,
		       const realnum *cndinv, realnum *fcnd);

/* macro wrappers around time-stepping functions: for performance reasons,
   if the inner loop is stride-1 then we use the stride-1 versions,
   which allow gcc (and possibly other compilers) to do additional
   optimizations, especially loop vectorization */

#define STEP_CURL(f, c, g1, g2, s1, s2, gv, dtdx, dsig, sig, kap, siginv, fu, dsigu, sigu, kapu, siginvu, dt, cnd, cndinv, fcnd) do { \
  if (LOOPS_ARE_STRIDE1(gv))						\
    step_curl_stride1(f, c, g1, g2, s1, s2, gv, dtdx, dsig, sig, kap, siginv, fu, dsigu, sigu, kapu, siginvu, dt, cnd, cndinv, fcnd); \
  else									\
    step_curl(f, c, g1, g2, s1, s2, gv, dtdx, dsig, sig, kap, siginv, fu, dsigu, sigu, kapu, siginvu, dt, cnd, cndinv, fcnd); \
} while (0)

#define STEP_UPDATE_EDHB(f, fc, gv, g, g1, g2, u, u1, u2, s, s1, s2, chi2, chi3, fw, dsigw, sigw, kapw) do { \
  if (LOOPS_ARE_STRIDE1(gv))						\
    step_update_EDHB_stride1(f, fc, gv, g, g1, g2, u, u1, u2, s, s1, s2, chi2, chi3, fw, dsigw, sigw, kapw); \
  else									\
    step_update_EDHB(f, fc, gv, g, g1, g2, u, u1, u2, s, s1, s2, chi2, chi3, fw, dsigw, sigw, kapw); \
} while (0)

#define STEP_BETA(f, c, g, gv, betadt, dsig, siginv, fu, dsigu, siginvu, cndinv, fcnd) do {	\
  if (LOOPS_ARE_STRIDE1(gv))						\
    step_beta_stride1(f, c,g, gv, betadt, dsig, siginv, fu, dsigu, siginvu, cndinv, fcnd);	\
  else									\
    step_beta(f, c,g, gv, betadt, dsig, siginv, fu, dsigu, siginvu, cndinv, fcnd);		\
} while (0)

} // namespace meep
