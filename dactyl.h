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

#include <complex>

#define MAXFLUXPLANES 20

class mat {
 public:
  double *eps, *invepser, *invepsez, *invepsep, a;
  double *Crez, *Crep, *Crhz, *Crhp;
  double *Cper, *Cpez, *Cphr, *Cphz;
  double *Czer, *Czep, *Czhr, *Czhp;
  long *pols;
  int nr, nz, numpols;
  int npmlr, npmlz; // Amount of pml
  double *freq[32], *width[32], height[32];
  const char *outdir;

  ~mat();
  mat(double eps(double r, double z),
      double rmax, double zmax, double a=1.0);
  void use_pml(int npmlr=16, int npmlz=16);

  void output_slices(const char *name);
  void set_output_directory(const char *name);
 private:
  void output_sigma_slice(const char *name);
};

class src {
 public:
  double freq, width, peaktime;
  double ez, ep, er;
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

  void get_fields(complex<double> *eigen, double *f, double *d,
                  int nbands, int n);
  int get_both_freqs(complex<double> *data1, complex<double> *data2, int n,
                     complex<double> *amps1, complex<double> *amps2, 
                     double *freqs, double *decays);
  int get_freqs(complex<double> *data, int n,
                complex<double> *amps, double *freqs, double *decays);
};

class fields {
 public:
  double *(hr[2]), *(hp[2]), *(hz[2]), *(er[2]), *(ep[2]), *(ez[2]);
  double *(hrp[2]), *(hpz[2]), *(hzr[2]), *(erp[2]), *(epz[2]), *(ezr[2]);
  double *(z_hrp[2][2]), *(z_hpz[2][2]), *(z_erp[2][2]), *(z_epz[2][2]);
  double *pol[32];
  double a, inva; // The "lattice constant" and its inverse!
  int nr, nz, numpols;
  int npmlr, npmlz; // Amount of pml
  int m, t;
  double k, cosknz, sinknz;
  bandsdata bands;
  src *e_sources, *h_sources;
  const mat *ma;
  const char *outdir;

  fields(const mat *, int m);
  void use_bloch(double kz);
  ~fields();

  void output_slices(const char *name);
  void output_real_imaginary_slices(const char *name);
  void step();

  void use_real_sources();
  void add_er_source(double freq, double width, double peaktime,
                     double cutoff, int z, double amp(double r));
  void add_ep_source(double freq, double width, double peaktime,
                     double cutoff, int z, double amp(double r));
  void add_ez_source(double freq, double width, double peaktime,
                     double cutoff, int z, double amp(double r));
  void add_hr_source(double freq, double width, double peaktime,
                     double cutoff, int z, double amp(double r));
  void add_hp_source(double freq, double width, double peaktime,
                     double cutoff, int z, double amp(double r));
  void add_hz_source(double freq, double width, double peaktime,
                     double cutoff, int z, double amp(double r));

  void output_point(FILE *, double r, double z, const char *name);

  void prepare_for_bands(int z, int ttot, double fmax, double qmin=1e300);
  void record_bands();
  double get_band(int n, double fmax, int maxbands=100);
  void output_bands(FILE *, const char *, int maxbands=100);
  void output_bands_and_modes(FILE *, const char *, int maxbands=100);
  double total_energy();
  double zflux(int ri, int ro, int z);
  double rflux(int zl, int zu, int r);
  void dft_flux();
  int add_zfluxplane(int ri, int ro, int z);
  int add_rfluxplane(int zl, int zu, int r);
  int set_frequency_range(double wl, double wu, double deltaw);
  void ttow(complex<double> field, double *retarget, double *imtarget, double time);
  void fluxw_output(FILE *outpf, char *header);
  void set_output_directory(const char *name);
 private: 
  double *(erw[2]), *(epw[2]), *(ezw[2]), *(hrw[2]), *(hpw[2]), *(hzw[2]);
  int iposmax, ifreqmax, nfreq, nzflux, *(nzfluxplane[MAXFLUXPLANES]);
  int nrflux, *(nrfluxplane[MAXFLUXPLANES]);
  double *freqs;
  void step_h_bulk();
  void step_h_pml();
  void step_h_boundaries();
  void step_h_source(const src *);
  void step_e_bulk();
  void step_e_pml();
  void step_e_boundaries();
  void step_e_source(const src *);
  void add_src_pt(int r, int z,
                  double Pr, double Pp, double Pz,
                  double freq, double width, double peaktime,
                  double cutoff, int is_h = 0);
  int setifreqmax_and_iposmax(int ifreq, int ipos);
  void out_bands(FILE *, const char *, int maxbands, int outmodes);
};

const double c = 0.5;

// The following is a utility function to parse the executable name use it
// to come up with a directory name, avoiding overwriting any existing
// directory, unless the source file hasn't changed.

const char *make_output_directory(const char *exename);
FILE *create_output_file(const char *dirname, const char *fname);
