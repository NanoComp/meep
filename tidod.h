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
#ifndef TIDOD_H
#define TIDOD_H

#include <complex>
#include "dactyl.h"

enum component_1d { Ex=0, Hy };

class polarizability_1d;
class polarization_1d;
class bandsdata_1d;
class grace;

class monitor_point_1d {
 public:
  monitor_point_1d();
  monitor_point_1d(double r, double z, const fields *f);
  ~monitor_point_1d();
  double z, t;
  complex<double> hy, ex;
  monitor_point_1d *next;

  complex<double> get_component(component_1d c) {
    return (c == Hy)?hy:ex;
  };

  // When called with only its first four arguments, fourier_transform
  // performs an FFT on its monitor points, putting the frequencies in f
  // and the amplitudes in a.  Yes, the frequencies are trivial and
  // redundant, but this saves you the risk of making a mistake in
  // converting your units.  Note also, that in this case f is always a
  // real number, although it's stored in a float.
  //
  // Note that in either case, fourier_transform assumes that the monitor
  // points are all equally spaced in time.
  void fourier_transform(component_1d w,
                         complex<double> **a, complex<double> **f, int *numout,
                         double fmin=0.0, double fmax=0.0, int maxbands=100);
  // harminv works much like fourier_transform, except that it is not yet
  // implemented.
  void harminv(component_1d w,
               complex<double> **a, double **f_re, double **f_im,
               int *numout, double fmin, double fmax,
               int maxbands);
};

class mat_1d {
 public:
  double *eps, *inveps, a;
  double *Czex, *Czhy;
  int nz;
  int npmlz; // Amount of pml
  polarizability_1d *pb;
  const char *outdir;

  ~mat_1d();
  mat_1d(double eps(double z), double zmax, double a=1.0);
  mat_1d(const mat_1d *);
  void make_average_eps();
  double use_integer_pml(int npmlz=16, double fmin=0.2);
  double use_pml(double pmlz=2.0, double fmin=0.2);

  void output_slices(const char *name = "");
  void set_output_directory(const char *name);
  void mix_with(const mat_1d *, double);

  void add_polarizability(double sigma(double), double omega, double gamma,
                          double delta_epsilon = 1.0);
  void add_plasma(double oneorzero(double), double omega_plasma, double gamma);
 private:
  void output_sigma_slice(const char *name);
  double pml_fmin;
};

class src_1d;
class fields_1d;

class fields_1d {
 public:
  double *(hy[2]), *(ex[2]);
  double *(backup_hy[2]);
  polarization_1d *pol, *olpol;
  double a, inva; // The "lattice constant" and its inverse!
  int nz;
  int npmlz; // Amount of pml
  int t, phasein_time;
  double k, cosknz, sinknz;
  complex<double> eiknz;
  bandsdata_1d *bands;
  src_1d *e_sources, *h_sources;
  const mat_1d *new_ma;
  mat_1d *ma;
  const char *outdir;
  double preferred_fmax;

  fields_1d(const mat_1d *);
  void use_bloch(double kz);
  ~fields_1d();

  void output_slices(const char *name = "");
  void energy_slices(const char *name = "");
  void output_real_imaginary_slices(const char *name = "");
  void step();
  inline double time() { return t*inva*c; };

  void use_real_sources();
  int find_last_source();
  // Note that the following plane source only works if m == 1.
  void add_source(component_1d whichf, double freq, double width, double peaktime,
                  double cutoff, double z, int is_continuous = 0);
  void add_continuous_source(component_1d whichf, double freq, double width, double peaktime,
                             double cutoff, double z);
  void initialize_with_nth(int n);
  void initialize_with_n(int n);
  void initialize_polarizations(polarization_1d *op=NULL, polarization_1d *np=NULL);
  int phase_in_material(const mat_1d *ma, double time);
  int is_phasing();

  monitor_point_1d *get_new_point(double z, monitor_point_1d *p=NULL);
  void get_point(monitor_point_1d *, double z);
  void output_point(FILE *, double z, const char *name);

  void prepare_for_bands(double total_time, double fmax=0,
                         double qmin=1e300, double frac_pow_min=0.0);
  void record_bands();
  complex<double> get_band(int n, int maxbands=100);
  void grace_bands(grace *, int maxbands=100);
  double energy_in_box(double zmin, double zmax);
  double thermo_energy_in_box(double zmin, double zmax);
  double electric_energy_in_box(double zmin, double zmax);
  double magnetic_energy_in_box(double zmin, double zmax);
  double total_energy() {return energy_in_box(0.0, nz*inva);};
  void set_output_directory(const char *name);
  void verbose(int v=1) { verbosity = v; }
 private: 
  int verbosity; // Turn on verbosity for debugging purposes...
  double *freqs;
  void find_source_z_position(double z, double shift, int *z1, int *z2,
                              complex<double> *amp1, complex<double> *amp2);
  void phase_material();
  void step_h();
  void step_h_source(const src_1d *);
  void step_e();
  void step_polarization_itself(polarization_1d *old = NULL, polarization_1d *newp = NULL);
  void half_step_polarization_energy(polarization_1d *old = NULL,
                                     polarization_1d *newp = NULL);
  void prepare_step_polarization_energy(polarization_1d *old = NULL,
                                        polarization_1d *newp = NULL);
  void step_e_polarization(polarization_1d *old = NULL, polarization_1d *newp = NULL);
  void step_e_source(const src_1d *);
  void add_src_pt(int z, complex<double> amp,
                  double freq, double width, double peaktime,
                  double cutoff, int is_h = 0, int is_continuous = 0);
  int cluster_some_bands_cleverly(double *tf, double *td, complex<double> *ta,
                                  int num_freqs, int fields_considered, int maxbands,
                                  complex<double> *fad, double *approx_power);
  void out_bands(FILE *, const char *, int maxbands, int outmodes);
  complex<double> *get_the_bands(int maxbands, double *approx_power = NULL);
  complex<double> *clever_cluster_bands(int maxbands, double *approx_power = NULL);
};

#endif
