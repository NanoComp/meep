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
#ifndef DACTYL_H
#define DACTYL_H

#include <stdio.h>

#include "vec.h"
#include "mympi.h"

const double c = 0.5;
const double pi = 3.141592653589793238462643383276L;

class polarizability;
class polarization;
class grace;

class mat_chunk {
 public:
  double *eps, a;
  double *inveps[10][5];
  double *C[5][10];
  double *Cdecay[5][10][5];
  volume v;
  geometric_volume gv;
  polarizability *pb;

  ~mat_chunk();
  mat_chunk(const volume &v, double eps(const vec &),
            const geometric_volume &vol_limit, int proc_num=0);
  mat_chunk(const mat_chunk *);
  void make_average_eps();
  void use_pml(direction, double dx, double boundary_loc);
  void update_pml_arrays();
  void update_Cdecay();

  void set_output_directory(const char *name);
  void mix_with(const mat_chunk *, double);

  void add_polarizability(double sigma(const vec &), double omega, double gamma,
                          double delta_epsilon = 1.0, double energy_saturation = 0.0);
  int n_proc() const { return the_proc; } // Says which proc owns me!
  int is_mine() const { return the_is_mine; }

  // monitor.cpp
  void interpolate_eps(const vec &loc, double val[8]) const;
  double max_eps() const;
 private:
  double pml_fmin;
  int the_proc;
  int the_is_mine;
};

class mat {
 public:
  mat_chunk **chunks;
  int num_chunks;
  int desired_num_chunks;
  volume v, user_volume;
  geometric_volume gv;
  symmetry S;
  const char *outdir;
  volume *effort_volumes;
  double *effort;
  int num_effort_volumes;

  ~mat();
  mat();
  mat(const volume &v, double eps(const vec &), int num_chunks = 0,
      const symmetry &s = identity());
  mat(const mat *);
  mat(const mat &);
  void add_to_effort_volumes(const volume &new_effort_volume, double extra_effort);
  void redefine_chunks(const int Nv, const volume *new_volumes, const int *procs);
  void optimize_volumes(int *Nv, volume *new_volumes, int *procs);
  void optimize_chunks();
  void choose_chunkdivision(const volume &v, double eps(const vec &),
                            int num_chunks = 1,
                            const symmetry &s = identity());
  void make_average_eps();
  void use_pml(direction d, boundary_side b, double dx, bool recalculate_chunks = true);
  void use_pml_everywhere(double dx, bool recalculate_chunks = true);

  void output_slices(const char *name = "") const;
  void output_slices(const geometric_volume &what, const char *name = "") const;
  void set_output_directory(const char *name);
  void mix_with(const mat *, double);

  void add_polarizability(double sigma(const vec &), double omega, double gamma,
                          double delta_epsilon = 1.0, double energy_saturation = 0.0);

  // monitor.cpp
  double get_eps(const vec &loc) const;
  double max_eps() const;
 private:
};

class src;
class bandsdata;
class fields_chunk;

class monitor_point {
 public:
  monitor_point();
  ~monitor_point();
  vec loc;
  double t;
  complex<double> f[10];
  monitor_point *next;

  complex<double> get_component(component);

  // When called with only its first four arguments, fourier_transform
  // performs an FFT on its monitor points, putting the frequencies in f
  // and the amplitudes in a.  Yes, the frequencies are trivial and
  // redundant, but this saves you the risk of making a mistake in
  // converting your units.  Note also, that in this case f is always a
  // real number, although it's stored in a float.
  //
  // Note that in either case, fourier_transform assumes that the monitor
  // points are all equally spaced in time.
  void fourier_transform(component w,
                         complex<double> **a, complex<double> **f, int *numout,
                         double fmin=0.0, double fmax=0.0, int maxbands=100);
  // harminv works much like fourier_transform, except that it is not yet
  // implemented.
  void harminv(component w,
               complex<double> **a, double **f_re, double **f_im,
               int *numout, double fmin, double fmax,
               int maxbands);
};

class partial_flux_plane;

class flux_plane {
 public:
  partial_flux_plane *partials;
  flux_plane *next;

  flux_plane(partial_flux_plane *);
  ~flux_plane();

  double flux();
};

enum in_or_out { Incoming=0, Outgoing=1 };

class fields_chunk {
 public:
  double *(f[10][2]);
  double *(f_backup[10][2]);
  double *(f_p_pml[10][2]);
  double *(f_m_pml[10][2]);
  double *(f_backup_p_pml[10][2]);
  double *(f_backup_m_pml[10][2]);

  double **zeroes[2]; // Holds pointers to metal points.
  int num_zeroes[2];
  double **(connections[2][2]);
  int num_connections[2][2];
  complex<double> *connection_phases[2];
  double *connection_factors[2];

  polarization *pol, *olpol;
  partial_flux_plane *fluxes;
  double a, inva; // The "lattice constant" and its inverse!
  volume v;
  geometric_volume gv;
  int m, is_real;
  bandsdata *bands;
  src *e_sources, *h_sources;
  const mat_chunk *new_ma;
  mat_chunk *ma;
  const char *outdir;

  fields_chunk(const mat_chunk *, const char *outdir, int m=0);
  ~fields_chunk();

  // step.cpp
  void step();
  double peek_field(component, const vec &);

  void use_real_fields();
  double find_last_source();
  // monitor.cpp
  complex<double> get_field(component, const ivec &) const;
  double get_eps(const ivec &iloc) const;
  complex<double> analytic_epsilon(double freq, const vec &) const;
  
  // slices.cpp
  double maxfieldmag(component) const;
  void output_eps_body(component, const symmetry &, int sn,
                       const geometric_volume &what, file *);
  complex<double> fields_chunk::field_mean(component c, bool abs_real,
                                           bool abs_imag) const;

  double electric_energy_in_box(const geometric_volume &, const symmetry &);
  double magnetic_energy_in_box(const geometric_volume &, const symmetry &);
  double thermo_energy_in_box(const geometric_volume &, const symmetry &);

  double backup_h();
  double restore_h();

  void set_output_directory(const char *name);
  void verbose(int v=1) { verbosity = v; }

  double count_volume(component);
  friend class fields;

  int n_proc() const { return ma->n_proc(); };
  int is_mine() const { return ma->is_mine(); };
  // boundaries.cpp
  void zero_metal(field_type);
  // fluxes.cpp
  partial_flux_plane *new_flux_plane(const vec &p1, const vec &p2);
  void update_fluxes();
 private: 
  int verbosity; // Turn on verbosity for debugging purposes...
  // fields.cpp
  void figure_out_step_plan();
  bool have_plus_deriv[10], have_minus_deriv[10];
  component plus_component[10], minus_component[10];
  direction plus_deriv_direction[10], minus_deriv_direction[10];
  int num_each_direction[3], stride_each_direction[3];
  // bands.cpp
  void record_bands(int tcount);
  // step.cpp
  void step_h_old();
  void step_e_old();
  void phase_in_material(const mat_chunk *ma);
  void phase_material(int phasein_time);
  void step_h();
  void step_h_source(const src *, double);
  void step_e();
  void step_polarization_itself(polarization *old = NULL, polarization *newp = NULL);
  void step_e_polarization(polarization *old = NULL, polarization *newp = NULL);
  void step_e_source(const src *, double);
  // polarization.cpp
  void prepare_step_polarization_energy(polarization *op = NULL, polarization *np = NULL);
  void half_step_polarization_energy(polarization *op = NULL, polarization *np = NULL);
  void update_polarization_saturation(polarization *op = NULL, polarization *np = NULL);
  // fields.cpp
  void alloc_f(component c);
  // monitor.cpp
  // sources.cpp

  // add_point_source returns 1 if the connections between chunks need to
  // be recalculated.  This allows us to avoid allocating TE or TM fields
  // until we know which is desired.
  int add_point_source(component whichf, double freq, double width, double peaktime,
                       double cutoff, const ivec &, complex<double> amp,
                       int is_continuous, double time);
  void add_indexed_source(component whichf, double freq, double width,
                          double peaktime, double cutoff, int theindex, 
                          complex<double> amp, int is_c, double time);
  // initialize.cpp
  void initialize_field(component, complex<double> f(const vec &));
  void initialize_polarizations(polarization *op=NULL, polarization *np=NULL);
  void initialize_with_nth_te(int n, double kz);
  void initialize_with_nth_tm(int n, double kz);
  // boundaries.cpp
  void alloc_extra_connections(field_type, in_or_out, int);
  // fluxes.cpp
  partial_flux_plane *nfp_1d(const vec &);
};

enum boundary_condition { Periodic=0, Metallic, Magnetic, None };
enum time_sink { Connecting, Stepping, Boundaries, MpiTime,
                 Slicing, Other };

class fields {
 public:
  int num_chunks;
  fields_chunk **chunks;
  flux_plane *fluxes;
  symmetry S;
  // The following is an array that is num_chunks by num_chunks.  Actually
  // it is two arrays, one for the imaginary and one for the real part.
  double **comm_blocks[2];
  // This is the same size as each comm_blocks array, and stores the sizes
  // of the comm blocks themselves.
  int *comm_sizes[2];
  int *comm_num_complex[2];
  int *comm_num_negate[2];
  double a, inva; // The "lattice constant" and its inverse!
  volume v, user_volume;
  geometric_volume gv;
  int m, t, phasein_time, is_real;
  complex<double> k[5], eikna[5];
  double coskna[5], sinkna[5];
  boundary_condition boundaries[2][5];
  bandsdata *bands;
  const char *outdir;
  // fields.cpp methods:
  fields(const mat *, int m=0);
  ~fields();
  void use_real_fields();
  // time.cpp
  double time_spent_on(time_sink);
  void print_times();
  // boundaries.cpp
  void set_boundary(boundary_side,direction,
                    boundary_condition, bool autoconnect=true,
                    complex<double> kcomponent=0.0);
  void use_metal_everywhere();
  void use_bloch(direction d, double k, bool autoconnect=true) {
    use_bloch(d, (complex<double>) k, autoconnect);
  }
  void use_bloch(direction, complex<double> kz, bool autoconnect=true);
  void use_bloch(const vec &k, bool autoconnect=true);
  vec lattice_vector(direction) const;
  // slices.cpp methods:
  void output_slices(const char *name = "");
  void output_slices(const geometric_volume &what, const char *name = "");
  void eps_slices(const char *name = "");
  void eps_slices(const vec &origin, const vec &xside, const vec &yside,
                  const double dx = 0.05, const char *name = "");
  void eps_slices(const geometric_volume &what, const char *name = "");
  void output_real_imaginary_slices(const char *name = "");
  void output_real_imaginary_slices(const geometric_volume &what,
                                    const char *name = "");
  double maxfieldmag_to_master(component) const;
  complex<double> optimal_phase_shift(component) const;
  // step.cpp methods:
  void step();
  void step_old(); // Step using old code.
  inline double time() const { return t*inva*c; };

  double find_last_source();
  void add_point_source(component whichf, double freq, double width, double peaktime,
                        double cutoff, const vec &, complex<double> amp = 1.0,
                        int is_continuous = 0);
  void initialize_field(component, complex<double> f(const vec &));
  void initialize_with_nth_te(int n);
  void initialize_with_nth_tm(int n);
  void initialize_with_n_te(int n);
  void initialize_with_n_tm(int n);
  void initialize_polarizations();
  int phase_in_material(const mat *ma, double time);
  int is_phasing();

  // monitor.cpp
  double get_eps(const vec &loc) const;
  void get_point(monitor_point *p, const vec &) const;
  monitor_point *get_new_point(const vec &, monitor_point *p=NULL) const;
  void output_point(file *, const vec &, const char *name);
  complex<double> analytic_epsilon(double freq, const vec &) const;
  
  void prepare_for_bands(const vec &, double end_time, double fmax=0,
                         double qmin=1e300, double frac_pow_min=0.0);
  void record_bands();
  complex<double> get_band(int n, int maxbands=100);
  void grace_bands(grace *, int maxbands=100);
  void output_bands(file *, const char *, int maxbands=100);
  // energy_and_flux.cpp
  double energy_in_box(const geometric_volume &);
  double electric_energy_in_box(const geometric_volume &);
  double magnetic_energy_in_box(const geometric_volume &);
  double thermo_energy_in_box(const geometric_volume &);
  double total_energy();
  double field_energy_in_box(const geometric_volume &);
  double field_energy();

  void set_output_directory(const char *name);
  void verbose(int v=1) { verbosity = v; }
  double count_volume(component);
  // fluxes.cpp
  flux_plane *add_flux_plane(const vec &, const vec &);
 private: 
  int verbosity; // Turn on verbosity for debugging purposes...
  unsigned long last_time;
  time_sink working_on, was_working_on;
  double times_spent[Other+1];
  // field.cpp
  bool have_component(component);
  // material.cpp
  double max_eps() const;
  // time.cpp
#ifdef WITH_TIMINGS
  void am_now_working_on(time_sink);
  void finished_working();
#else
  void am_now_working_on(time_sink) { return; }
  void finished_working() { return; }
#endif
  // boundaries.cpp
  void find_metals();
  void disconnect_chunks();
  void connect_chunks();
  void connect_the_chunks(); // Intended to be ultra-private...
  int is_metal(const ivec &);
  ivec ilattice_vector(direction) const;
  bool locate_component_point(component *, ivec *, complex<double> *) const;
  bool locate_point_in_user_volume(ivec *, complex<double> *phase) const;
  // step.cpp
  void phase_material();
  void step_h();
  void step_h_old();
  void step_h_source();
  void step_e();
  void step_e_old();
  void step_boundaries(field_type);
  void step_polarization_itself();
  void step_e_polarization();
  void step_e_source();
  void prepare_step_polarization_energy();
  void half_step_polarization_energy();
  void update_polarization_saturation();
  int cluster_some_bands_cleverly(double *tf, double *td, complex<double> *ta,
                                  int num_freqs, int fields_considered, int maxbands,
                                  complex<double> *fad, double *approx_power);
  void out_bands(file *, const char *, int maxbands);
  complex<double> *clever_cluster_bands(int maxbands, double *approx_power = NULL);
  // sources.cpp
  void add_point_source(component whichf, double freq, double width, double peaktime,
                        double cutoff, const ivec &p, complex<double> amp, int is_c);
  // slices.cpp
  void outline_chunks(file *name);
  bool has_eps_interface(vec *loc) const;
  complex<double> fields::field_mean(component c, bool abs_real = false,
                                     bool abs_imag = false) const;
  // energy_and_flux.cpp
  // fluxes.cpp
  void update_fluxes();
  // monitor.cpp
  complex<double> get_field(component c, const ivec &iloc) const;
  complex<double> get_field(component c, const vec &loc) const;
  double get_eps(const ivec &iloc) const;
};

class grace_point;
enum grace_type { XY, ERROR_BARS };

class grace {
 public:
  grace(const char *fname, const char *dirname = ".");
  ~grace();
  
  void new_set(grace_type t = XY);
  void new_curve();
  void set_legend(const char *);
  void set_range(double xmin, double xmax, double ymin, double ymax);
  void output_point(double x, double y,
                    double dy = -1.0, double extra = -1.0);
  void output_out_of_order(int n, double x, double y,
                           double dy = -1.0, double extra= -1.0);
 private:
  void flush_pts();
  file *f;
  char *fn, *dn;
  grace_point *pts;
  int set_num,sn;
};

// The following is a utility function to parse the executable name use it
// to come up with a directory name, avoiding overwriting any existing
// directory, unless the source file hasn't changed.

const char *make_output_directory(const char *exename, const char *jobname = NULL);
void trash_output_directory(const char *dirname);
file *create_output_file(const char *dirname, const char *fname);

// The following allows you to hit ctrl-C to tell your calculation to stop
// and clean up.
void deal_with_ctrl_c(int stop_now = 2);
// When a ctrl_c is called, the following variable (which starts with a
// zero value) is incremented.
extern int interrupt;

int do_harminv(complex<double> *data, int n, int sampling_rate, double a,
	       double fmin, double fmax, int maxbands,
	       complex<double> *amps, double *freq_re, double *freq_im, double *errors = NULL);

#endif
