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

enum component { Er=0, Ep, Ez, Hr, Hp, Hz };

const double c = 0.5;
const double pi = 3.141592653589793238462643383276L;

class polarizability;
class polarization;
class grace;

class mat {
 public:
  double *eps, *invepser, *invepsez, *invepsep, a;
  double *Crez, *Crep, *Crhz, *Crhp;
  double *Cper, *Cpez, *Cphr, *Cphz;
  double *Czer, *Czep, *Czhr, *Czhp;
  int nr, nz;
  int npmlr, npmlz; // Amount of pml
  polarizability *pb;
  const char *outdir;

  ~mat();
  mat(double eps(double r, double z),
      double rmax, double zmax, double a=1.0);
  mat(const mat *);
  void make_average_eps();
  void use_pml(int npmlr=16, int npmlz=16);

  void output_slices(const char *name = "");
  void set_output_directory(const char *name);
  void mix_with(const mat *, double);

  void add_polarizability(double sigma(double,double), double omega, double gamma,
                          double delta_epsilon = 1.0);
 private:
  void output_sigma_slice(const char *name);
};

class src;
class bandsdata;
class fields;
class weighted_flux_plane;

class flux_plane {
 public:
  double ymin, ymax, xconst;
  int is_rflux;
  int num_wf;
  double weights[2];
  int xpos[2];
  int verbosity;
  weighted_flux_plane *wf[2];
  flux_plane(double ymin, double ymax, double xconst, int is_rflux, double a);
  flux_plane(const flux_plane &fp);
  ~flux_plane();
  complex<double> flux(fields *f);
};

class a_field {
 public:
  a_field();
  ~a_field();
  a_field(complex<double>, complex<double>, complex<double>);
  void another_field(complex<double>, complex<double>, complex<double>);
  complex<double> r, p, z;
  a_field *next;
};

class monitor_point {
 public:
  monitor_point();
  monitor_point(double r, double z, const fields *f);
  ~monitor_point();
  double r, z, t;
  a_field h, e, *p;
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
  //
  // When fourier_transform is called with the other arguments as well,
  // harminv is run, with its output amplitudes and frequencies being
  // stored in a and f.
  void fourier_transform(component w,
                         complex<double> **a, complex<double> **f, int *numout,
                         double fmin=0.0, double fmax=0.0, int maxbands=100);
};

class fields {
 public:
  double *(hr[2]), *(hp[2]), *(hz[2]), *(er[2]), *(ep[2]), *(ez[2]);
  double *(backup_hr[2]), *(backup_hp[2]), *(backup_hz[2]);
  double *(hrp[2]), *(hpz[2]), *(hzr[2]), *(erp[2]), *(epz[2]), *(ezr[2]);
  double *(z_hrp[2][2]), *(z_hpz[2][2]), *(z_erp[2][2]), *(z_epz[2][2]);
  polarization *pol, *olpol;
  double a, inva; // The "lattice constant" and its inverse!
  int nr, nz;
  int npmlr, npmlz; // Amount of pml
  int m, t, phasein_time;
  double k, cosknz, sinknz;
  complex<double> eiknz;
  bandsdata *bands;
  src *e_sources, *h_sources;
  const mat *new_ma;
  mat *ma;
  const char *outdir;
  double preferred_fmax;

  fields(const mat *, int m);
  void use_bloch(double kz);
  ~fields();

  void output_slices(const char *name = "");
  void output_real_imaginary_slices(const char *name = "");
  void step();
  inline double time() { return t*inva*c; };

  void use_real_sources();
  int find_last_source();
  // Note that the following plane source only works if m == 1.
  void add_plane_source(double freq, double width, double peaktime,
                        double cutoff, double z, complex<double> amp(double r));
  void add_er_source(double freq, double width, double peaktime,
                     double cutoff, double z, complex<double> amp(double r));
  void add_ep_source(double freq, double width, double peaktime,
                     double cutoff, double z, complex<double> amp(double r));
  void add_ez_source(double freq, double width, double peaktime,
                     double cutoff, double z, complex<double> amp(double r));
  void add_hr_source(double freq, double width, double peaktime,
                     double cutoff, double z, complex<double> amp(double r));
  void add_hp_source(double freq, double width, double peaktime,
                     double cutoff, double z, complex<double> amp(double r));
  void add_hz_source(double freq, double width, double peaktime,
                     double cutoff, double z, complex<double> amp(double r));
  void initialize_with_nth_te(int n);
  void initialize_with_nth_tm(int n);
  void initialize_with_n_te(int n);
  void initialize_with_n_tm(int n);
  void initialize_polarizations(polarization *op=NULL, polarization *np=NULL);
  int phase_in_material(const mat *ma, double time);
  int is_phasing();

  void get_point(monitor_point *p, double r, double z);
  monitor_point *get_new_point(double r, double z, monitor_point *p=NULL);
  void output_point(FILE *, double r, double z, const char *name);

  flux_plane create_rflux_plane(double zmin, double zmax, double rconst);
  flux_plane create_zflux_plane(double rmin, double rmax, double zconst);
  complex<double> get_flux(flux_plane *fp);
  
  void prepare_for_bands(int z, int ttot, double fmax=0,
                         double qmin=1e300, double frac_pow_min=1e-7);
  void record_bands();
  complex<double> get_band(int n, int maxbands=100);
  void grace_bands(grace *, int maxbands=100);
  void output_bands(FILE *, const char *, int maxbands=100);
  void output_bands_and_modes(FILE *, const char *, int maxbands=100);
  double energy_in_box(double rmin, double rmax, double zmin, double zmax);
  double electric_energy_in_box(double rmin, double rmax, double zmin, double zmax);
  double magnetic_energy_in_box(double rmin, double rmax, double zmin, double zmax);
  double total_energy() {return energy_in_box(0.0, nr*inva, 0.0, nz*inva);};
  double zflux(int ri, int ro, int z);
  double rflux(int zl, int zu, int r);
  void dft_flux();
  int add_zfluxplane(int ri, int ro, int z);
  int add_rfluxplane(int zl, int zu, int r);
  int set_frequency_range(double wl, double wu, double deltaw);
  void ttow(complex<double> field, double *retarget, double *imtarget, double time);
  void fluxw_output(FILE *outpf, char *header);
  void set_output_directory(const char *name);
  void verbose(int v=1) { verbosity = v; }
 private: 
  int verbosity; // Turn on verbosity for debugging purposes...
  double *(erw[2]), *(epw[2]), *(ezw[2]), *(hrw[2]), *(hpw[2]), *(hzw[2]);
  int iposmax, ifreqmax, nfreq, nzflux, *(nzfluxplane[MAXFLUXPLANES]);
  int nrflux, *(nrfluxplane[MAXFLUXPLANES]);
  double *freqs;
  void find_source_z_position(double z, double shift, int *z1, int *z2,
                              complex<double> *amp1, complex<double> *amp2);
  void phase_material();
  void step_h_bulk();
  void step_h_pml();
  void step_h_boundaries();
  void step_h_source(const src *);
  void step_e_bulk();
  void step_e_pml();
  void step_e_boundaries();
  void step_polarization_itself(polarization *old = NULL, polarization *newp = NULL);
  void step_polarization_pml(polarization *old = NULL, polarization *newp = NULL);
  void step_e_polarization(polarization *old = NULL, polarization *newp = NULL);
  void step_e_pml_polarization(polarization *old = NULL, polarization *newp = NULL);
  void step_e_source(const src *);
  void add_src_pt(int r, int z,
                  complex<double> Pr, complex<double> Pp, complex<double> Pz,
                  double freq, double width, double peaktime,
                  double cutoff, int is_h = 0);
  int setifreqmax_and_iposmax(int ifreq, int ipos);
  void out_bands(FILE *, const char *, int maxbands, int outmodes);
  complex<double> *get_the_bands(int maxbands, double *approx_power = NULL);
};

class grace_point;
enum grace_type { XY, ERROR_BARS };

class grace {
 public:
  grace(const char *fname, const char *dirname = "");
  ~grace();
  
  void new_set(grace_type t = XY);
  void new_curve();
  void set_legend(const char *);
  void set_range(double xmin, double xmax, double ymin, double ymax);
  void output_point(double x, double y, double dy = -1.0);
  void output_out_of_order(int n, double x, double y, double dy = -1.0);
 private:
  void flush_pts();
  FILE *f;
  const char *fn, *dn;
  grace_point *pts;
  int set_num,sn;
};

// The following is a utility function to parse the executable name use it
// to come up with a directory name, avoiding overwriting any existing
// directory, unless the source file hasn't changed.

const char *make_output_directory(const char *exename, const char *jobname = NULL);
FILE *create_output_file(const char *dirname, const char *fname);

// The following allows you to hit ctrl-C to tell your calculation to stop
// and clean up.
void deal_with_ctrl_c(int stop_now = 2);
// When a ctrl_c is called, the following variable (which starts with a
// zero value) is incremented.
extern int interrupt;

int do_harminv(complex<double> *data, int n, int sampling_rate, double a,
	       double fmin, double fmax, int maxbands,
	       complex<double> *amps, double *freq_re, double *freq_im, double *errors = NULL);
