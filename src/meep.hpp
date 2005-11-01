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
#ifndef MEEP_H
#define MEEP_H

#include <stdio.h>
#include <math.h>

#include "meep/vec.hpp"
#include "meep/mympi.hpp"

namespace meep {

extern bool quiet; // if true, suppress all non-error messages from Meep

const double pi = 3.141592653589793238462643383276;

#ifdef INFINITY
const double infinity = INFINITY;
#else
const double infinity = 1e20; // this should be big enough for us
#endif

#ifdef NAN
const double nan = NAN;
#else
const double nan = -1.23454321e123; // ideally, a value never encountered in practice
#endif

class polarizability_identifier {
 public:
  double gamma, omeganot;
  double energy_saturation, saturated_sigma;
  bool operator==(const polarizability_identifier &);
};
class polarizability;
class polarization;
class grace;

// h5file.cpp: HDF5 file I/O.  Most users, if they use this
// class at all, will only use the constructor to open the file, and
// will otherwise use the fields::output_hdf5 functions.
class h5file {
public:
  typedef enum {
    READONLY, READWRITE, WRITE
  } access_mode;
  
  h5file(const char *filename_, access_mode m=READWRITE, bool parallel_=true);
  ~h5file(); // closes the files (and any open dataset)
  
  bool ok();
  
  double *read(const char *dataname, int *rank, int *dims, int maxrank);
  void write(const char *dataname, int rank, const int *dims, double *data,
	     bool single_precision = true);
  
  char *read(const char *dataname);
  void write(const char *dataname, const char *data);
  
  void create_data(const char *dataname, int rank, const int *dims,
		   bool append_data = false,
		   bool single_precision = true);
  void extend_data(const char *dataname, int rank, const int *dims);
  void create_or_extend_data(const char *dataname, int rank,
			     const int *dims,
			     bool append_data, bool single_precision);
  void write_chunk(int rank, const int *chunk_start, const int *chunk_dims,
		   double *data);
  void done_writing_chunks();
  
  void read_size(const char *dataname, int *rank, int *dims, int maxrank);
  void read_chunk(int rank, const int *chunk_start, const int *chunk_dims,
		  double *data);
  
  void remove();
  void remove_data(const char *dataname);
  
  const char *file_name() const { return filename; }

  void prevent_deadlock(); // hackery for exclusive mode
private:
  access_mode mode;
  char *filename;
  bool parallel;

  bool is_cur(const char *dataname);
  void unset_cur();
  void set_cur(const char *dataname, void *data_id);
  char *cur_dataname;
  bool cur_append_data;
  
  /* linked list to keep track of which datasets we are extending...
     this is necessary so that create_or_extend_data can know whether
     to create (overwrite) a dataset or extend it. */
  struct extending_s {
    int dindex;
    char *dataname;
    struct extending_s *next;
  } *extending;
  extending_s *get_extending(const char *dataname) const;
  
  /* store hid_t values as hid_t* cast to void*, so that
     files including meep.h don't need hdf5.h */
  void *id; /* file */
  void *cur_id; /* dataset, if any */

  void *get_id(); // get current (file) id, opening/creating file if needed
  void close_id();
};

/* This class is used to compute position-dependent material properties
   like the dielectric function, polarizability sigma, 
   nonlinearities, et cetera.  Simple cases of stateless functions are
   handled by canned subclasses below, but more complicated cases
   can be handled by creating a user-defined subclass of material_function.
   It is useful to group different properties into one class because
   it is likely that complicated implementations will share state between
   properties. */
class material_function {
     material_function(const material_function &ef) {} // prevent copying
 public:
     material_function() {}

     virtual ~material_function() {}

     /* Specify a restricted volume: all subsequent eps/sigma/etc
	calls will be for points inside gv, until the next set_volume. */
     virtual void set_volume(const geometric_volume &gv) {}
     virtual void unset_volume(void) {} // unrestrict the volume

     /* scalar dielectric function */
     virtual double eps(const vec &r) { return 1.0; }

     /* polarizability sigma function */
     virtual double sigma(const vec &r) { return 0.0; }
     /* specify polarizability used for subsequent calls to sigma(r) */
     virtual void set_polarizability(double omega, double gamma,
				     double delta_epsilon,
				     double energy_saturation) {}

     // Kerr coefficient
     
     virtual bool has_kerr() { return false; }
     virtual double kerr(const vec &r) { return 0.0; }

     // TODO: dielectric tensor, ...
};

class simple_material_function : public material_function {
     double (*f)(const vec &);
     
 public:
     simple_material_function(double (*func)(const vec &)) { f = func; }

     virtual ~simple_material_function() {}

     virtual double eps(const vec &r) { return f(r); }
     virtual double sigma(const vec &r) { return f(r); }
     virtual double kerr(const vec &r) { return f(r); }
};

class structure;

class structure_chunk {
 public:
  double a, Courant, dt; // res. a, Courant num., and timestep dt=Courant/a
  double *eps, *kerr[NUM_FIELD_COMPONENTS];
  double *inveps[NUM_FIELD_COMPONENTS][5];
  double *C[5][NUM_FIELD_COMPONENTS];
  double *Cdecay[5][NUM_FIELD_COMPONENTS][5];
  volume v;  // integer volume that could be bigger than non-overlapping gv below
  geometric_volume gv;
  polarizability *pb;

  ~structure_chunk();
  structure_chunk(const volume &v,
            const geometric_volume &vol_limit, double Courant, int proc_num);
  structure_chunk(const structure_chunk *);
  void set_epsilon(material_function &eps,
                   bool use_anisotropic_averaging, double minvol);
  void set_kerr(material_function &eps);
  void use_pml(direction, double dx, double boundary_loc, double strength);
  void update_pml_arrays();
  void update_Cdecay();

  void add_polarizability(double sigma(const vec &), double omega, double gamma,
                          double delta_epsilon = 1.0, double energy_saturation = 0.0);
  void add_polarizability(material_function &sigma, double omega, double gamma,
                          double delta_epsilon = 1.0, double energy_saturation = 0.0);

  void mix_with(const structure_chunk *, double);

  int n_proc() const { return the_proc; } // Says which proc owns me!
  int is_mine() const { return the_is_mine; }

  // monitor.cpp
  double get_inveps(component, direction, const ivec &iloc) const;
  double max_eps() const;
 private:
  double pml_fmin;
  int the_proc;
  int the_is_mine;
};

// linked list of descriptors for boundary regions (currently just for PML)
class boundary_region {
public:
  typedef enum { NOTHING_SPECIAL, PML } boundary_region_kind;
  
  boundary_region() :
    kind(NOTHING_SPECIAL), thickness(0.0), strength(1.0), d(NO_DIRECTION), side(Low), next(0) {}
  boundary_region(boundary_region_kind kind, double thickness, double strength, direction d, boundary_side side, boundary_region *next = 0) : kind(kind), thickness(thickness), strength(strength), d(d), side(side), next(next) {}

  boundary_region(const boundary_region &r) : kind(r.kind), thickness(r.thickness), strength(r.strength), d(r.d), side(r.side) { 
    next = r.next ? new boundary_region(*r.next) : 0;
  }

  ~boundary_region() { if (next) delete next; }
  
  void operator=(const boundary_region &r) {
    kind = r.kind; thickness = r.thickness; strength = r.strength;
    d = r.d; side = r.side;
    if (next) delete next;
    next = r.next ? new boundary_region(*r.next) : 0;
  }
  boundary_region operator+(const boundary_region &r0) const {
    boundary_region r(*this), *cur = &r;
    while (cur->next) cur = cur->next;
    cur->next = new boundary_region(r0);
    return r;
  }

  boundary_region operator*(double strength_mult) const {
    boundary_region r(*this), *cur = &r;
    while (cur) {
      cur->strength *= strength_mult;
      cur = cur->next;
    }
    return r;
  }

  void apply(structure *s) const;
  void apply(const structure *s, structure_chunk *sc) const;

private:
  boundary_region_kind kind;
  double thickness, strength;
  direction d;
  boundary_side side;
  boundary_region *next;
};

boundary_region pml(double thickness, direction d, boundary_side side);
boundary_region pml(double thickness, direction d);
boundary_region pml(double thickness);
#define no_pml() boundary_region()

class structure {
 public:
  structure_chunk **chunks;
  int num_chunks;
  volume v, user_volume;
  double a, Courant, dt; // res. a, Courant num., and timestep dt=Courant/a
  geometric_volume gv;
  symmetry S;
  const char *outdir;
  volume *effort_volumes;
  double *effort;
  int num_effort_volumes;

  ~structure();
  structure();
  structure(const volume &v, material_function &eps,
	    const boundary_region &br = boundary_region(),
	    const symmetry &s = meep::identity(),
	    int num_chunks = 0, double Courant = 0.5);
  structure(const volume &v, double eps(const vec &), 
	    const boundary_region &br = boundary_region(),
	    const symmetry &s = meep::identity(),
	    int num_chunks = 0, double Courant = 0.5);
  structure(const structure *);
  structure(const structure &);

  void set_materials(material_function &mat,
		     bool use_anisotropic_averaging=true, double minvol=0.0);
  void set_epsilon(material_function &eps,
                   bool use_anisotropic_averaging=true, double minvol=0.0);
  void set_epsilon(double eps(const vec &),
                   bool use_anisotropic_averaging=true, double minvol=0.0);
  void set_kerr(material_function &eps);
  void set_kerr(double eps(const vec &));
  polarizability_identifier
     add_polarizability(double sigma(const vec &), double omega, double gamma,
                  double delta_epsilon = 1.0, double energy_saturation = 0.0);
  polarizability_identifier
     add_polarizability(material_function &sigma, double omega, double gamma,
                  double delta_epsilon = 1.0, double energy_saturation = 0.0);

  void output_slices(const char *name="") const;
  void output_slices(const geometric_volume &what, const char *name="") const;

  void set_output_directory(const char *name);
  void mix_with(const structure *, double);

  bool equal_layout(const structure &) const;

  // monitor.cpp
  double get_inveps(component, direction, const ivec &origloc) const;
  double get_inveps(component, direction, const vec &loc) const;
  double get_eps(const vec &loc) const;
  double max_eps() const;

  friend class boundary_region;

 private:
  void use_pml(direction d, boundary_side b, double dx, double strength);
  void add_to_effort_volumes(const volume &new_effort_volume, 
			     double extra_effort);
  void choose_chunkdivision(const volume &v, int num_chunks,
			    const boundary_region &br, const symmetry &s);
  void check_chunks();
};

class src_vol;
class bandsdata;
class fields;
class fields_chunk;
class flux_vol;

// Time-dependence of a current source, intended to be overridden by
// subclasses.  current() and dipole() are be related by
// current = d(dipole)/dt (or rather, the finite-difference equivalent).
class src_time {
 public:
  src_time() { current_time = nan; current_current = 0.0; next = NULL; }
  virtual ~src_time() { delete next; }
  src_time(const src_time &t) { 
       current_time = t.current_time;
       current_current = t.current_current;
       current_dipole = t.current_dipole;
       if (t.next) next = t.next->clone(); else next = NULL;
  }
  
  complex<double> dipole() const { return current_dipole; }
  complex<double> current() const { return current_current; }
  void update(double time, double dt) {
    if (time != current_time) {
      current_dipole = dipole(time);
      current_current = current(time, dt);
      current_time = time;
    }
  }

  complex<double> current(double time, double dt) const { 
    return ((dipole(time + dt) - dipole(time)) / dt);
  }

  double last_time_max() { return last_time_max(0.0); }
  double last_time_max(double after);
  
  src_time *add_to(src_time *others, src_time **added) const;
  src_time *next;

  // subclasses should override these methods:
  virtual complex<double> dipole(double time) const { (void)time; return 0; }
  virtual double last_time() const { return 0.0; }
  virtual src_time *clone() const { return new src_time(*this); }
  virtual bool is_equal(const src_time &t) const { (void)t; return 1; }
  virtual complex<double> frequency() const { return 0.0; }

 private:
  double current_time;
  complex<double> current_dipole, current_current;
};

bool src_times_equal(const src_time &t1, const src_time &t2);

// Gaussian-envelope source with given frequency, width, peak-time, cutoff
class gaussian_src_time : public src_time {
 public:
  gaussian_src_time(double f, double fwidth, double s = 5.0);
  gaussian_src_time(double f, double w, double start_time, double end_time);
  virtual ~gaussian_src_time() {}

  virtual complex<double> dipole(double time) const;
  virtual double last_time() const { return peak_time + cutoff; };
  virtual src_time *clone() const { return new gaussian_src_time(*this); }
  virtual bool is_equal(const src_time &t) const;
  virtual complex<double> frequency() const { return freq; }

 private:
  double freq, width, peak_time, cutoff;
};

// Continuous (CW) source with (optional) slow turn-on and/or turn-off.
class continuous_src_time : public src_time {
 public:
  continuous_src_time(complex<double> f, double w = 0.0, 
		      double st = 0.0, double et = infinity,
		      double s = 3.0) : freq(f), width(w), start_time(st),
					end_time(et), slowness(s) {}
  virtual ~continuous_src_time() {}
  
  virtual complex<double> dipole(double time) const;
  virtual double last_time() const { return end_time; };
  virtual src_time *clone() const { return new continuous_src_time(*this); }
  virtual bool is_equal(const src_time &t) const;
  virtual complex<double> frequency() const { return freq; }
  
 private:
  complex<double> freq;
  double width, start_time, end_time, slowness;
};

// user-specified source function with start and end times
class custom_src_time : public src_time {
 public:
  custom_src_time(complex<double> (*func)(double t, void *), void *data,
		      double st = -infinity, double et = infinity)
    : func(func), data(data), start_time(st), end_time(et) {}
  virtual ~custom_src_time() {}
  
  virtual complex<double> dipole(double time) const { return func(time,data); }
  virtual double last_time() const { return end_time; };
  virtual src_time *clone() const { return new custom_src_time(*this); }
  virtual bool is_equal(const src_time &t) const;
  
 private:
  complex<double> (*func)(double t, void *);
  void *data;
  double start_time, end_time;
};

class monitor_point {
 public:
  monitor_point();
  ~monitor_point();
  vec loc;
  double t;
  complex<double> f[NUM_FIELD_COMPONENTS];
  monitor_point *next;

  complex<double> get_component(component);
  double poynting_in_direction(direction d);
  double poynting_in_direction(vec direction_v);

  // When called with only its first four arguments, fourier_transform
  // performs an FFT on its monitor points, putting the frequencies in f
  // and the amplitudes in a.  Yes, the frequencies are trivial and
  // redundant, but this saves you the risk of making a mistake in
  // converting your units.  Note also, that in this case f is always a
  // real number, although it's stored in a complex.
  //
  // Note that in either case, fourier_transform assumes that the monitor
  // points are all equally spaced in time.
  void fourier_transform(component w,
                         complex<double> **a, complex<double> **f, int *numout,
                         double fmin=0.0, double fmax=0.0, int maxbands=100);
  // harminv works much like fourier_transform, except that it is not yet
  // implemented.
  void harminv(component w,
               complex<double> **a, complex<double> **f,
               int *numout, double fmin, double fmax,
               int maxbands);
};

// dft.cpp
// this should normally only be created with fields::add_dft
class dft_chunk {
public:
  dft_chunk(fields_chunk *fc_,
	    ivec is_, ivec ie_,
	    vec s0_, vec s1_, vec e0_, vec e1_,
	    double dV0_, double dV1_,
	    complex<double> scale_,
	    component c_,
	    const void *data_);
  ~dft_chunk();
  
  void update_dft(double time);

  void scale_dft(complex<double> scale);

  void operator-=(const dft_chunk &chunk);

  // the frequencies to loop_in_chunks
  double omega_min, domega;
  int Nomega;

  component c; // component to DFT (possibly transformed by symmetry)

  int N; // number of spatial points (on epsilon grid)
  complex<double> *dft; // N x Nomega array of DFT values.

  struct dft_chunk *next_in_chunk; // per-fields_chunk list of DFT chunks
  struct dft_chunk *next_in_dft; // next for this particular DFT vol./component

private:
  // parameters passed from field_integrate:
  fields_chunk *fc;
  ivec is, ie;
  vec s0, s1, e0, e1;
  double dV0, dV1;
  complex<double> scale; // scale factor * phase from shift and symmetry

  // cache of exp(iwt) * scale, of length Nomega
  complex<double> *dft_phase;

  int avg1, avg2; // index offsets for average to get epsilon grid
};

void save_dft_hdf5(dft_chunk *dft_chunks, component c, h5file *file,
		   const char *dprefix = 0);
void load_dft_hdf5(dft_chunk *dft_chunks, component c, h5file *file,
		   const char *dprefix = 0);

// dft.cpp (normally created with fields::add_dft_flux)
class dft_flux {
public:
  dft_flux(const component cE_, const component cH_,
	   dft_chunk *E_, dft_chunk *H_,
	   double fmin, double fmax, int Nf);
  dft_flux(const dft_flux &f);

  double *flux();

  void save_hdf5(h5file *file, const char *dprefix = 0);
  void load_hdf5(h5file *file, const char *dprefix = 0);

  void operator-=(const dft_flux &fl) { if (E && fl.E) *E -= *fl.E; if (H && fl.H) *H -= *fl.H; }

  void save_hdf5(fields &f, const char *fname, const char *dprefix = 0);
  void load_hdf5(fields &f, const char *fname, const char *dprefix = 0);

  void scale_dfts(complex<double> scale);

  void remove() { delete E; delete H; E = H = 0; }

  double freq_min, dfreq;
  int Nfreq;
  dft_chunk *E, *H;
  component cE, cH;
};

enum in_or_out { Incoming=0, Outgoing };
enum connect_phase { CONNECT_PHASE = 0, CONNECT_NEGATE=1, CONNECT_COPY=2 };

class fields_chunk {
 public:
  double *f[NUM_FIELD_COMPONENTS][2];
  double *f_backup[NUM_FIELD_COMPONENTS][2];
  double *f_p_pml[NUM_FIELD_COMPONENTS][2];
  double *f_m_pml[NUM_FIELD_COMPONENTS][2];
  double *f_backup_p_pml[NUM_FIELD_COMPONENTS][2];
  double *f_backup_m_pml[NUM_FIELD_COMPONENTS][2];

  bool have_d_minus_p;
  double *d_minus_p[NUM_FIELD_COMPONENTS][2];

  dft_chunk *dft_chunks;

  double **zeroes[NUM_FIELD_TYPES]; // Holds pointers to metal points.
  int num_zeroes[NUM_FIELD_TYPES];
  double **connections[NUM_FIELD_TYPES][CONNECT_COPY+1][Outgoing+1];
  int num_connections[NUM_FIELD_TYPES][CONNECT_COPY+1][Outgoing+1];
  complex<double> *connection_phases[NUM_FIELD_TYPES];

  polarization *pol, *olpol;
  double a, Courant, dt; // res. a, Courant num., and timestep dt=Courant/a
  volume v;
  geometric_volume gv;
  double m, rshift;
  int is_real;
  bandsdata *bands;
  src_vol *e_sources, *h_sources;
  const structure_chunk *new_s;
  structure_chunk *s;
  const char *outdir;

  fields_chunk(const structure_chunk *, const char *outdir, double m);
  fields_chunk(const fields_chunk &);
  ~fields_chunk();

  // step.cpp
  double peek_field(component, const vec &);

  void use_real_fields();
  bool have_component(component c, bool is_complex = false) {
    switch (c) {
    case Dielectric:
      return !is_complex;
    default:
      return (f[c][0] && f[c][is_complex]);
    }
  }

  double last_source_time();
  // monitor.cpp
  complex<double> get_field(component, const ivec &) const;

  // for non-collective interpolation:
  geometric_volume get_field_gv(component) const;
  complex<double> get_field(component, const vec &) const;

  complex<double> get_polarization_field(const polarizability_identifier &p,
                                         component c, const ivec &iloc) const;
  double get_polarization_energy(const ivec &) const;
  double my_polarization_energy(const ivec &) const;
  double get_polarization_energy(const polarizability_identifier &, const ivec &) const;
  double my_polarization_energy(const polarizability_identifier &, const ivec &) const;
  double get_inveps(component, direction, const ivec &iloc) const;
  double get_eps(const ivec &iloc) const;
  complex<double> analytic_epsilon(double freq, const vec &) const;
  
  // slices.cpp
  double maxfieldmag(component) const;
  double maxpolenergy() const;
  double minpolenergy() const;
  void output_eps_body(component, const symmetry &, int sn,
                       const geometric_volume &what, file *,
                       complex<double> phase_shift);
  void output_eps_body(const polarizability_identifier &p,
                       component c, const symmetry &S, int sn,
                       const geometric_volume &what, file *out,
                       complex<double> phshift);
  complex<double> field_mean(component c, bool abs_real, bool abs_imag) const;

  void backup_h();
  void restore_h();

  void set_output_directory(const char *name);
  void verbose(int v=1) { verbosity = v; }

  double count_volume(component);
  friend class fields;

  int n_proc() const { return s->n_proc(); };
  int is_mine() const { return s->is_mine(); };
  // boundaries.cpp
  void zero_metal(field_type);
  // polarization.cpp
  void initialize_polarization_energy(const polarizability_identifier &,
                                      double energy(const vec &));

  // fields.cpp
  void alloc_f(component c);
  void remove_sources();
  void zero_fields();

 private: 
  int verbosity; // Turn on verbosity for debugging purposes...
  // fields.cpp
  void figure_out_step_plan();
  bool have_plus_deriv[NUM_FIELD_COMPONENTS], have_minus_deriv[NUM_FIELD_COMPONENTS];
  component plus_component[NUM_FIELD_COMPONENTS], minus_component[NUM_FIELD_COMPONENTS];
  direction plus_deriv_direction[NUM_FIELD_COMPONENTS],
            minus_deriv_direction[NUM_FIELD_COMPONENTS];
  int num_each_direction[3], stride_each_direction[3];
  int num_any_direction[5], stride_any_direction[5];
  // bands.cpp
  void record_bands(int tcount);
  // step.cpp
  void phase_in_material(const structure_chunk *s);
  void phase_material(int phasein_time);
  void step_h();
  void step_h_source(src_vol *);
  void step_d();
  void step_d_source(src_vol *);
  void update_e_from_d_prepare(void);
  void update_e_from_d_sources(void);
  void update_e_from_d_update(void);
  void update_e_from_d();
  void update_from_e();
  void calc_sources(double time);

  // initialize.cpp
  void initialize_field(component, complex<double> f(const vec &));
  void initialize_polarizations(polarization *op=NULL, polarization *np=NULL);
  void initialize_with_nth_te(int n, double kz);
  void initialize_with_nth_tm(int n, double kz);
  // boundaries.cpp
  void alloc_extra_connections(field_type, connect_phase, in_or_out, int);
  // dft.cpp
  void update_dfts(double timeE, double timeH);
};

enum boundary_condition { Periodic=0, Metallic, Magnetic, None };
enum time_sink { Connecting, Stepping, Boundaries, MpiTime,
                 FieldOutput, Other };

typedef void (*field_chunkloop)(fields_chunk *fc, int ichunk, component cgrid,
				ivec is, ivec ie,
				vec s0, vec s1, vec e0, vec e1,
				double dV0, double dV1,
				ivec shift, complex<double> shift_phase, 
				const symmetry &S, int sn,
				void *chunkloop_data);
typedef complex<double> (*field_function)(const complex<double> *fields,
					   const vec &loc,
					   void *integrand_data_);
typedef double (*field_rfunction)(const complex<double> *fields,
				   const vec &loc,
				   void *integrand_data_);

class fields {
 public:
  int num_chunks;
  fields_chunk **chunks;
  src_time *sources;
  flux_vol *fluxes;
  symmetry S;

  // The following is an array that is num_chunks by num_chunks.  Actually
  // it is two arrays, one for the imaginary and one for the real part.
  double **comm_blocks[NUM_FIELD_TYPES];
  // This is the same size as each comm_blocks array, and store the sizes
  // of the comm blocks themselves for each connection-phase type
  int *comm_sizes[NUM_FIELD_TYPES][CONNECT_COPY+1];
  int comm_size_tot(int f, int pair) const {
    int sum = 0; for (int ip=0; ip<3; ++ip) sum+=comm_sizes[f][ip][pair];
    return sum;
  }

  double a, dt; // The resolution a and timestep dt=Courant/a
  volume v, user_volume;
  geometric_volume gv;
  double m;
  int t, phasein_time, is_real;
  complex<double> k[5], eikna[5];
  double coskna[5], sinkna[5];
  boundary_condition boundaries[2][5];
  bandsdata *bands;
  char *outdir;
  // fields.cpp methods:
  fields(const structure *, double m=0);
  fields(const fields &);
  ~fields();
  bool equal_layout(const fields &f) const;
  void use_real_fields();
  void zero_fields();
  void remove_sources();
  void remove_fluxes();
  void reset();
  bool disable_sources; // set to true to turn off sources (w/o deleting)
  // time.cpp
  double time_spent_on(time_sink);
  void print_times();
  // boundaries.cpp
  void set_boundary(boundary_side,direction,
                    boundary_condition, complex<double> kcomponent=0.0);
  void use_bloch(direction d, double k) { use_bloch(d, (complex<double>) k); }
  void use_bloch(direction, complex<double> kz);
  void use_bloch(const vec &k);
  vec lattice_vector(direction) const;
  // slices.cpp methods:
  void output_slices(const char *name = "");
  void output_slices(const geometric_volume &what, const char *name = "");
  void eps_envelope(const char *name = "");
  void eps_envelope(const geometric_volume &what, const char *name = "");
  void eps_slices(const char *name = "");
  void eps_slices(const vec &origin, const vec &xside, const vec &yside,
                  const double dx = 0.05, const char *name = "");
  void eps_slices(const geometric_volume &what, const char *name = "");
  void eps_polarization_slice(const polarizability_identifier &, const char *name = "");
  void eps_polarization_slice(const polarizability_identifier &,
                              const geometric_volume &, const char *name = "");
  void eps_energy_slice(const polarizability_identifier &, const char *name = "");
  void eps_energy_slice(const polarizability_identifier &,
                        const geometric_volume &what, const char *name = "");
  void eps_energy_slice(const char *name = "");
  void eps_energy_slice(const geometric_volume &what, const char *name = "");
  void output_real_imaginary_slices(const char *name = "");
  void output_real_imaginary_slices(const geometric_volume &what,
                                    const char *name = "");

  // h5fields.cpp:
  // low-level function:
  void output_hdf5(h5file *file, const char *dataname,
		   int num_fields, const component *components,
		   field_function fun, int reim,
		   const geometric_volume &where,
		   void *fun_data_ = 0,
		   bool append_data = false,
		   bool single_precision = false);
  // higher-level functions
  void output_hdf5(const char *dataname,  // OUTPUT COMPLEX-VALUED FUNCTION
		   int num_fields, const component *components,
		   field_function fun,
		   const geometric_volume &where,
		   void *fun_data_ = 0,
		   h5file *file = 0,
		   bool append_data = false,
		   bool single_precision = false,
		   const char *prefix = 0);
  void output_hdf5(const char *dataname,  // OUTPUT REAL-VALUED FUNCTION
		   int num_fields, const component *components,
		   field_rfunction fun,
		   const geometric_volume &where,
		   void *fun_data_ = 0,
		   h5file *file = 0,
		   bool append_data = false,
		   bool single_precision = false,
		   const char *prefix = 0);
  void output_hdf5(component c,   // OUTPUT FIELD COMPONENT (or Dielectric)
		   const geometric_volume &where,
		   h5file *file = 0,
		   bool append_data = false,
		   bool single_precision = false,
		   const char *prefix = 0);
  void output_hdf5(derived_component c,   // OUTPUT DERIVED FIELD COMPONENT
		   const geometric_volume &where,
		   h5file *file = 0,
		   bool append_data = false,
		   bool single_precision = false,
		   const char *prefix = 0);
  h5file *open_h5file(const char *name, 
		      h5file::access_mode mode = h5file::WRITE,
		      const char *prefix = NULL, bool timestamp = false);
  const char *h5file_name(const char *name,
			  const char *prefix = NULL, bool timestamp = false);

  double maxfieldmag_to_master(component) const;
  double minpolenergy_to_master() const;
  double maxpolenergy_to_master() const;
  complex<double> optimal_phase_shift(component) const;
  // step.cpp methods:
  double last_step_output_wall_time;
  int last_step_output_t;
  void step();
  inline double time() const { return t*dt; };

  // cw_fields.cpp:
  bool solve_cw(double tol, int maxiters, complex<double> frequency, int L=2);
  bool solve_cw(double tol = 1e-8, int maxiters = 10000, int L=2);

  // sources.cpp:
  double last_source_time();
  void add_point_source(component c, double freq, double width, double peaktime,
                        double cutoff, const vec &, complex<double> amp = 1.0,
                        int is_continuous = 0);
  void add_point_source(component c, const src_time &src,
                        const vec &, complex<double> amp = 1.0);
  void add_volume_source(component c, const src_time &src,
			 const geometric_volume &, 
			 complex<double> A(const vec &),
			 complex<double> amp = 1.0);
  void add_volume_source(component c, const src_time &src,
			 const geometric_volume &, 
			 complex<double> amp = 1.0);
  void require_component(component c);

  // initialize.cpp:
  void initialize_field(component, complex<double> f(const vec &));
  void initialize_A(complex<double> A(component, const vec &), double freq);
  void initialize_with_nth_te(int n);
  void initialize_with_nth_tm(int n);
  void initialize_with_n_te(int n);
  void initialize_with_n_tm(int n);
  void initialize_polarizations();
  int phase_in_material(const structure *s, double time);
  int is_phasing();

  // loop_in_chunks.cpp
  void loop_in_chunks(field_chunkloop chunkloop, void *chunkloop_data,
		      const geometric_volume &where,
		      component cgrid = Dielectric,
		      bool use_symmetry = true,
		      bool snap_unit_dims = false);
  
  // integrate.cpp
  complex<double> integrate(int num_fields, const component *components,
			    field_function integrand,
			    const geometric_volume &where,
			    void *integrand_data_ = 0,
			    double *maxabs = 0);
  double max_abs(int num_fields, const component *components,
		 field_function integrand,
		 const geometric_volume &where,
		 void *integrand_data = 0);
  double max_abs(component c, const geometric_volume &where);
  
  // dft.cpp
  dft_chunk *add_dft(component c, const geometric_volume &where,
		     double freq_min, double freq_max, int Nfreq,
		     bool include_dV = true,
		     complex<double> weight = 1.0, dft_chunk *chunk_next = 0);
  dft_chunk *add_dft_pt(component c, const vec &where,
			double freq_min, double freq_max, int Nfreq);
  dft_chunk *add_dft(const geometric_volume_list *where,
		     double freq_min, double freq_max, int Nfreq,
		     bool include_dV = true);
  void update_dfts();
  dft_flux add_dft_flux(direction d, const geometric_volume &where,
			double freq_min, double freq_max, int Nfreq);
  dft_flux add_dft_flux_box(const geometric_volume &where,
			    double freq_min, double freq_max, int Nfreq);
  dft_flux add_dft_flux_plane(const geometric_volume &where,
			      double freq_min, double freq_max, int Nfreq);
  dft_flux add_dft_flux(const geometric_volume_list *where,
			double freq_min, double freq_max, int Nfreq);
  
  // monitor.cpp
  double get_inveps(component, direction, const vec &loc) const;
  double get_eps(const vec &loc) const;
  void get_point(monitor_point *p, const vec &) const;
  monitor_point *get_new_point(const vec &, monitor_point *p=NULL) const;
  complex<double> analytic_epsilon(double freq, const vec &) const;
  
  void prepare_for_bands(const vec &, double end_time, double fmax=0,
                         double qmin=1e300, double frac_pow_min=0.0);
  void record_bands();
  complex<double> get_band(int n, int maxbands=100);
  void grace_bands(grace *, int maxbands=100);
  void output_bands(file *, const char *, int maxbands=100);
  complex<double> get_field(int c, const vec &loc) const;
  complex<double> get_field(component c, const vec &loc) const;
  double get_field(derived_component c, const vec &loc) const;

  // energy_and_flux.cpp
  double energy_in_box(const geometric_volume &);
  double electric_energy_in_box(const geometric_volume &);
  double magnetic_energy_in_box(const geometric_volume &);
  double thermo_energy_in_box(const geometric_volume &);
  double total_energy();
  double field_energy_in_box(const geometric_volume &);
  double field_energy_in_box(component c, const geometric_volume &);
  double field_energy();
  double flux_in_box_wrongH(direction d, const geometric_volume &);
  double flux_in_box(direction d, const geometric_volume &);
  flux_vol *add_flux_vol(direction d, const geometric_volume &where);
  flux_vol *add_flux_plane(const geometric_volume &where);
  flux_vol *add_flux_plane(const vec &p1, const vec &p2);
  double electric_energy_max_in_box(const geometric_volume &where);
  double modal_volume_in_box(const geometric_volume &where);
  double electric_sqr_weighted_integral(double (*deps)(const vec &),
				       const geometric_volume &where);
  double electric_energy_weighted_integral(double (*f)(const vec &),
					   const geometric_volume &where);

  void set_output_directory(const char *name);
  void verbose(int v=1);
  double count_volume(component);
  // polarization.cpp
  void initialize_polarization_energy(const polarizability_identifier &,
                                      double energy(const vec &));

  // fields.cpp
  bool have_component(component);
  void set_rshift(double rshift);
  // material.cpp
  double max_eps() const;
  // step.cpp
  void force_consistency(field_type ft);

 private: 
  int verbosity; // Turn on verbosity for debugging purposes...
  double last_wall_time;
  time_sink working_on, was_working_on;
  double times_spent[Other+1];
  // time.cpp
  void am_now_working_on(time_sink);
  void finished_working();
  // boundaries.cpp
  bool chunk_connections_valid;
  void find_metals();
  void disconnect_chunks();
  void connect_chunks();
  void connect_the_chunks(); // Intended to be ultra-private...
  bool on_metal_boundary(const ivec &);
  ivec ilattice_vector(direction) const;
  bool locate_component_point(component *, ivec *, complex<double> *) const;
  bool locate_point_in_user_volume(ivec *, complex<double> *phase) const;
  void locate_volume_source_in_user_volume(const vec p1, const vec p2, vec newp1[8], vec newp2[8],
                                           complex<double> kphase[8], int &ncopies) const;
  // step.cpp
  void phase_material();
  void step_h();
  void step_h_source();
  void step_d();
  void step_d_source();
  void update_e_from_d();
  void update_from_e();
  void step_boundaries(field_type);
  void calc_sources(double tim);
  int cluster_some_bands_cleverly(double *tf, double *td, complex<double> *ta,
                                  int num_freqs, int fields_considered, int maxbands,
                                  complex<double> *fad, double *approx_power);
  void out_bands(file *, const char *, int maxbands);
  complex<double> *clever_cluster_bands(int maxbands, double *approx_power = NULL);
  // slices.cpp
  void outline_chunks(file *name);
  bool has_eps_interface(vec *loc) const;
  complex<double> field_mean(component c, bool abs_real = false,
                             bool abs_imag = false) const;
  // monitor.cpp
  complex<double> get_field(component c, const ivec &iloc) const;
  double get_polarization_energy(const ivec &) const;
  double get_polarization_energy(const vec &) const;
  double get_polarization_energy(const polarizability_identifier &, const ivec &) const;
  double get_polarization_energy(const polarizability_identifier &, const vec &) const;
  double get_inveps(component, direction, const ivec &iloc) const;
  double get_eps(const ivec &iloc) const;
};

class flux_vol {
 public:
  flux_vol(fields *f_, direction d_, const geometric_volume &where_) : where(where_) {
    f = f_; d = d_; cur_flux = cur_flux_half = 0; 
    next = f->fluxes; f->fluxes = this;
  }
  ~flux_vol() { delete next; }

  void update_half() { cur_flux_half = flux_wrongE(); 
                       if (next) next->update_half(); }
  void update() { cur_flux = (flux_wrongE() + cur_flux_half) * 0.5;
                  if (next) next->update(); }

  double flux() { return cur_flux; }

  flux_vol *next;
 private:
  double flux_wrongE() { return f->flux_in_box_wrongH(d, where); }
  fields *f;
  direction d;
  geometric_volume where;
  double cur_flux, cur_flux_half;
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

int do_harminv(complex<double> *data, int n, double dt,
	       double fmin, double fmax, int maxbands,
	       complex<double> *amps, double *freq_re, double *freq_im,
	       double *errors = NULL,
	       double spectral_density = 1.1, double Q_thresh = 50,
	       double rel_err_thresh = 1e20, double err_thresh = 0.01, 
	       double rel_amp_thresh = -1, double amp_thresh = -1);

} /* namespace meep */

#endif /* MEEP_H */
