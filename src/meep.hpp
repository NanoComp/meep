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
#ifndef MEEP_H
#define MEEP_H

#include <stdio.h>
#include <math.h>
#include <vector>


#include "meep/vec.hpp"
#include "meep/mympi.hpp"

namespace meep {

/* We use the type realnum for large arrays, e.g. the fields.
   For local variables and small arrays, we use double precision,
   but for things like the fields we can often get away with
   single precision (since the errors are not dominated by roundoff). 
   However, we will default to using double-precision for large 
   arrays, as the factor of two in memory and the moderate increase
   in speed currently don't seem worth the loss of precision. */
#define MEEP_SINGLE 0 // 1 for single precision, 0 for double
#if MEEP_SINGLE
typedef float realnum;
#else
typedef double realnum;
#endif

extern bool quiet; // if true, suppress all non-error messages from Meep

const double pi = 3.141592653589793238462643383276;

const double infinity = HUGE_VAL;

#ifdef NAN
const double nan = NAN;
#else
const double nan = -7.0415659787563146e103; // ideally, a value never encountered in practice
#endif

/* generic base class, only used by subclassing: represents susceptibility 
   polarizability vector P = chi(omega) W  (where W = E or H). */
class susceptibility {
public:
  susceptibility() { id = cur_id++; ntot = 0; next = NULL; 
    FOR_COMPONENTS(c) FOR_DIRECTIONS(d) {
      sigma[c][d] = NULL; trivial_sigma[c][d] = true; } }
  susceptibility(const susceptibility &s) { id = s.id; ntot = s.ntot; 
    next = NULL; FOR_COMPONENTS(c) FOR_DIRECTIONS(d) {
      sigma[c][d] = NULL; trivial_sigma[c][d] = true; } }
  virtual susceptibility *clone() const;
  virtual ~susceptibility() { 
    FOR_COMPONENTS(c) FOR_DIRECTIONS(d) delete[] sigma[c][d];
    delete next; }

  int get_id() const { return id; }
  bool operator==(const susceptibility &s) const { return id == s.id; };

  // update all of the internal polarization state given the W field
  // at the current time step, possibly the previous field W_prev, etc.
  virtual void update_P(realnum *W[NUM_FIELD_COMPONENTS][2], 
			realnum *W_prev[NUM_FIELD_COMPONENTS][2], 
			double dt, const grid_volume &gv,
			void *P_internal_data) const {
    (void) P; (void) W; (void) W_prev; (void) dt; (void) gv;
    (void) P_internal_data; // avoid warnings for unused params
  }

  // subtract all of the internal polarizations from the given f_minus_p
  // field.  Also given the fields array if it is needed for some reason.
  // Only update for ft fields.
  virtual void subtract_P(field_type ft,
			  realnum *f_minus_p[NUM_FIELD_COMPONENTS][2], 
			  void *P_internal_data) const {
    (void) ft; (void) f_minus_p; (void) P_internal_data;
  }

  // whether, for the given field W, Meep needs to allocate P[c]
  virtual bool needs_P(component c, int cmp,
		       realnum *W[NUM_FIELD_COMPONENTS][2]) const;

  // whether update_P will need the notowned part of W for this c
  // (which means that Meep will need to communicate it between chunks)
  virtual bool needs_W_notowned(component c,
				realnum *W[NUM_FIELD_COMPONENTS][2]) const;

  // whether update_P needs the W_prev field (from the previous timestep)
  virtual bool needs_W_prev() const { return false; }

  /* A susceptibility may be associated with any amount of internal
     data need to update the polarization field.  This includes the
     polarization field(s) itself.  It may also, for example, store
     the polarization field from previous timesteps, atomic-level
     populations, or other data.  These routines return the size of
     this internal-data array and initialize it. */
  virtual void* new_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
				  const grid_volume &gv) const { 
    (void) W; (void) gv; return 0;
  }
  virtual void delete_internal_data(void *data) const;
  virtual void init_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
			 double dt, const grid_volume &gv, void *data) const {
    (void) W; (void) dt; (void) gv; (void) data; }
  virtual void *copy_internal_data(void *data) const { (void)data; return 0; }

  /* The following methods are used in boundaries.cpp to set up any
     extra communications that may be necessary at chunk boundaries
     for the internal data of a susceptibility's polarization
     state. */

  /* the number of notowned fields/data in the internal data that
     are needed by update_P for the c Yee grid (note: we assume that we only
     have internal data for c's where we have external polarizations) */
  virtual int num_internal_notowned_needed(component c,
					   void *P_internal_data) const {
    (void) c; (void) P_internal_data; return 0; }
  /* the offset into the internal data of the n'th Yee-grid point in
     the c Yee grid for the inotowned internal field, where
     0 <= inotowned < size_internal_notowned_needed. */
  virtual realnum *internal_notowned_ptr(int inotowned, component c, int n,
					 void *P_internal_data) const {
    (void) inotowned; (void) n; (void) c; (void) P_internal_data; return 0; }
  
  /* same thing as above, except this gives (possibly complex)
     internal fields that need to be multiplied by the same phase
     factor as the fields at boundaries.  Note: we assume internal fields
     are complex if and only if !is_real (i.e. if EM fields are complex) */
  virtual int num_cinternal_notowned_needed(component c,
					   void *P_internal_data) const {
    (void) c; (void) P_internal_data; return 0; }
  // real/imaginary parts offsets for cmp = 0/1
  virtual realnum *cinternal_notowned_ptr(int inotowned, component c, int cmp, 
					  int n, 
					  void *P_internal_data) const {
    (void) inotowned; (void) n; (void) c; (void) cmp; (void) P_internal_data; 
    return 0; }
  
  susceptibility *next;
  int ntot;
  realnum *sigma[NUM_FIELD_COMPONENTS][5];

  /* trivial_sigma[c][d] is true only if *none* of the processes has a
     nontrivial sigma (c,d) component.  This differs, from sigma,
     which is non-NULL only if *this* process needs a nontrivial sigma
     (c,d).  Coordinated between processes at add_susceptibility, no
     communication elsewhere.  (We need this for boundary
     communcations between chunks, where one chunk might have sigma ==
     0 and the other != 0.) */
  bool trivial_sigma[NUM_FIELD_COMPONENTS][5];

private:
  static int cur_id; // unique id to assign to next susceptibility object
  int id; // id for this object and its clones, for comparison purposes
};

/* a Lorentzian susceptibility 
   \chi(\omega) = sigma * omega_0^2 / (\omega_0^2 - \omega^2 - i\gamma \omega)
  If no_omega_0_denominator is true, then we omit the omega_0^2 factor in the
  denominator to obtain a Drude model. */
class lorentzian_susceptibility : public susceptibility {
public:
  lorentzian_susceptibility(double omega_0, double gamma, bool no_omega_0_denominator = false) : omega_0(omega_0), gamma(gamma), no_omega_0_denominator(no_omega_0_denominator) {}
  virtual susceptibility *clone() const { return new lorentzian_susceptibility(*this); }
  virtual ~lorentzian_susceptibility() {}
  
  virtual void update_P(realnum *W[NUM_FIELD_COMPONENTS][2], 
			realnum *W_prev[NUM_FIELD_COMPONENTS][2], 
			double dt, const grid_volume &gv,
			void *P_internal_data) const;

  virtual void subtract_P(field_type ft,
			  realnum *f_minus_p[NUM_FIELD_COMPONENTS][2], 
			  void *P_internal_data) const;

  virtual void *new_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
				  const grid_volume &gv) const;
  virtual void init_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
			  double dt, const grid_volume &gv, void *data) const;
  virtual void *copy_internal_data(void *data) const;

  virtual int num_cinternal_notowned_needed(component c,
					    void *P_internal_data) const;
  virtual realnum *cinternal_notowned_ptr(int inotowned, component c, int cmp, 
					  int n, 
					  void *P_internal_data) const;
protected:
  double omega_0, gamma;
  bool no_omega_0_denominator;
};

/* like a Lorentzian susceptibility, but the polarization equation
   includes white noise with a specified amplitude */
class noisy_lorentzian_susceptibility : public lorentzian_susceptibility {
public:
  noisy_lorentzian_susceptibility(double noise_amp, double omega_0, double gamma, bool no_omega_0_denominator = false) : lorentzian_susceptibility(omega_0, gamma, no_omega_0_denominator), noise_amp(noise_amp) {}
  
  virtual susceptibility *clone() const { return new noisy_lorentzian_susceptibility(*this); }
  
  virtual void update_P(realnum *W[NUM_FIELD_COMPONENTS][2], 
			realnum *W_prev[NUM_FIELD_COMPONENTS][2], 
			double dt, const grid_volume &gv,
			void *P_internal_data) const;

protected:
  double noise_amp;
};

class multilevel_susceptibility : public susceptibility {
public:
  multilevel_susceptibility() : L(0), T(0), Gamma(0), N0(0), alpha(0), omega(0), gamma(0) {}
  multilevel_susceptibility(int L, int T,
			    const realnum *Gamma,
			    const realnum *N0,
			    const realnum *alpha,
			    const realnum *omega,
			    const realnum *gamma,
			    const realnum *sigmat);
  multilevel_susceptibility(const multilevel_susceptibility &from);
  virtual susceptibility *clone() const { return new multilevel_susceptibility(*this); }
  virtual ~multilevel_susceptibility();

  virtual void update_P(realnum *W[NUM_FIELD_COMPONENTS][2], 
			realnum *W_prev[NUM_FIELD_COMPONENTS][2], 
			double dt, const grid_volume &gv,
			void *P_internal_data) const;

  virtual void subtract_P(field_type ft,
			  realnum *f_minus_p[NUM_FIELD_COMPONENTS][2], 
			  void *P_internal_data) const;

  virtual void *new_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
				  const grid_volume &gv) const;
  virtual void init_internal_data(realnum *W[NUM_FIELD_COMPONENTS][2],
				  double dt, const grid_volume &gv, 
				  void *data) const;
  virtual void *copy_internal_data(void *data) const;
  virtual void delete_internal_data(void *data) const;

  virtual int num_cinternal_notowned_needed(component c,
					    void *P_internal_data) const;
  virtual realnum *cinternal_notowned_ptr(int inotowned, component c, int cmp, 
					  int n, 
					  void *P_internal_data) const;

  // always need notowned W and W_prev for E dot dP/dt terms
  virtual bool needs_W_notowned(component c,
				realnum *W[NUM_FIELD_COMPONENTS][2]) const {
    (void) c; (void) W;
    return true;
  }
  virtual bool needs_W_prev() const { return true; }

protected:
  int L; // number of atom levels
  int T; // number of optical transitions
  realnum *Gamma; // LxL matrix of relaxation rates Gamma[i*L+j] from i -> j
  realnum *N0; // L initial populations
  realnum *alpha; // LxT matrix of transition coefficients 1/omega
  realnum *omega; // T transition frequencies
  realnum *gamma; // T optical loss rates
  realnum *sigmat; // 5*T transition-specific sigma-diagonal factors
};

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
  
  realnum *read(const char *dataname, int *rank, int *dims, int maxrank);
  void write(const char *dataname, int rank, const int *dims, realnum *data,
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
		   realnum *data);
  void done_writing_chunks();
  
  void read_size(const char *dataname, int *rank, int *dims, int maxrank);
  void read_chunk(int rank, const int *chunk_start, const int *chunk_dims,
		  realnum *data);
  
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

typedef double (*pml_profile_func)(double u, void *func_data);

#define DEFAULT_SUBPIXEL_TOL 1e-4
#define DEFAULT_SUBPIXEL_MAXEVAL 100000

/* This class is used to compute position-dependent material properties
   like the dielectric function, permeability (mu), polarizability sigma, 
   nonlinearities, et cetera.  Simple cases of stateless functions are
   handled by canned subclasses below, but more complicated cases
   can be handled by creating a user-defined subclass of material_function.
   It is useful to group different properties into one class because
   it is likely that complicated implementations will share state between
   properties. */
class material_function {
  material_function(const material_function &ef) {(void)ef;} // prevent copying
public:
  material_function() {}
  virtual ~material_function() {}
  
  /* Specify a restricted grid_volume: all subsequent eps/sigma/etc
     calls will be for points inside v, until the next set_volume. */
  virtual void set_volume(const volume &v) {(void)v;}
  virtual void unset_volume(void) {} // unrestrict the grid_volume
  
  virtual double chi1p1(field_type ft, const vec &r) { (void)ft; (void)r; return 1.0; }

  /* scalar dielectric function */
  virtual double eps(const vec &r) { return chi1p1(E_stuff, r);  }

  /* scalar permeability function */
  virtual bool has_mu() { return false; } /* true if mu != 1 */
  virtual double mu(const vec &r) { return chi1p1(H_stuff, r);  }
  
  /* scalar conductivity function */
  virtual bool has_conductivity(component c) { (void)c; return false; }
  virtual double conductivity(component c, const vec &r) { 
    (void) c; (void)r; return 0.0; }  

  // fallback routine based on spherical quadrature
  vec normal_vector(field_type ft, const volume &v);

  /* Return c'th row of effective 1/(1+chi1) tensor in the given grid_volume v
     ... virtual so that e.g. libctl can override with more-efficient
     libctlgeom-based routines.  maxeval == 0 if no averaging desired. */
  virtual void eff_chi1inv_row(component c, double chi1inv_row[3],
			       const volume &v, 
			       double tol=DEFAULT_SUBPIXEL_TOL, 
			       int maxeval=DEFAULT_SUBPIXEL_MAXEVAL);
  
  /* polarizability sigma function: return c'th row of tensor */
  virtual void sigma_row(component c, double sigrow[3], const vec &r) {
    (void) c; (void) r; sigrow[0] = sigrow[1] = sigrow[2] = 0.0;
  }

  // Nonlinear susceptibilities
  virtual bool has_chi3(component c) { (void)c; return false; }
  virtual double chi3(component c, const vec &r) { (void)c; (void)r; return 0.0; }
  virtual bool has_chi2(component c) { (void)c; return false; }
  virtual double chi2(component c, const vec &r) { (void)c; (void)r; return 0.0; }
};

class simple_material_function : public material_function {
  double (*f)(const vec &);
  
public:
  simple_material_function(double (*func)(const vec &)) { f = func; }
  
  virtual ~simple_material_function() {}
  
  virtual double chi1p1(field_type ft, const vec &r) { (void)ft; return f(r); }
  virtual double eps(const vec &r) { return f(r); }
  virtual double mu(const vec &r) { return f(r); }
  virtual double conductivity(component c, const vec &r) { 
    (void)c; return f(r); }
  virtual void sigma_row(component c, double sigrow[3], const vec &r) {
    sigrow[0] = sigrow[1] = sigrow[2] = 0.0;
    sigrow[component_index(c)] = f(r);
  }
  virtual double chi3(component c, const vec &r) { (void)c; return f(r); }
  virtual double chi2(component c, const vec &r) { (void)c; return f(r); }
};

// in anisotropic_averaging.cpp:

// A simple polygon class, containing a number of 2D-points defining a polygon.
class simple_polygon {
public:
    // Default constructor creating polygon with zero points:
    simple_polygon() : _number_of_points(0), _polygon_points(0) {};

    /* Construct a simple_polygon from the double[2]-array points with length num_points.
     * (externally declared and initialized in this way: double (*points)[2] = new double[numpts][2];)
     * The matrix is copied internally, so can savely be deleted or reused after calling this constructor.
     */
    simple_polygon(const double (* const points)[2], std::size_t num_points);

    /* Construct a simple_polygon from the flattened points array matrix.
     * (row-major order: {x0, y0, x1, y1, ...})
     * dim1 is the number of points, dim2 must always be 2.
     * This constructor is needed for easier integration with the SWIG interface.
     * The matrix is copied internally, so can savely be deleted or reused after calling this constructor.
     */
    simple_polygon(const double* matrix, std::size_t dim1, std::size_t dim2 = 2);

    // copy-constructor:
    simple_polygon(const simple_polygon &pol);

    // swap function, needed for copy-and-swap in operator=
    friend void swap(simple_polygon& first, simple_polygon& second);

    // copy assignment operator:
    simple_polygon& operator=(simple_polygon other);

    virtual ~simple_polygon();

    std::size_t get_number_of_points() const { return _number_of_points; };
    double* &get_polygon_point(std::size_t index) const { return _polygon_points[index]; };

    // Clip the polygon by a rectangular boundary:
    virtual void clip_polygon(double left, double right, double bottom, double top);

    // Returns the absolute area of polygon.
    // For self-intersecting polygons, this will return the difference between
    // clockwise parts and anti-clockwise parts, not the real area.
    virtual double get_area();

private:
    std::size_t _number_of_points;
    double** _polygon_points;
};

/* A better polygon class. Like simple_polygon, but can additionally contain
 * inner polygons, i.e. empty polygonal holes of the parent polygon. */
class polygon : public simple_polygon{
public:
    /* Construct a polygon from the double[2]-array points with length num_points.
     * (externally declared and initialized in this way: double (*points)[2] = new double[numpts][2];)
     * The matrix is copied internally, so can savely be deleted or reused after calling this constructor. */
    polygon(const double (* const points)[2], std::size_t num_points) : simple_polygon(points, num_points) {};

    /* Construct a polygon from the flattened points array matrix.
     * (row-major order: {x0, y0, x1, y1, ...})
     * dim1 is the number of points, dim2 must always be 2.
//      * This constructor is needed for easier integration with the SWIG interface.
     * The matrix is copied internally, so can savely be deleted or reused after calling this constructor. */
    polygon(const double* matrix, std::size_t dim1, std::size_t dim2 = 2) : simple_polygon(matrix, dim1, dim2) {};

    // copy-constructor:
    polygon(const polygon &pol);

    // swap function, needed for copy-and-swap in operator=
    friend void swap(polygon& first, polygon& second);

    // copy assignment operator:
    polygon& operator=(polygon other);

    virtual ~polygon();

    // initialize with a simple_polygon:
    polygon(const simple_polygon &pol) : simple_polygon(pol) {};

    /* Add an inner polygon (empty polygonal hole of the parent polygon).
     * The parameters are the same than for the polygon constructor.
     * points is copied internally, so can savely be deleted or reused afterwards. */
    void add_inner_polygon(const double (* const points)[2], std::size_t num_points);
    /* Add an inner polygon (empty polygonal hole of the parent polygon).
     * The parameters are the same than for the polygon constructor.
     * matrix is copied internally, so can savely be deleted or reused afterwards. */
    void add_inner_polygon(const double* matrix, std::size_t dim1, std::size_t dim2 = 2);
    /* Add an inner polygon (empty polygonal hole of the parent polygon).
     * pol is copied internally, so can savely be deleted or reused afterwards. */
    void add_inner_polygon(const simple_polygon &pol);

    std::size_t get_number_inner_polygons() { return _inner_polygons.size(); };
    simple_polygon* get_inner_polygon(std::size_t index) { return _inner_polygons[index]; };

    // Clip the polygon by a rectangular boundary:
    virtual void clip_polygon(double left, double right, double bottom, double top);

    /* Returns the absolute area of polygon. This is the area of the parent
     * polygon minus the areas of all inner polygons.
     * For self-intersecting polygons, this will return the difference between
     * clockwise parts and anti-clockwise parts, not the real area. */
    virtual double get_area();

protected:
    std::vector<simple_polygon*> _inner_polygons;
};

/* A material stack defines a number of layers, each with its own thickness and
 * epsilon value, for use with material_function_for_polygons in 3D simulations.
 * These layers are stacked along the Z direction, while the polygons, which
 * are defined elsewhere, define on which polygonal areas on the X-Y plane
 * these stacks are located. */
class material_stack{
public:
    /* Define a material stack, consisting of number_of_layers layers, with
     * the bottom most layer starting at z=0. From bottom to top, the layers'
     * thicknesses and epsilon values are placed in arrays pointed to by
     * layer_thicknesses and layer_epsilons, respectively. If the topmost
     * layer ends before the end of the computational domain in z, the topmost
     * layer will be extended until the end of the computational domain.
     * Both arrays will be copied internally.
     */
    material_stack(const double* layer_thicknesses,
                   const double* layer_epsilons,
                   const int number_of_layers);

    // make a copy of an existing stack:
    material_stack(const material_stack &stack);
    ~material_stack();

    // return chi1p1 (epsilon) at z; The lowest layer starts at z=0.
    double chi1p1(const double &z) const;

    /**
     * Return the averaged normal component between lower_z and higher_z,
     * pointing from lower to higher eps. Return zero if there is no interface
     * between the two z values.
     */
    double normal(const double &lower_z, const double &higher_z) const;

    // Return averaged epsilon between lower_z and higher_z.
    double mean_eps(const double &lower_z, const double &higher_z) const;

    // Return averaged inverse epsilon between lower_z and higher_z.
    double mean_inveps(const double &lower_z, const double &higher_z) const;

    int get_number_of_layers() const {
        return _number_of_layers;
    };

protected:
    int _number_of_layers;
    // positions of upper boundaries of every layer:
    double* _high_z;
    // epsilon values of every layer:
    double* _epsilon;
};

/**
 * Defines one material type, which can be one epsilon value for 2D simulations
 * or one material_stack for 3D simulations, along with polygons where this
 * material type will be located on the X-Y plane.
 * Objects of this class never need to be created by the user. This class is
 * only used internally by material_function_for_polygons.
 */
class polygonal_material{
public:
    /**
     * Constructor for 3D simulations.
     * mat_stack defines the layer stack, it will be copied internally and can
     * be discarded afterwards.
     * block_size is the size of the internally used blocks, which usually is
     * half the pixel size, i.e. 0.5 * grid_volume.inva .
     * nx to nz are the number of these internal blocks for each cartesian
     * direction, usually twice the amount of pixels in each direction, i.e.
     * grid_volume.nx() * 2 etc.
     */
    polygonal_material(
        const material_stack* mat_stack,
        double block_size, int nx, int ny, int nz);

    /**
     * Constructor for 2D simulations.
     * epsilon is the dielectric epsilon value of the material.
     * block_size is the size of the internally used blocks, which usually is
     * half the pixel size, i.e. 0.5 * grid_volume.inva .
     * nx and nx are the number of these internal blocks for each cartesian
     * direction, usually twice the amount of pixels in each direction, i.e.
     * grid_volume.nx() * 2 etc.
     */
    polygonal_material(
        double epsilon, double block_size, int nx, int ny);

    ~polygonal_material();

    /**
     * Added polygon object will be copied internally, then clipped to the
     * computational cell */
    void add_polygon(const polygon& pol);

    double get_material_epsilon() const { return _epsilon; }
    const material_stack *get_material_stack() const
    { return _material_stack; }

    /**
     * Return true if the blocks between small_ivec and big_ivec have either
     * all an area of zero or all one. Only looks at ivecs' x and y
     * components. */
    bool is_trivial(const ivec &small_ivec, const ivec &big_ivec);

    /**
     * Return the surface vector of a straight interface
     * crossing the pixel (=2x2 ivec-pixels), centered at center.
     * If the interface is not straight, an approximate
     * surface vector is returned. */
    vec normal_vector(const ivec &center);

    /**
     * Return the portion of the block that is inside the material, i.e. the
     * portion inside any specified polygons, but outside inner polygons.
     * Ignores the z-component. */
    double get_area(const ivec &block)
    { return get_area(block.x(), block.y()); }
    /**
     * Return the portion of the block at (i, j) that is inside the material,
     * i.e. the portion inside any specified polygons, but outside inner
     * polygons. */
    double get_area(int i, int j);

protected:
    /**
     * General initializer for polygon_material. One of mat_stack or
     * epsilon must be 0, depending on whether it is 2D or 3D. */
    void init_polymat(
        const material_stack* mat_stack, double epsilon,
        double block_size, int nx, int ny, int nz);

    /**
     * Split all _polygons into blocks with size _blocksize x _blocksize, then
     * save the area they take up in each block to _areas. Afterwards, _polygons
     * will be empty and all split polygons will be deleted. */
    void update_splitting();

    /**
     * Split pol into (_nx x _ny) blocks and add (substract if inner_polygon)
     * the areas of each resulting part to protected field double** _areas
     * according to their (_nx x _ny) index.
     * col_arr and row_arr are only used as temporary buffers during splitting.
     * They must be allocated before the call, their items (the vectors)
     * will be emptied in the function (and any (x, y) points freed). */
    double split_polygon(
        simple_polygon* pol, bool inner_polygon,
        std::vector<double*>* col_arr, std::vector<double*>* row_arr);

    /**
     * Split the polygon in matrix (= array of (x, y) coordinates) along axis
     * (axis=0: split into columns; axis=1: split into rows) and place the
     * resulting polygons in arr (= array of vectors of (x, y) coordinates).
     * arr_size is the length of arr (= the number of cols/rows to split the
     * polygon into); matrix_len is the length of the input polygon.
     * The array arr must be initialized before call, its items (the vectors)
     * will be emptied in the function (and any (x, y) points freed).
     * It is important that the last point in matrix equals the first point.
     * The polygons in arr returned from this function will all have their last
     * points equal their first points. These polygons will all be scaled and
     * shifted along axis so that the point values in this direction will
     * be between 0 and 1.
     */
    void split_in_one_dimension(std::vector<double*>* arr, int arr_size,
                                double** matrix, int matrix_len, int axis);

    bool is_3D() const { return _nz != 0; };

    /**
     * Rotate vec v num times by pi/2 in anti-clockwise direction.
     * Negative num will rotate in clockwise direction.
     * If mirror_x_before or mirror_y_before is true, v will be flipped in x or
     * y direction, respectively, before rotation.
     * v will be changed in place and also returned.
     * If v is 3D, z-component will be kept unchanged. */
    vec rightanglerotate_2Dvec(
        vec v, int num, bool mirror_x_before = false,
        bool mirror_y_before = false) const;

    /**
     * The following three functions use geometric considerations to calculate
     * the (normalized, i.e. length=1) surface vector of an interface
     * intersecting two neigboring squares, each of area 1, such that both
     * squares are divided by the interface in two areas each. Needed for the
     * calculation are only the areas A1 and A2 (for square 1 and 2,
     * respectively) on one side of the intersecting interface, such that the
     * resulting vector point away from areas for which A1 and A2 are given.
     * The interface should be a straight line for correct results, otherwise
     * an approximately correct result is given. The three cases correspond to
     * different positions of the interface inside the squares. See the
     * implementations for details.
     * If the calculation is in 3D, the returned vec will also be 3D, with
     * a zero z component. */
    vec calc_grad_case1(const double A1, const double A2) const;
    vec calc_grad_case2(const double A1, const double A2) const;
    vec calc_grad_case3(const double A1, const double A2) const;

    /**
     * Return the normalized 2D surface vector of a flat material interface
     * intersecting a 4x4 grid of squares on the x-y plane. The single squares
     * each have a size of 1x1. Needed for the calculation are the areas the
     * material occupies in each of these squares. The area of one such square
     * might be 1 if the square is completely on one side of the interface
     * (completely inside material) or 0 if it is on the other side (completely
     * outside material), or a fraction inbetween if it is divided by the
     * interface. The input 'areas4x4' is a 4x4 matrix of these areas,
     * corresponding to the 4x4 grid of squares. The resulting vector will
     * point away from the material.
     * Mainly the areas in the center 2x2 grid of squares are relevant for the
     * surface vector. The additional squares outside are only needed to
     * determine the case of calculation, see calc_grad_case[1-3] above, where
     * it is not clear from the center 2x2 squares alone.
     * The interface should ideally be a straight line for exact results,
     * otherwise if an edge is present the algorithm will return approximately
     * correct values for the surface vector.
     * If the calculation is in 3D, the returned vec will also be 3D, with
     * a zero z component. */
    vec get_2D_gradient(double** areas4x4) const;

    /**
     * Return the (3D) surface vector of the material interface calculated
     * inside the pixel centered at pixel_center using the Lebedev quadrature
     * scheme, pointing away from the material (defined by polygons).
     * The 2D normal gradient2D (only x and y components regarded; averaged
     * inside the pixel) is needed to reconstruct the original polygon
     * interface from the occupied areas in the pixel. This is exact if the
     * original interface was a straight line, otherwise the gradient2D will
     * be an average already and so the interface used will be straightened.
     *
     * The use of this function does not really make sense in 2D, since there
     * the gradient2D is already the best solution. It is needed in 3D though
     * at pixels where there is both an interface in lateral as well as in z
     * direction.
     */
    vec get_gradient_from_sphere_quad(
        const ivec &pixel_center, const vec &gradient2D);

    bool splitting_done;

    /**
     * Polygons added by add_polygon will be added to this list. All these
     * polygons have the same eps (in 2D) or material stack (in 3D). */
    std::vector<polygon*> _polygons;

    // The material all polygons inside this class share:
    material_stack* _material_stack; // only used for 3D
    double _epsilon;                 // only used for 2D

    /**
     * Number of blocks to split polygons into. One block should have halve the
     * size of one pixel to allow calculating averaged eps for different
     * components (half pixel apart due to yee lattice). */
    int _nx, _ny, _nz;
    // The size of one such block ( = width = height = depth):
    double _block_size;

    /**
     * In this 2D array of areas, the 2D matrix corresponds to the _nx rows and
     * _ny columns the computational volume is split into. The areas are the
     * percentage of material in each block. Polygons only define the structure
     * on the x-y-plane, so the areas are not dependant on z. */
    double ** _areas;
};

/**
 * A material_function which allows defining the dielectric function for
 * structured layered media with polygons.
 * For 2D simulations, the dielectric function can be freely defined by adding
 * any amount of non-overlapping polygons and their epsilon value eps with
 * add_polygon(const polygon &pol, double eps);
 * For 3D simulations, instead of an epsilon value, all added polygons must be
 * assigned a material_stack, which defines the dielectric function along the
 * z axis, while the polygons describe the regions on the x-y-plane where the
 * stacks are located. Any amount of material stacks can be added with
 * add_material_stack.
 */
class material_function_for_polygons : public material_function {
public:
    /* creates a material function for polygons. thegv is the grid volume
     * object also used for the simulation.
     */
    material_function_for_polygons(const grid_volume &thegv);

    virtual ~material_function_for_polygons();

    /* Adds a material stack to the material function, for use in 3D case.
     * The stack consists of number_of_layers layers, with the bottom most
     * layer starting at z=0. From bottom to top, the layers' thicknesses and
     * epsilon values are placed in arrays pointed to by layer_thicknesses
     * and layer_epsilons, respectively. If the topmost layer ends before the
     * end of the computational domain in z, the topmost layer will be
     * extended until the end of the computational domain. Both arrays will
     * be copied internally. Returns mat_ID, which must be used for
     * add_polygon(polygon, mat_ID)
     */
    unsigned int add_material_stack(
        double* layer_thicknesses, double* layer_epsilons,
        int number_of_layers)
    {
        return add_material_stack(
            layer_thicknesses, number_of_layers, layer_epsilons, number_of_layers);
    };

    // Adds a material stack to the material function, for use in 3D case.
    // Returns mat_ID, which must be used for add_polygon(polygon, mat_ID)
    // This method is for easier integration with the SWIG interface:
    unsigned int add_material_stack(
        double* layer_thicknesses, int dim1,
        double* layer_epsilons, int dim2);

    // Adds a polygon with predefined stack (->mat_ID from add_material_stack)
    // to the simulation. Care must be taken that none of the polygons overlap.
    // Only for 3D simulation.
    // Added polygon object will be copied internally and can savely be deleted hereafter.
    void add_polygon(const polygon &pol, unsigned int mat_ID);

    // Adds a polygon with epsilon value eps to the simulation.
    // Care must be taken that none of the polygons overlap.
    // Added polygon object will be copied internally and can savely be deleted hereafter.
    void add_polygon(const polygon &pol, double eps);

    /**
     * Return c'th row of effective 1/(1+chi1) tensor in the given grid_volume v.
     * tol has no use. maxeval == 0 if only simple mean averaging desired,
     * otherwise, maxeval has no use. */
    virtual void eff_chi1inv_row(component c, double chi1inv_row[3],
                                 const volume &v,
                                 double tol = 0,
                                 int maxeval = 1);

    virtual double eps(const vec &r) { return mean_eps(r); };

    virtual double chi1p1(field_type ft, const vec &r) { (void)ft; return mean_eps(r); }

    // Return normal vector of interface averaged inside 1 pixel area (== 2x2 ivecs area) around center.
    vec normal_vector(const ivec &center);

    // Return mean epsilon of 1 pixel area (== 2x2 ivecs area) around center.
    double mean_eps(const ivec &center);
    double mean_eps(const vec &r)
    {   return mean_eps(round_to_block_edge(r)); };

    // Return mean inverse epsilon of 1 pixel area (== 2x2 ivecs area) around center.
    double mean_inveps(const ivec &center);

protected:

    /**
     * This rounds the vec v to the nearest ivec position, i.e. the nearest
     * border between neigboring blocks. This does not neccessarily return the
     * block index inside which the vec is. This is useful for vectors which
     * are already given on ivec positions. */
    ivec round_to_block_edge(const vec &v);

    bool is_3D() const { return nz != 0; };

    /**
     * Similar to add_material_stack, but for 2D, this adds a new item to
     * _polymats for every new epsilon value, and/or returns the index to
     * the corresponding polygon_material in _polymats. */
    unsigned int add_epsilon(double eps);

    // The items of this vector correspond to different epsilon/material_stacks;
    // All polygons of same eps/mat_stack are added to same vector item:
    std::vector<polygonal_material*> _polymats;

    // Number of blocks to split polygons into. One block should have halve the size of one pixel
    // to allow calculating averaged eps for different components (half pixel apart due to yee lattice).
    // The smooting diameter should stay 1 pixel (Meep default) or be an integer multiple thereof,
    // in which case material_function_for_polygons::eff_chi1inv_row() must be changed though.
    int nx, ny, nz;
    // the size of one such block ( = width = height = depth):
    double block_size;
};


class structure;

class structure_chunk {
 public:
  double a, Courant, dt; // res. a, Courant num., and timestep dt=Courant/a
  realnum *chi3[NUM_FIELD_COMPONENTS], *chi2[NUM_FIELD_COMPONENTS];
  realnum *chi1inv[NUM_FIELD_COMPONENTS][5];
  bool trivial_chi1inv[NUM_FIELD_COMPONENTS][5];
  realnum *conductivity[NUM_FIELD_COMPONENTS][5];
  realnum *condinv[NUM_FIELD_COMPONENTS][5]; // cache of 1/(1+conduct*dt/2)
  bool condinv_stale; // true if condinv needs to be recomputed
  double *sig[5], *kap[5], *siginv[5]; // conductivity array for uPML
  int sigsize[5]; // conductivity array size
  grid_volume gv;  // integer grid_volume that could be bigger than non-overlapping v below
  volume v;
  susceptibility *chiP[NUM_FIELD_TYPES]; // only E_stuff and H_stuff are used

  int refcount; // reference count of objects using this structure_chunk

  ~structure_chunk();
  structure_chunk(const grid_volume &gv,
            const volume &vol_limit, double Courant, int proc_num);
  structure_chunk(const structure_chunk *);
  void set_chi1inv(component c, material_function &eps,
                   bool use_anisotropic_averaging,
		   double tol, int maxeval);
  bool has_chi(component c, direction d) const;
  bool has_chisigma(component c, direction d) const;
  bool has_chi1inv(component c, direction d) const;
  void set_conductivity(component c, material_function &eps);
  void update_condinv();
  void set_chi3(component c, material_function &eps);
  void set_chi2(component c, material_function &eps);
  void use_pml(direction, double dx, double boundary_loc,
	       double Rasymptotic, double mean_stretch,
	       pml_profile_func pml_profile, void *pml_profile_data,
	       double pml_profile_integral, double pml_profile_integral_u);

  void add_susceptibility(material_function &sigma, field_type ft,
			  const susceptibility &sus);

  void mix_with(const structure_chunk *, double);

  int n_proc() const { return the_proc; } // Says which proc owns me!
  int is_mine() const { return the_is_mine; }

  void remove_susceptibilities();

  // monitor.cpp
  double get_chi1inv(component, direction, const ivec &iloc) const;
  double get_inveps(component c, direction d, const ivec &iloc) const {
    return get_chi1inv(c, d, iloc); }
  double max_eps() const;
 private:
  double pml_fmin;
  int the_proc;
  int the_is_mine;
};

double pml_quadratic_profile(double, void*);

// linked list of descriptors for boundary regions (currently just for PML)
class boundary_region {
public:
  typedef enum { NOTHING_SPECIAL, PML } boundary_region_kind;
  
  boundary_region() :
    kind(NOTHING_SPECIAL), thickness(0.0), Rasymptotic(1e-16), mean_stretch(1.0), pml_profile(NULL), pml_profile_data(NULL), pml_profile_integral(1.0), pml_profile_integral_u(1.0), d(NO_DIRECTION), side(Low), next(0) {}
  boundary_region(boundary_region_kind kind, double thickness, double Rasymptotic, double mean_stretch, pml_profile_func pml_profile, void* pml_profile_data, double pml_profile_integral, double pml_profile_integral_u, direction d, boundary_side side, boundary_region *next = 0) : kind(kind), thickness(thickness), Rasymptotic(Rasymptotic), mean_stretch(mean_stretch), pml_profile(pml_profile), pml_profile_data(pml_profile_data), pml_profile_integral(pml_profile_integral), pml_profile_integral_u(pml_profile_integral_u), d(d), side(side), next(next) {}

  boundary_region(const boundary_region &r) : kind(r.kind), thickness(r.thickness), Rasymptotic(r.Rasymptotic), mean_stretch(r.mean_stretch), pml_profile(r.pml_profile), pml_profile_data(r.pml_profile_data), pml_profile_integral(r.pml_profile_integral), pml_profile_integral_u(r.pml_profile_integral_u), d(r.d), side(r.side) { 
    next = r.next ? new boundary_region(*r.next) : 0;
  }

  ~boundary_region() { if (next) delete next; }
  
  void operator=(const boundary_region &r) {
    kind = r.kind; thickness = r.thickness; Rasymptotic = r.Rasymptotic; mean_stretch = r.mean_stretch;
    pml_profile = r.pml_profile; pml_profile_data = r.pml_profile_data;
    pml_profile_integral = r.pml_profile_integral;
    pml_profile_integral_u = r.pml_profile_integral_u;
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
      cur->Rasymptotic = pow(cur->Rasymptotic, strength_mult);
      cur = cur->next;
    }
    return r;
  }

  void apply(structure *s) const;
  void apply(const structure *s, structure_chunk *sc) const;
  bool check_ok(const grid_volume &gv) const;

private:
  boundary_region_kind kind;
  double thickness, Rasymptotic, mean_stretch;
  pml_profile_func pml_profile;
  void *pml_profile_data;
  double pml_profile_integral, pml_profile_integral_u;
  direction d;
  boundary_side side;
  boundary_region *next;
};

boundary_region pml(double thickness, direction d, boundary_side side,
		    double Rasymptotic = 1e-15, double mean_stretch = 1.0);
boundary_region pml(double thickness, direction d,
		    double Rasymptotic = 1e-15, double mean_stretch = 1.0);
boundary_region pml(double thickness,
		    double Rasymptotic = 1e-15, double mean_stretch = 1.0);
#define no_pml() boundary_region()

class structure {
 public:
  structure_chunk **chunks;
  int num_chunks;
  grid_volume gv, user_volume;
  double a, Courant, dt; // res. a, Courant num., and timestep dt=Courant/a
  volume v;
  symmetry S;
  const char *outdir;
  grid_volume *effort_volumes;
  double *effort;
  int num_effort_volumes;

  ~structure();
  structure();
  structure(const grid_volume &gv, material_function &eps,
	    const boundary_region &br = boundary_region(),
	    const symmetry &s = meep::identity(),
	    int num_chunks = 0, double Courant = 0.5,
	    bool use_anisotropic_averaging=false,
	    double tol=DEFAULT_SUBPIXEL_TOL,
	    int maxeval=DEFAULT_SUBPIXEL_MAXEVAL);
  structure(const grid_volume &gv, double eps(const vec &), 
	    const boundary_region &br = boundary_region(),
	    const symmetry &s = meep::identity(),
	    int num_chunks = 0, double Courant = 0.5,
	    bool use_anisotropic_averaging=false,
	    double tol=DEFAULT_SUBPIXEL_TOL,
	    int maxeval=DEFAULT_SUBPIXEL_MAXEVAL);
  structure(const structure *);
  structure(const structure &);

  void set_materials(material_function &mat,
		     bool use_anisotropic_averaging=true,
		     double tol=DEFAULT_SUBPIXEL_TOL,
		     int maxeval=DEFAULT_SUBPIXEL_MAXEVAL);
  void set_chi1inv(component c, material_function &eps,
                   bool use_anisotropic_averaging=true,
		   double tol=DEFAULT_SUBPIXEL_TOL,
		   int maxeval=DEFAULT_SUBPIXEL_MAXEVAL);
  bool has_chi(component c, direction d) const;
  void set_epsilon(material_function &eps,
                   bool use_anisotropic_averaging=true,
		   double tol=DEFAULT_SUBPIXEL_TOL,
		   int maxeval=DEFAULT_SUBPIXEL_MAXEVAL);
  void set_epsilon(double eps(const vec &),
                   bool use_anisotropic_averaging=true,
		   double tol=DEFAULT_SUBPIXEL_TOL,
		   int maxeval=DEFAULT_SUBPIXEL_MAXEVAL);
  void set_mu(material_function &eps,
	      bool use_anisotropic_averaging=true,
	      double tol=DEFAULT_SUBPIXEL_TOL,
	      int maxeval=DEFAULT_SUBPIXEL_MAXEVAL);
  void set_mu(double mu(const vec &),
	      bool use_anisotropic_averaging=true,
	      double tol=DEFAULT_SUBPIXEL_TOL,
	      int maxeval=DEFAULT_SUBPIXEL_MAXEVAL);
  void set_conductivity(component c, material_function &conductivity);
  void set_conductivity(component C, double conductivity(const vec &));
  void set_chi3(component c, material_function &eps);
  void set_chi3(material_function &eps);
  void set_chi3(double eps(const vec &));
  void set_chi2(component c, material_function &eps);
  void set_chi2(material_function &eps);
  void set_chi2(double eps(const vec &));
  
  void add_susceptibility(double sigma(const vec &), field_type c, const susceptibility &sus);
  void add_susceptibility(material_function &sigma, field_type c, const susceptibility &sus);
  void remove_susceptibilities();

  void set_output_directory(const char *name);
  void mix_with(const structure *, double);

  bool equal_layout(const structure &) const;
  void print_layout(void) const;

  // monitor.cpp
  double get_chi1inv(component, direction, const ivec &origloc) const;
  double get_chi1inv(component, direction, const vec &loc) const;
  double get_inveps(component c, direction d, const ivec &origloc) const {
    return get_chi1inv(c, d, origloc); }
  double get_inveps(component c, direction d, const vec &loc) const {
    return get_chi1inv(c, d, loc); }
  double get_eps(const vec &loc) const;
  double get_mu(const vec &loc) const;
  double max_eps() const;

  friend class boundary_region;

 private:
  void use_pml(direction d, boundary_side b, double dx);
  void add_to_effort_volumes(const grid_volume &new_effort_volume, 
			     double extra_effort);
  void choose_chunkdivision(const grid_volume &gv, int num_chunks,
			    const boundary_region &br, const symmetry &s);
  void check_chunks();
  void changing_chunks();
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
  // the following variable specifies whether the current
  // source is specified as a current or as an integrated
  // current (a dipole moment), if possible.  In the original Meep,
  // by default electric sources are integrated and magnetic
  // sources are not, but this may change.
  bool is_integrated;

  src_time() { is_integrated = true; current_time = nan; current_current = 0.0; next = NULL; }
  virtual ~src_time() { delete next; }
  src_time(const src_time &t) { 
       is_integrated = t.is_integrated;
       current_time = t.current_time;
       current_current = t.current_current;
       current_dipole = t.current_dipole;
       if (t.next) next = t.next->clone(); else next = NULL;
  }
  
  std::complex<double> dipole() const { return current_dipole; }
  std::complex<double> current() const { return current_current; }
  void update(double time, double dt) {
    if (time != current_time) {
      current_dipole = dipole(time);
      current_current = current(time, dt);
      current_time = time;
    }
  }

  // subclasses *can* override this method in order to specify the
  // current directly rather than as the derivative of dipole.
  // in that case you would probably ignore the dt argument.
  virtual std::complex<double> current(double time, double dt) const { 
    return ((dipole(time + dt) - dipole(time)) / dt);
  }

  double last_time_max() { return last_time_max(0.0); }
  double last_time_max(double after);
  
  src_time *add_to(src_time *others, src_time **added) const;
  src_time *next;

  // subclasses should override these methods:
  virtual std::complex<double> dipole(double time) const { (void)time; return 0; }
  virtual double last_time() const { return 0.0; }
  virtual src_time *clone() const { return new src_time(*this); }
  virtual bool is_equal(const src_time &t) const { (void)t; return 1; }
  virtual std::complex<double> frequency() const { return 0.0; }
  virtual void set_frequency(std::complex<double> f) { (void) f; }

 private:
  double current_time;
  std::complex<double> current_dipole, current_current;
};

bool src_times_equal(const src_time &t1, const src_time &t2);

// Gaussian-envelope source with given frequency, width, peak-time, cutoff
class gaussian_src_time : public src_time {
 public:
  gaussian_src_time(double f, double fwidth, double s = 5.0);
  gaussian_src_time(double f, double w, double start_time, double end_time);
  virtual ~gaussian_src_time() {}

  virtual std::complex<double> dipole(double time) const;
  virtual double last_time() const { return float(peak_time + cutoff); };
  virtual src_time *clone() const { return new gaussian_src_time(*this); }
  virtual bool is_equal(const src_time &t) const;
  virtual std::complex<double> frequency() const { return freq; }
  virtual void set_frequency(std::complex<double> f) { freq = real(f); }

 private:
  double freq, width, peak_time, cutoff;
};

// Continuous (CW) source with (optional) slow turn-on and/or turn-off.
class continuous_src_time : public src_time {
 public:
  continuous_src_time(std::complex<double> f, double w = 0.0, 
		      double st = 0.0, double et = infinity,
		      double s = 3.0) : freq(f), width(w), start_time(float(st)),
					end_time(float(et)), slowness(s) {}
  virtual ~continuous_src_time() {}
  
  virtual std::complex<double> dipole(double time) const;
  virtual double last_time() const { return end_time; };
  virtual src_time *clone() const { return new continuous_src_time(*this); }
  virtual bool is_equal(const src_time &t) const;
  virtual std::complex<double> frequency() const { return freq; }
  virtual void set_frequency(std::complex<double> f) { freq = f; }
  
 private:
  std::complex<double> freq;
  double width, start_time, end_time, slowness;
};

// user-specified source function with start and end times
class custom_src_time : public src_time {
 public:
  custom_src_time(std::complex<double> (*func)(double t, void *), void *data,
		  double st = -infinity, double et = infinity)
    : func(func), data(data), start_time(float(st)), end_time(float(et)) {}
  virtual ~custom_src_time() {}
  
  virtual std::complex<double> current(double time, double dt) const { 
    if (is_integrated) return src_time::current(time,dt);
    else return dipole(time);
  }
  virtual std::complex<double> dipole(double time) const { float rtime = float(time);
    if (rtime >= start_time && rtime <= end_time) return func(time,data); else return 0.0; }
  virtual double last_time() const { return end_time; };
  virtual src_time *clone() const { return new custom_src_time(*this); }
  virtual bool is_equal(const src_time &t) const;
  
 private:
  std::complex<double> (*func)(double t, void *);
  void *data;
  double start_time, end_time;
};

class monitor_point {
 public:
  monitor_point();
  ~monitor_point();
  vec loc;
  double t;
  std::complex<double> f[NUM_FIELD_COMPONENTS];
  monitor_point *next;

  std::complex<double> get_component(component);
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
                         std::complex<double> **a, std::complex<double> **f, int *numout,
                         double fmin=0.0, double fmax=0.0, int maxbands=100);
  // harminv works much like fourier_transform, except that it is not yet
  // implemented.
  void harminv(component w,
               std::complex<double> **a, std::complex<double> **f,
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
	    std::complex<double> scale_,
	    std::complex<double> extra_weight_,
	    component c_,
	    bool use_centered_grid,
            ivec shift_, const symmetry &S_, int sn_, int vc,
	    const void *data_);
  ~dft_chunk();
  
  void update_dft(double time);

  void scale_dft(std::complex<double> scale);

  void operator-=(const dft_chunk &chunk);

  // the frequencies to loop_in_chunks
  double omega_min, domega;
  int Nomega;

  component c; // component to DFT (possibly transformed by symmetry)

  int N; // number of spatial points (on epsilon grid)
  std::complex<realnum> *dft; // N x Nomega array of DFT values.

  struct dft_chunk *next_in_chunk; // per-fields_chunk list of DFT chunks
  struct dft_chunk *next_in_dft; // next for this particular DFT vol./component

  /* When computing things like -0.5*|E|^2 for the stress tensor,
     we cannot incorporate the minus sign into the scale factor
     because we only ever compute |scale|^2.  Thus, it is necessary
     to store an additional weight factor with the dft_chunk to record
     any additional negative or complex weight factor to be used
     in computations involving the fourier-transformed fields.  Because
     it is used in computations involving dft[...], it needs to be public. */
     std::complex<double> extra_weight;

  // parameters passed from field_integrate:
  fields_chunk *fc;
  ivec is, ie;
  vec s0, s1, e0, e1;
  double dV0, dV1;
  bool sqrt_dV_and_interp_weights;
  std::complex<double> scale; // scale factor * phase from shift and symmetry
  ivec shift;
  symmetry S; int sn;

  // cache of exp(iwt) * scale, of length Nomega
  std::complex<realnum> *dft_phase;

  int avg1, avg2; // index offsets for average to get epsilon grid

  int vc; // component descriptor from the original volume
};

void save_dft_hdf5(dft_chunk *dft_chunks, component c, h5file *file,
		   const char *dprefix = 0);
void load_dft_hdf5(dft_chunk *dft_chunks, component c, h5file *file,
		   const char *dprefix = 0);
void save_dft_hdf5(dft_chunk *dft_chunks, const char *name, h5file *file,
		   const char *dprefix = 0);
void load_dft_hdf5(dft_chunk *dft_chunks, const char *name, h5file *file,
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

  void save_hdf5(fields &f, const char *fname, const char *dprefix = 0,
		 const char *prefix = 0);
  void load_hdf5(fields &f, const char *fname, const char *dprefix = 0,
		 const char *prefix = 0);

  void scale_dfts(std::complex<double> scale);

  void remove();

  double freq_min, dfreq;
  int Nfreq;
  dft_chunk *E, *H;
  component cE, cH;
};

// stress.cpp (normally created with fields::add_dft_force)
class dft_force {
public:
  dft_force(dft_chunk *offdiag1_, dft_chunk *offdiag2_, dft_chunk *diag_,
	    double fmin, double fmax, int Nf);
  dft_force(const dft_force &f);

  double *force();

  void save_hdf5(h5file *file, const char *dprefix = 0);
  void load_hdf5(h5file *file, const char *dprefix = 0);

  void operator-=(const dft_force &fl);

  void save_hdf5(fields &f, const char *fname, const char *dprefix = 0,
		 const char *prefix = 0);
  void load_hdf5(fields &f, const char *fname, const char *dprefix = 0,
		 const char *prefix = 0);

  void scale_dfts(std::complex<double> scale);

  void remove();

  double freq_min, dfreq;
  int Nfreq;
  dft_chunk *offdiag1, *offdiag2, *diag;
};

// near2far.cpp (normally created with fields::add_dft_near2far)
class dft_near2far {
public:
  /* fourier tranforms of tangential E and H field components in a
     medium with the given scalar eps and mu */
  dft_near2far(dft_chunk *F,
               double fmin, double fmax, int Nf, double eps, double mu);
  dft_near2far(const dft_near2far &f);

  /* return an array (Ex,Ey,Ez,Hx,Hy,Hz) x Nfreq of the far fields at x */
  std::complex<double> *farfield(const vec &x);

  /* like farfield, but requires F to be Nfreq*6 preallocated array, and
     does *not* perform the reduction over processes...an MPI allreduce
     summation by the caller is required to get the final result ... used
     by other output routine to efficiently get far field on a grid of pts */
  void farfield_lowlevel(std::complex<double> *F, const vec &x);

  /* output far fields on a grid to an HDF5 file */
  void save_farfields(const char *fname, const char *prefix,
                      const volume &where, double resolution);

  void save_hdf5(h5file *file, const char *dprefix = 0);
  void load_hdf5(h5file *file, const char *dprefix = 0);

  void operator-=(const dft_near2far &fl);

  void save_hdf5(fields &f, const char *fname, const char *dprefix = 0,
		 const char *prefix = 0);
  void load_hdf5(fields &f, const char *fname, const char *dprefix = 0,
		 const char *prefix = 0);

  void scale_dfts(std::complex<double> scale);

  void remove();

  double freq_min, dfreq;
  int Nfreq;
  dft_chunk *F;
  double eps, mu;
};

/* Class to compute local-density-of-states spectra: the power spectrum
   P(omega) of the work done by the sources.  Specialized to handle only
   the case where all sources have the same time dependence, which greatly
   simplifies things because then we can do the spatial integral of E*J
   *first* and then do the Fourier transform, eliminating the need to
   store the Fourier transform per point or per current. */
class dft_ldos {
public:
  dft_ldos(double freq_min, double freq_max, int Nfreq);
  ~dft_ldos() { delete[] Fdft; delete[] Jdft; }

  void update(fields &f); // to be called after each timestep
  double *ldos() const; // returns array of Nomega values (after last timestep)
  std::complex<double> *F() const; // returns Fdft
  std::complex<double> *J() const; // returns Jdft

private:
  std::complex<realnum> *Fdft; // Nomega array of field * J*(x) DFT values
  std::complex<realnum> *Jdft; // Nomega array of J(t) DFT values
  double Jsum; // sum of |J| over all points
public:
  double omega_min, domega;
  int Nomega;
};

enum in_or_out { Incoming=0, Outgoing };
enum connect_phase { CONNECT_PHASE = 0, CONNECT_NEGATE=1, CONNECT_COPY=2 };

// data for each susceptibility
typedef struct polarization_state_s {
  void *data; // internal polarization data for the susceptibility
  const susceptibility *s;
  struct polarization_state_s *next; // linked list
} polarization_state;

class fields_chunk {
 public:
  realnum *f[NUM_FIELD_COMPONENTS][2]; // fields at current time

  // auxiliary fields needed for PML (at least in some components)
  realnum *f_u[NUM_FIELD_COMPONENTS][2]; // integrated from D/B
  realnum *f_w[NUM_FIELD_COMPONENTS][2]; // E/H integrated from these
  realnum *f_cond[NUM_FIELD_COMPONENTS][2]; // aux field for PML+conductivity

  /* sometimes, to synchronize the E and H fields, e.g. for computing
     flux at a given time, we need to timestep H by 1/2; in this case
     we save backup copies of (some of) the fields to resume timestepping */
  realnum *f_backup[NUM_FIELD_COMPONENTS][2];
  realnum *f_u_backup[NUM_FIELD_COMPONENTS][2];
  realnum *f_w_backup[NUM_FIELD_COMPONENTS][2];
  realnum *f_cond_backup[NUM_FIELD_COMPONENTS][2];

  // W (or E/H) field from prev. timestep, only stored if needed by update_pols
  realnum *f_w_prev[NUM_FIELD_COMPONENTS][2];

  // used to store D-P and B-P, e.g. when P implements dispersive media
  realnum *f_minus_p[NUM_FIELD_COMPONENTS][2];

  realnum *f_rderiv_int; // cache of helper field for 1/r d(rf)/dr derivative

  dft_chunk *dft_chunks;

  realnum **zeroes[NUM_FIELD_TYPES]; // Holds pointers to metal points.
  int num_zeroes[NUM_FIELD_TYPES];
  realnum **connections[NUM_FIELD_TYPES][CONNECT_COPY+1][Outgoing+1];
  int num_connections[NUM_FIELD_TYPES][CONNECT_COPY+1][Outgoing+1];
  std::complex<realnum> *connection_phases[NUM_FIELD_TYPES];

  int npol[NUM_FIELD_TYPES]; // only E_stuff and H_stuff are used
  polarization_state *pol[NUM_FIELD_TYPES]; // array of npol[i] polarization_state structures

  double a, Courant, dt; // res. a, Courant num., and timestep dt=Courant/a
  grid_volume gv;
  volume v;
  double m; // angular dependence in cyl. coords
  bool zero_fields_near_cylorigin; // fields=0 m pixels near r=0 for stability
  double beta;
  int is_real;
  bandsdata *bands;
  src_vol *sources[NUM_FIELD_TYPES];
  structure_chunk *new_s;
  structure_chunk *s;
  const char *outdir;

  fields_chunk(structure_chunk *, const char *outdir, double m,
	       double beta, bool zero_fields_near_cylorigin);

  fields_chunk(const fields_chunk &);
  ~fields_chunk();

  // step.cpp
  double peek_field(component, const vec &);

  void use_real_fields();
  bool have_component(component c, bool is_complex = false) {
    switch (c) {
    case Dielectric: case Permeability:
      return !is_complex;
    default:
      return (f[c][0] && f[c][is_complex]);
    }
  }

  double last_source_time();
  // monitor.cpp
  std::complex<double> get_field(component, const ivec &) const;

  // for non-collective interpolation:
  volume get_field_gv(component) const;
  std::complex<double> get_field(component, const vec &) const;

  double get_chi1inv(component, direction, const ivec &iloc) const;
  
  void backup_component(component c);
  void average_with_backup(component c);
  void restore_component(component c);
  
  void set_output_directory(const char *name);
  void verbose(int gv=1) { verbosity = gv; }

  double count_volume(component);
  friend class fields;

  int n_proc() const { return s->n_proc(); };
  int is_mine() const { return s->is_mine(); };
  // boundaries.cpp
  void zero_metal(field_type);
  bool needs_W_notowned(component c);
  // fields.cpp
  void remove_sources();
  void remove_susceptibilities();
  void zero_fields();

  // update_eh.cpp
  bool needs_W_prev(component c) const;
  bool update_eh(field_type ft, bool skip_w_components = false);

  bool alloc_f(component c);
  void figure_out_step_plan();

  void set_solve_cw_omega(std::complex<double> omega) {
    doing_solve_cw = true;
    solve_cw_omega = omega;
  }
  void unset_solve_cw_omega() {
    doing_solve_cw = false;
    solve_cw_omega = 0.0;
  }

 private: 
  // we set a flag during cw_solve to replace some
  // time-dependent stuff with the analogous frequency-domain operation
  bool doing_solve_cw; // true when inside solve_cw
  std::complex<double> solve_cw_omega; // current omega for solve_cw

  int verbosity; // Turn on verbosity for debugging purposes...
  // fields.cpp
  bool have_plus_deriv[NUM_FIELD_COMPONENTS], have_minus_deriv[NUM_FIELD_COMPONENTS];
  component plus_component[NUM_FIELD_COMPONENTS], minus_component[NUM_FIELD_COMPONENTS];
  direction plus_deriv_direction[NUM_FIELD_COMPONENTS],
            minus_deriv_direction[NUM_FIELD_COMPONENTS];
  // bands.cpp
  void record_bands(int tcount);
  // step.cpp
  void phase_in_material(structure_chunk *s);
  void phase_material(int phasein_time);
  bool step_db(field_type ft);
  void step_source(field_type ft, bool including_integrated);
  bool update_pols(field_type ft);
  void calc_sources(double time);

  // initialize.cpp
  void initialize_field(component, std::complex<double> f(const vec &));
  void initialize_with_nth_te(int n, double kz);
  void initialize_with_nth_tm(int n, double kz);
  // boundaries.cpp
  void alloc_extra_connections(field_type, connect_phase, in_or_out, int);
  // dft.cpp
  void update_dfts(double timeE, double timeH);

  void changing_structure();
};

enum boundary_condition { Periodic=0, Metallic, Magnetic, None };
enum time_sink { Connecting, Stepping, Boundaries, MpiTime,
                 FieldOutput, FourierTransforming, Other };

typedef void (*field_chunkloop)(fields_chunk *fc, int ichunk, component cgrid,
				ivec is, ivec ie,
				vec s0, vec s1, vec e0, vec e1,
				double dV0, double dV1,
				ivec shift, std::complex<double> shift_phase, 
				const symmetry &S, int sn,
				void *chunkloop_data);
typedef std::complex<double> (*field_function)(const std::complex<double> *fields,
					   const vec &loc,
					   void *integrand_data_);
typedef double (*field_rfunction)(const std::complex<double> *fields,
				   const vec &loc,
				   void *integrand_data_);

field_rfunction derived_component_func(derived_component c, const grid_volume &gv,
				       int &nfields, component cs[12]);

class fields {
 public:
  int num_chunks;
  fields_chunk **chunks;
  src_time *sources;
  flux_vol *fluxes;
  symmetry S;

  // The following is an array that is num_chunks by num_chunks.  Actually
  // it is two arrays, one for the imaginary and one for the real part.
  realnum **comm_blocks[NUM_FIELD_TYPES];
  // This is the same size as each comm_blocks array, and store the sizes
  // of the comm blocks themselves for each connection-phase type
  int *comm_sizes[NUM_FIELD_TYPES][CONNECT_COPY+1];
  int comm_size_tot(int f, int pair) const {
    int sum = 0; for (int ip=0; ip<3; ++ip) sum+=comm_sizes[f][ip][pair];
    return sum;
  }

  double a, dt; // The resolution a and timestep dt=Courant/a
  grid_volume gv, user_volume;
  volume v;
  double m;
  double beta;
  int t, phasein_time, is_real;
  std::complex<double> k[5], eikna[5];
  double coskna[5], sinkna[5];
  boundary_condition boundaries[2][5];
  bandsdata *bands;
  char *outdir;

  // fields.cpp methods:
  fields(structure *, double m=0, double beta=0,
	 bool zero_fields_near_cylorigin=true);
  fields(const fields &);
  ~fields();
  bool equal_layout(const fields &f) const;
  void use_real_fields();
  void zero_fields();
  void remove_sources();
  void remove_susceptibilities();
  void remove_fluxes();
  void reset();

  // time.cpp
  double time_spent_on(time_sink);
  void print_times();
  // boundaries.cpp
  void set_boundary(boundary_side,direction,boundary_condition);
  void use_bloch(direction d, double k) { use_bloch(d, (std::complex<double>) k); }
  void use_bloch(direction, std::complex<double> kz);
  void use_bloch(const vec &k);
  vec lattice_vector(direction) const;
  // update_eh.cpp
  void update_eh(field_type ft, bool skip_w_components = false);

  volume total_volume(void) const;

  // h5fields.cpp:
  // low-level function:
  void output_hdf5(h5file *file, const char *dataname,
		   int num_fields, const component *components,
		   field_function fun, void *fun_data_, int reim,
		   const volume &where,
		   bool append_data = false,
		   bool single_precision = false);
  // higher-level functions
  void output_hdf5(const char *dataname,  // OUTPUT COMPLEX-VALUED FUNCTION
		   int num_fields, const component *components,
		   field_function fun, void *fun_data_,
		   const volume &where,
		   h5file *file = 0,
		   bool append_data = false,
		   bool single_precision = false,
		   const char *prefix = 0,
		   bool real_part_only = false);
  void output_hdf5(const char *dataname,  // OUTPUT REAL-VALUED FUNCTION
		   int num_fields, const component *components,
		   field_rfunction fun, void *fun_data_,
		   const volume &where,
		   h5file *file = 0,
		   bool append_data = false,
		   bool single_precision = false,
		   const char *prefix = 0);
  void output_hdf5(component c,   // OUTPUT FIELD COMPONENT (or Dielectric)
		   const volume &where,
		   h5file *file = 0,
		   bool append_data = false,
		   bool single_precision = false,
		   const char *prefix = 0);
  void output_hdf5(derived_component c,   // OUTPUT DERIVED FIELD COMPONENT
		   const volume &where,
		   h5file *file = 0,
		   bool append_data = false,
		   bool single_precision = false,
		   const char *prefix = 0);
  h5file *open_h5file(const char *name, 
		      h5file::access_mode mode = h5file::WRITE,
		      const char *prefix = NULL, bool timestamp = false);
  const char *h5file_name(const char *name,
			  const char *prefix = NULL, bool timestamp = false);

  // step.cpp methods:
  double last_step_output_wall_time;
  int last_step_output_t;
  void step();

  // when comparing times, e.g. for source cutoffs, it
  // is useful to round to float to avoid gratuitous sensitivity
  // to floating-point roundoff error
  inline double round_time() const { return float(t*dt); };
  inline double time() const { return t*dt; };

  // cw_fields.cpp:
  bool solve_cw(double tol, int maxiters, std::complex<double> frequency, int L=2);
  bool solve_cw(double tol = 1e-8, int maxiters = 10000, int L=2);

  // sources.cpp:
  double last_source_time();
  void add_point_source(component c, double freq, double width, double peaktime,
                        double cutoff, const vec &, std::complex<double> amp = 1.0,
                        int is_continuous = 0);
  void add_point_source(component c, const src_time &src,
                        const vec &, std::complex<double> amp = 1.0);
  void add_volume_source(component c, const src_time &src,
			 const volume &, 
			 std::complex<double> A(const vec &),
			 std::complex<double> amp = 1.0);
  void add_volume_source(component c, const src_time &src,
			 const volume &, 
			 std::complex<double> amp = 1.0);
  void require_component(component c);

  // mpb.cpp
  void add_eigenmode_source(component c, const src_time &src,
			    direction d, const volume &where,
			    const volume &eig_vol,
			    int band_num, 
			    const vec &kpoint, bool match_frequency,
			    int parity,
			    double eig_resolution, double eigensolver_tol,
			    std::complex<double> amp,
			    std::complex<double> A(const vec &) = 0);

  // initialize.cpp:
  void initialize_field(component, std::complex<double> f(const vec &));
  void initialize_with_nth_te(int n);
  void initialize_with_nth_tm(int n);
  void initialize_with_n_te(int n);
  void initialize_with_n_tm(int n);
  int phase_in_material(const structure *s, double time);
  int is_phasing();

  // loop_in_chunks.cpp
  void loop_in_chunks(field_chunkloop chunkloop, void *chunkloop_data,
		      const volume &where,
		      component cgrid = Centered,
		      bool use_symmetry = true,
		      bool snap_unit_dims = false);
  
  // integrate.cpp
  std::complex<double> integrate(int num_fields, const component *components,
			    field_function fun, void *fun_data_,
			    const volume &where,
			    double *maxabs = 0);
  double integrate(int num_fields, const component *components,
		   field_rfunction fun, void *fun_data_,
		   const volume &where,
		   double *maxabs = 0);
  std::complex<double> integrate2(const fields &fields2,
			     int num_fields1,
			     const component *components1,
			     int num_fields2, 
			     const component *components2,
			     field_function integrand,
			     void *integrand_data_,
			     const volume &where,
			     double *maxabs = 0);
  double integrate2(const fields &fields2,
		    int num_fields1, const component *components1,
		    int num_fields2, const component *components2,
		    field_rfunction integrand,
		    void *integrand_data_,
		    const volume &where,
		    double *maxabs = 0);


  double max_abs(int num_fields, const component *components,
		 field_function fun, void *fun_data_,
		 const volume &where);
  double max_abs(int num_fields, const component *components,
		 field_rfunction fun, void *fun_data_,
		 const volume &where);

  double max_abs(int c, const volume &where);
  double max_abs(component c, const volume &where);
  double max_abs(derived_component c, const volume &where);
  
  // dft.cpp
  dft_chunk *add_dft(component c, const volume &where,
		     double freq_min, double freq_max, int Nfreq,
		     bool include_dV_and_interp_weights = true,
		     std::complex<double> weight = 1.0, dft_chunk *chunk_next = 0,
		     bool sqrt_dV_and_interp_weights = false,
		     std::complex<double> extra_weight = 1.0,
		     bool use_centered_grid = true, int vc = 0);
  dft_chunk *add_dft_pt(component c, const vec &where,
			double freq_min, double freq_max, int Nfreq);
  dft_chunk *add_dft(const volume_list *where,
		     double freq_min, double freq_max, int Nfreq,
		     bool include_dV = true);
  void update_dfts();
  dft_flux add_dft_flux(direction d, const volume &where,
			double freq_min, double freq_max, int Nfreq);
  dft_flux add_dft_flux_box(const volume &where,
			    double freq_min, double freq_max, int Nfreq);
  dft_flux add_dft_flux_plane(const volume &where,
			      double freq_min, double freq_max, int Nfreq);
  dft_flux add_dft_flux(const volume_list *where,
			double freq_min, double freq_max, int Nfreq);

  // stress.cpp
  dft_force add_dft_force(const volume_list *where,
			  double freq_min, double freq_max, int Nfreq);

  // near2far.cpp
  dft_near2far add_dft_near2far(const volume_list *where,
                                double freq_min, double freq_max, int Nfreq);
  // monitor.cpp
  double get_chi1inv(component, direction, const vec &loc) const;
  double get_inveps(component c, direction d, const vec &loc) const {
    return get_chi1inv(c, d, loc);
  }
  double get_eps(const vec &loc) const;
  double get_mu(const vec &loc) const;
  void get_point(monitor_point *p, const vec &) const;
  monitor_point *get_new_point(const vec &, monitor_point *p=NULL) const;
  
  void prepare_for_bands(const vec &, double end_time, double fmax=0,
                         double qmin=1e300, double frac_pow_min=0.0);
  void record_bands();
  std::complex<double> get_band(int n, int maxbands=100);
  void grace_bands(grace *, int maxbands=100);
  void output_bands(FILE *, const char *, int maxbands=100);
  std::complex<double> get_field(int c, const vec &loc) const;
  std::complex<double> get_field(component c, const vec &loc) const;
  double get_field(derived_component c, const vec &loc) const;

  // energy_and_flux.cpp
  void synchronize_magnetic_fields();
  void restore_magnetic_fields();
  double energy_in_box(const volume &);
  double electric_energy_in_box(const volume &);
  double magnetic_energy_in_box(const volume &);
  double thermo_energy_in_box(const volume &);
  double total_energy();
  double field_energy_in_box(const volume &);
  double field_energy_in_box(component c, const volume &);
  double field_energy();
  double flux_in_box_wrongH(direction d, const volume &);
  double flux_in_box(direction d, const volume &);
  flux_vol *add_flux_vol(direction d, const volume &where);
  flux_vol *add_flux_plane(const volume &where);
  flux_vol *add_flux_plane(const vec &p1, const vec &p2);
  double electric_energy_max_in_box(const volume &where);
  double modal_volume_in_box(const volume &where);
  double electric_sqr_weighted_integral(double (*deps)(const vec &),
				       const volume &where);
  double electric_energy_weighted_integral(double (*f)(const vec &),
					   const volume &where);

  void set_output_directory(const char *name);
  void verbose(int gv=1);
  double count_volume(component);
  // fields.cpp
  bool have_component(component);
  // material.cpp
  double max_eps() const;
  // step.cpp
  void step_boundaries(field_type);

  bool nosize_direction(direction d) const;
  direction normal_direction(const volume &where) const;

  // casimir.cpp
  std::complex<double> casimir_stress_dct_integral(direction dforce,
					      direction dsource,
					      double mx, double my, double mz,
					      field_type ft,
					      volume where,
					      bool is_bloch = false);

  void set_solve_cw_omega(std::complex<double> omega);
  void unset_solve_cw_omega();

 private: 
  int verbosity; // Turn on verbosity for debugging purposes...
  int synchronized_magnetic_fields; // count number of nested synchs
  double last_wall_time;
#define MEEP_TIMING_STACK_SZ 10
  time_sink working_on, was_working_on[MEEP_TIMING_STACK_SZ];
  double times_spent[Other+1];
  // fields.cpp
  void figure_out_step_plan();
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
  bool locate_point_in_user_volume(ivec *, std::complex<double> *phase) const;
  void locate_volume_source_in_user_volume(const vec p1, const vec p2, vec newp1[8], vec newp2[8],
                                           std::complex<double> kphase[8], int &ncopies) const;
  // mympi.cpp
  void boundary_communications(field_type);
  // step.cpp
  void phase_material();
  void step_db(field_type ft);
  void step_source(field_type ft, bool including_integrated = false);
  void update_pols(field_type ft);
  void calc_sources(double tim);
  int cluster_some_bands_cleverly(double *tf, double *td, std::complex<double> *ta,
                                  int num_freqs, int fields_considered, int maxbands,
                                  std::complex<double> *fad, double *approx_power);
  void out_bands(FILE *, const char *, int maxbands);
  std::complex<double> *clever_cluster_bands(int maxbands, double *approx_power = NULL);
public:
  // monitor.cpp
  std::complex<double> get_field(component c, const ivec &iloc) const;
  double get_chi1inv(component, direction, const ivec &iloc) const;
  // boundaries.cpp
  bool locate_component_point(component *, ivec *, std::complex<double> *) const;
};

class flux_vol {
 public:
  flux_vol(fields *f_, direction d_, const volume &where_) : where(where_) {
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
  volume where;
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
  FILE *f;
  char *fn, *dn;
  grace_point *pts;
  int set_num,sn;
};

// The following is a utility function to parse the executable name use it
// to come up with a directory name, avoiding overwriting any existing
// directory, unless the source file hasn't changed.

const char *make_output_directory(const char *exename, const char *jobname = NULL);
void trash_output_directory(const char *dirname);
FILE *create_output_file(const char *dirname, const char *fname);

// The following allows you to hit ctrl-C to tell your calculation to stop
// and clean up.
void deal_with_ctrl_c(int stop_now = 2);
// When a ctrl_c is called, the following variable (which starts with a
// zero value) is incremented.
extern int interrupt;

int do_harminv(std::complex<double> *data, int n, double dt,
	       double fmin, double fmax, int maxbands,
	       std::complex<double> *amps, double *freq_re, double *freq_im,
	       double *errors = NULL,
	       double spectral_density = 1.1, double Q_thresh = 50,
	       double rel_err_thresh = 1e20, double err_thresh = 0.01, 
	       double rel_amp_thresh = -1, double amp_thresh = -1);

std::complex<double> *make_casimir_gfunc(double T, double dt, double sigma, field_type ft,
				std::complex<double> (*eps_func)(std::complex<double> omega) = 0,
				double Tfft = 0);

std::complex<double> *make_casimir_gfunc_kz(double T, double dt, double sigma, field_type ft);

#if MEEP_SINGLE
// in mympi.cpp ... must be here in order to use realnum type
void broadcast(int from, realnum *data, int size);
#endif

// random number generation: random.cpp
void set_random_seed(unsigned long seed);
double uniform_random(double a, double b); // uniform random in [a,b]
double gaussian_random(double mean, double stddev); // normal random with given mean and stddev
int random_int(int a, int b); // uniform random in [a,b)

// Bessel function (in initialize.cpp)
double BesselJ(int m, double kr);

// analytical Green's functions (in near2far.cpp); upon return,
// EH[0..5] are set to the Ex,Ey,Ez,Hx,Hy,Hz field components at x
// from a c0 source of amplitude f0 at x0.
void green2d(std::complex<double> *EH, const vec &x,
             double freq, double eps, double mu,
             const vec &x0, component c0, std::complex<double> f0);
void green3d(std::complex<double> *EH, const vec &x,
             double freq, double eps, double mu,
             const vec &x0, component c0, std::complex<double> f0);

} /* namespace meep */

#endif /* MEEP_H */
