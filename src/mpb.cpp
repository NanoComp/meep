/* Copyright (C) 2005-2015 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "meep.hpp"
#include "config.h"

#ifdef HAVE_MPB
#  include <mpb.h>
#  ifndef SCALAR_COMPLEX
#    error Meep requires complex version of MPB
#  endif
#endif

using namespace std;

typedef complex<double> cdouble;

namespace meep {

/**************************************************************/
/* dummy versions of class methods for compiling without MPB  */
/**************************************************************/
#ifndef HAVE_MPB
void *fields::get_eigenmode(double &omega_src,
	     		    direction d, const volume &where,
			    const volume &eig_vol,
	    	            int band_num,
		            const vec &kpoint, bool match_frequency,
                            int parity,
                            double resolution, 
                            double eigensolver_tol) {

  (void) omega_src; (void) d; (void) where; (void) eig_vol; 
  (void) band_num;  (void) kpoint; (void) match_frequency; 
  (void) parity; (void) resolution; (void) eigensolver_tol;
  abort("Meep must be configured/compiled with MPB for get_eigenmode");
}

void fields::add_eigenmode_source(component c0, const src_time &src,
				  direction d, const volume &where,
				  const volume &eig_vol,
				  int band_num,
				  const vec &kpoint, bool match_frequency,
				  int parity,
				  double resolution, double eigensolver_tol,
				  complex<double> amp,
				  complex<double> A(const vec &)) {
  (void) c0; (void) src; (void) d; (void) where; (void) eig_vol; 
  (void) band_num;  (void) kpoint; (void) match_frequency; 
  (void) parity; (void) resolution; (void) eigensolver_tol;
  (void) amp; (void) A;
  abort("Meep must be configured/compiled with MPB for add_eigenmode_source");
}

cdouble fields::get_eigenmode_coefficient(dft_flux *flux,
                                          int num_freq,
                                          direction d,
                                          const volume &where,
                                          int band_num,
                                          kpoint_func k_func,
                                          void *k_func_data)
{ (void) flux; (void) num_freq; (void) d; (void) where;
  (void) band_num; (void) k_func; (void) k_func_data;
  abort("Meep must be configured/compiled with MPB for get_eigenmode_coefficient");
}



#else // HAVE_MPB

typedef struct {
  const double *s, *o;
  ndim dim;
  const fields *f;
} meep_mpb_eps_data;

static void meep_mpb_eps(symmetric_matrix *eps,
			 symmetric_matrix *eps_inv,
			 const mpb_real r[3],
			 void *eps_data_) {
  meep_mpb_eps_data *eps_data = (meep_mpb_eps_data *) eps_data_;
  const double *s = eps_data->s;
  const double *o = eps_data->o;
  vec p(eps_data->dim == D3 ?
	vec(o[0] + r[0] * s[0], o[1] + r[1] * s[1], o[2] + r[2] * s[2]) :
	(eps_data->dim == D2 ?
	 vec(o[0] + r[0] * s[0], o[1] + r[1] * s[1]) :
	 /* D1 */ vec(o[2] + r[2] * s[2])));
  const fields *f = eps_data->f;
  eps_inv->m00 = f->get_chi1inv(Ex, X, p);
  eps_inv->m11 = f->get_chi1inv(Ey, Y, p);
  eps_inv->m22 = f->get_chi1inv(Ez, Z, p);
  //  master_printf("eps_zz(%g,%g) = %g\n", p.x(), p.y(), 1/eps_inv->m00);
  ASSIGN_ESCALAR(eps_inv->m01, f->get_chi1inv(Ex, Y, p), 0);
  ASSIGN_ESCALAR(eps_inv->m02, f->get_chi1inv(Ex, Z, p), 0);
  ASSIGN_ESCALAR(eps_inv->m12, f->get_chi1inv(Ey, Z, p), 0);
  maxwell_sym_matrix_invert(eps, eps_inv);
}

/**************************************************************/
/* prototype for position-dependent amplitude function passed */
/* to add_volume_source                                       */
/**************************************************************/
typedef complex<double> (*amplitude_function)(const vec &);

// default implementation of amplitude_function
static complex<double> default_amp_func(const vec &pt) 
 {(void) pt; return 1.0;}

/*******************************************************************/
/* structure storing all data needed to compute position-dependent */
/* amplitude for eigenmode source (the fields of this structure    */
/* were formerly global variables)                                 */
/*******************************************************************/
typedef struct eigenmode_data
 { 
   maxwell_data *mdata;
   evectmatrix H;
   int n[3];
   int component;
   double s[3];
   vec center;
   amplitude_function amp_func;
   int band_num;
   double omega;
   double group_velocity;
 } eigenmode_data;

/*******************************************************************/
/* compute position-dependent amplitude for eigenmode source       */
/*  (similar to the routine formerly called meep_mpb_A)            */
/*******************************************************************/
static FILE *ESFieldFile=0; // "eigenmode-source log file"
complex<double> eigenmode_amplitude(const vec &p,
                                    eigenmode_data *edata) {
  
  if ( !edata || !(edata->mdata) )
   abort("%s:%i: internal error",__FILE__,__LINE__);
   
  maxwell_data *mdata         = edata->mdata;
  int *n                      = edata->n;
  int component               = edata->component % 3;
  double *s                   = edata->s;
  vec center                  = edata->center;
  amplitude_function amp_func = edata->amp_func;

  complex<mpb_real> *cdata = (complex<mpb_real> *) mdata->fft_data;
  const complex<mpb_real> *data = cdata + component;
  int nx = n[0];
  int ny = n[1];
  int nz = n[2];
  double r[3] = {0,0,0};
  vec p0(p - center);
  LOOP_OVER_DIRECTIONS(p.dim, d) r[d%3] = p0.in_direction(d) / s[d%3] + 0.5;
  double rx = r[0], ry = r[1], rz = r[2];

  /* linearly interpolate the amplitude from MPB at point p */

  int x, y, z, x2, y2, z2;
  double dx, dy, dz;

  /* get the point corresponding to r in the epsilon array grid: */
  x = int(rx * nx);
  y = int(ry * ny);
  z = int(rz * nz);

  /* get the difference between (x,y,z) and the actual point */
  dx = rx * nx - x;
  dy = ry * ny - y;
  dz = rz * nz - z;

  /* get the other closest point in the grid, with periodic boundaries: */
  x2 = (nx + (dx >= 0.0 ? x + 1 : x - 1)) % nx;
  y2 = (ny + (dy >= 0.0 ? y + 1 : y - 1)) % ny;
  z2 = (nz + (dz >= 0.0 ? z + 1 : z - 1)) % nz;
  x = x % nx; y = y % ny; z = z % nz;

  /* take abs(d{xyz}) to get weights for {xyz} and {xyz}2: */
  dx = fabs(dx);
  dy = fabs(dy);
  dz = fabs(dz);

  /* define a macro to give us data(x,y,z) on the grid,
     in row-major order (the order used by MPB): */
#define D(x,y,z) (data[(((x)*ny + (y))*nz + (z)) * 3])
  complex<mpb_real> ret;
  ret = (((D(x,y,z)*(1.0-dx) + D(x2,y,z)*dx) * (1.0-dy) +
	  (D(x,y2,z)*(1.0-dx) + D(x2,y2,z)*dx) * dy) * (1.0-dz) +
	 ((D(x,y,z2)*(1.0-dx) + D(x2,y,z2)*dx) * (1.0-dy) +
	  (D(x,y2,z2)*(1.0-dx) + D(x2,y2,z2)*dx) * dy) * dz);
#undef D


  cdouble amp = amp_func(p);
  if (ESFieldFile)
   fprintf(ESFieldFile,"%e %e %e %e %e %e %e\n",p.x(),p.y(),p.z(),real(ret),imag(ret),real(amp),imag(amp));

  return (complex<double>(double(real(ret)), double(imag(ret)))
	  * amp_func(p));
}

/***************************************************************/
/* entry point to eigenmode_amplitude with the right prototype */
/* for passage as the A parameter to add_volume_source         */
/***************************************************************/
static eigenmode_data *global_eigenmode_data=0;
static complex<double> meep_mpb_A(const vec &p)
 { return eigenmode_amplitude(p, global_eigenmode_data); }


/****************************************************************/
/* call MPB to get the band_numth eigenmode at freq omega_src.  */
/*                                                              */
/* this routine constitutes the first 75% of what was formerly  */
/* add_eigenmode_source; it has been split off as a separate    */
/* routine to allow it to be followed either by                 */
/*  (a) add_eigenmode_src()                                     */
/* or                                                           */
/*  (b) get_eigenmode_coefficient()                             */
/*                                                              */
/* the return value is an opaque pointer to an eigenmode_data   */
/* structure (needs to be opaque to allow compilation without   */
/* MPB, in which case maxwell_data and other types aren't       */
/* defined). this structure may then be passed to               */
/* eigenmode_amplitude (above) to compute values of eigenmode   */
/* magnetic field components at arbitrary points.               */
/*                                                              */
/* if want to compute *electric* field components of the        */
/* eigenmode, you can do that too, but you must first call      */
/* switch_eigenmode_data_to_electric_field() on the             */
/* eigenmode_data structure.                                    */
/****************************************************************/
void *fields::get_eigenmode(double &omega_src,
	     		    direction d, const volume &where,
			    const volume &eig_vol,
	    	            int band_num,
		            const vec &kpoint, bool match_frequency,
                            int parity,
                            double resolution, 
                            double eigensolver_tol) {

  /*--------------------------------------------------------------*/
  /*- part 1: preliminary setup for calling MPB  -----------------*/
  /*--------------------------------------------------------------*/

  if (resolution <= 0) resolution = 2 * gv.a; // default to twice resolution
  int n[3], local_N, N_start, alloc_N, mesh_size[3] = {1,1,1};
  mpb_real k[3] = {0,0,0}, kcart[3] = {0,0,0};
  double s[3] = {0,0,0}, o[3] = {0,0,0};
  mpb_real R[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  mpb_real G[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  mpb_real kdir[3] = {0,0,0};
  double kscale = 1.0;
  double match_tol = eigensolver_tol * 10;

  if (d == NO_DIRECTION || coordinate_mismatch(gv.dim, d))
    abort("invalid direction in add_eigenmode_source");
  if (where.dim != gv.dim || eig_vol.dim != gv.dim)
    abort("invalid volume dimensionality in add_eigenmode_source");

  if (!eig_vol.contains(where))
    abort("invalid grid_volume in add_eigenmode_source (WHERE must be in EIG_VOL)");

  switch (gv.dim) {
  case D3:
    o[0] = eig_vol.in_direction_min(X);
    o[1] = eig_vol.in_direction_min(Y);
    o[2] = eig_vol.in_direction_min(Z);
    s[0] = eig_vol.in_direction(X);
    s[1] = eig_vol.in_direction(Y);
    s[2] = eig_vol.in_direction(Z);
    k[0] = kpoint.in_direction(X);
    k[1] = kpoint.in_direction(Y);
    k[2] = kpoint.in_direction(Z);
    break;
  case D2:
    o[0] = eig_vol.in_direction_min(X);
    o[1] = eig_vol.in_direction_min(Y);
    s[0] = eig_vol.in_direction(X);
    s[1] = eig_vol.in_direction(Y);
    k[0] = kpoint.in_direction(X);
    k[1] = kpoint.in_direction(Y);
    // the following line was missing from the original mpb.cpp,
    // but I think it's needed! Consider a waveguide of 
    // constant (x,y) cross section with power flow in the z direction.
    k[2] = kpoint.in_direction(Z); 
    break;
  case D1:
    o[2] = eig_vol.in_direction_min(Z);
    s[2] = eig_vol.in_direction(Z);
    k[2] = kpoint.in_direction(Z);
    break;
  default:
    abort("unsupported dimensionality in add_eigenmode_source");
  }

  master_printf("KPOINT: %g, %g, %g\n", k[0], k[1], k[2]);

  // if match_frequency is true, all we need is a direction for k
  // and a crude guess for its value; we must supply this if k==0.
  if (match_frequency && k[0] == 0 && k[1] == 0 && k[2] == 0) {
    k[d-X] = omega_src * sqrt(get_eps(eig_vol.center()));
    master_printf("NEW KPOINT: %g, %g, %g\n", k[0], k[1], k[2]);
    if (s[d-X] > 0) {
      k[d-X] *= s[d-X]; // put k in G basis (inverted when we compute kcart)
      if (fabs(k[d-X]) > 0.4)  // ensure k is well inside the Brillouin zone
	k[d-X] = k[d-X] > 0 ? 0.4 : -0.4;
      master_printf("NEWER KPOINT: %g, %g, %g\n", k[0], k[1], k[2]);
    }
  }

  for (int i = 0; i < 3; ++i) {
    n[i] = int(resolution * s[i] + 0.5); if (n[i] == 0) n[i] = 1;
    R[i][i] = s[i] = s[i] == 0 ? 1 : s[i];
    G[i][i] = 1 / R[i][i]; // recip. latt. vectors / 2 pi
  }

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      kcart[i] += G[j][i] * k[j];
  double klen0 = sqrt(k[0]*k[0]+k[1]*k[1]+k[2]*k[2]);
  double klen = sqrt(kcart[0]*kcart[0]+kcart[1]*kcart[1]+kcart[2]*kcart[2]);
  if (klen == 0.0) {
    if (match_frequency) abort("need nonzero kpoint guess to match frequency");
    klen = 1;
  }
  kdir[0] = kcart[0] / klen;
  kdir[1] = kcart[1] / klen;
  kdir[2] = kcart[2] / klen;

  maxwell_data *mdata = create_maxwell_data(n[0], n[1], n[2],
					    &local_N, &N_start, &alloc_N,
					    band_num, band_num);
  if (local_N != n[0] * n[1] * n[2])
    abort("MPI version of MPB library not supported");

  set_maxwell_data_parity(mdata, parity);
  update_maxwell_data_k(mdata, k, G[0], G[1], G[2]);

  if (k[0] == 0 && k[1] == 0 && k[2] == 0) {
    evectmatrix H; H.p = band_num; H.c = 2;
    band_num -= maxwell_zero_k_num_const_bands(H, mdata);
    if (band_num == 0)
      abort("zero-frequency bands at k=0 are ill-defined");
  }

  meep_mpb_eps_data eps_data;
  eps_data.s = s; eps_data.o = o; eps_data.dim = gv.dim; eps_data.f = this;
  set_maxwell_dielectric(mdata, mesh_size, R, G, meep_mpb_eps,NULL, &eps_data);
  if (check_maxwell_dielectric(mdata, 0))
    abort("invalid dielectric function for MPB");

  evectmatrix H = create_evectmatrix(n[0] * n[1] * n[2], 2, band_num,
				     local_N, N_start, alloc_N);
  for (int i = 0; i < H.n * H.p; ++i) {
    ASSIGN_SCALAR(H.data[i], rand() * 1.0/RAND_MAX, rand() * 1.0/RAND_MAX);
  }

  mpb_real *eigvals = new mpb_real[band_num];
  int num_iters;
  evectmatrix W[3];
  for (int i = 0; i < 3; ++i)
    W[i] = create_evectmatrix(n[0] * n[1] * n[2], 2, band_num,
			      local_N, N_start, alloc_N);

  evectconstraint_chain *constraints = NULL;
  constraints = evect_add_constraint(constraints,
				     maxwell_parity_constraint,
				     (void *) mdata);
  if (k[0] == 0 && k[1] == 0 && k[2] == 0)
    constraints = evect_add_constraint(constraints,
				       maxwell_zero_k_constraint,
				       (void *) mdata);

  mpb_real knew[3]; for (int i = 0; i < 3; ++i) knew[i] = k[i];

  mpb_real vgrp; // Re( W[0]* (dTheta/dk) W[0] ) = group velocity

  /*--------------------------------------------------------------*/
  /*- part 2: newton iteration loop with call to MPB on each step */
  /*-         until eigenmode converged to requested tolerance    */
  /*--------------------------------------------------------------*/
  do {
    eigensolver(H, eigvals, maxwell_operator, (void *) mdata,
#if MPB_VERSION_MAJOR > 1 || (MPB_VERSION_MAJOR == 1 && MPB_VERSION_MINOR >= 6)
                NULL, NULL, /* eventually, we can support mu here */
#endif
		maxwell_preconditioner2, (void *) mdata,
		evectconstraint_chain_func,
		(void *) constraints,
		W, 3,
		eigensolver_tol, &num_iters,
		EIGS_DEFAULT_FLAGS |
		(am_master() && !quiet ? EIGS_VERBOSE : 0));
    if (!quiet)
      master_printf("MPB solved for omega_%d(%g,%g,%g) = %g after %d iters\n",
		    band_num, knew[0],knew[1],knew[2],
		    sqrt(eigvals[band_num-1]), num_iters);

    if (match_frequency) {
      // copy desired single eigenvector into scratch arrays
      evectmatrix_resize(&W[0], 1, 0);
      evectmatrix_resize(&W[1], 1, 0);
      for (int i = 0; i < H.n; ++i)
	W[0].data[i] = H.data[H.p-1 + i * H.p];

      // compute the group velocity in the k direction
      maxwell_ucross_op(W[0], W[1], mdata, kdir); // W[1] = (dTheta/dk) W[0]
      mpb_real vscratch; 
      evectmatrix_XtY_diag_real(W[0], W[1], &vgrp, &vscratch);
      vgrp /= sqrt(eigvals[band_num - 1]);

      // return to original size
      evectmatrix_resize(&W[0], band_num, 0);
      evectmatrix_resize(&W[1], band_num, 0);

      // update k via Newton step
      kscale = kscale - (sqrt(eigvals[band_num - 1]) - omega_src) / (vgrp*klen0);
      if (!quiet)
	master_printf("Newton step: group velocity v=%g, kscale=%g\n", vgrp, kscale);
      if (kscale < 0 || kscale > 100)
	abort("Newton solver not converging -- need a better starting kpoint");
      for (int i = 0; i < 3; ++i) knew[i] = k[i] * kscale;
      update_maxwell_data_k(mdata, knew, G[0], G[1], G[2]);
    }
  } while (match_frequency
	   && fabs(sqrt(eigvals[band_num - 1]) - omega_src) >
	   omega_src * match_tol);

  if (!match_frequency)
   omega_src = sqrt(eigvals[band_num - 1]);

  // cleanup temporary storage
  delete[] eigvals;
  evect_destroy_constraints(constraints);
  for (int i = 0; i < 3; ++i)
    destroy_evectmatrix(W[i]);

  /*--------------------------------------------------------------*/
  /*- part 3: do one stage of postprocessing to tabulate H-field  */
  /*-         components on the internal storage buffer in mdata  */
  /*-         (a similar step is performed subsequently in the    */
  /*-          calling routine to replace H-fields with E-fields  */
  /*--------------------------------------------------------------*/
  complex<mpb_real> *cdata = (complex<mpb_real> *) mdata->fft_data;

  maxwell_compute_h_from_H(mdata, H, (scalar_complex*)cdata, band_num - 1, 1);
  /* choose deterministic phase, maximizing power in real part;
     see fix_field_phase routine in MPB.*/
  {
    int i, N = mdata->fft_output_size * 3;
    double sq_sum0 = 0, sq_sum1 = 0, maxabs = 0.0;
    double theta;
    for (i = 0; i < N; ++i) {
      double a = real(cdata[i]), b = imag(cdata[i]);
      sq_sum0 += a*a - b*b;
      sq_sum1 += 2*a*b;
    }
    theta = 0.5 * atan2(-sq_sum1, sq_sum0);
    complex<mpb_real> phase(cos(theta), sin(theta));
    for (i = 0; i < N; ++i) {
      double r = fabs(real(cdata[i] * phase));
      if (r > maxabs) maxabs = r;
    }
    for (i = N-1; i >= 0 && fabs(real(cdata[i] * phase)) < 0.5 * maxabs; --i)
      ;
    if (real(cdata[i] * phase) < 0) phase = -phase;
    for (i = 0; i < N; ++i) cdata[i] *= phase;
    complex<mpb_real> *hdata = (complex<mpb_real> *) H.data;
    for (i = 0; i < H.n; ++i) hdata[i*H.p + (band_num-1)] *= phase;
  }

  /*--------------------------------------------------------------*/
  /*- part 4: initialize and return output data structures.       */
  /*--------------------------------------------------------------*/
  eigenmode_data *edata = new eigenmode_data;
  edata->mdata          = mdata;
  edata->H              = H;
  edata->n[0]           = n[0];
  edata->n[1]           = n[1];
  edata->n[2]           = n[2];
  edata->component      = Hx % 3;
  edata->s[0]           = s[0];
  edata->s[1]           = s[1];
  edata->s[2]           = s[2];
  edata->center         = eig_vol.center() - where.center();
  edata->amp_func       = default_amp_func;
  edata->band_num       = band_num;
  edata->omega          = omega_src;
  edata->group_velocity = (double) vgrp;
  return (void *)edata;
}

// the eigenmode_data structure returned by get_eigenmode 
// initially stores H-field data on its internal real-space grid;
// this routine switches that to E-field data.
void switch_eigenmode_data_to_electric_field(eigenmode_data *edata)
{
  maxwell_data *mdata      = edata->mdata;
  complex<mpb_real> *cdata = (complex<mpb_real> *)mdata->fft_data;
  evectmatrix H            = edata->H;
  int band_num             = edata->band_num;
  double omega             = edata->omega;

  maxwell_compute_d_from_H(mdata, H, (scalar_complex*)cdata, band_num - 1, 1);
  // d_from_H actually computes -omega*D (see mpb/src/maxwell/maxwell_op.c)
  double scale = -1.0 / omega;
  int N = mdata->fft_output_size * 3;
  for (int i = 0; i < N; ++i) 
   cdata[i] *= scale;

  maxwell_compute_e_from_d(mdata, (scalar_complex*)cdata, 1);
}

void destroy_eigenmode_data(eigenmode_data *edata)
{ 
  destroy_evectmatrix( edata->H  );
  destroy_maxwell_data( edata->mdata ); 
  delete edata;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
FILE *vfopen(const char *format, const char *mode, ...)
{
  va_list ap;
  char buffer[1000];
  va_start(ap,mode);
  vsnprintf(buffer,1000,format,ap);
  va_end(ap);

  FILE *f=fopen(buffer,mode);
  return f;
}

/***************************************************************/
/* call get_eigenmode() to solve for the specified eigenmode,  */
/* then call add_volume_source() to add current sources whose  */
/* radiated fields reproduce the eigenmode fields              */
/***************************************************************/
void fields::add_eigenmode_source(component c0, const src_time &src,
				  direction d, const volume &where,
				  const volume &eig_vol,
				  int band_num,
				  const vec &kpoint, bool match_frequency,
				  int parity,
				  double resolution, double eigensolver_tol,
				  complex<double> amp,
				  complex<double> A(const vec &)) {


  char *s=getenv("MEEP_PLOT_EIGENMODE_SOURCE_FIELDS");
  bool PlotESFields = (s && s[0]=='1');

  /*--------------------------------------------------------------*/
  /* step 1: call MPB to compute the eigenmode                    */
  /*--------------------------------------------------------------*/
  double omega_src = real(src.frequency());
  global_eigenmode_data
   =(eigenmode_data *)get_eigenmode(omega_src, d, where,
                                    eig_vol, band_num,
                                    kpoint, match_frequency,
                                    parity, resolution, 
                                    eigensolver_tol);

  global_eigenmode_data->amp_func = A ? A : default_amp_func;
  
  src_time *src_mpb = src.clone();
  if (!match_frequency)
    src_mpb->set_frequency(omega_src);

  if (is_D(c0)) c0 = direction_component(Ex, component_direction(c0));
  if (is_B(c0)) c0 = direction_component(Hx, component_direction(c0));

  /*--------------------------------------------------------------*/
  // step 2: add sources whose radiated field reproduces the      */
  //         the eigenmode                                        */
  /*--------------------------------------------------------------*/

  // step 2a: electric-current sources
  //           = nHat \times magnetic-field components

  // use principle of equivalence to obtain equivalent currents
  FOR_ELECTRIC_COMPONENTS(c)
    if (gv.has_field(c) && (c0 == Centered || c0 == c)
	&& component_direction(c) != d
	&& (gv.dim != D2 || !(parity & (EVEN_Z_PARITY | ODD_Z_PARITY))
	    || ((parity & EVEN_Z_PARITY) && !is_tm(c))
	    || ((parity & ODD_Z_PARITY) && is_tm(c)))) {
      // E current source = d x (eigenmode H)
      if ((d + 1) % 3 == component_direction(c) % 3) {
	global_eigenmode_data->component = (d + 2) % 3;
	if (PlotESFields) ESFieldFile=vfopen("/tmp/ESLog%i_%i_h%c_%g","w",my_rank(),band_num,'x' + global_eigenmode_data->component, omega_src);
	add_volume_source(c, *src_mpb, where, meep_mpb_A, -amp);
	if (ESFieldFile) fclose(ESFieldFile);
      }
      else {
	global_eigenmode_data->component= (d + 1) % 3;
	if (PlotESFields) ESFieldFile=vfopen("/tmp/ESLog%i_%i_h%c_%g","w",my_rank(),band_num,'x' + global_eigenmode_data->component, omega_src);
	add_volume_source(c, *src_mpb, where, meep_mpb_A, amp);
	if (ESFieldFile) fclose(ESFieldFile);
      }
      }

  // step 2b: post-processing step to replace H-field components
  //          with E-field components in the internal data buffer
  //          inside mdata; cf. Part 3 of get_eigenmode() above
  switch_eigenmode_data_to_electric_field(global_eigenmode_data);

  // step 2c: magnetic-current sources
  //           = - nHat \times electric-field components

  // use principle of equivalence to obtain equivalent currents
  FOR_MAGNETIC_COMPONENTS(c)
    if (gv.has_field(c) && (c0 == Centered || c0 == c)
	&& component_direction(c) != d
	&& (gv.dim != D2 || !(parity & (EVEN_Z_PARITY | ODD_Z_PARITY))
	    || ((parity & EVEN_Z_PARITY) && !is_tm(c))
	    || ((parity & ODD_Z_PARITY) && is_tm(c)))) {
      // H current source = - d x (eigenmode E)
      if ((d + 1) % 3 == component_direction(c) % 3) {
	global_eigenmode_data->component= (d + 2) % 3;
	if (PlotESFields) ESFieldFile=vfopen("/tmp/ESLog%i_%i_e%c_%g","w",my_rank(),band_num,'x' + global_eigenmode_data->component, omega_src);
	add_volume_source(c, *src_mpb, where, meep_mpb_A, amp);
	if (ESFieldFile) fclose(ESFieldFile);
	}
      else {
	global_eigenmode_data->component = (d + 1) % 3;
	if (PlotESFields) ESFieldFile=vfopen("/tmp/ESLog%i_%i_e%c_%g","w",my_rank(),band_num,'x' + global_eigenmode_data->component, omega_src);
	add_volume_source(c, *src_mpb, where, meep_mpb_A, -amp);
	if (ESFieldFile) fclose(ESFieldFile);
      }
    }

  delete src_mpb;
  destroy_eigenmode_data(global_eigenmode_data);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
void flux_output_hdf5(fields *f, dft_flux *flux, direction d,
                      int num_freq, char *file_base)
{
  //double omega = flux->freq_min + num_freq*flux->dfreq;

  h5file *file = open_h5file(dataname, h5file::WRITE, prefix, true);
  file->create_or_extend_data(dataname, rank, dims,
                              append_data, single_precision);

  /*--------------------------------------------------------------*/
  /*- this loop amounts to a "loop_in_dft_chunks()" function and  */
  /*- should maybe be promoted to a standalone function?          */
  /*--------------------------------------------------------------*/
  int Nfreq          = flux->Nfreq;
  for ( dft_chunk *E=flux->E, *H=flux->H; E && H;
        E=E->next_in_dft, H=H->next_in_dft
      )
   { 
     // create output files for E and H components
     char E_file_name[100];
     snprintf(E_file_name,100,"%s_%s.dat",file_base,component_name(E->c));
     FILE *E_file=fopen(E_file_name,"a");

     char H_file_name[100];
     snprintf(H_file_name,100,"%s_%s.dat",file_base,component_name(H->c));
     FILE *H_file=fopen(H_file_name,"a");

     // extract info from the current dft_chunk
     fields_chunk *fc = E->fc;
     ivec is          = E->is;
     ivec ie          = E->ie;
     vec s0           = E->s0;
     vec s1           = E->s1;
     vec e0           = E->e0;
     vec e1           = E->e1;
     double dV0       = E->dV0;
     double dV1       = E->dV1;
     ivec shift       = E->shift;
     symmetry S       = E->S;
     int sn           = E->sn;

     vec rshift(shift * (0.5*fc->gv.inva));

     // loop over all points in the current dft_chunk
     int chunk_idx = 0;
     LOOP_OVER_IVECS(fc->gv, is, ie, idx)
      { 
        // get the coordinates and integration weight for this grid point
        IVEC_LOOP_LOC(fc->gv, loc);
        loc = S.transform(loc, sn) + rshift;
        double w=IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2);

        // get the E and H field components at this grid point,
        //  dividing out any extra weight factors it may already contain
        cdouble E_flux = E->dft[ Nfreq*(chunk_idx++) + num_freq];
        if (E->include_dV_and_interp_weights)
         E_flux /= (E->sqrt_dV_and_interp_weights ? sqrt(w) : w);

        cdouble H_flux = H->dft[ Nfreq*(chunk_idx++) + num_freq];
        if (H->include_dV_and_interp_weights)
         H_flux /= (H->sqrt_dV_and_interp_weights ? sqrt(w) : w);

        // kinda byzantine: the second of the two E-field components
        // is stored with a minus sign, which we want to remove for 
        // our purposes; but *which* component is the second component 
        // depends on the direction 'd' that was used to create the dft_flux
        if (     (d==X && E->c==Ez)
             ||  (d==Y && E->c==Ex)
             ||  (d==R && E->c==Ez)
             ||  (d==P && E->c==Er)
             ||  (d==Z && f->gv.dim == Dcyl && E->c==Ep)
             ||  (d==Z && f->gv.dim != Dcyl && E->c==Ey)
           ) E_flux*= -1.0;
   
        fprintf(E_file,"%e %e %e %e %e\n",loc.x(), loc.y(), loc.z(), real(E_flux), imag(E_flux));
        fprintf(H_file,"%e %e %e %e %e\n",loc.x(), loc.y(), loc.z(), real(H_flux), imag(H_flux));

      }; // LOOP_OVER_IVECS

     fclose(E_file);
     fclose(H_file);

   }; // for ( dft_chunk *E=Echunks, *H=Hchunks ...

}
#endif

void output_flux(fields *f, dft_flux *flux, direction d,
                 int num_freq, char *file_base)
{
  //double omega = flux->freq_min + num_freq*flux->dfreq;

  /*--------------------------------------------------------------*/
  /*- this loop amounts to a "loop_in_dft_chunks()" function and  */
  /*- should maybe be promoted to a standalone function?          */
  /*--------------------------------------------------------------*/
  int Nfreq = flux->Nfreq;
  dft_chunk *EHList[2];
  EHList[0] = flux->E;
  EHList[1] = flux->H;
  for (int eh=0; eh<2; eh++)
   for (dft_chunk *EH=EHList[eh]; EH; EH=EH->next_in_dft)
    { 
     // create output files for E and H components
     char file_name[100];
     snprintf(file_name,100,"%s_%s.dat",file_base,component_name(EH->c));
     FILE *file=fopen(file_name,"a");

     // extract info from the current dft_chunk
     fields_chunk *fc = EH->fc;
     ivec is          = EH->is;
     ivec ie          = EH->ie;
     vec s0           = EH->s0;
     vec s1           = EH->s1;
     vec e0           = EH->e0;
     vec e1           = EH->e1;
     double dV0       = EH->dV0;
     double dV1       = EH->dV1;
     ivec shift       = EH->shift;
     symmetry S       = EH->S;
     int sn           = EH->sn;

     vec rshift(shift * (0.5*fc->gv.inva));

     // loop over all points in the current dft_chunk
     int chunk_idx = 0;
     LOOP_OVER_IVECS(fc->gv, is, ie, idx)
      { 
        // get the coordinates and integration weight for this grid point
        IVEC_LOOP_LOC(fc->gv, loc);
        loc = S.transform(loc, sn) + rshift;
        double w=IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2);

        // get the E and H field components at this grid point,
        //  dividing out any extra weight factors it may already contain
        cdouble flux = EH->dft[ Nfreq*(chunk_idx++) + num_freq];
        if (EH->include_dV_and_interp_weights)
         flux /= (EH->sqrt_dV_and_interp_weights ? sqrt(w) : w);

        // kinda byzantine: the second of the two E-field components
        // is stored with a minus sign, which we want to remove for 
        // our purposes; but *which* component is the second component 
        // depends on the direction 'd' that was used to create the dft_flux
        if (     (d==X && EH->c==Ez)
             ||  (d==Y && EH->c==Ex)
             ||  (d==R && EH->c==Ez)
             ||  (d==P && EH->c==Er)
             ||  (d==Z && f->gv.dim == Dcyl && EH->c==Ep)
             ||  (d==Z && f->gv.dim != Dcyl && EH->c==Ey)
           ) flux*= -1.0;
   
        fprintf(file,"%e %e %e %e %e\n",loc.x(), loc.y(), loc.z(), real(flux), imag(flux));

      }; // LOOP_OVER_IVECS
     fclose(file);

   }; // for (eh=0..1) for(dft_chunk *EH=EHList[eh]; EH; EH=EH->next_in_dft)

}

/***************************************************************/
/* call get_eigenmode() to solve for the specified eigenmode,  */
/* then call add_overlap_integral_contribution() multiple times*/
/* to sum all contributions to the numerator and denominator   */
/* of the eigenmode expansion coefficients.                    */
/***************************************************************/
cdouble fields::get_eigenmode_coefficient(dft_flux *flux,
                                          int num_freq,
                                          direction d,
                                          const volume &where,
                                          int band_num,
                                          kpoint_func k_func, void *k_func_data)
{
  char file_base[100];
  snprintf(file_base,100,"np%i_nb%i_nf%i",my_rank(),band_num,num_freq);
  output_flux(this, flux, d, num_freq, file_base);

  master_printf("Getting eigenmode coefficient (%i,%i)\n",num_freq, band_num);
#if 0
  /*--------------------------------------------------------------*/
  /* step 1: call MPB to compute the eigenmode                   -*/
  /*--------------------------------------------------------------*/
  double omega = flux->freq_min + num_freq*flux->dfreq;
  // call user's kpoint function if present
  vec kpoint(0.0, 0.0, 0.5); // TODO better default? 
  if (k_func) 
   kpoint=k_func(k_func_data, omega, band_num);

  bool match_frequency=true;
  int parity=0; 
  double resolution=a;
  double eigensolver_tol=1.0e-7;
  eigenmode_data *edata
   =(eigenmode_data *)get_eigenmode(omega, d, where, where,
                                    band_num, kpoint, match_frequency,
                                    parity, resolution, 
                                    eigensolver_tol);

  /*--------------------------------------------------------------*/
  /* step 2: sum contributions of all 4 surface-current cmpnents  */
  /*         to numerator and denominator of overlap integral     */
  /* num   = <caller's field | eigenmode>                         */
  /* denom = <eigenmode      | eigenmode>                         */
  /*--------------------------------------------------------------*/
  // step 2a: electric-current components 
  //            = nHat \times magnetic-field components
  cdouble numdenom[2]={0.0,0.0};

  FOR_ELECTRIC_COMPONENTS(c)
   {  
     if ( !(gv.has_field(c)) ) continue;
     // TODO restore parity check

     if ( (d+1)%3 == component_direction(c)%3 )
      { edata->component = (d+2)%3;
        add_overlap_integral_contribution(this, flux, d, num_freq, c, edata, -1.0, numdenom);
      }
     else if ( (d+2)%3 == component_direction(c)%3 )
      { edata->component = (d+1)%3;
        add_overlap_integral_contribution(this, flux, d, num_freq, c, edata, +1.0, numdenom);
      }
   };

  // step 2b: post-processing step to replace H-field components
  //          with E-field components in the internal data buffer
  //          inside mdata; cf. Part 3 of get_eigenmode() above
  switch_eigenmode_data_to_electric_field(edata);

  // step 2c: magnetic-current components 
  //            = -nHat \times electric-field components
  FOR_MAGNETIC_COMPONENTS(c)
   { 
     if ( !(gv.has_field(c)) ) continue;
     // TODO restore parity check

     if ( (d+1)%3 == component_direction(c)%3 )
      { edata->component = (d+2)%3;
        add_overlap_integral_contribution(this, flux, d, num_freq, c, edata, +1.0, numdenom);
      }
     else if ( (d+2)%3 == component_direction(c)%3 )
      { edata->component = (d+1)%3;
        add_overlap_integral_contribution(this, flux, d, num_freq, c, edata, -1.0, numdenom);
      }
   };

  destroy_eigenmode_data(edata);

  cdouble num=numdenom[0], denom=numdenom[1];
  if( denom==0.0 )
   { master_printf("**warning: denominator in get_eigenmode_coefficient**");
    return 0.0;
   };
  return num/denom;
#endif
return 0.0;
}

/***************************************************************/
/* get eigenmode coefficients for all frequencies in flux      */
/* and all band indices in the caller-populated bands array.   */
/*                                                             */
/* the array returned has length num_freqs x num_bands, with   */
/* the coefficient for frequency #nf, band #nb stored in slot  */ 
/* [ nb*num_freqs + nf ]                                       */
/***************************************************************/
std::vector<cdouble>
 fields::get_eigenmode_coefficients(dft_flux *flux, direction d,
                                    const volume &where,
                                    std::vector<int> bands,
                                    kpoint_func k_func, 
                                    void *k_func_data)
{ 
  int num_freqs = flux->Nfreq;
  int num_bands = bands.size();
  std::vector<cdouble> coeffs( num_freqs * num_bands );

  for(int nb=0; nb<num_bands; nb++)
   for(int nf=0; nf<num_freqs; nf++)
    coeffs[ nb*num_freqs + nf ] 
     = get_eigenmode_coefficient(flux, nf, d, where, bands[nb],
                                 k_func, k_func_data);

  return coeffs;
}

#endif // HAVE_MPB

} // namespace meep
