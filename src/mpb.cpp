/* Copyright (C) 2005-2009 Massachusetts Institute of Technology.
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

#include <stdlib.h>
#include "meep.hpp"
#include "config.h"

#ifdef HAVE_MPB
#  include <mpb.h>
#  ifndef SCALAR_COMPLEX
#    error Meep requires complex version of MPB
#  endif
#endif

namespace meep {

#ifdef HAVE_MPB

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
	vec(o[0] + r[0] * s[0], o[1] + r[1] * s[1], o[1] + r[1] * s[1]) :
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

static const complex<mpb_real> *meep_mpb_A_data = 0;
static const int *meep_mpb_A_n = 0;
static const double *meep_mpb_A_s = 0;
static int meep_mpb_A_component = 0;
static vec meep_mpb_A_center;
static complex<double> one(const vec &pt) {(void) pt; return 1.0;}
static complex<double> (*meep_mpb_A_A)(const vec &) = 0;
static complex<double> meep_mpb_A(const vec &p) {
  const complex<mpb_real> *data = meep_mpb_A_data + meep_mpb_A_component;
  int nx = meep_mpb_A_n[0];
  int ny = meep_mpb_A_n[1];
  int nz = meep_mpb_A_n[2];
  const double *s = meep_mpb_A_s;
  double r[3] = {0,0,0};
  vec p0(p - meep_mpb_A_center);
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

  return (complex<double>(double(real(ret)), double(imag(ret)))
	  * meep_mpb_A_A(p));
}

#endif /* HAVE_MPB */

void fields::add_eigenmode_source(component c0, const src_time &src,
				  direction d, const volume &where,
				  const volume &eig_vol,
				  int band_num, 
				  const vec &kpoint, bool match_frequency,
				  int parity,
				  double resolution, double eigensolver_tol,
				  complex<double> amp,
				  complex<double> A(const vec &)) {
#ifdef HAVE_MPB
  if (resolution <= 0) resolution = 2 * gv.a; // default to twice resolution
  int n[3], local_N, N_start, alloc_N, mesh_size[3] = {1,1,1};
  mpb_real k[3] = {0,0,0}, kcart[3] = {0,0,0};
  double s[3] = {0,0,0}, o[3] = {0,0,0};
  mpb_real R[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  mpb_real G[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  mpb_real kdir[3] = {0,0,0};
  double omega_src = real(src.frequency()), kscale = 1.0;
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
    break;
  case D1:
    o[2] = eig_vol.in_direction_min(Z);
    s[2] = eig_vol.in_direction(Z);
    k[2] = kpoint.in_direction(Z);
    break;
  default:
    abort("unsupported dimensionality in add_mpb_source");
  }
  for (int i = 0; i < 3; ++i) {
    n[i] = int(resolution * s[i] + 0.5); if (n[i] == 0) n[i] = 1;
    R[i][i] = s[i] = s[i] == 0 ? 1 : s[i];
    G[i][i] = 1 / R[i][i]; // recip. latt. vectors / 2 pi
  }
  
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      kcart[i] += G[j][i] * k[j];
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

  do {
    eigensolver(H, eigvals, maxwell_operator, (void *) mdata,
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
      mpb_real v, vscratch; // v = Re( W[0]* (dTheta/dk) W[0] ) = g. velocity
      evectmatrix_XtY_diag_real(W[0], W[1], &v, &vscratch);
      v /= sqrt(eigvals[band_num - 1]);

      // return to original size
      evectmatrix_resize(&W[0], band_num, 0);
      evectmatrix_resize(&W[1], band_num, 0);

      // update k via Newton step
      kscale = kscale - (sqrt(eigvals[band_num - 1]) - omega_src) 
	/ (v * abs(kpoint));
      if (!quiet)
	master_printf("Newton step: group velocity v=%g, kscale=%g\n",
		      v, kscale);
      if (kscale < 0 || kscale > 100)
	abort("Newton solver not converging -- need a better starting kpoint");
      for (int i = 0; i < 3; ++i) knew[i] = k[i] * kscale;
      update_maxwell_data_k(mdata, knew, G[0], G[1], G[2]);
    }
  } while (match_frequency 
	   && fabs(sqrt(eigvals[band_num - 1]) - omega_src) > 
	   omega_src * match_tol);
  
  evect_destroy_constraints(constraints);
  for (int i = 0; i < 3; ++i)
    destroy_evectmatrix(W[i]);
  
  src_time *src_mpb = src.clone();
  if (!match_frequency) src_mpb->set_frequency(sqrt(eigvals[band_num - 1]));
  
  complex<mpb_real> *cdata = (complex<mpb_real> *) mdata->fft_data;
  meep_mpb_A_s = s;
  meep_mpb_A_n = n;
  meep_mpb_A_data = cdata;
  meep_mpb_A_center = eig_vol.center() - where.center();
  meep_mpb_A_A = A ? A : one;

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

  if (is_D(c0)) c0 = direction_component(Ex, component_direction(c0));
  if (is_B(c0)) c0 = direction_component(Hx, component_direction(c0));

  // use principle of equivalence to obtain equivalent currents
  FOR_ELECTRIC_COMPONENTS(c) 
    if (gv.has_field(c) && (c0 == Centered || c0 == c)
	&& component_direction(c) != d
	&& (gv.dim != D2 || !(parity & (EVEN_Z_PARITY | ODD_Z_PARITY))
	    || ((parity & EVEN_Z_PARITY) && !is_tm(c))
	    || ((parity & ODD_Z_PARITY) && is_tm(c)))) {
      // E current source = d x (eigenmode H)
      if ((d + 1) % 3 == component_direction(c) % 3) {
	meep_mpb_A_component = (d + 2) % 3;
	add_volume_source(c, *src_mpb, where, meep_mpb_A, -amp);
      }
      else {
	meep_mpb_A_component = (d + 1) % 3;
	add_volume_source(c, *src_mpb, where, meep_mpb_A, amp);
      }
      }
  
  maxwell_compute_d_from_H(mdata, H, (scalar_complex*)cdata, band_num - 1, 1);
  { // d_from_H actually computes -omega*D (see mpb/src/maxwell/maxwell_op.c)
    double scale = -1.0 / omega_src;
    int N = mdata->fft_output_size * 3;
    for (int i = 0; i < N; ++i) cdata[i] *= scale;
  }
  maxwell_compute_e_from_d(mdata, (scalar_complex*)cdata, 1);
  // use principle of equivalence to obtain equivalent currents
  FOR_MAGNETIC_COMPONENTS(c) 
    if (gv.has_field(c) && (c0 == Centered || c0 == c)
	&& component_direction(c) != d
	&& (gv.dim != D2 || !(parity & (EVEN_Z_PARITY | ODD_Z_PARITY))
	    || ((parity & EVEN_Z_PARITY) && !is_tm(c))
	    || ((parity & ODD_Z_PARITY) && is_tm(c)))) {
      // H current source = - d x (eigenmode E)
      if ((d + 1) % 3 == component_direction(c) % 3) {
	meep_mpb_A_component = (d + 2) % 3;
	add_volume_source(c, *src_mpb, where, meep_mpb_A, amp);
	}
      else {
	meep_mpb_A_component = (d + 1) % 3;
	add_volume_source(c, *src_mpb, where, meep_mpb_A, -amp);
      }
    }
  
  delete src_mpb;
  destroy_evectmatrix(H);
  delete[] eigvals;
  destroy_maxwell_data(mdata);
#else /* !defined(HAVE_MPB) */
  abort("Meep must be configured/compiled with MPB for add_mpb_source");
#endif
}

} // namespace meep
