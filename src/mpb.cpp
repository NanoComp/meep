/* Copyright (C) 2005-2008 Massachusetts Institute of Technology.
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

#include "meep.hpp"
#include "config.h"

#ifdef HAVE_MPB
#define real mpb_real // avoid C++ conflict
#define SCALAR_COMPLEX 1 // complex version of MPB (not mpbi)
#include <mpb/maxwell.h>
#include <mpb/eigensolver.h>
#undef real
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
			 const double r[3],
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
  eps_inv->m01 = f->get_chi1inv(Ex, Y, p);
  eps_inv->m02 = f->get_chi1inv(Ex, Z, p);
  eps_inv->m12 = f->get_chi1inv(Ey, Z, p);
  maxwell_sym_matrix_invert(eps, eps_inv);
}

static const complex<double> *meep_mpb_A_data = 0;
static const int *meep_mpb_A_n = 0;
static const double *meep_mpb_A_s = 0;
static int meep_mpb_A_component = 0;
static vec meep_mpb_A_center;
static complex<double> meep_mpb_A(const vec &p) {
  const complex<double> *data = meep_mpb_A_data + meep_mpb_A_component;
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
  
  /* take abs(d{xyz}) to get weights for {xyz} and {xyz}2: */
  dx = fabs(dx);
  dy = fabs(dy);
  dz = fabs(dz);
  
  /* define a macro to give us data(x,y,z) on the grid,
     in row-major order (the order used by MPB): */
#define D(x,y,z) (data[(((x)*ny + (y))*nz + (z)) * 3])

  return(((D(x,y,z)*(1.0-dx) + D(x2,y,z)*dx) * (1.0-dy) +
	  (D(x,y2,z)*(1.0-dx) + D(x2,y2,z)*dx) * dy) * (1.0-dz) +
	 ((D(x,y,z2)*(1.0-dx) + D(x2,y,z2)*dx) * (1.0-dy) +
	  (D(x,y2,z2)*(1.0-dx) + D(x2,y2,z2)*dx) * dy) * dz);
  
#undef D
}

#endif /* HAVE_MPB */

void fields::add_eigenmode_source(const src_time &src,
				  const geometric_volume &where,
				  const geometric_volume &eig_vol,
				  int band_num, const vec &kpoint, int parity,
				  complex<double> amp,
				  double resolution,
				  double eigensolver_tol) {
#ifdef HAVE_MPB
  if (resolution <= 0) resolution = v.a;
  int n[3], local_N, N_start, alloc_N, mesh_size[3] = {1,1,1};
  double k[3] = {0,0,0}, s[3] = {0,0,0}, o[3] = {0,0,0};
  double R[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  double G[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

  if (!eig_vol.contains(where))
    abort("invalid volume in add_eigenmode_source (WHERE must be in EIG_VOL)");

  switch (v.dim) {
  case D3:
    o[0] = eig_vol.in_direction_min(X);
    o[1] = eig_vol.in_direction_min(Y);
    o[2] = eig_vol.in_direction_min(Z);
    s[0] = eig_vol.in_direction(X);
    s[1] = eig_vol.in_direction(X);
    s[2] = eig_vol.in_direction(X);
    k[0] = kpoint.in_direction(X);
    k[1] = kpoint.in_direction(Y);
    k[2] = kpoint.in_direction(Z);
    break;
  case D2:
    o[0] = eig_vol.in_direction_min(X);
    o[1] = eig_vol.in_direction_min(Y);
    s[0] = eig_vol.in_direction(X);
    s[1] = eig_vol.in_direction(X);
    k[0] = kpoint.in_direction(X);
    k[1] = kpoint.in_direction(Y);
    break;
  case D1:
    o[2] = eig_vol.in_direction_min(Z);
    s[2] = eig_vol.in_direction(X);
    k[2] = kpoint.in_direction(Z);
    break;
  default:
    abort("unsupported dimensionality in add_mpb_source");
  }
  for (int i = 0; i < 3; ++i) {
    n[i] = int(resolution * s[i] + 0.5);
    R[i][i] = s[i] = s[i] == 0 ? 1 : s[i];
    G[i][i] = 1 / R[i][i]; // recip. latt. vectors / 2 pi
  }
  
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
  eps_data.s = s; eps_data.o = o; eps_data.dim = v.dim; eps_data.f = this;
  set_maxwell_dielectric(mdata, mesh_size, R, G, meep_mpb_eps,NULL, &eps_data);
  if (check_maxwell_dielectric(mdata, 0))
    abort("invalid dielectric function for MPB");
  
  evectmatrix H = create_evectmatrix(n[0] * n[1] * n[2], 2, band_num,
				     local_N, N_start, alloc_N);
  
  double *eigvals = new double[band_num];
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
  
  eigensolver(H, eigvals, maxwell_operator, (void *) mdata,
	      maxwell_preconditioner2, (void *) mdata,
	      evectconstraint_chain_func,
	      (void *) constraints,
	      W, 3, 
	      eigensolver_tol, &num_iters, EIGS_DEFAULT_FLAGS);
  
  evect_destroy_constraints(constraints);
  for (int i = 0; i < 3; ++i)
    destroy_evectmatrix(W[i]);
  
  src_time *src_mpb = src.clone();
  src_mpb->set_frequency(sqrt(eigvals[band_num - 1]));
  
  complex<double> *cdata = (complex<double> *) mdata->fft_data;
  meep_mpb_A_s = s;
  meep_mpb_A_n = n;
  meep_mpb_A_data = cdata;
  meep_mpb_A_center = eig_vol.center() - where.center();
  
  maxwell_compute_h_from_H(mdata, H, (scalar_complex*)cdata, band_num - 1, 1);
  FOR_MAGNETIC_COMPONENTS(c) 
    if (v.has_field(c) && where.in_direction(component_direction(c))>0) {
      meep_mpb_A_component = component_direction(c) % 3;
      add_volume_source(c, *src_mpb, where, meep_mpb_A, amp);
    }
  
  maxwell_compute_d_from_H(mdata, H, (scalar_complex*)cdata, band_num - 1, 1);
  maxwell_compute_e_from_d(mdata, (scalar_complex*)cdata, 1);
  FOR_ELECTRIC_COMPONENTS(c) 
    if (v.has_field(c) && where.in_direction(component_direction(c))>0) {
      meep_mpb_A_component = component_direction(c) % 3;
      add_volume_source(c, *src_mpb, where, meep_mpb_A, amp);
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
