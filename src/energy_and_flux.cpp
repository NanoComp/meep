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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "meep.h"
#include "meep_internals.h"

namespace meep {

/* Energy calculation */

double fields::count_volume(component c) {
  double vol = 0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      vol += chunks[i]->count_volume(c);
  return sum_to_all(vol);
}

double fields_chunk::count_volume(component c) {
  double vol = 0;
  for (int i=0;i<v.ntot();i++)
    vol += v.dV(c,i).full_volume();
  return vol;
}

double fields::total_energy() {
  return energy_in_box(user_volume.surroundings());
}

double fields::field_energy() {
  return field_energy_in_box(user_volume.surroundings());
}

double fields::energy_in_box(const geometric_volume &where) {
  return thermo_energy_in_box(where) + field_energy_in_box(where);
}

double fields::field_energy_in_box(const geometric_volume &where) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->backup_h();
  step_h();
  step_boundaries(H_stuff);
  step_h_source();
  double next_step_magnetic_energy = magnetic_energy_in_box(where);
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->restore_h();

  return electric_energy_in_box(where) +
    0.5*next_step_magnetic_energy + 0.5*magnetic_energy_in_box(where);
}

struct dot_integrand_data {
  component c1, c2;
  long double sum;
};

static void dot_integrand(fields_chunk *fc, component cgrid,
			  ivec is, ivec ie,
			  vec s0, vec s1, vec e0, vec e1,
			  double dV0, double dV1,
			  ivec shift, complex<double> shift_phase,
			  const symmetry &S, int sn,
			  void *data_) {
  dot_integrand_data *data = (dot_integrand_data *) data_;

  (void) shift; // unused
  (void) shift_phase; // unused

  component c1 = S.transform(data->c1, -sn);
  component c2 = S.transform(data->c2, -sn);

  // We're integrating c1 * c2*, and we assume that
  // S.phase_shift(c1,sn) == S.phase_shift(c2,sn), so the phases all cancel

  // We also assume that cgrid is the same yee grid as that for c1 and c2,
  // so that no averaging is needed.
  (void) cgrid;

  for (int cmp = 0; cmp < 2; ++cmp)
    if (fc->f[c1][cmp] && fc->f[c2][cmp])
      LOOP_OVER_IVECS(fc->v, is, ie, idx)
	data->sum += IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2)
	  * fc->f[c1][cmp][idx] * fc->f[c2][cmp][idx];
}

double fields::field_energy_in_box(component c,
				   const geometric_volume &where) {
  if (coordinate_mismatch(v.dim, c))
    return 0.0;

  dot_integrand_data data;
  if (is_electric(c) || is_D(c)) {
    data.c1 = direction_component(Ex, component_direction(c));
    data.c2 = direction_component(Dx, component_direction(c));
  }
  else if (is_magnetic(c)) {
    data.c1 = data.c2 = direction_component(Hx, component_direction(c));
  }
  else
    abort("invalid field component in field_energy_in_box");

  data.sum = 0.0;
  integrate(dot_integrand, (void *) &data, where, c);
  return sum_to_all(data.sum) / (8*pi);
}

double fields::electric_energy_in_box(const geometric_volume &where) {
  long double sum = 0.0;
  FOR_ELECTRIC_COMPONENTS(c)
    sum += field_energy_in_box(c, where);
  return sum;
}

double fields::magnetic_energy_in_box(const geometric_volume &where) {
  long double sum = 0.0;
  FOR_MAGNETIC_COMPONENTS(c)
    sum += field_energy_in_box(c, where);
  return sum;
}

void fields_chunk::backup_h() {
  DOCMP FOR_MAGNETIC_COMPONENTS(c)
    if (f[c][cmp]) {
      if (f_backup[c][cmp] == NULL)
	f_backup[c][cmp] = new double[v.ntot()];
      if (f_backup_m_pml[c][cmp] == NULL)
	f_backup_m_pml[c][cmp] = new double[v.ntot()];
      if (f_backup_p_pml[c][cmp] == NULL)
	f_backup_p_pml[c][cmp] = new double[v.ntot()];
    }
  DOCMP FOR_MAGNETIC_COMPONENTS(c)
    if (f[c][cmp]) {
      memcpy(f_backup[c][cmp], f[c][cmp], v.ntot()*sizeof(double));
      memcpy(f_backup_p_pml[c][cmp], f_p_pml[c][cmp], v.ntot()*sizeof(double));
      memcpy(f_backup_m_pml[c][cmp], f_m_pml[c][cmp], v.ntot()*sizeof(double));
    }
}
  
void fields_chunk::restore_h() {
  DOCMP FOR_MAGNETIC_COMPONENTS(c)
    if (f[c][cmp]) {
      memcpy(f[c][cmp], f_backup[c][cmp], v.ntot()*sizeof(double));
      memcpy(f_p_pml[c][cmp], f_backup_p_pml[c][cmp], v.ntot()*sizeof(double));
      memcpy(f_m_pml[c][cmp], f_backup_m_pml[c][cmp], v.ntot()*sizeof(double));
    }
}

static void thermo_integrand(fields_chunk *fc, component cgrid,
			     ivec is, ivec ie,
			     vec s0, vec s1, vec e0, vec e1,
			     double dV0, double dV1,
			     ivec shift, complex<double> shift_phase,
			     const symmetry &S, int sn,
			     void *sum_) {
  long double *sum = (long double *) sum_;
  (void) shift; (void) shift_phase; (void) S; (void) sn; // unused
  if (fc->pol && fc->pol->energy[cgrid])
    LOOP_OVER_IVECS(fc->v, is, ie, idx)
      *sum += IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2)
	* fc->pol->energy[cgrid][idx];
}

double fields::thermo_energy_in_box(const geometric_volume &where) {
  long double sum = 0.0;
  FOR_ELECTRIC_COMPONENTS(c)
    if (!coordinate_mismatch(v.dim, c))
      integrate(thermo_integrand, (void *) &sum, where, c);
  return sum_to_all(sum);
}

struct flux_integrand_data {
  direction d; // flux direction
  component cE, cH; // components of E and H to get ExH in d
  long double sum;
};

static void flux_integrand(fields_chunk *fc, component cgrid,
			  ivec is, ivec ie,
			  vec s0, vec s1, vec e0, vec e1,
			  double dV0, double dV1,
			  ivec shift, complex<double> shift_phase,
			  const symmetry &S, int sn,
			  void *data_) {
  flux_integrand_data *data = (flux_integrand_data *) data_;

  (void) shift; // unused
  (void) shift_phase; // unused
  (void) cgrid; // == Dielectric

  component cE = S.transform(data->cE, -sn);
  component cH = S.transform(data->cH, -sn);

  // We're integrating Re[E * H*], and we assume that
  // S.phase_shift(cE,sn) == S.phase_shift(cH,sn), so the phases all cancel
  // ...except for at most an overall sign, which is fixed by checking
  // whether S.transform(data->d, -sn) is flipped, below.

  // offsets to average E and H components onto the Dielectric grid
  int oE1, oE2, oH1, oH2;
  fc->v.yee2diel_offsets(cE, oE1, oE2);
  fc->v.yee2diel_offsets(cH, oH1, oH2);

  long double sum = 0.0;

  for (int cmp = 0; cmp < 2; ++cmp)
    if (fc->f[cE][cmp] && fc->f[cH][cmp])
      LOOP_OVER_IVECS(fc->v, is, ie, idx) {
	double E, H;

	if (oE2)
	  E = 0.25 * (fc->f[cE][cmp][idx] + fc->f[cE][cmp][idx+oE1] +
		      fc->f[cE][cmp][idx+oE2] + fc->f[cE][cmp][idx+(oE1+oE2)]);
	else if (oE1)
	  E = 0.5 * (fc->f[cE][cmp][idx] + fc->f[cE][cmp][idx+oE1]);
	else
	  E = fc->f[cE][cmp][idx];

	if (oH2)
	  H = 0.25 * (fc->f[cH][cmp][idx] + fc->f[cH][cmp][idx+oH1] +
		      fc->f[cH][cmp][idx+oH2] + fc->f[cH][cmp][idx+(oH1+oH2)]);
	else if (oH1)
	  H = 0.5 * (fc->f[cH][cmp][idx] + fc->f[cH][cmp][idx+oH1]);
	else
	  H = fc->f[cH][cmp][idx];
	
	sum += IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2) * E * H;
    }

  data->sum += S.transform(data->d, -sn).flipped ? -sum : sum;
}

/* Compute ExH integral in box using current fields, ignoring fact
   that this E and H correspond to different times. */
double fields::flux_in_box_wrongH(direction d, const geometric_volume &where) {
  if (coordinate_mismatch(v.dim, d))
    return 0.0;

  component cE[2], cH[2];
  switch (d) {
  case X: cE[0] = Ey, cE[1] = Ez, cH[0] = Hz, cH[1] = Hy; break;
  case Y: cE[0] = Ez, cE[1] = Ex, cH[0] = Hx, cH[1] = Hz; break;
  case R: cE[0] = Ep, cE[1] = Ez, cH[0] = Hz, cH[1] = Hp; break;
  case P: cE[0] = Ez, cE[1] = Er, cH[0] = Hr, cH[1] = Hz; break;
  case Z:
    if (v.dim == Dcyl)
      cE[0] = Er, cE[1] = Ep, cH[0] = Hp, cH[1] = Hr;
    else
      cE[0] = Ex, cE[1] = Ey, cH[0] = Hy, cH[1] = Hx; 
    break;
  case NO_DIRECTION: abort("cannot get flux in NO_DIRECTION");
  }
  
  flux_integrand_data data;
  data.d = d;
  long double sum = 0.0;
  for (int i = 0; i < 2; ++i) {
    data.cE = cE[i]; data.cH = cH[i];
    data.sum = 0.0;
    integrate(flux_integrand, (void *) &data, where, Dielectric);
    sum += data.sum * (1 - 2*i);
  }
  return sum_to_all(sum) / (4*pi);
}

double fields::flux_in_box(direction d, const geometric_volume &where) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->backup_h();
  step_h();
  step_boundaries(H_stuff);
  step_h_source();
  double next_step_flux = flux_in_box_wrongH(d, where);
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->restore_h();
  return 0.5 * (next_step_flux + flux_in_box_wrongH(d, where));
}

flux_vol *fields::add_flux_vol(direction d, const geometric_volume &where) {
  if (where.dim != v.dim) abort("invalid dimensionality in add_flux_vol");
  if (d == NO_DIRECTION || coordinate_mismatch(v.dim, d))
    abort("invalid direction in add_flux_vol");
 return new flux_vol(this, d, where);
}

// As add_flux_vol, but infer direction from where (if possible)
flux_vol *fields::add_flux_plane(const geometric_volume &where) {
  return add_flux_vol(where.normal_direction(), where);
}

flux_vol *fields::add_flux_plane(const vec &p1, const vec &p2) {
  return add_flux_plane(geometric_volume(p1, p2));
}

} // namespace meep
