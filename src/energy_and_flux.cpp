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

double fields::energy_in_box(const geometric_volume &otherv) {
  return thermo_energy_in_box(otherv) + field_energy_in_box(otherv);
}

double fields::field_energy_in_box(const geometric_volume &otherv) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->backup_h();
  step_h();
  step_boundaries(H_stuff);
  step_h_source();
  double next_step_magnetic_energy = magnetic_energy_in_box(otherv);
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->restore_h();

  return electric_energy_in_box(otherv) +
    0.5*next_step_magnetic_energy + 0.5*magnetic_energy_in_box(otherv);
}

struct dot_integrand_data {
  component c1, c2;
  long double sum;
};

static void dot_integrand(fields_chunk *fc, component cgrid,
			  ivec is, ivec ie,
			  vec s0, vec s1, vec e0, vec e1,
			  double dV0, double dV1,
			  vec shift, complex<double> shift_phase,
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
      LOOP_OVER_IVECS(fc->v, is, ie, idx) {
        double w1 = IVEC_LOOP_WEIGHT(1);
	double dV = dV0 + dV1 * loop_i2;
	double w12 = w1 * IVEC_LOOP_WEIGHT(2) * dV;
	double w123 = w12 * IVEC_LOOP_WEIGHT(3);
	data->sum += w123 * fc->f[c1][cmp][idx] * fc->f[c2][cmp][idx];;
	(void) loop_is1; (void) loop_is2; (void) loop_is3; // unused
    }
}

double fields::field_energy_in_box(component c,
				   const geometric_volume &otherv) {
  if (coordinate_mismatch(v.dim, component_direction(c)))
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
  integrate(dot_integrand, (void *) &data, otherv, c);
  return sum_to_all(data.sum) / (8*pi);
}

double fields::electric_energy_in_box(const geometric_volume &otherv) {
  long double sum = 0.0;
  FOR_ELECTRIC_COMPONENTS(c)
    sum += field_energy_in_box(c, otherv);
  return sum;
}

double fields::magnetic_energy_in_box(const geometric_volume &otherv) {
  long double sum = 0.0;
  FOR_MAGNETIC_COMPONENTS(c)
    sum += field_energy_in_box(c, otherv);
  return sum;
}

double fields::thermo_energy_in_box(const geometric_volume &otherv) {
  double energy = 0.0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      energy += chunks[i]->thermo_energy_in_box(otherv, S);
  return sum_to_all(energy);
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

double fields_chunk::thermo_energy_in_box(const geometric_volume &otherv,
                                          const symmetry &S) {
  // FIXME this is buggy when either parallel or using symmetry.
  if (pol) {
    return pol->total_energy(otherv);
  } else {
    return 0.0;
  }
}

} // namespace meep
