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

double fields::electric_energy_in_box(const geometric_volume &otherv) {
  double energy = 0.0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      energy += chunks[i]->electric_energy_in_box(otherv, S);
  return sum_to_all(energy);
}

double fields::magnetic_energy_in_box(const geometric_volume &otherv) {
  double energy = 0.0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      energy += chunks[i]->magnetic_energy_in_box(otherv, S);
  return sum_to_all(energy);
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
        for (int i=0;i<v.ntot();i++) f_backup[c][cmp][i] = f[c][cmp][i];
        for (int i=0;i<v.ntot();i++) f_backup_p_pml[c][cmp][i] = f_p_pml[c][cmp][i];
        for (int i=0;i<v.ntot();i++) f_backup_m_pml[c][cmp][i] = f_m_pml[c][cmp][i];
      }
}

void fields_chunk::restore_h() {
  DOCMP FOR_MAGNETIC_COMPONENTS(c)
      if (f[c][cmp]) {
        for (int i=0;i<v.ntot();i++) f[c][cmp][i] = f_backup[c][cmp][i];
        for (int i=0;i<v.ntot();i++) f_p_pml[c][cmp][i] = f_backup_p_pml[c][cmp][i];
        for (int i=0;i<v.ntot();i++) f_m_pml[c][cmp][i] = f_backup_m_pml[c][cmp][i];
      }
}

double fields_chunk::electric_energy_in_box(const geometric_volume &otherv,
                                            const symmetry &S) {
  double energy = 0;
  DOCMP
    FOR_E_AND_D(c,dc)
      if (f[c][cmp])
        for (int i=0;i<v.ntot();i++) {
          const ivec p0 = v.iloc(c,i);
          for (int sn=0;sn<S.multiplicity();sn++)
            energy += (otherv & S.transform(gv & v.dV(p0),sn)).full_volume()*
              f[c][cmp][i]*f[dc][cmp][i];
        }
  return energy*(1.0/(8*pi));
}

double fields_chunk::magnetic_energy_in_box(const geometric_volume &otherv,
                                            const symmetry &S) {
  double energy = 0;
  DOCMP
    FOR_MAGNETIC_COMPONENTS(c)
      if (f[c][cmp])
        for (int i=0;i<v.ntot();i++) {
          const ivec p0 = v.iloc(c,i);
          for (int sn=0;sn<S.multiplicity();sn++)
            energy += (otherv & S.transform(gv & v.dV(p0),sn)).full_volume()*
              f[c][cmp][i]*f[c][cmp][i];
        }
  return energy*(1.0/(8*pi));
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
