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

#include "dactyl.h"
#include "dactyl_internals.h"

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
    vol += v.dV((component)c,i).intersection(v.dV((component)c,i));
  return vol;
}

double fields::total_energy() {
  return energy_in_box(user_volume);
}

double fields::field_energy() {
  return field_energy_in_box(user_volume);
}

double fields::energy_in_box(const volume &otherv) {
  return thermo_energy_in_box(otherv) + field_energy_in_box(otherv);
}

double fields::field_energy_in_box(const volume &otherv) {
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

double fields::electric_energy_in_box(const volume &otherv) {
  double energy = 0.0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      energy += chunks[i]->electric_energy_in_box(otherv, S);
  return sum_to_all(energy);
}

double fields::magnetic_energy_in_box(const volume &otherv) {
  double energy = 0.0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      energy += chunks[i]->magnetic_energy_in_box(otherv, S);
  return sum_to_all(energy);
}

double fields::thermo_energy_in_box(const volume &otherv) {
  double energy = 0.0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      energy += chunks[i]->thermo_energy_in_box(otherv, S);
  return sum_to_all(energy);
}

double fields_chunk::backup_h() {
  DOCMP {
    for (int c=0;c<10;c++)
      if (f[c][cmp] && is_magnetic((component)c)) {
        if (f_backup[c][cmp] == NULL)
          f_backup[c][cmp] = new double[v.ntot()];
        if (f_backup_pml[c][cmp] == NULL)
          f_backup_pml[c][cmp] = new double[v.ntot()];
      }
  }
  DOCMP {
    for (int c=0;c<10;c++)
      if (f[c][cmp] && is_magnetic((component)c)) {
        for (int i=0;i<v.ntot();i++) f_backup[c][cmp][i] = f[c][cmp][i];
        for (int i=0;i<v.ntot();i++) f_backup_pml[c][cmp][i] = f_pml[c][cmp][i];
      }
  }
}

double fields_chunk::restore_h() {
  DOCMP {
    for (int c=0;c<10;c++)
      if (f[c][cmp] && is_magnetic((component)c)) {
        for (int i=0;i<v.ntot();i++) f[c][cmp][i] = f_backup[c][cmp][i];
        for (int i=0;i<v.ntot();i++) f_pml[c][cmp][i] = f_backup_pml[c][cmp][i];
      }
  }
}

double fields_chunk::electric_energy_in_box(const volume &otherv,
                                            const symmetry &S) {
  double energy = 0;
  DOCMP
    FOR_ELECTRIC_COMPONENTS(c)
      if (f[c][cmp])
        for (int i=0;i<v.ntot();i++) {
          const ivec p0 = v.iloc((component)c,i);
          if (S.is_primitive(p0) && v.owns(p0)) {
            int num_times_mapped_to_self = 0;
            for (int sn=0;sn<S.multiplicity();sn++)
              if (S.transform(p0,sn)==p0) num_times_mapped_to_self++;
            for (int sn=0;sn<S.multiplicity();sn++) {
              const ivec pn = S.transform(p0,sn);
              if (otherv.owns(pn))
                // FIXME I need to rewrite this to deal with anisotropic
                // dielectric stuff.
                energy += otherv.intersection(v.dV(pn))*
                  f[c][cmp][i]*
                  (1./ma->inveps[c][component_direction(c)][i]*f[c][cmp][i])
                  /num_times_mapped_to_self;
            }
          }
        }
  return energy*(1.0/(8*pi));
}

double fields_chunk::magnetic_energy_in_box(const volume &otherv,
                                            const symmetry &S) {
  double energy = 0;
  DOCMP
    for (int c=0;c<10;c++)
      if (f[c][cmp] && is_magnetic((component)c))
        for (int i=0;i<v.ntot();i++) {
          const ivec p0 = v.iloc((component)c,i);
          if (S.is_primitive(p0) && v.owns(p0)) {
            int num_times_mapped_to_self = 0;
            for (int sn=0;sn<S.multiplicity();sn++)
              if (S.transform(p0,sn)==p0) num_times_mapped_to_self++;
            for (int sn=0;sn<S.multiplicity();sn++) {
              const ivec pn = S.transform(p0,sn);
              if ((pn!=p0 || sn==0) && otherv.owns(pn))
                energy += otherv.intersection(v.dV(pn))*
                  f[c][cmp][i]*f[c][cmp][i];
            }
          }
        }
  return energy*(1.0/(8*pi));
}

double fields_chunk::thermo_energy_in_box(const volume &otherv,
                                          const symmetry &S) {
  // FIXME this is buggy when either parallel or using symmetry.
  if (pol) {
    return pol->total_energy(otherv)/(4*pi);
  } else {
    return 0.0;
  }
}
