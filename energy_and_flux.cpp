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
    vol += v.dv((component)c,i);
  return vol;
}

double fields::total_energy() {
  return energy_in_box(v);
}

double fields::field_energy() {
  return field_energy_in_box(v);
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
      energy += chunks[i]->electric_energy_in_box(otherv);
  return sum_to_all(energy);
}

double fields::magnetic_energy_in_box(const volume &otherv) {
  double energy = 0.0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      energy += chunks[i]->magnetic_energy_in_box(otherv);
  return sum_to_all(energy);
}

double fields::thermo_energy_in_box(const volume &otherv) {
  double energy = 0.0;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      energy += chunks[i]->thermo_energy_in_box(otherv);
  return sum_to_all(energy);
}

double fields_chunk::backup_h() {
  DOCMP {
    for (int c=0;c<10;c++)
      if (v.has_field((component)c) && is_magnetic((component)c)) {
        if (f_backup[c][cmp] == NULL)
          f_backup[c][cmp] = new double[v.ntot()];
        if (f_backup_pml[c][cmp] == NULL)
          f_backup_pml[c][cmp] = new double[v.ntot()];
      }
  }
  DOCMP {
    for (int c=0;c<10;c++)
      if (v.has_field((component)c) && is_magnetic((component)c)) {
        for (int i=0;i<v.ntot();i++) f_backup[c][cmp][i] = f[c][cmp][i];
        for (int i=0;i<v.ntot();i++) f_backup_pml[c][cmp][i] = f_pml[c][cmp][i];
      }
  }
}

double fields_chunk::restore_h() {
  DOCMP {
    for (int c=0;c<10;c++)
      if (v.has_field((component)c) && is_magnetic((component)c)) {
        for (int i=0;i<v.ntot();i++) f[c][cmp][i] = f_backup[c][cmp][i];
        for (int i=0;i<v.ntot();i++) f_pml[c][cmp][i] = f_backup_pml[c][cmp][i];
      }
  }
}

double fields_chunk::electric_energy_in_box(const volume &otherv) {
  double energy = 0;
  DOCMP {
    for (int c=0;c<10;c++)
      if (v.has_field((component)c) && is_electric((component)c))
        for (int i=0;i<v.ntot();i++)
          if (otherv.owns(v.loc((component)c,i)))
            energy += otherv.intersection(v.dV((component)c,i))*
              f[c][cmp][i]*(1./ma->inveps[c][i]*f[c][cmp][i]);
  }
  return energy/(8*pi);
}

double fields_chunk::magnetic_energy_in_box(const volume &otherv) {
  double energy = 0;
  DOCMP {
    for (int c=0;c<10;c++)
      if (v.has_field((component)c) && is_magnetic((component)c))
        for (int i=0;i<v.ntot();i++) {
          vec here = v.loc((component)c,i);
          if (v.owns(here) && otherv.owns(here))
            energy += otherv.intersection(v.dV((component)c,i))*
              f[c][cmp][i]*f[c][cmp][i];
        }
  }
  return energy/(8*pi);
}

double fields_chunk::thermo_energy_in_box(const volume &otherv) {
  if (pol) {
    return pol->total_energy(otherv)/(4*pi);
  } else {
    return 0.0;
  }
}

complex<double> fields::flux(const vec &end1, const vec &end2) {
  switch (v.dim) {
  case d1: return oned_flux(end1); // assume end1 and end2 are identical
  }
  abort("Don't yet support flux planes in this dimensionality!\n");
}

complex<double> fields::oned_flux(const vec &loc) {
  monitor_point p;
  get_point(&p, loc);
  return p.get_component(Ex)*p.get_component(Hy)*(1.0/pi/8/c);
}
