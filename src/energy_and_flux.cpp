/* Copyright (C) 2005-2007 Massachusetts Institute of Technology
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

#include "meep.hpp"
#include "meep_internals.hpp"

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
  synchronize_magnetic_fields();
  double cur_step_magnetic_energy = magnetic_energy_in_box(where);
  restore_magnetic_fields();
  return electric_energy_in_box(where) + cur_step_magnetic_energy;
}

static complex<double> dot_integrand(const complex<double> *fields,
				     const vec &loc, void *data_)
{
  (void) loc; (void) data_; // unused;
  return real(conj(fields[0]) * fields[1]);
}

double fields::field_energy_in_box(component c,
				   const geometric_volume &where) {
  if (coordinate_mismatch(v.dim, c))
    return 0.0;

  component cs[2];
  if (is_electric(c) || is_D(c)) {
    cs[0] = direction_component(Ex, component_direction(c));
    cs[1] = direction_component(Dx, component_direction(c));
  }
  else if (is_magnetic(c) || is_B(c)) {
    cs[0] = direction_component(Hx, component_direction(c));
    cs[1] = direction_component(Bx, component_direction(c));
  }
  else
    abort("invalid field component in field_energy_in_box");

  return real(integrate(2, cs, dot_integrand, 0, where)) * 0.5;
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

void fields_chunk::backup_component(component c) {
  // if we have D-P, then f_prev stores D-P, not D, complicating matters
  if (is_D(c) && have_d_minus_p) abort("backup-D with D-P not implemented");

  DOCMP {
    if (f[c][cmp]) {
      if (!f_backup[c][cmp])
	f_backup[c][cmp] = new double[v.ntot()];  
      memcpy(f_backup[c][cmp], 
	     f_prev[c][cmp] ? f_prev[c][cmp] : f[c][cmp],
	     v.ntot()*sizeof(double));    
    }
  }
}
  
void fields_chunk::restore_component(component c) {
  DOCMP if (f_backup[c][cmp]) {
    if (f[c][cmp])
      memcpy(f[c][cmp], f_prev[c][cmp] ? f_prev[c][cmp] : f_backup[c][cmp], 
	     v.ntot()*sizeof(double));
    if (f_prev[c][cmp])
      memcpy(f_prev[c][cmp], f_backup[c][cmp], v.ntot()*sizeof(double));
  }
}

void fields_chunk::average_with_backup(component c) {
  DOCMP {
    double *fc = f[c][cmp];
    double *backup = f_prev[c][cmp] ? f_prev[c][cmp] : f_backup[c][cmp];
    if (fc && backup)
      for (int i = 0; i < v.ntot(); i++)
	fc[i] = 0.5 * (fc[i] + backup[i]);
  }
}

void fields::synchronize_magnetic_fields() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine()) {
      FOR_B_COMPONENTS(c) chunks[i]->backup_component(c);
      FOR_MAGNETIC_COMPONENTS(c) chunks[i]->backup_component(c);
    }
  step_b();
  step_b_source();
  step_boundaries(B_stuff);
  update_h_from_b();
  step_boundaries(H_stuff);
  for (int i=0;i<num_chunks;i++) 
    if (chunks[i]->is_mine()) {
      FOR_B_COMPONENTS(c) chunks[i]->average_with_backup(c);
      FOR_MAGNETIC_COMPONENTS(c) chunks[i]->average_with_backup(c);
    }
}

void fields::restore_magnetic_fields() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine()) {
      FOR_B_COMPONENTS(c) chunks[i]->restore_component(c);
      FOR_MAGNETIC_COMPONENTS(c) chunks[i]->restore_component(c);
    }
}

static void thermo_chunkloop(fields_chunk *fc, int ichunk, component cgrid,
			     ivec is, ivec ie,
			     vec s0, vec s1, vec e0, vec e1,
			     double dV0, double dV1,
			     ivec shift, complex<double> shift_phase,
			     const symmetry &S, int sn,
			     void *sum_) {
  long double *sum = (long double *) sum_;
  (void)shift; (void)shift_phase; (void)S; (void)sn; (void)ichunk; // unused
  for (polarization *pol = fc->pol; pol; pol = pol->next)
    if (pol->energy[cgrid])
      LOOP_OVER_IVECS(fc->v, is, ie, idx)
	*sum += IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2)
	* pol->energy[cgrid][idx];
}

double fields::thermo_energy_in_box(const geometric_volume &where) {
  long double sum = 0.0;
  FOR_ELECTRIC_COMPONENTS(c)
    if (!coordinate_mismatch(v.dim, c))
      loop_in_chunks(thermo_chunkloop, (void *) &sum, where, c);
  return sum_to_all(sum);
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
  
  long double sum = 0.0;
  for (int i = 0; i < 2; ++i) {
    component cs[2];
    cs[0] = cE[i]; cs[1] = cH[i];
    sum += real(integrate(2, cs, dot_integrand, 0, where)) * (1 - 2*i);
  }
  return sum;
}

double fields::flux_in_box(direction d, const geometric_volume &where) {
  synchronize_magnetic_fields();
  double cur_step_flux = flux_in_box_wrongH(d, where);
  restore_magnetic_fields();
  return cur_step_flux;
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

/************************************************************************/

/* Note that computation of modal volume by this definition is
   somewhat problematic computationally, because we need to compute
   max|D*E|, which requires averaging discontinuous functions.  Hence,
   except for the special case of 2d TM polarization, the computed
   value tends to have a large error bar if the maximum lies on a
   dielectric boundary as it commonly does. 

   A better method would be to average only continuous quantities in
   order to compute the fields on the Dielectric grid, but this
   is more expensive and requires us to know the boundary orientation, and
   does not seem worth the trouble at this point. */

static complex<double> dot3_max_integrand(const complex<double> *fields,
				      const vec &loc, void *data_)
{
  (void) loc; (void) data_; // unused;
  return (real(conj(fields[0]) * fields[3]) +
	  real(conj(fields[1]) * fields[4]) +
	  real(conj(fields[2]) * fields[5]));
}

double fields::electric_energy_max_in_box(const geometric_volume &where) {
  component cs[6];
  if (v.dim == Dcyl) {
    cs[0] = Er; cs[1] = Ep; cs[2] = Ez;
    cs[3+0] = Dr; cs[3+1] = Dp; cs[3+2] = Dz;
  }
  else {
    cs[0] = Ex; cs[1] = Ey; cs[2] = Ez;
    cs[3+0] = Dx; cs[3+1] = Dy; cs[3+2] = Dz;
  }
  
  return max_abs(6, cs, dot3_max_integrand, 0, where) * 0.5;
}

/* "modal" volume according to definition in:
      E. M. Purcell, Phys. Rev. B 69, 681 (1946).
    (based on spontaneous emission enhancement). */
double fields::modal_volume_in_box(const geometric_volume &where) {
  return electric_energy_in_box(where) / electric_energy_max_in_box(where);
}

/************************************************************************/

  /* compute integral f(x) * Re[conj(f1)*f2] * 0.5, which is useful for
     perturbation theory, etcetera, where f1 and f2 are two field components
     on the same Yee lattice (e.g. Hx and Hx or Ex and Dx). */

typedef double (*fx_func)(const vec &);

static complex<double> dot_fx_integrand(const complex<double> *fields,
					const vec &loc, void *data_) {
  fx_func fx = (fx_func) data_;
  return (real(conj(fields[0]) * fields[1]) * fx(loc));
}

/* computes integral of f(x) * |E|^2 / integral epsilon*|E|^2 */
double fields::electric_sqr_weighted_integral(double (*f)(const vec &),
					     const geometric_volume &where) {
  double sum = 0.0;
  FOR_ELECTRIC_COMPONENTS(c) 
    if (!coordinate_mismatch(v.dim, component_direction(c))) {
      component cs[2];
      cs[0] = cs[1] = direction_component(Ex, component_direction(c));
      sum += real(integrate(2, cs, dot_fx_integrand, (void *) f, where));
    }
  return sum * 0.5 / electric_energy_in_box(where);
}

/* computes integral of f(x) * epsilon*|E|^2 / integral epsilon*|E|^2 */
double fields::electric_energy_weighted_integral(double (*f)(const vec &),
					     const geometric_volume &where) {
  double sum = 0.0;
  FOR_ELECTRIC_COMPONENTS(c) 
    if (!coordinate_mismatch(v.dim, component_direction(c))) {
      component cs[2];
      cs[0] = direction_component(Ex, component_direction(c));
      cs[1] = direction_component(Dx, component_direction(c));
      sum += real(integrate(2, cs, dot_fx_integrand, (void *) f, where));
    }
  return sum * 0.5 / electric_energy_in_box(where);
}

} // namespace meep

