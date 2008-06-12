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

#include "meep.hpp"
#include "meep_internals.hpp"

#define RESTRICT

namespace meep {

static inline double it(int cmp, double *(f[2]), int ind) {
  return (f[1-cmp]) ? (1-2*cmp)*f[1-cmp][ind] : 0;
}

inline int rstart_0(const volume &v, double m) {
  return (int) max(0.0, m - (int)(v.origin_r()*v.a+0.5) - 1.0);
}
inline int rstart_1(const volume &v, double m) {
  return (int) max(1.0, m - (int)(v.origin_r()*v.a+0.5));
}

void fields::step_d() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_d();
}

void fields_chunk::step_d() {
  const volume v = this->v;
  bool have_pml = false;
  FOR_D_COMPONENTS(cc)
    if (s->sigsize[(component_direction(cc)+1)%3] > 1) have_pml = true;
  if (have_pml) backup_d();

  if (v.dim != Dcyl) {
    DOCMP FOR_D_COMPONENTS(cc)
      if (f[cc][cmp]) {
	const component c_p=plus_component[cc], c_m=minus_component[cc];
	const direction d_deriv_p = plus_deriv_direction[cc];
	const direction d_deriv_m = minus_deriv_direction[cc];
	const direction d_c = component_direction(cc);
	const bool have_p = have_plus_deriv[cc];
	const bool have_m = have_minus_deriv[cc];
	const bool have_pml = s->sigsize[(d_c+1)%3] > 1;
	const direction dsig = have_pml?(direction)((d_c+1)%3):NO_DIRECTION;
	const int stride_p = have_p?v.stride(d_deriv_p):0;
	const int stride_m = have_m?v.stride(d_deriv_m):0;
	RESTRICT const double *f_p = have_p?f[c_p][cmp]:NULL;
	RESTRICT const double *f_m = have_m?f[c_m][cmp]:NULL;
	RESTRICT double *the_f = f[cc][cmp];
	if (v.dim == D3) {
	  step_curl(the_f, cc, f_p, f_m, -stride_p, -stride_m, v, Courant, 
	dsig, have_pml?s->sig[dsig]:NULL, have_pml?s->siginv[dsig]:NULL);
	} else if (v.dim == D2) {
	  if (f[Ez][cmp]) { // TM
	    if (have_p && have_m)
	     step_curl(the_f, cc, f_p, f_m, -stride_p, -stride_m, v, Courant,
	dsig, have_pml?s->sig[dsig]:NULL, have_pml?s->siginv[dsig]:NULL);
	  } 
	  if (f[Hz][cmp]) { // TE
	    if (!have_p && have_m)
	      step_curl(the_f, cc, f_m, NULL, -stride_m, 0, v, -Courant, 
			dsig, have_pml?s->sig[dsig]:NULL, have_pml?s->siginv[dsig]:NULL);
	    else if (have_p && !have_m)
	      step_curl(the_f, cc, f_p, NULL, -stride_p, 0, v, Courant, 
			dsig, have_pml?s->sig[dsig]:NULL, have_pml?s->siginv[dsig]:NULL);
	  }
	} else if (v.dim == D1) {
	  if (!have_p && have_m)
	    step_curl(the_f, cc, f_m, NULL, -stride_m, 0, v, -Courant, 
        dsig, have_pml?s->sig[dsig]:NULL, have_pml?s->siginv[dsig]:NULL);
	}
      }
  } else if (v.dim == Dcyl) {
  } else {
    abort("Unsupported dimension.\n");
  }
}

} // namespace meep
