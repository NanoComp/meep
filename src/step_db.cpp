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

void fields::step_db(field_type ft) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_db(ft);
}

void fields_chunk::step_db(field_type ft) {
  if (ft != B_stuff && ft != D_stuff)
    abort("bug - step_db should only be called for B or D");

  bool have_pml = false;
  FOR_FT_COMPONENTS(ft, cc)
    if (s->sigsize[cycle_direction(v.dim,component_direction(cc),1)] > 1)
	have_pml = true;
  if (have_pml) FOR_FT_COMPONENTS(ft, cc) DOCMP
    if (f[cc][cmp]) {
      if (!f_prev[cc][cmp])
        f_prev[cc][cmp] = new double[v.ntot()];
      // update_e_from_d requires previous D-P, not D, if D-P is present
      memcpy(f_prev[cc][cmp],
	     f_minus_p[cc][cmp] ? f_minus_p[cc][cmp] : f[cc][cmp],
	     v.ntot()*sizeof(double));
    }

  if (v.dim != Dcyl) {
    DOCMP FOR_FT_COMPONENTS(ft, cc)
      if (f[cc][cmp]) {
	const component c_p=plus_component[cc], c_m=minus_component[cc];
	const direction d_deriv_p = plus_deriv_direction[cc];
	const direction d_deriv_m = minus_deriv_direction[cc];
	const direction d_c = component_direction(cc);
	const bool have_p = have_plus_deriv[cc];
	const bool have_m = have_minus_deriv[cc];
	const direction dsig0 = cycle_direction(v.dim,d_c,1);
	const bool have_pml = s->sigsize[dsig0] > 1;
	const direction dsig = have_pml ? dsig0 : NO_DIRECTION;
	int stride_p = have_p?v.stride(d_deriv_p):0;
	int stride_m = have_m?v.stride(d_deriv_m):0;
	RESTRICT const double *f_p = have_p?f[c_p][cmp]:NULL;
	RESTRICT const double *f_m = have_m?f[c_m][cmp]:NULL;
	RESTRICT double *the_f = f[cc][cmp];

	if (ft == D_stuff) { // strides are opposite sign for H curl
	  stride_p = -stride_p;
	  stride_m = -stride_m;
	}

	step_curl(the_f, cc, f_p, f_m, stride_p, stride_m, v, Courant, 
		  dsig, s->sig[dsig], s->siginv[dsig],
		  dt, s->conductivity[cc][d_c], s->condinv[cc][d_c]);
      }
  } else if (v.dim == Dcyl) {
  } else {
    abort("Unsupported dimension.\n");
  }
}

} // namespace meep
