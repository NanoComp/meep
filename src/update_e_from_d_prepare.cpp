/* Copyright (C) 2006 Massachusetts Institute of Technology
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

#include "meep_internals.hpp"

namespace meep {

void fields_chunk::update_e_from_d_prepare(void) {
  const int ntot = s->v.ntot();

  if (have_d_minus_p) {
    if (pol) {
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
	for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
	  if (is_real) for (int i = 0; i < ntot; ++i) {
	    np->energy[ec][i] = op->energy[ec][i] +
	      (0.5)*(np->P[ec][0][i] - op->P[ec][0][i])
              * f[ec][0][i];
	  }
	  else for (int i = 0; i < ntot; ++i) {
            np->energy[ec][i] = op->energy[ec][i] +
              (0.5)*(np->P[ec][0][i] - op->P[ec][0][i])
              * f[ec][0][i] +
              (0.5)*(np->P[ec][1][i] - op->P[ec][1][i])
              * f[ec][1][i];
          }
	}
	DOCMP {
	  for (int i=0;i<ntot;i++) {
	    double sum = f[dc][cmp][i];
            for (polarization *p = pol; p; p = p->next) {
              sum -= p->P[ec][cmp][i];
            }	  
            d_minus_p[ec][cmp][i] = sum;
	  }
	}
      }
    }
    else {
      FOR_E_AND_D(ec,dc) if (f[ec][0]) DOCMP
	memcpy(d_minus_p[ec][cmp], f[dc][cmp], ntot * sizeof(double));
    }
  }
}

} // namespace meep
