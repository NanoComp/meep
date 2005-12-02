/* Copyright (C) 2005 Massachusetts Institute of Technology
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

#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {

  inline double calc_nonlinear_inveps(const double Dsqr, const double e,
                                      const double alpha) {
    return 1.0 - (alpha*Dsqr)/(e + 3*(alpha*Dsqr));
  }

  void fields_chunk::update_e_from_d_update(void) {
  const int ntot = s->v.ntot();
#include "update_e_from_d_update.hpp"

  /* Do annoying special cases for r=0 in cylindrical coords.  Note
     that this only really matters for field output; the Ez and Ep
     components at r=0 don't usually affect the fields elsewhere
     because of the form of Maxwell's equations in cylindrical coords. */
  // (FIXME: handle Kerr case).
  if (v.dim == Dcyl && v.origin_r() == 0.0)
    DOCMP FOR_E_AND_D(ec,dc) if (f[ec][cmp] && ec != Er) {
      const int yee_idx = v.yee_index(ec);
      const int d_ec = component_direction(ec);
      const int sR = stride_any_direction[R];
      const double *D = have_d_minus_p ? d_minus_p[ec][cmp] : f[dc][cmp];
      for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
	const int i = yee_idx + iZ - sR;
	f[ec][cmp][i] = s->inveps[ec][d_ec][i] * D[i];
      }
    }
}

} // namespace meep
