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

  // FIXME: do annoying special cases for r=0 in cylindrical coords
  
}

} // namespace meep
