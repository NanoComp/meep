/* Copyright (C) 2005-2021 Massachusetts Institute of Technology
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

#include "material_data.hpp"

#include "meep/mympi.hpp"

namespace meep_geom {

medium_struct::medium_struct(double epsilon) :
      epsilon_diag{epsilon, epsilon, epsilon},
      epsilon_offdiag{},
      mu_diag{1, 1, 1},
      mu_offdiag{},
      E_susceptibilities(), H_susceptibilities(),
      E_chi2_diag{},
      E_chi3_diag{},
      H_chi2_diag{},
      H_chi3_diag{},
      D_conductivity_diag{},
      B_conductivity_diag{}
    {}

void medium_struct::check_offdiag_im_zero_or_abort() const {
  if (epsilon_offdiag.x.im != 0 || epsilon_offdiag.y.im != 0 ||
      epsilon_offdiag.z.im != 0 || mu_offdiag.x.im != 0 || mu_offdiag.y.im != 0 ||
      mu_offdiag.z.im != 0) {
    meep::abort("Found non-zero imaginary part of epsilon or mu offdiag.\n");
  }
}

}  // namespace meep_geom
