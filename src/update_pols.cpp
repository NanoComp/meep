/* Copyright (C) 2005-2009 Massachusetts Institute of Technology
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
#include "config.h"

namespace meep {

void fields::update_pols(field_type ft) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->update_pols(ft);
}

void fields_chunk::update_pols(field_type ft) {
  const int ntot = s->gv.ntot();
  polarization *pol = pols[ft];
  polarization *olpol = olpols[ft];

  DOCMP FOR_FT_COMPONENTS(ft, c) if (f[c][cmp])
    for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
      if (np->pb->ft != ft) abort("bug in update_pols");
      const double cn = 2 - op->pb->omeganot*op->pb->omeganot;
      const double co = 0.5 * op->pb->gamma - 1;
      const double funinv = 1.0 / (1 + 0.5*op->pb->gamma);
      const realnum * restrict fE = f_w[c][cmp] ? f_w[c][cmp] : f[c][cmp];
      const realnum * restrict npP = np->P[c][cmp], * restrict nps = np->s[c];
      realnum * restrict opP = op->P[c][cmp], * restrict npenergy = np->energy[c];
      if (npenergy)
	for (int i = 0; i < ntot; ++i) {
	  npenergy[i] += 0.5 * (npP[i] - opP[i]) * fE[i];
	  opP[i] = funinv * (cn * npP[i] + co * opP[i] + nps[i] * fE[i]);
	}
      else
	for (int i = 0; i < ntot; ++i)
	  opP[i] = funinv * (cn * npP[i] + co * opP[i] + nps[i] * fE[i]);
    }

  /* the old polarization is now the new polarization */
  olpols[ft] = pol;
  pols[ft] = olpol;
}

} // namespace meep
