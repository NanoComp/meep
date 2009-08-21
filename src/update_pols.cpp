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
#include <string.h>

#include "meep.hpp"
#include "meep_internals.hpp"
#include "config.h"

namespace meep {

void fields::update_pols(field_type ft) {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      if (chunks[i]->update_pols(ft))
	chunk_connections_valid = false;

  /* synchronize to avoid deadlocks if one process decides it needs
     to allocate E or H ... */
  chunk_connections_valid = and_to_all(chunk_connections_valid);
}

bool fields_chunk::update_pols(field_type ft) {
  bool allocated_fields = false;

  for (poldata *p = pol[ft]; p; p = p->next) {

  DOCMP FOR_FT_COMPONENTS(ft, c) if (f[c][cmp])
    for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
      if (np->pb->ft != ft) abort("bug in update_pols");
      const double cn = 2 - op->pb->omeganot*op->pb->omeganot;
      const double co = 0.5 * op->pb->gamma - 1;
      const double funinv = 1.0 / (1 + 0.5*op->pb->gamma);
      const realnum *fE = f[c][cmp];
      const realnum *npP = np->P[c][cmp], *nps = np->s[c];
      realnum *opP = op->P[c][cmp], *npenergy = np->energy[c];
      if (npenergy)
	for (int i = 0; i < ntot; ++i) {
	  npenergy[i] += 0.5 * (npP[i] - opP[i]) * fE[i];
	  opP[i] = funinv * (cn * npP[i] + co * opP[i] + nps[i] * fE[i]);
	}
	allocated_pol = true;
      }
    if (allocated_pol) {
      allocated_fields = true;
      if (p->data) { // TODO: warning or error message in this weird case?
	delete[] p->data; p->data = NULL;
      }
    }

    // Lazily allocate internal polarization data:
    if (!p->data) {
      p->ndata = p->s->num_internal_data(p->P, gv);
      if (p->ndata) {
	p->data = new realnum[p->ndata];
	p->s->init_internal_data(p->P, gv, p->data);
	allocated_fields = true;
      }
    }

    // Finally, timestep the polarizations:
    p->s->update_P(p->P, f, f_prev, dt, gv, p->data);
  }

  return allocated_fields;
}

} // namespace meep
