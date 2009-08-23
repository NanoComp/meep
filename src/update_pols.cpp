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

  realnum *w[NUM_FIELD_COMPONENTS][2];
  FOR_COMPONENTS(c) DOCMP2 w[c][cmp] = f_w[c][cmp] ? f_w[c][cmp] : f[c][cmp];

  for (polarization_state *p = pol[ft]; p; p = p->next) {

    // Lazily allocate polarizations P where needed:
    bool allocated_pol = false;
    FOR_FT_COMPONENTS(ft, c)
      if (!p->P[c][0] && p->s->needs_P(c, f)) {
	DOCMP {
	  p->P[c][cmp] = new realnum[gv.ntot()];
	  memset(p->P[c][cmp], 0, gv.ntot() * sizeof(realnum));
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
    p->s->update_P(p->P, w, f_w_prev, dt, gv, p->data);
  }

  return allocated_fields;
}

} // namespace meep
