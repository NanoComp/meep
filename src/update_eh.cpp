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

#include <string.h>

#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {
  
void fields::update_eh(field_type ft) {
  if (ft != E_stuff && ft != H_stuff) abort("update_eh only works with E/H");
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff; // for sources etc.
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine()) {
      src_vol *save_sources = chunks[i]->sources[ft2];
      if (disable_sources) chunks[i]->sources[ft2] = NULL; // temporary
      if (chunks[i]->update_eh(ft))
	chunk_connections_valid = false; // E/H allocated - reconnect chunks 
      chunks[i]->sources[ft2] = save_sources;
    }

  /* synchronize to avoid deadlocks if one process decides it needs
     to allocate E or H ... */
  chunk_connections_valid = and_to_all(chunk_connections_valid);
}

bool fields_chunk::update_eh(field_type ft) {
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff; // for sources etc.
  bool allocated_eh = false;

  bool have_int_sources = false;
  for (src_vol *sv = sources[ft2]; sv; sv = sv->next)
    if (sv->t->is_integrated) {
      have_int_sources = true;
      break;
    }

  FOR_FT_COMPONENTS(ft2, dc) DOCMP {
    if (f[dc][cmp] && (pols[ft] || have_int_sources)) {
      if (!f_minus_p[dc][cmp]) f_minus_p[dc][cmp] = new realnum[v.ntot()];
    }
    else if (f_minus_p[dc][cmp]) { // remove unneeded f_minus_p
      delete[] f_minus_p[dc][cmp];
      f_minus_p[dc][cmp] = 0;
    }
  }
  bool have_f_minus_p = false;
  FOR_FT_COMPONENTS(ft2, dc) if (f_minus_p[dc][0]) {
    have_f_minus_p = true;
    break;
  }

  const int ntot = s->v.ntot();

  //////////////////////////////////////////////////////////////////////////
  // First, initialize f_minus_p to D - P, if necessary

  if (have_f_minus_p) {
    if (pols[ft]) {
      FOR_FT_COMPONENTS(ft, ec) if (f[ec][0]) {
	for (polarization *np=pols[ft],*op=olpols[ft]; np; 
	     np=np->next, op=op->next) {
	  if (np->energy[ec] && op->energy[ec]) {
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
	}
	component dc = direction_component(first_field_component(ft2),
					   component_direction(ec));
	DOCMP {
	  for (int i=0;i<ntot;i++) {
	    double sum = f[dc][cmp][i];
            for (polarization *p = pols[ft]; p; p = p->next)
              sum -= p->P[ec][cmp][i];
            f_minus_p[dc][cmp][i] = sum;
	  }
	}
      }
    }
    else {
      FOR_FT_COMPONENTS(ft2, dc) if (f[dc][0]) DOCMP
	memcpy(f_minus_p[dc][cmp], f[dc][cmp], ntot * sizeof(realnum));
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Next, subtract time-integrated sources (i.e. polarizations, not currents)

  if (have_f_minus_p) {
    for (src_vol *sv = sources[ft2]; sv; sv = sv->next) {  
      if (sv->t->is_integrated && f[sv->c][0] && ft == type(sv->c)) {
	component c = field_type_component(ft2, sv->c);
	for (int j = 0; j < sv->npts; ++j) { 
	  const complex<double> A = sv->dipole(j);
	  DOCMP {
	    f_minus_p[c][cmp][sv->index[j]] -= 
	      (cmp) ? imag(A) :  real(A);
	  }
	}
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Finally, compute E = chi1inv * D
  
  realnum *dmp[NUM_FIELD_COMPONENTS][2];
  if (have_f_minus_p) {
    FOR_FT_COMPONENTS(ft2,dc) DOCMP2 dmp[dc][cmp] = f_minus_p[dc][cmp];
  } else {
    FOR_FT_COMPONENTS(ft2,dc) DOCMP2 dmp[dc][cmp] = f[dc][cmp];
  }

  DOCMP FOR_FT_COMPONENTS(ft,ec) if (f[ec][cmp]) {
    if (type(ec) != ft) abort("bug in FOR_FT_COMPONENTS");
    component dc = field_type_component(ft2, ec);
    const direction d_ec = component_direction(ec);
    const int s_ec = stride_any_direction[d_ec];
    const direction d_1 = cycle_direction(v.dim, d_ec, 1);
    const component dc_1 = direction_component(dc,d_1);
    const int s_1 = stride_any_direction[d_1];
    const direction d_2 = cycle_direction(v.dim, d_ec, 2);
    const component dc_2 = direction_component(dc,d_2);
    const int s_2 = stride_any_direction[d_2];

    direction dsig = d_2;
    direction dsigg = d_ec;
    direction dsig1 = d_1;
    direction dsig1inv = d_ec;
    direction dsig2 = d_2;
    direction dsig2inv = d_1;

    // lazily allocate any E/H fields that are needed (H==B initially)
    if (f[ec][cmp] == f[dc][cmp]
	&& (s->chi1inv[ec][d_ec] || have_f_minus_p
	    || s->sigsize[dsig] > 1
	    || s->sigsize[dsigg] > 1
	    || (s->sigsize[dsig1] > 1
		&& (s->chi1inv[ec][d_1] || s->chi1inv[ec][d_2])))) {
      f[ec][cmp] = new realnum[v.ntot()];
      memcpy(f[ec][cmp], f[dc][cmp], v.ntot() * sizeof(realnum));
      allocated_eh = true;
    }

    if (f[ec][cmp] != f[dc][cmp])
      step_update_EDHB(f[ec][cmp], ec, v, 
		       dmp[dc][cmp], dmp[dc_1][cmp], dmp[dc_2][cmp],
		       f_prev[dc][cmp], f_prev[dc_1][cmp], f_prev[dc_2][cmp],
		       s->chi1inv[ec][d_ec], dmp[dc_1][cmp]?s->chi1inv[ec][d_1]:NULL, dmp[dc_2][cmp]?s->chi1inv[ec][d_2]:NULL,
		       s_ec, s_1, s_2, s->chi2[ec], s->chi3[ec],
		       dsig, s->sig[dsig], s->siginv[dsig],
		       dsigg, s->sig[dsigg],
		       dsig1, s->sig[dsig1],
		       dsig1inv, s->sig[dsig1inv],
		       dsig2, s->sig[dsig2],
		       dsig2inv, s->sig[dsig2inv],
		       s->sigsize[dsig],s->sigsize[dsigg],s->sigsize[dsig1]);
  }

  /* Do annoying special cases for r=0 in cylindrical coords.  Note
     that this only really matters for field output; the Ez and Ep
     components at r=0 don't usually affect the fields elsewhere
     because of the form of Maxwell's equations in cylindrical coords. */
  // (FIXME: handle Kerr case?).
  if (v.dim == Dcyl && v.origin_r() == 0.0)
    DOCMP FOR_FT_COMPONENTS(ft,ec) if (f[ec][cmp] && (ec == Ep || ec == Ez
						      || ec == Hr)) {
      component dc = field_type_component(ft2, ec);
      if (f[ec][cmp] == f[dc][cmp]) continue;
      const int yee_idx = v.yee_index(ec);
      const int d_ec = component_direction(ec);
      const int sR = stride_any_direction[R];
      realnum *E = f[ec][cmp];
      const realnum *D = have_f_minus_p ? f_minus_p[dc][cmp] : f[dc][cmp];
      const realnum *chi1inv = s->chi1inv[ec][d_ec];
      if (chi1inv)
	for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
	  const int i = yee_idx + iZ - sR;
	  E[i] = chi1inv[i] * D[i];
	}
      else
	for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
	  const int i = yee_idx + iZ - sR;
	  E[i] = D[i];
	}
    }
  
  return allocated_eh;
}

} // namespace meep
