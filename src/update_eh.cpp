/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
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
#include <assert.h>

#include "meep.hpp"
#include "meep_internals.hpp"

using namespace std;

namespace meep {

void fields::update_eh(field_type ft, bool skip_w_components) {
  if (ft != E_stuff && ft != H_stuff) meep::abort("update_eh only works with E/H");

  // split the chunks' volume into subdomains for tiled execution of update_eh loop
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine() && changed_materials) {
      bool is_aniso = false;
      FOR_FT_COMPONENTS(ft, cc) {
        const direction d_c = component_direction(cc);
        const direction d_1 = cycle_direction(chunks[i]->gv.dim, d_c, 1);
        const direction d_2 = cycle_direction(chunks[i]->gv.dim, d_c, 2);
        if (chunks[i]->s->chi1inv[cc][d_1] && chunks[i]->s->chi1inv[cc][d_2]) {
          is_aniso = true;
          break;
        }
      }
      if (!chunks[i]->gvs_eh[ft].empty()) chunks[i]->gvs_eh[ft].clear();
      if (loop_tile_base_eh > 0 && is_aniso) {
        split_into_tiles(chunks[i]->gv, &chunks[i]->gvs_eh[ft], loop_tile_base_eh);
        check_tiles(chunks[i]->gv, chunks[i]->gvs_eh[ft]);
      }
      else { chunks[i]->gvs_eh[ft].push_back(chunks[i]->gv); }
    }

  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine())
      if (chunks[i]->update_eh(ft, skip_w_components)) {
        chunk_connections_valid = false; // E/H allocated - reconnect chunks
        assert(changed_materials);
      }
}

bool fields_chunk::needs_W_prev(component c) const {
  for (susceptibility *chiP = s->chiP[type(c)]; chiP; chiP = chiP->next)
    if (chiP->needs_W_prev()) return true;
  return false;
}

bool fields_chunk::update_eh(field_type ft, bool skip_w_components) {
  field_type ft2 = ft == E_stuff ? D_stuff : B_stuff; // for sources etc.
  bool allocated_eh = false;

  bool have_int_sources = false;
  if (!doing_solve_cw) {
    for (const src_vol &sv : sources[ft2]) {
      if (sv.t()->is_integrated) {
        have_int_sources = true;
        break;
      }
    }
  }

  FOR_FT_COMPONENTS(ft, ec) {
    component dc = field_type_component(ft2, ec);
    DOCMP {
      bool need_fmp = false;
      if (f[ec][cmp]) {
        need_fmp = have_int_sources;
        for (polarization_state *p = pol[ft]; p && !need_fmp; p = p->next)
          need_fmp = need_fmp || p->s->needs_P(ec, cmp, f);
      }
      if (need_fmp) {
        if (!f_minus_p[dc][cmp]) f_minus_p[dc][cmp] = new realnum[gv.ntot()];
      }
      else if (f_minus_p[dc][cmp]) { // remove unneeded f_minus_p
        delete[] f_minus_p[dc][cmp];
        f_minus_p[dc][cmp] = 0;
      }
    }
  }
  bool have_f_minus_p = false;
  FOR_FT_COMPONENTS(ft2, dc) {
    if (f_minus_p[dc][0]) {
      have_f_minus_p = true;
      break;
    }
  }

  const size_t ntot = s->gv.ntot();

  if (have_f_minus_p && doing_solve_cw)
    meep::abort("dispersive materials are not yet implemented for solve_cw");

  //////////////////////////////////////////////////////////////////////////
  // First, initialize f_minus_p to D - P, if necessary

  FOR_FT_COMPONENTS(ft, ec) if (f[ec][0]) {
    component dc = field_type_component(ft2, ec);
    DOCMP if (f_minus_p[dc][cmp]) {
      realnum *fmp = f_minus_p[dc][cmp];
      memcpy(fmp, f[dc][cmp], sizeof(realnum) * ntot);
    }
  }

  for (polarization_state *p = pol[ft]; p; p = p->next)
    if (p->data) p->s->subtract_P(ft, f_minus_p, p->data);

  //////////////////////////////////////////////////////////////////////////
  // Next, subtract time-integrated sources (i.e. polarizations, not currents)

  if (have_f_minus_p && !doing_solve_cw) {
    for (const src_vol &sv : sources[ft2]) {
      if (sv.t()->is_integrated && f[sv.c][0] && ft == type(sv.c)) {
        component c = field_type_component(ft2, sv.c);
        for (size_t j = 0; j < sv.num_points(); ++j) {
          const complex<double> A = sv.dipole(j);
          DOCMP { f_minus_p[c][cmp][sv.index_at(j)] -= (cmp) ? imag(A) : real(A); }
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Finally, compute E = chi1inv * D

  realnum *dmp[NUM_FIELD_COMPONENTS][2];
  FOR_FT_COMPONENTS(ft2, dc) DOCMP2 {
    dmp[dc][cmp] = f_minus_p[dc][cmp] ? f_minus_p[dc][cmp] : f[dc][cmp];
  }

  for (size_t i = 0; i < gvs_eh[ft].size(); ++i) {
    DOCMP FOR_FT_COMPONENTS(ft, ec) {
      if (f[ec][cmp]) {
        if (type(ec) != ft) meep::abort("bug in FOR_FT_COMPONENTS");
        component dc = field_type_component(ft2, ec);
        const direction d_ec = component_direction(ec);
        const ptrdiff_t s_ec = gv.stride(d_ec) * (ft == H_stuff ? -1 : +1);
        const direction d_1 = cycle_direction(gv.dim, d_ec, 1);
        const component dc_1 = direction_component(dc, d_1);
        const ptrdiff_t s_1 = gv.stride(d_1) * (ft == H_stuff ? -1 : +1);
        const direction d_2 = cycle_direction(gv.dim, d_ec, 2);
        const component dc_2 = direction_component(dc, d_2);
        const ptrdiff_t s_2 = gv.stride(d_2) * (ft == H_stuff ? -1 : +1);

        direction dsigw0 = d_ec;
        direction dsigw = s->sigsize[dsigw0] > 1 ? dsigw0 : NO_DIRECTION;

        // lazily allocate any E/H fields that are needed (H==B initially)
        if (i == 0 && f[ec][cmp] == f[dc][cmp] &&
            (s->chi1inv[ec][d_ec] || have_f_minus_p || dsigw != NO_DIRECTION)) {
          f[ec][cmp] = new realnum[gv.ntot()];
          memcpy(f[ec][cmp], f[dc][cmp], gv.ntot() * sizeof(realnum));
          allocated_eh = true;
        }

        // lazily allocate W auxiliary field
        if (i == 0 && !f_w[ec][cmp] && dsigw != NO_DIRECTION) {
          f_w[ec][cmp] = new realnum[gv.ntot()];
          memcpy(f_w[ec][cmp], f[ec][cmp], gv.ntot() * sizeof(realnum));
          if (needs_W_notowned(ec)) allocated_eh = true; // communication needed
        }

        // for solve_cw, when W exists we get W and E from special variables
        if (f_w[ec][cmp] && skip_w_components) continue;

        // save W field from this timestep in f_w_prev if needed by pols
        if (i == 0 && needs_W_prev(ec)) {
          if (!f_w_prev[ec][cmp]) f_w_prev[ec][cmp] = new realnum[gv.ntot()];
          memcpy(f_w_prev[ec][cmp], f_w[ec][cmp] ? f_w[ec][cmp] : f[ec][cmp],
                 sizeof(realnum) * gv.ntot());
        }

        if (f[ec][cmp] != f[dc][cmp])
          STEP_UPDATE_EDHB(f[ec][cmp], ec, gv, gvs_eh[ft][i].little_owned_corner0(ec),
                           gvs_eh[ft][i].big_corner(), dmp[dc][cmp], dmp[dc_1][cmp], dmp[dc_2][cmp],
                           s->chi1inv[ec][d_ec], dmp[dc_1][cmp] ? s->chi1inv[ec][d_1] : NULL,
                           dmp[dc_2][cmp] ? s->chi1inv[ec][d_2] : NULL, s_ec, s_1, s_2, s->chi2[ec],
                           s->chi3[ec], f_w[ec][cmp], dsigw, s->sig[dsigw], s->kap[dsigw]);
      }
    }
  }

  return allocated_eh;
}

} // namespace meep
