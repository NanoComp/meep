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

#include <algorithm>
#include <map>
#include <stdlib.h>
#include <complex>

#include "meep.hpp"
#include "meep/mympi.hpp"
#include "meep_internals.hpp"

#define UNUSED(x) (void)x // silence compiler warnings

using namespace std;

namespace meep {

namespace {

// Creates an optimized comms_sequence from a vector of comms_operations.
// Send operations are prioritized in descending order by the amount of data that is transferred.
comms_sequence optimize_comms_operations(const std::vector<comms_operation> &operations) {
  comms_sequence ret;
  std::map<int, size_t> send_size_by_my_chunk_idx;
  std::map<int, std::vector<comms_operation> > send_ops_by_my_chunk_idx;

  for (const auto &op : operations) {
    if (op.comm_direction == Incoming) {
      ret.receive_ops.push_back(op);
      continue;
    }

    // Group send operations by source chunk and accumulate the transfer size - excluding chunk
    // pairs that reside on the same processor.
    if (op.other_proc_id != my_rank()) {
      send_size_by_my_chunk_idx[op.my_chunk_idx] += op.transfer_size;
    }
    else {
      // Make sure that op.my_chunk_idx is represented in the map.
      send_size_by_my_chunk_idx[op.my_chunk_idx] += 0;
    }
    send_ops_by_my_chunk_idx[op.my_chunk_idx].push_back(op);
  }

  // Sort in descending order to prioritize large transfers.
  std::vector<std::pair<int, size_t> > send_op_sizes(send_size_by_my_chunk_idx.begin(),
                                                     send_size_by_my_chunk_idx.end());
  std::sort(send_op_sizes.begin(), send_op_sizes.end(),
            [](const std::pair<int, size_t> &a, const std::pair<int, size_t> &b) -> bool {
              return a.second > b.second;
            });

  // Assemble send operations.
  for (const auto &size_pair : send_op_sizes) {
    int my_chunk_idx = size_pair.first;
    const auto &ops_vector = send_ops_by_my_chunk_idx[my_chunk_idx];
    ret.send_ops.insert(std::end(ret.send_ops), std::begin(ops_vector), std::end(ops_vector));
  }
  return ret;
}

} // namespace

void fields::set_boundary(boundary_side b, direction d, boundary_condition cond) {
  if (boundaries[b][d] != cond) {
    boundaries[b][d] = cond;
    // we don't need to call sync_chunk_connections() since set_boundary()
    // should always be called on every process
    chunk_connections_valid = false;
  }
}

void fields::use_bloch(direction d, complex<double> kk) {
  k[d] = kk;
  for (int b = 0; b < 2; b++)
    set_boundary(boundary_side(b), d, Periodic);
  if (real(kk) * gv.num_direction(d) == 0.5 * a) // check b.z. edge exactly
    eikna[d] = -exp(-imag(kk) * ((2 * pi / a) * gv.num_direction(d)));
  else {
    const complex<double> I = complex<double>(0.0, 1.0);
    eikna[d] = exp(I * kk * ((2 * pi / a) * gv.num_direction(d)));
  }
  coskna[d] = real(eikna[d]);
  sinkna[d] = imag(eikna[d]);
  if (is_real && kk != 0.0) // FIXME: allow real phases (c.f. CONNECT_PHASE)
    meep::abort("Can't use real fields with bloch boundary conditions!\n");
  chunk_connections_valid = false; // FIXME: we don't always need to invalidate
}

void fields::use_bloch(const vec &k) {
  // Note that I allow a 1D k input when in cylindrical, since in that case
  // it is unambiguous.
  if (k.dim != gv.dim && !(k.dim == D1 && gv.dim == Dcyl))
    meep::abort("Aaaack, k has wrong dimensions!\n");
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    if (gv.has_boundary(Low, d) && d != R) use_bloch(d, k.in_direction(d));
  }
}

ivec fields::ilattice_vector(direction d) const {
  switch (user_volume.dim) {
    case D1: return ivec(2 * user_volume.nz());
    case Dcyl: return iveccyl(0, 2 * user_volume.nz()); // Only Z direction here
    case D2:
      switch (d) {
        case X: return ivec(user_volume.nx() * 2, 0);
        case Y: return ivec(0, user_volume.ny() * 2);
        case Z: // fall-thru
        case R: // fall-thru
        case P: // fall-thru
        case NO_DIRECTION: break;
      }
    case D3:
      switch (d) {
        case X: return ivec(user_volume.nx() * 2, 0, 0);
        case Y: return ivec(0, user_volume.ny() * 2, 0);
        case Z: return ivec(0, 0, user_volume.nz() * 2);
        case R: // fall-thru
        case P: // fall-thru
        case NO_DIRECTION: break;
      }
  }
  meep::abort("Aaack in ilattice_vector.\n");
  return ivec(0);
}

vec fields::lattice_vector(direction d) const { return gv[ilattice_vector(d)]; }

void fields::disconnect_chunks() {
  chunk_connections_valid = false;
  for (int i = 0; i < num_chunks; i++) {
    chunks[i]->connections_in.clear();
    chunks[i]->connections_out.clear();
    chunks[i]->connection_phases.clear();
  }
  FOR_FIELD_TYPES(ft) {
    for (int i = 0; i < num_chunks * num_chunks; i++) {
      delete[] comm_blocks[ft][i];
      comm_blocks[ft][i] = 0;
    }
    comms_sequence_for_field[ft].clear();
  }
  comm_sizes.clear();
}

// this should be called by any code that might set chunk_connections_valid = false,
// with the caveat that we need to be careful that we call it on all processes
void fields::sync_chunk_connections() {
  /* make sure all processes agree on chunk_connections_valid to avoid deadlocks
     when we eventually call connect_chunks */
  am_now_working_on(MpiAllTime);
  chunk_connections_valid = and_to_all(chunk_connections_valid);
  finished_working();
}

void fields::connect_chunks() {
  // might have invalidated connections in step_db, update_eh, or update_pols:
  if (changed_materials) sync_chunk_connections();

  if (!chunk_connections_valid) {
    am_now_working_on(Connecting);
    disconnect_chunks();
    find_metals();
    connect_the_chunks();
    finished_working();
    chunk_connections_valid = true;
  }
}

bool fields::on_metal_boundary(const ivec &here) {
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    if (user_volume.has_boundary(High, d) &&
        here.in_direction(d) == user_volume.big_corner().in_direction(d)) {
      if (boundaries[High][d] == Metallic) return true;
    }
    if (boundaries[Low][d] == Magnetic &&
        here.in_direction(d) == user_volume.little_corner().in_direction(d) + 1)
      return true;
    if (boundaries[Low][d] == Metallic &&
        here.in_direction(d) == user_volume.little_corner().in_direction(d))
      return true;
  }
  return false;
}

bool fields::locate_point_in_user_volume(ivec *there, complex<double> *phase) const {
  // Check if a translational symmetry is needed to bring the point in...
  if (!user_volume.owns(*there)) {
    LOOP_OVER_DIRECTIONS(gv.dim, d) {
      if (boundaries[High][d] == Periodic &&
          there->in_direction(d) <= user_volume.little_corner().in_direction(d)) {
        while (there->in_direction(d) <= user_volume.little_corner().in_direction(d)) {
          *there += ilattice_vector(d);
          *phase *= conj(eikna[d]);
        }
      }
      else if (boundaries[High][d] == Periodic &&
               there->in_direction(d) - ilattice_vector(d).in_direction(d) >
                   user_volume.little_corner().in_direction(d)) {
        while (there->in_direction(d) - ilattice_vector(d).in_direction(d) >
               user_volume.little_corner().in_direction(d)) {
          *there -= ilattice_vector(d);
          *phase *= eikna[d];
        }
      }
    }
  }
  return user_volume.owns(*there);
}

void fields::locate_volume_source_in_user_volume(const vec p1, const vec p2, vec newp1[8],
                                                 vec newp2[8], complex<double> kphase[8],
                                                 int &ncopies) const {
  // For periodic boundary conditions,
  // this function locates up to 8 translated copies of the initial grid_volume specified by (p1,p2)
  // First bring center of grid_volume inside
  ncopies = 1;
  newp1[0] = p1;
  newp2[0] = p2;
  kphase[0] = 1;
  vec cen = (newp1[0] + newp2[0]) * 0.5;
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    if (boundaries[High][d] == Periodic) {
      while (cen.in_direction(d) < gv.boundary_location(Low, d)) {
        newp1[0] += lattice_vector(d);
        newp2[0] += lattice_vector(d);
        kphase[0] *= conj(eikna[d]);
        cen = (newp1[0] + newp2[0]) * 0.5;
      }
      while (cen.in_direction(d) > gv.boundary_location(High, d)) {
        newp1[0] -= lattice_vector(d);
        newp2[0] -= lattice_vector(d);
        kphase[0] *= eikna[d];
        cen = (newp1[0] + newp2[0]) * 0.5;
      }
    }
  }

  // if grid_volume extends outside user_volume in any direction, we need to duplicate already
  // existing copies
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    if (boundaries[High][d] == Periodic) {
      if (newp1[0].in_direction(d) < gv.boundary_location(Low, d) ||
          newp2[0].in_direction(d) < gv.boundary_location(Low, d)) {
        for (int j = 0; j < ncopies; j++) {
          newp1[ncopies + j] = newp1[j] + lattice_vector(d);
          newp2[ncopies + j] = newp2[j] + lattice_vector(d);
          kphase[ncopies + j] = kphase[j] * conj(eikna[d]);
        }
        ncopies *= 2;
      }
      else if (newp1[0].in_direction(d) > gv.boundary_location(High, d) ||
               newp2[0].in_direction(d) > gv.boundary_location(High, d)) {
        for (int j = 0; j < ncopies; j++) {
          newp1[ncopies + j] = newp1[j] - lattice_vector(d);
          newp2[ncopies + j] = newp2[j] - lattice_vector(d);
          kphase[ncopies + j] = kphase[j] * eikna[d];
        }
        ncopies *= 2;
      }
    }
  }
}

bool fields::locate_component_point(component *c, ivec *there, complex<double> *phase) const {
  // returns true if this point and component exist in the user_volume.  If
  // that is the case, on return *c and *there store the component and
  // location of where the point actually is, and *phase determines holds
  // the phase needed to get the true field.  If the point is not located,
  // *c and *there will hold undefined values.

  // Check if nothing tricky is needed...
  *phase = 1.0;
  if (!locate_point_in_user_volume(there, phase)) return false;
  // Check if a rotation or inversion brings the point in...
  if (user_volume.owns(*there))
    for (int sn = 0; sn < S.multiplicity(); sn++) {
      const ivec here = S.transform(*there, sn);
      if (gv.owns(here)) {
        *there = here;
        *phase *= S.phase_shift(*c, sn);
        *c = direction_component(*c, S.transform(component_direction(*c), sn).d);
        return true;
      }
    }
  return false;
}

void fields_chunk::zero_metal(field_type ft) {
  for (size_t i = 0; i < num_zeroes[ft]; i++)
    *(zeroes[ft][i]) = 0.0;
}

void fields::find_metals() {
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine()) {
      const grid_volume vi = chunks[i]->gv;
      FOR_FIELD_TYPES(ft) {
        delete[] chunks[i]->zeroes[ft];
        // First electric components...
        chunks[i]->num_zeroes[ft] = 0;
        DOCMP FOR_COMPONENTS(c) {
          if (type(c) == ft && chunks[i]->f[c][cmp]) LOOP_OVER_VOL_OWNED(vi, c, n) {
              if (IVEC_LOOP_AT_BOUNDARY) { // todo: just loop over boundaries
                IVEC_LOOP_ILOC(vi, here);
                if (on_metal_boundary(here)) chunks[i]->num_zeroes[ft]++;
              }
            }
        }
        typedef realnum *realnum_ptr;
        chunks[i]->zeroes[ft] = new realnum_ptr[chunks[i]->num_zeroes[ft]];
        size_t num = 0;
        DOCMP FOR_COMPONENTS(c) {
          if (type(c) == ft && chunks[i]->f[c][cmp]) LOOP_OVER_VOL_OWNED(vi, c, n) {
              if (IVEC_LOOP_AT_BOUNDARY) { // todo: just loop over boundaries
                IVEC_LOOP_ILOC(vi, here);
                if (on_metal_boundary(here))
                  chunks[i]->zeroes[ft][num++] = chunks[i]->f[c][cmp] + n;
              }
            }
        }
      }
    }
}

bool fields_chunk::needs_W_notowned(component c) {
  for (susceptibility *chiP = s->chiP[type(c)]; chiP; chiP = chiP->next)
    if (chiP->needs_W_notowned(c, f)) return true;
  return false;
}

void fields::connect_the_chunks() {
  /* For some of the chunks, H==B, and we definitely don't need to
     send B between two such chunks.   We'll still send B when
     the recipient has H != B, since the recipient needs to get B
     from somewhere (although it could get it locally, in principle).
     When the sender has H != B, we'll skip sending B (we'll only send H)
     since we need to get the correct curl H in the E update.  This is
     a bit subtle since the non-owned B may be different from H even
     on an H==B chunk (true?), but since we don't use the non-owned B
     for anything(?) it shouldn't matter. */
  std::vector<int> B_redundant(num_chunks * 2 * 5);
  for (int i = 0; i < num_chunks; ++i)
    FOR_H_AND_B(hc, bc) {
      B_redundant[5 * (num_chunks + i) + bc - Bx] = chunks[i]->f[hc][0] == chunks[i]->f[bc][0];
    }
  am_now_working_on(MpiAllTime);
  and_to_all(B_redundant.data() + 5 * num_chunks, B_redundant.data(), 5 * num_chunks);
  finished_working();

  /* Figure out whether we need the notowned W field (== E/H in
     non-PML regions) in update_pols, e.g. if we have an anisotropic
     susceptibility.  In this case, we have an additional
     communication step where we communicate the notowned W.  Then,
     after updating the polarizations, we communicate the notowned E/H
     ... this does the E/H communication twice between non-PML regions
     and hence is somewhat wasteful, but greatly simplifies the case
     of communicating between a PML region (which has a separate W
     array) and a non-PML region (no separate W). */
  bool needs_W_notowned[NUM_FIELD_COMPONENTS];
  FOR_COMPONENTS(c) { needs_W_notowned[c] = false; }
  FOR_E_AND_H(c) {
    for (int i = 0; i < num_chunks; i++)
      needs_W_notowned[c] = needs_W_notowned[c] || chunks[i]->needs_W_notowned(c);
  }
  am_now_working_on(MpiAllTime);
  FOR_E_AND_H(c) { needs_W_notowned[c] = or_to_all(needs_W_notowned[c]); }
  finished_working();

  comm_sizes.clear();
  const size_t num_reals_per_voxel = is_real ? 1 : 2;
  for (int i = 0; i < num_chunks; i++) {
    // First count the border elements...
    const grid_volume vi = chunks[i]->gv;
    FOR_COMPONENTS(corig) {
      if (have_component(corig)) LOOP_OVER_VOL_NOTOWNED(vi, corig, n) {
          IVEC_LOOP_ILOC(vi, here);
          component c = corig;
          // We're looking at a border element...
          complex<double> thephase;
          if (locate_component_point(&c, &here, &thephase) && !on_metal_boundary(here))
            for (int j = 0; j < num_chunks; j++) {
              const std::pair<int, int> pair_j_to_i{j, i};
              if ((chunks[i]->is_mine() || chunks[j]->is_mine()) && chunks[j]->gv.owns(here) &&
                  !(is_B(corig) && is_B(c) && B_redundant[5 * i + corig - Bx] &&
                    B_redundant[5 * j + c - Bx])) {
                const connect_phase ip = thephase == 1.0
                                             ? CONNECT_COPY
                                             : (thephase == -1.0 ? CONNECT_NEGATE : CONNECT_PHASE);
                comm_sizes[{type(c), ip, pair_j_to_i}] += num_reals_per_voxel;

                if (needs_W_notowned[corig]) {
                  field_type f = is_electric(corig) ? WE_stuff : WH_stuff;
                  comm_sizes[{f, ip, pair_j_to_i}] += num_reals_per_voxel;
                }
                if (is_electric(corig) || is_magnetic(corig)) {
                  field_type f = is_electric(corig) ? PE_stuff : PH_stuff;
                  size_t ni = 0, cni = 0;
                  for (polarization_state *pi = chunks[i]->pol[type(corig)]; pi; pi = pi->next)
                    for (polarization_state *pj = chunks[j]->pol[type(c)]; pj; pj = pj->next)
                      if (*pi->s == *pj->s) {
                        if (pi->data && chunks[i]->is_mine()) {
                          ni += pi->s->num_internal_notowned_needed(corig, pi->data);
                          cni += pi->s->num_cinternal_notowned_needed(corig, pi->data);
                        }
                        else if (pj->data && chunks[j]->is_mine()) {
                          ni += pj->s->num_internal_notowned_needed(c, pj->data);
                          cni += pj->s->num_cinternal_notowned_needed(c, pj->data);
                        }
                      }
                  comm_sizes[{f, ip, pair_j_to_i}] += cni * num_reals_per_voxel;
                  comm_sizes[{f, CONNECT_COPY, pair_j_to_i}] += ni;
                }
              } // if is_mine and owns...
            }   // loop over j chunks
        }       // LOOP_OVER_VOL_NOTOWNED
    }           // FOR_COMPONENTS

    // Allocating comm blocks as we go...
    FOR_FIELD_TYPES(ft) {
      for (int j = 0; j < num_chunks; j++) {
        delete[] comm_blocks[ft][j + i * num_chunks];
        comm_blocks[ft][j + i * num_chunks] = new realnum[comm_size_tot(ft, {j, i})];
      }
    }
  } // loop over i chunks

  // Preallocate all connection vectors.
  for (const std::pair<const comms_key, size_t> &key_and_comm_size : comm_sizes) {
    const chunk_pair &pair_j_to_i = key_and_comm_size.first.pair;
    if (chunks[pair_j_to_i.first]->is_mine()) {
      chunks[pair_j_to_i.first]->connections_out[key_and_comm_size.first].reserve(
          key_and_comm_size.second);
    }
    if (chunks[pair_j_to_i.second]->is_mine()) {
      chunks[pair_j_to_i.second]->connections_in[key_and_comm_size.first].reserve(
          key_and_comm_size.second);
    }
  }

  // Next start setting up the connections...
  for (int i = 0; i < num_chunks; i++) {
    const grid_volume &vi = chunks[i]->gv;

    FOR_COMPONENTS(corig) {
      if (have_component(corig)) LOOP_OVER_VOL_NOTOWNED(vi, corig, n) {
          IVEC_LOOP_ILOC(vi, here);
          component c = corig;
          // We're looking at a border element...
          std::complex<double> thephase_double;
          if (locate_component_point(&c, &here, &thephase_double) && !on_metal_boundary(here)) {
            std::complex<realnum> thephase(thephase_double.real(), thephase_double.imag());
            for (int j = 0; j < num_chunks; j++) {
              const std::pair<int, int> pair_j_to_i{j, i};
              const bool i_is_mine = chunks[i]->is_mine();
              const bool j_is_mine = chunks[j]->is_mine();
              if (!i_is_mine && !j_is_mine) { continue; }

              auto push_back_phase = [this, &thephase, &pair_j_to_i](field_type f) {
                chunks[pair_j_to_i.second]
                    ->connection_phases[{f, CONNECT_PHASE, pair_j_to_i}]
                    .push_back(std::complex<realnum>(thephase.real(), thephase.imag()));
              };
              auto push_back_incoming_pointer = [this, &pair_j_to_i](field_type f, connect_phase ip,
                                                                     realnum *p) {
                chunks[pair_j_to_i.second]->connections_in[{f, ip, pair_j_to_i}].push_back(p);
              };
              auto push_back_outgoing_pointer = [this, &pair_j_to_i](field_type f, connect_phase ip,
                                                                     realnum *p) {
                chunks[pair_j_to_i.first]->connections_out[{f, ip, pair_j_to_i}].push_back(p);
              };

              if (chunks[j]->gv.owns(here) &&
                  !(is_B(corig) && is_B(c) && B_redundant[5 * i + corig - Bx] &&
                    B_redundant[5 * j + c - Bx])) {
                const connect_phase ip =
                    thephase == static_cast<realnum>(1.0)
                        ? CONNECT_COPY
                        : (thephase == static_cast<realnum>(-1.0) ? CONNECT_NEGATE : CONNECT_PHASE);
                const ptrdiff_t m = chunks[j]->gv.index(c, here);

                {
                  field_type f = type(c);
                  if (i_is_mine) {
                    if (ip == CONNECT_PHASE) { push_back_phase(f); }
                    DOCMP { push_back_incoming_pointer(f, ip, chunks[i]->f[corig][cmp] + n); }
                  }
                  if (j_is_mine) {
                    DOCMP { push_back_outgoing_pointer(f, ip, chunks[j]->f[c][cmp] + m); }
                  }
                }

                if (needs_W_notowned[corig]) {
                  field_type f = is_electric(corig) ? WE_stuff : WH_stuff;
                  if (i_is_mine) {
                    if (ip == CONNECT_PHASE) { push_back_phase(f); }
                    DOCMP {
                      push_back_incoming_pointer(f, ip,
                                                 (chunks[i]->f_w[corig][cmp]
                                                      ? chunks[i]->f_w[corig][cmp]
                                                      : chunks[i]->f[corig][cmp]) +
                                                     n);
                    }
                  }
                  if (j_is_mine) {
                    DOCMP {
                      push_back_outgoing_pointer(
                          f, ip,
                          (chunks[j]->f_w[c][cmp] ? chunks[j]->f_w[c][cmp] : chunks[j]->f[c][cmp]) +
                              m);
                    }
                  }
                }

                if (is_electric(corig) || is_magnetic(corig)) {
                  field_type f = is_electric(corig) ? PE_stuff : PH_stuff;
                  for (polarization_state *pi = chunks[i]->pol[type(corig)]; pi; pi = pi->next)
                    for (polarization_state *pj = chunks[j]->pol[type(c)]; pj; pj = pj->next)
                      if (*pi->s == *pj->s) {
                        polarization_state *po = NULL;
                        if (pi->data && chunks[i]->is_mine())
                          po = pi;
                        else if (pj->data && chunks[j]->is_mine())
                          po = pj;
                        if (po) {
                          const connect_phase iip = CONNECT_COPY;
                          const size_t ni = po->s->num_internal_notowned_needed(corig, po->data);
                          for (size_t k = 0; k < ni; ++k) {
                            if (i_is_mine) {
                              push_back_incoming_pointer(
                                  f, iip, po->s->internal_notowned_ptr(k, corig, n, pi->data));
                            }
                            if (j_is_mine) {
                              push_back_outgoing_pointer(
                                  f, iip, po->s->internal_notowned_ptr(k, c, m, pj->data));
                            }
                          }
                          const size_t cni = po->s->num_cinternal_notowned_needed(corig, po->data);
                          for (size_t k = 0; k < cni; ++k) {
                            if (i_is_mine) {
                              if (ip == CONNECT_PHASE) { push_back_phase(f); }

                              DOCMP {
                                push_back_incoming_pointer(
                                    f, ip,
                                    po->s->cinternal_notowned_ptr(k, corig, cmp, n, pi->data));
                              }
                            }
                            if (j_is_mine) {
                              DOCMP {
                                push_back_outgoing_pointer(
                                    f, ip, po->s->cinternal_notowned_ptr(k, c, cmp, m, pj->data));
                              }
                            }
                          }
                        }
                      }
                } // is_electric(corig)
              }   // if is_mine and owns...
            }     // loop over j chunks
          }       // here in user_volume
        }         // LOOP_OVER_VOL_NOTOWNED
    }             // FOR_COMPONENTS
  }               // loop over i chunks

  FOR_FIELD_TYPES(f) {

    // Calculate the sequence of sends and receives in advance.
    // Initiate receive operations as early as possible.
    std::unique_ptr<comms_manager> manager = create_comms_manager();
    std::vector<comms_operation> operations;
    std::vector<int> tagto(count_processors());

    for (int j = 0; j < num_chunks; j++) {
      for (int i = 0; i < num_chunks; i++) {
        const chunk_pair pair{j, i};
        const size_t comm_size = comm_size_tot(f, pair);
        if (!comm_size) continue;
        if (comm_size > manager->max_transfer_size()) {
          // MPI uses int for size to send/recv
          meep::abort("communications size too big for the current implementation");
        }
        const int pair_idx = j + i * num_chunks;

        if (chunks[j]->is_mine()) {
          operations.push_back(comms_operation{/*my_chunk_idx=*/j,
                                               /*other_chunk_idx=*/i,
                                               /*other_proc_id=*/chunks[i]->n_proc(),
                                               /*pair_idx=*/pair_idx,
                                               /*transfer_size=*/comm_size,
                                               /*comm_direction=*/Outgoing,
                                               /*tag=*/tagto[chunks[i]->n_proc()]++});
        }
        if (chunks[i]->is_mine()) {
          operations.push_back(comms_operation{/*my_chunk_idx=*/i,
                                               /*other_chunk_idx=*/j,
                                               /*other_proc_id=*/chunks[j]->n_proc(),
                                               /*pair_idx=*/pair_idx,
                                               /*transfer_size=*/comm_size,
                                               /*comm_direction=*/Incoming,
                                               /*tag=*/tagto[chunks[j]->n_proc()]++});
        }
      }
    }

    comms_sequence_for_field[f] = optimize_comms_operations(operations);
  }
}

} // namespace meep
