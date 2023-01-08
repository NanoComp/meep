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

#include <array>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "meep.hpp"
#include "meep_internals.hpp"

#include "config.h"

#define RESTRICT

using namespace std;

namespace meep {

void fields::step() {
  // however many times the fields have been synched, we want to restore now
  int save_synchronized_magnetic_fields = synchronized_magnetic_fields;
  if (synchronized_magnetic_fields) {
    synchronized_magnetic_fields = 1; // reset synchronization count
    restore_magnetic_fields();
  }

  am_now_working_on(Stepping);

  if (!t) {
    last_step_output_wall_time = wall_time();
    last_step_output_t = t;
  }
  if (verbosity > 0 && wall_time() > last_step_output_wall_time + MEEP_MIN_OUTPUT_TIME) {
    master_printf("on time step %d (time=%g), %g s/step\n", t, time(),
                  (wall_time() - last_step_output_wall_time) / (t - last_step_output_t));
    if (save_synchronized_magnetic_fields)
      master_printf("  (doing expensive timestepping of synched fields)\n");
    last_step_output_wall_time = wall_time();
    last_step_output_t = t;
  }

  phase_material();

  // update cached conductivity-inverse array, if needed
  for (int i = 0; i < num_chunks; i++)
    chunks[i]->s->update_condinv();

  calc_sources(time()); // for B sources
  {
    auto step_timer = with_timing_scope(FieldUpdateB);
    step_db(B_stuff);
  }
  step_source(B_stuff);
  {
    auto step_timer = with_timing_scope(BoundarySteppingB);
    step_boundaries(B_stuff);
  }
  calc_sources(time() + 0.5 * dt); // for integrated H sources
  {
    auto step_timer = with_timing_scope(FieldUpdateH);
    update_eh(H_stuff);
  }
  {
    auto step_timer = with_timing_scope(BoundarySteppingWH);
    step_boundaries(WH_stuff);
  }
  update_pols(H_stuff);
  {
    auto step_timer = with_timing_scope(BoundarySteppingPH);
    step_boundaries(PH_stuff);
  }
  {
    auto step_timer = with_timing_scope(BoundarySteppingH);
    step_boundaries(H_stuff);
  }

  if (fluxes) fluxes->update_half();

  calc_sources(time() + 0.5 * dt); // for D sources
  {
    auto step_timer = with_timing_scope(FieldUpdateD);
    step_db(D_stuff);
  }
  step_source(D_stuff);
  {
    auto step_timer = with_timing_scope(BoundarySteppingD);
    step_boundaries(D_stuff);
  }
  calc_sources(time() + dt); // for integrated E sources
  {
    auto step_timer = with_timing_scope(FieldUpdateE);
    update_eh(E_stuff);
  }
  {
    auto step_timer = with_timing_scope(BoundarySteppingWE);
    step_boundaries(WE_stuff);
  }
  update_pols(E_stuff);
  {
    auto step_timer = with_timing_scope(BoundarySteppingPE);
    step_boundaries(PE_stuff);
  }
  {
    auto step_timer = with_timing_scope(BoundarySteppingE);
    step_boundaries(E_stuff);
  }

  if (fluxes) fluxes->update();
  t += 1;
  update_dfts();
  finished_working();

  // re-synch magnetic fields if they were previously synchronized
  if (save_synchronized_magnetic_fields) {
    synchronize_magnetic_fields();
    synchronized_magnetic_fields = save_synchronized_magnetic_fields;
  }

  changed_materials = false; // any material changes were handled in connect_chunks()

  if (!std::isfinite(get_field(D_EnergyDensity, gv.center(), false)))
    meep::abort("simulation fields are NaN or Inf");
}

void fields::phase_material() {
  bool changed = false;
  if (is_phasing()) {
    CHUNK_OPENMP
    for (int i = 0; i < num_chunks; i++)
      if (chunks[i]->is_mine()) {
        chunks[i]->phase_material(phasein_time);
        changed = changed || chunks[i]->new_s;
      }
    phasein_time--;
    am_now_working_on(MpiAllTime);
    bool changed_mpi = or_to_all(changed);
    finished_working();
    if (changed_mpi) {
      calc_sources(time() + 0.5 * dt); // for integrated H sources
      update_eh(H_stuff);              // ensure H = 1/mu * B
      step_boundaries(H_stuff);
      calc_sources(time() + dt); // for integrated E sources
      update_eh(E_stuff);        // ensure E = 1/eps * D
      step_boundaries(E_stuff);
    }
  }
}

void fields_chunk::phase_material(int phasein_time) {
  if (new_s && phasein_time > 0) {
    changing_structure();
    s->mix_with(new_s, 1.0 / phasein_time);
  }
}

void fields::process_incoming_chunk_data(field_type ft, const chunk_pair &comm_pair) {
  am_now_working_on(Boundaries);
  int this_chunk_idx = comm_pair.second;
  const int pair_idx = chunk_pair_to_index(comm_pair);
  const realnum *pair_comm_block = static_cast<realnum *>(comm_blocks[ft][pair_idx]);

  {
    const comms_key key = {ft, CONNECT_PHASE, comm_pair};
    size_t num_transfers = get_comm_size(key) / 2; // Two realnums per complex
    if (num_transfers) {
      const std::complex<realnum> *pair_comm_block_complex =
          reinterpret_cast<const std::complex<realnum> *>(pair_comm_block);
      const std::vector<realnum *> &incoming_connection =
          chunks[this_chunk_idx]->connections_in.at(key);
      const std::vector<std::complex<realnum> > &connection_phase_for_ft =
          chunks[this_chunk_idx]->connection_phases[key];

      for (size_t n = 0; n < num_transfers; ++n) {
        std::complex<realnum> temp = connection_phase_for_ft[n] * pair_comm_block_complex[n];
        *(incoming_connection[2 * n]) = temp.real();
        *(incoming_connection[2 * n + 1]) = temp.imag();
      }
      pair_comm_block += 2 * num_transfers;
    }
  }

  {
    const comms_key key = {ft, CONNECT_NEGATE, comm_pair};
    const size_t num_transfers = get_comm_size(key);
    if (num_transfers) {
      const std::vector<realnum *> &incoming_connection =
          chunks[this_chunk_idx]->connections_in.at(key);
      for (size_t n = 0; n < num_transfers; ++n) {
        *(incoming_connection[n]) = -pair_comm_block[n];
      }
      pair_comm_block += num_transfers;
    }
  }

  {
    const comms_key key = {ft, CONNECT_COPY, comm_pair};
    const size_t num_transfers = get_comm_size(key);
    if (num_transfers) {
      const std::vector<realnum *> &incoming_connection =
          chunks[this_chunk_idx]->connections_in.at(key);
      for (size_t n = 0; n < num_transfers; ++n) {
        *(incoming_connection[n]) = pair_comm_block[n];
      }
    }
  }
  finished_working();
}

void fields::step_boundaries(field_type ft) {
  connect_chunks(); // re-connect if !chunk_connections_valid

  {
    // Initiate receive operations as early as possible.
    std::unique_ptr<comms_manager> manager = create_comms_manager();

    const auto &sequence = comms_sequence_for_field[ft];
    for (const comms_operation &op : sequence.receive_ops) {
      if (chunks[op.other_chunk_idx]->is_mine()) { continue; }
      chunk_pair comm_pair{op.other_chunk_idx, op.my_chunk_idx};
      comms_manager::receive_callback cb = [this, ft, comm_pair]() {
        process_incoming_chunk_data(ft, comm_pair);
      };
      manager->receive_real_async(comm_blocks[ft][op.pair_idx], static_cast<int>(op.transfer_size),
                                  op.other_proc_id, op.tag, cb);
    }

    // Do the metals first!
    for (int i = 0; i < num_chunks; i++)
      if (chunks[i]->is_mine()) chunks[i]->zero_metal(ft);

    // Copy outgoing data into buffers while following the predefined sequence of comms operations.
    // Trigger the asynchronous send immediately once the outgoing comms buffer has been filled.
    am_now_working_on(Boundaries);

    for (const comms_operation &op : sequence.send_ops) {
      const std::pair<int, int> comm_pair{op.my_chunk_idx, op.other_chunk_idx};
      const int pair_idx = op.pair_idx;

      realnum *outgoing_comm_block = comm_blocks[ft][pair_idx];
      for (connect_phase ip : all_connect_phases) {
        const comms_key key = {ft, ip, comm_pair};
        const size_t pair_comm_size = get_comm_size(key);
        if (pair_comm_size) {
          const std::vector<realnum *> &outgoing_connection =
              chunks[op.my_chunk_idx]->connections_out.at(key);
          for (size_t n = 0; n < pair_comm_size; ++n) {
            outgoing_comm_block[n] = *(outgoing_connection[n]);
          }
          outgoing_comm_block += pair_comm_size;
        }
      }
      if (chunks[op.other_chunk_idx]->is_mine()) { continue; }
      manager->send_real_async(comm_blocks[ft][pair_idx], static_cast<int>(op.transfer_size),
                               op.other_proc_id, op.tag);
    }

    // Process local transfers, which do not depend on a communication mechanism across nodes.
    for (const comms_operation &op : sequence.receive_ops) {
      if (chunks[op.other_chunk_idx]->is_mine()) {
        process_incoming_chunk_data(ft, {op.other_chunk_idx, op.my_chunk_idx});
      }
    }
    finished_working();

    am_now_working_on(MpiOneTime);
    // Let the communication manager drop out of scope to complete all outstanding requests.
    // As data is received, the installed callback handles copying the data from the comm buffer
    // back into the chunk field array.
  }
  finished_working();
}

void fields::step_source(field_type ft, bool including_integrated) {
  if (ft != D_stuff && ft != B_stuff) meep::abort("only step_source(D/B) is okay");
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine()) chunks[i]->step_source(ft, including_integrated);
}

void fields_chunk::step_source(field_type ft, bool including_integrated) {
  if (doing_solve_cw && !including_integrated) return;
  for (const src_vol &sv : sources[ft]) {
    component c = direction_component(first_field_component(ft), component_direction(sv.c));
    const realnum *cndinv = s->condinv[c][component_direction(sv.c)];
    if ((including_integrated || !sv.t()->is_integrated) && f[c][0] &&
        ((ft == D_stuff && is_electric(sv.c)) || (ft == B_stuff && is_magnetic(sv.c)))) {
      if (cndinv)
        for (size_t j = 0; j < sv.num_points(); j++) {
          const ptrdiff_t i = sv.index_at(j);
          const complex<double> A = sv.current(j) * dt * double(cndinv[i]);
          f[c][0][i] -= real(A);
          if (!is_real) f[c][1][i] -= imag(A);
        }
      else
        for (size_t j = 0; j < sv.num_points(); j++) {
          const complex<double> A = sv.current(j) * dt;
          const ptrdiff_t i = sv.index_at(j);
          f[c][0][i] -= real(A);
          if (!is_real) f[c][1][i] -= imag(A);
        }
    }
  }
}

void fields::calc_sources(double tim) {
  for (src_time *s = sources; s; s = s->next)
    s->update(tim, dt);
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine()) chunks[i]->calc_sources(tim);
}

void fields_chunk::calc_sources(double time) {
  (void)time; // unused;
}

} // namespace meep
