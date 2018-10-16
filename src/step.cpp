/* Copyright (C) 2005-2015 Massachusetts Institute of Technology
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

#include "multithreading.hpp"

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
  if (!quiet && wall_time() > last_step_output_wall_time + MIN_OUTPUT_TIME) {
    master_printf("on time step %d (time=%g), %g s/step\n", t, time(),
		  (wall_time() - last_step_output_wall_time) /
		  (t - last_step_output_t));
    if (save_synchronized_magnetic_fields)
      master_printf("  (doing expensive timestepping of synched fields)\n");
    last_step_output_wall_time = wall_time();
    last_step_output_t = t;
print_benchmarks();
  }

checkpoint("phase_material");
  phase_material();
checkpoint("calc_sources");

  // update cached conductivity-inverse array, if needed
  for (int i=0;i<num_chunks;i++) chunks[i]->s->update_condinv();
  calc_sources(time()); // for B sources
checkpoint("step_db B");
  step_db(B_stuff);
checkpoint("step_source B");
  step_source(B_stuff);
checkpoint("step_boundaries B");
  step_boundaries(B_stuff);
checkpoint("calc_sources 2");
  calc_sources(time() + 0.5*dt); // for integrated H sources
checkpoint("update_eh H");
  update_eh(H_stuff);
checkpoint("step_boundaries WH");
  step_boundaries(WH_stuff);
checkpoint("update_pols H");
  update_pols(H_stuff);
checkpoint("step_boundaries PH");
  step_boundaries(PH_stuff);
checkpoint("step_boundaries H");
  step_boundaries(H_stuff);

checkpoint("fluxes 1");
  if (fluxes) fluxes->update_half();

checkpoint("calc_sources 3");
  calc_sources(time() + 0.5*dt); // for D sources
checkpoint("step_db D");
  step_db(D_stuff);
checkpoint("step_source D");
  step_source(D_stuff);
checkpoint("step_boundaries D");
  step_boundaries(D_stuff);
checkpoint("calc_sources 4");
  calc_sources(time() + dt); // for integrated E sources
checkpoint("update_eh E");
  update_eh(E_stuff);
checkpoint("step_boundaries WE");
  step_boundaries(WE_stuff);
checkpoint("update_pols E");
  update_pols(E_stuff);
checkpoint("step_boundaries PE");
  step_boundaries(PE_stuff);
checkpoint("step_boundaries E");
  step_boundaries(E_stuff);

checkpoint("fluxes 2");
  if (fluxes) fluxes->update();
checkpoint("dft");
  t += 1;
  update_dfts();
checkpoint();

  finished_working();

  // re-synch magnetic fields if they were previously synchronized
  if (save_synchronized_magnetic_fields) {
    synchronize_magnetic_fields();
    synchronized_magnetic_fields = save_synchronized_magnetic_fields;
  }
}

double fields_chunk::peek_field(component c, const vec &where) {
  double w[8];
  ivec ilocs[8];
  gv.interpolate(c,where, ilocs, w);
  if (gv.contains(ilocs[0]) && f[c][0]) {
    double hello = 0.0;
    if (is_mine()) hello = f[c][0][gv.index(c,ilocs[0])];
    broadcast(n_proc(), &hello, 1);
    return hello;
  }
  //abort("Got no such %s field at %g %g!\n",
  //      component_name(c), gv[ilocs[0]].x(), gv[ilocs[0]].y());
  return 0.0;
}

void fields::phase_material() {
  bool changed = false;
  if (is_phasing()) {
    for (int i=0;i<num_chunks;i++)
      if (chunks[i]->is_mine()) {
      	chunks[i]->phase_material(phasein_time);
      	changed = changed || chunks[i]->new_s;
      }
    phasein_time--;
  }
  if (or_to_all(changed)) {
    calc_sources(time() + 0.5*dt); // for integrated H sources
    update_eh(H_stuff); // ensure H = 1/mu * B
    step_boundaries(H_stuff);
    calc_sources(time() + dt); // for integrated E sources
    update_eh(E_stuff); // ensure E = 1/eps * D
    step_boundaries(E_stuff);
  }
}

void fields_chunk::phase_material(int phasein_time) {
  if (new_s && phasein_time > 0) {
    changing_structure();
    s->mix_with(new_s, 1.0/phasein_time);
  }
}

void fields::step_boundaries(field_type ft) {
  connect_chunks(); // re-connect if !chunk_connections_valid

  am_now_working_on(MpiTime);

  // Do the metals first!
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine()) chunks[i]->zero_metal(ft);

  /* Note that the copying of data to/from buffers is order-sensitive,
     and must be kept consistent with the code in boundaries.cpp.
     In particular, we require that boundaries.cpp set up the connections
     array so that all of the connections for process i come before all
     of the connections for process i' for i < i'  */

  // First copy outgoing data to buffers...
  for (int j=0;j<num_chunks;j++)
    if (chunks[j]->is_mine()) {
      int wh[3] = {0,0,0};
      for (int i=0;i<num_chunks;i++) {
      	const int pair = j+i*num_chunks;
      	size_t n0 = 0;
      	for (int ip=0;ip<3;ip++) {
      	  for (size_t n=0;n<comm_sizes[ft][ip][pair];n++)
      	    comm_blocks[ft][pair][n0 + n] =
      	      *(chunks[j]->connections[ft][ip][Outgoing][wh[ip]++]);
      	  n0 += comm_sizes[ft][ip][pair];
      	}
      }
    }
  boundary_communications(ft);

  // Finally, copy incoming data to the fields themselves, multiplying phases:
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine()) {
      int wh[3] = {0,0,0};
      for (int j=0;j<num_chunks;j++) {
        const int pair = j+i*num_chunks;
      	connect_phase ip = CONNECT_PHASE;
        for (size_t n = 0; n < comm_sizes[ft][ip][pair]; n += 2, wh[ip] += 2) {
          const double phr = real(chunks[i]->connection_phases[ft][wh[ip]/2]);
          const double phi = imag(chunks[i]->connection_phases[ft][wh[ip]/2]);
          *(chunks[i]->connections[ft][ip][Incoming][wh[ip]]) =
            phr*comm_blocks[ft][pair][n] - phi*comm_blocks[ft][pair][n+1];
          *(chunks[i]->connections[ft][ip][Incoming][wh[ip]+1]) =
            phr*comm_blocks[ft][pair][n+1] + phi*comm_blocks[ft][pair][n];
        }
      	size_t n0 = comm_sizes[ft][ip][pair];
      	ip = CONNECT_NEGATE;
        for (size_t n = 0; n < comm_sizes[ft][ip][pair]; ++n)
          *(chunks[i]->connections[ft][ip][Incoming][wh[ip]++])
      	    = -comm_blocks[ft][pair][n0 + n];
      	n0 += comm_sizes[ft][ip][pair];
      	ip = CONNECT_COPY;
        for (size_t n = 0; n < comm_sizes[ft][ip][pair]; ++n)
          *(chunks[i]->connections[ft][ip][Incoming][wh[ip]++])
      	    = comm_blocks[ft][pair][n0 + n];
      }
    }

  finished_working();
}

void fields::step_source(field_type ft, bool including_integrated) {
  if (ft != D_stuff && ft != B_stuff) abort("only step_source(D/B) is okay");
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_source(ft, including_integrated);
}
void fields_chunk::step_source(field_type ft, bool including_integrated) {
  if (doing_solve_cw && !including_integrated) return;
  for (src_vol *sv = sources[ft]; sv; sv = sv->next) {
    component c = direction_component(first_field_component(ft),
				      component_direction(sv->c));
    const realnum *cndinv = s->condinv[c][component_direction(sv->c)];
    if ((including_integrated || !sv->t->is_integrated)	&& f[c][0]
	&& ((ft == D_stuff && is_electric(sv->c))
	    || (ft == B_stuff && is_magnetic(sv->c)))) {
      if (cndinv)
	for (size_t j=0; j<sv->npts; j++) {
	  const ptrdiff_t i = sv->index[j];
	  const complex<double> A = sv->current(j) * dt * double(cndinv[i]);
	  f[c][0][i] -= real(A);
	  if (!is_real) f[c][1][i] -= imag(A);
	}
      else
	for (size_t j=0; j<sv->npts; j++) {
	  const complex<double> A = sv->current(j) * dt;
	  const ptrdiff_t i = sv->index[j];
	  f[c][0][i] -= real(A);
	  if (!is_real) f[c][1][i] -= imag(A);
	}
    }
  }
}

void fields::calc_sources(double tim) {
  for (src_time *s = sources; s; s = s->next) s->update(tim, dt);
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->calc_sources(tim);
}

void fields_chunk::calc_sources(double time) {
  (void) time; // unused;
}

} // namespace meep
