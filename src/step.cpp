/* Copyright (C) 2003 Massachusetts Institute of Technology
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

#include "meep.h"
#include "meep_internals.h"

#include "config.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

#define RESTRICT

namespace meep {

void fields::step() {
  am_now_working_on(Stepping);
  phase_material();

  calc_sources(time() - 0.5 * dt); // for H sources

  //for (int i=0;i<num_chunks;i++)
  //  master_printf("Field is now %g\n", chunks[i]->peek_field(Ex,vec2d(1.55,0.6)));
  step_h();
  step_h_source();
  step_boundaries(H_stuff);
  // because step_boundaries overruns the timing stack...
  am_now_working_on(Stepping);

  if (fluxes) fluxes->update_half();

  calc_sources(time()); // for E sources

  step_d();
  step_boundaries(D_stuff);

  update_e_from_d();
  step_boundaries(E_stuff);

  // because step_boundaries overruns the timing stack...
  am_now_working_on(Stepping);

  update_from_e();
  step_boundaries(P_stuff);

  if (fluxes) fluxes->update();
  t += 1;
  update_dfts();
  finished_working();
}

double fields_chunk::peek_field(component c, const vec &where) {
  double w[8];
  ivec ilocs[8];
  v.interpolate(c,where, ilocs, w);
  if (v.contains(ilocs[0]) && f[c][0]) {
    double hello = 0.0;
    if (is_mine()) hello = f[c][0][v.index(c,ilocs[0])];
    broadcast(n_proc(), &hello, 1);
    return hello;
  }
  //abort("Got no such %s field at %g %g!\n",
  //      component_name(c), v[ilocs[0]].x(), v[ilocs[0]].y());
  return 0.0;
}

void fields::phase_material() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->phase_material(phasein_time);
  phasein_time--;
}

void fields_chunk::phase_material(int phasein_time) {
  if (new_s && phasein_time > 0) {
    // First multiply the electric field by epsilon...
    DOCMP {
      FOR_ELECTRIC_COMPONENTS(c)
        if (f[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f[c][cmp][i] /= s->inveps[c][component_direction(c)][i];
      FOR_ELECTRIC_COMPONENTS(c)
        if (f_p_pml[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f_p_pml[c][cmp][i] /= s->inveps[c][component_direction(c)][i];
      FOR_ELECTRIC_COMPONENTS(c)
        if (f_m_pml[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f_m_pml[c][cmp][i] /= s->inveps[c][component_direction(c)][i];
    }
    // Then change epsilon...
    s->mix_with(new_s, 1.0/phasein_time);
    // Finally divide the electric field by the new epsilon...
    DOCMP {
      FOR_ELECTRIC_COMPONENTS(c)
        if (f[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f[c][cmp][i] *= s->inveps[c][component_direction(c)][i];
      FOR_ELECTRIC_COMPONENTS(c)
        if (f_p_pml[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f_p_pml[c][cmp][i] *= s->inveps[c][component_direction(c)][i];
      FOR_ELECTRIC_COMPONENTS(c)
        if (f_m_pml[c][cmp])
          for (int i=0;i<v.ntot();i++)
            f_m_pml[c][cmp][i] *= s->inveps[c][component_direction(c)][i];
    }
  }
}

void fields::step_boundaries(field_type ft) {
  connect_chunks(); // re-connect if !chunk_connections_valid
  am_now_working_on(MpiTime);
  // Do the metals first!
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine()) chunks[i]->zero_metal(ft);
  // First copy outgoing data to buffers...
  int *wh = new int[num_chunks];
  for (int i=0;i<num_chunks;i++) wh[i] = 0;
  for (int i=0;i<num_chunks;i++)
    for (int j=0;j<num_chunks;j++)
      if (chunks[j]->is_mine()) {
        const int pair = j+i*num_chunks;
        for (int n=0;n<comm_sizes[ft][pair];n++)
          comm_blocks[ft][pair][n] =
            *(chunks[j]->connections[ft][Outgoing][wh[j]++]);
      }
  // Communicate the data around!
#if 0 // This is the blocking version, which should always be safe!
  for (int noti=0;noti<num_chunks;noti++)
    for (int j=0;j<num_chunks;j++) {
      const int i = (noti+j)%num_chunks;
      const int pair = j+i*num_chunks;
      DOCMP {
        send(chunks[j]->n_proc(), chunks[i]->n_proc(),
             comm_blocks[ft][pair], comm_sizes[ft][pair]);
      }
    }
#endif
#ifdef HAVE_MPI
  const int maxreq = num_chunks*num_chunks;
  MPI_Request *reqs = new MPI_Request[maxreq];
  MPI_Status *stats = new MPI_Status[maxreq];
  int reqnum = 0;
  int *tagto = new int[count_processors()];
  for (int i=0;i<count_processors();i++) tagto[i] = 0;
  for (int noti=0;noti<num_chunks;noti++)
    for (int j=0;j<num_chunks;j++) {
      const int i = (noti+j)%num_chunks;
      if (chunks[j]->n_proc() != chunks[i]->n_proc()) {
        const int pair = j+i*num_chunks;
        if (comm_sizes[ft][pair] > 0) {
          if (chunks[j]->is_mine())
            MPI_Isend(comm_blocks[ft][pair], comm_sizes[ft][pair],
                      MPI_DOUBLE, chunks[i]->n_proc(),
                      tagto[chunks[i]->n_proc()]++,
                      MPI_COMM_WORLD, &reqs[reqnum++]);
          if (chunks[i]->is_mine())
            MPI_Irecv(comm_blocks[ft][pair], comm_sizes[ft][pair],
                      MPI_DOUBLE, chunks[j]->n_proc(),
                      tagto[chunks[j]->n_proc()]++,
                      MPI_COMM_WORLD, &reqs[reqnum++]);
        }
      }
    }
  delete[] tagto;
  if (reqnum > maxreq) abort("Too many requests!!!\n");
  if (reqnum > 0) MPI_Waitall(reqnum, reqs, stats);
  delete[] reqs;
  delete[] stats;
#endif
  
  // Finally, copy incoming data to the fields themselves!
  for (int i=0;i<num_chunks;i++) {
    int wh = 0;
    if (chunks[i]->is_mine())
      for (int j=0;j<num_chunks;j++) {
        const int pair = j+i*num_chunks;
        int n;
        for (n=0;n<comm_num_complex[ft][pair];n+=2) {
          const double phr = real(chunks[i]->connection_phases[ft][wh/2]);
          const double phi = imag(chunks[i]->connection_phases[ft][wh/2]);
          *(chunks[i]->connections[ft][Incoming][wh]) =
            phr*comm_blocks[ft][pair][n] - phi*comm_blocks[ft][pair][n+1];
          *(chunks[i]->connections[ft][Incoming][wh+1]) =
            phr*comm_blocks[ft][pair][n+1] + phi*comm_blocks[ft][pair][n];
          wh += 2;
        }
        for (;n<comm_num_complex[ft][pair]+comm_num_negate[ft][pair];n++)
          *(chunks[i]->connections[ft][Incoming][wh++]) = -comm_blocks[ft][pair][n];
        for (;n<comm_sizes[ft][pair];n++)
          *(chunks[i]->connections[ft][Incoming][wh++]) = comm_blocks[ft][pair][n];
      }
  }
  delete[] wh;
  finished_working();
}

void fields::step_h_source() {
  const double tim = time();
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->step_h_source(chunks[i]->h_sources, tim);
}

void fields_chunk::step_h_source(src_vol *sv, double time) {
  if (sv == NULL) return;
  component c = sv->c;
  if (f[c][0] && is_magnetic(c))
    for (int j=0; j<sv->npts; j++) {
      const complex<double> A = sv->current(j);
      const int i = sv->index[j];
      f[c][0][i] += real(A);
      if (!is_real) f[c][1][i] += imag(A);
    }
  step_h_source(sv->next, time);
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
