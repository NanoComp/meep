/* Copyright (C) 2004 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "meep.h"
#include "meep_internals.h"
#include "h5io.h"

namespace meep {

struct dft_chunk_data { // for passing to field::integrate as void*
  double omega_min, domega;
  int Nomega;
  component c;
  dft_chunk *dft_chunks;
};

dft_chunk::dft_chunk(fields_chunk *fc_,
		     ivec is_, ivec ie_,
		     vec s0_, vec s1_, vec e0_, vec e1_,
		     double dV0_, double dV1_,
		     complex<double> shift_sym_phase_,
		     component c_,
		     const void *data_) {
  dft_chunk_data *data = (dft_chunk_data *) data_;
  if (!fc->f[c][0])
    abort("invalid fields_chunk/component combination in dft_chunk");
  
  fc = fc_;
  is = is_;
  ie = ie_;
  s0 = s0_;
  s1 = s1_;
  e0 = e0_;
  e1 = e1_;
  dV0 = dV0_;
  dV1 = dV1_;
  shift_sym_phase = shift_sym_phase_;
  c = c_;

  avg1 = avg2 = 0;
  LOOP_OVER_DIRECTIONS(fc->v.dim,d) {
    if (!fc->v.iyee_shift(c).in_direction(d)) {
      if (avg2) abort("weird yee shift for component %s", component_name(c));
      if (avg1) avg2 = fc->v.stride(d);
      else avg1 = fc->v.stride(d);
    }
  }
  
  omega_min = data->omega_min;
  domega = data->domega;
  Nomega = data->Nomega;
  dft_phase = new complex<double>[Nomega];
  
  LOOP_OVER_IVECS(fc->v, is, ie, idx) {
    (void) loop_is1; (void) loop_is2; (void) loop_is3; // unused
    N = loop_n1 * loop_n2 * loop_n3;
    goto stoploop; // we only need one loop iteration to get params
  }
 stoploop:
  dft = new complex<double>[N * Nomega];
  for (int i = 0; i < N * Nomega; ++i)
    dft[i] = 0.0;
  
  next_in_chunk = fc->dft_chunks;
  fc->dft_chunks = this;
  next_in_dft = data->dft_chunks;
}

dft_chunk::~dft_chunk() {
  delete[] dft;
  delete[] dft_phase;

  // delete from fields_chunk list
  dft_chunk *cur = fc->dft_chunks;
  while (cur && cur->next_in_chunk)
    cur = cur->next_in_chunk;
  if (cur && cur->next_in_chunk == this)
    cur->next_in_chunk = next_in_chunk;
  else if (fc->dft_chunks == this)
    fc->dft_chunks = next_in_chunk;
}

static void add_dft_integrand(fields_chunk *fc,
			      ivec is, ivec ie,
			      vec s0, vec s1, vec e0, vec e1,
			      double dV0, double dV1,
			      vec shift, complex<double> shift_phase,
			      const symmetry &S, int sn,
			      void *integrand_data)
{
  dft_chunk_data *data = (dft_chunk_data *) integrand_data;
  (void) shift; // unused
  
  component c = S.transform(data->c, -sn);
  if (!fc->f[c][0]) return; // this chunk doesn't have component c

  data->dft_chunks = new dft_chunk(fc,is,ie,s0,s1,e0,e1,dV0,dV1,
				   shift_phase * S.phase_shift(c, sn),
				   c, integrand_data);
}

dft_chunk *fields::add_dft(component c, const geometric_volume &where,
			   double freq_min, double freq_max, int Nfreq) {
  dft_chunk_data data;
  
  data.c = c;
  data.omega_min = freq_min * 2*pi;
  data.domega = Nfreq <= 1 ? 0.0 : 
    (freq_max * 2*pi - data.omega_min) / (Nfreq - 1);
  data.Nomega = Nfreq;
  data.dft_chunks = 0;
  
  integrate(add_dft_integrand, (void *) &data, where);

  return data.dft_chunks;
}

void fields::update_dfts() {
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine())
      chunks[i]->update_dfts(time(), time() - 0.5 * inva*c);
}

void fields_chunk::update_dfts(double timeE, double timeH) {
  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_chunk) {
    cur->update_dft(is_magnetic(cur->c) ? timeH : timeE);
  }
}

void dft_chunk::update_dft(double time) {
  for (int i = 0; i < Nomega; ++i)
    dft_phase[i] = polar(1.0, (omega_min + i*domega)*time) * shift_sym_phase;

  // get the integration weight along a given dimension
#define WEIGHT(i, n, dir) ((i > 1 && i < n - 2) ? 1.0 : (i == 0 ? s0.in_direction(direction(dir)) : (i == 1 ? s1.in_direction(direction(dir)) : i == n - 1 ? e0.in_direction(direction(dir)) : (i == n - 2 ? e1.in_direction(direction(dir)) : 1.0))))
  
  if (fc->f[c][1]) // complex fields
    LOOP_OVER_IVECS(fc->v, is, ie, idx) {
      // slightly evil use of loop vars defined within LOOP_OVER_IVECS...
      double w1 = WEIGHT(loop_i1, loop_n1, loop_d1);
      double dV = dV0 + dV1 * loop_i2;
      double w12 = w1 * WEIGHT(loop_i2, loop_n2, loop_d2) * dV;
      double w123 = w12 * WEIGHT(loop_i3, loop_n3, loop_d3);

      (void) loop_is1; (void) loop_is2; (void) loop_is3; // unused
      
      complex<double> f; // field value at epsilon point
      if (avg2 == 0) {
	if (avg1 == 0)
	  f = w123 * complex<double>(fc->f[c][0][idx], fc->f[c][1][idx]);
	else
	  f = (w123 * 0.5) * 
	    (complex<double>(fc->f[c][0][idx], fc->f[c][1][idx]) +
	     complex<double>(fc->f[c][0][idx+avg1], fc->f[c][1][idx+avg1]));
      }
      else { // avg1 != 0, avg2 != 0, avg2 != avg1
	f = (w123 * (1./3.)) * 
	  (complex<double>(fc->f[c][0][idx], fc->f[c][1][idx]) +
	   complex<double>(fc->f[c][0][idx+avg1], fc->f[c][1][idx+avg1]) +
	   complex<double>(fc->f[c][0][idx+avg2], fc->f[c][1][idx+avg2]));
      }
      
      int idx_dft = loop_i3 + loop_n3 * (loop_i2 + loop_n2 * loop_i1);
      for (int i = 0; i < Nomega; ++i)
	dft[Nomega * idx_dft + i] += dft_phase[i] * f;
    }
  else // real fields
    LOOP_OVER_IVECS(fc->v, is, ie, idx) {
      // slightly evil use of loop vars defined within LOOP_OVER_IVECS...
      double w1 = WEIGHT(loop_i1, loop_n1, loop_d1);
      double dV = dV0 + dV1 * loop_i2;
      double w12 = w1 * WEIGHT(loop_i2, loop_n2, loop_d2) * dV;
      double w123 = w12 * WEIGHT(loop_i3, loop_n3, loop_d3);
      
      (void) loop_is1; (void) loop_is2; (void) loop_is3; // unused
      
      double f; // field value at epsilon point
      if (avg2 == 0) {
	if (avg1 == 0)
	  f = w123 * fc->f[c][0][idx];
	else
	  f = (w123 * 0.5) * (fc->f[c][0][idx] + fc->f[c][0][idx+avg1]);
      }
      else { // avg1 != 0, avg2 != 0, avg2 != avg1
	f = (w123 * (1./3.)) * 
	  (fc->f[c][0][idx] + fc->f[c][0][idx+avg1] + fc->f[c][0][idx+avg2]);
      }
      
      int idx_dft = loop_i3 + loop_n3 * (loop_i2 + loop_n2 * loop_i1);
      for (int i = 0; i < Nomega; ++i)
	dft[Nomega * idx_dft + i] += dft_phase[i] * f;
    }

#undef WEIGHT
}

void dft_chunk::negate_dft() {
  for (int i = 0; i < N * Nomega; ++i)
    dft[i] = -dft[i];
}

static int dft_chunks_Ntotal(dft_chunk *dft_chunks, int *nchunks_p) {
  int n = 0, nchunks = 0;
  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_dft, ++nchunks)
    n += cur->N * cur->Nomega * 2;
  if (nchunks_p) *nchunks_p = max_to_all(nchunks);
  return sum_to_all(n);
}

static void dft_filename(char *s, int slen, component c, const char *outdir,
			 bool append_file, const char *prefix) {
  snprintf(s, slen, "%s/" "%s%s" "dft-%s.h5",
           outdir,
           prefix ? prefix : "", prefix && prefix[0] ? "-" : "",
           append_file ? "fields" : component_name(c));
}

void save_dft_hdf5(dft_chunk *dft_chunks, component c, const char *outdir,
		   bool append_file, const char *prefix) {
  int n, nchunks;
  n = dft_chunks_Ntotal(dft_chunks, &nchunks);

  const int buflen = 1024;
  char filename[buflen];
  dft_filename(filename, buflen, c, outdir, append_file, prefix);
  
  int ichunk = 0, istart = 0;
  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_dft, ++ichunk) {
    int Nchunk = cur->N * cur->Nomega * 2;
    h5io::write_chunk(filename, component_name(c),
		      1, &n,
		      (double *) cur->dft,
		      &istart, &Nchunk,
		      true, ichunk == 0,
		      false, -1, 
		      append_file, false);
    istart += Nchunk;
  }
  /* All processes need to call write_chunk in parallel, even if
     some processes have nothing to write. */
  for (; ichunk < nchunks; ++ichunk) {
    int Nchunk = 0;
    h5io::write_chunk(filename, component_name(c),
		      1, &n,
		      (double *) 0,
		      &istart, &Nchunk,
		      true, ichunk == 0,
		      false, -1, 
		      append_file, false);
  }
}

void load_dft_hdf5(dft_chunk *dft_chunks, component c, const char *outdir,
		   bool in_appended_file, const char *prefix) {
  int n, nchunks;
  n = dft_chunks_Ntotal(dft_chunks, &nchunks);

  const int buflen = 1024;
  char filename[buflen];
  dft_filename(filename, buflen, c, outdir, in_appended_file, prefix);
  
  int ichunk = 0, istart = 0;
  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_dft, ++ichunk) {
    int Nchunk = cur->N * cur->Nomega * 2;
    h5io::read_chunk(filename, component_name(c),
		     1, &n,
		     (double *) cur->dft,
		     &istart, &Nchunk,
		     true);
    istart += Nchunk;
  }
  /* All processes need to call read_chunk in parallel, even if
     some processes have nothing to write. */
  for (; ichunk < nchunks; ++ichunk) {
    int Nchunk = 0;
    h5io::read_chunk(filename, component_name(c),
		     1, &n,
		     (double *) 0,
		     &istart, &Nchunk,
		     true);
  }
}

} // namespace meep
