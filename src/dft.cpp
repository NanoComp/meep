/* Copyright (C) 2006 Massachusetts Institute of Technology.
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

#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {

struct dft_chunk_data { // for passing to field::loop_in_chunks as void*
  double omega_min, domega;
  int Nomega;
  component c;
  complex<double> weight;
  bool include_dV_and_interp_weights;
  dft_chunk *dft_chunks;
};

dft_chunk::dft_chunk(fields_chunk *fc_,
		     ivec is_, ivec ie_,
		     vec s0_, vec s1_, vec e0_, vec e1_,
		     double dV0_, double dV1_,
		     complex<double> scale_,
		     component c_,
		     const void *data_) {
  dft_chunk_data *data = (dft_chunk_data *) data_;
  if (!fc_->f[c_][0])
    abort("invalid fields_chunk/component combination in dft_chunk");
  
  fc = fc_;
  is = is_;
  ie = ie_;
  s0 = s0_;
  s1 = s1_;
  e0 = e0_;
  e1 = e1_;
  if (data->include_dV_and_interp_weights) {
    dV0 = dV0_;
    dV1 = dV1_;
  }
  else {
    /* this is for e.g. computing E x H, where we don't want to 
       multiply by the interpolation weights or the volume twice. */
    dV0 = 1;
    dV1 = 0;
    LOOP_OVER_DIRECTIONS(fc->v.dim, d) {
      s0.set_direction(d, 1.0);
      s1.set_direction(d, 1.0);
      e0.set_direction(d, 1.0);
      e1.set_direction(d, 1.0);
    }
  }
  scale = scale_ * data->weight;
  c = c_;

  fc->v.yee2diel_offsets(c, avg1, avg2);

  omega_min = data->omega_min;
  domega = data->domega;
  Nomega = data->Nomega;
  dft_phase = new complex<double>[Nomega];
  
  N = 1;
  LOOP_OVER_DIRECTIONS(is.dim, d)
    N *= (ie.in_direction(d) - is.in_direction(d)) / 2 + 1;
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
  while (cur && cur->next_in_chunk && cur->next_in_chunk != this)
    cur = cur->next_in_chunk;
  if (cur && cur->next_in_chunk == this)
    cur->next_in_chunk = next_in_chunk;
  else if (fc->dft_chunks == this)
    fc->dft_chunks = next_in_chunk;

  delete next_in_dft;
}

static void add_dft_chunkloop(fields_chunk *fc, int ichunk, component cgrid,
			      ivec is, ivec ie,
			      vec s0, vec s1, vec e0, vec e1,
			      double dV0, double dV1,
			      ivec shift, complex<double> shift_phase,
			      const symmetry &S, int sn,
			      void *chunkloop_data)
{
  dft_chunk_data *data = (dft_chunk_data *) chunkloop_data;
  (void) shift; // unused

  if (cgrid != Dielectric) abort("dft chunks should use the Dielectric grid");
  
  component c = S.transform(data->c, -sn);
  if (c >= NUM_FIELD_COMPONENTS || !fc->f[c][0])
       return; // this chunk doesn't have component c

  data->dft_chunks = new dft_chunk(fc,is,ie,s0,s1,e0,e1,dV0,dV1,
				   shift_phase * S.phase_shift(c, sn),
				   c, chunkloop_data);
}

dft_chunk *fields::add_dft(component c, const geometric_volume &where,
			   double freq_min, double freq_max, int Nfreq,
			   bool include_dV_and_interp_weights,
			   complex<double> weight, dft_chunk *chunk_next) {
  if (coordinate_mismatch(v.dim, c))
    return NULL;

  dft_chunk_data data;  
  data.c = c;
  data.omega_min = freq_min * 2*pi;
  data.domega = Nfreq <= 1 ? 0.0 : 
    (freq_max * 2*pi - data.omega_min) / (Nfreq - 1);
  data.Nomega = Nfreq;
  data.include_dV_and_interp_weights = include_dV_and_interp_weights;
  data.dft_chunks = chunk_next;
  data.weight = weight * (dt/sqrt(2*pi));
  loop_in_chunks(add_dft_chunkloop, (void *) &data, where);

  return data.dft_chunks;
}

dft_chunk *fields::add_dft(const geometric_volume_list *where,
			   double freq_min, double freq_max, int Nfreq,
			   bool include_dV_and_interp_weights) {
  dft_chunk *chunks = 0;
  while (where) {
    if (is_derived(where->c)) abort("derived_component invalid for dft");
    chunks = add_dft(component(where->c), where->gv,
		     freq_min, freq_max, Nfreq, include_dV_and_interp_weights,
		     where->weight, chunks);
    where = where->next;
  }
  return chunks;
}

dft_chunk *fields::add_dft_pt(component c, const vec &where,
			   double freq_min, double freq_max, int Nfreq) {
  return add_dft(c, where, freq_min, freq_max, Nfreq, false);
}

void fields::update_dfts() {
  am_now_working_on(FourierTransforming);
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine())
      chunks[i]->update_dfts(time(), time() - 0.5 * dt);
  finished_working();
}

void fields_chunk::update_dfts(double timeE, double timeH) {
  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_chunk) {
    cur->update_dft(is_magnetic(cur->c) ? timeH : timeE);
  }
}

void dft_chunk::update_dft(double time) {
  if (!fc->f[c][0]) return;

  for (int i = 0; i < Nomega; ++i)
    dft_phase[i] = polar(1.0, (omega_min + i*domega)*time) * scale;

  int numcmp = fc->f[c][1] ? 2 : 1;

  int idx_dft = 0;
  LOOP_OVER_IVECS(fc->v, is, ie, idx) {
    double w = IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2);
    double f[2]; // real/imag field value at epsilon point
    if (avg2)
      for (int cmp=0; cmp < numcmp; ++cmp)
	f[cmp] = (w * 0.25) * 
	  (fc->f[c][cmp][idx] + fc->f[c][cmp][idx+avg1]
	   + fc->f[c][cmp][idx+avg2] + fc->f[c][cmp][idx+(avg1+avg2)]);
    else if (avg1)
      for (int cmp=0; cmp < numcmp; ++cmp)
	f[cmp] = (w * 0.5) * (fc->f[c][cmp][idx] + fc->f[c][cmp][idx+avg1]);
    else
      for (int cmp=0; cmp < numcmp; ++cmp)
	f[cmp] = w * fc->f[c][cmp][idx];
    
    if (numcmp == 2) {
      complex<double> fc(f[0], f[1]);
      for (int i = 0; i < Nomega; ++i)
	dft[Nomega * idx_dft + i] += dft_phase[i] * fc;
    }
    else {
      double fr = f[0];
      for (int i = 0; i < Nomega; ++i)
	dft[Nomega * idx_dft + i] += dft_phase[i] * fr;
    }
    idx_dft++;
  }
}

void dft_chunk::scale_dft(complex<double> scale) {
  for (int i = 0; i < N * Nomega; ++i)
    dft[i] *= scale;
  if (next_in_dft)
    next_in_dft->scale_dft(scale);
}

void dft_chunk::operator-=(const dft_chunk &chunk) {
  if (c != chunk.c || N * Nomega != chunk.N * chunk.Nomega) abort("Mismatched chunks in dft_chunk::operator-=");

  for (int i = 0; i < N * Nomega; ++i)
    dft[i] -= chunk.dft[i];

  if (next_in_dft) {
    if (!chunk.next_in_dft) abort("Mismatched chunk lists in dft_chunk::operator-=");
    *next_in_dft -= *chunk.next_in_dft;
  }
}

static int dft_chunks_Ntotal(dft_chunk *dft_chunks, int *my_start) {
  int n = 0;
  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_dft)
    n += cur->N * cur->Nomega * 2;
  *my_start = partial_sum_to_all(n) - n; // sum(n) for processes before this
  return sum_to_all(n);
}

// Note: the file must have been created in parallel mode, typically via fields::open_h5file.
void save_dft_hdf5(dft_chunk *dft_chunks, component c, h5file *file,
		   const char *dprefix) {
  int istart;
  int n = dft_chunks_Ntotal(dft_chunks, &istart);

  char dataname[1024];
  snprintf(dataname, 1024, "%s%s" "%s_dft", 
	   dprefix ? dprefix : "", dprefix ? "_" : "", component_name(c));
  file->create_data(dataname, 1, &n);

  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_dft) {
    int Nchunk = cur->N * cur->Nomega * 2;
    file->write_chunk(1, &istart, &Nchunk, (double *) cur->dft);
    istart += Nchunk;
  }
  file->done_writing_chunks();
}

void load_dft_hdf5(dft_chunk *dft_chunks, component c, h5file *file,
		   const char *dprefix) {
  int istart;
  int n = dft_chunks_Ntotal(dft_chunks, &istart);

  char dataname[1024];
  snprintf(dataname, 1024, "%s%s" "%s_dft", 
	   dprefix ? dprefix : "", dprefix ? "_" : "", component_name(c));
  int file_rank, file_dims;
  file->read_size(dataname, &file_rank, &file_dims, 1);
  if (file_rank != 1 || file_dims != n)
    abort("incorrect dataset size (%d vs. %d) in load_dft_hdf5 %s:%s", file_dims, n, file->file_name(), dataname);
  
  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_dft) {
    int Nchunk = cur->N * cur->Nomega * 2;
    file->read_chunk(1, &istart, &Nchunk, (double *) cur->dft);
    istart += Nchunk;
  }
}

dft_flux::dft_flux(const component cE_, const component cH_,
		   dft_chunk *E_, dft_chunk *H_, 
		   double fmin, double fmax, int Nf)
{
  freq_min = fmin;
  Nfreq = Nf;
  dfreq = Nf <= 1 ? 0.0 : (fmax - fmin) / (Nf - 1);
  E = E_; H = H_;
  cE = cE_; cH = cH_;
}

dft_flux::dft_flux(const dft_flux &f) {
  freq_min = f.freq_min; Nfreq = f.Nfreq; dfreq = f.dfreq;
  E = f.E; H = f.H;
  cE = f.cE; cH = f.cH;
}

double *dft_flux::flux() {
  double *F = new double[Nfreq];
  for (int i = 0; i < Nfreq; ++i) F[i] = 0;
  for (dft_chunk *curE = E, *curH = H; curE && curH;
       curE = curE->next_in_dft, curH = curH->next_in_dft)
    for (int k = 0; k < curE->N; ++k)
      for (int i = 0; i < Nfreq; ++i)
	F[i] += real(curE->dft[k*Nfreq + i]
		     * conj(curH->dft[k*Nfreq + i]));
  double *Fsum = new double[Nfreq];
  sum_to_all(F, Fsum, Nfreq);
  delete[] F;
  return Fsum;
}

void dft_flux::save_hdf5(h5file *file, const char *dprefix) {
  save_dft_hdf5(E, cE, file, dprefix);
  file->prevent_deadlock(); // hackery
  save_dft_hdf5(H, cH, file, dprefix);
}

void dft_flux::load_hdf5(h5file *file, const char *dprefix) {
  load_dft_hdf5(E, cE, file, dprefix);
  file->prevent_deadlock(); // hackery
  load_dft_hdf5(H, cH, file, dprefix);
}

void dft_flux::save_hdf5(fields &f, const char *fname, const char *dprefix) {
  h5file *ff = f.open_h5file(fname);
  save_hdf5(ff, dprefix);
  delete ff;
}

void dft_flux::load_hdf5(fields &f, const char *fname, const char *dprefix) {
  h5file *ff = f.open_h5file(fname, h5file::READONLY);
  load_hdf5(ff, dprefix);
  delete ff;
}

void dft_flux::scale_dfts(complex<double> scale) {
  if (E) E->scale_dft(scale);
  if (H) H->scale_dft(scale);
}

dft_flux fields::add_dft_flux(const geometric_volume_list *where_,
			      double freq_min, double freq_max, int Nfreq) {
  dft_chunk *E = 0, *H = 0;
  component cE[2] = {Ex,Ey}, cH[2] = {Hy,Hx};

  geometric_volume_list *where = S.reduce(where_);
  while (where) {
    derived_component c = derived_component(where->c);
    if (coordinate_mismatch(v.dim, component_direction(c)))
      abort("coordinate-type mismatch in add_dft_flux");
    
    switch (c) {
    case Sx: cE[0] = Ey, cE[1] = Ez, cH[0] = Hz, cH[1] = Hy; break;
    case Sy: cE[0] = Ez, cE[1] = Ex, cH[0] = Hx, cH[1] = Hz; break;
    case Sr: cE[0] = Ep, cE[1] = Ez, cH[0] = Hz, cH[1] = Hp; break;
    case Sp: cE[0] = Ez, cE[1] = Er, cH[0] = Hr, cH[1] = Hz; break;
    case Sz:
      if (v.dim == Dcyl)
	cE[0] = Er, cE[1] = Ep, cH[0] = Hp, cH[1] = Hr;
      else
	cE[0] = Ex, cE[1] = Ey, cH[0] = Hy, cH[1] = Hx; 
      break;
    default: abort("invalid flux component!");
    }
    
    for (int i = 0; i < 2; ++i) {
      E = add_dft(cE[i], where->gv, freq_min, freq_max, Nfreq,
		  true, where->weight * double(1 - 2*i), E);
      H = add_dft(cH[i], where->gv, freq_min, freq_max, Nfreq,
		  false, 1.0, H);
    }
    
    where = where->next;
  }

  return dft_flux(cE[0], cH[0], E, H, freq_min, freq_max, Nfreq);
}

direction fields::normal_direction(const geometric_volume &where) const {
  direction d = where.normal_direction();
  if (d == NO_DIRECTION) {
    /* hack so that we still infer the normal direction correctly for
       volumes with empty dimensions */
    geometric_volume where_pad(where);
    LOOP_OVER_DIRECTIONS(where.dim, d1)
      if (nosize_direction(d1) && where.in_direction(d1) == 0.0)
	where_pad.set_direction_max(d1, where.in_direction_min(d1) + 0.1);
    d = where_pad.normal_direction();  
    if (d == NO_DIRECTION)
      abort("Could not determine normal direction for given volume.");
  }
  return d;
}

dft_flux fields::add_dft_flux(direction d, const geometric_volume &where,
			      double freq_min, double freq_max, int Nfreq) {
  if (d == NO_DIRECTION)
    d = normal_direction(where);
  geometric_volume_list gvl(where, direction_component(Sx, d));
  return add_dft_flux(&gvl, freq_min, freq_max, Nfreq);
}

dft_flux fields::add_dft_flux_box(const geometric_volume &where,
				  double freq_min, double freq_max, int Nfreq){
  geometric_volume_list *faces = 0;
  int ifaces = 0;
  LOOP_OVER_DIRECTIONS(where.dim, d)
    if (where.in_direction(d) > 0) {
      geometric_volume face(where);
      derived_component c = direction_component(Sx, d);
      face.set_direction_min(d, where.in_direction_max(d));
      faces = new geometric_volume_list(face, c, +1, faces);
      face.set_direction_min(d, where.in_direction_min(d));
      face.set_direction_max(d, where.in_direction_min(d));
      faces = new geometric_volume_list(face, c, -1, faces);
    }

  dft_flux flux = add_dft_flux(faces, freq_min, freq_max, Nfreq);
  delete faces;
  return flux;
}

dft_flux fields::add_dft_flux_plane(const geometric_volume &where,
			      double freq_min, double freq_max, int Nfreq) {
  return add_dft_flux(NO_DIRECTION, where, freq_min, freq_max, Nfreq);
}

} // namespace meep
