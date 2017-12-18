/* Copyright (C) 2005-2015 Massachusetts Institute of Technology.
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

using namespace std;

typedef complex<double> cdouble;

namespace meep {

struct dft_chunk_data { // for passing to field::loop_in_chunks as void*
  double omega_min, domega;
  int Nomega;
  component c;
  int vc;
  complex<double> weight, extra_weight;
  bool include_dV_and_interp_weights;
  bool sqrt_dV_and_interp_weights;
  dft_chunk *dft_chunks;
};

dft_chunk::dft_chunk(fields_chunk *fc_,
		     ivec is_, ivec ie_,
		     vec s0_, vec s1_, vec e0_, vec e1_,
		     double dV0_, double dV1_,
		     complex<double> extra_weight_,
		     complex<double> scale_,
		     component c_,
		     bool use_centered_grid,
                     ivec shift_, const symmetry &S_, int sn_, int vc_,
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
       multiply by the interpolation weights or the grid_volume twice. */
    dV0 = 1;
    dV1 = 0;
    LOOP_OVER_DIRECTIONS(fc->gv.dim, d) {
      s0.set_direction(d, 1.0);
      s1.set_direction(d, 1.0);
      e0.set_direction(d, 1.0);
      e1.set_direction(d, 1.0);
    }
  }
  /* an alternative way to avoid multipling by interpolation weights twice:
     multiply by square root of the weights */
  include_dV_and_interp_weights = data->include_dV_and_interp_weights;
  sqrt_dV_and_interp_weights = data->sqrt_dV_and_interp_weights;
  scale = scale_ * data->weight;
  extra_weight = extra_weight_;
  c = c_;

  if (use_centered_grid)
    fc->gv.yee2cent_offsets(c, avg1, avg2);
  else
    avg1 = avg2 = 0;

  shift = shift_;
  S = S_; sn = sn_;
  vc = vc_;

  omega_min = data->omega_min;
  domega = data->domega;
  Nomega = data->Nomega;
  dft_phase = new complex<realnum>[Nomega];
  
  N = 1;
  LOOP_OVER_DIRECTIONS(is.dim, d)
    N *= (ie.in_direction(d) - is.in_direction(d)) / 2 + 1;
  dft = new complex<realnum>[N * Nomega];
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
  if (cur == this)
    fc->dft_chunks = next_in_chunk;
  else {
    while (cur && cur->next_in_chunk && cur->next_in_chunk != this)
      cur = cur->next_in_chunk;
    if (cur && cur->next_in_chunk == this)
      cur->next_in_chunk = next_in_chunk;
  }
}

void dft_flux::remove()
{
  while (E) {
    dft_chunk *nxt = E->next_in_dft;
    delete E;
    E = nxt;
  }
  while (H) {
    dft_chunk *nxt = H->next_in_dft;
    delete H;
    H = nxt;
  }
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
  (void) ichunk; // unused

  component c = S.transform(data->c, -sn);
  if (c >= NUM_FIELD_COMPONENTS || !fc->f[c][0])
       return; // this chunk doesn't have component c

  data->dft_chunks = new dft_chunk(fc,is,ie,s0,s1,e0,e1,dV0,dV1,
				   data->extra_weight,
				   shift_phase * S.phase_shift(c, sn),
				   c, cgrid == Centered,
                                   shift, S, sn, data->vc,
				   chunkloop_data);
}

dft_chunk *fields::add_dft(component c, const volume &where,
			   double freq_min, double freq_max, int Nfreq,
			   bool include_dV_and_interp_weights,
			   complex<double> weight, dft_chunk *chunk_next,
			   bool sqrt_dV_and_interp_weights,
			   complex<double> extra_weight,
			   bool use_centered_grid, int vc) {
  if (coordinate_mismatch(gv.dim, c))
    return NULL;

  dft_chunk_data data;  
  data.c = c;
  data.vc = vc;
  if (Nfreq <= 1) freq_min = freq_max = (freq_min + freq_max) * 0.5;
  data.omega_min = freq_min * 2*pi;
  data.domega = Nfreq <= 1 ? 0.0 : 
    (freq_max * 2*pi - data.omega_min) / (Nfreq - 1);
  data.Nomega = Nfreq;
  data.include_dV_and_interp_weights = include_dV_and_interp_weights;
  data.sqrt_dV_and_interp_weights = sqrt_dV_and_interp_weights;
  data.dft_chunks = chunk_next;
  data.weight = weight * (dt/sqrt(2*pi));
  data.extra_weight = extra_weight;
  loop_in_chunks(add_dft_chunkloop, (void *) &data, where,
		 use_centered_grid ? Centered : c);

  return data.dft_chunks;
}

dft_chunk *fields::add_dft(const volume_list *where,
			   double freq_min, double freq_max, int Nfreq,
			   bool include_dV_and_interp_weights) {
  dft_chunk *chunks = 0;
  while (where) {
    if (is_derived(where->c)) abort("derived_component invalid for dft");
    chunks = add_dft(component(where->c), where->v,
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
  LOOP_OVER_IVECS(fc->gv, is, ie, idx) {
    double w = IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2);
    if (sqrt_dV_and_interp_weights) w = sqrt(w);
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
      complex<realnum> fc(f[0], f[1]);
      for (int i = 0; i < Nomega; ++i)
	dft[Nomega * idx_dft + i] += dft_phase[i] * fc;
    }
    else {
      realnum fr = f[0];
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
void save_dft_hdf5(dft_chunk *dft_chunks, const char *name, h5file *file,
		   const char *dprefix) {
  int istart;
  int n = dft_chunks_Ntotal(dft_chunks, &istart);

  char dataname[1024];
  snprintf(dataname, 1024, "%s%s" "%s_dft", 
	   dprefix ? dprefix : "", dprefix && dprefix[0] ? "_" : "", name);
  file->create_data(dataname, 1, &n);

  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_dft) {
    int Nchunk = cur->N * cur->Nomega * 2;
    file->write_chunk(1, &istart, &Nchunk, (realnum *) cur->dft);
    istart += Nchunk;
  }
  file->done_writing_chunks();
}

void save_dft_hdf5(dft_chunk *dft_chunks, component c, h5file *file,
		   const char *dprefix) {
  save_dft_hdf5(dft_chunks, component_name(c), file, dprefix);
}

void load_dft_hdf5(dft_chunk *dft_chunks, const char *name, h5file *file,
		   const char *dprefix) {
  int istart;
  int n = dft_chunks_Ntotal(dft_chunks, &istart);

  char dataname[1024];
  snprintf(dataname, 1024, "%s%s" "%s_dft", 
	   dprefix ? dprefix : "", dprefix && dprefix[0] ? "_" : "", name);
  int file_rank, file_dims;
  file->read_size(dataname, &file_rank, &file_dims, 1);
  if (file_rank != 1 || file_dims != n)
    abort("incorrect dataset size (%d vs. %d) in load_dft_hdf5 %s:%s", file_dims, n, file->file_name(), dataname);
  
  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_dft) {
    int Nchunk = cur->N * cur->Nomega * 2;
    file->read_chunk(1, &istart, &Nchunk, (realnum *) cur->dft);
    istart += Nchunk;
  }
}

void load_dft_hdf5(dft_chunk *dft_chunks, component c, h5file *file,
		   const char *dprefix) {
  load_dft_hdf5(dft_chunks, component_name(c), file, dprefix);
}

dft_flux::dft_flux(const component cE_, const component cH_,
		   dft_chunk *E_, dft_chunk *H_, 
		   double fmin, double fmax, int Nf)
{
  if (Nf <= 1) fmin = fmax = (fmin + fmax) * 0.5;
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

void dft_flux::save_hdf5(fields &f, const char *fname, const char *dprefix,
			 const char *prefix) {
  h5file *ff = f.open_h5file(fname, h5file::WRITE, prefix);
  save_hdf5(ff, dprefix);
  delete ff;
}

void dft_flux::load_hdf5(fields &f, const char *fname, const char *dprefix,
			 const char *prefix) {
  h5file *ff = f.open_h5file(fname, h5file::READONLY, prefix);
  load_hdf5(ff, dprefix);
  delete ff;
}

void dft_flux::scale_dfts(complex<double> scale) {
  if (E) E->scale_dft(scale);
  if (H) H->scale_dft(scale);
}

dft_flux fields::add_dft_flux(const volume_list *where_,
			      double freq_min, double freq_max, int Nfreq) {
  dft_chunk *E = 0, *H = 0;
  component cE[2] = {Ex,Ey}, cH[2] = {Hy,Hx};

  volume_list *where = S.reduce(where_);
  volume_list *where_save = where;
  while (where) {
    derived_component c = derived_component(where->c);
    if (coordinate_mismatch(gv.dim, component_direction(c)))
      abort("coordinate-type mismatch in add_dft_flux");
    
    switch (c) {
    case Sx: cE[0] = Ey, cE[1] = Ez, cH[0] = Hz, cH[1] = Hy; break;
    case Sy: cE[0] = Ez, cE[1] = Ex, cH[0] = Hx, cH[1] = Hz; break;
    case Sr: cE[0] = Ep, cE[1] = Ez, cH[0] = Hz, cH[1] = Hp; break;
    case Sp: cE[0] = Ez, cE[1] = Er, cH[0] = Hr, cH[1] = Hz; break;
    case Sz:
      if (gv.dim == Dcyl)
	cE[0] = Er, cE[1] = Ep, cH[0] = Hp, cH[1] = Hr;
      else
	cE[0] = Ex, cE[1] = Ey, cH[0] = Hy, cH[1] = Hx; 
      break;
    default: abort("invalid flux component!");
    }
    
    for (int i = 0; i < 2; ++i) {
      E = add_dft(cE[i], where->v, freq_min, freq_max, Nfreq,
		  true, where->weight * double(1 - 2*i), E);
      H = add_dft(cH[i], where->v, freq_min, freq_max, Nfreq,
		  false, 1.0, H);
    }
    
    where = where->next;
  }
  delete where_save;

  return dft_flux(cE[0], cH[0], E, H, freq_min, freq_max, Nfreq);
}

direction fields::normal_direction(const volume &where) const {
  direction d = where.normal_direction();
  if (d == NO_DIRECTION) {
    /* hack so that we still infer the normal direction correctly for
       volumes with empty dimensions */
    volume where_pad(where);
    LOOP_OVER_DIRECTIONS(where.dim, d1)
      if (nosize_direction(d1) && where.in_direction(d1) == 0.0)
	where_pad.set_direction_max(d1, where.in_direction_min(d1) + 0.1);
    d = where_pad.normal_direction();  
    if (d == NO_DIRECTION)
      abort("Could not determine normal direction for given grid_volume.");
  }
  return d;
}

dft_flux fields::add_dft_flux(direction d, const volume &where,
			      double freq_min, double freq_max, int Nfreq) {
  if (d == NO_DIRECTION)
    d = normal_direction(where);
  volume_list vl(where, direction_component(Sx, d));
  return add_dft_flux(&vl, freq_min, freq_max, Nfreq);
}

dft_flux fields::add_dft_flux_box(const volume &where,
				  double freq_min, double freq_max, int Nfreq){
  volume_list *faces = 0;
  LOOP_OVER_DIRECTIONS(where.dim, d)
    if (where.in_direction(d) > 0) {
      volume face(where);
      derived_component c = direction_component(Sx, d);
      face.set_direction_min(d, where.in_direction_max(d));
      faces = new volume_list(face, c, +1, faces);
      face.set_direction_min(d, where.in_direction_min(d));
      face.set_direction_max(d, where.in_direction_min(d));
      faces = new volume_list(face, c, -1, faces);
    }

  dft_flux flux = add_dft_flux(faces, freq_min, freq_max, Nfreq);
  delete faces;
  return flux;
}

dft_flux fields::add_dft_flux_plane(const volume &where,
			      double freq_min, double freq_max, int Nfreq) {
  return add_dft_flux(NO_DIRECTION, where, freq_min, freq_max, Nfreq);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef enum { OUTPUT_FLUX, OUTPUT_MODE, MODE_FLUX, MODE_MODE, NO_OP } flux_operation;

cdouble dft_chunk::do_flux_operation(int rank,
                                     direction *ds,
                                     ivec min_corner,
                                     h5file *file,
                                     double *buffer,
                                     int reim,
                                     void *mode1_data,
                                     component mode1_c,
                                     void *mode2_data,
                                     component mode2_c,
                                     int num_freq,
                                     double flux_sign)
{
   /*****************************************************************/
   /* compute the size of the chunk we own and its strides etc.     */
   /*****************************************************************/
   int start[3]={0,0,0}, count[3]={1,1,1};
   int offset[3]={0,0,0}, stride[3]={1,1,1};
   ivec isS = S.transform(is, sn) + shift;
   ivec ieS = S.transform(ie, sn) + shift;

   ivec permute(zero_ivec(fc->gv.dim));
   for (int i = 0; i < 3; ++i) 
    permute.set_direction(fc->gv.yucky_direction(i), i);
   permute = S.transform_unshifted(permute, sn);
   LOOP_OVER_DIRECTIONS(permute.dim, d)
    permute.set_direction(d, abs(permute.in_direction(d)));

   for (int i = 0; i < rank; ++i)
    { direction d = ds[i];
      int isd = isS.in_direction(d), ied = ieS.in_direction(d);
      start[i] = (min(isd, ied) - min_corner.in_direction(d)) / 2;
      count[i] = abs(ied - isd) / 2 + 1;
      if (ied < isd) offset[permute.in_direction(d)] = count[i] - 1;
    };

   for (int i = 0; i < rank; ++i)
    { direction d = ds[i];
      int j = permute.in_direction(d);
      for (int k = i + 1; k < rank; ++k)
      stride[j] *= count[k];
      offset[j] *= stride[j];
      if (offset[j]) stride[j] *= -1;
    };

   /***************************************************************/
   /* loop over all grid points in our piece of the volume        */
   /***************************************************************/
   vec rshift(shift * (0.5*fc->gv.inva));
   int chunk_idx = 0;
   cdouble integral=0.0;
   LOOP_OVER_IVECS(fc->gv, is, ie, idx)
    {
      IVEC_LOOP_LOC(fc->gv, loc);
      loc = S.transform(loc, sn) + rshift;
      double w = IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2);
      cdouble fluxval = flux_sign*dft[ Nomega*(chunk_idx++) + num_freq];
      if (include_dV_and_interp_weights)
       fluxval /= (sqrt_dV_and_interp_weights ? sqrt(w) : w);

      cdouble mode1val=0.0, mode2val=0.0;
      if (mode1_data)
       mode1val=eigenmode_amplitude(loc,mode1_data,mode1_c);
      if (mode2_data)
       mode2val=eigenmode_amplitude(loc,mode2_data,mode2_c);

      if (file)
       { int idx2 = ((((offset[0] + offset[1] + offset[2])
                                  + loop_i1 * stride[0])
                                  + loop_i2 * stride[1])
                                  + loop_i3 * stride[2]);
         cdouble val = (mode1_data ? mode1val : fluxval);
         buffer[idx2] = reim ? imag(val) : real(val);
       }
      else
       { if (mode2_data)
          integral += w*conj(mode1val)*mode2val;
         else
          integral += w*conj(mode1val)*fluxval;
       };

    }; // LOOP_OVER_IVECS(fc->gv, is, ie, idx)

  if (file)
   file->write_chunk(rank, start, count, buffer);

  return integral;

}

/***************************************************************/
/* flux_operation is an omnibus routine that serves as the     */
/* computational back end for the following routines:          */
/*  output_flux_fields()                                       */
/*  output_mode_fields()                                       */
/*  get_mode_flux_overlap()                                    */
/*  get_mode_mode_overlap()                                    */
/*                                                             */
/* This routine does one or two things depending on the input. */
/*                                                             */
/* (A) HDF5 file output (if HDF5FileName is non-null)          */
/*                                                             */
/*     (A1) If mode1_data is NULL, write all field components  */
/*          stored in flux (at all frequencies) to HDF5 file.  */
/*                                                             */
/*     (A2) If mode1_data is non_null, write all field         */
/*          components of the eigenmode field described by     */
/*          mode1_data to HDF5 file.                           */
/*                                                             */
/* (B) Computation of overlap integrals.                       */
/*                                                             */
/*     (B1) If mode1_data is non-NULL and mode2_data is NULL,  */
/*          compute and return the overlap integral between    */
/*          the fields described by flux (at the #num_freqth   */
/*          of the frequencies for which flux contains data)   */
/*          and the eigenmode field described by mode1_data.   */
/*                                                             */
/*     (B2) If flux is NULL and mode1_data, mode2_data are     */ 
/*          both non-NULL, compute and return the overlap      */ 
/*          integral between the eigenmode fields described    */ 
/*          by mode1_data and mode2_data.                      */
/***************************************************************/
cdouble fields::do_flux_operation(dft_flux *flux, const volume where,
                                  const char *HDF5FileName, void *mode1_data, void *mode2_data, int num_freq)
{ 
  /***************************************************************/
  /* look at input arguments to figure out what to do  ***********/
  /***************************************************************/
  flux_operation flux_ops[4];
  int num_ops=0;
  if (HDF5FileName!=0 && mode1_data==0)
   flux_ops[num_ops++] = OUTPUT_FLUX;
  if (HDF5FileName!=0 && mode1_data!=0)
   flux_ops[num_ops++] = OUTPUT_MODE;
  if (HDF5FileName==0 && mode1_data!=0 && mode2_data==0)
   flux_ops[num_ops++] = MODE_FLUX;
  if (HDF5FileName==0 && mode1_data!=0 && mode2_data!=0)
   flux_ops[num_ops++] = MODE_MODE;
  if (num_ops==0)
   abort("no operation specified for do_flux_operation");
  if (num_ops>1)
   abort("more than one operation specified for do_flux_operation");

  flux_operation flux_op=flux_ops[0];
  
  /***************************************************************/
  /* get statistics on the volume slice **************************/
  /***************************************************************/
  int bufsz=0;
  ivec min_corner = gv.round_vec(where.get_max_corner()) + one_ivec(gv.dim);
  ivec max_corner = gv.round_vec(where.get_min_corner()) - one_ivec(gv.dim);
  for (dft_chunk *E=flux->E; E; E=E->next_in_dft)
   {
     ivec isS = E->S.transform(E->is, E->sn) + E->shift;
     ivec ieS = E->S.transform(E->ie, E->sn) + E->shift;
     min_corner = min(min_corner, min(isS, ieS));
     max_corner = max(max_corner, max(isS, ieS));
     int this_bufsz=1;
     LOOP_OVER_DIRECTIONS(E->fc->gv.dim, d)
      this_bufsz *= (E->ie.in_direction(d) - E->is.in_direction(d)) / 2 + 1;
     bufsz = max(bufsz, this_bufsz);
   };
  max_corner = max_to_all(max_corner);
  min_corner = -max_to_all(-min_corner); // i.e., min_to_all

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int rank = 0, dims[3];
  direction ds[3];
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    if (rank >= 3) abort("too many dimensions in output_hdf5_flux");
    int n = (max_corner.in_direction(d)
	     - min_corner.in_direction(d)) / 2 + 1;
    if (n > 1) {
      ds[rank] = d;
      dims[rank++] = n;
    }
  };

  /***************************************************************/
  /* buffer for process-local contributions to HDF5 output files,*/
  /* like h5_output_data::buf in h5fields.cpp                    */
  /***************************************************************/
  realnum *buffer = 0;
  if (HDF5FileName)
   buffer = new realnum[bufsz];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  component cE[2]={Ex, Ey}, cH[2]={Hy, Hx};
  switch (normal_direction(where))
   { case X: cE[0] = Ey, cE[1] = Ez, cH[0] = Hz, cH[1] = Hy; break;
     case Y: cE[0] = Ez, cE[1] = Ex, cH[0] = Hx, cH[1] = Hz; break;
     case R: cE[0] = Ep, cE[1] = Ez, cH[0] = Hz, cH[1] = Hp; break;
     case P: cE[0] = Ez, cE[1] = Er, cH[0] = Hr, cH[1] = Hz; break;
     case Z:
      if (gv.dim == Dcyl)
	cE[0] = Er, cE[1] = Ep, cH[0] = Hp, cH[1] = Hr;
      else
	cE[0] = Ex, cE[1] = Ey, cH[0] = Hy, cH[1] = Hx; 
     break;
     default: abort("invalid flux component!");
   };

  component all_components[6] = {Ex, Ey, Ez, Hx, Hy, Hz};

  // tangential components and their 'conjugates'
  component tang_components[4], conj_components[4];
  tang_components[0] = cE[0]; conj_components[0] = cH[0];
  tang_components[1] = cE[1]; conj_components[1] = cH[1];
  tang_components[2] = cH[0]; conj_components[2] = cE[0];
  tang_components[3] = cH[1]; conj_components[3] = cE[1];

  int ncOverlap=2;
  if (flux_op==MODE_FLUX || flux_op==MODE_MODE)
   { // default ExHy - EyHx;
     // 0,1,2,3 ExHy, EyHx, HyEx, HxEy
     // 4       ExHy - EyHx + HyEx - HxEy
     // 5       HyEx - HxEy
     char *s=getenv("MEEP_OVERLAP_ALGORITHM");
     if (s && (s[0]>='0' && s[0]<='3') )
      { ncOverlap=1;
        int ic=s[0] - '0';
        tang_components[0] = tang_components[ic];
        conj_components[0] = conj_components[ic];
      }
     else if (s && s[0]=='4')
      ncOverlap=4;
     else if (s && s[0]=='5')
      { tang_components[0]=tang_components[2]; conj_components[0]=conj_components[2];
        tang_components[1]=tang_components[3]; conj_components[1]=conj_components[3];
      };
   };

  /***************************************************************/
  /* set up limits for loops over frequencies, components, etc.  */
  /*                                                             */
  /* which of the frequencies in the flux object will we use?    */
  /*  -- if writing flux fields to HDF5:   all frequencies       */
  /*  -- if computing <mode|flux> overlap: only freq #num_freq   */
  /*  -- otherwise:                        flux fields not used  */
  /*                                                             */
  /* which field components will we consider?                    */
  /*  -- if writing flux fields to HDF5: tangential cmpts only   */
  /*  -- if writing mode fields to HDF5: all components          */
  /*  -- if computing an overlap integral:                       */
  /*                                                             */
  /* loop over real/imaginary parts of field components? yes for */
  /*  -- if writing flux or mode fields to HDF5: yes             */
  /*  -- if computing overlap integral: no                       */
  /***************************************************************/
  int nf_min, nf_max, num_components, reim_max;
  switch(flux_op)
   { case OUTPUT_FLUX:
       nf_min=0; nf_max=flux->Nfreq-1; num_components=4; reim_max=1;
       break;
     case OUTPUT_MODE:
       nf_min=0; nf_max=0;             num_components=6; reim_max=1;
       break;
     case MODE_FLUX:  
       nf_min=nf_max=num_freq;         num_components=ncOverlap; reim_max=0;
       break;
     case MODE_MODE:  
       nf_min=nf_max=0;                num_components=ncOverlap; reim_max=0;
       break;
     case NO_OP:
       abort("%s:%i: internal error",__FILE__,__LINE__);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  bool append_data      = false;
  bool single_precision = false;
  bool First            = true;
  cdouble integral      = 0.0;
  for(int nf=nf_min; nf<=nf_max; nf++)
   for(int nc=0; nc<num_components; nc++)
    for(int reim=0; reim<=reim_max; reim++)
     { 
       char dataname[100];
       component c_flux=Ex, c_mode=Ex, c_mode2=Ex;
       if (flux_op==OUTPUT_FLUX)
        { 
          c_flux=tang_components[nc];
          snprintf(dataname,100,"%s_%i.%c",component_name(c_flux),nf,reim ? 'i' : 'r');
        }
       else if (flux_op==OUTPUT_MODE)
        { c_flux = tang_components[0];
          c_mode = all_components[nc];
          snprintf(dataname,100,"%s.%c",component_name(c_mode),reim ? 'i' : 'r');
        }
       else if (flux_op==MODE_FLUX)
        { c_mode=tang_components[nc];
          c_flux=conj_components[nc];
        }
       else if (flux_op==MODE_MODE)
        { c_mode =tang_components[nc];
          c_mode2=conj_components[nc];
        };

       h5file *file=0;
       if (HDF5FileName)
        { file = open_h5file(HDF5FileName, First ? h5file::WRITE : h5file::READWRITE, 0, false);
          First=false;
          file->create_or_extend_data(dataname, rank, dims, append_data, single_precision);
        };

       // the second component of the E field is stored with a
       // minus sign in the flux object, which we need to remove
       double flux_sign = 1.0;
       if ( (flux_op==OUTPUT_FLUX || flux_op==MODE_FLUX) && (c_flux==cE[1]) )
        flux_sign = -1.0;

       // 
       double integral_sign = ( ((nc%2)==1) ? -1.0 : 1.0);

       for (dft_chunk *EH = (c_flux>=Hx ? flux->H : flux->E); EH; EH=EH->next_in_dft)
        if (EH->c == c_flux)
         integral += integral_sign * EH->do_flux_operation(rank, ds, min_corner, file, buffer, reim,
                                                           mode1_data, c_mode, mode2_data, c_mode2, nf, flux_sign);

       if (file)
        { file->done_writing_chunks();
          file->prevent_deadlock(); // hackery
          delete file;
        };
     };

  if (buffer)
   delete[] buffer;

  return sum_to_all(integral);
  /***************************************************************/
  /* write the lower-left grid corner and grid spacing           */
  /* to the hdf5 file so that it contains enough information     */
  /* to recreate the coordinates of the grid points              */
  /***************************************************************/
#if 0
  if (am_master())
   { bool parallel=false, single_precision=false;
     char filename[100];
     snprintf(filename,100,"%s.h5",HDF5FileName);
     h5file file(filename, h5file::READWRITE, parallel);
     double xmin[3];
     for(int nd=0; nd<rank; nd++)
      xmin[nd] = where.in_direction_min( ds[nd] );

     dims[0]=rank;
     file.write("min_corner",1,dims,xmin,single_precision);

     dims[0]=1;
     file.write("inva",1,dims,&(gv.inva),single_precision);
   };
#endif

}

/***************************************************************/
/* entry points to flux_operation for various specific         */
/* calculations                                                */
/***************************************************************/
void fields::output_flux_fields(dft_flux *flux, const volume where, const char *HDF5FileName)
{ do_flux_operation(flux, where, HDF5FileName); }

void fields::output_mode_fields(void *mode_data, dft_flux *flux, const volume where, const char *HDF5FileName)
{ do_flux_operation(flux, where, HDF5FileName, mode_data); }
 
cdouble fields::get_mode_flux_overlap(void *mode_data, dft_flux *flux, int num_freq, const volume where)
{ return do_flux_operation(flux, where, 0, mode_data, 0, num_freq); }

cdouble fields::get_mode_mode_overlap(void *mode1_data, void *mode2_data, dft_flux *flux, const volume where)
{ return do_flux_operation(flux, where, 0, mode1_data, mode2_data); }

} // namespace meep
