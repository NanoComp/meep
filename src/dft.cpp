/* Copyright (C) 2005-2019 Massachusetts Institute of Technology.
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
#include <string.h>
#include <algorithm>

#include "meep.hpp"
#include "meep_internals.hpp"

using namespace std;

typedef complex<double> cdouble;

namespace meep {

struct dft_chunk_data { // for passing to field::loop_in_chunks as void*
  component c;
  int vc;
  double omega_min, domega;
  int Nomega;
  complex<double> stored_weight, extra_weight;
  double dt_factor;
  bool include_dV_and_interp_weights;
  bool sqrt_dV_and_interp_weights;
  bool empty_dim[5];
  dft_chunk *dft_chunks;
};

dft_chunk::dft_chunk(fields_chunk *fc_, ivec is_, ivec ie_, vec s0_, vec s1_, vec e0_, vec e1_,
                     double dV0_, double dV1_, component c_, bool use_centered_grid,
                     cdouble phase_factor, ivec shift_, const symmetry &S_, int sn_,
                     const void *data_) {
  dft_chunk_data *data = (dft_chunk_data *)data_;
  if (!fc_->f[c_][0]) abort("invalid fields_chunk/component combination in dft_chunk");

  fc = fc_;
  is = is_;
  ie = ie_;
  s0 = s0_;
  s1 = s1_;
  e0 = e0_;
  e1 = e1_;
  dV0 = dV0_;
  dV1 = dV1_;

  c = c_;

  if (use_centered_grid)
    fc->gv.yee2cent_offsets(c, avg1, avg2);
  else
    avg1 = avg2 = 0;

  stored_weight = data->stored_weight;
  extra_weight = data->extra_weight;
  scale = stored_weight * phase_factor * data->dt_factor;

  /* this is for e.g. computing E x H, where we don't want to
     multiply by the interpolation weights or the grid_volume twice. */
  include_dV_and_interp_weights = data->include_dV_and_interp_weights;

  /* an alternative way to avoid multipling by interpolation weights twice:
     multiply by square root of the weights */
  sqrt_dV_and_interp_weights = data->sqrt_dV_and_interp_weights;

  shift = shift_;
  S = S_;
  sn = sn_;
  vc = data->vc;

  omega_min = data->omega_min;
  domega = data->domega;
  Nomega = data->Nomega;
  dft_phase = new complex<realnum>[Nomega];

  N = 1;
  LOOP_OVER_DIRECTIONS(is.dim, d) { N *= (ie.in_direction(d) - is.in_direction(d)) / 2 + 1; }
  dft = new complex<realnum>[N * Nomega];
  for (size_t i = 0; i < N * Nomega; ++i)
    dft[i] = 0.0;
  for (int i = 0; i < 5; ++i)
    empty_dim[i] = data->empty_dim[i];

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
    if (cur && cur->next_in_chunk == this) cur->next_in_chunk = next_in_chunk;
  }
}

void dft_flux::remove() {
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

static void add_dft_chunkloop(fields_chunk *fc, int ichunk, component cgrid, ivec is, ivec ie,
                              vec s0, vec s1, vec e0, vec e1, double dV0, double dV1, ivec shift,
                              complex<double> shift_phase, const symmetry &S, int sn,
                              void *chunkloop_data) {
  dft_chunk_data *data = (dft_chunk_data *)chunkloop_data;
  (void)ichunk; // unused

  component c = S.transform(data->c, -sn);
  if (c >= NUM_FIELD_COMPONENTS || !fc->f[c][0]) return; // this chunk doesn't have component c

  data->dft_chunks =
      new dft_chunk(fc, is, ie, s0, s1, e0, e1, dV0, dV1, c, cgrid == Centered,
                    shift_phase * S.phase_shift(c, sn), shift, S, sn, chunkloop_data);
}

dft_chunk *fields::add_dft(component c, const volume &where, double freq_min, double freq_max,
                           int Nfreq, bool include_dV_and_interp_weights,
                           complex<double> stored_weight, dft_chunk *chunk_next,
                           bool sqrt_dV_and_interp_weights, complex<double> extra_weight,
                           bool use_centered_grid, int vc) {
  if (coordinate_mismatch(gv.dim, c)) return NULL;

  /* If you call add_dft before adding sources, it will do nothing
     since no fields will be found.   This is almost certainly not
     what the user wants. */
  if (!components_allocated)
    abort("allocate field components (by adding sources) before adding dft objects");
  if (!include_dV_and_interp_weights && sqrt_dV_and_interp_weights)
    abort("include_dV_and_interp_weights must be true for sqrt_dV_and_interp_weights=true in "
          "add_dft");

  dft_chunk_data data;
  data.c = c;
  data.vc = vc;
  if (Nfreq <= 1) freq_min = freq_max = (freq_min + freq_max) * 0.5;
  data.omega_min = freq_min * 2 * pi;
  data.domega = Nfreq <= 1 ? 0.0 : (freq_max * 2 * pi - data.omega_min) / (Nfreq - 1);
  data.Nomega = Nfreq;
  data.stored_weight = stored_weight;
  data.extra_weight = extra_weight;
  data.dt_factor = dt / sqrt(2.0 * pi);
  data.include_dV_and_interp_weights = include_dV_and_interp_weights;
  data.sqrt_dV_and_interp_weights = sqrt_dV_and_interp_weights;
  data.empty_dim[0] = data.empty_dim[1] = data.empty_dim[2] = data.empty_dim[3] = data.empty_dim[4] = false;
  LOOP_OVER_DIRECTIONS(where.dim, d) {
    data.empty_dim[d] = where.in_direction(d) == 0;
  }
  data.dft_chunks = chunk_next;
  loop_in_chunks(add_dft_chunkloop, (void *)&data, where, use_centered_grid ? Centered : c);

  return data.dft_chunks;
}

dft_chunk *fields::add_dft(const volume_list *where, double freq_min, double freq_max, int Nfreq,
                           bool include_dV_and_interp_weights) {
  dft_chunk *chunks = 0;
  while (where) {
    if (is_derived(where->c)) abort("derived_component invalid for dft");
    cdouble stored_weight = where->weight;
    chunks = add_dft(component(where->c), where->v, freq_min, freq_max, Nfreq,
                     include_dV_and_interp_weights, stored_weight, chunks);
    where = where->next;
  }
  return chunks;
}

dft_chunk *fields::add_dft_pt(component c, const vec &where, double freq_min, double freq_max,
                              int Nfreq) {
  return add_dft(c, where, freq_min, freq_max, Nfreq, false);
}

void fields::update_dfts() {
  am_now_working_on(FourierTransforming);
  for (int i = 0; i < num_chunks; i++)
    if (chunks[i]->is_mine()) chunks[i]->update_dfts(time(), time() - 0.5 * dt);
  finished_working();
}

void fields_chunk::update_dfts(double timeE, double timeH) {
  if (doing_solve_cw) return;
  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_chunk) {
    cur->update_dft(is_magnetic(cur->c) ? timeH : timeE);
  }
}

void dft_chunk::update_dft(double time) {
  if (!fc->f[c][0]) return;

  for (int i = 0; i < Nomega; ++i)
    dft_phase[i] = polar(1.0, (omega_min + i * domega) * time) * scale;

  int numcmp = fc->f[c][1] ? 2 : 1;

  size_t idx_dft = 0;
  LOOP_OVER_IVECS(fc->gv, is, ie, idx) {
    double w;
    if (include_dV_and_interp_weights) {
      w = IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2);
      if (sqrt_dV_and_interp_weights) w = sqrt(w);
    } else
      w = 1.0;
    double f[2]; // real/imag field value at epsilon point
    if (avg2)
      for (int cmp = 0; cmp < numcmp; ++cmp)
        f[cmp] = (w * 0.25) * (fc->f[c][cmp][idx] + fc->f[c][cmp][idx + avg1] +
                               fc->f[c][cmp][idx + avg2] + fc->f[c][cmp][idx + (avg1 + avg2)]);
    else if (avg1)
      for (int cmp = 0; cmp < numcmp; ++cmp)
        f[cmp] = (w * 0.5) * (fc->f[c][cmp][idx] + fc->f[c][cmp][idx + avg1]);
    else
      for (int cmp = 0; cmp < numcmp; ++cmp)
        f[cmp] = w * fc->f[c][cmp][idx];

    if (numcmp == 2) {
      complex<realnum> fc(f[0], f[1]);
      for (int i = 0; i < Nomega; ++i)
        dft[Nomega * idx_dft + i] += dft_phase[i] * fc;
    } else {
      realnum fr = f[0];
      for (int i = 0; i < Nomega; ++i)
        dft[Nomega * idx_dft + i] += dft_phase[i] * fr;
    }
    idx_dft++;
  }
}

void dft_chunk::scale_dft(complex<double> scale) {
  for (size_t i = 0; i < N * Nomega; ++i)
    dft[i] *= scale;
  if (next_in_dft) next_in_dft->scale_dft(scale);
}

void dft_chunk::operator-=(const dft_chunk &chunk) {
  if (c != chunk.c || N * Nomega != chunk.N * chunk.Nomega)
    abort("Mismatched chunks in dft_chunk::operator-=");

  for (size_t i = 0; i < N * Nomega; ++i)
    dft[i] -= chunk.dft[i];

  if (next_in_dft) {
    if (!chunk.next_in_dft) abort("Mismatched chunk lists in dft_chunk::operator-=");
    *next_in_dft -= *chunk.next_in_dft;
  }
}

size_t dft_chunks_Ntotal(dft_chunk *dft_chunks, size_t *my_start) {
  size_t n = 0;
  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_dft)
    n += cur->N * cur->Nomega * 2;
  *my_start = partial_sum_to_all(n) - n; // sum(n) for processes before this
  return sum_to_all(n);
}

// Note: the file must have been created in parallel mode, typically via fields::open_h5file.
void save_dft_hdf5(dft_chunk *dft_chunks, const char *name, h5file *file, const char *dprefix) {
  size_t istart;
  size_t n = dft_chunks_Ntotal(dft_chunks, &istart);

  char dataname[1024];
  snprintf(dataname, 1024,
           "%s%s"
           "%s_dft",
           dprefix ? dprefix : "", dprefix && dprefix[0] ? "_" : "", name);
  file->create_data(dataname, 1, &n);

  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_dft) {
    size_t Nchunk = cur->N * cur->Nomega * 2;
    file->write_chunk(1, &istart, &Nchunk, (realnum *)cur->dft);
    istart += Nchunk;
  }
  file->done_writing_chunks();
}

void save_dft_hdf5(dft_chunk *dft_chunks, component c, h5file *file, const char *dprefix) {
  save_dft_hdf5(dft_chunks, component_name(c), file, dprefix);
}

void load_dft_hdf5(dft_chunk *dft_chunks, const char *name, h5file *file, const char *dprefix) {
  size_t istart;
  size_t n = dft_chunks_Ntotal(dft_chunks, &istart);

  char dataname[1024];
  snprintf(dataname, 1024,
           "%s%s"
           "%s_dft",
           dprefix ? dprefix : "", dprefix && dprefix[0] ? "_" : "", name);
  int file_rank;
  size_t file_dims;
  file->read_size(dataname, &file_rank, &file_dims, 1);
  if (file_rank != 1 || file_dims != n)
    abort("incorrect dataset size (%zd vs. %zd) in load_dft_hdf5 %s:%s", file_dims, n,
          file->file_name(), dataname);

  for (dft_chunk *cur = dft_chunks; cur; cur = cur->next_in_dft) {
    size_t Nchunk = cur->N * cur->Nomega * 2;
    file->read_chunk(1, &istart, &Nchunk, (realnum *)cur->dft);
    istart += Nchunk;
  }
}

void load_dft_hdf5(dft_chunk *dft_chunks, component c, h5file *file, const char *dprefix) {
  load_dft_hdf5(dft_chunks, component_name(c), file, dprefix);
}

dft_flux::dft_flux(const component cE_, const component cH_, dft_chunk *E_, dft_chunk *H_,
                   double fmin, double fmax, int Nf, const volume &where_,
                   direction normal_direction_, bool use_symmetry_)
    : Nfreq(Nf), E(E_), H(H_), cE(cE_), cH(cH_), where(where_), normal_direction(normal_direction_),
      use_symmetry(use_symmetry_) {
  if (Nf <= 1) fmin = fmax = (fmin + fmax) * 0.5;
  freq_min = fmin;
  dfreq = Nf <= 1 ? 0.0 : (fmax - fmin) / (Nf - 1);
}

dft_flux::dft_flux(const dft_flux &f) : where(f.where) {
  freq_min = f.freq_min;
  Nfreq = f.Nfreq;
  dfreq = f.dfreq;
  E = f.E;
  H = f.H;
  cE = f.cE;
  cH = f.cH;
  normal_direction = f.normal_direction;
  use_symmetry = f.use_symmetry;
}

double *dft_flux::flux() {
  double *F = new double[Nfreq];
  for (int i = 0; i < Nfreq; ++i)
    F[i] = 0;
  for (dft_chunk *curE = E, *curH = H; curE && curH;
       curE = curE->next_in_dft, curH = curH->next_in_dft)
    for (size_t k = 0; k < curE->N; ++k)
      for (int i = 0; i < Nfreq; ++i)
        F[i] += real(curE->dft[k * Nfreq + i] * conj(curH->dft[k * Nfreq + i]));
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

void dft_flux::save_hdf5(fields &f, const char *fname, const char *dprefix, const char *prefix) {
  h5file *ff = f.open_h5file(fname, h5file::WRITE, prefix);
  save_hdf5(ff, dprefix);
  delete ff;
}

void dft_flux::load_hdf5(fields &f, const char *fname, const char *dprefix, const char *prefix) {
  h5file *ff = f.open_h5file(fname, h5file::READONLY, prefix);
  load_hdf5(ff, dprefix);
  delete ff;
}

void dft_flux::scale_dfts(complex<double> scale) {
  if (E) E->scale_dft(scale);
  if (H) H->scale_dft(scale);
}

dft_flux fields::add_dft_flux(const volume_list *where_, double freq_min, double freq_max,
                              int Nfreq, bool use_symmetry) {
  if (!where_) // handle empty list of volumes
    return dft_flux(Ex, Hy, NULL, NULL, freq_min, freq_max, Nfreq, v, NO_DIRECTION, use_symmetry);

  dft_chunk *E = 0, *H = 0;
  component cE[2] = {Ex, Ey}, cH[2] = {Hy, Hx};

  // the dft_flux object needs to store the (unreduced) volume for
  // mode-coefficient computation in mpb.cpp, but this only works
  // when the volume_list consists of a single volume, so it suffices
  // to store the first volume in the list.
  volume firstvol(where_->v);

  volume_list *where = use_symmetry ? S.reduce(where_) : new volume_list(where_);
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
      E = add_dft(cE[i], where->v, freq_min, freq_max, Nfreq, true,
                  where->weight * double(1 - 2 * i), E);
      H = add_dft(cH[i], where->v, freq_min, freq_max, Nfreq, false, 1.0, H);
    }

    where = where->next;
  }
  delete where_save;

  // if the volume list has only one entry, store its component's direction.
  // if the volume list has > 1 entry, store NO_DIRECTION.
  direction flux_dir = (where_->next ? NO_DIRECTION : component_direction(where_->c));
  return dft_flux(cE[0], cH[0], E, H, freq_min, freq_max, Nfreq, firstvol, flux_dir, use_symmetry);
}

dft_energy::dft_energy(dft_chunk *E_, dft_chunk *H_, dft_chunk *D_, dft_chunk *B_, double fmin,
                       double fmax, int Nf, const volume &where_)
    : Nfreq(Nf), E(E_), H(H_), D(D_), B(B_), where(where_) {
  if (Nf <= 1) fmin = fmax = (fmin + fmax) * 0.5;
  freq_min = fmin;
  dfreq = Nf <= 1 ? 0.0 : (fmax - fmin) / (Nf - 1);
}

dft_energy::dft_energy(const dft_energy &f) : where(f.where) {
  freq_min = f.freq_min;
  Nfreq = f.Nfreq;
  dfreq = f.dfreq;
  E = f.E;
  H = f.H;
  D = f.D;
  B = f.B;
}

double *dft_energy::electric() {
  double *F = new double[Nfreq];
  for (int i = 0; i < Nfreq; ++i)
    F[i] = 0;
  for (dft_chunk *curE = E, *curD = D; curE && curD;
       curE = curE->next_in_dft, curD = curD->next_in_dft)
    for (size_t k = 0; k < curE->N; ++k)
      for (int i = 0; i < Nfreq; ++i)
        F[i] += 0.5 * real(conj(curE->dft[k * Nfreq + i]) * curD->dft[k * Nfreq + i]);
  double *Fsum = new double[Nfreq];
  sum_to_all(F, Fsum, Nfreq);
  delete[] F;
  return Fsum;
}

double *dft_energy::magnetic() {
  double *F = new double[Nfreq];
  for (int i = 0; i < Nfreq; ++i)
    F[i] = 0;
  for (dft_chunk *curH = H, *curB = B; curH && curB;
       curH = curH->next_in_dft, curB = curB->next_in_dft)
    for (size_t k = 0; k < curH->N; ++k)
      for (int i = 0; i < Nfreq; ++i)
        F[i] += 0.5 * real(conj(curH->dft[k * Nfreq + i]) * curB->dft[k * Nfreq + i]);
  double *Fsum = new double[Nfreq];
  sum_to_all(F, Fsum, Nfreq);
  delete[] F;
  return Fsum;
}

double *dft_energy::total() {
  double *Fe = electric();
  double *Fm = magnetic();
  double *F = new double[Nfreq];
  for (int i = 0; i < Nfreq; ++i)
    F[i] = Fe[i] + Fm[i];
  delete[] Fe;
  delete[] Fm;
  return F;
}

dft_energy fields::add_dft_energy(const volume_list *where_, double freq_min, double freq_max,
                                  int Nfreq) {

  if (!where_) // handle empty list of volumes
    return dft_energy(NULL, NULL, NULL, NULL, freq_min, freq_max, Nfreq, v);

  dft_chunk *E = 0, *D = 0, *H = 0, *B = 0;
  volume firstvol(where_->v);
  volume_list *where = new volume_list(where_);
  volume_list *where_save = where;
  while (where) {
    LOOP_OVER_FIELD_DIRECTIONS(gv.dim, d) {
      E = add_dft(direction_component(Ex, d), where->v, freq_min, freq_max, Nfreq, true, 1.0, E);
      D = add_dft(direction_component(Dx, d), where->v, freq_min, freq_max, Nfreq, false, 1.0, D);
      H = add_dft(direction_component(Hx, d), where->v, freq_min, freq_max, Nfreq, true, 1.0, H);
      B = add_dft(direction_component(Bx, d), where->v, freq_min, freq_max, Nfreq, false, 1.0, B);
    }
    where = where->next;
  }
  delete where_save;

  return dft_energy(E, H, D, B, freq_min, freq_max, Nfreq, firstvol);
}

void dft_energy::save_hdf5(h5file *file, const char *dprefix) {
  save_dft_hdf5(E, "E", file, dprefix);
  file->prevent_deadlock(); // hackery
  save_dft_hdf5(D, "D", file, dprefix);
  file->prevent_deadlock(); // hackery
  save_dft_hdf5(H, "H", file, dprefix);
  file->prevent_deadlock(); // hackery
  save_dft_hdf5(B, "B", file, dprefix);
}

void dft_energy::load_hdf5(h5file *file, const char *dprefix) {
  load_dft_hdf5(E, "E", file, dprefix);
  file->prevent_deadlock(); // hackery
  load_dft_hdf5(D, "D", file, dprefix);
  file->prevent_deadlock(); // hackery
  load_dft_hdf5(H, "H", file, dprefix);
  file->prevent_deadlock(); // hackery
  load_dft_hdf5(B, "B", file, dprefix);
}

void dft_energy::save_hdf5(fields &f, const char *fname, const char *dprefix, const char *prefix) {
  h5file *ff = f.open_h5file(fname, h5file::WRITE, prefix);
  save_hdf5(ff, dprefix);
  delete ff;
}

void dft_energy::load_hdf5(fields &f, const char *fname, const char *dprefix, const char *prefix) {
  h5file *ff = f.open_h5file(fname, h5file::READONLY, prefix);
  load_hdf5(ff, dprefix);
  delete ff;
}

void dft_energy::scale_dfts(complex<double> scale) {
  if (E) E->scale_dft(scale);
  if (D) D->scale_dft(scale);
  if (H) H->scale_dft(scale);
  if (B) B->scale_dft(scale);
}

void dft_energy::remove() {
  while (E) {
    dft_chunk *nxt = E->next_in_dft;
    delete E;
    E = nxt;
  }
  while (D) {
    dft_chunk *nxt = D->next_in_dft;
    delete D;
    D = nxt;
  }
  while (H) {
    dft_chunk *nxt = H->next_in_dft;
    delete H;
    H = nxt;
  }
  while (B) {
    dft_chunk *nxt = B->next_in_dft;
    delete B;
    B = nxt;
  }
}

direction fields::normal_direction(const volume &where) const {
  direction d = where.normal_direction();
  if (d == NO_DIRECTION) {
    /* hack so that we still infer the normal direction correctly for
       volumes with empty dimensions */
    volume where_pad(where);
    LOOP_OVER_DIRECTIONS(where.dim, d1) {
      if (nosize_direction(d1) && where.in_direction(d1) == 0.0)
        where_pad.set_direction_max(d1, where.in_direction_min(d1) + 0.1);
    }
    d = where_pad.normal_direction();
    if (d == NO_DIRECTION) abort("Could not determine normal direction for given grid_volume.");
  }
  return d;
}

dft_flux fields::add_dft_flux(direction d, const volume &where, double freq_min, double freq_max,
                              int Nfreq, bool use_symmetry) {
  if (d == NO_DIRECTION) d = normal_direction(where);
  volume_list vl(where, direction_component(Sx, d));
  dft_flux flux = add_dft_flux(&vl, freq_min, freq_max, Nfreq, use_symmetry);
  flux.normal_direction = d;
  return flux;
}

dft_flux fields::add_mode_monitor(direction d, const volume &where, double freq_min,
                                  double freq_max, int Nfreq) {
  return add_dft_flux(d, where, freq_min, freq_max, Nfreq, /*use_symmetry=*/false);
}

dft_flux fields::add_dft_flux_box(const volume &where, double freq_min, double freq_max,
                                  int Nfreq) {
  volume_list *faces = 0;
  LOOP_OVER_DIRECTIONS(where.dim, d) {
    if (where.in_direction(d) > 0) {
      volume face(where);
      derived_component c = direction_component(Sx, d);
      face.set_direction_min(d, where.in_direction_max(d));
      faces = new volume_list(face, c, +1, faces);
      face.set_direction_min(d, where.in_direction_min(d));
      face.set_direction_max(d, where.in_direction_min(d));
      faces = new volume_list(face, c, -1, faces);
    }
  }

  dft_flux flux = add_dft_flux(faces, freq_min, freq_max, Nfreq);
  delete faces;
  return flux;
}

dft_flux fields::add_dft_flux_plane(const volume &where, double freq_min, double freq_max,
                                    int Nfreq) {
  return add_dft_flux(NO_DIRECTION, where, freq_min, freq_max, Nfreq);
}

dft_fields::dft_fields(dft_chunk *chunks_, double freq_min_, double freq_max_, int Nfreq_,
                       const volume &where_)
    : where(where_) {
  chunks = chunks_;
  freq_min = freq_min_;
  dfreq = Nfreq_ <= 1 ? 0.0 : (freq_max_ - freq_min_) / (Nfreq_ - 1);
  Nfreq = Nfreq_;
}

void dft_fields::scale_dfts(cdouble scale) { chunks->scale_dft(scale); }

void dft_fields::remove() {
  while (chunks) {
    dft_chunk *nxt = chunks->next_in_dft;
    delete chunks;
    chunks = nxt;
  }
}

dft_fields fields::add_dft_fields(component *components, int num_components, const volume where,
                                  double freq_min, double freq_max, int Nfreq) {
  bool include_dV_and_interp_weights = false;
  cdouble stored_weight = 1.0;
  dft_chunk *chunks = 0;
  for (int nc = 0; nc < num_components; nc++)
    chunks = add_dft(components[nc], where, freq_min, freq_max, Nfreq,
                     include_dV_and_interp_weights, stored_weight, chunks);

  return dft_fields(chunks, freq_min, freq_max, Nfreq, where);
}

/***************************************************************/
/* chunk-level processing for fields::process_dft_component.   */
/***************************************************************/
cdouble dft_chunk::process_dft_component(int rank, direction *ds, ivec min_corner, ivec max_corner,
                                         int num_freq, h5file *file, double *buffer, int reim,
                                         cdouble *field_array, void *mode1_data, void *mode2_data,
                                         component c_conjugate, bool retain_interp_weights) {
  /*****************************************************************/
  /* compute the size of the chunk we own and its strides etc.     */
  /*****************************************************************/
  size_t start[3] = {0, 0, 0};
  size_t file_count[3] = {1, 1, 1}, array_count[3] = {1, 1, 1};
  int file_offset[3] = {0, 0, 0};
  int file_stride[3] = {1, 1, 1};
  ivec isS = S.transform(is, sn) + shift;
  ivec ieS = S.transform(ie, sn) + shift;

  ivec permute(zero_ivec(fc->gv.dim));
  for (int i = 0; i < 3; ++i)
    permute.set_direction(fc->gv.yucky_direction(i), i);
  permute = S.transform_unshifted(permute, sn);
  LOOP_OVER_DIRECTIONS(permute.dim, d) { permute.set_direction(d, abs(permute.in_direction(d))); }

  for (int i = 0; i < rank; ++i) {
    direction d = ds[i];
    int isd = isS.in_direction(d), ied = ieS.in_direction(d);
    start[i] = (min(isd, ied) - min_corner.in_direction(d)) / 2;
    file_count[i] = abs(ied - isd) / 2 + 1;
    if (ied < isd) file_offset[permute.in_direction(d)] = file_count[i] - 1;
    array_count[i] = (max_corner.in_direction(d) - min_corner.in_direction(d)) / 2 + 1;
  }

  for (int i = 0; i < rank; ++i) {
    direction d = ds[i];
    int j = permute.in_direction(d);
    for (int k = i + 1; k < rank; ++k)
      file_stride[j] *= file_count[k];
    file_offset[j] *= file_stride[j];
    if (file_offset[j]) file_stride[j] *= -1;
  }

  /*****************************************************************/
  /* For collapsing empty dimensions, we want to retain interpolation
     weights for empty dimensions, but not interpolation weights for
     integration of edge pixels (for retain_interp_weights == true).
     All of the weights are stored in (s0, s1, e0, e1), so we make
     a copy of these with the weights for non-empty dimensions set to 1. */
  vec s0i(s0), s1i(s1), e0i(e0), e1i(e1);
  LOOP_OVER_DIRECTIONS(fc->gv.dim, d) {
    if (!empty_dim[d]) {
      s0i.set_direction(d, 1.0);
      s1i.set_direction(d, 1.0);
      e0i.set_direction(d, 1.0);
      e1i.set_direction(d, 1.0);
    }
  }

  /***************************************************************/
  /* loop over all grid points in our piece of the volume        */
  /***************************************************************/
  vec rshift(shift * (0.5 * fc->gv.inva));
  int chunk_idx = 0;
  cdouble integral = 0.0;
  LOOP_OVER_IVECS(fc->gv, is, ie, idx) {
    IVEC_LOOP_LOC(fc->gv, loc);
    loc = S.transform(loc, sn) + rshift;
    double w = IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2), interp_w = 1.0;
    cdouble dft_val = dft[Nomega * (chunk_idx++) + num_freq] / stored_weight;
    if (include_dV_and_interp_weights) dft_val /= (sqrt_dV_and_interp_weights ? sqrt(w) : w);

    if (retain_interp_weights)
      interp_w = IVEC_LOOP_WEIGHT(s0i, s1i, e0i, e1i, 1.0);

    cdouble mode1val = 0.0, mode2val = 0.0;
    if (mode1_data) mode1val = eigenmode_amplitude(mode1_data, loc, S.transform(c_conjugate, sn));
    if (mode2_data) mode2val = eigenmode_amplitude(mode2_data, loc, S.transform(c, sn));

    if (file) {
      int idx2 = ((((file_offset[0] + file_offset[1] + file_offset[2]) + loop_i1 * file_stride[0]) +
                   loop_i2 * file_stride[1]) +
                  loop_i3 * file_stride[2]);

      dft_val *= interp_w;

      cdouble val = (mode1_data ? mode1val : dft_val);
      buffer[idx2] = reim ? imag(val) : real(val);
    } else if (field_array) {
      IVEC_LOOP_ILOC(fc->gv, iloc);         // iloc <-- indices of parent point in Yee grid
      iloc = S.transform(iloc, sn) + shift; // iloc <-- indices of child point in Yee grid
      iloc -= min_corner;                   // iloc <-- 2*(indices of point in DFT array)

      // the index of point n1 or (n1,n2) or (n1,n2,n3) in a 1D, 2D, or 3D array is
      // (for a 1D array) n1
      // (for a 2D array) n2 + n1*N2
      // (for a 3D array) n3 + n2*N3 + n1*N2*N3
      // where NI = number of points in Ith direction.
      int idx2 = 0;
      for (int i = rank - 1, stride = 1; i >= 0; stride *= array_count[i--])
        idx2 += stride * (iloc.in_direction(ds[i]) / 2);
      field_array[idx2] = interp_w * dft_val;
    } else {
      mode1val = conj(mode1val); // conjugated inner product
      if (mode2_data)
        integral += w * mode1val * mode2val;
      else
        integral += w * mode1val * dft_val;
    }

  } // LOOP_OVER_IVECS(fc->gv, is, ie, idx)

  if (file) file->write_chunk(rank, start, file_count, buffer);

  return integral;
}

/***************************************************************/
/* low-level [actually intermediate-level, since it calls      */
/* dft_chunk::process_dft_component(), which is the true       */
/* low-level function] workhorse routine that forms the common */
/* backend for several operations involving DFT fields.        */
/*                                                             */
/* looks through the given collection of dft_chunks and        */
/* processes only those chunks that store component c, using   */
/* only data for frequency #num_freq.                          */
/*                                                             */
/* the meaning of 'processes' depends on the arguments:        */
/*                                                             */
/*  1. if HDF5FileName is non-null: write to the given HDF5    */
/*     file a new dataset describing either                    */
/*      (A) DFT field component c (if mode_data1 is null), or  */
/*      (B) mode field component c for the eigenmode described */
/*          by mode_data1 (if it is non-null)                  */
/*                                                             */
/*  2. if HDF5FileName is null but pfield_array is non-null:   */
/*     set *pfield_array equal to a newly allocated buffer     */
/*     populated on return with values of DFT field component  */
/*     c, equivalent to writing the data to HDF5 and reading   */
/*     it back into field_array                                */
/*                                                             */
/*  3. if both HDF5FileName and field_array are null: compute  */
/*     and return an  overlap integral between                 */
/*      (A) the DFT fields and the fields of the eigenmode     */
/*          described by mode_data1 (if mode_data2 is null)    */
/*      (B) the eigenmode fields described by mode_data1       */
/*          and the eigenmode fields described by mode_data2   */
/*          (if mode_data2 is non-null).                       */
/*     more specifically, the integral computed is             */
/*      < mode1_{c_conjugate} | dft_{c} >                      */
/*     in case (A) and                                         */
/*      < mode1_{c_conjugate} | mode2_{c} >                    */
/*     in case (B).                                            */
/*                                                             */
/* if where is non-null, only field components inside *where   */
/* are processed.                                              */
/***************************************************************/
cdouble fields::process_dft_component(dft_chunk **chunklists, int num_chunklists, int num_freq,
                                      component c, const char *HDF5FileName, cdouble **pfield_array,
                                      int *array_rank, int *array_dims, void *mode1_data,
                                      void *mode2_data, component c_conjugate,
                                      bool *first_component, bool retain_interp_weights) {
  /***************************************************************/
  /* get statistics on the volume slice **************************/
  /***************************************************************/
  volume *where = &v; // use full volume of fields
  size_t bufsz = 0;
  ivec min_corner = gv.round_vec(where->get_max_corner()) + one_ivec(gv.dim);
  ivec max_corner = gv.round_vec(where->get_min_corner()) - one_ivec(gv.dim);
  for (int ncl = 0; ncl < num_chunklists; ncl++)
    for (dft_chunk *chunk = chunklists[ncl]; chunk; chunk = chunk->next_in_dft) {
      if (chunk->c != c) continue;
      ivec isS = chunk->S.transform(chunk->is, chunk->sn) + chunk->shift;
      ivec ieS = chunk->S.transform(chunk->ie, chunk->sn) + chunk->shift;
      min_corner = min(min_corner, min(isS, ieS));
      max_corner = max(max_corner, max(isS, ieS));
      size_t this_bufsz = 1;
      LOOP_OVER_DIRECTIONS(chunk->fc->gv.dim, d) {
        this_bufsz *= (chunk->ie.in_direction(d) - chunk->is.in_direction(d)) / 2 + 1;
      }
      bufsz = max(bufsz, this_bufsz);
    }
  max_corner = max_to_all(max_corner);
  min_corner = -max_to_all(-min_corner); // i.e., min_to_all

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int rank = 0;
  size_t dims[3];
  direction ds[3];
  size_t array_size = 1;
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    if (rank >= 3) abort("too many dimensions in process_dft_component");
    size_t n = std::max(0, (max_corner.in_direction(d) - min_corner.in_direction(d)) / 2 + 1);

    if (n > 1) {
      ds[rank] = d;
      dims[rank++] = n;
      array_size *= n;
    }
  }
  if (array_rank) {
    *array_rank = rank;
    for (int d = 0; d < rank; d++)
      array_dims[d] = dims[d];
  }
  if (rank == 0) {
    if (pfield_array) *pfield_array = 0;
    return 0.0; // no chunks with the specified component on this processor
  }

  /***************************************************************/
  /* buffer for process-local contributions to HDF5 output files,*/
  /* like h5_output_data::buf in h5fields.cpp                    */
  /***************************************************************/
  realnum *buffer = 0;
  cdouble *field_array = 0;
  int reim_max = 0;
  if (HDF5FileName) {
    buffer = new realnum[bufsz];
    reim_max = 1;
  } else if (pfield_array)
    *pfield_array = field_array = (array_size ? new cdouble[array_size] : 0);

  bool append_data = false;
  bool single_precision = false;
  cdouble overlap = 0.0;
  for (int reim = 0; reim <= reim_max; reim++) {
    h5file *file = 0;
    if (HDF5FileName) {
      file = open_h5file(HDF5FileName, (*first_component) ? h5file::WRITE : h5file::READWRITE);
      *first_component = false;
      char dataname[100];
      snprintf(dataname, 100, "%s_%i.%c", component_name(c), num_freq, reim ? 'i' : 'r');
      file->create_or_extend_data(dataname, rank, dims, append_data, single_precision);
    }

    for (int ncl = 0; ncl < num_chunklists; ncl++)
      for (dft_chunk *chunk = chunklists[ncl]; chunk; chunk = chunk->next_in_dft)
        if (chunk->c == c)
          overlap +=
              chunk->process_dft_component(rank, ds, min_corner, max_corner, num_freq, file, buffer,
                                           reim, field_array, mode1_data, mode2_data, c_conjugate, retain_interp_weights);

    if (HDF5FileName) {
      file->done_writing_chunks();
      file->prevent_deadlock(); // hackery
      delete file;
    } else if (field_array) {
/***************************************************************/
/* repeatedly call sum_to_all to consolidate full field array  */
/* on all cores                                                */
/***************************************************************/
#define BUFSIZE 1 << 16 // use 64k buffer
      cdouble *buf = new cdouble[BUFSIZE];
      ptrdiff_t offset = 0;
      size_t remaining = array_size;
      while (remaining != 0) {
        size_t size = (remaining > BUFSIZE ? BUFSIZE : remaining);
        sum_to_all(field_array + offset, buf, size);
        memcpy(field_array + offset, buf, size * sizeof(cdouble));
        remaining -= size;
        offset += size;
      }
      delete[] buf;
    }
  } // for(int reim=0; reim<=reim_max; reim++)

  if (HDF5FileName)
    delete[] buffer;
  else
    overlap = sum_to_all(overlap);

  return overlap;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool increment(int n[3], int nMax[3], int rank) {
  for (rank--; rank >= 0; rank--)
    if (++n[rank] < nMax[rank])
      return false;
    else
      n[rank] = 0;
  return true;
}

cdouble *collapse_empty_dimensions(cdouble *array, int *rank, int dims[3], volume dft_volume) {
  /*--------------------------------------------------------------*/
  /*- detect empty dimensions and compute rank and strides for    */
  /*- collapsed array                                             */
  /*--------------------------------------------------------------*/
  int full_rank = *rank;
  if (full_rank == 0) return array;

  int reduced_rank = 0, reduced_dims[3], reduced_stride[3] = {1, 1, 1}, nd = 0;
  LOOP_OVER_DIRECTIONS(dft_volume.dim, d) {
    int dim = dims[nd++];
    if (dim == 0) continue;
    if (dft_volume.in_direction(d) == 0.0)
      reduced_stride[nd - 1] = 0; // degenerate dimension, to be collapsed
    else
      reduced_dims[reduced_rank++] = dim;
  }
  if (reduced_rank == full_rank) return array; // nothing to collapse

  /*--------------------------------------------------------------*/
  /*- set up strides into full and reduced arrays                 */
  /*--------------------------------------------------------------*/
  int stride[3] = {1, 1, 1}; // non-reduced array strides
  if (full_rank == 2)
    stride[0] = dims[1]; // rstride is already all set in this case
  else if (full_rank == 3) {
    stride[0] = dims[1] * dims[2];
    stride[1] = dims[2];
    if (reduced_stride[0] != 0)
      reduced_stride[0] = reduced_dims[1];
    else if (reduced_stride[1] != 0)
      reduced_stride[1] = reduced_dims[1];
    // else: two degenerate dimensions->reduced array is 1-diml, no strides needed
  }

  /*--------------------------------------------------------------*/
  /*- allocate reduced array and compress full array into it     -*/
  /*--------------------------------------------------------------*/
  int reduced_size = reduced_dims[0] * (reduced_rank == 2 ? reduced_dims[1] : 1);
  cdouble *reduced_array = new cdouble[reduced_size];
  if (!reduced_array) abort("%s:%i: out of memory (%i)", __FILE__, __LINE__, reduced_size);
  memset(reduced_array, 0, reduced_size * sizeof(cdouble));

  int n[3] = {0, 0, 0};
  do {
    int index = n[0] * stride[0] + n[1] * stride[1] + n[2] * stride[2];
    int rindex = n[0] * reduced_stride[0] + n[1] * reduced_stride[1] + n[2] * reduced_stride[2];
    reduced_array[rindex] += array[index];
  } while (!increment(n, dims, full_rank));

  *rank = reduced_rank;
  dims[0] = reduced_dims[0];
  if (reduced_rank == 2) dims[1] = reduced_dims[1];
  delete[] array;
  return reduced_array;
}

/***************************************************************/
/* routines for fetching arrays of dft fields                  */
/***************************************************************/
cdouble *fields::get_dft_array(dft_flux flux, component c, int num_freq, int *rank, int dims[3]) {
  dft_chunk *chunklists[2];
  chunklists[0] = flux.E;
  chunklists[1] = flux.H;
  cdouble *array;
  process_dft_component(chunklists, 2, num_freq, c, 0, &array, rank, dims);
  return collapse_empty_dimensions(array, rank, dims, flux.where);
}

cdouble *fields::get_dft_array(dft_force force, component c, int num_freq, int *rank, int dims[3]) {
  dft_chunk *chunklists[3];
  chunklists[0] = force.offdiag1;
  chunklists[1] = force.offdiag2;
  chunklists[2] = force.diag;
  cdouble *array;
  process_dft_component(chunklists, 3, num_freq, c, 0, &array, rank, dims);
  return collapse_empty_dimensions(array, rank, dims, force.where);
}

cdouble *fields::get_dft_array(dft_near2far n2f, component c, int num_freq, int *rank,
                               int dims[3]) {
  dft_chunk *chunklists[1];
  chunklists[0] = n2f.F;
  cdouble *array;
  process_dft_component(chunklists, 1, num_freq, c, 0, &array, rank, dims);
  return collapse_empty_dimensions(array, rank, dims, n2f.where);
}

cdouble *fields::get_dft_array(dft_fields fdft, component c, int num_freq, int *rank, int dims[3]) {
  dft_chunk *chunklists[1];
  chunklists[0] = fdft.chunks;
  cdouble *array;
  process_dft_component(chunklists, 1, num_freq, c, 0, &array, rank, dims);
  return collapse_empty_dimensions(array, rank, dims, fdft.where);
}

/***************************************************************/
/* wrapper around process_dft_component that writes HDF5       */
/* datasets for all components at all frequencies stored in    */
/* the given collection of DFT chunks                          */
/***************************************************************/
void fields::output_dft_components(dft_chunk **chunklists, int num_chunklists, volume dft_volume,
                                   const char *HDF5FileName) {
  int NumFreqs = 0;
  for (int nc = 0; nc < num_chunklists && NumFreqs == 0; nc++)
    if (chunklists[nc]) NumFreqs = chunklists[nc]->Nomega;

  // if the volume has zero thickness in one or more directions, the DFT
  // grid is two pixels thick in those directions, but we want the HDF5 output
  // to be just one pixel thick in those directions. solution: first get the
  // fields in array form (as get_dft_array), then collapse degenerate dimensions
  // and export the collapsed array to HDF5. in this case the max_to_all() below
  // is needed to make sure everybody agrees on how many frequencies there are,
  // because some processes' field chunks may have no overlap with dft_volume,
  // in which case those processes will think NumFreqs==0.
  bool have_empty_dims = false;
  LOOP_OVER_DIRECTIONS(dft_volume.dim, d) {
    if (dft_volume.in_direction(d) == 0.0) have_empty_dims = true;
  }

  h5file *file = 0;
  if (have_empty_dims && am_master()) {
    char filename[100];
    snprintf(filename, 100, "%s%s", HDF5FileName, strstr("%.h5", HDF5FileName) ? "" : ".h5");
    file = new h5file(filename, h5file::WRITE, false /*parallel*/);
  }
  if (have_empty_dims) NumFreqs = max_to_all(NumFreqs); // subtle!

  bool first_component = true;
  for (int num_freq = 0; num_freq < NumFreqs; num_freq++)
    FOR_E_AND_H(c) {
      if (!have_empty_dims) {
        process_dft_component(chunklists, num_chunklists, num_freq, c, HDF5FileName, 0, 0, 0, 0, 0,
                              Ex, &first_component);
      } else {
        cdouble *array = 0;
        int rank, dims[3];
        process_dft_component(chunklists, num_chunklists, num_freq, c, 0, &array, &rank, dims);
        if (rank > 0 && am_master()) {
          array = collapse_empty_dimensions(array, &rank, dims, dft_volume);
          if (rank == 0) abort("%s:%i: internal error", __FILE__, __LINE__);
          size_t stdims[2];
          stdims[0] = (size_t)dims[0];
          stdims[1] = (size_t)((rank > 1) ? dims[1] : 1);
          size_t array_size = stdims[0] * stdims[1];
          double *real_array = new double[array_size];
          if (!real_array) abort("%s:%i:out of memory(%lu)", __FILE__, __LINE__, array_size);
          for (int reim = 0; reim < 2; reim++) {
            for (size_t n = 0; n < array_size; n++)
              real_array[n] = (reim == 0 ? real(array[n]) : imag(array[n]));
            char dataname[100], filename[100];
            snprintf(dataname, 100, "%s_%i.%c", component_name(c), num_freq, reim ? 'i' : 'r');
            snprintf(filename, 100, "%s%s", HDF5FileName, strstr(".h5", HDF5FileName) ? "" : ".h5");
            bool single_precision = false;
            file->write(dataname, rank, stdims, real_array, single_precision);
          }
          delete[] real_array;
        }
        if (array) delete[] array;
      }
    }
  if (file) delete file;
}

void fields::output_dft(dft_flux flux, const char *HDF5FileName) {
  dft_chunk *chunklists[2];
  chunklists[0] = flux.E;
  chunklists[1] = flux.H;
  output_dft_components(chunklists, 2, flux.where, HDF5FileName);
}

void fields::output_dft(dft_force force, const char *HDF5FileName) {
  dft_chunk *chunklists[3];
  chunklists[0] = force.offdiag1;
  chunklists[1] = force.offdiag2;
  chunklists[2] = force.diag;
  output_dft_components(chunklists, 3, force.where, HDF5FileName);
}

void fields::output_dft(dft_near2far n2f, const char *HDF5FileName) {
  dft_chunk *chunklists[1];
  chunklists[0] = n2f.F;
  output_dft_components(chunklists, 1, n2f.where, HDF5FileName);
}

void fields::output_dft(dft_fields fdft, const char *HDF5FileName) {
  dft_chunk *chunklists[1];
  chunklists[0] = fdft.chunks;
  output_dft_components(chunklists, 1, fdft.where, HDF5FileName);
}

/***************************************************************/
/* does the same thing as output_dft(flux ...), but using      */
/* eigenmode fields instead of dft_flux fields.                */
/***************************************************************/
void fields::output_mode_fields(void *mode_data, dft_flux flux, const char *HDF5FileName) {
  h5file *file = open_h5file(HDF5FileName, h5file::WRITE);
  delete file;

  dft_chunk *chunklists[2];
  chunklists[0] = flux.E;
  chunklists[1] = flux.H;
  FOR_E_AND_H(c) { process_dft_component(chunklists, 2, 0, c, 0, 0, 0, 0, mode_data, 0, c); }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void fields::get_overlap(void *mode1_data, void *mode2_data, dft_flux flux, int num_freq,
                         cdouble overlaps[2]) {
  component cE[2], cH[2];
  switch (flux.normal_direction) {
    case X:
      cE[0] = Ey;
      cH[0] = Hz;
      cE[1] = Ez;
      cH[1] = Hy;
      break;
    case Y:
      cE[0] = Ez;
      cH[0] = Hx;
      cE[1] = Ex;
      cH[1] = Hz;
      break;
    case R:
      cE[0] = Ep;
      cH[0] = Hz;
      cE[1] = Ez;
      cH[1] = Hp;
      break;
    case P:
      cE[0] = Ez;
      cH[0] = Hr;
      cE[1] = Er;
      cH[1] = Hz;
      break;
    case Z:
      if (gv.dim == Dcyl)
        cE[0] = Er, cE[1] = Ep, cH[0] = Hp, cH[1] = Hr;
      else
        cE[0] = Ex, cE[1] = Ey, cH[0] = Hy, cH[1] = Hx;
      break;
    default: abort("invalid normal_direction in get_overlap");
  };

  dft_chunk *chunklists[2];
  chunklists[0] = flux.E;
  chunklists[1] = flux.H;
  cdouble ExHy = process_dft_component(chunklists, 2, num_freq, cE[0], 0, 0, 0, 0, mode1_data,
                                       mode2_data, cH[0]);
  cdouble EyHx = process_dft_component(chunklists, 2, num_freq, cE[1], 0, 0, 0, 0, mode1_data,
                                       mode2_data, cH[1]);
  cdouble HyEx = process_dft_component(chunklists, 2, num_freq, cH[0], 0, 0, 0, 0, mode1_data,
                                       mode2_data, cE[0]);
  cdouble HxEy = process_dft_component(chunklists, 2, num_freq, cH[1], 0, 0, 0, 0, mode1_data,
                                       mode2_data, cE[1]);
  overlaps[0] = ExHy - EyHx;
  overlaps[1] = HyEx - HxEy;
}

void fields::get_mode_flux_overlap(void *mode_data, dft_flux flux, int num_freq,
                                   std::complex<double> overlaps[2]) {
  get_overlap(mode_data, 0, flux, num_freq, overlaps);
}

void fields::get_mode_mode_overlap(void *mode1_data, void *mode2_data, dft_flux flux,
                                   std::complex<double> overlaps[2]) {
  get_overlap(mode1_data, mode2_data, flux, 0, overlaps);
}

} // namespace meep
