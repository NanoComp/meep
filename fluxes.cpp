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

#include "dactyl.h"
#include "dactyl_internals.h"

flux_plane *fields::add_flux_plane(const vec &p1, const vec &p2) {
  partial_flux_plane *hello_world = NULL;
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine()) {
      partial_flux_plane *temp = chunks[i]->new_flux_plane(p1,p2);
      if (temp) {
        temp->append(hello_world);
        hello_world = temp;
      }
    }
  flux_plane *nfp = new flux_plane(hello_world);
  nfp->next = fluxes;
  fluxes = nfp;
  return nfp;
}

partial_flux_plane *fields_chunk::new_flux_plane(const vec &p1, const vec &p2) {
  switch (v.dim) {
  case d1: return nfp_1d(p1);
  }
  return NULL;
}

partial_flux_plane *fields_chunk::nfp_1d(const vec &p1) {
  int inds[8];
  double w[8];
  v.interpolate(Hy,p1,inds,w);
  const int indhy = inds[0];
  v.interpolate(Ex,p1,inds,w);
  const int index = inds[0];
  const vec lochy = v.loc(Hy,indhy);
  const vec locex = v.loc(Ex,index);
  const double why = fabs(locex.z() - p1.z())/fabs(locex.z() - lochy.z());
  const double wex = 1.0 - why;
  partial_flux_plane *out = NULL;
  if (v.owns(v.iloc(Hy,indhy))) {
    partial_flux_plane *temp = new partial_flux_plane(this, 2);
    temp->weights[0] = temp->weights[1] = 0.5*why;
    temp->indH[0] = temp->indH[1] = indhy;
    temp->indE[0] = indhy;
    temp->indE[1] = indhy + 1;
    temp->cE = Ex;
    temp->cH = Hy;
    temp->next = out;
    temp->next_in_chunk = out;
    out = temp;
  }
  if (v.owns(v.iloc(Ex,index))) {
    partial_flux_plane *temp = new partial_flux_plane(this, 2);
    temp->weights[0] = temp->weights[1] = 0.5*wex;
    temp->indE[0] = temp->indE[1] = index;
    temp->indH[0] = index;
    temp->indH[1] = index - 1;
    temp->cE = Ex;
    temp->cH = Hy;
    temp->next = out;
    temp->next_in_chunk = out;
    out = temp;
  }
  if (out) {
    out->append_in_chunk(fluxes);
    fluxes = out;
  }
  return out;
}

void partial_flux_plane::append(partial_flux_plane *n) {
  partial_flux_plane *mover = this;
  while (mover->next) mover = mover->next;
  mover->next = n;
}

void partial_flux_plane::append_in_chunk(partial_flux_plane *n) {
  partial_flux_plane *mover = this;
  while (mover->next_in_chunk) mover = mover->next_in_chunk;
  mover->next_in_chunk = n;
}

partial_flux_plane::partial_flux_plane(fields_chunk *thech, int s) {
  numpts = s;
  next = next_in_chunk = NULL;
  weights = new double[s];
  f = thech;
  indE = new int[s];
  indH = new int[s];
  oldE[1] = NULL;
  for (int cmp=0;cmp<2;cmp++) {
    oldE[cmp] = new double[s];
    for (int i=0;i<s;i++) oldE[cmp][i] = 0.0;
  }
}

partial_flux_plane::~partial_flux_plane() {
  delete[] weights;
  delete[] indE;
  delete[] indH;
  for (int cmp=0;cmp<2;cmp++)
    delete[] oldE[cmp];
  delete next_in_chunk;
}

flux_plane::flux_plane(partial_flux_plane *p) {
  partials = p;
  next = NULL;
}

flux_plane::~flux_plane() {
  delete next;
}

double flux_plane::flux() {
  double fl = 0.0;
  if (partials) {
    const int is_real = partials->f->is_real;
    partial_flux_plane *mover = partials;
    do {
      DOCMP
        for (int n=0;n<mover->numpts;n++) {
          const double enew = mover->f->f[mover->cE][cmp][mover->indE[n]];
          const double eold = mover->oldE[cmp][n];
          fl += (enew + eold) * 0.5 *(1.0/(4.0*pi))*mover->weights[n]*
            mover->f->f[mover->cH][cmp][mover->indH[n]];
        }
      mover = mover->next;
    } while (mover);
  }
  return sum_to_all(fl);
}

void fields::update_fluxes() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->update_fluxes();
}

void fields_chunk::update_fluxes() {
  if (!fluxes) return;
  partial_flux_plane *mover = fluxes;
  do {
    DOCMP
      for (int i=0;i<mover->numpts;i++)
        mover->oldE[cmp][i] = f[mover->cE][cmp][mover->indE[i]];
    mover = mover->next_in_chunk;
  } while (mover);
}
