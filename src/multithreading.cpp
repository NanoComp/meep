/* Copyright (C) 2005-2020 Massachusetts Institute of Technology
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
#include <complex>
#include <map>
#include <limits.h>
#include <utility>
#include <vector>
#include <algorithm>

#include "multithreading.hpp"

using namespace std;

namespace meep {

/***************************************************************/
/* implementation of ivec_loop_counter                         */
/***************************************************************/
ivec_loop_counter::ivec_loop_counter(const grid_volume &gv, const ivec &_is, const ivec &_ie, int nt, int NT) {
  init(gv, _is, _ie, nt, NT);
}

void ivec_loop_counter::init(const grid_volume &gv, const ivec &_is, const ivec &_ie, int nt, int NT) {
  if (is==_is && ie==_ie && gvlc==gv.little_corner() && inva==gv.inva)
    return; // no need to reinitialize

  is   = _is;
  gvlc = gv.little_corner();
  ie   = _ie;
  inva = gv.inva;

  // tabulate 'active' dimensions, i.e. those in which loop is non-empty
  idx0 = 0;
  rank = 0;
  size_t num_iters = 1;
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    idx0 += gv.stride(d) * ((is - gvlc).in_direction(d)/2);
    ptrdiff_t count = (ie-is).in_direction(d)/2 + 1;
    if (count > 1) {
      idx_stride[rank] =  gv.stride(d);
      loop_dir[rank]   =  d;
      N[rank]          =  count;
      num_iters        *= count;
      rank++;
    }
  }
  N_inner = N[rank-1];
  idx_step = idx_stride[rank-1];

  if (NT<1) { nt=0; NT=1; }
  min_iter = (nt*num_iters) / NT;
  max_iter = ((nt+1)*num_iters) / NT;
}

/***************************************************************/
/* helper function for the 'start' routine that initializes    */
/* computation of 'k' index into PML sigma arrays for up to    */
/* two distinct directions                                     */
/***************************************************************/
void ivec_loop_counter::init_k(int nd, direction dsig) {
  if (dsig==NO_DIRECTION) {
    k0[nd]=INT_MAX;
    k_step[nd]=0;
  }
  else {
    k0[nd]=(is-gvlc).in_direction(dsig);
    for(int r=0; r<rank; r++)
      k_stride[nd][r] = (loop_dir[r]==dsig ? 2 : 0);
    k_step[nd]=k_stride[nd][rank-1];
  }
}

ptrdiff_t ivec_loop_counter::start(direction dsig1, direction dsig2) {
  complete = false;
  init_k(0, dsig1);
  init_k(1, dsig2);

  current_iter=min_iter;
  N_inner=0;
  return update_outer();
}

ptrdiff_t ivec_loop_counter::start(const grid_volume &gv, const ivec &_is, const ivec &_ie, int nt, int NT,
                                   direction dsig1, direction dsig2) {
  init(gv, _is, _ie, nt, NT);
  return start(dsig1, dsig2);
}

ptrdiff_t ivec_loop_counter::update_outer() {
  current_iter += N_inner;
  if (current_iter>=max_iter) {
    complete=true;
    return idx0;
  }
  ptrdiff_t idx = niter_to_narray(current_iter, n);
  N_inner = std::min( N[rank-1]-n[rank-1], max_iter-current_iter );
  idx_max = idx + N_inner*idx_step;
  return idx;
}

void ivec_loop_counter::advance() {
  if (++current_iter < max_iter)
    for(int r=rank-1; r>=0; n[r--]=0)
      if ( ++n[r] < N[r] )
        return;
  complete=true;
}

ptrdiff_t ivec_loop_counter::increment(int *k1, int *k2) {
  advance();
  get_k(k1,k2);
  return get_idx();
}

ptrdiff_t ivec_loop_counter::increment(size_t *k1, size_t *k2) {
  int ik1, ik2;
  ptrdiff_t idx=increment( k1 ? &ik1 : 0, k2 ? &ik2 : 0 );
  if (k1) *k1=ik1;
  if (k2) *k2=ik2;
  return idx;
}

// given the overall index of a loop iteration (between 0 and N[0]*...*N[rank-1]),
// resolve the loop indices in each dimension
ptrdiff_t ivec_loop_counter::niter_to_narray(size_t niter, size_t narray[3]) {
  ptrdiff_t idx=idx0;
  for(int r=rank-1; r>0; niter/=N[r--]) {
    narray[r] = niter % N[r];
    idx+=narray[r]*idx_stride[r];
  }
  narray[0]=niter;
  return idx+narray[0]*idx_stride[0];
}

ptrdiff_t ivec_loop_counter::get_idx(size_t *narray) {
  if (narray==0) narray=n;
  ptrdiff_t idx = idx0;
  for(int r=0; r<rank; r++)
    idx += narray[r]*idx_stride[r];
  return idx;
}

int ivec_loop_counter::get_k(int nd, size_t *narray) {
  if (k0[nd]==INT_MAX) return 0;
  if (narray==0) narray=n;
  int k=k0[nd];
  for(int r=0; r<rank; r++)
    k+=narray[r]*k_stride[nd][r];
  return k;
}

void ivec_loop_counter::get_k(int *k1, int *k2, size_t *narray) {
  if (k1) *k1=get_k(0, narray);
  if (k2) *k2=get_k(1, narray);
}

ivec ivec_loop_counter::get_iloc(size_t *narray) {
  if (narray==0) narray=n;
  ivec iloc=is;
  for(int r=0; r<rank; r++)
    iloc.set_direction(loop_dir[r], is.in_direction(loop_dir[r]) + 2*narray[r]);
  return iloc;
}

vec ivec_loop_counter::get_loc(size_t *narray) {
  ivec iloc=get_iloc(narray);
  vec loc(iloc.dim);
  LOOP_OVER_DIRECTIONS(iloc.dim,d)
    loc.set_direction(d, iloc.in_direction(d) * 0.5 * inva);
  return loc;
}

/******************************************************/
/* global variables related to threads plus a routine */
/* for setting their values via environment variable  */
/******************************************************/
vector<ivec_loop_counter> ilcs;
int meep_threads=0;
bool use_stride1=true;

void init_meep_threads() {
  char *s=getenv("MEEP_NUM_THREADS");
  if (s)
    sscanf(s,"%i",&meep_threads);
  if (verbosity > 0)
    master_printf("Using %i OpenMP threads.\n",meep_threads);

  s=getenv("MEEP_IGNORE_STRIDE1");
  use_stride1 = !s || s[0] != '1';
  if (verbosity > 0)
    master_printf("%s #pragma ivdep for stride-1 loops.\n",use_stride1 ? "Using" : "Not using");
}

void update_ilcs(const grid_volume &gv) {
  if (meep_threads==0 || ((int)ilcs.size())==meep_threads)
    return;
  ilcs.clear();
  ivec is(gv.dim);
  for(int nt=0; nt<meep_threads; nt++)
    ilcs.push_back(ivec_loop_counter(gv,is,is,nt,meep_threads));
}

} // namespace meep
