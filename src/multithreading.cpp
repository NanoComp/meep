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
#include <complex>
#include <map>
#include <limits.h>

#include "multithreading.hpp"

using namespace std;

namespace meep {

/***************************************************************************/
/* simple benchmarking: call checkpoint(__FILE__,__LINE) upon entering a   */
/* section of code, then call checkpoint() upon leaving it. periodically   */
/* call print_stats() to dump results to file and/or console.              */
/***************************************************************************/
typedef std::pair<string,double> snippet_time;
typedef std::map<string,double> snippet_time_map;
static snippet_time_map snippet_times;
typedef std::map<string,int> snippet_count_map;
static snippet_count_map snippet_counts;
size_t checkpoint(const char *name, int line)
{ 
  static string current_snippet;
  static double current_start_time=0.0;

  if (!current_snippet.empty() && snippet_times.find(current_snippet)!=snippet_times.end())
   { snippet_times[current_snippet]  += wall_time() - current_start_time;
     snippet_counts[current_snippet] ++;
   }
  current_snippet.clear();

  if (name==0)
   return 0;

  current_snippet = string(name) + (line==0 ? "" : ":" + to_string(line));
  if (snippet_times.find(current_snippet)==snippet_times.end())
   { snippet_times[current_snippet]=0.0;
     snippet_counts[current_snippet]=0;
   }
  current_start_time=wall_time();
  return 0;
}

void print_benchmarks()
{ for(snippet_time_map::iterator it=snippet_times.begin(); it!=snippet_times.end(); ++it)
   printf("%30s: %4ix  %.6es total   %.3es avg\n", it->first.c_str(),
           snippet_counts[it->first], it->second,it->second/snippet_counts[it->first]);
}

/***************************************************************/
/* implementation of ivec_loop_counter                         */
/***************************************************************/
ivec_loop_counter::ivec_loop_counter()
 { init(grid_volume(), ivec(), ivec()); }

ivec_loop_counter::ivec_loop_counter(grid_volume gv, ivec _is, ivec _ie, int nt, int NT)
{ init(gv, _is, _ie, nt, NT); } 

int ivec_loop_counter::init(grid_volume gv, ivec _is, ivec _ie, int nt, int NT)
{ 
  if (is==_is && ie==_ie && gvlc==gv.little_corner() && inva==gv.inva)
   return 0; // no need to reinitialize

  is   = _is;
  gvlc = gv.little_corner();
  ie   = _ie;
  inva = gv.inva;

  // tabulate 'active' dimensions, i.e. those in which loop is non-empty
  idx0 = 0;
  rank = 0;
  size_t num_iters = 1;
  LOOP_OVER_DIRECTIONS(gv.dim, d)
   { idx0 += gv.stride(d) * ((is - gvlc).in_direction(d)/2);
     ptrdiff_t count = (ie-is).in_direction(d)/2 + 1;
     if ( count > 1 )
      { idx_stride[rank] =  gv.stride(d);
        loop_dir[rank]   =  d;
        N[rank]          =  count;
        num_iters        *= count;
        rank++;
      }
   }
  N_inner = N[rank-1];

  if (NT<1) {nt=0; NT=1;}
  min_iter = (nt*num_iters) / NT;
  max_iter = ((nt+1)*num_iters) / NT;
  return 0;
}

/***************************************************************/
/* helper function for the 'start' routine that initializes    */
/* computation of 'k' index into PML sigma arrays for up to    */
/* two distinct directions                                     */
/***************************************************************/
void ivec_loop_counter::init_k(int nd, direction dsig)
{ 
  if (dsig==NO_DIRECTION)
   k0[nd]=INT_MAX;
  else
   { k0[nd]=(is-gvlc).in_direction(dsig);
     for(int r=0; r<rank; r++)
      k_stride[nd][r] = (loop_dir[r]==dsig ? 2 : 0);
   }
}

ptrdiff_t ivec_loop_counter::start(direction dsig1, int *k1, direction dsig2, int *k2)
{
  complete = false;
  N_inner=N[rank-1];
  current_iter=min_iter;
  size_t idx=niter_to_narray(current_iter, n);

  init_k(0, dsig1);
  init_k(1, dsig2);
  get_k(k1,k2);

  return idx;
}

ptrdiff_t ivec_loop_counter::start(size_t *k1, direction dsig1, size_t *k2, direction dsig2)
{ int ik1, ik2;
  ptrdiff_t idx = start(dsig1, k1 ? &ik1 : 0, dsig2, k2 ? &ik2 : 0);
  if (k1) *k1=ik1; 
  if (k2) *k2=ik2; 
  return idx;
}

void ivec_loop_counter::advance()
{ if (++current_iter < max_iter)
   for(int r=rank-1; r>=0; n[r--]=0)
    if ( ++n[r] < N[r] )
     return;
  complete=true;
  N_inner=0; // to break out of inner loop
}

ptrdiff_t ivec_loop_counter::increment(int *k1, int *k2)
{ advance();
  get_k(k1,k2);
  return get_idx();
}

ptrdiff_t ivec_loop_counter::increment(size_t *k1, size_t *k2)
{ int ik1, ik2;
  ptrdiff_t idx=increment( k1 ? &ik1 : 0, k2 ? &ik2 : 0 );
  if (k1) *k1=ik1;
  if (k2) *k2=ik2;
  return idx;
}

ptrdiff_t ivec_loop_counter::niter_to_narray(size_t niter, size_t narray[3])
{ for(int r=rank-1; r>=0; niter/=N[r--])
   narray[r] = niter%N[r];
  return get_idx(narray);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
ptrdiff_t ivec_loop_counter::get_idx(size_t *narray)
{ if (narray==0) narray=n;
  ptrdiff_t idx = idx0;
  for(int r=0; r<rank; r++)
   idx += narray[r]*idx_stride[r];
  return idx;
}

int ivec_loop_counter::get_k(int nd, size_t *narray)
{ if (narray==0) narray=n;
  if (k0[nd]==INT_MAX) return 0;
  int k=k0[nd];
  for(int r=0; r<rank; r++)
   k+=narray[r]*k_stride[nd][r];
  return k;
}

void ivec_loop_counter::get_k(int *k1, int *k2, size_t *narray)
{ if (k1) *k1=get_k(0, narray);
  if (k2) *k2=get_k(1, narray);
}

ivec ivec_loop_counter::get_iloc(size_t *narray)
{ if (narray==0) narray=n;
  ivec iloc=is;
  for(int r=0; r<rank; r++)
   iloc.set_direction(loop_dir[r], is.in_direction(loop_dir[r]) + 2*n[r]);
  return iloc;
}

vec ivec_loop_counter::get_loc(size_t *narray)
{ if(narray==0) narray=n;
  ivec iloc=get_iloc();
  vec loc(iloc.dim);
  LOOP_OVER_DIRECTIONS(iloc.dim,d)
   loc.set_direction(d, iloc.in_direction(d) * 0.5 * inva);
  return loc;
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
vector<ivec_loop_counter> ilc;
int meep_threads=0;
int thread_strategy=0;
void set_meep_threads(int num_threads, int new_strategy)
{ 
  meep_threads=num_threads;
  if (meep_threads != ilc.size()) 
   { ilc.clear();
     while(num_threads--)
      ilc.push_back(ivec_loop_counter());
printf("meep_threads=%i, ilc.size=%i\n",meep_threads, ilc.size());
   }
  if (new_strategy!=-1) thread_strategy=new_strategy;
}

} // namespace meep
