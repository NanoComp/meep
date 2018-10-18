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

#include <utility>
#include <map>
#include <vector>
#include <algorithm>

#include "multithreading.hpp"

using namespace std;

namespace meep {

/***************************************************************************/
/* simple benchmarking: call checkpoint(__FILE__,__LINE) upon entering a   */
/* section of code, then call checkpoint() upon leaving it. periodically   */
/* call print_benchmarks() to dump results to file and/or console.         */
/***************************************************************************/
std::map<string,double> snippet_times;
std::map<string,int> snippet_counts;
double t0=HUGE_VAL;
double benchmark_interval=5.0;
double next_benchmark_time=benchmark_interval;
bool checkpoint_started=false;
void checkpoint(const char *name, int line)
{ 
  if (meep_steps < 5) return;
  static string current_snippet;
  static double current_start_time=0.0;

  if (name==0 && line==-1)
   { snippet_times.clear();
     snippet_counts.clear();
     current_snippet.clear();
     t0=HUGE_VAL;
     benchmark_interval=5.0;
     next_benchmark_time=benchmark_interval;
     return;
   }
 
  if (t0==HUGE_VAL) t0=wall_time();

  if (!current_snippet.empty() && snippet_times.count(current_snippet) )
   { snippet_times[current_snippet]  += wall_time() - current_start_time;
     snippet_counts[current_snippet] ++;
   }
  current_snippet.clear();

  if (name==0) return;

  if (!strcmp(name,"step_generic.cpp"))
   current_snippet=string("s_g.cpp");
  else if (!strcmp(name,"step_generic_stride1.cpp"))
   current_snippet=string("s_g_s1.cpp");
  else if (!strcmp(name,"step_multithreaded.cpp"))
   current_snippet=string("s_m.cpp");
  else if (!strcmp(name,"step_multithreaded_stride1.cpp"))
   current_snippet=string("s_m_s1.cpp");
  else 
   current_snippet = string(name);
  if (line!=0) current_snippet += ":" + to_string(line);

  if (snippet_times.count(current_snippet)==0)
   { snippet_times[current_snippet]=0.0;
     snippet_counts[current_snippet]=0;
   }
  current_start_time=wall_time();
}

void print_benchmarks(double meep_time, char *FileName)
{
#ifndef BENCHMARK
  return;
#endif
  if (meep_time<next_benchmark_time) return;
  next_benchmark_time += benchmark_interval;

  double ttot=wall_time()-t0;
  FILE *f = FileName ? fopen(FileName,"a") : stdout;
  fprintf(f,"\n\nNT %i   STEPS %i   WALL(TOT) %e  WALL(PER) %e\n",meep_threads,meep_steps,ttot,ttot/meep_steps);
  for(auto it=snippet_times.begin(); it!=snippet_times.end(); ++it)
   fprintf(f,"NT %i %20s: %8ix  %.1e s tot (%.2f %%) %.3es avg\n",meep_threads,
            it->first.c_str(),snippet_counts[it->first],it->second,
            it->second/ttot,it->second/snippet_counts[it->first]);
  if (FileName) fclose(f);
}

/***************************************************************/
/* implementation of ivec_loop_counter                         */
/***************************************************************/
ivec_loop_counter::ivec_loop_counter(const grid_volume &gv, const ivec &_is, const ivec &_ie, int nt, int NT)
{ init(gv, _is, _ie, nt, NT); } 

void ivec_loop_counter::init(const grid_volume &gv, const ivec &_is, const ivec &_ie, int nt, int NT)
{ 
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
  N_inner  = N[rank-1];
  idx_step = idx_stride[rank-1];

  if (NT<1) {nt=0; NT=1;}
  min_iter = (nt*num_iters) / NT;
  max_iter = ((nt+1)*num_iters) / NT;
}

/***************************************************************/
/* helper function for the 'start' routine that initializes    */
/* computation of 'k' index into PML sigma arrays for up to    */
/* two distinct directions                                     */
/***************************************************************/
void ivec_loop_counter::init_k(int nd, direction dsig)
{ 
  if (dsig==NO_DIRECTION)
   { k0[nd]=INT_MAX;
     k_step[nd]=0;
   }
  else
   { k0[nd]=(is-gvlc).in_direction(dsig);
     for(int r=0; r<rank; r++)
      k_stride[nd][r] = (loop_dir[r]==dsig ? 2 : 0);
     k_step[nd]=k_stride[nd][rank-1];
   }
}

ptrdiff_t ivec_loop_counter::start(direction dsig1, direction dsig2)
{
  complete = false;
  init_k(0, dsig1);
  init_k(1, dsig2);

  current_iter=min_iter;
  N_inner=0;
  return update_outer();
}

ptrdiff_t ivec_loop_counter::start(const grid_volume &gv, const ivec &_is, const ivec &_ie, int nt, int NT,
                                   direction dsig1, direction dsig2)
{ init(gv, _is, _ie, nt, NT);
  return start(dsig1, dsig2);
}

ptrdiff_t ivec_loop_counter::update_outer()
{ 
  current_iter += N_inner;
  if (current_iter>=max_iter)
   { complete=true; return idx0; }
  ptrdiff_t idx = niter_to_narray(current_iter, n);
  N_inner = N[rank-1] - n[rank-1];
  if (current_iter + N_inner > max_iter)
   N_inner = max_iter - current_iter;
  idx_max = idx + N_inner*idx_step;
  return idx;
}

void ivec_loop_counter::advance()
{ if (++current_iter < max_iter)
   for(int r=rank-1; r>=0; n[r--]=0)
    if ( ++n[r] < N[r] )
     return;
  complete=true;
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
{ ptrdiff_t idx=idx0;
  for(int r=rank-1; r>0; niter/=N[r--])
   { narray[r] = niter % N[r];
     idx+=narray[r]*idx_stride[r];
   }
  narray[0]=niter;
  return idx+narray[0]*idx_stride[0];
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
{ if (k0[nd]==INT_MAX) return 0;
  if (narray==0) narray=n;
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
   iloc.set_direction(loop_dir[r], is.in_direction(loop_dir[r]) + 2*narray[r]);
  return iloc;
}

vec ivec_loop_counter::get_loc(size_t *narray)
{ ivec iloc=get_iloc(narray);
  vec loc(iloc.dim);
  LOOP_OVER_DIRECTIONS(iloc.dim,d)
   loc.set_direction(d, iloc.in_direction(d) * 0.5 * inva);
  return loc;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
vector<ivec_loop_counter> ilcs;
int meep_threads=0;
int thread_strategy=0;
bool use_stride1=false;
int meep_steps=0;
void set_meep_threads(int num_threads, int new_strategy)
{ 
  meep_threads=num_threads;
  if (new_strategy!=-1) thread_strategy=new_strategy;

  char *s=getenv("MEEP_USE_STRIDE1");
  use_stride1=(s && s[0]=='1');
  printf("%s #pragma ivdep for stride-1 loops.\n",use_stride1 ? "Using" : "Not using");

  s=getenv("MEEP_BENCHMARK_INTERVAL");
  if (s)
   { sscanf(s,"%le",&benchmark_interval);
     printf("Benchmarking at intervals of %g.\n",benchmark_interval);
   }
}

void update_ilcs(const grid_volume &gv)
{
  if (meep_threads==0 || ((int)ilcs.size())==meep_threads)
   return;
  ilcs.clear();
  ivec is(gv.dim);
  for(int nt=0; nt<meep_threads; nt++)
   ilcs.push_back(ivec_loop_counter(gv,is,is,nt,meep_threads));
}

} // namespace meep
