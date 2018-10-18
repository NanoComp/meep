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

#ifndef MULTITHREADING_H
#define MULTITHREADING_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex>

#include "meep.hpp"

using namespace std;
namespace meep { 

/**************************************************/
/* loop execution benchmarking                    */
/**************************************************/
void checkpoint(const char *code=0, int line=0);
void print_benchmarks(double meep_time, char *FileName=0);

#define BENCHMARK

#ifdef BENCHMARK
  #define CHECKPOINT(a,b) checkpoint(a,b);
#else
  #define CHECKPOINT(a,b)
#endif

// benchmarked versions of (non-multithreaded) grid-loop primitives
#define BLOOP_OVER_VOL_OWNED(gv, c, idx)  	\
  CHECKPOINT(__FILE__,__LINE__) 	  	\
  LOOP_OVER_VOL_OWNED(gv,c,idx)          

#define BLOOP_OVER_VOL_OWNED0(gv, c, idx) 	\
  CHECKPOINT(__FILE__,__LINE__);	  	\
  LOOP_OVER_VOL_OWNED0(gv,c,idx)

#define BS1LOOP_OVER_VOL_OWNED(gv, c, idx)	\
  CHECKPOINT(__FILE__,__LINE__) 	  	\
  S1LOOP_OVER_VOL_OWNED(gv,c,idx)          

#define BS1LOOP_OVER_VOL_OWNED0(gv, c, idx)	\
  CHECKPOINT(__FILE__,__LINE__);	  	\
  S1LOOP_OVER_VOL_OWNED0(gv,c,idx)


/**************************************************/
/* ivec_loop_counter is a class that replaces the */
/* LOOP_OVER_IVECS() macro and its derivatives to */
/* facilitate shared-memory loop parallelization. */
/**************************************************/
class ivec_loop_counter
 {
public:

    /***************************************************************/
    /* initialize to handle loops from _is to _ie in gv.           */
    /* if NT>1 and nt \in [0,NT), the loop is subdivided into NT   */
    /* equal chunks and the loop counter handles only chunk #nt.   */
    /***************************************************************/
    ivec_loop_counter(const grid_volume &gv, const ivec &_is, const ivec &_ie, int nt=0, int NT=1);
    void init(const grid_volume &gv, const ivec &_is, const ivec &_ie, int nt=0, int NT=1);

    /*--------------------------------------------------------------*/
    /* The start() and increment() routines implement a stateful    */
    /* paradigm in which the current state of the loop is stored in */
    /* internal class variables. In this case we have separate      */
    /* ivec_loop_counter structures for each thread, each covering  */
    /* a fraction of the total loop.                                */
    /*--------------------------------------------------------------*/

    /***************************************************************/
    /* Set the loop counter at the beginning of our chunk of the   */
    /* loop. The return value is 'idx' (index into field arrays) at*/
    /* the starting grid point. The optional direction arguments   */
    /* request computation of 'k' indices into PML sigma arrays for*/
    /* up to two directions. If k1,k2 are non-null they are filled */
    /* in with the corresponding k values at the starting grid pt. */
    /***************************************************************/
    ptrdiff_t start(direction dsig1=NO_DIRECTION, direction dsig2=NO_DIRECTION);
    ptrdiff_t start(const grid_volume &gv, const ivec &_is, const ivec &_ie, int nt, int NT,
                    direction dsig1=NO_DIRECTION, direction dsig2=NO_DIRECTION);

    ptrdiff_t update_outer();

    /***************************************************************/
    /* Advance the loop one iteration. The return value is idx.    */
    /* The internal class field 'complete' will be set to true if  */
    /* this was the final loop iteration.                          */
    /***************************************************************/
    ptrdiff_t increment(int *k1=0, int *k2=0);
    ptrdiff_t increment(size_t *k1, size_t *k2=0);
    ptrdiff_t operator++() { return increment(); }

    /*--------------------------------------------------------------*/
    /*- The niter_to_narray() routine implements an alternative,    */
    /*- fully stateless, model. In this case we will have only one  */
    /*- ivec_loop_counter structure used by *all* threads with no   */
    /*- thread-specific data stored in the structure. Threads call  */
    /*- niter_to_narray() with arbitrary values of niter in the     */
    /*- range [min_iter, max_iter) to compute narray, the loop      */
    /*- indices for that iteration (the return value is idx). Then  */
    /*- narray can be passed to routines like get_k, get_iloc, etc. */
    /*- to compute loop quantities for the given iteration.         */
    /*--------------------------------------------------------------*/
    ptrdiff_t niter_to_narray(size_t niter, size_t narray[3]);

    /***************************************************************/
    /* The following routines are used in both the state-ful model */
    /* (in which case they should be called with no arguments, to  */
    /* indicate that internal class variables should be used) and  */
    /* in the stateless approach, in which case the caller should  */
    /* supply an narray argument filled in by niter_to_narray.     */
    /***************************************************************/
    ptrdiff_t get_idx(size_t *narray=0);
    int get_k(int nd, size_t *narray=0);
    void get_k(int *k1, int *k2=0, size_t *narray=0);
    ivec get_iloc(size_t *narray=0);
    vec get_loc(size_t *narray=0);

// private

// routines intended for internal use
    void init_k(int nd, direction dsig);   // helper routine for start()
    void advance();                        // helper routine for increment()

// class data fields

    // fields used for both state-ful and stateless usage models
    ivec is, ie, gvlc; 
    double inva;
    direction loop_dir[3];
    ptrdiff_t idx_stride[3];
    size_t N[3];
    int rank;
    ptrdiff_t idx0;
    size_t min_iter, max_iter;
    int k0[2], k_stride[2][3]; 

    // fields used only for the state-ful model
    size_t n[3], N_inner;
    size_t current_iter;
    size_t idx_max, idx_step, k_step[2];
    bool complete;

 }; // class ivec_loop_counter

/**************************************************/
/* global array of ivec_loop_counters, one per    */
/* thread, initialized by set_meep_threads()      */
/**************************************************/
extern vector<ivec_loop_counter> ilcs;
extern int meep_threads;
extern int thread_strategy;
extern bool use_stride1;
void set_meep_threads(int num_threads=0, int new_strategy=-1);
void update_ilcs(const grid_volume &gv);

/**************************************************/
/* macros for multithreaded grid loops. The basic */
/* grid loop is PLOOP_OVER_IVECS. Versions with   */
/* _D and _DD suffixes are cases requiring one    */
/* or two 'k' indices into PML sigma arrays of    */
/* specific dimentions. S1PLOOP_ versions break   */
/* out the innermost loop dimension as a separate */
/* loop to allow compiler optimization using      */
/* #pragma ivdep (intel compilers) or             */
/* #pragma gcc ivdep (GNU compilers).             */
/**************************************************/
#define _NT meep_threads

#define PLOOP_OVER_IVECS(gv, is, ie, idx) 							\
 CHECKPOINT(__FILE__,__LINE__)									\
_Pragma("omp parallel for schedule(guided), num_threads(_NT)")					\
 for(int nt=0; nt<_NT; nt++)									\
  for(size_t idx=ilcs[nt].start(gv,is,ie,nt,_NT); !(ilcs[nt].complete); idx=++ilcs[nt])

#define PLOOP_OVER_IVECS_D(gv, is, ie, idx, d, k) 					\
 CHECKPOINT(__FILE__,__LINE__)								\
_Pragma("omp parallel for schedule(guided), num_threads(_NT)")				\
 for(int nt=0; nt<_NT; nt++)								\
  for(size_t idx=ilcs[nt].start(gv,is,ie,nt,_NT,d), k=ilcs[nt].get_k(0); !(ilcs[nt].complete); idx=ilcs[nt].increment(&k))

#define PLOOP_OVER_IVECS_DD(gv, is, ie, idx, d1, k1, d2, k2) 				\
 CHECKPOINT(__FILE__,__LINE__)								\
_Pragma("omp parallel for schedule(guided), num_threads(_NT)")				\
 for(int nt=0; nt<_NT; nt++)								\
  for(size_t idx=ilcs[nt].start(gv,is,ie,nt,_NT,d1,d2), k1=ilcs[nt].get_k(0), k2=ilcs[nt].get_k(1); !(ilcs[nt].complete); idx=ilcs[nt].increment(&k1,&k2))

#define S1PLOOP_OVER_IVECS(gv, is, ie, idx)						\
 CHECKPOINT(__FILE__,__LINE__)								\
_Pragma("omp parallel for schedule(guided), num_threads(_NT)")				\
 for(int nt=0; nt<_NT; nt++)								\
  for(size_t idx=ilcs[nt].start(gv,is,ie,nt,_NT); !(ilcs[nt].complete); idx=ilcs[nt].update_outer()) \
_Pragma(IVDEP)										\
   for(; idx<ilcs[nt].idx_max; idx+=ilcs[nt].idx_step)

#define S1PLOOP_OVER_IVECS_D(gv, is, ie, idx, d, k)					\
 CHECKPOINT(__FILE__,__LINE__)								\
_Pragma("omp parallel for schedule(guided), num_threads(_NT)")				\
 for(int nt=0; nt<_NT; nt++)								\
  for(size_t idx=ilcs[nt].start(gv,is,ie,nt,_NT, d); !(ilcs[nt].complete); idx=ilcs[nt].update_outer()) \
_Pragma(IVDEP)										\
   for(int k=ilcs[nt].get_k(0); idx<ilcs[nt].idx_max; idx+=ilcs[nt].idx_step, k+=ilcs[nt].k_step[0])

#define S1PLOOP_OVER_IVECS_DD(gv, is, ie, idx, d1, k1, d2, k2)				\
 CHECKPOINT(__FILE__,__LINE__)								\
_Pragma("omp parallel for schedule(guided), num_threads(_NT)")				\
 for(int nt=0; nt<_NT; nt++)								\
  for(size_t idx=ilcs[nt].start(gv,is,ie,nt,_NT, d1, d2); !(ilcs[nt].complete); idx=ilcs[nt].update_outer()) \
_Pragma(IVDEP)										\
   for(int k1=ilcs[nt].get_k(0), k2=ilcs[nt].get_k(1); idx<ilcs[nt].idx_max; idx+=ilcs[nt].idx_step, k1+=ilcs[nt].k_step[0], k2+=ilcs[nt].k_step[1])

#if 0
#define S1PLOOP_OVER_IVECS(gv, is, ie, idx)						\
 CHECKPOINT(__FILE__,__LINE__)								\
_Pragma("omp parallel for schedule(guided), num_threads(_NT)")				\
 for(int nt=0; nt<_NT; nt++)								\
  for(size_t idx=ilcs[nt].start(gv,is,ie,nt,_NT); !(ilcs[nt].complete); /*empty*/ )	\
_Pragma(IVDEP)										\
   for(size_t n_inner=0; n_inner<ilcs[nt].N_inner; n_inner++, idx=++ilcs[nt])

#define S1PLOOP_OVER_IVECS_D(gv, is, ie, idx, d, k)						\
 CHECKPOINT(__FILE__,__LINE__)									\
_Pragma("omp parallel for schedule(guided), num_threads(_NT)")					\
 for(int nt=0; nt<_NT; nt++)									\
  for(size_t k=0, idx=ilcs[nt].start(gv,is,ie,nt,_NT, d, &k); !(ilcs[nt].complete); /*empty*/) 	\
_Pragma(IVDEP)											\
   for(size_t n_inner=0; n_inner<ilcs[nt].N_inner; n_inner++, idx=ilcs[nt].incremen(&k))

#define S1PLOOP_OVER_IVECS_DD(gv, is, ie, idx, d1, k1, d2, k2)					\
 CHECKPOINT(__FILE__,__LINE__)									\
_Pragma("omp parallel for schedule(guided), num_threads(_NT)")					\
 for(int nt=0; nt<_NT; nt++)									\
  for(size_t k1=0, k2=0, idx=ilcs[nt].start(gv,is,ie,nt,_NT, d1, &k1, d2, &k2); !(ilcs[nt].complete); /*empty*/) \
_Pragma(IVDEP)											\
   for(size_t n_inner=0; n_inner<ilcs[nt].N_inner; n_inner++, idx=ilcs[nt].increment(&k1,&k2))
#endif

#define IVEC_PLOOP_LOC(loc)		\
 vec loc=ilcs[nt].get_loc();

#define IVEC_PLOOP_ILOC(iloc)		\
 ivec iloc=ilcs[nt].get_iloc();

} // namespace meep

#endif // MULTITHREADING_H
