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
    ivec_loop_counter(grid_volume gv, ivec _is, ivec _ie, int nt=0, int NT=1);
    ivec_loop_counter();
    int init(grid_volume gv, ivec _is, ivec _ie, int nt=0, int NT=1);

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
    ptrdiff_t start(direction dsig1=NO_DIRECTION, int *k1=0, direction dsig2=NO_DIRECTION, int *k2=0);
    ptrdiff_t start(size_t *k1, direction dsig1, size_t *k2=0, direction dsig2=NO_DIRECTION);
    ptrdiff_t start(direction dsig1, direction dsig2) { return start(dsig1,0,dsig2,0); }

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
    bool complete;

 }; // class ivec_loop_counter

/**************************************************/
/* global array of ivec_loop_counters, one per    */
/* thread, initialized by set_meep_threads()      */
/**************************************************/
extern vector<ivec_loop_counter> ilc;
extern int meep_threads;
extern int thread_strategy;
void set_meep_threads(int num_threads=0, int new_strategy=-1);

/**************************************************/
/* loop execution benchmarking                    */
/**************************************************/
size_t checkpoint(const char *code=0, int line=0);
void print_benchmarks();

#define BENCHMARK

#ifdef BENCHMARK
  #define CHECKPOINT(a,b) checkpoint(a,b);
#else
  #define CHECKPOINT(a,b)
#endif

// benchmarked versions of grid-loop primitives
#define BLOOP_OVER_VOL_OWNED(gv, c, idx)  \
  CHECKPOINT(__FILE__,__LINE__) 	  \
  LOOP_OVER_VOL_OWNED(gv,c,idx)          

#define BLOOP_OVER_VOL_OWNED0(gv, c, idx) \
  CHECKPOINT(__FILE__,__LINE__);	  \
  LOOP_OVER_VOL_OWNED0(gv,c,idx)

#define BS1LOOP_OVER_VOL_OWNED(gv, c, idx)  \
  CHECKPOINT(__FILE__,__LINE__) 	  \
  S1LOOP_OVER_VOL_OWNED(gv,c,idx)          

#define BS1LOOP_OVER_VOL_OWNED0(gv, c, idx) \
  CHECKPOINT(__FILE__,__LINE__);	  \
  S1LOOP_OVER_VOL_OWNED0(gv,c,idx)

} // namespace meep

#endif // MULTITHREADING_H
