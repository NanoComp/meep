/* Copyright (C) 2005-2009 Massachusetts Institute of Technology
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

#include <stdarg.h>
#include <string.h>
#include <stdlib.h>

#include "meep.hpp"
#include "config.h"

#ifdef HAVE_MPI
// undef'ing SEEK_* is needed for MPICH, possibly other MPI versions
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include <mpi.h>
#endif

#ifdef IGNORE_SIGFPE
#  include <signal.h>
#endif

#if defined(DEBUG) && defined(HAVE_FEENABLEEXCEPT)
#  ifndef _GNU_SOURCE
#    define _GNU_SOURCE 1
#  endif
#  include <fenv.h>
#  if !HAVE_DECL_FEENABLEEXCEPT
extern "C" int feenableexcept (int EXCEPTS);
#  endif
#endif

#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif
#ifdef HAVE_BSDGETTIMEOFDAY
#  ifndef HAVE_GETTIMEOFDAY
#    define gettimeofday BSDgettimeofday
#    define HAVE_GETTIMEOFDAY 1
#  endif
#endif

#define UNUSED(x) (void) x // silence compiler warnings

#define MPI_REALNUM (sizeof(realnum) == sizeof(double) ? MPI_DOUBLE:MPI_FLOAT)

namespace meep {

#ifdef HAVE_MPI
  static MPI_Comm mycomm = MPI_COMM_WORLD;
#endif

bool quiet = false; // defined in meep.h

initialize::initialize(int &argc, char** &argv) {
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  int major, minor;
  MPI_Get_version(&major, &minor);
  if (!quiet) master_printf("Using MPI version %d.%d, %d processes\n", 
			    major, minor, count_processors());
#else
  UNUSED(argc);
  UNUSED(argv);
#endif
#if defined(DEBUG_FP) && defined(HAVE_FEENABLEEXCEPT)
  feenableexcept(FE_INVALID | FE_OVERFLOW); //crash if NaN created, or overflow
#endif
#ifdef IGNORE_SIGFPE
  signal(SIGFPE, SIG_IGN);
#endif
  t_start = wall_time();
}

initialize::~initialize() {
  if (!quiet) master_printf("\nElapsed run time = %g s\n", elapsed_time());
#ifdef HAVE_MPI
  end_divide_parallel();
  MPI_Finalize();
#endif
}

double wall_time(void) {
#ifdef HAVE_MPI
  return MPI_Wtime();
#elif HAVE_GETTIMEOFDAY
  struct timeval tv;
  gettimeofday(&tv, 0);
  return(tv.tv_sec + tv.tv_usec * 1e-6);
#else
  return (clock() * 1.0 / CLOCKS_PER_SECOND);
#endif
}

void abort(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  fprintf(stderr, "meep: ");
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  if (fmt[strlen(fmt) - 1] != '\n') fputc('\n', stderr); // force newline
#ifdef HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  exit(1);
}

void send(int from, int to, double *data, int size) {
#ifdef HAVE_MPI
  if (from == to) return;
  if (size == 0) return;
  const int me = my_rank();
  if (from == me) MPI_Send(data, size, MPI_DOUBLE, to, 1, mycomm);
  MPI_Status stat;
  if (to == me) MPI_Recv(data, size, MPI_DOUBLE, from, 1, mycomm, &stat);
#else
  UNUSED(from);
  UNUSED(to);
  UNUSED(data);
  UNUSED(size);
#endif
}

#if MEEP_SINGLE
void broadcast(int from, realnum *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, size, MPI_FLOAT, from, mycomm);
#else
  UNUSED(from);
  UNUSED(data);
  UNUSED(size);
#endif
}
#endif

void broadcast(int from, double *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, size, MPI_DOUBLE, from, mycomm);
#else
  UNUSED(from);
  UNUSED(data);
  UNUSED(size);
#endif
}

void broadcast(int from, char *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, size, MPI_CHAR, from, mycomm);
#else
  UNUSED(from);
  UNUSED(data);
  UNUSED(size);
#endif
}

void broadcast(int from, complex<double> *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, 2*size, MPI_DOUBLE, from, mycomm);
#else
  UNUSED(from);
  UNUSED(data);
  UNUSED(size);
#endif
}

void broadcast(int from, int *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, size, MPI_INT, from, mycomm);
#else
  UNUSED(from);
  UNUSED(data);
  UNUSED(size);
#endif
}

complex<double> broadcast(int from, complex<double> data) {
#ifdef HAVE_MPI
  MPI_Bcast(&data, 2, MPI_DOUBLE, from, mycomm);
#else
  UNUSED(from);
#endif
  return data;
}

double broadcast(int from, double data) {
#ifdef HAVE_MPI
  MPI_Bcast(&data, 1, MPI_DOUBLE, from, mycomm);
#else
  UNUSED(from);
#endif
  return data;
}

int broadcast(int from, int data) {
#ifdef HAVE_MPI
  MPI_Bcast(&data, 1, MPI_INT, from, mycomm);
#else
  UNUSED(from);
#endif
  return data;
}

bool broadcast(int from, bool b) {
  return broadcast(from, (int) b);
}

double max_to_master(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Reduce(&in,&out,1,MPI_DOUBLE,MPI_MAX,0,mycomm);
#endif
  return out;
}

double max_to_all(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,1,MPI_DOUBLE,MPI_MAX,mycomm);
#endif
  return out;
}

int max_to_all(int in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,1,MPI_INT,MPI_MAX,mycomm);
#endif
  return out;
}

ivec max_to_all(const ivec &v) {
  int in[5], out[5];
  for (int i=0; i<5; ++i) in[i] = out[i] = v.in_direction(direction(i));
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,5,MPI_INT,MPI_MAX,mycomm);
#endif
  ivec vout(v.dim);
  for (int i=0; i<5; ++i) vout.set_direction(direction(i), out[i]);
  return vout;
}

double sum_to_master(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Reduce(&in,&out,1,MPI_DOUBLE,MPI_SUM,0,mycomm);
#endif
  return out;
}

double sum_to_all(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,1,MPI_DOUBLE,MPI_SUM,mycomm);
#endif
  return out;
}

void sum_to_all(const double *in, double *out, int size) {
#ifdef HAVE_MPI
  MPI_Allreduce((void*) in, out, size, MPI_DOUBLE,MPI_SUM,mycomm);
#else
  memcpy(out, in, sizeof(double) * size);
#endif
}

long double sum_to_all(long double in) {
  long double out = in;
#ifdef HAVE_MPI
  if (MPI_LONG_DOUBLE == MPI_DATATYPE_NULL)
    out = sum_to_all(double(in));
  else
    MPI_Allreduce(&in,&out,1,MPI_LONG_DOUBLE,MPI_SUM,mycomm);
#endif
  return out;
}

int sum_to_all(int in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,1,MPI_INT,MPI_SUM,mycomm);
#endif
  return out;
}

int partial_sum_to_all(int in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Scan(&in,&out,1,MPI_INT,MPI_SUM,mycomm);
#endif
  return out;
}

complex<double> sum_to_all(complex<double> in) {
  complex<double> out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,2,MPI_DOUBLE,MPI_SUM,mycomm);
#endif
  return out;
}

complex<long double> sum_to_all(complex<long double> in) {
  complex<long double> out = in;
#ifdef HAVE_MPI
  if (MPI_LONG_DOUBLE == MPI_DATATYPE_NULL) {
    complex<double> dout;
    dout = sum_to_all(complex<double>(double(in.real()), double(in.imag())));
    out = complex<long double>(dout.real(), dout.imag());
  }
  else
    MPI_Allreduce(&in,&out,2,MPI_LONG_DOUBLE,MPI_SUM,mycomm);
#endif
  return out;
}

bool or_to_all(bool in) {
  int in2 = in, out;
#ifdef HAVE_MPI
  MPI_Allreduce(&in2,&out,1,MPI_INT,MPI_LOR,mycomm);
#else
  out = in2;
#endif
  return (bool) out;
}

void or_to_all(const int *in, int *out, int size) {
#ifdef HAVE_MPI
  MPI_Allreduce((void*) in, out, size, MPI_INT,MPI_LOR,mycomm);
#else
  memcpy(out, in, sizeof(int) * size);
#endif
}

bool and_to_all(bool in) {
  int in2 = in, out;
#ifdef HAVE_MPI
  MPI_Allreduce(&in2,&out,1,MPI_INT,MPI_LAND,mycomm);
#else
  out = in2;
#endif
  return (bool) out;
}

void and_to_all(const int *in, int *out, int size) {
#ifdef HAVE_MPI
  MPI_Allreduce((void*) in, out, size, MPI_INT,MPI_LAND,mycomm);
#else
  memcpy(out, in, sizeof(int) * size);
#endif
}

void all_wait() {
#ifdef HAVE_MPI
  MPI_Barrier(mycomm);
#endif
}

int my_rank() {
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(mycomm, &rank);
  return rank;
#else
  return 0;
#endif
}

int count_processors() {
#ifdef HAVE_MPI
  int n;
  MPI_Comm_size(mycomm, &n);
  return n;
#else
  return 1;
#endif
}

void fields::boundary_communications(field_type ft) {
  // Communicate the data around!
#if 0 // This is the blocking version, which should always be safe!
  for (int noti=0;noti<num_chunks;noti++)
    for (int j=0;j<num_chunks;j++) {
      const int i = (noti+j)%num_chunks;
      const int pair = j+i*num_chunks;
      DOCMP {
        send(chunks[j]->n_proc(), chunks[i]->n_proc(),
             comm_blocks[ft][pair], comm_size_tot(ft,pair));
      }
    }
#endif
#ifdef HAVE_MPI
  const int maxreq = num_chunks*num_chunks;
  MPI_Request *reqs = new MPI_Request[maxreq];
  MPI_Status *stats = new MPI_Status[maxreq];
  int reqnum = 0;
  int *tagto = new int[count_processors()];
  for (int i=0;i<count_processors();i++) tagto[i] = 0;
  for (int noti=0;noti<num_chunks;noti++)
    for (int j=0;j<num_chunks;j++) {
      const int i = (noti+j)%num_chunks;
      const int pair = j+i*num_chunks;
      const int comm_size = comm_size_tot(ft,pair);
      if (comm_size > 0) {
	if (chunks[j]->is_mine() && !chunks[i]->is_mine())
	  MPI_Isend(comm_blocks[ft][pair], comm_size,
		    MPI_REALNUM, chunks[i]->n_proc(),
		    tagto[chunks[i]->n_proc()]++,
		    mycomm, &reqs[reqnum++]);
	if (chunks[i]->is_mine() && !chunks[j]->is_mine())
	  MPI_Irecv(comm_blocks[ft][pair], comm_size,
		    MPI_REALNUM, chunks[j]->n_proc(),
		    tagto[chunks[j]->n_proc()]++,
		    mycomm, &reqs[reqnum++]);
      }
    }
  delete[] tagto;
  if (reqnum > maxreq) abort("Too many requests!!!\n");
  if (reqnum > 0) MPI_Waitall(reqnum, reqs, stats);
  delete[] reqs;
  delete[] stats;
#else
  (void) ft; // unused
#endif
}

// IO Routines...

bool am_really_master() {
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return (rank == 0);
#else
  return true;
#endif
}

void master_printf(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  if (am_really_master()) { vprintf(fmt, ap); fflush(stdout); }
  va_end(ap);
}

static FILE *debf = NULL;

void debug_printf(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  if (debf == NULL) {
    char temp[50];
    snprintf(temp, 50, "debug_out_%d", my_rank());
    debf = fopen(temp,"w");
    if (!debf) abort("Unable to open debug output %s\n", temp);
  }
  vfprintf(debf, fmt, ap);
  fflush(debf);
  va_end(ap);
}

void master_fprintf(FILE *f, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  if (am_master()) { vfprintf(f, fmt, ap); fflush(f); }
  va_end(ap);
}
FILE *master_fopen(const char *name, const char *mode) {
  FILE *f = am_master() ? fopen(name, mode) : 0;

  /* other processes need to know if fopen returned zero, in order
     to abort if fopen failed.  If fopen was successfully, just return
     a random non-zero pointer (which is never used except to compare to zero)
     on non-master processes */
  if (broadcast(0, bool(f != 0)) && !am_master())
    f = (FILE *) name;
  return f;
}
void master_fclose(FILE *f) {
  if (am_master()) fclose(f);
}

/* The following functions bracket a "critical section," a region
   of code that should be executed by only one process at a time.

   They work by having each process wait for a message from the
   previous process before starting. 

   Each critical section is passed an integer "tag"...ideally, this
   should be a unique identifier for each critical section so that
   messages from different critical sections don't get mixed up
   somehow. */

void begin_critical_section(int tag)
{
#ifdef HAVE_MPI
     int process_rank;
     MPI_Comm_rank(mycomm, &process_rank);
     if (process_rank > 0) { /* wait for a message before continuing */
	  MPI_Status status;
	  int recv_tag = tag - 1; /* initialize to wrong value */
	  MPI_Recv(&recv_tag, 1, MPI_INT, process_rank - 1, tag, 
		   mycomm, &status);
	  if (recv_tag != tag) abort("invalid tag received in begin_critical_section");
     }
#else
     UNUSED(tag);
#endif
}

void end_critical_section(int tag)
{
#ifdef HAVE_MPI
     int process_rank, num_procs;
     MPI_Comm_rank(mycomm, &process_rank);
     MPI_Comm_size(mycomm, &num_procs);
     if (process_rank != num_procs - 1) { /* send a message to next process */
	  MPI_Send(&tag, 1, MPI_INT, process_rank + 1, tag, 
		   mycomm);
     }
#else
     UNUSED(tag);
#endif
}


/* Simple, somewhat hackish API to allow user to run multiple simulations
   in parallel in the same MPI job.  The user calls

   mygroup = divide_parallel_processes(numgroups);

   to divide all of the MPI processes into numgroups equal groups,
   and to return the index (from 0 to numgroups-1) of the current group.
   From this point on, all fields etc. that you create and all
   calls from mympi.cpp will only communicate within your group of
   processes.

   However, there are two calls that you can use to switch back to
   globally communication among all processes:

   begin_global_communications();
   ....do stuff....
   end_global_communications();

   It is important not to mix the two types; e.g. you cannot timestep
   a field created in the local group in global mode, or vice versa.
*/

int divide_parallel_processes(int numgroups)
{
#ifdef HAVE_MPI
  end_divide_parallel();
  if (numgroups > count_processors()) abort("numgroups > count_processors");
  int mygroup = (my_rank() * numgroups) / count_processors();
  MPI_Comm_split(MPI_COMM_WORLD, mygroup, my_rank(), &mycomm);
  return mygroup;
#else
  if (numgroups != 1) abort("cannot divide processes in non-MPI mode");
  return 0;
#endif
}

#ifdef HAVE_MPI
  static MPI_Comm mycomm_save = MPI_COMM_WORLD;
#endif

void begin_global_communications(void)
{
#ifdef HAVE_MPI
  mycomm_save = mycomm;
  mycomm = MPI_COMM_WORLD;
#endif
}

void end_global_communications(void)
{
#ifdef HAVE_MPI
  mycomm = mycomm_save;
#endif
}

void end_divide_parallel(void)
{
#ifdef HAVE_MPI
  if (mycomm != MPI_COMM_WORLD) MPI_Comm_free(&mycomm);
  if (mycomm_save != MPI_COMM_WORLD) MPI_Comm_free(&mycomm_save);
  mycomm = mycomm_save = MPI_COMM_WORLD;
#endif
}

} // namespace meep
