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

#include <stdarg.h>
#include <string.h>

#include "meep.h"
#include "config.h"

#ifdef HAVE_MPI
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

namespace meep {

initialize::initialize(int argc, char **argv) {
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  int major, minor;
  MPI_Get_version(&major, &minor);
  master_printf("Using MPI... version %d.%d\n", major, minor);
#else
  UNUSED(argc);
  UNUSED(argv);
#endif
#if defined(DEBUG) && defined(HAVE_FEENABLEEXCEPT)
  feenableexcept(FE_INVALID | FE_OVERFLOW); //crash if NaN created, or overflow
#endif
#ifdef IGNORE_SIGFPE
  signal(SIGFPE, SIG_IGN);
#endif
  t_start = wall_time();
}

initialize::~initialize() {
  master_printf("Elapsed run time = %g s\n", elapsed_time());
#ifdef HAVE_MPI
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
  if (from == me) MPI_Send(data, size, MPI_DOUBLE, to, 1, MPI_COMM_WORLD);
  MPI_Status stat;
  if (to == me) MPI_Recv(data, size, MPI_DOUBLE, from, 1, MPI_COMM_WORLD, &stat);
#else
  UNUSED(from);
  UNUSED(to);
  UNUSED(data);
  UNUSED(size);
#endif
}

void broadcast(int from, double *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, size, MPI_DOUBLE, from, MPI_COMM_WORLD);
#else
  UNUSED(from);
  UNUSED(data);
  UNUSED(size);
#endif
}

void broadcast(int from, char *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, size, MPI_CHAR, from, MPI_COMM_WORLD);
#else
  UNUSED(from);
  UNUSED(data);
  UNUSED(size);
#endif
}

void broadcast(int from, complex<double> *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, 2*size, MPI_DOUBLE, from, MPI_COMM_WORLD);
#else
  UNUSED(from);
  UNUSED(data);
  UNUSED(size);
#endif
}

void broadcast(int from, int *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, size, MPI_INT, from, MPI_COMM_WORLD);
#else
  UNUSED(from);
  UNUSED(data);
  UNUSED(size);
#endif
}

complex<double> broadcast(int from, complex<double> data) {
#ifdef HAVE_MPI
  MPI_Bcast(&data, 2, MPI_DOUBLE, from, MPI_COMM_WORLD);
#else
  UNUSED(from);
#endif
  return data;
}

double broadcast(int from, double data) {
#ifdef HAVE_MPI
  MPI_Bcast(&data, 1, MPI_DOUBLE, from, MPI_COMM_WORLD);
#else
  UNUSED(from);
#endif
  return data;
}

int broadcast(int from, int data) {
#ifdef HAVE_MPI
  MPI_Bcast(&data, 1, MPI_INT, from, MPI_COMM_WORLD);
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
  MPI_Reduce(&in,&out,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
#endif
  return out;
}

double max_to_all(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  return out;
}

int max_to_all(int in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
  return out;
}

ivec max_to_all(const ivec &v) {
  int in[5], out[5];
  for (int i=0; i<5; ++i) in[i] = out[i] = v.in_direction(direction(i));
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,5,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif
  ivec vout(v.dim);
  for (int i=0; i<5; ++i) vout.set_direction(direction(i), out[i]);
  return vout;
}

double sum_to_master(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Reduce(&in,&out,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#endif
  return out;
}

double sum_to_all(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
  return out;
}

void sum_to_all(const double *in, double *out, int size) {
  memcpy(out, in, sizeof(double) * size);
#ifdef HAVE_MPI
  MPI_Allreduce((void*) in, out, size, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
}

long double sum_to_all(long double in) {
  long double out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,1,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
  return out;
}

int sum_to_all(int in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  return out;
}

int partial_sum_to_all(int in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Scan(&in,&out,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  return out;
}

complex<double> sum_to_all(complex<double> in) {
  complex<double> out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
  return out;
}

complex<long double> sum_to_all(complex<long double> in) {
  complex<long double> out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,2,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
  return out;
}

bool or_to_all(bool in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,1,MPI_INT,MPI_LOR,MPI_COMM_WORLD);
#endif
  return (bool) out;
}

bool and_to_all(bool in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in,&out,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
#endif
  return (bool) out;
}

void all_wait() {
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

int my_rank() {
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
#else
  return 0;
#endif
}

int count_processors() {
#ifdef HAVE_MPI
  int n;
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  return n;
#else
  return 1;
#endif
}

// IO Routines...

void master_printf(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  if (am_master()) { vprintf(fmt, ap); fflush(stdout); }
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

// Scary file writing...

file *everyone_open_write(const char *name) {
#ifdef HAVE_MPI
  const int buflen = strlen(name)+1;
  char *buf = new char[buflen];
  strcpy(buf,name);
  MPI_File myf;
  MPI_File_open(MPI_COMM_WORLD, buf,
                MPI_MODE_WRONLY | MPI_MODE_CREATE,
                NULL, &myf);
  return (file *)myf;
#else
  return (file *)fopen(name,"w");
#endif
}

file *everyone_open_write(const char *filename, const char *directory) {
  char name[500];
  snprintf(name, 500, "%s/%s", directory, filename);
#ifdef HAVE_MPI
  const int buflen = strlen(name)+1;
  char *buf = new char[buflen];
  strcpy(buf,name);
  MPI_File myf;
  MPI_File_open(MPI_COMM_WORLD, buf,
                MPI_MODE_WRONLY | MPI_MODE_CREATE,
                NULL, &myf);
  return (file *)myf;
#else
  return (file *)fopen(name,"w");
#endif
}

void everyone_close(file *f) {
#ifdef HAVE_MPI
  MPI_File_close((MPI_File *)&f);
#else
  fclose((FILE *)f);
#endif
}

void i_flush(file *f) {
#ifdef HAVE_MPI
  // I don't know how to flush over MPI...
#else
  fflush((FILE *)f);
#endif
}

#ifdef HAVE_MPI
static const int buflen = 1 << 20;
static char *buf = NULL;
#endif

void i_fprintf(file *f, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
#ifdef HAVE_MPI
  if (!buf) buf = new char[buflen];
  if (!buf) abort("Can't allocate buffer in i_fprintf!\n");
  int written = vsnprintf(buf, buflen, fmt, ap);
  if (written <= 0 || written > buflen)
    abort("Aaack can't write that much in i_fprintf!\n");
  MPI_Status stat;
  MPI_File_write_shared((MPI_File)f, buf, written, MPI_CHAR, &stat);
#else
  vfprintf((FILE *)f, fmt, ap);
#endif
  va_end(ap);
}

void master_fprintf(file *f, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  if (my_rank() == 0) {
#ifdef HAVE_MPI
    if (!buf) buf = new char[buflen];
    if (!buf) abort("Can't allocate buffer in master_fprintf!\n");
    int written = vsnprintf(buf, buflen, fmt, ap);
    if (written <= 0 || written > buflen)
      abort("Aaack can't write that much in master_fprintf!\n");
    MPI_Status stat;
    MPI_File_write_shared((MPI_File)f, buf, written, MPI_CHAR, &stat);
#else
    vfprintf((FILE *)f, fmt, ap);
#endif
  }
  va_end(ap);
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
     MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
     if (process_rank > 0) { /* wait for a message before continuing */
	  MPI_Status status;
	  int recv_tag = tag - 1; /* initialize to wrong value */
	  MPI_Recv(&recv_tag, 1, MPI_INT, process_rank - 1, tag, 
		   MPI_COMM_WORLD, &status);
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
     MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
     MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
     if (process_rank != num_procs - 1) { /* send a message to next process */
	  MPI_Send(&tag, 1, MPI_INT, process_rank + 1, tag, 
		   MPI_COMM_WORLD);
     }
#else
     UNUSED(tag);
#endif
}

} // namespace meep
