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

#include "dactyl.h"
#include "config.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

initialize::initialize(int argc, char **argv) {
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  int major, minor;
  MPI_Get_version(&major, &minor);
  master_printf("Using MPI... version %d.%d\n", major, minor);
#endif
}

initialize::~initialize() {
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}

void abort(char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
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
#endif
}

void broadcast(int from, double *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, size, MPI_DOUBLE, from, MPI_COMM_WORLD);
#endif
}

void broadcast(int from, complex<double> *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, 2*size, MPI_DOUBLE, from, MPI_COMM_WORLD);
#endif
}

complex<double> broadcast(int from, complex<double> data) {
#ifdef HAVE_MPI
  MPI_Bcast(&data, 2, MPI_DOUBLE, from, MPI_COMM_WORLD);
#endif
  return data;
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
  if (am_master()) vprintf(fmt, ap);
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

void everyone_close(file *f) {
#ifdef HAVE_MPI
  MPI_File_close((MPI_File *)&f);
#else
  fclose((FILE *)f);
#endif
}

static const int buflen = 8192;
static char buf[buflen];

void i_fprintf(file *f, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
#ifdef HAVE_MPI
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
    int written = vsnprintf(buf, buflen, fmt, ap);
    if (written <= 0 || written > buflen)
      abort("Aaack can't write that much in i_fprintf!\n");
    MPI_Status stat;
    MPI_File_write_shared((MPI_File)f, buf, written, MPI_CHAR, &stat);
#else
    vfprintf((FILE *)f, fmt, ap);
#endif
  }
  va_end(ap);
}
