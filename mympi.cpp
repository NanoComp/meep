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

void initialize(int argc, char **argv) {
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  master_printf("Using MPI...\n");
#endif
}

void finished() {
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
  const int me = my_rank();
  MPI_Request req;
  MPI_Status stat;
  if (from == me) {
    MPI_Isend(data, size, MPI_DOUBLE, to, 1, MPI_COMM_WORLD, &req);
    MPI_Wait(&req, &stat);
  }
  if (to == me) {
    MPI_Irecv(data, size, MPI_DOUBLE, from, 1, MPI_COMM_WORLD, &req);
    MPI_Wait(&req, &stat);
  }
#endif
}

double max_to_master(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Reduce(&in,&out,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
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

void sync() {
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

void master_printf(char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  if (am_master()) vprintf(fmt, ap);
  va_end(ap);
}
