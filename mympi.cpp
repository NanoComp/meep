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

#include "dactyl.h"
#include "config.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

void initialize(int argc, char **argv) {
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  if (!my_rank()) printf("Using MPI!!!\n");
#endif
}

void finished() {
#ifdef HAVE_MPI
  MPI_Finalize();
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
