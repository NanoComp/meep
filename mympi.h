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
#ifndef MY_MPI_H
#define MY_MPI_H

// MPI helper routines!

class initialize {
 public:
  initialize(int argc, char **argv);
  ~initialize();
};

void abort(char *fmt, ...);
void all_wait();
int count_processors();
int my_rank();
inline int am_master() { return my_rank() == 0; };

void send(int from, int to, double *data, int size=1);
void broadcast(int from, double *data, int size);
void broadcast(int from, complex<double> *data, int size);
bool broadcast(int from, bool);
complex<double> broadcast(int from, complex<double> data);
double max_to_master(double); // Only returns the correct value to proc 0.
double max_to_all(double);
double sum_to_master(double); // Only returns the correct value to proc 0.
double sum_to_all(double);

// IO routines:
void master_printf(const char *fmt, ...);
void debug_printf(const char *fmt, ...);

// File is an abstract type to keep you from accidentally using it in an
// unsafe manner.
class file;

file *everyone_open_write(const char *);
void everyone_close(file *);
void i_flush(file *f);
void i_fprintf(file *, const char *fmt, ...);
void master_fprintf(file *, const char *fmt, ...);

#endif
