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
#ifndef MEEP_MY_MPI_H
#define MEEP_MY_MPI_H

#include <complex>
using namespace std;

namespace meep {

// MPI helper routines!

class initialize {
 public:
  initialize(int argc, char **argv);
  ~initialize();
};

#ifdef __GNUC__
#  define NORETURN_ATTR __attribute__((noreturn))
#  define PRINTF_ATTR(f,a) __attribute__((format(printf, f, a)))
#else
#  define NORETURN_ATTR 
#  define PRINTF_ATTR(f,a) 
#endif

void abort(const char *fmt, ...) NORETURN_ATTR PRINTF_ATTR(1,2);
void all_wait();
int count_processors();
int my_rank();
inline int am_master() { return my_rank() == 0; }

void send(int from, int to, double *data, int size=1);
void broadcast(int from, double *data, int size);
void broadcast(int from, char *data, int size);
void broadcast(int from, int *data, int size);
void broadcast(int from, complex<double> *data, int size);
complex<double> broadcast(int from, complex<double> data);
double broadcast(int from, double data);
int broadcast(int from, int data);
bool broadcast(int from, bool);
double max_to_master(double); // Only returns the correct value to proc 0.
double max_to_all(double);
int max_to_all(int);
double sum_to_master(double); // Only returns the correct value to proc 0.
double sum_to_all(double);
long double sum_to_all(long double);
int sum_to_all(int);
complex<double> sum_to_all(complex<double>);

// IO routines:
void master_printf(const char *fmt, ...) PRINTF_ATTR(1,2);
void debug_printf(const char *fmt, ...) PRINTF_ATTR(1,2);

// File is an abstract type to keep you from accidentally using it in an
// unsafe manner.
class file;

file *everyone_open_write(const char *);
file *everyone_open_write(const char *filename, const char *directory);
void everyone_close(file *);
void i_flush(file *f);
void i_fprintf(file *, const char *fmt, ...) PRINTF_ATTR(2,3);
void master_fprintf(file *, const char *fmt, ...) PRINTF_ATTR(2,3);

void begin_critical_section(int tag);
void end_critical_section(int tag);

} /* namespace meep */

#endif /* MEEP_MY_MPI_H */
