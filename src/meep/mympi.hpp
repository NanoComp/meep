// -*- C++ -*-
/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
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
#include <stddef.h>
#include <stdexcept>

namespace meep {

// MPI helper routines!

double wall_time(void);

class initialize {
public:
  initialize(int &argc, char **&argv);
  ~initialize();
  double elapsed_time() { return wall_time() - t_start; }

private:
  double t_start;
};

#ifdef __GNUC__
#define NORETURN_ATTR __attribute__((noreturn))
#define PRINTF_ATTR(f, a) __attribute__((format(printf, f, a)))
#else
#define NORETURN_ATTR
#define PRINTF_ATTR(f, a)
#endif

#ifdef SWIG
void abort(const char *fmt, ...) PRINTF_ATTR(1, 2);
#else
[[noreturn]] void abort(const char *fmt, ...) PRINTF_ATTR(1, 2);
#endif
void all_wait();
int count_processors();
int my_rank();
bool am_really_master();
inline int am_master() { return my_rank() == 0; }
bool with_mpi();

void send(int from, int to, double *data, int size = 1);
void broadcast(int from, float *data, int size);
void broadcast(int from, double *data, int size);
void broadcast(int from, char *data, int size);
void broadcast(int from, int *data, int size);
void broadcast(int from, size_t *data, int size);
void broadcast(int from, std::complex<double> *data, int size);
std::complex<double> broadcast(int from, std::complex<double> data);
double broadcast(int from, double data);
int broadcast(int from, int data);
bool broadcast(int from, bool);
double max_to_master(double); // Only returns the correct value to proc 0.
double max_to_all(double);
int max_to_all(int);
int min_to_all(int);
float sum_to_master(float);   // Only returns the correct value to proc 0.
double sum_to_master(double); // Only returns the correct value to proc 0.
double sum_to_all(double);
void sum_to_all(const float *in, float *out, int size);
void sum_to_all(const double *in, double *out, int size);
void sum_to_master(const float *in, float *out, int size);
void sum_to_master(const double *in, double *out, int size);
void sum_to_all(const float *in, double *out, int size);
void sum_to_all(const std::complex<float> *in, std::complex<double> *out, int size);
void sum_to_all(const std::complex<double> *in, std::complex<double> *out, int size);
void sum_to_all(const std::complex<float> *in, std::complex<float> *out, int size);
void sum_to_master(const std::complex<float> *in, std::complex<float> *out, int size);
void sum_to_master(const std::complex<double> *in, std::complex<double> *out, int size);
long double sum_to_all(long double);
std::complex<double> sum_to_all(std::complex<double> in);
std::complex<long double> sum_to_all(std::complex<long double> in);
int sum_to_all(int);
int partial_sum_to_all(int in);
size_t sum_to_all(size_t);
size_t partial_sum_to_all(size_t in);
void sum_to_all(const size_t *in, size_t *out, int size);
void sum_to_master(const size_t *in, size_t *out, int size);
bool or_to_all(bool in);
void or_to_all(const int *in, int *out, int size);
void bw_or_to_all(const size_t *in, size_t *out, int size);
bool and_to_all(bool in);
void and_to_all(const int *in, int *out, int size);

// IO routines:
void master_printf(const char *fmt, ...) PRINTF_ATTR(1, 2);
void master_printf_stderr(const char *fmt, ...) PRINTF_ATTR(1, 2);
void debug_printf(const char *fmt, ...) PRINTF_ATTR(1, 2);
void master_fprintf(FILE *f, const char *fmt, ...) PRINTF_ATTR(2, 3);
FILE *master_fopen(const char *name, const char *mode);
void master_fclose(FILE *f);
typedef void (*meep_printf_callback_func)(const char *s);
meep_printf_callback_func set_meep_printf_callback(meep_printf_callback_func func);
meep_printf_callback_func set_meep_printf_stderr_callback(meep_printf_callback_func func);

void begin_critical_section(int tag);
void end_critical_section(int tag);

int divide_parallel_processes(int numgroups);
void begin_global_communications(void);
void end_global_communications(void);
void end_divide_parallel(void);

int my_global_rank(void);

} /* namespace meep */

#endif /* MEEP_MY_MPI_H */
