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

#include <cstdlib>
#include <stdarg.h>
#include <string.h>

#include "meep.hpp"
#include "config.h"

#ifdef HAVE_MPI
#ifdef NEED_UNDEF_SEEK_FOR_MPI
// undef'ing SEEK_* is needed for MPICH, possibly other MPI versions
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#endif
#include <mpi.h>
#endif

#ifdef _OPENMP
#include "omp.h"
#else
#define omp_get_num_threads() (1)
#endif

#ifdef IGNORE_SIGFPE
#include <signal.h>
#endif

#if defined(DEBUG) && defined(HAVE_FEENABLEEXCEPT)
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <fenv.h>
#if !HAVE_DECL_FEENABLEEXCEPT
extern "C" int feenableexcept(int EXCEPTS);
#endif
#endif

#if HAVE_SYS_TIME_H
#include <sys/time.h>
#include <time.h>
#else
#include <time.h>
#endif
#ifdef HAVE_BSDGETTIMEOFDAY
#ifndef HAVE_GETTIMEOFDAY
#define gettimeofday BSDgettimeofday
#define HAVE_GETTIMEOFDAY 1
#endif
#endif

#if HAVE_IMMINTRIN_H
#include <immintrin.h>
#endif

#define UNUSED(x) (void)x // silence compiler warnings

#define MPI_REALNUM (sizeof(realnum) == sizeof(double) ? MPI_DOUBLE : MPI_FLOAT)

using namespace std;

namespace meep {

namespace {

#ifdef HAVE_MPI
MPI_Comm mycomm = MPI_COMM_WORLD;
#endif

// comms_manager implementation that uses MPI.
class mpi_comms_manager : public comms_manager {
public:
  mpi_comms_manager() {}
  ~mpi_comms_manager() override {
#ifdef HAVE_MPI
    int num_pending_requests = reqs.size();
    std::vector<int> completed_indices(num_pending_requests);
    while (num_pending_requests) {
      int num_completed_requests = 0;
      MPI_Waitsome(reqs.size(), reqs.data(), &num_completed_requests, completed_indices.data(),
                   MPI_STATUSES_IGNORE);
      for (int i = 0; i < num_completed_requests; ++i) {
        int request_idx = completed_indices[i];
        callbacks[request_idx]();
        reqs[request_idx] = MPI_REQUEST_NULL;
        --num_pending_requests;
      }
    }
#endif
  }

  void send_real_async(const void *buf, size_t count, int dest, int tag) override {
#ifdef HAVE_MPI
    reqs.emplace_back();
    callbacks.push_back(/*no-op*/ []{});
    MPI_Isend(buf, static_cast<int>(count), MPI_REALNUM, dest, tag, mycomm, &reqs.back());
#else
    (void)buf;
    (void)count;
    (void)dest;
    (void)tag;
#endif
  }

  void receive_real_async(void *buf, size_t count, int source, int tag,
                          const receive_callback &cb) override {
#ifdef HAVE_MPI
    reqs.emplace_back();
    callbacks.push_back(cb);
    MPI_Irecv(buf, static_cast<int>(count), MPI_REALNUM, source, tag, mycomm, &reqs.back());
#else
    (void)buf;
    (void)count;
    (void)source;
    (void)tag;
    (void)cb;
#endif
  }

#ifdef HAVE_MPI
  size_t max_transfer_size() const override { return std::numeric_limits<int>::max(); }
#endif

private:
#ifdef HAVE_MPI
  std::vector<MPI_Request> reqs;
#endif
  std::vector<receive_callback> callbacks;
};

} // namespace

std::unique_ptr<comms_manager> create_comms_manager() {
  return std::unique_ptr<comms_manager>(new mpi_comms_manager());
}

int verbosity = 1; // defined in meep.h

/* Set CPU to flush subnormal values to zero (if iszero == true).  This slightly
   reduces the range of floating-point numbers, but can greatly increase the speed
   in cases where subnormal values might arise (e.g. deep in the tails of
   exponentially decaying sources).

   See also meep#1708.

   code based on github.com/JuliaLang/julia/blob/master/src/processor_x86.cpp#L1087-L1104,
   which is free software under the GPL-compatible "MIT license" */
static void _set_zero_subnormals(bool iszero) {
#if HAVE_IMMINTRIN_H
  unsigned int flags =
      0x00008040; // assume a non-ancient processor with SSE2, supporting both FTZ and DAZ flags
  unsigned int state = _mm_getcsr();
  if (iszero)
    state |= flags;
  else
    state &= ~flags;
  _mm_setcsr(state);
#else
  (void)iszero; // unused
#endif
}
void set_zero_subnormals(bool iszero) {
  int n = omp_get_num_threads();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1)
#endif
  for (int i = 0; i < n; ++i)
    _set_zero_subnormals(iszero); // This has to be done in every thread for OpenMP.
}

void setup() {
  set_zero_subnormals(true);
#ifdef _OPENMP
  if (getenv("OMP_NUM_THREADS") == NULL) omp_set_num_threads(1);
#endif
}

initialize::initialize(int &argc, char **&argv) {
#ifdef HAVE_MPI
#ifdef _OPENMP
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  if (provided < MPI_THREAD_FUNNELED && omp_get_num_threads() > 1)
    abort("MPI does not support multi-threaded execution");
#else
  MPI_Init(&argc, &argv);
#endif
  int major, minor;
  MPI_Get_version(&major, &minor);
  if (verbosity > 0)
    master_printf("Using MPI version %d.%d, %d processes\n", major, minor, count_processors());
#else
  UNUSED(argc);
  UNUSED(argv);
#endif
#if defined(DEBUG_FP) && defined(HAVE_FEENABLEEXCEPT)
  feenableexcept(FE_INVALID | FE_OVERFLOW); // crash if NaN created, or overflow
#endif
#ifdef IGNORE_SIGFPE
  signal(SIGFPE, SIG_IGN);
#endif
  t_start = wall_time();
  setup();
}

initialize::~initialize() {
  if (verbosity > 0) master_printf("\nElapsed run time = %g s\n", elapsed_time());
#ifdef HAVE_MPI
  end_divide_parallel();
  MPI_Finalize();
#endif
}

double wall_time(void) {
#ifdef HAVE_MPI
  return MPI_Wtime();
#elif defined(_OPENMP)
  return omp_get_wtime();
#elif HAVE_GETTIMEOFDAY
  struct timeval tv;
  gettimeofday(&tv, 0);
  return (tv.tv_sec + tv.tv_usec * 1e-6);
#else
  return (clock() * 1.0 / CLOCKS_PER_SECOND);
#endif
}

[[noreturn]] void abort(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  char *s;
  vasprintf(&s, fmt, ap);
  va_end(ap);
  // Make a std::string to support older compilers (std::runtime_error(char *) was added in C++11)
  std::string error_msg(s);
  free(s);
#ifdef HAVE_MPI
  if (count_processors() == 1) { throw runtime_error("meep: " + error_msg); }
  fprintf(stderr, "meep: %s", error_msg.c_str());
  if (fmt[strlen(fmt) - 1] != '\n') fputc('\n', stderr); // force newline
  MPI_Abort(MPI_COMM_WORLD, 1);
  std::abort(); // Unreachable but MPI_Abort does not have the noreturn attribute.
#else
  throw runtime_error("meep: " + error_msg);
#endif
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

void broadcast(int from, float *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, size, MPI_FLOAT, from, mycomm);
#else
  UNUSED(from);
  UNUSED(data);
  UNUSED(size);
#endif
}

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
  MPI_Bcast(data, 2 * size, MPI_DOUBLE, from, mycomm);
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

void broadcast(int from, size_t *data, int size) {
#ifdef HAVE_MPI
  if (size == 0) return;
  MPI_Bcast(data, size, sizeof(size_t) == 4 ? MPI_UNSIGNED : MPI_UNSIGNED_LONG_LONG, from, mycomm);
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

bool broadcast(int from, bool b) { return broadcast(from, (int)b); }

double max_to_master(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Reduce(&in, &out, 1, MPI_DOUBLE, MPI_MAX, 0, mycomm);
#endif
  return out;
}

double max_to_all(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE, MPI_MAX, mycomm);
#endif
  return out;
}

int max_to_all(int in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_MAX, mycomm);
#endif
  return out;
}

int min_to_all(int in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_MIN, mycomm);
#endif
  return out;
}

ivec max_to_all(const ivec &pt) {
  int in[5], out[5];
  for (int i = 0; i < 5; ++i)
    in[i] = out[i] = pt.in_direction(direction(i));
#ifdef HAVE_MPI
  MPI_Allreduce(&in, &out, 5, MPI_INT, MPI_MAX, mycomm);
#endif
  ivec ptout(pt.dim);
  for (int i = 0; i < 5; ++i)
    ptout.set_direction(direction(i), out[i]);
  return ptout;
}

float sum_to_master(float in) {
  float out = in;
#ifdef HAVE_MPI
  MPI_Reduce(&in, &out, 1, MPI_FLOAT, MPI_SUM, 0, mycomm);
#endif
  return out;
}

double sum_to_master(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Reduce(&in, &out, 1, MPI_DOUBLE, MPI_SUM, 0, mycomm);
#endif
  return out;
}

double sum_to_all(double in) {
  double out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE, MPI_SUM, mycomm);
#endif
  return out;
}

void sum_to_all(const float *in, float *out, int size) {
#ifdef HAVE_MPI
  MPI_Allreduce((void *)in, out, size, MPI_FLOAT, MPI_SUM, mycomm);
#else
  memcpy(out, in, sizeof(float) * size);
#endif
}

void sum_to_all(const double *in, double *out, int size) {
#ifdef HAVE_MPI
  MPI_Allreduce((void *)in, out, size, MPI_DOUBLE, MPI_SUM, mycomm);
#else
  memcpy(out, in, sizeof(double) * size);
#endif
}

void sum_to_master(const float *in, float *out, int size) {
#ifdef HAVE_MPI
  MPI_Reduce((void *)in, out, size, MPI_FLOAT, MPI_SUM, 0, mycomm);
#else
  memcpy(out, in, sizeof(float) * size);
#endif
}

void sum_to_master(const double *in, double *out, int size) {
#ifdef HAVE_MPI
  MPI_Reduce((void *)in, out, size, MPI_DOUBLE, MPI_SUM, 0, mycomm);
#else
  memcpy(out, in, sizeof(double) * size);
#endif
}

void sum_to_all(const float *in, double *out, int size) {
  double *in2 = new double[size];
  for (int i = 0; i < size; ++i)
    in2[i] = in[i];
  sum_to_all(in2, out, size);
  delete[] in2;
}

void sum_to_all(const complex<double> *in, complex<double> *out, int size) {
  sum_to_all((const double *)in, (double *)out, 2 * size);
}

void sum_to_all(const complex<float> *in, complex<double> *out, int size) {
  sum_to_all((const float *)in, (double *)out, 2 * size);
}

void sum_to_all(const complex<float> *in, complex<float> *out, int size) {
  sum_to_all((const float *)in, (float *)out, 2 * size);
}

void sum_to_master(const complex<float> *in, complex<float> *out, int size) {
  sum_to_master((const float *)in, (float *)out, 2 * size);
}

void sum_to_master(const complex<double> *in, complex<double> *out, int size) {
  sum_to_master((const double *)in, (double *)out, 2 * size);
}

long double sum_to_all(long double in) {
  long double out = in;
#ifdef HAVE_MPI
  if (MPI_LONG_DOUBLE == MPI_DATATYPE_NULL)
    out = sum_to_all(double(in));
  else
    MPI_Allreduce(&in, &out, 1, MPI_LONG_DOUBLE, MPI_SUM, mycomm);
#endif
  return out;
}

int sum_to_all(int in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_SUM, mycomm);
#endif
  return out;
}

int partial_sum_to_all(int in) {
  int out = in;
#ifdef HAVE_MPI
  MPI_Scan(&in, &out, 1, MPI_INT, MPI_SUM, mycomm);
#endif
  return out;
}

size_t sum_to_all(size_t in) {
  size_t out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in, &out, 1, sizeof(size_t) == 4 ? MPI_UNSIGNED : MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                mycomm);
#endif
  return out;
}

void sum_to_all(const size_t *in, size_t *out, int size) {
#ifdef HAVE_MPI
  MPI_Allreduce((void *)in, out, size, sizeof(size_t) == 4 ? MPI_UNSIGNED : MPI_UNSIGNED_LONG_LONG,
                MPI_SUM, mycomm);
#else
  memcpy(out, in, sizeof(size_t) * size);
#endif
}

void sum_to_master(const size_t *in, size_t *out, int size) {
#ifdef HAVE_MPI
  MPI_Reduce((void *)in, out, size, sizeof(size_t) == 4 ? MPI_UNSIGNED : MPI_UNSIGNED_LONG_LONG,
             MPI_SUM, 0, mycomm);
#else
  memcpy(out, in, sizeof(size_t) * size);
#endif
}

size_t partial_sum_to_all(size_t in) {
  size_t out = in;
#ifdef HAVE_MPI
  MPI_Scan(&in, &out, 1, sizeof(size_t) == 4 ? MPI_UNSIGNED : MPI_UNSIGNED_LONG_LONG, MPI_SUM,
           mycomm);
#endif
  return out;
}

complex<double> sum_to_all(complex<double> in) {
  complex<double> out = in;
#ifdef HAVE_MPI
  MPI_Allreduce(&in, &out, 2, MPI_DOUBLE, MPI_SUM, mycomm);
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
    MPI_Allreduce(&in, &out, 2, MPI_LONG_DOUBLE, MPI_SUM, mycomm);
#endif
  return out;
}

bool or_to_all(bool in) {
  int in2 = in, out;
#ifdef HAVE_MPI
  MPI_Allreduce(&in2, &out, 1, MPI_INT, MPI_LOR, mycomm);
#else
  out = in2;
#endif
  return (bool)out;
}

void or_to_all(const int *in, int *out, int size) {
#ifdef HAVE_MPI
  MPI_Allreduce((void *)in, out, size, MPI_INT, MPI_LOR, mycomm);
#else
  memcpy(out, in, sizeof(int) * size);
#endif
}

void bw_or_to_all(const size_t *in, size_t *out, int size) {
#ifdef HAVE_MPI
  MPI_Allreduce((void *)in, out, size, sizeof(size_t) == 4 ? MPI_UNSIGNED : MPI_UNSIGNED_LONG_LONG,
                MPI_BOR, mycomm);
#else
  memcpy(out, in, sizeof(size_t) * size);
#endif
}

bool and_to_all(bool in) {
  int in2 = in, out;
#ifdef HAVE_MPI
  MPI_Allreduce(&in2, &out, 1, MPI_INT, MPI_LAND, mycomm);
#else
  out = in2;
#endif
  return (bool)out;
}

void and_to_all(const int *in, int *out, int size) {
#ifdef HAVE_MPI
  MPI_Allreduce((void *)in, out, size, MPI_INT, MPI_LAND, mycomm);
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

bool with_mpi() {
#ifdef HAVE_MPI
  return true;
#else
  return false;
#endif
}

// IO Routines...

bool am_really_master() { return (my_global_rank() == 0); }

static meep_printf_callback_func master_printf_callback = NULL;
static meep_printf_callback_func master_printf_stderr_callback = NULL;

meep_printf_callback_func set_meep_printf_callback(meep_printf_callback_func func) {
  meep_printf_callback_func old_func = master_printf_callback;
  master_printf_callback = func;
  return old_func;
}

meep_printf_callback_func set_meep_printf_stderr_callback(meep_printf_callback_func func) {
  meep_printf_callback_func old_func = master_printf_stderr_callback;
  master_printf_stderr_callback = func;
  return old_func;
}

static void _do_master_printf(FILE *output, meep_printf_callback_func callback, const char *fmt,
                              va_list ap) {
  if (am_really_master()) {
    if (callback) {
      char *s;
      vasprintf(&s, fmt, ap);
      callback(s);
      free(s);
    }
    else {
      vfprintf(output, fmt, ap);
      fflush(output);
    }
  }
  va_end(ap);
}

void master_printf(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  _do_master_printf(stdout, master_printf_callback, fmt, ap);
  va_end(ap);
}

void master_printf_stderr(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  _do_master_printf(stderr, master_printf_stderr_callback, fmt, ap);
  va_end(ap);
}

static FILE *debf = NULL;

void debug_printf(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  if (debf == NULL) {
    char temp[50];
    snprintf(temp, 50, "debug_out_%d", my_rank());
    debf = fopen(temp, "w");
    if (!debf) meep::abort("Unable to open debug output %s\n", temp);
  }
  vfprintf(debf, fmt, ap);
  fflush(debf);
  va_end(ap);
}

void master_fprintf(FILE *f, const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  if (am_master()) {
    vfprintf(f, fmt, ap);
    fflush(f);
  }
  va_end(ap);
}
FILE *master_fopen(const char *name, const char *mode) {
  FILE *f = am_master() ? fopen(name, mode) : 0;

  /* other processes need to know if fopen returned zero, in order
     to abort if fopen failed.  If fopen was successfully, just return
     a random non-zero pointer (which is never used except to compare to zero)
     on non-master processes */
  if (broadcast(0, bool(f != 0)) && !am_master()) f = (FILE *)name;
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

void begin_critical_section(int tag) {
#ifdef HAVE_MPI
  int process_rank;
  MPI_Comm_rank(mycomm, &process_rank);
  if (process_rank > 0) { /* wait for a message before continuing */
    MPI_Status status;
    int recv_tag = tag - 1; /* initialize to wrong value */
    MPI_Recv(&recv_tag, 1, MPI_INT, process_rank - 1, tag, mycomm, &status);
    if (recv_tag != tag) meep::abort("invalid tag received in begin_critical_section");
  }
#else
  UNUSED(tag);
#endif
}

void end_critical_section(int tag) {
#ifdef HAVE_MPI
  int process_rank, num_procs;
  MPI_Comm_rank(mycomm, &process_rank);
  MPI_Comm_size(mycomm, &num_procs);
  if (process_rank != num_procs - 1) { /* send a message to next process */
    MPI_Send(&tag, 1, MPI_INT, process_rank + 1, tag, mycomm);
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

int divide_parallel_processes(int numgroups) {
#ifdef HAVE_MPI
  end_divide_parallel();
  if (numgroups > count_processors()) meep::abort("numgroups > count_processors");
  int mygroup = (my_rank() * numgroups) / count_processors();
  MPI_Comm_split(MPI_COMM_WORLD, mygroup, my_rank(), &mycomm);
  return mygroup;
#else
  if (numgroups != 1) meep::abort("cannot divide processes in non-MPI mode");
  return 0;
#endif
}

#ifdef HAVE_MPI
static MPI_Comm mycomm_save = MPI_COMM_WORLD;
#endif

void begin_global_communications(void) {
#ifdef HAVE_MPI
  mycomm_save = mycomm;
  mycomm = MPI_COMM_WORLD;
#endif
}

void end_global_communications(void) {
#ifdef HAVE_MPI
  mycomm = mycomm_save;
  mycomm_save = MPI_COMM_WORLD;
#endif
}

void end_divide_parallel(void) {
#ifdef HAVE_MPI
  if (mycomm != MPI_COMM_WORLD) MPI_Comm_free(&mycomm);
  if (mycomm_save != MPI_COMM_WORLD) MPI_Comm_free(&mycomm_save);
  mycomm = mycomm_save = MPI_COMM_WORLD;
#endif
}

int my_global_rank() {
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
#else
  return 0;
#endif
}

} // namespace meep
