/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <cstdio>
#include <cstdlib>
#include <string.h>

#include "meep.hpp"

#define CHECK(condition, message)                                                                  \
  do {                                                                                             \
    if (!(condition)) {                                                                            \
      meep::abort("error on line %d of " __FILE__ ": " message "\n", __LINE__);                    \
    }                                                                                              \
  } while (0)

#include "config.h"

#ifdef HAVE_HDF5

/* don't use new HDF5 1.8 API (which isn't even fully documented yet, grrr) */
#define H5_USE_16_API 1

#include <hdf5.h>

/* HDF5 changed this datatype in their interfaces starting in version 1.6.4 */
#if H5_VERS_MAJOR > 1 || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6) ||                              \
    (H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 6 && H5_VERS_RELEASE > 3)
typedef hsize_t start_t;
#else
typedef hssize_t start_t;
#endif

#else
typedef int hid_t;
#endif

#define HID(x) (*((hid_t *)(x)))

/*****************************************************************************/
/* If we have the H5Pset_fapl_mpio function (which is available if HDF5 was
   compiled for MPI), then we can perform collective file i/o operations
   (e.g. all processes call H5Fcreate at the same time to create one file).

   If we don't, however, then we deal with it by having one process
   work with the file at a time: "exclusive" access.  The following
   macro helps us select different bits of code depending upon whether
   this is the case. */

#ifdef HAVE_H5PSET_MPI /* old name for this routine */
#define H5Pset_fapl_mpio H5Pset_mpi
#ifndef HAVE_H5PSET_FAPL_MPIO
#define HAVE_H5PSET_FAPL_MPIO 1
#endif
#endif

#if defined(HAVE_H5PSET_FAPL_MPIO) || !defined(HAVE_MPI)
#define IF_EXCLUSIVE(yes, no) no
#else
#define IF_EXCLUSIVE(yes, no) yes
static int h5io_critical_section_tag = 0;
#endif

/*****************************************************************************/
/* Normally, HDF5 prints out all sorts of error messages, e.g. if a dataset
   can't be found, in addition to returning an error code.  The following
   macro can be wrapped around code to temporarily suppress error messages. */

#define SUPPRESS_HDF5_ERRORS(statements)                                                           \
  {                                                                                                \
    H5E_auto_t xxxxx_err_func;                                                                     \
    void *xxxxx_err_func_data;                                                                     \
    H5Eget_auto(&xxxxx_err_func, &xxxxx_err_func_data);                                            \
    H5Eset_auto(NULL, NULL);                                                                       \
    { statements; }                                                                                \
    H5Eset_auto(xxxxx_err_func, xxxxx_err_func_data);                                              \
  }

/*****************************************************************************/

using namespace std;

namespace meep {

bool h5file::dataset_exists(const char *name) {
#if HAVE_HDF5
  hid_t data_id;
  SUPPRESS_HDF5_ERRORS(data_id = H5Dopen(HID(get_id()), name));
  if (data_id >= 0) H5Dclose(data_id);
  return (data_id >= 0);
#else
  return false;
#endif
}

// lazy file creation & locking
void *h5file::get_id() {
  if (HID(id) < 0) {
    if (parallel) all_wait();

#ifdef HAVE_HDF5
    hid_t access_props = H5Pcreate(H5P_FILE_ACCESS);
#ifdef HAVE_MPI
#ifdef HAVE_H5PSET_FAPL_MPIO
    if (parallel) H5Pset_fapl_mpio(access_props, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
    if (parallel) begin_critical_section(h5io_critical_section_tag);
#endif
#endif

    if (mode != WRITE || IF_EXCLUSIVE(parallel && !am_master(), 0))
      HID(id) = H5Fopen(filename, mode == READONLY ? H5F_ACC_RDONLY : H5F_ACC_RDWR, access_props);
    else
      HID(id) = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, access_props);

    H5Pclose(access_props);
#endif
  }
  return id;
}

// hackery: in some circumstances, for the exclusive-access mode
// we must close the id (i.e. the file) in order to prevent deadlock.
void h5file::prevent_deadlock() { IF_EXCLUSIVE(if (parallel) close_id(), (void)0); }

void h5file::close_id() {
  unset_cur();
  if (HID(id) >= 0)
    if (mode == WRITE) mode = READWRITE; // don't re-create on re-open
#ifdef HAVE_HDF5
  if (HID(id) >= 0) {
    H5Fclose(HID(id));
    IF_EXCLUSIVE(if (parallel) end_critical_section(h5io_critical_section_tag++), (void)0);
  }
#endif
  HID(id) = -1;
}

/* note: if parallel is true, then *all* processes must call this,
   and all processes will use I/O. */
h5file::h5file(const char *filename_, access_mode m, bool parallel_, bool local_) {
  cur_dataname = NULL;
  id = (void *)malloc(sizeof(hid_t));
  cur_id = (void *)malloc(sizeof(hid_t));
  HID(id) = -1;
  HID(cur_id) = -1;
  extending = 0;
  filename = new char[strlen(filename_) + 1];
  strcpy(filename, filename_);
  mode = m;

  if (parallel_ && local_) {
    meep::abort("Can not open h5file (%s) in both parallel and local mode.", filename);
  }
  parallel = parallel_;
  local = local_;
}

h5file::~h5file() {
  close_id();
  if (cur_dataname) free(cur_dataname); // allocated with realloc
  for (h5file::extending_s *cur = extending; cur;) {
    h5file::extending_s *next = cur->next;
    delete[] cur->dataname;
    delete cur;
    cur = next;
  }
  delete[] filename;
  free(cur_id);
  free(id);
}

bool h5file::ok() { return (HID(get_id()) >= 0); }

void h5file::remove() {
  close_id();
  if (mode == READWRITE) mode = WRITE; // now need to re-create file
  for (h5file::extending_s *cur = extending; cur;) {
    h5file::extending_s *next = cur->next;
    delete[] cur->dataname;
    delete cur;
    cur = next;
  }
  extending = 0;

  IF_EXCLUSIVE(if (parallel) all_wait(), (void)0);
  if (am_master() || local) {
    if (std::remove(filename)) meep::abort("error removing file %s", filename);
  }
}

h5file::extending_s *h5file::get_extending(const char *dataname) const {
  for (extending_s *cur = extending; cur; cur = cur->next)
    if (!strcmp(dataname, cur->dataname)) return cur;
  return NULL;
}

bool h5file::is_cur(const char *dataname) {
  return cur_dataname && !strcmp(cur_dataname, dataname);
}

void h5file::unset_cur() {
#ifdef HAVE_HDF5
  if (HID(cur_id) >= 0) H5Dclose(HID(cur_id));
#endif
  HID(cur_id) = -1;
  if (cur_dataname) cur_dataname[0] = 0;
}

void h5file::set_cur(const char *dataname, void *data_id) {
#ifdef HAVE_HDF5
  if (HID(cur_id) >= 0 && HID(cur_id) != HID(data_id)) H5Dclose(HID(cur_id));
#endif
  HID(cur_id) = HID(data_id);
  if (!is_cur(dataname)) {
    if (!cur_dataname || strlen(dataname) > strlen(cur_dataname))
      cur_dataname = (char *)realloc(cur_dataname, strlen(dataname) + 1);
    strcpy(cur_dataname, dataname);
  }
}

void h5file::read_size(const char *dataname, int *rank, size_t *dims, int maxrank) {
#ifdef HAVE_HDF5
  if (parallel || am_master() || local) {
    hid_t file_id = HID(get_id()), space_id, data_id;

    CHECK(file_id >= 0, "error opening HDF5 input file");

    if (is_cur(dataname))
      data_id = HID(cur_id);
    else {
      CHECK(dataset_exists(dataname), "missing dataset in HDF5 file");
      data_id = H5Dopen(file_id, dataname);
      set_cur(dataname, &data_id);
    }
    space_id = H5Dget_space(data_id);

    *rank = H5Sget_simple_extent_ndims(space_id);
    CHECK(*rank <= maxrank, "input array rank is too big");

    hsize_t *dims_copy = new hsize_t[*rank];
    hsize_t *maxdims = new hsize_t[*rank];
    H5Sget_simple_extent_dims(space_id, dims_copy, maxdims);
    for (int i = 0; i < *rank; ++i)
      dims[i] = dims_copy[i];
    delete[] maxdims;
    delete[] dims_copy;
    H5Sclose(space_id);
  }

  if (!parallel && !local) {
    *rank = broadcast(0, *rank);
    broadcast(0, dims, *rank);

    if (*rank == 1 && dims[0] == 1) *rank = 0;
  }
#endif
}

#ifdef HAVE_HDF5
/* check if the given name is a dataset in group_id, and if so set d
   to point to a char** with a copy of name. */

static herr_t find_dataset(hid_t group_id, const char *name, void *d) {
  char **dname = (char **)d;
  H5G_stat_t info;

  H5Gget_objinfo(group_id, name, 1, &info);
  if (info.type == H5G_DATASET) {
    *dname = new char[strlen(name) + 1];
    strcpy(*dname, name);
    return 1;
  }
  return 0;
}

#endif

#ifdef HAVE_HDF5
#define SIZE_T_H5T (sizeof(size_t) == 4 ? H5T_NATIVE_UINT32 : H5T_NATIVE_UINT64)
#else
#define SIZE_T_H5T 0
#endif

void *h5file::read(const char *dataname, int *rank, size_t *dims, int maxrank,
                   bool single_precision) {
#ifdef HAVE_HDF5
  void *data = 0;
  if (parallel || am_master() || local) {
    int i, N;
    hid_t file_id = HID(get_id()), space_id, data_id;

    CHECK(file_id >= 0, "error opening HDF5 input file");

    bool close_data_id = true;
    if (dataname) {
      if (is_cur(dataname)) {
        data_id = HID(cur_id);
        close_data_id = false;
      }
      else {
        if (!dataset_exists(dataname)) {
          meep::abort("missing dataset in HDF5 file: %s", dataname);
        }
        data_id = H5Dopen(file_id, dataname);
      }
    }
    else { // find first dataset in file
      char *dname = 0;
      CHECK(H5Giterate(file_id, "/", NULL, find_dataset, &dname) >= 0 && dname,
            "cannot find dataset in HDF5 file");
      if (is_cur(dname)) {
        data_id = HID(cur_id);
        close_data_id = false;
      }
      else { data_id = H5Dopen(file_id, dname); }
      delete[] dname;
    }
    space_id = H5Dget_space(data_id);

    *rank = H5Sget_simple_extent_ndims(space_id);
    CHECK(*rank <= maxrank, "input array rank is too big");

    hsize_t *dims_copy = new hsize_t[*rank];
    hsize_t *maxdims = new hsize_t[*rank];
    H5Sget_simple_extent_dims(space_id, dims_copy, maxdims);
    delete[] maxdims;
    for (N = 1, i = 0; i < *rank; ++i)
      N *= (dims[i] = dims_copy[i]);
    delete[] dims_copy;
    H5Sclose(space_id);

    if (single_precision)
      data = new float[N];
    else
      data = new double[N];
    H5Dread(data_id, single_precision ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, (void *)data);

    if (close_data_id) H5Dclose(data_id);
  }

  if (!parallel && !local) {
    *rank = broadcast(0, *rank);
    broadcast(0, dims, *rank);
    size_t N = 1;
    for (int i = 0; i < *rank; ++i)
      N *= dims[i];
    if (!am_master()) {
      if (single_precision)
        data = new float[N];
      else
        data = new double[N];
    }
    if (single_precision)
      broadcast(0, (float *)data, N);
    else
      broadcast(0, (double *)data, N);
  }

  if (*rank == 1 && dims[0] == 1) *rank = 0;

  return data;
#else
  return NULL;
#endif
}

char *h5file::read(const char *dataname) {
#ifdef HAVE_HDF5
  char *data = 0;
  int len = 0;
  if (parallel || am_master() || local) {
    hid_t file_id = HID(get_id()), space_id, data_id, type_id;

    CHECK(file_id >= 0, "error opening HDF5 input file");

    if (is_cur(dataname)) unset_cur();

    CHECK(dataset_exists(dataname), "missing dataset in HDF5 file");

    data_id = H5Dopen(file_id, dataname);
    space_id = H5Dget_space(data_id);
    type_id = H5Dget_type(data_id);

    CHECK(H5Sget_simple_extent_npoints(space_id) == 1,
          "expected single string in HDF5 file, but didn't get one");

    len = H5Tget_size(type_id);
    H5Tclose(type_id);
    type_id = H5Tcopy(H5T_C_S1);
    H5Tset_size(type_id, len);

    data = new char[len];
    H5Dread(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void *)data);

    H5Tclose(type_id);
    H5Sclose(space_id);
    H5Dclose(data_id);
  }

  if (!parallel && !local) {
    len = broadcast(0, len);
    if (!am_master()) data = new char[len];
    broadcast(0, data, len);
  }

  return data;
#else
  return NULL;
#endif
}

/*****************************************************************************/

/* Delete a dataset, if it exists.  In parallel mode, should be called
   by all processors. */
void h5file::remove_data(const char *dataname) {
#ifdef HAVE_HDF5
  hid_t file_id = HID(get_id());

  if (is_cur(dataname)) unset_cur();

  if (get_extending(dataname)) { // delete dataname from extending list
    extending_s *prev = 0, *cur = extending;
    for (; cur && strcmp(cur->dataname, dataname); cur = (prev = cur)->next)
      ;
    if (!cur) meep::abort("bug in remove_data: inconsistent get_extending");
    if (prev)
      prev->next = cur->next;
    else
      extending = cur->next;
    delete[] cur->dataname;
    delete cur;
  }

  if (dataset_exists(dataname)) {
    /* this is hackish ...need to pester HDF5 developers to make
       H5Gunlink a collective operation for parallel mode */
    if (!parallel || am_master() || local) {
      H5Gunlink(file_id, dataname); /* delete it */
      H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    }
    IF_EXCLUSIVE((void)0, if (parallel) all_wait());
  }
#endif
}

/* Create a dataset, for writing chunks etc.  Note that, in parallel mode,
   this should be called by *all* processors, even those not writing any
   data. */
void h5file::create_data(const char *dataname, int rank, const size_t *dims, bool append_data,
                         bool single_precision) {
#ifdef HAVE_HDF5
  int i;
  hid_t file_id = HID(get_id()), space_id, data_id;
  int rank1;

  CHECK(rank >= 0, "negative rank");

  // stupid HDF5 has problems with rank 0
  rank1 = (rank == 0 && !append_data) ? 1 : rank;

  CHECK(file_id >= 0, "error opening HDF5 output file");

  unset_cur();
  remove_data(dataname); // HDF5 gives error if we H5Dcreate existing dataset

  if (local || IF_EXCLUSIVE(!parallel || am_master(), 1)) {
    hsize_t *dims_copy = new hsize_t[rank1 + append_data];
    hsize_t *maxdims = new hsize_t[rank1 + append_data];
    hsize_t N = 1;
    for (i = 0; i < rank; ++i)
      N *= (maxdims[i] = dims_copy[i] = dims[i]);
    if (!rank) maxdims[0] = dims_copy[0] = 1;
    if (append_data) {
      dims_copy[rank1] = 1;
      maxdims[rank1] = H5S_UNLIMITED;
    }
    space_id = H5Screate_simple(rank1 + append_data, dims_copy, maxdims);
    delete[] maxdims;

    /* For unlimited datasets, we need to specify the size of the
       "chunks" in which the file data is allocated.  */
    hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
    if (append_data) {
      const int blocksize = 128;
      // make a chunk at least blocksize elements for efficiency
      dims_copy[rank1] = (blocksize + (N - 1)) / N;
      H5Pset_chunk(prop_id, rank1 + 1, dims_copy);
      dims_copy[rank1] = 1;
    }

    delete[] dims_copy;

    hid_t type_id = single_precision ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE;

    data_id = H5Dcreate(file_id, dataname, type_id, space_id, prop_id);
    if (data_id < 0) meep::abort("Error creating dataset");

    H5Pclose(prop_id);
  }
  else {
    data_id = H5Dopen(file_id, dataname);
    CHECK(data_id >= 0, "missing dataset for subsequent processor");
    space_id = H5Dget_space(data_id);

    CHECK(rank1 + append_data == H5Sget_simple_extent_ndims(space_id),
          "file data is inconsistent rank for subsequent processor");

    hsize_t *dims_copy = new hsize_t[rank1 + append_data];
    hsize_t *maxdims = new hsize_t[rank1 + append_data];
    H5Sget_simple_extent_dims(space_id, dims_copy, maxdims);
    CHECK(!append_data || maxdims[rank1] == H5S_UNLIMITED,
          "file data is missing unlimited dimension for append_data");
    delete[] maxdims;
    for (i = 0; i < rank; ++i)
      CHECK(dims[i] == dims_copy[i], "file data is inconsistent size for subsequent processor");
    if (rank < rank1) CHECK(dims_copy[0] == 1, "rank-0 data is incorrect size");

    delete[] dims_copy;
  }

  set_cur(dataname, &data_id);
  H5Sclose(space_id);

  if (append_data) {
    extending_s *cur = new extending_s;
    cur->dataname = new char[strlen(dataname) + 1];
    strcpy(cur->dataname, dataname);
    cur->dindex = 0;
    cur->next = extending;
    extending = cur;
  }
#else
  meep::abort("not compiled with HDF5, required for HDF5 output");
#endif
}

/* Assumed data already created with append_data == true, and is
   already open; extends it and increments cur_dindex.  Like
   create_data, this is a collective operation and must be called from
   all processes. */
void h5file::extend_data(const char *dataname, int rank, const size_t *dims) {
#ifdef HAVE_HDF5
  extending_s *cur = get_extending(dataname);
  CHECK(cur, "extend_data can only be called on extensible data");

  hid_t file_id = HID(get_id()), data_id;
  if (is_cur(dataname))
    data_id = HID(cur_id);
  else {
    data_id = H5Dopen(file_id, dataname);
    set_cur(dataname, &data_id);
  }
  hid_t space_id = H5Dget_space(data_id);

  CHECK(rank + 1 == H5Sget_simple_extent_ndims(space_id),
        "file data is inconsistent rank for subsequent extend_data");
  hsize_t *dims_copy = new hsize_t[rank + 1];
  hsize_t *maxdims = new hsize_t[rank + 1];
  H5Sget_simple_extent_dims(space_id, dims_copy, maxdims);
  CHECK(maxdims[rank] == H5S_UNLIMITED, "file data is missing unlimited dimension for extend_data");
  delete[] maxdims;
  for (int i = 0; i < rank; ++i)
    CHECK(dims[i] == dims_copy[i], "file data is inconsistent size for subsequent extend_data");

  H5Sclose(space_id);

  // Allocate more space along unlimited direction
  cur->dindex++;
  dims_copy[rank] = cur->dindex + 1;
  H5Dextend(data_id, dims_copy);

  delete[] dims_copy;

#else
  meep::abort("not compiled with HDF5, required for HDF5 output");
#endif
}

/* If append_data is true, dataname is the current dataset, and is
   extensible, then as extend_data; otherwise as create_data. */
void h5file::create_or_extend_data(const char *dataname, int rank, const size_t *dims,
                                   bool append_data, bool single_precision) {
  if (get_extending(dataname))
    extend_data(dataname, rank, dims);
  else
    create_data(dataname, rank, dims, append_data, single_precision);
}

/*****************************************************************************/

/* Write a chunk of data to dataset <dataname> in HDF5 file.
   The dataset has dimension dims[rank], and we are
   writing a chunk stored at <data> (row-major order) of size
   chunk_dims[rank], starting at chunk_start[rank].

   You *must* have already called create_data for the same dimensions
   (and extend_data, if necessary).

   In the special case of rank == 0 (writing a single datum), chunk_dims[0]
   should still be initialized to 1 (if the given process is writing data)
   or 0 (if it is not).

   This function does *not* need to be called on all CPUs (e.g. those
   that have no data can be skipped).
*/
static void _write_chunk(hid_t data_id, h5file::extending_s *cur, int rank,
                         const size_t *chunk_start, const size_t *chunk_dims, hid_t datatype,
                         void *data) {
#ifdef HAVE_HDF5
  int i;
  bool do_write = true;
  hid_t space_id, mem_space_id;
  int rank1;
  bool append_data = cur != NULL;
  int dindex = cur ? cur->dindex : 0;

  CHECK(data_id >= 0, "create_data must be called before write_chunk");

  CHECK(rank >= 0, "negative rank");
  CHECK(rank > 0 || chunk_dims[0] == 0 || chunk_dims[0] == 1, "invalid chunk_dims[0] for rank 0");

  // stupid HDF5 has problems with rank 0
  rank1 = (rank == 0 && !append_data) ? 1 : rank;

  space_id = H5Dget_space(data_id);

  /*******************************************************************/
  /* Before we can write the data to the data set, we must define
     the dimensions and "selections" of the arrays to be read & written: */

  start_t *start = new start_t[rank1 + append_data];
  hsize_t *count = new hsize_t[rank1 + append_data];

  size_t count_prod = 1;
  for (i = 0; i < rank; ++i) {
    start[i] = chunk_start[i];
    count[i] = chunk_dims[i];
    count_prod *= count[i];
  }
  if (!rank) {
    start[0] = 0;
    count[0] = chunk_dims[0]; // see comment at top
    count_prod *= count[0];
  }
  if (append_data) {
    start[rank1] = dindex;
    count[rank1] = 1;
  }

  if (count_prod > 0) {
    H5Sselect_hyperslab(space_id, H5S_SELECT_SET, start, NULL, count, NULL);
    mem_space_id = H5Screate_simple(!rank1 ? 1 : rank1, count, NULL);
    H5Sselect_all(mem_space_id);
  }
  else { /* this can happen on leftover processes in MPI */
    H5Sselect_none(space_id);
    mem_space_id = H5Scopy(space_id); /* can't create an empty space */
    H5Sselect_none(mem_space_id);
    do_write = false; /* HDF5 complains about empty dataspaces */
  }

  delete[] start;
  delete[] count;

  /*******************************************************************/
  /* Write the data, then free all the stuff we've allocated. */

  if (do_write) H5Dwrite(data_id, datatype, mem_space_id, space_id, H5P_DEFAULT, (void *)data);

  H5Sclose(mem_space_id);
  H5Sclose(space_id);
#else
  meep::abort("not compiled with HDF5, required for HDF5 output");
#endif
}

void h5file::write_chunk(int rank, const size_t *chunk_start, const size_t *chunk_dims,
                         float *data) {
  _write_chunk(HID(cur_id), get_extending(cur_dataname), rank, chunk_start, chunk_dims,
               H5T_NATIVE_FLOAT, data);
}

void h5file::write_chunk(int rank, const size_t *chunk_start, const size_t *chunk_dims,
                         double *data) {
  _write_chunk(HID(cur_id), get_extending(cur_dataname), rank, chunk_start, chunk_dims,
               H5T_NATIVE_DOUBLE, data);
}

void h5file::write_chunk(int rank, const size_t *chunk_start, const size_t *chunk_dims,
                         size_t *data) {
  _write_chunk(HID(cur_id), get_extending(cur_dataname), rank, chunk_start, chunk_dims, SIZE_T_H5T,
               (void *)data);
}

// collective call after completing all write_chunk calls
void h5file::done_writing_chunks() {
  /* hackery: in order to not deadlock when writing extensible datasets
     with a non-parallel version of HDF5, we need to close the file
     and release the lock after writing extensible chunks  ...here,
     I'm assuming(?) that non-extensible datasets will use different
     files, etcetera, for different timesteps.  All of this hackery
     goes away if we just use an MPI-compiled version of HDF5. */
  if (parallel && !local && cur_dataname && get_extending(cur_dataname))
    prevent_deadlock(); // closes id
}

void h5file::write(const char *dataname, int rank, const size_t *dims, void *data,
                   bool single_precision) {
  if (parallel || am_master() || local) {
    size_t *start = new size_t[rank + 1];
    for (int i = 0; i < rank; i++)
      start[i] = 0;
    create_data(dataname, rank, dims, false, single_precision);
    if (am_master()) {
      if (single_precision)
        write_chunk(rank, start, dims, (float *)data);
      else
        write_chunk(rank, start, dims, (double *)data);
    }
    done_writing_chunks();
    unset_cur();
    delete[] start;
  }
}

void h5file::write(const char *dataname, const char *data) {
#ifdef HAVE_HDF5
  if (local || IF_EXCLUSIVE(am_master(), parallel || am_master())) {
    hid_t file_id = HID(get_id()), type_id, data_id, space_id;

    CHECK(file_id >= 0, "error opening HDF5 output file");

    remove_data(dataname); // HDF5 gives error if we H5Dcreate existing dataset

    type_id = H5Tcopy(H5T_C_S1);
    H5Tset_size(type_id, strlen(data) + 1);
    space_id = H5Screate(H5S_SCALAR);

    data_id = H5Dcreate(file_id, dataname, type_id, space_id, H5P_DEFAULT);
    if (am_master()) H5Dwrite(data_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Sclose(space_id);
    H5Tclose(type_id);
    H5Dclose(data_id);
  }
#else
  meep::abort("not compiled with HDF5, required for HDF5 output");
#endif
}

/*****************************************************************************/

/* Inverse of write_chunk, above.  The caller must first get the
   total dataset's rank and dims first by calling read_size, above,
   (which also opens the dataset for reading). */
static void _read_chunk(hid_t data_id, int rank, const size_t *chunk_start,
                        const size_t *chunk_dims, hid_t datatype, void *data) {
#ifdef HAVE_HDF5
  bool do_read = true;
  int rank1;
  hid_t space_id, mem_space_id;

  CHECK(data_id >= 0, "read_size must be called before read_chunk");

  CHECK(rank >= 0, "negative rank");
  CHECK(rank > 0 || chunk_dims[0] == 0 || chunk_dims[0] == 1, "invalid chunk_dims[0] for rank 0");

  // stupid HDF5 has problems with rank 0
  rank1 = rank == 0 ? 1 : rank;

  space_id = H5Dget_space(data_id);

  /*******************************************************************/
  /* Before we can read the data from the data set, we must define
     the dimensions and "selections" of the arrays to be read & written: */

  start_t *start = new start_t[rank1];
  hsize_t *count = new hsize_t[rank1];

  size_t count_prod = 1;
  for (int i = 0; i < rank; ++i) {
    start[i] = chunk_start[i];
    count[i] = chunk_dims[i];
    count_prod *= count[i];
  }
  if (!rank) {
    start[0] = 0;
    count[0] = chunk_dims[0]; // see comment at top
    count_prod *= count[0];
  }

  if (count_prod > 0) {
    H5Sselect_hyperslab(space_id, H5S_SELECT_SET, start, NULL, count, NULL);
    mem_space_id = H5Screate_simple(rank1, count, NULL);
    H5Sselect_all(mem_space_id);
  }
  else { /* this can happen on leftover processes in MPI */
    H5Sselect_none(space_id);
    mem_space_id = H5Scopy(space_id); /* can't create an empty space */
    H5Sselect_none(mem_space_id);
    do_read = false; /* HDF5 complains about empty dataspaces */
  }

  delete[] count;
  delete[] start;

  /*******************************************************************/
  /* Read the data, then free all the stuff we've allocated. */

  if (do_read) H5Dread(data_id, datatype, mem_space_id, space_id, H5P_DEFAULT, (void *)data);

  H5Sclose(mem_space_id);
  H5Sclose(space_id);
#else
  meep::abort("not compiled with HDF5, required for HDF5 input");
#endif
}

void h5file::read_chunk(int rank, const size_t *chunk_start, const size_t *chunk_dims,
                        float *data) {
  _read_chunk(HID(cur_id), rank, chunk_start, chunk_dims, H5T_NATIVE_FLOAT, data);
}

void h5file::read_chunk(int rank, const size_t *chunk_start, const size_t *chunk_dims,
                        double *data) {
  _read_chunk(HID(cur_id), rank, chunk_start, chunk_dims, H5T_NATIVE_DOUBLE, data);
}

void h5file::read_chunk(int rank, const size_t *chunk_start, const size_t *chunk_dims,
                        size_t *data) {
  _read_chunk(HID(cur_id), rank, chunk_start, chunk_dims, SIZE_T_H5T, (void *)data);
}

} // namespace meep
