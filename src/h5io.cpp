/* Copyright (C) 2004 Massachusetts Institute of Technology
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


#include "h5io.h"

#include "mympi.h"
#define CHECK(condition, message) do { \
     if (!(condition))  { \
          abort("h5io error on line %d of " __FILE__ ": " \
                  message "\n", __LINE__); \
     } \
} while (0)

#include "config.h"

#ifdef HAVE_HDF5
#  include <hdf5.h>
#endif

namespace meep {

/*****************************************************************************/
/* If we have the H5Pset_fapl_mpio function (which is available if HDF5 was
   compiled for MPI), then we can perform collective file i/o operations
   (e.g. all processes call H5Fcreate at the same time to create one file).

   If we don't, however, then we deal with it by having one process
   work with the file at a time: "exclusive" access.  The following
   macro helps us select different bits of code depending upon whether
   this is the case. */

#ifdef HAVE_H5PSET_MPI /* old name for this routine */
#  define H5Pset_fapl_mpio H5Pset_mpi
#  ifndef HAVE_H5PSET_FAPL_MPIO
#    define HAVE_H5PSET_FAPL_MPIO 1
#  endif
#endif

#ifdef HAVE_H5PSET_FAPL_MPIO
#  define IF_EXCLUSIVE(yes,no) no
#else
#  define IF_EXCLUSIVE(yes,no) yes
#endif

static int matrixio_critical_section_tag = 0;

/*****************************************************************************/
/* Normally, HDF5 prints out all sorts of error messages, e.g. if a dataset
   can't be found, in addition to returning an error code.  The following
   macro can be wrapped around code to temporarily suppress error messages. */

#define SUPPRESS_HDF5_ERRORS(statements) { \
     H5E_auto_t xxxxx_err_func; \
     void *xxxxx_err_func_data; \
     H5Eget_auto(&xxxxx_err_func, &xxxxx_err_func_data); \
     H5Eset_auto(NULL, NULL); \
     { statements; } \
     H5Eset_auto(xxxxx_err_func, xxxxx_err_func_data); \
}

/*****************************************************************************/

#ifdef HAVE_HDF5
static bool dataset_exists(hid_t id, const char *name)
{
     hid_t data_id;
     SUPPRESS_HDF5_ERRORS(data_id = H5Dopen(id, name));
     if (data_id >= 0)
          H5Dclose(data_id);
     return (data_id >= 0);
}
#endif

double *h5io::read(const char *filename, const char *dataname,
		   int *rank, int *dims, int maxrank)
{
#ifdef HAVE_HDF5
     double *data;
     if (am_master()) {
	  int i, N;
	  hid_t file_id, space_id, data_id;
	  
	  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	  CHECK(file_id >= 0, "error opening HDF5 input file");
	  
	  if (!dataset_exists(file_id, dataname))
	       return NULL;
	  
	  data_id = H5Dopen(file_id, dataname);
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
	  
	  data = new double[N];
	  H5Dread(data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		  (void *) data);
	  
	  H5Dclose(data_id);
	  H5Fclose(file_id);
     }

     *rank = broadcast(0, *rank);
     broadcast(0, dims, *rank);
     int N = 1;
     for (int i = 0; i < *rank; ++i)
	  N *= dims[i];
     if (!am_master())
	  data = new double[N];
     broadcast(0, data, N);

     if (N == 1)
	  *rank = 0;

     return data;
#else
     return NULL;
#endif
}

/* Write a chunk of data to dataset <dataname> in HDF5 file
   <filename>.  The dataset has dimension dims[rank], and we are
   writing a chunk stored at <data> (row-major order) of size
   chunk_dims[rank], starting at chunk_start[rank].

   If append_data is true, then the file dataset has one extra
   dimension used for writing subsequent datasets, where the current
   dataset is at (nonnegative) location dindex in this dimension.

   If parallel is true, then all processes must call this routine
   simultaneously to write (non-overlapping) chunks in parallel, at
   the same dindex.  It should even be called from processes that have
   no data to write (i.e. chunk_dims[i] == 0).  If parallel is false,
   then this routine should only be called from a single process (with
   filesystem access).

   In the special case of rank == 0 (writing a single datum), chunk_dims[0]
   should still be initialized to 1 (if the given process is writing data)
   or 0 (if it is not).
   
   first_chunk should be set to true on the *first* chunk (or set of
   parallel chunks) written for *any* dindex.

   if single_precision is true, then the data is stored in the file
   using single precision; otherwise, double precision is used.

   If append_file is true, then the file can have pre-existing datasets
   with different names, which are not modified.  */
void h5io::write_chunk(const char *filename, const char *dataname,
		       int rank, const int *dims,
		       double *data,
		       const int *chunk_start, const int *chunk_dims,
		       bool append_data, int dindex,
		       bool parallel, bool first_chunk, 
		       bool single_precision,
		       bool append_file)
{
#ifdef HAVE_HDF5
     int i;
     bool do_write = true;
     hid_t file_id, space_id, mem_space_id, data_id;
     int rank1;

     CHECK(!append_data || dindex >= 0, "invalid dindex");
     CHECK(rank >= 0, "negative rank");
     CHECK(rank > 0 || chunk_dims[0] == 0 || chunk_dims[0] == 1,
	   "invalid chunk_dims[0] for rank 0");

     // stupid HDF5 has problems with rank 0
     rank1 = (rank == 0 && !append_data) ? 1 : rank;

     hid_t access_props = H5Pcreate (H5P_FILE_ACCESS);
#  if defined(HAVE_MPI) && defined(HAVE_H5PSET_FAPL_MPIO)
     if (parallel)
	  H5Pset_fapl_mpio(access_props, MPI_COMM_WORLD, MPI_INFO_NULL);
#  else
     if (parallel)
	  begin_critical_section(matrixio_critical_section_tag);
#  endif

     if (append_file || !first_chunk
	 || (IF_EXCLUSIVE(parallel && !am_master(), 1)))
	  file_id = H5Fopen(filename, H5F_ACC_RDWR, access_props);
     else
	  file_id = H5Fcreate(filename, H5F_ACC_TRUNC,
			      H5P_DEFAULT, access_props);

     H5Pclose(access_props);

     CHECK(file_id >= 0, "error opening HDF5 output file");     

     /* delete pre-existing datasets, or we'll have an error; I think
        we can only do this on the master process. (?) */
     if (first_chunk && dataset_exists(file_id, dataname)) {
          if (!parallel || am_master()) {
	       H5Gunlink(file_id, dataname);  /* delete it */
	       H5Fflush(file_id, H5F_SCOPE_GLOBAL);
	  }
	  IF_EXCLUSIVE((void) 0, if (parallel) all_wait());
     }
     
     if (first_chunk && (IF_EXCLUSIVE(!parallel || am_master(), 1))) {
	  hsize_t *dims_copy = new hsize_t[rank1 + append_data];
	  hsize_t *maxdims = new hsize_t[rank1 + append_data];
	  hsize_t N = 1;
	  for (i = 0; i < rank; ++i)
	       N *= (maxdims[i] = dims_copy[i] = dims[i]);
	  if (!rank)
	       maxdims[0] = dims_copy[0] = 1;
	  if (append_data) {
	       dims_copy[rank1] = 1;
	       maxdims[rank1] = H5S_UNLIMITED;
	  }
	  space_id = H5Screate_simple(rank1 + append_data, dims_copy, maxdims);
	  delete[] maxdims;

	  /* For unlimited datasets, we need to specify the size of the
	     "chunks" in which the file data is allocated.  We'll just
	     set the chunk size to the dataset dimensions (dims_copy). */
	  hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
	  if (append_data) {
	    const int blocksize = 128;
	    // make a chunk at least blocksize elements for efficiency
	    dims_copy[rank1] = (blocksize + (N - 1)) / N;
	    H5Pset_chunk(prop_id, rank1 + 1, dims_copy);
	    dims_copy[rank1] = 1;
	  }

	  delete[] dims_copy;
	  
	  hid_t type_id = 
	       single_precision ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE;

	  data_id = H5Dcreate(file_id, dataname, type_id, space_id, prop_id);

	  H5Pclose(prop_id);
     }
     else {
	  data_id = H5Dopen(file_id, dataname);
	  space_id = H5Dget_space(data_id);
	  
	  CHECK(rank1 + append_data == H5Sget_simple_extent_ndims(space_id),
		"file data is inconsistent rank for subsequent chunk");
	  
	  hsize_t *dims_copy = new hsize_t[rank1 + append_data];
	  hsize_t *maxdims = new hsize_t[rank1 + append_data];
	  H5Sget_simple_extent_dims(space_id, dims_copy, maxdims);
	  CHECK(!append_data || maxdims[rank1] == H5S_UNLIMITED,
		"file data is missing unlimited dimension for append_data");
	  delete[] maxdims;
	  for (i = 0; i < rank; ++i)
	       CHECK(dims[i] == (int) dims_copy[i],
		     "file data is inconsistent size for subsequent chunk");
	  if (rank < rank1)
	       CHECK(dims_copy[0] == 1, "rank-0 data is incorrect size");

	  // Allocate more space along unlimited direction, if needed:
	  if (append_data && dindex >= (int) dims_copy[rank1]) {
	       dims_copy[rank1] = dindex + 1;
	       H5Dextend(data_id, dims_copy);
	       H5Sclose(space_id); // I'm not sure if old space_id still valid
	       space_id = H5Dget_space(data_id);
	  }

	  delete[] dims_copy;
     }

     /*******************************************************************/
     /* Before we can write the data to the data set, we must define
	the dimensions and "selections" of the arrays to be read & written: */

     hssize_t *start = new hssize_t[rank1 + append_data];
     hsize_t *count = new hsize_t[rank1 + append_data];

     int count_prod = 1;
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
	  H5Sselect_hyperslab(space_id, H5S_SELECT_SET,
			      start, NULL, count, NULL);
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

     if (do_write)
	  H5Dwrite(data_id,
		   H5T_NATIVE_DOUBLE, mem_space_id, space_id, H5P_DEFAULT, 
		   (void *) data);

     H5Sclose(mem_space_id);
     H5Sclose(space_id);
     H5Fclose(file_id);

     IF_EXCLUSIVE(if (parallel)
		    end_critical_section(matrixio_critical_section_tag++),
		  (void) 0);
#endif
}

void h5io::write(const char *filename, const char *dataname,
		 double *data, int rank, const int *dims,
		 bool single_precision,
		 bool append_file)
{
  if (am_master()) {
    int *start = new int[rank + 1];
    for (int i = 0; i < rank; i++) start[i] = 0;
    h5io::write_chunk(filename, dataname, rank, dims, data, start, dims,
		      false,-1, false, true, single_precision, append_file);
    delete[] start;
  }
}

} // namespace meep
