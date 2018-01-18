/****** matrixio ******/

// For now this is copied from mpb/matrixio. It's currently linked statically with
// the mpb binary, but maybe we can add it to libmpb.so.

#define CHECK(condition, message) /* do nothing */
#define FNAME_SUFFIX ".h5"

#ifndef HAVE_H5PSET_FAPL_MPIO
static int matrixio_critical_section_tag = 0;
#endif

#ifdef HAVE_H5PSET_FAPL_MPIO
#  define IF_EXCLUSIVE(yes,no) no
#else
#  define IF_EXCLUSIVE(yes,no) yes
#endif

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

matrixio_id matrixio_create(const char *fname);
void fieldio_write_complex_field(scalar_complex *field,
                                 int rank,
                                 const int dims[3],
                                 const int local_dims[3],
                                 const int start[3],
                                 int which_component,
                                 int num_components,
                                 const mpb_real *kvector,
                                 matrixio_id file_id,
                                 int append,
                                 matrixio_id data_id[]);
void fieldio_write_real_vals(mpb_real *vals,
                             int rank,
                             const int dims[3],
                             const int local_dims[3],
                             const int start[3],
                             matrixio_id file_id,
                             int append,
                             const char *dataname,
                             matrixio_id *data_id);
matrixio_id matrixio_create_dataset(matrixio_id id, const char *name, const char *description,
                                    int rank, const int *dims);
void matrixio_write_real_data(matrixio_id data_id, const int *local_dims, const int *local_start,
                              int stride, mpb_real *data);
void matrixio_close_dataset(matrixio_id data_id);
void matrixio_write_data_attr(matrixio_id id, const char *name, const mpb_real *val,
                              int rank, const int *dims);
void matrixio_write_string_attr(matrixio_id id, const char *name, const char *val);
matrixio_id matrixio_open(const char *fname, int read_only);
void matrixio_close(matrixio_id id);
void evectmatrixio_readall_raw(const char *filename, evectmatrix a);
char *add_fname_suffix(const char *fname);
int matrixio_dataset_exists(matrixio_id id, const char *name);
void matrixio_dataset_delete(matrixio_id id, const char *name);
void write_attr(matrixio_id id, matrixio_id_ type_id, matrixio_id_ space_id,
                const char *name, const void *val);

static matrixio_id matrixio_create_(const char *fname, int parallel)
{
#if defined(HAVE_HDF5)
  char *new_fname;
  matrixio_id id;
  hid_t access_props;

  access_props = H5Pcreate (H5P_FILE_ACCESS);

#if defined(HAVE_MPI) && defined(HAVE_H5PSET_FAPL_MPIO)
  if (parallel) {
    CHECK(H5Pset_fapl_mpio(access_props, mpb_comm, MPI_INFO_NULL)
          >= 0, "error initializing MPI file access");
  }
#endif

  new_fname = add_fname_suffix(fname);

#ifdef HAVE_H5PSET_FAPL_MPIO
  id.id = H5Fcreate(new_fname, H5F_ACC_TRUNC, H5P_DEFAULT, access_props);
#else
  if (parallel) mpi_begin_critical_section(matrixio_critical_section_tag);
  if (mpi_is_master() || !parallel)
    id.id = H5Fcreate(new_fname, H5F_ACC_TRUNC,H5P_DEFAULT,access_props);
  else
    id.id = H5Fopen(new_fname, H5F_ACC_RDWR, access_props);
#endif
  id.parallel = parallel;

  CHECK(id.id >= 0, "error creating HDF output file");

  free(new_fname);

  H5Pclose(access_props);

  return id;
#else
  meep::master_fprintf(stderr, "matrixio: cannot output \"%s\" (compiled without HDF)\n", fname);
  {
    matrixio_id id = {0,0};
    return id;
  }
#endif
}

matrixio_id matrixio_create(const char *fname) {
  return matrixio_create_(fname, 1);
}

/* note that kvector here is given in the reciprocal basis
   ...data_id should be of length at 2*num_components */
void fieldio_write_complex_field(scalar_complex *field,
                                 int rank,
                                 const int dims[3],
                                 const int local_dims[3],
                                 const int start[3],
                                 int which_component, int num_components,
                                 const mpb_real *kvector,
                                 matrixio_id file_id,
                                 int append,
                                 matrixio_id data_id[]) {

  int i, j, k, component, ri_part;

  rank = dims[2] == 1 ? (dims[1] == 1 ? 1 : 2) : 3;

  if (kvector) {
    mpb_real s[3]; /* the step size between grid points dotted with k */

    for (i = 0; i < 3; ++i)
      s[i] = TWOPI * kvector[i] / dims[i];

    /* cache exp(ikx) along each of the directions, for speed */
    scalar_complex *phasex = (scalar_complex*)malloc(sizeof(scalar_complex) * local_dims[0]);
    scalar_complex *phasey = (scalar_complex*)malloc(sizeof(scalar_complex) * local_dims[1]);
    scalar_complex *phasez = (scalar_complex*)malloc(sizeof(scalar_complex) * local_dims[2]);

    for (i = 0; i < local_dims[0]; ++i) {
      mpb_real phase = s[0] * (i + start[0]);
      phasex[i].re = cos(phase);
      phasex[i].im = sin(phase);
    }
    for (j = 0; j < local_dims[1]; ++j) {
      mpb_real phase = s[1] * (j + start[1]);
      phasey[j].re = cos(phase);
      phasey[j].im = sin(phase);
    }
    for (k = 0; k < local_dims[2]; ++k) {
      mpb_real phase = s[2] * (k + start[2]);
      phasez[k].re = cos(phase);
      phasez[k].im = sin(phase);
    }

    /* Now, multiply field by exp(i k*r): */
    for (i = 0; i < local_dims[0]; ++i) {
      scalar_complex px = phasex[i];

      for (j = 0; j < local_dims[1]; ++j) {
        scalar_complex py;
        mpb_real re = phasey[j].re, im = phasey[j].im;
        py.re = px.re * re - px.im * im;
        py.im = px.re * im + px.im * re;

        for (k = 0; k < local_dims[2]; ++k) {
          int ijk = ((i*local_dims[1] + j)*local_dims[2] + k)*3;
          mpb_real p_re, p_im;
          mpb_real re = phasez[k].re, im = phasez[k].im;

          p_re = py.re * re - py.im * im;
          p_im = py.re * im + py.im * re;

          for (component = 0; component < 3; ++component) {
            int ijkc = ijk + component;
            re = field[ijkc].re; im = field[ijkc].im;
            field[ijkc].re = re * p_re - im * p_im;
            field[ijkc].im = im * p_re + re * p_im;
          }
        }
      }
    }

    free(phasez);
    free(phasey);
    free(phasex);
  }

  /* write hyperslabs for each field component: */
  for (component = 0; component < num_components; ++component) {
    if (component == which_component || which_component < 0) {
      for (ri_part = 0; ri_part < 2; ++ri_part) {
        char name[] = "x.i";
        name[0] = (num_components == 1 ? 'c' : 'x') + component;
        name[2] = ri_part ? 'i' : 'r';

        if (!append)
          data_id[component*2 + ri_part] = matrixio_create_dataset(file_id, name, NULL,
                                                                      rank, dims);

        matrixio_write_real_data(data_id[component*2 + ri_part], local_dims, start, 2 * num_components,
                                 ri_part ? &field[component].im : &field[component].re);
      }
    }
  }
}

void fieldio_write_real_vals(mpb_real *vals, int rank, const int dims[3], const int local_dims[3],
                             const int start[3], matrixio_id file_id, int append, const char *dataname,
                             matrixio_id *data_id) {

  rank = dims[2] == 1 ? (dims[1] == 1 ? 1 : 2) : 3;

  if (!append || data_id->id < 0)
    *data_id = matrixio_create_dataset(file_id, dataname, NULL, rank,dims);

  matrixio_write_real_data(*data_id,local_dims,start,1,vals);
}

matrixio_id matrixio_create_dataset(matrixio_id id, const char *name, const char *description,
                                    int rank, const int *dims) {

  matrixio_id data_id;
  data_id.id = 0;
  data_id.parallel = id.parallel;
#if defined(HAVE_HDF5)
  {
    int i;
    hid_t space_id, type_id;
    hsize_t *dims_copy;

    /* delete pre-existing datasets, or we'll have an error; I think
       we can only do this on the master process. (?) */
    if (matrixio_dataset_exists(id, name)) {
#ifdef HAVE_H5PSET_FAPL_MPIO /* H5Gunlink is collective */
      matrixio_dataset_delete(id, name);
#else
      if (mpi_is_master() || !id.parallel) {
        matrixio_dataset_delete(id, name);
        H5Fflush(id.id, H5F_SCOPE_GLOBAL);
      }
      IF_EXCLUSIVE(0,if (id.parallel) MPI_Barrier(mpb_comm));
#endif
    }

    CHECK(rank > 0, "non-positive rank");

    dims_copy = (hsize_t *)malloc(sizeof(hsize_t) * rank);
    for (i = 0; i < rank; ++i)
      dims_copy[i] = dims[i];

    space_id = H5Screate_simple(rank, dims_copy, NULL);

    free(dims_copy);

#if defined(SCALAR_SINGLE_PREC)
    type_id = H5T_NATIVE_FLOAT;
#elif defined(SCALAR_LONG_DOUBLE_PREC)
    type_id = H5T_NATIVE_LDOUBLE;
#else
    type_id = H5T_NATIVE_DOUBLE;
#endif

    /* Create the dataset.  Note that, on parallel machines, H5Dcreate
       should do the right thing; it is supposedly a collective operation. */
    IF_EXCLUSIVE(
    if (mpi_is_master() || !id.parallel)
      data_id.id = H5Dcreate(id.id,name,type_id,space_id,H5P_DEFAULT);
    else
      data_id.id = H5Dopen(id.id, name), data_id.id = H5Dcreate(id.id, name, type_id, space_id, H5P_DEFAULT));

    H5Sclose(space_id);  /* the dataset should have its own copy now */

    matrixio_write_string_attr(data_id, "description", description);
  }
#endif
  return data_id;
}

void matrixio_write_real_data(matrixio_id data_id, const int *local_dims, const int *local_start,
                              int stride, mpb_real *data) {
#if defined(HAVE_HDF5)
  int rank;
  hid_t space_id, type_id, mem_space_id;
  start_t *start;
  hsize_t *strides, *count, count_prod;
  int i;
  mpb_real *data_copy;
  int data_copy_stride = 1, free_data_copy = 0, do_write = 1;

  /*******************************************************************/
  /* Get dimensions of dataset */

  space_id = H5Dget_space(data_id.id);

  rank = H5Sget_simple_extent_ndims(space_id);

  hsize_t *dims = (hsize_t *)malloc(sizeof(hsize_t) * rank);
  hsize_t *maxdims = (hsize_t *)malloc(sizeof(hsize_t) * rank);

  H5Sget_simple_extent_dims(space_id, dims, maxdims);

  free(maxdims);

#if defined(SCALAR_SINGLE_PREC)
  type_id = H5T_NATIVE_FLOAT;
#elif defined(SCALAR_LONG_DOUBLE_PREC)
  type_id = H5T_NATIVE_LDOUBLE;
#else
  type_id = H5T_NATIVE_DOUBLE;
#endif

  /*******************************************************************/
  /* if stride > 1, make a contiguous copy; hdf5 is much faster
     in this case. */

  if (stride > 1) {
    int N = 1;
    for (i = 0; i < rank; ++i)
      N *= local_dims[i];
    data_copy = (mpb_real *)malloc(sizeof(mpb_real) * N);
    if (data_copy) {
      free_data_copy = 1;
      for (i = 0; i < (N & 3); ++i)
        data_copy[i] = data[i * stride];
      for (; i < N; i += 4) {
        mpb_real d0 = data[i * stride];
        mpb_real d1 = data[(i + 1) * stride];
        mpb_real d2 = data[(i + 2) * stride];
        mpb_real d3 = data[(i + 3) * stride];
        data_copy[i] = d0;
        data_copy[i+1] = d1;
        data_copy[i+2] = d2;
        data_copy[i+3] = d3;
      }
      CHECK(i == N, "bug in matrixio copy routine");
    }
    else {
      data_copy = data;
      data_copy_stride = stride;
    }
  }
  else
    data_copy = data;

  /*******************************************************************/
  /* Before we can write the data to the data set, we must define
     the dimensions and "selections" of the arrays to be read & written: */

  start = (start_t *)malloc(sizeof(start_t) * rank);
  strides = (hsize_t *)malloc(sizeof(hsize_t) * rank);
  count = (hsize_t *)malloc(sizeof(hsize_t) * rank);

  count_prod = 1;
  for (i = 0; i < rank; ++i) {
    start[i] = local_start[i];
    count[i] = local_dims[i];
    strides[i] = 1;
    count_prod *= count[i];
  }

  if (count_prod > 0) {
    H5Sselect_hyperslab(space_id, H5S_SELECT_SET, start, NULL, count, NULL);

    for (i = 0; i < rank; ++i)
         start[i] = 0;
    strides[rank - 1] = data_copy_stride;
    count[rank - 1] *= data_copy_stride;
    mem_space_id = H5Screate_simple(rank, count, NULL);
    count[rank - 1] = local_dims[rank - 1];
    H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET, start, data_copy_stride <= 1 ? NULL : strides, count, NULL);
  }
  else { /* this can happen on leftover processes in MPI */
    H5Sselect_none(space_id);
    mem_space_id = H5Scopy(space_id); /* can't create an empty space */
    H5Sselect_none(mem_space_id);
    do_write = 0;  /* HDF5 complains about empty dataspaces otherwise */
  }

  /*******************************************************************/
  /* Write the data, then free all the stuff we've allocated. */

  if (do_write)
    H5Dwrite(data_id.id, type_id, mem_space_id, space_id, H5P_DEFAULT, data_copy);

  if (free_data_copy)
    free(data_copy);
  H5Sclose(mem_space_id);
  free(count);
  free(strides);
  free(start);
  free(dims);
  H5Sclose(space_id);
#endif
}

void matrixio_close_dataset(matrixio_id data_id) {
#if defined(HAVE_HDF5)
  CHECK(H5Dclose(data_id.id) >= 0, "error closing HDF dataset");
#endif
}

void matrixio_write_data_attr(matrixio_id id, const char *name, const mpb_real *val,
                              int rank, const int *dims) {
#if defined(HAVE_HDF5)
  hid_t type_id;
  hid_t space_id;
  hsize_t *space_dims;
  int i;

  if (!val || !name || !name[0] || rank < 0 || !dims)
    return; /* don't try to create empty attributes */

#if defined(SCALAR_SINGLE_PREC)
  type_id = H5T_NATIVE_FLOAT;
#elif defined(SCALAR_LONG_DOUBLE_PREC)
  type_id = H5T_NATIVE_LDOUBLE;
#else
  type_id = H5T_NATIVE_DOUBLE;
#endif

  if (rank > 0) {
    space_dims = (hsize_t *)malloc(sizeof(hsize_t) * rank);
    for (i = 0; i < rank; ++i)
      space_dims[i] = dims[i];
    space_id = H5Screate_simple(rank, space_dims, NULL);
    free(space_dims);
  }
  else {
    space_id = H5Screate(H5S_SCALAR);
  }

  write_attr(id, type_id, space_id, name, val);
  H5Sclose(space_id);
#endif
}

void matrixio_write_string_attr(matrixio_id id, const char *name, const char *val) {

#if defined(HAVE_HDF5)
  hid_t type_id;
  hid_t space_id;

  if (!val || !name || !name[0] || !val[0])
    return; /* don't try to create empty attributes */

  type_id = H5Tcopy(H5T_C_S1);
  H5Tset_size(type_id, strlen(val) + 1);
  space_id = H5Screate(H5S_SCALAR);
  write_attr(id, type_id, space_id, name, val);
  H5Sclose(space_id);
  H5Tclose(type_id);
#endif
}

matrixio_id matrixio_open_(const char *fname, int read_only, int parallel) {
#if defined(HAVE_HDF5)
  char *new_fname;
  matrixio_id id;
  hid_t access_props;

  access_props = H5Pcreate (H5P_FILE_ACCESS);

#if defined(HAVE_MPI) && defined(HAVE_H5PSET_FAPL_MPIO)
  if (parallel)
    H5Pset_fapl_mpio(access_props, mpb_comm, MPI_INFO_NULL);
#endif

  new_fname = add_fname_suffix(fname);

  IF_EXCLUSIVE(if (parallel) mpi_begin_critical_section(matrixio_critical_section_tag),0);

  if (read_only)
    id.id = H5Fopen(new_fname, H5F_ACC_RDONLY, access_props);
  else
    id.id = H5Fopen(new_fname, H5F_ACC_RDWR, access_props);
  id.parallel = parallel;
  CHECK(id.id >= 0, "error opening HDF input file");

  free(new_fname);

  H5Pclose(access_props);

  return id;
#else
  CHECK(0, "no matrixio implementation is linked");
  {
    matrixio_id id = {0,0};
    return id;
  }
#endif
}

matrixio_id matrixio_open(const char *fname, int read_only) {
  return matrixio_open_(fname, read_only, 1);
}

void matrixio_close(matrixio_id id) {
#if defined(HAVE_HDF5)
  CHECK(H5Fclose(id.id) >= 0, "error closing HDF file");
  IF_EXCLUSIVE(if (id.parallel) mpi_end_critical_section(matrixio_critical_section_tag++),0);
#endif
}

void evectmatrixio_readall_raw(const char *filename, evectmatrix a) {
  int rank = 4, dims[4];
  matrixio_id file_id;

  dims[0] = a.N;
  dims[1] = a.c;
  dims[2] = a.p;
  dims[3] = SCALAR_NUMVALS;

  file_id = matrixio_open(filename, 1);

  CHECK(matrixio_read_real_data(file_id, "rawdata", &rank, dims, a.localN, a.Nstart, 1,
                                (mpb_real *) a.data), "error reading data set in file");

  matrixio_close(file_id);
}

char *add_fname_suffix(const char *fname) {
  int oldlen = strlen(fname);
  int suflen = strlen(FNAME_SUFFIX);

  CHECK(fname, "null filename!");

  char *new_fname = (char *)malloc(sizeof(char) * (oldlen + suflen + 1));

  strcpy(new_fname, fname);

  /* only add suffix if it is not already there: */
  if (strstr(new_fname, FNAME_SUFFIX) != new_fname + oldlen - suflen)
    strcat(new_fname, FNAME_SUFFIX);

  return new_fname;
}

int matrixio_dataset_exists(matrixio_id id, const char *name) {
#if defined(HAVE_HDF5)
  hid_t data_id;
  SUPPRESS_HDF5_ERRORS(data_id = H5Dopen(id.id, name));
  if (data_id >= 0)
    H5Dclose(data_id);
  return (data_id >= 0);
#else
  return 0;
#endif
}

void matrixio_dataset_delete(matrixio_id id, const char *name) {
#if defined(HAVE_HDF5)
  H5Gunlink(id.id, name);
#endif
}

/* Wrappers to write/read an attribute attached to id.  HDF5 attributes
   can *not* be attached to files, in which case we'll write/read it
   as an ordinary dataset.  Ugh. */

void write_attr(matrixio_id id, matrixio_id_ type_id, matrixio_id_ space_id,
                const char *name, const void *val) {

#if defined(HAVE_HDF5)
  hid_t attr_id;

#ifndef HAVE_H5PSET_FAPL_MPIO
  if (!mpi_is_master() && id.parallel)
    return; /* only one process should add attributes */
#else
  /* otherwise, the operations must be performed collectively */
#endif

  if (H5I_FILE == H5Iget_type(id.id)) {
    attr_id = H5Dcreate(id.id, name, type_id, space_id, H5P_DEFAULT);
    CHECK(attr_id >= 0, "error creating HDF attr");
    H5Dwrite(attr_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);
    H5Dclose(attr_id);
  }
  else {
    attr_id = H5Acreate(id.id, name, type_id, space_id, H5P_DEFAULT);
    CHECK(attr_id >= 0, "error creating HDF attr");
    H5Awrite(attr_id, type_id, val);
    H5Aclose(attr_id);
  }
#endif
}

/****** end matrixio ******/