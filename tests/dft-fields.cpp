/***************************************************************/
/* unit test for dft_fields and output_dft functionality       */
/***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <vector>

#include "meep.hpp"

#include "ctl-math.h"
#include "ctlgeom.h"

#include "meepgeom.hpp"

using namespace meep;

vector3 v3(double x, double y = 0.0, double z = 0.0) {
  vector3 v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}

/***************************************************************/
/* dummy material function needed to pass to structure( )      */
/* constructor as a placeholder before we can call             */
/* set_materials_from_geometry                                 */
/***************************************************************/
double dummy_eps(const vec &) { return 1.0; }

/***************************************************************/
/***************************************************************/
/***************************************************************/
void Run(bool Pulse, double resolution, std::complex<meep::realnum> **field_array = 0,
         int *array_rank = 0, size_t *array_dims = 0) {
  /***************************************************************/
  /* initialize geometry                                         */
  /***************************************************************/
  double n = 3.4;  // index of waveguide
  double w = 1.0;  // width of waveguide
  double r = 1.0;  // inner radius of ring
  double pad = 4;  // padding between waveguide and edge of PML
  double dpml = 2; // thickness of PML

  double sxy = 2.0 * (r + w + pad + dpml); // cell size
  geometry_lattice.size.x = sxy;
  geometry_lattice.size.y = sxy;
  geometry_lattice.size.z = 0.0;
  grid_volume gv = voltwo(sxy, sxy, resolution);
  gv.center_origin();
  symmetry sym = identity(); // mirror(Y, gv);
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  /***************************************************************/
  /* add objects                                                 */
  /***************************************************************/
  meep_geom::material_type dielectric = meep_geom::make_dielectric(n * n);
  geometric_object objects[2];
  vector3 v3zero = {0.0, 0.0, 0.0};
  vector3 zaxis = {0.0, 0.0, 1.0};
  objects[0] = make_cylinder(dielectric, v3zero, r + w, meep_geom::ENORMOUS, zaxis);
  objects[1] = make_cylinder(meep_geom::vacuum, v3zero, r, meep_geom::ENORMOUS, zaxis);
  geometric_object_list g = {2, objects};
  meep_geom::set_materials_from_geometry(&the_structure, g);
  fields f(&the_structure);
  f.step(); // single timestep to trigger internal initialization

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double fcen = 0.118;  // ; pulse center frequency
  double df = 0.1;      // ; df
  vec x0(r + 0.1, 0.0); // ; source location
  if (Pulse) {
    f.add_point_source(Ez, gaussian_src_time(fcen, df), x0);

    component components[6] = {Ex, Ey, Ez, Hx, Hy, Hz};
    dft_fields dftFields = f.add_dft_fields(components, 6, f.v, fcen, fcen, 1);
    dft_flux dftFlux = f.add_dft_flux(X, f.v, fcen, fcen, 1);

    while (f.round_time() < f.last_source_time() + 100.0)
      f.step();

    f.output_dft(dftFlux, "dft-flux");
    f.output_dft(dftFields, "dft-fields");

    *field_array = f.get_dft_array(dftFlux, Ez, 0, array_rank, array_dims);
  }
  else {
    f.add_point_source(Ez, continuous_src_time(fcen, df), x0);
    f.solve_cw(sizeof(realnum) == sizeof(float) ? 1e-5 : 1e-8, 10000, 10);
    h5file *file = f.open_h5file("cw-fields", h5file::WRITE, 0, false);
    f.output_hdf5(Ez, f.v, file);
    f.output_hdf5(Hx, f.v, file);
    f.output_hdf5(Hy, f.v, file);
    delete file;
  }
}

/***************************************************************/
/* return L2 norm of error normalized by average of L2 norms   */
/***************************************************************/
double compare_array_to_dataset(std::complex<meep::realnum> *field_array, int array_rank,
                                size_t *array_dims, const char *file, const char *name) {
  int file_rank;
  size_t file_dims[3];
  h5file f(file, h5file::READONLY, false);
  char dataname[100];
  snprintf(dataname, 100, "%s.r", name);
  double *rdata =
      (double *)f.read(dataname, &file_rank, file_dims, 2, false /* single_precision */);
  snprintf(dataname, 100, "%s.i", name);
  double *idata =
      (double *)f.read(dataname, &file_rank, file_dims, 2, false /* single_precision */);
  if (!rdata || !idata) return -1.0;
  if (file_rank != array_rank) return -1.0;
  for (int n = 0; n < file_rank; n++)
    if (file_dims[n] != array_dims[n]) return -1.0;

  double NormArray = 0.0, NormFile = 0.0, NormDelta = 0.0;
  for (size_t n = 0; n < file_dims[0] * file_dims[1]; n++) {
    std::complex<double> zArray = field_array[n];
    std::complex<double> zFile = std::complex<double>(rdata[n], idata[n]);
    NormArray += norm(zArray);
    NormFile += norm(zFile);
    NormDelta += norm(zArray - zFile);
  }
  NormArray = sqrt(NormArray);
  NormFile = sqrt(NormFile);
  NormDelta = sqrt(NormDelta);
  double RelErr = NormDelta / (0.5 * (NormArray + NormFile));
  return RelErr;
}

/***************************************************************/
/* compute the L2 norm of two complex-valued HDF5 datasets,    */
/* after (a) normalizing each data set by its maximum amplitude*/
/* and (b) compensating for a constant overall phase factor    */
/* between the datasets.                                       */
/***************************************************************/
double compare_complex_hdf5_datasets(const char *file1, const char *name1, const char *file2,
                                     const char *name2, int expected_rank = 2,
                                     double *max_dft = 0) {
  char dataname[100];

  // read dataset 1
  h5file f1(file1, h5file::READONLY, false);
  int rank1;
  size_t *dims1 = new size_t[expected_rank];
  snprintf(dataname, 100, "%s.r", name1);
  double *rdata1 =
      (double *)f1.read(dataname, &rank1, dims1, expected_rank, false /* single_precision */);
  snprintf(dataname, 100, "%s.i", name1);
  double *idata1 =
      (double *)f1.read(dataname, &rank1, dims1, expected_rank, false /* single_precision */);
  if (!rdata1 || !idata1) return -1.0;

  // read dataset 2
  h5file f2(file2, h5file::READONLY, false);
  int rank2;
  size_t *dims2 = new size_t[expected_rank];
  snprintf(dataname, 100, "%s.r", name2);
  double *rdata2 =
      (double *)f2.read(dataname, &rank2, dims2, expected_rank, false /* single_precision */);
  snprintf(dataname, 100, "%s.i", name2);
  double *idata2 =
      (double *)f2.read(dataname, &rank2, dims2, expected_rank, false /* single_precision */);
  if (!rdata2 || !idata2) return -1.0;

  // check same size
  bool same_size = (rank1 == rank2);
  for (int d = 0; same_size && d < rank1; d++)
    if (dims1[d] != dims2[d]) same_size = false;
  if (!same_size) return -1.0;

  // first pass to normalize each dataset to its maximum absolute magnitude;
  // we also note the phase difference between the datasets at their points
  // of maximum magnitude so we can compensate for this in the comparison below.
  size_t length = dims1[0];
  for (int d = 1; d < rank1; d++)
    length *= dims1[d];

  double max_abs1 = 0.0, max_abs2 = 0.0;
  double max_arg1 = 0.0, max_arg2 = 0.0;
  for (size_t n = 0; n < length; n++) {
    std::complex<double> z1 = std::complex<double>(rdata1[n], idata1[n]);
    if (abs(z1) > max_abs1) {
      max_abs1 = abs(z1);
      max_arg1 = arg(z1);
    }
    std::complex<double> z2 = std::complex<double>(rdata2[n], idata2[n]);
    if (abs(z2) > max_abs2) {
      max_abs2 = abs(z2);
      max_arg2 = arg(z2);
    }
  }
  *max_dft = max_abs1;

  // second pass to get L2 norm of difference between normalized data sets
  double norm1 = 0.0, norm2 = 0.0, normdiff = 0.0;
  std::complex<double> phase1 = exp(-std::complex<double>(0, 1) * max_arg1);
  std::complex<double> phase2 = exp(-std::complex<double>(0, 1) * max_arg2);
  for (size_t n = 0; n < length; n++) {
    std::complex<double> z1 = phase1 * std::complex<double>(rdata1[n], idata1[n]) / max_abs1;
    std::complex<double> z2 = phase2 * std::complex<double>(rdata2[n], idata2[n]) / max_abs2;
    norm1 += norm(z1);
    norm2 += norm(z2);
    normdiff += norm(z1 - z2);
  }
  norm1 = sqrt(norm1) / ((double)length);
  norm2 = sqrt(norm2) / ((double)length);
  normdiff = sqrt(normdiff) / ((double)length);

  master_printf("norm (DFT)=%e, norm(CW)=%e, norm(DFT-CT)=%e\n", norm1, norm2, normdiff);

  return normdiff / fmax(norm1, norm2);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) {
  initialize mpi(argc, argv);

  double resolution = 10.0;
  bool verbose = false;
  for (int narg = 1; narg < argc; narg++) {
    if (!strcasecmp(argv[narg], "--resolution")) {
      if (narg + 1 >= argc) meep::abort("--resolution requires an argument");
      sscanf(argv[narg + 1], "%le", &resolution);
      master_printf("Setting resolution=%e.\n", resolution);
      narg++;
    }
    else if (!strcasecmp(argv[narg], "--verbose"))
      verbose = true;
    else
      meep::abort("unknown argument %s", argv[narg]);
  }

  std::complex<meep::realnum> *field_array = 0;
  int array_rank;
  size_t array_dims[3];
  Run(true, resolution, &field_array, &array_rank, array_dims);
  Run(false, resolution);

  /* compare DFT field array to DFT HDF5 output */
  double L2ErrorArray =
      compare_array_to_dataset(field_array, array_rank, array_dims, "dft-fields.h5", "ez_0");
  if (verbose) master_printf("L2Error (array<-->file) = %e\n", L2ErrorArray);

  /* compare DFT fields to CW fields *****************************/
  double max_dft;
  double L2ErrorFile =
      compare_complex_hdf5_datasets("dft-fields.h5", "ez_0", "cw-fields.h5", "ez", 2, &max_dft);
  if (verbose) master_printf("L2Error (file<-->file) = %e\n", L2ErrorFile);

  bool unit_test = (argc == 1); // run unit-test checks if no command-line arguments
  if (unit_test) {
    if (L2ErrorFile == -1.0) // files couldn't be read or datasets had different sizes
    {
      master_printf("failed to compare data files");
      return -1;
    }

#define REF_MAX_DFT 6.740116e+00
    if ((max_dft - REF_MAX_DFT) > 1.0e-5 * REF_MAX_DFT) {
      master_printf("max dft amplitude=%e, should be %e\n", max_dft, REF_MAX_DFT);
      return -1;
    }

    if (L2ErrorFile > 1.0) {
      master_printf("L2 norm of file-file error=%e (should be <1)\n", L2ErrorFile);
      return -1;
    }

    if (L2ErrorArray == -1.0) {
      master_printf("failed to compare array to data file");
      return -1;
    }

    if (L2ErrorArray > 1.0) {
      master_printf("L2 norm of array-file error=%e (should be <1)\n", L2ErrorArray);
      return -1;
    }

    return 0;
  }
}
