/***************************************************************/
/* simple test for libmeepgeom, modeled after meep_test.ctl    */
/***************************************************************/
#include <memory>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <complex>

#include "meep.hpp"

#include "ctl-math.h"
#include "ctlgeom.h"

#include "meepgeom.hpp"

#ifndef DATADIR
#define DATADIR "./"
#endif

using namespace meep;

typedef std::complex<double> cdouble;

/***************************************************************/
/* return true if the datasets match, false if not             */
/***************************************************************/
bool compare_hdf5_datasets(const char *file1, const char *name1, const char *file2,
                           const char *name2, int expected_rank = 2, double rel_tol = 1.0e-4,
                           double abs_tol = sizeof(realnum) == sizeof(float) ? 1.0e-6 : 1.0e-8) {
  h5file f1(file1, h5file::READONLY, false);
  int rank1;
  std::vector<size_t> dims1(expected_rank);

  std::unique_ptr<realnum[]> data1;
  std::unique_ptr<realnum[]> data2;
  data1.reset(static_cast<realnum *>(
      f1.read(name1, &rank1, dims1.data(), expected_rank, sizeof(realnum) == sizeof(float))));
  if (!data1) return false;

  h5file f2(file2, h5file::READONLY, false);
  int rank2;
  std::vector<size_t> dims2(expected_rank);
  data2.reset(static_cast<realnum *>(
      f2.read(name2, &rank2, dims2.data(), expected_rank, sizeof(realnum) == sizeof(float))));
  if (!data2) return false;

  if (rank1 != expected_rank || rank2 != expected_rank) return false;

  size_t size = 1;
  for (int r = 0; r < expected_rank; r++) {
    if (dims1[r] != dims2[r]) return false;
    size *= dims1[r];
  };

  for (size_t n = 0; n < size; n++) {
    double d1 = data1[n], d2 = data2[n], diff = fabs(d1 - d2), max = fmax(fabs(d1), fabs(d2));
    if (diff > abs_tol || diff > max * rel_tol) return false;
  };

  return true;
}

/***************************************************************/
/* dummy material function needed to pass to structure( )      */
/* constructor as a placeholder before we can call             */
/* set_materials_from_geometry                                 */
/***************************************************************/
double dummy_eps(const vec &) { return 1.0; }

/***************************************************************/
/* usage: cyl-ellipsoid [ --polarization xx ]                  */
/*                      [ --eps_ref_file path/to/file ]        */
/*  where xx = S for TE polarization (default)                 */
/*             P for TM polarization                           */
/***************************************************************/
int main(int argc, char *argv[]) {
  initialize mpi(argc, argv);

  // simple argument parsing
  meep::component src_cmpt = Ez;
  std::string eps_ref_file = "cyl-ellipsoid-eps-ref.h5";
  for (int narg = 1; narg < argc; narg++) {
    if (argv[narg] && !strcmp(argv[narg], "--polarization")) {
      if (narg + 1 == argc)
        meep::abort("no option specified for --polarization");
      else if (!strcasecmp(argv[narg + 1], "S"))
        master_printf("Using S-polarization\n");
      else if (!strcasecmp(argv[narg + 1], "P")) {
        src_cmpt = Hz;
        master_printf("Using P-polarization\n");
      }
      else
        meep::abort("invalid --polarization %s", argv[narg + 1]);
    }
    else if (argv[narg] && !strcmp(argv[narg], "--eps_ref_file")) {
      if (narg + 1 == argc) meep::abort("no option specified for --eps_ref_file");
      eps_ref_file = argv[++narg];
    }
    else
      meep::abort("unrecognized command-line option %s", argv[narg]);
  };

  std::string eps_ref_path = DATADIR + eps_ref_file;

  //(set-param! resolution 100)
  double resolution = 100.0;

  // (set! geometry-lattice (make lattice (size 10 10 no-size)))
  // (set! pml-layers (list (make pml (thickness 1))))
  // (if (= src-cmpt Ez)
  //  (set! symmetries (list (make mirror-sym (direction X))
  //                         (make mirror-sym (direction Y)))))
  // (if (= src-cmpt Hz)
  //  (set! symmetries (list (make mirror-sym (direction X) (phase -1))
  //  (set! symmetries (list (make mirror-sym (direction Y) (phase -1))
  geometry_lattice.size.x = 10.0;
  geometry_lattice.size.y = 10.0;
  geometry_lattice.size.z = 0.0;
  grid_volume gv = voltwo(10.0, 10.0, resolution);
  gv.center_origin();
  symmetry sym = (src_cmpt == Ez) ? mirror(X, gv) + mirror(Y, gv) : -mirror(X, gv) - mirror(Y, gv);
  structure the_structure(gv, dummy_eps, pml(1.0), sym);

  // (set! geometry (list
  //    (make cylinder (center 0 0 0) (radius 3) (height infinity)
  //                   (material (make medium (index 3.5))))
  //    (make ellipsoid (center 0 0 0) (size 1 2 infinity)
  //                   (material air))))
  double n = 3.5; // index of refraction
  auto material_deleter = [](meep_geom::material_data *m) { meep_geom::material_free(m); };
  std::unique_ptr<meep_geom::material_data, decltype(material_deleter)> dielectric(
      meep_geom::make_dielectric(n * n), material_deleter);
  geometric_object objects[2];
  vector3 center = {0.0, 0.0, 0.0};
  double radius = 3.0;
  double height = meep_geom::ENORMOUS;
  vector3 xhat = {1.0, 0.0, 0.0};
  vector3 yhat = {0.0, 1.0, 0.0};
  vector3 zhat = {0.0, 0.0, 1.0};
  vector3 size = {1.0, 2.0, 1.0e20};
  objects[0] = make_cylinder(dielectric.get(), center, radius, height, zhat);
  objects[1] = make_ellipsoid(meep_geom::vacuum, center, xhat, yhat, zhat, size);
  geometric_object_list g = {2, objects};
  meep_geom::set_materials_from_geometry(&the_structure, g);

  // (set! sources (list (make source (src (make gaussian-src (frequency 1) (fwidth 0.1)))
  //  (center 0 0 0) (component src-cmpt))))
  fields f(&the_structure);
  double fcen = 1.0;
  double df = 0.1;
  gaussian_src_time src(fcen, df);
  vec src_point = vec(0.0, 0.0);
  vec src_size = vec(10.0, 10.0);
  f.add_point_source(src_cmpt, src, src_point);

  // first test: write permittivity to HDF5 file and
  // compare with contents of reference file
  if (am_really_master()) {
    f.output_hdf5(Dielectric, f.total_volume());
    bool status = compare_hdf5_datasets("eps-000000000.h5", "eps", eps_ref_path.c_str(), "eps");
    if (status)
      master_printf("Dielectric output test successful.\n");
    else
      meep::abort("Dielectric output error in cyl-ellipsoid-ll");
  };

  // (run-until 23 (at-beginning output-epsilon)
  //	           (at-end output-efield-z)
  //	           (at-end print-stuff))
  double duration = 23.0;
  double start_time = f.round_time();
  double stop_time = start_time + duration;
  while (f.round_time() < stop_time)
    f.step();

  // second test: compare field component at specified evaluation
  // point to reference values
  if (am_really_master()) {
// ref values obtained by running `meep cyl-ellipsoid.ctl`
#define REF_EZ -8.29555720049629e-5
#define REF_HZ -4.5623185899766e-5
#define RELTOL 0.05
    double ref_out_field = (src_cmpt == Ez) ? REF_EZ : REF_HZ;
    double out_field = real(f.get_field(src_cmpt, vec(4.13, 3.75)));
    double diff = fabs(out_field - ref_out_field);
    printf("field: %e\n", out_field);

    if (fabs(diff) <= RELTOL * fabs(ref_out_field))
      printf("Field component output test successful.");
    else
      meep::abort("field output error in cyl-ellipsoid-ll");
  };

  for (int n = 0; n < 2; n++) {
    geometric_object_destroy(objects[n]);
  }
  meep_geom::unset_default_material();

  return 0;
}
