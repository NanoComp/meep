/***************************************************************/
/* test for absorber functionality in libmeepgeom              */
/* modeled after absorber_1d.ctl                               */
/***************************************************************/
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
/* dummy material function needed to pass to structure( )      */
/* constructor as a placeholder before we can call             */
/* set_materials_from_geometry                                 */
/***************************************************************/
double dummy_eps(const vec &) { return 1.0; }

/***************************************************************/
/* usage: absorber-1d [ --pml ]                                */
/***************************************************************/
int main(int argc, char *argv[]) {
  initialize mpi(argc, argv);

  // simple argument parsing
  bool use_pml = false;
  bool verbose = false;
  for (int narg = 1; narg < argc; narg++) {
    if (argv[narg] && !strcmp(argv[narg], "--pml"))
      use_pml = true;
    else if (argv[narg] && !strcmp(argv[narg], "--verbose"))
      verbose = true;
    else
      meep::abort("unrecognized command-line option %s", argv[narg]);
  };
  if (verbose) master_printf("Using %s.\n", use_pml ? "pml" : "absorber");

  double resolution = 20.0;
  geometry_lattice.size.x = meep_geom::TINY;
  geometry_lattice.size.y = meep_geom::TINY;
  geometry_lattice.size.z = 10.0;
  grid_volume gv = volone(10.0, resolution);
  gv.center_origin();
  boundary_region br = use_pml ? pml(1.0) : boundary_region();
  structure the_structure(gv, dummy_eps, br);

  meep_geom::absorber_list alist = 0;
  if (!use_pml) {
    alist = meep_geom::create_absorber_list();
    meep_geom::add_absorbing_layer(alist, 1.0, Z);
  };

  geometric_object_list g = {0, 0};
  vector3 center = {0, 0, 0};
  meep_geom::set_materials_from_geometry(&the_structure, g, center,
                                         true,                     // use_anisotropic_averaging,
                                         DEFAULT_SUBPIXEL_TOL,     // tol
                                         DEFAULT_SUBPIXEL_MAXEVAL, // maxeval
                                         false,                    // ensure_periodicity
                                         meep_geom::vacuum, alist);

  if (alist) meep_geom::destroy_absorber_list(alist);

  // (set! sources (list (make source (src (make gaussian-src (frequency (/ 0.803)) (fwidth 0.1)))
  // (center 0 0 0) (component Ex))))
  fields f(&the_structure);
  double fcen = 1.0 / 0.803;
  double df = 0.1;
  gaussian_src_time src(fcen, df);
  vec x0(0.0);
  f.add_point_source(Ex, src, x0);

  double min_time = 50.0;
  double dt_output = 10.0;
  double next_output = dt_output;
  double start_time = f.time();
  double max_field = real(f.get_field(Ex, x0));
  double field = max_field;
  double f50 = 0.0;
  bool done = false;
  while (!done) {
    f.step();
    double t = f.time();
    field = real(f.get_field(Ex, x0));
    if (t >= 50.0 && f50 == 0.0) f50 = field;
    max_field = fmax(fabs(field), max_field);
    if (t >= next_output) {
      next_output += dt_output;
      if (verbose) master_printf("%.6e: %+.12e \n", t, field);
    };
    done = ((t - start_time) > min_time && (fabs(field) <= 1.0e-6 * fabs(max_field)));
  };
  double tFinal = f.time();
  double fFinal = real(f.get_field(Ex, x0));

  // test: compare {f50, tFinal, fFinal} to
  // their reference values
  double f50_ref = -5.090114e-01;
  double tFinal_ref = 9.195000e+01;
  double fFinal_ref = sizeof(realnum) == sizeof(float) ? 1.413336e-07 : 1.624782e-07;
  if (fabs(f50 - f50_ref) > 1.0e-6 * fabs(f50_ref) ||
      fabs(tFinal - tFinal_ref) > 1.0e-6 * fabs(tFinal_ref) ||
      fabs(fFinal - fFinal_ref) > 1.0e-6 * fabs(fFinal_ref)) {
    master_printf("{f50, tFinal, fFinal}={%e,%e,%e}\n", f50, tFinal, fFinal);
    master_printf(" should be:\n");
    master_printf("{f50, tFinal, fFinal}={%e,%e,%e}\n", f50_ref, tFinal_ref, fFinal_ref);
    meep::abort("Test failed.");
  }
  else if (verbose)
    master_printf("Test successful.\n");

  meep_geom::unset_default_material();

  return 0;
}
