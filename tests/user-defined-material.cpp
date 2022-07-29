/***************************************************************/
/* test of user-defined materials in libmeepgeom.              */
/* this code creates a user-defined material (defined by the   */
/* function my_material_func below) that reproduces the        */
/* material geometry of the cyl_ellipsoid example.             */
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

/***************************************************************/
/* geometric parameters ****************************************/
/***************************************************************/
#define R1X 0.5 // ellipsoid minor radius
#define R1Y 1.0 // ellipsoid major radius
#define R2 3.0  // cylinder radius

/***************************************************************/
/* parameters in (two-oscillator) Drude-Lorenz model of Ag     */
/***************************************************************/
#define EV_UM (1.0 / 1.23984193)
#define AG_WP (9.01 * EV_UM)
#define AG_F0 0.845

#define AG_FRQ0 1.0e-10
#define AG_GAM0 0.48 * EV_UM
#define AG_SIG0 (AG_F0 * AG_WP * AG_WP) / (AG_FRQ0 * AG_FRQ0)

#define AG_FRQ1 0.065
#define AG_GAM1 0.816 * EV_UM
#define AG_SIG1 (AG_F0 * AG_WP * AG_WP) / (AG_FRQ1 * AG_FRQ1)

using namespace meep;

typedef std::complex<double> cdouble;

/***************************************************************/
/* return true if the datasets match, false if not             */
/***************************************************************/
bool compare_hdf5_datasets(const char *file1, const char *name1, const char *file2,
                           const char *name2, int expected_rank = 2, double rel_tol = 1.0e-2) {
  h5file f1(file1, h5file::READONLY, false);
  int rank1;
  size_t *dims1 = new size_t[expected_rank];
  realnum *data1 =
      (realnum *)f1.read(name1, &rank1, dims1, expected_rank, sizeof(realnum) == sizeof(float));
  if (!data1) return false;

  h5file f2(file2, h5file::READONLY, false);
  int rank2;
  size_t *dims2 = new size_t[expected_rank];
  realnum *data2 =
      (realnum *)f2.read(name2, &rank2, dims2, expected_rank, sizeof(realnum) == sizeof(float));
  if (!data2) return false;

  if (rank1 != expected_rank || rank2 != expected_rank) return false;

  size_t size = 1;
  for (int r = 0; r < expected_rank; r++) {
    if (dims1[r] != dims2[r]) return false;
    size *= dims1[r];
  };

  double norm1 = 0.0, norm2 = 0.0, normDelta = 0.0;
  for (size_t n = 0; n < size; n++) {
    double d1 = data1[n], d2 = data2[n], delta = d1 - d2;
    norm1 += d1 * d1;
    norm2 += d2 * d2;
    normDelta += delta * delta;
  };
  norm1 = sqrt(norm1) / size;
  norm2 = sqrt(norm2) / size;
  normDelta = sqrt(normDelta) / size;

  double norm12 = fmax(sqrt(norm1), sqrt(norm2));
  if (normDelta > rel_tol * norm12) return false;
  return true;
}

/***************************************************************/
/* dummy material function needed to pass to structure( )      */
/* constructor as a placeholder before we can call             */
/* set_materials_from_geometry                                 */
/***************************************************************/
double dummy_eps(const vec &) { return 1.0; }

/***************************************************************/
/* material function that recreates the ellipsoid-in-cylinder  */
/* configuration of the cyl-ellipsoid sample code              */
/***************************************************************/
typedef struct my_material_func_data {
  double rxInner, ryInner, rOuter;
  bool with_susceptibility;
} my_material_func_data;

void my_material_func(vector3 p, void *user_data, meep_geom::medium_struct *m) {
  my_material_func_data *data = (my_material_func_data *)user_data;
  double rxInner = data->rxInner, rxInner2 = rxInner * rxInner;
  double ryInner = data->ryInner, ryInner2 = ryInner * ryInner;
  double rOuter = data->rOuter, rOuter2 = rOuter * rOuter;

  double x = p.x, x2 = x * x, y = p.y, y2 = y * y;

  // test for point inside inner ellipsoid
  bool innermost = ((x2 / rxInner2 + y2 / ryInner2) < 1.0);
  bool outermost = ((x * x + y * y) > rOuter2);
  bool in_middle = (!innermost && !outermost);

  // set permittivity
  double nn = in_middle ? 3.5 : 1.0;
  m->epsilon_diag.x = m->epsilon_diag.y = m->epsilon_diag.z = nn * nn;

  // add susceptibilities (two-oscillator model for Ag)
  if (in_middle && data->with_susceptibility) {
    m->E_susceptibilities.resize(2);

    m->E_susceptibilities[0].sigma_offdiag.x = 0.0;
    m->E_susceptibilities[0].sigma_offdiag.y = 0.0;
    m->E_susceptibilities[0].sigma_offdiag.z = 0.0;
    m->E_susceptibilities[0].sigma_diag.x = AG_SIG0;
    m->E_susceptibilities[0].sigma_diag.y = AG_SIG0;
    m->E_susceptibilities[0].sigma_diag.z = AG_SIG0;
    m->E_susceptibilities[0].frequency = AG_FRQ0;
    m->E_susceptibilities[0].gamma = AG_GAM0;
    m->E_susceptibilities[0].noise_amp = 0.0;
    m->E_susceptibilities[0].drude = true;
    m->E_susceptibilities[0].is_file = false;

    m->E_susceptibilities[1].sigma_offdiag.x = 0.0;
    m->E_susceptibilities[1].sigma_offdiag.y = 0.0;
    m->E_susceptibilities[1].sigma_offdiag.z = 0.0;
    m->E_susceptibilities[1].sigma_diag.x = AG_SIG1;
    m->E_susceptibilities[1].sigma_diag.y = AG_SIG1;
    m->E_susceptibilities[1].sigma_diag.z = AG_SIG1;
    m->E_susceptibilities[1].frequency = AG_FRQ1;
    m->E_susceptibilities[1].gamma = AG_GAM1;
    m->E_susceptibilities[1].noise_amp = 0.0;
    m->E_susceptibilities[1].drude = true;
    m->E_susceptibilities[1].is_file = false;
  }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) {
  initialize mpi(argc, argv);

  // simple argument parsing
  meep::component src_cmpt = Ez;
  std::string eps_ref_file = "cyl-ellipsoid-eps-ref.h5";
  bool with_susceptibility = true;
  for (int narg = 1; narg < argc; narg++) {
    if (argv[narg] && !strcmp(argv[narg], "--eps_ref_file")) {
      if (narg + 1 == argc) meep::abort("no option specified for --eps_ref_file");
      eps_ref_file = argv[++narg];
    }
    else if (argv[narg] && !strcmp(argv[narg], "--without_susceptibility"))
      with_susceptibility = false;
    else
      meep::abort("unrecognized command-line option %s", argv[narg]);
  };
  std::string eps_ref_path = DATADIR + eps_ref_file;

  double resolution = 100.0;

  geometry_lattice.size.x = 10.0;
  geometry_lattice.size.y = 10.0;
  geometry_lattice.size.z = 0.0;
  grid_volume gv = voltwo(10.0, 10.0, resolution);
  gv.center_origin();
  symmetry sym = (src_cmpt == Ez) ? mirror(X, gv) + mirror(Y, gv) : -mirror(X, gv) - mirror(Y, gv);
  structure the_structure(gv, dummy_eps, pml(1.0), sym);

  my_material_func_data data;
  data.with_susceptibility = with_susceptibility;
  data.rxInner = R1X;
  data.ryInner = R1Y;
  data.rOuter = R2;

  meep_geom::material_type my_material =
      meep_geom::make_user_material(my_material_func, (void *)&data, false);

  geometric_object_list g = {0, 0};
  vector3 center = {0, 0, 0};
  bool use_anisotropic_averaging = true;
  bool ensure_periodicity = true;
  meep_geom::set_materials_from_geometry(&the_structure, g, center, use_anisotropic_averaging,
                                         DEFAULT_SUBPIXEL_TOL, DEFAULT_SUBPIXEL_MAXEVAL,
                                         ensure_periodicity, my_material);

  fields f(&the_structure);

  f.output_hdf5(Dielectric, f.total_volume());
  bool status = compare_hdf5_datasets("eps-000000000.h5", "eps", eps_ref_path.c_str(), "eps");
  if (status)
    master_printf("user-defined-material test successful.\n");
  else
    meep::abort("user-defined-material test failed.\n");
  return 0;
}
