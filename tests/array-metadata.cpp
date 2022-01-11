/***************************************************************/
/***************************************************************/
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
using std::vector;

vector3 v3(double x = 0.0, double y = 0.0, double z = 0.0) {
  vector3 v;
  v.x = x;
  v.y = y;
  v.z = z;
  return v;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static ivec vec2diel_floor(const vec &pt, double a, const ivec &equal_shift) {
  ivec ipt(pt.dim);
  LOOP_OVER_DIRECTIONS(pt.dim, d) {
    ipt.set_direction(d, 1 + 2 * int(floor(pt.in_direction(d) * a - .5)));
    if (ipt.in_direction(d) == pt.in_direction(d))
      ipt.set_direction(d, ipt.in_direction(d) + equal_shift.in_direction(d));
  }
  return ipt;
}
static ivec vec2diel_ceil(const vec &pt, double a, const ivec &equal_shift) {
  ivec ipt(pt.dim);
  LOOP_OVER_DIRECTIONS(pt.dim, d) {
    ipt.set_direction(d, 1 + 2 * int(ceil(pt.in_direction(d) * a - .5)));
    if (ipt.in_direction(d) == pt.in_direction(d))
      ipt.set_direction(d, ipt.in_direction(d) + equal_shift.in_direction(d));
  }
  return ipt;
}
namespace meep {
void compute_boundary_weights(grid_volume gv, const volume &wherec, ivec &is, ivec &ie,
                              bool snap_empty_dims, vec &s0, vec &e0, vec &s1, vec &e1);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool equal_float(double *array1, double *array2, int N) {
  for (int n = 0; n < N; n++)
    if (((float)array1[n]) != ((float)array2[n])) return false;
  return true;
}

/***************************************************************/
/* check that the coordinates and weights computed from the    */
/* metadata match the correct values for all grid points in    */
/* where. return true for zero mismatches, false otherwise.    */
/*                                                             */
/* if the environment variable MEEP_ARRAY_METADATA_LOGFILE is  */
/* set, more detailed output is written to that file.          */
/***************************************************************/
bool test_array_metadata(meep::fields &f, const volume &where) {
  /***************************************************************/
  /* step 1: get coordinate grids and weights as reported by     */
  /* get_array_metadata                                          */
  /***************************************************************/
  size_t dims[3];
  direction dirs[3];
  int rank = f.get_array_slice_dimensions(where, dims, dirs);
  std::vector<double> xyzw = f.get_array_metadata(where);

  // convert to a more convenient format
  int offset = 0;
  size_t nxyz[3], nw = 1;
  vector<double> tics[3], weights;
  for (int i = 0; i < 3; ++i) {
    nxyz[i] = (size_t)xyzw[offset++];
    nw *= nxyz[i];
    for (size_t j = 0; j < nxyz[i]; ++j)
      tics[i].push_back(xyzw[offset++]);
  }
  for (size_t j = 0; j < nw; ++j)
    weights.push_back(xyzw[offset++]);

  size_t stride[3];
  stride[2] = 1;
  stride[1] = nxyz[2];
  stride[0] = nxyz[1] * nxyz[2];

  printf("Metadata: Rank=%i, dims=", rank);
  for (int r = 0; r < rank; r++)
    printf("%c %zu", r == 0 ? '{' : ',', dims[r]);
  printf("}, ");
  printf("xyz sizes={%zu, %zu, %zu}, ", nxyz[0], nxyz[1], nxyz[2]);
  printf("strides={%zu, %zu, %zu}\n", stride[0], stride[1], stride[2]);

  /***************************************************************/
  /* step 2: initialize loop over grid points in the volume via  */
  /*         standard libmeep looping primitives                 */
  /***************************************************************/
  component cgrid = Centered;
  grid_volume gv = f.gv;
  vec yee_c(gv.yee_shift(Centered) - gv.yee_shift(cgrid));
  ivec iyee_c(gv.iyee_shift(Centered) - gv.iyee_shift(cgrid));
  volume wherec(where + yee_c);
  ivec is(vec2diel_floor(wherec.get_min_corner(), gv.a, zero_ivec(gv.dim)));
  ivec ie(vec2diel_ceil(wherec.get_max_corner(), gv.a, zero_ivec(gv.dim)));

  ivec imin = gv.little_corner() + one_ivec(gv.dim), imax = gv.big_corner() - one_ivec(gv.dim);
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    if (is.in_direction(d) < imin.in_direction(d)) is.set_direction(d, imin.in_direction(d));
    if (ie.in_direction(d) > imax.in_direction(d)) ie.set_direction(d, imax.in_direction(d));
  }

  bool snap_empty_dims = true;
  vec s0(gv.dim), e0(gv.dim), s1(gv.dim), e1(gv.dim);
  // this initialization step seems to be necessary here to avoid winding
  // up with zero or undefined integration weights; I don't know why it
  // seems to be unnecessary for loop_in_chunks above.
  FOR_DIRECTIONS(d)
  if (!has_direction(gv.dim, d)) {
    s0.set_direction(d, 1.0);
    e0.set_direction(d, 1.0);
    s1.set_direction(d, 1.0);
    e1.set_direction(d, 1.0);
  }
  compute_boundary_weights(gv, wherec, is, ie, snap_empty_dims, s0, e0, s1, e1);

  // Determine integration "volumes" dV0 and dV1
  double dV0 = 1.0, dV1 = 0.0;
  LOOP_OVER_DIRECTIONS(gv.dim, d)
  if (wherec.in_direction(d) > 0.0) dV0 *= gv.inva;

  /***************************************************************/
  /* step 3: execute the loop and check that coordinates and     */
  /*         weights of each point as determined from the return */
  /*         values of get_array_metadata agree with those       */
  /*         determined by the libmeep loop primitives           */
  /***************************************************************/
  int num_points = 0, num_mismatches = 0;
  char *LogFileName = getenv("MEEP_ARRAY_METADATA_LOGFILE");
  FILE *LogFile = (LogFileName ? fopen(LogFileName, "w") : 0);
  LOOP_OVER_IVECS(gv, is, ie, idx) {
    // get the (correct) coordinates and weight for the current grid point,
    // or (for collapsed dimensions) the sum of the weights of the two
    // points from which we interpolate to get values at the array slice coordinate
    double xyzw_loop[4] = {0.0, 0.0, 0.0, 0.0};
    IVEC_LOOP_LOC(gv, loc);
    xyzw_loop[0] = has_direction(gv.dim, X) ? loc.x() : 0.0;
    xyzw_loop[1] = has_direction(gv.dim, Y) ? loc.y() : 0.0;
    xyzw_loop[2] = has_direction(gv.dim, Z) ? loc.z() : 0.0;
    xyzw_loop[3] = IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2);

    // coordinates and weight for current grid point according to metadata
    double xyzw_meta[4] = {0.0, 0.0, 0.0, 0.0};
    IVEC_LOOP_ILOC(gv, iloc);
    ivec two_n = iloc - is;
    int nx = 0, ny = 0, nz = 0, index = 0;
    if (has_direction(gv.dim, X)) {
      nx = two_n.in_direction(X) / 2;
      xyzw_meta[0] = tics[0][nx];
      index += nx * stride[0];
    }
    if (has_direction(gv.dim, Y)) {
      ny = two_n.in_direction(Y) / 2;
      xyzw_meta[1] = tics[1][ny];
      index += ny * stride[1];
    }
    if (has_direction(gv.dim, Z)) {
      nz = two_n.in_direction(Z) / 2;
      xyzw_meta[2] = tics[2][nz];
      index += nz * stride[2];
    }
    xyzw_meta[3] = weights[index];

    bool mismatch = !equal_float(xyzw_loop, xyzw_meta, 4);
    if (mismatch) num_mismatches++;
    if (LogFile) {
      fprintf(LogFile, "%i %i ", num_points++, mismatch ? 0 : 1);
      fprintf(LogFile, "%e %e %e %e ", xyzw_loop[0], xyzw_loop[1], xyzw_loop[2], xyzw_loop[3]);
      fprintf(LogFile, "%e %e %e %e ", xyzw_meta[0], xyzw_meta[1], xyzw_meta[2], xyzw_meta[3]);
      fprintf(LogFile, "\n");
    }

  } // LOOP_OVER_IVECS(gv, is, ie, idx)
  if (LogFile) fclose(LogFile);
  printf("%i/%i mismatches\n", num_mismatches, num_points);
  return (num_mismatches == 0);
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
int main(int argc, char *argv[]) {
  initialize mpi(argc, argv);

  /*--------------------------------------------------------------*/
  /*- set default geometric parameters ---------------------------*/
  /*--------------------------------------------------------------*/
  // size of computational cell
  double sx = 10.0;
  double sy = 5.0;
  double sz = 0.0;

  // corners of array volume
  double vxmin = -2.5, vxmax = -2.5;
  double vymin = -1.0, vymax = +3.0;
  double vzmin = 0.0, vzmax = 0.0;

  double res = 10.0;

  // double-valued command-line parameters
  vector<const char *> parm_name;
  vector<double *> parm_adrs;
  parm_name.push_back("--sx");
  parm_adrs.push_back(&sx);
  parm_name.push_back("--sy");
  parm_adrs.push_back(&sy);
  parm_name.push_back("--sz");
  parm_adrs.push_back(&sz);
  parm_name.push_back("--vxmin");
  parm_adrs.push_back(&vxmin);
  parm_name.push_back("--vymin");
  parm_adrs.push_back(&vymin);
  parm_name.push_back("--vzmin");
  parm_adrs.push_back(&vzmin);
  parm_name.push_back("--vxmax");
  parm_adrs.push_back(&vxmax);
  parm_name.push_back("--vymax");
  parm_adrs.push_back(&vymax);
  parm_name.push_back("--vzmax");
  parm_adrs.push_back(&vzmax);
  parm_name.push_back("--res");
  parm_adrs.push_back(&res);

  /*--------------------------------------------------------------*/
  /*- parse arguments --------------------------------------------*/
  /*--------------------------------------------------------------*/
  for (int narg = 1; narg < argc; narg++) {
    // process double-valued parameters
    size_t np;
    for (np = 0; np < parm_name.size(); np++)
      if (!strcasecmp(argv[narg], parm_name[np])) break;
    if (np == parm_name.size()) meep::abort("unknown command-line option %s", argv[narg]);
    if (narg + 1 == argc) meep::abort("no option specified for %s", argv[narg]);
    if (1 != sscanf(argv[narg + 1], "%le", parm_adrs[np]))
      meep::abort("invalid value %s specified for %s", argv[narg + 1], argv[narg]);
    printf("Setting %s=%e.\n", argv[narg], *(parm_adrs[np]));
    narg++;
  }

  /*--------------------------------------------------------------*/
  /*- initialize geometry ----------------------------------------*/
  /*--------------------------------------------------------------*/
  geometry_lattice.size.x = sx;
  geometry_lattice.size.y = sy;
  geometry_lattice.size.z = sz;
  grid_volume gv;
  if (sx == 0.0 && sy == 0.0)
    gv = vol1d(sz, res);
  else if (sz == 0.0)
    gv = vol2d(sx, sy, res);
  else
    gv = vol3d(sx, sy, sz, res);
  gv.center_origin();
  structure the_structure(gv, dummy_eps);

  meep_geom::material_type silicon = meep_geom::make_dielectric(12.0);
  geometric_object objects[1];
  vector3 origin = v3(0.0, 0.0, 0.0);
  vector3 wvg_size = v3(0.5 * sx, 0.5 * sy, 0.5 * sz);
  vector3 xhat = {1.0, 0.0, 0.0};
  vector3 yhat = {0.0, 1.0, 0.0};
  vector3 zhat = {0.0, 0.0, 1.0};
  objects[0] = make_block(silicon, origin, xhat, yhat, zhat, wvg_size);
  geometric_object_list g = {1, objects};
  meep_geom::set_materials_from_geometry(&the_structure, g);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  fields f(&the_structure);
  double fcen = 0.2;
  double df = 0.1;
  gaussian_src_time src(fcen, df);
  component src_cmpt = (gv.dim == D1 ? Ex : Ez);
  vec src_point = zero_vec(gv.dim);
  vec src_size = zero_vec(gv.dim);
  f.add_point_source(src_cmpt, src, src_point);

  vec vmin = zero_vec(gv.dim), vmax = zero_vec(gv.dim);
  if (has_direction(gv.dim, X)) {
    vmin.set_direction(X, vxmin);
    vmax.set_direction(X, vxmax);
  }
  if (has_direction(gv.dim, Y)) {
    vmin.set_direction(Y, vymin);
    vmax.set_direction(Y, vymax);
  }
  if (has_direction(gv.dim, Z)) {
    vmin.set_direction(Z, vzmin);
    vmax.set_direction(Z, vzmax);
  }
  volume slice(vmin, vmax);
  bool test_passed = test_array_metadata(f, slice);
  return test_passed ? 0 : 1;
}
