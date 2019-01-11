/***************************************************************/
/***************************************************************/
/***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <complex>

#include "ctl-math.h"
#include "ctlgeom.h"
#include "meep.hpp"
#include "meep/meepgeom.hpp"
#include "meep/vec.hpp"

using namespace meep;
using namespace std;

typedef std::complex<double> cdouble;

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

// Pretty=true  --> {x,y,z}
// Pretty=false -->  x y z
void fprint_vec(FILE *f, vec v, bool Pretty = false) {
  const char *s = Pretty ? "{" : " ";
  LOOP_OVER_DIRECTIONS(v.dim, d) {
    fprintf(f, "%s%e", s, v.in_direction(d));
    if (Pretty) s = ",";
  }
  fprintf(f, "%s", Pretty ? "} " : " ");
}

// Pretty=true  --> {x,y,z}
// Pretty=false -->  x y z
void fprint_ivec(FILE *f, ivec v, bool Pretty = false) {
  const char *s = Pretty ? "{" : " ";
  LOOP_OVER_DIRECTIONS(v.dim, d) {
    fprintf(f, "%s%i", s, v.in_direction(d));
    if (Pretty) s = ",";
  }
  fprintf(f, "%s", Pretty ? "} " : " ");
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void print_chunk_info(fields_chunk *fc, int ichunk, component cgrid, ivec is, ivec ie,
                             vec s0, vec s1, vec e0, vec e1, double dV0, double dV1, ivec shift,
                             complex<double> shift_phase, const symmetry &S, int sn, void *data_) {
  (void)s0;
  (void)s1;
  (void)e0;
  (void)e1;
  (void)dV0;
  (void)dV1;
  (void)shift_phase;
  (void)data_;

  int iproc = fc->n_proc();
  char FileName[100];
  snprintf(FileName, 100, "/tmp/ChunkInfo_%i_%i", iproc, ichunk);
  FILE *f = fopen(FileName, "w");

  // write some info on the chunk
  fprintf(f, "# chunk %i (process %i) \n", ichunk, iproc);
  fprintf(f, "# is={%i,%i,%i}\n", is.in_direction(X), is.in_direction(Y), is.in_direction(Z));
  fprintf(f, "# ie={%i,%i,%i}\n", ie.in_direction(X), ie.in_direction(Y), ie.in_direction(Z));
  fprintf(f, "# Ex=%s\n", component_name(S.transform(Ex, -sn)));
  fprintf(f, "# Ey=%s\n", component_name(S.transform(Ey, -sn)));
  fprintf(f, "# Hx=%s\n", component_name(S.transform(Hx, -sn)));
  fprintf(f, "# Hy=%s\n", component_name(S.transform(Hy, -sn)));
  fprintf(f, "# \n");
  fprintf(f, "# columns below: \n");
  fprintf(f, "# 1,2   process, chunk index\n");
  fprintf(f, "# 3-5   grid-point indices     (physical, i.e. before symmetry)\n");
  fprintf(f, "# 6-8   grid-point coordinates (physical, i.e. before symmetry)\n");
  fprintf(f, "# 9-11  grid-point indices     (logical,  i.e. after symmetry)\n");
  fprintf(f, "# 12-14 grid-point coordinates (logical,  i.e. after symmetry)\n");

  // loop over all grid points
  vec rshift(shift * (0.5 * fc->gv.inva));
  LOOP_OVER_IVECS(fc->gv, is, ie, idx) {
    fprintf(f, "%i %i ", iproc, ichunk);

    IVEC_LOOP_ILOC(fc->gv, iloc);
    fprint_ivec(f, iloc);
    iloc = S.transform(iloc, sn) + shift;
    fprint_ivec(f, iloc);

    IVEC_LOOP_LOC(fc->gv, loc);
    fprint_vec(f, loc);
    loc = S.transform(loc, sn) + rshift;
    fprint_vec(f, loc);

    fprintf(f, "\n");
  }
  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
structure create_structure(double resolution, int symmetries) {
  double sx = 8.0;   // size of cell in X direction
  double sy = 6.0;   // size of cell in Y direction
  double w = 3.0;    // width of waveguide
  double dpml = 1.0; // PML thickness

  geometry_lattice.size.x = sx;
  geometry_lattice.size.y = sy;
  geometry_lattice.size.z = 0.0;
  grid_volume gv = voltwo(sx, sy, resolution);
  gv.center_origin();
  symmetry S = (symmetries == 1 ? mirror(Y, gv)
                                : symmetries == 2 ? mirror(X, gv) + mirror(Y, gv) : symmetry());
  structure the_structure(gv, dummy_eps, pml(dpml), S);

  vector3 e1 = v3(1.0, 0.0, 0.0);
  vector3 e2 = v3(0.0, 1.0, 0.0);
  vector3 e3 = v3(0.0, 0.0, 1.0);
  meep_geom::material_type dielectric = meep_geom::make_dielectric(12.0);

  vector3 center = v3(0.0, 0.0, 0.0);
  vector3 size = v3(sx, w, 0.0);
  geometric_object objects[1];
  objects[0] = make_block(dielectric, center, e1, e2, e3, size);
  geometric_object_list g = {1, objects};
  meep_geom::set_materials_from_geometry(&the_structure, g);

  return the_structure;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) {
  initialize mpi(argc, argv);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double resolution = 5.0;
  int symmetries = 0; // {0,1,2} for none, Y-mirror, XY-mirror
  for (int narg = 1; narg < argc; narg++)
    if (!strcasecmp(argv[narg], "--resolution")) sscanf(argv[narg + 1], "%le", &resolution);
  for (int narg = 1; narg < argc - 1; narg++)
    if (!strcasecmp(argv[narg], "--symmetries")) sscanf(argv[narg + 1], "%i", &symmetries);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  structure s = create_structure(resolution, symmetries);
  fields f(&s);
  f.step();
  const char *prefix = "WriteChunkInfo";
  // f.output_hdf5(Dielectric, f.total_volume(), 0, false, true, prefix);

  f.loop_in_chunks(print_chunk_info, 0, f.total_volume());
}
