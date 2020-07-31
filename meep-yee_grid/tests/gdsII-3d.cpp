/***************************************************************/
/* example of a 3D geometry defined by GDSII file              */
/***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <complex>

#include "meep.hpp"

#include "ctl-math.h"
#include "ctlgeom.h"

#include "meepgeom.hpp"

using namespace meep;

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

/***************************************************************/
/* GDSII layers on which various geometric entities live       */
/***************************************************************/
#define GEOM_LAYER 0       //  computational cell
#define OXIDE_BULK_LAYER 1 //  oxide layer (bulk, i.e. oxide region)
#define OXIDE_VIA_LAYER 2  //  oxide layer (vias)
#define SILICON_LAYER 3    //  hexagon, rectangle
#define SLICE_LAYER 4      //  volumes for outputting epsilon slices

/***************************************************************/
/* layer thicknesses and materials *****************************/
/***************************************************************/
#define OXIDE_ZMIN 0.0
#define OXIDE_ZMAX 1.0
#define OXIDE_Z0 0.5 * (OXIDE_ZMIN + OXIDE_ZMAX)
#define OXIDE_EPS 2.2

#define SILICON_ZMIN OXIDE_ZMAX
#define SILICON_ZMAX 0.75
#define SILICON_Z0 0.5 * (SILICON_ZMIN + SILICON_ZMAX)
#define SILICON_EPS 12.0

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) {
  initialize mpi(argc, argv);

  const char *GDSIIFile = "gdsII-3d.gds";

  // set computational cell
  double dpml = 1.0;
  double resolution = 10.0;
  grid_volume gv = meep_geom::set_geometry_from_GDSII(resolution, GDSIIFile, GEOM_LAYER);
  structure the_structure(gv, dummy_eps, pml(dpml));

  // oxide layer, part 1: bulk of oxide layer
  meep_geom::material_type oxide = meep_geom::make_dielectric(OXIDE_EPS);
  geometric_object oxide_bulk_prism =
      meep_geom::get_GDSII_prism(oxide, GDSIIFile, OXIDE_BULK_LAYER, OXIDE_ZMIN, OXIDE_ZMAX);

  // oxide layer, part 2: via in oxide layer
  geometric_object oxide_via_prism = meep_geom::get_GDSII_prism(
      meep_geom::vacuum, GDSIIFile, OXIDE_VIA_LAYER, OXIDE_ZMIN, OXIDE_ZMAX);

  // silicon layer
  meep_geom::material_type silicon = meep_geom::make_dielectric(SILICON_EPS);
  geometric_object_list silicon_prisms =
      meep_geom::get_GDSII_prisms(silicon, GDSIIFile, SILICON_LAYER, SILICON_ZMIN, SILICON_ZMAX);

  // merge all prisms into a single geometric_object_list and instantiate meep geometry
  geometric_object_list all_prisms;
  all_prisms.num_items = 1 + 1 + silicon_prisms.num_items;
  all_prisms.items = new geometric_object[all_prisms.num_items];
  all_prisms.items[0] = oxide_bulk_prism;
  all_prisms.items[1] = oxide_via_prism;
  for (int n = 0; n < silicon_prisms.num_items; n++)
    all_prisms.items[2 + n] = silicon_prisms.items[n];
  meep_geom::set_materials_from_geometry(&the_structure, all_prisms);
  fields f(&the_structure);

  // define volumes for source and flux-monitor regions
  volume v1 =
      meep_geom::get_GDSII_volume(GDSIIFile, "yzplane", SLICE_LAYER, OXIDE_ZMIN, SILICON_ZMAX);
  volume v2 = meep_geom::get_GDSII_volume(GDSIIFile, "xyplane", SLICE_LAYER, OXIDE_Z0, OXIDE_Z0);
  volume v3 =
      meep_geom::get_GDSII_volume(GDSIIFile, "xyplane", SLICE_LAYER, SILICON_Z0, SILICON_Z0);

  f.step();
  f.output_hdf5(Dielectric, v1, 0, false, true, "v1");
  f.output_hdf5(Dielectric, v2, 0, false, true, "v2");
  f.output_hdf5(Dielectric, v3, 0, false, true, "v3");
}
