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

typedef std::complex<double> cdouble;

vector3 v3(double x, double y=0.0, double z=0.0)
{
  vector3 v;
  v.x=x; v.y=y; v.z=z;
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
void Run(bool Pulse)
{
  /***************************************************************/
  /* initialize geometry                                         */
  /***************************************************************/
  double n=3.4;     // index of waveguide
  double w=1.0;     // width of waveguide
  double r=1.0;     // inner radius of ring
  double pad=4;     // padding between waveguide and edge of PML
  double dpml=2;    // thickness of PML

  double resolution = 10.0;

  double sxy        = 2.0*(r+w+pad+dpml);  // cell size
  geometry_lattice.size.x=sxy;
  geometry_lattice.size.y=sxy;
  geometry_lattice.size.z=0.0;
  grid_volume gv = voltwo(sxy, sxy, resolution);
  gv.center_origin();
  symmetry sym=identity(); // mirror(Y, gv);
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  /***************************************************************/
  /* add objects                                                 */
  /***************************************************************/
  meep_geom::material_type dielectric = meep_geom::make_dielectric(n*n);
  geometric_object objects[2];
  vector3 v3zero = {0.0, 0.0, 0.0};
  vector3 zaxis  = {0.0, 0.0, 1.0};
  objects[0] = make_cylinder(dielectric, v3zero, r+w, ENORMOUS, zaxis);
  objects[1] = make_cylinder(meep_geom::vacuum,  v3zero, r, ENORMOUS, zaxis);
  geometric_object_list g={ 2, objects };
  meep_geom::set_materials_from_geometry(&the_structure, g);
  fields f(&the_structure);
  f.step(); // single timestep to trigger internal initialization
  f.output_hdf5(Dielectric, f.v);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double fcen = 0.15;  // ; pulse center frequency
  double df   = 0.1;   // ; df
  vec x0(r+0.1,0.0);   // ; source location
  if(Pulse)
   { 
     f.add_point_source(Ez, gaussian_src_time(fcen,df), x0);

     component components[6] = {Ex, Ey, Ez, Hx, Hy, Hz};
     dft_fields dftFields    = f.add_dft_fields(components, 6, f.v, fcen, fcen, 1);
     dft_flux dftFlux        = f.add_dft_flux(X, f.v, fcen, fcen, 1);

     while( f.round_time() < f.last_source_time() + 100.0 )
      f.step();

     f.output_dft(dftFlux,   "dft_flux");
     f.output_dft(dftFields, "dft_fields");
   }
  else
   { 
     f.add_point_source(Ez, continuous_src_time(fcen), x0);
     f.solve_cw();
     f.output_hdf5(Ez, f.v);
     f.output_hdf5(Hx, f.v);
     f.output_hdf5(Hy, f.v);
   }

  // success if we made it here
  // return 0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  initialize mpi(argc, argv);
  Run(true);
  Run(false);
}
