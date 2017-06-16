/***************************************************************/
/* C++ port of meep/examples/ring.ctl, using the               */
/* "low-level" meep C++ interface stack, which consists of     */
/* libmeep_geom + libctlgeom + libmeep                         */
/***************************************************************/
/*
; Calculating 2d ring-resonator modes, from the Meep tutorial.
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex>

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
int main(int argc, char *argv[])
{
  initialize mpi(argc, argv);

  double n=3.4;     // index of waveguide
  double w=1.0;     // width of waveguide
  double r=1.0;     // inner radius of ring

  double pad=4;     // padding between waveguide and edge of PML
  double dpml=2;    // thickness of PML

  double sxy        = 2.0*(r+w+pad+dpml);  // cell size
  double resolution = 10.0; 

  // (set-param! resolution 10)
  // (set! geometry-lattice (make lattice (size sxy sxy no-size)))
  grid_volume gv = voltwo(sxy, sxy, resolution);
  gv.center_origin();

  // (set! symmetries (list (make mirror-sym (direction Y))))
  symmetry sym=mirror(Y, gv);

  // (set! pml-layers (list (make pml (thickness dpml))))
  // ; exploit the mirror symmetry in structure+source:
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  // ; Create a ring waveguide by two overlapping cylinders - later objects
  // ; take precedence over earlier objects, so we put the outer cylinder first.
  // ; and the inner (air) cylinder second.
  // (set! geometry (list
  //	(make cylinder (center 0 0) (height infinity)
  //		(radius (+ r w)) (material (make dielectric (index n))))
  // 	(make cylinder (center 0 0) (height infinity)
  // 		(radius r) (material air))))
  material_type dielectric = meep_geom::make_dielectric(n*n);
  geometric_object objects[2];
  vector3 v3zero = {0.0, 0.0, 0.0};
  vector3 zaxis  = {0.0, 0.0, 1.0};
  objects[0] = make_cylinder(dielectric, v3zero, r+w, HUGE_VAL, zaxis);
  objects[1] = make_cylinder(meep_geom::air,  v3zero, r,   HUGE_VAL, zaxis);
  geometric_object_list g={ 2, objects }; 
  meep_geom::set_materials_from_geometry(&the_structure, g);
  fields f(&the_structure);

  // ; If we don't want to excite a specific mode symmetry, we can just
  // ; put a single point source at some arbitrary place, pointing in some
  // ; arbitrary direction.  We will only look for TM modes (E out of the plane).
  // (set! sources (list
  //              (make source
  //                (src (make gaussian-src (frequency fcen) (fwidth df)))
  //                (component Ez) (center (+ r 0.1) 0))))
  double fcen = 0.15;  // ; pulse center frequency
  double df   = 0.1;   // ; df
  gaussian_src_time src(fcen, df);
  volume v( vec(r+0.1,0.0), vec(0.0, 0.0) );
  f.add_volume_source(Ez, src, v);

  // (run-sources+ 300 
  // 	(at-beginning output-epsilon)
  // 	(after-sources (harminv Ez (vector3 (+ r 0.1)) fcen df)))
  double T = 300.0;
  double stop_time = f.last_source_time() + T;
  while( f.round_time() < stop_time )
   f.step();

  // TODO: translate call to harminv
  // int bands = do_harminv (... Ez, vec3(r+0.1), fcen, df) 

  // ; Output fields for one period at the end.  (If we output
  // ; at a single time, we might accidentally catch the Ez field
  // ; when it is almost zero and get a distorted view.)
  // (run-until (/ 1 fcen) (at-every (/ 1 fcen 20) output-efield-z))
  double DeltaT = 1.0/(20*fcen);
  double NextOutputTime = f.round_time() + DeltaT;
  while( f.round_time() < 1.0/fcen )
   { f.step();
     if ( f.round_time() >= NextOutputTime  )
      { f.output_hdf5(Ez, f.total_volume());
        NextOutputTime += DeltaT;
      };
   };

  // success if we made it here
  exit(0);

}
