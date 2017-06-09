/***************************************************************/
/* C++ port of meep/examples/bend-flux.ctl, using the          */
/* "low-level" meep C++ interface stack, which consists of     */
/* libmeep_geom + libctlgeom + libmeep                         */
/***************************************************************/

/*
; From the Meep tutorial: transmission around a 90-degree waveguide
; bend in 2d.
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex>

#include "meep.hpp"

#include "ctl-math.h"
#include "ctlgeom.h"

#include "meep_geom.hpp"

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
int main(int argv, char *argc[])
{

  double sx=16.0;       // size of cell in X direction
  double sy=32.0;       // size of cell in Y direction
  double pad=4.0;       // padding distance between waveguide and cell edge
  double w=1.0;         // width of waveguide
  double resolution=10; // (set-param! resolution 10)

  // (set! geometry-lattice (make lattice (size sx sy no-size)))
  // (set! pml-layers (list (make pml (thickness 1.0))))
  grid_volume gv = voltwo(sx, sy, resolution);
  structure the_structure(gv, dummy_eps, pml(1.0));

  double wvg_ycen = -0.5*(sy - w - 2.0*pad); // y center of horiz. wvg
  double wvg_xcen =  0.5*(sx - w - 2.0*pad); //x center of vert. wvg

   //  (set! geometry
   //      (if no-bend?
   //	  (list
   //	  (list
   //	   (make block 
   //	     (center wvg-xcen (* 0.5 pad))
   //	     (size w (- sy pad) infinity)
   //	     (material (make dielectric (epsilon 12)))))))
   //	   (make block 
   //	     (center 0 wvg-ycen)
   //	     (size infinity w infinity)
   //	     (material (make dielectric (epsilon 12)))))
   //	   (make block 
   //	     (center (* -0.5 pad) wvg-ycen)
   //	     (size (- sx pad) w infinity)
   //	     (material (make dielectric (epsilon 12))))
   vector e1=v3(1.0, 0.0, 0.0);
   vector e2=v3(0.0, 1.0, 0.0);
   vector e3=v3(0.0, 0.0, 1.0);

   material_type dielectric = make_dielectric(12.0);
   bool no_bend=false; // if true, have straight waveguide, not bend

   if (no_bend)
    {
      geometric_object objects[1];
      objects[0] = make_block(dielectric,
                              v3(0.0, wvg_ycen),
                              e1, e2, e3,
                              v3(HUGE_VAL, w, HUGE_VAL)
                             );
      geometric_object_list g={ objects, 1 };
      set_materials_from_geometry(the_structure, g);
    }
   else
    {
      geometric_object objects[3];
      objects[0] = make_block(dielectric,
                              v3(-0.5*pad, wvg_ycen),
                              e1, e2, e3,
                              v3(-sx*pad, w, HUGE_VAL)
                             );
   
      objects[1] = make_block(dielectric,
                              v3(0.0, wvg_ycen),
                              e1, e2, e3,
                              v3(-sx*pad, w, HUGE_VAL)
                             );
   
      objects[2] = make_block(dielectric,
                              v3(-0.5*pad, wvg_ycen),
                              e1, e2, e3,
                              v3(-sx*pad, w, HUGE_VAL)
                             );
   
      geometric_object_list g={ objects, 3 };
      set_materials_from_geometry(the_structure, g);
    };

  fields f(&the_structure);

  //  (set! sources (list
  //	       (make source 
  //		 (src (make gaussian-src (frequency fcen) (fwidth df)))
  //		 (component Ez)
  //		 (center (+ 1 (* -0.5 sx)) wvg-ycen)
  //		 (size 0 w))))
  //
  double fcen = 0.15;  // ; pulse center frequency
  double df   = 0.1;   // ; df
  gaussian_src_time src(fcen, df);
  f.add_volume_source(Ez, src, vec(1-0.5*sx, wvg-ycen), vec(0,w));

  //
  // didn't finish porting the rest yet
  //
  int nfreq=100; //  number of frequencies at which to compute flux

  //(define trans ; transmitted flux
  //      (add-flux fcen df nfreq
  //		(if no-bend?
  //		    (make flux-region
  //		      (center (- (/ sx 2) 1.5) wvg-ycen) (size 0 (* w 2)))
  //		    (make flux-region
  //		      (center wvg-xcen (- (/ sy 2) 1.5)) (size (* w 2) 0)))))

  //(define refl ; reflected flux
  //      (add-flux fcen df nfreq
  //		(make flux-region 
  //		  (center (+ (* -0.5 sx) 1.5) wvg-ycen) (size 0 (* w 2)))))
  //

  //; for normal run, load negated fields to subtract incident from refl. fields
  //(if (not no-bend?) (load-minus-flux "refl-flux" refl))
  //
  //(run-sources+ 
  // (stop-when-fields-decayed 50 Ez
  //			   (if no-bend? 
  //			       (vector3 (- (/ sx 2) 1.5) wvg-ycen)
  //			       (vector3 wvg-xcen (- (/ sy 2) 1.5)))
  //			   1e-3)
  // (at-beginning output-epsilon))

  //; for normalization run, save flux fields for refl. plane
  //(if no-bend? (save-flux "refl-flux" refl))
  //
  //(display-fluxes trans refl)

}
