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
void bend_flux(bool no_bend)
{
  double sx=16.0;       // size of cell in X direction
  double sy=32.0;       // size of cell in Y direction
  double pad=4.0;       // padding distance between waveguide and cell edge
  double w=1.0;         // width of waveguide
  double resolution=10; // (set-param! resolution 10)

  // (set! geometry-lattice (make lattice (size sx sy no-size)))
  // (set! pml-layers (list (make pml (thickness 1.0))))
  grid_volume gv = voltwo(sx, sy, resolution);
  gv.center_origin();
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
   vector3 e1=v3(1.0, 0.0, 0.0);
   vector3 e2=v3(0.0, 1.0, 0.0);
   vector3 e3=v3(0.0, 0.0, 1.0);

   material_type dielectric = meep_geom::make_dielectric(12.0);
   if (no_bend)
    {
      geometric_object objects[1];
      objects[0] = make_block(dielectric,
                              v3(0.0, wvg_ycen),
                              e1, e2, e3,
                              v3(HUGE_VAL, w, HUGE_VAL)
                             );
      geometric_object_list g={ 1, objects };
      meep_geom::set_materials_from_geometry(&the_structure, g);
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
   
      geometric_object_list g={ 3, objects };
      meep_geom::set_materials_from_geometry(&the_structure, g);
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
  volume v( vec(1.0-0.5*sx, wvg_ycen), vec(0.0,w) );
  f.add_volume_source(Ez, src, v);


  //(define trans ; transmitted flux
  //      (add-flux fcen df nfreq
  //		(if no-bend?
  //		    (make flux-region
  //		      (center (- (/ sx 2) 1.5) wvg-ycen) (size 0 (* w 2)))
  //		    (make flux-region
  //		      (center wvg-xcen (- (/ sy 2) 1.5)) (size (* w 2) 0)))))
  int nfreq=100; //  number of frequencies at which to compute flux
 
  volume *trans_volume=
   no_bend ? new volume(vec(0.5*sx-1.5, wvg_ycen), vec(0.0, 2.0*w))
           : new volume(vec(wvg_xcen, 0.5*sy-1.5), vec(2.0*w, 0.0));
  dft_flux trans=f.add_dft_flux( Z,
                                 *trans_volume,
                                 fcen-0.5*df, fcen+0.5*df, nfreq);
  
  //(define refl ; reflected flux
  //      (add-flux fcen df nfreq
  //		(make flux-region 
  //		  (center (+ (* -0.5 sx) 1.5) wvg-ycen) (size 0 (* w 2)))))
  //
  volume refl_volume( vec(-0.5*sx+1.5, wvg_ycen), vec(0.0,2.0*w));
  dft_flux refl=f.add_dft_flux(Z,
                               refl_volume, 
                               fcen-0.5*df, fcen+0.5*df, nfreq);

  //; for normal run, load negated fields to subtract incident from refl. fields
  //(if (not no-bend?) (load-minus-flux "refl-flux" refl))
  const char *dataname="refl-flux";
  if (!no_bend)
   { refl.load_hdf5(f, dataname);
     refl.scale_dfts(-1.0);
   };

  //(run-sources+ 
  // (stop-when-fields-decayed 50 Ez
  //			   (if no-bend? 
  //			       (vector3 (- (/ sx 2) 1.5) wvg-ycen)
  //			       (vector3 wvg-xcen (- (/ sy 2) 1.5)))
  //			   1e-3)
  // (at-beginning output-epsilon))
  vec eval_point = no_bend ? vec(0.5*sx-1.5, wvg_ycen)
                           : vec(wvg_xcen, 0.5*sy - 1.5);
  double DeltaT=50.0, NextCheckTime = f.round_time() + DeltaT;
  double Tol=1.0e-3;
  double max_abs=0.0;
  bool Done=false;
  do
   {  
     f.step();

     // manually check fields-decayed condition
     double absEz = abs(f.get_field(Ez, eval_point));
     max_abs = fmax(max_abs, absEz);
     if ( f.round_time() >= NextCheckTime )
      { NextCheckTime += DeltaT;
        if (max_abs > 0.0 && absEz<Tol*absEz)
         Done=true;
      };
     printf("%+.2e %+.2e \n",f.round_time(),absEz);

   } while(Done);

  //; for normalization run, save flux fields for refl. plane
  //(if no-bend? (save-flux "refl-flux" refl))
  if (no_bend)
   refl.save_hdf5(f, dataname);

  //printf("Trans: %e\n",
  //(display-fluxes trans refl)

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argv, char *argc[])
{
  (void )argv; // currently unused
  (void )argc;

  bend_flux(true);
  bend_flux(false);

}
