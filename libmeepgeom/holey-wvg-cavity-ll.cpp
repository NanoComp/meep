/***************************************************************/
/* C++ port of meep/examples/holey-wvg-cavity.ctl, using the   */
/* "low-level" meep C++ interface stack, which consists of     */
/* libmeep_geom + libctlgeom + libmeep.                        */
/*                                                             */
/* usage:                                                      */
/*  holey-wvg-cavity-ll [--compute-modes]                      */
/*                                                             */
/***************************************************************/

/*
; Meep Tutorial: TE transmission and reflection through a cavity
; formed by a periodic sequence of holes in a dielectric waveguide,
; with a defect formed by a larger spacing between one pair of holes.

; This structure is based on one analyzed in:
;    S. Fan, J. N. Winn, A. Devenyi, J. C. Chen, R. D. Meade, and
;    J. D. Joannopoulos, "Guided and defect modes in periodic dielectric
;    waveguides," J. Opt. Soc. Am. B, 12 (7), 1267-1272 (1995).
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
/* compute_modes == true  -> compute resonant modes            */
/*                  false -> compute transmission spectrum     */
/***************************************************************/
void holey_wvg_cavity(bool compute_modes=false)
{
  printf("Running holey_wvg_cavity to compute %s...\n", 
          compute_modes ? "resonant modes" : "transmission spectrum");

  // Some parameters to describe the geometry:
  double Eps=13.0;   // dielectric constant of waveguide
  double w=1.0;      // width of waveguide
  double r=0.36;     // radius of holes
  double d=1.4;      // defect spacing (ordinary spacing = 1)
  int N=3;           // number of holes on either side of defect

  // The cell dimensions
  double sy=6.0;     // size of cell in y direction (perpendicular to wvg.)
  double pad=2.0;    // padding between last hole and PML edge
  double dpml=1.0;   // PML thickness

  // compute-mode == true --> 
  //  (set! symmetries
  //   (list (make mirror-sym (direction Y) (phase -1))
  //         (make mirror-sym (direction X) (phase -1))))
  // compute-mode == false -->
  // (set! symmetries (list (make mirror-sym (direction Y) (phase -1))
  symmetry sym = -mirror(Y,gv);
  if (compute_modes)
   sym += -mirror(X,gv);

  // (define sx (+ (* 2 (+ pad dpml N)) d -1))
  // (set! geometry-lattice (make lattice (size sx sy no-size)))
  // (set! pml-layers (list (make pml (thickness dpml))))
  // (set-param! resolution 20)
  double sx = 2.0*(pad + dpml + N) + d -1.0; // size of cell in x dir
  double resolution=20.0;
  geometry_lattice.size.x=sx;
  geometry_lattice.size.y=sy;
  geometry_lattice.size.z=0.0;
  grid_volume gv = voltwo(sx, sy, resolution);
  gv.center_origin();
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  // (set! geometry
  //    (append ; combine lists of objects:
  //     (list (make block (center 0 0) (size infinity w infinity)
  // (material (make dielectric (epsilon eps)))))
  //      (geometric-object-duplicates (vector3 1 0) 0 (- N 1)
  // (make cylinder (center (/ d 2) 0) (radius r) (height infinity)
  // (material air)))
  // (geometric-object-duplicates (vector3 -1 0) 0 (- N 1)
  // (make cylinder (center (/ d -2) 0) (radius r) (height infinity)
  // (material air)))))
  material_type dielectric = meep_geom::make_dielectric(eps);
  geometric_object objects[7];
  vector3 origin = v3(0.0,   0.0,  0.0);
  vector3 xhat   = v3(1.0,   0.0,  0.0);
  vector3 yhat   = v3(0.0,   1.0,  0.0);
  vector3 zhat   = v3(0.0,   0.0,  1.0);
  vector3 center = v3(0.5*d, 0.0,  0.0);
  vector3 delta  = v3(1.0,   0.0,  0.0);
  vector3 size   = v3(ENORMOUS, w, ENORMOUS);

  objects[0] = make_block(dielectric, origin, xhat, yhat, zhat, size);
  objects[1] = make_cylinder(vacuum, center+0*delta, r, height, zhat);
  objects[2] = make_cylinder(vacuum, center+1*delta, r, height, zhat);
  objects[3] = make_cylinder(vacuum, center+2*delta, r, height, zhat);
  objects[4] = make_cylinder(vacuum, -center-0*delta, r, height, zhat);
  objects[5] = make_cylinder(vacuum, -center-1*delta, r, height, zhat);
  objects[6] = make_cylinder(vacuum, -center-2*delta, r, height, zhat);
  geometric_object_list g={ 7, objects };
  meep_geom::set_materials_from_geometry(&the_structure, g);
  fields f(&the_structure);

  //
  double fcen = 0.25;  // pulse center frequency
  double df   = 0.2;   // pulse width (in frequency)
  int nfreq   = 500;   // number of frequencies at which to compute flux
  
  if (compute_modes)
   { 
     //  (set! sources (list
     //  	(make source
     //		(src (make gaussian-src (frequency fcen) (fwidth df)))
     // 		(component Hz) (center 0 0))))
     gaussian_src_time src(fcen, df);
     component src_cmpt = Hz;
     f.add_point_source(src_cmpt, src, origin);

     // (run-sources+ 400
     //    (at-beginning output-epsilon)
     //    (after-sources (harminv Hz (vector3 0) fcen df)))
     f.output_hdf5(eps, f.total_volume());
     while( f.round_time() < f.last_source_time() )
      f.step();
     while( f.round_time() < f.last_source_time() + 400.0)
      f.step();

     // (run-until (/ 1 fcen) (at-every (/ 1 fcen 20) output-hfield-z))      
     double DeltaT = 1.0/(20*fcen);
     double NextOutputTime = f.round_time() + DeltaT;
     while( f.round_time() < 1.0/fcen )
      { f.step();
        if ( f.round_time() >= NextOutputTime  )
         { f.output_hdf5(Hz, f.total_volume());
           NextOutputTime += DeltaT;
         };
      };

   } // if (compute_modes)
  else
   { 
     // (set! sources (list
     //		     (make source
     //		       (src (make gaussian-src (frequency fcen) (fwidth df)))
     //		       (component Ey)
     //		       (center (+ dpml (* -0.5 sx)) 0)
     //		       (size 0 w))))
     gaussian_src_time src(fcen, df);
     component src_cmpt = Ey;
     vec src_center = vec(dpml - 0.5*sx, 0.0);
     vec src_size  = vec(0.0, w);
     volume src_vol(src_center, src_size);
     f.add_volume_source(src_cmpt, src, src_vol);
      
     // (define trans ; transmitted flux
     // (add-flux fcen df nfreq
     //	  (make flux-region
     // (center (- (* 0.5 sx) dpml 0.5) 0) (size 0 (* w 2)))))
      
     //    (run-sources+ (stop-when-fields-decayed 
     //		     50 Ey
     //		     (vector3 (- (* 0.5 sx) dpml 0.5) 0)
     //		     1e-3)
     //		    (at-beginning output-epsilon)
     //		    (during-sources
     //                     (in-volume (volume (center 0 0) (size sx 0))
     //                     (to-appended "hz-slice" (at-every 0.4 output-hfield-z)))))


     //      (display-fluxes trans) ; print out the flux spectrum
     //      )

   }; // if(compute_modes) ... else ... 

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  bool compute_modes=false;
  for(int narg=1; narg<argc; narg++)
   { if (!strcmp(argv[narg],"--compute-modes"))
      compute_modes=true;
     else
      meep::abort("unknown command-line option %s",argv[narg]);
   };

  holey_wvg_cavity(compute_modes);
   
}
