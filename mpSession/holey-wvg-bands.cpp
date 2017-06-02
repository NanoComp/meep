/***************************************************************/
/* C++ port of meep/examples/holey-wvg-bands.ctl               */
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
#include "mpSession.hpp"

int main(int argv, char *argc[])
{
  double eps = 13.0;  // dielectric constant of waveguide
  double w   = 1.2;   // width of waveguide
  double r   = 0.36;  // radius of holes

  int sx     = 12;    // size of cell in x direction
  int sy     = 12;    // size of cell in y direction (perpendicular to wvg.)
  int dpml   = 1;     // PML thickness (y direction only!)

  mpSession S;

  // (set! geometry-lattice (make lattice (size 1 sy no-size)))
  S.set_geometry_lattice(sx, sy);
  S.set_resolution=20;             // (set-param! resolution 20)

  // (set! geometry...
  //  (list (make block (center 0 0)
  //                    (size infinity w infinity)
  //                    (material (make dielectric (epsilon eps))))
  medium myMedium = make_medium(eps);
  vector3 center  = v3_zeros;
  vector3 size    = vector3_scale(w,v3_yaxisBar); // w*(inf 1 inf)
  S.add_block(myMedium,center,vector3_scale(w, v3_yaxisBar));

  //        (make cylinder (center 0 0) (radius r) (height infinity) (material air))))
  S.add_cylinder(air,center,r,inf);

  // (set! pml-layers (list (make pml (direction Y) (thickness dpml))))
  S.add_pml_layer(dpml, Y);

  double fcen=0.25;  // pulse center frequency
  double df=1.5;     // pulse freq. width: large df = short impulse
  vector3 src_point = vector3_scale(0.1234, v3_xaxis);
  S.add_gaussian_src(fcen, df, Hz, src_point);

  //(set! symmetries (list (make mirror-sym (direction Y) (phase -1))))
  S.add_mirror_symmetry(Y, -1.0);

  bool fixedkx=false; // if true, do run at this kx and get fields

  if (fixedkx)
   { 
     //  (run-sources+ 
     //  300 (at-beginning output-epsilon)
     //  (after-sources (harminv Hz (vector3 0.1234 0) fcen df)))

     S.add_step_func(output_epsilon, 0, AT_BEGINNING);
     harminv_data myhdata={.c=Hz, .pt=src_point, .fcen=fcen, .df=df};
     S.add_step_func(do_harminv, &myhdata, AFTER_SOURCES);

     S.run_until(300.0);

     // (run-until (/ 1 fcen) (at-every (/ 1 fcen 20) output-hfield-z)))
     S.add_step_func(outputhfield_z, 0, AT_EVERY, fcen/20.0);
     S.run_until(1.0/fcen);
  } // if fixedkx
 else
  {
    // (run-k-points 300 (interpolate k-interp (list (vector3 0) (vector3 0.5)))))
    int k_interp=19;            // k-points to interpolate, otherwise
    double kMin=0.0, kMax=0.5; 
    double kDelta=(kMax-kMin)/((double)(k_interp));
    for(int nk=0; nk<=kinterp; nk++)
     S.run_k_point(300.0, vector3_scale( nk*k_interp, v3_yaxis)
  };

}
